"""system.py -- the WHOLE system on one bus: clock lane + GF(2)/GF(3) state lanes.

Crystals (state vectors over their own field):
    bin  : GF(2)^24   binary Golay crystal, Floquet tick of order 22   (repo: truth/dmf/floquet.py)
    ter  : GF(3)^12   ternary crystal, schedule period 33              (ASSUMED realization, see below)
    quad : GF(2)^48   quadalic crystal, M = C2 o P47 o P23Q of order 24 (x -> -1/(2x+1) on P^1(F_47))

Wire format
    clock : 7 residues (channels 8,3,11,23,47 + redundant 53,59)
    A     : Z_6^24 = bin bits (GF(2) half) + ter trits (GF(3) half, first 12 positions)
    B     : Z_6^48 = quad bits (GF(2) half); GF(3) half unused (must be 0)
One Z_6 symbol = one bit and one trit (CRT on the fields).  Field-lane operations never mix.

State at tick t is  state0[powers[t mod period]]  where powers[k] is the cumulative gather
permutation after k ticks; the phase t mod period is READ FROM THE CLOCK RESIDUES, so jumping
to any tick costs the same regardless of t.  Storage = 22 + 33 + 24 = 79 permutations, not the
264 joint states.

ASSUMED: the exact ternary generators (P_mult, P_23^T, H_T schedule) are not in the extracted
corpus.  The realization used here (12-cycle; x -> 2x on 11 points; block 3-cycles; H_T every 11th
tick) reproduces the corpus numbers (base operator order 11, schedule period 33).  Replace
`ternary_step_perm` with the real generators when found; nothing else changes.
"""
from dataclasses import dataclass

import numpy as np

from .clock_bus import ClockBus
from .crt import crt_solve
from .zring import pack, unpack

CLOCK_NAME = {"bin": "bin_orbit", "ter": "ter_orbit", "quad": "quad_orbit"}
FIELD = {"bin": 2, "ter": 3, "quad": 2}
SIZE = {"bin": 24, "ter": 12, "quad": 48}


# ------------------------------------------------------------------ operators (gather permutations)
def binary_tick_perm():
    """Gather permutation of one Floquet tick: P23 -> P11 -> H_L -> S_inv (repo convention)."""
    inv = [0] + [pow(j, -1, 23) for j in range(1, 23)]
    v = list(range(24))
    w = v[:]; w[0] = v[22]; w[1:23] = v[0:22]; v = w
    w = v[:]
    for j in range(23):
        w[j] = v[(12 * j) % 23]
    v = w; w = v[:]
    for j in range(1, 23):
        w[j] = v[inv[j]]
    v = w; w = v[:]
    for j in range(1, 23):
        w[j] = v[(-inv[j]) % 23]
    return np.array(w, dtype=np.int64)


def quad_tick_perm():
    """Gather permutation of M = C2 o P47 o P23Q  (x -> -1/(2x+1) on P^1(F_47), inf = 47)."""
    scatter = [0] * 48
    for x in range(47):
        y = (2 * x + 1) % 47
        scatter[x] = 47 if y == 0 else (-pow(y, -1, 47)) % 47
    scatter[47] = 0
    gather = [0] * 48
    for i, j in enumerate(scatter):
        gather[j] = i
    return np.array(gather, dtype=np.int64)


def ternary_step_perm(t):
    """ASSUMED realization.  t is the 1-based tick index; H_T is applied on every 11th tick."""
    pcyc = [(i + 1) % 12 for i in range(12)]
    pmult = [(2 * i) % 11 for i in range(11)] + [11]
    A = np.array([pmult[pcyc[i]] for i in range(12)], dtype=np.int64)
    HT = np.zeros(12, dtype=np.int64)
    for b in range(4):
        for j in range(3):
            HT[3 * b + j] = 3 * b + (j + 1) % 3
    return A[HT] if t % 11 == 0 else A


# ------------------------------------------------------------------ a crystal
class Crystal:
    def __init__(self, name, n, step_perm, expected_period, assumed=False):
        self.name, self.n, self.assumed = name, n, assumed
        ident = np.arange(n)
        powers = [ident]
        for t in range(1, 10 * expected_period + 2):
            powers.append(powers[-1][step_perm(t)])          # c_t[j] = c_{t-1}[s_t[j]]
            if np.array_equal(powers[-1], ident):
                break
        self.period = len(powers) - 1
        if self.period != expected_period:
            raise AssertionError(f"{name}: orbit closes at {self.period}, expected {expected_period}")
        self.powers = np.array(powers[:-1])                  # powers[k], k = 0..period-1
        keys = {p.tobytes() for p in self.powers}
        assert len(keys) == self.period, "phases are not distinguishable (repeated permutation)"
        self._index = {p.tobytes(): k for k, p in enumerate(self.powers)}

    def apply(self, state, t):
        return np.asarray(state)[self.powers[t % self.period]]

    def locate_perm(self, perm):
        """Phase k with powers[k] == perm (unique), or None if perm is not on the orbit."""
        return self._index.get(np.asarray(perm, dtype=np.int64).tobytes())


@dataclass
class Wire:
    clock: tuple
    A: np.ndarray
    B: np.ndarray


class SystemBus:
    def __init__(self, clock=None):
        self.clock = clock or ClockBus()
        self.crystals = {
            "bin": Crystal("bin", 24, lambda t: binary_tick_perm(), 22),
            "ter": Crystal("ter", 12, ternary_step_perm, 33, assumed=True),
            "quad": Crystal("quad", 48, lambda t: quad_tick_perm(), 24),
        }
        self.orbit_periods = {k: c.period for k, c in self.crystals.items()}
        self._unified = None

    # ------------------------------------------------------------ states
    @staticmethod
    def default_state(seed=0):
        rng = np.random.default_rng(seed)
        return {"bin": rng.integers(0, 2, 24), "ter": rng.integers(0, 3, 12), "quad": rng.integers(0, 2, 48)}

    def state_at(self, t, state0, residues=None):
        """State at tick t. If residues are given, phases come from them (bus path)."""
        res = self.clock.encode(t) if residues is None else residues
        out = {}
        for k, c in self.crystals.items():
            phase = self.clock.phase(res, CLOCK_NAME[k])
            out[k] = np.asarray(state0[k])[c.powers[phase]]
        return out

    # ------------------------------------------------------------ wire
    def wire(self, t, state0):
        s = self.state_at(t, state0)
        ter_padded = np.concatenate([s["ter"], np.zeros(12, dtype=np.int64)])
        return Wire(self.clock.encode_full(t), pack(s["bin"], ter_padded), pack(s["quad"], np.zeros(48, dtype=np.int64)))

    def read_wire(self, w):
        """Decode a wire: corrected tick (or None), lane states, and padding-integrity flags."""
        t = self.clock.correct(w.clock)
        b2, b3 = unpack(w.A)
        q2, q3 = unpack(w.B)
        return {
            "t": t,
            "state": {"bin": b2, "ter": b3[:12], "quad": q2},
            "padding_ok": bool(not b3[12:].any() and not q3.any()),
        }

    def consistency_report(self, w, state0):
        """Whole-system check: is every lane where the (error-corrected) clock says it should be?"""
        r = self.read_wire(w)
        if r["t"] is None:
            return {"t": None, "clock_ok": False, "lanes": {}, "padding_ok": r["padding_ok"], "consistent": False}
        expected = self.state_at(r["t"], state0)
        lanes = {k: bool(np.array_equal(r["state"][k], expected[k])) for k in self.crystals}
        return {"t": r["t"], "clock_ok": True, "clock_was_corrupted": tuple(w.clock) != self.clock.encode_full(r["t"]),
                "lanes": lanes, "padding_ok": r["padding_ok"],
                "consistent": all(lanes.values()) and r["padding_ok"]}

    # ------------------------------------------------------------ phase / joint checks
    def phases_from_perms(self, perms):
        """Locate each crystal on its orbit from an observed label-permutation. Missing/foreign -> None."""
        return {k: self.crystals[k].locate_perm(p) for k, p in perms.items()}

    def joint_tick_mod(self, phases):
        """t mod 264 consistent with observed orbit phases {bin, ter, quad}, or None (= desync)."""
        pairs = [(phases[k], self.orbit_periods[k]) for k in phases if phases[k] is not None]
        if len(pairs) < len(phases):
            return None
        return crt_solve(pairs)

    # ------------------------------------------------------------ unified payload code
    @property
    def unified_code(self):
        if self._unified is None:
            from .unified_code import UnifiedCode24
            self._unified = UnifiedCode24()
        return self._unified
