import random
from math import gcd, log2

import numpy as np

from crt_bus.analysis import plugin_mi
from crt_bus.crt import lcm
from crt_bus.system import (SystemBus, binary_tick_perm, quad_tick_perm, ternary_step_perm)
from crt_bus.zring import unpack

S = SystemBus()


def naive(state0, t):
    cur = {k: v.copy() for k, v in state0.items()}
    for tick in range(1, t + 1):
        cur["bin"] = cur["bin"][binary_tick_perm()]
        cur["quad"] = cur["quad"][quad_tick_perm()]
        cur["ter"] = cur["ter"][ternary_step_perm(tick)]
    return cur


def test_orbits_close_at_22_33_24_with_distinguishable_phases():
    assert S.orbit_periods == {"bin": 22, "ter": 33, "quad": 24}
    assert S.crystals["ter"].assumed and not S.crystals["bin"].assumed


def test_jump_ahead_equals_simulation():
    s0 = S.default_state(1)
    rnd = random.Random(6)
    for t in [0, 1, 11, 22, 24, 33, 264] + [rnd.randrange(1, 700) for _ in range(12)]:
        a, b = S.state_at(t, s0), naive(s0, t)
        assert all(np.array_equal(a[k], b[k]) for k in a), t


def test_bus_path_equals_direct_path_for_huge_ticks():
    s0 = S.default_state(2)
    rnd = random.Random(7)
    for _ in range(300):
        t = rnd.randrange(S.clock.range)
        direct = {k: s0[k][c.powers[t % c.period]] for k, c in S.crystals.items()}
        via_bus = S.state_at(t, s0, residues=S.clock.encode(t))
        assert all(np.array_equal(direct[k], via_bus[k]) for k in direct)


def test_wire_roundtrip_uses_one_Z6_symbol_for_one_bit_and_one_trit():
    s0 = S.default_state(3)
    w = S.wire(12345, s0)
    r = S.read_wire(w)
    s = S.state_at(12345, s0)
    assert r["t"] == 12345 and r["padding_ok"]
    assert all(np.array_equal(r["state"][k], s[k]) for k in s)
    b2, b3 = unpack(w.A)
    assert np.array_equal(b2, s["bin"]) and np.array_equal(b3[:12], s["ter"])


def test_clock_corruption_is_corrected_and_lane_corruption_is_localized():
    s0 = S.default_state(4)
    rnd = random.Random(8)
    for _ in range(40):
        t = rnd.randrange(S.clock.range)
        w = S.wire(t, s0)
        assert S.consistency_report(w, s0)["consistent"]
        i = rnd.randrange(len(S.clock.all_moduli))
        m = S.clock.all_moduli[i]
        c = list(w.clock)
        c[i] = (c[i] + rnd.randrange(1, m)) % m
        w.clock = tuple(c)
        rep = S.consistency_report(w, s0)
        assert rep["t"] == t and rep["consistent"] and rep["clock_was_corrupted"]
    w = S.wire(777, s0)
    w.A = w.A.copy()
    w.A[3] = (w.A[3] + 3) % 6                              # flip one BIT in the binary lane
    rep = S.consistency_report(w, s0)
    assert rep["lanes"] == {"bin": False, "ter": True, "quad": True} and not rep["consistent"]


def test_padding_integrity_is_checked():
    s0 = S.default_state(5)
    w = S.wire(10, s0)
    w.B = w.B.copy()
    w.B[0] = (w.B[0] + 4) % 6                              # puts a trit in the unused GF(3) half
    assert not S.read_wire(w)["padding_ok"]


def test_joint_phase_check_detects_desync_and_recovers_tick_mod_264():
    rnd = random.Random(9)
    ident = {k: np.arange(len(c.powers[0])) for k, c in S.crystals.items()}
    for _ in range(200):
        t = rnd.randrange(10 ** 6)
        perms = {k: c.powers[t % c.period] for k, c in S.crystals.items()}
        ph = S.phases_from_perms(perms)
        assert ph == {k: t % c.period for k, c in S.crystals.items()}
        assert S.joint_tick_mod(ph) == (t % 264, 264)
    bad = {"bin": 1, "ter": 0, "quad": 0}                 # needs 1 = 0 (mod 11)
    assert S.joint_tick_mod(bad) is None
    assert S.phases_from_perms({"bin": np.arange(24)[::-1]})["bin"] is None


def test_mutual_information_law_on_the_actual_crystal_orbits():
    for a, b in (("bin", "ter"), ("bin", "quad"), ("ter", "quad")):
        ca, cb = S.crystals[a], S.crystals[b]
        L = lcm(ca.period, cb.period)
        xs = [ca.powers[t % ca.period].tobytes() for t in range(L)]
        ys = [cb.powers[t % cb.period].tobytes() for t in range(L)]
        law = log2(gcd(ca.period, cb.period))
        assert abs(plugin_mi(xs, ys) - law) < 1e-9
        # a coarse observable (one coordinate) can only carry LESS information (data processing)
        xs0 = [int(ca.powers[t % ca.period][0]) for t in range(L)]
        ys0 = [int(cb.powers[t % cb.period][0]) for t in range(L)]
        assert plugin_mi(xs0, ys0) <= law + 1e-9
