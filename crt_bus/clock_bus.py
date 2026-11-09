"""clock_bus.py -- the tick counter as a CRT residue bus.

Channels (8, 3, 11, 23, 47) are pairwise coprime; their product 285,384 = 2^3*3*11*23*47
equals the joint period of ALL six DMF clocks:

    orbit clocks   : binary 22 = 2*11     ternary 33 = 3*11     quadalic 24 = 8*3
    master clocks  : binary 46 = 2*23     ternary 759 = 3*11*23 quadalic 1081 = 23*47

A tick t is stored as five small residues.  Each clock reads its phase from just the
channels it needs (no extra counters).  Incrementing is carry-free (each channel +1 mod m).
Two redundant channels (53, 59), both larger than every information channel, make the
code MDS: any single corrupted channel is detected AND corrected.

Facts proven by tests (see tests/test_clock_bus.py):
  * phase(res, p) == t % p for every clock;
  * generalized CRT answers "when do these phases co-occur?" exactly;
  * only lcm/product of the phase space is reachable (1/66 of orbit-phase triples);
  * for two clocks sampled over a full joint period, MI = log2(gcd) bits.
"""
from math import gcd, log2

from .crt import crt_solve, decode_coprime, lcm, pairwise_coprime, prime_power_factors

DEFAULT_CHANNELS = (8, 3, 11, 23, 47)
DEFAULT_REDUNDANT = (53, 59)
DEFAULT_CLOCKS = {
    "bin_orbit": 22, "ter_orbit": 33, "quad_orbit": 24,
    "bin_master": 46, "ter_master": 759, "quad_master": 1081,
}


class ClockBus:
    def __init__(self, channels=DEFAULT_CHANNELS, redundant=DEFAULT_REDUNDANT, clocks=None):
        self.channels = tuple(channels)
        self.redundant = tuple(redundant)
        self.clocks = dict(DEFAULT_CLOCKS if clocks is None else clocks)
        self.all_moduli = self.channels + self.redundant
        if not pairwise_coprime(self.all_moduli):
            raise ValueError("channels (and redundant channels) must be pairwise coprime")
        self.range = 1
        for m in self.channels:
            self.range *= m
        self.total = 1
        for m in self.all_moduli:
            self.total *= m
        # which channel serves each prime-power part of each clock
        self._plan = {}
        for name, period in self.clocks.items():
            self._plan[name] = self._plan_for(period)

    # ---------------------------------------------------------------- reading phases
    def _plan_for(self, period):
        plan = []
        for q in prime_power_factors(period):
            idx = [i for i, m in enumerate(self.channels) if m % q == 0]
            if not idx:
                raise ValueError(f"period {period}: no channel is divisible by {q}")
            plan.append((idx[0], q))
        return plan

    def encode(self, t):
        return tuple(t % m for m in self.channels)

    def encode_full(self, t):
        return tuple(t % m for m in self.all_moduli)

    def increment(self, res, full=False):
        """Carry-free tick: every channel advances independently."""
        mods = self.all_moduli if full else self.channels
        return tuple((r + 1) % m for r, m in zip(res, mods))

    def phase(self, res, clock):
        """Phase of a named clock (or of an integer period) read from the residues only."""
        plan = self._plan[clock] if isinstance(clock, str) else self._plan_for(clock)
        return crt_solve([(res[i] % q, q) for i, q in plan])[0]

    def phases(self, res):
        return {name: self.phase(res, name) for name in self.clocks}

    # ---------------------------------------------------------------- queries
    @staticmethod
    def align(constraints):
        """constraints: {period: phase}. Returns (t0, L): all t = t0 (mod L) satisfy them; or None."""
        return crt_solve([(ph, per) for per, ph in constraints.items()])

    @staticmethod
    def consistent(constraints):
        return crt_solve([(ph, per) for per, ph in constraints.items()]) is not None

    @staticmethod
    def reachable_fraction(periods):
        """Fraction of the phase product space that can actually occur = lcm / product."""
        prod = 1
        for p in periods:
            prod *= p
        return lcm(*periods) / prod

    @staticmethod
    def mi_bits(n1, n2):
        """Mutual information (bits) between two clocks' phases over a full joint period."""
        return log2(gcd(n1, n2))

    # ---------------------------------------------------------------- fault tolerance
    def check_mds(self):
        """Any len(channels) of the moduli must multiply to >= range (=> distance r+1)."""
        from itertools import combinations
        k = len(self.channels)
        for combo in combinations(self.all_moduli, k):
            p = 1
            for m in combo:
                p *= m
            if p < self.range:
                return False
        return True

    def correct(self, res_full):
        """Return the tick (in [0, range)) consistent with all but at most one channel, else None."""
        x = decode_coprime(res_full, self.all_moduli) % self.total
        if x < self.range:
            return x
        cands = set()
        for i in range(len(self.all_moduli)):
            mods = [m for k, m in enumerate(self.all_moduli) if k != i]
            sub = [r for k, r in enumerate(res_full) if k != i]
            y = decode_coprime(sub, mods)
            if y < self.range:
                cands.add(y)
        return cands.pop() if len(cands) == 1 else None
