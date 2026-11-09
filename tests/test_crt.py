import random
from math import gcd

from crt_bus.crt import crt_merge, crt_solve, decode_coprime, egcd, lcm, pairwise_consistent


def test_egcd_identity():
    rnd = random.Random(1)
    for _ in range(500):
        a, b = rnd.randrange(1, 10_000), rnd.randrange(1, 10_000)
        g, x, y = egcd(a, b)
        assert g == gcd(a, b) and a * x + b * y == g


def test_merge_matches_brute_force_including_noncoprime():
    rnd = random.Random(2)
    for _ in range(400):
        n1, n2 = rnd.randrange(1, 40), rnd.randrange(1, 40)
        a1, a2 = rnd.randrange(n1), rnd.randrange(n2)
        L = lcm(n1, n2)
        brute = [t for t in range(L) if t % n1 == a1 and t % n2 == a2]
        got = crt_merge(a1, n1, a2, n2)
        if brute:
            assert got == (brute[0], L) and len(brute) == 1
        else:
            assert got is None


def test_solvable_iff_pairwise_consistent():
    rnd = random.Random(3)
    for _ in range(300):
        ms = [rnd.randrange(2, 30) for _ in range(3)]
        pairs = [(rnd.randrange(m), m) for m in ms]
        L = lcm(*ms)
        solvable = any(all(t % m == a for a, m in pairs) for t in range(L))
        assert (crt_solve(pairs) is not None) == solvable == pairwise_consistent(pairs)


def test_decode_coprime_roundtrip():
    mods = (8, 3, 11, 23, 47)
    for t in range(0, 285384, 997):
        assert decode_coprime([t % m for m in mods], mods) == t
