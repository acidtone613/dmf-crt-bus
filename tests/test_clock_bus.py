import random
from math import gcd, log2

from crt_bus.analysis import plugin_mi
from crt_bus.clock_bus import DEFAULT_CLOCKS, ClockBus
from crt_bus.crt import crt_solve, lcm


def test_all_six_clocks_read_from_five_residues():
    b, rnd = ClockBus(), random.Random(4)
    assert b.range == 285384 == lcm(*DEFAULT_CLOCKS.values())
    for _ in range(3000):
        t = rnd.randrange(b.range)
        res = b.encode(t)
        for name, p in DEFAULT_CLOCKS.items():
            assert b.phase(res, name) == t % p


def test_increment_is_carry_free_and_correct():
    b = ClockBus()
    res = b.encode(0)
    for t in range(1, 1000):
        res = b.increment(res)
        assert res == b.encode(t)
    assert b.increment(b.encode(b.range - 1)) == b.encode(0)      # wraps at the joint period


def test_alignment_queries_match_brute_force():
    b = ClockBus()
    assert b.align({46: 23, 759: 0, 1081: 0}) == (35673, 71346)
    hits = [t for t in range(71346) if t % 46 == 23 and t % 759 == 0 and t % 1081 == 0]
    assert hits == [35673]
    assert b.align({46: 0, 759: 0, 1081: 0}) == (0, 71346)
    assert b.align({46: 1, 759: 0}) is None                       # needs 1 = 0 (mod 23)
    assert not any(t % 46 == 1 and t % 759 == 0 for t in range(1518))


def test_only_one_in_66_orbit_phase_triples_is_reachable():
    ok = sum(crt_solve([(i, 22), (j, 33), (k, 24)]) is not None
             for i in range(22) for j in range(33) for k in range(24))
    assert ok == 264 and 22 * 33 * 24 == 17424
    assert abs(ClockBus.reachable_fraction([22, 33, 24]) - 1 / 66) < 1e-12


def test_mutual_information_between_clocks_is_log2_gcd():
    for a, b in ((22, 33), (22, 24), (33, 24), (46, 759), (759, 1081), (47, 23)):
        L = lcm(a, b)
        xs, ys = [t % a for t in range(L)], [t % b for t in range(L)]
        assert abs(plugin_mi(xs, ys) - log2(gcd(a, b))) < 1e-9
        assert abs(ClockBus.mi_bits(a, b) - log2(gcd(a, b))) < 1e-12


def test_redundant_channels_are_mds_and_correct_any_single_error():
    b, rnd = ClockBus(), random.Random(5)
    assert b.check_mds()
    for _ in range(40):
        t = rnd.randrange(b.range)
        full = list(b.encode_full(t))
        assert b.correct(tuple(full)) == t                        # no error
        for i, m in enumerate(b.all_moduli):
            for delta in range(1, m):                             # EVERY wrong residue in EVERY channel
                bad = full[:]
                bad[i] = (bad[i] + delta) % m
                assert b.correct(tuple(bad)) == t


def test_values_beyond_range_are_flagged_not_wrapped():
    b = ClockBus()
    over = tuple((b.range + 5000) % m for m in b.all_moduli)
    assert b.correct(over) is None
