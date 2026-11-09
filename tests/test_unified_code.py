import itertools

import numpy as np

from _compat import skip
from crt_bus.unified_code import (INF, UnifiedCode24, automorphism_generators, automorphism_orders,
                                  build_golay24, build_xqr23, element_from_word)
from crt_bus.zring import ZMonomial, pack, unpack

U = UnifiedCode24()
EXPECTED_ORDERS = {(1, 1), (2, 1), (3, 3), (4, 2), (6, 3), (8, 4), (11, 11), (12, 6), (22, 11), (23, 23), (24, 12), (46, 23)}


def _rank_mod(G, p):
    from crt_bus.unified_code import _rref_mod
    return len(_rref_mod(G, p)[1])


def test_both_codes_are_self_dual_rank_12_length_24():
    for G, p in ((U.G2, 2), (U.G3, 3)):
        assert G.shape == (12, 24) and _rank_mod(G, p) == 12
        assert not ((G @ G.T) % p).any()                  # self-orthogonal + dimension 12 => self-dual


def test_weight_distributions_and_minimum_distances():
    m2 = np.array([[(m >> i) & 1 for i in range(12)] for m in range(4096)])
    w2 = np.bincount(((m2 @ U.G2) % 2).sum(axis=1))
    assert {w: int(c) for w, c in enumerate(w2) if c} == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    idx = np.arange(3 ** 12)
    m3 = np.stack([(idx // 3 ** k) % 3 for k in range(12)], axis=1)
    w3 = np.bincount(((m3 @ U.G3) % 3 != 0).sum(axis=1))
    assert {w: int(c) for w, c in enumerate(w3) if c} == {0: 1, 9: 4048, 12: 61824, 15: 242880, 18: 198352, 21: 24288, 24: 48}


def test_psl23_generators_preserve_BOTH_lanes_as_one_Z6_map():
    for name, g in automorphism_generators().items():
        assert U.preserves(g), name


def test_the_literal_inverse_maps_are_NOT_automorphisms():
    inv = [0] + [pow(i, -1, 23) for i in range(1, 23)] + [INF]
    refl = [0] + [(-i) % 23 for i in range(1, 23)] + [INF]
    assert not U.preserves(ZMonomial(inv))
    assert not U.preserves(ZMonomial(refl))


def test_signs_are_needed_on_the_ternary_lane_for_S():
    S = automorphism_generators()["S"]
    unsigned = ZMonomial(S.perm)
    zero = np.zeros(24, dtype=np.int64)
    bin_ok = all(not U.syndromes(unsigned.apply(pack(r, zero)))[0].any() for r in U.G2)
    ter_ok = all(not U.syndromes(unsigned.apply(pack(zero, r)))[1].any() for r in U.G3)
    assert bin_ok and not ter_ok                          # pure permutation works for bits only


def test_automorphism_orders_include_22_24_46_via_the_ternary_double_cover():
    seen = automorphism_orders(maxlen=5)
    assert set(seen) == EXPECTED_ORDERS
    for gf3, gf2 in ((22, 11), (24, 12), (46, 23)):
        assert (gf3, gf2) in seen and gf3 == 2 * gf2          # the factor 2 lives in the trit lane only
        el = element_from_word(seen[(gf3, gf2)])
        assert el.lane_orders() == (gf2, gf3) and U.preserves(el)


def test_signed_group_has_order_12144_and_binary_quotient_6072():
    try:
        from sympy.combinatorics import Permutation, PermutationGroup
    except ImportError:
        skip("sympy not installed")

    def pts(g, signed):
        img = [0] * 48
        for i in range(24):
            j = int(g.perm[i])
            plus = (not signed) or g.units[i] == 1
            img[i] = j if plus else 24 + j
            img[24 + i] = 24 + j if plus else j
        return Permutation(img)
    gens = automorphism_generators()
    assert PermutationGroup([pts(g, True) for g in gens.values()]).order() == 12144
    assert PermutationGroup([Permutation([int(x) for x in g.perm]) for g in gens.values()]).order() == 6072


def test_encode_decode_corrects_3_bit_errors_and_4_trit_errors_at_independent_positions():
    rnd = np.random.default_rng(11)
    for _ in range(25):
        m2, m3 = rnd.integers(0, 2, 12), rnd.integers(0, 3, 12)
        w = U.encode(m2, m3)
        assert U.is_codeword(w)
        bad = w.copy()
        pos = rnd.permutation(24)
        for p in pos[:3]:
            bad[p] = (bad[p] + 3) % 6                     # flips the GF(2) part only
        for p in pos[3:7]:
            bad[p] = (bad[p] + 4 * rnd.integers(1, 3)) % 6   # changes the GF(3) part only
        r = U.decode(bad)
        assert r.ok and np.array_equal(r.m2, m2) and np.array_equal(r.m3, m3) and np.array_equal(r.word, w)


def test_beyond_radius_is_reported_not_silently_accepted():
    rnd = np.random.default_rng(12)
    for _ in range(25):
        w = U.encode(rnd.integers(0, 2, 12), rnd.integers(0, 3, 12))
        bad = w.copy()
        for p in rnd.permutation(24)[:4]:
            bad[p] = (bad[p] + 3) % 6                     # 4 bit errors: extended Golay detects, cannot correct
        assert not U.decode(bad).ok


def test_syndrome_tables_have_no_collisions_up_to_the_decoding_radius():
    assert len(U.table2) == 1 + 24 + 276 + 2024
    assert len(U.table3) == 1 + 48 + 1104 + 16192 + 170016
