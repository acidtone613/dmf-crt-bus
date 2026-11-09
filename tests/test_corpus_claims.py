"""Audit of specific DMF corpus claims (see FINDINGS.md).  Each test states what the corpus says and what is true."""
from itertools import combinations
from math import lcm

import numpy as np

from _compat import skip
from crt_bus.system import binary_tick_perm, quad_tick_perm, ternary_step_perm


def compose(*fs):                      # f1 o f2 o ... : apply the LAST one first (maps point -> point)
    out = list(range(len(fs[0])))
    for f in reversed(fs):
        out = [f[x] for x in out]
    return out


def order(p):
    p, ident, k = list(p), list(range(len(p))), 1
    cur = list(p)
    while cur != ident:
        cur = [p[x] for x in cur]
        k += 1
    return k


P23 = [(i + 1) % 23 if i < 23 else 23 for i in range(24)]
P11 = [(2 * i) % 23 if i < 23 else 23 for i in range(24)]
HL = [0] + [pow(i, -1, 23) for i in range(1, 23)] + [23]
SINV = [0] + [pow((-i) % 23, -1, 23) for i in range(1, 23)] + [23]
R23 = [0] + [(-i) % 23 for i in range(1, 23)] + [23]
SMOB = [23] + [(-pow(i, -1, 23)) % 23 for i in range(1, 23)] + [0]


def test_joint_period_arithmetic():
    assert lcm(22, 33, 24) == 264 and lcm(46, 759, 1081) == 71346 and lcm(22, 33, 24, 46, 759, 1081) == 285384
    assert 270 * 264 == 71280 != 71346 and 71346 % 264 == 66        # the 'Time Spiral' line is wrong
    assert 46 * 759 * 1081 == 37742034


def test_rc207_three_factor_composite_has_order_21_not_22_but_the_repo_four_factor_tick_has_22():
    assert order(compose(HL, P23, P11)) == 21
    assert order(compose(SINV, HL, P11, P23)) == 22
    tick = binary_tick_perm()
    cycles = {len(c) for c in _cycles(tick)}
    assert order(list(tick)) == 22 and sorted(len(c) for c in _cycles(tick)) == [1, 1, 22]
    assert sorted(i for i in range(24) if tick[i] == i) == [7, 23]


def _cycles(p):
    seen, out = set(), []
    for i in range(len(p)):
        if i not in seen:
            c, j = [], i
            while j not in seen:
                seen.add(j)
                c.append(j)
                j = p[j]
            out.append(c)
    return out


def test_d23_claim_reflection_vs_inversion():
    assert order(compose(P23, R23)) == 2        # true dihedral reflection (what engine.py uses)
    assert order(compose(P23, HL)) == 23        # the inversion map is NOT a reflection; no element of order 46 exists


def test_quadalic_composite_orders_depend_on_composition_order():
    p47 = [(i + 1) % 47 for i in range(47)] + [47]
    p23q = [(2 * i) % 47 for i in range(47)] + [47]
    c2 = [47] + [(-pow(i, -1, 47)) % 47 for i in range(1, 47)] + [0]
    assert order(compose(p47, p23q)) == 23
    assert order(compose(c2, p47, p23q)) == 24          # P23^Q applied first (corpus RC-267)
    assert order(compose(p23q, p47, c2)) == 4           # the other composition order
    scatter = compose(c2, p47, p23q)
    gather = [0] * 48
    for i, j in enumerate(scatter):
        gather[j] = i
    assert list(quad_tick_perm()) == gather


def test_order_33_759_805_cannot_be_permutations_of_12_points():
    def partitions(n, m=None):
        m = n if m is None else m
        if n == 0:
            yield []
            return
        for k in range(min(n, m), 0, -1):
            for rest in partitions(n - k, k):
                yield [k] + rest
    orders = {lcm(*p) for p in partitions(12)}
    assert max(orders) == 60
    assert not ({23, 33, 69, 759, 805} & orders)
    assert 5 + 7 == 12                                    # 805 = 5*7*23: 5- and 7-cycles fit, the 23 is a separate clock


def test_ternary_assumed_schedule_matches_corpus_numbers():
    pcyc = [(i + 1) % 12 for i in range(12)]
    pmult = [(2 * i) % 11 for i in range(11)] + [11]
    base = [pmult[pcyc[i]] for i in range(12)]
    assert order(base) == 11
    c = np.arange(12)
    for t in range(1, 34):
        c = c[ternary_step_perm(t)]
    assert np.array_equal(c, np.arange(12))
    c = np.arange(12)
    for t in range(1, 33):
        c = c[ternary_step_perm(t)]
        assert not np.array_equal(c, np.arange(12))      # first return is exactly 33


def test_rc232_subgroup_data_is_the_wreath_product_C23_wr_C3():
    n, counts = 23, {1: 0, 3: 0, 23: 0, 69: 0}
    for k in range(3):
        for v0 in range(n):
            for v1 in range(n):
                for v2 in range(n):
                    if k == 0:
                        o = 1 if (v0, v1, v2) == (0, 0, 0) else 23
                    else:
                        o = 3 if (v0 + v1 + v2) % n == 0 else 69
                    counts[o] += 1
    assert 3 * 23 ** 3 == 36501 and counts[1] == 1 and counts[3] == 1058 and counts[23] == 12166 and counts[69] == 23276


def test_6072_needs_the_mobius_S_and_the_literal_H_L_breaks_it():
    try:
        from sympy.combinatorics import Permutation, PermutationGroup
    except ImportError:
        skip("sympy not installed")
    G = lambda *gs: PermutationGroup([Permutation(g) for g in gs]).order()
    assert G(P23, P11, SMOB) == 6072
    assert G(P23, P11, SMOB, HL) > 10 ** 20                # A_24, not 6072
    assert G(P23, P11, SINV, HL) > 10 ** 20                # the literal generators generate S_23
    p47 = [(i + 1) % 47 for i in range(47)] + [47]
    p23q = [(2 * i) % 47 for i in range(47)] + [47]
    c2 = [47] + [(-pow(i, -1, 47)) % 47 for i in range(1, 47)] + [0]
    assert G(p47, p23q, c2) == 51888
    assert G(p47, p23q) == 1081
