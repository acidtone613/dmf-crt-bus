import itertools

import numpy as np

from crt_bus.zring import UNITS, ZMonomial, pack, unpack


def test_pack_is_a_bijection_Z6_equals_GF2_x_GF3():
    seen = {int(pack(np.array(a), np.array(b))) for a in (0, 1) for b in (0, 1, 2)}
    assert seen == set(range(6))
    for z in range(6):
        b2, b3 = unpack(np.array(z))
        assert int(pack(b2, b3)) == z


def test_ring_operations_act_on_both_fields_at_once():
    for x, y in itertools.product(range(6), repeat=2):
        (x2, x3), (y2, y3) = unpack(np.array(x)), unpack(np.array(y))
        s2, s3 = unpack(np.array((x + y) % 6))
        p2, p3 = unpack(np.array((x * y) % 6))
        assert (int(s2), int(s3)) == ((int(x2) + int(y2)) % 2, (int(x3) + int(y3)) % 3)
        assert (int(p2), int(p3)) == ((int(x2) * int(y2)) % 2, (int(x3) * int(y3)) % 3)


def test_units_are_plus_minus_one_and_sign_is_invisible_to_the_bit_lane():
    assert UNITS == (1, 5)
    b2, b3 = unpack(np.array(5))
    assert (int(b2), int(b3)) == (1, 2)          # 5 = (1 in GF(2), -1 in GF(3))
    for z in range(6):
        b2p, _ = unpack(np.array((5 * z) % 6))
        b2o, _ = unpack(np.array(z))
        assert int(b2p) == int(b2o)              # multiplying by -1 does not change the GF(2) part


def test_zmonomial_orders_and_lane_orders():
    g = ZMonomial([1, 2, 0], [1, 5, 1])          # 3-cycle with one sign flip
    assert g.order() == 6 and g.lane_orders() == (3, 6)
    assert ZMonomial([0, 1], [5, 5]).lane_orders() == (1, 2)   # global -1: order 2 on trits, 1 on bits
