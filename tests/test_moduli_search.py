import random
from itertools import product as iproduct

from crt_bus.moduli_search import (_smallest_prime, brute_force_choose, candidate_options, choose_moduli,
                                   verify_choice)

DMF_PERIODS = (22, 33, 24, 46, 759, 1081)


def _bf_key(R, bound, den, objective):
    opts, forced = candidate_options(bound, den)
    best = None
    for ch in iproduct(*[[None] + opts[p] for p in sorted(opts)]):
        c = [q for q in ch if q]
        if not c:
            continue
        pr = 1
        for q in c:
            pr *= q
        if pr < R:
            continue
        k = {"channels": (len(c), pr), "bits": (pr, len(c)), "width": (max(c), len(c), pr)}[objective]
        if best is None or k < best[0]:
            best = (k, tuple(sorted(c)))
    return best


def test_search_matches_exhaustive_enumeration_on_random_small_specs():
    rnd = random.Random(10)
    checked = 0
    for _ in range(80):
        R, bound = rnd.randrange(2, 3000), rnd.choice([8, 12, 16, 20, 24])
        den = tuple(rnd.sample([2, 3, 5, 7], rnd.choice([0, 0, 1])))
        for obj in ("channels", "bits", "width"):
            try:
                got = choose_moduli(R, bound=bound, denominators=den, objective=obj)["moduli"]
            except ValueError:
                got = None
            best = _bf_key(R, bound, den, obj)
            assert got == (best[1] if best else None), (R, bound, den, obj)
            checked += 1
    assert checked == 240


def test_brute_force_reference_agrees_for_channels_objective():
    for R in (10, 100, 1000, 2500):
        a = choose_moduli(R, bound=20)["moduli"]
        assert a == brute_force_choose(R, 20, objective="channels")[0]


def test_dmf_clock_requirement_reproduces_the_hand_built_bus():
    r = choose_moduli(285384, bound=64, must_cover=DMF_PERIODS, redundancy=2)
    assert r["moduli"] == (3, 8, 11, 23, 47) and r["redundant"] == (53, 59) and r["product"] == 285384
    assert verify_choice(r["moduli"], 285384, must_cover=DMF_PERIODS, redundant=r["redundant"]) == []


def test_denominators_are_made_whole_by_avoiding_their_primes():
    r = choose_moduli(10 ** 6, bound=64, denominators=(2, 3))     # so that 1/2 and 1/3 exist mod every channel
    assert all(m % 2 and m % 3 for m in r["moduli"])
    assert verify_choice(r["moduli"], 10 ** 6, denominators=(2, 3)) == []
    assert verify_choice((6, 7), 42, denominators=(2,)) != []     # mod 6 there is no 1/2


def test_ntt_channels_need_primes_congruent_1_mod_n():
    r = choose_moduli(10 ** 12, bound=20000, root_orders=(256,))
    assert all((m - 1) % 256 == 0 for m in r["moduli"])
    assert verify_choice(r["moduli"], 10 ** 12, root_orders=(256,)) == []
    assert verify_choice((3329, 7681), 1000, root_orders=(256,)) == []   # real lattice-crypto moduli qualify
    assert verify_choice((3329, 7681), 1000, root_orders=(512,)) != []   # 3329 has no 512th root (3328 = 2^8 * 13)


def test_verifier_catches_each_kind_of_violation():
    assert any("coprime" in v for v in verify_choice((6, 9), 10))
    assert any("product" in v for v in verify_choice((3, 5), 100))
    assert any("readable" in v for v in verify_choice((3, 5), 10, must_cover=(22,)))
    assert any("MDS" in v for v in verify_choice((3, 5), 10, redundant=(7, 2)))


def test_width_objective_trades_channels_for_narrowness():
    narrow = choose_moduli(285384, bound=64, objective="width")
    few = choose_moduli(285384, bound=64, objective="channels")
    assert narrow["max"] <= few["max"] and few["count"] <= narrow["count"]
    assert _smallest_prime(49) == 7
