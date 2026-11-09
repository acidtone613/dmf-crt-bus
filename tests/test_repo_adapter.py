import numpy as np

from _compat import skip
from crt_bus.analysis import plugin_mi
from crt_bus.repo_adapter import compare_with_repo, load_repo
from crt_bus.system import SystemBus
from math import gcd, log2


def test_package_matches_the_real_dmf_repo():
    dmf = load_repo()
    if dmf is None:
        skip("DMF repo not found (set DMF_TRUTH_PATH to the 'truth' folder)")
    res = compare_with_repo(dmf)
    assert all(res.values()), res


def test_mi_law_on_the_repos_real_binary_floquet_tick():
    dmf = load_repo()
    if dmf is None:
        skip("DMF repo not found (set DMF_TRUTH_PATH to the 'truth' folder)")
    S = SystemBus()
    # real binary clock: iterate the repo's own tick on a label vector for one full joint period
    v, states = np.arange(24), []
    for _ in range(66):
        states.append(v.tobytes())
        v = dmf.floquet.apply_tick_vector(v.copy())
    ter = S.crystals["ter"]
    xs = states
    ys = [ter.powers[t % 33].tobytes() for t in range(66)]
    assert abs(plugin_mi(xs, ys) - log2(gcd(22, 33))) < 1e-9
