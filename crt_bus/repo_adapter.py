"""repo_adapter.py -- optional cross-check against an external reference tree.

If a reference tree providing the `dmf` package (the "truth/dmf" modules) is
available (set DMF_TRUTH_PATH to that folder, or place it next to this
package), `compare_with_repo()` verifies that this package's operators and
codes are IDENTICAL to the reference's. If it is not found, everything else
in the package still works (it bundles its own copies).
"""
import importlib
import os
import sys

import numpy as np

from .system import binary_tick_perm
from .unified_code import build_golay24, build_xqr23


def find_truth_path(hint=None):
    here = os.path.dirname(os.path.abspath(__file__))
    cands = [hint, os.environ.get("DMF_TRUTH_PATH"),
             os.path.join(here, "..", "..", "truth"), os.path.join(here, "..", "truth"),
             os.path.join(os.getcwd(), "truth"), os.path.join(os.getcwd(), "..", "truth")]
    for c in cands:
        if c and os.path.isfile(os.path.join(c, "dmf", "__init__.py")):
            return os.path.abspath(c)
    return None


def load_repo(hint=None):
    """Return the imported `dmf` package or None."""
    path = find_truth_path(hint)
    if path is None:
        return None
    if path not in sys.path:
        sys.path.insert(0, path)
    try:
        return importlib.import_module("dmf")
    except Exception:
        return None


def compare_with_repo(dmf):
    out = {}
    out["golay_G24_equal"] = bool(np.array_equal(dmf.golay.G24, build_golay24()))
    _, _, G, _ = dmf.ternary.xqr23_code()
    mine = build_xqr23()
    out["xqr23_generator_equal"] = bool(np.array_equal(mine[:, :23], G)
                                        and np.array_equal(mine[:, 23], (-G.sum(axis=1)) % 3))
    ident = np.arange(24)
    out["binary_tick_perm_equal"] = bool(np.array_equal(dmf.floquet.apply_tick_vector(ident.copy()),
                                                        binary_tick_perm()))
    out["floquet_period_is_22"] = dmf.floquet.FLOQUET_PERIOD == 22
    return out
