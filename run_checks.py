"""run_checks.py -- run every test with or without pytest.   python3 run_checks.py"""
import importlib
import os
import sys
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tests"))


def main():
    passed = skipped = failed = 0
    for fn in sorted(os.listdir(os.path.join(HERE, "tests"))):
        if not (fn.startswith("test_") and fn.endswith(".py")):
            continue
        mod = importlib.import_module(fn[:-3])
        for name in sorted(n for n in dir(mod) if n.startswith("test_")):
            t0 = time.time()
            try:
                getattr(mod, name)()
                passed += 1
                status = "PASS"
            except unittest.SkipTest as e:
                skipped += 1
                status = f"SKIP ({e})"
            except Exception as e:                                  # noqa: BLE001
                failed += 1
                status = f"FAIL: {type(e).__name__}: {e}"
            print(f"{status:8s} {fn[:-3]}::{name}  [{time.time()-t0:.1f}s]")
    print(f"\n{passed} passed, {skipped} skipped, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
