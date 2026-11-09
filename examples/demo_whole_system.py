"""A guided tour of the whole-system bus.   python3 examples/demo_whole_system.py"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from crt_bus.clock_bus import ClockBus
from crt_bus.moduli_search import choose_moduli, verify_choice
from crt_bus.system import SystemBus
from crt_bus.unified_code import automorphism_orders

bus = SystemBus()
clk = bus.clock
print("1. CLOCK LANE  channels", clk.channels, "-> range", clk.range)
t = 123_456
print("   tick", t, "-> residues", clk.encode(t), "-> phases", clk.phases(clk.encode(t)))
print("   when do binary@23, ternary@0, quadalic@0 (master clocks) align?", ClockBus.align({46: 23, 759: 0, 1081: 0}))
print("   can binary@1 and ternary@0 (masters) ever co-occur?", ClockBus.align({46: 1, 759: 0}))
print("   MI(binary orbit, ternary orbit) =", round(ClockBus.mi_bits(22, 33), 4), "bits = log2 gcd(22,33)")

print("\n2. FIELD LANES (one Z_6 symbol = one bit + one trit)")
s0 = bus.default_state(seed=42)
w = bus.wire(t, s0)
print("   wire A (24 symbols):", w.A[:12], "...")
print("   consistency report:", {k: v for k, v in bus.consistency_report(w, s0).items() if k != "lanes"})
bad = list(w.clock); bad[3] = (bad[3] + 9) % 23; w.clock = tuple(bad)          # corrupt the 23-channel
rep = bus.consistency_report(w, s0)
print("   after corrupting one clock channel: tick still recovered =", rep["t"], "| corrected =", rep["clock_was_corrupted"])

print("\n3. UNIFIED CODE  binary Golay [24,12,8] x ternary XQR23 [24,12,9] as one Z_6 code")
U = bus.unified_code
rng = np.random.default_rng(0)
m2, m3 = rng.integers(0, 2, 12), rng.integers(0, 3, 12)
cw = U.encode(m2, m3)
noisy = cw.copy()
for p in (0, 5, 9): noisy[p] = (noisy[p] + 3) % 6            # 3 bit errors
for p in (1, 2, 13, 20): noisy[p] = (noisy[p] + 4) % 6       # 4 trit errors, other positions
r = U.decode(noisy)
print("   7 symbol errors (3 on the bit lane, 4 on the trit lane) corrected:", r.ok,
      "| message recovered:", bool(np.array_equal(r.m2, m2) and np.array_equal(r.m3, m3)))
orders = automorphism_orders()
print("   symmetry orders (GF(3)-lane, GF(2)-lane):", sorted(orders))

print("\n4. MODULE SELECTION")
r = choose_moduli(285384, bound=64, must_cover=(22, 33, 24, 46, 759, 1081), redundancy=2)
print("   target clocks ->", r["moduli"], "+ redundant", r["redundant"], "| violations:",
      verify_choice(r["moduli"], 285384, must_cover=(22, 33, 24, 46, 759, 1081), redundant=r["redundant"]))
print("   10^12 range with NTT size 256 ->", choose_moduli(10 ** 12, bound=20000, root_orders=(256,))["moduli"])
