# Corpus audit (from the DMF `truth` extraction) — what holds and what does not

Method: every claim was recomputed from the repo's own operators/codes (`truth/dmf`) with exact arithmetic.
Scripts: `run_checks.py` (this package) reproduces every item marked [T].

## Holds
* [T] CRT/LCM arithmetic: all pairwise and triple LCMs, gcd grid and event counts; LCM(22,33,24)=264; LCM(46,759,1081)=71,346.
* [T] Quadalic crystal: P47, P23^Q, C2 have orders 47, 23, 2; semidirect relation; ⟨P47,P23^Q⟩=1081; M4=P47∘P23^Q has order 23; M=C2∘P47∘P23^Q has order **24 only with P23^Q applied first** (the other composition order gives 4); ⟨P47,P23^Q,C2⟩ = PSL(2,47), order 51,888.
* [T] Binary Floquet tick (4 factors, `floquet.py`) has order 22: one 22-cycle, fixed points 7 and 23.
* [T] RC-232 subgroup data for ⟨H_T,P23_T⟩ (order 36,501; element orders 1/3/23/69 with counts 1/1,058/12,166/23,276) is exactly C₂₃ ≀ C₃. So P23_T acts as three independent Z_23 phases, **not** on the 12 coordinates. The 805 = 5·7·23 fits a 12-point part (5+7) times a 23-clock.
* [T] 759 = LCM(23, 33): a clock-times-schedule period, not an element order. Order 33 cannot be a permutation of 12 points (needs 3+11 = 14); under the assumed realization the schedule "H_T every 11th tick" has period 33 with base operator order 11.
* [T] Orbit-level CRT Z_264 ≅ Z_8 × Z_3 × Z_11 and group-level prime-factor CRT (injective map Z_71346 → Z_46 × Z_759 × Z_1081).

## Does not hold (as stated)
* "71,346 = 270 × 264 = LCM(22,33,24)" (Time Spiral): 270×264 = 71,280; 264 does not divide 71,346 (remainder 66); LCM(22,33,24)=264. The joint period of *all six* clocks is 285,384 = 2³·3·11·23·47.
* "M_U = H_L∘P23∘P11 has order 22" (RC-207): that 3-factor composite has order **21**. The 4-factor composite with S_inv (`floquet.py`) has order 22.
* "D23 = P23∘H_L has order 46": with H_L the inversion map i→1/i the composite has order 23; no element of D₂₃ has order 46. (`engine.py` correctly uses the reflection i→−i for D₂₃.)
* "⟨P23,P11,S,H_L⟩ = 6072" (RC-146): not reproducible with the literal H_L (see README item 3). `automorphism.py` itself notes the S matrix (RC-100) is not in the extraction.
* The literal H_L, S_inv and R23: **not automorphisms** of the extended Golay code (tested; P23, P11 and the Möbius S are).
* RC-139b's listed permutation F has order 12 (cycle type 1,3,4⁵), not 22: an older operator, superseded.
* Arithmetic slip in a dossier: 46·759·1081 = 37,742,034 (not 37,736,000).

## New (constructive)
* Both length-24 crystal codes (binary Golay, ternary XQR23) are extended QR codes at p = 23 and form one Z_6-linear code (`unified_code.py`); its symmetry group is the signed 2·PSL(2,23) of order 12,144.
* Orders 22, 24, 46 occur in that group (lane orders (22,11), (24,12), (46,23)), i.e. the clock periods that PSL(2,23) cannot supply appear in the double cover.
* MI between two clocks over a full joint period = log₂ gcd; fraction of reachable joint phases = lcm/product.
* Channels (3, 8, 11, 23, 47) are the minimal prime-power channel set covering all six clock periods.

## Suggested corrections to the corpus text
1. Replace the Time Spiral line with "joint period of orbit clocks = 264; of all six clocks = 285,384; 71,346 is the LCM of the three master clocks only".
2. Replace RC-207's 3-factor composite by the 4-factor composite of `floquet.py`, or state order 21.
3. Define S (RC-100) explicitly (the Möbius map x→−1/x swapping 0 and ∞ gives PSL(2,23) = 6,072) and say whether H_L belongs to the group.
4. State that 33 and 759 are schedule/clock periods, and put the exact ternary generators in the repo.
