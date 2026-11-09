# DMF CRT Bus

A small, exact, fully-tested library that puts the **whole DMF system** on one Chinese-Remainder-Theorem bus:
("DMF" is the parent research corpus this grew out of; this library is self-contained and depends on none of it.)

| Lane | What it carries | CRT fact it uses |
|---|---|---|
| **Clock lane** | the tick counter, as residues mod (8, 3, 11, 23, 47) plus redundant (53, 59) | Z_285384 ≅ Z_8 × Z_3 × Z_11 × Z_23 × Z_47; generalized CRT for clocks that share a prime |
| **Field lanes** | GF(2) state (binary + quadalic crystals) and GF(3) state (ternary crystal) | **Z_6 ≅ GF(2) × GF(3)**: one ring symbol = one bit **and** one trit |
| **Code lane** | the binary Golay [24,12,8] and ternary QR [24,12,9] codes as **one** Z_6-linear code | both are extended QR codes at the prime 23, so they share the symmetry group PSL(2,23) (signed, on the trit side) |

Everything is integer arithmetic and every claim below is asserted by a test. Nothing here is cryptanalysis, and nothing
depends on the removed physics-layer claims of the DMF corpus.

## Quick start
```
pip install numpy            # sympy / pytest optional
python3 run_checks.py        # 50 tests (48 run standalone, 2 need the repo), a few seconds
python3 examples/demo_whole_system.py
```
To cross-check against the real repo, point at its `truth` folder:
`DMF_TRUTH_PATH=/path/to/truth python3 run_checks.py` (without it, the two repo tests are skipped; the package bundles its own copies of the operators and codes).

## Package map (`crt_bus/`)
| Module | Purpose | Key API |
|---|---|---|
| `crt.py` | exact number theory | `crt_merge`, `crt_solve` (non-coprime OK), `pairwise_consistent`, `decode_coprime` |
| `clock_bus.py` | tick counter as residues | `ClockBus.encode / increment / phase / phases / align / correct / mi_bits / reachable_fraction` |
| `zring.py` | Z_6 = GF(2)×GF(3) packing | `pack`, `unpack`, `ZMonomial` (signed permutations; `.order()`, `.lane_orders()`) |
| `unified_code.py` | Golay ⊗ XQR23 as one Z_6 code | `UnifiedCode24.encode / decode / preserves`, `automorphism_generators`, `automorphism_orders` |
| `system.py` | whole-system bus | `SystemBus.state_at / wire / read_wire / consistency_report / joint_tick_mod` |
| `moduli_search.py` | fewest/narrowest channels for a computation | `choose_moduli`, `verify_choice` (exact verifier), `brute_force_choose` |
| `repo_adapter.py` | cross-check vs the real `truth/dmf` | `load_repo`, `compare_with_repo` |
| `analysis.py` | plug-in mutual information | `plugin_mi` |

## What is verified (test → fact)
* **Six clocks from five residues.** `phase(residues, p) == t % p` for periods 22, 33, 24, 46, 759, 1081; the joint period is 285,384 = 2³·3·11·23·47. Increment is carry-free.
* **Alignment queries.** binary@23 + ternary@0 + quadalic@0 (master clocks) first coincide at t = 35,673 and then every 71,346; binary@1 with ternary@0 can never co-occur (needs 1 ≡ 0 mod 23). Pairwise consistency ⇔ solvable.
* **Desync detector.** Only 264 of the 17,424 orbit-phase triples (1 in 66) are reachable.
* **Mutual-information law.** Over a full joint period, two clocks have MI = **log₂ gcd** of their periods (22,33 → 3.459 bits; coprime → 0). Verified on the actual crystal orbits and on the repo's real binary Floquet tick; coarse observables carry less.
* **Fault tolerance.** With redundant channels 53 and 59, *every* single-channel corruption of *every* residue value is corrected (MDS: any 5 of the 7 moduli cover the range); out-of-range values are flagged, not wrapped.
* **Whole-system wire.** Binary, ternary and quadalic states share Z_6 symbols; corrupted clock channels are corrected; a flipped bit is localized to its lane; non-zero padding is detected.
* **Unified Z_6 code.** Golay [24,12,8] (weights 0/8/12/16/24: 1/759/2576/759/1) and XQR23 [24,12,9] (weights 9…24: 4048/61824/242880/198352/24288/48) are both self-dual and both preserved by P23, P11 and S: x → −1/x (0↔∞). S is a plain permutation on bits but a **signed** permutation on trits (ε₀ = −1, ε∞ = +1, εᵢ = Legendre(i)); signs are units of Z_6, and −1 = +1 in GF(2). The signed group has order **12,144 = 2·|PSL(2,23)|** and contains the global −1; the binary lane sees only PSL(2,23) (6,072).
* **Joint decoding.** 3 bit-lane errors and 4 trit-lane errors at *independent* positions (7 symbol errors) are corrected together; 4 bit-lane errors are reported as uncorrectable, never accepted silently.
* **Clock orders 22, 24, 46 from the double cover.** PSL(2,23) has no elements of these orders, but the signed group acting on the Z_6 code does: (GF(3)-lane order, GF(2)-lane order) = (22, 11), (24, 12), (46, 23). The factor 2 lives entirely in the trit lane. *This shows such elements exist as code automorphisms; it does not by itself prove the corpus clocks are these elements.*
* **Module selection.** The exhaustive search reproduces the hand-built bus from the requirement "all six DMF periods readable": (3, 8, 11, 23, 47) + redundant (53, 59). It agrees with brute force on 240 random small specs; it also finds NTT-friendly channels (primes ≡ 1 mod 256).

## Assumptions and open items (read before relying on this)
1. **Ternary generators are assumed.** The exact P_mult, P_23^T and H_T schedule are not in the extracted corpus. `system.ternary_step_perm` uses a natural realization (12-cycle; x→2x on 11 points; block 3-cycles; H_T every 11th tick) that reproduces the corpus numbers (base operator of order 11, schedule period 33). Swap in the real generators; nothing else changes. Crystals carry an `assumed` flag.
2. **The repo's literal binary Floquet tick is *not* a Golay automorphism** (H_L: i→1/i and S_inv fix 0 and ∞; they do not preserve the code). Its period 22 is real (single 22-cycle, fixed points 7 and 23) but the state does not stay inside the code. The code-preserving alternative is the signed group in `unified_code.py`.
3. **The "6072" claim** (⟨P23,P11,S,H_L⟩) is not reproducible with H_L as defined: ⟨P23,P11,S_Möbius⟩ is exactly 6,072 (= PSL(2,23)); adding the literal H_L gives A₂₄.
4. Z_4 is **not** GF(4), and 2 and 4 are not coprime, so CRT cannot join a "quaternary" part to the binary part; GF(4) is an extension of GF(2). F₄ in DMF is the hexacode only.
5. The unified Z_6 code lives on 24 coordinates (binary Golay + ternary XQR23). The ternary *Golay* [12,6,6] crystal is a different, length-12 object and is carried on the ternary state lane, not inside this code.

## Reuse ideas
* **Reward/verifier for learned search.** `moduli_search.verify_choice` is an exact checker and `brute_force_choose` gives ground truth at small sizes, so a learned policy that proposes channel sets can be scored with no ambiguity.
* **Environment generator.** Any pairwise-coprime clock set, any pair of fields (p, q) with Z_pq, any list of denominators or root-of-unity orders gives a new instance with known optimum.
* **Template.** To bus a new periodic subsystem: factor its period into prime powers, add channels for new primes, and `ClockBus` reads it out; `crt_solve` answers alignment and desync questions across clocks that share primes.
* **Z_6 trick.** Any two coprime-characteristic lanes (GF(2)/GF(3), or GF(p)/GF(q)) can share one ring symbol; field operations never mix.

See `FINDINGS.md` for the corpus audit that motivated this.
