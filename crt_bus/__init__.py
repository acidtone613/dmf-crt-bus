"""crt_bus -- a Chinese-Remainder-Theorem bus for a whole system of periodic processes.

Lanes
  clock lane   : tick counter as residues (8, 3, 11, 23, 47) + redundant (53, 59)
  field lanes  : GF(2) and GF(3) state packed into one ring, Z_6 = GF(2) x GF(3)
  code lane    : unified Z_6 code = binary Golay [24,12,8] (x) ternary QR [24,12,9]
See README.md for the map of the package.
"""
__version__ = "1.0.0"
