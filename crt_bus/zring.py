"""zring.py -- CRT applied to the FIELDS: Z_6 = GF(2) x GF(3).

Because gcd(2, 3) = 1, the ring Z_6 splits as  Z_6  ~  GF(2) x GF(3).
One Z_6 symbol therefore carries one bit AND one trit at the same time, and
ring addition / multiplication act on both at once:

    unpack(x + y) == (b2 xor b2', (b3 + b3') mod 3)
    unpack(x * y) == (b2 and b2', (b3 * b3') mod 3)

The unit group of Z_6 is {1, 5} = {+1, -1}.  Under the CRT split:
    +1  ->  (1 in GF(2), +1 in GF(3))
    -1  ->  (1 in GF(2), -1 in GF(3))      (-1 == +1 in GF(2): the sign is INVISIBLE to the bit lane)
So signed permutations (monomial maps) of a ternary code are exactly Z_6-monomial maps,
and they act on the binary lane as plain permutations.  See unified_code.py.

NOTE: Z_4 is NOT GF(4), and 2 and 4 are not coprime, so binary and "quaternary" parts
cannot be joined by CRT.  GF(4) is a field EXTENSION of GF(2) (a different mechanism).
"""
import numpy as np

E2 = 3   # idempotent: 3 = 1 (mod 2), 3 = 0 (mod 3)
E3 = 4   # idempotent: 4 = 0 (mod 2), 4 = 1 (mod 3)
UNITS = (1, 5)


def pack(b2, b3):
    """(GF(2) array, GF(3) array) -> Z_6 array of the same shape."""
    b2 = np.asarray(b2, dtype=np.int64)
    b3 = np.asarray(b3, dtype=np.int64)
    if b2.shape != b3.shape:
        raise ValueError("lanes must have the same shape (pad the shorter lane with zeros)")
    if ((b2 < 0) | (b2 > 1)).any() or ((b3 < 0) | (b3 > 2)).any():
        raise ValueError("lane values out of range: GF(2) in {0,1}, GF(3) in {0,1,2}")
    return (E2 * b2 + E3 * b3) % 6


def unpack(z):
    """Z_6 array -> (GF(2) array, GF(3) array)."""
    z = np.asarray(z, dtype=np.int64) % 6
    return z % 2, z % 3


def sign_to_unit(s):
    """+1 -> 1, -1 -> 5 (the units of Z_6)."""
    return 1 if s == 1 else 5


class ZMonomial:
    """A Z_6-monomial map: out[perm[i]] = units[i] * word[i]  (mod 6).

    `units` entries must be in {1, 5}.  Maps act on binary lanes as permutations and on
    ternary lanes as signed permutations.
    """

    def __init__(self, perm, units=None):
        self.perm = np.asarray(perm, dtype=np.int64)
        n = len(self.perm)
        self.units = np.ones(n, dtype=np.int64) if units is None else np.asarray(units, dtype=np.int64)
        if sorted(self.perm.tolist()) != list(range(n)):
            raise ValueError("perm is not a permutation")
        if not set(self.units.tolist()) <= set(UNITS):
            raise ValueError("units must be 1 or 5 (the units of Z_6)")

    def __len__(self):
        return len(self.perm)

    def apply(self, word):
        word = np.asarray(word, dtype=np.int64)
        out = np.zeros_like(word)
        out[..., self.perm] = (self.units * word) % 6
        return out

    def then(self, other):
        """The map 'apply self, then other'."""
        n = len(self)
        perm = np.empty(n, dtype=np.int64)
        units = np.empty(n, dtype=np.int64)
        for i in range(n):
            j = self.perm[i]
            perm[i] = other.perm[j]
            units[i] = (self.units[i] * other.units[j]) % 6
        return ZMonomial(perm, units)

    def identity(n):  # noqa: N805  (static helper)
        return ZMonomial(np.arange(n))
    identity = staticmethod(identity)

    def order(self, limit=10_000):
        cur, ident = self, ZMonomial.identity(len(self))
        k = 1
        while not (np.array_equal(cur.perm, ident.perm) and np.array_equal(cur.units, ident.units)):
            cur = cur.then(self)
            k += 1
            if k > limit:
                raise RuntimeError("order exceeds limit")
        return k

    def lane_orders(self):
        """(order on the GF(2) lane, order on the GF(3) lane)."""
        binary_part = ZMonomial(self.perm)                 # drop signs: what GF(2) sees
        return binary_part.order(), self.order()
