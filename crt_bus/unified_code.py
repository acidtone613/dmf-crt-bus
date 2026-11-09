"""unified_code.py -- the GF(2) and GF(3) codes here as ONE Z_6-linear code.

Both crystals' length-24 codes are extended quadratic-residue codes at the SAME prime 23:

    binary Golay    [24,12,8]_2   (2 is a QR mod 23)
    ternary QR      [24,12,9]_3   (3 is a QR mod 23)   -- 'XQR23' in this package

Coordinates are the projective line P^1(F_23) = {0..22, inf}; inf is index 23.
By CRT (Z_6 = GF(2) x GF(3)), the pair is one Z_6-linear code  C6 = CRT(Golay, XQR23)
of length 24:  a codeword is pack(c2, c3) with c2 in Golay, c3 in XQR23.

SYMMETRY (verified by tests): both codes are preserved by
    P23 : i -> i+1 (mod 23)           P11 : i -> 2i (mod 23)
    S   : i -> -1/i  (0 <-> inf)      -- a pure permutation for the binary code, but a
                                         SIGNED permutation for the ternary code, with signs
                                         eps_0 = -1, eps_inf = +1, eps_i = Legendre(i).
Signs are units of Z_6 (+1 -> 1, -1 -> 5), and -1 == +1 in GF(2), so the same Z_6 map acts as
a permutation on the bit lane and a signed permutation on the trit lane.  The group generated
has order 12,144 = 2 * |PSL(2,23)|; it contains the global sign -1 (S^2 = -1).  The binary lane
only sees the quotient PSL(2,23) (order 6,072).

CONSEQUENCE FOR CLOCKS: PSL(2,23) has no element of order 22, 24 or 46, but the double cover
acting on C6 does (GF(3)-lane orders 22, 24, 46 with GF(2)-lane orders 11, 12, 23).  The factor
2 lives entirely in the trit lane.  See `automorphism_orders`.

Decoding: because both codes are self-dual, the parity-check matrix equals the generator
matrix.  Syndrome tables correct up to 3 errors on the bit lane and up to 4 on the trit lane
(independently); an error pattern beyond that is reported as uncorrectable, never silently
accepted.
"""
import itertools
from dataclasses import dataclass

import numpy as np

from .zring import ZMonomial, pack, sign_to_unit, unpack

N = 24
INF = 23
GASC = [1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1]      # g(x) = x^11+x^10+x^6+x^5+x^4+x^2+1, ascending


# ------------------------------------------------------------------ mod-p linear algebra
def _rref_mod(M, p):
    M = np.array(M, dtype=np.int64) % p
    rows, cols = M.shape
    piv, r = [], 0
    for c in range(cols):
        if r >= rows:
            break
        pr = next((i for i in range(r, rows) if M[i, c] % p), None)
        if pr is None:
            continue
        M[[r, pr]] = M[[pr, r]]
        M[r] = (M[r] * pow(int(M[r, c]), -1, p)) % p
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        piv.append(c)
        r += 1
    return M, piv


def _inv_mod(A, p):
    n = A.shape[0]
    aug = np.hstack([A % p, np.eye(n, dtype=np.int64)])
    R, piv = _rref_mod(aug, p)
    assert piv[:n] == list(range(n)), "matrix is singular mod p"
    return R[:, n:]


# ------------------------------------------------------------------ the two codes
def build_golay24():
    """12x24 generator matrix of the extended binary Golay code."""
    seed = GASC + [0] * 11
    G23 = np.array([[seed[(j - i) % 23] for j in range(23)] for i in range(12)], dtype=np.int64)
    return np.hstack([G23, (G23.sum(axis=1) % 2).reshape(-1, 1)])


def build_xqr23():
    """12x24 generator matrix of the extended ternary QR code [24,12,9]_3."""
    qr = sorted({(x * x) % 23 for x in range(1, 23)})
    q = [0] * 23
    for s in qr:
        q[s] = 1
    rows = []
    for i in range(12):
        row = [0] * 23
        for j, c in enumerate(q):
            row[(j + i) % 23] = (row[(j + i) % 23] + c) % 3
        rows.append(row)
    G23 = np.array(rows, dtype=np.int64)
    return np.hstack([G23, (-G23.sum(axis=1) % 3).reshape(-1, 1)])


def legendre23(i):
    i %= 23
    return 0 if i == 0 else (1 if pow(i, 11, 23) == 1 else -1)


def automorphism_generators():
    """P23, P11, S as Z_6-monomial maps on length-24 words."""
    p23 = [(i + 1) % 23 if i < 23 else INF for i in range(N)]
    p11 = [(2 * i) % 23 if i < 23 else INF for i in range(N)]
    s = [INF if i == 0 else (0 if i == INF else (-pow(i, -1, 23)) % 23) for i in range(N)]
    eps = [-1] + [legendre23(i) for i in range(1, 23)] + [1]
    return {
        "P23": ZMonomial(p23),
        "P11": ZMonomial(p11),
        "S": ZMonomial(s, [sign_to_unit(e) for e in eps]),
    }


def _syndrome_table(H, p, tmax):
    """syndrome-int -> ((positions), (values)) for all error patterns of weight <= tmax."""
    k, n = H.shape
    wts = p ** np.arange(k, dtype=np.int64)
    table = {0: ((), ())}
    for w in range(1, tmax + 1):
        combos = np.array(list(itertools.combinations(range(n), w)), dtype=np.int64)
        vals = np.array(list(itertools.product(range(1, p), repeat=w)), dtype=np.int64)
        Hc = H.T[combos]                                   # (N, w, k)
        S = np.einsum("nwk,vw->nvk", Hc, vals) % p         # (N, V, k)
        keys = S @ wts                                     # (N, V)
        for ni in range(len(combos)):
            pos = tuple(int(x) for x in combos[ni])
            for vi in range(len(vals)):
                key = int(keys[ni, vi])
                if key in table:
                    raise AssertionError("syndrome collision: radius exceeds (d-1)/2")
                table[key] = (pos, tuple(int(x) for x in vals[vi]))
    return table


@dataclass
class DecodeResult:
    word: object            # corrected Z_6 word, or None if uncorrectable
    m2: object              # 12 message bits (or None)
    m3: object              # 12 message trits (or None)
    errors2: list           # [(position, value)] corrected on the bit lane
    errors3: list           # [(position, value)] corrected on the trit lane
    ok: bool


class UnifiedCode24:
    T2, T3 = 3, 4           # correction radii: floor((d-1)/2) for d = 8 and d = 9

    def __init__(self):
        self.G2 = build_golay24()
        self.G3 = build_xqr23()
        self._tab2 = None
        self._tab3 = None
        self._w2 = 2 ** np.arange(12, dtype=np.int64)
        self._w3 = 3 ** np.arange(12, dtype=np.int64)
        self._info2, self._Ti2 = self._info(self.G2, 2)
        self._info3, self._Ti3 = self._info(self.G3, 3)

    @staticmethod
    def _info(G, p):
        _, piv = _rref_mod(G, p)
        return piv, _inv_mod(G[:, piv], p)

    # tables are built on first use (the ternary one has ~187k entries)
    @property
    def table2(self):
        if self._tab2 is None:
            self._tab2 = _syndrome_table(self.G2, 2, self.T2)
        return self._tab2

    @property
    def table3(self):
        if self._tab3 is None:
            self._tab3 = _syndrome_table(self.G3, 3, self.T3)
        return self._tab3

    # ------------------------------------------------------------ encode / test
    def encode(self, m2, m3):
        c2 = (np.asarray(m2, dtype=np.int64) @ self.G2) % 2
        c3 = (np.asarray(m3, dtype=np.int64) @ self.G3) % 3
        return pack(c2, c3)

    def syndromes(self, word):
        r2, r3 = unpack(word)
        return (self.G2 @ r2) % 2, (self.G3 @ r3) % 3     # self-dual => H == G

    def is_codeword(self, word):
        s2, s3 = self.syndromes(word)
        return not s2.any() and not s3.any()

    # ------------------------------------------------------------ decode
    def decode(self, word):
        r2, r3 = unpack(word)
        s2, s3 = self.syndromes(word)
        e2 = self.table2.get(int(s2 @ self._w2))
        e3 = self.table3.get(int(s3 @ self._w3))
        if e2 is None or e3 is None:
            return DecodeResult(None, None, None, [], [], False)
        c2, c3 = r2.copy(), r3.copy()
        for pos, val in zip(*e2):
            c2[pos] = (c2[pos] - val) % 2
        for pos, val in zip(*e3):
            c3[pos] = (c3[pos] - val) % 3
        m2 = (c2[self._info2] @ self._Ti2) % 2
        m3 = (c3[self._info3] @ self._Ti3) % 3
        return DecodeResult(pack(c2, c3), m2, m3, list(zip(*e2)), list(zip(*e3)), True)

    # ------------------------------------------------------------ symmetry
    def preserves(self, g):
        """Does the Z_6-monomial map g map the code to itself? (checked on generator rows)."""
        zero = np.zeros(N, dtype=np.int64)
        for row in self.G2:
            img = g.apply(pack(row, zero))
            s2, _ = self.syndromes(img)
            if s2.any() or (img % 3).any():
                return False
        for row in self.G3:
            img = g.apply(pack(zero, row))
            _, s3 = self.syndromes(img)
            if s3.any() or (img % 2).any():
                return False
        return True


def automorphism_orders(maxlen=5):
    """{(order on GF(3) lane = order on Z_6 word, order on GF(2) lane): example word}.

    Words are products of the generators 'a'=P23, 'b'=P11, 's'=S (applied left to right).
    """
    gens = automorphism_generators()
    letters = {"a": gens["P23"], "b": gens["P11"], "s": gens["S"]}
    seen = {}
    for L in range(1, maxlen + 1):
        for w in itertools.product("abs", repeat=L):
            el = letters[w[0]]
            for ch in w[1:]:
                el = el.then(letters[ch])
            key = el.lane_orders()[::-1]          # (GF(3) lane, GF(2) lane)
            seen.setdefault(key, "".join(w))
    return seen


def element_from_word(word):
    gens = automorphism_generators()
    letters = {"a": gens["P23"], "b": gens["P11"], "s": gens["S"]}
    el = letters[word[0]]
    for ch in word[1:]:
        el = el.then(letters[ch])
    return el
