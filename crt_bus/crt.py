"""crt.py -- number-theoretic core of the CRT bus.

Everything else in this package is built from three facts:

1. Chinese Remainder Theorem (coprime moduli):  Z_(m1*m2*...) is isomorphic to
   Z_m1 x Z_m2 x ...   Addition and multiplication act on each factor independently,
   with no carries between factors.
2. Generalized CRT (NON-coprime moduli): the system  t = a1 (mod n1), t = a2 (mod n2)
   is solvable  iff  a1 = a2 (mod gcd(n1, n2)).  Then the solution is unique mod lcm(n1, n2).
   This is the right tool when clocks share a prime (the bundled clocks all share 23).
3. For several congruences, pairwise consistency is necessary AND sufficient.

All functions are pure integer arithmetic (exact, no floating point).
"""
from math import gcd


def egcd(a, b):
    """Return (g, x, y) with a*x + b*y == g == gcd(a, b). Iterative (no recursion limit)."""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def lcm(*xs):
    out = 1
    for x in xs:
        out = out // gcd(out, x) * x
    return out


def crt_merge(a1, n1, a2, n2):
    """Merge t = a1 (mod n1) and t = a2 (mod n2) for ANY positive n1, n2.

    Returns (t0, L) with L = lcm(n1, n2) and 0 <= t0 < L, or None if inconsistent.
    """
    g, p, _ = egcd(n1, n2)
    if (a2 - a1) % g:
        return None
    L = n1 // g * n2
    k = (((a2 - a1) // g) * p) % (n2 // g)
    return (a1 + n1 * k) % L, L


def crt_solve(pairs):
    """Solve a list of (residue, modulus) pairs. Returns (t0, lcm) or None."""
    a, n = 0, 1
    for ai, ni in pairs:
        r = crt_merge(a, n, ai % ni, ni)
        if r is None:
            return None
        a, n = r
    return a, n


def pairwise_consistent(pairs):
    """True iff every two congruences agree modulo the gcd of their moduli (== solvable)."""
    pairs = list(pairs)
    return all((a - b) % gcd(m, n) == 0
               for i, (a, m) in enumerate(pairs) for (b, n) in pairs[i + 1:])


def decode_coprime(residues, moduli):
    """Reconstruct x in [0, prod(moduli)) from residues for PAIRWISE COPRIME moduli."""
    M = 1
    for m in moduli:
        M *= m
    x = 0
    for r, m in zip(residues, moduli):
        Mi = M // m
        x += r * Mi * pow(Mi, -1, m)
    return x % M


def pairwise_coprime(moduli):
    moduli = list(moduli)
    return all(gcd(a, b) == 1 for i, a in enumerate(moduli) for b in moduli[i + 1:])


def prime_power_factors(n):
    """[2**a, 3**b, ...] -- the prime-power parts of n (n >= 1)."""
    out, p = [], 2
    while n > 1:
        if n % p == 0:
            q = 1
            while n % p == 0:
                n //= p
                q *= p
            out.append(q)
        p += 1
    return out
