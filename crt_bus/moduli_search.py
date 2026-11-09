"""moduli_search.py -- choose the FEWEST / NARROWEST CRT channels that keep a computation whole.

Problem (exact, checkable):  given a computation made of +, -, * (and optionally division by known
constants and roots of unity) whose final result lies in [0, R), find pairwise-coprime moduli with

  1. product >= R                       (the result fits; intermediates may wrap -- ring ops commute with reduction)
  2. every denominator d invertible      (gcd(d, m) = 1 for every channel m)  -> fractions like 1/2 are 'whole'
  3. required roots of unity exist       (for each n in root_orders: n | m-1, m prime)  -> NTT-friendly channels
  4. every period in must_cover readable (each prime-power part q of the period divides some channel)
  5. optional redundancy r               (r extra prime channels larger than all others; MDS => corrects floor(r/2)... see ClockBus)

Objectives:  "channels" (fewest channels, then smallest product)
             "width"    (smallest largest-channel, then fewest channels)   -- hardware width
             "bits"     (smallest product, then fewest channels)

`verify_choice` is the exact verifier; `choose_moduli` is an exhaustive branch-and-bound search whose
optimum is cross-checked against `brute_force_choose` in the tests.  A learned search policy can use
`verify_choice` as its reward and `brute_force_choose` as ground truth at small sizes.

Two meanings of 'whole' pull opposite ways:  in positional base b a fraction terminates iff its
denominator's primes divide b; modulo m a fraction is a whole number iff its denominator is COPRIME to m.
"""
from itertools import product as iproduct
from math import gcd

from .crt import pairwise_coprime, prime_power_factors


def _primes_upto(n):
    sieve = [True] * (n + 1)
    out = []
    for i in range(2, n + 1):
        if sieve[i]:
            out.append(i)
            for j in range(i * i, n + 1, i):
                sieve[j] = False
    return out


def _smallest_prime(q):
    p = 2
    while q % p:
        p += 1
    return p


def candidate_options(bound, denominators=(), root_orders=(), must_cover=()):
    """({prime: [allowed channel values]}, forced_primes)."""
    need = {}
    for period in must_cover:
        for q in prime_power_factors(period):
            p = _smallest_prime(q)
            need[p] = max(need.get(p, 1), q)
    opts = {}
    for p in _primes_upto(bound):
        blocked = any(d % p == 0 for d in denominators)
        if blocked:
            if p in need:
                raise ValueError(f"must_cover needs prime {p} but a required denominator is divisible by {p}")
            continue
        powers, q = [], p
        while q <= bound:
            powers.append(q)
            q *= p
        if root_orders:
            powers = [p] if all((p - 1) % n == 0 for n in root_orders) else []
        powers = [q for q in powers if q >= need.get(p, 1)]
        if powers:
            opts[p] = powers
        elif p in need:
            raise ValueError(f"no allowed channel for required prime {p} within bound {bound}")
    forced = {p for p in need if p in opts}
    missing = set(need) - set(opts)
    if missing:
        raise ValueError(f"required primes {sorted(missing)} unavailable")
    return opts, forced


def _key(objective, chosen, prod):
    if objective == "channels":
        return (len(chosen), prod)
    if objective == "bits":
        return (prod, len(chosen))
    if objective == "width":
        return (max(chosen), len(chosen), prod)
    raise ValueError("objective must be 'channels', 'width' or 'bits'")


def _search(opts, forced, R, objective):
    primes = sorted(opts, key=lambda p: -max(opts[p]))
    best = {"key": None, "set": None}
    top = [max(opts[p]) for p in primes]
    suffix_ub = [1] * (len(primes) + 1)
    for i in range(len(primes) - 1, -1, -1):
        suffix_ub[i] = suffix_ub[i + 1] * top[i]
    forced_after = [0] * (len(primes) + 1)
    for i in range(len(primes) - 1, -1, -1):
        forced_after[i] = forced_after[i + 1] + (1 if primes[i] in forced else 0)

    def dfs(i, chosen, prod, n_forced_taken):
        if n_forced_taken == len(forced) and prod >= R:
            k = _key(objective, chosen, prod)
            if best["key"] is None or k < best["key"]:
                best["key"], best["set"] = k, tuple(sorted(chosen))
            return                                    # adding channels only worsens every objective
        if i == len(primes):
            return
        if prod * suffix_ub[i] < R:
            return
        if best["key"] is not None and objective == "channels" and len(chosen) + 1 > best["key"][0] and prod < R:
            return
        p = primes[i]
        for q in sorted(opts[p], reverse=True):
            dfs(i + 1, chosen + [q], prod * q, n_forced_taken + (1 if p in forced else 0))
        if p not in forced:
            dfs(i + 1, chosen, prod, n_forced_taken)

    dfs(0, [], 1, 0)
    return best["set"], best["key"]


def verify_choice(moduli, R, denominators=(), root_orders=(), must_cover=(), redundant=()):
    """Exact verifier: returns a list of violated requirements (empty list = valid)."""
    bad = []
    allm = list(moduli) + list(redundant)
    if not pairwise_coprime(allm):
        bad.append("moduli are not pairwise coprime")
    prod = 1
    for m in moduli:
        prod *= m
    if prod < R:
        bad.append(f"product {prod} < range {R}")
    for d in denominators:
        for m in allm:
            if gcd(d, m) != 1:
                bad.append(f"denominator {d} not invertible mod {m}")
    for n in root_orders:
        for m in allm:
            if (m - 1) % n or any(m % p == 0 for p in range(2, int(m ** 0.5) + 1)):
                bad.append(f"channel {m} lacks a primitive {n}-th root of unity (need prime m with n | m-1)")
    for period in must_cover:
        for q in prime_power_factors(period):
            if not any(m % q == 0 for m in moduli):
                bad.append(f"period {period}: prime-power part {q} not readable from any channel")
    if redundant:
        k = len(moduli)
        from itertools import combinations
        for combo in combinations(allm, k):
            p = 1
            for m in combo:
                p *= m
            if p < prod:
                bad.append("redundant set is not MDS (some %d channels cover less than the range)" % k)
                break
    return bad


def choose_moduli(range_needed, bound=64, denominators=(), root_orders=(), must_cover=(),
                  redundancy=0, objective="channels"):
    """Best channel set under the objective; returns a dict, or raises ValueError if infeasible."""
    R = range_needed
    if objective == "width":
        allopts, forced = candidate_options(bound, denominators, root_orders, must_cover)
        values = sorted({q for qs in allopts.values() for q in qs})
        result = None
        for W in values:
            opts = {p: [q for q in qs if q <= W] for p, qs in allopts.items()}
            opts = {p: qs for p, qs in opts.items() if qs}
            if any(p not in opts for p in forced):
                continue
            s, k = _search(opts, forced, R, "width")
            if s is not None:
                result = s
                break
    else:
        opts, forced = candidate_options(bound, denominators, root_orders, must_cover)
        result, _ = _search(opts, forced, R, objective)
    if result is None:
        raise ValueError("infeasible: no channel set within the bound covers the range")
    prod = 1
    for m in result:
        prod *= m
    red = ()
    if redundancy:
        cand, p = [], max(result) + 1
        while len(cand) < redundancy:
            if all(p % d for d in range(2, int(p ** 0.5) + 1)) and p > 1 \
                    and all(gcd(p, m) == 1 for m in result) and all(gcd(p, d) == 1 for d in denominators) \
                    and (not root_orders or all((p - 1) % n == 0 for n in root_orders)):
                cand.append(p)
            p += 1
        red = tuple(cand)
    return {"moduli": result, "redundant": red, "product": prod, "count": len(result), "max": max(result)}


def brute_force_choose(range_needed, bound, denominators=(), root_orders=(), must_cover=(), objective="channels"):
    """Reference optimum by exhaustive enumeration (use only for small bounds)."""
    opts, forced = candidate_options(bound, denominators, root_orders, must_cover)
    primes = sorted(opts)
    best, best_set = None, None
    for choice in iproduct(*[[None] + opts[p] for p in primes]):
        chosen = [q for q in choice if q is not None]
        if not chosen:
            continue
        if any(all(c is None or _smallest_prime(c) != p for c in choice) for p in forced):
            continue
        prod = 1
        for q in chosen:
            prod *= q
        if prod < range_needed:
            continue
        k = _key(objective, chosen, prod)
        if best is None or k < best:
            best, best_set = k, tuple(sorted(chosen))
    return best_set, best
