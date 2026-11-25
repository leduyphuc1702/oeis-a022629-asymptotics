"""Verification for "Asymptotics of the coefficients of prod (1 + k x^k)" (OEIS A022629).

1. Exact coefficients by dynamic programming up to NEXACT, compared with the OEIS data
   and with the b-file (downloaded once to b022629.txt).
2. The bounds log K! <= log a(n) <= log K! + O(sqrt(n)/log n) of Theorem 1, and the
   conjectured main term sqrt(n/2)(log(2n) - 2).
3. The saddle-point formula of Theorem 2: a(n) ~ exp(Phi_n(s_n)) / sqrt(2 pi V(s_n)),
   evaluated at s_n with m(s_n) = n, up to n = 10000 using the b-file.

Usage: python3 verify_a022629.py [NEXACT] [--no-download]
"""
import os
import sys
import urllib.request
from math import exp, lgamma, log, log1p, pi, sqrt, isqrt

DATA = [1, 1, 2, 5, 7, 15, 25, 43, 64, 120, 186, 288, 463, 695, 1105, 1728, 2525, 3741, 5775,
        8244, 12447, 18302, 26424, 37827, 54729, 78330, 111184, 159538, 225624, 315415,
        444708, 618666, 858165, 1199701, 1646076, 2288961, 3150951, 4303995, 5870539,
        8032571, 10881794, 14749051, 19992626]

args = [a for a in sys.argv[1:] if not a.startswith('--')]
NEXACT = int(args[0]) if args else 3000

# ------------------------------------------------------------ 1. coefficients
a = [1] + [0] * NEXACT
for k in range(1, NEXACT + 1):
    for i in range(NEXACT, k - 1, -1):
        a[i] += k * a[i - k]
assert a[:len(DATA)] == DATA
print("1. DP coefficients a(0..%d); first %d agree with the OEIS data" % (NEXACT, len(DATA)))

bfile = 'b022629.txt'
if not os.path.exists(bfile) and '--no-download' not in sys.argv:
    req = urllib.request.Request('https://oeis.org/A022629/b022629.txt', headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=300) as f:
        data = f.read()
    with open(bfile + '.part', 'wb') as g:
        g.write(data)
    os.replace(bfile + '.part', bfile)
B = {}
if os.path.exists(bfile):
    for line in open(bfile):
        p = line.split()
        if len(p) == 2 and not line.startswith('#'):
            B[int(p[0])] = int(p[1])
    common = [n for n in range(NEXACT + 1) if n in B]
    assert all(a[n] == B[n] for n in common)
    print("   b-file: %d terms (n <= %d); all %d terms with n <= %d agree with the DP"
          % (len(B), max(B), len(common), NEXACT))


def coef(n):
    return a[n] if n <= NEXACT else B[n]


def logint(x):
    sh = max(x.bit_length() - 900, 0)
    return log(x >> sh) + sh * log(2)


# ------------------------------------------------------------ 2./3. asymptotics
def kmax(n):
    K = (isqrt(8 * n + 1) - 1) // 2
    assert K * (K + 1) // 2 <= n < (K + 1) * (K + 2) // 2
    return K


def sums(s):
    """log Z(s), m(s), V(s) for Z(s) = prod (1 + k e^{-sk})."""
    lz = m = v = 0.0
    k = 1
    while True:
        u = log(k) - s * k
        if u < -60 and k > 1 / s:
            break
        lz += max(u, 0.0) + log1p(exp(-abs(u)))
        p = 1 / (1 + exp(-u))
        m += k * p
        v += k * k * p * (1 - p)
        k += 1
    return lz, m, v


def saddle(n):
    lo, hi = 1e-6, 5.0
    for _ in range(200):
        mid = sqrt(lo * hi)
        if sums(mid)[1] > n:
            lo = mid
        else:
            hi = mid
    s = sqrt(lo * hi)
    lz, m, v = sums(s)
    return s, s * n + lz, v


print("2./3. n, K, log a(n), log K!, main term, gap/(sqrt(n)/log n), saddle error")
print("      saddle error = log a(n) - (Phi_n(s_n) - log(2 pi V)/2)")
rows = []
for n in [100, 200, 500, 1000, 2000, 3000, 5000, 7000, 10000]:
    if n > NEXACT and n not in B:
        continue
    K = kmax(n)
    la = logint(coef(n))
    lk = lgamma(K + 1)
    main = sqrt(n / 2) * (log(2 * n) - 2)
    s, phi, v = saddle(n)
    err = la - (phi - 0.5 * log(2 * pi * v))
    gap = (la - lk) / (sqrt(n) / log(n))
    assert la >= lk
    rows.append((n, K, la, lk, main, gap, phi, err))
    print("  %6d %4d %11.4f %11.4f %11.4f %8.4f %11.4f %10.6f" % (n, K, la, lk, main, gap, phi, err))

errs = [r[-1] for r in rows]
assert all(abs(e) < 0.02 for e in errs)
assert all(abs(x) > abs(y) for x, y in zip(errs, errs[1:]))  # saddle error decreases to 0
assert all(r[5] < 6 for r in rows)  # gap = O(sqrt(n)/log n)
print("ALL CHECKS PASSED")
