#!/usr/bin/env python3
"""Generates src/mcfarland_normal.h: the tables for McFarland's modified
ziggurat for the standard normal (rng_normal(method = "mcfarland")).

C. D. McFarland (2016), "A modified ziggurat algorithm for generating
exponentially and normally distributed pseudorandom numbers", Journal of
Statistical Computation and Simulation 86(7), 1281-1294,
doi:10.1080/00949655.2015.1060234. Reference implementation (MIT):
https://github.com/cd-mcfarland/fast_prng -- this follows its
create_layers.py construction, computed with mpmath at 60 digits instead of
numpy longdouble, so the result does not depend on the platform's long
double, and written as hex floats and exact integers.

Construction, with f(x) = exp(-x^2/2) on x >= 0 (area sqrt(pi/2)):

  * 256 slots of equal volume v = sqrt(pi/2) / 256.
  * Layers i = 0, 1, ...: the rectangle [0, X_i] x [Y_{i-1}, Y_i] with
    Y_i = f(X_i), Y_{-1} = 0 and X_i (f(X_i) - Y_{i-1}) = v, taking the long
    flat solution; stop when no further layer fits. That gives N = 253
    layers; X_N = 0 and Y_N = 1 close the tables. Every point of a layer is
    under f, so a draw in layer i is X_i * u, accepted outright.
  * The 3v left over is the tail (x > X_0, region j = 0) and the overhang
    above each layer (region j = 1..N: X_j <= x <= X_{j-1}, Y_{j-1} <= y
    <= f(x)), chosen by Walker/Vose alias sampling over 256 cells.
  * f is concave below x = 1 and convex above it. J_INFLECTION is the
    overhang straddling x = 1. Above the inflection index (x < 1) the curve
    lies above the chord and a point below the chord is accepted at once;
    below it (x > 1) the curve lies under the chord and a point above the
    chord is reflected. E_CONVEX / E_CONCAVE bound the largest distance
    between curve and chord in any overhang of each kind, in units of that
    overhang's height scaled by 2^63, so most overhang draws skip exp().

Output scaling: X and Y are divided by 2^63 (exact), so X[i] * (double)w
for a signed 64-bit word w is a uniform point of layer i with a sign, and
an overhang coordinate is v[j] * 2^63 + u * (v[j-1] - v[j]) for u in
[0, 2^63). Alias thresholds are 56-bit: region j is kept when
(word >> 8) < T[j], otherwise ALIAS[j] is used; the low 8 bits pick j.

Run from the package root (needs mpmath):
    python3 tools/generate-mcfarland-tables.py
"""
import sys
from mpmath import mp, mpf, exp, sqrt, pi, erf, erfc, findroot, floor, ceil

mp.dps = 60
SLOTS = 256
TWO63 = mpf(2) ** 63
TWO56 = 2 ** 56
# Headroom on the early-exit bounds, in units of 2^-63 of an overhang's
# height: the bounds are exact for the real tables; this covers the double
# rounding of the stored tables and of the interpolation, so a shortcut
# never decides a point the exp() test would decide the other way by more
# than rounding. 2^20 units is 1.1e-13 of the height.
E_SLACK = 2 ** 20


def f(x):
    return exp(-x * x / 2)


def area(a, b):  # integral of f over [a, b]
    return sqrt(pi / 2) * (erf(b / sqrt(2)) - erf(a / sqrt(2)))


def bisect(g, lo, hi, iters=260):
    glo = g(lo)
    for _ in range(iters):
        mid = (lo + hi) / 2
        gm = g(mid)
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return (lo + hi) / 2


v = sqrt(pi / 2) / SLOTS

# --- layers -----------------------------------------------------------------
X = []
yprev = mpf(0)
xprev = mpf(10)
while True:
    h = lambda x: x * (f(x) - yprev) - v
    if yprev == 0:
        # x f(x) peaks at x = 1; the long flat layer is the root above it
        lo = mpf(1)
    else:
        # h peaks where f(x)(1 - x^2) = yprev, on (0, xprev)
        dh = lambda x: f(x) * (1 - x * x) - yprev
        lo = bisect(dh, mpf(0), xprev)
    if h(lo) <= 0:
        break
    xi = bisect(h, lo, xprev)
    X.append(xi)
    yprev = f(xi)
    xprev = xi
N = len(X)
X.append(mpf(0))
Y = [f(x) for x in X]
assert N == 253, N
for i in range(N):  # every layer has volume v
    assert abs(X[i] * (Y[i] - (Y[i - 1] if i else 0)) - v) < mpf(10) ** -50

# --- regions and alias table -------------------------------------------------
V = [sqrt(pi / 2) * erfc(X[0] / sqrt(2)) / 2 * 2]  # tail: integral of f over [X_0, inf)
V[0] = sqrt(pi / 2) * erfc(X[0] / sqrt(2))
for j in range(1, N + 1):
    V.append(area(X[j], X[j - 1]) - Y[j - 1] * (X[j - 1] - X[j]))
assert all(x > 0 for x in V)
assert abs(sum(V) - (SLOTS - N) * v) < mpf(10) ** -45, sum(V) / v
V += [mpf(0)] * (SLOTS - len(V))
total = sum(V)
p = [x * SLOTS / total for x in V]  # mean 1

# Vose's alias method; keep probability prob[j], else alias[j]
prob = [mpf(0)] * SLOTS
alias = list(range(SLOTS))
small = [j for j in range(SLOTS) if p[j] < 1]
large = [j for j in range(SLOTS) if p[j] >= 1]
q = p[:]
while small and large:
    s = small.pop()
    l = large.pop()
    prob[s] = q[s]
    alias[s] = l
    q[l] = q[l] - (1 - q[s])
    (small if q[l] < 1 else large).append(l)
for j in large + small:
    prob[j] = mpf(1)
    alias[j] = j
T = [min(TWO56, int(floor(prob[j] * TWO56 + mpf(1) / 2))) for j in range(SLOTS)]
# the realised probabilities of each region
real = [mpf(0)] * SLOTS
for j in range(SLOTS):
    real[j] += mpf(T[j]) / TWO56 / SLOTS
    real[alias[j]] += (1 - mpf(T[j]) / TWO56) / SLOTS
for j in range(SLOTS):
    assert abs(real[j] - V[j] / total) < mpf(2) ** -54, (j, real[j], V[j] / total)

# --- inflection and early-exit bounds ----------------------------------------
J_INFL = next(j for j in range(1, N + 1) if X[j] < 1 < X[j - 1])


def max_gap(j, curve_above):
    """max over the overhang of |f - chord| / (Y_j - Y_{j-1})"""
    a, b = X[j], X[j - 1]
    slope = (Y[j - 1] - Y[j]) / (b - a)
    # extremum of f - chord where f'(x) = -x f(x) = slope
    g = lambda x: -x * f(x) - slope
    x = bisect(g, a, b)
    chord = Y[j] + (x - a) * slope
    gap = (f(x) - chord) if curve_above else (chord - f(x))
    assert gap > 0, (j, gap)
    return gap / (Y[j] - Y[j - 1])


e_convex = max(max_gap(j, True) for j in range(J_INFL + 1, N + 1))
e_concave = max(max_gap(j, False) for j in range(1, J_INFL))
E_CONVEX = int(ceil(e_convex * TWO63)) + E_SLACK
E_CONCAVE = int(ceil(e_concave * TWO63)) + E_SLACK

# --- output -------------------------------------------------------------------


def hexd(x):
    return float(x).hex()


def rows(vals, per):
    return ",\n".join("    " + ", ".join(vals[i:i + per]) for i in range(0, len(vals), per))


x0 = X[0]
out = []
out.append("""/* Generated by tools/generate-mcfarland-tables.py -- do not edit.
 *
 * Tables for McFarland's modified ziggurat for the standard normal,
 * rng_normal(method = "mcfarland"): C. D. McFarland (2016), J. Stat.
 * Comput. Simul. 86(7), 1281-1294, doi:10.1080/00949655.2015.1060234;
 * construction from https://github.com/cd-mcfarland/fast_prng (MIT),
 * computed here with mpmath at 60 digits. See the script for the layout.
 *
 * zurand_mcf_x[i], zurand_mcf_y[i]: layer edge X_i and f(X_i) = exp(-X_i^2/2),
 *   both divided by 2^63; index ZURAND_MCF_LAYERS closes the tables with
 *   X = 0, Y = 1 (scaled).
 * zurand_mcf_alias_t[j] / zurand_mcf_alias_j[j]: alias sampling over 256
 *   cells; keep region j when (word >> 8) < t[j], otherwise take alias_j[j].
 *   Region 0 is the tail beyond X_0, region j >= 1 the overhang above
 *   layer j - 1.
 */
#ifndef ZURAND_MCFARLAND_NORMAL_H
#define ZURAND_MCFARLAND_NORMAL_H
""")
out.append("#define ZURAND_MCF_LAYERS %d" % N)
out.append("#define ZURAND_MCF_J_INFLECTION %d" % J_INFL)
out.append("/* max curve-chord distance in convex (x < 1) / concave (x > 1) overhangs,")
out.append(" * in units of 2^-63 of the overhang height, plus %d units of headroom */" % E_SLACK)
out.append("#define ZURAND_MCF_E_CONVEX UINT64_C(%d)   /* %.6f */" % (E_CONVEX, float(e_convex)))
out.append("#define ZURAND_MCF_E_CONCAVE UINT64_C(%d)   /* %.6f */" % (E_CONCAVE, float(e_concave)))
out.append("static const double zurand_mcf_x0 = %s;       /* X_0 = %s */" % (hexd(x0), mp.nstr(x0, 20)))
out.append("static const double zurand_mcf_inv_x0 = %s;   /* 1 / X_0 */" % hexd(1 / x0))
out.append("")
out.append("static const double zurand_mcf_x[%d] = {\n%s\n};" % (N + 1, rows([hexd(x / TWO63) for x in X], 4)))
out.append("static const double zurand_mcf_y[%d] = {\n%s\n};" % (N + 1, rows([hexd(y / TWO63) for y in Y], 4)))
out.append("static const uint64_t zurand_mcf_alias_t[%d] = {\n%s\n};" % (
    SLOTS, rows(["UINT64_C(0x%014x)" % t for t in T], 4)))
out.append("static const unsigned char zurand_mcf_alias_j[%d] = {\n%s\n};" % (
    SLOTS, rows(["%3d" % a for a in alias], 16)))
out.append("\n#endif")

path = "src/mcfarland_normal.h"
with open(path, "w") as fh:
    fh.write("\n".join(out) + "\n")

sys.stderr.write("layers %d, inflection overhang %d, X_0 %s\n" % (N, J_INFL, mp.nstr(x0, 17)))
sys.stderr.write("regions: tail %.4f%%, convex %.4f%%, inflection %.4f%%, concave %.4f%% of the 3/256 edge\n" % (
    100 * float(V[0] / total), 100 * float(sum(V[J_INFL + 1:N + 1]) / total),
    100 * float(V[J_INFL] / total), 100 * float(sum(V[1:J_INFL]) / total)))
sys.stderr.write("E_convex %.6f  E_concave %.6f\n" % (float(e_convex), float(e_concave)))
sys.stderr.write("wrote %s\n" % path)
