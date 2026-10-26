#!/usr/bin/env python3
"""Cold, from-scratch verifier for C_1d (signed third autocorrelation) upper-bound certificates.

Problem (Tao optimizationproblems 1d): for real f supported on [-1/4,1/4],
    C_1d <= max_t |(f*f)(t)| / (int f)^2 .

Reduction used here (derived independently, see --selftest which checks it by brute force):
  f = sum_i v_i 1_[x_i, x_i+h),  x_i = -1/4 + i h,  h = 1/(2n),  i = 0..n-1.
  1_[a,a+h) * 1_[b,b+h) is the tent with knots a+b, a+b+h, a+b+2h and peak value h.
  Hence f*f is continuous, piecewise linear, supported on [-1/2,1/2], with knots at -1/2 + m h
  (m = 0..2n), and value  (f*f)(-1/2 + m h) = h * c_{m-1},  c_k = sum_{i+j=k} v_i v_j
  (c_{-1} = c_{2n-1} = 0).  |piecewise linear| is maximal at a knot, so
      max_t |f*f| = h * max_{0<=k<=2n-2} |c_k|,   int f = h * S,  S = sum_i v_i
      ratio = h max|c| / (h S)^2 = 2n * max|c_k| / S^2.

Exact computation: c_k via 16-bit limbs and numpy.convolve on int64 arrays (bounds asserted),
recombined with Python ints; cross-checked lag by lag against an independent C __int128 loop
(xcheck.c, same directory) and by direct Python big-int sums at selected lags; float64 FFT sanity.

Usage:  python3 -I cold_c1d.py CERT.json RECORD [--claim CLAIM] [--no-c] [--key v_int]
        python3 -I cold_c1d.py --selftest
RECORD / CLAIM are fractions ("num/den") or decimal strings, parsed exactly by fractions.Fraction.
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile
import time
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
B = 16  # limb width in bits


def dec(fr, digits=30):
    getcontext().prec = digits + 10
    q = Decimal(fr.numerator) / Decimal(fr.denominator)
    return format(q, "." + str(digits) + "f")


def make_limbs(vals, nl):
    a = np.array(vals, dtype=np.int64)
    out, r = [], a.copy()
    mask = (1 << B) - 1
    for _ in range(nl - 1):
        out.append(r & mask)
        r = r >> B  # arithmetic shift = floor division, top limb stays signed
    out.append(r)
    # exact reconstruction check in Python ints
    rec = [0] * len(vals)
    for j, d in enumerate(out):
        for i, x in enumerate(d.tolist()):
            rec[i] += x << (B * j)
    assert rec == list(vals), "limb reconstruction failed"
    for d in out:
        assert int(np.abs(d).max()) <= (1 << B)
    return out


def exact_autoconv(vals):
    """c_k = sum_{i+j=k} v_i v_j for k=0..2n-2 as Python ints, by limb convolution."""
    n = len(vals)
    maxbits = max(abs(x) for x in vals).bit_length()
    nl = max(1, -(-maxbits // B))
    limbs = make_limbs(vals, nl)
    # each convolution entry is a sum of <= n terms with |limb|<=2^B; a combined slot adds <= nl pair terms (x2 for symmetry)
    bound = 2 * nl * n * (1 << (2 * B))
    assert bound < (1 << 62), "int64 overflow risk: %d" % bound
    slots = [None] * (2 * nl - 1)
    for a in range(nl):
        for b in range(a, nl):
            t = np.convolve(limbs[a], limbs[b])
            if a != b:
                t = 2 * t
            s = a + b
            slots[s] = t if slots[s] is None else slots[s] + t
    res = [0] * (2 * n - 1)
    for s, arr in enumerate(slots):
        for k, x in enumerate(arr.tolist()):
            res[k] += x << (B * s)
    return res


def run_c(mode, arrays, extra, workdir):
    """Compile xcheck.c (own code) and run it; returns the exact integers it prints."""
    src = os.path.join(HERE, "xcheck.c")
    exe = os.path.join(workdir, "xcheck")
    if not os.path.exists(exe):
        subprocess.run(["clang", "-O2", "-o", exe, src], check=True)
    paths = []
    for idx, arr in enumerate(arrays):
        p = os.path.join(workdir, "in%d.txt" % idx)
        with open(p, "w") as f:
            f.write("\n".join(str(x) for x in arr) + "\n")
        paths.append(p)
    out = os.path.join(workdir, "out.txt")
    subprocess.run([exe, mode] + paths + [out] + extra, check=True)
    with open(out) as f:
        return [int(line) for line in f.read().split()]


def check_hygiene(vals, n):
    assert isinstance(vals, list), "not a list"
    assert len(vals) == n, "declared n=%d but len=%d" % (n, len(vals))
    assert all(type(x) is int for x in vals), "non-int entry (float/bool/str?)"
    S = sum(vals)
    assert S != 0, "integral is zero"
    return S


def evaluate(vals, n, use_c=True, workdir=None):
    S = check_hygiene(vals, n)
    t0 = time.time()
    c = exact_autoconv(vals)
    print("numpy-limb exact autoconvolution done in %.1fs, len(c)=%d" % (time.time() - t0, len(c)))
    m = max(abs(x) for x in c)
    kstar = [k for k, x in enumerate(c) if abs(x) == m]
    print("max|c_k| = %d at lag(s) %s (sign %s)" % (m, kstar, [("+" if c[k] > 0 else "-") for k in kstar]))
    # direct big-int spot checks: argmax, top 10 |c|, 200 random lags
    order = sorted(range(len(c)), key=lambda k: -abs(c[k]))[:10]
    rnd = random.Random(12345)
    lags = set(order) | set(rnd.sample(range(len(c)), 200)) | {0, 1, len(c) - 1, n - 1}
    for k in sorted(lags):
        d = sum(vals[i] * vals[k - i] for i in range(max(0, k - n + 1), min(k, n - 1) + 1))
        assert d == c[k], "direct sum mismatch at lag %d" % k
    print("direct Python big-int sums agree at %d lags (argmax, top-10, 200 random, edges)" % len(lags))
    if use_c:
        t0 = time.time()
        maxabs = max(abs(x) for x in vals)
        assert n * maxabs * maxabs < (1 << 126) and maxabs < (1 << 62)
        cc = run_c("conv", [vals], [], workdir)
        assert len(cc) == len(c) and cc == c, "C __int128 loop disagrees with numpy-limb result"
        print("C __int128 direct loop agrees with numpy-limb result at ALL %d lags (%.1fs)" % (len(c), time.time() - t0))
    ratio = Fraction(2 * n * m, S * S)
    # float64 FFT sanity
    N = 1 << (2 * n - 1).bit_length()
    x = np.array(vals, dtype=np.float64)
    cf = np.fft.irfft(np.fft.rfft(x, N) ** 2, N)[: 2 * n - 1]
    rf = 2 * n * float(np.abs(cf).max()) / float(x.sum()) ** 2
    print("float64 FFT ratio = %.12f ; exact ratio = %.12f ; rel diff = %.2e" % (rf, float(ratio), abs(rf - float(ratio)) / float(ratio)))
    assert abs(rf - float(ratio)) / float(ratio) < 1e-9, "float sanity check failed"
    return ratio, S, m


def selftest():
    """Brute-force check of the reduction on small random instances, exact rationals, many non-knot t."""
    rnd = random.Random(7)
    for trial in range(40):
        n = rnd.randint(1, 7)
        v = [rnd.randint(-9, 9) for _ in range(n)]
        if sum(v) == 0:
            v[0] += 1
        h = Fraction(1, 2 * n)
        xs = [Fraction(-1, 4) + i * h for i in range(n)]

        def F(t):  # (f*f)(t) = sum_ij v_i v_j |cell_i  ∩  (t - cell_j)|
            tot = Fraction(0)
            for i in range(n):
                for j in range(n):
                    lo = max(xs[i], t - xs[j] - h)
                    hi = min(xs[i] + h, t - xs[j])
                    if hi > lo:
                        tot += v[i] * v[j] * (hi - lo)
            return tot

        c = [sum(v[i] * v[k - i] for i in range(max(0, k - n + 1), min(k, n - 1) + 1)) for k in range(2 * n - 1)]
        for m in range(0, 2 * n + 1):
            t = Fraction(-1, 2) + m * h
            expect = h * (c[m - 1] if 1 <= m <= 2 * n - 1 else 0)
            assert F(t) == expect, ("knot formula", n, v, m)
        sample_max = Fraction(0)
        for m in range(-2, 2 * n + 2):
            for s in (Fraction(0), Fraction(1, 7), Fraction(1, 3), Fraction(1, 2), Fraction(5, 7), Fraction(1)):
                t = Fraction(-1, 2) + (m + s) * h
                val = F(t)
                lin = (1 - s) * F(Fraction(-1, 2) + m * h) + s * F(Fraction(-1, 2) + (m + 1) * h)
                assert val == lin, "not linear between knots"
                sample_max = max(sample_max, abs(val))
        assert sample_max == h * max(abs(x) for x in c), "sup over samples != knot max"
        # outside support is zero
        assert F(Fraction(-1, 2) - h) == 0 and F(Fraction(1, 2) + h) == 0
        S = sum(v)
        sup_ratio = sample_max / (h * S) ** 2
        assert sup_ratio == Fraction(2 * n * max(abs(x) for x in c), S * S)
    print("selftest OK: knot formula, piecewise linearity, zero outside support, ratio = 2n max|c|/S^2 (40 random cases)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert", nargs="?")
    ap.add_argument("record", nargs="?", help="record to beat (upper bound; smaller is better)")
    ap.add_argument("--claim", default=None)
    ap.add_argument("--key", default="v_int")
    ap.add_argument("--no-c", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    cert = json.load(open(a.cert))
    n, vals = cert["n"], cert[a.key]
    print("certificate:", a.cert, "| n =", n, "| entries =", len(vals))
    with tempfile.TemporaryDirectory() as wd:
        ratio, S, m = evaluate(vals, n, use_c=not a.no_c, workdir=wd)
    print("S = sum v =", S)
    print("max|c_k| =", m)
    print("EXACT ratio = %d / %d  (reduced: %d / %d)" % (2 * n * m, S * S, ratio.numerator, ratio.denominator))
    print("ratio = %s" % dec(ratio))
    if a.claim:
        cl = Fraction(a.claim)
        print("claim   = %s ; exact equality with claim: %s" % (dec(cl), ratio == cl))
    if a.record:
        rec = Fraction(a.record)
        print("record  = %s" % dec(rec))
        print("ratio < record (strictly better upper bound): %s" % (ratio < rec))
        print("record - ratio = %s  (= %s as float)" % (dec(rec - ratio, 20), float(rec - ratio)))


if __name__ == "__main__":
    main()
