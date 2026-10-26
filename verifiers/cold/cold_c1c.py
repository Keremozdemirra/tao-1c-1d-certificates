#!/usr/bin/env python3
"""Cold, from-scratch verifier for C_1c (Barnard-Steinerberger) lower-bound certificates.

Problem (Tao optimizationproblems 1c): for non-negative f in L^1(R),
    min_{0<=t<=1} int f(x) f(x+t) dx <= C_1c (int f)^2,
so any admissible f gives the lower bound  C_1c >= min_{0<=t<=1} A(t) / (int f)^2,  A(t) = int f(x)f(x+t)dx.

Reduction (derived independently; --selftest checks it by brute force with exact rationals):
  f = sum_{i=0}^{N-1} w_i 1_[i/K,(i+1)/K), w_i >= 0.  The cross-correlation of two cells i and j is the
  tent  (1/K) max(0, 1 - |j - i - K t|).  So, with R_d = sum_i w_i w_{i+d} (d >= 0),
      A(t) = (1/K) sum_d R_|d| max(0, 1 - |d - K t|)      (t >= 0 uses only d >= 0 terms),
  and for K t in [k, k+1]:  A(t) = (1/K)[ (1-s) R_k + s R_{k+1} ],  s = K t - k.
  A is piecewise linear with knots exactly at k/K; t = 1 is the knot k = K, so the min over [0,1]
  is attained at a knot:  min_{0<=t<=1} A = (1/K) min_{0<=k<=K} R_k.   int f = S/K, S = sum w_i.
      score = K * min_{0<=k<=K} R_k / S^2 .
  (The support [0, N/K) may be longer than 1; lags k <= K only are needed because t <= 1.)

Exact computation: R_k via 21-bit limbs and numpy int64 dot products (bounds asserted), recombined with
Python ints; cross-checked at every lag 0..K against an independent C __int128 loop (xcheck.c) and by
direct Python big-int sums at selected lags; float64 FFT sanity check.

Usage:  python3 -I cold_c1c.py CERT.json RECORD [--claim CLAIM] [--no-c] [--no-numpy]
        python3 -I cold_c1c.py --selftest
RECORD / CLAIM are fractions ("num/den") or decimal strings, parsed exactly by fractions.Fraction.
RECORD is a lower bound: a certificate beats it when score > RECORD.
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
B = 21  # limb width in bits


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
        r = r >> B
    out.append(r)
    rec = [0] * len(vals)
    for j, d in enumerate(out):
        for i, x in enumerate(d.tolist()):
            rec[i] += x << (B * j)
    assert rec == list(vals), "limb reconstruction failed"
    for d in out:
        assert int(np.abs(d).max()) <= (1 << B)
    return out


def exact_autocorr_numpy(vals, maxlag):
    """R_k = sum_i w_i w_{i+k}, k=0..maxlag, via limbs and int64 dot products."""
    N = len(vals)
    maxbits = max(abs(x) for x in vals).bit_length()
    nl = max(1, -(-maxbits // B))
    limbs = make_limbs(vals, nl)
    assert N * (1 << (2 * B)) < (1 << 62), "int64 overflow risk in a single dot product"
    res = []
    t0 = time.time()
    for k in range(maxlag + 1):
        tot = 0
        for a in range(nl):
            la = limbs[a][: N - k]
            for b in range(nl):
                tot += int(np.dot(la, limbs[b][k:])) << (B * (a + b))
        res.append(tot)
        if k and k % 20000 == 0:
            print("   numpy lag %d / %d  (%.0fs)" % (k, maxlag, time.time() - t0), flush=True)
    return res


def run_c(mode, arrays, extra, workdir):
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


def evaluate(vals, K, use_c=True, use_np=True, workdir=None):
    N = len(vals)
    assert type(K) is int and K > 0
    assert all(type(x) is int for x in vals), "non-int entry (float/bool/str?)"
    assert min(vals) >= 0, "negative weight"
    S = sum(vals)
    assert S > 0
    assert N >= K + 1, "support must exceed K cells or lag K has zero correlation"
    print("N = %d cells, K = %d (support [0, %d/%d) = [0, %.6f)), min w = %d, max w = %d, zeros = %d"
          % (N, K, N, K, N / K, min(vals), max(vals), sum(1 for x in vals if x == 0)))
    maxabs = max(vals)
    assert N * maxabs * maxabs < (1 << 126) and maxabs < (1 << 62)
    R = None
    if use_np:
        t0 = time.time()
        R = exact_autocorr_numpy(vals, K)
        print("numpy-limb exact autocorrelation lags 0..K done in %.1fs" % (time.time() - t0))
    if use_c:
        t0 = time.time()
        Rc = run_c("acorr", [vals], [str(K)], workdir)
        print("C __int128 autocorrelation lags 0..K done in %.1fs" % (time.time() - t0))
        if R is not None:
            assert Rc == R, "C and numpy-limb results disagree"
            print("C __int128 loop agrees with numpy-limb result at ALL %d lags" % (K + 1))
        else:
            R = Rc
    assert R is not None
    mn = min(R)
    kst = [k for k, x in enumerate(R) if x == mn]
    print("min_k R_k = %d at lag(s) %s" % (mn, kst))
    # direct Python big-int spot checks
    rnd = random.Random(2468)
    order = sorted(range(K + 1), key=lambda k: R[k])[:12]
    lags = set(order) | set(rnd.sample(range(K + 1), 60)) | {0, 1, K - 1, K}
    for k in sorted(lags):
        d = sum(vals[i] * vals[i + k] for i in range(N - k))
        assert d == R[k], "direct sum mismatch at lag %d" % k
    print("direct Python big-int sums agree at %d lags (argmin, 12 smallest, 60 random, edges)" % len(lags))
    print("R_0 = %d ; R_K = %d ; smallest 5 (lag:R/S^2*K) = %s" % (
        R[0], R[K], [(k, float(Fraction(K * R[k], S * S))) for k in order[:5]]))
    score = Fraction(K * mn, S * S)
    # float64 FFT sanity
    N2 = 1 << (2 * N - 1).bit_length()
    x = np.array(vals, dtype=np.float64)
    ac = np.fft.irfft(np.abs(np.fft.rfft(x, N2)) ** 2, N2)[: K + 1]
    sf = K * float(ac.min()) / float(x.sum()) ** 2
    print("float64 FFT score = %.12f ; exact score = %.12f ; rel diff = %.2e" % (sf, float(score), abs(sf - float(score)) / float(score)))
    assert abs(sf - float(score)) / float(score) < 1e-9, "float sanity check failed"
    return score, S, mn


def selftest():
    """Brute-force check of the reduction on small random instances with exact rationals."""
    rnd = random.Random(11)
    for trial in range(40):
        K = rnd.randint(1, 6)
        N = K + rnd.randint(1, 4)
        w = [rnd.randint(0, 9) for _ in range(N)]
        if sum(w) == 0:
            w[0] = 1
        h = Fraction(1, K)
        xs = [i * h for i in range(N)]

        def A(t):  # int f(x) f(x+t) dx = sum_ij w_i w_j |cell_i ∩ (cell_j - t)|
            tot = Fraction(0)
            for i in range(N):
                for j in range(N):
                    lo = max(xs[i], xs[j] - t)
                    hi = min(xs[i] + h, xs[j] + h - t)
                    if hi > lo:
                        tot += w[i] * w[j] * (hi - lo)
            return tot

        R = [sum(w[i] * w[i + k] for i in range(N - k)) for k in range(N)]
        for k in range(0, K + 1):
            assert A(Fraction(k, K)) == Fraction(R[k] if k < N else 0, K), ("knot", w, K, k)
        sample_min = None
        for k in range(0, K):
            for s in (Fraction(0), Fraction(1, 7), Fraction(1, 3), Fraction(1, 2), Fraction(5, 7), Fraction(1)):
                t = (k + s) * h
                val = A(t)
                lin = (1 - s) * A(Fraction(k, K)) + s * A(Fraction(k + 1, K))
                assert val == lin, "not linear between knots"
                sample_min = val if sample_min is None else min(sample_min, val)
        sample_min = min(sample_min, A(Fraction(1)))
        assert sample_min == Fraction(min(R[: K + 1]), K), "min over samples != knot min"
        S = sum(w)
        assert sample_min / (Fraction(S, K)) ** 2 == Fraction(K * min(R[: K + 1]), S * S)
    print("selftest OK: knot formula A(k/K)=R_k/K, piecewise linearity, min over t in [0,1] = knot min, score = K min R / S^2 (40 random cases)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert", nargs="?")
    ap.add_argument("record", nargs="?", help="record to beat (lower bound; larger is better)")
    ap.add_argument("--claim", default=None)
    ap.add_argument("--key", default="weights_int")
    ap.add_argument("--no-c", action="store_true")
    ap.add_argument("--no-numpy", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    cert = json.load(open(a.cert))
    K, vals = cert["K"], cert[a.key]
    print("certificate:", a.cert, "| K =", K, "| entries =", len(vals))
    with tempfile.TemporaryDirectory() as wd:
        score, S, mn = evaluate(vals, K, use_c=not a.no_c, use_np=not a.no_numpy, workdir=wd)
    print("S = sum w =", S)
    print("min R_k  =", mn)
    print("EXACT score = %d / %d  (reduced: %d / %d)" % (K * mn, S * S, score.numerator, score.denominator))
    print("score = %s" % dec(score))
    if a.claim:
        cl = Fraction(a.claim)
        print("claim = %s ; exact equality with claim: %s" % (dec(cl), score == cl))
    if a.record:
        rec = Fraction(a.record)
        print("record = %s" % dec(rec))
        print("score > record (strictly better lower bound): %s" % (score > rec))
        print("score - record = %s  (= %s as float)" % (dec(score - rec, 20), float(score - rec)))


if __name__ == "__main__":
    main()
