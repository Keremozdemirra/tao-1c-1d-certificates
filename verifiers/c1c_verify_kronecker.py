"""Exact score of a C_1c certificate (Barnard-Steinerberger autocorrelation constant, lower bound).

f = sum_i w_i 1_[i/K, (i+1)/K) with non-negative integer weights w_i. The autocorrelation
g(t) = int f(x) f(x+t) dx is linear between consecutive multiples of 1/K, with g(k/K) = A_k / K,
A_k = sum_i w_i w_{i+k}, and int f = sum(w) / K. So

    min_{0<=t<=1} g(t) / (int f)^2 = K * min_{0<=k<=K} A_k / (sum w)^2      (exact rational)

All A_k come from one big-integer product (Kronecker substitution of w and reversed w).
Standard library only.

Usage: python3 -I verify_kronecker.py certificate.json [record_numerator/record_denominator]
"""
import json
import sys
from fractions import Fraction


def pack(a, nbytes):
    buf = bytearray(len(a) * nbytes)
    for i, x in enumerate(a):
        if x:
            buf[i * nbytes:(i + 1) * nbytes] = x.to_bytes(nbytes, "little")
    return int.from_bytes(buf, "little")


def main():
    cert = json.load(open(sys.argv[1]))
    K = cert["K"]
    w = cert["weights_int"]
    assert all(isinstance(x, int) and x >= 0 for x in w), "weights must be non-negative integers"
    n = len(w)
    nbytes = (2 * max(w).bit_length() + n.bit_length() + 8) // 8 + 1
    prod = pack(w, nbytes) * pack(w[::-1], nbytes)
    raw = prod.to_bytes((2 * n - 1) * nbytes, "little")
    # coefficient of index n-1+k in w * reverse(w) is A_k
    A = [int.from_bytes(raw[(n - 1 + k) * nbytes:(n + k) * nbytes], "little") for k in range(min(K, n - 1) + 1)]
    A += [0] * (K + 1 - len(A))
    s = sum(w)
    assert A[0] == sum(x * x for x in w)
    mn = min(A)
    arg = [k for k, x in enumerate(A) if x == mn]
    score = Fraction(K * mn, s * s)
    print("K =", K, " cells =", n, " support = [0, %d/%d)" % (n, K))
    print("sum w =", s)
    print("min A_k =", mn, "at k =", arg[:5], "(%d lags attain it)" % len(arg))
    print("score = %d / %d" % (score.numerator, score.denominator))
    digits = (score.numerator * 10 ** 30) // score.denominator
    print("score = %d.%s (truncated)" % (digits // 10 ** 30, str(digits % 10 ** 30).zfill(30)))
    if len(sys.argv) > 2:
        num, den = sys.argv[2].split("/")
        rec = Fraction(int(num), int(den))
        print("record =", float(rec))
        print("beats record:", score > rec, " margin =", float(score - rec))


if __name__ == "__main__":
    main()
