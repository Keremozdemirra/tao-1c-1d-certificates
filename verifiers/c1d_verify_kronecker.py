"""Exact score of a C_1d certificate (signed third autocorrelation constant).

f is a step function with n equal cells on [-1/4, 1/4] and integer cell values v_i.
f*f is piecewise linear with knots at multiples of h = 1/(2n), so max_t |f*f(t)| is a knot
value h * |c_m| with c = v*v (full discrete autoconvolution). With int f = h * sum(v):

    max_t |f*f| / (int f)^2 = 2n * max_m |c_m| / (sum v)^2      (exact rational)

c is computed exactly with Python integers via Kronecker substitution: v = p - q with p, q >= 0,
c = p*p - 2*p*q + q*q, each product one big-integer multiplication. Standard library only.

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


def unpack(x, nbytes, count):
    raw = x.to_bytes(count * nbytes, "little")
    return [int.from_bytes(raw[i * nbytes:(i + 1) * nbytes], "little") for i in range(count)]


def autoconvolution(v):
    n = len(v)
    p = [x if x > 0 else 0 for x in v]
    q = [-x if x < 0 else 0 for x in v]
    w = max(max(abs(x) for x in v).bit_length(), 1)
    # every coefficient of p*p, p*q, q*q is < n * 2^(2w), so this digit width never carries
    nbytes = (2 * w + n.bit_length() + 8) // 8 + 1
    P, Q = pack(p, nbytes), pack(q, nbytes)
    m = 2 * n - 1
    pp = unpack(P * P, nbytes, m)
    pq = unpack(P * Q, nbytes, m)
    qq = unpack(Q * Q, nbytes, m)
    return [a - 2 * b + c for a, b, c in zip(pp, pq, qq)]


def main():
    cert = json.load(open(sys.argv[1]))
    v = cert["v_int"]
    assert all(isinstance(x, int) for x in v), "cell values must be integers"
    n = len(v)
    assert cert.get("n", n) == n
    c = autoconvolution(v)
    s = sum(v)
    assert s != 0
    # consistency checks: sum of all c_m equals s^2, and c_0 = v_0^2
    assert sum(c) == s * s and c[0] == v[0] * v[0] and c[-1] == v[-1] * v[-1]
    mx = max(abs(x) for x in c)
    arg = [i for i, x in enumerate(c) if abs(x) == mx]
    ratio = Fraction(2 * n * mx, s * s)
    print("n =", n)
    print("sum v =", s)
    print("max |c_m| =", mx, "at m =", arg[:5], "(%d knots attain it)" % len(arg))
    print("ratio = %d / %d" % (ratio.numerator, ratio.denominator))
    digits = (ratio.numerator * 10 ** 30) // ratio.denominator
    print("ratio = %d.%s (truncated)" % (digits // 10 ** 30, str(digits % 10 ** 30).zfill(30)))
    if len(sys.argv) > 2:
        num, den = sys.argv[2].split("/")
        rec = Fraction(int(num), int(den))
        print("record =", float(rec))
        print("beats record:", ratio < rec, " margin =", float(rec - ratio))


if __name__ == "__main__":
    main()
