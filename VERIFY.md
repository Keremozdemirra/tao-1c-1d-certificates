# Certificates: C_1c and C_1d

Two step-function certificates for constants on Terence Tao's optimization-constants page
(https://github.com/teorth/optimizationproblems, pages `constants/1c.md` and `1d.md`), frozen 2026-10-26.
The Erdős minimum overlap (C_1b) certificate is kept apart in `../release-1b/`.
Every number below is an exact rational computed from the integer data in the JSON files.

| File | Constant, direction | Exact value | Newest earlier claim (seed) | Margin |
|---|---|---|---|---|
| `c1d_n40960_1.4556234448.json` | C_1d, upper bound | 279676466330630562685427712000 / 192135175720687825901095191121 = 1.45562344470022455281... | 76358168487433178963290071040 / 52455219864513498497744894329 = 1.45568293650581519488... (open PR #269 by echul-seal, 26 Oct 2026; this was the seed) | 5.95e-5 lower |
| `c1c_K98304_0.4103454241.json` | C_1c, lower bound | 2574880088097860246697246720 / 6274908738633606834153301849 = 0.41034542418836093382... | 2594500369901311586790309888 / 6322755819680014072464230569 = 0.41034328129923190079... (open PR #270 by V. Patankar, 26 Oct 2026; this was the seed) | 2.14e-6 higher |

## Formats and scores

- **C_1d** (`n`, `v_int`): f has n equal cells on [-1/4, 1/4] with values v_i (signed integers). f*f is piecewise
  linear with knots at multiples of h = 1/(2n), so max_t |f*f| / (int f)^2 = 2n max_m |c_m| / (sum v)^2 with
  c = v*v (full discrete autoconvolution).
- **C_1c** (`K`, `weights_int`): f = sum_i w_i 1_[i/K,(i+1)/K) with positive integer weights (98368 cells, so the
  support is [0, 98368/98304)). The autocorrelation is linear between multiples of 1/K, so
  min_{0<=t<=1} int f(x) f(x+t) dx / (int f)^2 = K min_{0<=k<=K} A_k / (sum w)^2, A_k = sum_i w_i w_{i+k}.

## Re-running the checks

Standard library Python (3.9+) and a C compiler; the cold verifiers also use numpy.

```
python3 -I verifiers/c1d_verify_kronecker.py c1d_n40960_1.4556234448.json 76358168487433178963290071040/52455219864513498497744894329
python3 -I verifiers/c1c_verify_kronecker.py c1c_K98304_0.4103454241.json 2594500369901311586790309888/6322755819680014072464230569

cc -O2 -o c1d_direct verifiers/c1d_verify_direct.c && python3 -I verifiers/to_txt.py c1d_n40960_1.4556234448.json c1d.txt && ./c1d_direct c1d.txt
cc -O2 -o c1c_direct verifiers/c1c_verify_direct.c && python3 -I verifiers/to_txt.py c1c_K98304_0.4103454241.json c1c.txt && ./c1c_direct c1c.txt

python3 -I verifiers/cold/cold_c1d.py c1d_n40960_1.4556234448.json 76358168487433178963290071040/52455219864513498497744894329
python3 -I verifiers/cold/cold_c1c.py c1c_K98304_0.4103454241.json 2594500369901311586790309888/6322755819680014072464230569
```

`*_verify_kronecker.py` gets every lag from one big-integer product; `*_verify_direct.c` sums every lag in
128-bit integers. The `cold/` verifiers were written by a separate agent from the problem statements alone,
without access to the search code or to the other verifiers (int64-limb numpy convolution, cross-checked
against their own `__int128` loop in `xcheck.c`). All three give the fractions in the table. Certified bounds, rounded outward: C_1d <= 1.4556234448, C_1c >= 0.4103454241.

`SHA256SUMS` covers every file in this directory.
