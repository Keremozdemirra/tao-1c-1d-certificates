# Step-function certificates for Tao's constants C_1c and C_1d

Two exact certificates for constants in
[teorth/optimizationproblems](https://github.com/teorth/optimizationproblems):

| Constant | Bound | File |
|---|---|---|
| C_1d (upper bound) | 1.4556234448 | `c1d_n40960_1.4556234448.json` |
| C_1c (lower bound) | 0.4103454241 | `c1c_K98304_0.4103454241.json` |

Exact values, seeds and verification commands are in `VERIFY.md`; file
hashes in `SHA256SUMS`. The C_1d certificate refines the 10240-cell
certificate of teorth/optimizationproblems PR #269 (echul-seal), and the C_1c
certificate refines the K = 98304 certificate of PR #270 (V. Patankar).

Found and prepared with Claude (Anthropic) under the direction of Kerem
Özdemir. No human mathematician has checked the argument independently.

Certificates: CC BY 4.0. Verifier code: MIT.
