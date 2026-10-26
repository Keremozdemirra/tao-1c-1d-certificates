# Step-function certificates for Tao's constants C_1c and C_1d

Two exact certificates for constants in
[teorth/optimizationproblems](https://github.com/teorth/optimizationproblems):

| Constant | Bound | File |
|---|---|---|
| C_1d (upper bound) | 1.4556234448 | `c1d_n40960_1.4556234448.json` |
| C_1c (lower bound) | 0.4103454241, below the bound of Chen and Gu (see the note) | `c1c_K98304_0.4103454241.json` |

Exact values, seeds and verification commands are in `VERIFY.md`; file
hashes in `SHA256SUMS`. The C_1d certificate refines the 10240-cell
certificate of teorth/optimizationproblems PR #269 (echul-seal), and the C_1c
certificate refines the K = 98304 certificate of PR #270 (V. Patankar).

## Note on C_1c (26 October 2026)

The C_1c certificate is not a record. S. Chen and E. Gu, *An improved lower bound for the
minimum-autocorrelation constant*, arXiv:2610.18888 (8 October 2026), prove C_1c >= 0.4103779774,
above the 0.4103454241 of this certificate. We found their paper after opening the pull request
and withdrew its C_1c part. The certificate stays here so that its value can still be checked.
The C_1d certificate is not affected.

Found and prepared with Claude (Anthropic) under the direction of Kerem
Özdemir. No human mathematician has checked the argument independently.

Certificates and documentation: CC BY 4.0 (`LICENSE-data.md`). Verifier code: MIT (`LICENSE`).
