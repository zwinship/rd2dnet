# Overnight decisions (2026-10-19)

Each entry: the decision, the alternatives, and why. All are reversible.

1. **Working branch `stata-testing`** for the Stata fixes, merged through a PR once CI is green.

2. **Test expectation for evaluation point 1 changed (Stata test suite).** The suite expected
   eval 1 of the mock city to have status 4 ("no bandwidth") with automatic bandwidths.
   Since 0.2.0 that point is fitted in `s` alone (basis 2) and estimated, in Python too
   (`basis = "s"`, status ok). The check now asserts status 0 and basis 2.
   Alternative: keep the old expectation and change the code; rejected because the
   0.2.0 behaviour is intended and matches Python and R.

3. **Help example for segment networks made runnable.** `help rdnetset` used
   `import delimited roads.csv`, a file nobody has. The example now builds the segment
   frame from the mock city's streets with `frlink`/`frget` and mentions the CSV route
   in the text. Alternative: ship a `roads.csv`; rejected to keep the package small.

4. **R tooling.** Rtools is not installed on this machine, so the R package (C++ code) cannot
   be built locally. I tried a silent per-user install of Rtools 4.5 from CRAN; the
   installer asked for administrator rights and is waiting behind a Windows elevation
   prompt (processes `rtools45.exe`/`rtools45.tmp`, which I cannot stop without
   elevation). Cancel that prompt in the morning, or install Rtools if you want local R
   checks. R is verified through CI (R CMD check --as-cran on five platforms) and,
   for cross-language numbers, through a GitHub Actions job that runs R on the same data.

5. **s-only fallback triggered by pilot fits (needs your decision).** On a jittered 100 m
   grid (parity case 05), points are fitted in `s` alone although the order-p (s, o) fit
   is identified: only 2-3 streets cross the boundary within the bandwidth, so the
   cubic bias-pilot fit in `o` is rank-deficient, `rdbwnet` marks (s, o) "singular" and
   falls back. Options: keep (documented behaviour), or reduce only the pilot's degree
   in `o`. Left unchanged.

6. **Parity harness on its own branch (`parity-harness`)**; generated cases are not
   committed (`parity/.gitignore`); CI regenerates them.
