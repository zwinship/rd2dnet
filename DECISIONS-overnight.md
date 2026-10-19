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

7. **Exact distances on directed networks and with transition costs: split edges at the
   locations.** Alternatives: (a) state-level "arc entry" distances in the engine (no graph
   growth, but a rewrite of every distance routine in three languages, including the Conley
   kernel); (b) split only for these metrics (chosen). Undirected plain-cost networks, the
   common case, are untouched. Cost: the network grows by one node per distinct location for
   directed/transition metrics. `coords.network` therefore carries `rdnet_point_*` nodes.

8. **Covariate side weights 1/h_side² (rd2d's choice)** applied in all languages, also in the
   Conley auxiliary fit and for `pooled`. Changes estimates only when h0 != h1.

9. **Leverage clip 1e-8 (rd2d) instead of 1e-12.**

10. **Parity tolerance 1e-7** instead of 1e-8: automatic bandwidths differ by up to 3.5e-8
    between Python and R from rounding in near-ill-conditioned fits.

11. **AppStream sign-in.** The Stata session expired around 06:30 and the link now asks for an
    AWS sign-in; per the rules I did not sign in. Stata changes after commit 253b969 are
    untested in Stata. PR #5 is merged only if CI is green, on the basis that main is
    unpublished and the Stata port has never been released; please run the Stata suite
    before any release.

## Waiting for your decision (from the audits)

- **Near-singular q-fits:** (a) flag points whose order-q fit has leverage 1 or a condition
  number above a threshold (new status "unstable"); (b) when only the q-fit is singular in
  (s, o), keep the p-fit in (s, o) and compute the bias correction in s alone, instead of
  refitting everything in s alone. I recommend both; not done because it changes results.
- **Density bandwidth selector** is not rddensity's (normal-reference pilots, `comb` default).
  rd2dnet uses rd2d-style constants and `common = min(hl, hr)`. Together with the
  over-rejection on round zones, consider porting rddensity's selector.
- **min_obs/bwcheck default:** rd2d's `rd2d()` uses 50 + p + 1; rd2dnet uses 20 + #terms(q).
- **Cluster SEs:** rdnet uses rd2d's `vce="hc0"` + CR1; rd2d's default (`hc1` + cluster) is
  about 1.7% larger. Also clustering under-covers in several Monte Carlo runs; Conley does well.
- **Uniform band:** rd2d uses repp=1000, a fixed seed and type-2 quantiles; rd2dnet uses
  repp=2000, `seed=None`. Consider a fixed default seed for reproducibility.
- **Off-network snapping cost for travel-time metrics** (currently free).
