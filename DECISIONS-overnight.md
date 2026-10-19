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

4. **R tooling.** Rtools is not installed on this machine, so the R package (which has
   C++ code) cannot be built locally at first. I installed Rtools 4.5 from CRAN into
   the user folder (`C:\Users\zwins\rtools45`) to run R tests locally; CI (R CMD check
   on five platforms) remains the authority.
