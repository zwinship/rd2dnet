# Cross-language parity

Random designs are run in Python (the reference), R and Stata, and every coordinate,
estimate, standard error, bandwidth and status must agree to 1e-7 (relative to
max(1, |Python value|)).

The designs cover:

- networks: jittered grids, rivers with tributaries, ring boundaries, directed
  networks with `out` and `in`;
- costs: transition costs and an edge-cost column;
- estimation: fuzzy designs, covariates, weights, clusters, Conley, derivatives, `p = 2`,
  every kernel and most `bwselect` options;
- the density test and `rdbwnet`.

Node positions are jittered and lengths are not multiples of each other. This avoids
exact ties between equidistant boundary nodes, which Python breaks differently from R
and Stata.

```bash
python parity/make_cases.py --n 40 --seed 1      # cases + Python results in parity/cases/
Rscript parity/run_r.R parity/cases                # R results (needs the R package installed)
python parity/compare.py parity/cases --prefix r   # compare R with Python
```

The `parity` GitHub Actions job runs these three steps on every pull request that
touches `src/`, `R-package/` or `parity/`. When it fails, the cases are uploaded as an
artifact.

For Stata, put the cases where Stata can read them (a clone, or a copy of
`parity/cases` next to `parity/run_stata.do`). Then, from the folder that contains
`parity/`, with rd2dnet installed or on the ado-path:

```stata
do parity/run_stata.do
```

It prints one line per case and ends with `PARITY PASSED` or a list of the differences.
`parity/fetch_stata.do <sha>` downloads `run_stata.do` and the cases from GitHub at a
commit. The cases are not committed on `main`, so push them to a scratch branch first
with `git add -f parity/cases`.
