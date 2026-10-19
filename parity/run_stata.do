* Cross-language parity: run every case in parity/cases with the Stata port and
* compare with the Python reference results written by make_cases.py.
*
* Run from the directory that contains parity/ (the repository root, or the
* folder made by parity/fetch_stata.do), with rd2dnet on the ado-path:
*     do parity/run_stata.do
* Prints one line per case and a summary; differences are relative to
* max(1, |Python value|), and status codes must agree exactly.

version 16
set more off
global par_tol 1e-8
global par_worst 0
global par_bad
global par_err

capture program drop rdn_par_note
program define rdn_par_note
    args case label diff
    local diff = `diff'
    if `diff' >= . local diff 1e300
    if `diff' > $par_worst global par_worst = `diff'
    if `diff' > $par_tol {
        di as err "  case `case' `label': rel diff " %9.2e `diff'
        global par_bad `"$par_bad"' `"`case' `label' `: display %9.2e `diff''"'
    }
end

capture program drop rdn_par_load
program define rdn_par_load
    args dir
    import delimited using "`dir'/obs.csv", clear asdouble case(preserve) varnames(1)
    foreach f in nodes edges bnd rgn pr {
        capture frame drop `f'
    }
    frame create nodes
    frame nodes: import delimited using "`dir'/nodes.csv", clear asdouble case(preserve) varnames(1)
    frame create edges
    frame edges: import delimited using "`dir'/edges.csv", clear asdouble case(preserve) varnames(1) stringcols(5)
    frame create bnd
    frame bnd {
        import delimited using "`dir'/boundary.csv", clear asdouble varnames(1)
        mkmat x y, matrix(BL)
    }
    capture confirm file "`dir'/region.csv"
    if !_rc {
        frame create rgn
        frame rgn {
            import delimited using "`dir'/region.csv", clear asdouble varnames(1)
            mkmat x y, matrix(RG)
        }
    }
end

capture program drop rdn_par_coords
program define rdn_par_coords
    args dir case
    tempvar d
    quietly gen double `d' = abs(_rdn_s - py_s) / max(1, abs(py_s))
    quietly replace `d' = 0 if _rdn_s >= . & py_s >= .
    quietly replace `d' = 1e300 if `d' >= .
    quietly summarize `d', meanonly
    rdn_par_note `case' "coords s" r(max)
    quietly count if _rdn_treated != py_t
    rdn_par_note `case' "coords treated (count)" r(N)
    quietly replace `d' = abs(_rdn_arc - py_arc) / max(1, abs(py_arc))
    quietly replace `d' = 0 if _rdn_arc >= . & py_arc >= .
    quietly replace `d' = 1e300 if `d' >= .
    quietly summarize `d', meanonly
    rdn_par_note `case' "coords arc" r(max)
end

capture program drop rdn_par_ref
program define rdn_par_ref
    * loads a Python reference file into matrix PYREF
    args file vars
    capture frame drop pr
    frame create pr
    frame pr {
        import delimited using "`file'", clear asdouble varnames(1)
        mkmat `vars', matrix(PYREF)
    }
end

capture program drop rdn_par_res
program define rdn_par_res
    args file case label
    rdn_par_ref "`file'" "h0 h1 n0 n1 tau se tau_bc se_rb status basis"
    tempname R P
    matrix `R' = e(results)
    matrix `P' = PYREF
    local cols 4 5 6 7 8 9 10 11
    local J = rowsof(`P')
    if rowsof(`R') != `J' {
        rdn_par_note `case' "`label' rows" 1e300
        exit
    }
    local nst 0
    local worst 0
    forvalues j = 1/`J' {
        if `R'[`j', 15] != `P'[`j', 9] local ++nst
        else if `P'[`j', 9] == 0 {
            forvalues k = 1/8 {
                local c : word `k' of `cols'
                local d = abs(`R'[`j', `c'] - `P'[`j', `k']) / max(1, abs(`P'[`j', `k']))
                if `d' >= . local d 1e300
                if `d' > `worst' local worst `d'
            }
            if `R'[`j', 22] != `P'[`j', 10] local ++nst
        }
    }
    rdn_par_note `case' "`label' status/basis mismatches" `nst'
    rdn_par_note `case' "`label' values" `worst'
end

capture program drop rdn_par_dens
program define rdn_par_dens
    args file case label
    local names diff diff_bc se h_left h_right
    tempname V
    matrix `V' = (r(diff), r(diff_bc), r(se), r(h_left), r(h_right))
    rdn_par_ref "`file'" "`names'"
    local worst 0
    forvalues k = 1/5 {
        local d = abs(`V'[1, `k'] - PYREF[1, `k']) / max(1, abs(PYREF[1, `k']))
        if `d' >= . local d 1e300
        if `d' > `worst' local worst `d'
    }
    rdn_par_note `case' "`label' values" `worst'
end

capture program drop rdn_par_bw
program define rdn_par_bw
    args file case label
    tempname B
    matrix `B' = r(bws)
    rdn_par_ref "`file'" "h0 h1 status"
    local J = rowsof(PYREF)
    local nst 0
    local worst 0
    forvalues j = 1/`J' {
        if `B'[`j', 3] != PYREF[`j', 3] local ++nst
        else if PYREF[`j', 3] == 0 {
            forvalues k = 1/2 {
                local c = 20 + `k'
                local d = abs(`B'[`j', `c'] - PYREF[`j', `k']) / max(1, abs(PYREF[`j', `k']))
                if `d' >= . local d 1e300
                if `d' > `worst' local worst `d'
            }
        }
    }
    rdn_par_note `case' "`label' status mismatches" `nst'
    rdn_par_note `case' "`label' values" `worst'
end

file open fh using "parity/cases/cases.txt", read text
file read fh line
local cases
while r(eof) == 0 {
    local cases `cases' `line'
    file read fh line
}
file close fh

foreach cs of local cases {
    di as txt _n "{hline 60}" _n "parity case `cs'"
    local before $par_worst
    global par_worst 0
    capture noisily do "parity/cases/`cs'/run.do"
    if _rc {
        di as err "  case `cs': error rc = " _rc
        global par_err `"$par_err `cs'(`=_rc')"'
    }
    di as txt "  case `cs': max rel diff " %9.2e $par_worst
    global par_worst = max(`before', $par_worst)
}

di _n "{hline 60}"
di as txt "worst relative difference: " %9.2e $par_worst
if `"$par_bad"' == "" & `"$par_err"' == "" di as res "PARITY PASSED"
else {
    foreach b of global par_bad {
        di as err "  FAIL `b'"
    }
    if `"$par_err"' != "" di as err "  errors in cases: $par_err"
}
