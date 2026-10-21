*! version 0.2.0  rd2dnet
*! rdnetmanip: manipulation (sorting) test at the boundary (Bugni, Crippa and Restrepo 2026)
program define rdnetmanip, rclass
    version 16
    syntax [if] [in] [, q(numlist integer >1) KMax(integer -1) REPs(integer 10000) Level(cilevel) ///
        SEED(integer 0) ]

    capture confirm variable _rdn_id
    if _rc {
        di as err "network coordinates not found; run {cmd:rdnetset} first"
        exit 301
    }
    marksample touse, novarlist
    markout `touse' _rdn_id
    local qlist `q'
    if `kmax' < 1 & `kmax' != -1 {
        di as err "kmax() must be at least 1"
        exit 198
    }
    if `kmax' < 0 local kmax
    rd2dnet_mata
    mata: rdn_manip_entry()
    tempname P
    matrix `P' = r(path)
    di
    di as txt "rd2dnet manipulation test (Bugni, Crippa and Restrepo)"
    di as txt "  locations nearest the boundary: q = " as res r(q) as txt ", within " ///
        as res %9.4g r(distance) as txt " of it; " as res r(n_treated) as txt " treated"
    di as txt "  S = " as res %7.3f r(statistic) as txt "   p-value = " as res %6.4f r(pvalue) ///
        as txt "   critical value (`level'%) = " as res %6.3f r(crit) as txt "   draws = " as res r(reps)
    if rowsof(`P') > 1 {
        di as txt "  {hline 44}"
        di as txt "       q" _col(14) "statistic" _col(26) "p-value" _col(36) "distance"
        forvalues k = 1/`=rowsof(`P')' {
            di as res %8.0f `P'[`k', 1] _col(14) %9.3f `P'[`k', 2] _col(26) %7.4f `P'[`k', 3] ///
                _col(36) %8.4g `P'[`k', 6]
        }
    }
    return add
end
