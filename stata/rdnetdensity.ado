*! version 0.1.0  rd2dnet
*! rdnetdensity: density (sorting) test at a network boundary (port of rddensity)
program define rdnetdensity, rclass
    version 16
    syntax [if] [in] [aweight] [, EVALindex(integer -1) WINDow(real -1) h(numlist max=2 >0) ///
        p(integer 2) q(integer -1) NORMalize(string) KERnel(string) BWSelect(string)          ///
        SCALEREGul(real 1) ]

    capture confirm variable _rdn_id
    if _rc {
        di as err "network coordinates not found; run {cmd:rdnetset} first"
        exit 301
    }
    if `evalindex' >= 0 & `window' < 0 {
        di as err "give window() together with evalindex()"
        exit 198
    }
    if "`kernel'" == "" local kernel triangular
    local normalize = lower("`normalize'")
    local bwselect = lower("`bwselect'")
    marksample touse, novarlist
    markout `touse' _rdn_id
    local wvar
    if "`weight'" != "" {
        tempvar w
        quietly gen double `w' `exp' if `touse'
        markout `touse' `w'
        local wvar `w'
    }
    if `evalindex' < 0 local evalindex
    if `window' < 0 local window
    if `q' < 0 local q
    local hlist `h'
    rd2dnet_mata
    mata: rdn_density_entry()
    di
    di as txt "rd2dnet density test (normalize = " as res "`r(normalize)'" as txt ")"
    di as txt "  jump f+ - f- (order p)  = " as res %10.4g r(diff)
    di as txt "  bias-corrected jump     = " as res %10.4g r(diff_bc)
    di as txt "  jackknife se            = " as res %10.4g r(se)
    di as txt "  z = " as res %7.3f r(z) as txt "   p-value = " as res %6.4f r(pvalue) ///
        as txt "   log(f+/f-) = " as res %7.3f r(theta)
    di as txt "  observations (control/treated) = " as res r(N_left) as txt " / " as res r(N_right) ///
        as txt ", within bandwidth " as res r(N_h_left) as txt " / " as res r(N_h_right)
    di as txt "  bandwidths (control/treated)   = " as res %9.1f r(h_left) as txt " / " as res %9.1f r(h_right)
end
