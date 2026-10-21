*! version 0.2.0  rd2dnet
*! rdbwnet: data-driven bandwidths for rdnet (port of rd2d's rdbw2d plug-in)
program define rdbwnet, rclass
    version 16
    syntax varname(numeric) [if] [in] [aweight] [, p(integer 1)                         ///
        DERiv(numlist min=2 max=2 integer >=0) TANGvec(numlist min=2 max=2)              ///
        COVs(varlist numeric) KERnel(string) BWSelect(string) vce(string)                 ///
        CLuster(varname) FUZzy(varname numeric) BWCHeck(integer -1)                      ///
        SCALEREGul(real 1) SCALEBIAScrct(real 1) ]

    capture confirm variable _rdn_id
    if _rc {
        di as err "network coordinates not found; run {cmd:rdnetset} first"
        exit 301
    }
    if "`deriv'" != "" & "`tangvec'" != "" {
        di as err "give either deriv() or tangvec(), not both"
        exit 198
    }
    if "`kernel'" == "" local kernel triangular
    if "`bwselect'" == "" local bwselect mserd
    local vcetype = lower("`vce'")
    if "`vcetype'" == "" local vcetype hc3
    if !inlist("`vcetype'", "hc0", "hc1", "hc2", "hc3") {
        di as err "vce() must be hc0, hc1, hc2 or hc3 (use cluster() for clustering)"
        exit 198
    }

    marksample touse
    markout `touse' _rdn_id `covs' `fuzzy'
    if "`cluster'" != "" markout `touse' `cluster', strok
    local wvar
    if "`weight'" != "" {
        tempvar w
        quietly gen double `w' `exp' if `touse'
        markout `touse' `w'
        local wvar `w'
    }
    local clvar
    if "`cluster'" != "" {
        tempvar clcode
        quietly egen long `clcode' = group(`cluster') if `touse'
        local clvar `clcode'
    }
    quietly count if `touse'
    if r(N) == 0 error 2000

    local yvar `varlist'
    if `bwcheck' < 0 local bwcheck
    tempname bwmat
    local bwmat `bwmat'
    rd2dnet_mata
    mata: rdn_rdbwnet_entry()

    local J = rowsof(`bwmat')
    di
    di as txt "rd2dnet bandwidth selection (`bwselect', p = `p')"
    di as txt "{hline 72}"
    di as txt "  eval" _col(10) "      arc" _col(22) "       h0" _col(34) "       h1" ///
        _col(46) "    n0" _col(54) "    n1" _col(62) "status"
    di as txt "{hline 72}"
    forvalues j = 1/`J' {
        local st = `bwmat'[`j', 3]
        local stl = cond(`st' == 0, "ok", cond(`st' == 1, "masked", cond(`st' == 2, "insufficient", "singular")))
        di as txt %6.0f `j' _col(10) as res %9.1f `bwmat'[`j', 2] _col(22) %9.2f `bwmat'[`j', 21] ///
            _col(34) %9.2f `bwmat'[`j', 22] _col(46) %6.0f `bwmat'[`j', 7] _col(54) %6.0f `bwmat'[`j', 8] ///
            _col(62) as txt "`stl'"
    }
    di as txt "{hline 72}"
    return matrix bws = `bwmat'
    return local bwselect `bwselect'
    return local kernel `kernel'
end
