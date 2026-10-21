*! version 0.2.0  rd2dnet
*! rdnet: boundary average treatment effect curve on a street network
program define rdnet, eclass
    version 16
    if replay() {
        if "`e(cmd)'" != "rdnet" error 301
        syntax [, Level(cilevel)]
        Display
        exit
    }
    syntax varname(numeric) [if] [in] [aweight] [,                                      ///
        h(numlist >0) HMATrix(name) p(integer 1) q(integer -1)                            ///
        DERiv(numlist min=2 max=2 integer >=0) TANGvec(numlist min=2 max=2)              ///
        COVs(varlist numeric) KERnel(string) BWSelect(string) vce(string)                 ///
        FUZzy(varname numeric) POOLed Level(cilevel) NOCBands REPP(integer 2000)          ///
        MINobs(integer -1) SCALEREGul(real 3) SCALEBIAScrct(real 1)                       ///
        CONLEYRESIDScale(real 3) CONLEYRESIDOrder(integer -1) SEED(integer -1) ]

    capture confirm variable _rdn_id
    if _rc {
        di as err "network coordinates not found; run {cmd:rdnetset} first"
        exit 301
    }
    if "`h'" != "" & "`hmatrix'" != "" {
        di as err "give h() or hmatrix(), not both"
        exit 198
    }
    if "`deriv'" != "" & "`tangvec'" != "" {
        di as err "give either deriv() or tangvec(), not both"
        exit 198
    }
    if "`kernel'" == "" local kernel triangular
    if "`bwselect'" == "" local bwselect mserd
    local bwselect = lower("`bwselect'")
    if !inlist("`bwselect'", "mserd", "msetwo", "imserd", "imsetwo", "cerrd", "certwo", "icerrd", "icertwo") {
        di as err "bwselect() must be mserd, msetwo, imserd, imsetwo, cerrd, certwo, icerrd or icertwo"
        exit 198
    }

    * vce(hc0|hc1|hc2|hc3 | cluster clustvar | conley # [network|euclidean])
    local vcetype hc3
    local cutoff
    local euclid 0
    local clorig
    if `"`vce'"' != "" {
        gettoken vt rest : vce
        local vt = lower("`vt'")
        if inlist("`vt'", "hc0", "hc1", "hc2", "hc3") {
            local vcetype `vt'
        }
        else if "`vt'" == "cluster" {
            local vcetype cluster
            gettoken clorig rest : rest
            confirm variable `clorig'
        }
        else if "`vt'" == "conley" {
            local vcetype conley
            gettoken cutoff rest : rest
            capture confirm number `cutoff'
            if _rc | "`cutoff'" == "" {
                di as err "vce(conley #): give a positive distance cutoff"
                exit 198
            }
            if `cutoff' <= 0 {
                di as err "vce(conley #): the cutoff must be positive"
                exit 198
            }
            local dt = lower(strtrim("`rest'"))
            if "`dt'" == "euclidean" local euclid 1
            else if !inlist("`dt'", "", "network") {
                di as err "vce(conley # [network|euclidean])"
                exit 198
            }
        }
        else {
            di as err "vce() must be hc0, hc1, hc2, hc3, cluster clustvar, or conley # [network|euclidean]"
            exit 198
        }
    }

    marksample touse
    markout `touse' _rdn_id `covs' `fuzzy'
    if "`clorig'" != "" markout `touse' `clorig', strok
    local wvar
    if "`weight'" != "" {
        tempvar w
        quietly gen double `w' `exp' if `touse'
        markout `touse' `w'
        quietly count if `touse' & `w' < 0
        if r(N) > 0 {
            di as err "weights must be non-negative"
            exit 402
        }
        local wvar `w'
    }
    local clvar
    if "`clorig'" != "" {
        tempvar clcode
        quietly egen long `clcode' = group(`clorig') if `touse'
        local clvar `clcode'
    }
    quietly count if `touse'
    if r(N) == 0 error 2000

    * options for Mata
    local yvar `varlist'
    local hlist `h'
    local hmat `hmatrix'
    if `q' < 0 local q
    if `minobs' < 0 local minobs
    if `conleyresidorder' < 0 local rorder
    else local rorder `conleyresidorder'
    if `seed' < 0 local seed
    local rscale `conleyresidscale'
    local pooled = ("`pooled'" != "")
    local cbands = ("`nocbands'" == "")
    local fuzzy `fuzzy'
    tempname resmat covmat bwmat
    local resmat `resmat'
    local covmat `covmat'
    local bwmat `bwmat'

    rd2dnet_mata
    mata: rdn_rdnet_entry()

    * e(b), e(V): bias-corrected estimates at the estimated evaluation points
    local J = rowsof(`resmat')
    local names
    local rows
    forvalues j = 1/`J' {
        if `resmat'[`j', 15] == 0 {
            local names `names' eval`j'
            local rows `rows' `j'
        }
    }
    local K : word count `rows'
    tempname b V
    if `K' > 0 {
        matrix `b' = J(1, `K', .)
        matrix `V' = J(`K', `K', .)
        forvalues k = 1/`K' {
            local jk : word `k' of `rows'
            matrix `b'[1, `k'] = `resmat'[`jk', 10]
            forvalues l = 1/`K' {
                local jl : word `l' of `rows'
                matrix `V'[`k', `l'] = `covmat'[`jk', `jl']
            }
        }
        matrix colnames `b' = `names'
        matrix colnames `V' = `names'
        matrix rownames `V' = `names'
        ereturn post `b' `V', esample(`touse') obs(`nobs') depname(`varlist')
    }
    else {
        ereturn post, esample(`touse') obs(`nobs') depname(`varlist')
    }
    ereturn matrix results = `resmat'
    capture confirm matrix `covmat'
    if !_rc ereturn matrix cov = `covmat'
    capture confirm matrix `bwmat'
    if !_rc ereturn matrix bandwidths = `bwmat'
    ereturn scalar N_eval = `J'
    ereturn scalar N_ok = `K'
    ereturn scalar level = `level'
    ereturn scalar p = `p'
    ereturn scalar q = `q'
    ereturn scalar min_obs = `minobs'
    if "`crit'" != "" & "`crit'" != "." ereturn scalar crit = `crit'
    if "`cutoff'" != "" ereturn scalar cutoff = `cutoff'
    ereturn local kernel `kname'
    ereturn local vce `vcetype'
    if "`clorig'" != "" ereturn local clustvar `clorig'
    if "`vcetype'" == "conley" ereturn local conley_distance = cond(`euclid', "euclidean", "network")
    ereturn local bwselect = cond("`h'`hmatrix'" == "", "`bwselect'", "manual")
    ereturn local fuzzy `fuzzy'
    ereturn local covs `covs'
    ereturn local deriv `deriv'
    ereturn local tangvec `tangvec'
    ereturn local pooled = cond(`pooled', "pooled", "")
    local costname : char _dta[rdnet_cost]
    ereturn local cost `costname'
    ereturn local depvar `varlist'
    ereturn local title "Network boundary discontinuity estimates"
    ereturn local estat_cmd rdnet_estat
    ereturn local cmd rdnet

    Display
end

program define Display
    tempname R
    matrix `R' = e(results)
    local J = rowsof(`R')
    local lev = e(level)
    local design = cond("`e(fuzzy)'" != "", "fuzzy", "sharp")
    di
    di as txt "rd2dnet: network boundary discontinuity estimates (`design')"
    di as txt "  distance = " as res "`e(cost)'" as txt "   kernel = " as res "`e(kernel)'" ///
        as txt "   p = " as res e(p) as txt "   q = " as res e(q) as txt "   vce = " as res "`e(vce)'"
    di as txt "  observations = " as res e(N) as txt "   evaluation points = " as res `J' ///
        as txt "   estimated = " as res e(N_ok)
    if e(crit) < . {
        di as txt "  uniform band critical value = " as res %6.3f e(crit)
    }
    di as txt "{hline 100}"
    di as txt "  eval" _col(9) "      arc" _col(19) "     h0" _col(28) "     h1" _col(37) "   n0" ///
        _col(44) "   n1" _col(51) "     tau" _col(61) "  tau_bc" _col(71) "   se_rb" ///
        _col(81) "[`lev'% robust CI]" _col(100) ""
    di as txt "{hline 100}"
    forvalues j = 1/`J' {
        local st = `R'[`j', 15]
        local stl = cond(`st' == 0, "", cond(`st' == 1, "masked", cond(`st' == 2, "insufficient", ///
            cond(`st' == 3, "singular", cond(`st' == 4, "no bandwidth", ///
            cond(`st' == 5, "weak first stage", "few clusters"))))))
        if `st' == 0 {
            di as txt %6.0f `j' _col(9) as res %9.1f `R'[`j', 2] _col(19) %7.1f `R'[`j', 4] ///
                _col(28) %7.1f `R'[`j', 5] _col(37) %5.0f `R'[`j', 6] _col(44) %5.0f `R'[`j', 7] ///
                _col(51) %8.3f `R'[`j', 8] _col(61) %8.3f `R'[`j', 10] _col(71) %8.3f `R'[`j', 11] ///
                _col(81) %8.3f `R'[`j', 12] _col(90) %8.3f `R'[`j', 13]
        }
        else {
            di as txt %6.0f `j' _col(9) as res %9.1f `R'[`j', 2] _col(19) as txt "(`stl')"
        }
    }
    di as txt "{hline 100}"
    di as txt "  tau: order-p estimate; tau_bc, se_rb and the CI: robust bias-corrected (order q)."
    tempname BS
    matrix `BS' = `R'[1..., 22]
    local n_s = 0
    forvalues j = 1/`J' {
        if `BS'[`j', 1] == 2 local ++n_s
    }
    if `n_s' > 0 {
        di as txt "  note: " as res `n_s' as txt " points fitted in s alone (few boundary crossings nearby)."
    }
    di as txt "  e(b) and e(V) hold tau_bc and its joint covariance at the estimated points."
end
