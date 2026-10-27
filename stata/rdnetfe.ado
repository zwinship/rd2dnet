*! version 0.2.0  rd2dnet
*! rdnetfe: pooled nearest-crossing BATE with crossing fixed effects and
*! propensity weights (run rdnetbins first).
program define rdnetfe, eclass
    version 16
    if replay() {
        if "`e(cmd)'" != "rdnetfe" error 301
        Display
        exit
    }
    syntax varname(numeric) [if] [in] [aweight] [,                                      ///
        h(numlist >0) HMATrix(name) p(integer 1) q(integer -1) KERnel(string)            ///
        BWSelect(string) vce(string) COVs(varlist numeric) NOIPW NOBINSlopes             ///
        Level(cilevel) MINobs(integer 1) SCALEREGul(real 3) CR1 ]

    capture confirm variable _rdn_bin
    if _rc {
        di as err "nearest-crossing bins not found; run {cmd:rdnetset} and {cmd:rdnetbins} first"
        exit 301
    }
    if "`h'" != "" & "`hmatrix'" != "" {
        di as err "give h() or hmatrix(), not both"
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
    markout `touse' _rdn_id `covs'
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
    if "`q'" == "" local qv = `p' + 1
    else local qv `q'
    local noipw = ("`noipw'" != "")
    local nobinslopes = ("`nobinslopes'" != "")
    local fuzzy

    rd2dnet_mata
    mata: rdn_fe_entry()

    tempname b V B
    matrix `B' = r(bins)
    matrix `b' = (r(tau_bc))
    matrix `V' = (r(se_rb)^2)
    matrix colnames `b' = BATE
    matrix colnames `V' = BATE
    matrix rownames `V' = BATE
    local tau = r(tau)
    local se = r(se)
    local tau_bc = r(tau_bc)
    local se_rb = r(se_rb)
    local ci_lo = r(ci_lo)
    local ci_hi = r(ci_hi)
    local pv = r(pvalue)
    local df = r(df)
    local nused = r(N_bins_used)
    local nfit = r(N)
    ereturn post `b' `V', esample(`touse') obs(`nobs') depname(`varlist')
    ereturn scalar tau = `tau'
    ereturn scalar se = `se'
    ereturn scalar tau_bc = `tau_bc'
    ereturn scalar se_rb = `se_rb'
    ereturn scalar ci_lo = `ci_lo'
    ereturn scalar ci_hi = `ci_hi'
    ereturn scalar pvalue = `pv'
    if `df' < . ereturn scalar df = `df'
    ereturn scalar N_bins = rowsof(`B')
    ereturn scalar N_bins_used = `nused'
    ereturn scalar N_fit = `nfit'
    ereturn scalar level = `level'
    ereturn scalar p = `p'
    ereturn scalar q = `qv'
    ereturn matrix bins = `B'
    ereturn local terms `terms'
    ereturn local ipw = cond(`noipw', "", "ipw")
    ereturn local vce `vcetype'
    if "`vcetype'" == "cluster" ereturn local clustvce = cond("`cr1'" != "", "cr1", "cr2")
    if "`clorig'" != "" ereturn local clustvar `clorig'
    if "`cutoff'" != "" ereturn scalar cutoff = `cutoff'
    ereturn local kernel `kernel'
    ereturn local depvar `varlist'
    ereturn local title "Pooled nearest-crossing BATE (crossing fixed effects)"
    ereturn local cmd rdnetfe
    Display
end

program define Display
    local lev = e(level)
    di
    di as txt "rd2dnet: pooled nearest-crossing BATE with crossing fixed effects"
    di as txt "  p = " as res e(p) as txt "   q = " as res e(q) as txt "   vce = " as res "`e(vce)'" ///
        as txt "   bins used = " as res e(N_bins_used) as txt " of " as res e(N_bins) ///
        as txt "   observations in the fit = " as res e(N_fit)
    di as txt "  propensity weights: " as res cond("`e(ipw)'" != "", "yes (density-weighted BATE)", "no")
    di as txt "  slope terms: " as res cond("`e(terms)'" != "", "`e(terms)'", "none")
    di as txt "{hline 78}"
    di as txt "       tau        se    tau_bc     se_rb    [`lev'% robust CI]      p-value"
    di as res %10.4f e(tau) %10.4f e(se) %10.4f e(tau_bc) %10.4f e(se_rb) "  " ///
        %9.4f e(ci_lo) %9.4f e(ci_hi) "   " %9.4f e(pvalue)
    di as txt "{hline 78}"
    di as txt "  tau: order-p estimate; tau_bc, se_rb, the CI and p-value: bias-corrected (order q)."
end
