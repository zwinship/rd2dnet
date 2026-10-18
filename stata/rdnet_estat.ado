*! version 0.2.0  rd2dnet
*! estat after rdnet: aggregate effects (WBATE, LBATE)
program define rdnet_estat, rclass
    version 16
    if "`e(cmd)'" != "rdnet" error 301
    gettoken sub rest : 0, parse(" ,")
    local sub = lower("`sub'")
    if "`sub'" == "wbate" {
        Wbate `rest'
    }
    else if "`sub'" == "lbate" {
        Lbate `rest'
    }
    else {
        estat_default `0'
    }
    return add
end

program define Wbate, rclass
    syntax , Weights(string)
    tempname R wmat
    matrix `R' = e(results)
    local J = rowsof(`R')
    capture confirm matrix `weights'
    if !_rc {
        matrix `wmat' = `weights'
    }
    else if lower("`weights'") == "equal" {
        matrix `wmat' = J(`J', 1, 0)
        forvalues j = 1/`J' {
            if `R'[`j', 15] == 0 matrix `wmat'[`j', 1] = 1
        }
    }
    else {
        capture numlist "`weights'"
        if _rc {
            di as err "weights() must be a matrix name, a list of numbers, or equal"
            exit 198
        }
        local ws `r(numlist)'
        matrix `wmat' = J(`: word count `ws'', 1, 0)
        local k 0
        foreach x of local ws {
            local ++k
            matrix `wmat'[`k', 1] = `x'
        }
    }
    local wmat `wmat'
    rd2dnet_mata
    mata: rdn_wbate_entry()
    local lev = e(level)
    di
    di as txt "Weighted boundary average treatment effect (WBATE)"
    di as txt "  tau (order p)    = " as res %9.4f r(tau)
    di as txt "  tau_bc           = " as res %9.4f r(tau_bc)
    di as txt "  robust se        = " as res %9.4f r(se_rb)
    di as txt "  `lev'% robust CI   = [" as res %9.4f r(ci_lo) as txt ", " as res %9.4f r(ci_hi) as txt "]"
    di as txt "  p-value          = " as res %9.4f r(pvalue)
    return add
end

program define Lbate, rclass
    syntax
    rd2dnet_mata
    mata: rdn_lbate_entry()
    di
    di as txt "Largest boundary average treatment effect (LBATE)"
    di as txt "  max tau (order p) = " as res %9.4f r(tau)
    di as txt "  max tau_bc        = " as res %9.4f r(tau_bc)
    di as txt "  uniform CI        = [" as res %9.4f r(ci_lo) as txt ", " as res %9.4f r(ci_hi) as txt "]"
    return add
end
