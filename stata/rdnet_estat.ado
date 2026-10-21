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
        * equal weight per distinct place (node_arc) among the estimated points;
        * points that snapped to the same place share its weight
        rd2dnet_mata
        mata: st_matrix(st_local("wmat"), rdn_place_weights(st_matrix(st_local("R"))))
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
    syntax [, MIN ABSmax]
    if "`min'" != "" & "`absmax'" != "" {
        di as err "give min or absmax, not both"
        exit 198
    }
    local kind = cond("`min'" != "", "min", cond("`absmax'" != "", "absmax", "max"))
    rd2dnet_mata
    mata: rdn_lbate_entry()
    di
    local what = cond("`kind'" == "min", "Smallest", cond("`kind'" == "absmax", "Largest absolute", "Largest"))
    local f = cond("`kind'" == "min", "min", cond("`kind'" == "absmax", "max|.|", "max"))
    di as txt "`what' boundary average treatment effect (LBATE)"
    di as txt "  `f' tau (order p)" _col(22) "= " as res %9.4f r(tau)
    di as txt "  `f' tau_bc" _col(22) "= " as res %9.4f r(tau_bc)
    di as txt "  uniform CI" _col(22) "= [" as res %9.4f r(ci_lo) as txt ", " as res %9.4f r(ci_hi) as txt "]"
    return add
end
