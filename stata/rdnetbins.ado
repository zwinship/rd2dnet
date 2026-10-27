*! version 0.2.0  rd2dnet
*! rdnetbins: turn the network coordinates from rdnetset into nearest-crossing
*! bins (an additional design; rerun rdnetset for the default design).
program define rdnetbins, rclass
    version 16
    syntax [, NOTOPology NOSIMilarity NOBRanches HS(string) KERNELS(string) ]

    capture confirm variable _rdn_id
    if _rc {
        di as err "network coordinates not found; run {cmd:rdnetset} first"
        exit 301
    }
    if "`kernels'" == "" local kernels triangular

    foreach v in _rdn_bin _rdn_branch _rdn_S {
        capture drop `v'
    }
    quietly gen long _rdn_bin = .
    quietly gen double _rdn_branch = .
    quietly gen double _rdn_S = .
    label variable _rdn_bin "rdnetbins: nearest boundary crossing (bin)"
    label variable _rdn_branch "rdnetbins: branch (edge where the path from the crossing leaves a junction; 0 trunk)"
    label variable _rdn_S "rdnetbins: edge dissimilarity to the crossing's segment"

    * options for Mata
    local notopology = ("`notopology'" != "")
    local nosimilarity = ("`nosimilarity'" != "")
    local nobranches = ("`nobranches'" != "")
    local hs `hs'
    local kernels `kernels'

    rd2dnet_mata
    mata: rdn_bins_entry()

    tempname B
    matrix `B' = r(bins)
    local M = r(N_bins)
    local hS = r(h_S)
    local rank = r(rank)
    char _dta[rdnet_design] bins

    di as txt "rd2dnet nearest-crossing bins"
    di as txt "  bins (crossings reached first by some observation): " as res `M'
    if `hS' < . di as txt "  similarity bandwidth h_S = " as res %8.4f `hS' as txt "  (features of rank " as res `rank' as txt ")"
    else di as txt "  no similarity kernel (h_S = inf)"
    di as txt "  variables created: _rdn_bin _rdn_branch _rdn_S"
    di as txt "  rdnet now estimates one effect per bin; estat wbate, weights(density|gls) aggregates;"
    di as txt "  rdnetfe pools the bins. Rerun rdnetset for the default design."

    return scalar N_bins = `M'
    return scalar h_S = `hS'
    return scalar rank = `rank'
    return matrix bins = `B'
end
