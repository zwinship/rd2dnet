*! version 0.2.0  rd2dnet
*! rdnetset: declare the street network and boundary, and compute network
*! coordinates for the observations in memory.
program define rdnetset, rclass
    version 16
    syntax varlist(min=2 max=2 numeric) [if] [in], edges(name)                            ///
        [ nodes(name) boundary(name) ECOORDs(namelist min=4 max=4) JOINtol(real 1e-6) nodeid(name) nodex(name) nodey(name) from(name) to(name) cost(name) length(name)   ///
          DIRected direction(string) treatedregion(name) treatedside(string)                 ///
          treated(varname numeric) maskarcs(name) maskpolygons(name) classvar(name)          ///
          penalty(real 0) penaltymatrix(name) connectors(string) snapscale(real -1)          ///
          tol(real 1e-6) neval(integer 25) evalpoints(name) NOSPLIT snapnode ]

    gettoken xvar yvar : varlist
    local yvar = strtrim("`yvar'")

    * defaults for variable names in the node and edge frames
    if "`nodeid'" == "" local nodeid id
    if "`nodex'" == ""  local nodex x
    if "`nodey'" == ""  local nodey y
    if "`from'" == ""   local from from
    if "`to'" == ""     local to to
    if "`cost'" == ""   local cost length
    if "`length'" == "" local length length
    if "`direction'" == "" local direction both
    local direction = lower("`direction'")
    if !inlist("`direction'", "both", "out", "in") {
        di as err "direction() must be both, out or in"
        exit 198
    }
    local treatedside = lower("`treatedside'")
    if "`treatedregion'" != "" & "`treatedside'" != "" {
        di as err "give only one of treatedregion() or treatedside(left|right)"
        exit 198
    }
    if "`treatedregion'" == "" & "`treatedside'" == "" & "`treated'" == "" {
        di as err "give treatedregion(), treatedside(left|right), or treated() to infer the treated side"
        exit 198
    }
    if "`boundary'" == "" & "`evalpoints'" == "" {
        di as err "give boundary(), or evalpoints() (a J x 2 matrix of points along the boundary, as in rd2d)"
        exit 198
    }
    if "`ecoords'" == "" local ecoords x1 y1 x2 y2
    if "`treatedside'" != "" & !inlist("`treatedside'", "left", "right") {
        di as err "treatedside() must be left or right"
        exit 198
    }
    if `neval' < 1 {
        di as err "neval() must be a positive integer"
        exit 198
    }
    foreach f in nodes edges {
        if "``f''" == "" continue
        capture frame ``f'': quietly count
        if _rc {
            di as err "frame ``f'' not found"
            exit 111
        }
    }
    foreach m in boundary treatedregion maskarcs maskpolygons penaltymatrix evalpoints {
        if "``m''" != "" confirm matrix ``m''
    }

    marksample touse
    if "`treated'" != "" markout `touse' `treated'
    quietly count if `touse'
    if r(N) == 0 error 2000

    * rdnetset owns these variables: they are replaced on every call
    foreach v in _rdn_id _rdn_s _rdn_treated _rdn_arc {
        capture drop `v'
    }
    quietly gen long _rdn_id = sum(`touse') if `touse'
    quietly gen double _rdn_s = .
    quietly gen byte _rdn_treated = .
    quietly gen double _rdn_arc = .
    label variable _rdn_id "rdnetset: observation id"
    label variable _rdn_s "rdnetset: signed network distance to the boundary"
    label variable _rdn_treated "rdnetset: treated side"
    label variable _rdn_arc "rdnetset: where the observation reaches the boundary (arc length)"

    * options for Mata
    local fnodes `nodes'
    local fedges `edges'
    local fromv `from'
    local tov `to'
    local costv `cost'
    local lenv `length'
    local classv `classvar'
    local directed = ("`directed'" != "")
    local bmat `boundary'
    local rmat `treatedregion'
    local side `treatedside'
    local tvar `treated'
    local marcs `maskarcs'
    local mpoly `maskpolygons'
    local pmat `penaltymatrix'
    local emat `evalpoints'
    local nosplit = ("`nosplit'" != "")
    local jointol `jointol' 
    local snapnode = ("`snapnode'" != "")

    rd2dnet_mata
    mata: rdn_setup()

    local N = r(N)
    local nt = r(n_treated)
    local nloc = r(n_locations)
    local nev = r(n_eval)
    local nmask = r(n_masked)
    local nun = r(n_unreachable)
    local ncr = r(n_crossings)
    local nn = r(n_nodes)
    local ne = r(n_edges)
    local blen = r(length)
    tempname E
    matrix `E' = r(eval)

    char _dta[rdnet_set] 1
    char _dta[rdnet_cost] `cost'
    char _dta[rdnet_neval] `nev'

    di as txt "rd2dnet network coordinates"
    di as txt "  network: " as res `nn' as txt " nodes, " as res `ne' as txt " edges" ///
        as txt " (" as res `ncr' as txt " split at the boundary); distance = " as res "`cost'"
    di as txt "  observations: " as res `N' as txt " (" as res `nt' as txt " treated), " ///
        as res `nloc' as txt " distinct locations"
    di as txt "  evaluation points: " as res `nev' as txt " (" as res `nmask' as txt " masked)"
    if `nun' > 0 {
        di as txt "  note: " as res `nun' as txt " observations have no network path to the boundary; _rdn_s is missing for them"
    }
    di as txt "  variables created: _rdn_id _rdn_s _rdn_treated _rdn_arc"

    return scalar N = `N'
    return scalar N_treated = `nt'
    return scalar N_locations = `nloc'
    return scalar N_eval = `nev'
    return scalar N_masked = `nmask'
    return scalar N_unreachable = `nun'
    return scalar N_crossings = `ncr'
    return scalar N_nodes = `nn'
    return scalar N_edges = `ne'
    return scalar length = `blen'
    return matrix eval = `E'
end
