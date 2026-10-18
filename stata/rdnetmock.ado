*! version 0.2.0  rd2dnet
*! rdnetmock: load the bundled mock city (works offline)
program define rdnetmock
    version 16
    syntax [, clear]
    if "`clear'" == "" & c(changed) {
        di as err "no; data in memory would be lost (use the clear option)"
        exit 4
    }
    quietly findfile rdnet_mock_obs.dta
    local fobs `"`r(fn)'"'
    quietly findfile rdnet_mock_nodes.dta
    local fnodes `"`r(fn)'"'
    quietly findfile rdnet_mock_edges.dta
    local fedges `"`r(fn)'"'
    use `"`fobs'"', clear
    foreach f in rdn_nodes rdn_edges {
        capture frame drop `f'
        frame create `f'
    }
    frame rdn_nodes: use `"`fnodes'"', clear
    frame rdn_edges: use `"`fedges'"', clear
    * boundary along 60th Street (y = 1600); treated zone to the south; park masked
    matrix rdn_line = (-400, 1600 \ 2600, 1600)
    matrix rdn_region = (-100000, -100000 \ 100000, -100000 \ 100000, 1600 \ -100000, 1600)
    matrix rdn_park = (750, 1600 \ 1750, 1600 \ 1750, 2800 \ 750, 2800)
    di as txt "mock city loaded: " as res _N as txt " parcels in memory; frames " ///
        as res "rdn_nodes" as txt " (" as res "450" as txt " intersections) and " ///
        as res "rdn_edges" as txt " (streets);"
    di as txt "matrices " as res "rdn_line rdn_region rdn_park" as txt " hold the boundary, treated zone and park. Next:"
    di as txt `"  {cmd:rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)}"'
end
