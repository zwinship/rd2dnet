* Installation check for rd2dnet: install with -net install- at a given commit,
* then run every help-file example and the postestimation commands from an
* empty directory, and list what is left behind in memory.
*
* usage:  do https://raw.githubusercontent.com/zwinship/rd2dnet/<sha>/stata/test/install_check.do <sha>

args sha
if "`sha'" == "" local sha main
version 16
clear all
set more off
capture mata: mata drop rdn_*()
capture program drop _all
* drop development copies from the ado-path (rd2dnet_tests.do adds one)
tokenize `"`c(adopath)'"', parse(";")
local i 1
while `"``i''"' != "" {
    if strpos(`"``i''"', "rd2dnet") capture adopath - `"``i''"'
    local ++i
}
adopath
capture ado uninstall rd2dnet
net install rd2dnet, from("https://raw.githubusercontent.com/zwinship/rd2dnet/`sha'/stata") replace

local home : env USERPROFILE
if "`home'" == "" local home : env HOME
cd "`home'"
capture mkdir rd2dnet_empty
cd rd2dnet_empty
local files : dir . files "*"
if `"`files'"' != "" di as err "note: the working directory is not empty: `files'"

* everything must come from the installed copy
which rdnetset
which rdnet
findfile lrd2dnet.mata
findfile rdnet_mock_obs.dta

global rdn_ifails 0
global rdn_ifaillist
capture program drop rdn_try
program define rdn_try
    * runs one command and records failures
    capture noisily `0'
    if _rc {
        global rdn_ifails = $rdn_ifails + 1
        global rdn_ifaillist `"$rdn_ifaillist"' `"`0' (rc `=_rc')"'
    }
end

* help rd2dnet
rdn_try rdnetmock, clear
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
rdn_try rdnet outcome [aw=far], vce(cluster block)
rdn_try estat wbate, weights(equal)
rdn_try rdnetdensity

* help rdnetset
rdn_try rdnetmock, clear
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
rdn_try matrix B = (100, 1600 \ 1000, 1600 \ 2000, 1600 \ 2900, 1600)
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) evalpoints(B) treated(treated)
rdn_try frame copy rdn_edges roads, replace
rdn_try frame roads: frlink m:1 from, frame(rdn_nodes id)
rdn_try frame roads: frget x1 = x y1 = y, from(rdn_nodes)
rdn_try frame roads: drop rdn_nodes
rdn_try frame roads: frlink m:1 to, frame(rdn_nodes id)
rdn_try frame roads: frget x2 = x y2 = y, from(rdn_nodes)
rdn_try rdnetset x y, edges(roads) ecoords(x1 y1 x2 y2) jointol(0.5) boundary(rdn_line) treatedside(right)
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) cost(travel_time) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) classvar(highway) penalty(500) connectors(motorway_link) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)

* help rdnet
rdn_try rdnetmock, clear
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
rdn_try rdnet outcome, seed(1)
rdn_try estat wbate, weights(equal)
rdn_try estat lbate
rdn_try rdnet
rdn_try test eval3 = eval4
rdn_try lincom eval3 - eval4
rdn_try estimates store A
rdn_try rdnet outcome [aw=far], h(600) vce(conley 150)
rdn_try rdnet outcome, h(800) p(2) deriv(1 0)
rdn_try generate lnfar = ln(far)
rdn_try rdnet outcome, h(600) covs(lnfar) vce(cluster block)
rdn_try rdnet outcome if abs(_rdn_s) >= 40, h(700)
rdn_try estat summarize

* help rdbwnet
rdn_try rdnetmock, clear
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
rdn_try rdbwnet outcome, bwselect(msetwo)
rdn_try matrix list r(bws)

* help rdnetdensity
rdn_try rdnetmock, clear
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)
rdn_try rdnetdensity
rdn_try rdnetdensity, normalize(none)

* help rdnetmock
rdn_try rdnetmock, clear
rdn_try rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)

* the help files open
foreach h in rd2dnet rdnetset rdnet rdbwnet rdnetdensity rdnetmock {
    rdn_try findfile `h'.sthlp
}

* what is left in memory (documented: frames rdn_nodes rdn_edges roads,
* matrices rdn_line rdn_region rdn_park B, Mata externals rdn_C rdn_DIM)
di _n as txt "{hline 70}" _n "left behind" _n "{hline 70}"
frames dir
matrix dir
macro dir
mata: direxternal("*")
local nw : dir . files "*"
di as txt "files written to the working directory: `nw'"

di _n "{hline 70}"
if $rdn_ifails == 0 di as res "INSTALL CHECK PASSED"
else {
    di as err "install check: $rdn_ifails commands failed"
    foreach f of global rdn_ifaillist {
        di as err "  FAILED: `f'"
    }
}
