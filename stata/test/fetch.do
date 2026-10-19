* Download the Stata files of rd2dnet at a given commit into ./rd2dnet_test
* usage:  do https://raw.githubusercontent.com/zwinship/rd2dnet/<sha>/stata/test/fetch.do <sha>
args sha
if "`sha'" == "" local sha main
capture mkdir rd2dnet_test
cd rd2dnet_test
local base "https://raw.githubusercontent.com/zwinship/rd2dnet/`sha'/stata"
foreach f in lrd2dnet.mata rd2dnet_mata.ado rdnetset.ado rdnet.ado rdnet_estat.ado rdbwnet.ado rdnetdensity.ado rdnetmock.ado rd2dnet.sthlp rdnetset.sthlp rdnet.sthlp rdbwnet.sthlp rdnetdensity.sthlp rdnetmock.sthlp rdnet_mock_obs.dta rdnet_mock_nodes.dta rdnet_mock_edges.dta {
    copy "`base'/`f'" "`f'", replace
}
capture mkdir test
foreach f in rd2dnet_tests.do river_obs.dta river_segs.dta {
    copy "`base'/test/`f'" "test/`f'", replace
}
di as txt "fetched rd2dnet Stata files at `sha' into `c(pwd)'"
