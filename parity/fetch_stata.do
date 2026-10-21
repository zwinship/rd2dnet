* Download the parity cases and runner at a commit into ~/rd2dnet_parity
*   do https://raw.githubusercontent.com/zwinship/rd2dnet/<sha>/parity/fetch_stata.do <sha>
* then (with rd2dnet installed or on the ado-path):  do parity/run_stata.do
args sha
local home : env USERPROFILE
if "`home'" == "" local home : env HOME
cd "`home'"
capture mkdir rd2dnet_parity
cd rd2dnet_parity
capture mkdir parity
capture mkdir parity/cases
local base "https://raw.githubusercontent.com/zwinship/rd2dnet/`sha'/parity"
copy "`base'/run_stata.do" "parity/run_stata.do", replace
copy "`base'/cases/cases.txt" "parity/cases/cases.txt", replace
file open fh using "parity/cases/cases.txt", read text
file read fh line
while r(eof) == 0 {
    local cs = strtrim("`line'")
    capture mkdir "parity/cases/`cs'"
    foreach f in run.do nodes.csv edges.csv obs.csv boundary.csv py_coords.csv py_res_1.csv py_res_2.csv py_res_3.csv py_dens_1.csv py_bw_1.csv py_agg_1.csv py_manip_1.csv {
        quietly copy "`base'/cases/`cs'/`f'" "parity/cases/`cs'/`f'", replace
    }
    capture quietly copy "`base'/cases/`cs'/region.csv" "parity/cases/`cs'/region.csv", replace
    di as txt "case `cs' fetched"
    file read fh line
}
file close fh
