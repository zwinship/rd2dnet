* rd2dnet (Stata) test suite
*
* Run from this directory:   cd stata/test   then   do rd2dnet_tests.do
*
* Reference values come from the Python package rd2dnet 0.2.0 and the R
* package (which matches Python to about 1e-12 on this data and uses the same
* tie-breaking rule as the Stata port). Every check prints "ok" or "FAIL",
* errors inside a block are caught, and a summary is printed at the end.

version 16
clear all
set more off
adopath ++ "`c(pwd)'/.."
capture mata: mata drop rdn_*()

global rdn_fails 0
global rdn_errors 0
* failure labels in numbered globals rdn_fail1, rdn_fail2, ... and rdn_err1, ...

capture program drop rdn_expect
program define rdn_expect
    args got want tol label
    local g = `got'
    local d = abs(`g' - (`want')) / max(1, abs(`want'))
    if `g' >= . | !(`d' <= `tol') {
        di as err "FAIL  `label': got " %21.15g `g' "  want " %21.15g `want' "  (rel diff " %9.2e `d' ")"
        global rdn_fails = $rdn_fails + 1
        global rdn_fail$rdn_fails `"`label' (got `: display %12.0g `g'')"'
    }
    else di as txt "ok    `label'"
end

capture program drop rdn_expect_rows
program define rdn_expect_rows
    * compares matrix `mat' column `col' rows 3..6 with four reference values
    args mat col refs tol label
    forvalues k = 1/4 {
        local j = `k' + 2
        local want : word `k' of `refs'
        rdn_expect `mat'[`j',`col'] `want' `tol' "`label' [eval`j']"
    }
end

capture program drop rdn_block
program define rdn_block
    * runs a block (a program name) and counts errors instead of stopping
    args name
    di as txt _n "{hline 70}" _n "block: `name'" _n "{hline 70}"
    capture noisily `name'
    if _rc {
        di as err "ERROR in block `name' (rc = " _rc ")"
        global rdn_errors = $rdn_errors + 1
        global rdn_err$rdn_errors `"`name' (rc `=_rc')"'
    }
end

* --------------------------------------------------------------------------
capture program drop t_setup
program define t_setup
    rdnetmock, clear
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) ///
        treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
    rdn_expect r(N) 3000 0 "N"
    rdn_expect r(N_eval) 21 0 "N_eval"
    rdn_expect r(N_masked) 7 0 "N_masked"
    rdn_expect r(N_nodes) 450 0 "nodes (nothing to split)"
    matrix E = r(eval)
    local arcs 100 240 380 520
    forvalues j = 1/4 {
        rdn_expect E[`j',2] `: word `j' of `arcs'' 1e-10 "eval arc `j'"
    }
    local s 829.087 1649.546 -1353.994 -1783.17 1648.397 1653.963
    local a 1900 400 2400 2150 2400 2650
    forvalues i = 1/6 {
        rdn_expect _rdn_s[`i'] `: word `i' of `s'' 1e-10 "s[`i']"
        rdn_expect _rdn_arc[`i'] `: word `i' of `a'' 1e-10 "crossing arc[`i']"
    }
    quietly count if (_rdn_s >= 0) != _rdn_treated & _rdn_s < .
    rdn_expect r(N) 0 0 "sign of s matches treatment"
    generate double lnfar = ln(far)
    generate double D = (far > 3) * 0.5 + _rdn_treated * 0.4
end

capture program drop t_hc1
program define t_hc1
    rdnet outcome, h(600) vce(hc1) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 8  "2.2879700352572083 2.226131870704858 2.2374576858442436 2.216318573301832" 1e-8 "hc1 tau"
    rdn_expect_rows R 9  "0.6423910522468664 0.519468183514787 0.3817868050575414 0.450498262570776" 1e-8 "hc1 se"
    rdn_expect_rows R 10 "2.5630875033719818 2.671419704490127 2.338716529414876 2.3138693971233018" 1e-8 "hc1 tau_bc"
    rdn_expect_rows R 11 "0.9542401673002789 0.6640019263539401 0.6604308274435798 0.7795185876927633" 1e-8 "hc1 se_rb"
    matrix V = e(V)
    rdn_expect sqrt(V[1,1]) R[3,11] 1e-10 "e(V) diagonal = se_rb^2"
    rdn_expect R[9,15] 1 0 "masked point has status 1"
    rdn_expect (e(crit)>1.96) 1 0 "uniform critical value > 1.96"
    rdn_expect (R[3,17]<=R[3,12]+1e-12) 1 0 "band contains the pointwise CI"
    estat wbate, weights(equal)
    rdn_expect r(tau_bc) 1.13552840803952 1e-8 "WBATE tau_bc (equal weights)"
    rdn_expect r(se_rb) 0.35252971060533 1e-8 "WBATE se"
    estat lbate
    rdn_expect (r(ci_lo)<=r(tau_bc)&r(tau_bc)<=r(ci_hi)) 1 0 "LBATE interval contains max tau_bc"
    estat lbate, min
    rdn_expect r(tau_bc) -0.09516753502656217 1e-8 "LBATE min tau_bc"
    estat lbate, absmax
    rdn_expect r(tau_bc) 2.6714197044901438 1e-8 "LBATE absmax tau_bc"
    * equal WBATE weights per distinct place: the two highway-node points share one
    * (compared with explicit weights; at this bandwidth a few points equidistant
    * from two boundary nodes enter, and Python breaks those ties differently)
    rdnet outcome, h(2000) vce(hc1) seed(0)
    estat wbate, weights(equal)
    local wpl = r(tau_bc)
    matrix RW = e(results)
    matrix WW = J(rowsof(RW), 1, 0)
    forvalues j = 1/`=rowsof(RW)' {
        if RW[`j', 15] == 0 & `j' != 2 matrix WW[`j', 1] = 1
    }
    estat wbate, weights(WW)
    rdn_expect `wpl' r(tau_bc) 1e-12 "WBATE equal weights per place"
    lincom eval3 - eval4
end

capture program drop t_vce
program define t_vce
    rdnet outcome, h(600) vce(conley 150) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 11 "0.8979338736065969 0.6601068377141024 0.5275681466518771 0.6208458892927217" 1e-8 "conley network se_rb"
    rdnet outcome, h(600) vce(conley 150 euclidean) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 11 "0.748630568292863 0.602154181109756 0.449459055690671 0.491291224104242" 1e-8 "conley euclidean se_rb"
    rdnet outcome, h(600) vce(cluster block) cr1 seed(0)
    matrix R = e(results)
    rdn_expect_rows R 11 "0.3221750701567485 0.4308191219209627 0.45480015406506874 0.45409127249969544" 1e-8 "cluster CR1 se_rb"
    rdnet outcome, h(600) vce(cluster block) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 11 "0.8561412149244739 0.8145564381252992 0.6044838146114162 0.5914486569282575" 1e-8 "cluster CR2 se_rb"
    rdn_expect_rows R 23 "2.1342578317476404 2.608646673715296 3.6268505432025187 2.954971141665643" 1e-8 "cluster CR2 df"
    rdn_expect_rows R 12 "-0.9072022435686362 -0.1553808739725122 0.5900578237239358 0.41528505568326657" 1e-8 "cluster CR2 ci_lo"
    rdnet outcome, h(600) vce(hc3) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 9 "0.722066730205256 0.575806269535848 0.403330350144766 0.481835465352973" 1e-8 "hc3 se"
    * a vanishing Conley cutoff reproduces HC0
    rdnet outcome, h(600) vce(hc0) seed(0)
    matrix A = e(results)
    rdnet outcome, h(600) vce(conley 1e-9) conleyresidscale(1) seed(0)
    matrix B = e(results)
    rdn_expect B[3,11] A[3,11] 1e-8 "conley(1e-9) = hc0"
end

capture program drop t_models
program define t_models
    rdnet outcome, h(600) pooled seed(0)
    matrix R = e(results)
    rdn_expect_rows R 8 "2.33711903774262 2.28628744363469 2.24378223855976 2.20960952447098" 1e-8 "pooled tau"
    rdnet outcome [aw=far], h(500 700) vce(hc1) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 8  "2.39052455698188 2.60204927186686 2.44069271559587 2.4799033928382" 1e-8 "weights tau"
    rdn_expect_rows R 11 "0.780143295867975 0.507737602420415 0.494117067302827 0.577625222741211" 1e-8 "weights se_rb"
    rdnet outcome, h(600) covs(lnfar) vce(hc1) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 8 "2.26211114648349 2.2266580452012 2.23727453004989 2.20988779475705" 1e-8 "covs tau"
    rdn_expect_rows R 9 "0.647342414736302 0.519161733126116 0.382049751911618 0.448214407246662" 1e-8 "covs se"
    rdnet outcome, h(700) fuzzy(D) vce(hc1) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 8  "7.22528138668698 5.73956721098894 5.87435849102914 5.24722077290505" 1e-8 "fuzzy tau"
    rdn_expect_rows R 20 "0.292335873427158 0.38559002377044 0.379865950971725 0.418673765634722" 1e-8 "fuzzy first stage"
    rdn_expect_rows R 21 "0.121656177554145 0.101418775335716 0.073703401650608 0.0757686261226911" 1e-8 "fuzzy first-stage se"
    rdnet outcome, h(800) p(2) deriv(1 0) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 8 "-0.000429057343545343 0.000998060571910066 -4.63876866083853e-05 0.000418190439912467" 1e-6 "kink tau"
    rdnet outcome, h(800) p(2) tangvec(1 0) seed(0)
    matrix T = e(results)
    rdn_expect T[3,8] R[3,8] 1e-8 "tangvec(1 0) = deriv(1 0)"
end

capture program drop t_bandwidth
program define t_bandwidth
    rdnet outcome, vce(hc1) seed(0)
    matrix R = e(results)
    rdn_expect_rows R 4  "645.51930767128 730.659206427927 560.790448573427 821.495489610758" 1e-6 "auto h"
    rdn_expect_rows R 10 "2.54273718655865 2.42473790752563 2.41529140108223 2.21032099845553" 1e-6 "auto tau_bc"
    * eval 1 (the isolated highway node) is estimated; the crossings within its wide
    * bandwidth identify the linear term in o, so the basis is (s, o), as in Python
    rdn_expect R[1,15] 0 0 "eval 1 estimated"
    rdn_expect R[1,22] 1 0 "eval 1 basis = (s, o)"
    rdbwnet outcome, bwselect(msetwo) vce(hc1)
    matrix B = r(bws)
    rdn_expect_rows B 21 "1010.98975735317 750.620801408354 584.100691648631 866.851382283764" 1e-6 "msetwo h0"
    rdn_expect_rows B 22 "804.073938664862 849.151186796829 685.342258189746 879.02354935806" 1e-6 "msetwo h1"
    foreach sel in imserd imsetwo cerrd certwo {
        rdbwnet outcome, bwselect(`sel')
    }
    rdbwnet outcome, covs(lnfar) fuzzy(D) cluster(block)
    rdnet outcome, bwselect(msetwo) nocbands
end

capture program drop t_density
program define t_density
    * default: straight-line distance to the boundary (reference values from Python)
    rdnetdensity
    rdn_expect r(diff_bc) -1.16812813223094e-05 1e-6 "density (euclidean) diff_bc"
    rdn_expect r(se) 0.00010308279614672041 1e-6 "density (euclidean) se"
    rdn_expect r(h_left) 392.21811115656584 1e-6 "density (euclidean) h"
    rdnetdensity, normalize(network)
    rdn_expect r(diff_bc) 3.370389666344076e-06 1e-6 "density diff_bc"
    rdn_expect r(se) 3.200746336495476e-06 1e-6 "density se"
    rdn_expect r(h_left) 12249.64743015767 1e-6 "density h"
    rdnetdensity, normalize(none)
    rdn_expect r(diff_bc) -5.152889755147277e-05 1e-6 "density (none) diff_bc"
    rdn_expect r(se) 9.395949397424248e-05 1e-6 "density (none) se"
    rdn_expect r(h_left) 326.7653651609291 1e-6 "density (none) h"
    rdnetdensity, normalize(network) bwselect(each)
    rdn_expect r(diff_bc) 3.24605351338202e-06 1e-6 "density (each) diff_bc"
    rdn_expect r(h_left) 12409.2089484406 1e-6 "density (each) h_left"
    rdn_expect r(h_right) 12249.6474300536 1e-6 "density (each) h_right"
    rdnetdensity, normalize(network) evalindex(7) window(600) h(400)
    rdnetdensity, evalindex(7) window(600) h(400)
    * manipulation test: same clusters and draws as Python (reference values from Python)
    rdnetmanip, q(50 200) reps(2000)
    matrix MP = r(path)
    rdn_expect MP[1,2] 1.1666666666666667 1e-12 "manipulation S (q = 50)"
    rdn_expect MP[1,3] 0.2375 1e-12 "manipulation p (q = 50)"
    rdn_expect MP[2,3] 0.689 1e-12 "manipulation p (q = 200)"
    rdn_expect MP[2,5] 100 0 "manipulation n_treated (q = 200)"
end

capture program drop t_metrics
program define t_metrics
    rdnetmock, clear
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) cost(travel_time) boundary(rdn_line) ///
        treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
    local s 80.508528 111.763714909091 -121.85946 -183.33648 150.969168 156.090672
    forvalues i = 1/6 {
        rdn_expect _rdn_s[`i'] `: word `i' of `s'' 1e-9 "travel time s[`i']"
    }
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) classvar(highway) penalty(500) ///
        connectors(motorway_link) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
    local s 1329.087 2149.546 -1353.994 -2283.17 2148.397 2153.963
    forvalues i = 1/6 {
        rdn_expect _rdn_s[`i'] `: word `i' of `s'' 1e-9 "transition s[`i']"
    }
    matrix P = (0, 1000 \ 10, 0)
    matrix rownames P = motorway residential
    matrix colnames P = motorway residential
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) classvar(highway) penaltymatrix(P) ///
        connectors(motorway_link) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)
end

capture program drop t_options
program define t_options
    rdnetmock, clear
    * boundary given by side instead of region; masks by arc interval; custom evaluation points
    matrix A = (1150, 2150)
    matrix EP = (500 \ 900 \ 1300 \ 2300)
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedside(right) ///
        maskarcs(A) evalpoints(EP)
    rdn_expect r(N_eval) 4 0 "custom evaluation points"
    rdnet outcome, h(600)
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) snapnode
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) ///
        directed direction(both)
    rdnetset x y if y < 2500, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region)
    rdnet outcome, h(600) kernel(epanechnikov)
    rdnet outcome, h(600) kernel(uniform) vce(hc2)
    rdnet outcome, h(600) kernel(gaussian)
    * a boundary that cuts streets between nodes: edges are split
    matrix L2 = (-400, 1640 \ 2600, 1640)
    matrix R2 = (-100000, -100000 \ 100000, -100000 \ 100000, 1640 \ -100000, 1640)
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(L2) treatedregion(R2)
    rdn_expect (r(N_crossings)>0) 1 0 "edges split at a mid-block boundary"
    rdnet outcome, h(600)
end

capture program drop t_entry
program define t_entry
    * a network given as segments with end coordinates (no node frame)
    rdnetmock, clear
    capture frame drop segs
    frame copy rdn_edges segs
    frame segs {
        frlink m:1 from, frame(rdn_nodes id)
        frget x1 = x y1 = y, from(rdn_nodes)
        drop rdn_nodes
        frlink m:1 to, frame(rdn_nodes id)
        frget x2 = x y2 = y, from(rdn_nodes)
        drop rdn_nodes from to
    }
    rdnetset x y, edges(segs) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
    rdn_expect r(N_nodes) 450 0 "segments: nodes joined at shared ends"
    local s 829.087 1649.546 -1353.994 -1783.17 1648.397 1653.963
    forvalues i = 1/6 {
        rdn_expect _rdn_s[`i'] `: word `i' of `s'' 1e-10 "segments: s[`i']"
    }
    * rd2d-style entry: boundary through the evaluation points, side from treated()
    matrix B = (100, 1600 \ 1000, 1600 \ 2000, 1600 \ 2900, 1600)
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) evalpoints(B) treated(treated)
    rdn_expect r(N_eval) 4 0 "rd2d-style: one evaluation point per row of B"
    quietly count if _rdn_treated != treated
    rdn_expect r(N) 0 0 "rd2d-style: _rdn_treated = treated"
    rdnet outcome, h(900)
    * density test with repeated measurements
    rdnetmock, clear
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)
    rdnetdensity
    local p1 = r(pvalue)
    expand 3
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)
    rdnetdensity, unique
    rdn_expect r(pvalue) `p1' 1e-8 "density: unique locations = one row per location"
end

capture program drop t_river
program define t_river
    * three parallel rivers crossing the border: no variation along the boundary
    * near each crossing, so the fit is in s alone (reference values from Python)
    use river_obs, clear
    capture frame drop rivers
    frame create rivers
    frame rivers: use river_segs, clear
    matrix L = (0, -4000 \ 0, 4000)
    matrix B = (0, -3000 \ 0, 0 \ 0, 3000)
    rdnetset x y, edges(rivers) boundary(L) treatedside(right) evalpoints(B)
    local s -538.52 1974.13 102.73 -4022.22
    forvalues i = 1/4 {
        rdn_expect _rdn_s[`i'] `: word `i' of `s'' 1e-9 "river s[`i']"
    }
    rdnet y_out, h(1500) vce(hc1) seed(0)
    matrix R = e(results)
    local tau -1.3986516707699939 -1.2344053207051944 -1.4911025736082353
    local tbc -1.23500711932145 -0.9986365254587386 -1.2862431836128785
    local srb 0.2984688242880504 0.29074350841404983 0.20534562299776382
    forvalues j = 1/3 {
        rdn_expect R[`j',15] 0 0 "river status ok [eval`j']"
        rdn_expect R[`j',22] 2 0 "river basis = s alone [eval`j']"
        rdn_expect R[`j',8] `: word `j' of `tau'' 1e-8 "river tau [eval`j']"
        rdn_expect R[`j',10] `: word `j' of `tbc'' 1e-8 "river tau_bc [eval`j']"
        rdn_expect R[`j',11] `: word `j' of `srb'' 1e-8 "river se_rb [eval`j']"
    }
    rdbwnet y_out
end

capture program drop t_floor
program define t_floor
    * bandwidths at the selector's floor (the distance to the bwcheck-th nearest
    * observation) must leave that observation inside the kernel, so points
    * are estimable rather than "insufficient"
    rdnetmock, clear
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) ///
        treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
    rdbwnet outcome
    matrix W = r(bws)
    local J = rowsof(W)
    matrix HF = J(`J', 2, 1000)
    forvalues j = 1/`J' {
        if W[`j', 3] == 0 {
            local f = max(W[`j', 9], W[`j', 10])
            matrix HF[`j', 1] = `f'
            matrix HF[`j', 2] = `f'
        }
    }
    rdnet outcome, hmatrix(HF) seed(0)
    matrix R = e(results)
    local bad 0
    forvalues j = 1/`J' {
        if W[`j', 3] == 0 & R[`j', 15] == 2 local ++bad
    }
    rdn_expect `bad' 0 0 "no point at the bandwidth floor is insufficient"
    * two numbers in h() are (h_control h_treated), whatever the number of points
    matrix EP = (600 \ 2600)
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) evalpoints(EP)
    rdnet outcome, h(800 1000)
    matrix R = e(results)
    rdn_expect R[2,4] 800 0 "h(800 1000): h0 = 800"
    rdn_expect R[2,5] 1000 0 "h(800 1000): h1 = 1000"
end

capture program drop t_pointnodes
program define t_pointnodes
    * a location on a one-way edge is reached in the edge's direction
    * (X on the boundary; the one-way edge N -> M is entered at N via X -> R -> N)
    capture frame drop pn_nodes
    capture frame drop pn_edges
    frame create pn_nodes
    frame pn_nodes {
        set obs 7
        gen id = _n
        gen x = 0
        gen y = 0
        replace y = -100 in 2
        replace y = -200 in 3
        replace x = 300 in 4
        replace y = -100 in 4
        replace y = 100 in 5
        replace x = -100 in 6
        replace x = 100 in 7
    }
    frame create pn_edges
    frame pn_edges {
        set obs 7
        gen from = .
        gen to = .
        local k 0
        foreach e in "1 2" "3 2" "1 4" "4 3" "1 5" "6 1" "1 7" {
            local ++k
            replace from = real(word("`e'", 1)) in `k'
            replace to = real(word("`e'", 2)) in `k'
        }
        gen double length = 100
        replace length = sqrt(300^2 + 100^2) if (from == 1 & to == 4) | (from == 4 & to == 3)
    }
    clear
    set obs 2
    gen double x = 0
    gen double y = cond(_n == 1, -110, -150)
    matrix PL = (-100, 0 \ 100, 0)
    matrix PE = (100)
    rdnetset x y, nodes(pn_nodes) edges(pn_edges) boundary(PL) treatedside(left) evalpoints(PE) ///
        directed direction(out)
    local w1 = -(2 * sqrt(300^2 + 100^2) + 90)
    local w2 = -(2 * sqrt(300^2 + 100^2) + 50)
    rdn_expect _rdn_s[1] `w1' 1e-12 "one-way edge: s[1]"
    rdn_expect _rdn_s[2] `w2' 1e-12 "one-way edge: s[2]"
    * direction(flow): the shorter directed path either way; from the points on
    * N -> M only the downstream end M is reachable, and M is a dead end, so the
    * path is the one from the boundary (as with direction(out))
    rdnetset x y, nodes(pn_nodes) edges(pn_edges) boundary(PL) treatedside(left) evalpoints(PE) ///
        directed direction(flow)
    rdn_expect _rdn_s[1] `w1' 1e-12 "flow: s[1]"
    rdn_expect _rdn_s[2] `w2' 1e-12 "flow: s[2]"
    * the transition penalty applies to the edge a location sits on: points on the
    * primary edge P-Q reach the boundary along primary streets via Z, not via X
    frame pn_nodes {
        clear
        set obs 7
        gen id = _n
        gen x = cond(_n >= 5, 300, 0)
        gen y = 0
        replace y = -200 in 1
        replace y = -100 in 2
        replace y = 100 in 4
        replace y = -100 in 5
        replace y = 100 in 7
    }
    frame pn_edges {
        clear
        set obs 7
        gen from = .
        gen to = .
        gen str12 highway = "primary"
        local k 0
        foreach e in "1 2" "2 3" "3 4" "2 5" "5 6" "6 7" "3 6" {
            local ++k
            replace from = real(word("`e'", 1)) in `k'
            replace to = real(word("`e'", 2)) in `k'
        }
        replace highway = "residential" in 2/3
        gen double length = cond(inlist(_n, 4, 7), 300, 100)
    }
    clear
    set obs 3
    gen double x = 0
    gen double y = -150
    replace y = -199.9 in 2
    replace y = -200 in 3
    matrix PL = (-50, 0 \ 350, 0)
    matrix PE = (50)
    rdnetset x y, nodes(pn_nodes) edges(pn_edges) boundary(PL) treatedside(left) evalpoints(PE) ///
        classvar(highway) penalty(1000)
    rdn_expect _rdn_s[1] -450 1e-12 "transition, own edge class: s[1]"
    rdn_expect _rdn_s[2] -499.9 1e-12 "transition, own edge class: s[2]"
    rdn_expect _rdn_s[3] -500 1e-12 "transition, own edge class: s[3]"
end

capture program drop t_round2
program define t_round2
    * transition penalties with direction(in) follow the order of travel:
    * L -> A residential, A -> B primary (crossing the boundary at x = 150)
    capture frame drop r2_nodes
    capture frame drop r2_edges
    frame create r2_nodes
    frame r2_nodes {
        set obs 4
        gen id = _n
        gen x = 100 * (_n - 1)
        gen y = 0
    }
    frame create r2_edges
    frame r2_edges {
        set obs 4
        gen from = .
        gen to = .
        local k 0
        foreach e in "1 2" "2 3" "4 3" "3 4" {
            local ++k
            replace from = real(word("`e'", 1)) in `k'
            replace to = real(word("`e'", 2)) in `k'
        }
        gen str12 highway = cond(_n == 1, "residential", "primary")
        gen double length = 100
    }
    matrix P = (0, 50 \ 0, 0)
    matrix rownames P = residential primary
    matrix colnames P = residential primary
    clear
    set obs 2
    gen double x = cond(_n == 1, 10, 250)
    gen double y = 0
    matrix PL = (150, -100 \ 150, 100)
    matrix PE = (100)
    foreach d in in flow {
        rdnetset x y, nodes(r2_nodes) edges(r2_edges) boundary(PL) treatedside(right) evalpoints(PE) ///
            directed direction(`d') classvar(highway) penaltymatrix(P)
        rdn_expect _rdn_s[1] -190 1e-12 "transition, direction(`d'): s[1]"
    }
    * fuzzy points with a first stage near zero get status 5 (weak first stage)
    rdnetmock, clear
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)
    set seed 1
    gen byte d = runiform() < 0.5
    rdnet outcome, h(600) fuzzy(d)
    mata: st_numscalar("r2_weak", sum(st_matrix("e(results)")[., 15] :== 5))
    mata: st_numscalar("r2_bad", sum((st_matrix("e(results)")[., 15] :== 0) :& (abs(st_matrix("e(results)")[., 20]) :< 2 :* st_matrix("e(results)")[., 21])))
    rdn_expect (scalar(r2_weak)>0) 1 0 "weak first stage flagged"
    rdn_expect scalar(r2_bad) 0 0 "no ok point with |fs| < 2 fs_se"
    rcof "noisily rdnet outcome, h(600) tangvec(0 0)" == 198
    rcof "noisily rdnetmanip, kmax(0)" == 198
end

capture program drop t_round3
program define t_round3
    * an edge crossing the boundary twice (both ends on the control side) is split twice
    capture frame drop r3_nodes
    capture frame drop r3_edges
    frame create r3_nodes
    frame r3_nodes {
        set obs 4
        gen id = _n
        gen x = 20
        gen y = 45
        replace x = 80 in 2
        replace x = 0 in 3/4
        replace y = 0 in 3
        replace y = 100 in 4
    }
    frame create r3_edges
    frame r3_edges {
        set obs 3
        gen from = cond(_n == 2, 3, 1)
        gen to = cond(_n == 1, 2, cond(_n == 2, 4, 3))
        gen double length = cond(_n == 1, 60, cond(_n == 2, 100, sqrt(20^2 + 45^2)))
    }
    clear
    set obs 2
    gen double x = cond(_n == 1, 45, 30)
    gen double y = 46
    matrix NL = (-10, 50 \ 40, 50 \ 40, 40 \ 60, 40 \ 60, 50 \ 110, 50)
    matrix PE = (55)
    rdnetset x y, nodes(r3_nodes) edges(r3_edges) boundary(NL) treatedside(left) evalpoints(PE)
    rdn_expect _rdn_s[1] 6 1e-9 "double crossing: s[1]"
    rdn_expect _rdn_s[2] -11 1e-9 "double crossing: s[2]"
    * evaluation points on a ring that follows streets stay where they were asked for
    frame r3_nodes {
        clear
        set obs 81
        gen id = _n
        gen x = 25 * mod(_n - 1, 9)
        gen y = 25 * floor((_n - 1) / 9)
    }
    frame r3_edges {
        clear
        set obs 144
        gen from = .
        gen to = .
        local k 0
        forvalues i = 1/81 {
            if mod(`i' - 1, 9) < 8 {
                local ++k
                replace from = `i' in `k'
                replace to = `i' + 1 in `k'
            }
            if `i' <= 72 {
                local ++k
                replace from = `i' in `k'
                replace to = `i' + 9 in `k'
            }
        }
        gen double length = 25
    }
    clear
    set seed 1
    set obs 500
    gen double x = runiform(5, 195)
    gen double y = runiform(5, 195)
    matrix RG = (62.5, 50 \ 150, 50 \ 150, 150 \ 50, 150 \ 50, 50 \ 62.5, 50)
    matrix PE = (237.5, 137.5, 337.5, 77.5, 397.5)
    rdnetset x y, nodes(r3_nodes) edges(r3_edges) boundary(RG) treatedregion(RG) evalpoints(PE) tol(0.5)
    matrix E = r(eval)
    forvalues j = 1/5 {
        rdn_expect E[`j',3] PE[1,`j'] 1e-9 "ring eval point `j' stays in place"
    }
end

* --------------------------------------------------------------------------
* Nearest-crossing bins (rdnetbins, rdnetfe, estat wbate density|gls).
* Reference values from Python rd2dnet (bin_coordinates on the mock city).
capture program drop t_bins
program define t_bins
    rdnetmock, clear
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) ///
        treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
    rdnetbins
    rdn_expect r(N_bins) 11 0 "bins: one per crossing reached"
    rdn_expect r(rank) 3 0 "bins: rank of the topology features (no triangles on a grid)"
    rdn_expect r(h_S) 3.953410654928521 1e-12 "bins: rule-of-thumb h_S"
    matrix B = r(bins)
    rdn_expect B[2,3] 214 0 "bins: n0 of bin 2"
    rdn_expect B[2,4] 187 0 "bins: n1 of bin 2"
    rdn_expect B[3,5] 6 0 "bins: branches of bin 3"
    rdn_expect B[4,6] 1 0 "bins: bin 4 masked"
    local S 1.807995690261959 2.325535786188288 1.1542993822150858 1.6480573193335377 2.241610069418619 1.8462591826814174
    local bn 7 1 9 8 9 10
    forvalues i = 1/6 {
        rdn_expect _rdn_S[`i'] `: word `i' of `S'' 1e-9 "bins: S[`i']"
        rdn_expect _rdn_bin[`i'] `: word `i' of `bn'' 0 "bins: bin[`i']"
    }
    rdnet outcome, h(700) minobs(10)
    rdn_expect ("`e(design)'" == "bins") 1 0 "bins: e(design)"
    matrix R = e(results)
    local rows 1 2 3 9 10 11
    local tbc 2.6059830764548826 1.7448275742516373 2.6104339718301404 0.7479398969433454 -0.44696316640305556 0.2975564306702405
    local srb 1.1296253005948182 0.8552830518211703 1.5565437316143036 0.9590629182735356 1.0286060688798866 1.4690102717534774
    forvalues k = 1/6 {
        local j : word `k' of `rows'
        rdn_expect R[`j',10] `: word `k' of `tbc'' 1e-7 "bins: tau_bc bin `j'"
        rdn_expect R[`j',11] `: word `k' of `srb'' 1e-7 "bins: se_rb bin `j'"
    }
    rdn_expect R[5,15] 1 0 "bins: masked bin"
    estat wbate, weights(density)
    rdn_expect r(tau_bc) 1.2399885813481475 1e-7 "bins: density WBATE"
    rdn_expect r(se_rb) 0.4857654141757674 1e-7 "bins: density WBATE se"
    estat wbate, weights(gls)
    rdn_expect r(tau_bc) 1.1925819893045733 1e-7 "bins: GLS WBATE"
    rdn_expect r(se_rb) 0.4445901952843397 1e-7 "bins: GLS WBATE se"
    rdnetfe outcome, h(700)
    rdn_expect e(tau) 1.4329741053916103 1e-7 "rdnetfe: tau"
    rdn_expect e(se) 0.2693587968689204 1e-7 "rdnetfe: se"
    rdn_expect e(tau_bc) 1.1824789662250355 1e-7 "rdnetfe: tau_bc"
    rdn_expect e(se_rb) 0.48376984245544763 1e-7 "rdnetfe: se_rb"
    rdnetfe outcome, h(700) vce(cluster block)
    rdn_expect e(se_rb) 0.6089516068407178 1e-7 "rdnetfe: CR2 se_rb"
    rdn_expect e(df) 8.331857443318553 1e-7 "rdnetfe: CR2 degrees of freedom"
    * rdnetset brings back the default design
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) ///
        treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
    rdnet outcome, h(600)
    rdn_expect ("`e(design)'" == "location") 1 0 "default design after rdnetset"
end

* Branch slopes recover the jump exactly on a comb: vertical streets crossing
* y = 0 (treated below), each splitting 200 m from the boundary into a
* straight continuation and a spur (mirrors tests/test_bins.py).
capture program drop t_bins_comb
program define t_bins_comb
    capture frame drop cnodes
    capture frame drop cedges
    frame create cnodes str8 id double(x y)
    frame create cedges str8(from to) double(length lanes)
    forvalues k = 0/2 {
        local x = `k' * 2000
        foreach sg in 1 -1 {
            local t = cond(`sg' > 0, "n", "s")
            frame post cnodes ("`t'j`k'") (`x') (`sg' * 200)
            frame post cnodes ("`t'e`k'") (`x') (`sg' * 600)
            frame post cnodes ("`t'p`k'") (`x' + 300) (`sg' * 500)
            frame post cedges ("`t'j`k'") ("`t'e`k'") (400) (2)
            frame post cedges ("`t'j`k'") ("`t'p`k'") (sqrt(2) * 300) (1)
        }
        frame post cedges ("nj`k'") ("sj`k'") (400) (4)
    }
    forvalues k = 0/1 {
        local k1 = `k' + 1
        frame post cedges ("ne`k'") ("ne`k1'") (2000) (1)
        frame post cedges ("se`k'") ("se`k1'") (2000) (1)
    }
    clear
    quietly set obs 540
    generate long i = _n - 1
    generate int k = floor(i / 180)
    generate int sg = cond(mod(floor(i / 90), 2) == 0, 1, -1)
    generate int seg = mod(floor(i / 30), 3)
    generate double t = (mod(i, 30) + 0.5) / 30
    generate double x = k * 2000 + cond(seg == 2, 300 * t, 0)
    generate double y = sg * cond(seg == 0, 199 * t + 0.5, cond(seg == 1, 200 + 400 * t, 200 + 300 * t))
    matrix L = (-500, 0 \ 4500, 0)
    rdnetset x y, nodes(cnodes) edges(cedges) boundary(L) treatedside(right) features(lanes) neval(3)
    rdnetbins, notopology
    rdn_expect r(N_bins) 3 0 "comb: three bins"
    rdn_expect r(rank) 1 0 "comb: one feature"
    quietly count if _rdn_bin != k + 1
    rdn_expect r(N) 0 0 "comb: bins are the nearest crossings"
    quietly count if seg == 0 & _rdn_branch != 0
    rdn_expect r(N) 0 0 "comb: trunk before the junction"
    quietly count if seg > 0 & _rdn_branch == 0
    rdn_expect r(N) 0 0 "comb: branches beyond the junction"
    tempvar nb
    bysort _rdn_bin _rdn_treated seg: egen `nb' = nvals(_rdn_branch)
    quietly summarize `nb'
    rdn_expect r(max) 1 0 "comb: one branch per segment and side"
    frame cedges: quietly summarize lanes
    local sd = r(sd)
    generate double S_want = cond(seg == 0, 0, cond(seg == 1, 2, 3)) / `sd'
    generate double dS = abs(_rdn_S - S_want)
    quietly summarize dS
    rdn_expect r(max) 0 1e-10 "comb: S is the Mahalanobis distance of the lanes"
    egen br = group(_rdn_branch)
    generate double slope = cond(_rdn_branch == 0, 0.002, 0.001 * (br + 2) * cond(mod(br, 2), 1, -1))
    generate double yy = cond(_rdn_bin == 1, 5, cond(_rdn_bin == 2, -1, 0.5)) + _rdn_bin * _rdn_treated + slope * _rdn_s
    rdnetbins, notopology hs(inf)
    rdnet yy, h(1500) minobs(5)
    matrix R = e(results)
    forvalues j = 1/3 {
        rdn_expect R[`j',8] `j' 1e-8 "comb: tau bin `j' with branch slopes"
        rdn_expect R[`j',10] `j' 1e-8 "comb: tau_bc bin `j' with branch slopes"
    }
    rdnetbins, nosimilarity nobranches
    rdnet yy, h(1500) minobs(5)
    matrix R = e(results)
    local dev = max(abs(R[1,8] - 1), abs(R[2,8] - 2), abs(R[3,8] - 3))
    rdn_expect (`dev' > 0.05) 1 0 "comb: biased without branch slopes"
    * rdnetfe with p = 0 and no slopes: the IPW-weighted fixed effects average
    * the bin effects by kernel mass
    generate double y0 = cond(_rdn_bin == 1, 5, cond(_rdn_bin == 2, -1, 0.5)) + cond(_rdn_bin == 3, 4, _rdn_bin) * _rdn_treated
    rdnetfe y0, h(500) p(0) q(0) nobinslopes
    generate double W = max(1 - abs(_rdn_s) / 500, 0)
    quietly summarize W if _rdn_bin == 1
    local m1 = r(sum)
    quietly summarize W if _rdn_bin == 2
    local m2 = r(sum)
    quietly summarize W if _rdn_bin == 3
    local m3 = r(sum)
    rdn_expect e(tau) (`m1' + 2 * `m2' + 4 * `m3') / (`m1' + `m2' + `m3') 1e-10 "rdnetfe: mass-weighted average of the bin effects"
end

capture program drop t_errors
program define t_errors
    rdnetmock, clear
    rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)
    rcof "noisily rdnet outcome, vce(bogus)" == 198
    rcof "noisily rdnet outcome, deriv(2 0)" != 0
    rcof "noisily rdnet outcome, deriv(1 0) tangvec(1 0)" == 198
    rcof "noisily rdnet outcome, bwselect(bogus)" == 198
    rcof "noisily rdnet outcome, vce(conley -1)" == 198
    rcof "noisily rdnetdensity, evalindex(3)" == 198
    rcof "noisily estat wbate, weights(1 2 3)" != 0
    rcof `"noisily rdnetset x y, nodes(nope) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region)"' == 111
    rcof `"noisily rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line)"' == 198
    rcof `"noisily rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) treatedside(left)"' == 198
    rcof "noisily rdnet outcome, h(600) vce(cluster treated)" == 0
    matrix R = e(results)
    rdn_expect R[3,15] 6 0 "one cluster per side: status 6 (few clusters)"
end

rdn_block t_setup
rdn_block t_hc1
rdn_block t_vce
rdn_block t_models
rdn_block t_bandwidth
rdn_block t_density
rdn_block t_metrics
rdn_block t_options
rdn_block t_entry
rdn_block t_river
rdn_block t_floor
rdn_block t_pointnodes
rdn_block t_errors
rdn_block t_round2
rdn_block t_round3
rdn_block t_bins
rdn_block t_bins_comb

di _n "{hline 70}"
if $rdn_fails == 0 & $rdn_errors == 0 {
    di as res "ALL TESTS PASSED"
}
else {
    di as err "failed checks: $rdn_fails   blocks with errors: $rdn_errors"
    forvalues k = 1/$rdn_fails {
        di as err `"  FAIL ${rdn_fail`k'}"'
    }
    forvalues k = 1/$rdn_errors {
        di as err `"  ERROR ${rdn_err`k'}"'
    }
}
