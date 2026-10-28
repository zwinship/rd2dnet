* parity case 24 (grid)
rdn_par_load "parity/cases/24"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) features(lanes) treatedside(right) directed direction(out)
rdn_par_coords "parity/cases/24" 24
rdnet y_out, kernel(epanechnikov) p(1) seed(0) bwselect(imsetwo) vce(conley 141.7)
rdn_par_res "parity/cases/24/py_res_1.csv" 24 rdnet1
rdn_par_agg "parity/cases/24/py_agg_1.csv" 24 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc3) bwselect(mserd) deriv(1 0)
rdn_par_res "parity/cases/24/py_res_2.csv" 24 rdnet2
rdnet y_out, kernel(epanechnikov) p(2) seed(0) vce(hc0) h(421.6 318.8)
rdn_par_res "parity/cases/24/py_res_3.csv" 24 rdnet3
rdnetdensity, normalize(euclidean) bwselect(each)
rdn_par_dens "parity/cases/24/py_dens_1.csv" 24 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/24/py_bw_1.csv" 24 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/24/py_manip_1.csv" 24 manipulation
rdn_par_bins "parity/cases/24" 24 350.0
