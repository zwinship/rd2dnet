* parity case 18 (ring)
rdn_par_load "parity/cases/18"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) features(lanes) treatedregion(RG) cost(cost2)
rdn_par_coords "parity/cases/18" 18
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc3) bwselect(mserd)
rdn_par_res "parity/cases/18/py_res_1.csv" 18 rdnet1
rdn_par_agg "parity/cases/18/py_agg_1.csv" 18 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) h(420.4) vce(conley 109.0)
rdn_par_res "parity/cases/18/py_res_2.csv" 18 rdnet2
rdnet y_out, kernel(uniform) p(2) seed(0) vce(hc2) bwselect(msetwo)
rdn_par_res "parity/cases/18/py_res_3.csv" 18 rdnet3
rdnetdensity, normalize(none) bwselect(common)
rdn_par_dens "parity/cases/18/py_dens_1.csv" 18 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/18/py_bw_1.csv" 18 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/18/py_manip_1.csv" 18 manipulation
rdn_par_bins "parity/cases/18" 18 350.0
