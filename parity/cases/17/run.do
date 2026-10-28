* parity case 17 (grid)
rdn_par_load "parity/cases/17"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) features(lanes) treatedside(left) cost(cost2)
rdn_par_coords "parity/cases/17" 17
rdnet y_out, kernel(epanechnikov) p(1) seed(0) h(308.0) vce(conley 179.5)
rdn_par_res "parity/cases/17/py_res_1.csv" 17 rdnet1
rdn_par_agg "parity/cases/17/py_agg_1.csv" 17 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) bwselect(msetwo) vce(conley 94.9)
rdn_par_res "parity/cases/17/py_res_2.csv" 17 rdnet2
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc1) bwselect(imserd) covs(z1)
rdn_par_res "parity/cases/17/py_res_3.csv" 17 rdnet3
rdnetdensity, normalize(none) bwselect(common)
rdn_par_dens "parity/cases/17/py_dens_1.csv" 17 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/17/py_bw_1.csv" 17 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/17/py_manip_1.csv" 17 manipulation
rdn_par_bins "parity/cases/17" 17 350.0
