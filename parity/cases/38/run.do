* parity case 38 (ring)
rdn_par_load "parity/cases/38"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) features(lanes) treatedregion(RG) cost(cost2)
rdn_par_coords "parity/cases/38" 38
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc2) bwselect(imserd) fuzzy(D)
rdn_par_res "parity/cases/38/py_res_1.csv" 38 rdnet1
rdn_par_agg "parity/cases/38/py_agg_1.csv" 38 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc1) h(358.9) fuzzy(D)
rdn_par_res "parity/cases/38/py_res_2.csv" 38 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) h(429.7) vce(cluster cl)
rdn_par_res "parity/cases/38/py_res_3.csv" 38 rdnet3
rdnetdensity, normalize(none) bwselect(common)
rdn_par_dens "parity/cases/38/py_dens_1.csv" 38 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/38/py_bw_1.csv" 38 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/38/py_manip_1.csv" 38 manipulation
rdn_par_bins "parity/cases/38" 38 350.0
