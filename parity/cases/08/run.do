* parity case 8 (grid)
rdn_par_load "parity/cases/08"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedregion(RG)
rdn_par_coords "parity/cases/08" 8
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc2) bwselect(imserd) fuzzy(D)
rdn_par_res "parity/cases/08/py_res_1.csv" 8 rdnet1
rdn_par_agg "parity/cases/08/py_agg_1.csv" 8 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) h(416.1 398.8) fuzzy(D)
rdn_par_res "parity/cases/08/py_res_2.csv" 8 rdnet2
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc1) h(443.6) fuzzy(D)
rdn_par_res "parity/cases/08/py_res_3.csv" 8 rdnet3
rdnetdensity, normalize(euclidean) bwselect(common)
rdn_par_dens "parity/cases/08/py_dens_1.csv" 8 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/08/py_bw_1.csv" 8 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/08/py_manip_1.csv" 8 manipulation
