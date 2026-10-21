* parity case 34 (ring)
rdn_par_load "parity/cases/34"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) treatedregion(RG)
rdn_par_coords "parity/cases/34" 34
rdnet y_out, kernel(epanechnikov) p(1) seed(0) h(387.4) vce(conley 114.6)
rdn_par_res "parity/cases/34/py_res_1.csv" 34 rdnet1
rdn_par_agg "parity/cases/34/py_agg_1.csv" 34 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) bwselect(mserd)
rdn_par_res "parity/cases/34/py_res_2.csv" 34 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) h(254.1) vce(cluster cl)
rdn_par_res "parity/cases/34/py_res_3.csv" 34 rdnet3
rdnetdensity, normalize(none) bwselect(each)
rdn_par_dens "parity/cases/34/py_dens_1.csv" 34 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/34/py_bw_1.csv" 34 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/34/py_manip_1.csv" 34 manipulation
