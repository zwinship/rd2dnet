* parity case 3 (river)
rdn_par_load "parity/cases/03"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) treatedside(right)
rdn_par_coords "parity/cases/03" 3
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc2) bwselect(certwo)
rdn_par_res "parity/cases/03/py_res_1.csv" 3 rdnet1
rdn_par_agg "parity/cases/03/py_agg_1.csv" 3 aggregates
rdnet y_out, kernel(epanechnikov) p(2) seed(0) vce(hc2) bwselect(certwo)
rdn_par_res "parity/cases/03/py_res_2.csv" 3 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) bwselect(certwo) vce(cluster cl)
rdn_par_res "parity/cases/03/py_res_3.csv" 3 rdnet3
rdnetdensity, normalize(euclidean) bwselect(common)
rdn_par_dens "parity/cases/03/py_dens_1.csv" 3 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/03/py_bw_1.csv" 3 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/03/py_manip_1.csv" 3 manipulation
