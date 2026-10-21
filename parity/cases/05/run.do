* parity case 5 (grid)
rdn_par_load "parity/cases/05"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedside(left)
rdn_par_coords "parity/cases/05" 5
rdnet y_out, kernel(triangular) p(1) seed(0) bwselect(imserd) vce(cluster cl)
rdn_par_res "parity/cases/05/py_res_1.csv" 5 rdnet1
rdn_par_agg "parity/cases/05/py_agg_1.csv" 5 aggregates
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc2) h(323.1 354.8) deriv(1 0)
rdn_par_res "parity/cases/05/py_res_2.csv" 5 rdnet2
rdnet y_out, kernel(epanechnikov) p(2) seed(0) vce(hc0) h(258.0 415.5)
rdn_par_res "parity/cases/05/py_res_3.csv" 5 rdnet3
rdnetdensity, normalize(network) bwselect(common)
rdn_par_dens "parity/cases/05/py_dens_1.csv" 5 density1
rdbwnet y_out, bwselect(mserd)
rdn_par_bw "parity/cases/05/py_bw_1.csv" 5 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/05/py_manip_1.csv" 5 manipulation
