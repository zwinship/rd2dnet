* parity case 23 (river)
rdn_par_load "parity/cases/23"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) treatedside(right)
rdn_par_coords "parity/cases/23" 23
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc3) h(293.1 378.8) covs(z1)
rdn_par_res "parity/cases/23/py_res_1.csv" 23 rdnet1
rdn_par_agg "parity/cases/23/py_agg_1.csv" 23 aggregates
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc1) h(440.5 343.5) fuzzy(D)
rdn_par_res "parity/cases/23/py_res_2.csv" 23 rdnet2
rdnet y_out, kernel(epanechnikov) p(1) seed(0) bwselect(imsetwo) vce(cluster cl)
rdn_par_res "parity/cases/23/py_res_3.csv" 23 rdnet3
rdnetdensity, normalize(network) bwselect(each)
rdn_par_dens "parity/cases/23/py_dens_1.csv" 23 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/23/py_bw_1.csv" 23 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/23/py_manip_1.csv" 23 manipulation
