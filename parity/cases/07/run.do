* parity case 7 (river)
rdn_par_load "parity/cases/07"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedside(right)
rdn_par_coords "parity/cases/07" 7
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc0) bwselect(mserd) covs(z1)
rdn_par_res "parity/cases/07/py_res_1.csv" 7 rdnet1
rdn_par_agg "parity/cases/07/py_agg_1.csv" 7 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc0) h(412.6 296.1)
rdn_par_res "parity/cases/07/py_res_2.csv" 7 rdnet2
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) h(306.8)
rdn_par_res "parity/cases/07/py_res_3.csv" 7 rdnet3
rdnetdensity, normalize(euclidean) bwselect(each)
rdn_par_dens "parity/cases/07/py_dens_1.csv" 7 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/07/py_bw_1.csv" 7 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/07/py_manip_1.csv" 7 manipulation
