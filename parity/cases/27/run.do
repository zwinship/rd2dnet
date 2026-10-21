* parity case 27 (river)
rdn_par_load "parity/cases/27"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedside(right)
rdn_par_coords "parity/cases/27" 27
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc3) h(428.1)
rdn_par_res "parity/cases/27/py_res_1.csv" 27 rdnet1
rdn_par_agg "parity/cases/27/py_agg_1.csv" 27 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc1) bwselect(imsetwo) covs(z1)
rdn_par_res "parity/cases/27/py_res_2.csv" 27 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc0) h(362.5) covs(z1)
rdn_par_res "parity/cases/27/py_res_3.csv" 27 rdnet3
rdnetdensity, normalize(network) bwselect(each)
rdn_par_dens "parity/cases/27/py_dens_1.csv" 27 density1
rdbwnet y_out, bwselect(mserd)
rdn_par_bw "parity/cases/27/py_bw_1.csv" 27 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/27/py_manip_1.csv" 27 manipulation
