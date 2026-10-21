* parity case 19 (river)
rdn_par_load "parity/cases/19"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) treatedside(right)
rdn_par_coords "parity/cases/19" 19
rdnet y_out, kernel(uniform) p(1) seed(0) bwselect(imserd) vce(conley 127.1)
rdn_par_res "parity/cases/19/py_res_1.csv" 19 rdnet1
rdn_par_agg "parity/cases/19/py_agg_1.csv" 19 aggregates
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc1) h(405.9)
rdn_par_res "parity/cases/19/py_res_2.csv" 19 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc1) h(428.3)
rdn_par_res "parity/cases/19/py_res_3.csv" 19 rdnet3
rdnetdensity, normalize(network) bwselect(each)
rdn_par_dens "parity/cases/19/py_dens_1.csv" 19 density1
rdbwnet y_out, bwselect(mserd)
rdn_par_bw "parity/cases/19/py_bw_1.csv" 19 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/19/py_manip_1.csv" 19 manipulation
