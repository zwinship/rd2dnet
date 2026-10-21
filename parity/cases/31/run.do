* parity case 31 (river)
rdn_par_load "parity/cases/31"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedside(right)
rdn_par_coords "parity/cases/31" 31
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) h(284.5) fuzzy(D)
rdn_par_res "parity/cases/31/py_res_1.csv" 31 rdnet1
rdn_par_agg "parity/cases/31/py_agg_1.csv" 31 aggregates
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc0) h(302.2 266.0)
rdn_par_res "parity/cases/31/py_res_2.csv" 31 rdnet2
rdnet y_out [aw=w], kernel(epanechnikov) p(1) seed(0) vce(hc1) bwselect(mserd)
rdn_par_res "parity/cases/31/py_res_3.csv" 31 rdnet3
rdnetdensity, normalize(none) bwselect(common)
rdn_par_dens "parity/cases/31/py_dens_1.csv" 31 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/31/py_bw_1.csv" 31 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/31/py_manip_1.csv" 31 manipulation
