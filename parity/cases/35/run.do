* parity case 35 (river)
rdn_par_load "parity/cases/35"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) treatedside(right)
rdn_par_coords "parity/cases/35" 35
rdnet y_out [aw=w], kernel(uniform) p(1) seed(0) vce(hc2) h(251.7)
rdn_par_res "parity/cases/35/py_res_1.csv" 35 rdnet1
rdn_par_agg "parity/cases/35/py_agg_1.csv" 35 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc3) h(292.8)
rdn_par_res "parity/cases/35/py_res_2.csv" 35 rdnet2
rdnet y_out [aw=w], kernel(uniform) p(1) seed(0) vce(hc3) bwselect(imsetwo)
rdn_par_res "parity/cases/35/py_res_3.csv" 35 rdnet3
rdnetdensity, normalize(none) bwselect(each)
rdn_par_dens "parity/cases/35/py_dens_1.csv" 35 density1
rdbwnet y_out, bwselect(mserd)
rdn_par_bw "parity/cases/35/py_bw_1.csv" 35 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/35/py_manip_1.csv" 35 manipulation
