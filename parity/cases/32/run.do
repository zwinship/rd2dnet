* parity case 32 (grid)
rdn_par_load "parity/cases/32"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) treatedregion(RG)
rdn_par_coords "parity/cases/32" 32
rdnet y_out, kernel(epanechnikov) p(1) seed(0) h(410.4) vce(cluster cl)
rdn_par_res "parity/cases/32/py_res_1.csv" 32 rdnet1
rdn_par_agg "parity/cases/32/py_agg_1.csv" 32 aggregates
rdnet y_out, kernel(uniform) p(2) seed(0) vce(hc1) h(417.6)
rdn_par_res "parity/cases/32/py_res_2.csv" 32 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) h(431.8) vce(cluster cl)
rdn_par_res "parity/cases/32/py_res_3.csv" 32 rdnet3
rdnetdensity, normalize(euclidean) bwselect(each)
rdn_par_dens "parity/cases/32/py_dens_1.csv" 32 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/32/py_bw_1.csv" 32 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/32/py_manip_1.csv" 32 manipulation
