* parity case 6 (ring)
rdn_par_load "parity/cases/06"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) treatedregion(RG)
rdn_par_coords "parity/cases/06" 6
rdnet y_out [aw=w], kernel(epanechnikov) p(1) seed(0) vce(hc0) bwselect(imsetwo)
rdn_par_res "parity/cases/06/py_res_1.csv" 6 rdnet1
rdn_par_agg "parity/cases/06/py_agg_1.csv" 6 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc2) h(272.0)
rdn_par_res "parity/cases/06/py_res_2.csv" 6 rdnet2
rdnet y_out, kernel(epanechnikov) p(2) seed(0) vce(hc1) bwselect(cerrd)
rdn_par_res "parity/cases/06/py_res_3.csv" 6 rdnet3
rdnetdensity, normalize(euclidean) bwselect(common)
rdn_par_dens "parity/cases/06/py_dens_1.csv" 6 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/06/py_bw_1.csv" 6 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/06/py_manip_1.csv" 6 manipulation
