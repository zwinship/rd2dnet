* parity case 22 (ring)
rdn_par_load "parity/cases/22"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) features(lanes) treatedregion(RG) cost(cost2)
rdn_par_coords "parity/cases/22" 22
rdnet y_out [aw=w], kernel(uniform) p(1) seed(0) vce(hc3) h(357.6)
rdn_par_res "parity/cases/22/py_res_1.csv" 22 rdnet1
rdn_par_agg "parity/cases/22/py_agg_1.csv" 22 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc1) h(251.3) deriv(1 0)
rdn_par_res "parity/cases/22/py_res_2.csv" 22 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc3) h(371.7)
rdn_par_res "parity/cases/22/py_res_3.csv" 22 rdnet3
rdnetdensity, normalize(none) bwselect(each)
rdn_par_dens "parity/cases/22/py_dens_1.csv" 22 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/22/py_bw_1.csv" 22 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/22/py_manip_1.csv" 22 manipulation
rdn_par_bins "parity/cases/22" 22 350.0
