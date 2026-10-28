* parity case 29 (grid)
rdn_par_load "parity/cases/29"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) features(lanes) treatedregion(RG)
rdn_par_coords "parity/cases/29" 29
rdnet y_out [aw=w], kernel(triangular) p(1) seed(0) vce(hc2) bwselect(cerrd)
rdn_par_res "parity/cases/29/py_res_1.csv" 29 rdnet1
rdn_par_agg "parity/cases/29/py_agg_1.csv" 29 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc2) h(261.7)
rdn_par_res "parity/cases/29/py_res_2.csv" 29 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) bwselect(cerrd) vce(conley 153.8)
rdn_par_res "parity/cases/29/py_res_3.csv" 29 rdnet3
rdnetdensity, normalize(none) bwselect(common)
rdn_par_dens "parity/cases/29/py_dens_1.csv" 29 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/29/py_bw_1.csv" 29 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/29/py_manip_1.csv" 29 manipulation
rdn_par_bins "parity/cases/29" 29 350.0
