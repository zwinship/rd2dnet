* parity case 2 (ring)
rdn_par_load "parity/cases/02"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) features(lanes) treatedregion(RG)
rdn_par_coords "parity/cases/02" 2
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc2) bwselect(imserd) fuzzy(D)
rdn_par_res "parity/cases/02/py_res_1.csv" 2 rdnet1
rdn_par_agg "parity/cases/02/py_agg_1.csv" 2 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) bwselect(cerrd) covs(z1)
rdn_par_res "parity/cases/02/py_res_2.csv" 2 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) h(300.1) vce(conley 157.2)
rdn_par_res "parity/cases/02/py_res_3.csv" 2 rdnet3
rdnetdensity, normalize(none) bwselect(each)
rdn_par_dens "parity/cases/02/py_dens_1.csv" 2 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/02/py_bw_1.csv" 2 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/02/py_manip_1.csv" 2 manipulation
rdn_par_bins "parity/cases/02" 2 350.0
