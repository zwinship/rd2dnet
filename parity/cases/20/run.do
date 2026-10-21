* parity case 20 (grid)
rdn_par_load "parity/cases/20"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) treatedregion(RG) classvar(cls) penalty(87.4) connectors(c)
rdn_par_coords "parity/cases/20" 20
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) h(262.5 447.0) fuzzy(D)
rdn_par_res "parity/cases/20/py_res_1.csv" 20 rdnet1
rdn_par_agg "parity/cases/20/py_agg_1.csv" 20 aggregates
rdnet y_out [aw=w], kernel(epanechnikov) p(1) seed(0) vce(hc0) h(291.7)
rdn_par_res "parity/cases/20/py_res_2.csv" 20 rdnet2
rdnet y_out [aw=w], kernel(epanechnikov) p(1) seed(0) vce(hc2) h(356.3 392.2)
rdn_par_res "parity/cases/20/py_res_3.csv" 20 rdnet3
rdnetdensity, normalize(none) bwselect(each)
rdn_par_dens "parity/cases/20/py_dens_1.csv" 20 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/20/py_bw_1.csv" 20 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/20/py_manip_1.csv" 20 manipulation
