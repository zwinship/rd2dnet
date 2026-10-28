* parity case 26 (ring)
rdn_par_load "parity/cases/26"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) features(lanes) treatedregion(RG) classvar(cls) penalty(171.7) connectors(c)
rdn_par_coords "parity/cases/26" 26
rdnet y_out [aw=w], kernel(epanechnikov) p(1) seed(0) vce(hc0) h(344.9)
rdn_par_res "parity/cases/26/py_res_1.csv" 26 rdnet1
rdn_par_agg "parity/cases/26/py_agg_1.csv" 26 aggregates
rdnet y_out [aw=w], kernel(uniform) p(1) seed(0) vce(hc2) h(350.4)
rdn_par_res "parity/cases/26/py_res_2.csv" 26 rdnet2
rdnet y_out, kernel(epanechnikov) p(2) seed(0) vce(hc3) bwselect(imsetwo)
rdn_par_res "parity/cases/26/py_res_3.csv" 26 rdnet3
rdnetdensity, normalize(euclidean) bwselect(each)
rdn_par_dens "parity/cases/26/py_dens_1.csv" 26 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/26/py_bw_1.csv" 26 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/26/py_manip_1.csv" 26 manipulation
rdn_par_bins "parity/cases/26" 26 350.0
