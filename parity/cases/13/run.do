* parity case 13 (grid)
rdn_par_load "parity/cases/13"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedregion(RG)
rdn_par_coords "parity/cases/13" 13
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) h(409.9 285.9) fuzzy(D)
rdn_par_res "parity/cases/13/py_res_1.csv" 13 rdnet1
rdn_par_agg "parity/cases/13/py_agg_1.csv" 13 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) bwselect(imsetwo) vce(conley 154.1)
rdn_par_res "parity/cases/13/py_res_2.csv" 13 rdnet2
rdnet y_out, kernel(epanechnikov) p(1) seed(0) bwselect(cerrd) vce(conley 114.4)
rdn_par_res "parity/cases/13/py_res_3.csv" 13 rdnet3
rdnetdensity, normalize(network) bwselect(common)
rdn_par_dens "parity/cases/13/py_dens_1.csv" 13 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/13/py_bw_1.csv" 13 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/13/py_manip_1.csv" 13 manipulation
