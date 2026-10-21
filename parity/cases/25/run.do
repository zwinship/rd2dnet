* parity case 25 (grid)
rdn_par_load "parity/cases/25"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedside(right)
rdn_par_coords "parity/cases/25" 25
rdnet y_out, kernel(triangular) p(1) seed(0) bwselect(imsetwo) vce(conley 136.8)
rdn_par_res "parity/cases/25/py_res_1.csv" 25 rdnet1
rdn_par_agg "parity/cases/25/py_agg_1.csv" 25 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc1) bwselect(imsetwo) fuzzy(D)
rdn_par_res "parity/cases/25/py_res_2.csv" 25 rdnet2
rdnet y_out, kernel(triangular) p(2) seed(0) vce(hc0) bwselect(mserd)
rdn_par_res "parity/cases/25/py_res_3.csv" 25 rdnet3
rdnetdensity, normalize(none) bwselect(common)
rdn_par_dens "parity/cases/25/py_dens_1.csv" 25 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/25/py_bw_1.csv" 25 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/25/py_manip_1.csv" 25 manipulation
