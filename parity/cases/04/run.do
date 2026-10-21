* parity case 4 (grid)
rdn_par_load "parity/cases/04"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedside(left)
rdn_par_coords "parity/cases/04" 4
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc2) h(368.2 274.4)
rdn_par_res "parity/cases/04/py_res_1.csv" 4 rdnet1
rdn_par_agg "parity/cases/04/py_agg_1.csv" 4 aggregates
rdnet y_out, kernel(uniform) p(1) seed(0) bwselect(cerrd) vce(cluster cl)
rdn_par_res "parity/cases/04/py_res_2.csv" 4 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) bwselect(msetwo) vce(cluster cl)
rdn_par_res "parity/cases/04/py_res_3.csv" 4 rdnet3
rdnetdensity, normalize(network) bwselect(common)
rdn_par_dens "parity/cases/04/py_dens_1.csv" 4 density1
rdbwnet y_out, bwselect(mserd)
rdn_par_bw "parity/cases/04/py_bw_1.csv" 4 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/04/py_manip_1.csv" 4 manipulation
