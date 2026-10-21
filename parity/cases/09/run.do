* parity case 9 (grid)
rdn_par_load "parity/cases/09"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) treatedside(left) classvar(cls) penalty(50.4) connectors(c)
rdn_par_coords "parity/cases/09" 9
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc1) h(407.4) covs(z1)
rdn_par_res "parity/cases/09/py_res_1.csv" 9 rdnet1
rdn_par_agg "parity/cases/09/py_agg_1.csv" 9 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc1) h(265.4) covs(z1)
rdn_par_res "parity/cases/09/py_res_2.csv" 9 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc2) bwselect(cerrd) covs(z1)
rdn_par_res "parity/cases/09/py_res_3.csv" 9 rdnet3
rdnetdensity, normalize(network) bwselect(common)
rdn_par_dens "parity/cases/09/py_dens_1.csv" 9 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/09/py_bw_1.csv" 9 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/09/py_manip_1.csv" 9 manipulation
