* parity case 1 (grid)
rdn_par_load "parity/cases/01"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) treatedside(left) classvar(cls) penalty(86.0) connectors(c)
rdn_par_coords "parity/cases/01" 1
rdnet y_out, kernel(uniform) p(1) seed(0) h(415.5 331.8) vce(conley 146.0)
rdn_par_res "parity/cases/01/py_res_1.csv" 1 rdnet1
rdn_par_agg "parity/cases/01/py_agg_1.csv" 1 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc0) bwselect(certwo) covs(z1)
rdn_par_res "parity/cases/01/py_res_2.csv" 1 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) bwselect(mserd) vce(conley 134.4)
rdn_par_res "parity/cases/01/py_res_3.csv" 1 rdnet3
rdnetdensity, normalize(none) bwselect(common)
rdn_par_dens "parity/cases/01/py_dens_1.csv" 1 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/01/py_bw_1.csv" 1 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/01/py_manip_1.csv" 1 manipulation
