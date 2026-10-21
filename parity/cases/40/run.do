* parity case 40 (grid)
rdn_par_load "parity/cases/40"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) treatedside(left) directed direction(flow) classvar(cls) penalty(220.1) connectors(c)
rdn_par_coords "parity/cases/40" 40
rdnet y_out, kernel(triangular) p(1) seed(0) bwselect(mserd) vce(conley 194.0)
rdn_par_res "parity/cases/40/py_res_1.csv" 40 rdnet1
rdn_par_agg "parity/cases/40/py_agg_1.csv" 40 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) bwselect(imserd) vce(conley 172.4)
rdn_par_res "parity/cases/40/py_res_2.csv" 40 rdnet2
rdnet y_out, kernel(triangular) p(2) seed(0) vce(hc0) bwselect(imserd)
rdn_par_res "parity/cases/40/py_res_3.csv" 40 rdnet3
rdnetdensity, normalize(network) bwselect(common)
rdn_par_dens "parity/cases/40/py_dens_1.csv" 40 density1
rdbwnet y_out, bwselect(mserd)
rdn_par_bw "parity/cases/40/py_bw_1.csv" 40 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/40/py_manip_1.csv" 40 manipulation
