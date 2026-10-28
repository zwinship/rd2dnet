* parity case 12 (grid)
rdn_par_load "parity/cases/12"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) features(lanes) treatedregion(RG) classvar(cls) penalty(156.9) connectors(c)
rdn_par_coords "parity/cases/12" 12
rdnet y_out, kernel(uniform) p(2) seed(0) vce(hc2) h(317.8 252.6)
rdn_par_res "parity/cases/12/py_res_1.csv" 12 rdnet1
rdn_par_agg "parity/cases/12/py_agg_1.csv" 12 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) h(427.0 318.6) deriv(1 0)
rdn_par_res "parity/cases/12/py_res_2.csv" 12 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) bwselect(msetwo) vce(conley 140.6)
rdn_par_res "parity/cases/12/py_res_3.csv" 12 rdnet3
rdnetdensity, normalize(network) bwselect(common)
rdn_par_dens "parity/cases/12/py_dens_1.csv" 12 density1
rdbwnet y_out, bwselect(msetwo)
rdn_par_bw "parity/cases/12/py_bw_1.csv" 12 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/12/py_manip_1.csv" 12 manipulation
rdn_par_bins "parity/cases/12" 12 350.0
