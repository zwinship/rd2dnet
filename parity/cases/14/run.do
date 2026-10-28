* parity case 14 (ring)
rdn_par_load "parity/cases/14"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) features(lanes) treatedregion(RG) cost(cost2)
rdn_par_coords "parity/cases/14" 14
rdnet y_out, kernel(triangular) p(2) seed(0) vce(hc0) h(284.1)
rdn_par_res "parity/cases/14/py_res_1.csv" 14 rdnet1
rdn_par_agg "parity/cases/14/py_agg_1.csv" 14 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc2) h(400.1) covs(z1)
rdn_par_res "parity/cases/14/py_res_2.csv" 14 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) bwselect(msetwo) vce(cluster cl)
rdn_par_res "parity/cases/14/py_res_3.csv" 14 rdnet3
rdnetdensity, normalize(network) bwselect(each)
rdn_par_dens "parity/cases/14/py_dens_1.csv" 14 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/14/py_bw_1.csv" 14 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/14/py_manip_1.csv" 14 manipulation
rdn_par_bins "parity/cases/14" 14 350.0
