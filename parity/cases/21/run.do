* parity case 21 (grid)
rdn_par_load "parity/cases/21"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) features(lanes) treatedregion(RG) directed direction(flow) cost(cost2)
rdn_par_coords "parity/cases/21" 21
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc3) bwselect(certwo) covs(z1)
rdn_par_res "parity/cases/21/py_res_1.csv" 21 rdnet1
rdn_par_agg "parity/cases/21/py_agg_1.csv" 21 aggregates
rdnet y_out [aw=w], kernel(epanechnikov) p(1) seed(0) vce(hc1) h(257.4)
rdn_par_res "parity/cases/21/py_res_2.csv" 21 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc0) bwselect(imsetwo) covs(z1)
rdn_par_res "parity/cases/21/py_res_3.csv" 21 rdnet3
rdnetdensity, normalize(network) bwselect(each)
rdn_par_dens "parity/cases/21/py_dens_1.csv" 21 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/21/py_bw_1.csv" 21 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/21/py_manip_1.csv" 21 manipulation
rdn_par_bins "parity/cases/21" 21 350.0
