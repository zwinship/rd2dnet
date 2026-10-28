* parity case 28 (grid)
rdn_par_load "parity/cases/28"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(9) features(lanes) treatedregion(RG)
rdn_par_coords "parity/cases/28" 28
rdnet y_out, kernel(triangular) p(2) seed(0) vce(hc2) bwselect(certwo)
rdn_par_res "parity/cases/28/py_res_1.csv" 28 rdnet1
rdn_par_agg "parity/cases/28/py_agg_1.csv" 28 aggregates
rdnet y_out, kernel(uniform) p(2) seed(0) vce(hc1) bwselect(msetwo)
rdn_par_res "parity/cases/28/py_res_2.csv" 28 rdnet2
rdnet y_out, kernel(triangular) p(1) seed(0) vce(hc1) h(327.2) covs(z1)
rdn_par_res "parity/cases/28/py_res_3.csv" 28 rdnet3
rdnetdensity, normalize(euclidean) bwselect(each)
rdn_par_dens "parity/cases/28/py_dens_1.csv" 28 density1
rdbwnet y_out, bwselect(mserd)
rdn_par_bw "parity/cases/28/py_bw_1.csv" 28 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/28/py_manip_1.csv" 28 manipulation
rdn_par_bins "parity/cases/28" 28 350.0
