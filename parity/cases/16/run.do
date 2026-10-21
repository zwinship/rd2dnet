* parity case 16 (grid)
rdn_par_load "parity/cases/16"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) treatedside(left)
rdn_par_coords "parity/cases/16" 16
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc2) h(258.4 429.8) fuzzy(D)
rdn_par_res "parity/cases/16/py_res_1.csv" 16 rdnet1
rdn_par_agg "parity/cases/16/py_agg_1.csv" 16 aggregates
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc2) h(389.2) deriv(1 0)
rdn_par_res "parity/cases/16/py_res_2.csv" 16 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) vce(hc2) bwselect(cerrd) covs(z1)
rdn_par_res "parity/cases/16/py_res_3.csv" 16 rdnet3
rdnetdensity, normalize(euclidean) bwselect(common)
rdn_par_dens "parity/cases/16/py_dens_1.csv" 16 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/16/py_bw_1.csv" 16 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/16/py_manip_1.csv" 16 manipulation
