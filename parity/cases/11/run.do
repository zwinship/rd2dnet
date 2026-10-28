* parity case 11 (river)
rdn_par_load "parity/cases/11"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) features(lanes) treatedside(right)
rdn_par_coords "parity/cases/11" 11
rdnet y_out, kernel(triangular) p(1) seed(0) h(404.2) vce(conley 190.6)
rdn_par_res "parity/cases/11/py_res_1.csv" 11 rdnet1
rdn_par_agg "parity/cases/11/py_agg_1.csv" 11 aggregates
rdnet y_out, kernel(triangular) p(1) seed(0) h(432.7) vce(cluster cl)
rdn_par_res "parity/cases/11/py_res_2.csv" 11 rdnet2
rdnet y_out, kernel(uniform) p(1) seed(0) h(327.5 254.5) vce(conley 91.4)
rdn_par_res "parity/cases/11/py_res_3.csv" 11 rdnet3
rdnetdensity, normalize(euclidean) bwselect(each)
rdn_par_dens "parity/cases/11/py_dens_1.csv" 11 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/11/py_bw_1.csv" 11 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/11/py_manip_1.csv" 11 manipulation
rdn_par_bins "parity/cases/11" 11 350.0
