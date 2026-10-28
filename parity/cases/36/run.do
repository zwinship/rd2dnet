* parity case 36 (grid)
rdn_par_load "parity/cases/36"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) features(lanes) treatedregion(RG)
rdn_par_coords "parity/cases/36" 36
rdnet y_out, kernel(uniform) p(1) seed(0) bwselect(cerrd) vce(conley 106.7)
rdn_par_res "parity/cases/36/py_res_1.csv" 36 rdnet1
rdn_par_agg "parity/cases/36/py_agg_1.csv" 36 aggregates
rdnet y_out [aw=w], kernel(triangular) p(1) seed(0) vce(hc3) h(435.3)
rdn_par_res "parity/cases/36/py_res_2.csv" 36 rdnet2
rdnet y_out, kernel(triangular) p(2) seed(0) vce(hc3) bwselect(imsetwo)
rdn_par_res "parity/cases/36/py_res_3.csv" 36 rdnet3
rdnetdensity, normalize(none) bwselect(each)
rdn_par_dens "parity/cases/36/py_dens_1.csv" 36 density1
rdbwnet y_out, bwselect(imserd)
rdn_par_bw "parity/cases/36/py_bw_1.csv" 36 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/36/py_manip_1.csv" 36 manipulation
rdn_par_bins "parity/cases/36" 36 350.0
