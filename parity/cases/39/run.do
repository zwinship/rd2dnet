* parity case 39 (river)
rdn_par_load "parity/cases/39"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(5) features(lanes) treatedside(right)
rdn_par_coords "parity/cases/39" 39
rdnet y_out [aw=w], kernel(triangular) p(1) seed(0) vce(hc2) bwselect(mserd)
rdn_par_res "parity/cases/39/py_res_1.csv" 39 rdnet1
rdn_par_agg "parity/cases/39/py_agg_1.csv" 39 aggregates
rdnet y_out, kernel(epanechnikov) p(2) seed(0) vce(hc1) bwselect(imsetwo)
rdn_par_res "parity/cases/39/py_res_2.csv" 39 rdnet2
rdnet y_out [aw=w], kernel(epanechnikov) p(1) seed(0) vce(hc3) bwselect(imsetwo)
rdn_par_res "parity/cases/39/py_res_3.csv" 39 rdnet3
rdnetdensity, normalize(network) bwselect(each)
rdn_par_dens "parity/cases/39/py_dens_1.csv" 39 density1
rdbwnet y_out, bwselect(mserd)
rdn_par_bw "parity/cases/39/py_bw_1.csv" 39 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/39/py_manip_1.csv" 39 manipulation
rdn_par_bins "parity/cases/39" 39 350.0
