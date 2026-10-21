* parity case 33 (grid)
rdn_par_load "parity/cases/33"
rdnetset x y, nodes(nodes) edges(edges) boundary(BL) neval(7) treatedregion(RG) cost(cost2)
rdn_par_coords "parity/cases/33" 33
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc3) bwselect(certwo) deriv(1 0)
rdn_par_res "parity/cases/33/py_res_1.csv" 33 rdnet1
rdn_par_agg "parity/cases/33/py_agg_1.csv" 33 aggregates
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc0) bwselect(mserd)
rdn_par_res "parity/cases/33/py_res_2.csv" 33 rdnet2
rdnet y_out, kernel(epanechnikov) p(1) seed(0) vce(hc2) bwselect(imserd) covs(z1)
rdn_par_res "parity/cases/33/py_res_3.csv" 33 rdnet3
rdnetdensity, normalize(euclidean) bwselect(each)
rdn_par_dens "parity/cases/33/py_dens_1.csv" 33 density1
rdbwnet y_out, bwselect(cerrd)
rdn_par_bw "parity/cases/33/py_bw_1.csv" 33 rdbwnet1
rdnetmanip, q(60) reps(500)
rdn_par_manip "parity/cases/33/py_manip_1.csv" 33 manipulation
