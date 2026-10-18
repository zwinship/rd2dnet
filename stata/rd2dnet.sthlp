{smcl}
{* *! version 0.2.0}{...}
{title:Title}

{phang}
{bf:rd2dnet} {hline 2} Boundary discontinuity designs on street networks


{title:Description}

{pstd}
{cmd:rd2dnet} adapts the location-based boundary discontinuity framework of
{bf:rd2d} (Cattaneo, Titiunik and Yu) from Euclidean coordinates to
network-native coordinates and user-defined network metrics. It is the Stata
version of the Python and R packages of the same name
({browse "https://github.com/zwinship/rd2dnet"}).

{pstd}
{bf:Built on rd2d.} The estimands, side-specific local polynomial estimation,
robust bias-corrected inference, uniform bands and the bandwidth selector (a
port of {cmd:rdbw2d}'s direct plug-in) follow the work of Matias D. Cattaneo,
Rocio Titiunik and Ruiqi Rae Yu. The density test follows Cattaneo, Jansson
and Ma's {cmd:rddensity}. Please cite the rd2d papers listed below when you use
{cmd:rd2dnet}. {cmd:rd2dnet} is an independent project, not affiliated with or
endorsed by the rd2d authors.

{pstd}
Commands:

{p2colset 5 22 24 2}{...}
{p2col:{helpb rdnetset}}declare the network and boundary; compute network coordinates{p_end}
{p2col:{helpb rdnet}}estimate the boundary average treatment effect curve{p_end}
{p2col:{helpb rdnet##estat:estat wbate}}weighted boundary average treatment effect{p_end}
{p2col:{helpb rdnet##estat:estat lbate}}largest boundary average treatment effect{p_end}
{p2col:{helpb rdbwnet}}data-driven bandwidths{p_end}
{p2col:{helpb rdnetdensity}}density (sorting) test at the boundary{p_end}
{p2col:{helpb rdnetmock}}load the bundled mock city{p_end}
{p2colreset}{...}


{title:Method}

{pstd}
For observation {it:i} and an evaluation point {it:b} on the boundary,
{cmd:rd2dnet} replaces rd2d's coordinates (X1 - b1, X2 - b2) with

{p2colset 9 18 20 2}{...}
{p2col:{it:s_i}}the signed network distance to the boundary (negative on the control side){p_end}
{p2col:{it:o_i(b)}}the signed offset along the boundary from {it:b} to where {it:i} reaches it{p_end}
{p2colreset}{...}

{pstd}
and the kernel distance with {it:d_i(b)}, the shortest-path distance from
{it:b} to {it:i}. Distance is whatever the edge cost variable measures (metres,
seconds, or any non-negative dissimilarity), optionally with penalties for
switching road class. The treatment side is always geographic.


{title:Example}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)}{p_end}
{phang2}{cmd:. rdnet outcome [aw=far], vce(cluster block)}{p_end}
{phang2}{cmd:. estat wbate, weights(equal)}{p_end}
{phang2}{cmd:. rdnetdensity}{p_end}


{title:Requirements}

{pstd}
Stata 16 or newer (frames). The computations are written in Mata; the library
{cmd:lrd2dnet.mata} is compiled automatically the first time a command runs in
a session.


{title:References}

{phang}
Cattaneo, M. D., R. Titiunik and R. R. Yu. 2026. rd2d: Causal Inference in
Boundary Discontinuity Designs. Working paper, arXiv:2505.07989.

{phang}
Cattaneo, M. D., R. Titiunik and R. R. Yu. 2026. Estimation and Inference in
Boundary Discontinuity Designs: Distance-Based Methods. {it:Journal of Econometrics} 256(A): 106266.

{phang}
Cattaneo, M. D., R. Titiunik and R. R. Yu. 2026. Estimation and Inference in
Boundary Discontinuity Designs: Location-Based Methods. {it:Journal of the American Statistical Association}, revise and resubmit.

{phang}
Cattaneo, M. D., M. Jansson and X. Ma. 2020. Simple Local Polynomial Density
Estimators. {it:Journal of the American Statistical Association} 115(531): 1449-1455.

{phang}
Calonico, S., M. D. Cattaneo and R. Titiunik. 2014. Robust Nonparametric
Confidence Intervals for Regression-Discontinuity Designs. {it:Econometrica} 82(6): 2295-2326.

{phang}
Conley, T. G. 1999. GMM Estimation with Cross Sectional Dependence.
{it:Journal of Econometrics} 92(1): 1-45.


{title:Author}

{pstd}
Zach Winship, zachtwinship@gmail.com{break}
{browse "https://github.com/zwinship/rd2dnet"}


{title:License}

{pstd}
GPL-3, the same license as rd2d and rddensity, from which parts of the
implementation are ported.
