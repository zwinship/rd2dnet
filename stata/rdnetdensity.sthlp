{smcl}
{* *! version 0.2.0}{...}
{title:Title}

{phang}
{bf:rdnetdensity} {hline 2} Density (sorting) test at a network boundary


{title:Syntax}

{p 8 16 2}
{cmd:rdnetdensity} {ifin} [{cmd:aweight=}{it:w}] [{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt norm:alize(network|none)}}position in cumulative street length (default) or in network distance{p_end}
{synopt:{opt h(# [#])}}bandwidth, or (h_control h_treated); selected when omitted{p_end}
{synopt:{opt bws:elect(common|each)}}one bandwidth for both sides (default) or one per side{p_end}
{synopt:{opt p(#)}}order of the CDF fit (density = its slope); default 2{p_end}
{synopt:{opt q(#)}}order of the bias-corrected fit; default p + 1{p_end}
{synopt:{opt eval:index(#)} {opt wind:ow(#)}}restrict to observations reaching the boundary within # (arc length) of evaluation point #{p_end}
{synopt:{opt ker:nel(name)}}kernel; default {cmd:triangular}{p_end}
{synopt:{opt scalereg:ul(#)}}weight on the variance of the bias estimate; default 1{p_end}
{synopt:{opt uniq:ue}}count each location once (for repeated measurements: panels, repeat sales, sensor readings){p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:rdnetdensity} tests for a discontinuity in the density of observations at
the boundary with the local polynomial density estimator of Cattaneo, Jansson
and Ma (as in {cmd:rddensity}), with jackknife standard errors and
bias-corrected inference, applied to the signed network distance.

{pstd}
With {cmd:normalize(network)} (the default) each observation's distance is
re-expressed as cumulative street length from the boundary, so the test asks
whether there are more observations per metre of street on one side. Jumps in
street supply at the boundary (a park frontage, a highway) then do not pass
for sorting. {cmd:normalize(none)} tests the density of the network distance
itself, as {cmd:rddensity} would. Masked stretches of the boundary are always
excluded.

{pstd}
With repeated measurements, every row counts as an observation and repeated
rows at one place look like a density spike; option {opt unique} tests the
density of locations instead.


{title:Stored results}

{synoptset 20 tabbed}{...}
{synopt:{cmd:r(diff)}}density jump f+ - f- (order p){p_end}
{synopt:{cmd:r(diff_bc)}}bias-corrected jump{p_end}
{synopt:{cmd:r(se)}, {cmd:r(z)}, {cmd:r(pvalue)}}jackknife standard error, z statistic, p-value{p_end}
{synopt:{cmd:r(theta)}}log(f+/f-){p_end}
{synopt:{cmd:r(f_left)}, {cmd:r(f_right)}}densities on each side{p_end}
{synopt:{cmd:r(h_left)}, {cmd:r(h_right)}}bandwidths{p_end}
{synopt:{cmd:r(N_left)}, {cmd:r(N_right)}, {cmd:r(N_h_left)}, {cmd:r(N_h_right)}}sample sizes (all, within the bandwidth){p_end}
{synopt:{cmd:r(normalize)}}normalisation used{p_end}


{title:Example}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)}{p_end}
{phang2}{cmd:. rdnetdensity}{p_end}
{phang2}{cmd:. rdnetdensity, normalize(none)}{p_end}


{title:Also see}

{psee}
{helpb rd2dnet}, {helpb rdnetset}
