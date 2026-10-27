{smcl}
{* *! version 0.2.0}{...}
{viewerjumpto "Syntax" "rdnet##syntax"}{...}
{viewerjumpto "Description" "rdnet##description"}{...}
{viewerjumpto "Options" "rdnet##options"}{...}
{viewerjumpto "Postestimation" "rdnet##estat"}{...}
{viewerjumpto "Stored results" "rdnet##results"}{...}
{viewerjumpto "Examples" "rdnet##examples"}{...}
{title:Title}

{phang}
{bf:rdnet} {hline 2} Boundary discontinuity estimation on a street network


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:rdnet} {depvar} {ifin} [{cmd:aweight=}{it:intensity}] [{cmd:,} {it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Bandwidth}
{synopt:{opt h(numlist)}}one bandwidth, (h_control h_treated) (two numbers always mean this), or one per evaluation point{p_end}
{synopt:{opt hmat:rix(matname)}}J x 2 matrix of (control, treated) bandwidths{p_end}
{synopt:{opt bws:elect(method)}}selector when {opt h()} is omitted: {cmd:mserd} (default), {cmd:msetwo}, {cmd:imserd}, {cmd:imsetwo}, {cmd:cerrd}, {cmd:certwo}, {cmd:icerrd}, {cmd:icertwo}{p_end}
{synopt:{opt scalereg:ul(#)}}bandwidth regularisation weight; default 3{p_end}
{synopt:{opt scalebias:crct(#)}}weight on the higher-order bias term; default 1{p_end}

{syntab:Model}
{synopt:{opt p(#)}}polynomial order for estimation; default 1{p_end}
{synopt:{opt q(#)}}order for bias-corrected inference; default p + 1{p_end}
{synopt:{opt der:iv(a c)}}jump in d^(a+c) mu / ds^a do^c (e.g. {cmd:deriv(1 0)}: kink across the boundary){p_end}
{synopt:{opt tang:vec(ts to)}}jump in the directional derivative ts d/ds + to d/do{p_end}
{synopt:{opt cov:s(varlist)}}covariates with coefficients common to both sides{p_end}
{synopt:{opt fuz:zy(varname)}}treatment take-up: fuzzy design{p_end}
{synopt:{opt pool:ed}}common-slope model y ~ 1 + T + poly(s, o){p_end}
{synopt:{opt ker:nel(name)}}{cmd:triangular} (default), {cmd:epanechnikov}, {cmd:uniform}, {cmd:gaussian}{p_end}
{synopt:{opt min:obs(#)}}minimum local observations per side; default 20 + the number of order-q terms{p_end}

{syntab:Standard errors and inference}
{synopt:{opt vce(vcetype)}}{cmd:hc0}, {cmd:hc1}, {cmd:hc2}, {cmd:hc3} (default), {cmd:cluster} {it:clustvar}, or {cmd:conley} {it:#} [{cmd:network}|{cmd:euclidean}]{p_end}
{synopt:{opt cr1}}with {cmd:vce(cluster)}: rd2d's CR1 cluster variance with normal critical values instead of the default CR2 with t critical values{p_end}
{synopt:{opt conleyresids:cale(#)}}auxiliary fit for Conley residuals at # x h; default 3{p_end}
{synopt:{opt conleyresido:rder(#)}}order of that auxiliary fit; default q{p_end}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:c(level)}{p_end}
{synopt:{opt nocb:ands}}do not add the uniform band to {cmd:e(results)}{p_end}
{synopt:{opt repp(#)}}draws for the uniform critical value; default 2000{p_end}
{synopt:{opt seed(#)}}seed for those draws; default 0, so bands are reproducible; a negative value draws from the session's stream (its RNG state is restored){p_end}
{synoptline}
{p 4 6 2}{cmd:rdnetset} must be run first. {cmd:aweight}s are intensity weights (e.g. floor area) multiplying the kernel weights.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:rdnet} estimates the boundary average treatment effect curve (BATEC) at
each evaluation point along the boundary, carrying rd2d's location-based
estimator over to network coordinates: at evaluation point {it:b} it fits,
separately on each side, a local polynomial of order {opt p()} in
({it:s}, {it:o(b)}) with kernel weights K({it:d(b)}/h). The estimate is the
difference of the two intercepts (the ratio of the outcome and take-up jumps
in fuzzy designs). Robust bias-corrected inference uses the order-{opt q()}
fit. A uniform confidence band across evaluation points uses the joint
covariance of the estimates.

{pstd}
Conley standard errors are a Bartlett spatial HAC over network (or Euclidean)
distance with residuals from a wider auxiliary fit ({opt conleyresidscale()});
residuals from the estimation fit itself absorb spatially correlated shocks
and make the HAC variance too small.

{pstd}
Cluster-robust errors are by default the bias-reduced CR2 variance of Bell and
McCaffrey (2002) with their Satterthwaite degrees of freedom (column {cmd:df} of
{cmd:e(results)}); intervals, p-values, the uniform band and {cmd:estat wbate}
use t critical values. A local fit puts most of its weight on a few clusters
near the point, whose residuals absorb their own cluster effect: in Monte Carlo
runs rd2d's CR1 (option {opt cr1}, factor (n - 1)/(n - k) x g/(g - 1), normal
critical values) covered 0.66-0.91 with 6 to 100 clusters per side, CR2
0.955-0.97.


{marker estat}{...}
{title:Postestimation}

{p 8 16 2}
{cmd:estat wbate}{cmd:,} {opt w:eights}({it:matname} | {it:numlist} | {cmd:equal})

{p 8 16 2}
{cmd:estat lbate} [{cmd:,} {cmd:min} | {cmd:absmax}]

{pstd}
{cmd:estat wbate} reports the weighted boundary average treatment effect
w'tau (rd2d's WBATE). Weights have one entry per evaluation point and are
normalised to sum to one; points that were not estimated must get weight 0
(points with a weak first stage may be given weight).
{cmd:weights(equal)} gives every distinct place of the estimated points equal
weight (points that snapped to the same place, {cmd:node_arc}, share it).
{cmd:weights(density)} weights the estimated points by their local sample n0 + n1
(after {helpb rdnetbins}: the density-weighted BATE), and {cmd:weights(gls)} uses the
minimum-variance weights Sigma^-1 1 / (1' Sigma^-1 1) (a note says when some are
negative). The standard error sqrt(w' Sigma w) uses the joint covariance {cmd:e(cov)}.

{pstd}
{cmd:estat lbate} reports the largest effect max_j tau(b_j) (rd2d's LBATE)
with the interval [max_j(tau_bc - c se), max_j(tau_bc + c se)], c the uniform
critical value. With {cmd:min}, the smallest effect, [min_j lo_j, min_j hi_j];
with {cmd:absmax}, the largest absolute effect,
[max_j dist(0, [lo_j, hi_j]), max_j max(|lo_j|, |hi_j|)].

{pstd}
Because {cmd:e(b)} and {cmd:e(V)} hold the bias-corrected estimates and their
joint covariance, {cmd:test} and {cmd:lincom} work directly, e.g.
{cmd:lincom eval5 - eval10}.


{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}observations{p_end}
{synopt:{cmd:e(N_eval)}}evaluation points{p_end}
{synopt:{cmd:e(N_ok)}}estimated evaluation points{p_end}
{synopt:{cmd:e(crit)}}uniform critical value{p_end}
{synopt:{cmd:e(level)}, {cmd:e(p)}, {cmd:e(q)}, {cmd:e(min_obs)}, {cmd:e(cutoff)}}settings{p_end}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:rdnet}{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(clustvce)} (cr2 or cr1), {cmd:e(kernel)}, {cmd:e(bwselect)}, {cmd:e(cost)}, {cmd:e(fuzzy)}, {cmd:e(covs)}}settings{p_end}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}tau_bc and its covariance at the estimated points{p_end}
{synopt:{cmd:e(results)}}one row per evaluation point: eval arc node_arc h0 h1 n0 n1 tau se tau_bc se_rb ci_lo ci_hi pvalue status h cb_lo cb_hi itt fs fs_se basis df (CR2 degrees of freedom){p_end}
{synopt:{cmd:e(cov)}}J x J covariance of tau_bc (missing for points not estimated){p_end}
{synopt:{cmd:e(bandwidths)}}selected bandwidths (see {helpb rdbwnet}), when {opt h()} is omitted{p_end}

{pstd}
Status codes in {cmd:e(results)}: 0 ok, 1 masked, 2 insufficient
observations, 3 singular local design, 4 no bandwidth, 5 weak first stage,
6 few clusters (fewer than two clusters on a side of the local sample). In fuzzy
designs a point whose first stage is within two standard errors of zero
(|fs| < 2 fs_se) has status 5: its estimates are shown (marked w) and kept in
{cmd:e(results)}, but not in {cmd:e(b)}; the delta-method
interval for itt/fs is unreliable, and {cmd:weights(equal)} and {cmd:estat lbate}
leave it out.

{pstd}
The local fit is a polynomial in (s, o). Where the network meets the boundary
at only a few points (rivers, rail lines, highways crossing a border), o hardly
varies near an evaluation point and that fit is singular; the point is then
fitted in s alone, and the {cmd:basis} column of {cmd:e(results)} is 2 instead
of 1. Put one evaluation point at each crossing ({opt evalpoints()} in
{helpb rdnetset}).


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)}{p_end}
{phang2}{cmd:. rdnet outcome, seed(1)}{p_end}
{phang2}{cmd:. estat wbate, weights(equal)}{p_end}
{phang2}{cmd:. estat lbate}{p_end}
{phang2}{cmd:. rdnet outcome [aw=far], h(600) vce(conley 150)}{p_end}
{phang2}{cmd:. rdnet outcome, h(800) p(2) deriv(1 0)}{p_end}
{phang2}{cmd:. generate lnfar = ln(far)}{p_end}
{phang2}{cmd:. rdnet outcome, h(600) covs(lnfar) vce(cluster block)}{p_end}

{pstd}Donut: drop observations within 40 m of the boundary{p_end}
{phang2}{cmd:. rdnet outcome if abs(_rdn_s) >= 40, h(700)}{p_end}


{title:Also see}

{psee}
{helpb rd2dnet}, {helpb rdnetset}, {helpb rdbwnet}, {helpb rdnetdensity}
