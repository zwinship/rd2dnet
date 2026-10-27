{smcl}
{* *! version 0.2.0}{...}
{title:Title}

{phang}
{bf:rdnetfe} {hline 2} Pooled nearest-crossing BATE with crossing fixed effects


{title:Syntax}

{p 8 16 2}
{cmd:rdnetfe} {it:depvar} {ifin} [{it:aweight}] [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt h(numlist)}}bandwidth: one number, (h_control h_treated), or one per bin; default: selected per bin{p_end}
{synopt:{opt hmat:rix(name)}}bins x 2 matrix of (control, treated) bandwidths{p_end}
{synopt:{opt p(#)}}polynomial order in R; default 1{p_end}
{synopt:{opt q(#)}}order for the bias-corrected estimate; default p + 1{p_end}
{synopt:{opt ker:nel(string)}}kernel K1; default triangular{p_end}
{synopt:{opt bws:elect(string)}}selector when h() is not given; default mserd{p_end}
{synopt:{opt vce(string)}}{cmd:hc0}-{cmd:hc3} (default {cmd:hc3}), {cmd:cluster} {it:clustvar} (CR2; add {opt cr1} for CR1), or {cmd:conley} {it:#} [{cmd:network}|{cmd:euclidean}]{p_end}
{synopt:{opt cov:s(varlist)}}covariates{p_end}
{synopt:{opt noipw}}plain fixed-effects weights instead of propensity weights{p_end}
{synopt:{opt nobins:lopes}}do not interact the slopes with the crossings' features{p_end}
{synopt:{opt min:obs(#)}}distinct locations per side for a bin to enter; default 1{p_end}
{synopt:{opt l:evel(#)}}confidence level{p_end}
{synoptline}


{title:Description}

{pstd}
Run {cmd:rdnetset} and {cmd:rdnetbins} first. {cmd:rdnetfe} is for boundaries
that cross the network at many points with few observations at each. It fits
one weighted regression over the members of all bins:

{p 8 8 2}
y = a_m + T_B T + poly_p(R) + T poly_p(R) + R Z_m'th1 + T R Z_m'th2 + d1 R S + d2 T R S + covs

{pstd}
Here T is the 0/1 treatment indicator and R = s/h, with each side's bandwidth.
a_m is one intercept per bin, so the jump is identified only from comparisons
within bins. Z_m are the bin's crossing features, standardised across bins, and
S is the edge dissimilarity. The weights are the product kernel times, unless
{opt noipw} is given, 1/p_m for treated and 1/(1 - p_m) for control members.
p_m is the kernel-weighted treated share of bin m.

{pstd}
Plain fixed-effects regression weights a bin by p_m(1 - p_m) times its mass.
After the reweighting, both sides of every bin carry the bin's kernel mass, so
T_B weights bins by their local sample: it estimates the density-weighted BATE.
{cmd:e(tau)} is the order-p estimate. {cmd:e(tau_bc)}, {cmd:e(se_rb)}, the
interval and the p-value come from the order-q fit, as in {cmd:rdnet}. With
{cmd:vce(conley ...)} the errors use the pooled fit's own residuals.


{title:Examples}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)}{p_end}
{phang2}{cmd:. rdnetbins}{p_end}
{phang2}{cmd:. rdnetfe outcome, h(700)}{p_end}
{phang2}{cmd:. rdnetfe outcome, h(700) vce(cluster block)}{p_end}


{title:Stored results}

{pstd}
{cmd:rdnetfe} stores in {cmd:e()}:

{synoptset 20 tabbed}{...}
{synopt:{cmd:e(b), e(V)}}tau_bc and se_rb^2 (column BATE){p_end}
{synopt:{cmd:e(tau), e(se)}}order-p estimate and its standard error{p_end}
{synopt:{cmd:e(tau_bc), e(se_rb)}}bias-corrected estimate and robust standard error{p_end}
{synopt:{cmd:e(ci_lo), e(ci_hi), e(pvalue)}}robust interval and p-value{p_end}
{synopt:{cmd:e(df)}}CR2 degrees of freedom (cluster errors){p_end}
{synopt:{cmd:e(N_bins), e(N_bins_used), e(N_fit)}}bins, bins with both sides, observations in the fit{p_end}
{synopt:{cmd:e(bins)}}per bin: bin, arc, h0, h1, n0, n1, phat, used{p_end}
{synopt:{cmd:e(terms)}}slope terms in the regression{p_end}


{title:Also see}

{psee}
{help rdnetbins}, {help rdnet}, {help rd2dnet}
{p_end}
