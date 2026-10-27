{smcl}
{* *! version 0.2.0}{...}
{title:Title}

{phang}
{bf:rdnetbins} {hline 2} Nearest-crossing bins: an additional design for network boundary discontinuities


{title:Syntax}

{p 8 16 2}
{cmd:rdnetbins} [{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt notop:ology}}do not use the topology features (only those in {cmd:rdnetset}'s {opt features()}){p_end}
{synopt:{opt nosim:ilarity}}no edge similarity: plain bins (no kernel in S, no R*S slopes){p_end}
{synopt:{opt nobr:anches}}no branch-specific slopes{p_end}
{synopt:{opt minbr:anch(#)}}a branch gets its own slope only with at least # distinct locations on its side of the local sample; default 10 (smaller branches share the reference slope: the trunk's, or else the largest branch's){p_end}
{synopt:{opt hs(rot|#|inf)}}similarity bandwidth h_S; default {cmd:rot}; {cmd:inf} keeps the R*S slopes without the kernel{p_end}
{synopt:{opt kernels(string)}}similarity kernel K2: {cmd:triangular} (default), {cmd:epanechnikov}, {cmd:uniform} or {cmd:gaussian}{p_end}
{synoptline}


{title:Description}

{pstd}
The default design of {cmd:rdnet} estimates the effect at points along the
boundary from every observation near each point, with a polynomial in the
signed distance s and the offset o along the boundary. {cmd:rdnetbins} switches
the coordinates in memory to an additional design, which compares each
observation only with others that reach the boundary at the same crossing:

{phang2}
1. {it:Bins.} Every observation belongs to the boundary crossing nearest to it
along the network (variable {cmd:_rdn_bin}); R = s is its signed distance to that
crossing. Each bin becomes one evaluation point, and its local sample is the
bin's members only. Where a street runs along the boundary, each of its nodes is
a crossing.

{phang2}
2. {it:Edge similarity.} Each street segment has a feature vector: the edge
variables named in {cmd:rdnetset}'s {opt features()} plus its topology on the
network as given ({it:n_adjacent}: segments sharing an end; {it:n_second}:
segments two steps away; {it:degree_max}: the larger end degree;
{it:clustering}: the mean clustering coefficient of its ends). S
({cmd:_rdn_S}) is the Mahalanobis distance between an observation's segment and
the crossing's segment, with the covariance of the features over all segments.

{phang2}
3. {it:Bin fits.} After {cmd:rdnetbins}, {cmd:rdnet} fits each bin with the
product kernel K1(R/h_R) K2(S/h_S). The jump is common to the bin. The slopes in
R differ by side and by branch ({cmd:_rdn_branch}: the edge where the path from
the crossing first passes a junction; 0 before it). R*S enters with a slope per
side. Every {cmd:rdnet} option works: bias correction, the h_R selector (per bin,
with the K2 weights, without the branch and S slopes), cluster and Conley errors,
fuzzy designs and covariates.

{phang2}
4. {it:Aggregation.} {cmd:estat wbate, weights(density)} weights bins by their
local sample (the density-weighted BATE); {cmd:estat wbate, weights(gls)} uses the
minimum-variance weights cov^-1 1 / (1' cov^-1 1). Both use se = sqrt(w' cov w);
with {cmd:vce(conley ...)} the covariance includes the correlation between
neighbouring observations in different bins.

{pstd}
{cmd:rdnetfe} pools all bins in one regression with crossing fixed effects, for
boundaries with many sparsely populated crossings.

{pstd}
{it:Choosing h_S.} {cmd:hs(rot)} is sqrt(2 chi2_k(0.95)), k the number of
independent features: if two segments' features were independent draws, it
would keep all but the 5% most dissimilar pairs, with weights falling linearly
in S. It is scale-free but not optimal in any sense, so report the BATE over a
few values (see the example). A data-driven selector for h_S is left for future
research: S is not a running variable extrapolated to a cutoff, so the
MSE-optimal plug-in theory behind {cmd:rdbwnet} does not cover it.

{pstd}
Run {cmd:rdnetbins} again to change its options (h_S, for example). Run
{cmd:rdnetset} again to return to the default design.


{title:Examples}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) features(speed_kph)}{p_end}
{phang2}{cmd:. rdnetbins}{p_end}
{phang2}{cmd:. rdnet outcome, h(700) minobs(10)}{p_end}
{phang2}{cmd:. estat wbate, weights(density)}{p_end}
{phang2}{cmd:. estat wbate, weights(gls)}{p_end}
{phang2}{cmd:. rdnetfe outcome, h(700)}{p_end}

{pstd}
The BATE over a range of h_S (half, one and two times the rule of thumb, and no kernel):

{phang2}{cmd:. rdnetbins}{p_end}
{phang2}{cmd:. local h0 = r(h_S)}{p_end}
{phang2}{cmd:. foreach f in 0.5 1 2 {c -(}}{p_end}
{phang2}{cmd:.     rdnetbins, hs(`=`f' * `h0'')}{p_end}
{phang2}{cmd:.     quietly rdnet outcome, h(700) minobs(10)}{p_end}
{phang2}{cmd:.     estat wbate, weights(density)}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. rdnetbins, hs(inf)}{p_end}


{title:Stored results}

{pstd}
{cmd:rdnetbins} stores in {cmd:r()}:

{synoptset 20 tabbed}{...}
{synopt:{cmd:r(N_bins)}}number of bins{p_end}
{synopt:{cmd:r(h_S)}}similarity bandwidth (missing: no kernel){p_end}
{synopt:{cmd:r(rank)}}rank of the feature covariance{p_end}
{synopt:{cmd:r(bins)}}one row per bin: bin, arc, n0, n1, n_branches, masked, then the crossing's features{p_end}

{pstd}
It creates {cmd:_rdn_bin}, {cmd:_rdn_branch} and {cmd:_rdn_S}.


{title:Also see}

{psee}
{help rdnetfe}, {help rdnet}, {help rdnetset}, {help rd2dnet}
{p_end}
