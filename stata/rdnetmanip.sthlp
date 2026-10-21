{smcl}
{* *! version 0.2.0}{...}
{title:Title}

{phang}
{bf:rdnetmanip} {hline 2} Manipulation (sorting) test at the boundary (Bugni, Crippa and Restrepo)


{title:Syntax}

{p 8 16 2}
{cmd:rdnetmanip} {ifin} [{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt q(numlist)}}locations nearest the boundary; several values give a path; default 200 (a tenth of the usable locations when there are fewer than 2,000){p_end}
{synopt:{opt km:ax(#)}}largest number of clusters along the boundary; default q{p_end}
{synopt:{opt rep:s(#)}}relabelling draws for the null distribution; default 10000{p_end}
{synopt:{opt l:evel(#)}}confidence level for the reported critical value{p_end}
{synopt:{opt seed(#)}}seed of the relabelling draws; default 0{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:rdnetmanip} implements the manipulation test of Bugni, Crippa and Restrepo
(2026). Without sorting across the boundary, the locations closest to it are
about equally likely to lie on either side, anywhere along the boundary. The
q locations nearest the boundary (straight-line distance in the plane: sorting
is about where units locate, not about the network) are projected onto it and
grouped into k contiguous clusters for every k from 1 to {opt kmax()}; the
statistic is the largest standardised imbalance max sqrt(n) |share treated -
1/2|, and its null distribution comes from relabelling the q locations as
independent Bernoulli(1/2) draws with the clusters held fixed.

{pstd}
The clusters are the exact optimal 1-D k-means partitions of the positions along
the boundary (a closed ring is cut at its largest gap), and the draws come from
the Park-Miller generator, so the Stata, Python and R versions give the same
answer. Locations count once; locations on the boundary line, in masked
stretches, or projecting onto an end of an open boundary are left out.

{pstd}
The test has power against sorting that offsets along the boundary (one side
denser on one stretch, the other on another), which a density test of the
distance to the boundary cannot detect; against sorting in one direction
everywhere, {helpb rdnetdensity} is more powerful. Report several values of
q, e.g. {cmd:q(50 100 200)}. The test is not valid for observations that lie on
or next to the network itself (stations along a river, points on streets, GPS
fixes); a warning is shown then.


{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(statistic)}}S (for the first q){p_end}
{synopt:{cmd:r(pvalue)}}randomisation p-value{p_end}
{synopt:{cmd:r(crit)}}critical value at {opt level()}{p_end}
{synopt:{cmd:r(q)}}number of locations used{p_end}
{synopt:{cmd:r(n_treated)}}of which on the treated side{p_end}
{synopt:{cmd:r(distance)}}distance to the boundary of the farthest location used{p_end}
{synopt:{cmd:r(k_max)}}largest number of clusters{p_end}
{synopt:{cmd:r(reps)}}draws{p_end}
{synopt:{cmd:r(N_locations)}}usable locations{p_end}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(path)}}one row per q: q, statistic, pvalue, crit, n_treated, distance, k_max{p_end}


{title:Examples}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)}{p_end}
{phang2}{cmd:. rdnetmanip, q(50 100 200) reps(2000)}{p_end}


{title:Reference}

{phang}
Bugni, F. A., F. Crippa and D. Restrepo. 2026. Manipulation testing in boundary
discontinuity designs. arXiv:2609.00350.


{title:Also see}

{psee}
{helpb rd2dnet}, {helpb rdnetdensity}, {helpb rdnetset}
