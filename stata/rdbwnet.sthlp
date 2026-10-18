{smcl}
{* *! version 0.1.0}{...}
{title:Title}

{phang}
{bf:rdbwnet} {hline 2} Bandwidth selection for network boundary discontinuity designs


{title:Syntax}

{p 8 16 2}
{cmd:rdbwnet} {depvar} {ifin} [{cmd:aweight=}{it:intensity}] [{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt p(#)}}polynomial order; default 1{p_end}
{synopt:{opt bws:elect(method)}}{cmd:mserd} (default), {cmd:msetwo}, {cmd:imserd}, {cmd:imsetwo}, {cmd:cerrd}, {cmd:certwo}, {cmd:icerrd}, {cmd:icertwo}{p_end}
{synopt:{opt der:iv(a c)}, {opt tang:vec(ts to)}}derivative targets, as in {helpb rdnet}{p_end}
{synopt:{opt cov:s(varlist)}}covariates partialled out with common coefficients{p_end}
{synopt:{opt fuz:zy(varname)}}take-up; uses the linearised outcome{p_end}
{synopt:{opt vce(hc0|hc1|hc2|hc3)}}variance used inside the selector; default {cmd:hc1}{p_end}
{synopt:{opt cl:uster(varname)}}cluster-robust variance inside the selector{p_end}
{synopt:{opt ker:nel(name)}}kernel; default {cmd:triangular}{p_end}
{synopt:{opt bwch:eck(#)}}minimum bandwidth: distance to the #-th nearest observation per side; default 20 + the number of order-(p+1) terms{p_end}
{synopt:{opt scalereg:ul(#)}}regularisation weight; default 1{p_end}
{synopt:{opt scalebias:crct(#)}}higher-order bias weight; default 1{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:rdbwnet} selects MSE- or coverage-error-optimal bandwidths (in the units
of the network distance) at each evaluation point, porting rd2d's direct
plug-in selector ({cmd:rdbw2d}) to network coordinates: a rule-of-thumb pilot,
a separate pilot for the bias, regularisation and a higher-order bias term.
{cmd:mserd} uses one bandwidth for both sides, {cmd:msetwo} one per side, and
the {cmd:i*} variants one bandwidth for all evaluation points. The {cmd:cer*}
variants shrink the MSE-optimal bandwidth by n^(1/(2p+4) - 1/(p+4)).
{helpb rdnet} calls it when {opt h()} is omitted.


{title:Stored results}

{synoptset 20 tabbed}{...}
{synopt:{cmd:r(bws)}}one row per evaluation point: eval arc status dn bn0 bn1 n0 n1 hmin0 hmin1 hmax0 hmax1 V0 B0 Reg10 Reg20 V1 B1 Reg11 Reg21 h0 h1 h{p_end}
{p 4 6 2}Status: 0 ok, 1 masked, 2 insufficient, 3 singular. Suffix 0 is the control side, 1 the treated side.{p_end}


{title:Example}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)}{p_end}
{phang2}{cmd:. rdbwnet outcome, bwselect(msetwo)}{p_end}
{phang2}{cmd:. matrix list r(bws)}{p_end}


{title:Also see}

{psee}
{helpb rd2dnet}, {helpb rdnet}
