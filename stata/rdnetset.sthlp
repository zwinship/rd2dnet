{smcl}
{* *! version 0.2.0}{...}
{viewerjumpto "Syntax" "rdnetset##syntax"}{...}
{viewerjumpto "Description" "rdnetset##description"}{...}
{viewerjumpto "Options" "rdnetset##options"}{...}
{viewerjumpto "Stored results" "rdnetset##results"}{...}
{viewerjumpto "Examples" "rdnetset##examples"}{...}
{title:Title}

{phang}
{bf:rdnetset} {hline 2} Declare a street network and boundary and compute network coordinates


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:rdnetset} {it:xvar} {it:yvar} {ifin}{cmd:,}
{opt edges(framename)} [{opt nodes(framename)}]
{c -(}{opt boundary(matname)} | {opt evalpoints(matname)}{c )-}
[{opt treatedregion(matname)} | {opt treatedside(left|right)} | {opt treated(varname)}]
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Network}
{synopt:{opt edges(framename)}}frame with one row per street segment (required){p_end}
{synopt:{opt nodes(framename)}}frame with one row per network node; without it, {opt edges()} holds segment end coordinates{p_end}
{synopt:{opt ecoord:s(x1 y1 x2 y2)}}segment end coordinates in {opt edges()} when there is no {opt nodes()}; default {cmd:x1 y1 x2 y2}{p_end}
{synopt:{opt join:tol(#)}}segment ends closer than this become one node; default 1e-6{p_end}

{syntab:Boundary and treated side}
{synopt:{opt boundary(matname)}}k x 2 matrix of boundary vertices (projected coordinates){p_end}
{synopt:{opt evalpoints(matname)}}evaluation points; without {opt boundary()}, a J x 2 matrix of points along the boundary, in order, through which the boundary is drawn (as in rd2d){p_end}
{synopt:{opt treatedregion(matname)}}polygon containing the treated side{p_end}
{synopt:{opt treatedside(left|right)}}treated side relative to the boundary's direction{p_end}
{synopt:{opt treated(varname)}}treatment indicator; also infers the treated side when neither option above is given{p_end}

{syntab:Network variables}
{synopt:{opt nodeid(varname)}}node id in the node frame; default {cmd:id}{p_end}
{synopt:{opt nodex(varname)} {opt nodey(varname)}}node coordinates; default {cmd:x} and {cmd:y}{p_end}
{synopt:{opt from(varname)} {opt to(varname)}}edge endpoints in the edge frame; default {cmd:from} and {cmd:to}{p_end}
{synopt:{opt cost(varname)}}edge cost that defines distance; default {cmd:length}{p_end}
{synopt:{opt length(varname)}}geometric edge length; default {cmd:length}{p_end}
{synopt:{opt dir:ected}}edges are one-way ({it:from} -> {it:to}){p_end}
{synopt:{opt direction(both|out|in)}}how directed edges are traversed; default {cmd:both}{p_end}

{syntab:Road-class transition penalties}
{synopt:{opt classvar(varname)}}edge class (e.g. OSM {cmd:highway}); turns penalties on{p_end}
{synopt:{opt penalty(#)}}penalty for switching between different classes; default 0{p_end}
{synopt:{opt penaltymatrix(matname)}}class-by-class penalties (row and column names are classes){p_end}
{synopt:{opt connectors(string)}}transparent classes (ramps, links), separated by spaces{p_end}

{syntab:Evaluation points and masks}
{synopt:{opt maskarcs(matname)}}m x 2 matrix of masked arc-length intervals{p_end}
{synopt:{opt maskpolygons(matname)}}k x 2 (one polygon) or k x 3 (x, y, polygon id) masking polygons{p_end}
{synopt:{opt neval(#)}}number of equally spaced evaluation points; default 25{p_end}
{synopt:{opt tol(#)}}distance (in coordinate units) at which a node lies on the boundary; default 1e-6. Increase it (e.g. 0.5 metres) when the boundary follows a street{p_end}

{syntab:Snapping}
{synopt:{opt snapscale(#)}}cost per unit of off-network snapping distance; default 1 if {cmd:cost(length)}, else 0{p_end}
{synopt:{opt snapnode}}snap observations to the nearest node instead of the nearest street{p_end}
{synopt:{opt nosplit}}do not split edges where they cross the boundary{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:rdnetset} declares the street network (two frames: nodes and edges) and
the boundary, snaps each observation at ({it:xvar}, {it:yvar}) to the network
on its own side of the boundary, and computes its network coordinates. All
coordinates must be in the same projected (planar) units.

{pstd}
Distance is the sum of {opt cost()} along the cheapest path: use lengths for
metres, travel times for seconds, or any non-negative edge dissimilarity you
build in Stata (for example {cmd:length + 50*abs(landuse change)}). With
{opt classvar()}, switching road class costs {opt penalty()} (or the entry of
{opt penaltymatrix()}); classes listed in {opt connectors()} are transparent,
so highway -> ramp -> highway is cheap while highway -> ramp -> local street is not.

{pstd}
Edges that cross the boundary are split at the crossing so that every path
between the two sides passes through a node on the boundary. Edges are treated
as straight segments between their nodes.

{pstd}
{cmd:rdnetset} creates the variables {cmd:_rdn_id}, {cmd:_rdn_s} (signed
network distance to the boundary), {cmd:_rdn_treated} and {cmd:_rdn_arc}
(where the observation reaches the boundary), and keeps the network
coordinates in memory for {helpb rdnet}, {helpb rdbwnet} and
{helpb rdnetdensity}. They are not saved with the data: rerun {cmd:rdnetset}
after loading the data in a new session. Repeated locations (panels, sensors)
are processed once.


{pstd}
{bf:Data entry as in rd2d.} rd2d takes an outcome, locations, a treatment
indicator and points b along the boundary. Give the same here: {it:xvar}
{it:yvar}, {opt treated()} and {opt evalpoints()} (a J x 2 matrix of the
points b, in order along the boundary). The boundary is the line through
them and the treated side is the side where the treated observations lie.

{pstd}
{bf:Networks from GIS layers.} A road (river, rail) layer exported as one row
per segment with its end coordinates needs no node table: put it in a frame
with variables {cmd:x1 y1 x2 y2} (and any cost variables) and omit
{opt nodes()}. Segment ends within {opt jointol()} of each other become one
node. Lines must be split into straight segments at their vertices and at
junctions.

{pstd}
{bf:Directed networks.} With {opt directed} and {opt direction(out)},
distances run from the boundary to each observation along edge directions (on
a river drawn downstream, observations upstream of the boundary cannot be
reached); {opt direction(in)} measures distances from each observation to the
boundary.

{pstd}
{bf:Closed boundaries.} When the first and last vertices of {opt boundary()}
coincide (a zone outline), offsets along the boundary go the shorter way round
and the default evaluation points go around the whole ring.


{marker options}{...}
{title:Options}

{phang}
{opt treatedregion(matname)} gives the treated zone as a polygon: a k x 2
matrix of vertices, or k x 3 with a ring id in the third column (holes and
several parts use the even-odd rule). Observations inside it, or within
{opt tol()} of it, are treated.

{phang}
{opt maskarcs()} and {opt maskpolygons()} mark stretches of the boundary with
no mass on one side (a park, a river). Evaluation points there are not
estimated, and observations reaching the boundary there are left out of the
density test.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:rdnetset} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}observations{p_end}
{synopt:{cmd:r(N_treated)}}treated observations{p_end}
{synopt:{cmd:r(N_locations)}}distinct locations{p_end}
{synopt:{cmd:r(N_eval)}}evaluation points{p_end}
{synopt:{cmd:r(N_masked)}}masked evaluation points{p_end}
{synopt:{cmd:r(N_unreachable)}}observations with no network path to the boundary{p_end}
{synopt:{cmd:r(N_crossings)}}edges split at the boundary{p_end}
{synopt:{cmd:r(N_nodes)}, {cmd:r(N_edges)}}network size after splitting{p_end}
{synopt:{cmd:r(length)}}length of the boundary line{p_end}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(eval)}}evaluation points: eval, arc, node_arc, x, y, masked{p_end}


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)}{p_end}

{pstd}rd2d-style entry: boundary through the evaluation points, treated side from the treatment variable:{p_end}
{phang2}{cmd:. matrix B = (100, 1600 \ 1000, 1600 \ 2000, 1600 \ 2900, 1600)}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) evalpoints(B) treated(treated)}{p_end}

{pstd}A road layer of segments with end coordinates (no node table):{p_end}
{phang2}{cmd:. frame create roads}{p_end}
{phang2}{cmd:. frame roads: import delimited roads.csv}{p_end}
{phang2}{cmd:. rdnetset x y, edges(roads) ecoords(x1 y1 x2 y2) jointol(0.5) boundary(line) treatedside(left)}{p_end}

{pstd}Travel time instead of length:{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) cost(travel_time) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)}{p_end}

{pstd}Switching road class costs 500 m; ramps are transparent:{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) classvar(highway) penalty(500) connectors(motorway_link) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park)}{p_end}


{title:Also see}

{psee}
{helpb rd2dnet}, {helpb rdnet}, {helpb rdbwnet}, {helpb rdnetdensity}
