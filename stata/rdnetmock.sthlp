{smcl}
{* *! version 0.1.0}{...}
{title:Title}

{phang}
{bf:rdnetmock} {hline 2} Load the bundled mock city


{title:Syntax}

{p 8 16 2}
{cmd:rdnetmock} [{cmd:,} {opt clear}]


{title:Description}

{pstd}
{cmd:rdnetmock} loads a small synthetic Midtown-Manhattan-style city that ships
with {cmd:rd2dnet}, so examples work offline. It is the same data as the Python
and R packages' mock city.

{pstd}
The 3,000 simulated parcels are loaded into memory (variables {cmd:x},
{cmd:y}, {cmd:treated}, {cmd:far}, {cmd:block}, {cmd:outcome}). The street
network is loaded into the frames {cmd:rdn_nodes} ({cmd:id}, {cmd:x},
{cmd:y}, {cmd:far}) and {cmd:rdn_edges} ({cmd:from}, {cmd:to}, {cmd:length},
{cmd:highway}, {cmd:speed_kph}, {cmd:travel_time}). The matrices
{cmd:rdn_line} (the boundary, 60th Street at y = 1600), {cmd:rdn_region} (the
treated zone, south of the boundary) and {cmd:rdn_park} (a park on the
boundary, to be masked) are defined.

{pstd}
Layout (metres): 11 north-south avenues 250 m apart, 41 east-west streets
80 m apart, a west-side highway joined to the grid by ramps, and a park north
of the boundary between avenues 3 and 7. The outcome is
0.0008(y - 1600) + 0.6cos(x/700) + 0.2ln(far) + tau(arc)T + N(0, 1), with
tau(arc) = 1.5 + 0.75sin(2 pi (arc - 400)/2500).


{title:Example}

{phang2}{cmd:. rdnetmock, clear}{p_end}
{phang2}{cmd:. rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)}{p_end}


{title:Also see}

{psee}
{helpb rd2dnet}, {helpb rdnetset}
