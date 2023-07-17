# Gravi3D
A 3D forward modelling software using gravity data, to resolve the geometry of subsurface objects.

Gravi3D is a program written in Python 3 for the purpose to estimate the 3-D geometry of an arbitrarily shaped subsurface body (e.g., valley overdeepenings filled with Quaternary sediments), using its gravity effect on measurement points (MPs) organized either in a profile or in a grid by forward modelling. 
The calculation of the gravity effect of an arbitrarily shaped body requires the integration of the gravity effect by an infinitesimally small volume over the full extension of the body in all three space dimensions. As there exists analytical solutions to such 3D integration for specific and regularly shaped bodies only (i.e., spheres, prisms, layers of polygonal shape, etc), the arbitrarily shaped body is approximated by a series of such regularly shaped bodies and the overall gravity effect denotes the weighted sum over all effects of these bodies.

Gravi3D comprises of two routines solving for two of those geometries. The first routine is called PRISMA and is based on the analytical solution for a prismatic body originally proposed by Nagy (1966). The second routine is called BGPoly (Bouguer Gravity Polygons) and it is based on the analytical solution for a thin layer of polygonal shape as proposed by Talwani and Ewing (1960).

The geometry of overdeepenings has been obtained with Gravi3D (Bandou et al. 2022) thanks to a workflow developed specifically for this use (Bandou, 2023)

References:

Bandou, D., Schlunegger, F., Kissling, E., Marti, U., Schwenk, M., Schläfli, P., Douillet, G., Mair, D., 2022. Three-dimensional gravity modelling of a Quaternary overdeepening fill in the Bern area of Switzerland discloses two stages of glacial carving. Scientific Reports 12, 1441. doi:10.1038/s41598-022-04830-x.

Bandou, D., 2023. Overdeepenings in the Bern region, Switzerland: Understanding their formation processes with 3D gravity forward modelling. Bern, Switzerland.

Nagy, D., 1966. The gravitational attraction of a right rectangular prism. Geophysics XXXI, 362-371. doi:10.1190/1.1439779.

Talwani, M., Ewing, M., 1960. Rapid computation of gravitational attraction of three-dimensional bodies of arbitrary shape. GEOPHYSICS 25, 203-225. doi:10.1190/1.1438687.
