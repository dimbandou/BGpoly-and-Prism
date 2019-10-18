# BGpoly-and-Prism
BGpoly and Prism program repository

This code purpose is to model the geometry of a subsurface body, using its gravity effect on measurement points (MPs);
It can be used with a grid of points or with a set of MPs, depending on the user input data.

To do so it allows the user to chose between two methods :
    - The 'BGPoly' method, from Talwani & Ewing 1960
    - The 'Prisma' method, from Nagy 1966

The first method will allow the user to use height lines corresponding to the geometry of the body.
Calculation of the gravity effect is done by summing the effect of each height line on the MPs,
and by interpolating between each height line.
Require the density for each height line.

the second method let the user calculate the gravity effect on the MPs, by approximating the shape of the subsurface bodies,
with prisms. Calculation of the gravity effect is done by summing the effect of each prisms on the MPs.
This method doesn't involve interpolation.
Each prism is made of 8 points, in a (x,y,z) cartesian system. 
Require the lithology to be known and the density of each layer.
