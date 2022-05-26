# alpha-damping-example

Viscous flux discretization

Historically, the discretization of the viscous terms in the Navier-Stokes equations has received less attention in the literature than the discretization of the inviscid terms. Much of the research is focused on the discretization of the inviscid terms especially to reduce spurious oscillations near discontinuities. The importance of high-frequency damping in high-order conservative finite-difference schemes for viscous terms in the Navier-Stokes equations.

Fourth-order finite-difference viscous scheme for viscous fluxes and MP5 for inviscid fluxes.

https://arxiv.org/abs/2205.01034

In this paper I have published 6th-order scheme
https://doi.org/10.1016/j.jcp.2022.111195


A poor viscous scheme can lead to oscillations and failure of the simulations (and it does lot more...will be explained in future papers)

This code kind of reproduces the resutls in the above (first) paper but it may or may not work if you change the parameters or the test case. Its upto you.

Sainadh Chamarthi
