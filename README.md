# FVLT_PDE
Finite Volume code for PDE of living tissues

This repository contains a serie of MATLAB codes to simulate PDE models encountered in the study of living tssues.
The finite volume framework is used to discretize space.
Multiple choices are available to compute the convective fluxes: from the simple upwind method to TVD or MUSCL methods to gain high-order accuracy in space.

Development of the code is still in process.
At the moment, the code has 1D examples for the following equation:
- Cahn-Hilliard equation
- Keller-Segel model
- Porous medium equation
- Coupled system of porous medium equations
- Simple hyperbolic transport