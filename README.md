# Shallow_water_FV
Solve the 2D shallow water equations with a finite volume method based on a van Leer Q-scheme.

This program written in Fortran solves the two-dimensional shallow water equations with both the bed slope and bed friction.
The numerical discretization is based on a finite volume method with a upwind van Leer Q-scheme used for the flux terms and for the geometrical source term.
A third order Total Variation Diminishing Runge-Kutta and a fourth order explicit Runge-Kutta time integration methods are implemented.
Steady and unsteady computations can be simulated by switching on/off a local time stepping calculation.
Further information can be found in the documentation [shallow_water.pdf](Doc/shallow_water.pdf).

The code can be compiled by cmake as the necessary files ([CMakeLists.txt](CMakeLists.txt) and [CMake.config](CMake.config)) are provided.
The user has to modify the content of [CMake.config](CMake.config) to his/her own configuration.

Hereunder an example of a supercritical flow inside a converging channel:
![supercritical flow inside a converging channel](https://github.com/xavierdechamps/Shallow_water_FV/blob/main/Doc/pics/supercritical_symmetrical_contraction_2Dsol.png)
