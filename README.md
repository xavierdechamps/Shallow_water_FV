# Shallow_water_FV
A Fortran program that solves the 2D shallow water equations with a finite volume method based on a van Leer Q-scheme.

This program solves the two-dimensional shallow water equations with both the bed slope and bed+wall friction.
The numerical discretization is based on a finite volume method with a upwind van Leer Q-scheme used for the flux terms and for the geometrical source term.
A third order Total Variation Diminishing Runge-Kutta and a fourth order explicit Runge-Kutta time integration methods are implemented.
Steady and unsteady computations can be simulated by switching on/off a local time stepping calculation.
Further information can be found in the documentation [shallow_water.pdf](Doc/shallow_water.pdf).

## Organization

* [**Doc**](https://github.com/xavierdechamps/Shallow_water_FV/tree/main/Doc): this folder contains the PDF file that explains the theory and shows the results from the program on several test cases.
* [**Examples**](https://github.com/xavierdechamps/Shallow_water_FV/tree/main/Examples): this folder contains some test cases that are helpful to get used to the program. A PDF file explains the tests cases and the operations to be done to perform computations.
* [**SRC**](https://github.com/xavierdechamps/Shallow_water_FV/tree/main/SRC): this folder contains the source code

## Get the code

You can get the latest code by cloning the master branch:

```
git clone https://github.com/xavierdechamps/Shallow_water_FV.git
```
or by downloading it as a [zip file](https://github.com/xavierdechamps/Shallow_water_FV/archive/refs/heads/main.zip).

## Build the code

The code can be compiled by CMake as the necessary files ([CMakeLists.txt](CMakeLists.txt) and [CMake.config](CMake.config)) are provided.
The user can modify the content of [CMake.config](CMake.config) to his/her own configuration, in particular the parameters Windows, Have_SigWatch, Have_Gnuplot and Have_OpenMP.
The parameter MANGLING is required to link with the external library sigwatch. Depending on you compiler you may choose between uppercase/lowercase with addition of a trailer "_" or not. 
On Linux you can get the name mangling with the command "nm libsigwatch.a". On Windows you get the name mangling with the command "dumpbin /ALL sigwatch.lib".
The parallelism inside the subroutine [flux.f90](SRC/flux.f90) with OpenMP can be tuned by setting the number of threads to use (environment variable OMP_NUM_THREADS).

```
mkdir build
cd build
cmake ..
make
```

Hereunder an example of a supercritical flow inside a converging channel:
![supercritical flow inside a converging channel](https://github.com/xavierdechamps/Shallow_water_FV/blob/main/Doc/pics/supercritical_symmetrical_contraction_2Dsol.png)
