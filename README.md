# SPECTRAL-3D

This README contains basic instructions for building the `Spectral-3D` code and installing the 2DECOMP&FFT library.

## General Description

Spectral code for 3D simulation of turbulent Newtonian incompressible flows in a periodic parallelepipedal domain 
of equally spaced points. It solves the PDEs using the pseudo-spectral approach. In short, it projects the PDEs in the Fourier space transforming the equations 
from differential to algebraic. 
\
It makes extensive use of the Discrete Fast Fourier Transform with the fftw 
implementation in an MPI parallel environment. The domain is decomposed in pencil aligned along the x-direction in the physical space and along the 
z-direction in the Fourier space. 
\
The flow supported in the current version are:
- Homogeneus Isotropic Turbulence (HIT) with and without forcing;
- Temporally evolving Planar Jets.

The executable `run3d` can perform simulation in single or double precision. All inputs are taken form the `THI_PJET` file.


```
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*    A NAVIER-STOKES SOLVER FOR 3D FULLY PERIODIC COMPUTATIONS        *
*    USING MPI/2DECOMP WITH X (Real) TO Z (Spectral) PENCIL           *
*    ARRANGEMENTS IN A PARALLELEPIPEDAL DOMAIN OF EQUI-SPACED POINTS  *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *    

 Example of a X - pencil decomposition with P_row=2, P_col=3 (nprocess=2x3=6)

                          PHYSICAL SPACE               :                       FOURIER SPACE
                 ----------------------------------    :              ----------------------------------
                /           /          /          /|   :             /           /          /          /|
               /           /          /          / |   :            /           /          /          / |
              /           /          /          /  |   :           /           /          /          /  |
             /           /          /          /   |   :          /           /          /          /   |
            /           /          /          /    |   :         /           /          /          /    |
           /           /          /          /     |   :        /           /          /          /     |
    Y  ^  /           /          /          /     /|   : Y  ^  /           /          /          /     /|
       | /           /          /          /     / |   :    | /           /          /          /     / |
       |/           /          /          /     /  |   :    |/           /          /          /     /  |
       |----------------------------------     /   |   :    |----------------------------------     /   |
       |           |          |          |    /    |   :    |           |          |          |    /    |
       |           |          |          |   /     |   :    |           |          |          |   /     |
       |     3     |     4    |     5    |  /     /    :    |     3     |     4    |     5    |  /     /
       |           |          |          | /     /     :    |           |          |          | /     /
       |-----------|----------|-----------/     /      :    |-----------|----------|-----------/     /
       |    /      |          |          |     /       :    |   /       |          |          |     /
       |   /       |          |          |    /        :    |  /        |          |          |    /
       |  /   0    |     1    |     2    |   /         :    | /   0     |     1    |     2    |   /
       | /         |          |          |  /          :    |/          |          |          |  /
       X___________|__________|__________|_/____>      :    Z___________|__________|__________|_/____>
    X                                             Z    :  Z                                             X
```
## Project Structure

The project has the following structure:

- **README.md**: The file you are reading now, which provides information about the project.
- **build.conf**: The file with compilation options for the `Spectral-3D` code (virtually the only file that matters for normal use).
- **src/**: Contains the source code of the project.
  - **main.f90**: Main Fortran file of the program.
  - **m_\*.f90**: Modules files with the needed procedures.
- **bin/**: Where the executable is being built, together the input `THI_PJET` file.
- **dependencies/**: Contains the additional libraries needed.
  - **2decomp-fft**: Provides access to a pencil decompositions as well as FFTs. 
    - **Makefile**: The file with the compilation options for the 2DECOMP&FFT library (virtually the only file that matters for normal use).

## Building

The code build is driven by `make` compilation. 

To get the code from Github:
```
git clone --recurse-submodules https://github.com/fortaia/Spectral-3D.git
```
and the code will be downloaded together with its dependencies. Key compiling parameters can be modified in 
the `build.conf` file (for Spectral-3D) and `dependencies/2decomp-fft/Makefile` (for 2DECOMP&FFT library).
Once the configuration files mentioned are ready first compile the library:
```
make libs
```
and then the code:
```
make
```
the executable `run3d` will be placed in the `bin` folder together with the input file `THI_PJET`.
If a clean of `Spectral-3D` is needed:
```
make clean
```
will remove all the compilation files. Occasionally, a clean of the 2DECOMP&FFT library might be
necessary:
```
make libsclean
```
## Running
The `Spectral-3D` code makes extensive use of MPI. To run the executable navigate to the folder where the `run3d` 
file is and use the appropriate MPI executor. A typical command is:
```
mpirun -np $number_of_cores_needed ./run3d
```
note that the input file `THI_PJET` must be in the same folder of the executable. All the simulation parameters are 
specified in the `THI_PJET` file. The file is extensively commented and it is meant to be self-explanatory
regarding the several options available.

## Contributing
The recommended strategy to contribute is to start with a [discussion](https://github.com/fortaia/Spectral-3D.git/discussions) 
or to pick an existing issue. To modify or experiment with the code, fork the `Spectral-3D` github repository and 
commit changes in a dedicated branch of your fork. When the modification is ready for review, one can 
open a pull request. Pull requests must be focused, small, coherent and have a detailed description.
The longer the pull request, the harder the review. Please empathise with your fellow contributors who are going to spend time reviewing your code.


