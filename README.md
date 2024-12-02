# SPECTRAL FULLY PERIODIC CODE

This README contains basic instructions for building the code and installing the 2DECOMP&FFT library.

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

## Building

The build system is driven either by `CMake` or `Makefile`. However, the `Makefile` uses some CMake command line so it 
must be installed.

### 2decomp&FFT
Here a default procedure is described, for more information on how to build the 2DECOMP&FFT library 
please see [INSTALL.md](INSTALL.md) in the 2decomp-fft folder.

To download the library:
```
git clone https://github.com/2decomp-fft/2decomp-fft.git                                                                                                                                  
```
and an example of compilation and building could be:
```
cmake -S . -B ./build -DEVEN=ON -DOVERWRITE=ON -DDOUBLE_PREC=ON -DFFT_Choice=fftw_f03 -DFFTW_ROOT=/usr/local
```
Additional flags (that should be activated by default):
```
 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=./build/opt
```
in the directoy where the code has been build (e.g. ./build) type the last 2 commands to complete the process:
``` 
cmake --build .
cmake --install .
```
The default location for `libdecomp2d.a` is `$path_to_build_directory/opt/lib` unless the variable `CMAKE_INSTALL_PREFIX` is modified.
The module files generated by the build process will similarly be installed to `$path_to_build_directory/opt/include`, users of the library should add this to the include paths for their program.

Occasionally a clean build is required, this can be performed by running
```
$ cmake --build $path_to_build_directory --target clean
```
### Spectral code

Completed the 2DECOMP&FFT library install to build the code ... 

### Contributing

If you would like to contribute to the development of the 2DECOMP&FFT library or report a bug please refer to 
the [Contributing section](Contribute.md)
