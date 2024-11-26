cmake_minimum_required(VERSION 3.20)
project(Nuovo_2decomp Fortran)
enable_language(Fortran)

# Set the version of the project
set(PROJECT_VERSION "1.0")

# Find MPI
find_package(MPI REQUIRED)
# Set the Fortran compiler to the MPI wrapper
set(CMAKE_Fortran_COMPILER mpif90)

# Set paths to FFTW manually
set(FFTW_INCLUDES "/usr/include")
set(FFTW_LIBRARIES "/usr/lib/x86_64-linux-gnu/libfftw3.so;/usr/lib/x86_64-linux-gnu/libfftw3f.so")

# Add FFTW include directories
include_directories(${FFTW_INCLUDES})

# Set compiler flags for Fortran (you can adjust optimization/debugging flags here)
set(CMAKE_Fortran_FLAGS "-O0 -cpp -Wall -fbacktrace -fbounds-check -ffree-line-length-none")
#set(CMAKE_Fortran_FLAGS "-O3 -cpp -ffree-line-length-none")

# Set paths for 2decomp&fft library
set(F_OBJ ./2decomp_fft/build/src/decomp_2d_constants.f90
        ./2decomp_fft/build/src/decomp_2d_mpi.f90
        ./2decomp_fft/build/src/profiler_none.f90
        ./2decomp_fft/build/src/factor.f90
        ./2decomp_fft/build/src/decomp_2d.f90
        ./2decomp_fft/build/src/glassman.f90
        ./2decomp_fft/build/src/fft_fftw3.f90
        ./2decomp_fft/build/src/fft_log.f90
        ./2decomp_fft/build/src/interp.f90
        ./2decomp_fft/build/src/io_utilities.f90
        ./2decomp_fft/build/src/io_adios_none.f90
        ./2decomp_fft/build/src/io_object_mpi.f90
        ./2decomp_fft/build/src/io_mpi.f90
        ./2decomp_fft/build/src/log.f90
        ./2decomp_fft/build/src/transpose_x_to_y.f90
        ./2decomp_fft/build/src/transpose_y_to_x.f90
        ./2decomp_fft/build/src/transpose_y_to_z.f90
        ./2decomp_fft/build/src/transpose_z_to_y.f90)
add_library(2decomp_fft STATIC ${F_OBJ})
target_compile_definitions(2decomp_fft PUBLIC -DEVEN -DDOUBLE_PREC)
target_link_libraries(2decomp_fft -L2decomp_fft/build/opt/lib -L2decomp_fft -L/usr/lib/x86_64-linux-gnu -lfftw3f -lfftw3)

# Create an executable target for your program
set(SOURCES
        m_glob_params.f90
        m_aux_phys.f90
        m_aux_spect.f90
        m_utils.f90
        m_io.f90
        m_initial_cond.f90
        m_equations_terms.f90
        m_core_calc.f90
        main.f90
)
add_executable(run3d ${SOURCES})

# Link against MPI, FFTW, and 2decomp&fft libraries
target_link_libraries(run3d MPI::MPI_Fortran -L/usr/lib/x86_64-linux-gnu -lfftw3 -lfftw3f fftw3 fftw3f ${FFTW_LIBRARIES} 2decomp_fft)
