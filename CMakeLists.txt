cmake_minimum_required(VERSION 3.20)
project(Nuovo_2decomp Fortran)
enable_language(Fortran)
#
### MPI Start ##
#find_package(MPI REQUIRED)
#if(NOT MPI_Fortran_FOUND)
#    message(FATAL_ERROR "Could not find Fortran MPI.  Please set MPI_Fortran_COMPILER to point to the mpif90 wrapper.")
#endif()
#set(CMAKE_Fortran_COMPILER "mpif90")
#include_directories(${MPI_Fortran_INCLUDE_PATH})
#add_compile_options(${MPI_Fortran_COMPILE_FLAGS})
#link_directories(${MPI_Fortran_LIBRARIES})
### MPI End ##
#
#set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -cpp -fbounds-check ")#-Wall
##set (CMAKE_Fortran_FLAGS_DEBUG   "-fallow-argument-mismatch -Wall -cpp -O0 -g3 -fbounds-check -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -fopenmp -Waliasing")
#
### 2decomp Start ##
##set(F_OBJ ./2decomp_fft/src/decomp_2d.f90
##          ./2decomp_fft/src/io_mpi.f90
##          ./2decomp_fft/src/glassman.f90
##          ./2decomp_fft/src/fft_fftw3.f90)
##
##add_library(2decomp_fft STATIC ${F_OBJ})
##
##target_compile_definitions(2decomp_fft PUBLIC -DOVERWRITE -DEVEN -DDOUBLE_PREC)
##target_include_directories(2decomp_fft PUBLIC 2decomp_fft/build/opt/include)
###target_compile_options(2decomp_fft PUBLIC  -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -fopenmp -Waliasing)
##target_link_libraries(2decomp_fft PUBLIC 2decomp_fft/build/opt/lib)
##
#
### 2decomp End ##
#
#####################################################
#set(TARG "Nuovo_2decomp")
#####################################################
#set(SOURCE_FILES main.f90)
#
#if(${TARG} STREQUAL "run3d")
#    set(SOURCE_FILES_loc
#            shear.f90
#            rkuvws.f90
#            rkuvws_visc.f90
#            forcing_pv.f90
#            init_phys.f90
#            init_pjet_phys.f90
#            cut_sft_phys.f90
#            cont.f90
#            carreaustats.f90)
#elseif(${TARG} STREQUAL "post3d")
#    set(SOURCE_FILES_loc
#            post.f90
#            equilibrium.f90)
#    set(HDF5_DIR 2dechdf5_link)
#    include_directories(${HDF5_DIR}/include)
#    include_directories(${HDF5_DIR}/lib)
#    #set(FortHDF5_DIR 2dechdf5_fort_link)
#    #include_directories(${FortHDF5_DIR})
#    link_directories(${HDF5_DIR}/lib)
#    set(libreria hdf5hl_fortran
#            hdf5_fortran)
#endif()
#
#add_executable(${TARG} ${SOURCE_FILES_loc} ${SOURCE_FILES})
#target_link_libraries(${TARG} PUBLIC 2decomp_fft ${MPI_Fortran_LIBRARIES})

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
#set(2DECOMP_FFT_DIR "/home/marseo/Programmi/CLion/Progetti/Nuovo_2decomp/2decomp_fft/build/opt")
#include_directories(${2DECOMP_FFT_DIR}/include)
#link_directories(${2DECOMP_FFT_DIR}/lib)
#set(F_OBJ ./2decomp_fft/build/src/decomp_2d.f90
#          ./2decomp_fft/src/io.f90
#         ./2decomp_fft/src/glassman.f90
#         ./2decomp_fft/src/fft_fftw3.f90)
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
#target_compile_definitions(2decomp_fft PUBLIC -DOVERWRITE -DEVEN)# -DHDF5_W)#-DDOUBLE_PREC
target_compile_definitions(2decomp_fft PUBLIC -DEVEN -DOVERWRITE -DDOUBLE_PREC)
#target_compile_options(2decomp_fft PUBLIC -g -cpp -c -g3 -fbounds-check)# -Wall -Wextra -Wpedantic)
#target_compile_options(2decomp_fft PUBLIC  -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underf    low -finit-real=nan -fopenmp -Waliasing)
target_link_libraries(2decomp_fft -L2decomp_fft/build/opt/lib -L2decomp_fft -L/usr/lib/x86_64-linux-gnu -lfftw3f -lfftw3)
##########################

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
#target_link_libraries(run3d MPI::MPI_Fortran -L/usr/lib/x86_64-linux-gnu -lfftw3 fftw3 fftw3f ${FFTW_LIBRARIES} ${2DECOMP_FFT_DIR}/lib/libdecomp2d.a)
target_link_libraries(run3d MPI::MPI_Fortran -L/usr/lib/x86_64-linux-gnu -lfftw3 -lfftw3f fftw3 fftw3f ${FFTW_LIBRARIES} 2decomp_fft)