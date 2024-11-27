program Spectral_3D
    
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !       A NAVIER-STOKES SOLVER FOR 3D FULLY PERIODIC COMPUTATIONS      *
    !       USING MPI/2DECOMP WITH X (Real) TO Z (Spectral) PENCIL         *
    !       ARRANGEMENTS                                                   *
    !                                                                      *
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    ! ======================================================================
    !
    !       nrank numbers (0=master) - PROCESS ORGANIZATION
    !
    !       Example of a X - pencil decomposition
    !       with P_row=2, P_col=3 (nprocess=2x3=6)
    !
    ! ======================================================================
    !
    !                          PHYSICAL SPACE
    !
    !                 ----------------------------------
    !                /           /          /          /|
    !               /           /          /          / |
    !              /           /          /          /  |
    !             /           /          /          /   |
    !            /           /          /          /    |
    !           /           /          /          /     |
    !    Y  ^  /           /          /          /     /|
    !       | /           /          /          /     / |
    !       |/           /          /          /     /  |
    !       |----------------------------------     /   |
    !       |           |          |          |    /    |
    !       |           |          |          |   /     |
    !       |     3     |     4    |     5    |  /     /
    !       |           |          |          | /     /
    !       |-----------|----------|-----------/     /
    !       |    /      |          |          |     /
    !       |   /       |          |          |    /
    !       |  /   0    |     1    |     2    |   /
    !       | /         |          |          |  /
    !       X___________|__________|__________|_/____>
    !    X                                             Z
    !
    !
    !
    !                          FOURIER SPACE
    !
    !                 ----------------------------------
    !                /           /          /          /|
    !               /           /          /          / |
    !              /           /          /          /  |
    !             /           /          /          /   |
    !            /           /          /          /    |
    !           /           /          /          /     |
    !    Y  ^  /           /          /          /     /|
    !       | /           /          /          /     / |
    !       |/           /          /          /     /  |
    !       |----------------------------------     /   |
    !       |           |          |          |    /    |
    !       |           |          |          |   /     |
    !       |     3     |     4    |     5    |  /     /
    !       |           |          |          | /     /
    !       |-----------|----------|-----------/     /
    !       |           |          |          |     /
    !       |           |          |          |    /
    !       |     0     |     1    |     2    |   /
    !       |           |          |          |  /
    !       Z___________|__________|__________|_/____>
    !    Z                                             X
    !
    !
    ! ///////////////////////////////////////////////////// !

    use decomp_2d_constants
    use decomp_2d_mpi
    use decomp_2d_fft
    use decomp_2d
    use MPI

    use m_glob_params
    use m_io, only: Read_THI_PJET, Write_3d_velocities
    use m_initial_cond
    use m_core_calc, only: Core_calc

    implicit none
    
!!!!!!!!!!!!!!!!!!! Variable Declaration START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !****** 2DECOMP *******
    ! LOCAL GRID SIZE (NOTE: row*col = numprocs when running.)
    ! NOTE: 0*0 SELF SELECTS (MOST LIKELY) BEST GRID - AUTO-TUNNING
    ! INTEGER :: P_row , P_col
    type(decomp_info), pointer :: ph => null(), sp => null()

    ! LOCAL GRID SIZE (NOTE: row*col = numprocs when running.)
    ! NOTE: 0*0 SELF SELECTS (MOST LIKELY) BEST GRID - AUTO-TUNNING
    ! INTEGER :: P_row , P_col
    logical, dimension(3) :: periodic_bc
    
!!!!!!!!!!!!!!!!!!! Variable Declaration END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! MPI INITIALIZATION
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierr)
    
    if (nrank==0) then
        print *, " ------------------"
        print *, "|   MPI Started    |"
        print *, " ------------------"
    endif
    
    ! ----- Read THI_PJET -------
    call Read_THI_PJET
    
    if (nrank==0) then
      print *, " ----------------------"
      print *, "| Initialising 2DECOMP |"
      print *, " ----------------------"
    endif
    
    ! 2DECOMP initialization
    periodic_bc=.true.
    ! Initialize and use the 2decomp_fft library
    call decomp_2d_init(N1G, N2G, N3G, P_row, P_col, periodic_BC)
    
    ! Initialises 2decomp ffts (size of local 'pencils' for physical and spectral (SLAVES) arrays)
    call decomp_2d_fft_init(PHYSICAL_IN_X)
    ! Pointers for pencils dimensions
    ph => decomp_2d_fft_get_ph()
    sp => decomp_2d_fft_get_sp()
    
    
    if (nrank==0) then
        
        write(*,*)
        write(*,*) '////////////////////////////////////////////////////'
        write(*,*)
        write(*,*) '  A 3D PSEUDO-SPECTRAL CODE FOR                     '
        write(*,*) '  NAVIER-STOKES SIMULATIONS                         '
        write(*,*)
        write(*,*) '  PARALLELIZED WITH:  MPI AND 2DECOMP'
        write(*,*)
        write(*,'("   NUMBER OF PROCS: ",I5)') nproc
        write(*,*)
        write(*, '("   P_ROW x P_COL = ",I5,"  x",I5)') P_ROW, P_COL
        write(*,*)
        write(*, '("   GRID SIZE : ",I5,"  x  ",I5,"  x  ",I5)') N1G,N2G,N3G
        write(*,*)
        write(*,*) '////////////////////////////////////////////////////'
    
    endif
    
    ! Initial Conditions: where u,v,w: cmplx and ur,vr,wr: real are allocated
    call Initial_cond(ph,sp)
    
    ! Calculations
    call Core_calc(ph,sp)
    
    ! Writing final files
    call Write_3d_velocities(trim(CURDAT),u,v,w,sp)
    
    if (nrank==0) then
        
        write(*,*) '----------------------------------------------------'
        write(*,*)
        write(*,*) '  RUN COMPLETED SUCCESSUFULLY                       '
        write(*,*)
        write(*,*) '----------------------------------------------------'
    
    endif
    
    deallocate(u,v,w)
    nullify(ph)
    nullify(sp)
    call decomp_2d_fft_finalize
    call decomp_2d_finalize

    call MPI_FINALIZE(ierr)
end program
