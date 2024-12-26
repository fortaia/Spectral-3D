!
! MODULE purpose: Group input/output procedures
!
module m_io
  
  implicit none
  
  public
  
  ! Generic interface to handle multiple data types
  interface Read_3d_field
    module procedure Read_3d_field_cmplx
    module procedure Read_3d_field_real
  end interface Read_3d_field

contains
  !*************************************************************
  !
  ! Read values from THI_PJET file
  !
  !*************************************************************
  subroutine Read_THI_PJET()
    !
    !*************************************************************
    !
    use m_glob_params
    use decomp_2d_mpi, only: nrank, mytype, real_type
    use MPI
    
    implicit none
    
    character(50) :: alist
    integer :: nnf,errorcode
    real(mytype) :: dj2
    real(mytype), parameter :: TWOPI=6.28318530717958647692528676655900
    logical :: file_exists
    
    !---------   Calculations start ------------------
    if (nrank==0) then
      
      ! Check if the file exists
      inquire(file="THI_PJET", exist=file_exists)
      
      if (.not. file_exists) then
        print *, "THI_PJET NOT FOUND."
        call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
      end if
      print *, ' '
      print *, '** Reading THI_PJET **'
      print *, ' '
      nnf = 55
      OPEN (nnf, FILE = 'THI_PJET', FORM = 'FORMATTED', STATUS = 'OLD')
      
      !     PARAMETERS FIXING THE NATURE OF THE BOUNDARY CONDITIONS
      read(nnf, *) NBOUND
      
      !     NUMBER OF PENCILS ALONG Y AND Z
      !     P_row: number of pencils along y (Nmytypeocs=P_row * P_col)
      read(nnf, *) P_row
      !     P_col: number of pencils along z (Nmytypeocs=P_row * P_col)
      read(nnf, *) P_col
      !     Random Generator Option
      read(nnf, *) ROPTION
      
      !     NUMBER OF (GLOBAL) PHYSICAL GRID POINTS ALONG X,Y, AND Z
      read(nnf, *) N1G
      read(nnf, *) N2G
      read(nnf, *) N3G
      
      !     TRUNCATION OUTSIDE OF A SPHERE OF RADIUS (ALPHA0*NMAX/2)
      read(nnf, *) ALPHA0
      !     NUMBER OF TIMES THE RANDOM NUMBER GENERATOR IS CALLED
      read(nnf, *) NRANDCALL
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     NRKU: ORDER OF Runge-KuttA SCHEME
      read(nnf, *) NRKU
      !     TIME STEP
      read(nnf, *) DT
      !     MAXIMUM NUMBER OF ITERATIONS
      read(nnf, *) NTMAX
      !     ICFL: IF(ICFL==1), THE TIME STEP IS
      !     COMPUTED BY MEANS OF A CFL CONDITION
      read(nnf, *) ICFL
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     IF ISIMU = 0 HOMOGENEOUS ISOTROPIC TURBULENCE
      !     IF ISIMU = 1 TEMPORAL PLANE JET
      !     IF ISIMU = 2 SHEAR-FREE TURBULENCE
      !     IF ISIMU = 3 TEMPORAL PLANE WAKE
      !     IF ISIMU = 4 TEMPORAL ROUND JET
      read(nnf, *) ISIMU
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     IF ISCAR = 1  PASSIVE SCALAR  COMPUTATION TURNED ON
      !     IF ISCAR = 2  PASSIVE SCALAR USED TO COMPUTE Ksgs
      read(nnf, *) ISCAR
      !     SCHMIDT NUMBER
      read(nnf, *) SCH
      !     RNU MOLECULAR VISCOSITY
      read(nnf, *) RNU
      !     IF (ICHOPA==1) S.G.S MODELLING (Lambalais)
      read(nnf, *) ICHOPA
      !     IF(ISTRUC==1): Structure function SGS
      read(nnf, *) ISTRUC
      !     IF (ISMAG==1): Smagorinsky SGS
      read(nnf, *) ISMAG
      read(nnf, *) IHVISC
      !     NEXPH       HYPER-VISCOSITY ORDER
      !    (The Laplacian diffusive operater is iterated NEXPH times)
      read(nnf, *) NEXPH
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     NAME OF THE CURRENT DATASET DISPOSED
      read (nnf, *) CURDAT
      !     ACQUIRE DATA FROM mytypeEVIOUS RUN TO INITIALIZE PRESENT RUN
      read (nnf, *) RESDAT
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     IF JINIT=1 : ACQUIRE DATA FROM PREVIOUS RUN
      !     IF JINIT=0 : INITIAL CONDITIONS DEFINED BY SUB INIT
      read(nnf, *) JINIT
      !     readING OF RECORD NUMBER NI OF CURDAT FOR POST-PROCESSING
      read(nnf, *) NI
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     IF JWRITE=1: THE FIELDS ARE WRITTEN AT EACH JSAVE TIME STEPS
      !     WARNING: JSAVE MUST DIVIDE NAVG
      read(nnf, *) JWRITE
      read(nnf, *) JSAVE
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     WRITE OF STATISTICS AT EACH NOUT TIME STEPS
      read(nnf, *) NOUT
      !     WRITE OF SPECTRA  AT EACH NOUTT  TIME STEPS
      read(nnf, *) NOUTT
      !     IF ImytypeESC==1 COMPUTE PRESSURE STATISTICS
      read(nnf, *) ImytypeESC
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     IPDF=1: PDF COMPUTED
      read (nnf, *) IPDF
      !     NOPDF  : PDF COMPUTED EVERY NOPDF TIME-STEP IF ICFL==0
      !              PDF COMPUTED EVERY NOPDF Tret      IF ICFL==1
      read (nnf, *) NOPDF
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     PARAMETERS USED IN INITIALISING FIELDS
      read(nnf, *) C1
      read(nnf, *) NSLOPE1
      read(nnf, *) C4
      read(nnf, *) NSLOPE4
      !     TOTE IS THE RMS VALUE OF ANY OF THE THREE VELOCITY
      !     COMPONENTS
      read(nnf, *) TOTE
      !     TOTS IS THE VARIANCE OF THE TEMPERATURE
      !     (OR PASSIVE SCALAR)
      read(nnf, *) TOTS
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     read THE PARAMETERS FOR THE SCALAR FORCING FUNCTION
      !     IFORCE: if =1 scalar forcing turned on using Alvelius method (Phys. Fluids, 1999)
      !     IFORCE: if =2 scalar forcing turned on using mean scalar gradient (Overholt and Pope, Phys. Fluids, 1996)
      !     read THE PARAMETERS FOR THE FORCING FUNCTION
      !     IFORCE: if =1 forcing turned on using Alvelius method (Phys. Fluids, 1999)
      !     IFORCE: if =2,3,4... (other methods not coded)
      read(nnf, *) IFORCE
      !     X0: Forcing wavenumber
      read(nnf, *) X0
      !     TOTF: FORCING  INTENSITY
      read(nnf, *) TOTF
      !     WIDTH_F: NUMBER OF MODES WHERE VELOCITY FORCE IS CONCENTRATED
      read(nnf, *) WIDTH_F
      !     read THE PARAMETERS FOR THE SCALAR FORCING FUNCTION
      !     IFORCE: if =1 scalar forcing turned on using Alvelius method (Phys. Fluids, 1999)
      !     IFORCE: if =2 scalar forcing turned on using mean scalar gradient (Overholt and Pope, Phys. Fluids, 1996)
      read(nnf, *) ISFORCE
      !     X0S: Scalar Forcing wavenumber
      read(nnf, *) X0S
      !     TOTFS: SCALAR FORCING INTENSITY
      read(nnf, *) TOTFS
      !     WIDTH_FS: NUMBER OF MODES WHERE SCALAR FORCE IS CONCENTRATED
      read(nnf, *) WIDTH_F_SCAL
      !     BETA_SCAL: MEAN SCALAR GRADIENT ALONG Y-DIRECTION
      read(nnf, *) BETA_SCAL
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     read THE PARAMETERS FOR THE PLAN JET
      !     (ONLY IF ISIMU==1)
      !     REYNOLDS NUMBER FOR THE PLANE JET
      !    (ONLY IN PLANE JET COMPUTATIONS IT IS read HERE)
      read(nnf, *) REYH
      !     A0X IS THE RATIO LX/LY
      read(nnf, *) A0X
      !     A0Z IS THE RATIO LX/LZ
      read(nnf, *) A0Z
      !     BLY IS THE BOX SIZE (ALONG Y DIRECTION) IN H/2 UNITS
      read(nnf, *) BLY
      !     HUTHETA IS THE RATIO OF H (SLOT-WIDTH) TO THE INIIAL MOMENTUM
      !     THICKNESS
      !     H/THETA FOR THE INITIAL VELOCITY PROFILE
      read(nnf, *) HUTHETA
      !     U2 IS THE MAXIMUM INITIAL JET VELOCITY
      read(nnf, *) U2
      !     U1 IS THE MINIMUM INITIAL JET VELOCITY
      read(nnf, *) U1
      !     SHTHETA IS THE RATIO OF H (SLOT-WIDTH) TO THE INIIAL SCALAR
      !     MOMENTUM THICKNESS
      !     H/THETA FOR THE INITIAL SCALAR PROFILE
      read(nnf, *) HSTHETA
      !     T2 IS THE MAXIMUM INITIAL SCALAR FIELD
      read(nnf, *) T2
      !     T1 IS THE MINIMUM INITIAL SCALAR FIELD
      read(nnf, *) T1
      !     COMPUTE MOLECULAR VISCOSITY (RNU)
      !    (ONLY IN PLANE JET COMPUTATIONS IT IS COMPUTED AGAIN HERE)
      if (ISIMU==1.OR.ISIMU==3.OR.ISIMU==4) then
        !     (H IS THE EQUIVALENT TO THE ROUND JET DIAMETER)
        dj2 = TWOPI / BLY
        !     HUTHETA IS THE RATIO OF H (SLOT-WIDTH) TO THE INIIAL MOMENTUM
        !     THICKNESS H/THETA FOR THE INITIAL VELOCITY mytypeOFILE
        RNU = (U2 - U1) * dj2 / REYH
      end if
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     NJJ IS THE SIZE OF THE DOMAIN SET TO ZERO IN RDT SIMUATIONS
      read(nnf, *) NJJ ! NOT ACTIVE(!). USING BLY, HUTHETA and SHTHETA ABOVE
      !     NNRDT=1 Start RDT simulation from THI; =2 Re-Start RDT
      read(nnf, *) NNRDT
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     IPARTIC:
      !     IPARTIC=1: COMPUTE PARTICLE TRACKING
      !    (IPARTIC=0: NO PARTICLE TRACKING)
      read(nnf, *) IPARTIC
      !     INTERPOL_TYPE: 0=EXACT; 1=LINEAR; 2=3CUB_FFT; 3=3CUB_DIFF; 4=TS13
      read(nnf, *) INTERPOL_TYPE
      !     NPARTIC IS THE TOTAL NUMBER OF PARTICLES (GLOBAL)
      read(nnf, *) NPARTIC_G
      !     PARTIC_LOC IS THE INITIAL POSITION OF PARTICLES
      !    (PARTIC_LOC=1: ALL DOMAIN
      !     PARTIC_LOC=2: SHEAR LAYER
      !     PARTIC_LOC=3: OUTSIDE SHEAR LAYER)
      read(nnf, *) PARTIC_LOC
      !     NPART_VAR: VARIABLES TO WRITE ON PARTICLES DURING RUN (=0 usual i.e. x,y,z,u,v,w; NE.0 e.g.3 read 3 variables from THI_PJET_partic
      read(nnf, *) NPART_VAR
      !     WRITE PARTIC DATA EVERY NOUTPART TIME
      read(nnf, *) NOUTPART
      !     (=0: INITIAL CONDITIONS DEFINED BY PREPART; =1 RESTART)
      !     IF JPINIT=1 : ACQUIRE DATA FROM PREVIOUS RUN
      !     IF JPINIT=0 : INITIAL CONDITIONS DEFINED BY PREPART
      read(nnf, *) JPINIT
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     IVISC: =0 NEWTONIAN; =1 VISCOELASTIC
      read(nnf, *) IVISC
      !     LPOLY/AAA  : LENGTH OF POLYMER MOLECULE (IF IVISC==3: JOHNSON-SEGALMAN LPOLY=aaa)
      !     aaa=1 : OLDROYD-B
      !     aaa=-1: OLDROYD-A
      !     aaa=0 : Co-rotational model
      read(nnf, *) LPOLY
      !     RLTIME (LAMBDA) : RELAXATION TIME
      read(nnf, *) RLTIME
      !     POLY_BETA : POLYMER CONCENTRATION
      read(nnf, *) POLY_BETA
      !     LAMBDA CARREAU MODEL
      read (nnf, *) LAMPOWER
      !     NPOWER CARREAU MODEL
      read (nnf, *) NPOWER
      !     APOWER CARREAU MODEL
      read (nnf, *) APOWER
      !     NUZERO CARREAU MODEL
      read (nnf, *) NUZERO
      !     NUINFTY CARREAU MODEL
      read (nnf, *) NUINFTY
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     For Buiding Animation
      read(nnf, *) BUILDANIMATION
      read(nnf, *) TIMESTART
      read(nnf, *) TIMESKIP
      read(nnf, *) NUMBERFRAMES
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     For Post mytypeocessing Writing
      read(nnf, *) WRITEBINARY
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     FOR WRITTING PARTIAL DATA FOR FURTHER PROCESSING
      !     NDIMS: =0 NO PARTIAL OUTPUT FIELD; =1 1D; =2 2D; =3 3D
      read(nnf, *) NDIMS
      !     ISTART : I START INDEX FOR PARTIAL FIELD
      read(nnf, *) ISTART
      !     IEND   : I END INDEX FOR PARTIAL FIELD
      read(nnf, *) IEND
      !     JSTART : J START INDEX FOR PARTIAL FIELD
      read(nnf, *) JSTART
      !     JEND   : J END INDEX FOR PARTIAL FIELD
      read(nnf, *) JEND
      !     KSTART : K START INDEX FOR PARTIAL FIELD
      read(nnf, *) KSTART
      !     KEND   : K END INDEX FOR PARTIAL FIELD
      read(nnf, *) KEND
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(nnf, *) ALIST
      !   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !     IMHD: =1: MHD TURNED ON (MAGNETOHYDRODYNAMICS AT LOW MAGNETIC REYNOLDS NUMBER)
      read(nnf, *) IMHD
      !     INITIAL GUESS FOR INTEGRAL SCALE (FROM PREVIOUS DNS WITH IMHD=0)
      read(nnf, *) L0_E
      !     INITIAL GUESS FOR ROOT-MEAN-SQUARE VELOCITY (FROM PREVIOUS DNS WITH WITH IMHD=0)
      read(nnf, *) U0_E
      !     MAGNETIC INTERACTION NUMBER
      read(nnf, *) N0_E
      !     ELECTRICAL CONDUTIVITY
      read(nnf, *) SIGMA_E
      
      ! @@@@@@@@@@@@@@@@@@@@ CHECKS @@@@@@@@@@@@@@@@@
      !     CHECKS FOR SHEAR-FREE TURBULENCE
      IF ((ISIMU==1.OR.ISIMU==2.OR.ISIMU==3.OR.ISIMU==4).AND. &
              (IFORCE==1.OR.ISFORCE==1)) THEN
        IF (nrank==0) THEN
          print *, ' ERROR IN FILE THI_PJET FOR ISIMU=1,2,3 OR 4'
          print *, ' FORCING IS ACTIVE (IFORCE OR ISFORCE)'
          print *, ' FORTRAN STOP'
        ENDIF
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
      ENDIF
      
      !Checks domain SIZE
      if (P_row/=0.AND.P_col/=0) then
        if(mod(N2G, P_row)/=0.OR.mod(N3G, P_col)/=0)THEN
          print *, ' ERROR IN FILE THI_PJET'
          print *, ' N2G AND N3G SHOULD BE DIVISIBLE BY P_ROW AND P_COL, RESPECTEVELY'
          print *, ' FORTRAN STOP'
          
          call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
        endif
      endif
    
    endif
    
    call MPI_BCAST(NBOUND, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(N1G, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(N2G, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(N3G, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(P_ROW, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(P_COL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(ROPTION, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(ALPHA0, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NRANDCALL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NRKU, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(DT, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NTMAX, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(ICFL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ISIMU, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ISCAR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(SCH, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(RNU, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ICHOPA, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(ISTRUC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ISMAG, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(IHVISC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NEXPH, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(CURDAT, 20, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(RESDAT, 20, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(JINIT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NI, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(JWRITE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(JSAVE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NOUT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NOUTT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(ImytypeESC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(IPDF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NOPDF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(c1, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nslope1, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(c4, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nslope4, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(TOTE, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(TOTS, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(IFORCE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(X0, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(TOTF, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(WIDTH_F, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ISFORCE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(X0S, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(TOTFS, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(WIDTH_F_SCAL, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(BETA_SCAL, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(REYH, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(A0X, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(A0Z, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(BLY, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(HUTHETA, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(U2, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(U1, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(HSTHETA, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(T2, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(T1, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(NJJ, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NNRDT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(IPARTIC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(INTERPOL_TYPE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NPARTIC_G, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(PARTIC_LOC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NPART_VAR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(NOUTPART, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(JPINIT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(IVISC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(LPOLY, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(RLTIME, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(POLY_BETA, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(LAMPOWER, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NPOWER, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(APOWER, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NUZERO, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NUINFTY, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(BUILDANIMATION, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(TIMESTART, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(TIMESKIP, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(NUMBERFRAMES, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(WRITEBINARY, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(NDIMS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ISTART, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(IEND, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(JSTART, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(JEND, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(KSTART, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(KEND, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(IMHD, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(L0_E, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(U0_E, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(N0_E, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(SIGMA_E, 1, real_type, 0, MPI_COMM_WORLD, ierr)
  
  end subroutine Read_THI_PJET
  
  !*************************************************************
  !
  ! Read 3d array file: REAL
  !
  !*************************************************************
  subroutine Read_3d_field_real(var,filename,ph)
    !
    !*************************************************************
    !
    use decomp_2d_io, only: decomp_2d_read_one
    use decomp_2d_mpi, only: mytype,nrank,decomp_2d_abort
    use decomp_2d, only: decomp_info
    
    real(mytype), dimension(:,:,:) :: var
    character(len=*), intent(IN) :: filename
    logical :: file_exists
    type(decomp_info), pointer :: ph
    
    !---------   Calculations start ------------------
    
    if (nrank == 0) then
      
      inquire (file=filename, exist=file_exists)
      
      if (.not. file_exists) then
        call decomp_2d_abort(1, "ERROR: data -> "//filename//" NOT found in folder!")
      end if
      
      print *, " Reading file : "//filename
    end if
    
    call decomp_2d_read_one(1, var, filename,opt_decomp=ph)
  
  end subroutine Read_3d_field_real
  !*************************************************************
  !
  ! Read 3d array file: COMPLEX
  !
  !*************************************************************
  subroutine Read_3d_field_cmplx(var,filename,sp)
    !
    !*************************************************************
    !
    use decomp_2d_io, only: decomp_2d_read_one
    use decomp_2d_mpi, only: mytype,nrank,decomp_2d_abort
    use decomp_2d, only: decomp_info
    
    complex(mytype), dimension(:,:,:) :: var
    character(len=*), intent(IN) :: filename
    logical :: file_exists
    type(decomp_info), pointer :: sp
    !---------   Subroutine start ------------------
    
    if (nrank == 0) then
      
      inquire (file=filename, exist=file_exists)
      
      if (.not. file_exists) then
        call decomp_2d_abort(1, "ERROR: data -> "//filename//" NOT found in folder!")
      end if
      
      print *, " Reading file : "//filename
    
    end if
    
    call decomp_2d_read_one(3,var,filename,opt_decomp=sp)
  
  end subroutine Read_3d_field_cmplx
  
  subroutine Write_3d_velocities(file_init,v1,v2,v3,sp)
    
    use decomp_2d_io, only: decomp_2d_io_init,decomp_2d_write_one,decomp_2d_io_fin
    use decomp_2d_mpi, only: mytype
    use decomp_2d, only: decomp_info
    
    complex(mytype), dimension(:,:,:) :: v1,v2,v3
    character(len=*), intent(IN) :: file_init
    character(len=len(file_init)+3) :: filename
    type(decomp_info), pointer :: sp
    !---------   Calculations start ------------------
    
    call decomp_2d_io_init()
    filename=trim(file_init//".ux")
    call decomp_2d_write_one(3,v1,filename,opt_decomp=sp)
    filename=trim(file_init//".uy")
    call decomp_2d_write_one(3,v2,filename,opt_decomp=sp)
    filename=trim(file_init//".uz")
    call decomp_2d_write_one(3,v3,filename,opt_decomp=sp)
    call decomp_2d_io_fin
  
  end subroutine Write_3d_velocities
  !
  ! Write store*.info file -> USE ONLY 1 PROCESSOR
  !
  subroutine Write_store_file(filename,time,tref,dt)
    
    use decomp_2d_mpi, only: mytype
    
    character(len=*), intent(IN) :: filename
    real(mytype), intent(IN) :: time,tref,dt
    ! -------- Start subroutine --------------
    
    ! store simulation data, T, Tret
    open(594,file=filename,status='unknown',form='formatted')
    write(594,*) time, tref, dt
    close(594)
  
  end subroutine Write_store_file
  !
  ! Read store.info file to restart simulation -> ALL read the same store file
  !
  subroutine Read_store_file(time,tref,dt)
    
    use decomp_2d_mpi, only: mytype,nrank, decomp_2d_abort
    
    real(mytype), intent(OUT) :: time,tref,dt
    integer :: fh
    logical :: file_exists
    
    
    inquire (file='store.info', exist=file_exists)
    
    if (.not. file_exists) then
      call decomp_2d_abort(1, "ERROR: store.info file NOT found in folder!")
    end if
    
    fh=437+nrank
    
    open(fh,file='store.info',status='old',form='formatted')
    read(fh,*) time,tref,dt
    close(fh)
  
  end subroutine Read_store_file
  !
  ! Write 7 - Arrays MAX to an ascii file
  !
  subroutine Write_ascii(filename,arr_1,arr_2,arr_3,arr_4,arr_5,arr_6,arr_7)
    
    use decomp_2d_mpi, only: mytype, decomp_2d_abort
    
    real(mytype), dimension(:), intent(IN) :: arr_1
    real(mytype), dimension(:), intent(IN), optional :: arr_2,arr_3,arr_4,arr_5,arr_6,arr_7
    character(len=*), intent(IN) :: filename
    
    integer :: fh, stat, i
    logical :: file_exist
    ! ------- Start subroutine -----------------
    
    fh=623
    
    ! Check if the file exists
    inquire(file=filename, exist=file_exist)
    
    ! Open the file for appending if it exists, otherwise create a new file
    if (file_exist) then
      ! Open the file in append mode
      open(newunit=fh, file=filename, status="old", action="write", position="append", iostat=stat)
      if (stat /= 0) then
        call decomp_2d_abort(1, "ERROR opening file: "//filename//" for appending.")
      end if
      print *, "Appending to file: ", filename
    else
      ! Create a new file
      open(newunit=fh, file=filename, status="new", action="write", iostat=stat)
      if (stat /= 0) then
        call decomp_2d_abort(1, "ERROR creating file: "//filename)
      end if
      write(*,*) "                 "
      print *, "---> Creating new file: ", filename
      write(*,*) "                 "
    end if
    
    ! Write data to the file
    do i=1,size(arr_1)
      
      if(present(arr_7)) then
        write(fh, *) arr_1(i),arr_2(i),arr_3(i),arr_4(i),arr_5(i),arr_6(i),arr_7(i)
      elseif(present(arr_6)) then
        write(fh, *) arr_1(i),arr_2(i),arr_3(i),arr_4(i),arr_5(i),arr_6(i)
      elseif(present(arr_5)) then
        write(fh, *) arr_1(i),arr_2(i),arr_3(i),arr_4(i),arr_5(i)
      elseif(present(arr_4)) then
        write(fh, *) arr_1(i),arr_2(i),arr_3(i),arr_4(i)
      elseif(present(arr_3)) then
        write(fh, *) arr_1(i),arr_2(i),arr_3(i)
      elseif(present(arr_2)) then
        write(fh, *) arr_1(i),arr_2(i)
      else
        write(fh, *) arr_1(i)
      end if
    
    end do
    
    ! Close the file
    close(fh)
  
  end subroutine Write_ascii
  
  !
  ! Write messages to screen after each iteration
  subroutine Write_screen_forcing(u,v,w,fs_x,fs_y,fs_z,sp)
    
    use decomp_2d_mpi, only: mytype,nrank
    use decomp_2d, only: decomp_info
    use m_aux_spect, only: Spherical_mult
    
    complex(mytype), dimension(:,:,:), intent(IN) :: u,v,w,fs_x,fs_y,fs_z
    real(mytype) :: fx_sq,fy_sq,fz_sq,fx_u,fy_v,fz_w,u_sq,v_sq,w_sq
    
    type(decomp_info), pointer :: sp
    !--------- Start function ----------------
    
    fx_sq = Spherical_mult(fs_x,fs_x,sp)
    fy_sq = Spherical_mult(fs_y,fs_y,sp)
    fz_sq = Spherical_mult(fs_z,fs_z,sp)
    
    if(nrank==0) then
      write(*,9) '<fx*fx>','<fy*fy>','<fz*fz>'
      write(*,11)   fx_sq   ,  fy_sq   ,  fz_sq
    end if
    
    fx_u = Spherical_mult(fs_x,u,sp)
    fy_v = Spherical_mult(fs_y,v,sp)
    fz_w = Spherical_mult(fs_z,w,sp)
    
    u_sq = Spherical_mult(u,u,sp)
    v_sq = Spherical_mult(v,v,sp)
    w_sq = Spherical_mult(w,w,sp)
    
    fx_u = fx_u/sqrt(fx_sq*u_sq)
    fy_v = fy_v/sqrt(fy_sq*v_sq)
    fz_w = fz_w/sqrt(fz_sq*w_sq)
    
    if(nrank==0) then
      write(*,10) '<fx*u>','<fy*v>','<fz*w>'
      write(*,11)   fx_u  , fy_v   , fz_w
      write(*,10) '<u*u>','<v*v>','<w*w>'
      write(*,11)   u_sq  , v_sq , w_sq
      write(*,*) '======================================================='
    end if
    
    ! ==========================================================
    9 FORMAT (4(2x,A12))
    10 FORMAT (/,4(2x,A12))
    11 FORMAT (4(2x,1PE12.4))

  end subroutine Write_screen_forcing

end module m_io
