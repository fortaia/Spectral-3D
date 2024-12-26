!
! MODULE purpose: Set the initial conditions
!
module m_initial_cond
  
  use MPI
  
  use m_glob_params
  
  use decomp_2d_constants
  use decomp_2d_io
  use decomp_2d_mpi
  use decomp_2d_fft
  use decomp_2d
  
  implicit none
  
  private
  
  public :: Initial_cond
  
  !*** VELOCITY FIELDS ***
  ! In physical space
  real(mytype), allocatable, dimension (:,:,:), public :: ur,vr,wr
  ! In Spectral space
  complex(mytype), allocatable, dimension (:,:,:), public :: u,v,w

contains
  !
  ! FIELDS INITIALIZATION - MAIN PROCEDURE
  !
  subroutine Initial_cond(ph,sp)
    
    use m_io
    
    type(decomp_info), pointer :: ph,sp
    real(mytype), parameter :: TWOPI=6.28318530717958647692528676655900
    ! ------------ Start subroutine ------------------
    
    ! input is X-pencil data
    call alloc_x(ur, ph, .true.)
    call alloc_x(vr, ph, .true.)
    call alloc_x(wr, ph, .true.)
    
    ! output is Z-pencil data
    call alloc_z(u, sp, .true.)
    call alloc_z(v, sp, .true.)
    call alloc_z(w, sp, .true.)
    
    select case (JINIT)
      ! New simulation - Initialize fields
    case (0)
      
      if (nrank==0) then
        print *, "                       "
        print *, "** Fields Initialization **"
      endif
      
      if(ISIMU==0) then
        
        ! Initializing the field as a random noise with imposed spectrum
        call Noise(ph,sp)
        
      else
        !Calculating noise
        call Noise(ph,sp)
        
        !Mean field convolution (physical space)
        call IC_mean_field(ph)
        
        !Transform into spectral space
        call decomp_2d_fft_3d(ur,u)
        call decomp_2d_fft_3d(vr,v)
        call decomp_2d_fft_3d(wr,w)
        !DFT Normalization
        u=u/real(N1G*N2G*N3G, mytype)
        v=v/real(N1G*N2G*N3G, mytype)
        w=w/real(N1G*N2G*N3G, mytype)
      
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Restaring from an existing simulation
    case(1)
      
      if (nrank==0) then
        print *, " "
        print *, "      RESTART          "
        print *, "--------------------------------------"
      endif
      
      ! Reading Starting fields
      call Read_3d_field(u,trim(RESDAT)//".ux",sp)
      call Read_3d_field(v,trim(RESDAT)//".uy",sp)
      call Read_3d_field(w,trim(RESDAT)//".uz",sp)
      
      ! Error case handling
    case default
      
      call decomp_2d_abort(__FILE__,__LINE__,JINIT,&
              "ERROR - JINIT value not allowed")
    end select
  
  end subroutine Initial_cond
  
  ! Generation of the mean profile
  subroutine IC_mean_field(ph)
    
    implicit none
    
    type(decomp_info), pointer :: ph
    integer :: j,k
    real(mytype) :: y
    ! *********************************************************************
    !   TOTAL VELOCITY = MEAN + FLUCTUATIONS
    ! *********************************************************************
    
    !---------------------------------------------------
    !   Plane jet profile : hyperbolic tangent formula |
    !---------------------------------------------------
    
    do k=ph%xst(3),ph%xen(3)
      do j=ph%xst(2),ph%xen(2)
        
        ! y is the global coordinate
        y=real(bly,mytype)*(real(j,mytype)-1._mytype)/real(N2G,mytype)-real(bly,mytype)*0.5_mytype
        
        !Stanley and Sarkar:
        if (y<0._mytype) then
          ur(:,j,k)=0.5_mytype*( (U2+U1)+(U2-U1)*tanh(0.5_mytype*(y+0.5_mytype)*HUTHETA) ) + ur(:,j,k)
        elseif (y>=0._mytype) then
          ur(:,j,k)=0.5_mytype*( (U2+U1)+(U2-U1)*tanh(0.5_mytype*(0.5_mytype-y)*HUTHETA) ) + ur(:,j,k)
        endif
      
      end do
    enddo
    
    vr=0._mytype
    wr=0._mytype
  
  
  end subroutine IC_mean_field
  
  ! Random noise to convolute with the mean profile
  subroutine Noise(ph,sp)
    
    use m_utils, only: Random_number_gen_init
    use m_aux_phys, only: Variance_ph
    use m_aux_spect, only : Wavenum
    
    integer :: i,j,k
    real(mytype), parameter :: TWOPI=6.28318530717958647692528676655900
    real(mytype) :: s1,s2,wrk1,wrk2,wrk3,wrk1b,wrk2b,wrk3b,alpha_rnd,rx,ry,rz
    real(mytype) :: uprime2,vprime2,wprime2,specvel,k_norm,k_proj
    real(mytype), allocatable, dimension (:) :: wavenumG1,wavenumG2,wavenumG3
    type(decomp_info), pointer :: ph,sp
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  THE FOURIER VELOCITY COMPONENTS AT K     ->
    !  ARE GENERATED IN THE PLANE ORTHOGONAL TO K
    !      -> ->
    !     (e1,e2) WHERE:
    !                   ->      (Ky , -Kx , 0)
    !                   e1 = -------------------
    !                        SQRT[Ky*Ky + Kx*Kx]
    !                        ->  ->
    !                   ->   K ^ e1
    !                   e2 = -------
    !                          K
    !
    !  [x]   [ sin(theta)cos(phi) cos(theta)cos(phi) -sin(phi) ] [   r   ]
    !  |y| = | sin(theta)sin(phi) cos(theta)sin(phi)  cos(phi) | | theta |
    !  [z]   [ cos(theta)             -sin(theta)       0      ] [  phi  ]
    !
    !          ->
    !          ^ ->  ->  ->     ->  ->
    !          U(K)= Ure(K) + i Uim(K)
    !
    !                ->  ->              ->            ->
    !                Ure(K) = { cos(Phi) e1 + sin(Phi) e2} *
    !                         { SQRT[E(K,0)] | K }
    !
    !                Uim(K) =  Ure(K) .
    !
    !   AND     E(K,0)  :  INITIAL KINETIC-ENERGY SPECTRUM .
    !               ->
    !           Phi(K)   :  RANDOM  NUMBER IN [-PI, PI].
    ! ========================================================
    ! Used constants
    
    s1 = (NSLOPE1-2.0_mytype) / 2.0_mytype
    if (ISIMU==0)then
      s2 = -NSLOPE1/(4.0_mytype*C1*C1)
    else
      s2 = -NSLOPE1/(4.0_mytype*C1*C1)
    end if
    
    
    call Wavenum(wavenumG1,wavenumG2,wavenumG3)
    call Random_number_gen_init()
    
    do k=sp%zst(3),sp%zen(3)
      wrk3=wavenumG3(k)
      wrk3b=wavenumG3(k)*wavenumG3(k)
      do j=sp%zst(2),sp%zen(2)
        wrk2=wavenumG2(j)
        wrk2b=wavenumG2(j)*wavenumG2(j)
        do i=sp%zst(1),sp%zen(1)
          wrk1=wavenumG1(i)
          wrk1b=wavenumG1(i)*wavenumG1(i)
          !                                avoid zero-division
          k_norm=sqrt(wrk1b+wrk2b+wrk3b) + 1.e-10
          k_proj=sqrt(wrk1b+wrk2b)       + 1.e-10
          
          ! Generate random number
          call random_number(alpha_rnd)
          alpha_rnd=TWOPI*alpha_rnd
          
          !             ->                 ->
          !            |k|**s * e**(-s/2 *(|k|/k_peak)**2)
          specvel = k_norm**s1 * exp( s2 * k_norm*k_norm )
          
          !                         cos(theta)cos(phi)          -               sin(phi)
          rx =(sin(alpha_rnd) * wrk1/k_proj * wrk3/k_norm  - cos(alpha_rnd) * wrk2/k_proj ) * specvel
          !                       cos(theta)sin(phi)            +               cos(phi)
          ry =(sin(alpha_rnd) * wrk2/k_proj * wrk3/k_norm  + cos(alpha_rnd) * wrk1/k_proj ) * specvel
          !                          -sin(theta)                +   0
          rz =(sin(alpha_rnd) *   (- k_proj/k_norm) ) * specvel
          
          u(i,j,k) = cmplx(rx,rx,kind=mytype)
          v(i,j,k) = cmplx(ry,ry,kind=mytype)
          w(i,j,k) = cmplx(rz,rz,kind=mytype)
        
        enddo
      enddo
    enddo
    
    !  =============================================================
    !             SPECTRAL  VELOCITY  NORMALISATION
    !
    !   FOR EACH DIRECTION i,
    !   ui <---- ui * FACUVW, WHERE
    !
    !     FACUVW = 1 / SQRT[(< ux*ux >+< uy*uy >+< uz*uz >) / 3TOTE]
    !
    !   SO ,  < ui*ui > = TOTE .
    !   AND  Ec(t=0) = 1/2 * < ux*ux+uy*uy+uz*uz >
    !                = 3/2 * TOTE
    !                = 1/2 Vchar**2
    !
    !   Velocity scale used for normalisation
    !        Vnorm = 1. / normaliz
    !
    !  =============================================================
    
    call decomp_2d_fft_3d(u,ur)
    call decomp_2d_fft_3d(v,vr)
    call decomp_2d_fft_3d(w,wr)
    
    uprime2=Variance_ph(ur)
    vprime2=Variance_ph(vr)
    wprime2=Variance_ph(wr)
    
    wrk1 = sqrt(3.0_mytype*TOTE / (uprime2+vprime2+wprime2))
    
    if(ISIMU==0)then
      u=u*wrk1
      v=v*wrk1
      w=w*wrk1
      return
    end if
    ur=wrk1*ur
    vr=wrk1*vr
    wr=wrk1*wr
    
    ! variance(const*u)=const**2*variance(u)
    wrk1= wrk1**2
    uprime2=uprime2*wrk1
    vprime2=vprime2*wrk1
    wprime2=wprime2*wrk1
    
    call noise_weight(ur,ph)
    call noise_weight(vr,ph)
    call noise_weight(wr,ph)
    
    !--maximum amplitude normalization
    CALL MPI_ALLREDUCE(maxval(ur),wrk1b,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(maxval(vr),wrk2b,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(maxval(wr),wrk3b,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    
    wrk1=TOTE/wrk1b
    wrk2=TOTE/wrk2b
    wrk3=TOTE/wrk3b
    !---Normalsing with prescribed noise amplitude
    ur=wrk1*ur
    vr=wrk1*vr
    wr=wrk1*wr
  
  end subroutine Noise
  ! Random noise localization in physical space
  subroutine Noise_weight(field,ph)
    
    integer :: j,k
    real(mytype) :: y
    real(mytype), dimension(:,:,:) :: field
    type(decomp_info), pointer :: ph
    
    do k=ph%xst(3),ph%xen(3)
      do j=ph%xst(2),ph%xen(2)
        
        ! y is the global coordinate
        y=real(BLY,mytype)*(real(j,mytype)-1._mytype)/real(N2G,mytype)-real(BLY,mytype)*0.5_mytype
        y=abs(y)
        ! Stanley and Sarkar:
        if (y<0.3_mytype) then
          field(:,j,k)=0.75_mytype*field(:,j,k)
        elseif (y>0.7_mytype) then
          field(:,j,k)=0.0_mytype
          !        elseif (y>=0.3_mytype.and.y<=0.7_mytype) then
          !          field(:,j,k)=1.0_mytype*field(:,j,k)
          !          (commented to emphasize the meaning but computationally redundant)
        endif
      
      enddo
    enddo
  
  end subroutine Noise_weight

end module m_initial_cond