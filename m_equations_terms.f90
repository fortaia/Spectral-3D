module m_equations_terms
  
  use decomp_2d_constants, only: mytype
  use decomp_2d, only: decomp_info
  
  implicit none
  
  private
  
  public :: Continuity, Dissipative, NonLin, Forcing, RHS
  
contains
  
  !
  ! Continuity equation
  !
  subroutine Continuity(u,v,w,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
    ! ================================================
    ! THIS ROUTINE GENERATES A DIVERGENCE-FREE SPECTRAL
    ! VELOCITY FIELD, FROM THE INPUT FIELD.
    !
    ! INPUT: (u,v,w) SPECTRAL VELOCITY COMPONENTS
    !
    !                                | ->      -> |
    !       ->                 ->    | k.FFTD{ u }|  ->
    ! FFTD{ u }  <-----  FFTD{ u } - | -----------|  k .
    !                                |     k**2   |
    !
    !           ->   ->
    ! FFTD{ u }( 0 ) IS NOT MODIFIED HERE, BUT
    ! IN ROUTINE trunc.
    ! ================================================
    type(decomp_info), pointer :: sp
    
    integer :: i,j,k,start
    real(mytype) :: wave_sq,wv1,wv2,wv3,sq_w2,sq_w3
    real(mytype), dimension(:), intent(IN) :: wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3
    complex(mytype), dimension(sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(INOUT) :: u,v,w
    !!!!!!!!!!!!!!!!!!! Variable Declaration END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    start=sp%zst(1)
    
    if(sp%zst(1)==1 .and. sp%zst(2)==1 .and. sp%zst(3)==1) start=2
    
    do k=sp%zst(3),sp%zen(3)
      sq_w3=sq_wnumG3(k)
      wv3=wavenumG3(k)
      do j=sp%zst(2),sp%zen(2)
        sq_w2=sq_wnumG2(j)
        wv2=wavenumG2(j)
        do i=start,sp%zen(1)
          wv1=wavenumG1(i)
          ! k**2
          wave_sq = (sq_wnumG1(i)+ sq_w2 + sq_w3)
          ! 1/k**2
          wave_sq=1.0_mytype/wave_sq
          
          ! Correcting the velocity field: es. u - k_1/|k|*(k*u)
          u(i,j,k) = u(i,j,k) - wv1*(wv1*u(i,j,k) +wv2*v(i,j,k)+ wv3*w(i,j,k))*wave_sq
          v(i,j,k) = v(i,j,k) - wv2*(wv1*u(i,j,k) +wv2*v(i,j,k)+ wv3*w(i,j,k))*wave_sq
          w(i,j,k) = w(i,j,k) - wv3*(wv1*u(i,j,k) +wv2*v(i,j,k)+ wv3*w(i,j,k))*wave_sq
        
        enddo
        start=sp%zst(1)
      enddo
    enddo
    
    return
  
  end subroutine Continuity
  
  !
  ! Dissipative TERM
  !
  subroutine Dissipative(u,v,w,sq_wnumG1,sq_wnumG2,sq_wnumG3,rhs_x,rhs_y,rhs_z,sp)
    !
    ! *****************************************
    !
    ! ====================================================================
    ! COMPUTATION OF THE VISCOUS TERMS OF THE NAVIER-STOKES EQUATIONS
    !
    !
    !   NAVIER-STOKES:
    !
    !   (Physical space)               (Spectral space)
    !
    !   d        d                       2 ^
    !   -  [  nu - (u_i)]  -> - nu |Kappa| u_i
    !   dx       dx
    !     j        j
    !
    !              2  2  2  2
    ! WHERE: |Kappa|=k1+k2+k3
    !
    ! ====================================================================
    use m_glob_params, only: RNU
    
    type(decomp_info), pointer :: sp
    
    integer :: i,j,k
    real(mytype) :: visc_ksq,dummy2,dummy3
    real(mytype), dimension (:), intent(IN) :: sq_wnumG1,sq_wnumG2,sq_wnumG3
    complex(mytype), dimension (sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(IN) :: u,v,w
    complex(mytype), dimension (sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(INOUT) :: rhs_x,rhs_y,rhs_z
    
    !------------ Start subroutine ----------------------------
    
    do k=sp%zst(3),sp%zen(3)
      dummy3=sq_wnumG3(k)
      do j=sp%zst(2),sp%zen(2)
        dummy2=sq_wnumG2(j)
        do i=sp%zst(1),sp%zen(1)
          
          visc_ksq= -RNU *(sq_wnumG1(i) + dummy2 + dummy3)
          
          rhs_x(i,j,k) =  visc_ksq*u(i,j,k)
          rhs_y(i,j,k) =  visc_ksq*v(i,j,k)
          rhs_z(i,j,k) =  visc_ksq*w(i,j,k)
        
        enddo
      enddo
    enddo
  
  end subroutine Dissipative
  
  !
  ! Non-linear (convective) TERM
  !
  subroutine NonLin(u,v,w,ur,vr,wr,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph,&
                    wavenumG1,wavenumG2,wavenumG3,normaliz,trunc_index,ph,sp)
    
    use m_aux_spect, only: Calc_vorticity_sp
    use m_utils, only: Dealiasing_trunc
    
    use decomp_2d_fft
    
    type(decomp_info), pointer :: ph,sp
    
    integer :: i,j,k
    integer, allocatable, dimension(:,:), intent(IN) :: trunc_index
    real(mytype) :: dummy1
    real(mytype), intent(IN) :: normaliz
    real(mytype), dimension (:), intent(IN) :: wavenumG1,wavenumG2,wavenumG3
    ! In Physical space: Working arrays
    real(mytype), dimension(ph%xst(1):,ph%xst(2):,ph%xst(3):), intent(INOUT) :: ur,vr,wr
    real(mytype), dimension(ph%xst(1):,ph%xst(2):,ph%xst(3):) :: work_xph,work_yph,work_zph
    ! In Spectral space
    complex(mytype), dimension(sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(INOUT) :: u,v,w
    complex(mytype), dimension(sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(INOUT) :: work_xsp1,work_ysp1,work_zsp1
    
    !  ======================================================
    !  USING PSEUDO-SPECTRAL METHODS, THIS ROUTINE COMPUTES
    !  THE NONLINEAR TERMS OF THE NAVIER-STOKES Eqs.
    !
    !                 -->  ->            ->    -->      ->
    !  OUTPUT :       rhs( k ) =  FFTD{  u  ^  rot(u)}( k )
    !
    !                       ->               ->
    !                 WHERE u  = FFTI[ FFTD{ u }]
    !
    !                  -->                   -->
    !              &   rot(u)  = FFTI[ FFTD{rot(u)} ]
    !
    !            ^ ^ ^
    !  INPUT :  (u,v,w)   SPECTRAL VELOCITY COMPONENTS
    !
    !           (ur,vr,wr) WILL CONTAIN THE VELOCITY
    !                      COMPONENTS IN THE PHYSICAL SPACE.
    !
    !   ################################################
    !
    !        -->         ->         ->
    !  FFTD{ rot(u)} = i k  ^ FFTD{ u }
    !
    !                         ^     ^
    !                    | k1*w - k3*v
    !        -->         |
    !                    |    ^      ^
    !  FFTD{ rot(u)} = i | k3*u - k1*w
    !                    |
    !                    |    ^      ^
    !                    | k1*v - k2*u
    !                    |
    !
    !  ======================================================
    ! Velocity in physical space
    
    work_xsp1=u
    call decomp_2d_fft_3d(work_xsp1,ur)
    work_xsp1=v
    call decomp_2d_fft_3d(work_xsp1,vr)
    work_xsp1=w
    call decomp_2d_fft_3d(work_xsp1,wr)
    !
    !  FIRST STEP : calculate curl of u
    !                    ->          ->
    !  (tnx,tny,tnz) = i k  ^  FFTD{ u }
    
    call Calc_vorticity_sp(u,v,w,wavenumG1,wavenumG2,wavenumG3,work_xsp1,work_ysp1,work_zsp1,sp)
    
    !
    ! SECOND STEP: Put all in the physical space
    !                          -->                   ->
    !   (rhs_x,rhs_y,rhs_z) =  rot(u) ; (ur,vr,wr) = u
    
    call decomp_2d_fft_3d(work_xsp1,work_xph)
    call decomp_2d_fft_3d(work_ysp1,work_yph)
    call decomp_2d_fft_3d(work_zsp1,work_zph)
    
    ! THIRD STEP: do the cross product
    !                 ->   -->
    ! (tnx,tny,tnz) = u  ^ rot(u)
    
    do k=ph%xst(3),ph%xen(3)
      do j=ph%xst(2),ph%xen(2)
        do i=ph%xst(1),ph%xen(1)
          
          dummy1          = vr(i,j,k)*work_zph(i,j,k) - wr(i,j,k)*work_yph(i,j,k)
          work_zph(i,j,k) = wr(i,j,k)*work_xph(i,j,k) - ur(i,j,k)*work_zph(i,j,k)
          work_yph(i,j,k) = ur(i,j,k)*work_yph(i,j,k) - vr(i,j,k)*work_xph(i,j,k)
          work_xph(i,j,k) = dummy1
          dummy1          = work_yph(i,j,k)
          work_yph(i,j,k) = work_zph(i,j,k)
          work_zph(i,j,k) = dummy1
        
        enddo
      enddo
    enddo
    
    
    !FORTH STEP: Back to Spectral space
    !                       ->   -->
    ! (tnx,tny,tnz) = FFTD {u  ^ rot(u)}
    call decomp_2d_fft_3d(work_xph,work_xsp1)
    call decomp_2d_fft_3d(work_yph,work_ysp1)
    call decomp_2d_fft_3d(work_zph,work_zsp1)
    
    call Dealiasing_trunc(work_xsp1, work_ysp1, work_zsp1, trunc_index,sp)
    
    work_xsp1 = work_xsp1/normaliz
    work_ysp1 = work_ysp1/normaliz
    work_zsp1 = work_zsp1/normaliz
  
  end subroutine NonLin
  !
  ! Forcing HIT TERM
  !
  subroutine Forcing(u,v,w,fs_x,fs_y,fs_z,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
    
    ! =================================================
    ! Sepand Ossia -  February 2000
    ! Ref: "RANDOM FORCING OF 3D HOMOGENEOUS TURBULENCE"
    ! Phys. Fluids , Vol 11(7) , 1999, pp 1880-1889.
    !
    ! The forcing scheme provides a divergence-free force
    ! which is uncorrelated in time with the velocity.
    ! =================================================
    use m_glob_params
    use m_utils, only: Random_number_gen
    use m_aux_spect, only: Spherical_mult
    
    use decomp_2d_mpi, only: mytype,nrank
    
    integer :: i,j,k,skip
    real(mytype), parameter :: TWOPI=6.28318530717958647692528676655900
    real(mytype) :: w1,w2,w3,w1_s,w2_s,w3_s,normaliz_fs,numerator,denominator,mod_const,k_norm,k_proj
    real(mytype) :: e2(3),e1(2),psi_rnd,two_phi_rnd,theta1,theta2
    real(mytype), dimension (:), intent(IN) :: wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3
    complex(mytype) :: xsi1,xsi2,a_rnd,b_rnd
    complex(mytype), dimension(:,:,:), intent(IN) :: u,v,w
    complex(mytype), dimension(:,:,:), intent(INOUT) :: fs_x,fs_y,fs_z
    type(decomp_info), pointer :: sp
    
    ! ------------ Start subroutine ------------------
    mod_const=-1.0_mytype
    skip=0
    if(nrank==0) skip=1
    
    do k=sp%zst(3),sp%zen(3)
      w3=wavenumG3(k)
      w3_s=sq_wnumG3(k)
      do j=sp%zst(2),sp%zen(2)
        w2=wavenumG2(j)
        w2_s=sq_wnumG2(j)
        do i=sp%zst(1)+skip,sp%zen(1)
          w1=wavenumG1(i)
          w1_s=sq_wnumG1(i)
          !                                avoid zero-division
          k_norm=sqrt(w1_s+w2_s+w3_s) + 1.e-10
          if (k_norm > X0) cycle
          
          k_proj=sqrt(w1_s+w2_s)      + 1.e-10
          
          e2(1)=w1*w3/(k_norm*k_proj)
          e2(2)=w2*w3/(k_norm*k_proj)
          e2(3)=-(k_proj/k_norm)
          
          e1(1)= w2/k_proj
          e1(2)=-w1/k_proj
          
          mod_const = exp(-0.5_mytype / WIDTH_F * (k_norm-X0)*(k_norm-X0))/k_norm
          
          psi_rnd = TWOPI*Random_number_gen()
          two_phi_rnd = TWOPI*Random_number_gen()
          
          xsi1=cmplx(real(u(i,j,k))*e1(1)+real(v(i,j,k))*e1(2),&
                  imag(u(i,j,k))*e1(1)+imag(v(i,j,k))*e1(2), mytype)
          
          xsi2=cmplx(real(u(i,j,k))*e2(1)+real(v(i,j,k))*e2(2)+real(w(i,j,k))*e2(3),&
                  imag(u(i,j,k))*e2(1)+imag(v(i,j,k))*e2(2)+imag(w(i,j,k))*e2(3), mytype)
          
          numerator= sin(two_phi_rnd)*real(xsi1)+ cos(two_phi_rnd)*(sin(psi_rnd)*imag(xsi2) + cos(psi_rnd)*real(xsi2))
          denominator= -sin(two_phi_rnd)*imag(xsi1)+ cos(two_phi_rnd)*(sin(psi_rnd)*real(xsi2) -cos(psi_rnd)*imag(xsi2))
          
          if(numerator==0._mytype.and.denominator==numerator) cycle
          
          theta1=atan(numerator/(denominator+1.e-20))
          theta2=psi_rnd+theta1
          
          a_rnd = cmplx(cos(theta1),sin(theta1), mytype)*mod_const*sin(two_phi_rnd)
          b_rnd = cmplx(cos(theta2),sin(theta2), mytype)*mod_const*cos(two_phi_rnd)
          
          fs_x(i,j,k) = cmplx(real(a_rnd)*e1(1)+real(b_rnd)*e2(1), imag(a_rnd)*e1(1)+imag(b_rnd)*e2(1), mytype)
          fs_y(i,j,k) = cmplx(real(a_rnd)*e1(2)+real(b_rnd)*e2(2), imag(a_rnd)*e1(2)+imag(b_rnd)*e2(2), mytype)
          fs_z(i,j,k) = cmplx(real(b_rnd)*e2(3)                  , imag(b_rnd)*e2(3), mytype)
        
        enddo
        skip=0
      enddo
    enddo
    
    normaliz_fs = Spherical_mult(fs_x,fs_x,sp)+Spherical_mult(fs_y,fs_y,sp)+Spherical_mult(fs_z,fs_z,sp)
    
    !   Normalization such that
    !     0.5 * < fi*fi > = Totf / Dt
    
    normaliz_fs = sqrt( (2.0_mytype*(TOTF/DT)) / normaliz_fs)
    
    fs_x = fs_x * normaliz_fs
    fs_y = fs_y * normaliz_fs
    fs_z = fs_z * normaliz_fs
  
  end subroutine Forcing
  
  !
  ! Non-linear (convective) + Dissipative TERM
  !
  subroutine RHS(u,v,w,ur,vr,wr,rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph, &
                 wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,normaliz,trunc_index,ph,sp)
    
    use m_glob_params, only: RNU
    use m_utils, only: Dealiasing_trunc
    
    use decomp_2d_fft
    
    type(decomp_info), pointer :: ph,sp
    
    integer :: i,j,k
    integer, allocatable, dimension(:,:), intent(IN) :: trunc_index
    real(mytype) :: dummy1,dummy2,dummy3,Ru,Rv,Rw,Iu,Iv,Iw,sq_dummy2,sq_dummy3,visc_ksq
    real(mytype), intent(IN) :: normaliz
    real(mytype), dimension (:), intent(IN) :: wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3
    ! In Physical space: Working arrays
    real(mytype), dimension (ph%xst(1):,ph%xst(2):,ph%xst(3):), intent(INOUT) :: ur,vr,wr
    real(mytype), dimension (ph%xst(1):,ph%xst(2):,ph%xst(3):), intent(INOUT) :: work_xph,work_yph,work_zph
    ! In Spectral space
    complex(mytype), dimension(sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(INOUT) :: u,v,w
    complex(mytype), dimension (sp%zst(1):,sp%zst(2):,sp%zst(3):) :: rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1
    
    !  ======================================================
    !  USING PSEUDO-SPECTRAL METHODS, THIS ROUTINE COMPUTES
    !  THE RIGHT HAND SIDE TERMS OF THE NAVIER-STOKES Eqs.
    !
    !  FIRST STEP : calculate curl of u +  DISSIPATIVE TERM
    !                    ->          ->
    !  (tnx,tny,tnz) = i k  ^  FFTD{ u }
    
    do k=sp%zst(3),sp%zen(3)
      dummy3=wavenumG3(k)
      sq_dummy3=sq_wnumG3(k)
      do j=sp%zst(2),sp%zen(2)
        dummy2=wavenumG2(j)
        sq_dummy2=sq_wnumG2(j)
        do i=sp%zst(1),sp%zen(1)
          dummy1=wavenumG1(i)
          Ru=REAL(u(i,j,k));  Rv=REAL(v(i,j,k));  Rw=REAL(w(i,j,k))
          Iu=AIMAG(u(i,j,k)); Iv=AIMAG(v(i,j,k)); Iw=AIMAG(w(i,j,k))
          
          ! Calculate real and imaginary parts directly to avoid complex multiplication overhead
          work_xsp1(i,j,k) = cmplx(-Iw * dummy2 + Iv * dummy3, &
                  Rw * dummy2 - Rv * dummy3, mytype)
          
          work_ysp1(i,j,k) = cmplx(-Iu * dummy3 + Iw * dummy1, &
                  Ru * dummy3 - Rw * dummy1, mytype)
          
          work_zsp1(i,j,k) = cmplx(-Iv * dummy1 + Iu * dummy2, &
                  Rv * dummy1 - Ru * dummy2, mytype)
          
          visc_ksq = -RNU * (sq_wnumG1(i)+ sq_dummy2 + sq_dummy3)
          
          rhs_x(i,j,k) =  visc_ksq*u(i,j,k)
          rhs_y(i,j,k) =  visc_ksq*v(i,j,k)
          rhs_z(i,j,k) =  visc_ksq*w(i,j,k)
        
        enddo
      enddo
    enddo
    
    !
    ! SECOND STEP: Put all in the physical space
    !                          -->                          ->
    !   (rhs_x,rhs_y,rhs_z) =  rot(u) ; whereas (ur,vr,wr) = u is known
    
    call decomp_2d_fft_3d(work_xsp1,work_xph)
    call decomp_2d_fft_3d(work_ysp1,work_yph)
    call decomp_2d_fft_3d(work_zsp1,work_zph)
    
    ! THIRD STEP: do the cross product
    !                 ->   -->
    ! (tnx,tny,tnz) = u  ^ rot(u)
    
    do k=ph%xst(3),ph%xen(3)
      do j=ph%xst(2),ph%xen(2)
        do i=ph%xst(1),ph%xen(1)
          
          dummy1          = vr(i,j,k)*work_zph(i,j,k) - wr(i,j,k)*work_yph(i,j,k)
          work_zph(i,j,k) = wr(i,j,k)*work_xph(i,j,k) - ur(i,j,k)*work_zph(i,j,k)
          work_yph(i,j,k) = ur(i,j,k)*work_yph(i,j,k) - vr(i,j,k)*work_xph(i,j,k)
          work_xph(i,j,k) = dummy1
          dummy1          = work_yph(i,j,k)
          work_yph(i,j,k) = work_zph(i,j,k)
          work_zph(i,j,k) = dummy1
        
        enddo
      enddo
    enddo
    
    !FORTH STEP: Back to Spectral space
    !                       ->   -->
    ! (tnx,tny,tnz) = FFTD {u  ^ rot(u)}
    call decomp_2d_fft_3d(work_xph,work_xsp1)
    call decomp_2d_fft_3d(work_yph,work_ysp1)
    call decomp_2d_fft_3d(work_zph,work_zsp1)
    
    call dealiasing_trunc(work_xsp1, work_ysp1, work_zsp1, trunc_index,sp)
    
    work_xsp1 = work_xsp1/normaliz
    work_ysp1 = work_ysp1/normaliz
    work_zsp1 = work_zsp1/normaliz
    
    rhs_x=rhs_x+work_xsp1
    rhs_y=rhs_y+work_ysp1
    rhs_z=rhs_z+work_zsp1
  
  end subroutine RHS
  
end module m_equations_terms