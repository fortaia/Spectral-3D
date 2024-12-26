!
! MODULE purpose: Group useful procedures for calculations in spectral space
!
module m_aux_spect
  
  implicit none
  
  private
  
  public :: Wavenum, Shell_average, Calc_kin_en_sp, Spherical_mult, Calc_vorticity_sp

contains
  !
  ! Calculates the wavenumbers in each pencil
  !
  subroutine Wavenum(wavenumG1,wavenumG2,wavenumG3)
    
    use decomp_2d_mpi
    use m_glob_params
    
    implicit none
    
    integer :: i,j,k,cnt
    real(mytype), parameter :: TWOPI=6.28318530717958647692528676655900
    real(mytype) :: FX,FY,FZ,HX,HY,HZ
    
    real(mytype), allocatable, dimension (:), intent(OUT) :: wavenumG1,wavenumG2,wavenumG3
    
    ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    allocate(wavenumG1(N1G/2+1),wavenumG2(N2G),wavenumG3(N3G))
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (NBOUND==0) then
      
      if (ISIMU==0.OR.ISIMU==2) then
        ! HOMOGENEOUS ISOTROPIC TURBULENCE
        ! WITH BOX WITH EQUAL SIZES (A0=1.)
        
        HX = TWOPI ! Box side length in x
        HY = TWOPI ! Box side length in y
        HZ = TWOPI ! Box side length in z
        
        FX = 1._mytype ! TWOPI/HX
        FY = 1._mytype
        FZ = 1._mytype
      
      elseif (ISIMU==1.OR.ISIMU==3.OR.ISIMU==4) then
        ! TEMPORAL TURBULENT PLANE JET OR TEMPORAL TURBULENT PLANE WAKE
        ! WITH ASPECT RATIOS A0X=LX/LY AND A0Z=LZ/LY
        
        HX = TWOPI * A0X ! Box side length in x
        HY = TWOPI       ! Box side length in y
        HZ = TWOPI * A0Z ! Box side length in z
        
        FX = A0X        !TWOPI/HX
        FY = 1._mytype  !TWOPI/HY
        FZ = A0Z        !TWOPI/HZ
      
      endif
    
    else
      
      HX = TWOPI * 0.5_mytype/ A0X
      HY = TWOPI * 0.5_mytype
      HZ = TWOPI * 0.5_mytype/ A0Z
      
      FX = TWOPI * 0.5_mytype/ HX
      FY = TWOPI * 0.5_mytype/ HY
      FZ = TWOPI * 0.5_mytype/ HZ
    
    endif
    
    !======
    ! GENERATING X-DIRECTION WAVENUMBERS
    !======
    do i=1,N1G/2+1
      wavenumG1(i)  = real(i-1,mytype) * FX
    enddo
    !======
    ! GENERATING Y-DIRECTION WAVENUMBERS
    !======
    if (N2G==1) then
      wavenumG2(1) = 0._mytype
    else
      do j=1,N2G
        cnt = j-1
        if (j>N2G/2)  cnt = cnt-N2G
        !IF (j==N2G/2+1) jj = 0 not clear why
        wavenumG2(j) = cnt * FY
      enddo
    endif
    !======
    ! GENERATING Z-DIRECTION WAVENUMBERS
    !======
    if (NBOUND==0) then
      do k=1,N3G
        cnt = k-1
        if (k>N3G/2)  cnt = cnt-n3G
!        if (K.EQ.N3hp) KK = 0 not clear why
        wavenumG3(k) = cnt * FZ
      enddo
    else
      do k=1,N3G
        wavenumG3(k) = real(k-1,mytype) * FZ
      enddo
    endif
    ! SETTING THE USELESS WAVENUMBERS TO ZERO. THESE
    ! WAVENUMBERS ARE ONLY GENERATED SO THAT THE FIELDS
    ! FULFILL THE SPECIFIC DIMENSIONS OF THE FFT.
    
!    FAC3G(n3p) = 0D0 NOT CLEAR WHY
    return
    
  end subroutine Wavenum
  !
  ! Calculates shell averages of a field
  !
  subroutine Shell_average(in_var,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp,count,out_var)
    
    use m_glob_params, only: N2G
    use decomp_2d_mpi, only: mytype
    use decomp_2d, only: decomp_info
    
    type(decomp_info), pointer :: sp
    
    integer :: i,j,k,wave
    integer, dimension(N2G/2+1), intent(OUT) :: count
    real(mytype) :: wave_radius,sq_w2,sq_w3,symm_factor
    real(mytype), contiguous, dimension(:), intent(IN) :: sq_wnumG1,sq_wnumG2,sq_wnumG3
    real(mytype), dimension(N2G/2+1), intent(OUT) :: out_var
    complex(mytype), dimension(sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(IN) :: in_var
    
    !---------   Subroutine start  ------------------
    
    count = 0
    out_var = 0.0_mytype
    
    do k=sp%zst(3),sp%zen(3)
      sq_w3=sq_wnumG3(k)
      do j=sp%zst(2),sp%zen(2)
        sq_w2=sq_wnumG2(j)
        do i=sp%zst(1),sp%zen(1)
          
          ! ||k||
          wave_radius = sqrt(sq_wnumG1(i)+ sq_w2 + sq_w3)
          symm_factor = 2.0_mytype
          if(i==1) symm_factor = 1.0_mytype
          ! Round to closest integer (0.5 rounded to 0)
          wave = nint(wave_radius)
          if (wave>N2G/2) cycle
          ! Accumulating quantities in array
          out_var(wave+1) = out_var(wave+1) + real(in_var(i,j,k)) * symm_factor
          count(wave+1) = count(wave+1) + 1 * int(symm_factor)
          
        enddo
      enddo
    enddo
  
  end subroutine Shell_average
  !
  ! Collective square in spectral space whole sphere
  !
  function Spherical_mult(field1,field2,sp) result(res_sph)

    use decomp_2d_constants, only: mytype, real_type
    use decomp_2d
    use MPI
    
    type(decomp_info), pointer :: sp
    
    integer :: ierr
    real(mytype) :: res_loc, res_sph
    complex(mytype), dimension(:,:,:), intent(IN) :: field1,field2
    
    !---------   Subroutine start  ------------------
    
    if (sp%zst(1)==1) then
      
      res_loc = 2.0_mytype*sum( real(field1(2:sp%zen(1),:,:))*real(field2(2:sp%zen(1),:,:)) + &
                                imag(field1(2:sp%zen(1),:,:))*imag(field2(2:sp%zen(1),:,:)) ) + &
                sum( real(field1(1,:,:))*real(field2(1,:,:)) + imag(field1(1,:,:))*imag(field2(1,:,:)) )
    
    else
      
      res_loc = 2.0_mytype*sum( real(field1)*real(field2) + imag(field1)*imag(field2) )
      
    end if
    
    CALL MPI_ALLREDUCE(res_loc,res_sph,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)

  end function Spherical_mult
  !
  ! Calculates kinetic energy in spectral space
  !
  subroutine Calc_kin_en_sp(u,v,w,kin_en)
    
    use decomp_2d_mpi, only: mytype
    
    complex(mytype), contiguous, dimension(:,:,:), intent(IN) :: u,v,w
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: kin_en
    
    !---------   Subroutine start  ------------------
    
    kin_en%im = 0.0_mytype
    kin_en%re = 0.5_mytype * (u%re*u%re + u%im*u%im + &
            v%re*v%re + v%im*v%im + &
            w%re*w%re + w%im*w%im )
  end subroutine Calc_kin_en_sp
  !
  ! Calculates vorticity in spectral space -> OUT var: vortx,vorty,vortz
  !
  subroutine Calc_vorticity_sp(u,v,w,wavenumG1,wavenumG2,wavenumG3,vortx,vorty,vortz,sp)
    
    use decomp_2d_constants, only: mytype
    use decomp_2d, only: decomp_info
    
    type(decomp_info), pointer :: sp
    
    integer :: i,j,k
    real(mytype) :: dummy1,dummy2,dummy3,Ru,Rv,Rw,Iu,Iv,Iw
    real(mytype), dimension(:), intent(IN) :: wavenumG1,wavenumG2,wavenumG3
    complex(mytype), dimension(sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(IN) :: u,v,w
    complex(mytype), dimension(sp%zst(1):,sp%zst(2):,sp%zst(3):), intent(INOUT) :: vortx,vorty,vortz
    
    !---------   Subroutine start  ------------------
    !                                ->   ->
    ! Vorticity in spectral space = i k ^ u
    !
    
    do k=sp%zst(3),sp%zen(3)
      dummy3=wavenumG3(k)
      do j=sp%zst(2),sp%zen(2)
        dummy2=wavenumG2(j)
        do i=sp%zst(1),sp%zen(1)
          dummy1=wavenumG1(i)
          
          Ru=u(i,j,k)%RE; Rv=v(i,j,k)%RE; Rw=w(i,j,k)%RE
          Iu=u(i,j,k)%IM; Iv=v(i,j,k)%IM; Iw=w(i,j,k)%IM
          
          ! Calculate real and imaginary parts directly to avoid complex multiplication overhead
          vortx(i,j,k) = cmplx(-Iw * dummy2 + Iv * dummy3, Rw * dummy2 - Rv * dummy3, mytype)
          
          vorty(i,j,k) = cmplx(-Iu * dummy3 + Iw * dummy1, Ru * dummy3 - Rw * dummy1, mytype)
          
          vortz(i,j,k) = cmplx(-Iv * dummy1 + Iu * dummy2, Rv * dummy1 - Ru * dummy2, mytype)
        
        enddo
      enddo
    enddo
    
  end subroutine Calc_vorticity_sp

end module m_aux_spect