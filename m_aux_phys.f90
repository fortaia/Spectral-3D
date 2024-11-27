module m_aux_phys
  
  use m_glob_params
  use decomp_2d_constants
  
  implicit none
  
  private
  
  public :: Variance_ph, Stats_quant_init
  
  ! Generic interface to handle multiple data types
  interface Variance_ph
    module procedure Variance_ph_1field
    module procedure Variance_ph_2fields
  end interface
  
contains
  
  ! Calulate variance var(u)=mean(u**2)-mean(u)**2
  function Variance_ph_1field(field1) result(variance)
    
    use MPI
    
    real(mytype), dimension(:,:,:), intent(IN) :: field1
    real(mytype) :: variance
    real(mytype), dimension(2) :: loc_sum,gl_sum
    
    ! COMPUTE MEAN AND TAKE MEAN OUT
    loc_sum(1) = sum(field1*field1)
    loc_sum(2) = sum(field1)

    call MPI_ALLREDUCE(loc_sum,gl_sum,2,real_type,MPI_SUM, MPI_COMM_WORLD,ierr)
    
    gl_sum(2) = gl_sum(2)/real(N1G*N2G*N3G,mytype)
    gl_sum(1) = gl_sum(1)/real(N1G*N2G*N3G,mytype)
    
    variance = gl_sum(1) - gl_sum(2)*gl_sum(2)
    
    return
    
  end function Variance_ph_1field
  
  ! Calulate variance var(uv)=mean(u*v)-mean(u)*mean(v)
  function Variance_ph_2fields(field1,field2) result(variance)
    
    use MPI
    
    real(mytype), dimension(:,:,:), intent(IN) :: field1,field2
    real(mytype) :: variance
    real(mytype), dimension(3) :: loc_sum,gl_sum
    
    ! COMPUTE MEAN AND TAKE MEAN OUT
    loc_sum(1) = sum(field1)
    loc_sum(2) = sum(field2)
    loc_sum(3) = sum(field1*field2)
    
    call MPI_ALLREDUCE(loc_sum,gl_sum,2,real_type,MPI_SUM, MPI_COMM_WORLD,ierr)
    
    gl_sum(1)=gl_sum(1)/real(N1G*N2G*N3G,mytype)
    gl_sum(2)=gl_sum(2)/real(N1G*N2G*N3G,mytype)
    
    variance=gl_sum(3)/real(N1G*N2G*N3G,mytype) - gl_sum(1)*gl_sum(2)
    
    return
  
  end function Variance_ph_2fields
  
  subroutine Stats_quant_init(ur,vr,wr,Tref)
    
    real(mytype), dimension(:,:,:), intent(IN) :: ur,vr,wr
    real(mytype), parameter :: TWOPI=6.28318530717958647692528676655900
    real(mytype) :: uprime2,vprime2,wprime2,Vchar
    real(mytype), intent(OUT) :: Tref
    
    !  =============================================================
    !   TIME NORMALISATION: Instead of using the initial integral
    !   scale Linteg(t=0) we used the box size:
    !
    !   Tret   = Linteg(t=0) / Vchar
    !
    !   Characteristic velocity  scale used in
    !   literature relative to isotropic turbulence:
    !        Vchar = SQRT(< ux*ux+uy*uy+uz*uz >) = SQRT(3*TOTE)
    !  =============================================================
    
    uprime2=variance_ph(ur)
    vprime2=variance_ph(vr)
    wprime2=variance_ph(wr)
    
    Vchar = sqrt(uprime2 + vprime2 + wprime2)
    Tref  = TWOPI / Vchar
    !    Uturb = Vchar/dSQRT(3D0)
    !    Reynolds = (Uturb * Linteg )/ (RNU+1.E-20)
    !------------------------------------------------------------------------------
    !    ! ==========================================================
    !    IF (Ntsr.Eq.0) THEN
    !
    !      Vchar= (Totf / X0)**(1D0/3D0)
    !
    !      !           LARGE-EDDY TRUNOVER TIME
    !      Tret = 1D0/(X0 * Vchar)
    !
    !      IF (Rnu.NE.0D0) THEN
    !        Re_kf = Totf**(1D0/3D0)   * X0**(-4D0/3D0) / Rnu
    !        Etha  = 1D0/ (Re_kf**(3D0/4D0) * X0)
    !      ELSE
    !        Re_kf = 1.0E20
    !        Etha  = 1.0E-20
    !      ENDIF
    !
    !      !    Minimum viscosity for DNS
    !
    !      Quality   = 1.5D0
    !      Etha_dns  = Quality / Xkc
    !      Re_kf_dns     = 1D0/(X0*Etha_dns)**(4D0/3D0)
    !      !            Rnu_dns       = Totf**(1./3.)*X0**(-4./3.) / Re_kf_dns
    !      Rnu_dns       = Totf**(1D0/3D0) * (Quality / Xkc)**(4D0/3D0)
    !      !            Re_lambda_dns = dsqrt((2./3.)*Ekin) * Lambda/ Rnu_dns
    !
    !      ! DETERMINATION OF max||div(f)|| = max||k.f||
    !
    !      DIV=-1.E+35
    !
    !      DO k=1,nCdim(3)
    !        DO j=1,nCdim(2)
    !          DO i=1,nCdim(1),2
    !            ip=i+1
    !
    !            DIVTEMP = ( fac1l(i)*fs_x(i,j,k)+       &
    !                    fac2l(j)*fs_y(i,j,k)+       &
    !                    fac3l(k)*fs_z(i,j,k)  )**2
    !
    !            DIVTEMP =   DIVTEMP +                 &
    !                    ( fac1l(i)*fs_x(ip,j,k)+         &
    !                            fac2l(j)*fs_y(ip,j,k)+         &
    !                            fac3l(k)*fs_z(ip,j,k) )**2
    !
    !            DIVTEMP = DSQRT(DIVTEMP)
    !
    !            IF (DIVTEMP.GT.DIV) DIV = DIVTEMP
    !
    !          ENDDO
    !        ENDDO
    !      ENDDO
    !
    !#ifdef MPI30
    !  CALL MPI_IALLREDUCE(DIV,DIVTEMP,1,real_type,MPI_MAX, &
    !                      world,handle,ierr)
    !  CALL MPI_WAIT(handle,status, ierr)
    !#else
    !      CALL MPI_ALLREDUCE(DIV,DIVTEMP,1,real_type,MPI_MAX, &
    !              world,ierr)
    !#endif
    !
    !      DIV = DIVTEMP
    !
    !
    !      ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !      if(myid.eq.0)then
    !
    !        WRITE(*,100) 'INITIAL FORCING FIELD'
    !        WRITE(*,110)
    !        WRITE(*,120) 'Truncated Field    ',UXX1,UYY1,UZZ1
    !        !       WRITE(*,120) 'K-normal Projection',UXX2,UYY2,UZZ2
    !        WRITE(*,120) 'Normalised Field   ',FXSQ,FYSQ,FZSQ
    !        WRITE(*,130)  'Totf'  , 'Vchar', 'Tret'
    !        WRITE(*,135)   Totf   ,  Vchar ,  Tret
    !        WRITE(*,130)  'Div(f)', 'Re_kf', 'Rnu_dns'
    !        WRITE(*,135)   Div    ,  Re_kf ,  Rnu_dns
    !
    !        OPEN(nwf2,file=curdat(1:len_trim(curdat))//'.st',access='append')
    !        WRITE(nwf2,100) 'INITIAL FORCING FIELD'
    !        WRITE(nwf2,110)
    !        WRITE(nwf2,120) 'Truncated Field    ',UXX1,UYY1,UZZ1
    !        !       WRITE(*,120) 'K-normal Projection',UXX2,UYY2,UZZ2
    !        WRITE(nwf2,120) 'Normalised Field   ',FXSQ,FYSQ,FZSQ
    !        WRITE(nwf2,130)  'Totf'  , 'Vchar', 'Tret'
    !        WRITE(nwf2,135)   Totf   ,  Vchar ,  Tret
    !        WRITE(nwf2,130)  'Div(f)', 'Re_kf', 'Rnu_dns'
    !        WRITE(nwf2,135)   Div    ,  Re_kf ,  Rnu_dns
    !        Close(nwf2)
    !
    !        !      WRITE(NWF1,100) 'INITIAL FORCING FIELD'
    !        !      WRITE(NWF1,110)
    !        !      WRITE(NWF1,120) 'Truncated Field    ',UXX1,UYY1,UZZ1
    !        !cc      WRITE(NWF1,120) 'K-normal Projection',UXX2,UYY2,UZZ2
    !        !      WRITE(NWF1,120) 'Normalised Field   ',FXSQ,FYSQ,FZSQ
    !        !      WRITE(NWF1,130)  'Totf'  , 'Vchar', 'Tret'
    !        !      WRITE(NWF1,135)   Totf   ,  Vchar ,  Tret
    !        !      WRITE(NWF1,130)  'Div(f)', 'Re_kf', 'Rnu_dns'
    !        !      WRITE(NWF1,135)   Div    ,  Re_kf ,  Rnu_dns
    !
    !      endif
    !
    !      ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    !      !      CALL specinetiq                                      &
    !      !        (n1d,n2d,n3d,n1,n2,n3,n3b,n3p,nspec,              &
    !      !         fs_x,fs_y,fs_z,xk,energy,ux,uy,uz,fsq1,fsq2,fsq3,f3)
    !      !       CALL specinetiq_PV                        &
    !      !           (n1d,n2d,n3d,n1,n2,n3,n3b,n3p,nspec,  &
    !      !            fs_x,fs_y,fs_z,xk,energy,ux,uy,uz,fsq1, &
    !      !            fsq2,fsq3,f3, 2)
    !      call  SPECINETIQ_PV_N(nspec, &
    !              fs_x,fs_y,fs_z,xk,energy,ux,uy,uz,fsq1L,fsq2L,fsq3L,f3, sp,2)
    !      energy = energy * dt
    !
    !      if(myid.eq.0) then
    !        name =   'force.dat'
    !        OPEN (70,FILE = name,Form = 'FORMATTED')
    !
    !        DO J=1,NSPEC
    !          WRITE(70,9000)  J,ENERGY(J)
    !        ENDDO
    !        9000 FORMAT(I3,3x,1PE13.5)
    !        CLOSE(70)
    !
    !      endif
    !      !       CALL specwrite (energy,nspec,name)
    !
    !    ENDIF
  end subroutine Stats_quant_init
  
  
end module m_aux_phys