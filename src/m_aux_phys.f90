!
! MODULE purpose: Group useful procedures for calculations in physical space
!
module m_aux_phys
  
  use m_glob_params
  use decomp_2d_constants, only: real_type
  
  implicit none
  
  private
  
  public :: Variance_ph, Stats_quant_init, Calc_kin_en_ph
  
  ! Generic interface to handle multiple data types
  interface Variance_ph
    module procedure Variance_ph_1field
    module procedure Variance_ph_2fields
  end interface

contains
  
  ! Calculate variance var(u)=mean(u**2)-mean(u)**2
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
  
  ! Calculate variance var(uv)=mean(u*v)-mean(u)*mean(v)
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
  
  ! Calculate reference quantities at the beginning of a new simulation
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
    
    uprime2=Variance_ph(ur)
    vprime2=Variance_ph(vr)
    wprime2=Variance_ph(wr)
    
    if (ISIMU==0)then
      Vchar = (TOTF / X0) ** 0.3333333333333333333333333333333_mytype
      Tref = 1.0_mytype/(X0 * Vchar) !  LARGE-EDDY TURNOVER TIME
    else
      Vchar = sqrt(uprime2 + vprime2 + wprime2)
      Tref  = TWOPI / Vchar
    end if
  
  end subroutine Stats_quant_init
  
  subroutine Calc_kin_en_ph(ur,vr,wr,normaliz)
    
    use decomp_2d_mpi, only: mytype, real_type, nrank
    use MPI
    
    real(mytype), contiguous, dimension(:,:,:), intent(IN) :: ur,vr,wr
    
    real(mytype) :: kin_en_loc, kin_en, normaliz
    
    !---------   Subroutine start  ------------------
    
    kin_en_loc=sum(ur*ur)+sum(vr*vr)+sum(wr*wr)
    call MPI_Reduce(kin_en_loc, kin_en, 1, real_type, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    if (nrank==0) print*, kin_en/(2.0_mytype*normaliz)
  
  end subroutine Calc_kin_en_ph

end module m_aux_phys