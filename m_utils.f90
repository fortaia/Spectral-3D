module m_utils
  
  implicit none
  
  private
  
  public :: Dealiasing_init, Dealiasing_trunc, Cfl, Random_number_gen
  
  contains
    
    ! Aliasing removal - TRUNCATION (2/3 rule) - Creating the array of points to truncate
    ! Covering 2 scenarios in each pencil: 1- Truncate all; 2- Truncate some
    !
    subroutine Dealiasing_init(sq_wnumG1,sq_wnumG2,sq_wnumG3,trunc_index,sp)
      
      use m_glob_params
      use decomp_2d_constants, only: mytype
      use decomp_2d, only: decomp_info
      
      integer :: i,j,k
      integer(8) :: count
      integer, allocatable, dimension(:,:), intent(OUT) :: trunc_index
      real(mytype), allocatable, dimension (:), intent(IN) :: sq_wnumG1,sq_wnumG2,sq_wnumG3
      real(mytype) :: kmax,kmax2,wave_sq,sq_w2,sq_w3
      
      type(decomp_info), pointer :: sp
      
      !---------   Calculations start ------------------
      ! kmax - maximum wavenumber allowed
      kmax = ( real(N2G,mytype)/2._mytype * ALPHA0 )
      kmax2=kmax**2
      
      ! Deal with the case in which ALL points must be truncated
      ! Check whether the verteces of the pencil have associated |k|**2 > kmax2
      if(sq_wnumG1(sp%zst(1))+sq_wnumG2(sp%zst(2))+sq_wnumG3(sp%zst(3))>kmax2.and.&
              sq_wnumG1(sp%zst(1))+sq_wnumG2(sp%zen(2))+sq_wnumG3(sp%zst(3))>kmax2 ) return
      
      ! Deal with the case in which SOME points must be truncated
      count=0
      ! Count how many points to truncate
      do k=sp%zst(3),int(kmax)+1
        sq_w3=sq_wnumG3(k)
        do j=sp%zst(2),sp%zen(2)
          sq_w2=sq_wnumG2(j)
          do i=sp%zst(1),sp%zen(1)
            ! k**2
            wave_sq = (sq_wnumG1(i)+ sq_w2 + sq_w3)
            
            if(wave_sq>kmax2) then
              count=count+1
              exit
            end if
          
          enddo
        enddo
      enddo
      
      ! Store indexes of points to truncate
      allocate(trunc_index(count,3))
      
      count=0
      do k=sp%zst(3),int(kmax)+1
        sq_w3=sq_wnumG3(k)
        do j=sp%zst(2),sp%zen(2)
          sq_w2=sq_wnumG2(j)
          do i=sp%zst(1),sp%zen(1)
            ! k**2
            wave_sq = (sq_wnumG1(i)+ sq_w2 + sq_w3)
            
            if(wave_sq>kmax2) then
              count=count+1
              trunc_index(count,:)=[i,j,k]
              exit
            end if
          
          enddo
        enddo
      enddo
    
    end subroutine dealiasing_init
    !
    ! Aliasing removal - TRUNCATION (2/3 rule) - Setting values to 0
    !
    subroutine Dealiasing_trunc(field1, field2, field3, trunc_index,sp)
      
      use m_glob_params
      
      use decomp_2d_constants, only: mytype
      use decomp_2d, only: decomp_info
      
      integer :: i,j,k,idx_kmax,start_int,end_int
      integer(8) :: count
      integer, allocatable, dimension(:,:), intent(IN) :: trunc_index
      complex(mytype), allocatable, dimension(:,:,:), intent(INOUT) :: field1,field2,field3
      type(decomp_info), pointer :: sp
      
      !------------ Start subroutine -------------------
      
      ! If a pencil must be completely truncated
      if(.not. allocated(trunc_index)) then
        field1=(0.0_mytype,0.0_mytype)
        field2=(0.0_mytype,0.0_mytype)
        field3=(0.0_mytype,0.0_mytype)
        return
      end if
      
      idx_kmax=size(trunc_index,1)
      
      ! Truncating i-line by line (along - x)
      
      do count=1,idx_kmax
        i=trunc_index(count,1)
        j=trunc_index(count,2)
        k=trunc_index(count,3)
        ! In the first N3G half
        field1(i:sp%zen(1),j,k)=(0.0_mytype,0.0_mytype)
        field2(i:sp%zen(1),j,k)=(0.0_mytype,0.0_mytype)
        field3(i:sp%zen(1),j,k)=(0.0_mytype,0.0_mytype)
        ! In the second N3G half
        end_int=N3G+1-k
        field1(i:sp%zen(1),j,end_int)=(0.0_mytype,0.0_mytype)
        field2(i:sp%zen(1),j,end_int)=(0.0_mytype,0.0_mytype)
        field3(i:sp%zen(1),j,end_int)=(0.0_mytype,0.0_mytype)
      enddo
      
      !Truncating contiguos "Central Portion"
      idx_kmax = int( real(N2G,mytype)/2.0_mytype * ALPHA0 ) + 1
      start_int=idx_kmax+1
      end_int=N3G-idx_kmax
      
      field1(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),start_int:end_int)=(0.0_mytype,0.0_mytype)
      field2(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),start_int:end_int)=(0.0_mytype,0.0_mytype)
      field3(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),start_int:end_int)=(0.0_mytype,0.0_mytype)
    
    end subroutine dealiasing_trunc
    
    ! CFL - Calculation
    function Cfl(ur,vr,wr) result(dt_cfl)
      
      use m_glob_params
      use decomp_2d_mpi, only: mytype,nrank,real_type
      
      use MPI
      
      integer :: errorcode
      real(mytype), parameter :: TWOPI=6.28318530717958647692528676655900
      real(mytype) :: cfl_conv, cfl_visc, u_lmax, u_gmax, deltax, dt_conv, dt_visc, dt_cfl_loc, dt_cfl
      real(mytype), allocatable, dimension(:,:,:), intent(IN) :: ur,vr,wr
      
      cfl_conv=0.6_mytype
      cfl_visc=0.5_mytype
      
      ! Calculate local MAX : why not * 1/3 ZZZZ
      u_lmax = sqrt(MAXVAL(( ur*ur + vr*vr +  wr*wr )))
      ! Calculate global MAX
      CALL MPI_ALLREDUCE(u_lmax,u_gmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
      
      !-----ISOTROPIC TURBULENCE, TEMPORAL PLANE JET AND WAKE AND ROUND JET, OR RDT SIMULATIONS
      if (ISIMU==0.or.ISIMU==2) then
        deltax  = TWOPI / real(N1G,mytype)
      elseif (ISIMU==1.OR.ISIMU==3.OR.ISIMU==4) then
        deltax  = A0X*TWOPI / real(N1G,mytype)
      endif
      
      dt_conv = (cfl_conv*deltax) / u_gmax
      
      if (RNU==0D0) THEN
        dt_visc = 1.E30
      else
        dt_visc = (cfl_visc*deltax*deltax) / RNU
      endif
      
      dt_cfl_loc = min(dt_conv,dt_visc)
      
      CALL MPI_ALLREDUCE(dt_cfl_loc,dt_cfl,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
      
      if(nrank==0) then
        if (dt<1.E-10) then
          
          print *,'                                        '
          print *,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          print *,'    TIME-STEP TOO SMALL: DIVERGENCE'
          print *,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          
          call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
        
        endif
      endif
      
      return
    
    end function cfl
    
    ! Calculate 1 random numbers
    function Random_number_gen() result(a_rand)
      
      use m_glob_params, only:ROPTION
      use decomp_2d_constants, only: mytype
      
      integer :: time, n
      integer, dimension(:), allocatable :: seed
      real(mytype) :: a_rand
      
      ! option = 0 for testing
      ! option = 1 direct paralel random with myid as seed.
      
      if(ROPTION==0) then
        
        ! not random at all, but used for non-random testing.
        !      aux(:,:,:)    = 0.4D0
        a_rand=0.4_mytype
      endif
      
      if(ROPTION==1) then
        
        ! Get the required seed size
        call RANDOM_SEED(size = n)
        allocate(seed(n))
        
        call SYSTEM_CLOCK(count=time)
        seed = time
        CALL RANDOM_SEED(put=seed)
        
        CALL RANDOM_NUMBER(a_rand)
      
      endif
      return
    
    end function Random_number_gen
  
end module m_utils