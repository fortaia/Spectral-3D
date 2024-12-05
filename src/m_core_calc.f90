module m_core_calc
  
  use m_glob_params
  use m_initial_cond, only: ur,vr,wr,u,v,w
  use m_aux_spect
  use m_utils
  use m_equations_terms
  
  use decomp_2d_io
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d_fft
  use decomp_2d
  
  use MPI
  
  implicit none
  
  private
  
  public :: Core_calc
  
  real(mytype), allocatable, dimension (:), private :: wavenumG1,wavenumG2,wavenumG3
  real(mytype), allocatable, dimension (:), private :: sq_wnumG1,sq_wnumG2,sq_wnumG3
  real(mytype), private :: normaliz

contains
  
  subroutine Core_calc(ph,sp)
    
    use m_aux_phys, only: Stats_quant_init
    use m_io, only: Write_3d_velocities, Write_store_file, Read_store_file,Write_screen_forcing
    
    integer :: iter
    real(mytype) :: const_rk3(5), time, t_start, t_end, t_ave, tref
    character(len=5) :: iter_string
    integer(8), dimension(:), allocatable :: linear_index
    integer, allocatable, dimension(:,:) :: trunc_index
    real(mytype), dimension(:,:), allocatable :: forc_init
    ! In Physical space: Working arrays
    real(mytype), allocatable, dimension (:,:,:) :: work_xph,work_yph,work_zph
    ! In Spectral space: Working arrays
    complex(mytype), allocatable, dimension (:,:,:) :: rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1
    complex(mytype), allocatable, dimension (:,:,:) :: fs_x,fs_y,fs_z
    
    type(decomp_info), pointer :: ph,sp
    !------------- Start subroutine -----------------------------
    
    if (nrank==0) then
      
      write(*,*) '----------------------------------------------------'
      write(*,*)
      write(*,*) '  Starting CALCULATIONS                             '
      write(*,*)
      write(*,*) '----------------------------------------------------'
    endif
    
    call wavenum(wavenumG1,wavenumG2,wavenumG3)
    
    normaliz = real(N1G*N2G*N3G,mytype)
    
    ! Calculating useful constants before the loop
    allocate(sq_wnumG1(size(wavenumG1)),sq_wnumG2(size(wavenumG2)),sq_wnumG3(size(wavenumG3)))
    sq_wnumG1=wavenumG1**2;sq_wnumG2=wavenumG2**2;sq_wnumG3=wavenumG3**2
    
    ! 3-steps Runge-Kutta coefficients
    const_rk3=[dt/3._mytype,-5._mytype/9._mytype,15._mytype/16._mytype * dt,- 153._mytype/128._mytype,8._mytype/15._mytype*dt]
    
    ! Determine points to truncate
    call Dealiasing_init(sq_wnumG1,sq_wnumG2,sq_wnumG3,trunc_index,sp)
    
    !Allocations needed for the working arrays
    call alloc_z(rhs_x,sp,.true.)
    call alloc_z(rhs_y,sp,.true.)
    call alloc_z(rhs_z,sp,.true.)
    call alloc_z(work_xsp1,sp,.true.)
    call alloc_z(work_ysp1,sp,.true.)
    call alloc_z(work_zsp1,sp,.true.)
    call alloc_x(work_xph,ph,.true.)
    call alloc_x(work_yph,ph,.true.)
    call alloc_x(work_zph,ph,.true.)
    
    ! In case of active forcing
    if(IFORCE==1) then
      call alloc_z(fs_x,sp,.true.)
      call alloc_z(fs_y,sp,.true.)
      call alloc_z(fs_z,sp,.true.)
      
      fs_x = (0.0_mytype,0.0_mytype)
      fs_y = (0.0_mytype,0.0_mytype)
      fs_z = (0.0_mytype,0.0_mytype)
      
      call Forcing_init(wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,linear_index,forc_init,sp)
    end if
    
    ! In case of a new run apply continuity, dealiasing, and compute ref. quant.
    if(JINIT==0) then
      
      call Dealiasing_trunc(u, v, w, trunc_index,sp)
      call Continuity(u,v,w,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
      ! Writing initial fields
      call Write_3d_velocities(trim(RESDAT)//"_00000",u,v,w,sp)
      ! Calculating initial variables
      work_xsp1=u
      call decomp_2d_fft_3d(work_xsp1,ur)
      work_xsp1=v
      call decomp_2d_fft_3d(work_xsp1,vr)
      work_xsp1=w
      call decomp_2d_fft_3d(work_xsp1,wr)
      
      call Stats_quant_init(ur,vr,wr,tref)
      
      time = 0.0_mytype
    else
      ! In case of restart load previuosly calculated quantities
      call Read_store_file(time,tref,dt)
      if(IFORCE==1) call Random_number_gen_init()
    endif
    
    call Stat_run(work_xsp1,0,time,sp)
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    t_ave=0.0_mytype
    
    ! -------------------------- CORE CYCLE ** START -------------------------------
    do iter=1,NTMAX
      ! Time start to compute iteration time
      t_start = MPI_WTIME()
      
      work_xsp1=u
      call decomp_2d_fft_3d(work_xsp1,ur)
      work_xsp1=v
      call decomp_2d_fft_3d(work_xsp1,vr)
      work_xsp1=w
      call decomp_2d_fft_3d(work_xsp1,wr)
      
      ! Time step controlled by CFL
      if(icfl==1) then
        dt=Cfl(ur,vr,wr)
        ! 3-steps Runge-Kutta coefficients
        const_rk3(1) = dt* 0.3333333333333333333333333333333333333333333333333333333333333333_mytype !1/3
        const_rk3(3) = dt* 0.9375_mytype !15/16
        const_rk3(5) = dt* 0.5333333333333333333333333333333333333333333333333333333333333333_mytype !8/15
      end if
      
      ! Time Advancement - CORE CALCULATIONS
      if(IFORCE /= 1) then
        ! NO FORCING
        call Rk3(rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph,const_rk3,trunc_index,ph,sp)
      else
        ! ACTIVE FORCING
        call Rk3_f(rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph,fs_x,fs_y,fs_z,&
                const_rk3,trunc_index,linear_index,forc_init,ph,sp)
      end if
      
      ! Calculating and writing statistics during time advancement
      if (mod(iter,NOUT)==0) call Stat_run(work_xsp1,iter,time,sp)
      ! Calculating and writing planes during time advancement
      if (mod(iter,NOUTT)==0) then
        if(ISIMU/=0) then
          ! Print U-velocity
          work_xsp1=u
          call decomp_2d_fft_3d(work_xsp1,ur)
          write(iter_string, '(I5.5)') NTMAX
          call decomp_2d_write_plane(1, ur, 'plane_U_'//iter_string//'.raw',&
                  opt_iplane=N3G / 2, opt_decomp=ph, opt_reduce_prec=.false.)
          
          ! Calculate and print enstrophy
          call Calc_vorticity_sp(u,v,w,wavenumG1,wavenumG2,wavenumG3,work_xsp1,work_ysp1,work_zsp1,sp)
          
          call decomp_2d_fft_3d(work_xsp1,work_xph)
          call decomp_2d_fft_3d(work_ysp1,work_yph)
          call decomp_2d_fft_3d(work_zsp1,work_zph)
          
          work_xph=(work_xph*work_xph + work_yph*work_yph + work_zph*work_zph) * 0.5_mytype
          call decomp_2d_write_plane(1, ur, 'plane_enstr_'//iter_string//'.raw',&
                  opt_iplane=N3G / 2, opt_decomp=ph, opt_reduce_prec=.false.)
        end if
      end if
      ! Writing 3D fields at prescribed iteration
      if (mod(iter,JSAVE)==0.and.iter/=NTMAX) then
        write(iter_string, '(I5.5)') iter
        call Write_3d_velocities(trim(RESDAT)//"_"//iter_string,u,v,w,sp)
        if(nrank==0) call Write_store_file('store_'//iter_string//'.info',time,tref,dt)
      end if
      
      time = time + dt
      
      t_end = MPI_WTIME() - t_start
      t_ave = t_ave + t_end
      
      ! Write to screen
      if (nrank==0) then
        
        write(*,110) time, time/tref
        write(*,'(A13,I5,2x,A20,F7.3)') 'Iteration = ',iter,'|| Iteration time =', t_end
        write(*,'(A20,F7.3)') 'Average iter. time:', t_ave/iter
        write(*,*) '---------------------------------------------------'
        
        110 FORMAT(1x,'Time = ',1PG11.4,2x,'|| T/Tref = ', 1PG11.4)
      
      endif
      
      if(IFORCE==1) then
        call Write_screen_forcing(u,v,w,fs_x,fs_y,fs_z,sp)
      end if
    
    end do
    ! -------------------------- CORE CYCLE ** END -------------------------------

    if (nrank==0) then
      
      call Write_store_file("store.info",time,tref,dt)
      
      write(*,*) '---------------------------------------------------'
      write(*,*)
      write(*,*) '  ENDING CALCULATIONS                              '
      write(*,*)
      write(*,*) '---------------------------------------------------'
    
    endif
    
    !--------- Writing planes for all ISIMU --------------------------------
    ! Writing 2 planes of U-velocity and Enstrophy
    ! Print U-velocity
    work_xsp1=u
    call decomp_2d_fft_3d(work_xsp1,ur)
    
    write(iter_string, '(I5.5)') NTMAX
    call decomp_2d_write_plane(1, ur, 'plane_U_'//iter_string//'.raw',&
            opt_iplane=N3G / 2, opt_decomp=ph, opt_reduce_prec=.false.)
    
    ! Calculate and print enstrophy
    call Calc_vorticity_sp(u,v,w,wavenumG1,wavenumG2,wavenumG3,work_xsp1,work_ysp1,work_zsp1,sp)
    
    call decomp_2d_fft_3d(work_xsp1,work_xph)
    call decomp_2d_fft_3d(work_ysp1,work_yph)
    call decomp_2d_fft_3d(work_zsp1,work_zph)
    
    work_xph=(work_xph*work_xph + work_yph*work_yph + work_zph*work_zph) * 0.5_mytype
    
    call decomp_2d_write_plane(1, ur, 'plane_enstr_'//iter_string//'.raw',&
            opt_iplane=N3G / 2, opt_decomp=ph, opt_reduce_prec=.false.)
    ! ---------------------------------------------------------------------------
    
    ! DEALLOCATING arrays
    deallocate(ur,vr,wr,rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph)
    deallocate(wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3)
    if (IFORCE==1) deallocate(fs_x,fs_y,fs_z)
    
  end subroutine Core_calc
  
  subroutine Stat_run(work_xsp1,iter,time,sp)
    
    use m_io, only: Write_ascii
    
    integer,intent(IN) :: iter
    integer, dimension(N2G/2+1) :: count,count_g
    real(mytype),intent(IN) :: time
    real(mytype) :: kin_en, diss, u_rms_sq, rey_lambda, eta, l_int, kmax_eta, two_third, PI_half
    real(mytype), dimension(N2G/2+1) :: sh_ave,en_spect
    complex(mytype), dimension (:,:,:) :: work_xsp1
    character(len=5) :: iter_string
    type(decomp_info), pointer :: sp
    
    !------- Start procedure -------------------
    call Calc_kin_en_sp(u,v,w,work_xsp1)
    call Shell_average(work_xsp1,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp,count,sh_ave)
    
    call MPI_Reduce(count, count_g, N2G/2+1, MPI_INTEGER, MPI_SUM, nproc-1, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(sh_ave, en_spect, N2G/2+1, real_type, MPI_SUM, nproc-1, MPI_COMM_WORLD, ierr)
    
    if(nrank==nproc-1) then
      en_spect=en_spect/real(count_g,mytype)
      write(iter_string, '(I5.5)') iter
      call Write_ascii('k_Ek_'//iter_string//'.dat',wavenumG1,en_spect)
      
      ! Average kinetic energy : approximating integral with trapz
      kin_en=sum(en_spect)-0.5_mytype*(en_spect(1)+en_spect(N2G/2+1))
      ! Average Dissipation : approximating integral with trapz (note that wavenumG1(1)=0)
      diss = 2.0_mytype*RNU*(sum(en_spect*wavenumG1*wavenumG1)-0.5_mytype*en_spect(N2G/2+1)*wavenumG1(N2G/2+1)**2)
      ! Ave. Reynolds num. (based on Taylor length= sqrt(15*RNU/diss)*u)
      ! u**2=(2/3 * kin_en);
      two_third = 0.6666666666666666666666666666666666_mytype
      u_rms_sq = kin_en*two_third
      rey_lambda = sqrt(15.0_mytype/(diss*RNU))*u_rms_sq
      ! Integral scale
      PI_half = 1.5707963267948966192313216916397514420985846996875529104874722961_mytype
      l_int = (sum(en_spect(2:N2G/2+1)/wavenumG1(2:N2G/2+1))-0.5_mytype*en_spect(N2G/2+1)/wavenumG1(N2G/2+1))/&
              u_rms_sq * PI_half
      ! Average Kolmogorv scale
      eta = (RNU**3/diss)**0.25_mytype
      ! Average resolution
      kmax_eta = eta * wavenumG1(N2G/2+1) * X0
      call Write_ascii('Time_Ekin_Diss_ReyLam_L_KmaxEta.dat',[time],[kin_en],[diss],[rey_lambda],[l_int],[kmax_eta])
      
    end if
  
  end subroutine Stat_run
  
  !
  ! 3 step Runge-Kutta explicit time advancement
  !
  subroutine Rk3(rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph,const_rk3,trunc_index,ph,sp)
    
    integer, allocatable, dimension(:,:), intent(IN) :: trunc_index
    real(mytype), intent(IN) :: const_rk3(5)
    ! In Physical space: Working arrays
    real(mytype), dimension (:,:,:) :: work_xph,work_yph,work_zph
    ! In Spectral space
    complex(mytype), dimension (:,:,:) :: rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1
    
    type(decomp_info), pointer :: ph,sp
    
    !-----------------------------------------------------------
    ! 1st - step
    call RHS(u,v,w,ur,vr,wr,rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph, &
            wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,normaliz,trunc_index,ph,sp)
    
    u = u + const_rk3(1) * rhs_x
    v = v + const_rk3(1) * rhs_y
    w = w + const_rk3(1) * rhs_z
    
    call Continuity(u,v,w,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
    ! 2nd - step
    call Dissipative(u,v,w,sq_wnumG1,sq_wnumG2,sq_wnumG3,work_xsp1,work_ysp1,work_zsp1,sp)
    
    rhs_x = const_rk3(2) * rhs_x + work_xsp1 * dt
    rhs_y = const_rk3(2) * rhs_y + work_ysp1 * dt
    rhs_z = const_rk3(2) * rhs_z + work_zsp1 * dt
    
    call NonLin(u,v,w,ur,vr,wr,work_xph,work_yph,work_zph,work_xsp1,work_ysp1,work_zsp1,&
            wavenumG1,wavenumG2,wavenumG3,normaliz,trunc_index,ph,sp)
    
    rhs_x = rhs_x + work_xsp1 * dt
    rhs_y = rhs_y + work_ysp1 * dt
    rhs_z = rhs_z + work_zsp1 * dt
    
    u = u + const_rk3(3) * rhs_x
    v = v + const_rk3(3) * rhs_y
    w = w + const_rk3(3) * rhs_z
    
    call Continuity(u,v,w,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
    ! 3rd - step
    call Dissipative(u,v,w,sq_wnumG1,sq_wnumG2,sq_wnumG3,work_xsp1,work_ysp1,work_zsp1,sp)
    
    rhs_x = const_rk3(4) * rhs_x + work_xsp1 * dt
    rhs_y = const_rk3(4) * rhs_y + work_ysp1 * dt
    rhs_z = const_rk3(4) * rhs_z + work_zsp1 * dt
    
    call NonLin(u,v,w,ur,vr,wr,work_xph,work_yph,work_zph,work_xsp1,work_ysp1,work_zsp1,&
            wavenumG1,wavenumG2,wavenumG3,normaliz,trunc_index,ph,sp)
    
    rhs_x = rhs_x + work_xsp1 * dt
    rhs_y = rhs_y + work_ysp1 * dt
    rhs_z = rhs_z + work_zsp1 * dt
    
    u = u + const_rk3(5) * rhs_x
    v = v + const_rk3(5) * rhs_y
    w = w + const_rk3(5) * rhs_z
    
    call Continuity(u,v,w,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
  
  end subroutine Rk3
  
  !
  ! 3 step Runge-Kutta explicit time advancement with forcing
  !
  subroutine Rk3_f(rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph,fs_x,fs_y,fs_z,&
          const_rk3,trunc_index,linear_index,forc_init,ph,sp)
    
    integer(8), dimension(:), allocatable, intent(IN) :: linear_index
    integer, allocatable, dimension(:,:), intent(IN) :: trunc_index
    real(mytype), dimension(:,:), allocatable, intent(IN) :: forc_init
    real(mytype), intent(IN) :: const_rk3(5)
    ! In Physical space: Working arrays
    real(mytype), dimension (:,:,:) :: work_xph,work_yph,work_zph
    ! In Spectral space
    complex(mytype), dimension (:,:,:) :: rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1
    complex(mytype), dimension (:,:,:), intent(INOUT) :: fs_x,fs_y,fs_z
    type(decomp_info), pointer :: ph,sp
    
    !-----------------------------------------------------------
    ! 1st - step
    call RHS(u,v,w,ur,vr,wr,rhs_x,rhs_y,rhs_z,work_xsp1,work_ysp1,work_zsp1,work_xph,work_yph,work_zph, &
            wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,normaliz,trunc_index,ph,sp)
    
    call Forcing(u,v,w,fs_x,fs_y,fs_z,linear_index,forc_init,sp)
    
    u = u + const_rk3(1) * (rhs_x + fs_x)
    v = v + const_rk3(1) * (rhs_y + fs_y)
    w = w + const_rk3(1) * (rhs_z + fs_z)
    
    call Continuity(u,v,w,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
    
    ! 2nd - step
    call Dissipative(u,v,w,sq_wnumG1,sq_wnumG2,sq_wnumG3,work_xsp1,work_xsp1,work_xsp1,sp)
    
    rhs_x = const_rk3(2) * (rhs_x + fs_x) + work_xsp1 * dt
    rhs_y = const_rk3(2) * (rhs_y + fs_y) + work_ysp1 * dt
    rhs_z = const_rk3(2) * (rhs_z + fs_z) + work_zsp1 * dt
    
    call NonLin(u,v,w,ur,vr,wr,work_xph,work_yph,work_zph,work_xsp1,work_ysp1,work_zsp1,&
            wavenumG1,wavenumG2,wavenumG3,normaliz,trunc_index,ph,sp)
!                      nonlin  + forc
    rhs_x = rhs_x +  (work_xsp1 + fs_x) * dt
    rhs_y = rhs_y +  (work_ysp1 + fs_y) * dt
    rhs_z = rhs_z +  (work_zsp1 + fs_z) * dt
    
    u = u + const_rk3(3) * rhs_x
    v = v + const_rk3(3) * rhs_y
    w = w + const_rk3(3) * rhs_z
    
    call Continuity(u,v,w,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
    
    ! 3rd - step
    call Dissipative(u,v,w,sq_wnumG1,sq_wnumG2,sq_wnumG3,work_xsp1,work_xsp1,work_xsp1,sp)
    
    rhs_x = const_rk3(4) * rhs_x + work_xsp1 * dt
    rhs_y = const_rk3(4) * rhs_y + work_ysp1 * dt
    rhs_z = const_rk3(4) * rhs_z + work_zsp1 * dt
    
    call NonLin(u,v,w,ur,vr,wr,work_xph,work_yph,work_zph,work_xsp1,work_ysp1,work_zsp1,&
            wavenumG1,wavenumG2,wavenumG3,normaliz,trunc_index,ph,sp)
    !                  nonlin  + forc
    rhs_x = rhs_x +  (work_xsp1 + fs_x) * dt
    rhs_y = rhs_y +  (work_ysp1 + fs_y) * dt
    rhs_z = rhs_z +  (work_zsp1 + fs_z) * dt
    
    u = u + const_rk3(5) * rhs_x
    v = v + const_rk3(5) * rhs_y
    w = w + const_rk3(5) * rhs_z
    
    call Continuity(u,v,w,wavenumG1,wavenumG2,wavenumG3,sq_wnumG1,sq_wnumG2,sq_wnumG3,sp)
  
  end subroutine Rk3_f

end module m_core_calc