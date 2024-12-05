! Module for global parameters.

module m_glob_params
    
    use decomp_2d_constants, only: mytype
    
    implicit none
    
    public
    
    !***** THI_PJET *****
    integer ::N1G,N2G,N3G,P_ROW,P_COL
    integer::NBOUND,ROPTION,NRANDCALL
    real(mytype):: ALPHA0
    integer::NRKU
    real(mytype)::DT
    integer::NTMAX,ICFL
    integer::ISIMU,ISCAR
    real(mytype)::SCH,RNU
    integer::ICHOPA,ISTRUC,ISMAG,IHVISC,NEXPH
    character(len=20)::CURDAT,RESDAT
    integer::JINIT,NI
    integer::JWRITE,JSAVE,NOUT,NOUTT
    integer::IMYTYPEESC,IPDF,NOPDF
    real(mytype)::C1,NSLOPE1,C4,NSLOPE4,TOTE,TOTS
    integer::IFORCE
    real(mytype)::X0,TOTF,WIDTH_F
    integer::ISFORCE
    real(mytype)::X0S,TOTFS,WIDTH_F_SCAL,BETA_SCAL
    real(mytype)::REYH,A0X,A0Y,A0Z,BLY,HUTHETA,U2,U1,HSTHETA,T2,T1
    integer::NJJ,NNRDT
    integer::IPARTIC,INTERPOL_TYPE,NPARTIC_G,PARTIC_LOC,NPART_VAR,NOUTPART,JPINIT
    integer::IVISC
    real(mytype)::LPOLY,RLTIME,POLY_BETA,LAMPOWER,NPOWER,APOWER,NUZERO,NUINFTY
    integer::BUILDANIMATION,TIMESTART,TIMESKIP,NUMBERFRAMES
    integer::WRITEBINARY,NDIMS,ISTART,IEND,JSTART,JEND,KSTART,KEND
    integer::IMHD
    real(mytype)::L0_E,U0_E,N0_E,SIGMA_E


    !***** MPI *****
    integer:: IERR

end module m_glob_params