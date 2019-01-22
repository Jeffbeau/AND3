


PGF90 (Version     10.5)          09/23/2010  09:49:03      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir /usr/local/include
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /usr/local/include
          -idir /usr/include/mpich2-x86_64
          -idir /usr/include/mpich2-x86_64
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: zdf_oce.F90

(    1) MODULE zdf_oce
(    2)    !!======================================================================
(    3)    !!              ***  MODULE  zdf_oce  ***
(    4)    !! Ocean physics : define vertical mixing variables
(    5)    !!=====================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf_oce.F90,v 1.4 2005/09/02 15:02:47 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !!----------------------------------------------------------------------
(   12)    !!   zdf_init    : initialization, namelist read, and parameters control
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE par_oce         ! mesh and scale factors
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! * Share Module variables
(   21)    LOGICAL, PARAMETER, PUBLIC ::    &   !:
(   22) #if defined key_zdfcst   ||   defined key_esopa
(   23)       lk_zdfcst        = .TRUE.         !: constant vertical mixing flag
(   24) #else
(   25)       lk_zdfcst        = .FALSE.        !: constant vertical mixing flag
(   26) #endif
(   27)    LOGICAL, PUBLIC ::    &   !:
(   28)       ln_zdfevd        = .TRUE.  ,   &  !: convection: enhanced vertical diffusion flag
(   29)       ln_zdfnpc        = .FALSE.        !: convection: non-penetrative convection flag
(   30) 
(   31)    LOGICAL, PUBLIC ::    &   !:
(   32)       l_trazdf_exp     = .FALSE. ,   &  !: ???
(   33)       l_trazdf_imp     = .FALSE. ,   &  !: 
(   34)       l_dynzdf_exp     = .FALSE. ,   &  !: 
(   35)       l_dynzdf_imp     = .TRUE.  ,   &  !:
(   36)       l_dynzdf_imp_tsk = .FALSE.        !:
(   37) 
(   38)    INTEGER, PUBLIC ::    & !!: namzdf:  vertical diffusion
(   39)       n_zdfexp = 3    ,  &  !: number of sub-time step (explicit time stepping)
(   40)       nevdm    = 1          !: =0/1 flag to apply enhanced avm or not
(   41)  
(   42)    REAL(wp), PUBLIC ::   & !!: namzdf   vertical diffusion
(   43)       avm0  = 1.e-4_wp,  &  !: vertical eddy viscosity (m2/s)
(   44)       avt0  = 1.e-5_wp,  &  !: vertical eddy diffusivity (m2/s)
(   45)       avevd = 1._wp         !: vertical eddy coeff. for enhanced vert. diff. (m2/s)
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 2

(   47)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   48)       avmu,              &  !: vertical viscosity coeff. at uw-, vw-points
(   49)       avmv,              &  !: vertical viscosity coeff. at uw-, vw-points
(   50)       avt ,              &  !: vertical diffusivity coeff. at w-point
(   51)       avt_evd,           &  !: convection: enhanced vertical diffusivity coeff. at w-point
(   52)       avmu_evd              !: convection: enhanced vertical viscosity   coeff. at w-point
(   53)  
(   54)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &   !:
(   55)       avmb, avtb            !: background profile of avm and avt
(   56)  
(   57)    !!======================================================================
(   58) END MODULE zdf_oce




















































PGF90 (Version     10.5)          09/23/2010  09:49:03      page 3

(    1) # 1 "zdf_oce.F90"
(    1) MODULE zdf_oce
(    2)    !!======================================================================
(    3)    !!              ***  MODULE  zdf_oce  ***
(    4)    !! Ocean physics : define vertical mixing variables
(    5)    !!=====================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf_oce.F90,v 1.4 2005/09/02 15:02:47 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !!----------------------------------------------------------------------
(   12)    !!   zdf_init    : initialization, namelist read, and parameters control
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE par_oce         ! mesh and scale factors
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! * Share Module variables
(   21)    LOGICAL, PARAMETER, PUBLIC ::    &   !:
(   25)       lk_zdfcst        = .FALSE.        !: constant vertical mixing flag
(   26) # 27
(   27)    LOGICAL, PUBLIC ::    &   !:
(   28)       ln_zdfevd        = .TRUE.  ,   &  !: convection: enhanced vertical diffusion flag
(   29)       ln_zdfnpc        = .FALSE.        !: convection: non-penetrative convection flag
(   30) 
(   31)    LOGICAL, PUBLIC ::    &   !:
(   32)       l_trazdf_exp     = .FALSE. ,   &  !: ???
(   33)       l_trazdf_imp     = .FALSE. ,   &  !: 
(   34)       l_dynzdf_exp     = .FALSE. ,   &  !: 
(   35)       l_dynzdf_imp     = .TRUE.  ,   &  !:
(   36)       l_dynzdf_imp_tsk = .FALSE.        !:
(   37) 
(   38)    INTEGER, PUBLIC ::    & !!: namzdf:  vertical diffusion
(   39)       n_zdfexp = 3    ,  &  !: number of sub-time step (explicit time stepping)
(   40)       nevdm    = 1          !: =0/1 flag to apply enhanced avm or not
(   41)  
(   42)    REAL(wp), PUBLIC ::   & !!: namzdf   vertical diffusion
(   43)       avm0  = 1.e-4_wp,  &  !: vertical eddy viscosity (m2/s)
(   44)       avt0  = 1.e-5_wp,  &  !: vertical eddy diffusivity (m2/s)
(   45)       avevd = 1._wp         !: vertical eddy coeff. for enhanced vert. diff. (m2/s)
(   46) 
(   47)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   48)       avmu,              &  !: vertical viscosity coeff. at uw-, vw-points
(   49)       avmv,              &  !: vertical viscosity coeff. at uw-, vw-points
(   50)       avt ,              &  !: vertical diffusivity coeff. at w-point
(   51)       avt_evd,           &  !: convection: enhanced vertical diffusivity coeff. at w-point
(   52)       avmu_evd              !: convection: enhanced vertical viscosity   coeff. at w-point
(   53)  
(   54)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &   !:
(   55)       avmb, avtb            !: background profile of avm and avt
(   56)  
(   57)    !!======================================================================
(   58) END MODULE zdf_oce
