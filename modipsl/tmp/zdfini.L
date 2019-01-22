


PGF90 (Version     10.5)          09/23/2010  09:50:22      page 1

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

Filename: zdfini.F90

(    1) MODULE zdfini
(    2)    !!======================================================================
(    3)    !!              ***  MODULE  zdfini  ***
(    4)    !! Ocean physics : define vertical mixing variables
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   zdf_init    : initialization, namelist read, and parameters control
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE par_oce         ! mesh and scale factors
(   12)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   13)    USE ldfdyn_oce      ! ocean dynamics lateral physics
(   14)    USE zdf_oce         ! TKE vertical mixing          
(   15)    USE lib_mpp         ! distribued memory computing
(   16)    USE zdftke          ! TKE vertical mixing  
(   17)    USE zdfkpp          ! KPP vertical mixing          
(   18)    USE zdfddm          ! double diffusion mixing      
(   19)    USE zdfevd          ! enhanced vertical diffusion  
(   20)    USE zdfric          ! Richardson vertical mixing   
(   21)    USE tranpc          ! convection: non penetrative adjustment
(   22)    USE ldfslp          ! iso-neutral slopes
(   23) 
(   24)    USE in_out_manager  ! I/O manager
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! *  Routine accessibility
(   30)    PUBLIC zdf_init          ! routine called by opa.F90
(   31)    !!----------------------------------------------------------------------
(   32)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfini.F90,v 1.4 2005/09/22 10:57:33 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36)    
(   37) CONTAINS
(   38) 
(   39)    SUBROUTINE zdf_init
(   40)       !!----------------------------------------------------------------------
(   41)       !!                  ***  ROUTINE zdf_init  ***
(   42)       !! 
(   43)       !! ** Purpose :   initializations of the vertical ocean physics
(   44)       !!
(   45)       !! ** Method  :   Read namelist namzdf, control logicals 
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 2

(   47)       !! History :
(   48)       !!        !  97-06  (G. Madec)  Original code from inimix
(   49)       !!   8.5  !  02-08  (G. Madec)  F90 : free form
(   50)       !!   9.0  !  05-06  (C. Ethe) KPP parameterization
(   51)       !!----------------------------------------------------------------------
(   52)       !! * Local declarations
(   53)       INTEGER ::   ioptio       ! temporary scalar
(   54)       LOGICAL ::           & !!! namzdf: vertical diffusion
(   55)          ln_zdfexp = .FALSE.   ! explicit vertical diffusion scheme flag
(   56) 
(   57)       !! * Namelist
(   58)       NAMELIST/namzdf/ ln_zdfevd, ln_zdfnpc,   &
(   59)          &             avm0, avt0, avevd, nevdm, ln_zdfexp, n_zdfexp
(   60)       !!----------------------------------------------------------------------
(   61)       !!  OPA 9.0, LODYC-IPSL (2003)
(   62)       !!----------------------------------------------------------------------
(   63) 
(   64)       ! Read namzdf namelist : vertical mixing parameters
(   65)       ! --------------------
(   66)       REWIND( numnam )
(   67)       READ  ( numnam, namzdf )
(   68) 
(   69)       ! Parameter print
(   70)       ! ---------------
(   71)       IF(lwp) THEN
(   72)          WRITE(numout,*)
(   73)          WRITE(numout,*) 'zdf_init: vertical physics'
(   74)          WRITE(numout,*) '~~~~~~~~'
(   75)          WRITE(numout,*) '          Namelist namzdf : set vertical mixing mixing parameters'
(   76)          WRITE(numout,*) '             enhanced vertical diffusion      ln_zdfevd = ', ln_zdfevd
(   77)          WRITE(numout,*) '             non-penetrative convection       ln_zdfnpc = ', ln_zdfnpc
(   78)          WRITE(numout,*) '             vertical eddy viscosity             avm0   = ', avm0
(   79)          WRITE(numout,*) '             vertical eddy diffusivity           avt0   = ', avt0
(   80)          WRITE(numout,*) '             vertical coefficient for evd        avevd  = ', avevd
(   81)          WRITE(numout,*) '                applied on momentum (=1/0)       nevdm  = ', nevdm
(   82)          WRITE(numout,*) '             time splitting / backward scheme ln_zdfexp = ', ln_zdfexp
(   83)          WRITE(numout,*) '             number of time step               n_zdfexp = ', n_zdfexp
(   84)       ENDIF
(   85) 
(   86)       ! Parameter & logicals controls
(   87)       ! -----------------------------
(   88)       ! ... vertical mixing
(   89)       ! time stepping scheme (N.B. TKE or KPP schemes => force the use of implicit scheme)
(   90)       IF( ( ln_zdfexp .AND. .NOT.lk_zdftke ) .OR. ( ln_zdfexp .AND. .NOT.lk_zdfkpp ) ) THEN  
(   91)          l_trazdf_exp = .TRUE.           ! use explicit scheme
(   92)          l_trazdf_imp = .FALSE.
(   93)          l_dynzdf_exp = .TRUE.           ! use explicit scheme
(   94)          l_dynzdf_imp = .FALSE.
(   95)       ELSE
(   96)          l_trazdf_exp = .FALSE.          ! use implicit scheme
(   97)          l_trazdf_imp = .TRUE. 
(   98)          l_dynzdf_exp = .FALSE.          ! use implicit scheme
(   99)          l_dynzdf_imp = .TRUE. 
(  100)       ENDIF
(  101)       IF( l_trazdf_iso .OR. l_trazdf_iso_vo ) THEN  
(  102)          l_trazdf_exp = .FALSE.          ! iso-neutral diffusion : 
(  103)          l_trazdf_imp = .FALSE.          ! implicit scheme included in iso-neutral routine
(  104)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 3

(  105)       IF( l_dynldf_iso ) THEN  
(  106)          l_dynzdf_exp = .FALSE.          ! iso-neutral diffusion :
(  107)          l_dynzdf_imp = .FALSE.          ! implicit scheme included in iso-neutral routine
(  108)       ENDIF
(  109) #if defined key_autotasking
(  110)       IF( l_dynzdf_imp ) THEN
(  111)          l_dynzdf_imp     = .FALSE.
(  112)          l_dynzdf_imp_tsk = .TRUE.
(  113)       ENDIF
(  114) #else
(  115)       l_dynzdf_imp_tsk = .FALSE.
(  116) #endif
(  117)       IF( lk_esopa  ) THEN
(  118)          l_trazdf_exp = .TRUE.           ! esopa: use all options
(  119)          l_trazdf_imp = .TRUE.
(  120)          l_dynzdf_exp     = .TRUE.           ! esopa: use all options
(  121)          l_dynzdf_imp     = .TRUE.
(  122)          l_dynzdf_imp_tsk = .TRUE.
(  123)       ENDIF
(  124) 
(  125)       IF(lwp) WRITE(numout,*)
(  126)       IF(lwp) WRITE(numout,*) '          vertical mixing option :'
(  127)       ioptio = 0
(  128)       IF( lk_zdfcst ) THEN
(  129)          IF(lwp) WRITE(numout,*) '             constant eddy diffusion coef.'
(  130)          ioptio = ioptio+1
(  131)       ENDIF
(  132)       IF( lk_zdfric ) THEN
(  133)          IF(lwp) WRITE(numout,*) '             Richardson dependent eddy coef.'
(  134)          ioptio = ioptio+1
(  135)       ENDIF
(  136)       IF( lk_zdftke ) THEN
(  137)          IF(lwp) WRITE(numout,*) '             TKE dependent eddy coef.'
(  138)          ioptio = ioptio+1
(  139)       ENDIF
(  140)       IF( lk_zdfkpp ) THEN
(  141)          IF(lwp) WRITE(numout,*) '             KPP dependent eddy coef.'
(  142)          ioptio = ioptio+1
(  143)       ENDIF
(  144)       IF( ioptio == 0 .OR. ioptio > 1 .AND. .NOT. lk_esopa ) THEN
(  145)           IF(lwp) WRITE(numout,cform_err)
(  146)           IF(lwp) WRITE(numout,*) ' one and only one vertical diffusion option has to be defined '
(  147)           nstop = nstop + 1
(  148)       ENDIF
(  149) 
(  150)       ! ... Convection
(  151)       IF(lwp) WRITE(numout,*)
(  152)       IF(lwp) WRITE(numout,*) '          convection :'
(  153)       ioptio = 0
(  154)       IF( ln_zdfnpc ) THEN
(  155)          IF(lwp) WRITE(numout,*) '             use non penetrative convective scheme'
(  156)          ioptio = ioptio+1
(  157)       ENDIF
(  158)       IF( ln_zdfevd ) THEN
(  159)          IF(lwp) WRITE(numout,*) '             use enhanced vertical dif. scheme'
(  160)          ioptio = ioptio+1
(  161)       ENDIF
(  162)       IF( lk_zdftke ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 4

(  163)          IF(lwp) WRITE(numout,*) '             use the 1.5 turbulent closure'
(  164)       ENDIF
(  165)       IF( lk_zdfkpp ) THEN
(  166)          IF(lwp) WRITE(numout,*) '             use the KPP closure scheme'
(  167)          IF(lk_mpp) THEN
(  168)             IF(lwp) WRITE(numout,cform_err)
(  169)             IF(lwp) WRITE(numout,*) '             The KPP scheme is not ready to run in MPI'
(  170)          ENDIF
(  171)       ENDIF
(  172)       IF ( ioptio > 1 .AND. .NOT. lk_esopa ) THEN
(  173)           IF(lwp) WRITE(numout,cform_err)
(  174)           IF(lwp) WRITE(numout,*) ' chose between ln_zdfnpc'
(  175)           IF(lwp) WRITE(numout,*) '           and ln_zdfevd'
(  176)           nstop = nstop + 1
(  177)       ENDIF
(  178)       IF( ioptio == 0 .AND. .NOT. lk_zdftke ) THEN
(  179)          IF(lwp) WRITE(numout,cform_err)
(  180)          IF(lwp) WRITE(numout,*) ' except for TKE scheme, a convection scheme is'
(  181)          IF(lwp) WRITE(numout,*) ' required: ln_zdfevd or ln_zdfnpc logicals'
(  182)          nstop = nstop + 1
(  183)       ENDIF
(  184) 
(  185)    END SUBROUTINE zdf_init
(  186) 
(  187)    !!======================================================================
(  188) END MODULE zdfini






































PGF90 (Version     10.5)          09/23/2010  09:50:22      page 5

(    1) # 1 "zdfini.F90"
(    1) MODULE zdfini
(    2)    !!======================================================================
(    3)    !!              ***  MODULE  zdfini  ***
(    4)    !! Ocean physics : define vertical mixing variables
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   zdf_init    : initialization, namelist read, and parameters control
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE par_oce         ! mesh and scale factors
(   12)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   13)    USE ldfdyn_oce      ! ocean dynamics lateral physics
(   14)    USE zdf_oce         ! TKE vertical mixing          
(   15)    USE lib_mpp         ! distribued memory computing
(   16)    USE zdftke          ! TKE vertical mixing  
(   17)    USE zdfkpp          ! KPP vertical mixing          
(   18)    USE zdfddm          ! double diffusion mixing      
(   19)    USE zdfevd          ! enhanced vertical diffusion  
(   20)    USE zdfric          ! Richardson vertical mixing   
(   21)    USE tranpc          ! convection: non penetrative adjustment
(   22)    USE ldfslp          ! iso-neutral slopes
(   23) 
(   24)    USE in_out_manager  ! I/O manager
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! *  Routine accessibility
(   30)    PUBLIC zdf_init          ! routine called by opa.F90
(   31)    !!----------------------------------------------------------------------
(   32)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfini.F90,v 1.4 2005/09/22 10:57:33 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36)    
(   37) CONTAINS
(   38) 
(   39)    SUBROUTINE zdf_init
(   40)       !!----------------------------------------------------------------------
(   41)       !!                  ***  ROUTINE zdf_init  ***
(   42)       !! 
(   43)       !! ** Purpose :   initializations of the vertical ocean physics
(   44)       !!
(   45)       !! ** Method  :   Read namelist namzdf, control logicals 
(   46)       !!
(   47)       !! History :
(   48)       !!        !  97-06  (G. Madec)  Original code from inimix
(   49)       !!   8.5  !  02-08  (G. Madec)  F90 : free form
(   50)       !!   9.0  !  05-06  (C. Ethe) KPP parameterization
(   51)       !!----------------------------------------------------------------------
(   52)       !! * Local declarations
(   53)       INTEGER ::   ioptio       ! temporary scalar
(   54)       LOGICAL ::           & !!! namzdf: vertical diffusion
(   55)          ln_zdfexp = .FALSE.   ! explicit vertical diffusion scheme flag
(   56) 
(   57)       !! * Namelist






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 6

(   58)       NAMELIST/namzdf/ ln_zdfevd, ln_zdfnpc,   &
(   59)          &             avm0, avt0, avevd, nevdm, ln_zdfexp, n_zdfexp
(   60)       !!----------------------------------------------------------------------
(   61)       !!  OPA 9.0, LODYC-IPSL (2003)
(   62)       !!----------------------------------------------------------------------
(   63) 
(   64)       ! Read namzdf namelist : vertical mixing parameters
(   65)       ! --------------------
(   66)       REWIND( numnam )
(   67)       READ  ( numnam, namzdf )
(   68) 
(   69)       ! Parameter print
(   70)       ! ---------------
(   71)       IF(lwp) THEN
(   72)          WRITE(numout,*)
(   73)          WRITE(numout,*) 'zdf_init: vertical physics'
(   74)          WRITE(numout,*) '~~~~~~~~'
(   75)          WRITE(numout,*) '          Namelist namzdf : set vertical mixing mixing parameters'
(   76)          WRITE(numout,*) '             enhanced vertical diffusion      ln_zdfevd = ', ln_zdfevd
(   77)          WRITE(numout,*) '             non-penetrative convection       ln_zdfnpc = ', ln_zdfnpc
(   78)          WRITE(numout,*) '             vertical eddy viscosity             avm0   = ', avm0
(   79)          WRITE(numout,*) '             vertical eddy diffusivity           avt0   = ', avt0
(   80)          WRITE(numout,*) '             vertical coefficient for evd        avevd  = ', avevd
(   81)          WRITE(numout,*) '                applied on momentum (=1/0)       nevdm  = ', nevdm
(   82)          WRITE(numout,*) '             time splitting / backward scheme ln_zdfexp = ', ln_zdfexp
(   83)          WRITE(numout,*) '             number of time step               n_zdfexp = ', n_zdfexp
(   84)       ENDIF
(   85) 
(   86)       ! Parameter & logicals controls
(   87)       ! -----------------------------
(   88)       ! ... vertical mixing
(   89)       ! time stepping scheme (N.B. TKE or KPP schemes => force the use of implicit scheme)
(   90)       IF( ( ln_zdfexp .AND. .NOT.lk_zdftke ) .OR. ( ln_zdfexp .AND. .NOT.lk_zdfkpp ) ) THEN  
(   91)          l_trazdf_exp = .TRUE.           ! use explicit scheme
(   92)          l_trazdf_imp = .FALSE.
(   93)          l_dynzdf_exp = .TRUE.           ! use explicit scheme
(   94)          l_dynzdf_imp = .FALSE.
(   95)       ELSE
(   96)          l_trazdf_exp = .FALSE.          ! use implicit scheme
(   97)          l_trazdf_imp = .TRUE. 
(   98)          l_dynzdf_exp = .FALSE.          ! use implicit scheme
(   99)          l_dynzdf_imp = .TRUE. 
(  100)       ENDIF
(  101)       IF( l_trazdf_iso .OR. l_trazdf_iso_vo ) THEN  
(  102)          l_trazdf_exp = .FALSE.          ! iso-neutral diffusion : 
(  103)          l_trazdf_imp = .FALSE.          ! implicit scheme included in iso-neutral routine
(  104)       ENDIF
(  105)       IF( l_dynldf_iso ) THEN  
(  106)          l_dynzdf_exp = .FALSE.          ! iso-neutral diffusion :
(  107)          l_dynzdf_imp = .FALSE.          ! implicit scheme included in iso-neutral routine
(  108)       ENDIF
(  109) # 115
(  115)       l_dynzdf_imp_tsk = .FALSE.
(  116) # 117
(  117)       IF( lk_esopa  ) THEN
(  118)          l_trazdf_exp = .TRUE.           ! esopa: use all options
(  119)          l_trazdf_imp = .TRUE.
(  120)          l_dynzdf_exp     = .TRUE.           ! esopa: use all options






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 7

(  121)          l_dynzdf_imp     = .TRUE.
(  122)          l_dynzdf_imp_tsk = .TRUE.
(  123)       ENDIF
(  124) 
(  125)       IF(lwp) WRITE(numout,*)
(  126)       IF(lwp) WRITE(numout,*) '          vertical mixing option :'
(  127)       ioptio = 0
(  128)       IF( lk_zdfcst ) THEN
(  129)          IF(lwp) WRITE(numout,*) '             constant eddy diffusion coef.'
(  130)          ioptio = ioptio+1
(  131)       ENDIF
(  132)       IF( lk_zdfric ) THEN
(  133)          IF(lwp) WRITE(numout,*) '             Richardson dependent eddy coef.'
(  134)          ioptio = ioptio+1
(  135)       ENDIF
(  136)       IF( lk_zdftke ) THEN
(  137)          IF(lwp) WRITE(numout,*) '             TKE dependent eddy coef.'
(  138)          ioptio = ioptio+1
(  139)       ENDIF
(  140)       IF( lk_zdfkpp ) THEN
(  141)          IF(lwp) WRITE(numout,*) '             KPP dependent eddy coef.'
(  142)          ioptio = ioptio+1
(  143)       ENDIF
(  144)       IF( ioptio == 0 .OR. ioptio > 1 .AND. .NOT. lk_esopa ) THEN
(  145)           IF(lwp) WRITE(numout,cform_err)
(  146)           IF(lwp) WRITE(numout,*) ' one and only one vertical diffusion option has to be defined '
(  147)           nstop = nstop + 1
(  148)       ENDIF
(  149) 
(  150)       ! ... Convection
(  151)       IF(lwp) WRITE(numout,*)
(  152)       IF(lwp) WRITE(numout,*) '          convection :'
(  153)       ioptio = 0
(  154)       IF( ln_zdfnpc ) THEN
(  155)          IF(lwp) WRITE(numout,*) '             use non penetrative convective scheme'
(  156)          ioptio = ioptio+1
(  157)       ENDIF
(  158)       IF( ln_zdfevd ) THEN
(  159)          IF(lwp) WRITE(numout,*) '             use enhanced vertical dif. scheme'
(  160)          ioptio = ioptio+1
(  161)       ENDIF
(  162)       IF( lk_zdftke ) THEN
(  163)          IF(lwp) WRITE(numout,*) '             use the 1.5 turbulent closure'
(  164)       ENDIF
(  165)       IF( lk_zdfkpp ) THEN
(  166)          IF(lwp) WRITE(numout,*) '             use the KPP closure scheme'
(  167)          IF(lk_mpp) THEN
(  168)             IF(lwp) WRITE(numout,cform_err)
(  169)             IF(lwp) WRITE(numout,*) '             The KPP scheme is not ready to run in MPI'
(  170)          ENDIF
(  171)       ENDIF
(  172)       IF ( ioptio > 1 .AND. .NOT. lk_esopa ) THEN
(  173)           IF(lwp) WRITE(numout,cform_err)
(  174)           IF(lwp) WRITE(numout,*) ' chose between ln_zdfnpc'
(  175)           IF(lwp) WRITE(numout,*) '           and ln_zdfevd'
(  176)           nstop = nstop + 1
(  177)       ENDIF
(  178)       IF( ioptio == 0 .AND. .NOT. lk_zdftke ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 8

(  179)          IF(lwp) WRITE(numout,cform_err)
(  180)          IF(lwp) WRITE(numout,*) ' except for TKE scheme, a convection scheme is'
(  181)          IF(lwp) WRITE(numout,*) ' required: ln_zdfevd or ln_zdfnpc logicals'
(  182)          nstop = nstop + 1
(  183)       ENDIF
(  184) 
(  185)    END SUBROUTINE zdf_init
(  186) 
(  187)    !!======================================================================
(  188) END MODULE zdfini
