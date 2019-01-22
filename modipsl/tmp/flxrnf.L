


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

Filename: flxrnf.F90

(    1) !!DB: 2009.08.31 -- deleted ORCA
(    2) MODULE flxrnf
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE  flxrnf  ***
(    5)    !! Ocean forcing:  runoff
(    6)    !!=====================================================================
(    7) !#if defined key_orca_r05
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_orca_r05'                               ORCA R05 configuration
(   10)    !!----------------------------------------------------------------------
(   11) !#  include "flxrnf_ORCA_R05.h90"
(   12) !#else
(   13)    !!----------------------------------------------------------------------
(   14)    !!   Default option                                     Standard runoffs
(   15)    !!----------------------------------------------------------------------
(   16) 
(   17)    !!----------------------------------------------------------------------
(   18)    !!   flx_rnf      : monthly runoff read in a NetCDF file
(   19)    !!----------------------------------------------------------------------
(   20)    !! * Modules used
(   21)    USE dom_oce         ! ocean space and time domain
(   22)    USE phycst          ! physical constants
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE daymod          ! calendar
(   25)    USE ioipsl          ! NetCDF IPSL library
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Routine accessibility
(   31)    PUBLIC flx_rnf          ! routine call in step module
(   32) 
(   33)    !! * Shared module variables
(   34)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   35)       runoff,           &  !: monthly runoff (kg/m2/s)
(   36)       upsadv,           &  !: mixed adv scheme in straits vicinity (hori.)
(   37)       upsrnfh              !: mixed adv scheme in runoffs vicinity (hori.)
(   38)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
(   39)       upsrnfz              !: mixed adv scheme in runoffs vicinity (vert.)
(   40)    INTEGER, PUBLIC ::   &  !:
(   41)       nrunoff =  0 ,    &  !: runoff option (namelist)
(   42)       nrnf1, nrnf2         !: first and second record used
(   43) 
(   44)    !! * Module variable
(   45)    REAL(wp), DIMENSION(jpi,jpj,2) ::   &  !:
(   46)       rnfdta               !: monthly runoff data array (kg/m2/s)






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 2

(   47)    !!----------------------------------------------------------------------
(   48)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   49)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxrnf.F90,v 1.9 2006/03/09 17:22:02 opalod Exp $ 
(   50)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   51)    !!----------------------------------------------------------------------
(   52) 
(   53) CONTAINS
(   54) 
(   55)    SUBROUTINE flx_rnf( kt )
(   56)       !!----------------------------------------------------------------------
(   57)       !!                  ***  ROUTINE flx_rnf  ***
(   58)       !!       
(   59)       !! ** Purpose :   Introduce a climatological run off forcing
(   60)       !!
(   61)       !! ** Method :
(   62)       !!      Initialze each mouth of river with a monthly climatology 
(   63)       !!      provided from different data.
(   64)       !!     C a u t i o n : upward water flux, runoff is negative
(   65)       !!                     set at the last loop of the routine
(   66)       !!
(   67)       !! ** Action :
(   68)       !!
(   69)       !! References : 
(   70)       !!       J. D. Milliman and R. H. Meade, 1983 : world-wide delivery
(   71)       !!          of river sediment to the oceans, journal of geology vol 91
(   72)       !!          pp 1-21.
(   73)       !!       G. L. Russell and J. R. Miller, 1990 : global river runoff
(   74)       !!          calculated from a global atmospheric general circulation
(   75)       !!          model, journal of hydrology, 117(1990), pp 241-254.
(   76)       !!       F. Van Der Leeden, Troise F. L., Todd D. K. : the water
(   77)       !!          encyclopedia, second edition, lewis publishers.
(   78)       !!       J. W. Weatherly, J. E. Walsh : The effects of precipitation
(   79)       !!          and river runoff in a coupled ice-ocean model of Arctic
(   80)       !!          Climate dynamics 1996 12:785,798
(   81)       !!       Jacobs et al. 1992. J. Glaciol. 38 (130) 375-387.
(   82)       !!
(   83)       !! History :
(   84)       !!        !  94-10  (G.Madec, M. Pontaud, M. Imbard)  Original code
(   85)       !!        !  97-03  (G.Madec)  time dependent version
(   86)       !!        !  98-06  (J.M. Molines)  exact computation of zxy 
(   87)       !!                         for months that are not 30 days
(   88)       !!        !  98-07  (M. Imbard)  ORCA and mpp option
(   89)       !!        !  99-08  (J.P. Boulanger H.L.Ayina)  New rivers and 
(   90)       !!                         values given in m3/s 
(   91)       !!        !  00-04  (G. Madec, K. Roberts) add antarctica ice discharge.
(   92)       !!        !  00-11  (R. Hordoir, E. Durand)  NetCDF FORMAT
(   93)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   94)       !!----------------------------------------------------------------------
(   95)       !! * arguments
(   96)       INTEGER, INTENT( in  ) ::   kt       ! ocean time step
(   97) 
(   98)       !! * Local declarations
(   99) # if ! defined key_coupled
(  100)       INTEGER  ::   ji, jj                 ! dummy loop indices
(  101)       INTEGER ::   &
(  102)          i15 , imois , iman,            &  ! temporary integers
(  103)          idbd, idmeom                      !    "          "
(  104)       REAL(wp) ::   zxy






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 3

(  105) # endif
(  106)       CHARACTER (len=32) ::   &
(  107)          clname                            ! monthly runoff filename
(  108)       INTEGER, PARAMETER :: jpmois = 12
(  109)       INTEGER  ::   ipi, ipj, ipk          ! temporary integers
(  110)       INTEGER  ::   ii0, ii1, ij0, ij1     !    "          "
(  111)       INTEGER, DIMENSION(jpmois) ::     &
(  112)          istep                             ! temporary workspace
(  113)       REAL(wp) ::   zdate0, zdt            ! temporary scalars
(  114)       REAL(wp), DIMENSION(jpk) ::       &
(  115)          zlev                              ! temporary workspace
(  116)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  117)          zlon, zlat,                    &  ! temporary workspace
(  118)          zcoefr                            ! coeff of advection link to runoff
(  119)       !!----------------------------------------------------------------------
(  120)          clname = 'runoff_1m_nomask'       ! monthly runoff filename
(  121)       
(  122)       IF( kt == nit000 ) THEN
(  123) 
(  124)          SELECT CASE ( nrunoff )
(  125) 
(  126)          CASE ( 0 )
(  127)             IF(lwp) WRITE(numout,*)
(  128)             IF(lwp) WRITE(numout,*) 'flx_rnf : No runoff in this simulation (nrunoff=0)'
(  129)             IF(lwp) WRITE(numout,*) '~~~~~~~'
(  130)             
(  131)          CASE ( 1 )
(  132)             IF(lwp) WRITE(numout,*)
(  133)             IF(lwp) WRITE(numout,*) 'flx_rnf : monthly runoff (nrunoff=1)'
(  134)             IF(lwp) WRITE(numout,*) '~~~~~~~'
(  135) 
(  136)          CASE ( 2 )
(  137)             IF(lwp) WRITE(numout,*)
(  138)             IF(lwp) WRITE(numout,*) 'flx_rnf : monthly runoff with upsteam advection'
(  139)             IF(lwp) WRITE(numout,*) '~~~~~~~   in the vicinity of river mouths (nrunoff=2)'
(  140) 
(  141)          CASE DEFAULT
(  142)             IF(lwp) WRITE(numout,cform_err)
(  143)             IF(lwp) WRITE(numout,*) ' Error nrunoff = ', nrunoff, ' /= 0, 1 or 2'
(  144)             nstop = nstop + 1
(  145) 
(  146)          END SELECT
(  147) 
(  148)          ! Set runoffs and upstream coeff to zero
(  149)          runoff (:,:) = 0.e0
(  150)          upsrnfh(:,:) = 0.e0
(  151)          upsrnfz(:)   = 0.e0 
(  152)          upsadv (:,:) = 0.e0
(  153) 
(  154)       ENDIF
(  155) 
(  156) 
(  157)       ! 1. Initialization
(  158)       ! -----------------
(  159) 
(  160)       IF( nrunoff == 1 .OR. nrunoff == 2 ) THEN
(  161) # if ! defined key_coupled
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 4

(  163)          ! year, month, day
(  164)          i15   = nday / 16
(  165)          imois = nmonth + i15 - 1
(  166)          IF( imois == 0 )   imois = jpmois
(  167)          ! Number of days in the month
(  168)          IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
(  169)             idbd = nbiss(imois)
(  170)          ELSEIF( nleapy > 1 ) THEN
(  171)             idbd = nleapy
(  172)          ELSE
(  173)             idbd = nobis(imois)
(  174)          ENDIF
(  175)          ! Number of days between imois, 15 and the end of month
(  176)          idmeom = idbd - 15
(  177) # endif
(  178)          ipi = jpiglo
(  179)          ipj = jpjglo
(  180)          ipk = jpk
(  181)          zdt = rdt
(  182)          
(  183)          ! Open file
(  184) 
(  185)          IF( kt == nit000 ) THEN
(  186)             iman = jpmois
(  187)             CALL flinopen( clname, mig(1), nlci, mjg(1), nlcj,    &
(  188)                &           .false., ipi, ipj, ipk, zlon,        &
(  189)                &           zlat, zlev, iman, istep, zdate0,   &
(  190)                &           zdt, numrnf )
(  191)             !   Title, dimensions and tests
(  192) # if ! defined key_coupled
(  193)             IF( iman /= jpmois ) THEN
(  194)                IF(lwp) WRITE(numout,*)
(  195)                IF(lwp) WRITE(numout,*) 'problem with time coordinates'
(  196)                IF(lwp) WRITE(numout,*) ' iman ', iman, ' jpmois ', jpmois
(  197)                nstop = nstop + 1
(  198)             ENDIF
(  199)             IF(lwp) WRITE(numout,*) iman, istep, zdate0, rdt, numrnf
(  200)             IF(lwp) WRITE(numout,*) 'numrnf=', numrnf
(  201)             IF(lwp) WRITE(numout,*) 'jpmois=', jpmois
(  202)             IF(lwp) WRITE(numout,*) 'zdt=', zdt
(  203) # endif
(  204)             IF(ipi /= jpidta .AND. ipj /= jpjdta .AND. ipk /= 1) THEN
(  205)                IF(lwp)WRITE(numout,*) ' '
(  206)                IF(lwp)WRITE(numout,*) 'problem with dimensions'
(  207)                IF(lwp)WRITE(numout,*) ' ipi ', ipi, ' jpidta ', jpidta
(  208)                IF(lwp)WRITE(numout,*) ' ipj ', ipj, ' jpjdta ', jpjdta
(  209)                IF(lwp)WRITE(numout,*) ' ipk ', ipk, ' =? 1'
(  210)                nstop = nstop + 1
(  211)             ENDIF
(  212)             IF(lwp)WRITE(numout,*) 'ipi=', ipi, ' ipj=', ipj, ' ipk=', ipk
(  213)          ENDIF
(  214)          
(  215) # if ! defined key_coupled
(  216) 
(  217)          ! 2. Read monthly file of runoff
(  218)          ! ------------------------------
(  219) 
(  220)          IF( kt == nit000 .OR. imois /= nrnf1 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 5

(  221) 
(  222)             ! Calendar computation for interpolation
(  223)             !     nrnf1 number of the first array record used in the simulation
(  224)             !     nrnf2 number of the last  array record
(  225) 
(  226)             iman = jpmois
(  227)             nrnf1 = imois
(  228)             nrnf2 = nrnf1 + 1
(  229)             nrnf1 = MOD( nrnf1, iman )
(  230)             IF( nrnf1 == 0 ) nrnf1 = iman
(  231)             nrnf2 = MOD( nrnf2, iman )
(  232)             IF( nrnf2 == 0 ) nrnf2 = iman
(  233)             
(  234)             IF(lwp) THEN
(  235)                WRITE(numout,*)
(  236)                WRITE(numout,*) ' runoff monthly field'
(  237)                WRITE(numout,*) ' --------------------'
(  238)                WRITE(numout,*) ' NetCDF format'
(  239)                WRITE(numout,*)
(  240)                WRITE(numout,*) 'first array record used nrnf1 ',nrnf1
(  241)                WRITE(numout,*) 'last  array record used nrnf2 ',nrnf2
(  242)                WRITE(numout,*)
(  243)             ENDIF
(  244)             
(  245)             ! Read monthly runoff data in kg/m2/s
(  246) !ibug
(  247)             IF( kt == nit000 )   rnfdta(:,:,:) = 0.e0
(  248) !ibug
(  249)             CALL flinget( numrnf, 'sorunoff', jpidta, jpjdta, 1, jpmois   &
(  250)                &        , nrnf1, nrnf1, mig(1), nlci, mjg(1), nlcj, rnfdta(1:nlci,1:nlcj,1) )
(  251)             CALL flinget( numrnf, 'sorunoff', jpidta, jpjdta, 1, jpmois   &
(  252)                &        , nrnf2, nrnf2, mig(1), nlci, mjg(1), nlcj, rnfdta(1:nlci,1:nlcj,2) )
(  253) 
(  254)             IF(lwp) WRITE(numout,*)
(  255)             IF(lwp) WRITE(numout,*) ' read runoff field ok'
(  256)             IF(lwp) WRITE(numout,*)
(  257) 
(  258)          ENDIF
(  259) 
(  260)          ! Linear interpolation and conversion in upward water flux
(  261)          ! C a u t i o n : runoff is negative and in kg/m2/s 
(  262) 
(  263)          zxy = FLOAT( nday + idmeom - idbd * i15 ) / idbd
(  264) !byoung for steady state
(  265)          zxy = 0.0
(  266) 
(  267)          runoff(:,:) = -( ( 1.e0 - zxy ) * rnfdta(:,:,1) + zxy * rnfdta(:,:,2) )
(  268) 
(  269)          
(  270) # endif
(  271) 
(  272)       ENDIF
(  273) 
(  274) 
(  275)       ! 3. Mixed advection scheme 
(  276)       ! -------------------------
(  277) 
(  278)       IF( nrunoff == 2 .AND. kt == nit000 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 6

(  279) 
(  280)          ! Upstream and centered scheme in the vicinity of river mouths
(  281) 
(  282)          !  Creates the array coef that contains the coefficient to affect to
(  283)          !  the upstream scheme. advection scheme will be:
(  284)          !  coefr * upstream + (1- coefr) centered
(  285)          !  coefr must be between 0 and 1.
(  286) !ibug
(  287)          zcoefr(:,:) = 0.e0
(  288) !ibug
(  289) 
(  290)          CALL flinget( numrnf, 'socoefr', jpidta, jpjdta, 1, jpmois, nrnf1,   &
(  291)             &          nrnf1, mig(1), nlci, mjg(1), nlcj, zcoefr(1:nlci,1:nlcj) )
(  292) 
(  293)          IF(lwp) WRITE(numout,*)
(  294)          IF(lwp) WRITE(numout,*) ' read coefr for advection ok'
(  295)          IF(lwp) WRITE(numout,*)
(  296)          
(  297)          upsrnfh(:,:) = zcoefr(:,:)
(  298)          upsrnfz(:)   = 0.e0
(  299)          upsrnfz(1)   = 1.0
(  300)          upsrnfz(2)   = 1.0
(  301)          upsrnfz(3)   = 0.5
(  302)          upsrnfz(4)   = 0.25
(  303)          upsrnfz(5)   = 0.125
(  304)          
(  305) 
(  306)       ENDIF
(  307) 
(  308)       ! Upstream and centered scheme in the vicinity of some straits
(  309) 
(  310)       IF( kt == nit000 ) THEN 
(  311) 
(  312) !         IF( cp_cfg == "orca" ) THEN
(  313) 
(  314)       ENDIF
(  315)      
(  316)       ! 4. Closing all files
(  317)       ! --------------------
(  318) 
(  319)       IF( kt == nitend .AND. nrunoff >= 1 )   CALL flinclo( numrnf )
(  320) 
(  321)    END SUBROUTINE flx_rnf
(  322) 
(  323) !#endif
(  324)    !!======================================================================
(  325) END MODULE flxrnf

















PGF90 (Version     10.5)          09/23/2010  09:49:03      page 7

(    1) # 1 "flxrnf.F90"
(    1) !!DB: 2009.08.31 -- deleted ORCA
(    2) MODULE flxrnf
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE  flxrnf  ***
(    5)    !! Ocean forcing:  runoff
(    6)    !!=====================================================================
(    7) !#if defined key_orca_r05
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_orca_r05'                               ORCA R05 configuration
(   10)    !!----------------------------------------------------------------------
(   11) !#  include "flxrnf_ORCA_R05.h90"
(   12) !#else
(   13)    !!----------------------------------------------------------------------
(   14)    !!   Default option                                     Standard runoffs
(   15)    !!----------------------------------------------------------------------
(   16) 
(   17)    !!----------------------------------------------------------------------
(   18)    !!   flx_rnf      : monthly runoff read in a NetCDF file
(   19)    !!----------------------------------------------------------------------
(   20)    !! * Modules used
(   21)    USE dom_oce         ! ocean space and time domain
(   22)    USE phycst          ! physical constants
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE daymod          ! calendar
(   25)    USE ioipsl          ! NetCDF IPSL library
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Routine accessibility
(   31)    PUBLIC flx_rnf          ! routine call in step module
(   32) 
(   33)    !! * Shared module variables
(   34)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   35)       runoff,           &  !: monthly runoff (kg/m2/s)
(   36)       upsadv,           &  !: mixed adv scheme in straits vicinity (hori.)
(   37)       upsrnfh              !: mixed adv scheme in runoffs vicinity (hori.)
(   38)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
(   39)       upsrnfz              !: mixed adv scheme in runoffs vicinity (vert.)
(   40)    INTEGER, PUBLIC ::   &  !:
(   41)       nrunoff =  0 ,    &  !: runoff option (namelist)
(   42)       nrnf1, nrnf2         !: first and second record used
(   43) 
(   44)    !! * Module variable
(   45)    REAL(wp), DIMENSION(jpi,jpj,2) ::   &  !:
(   46)       rnfdta               !: monthly runoff data array (kg/m2/s)
(   47)    !!----------------------------------------------------------------------
(   48)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   49)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxrnf.F90,v 1.9 2006/03/09 17:22:02 opalod Exp $ 
(   50)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   51)    !!----------------------------------------------------------------------
(   52) 
(   53) CONTAINS
(   54) 
(   55)    SUBROUTINE flx_rnf( kt )
(   56)       !!----------------------------------------------------------------------
(   57)       !!                  ***  ROUTINE flx_rnf  ***






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 8

(   58)       !!       
(   59)       !! ** Purpose :   Introduce a climatological run off forcing
(   60)       !!
(   61)       !! ** Method :
(   62)       !!      Initialze each mouth of river with a monthly climatology 
(   63)       !!      provided from different data.
(   64)       !!     C a u t i o n : upward water flux, runoff is negative
(   65)       !!                     set at the last loop of the routine
(   66)       !!
(   67)       !! ** Action :
(   68)       !!
(   69)       !! References : 
(   70)       !!       J. D. Milliman and R. H. Meade, 1983 : world-wide delivery
(   71)       !!          of river sediment to the oceans, journal of geology vol 91
(   72)       !!          pp 1-21.
(   73)       !!       G. L. Russell and J. R. Miller, 1990 : global river runoff
(   74)       !!          calculated from a global atmospheric general circulation
(   75)       !!          model, journal of hydrology, 117(1990), pp 241-254.
(   76)       !!       F. Van Der Leeden, Troise F. L., Todd D. K. : the water
(   77)       !!          encyclopedia, second edition, lewis publishers.
(   78)       !!       J. W. Weatherly, J. E. Walsh : The effects of precipitation
(   79)       !!          and river runoff in a coupled ice-ocean model of Arctic
(   80)       !!          Climate dynamics 1996 12:785,798
(   81)       !!       Jacobs et al. 1992. J. Glaciol. 38 (130) 375-387.
(   82)       !!
(   83)       !! History :
(   84)       !!        !  94-10  (G.Madec, M. Pontaud, M. Imbard)  Original code
(   85)       !!        !  97-03  (G.Madec)  time dependent version
(   86)       !!        !  98-06  (J.M. Molines)  exact computation of zxy 
(   87)       !!                         for months that are not 30 days
(   88)       !!        !  98-07  (M. Imbard)  ORCA and mpp option
(   89)       !!        !  99-08  (J.P. Boulanger H.L.Ayina)  New rivers and 
(   90)       !!                         values given in m3/s 
(   91)       !!        !  00-04  (G. Madec, K. Roberts) add antarctica ice discharge.
(   92)       !!        !  00-11  (R. Hordoir, E. Durand)  NetCDF FORMAT
(   93)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   94)       !!----------------------------------------------------------------------
(   95)       !! * arguments
(   96)       INTEGER, INTENT( in  ) ::   kt       ! ocean time step
(   97) 
(   98)       !! * Local declarations
(   99) # 100
(  100)       INTEGER  ::   ji, jj                 ! dummy loop indices
(  101)       INTEGER ::   &
(  102)          i15 , imois , iman,            &  ! temporary integers
(  103)          idbd, idmeom                      !    "          "
(  104)       REAL(wp) ::   zxy
(  105) # 106
(  106)       CHARACTER (len=32) ::   &
(  107)          clname                            ! monthly runoff filename
(  108)       INTEGER, PARAMETER :: jpmois = 12
(  109)       INTEGER  ::   ipi, ipj, ipk          ! temporary integers
(  110)       INTEGER  ::   ii0, ii1, ij0, ij1     !    "          "
(  111)       INTEGER, DIMENSION(jpmois) ::     &
(  112)          istep                             ! temporary workspace
(  113)       REAL(wp) ::   zdate0, zdt            ! temporary scalars
(  114)       REAL(wp), DIMENSION(jpk) ::       &
(  115)          zlev                              ! temporary workspace






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 9

(  116)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  117)          zlon, zlat,                    &  ! temporary workspace
(  118)          zcoefr                            ! coeff of advection link to runoff
(  119)       !!----------------------------------------------------------------------
(  120)          clname = 'runoff_1m_nomask'       ! monthly runoff filename
(  121)       
(  122)       IF( kt == nit000 ) THEN
(  123) 
(  124)          SELECT CASE ( nrunoff )
(  125) 
(  126)          CASE ( 0 )
(  127)             IF(lwp) WRITE(numout,*)
(  128)             IF(lwp) WRITE(numout,*) 'flx_rnf : No runoff in this simulation (nrunoff=0)'
(  129)             IF(lwp) WRITE(numout,*) '~~~~~~~'
(  130)             
(  131)          CASE ( 1 )
(  132)             IF(lwp) WRITE(numout,*)
(  133)             IF(lwp) WRITE(numout,*) 'flx_rnf : monthly runoff (nrunoff=1)'
(  134)             IF(lwp) WRITE(numout,*) '~~~~~~~'
(  135) 
(  136)          CASE ( 2 )
(  137)             IF(lwp) WRITE(numout,*)
(  138)             IF(lwp) WRITE(numout,*) 'flx_rnf : monthly runoff with upsteam advection'
(  139)             IF(lwp) WRITE(numout,*) '~~~~~~~   in the vicinity of river mouths (nrunoff=2)'
(  140) 
(  141)          CASE DEFAULT
(  142)             IF(lwp) WRITE(numout,cform_err)
(  143)             IF(lwp) WRITE(numout,*) ' Error nrunoff = ', nrunoff, ' /= 0, 1 or 2'
(  144)             nstop = nstop + 1
(  145) 
(  146)          END SELECT
(  147) 
(  148)          ! Set runoffs and upstream coeff to zero
(  149)          runoff (:,:) = 0.e0
(  150)          upsrnfh(:,:) = 0.e0
(  151)          upsrnfz(:)   = 0.e0 
(  152)          upsadv (:,:) = 0.e0
(  153) 
(  154)       ENDIF
(  155) 
(  156) 
(  157)       ! 1. Initialization
(  158)       ! -----------------
(  159) 
(  160)       IF( nrunoff == 1 .OR. nrunoff == 2 ) THEN
(  161) 
(  162) # 163
(  163)          ! year, month, day
(  164)          i15   = nday / 16
(  165)          imois = nmonth + i15 - 1
(  166)          IF( imois == 0 )   imois = jpmois
(  167)          ! Number of days in the month
(  168)          IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
(  169)             idbd = nbiss(imois)
(  170)          ELSEIF( nleapy > 1 ) THEN
(  171)             idbd = nleapy
(  172)          ELSE
(  173)             idbd = nobis(imois)






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 10

(  174)          ENDIF
(  175)          ! Number of days between imois, 15 and the end of month
(  176)          idmeom = idbd - 15
(  177) # 178
(  178)          ipi = jpiglo
(  179)          ipj = jpjglo
(  180)          ipk = jpk
(  181)          zdt = rdt
(  182)          
(  183)          ! Open file
(  184) 
(  185)          IF( kt == nit000 ) THEN
(  186)             iman = jpmois
(  187)             CALL flinopen( clname, mig(1), nlci, mjg(1), nlcj,    &
(  188)                &           .false., ipi, ipj, ipk, zlon,        &
(  189)                &           zlat, zlev, iman, istep, zdate0,   &
(  190)                &           zdt, numrnf )
(  191)             !   Title, dimensions and tests
(  192) # 193
(  193)             IF( iman /= jpmois ) THEN
(  194)                IF(lwp) WRITE(numout,*)
(  195)                IF(lwp) WRITE(numout,*) 'problem with time coordinates'
(  196)                IF(lwp) WRITE(numout,*) ' iman ', iman, ' jpmois ', jpmois
(  197)                nstop = nstop + 1
(  198)             ENDIF
(  199)             IF(lwp) WRITE(numout,*) iman, istep, zdate0, rdt, numrnf
(  200)             IF(lwp) WRITE(numout,*) 'numrnf=', numrnf
(  201)             IF(lwp) WRITE(numout,*) 'jpmois=', jpmois
(  202)             IF(lwp) WRITE(numout,*) 'zdt=', zdt
(  203) # 204
(  204)             IF(ipi /= jpidta .AND. ipj /= jpjdta .AND. ipk /= 1) THEN
(  205)                IF(lwp)WRITE(numout,*) ' '
(  206)                IF(lwp)WRITE(numout,*) 'problem with dimensions'
(  207)                IF(lwp)WRITE(numout,*) ' ipi ', ipi, ' jpidta ', jpidta
(  208)                IF(lwp)WRITE(numout,*) ' ipj ', ipj, ' jpjdta ', jpjdta
(  209)                IF(lwp)WRITE(numout,*) ' ipk ', ipk, ' =? 1'
(  210)                nstop = nstop + 1
(  211)             ENDIF
(  212)             IF(lwp)WRITE(numout,*) 'ipi=', ipi, ' ipj=', ipj, ' ipk=', ipk
(  213)          ENDIF
(  214)          
(  215) 
(  216) # 217
(  217)          ! 2. Read monthly file of runoff
(  218)          ! ------------------------------
(  219) 
(  220)          IF( kt == nit000 .OR. imois /= nrnf1 ) THEN
(  221) 
(  222)             ! Calendar computation for interpolation
(  223)             !     nrnf1 number of the first array record used in the simulation
(  224)             !     nrnf2 number of the last  array record
(  225) 
(  226)             iman = jpmois
(  227)             nrnf1 = imois
(  228)             nrnf2 = nrnf1 + 1
(  229)             nrnf1 = MOD( nrnf1, iman )
(  230)             IF( nrnf1 == 0 ) nrnf1 = iman
(  231)             nrnf2 = MOD( nrnf2, iman )






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 11

(  232)             IF( nrnf2 == 0 ) nrnf2 = iman
(  233)             
(  234)             IF(lwp) THEN
(  235)                WRITE(numout,*)
(  236)                WRITE(numout,*) ' runoff monthly field'
(  237)                WRITE(numout,*) ' --------------------'
(  238)                WRITE(numout,*) ' NetCDF format'
(  239)                WRITE(numout,*)
(  240)                WRITE(numout,*) 'first array record used nrnf1 ',nrnf1
(  241)                WRITE(numout,*) 'last  array record used nrnf2 ',nrnf2
(  242)                WRITE(numout,*)
(  243)             ENDIF
(  244)             
(  245)             ! Read monthly runoff data in kg/m2/s
(  246) !ibug
(  247)             IF( kt == nit000 )   rnfdta(:,:,:) = 0.e0
(  248) !ibug
(  249)             CALL flinget( numrnf, 'sorunoff', jpidta, jpjdta, 1, jpmois   &
(  250)                &        , nrnf1, nrnf1, mig(1), nlci, mjg(1), nlcj, rnfdta(1:nlci,1:nlcj,1) )
(  251)             CALL flinget( numrnf, 'sorunoff', jpidta, jpjdta, 1, jpmois   &
(  252)                &        , nrnf2, nrnf2, mig(1), nlci, mjg(1), nlcj, rnfdta(1:nlci,1:nlcj,2) )
(  253) 
(  254)             IF(lwp) WRITE(numout,*)
(  255)             IF(lwp) WRITE(numout,*) ' read runoff field ok'
(  256)             IF(lwp) WRITE(numout,*)
(  257) 
(  258)          ENDIF
(  259) 
(  260)          ! Linear interpolation and conversion in upward water flux
(  261)          ! C a u t i o n : runoff is negative and in kg/m2/s 
(  262) 
(  263)          zxy = FLOAT( nday + idmeom - idbd * i15 ) / idbd
(  264) !byoung for steady state
(  265)          zxy = 0.0
(  266) 
(  267)          runoff(:,:) = -( ( 1.e0 - zxy ) * rnfdta(:,:,1) + zxy * rnfdta(:,:,2) )
(  268) 
(  269)          
(  270) 
(  271) # 272
(  272)       ENDIF
(  273) 
(  274) 
(  275)       ! 3. Mixed advection scheme 
(  276)       ! -------------------------
(  277) 
(  278)       IF( nrunoff == 2 .AND. kt == nit000 ) THEN
(  279) 
(  280)          ! Upstream and centered scheme in the vicinity of river mouths
(  281) 
(  282)          !  Creates the array coef that contains the coefficient to affect to
(  283)          !  the upstream scheme. advection scheme will be:
(  284)          !  coefr * upstream + (1- coefr) centered
(  285)          !  coefr must be between 0 and 1.
(  286) !ibug
(  287)          zcoefr(:,:) = 0.e0
(  288) !ibug
(  289) 






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 12

(  290)          CALL flinget( numrnf, 'socoefr', jpidta, jpjdta, 1, jpmois, nrnf1,   &
(  291)             &          nrnf1, mig(1), nlci, mjg(1), nlcj, zcoefr(1:nlci,1:nlcj) )
(  292) 
(  293)          IF(lwp) WRITE(numout,*)
(  294)          IF(lwp) WRITE(numout,*) ' read coefr for advection ok'
(  295)          IF(lwp) WRITE(numout,*)
(  296)          
(  297)          upsrnfh(:,:) = zcoefr(:,:)
(  298)          upsrnfz(:)   = 0.e0
(  299)          upsrnfz(1)   = 1.0
(  300)          upsrnfz(2)   = 1.0
(  301)          upsrnfz(3)   = 0.5
(  302)          upsrnfz(4)   = 0.25
(  303)          upsrnfz(5)   = 0.125
(  304)          
(  305) 
(  306)       ENDIF
(  307) 
(  308)       ! Upstream and centered scheme in the vicinity of some straits
(  309) 
(  310)       IF( kt == nit000 ) THEN 
(  311) 
(  312) !         IF( cp_cfg == "orca" ) THEN
(  313) 
(  314)       ENDIF
(  315)      
(  316)       ! 4. Closing all files
(  317)       ! --------------------
(  318) 
(  319)       IF( kt == nitend .AND. nrunoff >= 1 )   CALL flinclo( numrnf )
(  320) 
(  321)    END SUBROUTINE flx_rnf
(  322) 
(  323) !#endif
(  324)    !!======================================================================
(  325) END MODULE flxrnf
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (flxrnf.F90: 263)
