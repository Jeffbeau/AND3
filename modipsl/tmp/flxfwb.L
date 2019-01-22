


PGF90 (Version     10.5)          09/23/2010  09:49:49      page 1

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

Filename: flxfwb.F90

(    1) MODULE flxfwb
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  flxfwb  ***
(    4)    !! Ocean fluxes   : domain averaged freshwater budget
(    5)    !!======================================================================
(    6) #if ! defined key_dynspg_rl
(    7)    !!----------------------------------------------------------------------
(    8)    !!   NOT 'key_dynspg_rl'                        Free surface formulation
(    9)    !!----------------------------------------------------------------------
(   10)    !!   flx_fwb      : freshwater budget for global ocean configurations
(   11)    !!                  in free surface and forced mode
(   12)    !!   flx_fwb_init : freshwater budget for global ocean configurations
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE cpl_oce         ! coupled atmosphere/ocean
(   18)    USE phycst          ! physical constants
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lib_mpp         ! distribued memory computing library
(   21)    USE flxrnf          ! ocean runoffs
(   22)    USE ocesbc          ! ocean surface boudaries conditions
(   23)    USE blk_oce
(   24)    USE flxblk          ! bulk formulea
(   25)    USE daymod          ! calendar
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Routine accessibility
(   31)    PUBLIC flx_fwb      ! routine called by step
(   32)    PUBLIC flx_fwb_init ! routine called by opa
(   33) 
(   34)    !! * Share module variables
(   35)    LOGICAL, PUBLIC ::    & !!: * namelist *
(   36)       ln_fwb = .TRUE.       !: Flag to activate the fwb computation
(   37)    REAL(wp), PUBLIC ::   &  !:
(   38)       a_fwb_b  ,         &  !: annual domain averaged freshwater budget
(   39)       a_fwb                 !: for 2 year before (_b) and before year.
(   40) 
(   41)    !! * Module variables
(   42)    REAL(wp) ::   &
(   43)       a_emp   ,  & ! domain averaged evaporation minus precipitation
(   44)       a_precip,  & ! domain averaged precipitation
(   45)       a_rnf   ,  & ! domain averaged runoff
(   46)       a_sshb  ,  & ! domain averaged sea surface heigh at nit000






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 2

(   47)       a_sshend,  & ! domain averaged sea surface heigh at nitend
(   48)       a_sal000,  & ! domain averaged ocean salinity at nit000 (before time step)
(   49)       a_salend,  & ! domain averaged ocean salinity at nitend (now time step)
(   50)       a_aminus,  & ! 
(   51)       a_aplus
(   52)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   53)          e1e2_i    ! area of the interior domain (e1t*e2t*tmask_i)
(   54) 
(   55)    !! * Substitutions
(   56) #  include "domzgr_substitute.h90"
(   57) #  include "vectopt_loop_substitute.h90"
(   58)    !!----------------------------------------------------------------------
(   59)    !!----------------------------------------------------------------------
(   60)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   61)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxfwb.F90,v 1.6 2005/12/21 10:46:40 opalod Exp $ 
(   62)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   63)    !!----------------------------------------------------------------------
(   64) CONTAINS
(   65) 
(   66)    SUBROUTINE flx_fwb( kt )
(   67)       !!---------------------------------------------------------------------
(   68)       !!                  ***  ROUTINE flx_fwb  ***
(   69)       !!
(   70)       !! ** Purpose :
(   71)       !!
(   72)       !! ** Method  :
(   73)       !!	
(   74)       !! History :
(   75)       !!   8.2  !  01-02  (E. Durand)  Original code
(   76)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   77)       !!----------------------------------------------------------------------
(   78)       !! * Arguments
(   79)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   80) 
(   81)       !! * Local declarations
(   82)       INTEGER  ::   ji, jj, jk        ! dummy loop indices
(   83)       INTEGER  ::   inum = 11         ! temporary logical unit
(   84)       INTEGER  ::   ikty              ! 
(   85)       REAL(wp) ::   &
(   86)          zarea, zvol, zwei,       &  ! temporary scalars
(   87)          zsm0, zempnew                !    "         "
(   88)       !!----------------------------------------------------------------------
(   89) 
(   90)       ! Mean global salinity
(   91)       zsm0 = 34.72654
(   92) 
(   93)       ! To compute emp mean value mean emp
(   94) 
(   95)       IF( kt == nit000 ) THEN
(   96)          IF(lwp) THEN
(   97)             WRITE(numout,*)
(   98)             WRITE(numout,*) 'flx_fwb : FreshWater Budget correction'
(   99)             WRITE(numout,*) '~~~~~~~'
(  100)          ENDIF
(  101) 
(  102)          a_emp    = 0.e0
(  103)          a_precip = 0.e0
(  104)          a_rnf    = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 3

(  105)          a_sshb   = 0.e0   ! averaged sea surface heigh at nit000 
(  106)          a_sal000 = 0.e0   ! averaged ocean salinity at nit000
(  107)          a_aminus = 0.e0
(  108)          a_aplus  = 0.e0
(  109)          
(  110)          e1e2_i(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:)
(  111)          a_sshb = SUM( e1e2_i(:,:) * sshn(:,:) )
(  112)          IF( lk_mpp )   CALL  mpp_sum( a_sshb   )   ! sum over all the global domain
(  113)          DO jk = 1, jpkm1
(  114)             DO jj = 2, jpjm1
(  115)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  116)                   zwei  = e1e2_i(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)
(  117)                   a_sal000 = a_sal000 + ( sb(ji,jj,jk) - zsm0 ) * zwei
(  118)                END DO
(  119)             END DO
(  120)          END DO
(  121)          IF( lk_mpp )   CALL  mpp_sum( a_sal000 )   ! sum over the global domain
(  122) 
(  123)       ENDIF
(  124)       
(  125)       ! cumulate surface freshwater budget at each time-step
(  126)       ! --------------------------------------====----------
(  127)       a_emp    = SUM( e1e2_i(:,:) * emp   (:,:) )
(  128)       IF( lk_mpp )   CALL  mpp_sum( a_emp    )   ! sum over the global domain
(  129) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  130)       a_precip = SUM( e1e2_i(:,:) * watm  (:,:) )
(  131)       IF( lk_mpp )   CALL  mpp_sum( a_precip )   ! sum over the global domain
(  132) #endif
(  133)       a_rnf    = SUM( e1e2_i(:,:) * runoff(:,:) )
(  134)       IF( lk_mpp )   CALL  mpp_sum( a_rnf    )   ! sum over the global domain
(  135) 
(  136)       IF( aminus /= 0.e0 ) a_aminus = a_aminus + ( MIN( aplus, aminus ) / aminus )
(  137)       IF( aplus  /= 0.e0 ) a_aplus  = a_aplus  + ( MIN( aplus, aminus ) / aplus  )
(  138) 
(  139) 
(  140)       ! Update empold if new year start
(  141)       ikty = 365 * 86400 / rdttra(1)    !!bug  use of 365 days leap year or 360d year !!!!!!!
(  142)       IF( MOD( kt, ikty ) == 0 ) THEN
(  143)          zarea    = SUM( e1e2_i(:,:)            )
(  144)          IF( lk_mpp )   CALL  mpp_sum( zarea    )   ! sum over the global domain
(  145)          a_fwb_b = a_fwb
(  146)          a_fwb   = SUM( e1e2_i(:,:) * sshn(:,:) )
(  147)          IF( lk_mpp )   CALL  mpp_sum( a_fwb    )   ! sum over the global domain
(  148) 
(  149)          a_fwb   = a_fwb * 1.e+3 / ( zarea * 86400. * 365. )    ! convert in Kg/m3/s = mm/s
(  150)          !                                                      !!bug 365d year 
(  151)          empold =  a_fwb                 ! current year freshwater budget correction
(  152)          !                               ! estimate from the previous year budget
(  153)       ENDIF
(  154) 
(  155) 
(  156)       IF( kt == nitend ) THEN
(  157)          zvol  = 0.e0
(  158)          zempnew = 0.e0
(  159)          ! Mean sea level at nitend
(  160)          a_sshend = SUM( e1e2_i(:,:) * sshn(:,:) )
(  161)          zarea    = SUM( e1e2_i(:,:)             )
(  162)          






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 4

(  163)          a_salend = 0.e0
(  164)          DO jk = 1, jpkm1   
(  165)             DO jj = 2, jpjm1
(  166)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  167)                   zwei  = e1e2_i(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)
(  168)                   a_salend = a_salend + ( sn(ji,jj,jk) - zsm0 ) * zwei
(  169)                   zvol  = zvol  + zwei
(  170)                END DO
(  171)             END DO
(  172)          END DO
(  173)          IF( lk_mpp ) THEN
(  174)             CALL  mpp_sum( zarea    )   ! sums over all the global domain
(  175)             CALL  mpp_sum( a_sshend )        
(  176)             CALL  mpp_sum( a_salend )        
(  177)          ENDIF
(  178) 
(  179)          a_aminus = a_aminus / ( nitend - nit000 + 1 )
(  180)          a_aplus  = a_aplus  / ( nitend - nit000 + 1 )
(  181) 
(  182)          ! Conversion in m3
(  183)          a_emp    = a_emp    * rdttra(1) * 1.e-3 
(  184)          a_precip = a_precip * rdttra(1) * 1.e-3 / rday
(  185)          a_rnf    = a_rnf    * rdttra(1) * 1.e-3
(  186)          
(  187)       ENDIF
(  188) 
(  189) 
(  190)       ! Ecriture des diagnostiques 
(  191)       ! --------------------------
(  192) 
(  193)       IF( kt == nitend ) THEN
(  194) 
(  195)          OPEN( inum, FILE='EMPave.dat' )
(  196)          WRITE(inum, "(24X,I8,2ES24.16)" ) nyear, a_fwb_b, a_fwb
(  197)          WRITE(inum,*)
(  198)          WRITE(inum,*)    'Net freshwater budget '
(  199)          WRITE(inum,9010) '  emp    = ', a_emp   , ' m3 =', a_emp   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  200)          WRITE(inum,9010) '  precip = ', a_precip, ' m3 =', a_precip/((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  201)          WRITE(inum,9010) '  a_rnf  = ', a_rnf   , ' m3 =', a_rnf   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  202)          WRITE(inum,*)
(  203)          WRITE(inum,9010) '  zarea =',zarea
(  204)          WRITE(inum,9010) '  zvol  =',zvol
(  205)          WRITE(inum,*)
(  206)          WRITE(inum,*)    'Mean sea level : '
(  207)          WRITE(inum,9010) '  at nit000 = ',a_sshb           ,' m3 '
(  208)          WRITE(inum,9010) '  at nitend = ',a_sshend         ,' m3 '
(  209)          WRITE(inum,9010) '  diff      = ',(a_sshend-a_sshb),' m3 =', (a_sshend-a_sshb)/((nitend-nit000+1)*rdt) * 1.e-6,' Sv
(  210)          WRITE(inum,9020) '  mean sea level elevation    =', a_sshend/zarea,' m'
(  211)          WRITE(inum,*)
(  212)          WRITE(inum,*)    'Anomaly of salinity content : '
(  213)          WRITE(inum,9010) '  at nit000 = ', a_sal000           ,' psu.m3 '
(  214)          WRITE(inum,9010) '  at nitend = ', a_salend           ,' psu.m3 '
(  215)          WRITE(inum,9010) '  diff      = ', a_salend - a_sal000,' psu.m3'
(  216)          WRITE(inum,*)
(  217)          WRITE(inum,*)    'Mean salinity : '
(  218)          WRITE(inum,9020) '  at nit000 = ',  a_sal000/zvol+zsm0     ,' psu '
(  219)          WRITE(inum,9020) '  at nitend = ',  a_salend/zvol+zsm0     ,' psu '
(  220)          WRITE(inum,9020) '  diff      = ', (a_salend-a_sal000)/zvol,' psu'






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 5

(  221)          WRITE(inum,9020) '  S-SLevitus= ',  a_salend          /zvol,' psu'
(  222)          WRITE(inum,*)
(  223)          WRITE(inum,*)    'Coeff : '
(  224)          WRITE(inum,9030) '  Alpha+   =  ', a_aplus
(  225)          WRITE(inum,9030) '  Alpha-   =  ', a_aminus
(  226)          WRITE(inum,*)
(  227)       ENDIF
(  228) 
(  229)  9006 FORMAT(1X,A,ES24.16)
(  230)  9010 FORMAT(1X,A,ES12.5,A,F10.5,A)
(  231)  9020 FORMAT(1X,A,F10.5,A)
(  232)  9030 FORMAT(1X,A,F8.2,A)
(  233) 
(  234)    END SUBROUTINE flx_fwb
(  235) 
(  236) 
(  237)    SUBROUTINE flx_fwb_init
(  238)       !!---------------------------------------------------------------------
(  239)       !!                ***  ROUTINE flx_fwb_init  ***
(  240)       !!
(  241)       !! ** Purpose :
(  242)       !!
(  243)       !! ** Method  :
(  244)       !!  
(  245)       !! History :
(  246)       !!   9.0  !  03-09  (G. Madec)  Original code
(  247)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  248)       !!----------------------------------------------------------------------
(  249)       !! * Local declarations
(  250)       LOGICAL ::   llbon
(  251)       CHARACTER (len=32) ::   &
(  252)          clname = 'EMPave_old.dat'
(  253)       INTEGER ::   inum = 11         ! temporary logical unit
(  254)       INTEGER ::   iyear
(  255) 
(  256)       NAMELIST/namfwb/ ln_fwb 
(  257)       !!----------------------------------------------------------------------
(  258)          
(  259)       ! Read Namelist namfwb : freshWater Budget correction
(  260)       ! --------------------
(  261)       REWIND( numnam )
(  262)       READ  ( numnam, namfwb )            
(  263)       
(  264)       ! Parameter control and print
(  265)       ! ---------------------------
(  266)       ! Control print
(  267)       IF(lwp) THEN
(  268)          WRITE(numout,*)
(  269)          WRITE(numout,*) 'flx_fwb_init : FreshWater Budget correction'
(  270)          WRITE(numout,*) '~~~~~~~~~~~~'
(  271)          WRITE(numout,*) '               Namelist namfwb : set fwb parameters'
(  272)          WRITE(numout,*) '                  use or not fwb correction      ln_fwb   = ', ln_fwb
(  273)       ENDIF
(  274)       ! Option consistency
(  275) #if defined key_dynspg_rl
(  276)       IF(lwp) WRITE '               Rigid-lid option, fwb correction is useless, but valid'
(  277) #endif
(  278)       IF( lk_cpl ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 6

(  279)          IF(lwp) WRITE(numout,*) '               Coupled option, fwb correction is a flux correction ! '
(  280)          IF(lwp) WRITE(numout,*) '               ln_fwb = .FALSE. is recommanded'
(  281)       ENDIF
(  282) 
(  283)       !                                        ! ==============================
(  284)       IF( ln_fwb ) THEN                        !  Freshwater budget correction 
(  285)          !                                     ! ==============================
(  286)          ! Read the corrective factor on precipitations (empold)
(  287)          INQUIRE( FILE=clname, EXIST=llbon )
(  288)          IF( llbon ) THEN
(  289)             OPEN ( inum, FILE=clname)
(  290)             READ ( inum, "(24X,I8,2ES24.16)" ) iyear, a_fwb_b, a_fwb
(  291)             CLOSE( inum )
(  292)             empold = a_fwb                  ! current year freshwater budget correction
(  293)             !                               ! estimate from the previous year budget
(  294)             IF(lwp)WRITE(numout,*)
(  295)             IF(lwp)WRITE(numout,*)'flx_fwb_init : year = ',iyear  , ' freshwater budget correction = ', empold
(  296)             IF(lwp)WRITE(numout,*)'               year = ',iyear-1, ' freshwater budget read       = ', a_fwb
(  297)             IF(lwp)WRITE(numout,*)'               year = ',iyear-2, ' freshwater budget read       = ', a_fwb_b
(  298)          ELSE
(  299)             IF(lwp) WRITE(numout,*)
(  300)             IF(lwp) WRITE(numout,*)'flx_fwb_init : unable to read the file', clname
(  301)             nstop = nstop + 1
(  302)          ENDIF
(  303)          !                                    ! ==============================
(  304)       ELSE                                    !      NO  budget correction 
(  305)          !                                    ! ==============================
(  306)          IF(lwp) WRITE(numout,*)
(  307)          IF(lwp) WRITE(numout,*)'flx_fwb_init : NO freshwater budget correction'
(  308)          empold = 0.e0
(  309)       ENDIF
(  310) 
(  311)    END SUBROUTINE flx_fwb_init
(  312) 
(  313) #else
(  314)    !!----------------------------------------------------------------------
(  315)    !!   Default case :                       
(  316)    !!----------------------------------------------------------------------
(  317)    USE in_out_manager
(  318)    LOGICAL, PUBLIC ::   ln_fwb = .FALSE.   !: no fwb forced
(  319) CONTAINS
(  320)    SUBROUTINE flx_fwb( kt )                ! dummy routine
(  321)       if(lwp) WRITE(numout,*) 'flx_fwb: You should not have seen this print! error?', kt
(  322)    END SUBROUTINE flx_fwb
(  323)    SUBROUTINE flx_fwb_init
(  324)    END SUBROUTINE flx_fwb_init
(  325) #endif
(  326)    !!======================================================================
(  327) END MODULE flxfwb















PGF90 (Version     10.5)          09/23/2010  09:49:49      page 7

(    1) # 1 "flxfwb.F90"
(    1) MODULE flxfwb
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  flxfwb  ***
(    4)    !! Ocean fluxes   : domain averaged freshwater budget
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   NOT 'key_dynspg_rl'                        Free surface formulation
(    9)    !!----------------------------------------------------------------------
(   10)    !!   flx_fwb      : freshwater budget for global ocean configurations
(   11)    !!                  in free surface and forced mode
(   12)    !!   flx_fwb_init : freshwater budget for global ocean configurations
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE cpl_oce         ! coupled atmosphere/ocean
(   18)    USE phycst          ! physical constants
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lib_mpp         ! distribued memory computing library
(   21)    USE flxrnf          ! ocean runoffs
(   22)    USE ocesbc          ! ocean surface boudaries conditions
(   23)    USE blk_oce
(   24)    USE flxblk          ! bulk formulea
(   25)    USE daymod          ! calendar
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Routine accessibility
(   31)    PUBLIC flx_fwb      ! routine called by step
(   32)    PUBLIC flx_fwb_init ! routine called by opa
(   33) 
(   34)    !! * Share module variables
(   35)    LOGICAL, PUBLIC ::    & !!: * namelist *
(   36)       ln_fwb = .TRUE.       !: Flag to activate the fwb computation
(   37)    REAL(wp), PUBLIC ::   &  !:
(   38)       a_fwb_b  ,         &  !: annual domain averaged freshwater budget
(   39)       a_fwb                 !: for 2 year before (_b) and before year.
(   40) 
(   41)    !! * Module variables
(   42)    REAL(wp) ::   &
(   43)       a_emp   ,  & ! domain averaged evaporation minus precipitation
(   44)       a_precip,  & ! domain averaged precipitation
(   45)       a_rnf   ,  & ! domain averaged runoff
(   46)       a_sshb  ,  & ! domain averaged sea surface heigh at nit000
(   47)       a_sshend,  & ! domain averaged sea surface heigh at nitend
(   48)       a_sal000,  & ! domain averaged ocean salinity at nit000 (before time step)
(   49)       a_salend,  & ! domain averaged ocean salinity at nitend (now time step)
(   50)       a_aminus,  & ! 
(   51)       a_aplus
(   52)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   53)          e1e2_i    ! area of the interior domain (e1t*e2t*tmask_i)
(   54) 
(   55)    !! * Substitutions
(   56) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 8

(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 57 "flxfwb.F90"
(   57) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 58 "flxfwb.F90"
(   58) # 58
(   58)    !!----------------------------------------------------------------------
(   59)    !!----------------------------------------------------------------------
(   60)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   61)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxfwb.F90,v 1.6 2005/12/21 10:46:40 opalod Exp $ 
(   62)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   63)    !!----------------------------------------------------------------------
(   64) CONTAINS
(   65) 
(   66)    SUBROUTINE flx_fwb( kt )
(   67)       !!---------------------------------------------------------------------
(   68)       !!                  ***  ROUTINE flx_fwb  ***
(   69)       !!
(   70)       !! ** Purpose :
(   71)       !!
(   72)       !! ** Method  :
(   73)       !!	
(   74)       !! History :
(   75)       !!   8.2  !  01-02  (E. Durand)  Original code
(   76)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   77)       !!----------------------------------------------------------------------
(   78)       !! * Arguments
(   79)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   80) 
(   81)       !! * Local declarations
(   82)       INTEGER  ::   ji, jj, jk        ! dummy loop indices
(   83)       INTEGER  ::   inum = 11         ! temporary logical unit
(   84)       INTEGER  ::   ikty              ! 
(   85)       REAL(wp) ::   &
(   86)          zarea, zvol, zwei,       &  ! temporary scalars
(   87)          zsm0, zempnew                !    "         "






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 9

(   88)       !!----------------------------------------------------------------------
(   89) 
(   90)       ! Mean global salinity
(   91)       zsm0 = 34.72654
(   92) 
(   93)       ! To compute emp mean value mean emp
(   94) 
(   95)       IF( kt == nit000 ) THEN
(   96)          IF(lwp) THEN
(   97)             WRITE(numout,*)
(   98)             WRITE(numout,*) 'flx_fwb : FreshWater Budget correction'
(   99)             WRITE(numout,*) '~~~~~~~'
(  100)          ENDIF
(  101) 
(  102)          a_emp    = 0.e0
(  103)          a_precip = 0.e0
(  104)          a_rnf    = 0.e0
(  105)          a_sshb   = 0.e0   ! averaged sea surface heigh at nit000 
(  106)          a_sal000 = 0.e0   ! averaged ocean salinity at nit000
(  107)          a_aminus = 0.e0
(  108)          a_aplus  = 0.e0
(  109)          
(  110)          e1e2_i(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:)
(  111)          a_sshb = SUM( e1e2_i(:,:) * sshn(:,:) )
(  112)          IF( lk_mpp )   CALL  mpp_sum( a_sshb   )   ! sum over all the global domain
(  113)          DO jk = 1, jpkm1
(  114)             DO jj = 2, jpjm1
(  115)                DO ji = 2, jpim1   ! vector opt.
(  116)                   zwei  = e1e2_i(ji,jj) * e3t(jk) * tmask(ji,jj,jk)
(  117)                   a_sal000 = a_sal000 + ( sb(ji,jj,jk) - zsm0 ) * zwei
(  118)                END DO
(  119)             END DO
(  120)          END DO
(  121)          IF( lk_mpp )   CALL  mpp_sum( a_sal000 )   ! sum over the global domain
(  122) 
(  123)       ENDIF
(  124)       
(  125)       ! cumulate surface freshwater budget at each time-step
(  126)       ! --------------------------------------====----------
(  127)       a_emp    = SUM( e1e2_i(:,:) * emp   (:,:) )
(  128)       IF( lk_mpp )   CALL  mpp_sum( a_emp    )   ! sum over the global domain
(  129) # 130
(  130)       a_precip = SUM( e1e2_i(:,:) * watm  (:,:) )
(  131)       IF( lk_mpp )   CALL  mpp_sum( a_precip )   ! sum over the global domain
(  132) # 133
(  133)       a_rnf    = SUM( e1e2_i(:,:) * runoff(:,:) )
(  134)       IF( lk_mpp )   CALL  mpp_sum( a_rnf    )   ! sum over the global domain
(  135) 
(  136)       IF( aminus /= 0.e0 ) a_aminus = a_aminus + ( MIN( aplus, aminus ) / aminus )
(  137)       IF( aplus  /= 0.e0 ) a_aplus  = a_aplus  + ( MIN( aplus, aminus ) / aplus  )
(  138) 
(  139) 
(  140)       ! Update empold if new year start
(  141)       ikty = 365 * 86400 / rdttra(1)    !!bug  use of 365 days leap year or 360d year !!!!!!!
(  142)       IF( MOD( kt, ikty ) == 0 ) THEN
(  143)          zarea    = SUM( e1e2_i(:,:)            )
(  144)          IF( lk_mpp )   CALL  mpp_sum( zarea    )   ! sum over the global domain
(  145)          a_fwb_b = a_fwb






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 10

(  146)          a_fwb   = SUM( e1e2_i(:,:) * sshn(:,:) )
(  147)          IF( lk_mpp )   CALL  mpp_sum( a_fwb    )   ! sum over the global domain
(  148) 
(  149)          a_fwb   = a_fwb * 1.e+3 / ( zarea * 86400. * 365. )    ! convert in Kg/m3/s = mm/s
(  150)          !                                                      !!bug 365d year 
(  151)          empold =  a_fwb                 ! current year freshwater budget correction
(  152)          !                               ! estimate from the previous year budget
(  153)       ENDIF
(  154) 
(  155) 
(  156)       IF( kt == nitend ) THEN
(  157)          zvol  = 0.e0
(  158)          zempnew = 0.e0
(  159)          ! Mean sea level at nitend
(  160)          a_sshend = SUM( e1e2_i(:,:) * sshn(:,:) )
(  161)          zarea    = SUM( e1e2_i(:,:)             )
(  162)          
(  163)          a_salend = 0.e0
(  164)          DO jk = 1, jpkm1   
(  165)             DO jj = 2, jpjm1
(  166)                DO ji = 2, jpim1   ! vector opt.
(  167)                   zwei  = e1e2_i(ji,jj) * e3t(jk) * tmask(ji,jj,jk)
(  168)                   a_salend = a_salend + ( sn(ji,jj,jk) - zsm0 ) * zwei
(  169)                   zvol  = zvol  + zwei
(  170)                END DO
(  171)             END DO
(  172)          END DO
(  173)          IF( lk_mpp ) THEN
(  174)             CALL  mpp_sum( zarea    )   ! sums over all the global domain
(  175)             CALL  mpp_sum( a_sshend )        
(  176)             CALL  mpp_sum( a_salend )        
(  177)          ENDIF
(  178) 
(  179)          a_aminus = a_aminus / ( nitend - nit000 + 1 )
(  180)          a_aplus  = a_aplus  / ( nitend - nit000 + 1 )
(  181) 
(  182)          ! Conversion in m3
(  183)          a_emp    = a_emp    * rdttra(1) * 1.e-3 
(  184)          a_precip = a_precip * rdttra(1) * 1.e-3 / rday
(  185)          a_rnf    = a_rnf    * rdttra(1) * 1.e-3
(  186)          
(  187)       ENDIF
(  188) 
(  189) 
(  190)       ! Ecriture des diagnostiques 
(  191)       ! --------------------------
(  192) 
(  193)       IF( kt == nitend ) THEN
(  194) 
(  195)          OPEN( inum, FILE='EMPave.dat' )
(  196)          WRITE(inum, "(24X,I8,2ES24.16)" ) nyear, a_fwb_b, a_fwb
(  197)          WRITE(inum,*)
(  198)          WRITE(inum,*)    'Net freshwater budget '
(  199)          WRITE(inum,9010) '  emp    = ', a_emp   , ' m3 =', a_emp   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  200)          WRITE(inum,9010) '  precip = ', a_precip, ' m3 =', a_precip/((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  201)          WRITE(inum,9010) '  a_rnf  = ', a_rnf   , ' m3 =', a_rnf   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  202)          WRITE(inum,*)
(  203)          WRITE(inum,9010) '  zarea =',zarea






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 11

(  204)          WRITE(inum,9010) '  zvol  =',zvol
(  205)          WRITE(inum,*)
(  206)          WRITE(inum,*)    'Mean sea level : '
(  207)          WRITE(inum,9010) '  at nit000 = ',a_sshb           ,' m3 '
(  208)          WRITE(inum,9010) '  at nitend = ',a_sshend         ,' m3 '
(  209)          WRITE(inum,9010) '  diff      = ',(a_sshend-a_sshb),' m3 =', (a_sshend-a_sshb)/((nitend-nit000+1)*rdt) * 1.e-6,' Sv'
(  210)          WRITE(inum,9020) '  mean sea level elevation    =', a_sshend/zarea,' m'
(  211)          WRITE(inum,*)
(  212)          WRITE(inum,*)    'Anomaly of salinity content : '
(  213)          WRITE(inum,9010) '  at nit000 = ', a_sal000           ,' psu.m3 '
(  214)          WRITE(inum,9010) '  at nitend = ', a_salend           ,' psu.m3 '
(  215)          WRITE(inum,9010) '  diff      = ', a_salend - a_sal000,' psu.m3'
(  216)          WRITE(inum,*)
(  217)          WRITE(inum,*)    'Mean salinity : '
(  218)          WRITE(inum,9020) '  at nit000 = ',  a_sal000/zvol+zsm0     ,' psu '
(  219)          WRITE(inum,9020) '  at nitend = ',  a_salend/zvol+zsm0     ,' psu '
(  220)          WRITE(inum,9020) '  diff      = ', (a_salend-a_sal000)/zvol,' psu'
(  221)          WRITE(inum,9020) '  S-SLevitus= ',  a_salend          /zvol,' psu'
(  222)          WRITE(inum,*)
(  223)          WRITE(inum,*)    'Coeff : '
(  224)          WRITE(inum,9030) '  Alpha+   =  ', a_aplus
(  225)          WRITE(inum,9030) '  Alpha-   =  ', a_aminus
(  226)          WRITE(inum,*)
(  227)       ENDIF
(  228) 
(  229)  9006 FORMAT(1X,A,ES24.16)
(  230)  9010 FORMAT(1X,A,ES12.5,A,F10.5,A)
(  231)  9020 FORMAT(1X,A,F10.5,A)
(  232)  9030 FORMAT(1X,A,F8.2,A)
(  233) 
(  234)    END SUBROUTINE flx_fwb
(  235) 
(  236) 
(  237)    SUBROUTINE flx_fwb_init
(  238)       !!---------------------------------------------------------------------
(  239)       !!                ***  ROUTINE flx_fwb_init  ***
(  240)       !!
(  241)       !! ** Purpose :
(  242)       !!
(  243)       !! ** Method  :
(  244)       !!  
(  245)       !! History :
(  246)       !!   9.0  !  03-09  (G. Madec)  Original code
(  247)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  248)       !!----------------------------------------------------------------------
(  249)       !! * Local declarations
(  250)       LOGICAL ::   llbon
(  251)       CHARACTER (len=32) ::   &
(  252)          clname = 'EMPave_old.dat'
(  253)       INTEGER ::   inum = 11         ! temporary logical unit
(  254)       INTEGER ::   iyear
(  255) 
(  256)       NAMELIST/namfwb/ ln_fwb 
(  257)       !!----------------------------------------------------------------------
(  258)          
(  259)       ! Read Namelist namfwb : freshWater Budget correction
(  260)       ! --------------------
(  261)       REWIND( numnam )






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 12

(  262)       READ  ( numnam, namfwb )            
(  263)       
(  264)       ! Parameter control and print
(  265)       ! ---------------------------
(  266)       ! Control print
(  267)       IF(lwp) THEN
(  268)          WRITE(numout,*)
(  269)          WRITE(numout,*) 'flx_fwb_init : FreshWater Budget correction'
(  270)          WRITE(numout,*) '~~~~~~~~~~~~'
(  271)          WRITE(numout,*) '               Namelist namfwb : set fwb parameters'
(  272)          WRITE(numout,*) '                  use or not fwb correction      ln_fwb   = ', ln_fwb
(  273)       ENDIF
(  274)       ! Option consistency
(  275) # 278
(  278)       IF( lk_cpl ) THEN
(  279)          IF(lwp) WRITE(numout,*) '               Coupled option, fwb correction is a flux correction ! '
(  280)          IF(lwp) WRITE(numout,*) '               ln_fwb = .FALSE. is recommanded'
(  281)       ENDIF
(  282) 
(  283)       !                                        ! ==============================
(  284)       IF( ln_fwb ) THEN                        !  Freshwater budget correction 
(  285)          !                                     ! ==============================
(  286)          ! Read the corrective factor on precipitations (empold)
(  287)          INQUIRE( FILE=clname, EXIST=llbon )
(  288)          IF( llbon ) THEN
(  289)             OPEN ( inum, FILE=clname)
(  290)             READ ( inum, "(24X,I8,2ES24.16)" ) iyear, a_fwb_b, a_fwb
(  291)             CLOSE( inum )
(  292)             empold = a_fwb                  ! current year freshwater budget correction
(  293)             !                               ! estimate from the previous year budget
(  294)             IF(lwp)WRITE(numout,*)
(  295)             IF(lwp)WRITE(numout,*)'flx_fwb_init : year = ',iyear  , ' freshwater budget correction = ', empold
(  296)             IF(lwp)WRITE(numout,*)'               year = ',iyear-1, ' freshwater budget read       = ', a_fwb
(  297)             IF(lwp)WRITE(numout,*)'               year = ',iyear-2, ' freshwater budget read       = ', a_fwb_b
(  298)          ELSE
(  299)             IF(lwp) WRITE(numout,*)
(  300)             IF(lwp) WRITE(numout,*)'flx_fwb_init : unable to read the file', clname
(  301)             nstop = nstop + 1
(  302)          ENDIF
(  303)          !                                    ! ==============================
(  304)       ELSE                                    !      NO  budget correction 
(  305)          !                                    ! ==============================
(  306)          IF(lwp) WRITE(numout,*)
(  307)          IF(lwp) WRITE(numout,*)'flx_fwb_init : NO freshwater budget correction'
(  308)          empold = 0.e0
(  309)       ENDIF
(  310) 
(  311)    END SUBROUTINE flx_fwb_init
(  312) 
(  313) # 326
(  326)    !!======================================================================
(  327) END MODULE flxfwb
