


PGF90 (Version     10.5)          09/23/2010  09:49:38      page 1

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

Filename: obcdyn_bt.F90

(    1) MODULE obcdyn_bt
(    2) #if ( defined key_dynspg_ts || defined key_dynspg_exp ) & defined key_obc
(    3)    !!=================================================================================
(    4)    !!                       ***  MODULE  obcdyn_bt  ***
(    5)    !! Ocean dynamics:   Radiation/prescription of sea surface heights
(    6)    !!                   on each open boundary
(    7)    !!=================================================================================
(    8) 
(    9)    !!---------------------------------------------------------------------------------
(   10)    !!   obc_dyn_bt        : call the subroutine for each open boundary
(   11)    !!   obc_dyn_bt_east   : Flather's algorithm at the east open boundary
(   12)    !!   obc_dyn_bt_west   : Flather's algorithm at the west open boundary
(   13)    !!   obc_dyn_bt_north  : Flather's algorithm at the north open boundary
(   14)    !!   obc_dyn_bt_south  : Flather's algorithm at the south open boundary
(   15)    !!----------------------------------------------------------------------------------
(   16) 
(   17)    !!----------------------------------------------------------------------------------
(   18)    !! * Modules used
(   19)    USE oce             ! ocean dynamics and tracers 
(   20)    USE dom_oce         ! ocean space and time domain
(   21)    USE phycst          ! physical constants
(   22)    USE obc_oce         ! ocean open boundary conditions
(   23)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   24)    USE lib_mpp         ! distributed memory computing
(   25)    USE obccli          ! ocean open boundary conditions: climatology
(   26)    USE in_out_manager  ! I/O manager
(   27)    USE dynspg_oce      ! surface pressure gradient     (free surface with time-splitting)
(   28) 
(   29)    IMPLICIT NONE
(   30)    PRIVATE
(   31) 
(   32)    !! * Accessibility
(   33)    PUBLIC obc_dyn_bt  ! routine called in dynnxt (explicit free surface case)
(   34) 
(   35)    !!---------------------------------------------------------------------------------
(   36)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdyn_bt.F90,v 1.1 2005/12/28 09:25:07 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE obc_dyn_bt ( kt )
(   44)       !!------------------------------------------------------------------------------
(   45)       !!                      SUBROUTINE obc_dyn_bt
(   46)       !!                     ***********************






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 2

(   47)       !! ** Purpose :
(   48)       !!      Apply Flather's algorithm at open boundaries for the explicit
(   49)       !!      free surface case and free surface case with time-splitting
(   50)       !!
(   51)       !!      This routine is called in dynnxt.F routine and updates ua, va and sshn. 
(   52)       !!
(   53)       !!      The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north, 
(   54)       !!      and/or lp_obc_south allow the user to determine which boundary is an
(   55)       !!      open one (must be done in the param_obc.h90 file).
(   56)       !!
(   57)       !! ** Reference : 
(   58)       !!         Flather, R. A., 1976, Mem. Soc. R. Sci. Liege, Ser. 6, 10, 141-164
(   59)       !!
(   60)       !! History :
(   61)       !!   9.0  !  05-12  (V. Garnier) original 
(   62)       !!----------------------------------------------------------------------
(   63)       !! * Arguments
(   64)       INTEGER, INTENT( in ) ::   kt
(   65) 
(   66)       !!----------------------------------------------------------------------
(   67) 
(   68)       IF( lp_obc_east  )   CALL obc_dyn_bt_east 
(   69)       IF( lp_obc_west  )   CALL obc_dyn_bt_west 
(   70)       IF( lp_obc_north )   CALL obc_dyn_bt_north
(   71)       IF( lp_obc_south )   CALL obc_dyn_bt_south
(   72) 
(   73)       IF( lk_mpp ) THEN
(   74)          IF( kt >= nit000+3 .AND. ln_rstart ) THEN
(   75)             CALL lbc_lnk( sshb, 'T',  1. )
(   76)             CALL lbc_lnk( ub  , 'U', -1. )
(   77)             CALL lbc_lnk( vb  , 'V', -1. )
(   78)          END IF
(   79)          CALL lbc_lnk( sshn, 'T',  1. )
(   80)          CALL lbc_lnk( ua  , 'U', -1. )
(   81)          CALL lbc_lnk( va  , 'V', -1. )
(   82)       ENDIF
(   83) 
(   84)    END SUBROUTINE obc_dyn_bt
(   85) 
(   86) # if defined key_dynspg_exp
(   87)    SUBROUTINE obc_dyn_bt_east 
(   88)       !!------------------------------------------------------------------------------
(   89)       !!                  ***  SUBROUTINE obc_dyn_bt_east  ***
(   90)       !!              
(   91)       !! ** Purpose :
(   92)       !!      Apply Flather's algorithm on east OBC velocities ua, va 
(   93)       !!      Fix sea surface height (sshn) on east open boundary
(   94)       !!      The logical lfbceast must be .TRUE.
(   95)       !!
(   96)       !!  History :
(   97)       !!   9.0  !  05-12  (V. Garnier) original
(   98)       !!------------------------------------------------------------------------------
(   99)       !! * Local declaration
(  100)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  101)       REAL(wp) ::   z05cx, ztau, zin
(  102)       !!------------------------------------------------------------------------------
(  103) 
(  104)       DO ji = nie0, nie1






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 3

(  105)          DO jk = 1, jpkm1
(  106)             DO jj = 1, jpj
(  107)                ua(ji,jj,jk) = ua(ji,jj,jk) + sqrt( grav*hur (ji,jj) )               &
(  108)                   &                      * ( ( sshn(ji,jj) + sshn(ji+1,jj) ) * 0.5  &
(  109)                   &                          - sshfoe(jj) ) * uemsk(jj,jk)
(  110)             END DO
(  111)          END DO
(  112)       END DO
(  113)       DO ji = nie0p1, nie1p1
(  114)          DO jj = 1, jpj
(  115)             sshn(ji,jj) = sshn(ji,jj) * (1.-temsk(jj,1)) + temsk(jj,1)*sshfoe(jj)
(  116)          END DO
(  117)       END DO
(  118) 
(  119)    END SUBROUTINE obc_dyn_bt_east
(  120) 
(  121) 
(  122)    SUBROUTINE obc_dyn_bt_west 
(  123)       !!------------------------------------------------------------------------------
(  124)       !!                  ***  SUBROUTINE obc_dyn_bt_west  ***
(  125)       !!                  
(  126)       !! ** Purpose :
(  127)       !!      Apply Flather algorithm on west OBC velocities ua, va
(  128)       !!      Fix sea surface height (sshn) on west open boundary
(  129)       !!      The logical lfbcwest must be .TRUE.
(  130)       !!
(  131)       !!  History :
(  132)       !!   9.0  !  05-12  (V. Garnier) original
(  133)       !!------------------------------------------------------------------------------
(  134)       !! * Local declaration
(  135)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  136)       REAL(wp) ::   z05cx, ztau, zin
(  137)       !!------------------------------------------------------------------------------
(  138) 
(  139)       DO ji = niw0, niw1
(  140)          DO jk = 1, jpkm1
(  141)             DO jj = 1, jpj
(  142)                ua(ji,jj,jk) = ua(ji,jj,jk) - sqrt( grav*hur (ji,jj) )               &
(  143)                   &                      * ( ( sshn(ji,jj) + sshn(ji+1,jj) ) * 0.5  &
(  144)                   &                          - sshfow(jj) ) * uwmsk(jj,jk)
(  145)             END DO
(  146)          END DO
(  147)          DO jj = 1, jpj
(  148)             sshn(ji,jj) = sshn(ji,jj) * (1.-twmsk(jj,1)) + twmsk(jj,1)*sshfow(jj)
(  149)          END DO
(  150)       END DO
(  151) 
(  152)    END SUBROUTINE obc_dyn_bt_west
(  153) 
(  154)    SUBROUTINE obc_dyn_bt_north 
(  155)       !!------------------------------------------------------------------------------
(  156)       !!                ***  SUBROUTINE obc_dyn_bt_north  ***
(  157)       !!
(  158)       !! ** Purpose :
(  159)       !!      Apply Flather algorithm on north OBC velocities ua, va
(  160)       !!      Fix sea surface height (sshn) on north open boundary
(  161)       !!      The logical lfbcnorth must be .TRUE.
(  162)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 4

(  163)       !!  History :
(  164)       !!   9.0  !  05-12  (V. Garnier) original
(  165)       !!------------------------------------------------------------------------------
(  166)       !! * Local declaration
(  167)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  168)       REAL(wp) ::   z05cx, ztau, zin
(  169)       !!------------------------------------------------------------------------------
(  170) 
(  171)       DO jj = njn0, njn1
(  172)          DO jk = 1, jpkm1
(  173)             DO ji = 1, jpi
(  174)                va(ji,jj,jk) = va(ji,jj,jk) + sqrt( grav*hvr (ji,jj) )               &
(  175)                   &                      * ( ( sshn(ji,jj) + sshn(ji,jj+1) ) * 0.5  &
(  176)                   &                          - sshfon(ji) ) * vnmsk(ji,jk)
(  177)             END DO
(  178)          END DO
(  179)       END DO
(  180)       DO jj = njn0p1, njn1p1
(  181)          DO ji = 1, jpi
(  182)             sshn(ji,jj)= sshn(ji,jj) * (1.-tnmsk(ji,1)) + sshfon(ji)*tnmsk(ji,1)
(  183)          END DO
(  184)       END DO
(  185) 
(  186)    END SUBROUTINE obc_dyn_bt_north
(  187) 
(  188)    SUBROUTINE obc_dyn_bt_south 
(  189)       !!------------------------------------------------------------------------------
(  190)       !!                ***  SUBROUTINE obc_dyn_bt_south  ***
(  191)       !!                    
(  192)       !! ** Purpose :
(  193)       !!      Apply Flather algorithm on south OBC velocities ua, va
(  194)       !!      Fix sea surface height (sshn) on south open boundary
(  195)       !!      The logical lfbcsouth must be .TRUE.
(  196)       !!
(  197)       !!  History :
(  198)       !!   9.0  !  05-12  (V. Garnier) original
(  199)       !!------------------------------------------------------------------------------
(  200)       !! * Local declaration
(  201)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  202)       REAL(wp) ::   z05cx, ztau, zin
(  203) 
(  204)       !!------------------------------------------------------------------------------
(  205) 
(  206)       DO jj = njs0, njs1
(  207)          DO jk = 1, jpkm1
(  208)             DO ji = 1, jpi
(  209)                va(ji,jj,jk) = va(ji,jj,jk) - sqrt( grav*hvr (ji,jj) )               &
(  210)                   &                       * ( ( sshn(ji,jj) + sshn(ji,jj+1) ) * 0.5 &
(  211)                   &                           - sshfos(ji) ) * vsmsk(ji,jk)
(  212)             END DO
(  213)          END DO
(  214)          DO ji = 1, jpi
(  215)             sshn(ji,jj)= sshn(ji,jj) * (1.-tsmsk(ji,1)) + tsmsk(ji,1) * sshfos(ji)
(  216)          END DO
(  217)       END DO
(  218) 
(  219)    END SUBROUTINE obc_dyn_bt_south
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 5

(  221) # elif defined key_dynspg_ts
(  222) 
(  223)    SUBROUTINE obc_dyn_bt_east 
(  224)       !!------------------------------------------------------------------------------
(  225)       !!                  ***  SUBROUTINE obc_dyn_bt_east  ***
(  226)       !!
(  227)       !! ** Purpose :
(  228)       !!      Apply Flather's algorithm on east OBC velocities ua, va
(  229)       !!      Fix sea surface height (sshn) on east open boundary
(  230)       !!      The logical lfbceast must be .TRUE.
(  231)       !!
(  232)       !!  History :
(  233)       !!   9.0  !  05-12  (V. Garnier) original
(  234)       !!------------------------------------------------------------------------------
(  235)       !! * Local declaration
(  236)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  237)       REAL(wp) ::   z05cx, ztau, zin
(  238)       !!------------------------------------------------------------------------------
(  239) 
(  240)       DO ji = nie0, nie1
(  241)          DO jk = 1, jpkm1
(  242)             DO jj = 1, jpj
(  243) !ylu               ua(ji,jj,jk) = ( ua(ji,jj,jk) + sshfoe_b(ji,jj) ) * uemsk(jj,jk)
(  244) !ylu               ua(ji,jj,jk) = ( ua(ji,jj,jk) + ua_e(ji,jj) * hur(ji,jj) ) * uemsk(jj,jk)
(  245)             END DO
(  246)          END DO
(  247)       END DO
(  248)       DO ji = nie0p1, nie1p1
(  249)          DO jj = 1, jpj
(  250)             sshn(ji,jj) = sshn(ji,jj) * (1.-temsk(jj,1)) + temsk(jj,1)*sshn_b(ji,jj)
(  251)          END DO
(  252)       END DO
(  253) 
(  254)    END SUBROUTINE obc_dyn_bt_east
(  255) 
(  256)    SUBROUTINE obc_dyn_bt_west 
(  257)       !!------------------------------------------------------------------------------
(  258)       !!                  ***  SUBROUTINE obc_dyn_bt_west  ***
(  259)       !!
(  260)       !! ** Purpose :
(  261)       !! ** Purpose :
(  262)       !!      Apply Flather algorithm on west OBC velocities ua, va
(  263)       !!      Fix sea surface height (sshn) on west open boundary
(  264)       !!      The logical lfbcwest must be .TRUE.
(  265)       !!
(  266)       !!  History :
(  267)       !!   9.0  !  05-12  (V. Garnier) original
(  268)       !!------------------------------------------------------------------------------
(  269)       !! * Local declaration
(  270)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  271)       REAL(wp) ::   z05cx, ztau, zin
(  272)       !!------------------------------------------------------------------------------
(  273) 
(  274)       DO ji = niw0, niw1
(  275)          DO jk = 1, jpkm1
(  276)             DO jj = 1, jpj
(  277) !ylu               ua(ji,jj,jk) = ( ua(ji,jj,jk) + sshfow_b(ji,jj) ) * uwmsk(jj,jk)
(  278) !ylu               ua(ji,jj,jk) = ( ua(ji,jj,jk) + ua_e(ji,jj) * hur(ji,jj) ) * uwmsk(jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 6

(  279)             END DO
(  280)          END DO
(  281)          DO jj = 1, jpj
(  282)             sshn(ji,jj) = sshn(ji,jj) * (1.-twmsk(jj,1)) + twmsk(jj,1)*sshn_b(ji,jj)
(  283)          END DO
(  284)       END DO
(  285) 
(  286)    END SUBROUTINE obc_dyn_bt_west
(  287) 
(  288)    SUBROUTINE obc_dyn_bt_north 
(  289)       !!------------------------------------------------------------------------------
(  290)       !!                     SUBROUTINE obc_dyn_bt_north
(  291)       !!                    *************************
(  292)       !! ** Purpose :
(  293)       !!      Apply Flather algorithm on north OBC velocities ua, va
(  294)       !!      Fix sea surface height (sshn) on north open boundary
(  295)       !!      The logical lfbcnorth must be .TRUE.
(  296)       !!
(  297)       !!  History :
(  298)       !!   9.0  !  05-12  (V. Garnier) original
(  299)       !!------------------------------------------------------------------------------
(  300)       !! * Local declaration
(  301)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  302)       REAL(wp) ::   z05cx, ztau, zin
(  303)       !!------------------------------------------------------------------------------
(  304) 
(  305)       DO jj = njn0, njn1
(  306)          DO jk = 1, jpkm1
(  307)             DO ji = 1, jpi
(  308) !ylu               va(ji,jj,jk) = ( va(ji,jj,jk) + sshfon_b(ji,jj) ) * vnmsk(jj,jk)
(  309) !ylu               va(ji,jj,jk) = ( va(ji,jj,jk) + va_e(ji,jj) * hvr(ji,jj) ) * vnmsk(jj,jk)
(  310)             END DO
(  311)          END DO
(  312)       END DO
(  313)       DO jj = njn0p1, njn1p1
(  314)          DO ji = 1, jpi
(  315)             sshn(ji,jj)= sshn(ji,jj) * (1.-tnmsk(ji,1)) + sshn_b(ji,jj)*tnmsk(ji,1)
(  316)          END DO
(  317)       END DO
(  318) 
(  319)    END SUBROUTINE obc_dyn_bt_north
(  320) 
(  321)    SUBROUTINE obc_dyn_bt_south 
(  322)       !!------------------------------------------------------------------------------
(  323)       !!                     SUBROUTINE obc_dyn_bt_south
(  324)       !!                    *************************
(  325)       !! ** Purpose :
(  326)       !!      Apply Flather algorithm on south OBC velocities ua, va
(  327)       !!      Fix sea surface height (sshn) on south open boundary
(  328)       !!      The logical lfbcsouth must be .TRUE.
(  329)       !!
(  330)       !!  History :
(  331)       !!   9.0  !  05-12  (V. Garnier) original
(  332)       !!------------------------------------------------------------------------------
(  333)       !! * Local declaration
(  334)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  335)       REAL(wp) ::   z05cx, ztau, zin
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 7

(  337)       !!------------------------------------------------------------------------------
(  338) 
(  339)       DO jj = njs0, njs1
(  340)          DO jk = 1, jpkm1
(  341)             DO ji = 1, jpi
(  342) !ylu               va(ji,jj,jk) = ( va(ji,jj,jk) + sshfos_b(ji,jj) ) * vsmsk(jj,jk)
(  343) !ylu               va(ji,jj,jk) = ( va(ji,jj,jk) + va_e(ji,jj)* hvr(ji,jj)) * vsmsk(jj,jk)
(  344)             END DO
(  345)          END DO
(  346)          DO ji = 1, jpi
(  347)             sshn(ji,jj)= sshn(ji,jj) * (1.-tsmsk(ji,1)) + tsmsk(ji,1) * sshn_b(ji,jj)
(  348)          END DO
(  349)       END DO
(  350) 
(  351)    END SUBROUTINE obc_dyn_bt_south
(  352) 
(  353) # endif
(  354) #else
(  355)    !!=================================================================================
(  356)    !!                       ***  MODULE  obcdyn_bt  ***
(  357)    !! Ocean dynamics:   Radiation of velocities on each open boundary
(  358)    !!=================================================================================
(  359) CONTAINS
(  360) 
(  361)    SUBROUTINE obc_dyn_bt
(  362)                               ! No open boundaries ==> empty routine
(  363)    END SUBROUTINE obc_dyn_bt
(  364) #endif
(  365) 
(  366) END MODULE obcdyn_bt


































PGF90 (Version     10.5)          09/23/2010  09:49:38      page 8

(    1) # 1 "obcdyn_bt.F90"
(    1) MODULE obcdyn_bt
(    2) # 3
(    3)    !!=================================================================================
(    4)    !!                       ***  MODULE  obcdyn_bt  ***
(    5)    !! Ocean dynamics:   Radiation/prescription of sea surface heights
(    6)    !!                   on each open boundary
(    7)    !!=================================================================================
(    8) 
(    9)    !!---------------------------------------------------------------------------------
(   10)    !!   obc_dyn_bt        : call the subroutine for each open boundary
(   11)    !!   obc_dyn_bt_east   : Flather's algorithm at the east open boundary
(   12)    !!   obc_dyn_bt_west   : Flather's algorithm at the west open boundary
(   13)    !!   obc_dyn_bt_north  : Flather's algorithm at the north open boundary
(   14)    !!   obc_dyn_bt_south  : Flather's algorithm at the south open boundary
(   15)    !!----------------------------------------------------------------------------------
(   16) 
(   17)    !!----------------------------------------------------------------------------------
(   18)    !! * Modules used
(   19)    USE oce             ! ocean dynamics and tracers 
(   20)    USE dom_oce         ! ocean space and time domain
(   21)    USE phycst          ! physical constants
(   22)    USE obc_oce         ! ocean open boundary conditions
(   23)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   24)    USE lib_mpp         ! distributed memory computing
(   25)    USE obccli          ! ocean open boundary conditions: climatology
(   26)    USE in_out_manager  ! I/O manager
(   27)    USE dynspg_oce      ! surface pressure gradient     (free surface with time-splitting)
(   28) 
(   29)    IMPLICIT NONE
(   30)    PRIVATE
(   31) 
(   32)    !! * Accessibility
(   33)    PUBLIC obc_dyn_bt  ! routine called in dynnxt (explicit free surface case)
(   34) 
(   35)    !!---------------------------------------------------------------------------------
(   36)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdyn_bt.F90,v 1.1 2005/12/28 09:25:07 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE obc_dyn_bt ( kt )
(   44)       !!------------------------------------------------------------------------------
(   45)       !!                      SUBROUTINE obc_dyn_bt
(   46)       !!                     ***********************
(   47)       !! ** Purpose :
(   48)       !!      Apply Flather's algorithm at open boundaries for the explicit
(   49)       !!      free surface case and free surface case with time-splitting
(   50)       !!
(   51)       !!      This routine is called in dynnxt.F routine and updates ua, va and sshn. 
(   52)       !!
(   53)       !!      The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north, 
(   54)       !!      and/or lp_obc_south allow the user to determine which boundary is an
(   55)       !!      open one (must be done in the param_obc.h90 file).
(   56)       !!
(   57)       !! ** Reference : 






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 9

(   58)       !!         Flather, R. A., 1976, Mem. Soc. R. Sci. Liege, Ser. 6, 10, 141-164
(   59)       !!
(   60)       !! History :
(   61)       !!   9.0  !  05-12  (V. Garnier) original 
(   62)       !!----------------------------------------------------------------------
(   63)       !! * Arguments
(   64)       INTEGER, INTENT( in ) ::   kt
(   65) 
(   66)       !!----------------------------------------------------------------------
(   67) 
(   68)       IF( lp_obc_east  )   CALL obc_dyn_bt_east 
(   69)       IF( lp_obc_west  )   CALL obc_dyn_bt_west 
(   70)       IF( lp_obc_north )   CALL obc_dyn_bt_north
(   71)       IF( lp_obc_south )   CALL obc_dyn_bt_south
(   72) 
(   73)       IF( lk_mpp ) THEN
(   74)          IF( kt >= nit000+3 .AND. ln_rstart ) THEN
(   75)             CALL lbc_lnk( sshb, 'T',  1. )
(   76)             CALL lbc_lnk( ub  , 'U', -1. )
(   77)             CALL lbc_lnk( vb  , 'V', -1. )
(   78)          END IF
(   79)          CALL lbc_lnk( sshn, 'T',  1. )
(   80)          CALL lbc_lnk( ua  , 'U', -1. )
(   81)          CALL lbc_lnk( va  , 'V', -1. )
(   82)       ENDIF
(   83) 
(   84)    END SUBROUTINE obc_dyn_bt
(   85) 
(   86) 
(   87) # 223
(  223)    SUBROUTINE obc_dyn_bt_east 
(  224)       !!------------------------------------------------------------------------------
(  225)       !!                  ***  SUBROUTINE obc_dyn_bt_east  ***
(  226)       !!
(  227)       !! ** Purpose :
(  228)       !!      Apply Flather's algorithm on east OBC velocities ua, va
(  229)       !!      Fix sea surface height (sshn) on east open boundary
(  230)       !!      The logical lfbceast must be .TRUE.
(  231)       !!
(  232)       !!  History :
(  233)       !!   9.0  !  05-12  (V. Garnier) original
(  234)       !!------------------------------------------------------------------------------
(  235)       !! * Local declaration
(  236)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  237)       REAL(wp) ::   z05cx, ztau, zin
(  238)       !!------------------------------------------------------------------------------
(  239) 
(  240)       DO ji = nie0, nie1
(  241)          DO jk = 1, jpkm1
(  242)             DO jj = 1, jpj
(  243) !ylu               ua(ji,jj,jk) = ( ua(ji,jj,jk) + sshfoe_b(ji,jj) ) * uemsk(jj,jk)
(  244) !ylu               ua(ji,jj,jk) = ( ua(ji,jj,jk) + ua_e(ji,jj) * hur(ji,jj) ) * uemsk(jj,jk)
(  245)             END DO
(  246)          END DO
(  247)       END DO
(  248)       DO ji = nie0p1, nie1p1
(  249)          DO jj = 1, jpj
(  250)             sshn(ji,jj) = sshn(ji,jj) * (1.-temsk(jj,1)) + temsk(jj,1)*sshn_b(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 10

(  251)          END DO
(  252)       END DO
(  253) 
(  254)    END SUBROUTINE obc_dyn_bt_east
(  255) 
(  256)    SUBROUTINE obc_dyn_bt_west 
(  257)       !!------------------------------------------------------------------------------
(  258)       !!                  ***  SUBROUTINE obc_dyn_bt_west  ***
(  259)       !!
(  260)       !! ** Purpose :
(  261)       !! ** Purpose :
(  262)       !!      Apply Flather algorithm on west OBC velocities ua, va
(  263)       !!      Fix sea surface height (sshn) on west open boundary
(  264)       !!      The logical lfbcwest must be .TRUE.
(  265)       !!
(  266)       !!  History :
(  267)       !!   9.0  !  05-12  (V. Garnier) original
(  268)       !!------------------------------------------------------------------------------
(  269)       !! * Local declaration
(  270)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  271)       REAL(wp) ::   z05cx, ztau, zin
(  272)       !!------------------------------------------------------------------------------
(  273) 
(  274)       DO ji = niw0, niw1
(  275)          DO jk = 1, jpkm1
(  276)             DO jj = 1, jpj
(  277) !ylu               ua(ji,jj,jk) = ( ua(ji,jj,jk) + sshfow_b(ji,jj) ) * uwmsk(jj,jk)
(  278) !ylu               ua(ji,jj,jk) = ( ua(ji,jj,jk) + ua_e(ji,jj) * hur(ji,jj) ) * uwmsk(jj,jk)
(  279)             END DO
(  280)          END DO
(  281)          DO jj = 1, jpj
(  282)             sshn(ji,jj) = sshn(ji,jj) * (1.-twmsk(jj,1)) + twmsk(jj,1)*sshn_b(ji,jj)
(  283)          END DO
(  284)       END DO
(  285) 
(  286)    END SUBROUTINE obc_dyn_bt_west
(  287) 
(  288)    SUBROUTINE obc_dyn_bt_north 
(  289)       !!------------------------------------------------------------------------------
(  290)       !!                     SUBROUTINE obc_dyn_bt_north
(  291)       !!                    *************************
(  292)       !! ** Purpose :
(  293)       !!      Apply Flather algorithm on north OBC velocities ua, va
(  294)       !!      Fix sea surface height (sshn) on north open boundary
(  295)       !!      The logical lfbcnorth must be .TRUE.
(  296)       !!
(  297)       !!  History :
(  298)       !!   9.0  !  05-12  (V. Garnier) original
(  299)       !!------------------------------------------------------------------------------
(  300)       !! * Local declaration
(  301)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  302)       REAL(wp) ::   z05cx, ztau, zin
(  303)       !!------------------------------------------------------------------------------
(  304) 
(  305)       DO jj = njn0, njn1
(  306)          DO jk = 1, jpkm1
(  307)             DO ji = 1, jpi
(  308) !ylu               va(ji,jj,jk) = ( va(ji,jj,jk) + sshfon_b(ji,jj) ) * vnmsk(jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 11

(  309) !ylu               va(ji,jj,jk) = ( va(ji,jj,jk) + va_e(ji,jj) * hvr(ji,jj) ) * vnmsk(jj,jk)
(  310)             END DO
(  311)          END DO
(  312)       END DO
(  313)       DO jj = njn0p1, njn1p1
(  314)          DO ji = 1, jpi
(  315)             sshn(ji,jj)= sshn(ji,jj) * (1.-tnmsk(ji,1)) + sshn_b(ji,jj)*tnmsk(ji,1)
(  316)          END DO
(  317)       END DO
(  318) 
(  319)    END SUBROUTINE obc_dyn_bt_north
(  320) 
(  321)    SUBROUTINE obc_dyn_bt_south 
(  322)       !!------------------------------------------------------------------------------
(  323)       !!                     SUBROUTINE obc_dyn_bt_south
(  324)       !!                    *************************
(  325)       !! ** Purpose :
(  326)       !!      Apply Flather algorithm on south OBC velocities ua, va
(  327)       !!      Fix sea surface height (sshn) on south open boundary
(  328)       !!      The logical lfbcsouth must be .TRUE.
(  329)       !!
(  330)       !!  History :
(  331)       !!   9.0  !  05-12  (V. Garnier) original
(  332)       !!------------------------------------------------------------------------------
(  333)       !! * Local declaration
(  334)       INTEGER ::   ji, jj, jk ! dummy loop indices
(  335)       REAL(wp) ::   z05cx, ztau, zin
(  336) 
(  337)       !!------------------------------------------------------------------------------
(  338) 
(  339)       DO jj = njs0, njs1
(  340)          DO jk = 1, jpkm1
(  341)             DO ji = 1, jpi
(  342) !ylu               va(ji,jj,jk) = ( va(ji,jj,jk) + sshfos_b(ji,jj) ) * vsmsk(jj,jk)
(  343) !ylu               va(ji,jj,jk) = ( va(ji,jj,jk) + va_e(ji,jj)* hvr(ji,jj)) * vsmsk(jj,jk)
(  344)             END DO
(  345)          END DO
(  346)          DO ji = 1, jpi
(  347)             sshn(ji,jj)= sshn(ji,jj) * (1.-tsmsk(ji,1)) + tsmsk(ji,1) * sshn_b(ji,jj)
(  348)          END DO
(  349)       END DO
(  350) 
(  351)    END SUBROUTINE obc_dyn_bt_south
(  352) 
(  353) 
(  354) # 366
(  366) END MODULE obcdyn_bt
