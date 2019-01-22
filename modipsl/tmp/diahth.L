


PGF90 (Version     10.5)          09/23/2010  09:49:28      page 1

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

Filename: diahth.F90

(    1) MODULE diahth
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  diahth  ***
(    4)    !! Ocean diagnostics: thermocline and 20 degree depth
(    5)    !!======================================================================
(    6) #if   defined key_diahth   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_diahth' :                              thermocline depth diag.
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dia_hth      : Compute diagnostics associated with the thermocline
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE phycst          ! physical constants
(   16)    USE in_out_manager  ! I/O manager
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC dia_hth    ! routine called by step.F90
(   23) 
(   24)    !! * Shared module variables
(   25)    LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .TRUE.   !: thermocline-20d depths flag
(   26)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   27)       hth  ,      &  !: depth of the max vertical temperature gradient (m)
(   28)       hd20 ,      &  !: depth of 20 C isotherm (m)
(   29)       hd28 ,      &  !: depth of 28 C isotherm (m)
(   30)       htc3           !: heat content of first 300 m
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diahth.F90,v 1.3 2005/03/27 18:34:55 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42)    SUBROUTINE dia_hth( kt )
(   43)       !!---------------------------------------------------------------------
(   44)       !!                  ***  ROUTINE dia_hth  ***
(   45)       !!
(   46)       !! ** Purpose :






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 2

(   47)       !!      Computes the depth of strongest vertical temperature gradient
(   48)       !!      Computes the depth of the 20 degree isotherm
(   49)       !!      Computes the depth of the 28 degree isotherm
(   50)       !!      Computes the heat content of first 300 m
(   51)       !!
(   52)       !! ** Method : 
(   53)       !!
(   54)       !! History :
(   55)       !!        !  94-09  (J.-P. Boulanger)  Original code
(   56)       !!        !  96-11  (E. Guilyardi)  OPA8 
(   57)       !!        !  97-08  (G. Madec)  optimization
(   58)       !!        !  99-07  (E. Guilyardi)  hd28 + heat content 
(   59)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   60)       !!-------------------------------------------------------------------
(   61)       !! * Arguments
(   62)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   63) 
(   64)       !! * Local declarations
(   65)       INTEGER :: ji, jj, jk         ! dummy loop arguments
(   66)       INTEGER :: iid, iif, ilevel   ! temporary integers
(   67)       INTEGER, DIMENSION(jpi) ::   idepth
(   68)       INTEGER, DIMENSION(jpi,jpj) ::   ikc
(   69) 
(   70)       REAL(wp) :: zd, zmoy              ! temporary scalars
(   71)       REAL(wp), DIMENSION(jpi) ::   zmax
(   72)       REAL(wp), DIMENSION(jpi,jpk) ::   zdzt
(   73)       !!----------------------------------------------------------------------
(   74) 
(   75)       IF( kt == nit000 ) THEN
(   76)          IF(lwp) WRITE(numout,*)
(   77)          IF(lwp) WRITE(numout,*) 'dia_hth : diagnostics of the thermocline depth'
(   78)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   79)          IF(lwp) WRITE(numout,*)
(   80)       ENDIF
(   81) 
(   82) 
(   83)       ! -------------------------- !
(   84)       !  Depth of the thermocline  !
(   85)       ! -------------------------- !
(   86)       ! The depth of the thermocline is defined as the depth of the 
(   87)       ! strongest vertical temperature gradient
(   88)       
(   89)       DO jj = 1, jpj
(   90)          
(   91)          ! vertical gradient of temperature
(   92)          DO jk = 2, jpkm1
(   93)             zdzt(:,jk) = ( tn(:,jj,jk-1) - tn(:,jj,jk) ) / fse3w(:,jj,jk) * tmask(:,jj,jk)
(   94)          END DO
(   95)          
(   96)          ! search the level of maximum vertical temperature gradient
(   97)          zmax  (:) = 0.e0
(   98)          idepth(:) = 1
(   99)          DO jk = jpkm1, 2, -1
(  100)             DO ji = 1, jpi
(  101)                IF( zdzt(ji,jk) > zmax(ji) ) THEN
(  102)                   zmax  (ji) = zdzt(ji,jk)
(  103)                   idepth(ji) = jk
(  104)                ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 3

(  105)             END DO
(  106)          END DO
(  107) 
(  108)          ! depth of the thermocline
(  109)          DO ji = 1, jpi
(  110)             hth(ji,jj) = fsdepw(ji,jj,idepth(ji))
(  111)          END DO
(  112)          
(  113)       END DO
(  114) 
(  115) 
(  116)       ! ----------------------- !
(  117)       !  Depth of 20C isotherm  !
(  118)       ! ----------------------- !
(  119) 
(  120)       ! initialization to the number of ocean w-point mbathy
(  121)       ! (cf dommsk, minimum value: 1)
(  122)       ikc(:,:) = 1
(  123) 
(  124)       ! search the depth of 20 degrees isotherm
(  125)       ! ( starting from the top, last level above 20C, if not exist, = 1)
(  126)       DO jk = 1, jpkm1
(  127)          DO jj = 1, jpj
(  128)             DO ji = 1, jpi
(  129)                IF( tn(ji,jj,jk) >= 20. ) ikc(ji,jj) = jk
(  130)             END DO
(  131)          END DO
(  132)       END DO
(  133)       
(  134)       ! Depth of 20C isotherm
(  135)       DO jj = 1, jpj
(  136)          DO ji = 1, jpi
(  137)             iid = ikc(ji,jj)
(  138)             iif = mbathy(ji,jj)
(  139)             IF( iid /= 1 ) THEN 
(  140)                ! linear interpolation
(  141)                zd =  fsdept(ji,jj,iid)   &
(  142)                   + (    fsdept(ji,jj,iid+1) - fsdept(ji,jj,iid) )   &
(  143)                   * ( 20.*tmask(ji,jj,iid+1) -     tn(ji,jj,iid) )   &
(  144)                   / (        tn(ji,jj,iid+1) -     tn(ji,jj,iid)    &
(  145)                   + (1.-tmask(ji,jj,1))                       )
(  146)                ! bound by the ocean depth, minimum value, first T-point depth
(  147)                hd20(ji,jj) = MIN( zd*tmask(ji,jj,1), fsdepw(ji,jj,iif))
(  148)             ELSE 
(  149)                hd20(ji,jj)=0.
(  150)             ENDIF
(  151)          END DO
(  152)       END DO
(  153) 
(  154)       ! ----------------------- !
(  155)       !  Depth of 28C isotherm  ! 
(  156)       ! ----------------------- !
(  157)       
(  158)       ! initialization to the number of ocean w-point mbathy
(  159)       ! (cf dommsk, minimum value: 1)
(  160)       ikc(:,:) = 1
(  161)       
(  162)       ! search the depth of 28 degrees isotherm






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 4

(  163)       ! ( starting from the top, last level above 28C, if not exist, = 1)
(  164)       DO jk = 1, jpkm1
(  165)          DO jj = 1, jpj
(  166)             DO ji = 1, jpi
(  167)                IF( tn(ji,jj,jk) >= 28. ) ikc(ji,jj) = jk
(  168)             END DO
(  169)          END DO
(  170)       END DO
(  171)       
(  172)       ! Depth of 28C isotherm
(  173)       DO jj = 1, jpj
(  174)          DO ji = 1, jpi
(  175)             iid = ikc(ji,jj)
(  176)             iif = mbathy(ji,jj)
(  177)             IF( iid /= 1 ) THEN 
(  178)                ! linear interpolation
(  179)                zd =  fsdept(ji,jj,iid)   &
(  180)                   + (    fsdept(ji,jj,iid+1) - fsdept(ji,jj,iid) )   &
(  181)                   * ( 28.*tmask(ji,jj,iid+1) -     tn(ji,jj,iid) )   &
(  182)                   / (        tn(ji,jj,iid+1) -     tn(ji,jj,iid)    &
(  183)                   + ( 1. - tmask(ji,jj,1) )  )
(  184)                ! bound by the ocean depth, minimum value, first T-point depth
(  185)                hd28(ji,jj) = MIN( zd*tmask(ji,jj,1), fsdepw(ji,jj,iif) )
(  186)             ELSE 
(  187)                hd28(ji,jj) = 0.
(  188)             ENDIF
(  189)          END DO
(  190)       END DO
(  191) 
(  192)       ! ----------------------------------------- !
(  193)       !  Heat content of first 300 m (18 levels)  !
(  194)       ! ----------------------------------------- !
(  195) 
(  196)       htc3(:,:) = 0.e0
(  197)       ilevel = 18
(  198)       zmoy = rau0 * rcp * 0.5
(  199)       
(  200)       ! intregrate tn from surface to klevel
(  201) 
(  202)       DO jk = 1, ilevel
(  203)                htc3(:,:) = htc3(:,:)   &
(  204)                          + zmoy * ( tn(:,:,jk) + tn(:,:,jk+1) ) * fse3w(:,:,jk) * tmask(:,:,jk)
(  205)       END DO
(  206) 
(  207)    END SUBROUTINE dia_hth
(  208) 
(  209) #else
(  210)    !!----------------------------------------------------------------------
(  211)    !!   Default option :                                       Empty module
(  212)    !!----------------------------------------------------------------------
(  213)    USE in_out_manager
(  214)    LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .FALSE.  !: thermocline-20d depths flag
(  215) CONTAINS
(  216)    SUBROUTINE dia_hth( kt )         ! Empty routine
(  217)       if(lwp) WRITE(numout,*) 'dia_hth: You should not have seen this print! error?', kt
(  218)    END SUBROUTINE dia_hth
(  219) #endif
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 5

(  221)    !!======================================================================
(  222) END MODULE diahth






























































PGF90 (Version     10.5)          09/23/2010  09:49:28      page 6

(    1) # 1 "diahth.F90"
(    1) MODULE diahth
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  diahth  ***
(    4)    !! Ocean diagnostics: thermocline and 20 degree depth
(    5)    !!======================================================================
(    6) # 210
(  210)    !!----------------------------------------------------------------------
(  211)    !!   Default option :                                       Empty module
(  212)    !!----------------------------------------------------------------------
(  213)    USE in_out_manager
(  214)    LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .FALSE.  !: thermocline-20d depths flag
(  215) CONTAINS
(  216)    SUBROUTINE dia_hth( kt )         ! Empty routine
(  217)       if(lwp) WRITE(numout,*) 'dia_hth: You should not have seen this print! error?', kt
(  218)    END SUBROUTINE dia_hth
(  219) 
(  220) # 221
(  221)    !!======================================================================
(  222) END MODULE diahth
