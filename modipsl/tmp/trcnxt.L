


PGF90 (Version     10.5)          09/23/2010  09:49:21      page 1

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

Filename: trcnxt.F90

(    1) MODULE trcnxt
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trcnxt  ***
(    4)    !! Ocean passive tracers:  time stepping on passives tracers
(    5)    !!======================================================================
(    6) #if defined key_passivetrc   
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_nxt     : time stepping on passive tracers
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce_trc         ! ocean dynamics and tracers variables
(   12)    USE trc             ! ocean passive tracers variables
(   13)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   14) !!DB
(   15) #ifdef key_BGCM_01
(   16) !   USE bgcm_01_initrc
(   17)    USE lib_bgcm_01
(   18) #else
(   19)    USE trctrp_lec      ! passive tracers transport
(   20) #endif
(   21)    USE prtctl_trc      ! Print control for debbuging
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC trc_nxt          ! routine called by step.F90
(   28)    !!----------------------------------------------------------------------
(   29)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcnxt.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36)    SUBROUTINE trc_nxt( kt )
(   37)       !!----------------------------------------------------------------------
(   38)       !!                   ***  ROUTINE trcnxt  ***
(   39)       !!
(   40)       !! ** Purpose :   Compute the passive tracers fields at the 
(   41)       !!      next time-step from their temporal trends and swap the fields.
(   42)       !! 
(   43)       !! ** Method  :   Apply lateral boundary conditions on (ua,va) through 
(   44)       !!      call to lbc_lnk routine
(   45)       !!   default:
(   46)       !!      arrays swap






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 2

(   47)       !!         (trn) = (tra) ; (tra) = (0,0)
(   48)       !!         (trb) = (trn) 
(   49)       !!
(   50)       !!   For Arakawa or TVD Scheme : 
(   51)       !!      A Asselin time filter applied on now tracers (trn) to avoid
(   52)       !!      the divergence of two consecutive time-steps and tr arrays
(   53)       !!      to prepare the next time_step:
(   54)       !!         (trb) = (trn) + atfp [ (trb) + (tra) - 2 (trn) ]
(   55)       !!         (trn) = (tra) ; (tra) = (0,0)
(   56)       !!
(   57)       !!
(   58)       !! ** Action  : - update trb, trn
(   59)       !!
(   60)       !! History :
(   61)       !!   7.0  !  91-11  (G. Madec)  Original code
(   62)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   63)       !!        !  95-02  (M. Levy)   passive tracers 
(   64)       !!        !  96-02  (G. Madec & M. Imbard)  opa release 8.0
(   65)       !!   8.0  !  96-04  (A. Weaver)  Euler forward step
(   66)       !!   8.2  !  99-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
(   67)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   68)       !!        !  02-11  (C. Talandier, A-M Treguier) Open boundaries
(   69)       !!   9.0  !  04-03  (C. Ethe) passive tracers 
(   70)       !!----------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   73)       !! * Local declarations
(   74)       INTEGER ::   ji, jj, jk,jn   ! dummy loop indices
(   75)       REAL(wp) ::   zfact, ztra    ! temporary scalar
(   76)       CHARACTER (len=22) :: charout
(   77)       !!----------------------------------------------------------------------
(   78) 
(   79)       IF( kt == nittrc000 .AND. lwp ) THEN
(   80)          WRITE(numout,*)
(   81)          WRITE(numout,*) 'trc_nxt : time stepping on passive tracers'
(   82)       ENDIF
(   83) 
(   84) 
(   85)       DO jn = 1, jptra
(   86) 
(   87)          ! 0. Lateral boundary conditions on tra (T-point, unchanged sign)
(   88)          ! ---------------------------------============
(   89)          CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )   
(   90)          
(   91)          !                                                ! ===============
(   92)          DO jk = 1, jpk                                   ! Horizontal slab
(   93)             !                                             ! ===============
(   94)             ! 1. Leap-frog scheme (only in explicit case, otherwise the 
(   95)             ! -------------------  time stepping is already done in trczdf)
(   96)             IF( l_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) THEN
(   97)                zfact = 2. * rdttra(jk) * FLOAT(ndttrc) 
(   98)                IF( neuler == 0 .AND. kt == nittrc000 ) zfact = rdttra(jk) * FLOAT(ndttrc) 
(   99)                tra(:,:,jk,jn) = ( trb(:,:,jk,jn) + zfact * tra(:,:,jk,jn) ) * tmask(:,:,jk)
(  100)             ENDIF
(  101) 
(  102)          END DO
(  103) 
(  104) #if defined key_obc






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 3

(  105) !!DB 
(  106)          if( kt == nittrc000 .AND. lwp ) then
(  107)             IF(lwp) WRITE(numout,*) 'DB:  Passive tracer Open Boundary conditions applied  '
(  108)          endif
(  109) !!DB (see bgcm_01_initrc) 
(  110)          call update_boundary_vals
(  111) !!DB: as part of OBC procedure call   SUBROUTINE bgcm_N_obc ( kt )
(  112) !!NB: Hardwired to jn = 1 --> N 
(  113) !!     (although really only need/want this called once per dt so jn=1 is a safe choice)
(  114)          if(jn == 1) then
(  115)             call bgcm_N_obc ( kt )
(  116)          endif
(  117) 
(  118) #endif
(  119) 
(  120) 
(  121)          DO jk = 1, jpk  
(  122) 
(  123)             ! 2. Time filter and swap of arrays
(  124)             ! ---------------------------------
(  125)             IF ( ln_trcadv_cen2 .OR. ln_trcadv_tvd  ) THEN
(  126) 
(  127)                IF( neuler == 0 .AND. kt == nittrc000 ) THEN
(  128)                   DO jj = 1, jpj
(  129)                      DO ji = 1, jpi
(  130)                         trb(ji,jj,jk,jn) = trn(ji,jj,jk,jn)
(  131)                         trn(ji,jj,jk,jn) = tra(ji,jj,jk,jn)
(  132)                         tra(ji,jj,jk,jn) = 0.
(  133)                      END DO
(  134)                   END DO
(  135) 
(  136)                ELSE
(  137)                   DO jj = 1, jpj
(  138)                      DO ji = 1, jpi
(  139)                         trb(ji,jj,jk,jn) = atfp  * ( trb(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) + atfp1 * trn(ji,jj,jk,jn)
(  140)                         trn(ji,jj,jk,jn) = tra(ji,jj,jk,jn)
(  141)                         tra(ji,jj,jk,jn) = 0.
(  142)                      END DO
(  143)                   END DO
(  144) 
(  145)                ENDIF
(  146) 
(  147)             ELSE
(  148) !  case of smolar scheme or muscl
(  149)                DO jj = 1, jpj
(  150)                   DO ji = 1, jpi
(  151)                      trb(ji,jj,jk,jn) = tra(ji,jj,jk,jn)
(  152)                      trn(ji,jj,jk,jn) = tra(ji,jj,jk,jn)
(  153)                      tra(ji,jj,jk,jn) = 0.
(  154)                   END DO
(  155)                END DO
(  156) 
(  157)             ENDIF
(  158)             !                                             ! ===============
(  159)          END DO                                           !   End of slab
(  160)          !                                                ! ===============
(  161)       END DO
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 4

(  163)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  164)          WRITE(charout, FMT="('nxt')")
(  165)          CALL prt_ctl_trc_info(charout)
(  166)          CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
(  167)       ENDIF
(  168) 
(  169) 
(  170)    END SUBROUTINE trc_nxt
(  171) 
(  172) #else
(  173)    !!----------------------------------------------------------------------
(  174)    !!   Default option                                         Empty module
(  175)    !!----------------------------------------------------------------------
(  176) CONTAINS
(  177)    SUBROUTINE trc_nxt( kt )  
(  178)       INTEGER, INTENT(in) :: kt
(  179) !      WRITE(*,*) 'trc_nxt: You should not have seen this print! error?', kt
(  180)    END SUBROUTINE trc_nxt
(  181) #endif
(  182)    !!======================================================================
(  183) END MODULE trcnxt











































PGF90 (Version     10.5)          09/23/2010  09:49:21      page 5

(    1) # 1 "trcnxt.F90"
(    1) MODULE trcnxt
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trcnxt  ***
(    4)    !! Ocean passive tracers:  time stepping on passives tracers
(    5)    !!======================================================================
(    6) # 173
(  173)    !!----------------------------------------------------------------------
(  174)    !!   Default option                                         Empty module
(  175)    !!----------------------------------------------------------------------
(  176) CONTAINS
(  177)    SUBROUTINE trc_nxt( kt )  
(  178)       INTEGER, INTENT(in) :: kt
(  179) !      WRITE(*,*) 'trc_nxt: You should not have seen this print! error?', kt
(  180)    END SUBROUTINE trc_nxt
(  181) # 182
(  182)    !!======================================================================
(  183) END MODULE trcnxt
