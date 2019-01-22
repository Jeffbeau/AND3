


PGF90 (Version     10.5)          09/23/2010  09:50:12      page 1

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

Filename: traadv_tvd.F90

(    1) MODULE traadv_tvd
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_tvd  ***
(    4)    !! Ocean active tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_tvd  : update the tracer trend with the horizontal
(    9)    !!                  and vertical advection trends using a TVD scheme
(   10)    !!   nonosc       : compute monotonic tracer fluxes by a nonoscillatory
(   11)    !!                  algorithm 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and active tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE trdmod          ! ocean active tracers trends 
(   17)    USE trdmod_oce      ! ocean variables trends
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   20)    USE trabbl          ! Advective term of BBL
(   21)    USE lib_mpp
(   22)    USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
(   23)    USE diaptr          ! poleward transport diagnostics
(   24)    USE prtctl          ! Print control
(   25) #ifdef key_RIVER_INPUT
(   26) !!DB
(   27)    USE rivers
(   28) #endif
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33) 
(   34)    !! * Accessibility
(   35)    PUBLIC tra_adv_tvd    ! routine called by step.F90
(   36) 
(   37)    !! * Substitutions
(   38) #  include "domzgr_substitute.h90"
(   39) #  include "vectopt_loop_substitute.h90"
(   40)    !!----------------------------------------------------------------------
(   41)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   42)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_tvd.F90,v 1.12 2006/03/20 16:34:20 opalod Exp $ 
(   43)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   44)    !!----------------------------------------------------------------------
(   45) 
(   46) CONTAINS






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 2

(   47) 
(   48)    SUBROUTINE tra_adv_tvd( kt )
(   49)       !!----------------------------------------------------------------------
(   50)       !!                  ***  ROUTINE tra_adv_tvd  ***
(   51)       !! 
(   52)       !! **  Purpose :   Compute the now trend due to total advection of 
(   53)       !!       tracers and add it to the general trend of tracer equations
(   54)       !!
(   55)       !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
(   56)       !!       corrected flux (monotonic correction)
(   57)       !!       note: - this advection scheme needs a leap-frog time scheme
(   58)       !!
(   59)       !! ** Action : - update (ta,sa) with the now advective tracer trends
(   60)       !!             - save the trends in (ttrdh,strdh) ('key_trdtra')
(   61)       !!
(   62)       !! History :
(   63)       !!        !  95-12  (L. Mortier)  Original code
(   64)       !!        !  00-01  (H. Loukos)  adapted to ORCA 
(   65)       !!        !  00-10  (MA Foujols E.Kestenare)  include file not routine
(   66)       !!        !  00-12  (E. Kestenare M. Levy)  fix bug in trtrd indexes
(   67)       !!        !  01-07  (E. Durand G. Madec)  adaptation to ORCA config
(   68)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   69)       !!   9.0  !  04-01  (A. de Miranda, G. Madec, J.M. Molines ): advective bbl
(   70)       !!    "   !  08-04  (S. Cravatte) add the i-, j- & k- trends computation
(   71)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   72)       !!----------------------------------------------------------------------
(   73)       !! * Modules used
(   74)       USE trdmod_oce         , ztay => tladj,  &  ! use tladj latter
(   75)          &                     zsay => sladj,  &  ! use sladj latter
(   76)          &                     ztaz => tladi,  &  ! use ua as workspace
(   77)          &                     zsaz => sladi      ! use ua as workspace
(   78) #if defined key_trabbl_adv
(   79)       USE oce                , zun => ua,  &  ! use ua as workspace
(   80)          &                     zvn => va      ! use va as workspace
(   81) 
(   82)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
(   83) #else
(   84)       USE oce                , zun => un,  &  ! When no bbl, zun == un
(   85)                                zvn => vn,  &  !             zvn == vn
(   86)                                zwn => wn      !             zwn == wn
(   87) #endif
(   88) 
(   89) 
(   90) 
(   91)       !! * Arguments
(   92)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   93) 
(   94)       !! * Local declarations
(   95)       INTEGER  ::   ji, jj, jk              ! dummy loop indices
(   96)       REAL(wp) ::                        &  ! temporary scalar
(   97)          ztai, ztaj, ztak,               &  !    "         "   
(   98)          zsai, zsaj, zsak                   !    "         "   
(   99)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
(  100)          zti, ztu, ztv, ztw,             &  ! temporary workspace
(  101)          zsi, zsu, zsv, zsw,             &  !    "           "
(  102)          ztdta, ztdsa                       !    "           "
(  103)       REAL(wp) ::   &
(  104)          z2dtt, zbtr, zeu, zev, zew, z2, &  ! temporary scalar






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 3

(  105)          zfp_ui, zfp_vj, zfp_wk,         &  !    "         "
(  106)          zfm_ui, zfm_vj, zfm_wk             !    "         "
(  107)       !!----------------------------------------------------------------------
(  108) !!ylu DB -- removed the delta_adv code as it did not work for us
(  109) !      REAL(wp), DIMENSION (jpi,jpj) ::   &
(  110) !         delta_adv  ! temporary workspace
(  111) 
(  112) #ifdef key_RIVER_INPUT
(  113) !!DB: 2009.07.30 -- SLE runoff related (see also dynspg_ts.)
(  114) !!    Put a u,v + a S,T at head to SLE
(  115)       call river_SLE_01(kt)
(  116) #endif
(  117) !!----------------------------------------------------------------------
(  118) 
(  119)       zti(:,:,:) = 0.e0   ;   zsi(:,:,:) = 0.e0
(  120) 
(  121)       IF( kt == nit000 .AND. lwp ) THEN
(  122)          WRITE(numout,*)
(  123)          WRITE(numout,*) 'tra_adv_tvd : TVD advection scheme'
(  124)          WRITE(numout,*) '~~~~~~~~~~~'
(  125)       ENDIF
(  126) 
(  127)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  128)          z2=1.
(  129)       ELSE
(  130)          z2=2.
(  131)       ENDIF
(  132) 
(  133)       ! Save ta and sa trends
(  134)       IF( l_trdtra )   THEN
(  135)          ztdta(:,:,:) = ta(:,:,:) 
(  136)          ztdsa(:,:,:) = sa(:,:,:) 
(  137)          l_adv = 'tvd'
(  138)       ENDIF
(  139) 
(  140) #if defined key_trabbl_adv
(  141)       ! Advective Bottom boundary layer: add the velocity
(  142)       ! -------------------------------------------------
(  143)       zun(:,:,:) = un (:,:,:) - u_bbl(:,:,:)
(  144)       zvn(:,:,:) = vn (:,:,:) - v_bbl(:,:,:)
(  145)       zwn(:,:,:) = wn (:,:,:) + w_bbl(:,:,:)
(  146) #endif
(  147) 
(  148)       ! 1. Bottom value : flux set to zero
(  149)       ! ---------------
(  150)       ztu(:,:,jpk) = 0.e0   ;   zsu(:,:,jpk) = 0.e0
(  151)       ztv(:,:,jpk) = 0.e0   ;   zsv(:,:,jpk) = 0.e0
(  152)       ztw(:,:,jpk) = 0.e0   ;   zsw(:,:,jpk) = 0.e0
(  153)       zti(:,:,jpk) = 0.e0   ;   zsi(:,:,jpk) = 0.e0
(  154) 
(  155) 
(  156)       ! 2. upstream advection with initial mass fluxes & intermediate update
(  157)       ! --------------------------------------------------------------------
(  158)       ! upstream tracer flux in the i and j direction
(  159)       DO jk = 1, jpkm1
(  160)          DO jj = 1, jpjm1
(  161)             DO ji = 1, fs_jpim1   ! vector opt.
(  162)                zeu = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 4

(  163)                zev = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  164)                ! upstream scheme
(  165)                zfp_ui = zeu + ABS( zeu )
(  166)                zfm_ui = zeu - ABS( zeu )
(  167)                zfp_vj = zev + ABS( zev )
(  168)                zfm_vj = zev - ABS( zev )
(  169)                ztu(ji,jj,jk) = zfp_ui * tb(ji,jj,jk) + zfm_ui * tb(ji+1,jj  ,jk)
(  170)                ztv(ji,jj,jk) = zfp_vj * tb(ji,jj,jk) + zfm_vj * tb(ji  ,jj+1,jk)
(  171)                zsu(ji,jj,jk) = zfp_ui * sb(ji,jj,jk) + zfm_ui * sb(ji+1,jj  ,jk)
(  172)                zsv(ji,jj,jk) = zfp_vj * sb(ji,jj,jk) + zfm_vj * sb(ji  ,jj+1,jk)
(  173)             END DO
(  174)          END DO
(  175)       END DO
(  176) 
(  177)       ! upstream tracer flux in the k direction
(  178)       ! Surface value
(  179)       IF( lk_dynspg_rl ) THEN				! rigid lid : flux set to zero
(  180)          ztw(:,:,1) = 0.e0
(  181)          zsw(:,:,1) = 0.e0
(  182)       ELSE						! free surface
(  183)          DO jj = 1, jpj
(  184)             DO ji = 1, jpi
(  185)                zew = e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,1)
(  186)                ztw(ji,jj,1) = zew * tb(ji,jj,1)
(  187)                zsw(ji,jj,1) = zew * sb(ji,jj,1)
(  188)             END DO
(  189)          END DO
(  190)       ENDIF
(  191) 
(  192)       ! Interior value
(  193)       DO jk = 2, jpkm1
(  194)          DO jj = 1, jpj
(  195)             DO ji = 1, jpi
(  196)                zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
(  197)                zfp_wk = zew + ABS( zew )
(  198)                zfm_wk = zew - ABS( zew )
(  199)                ztw(ji,jj,jk) = zfp_wk * tb(ji,jj,jk) + zfm_wk * tb(ji,jj,jk-1)
(  200)                zsw(ji,jj,jk) = zfp_wk * sb(ji,jj,jk) + zfm_wk * sb(ji,jj,jk-1)
(  201)             END DO
(  202)          END DO
(  203)       END DO
(  204) 
(  205)       ! total advective trend
(  206)       DO jk = 1, jpkm1
(  207)          DO jj = 2, jpjm1
(  208)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  209)                zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  210) 
(  211)                ! i- j- horizontal & k- vertical advective trends
(  212)                ztai = - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  ) ) * zbtr
(  213)                ztaj = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr
(  214)                ztak = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr
(  215)                zsai = - ( zsu(ji,jj,jk) - zsu(ji-1,jj  ,jk  ) ) * zbtr
(  216)                zsaj = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr
(  217)                zsak = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr
(  218) 
(  219)                ! total intermediate advective trends
(  220)                zti(ji,jj,jk) = ztai + ztaj + ztak






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 5

(  221)                zsi(ji,jj,jk) = zsai + zsaj + zsak
(  222)             END DO
(  223)          END DO
(  224)       END DO
(  225) 
(  226) 
(  227)      ! Save the intermediate vertical & j- horizontal advection trends
(  228)      IF( l_trdtra )   THEN
(  229)          DO jk = 1, jpkm1
(  230)             DO jj = 2, jpjm1
(  231)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  232)                   zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  233)                   ztay(ji,jj,jk) = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr
(  234)                   zsay(ji,jj,jk) = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr
(  235)                   ztaz(ji,jj,jk) = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr
(  236)                   zsaz(ji,jj,jk) = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr
(  237)                END DO
(  238)             END DO
(  239)          END DO
(  240)       ENDIF
(  241) 
(  242)       ! update and guess with monotonic sheme
(  243)       DO jk = 1, jpkm1
(  244)          z2dtt = z2 * rdttra(jk)
(  245)          DO jj = 2, jpjm1
(  246)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  247)                ta(ji,jj,jk) =  ta(ji,jj,jk) + zti(ji,jj,jk)
(  248)                sa(ji,jj,jk) =  sa(ji,jj,jk) + zsi(ji,jj,jk)
(  249)                zti (ji,jj,jk) = ( tb(ji,jj,jk) + z2dtt * zti(ji,jj,jk) ) * tmask(ji,jj,jk)
(  250)                zsi (ji,jj,jk) = ( sb(ji,jj,jk) + z2dtt * zsi(ji,jj,jk) ) * tmask(ji,jj,jk)
(  251)             END DO
(  252)          END DO
(  253)       END DO
(  254) 
(  255)       ! Lateral boundary conditions on zti, zsi   (unchanged sign)
(  256)       CALL lbc_lnk( zti, 'T', 1. )
(  257)       CALL lbc_lnk( zsi, 'T', 1. )
(  258) 
(  259) 
(  260)       ! 3. antidiffusive flux : high order minus low order
(  261)       ! --------------------------------------------------
(  262)       ! antidiffusive flux on i and j
(  263)       DO jk = 1, jpkm1
(  264)          DO jj = 1, jpjm1
(  265)             DO ji = 1, fs_jpim1   ! vector opt.
(  266)                zeu = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  267)                zev = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  268)                ztu(ji,jj,jk) = zeu * ( tn(ji,jj,jk) + tn(ji+1,jj,jk) ) - ztu(ji,jj,jk)
(  269)                zsu(ji,jj,jk) = zeu * ( sn(ji,jj,jk) + sn(ji+1,jj,jk) ) - zsu(ji,jj,jk)
(  270)                ztv(ji,jj,jk) = zev * ( tn(ji,jj,jk) + tn(ji,jj+1,jk) ) - ztv(ji,jj,jk)
(  271)                zsv(ji,jj,jk) = zev * ( sn(ji,jj,jk) + sn(ji,jj+1,jk) ) - zsv(ji,jj,jk)
(  272)             END DO
(  273)          END DO
(  274)       END DO
(  275)       
(  276)       ! antidiffusive flux on k
(  277)       ! Surface value
(  278)       ztw(:,:,1) = 0.






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 6

(  279)       zsw(:,:,1) = 0.
(  280) 
(  281)       ! Interior value
(  282)       DO jk = 2, jpkm1
(  283)          DO jj = 1, jpj
(  284)             DO ji = 1, jpi
(  285)                zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
(  286)                ztw(ji,jj,jk) = zew * ( tn(ji,jj,jk) + tn(ji,jj,jk-1) ) - ztw(ji,jj,jk)
(  287)                zsw(ji,jj,jk) = zew * ( sn(ji,jj,jk) + sn(ji,jj,jk-1) ) - zsw(ji,jj,jk)
(  288)             END DO
(  289)          END DO
(  290)       END DO
(  291) 
(  292)       ! Lateral bondary conditions
(  293)       CALL lbc_lnk( ztu, 'U', -1. )   ;   CALL lbc_lnk( zsu, 'U', -1. )
(  294)       CALL lbc_lnk( ztv, 'V', -1. )   ;   CALL lbc_lnk( zsv, 'V', -1. )
(  295)       CALL lbc_lnk( ztw, 'W',  1. )   ;   CALL lbc_lnk( zsw, 'W',  1. )
(  296) 
(  297)       ! 4. monotonicity algorithm
(  298)       ! -------------------------
(  299)       CALL nonosc( tb, ztu, ztv, ztw, zti, z2 )
(  300)       CALL nonosc( sb, zsu, zsv, zsw, zsi, z2 )
(  301) 
(  302) 
(  303)       ! 5. final trend with corrected fluxes
(  304)       ! ------------------------------------
(  305)       DO jk = 1, jpkm1
(  306)          DO jj = 2, jpjm1
(  307)             DO ji = fs_2, fs_jpim1   ! vector opt.  
(  308)                zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  309)                ! i- j- horizontal & k- vertical advective trends
(  310)                ztai = - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  )) * zbtr
(  311)                ztaj = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  )) * zbtr
(  312)                ztak = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1)) * zbtr
(  313)                zsai = - ( zsu(ji,jj,jk) - zsu(ji-1,jj  ,jk  )) * zbtr
(  314)                zsaj = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  )) * zbtr
(  315)                zsak = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1)) * zbtr
(  316) 
(  317)                ! add them to the general tracer trends
(  318)                ta(ji,jj,jk) = ta(ji,jj,jk) + ztai + ztaj + ztak
(  319)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsai + zsaj + zsak
(  320)             END DO
(  321)          END DO
(  322)       END DO
(  323) 
(  324)       ! save the advective trends for diagnostic
(  325)       ! tracers trends
(  326)       IF( l_trdtra )   THEN
(  327)          ! Compute the final vertical & j- horizontal advection trends
(  328)          DO jk = 1, jpkm1
(  329)             DO jj = 2, jpjm1
(  330)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  331)                   zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  332)                   ztay(ji,jj,jk) = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr   &
(  333)                      &             + ztay(ji,jj,jk) 
(  334)                   zsay(ji,jj,jk) = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr   &
(  335)                      &             + zsay(ji,jj,jk) 
(  336)                   ztaz(ji,jj,jk) = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr   &






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 7

(  337)                      &             + ztaz(ji,jj,jk) 
(  338)                   zsaz(ji,jj,jk) = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr   &
(  339)                      &             + zsaz(ji,jj,jk) 
(  340)                END DO
(  341)             END DO
(  342)          END DO
(  343) 
(  344)          ! horizontal advection: 
(  345)          ! make the difference between the new trends ta()/sa() and the 
(  346)          ! previous one ztdta()/ztdsa() to have the total advection trends 
(  347)          ! to which we substract the vertical trends ztaz()/zsaz()
(  348)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - ztaz(:,:,:)
(  349)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - zsaz(:,:,:)
(  350) 
(  351)          ! Add the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
(  352)          ztdta(:,:,:) = ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:)
(  353)          ztdsa(:,:,:) = ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)
(  354) 
(  355)          CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)
(  356) 
(  357)          ! vertical advection: 
(  358)          ! Substract the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
(  359)          ztaz(:,:,:) = ztaz(:,:,:) - tn(:,:,:) * hdivn(:,:,:)
(  360)          zsaz(:,:,:) = zsaz(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  361) 
(  362)          CALL trd_mod(ztaz, zsaz, jpttdzad, 'TRA', kt)
(  363) 
(  364)       ENDIF
(  365) 
(  366)       IF(ln_ctl) THEN
(  367)          CALL prt_ctl(tab3d_1=ta, clinfo1=' tvd adv  - Ta: ', mask1=tmask, &
(  368)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  369)       ENDIF
(  370) 
(  371)       ! "zonal" mean advective heat and salt transport
(  372)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  373)          pht_adv(:) = ptr_vj( ztv(:,:,:) )
(  374)          pst_adv(:) = ptr_vj( zsv(:,:,:) )
(  375)       ENDIF
(  376) 
(  377)    END SUBROUTINE tra_adv_tvd
(  378) 
(  379) 
(  380)    SUBROUTINE nonosc( pbef, paa, pbb, pcc, paft, prdt )
(  381)       !!---------------------------------------------------------------------
(  382)       !!                    ***  ROUTINE nonosc  ***
(  383)       !!     
(  384)       !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
(  385)       !!       scheme and the before field by a nonoscillatory algorithm 
(  386)       !!
(  387)       !! **  Method  :   ... ???
(  388)       !!       warning : pbef and paft must be masked, but the boundaries
(  389)       !!       conditions on the fluxes are not necessary zalezak (1979)
(  390)       !!       drange (1995) multi-dimensional forward-in-time and upstream-
(  391)       !!       in-space based differencing for fluid
(  392)       !!
(  393)       !! History :
(  394)       !!        !  97-04  (L. Mortier) Original code






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 8

(  395)       !!        !  00-02  (H. Loukos)  rewritting for opa8
(  396)       !!        !  00-10  (M.A Foujols, E. Kestenare)  lateral b.c.
(  397)       !!        !  01-03  (E. Kestenare)  add key_passivetrc
(  398)       !!        !  01-07  (E. Durand G. Madec)  adapted for T & S
(  399)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  400)       !!----------------------------------------------------------------------
(  401)       !! * Arguments
(  402)       REAL(wp), INTENT( in ) ::   &
(  403)          prdt                               ! ???
(  404)       REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT( inout ) ::   &
(  405)          pbef,                            & ! before field
(  406)          paft,                            & ! after field
(  407)          paa,                             & ! monotonic flux in the i direction
(  408)          pbb,                             & ! monotonic flux in the j direction
(  409)          pcc                                ! monotonic flux in the k direction
(  410) 
(  411)       !! * Local declarations
(  412)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(  413)       INTEGER ::   ikm1
(  414)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   zbetup, zbetdo
(  415)       REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn, z2dtt
(  416)       !!----------------------------------------------------------------------
(  417) 
(  418)       zbig = 1.e+40
(  419)       zrtrn = 1.e-15
(  420)       zbetup(:,:,:) = 0.e0   ;   zbetdo(:,:,:) = 0.e0
(  421) 
(  422)       ! Search local extrema
(  423)       ! --------------------
(  424)       ! large negative value (-zbig) inside land
(  425)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
(  426)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
(  427)       ! search maximum in neighbourhood
(  428)       DO jk = 1, jpkm1
(  429)          ikm1 = MAX(jk-1,1)
(  430)          DO jj = 2, jpjm1
(  431)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  432)                zbetup(ji,jj,jk) = MAX(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
(  433)                   &                     pbef(ji-1,jj  ,jk  ), pbef(ji+1,jj  ,jk  ),   &
(  434)                   &                     paft(ji-1,jj  ,jk  ), paft(ji+1,jj  ,jk  ),   &
(  435)                   &                     pbef(ji  ,jj-1,jk  ), pbef(ji  ,jj+1,jk  ),   &
(  436)                   &                     paft(ji  ,jj-1,jk  ), paft(ji  ,jj+1,jk  ),   &
(  437)                   &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
(  438)                   &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
(  439)             END DO
(  440)          END DO
(  441)       END DO
(  442)       ! large positive value (+zbig) inside land
(  443)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
(  444)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
(  445)       ! search minimum in neighbourhood
(  446)       DO jk = 1, jpkm1
(  447)          ikm1 = MAX(jk-1,1)
(  448)          DO jj = 2, jpjm1
(  449)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  450)                zbetdo(ji,jj,jk) = MIN(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
(  451)                   &                     pbef(ji-1,jj  ,jk  ), pbef(ji+1,jj  ,jk  ),   &
(  452)                   &                     paft(ji-1,jj  ,jk  ), paft(ji+1,jj  ,jk  ),   &






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 9

(  453)                   &                     pbef(ji  ,jj-1,jk  ), pbef(ji  ,jj+1,jk  ),   &
(  454)                   &                     paft(ji  ,jj-1,jk  ), paft(ji  ,jj+1,jk  ),   &
(  455)                   &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
(  456)                   &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
(  457)             END DO
(  458)          END DO
(  459)       END DO
(  460) 
(  461)       ! restore masked values to zero
(  462)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:)
(  463)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:)
(  464)  
(  465) 
(  466)       ! 2. Positive and negative part of fluxes and beta terms
(  467)       ! ------------------------------------------------------
(  468) 
(  469)       DO jk = 1, jpkm1
(  470)          z2dtt = prdt * rdttra(jk)
(  471)          DO jj = 2, jpjm1
(  472)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  473)                ! positive & negative part of the flux
(  474)                zpos = MAX( 0., paa(ji-1,jj  ,jk  ) ) - MIN( 0., paa(ji  ,jj  ,jk  ) )   &
(  475)                   & + MAX( 0., pbb(ji  ,jj-1,jk  ) ) - MIN( 0., pbb(ji  ,jj  ,jk  ) )   &
(  476)                   & + MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )
(  477)                zneg = MAX( 0., paa(ji  ,jj  ,jk  ) ) - MIN( 0., paa(ji-1,jj  ,jk  ) )   &
(  478)                   & + MAX( 0., pbb(ji  ,jj  ,jk  ) ) - MIN( 0., pbb(ji  ,jj-1,jk  ) )   &
(  479)                   & + MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )
(  480)                ! up & down beta terms
(  481)                zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) / z2dtt
(  482)                zbetup(ji,jj,jk) = ( zbetup(ji,jj,jk) - paft(ji,jj,jk) ) / (zpos+zrtrn) * zbt
(  483)                zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zbetdo(ji,jj,jk) ) / (zneg+zrtrn) * zbt
(  484)             END DO
(  485)          END DO
(  486)       END DO
(  487) 
(  488)       ! lateral boundary condition on zbetup & zbetdo   (unchanged sign)
(  489)       CALL lbc_lnk( zbetup, 'T', 1. )
(  490)       CALL lbc_lnk( zbetdo, 'T', 1. )
(  491) 
(  492) 
(  493)       ! 3. monotonic flux in the i & j direction (paa & pbb)
(  494)       ! ----------------------------------------
(  495)       DO jk = 1, jpkm1
(  496)          DO jj = 2, jpjm1
(  497)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  498)                za = MIN( 1.e0, zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
(  499)                zb = MIN( 1.e0, zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
(  500)                zc = 0.5 * ( 1.e0 + SIGN( 1.e0, paa(ji,jj,jk) ) )
(  501)                paa(ji,jj,jk) = paa(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
(  502) 
(  503)                za = MIN( 1.e0, zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )
(  504)                zb = MIN( 1.e0, zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )
(  505)                zc = 0.5 * ( 1.e0 + SIGN( 1.e0, pbb(ji,jj,jk) ) )
(  506)                pbb(ji,jj,jk) = pbb(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
(  507)             END DO
(  508)          END DO
(  509)       END DO
(  510) 






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 10

(  511) 
(  512)       ! monotonic flux in the k direction, i.e. pcc
(  513)       ! -------------------------------------------
(  514)       DO jk = 2, jpkm1
(  515)          DO jj = 2, jpjm1
(  516)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  517) 
(  518)                za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj,jk-1) )
(  519)                zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj,jk-1) )
(  520)                zc = 0.5 * ( 1.e0 + SIGN( 1.e0, pcc(ji,jj,jk) ) )
(  521)                pcc(ji,jj,jk) = pcc(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
(  522)             END DO
(  523)          END DO
(  524)       END DO
(  525) 
(  526)       ! lateral boundary condition on paa, pbb, pcc
(  527)       CALL lbc_lnk( paa, 'U', -1. )      ! changed sign
(  528)       CALL lbc_lnk( pbb, 'V', -1. )      ! changed sign
(  529)       CALL lbc_lnk( pcc, 'W',  1. )      ! NO changed sign
(  530) 
(  531)    END SUBROUTINE nonosc
(  532) 
(  533)    !!======================================================================
(  534) END MODULE traadv_tvd








































PGF90 (Version     10.5)          09/23/2010  09:50:12      page 11

(    1) # 1 "traadv_tvd.F90"
(    1) MODULE traadv_tvd
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_tvd  ***
(    4)    !! Ocean active tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_tvd  : update the tracer trend with the horizontal
(    9)    !!                  and vertical advection trends using a TVD scheme
(   10)    !!   nonosc       : compute monotonic tracer fluxes by a nonoscillatory
(   11)    !!                  algorithm 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and active tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE trdmod          ! ocean active tracers trends 
(   17)    USE trdmod_oce      ! ocean variables trends
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   20)    USE trabbl          ! Advective term of BBL
(   21)    USE lib_mpp
(   22)    USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
(   23)    USE diaptr          ! poleward transport diagnostics
(   24)    USE prtctl          ! Print control
(   25) # 26
(   26) !!DB
(   27)    USE rivers
(   28) 
(   29) # 30
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33) 
(   34)    !! * Accessibility
(   35)    PUBLIC tra_adv_tvd    ! routine called by step.F90
(   36) 
(   37)    !! * Substitutions
(   38) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
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
(   48) # 39 "traadv_tvd.F90"
(   39) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 12

(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 40 "traadv_tvd.F90"
(   40) # 40
(   40)    !!----------------------------------------------------------------------
(   41)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   42)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_tvd.F90,v 1.12 2006/03/20 16:34:20 opalod Exp $ 
(   43)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   44)    !!----------------------------------------------------------------------
(   45) 
(   46) CONTAINS
(   47) 
(   48)    SUBROUTINE tra_adv_tvd( kt )
(   49)       !!----------------------------------------------------------------------
(   50)       !!                  ***  ROUTINE tra_adv_tvd  ***
(   51)       !! 
(   52)       !! **  Purpose :   Compute the now trend due to total advection of 
(   53)       !!       tracers and add it to the general trend of tracer equations
(   54)       !!
(   55)       !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
(   56)       !!       corrected flux (monotonic correction)
(   57)       !!       note: - this advection scheme needs a leap-frog time scheme
(   58)       !!
(   59)       !! ** Action : - update (ta,sa) with the now advective tracer trends
(   60)       !!             - save the trends in (ttrdh,strdh) ('key_trdtra')
(   61)       !!
(   62)       !! History :
(   63)       !!        !  95-12  (L. Mortier)  Original code
(   64)       !!        !  00-01  (H. Loukos)  adapted to ORCA 
(   65)       !!        !  00-10  (MA Foujols E.Kestenare)  include file not routine
(   66)       !!        !  00-12  (E. Kestenare M. Levy)  fix bug in trtrd indexes
(   67)       !!        !  01-07  (E. Durand G. Madec)  adaptation to ORCA config
(   68)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   69)       !!   9.0  !  04-01  (A. de Miranda, G. Madec, J.M. Molines ): advective bbl
(   70)       !!    "   !  08-04  (S. Cravatte) add the i-, j- & k- trends computation
(   71)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   72)       !!----------------------------------------------------------------------
(   73)       !! * Modules used
(   74)       USE trdmod_oce         , ztay => tladj,  &  ! use tladj latter
(   75)          &                     zsay => sladj,  &  ! use sladj latter
(   76)          &                     ztaz => tladi,  &  ! use ua as workspace
(   77)          &                     zsaz => sladi      ! use ua as workspace
(   78) # 84
(   84)       USE oce                , zun => un,  &  ! When no bbl, zun == un
(   85)                                zvn => vn,  &  !             zvn == vn
(   86)                                zwn => wn      !             zwn == wn
(   87) 
(   88) 
(   89) 
(   90) # 91
(   91)       !! * Arguments
(   92)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 13

(   93) 
(   94)       !! * Local declarations
(   95)       INTEGER  ::   ji, jj, jk              ! dummy loop indices
(   96)       REAL(wp) ::                        &  ! temporary scalar
(   97)          ztai, ztaj, ztak,               &  !    "         "   
(   98)          zsai, zsaj, zsak                   !    "         "   
(   99)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
(  100)          zti, ztu, ztv, ztw,             &  ! temporary workspace
(  101)          zsi, zsu, zsv, zsw,             &  !    "           "
(  102)          ztdta, ztdsa                       !    "           "
(  103)       REAL(wp) ::   &
(  104)          z2dtt, zbtr, zeu, zev, zew, z2, &  ! temporary scalar
(  105)          zfp_ui, zfp_vj, zfp_wk,         &  !    "         "
(  106)          zfm_ui, zfm_vj, zfm_wk             !    "         "
(  107)       !!----------------------------------------------------------------------
(  108) !!ylu DB -- removed the delta_adv code as it did not work for us
(  109) !      REAL(wp), DIMENSION (jpi,jpj) ::   &
(  110) !         delta_adv  ! temporary workspace
(  111) 
(  112) # 113
(  113) !!DB: 2009.07.30 -- SLE runoff related (see also dynspg_ts.)
(  114) !!    Put a u,v + a S,T at head to SLE
(  115)       call river_SLE_01(kt)
(  116) # 117
(  117) !!----------------------------------------------------------------------
(  118) 
(  119)       zti(:,:,:) = 0.e0   ;   zsi(:,:,:) = 0.e0
(  120) 
(  121)       IF( kt == nit000 .AND. lwp ) THEN
(  122)          WRITE(numout,*)
(  123)          WRITE(numout,*) 'tra_adv_tvd : TVD advection scheme'
(  124)          WRITE(numout,*) '~~~~~~~~~~~'
(  125)       ENDIF
(  126) 
(  127)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  128)          z2=1.
(  129)       ELSE
(  130)          z2=2.
(  131)       ENDIF
(  132) 
(  133)       ! Save ta and sa trends
(  134)       IF( l_trdtra )   THEN
(  135)          ztdta(:,:,:) = ta(:,:,:) 
(  136)          ztdsa(:,:,:) = sa(:,:,:) 
(  137)          l_adv = 'tvd'
(  138)       ENDIF
(  139) 
(  140) 
(  141) # 148
(  148)       ! 1. Bottom value : flux set to zero
(  149)       ! ---------------
(  150)       ztu(:,:,jpk) = 0.e0   ;   zsu(:,:,jpk) = 0.e0
(  151)       ztv(:,:,jpk) = 0.e0   ;   zsv(:,:,jpk) = 0.e0
(  152)       ztw(:,:,jpk) = 0.e0   ;   zsw(:,:,jpk) = 0.e0
(  153)       zti(:,:,jpk) = 0.e0   ;   zsi(:,:,jpk) = 0.e0
(  154) 
(  155) 
(  156)       ! 2. upstream advection with initial mass fluxes & intermediate update






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 14

(  157)       ! --------------------------------------------------------------------
(  158)       ! upstream tracer flux in the i and j direction
(  159)       DO jk = 1, jpkm1
(  160)          DO jj = 1, jpjm1
(  161)             DO ji = 1, jpim1   ! vector opt.
(  162)                zeu = 0.5 * e2u(ji,jj) * e3t(jk) * zun(ji,jj,jk)
(  163)                zev = 0.5 * e1v(ji,jj) * e3t(jk) * zvn(ji,jj,jk)
(  164)                ! upstream scheme
(  165)                zfp_ui = zeu + ABS( zeu )
(  166)                zfm_ui = zeu - ABS( zeu )
(  167)                zfp_vj = zev + ABS( zev )
(  168)                zfm_vj = zev - ABS( zev )
(  169)                ztu(ji,jj,jk) = zfp_ui * tb(ji,jj,jk) + zfm_ui * tb(ji+1,jj  ,jk)
(  170)                ztv(ji,jj,jk) = zfp_vj * tb(ji,jj,jk) + zfm_vj * tb(ji  ,jj+1,jk)
(  171)                zsu(ji,jj,jk) = zfp_ui * sb(ji,jj,jk) + zfm_ui * sb(ji+1,jj  ,jk)
(  172)                zsv(ji,jj,jk) = zfp_vj * sb(ji,jj,jk) + zfm_vj * sb(ji  ,jj+1,jk)
(  173)             END DO
(  174)          END DO
(  175)       END DO
(  176) 
(  177)       ! upstream tracer flux in the k direction
(  178)       ! Surface value
(  179)       IF( lk_dynspg_rl ) THEN				! rigid lid : flux set to zero
(  180)          ztw(:,:,1) = 0.e0
(  181)          zsw(:,:,1) = 0.e0
(  182)       ELSE						! free surface
(  183)          DO jj = 1, jpj
(  184)             DO ji = 1, jpi
(  185)                zew = e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,1)
(  186)                ztw(ji,jj,1) = zew * tb(ji,jj,1)
(  187)                zsw(ji,jj,1) = zew * sb(ji,jj,1)
(  188)             END DO
(  189)          END DO
(  190)       ENDIF
(  191) 
(  192)       ! Interior value
(  193)       DO jk = 2, jpkm1
(  194)          DO jj = 1, jpj
(  195)             DO ji = 1, jpi
(  196)                zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
(  197)                zfp_wk = zew + ABS( zew )
(  198)                zfm_wk = zew - ABS( zew )
(  199)                ztw(ji,jj,jk) = zfp_wk * tb(ji,jj,jk) + zfm_wk * tb(ji,jj,jk-1)
(  200)                zsw(ji,jj,jk) = zfp_wk * sb(ji,jj,jk) + zfm_wk * sb(ji,jj,jk-1)
(  201)             END DO
(  202)          END DO
(  203)       END DO
(  204) 
(  205)       ! total advective trend
(  206)       DO jk = 1, jpkm1
(  207)          DO jj = 2, jpjm1
(  208)             DO ji = 2, jpim1   ! vector opt.
(  209)                zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * e3t(jk) )
(  210) 
(  211)                ! i- j- horizontal & k- vertical advective trends
(  212)                ztai = - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  ) ) * zbtr
(  213)                ztaj = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr
(  214)                ztak = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 15

(  215)                zsai = - ( zsu(ji,jj,jk) - zsu(ji-1,jj  ,jk  ) ) * zbtr
(  216)                zsaj = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr
(  217)                zsak = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr
(  218) 
(  219)                ! total intermediate advective trends
(  220)                zti(ji,jj,jk) = ztai + ztaj + ztak
(  221)                zsi(ji,jj,jk) = zsai + zsaj + zsak
(  222)             END DO
(  223)          END DO
(  224)       END DO
(  225) 
(  226) 
(  227)      ! Save the intermediate vertical & j- horizontal advection trends
(  228)      IF( l_trdtra )   THEN
(  229)          DO jk = 1, jpkm1
(  230)             DO jj = 2, jpjm1
(  231)                DO ji = 2, jpim1   ! vector opt.
(  232)                   zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * e3t(jk) )
(  233)                   ztay(ji,jj,jk) = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr
(  234)                   zsay(ji,jj,jk) = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr
(  235)                   ztaz(ji,jj,jk) = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr
(  236)                   zsaz(ji,jj,jk) = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr
(  237)                END DO
(  238)             END DO
(  239)          END DO
(  240)       ENDIF
(  241) 
(  242)       ! update and guess with monotonic sheme
(  243)       DO jk = 1, jpkm1
(  244)          z2dtt = z2 * rdttra(jk)
(  245)          DO jj = 2, jpjm1
(  246)             DO ji = 2, jpim1   ! vector opt.
(  247)                ta(ji,jj,jk) =  ta(ji,jj,jk) + zti(ji,jj,jk)
(  248)                sa(ji,jj,jk) =  sa(ji,jj,jk) + zsi(ji,jj,jk)
(  249)                zti (ji,jj,jk) = ( tb(ji,jj,jk) + z2dtt * zti(ji,jj,jk) ) * tmask(ji,jj,jk)
(  250)                zsi (ji,jj,jk) = ( sb(ji,jj,jk) + z2dtt * zsi(ji,jj,jk) ) * tmask(ji,jj,jk)
(  251)             END DO
(  252)          END DO
(  253)       END DO
(  254) 
(  255)       ! Lateral boundary conditions on zti, zsi   (unchanged sign)
(  256)       CALL lbc_lnk( zti, 'T', 1. )
(  257)       CALL lbc_lnk( zsi, 'T', 1. )
(  258) 
(  259) 
(  260)       ! 3. antidiffusive flux : high order minus low order
(  261)       ! --------------------------------------------------
(  262)       ! antidiffusive flux on i and j
(  263)       DO jk = 1, jpkm1
(  264)          DO jj = 1, jpjm1
(  265)             DO ji = 1, jpim1   ! vector opt.
(  266)                zeu = 0.5 * e2u(ji,jj) * e3t(jk) * zun(ji,jj,jk)
(  267)                zev = 0.5 * e1v(ji,jj) * e3t(jk) * zvn(ji,jj,jk)
(  268)                ztu(ji,jj,jk) = zeu * ( tn(ji,jj,jk) + tn(ji+1,jj,jk) ) - ztu(ji,jj,jk)
(  269)                zsu(ji,jj,jk) = zeu * ( sn(ji,jj,jk) + sn(ji+1,jj,jk) ) - zsu(ji,jj,jk)
(  270)                ztv(ji,jj,jk) = zev * ( tn(ji,jj,jk) + tn(ji,jj+1,jk) ) - ztv(ji,jj,jk)
(  271)                zsv(ji,jj,jk) = zev * ( sn(ji,jj,jk) + sn(ji,jj+1,jk) ) - zsv(ji,jj,jk)
(  272)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 16

(  273)          END DO
(  274)       END DO
(  275)       
(  276)       ! antidiffusive flux on k
(  277)       ! Surface value
(  278)       ztw(:,:,1) = 0.
(  279)       zsw(:,:,1) = 0.
(  280) 
(  281)       ! Interior value
(  282)       DO jk = 2, jpkm1
(  283)          DO jj = 1, jpj
(  284)             DO ji = 1, jpi
(  285)                zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
(  286)                ztw(ji,jj,jk) = zew * ( tn(ji,jj,jk) + tn(ji,jj,jk-1) ) - ztw(ji,jj,jk)
(  287)                zsw(ji,jj,jk) = zew * ( sn(ji,jj,jk) + sn(ji,jj,jk-1) ) - zsw(ji,jj,jk)
(  288)             END DO
(  289)          END DO
(  290)       END DO
(  291) 
(  292)       ! Lateral bondary conditions
(  293)       CALL lbc_lnk( ztu, 'U', -1. )   ;   CALL lbc_lnk( zsu, 'U', -1. )
(  294)       CALL lbc_lnk( ztv, 'V', -1. )   ;   CALL lbc_lnk( zsv, 'V', -1. )
(  295)       CALL lbc_lnk( ztw, 'W',  1. )   ;   CALL lbc_lnk( zsw, 'W',  1. )
(  296) 
(  297)       ! 4. monotonicity algorithm
(  298)       ! -------------------------
(  299)       CALL nonosc( tb, ztu, ztv, ztw, zti, z2 )
(  300)       CALL nonosc( sb, zsu, zsv, zsw, zsi, z2 )
(  301) 
(  302) 
(  303)       ! 5. final trend with corrected fluxes
(  304)       ! ------------------------------------
(  305)       DO jk = 1, jpkm1
(  306)          DO jj = 2, jpjm1
(  307)             DO ji = 2, jpim1   ! vector opt.  
(  308)                zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * e3t(jk) )
(  309)                ! i- j- horizontal & k- vertical advective trends
(  310)                ztai = - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  )) * zbtr
(  311)                ztaj = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  )) * zbtr
(  312)                ztak = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1)) * zbtr
(  313)                zsai = - ( zsu(ji,jj,jk) - zsu(ji-1,jj  ,jk  )) * zbtr
(  314)                zsaj = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  )) * zbtr
(  315)                zsak = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1)) * zbtr
(  316) 
(  317)                ! add them to the general tracer trends
(  318)                ta(ji,jj,jk) = ta(ji,jj,jk) + ztai + ztaj + ztak
(  319)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsai + zsaj + zsak
(  320)             END DO
(  321)          END DO
(  322)       END DO
(  323) 
(  324)       ! save the advective trends for diagnostic
(  325)       ! tracers trends
(  326)       IF( l_trdtra )   THEN
(  327)          ! Compute the final vertical & j- horizontal advection trends
(  328)          DO jk = 1, jpkm1
(  329)             DO jj = 2, jpjm1
(  330)                DO ji = 2, jpim1   ! vector opt.






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 17

(  331)                   zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * e3t(jk) )
(  332)                   ztay(ji,jj,jk) = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr   &
(  333)                      &             + ztay(ji,jj,jk) 
(  334)                   zsay(ji,jj,jk) = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr   &
(  335)                      &             + zsay(ji,jj,jk) 
(  336)                   ztaz(ji,jj,jk) = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr   &
(  337)                      &             + ztaz(ji,jj,jk) 
(  338)                   zsaz(ji,jj,jk) = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr   &
(  339)                      &             + zsaz(ji,jj,jk) 
(  340)                END DO
(  341)             END DO
(  342)          END DO
(  343) 
(  344)          ! horizontal advection: 
(  345)          ! make the difference between the new trends ta()/sa() and the 
(  346)          ! previous one ztdta()/ztdsa() to have the total advection trends 
(  347)          ! to which we substract the vertical trends ztaz()/zsaz()
(  348)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - ztaz(:,:,:)
(  349)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - zsaz(:,:,:)
(  350) 
(  351)          ! Add the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
(  352)          ztdta(:,:,:) = ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:)
(  353)          ztdsa(:,:,:) = ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)
(  354) 
(  355)          CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)
(  356) 
(  357)          ! vertical advection: 
(  358)          ! Substract the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
(  359)          ztaz(:,:,:) = ztaz(:,:,:) - tn(:,:,:) * hdivn(:,:,:)
(  360)          zsaz(:,:,:) = zsaz(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  361) 
(  362)          CALL trd_mod(ztaz, zsaz, jpttdzad, 'TRA', kt)
(  363) 
(  364)       ENDIF
(  365) 
(  366)       IF(ln_ctl) THEN
(  367)          CALL prt_ctl(tab3d_1=ta, clinfo1=' tvd adv  - Ta: ', mask1=tmask, &
(  368)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  369)       ENDIF
(  370) 
(  371)       ! "zonal" mean advective heat and salt transport
(  372)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  373)          pht_adv(:) = ptr_vj( ztv(:,:,:) )
(  374)          pst_adv(:) = ptr_vj( zsv(:,:,:) )
(  375)       ENDIF
(  376) 
(  377)    END SUBROUTINE tra_adv_tvd
(  378) 
(  379) 
(  380)    SUBROUTINE nonosc( pbef, paa, pbb, pcc, paft, prdt )
(  381)       !!---------------------------------------------------------------------
(  382)       !!                    ***  ROUTINE nonosc  ***
(  383)       !!     
(  384)       !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
(  385)       !!       scheme and the before field by a nonoscillatory algorithm 
(  386)       !!
(  387)       !! **  Method  :   ... ???
(  388)       !!       warning : pbef and paft must be masked, but the boundaries






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 18

(  389)       !!       conditions on the fluxes are not necessary zalezak (1979)
(  390)       !!       drange (1995) multi-dimensional forward-in-time and upstream-
(  391)       !!       in-space based differencing for fluid
(  392)       !!
(  393)       !! History :
(  394)       !!        !  97-04  (L. Mortier) Original code
(  395)       !!        !  00-02  (H. Loukos)  rewritting for opa8
(  396)       !!        !  00-10  (M.A Foujols, E. Kestenare)  lateral b.c.
(  397)       !!        !  01-03  (E. Kestenare)  add key_passivetrc
(  398)       !!        !  01-07  (E. Durand G. Madec)  adapted for T & S
(  399)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  400)       !!----------------------------------------------------------------------
(  401)       !! * Arguments
(  402)       REAL(wp), INTENT( in ) ::   &
(  403)          prdt                               ! ???
(  404)       REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT( inout ) ::   &
(  405)          pbef,                            & ! before field
(  406)          paft,                            & ! after field
(  407)          paa,                             & ! monotonic flux in the i direction
(  408)          pbb,                             & ! monotonic flux in the j direction
(  409)          pcc                                ! monotonic flux in the k direction
(  410) 
(  411)       !! * Local declarations
(  412)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(  413)       INTEGER ::   ikm1
(  414)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   zbetup, zbetdo
(  415)       REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn, z2dtt
(  416)       !!----------------------------------------------------------------------
(  417) 
(  418)       zbig = 1.e+40
(  419)       zrtrn = 1.e-15
(  420)       zbetup(:,:,:) = 0.e0   ;   zbetdo(:,:,:) = 0.e0
(  421) 
(  422)       ! Search local extrema
(  423)       ! --------------------
(  424)       ! large negative value (-zbig) inside land
(  425)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
(  426)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
(  427)       ! search maximum in neighbourhood
(  428)       DO jk = 1, jpkm1
(  429)          ikm1 = MAX(jk-1,1)
(  430)          DO jj = 2, jpjm1
(  431)             DO ji = 2, jpim1   ! vector opt.
(  432)                zbetup(ji,jj,jk) = MAX(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
(  433)                   &                     pbef(ji-1,jj  ,jk  ), pbef(ji+1,jj  ,jk  ),   &
(  434)                   &                     paft(ji-1,jj  ,jk  ), paft(ji+1,jj  ,jk  ),   &
(  435)                   &                     pbef(ji  ,jj-1,jk  ), pbef(ji  ,jj+1,jk  ),   &
(  436)                   &                     paft(ji  ,jj-1,jk  ), paft(ji  ,jj+1,jk  ),   &
(  437)                   &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
(  438)                   &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
(  439)             END DO
(  440)          END DO
(  441)       END DO
(  442)       ! large positive value (+zbig) inside land
(  443)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
(  444)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
(  445)       ! search minimum in neighbourhood
(  446)       DO jk = 1, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 19

(  447)          ikm1 = MAX(jk-1,1)
(  448)          DO jj = 2, jpjm1
(  449)             DO ji = 2, jpim1   ! vector opt.
(  450)                zbetdo(ji,jj,jk) = MIN(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
(  451)                   &                     pbef(ji-1,jj  ,jk  ), pbef(ji+1,jj  ,jk  ),   &
(  452)                   &                     paft(ji-1,jj  ,jk  ), paft(ji+1,jj  ,jk  ),   &
(  453)                   &                     pbef(ji  ,jj-1,jk  ), pbef(ji  ,jj+1,jk  ),   &
(  454)                   &                     paft(ji  ,jj-1,jk  ), paft(ji  ,jj+1,jk  ),   &
(  455)                   &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
(  456)                   &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
(  457)             END DO
(  458)          END DO
(  459)       END DO
(  460) 
(  461)       ! restore masked values to zero
(  462)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:)
(  463)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:)
(  464)  
(  465) 
(  466)       ! 2. Positive and negative part of fluxes and beta terms
(  467)       ! ------------------------------------------------------
(  468) 
(  469)       DO jk = 1, jpkm1
(  470)          z2dtt = prdt * rdttra(jk)
(  471)          DO jj = 2, jpjm1
(  472)             DO ji = 2, jpim1   ! vector opt.
(  473)                ! positive & negative part of the flux
(  474)                zpos = MAX( 0., paa(ji-1,jj  ,jk  ) ) - MIN( 0., paa(ji  ,jj  ,jk  ) )   &
(  475)                   & + MAX( 0., pbb(ji  ,jj-1,jk  ) ) - MIN( 0., pbb(ji  ,jj  ,jk  ) )   &
(  476)                   & + MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )
(  477)                zneg = MAX( 0., paa(ji  ,jj  ,jk  ) ) - MIN( 0., paa(ji-1,jj  ,jk  ) )   &
(  478)                   & + MAX( 0., pbb(ji  ,jj  ,jk  ) ) - MIN( 0., pbb(ji  ,jj-1,jk  ) )   &
(  479)                   & + MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )
(  480)                ! up & down beta terms
(  481)                zbt = e1t(ji,jj) * e2t(ji,jj) * e3t(jk) / z2dtt
(  482)                zbetup(ji,jj,jk) = ( zbetup(ji,jj,jk) - paft(ji,jj,jk) ) / (zpos+zrtrn) * zbt
(  483)                zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zbetdo(ji,jj,jk) ) / (zneg+zrtrn) * zbt
(  484)             END DO
(  485)          END DO
(  486)       END DO
(  487) 
(  488)       ! lateral boundary condition on zbetup & zbetdo   (unchanged sign)
(  489)       CALL lbc_lnk( zbetup, 'T', 1. )
(  490)       CALL lbc_lnk( zbetdo, 'T', 1. )
(  491) 
(  492) 
(  493)       ! 3. monotonic flux in the i & j direction (paa & pbb)
(  494)       ! ----------------------------------------
(  495)       DO jk = 1, jpkm1
(  496)          DO jj = 2, jpjm1
(  497)             DO ji = 2, jpim1   ! vector opt.
(  498)                za = MIN( 1.e0, zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
(  499)                zb = MIN( 1.e0, zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
(  500)                zc = 0.5 * ( 1.e0 + SIGN( 1.e0, paa(ji,jj,jk) ) )
(  501)                paa(ji,jj,jk) = paa(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
(  502) 
(  503)                za = MIN( 1.e0, zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )
(  504)                zb = MIN( 1.e0, zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 20

(  505)                zc = 0.5 * ( 1.e0 + SIGN( 1.e0, pbb(ji,jj,jk) ) )
(  506)                pbb(ji,jj,jk) = pbb(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
(  507)             END DO
(  508)          END DO
(  509)       END DO
(  510) 
(  511) 
(  512)       ! monotonic flux in the k direction, i.e. pcc
(  513)       ! -------------------------------------------
(  514)       DO jk = 2, jpkm1
(  515)          DO jj = 2, jpjm1
(  516)             DO ji = 2, jpim1   ! vector opt.
(  517) 
(  518)                za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj,jk-1) )
(  519)                zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj,jk-1) )
(  520)                zc = 0.5 * ( 1.e0 + SIGN( 1.e0, pcc(ji,jj,jk) ) )
(  521)                pcc(ji,jj,jk) = pcc(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
(  522)             END DO
(  523)          END DO
(  524)       END DO
(  525) 
(  526)       ! lateral boundary condition on paa, pbb, pcc
(  527)       CALL lbc_lnk( paa, 'U', -1. )      ! changed sign
(  528)       CALL lbc_lnk( pbb, 'V', -1. )      ! changed sign
(  529)       CALL lbc_lnk( pcc, 'W',  1. )      ! NO changed sign
(  530) 
(  531)    END SUBROUTINE nonosc
(  532) 
(  533)    !!======================================================================
(  534) END MODULE traadv_tvd
