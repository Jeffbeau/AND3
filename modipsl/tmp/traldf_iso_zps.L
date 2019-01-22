


PGF90 (Version     10.5)          09/23/2010  09:50:15      page 1

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

Filename: traldf_iso_zps.F90

(    1) MODULE traldf_iso_zps
(    2)    !!==============================================================================
(    3)    !!                   ***  MODULE  traldf_iso_zps  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) #if ( defined key_ldfslp   &&   defined key_partial_steps )   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'               slope of the lateral diffusive direction
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_ldf_iso_zps : update the tracer trend with the horizontal 
(   11)    !!                     component of a iso-neutral laplacian operator
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and active tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   17)    USE trdmod          ! ocean active tracers trends 
(   18)    USE trdmod_oce      ! ocean variables trends
(   19)    USE zdf_oce         ! ocean vertical physics
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE ldfslp          ! iso-neutral slopes
(   22)    USE diaptr          ! poleward transport diagnostics
(   23)    USE prtctl          ! Print control
(   24) 
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Accessibility
(   30)    PUBLIC tra_ldf_iso_zps  ! routine called by step.F90
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34) #  include "ldftra_substitute.h90"
(   35) #  include "ldfeiv_substitute.h90"
(   36) #  include "vectopt_loop_substitute.h90"
(   37)    !!----------------------------------------------------------------------
(   38)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_iso_zps.F90,v 1.9 2006/03/20 16:52:22 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------
(   42) 
(   43) CONTAINS
(   44) 
(   45)    SUBROUTINE tra_ldf_iso_zps( kt )
(   46)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 2

(   47)       !!                  ***  ROUTINE tra_ldf_iso_zps  ***
(   48)       !!
(   49)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   50)       !!      trend and add it to the general trend of tracer equation.
(   51)       !!
(   52)       !! ** Method  :   The horizontal component of the lateral diffusive trends 
(   53)       !!      is provided by a 2nd order operator rotated along neural or geopo-
(   54)       !!      tential surfaces to which an eddy induced advection can be added
(   55)       !!      It is computed using before fields (forward in time) and isopyc-
(   56)       !!      nal or geopotential slopes computed in routine ldfslp.
(   57)       !!
(   58)       !!      horizontal fluxes associated with the rotated lateral mixing:
(   59)       !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
(   60)       !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
(   61)       !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
(   62)       !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
(   63)       !!      add horizontal Eddy Induced advective fluxes (lk_traldf_eiv=T):
(   64)       !!         zftu = zftu - dk-1[ aht e2u mi(wslpi) ] mi( tb ) 
(   65)       !!         zftv = zftv - dk-1[ aht e1v mj(wslpj) ] mj( tb ) 
(   66)       !!      take the horizontal divergence of the fluxes:
(   67)       !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
(   68)       !!      Add this trend to the general trend (ta,sa):
(   69)       !!         ta = ta + difft
(   70)       !!
(   71)       !!      'key_trdtra' defined: the trend is saved for diagnostics.
(   72)       !!
(   73)       !!      macro-tasked on horizontal slab (jk-loop).
(   74)       !!
(   75)       !! ** Action :
(   76)       !!         Update (ta,sa) arrays with the before along level biharmonic
(   77)       !!      mixing trend.
(   78)       !!         Save in (ztdta,ztdsa) arrays the trends if 'key_trdtra' defined
(   79)       !!
(   80)       !! History :
(   81)       !!        !  94-08  (G. Madec, M. Imbard)
(   82)       !!        !  97-05  (G. Madec)  split into traldf and trazdf
(   83)       !!   8.5  !  02-08  (G. Madec)  Free form, F90
(   84)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   85)       !!----------------------------------------------------------------------
(   86)       !! * Modules used
(   87)       USE oce           , zftu => ua,  &  ! use ua as workspace
(   88)          &                zfsu => va      ! use va as workspace
(   89) 
(   90)       !! * Arguments
(   91)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   92) 
(   93)       !! * Local declarations
(   94)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(   95)       INTEGER ::   iku, ikv               ! temporary integer
(   96)       REAL(wp) ::   &
(   97)          zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
(   98)          zmsku, zmskv, zbtr, zta, zsa     !    "           "
(   99)       REAL(wp), DIMENSION(jpi,jpj) ::   & ! temporary workspace
(  100)          zdkt , zdk1t, zdks , zdk1s       !    "           "
(  101)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   & 
(  102)          zftv, zgtbu, zgtbv,                &  ! temporary workspace
(  103)          zfsv, zgsbu, zgsbv,                &  !    "           "
(  104)          ztdta, ztdsa






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 3

(  105)          
(  106) #if defined key_traldf_eiv 
(  107)       REAL(wp) ::   &
(  108)          zcg1, zcg2, zuwk, zvwk,            &  ! temporary scalars
(  109)          zuwk1, zvwk1                          !    "           "
(  110)       REAL(wp), DIMENSION(jpi,jpj) ::       &  ! temporary workspace
(  111)          zftug, zftvg, zfsug, zfsvg            !     "        "     
(  112) #endif
(  113)       !!----------------------------------------------------------------------
(  114) 
(  115)       IF( kt == nit000 ) THEN
(  116)          IF(lwp) WRITE(numout,*)
(  117)          IF(lwp) WRITE(numout,*) 'tra_ldf_iso_zps : iso neutral laplacian diffusion in '
(  118)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   z-coordinates with partial steps'
(  119) #if defined key_diaeiv
(  120)          u_eiv(:,:,:) = 0.e0
(  121)          v_eiv(:,:,:) = 0.e0
(  122) #endif
(  123)       ENDIF
(  124) 
(  125)       ! Save ta and sa trends
(  126)       IF( l_trdtra )   THEN
(  127)          ztdta(:,:,:) = ta(:,:,:) 
(  128)          ztdsa(:,:,:) = sa(:,:,:) 
(  129)       ENDIF
(  130) 
(  131)       ! Horizontal temperature and salinity gradient 
(  132)       DO jk = 1, jpk
(  133)          DO jj = 1, jpj-1
(  134)             DO ji = 1, fs_jpim1   ! vector opt.
(  135)                zgtbu(ji,jj,jk) = tmask(ji,jj,jk) * ( tb(ji+1,jj  ,jk) - tb(ji,jj,jk) )
(  136)                zgsbu(ji,jj,jk) = tmask(ji,jj,jk) * ( sb(ji+1,jj  ,jk) - sb(ji,jj,jk) )
(  137)                zgtbv(ji,jj,jk) = tmask(ji,jj,jk) * ( tb(ji  ,jj+1,jk) - tb(ji,jj,jk) )
(  138)                zgsbv(ji,jj,jk) = tmask(ji,jj,jk) * ( sb(ji  ,jj+1,jk) - sb(ji,jj,jk) )
(  139)             END DO
(  140)          END DO
(  141)       END DO
(  142)       ! partial steps correction at the last level 
(  143)       DO jj = 1, jpj-1
(  144)          DO ji = 1, jpi-1
(  145)             ! last level
(  146)             iku = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
(  147)             ikv = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
(  148)             zgtbu(ji,jj,iku) = gtu(ji,jj) 
(  149)             zgsbu(ji,jj,iku) = gsu(ji,jj)               
(  150)             zgtbv(ji,jj,ikv) = gtv(ji,jj) 
(  151)             zgsbv(ji,jj,ikv) = gsv(ji,jj)               
(  152)          END DO
(  153)       END DO
(  154)       
(  155)       !                                                ! ===============
(  156)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  157)          !                                             ! ===============
(  158)          ! 1. Vertical tracer gradient at level jk and jk+1
(  159)          ! ------------------------------------------------
(  160)          ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
(  161) 
(  162)          zdk1t(:,:) = ( tb(:,:,jk) - tb(:,:,jk+1) ) * tmask(:,:,jk+1)






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 4

(  163)          zdk1s(:,:) = ( sb(:,:,jk) - sb(:,:,jk+1) ) * tmask(:,:,jk+1)
(  164) 
(  165)          IF( jk == 1 ) THEN
(  166)             zdkt(:,:) = zdk1t(:,:)
(  167)             zdks(:,:) = zdk1s(:,:)
(  168)          ELSE
(  169)             zdkt(:,:) = ( tb(:,:,jk-1) - tb(:,:,jk) ) * tmask(:,:,jk)
(  170)             zdks(:,:) = ( sb(:,:,jk-1) - sb(:,:,jk) ) * tmask(:,:,jk)
(  171)          ENDIF
(  172) 
(  173) 
(  174)          ! 2. Horizontal fluxes
(  175)          ! --------------------
(  176) 
(  177)          DO jj = 1 , jpjm1
(  178)             DO ji = 1, fs_jpim1   ! vector opt.
(  179)                zabe1 = ( fsahtu(ji,jj,jk) + ahtb0 ) * e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
(  180)                zabe2 = ( fsahtv(ji,jj,jk) + ahtb0 ) * e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
(  181) 
(  182)                zmsku = 1. / MAX(  tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
(  183)                                 + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )
(  184) 
(  185)                zmskv = 1. / MAX(  tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
(  186)                                 + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )
(  187) 
(  188)                zcof1 = -fsahtu(ji,jj,jk) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
(  189)                zcof2 = -fsahtv(ji,jj,jk) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv
(  190) 
(  191)                zftu(ji,jj,jk) = umask(ji,jj,jk) * (  zabe1 * zgtbu(ji,jj,jk)   &
(  192)                   &                                + zcof1 * (  zdkt (ji+1,jj) + zdk1t(ji,jj)      &
(  193)                   &                                           + zdk1t(ji+1,jj) + zdkt (ji,jj)  )  )
(  194)                zftv(ji,jj,jk) = vmask(ji,jj,jk) * (  zabe2 * zgtbv(ji,jj,jk)   &
(  195)                   &                                + zcof2 * (  zdkt (ji,jj+1) + zdk1t(ji,jj)      &
(  196)                   &                                           + zdk1t(ji,jj+1) + zdkt (ji,jj)  )  )
(  197)                zfsu(ji,jj,jk) = umask(ji,jj,jk) * (  zabe1 * zgsbu(ji,jj,jk)   &
(  198)                   &                                + zcof1 * (  zdks (ji+1,jj) + zdk1s(ji,jj)      &
(  199)                   &                                           + zdk1s(ji+1,jj) + zdks (ji,jj)  )  )
(  200)                zfsv(ji,jj,jk) = vmask(ji,jj,jk) * (  zabe2 * zgsbv(ji,jj,jk)   &
(  201)                   &                                + zcof2 * (  zdks (ji,jj+1) + zdk1s(ji,jj)      &
(  202)                   &                                           + zdk1s(ji,jj+1) + zdks (ji,jj)  )  )
(  203)             END DO
(  204)          END DO
(  205) 
(  206) #if defined key_traldf_eiv
(  207)                                         ! ---------------------------------------!
(  208)                                         ! Eddy induced vertical advective fluxes !
(  209)                                         ! ---------------------------------------!
(  210)             DO jj = 1, jpjm1
(  211)                DO ji = 1, fs_jpim1   ! vector opt.
(  212)                   zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj  ,jk  ) ) * fsaeiu(ji,jj,jk  ) * umask(ji,jj,jk  )
(  213)                   zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj  ,jk+1) ) * fsaeiu(ji,jj,jk+1) * umask(ji,jj,jk+1)
(  214)                   zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji  ,jj+1,jk  ) ) * fsaeiv(ji,jj,jk  ) * vmask(ji,jj,jk  )
(  215)                   zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji  ,jj+1,jk+1) ) * fsaeiv(ji,jj,jk+1) * vmask(ji,jj,jk+1)
(  216) 
(  217)                   zcg1= -0.25 * e2u(ji,jj) * umask(ji,jj,jk) * ( zuwk-zuwk1 )
(  218)                   zcg2= -0.25 * e1v(ji,jj) * vmask(ji,jj,jk) * ( zvwk-zvwk1 )
(  219) 
(  220)                   zftug(ji,jj) = zcg1 * ( tb(ji+1,jj,jk) + tb(ji,jj,jk) )






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 5

(  221)                   zftvg(ji,jj) = zcg2 * ( tb(ji,jj+1,jk) + tb(ji,jj,jk) )
(  222)                   zfsug(ji,jj) = zcg1 * ( sb(ji+1,jj,jk) + sb(ji,jj,jk) )
(  223)                   zfsvg(ji,jj) = zcg2 * ( sb(ji,jj+1,jk) + sb(ji,jj,jk) )
(  224) 
(  225)                   zftu(ji,jj,jk) = zftu(ji,jj,jk) + zftug(ji,jj)
(  226)                   zftv(ji,jj,jk) = zftv(ji,jj,jk) + zftvg(ji,jj)
(  227)                   zfsu(ji,jj,jk) = zfsu(ji,jj,jk) + zfsug(ji,jj)
(  228)                   zfsv(ji,jj,jk) = zfsv(ji,jj,jk) + zfsvg(ji,jj)
(  229) #   if defined key_diaeiv
(  230)                   u_eiv(ji,jj,jk) = -2. * zcg1 / ( e2u(ji,jj) * fse3u(ji,jj,jk) )
(  231)                   v_eiv(ji,jj,jk) = -2. * zcg2 / ( e1v(ji,jj) * fse3v(ji,jj,jk) )
(  232) #   endif
(  233)                END DO
(  234)             END DO
(  235) #endif
(  236) 
(  237)          ! II.4 Second derivative (divergence) and add to the general trend
(  238)          ! ----------------------------------------------------------------
(  239) 
(  240)          DO jj = 2 , jpjm1
(  241)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  242)                zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  243)                zta = zbtr * ( zftu(ji,jj,jk) - zftu(ji-1,jj,jk) + zftv(ji,jj,jk) - zftv(ji,jj-1,jk)  )
(  244)                zsa = zbtr * ( zfsu(ji,jj,jk) - zfsu(ji-1,jj,jk) + zfsv(ji,jj,jk) - zfsv(ji,jj-1,jk)  )
(  245)                ta (ji,jj,jk) = ta (ji,jj,jk) + zta
(  246)                sa (ji,jj,jk) = sa (ji,jj,jk) + zsa
(  247)             END DO
(  248)          END DO
(  249)          !                                          ! ===============
(  250)       END DO                                        !   End of slab  
(  251)       !                                             ! ===============
(  252) 
(  253)       ! save the trends for diagnostic
(  254)       ! save the horizontal diffusive trends
(  255)       IF( l_trdtra )   THEN
(  256) #   if defined key_traldf_eiv
(  257)          DO jk = 1 , jpkm1
(  258)             DO jj = 2 , jpjm1
(  259)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  260)                   zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  261)                   tladi(ji,jj,jk) = ( zftug(ji,jj) - zftug(ji-1,jj  ) ) * zbtr
(  262)                   tladj(ji,jj,jk) = ( zftvg(ji,jj) - zftvg(ji  ,jj-1) ) * zbtr
(  263)                   sladi(ji,jj,jk) = ( zfsug(ji,jj) - zfsug(ji-1,jj  ) ) * zbtr
(  264)                   sladj(ji,jj,jk) = ( zfsvg(ji,jj) - zfsvg(ji  ,jj-1) ) * zbtr
(  265)                END DO
(  266)             END DO
(  267)          END DO
(  268) #   else
(  269)          tladi(:,:,:) = 0.e0
(  270)          tladj(:,:,:) = 0.e0
(  271)          sladi(:,:,:) = 0.e0
(  272)          sladj(:,:,:) = 0.e0
(  273) #   endif
(  274) 
(  275)          ! Substract the eddy induced velocity for T/S
(  276)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tladi(:,:,:) - tladj(:,:,:) 
(  277)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sladi(:,:,:) - sladj(:,:,:) 
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 6

(  279)          CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
(  280)       ENDIF
(  281) 
(  282)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  283)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  284)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  285)       ENDIF
(  286) 
(  287) 
(  288)       !!bug  no separation of diff iso and eiv
(  289)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  290)          ! "zonal" mean lateral diffusive heat and salt transports
(  291)          pht_ldf(:) = ptr_vj( zftv(:,:,:) )
(  292)          pst_ldf(:) = ptr_vj( zfsv(:,:,:) )
(  293)          ! "zonal" mean lateral eddy induced velocity heat and salt transports
(  294) #if defined key_diaeiv
(  295)          pht_eiv(:) = ptr_vj( zftv(:,:,:) )
(  296)          pst_eiv(:) = ptr_vj( zfsv(:,:,:) )
(  297) #endif
(  298)       ENDIF
(  299) 
(  300)    END SUBROUTINE tra_ldf_iso_zps
(  301) 
(  302) #else
(  303)    !!----------------------------------------------------------------------
(  304)    !!   default option :   Dummy code   NO rotation of the diffusive tensor
(  305)    !!----------------------------------------------------------------------
(  306) CONTAINS
(  307)    SUBROUTINE tra_ldf_iso_zps( kt )               ! Empty routine
(  308) !      WRITE(*,*) 'tra_ldf_iso_zps: You should not have seen this print! error?', kt
(  309)    END SUBROUTINE tra_ldf_iso_zps
(  310) #endif
(  311) 
(  312)    !!==============================================================================
(  313) END MODULE traldf_iso_zps





























PGF90 (Version     10.5)          09/23/2010  09:50:15      page 7

(    1) # 1 "traldf_iso_zps.F90"
(    1) MODULE traldf_iso_zps
(    2)    !!==============================================================================
(    3)    !!                   ***  MODULE  traldf_iso_zps  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) # 303
(  303)    !!----------------------------------------------------------------------
(  304)    !!   default option :   Dummy code   NO rotation of the diffusive tensor
(  305)    !!----------------------------------------------------------------------
(  306) CONTAINS
(  307)    SUBROUTINE tra_ldf_iso_zps( kt )               ! Empty routine
(  308) !      WRITE(*,*) 'tra_ldf_iso_zps: You should not have seen this print! error?', kt
(  309)    END SUBROUTINE tra_ldf_iso_zps
(  310) 
(  311) # 312
(  312)    !!==============================================================================
(  313) END MODULE traldf_iso_zps
