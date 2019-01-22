


PGF90 (Version     10.5)          09/23/2010  09:50:14      page 1

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

Filename: traldf_iso.F90

(    1) MODULE traldf_iso
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  traldf_iso  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_ldfslp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_ldf_iso : update the tracer trend with the horizontal component
(   11)    !!                 of iso neutral laplacian operator or horizontal 
(   12)    !!                 laplacian operator in s-coordinate
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables
(   17)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   18)    USE trdmod          ! ocean active tracers trends 
(   19)    USE trdmod_oce      ! ocean variables trends
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE ldfslp          ! iso-neutral slopes
(   22)    USE diaptr          ! poleward transport diagnostics
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC tra_ldf_iso  ! routine called by step.F90
(   30) 
(   31)    !! * Substitutions
(   32) #  include "domzgr_substitute.h90"
(   33) #  include "ldftra_substitute.h90"
(   34) #  include "ldfeiv_substitute.h90"
(   35) #  include "vectopt_loop_substitute.h90"
(   36)    !!----------------------------------------------------------------------
(   37)    !!----------------------------------------------------------------------
(   38)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_iso.F90,v 1.10 2006/03/20 16:52:22 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE tra_ldf_iso( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE tra_ldf_iso  ***






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 2

(   47)       !! 
(   48)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   49)       !!      trend and add it to the general trend of tracer equation.
(   50)       !!
(   51)       !! ** Method  :   The horizontal component of the lateral diffusive trends 
(   52)       !!      is provided by a 2nd order operator rotated along neural or geopo-
(   53)       !!      tential surfaces to which an eddy induced advection can be added
(   54)       !!      It is computed using before fields (forward in time) and isopyc-
(   55)       !!      nal or geopotential slopes computed in routine ldfslp.
(   56)       !!
(   57)       !!      horizontal fluxes associated with the rotated lateral mixing:
(   58)       !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
(   59)       !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
(   60)       !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
(   61)       !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
(   62)       !!      add horizontal Eddy Induced advective fluxes (lk_traldf_eiv=T):
(   63)       !!         zftu = zftu - dk-1[ aht e2u mi(wslpi) ] mi( tb ) 
(   64)       !!         zftv = zftv - dk-1[ aht e1v mj(wslpj) ] mj( tb ) 
(   65)       !!      take the horizontal divergence of the fluxes:
(   66)       !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
(   67)       !!      Add this trend to the general trend (ta,sa):
(   68)       !!         ta = ta + difft
(   69)       !!
(   70)       !! ** Action  : - Update (ta,sa) arrays with the before isopycnal or
(   71)       !!                geopotential s-coord harmonic mixing trend.
(   72)       !!              - Save the trends in (ztdta,ztdsa) ('key_trdtra')
(   73)       !!
(   74)       !! History :
(   75)       !!        !  94-08  (G. Madec, M. Imbard)
(   76)       !!        !  97-05  (G. Madec)  split into traldf and trazdf
(   77)       !!   8.5  !  02-08  (G. Madec)  Free form, F90
(   78)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   79)       !!----------------------------------------------------------------------
(   80)       !! * Modules used
(   81)       USE oce           , zftu => ua,  &  ! use ua as workspace
(   82)          &                zfsu => va      ! use va as workspace
(   83) 
(   84)       !! * Arguments
(   85)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   86) 
(   87)       !! * Local declarations
(   88)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(   89)       REAL(wp) ::   &
(   90)          zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
(   91) #if defined key_traldf_eiv
(   92)          zcg1, zcg2, zuwk, zvwk,       &
(   93)          zuwk1, zvwk1,                 &
(   94) #endif
(   95)          zmsku, zmskv, zbtr, zta, zsa
(   96)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   97)          zdkt, zdk1t,      &  ! workspace
(   98)          zdks, zdk1s
(   99) #if defined key_traldf_eiv
(  100)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  101)          zftug, zftvg,                  &
(  102)          zfsug, zfsvg
(  103) #endif
(  104)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 3

(  105)          zftv, zfsv,                        &  ! workspace
(  106)          ztdta, ztdsa
(  107)       !!----------------------------------------------------------------------
(  108) 
(  109)       IF( kt == nit000 ) THEN
(  110)          IF(lwp) WRITE(numout,*)
(  111)          IF(lwp) WRITE(numout,*) 'tra_ldf_iso : iso neutral lateral diffusion or'
(  112)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   horizontal laplacian diffusion in s-coordinate'
(  113) #if defined key_diaeiv
(  114)          u_eiv(:,:,:) = 0.e0
(  115)          v_eiv(:,:,:) = 0.e0
(  116) #endif
(  117)       ENDIF
(  118) 
(  119)       ! Save ta and sa trends
(  120)       IF( l_trdtra )   THEN
(  121)          ztdta(:,:,:) = ta(:,:,:) 
(  122)          ztdsa(:,:,:) = sa(:,:,:) 
(  123)       ENDIF
(  124) 
(  125)       !                                                ! ===============
(  126)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  127)          !                                             ! ===============
(  128)          ! 1. Vertical tracer gradient at level jk and jk+1
(  129)          ! ------------------------------------------------
(  130)          ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
(  131) 
(  132)          zdk1t(:,:) = ( tb(:,:,jk) - tb(:,:,jk+1) ) * tmask(:,:,jk+1)
(  133)          zdk1s(:,:) = ( sb(:,:,jk) - sb(:,:,jk+1) ) * tmask(:,:,jk+1)
(  134) 
(  135)          IF( jk == 1 ) THEN
(  136)             zdkt(:,:) = zdk1t(:,:)
(  137)             zdks(:,:) = zdk1s(:,:)
(  138)          ELSE
(  139)             zdkt(:,:) = ( tb(:,:,jk-1) - tb(:,:,jk) ) * tmask(:,:,jk)
(  140)             zdks(:,:) = ( sb(:,:,jk-1) - sb(:,:,jk) ) * tmask(:,:,jk)
(  141)          ENDIF
(  142) 
(  143) 
(  144)          ! 2. Horizontal fluxes
(  145)          ! --------------------
(  146) 
(  147)          DO jj = 1 , jpjm1
(  148)             DO ji = 1, fs_jpim1   ! vector opt.
(  149)                zabe1 = ( fsahtu(ji,jj,jk) + ahtb0 ) * e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
(  150)                zabe2 = ( fsahtv(ji,jj,jk) + ahtb0 ) * e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
(  151) 
(  152)                zmsku = 1. / MAX(   tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
(  153)                                  + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )
(  154) 
(  155)                zmskv = 1. / MAX(   tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
(  156)                                  + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )
(  157) 
(  158)                zcof1 = -fsahtu(ji,jj,jk) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
(  159)                zcof2 = -fsahtv(ji,jj,jk) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv
(  160) 
(  161)                zftu(ji,jj,jk) = umask(ji,jj,jk) * (   zabe1 * (   tb(ji+1,jj,jk) - tb(ji,jj,jk)  )   &
(  162)                   &                              + zcof1 * (   zdkt (ji+1,jj) + zdk1t(ji,jj)      &






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 4

(  163)                   &                                          + zdk1t(ji+1,jj) + zdkt (ji,jj)  )  )
(  164) 
(  165)                zftv(ji,jj,jk) = vmask(ji,jj,jk) * (   zabe2 * (   tb(ji,jj+1,jk) - tb(ji,jj,jk)  )   &
(  166)                   &                              + zcof2 * (   zdkt (ji,jj+1) + zdk1t(ji,jj)      &
(  167)                   &                                          + zdk1t(ji,jj+1) + zdkt (ji,jj)  )  )
(  168) 
(  169)                zfsu(ji,jj,jk) = umask(ji,jj,jk) * (   zabe1 * (   sb(ji+1,jj,jk) - sb(ji,jj,jk)  )   &
(  170)                   &                              + zcof1 * (   zdks (ji+1,jj) + zdk1s(ji,jj)      &
(  171)                   &                              + zdk1s(ji+1,jj) + zdks (ji,jj)  )  )
(  172) 
(  173)                zfsv(ji,jj,jk) = vmask(ji,jj,jk) * (   zabe2 * (   sb(ji,jj+1,jk) - sb(ji,jj,jk)  )   &
(  174)                   &                              + zcof2 * (   zdks (ji,jj+1) + zdk1s(ji,jj)      &
(  175)                   &                              + zdk1s(ji,jj+1) + zdks (ji,jj)  )  )
(  176)             END DO
(  177)          END DO
(  178) 
(  179) #   if defined key_traldf_eiv
(  180)          !                              ! ---------------------------------------!
(  181)          !                              ! Eddy induced vertical advective fluxes !
(  182)          !                              ! ---------------------------------------!
(  183)             DO jj = 1, jpjm1
(  184)                DO ji = 1, fs_jpim1   ! vector opt.
(  185)                   zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj,jk  ) ) * fsaeiu(ji,jj,jk  ) * umask(ji,jj,jk  )
(  186)                   zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj,jk+1) ) * fsaeiu(ji,jj,jk+1) * umask(ji,jj,jk+1)
(  187)                   zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji,jj+1,jk  ) ) * fsaeiv(ji,jj,jk  ) * vmask(ji,jj,jk  )
(  188)                   zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji,jj+1,jk+1) ) * fsaeiv(ji,jj,jk+1) * vmask(ji,jj,jk+1)
(  189)    
(  190)                   zcg1= -0.25 * e2u(ji,jj) * umask(ji,jj,jk) * ( zuwk-zuwk1 )
(  191)                   zcg2= -0.25 * e1v(ji,jj) * vmask(ji,jj,jk) * ( zvwk-zvwk1 )
(  192)    
(  193)                   zftug(ji,jj) = zcg1 * ( tb(ji+1,jj,jk) + tb(ji,jj,jk) )
(  194)                   zftvg(ji,jj) = zcg2 * ( tb(ji,jj+1,jk) + tb(ji,jj,jk) )
(  195)                   zfsug(ji,jj) = zcg1 * ( sb(ji+1,jj,jk) + sb(ji,jj,jk) )
(  196)                   zfsvg(ji,jj) = zcg2 * ( sb(ji,jj+1,jk) + sb(ji,jj,jk) )
(  197)    
(  198)                   zftu(ji,jj,jk) = zftu(ji,jj,jk) + zftug(ji,jj)
(  199)                   zftv(ji,jj,jk) = zftv(ji,jj,jk) + zftvg(ji,jj)
(  200)                   zfsu(ji,jj,jk) = zfsu(ji,jj,jk) + zfsug(ji,jj)
(  201)                   zfsv(ji,jj,jk) = zfsv(ji,jj,jk) + zfsvg(ji,jj)
(  202) #   if defined key_diaeiv
(  203)                   u_eiv(ji,jj,jk) = -2. * zcg1 / ( e2u(ji,jj) * fse3u(ji,jj,jk) )
(  204)                   v_eiv(ji,jj,jk) = -2. * zcg2 / ( e1v(ji,jj) * fse3v(ji,jj,jk) )
(  205) #   endif
(  206)                END DO
(  207)             END DO
(  208) #   endif
(  209) 
(  210)          ! II.4 Second derivative (divergence) and add to the general trend
(  211)          ! ----------------------------------------------------------------
(  212) 
(  213)          DO jj = 2 , jpjm1
(  214)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  215)                zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  216)                zta = zbtr * (  zftu(ji,jj,jk) - zftu(ji-1,jj  ,jk)   &
(  217)                   &          + zftv(ji,jj,jk) - zftv(ji  ,jj-1,jk)  )
(  218)                zsa = zbtr * (  zfsu(ji,jj,jk) - zfsu(ji-1,jj  ,jk)   &
(  219)                   &          + zfsv(ji,jj,jk) - zfsv(ji  ,jj-1,jk)  )
(  220)                ta (ji,jj,jk) = ta (ji,jj,jk) + zta






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 5

(  221)                sa (ji,jj,jk) = sa (ji,jj,jk) + zsa
(  222)             END DO
(  223)          END DO
(  224)          !                                          ! ===============
(  225)       END DO                                        !   End of slab  
(  226)       !                                             ! ===============
(  227) 
(  228)       ! save the trends for diagnostic
(  229)       ! save the horizontal diffusive trends
(  230)       IF( l_trdtra )   THEN
(  231) #   if defined key_traldf_eiv
(  232)          DO jk = 1 , jpkm1
(  233)             DO jj = 2 , jpjm1
(  234)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  235)                   zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  236)                   tladi(ji,jj,jk) = ( zftug(ji,jj) - zftug(ji-1,jj  ) ) * zbtr
(  237)                   tladj(ji,jj,jk) = ( zftvg(ji,jj) - zftvg(ji  ,jj-1) ) * zbtr
(  238)                   sladi(ji,jj,jk) = ( zfsug(ji,jj) - zfsug(ji-1,jj  ) ) * zbtr
(  239)                   sladj(ji,jj,jk) = ( zfsvg(ji,jj) - zfsvg(ji  ,jj-1) ) * zbtr
(  240)                END DO
(  241)             END DO
(  242)          END DO
(  243) #   else
(  244)          tladi(:,:,:) = 0.e0
(  245)          tladj(:,:,:) = 0.e0
(  246)          sladi(:,:,:) = 0.e0
(  247)          sladj(:,:,:) = 0.e0
(  248) #   endif
(  249)         
(  250)          ! Substract the eddy induced velocity for T/S
(  251)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tladi(:,:,:) - tladj(:,:,:) 
(  252)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sladi(:,:,:) - sladj(:,:,:) 
(  253) 
(  254)          CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
(  255)       ENDIF
(  256) 
(  257)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  258)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  259)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  260)       ENDIF
(  261) 
(  262)       !!bug  no separation of diff iso and eiv
(  263)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  264)          ! "zonal" mean lateral diffusive heat and salt transports
(  265)          pht_ldf(:) = ptr_vj( zftv(:,:,:) )
(  266)          pst_ldf(:) = ptr_vj( zfsv(:,:,:) )
(  267)          ! "zonal" mean lateral eddy induced velocity heat and salt transports
(  268) #if defined key_diaeiv
(  269)          pht_eiv(:) = ptr_vj( zftv(:,:,:) )
(  270)          pst_eiv(:) = ptr_vj( zfsv(:,:,:) )
(  271) #endif
(  272)       ENDIF
(  273) 
(  274)    END SUBROUTINE tra_ldf_iso
(  275) 
(  276) #else
(  277)    !!----------------------------------------------------------------------
(  278)    !!   Dummy module :             No rotation of the lateral mixing tensor






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 6

(  279)    !!----------------------------------------------------------------------
(  280) CONTAINS
(  281)    SUBROUTINE tra_ldf_iso( kt )               ! Empty routine
(  282) !      WRITE(*,*) 'tra_ldf_iso: You should not have seen this print! error?', kt
(  283)    END SUBROUTINE tra_ldf_iso
(  284) #endif
(  285) 
(  286)    !!==============================================================================
(  287) END MODULE traldf_iso























































PGF90 (Version     10.5)          09/23/2010  09:50:14      page 7

(    1) # 1 "traldf_iso.F90"
(    1) MODULE traldf_iso
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  traldf_iso  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_ldf_iso : update the tracer trend with the horizontal component
(   11)    !!                 of iso neutral laplacian operator or horizontal 
(   12)    !!                 laplacian operator in s-coordinate
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables
(   17)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   18)    USE trdmod          ! ocean active tracers trends 
(   19)    USE trdmod_oce      ! ocean variables trends
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE ldfslp          ! iso-neutral slopes
(   22)    USE diaptr          ! poleward transport diagnostics
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC tra_ldf_iso  ! routine called by step.F90
(   30) 
(   31)    !! * Substitutions
(   32) # 1 "./domzgr_substitute.h90"
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
(   48) # 33 "traldf_iso.F90"
(   33) # 1 "./ldftra_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** ldftra_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_substitute.h90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 8

(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19
(   19) !   'key_traldf_c3d' :                 aht: 3D coefficient
(   20) # 34 "traldf_iso.F90"
(   34) # 1 "./ldfeiv_substitute.h90"
(    1) # 40 "./ldfeiv_substitute.h90"
(   40) !   Default option :                         eiv: dummy variables
(   41) # 35 "traldf_iso.F90"
(   35) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 36 "traldf_iso.F90"
(   36) # 36
(   36)    !!----------------------------------------------------------------------
(   37)    !!----------------------------------------------------------------------
(   38)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_iso.F90,v 1.10 2006/03/20 16:52:22 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE tra_ldf_iso( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE tra_ldf_iso  ***
(   47)       !! 
(   48)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   49)       !!      trend and add it to the general trend of tracer equation.
(   50)       !!
(   51)       !! ** Method  :   The horizontal component of the lateral diffusive trends 
(   52)       !!      is provided by a 2nd order operator rotated along neural or geopo-
(   53)       !!      tential surfaces to which an eddy induced advection can be added
(   54)       !!      It is computed using before fields (forward in time) and isopyc-
(   55)       !!      nal or geopotential slopes computed in routine ldfslp.
(   56)       !!
(   57)       !!      horizontal fluxes associated with the rotated lateral mixing:
(   58)       !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
(   59)       !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
(   60)       !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
(   61)       !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
(   62)       !!      add horizontal Eddy Induced advective fluxes (lk_traldf_eiv=T):
(   63)       !!         zftu = zftu - dk-1[ aht e2u mi(wslpi) ] mi( tb ) 
(   64)       !!         zftv = zftv - dk-1[ aht e1v mj(wslpj) ] mj( tb ) 
(   65)       !!      take the horizontal divergence of the fluxes:
(   66)       !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
(   67)       !!      Add this trend to the general trend (ta,sa):
(   68)       !!         ta = ta + difft
(   69)       !!
(   70)       !! ** Action  : - Update (ta,sa) arrays with the before isopycnal or






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 9

(   71)       !!                geopotential s-coord harmonic mixing trend.
(   72)       !!              - Save the trends in (ztdta,ztdsa) ('key_trdtra')
(   73)       !!
(   74)       !! History :
(   75)       !!        !  94-08  (G. Madec, M. Imbard)
(   76)       !!        !  97-05  (G. Madec)  split into traldf and trazdf
(   77)       !!   8.5  !  02-08  (G. Madec)  Free form, F90
(   78)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   79)       !!----------------------------------------------------------------------
(   80)       !! * Modules used
(   81)       USE oce           , zftu => ua,  &  ! use ua as workspace
(   82)          &                zfsu => va      ! use va as workspace
(   83) 
(   84)       !! * Arguments
(   85)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   86) 
(   87)       !! * Local declarations
(   88)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(   89)       REAL(wp) ::   &
(   90)          zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
(   95)          zmsku, zmskv, zbtr, zta, zsa
(   96)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   97)          zdkt, zdk1t,      &  ! workspace
(   98)          zdks, zdk1s
(   99) # 104
(  104)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  105)          zftv, zfsv,                        &  ! workspace
(  106)          ztdta, ztdsa
(  107)       !!----------------------------------------------------------------------
(  108) 
(  109)       IF( kt == nit000 ) THEN
(  110)          IF(lwp) WRITE(numout,*)
(  111)          IF(lwp) WRITE(numout,*) 'tra_ldf_iso : iso neutral lateral diffusion or'
(  112)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   horizontal laplacian diffusion in s-coordinate'
(  113) # 117
(  117)       ENDIF
(  118) 
(  119)       ! Save ta and sa trends
(  120)       IF( l_trdtra )   THEN
(  121)          ztdta(:,:,:) = ta(:,:,:) 
(  122)          ztdsa(:,:,:) = sa(:,:,:) 
(  123)       ENDIF
(  124) 
(  125)       !                                                ! ===============
(  126)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  127)          !                                             ! ===============
(  128)          ! 1. Vertical tracer gradient at level jk and jk+1
(  129)          ! ------------------------------------------------
(  130)          ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
(  131) 
(  132)          zdk1t(:,:) = ( tb(:,:,jk) - tb(:,:,jk+1) ) * tmask(:,:,jk+1)
(  133)          zdk1s(:,:) = ( sb(:,:,jk) - sb(:,:,jk+1) ) * tmask(:,:,jk+1)
(  134) 
(  135)          IF( jk == 1 ) THEN
(  136)             zdkt(:,:) = zdk1t(:,:)
(  137)             zdks(:,:) = zdk1s(:,:)
(  138)          ELSE
(  139)             zdkt(:,:) = ( tb(:,:,jk-1) - tb(:,:,jk) ) * tmask(:,:,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 10

(  140)             zdks(:,:) = ( sb(:,:,jk-1) - sb(:,:,jk) ) * tmask(:,:,jk)
(  141)          ENDIF
(  142) 
(  143) 
(  144)          ! 2. Horizontal fluxes
(  145)          ! --------------------
(  146) 
(  147)          DO jj = 1 , jpjm1
(  148)             DO ji = 1, jpim1   ! vector opt.
(  149)                zabe1 = ( ahtu(ji,jj,jk) + ahtb0 ) * e2u(ji,jj) * e3t(jk) / e1u(ji,jj)
(  150)                zabe2 = ( ahtv(ji,jj,jk) + ahtb0 ) * e1v(ji,jj) * e3t(jk) / e2v(ji,jj)
(  151) 
(  152)                zmsku = 1. / MAX(   tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
(  153)                                  + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )
(  154) 
(  155)                zmskv = 1. / MAX(   tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
(  156)                                  + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )
(  157) 
(  158)                zcof1 = -ahtu(ji,jj,jk) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
(  159)                zcof2 = -ahtv(ji,jj,jk) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv
(  160) 
(  161)                zftu(ji,jj,jk) = umask(ji,jj,jk) * (   zabe1 * (   tb(ji+1,jj,jk) - tb(ji,jj,jk)  )   &
(  162)                   &                              + zcof1 * (   zdkt (ji+1,jj) + zdk1t(ji,jj)      &
(  163)                   &                                          + zdk1t(ji+1,jj) + zdkt (ji,jj)  )  )
(  164) 
(  165)                zftv(ji,jj,jk) = vmask(ji,jj,jk) * (   zabe2 * (   tb(ji,jj+1,jk) - tb(ji,jj,jk)  )   &
(  166)                   &                              + zcof2 * (   zdkt (ji,jj+1) + zdk1t(ji,jj)      &
(  167)                   &                                          + zdk1t(ji,jj+1) + zdkt (ji,jj)  )  )
(  168) 
(  169)                zfsu(ji,jj,jk) = umask(ji,jj,jk) * (   zabe1 * (   sb(ji+1,jj,jk) - sb(ji,jj,jk)  )   &
(  170)                   &                              + zcof1 * (   zdks (ji+1,jj) + zdk1s(ji,jj)      &
(  171)                   &                              + zdk1s(ji+1,jj) + zdks (ji,jj)  )  )
(  172) 
(  173)                zfsv(ji,jj,jk) = vmask(ji,jj,jk) * (   zabe2 * (   sb(ji,jj+1,jk) - sb(ji,jj,jk)  )   &
(  174)                   &                              + zcof2 * (   zdks (ji,jj+1) + zdk1s(ji,jj)      &
(  175)                   &                              + zdk1s(ji,jj+1) + zdks (ji,jj)  )  )
(  176)             END DO
(  177)          END DO
(  178) 
(  179) 
(  180) # 210
(  210)          ! II.4 Second derivative (divergence) and add to the general trend
(  211)          ! ----------------------------------------------------------------
(  212) 
(  213)          DO jj = 2 , jpjm1
(  214)             DO ji = 2, jpim1   ! vector opt.
(  215)                zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*e3t(jk) )
(  216)                zta = zbtr * (  zftu(ji,jj,jk) - zftu(ji-1,jj  ,jk)   &
(  217)                   &          + zftv(ji,jj,jk) - zftv(ji  ,jj-1,jk)  )
(  218)                zsa = zbtr * (  zfsu(ji,jj,jk) - zfsu(ji-1,jj  ,jk)   &
(  219)                   &          + zfsv(ji,jj,jk) - zfsv(ji  ,jj-1,jk)  )
(  220)                ta (ji,jj,jk) = ta (ji,jj,jk) + zta
(  221)                sa (ji,jj,jk) = sa (ji,jj,jk) + zsa
(  222)             END DO
(  223)          END DO
(  224)          !                                          ! ===============
(  225)       END DO                                        !   End of slab  
(  226)       !                                             ! ===============






PGF90 (Version     10.5)          09/23/2010  09:50:14      page 11

(  227) 
(  228)       ! save the trends for diagnostic
(  229)       ! save the horizontal diffusive trends
(  230)       IF( l_trdtra )   THEN
(  231) # 244
(  244)          tladi(:,:,:) = 0.e0
(  245)          tladj(:,:,:) = 0.e0
(  246)          sladi(:,:,:) = 0.e0
(  247)          sladj(:,:,:) = 0.e0
(  248) # 249
(  249)         
(  250)          ! Substract the eddy induced velocity for T/S
(  251)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tladi(:,:,:) - tladj(:,:,:) 
(  252)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sladi(:,:,:) - sladj(:,:,:) 
(  253) 
(  254)          CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
(  255)       ENDIF
(  256) 
(  257)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  258)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  259)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  260)       ENDIF
(  261) 
(  262)       !!bug  no separation of diff iso and eiv
(  263)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  264)          ! "zonal" mean lateral diffusive heat and salt transports
(  265)          pht_ldf(:) = ptr_vj( zftv(:,:,:) )
(  266)          pst_ldf(:) = ptr_vj( zfsv(:,:,:) )
(  267)          ! "zonal" mean lateral eddy induced velocity heat and salt transports
(  268) # 272
(  272)       ENDIF
(  273) 
(  274)    END SUBROUTINE tra_ldf_iso
(  275) 
(  276) 
(  277) # 286
(  286)    !!==============================================================================
(  287) END MODULE traldf_iso
