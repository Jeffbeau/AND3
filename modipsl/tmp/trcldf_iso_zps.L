


PGF90 (Version     10.5)          09/23/2010  09:49:20      page 1

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

Filename: trcldf_iso_zps.F90

(    1) MODULE trcldf_iso_zps
(    2)    !!==============================================================================
(    3)    !!                   ***  MODULE  trcldf_iso_zps  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) #if key_passivetrc && ( defined key_ldfslp   &&   defined key_partial_steps )
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'               slope of the lateral diffusive direction
(    9)    !!----------------------------------------------------------------------
(   10)    !!   trc_ldf_iso_zps : update the tracer trend with the horizontal 
(   11)    !!                     component of a iso-neutral laplacian operator
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce_trc             ! ocean dynamics and active tracers variables
(   15)    USE trc                 ! ocean passive tracers variables
(   16)    USE prtctl_trc          ! Print control for debbuging
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC trc_ldf_iso_zps  ! routine called by step.F90
(   23) 
(   24)    !! * Substitutions
(   25) #  include "passivetrc_substitute.h90"
(   26)    !!----------------------------------------------------------------------
(   27)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_iso_zps.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!----------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE trc_ldf_iso_zps( kt )
(   35)       !!----------------------------------------------------------------------
(   36)       !!                  ***  ROUTINE trc_ldf_iso_zps  ***
(   37)       !!
(   38)       !! ** Purpose :   Compute the before horizontal tracer  diffusive 
(   39)       !!      trend and add it to the general trend of tracer equation.
(   40)       !!
(   41)       !! ** Method  :   The horizontal component of the lateral diffusive trends 
(   42)       !!      is provided by a 2nd order operator rotated along neural or geopo-
(   43)       !!      tential surfaces to which an eddy induced advection can be added
(   44)       !!      It is computed using before fields (forward in time) and isopyc-
(   45)       !!      nal or geopotential slopes computed in routine ldfslp.
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 2

(   47)       !!      horizontal fluxes associated with the rotated lateral mixing:
(   48)       !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
(   49)       !!               - aht       e2u*uslp    dk[ mi(mk(trb)) ]
(   50)       !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
(   51)       !!               - aht       e2u*vslp    dk[ mj(mk(trb)) ]
(   52)       !!      add horizontal Eddy Induced advective fluxes (lk_traldf_eiv=T):
(   53)       !!         zftu = zftu - dk-1[ aht e2u mi(wslpi) ] mi( trb ) 
(   54)       !!         zftv = zftv - dk-1[ aht e1v mj(wslpj) ] mj( trb ) 
(   55)       !!      take the horizontal divergence of the fluxes:
(   56)       !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
(   57)       !!      Add this trend to the general trend tra :
(   58)       !!         tra = tra + difft
(   59)       !!
(   60)       !!      'key_trdtra' defined: the trend is saved for diagnostics.
(   61)       !!
(   62)       !!      macro-tasked on horizontal slab (jk-loop).
(   63)       !!
(   64)       !! ** Action :
(   65)       !!         Update tra arrays with the before along level biharmonic
(   66)       !!      mixing trend.
(   67)       !!         Save in trtrd arrays the trends if 'key_trc_diatrd' defined
(   68)       !!
(   69)       !! History :
(   70)       !!        !  94-08  (G. Madec, M. Imbard)
(   71)       !!        !  97-05  (G. Madec)  split into traldf and trazdf
(   72)       !!   8.5  !  02-08  (G. Madec)  Free form, F90
(   73)       !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
(   74)       !!----------------------------------------------------------------------
(   75)       !! * Modules used
(   76)       USE oce_trc       , zftu => ua,  &  ! use ua as workspace
(   77)          &                zfsu => va      ! use va as workspace
(   78) 
(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   81) 
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk,jn          ! dummy loop indices
(   84)       INTEGER ::   iku, ikv               ! temporary integer
(   85)       REAL(wp) ::   &
(   86)          zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
(   87)          zmsku, zmskv, zbtr, ztra,     &
(   88)          ztagu, ztagv
(   89) 
(   90)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   91)          zdkt , zdk1t                     ! temporary workspace
(   92) 
(   93)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   & 
(   94)          zftv, zgtbu, zgtbv               ! temporary workspace
(   95) 
(   96) #if defined key_trcldf_eiv
(   97)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   98)          zftug, zftvg                    ! temporary workspace
(   99) 
(  100)       REAL(wp) ::   &
(  101)          zuwk, zvwk,                   &
(  102)          zuwk1, zvwk1,                 &
(  103)          zcg1,zcg2
(  104) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 3

(  105)       CHARACTER (len=22) :: charout
(  106)       !!----------------------------------------------------------------------
(  107) 
(  108)       IF( kt == nittrc000 ) THEN
(  109)          IF(lwp) WRITE(numout,*)
(  110)          IF(lwp) WRITE(numout,*) 'trc_ldf_iso_zps : iso neutral laplacian diffusion in '
(  111)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   z-coordinates with partial steps'
(  112) #if defined key_trcldf_eiv && defined key_diaeiv
(  113)          u_trc_eiv(:,:,:) = 0.e0
(  114)          v_trc_eiv(:,:,:) = 0.e0
(  115) #endif
(  116)       ENDIF
(  117) 
(  118)       DO jn = 1, jptra
(  119) 
(  120)          ztagu = 0.e0
(  121)          ztagv = 0.e0
(  122) 
(  123)          ! Horizontal passive tracer gradient 
(  124)          DO jk = 1, jpk
(  125)             DO jj = 1, jpj-1
(  126)                DO ji = 1, fs_jpim1   ! vector opt.
(  127)                   zgtbu(ji,jj,jk) = tmask(ji,jj,jk) * ( trb(ji+1,jj  ,jk,jn) - trb(ji,jj,jk,jn) )
(  128)                   zgtbv(ji,jj,jk) = tmask(ji,jj,jk) * ( trb(ji  ,jj+1,jk,jn) - trb(ji,jj,jk,jn) )
(  129)                END DO
(  130)             END DO
(  131)          END DO
(  132)          ! partial steps correction at the last level 
(  133)          DO jj = 1, jpj-1
(  134)             DO ji = 1, jpi-1
(  135)                ! last level
(  136)                iku = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
(  137)                ikv = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
(  138)                zgtbu(ji,jj,iku) = gtru(ji,jj,jn)               
(  139)                zgtbv(ji,jj,ikv) = gtrv(ji,jj,jn)               
(  140)             END DO
(  141)          END DO
(  142) 
(  143)          !                                                ! ===============
(  144)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  145)             !                                             ! ===============
(  146)             ! 1. Vertical tracer gradient at level jk and jk+1
(  147)             ! ------------------------------------------------
(  148)             ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
(  149) 
(  150)             zdk1t(:,:) = ( trb(:,:,jk,jn) - trb(:,:,jk+1,jn) ) * tmask(:,:,jk+1)
(  151) 
(  152)             IF( jk == 1 ) THEN
(  153)                zdkt(:,:) = zdk1t(:,:)
(  154)             ELSE
(  155)                zdkt(:,:) = ( trb(:,:,jk-1,jn) - trb(:,:,jk,jn) ) * tmask(:,:,jk)
(  156)             ENDIF
(  157) 
(  158) 
(  159)             ! 2. Horizontal fluxes
(  160)             ! --------------------
(  161) 
(  162)             DO jj = 1 , jpjm1






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 4

(  163)                DO ji = 1, fs_jpim1   ! vector opt.
(  164)                   zabe1 = ( fsahtru(ji,jj,jk) + ahtrb0 ) * e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
(  165)                   zabe2 = ( fsahtrv(ji,jj,jk) + ahtrb0 ) * e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
(  166) 
(  167)                   zmsku = 1. / MAX(  tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
(  168)                      + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )
(  169) 
(  170)                   zmskv = 1. / MAX(  tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
(  171)                      + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )
(  172) 
(  173)                   zcof1 = -fsahtru(ji,jj,jk) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
(  174)                   zcof2 = -fsahtrv(ji,jj,jk) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv
(  175) 
(  176)                   zftu(ji,jj,jk) = umask(ji,jj,jk) * (  zabe1 * zgtbu(ji,jj,jk)   &
(  177)                      &                                + zcof1 * (  zdkt (ji+1,jj) + zdk1t(ji,jj)      &
(  178)                      &                                           + zdk1t(ji+1,jj) + zdkt (ji,jj)  )  )
(  179)                   zftv(ji,jj,jk) = vmask(ji,jj,jk) * (  zabe2 * zgtbv(ji,jj,jk)   &
(  180)                      &                                + zcof2 * (  zdkt (ji,jj+1) + zdk1t(ji,jj)      &
(  181)                      &                                           + zdk1t(ji,jj+1) + zdkt (ji,jj)  )  )
(  182)                END DO
(  183)             END DO
(  184) 
(  185) # if defined key_trcldf_eiv
(  186)             ! ---------------------------------------!
(  187)             ! Eddy induced vertical advective fluxes !
(  188)             ! ---------------------------------------!
(  189)             DO jj = 1, jpjm1
(  190)                DO ji = 1, fs_jpim1   ! vector opt.
(  191)                   zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj  ,jk  ) ) * fsaeitru(ji,jj,jk  ) * umask(ji,jj,jk  )
(  192)                   zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj  ,jk+1) ) * fsaeitru(ji,jj,jk+1) * umask(ji,jj,jk+1)
(  193)                   zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji  ,jj+1,jk  ) ) * fsaeitrv(ji,jj,jk  ) * vmask(ji,jj,jk  )
(  194)                   zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji  ,jj+1,jk+1) ) * fsaeitrv(ji,jj,jk+1) * vmask(ji,jj,jk+1)
(  195) 
(  196)                   zcg1= -0.25 * e2u(ji,jj) * umask(ji,jj,jk) * ( zuwk-zuwk1 )
(  197)                   zcg2= -0.25 * e1v(ji,jj) * vmask(ji,jj,jk) * ( zvwk-zvwk1 )
(  198) 
(  199)                   zftug(ji,jj) = zcg1 * ( trb(ji+1,jj,jk,jn) + trb(ji,jj,jk,jn) )
(  200)                   zftvg(ji,jj) = zcg2 * ( trb(ji,jj+1,jk,jn) + trb(ji,jj,jk,jn) )
(  201) 
(  202)                   zftu(ji,jj,jk) = zftu(ji,jj,jk) + zftug(ji,jj)
(  203)                   zftv(ji,jj,jk) = zftv(ji,jj,jk) + zftvg(ji,jj)
(  204) 
(  205) #   if defined key_diaeiv
(  206)                   u_trc_eiv(ji,jj,jk) = -2. * zcg1 / ( e2u(ji,jj) * fse3u(ji,jj,jk) )
(  207)                   v_trc_eiv(ji,jj,jk) = -2. * zcg2 / ( e1v(ji,jj) * fse3v(ji,jj,jk) )
(  208) #   endif
(  209)                END DO
(  210)             END DO
(  211) # endif
(  212) 
(  213)             ! II.4 Second derivative (divergence) and add to the general trend
(  214)             ! ----------------------------------------------------------------
(  215) 
(  216)             DO jj = 2 , jpjm1
(  217)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  218)                   zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  219)                   ztra = zbtr * ( zftu(ji,jj,jk) - zftu(ji-1,jj,jk) + zftv(ji,jj,jk) - zftv(ji,jj-1,jk)  )
(  220)                   tra (ji,jj,jk,jn) = tra (ji,jj,jk,jn) + ztra






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 5

(  221) #if defined key_trc_diatrd
(  222)                   IF (luttrd(jn)) trtrd (ji,jj,jk,ikeep(jn),4) = zbtr * ( zftu(ji,jj,jk) - zftu(ji-1,  jj,jk) )
(  223)                   IF (luttrd(jn)) trtrd (ji,jj,jk,ikeep(jn),5) = zbtr * ( zftv(ji,jj,jk) - zftv(ji  ,jj-1,jk) )
(  224) #endif
(  225)                END DO
(  226)             END DO
(  227) #if defined key_trc_diatrd
(  228) #   if defined key_trcldf_eiv
(  229)             DO jj = 2 , jpjm1
(  230)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  231)                   zbtr= 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  232)                   ztagu = ( zftug(ji,jj) - zftug(ji-1,jj  ) ) * zbtr
(  233)                   ztagv = ( zftvg(ji,jj) - zftvg(ji  ,jj-1) ) * zbtr
(  234)                   IF (luttrd(jn)) trtrd (ji,jj,jk,ikeep(jn),4) = trtrd(ji,jj,jk,ikeep(jn),4) - ztagu
(  235)                   IF (luttrd(jn)) trtrd (ji,jj,jk,ikeep(jn),5) = trtrd(ji,jj,jk,ikeep(jn),5) - ztagv
(  236)                END DO
(  237)             END DO
(  238) #   endif
(  239) #endif
(  240)             !                                          ! ===============
(  241)          END DO                                        !   End of slab  
(  242)          !                                             ! ===============
(  243)       END DO
(  244) 
(  245)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  246)          WRITE(charout, FMT="('ldf - iso/zps')")
(  247)          CALL prt_ctl_trc_info(charout)
(  248)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  249)       ENDIF
(  250) 
(  251)    END SUBROUTINE trc_ldf_iso_zps
(  252) 
(  253) #else
(  254)    !!----------------------------------------------------------------------
(  255)    !!   default option :   Dummy code   NO rotation of the diffusive tensor
(  256)    !!----------------------------------------------------------------------
(  257) CONTAINS
(  258)    SUBROUTINE trc_ldf_iso_zps( kt )               ! Empty routine
(  259)       INTEGER, INTENT(in) :: kt
(  260) !      WRITE(*,*) 'trc_ldf_iso_zps: You should not have seen this print! error?', kt
(  261)    END SUBROUTINE trc_ldf_iso_zps
(  262) #endif
(  263) 
(  264)    !!==============================================================================
(  265) END MODULE trcldf_iso_zps



















PGF90 (Version     10.5)          09/23/2010  09:49:20      page 6

(    1) # 1 "trcldf_iso_zps.F90"
(    1) MODULE trcldf_iso_zps
(    2)    !!==============================================================================
(    3)    !!                   ***  MODULE  trcldf_iso_zps  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) # 254
(  254)    !!----------------------------------------------------------------------
(  255)    !!   default option :   Dummy code   NO rotation of the diffusive tensor
(  256)    !!----------------------------------------------------------------------
(  257) CONTAINS
(  258)    SUBROUTINE trc_ldf_iso_zps( kt )               ! Empty routine
(  259)       INTEGER, INTENT(in) :: kt
(  260) !      WRITE(*,*) 'trc_ldf_iso_zps: You should not have seen this print! error?', kt
(  261)    END SUBROUTINE trc_ldf_iso_zps
(  262) 
(  263) # 264
(  264)    !!==============================================================================
(  265) END MODULE trcldf_iso_zps
