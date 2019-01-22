


PGF90 (Version     10.5)          09/23/2010  09:49:37      page 1

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

Filename: dynldf_iso.F90

(    1) MODULE dynldf_iso
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  dynldf_iso  ***
(    4)    !! Ocean dynamics:  lateral viscosity trend
(    5)    !!======================================================================
(    6) #if defined key_ldfslp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                      slopes of the direction of mixing
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dyn_ldf_iso  : update the momentum trend with the horizontal part
(   11)    !!                  of the lateral diffusion using isopycnal or horizon-
(   12)    !!                  tal s-coordinate laplacian operator.
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE ldfdyn_oce      ! ocean dynamics lateral physics
(   18)    USE ldftra_oce      ! ocean tracer   lateral physics
(   19)    USE zdf_oce         ! ocean vertical physics
(   20)    USE trdmod          ! ocean dynamics trends 
(   21)    USE trdmod_oce      ! ocean variables trends
(   22)    USE ldfslp          ! iso-neutral slopes 
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC dyn_ldf_iso           ! called by step.F90
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34) #  include "ldfdyn_substitute.h90"
(   35) #  include "vectopt_loop_substitute.h90"
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_iso.F90,v 1.6 2005/09/02 15:45:24 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE dyn_ldf_iso( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                     ***  ROUTINE dyn_ldf_iso  ***






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 2

(   47)       !!                       
(   48)       !! ** Purpose :   Compute the before trend of the horizontal part of the
(   49)       !!      lateral momentum diffusion and add it to the general trend of
(   50)       !!      momentum equation.
(   51)       !!
(   52)       !! ** Method :
(   53)       !!        The horizontal component of the lateral diffusive trends on
(   54)       !!      momentum is provided by a 2nd order operator rotated along neu-
(   55)       !!      tral or geopotential surfaces (in s-coordinates).
(   56)       !!      It is computed using before fields (forward in time) and isopyc-
(   57)       !!      nal or geopotential slopes computed in routine ldfslp or inildf.
(   58)       !!      Here, u and v components are considered as 2 independent scalar
(   59)       !!      fields. Therefore, the property of splitting divergent and rota-
(   60)       !!      tional part of the flow of the standard, z-coordinate laplacian
(   61)       !!      momentum diffusion is lost.
(   62)       !!      horizontal fluxes associated with the rotated lateral mixing:
(   63)       !!      u-component:
(   64)       !!         ziut = ( ahmt + ahmb0 ) e2t * e3t / e1t  di[ ub ]
(   65)       !!               -      ahmt       e2t * mi-1(uslp) dk[ mi(mk(ub)) ]
(   66)       !!         zjuf = ( ahmf + ahmb0 ) e1f * e3f / e2f  dj[ ub ]
(   67)       !!               -      ahmf       e1f * mi(vslp)   dk[ mj(mk(ub)) ]
(   68)       !!      v-component:
(   69)       !!         zivf = ( ahmf + ahmb0 ) e2t * e3t / e1t  di[ vb ]
(   70)       !!               -      ahmf       e2t * mj(uslp)   dk[ mi(mk(vb)) ]
(   71)       !!         zjvt = ( ahmt + ahmb0 ) e1f * e3f / e2f  dj[ ub ]
(   72)       !!               -      ahmt       e1f * mj-1(vslp) dk[ mj(mk(vb)) ]
(   73)       !!      take the horizontal divergence of the fluxes:
(   74)       !!         diffu = 1/(e1u*e2u*e3u) {  di  [ ziut ] + dj-1[ zjuf ]  }
(   75)       !!         diffv = 1/(e1v*e2v*e3v) {  di-1[ zivf ] + dj  [ zjvt ]  }
(   76)       !!      Add this trend to the general trend (ua,va):
(   77)       !!         ua = ua + diffu
(   78)       !!      'key_trddyn' defined: the trends are saved for diagnostics.
(   79)       !!
(   80)       !! ** Action :
(   81)       !!        Update (ua,va) arrays with the before geopotential biharmonic
(   82)       !!      mixing trend.
(   83)       !!        Save in (uldftrd,vldftrd) arrays the trends if 'key_trddyn' defined
(   84)       !!
(   85)       !! History :
(   86)       !!   8.0  !  97-07  (G. Madec)  Original code
(   87)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   88)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   89)       !!----------------------------------------------------------------------
(   90)       !! * Modules used     
(   91)       USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
(   92)                          ztdva => sa      ! use sa as 3D workspace   
(   93)       !! * Arguments
(   94)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   95) 
(   96)       !! * Local declarations
(   97)       INTEGER  ::   ji, jj, jk            ! dummy loop indices
(   98)       REAL(wp) ::   &
(   99)          zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
(  100)          zmskt, zmskf, zbu, zbv,       &
(  101)          zuah, zvah
(  102)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  103)          ziut, zjuf, zjvt, zivf,        & ! temporary workspace
(  104)          zdku, zdk1u, zdkv, zdk1v






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 3

(  105)       !!----------------------------------------------------------------------
(  106) 
(  107)       IF( kt == nit000 ) THEN
(  108)          IF(lwp) WRITE(numout,*)
(  109)          IF(lwp) WRITE(numout,*) 'dyn_ldf_iso : iso-neutral laplacian diffusive operator or '
(  110)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate horizontal diffusive operator'
(  111)       ENDIF
(  112) 
(  113)       ! Save ua and va trends
(  114)       IF( l_trddyn )   THEN
(  115)          ztdua(:,:,:) = ua(:,:,:) 
(  116)          ztdva(:,:,:) = va(:,:,:) 
(  117)       ENDIF
(  118)       !                                                ! ===============
(  119)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  120)          !                                             ! ===============
(  121) 
(  122)          ! Vertical u- and v-shears at level jk and jk+1
(  123)          ! ---------------------------------------------
(  124)          ! surface boundary condition: zdku(jk=1)=zdku(jk=2)
(  125)          !                             zdkv(jk=1)=zdkv(jk=2)
(  126) 
(  127)          zdk1u(:,:) = ( ub(:,:,jk) -ub(:,:,jk+1) ) * umask(:,:,jk+1)
(  128)          zdk1v(:,:) = ( vb(:,:,jk) -vb(:,:,jk+1) ) * vmask(:,:,jk+1)
(  129) 
(  130)          IF( jk == 1 ) THEN
(  131)             zdku(:,:) = zdk1u(:,:)
(  132)             zdkv(:,:) = zdk1v(:,:)
(  133)          ELSE
(  134)             zdku(:,:) = ( ub(:,:,jk-1) - ub(:,:,jk) ) * umask(:,:,jk)
(  135)             zdkv(:,:) = ( vb(:,:,jk-1) - vb(:,:,jk) ) * vmask(:,:,jk)
(  136)          ENDIF
(  137) 
(  138)          !                               -----f-----
(  139)          ! Horizontal fluxes on U             |  
(  140)          ! --------------------===        t   u   t
(  141)          !                                    |  
(  142)          ! i-flux at t-point             -----f-----
(  143) 
(  144)          DO jj = 2, jpjm1
(  145)             DO ji = fs_2, jpi   ! vector opt.
(  146)                zabe1 = ( fsahmt(ji,jj,jk) + ahmb0 )   &
(  147) #if defined key_partial_steps
(  148)                      * e2t(ji,jj) * MIN( fse3u(ji,jj,jk), fse3u(ji-1, jj,jk) ) / e1t(ji,jj)
(  149) #else
(  150)                      * e2t(ji,jj) * fse3t(ji,jj,jk) / e1t(ji,jj)
(  151) #endif
(  152) 
(  153)                zmskt = 1./MAX(  umask(ji-1,jj,jk  )+umask(ji,jj,jk+1)   &
(  154)                               + umask(ji-1,jj,jk+1)+umask(ji,jj,jk  ), 1. )
(  155) 
(  156)                zcof1 = - aht0 * e2t(ji,jj) * zmskt   &
(  157)                      * 0.5  * ( uslp(ji-1,jj,jk) + uslp(ji,jj,jk) )
(  158) 
(  159)                ziut(ji,jj) = tmask(ji,jj,jk) *   &
(  160)                            (  zabe1 * ( ub(ji,jj,jk) - ub(ji-1,jj,jk) )   &
(  161)                             + zcof1 * ( zdku (ji,jj) + zdk1u(ji-1,jj)     &
(  162)                                        +zdk1u(ji,jj) + zdku (ji-1,jj) )  )






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 4

(  163)             END DO
(  164)          END DO
(  165) 
(  166)          ! j-flux at f-point
(  167)          DO jj = 1, jpjm1
(  168)             DO ji = 1, fs_jpim1   ! vector opt.
(  169)                zabe2 = ( fsahmf(ji,jj,jk) + ahmb0 )   &
(  170)                      * e1f(ji,jj) * fse3f(ji,jj,jk) / e2f(ji,jj)
(  171) 
(  172)                zmskf = 1./MAX(  umask(ji,jj+1,jk  )+umask(ji,jj,jk+1)   &
(  173)                               + umask(ji,jj+1,jk+1)+umask(ji,jj,jk  ), 1. )
(  174) 
(  175)                zcof2 = - aht0 * e1f(ji,jj) * zmskf   &
(  176)                      * 0.5  * ( vslp(ji+1,jj,jk) + vslp(ji,jj,jk) )
(  177) 
(  178)                zjuf(ji,jj) = fmask(ji,jj,jk) *   &
(  179)                            (  zabe2 * ( ub(ji,jj+1,jk) - ub(ji,jj,jk) )   &
(  180)                             + zcof2 * ( zdku (ji,jj+1) + zdk1u(ji,jj)     &
(  181)                                        +zdk1u(ji,jj+1) + zdku (ji,jj) )  )
(  182)             END DO
(  183)          END DO
(  184) 
(  185)          !                                |   t   |
(  186)          ! Horizontal fluxes on V         |       |
(  187)          ! --------------------===        f---v---f
(  188)          !                                |       |
(  189)          ! i-flux at f-point              |   t   |
(  190) 
(  191)          DO jj = 2, jpjm1
(  192)             DO ji = 1, fs_jpim1   ! vector opt.
(  193)                zabe1 = ( fsahmf(ji,jj,jk) + ahmb0 )   &
(  194)                      * e2f(ji,jj) * fse3f(ji,jj,jk) / e1f(ji,jj)
(  195) 
(  196)                zmskf = 1./MAX(  vmask(ji+1,jj,jk  )+vmask(ji,jj,jk+1)   &
(  197)                               + vmask(ji+1,jj,jk+1)+vmask(ji,jj,jk  ), 1. )
(  198) 
(  199)                zcof1 = - aht0 * e2f(ji,jj) * zmskf   &
(  200)                      * 0.5 * ( uslp(ji,jj+1,jk) + uslp(ji,jj,jk) )
(  201) 
(  202)                zivf(ji,jj) = fmask(ji,jj,jk) *   &
(  203)                            (  zabe1 * ( vb(ji+1,jj,jk) - vb(ji,jj,jk) )   &
(  204)                             + zcof1 * ( zdkv (ji,jj) + zdk1v(ji+1,jj)     &
(  205)                                        +zdk1v(ji,jj) + zdkv (ji+1,jj) )  )
(  206)             END DO
(  207)          END DO
(  208) 
(  209)          ! j-flux at t-point
(  210)          DO jj = 2, jpj
(  211)             DO ji = 1, fs_jpim1   ! vector opt.
(  212)                zabe2 = ( fsahmt(ji,jj,jk) + ahmb0 )   &
(  213) #if defined key_partial_steps
(  214)                      * e1t(ji,jj) * MIN( fse3v(ji,jj,jk), fse3v(ji, jj-1, jk) ) / e2t(ji,jj)
(  215) #else
(  216)                      * e1t(ji,jj) * fse3t(ji,jj,jk) / e2t(ji,jj)
(  217) #endif
(  218) 
(  219)                zmskt = 1./MAX(  vmask(ji,jj-1,jk  )+vmask(ji,jj,jk+1)   &
(  220)                               + vmask(ji,jj-1,jk+1)+vmask(ji,jj,jk  ), 1. )






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 5

(  221) 
(  222)                zcof2 = - aht0 * e1t(ji,jj) * zmskt   &
(  223)                      * 0.5 * ( vslp(ji,jj-1,jk) + vslp(ji,jj,jk) )
(  224) 
(  225)                zjvt(ji,jj) = tmask(ji,jj,jk) *   &
(  226)                            (  zabe2 * ( vb(ji,jj,jk) - vb(ji,jj-1,jk) )   &
(  227)                             + zcof2 * ( zdkv (ji,jj-1) + zdk1v(ji,jj)     &
(  228)                                        +zdk1v(ji,jj-1) + zdkv (ji,jj) )  )
(  229)             END DO
(  230)          END DO
(  231) 
(  232) 
(  233)          ! Second derivative (divergence) and add to the general trend
(  234)          ! -----------------------------------------------------------
(  235) 
(  236)          DO jj = 2, jpjm1
(  237)             DO ji = 2, jpim1          !! Question vectop possible??? !!bug
(  238)                ! volume elements
(  239)                zbu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
(  240)                zbv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
(  241)                ! horizontal component of isopycnal momentum diffusive trends
(  242)                zuah =( ziut (ji+1,jj) - ziut (ji,jj  ) +   &
(  243)                        zjuf (ji  ,jj) - zjuf (ji,jj-1)  ) / zbu
(  244)                zvah =( zivf (ji,jj  ) - zivf (ji-1,jj) +   &
(  245)                        zjvt (ji,jj+1) - zjvt (ji,jj  )  ) / zbv
(  246)                ! add the trends to the general trends
(  247)                ua (ji,jj,jk) = ua (ji,jj,jk) + zuah
(  248)                va (ji,jj,jk) = va (ji,jj,jk) + zvah
(  249)             END DO
(  250)          END DO
(  251)          !                                             ! ===============
(  252)       END DO                                           !   End of slab
(  253)       !                                                ! ===============
(  254) 
(  255)       ! save the lateral diffusion trends for diagnostic
(  256)       ! momentum trends will be saved in dynzdf_iso.F90
(  257)       IF( l_trddyn )   THEN
(  258)          uldftrd(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  259)          vldftrd(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  260)       ENDIF
(  261) 
(  262)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  263)          CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
(  264)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  265)       ENDIF
(  266) 
(  267)    END SUBROUTINE dyn_ldf_iso
(  268) 
(  269) # else
(  270)    !!----------------------------------------------------------------------
(  271)    !!   Dummy module                           NO rotation of mixing tensor
(  272)    !!----------------------------------------------------------------------
(  273)    USE in_out_manager
(  274) CONTAINS
(  275)    SUBROUTINE dyn_ldf_iso( kt )               ! Empty routine
(  276)       if(lwp) WRITE(numout,*) 'dyn_ldf_iso: You should not have seen this print! error?', kt
(  277)    END SUBROUTINE dyn_ldf_iso
(  278) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 6

(  279) 
(  280)    !!======================================================================
(  281) END MODULE dynldf_iso





























































PGF90 (Version     10.5)          09/23/2010  09:49:37      page 7

(    1) # 1 "dynldf_iso.F90"
(    1) MODULE dynldf_iso
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  dynldf_iso  ***
(    4)    !! Ocean dynamics:  lateral viscosity trend
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                      slopes of the direction of mixing
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dyn_ldf_iso  : update the momentum trend with the horizontal part
(   11)    !!                  of the lateral diffusion using isopycnal or horizon-
(   12)    !!                  tal s-coordinate laplacian operator.
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE ldfdyn_oce      ! ocean dynamics lateral physics
(   18)    USE ldftra_oce      ! ocean tracer   lateral physics
(   19)    USE zdf_oce         ! ocean vertical physics
(   20)    USE trdmod          ! ocean dynamics trends 
(   21)    USE trdmod_oce      ! ocean variables trends
(   22)    USE ldfslp          ! iso-neutral slopes 
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC dyn_ldf_iso           ! called by step.F90
(   31) 
(   32)    !! * Substitutions
(   33) # 1 "./domzgr_substitute.h90"
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
(   48) # 34 "dynldf_iso.F90"
(   34) # 1 "./ldfdyn_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  ldfdyn_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsahm., the lateral eddy viscosity coeff. 
(    5)    !!      with a constant, or 1D, or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 8

(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn_substitute.h90,v 1.3 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19
(   19) !   ' key_dynldf_c3d' :                  3D coefficient
(   20) # 35 "dynldf_iso.F90"
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
(   12) # 36 "dynldf_iso.F90"
(   36) # 36
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_iso.F90,v 1.6 2005/09/02 15:45:24 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE dyn_ldf_iso( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                     ***  ROUTINE dyn_ldf_iso  ***
(   47)       !!                       
(   48)       !! ** Purpose :   Compute the before trend of the horizontal part of the
(   49)       !!      lateral momentum diffusion and add it to the general trend of
(   50)       !!      momentum equation.
(   51)       !!
(   52)       !! ** Method :
(   53)       !!        The horizontal component of the lateral diffusive trends on
(   54)       !!      momentum is provided by a 2nd order operator rotated along neu-
(   55)       !!      tral or geopotential surfaces (in s-coordinates).
(   56)       !!      It is computed using before fields (forward in time) and isopyc-
(   57)       !!      nal or geopotential slopes computed in routine ldfslp or inildf.
(   58)       !!      Here, u and v components are considered as 2 independent scalar
(   59)       !!      fields. Therefore, the property of splitting divergent and rota-
(   60)       !!      tional part of the flow of the standard, z-coordinate laplacian
(   61)       !!      momentum diffusion is lost.
(   62)       !!      horizontal fluxes associated with the rotated lateral mixing:
(   63)       !!      u-component:
(   64)       !!         ziut = ( ahmt + ahmb0 ) e2t * e3t / e1t  di[ ub ]
(   65)       !!               -      ahmt       e2t * mi-1(uslp) dk[ mi(mk(ub)) ]
(   66)       !!         zjuf = ( ahmf + ahmb0 ) e1f * e3f / e2f  dj[ ub ]
(   67)       !!               -      ahmf       e1f * mi(vslp)   dk[ mj(mk(ub)) ]
(   68)       !!      v-component:
(   69)       !!         zivf = ( ahmf + ahmb0 ) e2t * e3t / e1t  di[ vb ]
(   70)       !!               -      ahmf       e2t * mj(uslp)   dk[ mi(mk(vb)) ]
(   71)       !!         zjvt = ( ahmt + ahmb0 ) e1f * e3f / e2f  dj[ ub ]
(   72)       !!               -      ahmt       e1f * mj-1(vslp) dk[ mj(mk(vb)) ]
(   73)       !!      take the horizontal divergence of the fluxes:






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 9

(   74)       !!         diffu = 1/(e1u*e2u*e3u) {  di  [ ziut ] + dj-1[ zjuf ]  }
(   75)       !!         diffv = 1/(e1v*e2v*e3v) {  di-1[ zivf ] + dj  [ zjvt ]  }
(   76)       !!      Add this trend to the general trend (ua,va):
(   77)       !!         ua = ua + diffu
(   78)       !!      'key_trddyn' defined: the trends are saved for diagnostics.
(   79)       !!
(   80)       !! ** Action :
(   81)       !!        Update (ua,va) arrays with the before geopotential biharmonic
(   82)       !!      mixing trend.
(   83)       !!        Save in (uldftrd,vldftrd) arrays the trends if 'key_trddyn' defined
(   84)       !!
(   85)       !! History :
(   86)       !!   8.0  !  97-07  (G. Madec)  Original code
(   87)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   88)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   89)       !!----------------------------------------------------------------------
(   90)       !! * Modules used     
(   91)       USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
(   92)                          ztdva => sa      ! use sa as 3D workspace   
(   93)       !! * Arguments
(   94)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   95) 
(   96)       !! * Local declarations
(   97)       INTEGER  ::   ji, jj, jk            ! dummy loop indices
(   98)       REAL(wp) ::   &
(   99)          zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
(  100)          zmskt, zmskf, zbu, zbv,       &
(  101)          zuah, zvah
(  102)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  103)          ziut, zjuf, zjvt, zivf,        & ! temporary workspace
(  104)          zdku, zdk1u, zdkv, zdk1v
(  105)       !!----------------------------------------------------------------------
(  106) 
(  107)       IF( kt == nit000 ) THEN
(  108)          IF(lwp) WRITE(numout,*)
(  109)          IF(lwp) WRITE(numout,*) 'dyn_ldf_iso : iso-neutral laplacian diffusive operator or '
(  110)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate horizontal diffusive operator'
(  111)       ENDIF
(  112) 
(  113)       ! Save ua and va trends
(  114)       IF( l_trddyn )   THEN
(  115)          ztdua(:,:,:) = ua(:,:,:) 
(  116)          ztdva(:,:,:) = va(:,:,:) 
(  117)       ENDIF
(  118)       !                                                ! ===============
(  119)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  120)          !                                             ! ===============
(  121) 
(  122)          ! Vertical u- and v-shears at level jk and jk+1
(  123)          ! ---------------------------------------------
(  124)          ! surface boundary condition: zdku(jk=1)=zdku(jk=2)
(  125)          !                             zdkv(jk=1)=zdkv(jk=2)
(  126) 
(  127)          zdk1u(:,:) = ( ub(:,:,jk) -ub(:,:,jk+1) ) * umask(:,:,jk+1)
(  128)          zdk1v(:,:) = ( vb(:,:,jk) -vb(:,:,jk+1) ) * vmask(:,:,jk+1)
(  129) 
(  130)          IF( jk == 1 ) THEN
(  131)             zdku(:,:) = zdk1u(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 10

(  132)             zdkv(:,:) = zdk1v(:,:)
(  133)          ELSE
(  134)             zdku(:,:) = ( ub(:,:,jk-1) - ub(:,:,jk) ) * umask(:,:,jk)
(  135)             zdkv(:,:) = ( vb(:,:,jk-1) - vb(:,:,jk) ) * vmask(:,:,jk)
(  136)          ENDIF
(  137) 
(  138)          !                               -----f-----
(  139)          ! Horizontal fluxes on U             |  
(  140)          ! --------------------===        t   u   t
(  141)          !                                    |  
(  142)          ! i-flux at t-point             -----f-----
(  143) 
(  144)          DO jj = 2, jpjm1
(  145)             DO ji = 2, jpi   ! vector opt.
(  146)                zabe1 = ( ahm1(ji,jj,jk) + ahmb0 )   &
(  150)                      * e2t(ji,jj) * e3t(jk) / e1t(ji,jj)
(  151) 
(  152) # 153
(  153)                zmskt = 1./MAX(  umask(ji-1,jj,jk  )+umask(ji,jj,jk+1)   &
(  154)                               + umask(ji-1,jj,jk+1)+umask(ji,jj,jk  ), 1. )
(  155) 
(  156)                zcof1 = - aht0 * e2t(ji,jj) * zmskt   &
(  157)                      * 0.5  * ( uslp(ji-1,jj,jk) + uslp(ji,jj,jk) )
(  158) 
(  159)                ziut(ji,jj) = tmask(ji,jj,jk) *   &
(  160)                            (  zabe1 * ( ub(ji,jj,jk) - ub(ji-1,jj,jk) )   &
(  161)                             + zcof1 * ( zdku (ji,jj) + zdk1u(ji-1,jj)     &
(  162)                                        +zdk1u(ji,jj) + zdku (ji-1,jj) )  )
(  163)             END DO
(  164)          END DO
(  165) 
(  166)          ! j-flux at f-point
(  167)          DO jj = 1, jpjm1
(  168)             DO ji = 1, jpim1   ! vector opt.
(  169)                zabe2 = ( ahm2(ji,jj,jk) + ahmb0 )   &
(  170)                      * e1f(ji,jj) * e3t(jk) / e2f(ji,jj)
(  171) 
(  172)                zmskf = 1./MAX(  umask(ji,jj+1,jk  )+umask(ji,jj,jk+1)   &
(  173)                               + umask(ji,jj+1,jk+1)+umask(ji,jj,jk  ), 1. )
(  174) 
(  175)                zcof2 = - aht0 * e1f(ji,jj) * zmskf   &
(  176)                      * 0.5  * ( vslp(ji+1,jj,jk) + vslp(ji,jj,jk) )
(  177) 
(  178)                zjuf(ji,jj) = fmask(ji,jj,jk) *   &
(  179)                            (  zabe2 * ( ub(ji,jj+1,jk) - ub(ji,jj,jk) )   &
(  180)                             + zcof2 * ( zdku (ji,jj+1) + zdk1u(ji,jj)     &
(  181)                                        +zdk1u(ji,jj+1) + zdku (ji,jj) )  )
(  182)             END DO
(  183)          END DO
(  184) 
(  185)          !                                |   t   |
(  186)          ! Horizontal fluxes on V         |       |
(  187)          ! --------------------===        f---v---f
(  188)          !                                |       |
(  189)          ! i-flux at f-point              |   t   |
(  190) 
(  191)          DO jj = 2, jpjm1
(  192)             DO ji = 1, jpim1   ! vector opt.






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 11

(  193)                zabe1 = ( ahm2(ji,jj,jk) + ahmb0 )   &
(  194)                      * e2f(ji,jj) * e3t(jk) / e1f(ji,jj)
(  195) 
(  196)                zmskf = 1./MAX(  vmask(ji+1,jj,jk  )+vmask(ji,jj,jk+1)   &
(  197)                               + vmask(ji+1,jj,jk+1)+vmask(ji,jj,jk  ), 1. )
(  198) 
(  199)                zcof1 = - aht0 * e2f(ji,jj) * zmskf   &
(  200)                      * 0.5 * ( uslp(ji,jj+1,jk) + uslp(ji,jj,jk) )
(  201) 
(  202)                zivf(ji,jj) = fmask(ji,jj,jk) *   &
(  203)                            (  zabe1 * ( vb(ji+1,jj,jk) - vb(ji,jj,jk) )   &
(  204)                             + zcof1 * ( zdkv (ji,jj) + zdk1v(ji+1,jj)     &
(  205)                                        +zdk1v(ji,jj) + zdkv (ji+1,jj) )  )
(  206)             END DO
(  207)          END DO
(  208) 
(  209)          ! j-flux at t-point
(  210)          DO jj = 2, jpj
(  211)             DO ji = 1, jpim1   ! vector opt.
(  212)                zabe2 = ( ahm1(ji,jj,jk) + ahmb0 )   &
(  216)                      * e1t(ji,jj) * e3t(jk) / e2t(ji,jj)
(  217) 
(  218) # 219
(  219)                zmskt = 1./MAX(  vmask(ji,jj-1,jk  )+vmask(ji,jj,jk+1)   &
(  220)                               + vmask(ji,jj-1,jk+1)+vmask(ji,jj,jk  ), 1. )
(  221) 
(  222)                zcof2 = - aht0 * e1t(ji,jj) * zmskt   &
(  223)                      * 0.5 * ( vslp(ji,jj-1,jk) + vslp(ji,jj,jk) )
(  224) 
(  225)                zjvt(ji,jj) = tmask(ji,jj,jk) *   &
(  226)                            (  zabe2 * ( vb(ji,jj,jk) - vb(ji,jj-1,jk) )   &
(  227)                             + zcof2 * ( zdkv (ji,jj-1) + zdk1v(ji,jj)     &
(  228)                                        +zdk1v(ji,jj-1) + zdkv (ji,jj) )  )
(  229)             END DO
(  230)          END DO
(  231) 
(  232) 
(  233)          ! Second derivative (divergence) and add to the general trend
(  234)          ! -----------------------------------------------------------
(  235) 
(  236)          DO jj = 2, jpjm1
(  237)             DO ji = 2, jpim1          !! Question vectop possible??? !!bug
(  238)                ! volume elements
(  239)                zbu = e1u(ji,jj) * e2u(ji,jj) * e3t(jk)
(  240)                zbv = e1v(ji,jj) * e2v(ji,jj) * e3t(jk)
(  241)                ! horizontal component of isopycnal momentum diffusive trends
(  242)                zuah =( ziut (ji+1,jj) - ziut (ji,jj  ) +   &
(  243)                        zjuf (ji  ,jj) - zjuf (ji,jj-1)  ) / zbu
(  244)                zvah =( zivf (ji,jj  ) - zivf (ji-1,jj) +   &
(  245)                        zjvt (ji,jj+1) - zjvt (ji,jj  )  ) / zbv
(  246)                ! add the trends to the general trends
(  247)                ua (ji,jj,jk) = ua (ji,jj,jk) + zuah
(  248)                va (ji,jj,jk) = va (ji,jj,jk) + zvah
(  249)             END DO
(  250)          END DO
(  251)          !                                             ! ===============
(  252)       END DO                                           !   End of slab
(  253)       !                                                ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 12

(  254) 
(  255)       ! save the lateral diffusion trends for diagnostic
(  256)       ! momentum trends will be saved in dynzdf_iso.F90
(  257)       IF( l_trddyn )   THEN
(  258)          uldftrd(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  259)          vldftrd(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  260)       ENDIF
(  261) 
(  262)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  263)          CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
(  264)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  265)       ENDIF
(  266) 
(  267)    END SUBROUTINE dyn_ldf_iso
(  268) 
(  269) 
(  270) # 280
(  280)    !!======================================================================
(  281) END MODULE dynldf_iso
