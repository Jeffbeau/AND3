


PGF90 (Version     10.5)          09/23/2010  09:49:46      page 1

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

Filename: dynzdf_iso.F90

(    1) MODULE dynzdf_iso
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  dynzdf_iso  ***
(    4)    !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
(    5)    !!==============================================================================
(    6) #if   defined key_ldfslp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                          rotation of the mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dyn_zdf_iso  : update the momentum trend with the vertical diffusion
(   11)    !!                  (vertical mixing + vertical component of lateral
(   12)    !!                  mixing) (rotated lateral operator case)
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE phycst          ! physical constants
(   18)    USE zdf_oce         ! ocean vertical physics
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE taumod          ! surface ocean stress
(   21)    USE trdmod          ! ocean dynamics trends 
(   22)    USE trdmod_oce      ! ocean variables trends
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC dyn_zdf_iso    ! called by step.F90
(   30) 
(   31)    !! * Substitutions
(   32) #  include "domzgr_substitute.h90"
(   33) #  include "vectopt_loop_substitute.h90"
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_iso.F90,v 1.6 2005/09/02 15:45:24 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42)    SUBROUTINE dyn_zdf_iso( kt )
(   43)       !!----------------------------------------------------------------------
(   44)       !!                  ***  ROUTINE dyn_zdf_iso  ***
(   45)       !!                   
(   46)       !! ** Purpose :






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 2

(   47)       !!         Compute the vertical momentum trend due to both vertical and
(   48)       !!      lateral mixing (only for second order lateral operator, for
(   49)       !!      fourth order it is already computed and add to the general trend
(   50)       !!      in dynldf.F) and the surface forcing, and add it to the general
(   51)       !!      trend of the momentum equations.
(   52)       !!
(   53)       !! ** Method :
(   54)       !!         The vertical component of the lateral diffusive trends is
(   55)       !!      provided by a 2nd order operator rotated along neural or geopo-
(   56)       !!      tential surfaces to which an eddy induced advection can be added
(   57)       !!      It is computed using before fields (forward in time) and isopyc-
(   58)       !!      nal or geopotential slopes computed in routine ldfslp.
(   59)       !!
(   60)       !!      First part: vertical trends associated with the lateral mixing
(   61)       !!      ==========  (excluding the vertical flux proportional to dk[U] )
(   62)       !!      vertical fluxes associated with the rotated lateral mixing:
(   63)       !!         zfuw =-ahm {  e2t*mi(wslpi) di[ mi(mk(ub)) ]
(   64)       !!                     + e1t*mj(wslpj) dj[ mj(mk(ub)) ]  }
(   65)       !!      update and save in zavt the vertical eddy viscosity coefficient:
(   66)       !!         avmu = avmu + mi(wslpi)^2 + mj(wslj)^2
(   67)       !!      take the horizontal divergence of the fluxes:
(   68)       !!         diffu = 1/(e1u*e2u*e3u) dk[ zfuw ]
(   69)       !!      Add this trend to the general trend (ta,sa):
(   70)       !!         ua = ua + difft
(   71)       !!
(   72)       !!      Second part: vertical trend associated with the vertical physics
(   73)       !!      ===========  (including the vertical flux proportional to dk[U]
(   74)       !!                    associated with the lateral mixing, through the
(   75)       !!                    update of avmu)
(   76)       !!      The vertical diffusion of momentum is given by:
(   77)       !!             diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
(   78)       !!      using a backward (implicit) time stepping.
(   79)       !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
(   80)       !!      Add this trend to the general trend ua :
(   81)       !!         ua = ua + dz( avmu dz(u) )
(   82)       !!
(   83)       !!      'key_trddyn' defined: trend saved for further diagnostics.
(   84)       !!
(   85)       !!      macro-tasked on vertical slab (jj-loop)
(   86)       !!
(   87)       !! ** Action : - Update (ua,va) arrays with the after vertical diffusive
(   88)       !!               mixing trend.
(   89)       !!             - Save the trends in (ztdua,ztdva) ('key_trddyn')
(   90)       !!
(   91)       !! History :
(   92)       !!        !  90-10  (B. Blanke)  Original code
(   93)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   94)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   95)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   96)       !!---------------------------------------------------------------------
(   97)       !! * Modules used
(   98)       USE ldfslp    , ONLY : wslpi, wslpj
(   99)       USE ldftra_oce, ONLY : aht0
(  100)       USE oce, ONLY :    ztdua => ta,        & ! use ta as 3D workspace   
(  101)                          ztdva => sa           ! use sa as 3D workspace   
(  102) 
(  103)       !! * Arguments
(  104)       INTEGER, INTENT( in ) ::   kt            ! ocean time-step index






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 3

(  105)       
(  106)       !! * Local declarations
(  107)       INTEGER ::   ji, jj, jk                  ! dummy loop indices
(  108)       INTEGER ::   &
(  109)          ikst, ikenm2, ikstp1,               & ! temporary integers
(  110)          ikbu, ikbum1 , ikbv, ikbvm1           !    "        "      
(  111)       REAL(wp) ::   &
(  112)          zrau0r, z2dt,                       & ! temporary scalars
(  113)          z2dtf, zcoef, zzws
(  114)       REAL(wp) ::   &
(  115)          zcoef0, zcoef3, zcoef4, zbu, zbv, zmkt, zmkf,   &
(  116)          zuav, zvav, zuwslpi, zuwslpj, zvwslpi, zvwslpj
(  117)       REAL(wp), DIMENSION(jpi,jpk) ::        &
(  118)          zwx, zwy, zwz,                      & ! workspace arrays
(  119)          zwd, zws, zwi, zwt,                 & !    "        "
(  120)          zfuw, zdiu, zdju, zdj1u,            & !    "        "
(  121)          zfvw, zdiv, zdjv, zdj1v
(  122)       REAL(wp), DIMENSION(jpi,jpj) ::        &
(  123)          ztsx, ztsy, ztbx, ztby                ! temporary workspace arrays
(  124)       !!----------------------------------------------------------------------
(  125) 
(  126)       IF( kt == nit000 ) THEN
(  127)          IF(lwp) WRITE(numout,*)
(  128)          IF(lwp) WRITE(numout,*) 'dyn_zdf_iso : vertical momentum diffusion isopycnal operator'
(  129)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(  130)       ENDIF
(  131) 
(  132)       ! 0. Local constant initialization
(  133)       ! --------------------------------
(  134)       
(  135)       ! inverse of the reference density
(  136)       zrau0r = 1. / rau0
(  137)       ! Leap-frog environnement
(  138)       z2dt = 2. * rdt
(  139)       ! workspace arrays
(  140)       ztsx(:,:)   = 0.e0
(  141)       ztsy(:,:)   = 0.e0 
(  142)       ztbx(:,:)   = 0.e0
(  143)       ztby(:,:)   = 0.e0
(  144)       ! Euler time stepping when starting from rest
(  145)       IF ( neuler == 0 .AND. kt == nit000 ) z2dt = rdt
(  146) 
(  147)       ! Save ua and va trends
(  148)       IF( l_trddyn )   THEN
(  149)          ztdua(:,:,:) = ua(:,:,:) 
(  150)          ztdva(:,:,:) = va(:,:,:) 
(  151)       ENDIF
(  152) 
(  153)       !                                                ! ===============
(  154)       DO jj = 2, jpjm1                                 !  Vertical slab
(  155)          !                                             ! ===============
(  156) 
(  157)   
(  158)          ! I. vertical trends associated with the lateral mixing
(  159)          ! =====================================================
(  160)          !  (excluding the vertical flux proportional to dk[t]
(  161) 
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 4

(  163)          ! I.1 horizontal momentum gradient
(  164)          ! --------------------------------
(  165) 
(  166)          DO jk = 1, jpk
(  167)             DO ji = 2, jpi
(  168)                ! i-gradient of u at jj
(  169)                zdiu (ji,jk) = tmask(ji,jj  ,jk) * ( ub(ji,jj  ,jk) - ub(ji-1,jj  ,jk) )
(  170)                ! j-gradient of u and v at jj
(  171)                zdju (ji,jk) = fmask(ji,jj  ,jk) * ( ub(ji,jj+1,jk) - ub(ji  ,jj  ,jk) )
(  172)                zdjv (ji,jk) = tmask(ji,jj  ,jk) * ( vb(ji,jj  ,jk) - vb(ji  ,jj-1,jk) )
(  173)                ! j-gradient of u and v at jj+1
(  174)                zdj1u(ji,jk) = fmask(ji,jj-1,jk) * ( ub(ji,jj  ,jk) - ub(ji  ,jj-1,jk) )
(  175)                zdj1v(ji,jk) = tmask(ji,jj+1,jk) * ( vb(ji,jj+1,jk) - vb(ji  ,jj  ,jk) )
(  176)             END DO
(  177)          END DO
(  178)          DO jk = 1, jpk
(  179)             DO ji = 1, jpim1
(  180)                ! i-gradient of v at jj
(  181)                zdiv (ji,jk) = fmask(ji,jj  ,jk) * ( vb(ji+1,jj,jk) - vb(ji  ,jj  ,jk) )
(  182)             END DO
(  183)          END DO
(  184) 
(  185) 
(  186)          ! I.2 Vertical fluxes
(  187)          ! -------------------
(  188) 
(  189)          ! Surface and bottom vertical fluxes set to zero
(  190)          DO ji = 1, jpi
(  191)             zfuw(ji, 1 ) = 0.e0
(  192)             zfvw(ji, 1 ) = 0.e0
(  193)             zfuw(ji,jpk) = 0.e0
(  194)             zfvw(ji,jpk) = 0.e0
(  195)          END DO
(  196) 
(  197)          ! interior (2=<jk=<jpk-1) on U field
(  198)          DO jk = 2, jpkm1
(  199)             DO ji = 2, jpim1
(  200)                zcoef0= 0.5 * aht0 * umask(ji,jj,jk)
(  201) 
(  202)                zuwslpi = zcoef0 * ( wslpi(ji+1,jj,jk) + wslpi(ji,jj,jk) )
(  203)                zuwslpj = zcoef0 * ( wslpj(ji+1,jj,jk) + wslpj(ji,jj,jk) )
(  204) 
(  205)                zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji+1,jj,jk-1)   &
(  206)                              + tmask(ji,jj,jk  )+tmask(ji+1,jj,jk  ), 1. )
(  207)                zmkf = 1./MAX(  fmask(ji,jj-1,jk-1)+fmask(ji,jj,jk-1)   &
(  208)                              + fmask(ji,jj-1,jk  )+fmask(ji,jj,jk  ), 1. )
(  209) 
(  210)                zcoef3 = - e2u(ji,jj) * zmkt * zuwslpi
(  211)                zcoef4 = - e1u(ji,jj) * zmkf * zuwslpj
(  212)                ! vertical flux on u field
(  213)                zfuw(ji,jk) = zcoef3 * ( zdiu (ji,jk-1) + zdiu (ji+1,jk-1)     &
(  214)                                        +zdiu (ji,jk  ) + zdiu (ji+1,jk  ) )   &
(  215)                            + zcoef4 * ( zdj1u(ji,jk-1) + zdju (ji  ,jk-1)     &
(  216)                                        +zdj1u(ji,jk  ) + zdju (ji  ,jk  ) )
(  217)                ! update avmu (add isopycnal vertical coefficient to avmu)
(  218)                avmu(ji,jj,jk) = avmu(ji,jj,jk) + ( zuwslpi * zuwslpi          &
(  219)                                                  + zuwslpj * zuwslpj ) / aht0
(  220)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 5

(  221)          END DO
(  222) 
(  223)          ! interior (2=<jk=<jpk-1) on V field
(  224)          DO jk = 2, jpkm1
(  225)             DO ji = 2, jpim1
(  226)                zcoef0= 0.5 * aht0 * vmask(ji,jj,jk)
(  227) 
(  228)                zvwslpi = zcoef0 * ( wslpi(ji,jj+1,jk) + wslpi(ji,jj,jk) )
(  229)                zvwslpj = zcoef0 * ( wslpj(ji,jj+1,jk) + wslpj(ji,jj,jk) )
(  230) 
(  231)                zmkf = 1./MAX(  fmask(ji-1,jj,jk-1)+fmask(ji,jj,jk-1)   &
(  232)                              + fmask(ji-1,jj,jk  )+fmask(ji,jj,jk  ), 1. )
(  233)                zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji,jj+1,jk-1)   &
(  234)                              + tmask(ji,jj,jk  )+tmask(ji,jj+1,jk  ), 1. )
(  235) 
(  236)                zcoef3 = - e2v(ji,jj) * zmkf * zvwslpi
(  237)                zcoef4 = - e1v(ji,jj) * zmkt * zvwslpj
(  238)                ! vertical flux on v field
(  239)                zfvw(ji,jk) = zcoef3 * ( zdiv (ji,jk-1) + zdiv (ji-1,jk-1)     &
(  240)                                        +zdiv (ji,jk  ) + zdiv (ji-1,jk  ) )   &
(  241)                            + zcoef4 * ( zdjv (ji,jk-1) + zdj1v(ji  ,jk-1)     &
(  242)                                        +zdjv (ji,jk  ) + zdj1v(ji  ,jk  ) )
(  243)                ! update avmv (add isopycnal vertical coefficient to avmv)
(  244)                avmv(ji,jj,jk) = avmv(ji,jj,jk) + ( zvwslpi * zvwslpi          &
(  245)                                                  + zvwslpj * zvwslpj ) / aht0
(  246)             END DO
(  247)          END DO
(  248) 
(  249) 
(  250)          ! I.3 Divergence of vertical fluxes added to the general tracer trend
(  251)          ! -------------------------------------------------------------------
(  252) 
(  253)          DO jk = 1, jpkm1
(  254)             DO ji = 2, jpim1
(  255)                ! volume elements
(  256)                zbu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
(  257)                zbv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
(  258)                ! part of the k-component of isopycnal momentum diffusive trends
(  259)                zuav = ( zfuw(ji,jk) - zfuw(ji,jk+1) ) / zbu
(  260)                zvav = ( zfvw(ji,jk) - zfvw(ji,jk+1) ) / zbv
(  261)                ! add the trends to the general trends
(  262)                ua(ji,jj,jk) = ua(ji,jj,jk) + zuav
(  263)                va(ji,jj,jk) = va(ji,jj,jk) + zvav
(  264)             END DO
(  265)          END DO
(  266)          !                                             ! ===============
(  267)       END DO                                           !   End of slab
(  268)       !                                                ! ===============
(  269)       IF( l_trddyn )   THEN
(  270)          ! save these trends in addition to the lateral diffusion one for diagnostics
(  271)          uldftrd(:,:,:) = uldftrd(:,:,:) + ua(:,:,:) - ztdua(:,:,:)
(  272)          vldftrd(:,:,:) = vldftrd(:,:,:) + va(:,:,:) - ztdva(:,:,:)
(  273) 
(  274)          ! save new trends ua and va
(  275)          ztdua(:,:,:) = ua(:,:,:) 
(  276)          ztdva(:,:,:) = va(:,:,:) 
(  277)       ENDIF
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 6

(  279)       !                                                ! ===============
(  280)       DO jj = 2, jpjm1                                 !  Vertical slab
(  281)          !                                             ! ===============
(  282)          ! 1. Vertical diffusion on u
(  283)          ! ---------------------------
(  284) 
(  285)          ! 1.0 Matrix and second member construction
(  286)          ! bottom boundary condition: only zws must be masked as avmu can take
(  287)          ! non zero value at the ocean bottom depending on the bottom friction
(  288)          ! used (see zdfmix.F)
(  289)          DO jk = 1, jpkm1
(  290)             DO ji = 2, jpim1
(  291)                zcoef = - z2dt / fse3u(ji,jj,jk)
(  292)                zwi(ji,jk) = zcoef * avmu(ji,jj,jk  ) / fse3uw(ji,jj,jk  )
(  293)                zzws       = zcoef * avmu(ji,jj,jk+1) / fse3uw(ji,jj,jk+1)
(  294)                zws(ji,jk) = zzws * umask(ji,jj,jk+1)
(  295)                zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
(  296)                zwy(ji,jk) = ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)
(  297)             END DO  
(  298)          END DO  
(  299) 
(  300)          ! 1.1 Surface boudary conditions
(  301)          DO ji = 2, jpim1
(  302)             z2dtf = z2dt / ( fse3u(ji,jj,1)*rau0 )
(  303)             zwi(ji,1) = 0.
(  304)             zwd(ji,1) = 1. - zws(ji,1)
(  305)             zwy(ji,1) = zwy(ji,1) + z2dtf * taux(ji,jj)
(  306)          END DO  
(  307) 
(  308)          ! 1.2 Matrix inversion starting from the first level
(  309)          ikst = 1
(  310) #include "zdf.matrixsolver.h90"
(  311) 
(  312)          ! 1.3 Normalization to obtain the general momentum trend ua
(  313)          DO jk = 1, jpkm1
(  314)             DO ji = 2, jpim1
(  315)                ua(ji,jj,jk) = ( zwx(ji,jk) - ub(ji,jj,jk) ) / z2dt
(  316)             END DO  
(  317)          END DO  
(  318) 
(  319)          ! 1.4 diagnose surface and bottom momentum fluxes
(  320)          DO ji = 2, jpim1
(  321)             ! save the surface forcing momentum fluxes
(  322)             ztsx(ji,jj) = taux(ji,jj) / ( fse3u(ji,jj,1)*rau0 )
(  323)             ! save bottom friction momentum fluxes
(  324)             ikbu   = min( mbathy(ji+1,jj), mbathy(ji,jj) )
(  325)             ikbum1 = max( ikbu-1, 1 )
(  326)             ztbx(ji,jj) = - avmu(ji,jj,ikbu) * zwx(ji,ikbum1)   &
(  327)                             / ( fse3u(ji,jj,ikbum1)*fse3uw(ji,jj,ikbu) )
(  328)             ! subtract surface forcing and bottom friction trend from vertical
(  329)             ! diffusive momentum trend
(  330)             ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
(  331)             ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
(  332)          END DO
(  333) 
(  334)          ! 2. Vertical diffusion on v
(  335)          ! ---------------------------
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 7

(  337)          ! 2.0 Matrix and second member construction
(  338)          ! bottom boundary condition: only zws must be masked as avmv can take
(  339)          ! non zero value at the ocean bottom depending on the bottom friction
(  340)          ! used (see zdfmix.F)
(  341)          DO jk = 1, jpkm1
(  342)             DO ji = 2, jpim1
(  343)                zcoef = -z2dt/fse3v(ji,jj,jk)
(  344)                zwi(ji,jk) = zcoef * avmv(ji,jj,jk  ) / fse3vw(ji,jj,jk  )
(  345)                zzws       = zcoef * avmv(ji,jj,jk+1) / fse3vw(ji,jj,jk+1)
(  346)                zws(ji,jk) =  zzws * vmask(ji,jj,jk+1)
(  347)                zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
(  348)                zwy(ji,jk) = vb(ji,jj,jk) + z2dt * va(ji,jj,jk)
(  349)             END DO  
(  350)          END DO  
(  351) 
(  352)          ! 2.1 Surface boudary conditions
(  353)          DO ji = 2, jpim1
(  354)             z2dtf = z2dt / ( fse3v(ji,jj,1)*rau0 )
(  355)             zwi(ji,1) = 0.e0
(  356)             zwd(ji,1) = 1. - zws(ji,1)
(  357)             zwy(ji,1) = zwy(ji,1) + z2dtf * tauy(ji,jj)
(  358)          END DO  
(  359) 
(  360)          ! 2.2 Matrix inversion starting from the first level
(  361)          ikst = 1
(  362) #include "zdf.matrixsolver.h90"
(  363) 
(  364)          ! 2.3 Normalization to obtain the general momentum trend va
(  365)          DO jk = 1, jpkm1
(  366)             DO ji = 2, jpim1
(  367)                va(ji,jj,jk) = ( zwx(ji,jk) - vb(ji,jj,jk) ) / z2dt
(  368)             END DO    
(  369)          END DO   
(  370) 
(  371)          ! 2.4 diagnose surface and bottom momentum fluxes
(  372)          DO ji = 2, jpim1
(  373)             ! save the surface forcing momentum fluxes
(  374)             ztsy(ji,jj) = tauy(ji,jj) / ( fse3v(ji,jj,1)*rau0 )
(  375)             ! save bottom friction momentum fluxes
(  376)             ikbv   = min( mbathy(ji,jj+1), mbathy(ji,jj) )
(  377)             ikbvm1 = max( ikbv-1, 1 )
(  378)             ztby(ji,jj) = - avmv(ji,jj,ikbv) * zwx(ji,ikbvm1)   &
(  379)                             / ( fse3v(ji,jj,ikbvm1)*fse3vw(ji,jj,ikbv) )
(  380)             ! subtract surface forcing and bottom friction trend from vertical
(  381)             ! diffusive momentum trend
(  382)             ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
(  383)             ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
(  384)          END DO
(  385)          !                                             ! ===============
(  386)       END DO                                           !   End of slab
(  387)       !                                                ! ===============
(  388) 
(  389)       ! save the vertical diffusive trends for diagnostic
(  390)       ! momentum trends
(  391)       IF( l_trddyn )  THEN 
(  392)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  393)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 8

(  395)          CALL trd_mod(uldftrd, vldftrd, jpdtdldf, 'DYN', kt)
(  396)          CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
(  397)          ztdua(:,:,:) = 0.e0
(  398)          ztdva(:,:,:) = 0.e0
(  399)          ztdua(:,:,1) = ztsx(:,:)
(  400)          ztdva(:,:,1) = ztsy(:,:)
(  401)          CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
(  402)          ztdua(:,:,:) = 0.e0
(  403)          ztdva(:,:,:) = 0.e0
(  404)          ztdua(:,:,1) = ztbx(:,:)
(  405)          ztdva(:,:,1) = ztby(:,:)
(  406)          CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
(  407)       ENDIF
(  408) 
(  409)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  410)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
(  411)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  412)       ENDIF
(  413) 
(  414)    END SUBROUTINE dyn_zdf_iso
(  415) 
(  416) #else
(  417)    !!----------------------------------------------------------------------
(  418)    !!   Dummy module                       NO rotation of the mixing tensor
(  419)    !!----------------------------------------------------------------------
(  420)    USE in_out_manager
(  421) CONTAINS
(  422)    SUBROUTINE dyn_zdf_iso( kt )                        ! Dummy routine
(  423)       if(lwp) WRITE(numout,*) 'dyn_zdf_iso: You should not have seen this print! error?', kt
(  424)    END SUBROUTINE dyn_zdf_iso
(  425) #endif
(  426) 
(  427)    !!==============================================================================
(  428) END MODULE dynzdf_iso






























PGF90 (Version     10.5)          09/23/2010  09:49:46      page 9

(    1) # 1 "dynzdf_iso.F90"
(    1) MODULE dynzdf_iso
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  dynzdf_iso  ***
(    4)    !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
(    5)    !!==============================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                          rotation of the mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dyn_zdf_iso  : update the momentum trend with the vertical diffusion
(   11)    !!                  (vertical mixing + vertical component of lateral
(   12)    !!                  mixing) (rotated lateral operator case)
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE phycst          ! physical constants
(   18)    USE zdf_oce         ! ocean vertical physics
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE taumod          ! surface ocean stress
(   21)    USE trdmod          ! ocean dynamics trends 
(   22)    USE trdmod_oce      ! ocean variables trends
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC dyn_zdf_iso    ! called by step.F90
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
(   48) # 33 "dynzdf_iso.F90"
(   33) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 10

(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 34 "dynzdf_iso.F90"
(   34) # 34
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_iso.F90,v 1.6 2005/09/02 15:45:24 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42)    SUBROUTINE dyn_zdf_iso( kt )
(   43)       !!----------------------------------------------------------------------
(   44)       !!                  ***  ROUTINE dyn_zdf_iso  ***
(   45)       !!                   
(   46)       !! ** Purpose :
(   47)       !!         Compute the vertical momentum trend due to both vertical and
(   48)       !!      lateral mixing (only for second order lateral operator, for
(   49)       !!      fourth order it is already computed and add to the general trend
(   50)       !!      in dynldf.F) and the surface forcing, and add it to the general
(   51)       !!      trend of the momentum equations.
(   52)       !!
(   53)       !! ** Method :
(   54)       !!         The vertical component of the lateral diffusive trends is
(   55)       !!      provided by a 2nd order operator rotated along neural or geopo-
(   56)       !!      tential surfaces to which an eddy induced advection can be added
(   57)       !!      It is computed using before fields (forward in time) and isopyc-
(   58)       !!      nal or geopotential slopes computed in routine ldfslp.
(   59)       !!
(   60)       !!      First part: vertical trends associated with the lateral mixing
(   61)       !!      ==========  (excluding the vertical flux proportional to dk[U] )
(   62)       !!      vertical fluxes associated with the rotated lateral mixing:
(   63)       !!         zfuw =-ahm {  e2t*mi(wslpi) di[ mi(mk(ub)) ]
(   64)       !!                     + e1t*mj(wslpj) dj[ mj(mk(ub)) ]  }
(   65)       !!      update and save in zavt the vertical eddy viscosity coefficient:
(   66)       !!         avmu = avmu + mi(wslpi)^2 + mj(wslj)^2
(   67)       !!      take the horizontal divergence of the fluxes:
(   68)       !!         diffu = 1/(e1u*e2u*e3u) dk[ zfuw ]
(   69)       !!      Add this trend to the general trend (ta,sa):
(   70)       !!         ua = ua + difft
(   71)       !!
(   72)       !!      Second part: vertical trend associated with the vertical physics
(   73)       !!      ===========  (including the vertical flux proportional to dk[U]
(   74)       !!                    associated with the lateral mixing, through the
(   75)       !!                    update of avmu)
(   76)       !!      The vertical diffusion of momentum is given by:
(   77)       !!             diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
(   78)       !!      using a backward (implicit) time stepping.
(   79)       !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
(   80)       !!      Add this trend to the general trend ua :
(   81)       !!         ua = ua + dz( avmu dz(u) )
(   82)       !!
(   83)       !!      'key_trddyn' defined: trend saved for further diagnostics.
(   84)       !!
(   85)       !!      macro-tasked on vertical slab (jj-loop)
(   86)       !!
(   87)       !! ** Action : - Update (ua,va) arrays with the after vertical diffusive






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 11

(   88)       !!               mixing trend.
(   89)       !!             - Save the trends in (ztdua,ztdva) ('key_trddyn')
(   90)       !!
(   91)       !! History :
(   92)       !!        !  90-10  (B. Blanke)  Original code
(   93)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   94)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   95)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   96)       !!---------------------------------------------------------------------
(   97)       !! * Modules used
(   98)       USE ldfslp    , ONLY : wslpi, wslpj
(   99)       USE ldftra_oce, ONLY : aht0
(  100)       USE oce, ONLY :    ztdua => ta,        & ! use ta as 3D workspace   
(  101)                          ztdva => sa           ! use sa as 3D workspace   
(  102) 
(  103)       !! * Arguments
(  104)       INTEGER, INTENT( in ) ::   kt            ! ocean time-step index
(  105)       
(  106)       !! * Local declarations
(  107)       INTEGER ::   ji, jj, jk                  ! dummy loop indices
(  108)       INTEGER ::   &
(  109)          ikst, ikenm2, ikstp1,               & ! temporary integers
(  110)          ikbu, ikbum1 , ikbv, ikbvm1           !    "        "      
(  111)       REAL(wp) ::   &
(  112)          zrau0r, z2dt,                       & ! temporary scalars
(  113)          z2dtf, zcoef, zzws
(  114)       REAL(wp) ::   &
(  115)          zcoef0, zcoef3, zcoef4, zbu, zbv, zmkt, zmkf,   &
(  116)          zuav, zvav, zuwslpi, zuwslpj, zvwslpi, zvwslpj
(  117)       REAL(wp), DIMENSION(jpi,jpk) ::        &
(  118)          zwx, zwy, zwz,                      & ! workspace arrays
(  119)          zwd, zws, zwi, zwt,                 & !    "        "
(  120)          zfuw, zdiu, zdju, zdj1u,            & !    "        "
(  121)          zfvw, zdiv, zdjv, zdj1v
(  122)       REAL(wp), DIMENSION(jpi,jpj) ::        &
(  123)          ztsx, ztsy, ztbx, ztby                ! temporary workspace arrays
(  124)       !!----------------------------------------------------------------------
(  125) 
(  126)       IF( kt == nit000 ) THEN
(  127)          IF(lwp) WRITE(numout,*)
(  128)          IF(lwp) WRITE(numout,*) 'dyn_zdf_iso : vertical momentum diffusion isopycnal operator'
(  129)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(  130)       ENDIF
(  131) 
(  132)       ! 0. Local constant initialization
(  133)       ! --------------------------------
(  134)       
(  135)       ! inverse of the reference density
(  136)       zrau0r = 1. / rau0
(  137)       ! Leap-frog environnement
(  138)       z2dt = 2. * rdt
(  139)       ! workspace arrays
(  140)       ztsx(:,:)   = 0.e0
(  141)       ztsy(:,:)   = 0.e0 
(  142)       ztbx(:,:)   = 0.e0
(  143)       ztby(:,:)   = 0.e0
(  144)       ! Euler time stepping when starting from rest
(  145)       IF ( neuler == 0 .AND. kt == nit000 ) z2dt = rdt






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 12

(  146) 
(  147)       ! Save ua and va trends
(  148)       IF( l_trddyn )   THEN
(  149)          ztdua(:,:,:) = ua(:,:,:) 
(  150)          ztdva(:,:,:) = va(:,:,:) 
(  151)       ENDIF
(  152) 
(  153)       !                                                ! ===============
(  154)       DO jj = 2, jpjm1                                 !  Vertical slab
(  155)          !                                             ! ===============
(  156) 
(  157)   
(  158)          ! I. vertical trends associated with the lateral mixing
(  159)          ! =====================================================
(  160)          !  (excluding the vertical flux proportional to dk[t]
(  161) 
(  162) 
(  163)          ! I.1 horizontal momentum gradient
(  164)          ! --------------------------------
(  165) 
(  166)          DO jk = 1, jpk
(  167)             DO ji = 2, jpi
(  168)                ! i-gradient of u at jj
(  169)                zdiu (ji,jk) = tmask(ji,jj  ,jk) * ( ub(ji,jj  ,jk) - ub(ji-1,jj  ,jk) )
(  170)                ! j-gradient of u and v at jj
(  171)                zdju (ji,jk) = fmask(ji,jj  ,jk) * ( ub(ji,jj+1,jk) - ub(ji  ,jj  ,jk) )
(  172)                zdjv (ji,jk) = tmask(ji,jj  ,jk) * ( vb(ji,jj  ,jk) - vb(ji  ,jj-1,jk) )
(  173)                ! j-gradient of u and v at jj+1
(  174)                zdj1u(ji,jk) = fmask(ji,jj-1,jk) * ( ub(ji,jj  ,jk) - ub(ji  ,jj-1,jk) )
(  175)                zdj1v(ji,jk) = tmask(ji,jj+1,jk) * ( vb(ji,jj+1,jk) - vb(ji  ,jj  ,jk) )
(  176)             END DO
(  177)          END DO
(  178)          DO jk = 1, jpk
(  179)             DO ji = 1, jpim1
(  180)                ! i-gradient of v at jj
(  181)                zdiv (ji,jk) = fmask(ji,jj  ,jk) * ( vb(ji+1,jj,jk) - vb(ji  ,jj  ,jk) )
(  182)             END DO
(  183)          END DO
(  184) 
(  185) 
(  186)          ! I.2 Vertical fluxes
(  187)          ! -------------------
(  188) 
(  189)          ! Surface and bottom vertical fluxes set to zero
(  190)          DO ji = 1, jpi
(  191)             zfuw(ji, 1 ) = 0.e0
(  192)             zfvw(ji, 1 ) = 0.e0
(  193)             zfuw(ji,jpk) = 0.e0
(  194)             zfvw(ji,jpk) = 0.e0
(  195)          END DO
(  196) 
(  197)          ! interior (2=<jk=<jpk-1) on U field
(  198)          DO jk = 2, jpkm1
(  199)             DO ji = 2, jpim1
(  200)                zcoef0= 0.5 * aht0 * umask(ji,jj,jk)
(  201) 
(  202)                zuwslpi = zcoef0 * ( wslpi(ji+1,jj,jk) + wslpi(ji,jj,jk) )
(  203)                zuwslpj = zcoef0 * ( wslpj(ji+1,jj,jk) + wslpj(ji,jj,jk) )






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 13

(  204) 
(  205)                zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji+1,jj,jk-1)   &
(  206)                              + tmask(ji,jj,jk  )+tmask(ji+1,jj,jk  ), 1. )
(  207)                zmkf = 1./MAX(  fmask(ji,jj-1,jk-1)+fmask(ji,jj,jk-1)   &
(  208)                              + fmask(ji,jj-1,jk  )+fmask(ji,jj,jk  ), 1. )
(  209) 
(  210)                zcoef3 = - e2u(ji,jj) * zmkt * zuwslpi
(  211)                zcoef4 = - e1u(ji,jj) * zmkf * zuwslpj
(  212)                ! vertical flux on u field
(  213)                zfuw(ji,jk) = zcoef3 * ( zdiu (ji,jk-1) + zdiu (ji+1,jk-1)     &
(  214)                                        +zdiu (ji,jk  ) + zdiu (ji+1,jk  ) )   &
(  215)                            + zcoef4 * ( zdj1u(ji,jk-1) + zdju (ji  ,jk-1)     &
(  216)                                        +zdj1u(ji,jk  ) + zdju (ji  ,jk  ) )
(  217)                ! update avmu (add isopycnal vertical coefficient to avmu)
(  218)                avmu(ji,jj,jk) = avmu(ji,jj,jk) + ( zuwslpi * zuwslpi          &
(  219)                                                  + zuwslpj * zuwslpj ) / aht0
(  220)             END DO
(  221)          END DO
(  222) 
(  223)          ! interior (2=<jk=<jpk-1) on V field
(  224)          DO jk = 2, jpkm1
(  225)             DO ji = 2, jpim1
(  226)                zcoef0= 0.5 * aht0 * vmask(ji,jj,jk)
(  227) 
(  228)                zvwslpi = zcoef0 * ( wslpi(ji,jj+1,jk) + wslpi(ji,jj,jk) )
(  229)                zvwslpj = zcoef0 * ( wslpj(ji,jj+1,jk) + wslpj(ji,jj,jk) )
(  230) 
(  231)                zmkf = 1./MAX(  fmask(ji-1,jj,jk-1)+fmask(ji,jj,jk-1)   &
(  232)                              + fmask(ji-1,jj,jk  )+fmask(ji,jj,jk  ), 1. )
(  233)                zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji,jj+1,jk-1)   &
(  234)                              + tmask(ji,jj,jk  )+tmask(ji,jj+1,jk  ), 1. )
(  235) 
(  236)                zcoef3 = - e2v(ji,jj) * zmkf * zvwslpi
(  237)                zcoef4 = - e1v(ji,jj) * zmkt * zvwslpj
(  238)                ! vertical flux on v field
(  239)                zfvw(ji,jk) = zcoef3 * ( zdiv (ji,jk-1) + zdiv (ji-1,jk-1)     &
(  240)                                        +zdiv (ji,jk  ) + zdiv (ji-1,jk  ) )   &
(  241)                            + zcoef4 * ( zdjv (ji,jk-1) + zdj1v(ji  ,jk-1)     &
(  242)                                        +zdjv (ji,jk  ) + zdj1v(ji  ,jk  ) )
(  243)                ! update avmv (add isopycnal vertical coefficient to avmv)
(  244)                avmv(ji,jj,jk) = avmv(ji,jj,jk) + ( zvwslpi * zvwslpi          &
(  245)                                                  + zvwslpj * zvwslpj ) / aht0
(  246)             END DO
(  247)          END DO
(  248) 
(  249) 
(  250)          ! I.3 Divergence of vertical fluxes added to the general tracer trend
(  251)          ! -------------------------------------------------------------------
(  252) 
(  253)          DO jk = 1, jpkm1
(  254)             DO ji = 2, jpim1
(  255)                ! volume elements
(  256)                zbu = e1u(ji,jj) * e2u(ji,jj) * e3t(jk)
(  257)                zbv = e1v(ji,jj) * e2v(ji,jj) * e3t(jk)
(  258)                ! part of the k-component of isopycnal momentum diffusive trends
(  259)                zuav = ( zfuw(ji,jk) - zfuw(ji,jk+1) ) / zbu
(  260)                zvav = ( zfvw(ji,jk) - zfvw(ji,jk+1) ) / zbv
(  261)                ! add the trends to the general trends






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 14

(  262)                ua(ji,jj,jk) = ua(ji,jj,jk) + zuav
(  263)                va(ji,jj,jk) = va(ji,jj,jk) + zvav
(  264)             END DO
(  265)          END DO
(  266)          !                                             ! ===============
(  267)       END DO                                           !   End of slab
(  268)       !                                                ! ===============
(  269)       IF( l_trddyn )   THEN
(  270)          ! save these trends in addition to the lateral diffusion one for diagnostics
(  271)          uldftrd(:,:,:) = uldftrd(:,:,:) + ua(:,:,:) - ztdua(:,:,:)
(  272)          vldftrd(:,:,:) = vldftrd(:,:,:) + va(:,:,:) - ztdva(:,:,:)
(  273) 
(  274)          ! save new trends ua and va
(  275)          ztdua(:,:,:) = ua(:,:,:) 
(  276)          ztdva(:,:,:) = va(:,:,:) 
(  277)       ENDIF
(  278) 
(  279)       !                                                ! ===============
(  280)       DO jj = 2, jpjm1                                 !  Vertical slab
(  281)          !                                             ! ===============
(  282)          ! 1. Vertical diffusion on u
(  283)          ! ---------------------------
(  284) 
(  285)          ! 1.0 Matrix and second member construction
(  286)          ! bottom boundary condition: only zws must be masked as avmu can take
(  287)          ! non zero value at the ocean bottom depending on the bottom friction
(  288)          ! used (see zdfmix.F)
(  289)          DO jk = 1, jpkm1
(  290)             DO ji = 2, jpim1
(  291)                zcoef = - z2dt / e3t(jk)
(  292)                zwi(ji,jk) = zcoef * avmu(ji,jj,jk  ) / e3w(jk  )
(  293)                zzws       = zcoef * avmu(ji,jj,jk+1) / e3w(jk+1)
(  294)                zws(ji,jk) = zzws * umask(ji,jj,jk+1)
(  295)                zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
(  296)                zwy(ji,jk) = ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)
(  297)             END DO  
(  298)          END DO  
(  299) 
(  300)          ! 1.1 Surface boudary conditions
(  301)          DO ji = 2, jpim1
(  302)             z2dtf = z2dt / ( e3t(1)*rau0 )
(  303)             zwi(ji,1) = 0.
(  304)             zwd(ji,1) = 1. - zws(ji,1)
(  305)             zwy(ji,1) = zwy(ji,1) + z2dtf * taux(ji,jj)
(  306)          END DO  
(  307) 
(  308)          ! 1.2 Matrix inversion starting from the first level
(  309)          ikst = 1
(  310) # 1 "./zdf.matrixsolver.h90"
(    1) !! 
(    2) !!         ZDF.MATRIXSOLVER
(    3) !!       ********************
(    4) !! 
(    5) !! Matrix inversion
(    6) !!----------------------------------------------------------------------
(    7) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(    8) !
(    9) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 15

(   10) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(   11) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(   12) !        (        ...               )( ...  ) ( ...  )
(   13) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(   14) !
(   15) !   m is decomposed in the product of an upper and lower triangular
(   16) !   matrix
(   17) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(   18) !   The second member is in 2d array zwy
(   19) !   The solution is in 2d array zwx
(   20) !   The 2d arry zwt and zwz are work space arrays
(   21) !
(   22) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(   23) !   the resolution of tke matrix where surface tke value is prescribed
(   24) !   so that ikstrt=2.
(   25) !!----------------------------------------------------------------------
(   26) !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   27) !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf.matrixsolver.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   28) !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   29) !!----------------------------------------------------------------------
(   30) 
(   31)       ikstp1 = ikst + 1
(   32)       ikenm2 = jpk - 2
(   33)       DO ji = 2, jpim1
(   34)          zwt(ji,ikst) = zwd(ji,ikst)
(   35)       END DO
(   36)       DO jk = ikstp1, jpkm1
(   37)          DO ji = 2, jpim1
(   38)             zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(   39)          END DO
(   40)       END DO
(   41)       DO ji = 2, jpim1
(   42)          zwz(ji,ikst) = zwy(ji,ikst)
(   43)       END DO
(   44)       DO jk = ikstp1, jpkm1
(   45)          DO ji = 2, jpim1
(   46)             zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(   47)          END DO
(   48)       END DO
(   49)       DO ji = 2, jpim1
(   50)          zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(   51)       END DO
(   52)       DO jk = ikenm2, ikst, -1
(   53)          DO ji = 2, jpim1
(   54)             zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(   55)          END DO
(   56)       END DO
(   57)       
(   58) # 311 "dynzdf_iso.F90"
(  311) 
(  312) # 312
(  312)          ! 1.3 Normalization to obtain the general momentum trend ua
(  313)          DO jk = 1, jpkm1
(  314)             DO ji = 2, jpim1
(  315)                ua(ji,jj,jk) = ( zwx(ji,jk) - ub(ji,jj,jk) ) / z2dt
(  316)             END DO  
(  317)          END DO  
(  318) 






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 16

(  319)          ! 1.4 diagnose surface and bottom momentum fluxes
(  320)          DO ji = 2, jpim1
(  321)             ! save the surface forcing momentum fluxes
(  322)             ztsx(ji,jj) = taux(ji,jj) / ( e3t(1)*rau0 )
(  323)             ! save bottom friction momentum fluxes
(  324)             ikbu   = min( mbathy(ji+1,jj), mbathy(ji,jj) )
(  325)             ikbum1 = max( ikbu-1, 1 )
(  326)             ztbx(ji,jj) = - avmu(ji,jj,ikbu) * zwx(ji,ikbum1)   &
(  327)                             / ( e3t(ikbum1)*e3w(ikbu) )
(  328)             ! subtract surface forcing and bottom friction trend from vertical
(  329)             ! diffusive momentum trend
(  330)             ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
(  331)             ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
(  332)          END DO
(  333) 
(  334)          ! 2. Vertical diffusion on v
(  335)          ! ---------------------------
(  336) 
(  337)          ! 2.0 Matrix and second member construction
(  338)          ! bottom boundary condition: only zws must be masked as avmv can take
(  339)          ! non zero value at the ocean bottom depending on the bottom friction
(  340)          ! used (see zdfmix.F)
(  341)          DO jk = 1, jpkm1
(  342)             DO ji = 2, jpim1
(  343)                zcoef = -z2dt/e3t(jk)
(  344)                zwi(ji,jk) = zcoef * avmv(ji,jj,jk  ) / e3w(jk  )
(  345)                zzws       = zcoef * avmv(ji,jj,jk+1) / e3w(jk+1)
(  346)                zws(ji,jk) =  zzws * vmask(ji,jj,jk+1)
(  347)                zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
(  348)                zwy(ji,jk) = vb(ji,jj,jk) + z2dt * va(ji,jj,jk)
(  349)             END DO  
(  350)          END DO  
(  351) 
(  352)          ! 2.1 Surface boudary conditions
(  353)          DO ji = 2, jpim1
(  354)             z2dtf = z2dt / ( e3t(1)*rau0 )
(  355)             zwi(ji,1) = 0.e0
(  356)             zwd(ji,1) = 1. - zws(ji,1)
(  357)             zwy(ji,1) = zwy(ji,1) + z2dtf * tauy(ji,jj)
(  358)          END DO  
(  359) 
(  360)          ! 2.2 Matrix inversion starting from the first level
(  361)          ikst = 1
(  362) # 1 "./zdf.matrixsolver.h90"
(    1) !! 
(    2) !!         ZDF.MATRIXSOLVER
(    3) !!       ********************
(    4) !! 
(    5) !! Matrix inversion
(    6) !!----------------------------------------------------------------------
(    7) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(    8) !
(    9) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(   10) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(   11) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(   12) !        (        ...               )( ...  ) ( ...  )
(   13) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(   14) !






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 17

(   15) !   m is decomposed in the product of an upper and lower triangular
(   16) !   matrix
(   17) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(   18) !   The second member is in 2d array zwy
(   19) !   The solution is in 2d array zwx
(   20) !   The 2d arry zwt and zwz are work space arrays
(   21) !
(   22) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(   23) !   the resolution of tke matrix where surface tke value is prescribed
(   24) !   so that ikstrt=2.
(   25) !!----------------------------------------------------------------------
(   26) !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   27) !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf.matrixsolver.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   28) !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   29) !!----------------------------------------------------------------------
(   30) 
(   31)       ikstp1 = ikst + 1
(   32)       ikenm2 = jpk - 2
(   33)       DO ji = 2, jpim1
(   34)          zwt(ji,ikst) = zwd(ji,ikst)
(   35)       END DO
(   36)       DO jk = ikstp1, jpkm1
(   37)          DO ji = 2, jpim1
(   38)             zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(   39)          END DO
(   40)       END DO
(   41)       DO ji = 2, jpim1
(   42)          zwz(ji,ikst) = zwy(ji,ikst)
(   43)       END DO
(   44)       DO jk = ikstp1, jpkm1
(   45)          DO ji = 2, jpim1
(   46)             zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(   47)          END DO
(   48)       END DO
(   49)       DO ji = 2, jpim1
(   50)          zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(   51)       END DO
(   52)       DO jk = ikenm2, ikst, -1
(   53)          DO ji = 2, jpim1
(   54)             zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(   55)          END DO
(   56)       END DO
(   57)       
(   58) # 363 "dynzdf_iso.F90"
(  363) 
(  364) # 364
(  364)          ! 2.3 Normalization to obtain the general momentum trend va
(  365)          DO jk = 1, jpkm1
(  366)             DO ji = 2, jpim1
(  367)                va(ji,jj,jk) = ( zwx(ji,jk) - vb(ji,jj,jk) ) / z2dt
(  368)             END DO    
(  369)          END DO   
(  370) 
(  371)          ! 2.4 diagnose surface and bottom momentum fluxes
(  372)          DO ji = 2, jpim1
(  373)             ! save the surface forcing momentum fluxes
(  374)             ztsy(ji,jj) = tauy(ji,jj) / ( e3t(1)*rau0 )
(  375)             ! save bottom friction momentum fluxes






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 18

(  376)             ikbv   = min( mbathy(ji,jj+1), mbathy(ji,jj) )
(  377)             ikbvm1 = max( ikbv-1, 1 )
(  378)             ztby(ji,jj) = - avmv(ji,jj,ikbv) * zwx(ji,ikbvm1)   &
(  379)                             / ( e3t(ikbvm1)*e3w(ikbv) )
(  380)             ! subtract surface forcing and bottom friction trend from vertical
(  381)             ! diffusive momentum trend
(  382)             ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
(  383)             ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
(  384)          END DO
(  385)          !                                             ! ===============
(  386)       END DO                                           !   End of slab
(  387)       !                                                ! ===============
(  388) 
(  389)       ! save the vertical diffusive trends for diagnostic
(  390)       ! momentum trends
(  391)       IF( l_trddyn )  THEN 
(  392)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  393)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  394) 
(  395)          CALL trd_mod(uldftrd, vldftrd, jpdtdldf, 'DYN', kt)
(  396)          CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
(  397)          ztdua(:,:,:) = 0.e0
(  398)          ztdva(:,:,:) = 0.e0
(  399)          ztdua(:,:,1) = ztsx(:,:)
(  400)          ztdva(:,:,1) = ztsy(:,:)
(  401)          CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
(  402)          ztdua(:,:,:) = 0.e0
(  403)          ztdva(:,:,:) = 0.e0
(  404)          ztdua(:,:,1) = ztbx(:,:)
(  405)          ztdva(:,:,1) = ztby(:,:)
(  406)          CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
(  407)       ENDIF
(  408) 
(  409)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  410)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
(  411)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  412)       ENDIF
(  413) 
(  414)    END SUBROUTINE dyn_zdf_iso
(  415) 
(  416) 
(  417) # 427
(  427)    !!==============================================================================
(  428) END MODULE dynzdf_iso
