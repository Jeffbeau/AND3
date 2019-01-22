


PGF90 (Version     10.5)          09/23/2010  09:49:45      page 1

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

Filename: dynzdf_imp_atsk.F90

(    1) MODULE dynzdf_imp_atsk
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  dynzdf_imp_atsk  ***
(    4)    !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_zdf_imp_tsk : update the momentum trend with the vertical
(    9)    !!                     diffusion using an implicit time-stepping and
(   10)    !!                     j-k-i loops.
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE phycst          ! physical constants
(   16)    USE zdf_oce         ! ocean vertical physics
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE taumod          ! surface ocean stress
(   19)    USE trdmod          ! ocean dynamics trends 
(   20)    USE trdmod_oce      ! ocean variables trends
(   21)    USE prtctl          ! Print control
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC dyn_zdf_imp_tsk     ! called by step.F90
(   28) 
(   29)    !! * Substitutions
(   30) #  include "domzgr_substitute.h90"
(   31) #  include "vectopt_loop_substitute.h90"
(   32)    !!----------------------------------------------------------------------
(   33)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   34)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_imp_atsk.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   35)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   36)    !!----------------------------------------------------------------------
(   37) 
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE dyn_zdf_imp_tsk( kt )
(   41)       !!----------------------------------------------------------------------
(   42)       !!                  ***  ROUTINE dyn_zdf_imp_tsk  ***
(   43)       !!                   
(   44)       !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
(   45)       !!      and the surface forcing, and add it to the general trend of 
(   46)       !!      the momentum equations.






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 2

(   47)       !!
(   48)       !! ** Method  :   The vertical momentum mixing trend is given by :
(   49)       !!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
(   50)       !!      backward time stepping
(   51)       !!      Surface boundary conditions: wind stress input
(   52)       !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
(   53)       !!      Add this trend to the general trend ua :
(   54)       !!         ua = ua + dz( avmu dz(u) )
(   55)       !!
(   56)       !! ** Action : - Update (ua,va) arrays with the after vertical diffusive
(   57)       !!               mixing trend.
(   58)       !!             - Save the trends in (ztdua,ztdva) ('l_trddyn')
(   59)       !!
(   60)       !! History :
(   61)       !!   8.5  !  02-08  (G. Madec)  auto-tasking option
(   62)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   63)       !!---------------------------------------------------------------------
(   64)       !! * Modules used     
(   65)       USE oce, ONLY :    ztdua => ta,  & ! use ta as 3D workspace   
(   66)                          ztdva => sa     ! use sa as 3D workspace   
(   67) 
(   68)       !! * Arguments
(   69)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   70) 
(   71)       !! * Local declarations
(   72)       INTEGER ::   ji, jj, jk            ! dummy loop indices
(   73)       INTEGER ::   &
(   74)          ikst, ikenm2, ikstp1,         & ! temporary integers
(   75)          ikbu, ikbum1, ikbv, ikbvm1      !    "         "      
(   76)       REAL(wp) ::   &
(   77)          zrau0r, z2dt,                 & !temporary scalars
(   78)          z2dtf, zcoef, zzws
(   79)       REAL(wp), DIMENSION(jpi,jpk) ::  &
(   80)          zwx, zwy, zwz,                & ! workspace
(   81)          zwd, zws, zwi, zwt
(   82)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   83)          ztsx, ztsy, ztbx, ztby          ! temporary workspace arrays
(   84)       !!----------------------------------------------------------------------
(   85) 
(   86)       IF( kt == nit000 ) THEN
(   87)          IF(lwp) WRITE(numout,*)
(   88)          IF(lwp) WRITE(numout,*) 'dyn_zdf_imp_tsk : vertical momentum diffusion implicit operator'
(   89)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   auto-task case (j-k-i loop)'
(   90)       ENDIF
(   91) 
(   92) 
(   93)       ! 0. Local constant initialization
(   94)       ! --------------------------------
(   95)       zrau0r = 1. / rau0      ! inverse of the reference density
(   96)       z2dt   = 2. * rdt       ! Leap-frog environnement
(   97)       ztsx(:,:)   = 0.e0
(   98)       ztsy(:,:)   = 0.e0 
(   99)       ztbx(:,:)   = 0.e0
(  100)       ztby(:,:)   = 0.e0
(  101)       ! Euler time stepping when starting from rest
(  102)       IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt
(  103) 
(  104)       ! Save ua and va trends






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 3

(  105)       IF( l_trddyn )   THEN
(  106)          ztdua(:,:,:) = ua(:,:,:) 
(  107)          ztdva(:,:,:) = va(:,:,:) 
(  108)       ENDIF
(  109) 
(  110)       !                                                ! ===============
(  111)       DO jj = 2, jpjm1                                 !  Vertical slab
(  112)          !                                             ! ===============
(  113)          ! 1. Vertical diffusion on u
(  114)          ! ---------------------------
(  115) 
(  116)          ! Matrix and second member construction
(  117)          ! bottom boundary condition: only zws must be masked as avmu can take
(  118)          ! non zero value at the ocean bottom depending on the bottom friction
(  119)          ! used (see zdfmix.F)
(  120)          DO jk = 1, jpkm1
(  121)             DO ji = 2, jpim1
(  122)                zcoef = - z2dt / fse3u(ji,jj,jk)
(  123)                zwi(ji,jk) = zcoef * avmu(ji,jj,jk  ) / fse3uw(ji,jj,jk  )
(  124)                zzws       = zcoef * avmu(ji,jj,jk+1) / fse3uw(ji,jj,jk+1)
(  125)                zws(ji,jk) = zzws * umask(ji,jj,jk+1)
(  126)                zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
(  127)                zwy(ji,jk) = ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)
(  128)             END DO
(  129)          END DO
(  130) 
(  131)          ! Surface boudary conditions
(  132)          DO ji = 2, jpim1
(  133)             z2dtf = z2dt / ( fse3u(ji,jj,1)*rau0 )
(  134)             zwi(ji,1) = 0.
(  135)             zwd(ji,1) = 1. - zws(ji,1)
(  136)             zwy(ji,1) = zwy(ji,1) + z2dtf * taux(ji,jj)
(  137)          END DO
(  138) 
(  139)          ! Matrix inversion starting from the first level
(  140)          ikst = 1
(  141) !!----------------------------------------------------------------------
(  142) !!         ZDF.MATRIXSOLVER
(  143) !!       ********************
(  144) !!----------------------------------------------------------------------
(  145) !! Matrix inversion
(  146) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  147) !
(  148) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  149) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  150) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  151) !        (        ...               )( ...  ) ( ...  )
(  152) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  153) !
(  154) !   m is decomposed in the product of an upper and lower triangular
(  155) !   matrix
(  156) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  157) !   The second member is in 2d array zwy
(  158) !   The solution is in 2d array zwx
(  159) !   The 2d arry zwt and zwz are work space arrays
(  160) !
(  161) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(  162) !   the resolution of tke matrix where surface tke value is prescribed






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 4

(  163) !   so that ikstrt=2.
(  164) !!----------------------------------------------------------------------
(  165) 
(  166)          ikstp1 = ikst + 1
(  167)          ikenm2 = jpk - 2
(  168)          DO ji = 2, jpim1
(  169)             zwt(ji,ikst) = zwd(ji,ikst)
(  170)          END DO
(  171)          DO jk = ikstp1, jpkm1
(  172)             DO ji = 2, jpim1
(  173)                zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(  174)             END DO
(  175)          END DO
(  176)          DO ji = 2, jpim1
(  177)             zwz(ji,ikst) = zwy(ji,ikst)
(  178)          END DO
(  179)          DO jk = ikstp1, jpkm1
(  180)             DO ji = 2, jpim1
(  181)                zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(  182)             END DO
(  183)          END DO
(  184)          DO ji = 2, jpim1
(  185)             zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(  186)          END DO
(  187)          DO jk = ikenm2, ikst, -1
(  188)             DO ji = 2, jpim1
(  189)             zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(  190)             END DO
(  191)          END DO
(  192)       
(  193)          ! Normalization to obtain the general momentum trend ua
(  194)          DO jk = 1, jpkm1
(  195)             DO ji = 2, jpim1
(  196)                ua(ji,jj,jk) = ( zwx(ji,jk) - ub(ji,jj,jk) ) / z2dt
(  197)             END DO
(  198)          END DO
(  199) 
(  200)          ! diagnose surface and bottom momentum fluxes
(  201)          ! for trends diagnostics
(  202)          DO ji = 2, jpim1
(  203)             ! save the surface forcing momentum fluxes
(  204)             ztsx(ji,jj) = taux(ji,jj) / ( fse3u(ji,jj,1)*rau0 )
(  205)             ! save bottom friction momentum fluxes
(  206)             ikbu   = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
(  207)             ikbum1 = MAX( ikbu-1, 1 )
(  208)             ztbx(ji,jj) = - avmu(ji,jj,ikbu) * zwx(ji,ikbum1)   &
(  209)                / ( fse3u(ji,jj,ikbum1)*fse3uw(ji,jj,ikbu) )
(  210)             ! subtract surface forcing and bottom friction trend from vertical
(  211)             ! diffusive momentum trend
(  212)             ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
(  213)             ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
(  214)          END DO
(  215) 
(  216)          ! 2. Vertical diffusion on v
(  217)          ! ---------------------------
(  218) 
(  219)          ! Matrix and second member construction
(  220)          ! bottom boundary condition: only zws must be masked as avmv can take






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 5

(  221)          ! non zero value at the ocean bottom depending on the bottom friction
(  222)          ! used (see zdfmix.F)
(  223)          DO jk = 1, jpkm1
(  224)             DO ji = 2, jpim1
(  225)                zcoef = -z2dt/fse3v(ji,jj,jk)
(  226)                zwi(ji,jk) = zcoef * avmv(ji,jj,jk  ) / fse3vw(ji,jj,jk  )
(  227)                zzws       = zcoef * avmv(ji,jj,jk+1) / fse3vw(ji,jj,jk+1)
(  228)                zws(ji,jk) =  zzws * vmask(ji,jj,jk+1)
(  229)                zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
(  230)                zwy(ji,jk) = vb(ji,jj,jk) + z2dt * va(ji,jj,jk)
(  231)             END DO
(  232)          END DO
(  233) 
(  234)          ! Surface boudary conditions
(  235)          DO ji = 2, jpim1
(  236)             z2dtf = z2dt / ( fse3v(ji,jj,1)*rau0 )
(  237)             zwi(ji,1) = 0.e0
(  238)             zwd(ji,1) = 1. - zws(ji,1)
(  239)             zwy(ji,1) = zwy(ji,1) + z2dtf * tauy(ji,jj)
(  240)          END DO
(  241) 
(  242)          ! Matrix inversion starting from the first level
(  243)          ikst = 1
(  244) !!----------------------------------------------------------------------
(  245) !!         ZDF.MATRIXSOLVER
(  246) !!       ********************
(  247) !!----------------------------------------------------------------------
(  248) !! Matrix inversion
(  249) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  250) !
(  251) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  252) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  253) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  254) !        (        ...               )( ...  ) ( ...  )
(  255) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  256) !
(  257) !   m is decomposed in the product of an upper and lower triangular
(  258) !   matrix
(  259) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  260) !   The second member is in 2d array zwy
(  261) !   The solution is in 2d array zwx
(  262) !   The 2d arry zwt and zwz are work space arrays
(  263) !
(  264) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(  265) !   the resolution of tke matrix where surface tke value is prescribed
(  266) !   so that ikstrt=2.
(  267) !!----------------------------------------------------------------------
(  268) 
(  269)          ikstp1 = ikst + 1
(  270)          ikenm2 = jpk - 2
(  271)          DO ji = 2, jpim1
(  272)             zwt(ji,ikst) = zwd(ji,ikst)
(  273)          END DO
(  274)          DO jk = ikstp1, jpkm1
(  275)             DO ji = 2, jpim1
(  276)                zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(  277)             END DO
(  278)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 6

(  279)          DO ji = 2, jpim1
(  280)             zwz(ji,ikst) = zwy(ji,ikst)
(  281)          END DO
(  282)          DO jk = ikstp1, jpkm1
(  283)             DO ji = 2, jpim1
(  284)                zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(  285)             END DO
(  286)          END DO
(  287)          DO ji = 2, jpim1
(  288)             zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(  289)          END DO
(  290)          DO jk = ikenm2, ikst, -1
(  291)             DO ji = 2, jpim1
(  292)                zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(  293)             END DO
(  294)          END DO
(  295)       
(  296)          ! Normalization to obtain the general momentum trend va
(  297)          DO jk = 1, jpkm1
(  298)             DO ji = 2, jpim1
(  299)                va(ji,jj,jk) = ( zwx(ji,jk) - vb(ji,jj,jk) ) / z2dt
(  300)             END DO
(  301)          END DO
(  302) 
(  303)          ! diagnose surface and bottom momentum fluxes
(  304)          ! for trends diagnostics
(  305)          DO ji = 2, jpim1
(  306)             ! save the surface forcing momentum fluxes
(  307)          ztsy(ji,jj) = tauy(ji,jj) / ( fse3v(ji,jj,1)*rau0 )
(  308)             ! save bottom friction momentum fluxes
(  309)             ikbv   = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
(  310)             ikbvm1 = MAX( ikbv-1, 1 )
(  311)          ztby(ji,jj) = - avmv(ji,jj,ikbv) * zwx(ji,ikbvm1)   &
(  312)                / ( fse3v(ji,jj,ikbvm1)*fse3vw(ji,jj,ikbv) )
(  313)             ! subtract surface forcing and bottom friction trend from vertical
(  314)             ! diffusive momentum trend
(  315)             ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
(  316)             ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
(  317)          END DO
(  318)          !                                             ! ===============
(  319)       END DO                                           !   End of slab
(  320)       !                                                ! ===============
(  321) 
(  322)       ! save the vertical diffusion trends for diagnostic
(  323)       ! momentum trends
(  324)       IF( l_trddyn )  THEN 
(  325)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  326)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  327) 
(  328)          CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
(  329)          ztdua(:,:,:) = 0.e0
(  330)          ztdva(:,:,:) = 0.e0
(  331)          ztdua(:,:,1) = ztsx(:,:)
(  332)          ztdva(:,:,1) = ztsy(:,:)
(  333)          CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
(  334)          ztdua(:,:,:) = 0.e0
(  335)          ztdva(:,:,:) = 0.e0
(  336)          ztdua(:,:,1) = ztbx(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 7

(  337)          ztdva(:,:,1) = ztby(:,:)
(  338)          CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
(  339)       ENDIF
(  340) 
(  341)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  342)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
(  343)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  344)       ENDIF
(  345) 
(  346)    END SUBROUTINE dyn_zdf_imp_tsk
(  347) 
(  348)    !!==============================================================================
(  349) END MODULE dynzdf_imp_atsk



















































PGF90 (Version     10.5)          09/23/2010  09:49:45      page 8

(    1) # 1 "dynzdf_imp_atsk.F90"
(    1) MODULE dynzdf_imp_atsk
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  dynzdf_imp_atsk  ***
(    4)    !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_zdf_imp_tsk : update the momentum trend with the vertical
(    9)    !!                     diffusion using an implicit time-stepping and
(   10)    !!                     j-k-i loops.
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE phycst          ! physical constants
(   16)    USE zdf_oce         ! ocean vertical physics
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE taumod          ! surface ocean stress
(   19)    USE trdmod          ! ocean dynamics trends 
(   20)    USE trdmod_oce      ! ocean variables trends
(   21)    USE prtctl          ! Print control
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC dyn_zdf_imp_tsk     ! called by step.F90
(   28) 
(   29)    !! * Substitutions
(   30) # 1 "./domzgr_substitute.h90"
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
(   48) # 31 "dynzdf_imp_atsk.F90"
(   31) # 1 "./vectopt_loop_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 9

(   12) # 32 "dynzdf_imp_atsk.F90"
(   32) # 32
(   32)    !!----------------------------------------------------------------------
(   33)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   34)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_imp_atsk.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   35)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   36)    !!----------------------------------------------------------------------
(   37) 
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE dyn_zdf_imp_tsk( kt )
(   41)       !!----------------------------------------------------------------------
(   42)       !!                  ***  ROUTINE dyn_zdf_imp_tsk  ***
(   43)       !!                   
(   44)       !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
(   45)       !!      and the surface forcing, and add it to the general trend of 
(   46)       !!      the momentum equations.
(   47)       !!
(   48)       !! ** Method  :   The vertical momentum mixing trend is given by :
(   49)       !!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
(   50)       !!      backward time stepping
(   51)       !!      Surface boundary conditions: wind stress input
(   52)       !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
(   53)       !!      Add this trend to the general trend ua :
(   54)       !!         ua = ua + dz( avmu dz(u) )
(   55)       !!
(   56)       !! ** Action : - Update (ua,va) arrays with the after vertical diffusive
(   57)       !!               mixing trend.
(   58)       !!             - Save the trends in (ztdua,ztdva) ('l_trddyn')
(   59)       !!
(   60)       !! History :
(   61)       !!   8.5  !  02-08  (G. Madec)  auto-tasking option
(   62)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   63)       !!---------------------------------------------------------------------
(   64)       !! * Modules used     
(   65)       USE oce, ONLY :    ztdua => ta,  & ! use ta as 3D workspace   
(   66)                          ztdva => sa     ! use sa as 3D workspace   
(   67) 
(   68)       !! * Arguments
(   69)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   70) 
(   71)       !! * Local declarations
(   72)       INTEGER ::   ji, jj, jk            ! dummy loop indices
(   73)       INTEGER ::   &
(   74)          ikst, ikenm2, ikstp1,         & ! temporary integers
(   75)          ikbu, ikbum1, ikbv, ikbvm1      !    "         "      
(   76)       REAL(wp) ::   &
(   77)          zrau0r, z2dt,                 & !temporary scalars
(   78)          z2dtf, zcoef, zzws
(   79)       REAL(wp), DIMENSION(jpi,jpk) ::  &
(   80)          zwx, zwy, zwz,                & ! workspace
(   81)          zwd, zws, zwi, zwt
(   82)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   83)          ztsx, ztsy, ztbx, ztby          ! temporary workspace arrays
(   84)       !!----------------------------------------------------------------------
(   85) 
(   86)       IF( kt == nit000 ) THEN
(   87)          IF(lwp) WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 10

(   88)          IF(lwp) WRITE(numout,*) 'dyn_zdf_imp_tsk : vertical momentum diffusion implicit operator'
(   89)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   auto-task case (j-k-i loop)'
(   90)       ENDIF
(   91) 
(   92) 
(   93)       ! 0. Local constant initialization
(   94)       ! --------------------------------
(   95)       zrau0r = 1. / rau0      ! inverse of the reference density
(   96)       z2dt   = 2. * rdt       ! Leap-frog environnement
(   97)       ztsx(:,:)   = 0.e0
(   98)       ztsy(:,:)   = 0.e0 
(   99)       ztbx(:,:)   = 0.e0
(  100)       ztby(:,:)   = 0.e0
(  101)       ! Euler time stepping when starting from rest
(  102)       IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt
(  103) 
(  104)       ! Save ua and va trends
(  105)       IF( l_trddyn )   THEN
(  106)          ztdua(:,:,:) = ua(:,:,:) 
(  107)          ztdva(:,:,:) = va(:,:,:) 
(  108)       ENDIF
(  109) 
(  110)       !                                                ! ===============
(  111)       DO jj = 2, jpjm1                                 !  Vertical slab
(  112)          !                                             ! ===============
(  113)          ! 1. Vertical diffusion on u
(  114)          ! ---------------------------
(  115) 
(  116)          ! Matrix and second member construction
(  117)          ! bottom boundary condition: only zws must be masked as avmu can take
(  118)          ! non zero value at the ocean bottom depending on the bottom friction
(  119)          ! used (see zdfmix.F)
(  120)          DO jk = 1, jpkm1
(  121)             DO ji = 2, jpim1
(  122)                zcoef = - z2dt / e3t(jk)
(  123)                zwi(ji,jk) = zcoef * avmu(ji,jj,jk  ) / e3w(jk  )
(  124)                zzws       = zcoef * avmu(ji,jj,jk+1) / e3w(jk+1)
(  125)                zws(ji,jk) = zzws * umask(ji,jj,jk+1)
(  126)                zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
(  127)                zwy(ji,jk) = ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)
(  128)             END DO
(  129)          END DO
(  130) 
(  131)          ! Surface boudary conditions
(  132)          DO ji = 2, jpim1
(  133)             z2dtf = z2dt / ( e3t(1)*rau0 )
(  134)             zwi(ji,1) = 0.
(  135)             zwd(ji,1) = 1. - zws(ji,1)
(  136)             zwy(ji,1) = zwy(ji,1) + z2dtf * taux(ji,jj)
(  137)          END DO
(  138) 
(  139)          ! Matrix inversion starting from the first level
(  140)          ikst = 1
(  141) !!----------------------------------------------------------------------
(  142) !!         ZDF.MATRIXSOLVER
(  143) !!       ********************
(  144) !!----------------------------------------------------------------------
(  145) !! Matrix inversion






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 11

(  146) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  147) !
(  148) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  149) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  150) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  151) !        (        ...               )( ...  ) ( ...  )
(  152) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  153) !
(  154) !   m is decomposed in the product of an upper and lower triangular
(  155) !   matrix
(  156) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  157) !   The second member is in 2d array zwy
(  158) !   The solution is in 2d array zwx
(  159) !   The 2d arry zwt and zwz are work space arrays
(  160) !
(  161) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(  162) !   the resolution of tke matrix where surface tke value is prescribed
(  163) !   so that ikstrt=2.
(  164) !!----------------------------------------------------------------------
(  165) 
(  166)          ikstp1 = ikst + 1
(  167)          ikenm2 = jpk - 2
(  168)          DO ji = 2, jpim1
(  169)             zwt(ji,ikst) = zwd(ji,ikst)
(  170)          END DO
(  171)          DO jk = ikstp1, jpkm1
(  172)             DO ji = 2, jpim1
(  173)                zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(  174)             END DO
(  175)          END DO
(  176)          DO ji = 2, jpim1
(  177)             zwz(ji,ikst) = zwy(ji,ikst)
(  178)          END DO
(  179)          DO jk = ikstp1, jpkm1
(  180)             DO ji = 2, jpim1
(  181)                zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(  182)             END DO
(  183)          END DO
(  184)          DO ji = 2, jpim1
(  185)             zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(  186)          END DO
(  187)          DO jk = ikenm2, ikst, -1
(  188)             DO ji = 2, jpim1
(  189)             zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(  190)             END DO
(  191)          END DO
(  192)       
(  193)          ! Normalization to obtain the general momentum trend ua
(  194)          DO jk = 1, jpkm1
(  195)             DO ji = 2, jpim1
(  196)                ua(ji,jj,jk) = ( zwx(ji,jk) - ub(ji,jj,jk) ) / z2dt
(  197)             END DO
(  198)          END DO
(  199) 
(  200)          ! diagnose surface and bottom momentum fluxes
(  201)          ! for trends diagnostics
(  202)          DO ji = 2, jpim1
(  203)             ! save the surface forcing momentum fluxes






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 12

(  204)             ztsx(ji,jj) = taux(ji,jj) / ( e3t(1)*rau0 )
(  205)             ! save bottom friction momentum fluxes
(  206)             ikbu   = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
(  207)             ikbum1 = MAX( ikbu-1, 1 )
(  208)             ztbx(ji,jj) = - avmu(ji,jj,ikbu) * zwx(ji,ikbum1)   &
(  209)                / ( e3t(ikbum1)*e3w(ikbu) )
(  210)             ! subtract surface forcing and bottom friction trend from vertical
(  211)             ! diffusive momentum trend
(  212)             ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
(  213)             ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
(  214)          END DO
(  215) 
(  216)          ! 2. Vertical diffusion on v
(  217)          ! ---------------------------
(  218) 
(  219)          ! Matrix and second member construction
(  220)          ! bottom boundary condition: only zws must be masked as avmv can take
(  221)          ! non zero value at the ocean bottom depending on the bottom friction
(  222)          ! used (see zdfmix.F)
(  223)          DO jk = 1, jpkm1
(  224)             DO ji = 2, jpim1
(  225)                zcoef = -z2dt/e3t(jk)
(  226)                zwi(ji,jk) = zcoef * avmv(ji,jj,jk  ) / e3w(jk  )
(  227)                zzws       = zcoef * avmv(ji,jj,jk+1) / e3w(jk+1)
(  228)                zws(ji,jk) =  zzws * vmask(ji,jj,jk+1)
(  229)                zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
(  230)                zwy(ji,jk) = vb(ji,jj,jk) + z2dt * va(ji,jj,jk)
(  231)             END DO
(  232)          END DO
(  233) 
(  234)          ! Surface boudary conditions
(  235)          DO ji = 2, jpim1
(  236)             z2dtf = z2dt / ( e3t(1)*rau0 )
(  237)             zwi(ji,1) = 0.e0
(  238)             zwd(ji,1) = 1. - zws(ji,1)
(  239)             zwy(ji,1) = zwy(ji,1) + z2dtf * tauy(ji,jj)
(  240)          END DO
(  241) 
(  242)          ! Matrix inversion starting from the first level
(  243)          ikst = 1
(  244) !!----------------------------------------------------------------------
(  245) !!         ZDF.MATRIXSOLVER
(  246) !!       ********************
(  247) !!----------------------------------------------------------------------
(  248) !! Matrix inversion
(  249) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  250) !
(  251) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  252) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  253) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  254) !        (        ...               )( ...  ) ( ...  )
(  255) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  256) !
(  257) !   m is decomposed in the product of an upper and lower triangular
(  258) !   matrix
(  259) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  260) !   The second member is in 2d array zwy
(  261) !   The solution is in 2d array zwx






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 13

(  262) !   The 2d arry zwt and zwz are work space arrays
(  263) !
(  264) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(  265) !   the resolution of tke matrix where surface tke value is prescribed
(  266) !   so that ikstrt=2.
(  267) !!----------------------------------------------------------------------
(  268) 
(  269)          ikstp1 = ikst + 1
(  270)          ikenm2 = jpk - 2
(  271)          DO ji = 2, jpim1
(  272)             zwt(ji,ikst) = zwd(ji,ikst)
(  273)          END DO
(  274)          DO jk = ikstp1, jpkm1
(  275)             DO ji = 2, jpim1
(  276)                zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(  277)             END DO
(  278)          END DO
(  279)          DO ji = 2, jpim1
(  280)             zwz(ji,ikst) = zwy(ji,ikst)
(  281)          END DO
(  282)          DO jk = ikstp1, jpkm1
(  283)             DO ji = 2, jpim1
(  284)                zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(  285)             END DO
(  286)          END DO
(  287)          DO ji = 2, jpim1
(  288)             zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(  289)          END DO
(  290)          DO jk = ikenm2, ikst, -1
(  291)             DO ji = 2, jpim1
(  292)                zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(  293)             END DO
(  294)          END DO
(  295)       
(  296)          ! Normalization to obtain the general momentum trend va
(  297)          DO jk = 1, jpkm1
(  298)             DO ji = 2, jpim1
(  299)                va(ji,jj,jk) = ( zwx(ji,jk) - vb(ji,jj,jk) ) / z2dt
(  300)             END DO
(  301)          END DO
(  302) 
(  303)          ! diagnose surface and bottom momentum fluxes
(  304)          ! for trends diagnostics
(  305)          DO ji = 2, jpim1
(  306)             ! save the surface forcing momentum fluxes
(  307)          ztsy(ji,jj) = tauy(ji,jj) / ( e3t(1)*rau0 )
(  308)             ! save bottom friction momentum fluxes
(  309)             ikbv   = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
(  310)             ikbvm1 = MAX( ikbv-1, 1 )
(  311)          ztby(ji,jj) = - avmv(ji,jj,ikbv) * zwx(ji,ikbvm1)   &
(  312)                / ( e3t(ikbvm1)*e3w(ikbv) )
(  313)             ! subtract surface forcing and bottom friction trend from vertical
(  314)             ! diffusive momentum trend
(  315)             ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
(  316)             ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
(  317)          END DO
(  318)          !                                             ! ===============
(  319)       END DO                                           !   End of slab






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 14

(  320)       !                                                ! ===============
(  321) 
(  322)       ! save the vertical diffusion trends for diagnostic
(  323)       ! momentum trends
(  324)       IF( l_trddyn )  THEN 
(  325)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  326)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  327) 
(  328)          CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
(  329)          ztdua(:,:,:) = 0.e0
(  330)          ztdva(:,:,:) = 0.e0
(  331)          ztdua(:,:,1) = ztsx(:,:)
(  332)          ztdva(:,:,1) = ztsy(:,:)
(  333)          CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
(  334)          ztdua(:,:,:) = 0.e0
(  335)          ztdva(:,:,:) = 0.e0
(  336)          ztdua(:,:,1) = ztbx(:,:)
(  337)          ztdva(:,:,1) = ztby(:,:)
(  338)          CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
(  339)       ENDIF
(  340) 
(  341)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  342)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
(  343)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  344)       ENDIF
(  345) 
(  346)    END SUBROUTINE dyn_zdf_imp_tsk
(  347) 
(  348)    !!==============================================================================
(  349) END MODULE dynzdf_imp_atsk
