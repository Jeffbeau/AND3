


PGF90 (Version     10.5)          09/23/2010  09:50:11      page 1

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

Filename: traadv_muscl2.F90

(    1) MODULE traadv_muscl2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_muscl2  ***
(    4)    !! Ocean active tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_muscl2 : update the tracer trend with the horizontal
(    9)    !!                    and vertical advection trends using MUSCL2 scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE trdmod          ! ocean active tracers trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   18)    USE trabbl          ! tracers: bottom boundary layer
(   19)    USE lib_mpp
(   20)    USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
(   21)    USE diaptr          ! poleward transport diagnostics
(   22)    USE prtctl          ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Accessibility
(   28)    PUBLIC tra_adv_muscl2        ! routine called by step.F90
(   29) 
(   30)    !! * Substitutions
(   31) #  include "domzgr_substitute.h90"
(   32) #  include "vectopt_loop_substitute.h90"
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_muscl2.F90,v 1.10 2005/12/28 09:25:10 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE tra_adv_muscl2( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                   ***  ROUTINE tra_adv_muscl2  ***
(   44)       !!
(   45)       !! ** Purpose :   Compute the now trend due to total advection of T and 
(   46)       !!      S using a MUSCL scheme (Monotone Upstream-centered Scheme for






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 2

(   47)       !!      Conservation Laws) and add it to the general tracer trend.
(   48)       !!
(   49)       !! ** Method  : MUSCL scheme plus centered scheme at ocean boundaries
(   50)       !!
(   51)       !! ** Action  : - update (ta,sa) with the now advective tracer trends
(   52)       !!              - save trends in (ztdta,ztdsa) ('key_trdtra')
(   53)       !!
(   54)       !! References :                
(   55)       !!      Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
(   56)       !!	IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
(   57)       !!
(   58)       !! History :
(   59)       !!        !  06-00  (A.Estublier)  for passive tracers
(   60)       !!        !  01-08  (E.Durand G.Madec)  adapted for T & S
(   61)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   62)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   63)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   64)       !!----------------------------------------------------------------------
(   65)       !! * modules used
(   66) #if defined key_trabbl_adv
(   67)       USE oce                , zun => ua,  &  ! use ua as workspace
(   68)          &                     zvn => va      ! use va as workspace
(   69)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
(   70) #else
(   71)       USE oce                , zun => un,  &  ! When no bbl, zun == un
(   72)                                zvn => vn,  &  !              zvn == vn
(   73)                                zwn => wn      !              zwn == wn
(   74) #endif
(   75) 
(   76)       !! * Arguments
(   77)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   78) 
(   79)       !! * Local declarations
(   80)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(   81)       REAL(wp) ::   &
(   82)          zu, zv, zw, zeu, zev,           &  
(   83)          zew, zbtr, zstep,               &
(   84)          z0u, z0v, z0w,                  &
(   85)          zzt1, zzt2, zalpha,             &
(   86)          zzs1, zzs2, z2,                 &
(   87)          zta, zsa
(   88)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
(   89)          zt1, zt2, ztp1, ztp2,   &
(   90)          zs1, zs2, zsp1, zsp2,   &
(   91)          ztdta, ztdsa
(   92)       !!----------------------------------------------------------------------
(   93) 
(   94)       IF( kt == nit000 .AND. lwp ) THEN
(   95)          WRITE(numout,*)
(   96)          WRITE(numout,*) 'tra_adv_muscl2 : MUSCL2 advection scheme'
(   97)          WRITE(numout,*) '~~~~~~~~~~~~~~~'
(   98)       ENDIF
(   99) 
(  100)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  101)           z2=1.
(  102)       ELSE
(  103)           z2=2.
(  104)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 3

(  105) 
(  106)       ! Save ta and sa trends
(  107)       IF( l_trdtra )   THEN
(  108)          ztdta(:,:,:) = ta(:,:,:) 
(  109)          ztdsa(:,:,:) = sa(:,:,:) 
(  110)          l_adv = 'mu2'
(  111)       ENDIF
(  112) 
(  113) #if defined key_trabbl_adv
(  114)       ! Advective bottom boundary layer
(  115)       ! -------------------------------
(  116)       zun(:,:,:) = un (:,:,:) - u_bbl(:,:,:)
(  117)       zvn(:,:,:) = vn (:,:,:) - v_bbl(:,:,:)
(  118)       zwn(:,:,:) = wn (:,:,:) + w_bbl( :,:,:)
(  119) #endif
(  120) 
(  121) 
(  122)       ! I. Horizontal advective fluxes
(  123)       ! ------------------------------
(  124) 
(  125)       ! first guess of the slopes
(  126)       ! interior values
(  127)       DO jk = 1, jpkm1
(  128)          DO jj = 1, jpjm1      
(  129)             DO ji = 1, fs_jpim1   ! vector opt.
(  130)                zt1(ji,jj,jk) = umask(ji,jj,jk) * ( tb(ji+1,jj,jk) - tb(ji,jj,jk) )
(  131)                zs1(ji,jj,jk) = umask(ji,jj,jk) * ( sb(ji+1,jj,jk) - sb(ji,jj,jk) )
(  132)                zt2(ji,jj,jk) = vmask(ji,jj,jk) * ( tb(ji,jj+1,jk) - tb(ji,jj,jk) )
(  133)                zs2(ji,jj,jk) = vmask(ji,jj,jk) * ( sb(ji,jj+1,jk) - sb(ji,jj,jk) )
(  134)             END DO
(  135)          END DO
(  136)       END DO
(  137)       ! bottom values
(  138)       zt1(:,:,jpk) = 0.e0    ;    zt2(:,:,jpk) = 0.e0
(  139)       zs1(:,:,jpk) = 0.e0    ;    zs2(:,:,jpk) = 0.e0
(  140) 
(  141)       ! lateral boundary conditions on zt1, zt2 ; zs1, zs2   (changed sign)
(  142)       CALL lbc_lnk( zt1, 'U', -1. )   ;   CALL lbc_lnk( zs1, 'U', -1. )
(  143)       CALL lbc_lnk( zt2, 'V', -1. )   ;   CALL lbc_lnk( zs2, 'V', -1. )
(  144) 
(  145)       ! Slopes
(  146)       ! interior values
(  147)       DO jk = 1, jpkm1
(  148)          DO jj = 2, jpj
(  149)             DO ji = fs_2, jpi   ! vector opt.
(  150)                ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji-1,jj  ,jk) )   &
(  151)                   &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji-1,jj  ,jk) ) )
(  152)                zsp1(ji,jj,jk) =                    ( zs1(ji,jj,jk) + zs1(ji-1,jj  ,jk) )   &
(  153)                   &           * ( 0.25 + SIGN( 0.25, zs1(ji,jj,jk) * zs1(ji-1,jj  ,jk) ) )
(  154)                ztp2(ji,jj,jk) =                    ( zt2(ji,jj,jk) + zt2(ji  ,jj-1,jk) )   &
(  155)                   &           * ( 0.25 + SIGN( 0.25, zt2(ji,jj,jk) * zt2(ji  ,jj-1,jk) ) )
(  156)                zsp2(ji,jj,jk) =                    ( zs2(ji,jj,jk) + zs2(ji  ,jj-1,jk) )   &
(  157)                   &           * ( 0.25 + SIGN( 0.25, zs2(ji,jj,jk) * zs2(ji  ,jj-1,jk) ) )
(  158)             END DO
(  159)          END DO
(  160)       END DO
(  161)       ! bottom values
(  162)       ztp1(:,:,jpk) = 0.e0    ;    ztp2(:,:,jpk) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 4

(  163)       zsp1(:,:,jpk) = 0.e0    ;    zsp2(:,:,jpk) = 0.e0
(  164) 
(  165)       ! Slopes limitation
(  166)       DO jk = 1, jpkm1
(  167)          DO jj = 2, jpj
(  168)             DO ji = fs_2, jpi   ! vector opt.
(  169)                ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  170)                   &           * MIN(    ABS( ztp1(ji  ,jj,jk) ),   &
(  171)                   &                  2.*ABS( zt1 (ji-1,jj,jk) ),   &
(  172)                   &                  2.*ABS( zt1 (ji  ,jj,jk) ) )
(  173)                zsp1(ji,jj,jk) = SIGN( 1., zsp1(ji,jj,jk) )   &
(  174)                   &           * MIN(    ABS( zsp1(ji  ,jj,jk) ),   &
(  175)                   &                  2.*ABS( zs1 (ji-1,jj,jk) ),   &
(  176)                   &                  2.*ABS( zs1 (ji  ,jj,jk) ) )
(  177)                ztp2(ji,jj,jk) = SIGN( 1., ztp2(ji,jj,jk) )   &
(  178)                   &           * MIN(    ABS( ztp2(ji,jj  ,jk) ),   &
(  179)                   &                  2.*ABS( zt2 (ji,jj-1,jk) ),   &
(  180)                   &                  2.*ABS( zt2 (ji,jj  ,jk) ) )
(  181)                zsp2(ji,jj,jk) = SIGN( 1., zsp2(ji,jj,jk) )   &
(  182)                   &           * MIN(    ABS( zsp2(ji,jj  ,jk) ),   &
(  183)                   &                  2.*ABS( zs2 (ji,jj-1,jk) ),   &
(  184)                   &                  2.*ABS( zs2 (ji,jj  ,jk) ) )
(  185)             END DO
(  186)          END DO
(  187)       END DO        
(  188) 
(  189)       ! Advection terms
(  190)       ! interior values
(  191)       DO jk = 1, jpkm1
(  192)          zstep  = z2 * rdttra(jk)
(  193)          DO jj = 2, jpjm1      
(  194)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  195)                ! volume fluxes
(  196) #if defined key_s_coord || defined key_partial_steps
(  197)                zeu = e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  198)                zev = e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  199) #else
(  200)                zeu = e2u(ji,jj) * zun(ji,jj,jk)
(  201)                zev = e1v(ji,jj) * zvn(ji,jj,jk)
(  202) #endif
(  203)                ! MUSCL fluxes
(  204)                z0u = SIGN( 0.5, zun(ji,jj,jk) )            
(  205)                zalpha = 0.5 - z0u
(  206)                zu  = z0u - 0.5 * zun(ji,jj,jk) * zstep / e1u(ji,jj)
(  207)                zzt1 = tb(ji+1,jj,jk) + zu*ztp1(ji+1,jj,jk)
(  208)                zzt2 = tb(ji  ,jj,jk) + zu*ztp1(ji  ,jj,jk)
(  209)                zzs1 = sb(ji+1,jj,jk) + zu*zsp1(ji+1,jj,jk)
(  210)                zzs2 = sb(ji  ,jj,jk) + zu*zsp1(ji  ,jj,jk)
(  211)                zt1(ji,jj,jk) = zeu * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  212)                zs1(ji,jj,jk) = zeu * ( zalpha * zzs1 + (1.-zalpha) * zzs2 )
(  213) 
(  214)                z0v = SIGN( 0.5, zvn(ji,jj,jk) )            
(  215)                zalpha = 0.5 - z0v
(  216)                zv  = z0v - 0.5 * zvn(ji,jj,jk) * zstep / e2v(ji,jj)
(  217)                zzt1 = tb(ji,jj+1,jk) + zv*ztp2(ji,jj+1,jk)
(  218)                zzt2 = tb(ji,jj  ,jk) + zv*ztp2(ji,jj  ,jk)
(  219)                zzs1 = sb(ji,jj+1,jk) + zv*zsp2(ji,jj+1,jk)
(  220)                zzs2 = sb(ji,jj  ,jk) + zv*zsp2(ji,jj  ,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 5

(  221)                zt2(ji,jj,jk) = zev * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  222)                zs2(ji,jj,jk) = zev * ( zalpha * zzs1 + (1.-zalpha) * zzs2 )
(  223)             END DO
(  224)          END DO
(  225)       END DO
(  226) 
(  227)       !!!!  centered scheme at lateral b.C. if off-shore velocity
(  228)       DO jk = 1, jpkm1
(  229)         DO jj = 2, jpjm1
(  230)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  231) #if defined key_s_coord || defined key_partial_steps
(  232)                zev = e1v(ji,jj) * fse3v(ji,jj,jk)
(  233)                IF( umask(ji,jj,jk) == 0. ) THEN
(  234)                   IF( zun(ji+1,jj,jk) > 0. .AND. ji /= jpi ) THEN
(  235)                      zt1(ji+1,jj,jk) = e2u(ji+1,jj)* fse3u(ji+1,jj,jk)   &
(  236)                         &            * zun(ji+1,jj,jk) * ( tb(ji+1,jj,jk) + tb(ji+2,jj,jk) ) * 0.5
(  237)                      zs1(ji+1,jj,jk) = e2u(ji+1,jj)* fse3u(ji+1,jj,jk)   &
(  238)                         &            * zun(ji+1,jj,jk) * ( sb(ji+1,jj,jk) + sb(ji+2,jj,jk) ) * 0.5
(  239)                   ENDIF
(  240)                   IF( zun(ji-1,jj,jk) < 0. ) THEN
(  241)                      zt1(ji-1,jj,jk) = e2u(ji-1,jj)* fse3u(ji-1,jj,jk)   &
(  242)                         &            * zun(ji-1,jj,jk) * ( tb(ji-1,jj,jk) + tb(ji  ,jj,jk) ) * 0.5
(  243)                      zs1(ji-1,jj,jk) = e2u(ji-1,jj)* fse3u(ji-1,jj,jk)   &
(  244)                         &            * zun(ji-1,jj,jk) * ( sb(ji-1,jj,jk) + sb(ji  ,jj,jk) ) * 0.5
(  245)                   ENDIF
(  246)                ENDIF
(  247)                IF( vmask(ji,jj,jk) == 0. ) THEN
(  248)                   IF( zvn(ji,jj+1,jk) > 0. .AND. jj /= jpj ) THEN
(  249)                      zt2(ji,jj+1,jk) = e1v(ji,jj+1) * fse3v(ji,jj+1,jk)   &
(  250)                         &            * zvn(ji,jj+1,jk) * ( tb(ji,jj+1,jk) + tb(ji,jj+2,jk) ) * 0.5
(  251)                      zs2(ji,jj+1,jk) = e1v(ji,jj+1) * fse3v(ji,jj+1,jk)   &
(  252)                         &            * zvn(ji,jj+1,jk) * ( sb(ji,jj+1,jk) + sb(ji,jj+2,jk) ) * 0.5
(  253)                   ENDIF
(  254)                   IF( zvn(ji,jj-1,jk) < 0. ) THEN
(  255)                      zt2(ji,jj-1,jk) = e1v(ji,jj-1)* fse3v(ji,jj-1,jk)   &
(  256)                         &            * zvn(ji,jj-1,jk) * ( tb(ji,jj-1,jk) + tb(ji  ,jj,jk) ) * 0.5
(  257)                      zs2(ji,jj-1,jk) = e1v(ji,jj-1)* fse3v(ji,jj-1,jk)   &
(  258)                         &            * zvn(ji,jj-1,jk) * ( sb(ji,jj-1,jk) + sb(ji  ,jj,jk) ) * 0.5
(  259)                   ENDIF
(  260)                ENDIF
(  261) 
(  262) #else
(  263)                IF( umask(ji,jj,jk) == 0. ) THEN
(  264)                   IF( zun(ji+1,jj,jk) > 0. .AND. ji /= jpi ) THEN
(  265)                      zt1(ji+1,jj,jk) = e2u(ji+1,jj) * zun(ji+1,jj,jk) * ( tb(ji+1,jj,jk) + tb(ji+2,jj,jk) ) * 0.5
(  266)                      zs1(ji+1,jj,jk) = e2u(ji+1,jj) * zun(ji+1,jj,jk) * ( sb(ji+1,jj,jk) + sb(ji+2,jj,jk) ) * 0.5
(  267)                   ENDIF
(  268)                   IF( zun(ji-1,jj,jk) < 0. ) THEN
(  269)                      zt1(ji-1,jj,jk) = e2u(ji-1,jj) * zun(ji-1,jj,jk) * ( tb(ji-1,jj,jk) + tb(ji  ,jj,jk) ) * 0.5
(  270)                      zs1(ji-1,jj,jk) = e2u(ji-1,jj) * zun(ji-1,jj,jk) * ( sb(ji-1,jj,jk) + sb(ji  ,jj,jk) ) * 0.5
(  271)                   ENDIF
(  272)                ENDIF
(  273)                IF( vmask(ji,jj,jk) == 0. ) THEN
(  274)                   IF( zvn(ji,jj+1,jk) > 0. .AND. jj /= jpj ) THEN
(  275)                      zt2(ji,jj+1,jk) = e1v(ji,jj+1) * zvn(ji,jj+1,jk) * ( tb(ji,jj+1,jk) + tb(ji,jj+2,jk) ) * 0.5
(  276)                      zs2(ji,jj+1,jk) = e1v(ji,jj+1) * zvn(ji,jj+1,jk) * ( sb(ji,jj+1,jk) + sb(ji,jj+2,jk) ) * 0.5
(  277)                   ENDIF
(  278)                   IF( zvn(ji,jj-1,jk) < 0. ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 6

(  279)                      zt2(ji,jj-1,jk) = e1v(ji,jj-1) * zvn(ji,jj-1,jk) * ( tb(ji,jj-1,jk) + tb(ji  ,jj,jk) ) * 0.5
(  280)                      zs2(ji,jj-1,jk) = e1v(ji,jj-1) * zvn(ji,jj-1,jk) * ( sb(ji,jj-1,jk) + sb(ji  ,jj,jk) ) * 0.5
(  281)                   ENDIF
(  282)                ENDIF
(  283) #endif
(  284)             END DO
(  285)          END DO
(  286)       END DO
(  287) 
(  288)       ! lateral boundary conditions on zt1, zt2 ; zs1, zs2   (changed sign)
(  289)       CALL lbc_lnk( zt1, 'U', -1. )   ;   CALL lbc_lnk( zs1, 'U', -1. ) 
(  290)       CALL lbc_lnk( zt2, 'V', -1. )   ;   CALL lbc_lnk( zs2, 'V', -1. )
(  291) 
(  292)       ! Save MUSCL fluxes to compute i- & j- horizontal 
(  293)       ! advection trends in the MLD
(  294)       IF( l_trdtra )   THEN
(  295)          ! save i- terms
(  296)          tladi(:,:,:) = zt1(:,:,:) 
(  297)          sladi(:,:,:) = zs1(:,:,:) 
(  298)          ! save j- terms
(  299)          tladj(:,:,:) = zt2(:,:,:) 
(  300)          sladj(:,:,:) = zs2(:,:,:) 
(  301)       ENDIF
(  302) 
(  303)       ! Compute & add the horizontal advective trend
(  304) 
(  305)       DO jk = 1, jpkm1
(  306)          DO jj = 2, jpjm1      
(  307)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  308) #if defined key_s_coord || defined key_partial_steps
(  309)                zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  310) #else
(  311)                zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  312) #endif
(  313)                ! horizontal advective trends
(  314)                zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk  )   &
(  315)                   &           + zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk  ) )
(  316)                zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji-1,jj  ,jk  )   &
(  317)                   &           + zs2(ji,jj,jk) - zs2(ji  ,jj-1,jk  ) ) 
(  318)                ! add it to the general tracer trends
(  319)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  320)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  321)             END DO
(  322)         END DO
(  323)       END DO        
(  324) 
(  325)       ! Save the horizontal advective trends for diagnostic
(  326) 
(  327)       IF( l_trdtra )   THEN
(  328)          ! Recompute the hoizontal advection zta & zsa trends computed 
(  329)          ! at the step 2. above in making the difference between the new 
(  330)          ! trends and the previous one ta()/sa - ztdta()/ztdsa() and add
(  331)          ! the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
(  332)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:) 
(  333)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)
(  334) 
(  335)          CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 7

(  337)          ! Save the new ta and sa trends
(  338)          ztdta(:,:,:) = ta(:,:,:) 
(  339)          ztdsa(:,:,:) = sa(:,:,:) 
(  340) 
(  341)       ENDIF
(  342) 
(  343)       IF(ln_ctl) THEN
(  344)          CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl2 had  - Ta: ', mask1=tmask, &
(  345)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  346)       ENDIF
(  347) 
(  348)       ! "zonal" mean advective heat and salt transport
(  349)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  350) # if defined key_s_coord || defined key_partial_steps
(  351)          pht_adv(:) = ptr_vj( zt2(:,:,:) )
(  352)          pst_adv(:) = ptr_vj( zs2(:,:,:) )
(  353) # else
(  354)          DO jk = 1, jpkm1
(  355)             DO jj = 2, jpjm1
(  356)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  357)                  zt2(ji,jj,jk) = zt2(ji,jj,jk) * fse3v(ji,jj,jk)
(  358)                  zs2(ji,jj,jk) = zs2(ji,jj,jk) * fse3v(ji,jj,jk)
(  359)                END DO
(  360)             END DO
(  361)          END DO
(  362)          pht_adv(:) = ptr_vj( zt2(:,:,:) )
(  363)          pst_adv(:) = ptr_vj( zs2(:,:,:) )
(  364) # endif
(  365)       ENDIF
(  366) 
(  367)       ! II. Vertical advective fluxes
(  368)       ! -----------------------------
(  369)       
(  370)       ! First guess of the slope
(  371)       ! interior values
(  372)       DO jk = 2, jpkm1
(  373)          zt1(:,:,jk) = tmask(:,:,jk) * ( tb(:,:,jk-1) - tb(:,:,jk) )
(  374)          zs1(:,:,jk) = tmask(:,:,jk) * ( sb(:,:,jk-1) - sb(:,:,jk) )
(  375)       END DO
(  376)       ! surface & bottom boundary conditions
(  377)       zt1 (:,:, 1 ) = 0.e0    ;    zt1 (:,:,jpk) = 0.e0
(  378)       zs1 (:,:, 1 ) = 0.e0    ;    zs1 (:,:,jpk) = 0.e0
(  379) 
(  380)       ! Slopes
(  381)       DO jk = 2, jpkm1
(  382)          DO jj = 1, jpj
(  383)             DO ji = 1, jpi
(  384)                ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji,jj,jk+1) )   &
(  385)                   &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji,jj,jk+1) ) )
(  386)                zsp1(ji,jj,jk) =                    ( zs1(ji,jj,jk) + zs1(ji,jj,jk+1) )   &
(  387)                   &           * ( 0.25 + SIGN( 0.25, zs1(ji,jj,jk) * zs1(ji,jj,jk+1) ) )
(  388)             END DO
(  389)          END DO
(  390)       END DO
(  391) 
(  392)       ! Slopes limitation
(  393)       ! interior values
(  394)       DO jk = 2, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 8

(  395)          DO jj = 1, jpj
(  396)             DO ji = 1, jpi
(  397)                ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  398)                   &           * MIN(    ABS( ztp1(ji,jj,jk  ) ),   &
(  399)                   &                  2.*ABS( zt1 (ji,jj,jk+1) ),   &
(  400)                   &                  2.*ABS( zt1 (ji,jj,jk  ) ) )
(  401)                zsp1(ji,jj,jk) = SIGN( 1., zsp1(ji,jj,jk) )   &
(  402)                   &           * MIN(    ABS( zsp1(ji,jj,jk  ) ),   &
(  403)                   &                  2.*ABS( zs1 (ji,jj,jk+1) ),   &
(  404)                   &                  2.*ABS( zs1 (ji,jj,jk  ) ) )
(  405)             END DO
(  406)          END DO
(  407)       END DO
(  408)       ! surface values
(  409)       ztp1(:,:,1) = 0.e0
(  410)       zsp1(:,:,1) = 0.e0
(  411) 
(  412)       ! vertical advective flux
(  413)       ! interior values
(  414)       DO jk = 1, jpkm1
(  415)          zstep  = z2 * rdttra(jk)
(  416)          DO jj = 2, jpjm1      
(  417)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  418)                zew = zwn(ji,jj,jk+1)
(  419)                z0w = SIGN( 0.5, zwn(ji,jj,jk+1) )
(  420)                zalpha = 0.5 + z0w
(  421)                zw  = z0w - 0.5 * zwn(ji,jj,jk+1)*zstep / fse3w(ji,jj,jk+1)
(  422)                zzt1 = tb(ji,jj,jk+1) + zw*ztp1(ji,jj,jk+1)
(  423)                zzt2 = tb(ji,jj,jk  ) + zw*ztp1(ji,jj,jk  )
(  424)                zzs1 = sb(ji,jj,jk+1) + zw*zsp1(ji,jj,jk+1)
(  425)                zzs2 = sb(ji,jj,jk  ) + zw*zsp1(ji,jj,jk  )
(  426)                zt1(ji,jj,jk+1) = zew * ( zalpha * zzt1 + (1.-zalpha)*zzt2 )
(  427)                zs1(ji,jj,jk+1) = zew * ( zalpha * zzs1 + (1.-zalpha)*zzs2 )
(  428)             END DO
(  429)          END DO
(  430)       END DO
(  431)       DO jk = 2, jpkm1
(  432)         DO jj = 2, jpjm1
(  433)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  434)                IF( tmask(ji,jj,jk+1) == 0. ) THEN
(  435)                   IF( zwn(ji,jj,jk) > 0. ) THEN
(  436)                      zt1(ji,jj,jk) = zwn(ji,jj,jk) * ( tb(ji,jj,jk-1) + tb(ji,jj,jk) ) * 0.5
(  437)                      zs1(ji,jj,jk) = zwn(ji,jj,jk) * ( sb(ji,jj,jk-1) + sb(ji,jj,jk) ) * 0.5
(  438)                   ENDIF
(  439)                ENDIF
(  440)             END DO
(  441)          END DO
(  442)       END DO
(  443) 
(  444)       ! surface values
(  445)       IF( lk_dynspg_rl ) THEN                           ! rigid lid : flux set to zero
(  446)          zt1(:,:, 1 ) = 0.e0
(  447)          zs1(:,:, 1 ) = 0.e0
(  448)       ELSE                                              ! free surface
(  449)          zt1(:,:, 1 ) = zwn(:,:,1) * tb(:,:,1)
(  450)          zs1(:,:, 1 ) = zwn(:,:,1) * sb(:,:,1)
(  451)       ENDIF
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 9

(  453)       ! bottom values
(  454)       zt1(:,:,jpk) = 0.e0
(  455)       zs1(:,:,jpk) = 0.e0
(  456) 
(  457) 
(  458)       ! Compute & add the vertical advective trend
(  459) 
(  460)       DO jk = 1, jpkm1
(  461)          DO jj = 2, jpjm1      
(  462)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  463)                zbtr = 1. / fse3t(ji,jj,jk)
(  464)                ! horizontal advective trends
(  465)                zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji,jj,jk+1) )
(  466)                zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji,jj,jk+1) )
(  467)                ! add it to the general tracer trends
(  468)                ta(ji,jj,jk) =  ta(ji,jj,jk) + zta
(  469)                sa(ji,jj,jk) =  sa(ji,jj,jk) + zsa
(  470)             END DO
(  471)          END DO
(  472)       END DO
(  473) 
(  474)       ! Save the vertical advective trends for diagnostic
(  475) 
(  476)       IF( l_trdtra )   THEN
(  477)          ! Recompute the vertical advection zta & zsa trends computed 
(  478)          ! at the step 2. above in making the difference between the new 
(  479)          ! trends and the previous one: ta()/sa - ztdta()/ztdsa() and substract
(  480)          ! the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
(  481)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tn(:,:,:) * hdivn(:,:,:) 
(  482)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  483) 
(  484)          CALL trd_mod(ztdta, ztdsa, jpttdzad, 'TRA', kt)
(  485)       ENDIF
(  486) 
(  487)       IF(ln_ctl) THEN
(  488)          CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl2 zad  - Ta: ', mask1=tmask, &
(  489)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  490) 
(  491)       ENDIF
(  492) 
(  493)    END SUBROUTINE tra_adv_muscl2
(  494) 
(  495)    !!======================================================================
(  496) END MODULE traadv_muscl2




















PGF90 (Version     10.5)          09/23/2010  09:50:11      page 10

(    1) # 1 "traadv_muscl2.F90"
(    1) MODULE traadv_muscl2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_muscl2  ***
(    4)    !! Ocean active tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_muscl2 : update the tracer trend with the horizontal
(    9)    !!                    and vertical advection trends using MUSCL2 scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE trdmod          ! ocean active tracers trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   18)    USE trabbl          ! tracers: bottom boundary layer
(   19)    USE lib_mpp
(   20)    USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
(   21)    USE diaptr          ! poleward transport diagnostics
(   22)    USE prtctl          ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Accessibility
(   28)    PUBLIC tra_adv_muscl2        ! routine called by step.F90
(   29) 
(   30)    !! * Substitutions
(   31) # 1 "./domzgr_substitute.h90"
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
(   48) # 32 "traadv_muscl2.F90"
(   32) # 1 "./vectopt_loop_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 11

(   11)    !!----------------------------------------------------------------------
(   12) # 33 "traadv_muscl2.F90"
(   33) # 33
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_muscl2.F90,v 1.10 2005/12/28 09:25:10 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE tra_adv_muscl2( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                   ***  ROUTINE tra_adv_muscl2  ***
(   44)       !!
(   45)       !! ** Purpose :   Compute the now trend due to total advection of T and 
(   46)       !!      S using a MUSCL scheme (Monotone Upstream-centered Scheme for
(   47)       !!      Conservation Laws) and add it to the general tracer trend.
(   48)       !!
(   49)       !! ** Method  : MUSCL scheme plus centered scheme at ocean boundaries
(   50)       !!
(   51)       !! ** Action  : - update (ta,sa) with the now advective tracer trends
(   52)       !!              - save trends in (ztdta,ztdsa) ('key_trdtra')
(   53)       !!
(   54)       !! References :                
(   55)       !!      Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
(   56)       !!	IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
(   57)       !!
(   58)       !! History :
(   59)       !!        !  06-00  (A.Estublier)  for passive tracers
(   60)       !!        !  01-08  (E.Durand G.Madec)  adapted for T & S
(   61)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   62)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   63)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   64)       !!----------------------------------------------------------------------
(   65)       !! * modules used
(   66) # 71
(   71)       USE oce                , zun => un,  &  ! When no bbl, zun == un
(   72)                                zvn => vn,  &  !              zvn == vn
(   73)                                zwn => wn      !              zwn == wn
(   74) 
(   75) # 76
(   76)       !! * Arguments
(   77)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   78) 
(   79)       !! * Local declarations
(   80)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(   81)       REAL(wp) ::   &
(   82)          zu, zv, zw, zeu, zev,           &  
(   83)          zew, zbtr, zstep,               &
(   84)          z0u, z0v, z0w,                  &
(   85)          zzt1, zzt2, zalpha,             &
(   86)          zzs1, zzs2, z2,                 &
(   87)          zta, zsa
(   88)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
(   89)          zt1, zt2, ztp1, ztp2,   &
(   90)          zs1, zs2, zsp1, zsp2,   &
(   91)          ztdta, ztdsa






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 12

(   92)       !!----------------------------------------------------------------------
(   93) 
(   94)       IF( kt == nit000 .AND. lwp ) THEN
(   95)          WRITE(numout,*)
(   96)          WRITE(numout,*) 'tra_adv_muscl2 : MUSCL2 advection scheme'
(   97)          WRITE(numout,*) '~~~~~~~~~~~~~~~'
(   98)       ENDIF
(   99) 
(  100)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  101)           z2=1.
(  102)       ELSE
(  103)           z2=2.
(  104)       ENDIF
(  105) 
(  106)       ! Save ta and sa trends
(  107)       IF( l_trdtra )   THEN
(  108)          ztdta(:,:,:) = ta(:,:,:) 
(  109)          ztdsa(:,:,:) = sa(:,:,:) 
(  110)          l_adv = 'mu2'
(  111)       ENDIF
(  112) 
(  113) 
(  114) 
(  115) # 122
(  122)       ! I. Horizontal advective fluxes
(  123)       ! ------------------------------
(  124) 
(  125)       ! first guess of the slopes
(  126)       ! interior values
(  127)       DO jk = 1, jpkm1
(  128)          DO jj = 1, jpjm1      
(  129)             DO ji = 1, jpim1   ! vector opt.
(  130)                zt1(ji,jj,jk) = umask(ji,jj,jk) * ( tb(ji+1,jj,jk) - tb(ji,jj,jk) )
(  131)                zs1(ji,jj,jk) = umask(ji,jj,jk) * ( sb(ji+1,jj,jk) - sb(ji,jj,jk) )
(  132)                zt2(ji,jj,jk) = vmask(ji,jj,jk) * ( tb(ji,jj+1,jk) - tb(ji,jj,jk) )
(  133)                zs2(ji,jj,jk) = vmask(ji,jj,jk) * ( sb(ji,jj+1,jk) - sb(ji,jj,jk) )
(  134)             END DO
(  135)          END DO
(  136)       END DO
(  137)       ! bottom values
(  138)       zt1(:,:,jpk) = 0.e0    ;    zt2(:,:,jpk) = 0.e0
(  139)       zs1(:,:,jpk) = 0.e0    ;    zs2(:,:,jpk) = 0.e0
(  140) 
(  141)       ! lateral boundary conditions on zt1, zt2 ; zs1, zs2   (changed sign)
(  142)       CALL lbc_lnk( zt1, 'U', -1. )   ;   CALL lbc_lnk( zs1, 'U', -1. )
(  143)       CALL lbc_lnk( zt2, 'V', -1. )   ;   CALL lbc_lnk( zs2, 'V', -1. )
(  144) 
(  145)       ! Slopes
(  146)       ! interior values
(  147)       DO jk = 1, jpkm1
(  148)          DO jj = 2, jpj
(  149)             DO ji = 2, jpi   ! vector opt.
(  150)                ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji-1,jj  ,jk) )   &
(  151)                   &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji-1,jj  ,jk) ) )
(  152)                zsp1(ji,jj,jk) =                    ( zs1(ji,jj,jk) + zs1(ji-1,jj  ,jk) )   &
(  153)                   &           * ( 0.25 + SIGN( 0.25, zs1(ji,jj,jk) * zs1(ji-1,jj  ,jk) ) )
(  154)                ztp2(ji,jj,jk) =                    ( zt2(ji,jj,jk) + zt2(ji  ,jj-1,jk) )   &
(  155)                   &           * ( 0.25 + SIGN( 0.25, zt2(ji,jj,jk) * zt2(ji  ,jj-1,jk) ) )






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 13

(  156)                zsp2(ji,jj,jk) =                    ( zs2(ji,jj,jk) + zs2(ji  ,jj-1,jk) )   &
(  157)                   &           * ( 0.25 + SIGN( 0.25, zs2(ji,jj,jk) * zs2(ji  ,jj-1,jk) ) )
(  158)             END DO
(  159)          END DO
(  160)       END DO
(  161)       ! bottom values
(  162)       ztp1(:,:,jpk) = 0.e0    ;    ztp2(:,:,jpk) = 0.e0
(  163)       zsp1(:,:,jpk) = 0.e0    ;    zsp2(:,:,jpk) = 0.e0
(  164) 
(  165)       ! Slopes limitation
(  166)       DO jk = 1, jpkm1
(  167)          DO jj = 2, jpj
(  168)             DO ji = 2, jpi   ! vector opt.
(  169)                ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  170)                   &           * MIN(    ABS( ztp1(ji  ,jj,jk) ),   &
(  171)                   &                  2.*ABS( zt1 (ji-1,jj,jk) ),   &
(  172)                   &                  2.*ABS( zt1 (ji  ,jj,jk) ) )
(  173)                zsp1(ji,jj,jk) = SIGN( 1., zsp1(ji,jj,jk) )   &
(  174)                   &           * MIN(    ABS( zsp1(ji  ,jj,jk) ),   &
(  175)                   &                  2.*ABS( zs1 (ji-1,jj,jk) ),   &
(  176)                   &                  2.*ABS( zs1 (ji  ,jj,jk) ) )
(  177)                ztp2(ji,jj,jk) = SIGN( 1., ztp2(ji,jj,jk) )   &
(  178)                   &           * MIN(    ABS( ztp2(ji,jj  ,jk) ),   &
(  179)                   &                  2.*ABS( zt2 (ji,jj-1,jk) ),   &
(  180)                   &                  2.*ABS( zt2 (ji,jj  ,jk) ) )
(  181)                zsp2(ji,jj,jk) = SIGN( 1., zsp2(ji,jj,jk) )   &
(  182)                   &           * MIN(    ABS( zsp2(ji,jj  ,jk) ),   &
(  183)                   &                  2.*ABS( zs2 (ji,jj-1,jk) ),   &
(  184)                   &                  2.*ABS( zs2 (ji,jj  ,jk) ) )
(  185)             END DO
(  186)          END DO
(  187)       END DO        
(  188) 
(  189)       ! Advection terms
(  190)       ! interior values
(  191)       DO jk = 1, jpkm1
(  192)          zstep  = z2 * rdttra(jk)
(  193)          DO jj = 2, jpjm1      
(  194)             DO ji = 2, jpim1   ! vector opt.
(  195)                ! volume fluxes
(  196) # 200
(  200)                zeu = e2u(ji,jj) * zun(ji,jj,jk)
(  201)                zev = e1v(ji,jj) * zvn(ji,jj,jk)
(  202) # 203
(  203)                ! MUSCL fluxes
(  204)                z0u = SIGN( 0.5, zun(ji,jj,jk) )            
(  205)                zalpha = 0.5 - z0u
(  206)                zu  = z0u - 0.5 * zun(ji,jj,jk) * zstep / e1u(ji,jj)
(  207)                zzt1 = tb(ji+1,jj,jk) + zu*ztp1(ji+1,jj,jk)
(  208)                zzt2 = tb(ji  ,jj,jk) + zu*ztp1(ji  ,jj,jk)
(  209)                zzs1 = sb(ji+1,jj,jk) + zu*zsp1(ji+1,jj,jk)
(  210)                zzs2 = sb(ji  ,jj,jk) + zu*zsp1(ji  ,jj,jk)
(  211)                zt1(ji,jj,jk) = zeu * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  212)                zs1(ji,jj,jk) = zeu * ( zalpha * zzs1 + (1.-zalpha) * zzs2 )
(  213) 
(  214)                z0v = SIGN( 0.5, zvn(ji,jj,jk) )            
(  215)                zalpha = 0.5 - z0v
(  216)                zv  = z0v - 0.5 * zvn(ji,jj,jk) * zstep / e2v(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 14

(  217)                zzt1 = tb(ji,jj+1,jk) + zv*ztp2(ji,jj+1,jk)
(  218)                zzt2 = tb(ji,jj  ,jk) + zv*ztp2(ji,jj  ,jk)
(  219)                zzs1 = sb(ji,jj+1,jk) + zv*zsp2(ji,jj+1,jk)
(  220)                zzs2 = sb(ji,jj  ,jk) + zv*zsp2(ji,jj  ,jk)
(  221)                zt2(ji,jj,jk) = zev * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  222)                zs2(ji,jj,jk) = zev * ( zalpha * zzs1 + (1.-zalpha) * zzs2 )
(  223)             END DO
(  224)          END DO
(  225)       END DO
(  226) 
(  227)       !!!!  centered scheme at lateral b.C. if off-shore velocity
(  228)       DO jk = 1, jpkm1
(  229)         DO jj = 2, jpjm1
(  230)             DO ji = 2, jpim1   ! vector opt.
(  231) # 263
(  263)                IF( umask(ji,jj,jk) == 0. ) THEN
(  264)                   IF( zun(ji+1,jj,jk) > 0. .AND. ji /= jpi ) THEN
(  265)                      zt1(ji+1,jj,jk) = e2u(ji+1,jj) * zun(ji+1,jj,jk) * ( tb(ji+1,jj,jk) + tb(ji+2,jj,jk) ) * 0.5
(  266)                      zs1(ji+1,jj,jk) = e2u(ji+1,jj) * zun(ji+1,jj,jk) * ( sb(ji+1,jj,jk) + sb(ji+2,jj,jk) ) * 0.5
(  267)                   ENDIF
(  268)                   IF( zun(ji-1,jj,jk) < 0. ) THEN
(  269)                      zt1(ji-1,jj,jk) = e2u(ji-1,jj) * zun(ji-1,jj,jk) * ( tb(ji-1,jj,jk) + tb(ji  ,jj,jk) ) * 0.5
(  270)                      zs1(ji-1,jj,jk) = e2u(ji-1,jj) * zun(ji-1,jj,jk) * ( sb(ji-1,jj,jk) + sb(ji  ,jj,jk) ) * 0.5
(  271)                   ENDIF
(  272)                ENDIF
(  273)                IF( vmask(ji,jj,jk) == 0. ) THEN
(  274)                   IF( zvn(ji,jj+1,jk) > 0. .AND. jj /= jpj ) THEN
(  275)                      zt2(ji,jj+1,jk) = e1v(ji,jj+1) * zvn(ji,jj+1,jk) * ( tb(ji,jj+1,jk) + tb(ji,jj+2,jk) ) * 0.5
(  276)                      zs2(ji,jj+1,jk) = e1v(ji,jj+1) * zvn(ji,jj+1,jk) * ( sb(ji,jj+1,jk) + sb(ji,jj+2,jk) ) * 0.5
(  277)                   ENDIF
(  278)                   IF( zvn(ji,jj-1,jk) < 0. ) THEN
(  279)                      zt2(ji,jj-1,jk) = e1v(ji,jj-1) * zvn(ji,jj-1,jk) * ( tb(ji,jj-1,jk) + tb(ji  ,jj,jk) ) * 0.5
(  280)                      zs2(ji,jj-1,jk) = e1v(ji,jj-1) * zvn(ji,jj-1,jk) * ( sb(ji,jj-1,jk) + sb(ji  ,jj,jk) ) * 0.5
(  281)                   ENDIF
(  282)                ENDIF
(  283) # 284
(  284)             END DO
(  285)          END DO
(  286)       END DO
(  287) 
(  288)       ! lateral boundary conditions on zt1, zt2 ; zs1, zs2   (changed sign)
(  289)       CALL lbc_lnk( zt1, 'U', -1. )   ;   CALL lbc_lnk( zs1, 'U', -1. ) 
(  290)       CALL lbc_lnk( zt2, 'V', -1. )   ;   CALL lbc_lnk( zs2, 'V', -1. )
(  291) 
(  292)       ! Save MUSCL fluxes to compute i- & j- horizontal 
(  293)       ! advection trends in the MLD
(  294)       IF( l_trdtra )   THEN
(  295)          ! save i- terms
(  296)          tladi(:,:,:) = zt1(:,:,:) 
(  297)          sladi(:,:,:) = zs1(:,:,:) 
(  298)          ! save j- terms
(  299)          tladj(:,:,:) = zt2(:,:,:) 
(  300)          sladj(:,:,:) = zs2(:,:,:) 
(  301)       ENDIF
(  302) 
(  303)       ! Compute & add the horizontal advective trend
(  304) 
(  305)       DO jk = 1, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 15

(  306)          DO jj = 2, jpjm1      
(  307)             DO ji = 2, jpim1   ! vector opt.
(  308) # 311
(  311)                zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  312) # 313
(  313)                ! horizontal advective trends
(  314)                zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk  )   &
(  315)                   &           + zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk  ) )
(  316)                zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji-1,jj  ,jk  )   &
(  317)                   &           + zs2(ji,jj,jk) - zs2(ji  ,jj-1,jk  ) ) 
(  318)                ! add it to the general tracer trends
(  319)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  320)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  321)             END DO
(  322)         END DO
(  323)       END DO        
(  324) 
(  325)       ! Save the horizontal advective trends for diagnostic
(  326) 
(  327)       IF( l_trdtra )   THEN
(  328)          ! Recompute the hoizontal advection zta & zsa trends computed 
(  329)          ! at the step 2. above in making the difference between the new 
(  330)          ! trends and the previous one ta()/sa - ztdta()/ztdsa() and add
(  331)          ! the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
(  332)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:) 
(  333)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)
(  334) 
(  335)          CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)
(  336) 
(  337)          ! Save the new ta and sa trends
(  338)          ztdta(:,:,:) = ta(:,:,:) 
(  339)          ztdsa(:,:,:) = sa(:,:,:) 
(  340) 
(  341)       ENDIF
(  342) 
(  343)       IF(ln_ctl) THEN
(  344)          CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl2 had  - Ta: ', mask1=tmask, &
(  345)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  346)       ENDIF
(  347) 
(  348)       ! "zonal" mean advective heat and salt transport
(  349)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  350) # 354
(  354)          DO jk = 1, jpkm1
(  355)             DO jj = 2, jpjm1
(  356)                DO ji = 2, jpim1   ! vector opt.
(  357)                  zt2(ji,jj,jk) = zt2(ji,jj,jk) * e3t(jk)
(  358)                  zs2(ji,jj,jk) = zs2(ji,jj,jk) * e3t(jk)
(  359)                END DO
(  360)             END DO
(  361)          END DO
(  362)          pht_adv(:) = ptr_vj( zt2(:,:,:) )
(  363)          pst_adv(:) = ptr_vj( zs2(:,:,:) )
(  364) # 365
(  365)       ENDIF
(  366) 
(  367)       ! II. Vertical advective fluxes
(  368)       ! -----------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 16

(  369)       
(  370)       ! First guess of the slope
(  371)       ! interior values
(  372)       DO jk = 2, jpkm1
(  373)          zt1(:,:,jk) = tmask(:,:,jk) * ( tb(:,:,jk-1) - tb(:,:,jk) )
(  374)          zs1(:,:,jk) = tmask(:,:,jk) * ( sb(:,:,jk-1) - sb(:,:,jk) )
(  375)       END DO
(  376)       ! surface & bottom boundary conditions
(  377)       zt1 (:,:, 1 ) = 0.e0    ;    zt1 (:,:,jpk) = 0.e0
(  378)       zs1 (:,:, 1 ) = 0.e0    ;    zs1 (:,:,jpk) = 0.e0
(  379) 
(  380)       ! Slopes
(  381)       DO jk = 2, jpkm1
(  382)          DO jj = 1, jpj
(  383)             DO ji = 1, jpi
(  384)                ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji,jj,jk+1) )   &
(  385)                   &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji,jj,jk+1) ) )
(  386)                zsp1(ji,jj,jk) =                    ( zs1(ji,jj,jk) + zs1(ji,jj,jk+1) )   &
(  387)                   &           * ( 0.25 + SIGN( 0.25, zs1(ji,jj,jk) * zs1(ji,jj,jk+1) ) )
(  388)             END DO
(  389)          END DO
(  390)       END DO
(  391) 
(  392)       ! Slopes limitation
(  393)       ! interior values
(  394)       DO jk = 2, jpkm1
(  395)          DO jj = 1, jpj
(  396)             DO ji = 1, jpi
(  397)                ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  398)                   &           * MIN(    ABS( ztp1(ji,jj,jk  ) ),   &
(  399)                   &                  2.*ABS( zt1 (ji,jj,jk+1) ),   &
(  400)                   &                  2.*ABS( zt1 (ji,jj,jk  ) ) )
(  401)                zsp1(ji,jj,jk) = SIGN( 1., zsp1(ji,jj,jk) )   &
(  402)                   &           * MIN(    ABS( zsp1(ji,jj,jk  ) ),   &
(  403)                   &                  2.*ABS( zs1 (ji,jj,jk+1) ),   &
(  404)                   &                  2.*ABS( zs1 (ji,jj,jk  ) ) )
(  405)             END DO
(  406)          END DO
(  407)       END DO
(  408)       ! surface values
(  409)       ztp1(:,:,1) = 0.e0
(  410)       zsp1(:,:,1) = 0.e0
(  411) 
(  412)       ! vertical advective flux
(  413)       ! interior values
(  414)       DO jk = 1, jpkm1
(  415)          zstep  = z2 * rdttra(jk)
(  416)          DO jj = 2, jpjm1      
(  417)             DO ji = 2, jpim1   ! vector opt.
(  418)                zew = zwn(ji,jj,jk+1)
(  419)                z0w = SIGN( 0.5, zwn(ji,jj,jk+1) )
(  420)                zalpha = 0.5 + z0w
(  421)                zw  = z0w - 0.5 * zwn(ji,jj,jk+1)*zstep / e3w(jk+1)
(  422)                zzt1 = tb(ji,jj,jk+1) + zw*ztp1(ji,jj,jk+1)
(  423)                zzt2 = tb(ji,jj,jk  ) + zw*ztp1(ji,jj,jk  )
(  424)                zzs1 = sb(ji,jj,jk+1) + zw*zsp1(ji,jj,jk+1)
(  425)                zzs2 = sb(ji,jj,jk  ) + zw*zsp1(ji,jj,jk  )
(  426)                zt1(ji,jj,jk+1) = zew * ( zalpha * zzt1 + (1.-zalpha)*zzt2 )






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 17

(  427)                zs1(ji,jj,jk+1) = zew * ( zalpha * zzs1 + (1.-zalpha)*zzs2 )
(  428)             END DO
(  429)          END DO
(  430)       END DO
(  431)       DO jk = 2, jpkm1
(  432)         DO jj = 2, jpjm1
(  433)             DO ji = 2, jpim1   ! vector opt.
(  434)                IF( tmask(ji,jj,jk+1) == 0. ) THEN
(  435)                   IF( zwn(ji,jj,jk) > 0. ) THEN
(  436)                      zt1(ji,jj,jk) = zwn(ji,jj,jk) * ( tb(ji,jj,jk-1) + tb(ji,jj,jk) ) * 0.5
(  437)                      zs1(ji,jj,jk) = zwn(ji,jj,jk) * ( sb(ji,jj,jk-1) + sb(ji,jj,jk) ) * 0.5
(  438)                   ENDIF
(  439)                ENDIF
(  440)             END DO
(  441)          END DO
(  442)       END DO
(  443) 
(  444)       ! surface values
(  445)       IF( lk_dynspg_rl ) THEN                           ! rigid lid : flux set to zero
(  446)          zt1(:,:, 1 ) = 0.e0
(  447)          zs1(:,:, 1 ) = 0.e0
(  448)       ELSE                                              ! free surface
(  449)          zt1(:,:, 1 ) = zwn(:,:,1) * tb(:,:,1)
(  450)          zs1(:,:, 1 ) = zwn(:,:,1) * sb(:,:,1)
(  451)       ENDIF
(  452) 
(  453)       ! bottom values
(  454)       zt1(:,:,jpk) = 0.e0
(  455)       zs1(:,:,jpk) = 0.e0
(  456) 
(  457) 
(  458)       ! Compute & add the vertical advective trend
(  459) 
(  460)       DO jk = 1, jpkm1
(  461)          DO jj = 2, jpjm1      
(  462)             DO ji = 2, jpim1   ! vector opt.
(  463)                zbtr = 1. / e3t(jk)
(  464)                ! horizontal advective trends
(  465)                zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji,jj,jk+1) )
(  466)                zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji,jj,jk+1) )
(  467)                ! add it to the general tracer trends
(  468)                ta(ji,jj,jk) =  ta(ji,jj,jk) + zta
(  469)                sa(ji,jj,jk) =  sa(ji,jj,jk) + zsa
(  470)             END DO
(  471)          END DO
(  472)       END DO
(  473) 
(  474)       ! Save the vertical advective trends for diagnostic
(  475) 
(  476)       IF( l_trdtra )   THEN
(  477)          ! Recompute the vertical advection zta & zsa trends computed 
(  478)          ! at the step 2. above in making the difference between the new 
(  479)          ! trends and the previous one: ta()/sa - ztdta()/ztdsa() and substract
(  480)          ! the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
(  481)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tn(:,:,:) * hdivn(:,:,:) 
(  482)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  483) 
(  484)          CALL trd_mod(ztdta, ztdsa, jpttdzad, 'TRA', kt)






PGF90 (Version     10.5)          09/23/2010  09:50:11      page 18

(  485)       ENDIF
(  486) 
(  487)       IF(ln_ctl) THEN
(  488)          CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl2 zad  - Ta: ', mask1=tmask, &
(  489)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  490) 
(  491)       ENDIF
(  492) 
(  493)    END SUBROUTINE tra_adv_muscl2
(  494) 
(  495)    !!======================================================================
(  496) END MODULE traadv_muscl2
