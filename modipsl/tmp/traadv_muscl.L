


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

Filename: traadv_muscl.F90

(    1) MODULE traadv_muscl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_muscl  ***
(    4)    !! Ocean active tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_muscl : update the tracer trend with the horizontal
(    9)    !!                   and vertical advection trends using MUSCL scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE trdmod          ! ocean active tracers trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   18)    USE trabbl          ! tracers: bottom boundary layer
(   19)    USE lib_mpp         ! distribued memory computing
(   20)    USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
(   21)    USE diaptr          ! poleward transport diagnostics
(   22)    USE prtctl          ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Accessibility
(   28)    PUBLIC tra_adv_muscl  ! routine called by step.F90
(   29) 
(   30)    !! * Substitutions
(   31) #  include "domzgr_substitute.h90"
(   32) #  include "vectopt_loop_substitute.h90"
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_muscl.F90,v 1.10 2005/12/28 09:25:10 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE tra_adv_muscl( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                    ***  ROUTINE tra_adv_muscl  ***
(   44)       !!
(   45)       !! ** Purpose :   Compute the now trend due to total advection of T and 
(   46)       !!      S using a MUSCL scheme (Monotone Upstream-centered Scheme for






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 2

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
(   89)          zt1, zt2, ztp1, ztp2,               &
(   90)          zs1, zs2, zsp1, zsp2,               &
(   91)          ztdta, ztdsa
(   92)       !!----------------------------------------------------------------------
(   93) 
(   94)       IF( kt == nit000 .AND. lwp ) THEN
(   95)          WRITE(numout,*)
(   96)          WRITE(numout,*) 'tra_adv : MUSCL advection scheme'
(   97)          WRITE(numout,*) '~~~~~~~'
(   98)       ENDIF
(   99) 
(  100)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  101)           z2=1.
(  102)       ELSE
(  103)           z2=2.
(  104)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 3

(  105) 
(  106)       ! Save ta and sa trends
(  107)       IF( l_trdtra )   THEN
(  108)          ztdta(:,:,:) = ta(:,:,:) 
(  109)          ztdsa(:,:,:) = sa(:,:,:) 
(  110)          l_adv = 'mus'
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






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 4

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






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 5

(  221)                zt2(ji,jj,jk) = zev * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  222)                zs2(ji,jj,jk) = zev * ( zalpha * zzs1 + (1.-zalpha) * zzs2 )
(  223)             END DO
(  224)          END DO
(  225)       END DO
(  226) 
(  227)       ! lateral boundary conditions on zt1, zt2 ; zs1, zs2   (changed sign)
(  228)       CALL lbc_lnk( zt1, 'U', -1. )   ;   CALL lbc_lnk( zs1, 'U', -1. ) 
(  229)       CALL lbc_lnk( zt2, 'V', -1. )   ;   CALL lbc_lnk( zs2, 'V', -1. )
(  230) 
(  231)       ! Save MUSCL fluxes to compute i- & j- horizontal 
(  232)       ! advection trends in the MLD
(  233)       IF( l_trdtra )   THEN
(  234)          ! save i- terms
(  235)          tladi(:,:,:) = zt1(:,:,:) 
(  236)          sladi(:,:,:) = zs1(:,:,:) 
(  237)          ! save j- terms
(  238)          tladj(:,:,:) = zt2(:,:,:) 
(  239)          sladj(:,:,:) = zs2(:,:,:) 
(  240)       ENDIF
(  241) 
(  242)       ! Compute & add the horizontal advective trend
(  243) 
(  244)       DO jk = 1, jpkm1
(  245)          DO jj = 2, jpjm1      
(  246)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  247) #if defined key_s_coord || defined key_partial_steps
(  248)                zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  249) #else
(  250)                zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  251) #endif
(  252)                ! horizontal advective trends
(  253)                zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk  )   &
(  254)                   &           + zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk  ) )
(  255)                zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji-1,jj  ,jk  )   &
(  256)                   &           + zs2(ji,jj,jk) - zs2(ji  ,jj-1,jk  ) ) 
(  257)                ! add it to the general tracer trends
(  258)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  259)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  260)             END DO
(  261)         END DO
(  262)       END DO        
(  263) 
(  264)       ! Save the horizontal advective trends for diagnostic
(  265) 
(  266)       IF( l_trdtra )   THEN
(  267)          ! Recompute the horizontal advection zta & zsa trends computed 
(  268)          ! at the step 2. above in making the difference between the new 
(  269)          ! trends and the previous one ta()/sa - ztdta()/ztdsa() and add
(  270)          ! the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
(  271)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:) 
(  272)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)
(  273) 
(  274)          CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)
(  275) 
(  276)          ! Save the new ta and sa trends
(  277)          ztdta(:,:,:) = ta(:,:,:) 
(  278)          ztdsa(:,:,:) = sa(:,:,:) 






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 6

(  279) 
(  280)       ENDIF
(  281) 
(  282)       IF(ln_ctl) THEN
(  283)          CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl had  - Ta: ', mask1=tmask ,  &
(  284)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra' )
(  285)       ENDIF
(  286) 
(  287)       ! "zonal" mean advective heat and salt transport
(  288)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  289) # if defined key_s_coord || defined key_partial_steps
(  290)          pht_adv(:) = ptr_vj( zt2(:,:,:) )
(  291)          pst_adv(:) = ptr_vj( zs2(:,:,:) )
(  292) # else
(  293)          DO jk = 1, jpkm1
(  294)             DO jj = 2, jpjm1
(  295)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  296)                  zt2(ji,jj,jk) = zt2(ji,jj,jk) * fse3v(ji,jj,jk)
(  297)                  zs2(ji,jj,jk) = zs2(ji,jj,jk) * fse3v(ji,jj,jk)
(  298)                END DO
(  299)             END DO
(  300)          END DO
(  301)          pht_adv(:) = ptr_vj( zt2(:,:,:) )
(  302)          pst_adv(:) = ptr_vj( zs2(:,:,:) )
(  303) # endif
(  304)       ENDIF
(  305) 
(  306)       ! II. Vertical advective fluxes
(  307)       ! -----------------------------
(  308)       
(  309)       ! First guess of the slope
(  310)       ! interior values
(  311)       DO jk = 2, jpkm1
(  312)          zt1(:,:,jk) = tmask(:,:,jk) * ( tb(:,:,jk-1) - tb(:,:,jk) )
(  313)          zs1(:,:,jk) = tmask(:,:,jk) * ( sb(:,:,jk-1) - sb(:,:,jk) )
(  314)       END DO
(  315)       ! surface & bottom boundary conditions
(  316)       zt1 (:,:, 1 ) = 0.e0    ;    zt1 (:,:,jpk) = 0.e0
(  317)       zs1 (:,:, 1 ) = 0.e0    ;    zs1 (:,:,jpk) = 0.e0
(  318) 
(  319)       ! Slopes
(  320)       DO jk = 2, jpkm1
(  321)          DO jj = 1, jpj
(  322)             DO ji = 1, jpi
(  323)                ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji,jj,jk+1) )   &
(  324)                   &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji,jj,jk+1) ) )
(  325)                zsp1(ji,jj,jk) =                    ( zs1(ji,jj,jk) + zs1(ji,jj,jk+1) )   &
(  326)                   &           * ( 0.25 + SIGN( 0.25, zs1(ji,jj,jk) * zs1(ji,jj,jk+1) ) )
(  327)             END DO
(  328)          END DO
(  329)       END DO
(  330) 
(  331)       ! Slopes limitation
(  332)       ! interior values
(  333)       DO jk = 2, jpkm1
(  334)          DO jj = 1, jpj
(  335)             DO ji = 1, jpi
(  336)                ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 7

(  337)                   &           * MIN(    ABS( ztp1(ji,jj,jk  ) ),   &
(  338)                   &                  2.*ABS( zt1 (ji,jj,jk+1) ),   &
(  339)                   &                  2.*ABS( zt1 (ji,jj,jk  ) ) )
(  340)                zsp1(ji,jj,jk) = SIGN( 1., zsp1(ji,jj,jk) )   &
(  341)                   &           * MIN(    ABS( zsp1(ji,jj,jk  ) ),   &
(  342)                   &                  2.*ABS( zs1 (ji,jj,jk+1) ),   &
(  343)                   &                  2.*ABS( zs1 (ji,jj,jk  ) ) )
(  344)             END DO
(  345)          END DO
(  346)       END DO
(  347)       ! surface values
(  348)       ztp1(:,:,1) = 0.e0
(  349)       zsp1(:,:,1) = 0.e0
(  350) 
(  351)       ! vertical advective flux
(  352)       ! interior values
(  353)       DO jk = 1, jpkm1
(  354)          zstep  = z2 * rdttra(jk)
(  355)          DO jj = 2, jpjm1      
(  356)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  357)                zew = zwn(ji,jj,jk+1)
(  358)                z0w = SIGN( 0.5, zwn(ji,jj,jk+1) )
(  359)                zalpha = 0.5 + z0w
(  360)                zw  = z0w - 0.5 * zwn(ji,jj,jk+1)*zstep / fse3w(ji,jj,jk+1)
(  361)                zzt1 = tb(ji,jj,jk+1) + zw*ztp1(ji,jj,jk+1)
(  362)                zzt2 = tb(ji,jj,jk  ) + zw*ztp1(ji,jj,jk  )
(  363)                zzs1 = sb(ji,jj,jk+1) + zw*zsp1(ji,jj,jk+1)
(  364)                zzs2 = sb(ji,jj,jk  ) + zw*zsp1(ji,jj,jk  )
(  365)                zt1(ji,jj,jk+1) = zew * ( zalpha * zzt1 + (1.-zalpha)*zzt2 )
(  366)                zs1(ji,jj,jk+1) = zew * ( zalpha * zzs1 + (1.-zalpha)*zzs2 )
(  367)             END DO
(  368)          END DO
(  369)       END DO
(  370)       ! surface values
(  371)       IF( lk_dynspg_rl ) THEN                           ! rigid lid : flux set to zero
(  372)          zt1(:,:, 1 ) = 0.e0
(  373)          zs1(:,:, 1 ) = 0.e0
(  374)       ELSE                                              ! free surface
(  375)          zt1(:,:, 1 ) = zwn(:,:,1) * tb(:,:,1)
(  376)          zs1(:,:, 1 ) = zwn(:,:,1) * sb(:,:,1)
(  377)       ENDIF
(  378) 
(  379)       ! bottom values
(  380)       zt1(:,:,jpk) = 0.e0
(  381)       zs1(:,:,jpk) = 0.e0
(  382) 
(  383) 
(  384)       ! Compute & add the vertical advective trend
(  385) 
(  386)       DO jk = 1, jpkm1
(  387)          DO jj = 2, jpjm1      
(  388)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  389)                zbtr = 1. / fse3t(ji,jj,jk)
(  390)                ! horizontal advective trends
(  391)                zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji,jj,jk+1) )
(  392)                zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji,jj,jk+1) )
(  393)                ! add it to the general tracer trends
(  394)                ta(ji,jj,jk) =  ta(ji,jj,jk) + zta






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 8

(  395)                sa(ji,jj,jk) =  sa(ji,jj,jk) + zsa
(  396)             END DO
(  397)          END DO
(  398)       END DO
(  399) 
(  400)       ! Save the vertical advective trends for diagnostic
(  401) 
(  402)       IF( l_trdtra )   THEN
(  403)          ! Recompute the vertical advection zta & zsa trends computed 
(  404)          ! at the step 2. above in making the difference between the new 
(  405)          ! trends and the previous one: ta()/sa - ztdta()/ztdsa() and substract
(  406)          ! the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
(  407)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tn(:,:,:) * hdivn(:,:,:) 
(  408)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  409) 
(  410)          CALL trd_mod(ztdta, ztdsa, jpttdzad, 'TRA', kt)
(  411)       ENDIF
(  412) 
(  413)       IF(ln_ctl) THEN
(  414)          CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl zad  - Ta: ', mask1=tmask ,  &
(  415)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  416)       ENDIF
(  417) 
(  418)    END SUBROUTINE tra_adv_muscl
(  419) 
(  420)    !!======================================================================
(  421) END MODULE traadv_muscl





































PGF90 (Version     10.5)          09/23/2010  09:50:12      page 9

(    1) # 1 "traadv_muscl.F90"
(    1) MODULE traadv_muscl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_muscl  ***
(    4)    !! Ocean active tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_muscl : update the tracer trend with the horizontal
(    9)    !!                   and vertical advection trends using MUSCL scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE trdmod          ! ocean active tracers trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   18)    USE trabbl          ! tracers: bottom boundary layer
(   19)    USE lib_mpp         ! distribued memory computing
(   20)    USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
(   21)    USE diaptr          ! poleward transport diagnostics
(   22)    USE prtctl          ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Accessibility
(   28)    PUBLIC tra_adv_muscl  ! routine called by step.F90
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
(   48) # 32 "traadv_muscl.F90"
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






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 10

(   11)    !!----------------------------------------------------------------------
(   12) # 33 "traadv_muscl.F90"
(   33) # 33
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_muscl.F90,v 1.10 2005/12/28 09:25:10 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE tra_adv_muscl( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                    ***  ROUTINE tra_adv_muscl  ***
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
(   89)          zt1, zt2, ztp1, ztp2,               &
(   90)          zs1, zs2, zsp1, zsp2,               &
(   91)          ztdta, ztdsa






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 11

(   92)       !!----------------------------------------------------------------------
(   93) 
(   94)       IF( kt == nit000 .AND. lwp ) THEN
(   95)          WRITE(numout,*)
(   96)          WRITE(numout,*) 'tra_adv : MUSCL advection scheme'
(   97)          WRITE(numout,*) '~~~~~~~'
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
(  110)          l_adv = 'mus'
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






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 12

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






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 13

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
(  227)       ! lateral boundary conditions on zt1, zt2 ; zs1, zs2   (changed sign)
(  228)       CALL lbc_lnk( zt1, 'U', -1. )   ;   CALL lbc_lnk( zs1, 'U', -1. ) 
(  229)       CALL lbc_lnk( zt2, 'V', -1. )   ;   CALL lbc_lnk( zs2, 'V', -1. )
(  230) 
(  231)       ! Save MUSCL fluxes to compute i- & j- horizontal 
(  232)       ! advection trends in the MLD
(  233)       IF( l_trdtra )   THEN
(  234)          ! save i- terms
(  235)          tladi(:,:,:) = zt1(:,:,:) 
(  236)          sladi(:,:,:) = zs1(:,:,:) 
(  237)          ! save j- terms
(  238)          tladj(:,:,:) = zt2(:,:,:) 
(  239)          sladj(:,:,:) = zs2(:,:,:) 
(  240)       ENDIF
(  241) 
(  242)       ! Compute & add the horizontal advective trend
(  243) 
(  244)       DO jk = 1, jpkm1
(  245)          DO jj = 2, jpjm1      
(  246)             DO ji = 2, jpim1   ! vector opt.
(  247) # 250
(  250)                zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  251) # 252
(  252)                ! horizontal advective trends
(  253)                zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk  )   &
(  254)                   &           + zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk  ) )
(  255)                zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji-1,jj  ,jk  )   &
(  256)                   &           + zs2(ji,jj,jk) - zs2(ji  ,jj-1,jk  ) ) 
(  257)                ! add it to the general tracer trends
(  258)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  259)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  260)             END DO
(  261)         END DO
(  262)       END DO        
(  263) 
(  264)       ! Save the horizontal advective trends for diagnostic
(  265) 
(  266)       IF( l_trdtra )   THEN
(  267)          ! Recompute the horizontal advection zta & zsa trends computed 
(  268)          ! at the step 2. above in making the difference between the new 
(  269)          ! trends and the previous one ta()/sa - ztdta()/ztdsa() and add
(  270)          ! the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
(  271)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:) 
(  272)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)
(  273) 
(  274)          CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)
(  275) 
(  276)          ! Save the new ta and sa trends






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 14

(  277)          ztdta(:,:,:) = ta(:,:,:) 
(  278)          ztdsa(:,:,:) = sa(:,:,:) 
(  279) 
(  280)       ENDIF
(  281) 
(  282)       IF(ln_ctl) THEN
(  283)          CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl had  - Ta: ', mask1=tmask ,  &
(  284)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra' )
(  285)       ENDIF
(  286) 
(  287)       ! "zonal" mean advective heat and salt transport
(  288)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  289) # 293
(  293)          DO jk = 1, jpkm1
(  294)             DO jj = 2, jpjm1
(  295)                DO ji = 2, jpim1   ! vector opt.
(  296)                  zt2(ji,jj,jk) = zt2(ji,jj,jk) * e3t(jk)
(  297)                  zs2(ji,jj,jk) = zs2(ji,jj,jk) * e3t(jk)
(  298)                END DO
(  299)             END DO
(  300)          END DO
(  301)          pht_adv(:) = ptr_vj( zt2(:,:,:) )
(  302)          pst_adv(:) = ptr_vj( zs2(:,:,:) )
(  303) # 304
(  304)       ENDIF
(  305) 
(  306)       ! II. Vertical advective fluxes
(  307)       ! -----------------------------
(  308)       
(  309)       ! First guess of the slope
(  310)       ! interior values
(  311)       DO jk = 2, jpkm1
(  312)          zt1(:,:,jk) = tmask(:,:,jk) * ( tb(:,:,jk-1) - tb(:,:,jk) )
(  313)          zs1(:,:,jk) = tmask(:,:,jk) * ( sb(:,:,jk-1) - sb(:,:,jk) )
(  314)       END DO
(  315)       ! surface & bottom boundary conditions
(  316)       zt1 (:,:, 1 ) = 0.e0    ;    zt1 (:,:,jpk) = 0.e0
(  317)       zs1 (:,:, 1 ) = 0.e0    ;    zs1 (:,:,jpk) = 0.e0
(  318) 
(  319)       ! Slopes
(  320)       DO jk = 2, jpkm1
(  321)          DO jj = 1, jpj
(  322)             DO ji = 1, jpi
(  323)                ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji,jj,jk+1) )   &
(  324)                   &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji,jj,jk+1) ) )
(  325)                zsp1(ji,jj,jk) =                    ( zs1(ji,jj,jk) + zs1(ji,jj,jk+1) )   &
(  326)                   &           * ( 0.25 + SIGN( 0.25, zs1(ji,jj,jk) * zs1(ji,jj,jk+1) ) )
(  327)             END DO
(  328)          END DO
(  329)       END DO
(  330) 
(  331)       ! Slopes limitation
(  332)       ! interior values
(  333)       DO jk = 2, jpkm1
(  334)          DO jj = 1, jpj
(  335)             DO ji = 1, jpi
(  336)                ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  337)                   &           * MIN(    ABS( ztp1(ji,jj,jk  ) ),   &






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 15

(  338)                   &                  2.*ABS( zt1 (ji,jj,jk+1) ),   &
(  339)                   &                  2.*ABS( zt1 (ji,jj,jk  ) ) )
(  340)                zsp1(ji,jj,jk) = SIGN( 1., zsp1(ji,jj,jk) )   &
(  341)                   &           * MIN(    ABS( zsp1(ji,jj,jk  ) ),   &
(  342)                   &                  2.*ABS( zs1 (ji,jj,jk+1) ),   &
(  343)                   &                  2.*ABS( zs1 (ji,jj,jk  ) ) )
(  344)             END DO
(  345)          END DO
(  346)       END DO
(  347)       ! surface values
(  348)       ztp1(:,:,1) = 0.e0
(  349)       zsp1(:,:,1) = 0.e0
(  350) 
(  351)       ! vertical advective flux
(  352)       ! interior values
(  353)       DO jk = 1, jpkm1
(  354)          zstep  = z2 * rdttra(jk)
(  355)          DO jj = 2, jpjm1      
(  356)             DO ji = 2, jpim1   ! vector opt.
(  357)                zew = zwn(ji,jj,jk+1)
(  358)                z0w = SIGN( 0.5, zwn(ji,jj,jk+1) )
(  359)                zalpha = 0.5 + z0w
(  360)                zw  = z0w - 0.5 * zwn(ji,jj,jk+1)*zstep / e3w(jk+1)
(  361)                zzt1 = tb(ji,jj,jk+1) + zw*ztp1(ji,jj,jk+1)
(  362)                zzt2 = tb(ji,jj,jk  ) + zw*ztp1(ji,jj,jk  )
(  363)                zzs1 = sb(ji,jj,jk+1) + zw*zsp1(ji,jj,jk+1)
(  364)                zzs2 = sb(ji,jj,jk  ) + zw*zsp1(ji,jj,jk  )
(  365)                zt1(ji,jj,jk+1) = zew * ( zalpha * zzt1 + (1.-zalpha)*zzt2 )
(  366)                zs1(ji,jj,jk+1) = zew * ( zalpha * zzs1 + (1.-zalpha)*zzs2 )
(  367)             END DO
(  368)          END DO
(  369)       END DO
(  370)       ! surface values
(  371)       IF( lk_dynspg_rl ) THEN                           ! rigid lid : flux set to zero
(  372)          zt1(:,:, 1 ) = 0.e0
(  373)          zs1(:,:, 1 ) = 0.e0
(  374)       ELSE                                              ! free surface
(  375)          zt1(:,:, 1 ) = zwn(:,:,1) * tb(:,:,1)
(  376)          zs1(:,:, 1 ) = zwn(:,:,1) * sb(:,:,1)
(  377)       ENDIF
(  378) 
(  379)       ! bottom values
(  380)       zt1(:,:,jpk) = 0.e0
(  381)       zs1(:,:,jpk) = 0.e0
(  382) 
(  383) 
(  384)       ! Compute & add the vertical advective trend
(  385) 
(  386)       DO jk = 1, jpkm1
(  387)          DO jj = 2, jpjm1      
(  388)             DO ji = 2, jpim1   ! vector opt.
(  389)                zbtr = 1. / e3t(jk)
(  390)                ! horizontal advective trends
(  391)                zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji,jj,jk+1) )
(  392)                zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji,jj,jk+1) )
(  393)                ! add it to the general tracer trends
(  394)                ta(ji,jj,jk) =  ta(ji,jj,jk) + zta
(  395)                sa(ji,jj,jk) =  sa(ji,jj,jk) + zsa






PGF90 (Version     10.5)          09/23/2010  09:50:12      page 16

(  396)             END DO
(  397)          END DO
(  398)       END DO
(  399) 
(  400)       ! Save the vertical advective trends for diagnostic
(  401) 
(  402)       IF( l_trdtra )   THEN
(  403)          ! Recompute the vertical advection zta & zsa trends computed 
(  404)          ! at the step 2. above in making the difference between the new 
(  405)          ! trends and the previous one: ta()/sa - ztdta()/ztdsa() and substract
(  406)          ! the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
(  407)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tn(:,:,:) * hdivn(:,:,:) 
(  408)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  409) 
(  410)          CALL trd_mod(ztdta, ztdsa, jpttdzad, 'TRA', kt)
(  411)       ENDIF
(  412) 
(  413)       IF(ln_ctl) THEN
(  414)          CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl zad  - Ta: ', mask1=tmask ,  &
(  415)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  416)       ENDIF
(  417) 
(  418)    END SUBROUTINE tra_adv_muscl
(  419) 
(  420)    !!======================================================================
(  421) END MODULE traadv_muscl
