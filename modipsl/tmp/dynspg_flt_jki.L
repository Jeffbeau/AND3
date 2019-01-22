


PGF90 (Version     10.5)          09/23/2010  09:49:42      page 1

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

Filename: dynspg_flt_jki.F90

(    1) MODULE dynspg_flt_jki
(    2)    !!======================================================================
(    3)    !!                  ***  MODULE  dynspg_flt_jki  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) #if ( defined key_dynspg_flt && defined key_autotasking )   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_dynspg_flt' & 'key_autotasking'          filtered free surface
(    9)    !!                                                   j-k-i loop (j-slab)
(   10)    !!----------------------------------------------------------------------
(   11)    !!   dyn_spg_flt_jki : Update the momentum trend with the surface pressure
(   12)    !!                     gradient for the free surf. constant volume case
(   13)    !!                     with auto-tasking (j-slab) (no vectior opt.)
(   14)    !!----------------------------------------------------------------------
(   15)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   16)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_flt_jki.F90,v 1.3 2006/03/20 17:25:49 opalod Exp $ 
(   17)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   18)    !!----------------------------------------------------------------------
(   19)    !! * Modules used
(   20)    USE oce             ! ocean dynamics and tracers 
(   21)    USE dom_oce         ! ocean space and time domain
(   22)    USE zdf_oce         ! ocean vertical physics
(   23)    USE phycst          ! physical constant
(   24)    USE ocesbc          ! Ocean Surface Boundary condition
(   25)    USE flxrnf          ! ocean runoffs
(   26)    USE sol_oce         ! ocean elliptic solver
(   27)    USE solpcg          ! preconditionned conjugate gradient solver
(   28)    USE solsor          ! Successive Over-relaxation solver
(   29)    USE solfet          ! FETI solver
(   30)    USE solsor_e        ! Successive Over-relaxation solver with MPP optimization
(   31)    USE obc_oce         ! Lateral open boundary condition
(   32)    USE obcdyn          ! ocean open boundary condition (obc_dyn routines)
(   33)    USE obcvol          ! ocean open boundary condition (obc_vol routines)
(   34)    USE lib_mpp         ! distributed memory computing library
(   35)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   36)    USE cla_dynspg      ! cross land advection
(   37)    USE prtctl          ! Print control
(   38)    USE in_out_manager  ! I/O manager
(   39)    USE solmat          ! matrix construction for elliptic solvers
(   40)    USE agrif_opa_interp
(   41) 
(   42)    IMPLICIT NONE
(   43)    PRIVATE
(   44) 
(   45)    !! * Accessibility
(   46)    PUBLIC dyn_spg_flt_jki        ! routine called by step.F90






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 2

(   47) 
(   48)    !! * Substitutions
(   49) #  include "domzgr_substitute.h90"
(   50)    !!----------------------------------------------------------------------
(   51)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   52)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_flt_jki.F90,v 1.3 2006/03/20 17:25:49 opalod Exp $ 
(   53)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   54)    !!----------------------------------------------------------------------
(   55) 
(   56) CONTAINS
(   57) 
(   58)    SUBROUTINE dyn_spg_flt_jki( kt, kindic )
(   59)       !!----------------------------------------------------------------------
(   60)       !!                  ***  routine dyn_spg_flt_jki  ***
(   61)       !!
(   62)       !! ** Purpose :   Compute the now trend due to the surface pressure 
(   63)       !!      gradient for free surface formulation with a constant ocean
(   64)       !!      volume case, add it to the general trend of momentum equation.
(   65)       !!
(   66)       !! ** Method  :   Free surface formulation. The surface pressure gradient
(   67)       !!      is given by:
(   68)       !!         spgu = 1/rau0 d/dx(ps) =  1/e1u di( etn + rnu btda )
(   69)       !!         spgv = 1/rau0 d/dy(ps) =  1/e2v dj( etn + rnu btda )
(   70)       !!      where etn is the free surface elevation and btda is the after
(   71)       !!      of the free surface elevation
(   72)       !!       -1- compute the after sea surface elevation from the cinematic
(   73)       !!      surface boundary condition:
(   74)       !!              zssha = sshb + 2 rdt ( wn - emp )
(   75)       !!           Time filter applied on now sea surface elevation to avoid 
(   76)       !!      the divergence of two consecutive time-steps and swap of free
(   77)       !!      surface arrays to start the next time step:
(   78)       !!              sshb = sshn + atfp * [ sshb + zssha - 2 sshn ]
(   79)       !!              sshn = zssha
(   80)       !!       -2- evaluate the surface presure trend (including the addi-
(   81)       !!      tional force) in three steps:
(   82)       !!        a- compute the right hand side of the elliptic equation:
(   83)       !!            gcb = 1/(e1t e2t) [ di(e2u spgu) + dj(e1v spgv) ]
(   84)       !!         where (spgu,spgv) are given by:
(   85)       !!            spgu = vertical sum[ e3u (ub+ 2 rdt ua ) ]
(   86)       !!                 - grav 2 rdt hu /e1u di[sshn + emp]
(   87)       !!            spgv = vertical sum[ e3v (vb+ 2 rdt va) ]
(   88)       !!                 - grav 2 rdt hv /e2v dj[sshn + emp]
(   89)       !!         and define the first guess from previous computation :
(   90)       !!            zbtd = btda
(   91)       !!            btda = 2 zbtd - btdb
(   92)       !!            btdb = zbtd
(   93)       !!        b- compute the relative accuracy to be reached by the
(   94)       !!         iterative solver
(   95)       !!        c- apply the solver by a call to sol... routine
(   96)       !!       -3- compute and add the free surface pressure gradient inclu-
(   97)       !!      ding the additional force used to stabilize the equation.
(   98)       !!      several slabs used ('key-autotasking')
(   99)       !!
(  100)       !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
(  101)       !!
(  102)       !! References :
(  103)       !!      Roullet and Madec 1999, JGR.
(  104)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 3

(  105)       !! History :
(  106)       !!        !  98-05 (G. Roullet)  Original code
(  107)       !!        !  98-10 (G. Madec, M. Imbard)  release 8.2
(  108)       !!   8.5  !  02-08 (G. Madec)  F90: Free form and module
(  109)       !!        !  02-11 (C. Talandier, A-M Treguier) Open boundaries
(  110)       !!   9.0  !  04-08 (C. Talandier) New trends organization
(  111)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  112)       !!---------------------------------------------------------------------
(  113)       !! * Arguments
(  114)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(  115)       INTEGER, INTENT( out ) ::   kindic     ! solver convergence flag
(  116)                                              ! if the solver doesn t converge
(  117)                                              ! the flag is < 0
(  118)       !! * Local declarations
(  119)       INTEGER  ::   ji, jj, jk               ! dummy loop indices
(  120)       REAL(wp) ::   &              ! temporary scalars
(  121)          z2dt, z2dtg, zraur, znugdt, znurau,   &
(  122)          zssha, zgcb, zbtd, ztdgu, ztdgv, zgwgt
(  123)       !!----------------------------------------------------------------------
(  124) 
(  125)       IF( kt == nit000 ) THEN
(  126)          IF(lwp) WRITE(numout,*)
(  127)          IF(lwp) WRITE(numout,*) 'dyn_spg_flt_jki : surface pressure gradient trend'
(  128)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~  (free surface constant volume, autotasking case)'
(  129) 
(  130)          ! set to zero free surface specific arrays 
(  131)          spgu(:,:) = 0.e0      ! surface pressure gradient (i-direction)
(  132)          spgv(:,:) = 0.e0      ! surface pressure gradient (j-direction)
(  133)       ENDIF
(  134) 
(  135)       ! 0. Local constant initialization
(  136)       ! --------------------------------
(  137)       ! time step: leap-frog
(  138)       z2dt = 2. * rdt
(  139)       ! time step: Euler if restart from rest
(  140)       IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt
(  141)       IF( neuler == 0 .AND. kt == nit000+1 ) CALL sol_mat(kt)
(  142)       ! coefficients
(  143)       z2dtg  = grav * z2dt
(  144)       zraur  = 1. / rauw
(  145)       znugdt =  rnu * grav * z2dt
(  146)       znurau =  znugdt * zraur
(  147) 
(  148)       !                                                ! ===============
(  149)       DO jj = 2, jpjm1                                 !  Vertical slab
(  150)          !                                             ! ===============
(  151)          ! Surface pressure gradient (now)
(  152)          DO ji = 2, jpim1
(  153)             spgu(ji,jj) = - grav * ( sshn(ji+1,jj  ) - sshn(ji,jj) ) / e1u(ji,jj)
(  154)             spgv(ji,jj) = - grav * ( sshn(ji  ,jj+1) - sshn(ji,jj) ) / e2v(ji,jj)
(  155)          END DO 
(  156) 
(  157)          !  Add the surface pressure trend to the general trend
(  158)          DO jk = 1, jpkm1
(  159)             DO ji = 2, jpim1
(  160)                ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
(  161)                va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
(  162)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 4

(  163)          END DO
(  164) 
(  165)          ! Evaluate the masked next velocity (effect of the additional force not included)
(  166)          DO jk = 1, jpkm1
(  167)             DO ji = 2, jpim1
(  168)                ua(ji,jj,jk) = ( ub(ji,jj,jk) + z2dt * ua(ji,jj,jk) ) * umask(ji,jj,jk)
(  169)                va(ji,jj,jk) = ( vb(ji,jj,jk) + z2dt * va(ji,jj,jk) ) * vmask(ji,jj,jk)
(  170)             END DO
(  171)          END DO
(  172) 
(  173)          !                                             ! ===============
(  174)       END DO                                           !   End of slab
(  175)       !                                                ! ===============
(  176) 
(  177) #if defined key_obc
(  178)          ! Update velocities on each open boundary with the radiation algorithm
(  179)          CALL obc_dyn( kt )
(  180)          ! Correction of the barotropic componant velocity to control the volume of the system
(  181)          CALL obc_vol( kt )
(  182) #endif
(  183) #if defined key_agrif
(  184)       ! Update velocities on each coarse/fine interfaces
(  185) 
(  186)       CALL Agrif_dyn( kt )
(  187) #endif
(  188) 
(  189) !!DB: deleted ORCA
(  190) !#if defined key_orca_r2
(  191) !      IF( n_cla == 1 )   CALL dyn_spg_cla( kt )      ! Cross Land Advection (Update (ua,va))
(  192) !#endif
(  193) 
(  194)       !                                                ! ===============
(  195)       DO jj = 2, jpjm1                                 !  Vertical slab
(  196)          !                                             ! ===============
(  197) 
(  198)          ! 2. compute the next vertically averaged velocity
(  199)          ! ------------------------------------------------
(  200)          !     (effect of the additional force not included)
(  201)          ! initialize to zero
(  202)          DO ji = 2, jpim1
(  203)             spgu(ji,jj) = 0.e0
(  204)             spgv(ji,jj) = 0.e0
(  205)          END DO
(  206) 
(  207)          ! vertical sum
(  208)          DO jk = 1, jpk
(  209)             DO ji = 2, jpim1
(  210)                spgu(ji,jj) = spgu(ji,jj) + fse3u(ji,jj,jk) * ua(ji,jj,jk)
(  211)                spgv(ji,jj) = spgv(ji,jj) + fse3v(ji,jj,jk) * va(ji,jj,jk)
(  212)             END DO
(  213)          END DO
(  214) 
(  215)          ! transport: multiplied by the horizontal scale factor
(  216)          DO ji = 2, jpim1
(  217)             spgu(ji,jj) = spgu(ji,jj) * e2u(ji,jj)
(  218)             spgv(ji,jj) = spgv(ji,jj) * e1v(ji,jj)
(  219)          END DO
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 5

(  221)          !                                             ! ===============
(  222)       END DO                                           !   End of slab
(  223)       !                                                ! ===============
(  224) 
(  225)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  226)       
(  227)       ! Boundary conditions on (spgu,spgv)
(  228)       CALL lbc_lnk( spgu, 'U', -1. )
(  229)       CALL lbc_lnk( spgv, 'V', -1. )
(  230) 
(  231)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  232) 
(  233)       ! 3. Right hand side of the elliptic equation and first guess
(  234)       ! -----------------------------------------------------------
(  235)       DO jj = 2, jpjm1
(  236)          DO ji = 2, jpim1
(  237)             ! Divergence of the after vertically averaged velocity
(  238)             zgcb =  spgu(ji,jj) - spgu(ji-1,jj)   &
(  239)                   + spgv(ji,jj) - spgv(ji,jj-1)
(  240)             gcb(ji,jj) = gcdprc(ji,jj) * zgcb
(  241)             ! First guess of the after barotropic transport divergence
(  242)             zbtd = gcx(ji,jj)
(  243)             gcx (ji,jj) = 2. * zbtd   - gcxb(ji,jj)
(  244)             gcxb(ji,jj) =      zbtd
(  245)          END DO
(  246)       END DO
(  247)       ! applied the lateral boundary conditions
(  248)       IF( nsolv == 4)   CALL lbc_lnk_e( gcb, c_solver_pt, 1. )  
(  249) 
(  250) #if defined key_agrif
(  251)       IF( .NOT. AGRIF_ROOT() ) THEN
(  252)          ! add contribution of gradient of after barotropic transport divergence 
(  253)         IF( (nbondi == -1) .OR. (nbondi == 2) ) gcb(3,:) = gcb(3,:) &
(  254)            &            -znugdt * z2dt*laplacu(2,:)*gcdprc(3,:)*hu(2,:)*e2u(2,:)
(  255)         IF( (nbondi ==  1) .OR. (nbondi == 2) )  gcb(nlci-2,:) = gcb(nlci-2,:) &
(  256)            &           +znugdt * z2dt*laplacu(nlci-2,:)*gcdprc(nlci-2,:)*hu(nlci-2,:)*e2u(nlci-2,:)
(  257)         IF( (nbondj == -1) .OR. (nbondj == 2) ) gcb(:,3) = gcb(:,3) &
(  258)            &           -znugdt * z2dt*laplacv(:,2)*gcdprc(:,3)*hv(:,2)*e1v(:,2)
(  259)         IF( (nbondj ==  1) .OR. (nbondj == 2) )  gcb(:,nlcj-2) = gcb(:,nlcj-2) &
(  260)            &           +znugdt * z2dt*laplacv(:,nlcj-2)*gcdprc(:,nlcj-2)*hv(:,nlcj-2)*e1v(:,nlcj-2)
(  261)       ENDIF
(  262) #endif
(  263) 
(  264)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  265)       
(  266)       ! 4. Relative precision (computation on one processor)
(  267)       ! ---------------------
(  268)       rnorme =0.
(  269)       DO jj = 1, jpj
(  270)          DO ji = 1, jpi
(  271)             zgwgt  = gcdmat(ji,jj) * gcb(ji,jj)
(  272)             rnorme = rnorme + gcb(ji,jj) * zgwgt * bmask(ji,jj)
(  273)          END DO
(  274)       END DO
(  275)       IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  276) 
(  277)       epsr = eps * eps * rnorme
(  278)       ncut = 0






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 6

(  279)       ! if rnorme is 0, the solution is 0, the solver isn't called
(  280)       IF( rnorme == 0.e0 ) THEN
(  281)          gcx(:,:) = 0.e0
(  282)          res   = 0.e0
(  283)          niter = 0
(  284)          ncut  = 999
(  285)       ENDIF
(  286)       
(  287)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  288)       
(  289)       ! 5. Evaluate the next transport divergence
(  290)       ! -----------------------------------------
(  291)       !    Iterarive solver for the elliptic equation (except IF sol.=0)
(  292)       !    (output in gcx with boundary conditions applied)
(  293)       kindic = 0
(  294)       IF( ncut == 0 ) THEN
(  295)          IF( nsolv == 1 ) THEN         ! diagonal preconditioned conjuguate gradient
(  296)             CALL sol_pcg( kindic )
(  297)          ELSEIF( nsolv == 2 ) THEN     ! successive-over-relaxation
(  298)             CALL sol_sor( kindic )
(  299)          ELSEIF( nsolv == 3 ) THEN     ! FETI solver
(  300)             CALL sol_fet( kindic )
(  301)          ELSEIF( nsolv == 4 ) THEN     ! successive-over-relaxation with extra outer halo
(  302)             CALL sol_sor_e( kindic )
(  303)          ELSE                          ! e r r o r in nsolv namelist parameter
(  304)             IF(lwp) WRITE(numout,cform_err)
(  305)             IF(lwp) WRITE(numout,*) ' dyn_spg_flt_jki : e r r o r, nsolv = 1, 2, 3 or 4'
(  306)             IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~                not = ', nsolv
(  307)             nstop = nstop + 1
(  308)          ENDIF
(  309)       ENDIF
(  310) 
(  311)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  312) 
(  313) !CDIR PARALLEL DO
(  314) !$OMP PARALLEL DO
(  315)       !                                                ! ===============
(  316)       DO jj = 2, jpjm1                                 !  Vertical slab
(  317)          !                                             ! ===============
(  318)          
(  319)          ! 6. Transport divergence gradient multiplied by z2dt
(  320)          ! -----------------------------------------------====
(  321)          DO ji = 2, jpim1
(  322)             ! trend of Transport divergence gradient
(  323)             ztdgu = znugdt * (gcx(ji+1,jj  ) - gcx(ji,jj) ) / e1u(ji,jj)
(  324)             ztdgv = znugdt * (gcx(ji  ,jj+1) - gcx(ji,jj) ) / e2v(ji,jj)
(  325)             ! multiplied by z2dt
(  326) #if defined key_obc
(  327)             ! caution : grad D = 0 along open boundaries
(  328)             spgu(ji,jj) = z2dt * ztdgu * obcumask(ji,jj)
(  329)             spgv(ji,jj) = z2dt * ztdgv * obcvmask(ji,jj)
(  330) #else
(  331)             spgu(ji,jj) = z2dt * ztdgu
(  332)             spgv(ji,jj) = z2dt * ztdgv
(  333) #endif
(  334)          END DO
(  335)          
(  336) #if defined key_agrif      






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 7

(  337)       IF( .NOT. Agrif_Root() ) THEN
(  338)          ! caution : grad D (fine) = grad D (coarse) at coarse/fine interface
(  339)          IF( (nbondi == -1) .OR. (nbondi == 2) ) spgu(2,:) = znugdt * z2dt * laplacu(2,:) * umask(2,:,1)
(  340)          IF( (nbondi ==  1) .OR. (nbondi == 2) ) spgu(nlci-2,:) = znugdt * z2dt * laplacu(nlci-2,:) * umask(nlci-2,:,1)
(  341)          IF( (nbondj == -1) .OR. (nbondj == 2) ) spgv(:,2) = znugdt * z2dt * laplacv(:,2) * vmask(:,2,1)
(  342)          IF( (nbondj ==  1) .OR. (nbondj == 2) ) spgv(:,nlcj-2) = znugdt * z2dt * laplacv(:,nlcj-2) * vmask(:,nlcj-2,1)
(  343)       ENDIF
(  344) #endif
(  345)       
(  346)          ! 7.  Add the trends multiplied by z2dt to the after velocity
(  347)          ! -----------------------------------------------------------
(  348)          !     ( c a u t i o n : (ua,va) here are the after velocity not the
(  349)          !                       trend, the leap-frog time stepping will not
(  350)          !                       be done in dynnxt.F routine)
(  351)          DO jk = 1, jpkm1
(  352)             DO ji = 2, jpim1
(  353)                ua(ji,jj,jk) = (ua(ji,jj,jk) + spgu(ji,jj)) * umask(ji,jj,jk)
(  354)                va(ji,jj,jk) = (va(ji,jj,jk) + spgv(ji,jj)) * vmask(ji,jj,jk)
(  355)             END DO
(  356)          END DO
(  357) 
(  358)          ! 8. Sea surface elevation time stepping
(  359)          ! --------------------------------------
(  360)          ! Euler (forward) time stepping, no time filter
(  361)          IF( neuler == 0 .AND. kt == nit000 ) THEN
(  362)             DO ji = 1, jpi
(  363)                ! after free surface elevation
(  364)                zssha = sshb(ji,jj) + rdt * ( wn(ji,jj,1) - emp(ji,jj) * zraur ) * tmask(ji,jj,1)
(  365)                ! swap of arrays
(  366)                sshb(ji,jj) = sshn(ji,jj)
(  367)                sshn(ji,jj) = zssha
(  368)             END DO
(  369)          ELSE
(  370)             ! Leap-frog time stepping and time filter
(  371)             DO ji = 1, jpi
(  372)                ! after free surface elevation
(  373)                zssha = sshb(ji,jj) + z2dt * ( wn(ji,jj,1) - emp(ji,jj) * zraur ) * tmask(ji,jj,1)
(  374)                ! time filter and array swap
(  375)                sshb(ji,jj) = atfp * ( sshb(ji,jj) + zssha ) + atfp1 * sshn(ji,jj)
(  376)                sshn(ji,jj) = zssha
(  377)             END DO
(  378)          ENDIF
(  379)          !                                             ! ===============
(  380)       END DO                                           !   End of slab
(  381)       !                                                ! ===============
(  382) 
(  383)       !Boundary conditions on sshn
(  384)       CALL lbc_lnk( sshn, 'T', 1. )
(  385) 
(  386)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  387)          CALL prt_ctl( tab3d_1=ua  , clinfo1=' spg  - Ua : ', mask1=umask,   &
(  388)             &          tab3d_2=va  , clinfo2='        Va : ', mask2=vmask )
(  389)          CALL prt_ctl( tab2d_1=sshn, clinfo1=' spg  - ssh: ', mask1=tmask)
(  390)       ENDIF
(  391) 
(  392) 
(  393)    END SUBROUTINE dyn_spg_flt_jki
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 8

(  395) #else
(  396)    !!----------------------------------------------------------------------
(  397)    !!   Default case :                                         Empty module
(  398)    !!----------------------------------------------------------------------
(  399)    USE in_out_manager
(  400) CONTAINS
(  401)    SUBROUTINE dyn_spg_flt_jki( kt, kindic )      ! Empty module
(  402)       if(lwp) WRITE(numout,*) 'dyn_spg_flt_jki: You should not have seen this print! error?', kt, kindic
(  403)    END SUBROUTINE dyn_spg_flt_jki
(  404) #endif
(  405)    
(  406)    !!======================================================================
(  407) END MODULE dynspg_flt_jki



















































PGF90 (Version     10.5)          09/23/2010  09:49:42      page 9

(    1) # 1 "dynspg_flt_jki.F90"
(    1) MODULE dynspg_flt_jki
(    2)    !!======================================================================
(    3)    !!                  ***  MODULE  dynspg_flt_jki  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) # 396
(  396)    !!----------------------------------------------------------------------
(  397)    !!   Default case :                                         Empty module
(  398)    !!----------------------------------------------------------------------
(  399)    USE in_out_manager
(  400) CONTAINS
(  401)    SUBROUTINE dyn_spg_flt_jki( kt, kindic )      ! Empty module
(  402)       if(lwp) WRITE(numout,*) 'dyn_spg_flt_jki: You should not have seen this print! error?', kt, kindic
(  403)    END SUBROUTINE dyn_spg_flt_jki
(  404) # 405
(  405)    
(  406)    !!======================================================================
(  407) END MODULE dynspg_flt_jki
