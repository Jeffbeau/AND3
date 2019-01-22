


PGF90 (Version     10.5)          09/23/2010  09:49:43      page 1

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

Filename: dynspg_flt.F90

(    1) MODULE dynspg_flt
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_flt  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) #if ( defined key_dynspg_flt && ! defined key_autotasking )   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_dynspg_flt'                              filtered free surface
(    9)    !!   NOT 'key_autotasking'                      k-j-i loop (vector opt.)
(   10)    !!----------------------------------------------------------------------
(   11)    !!   dyn_spg_flt  : update the momentum trend with the surface pressure
(   12)    !!                  gradient in the filtered free surface case with
(   13)    !!                  vector optimization
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers 
(   17)    USE dom_oce         ! ocean space and time domain 
(   18)    USE zdf_oce         ! ocean vertical physics
(   19)    USE phycst          ! physical constants
(   20)    USE ocesbc          ! ocean surface boundary condition
(   21)    USE flxrnf          ! ocean runoffs
(   22)    USE sol_oce         ! ocean elliptic solver
(   23)    USE solpcg          ! preconditionned conjugate gradient solver
(   24)    USE solsor          ! Successive Over-relaxation solver
(   25)    USE solfet          ! FETI solver
(   26)    USE solsor_e        ! Successive Over-relaxation solver with MPP optimization
(   27)    USE obc_oce         ! Lateral open boundary condition
(   28)    USE obcdyn          ! ocean open boundary condition (obc_dyn routines)
(   29)    USE obcvol          ! ocean open boundary condition (obc_vol routines)
(   30)    USE lib_mpp         ! distributed memory computing library
(   31)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   32)    USE cla_dynspg      ! cross land advection
(   33)    USE prtctl          ! Print control
(   34)    USE in_out_manager  ! I/O manager
(   35)    USE solmat          ! matrix construction for elliptic solvers
(   36)    USE agrif_opa_interp
(   37) 
(   38)    IMPLICIT NONE
(   39)    PRIVATE
(   40) 
(   41)    !! * Accessibility
(   42)    PUBLIC dyn_spg_flt  ! routine called by step.F90
(   43) 
(   44)    !! * Substitutions
(   45) #  include "domzgr_substitute.h90"
(   46) #  include "vectopt_loop_substitute.h90"






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 2

(   47)    !!----------------------------------------------------------------------
(   48)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   49)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_flt.F90,v 1.4 2006/03/20 17:25:49 opalod Exp $ 
(   50)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   51)    !!----------------------------------------------------------------------
(   52) 
(   53) CONTAINS
(   54) 
(   55)    SUBROUTINE dyn_spg_flt( kt, kindic )
(   56)       !!----------------------------------------------------------------------
(   57)       !!                  ***  routine dyn_spg_flt  ***
(   58)       !!
(   59)       !! ** Purpose :   Compute the now trend due to the surface pressure 
(   60)       !!      gradient in case of filtered free surface formulation  and add
(   61)       !!      it to the general trend of momentum equation.
(   62)       !!
(   63)       !! ** Method  :   Filtered free surface formulation. The surface
(   64)       !!      pressure gradient is given by:
(   65)       !!         spgu = 1/rau0 d/dx(ps) =  1/e1u di( sshn + rnu btda )
(   66)       !!         spgv = 1/rau0 d/dy(ps) =  1/e2v dj( sshn + rnu btda )
(   67)       !!      where sshn is the free surface elevation and btda is the after
(   68)       !!      of the free surface elevation
(   69)       !!       -1- compute the after sea surface elevation from the kinematic
(   70)       !!      surface boundary condition:
(   71)       !!              zssha = sshb + 2 rdt ( wn - emp )
(   72)       !!           Time filter applied on now sea surface elevation to avoid 
(   73)       !!      the divergence of two consecutive time-steps and swap of free
(   74)       !!      surface arrays to start the next time step:
(   75)       !!              sshb = sshn + atfp * [ sshb + zssha - 2 sshn ]
(   76)       !!              sshn = zssha
(   77)       !!       -2- evaluate the surface presure trend (including the addi-
(   78)       !!      tional force) in three steps:
(   79)       !!        a- compute the right hand side of the elliptic equation:
(   80)       !!            gcb = 1/(e1t e2t) [ di(e2u spgu) + dj(e1v spgv) ]
(   81)       !!         where (spgu,spgv) are given by:
(   82)       !!            spgu = vertical sum[ e3u (ub+ 2 rdt ua ) ]
(   83)       !!                 - grav 2 rdt hu /e1u di[sshn + emp]
(   84)       !!            spgv = vertical sum[ e3v (vb+ 2 rdt va) ]
(   85)       !!                 - grav 2 rdt hv /e2v dj[sshn + emp]
(   86)       !!         and define the first guess from previous computation :
(   87)       !!            zbtd = btda
(   88)       !!            btda = 2 zbtd - btdb
(   89)       !!            btdb = zbtd
(   90)       !!        b- compute the relative accuracy to be reached by the
(   91)       !!         iterative solver
(   92)       !!        c- apply the solver by a call to sol... routine
(   93)       !!       -3- compute and add the free surface pressure gradient inclu-
(   94)       !!      ding the additional force used to stabilize the equation.
(   95)       !!
(   96)       !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
(   97)       !!
(   98)       !! References :
(   99)       !!      Roullet and Madec 1999, JGR.
(  100)       !!
(  101)       !! History :
(  102)       !!        !  98-05 (G. Roullet)  Original code
(  103)       !!        !  98-10 (G. Madec, M. Imbard)  release 8.2
(  104)       !!   8.5  !  02-08 (G. Madec)  F90: Free form and module






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 3

(  105)       !!        !  02-11 (C. Talandier, A-M Treguier) Open boundaries
(  106)       !!   9.0  !  04-08 (C. Talandier) New trends organization
(  107)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  108)       !!---------------------------------------------------------------------
(  109)       !! * Arguments
(  110)       INTEGER, INTENT( in )  ::   kt         ! ocean time-step index
(  111)       INTEGER, INTENT( out ) ::   kindic     ! solver convergence flag
(  112)                                              ! if the solver doesn t converge
(  113)                                              ! the flag is < 0
(  114)       !! * Local declarations
(  115)       INTEGER  ::   ji, jj, jk               ! dummy loop indices
(  116)       REAL(wp) ::                         & 
(  117)          z2dt, z2dtg, zraur, znugdt,      &  ! temporary scalars
(  118)          znurau, zssha, zgcb, zbtd,       &  !   "          "
(  119)          ztdgu, ztdgv                        !   "          "
(  120)       !!----------------------------------------------------------------------
(  121) 
(  122)       IF( kt == nit000 ) THEN
(  123)          IF(lwp) WRITE(numout,*)
(  124)          IF(lwp) WRITE(numout,*) 'dyn_spg_flt : surface pressure gradient trend'
(  125)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   (free surface constant volume case)'
(  126)        
(  127)          ! set to zero free surface specific arrays
(  128)          spgu(:,:) = 0.e0                     ! surface pressure gradient (i-direction)
(  129)          spgv(:,:) = 0.e0                     ! surface pressure gradient (j-direction)
(  130)       ENDIF
(  131) 
(  132)       ! Local constant initialization
(  133)       z2dt = 2. * rdt                                       ! time step: leap-frog
(  134)       IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt       ! time step: Euler if restart from rest
(  135)       IF( neuler == 0 .AND. kt == nit000+1 ) CALL sol_mat(kt)
(  136)       z2dtg  = grav * z2dt
(  137)       zraur  = 1. / rauw
(  138)       znugdt =  rnu * grav * z2dt
(  139)       znurau =  znugdt * zraur
(  140) 
(  141)       ! Surface pressure gradient (now)
(  142)       DO jj = 2, jpjm1
(  143)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  144)             spgu(ji,jj) = - grav * ( sshn(ji+1,jj) - sshn(ji,jj) ) / e1u(ji,jj)
(  145)             spgv(ji,jj) = - grav * ( sshn(ji,jj+1) - sshn(ji,jj) ) / e2v(ji,jj)
(  146)          END DO 
(  147)       END DO 
(  148) 
(  149)       ! Add the surface pressure trend to the general trend
(  150)       DO jk = 1, jpkm1
(  151)          DO jj = 2, jpjm1
(  152)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  153)                ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
(  154)                va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
(  155)             END DO
(  156)          END DO
(  157)       END DO
(  158) 
(  159)       ! Evaluate the masked next velocity (effect of the additional force not included)
(  160)       DO jk = 1, jpkm1
(  161)          DO jj = 2, jpjm1
(  162)             DO ji = fs_2, fs_jpim1   ! vector opt.






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 4

(  163)                ua(ji,jj,jk) = ( ub(ji,jj,jk) + z2dt * ua(ji,jj,jk) ) * umask(ji,jj,jk)
(  164)                va(ji,jj,jk) = ( vb(ji,jj,jk) + z2dt * va(ji,jj,jk) ) * vmask(ji,jj,jk)
(  165)             END DO
(  166)          END DO
(  167)       END DO
(  168) 
(  169) #if defined key_obc
(  170)       ! Update velocities on each open boundary with the radiation algorithm
(  171)       CALL obc_dyn( kt )
(  172)       ! Correction of the barotropic componant velocity to control the volume of the system
(  173)       CALL obc_vol( kt )
(  174) #endif
(  175) #if defined key_agrif
(  176)       ! Update velocities on each coarse/fine interfaces
(  177) 
(  178)       CALL Agrif_dyn( kt )
(  179) #endif
(  180) 
(  181) !!DB: deleted ORCA
(  182) !#if defined key_orca_r2
(  183) !      IF( n_cla == 1 )   CALL dyn_spg_cla( kt )      ! Cross Land Advection (update (ua,va))
(  184) !#endif
(  185) 
(  186)       ! compute the next vertically averaged velocity (effect of the additional force not included)
(  187)       ! ---------------------------------------------
(  188)       DO jj = 2, jpjm1
(  189)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  190)             spgu(ji,jj) = 0.e0
(  191)             spgv(ji,jj) = 0.e0
(  192)          END DO
(  193)       END DO
(  194) 
(  195)       ! vertical sum
(  196) !CDIR NOLOOPCHG
(  197)       IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
(  198)          DO jk = 1, jpkm1
(  199)             DO ji = 1, jpij
(  200)                spgu(ji,1) = spgu(ji,1) + fse3u(ji,1,jk) * ua(ji,1,jk)
(  201)                spgv(ji,1) = spgv(ji,1) + fse3v(ji,1,jk) * va(ji,1,jk)
(  202)             END DO
(  203)          END DO
(  204)       ELSE                        ! No  vector opt.
(  205)          DO jk = 1, jpkm1
(  206)             DO jj = 2, jpjm1
(  207)                DO ji = 2, jpim1
(  208)                   spgu(ji,jj) = spgu(ji,jj) + fse3u(ji,jj,jk) * ua(ji,jj,jk)
(  209)                   spgv(ji,jj) = spgv(ji,jj) + fse3v(ji,jj,jk) * va(ji,jj,jk)
(  210)                END DO
(  211)             END DO
(  212)          END DO
(  213)       ENDIF
(  214) 
(  215)       ! transport: multiplied by the horizontal scale factor
(  216)       DO jj = 2, jpjm1
(  217)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  218)             spgu(ji,jj) = spgu(ji,jj) * e2u(ji,jj)
(  219)             spgv(ji,jj) = spgv(ji,jj) * e1v(ji,jj)
(  220)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 5

(  221)       END DO
(  222) 
(  223)       ! Boundary conditions on (spgu,spgv)
(  224)       CALL lbc_lnk( spgu, 'U', -1. )
(  225)       CALL lbc_lnk( spgv, 'V', -1. )
(  226) 
(  227)       ! Right hand side of the elliptic equation and first guess
(  228)       ! -----------------------------------------------------------
(  229)       DO jj = 2, jpjm1
(  230)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  231)             ! Divergence of the after vertically averaged velocity
(  232)             zgcb =  spgu(ji,jj) - spgu(ji-1,jj)   &
(  233)                   + spgv(ji,jj) - spgv(ji,jj-1)
(  234)             gcb(ji,jj) = gcdprc(ji,jj) * zgcb
(  235)             ! First guess of the after barotropic transport divergence
(  236)             zbtd = gcx(ji,jj)
(  237)             gcx (ji,jj) = 2. * zbtd   - gcxb(ji,jj)
(  238)             gcxb(ji,jj) =      zbtd
(  239)          END DO
(  240)       END DO
(  241)       ! applied the lateral boundary conditions
(  242)       IF( nsolv == 4 )   CALL lbc_lnk_e( gcb, c_solver_pt, 1. )   
(  243) 
(  244) #if defined key_agrif
(  245)       IF( .NOT. AGRIF_ROOT() ) THEN
(  246)          ! add contribution of gradient of after barotropic transport divergence 
(  247)          IF( (nbondi == -1) .OR. (nbondi == 2) ) gcb(3,:) = gcb(3,:) &
(  248)             &            -znugdt * z2dt*laplacu(2,:)*gcdprc(3,:)*hu(2,:)*e2u(2,:)
(  249)          IF( (nbondi ==  1) .OR. (nbondi == 2) )  gcb(nlci-2,:) = gcb(nlci-2,:) &
(  250)             &           +znugdt * z2dt*laplacu(nlci-2,:)*gcdprc(nlci-2,:)*hu(nlci-2,:)*e2u(nlci-2,:)
(  251)          IF( (nbondj == -1) .OR. (nbondj == 2) ) gcb(:,3) = gcb(:,3) &
(  252)             &           -znugdt * z2dt*laplacv(:,2)*gcdprc(:,3)*hv(:,2)*e1v(:,2)
(  253)          IF( (nbondj ==  1) .OR. (nbondj == 2) )  gcb(:,nlcj-2) = gcb(:,nlcj-2) &
(  254)             &           +znugdt * z2dt*laplacv(:,nlcj-2)*gcdprc(:,nlcj-2)*hv(:,nlcj-2)*e1v(:,nlcj-2)
(  255)       ENDIF
(  256) #endif
(  257) 
(  258) 
(  259)       ! Relative precision (computation on one processor)
(  260)       ! ------------------
(  261)       rnorme =0.
(  262)       rnorme = SUM( gcb(1:jpi,1:jpj) * gcdmat(1:jpi,1:jpj) * gcb(1:jpi,1:jpj) * bmask(:,:) )
(  263)       IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  264) 
(  265)       epsr = eps * eps * rnorme
(  266)       ncut = 0
(  267)       ! if rnorme is 0, the solution is 0, the solver isn't called
(  268)       IF( rnorme == 0.e0 ) THEN
(  269)          gcx(:,:) = 0.e0
(  270)          res   = 0.e0
(  271)          niter = 0
(  272)          ncut  = 999
(  273)       ENDIF
(  274) 
(  275)       ! Evaluate the next transport divergence
(  276)       ! --------------------------------------
(  277)       !    Iterarive solver for the elliptic equation (except IF sol.=0)
(  278)       !    (output in gcx with boundary conditions applied)






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 6

(  279)       kindic = 0
(  280)       IF( ncut == 0 ) THEN
(  281)          IF( nsolv == 1 ) THEN         ! diagonal preconditioned conjuguate gradient
(  282)             CALL sol_pcg( kindic )
(  283)          ELSEIF( nsolv == 2 ) THEN     ! successive-over-relaxation
(  284)             CALL sol_sor( kindic )
(  285)          ELSEIF( nsolv == 3 ) THEN     ! FETI solver
(  286)             CALL sol_fet( kindic )
(  287)          ELSEIF( nsolv == 4 ) THEN     ! successive-over-relaxation with extra outer halo
(  288)             CALL sol_sor_e( kindic )
(  289)          ELSE                          ! e r r o r in nsolv namelist parameter
(  290)             IF(lwp) WRITE(numout,cform_err)
(  291)             IF(lwp) WRITE(numout,*) ' dyn_spg_flt : e r r o r, nsolv = 1, 2 ,3 or 4'
(  292)             IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~                not = ', nsolv
(  293)             nstop = nstop + 1
(  294)          ENDIF
(  295)       ENDIF
(  296) 
(  297)       ! Transport divergence gradient multiplied by z2dt
(  298)       ! --------------------------------------------====
(  299)       DO jj = 2, jpjm1
(  300)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  301)             ! trend of Transport divergence gradient
(  302)             ztdgu = znugdt * (gcx(ji+1,jj  ) - gcx(ji,jj) ) / e1u(ji,jj)
(  303)             ztdgv = znugdt * (gcx(ji  ,jj+1) - gcx(ji,jj) ) / e2v(ji,jj)
(  304)             ! multiplied by z2dt
(  305) #if defined key_obc
(  306)             ! caution : grad D = 0 along open boundaries
(  307)             spgu(ji,jj) = z2dt * ztdgu * obcumask(ji,jj)
(  308)             spgv(ji,jj) = z2dt * ztdgv * obcvmask(ji,jj)
(  309) #else
(  310)             spgu(ji,jj) = z2dt * ztdgu
(  311)             spgv(ji,jj) = z2dt * ztdgv
(  312) #endif
(  313)          END DO
(  314)       END DO
(  315) 
(  316) #if defined key_agrif      
(  317)       IF( .NOT. Agrif_Root() ) THEN
(  318)          ! caution : grad D (fine) = grad D (coarse) at coarse/fine interface
(  319)          IF( (nbondi == -1) .OR. (nbondi == 2) ) spgu(2,:) = znugdt * z2dt * laplacu(2,:) * umask(2,:,1)
(  320)          IF( (nbondi ==  1) .OR. (nbondi == 2) ) spgu(nlci-2,:) = znugdt * z2dt * laplacu(nlci-2,:) * umask(nlci-2,:,1)
(  321)          IF( (nbondj == -1) .OR. (nbondj == 2) ) spgv(:,2) = znugdt * z2dt * laplacv(:,2) * vmask(:,2,1)
(  322)          IF( (nbondj ==  1) .OR. (nbondj == 2) ) spgv(:,nlcj-2) = znugdt * z2dt * laplacv(:,nlcj-2) * vmask(:,nlcj-2,1)
(  323)       ENDIF
(  324) #endif      
(  325)       ! 7.  Add the trends multiplied by z2dt to the after velocity
(  326)       ! -----------------------------------------------------------
(  327)       !     ( c a u t i o n : (ua,va) here are the after velocity not the
(  328)       !                       trend, the leap-frog time stepping will not
(  329)       !                       be done in dynnxt.F routine)
(  330)       DO jk = 1, jpkm1
(  331)          DO jj = 2, jpjm1
(  332)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  333)                ua(ji,jj,jk) = (ua(ji,jj,jk) + spgu(ji,jj)) * umask(ji,jj,jk)
(  334)                va(ji,jj,jk) = (va(ji,jj,jk) + spgv(ji,jj)) * vmask(ji,jj,jk)
(  335)             END DO
(  336)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 7

(  337)       END DO
(  338) 
(  339) 
(  340)       ! Sea surface elevation time stepping
(  341)       ! -----------------------------------
(  342)       IF( neuler == 0 .AND. kt == nit000 ) THEN      ! Euler (forward) time stepping, no time filter
(  343)          DO jj = 1, jpj
(  344)             DO ji = 1, jpi
(  345)                ! after free surface elevation
(  346)                zssha = sshb(ji,jj) + rdt * ( wn(ji,jj,1) - emp(ji,jj) * zraur ) * tmask(ji,jj,1)
(  347)                ! swap of arrays
(  348)                sshb(ji,jj) = sshn(ji,jj)
(  349)                sshn(ji,jj) = zssha
(  350)             END DO
(  351)          END DO
(  352)       ELSE                                           ! Leap-frog time stepping and time filter
(  353)          DO jj = 1, jpj
(  354)             DO ji = 1, jpi
(  355)                ! after free surface elevation
(  356)                zssha = sshb(ji,jj) + z2dt * ( wn(ji,jj,1) - emp(ji,jj) * zraur ) * tmask(ji,jj,1)
(  357)                ! time filter and array swap
(  358)                sshb(ji,jj) = atfp * ( sshb(ji,jj) + zssha ) + atfp1 * sshn(ji,jj)
(  359)                sshn(ji,jj) = zssha
(  360)             END DO
(  361)          END DO
(  362)       ENDIF
(  363) 
(  364)       !                       ! print sum trends (used for debugging)
(  365)       IF(ln_ctl)   CALL prt_ctl( tab2d_1=sshn, clinfo1=' spg  - ssh: ', mask1=tmask )
(  366) 
(  367)    END SUBROUTINE dyn_spg_flt
(  368) 
(  369) #else
(  370)    !!----------------------------------------------------------------------
(  371)    !!   Default case :   Empty module   No standart free surface cst volume
(  372)    !!----------------------------------------------------------------------
(  373)    USE in_out_manager
(  374) CONTAINS
(  375)    SUBROUTINE dyn_spg_flt( kt, kindic )       ! Empty routine
(  376)       if(lwp) WRITE(numout,*) 'dyn_spg_flt: You should not have seen this print! error?', kt, kindic
(  377)    END SUBROUTINE dyn_spg_flt
(  378) #endif
(  379)    
(  380)    !!======================================================================
(  381) END MODULE dynspg_flt



















PGF90 (Version     10.5)          09/23/2010  09:49:43      page 8

(    1) # 1 "dynspg_flt.F90"
(    1) MODULE dynspg_flt
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_flt  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) # 370
(  370)    !!----------------------------------------------------------------------
(  371)    !!   Default case :   Empty module   No standart free surface cst volume
(  372)    !!----------------------------------------------------------------------
(  373)    USE in_out_manager
(  374) CONTAINS
(  375)    SUBROUTINE dyn_spg_flt( kt, kindic )       ! Empty routine
(  376)       if(lwp) WRITE(numout,*) 'dyn_spg_flt: You should not have seen this print! error?', kt, kindic
(  377)    END SUBROUTINE dyn_spg_flt
(  378) # 379
(  379)    
(  380)    !!======================================================================
(  381) END MODULE dynspg_flt
