


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

Filename: dynspg_rl.F90

(    1) MODULE dynspg_rl
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  dynspg_rl  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) #if   defined key_dynspg_rl   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_dynspg_rl'                                           rigid lid
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dyn_spg_rl   : update the momentum trend with the surface pressure
(   11)    !!                  for the rigid-lid case.
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE phycst          ! physical constants
(   17)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   18)    USE ldfdyn_oce      ! ocean dynamics: lateral physics
(   19)    USE zdf_oce         ! ocean vertical physics
(   20)    USE sol_oce         ! ocean elliptic solver
(   21)    USE solpcg          ! preconditionned conjugate gradient solver
(   22)    USE solsor          ! Successive Over-relaxation solver
(   23)    USE solfet          ! FETI solver
(   24)    USE solsor_e        ! Successive Over-relaxation solver with MPP optimization
(   25)    USE solisl          ! ???
(   26)    USE obc_oce         ! Lateral open boundary condition
(   27)    USE lib_mpp         ! distributed memory computing library
(   28)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   29)    USE in_out_manager  ! I/O manager
(   30) 
(   31)    IMPLICIT NONE
(   32)    PRIVATE
(   33) 
(   34)    !! * Accessibility
(   35)    PUBLIC dyn_spg_rl   ! called by step.F90
(   36) 
(   37)    !! * Substitutions
(   38) #  include "domzgr_substitute.h90"
(   39) #  include "vectopt_loop_substitute.h90"
(   40) #  include "obc_vectopt_loop_substitute.h90"
(   41)    !!----------------------------------------------------------------------
(   42)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   43)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_rl.F90,v 1.8 2005/12/21 10:46:38 opalod Exp $ 
(   44)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   45)    !!----------------------------------------------------------------------
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 2

(   47) CONTAINS
(   48) 
(   49)    SUBROUTINE dyn_spg_rl( kt, kindic )
(   50)       !!----------------------------------------------------------------------
(   51)       !!                  ***  routine dyn_spg_rl  ***
(   52)       !!
(   53)       !! ** Purpose :   Compute the now trend due to the surface pressure
(   54)       !!      gradient for the rigid-lid case,  add it to the general trend of 
(   55)       !!      momentum equation.
(   56)       !!
(   57)       !! ** Method  :   Rigid-lid appromimation: the surface pressure gradient 
(   58)       !!      is given by:
(   59)       !!         spgu = 1/rau0 d/dx(ps) = Mu + 1/(hu e2u) dj-1(bsfd)
(   60)       !!         spgv = 1/rau0 d/dy(ps) = Mv - 1/(hv e1v) di-1(bsfd)
(   61)       !!      where (Mu,Mv) is the vertically averaged momentum trend (i.e.
(   62)       !!      the vertical ponderated sum of the general momentum trend),
(   63)       !!      and bsfd is the barotropic streamfunction trend.
(   64)       !!      The trend is computed as follows:
(   65)       !!         -1- compute the vertically averaged momentum trend (Mu,Mv)
(   66)       !!         -2- compute the barotropic streamfunction trend by solving an
(   67)       !!      ellipic equation using a diagonal preconditioned conjugate
(   68)       !!      gradient or a successive-over-relaxation method (depending
(   69)       !!      on nsolv, a namelist parameter).
(   70)       !!         -3- add to bsfd the island trends if lk_isl=T
(   71)       !!         -4- compute the after streamfunction is for further diagnos-
(   72)       !!      tics using a leap-frog scheme.
(   73)       !!         -5- add the momentum trend associated with the surface pres-
(   74)       !!      sure gradient to the general trend.
(   75)       !!
(   76)       !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
(   77)       !!
(   78)       !! References :
(   79)       !!      Madec et al. 1988, ocean modelling, issue 78, 1-6.
(   80)       !!
(   81)       !! History :
(   82)       !!        !  96-05  (G. Madec, M. Imbard, M. Guyon)  rewitting in 1
(   83)       !!                  routine, without pointers, and with the same matrix
(   84)       !!                  for sor and pcg, mpp exchange, and symmetric conditions
(   85)       !!        !  96-07  (A. Weaver)  Euler forward step
(   86)       !!        !  96-11  (A. Weaver)  correction to preconditioning:
(   87)       !!        !  98-02  (M. Guyon)  FETI method
(   88)       !!        !  98-05  (G. Roullet)  free surface
(   89)       !!        !  97-09  (J.-M. Molines)  Open boundaries
(   90)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   91)       !!        !  02-11  (C. Talandier, A-M Treguier) Open boundaries
(   92)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   93)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   94)       !!---------------------------------------------------------------------
(   95)       !! * Arguments
(   96)       INTEGER, INTENT( in  ) ::   kt       ! ocean time-step index
(   97)       INTEGER, INTENT( out ) ::   kindic   ! solver flag, take a negative value
(   98)       !                                    ! when the solver doesnot converge
(   99)       !! * Local declarations
(  100)       INTEGER ::   ji, jj, jk    ! dummy loop indices
(  101)       REAL(wp) ::  zbsfa, zgcx, z2dt
(  102) # if defined key_obc
(  103)       INTEGER ::   ip, ii, ij
(  104)       INTEGER ::   iii, ijj, jip, jnic






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 3

(  105)       INTEGER ::   it, itm, itm2, ib, ibm, ibm2
(  106)       REAL(wp) ::   z2dtr
(  107) # endif
(  108)       !!----------------------------------------------------------------------
(  109)       
(  110)       IF( kt == nit000 ) THEN
(  111)          IF(lwp) WRITE(numout,*)
(  112)          IF(lwp) WRITE(numout,*) 'dyn_spg_rl : surface pressure gradient trend'
(  113)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
(  114) 
(  115)          ! set to zero rigid-lid specific arrays
(  116)          spgu(:,:) = 0.e0      ! surface pressure gradient (i-direction) 
(  117)          spgv(:,:) = 0.e0      ! surface pressure gradient (j-direction)
(  118)       ENDIF
(  119) 
(  120)       ! 0. Initializations:
(  121)       ! -------------------
(  122) # if defined key_obc
(  123)       ! space index on boundary arrays
(  124)       ib = 1
(  125)       ibm = 2
(  126)       ibm2 = 3
(  127)       ! time index on boundary arrays
(  128)       it = 1
(  129)       itm = 2
(  130)       itm2 = 3
(  131) # endif
(  132) 
(  133)       !                                                ! ===============
(  134)       DO jj = 2, jpjm1                                 !  Vertical slab
(  135)          !                                             ! ===============
(  136) 
(  137)          ! 1. Vertically averaged momentum trend
(  138)          ! -------------------------------------
(  139)          ! initialization to zero
(  140)          spgu(:,jj) = 0.
(  141)          spgv(:,jj) = 0.
(  142)          ! vertical sum
(  143)          DO jk = 1, jpk
(  144)             DO ji = 2, jpim1
(  145)                spgu(ji,jj) = spgu(ji,jj) + ua(ji,jj,jk) * fse3u(ji,jj,jk) * umask(ji,jj,jk)
(  146)                spgv(ji,jj) = spgv(ji,jj) + va(ji,jj,jk) * fse3v(ji,jj,jk) * vmask(ji,jj,jk)
(  147)             END DO 
(  148)          END DO 
(  149)          ! divide by the depth
(  150)          spgu(:,jj) = spgu(:,jj) * hur(:,jj)
(  151)          spgv(:,jj) = spgv(:,jj) * hvr(:,jj)
(  152) 
(  153)          !                                             ! ===============
(  154)       END DO                                           !   End of slab
(  155)       !                                                ! ===============
(  156) 
(  157)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  158) 
(  159)       ! Boundary conditions on (spgu,spgv)
(  160)       CALL lbc_lnk( spgu, 'U', -1. )
(  161)       CALL lbc_lnk( spgv, 'V', -1. )
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 4

(  163)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  164) 
(  165)       ! 2. Barotropic streamfunction trend (bsfd)
(  166)       ! ----------------------------------
(  167) 
(  168)       ! Curl of the vertically averaged velocity
(  169)       DO jj = 2, jpjm1
(  170)          DO ji = 2, jpim1
(  171)             gcb(ji,jj) = -gcdprc(ji,jj)   &
(  172)                        *(  ( e2v(ji+1,jj  )*spgv(ji+1,jj  ) - e2v(ji,jj)*spgv(ji,jj) )   &
(  173)                           -( e1u(ji  ,jj+1)*spgu(ji  ,jj+1) - e1u(ji,jj)*spgu(ji,jj) )   ) 
(  174)          END DO
(  175)       END DO
(  176) 
(  177) # if defined key_obc
(  178)       ! Open boundary contribution
(  179)       DO jj = 2, jpjm1
(  180)          DO ji = 2, jpim1
(  181)            gcb(ji,jj) = gcb(ji,jj) - gcdprc(ji,jj) * gcbob(ji,jj)
(  182)          END DO
(  183)       END DO
(  184) # else
(  185)       ! No open boundary contribution
(  186) # endif
(  187) 
(  188)       ! First guess using previous solution of the elliptic system and
(  189)       ! not bsfd since the system is solved with 0 as coastal boundary
(  190)       ! condition. Also include a swap array (gcx,gxcb)
(  191)       DO jj = 2, jpjm1
(  192)          DO ji = 2, jpim1
(  193)             zgcx        = gcx(ji,jj)
(  194)             gcx (ji,jj) = 2.*zgcx - gcxb(ji,jj)
(  195)             gcxb(ji,jj) =    zgcx
(  196)          END DO
(  197)       END DO
(  198)       ! applied the lateral boundary conditions
(  199)       IF( nsolv == 4)   CALL lbc_lnk_e( gcb, c_solver_pt, 1. )   
(  200) 
(  201)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  202) 
(  203)       ! Relative precision (computation on one processor)
(  204)       rnorme = 0.e0
(  205)       rnorme = SUM( gcb(1:nlci,1:nlcj) * gcdmat(1:nlci,1:nlcj) * gcb(1:nlci,1:nlcj) * bmask(:,:) )
(  206)       IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  207) 
(  208)       epsr = eps*eps*rnorme
(  209)       ncut = 0
(  210)       ! if rnorme is 0, the solution is 0, the solver isn't called
(  211)       IF( rnorme == 0.e0 ) THEN
(  212)          bsfd (:,:) = 0.e0
(  213)          res   = 0.e0
(  214)          niter = 0
(  215)          ncut  = 999
(  216)       ENDIF
(  217) 
(  218)       kindic = 0
(  219)       ! solve the bsf system  ===> solution in gcx array
(  220)       IF( ncut == 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 5

(  221)          SELECT CASE ( nsolv )
(  222)          CASE ( 1 )                    ! diagonal preconditioned conjuguate gradient
(  223)             CALL sol_pcg( kindic )
(  224)          CASE( 2 )                     ! successive-over-relaxation
(  225)             CALL sol_sor( kindic )
(  226)          CASE( 3 )                     ! FETI solver
(  227)             CALL sol_fet( kindic )
(  228)          CASE( 4 )                     ! successive-over-relaxation with extra outer halo
(  229)             CALL sol_sor_e( kindic )
(  230)          CASE DEFAULT                  ! e r r o r in nsolv namelist parameter
(  231)             IF(lwp) WRITE(numout,cform_err)
(  232)             IF(lwp) WRITE(numout,*) ' dyn_spg_rl : e r r o r, nsolv = 1, 2 ,3 or 4'
(  233)             IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~                not = ', nsolv
(  234)             nstop = nstop + 1
(  235)          END SELECT
(  236)       ENDIF
(  237) 
(  238) 
(  239)       ! bsf trend update
(  240)       ! ----------------
(  241) 
(  242)       bsfd(1:nlci,1:nlcj) = gcx(1:nlci,1:nlcj)
(  243) 
(  244)       
(  245)       ! update bsf trend with islands trend
(  246)       ! -----------------------------------
(  247) 
(  248)       IF( lk_isl )   CALL isl_dyn_spg         ! update bsfd
(  249) 
(  250) 
(  251) # if defined key_obc
(  252)       ! Compute bsf trend for OBC points (not masked)
(  253) 
(  254)       IF( lp_obc_east ) THEN
(  255)       ! compute bsf trend at the boundary from bsfeob, computed in obc_spg
(  256)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  257)          z2dtr = 1. / rdt
(  258)       ELSE
(  259)          z2dtr = 1. / (2. * rdt )
(  260)       ENDIF
(  261)       ! (jped,jpefm1),nieob
(  262)       DO ji = fs_nie0, fs_nie1   ! vector opt.
(  263)          DO jj = nje0m1, nje1m1
(  264)             bsfd(ji,jj) = ( bsfeob(jj) - bsfb(ji,jj) ) * z2dtr
(  265)          END DO
(  266)       END DO
(  267)       ENDIF
(  268) 
(  269)       IF( lp_obc_west ) THEN
(  270)       ! compute bsf trend at the boundary from bsfwob, computed in obc_spg
(  271)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  272)          z2dtr = 1. / rdt
(  273)       ELSE
(  274)          z2dtr = 1. / ( 2. * rdt )
(  275)       ENDIF
(  276)       ! (jpwd,jpwfm1),niwob
(  277)       DO ji = fs_niw0, fs_niw1   ! vector opt.
(  278)          DO jj = njw0m1, njw1m1






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 6

(  279)             bsfd(ji,jj) = ( bsfwob(jj) - bsfb(ji,jj) ) * z2dtr
(  280)          END DO
(  281)       END DO
(  282)       ENDIF
(  283) 
(  284)       IF( lp_obc_north ) THEN
(  285)       ! compute bsf trend at the boundary from bsfnob, computed in obc_spg
(  286)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  287)          z2dtr = 1. / rdt
(  288)       ELSE
(  289)          z2dtr = 1. / ( 2. * rdt )
(  290)       ENDIF
(  291)       ! njnob,(jpnd,jpnfm1)
(  292)       DO jj = fs_njn0, fs_njn1   ! vector opt.
(  293)          DO ji = nin0m1, nin1m1
(  294)             bsfd(ji,jj) = ( bsfnob(ji) - bsfb(ji,jj) ) * z2dtr
(  295)          END DO
(  296)       END DO
(  297)       ENDIF
(  298) 
(  299)       IF( lp_obc_south ) THEN
(  300)       ! compute bsf trend at the boundary from bsfsob, computed in obc_spg
(  301)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  302)          z2dtr = 1. / rdt
(  303)       ELSE
(  304)          z2dtr = 1. / ( 2. * rdt )
(  305)       ENDIF
(  306)       ! njsob,(jpsd,jpsfm1)
(  307)       DO jj = fs_njs0, fs_njs1   ! vector opt.
(  308)          DO ji = nis0m1, nis1m1
(  309)             bsfd(ji,jj) = ( bsfsob(ji) - bsfb(ji,jj) ) * z2dtr
(  310)          END DO
(  311)       END DO
(  312)       ENDIF
(  313) 
(  314)       ! compute bsf trend for isolated coastline points
(  315) 
(  316)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  317)          z2dtr = 1. / rdt
(  318)       ELSE
(  319)          z2dtr = 1. /( 2. * rdt )
(  320)       ENDIF
(  321) 
(  322)       IF( nbobc > 1 ) THEN
(  323)          DO jnic = 1,nbobc - 1
(  324)             ip = mnic(0,jnic)
(  325)             DO jip = 1,ip
(  326)                ii = miic(jip,0,jnic)
(  327)                ij = mjic(jip,0,jnic)
(  328)                IF( ii >= 1 + nimpp - 1 .AND. ii <= jpi + nimpp -1 .AND. &
(  329)                    ij >= 1 + njmpp - 1 .AND. ij <= jpj + njmpp -1 ) THEN
(  330)                   iii = ii - nimpp + 1
(  331)                   ijj = ij - njmpp + 1
(  332)                   bsfd(iii,ijj) = ( bsfic(jnic) - bsfb(iii,ijj) ) * z2dtr
(  333)                ENDIF
(  334)             END DO
(  335)          END DO
(  336)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 7

(  337) # endif
(  338) 
(  339)       ! 4. Barotropic stream function and array swap
(  340)       ! --------------------------------------------
(  341) 
(  342)       ! Leap-frog time scheme, time filter and array swap
(  343)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  344)          ! Euler time stepping (first time step, starting from rest)
(  345)          z2dt = rdt
(  346)          DO jj = 1, jpj
(  347)             DO ji = 1, jpi
(  348)                zbsfa       = bsfb(ji,jj) + z2dt * bsfd(ji,jj)
(  349)                bsfb(ji,jj) = bsfn(ji,jj)
(  350)                bsfn(ji,jj) = zbsfa 
(  351)             END DO
(  352)          END DO
(  353)       ELSE
(  354)          ! Leap-frog time stepping - Asselin filter
(  355)          z2dt = 2.*rdt
(  356)          DO jj = 1, jpj
(  357)             DO ji = 1, jpi
(  358)                zbsfa       = bsfb(ji,jj) + z2dt * bsfd(ji,jj)
(  359)                bsfb(ji,jj) = atfp * ( bsfb(ji,jj) + zbsfa ) + atfp1 * bsfn(ji,jj)
(  360)                bsfn(ji,jj) = zbsfa
(  361)             END DO
(  362)          END DO
(  363)       ENDIF
(  364) 
(  365) # if defined key_obc
(  366)       ! Swap of boundary arrays
(  367)       IF( lp_obc_east ) THEN
(  368)       ! (jped,jpef),nieob
(  369)       IF( kt < nit000+3 .AND. .NOT.ln_rstart ) THEN
(  370)          DO jj = nje0m1, nje1
(  371)             ! fields itm2 <== itm
(  372)             bebnd(jj,ib  ,itm2) = bebnd(jj,ib  ,itm)
(  373)             bebnd(jj,ibm ,itm2) = bebnd(jj,ibm ,itm)
(  374)             bebnd(jj,ibm2,itm2) = bebnd(jj,ibm2,itm)
(  375)             bebnd(jj,ib  ,itm ) = bebnd(jj,ib  ,it )
(  376)          END DO
(  377)       ELSE
(  378)          ! fields itm <== it  plus time filter at the boundary
(  379)          DO ji = fs_nie0, fs_nie1   ! vector opt.
(  380)             DO jj = nje0m1, nje1
(  381)                bebnd(jj,ib  ,itm2) = bebnd(jj,ib  ,itm)
(  382)                bebnd(jj,ibm ,itm2) = bebnd(jj,ibm ,itm)
(  383)                bebnd(jj,ibm2,itm2) = bebnd(jj,ibm2,itm)
(  384)                bebnd(jj,ib  ,itm ) = atfp * ( bebnd(jj,ib,itm) + bsfn(ji,jj) ) + atfp1 * bebnd(jj,ib,it)
(  385)                bebnd(jj,ibm ,itm ) = bebnd(jj,ibm ,it )
(  386)                bebnd(jj,ibm2,itm ) = bebnd(jj,ibm2,it )
(  387)             END DO
(  388)          END DO
(  389)       ENDIF
(  390)       ! fields it <== now (kt+1)
(  391)       DO ji = fs_nie0, fs_nie1   ! vector opt.
(  392)          DO jj = nje0m1, nje1
(  393)             bebnd(jj,ib  ,it  ) = bsfn (ji  ,jj)
(  394)             bebnd(jj,ibm ,it  ) = bsfn (ji-1,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 8

(  395)             bebnd(jj,ibm2,it  ) = bsfn (ji-2,jj)
(  396)          END DO
(  397)       END DO
(  398)       IF( lk_mpp )   CALL mppobc( bebnd, jpjed, jpjef, jpieob, 3*3, 2, jpj )
(  399)       ENDIF
(  400) 
(  401)       IF( lp_obc_west ) THEN
(  402)       ! (jpwd,jpwf),niwob
(  403)       IF( kt < nit000+3 .AND. .NOT.ln_rstart ) THEN
(  404)          DO jj = njw0m1, njw1
(  405)             ! fields itm2 <== itm
(  406)             bwbnd(jj,ib  ,itm2) = bwbnd(jj,ib  ,itm)
(  407)             bwbnd(jj,ibm ,itm2) = bwbnd(jj,ibm ,itm)
(  408)             bwbnd(jj,ibm2,itm2) = bwbnd(jj,ibm2,itm)
(  409)             bwbnd(jj,ib  ,itm ) = bwbnd(jj,ib  ,it )
(  410)          END DO
(  411)       ELSE
(  412)          DO ji = fs_niw0, fs_niw1   ! Vector opt.
(  413)             DO jj = njw0m1, njw1
(  414)                bwbnd(jj,ib  ,itm2) = bwbnd(jj,ib  ,itm)
(  415)                bwbnd(jj,ibm ,itm2) = bwbnd(jj,ibm ,itm)
(  416)                bwbnd(jj,ibm2,itm2) = bwbnd(jj,ibm2,itm)
(  417)                ! fields itm <== it  plus time filter at the boundary
(  418)                bwbnd(jj,ib  ,itm ) = atfp * ( bwbnd(jj,ib,itm) + bsfn(ji,jj) ) + atfp1 * bwbnd(jj,ib,it)
(  419)                bwbnd(jj,ibm ,itm ) = bwbnd(jj,ibm ,it )
(  420)                bwbnd(jj,ibm2,itm ) = bwbnd(jj,ibm2,it )
(  421)             END DO
(  422)          END DO
(  423)       ENDIF
(  424)       ! fields it <== now (kt+1)
(  425)       DO ji = fs_niw0, fs_niw1   ! Vector opt.
(  426)          DO jj = njw0m1, njw1
(  427)             bwbnd(jj,ib  ,it  ) = bsfn (ji  ,jj)
(  428)             bwbnd(jj,ibm ,it  ) = bsfn (ji+1,jj)
(  429)             bwbnd(jj,ibm2,it  ) = bsfn (ji+2,jj)
(  430)          END DO
(  431)       END DO
(  432)       IF( lk_mpp )   CALL mppobc( bwbnd, jpjwd, jpjwf, jpiwob, 3*3, 2, jpj )
(  433)       ENDIF
(  434) 
(  435)       IF( lp_obc_north ) THEN
(  436)    ! njnob,(jpnd,jpnf)
(  437)       IF( kt < nit000 + 3 .AND. .NOT.ln_rstart ) THEN
(  438)          DO ji = nin0m1, nin1
(  439)             ! fields itm2 <== itm
(  440)             bnbnd(ji,ib  ,itm2) = bnbnd(ji,ib  ,itm)
(  441)             bnbnd(ji,ibm ,itm2) = bnbnd(ji,ibm ,itm)
(  442)             bnbnd(ji,ibm2,itm2) = bnbnd(ji,ibm2,itm)
(  443)             bnbnd(ji,ib  ,itm ) = bnbnd(ji,ib  ,it )
(  444)          END DO
(  445)       ELSE
(  446)          DO jj = fs_njn0, fs_njn1   ! Vector opt.
(  447)             DO ji = nin0m1, nin1
(  448)                bnbnd(ji,ib  ,itm2) = bnbnd(ji,ib  ,itm)
(  449)                bnbnd(ji,ibm ,itm2) = bnbnd(ji,ibm ,itm)
(  450)                bnbnd(ji,ibm2,itm2) = bnbnd(ji,ibm2,itm)
(  451)                ! fields itm <== it  plus time filter at the boundary
(  452)                bnbnd(jj,ib  ,itm ) = atfp * ( bnbnd(jj,ib,itm) + bsfn(ji,jj) ) + atfp1 * bnbnd(jj,ib,it)






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 9

(  453)                bnbnd(ji,ibm ,itm ) = bnbnd(ji,ibm ,it )
(  454)                bnbnd(ji,ibm2,itm ) = bnbnd(ji,ibm2,it )
(  455)             END DO
(  456)           END DO
(  457)       ENDIF
(  458)       ! fields it <== now (kt+1)
(  459)       DO jj = fs_njn0, fs_njn1   ! Vector opt.
(  460)          DO ji = nin0m1, nin1
(  461)             bnbnd(ji,ib  ,it  ) = bsfn (ji,jj  )
(  462)             bnbnd(ji,ibm ,it  ) = bsfn (ji,jj-1)
(  463)             bnbnd(ji,ibm2,it  ) = bsfn (ji,jj-2)
(  464)          END DO
(  465)       END DO
(  466)       IF( lk_mpp )   CALL mppobc( bnbnd, jpind, jpinf, jpjnob, 3*3, 1, jpi )
(  467)       ENDIF
(  468) 
(  469)       IF( lp_obc_south ) THEN
(  470)          ! njsob,(jpsd,jpsf)
(  471)          IF( kt < nit000+3 .AND. .NOT.ln_rstart ) THEN
(  472)             DO ji = nis0m1, nis1
(  473)                ! fields itm2 <== itm
(  474)                bsbnd(ji,ib  ,itm2) = bsbnd(ji,ib  ,itm)
(  475)                bsbnd(ji,ibm ,itm2) = bsbnd(ji,ibm ,itm)
(  476)                bsbnd(ji,ibm2,itm2) = bsbnd(ji,ibm2,itm)
(  477)                bsbnd(ji,ib  ,itm ) = bsbnd(ji,ib  ,it )
(  478)             END DO
(  479)          ELSE
(  480)             DO jj = fs_njs0, fs_njs1   ! vector opt.
(  481)                DO ji = nis0m1, nis1
(  482)                   bsbnd(ji,ib  ,itm2) = bsbnd(ji,ib  ,itm)
(  483)                   bsbnd(ji,ibm ,itm2) = bsbnd(ji,ibm ,itm)
(  484)                   bsbnd(ji,ibm2,itm2) = bsbnd(ji,ibm2,itm)
(  485)                   ! fields itm <== it  plus time filter at the boundary
(  486)                   bsbnd(jj,ib  ,itm ) = atfp * ( bsbnd(jj,ib,itm) + bsfn(ji,jj) ) + atfp1 * bsbnd(jj,ib,it)
(  487)                   bsbnd(ji,ibm ,itm ) = bsbnd(ji,ibm ,it )
(  488)                   bsbnd(ji,ibm2,itm ) = bsbnd(ji,ibm2,it )
(  489)                END DO
(  490)             END DO
(  491)          ENDIF
(  492)          DO jj = fs_njs0, fs_njs1   ! vector opt.
(  493)             DO ji = nis0m1, nis1 
(  494)                ! fields it <== now (kt+1)
(  495)                bsbnd(ji,ib  ,it  ) = bsfn (ji,jj  )
(  496)                bsbnd(ji,ibm ,it  ) = bsfn (ji,jj+1)
(  497)                bsbnd(ji,ibm2,it  ) = bsfn (ji,jj+2)
(  498)             END DO
(  499)          END DO
(  500)          IF( lk_mpp )   CALL mppobc( bsbnd, jpisd, jpisf, jpjsob, 3*3, 1, jpi )
(  501)       ENDIF
(  502) # endif
(  503)       !
(  504)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  505)       
(  506)       !  add the surface pressure trend to the general trend
(  507)       ! -----------------------------------------------------
(  508)       
(  509)       DO jj = 2, jpjm1
(  510) 






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 10

(  511)          ! update the surface pressure gradient with the barotropic trend
(  512)          DO ji = 2, jpim1
(  513)             spgu(ji,jj) = spgu(ji,jj) + hur(ji,jj) / e2u(ji,jj) * ( bsfd(ji,jj) - bsfd(ji  ,jj-1) )
(  514)             spgv(ji,jj) = spgv(ji,jj) - hvr(ji,jj) / e1v(ji,jj) * ( bsfd(ji,jj) - bsfd(ji-1,jj  ) )
(  515)          END DO
(  516)          ! add the surface pressure gradient trend to the general trend
(  517)          DO jk = 1, jpkm1
(  518)             DO ji = 2, jpim1
(  519)                ua(ji,jj,jk) = ua(ji,jj,jk) - spgu(ji,jj)
(  520)                va(ji,jj,jk) = va(ji,jj,jk) - spgv(ji,jj)
(  521)             END DO
(  522)          END DO
(  523) 
(  524)       END DO
(  525) 
(  526)    END SUBROUTINE dyn_spg_rl
(  527) 
(  528) #else
(  529)    !!----------------------------------------------------------------------
(  530)    !!   'key_dynspg_rl'                                        NO rigid lid
(  531)    !!----------------------------------------------------------------------
(  532)    USE in_out_manager
(  533) CONTAINS
(  534)    SUBROUTINE dyn_spg_rl( kt, kindic )       ! Empty routine
(  535)       if(lwp) WRITE(numout,*) 'dyn_spg_rl: You should not have seen this print! error?', kt, kindic
(  536)    END SUBROUTINE dyn_spg_rl
(  537) #endif
(  538)    
(  539)    !!======================================================================
(  540) END MODULE dynspg_rl


































PGF90 (Version     10.5)          09/23/2010  09:49:43      page 11

(    1) # 1 "dynspg_rl.F90"
(    1) MODULE dynspg_rl
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  dynspg_rl  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) # 529
(  529)    !!----------------------------------------------------------------------
(  530)    !!   'key_dynspg_rl'                                        NO rigid lid
(  531)    !!----------------------------------------------------------------------
(  532)    USE in_out_manager
(  533) CONTAINS
(  534)    SUBROUTINE dyn_spg_rl( kt, kindic )       ! Empty routine
(  535)       if(lwp) WRITE(numout,*) 'dyn_spg_rl: You should not have seen this print! error?', kt, kindic
(  536)    END SUBROUTINE dyn_spg_rl
(  537) # 538
(  538)    
(  539)    !!======================================================================
(  540) END MODULE dynspg_rl
