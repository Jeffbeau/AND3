


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

Filename: solmat.F90

(    1) MODULE solmat
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  solmat  ***
(    4)    !! solver       : construction of the matrix 
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   sol_mat       : Construction of the matrix of used by the elliptic solvers
(    9)    !!   fetsch        :
(   10)    !!   fetmat        :
(   11)    !!   fetstr        :
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and active tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE sol_oce         ! ocean solver
(   17)    USE phycst          ! physical constants
(   18)    USE obc_oce         ! ocean open boundary conditions
(   19)    USE lbclnk          ! lateral boudary conditions
(   20)    USE lib_mpp         ! distributed memory computing
(   21)    USE in_out_manager  ! I/O manager
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC sol_mat     ! routine called by inisol.F90
(   28)    !!----------------------------------------------------------------------
(   29)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solmat.F90,v 1.11 2006/03/20 17:27:14 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36)    SUBROUTINE sol_mat( kt )
(   37)       !!----------------------------------------------------------------------
(   38)       !!                  ***  ROUTINE sol_mat  ***
(   39)       !!
(   40)       !! ** Purpose :   Construction of the matrix of used by the elliptic 
(   41)       !!      solvers (either sor, pcg or feti methods).
(   42)       !!
(   43)       !! ** Method  :   The matrix depends on the type of free surface:
(   44)       !!       * lk_dynspg_rl=T: rigid lid formulation
(   45)       !!      The matrix is built for the barotropic stream function system.
(   46)       !!      a diagonal preconditioning matrix is also defined.






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 2

(   47)       !!       * lk_dynspg_flt=T: free surface formulation
(   48)       !!      The matrix is built for the divergence of the transport system
(   49)       !!      a diagonal preconditioning matrix is also defined.
(   50)       !!        Note that for feti solver (nsolv=3) a specific initialization 
(   51)       !!      is required (call to fetstr.F) for memory allocation and inter-
(   52)       !!      face definition.
(   53)       !! 
(   54)       !! ** Action  : - gcp    : extra-diagonal elements of the matrix
(   55)       !!              - gcdmat : preconditioning matrix (diagonal elements)
(   56)       !!              - gcdprc : inverse of the preconditioning matrix
(   57)       !!
(   58)       !! History :
(   59)       !!   1.0  !  88-04  (G. Madec)  Original code
(   60)       !!        !  91-11  (G. Madec)
(   61)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   62)       !!        !  93-06  (M. Guyon)  suppress pointers
(   63)       !!        !  96-05  (G. Madec)  merge sor and pcg formulations
(   64)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(   65)       !!        !  98-02  (M. Guyon)  FETI method
(   66)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   67)       !!        !  02-11  (C. Talandier, A-M. Treguier) Free surface & Open boundaries
(   68)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   69)       !!----------------------------------------------------------------------
(   70)       !! * Arguments
(   71)       INTEGER, INTENT(in) :: kt
(   72) 
(   73)       !! * Local declarations
(   74)       INTEGER ::   ji, jj                    ! dummy loop indices
(   75)       INTEGER ::   ii, ij, iiend, ijend      ! temporary integers
(   76)       REAL(wp) ::   zcoefs, zcoefw, zcoefe, zcoefn  ! temporary scalars
(   77)       REAL(wp) ::   z2dt, zcoef
(   78)       !!----------------------------------------------------------------------
(   79) 
(   80)       ! FETI method ( nsolv = 3)
(   81)       ! memory allocation and interface definition for the solver
(   82) 
(   83)       IF( nsolv == 3 )   CALL fetstr
(   84) 
(   85)       
(   86)       ! 1. Construction of the matrix
(   87)       ! -----------------------------
(   88)       
(   89)       ! initialize to zero
(   90)       zcoef = 0.e0
(   91)       gcp(:,:,1) = 0.e0
(   92)       gcp(:,:,2) = 0.e0
(   93)       gcp(:,:,3) = 0.e0
(   94)       gcp(:,:,4) = 0.e0
(   95)       
(   96)       gcdprc(:,:) = 0.e0
(   97)       gcdmat(:,:) = 0.e0
(   98)       
(   99)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  100)          z2dt = rdt
(  101)       ELSE
(  102)          z2dt = 2. * rdt
(  103)       ENDIF
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 3

(  105) #if defined key_dynspg_flt && ! defined key_obc
(  106) !!cr      IF( lk_dynspg_flt .AND. .NOT.lk_obc ) THEN   !bug missing lk_dynspg_flt_atsk
(  107) 
(  108)       ! defined the coefficients for free surface elliptic system
(  109) 
(  110)       DO jj = 2, jpjm1
(  111)          DO ji = 2, jpim1
(  112)             zcoef = z2dt * z2dt * grav * rnu * bmask(ji,jj)
(  113)             zcoefs = -zcoef * hv(ji  ,jj-1) * e1v(ji  ,jj-1) / e2v(ji  ,jj-1)    ! south coefficient
(  114)             zcoefw = -zcoef * hu(ji-1,jj  ) * e2u(ji-1,jj  ) / e1u(ji-1,jj  )    ! west coefficient
(  115)             zcoefe = -zcoef * hu(ji  ,jj  ) * e2u(ji  ,jj  ) / e1u(ji  ,jj  )    ! east coefficient
(  116)             zcoefn = -zcoef * hv(ji  ,jj  ) * e1v(ji  ,jj  ) / e2v(ji  ,jj  )    ! north coefficient
(  117)             gcp(ji,jj,1) = zcoefs
(  118)             gcp(ji,jj,2) = zcoefw
(  119)             gcp(ji,jj,3) = zcoefe
(  120)             gcp(ji,jj,4) = zcoefn
(  121)             gcdmat(ji,jj) = e1t(ji,jj) * e2t(ji,jj) * bmask(ji,jj)    &          ! diagonal coefficient
(  122)                &          - zcoefs -zcoefw -zcoefe -zcoefn
(  123)          END DO
(  124)       END DO
(  125)       
(  126) #  elif defined key_dynspg_flt && defined key_obc
(  127) !!cr      ELSEIF( lk_dynspg_flt .AND. lk_obc ) THEN     !bug missing lk_dynspg_flt_atsk 
(  128) 
(  129)       !   defined gcdmat in the case of open boundaries
(  130) 
(  131)       DO jj = 2, jpjm1
(  132)          DO ji = 2, jpim1
(  133)             zcoef = z2dt * z2dt * grav * rnu * bmask(ji,jj)
(  134)             !  south coefficient
(  135)             IF( lp_obc_south .AND. ( jj == njs0p1 ) ) THEN
(  136)                zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)*(1.-vsmsk(ji,1))
(  137)             ELSE
(  138)                zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)
(  139)             END IF
(  140)             gcp(ji,jj,1) = zcoefs
(  141) 
(  142)             !  west coefficient
(  143)             IF( lp_obc_west  .AND. ( ji == niw0p1 ) ) THEN
(  144)                zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)*(1.-uwmsk(jj,1))
(  145)             ELSE
(  146)                zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)
(  147)             END IF
(  148)             gcp(ji,jj,2) = zcoefw
(  149) 
(  150)             !   east coefficient
(  151)             IF( lp_obc_east  .AND. ( ji == nie0 ) ) THEN
(  152)                zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)*(1.-uemsk(jj,1))
(  153)             ELSE
(  154)                zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)
(  155)             END IF
(  156)             gcp(ji,jj,3) = zcoefe
(  157) 
(  158)             !   north coefficient
(  159)             IF( lp_obc_north .AND. ( jj == njn0 ) ) THEN
(  160)                zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)*(1.-vnmsk(ji,1))
(  161)             ELSE
(  162)                zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 4

(  163)             END IF
(  164)             gcp(ji,jj,4) = zcoefn
(  165) 
(  166)             ! diagonal coefficient
(  167)             gcdmat(ji,jj) = e1t(ji,jj)*e2t(ji,jj)*bmask(ji,jj) &
(  168)                             - zcoefs -zcoefw -zcoefe -zcoefn
(  169)          END DO
(  170)       END DO
(  171) 
(  172) #  else
(  173) !!cr      ELSE
(  174) 
(  175)       !   defined the coefficients for bsf elliptic system
(  176)       
(  177)       DO jj = 2, jpjm1
(  178)          DO ji = 2, jpim1
(  179)             zcoefs = -hur(ji  ,jj  ) * e1u(ji  ,jj  ) / e2u(ji  ,jj  ) * bmask(ji,jj)   ! south coefficient
(  180)             zcoefw = -hvr(ji  ,jj  ) * e2v(ji  ,jj  ) / e1v(ji  ,jj  ) * bmask(ji,jj)   ! west coefficient
(  181)             zcoefe = -hvr(ji+1,jj  ) * e2v(ji+1,jj  ) / e1v(ji+1,jj  ) * bmask(ji,jj)   ! east coefficient
(  182)             zcoefn = -hur(ji  ,jj+1) * e1u(ji  ,jj+1) / e2u(ji  ,jj+1) * bmask(ji,jj)   ! north coefficient
(  183)             gcp(ji,jj,1) = zcoefs
(  184)             gcp(ji,jj,2) = zcoefw
(  185)             gcp(ji,jj,3) = zcoefe
(  186)             gcp(ji,jj,4) = zcoefn
(  187)             gcdmat(ji,jj) = -zcoefs -zcoefw -zcoefe -zcoefn                             ! diagonal coefficient
(  188)          END DO
(  189)       END DO
(  190)       
(  191) !!cr  ENDIF
(  192) #endif
(  193) #if defined key_agrif
(  194)        IF (.NOT.AGRIF_ROOT()) THEN
(  195)        
(  196)        IF ( (nbondi == -1)  .OR. (nbondi == 2) ) bmask(2,:)=0.
(  197)        IF ( (nbondi ==  1)  .OR. (nbondi == 2) ) bmask(nlci-1,:)=0.
(  198)        IF ( (nbondj == -1)  .OR. (nbondj == 2) ) bmask(:,2)=0.
(  199)        IF ( (nbondj ==  1)  .OR. (nbondj == 2) ) bmask(:,nlcj-1)=0.
(  200) 
(  201)       DO jj = 2, jpjm1
(  202)          DO ji = 2, jpim1
(  203)             zcoef = z2dt * z2dt * grav * rnu * bmask(ji,jj)
(  204)             !  south coefficient
(  205)             IF( ((nbondj == -1)  .OR. (nbondj == 2)) .AND. ( jj == 3 ) ) THEN
(  206)                zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)*(1.-vmask(ji,jj-1,1))
(  207)             ELSE
(  208)                zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)
(  209)             END IF
(  210)             gcp(ji,jj,1) = zcoefs
(  211) 
(  212)             !  west coefficient
(  213) 	    IF( ( (nbondi == -1)  .OR. (nbondi == 2) ) .AND. ( ji == 3 )  ) THEN
(  214)                zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)*(1.-umask(ji-1,jj,1))
(  215)             ELSE
(  216)                zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)
(  217)             END IF
(  218)             gcp(ji,jj,2) = zcoefw
(  219) 
(  220)             !   east coefficient






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 5

(  221)             IF( ((nbondi == 1)  .OR. (nbondi == 2)) .AND. ( ji == nlci-2 ) ) THEN
(  222)                zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)*(1.-umask(ji,jj,1))
(  223)             ELSE
(  224)                zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)
(  225)             END IF
(  226)             gcp(ji,jj,3) = zcoefe
(  227) 
(  228)             !   north coefficient
(  229)             IF( ((nbondj == 1)  .OR. (nbondj == 2)) .AND. ( jj == nlcj-2 ) ) THEN
(  230)                zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)*(1.-vmask(ji,jj,1))
(  231)             ELSE
(  232)                zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)
(  233)             END IF
(  234)             gcp(ji,jj,4) = zcoefn
(  235) 
(  236)             ! diagonal coefficient
(  237)             gcdmat(ji,jj) = e1t(ji,jj)*e2t(ji,jj)*bmask(ji,jj) &
(  238)                             - zcoefs -zcoefw -zcoefe -zcoefn
(  239)          END DO
(  240)       END DO
(  241)       
(  242)        ENDIF
(  243) #endif
(  244) 
(  245)       ! 2. Boundary conditions 
(  246)       ! ----------------------
(  247)       
(  248)       ! Cyclic east-west boundary conditions
(  249)       !     ji=2 is the column east of ji=jpim1 and reciprocally,
(  250)       !     ji=jpim1 is the column west of ji=2
(  251)       !     all the coef are already set to zero as bmask is initialized to
(  252)       !     zero for ji=1 and ji=jpj in dommsk.
(  253)       
(  254)       ! Symetrical conditions
(  255)       ! free surface: no specific action
(  256)       ! bsf system: n-s gradient of bsf = 0 along j=2 (perhaps a bug !!!!!!)
(  257)       ! the diagonal coefficient of the southern grid points must be modify to
(  258)       ! account for the existence of the south symmetric bassin.
(  259)       
(  260) !!cr      IF( .NOT.lk_dynspg_flt ) THEN   !bug missing lk_dynspg_flt_atsk
(  261) #if ! defined key_dynspg_flt
(  262)       IF( nperio == 2 ) THEN
(  263)          DO ji = 1, jpi
(  264)             IF( bmask(ji,2) /= 0.e0 ) THEN
(  265)                zcoefs = - hur(ji,2)*e1u(ji,2)/e2u(ji,2)
(  266)                gcdmat(ji,2) = gcdmat(ji,2) - zcoefs
(  267)             ENDIF
(  268)          END DO
(  269)       ENDIF
(  270) !!cr      ENDIF
(  271) #endif
(  272)       
(  273)       ! North fold boundary condition
(  274)       !     all the coef are already set to zero as bmask is initialized to
(  275)       !     zero on duplicated lignes and portion of lignes
(  276)       
(  277)       ! 3. Preconditioned matrix
(  278)       ! ------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 6

(  279)       
(  280)       IF( nsolv /= 3 ) THEN
(  281)          
(  282)          ! SOR and PCG solvers
(  283)          DO jj = 1, jpj
(  284)             DO ji = 1, jpi
(  285)                IF( bmask(ji,jj) /= 0.e0 )   gcdprc(ji,jj) = 1.e0 / gcdmat(ji,jj)
(  286)             END DO
(  287)          END DO
(  288)          
(  289)          gcp(:,:,1) = gcp(:,:,1) * gcdprc(:,:)
(  290)          gcp(:,:,2) = gcp(:,:,2) * gcdprc(:,:)
(  291)          gcp(:,:,3) = gcp(:,:,3) * gcdprc(:,:)
(  292)          gcp(:,:,4) = gcp(:,:,4) * gcdprc(:,:)
(  293)          IF( ( nsolv == 2 ) .OR. ( nsolv == 4 ) )  gccd(:,:) = sor * gcp(:,:,2)
(  294) 
(  295)          IF( nsolv == 4 ) THEN
(  296)             CALL lbc_lnk_e( gcp   (:,:,1), c_solver_pt, 1. )   ! lateral boundary conditions
(  297)             CALL lbc_lnk_e( gcp   (:,:,2), c_solver_pt, 1. )   ! lateral boundary conditions
(  298)             CALL lbc_lnk_e( gcp   (:,:,3), c_solver_pt, 1. )   ! lateral boundary conditions
(  299)             CALL lbc_lnk_e( gcp   (:,:,4), c_solver_pt, 1. )   ! lateral boundary conditions
(  300)             CALL lbc_lnk_e( gcdprc(:,:)  , c_solver_pt, 1. )   ! lateral boundary conditions
(  301)             CALL lbc_lnk_e( gcdmat(:,:)  , c_solver_pt, 1. )   ! lateral boundary conditions         
(  302)             IF( npolj /= 0 ) CALL sol_exd( gcp , c_solver_pt ) ! switch northernelements
(  303)          END IF
(  304) 
(  305)       ELSE
(  306)          
(  307)          ! FETI method
(  308)          ! if feti solver : gcdprc is a mask for the non-overlapping
(  309)          !   data structuring
(  310)          
(  311)          DO jj = 1, jpj
(  312)             DO ji = 1, jpi
(  313)                IF( bmask(ji,jj) /= 0.e0 ) THEN
(  314)                   gcdprc(ji,jj) = 1.e0
(  315)                ELSE
(  316)                   gcdprc(ji,jj) = 0.e0
(  317)                ENDIF
(  318)             END DO
(  319)          END DO
(  320)          
(  321)          ! so "common" line & "common" column have to be !=0 except on global
(  322)          !   domain boundaries
(  323)          ! pbs with nbondi if nperio != 2 ?
(  324)          !   ii = nldi-1
(  325)          ! pb with nldi value if jperio==1 : nbondi modifyed at the end
(  326)          !   of inimpp.F => pb
(  327)          ! pb with periodicity conditions : iiend, ijend
(  328)          
(  329)          ijend = nlej
(  330)          iiend = nlei
(  331)          IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 )   iiend = nlci - jpreci
(  332)          ii = jpreci
(  333)          
(  334)          ! case number 1
(  335)          
(  336)          IF( nbondi /= -1 .AND. nbondi /= 2 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 7

(  337)             DO jj = 1, ijend
(  338)                IF( fmask(ii,jj,1) == 1. ) gcdprc(ii,jj) = 1.
(  339)             END DO
(  340)          ENDIF
(  341)          
(  342)          ! case number 2
(  343)          
(  344)          IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
(  345)             DO jj = 1, ijend
(  346)                IF( fmask(ii,jj,1) == 1. ) gcdprc(ii,jj) = 1.
(  347)             END DO
(  348)          ENDIF
(  349)          
(  350)          !      ij = nldj-1
(  351)          ! pb with nldi value if jperio==1 : nbondi modifyed at the end
(  352)          !   of inimpp.F => pb, here homogeneisation...
(  353)          
(  354)          ij = jprecj
(  355)          IF( nbondj /= -1 .AND. nbondj /= 2 ) THEN
(  356)             DO ji = 1, iiend
(  357)                IF( fmask(ji,ij,1) == 1. ) gcdprc(ji,ij) = 1.
(  358)             END DO
(  359)          ENDIF
(  360)       ENDIF
(  361)       
(  362)       
(  363)       ! 4. Initialization the arrays used in pcg
(  364)       ! ----------------------------------------
(  365)       gcx  (:,:) = 0.e0
(  366)       gcxb (:,:) = 0.e0
(  367)       gcb  (:,:) = 0.e0
(  368)       gcr  (:,:) = 0.e0
(  369)       gcdes(:,:) = 0.e0
(  370)       gccd (:,:) = 0.e0
(  371)       
(  372)       ! FETI method
(  373)       IF( nsolv == 3 ) THEN
(  374)          CALL fetmat       ! Matrix treatment : Neumann condition, inverse computation
(  375)          CALL fetsch       ! data framework for the Schur Dual solver
(  376)       ENDIF
(  377)       
(  378)    END SUBROUTINE sol_mat
(  379) 
(  380) 
(  381)    SUBROUTINE sol_exd( pt3d, cd_type )
(  382)       !!----------------------------------------------------------------------
(  383)       !!                  ***  routine sol_exd  ***
(  384)       !!                  
(  385)       !! ** Purpose :   Reorder gcb coefficient on the extra outer  halo 
(  386)       !!                at north fold in case of T or F pivot
(  387)       !!
(  388)       !! ** Method  :   Perform a circular permutation of the coefficients on 
(  389)       !!                the total area strictly above the pivot point,
(  390)       !!                and on the semi-row of the pivot point   
(  391)       !!                
(  392)       !! History :
(  393)       !!   9.0  !  05-09  (R. Benshila)  original routine
(  394)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 8

(  395)       !! * Arguments
(  396)       CHARACTER(len=1) , INTENT( in ) ::   &
(  397)          cd_type       ! define the nature of pt2d array grid-points
(  398)          !             !  = T , U , V , F , W 
(  399)          !             !  = S : T-point, north fold treatment
(  400)          !             !  = G : F-point, north fold treatment
(  401)          !             !  = I : sea-ice velocity at F-point with index shift
(  402)       REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4), INTENT( inout ) ::   &
(  403)          pt3d          ! 2D array on which the boundary condition is applied
(  404) 
(  405)       !! * Local variables
(  406)       INTEGER  ::   ji, jk      ! dummy loop indices
(  407)       INTEGER  ::   iloc                ! temporary integers
(  408)       REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4) ::   &
(  409)          ztab          ! 2D array on which the boundary condition is applied
(  410)       !!----------------------------------------------------------------------
(  411) 
(  412)       ztab = pt3d
(  413) 
(  414)       ! north fold treatment
(  415)       ! -----------------------
(  416)   
(  417)       SELECT CASE ( npolj )
(  418)          
(  419)          CASE ( 3 , 4 )   !  T pivot
(  420)          iloc = jpiglo/2 +1 
(  421)             
(  422)             SELECT CASE ( cd_type )
(  423)   
(  424)             CASE ( 'T', 'S', 'U', 'W' )
(  425)                DO jk =1, 4
(  426)                   DO ji = 1-jpr2di, nlci+jpr2di
(  427)                      pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
(  428)                   ENDDO
(  429)                ENDDO
(  430) 
(  431)               DO jk =1, 4
(  432)                   DO ji = nlci+jpr2di, 1-jpr2di,  -1
(  433)                      IF( ( ji .LT. mi0(iloc) .AND. mi0(iloc) /= 1 ) &
(  434)                        & .OR. ( mi0(iloc) == jpi+1 ) ) EXIT
(  435)                      pt3d(ji,nlcj-1,jk) = ztab(ji,nlcj-1,jk+3-2*MOD(jk+3,4))
(  436)                   ENDDO
(  437)                ENDDO
(  438) 
(  439)             CASE ( 'F' ,'G' , 'I', 'V' )
(  440)                DO jk =1, 4
(  441)                   DO ji = 1-jpr2di, nlci+jpr2di
(  442)                      pt3d(ji,nlcj-1:nlcj+jpr2dj,jk) = ztab(ji,nlcj-1:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
(  443)                   ENDDO
(  444)                ENDDO
(  445) 
(  446)             END SELECT   ! cd_type
(  447)   
(  448)          CASE ( 5 , 6 )                 ! F pivot
(  449)           iloc=jpiglo/2
(  450) 
(  451)             SELECT CASE (cd_type )
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 9

(  453)             CASE ( 'T'  ,'S', 'U', 'W')
(  454)                DO jk =1, 4
(  455)                   DO ji = 1-jpr2di, nlci+jpr2di
(  456)                      pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
(  457)                   ENDDO
(  458)                ENDDO
(  459) 
(  460)             CASE ( 'F' ,'G' , 'I', 'V' )
(  461)                DO jk =1, 4
(  462)                   DO ji = 1-jpr2di, nlci+jpr2di
(  463)                      pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
(  464)                   ENDDO
(  465)                ENDDO
(  466)                DO jk =1, 4
(  467)                   DO ji = nlci+jpr2di, 1-jpr2di,  -1
(  468)                     IF ( ( ji .LT. mi0(iloc) .AND. mi0(iloc) /= 1 ) &
(  469)                        & .OR. ( mi0(iloc) == jpi+1 ) ) EXIT
(  470)                     pt3d(ji,nlcj-1,jk) = ztab(ji,nlcj-1,jk+3-2*MOD(jk+3,4))
(  471)                   ENDDO
(  472)                ENDDO
(  473) 
(  474)             END SELECT   ! cd_type
(  475) 
(  476)          END SELECT   ! npolj
(  477)   
(  478)    END SUBROUTINE sol_exd
(  479) 
(  480) #if defined key_feti
(  481) 
(  482)    SUBROUTINE fetstr
(  483)       !!---------------------------------------------------------------------
(  484)       !!                  ***  ROUTINE fetstr  ***
(  485)       !!               
(  486)       !! ** Purpose :   Construction of the matrix of the barotropic stream 
(  487)       !!      function system.
(  488)       !!      Finite Elements Tearing & Interconnecting (FETI) approach
(  489)       !!      Memory allocation and interface definition for the solver
(  490)       !!
(  491)       !! ** Method :
(  492)       !!
(  493)       !! References :
(  494)       !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
(  495)       !!      A domain decomposition solver to compute the barotropic 
(  496)       !!      component of an OGCM in the parallel processing field.
(  497)       !!      Ocean Modelling, issue 105, december 94.
(  498)       !!
(  499)       !! History :
(  500)       !!        !  98-02 (M. Guyon)  Original code
(  501)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  502)       !!----------------------------------------------------------------------
(  503)       !! * Modules used
(  504)       USE lib_feti            ! feti librairy
(  505)       !! * Local declarations
(  506)       INTEGER ::   iiend, ijend, iperio   ! temporary integers
(  507)       !!---------------------------------------------------------------------
(  508)       
(  509)       
(  510)       ! Preconditioning technics of the Dual Schur Operator






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 10

(  511)       ! <= definition of the Coarse Grid solver
(  512)       ! <= dimension of the nullspace of the local operators
(  513)       ! <= Neumann boundaries conditions
(  514)       
(  515)       ! 0. Initializations
(  516)       ! ------------------
(  517)       
(  518)       ndkerep = 1
(  519) 
(  520)       ! initialization of the superstructures management
(  521) 
(  522)       malxm = 1
(  523)       malim = 1
(  524) 
(  525)       ! memory space for the pcpg associated with the FETI dual formulation
(  526)       ! ndkerep is associated to the list of rigid modes, 
(  527)       ! ndkerep == 1 because the Dual Operator
(  528)       ! is a first order operator due to SPG elliptic Operator is a
(  529)       ! second order operator
(  530) 
(  531)       nim = 50
(  532)       nim = nim + ndkerep
(  533)       nim = nim + 2*jpi + 2*jpj
(  534)       nim = nim + jpi*jpj
(  535) 
(  536)       nxm = 33
(  537)       nxm = nxm + 4*jpnij
(  538)       nxm = nxm + 19*(jpi+jpj)
(  539)       nxm = nxm + 13*jpi*jpj
(  540)       nxm = nxm + jpi*jpi*jpj
(  541) 
(  542)       ! krylov space memory
(  543) 
(  544)       iperio = 0
(  545)       IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6) iperio = 1
(  546)       nxm = nxm + 3*(jpnij-jpni)*jpi
(  547)       nxm = nxm + 3*(jpnij-jpnj+iperio)*jpj
(  548)       nxm = nxm + 2*(jpi+jpj)*(jpnij-jpni)*jpi
(  549)       nxm = nxm + 2*(jpi+jpj)*(jpnij-jpnj+iperio)*jpj
(  550) 
(  551)       ! Resolution with the Schur dual Method ( frontal and local solver by
(  552)       ! blocks
(  553)       ! Case with a local symetrical matrix
(  554)       ! The local matrix is stored in a multi-column form
(  555)       ! The total number of nodes for this subdomain is named "noeuds"
(  556) 
(  557)       noeuds = jpi*jpj
(  558)       nifmat = jpi-1
(  559)       njfmat = jpj-1
(  560)       nelem = nifmat*njfmat
(  561)       npe = 4
(  562)       nmorse = 5*noeuds
(  563)       
(  564)       ! 1. mesh building
(  565)       ! ----------------
(  566)       
(  567)       ! definition of specific information for a subdomain
(  568)       !  narea           : subdomain number = processor number +1 






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 11

(  569)       !  ninterf         : neighbour subdomain number
(  570)       !  nni             : interface point number
(  571)       !  ndvois array    : neighbour subdomain list 
(  572)       !  maplistin array : node pointer at each interface
(  573)       !  maplistin array : concatened list of interface nodes
(  574)       
(  575)       !  messag coding is necessary by interface type for avoid collision
(  576)       !  if nperio == 1 
(  577)       
(  578)       !  lint  array     : indoor interface list / type
(  579)       !  lext  array     : outdoor interface list / type
(  580)       
(  581)       !  domain with jpniXjpnj subdomains
(  582)       
(  583)       CALL feti_inisub(nifmat,njfmat,nbondi,nbondj,nperio,   &
(  584)           nbsw,nbnw,nbse,nbne,ninterf,ninterfc,nni,nnic)
(  585) 
(  586)       CALL feti_creadr(malim,malimax,nim,3*ninterf ,mandvois ,'ndvois' )
(  587)       CALL feti_creadr(malim,malimax,nim,3*ninterfc,mandvoisc,'ndvoisc')
(  588)       CALL feti_creadr(malim,malimax,nim,ninterfc+1,maplistin,'plistin')
(  589)       CALL feti_creadr(malim,malimax,nim,nnic      ,malistin ,'listin' )
(  590) 
(  591)       ! pb with periodicity conditions : iiend, ijend
(  592) 
(  593)       ijend = nlej
(  594)       iiend = nlei
(  595)       IF (jperio == 1) iiend = nlci - jpreci
(  596) 
(  597)       CALL feti_subound(nifmat,njfmat,nldi,iiend,nldj,ijend,   &
(  598)           narea,nbondi,nbondj,nperio,   &
(  599)           ninterf,ninterfc,   &
(  600)           nowe,noea,noso,nono,   &
(  601)           nbsw,nbnw,nbse,nbne,   &
(  602)           npsw,npnw,npse,npne,   &
(  603)           mfet(mandvois),mfet(mandvoisc),   &
(  604)           mfet(maplistin),nnic,mfet(malistin) )
(  605) 
(  606)    END SUBROUTINE fetstr
(  607) 
(  608) 
(  609)    SUBROUTINE fetmat
(  610)       !!---------------------------------------------------------------------
(  611)       !!                   ***  ROUTINE fetmat  ***
(  612)       !!            
(  613)       !! ** Purpose :   Construction of the matrix of the barotropic stream 
(  614)       !!      function system.
(  615)       !!      Finite Elements Tearing & Interconnecting (FETI) approach
(  616)       !!      Matrix treatment : Neumann condition, inverse computation
(  617)       !!
(  618)       !! ** Method :
(  619)       !!
(  620)       !! References :
(  621)       !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
(  622)       !!      A domain decomposition solver to compute the barotropic 
(  623)       !!      component of an OGCM in the parallel processing field.
(  624)       !!      Ocean Modelling, issue 105, december 94.
(  625)       !!
(  626)       !! History :






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 12

(  627)       !!        !  98-02 (M. Guyon)  Original code
(  628)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  629)       !!----------------------------------------------------------------------
(  630)       !! * Modules used
(  631)       USE lib_feti            ! feti librairy
(  632)       !! * Local declarations
(  633)       INTEGER ::   ji, jj, jk, jl
(  634)       INTEGER ::   iimask(jpi,jpj)
(  635)       INTEGER ::   iiend, ijend
(  636)       REAL(wp) ::   zres, zres2, zdemi
(  637)       !!---------------------------------------------------------------------
(  638) 
(  639)       ! Matrix computation
(  640)       ! ------------------
(  641) 
(  642)       CALL feti_creadr(malxm,malxmax,nxm,nmorse,maan,'matrice a')
(  643) 
(  644)       nnitot = nni
(  645) 
(  646)       CALL mpp_sum( nnitot, 1, numit0ete )
(  647)       CALL feti_creadr(malxm,malxmax,nxm,npe*npe,maae,'ae')
(  648) 
(  649)       ! initialisation of the local barotropic matrix
(  650)       ! local boundary conditions on the halo
(  651) 
(  652)       CALL lbc_lnk( gcp(:,:,1), 'F', 1)
(  653)       CALL lbc_lnk( gcp(:,:,2), 'F', 1) 
(  654)       CALL lbc_lnk( gcp(:,:,3), 'F', 1) 
(  655)       CALL lbc_lnk( gcp(:,:,4), 'F', 1) 
(  656)       CALL lbc_lnk( gcdmat    , 'T', 1)
(  657) 
(  658)       ! Neumann conditions
(  659)       ! initialisation of the integer Neumann Mask
(  660) 
(  661)       CALL feti_iclr(jpi*jpj,iimask)
(  662)       DO jj = 1, jpj
(  663)          DO ji = 1, jpi
(  664)             iimask(ji,jj) = INT( gcdprc(ji,jj) )
(  665)          END DO
(  666)       END DO
(  667) 
(  668)       ! regularization of the local matrix
(  669) 
(  670)       DO jj = 1, jpj
(  671)          DO ji = 1, jpi
(  672)             gcdmat(ji,jj) = gcdmat(ji,jj) * gcdprc(ji,jj) + 1. - gcdprc(ji,jj)
(  673)          END DO
(  674)       END DO
(  675) 
(  676)       DO jk = 1, 4
(  677)          DO jj = 1, jpj
(  678)             DO ji = 1, jpi
(  679)                gcp(ji,jj,jk) = gcp(ji,jj,jk) * gcdprc(ji,jj)
(  680)             END DO
(  681)          END DO
(  682)       END DO
(  683)       
(  684)       ! implementation of the west, east, north & south Neumann conditions






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 13

(  685) 
(  686)       zdemi  = 0.5
(  687) 
(  688)       ! pb with periodicity conditions : iiend, ijend
(  689) 
(  690)       ijend = nlej
(  691)       iiend = nlei
(  692)       IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) iiend = nlci - jpreci
(  693) 
(  694)       IF( nbondi == 2 .AND. (nperio /= 1 .OR. nperio /= 4 .OR. nperio == 6) ) THEN
(  695) 
(  696)          ! with the periodicity : no east/west interface if nbondi = 2
(  697)          ! and nperio != 1
(  698) 
(  699)       ELSE 
(  700)          ! west
(  701)          IF( nbondi /= -1 ) THEN 
(  702)             DO jj = 1, jpj
(  703)                IF( iimask(1,jj) /= 0 ) THEN
(  704)                   gcp(1,jj,2) = 0.e0
(  705)                   gcp(1,jj,1) = zdemi * gcp(1,jj,1)
(  706)                   gcp(1,jj,4) = zdemi * gcp(1,jj,4)
(  707)                ENDIF
(  708)             END DO
(  709)             DO jj = 1, jpj
(  710)                IF( iimask(1,jj) /= 0 ) THEN
(  711)                   gcdmat(1,jj) = - ( gcp(1,jj,1) + gcp(1,jj,2) + gcp(1,jj,3) + gcp(1,jj,4) )
(  712)                ENDIF
(  713)             END DO
(  714)          ENDIF
(  715)          ! east
(  716)          IF( nbondi /= 1 ) THEN 
(  717)             DO jj = 1, jpj
(  718)                IF( iimask(iiend,jj) /= 0 ) THEN
(  719)                   gcp(iiend,jj,3) = 0.e0
(  720)                   gcp(iiend,jj,1) = zdemi * gcp(iiend,jj,1)
(  721)                   gcp(iiend,jj,4) = zdemi * gcp(iiend,jj,4)
(  722)                ENDIF
(  723)             END DO
(  724)             DO jj = 1, jpj
(  725)                IF( iimask(iiend,jj) /= 0 ) THEN
(  726)                   gcdmat(iiend,jj) = - ( gcp(iiend,jj,1) + gcp(iiend,jj,2)   &
(  727)                                        + gcp(iiend,jj,3) + gcp(iiend,jj,4) )
(  728)                ENDIF
(  729)             END DO
(  730)          ENDIF
(  731)       ENDIF
(  732) 
(  733)       ! south
(  734)       IF( nbondj /= -1 .AND. nbondj /= 2 ) THEN 
(  735)          DO ji = 1, jpi
(  736)             IF( iimask(ji,1) /= 0 ) THEN
(  737)                gcp(ji,1,1) = 0.e0
(  738)                gcp(ji,1,2) = zdemi * gcp(ji,1,2)
(  739)                gcp(ji,1,3) = zdemi * gcp(ji,1,3)
(  740)             ENDIF
(  741)          END DO
(  742)          DO ji = 1, jpi






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 14

(  743)             IF( iimask(ji,1) /= 0 ) THEN
(  744)                gcdmat(ji,1) = - ( gcp(ji,1,1) + gcp(ji,1,2) + gcp(ji,1,3) + gcp(ji,1,4) )
(  745)             ENDIF
(  746)          END DO
(  747)       ENDIF
(  748)       
(  749)       ! north
(  750)       IF( nbondj /= 1 .AND. nbondj /= 2 ) THEN 
(  751)          DO ji = 1, jpi
(  752)             IF( iimask(ji,ijend) /= 0 ) THEN
(  753)                gcp(ji,ijend,4) = 0.e0
(  754)                gcp(ji,ijend,2) = zdemi * gcp(ji,ijend,2) 
(  755)                gcp(ji,ijend,3) = zdemi * gcp(ji,ijend,3)
(  756)             ENDIF
(  757)          END DO
(  758)          DO ji = 1, jpi
(  759)             IF( iimask(ji,ijend) /= 0 ) THEN
(  760)                gcdmat(ji,ijend) = - ( gcp(ji,ijend,1) + gcp(ji,ijend,2)   &
(  761)                                     + gcp(ji,ijend,3) + gcp(ji,ijend,4) )
(  762)             ENDIF
(  763)          END DO
(  764)       ENDIF
(  765) 
(  766)       ! matrix terms are  saved in FETI solver arrays
(  767)       CALL feti_vmov(noeuds,gcp(1,1,1),wfeti(maan))
(  768)       CALL feti_vmov(noeuds,gcp(1,1,2),wfeti(maan+noeuds))
(  769)       CALL feti_vmov(noeuds,gcdmat,wfeti(maan+2*noeuds))
(  770)       CALL feti_vmov(noeuds,gcp(1,1,3),wfeti(maan+3*noeuds))
(  771)       CALL feti_vmov(noeuds,gcp(1,1,4),wfeti(maan+4*noeuds))
(  772) 
(  773)       ! construction of Dirichlet liberty degrees array
(  774)       CALL feti_subdir(nifmat,njfmat,noeuds,ndir,iimask)
(  775)       CALL feti_creadr(malim,malimax,nim,ndir,malisdir,'lisdir')
(  776)       CALL feti_listdir(jpi,jpj,iimask,ndir,mfet(malisdir))
(  777) 
(  778)       ! stop onto  matrix term for Dirichlet conditions
(  779)       CALL feti_blomat(nifmat+1,njfmat+1,wfeti(maan),ndir,mfet(malisdir))
(  780) 
(  781)       ! reservation of factorized diagonal blocs and temporary array for
(  782)       ! factorization
(  783)       npblo = (njfmat+1) * (nifmat+1) * (nifmat+1)
(  784)       ndimax = nifmat+1
(  785) 
(  786)       CALL feti_creadr(malxm,malxmax,nxm,npblo,mablo,'blo')
(  787)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,madia,'dia')
(  788)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,mav,'v')
(  789)       CALL feti_creadr(malxm,malxmax,nxm,ndimax*ndimax,mautil,'util')
(  790) 
(  791)       ! stoping the rigid modes
(  792) 
(  793)       ! the number of rigid modes =< Max [dim(Ker(Ep))]
(  794)       !                                p=1,Np
(  795) 
(  796)       CALL feti_creadr(malim,malimax,nim,ndkerep,malisblo,'lisblo')
(  797) 
(  798)       ! Matrix factorization
(  799) 
(  800)       CALL feti_front(noeuds,nifmat+1,njfmat+1,wfeti(maan),npblo,   &






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 15

(  801)           wfeti(mablo),wfeti(madia),   &
(  802)           wfeti(mautil),wfeti(mav),ndlblo,mfet(malisblo),ndkerep)
(  803)       CALL feti_prext(noeuds,wfeti(madia))
(  804) 
(  805)       ! virtual dealloc => we have to see for a light f90 version
(  806)       ! the super structure is removed to clean the coarse grid
(  807)       ! solver structure
(  808) 
(  809)       malxm = madia
(  810)       CALL feti_vclr(noeuds,wfeti(madia))
(  811)       CALL feti_vclr(noeuds,wfeti(mav))
(  812)       CALL feti_vclr(ndimax*ndimax,wfeti(mautil))
(  813) 
(  814)       ! ndlblo is the dimension of the local nullspace .=<. the size of the 
(  815)       ! memory of the superstructure associated to the nullspace : ndkerep
(  816)       ! ndkerep is introduced to avoid messages "out of bounds" when memory
(  817)       ! is checked
(  818) 
(  819)       ! copy matrix for Dirichlet condition
(  820) 
(  821)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,miax,'x')
(  822)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,may,'y')
(  823)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,maz,'z')
(  824) 
(  825)       ! stoping the rigid modes
(  826) 
(  827)       ! ndlblo is the dimension of the local nullspace .=<. the size of the 
(  828)       ! memory of the superstructure associated to the nullspace : ndkerep
(  829)       ! ndkerep is introduced to avoid messages "out of bounds" when memory
(  830)       ! is checked
(  831) 
(  832)       CALL feti_creadr(malxm,malxmax,nxm,ndkerep*noeuds,mansp,'nsp')
(  833)       CALL feti_blomat1(nifmat+1,njfmat+1,wfeti(maan),ndlblo,   &
(  834)           mfet(malisblo),wfeti(mansp))      
(  835) 
(  836)       ! computation of operator kernel
(  837) 
(  838)       CALL feti_nullsp(noeuds,nifmat+1,njfmat+1,npblo,wfeti(mablo),   &
(  839)           wfeti(maan),ndlblo,mfet(malisblo),wfeti(mansp),   &
(  840)           wfeti(maz))
(  841) 
(  842)       ! test of the factorisation onto each sub domain
(  843) 
(  844)       CALL feti_init(noeuds,wfeti(may))
(  845)       CALL feti_blodir(noeuds,wfeti(may),ndir,mfet(malisdir))
(  846)       CALL feti_blodir(noeuds,wfeti(may),ndlblo,mfet(malisblo))
(  847)       CALL feti_vclr(noeuds,wfeti(miax))
(  848)       CALL feti_resloc(noeuds,nifmat+1,njfmat+1,wfeti(maan),npblo,   &
(  849)           wfeti(mablo),wfeti(may),wfeti(miax),wfeti(maz)) 
(  850)       CALL feti_proax(noeuds,nifmat+1,njfmat+1,wfeti(maan),wfeti(miax),   &
(  851)           wfeti(maz))
(  852)       CALL feti_blodir(noeuds,wfeti(maz),ndlblo,mfet(malisblo))
(  853)       CALL feti_vsub(noeuds,wfeti(may),wfeti(maz),wfeti(maz))
(  854) 
(  855)       zres2 = 0.e0
(  856)       DO jl = 1, noeuds
(  857)          zres2 = zres2 + wfeti(may+jl-1) * wfeti(may+jl-1)
(  858)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 16

(  859)       CALL mpp_sum(zres2,1,zres)
(  860) 
(  861)       res2 = 0.e0
(  862)       DO jl = 1, noeuds
(  863)          res2 = res2 + wfeti(maz+jl-1) * wfeti(maz+jl-1)
(  864)       END DO 
(  865)       res2 = res2 / zres2
(  866)       CALL mpp_sum(res2,1,zres)
(  867) 
(  868)       res2 = SQRT(res2)
(  869)       IF(lwp) WRITE(numout,*) 'global residu : sqrt((Ax-b,Ax-b)/(b.b)) =', res2
(  870) 
(  871)       IF( res2 > (eps/100.) ) THEN 
(  872)          IF(lwp) WRITE (numout,*) 'eps is :',eps
(  873)          IF(lwp) WRITE (numout,*) 'factorized matrix precision :',res2
(  874)          STOP 
(  875)       ENDIF
(  876) 
(  877)    END SUBROUTINE fetmat
(  878) 
(  879) 
(  880)    SUBROUTINE fetsch
(  881)       !!---------------------------------------------------------------------
(  882)       !!                  ***  ROUTINE fetsch  ***
(  883)       !!        
(  884)       !! ** Purpose :
(  885)       !!     Construction of the matrix of the barotropic stream function
(  886)       !!     system.
(  887)       !!     Finite Elements Tearing & Interconnecting (FETI) approach
(  888)       !!     Data framework for the Schur Dual solve
(  889)       !!
(  890)       !! ** Method :
(  891)       !!
(  892)       !! References :
(  893)       !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
(  894)       !!      A domain decomposition solver to compute the barotropic 
(  895)       !!      component of an OGCM in the parallel processing field.
(  896)       !!      Ocean Modelling, issue 105, december 94.
(  897)       !!
(  898)       !! History :
(  899)       !!        !  98-02 (M. Guyon)  Original code
(  900)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  901)       !!----------------------------------------------------------------------
(  902)       !! * Modules used
(  903)       USE lib_feti            ! feti librairy
(  904)       !! * Local declarations
(  905)       !!---------------------------------------------------------------------
(  906) 
(  907)       !  computing weights for the conform construction
(  908) 
(  909)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,mapoids ,'poids' )
(  910)       CALL feti_creadr(malxm,malxmax,nxm,nnic  ,mabufin ,'bufin' )
(  911)       CALL feti_creadr(malxm,malxmax,nxm,nnic  ,mabufout,'bufout')
(  912) 
(  913) !!    CALL feti_poids(ninterfc,mfet(mandvoisc),mfet(maplistin),nnic,   &
(  914) !!        mfet(malistin),narea,noeuds,wfeti(mapoids),wfeti(mabufin),   &
(  915) !!        wfeti(mabufout) )
(  916)       CALL feti_poids(ninterfc,                                nnic,   &






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 17

(  917)           mfet(malistin),      noeuds,wfeti(mapoids)                )
(  918) 
(  919) 
(  920)       ! Schur dual arrays
(  921)  
(  922)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,mabitw,'bitw') 
(  923)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,mautilu,'utilu') 
(  924)       CALL feti_creadr(malxm,malxmax,nxm,nni,malambda,'lambda') 
(  925)       CALL feti_creadr(malxm,malxmax,nxm,nni,mag,'g') 
(  926)       CALL feti_creadr(malxm,malxmax,nxm,nni,mapg,'pg') 
(  927)       CALL feti_creadr(malxm,malxmax,nxm,nni,mamg,'mg') 
(  928)       CALL feti_creadr(malxm,malxmax,nxm,nni,maw,'w') 
(  929)       CALL feti_creadr(malxm,malxmax,nxm,nni,madw,'dw')
(  930) 
(  931)       !  coarse grid solver dimension and arrays
(  932) 
(  933)       nitmaxete = ndlblo
(  934)       CALL  mpp_sum(nitmaxete,1,numit0ete)
(  935) 
(  936)       nitmaxete = nitmaxete + 1
(  937)       CALL feti_creadr(malxm,malxmax,nxm,ndkerep,maxnul,'xnul')
(  938)       CALL feti_creadr(malxm,malxmax,nxm,ndkerep,maynul,'ynul')
(  939)       CALL feti_creadr(malxm,malxmax,nxm,ndkerep,maeteg,'eteg')
(  940)       CALL feti_creadr(malxm,malxmax,nxm,ndkerep,maeteag,'eteag')
(  941)       CALL feti_creadr(malxm,malxmax,nxm,ndkerep*nitmaxete,maeted,'eted')
(  942)       CALL feti_creadr(malxm,malxmax,nxm,ndkerep*nitmaxete,maetead,'etead')
(  943)       CALL feti_creadr(malxm,malxmax,nxm,nitmaxete,maeteadd,'eteadd')
(  944)       CALL feti_creadr(malxm,malxmax,nxm,nitmaxete,maetegamm,'etegamm')
(  945)       CALL feti_creadr(malxm,malxmax,nxm,nni,maetew,'etew') 
(  946)       CALL feti_creadr(malxm,malxmax,nxm,noeuds,maetev,'etev') 
(  947) 
(  948)       ! construction of semi interface arrays
(  949) 
(  950)       CALL feti_creadr(malim,malimax,nim,ninterf+1,maplistih,'plistih')
(  951) !!    CALL feti_halfint(ninterf,mfet(mandvois),mfet(maplistin),nni,   &
(  952) !!        mfet(maplistih),nnih,narea)
(  953)       CALL feti_halfint(ninterf,mfet(mandvois),mfet(maplistin),       & 
(  954)           mfet(maplistih),nnih      )
(  955) 
(  956)       CALL feti_creadr(malxm,malxmax,nxm,nnih,magh,'gh') 
(  957) 
(  958)       ! Schur Dual Method
(  959) 
(  960)       nmaxd = nnitot / 2
(  961) 
(  962)       !  computation of the remain array for descent directions
(  963) 
(  964)       nmaxd = min0(nmaxd,(nxm-nitmaxete-malxm)/(2*nnih+3))
(  965)       CALL mpp_min(nmaxd,1,numit0ete)
(  966) 
(  967)       nitmax = nnitot/2
(  968)       epsilo = eps
(  969)       ntest = 0
(  970) 
(  971)       ! Krylov space construction
(  972) 
(  973)       CALL feti_creadr(malxm,malxmax,nxm,nnih*nmaxd,mawj,'wj') 
(  974)       CALL feti_creadr(malxm,malxmax,nxm,nnih*nmaxd,madwj,'dwj') 






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 18

(  975)       CALL feti_creadr(malxm,malxmax,nxm,nmaxd,madwwj,'dwwj') 
(  976)       CALL feti_creadr(malxm,malxmax,nxm,nmaxd,magamm,'gamm') 
(  977)       CALL feti_creadr(malxm,malxmax,nxm,max0(nmaxd,nitmaxete),mawork,'work')
(  978)       mjj0 = 0 
(  979)       numit0ete = 0 
(  980) 
(  981)    END SUBROUTINE fetsch
(  982) 
(  983) #else
(  984)    SUBROUTINE fetstr                 ! Empty routine
(  985)    END SUBROUTINE fetstr
(  986)    SUBROUTINE fetmat                 ! Empty routine
(  987)    END SUBROUTINE fetmat
(  988)    SUBROUTINE fetsch                 ! Empty routine
(  989)    END SUBROUTINE fetsch
(  990) #endif
(  991) 
(  992)    !!======================================================================
(  993) END MODULE solmat













































PGF90 (Version     10.5)          09/23/2010  09:49:42      page 19

(    1) # 1 "solmat.F90"
(    1) MODULE solmat
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  solmat  ***
(    4)    !! solver       : construction of the matrix 
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   sol_mat       : Construction of the matrix of used by the elliptic solvers
(    9)    !!   fetsch        :
(   10)    !!   fetmat        :
(   11)    !!   fetstr        :
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and active tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE sol_oce         ! ocean solver
(   17)    USE phycst          ! physical constants
(   18)    USE obc_oce         ! ocean open boundary conditions
(   19)    USE lbclnk          ! lateral boudary conditions
(   20)    USE lib_mpp         ! distributed memory computing
(   21)    USE in_out_manager  ! I/O manager
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC sol_mat     ! routine called by inisol.F90
(   28)    !!----------------------------------------------------------------------
(   29)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solmat.F90,v 1.11 2006/03/20 17:27:14 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36)    SUBROUTINE sol_mat( kt )
(   37)       !!----------------------------------------------------------------------
(   38)       !!                  ***  ROUTINE sol_mat  ***
(   39)       !!
(   40)       !! ** Purpose :   Construction of the matrix of used by the elliptic 
(   41)       !!      solvers (either sor, pcg or feti methods).
(   42)       !!
(   43)       !! ** Method  :   The matrix depends on the type of free surface:
(   44)       !!       * lk_dynspg_rl=T: rigid lid formulation
(   45)       !!      The matrix is built for the barotropic stream function system.
(   46)       !!      a diagonal preconditioning matrix is also defined.
(   47)       !!       * lk_dynspg_flt=T: free surface formulation
(   48)       !!      The matrix is built for the divergence of the transport system
(   49)       !!      a diagonal preconditioning matrix is also defined.
(   50)       !!        Note that for feti solver (nsolv=3) a specific initialization 
(   51)       !!      is required (call to fetstr.F) for memory allocation and inter-
(   52)       !!      face definition.
(   53)       !! 
(   54)       !! ** Action  : - gcp    : extra-diagonal elements of the matrix
(   55)       !!              - gcdmat : preconditioning matrix (diagonal elements)
(   56)       !!              - gcdprc : inverse of the preconditioning matrix
(   57)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 20

(   58)       !! History :
(   59)       !!   1.0  !  88-04  (G. Madec)  Original code
(   60)       !!        !  91-11  (G. Madec)
(   61)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   62)       !!        !  93-06  (M. Guyon)  suppress pointers
(   63)       !!        !  96-05  (G. Madec)  merge sor and pcg formulations
(   64)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(   65)       !!        !  98-02  (M. Guyon)  FETI method
(   66)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   67)       !!        !  02-11  (C. Talandier, A-M. Treguier) Free surface & Open boundaries
(   68)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   69)       !!----------------------------------------------------------------------
(   70)       !! * Arguments
(   71)       INTEGER, INTENT(in) :: kt
(   72) 
(   73)       !! * Local declarations
(   74)       INTEGER ::   ji, jj                    ! dummy loop indices
(   75)       INTEGER ::   ii, ij, iiend, ijend      ! temporary integers
(   76)       REAL(wp) ::   zcoefs, zcoefw, zcoefe, zcoefn  ! temporary scalars
(   77)       REAL(wp) ::   z2dt, zcoef
(   78)       !!----------------------------------------------------------------------
(   79) 
(   80)       ! FETI method ( nsolv = 3)
(   81)       ! memory allocation and interface definition for the solver
(   82) 
(   83)       IF( nsolv == 3 )   CALL fetstr
(   84) 
(   85)       
(   86)       ! 1. Construction of the matrix
(   87)       ! -----------------------------
(   88)       
(   89)       ! initialize to zero
(   90)       zcoef = 0.e0
(   91)       gcp(:,:,1) = 0.e0
(   92)       gcp(:,:,2) = 0.e0
(   93)       gcp(:,:,3) = 0.e0
(   94)       gcp(:,:,4) = 0.e0
(   95)       
(   96)       gcdprc(:,:) = 0.e0
(   97)       gcdmat(:,:) = 0.e0
(   98)       
(   99)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  100)          z2dt = rdt
(  101)       ELSE
(  102)          z2dt = 2. * rdt
(  103)       ENDIF
(  104) 
(  105) # 173
(  173) !!cr      ELSE
(  174) 
(  175)       !   defined the coefficients for bsf elliptic system
(  176)       
(  177)       DO jj = 2, jpjm1
(  178)          DO ji = 2, jpim1
(  179)             zcoefs = -hur(ji  ,jj  ) * e1u(ji  ,jj  ) / e2u(ji  ,jj  ) * bmask(ji,jj)   ! south coefficient
(  180)             zcoefw = -hvr(ji  ,jj  ) * e2v(ji  ,jj  ) / e1v(ji  ,jj  ) * bmask(ji,jj)   ! west coefficient
(  181)             zcoefe = -hvr(ji+1,jj  ) * e2v(ji+1,jj  ) / e1v(ji+1,jj  ) * bmask(ji,jj)   ! east coefficient
(  182)             zcoefn = -hur(ji  ,jj+1) * e1u(ji  ,jj+1) / e2u(ji  ,jj+1) * bmask(ji,jj)   ! north coefficient






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 21

(  183)             gcp(ji,jj,1) = zcoefs
(  184)             gcp(ji,jj,2) = zcoefw
(  185)             gcp(ji,jj,3) = zcoefe
(  186)             gcp(ji,jj,4) = zcoefn
(  187)             gcdmat(ji,jj) = -zcoefs -zcoefw -zcoefe -zcoefn                             ! diagonal coefficient
(  188)          END DO
(  189)       END DO
(  190)       
(  191) !!cr  ENDIF
(  192) 
(  193) # 245
(  245)       ! 2. Boundary conditions 
(  246)       ! ----------------------
(  247)       
(  248)       ! Cyclic east-west boundary conditions
(  249)       !     ji=2 is the column east of ji=jpim1 and reciprocally,
(  250)       !     ji=jpim1 is the column west of ji=2
(  251)       !     all the coef are already set to zero as bmask is initialized to
(  252)       !     zero for ji=1 and ji=jpj in dommsk.
(  253)       
(  254)       ! Symetrical conditions
(  255)       ! free surface: no specific action
(  256)       ! bsf system: n-s gradient of bsf = 0 along j=2 (perhaps a bug !!!!!!)
(  257)       ! the diagonal coefficient of the southern grid points must be modify to
(  258)       ! account for the existence of the south symmetric bassin.
(  259)       
(  260) !!cr      IF( .NOT.lk_dynspg_flt ) THEN   !bug missing lk_dynspg_flt_atsk
(  261) # 262
(  262)       IF( nperio == 2 ) THEN
(  263)          DO ji = 1, jpi
(  264)             IF( bmask(ji,2) /= 0.e0 ) THEN
(  265)                zcoefs = - hur(ji,2)*e1u(ji,2)/e2u(ji,2)
(  266)                gcdmat(ji,2) = gcdmat(ji,2) - zcoefs
(  267)             ENDIF
(  268)          END DO
(  269)       ENDIF
(  270) !!cr      ENDIF
(  271) # 272
(  272)       
(  273)       ! North fold boundary condition
(  274)       !     all the coef are already set to zero as bmask is initialized to
(  275)       !     zero on duplicated lignes and portion of lignes
(  276)       
(  277)       ! 3. Preconditioned matrix
(  278)       ! ------------------------
(  279)       
(  280)       IF( nsolv /= 3 ) THEN
(  281)          
(  282)          ! SOR and PCG solvers
(  283)          DO jj = 1, jpj
(  284)             DO ji = 1, jpi
(  285)                IF( bmask(ji,jj) /= 0.e0 )   gcdprc(ji,jj) = 1.e0 / gcdmat(ji,jj)
(  286)             END DO
(  287)          END DO
(  288)          
(  289)          gcp(:,:,1) = gcp(:,:,1) * gcdprc(:,:)
(  290)          gcp(:,:,2) = gcp(:,:,2) * gcdprc(:,:)
(  291)          gcp(:,:,3) = gcp(:,:,3) * gcdprc(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 22

(  292)          gcp(:,:,4) = gcp(:,:,4) * gcdprc(:,:)
(  293)          IF( ( nsolv == 2 ) .OR. ( nsolv == 4 ) )  gccd(:,:) = sor * gcp(:,:,2)
(  294) 
(  295)          IF( nsolv == 4 ) THEN
(  296)             CALL lbc_lnk_e( gcp   (:,:,1), c_solver_pt, 1. )   ! lateral boundary conditions
(  297)             CALL lbc_lnk_e( gcp   (:,:,2), c_solver_pt, 1. )   ! lateral boundary conditions
(  298)             CALL lbc_lnk_e( gcp   (:,:,3), c_solver_pt, 1. )   ! lateral boundary conditions
(  299)             CALL lbc_lnk_e( gcp   (:,:,4), c_solver_pt, 1. )   ! lateral boundary conditions
(  300)             CALL lbc_lnk_e( gcdprc(:,:)  , c_solver_pt, 1. )   ! lateral boundary conditions
(  301)             CALL lbc_lnk_e( gcdmat(:,:)  , c_solver_pt, 1. )   ! lateral boundary conditions         
(  302)             IF( npolj /= 0 ) CALL sol_exd( gcp , c_solver_pt ) ! switch northernelements
(  303)          END IF
(  304) 
(  305)       ELSE
(  306)          
(  307)          ! FETI method
(  308)          ! if feti solver : gcdprc is a mask for the non-overlapping
(  309)          !   data structuring
(  310)          
(  311)          DO jj = 1, jpj
(  312)             DO ji = 1, jpi
(  313)                IF( bmask(ji,jj) /= 0.e0 ) THEN
(  314)                   gcdprc(ji,jj) = 1.e0
(  315)                ELSE
(  316)                   gcdprc(ji,jj) = 0.e0
(  317)                ENDIF
(  318)             END DO
(  319)          END DO
(  320)          
(  321)          ! so "common" line & "common" column have to be !=0 except on global
(  322)          !   domain boundaries
(  323)          ! pbs with nbondi if nperio != 2 ?
(  324)          !   ii = nldi-1
(  325)          ! pb with nldi value if jperio==1 : nbondi modifyed at the end
(  326)          !   of inimpp.F => pb
(  327)          ! pb with periodicity conditions : iiend, ijend
(  328)          
(  329)          ijend = nlej
(  330)          iiend = nlei
(  331)          IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 )   iiend = nlci - jpreci
(  332)          ii = jpreci
(  333)          
(  334)          ! case number 1
(  335)          
(  336)          IF( nbondi /= -1 .AND. nbondi /= 2 ) THEN
(  337)             DO jj = 1, ijend
(  338)                IF( fmask(ii,jj,1) == 1. ) gcdprc(ii,jj) = 1.
(  339)             END DO
(  340)          ENDIF
(  341)          
(  342)          ! case number 2
(  343)          
(  344)          IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
(  345)             DO jj = 1, ijend
(  346)                IF( fmask(ii,jj,1) == 1. ) gcdprc(ii,jj) = 1.
(  347)             END DO
(  348)          ENDIF
(  349)          






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 23

(  350)          !      ij = nldj-1
(  351)          ! pb with nldi value if jperio==1 : nbondi modifyed at the end
(  352)          !   of inimpp.F => pb, here homogeneisation...
(  353)          
(  354)          ij = jprecj
(  355)          IF( nbondj /= -1 .AND. nbondj /= 2 ) THEN
(  356)             DO ji = 1, iiend
(  357)                IF( fmask(ji,ij,1) == 1. ) gcdprc(ji,ij) = 1.
(  358)             END DO
(  359)          ENDIF
(  360)       ENDIF
(  361)       
(  362)       
(  363)       ! 4. Initialization the arrays used in pcg
(  364)       ! ----------------------------------------
(  365)       gcx  (:,:) = 0.e0
(  366)       gcxb (:,:) = 0.e0
(  367)       gcb  (:,:) = 0.e0
(  368)       gcr  (:,:) = 0.e0
(  369)       gcdes(:,:) = 0.e0
(  370)       gccd (:,:) = 0.e0
(  371)       
(  372)       ! FETI method
(  373)       IF( nsolv == 3 ) THEN
(  374)          CALL fetmat       ! Matrix treatment : Neumann condition, inverse computation
(  375)          CALL fetsch       ! data framework for the Schur Dual solver
(  376)       ENDIF
(  377)       
(  378)    END SUBROUTINE sol_mat
(  379) 
(  380) 
(  381)    SUBROUTINE sol_exd( pt3d, cd_type )
(  382)       !!----------------------------------------------------------------------
(  383)       !!                  ***  routine sol_exd  ***
(  384)       !!                  
(  385)       !! ** Purpose :   Reorder gcb coefficient on the extra outer  halo 
(  386)       !!                at north fold in case of T or F pivot
(  387)       !!
(  388)       !! ** Method  :   Perform a circular permutation of the coefficients on 
(  389)       !!                the total area strictly above the pivot point,
(  390)       !!                and on the semi-row of the pivot point   
(  391)       !!                
(  392)       !! History :
(  393)       !!   9.0  !  05-09  (R. Benshila)  original routine
(  394)       !!----------------------------------------------------------------------
(  395)       !! * Arguments
(  396)       CHARACTER(len=1) , INTENT( in ) ::   &
(  397)          cd_type       ! define the nature of pt2d array grid-points
(  398)          !             !  = T , U , V , F , W 
(  399)          !             !  = S : T-point, north fold treatment
(  400)          !             !  = G : F-point, north fold treatment
(  401)          !             !  = I : sea-ice velocity at F-point with index shift
(  402)       REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4), INTENT( inout ) ::   &
(  403)          pt3d          ! 2D array on which the boundary condition is applied
(  404) 
(  405)       !! * Local variables
(  406)       INTEGER  ::   ji, jk      ! dummy loop indices
(  407)       INTEGER  ::   iloc                ! temporary integers






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 24

(  408)       REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4) ::   &
(  409)          ztab          ! 2D array on which the boundary condition is applied
(  410)       !!----------------------------------------------------------------------
(  411) 
(  412)       ztab = pt3d
(  413) 
(  414)       ! north fold treatment
(  415)       ! -----------------------
(  416)   
(  417)       SELECT CASE ( npolj )
(  418)          
(  419)          CASE ( 3 , 4 )   !  T pivot
(  420)          iloc = jpiglo/2 +1 
(  421)             
(  422)             SELECT CASE ( cd_type )
(  423)   
(  424)             CASE ( 'T', 'S', 'U', 'W' )
(  425)                DO jk =1, 4
(  426)                   DO ji = 1-jpr2di, nlci+jpr2di
(  427)                      pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
(  428)                   ENDDO
(  429)                ENDDO
(  430) 
(  431)               DO jk =1, 4
(  432)                   DO ji = nlci+jpr2di, 1-jpr2di,  -1
(  433)                      IF( ( ji .LT. mi0(iloc) .AND. mi0(iloc) /= 1 ) &
(  434)                        & .OR. ( mi0(iloc) == jpi+1 ) ) EXIT
(  435)                      pt3d(ji,nlcj-1,jk) = ztab(ji,nlcj-1,jk+3-2*MOD(jk+3,4))
(  436)                   ENDDO
(  437)                ENDDO
(  438) 
(  439)             CASE ( 'F' ,'G' , 'I', 'V' )
(  440)                DO jk =1, 4
(  441)                   DO ji = 1-jpr2di, nlci+jpr2di
(  442)                      pt3d(ji,nlcj-1:nlcj+jpr2dj,jk) = ztab(ji,nlcj-1:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
(  443)                   ENDDO
(  444)                ENDDO
(  445) 
(  446)             END SELECT   ! cd_type
(  447)   
(  448)          CASE ( 5 , 6 )                 ! F pivot
(  449)           iloc=jpiglo/2
(  450) 
(  451)             SELECT CASE (cd_type )
(  452) 
(  453)             CASE ( 'T'  ,'S', 'U', 'W')
(  454)                DO jk =1, 4
(  455)                   DO ji = 1-jpr2di, nlci+jpr2di
(  456)                      pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
(  457)                   ENDDO
(  458)                ENDDO
(  459) 
(  460)             CASE ( 'F' ,'G' , 'I', 'V' )
(  461)                DO jk =1, 4
(  462)                   DO ji = 1-jpr2di, nlci+jpr2di
(  463)                      pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
(  464)                   ENDDO
(  465)                ENDDO






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 25

(  466)                DO jk =1, 4
(  467)                   DO ji = nlci+jpr2di, 1-jpr2di,  -1
(  468)                     IF ( ( ji .LT. mi0(iloc) .AND. mi0(iloc) /= 1 ) &
(  469)                        & .OR. ( mi0(iloc) == jpi+1 ) ) EXIT
(  470)                     pt3d(ji,nlcj-1,jk) = ztab(ji,nlcj-1,jk+3-2*MOD(jk+3,4))
(  471)                   ENDDO
(  472)                ENDDO
(  473) 
(  474)             END SELECT   ! cd_type
(  475) 
(  476)          END SELECT   ! npolj
(  477)   
(  478)    END SUBROUTINE sol_exd
(  479) 
(  480) # 984
(  984)    SUBROUTINE fetstr                 ! Empty routine
(  985)    END SUBROUTINE fetstr
(  986)    SUBROUTINE fetmat                 ! Empty routine
(  987)    END SUBROUTINE fetmat
(  988)    SUBROUTINE fetsch                 ! Empty routine
(  989)    END SUBROUTINE fetsch
(  990) 
(  991) # 992
(  992)    !!======================================================================
(  993) END MODULE solmat
