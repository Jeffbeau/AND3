


PGF90 (Version     10.5)          09/23/2010  09:49:26      page 1

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

Filename: solsor.F90

(    1) MODULE solsor
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  solsor  ***
(    4)    !! Ocean solver :  Successive Over-Relaxation solver
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   sol_sor     : Red-Black Successive Over-Relaxation solver
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE zdf_oce         ! ocean vertical physics variables
(   14)    USE sol_oce         ! solver variables
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE lib_mpp         ! distributed memory computing
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Routine accessibility
(   23)    PUBLIC sol_sor              ! ???
(   24) 
(   25)    !!----------------------------------------------------------------------
(   26)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solsor.F90,v 1.8 2005/12/21 10:46:42 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32)       
(   33)    SUBROUTINE sol_sor( kindic )
(   34)       !!----------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE sol_sor  ***
(   36)       !!                 
(   37)       !! ** Purpose :   Solve the ellipic equation for the barotropic stream 
(   38)       !!      function system (lk_dynspg_rl=T) or the transport divergence 
(   39)       !!      system (lk_dynspg_flt=T) using a red-black successive-over-
(   40)       !!      relaxation method.
(   41)       !!       In the former case, the barotropic stream function trend has a
(   42)       !!     zero boundary condition along all coastlines (i.e. continent
(   43)       !!     as well as islands) while in the latter the boundary condition
(   44)       !!     specification is not required.
(   45)       !!
(   46)       !! ** Method  :   Successive-over-relaxation method using the red-black 






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 2

(   47)       !!      technique. The former technique used was not compatible with 
(   48)       !!      the north-fold boundary condition used in orca configurations.
(   49)       !!
(   50)       !! References :
(   51)       !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
(   52)       !!
(   53)       !! History :
(   54)       !!        !  90-10  (G. Madec)  Original code
(   55)       !!        !  91-11  (G. Madec)
(   56)       !!   7.1  !  93-04  (G. Madec)  time filter
(   57)       !!        !  96-05  (G. Madec)  merge sor and pcg formulations
(   58)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(   59)       !!   9.0  !  03-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
(   60)       !!----------------------------------------------------------------------
(   61)       !! * Arguments
(   62)       INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
(   63)       !                                      ! gence is not reached: the model is
(   64)       !                                      ! stopped in step
(   65)       !                                      ! set to zero before the call of solsor
(   66)       !! * Local declarations
(   67)       INTEGER  ::   ji, jj, jn               ! dummy loop indices
(   68)       INTEGER  ::   ishift
(   69)       REAL(wp) ::   ztmp, zres, zres2
(   70) 
(   71)       INTEGER  ::   ijmppodd, ijmppeven
(   72)       !!----------------------------------------------------------------------
(   73)       
(   74)       ijmppeven = MOD(nimpp+njmpp  ,2)
(   75)       ijmppodd  = MOD(nimpp+njmpp+1,2)
(   76)       !                                                       ! ==============
(   77)       DO jn = 1, nmax                                         ! Iterative loop 
(   78)          !                                                    ! ==============
(   79) 
(   80)          CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! applied the lateral boundary conditions
(   81)          
(   82)          ! Residus
(   83)          ! -------
(   84) 
(   85)          ! Guess black update
(   86)          DO jj = 2, jpjm1
(   87)             ishift = MOD( jj-ijmppodd, 2 )
(   88)             DO ji = 2+ishift, jpim1, 2
(   89)                ztmp =                  gcb(ji  ,jj  )   &
(   90)                   &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(   91)                   &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(   92)                   &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(   93)                   &   - gcp(ji,jj,4) * gcx(ji  ,jj+1)
(   94)                ! Estimate of the residual
(   95)                zres = ztmp - gcx(ji,jj)
(   96)                gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
(   97)                ! Guess update
(   98)                gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
(   99)             END DO
(  100)          END DO
(  101) 
(  102)          CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! applied the lateral boubary conditions
(  103) 
(  104)          ! Guess red update






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 3

(  105)          DO jj = 2, jpjm1
(  106)             ishift = MOD( jj-ijmppeven, 2 )
(  107)             DO ji = 2+ishift, jpim1, 2
(  108)                ztmp =                  gcb(ji  ,jj  )   &
(  109)                   &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(  110)                   &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(  111)                   &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(  112)                   &   - gcp(ji,jj,4) * gcx(ji  ,jj+1) 
(  113)                ! Estimate of the residual
(  114)                zres = ztmp - gcx(ji,jj)
(  115)                gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
(  116)                ! Guess update
(  117)                gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
(  118)             END DO
(  119)          END DO
(  120) 
(  121)          ! test of convergence
(  122)          IF ( jn > nmin .AND. MOD( jn-nmin, nmod ) == 0 ) then
(  123) 
(  124)             SELECT CASE ( nsol_arp )
(  125)             CASE ( 0 )                 ! absolute precision (maximum value of the residual)
(  126)                zres2 = MAXVAL( gcr(2:jpim1,2:jpjm1) )
(  127)                IF( lk_mpp )   CALL mpp_max( zres2 )   ! max over the global domain
(  128)                ! test of convergence
(  129)                IF( zres2 < resmax .OR. jn == nmax ) THEN
(  130)                   res = SQRT( zres2 )
(  131)                   niter = jn
(  132)                   ncut = 999
(  133)                ENDIF
(  134)             CASE ( 1 )                 ! relative precision
(  135)                rnorme = SUM( gcr(2:jpim1,2:jpjm1) )
(  136)                IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  137)                ! test of convergence
(  138)                IF( rnorme < epsr .OR. jn == nmax ) THEN
(  139)                   res = SQRT( rnorme )
(  140)                   niter = jn
(  141)                   ncut = 999
(  142)                ENDIF
(  143)             END SELECT
(  144)          
(  145)          !****
(  146) 9300     FORMAT('          niter :',i4,' res :',e20.10,' b :',e20.10)
(  147)          !****
(  148)          
(  149)          ENDIF
(  150)          ! indicator of non-convergence or explosion
(  151)          IF( jn == nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
(  152)          IF( ncut == 999 ) GOTO 999
(  153)          
(  154)          !                                                 ! =====================
(  155)       END DO                                               ! END of iterative loop
(  156)       !                                                    ! =====================
(  157)       
(  158) 999   CONTINUE
(  159)       
(  160)       
(  161)       !  Output in gcx
(  162)       !  -------------






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 4

(  163) 
(  164)       CALL lbc_lnk( gcx, c_solver_pt, 1. )    ! boundary conditions
(  165) 
(  166)       
(  167)    END SUBROUTINE sol_sor
(  168) 
(  169)    !!=====================================================================
(  170) END MODULE solsor
























































PGF90 (Version     10.5)          09/23/2010  09:49:26      page 5

(    1) # 1 "solsor.F90"
(    1) MODULE solsor
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  solsor  ***
(    4)    !! Ocean solver :  Successive Over-Relaxation solver
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   sol_sor     : Red-Black Successive Over-Relaxation solver
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE zdf_oce         ! ocean vertical physics variables
(   14)    USE sol_oce         ! solver variables
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE lib_mpp         ! distributed memory computing
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Routine accessibility
(   23)    PUBLIC sol_sor              ! ???
(   24) 
(   25)    !!----------------------------------------------------------------------
(   26)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solsor.F90,v 1.8 2005/12/21 10:46:42 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32)       
(   33)    SUBROUTINE sol_sor( kindic )
(   34)       !!----------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE sol_sor  ***
(   36)       !!                 
(   37)       !! ** Purpose :   Solve the ellipic equation for the barotropic stream 
(   38)       !!      function system (lk_dynspg_rl=T) or the transport divergence 
(   39)       !!      system (lk_dynspg_flt=T) using a red-black successive-over-
(   40)       !!      relaxation method.
(   41)       !!       In the former case, the barotropic stream function trend has a
(   42)       !!     zero boundary condition along all coastlines (i.e. continent
(   43)       !!     as well as islands) while in the latter the boundary condition
(   44)       !!     specification is not required.
(   45)       !!
(   46)       !! ** Method  :   Successive-over-relaxation method using the red-black 
(   47)       !!      technique. The former technique used was not compatible with 
(   48)       !!      the north-fold boundary condition used in orca configurations.
(   49)       !!
(   50)       !! References :
(   51)       !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
(   52)       !!
(   53)       !! History :
(   54)       !!        !  90-10  (G. Madec)  Original code
(   55)       !!        !  91-11  (G. Madec)
(   56)       !!   7.1  !  93-04  (G. Madec)  time filter
(   57)       !!        !  96-05  (G. Madec)  merge sor and pcg formulations






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 6

(   58)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(   59)       !!   9.0  !  03-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
(   60)       !!----------------------------------------------------------------------
(   61)       !! * Arguments
(   62)       INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
(   63)       !                                      ! gence is not reached: the model is
(   64)       !                                      ! stopped in step
(   65)       !                                      ! set to zero before the call of solsor
(   66)       !! * Local declarations
(   67)       INTEGER  ::   ji, jj, jn               ! dummy loop indices
(   68)       INTEGER  ::   ishift
(   69)       REAL(wp) ::   ztmp, zres, zres2
(   70) 
(   71)       INTEGER  ::   ijmppodd, ijmppeven
(   72)       !!----------------------------------------------------------------------
(   73)       
(   74)       ijmppeven = MOD(nimpp+njmpp  ,2)
(   75)       ijmppodd  = MOD(nimpp+njmpp+1,2)
(   76)       !                                                       ! ==============
(   77)       DO jn = 1, nmax                                         ! Iterative loop 
(   78)          !                                                    ! ==============
(   79) 
(   80)          CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! applied the lateral boundary conditions
(   81)          
(   82)          ! Residus
(   83)          ! -------
(   84) 
(   85)          ! Guess black update
(   86)          DO jj = 2, jpjm1
(   87)             ishift = MOD( jj-ijmppodd, 2 )
(   88)             DO ji = 2+ishift, jpim1, 2
(   89)                ztmp =                  gcb(ji  ,jj  )   &
(   90)                   &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(   91)                   &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(   92)                   &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(   93)                   &   - gcp(ji,jj,4) * gcx(ji  ,jj+1)
(   94)                ! Estimate of the residual
(   95)                zres = ztmp - gcx(ji,jj)
(   96)                gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
(   97)                ! Guess update
(   98)                gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
(   99)             END DO
(  100)          END DO
(  101) 
(  102)          CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! applied the lateral boubary conditions
(  103) 
(  104)          ! Guess red update
(  105)          DO jj = 2, jpjm1
(  106)             ishift = MOD( jj-ijmppeven, 2 )
(  107)             DO ji = 2+ishift, jpim1, 2
(  108)                ztmp =                  gcb(ji  ,jj  )   &
(  109)                   &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(  110)                   &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(  111)                   &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(  112)                   &   - gcp(ji,jj,4) * gcx(ji  ,jj+1) 
(  113)                ! Estimate of the residual
(  114)                zres = ztmp - gcx(ji,jj)
(  115)                gcr(ji,jj) = zres * gcdmat(ji,jj) * zres






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 7

(  116)                ! Guess update
(  117)                gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
(  118)             END DO
(  119)          END DO
(  120) 
(  121)          ! test of convergence
(  122)          IF ( jn > nmin .AND. MOD( jn-nmin, nmod ) == 0 ) then
(  123) 
(  124)             SELECT CASE ( nsol_arp )
(  125)             CASE ( 0 )                 ! absolute precision (maximum value of the residual)
(  126)                zres2 = MAXVAL( gcr(2:jpim1,2:jpjm1) )
(  127)                IF( lk_mpp )   CALL mpp_max( zres2 )   ! max over the global domain
(  128)                ! test of convergence
(  129)                IF( zres2 < resmax .OR. jn == nmax ) THEN
(  130)                   res = SQRT( zres2 )
(  131)                   niter = jn
(  132)                   ncut = 999
(  133)                ENDIF
(  134)             CASE ( 1 )                 ! relative precision
(  135)                rnorme = SUM( gcr(2:jpim1,2:jpjm1) )
(  136)                IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  137)                ! test of convergence
(  138)                IF( rnorme < epsr .OR. jn == nmax ) THEN
(  139)                   res = SQRT( rnorme )
(  140)                   niter = jn
(  141)                   ncut = 999
(  142)                ENDIF
(  143)             END SELECT
(  144)          
(  145)          !****
(  146) 9300     FORMAT('          niter :',i4,' res :',e20.10,' b :',e20.10)
(  147)          !****
(  148)          
(  149)          ENDIF
(  150)          ! indicator of non-convergence or explosion
(  151)          IF( jn == nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
(  152)          IF( ncut == 999 ) GOTO 999
(  153)          
(  154)          !                                                 ! =====================
(  155)       END DO                                               ! END of iterative loop
(  156)       !                                                    ! =====================
(  157)       
(  158) 999   CONTINUE
(  159)       
(  160)       
(  161)       !  Output in gcx
(  162)       !  -------------
(  163) 
(  164)       CALL lbc_lnk( gcx, c_solver_pt, 1. )    ! boundary conditions
(  165) 
(  166)       
(  167)    END SUBROUTINE sol_sor
(  168) 
(  169)    !!=====================================================================
(  170) END MODULE solsor
