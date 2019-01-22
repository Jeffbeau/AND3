


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

Filename: solsor_e.F90

(    1) MODULE solsor_e
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
(   23)    PUBLIC sol_sor_e              ! ???
(   24) 
(   25)    !!----------------------------------------------------------------------
(   26)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solsor_e.F90,v 1.2 2005/12/21 10:46:42 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32)       
(   33)    SUBROUTINE sol_sor_e( kindic )
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
(   45)       !!       This routine provides a MPI optimization to the existing solsor
(   46)       !!     by reducing the number of call to lbc.






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 2

(   47)       !! 
(   48)       !! ** Method  :   Successive-over-relaxation method using the red-black 
(   49)       !!      technique. The former technique used was not compatible with 
(   50)       !!      the north-fold boundary condition used in orca configurations.
(   51)       !!      Compared to the classical sol_sor, this routine provides a 
(   52)       !!      mpp optimization by reducing the number of calls to lnc_lnk
(   53)       !!      The solution is computed on a larger area and the boudary
(   54)       !!      conditions only when the inside domain is reached.
(   55)       !! 
(   56)       !! References :
(   57)       !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
(   58)       !!      Beare and Stevens 1997 Ann. Geophysicae 15, 1369-1377
(   59)       !!
(   60)       !! History :
(   61)       !!        !  90-10  (G. Madec)  Original code
(   62)       !!        !  91-11  (G. Madec)
(   63)       !!   7.1  !  93-04  (G. Madec)  time filter
(   64)       !!        !  96-05  (G. Madec)  merge sor and pcg formulations
(   65)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(   66)       !!   9.0  !  03-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
(   67)       !!   9.0  !  05-09  (R. Benshila, G. Madec)  MPI optimization
(   68)       !!----------------------------------------------------------------------
(   69)       !! * Arguments
(   70)       INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
(   71)       !                                      ! gence is not reached: the model is
(   72)       !                                      ! stopped in step
(   73)       !                                      ! set to zero before the call of solsor
(   74)       !! * Local declarations
(   75)       INTEGER  ::   ji, jj, jn               ! dummy loop indices
(   76)       INTEGER  ::   ishift, icount
(   77)       REAL(wp) ::   ztmp, zres, zres2
(   78) 
(   79)       INTEGER  ::   ijmppodd, ijmppeven
(   80)       INTEGER  ::   ijpr2d
(   81)       !!----------------------------------------------------------------------
(   82)       
(   83)       ijmppeven = MOD(nimpp+njmpp+jpr2di+jpr2dj,2)
(   84)       ijmppodd  = MOD(nimpp+njmpp+jpr2di+jpr2dj+1,2)
(   85)       ijpr2d = MAX(jpr2di,jpr2dj)
(   86)       icount = 0
(   87)       !                                                       ! ==============
(   88)       DO jn = 1, nmax                                         ! Iterative loop 
(   89)          !                                                    ! ==============
(   90) 
(   91)          ! applied the lateral boundary conditions
(   92)          IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e( gcx, c_solver_pt, 1. )   
(   93)         
(   94)          ! Residus
(   95)          ! -------
(   96) 
(   97)          ! Guess black update
(   98)          DO jj = 2-jpr2dj, nlcj-1+jpr2dj
(   99)             ishift = MOD( jj-ijmppodd-jpr2dj, 2 )
(  100)             DO ji = 2-jpr2di+ishift, nlci-1+jpr2di, 2
(  101)                ztmp =                  gcb(ji  ,jj  )   &
(  102)                   &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(  103)                   &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(  104)                   &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 3

(  105)                   &   - gcp(ji,jj,4) * gcx(ji  ,jj+1)
(  106)                ! Estimate of the residual
(  107)                zres = ztmp - gcx(ji,jj)
(  108)                gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
(  109)                ! Guess update
(  110)                gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
(  111)             END DO
(  112)          END DO
(  113)          icount = icount + 1 
(  114)  
(  115)          ! applied the lateral boundary conditions
(  116)          IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e( gcx, c_solver_pt, 1. )  
(  117) 
(  118)          ! Guess red update
(  119)          DO jj = 2-jpr2dj, nlcj-1+jpr2dj
(  120)             ishift = MOD( jj-ijmppeven-jpr2dj, 2 )
(  121)             DO ji = 2-jpr2di+ishift, nlci-1+jpr2di, 2
(  122)                ztmp =                  gcb(ji  ,jj  )   &
(  123)                   &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(  124)                   &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(  125)                   &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(  126)                   &   - gcp(ji,jj,4) * gcx(ji  ,jj+1) 
(  127)                ! Estimate of the residual
(  128)                zres = ztmp - gcx(ji,jj)
(  129)                gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
(  130)                ! Guess update
(  131)                gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
(  132)             END DO
(  133)          END DO
(  134)          icount = icount + 1
(  135) 
(  136)          ! test of convergence
(  137)          IF ( jn > nmin .AND. MOD( jn-nmin, nmod ) == 0 ) then
(  138) 
(  139)             SELECT CASE ( nsol_arp )
(  140)             CASE ( 0 )                 ! absolute precision (maximum value of the residual)
(  141)                zres2 = MAXVAL( gcr(2:nlci-1,2:nlcj-1) )
(  142)                IF( lk_mpp )   CALL mpp_max( zres2 )   ! max over the global domain
(  143)                ! test of convergence
(  144)                IF( zres2 < resmax .OR. jn == nmax ) THEN
(  145)                   res = SQRT( zres2 )
(  146)                   niter = jn
(  147)                   ncut = 999
(  148)                ENDIF
(  149)             CASE ( 1 )                 ! relative precision
(  150)                rnorme = SUM( gcr(2:nlci-1,2:nlcj-1) )
(  151)                IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  152)                ! test of convergence
(  153)                IF( rnorme < epsr .OR. jn == nmax ) THEN
(  154)                   res = SQRT( rnorme )
(  155)                   niter = jn
(  156)                   ncut = 999
(  157)                ENDIF
(  158)             END SELECT
(  159)          
(  160)          !****
(  161) 9300     FORMAT('          niter :',i4,' res :',e20.10,' b :',e20.10)
(  162)          !****






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 4

(  163)          
(  164)          ENDIF
(  165)          ! indicator of non-convergence or explosion
(  166)          IF( jn == nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
(  167)          IF( ncut == 999 ) GOTO 999
(  168)          
(  169)          !                                                 ! =====================
(  170)       END DO                                               ! END of iterative loop
(  171)       !                                                    ! =====================
(  172)       
(  173) 999   CONTINUE
(  174)       
(  175)       
(  176)       !  Output in gcx
(  177)       !  -------------
(  178) 
(  179)       CALL lbc_lnk_e( gcx, c_solver_pt, 1. )    ! boundary conditions
(  180) 
(  181)       
(  182)    END SUBROUTINE sol_sor_e
(  183) 
(  184)    !!=====================================================================
(  185) END MODULE solsor_e









































PGF90 (Version     10.5)          09/23/2010  09:49:42      page 5

(    1) # 1 "solsor_e.F90"
(    1) MODULE solsor_e
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
(   23)    PUBLIC sol_sor_e              ! ???
(   24) 
(   25)    !!----------------------------------------------------------------------
(   26)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solsor_e.F90,v 1.2 2005/12/21 10:46:42 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32)       
(   33)    SUBROUTINE sol_sor_e( kindic )
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
(   45)       !!       This routine provides a MPI optimization to the existing solsor
(   46)       !!     by reducing the number of call to lbc.
(   47)       !! 
(   48)       !! ** Method  :   Successive-over-relaxation method using the red-black 
(   49)       !!      technique. The former technique used was not compatible with 
(   50)       !!      the north-fold boundary condition used in orca configurations.
(   51)       !!      Compared to the classical sol_sor, this routine provides a 
(   52)       !!      mpp optimization by reducing the number of calls to lnc_lnk
(   53)       !!      The solution is computed on a larger area and the boudary
(   54)       !!      conditions only when the inside domain is reached.
(   55)       !! 
(   56)       !! References :
(   57)       !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 6

(   58)       !!      Beare and Stevens 1997 Ann. Geophysicae 15, 1369-1377
(   59)       !!
(   60)       !! History :
(   61)       !!        !  90-10  (G. Madec)  Original code
(   62)       !!        !  91-11  (G. Madec)
(   63)       !!   7.1  !  93-04  (G. Madec)  time filter
(   64)       !!        !  96-05  (G. Madec)  merge sor and pcg formulations
(   65)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(   66)       !!   9.0  !  03-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
(   67)       !!   9.0  !  05-09  (R. Benshila, G. Madec)  MPI optimization
(   68)       !!----------------------------------------------------------------------
(   69)       !! * Arguments
(   70)       INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
(   71)       !                                      ! gence is not reached: the model is
(   72)       !                                      ! stopped in step
(   73)       !                                      ! set to zero before the call of solsor
(   74)       !! * Local declarations
(   75)       INTEGER  ::   ji, jj, jn               ! dummy loop indices
(   76)       INTEGER  ::   ishift, icount
(   77)       REAL(wp) ::   ztmp, zres, zres2
(   78) 
(   79)       INTEGER  ::   ijmppodd, ijmppeven
(   80)       INTEGER  ::   ijpr2d
(   81)       !!----------------------------------------------------------------------
(   82)       
(   83)       ijmppeven = MOD(nimpp+njmpp+jpr2di+jpr2dj,2)
(   84)       ijmppodd  = MOD(nimpp+njmpp+jpr2di+jpr2dj+1,2)
(   85)       ijpr2d = MAX(jpr2di,jpr2dj)
(   86)       icount = 0
(   87)       !                                                       ! ==============
(   88)       DO jn = 1, nmax                                         ! Iterative loop 
(   89)          !                                                    ! ==============
(   90) 
(   91)          ! applied the lateral boundary conditions
(   92)          IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e( gcx, c_solver_pt, 1. )   
(   93)         
(   94)          ! Residus
(   95)          ! -------
(   96) 
(   97)          ! Guess black update
(   98)          DO jj = 2-jpr2dj, nlcj-1+jpr2dj
(   99)             ishift = MOD( jj-ijmppodd-jpr2dj, 2 )
(  100)             DO ji = 2-jpr2di+ishift, nlci-1+jpr2di, 2
(  101)                ztmp =                  gcb(ji  ,jj  )   &
(  102)                   &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(  103)                   &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(  104)                   &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(  105)                   &   - gcp(ji,jj,4) * gcx(ji  ,jj+1)
(  106)                ! Estimate of the residual
(  107)                zres = ztmp - gcx(ji,jj)
(  108)                gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
(  109)                ! Guess update
(  110)                gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
(  111)             END DO
(  112)          END DO
(  113)          icount = icount + 1 
(  114)  
(  115)          ! applied the lateral boundary conditions






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 7

(  116)          IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e( gcx, c_solver_pt, 1. )  
(  117) 
(  118)          ! Guess red update
(  119)          DO jj = 2-jpr2dj, nlcj-1+jpr2dj
(  120)             ishift = MOD( jj-ijmppeven-jpr2dj, 2 )
(  121)             DO ji = 2-jpr2di+ishift, nlci-1+jpr2di, 2
(  122)                ztmp =                  gcb(ji  ,jj  )   &
(  123)                   &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(  124)                   &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(  125)                   &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(  126)                   &   - gcp(ji,jj,4) * gcx(ji  ,jj+1) 
(  127)                ! Estimate of the residual
(  128)                zres = ztmp - gcx(ji,jj)
(  129)                gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
(  130)                ! Guess update
(  131)                gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
(  132)             END DO
(  133)          END DO
(  134)          icount = icount + 1
(  135) 
(  136)          ! test of convergence
(  137)          IF ( jn > nmin .AND. MOD( jn-nmin, nmod ) == 0 ) then
(  138) 
(  139)             SELECT CASE ( nsol_arp )
(  140)             CASE ( 0 )                 ! absolute precision (maximum value of the residual)
(  141)                zres2 = MAXVAL( gcr(2:nlci-1,2:nlcj-1) )
(  142)                IF( lk_mpp )   CALL mpp_max( zres2 )   ! max over the global domain
(  143)                ! test of convergence
(  144)                IF( zres2 < resmax .OR. jn == nmax ) THEN
(  145)                   res = SQRT( zres2 )
(  146)                   niter = jn
(  147)                   ncut = 999
(  148)                ENDIF
(  149)             CASE ( 1 )                 ! relative precision
(  150)                rnorme = SUM( gcr(2:nlci-1,2:nlcj-1) )
(  151)                IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  152)                ! test of convergence
(  153)                IF( rnorme < epsr .OR. jn == nmax ) THEN
(  154)                   res = SQRT( rnorme )
(  155)                   niter = jn
(  156)                   ncut = 999
(  157)                ENDIF
(  158)             END SELECT
(  159)          
(  160)          !****
(  161) 9300     FORMAT('          niter :',i4,' res :',e20.10,' b :',e20.10)
(  162)          !****
(  163)          
(  164)          ENDIF
(  165)          ! indicator of non-convergence or explosion
(  166)          IF( jn == nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
(  167)          IF( ncut == 999 ) GOTO 999
(  168)          
(  169)          !                                                 ! =====================
(  170)       END DO                                               ! END of iterative loop
(  171)       !                                                    ! =====================
(  172)       
(  173) 999   CONTINUE






PGF90 (Version     10.5)          09/23/2010  09:49:42      page 8

(  174)       
(  175)       
(  176)       !  Output in gcx
(  177)       !  -------------
(  178) 
(  179)       CALL lbc_lnk_e( gcx, c_solver_pt, 1. )    ! boundary conditions
(  180) 
(  181)       
(  182)    END SUBROUTINE sol_sor_e
(  183) 
(  184)    !!=====================================================================
(  185) END MODULE solsor_e
