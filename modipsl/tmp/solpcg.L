


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

Filename: solpcg.F90

(    1) MODULE solpcg
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  solfet
(    4)    !! Ocean solver :  preconditionned conjugate gradient solver
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   sol_pcg    : preconditionned conjugate gradient solver
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE sol_oce         ! ocean solver variables
(   14)    USE lib_mpp         ! distributed memory computing
(   15)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   16)    USE in_out_manager  ! I/O manager
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC sol_pcg              ! ???
(   23) 
(   24)    !! * Substitutions
(   25) #  include "vectopt_loop_substitute.h90"
(   26)    !!----------------------------------------------------------------------
(   27)    !!----------------------------------------------------------------------
(   28)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   29)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solpcg.F90,v 1.4 2005/12/21 10:46:42 opalod Exp $ 
(   30)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   31)    !!----------------------------------------------------------------------
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE sol_pcg( kindic )
(   35)       !!----------------------------------------------------------------------
(   36)       !!                  ***  ROUTINE sol_pcg  ***
(   37)       !!                    
(   38)       !! ** Purpose :   Solve the ellipic equation for the barotropic stream 
(   39)       !!      function system (lk_dynspg_rl=T) or the transport divergence 
(   40)       !!      system (lk_dynspg_flt=T) using a diagonal preconditionned
(   41)       !!      conjugate gradient method.
(   42)       !!      In the former case, the barotropic stream function trend has a
(   43)       !!      zero boundary condition along all coastlines (i.e. continent
(   44)       !!      as well as islands) while in the latter the boundary condition
(   45)       !!      specification is not required.
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 2

(   47)       !! ** Method  :   Diagonal preconditionned conjugate gradient method.
(   48)       !!      the algorithm is multitasked. (case of 5 points matrix)
(   49)       !!      define              pa  = q^-1 * a
(   50)       !!                        pgcb  = q^-1 * gcb
(   51)       !!                 < . ; . >_q  = ( . )^t q ( . )
(   52)       !!      where q is the preconditioning matrix = diagonal matrix of the
(   53)       !!                                              diagonal elements of a
(   54)       !!      Initialization:
(   55)       !!         x(o) = gcx
(   56)       !!         r(o) = d(o) = pgcb - pa.x(o)
(   57)       !!         rr(o)= < r(o) , r(o) >_q
(   58)       !!      Iteration n   :
(   59)       !!         z(n)   = pa.d(n)
(   60)       !!         alp(n) = rr(n) / < z(n) , d(n) >_q
(   61)       !!         x(n+1) = x(n) + alp(n) d(n)
(   62)       !!         r(n+1) = r(n) - alp(n) z(n)
(   63)       !!         rr(n+1)= < r(n+1) , r(n+1) >_q
(   64)       !!         bet(n) = rr(n+1) / rr(n)
(   65)       !!         r(n+1) = r(n+1) + bet(n+1) d(n)
(   66)       !!      Convergence test :
(   67)       !!         rr(n+1) / < gcb , gcb >_q   =< epsr
(   68)       !!
(   69)       !! ** Action : - niter  : solver number of iteration done
(   70)       !!             - res    : solver residu reached
(   71)       !!             - gcx()  : solution of the elliptic system
(   72)       !!
(   73)       !! References :
(   74)       !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
(   75)       !!
(   76)       !! History :
(   77)       !!        !  90-10  (G. Madec)  Original code
(   78)       !!        !  91-11  (G. Madec)
(   79)       !!        !  93-04  (M. Guyon)  loops and suppress pointers
(   80)       !!        !  95-09  (M. Imbard, J. Escobar)  mpp exchange 
(   81)       !!        !  96-05  (G. Madec)  merge sor and pcg formulations
(   82)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(   83)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   84)       !!----------------------------------------------------------------------
(   85)       !! * Arguments
(   86)       INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
(   87)       !                                      ! gence is not reached: the model is
(   88)       !                                      ! stopped in step
(   89)       !                                      ! set to zero before the call of solpcg
(   90) 
(   91)       !! * Local declarations
(   92)       INTEGER ::   ji, jj, jn                ! dummy loop indices
(   93)       REAL(wp) ::   zgcad                    ! temporary scalars
(   94)       !!----------------------------------------------------------------------
(   95) 
(   96)       !                                                !================
(   97)       DO jn = 1, nmax                                  ! Iterative loop
(   98)          !                                             !================
(   99) 
(  100)          IF( jn == 1 ) THEN           ! Initialization of the algorithm
(  101) 
(  102)             CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! lateral boundary condition
(  103) 
(  104)             ! gcr   = gcb-a.gcx






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 3

(  105)             ! gcdes = gsr
(  106)             DO jj = 2, jpjm1
(  107)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  108)                   zgcad = bmask(ji,jj) * ( gcb(ji,jj  ) -                gcx(ji  ,jj  )   &
(  109)                      &                                  - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(  110)                      &                                  - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(  111)                      &                                  - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(  112)                      &                                  - gcp(ji,jj,4) * gcx(ji  ,jj+1)   )
(  113)                   gcr  (ji,jj) = zgcad
(  114)                   gcdes(ji,jj) = zgcad
(  115)                END DO
(  116)             END DO
(  117)             
(  118)             rnorme = SUM(  gcr(:,:) * gcdmat(:,:) * gcr(:,:)  )
(  119)             IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  120)             rr = rnorme
(  121) 
(  122)          ENDIF
(  123)         
(  124)          !                             ! Algorithm
(  125)         
(  126)          CALL lbc_lnk( gcdes, c_solver_pt, 1. )   ! lateral boundary condition
(  127)         
(  128)          ! ... gccd = matrix . gcdes
(  129)          DO jj = 2, jpjm1
(  130)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  131)                gccd(ji,jj) = bmask(ji,jj)*( gcdes(ji,jj)   &
(  132)                   &        +gcp(ji,jj,1)*gcdes(ji,jj-1)+gcp(ji,jj,2)*gcdes(ji-1,jj)   &
(  133)                   &        +gcp(ji,jj,4)*gcdes(ji,jj+1)+gcp(ji,jj,3)*gcdes(ji+1,jj)   )
(  134)             END DO
(  135)          END DO
(  136)  
(  137)          ! alph = (gcr,gcr)/(gcdes,gccd)
(  138)          radd = SUM(  gcdes(:,:) * gcdmat(:,:) * gccd(:,:)  )
(  139)          IF( lk_mpp )   CALL mpp_sum( radd )   ! sum over the global domain
(  140)          alph = rr / radd
(  141)          
(  142)          ! gcx = gcx + alph * gcdes
(  143)          ! gcr = gcr - alph * gccd
(  144)          DO jj = 2, jpjm1
(  145)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  146)                gcx(ji,jj) = bmask(ji,jj) * ( gcx(ji,jj) + alph * gcdes(ji,jj) )
(  147)                gcr(ji,jj) = bmask(ji,jj) * ( gcr(ji,jj) - alph * gccd (ji,jj) )
(  148)             END DO
(  149)          END DO
(  150)         
(  151)          ! rnorme = (gcr,gcr)
(  152)          rnorme = SUM(  gcr(:,:) * gcdmat(:,:) * gcr(:,:)  )
(  153)          IF( lk_mpp )   CALL  mpp_sum( rnorme )   ! sum over the global domain
(  154)         
(  155)          ! test of convergence
(  156)          IF( rnorme < epsr .OR. jn == nmax ) THEN
(  157)             res = SQRT( rnorme )
(  158)             niter = jn
(  159)             ncut = 999
(  160)          ENDIF
(  161)         
(  162)          ! beta = (rk+1,rk+1)/(rk,rk)






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 4

(  163)          beta = rnorme / rr
(  164)          rr   = rnorme
(  165) 
(  166)          ! indicator of non-convergence or explosion
(  167)          IF( jn == nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
(  168)          IF( ncut == 999 ) GOTO 999
(  169) 
(  170)          ! gcdes = gcr + beta * gcdes
(  171)          DO jj = 2, jpjm1
(  172)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  173)                gcdes(ji,jj) = bmask(ji,jj)*( gcr(ji,jj) + beta * gcdes(ji,jj) )
(  174)             END DO
(  175)          END DO
(  176)         
(  177)          !                                             !================
(  178)       END DO                                           !    End Loop
(  179)       !                                                !================
(  180)      
(  181) 999   CONTINUE
(  182)      
(  183)      
(  184)       ! Output in gcx with lateral b.c. applied
(  185)       ! ---------------------------------------
(  186)      
(  187)       CALL lbc_lnk( gcx, c_solver_pt, 1. )
(  188)      
(  189)    END SUBROUTINE sol_pcg
(  190) 
(  191)    !!=====================================================================
(  192) END MODULE solpcg


































PGF90 (Version     10.5)          09/23/2010  09:49:26      page 5

(    1) # 1 "solpcg.F90"
(    1) MODULE solpcg
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  solfet
(    4)    !! Ocean solver :  preconditionned conjugate gradient solver
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   sol_pcg    : preconditionned conjugate gradient solver
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE sol_oce         ! ocean solver variables
(   14)    USE lib_mpp         ! distributed memory computing
(   15)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   16)    USE in_out_manager  ! I/O manager
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC sol_pcg              ! ???
(   23) 
(   24)    !! * Substitutions
(   25) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 26 "solpcg.F90"
(   26) # 26
(   26)    !!----------------------------------------------------------------------
(   27)    !!----------------------------------------------------------------------
(   28)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   29)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solpcg.F90,v 1.4 2005/12/21 10:46:42 opalod Exp $ 
(   30)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   31)    !!----------------------------------------------------------------------
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE sol_pcg( kindic )
(   35)       !!----------------------------------------------------------------------
(   36)       !!                  ***  ROUTINE sol_pcg  ***
(   37)       !!                    
(   38)       !! ** Purpose :   Solve the ellipic equation for the barotropic stream 
(   39)       !!      function system (lk_dynspg_rl=T) or the transport divergence 
(   40)       !!      system (lk_dynspg_flt=T) using a diagonal preconditionned
(   41)       !!      conjugate gradient method.
(   42)       !!      In the former case, the barotropic stream function trend has a
(   43)       !!      zero boundary condition along all coastlines (i.e. continent
(   44)       !!      as well as islands) while in the latter the boundary condition






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 6

(   45)       !!      specification is not required.
(   46)       !!
(   47)       !! ** Method  :   Diagonal preconditionned conjugate gradient method.
(   48)       !!      the algorithm is multitasked. (case of 5 points matrix)
(   49)       !!      define              pa  = q^-1 * a
(   50)       !!                        pgcb  = q^-1 * gcb
(   51)       !!                 < . ; . >_q  = ( . )^t q ( . )
(   52)       !!      where q is the preconditioning matrix = diagonal matrix of the
(   53)       !!                                              diagonal elements of a
(   54)       !!      Initialization:
(   55)       !!         x(o) = gcx
(   56)       !!         r(o) = d(o) = pgcb - pa.x(o)
(   57)       !!         rr(o)= < r(o) , r(o) >_q
(   58)       !!      Iteration n   :
(   59)       !!         z(n)   = pa.d(n)
(   60)       !!         alp(n) = rr(n) / < z(n) , d(n) >_q
(   61)       !!         x(n+1) = x(n) + alp(n) d(n)
(   62)       !!         r(n+1) = r(n) - alp(n) z(n)
(   63)       !!         rr(n+1)= < r(n+1) , r(n+1) >_q
(   64)       !!         bet(n) = rr(n+1) / rr(n)
(   65)       !!         r(n+1) = r(n+1) + bet(n+1) d(n)
(   66)       !!      Convergence test :
(   67)       !!         rr(n+1) / < gcb , gcb >_q   =< epsr
(   68)       !!
(   69)       !! ** Action : - niter  : solver number of iteration done
(   70)       !!             - res    : solver residu reached
(   71)       !!             - gcx()  : solution of the elliptic system
(   72)       !!
(   73)       !! References :
(   74)       !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
(   75)       !!
(   76)       !! History :
(   77)       !!        !  90-10  (G. Madec)  Original code
(   78)       !!        !  91-11  (G. Madec)
(   79)       !!        !  93-04  (M. Guyon)  loops and suppress pointers
(   80)       !!        !  95-09  (M. Imbard, J. Escobar)  mpp exchange 
(   81)       !!        !  96-05  (G. Madec)  merge sor and pcg formulations
(   82)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(   83)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   84)       !!----------------------------------------------------------------------
(   85)       !! * Arguments
(   86)       INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
(   87)       !                                      ! gence is not reached: the model is
(   88)       !                                      ! stopped in step
(   89)       !                                      ! set to zero before the call of solpcg
(   90) 
(   91)       !! * Local declarations
(   92)       INTEGER ::   ji, jj, jn                ! dummy loop indices
(   93)       REAL(wp) ::   zgcad                    ! temporary scalars
(   94)       !!----------------------------------------------------------------------
(   95) 
(   96)       !                                                !================
(   97)       DO jn = 1, nmax                                  ! Iterative loop
(   98)          !                                             !================
(   99) 
(  100)          IF( jn == 1 ) THEN           ! Initialization of the algorithm
(  101) 
(  102)             CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! lateral boundary condition






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 7

(  103) 
(  104)             ! gcr   = gcb-a.gcx
(  105)             ! gcdes = gsr
(  106)             DO jj = 2, jpjm1
(  107)                DO ji = 2, jpim1   ! vector opt.
(  108)                   zgcad = bmask(ji,jj) * ( gcb(ji,jj  ) -                gcx(ji  ,jj  )   &
(  109)                      &                                  - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
(  110)                      &                                  - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
(  111)                      &                                  - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
(  112)                      &                                  - gcp(ji,jj,4) * gcx(ji  ,jj+1)   )
(  113)                   gcr  (ji,jj) = zgcad
(  114)                   gcdes(ji,jj) = zgcad
(  115)                END DO
(  116)             END DO
(  117)             
(  118)             rnorme = SUM(  gcr(:,:) * gcdmat(:,:) * gcr(:,:)  )
(  119)             IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  120)             rr = rnorme
(  121) 
(  122)          ENDIF
(  123)         
(  124)          !                             ! Algorithm
(  125)         
(  126)          CALL lbc_lnk( gcdes, c_solver_pt, 1. )   ! lateral boundary condition
(  127)         
(  128)          ! ... gccd = matrix . gcdes
(  129)          DO jj = 2, jpjm1
(  130)             DO ji = 2, jpim1   ! vector opt.
(  131)                gccd(ji,jj) = bmask(ji,jj)*( gcdes(ji,jj)   &
(  132)                   &        +gcp(ji,jj,1)*gcdes(ji,jj-1)+gcp(ji,jj,2)*gcdes(ji-1,jj)   &
(  133)                   &        +gcp(ji,jj,4)*gcdes(ji,jj+1)+gcp(ji,jj,3)*gcdes(ji+1,jj)   )
(  134)             END DO
(  135)          END DO
(  136)  
(  137)          ! alph = (gcr,gcr)/(gcdes,gccd)
(  138)          radd = SUM(  gcdes(:,:) * gcdmat(:,:) * gccd(:,:)  )
(  139)          IF( lk_mpp )   CALL mpp_sum( radd )   ! sum over the global domain
(  140)          alph = rr / radd
(  141)          
(  142)          ! gcx = gcx + alph * gcdes
(  143)          ! gcr = gcr - alph * gccd
(  144)          DO jj = 2, jpjm1
(  145)             DO ji = 2, jpim1   ! vector opt.
(  146)                gcx(ji,jj) = bmask(ji,jj) * ( gcx(ji,jj) + alph * gcdes(ji,jj) )
(  147)                gcr(ji,jj) = bmask(ji,jj) * ( gcr(ji,jj) - alph * gccd (ji,jj) )
(  148)             END DO
(  149)          END DO
(  150)         
(  151)          ! rnorme = (gcr,gcr)
(  152)          rnorme = SUM(  gcr(:,:) * gcdmat(:,:) * gcr(:,:)  )
(  153)          IF( lk_mpp )   CALL  mpp_sum( rnorme )   ! sum over the global domain
(  154)         
(  155)          ! test of convergence
(  156)          IF( rnorme < epsr .OR. jn == nmax ) THEN
(  157)             res = SQRT( rnorme )
(  158)             niter = jn
(  159)             ncut = 999
(  160)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 8

(  161)         
(  162)          ! beta = (rk+1,rk+1)/(rk,rk)
(  163)          beta = rnorme / rr
(  164)          rr   = rnorme
(  165) 
(  166)          ! indicator of non-convergence or explosion
(  167)          IF( jn == nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
(  168)          IF( ncut == 999 ) GOTO 999
(  169) 
(  170)          ! gcdes = gcr + beta * gcdes
(  171)          DO jj = 2, jpjm1
(  172)             DO ji = 2, jpim1   ! vector opt.
(  173)                gcdes(ji,jj) = bmask(ji,jj)*( gcr(ji,jj) + beta * gcdes(ji,jj) )
(  174)             END DO
(  175)          END DO
(  176)         
(  177)          !                                             !================
(  178)       END DO                                           !    End Loop
(  179)       !                                                !================
(  180)      
(  181) 999   CONTINUE
(  182)      
(  183)      
(  184)       ! Output in gcx with lateral b.c. applied
(  185)       ! ---------------------------------------
(  186)      
(  187)       CALL lbc_lnk( gcx, c_solver_pt, 1. )
(  188)      
(  189)    END SUBROUTINE sol_pcg
(  190) 
(  191)    !!=====================================================================
(  192) END MODULE solpcg
