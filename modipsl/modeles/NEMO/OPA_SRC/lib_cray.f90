!  Cray subroutines or functions used by OPA model and possibly 
!  not found on other platforms.
!
!  check their existence
!  
!  sdot
!  wheneq
!  saxpy
!  isrchne
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_cray.f90,v 1.3 2005/03/27 18:34:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
!---------------------------------------------------------
     FUNCTION sdot( I, X, J, Y, K )
        DIMENSION X(1), Y(1)
        SDOT = 0.
        DO N = 1, I
        SDOT = SDOT + X(1+(N-1)*J) * Y(1+(N-1)*K)
        END DO
     END FUNCTION sdot
!---------------------------------------------------------
     SUBROUTINE wheneq ( i, x, j, t, ind, nn )
        IMPLICIT NONE

        INTEGER , INTENT (  in ) :: i, j
        INTEGER , INTENT ( out ) :: nn
        REAL    , INTENT (  in ), DIMENSION (1+(i-1)*j) :: x
        REAL    , INTENT (  in ) :: t
        INTEGER , INTENT ( out ), DIMENSION (1+(i-1)*j) :: ind
        INTEGER :: n, k
        nn = 0
        DO n = 1, i
          k = 1 + (n-1) * j
          IF ( x ( k) == t ) THEN 
              nn = nn + 1
              ind (nn) = k
          ENDIF
        END DO 

     END SUBROUTINE wheneq
!---------------------------------------------------------
     SUBROUTINE saxpy( I, A, X, J, Y, K )
        DIMENSION X(1),Y(1)
        DO N = 1, I
           Y(1+(N-1)*K)=A*X(1+(N-1)*J)+Y(1+(N-1)*K)
        END DO
     END SUBROUTINE saxpy
!---------------------------------------------------------
     FUNCTION isrchne( K, X, I, B )
        DIMENSION X(1)
        DO N = 1, K
           IF( X(1+(N-1)*I) /= B ) THEN
              ISRCHNE = N
              RETURN
           ELSE
              ISRCHNE = N + 1
           ENDIF
        END DO
     END FUNCTION isrchne
