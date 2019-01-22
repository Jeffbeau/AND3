


PGF90 (Version     10.5)          09/23/2010  18:36:34      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /distrib/local/netcdf/pgf/include/
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: lib_cray.f90

(    1) !  Cray subroutines or functions used by OPA model and possibly 
(    2) !  not found on other platforms.
(    3) !
(    4) !  check their existence
(    5) !  
(    6) !  sdot
(    7) !  wheneq
(    8) !  saxpy
(    9) !  isrchne
(   10)    !!----------------------------------------------------------------------
(   11)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   12)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_cray.f90,v 1.3 2005/03/27 18:34:47 opalod Exp $ 
(   13)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   14)    !!----------------------------------------------------------------------
(   15) !---------------------------------------------------------
(   16)      FUNCTION sdot( I, X, J, Y, K )
(   17)         DIMENSION X(1), Y(1)
(   18)         SDOT = 0.
(   19)         DO N = 1, I
(   20)         SDOT = SDOT + X(1+(N-1)*J) * Y(1+(N-1)*K)
(   21)         END DO
(   22)      END FUNCTION sdot
(   23) !---------------------------------------------------------
(   24)      SUBROUTINE wheneq ( i, x, j, t, ind, nn )
(   25)         IMPLICIT NONE
(   26) 
(   27)         INTEGER , INTENT (  in ) :: i, j
(   28)         INTEGER , INTENT ( out ) :: nn
(   29)         REAL    , INTENT (  in ), DIMENSION (1+(i-1)*j) :: x
(   30)         REAL    , INTENT (  in ) :: t
(   31)         INTEGER , INTENT ( out ), DIMENSION (1+(i-1)*j) :: ind
(   32)         INTEGER :: n, k
(   33)         nn = 0
(   34)         DO n = 1, i
(   35)           k = 1 + (n-1) * j
(   36)           IF ( x ( k) == t ) THEN 
(   37)               nn = nn + 1
(   38)               ind (nn) = k
(   39)           ENDIF
(   40)         END DO 
(   41) 
(   42)      END SUBROUTINE wheneq
(   43) !---------------------------------------------------------
(   44)      SUBROUTINE saxpy( I, A, X, J, Y, K )
(   45)         DIMENSION X(1),Y(1)
(   46)         DO N = 1, I
(   47)            Y(1+(N-1)*K)=A*X(1+(N-1)*J)+Y(1+(N-1)*K)
(   48)         END DO
(   49)      END SUBROUTINE saxpy






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 2

(   50) !---------------------------------------------------------
(   51)      FUNCTION isrchne( K, X, I, B )
(   52)         DIMENSION X(1)
(   53)         DO N = 1, K
(   54)            IF( X(1+(N-1)*I) /= B ) THEN
(   55)               ISRCHNE = N
(   56)               RETURN
(   57)            ELSE
(   58)               ISRCHNE = N + 1
(   59)            ENDIF
(   60)         END DO
(   61)      END FUNCTION isrchne
