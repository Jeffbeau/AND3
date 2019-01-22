


PGF90 (Version     10.5)          09/23/2010  18:36:34      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /distrib/local/netcdf/pgf/include/
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: par_kind.F90

(    1) MODULE par_kind
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE par_kind  ***
(    4)    !! Ocean :  define the kind of real for the whole model
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5   02/06  (G. Madec)  Original code
(    8)    !!----------------------------------------------------------------------
(    9)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/par_kind.F90,v 1.3 2005/03/27 18:34:48 opalod Exp $ 
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)    !!----------------------------------------------------------------------
(   13) 
(   14)    IMPLICIT NONE
(   15)    PRIVATE
(   16) 
(   17)    INTEGER, PUBLIC, PARAMETER ::    &  !:
(   18)       jpbyt   = 8       ,           &  !: real size for mpp communications
(   19)       jpbytda = 4       ,           &  !: real size in input data files 4 or 8
(   20)       jpbi3e  = 4                      !: real size for T3E
(   21) 
(   22)    ! Number model from which the SELECTED_*_KIND are requested:
(   23)    !             4 byte REAL       8 byte REAL
(   24)    ! CRAY:           -            precision = 13
(   25)    !                              exponent = 2465
(   26)    ! IEEE:      precision = 6     precision = 15
(   27)    !            exponent = 37     exponent = 307
(   28) 
(   29)    INTEGER, PUBLIC, PARAMETER ::        &  !: Floating point section
(   30)       sp = SELECTED_REAL_KIND( 6, 37),  &  !: single precision (real 4)
(   31)       dp = SELECTED_REAL_KIND(12,307),  &  !: double precision (real 8)
(   32)       wp = dp                              !: working precision
(   33) 
(   34)    INTEGER, PUBLIC, PARAMETER ::        &  !: Integer section
(   35)       i4 = SELECTED_INT_KIND(9) ,       &  !: single precision (integer 4)
(   36)       i8 = SELECTED_INT_KIND(14)           !: double precision (integer 8)
(   37) 
(   38) !!----------------------------------------------------------------------
(   39) !!----------------------------------------------------------------------
(   40)    PUBLIC SIGN, DDPDD
(   41) 
(   42)    INTERFACE SIGN
(   43)       MODULE PROCEDURE SIGN_SCALAR, SIGN_ARRAY_1D, SIGN_ARRAY_2D, SIGN_ARRAY_3D
(   44)    END INTERFACE
(   45) 
(   46)    INTERFACE DDPDD
(   47)       MODULE PROCEDURE DDPDD_S
(   48)    END INTERFACE
(   49) 






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 2

(   50) 
(   51) !!----------------------------------------------------------------------
(   52) ! Fortran 95 Standard has changed the behaviour of the SIGN(a,b) intrinsic function
(   53) ! Since the OPA model (particularly LIM) relies on the Fortran 90 behaviour
(   54) ! this function is introduced here to reproduce the Fortran 90 results.
(   55) ! Michael Dunphy 2007-Apr DunphyM@mar.dfo-mpo.gc.ca
(   56) ! Thanks to Zeliang Wang for the original scalar implementation of this function
(   57) !!----------------------------------------------------------------------
(   58) 
(   59)    CONTAINS
(   60) 
(   61)    FUNCTION SIGN_SCALAR(a,b)
(   62)       REAL(wp) :: a,b          ! input
(   63)       REAL(wp) :: SIGN_SCALAR  ! result
(   64) 
(   65)       IF ( b .GE. 0.e0) THEN
(   66)          SIGN_SCALAR = ABS(a)
(   67)       ELSE
(   68)          SIGN_SCALAR =-ABS(a)
(   69)       ENDIF
(   70) 
(   71)       RETURN
(   72) 
(   73)    END FUNCTION SIGN_SCALAR
(   74) 
(   75)    FUNCTION SIGN_ARRAY_1D(a,b) 
(   76)       IMPLICIT NONE
(   77)       REAL(wp) :: a(:),b(:)      ! input
(   78)       REAL(wp) :: SIGN_ARRAY_1D(SIZE(a,1))  ! result
(   79) 
(   80)       WHERE ( b .GE. 0.e0 )
(   81)          SIGN_ARRAY_1D = ABS(a)
(   82)       ELSEWHERE
(   83)          SIGN_ARRAY_1D =-ABS(a)
(   84)       END WHERE
(   85) 
(   86)       RETURN
(   87) 
(   88)    END FUNCTION SIGN_ARRAY_1D
(   89) 
(   90)    FUNCTION SIGN_ARRAY_2D(a,b) 
(   91)       IMPLICIT NONE
(   92)       REAL(wp) :: a(:,:),b(:,:)      ! input
(   93)       REAL(wp) :: SIGN_ARRAY_2D(SIZE(a,1),SIZE(a,2))  ! result
(   94) 
(   95)       WHERE ( b .GE. 0.e0 )
(   96)          SIGN_ARRAY_2D = ABS(a)
(   97)       ELSEWHERE
(   98)          SIGN_ARRAY_2D =-ABS(a)
(   99)       END WHERE
(  100) 
(  101)       RETURN
(  102) 
(  103)    END FUNCTION SIGN_ARRAY_2D
(  104) 
(  105)    FUNCTION SIGN_ARRAY_3D(a,b) 
(  106)       IMPLICIT NONE
(  107)       REAL(wp) :: a(:,:,:),b(:,:,:)      ! input






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 3

(  108)       REAL(wp) :: SIGN_ARRAY_3D(SIZE(a,1),SIZE(a,2),SIZE(a,3))  ! result
(  109) 
(  110)       WHERE ( b .GE. 0.e0 )
(  111)          SIGN_ARRAY_3D = ABS(a)
(  112)       ELSEWHERE
(  113)          SIGN_ARRAY_3D =-ABS(a)
(  114)       END WHERE
(  115) 
(  116)       RETURN
(  117) 
(  118)    END FUNCTION SIGN_ARRAY_3D
(  119) 
(  120)    !!---------------------------------------------------------------------
(  121)    !!   Routine DDPDD_S, performs double-double summation (scalar)
(  122)    !!
(  123)    !!   Modification of original codes written by David H. Bailey
(  124)    !!   This subroutine computes ddb(i) = dda(i)+ddb(i)
(  125)    !!---------------------------------------------------------------------
(  126)       SUBROUTINE DDPDD_S (dda, ddb, len, itype)
(  127)       IMPLICIT NONE
(  128) 
(  129)     !! * Arguments
(  130)       INTEGER, INTENT(in)     :: len, itype
(  131)       COMPLEX, INTENT(in)     :: dda
(  132)       COMPLEX, INTENT(inout)  :: ddb
(  133) 
(  134)     !! * Local variables
(  135)       REAL :: e, t1, t2  ! local work variables
(  136) 
(  137)    ! Compute dda + ddb using Knuth's trick.
(  138)         t1 = real(dda) + real(ddb)
(  139)         e = t1 - real(dda)
(  140)         t2 = ((real(ddb) - e) + (real(dda) - (t1 - e))) &
(  141)              +imag(dda) + imag(ddb)
(  142) 
(  143)    ! The result is t1 + t2, after normalization.
(  144)         ddb = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
(  145)       END SUBROUTINE DDPDD_S
(  146)    !!---------------------------------------------------------------------
(  147) 
(  148) 
(  149) END MODULE par_kind






















PGF90 (Version     10.5)          09/23/2010  18:36:34      page 4

(    1) # 1 "par_kind.F90"
(    1) MODULE par_kind
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE par_kind  ***
(    4)    !! Ocean :  define the kind of real for the whole model
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5   02/06  (G. Madec)  Original code
(    8)    !!----------------------------------------------------------------------
(    9)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/par_kind.F90,v 1.3 2005/03/27 18:34:48 opalod Exp $ 
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)    !!----------------------------------------------------------------------
(   13) 
(   14)    IMPLICIT NONE
(   15)    PRIVATE
(   16) 
(   17)    INTEGER, PUBLIC, PARAMETER ::    &  !:
(   18)       jpbyt   = 8       ,           &  !: real size for mpp communications
(   19)       jpbytda = 4       ,           &  !: real size in input data files 4 or 8
(   20)       jpbi3e  = 4                      !: real size for T3E
(   21) 
(   22)    ! Number model from which the SELECTED_*_KIND are requested:
(   23)    !             4 byte REAL       8 byte REAL
(   24)    ! CRAY:           -            precision = 13
(   25)    !                              exponent = 2465
(   26)    ! IEEE:      precision = 6     precision = 15
(   27)    !            exponent = 37     exponent = 307
(   28) 
(   29)    INTEGER, PUBLIC, PARAMETER ::        &  !: Floating point section
(   30)       sp = SELECTED_REAL_KIND( 6, 37),  &  !: single precision (real 4)
(   31)       dp = SELECTED_REAL_KIND(12,307),  &  !: double precision (real 8)
(   32)       wp = dp                              !: working precision
(   33) 
(   34)    INTEGER, PUBLIC, PARAMETER ::        &  !: Integer section
(   35)       i4 = SELECTED_INT_KIND(9) ,       &  !: single precision (integer 4)
(   36)       i8 = SELECTED_INT_KIND(14)           !: double precision (integer 8)
(   37) 
(   38) !!----------------------------------------------------------------------
(   39) !!----------------------------------------------------------------------
(   40)    PUBLIC SIGN, DDPDD
(   41) 
(   42)    INTERFACE SIGN
(   43)       MODULE PROCEDURE SIGN_SCALAR, SIGN_ARRAY_1D, SIGN_ARRAY_2D, SIGN_ARRAY_3D
(   44)    END INTERFACE
(   45) 
(   46)    INTERFACE DDPDD
(   47)       MODULE PROCEDURE DDPDD_S
(   48)    END INTERFACE
(   49) 
(   50) 
(   51) !!----------------------------------------------------------------------
(   52) ! Fortran 95 Standard has changed the behaviour of the SIGN(a,b) intrinsic function
(   53) ! Since the OPA model (particularly LIM) relies on the Fortran 90 behaviour
(   54) ! this function is introduced here to reproduce the Fortran 90 results.
(   55) ! Michael Dunphy 2007-Apr DunphyM@mar.dfo-mpo.gc.ca
(   56) ! Thanks to Zeliang Wang for the original scalar implementation of this function
(   57) !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 5

(   58) 
(   59)    CONTAINS
(   60) 
(   61)    FUNCTION SIGN_SCALAR(a,b)
(   62)       REAL(wp) :: a,b          ! input
(   63)       REAL(wp) :: SIGN_SCALAR  ! result
(   64) 
(   65)       IF ( b .GE. 0.e0) THEN
(   66)          SIGN_SCALAR = ABS(a)
(   67)       ELSE
(   68)          SIGN_SCALAR =-ABS(a)
(   69)       ENDIF
(   70) 
(   71)       RETURN
(   72) 
(   73)    END FUNCTION SIGN_SCALAR
(   74) 
(   75)    FUNCTION SIGN_ARRAY_1D(a,b) 
(   76)       IMPLICIT NONE
(   77)       REAL(wp) :: a(:),b(:)      ! input
(   78)       REAL(wp) :: SIGN_ARRAY_1D(SIZE(a,1))  ! result
(   79) 
(   80)       WHERE ( b .GE. 0.e0 )
(   81)          SIGN_ARRAY_1D = ABS(a)
(   82)       ELSEWHERE
(   83)          SIGN_ARRAY_1D =-ABS(a)
(   84)       END WHERE
(   85) 
(   86)       RETURN
(   87) 
(   88)    END FUNCTION SIGN_ARRAY_1D
(   89) 
(   90)    FUNCTION SIGN_ARRAY_2D(a,b) 
(   91)       IMPLICIT NONE
(   92)       REAL(wp) :: a(:,:),b(:,:)      ! input
(   93)       REAL(wp) :: SIGN_ARRAY_2D(SIZE(a,1),SIZE(a,2))  ! result
(   94) 
(   95)       WHERE ( b .GE. 0.e0 )
(   96)          SIGN_ARRAY_2D = ABS(a)
(   97)       ELSEWHERE
(   98)          SIGN_ARRAY_2D =-ABS(a)
(   99)       END WHERE
(  100) 
(  101)       RETURN
(  102) 
(  103)    END FUNCTION SIGN_ARRAY_2D
(  104) 
(  105)    FUNCTION SIGN_ARRAY_3D(a,b) 
(  106)       IMPLICIT NONE
(  107)       REAL(wp) :: a(:,:,:),b(:,:,:)      ! input
(  108)       REAL(wp) :: SIGN_ARRAY_3D(SIZE(a,1),SIZE(a,2),SIZE(a,3))  ! result
(  109) 
(  110)       WHERE ( b .GE. 0.e0 )
(  111)          SIGN_ARRAY_3D = ABS(a)
(  112)       ELSEWHERE
(  113)          SIGN_ARRAY_3D =-ABS(a)
(  114)       END WHERE
(  115) 






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 6

(  116)       RETURN
(  117) 
(  118)    END FUNCTION SIGN_ARRAY_3D
(  119) 
(  120)    !!---------------------------------------------------------------------
(  121)    !!   Routine DDPDD_S, performs double-double summation (scalar)
(  122)    !!
(  123)    !!   Modification of original codes written by David H. Bailey
(  124)    !!   This subroutine computes ddb(i) = dda(i)+ddb(i)
(  125)    !!---------------------------------------------------------------------
(  126)       SUBROUTINE DDPDD_S (dda, ddb, len, itype)
(  127)       IMPLICIT NONE
(  128) 
(  129)     !! * Arguments
(  130)       INTEGER, INTENT(in)     :: len, itype
(  131)       COMPLEX, INTENT(in)     :: dda
(  132)       COMPLEX, INTENT(inout)  :: ddb
(  133) 
(  134)     !! * Local variables
(  135)       REAL :: e, t1, t2  ! local work variables
(  136) 
(  137)    ! Compute dda + ddb using Knuth's trick.
(  138)         t1 = real(dda) + real(ddb)
(  139)         e = t1 - real(dda)
(  140)         t2 = ((real(ddb) - e) + (real(dda) - (t1 - e))) &
(  141)              +imag(dda) + imag(ddb)
(  142) 
(  143)    ! The result is t1 + t2, after normalization.
(  144)         ddb = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
(  145)       END SUBROUTINE DDPDD_S
(  146)    !!---------------------------------------------------------------------
(  147) 
(  148) 
(  149) END MODULE par_kind
