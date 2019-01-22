MODULE par_kind
   !!======================================================================
   !!                   ***  MODULE par_kind  ***
   !! Ocean :  define the kind of real for the whole model
   !!======================================================================
   !! History :
   !!   8.5   02/06  (G. Madec)  Original code
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/par_kind.F90,v 1.3 2005/03/27 18:34:48 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE

   INTEGER, PUBLIC, PARAMETER ::    &  !:
      jpbyt   = 8       ,           &  !: real size for mpp communications
      jpbytda = 4       ,           &  !: real size in input data files 4 or 8
      jpbi3e  = 4                      !: real size for T3E

   ! Number model from which the SELECTED_*_KIND are requested:
   !             4 byte REAL       8 byte REAL
   ! CRAY:           -            precision = 13
   !                              exponent = 2465
   ! IEEE:      precision = 6     precision = 15
   !            exponent = 37     exponent = 307

   INTEGER, PUBLIC, PARAMETER ::        &  !: Floating point section
      sp = SELECTED_REAL_KIND( 6, 37),  &  !: single precision (real 4)
      dp = SELECTED_REAL_KIND(12,307),  &  !: double precision (real 8)
      wp = dp                              !: working precision

   INTEGER, PUBLIC, PARAMETER ::        &  !: Integer section
      i4 = SELECTED_INT_KIND(9) ,       &  !: single precision (integer 4)
      i8 = SELECTED_INT_KIND(14)           !: double precision (integer 8)

!!----------------------------------------------------------------------
!!----------------------------------------------------------------------
   PUBLIC SIGN, DDPDD

   INTERFACE SIGN
      MODULE PROCEDURE SIGN_SCALAR, SIGN_ARRAY_1D, SIGN_ARRAY_2D, SIGN_ARRAY_3D
   END INTERFACE

   INTERFACE DDPDD
      MODULE PROCEDURE DDPDD_S
   END INTERFACE


!!----------------------------------------------------------------------
! Fortran 95 Standard has changed the behaviour of the SIGN(a,b) intrinsic function
! Since the OPA model (particularly LIM) relies on the Fortran 90 behaviour
! this function is introduced here to reproduce the Fortran 90 results.
! Michael Dunphy 2007-Apr DunphyM@mar.dfo-mpo.gc.ca
! Thanks to Zeliang Wang for the original scalar implementation of this function
!!----------------------------------------------------------------------

   CONTAINS

   FUNCTION SIGN_SCALAR(a,b)
      REAL(wp) :: a,b          ! input
      REAL(wp) :: SIGN_SCALAR  ! result

      IF ( b .GE. 0.e0) THEN
         SIGN_SCALAR = ABS(a)
      ELSE
         SIGN_SCALAR =-ABS(a)
      ENDIF

      RETURN

   END FUNCTION SIGN_SCALAR

   FUNCTION SIGN_ARRAY_1D(a,b) 
      IMPLICIT NONE
      REAL(wp) :: a(:),b(:)      ! input
      REAL(wp) :: SIGN_ARRAY_1D(SIZE(a,1))  ! result

      WHERE ( b .GE. 0.e0 )
         SIGN_ARRAY_1D = ABS(a)
      ELSEWHERE
         SIGN_ARRAY_1D =-ABS(a)
      END WHERE

      RETURN

   END FUNCTION SIGN_ARRAY_1D

   FUNCTION SIGN_ARRAY_2D(a,b) 
      IMPLICIT NONE
      REAL(wp) :: a(:,:),b(:,:)      ! input
      REAL(wp) :: SIGN_ARRAY_2D(SIZE(a,1),SIZE(a,2))  ! result

      WHERE ( b .GE. 0.e0 )
         SIGN_ARRAY_2D = ABS(a)
      ELSEWHERE
         SIGN_ARRAY_2D =-ABS(a)
      END WHERE

      RETURN

   END FUNCTION SIGN_ARRAY_2D

   FUNCTION SIGN_ARRAY_3D(a,b) 
      IMPLICIT NONE
      REAL(wp) :: a(:,:,:),b(:,:,:)      ! input
      REAL(wp) :: SIGN_ARRAY_3D(SIZE(a,1),SIZE(a,2),SIZE(a,3))  ! result

      WHERE ( b .GE. 0.e0 )
         SIGN_ARRAY_3D = ABS(a)
      ELSEWHERE
         SIGN_ARRAY_3D =-ABS(a)
      END WHERE

      RETURN

   END FUNCTION SIGN_ARRAY_3D

   !!---------------------------------------------------------------------
   !!   Routine DDPDD_S, performs double-double summation (scalar)
   !!
   !!   Modification of original codes written by David H. Bailey
   !!   This subroutine computes ddb(i) = dda(i)+ddb(i)
   !!---------------------------------------------------------------------
      SUBROUTINE DDPDD_S (dda, ddb, len, itype)
      IMPLICIT NONE

    !! * Arguments
      INTEGER, INTENT(in)     :: len, itype
      COMPLEX, INTENT(in)     :: dda
      COMPLEX, INTENT(inout)  :: ddb

    !! * Local variables
      REAL :: e, t1, t2  ! local work variables

   ! Compute dda + ddb using Knuth's trick.
        t1 = real(dda) + real(ddb)
        e = t1 - real(dda)
        t2 = ((real(ddb) - e) + (real(dda) - (t1 - e))) &
             +imag(dda) + imag(ddb)

   ! The result is t1 + t2, after normalization.
        ddb = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
      END SUBROUTINE DDPDD_S
   !!---------------------------------------------------------------------


END MODULE par_kind
