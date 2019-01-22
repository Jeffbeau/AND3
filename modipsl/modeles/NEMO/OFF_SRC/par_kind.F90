MODULE par_kind
   !!======================================================================
   !!                   ***  MODULE par_kind  ***
   !! Ocean :  define the kind of real for the whole model
   !!======================================================================
   !! History :
   !!   8.5   02/06  (G. Madec)  Original code
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/par_kind.F90,v 1.1.1.1 2005/11/14 10:41:07 opalod Exp $ 
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
END MODULE par_kind
