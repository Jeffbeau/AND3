
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zslow.F,v 1.2 2005/11/14 16:42:42 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
! CDIR$ LIST
      SUBROUTINE p4zslow(kt)
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC          ROUTINE p4zslow : PISCES MODEL
CCC          *****************************
CCC
CCC     PURPOSE.
CCC     --------
CCC          *P4ZSLOW* CALL VARIOUS ROUTINES ON A DAILY BASIS
CCC
CCC
CC     EXTERNALS.
CC     ----------
CC          p4zche, p4zint
CC
CC   MODIFICATIONS:
CC   --------------
CC      NEW : 2004 O. Aumont
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
! CDIR$ nolist
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
! CDIR$ list
CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER iyy, imm, idd, kt
C
      iyy = ndastp/10000
      imm = (ndastp - iyy*10000)/100
      idd = (ndastp - iyy*10000 - imm*100)

       IF (ndayflxtr /= idd) THEN
           ndayflxtr = idd
C
C* 1. COMPUTATION OF THE CHEMICAL CONSTANTS
C  ----------------------------------------
C
       CALL p4zche
C
C  2. COMPUTATION OF VARIOUS RATES FOR BIOGEOCHEMISTRY
C  ---------------------------------------------------
C
       CALL p4zint(kt)
C
      ENDIF
C      
#endif
C
      RETURN
      END
