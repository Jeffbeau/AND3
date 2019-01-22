
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zprg.F,v 1.6 2005/11/14 16:42:41 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
! CDIR$ LIST
       SUBROUTINE p4zprg(kt)
CCC---------------------------------------------------------------------
CCC
CCC           ROUTINE p4zprg : PISCES MODEL
CCC           *****************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Call Biological sources and sinks subroutines
CCC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC              all the common defined in opa
CC
CC
CC   OUTPUT :                   : no
CC   ------
CC
CC   WORKSPACE :
CC   ---------
CC
CC   EXTERNAL :
CC   --------
CC      p4zche, p4zint, p4zlys, p4zbio, p4zsed, p4zflx
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : O. AUMONT (2004)
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
! CDIR$ NOLIST
      USE trp_trc
      USE sms
CC
      USE oce_trc
      USE lbclnk
      USE lib_mpp
CC
      IMPLICIT NONE
! CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================

      INTEGER kt
#if defined key_passivetrc && defined key_trc_pisces
      INTEGER jnt, jn

C
C this part is without macrotasking coding
C
C Call an intermediate routine that in turns, calls chemistry
C and another routine on a daily basis
C -----------------------------------------------------------
C
      CALL p4zslow(kt)

C......................................................................
C
C Compute soft tissue production (POC)
C ------------------------------------

       do jnt=1,nrdttrc
C
         CALL p4zbio

C
C......................................................................
C
C Compute soft tissue remineralisation
C ------------------------------------
C
         CALL p4zsed

C
          trb=trn
        END DO
C
C......................................................................
C
C Compute CaCO3 saturation
C ------------------------
C
      CALL p4zlys

C
C......................................................................
C
C Compute surface fluxes
C ----------------------
C
      CALL p4zflx


      DO jn=1 , jptra
        CALL lbc_lnk(trn(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk(trb(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk(tra(:,:,:,jn), 'T', 1. )
      END DO

C
C......................................................................
C
#endif
C
      RETURN
      END
