
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zint.F,v 1.6 2005/11/14 16:42:40 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
! CDIR$ LIST
      SUBROUTINE p4zint(kt)
#if defined key_passivetrc && defined key_trc_pisces
CCC
CCC
CCC       ROUTINE p4zint : PISCES MODEL
CCC       *****************************
CCC
CC
CC   PURPOSE :
CC   ---------
CC           *P4ZINT* INTERPOLATION AND COMPUTATION OF 
CC                    VARIOUS ACCESSORY FIELDS
CC   INPUT :
CC   -----
CC      argument
CC              kt              : time step
CC
CC   EXTERNAL :
CC   ----------
CC          NONE
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 2004    O. Aumont
CC ----------------------------------------------------------------
CC parameters and commons
CC ======================
! CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
! CDIR$ LIST
CC-----------------------------------------------------------------
CC------
CC local declarations
CC ==================
C
      INTEGER kt
      INTEGER ji, jj
      INTEGER iman
      INTEGER nspyr,nvit1t,nvit2t
      REAL zpdtan, zman, zpdtmo, zdemi
      REAL zt, zdum
C
C
      zpdtan = raass / rdt
      nspyr  = nint(zpdtan)
      zman   = 12.
      iman   = 12
      zpdtmo = zpdtan / zman
      zdemi  = zpdtmo / 2.
      zt     = ( float ( kt) + zdemi) / zpdtmo
      
C  recherche de l'indice des enregistrements
C  du modele dynamique encadrant le pas de temps kt.
C  --------------------------------------------------
C
      xtvit = zt - float(int ( zt))
      nvit1t = int (zt)
      nvit2t = nvit1t+1
      nvit1t = MOD ( nvit1t, iman)
      IF ( nvit1t .EQ. 0 ) nvit1t = iman
      nvit2t = MOD ( nvit2t, iman)
      IF ( nvit2t .EQ. 0 ) nvit2t = iman
C
C Interpolation of dust deposition
C --------------------------------
C
         dust(:,:) = (1.-xtvit)*dustmo(:,:,nvit1t)
     $            +xtvit*dustmo(:,:,nvit2t)
C
C
C Computation of phyto and zoo metabolic rate
C -------------------------------------------
C

         Tgfunc(:,:,:) = exp(0.063913*tn(:,:,:))
         Tgfunc2(:,:,:) = exp(0.07608*tn(:,:,:))
C
C      Computation of the silicon dependant half saturation
C      constant for silica uptake
C       ---------------------------------------------------
C
        DO ji=1,jpi
          DO jj=1,jpj
          zdum=trn(ji,jj,1,jpsil)**2
          xksimax(ji,jj) = max(xksimax(ji,jj),(1.+7.*zdum
     &      /(xksi2*xksi2*25.+zdum))*1E-6)
          END DO
        END DO
C
        IF (nday_year.EQ.365) THEN
           xksi=xksimax
           xksimax=0.
        ENDIF
C
#endif
C
      RETURN
      END
