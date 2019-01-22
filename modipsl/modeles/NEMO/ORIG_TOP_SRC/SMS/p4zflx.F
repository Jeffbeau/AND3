
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zflx.F,v 1.6 2005/11/14 16:42:40 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
CDIR$ LIST
      SUBROUTINE p4zflx
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC          ROUTINE p4zflx : PISCES MODEL
CCC          *****************************
CCC
CCC
CC     PURPOSE.
CC     --------
CC          *P4ZFLX* CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
CC
CC     EXTERNALS.
CC     ----------
CC          NONE.
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1988-07 E. MAIER-REIMER      MPI HAMBURG
CC      additions     : 1998    O. Aumont
CC      modifications : 1999    C. Le Quere
CC      modifications : 2004    O. Aumont
CC     -----------------------------------------------------------------
CC  parameters and commons
CC ======================
CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
#include "domzgr_substitute.h90"
CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER nspyr, ji, jj, krorr
      REAL zpdtan
      REAL kgco2(jpi,jpj),kgo2(jpi,jpj)
      REAL ttc, ws
      REAL fld, flu, oxy16, flu16, zfact
      REAL zph,ah2,zbot,zdic,zalk,schmitto2, zalka
      REAL schmittco2
C
C
C  1. ASSIGNATION TO EXPONENTS IN THE LISS AND MERLIVAT
C     FORMULATION OF THE GAS EXCHANGE RATE
c -----------------------------------------------------
C
      zpdtan = raass / rdt
      nspyr  = nint(zpdtan)
C
C* 1.1 SURFACE CHEMISTRY (PCO2 AND [H+] IN
C     SURFACE LAYER); THE RESULT OF THIS CALCULATION
C     IS USED TO COMPUTE AIR-SEA FLUX OF CO2
C ---------------------------------------------------
C
      DO krorr = 1,10
C
        DO jj = 1,jpj
          DO ji = 1,jpi
C
C* 1.2 DUMMY VARIABLES FOR DIC, H+, AND BORATE
C --------------------------------------------
C
        zbot = borat(ji,jj,1)
        zfact = rhop(ji,jj,1)/1000.+rtrn
        zdic  = trn(ji,jj,1,jpdic)/zfact
        zph = max(hi(ji,jj,1),1.E-10)/zfact
        zalka = trn(ji,jj,1,jptal)/zfact
C
C* 1.3 CALCULATE [ALK]([CO3--], [HCO3-])
C ------------------------------------
C
        zalk=zalka-
     &        (akw3(ji,jj,1)/zph-zph+zbot/(1.+zph/akb3(ji,jj,1)))
C
C* 1.4 CALCULATE [H+] AND [H2CO3]
C -----------------------------------------
C
         ah2=sqrt((zdic-zalk)**2+4*(zalk*ak23(ji,jj,1)
     &     /ak13(ji,jj,1))*(2*zdic-zalk))
        ah2=0.5*ak13(ji,jj,1)/zalk*((zdic-zalk)+ah2)
        h2co3(ji,jj) = (2*zdic-zalk)/(2.+ak13(ji,jj,1)/ah2)*zfact
        hi(ji,jj,1)  = ah2*zfact
          END DO
        END DO
      END DO
C
C
C 2. COMPUTE FLUXES
C --------------
C
C 2.1 FIRST COMPUTE GAS EXCHANGE COEFFICIENTS
C -------------------------------------------
C
      DO jj = 1,jpj
        DO ji = 1,jpi
C
          ttc = min(35.,tn(ji,jj,1))
          schmittco2=2073.1-125.62*ttc+3.6276*ttc**2
     &      -0.043126*ttc**3
          ws=vatm(ji,jj)
C
C 2.2 COMPUTE GAS EXCHANGE FOR CO2
C --------------------------------
C
          kgco2(ji,jj) = (0.3*ws*ws + 2.5*(0.5246+ttc*(0.016256+
     &      ttc*0.00049946)))*sqrt(660./schmittco2)
C
C 2.3 CONVERT TO m/s, and apply sea-ice cover
C -----------------------------------------------------
C
          kgco2(ji,jj) = kgco2(ji,jj)/(100.*3600.)
     &      *(1-freeze(ji,jj))*tmask(ji,jj,1)
C
         END DO
       END DO
C
C 2.5 COMPUTE GAS EXCHANGE COEFFICIENT FO O2 FROM
C      Waninkhof EQUATIONS
C -----------------------------------------------
C
       DO jj = 1,jpj
         DO ji = 1,jpi
C
          ws = vatm(ji,jj)
          schmitto2 = 1953.4-128.0*ttc+3.9918*ttc**2
     &      -0.050091*ttc**3

          kgo2(ji,jj) = (0.3*ws*ws + 2.5*(0.5246+ttc*(0.016256+
     &      ttc*0.00049946)))*sqrt(660./schmitto2)
C
C CONVERT TO m/s AND APPLY SEA ICE COVER
C -------------------------------------
C
          kgo2(ji,jj) = kgo2(ji,jj)/(100.*3600.)
     $      *(1-freeze(ji,jj))*tmask(ji,jj,1)
C
         ENDDO
       ENDDO
C
       DO jj = 1,jpj
         DO ji = 1,jpi
C
C Compute CO2 flux for the sea and air
C ------------------------------------
C
          fld = atcco2*tmask(ji,jj,1)*chemc(ji,jj,3)*kgco2(ji,jj)
          flu = h2co3(ji,jj)*tmask(ji,jj,1)*kgco2(ji,jj)
          tra(ji,jj,1,jpdic)= tra(ji,jj,1,jpdic)+(fld-flu)
     &      /fse3t(ji,jj,1)
C
C Compute O2 flux 
C ---------------
C
          oxy16 = trn(ji,jj,1,jpoxy)
          flu16 = (atcox*chemc(ji,jj,2)-oxy16)*kgo2(ji,jj)
          tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy)+flu16
     &      /fse3t(ji,jj,1)
C
C Save diagnostics
C ----------------
C
#    if defined key_trc_diaadd
          trc2d(ji,jj,1) = (fld-flu)*1000.
          trc2d(ji,jj,2) = flu16*1000.
          trc2d(ji,jj,3) = kgco2(ji,jj)
          trc2d(ji,jj,4) = atcco2-h2co3(ji,jj)/(chemc(ji,jj,1)+rtrn)
#    endif
C
        END DO
      END DO
C
#endif
      RETURN
      END
