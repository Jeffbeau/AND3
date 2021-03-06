
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zlim.F,v 1.2 2005/11/14 16:42:40 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
CDIR$ LIST
      SUBROUTINE p4zlim
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC             ROUTINE p4zlim : PISCES MODEL
CCC             *****************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC         Compute the co-limitations by the various nutrients
CCC         for the various phytoplankton species
CCC
CC   INPUT :
CC   -----
CC      argument
CC              None
CC      common
CC              all the common defined in opa
CC
CC
CC   OUTPUT :                   : no
CC   ------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : O. Aumont (2004)
CC----------------------------------------------------------------------
CC parameters and commons
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
      INTEGER ji, jj, jk
      REAL xlim1,xlim2,xlim3,xlim4
      REAL xconctemp,xconctemp2,xconctempn,xconctempn2
C
C  Tuning of the iron concentration to a minimum
C  level that is set to the detection limit
C  -------------------------------------
C

        trn(:,:,:,jpfer)=max(trn(:,:,:,jpfer),1.E-11)
C
C  Computation of a variable Ks for iron on diatoms
C  taking into account that increasing biomass is
C  made of generally bigger cells
C  ------------------------------------------------
C
        DO jk=1,jpkm1
          DO jj=1,jpj
            DO ji=1,jpi
        xconctemp=max(0.,trn(ji,jj,jk,jpdia)-5E-7)
        xconctemp2=min(5.E-7,trn(ji,jj,jk,jpdia))
        xconctempn=max(0.,trn(ji,jj,jk,jpphy)-1E-6)
        xconctempn2=min(1.E-6,trn(ji,jj,jk,jpphy))
        concdfe(ji,jj,jk)=(xconctemp2*conc3+0.4E-9*
     .    xconctemp)/(xconctemp2+xconctemp+rtrn)
        concdfe(ji,jj,jk)=max(conc3,concdfe(ji,jj,jk))
        concnfe(ji,jj,jk)=(xconctempn2*conc2+0.08E-9*
     .    xconctempn)/(xconctempn2+xconctempn+rtrn)
        concnfe(ji,jj,jk)=max(conc2,concnfe(ji,jj,jk))
            END DO
          END DO
        END DO
C
        DO jk = 1,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
C    
C      Michaelis-Menten Limitation term for nutrients
C      Small flagellates
C      -----------------------------------------------
C
        xnanono3(ji,jj,jk)=trn(ji,jj,jk,jpno3)*concnnh4
     &      /(conc0*concnnh4+concnnh4*trn(ji,jj,jk,jpno3)+
     &        conc0*trn(ji,jj,jk,jpnh4))
        xnanonh4(ji,jj,jk)=trn(ji,jj,jk,jpnh4)*conc0
     &      /(conc0*concnnh4+concnnh4*trn(ji,jj,jk,jpno3)+
     &        conc0*trn(ji,jj,jk,jpnh4))
        xlim1=xnanono3(ji,jj,jk)+xnanonh4(ji,jj,jk)
        xlim2=trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+concnnh4)
        xlim3=trn(ji,jj,jk,jpfer)/(trn(ji,jj,jk,jpfer)
     &     +concnfe(ji,jj,jk))
        xlimphy(ji,jj,jk)=min(xlim1,xlim2,xlim3)
        xlim1=trn(ji,jj,jk,jpnh4)/(concnnh4+trn(ji,jj,jk,jpnh4))
        xlim3=trn(ji,jj,jk,jpfer)/(trn(ji,jj,jk,jpfer)+conc2)
        xlim4=trn(ji,jj,jk,jpdoc)/(trn(ji,jj,jk,jpdoc)+xkdoc2)
        xlimbac(ji,jj,jk)=min(xlim1,xlim2,xlim3)*xlim4
C
            END DO
          END DO
        END DO
C
        DO jk = 1,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
C
C   Michaelis-Menten Limitation term for nutrients
C   Diatoms
C   ----------------------------------------------
C
        xdiatno3(ji,jj,jk)=trn(ji,jj,jk,jpno3)*concdnh4
     &      /(conc1*concdnh4+concdnh4*trn(ji,jj,jk,jpno3)+
     &        conc1*trn(ji,jj,jk,jpnh4))
        xdiatnh4(ji,jj,jk)=trn(ji,jj,jk,jpnh4)*conc1
     &      /(conc1*concdnh4+concdnh4*trn(ji,jj,jk,jpno3)+
     &        conc1*trn(ji,jj,jk,jpnh4))

        xlim1=xdiatno3(ji,jj,jk)+xdiatnh4(ji,jj,jk)
        xlim2=trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+concdnh4)
        xlim3=trn(ji,jj,jk,jpsil)/(trn(ji,jj,jk,jpsil)+xksi(ji,jj))
        xlim4=trn(ji,jj,jk,jpfer)/(trn(ji,jj,jk,jpfer)
     &      +concdfe(ji,jj,jk))
        xlimdia(ji,jj,jk)=min(xlim1,xlim2,xlim3,xlim4)
        xlim4=trn(ji,jj,jk,jpfer)/(trn(ji,jj,jk,jpfer)+conc3)
        xlimdia2(ji,jj,jk)=min(xlim1,xlim2,xlim3,xlim4)
C
            END DO
          END DO
        END DO
C
#endif
      RETURN
      END
