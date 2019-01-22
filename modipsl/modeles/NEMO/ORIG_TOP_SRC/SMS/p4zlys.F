
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zlys.F,v 1.6 2005/11/14 16:42:40 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
CDIR$ LIST
       SUBROUTINE p4zlys
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC        ROUTINE p4zlys : PISCES MODEL
CCC        *****************************
CCC
CCC
CCC     PURPOSE.
CCC     --------
CCC          *P4ZLYS*  CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
CCC                    COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
CCC                    OF CACO3 TO THE CACO3 SEDIMENT POOL.
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
CC ---------------------------------------------------------------------------
CC parameters and commons
CC ======================
CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk, jn
      REAL zbot, zalk, zdic, zph, remco3, ah2
      REAL delco3, excess, dispot, zfact, zalka
C
C
C* 1.1  BEGIN OF ITERATION
C ------------------------
C
      DO jn = 1,5
C
C* 1.2  COMPUTE [CO3--] and [H+] CONCENTRATIONS
C -------------------------------------------
C
      DO jk = 1,jpkm1
        DO jj=1,jpj
          DO ji = 1, jpi
C
C* 1.3  SET DUMMY VARIABLE FOR TOTAL BORATE
C -----------------------------------------
C
        zbot = borat(ji,jj,jk)
        zfact=rhop(ji,jj,jk)/1000.+rtrn
C
C* 1.4  SET DUMMY VARIABLE FOR [H+]
C ---------------------------------
C
        zph = hi(ji,jj,jk)*tmask(ji,jj,jk)/zfact
     &    +(1.-tmask(ji,jj,jk))*1.e-9
C
C* 1.5  SET DUMMY VARIABLE FOR [SUM(CO2)]GIVEN 
C -------------------------------------------
C
        zdic=trn(ji,jj,jk,jpdic)/zfact
        zalka=trn(ji,jj,jk,jptal)/zfact
C
C* 1.6 CALCULATE [ALK]([CO3--], [HCO3-])
C ------------------------------------
C
        zalk=zalka-(akw3(ji,jj,jk)/zph-zph
     &     +zbot/(1.+zph/akb3(ji,jj,jk)))
C
C* 2.10 CALCULATE [H+] and [CO3--]
C -----------------------------------------
C
        ah2=sqrt((zdic-zalk)*(zdic-zalk)+
     &     4.*(zalk*ak23(ji,jj,jk)/ak13(ji,jj,jk))
     &     *(2*zdic-zalk))
C
        ah2=0.5*ak13(ji,jj,jk)/zalk*((zdic-zalk)+ah2)
        co3(ji,jj,jk) = zalk/(2.+ah2/ak23(ji,jj,jk))*zfact

        hi(ji,jj,jk)  = ah2*zfact
C
          ENDDO
        ENDDO
      END DO
C
      END DO 
C
C     ---------------------------------------------------------
C*    2. CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
C        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
C        MGCO3)
C     ---------------------------------------------------------
C
      DO jk = 1,jpkm1
        DO jj = 1,jpj
          DO ji = 1, jpi
C
C* 2.1  DEVIATION OF [CO3--] FROM SATURATION VALUE
C ------------------------------------------------
C
            delco3 = co3(ji,jj,jk)-aksp(ji,jj,jk)/calcon
C
C* 2.2  SET DEGREE OF UNDER-/SUPERSATURATION
C ------------------------------------------
C
            excess = max(0.,-delco3)
C
C* 2.3  AMOUNT CACO3 (12C) THAT RE-ENTERS SOLUTION
C       (ACCORDING TO THIS FORMULATION ALSO SOME PARTICULATE
C       CACO3 GETS DISSOLVED EVEN IN THE CASE OF OVERSATURATION)
C --------------------------------------------------------------
C
            dispot = trn(ji,jj,jk,jpcal)*
     &        excess/(dispo0+excess)
#    if defined key_off_degrad
     &        *facvol(ji,jj,jk)
#    endif
C
C* 2.4  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
C       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
C -----------------------------------------------------------
C
            remco3=dispot/rmoss
            co3(ji,jj,jk) = co3(ji,jj,jk)+
     &        remco3*rfact
            tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal)+
     &        2.*remco3
            tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal)-
     &        remco3
            tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic)+
     &        remco3
C
          ENDDO
        ENDDO
      END DO

#    if defined key_trc_dia3d
         trc3d(:,:,:,1) = rhop(:,:,:)
         trc3d(:,:,:,2) = co3(:,:,:)
         trc3d(:,:,:,3) = aksp(:,:,:)/calcon
#    endif

C
#endif
      RETURN
      END
