
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zrem.F,v 1.6 2005/11/14 16:42:41 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
! CDIR$ LIST
      SUBROUTINE p4zrem
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC          ROUTINE p4zrem : PISCES MODEL
CCC          *****************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC         Compute remineralization/scavenging of organic compounds
CCC
CC   INPUT :
CC   -----
CC      common
CC              all the common defined in opa
CC
CC
CC   OUTPUT :                   : no
CC   ------
CC
CC   EXTERNAL :
CC   --------
CC            None
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : 2004 - O. Aumont 
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
! CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
#include "domzgr_substitute.h90"
! CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj, jk
      REAL remip,remik,xlam1b
      REAL xkeq,xfeequi,siremin
      REAL zsatur,zsatur2,znusil,zdepbac(jpi,jpj,jpk)
      REAL zlamfac,zstep,fesatur(jpi,jpj,jpk)
C
C      Time step duration for the biology
C
       zstep=rfact2/rjjss
C
C      Computation of the mean phytoplankton concentration as
C      a crude estimate of the bacterial biomass
C      --------------------------------------------------
C
        DO jk=1,12
         zdepbac(:,:,jk)=min(0.7*(trn(:,:,jk,jpzoo)+2*trn(:,:,jk,jpmes))
     &     ,4E-6)
        END DO
C
C      Vertical decay of the bacterial activity
C      ----------------------------------------
C
         do jk=13,jpk
           do jj=1,jpj
             do ji=1,jpi
         zdepbac(ji,jj,jk)=min(1.,fsdept(ji,jj,12)/fsdept(ji,jj,jk))
     &      *zdepbac(ji,jj,12)
             end do
           end do
         end do

         DO jk = 1,jpkm1
           DO jj = 1,jpj
             DO ji = 1,jpi
C
C    DENITRIFICATION FACTOR COMPUTED FROM O2 LEVELS
C    ----------------------------------------------
C
          nitrfac(ji,jj,jk)=
     &      max(0.,0.4*(6.E-6-trn(ji,jj,jk,jpoxy))/(oxymin+
     &      trn(ji,jj,jk,jpoxy)))
             END DO
           END DO
         END DO

          nitrfac(:,:,:)=min(1.,nitrfac(:,:,:))

         DO jk = 1,jpkm1
           DO jj = 1,jpj
             DO ji = 1,jpi
C
C     DOC ammonification. Depends on depth, phytoplankton biomass
C     and a limitation term which is supposed to be a parameterization
C     of the bacterial activity. 
C     ----------------------------------------------------------------
C
         remik = xremik*zstep/1E-6*xlimbac(ji,jj,jk)
     &     *zdepbac(ji,jj,jk)
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif
         remik=max(remik,5.5E-4*zstep)
C
C     Ammonification in oxic waters with oxygen consumption
C     -----------------------------------------------------
C
         olimi(ji,jj,jk)=min((trn(ji,jj,jk,jpoxy)-rtrn)/o2ut,
     &     remik*(1.-nitrfac(ji,jj,jk))*trn(ji,jj,jk,jpdoc)) 
C
C     Ammonification in suboxic waters with denitrification
C     -------------------------------------------------------
C
         denitr(ji,jj,jk)=min((trn(ji,jj,jk,jpno3)-rtrn)/rdenit,
     &     remik*nitrfac(ji,jj,jk)*trn(ji,jj,jk,jpdoc))
             END DO
           END DO
         END DO
C
         olimi(:,:,:)=max(0.,olimi(:,:,:))
         denitr(:,:,:)=max(0.,denitr(:,:,:))
C
         DO jk = 1,jpkm1
           DO jj = 1,jpj
             DO ji = 1,jpi
C
C    NH4 nitrification to NO3. Ceased for oxygen concentrations
C    below 2 umol/L. Inhibited at strong light 
C    ----------------------------------------------------------
C
         onitr(ji,jj,jk)=nitrif*zstep*trn(ji,jj,jk,jpnh4)/(1.
     &     +emoy(ji,jj,jk))*(1.-nitrfac(ji,jj,jk))
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif
             END DO
           END DO
         END DO

         DO jk = 1,jpkm1
           DO jj = 1,jpj
             DO ji = 1,jpi
C
C    Bacterial uptake of iron. No iron is available in DOC. So
C    Bacteries are obliged to take up iron from the water. Some
C    studies (especially at Papa) have shown this uptake to be
C    significant
C    ----------------------------------------------------------
C
         xbactfer(ji,jj,jk)=15E-6*rfact2*4.*0.4*prmax(ji,jj,jk)
     &     *(xlimphy(ji,jj,jk)*zdepbac(ji,jj,jk))**2
     &     /(xkgraz2+zdepbac(ji,jj,jk))
     &     *(0.5+sign(0.5,trn(ji,jj,jk,jpfer)-2E-11))
C
             END DO
           END DO
         END DO
C
         DO jk = 1,jpkm1
           DO jj = 1,jpj
             DO ji = 1,jpi
C
C    POC disaggregation by turbulence and bacterial activity. 
C    -------------------------------------------------------------
C
         remip=xremip*zstep*tgfunc(ji,jj,jk)*(1.-0.5*nitrfac(ji,jj,jk))
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif
C
C    POC disaggregation rate is reduced in anoxic zone as shown by
C    sediment traps data. In oxic area, the exponent of the martin's
C    law is around -0.87. In anoxic zone, it is around -0.35. This
C    means a disaggregation constant about 0.5 the value in oxic zones
C    -----------------------------------------------------------------
C
         orem(ji,jj,jk)=remip*trn(ji,jj,jk,jppoc)
         orem2(ji,jj,jk)=remip*trn(ji,jj,jk,jpgoc)
         ofer(ji,jj,jk)=remip*trn(ji,jj,jk,jpsfe)
         ofer2(ji,jj,jk)=remip*trn(ji,jj,jk,jpbfe)
C
             END DO
           END DO
         END DO

         DO jk = 1,jpkm1
           DO jj = 1,jpj
             DO ji = 1,jpi
C
C     Remineralization rate of BSi depedant on T and saturation
C     ---------------------------------------------------------
C
         zsatur=(sio3eq(ji,jj,jk)-trn(ji,jj,jk,jpsil))/
     &     (sio3eq(ji,jj,jk)+rtrn)
         zsatur=max(rtrn,zsatur)
         zsatur2=zsatur*(1.+tn(ji,jj,jk)/400.)**4
         znusil=0.225*(1.+tn(ji,jj,jk)/15.)*zsatur+0.775*zsatur2**9

         siremin=xsirem*zstep*znusil
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif
C
         osil(ji,jj,jk)=siremin*trn(ji,jj,jk,jpdsi)
             END DO
           END DO
         END DO
C
         fesatur(:,:,:)=0.6E-9
C
         DO jk = 1,jpkm1
           DO jj = 1,jpj
             DO ji = 1,jpi
C
C     scavenging rate of iron. this scavenging rate depends on the
C     load in particles on which they are adsorbed. The
C     parameterization has been taken from studies on Th
C     ------------------------------------------------------------
C
         xkeq=fekeq(ji,jj,jk)
         xfeequi=(-(1.+fesatur(ji,jj,jk)*xkeq-xkeq*trn(ji,jj,jk,jpfer))+
     &     sqrt((1.+fesatur(ji,jj,jk)*xkeq-xkeq*trn(ji,jj,jk,jpfer))**2
     &     +4.*trn(ji,jj,jk,jpfer)*xkeq))/(2.*xkeq)

         xlam1b=3E-5+xlam1*(trn(ji,jj,jk,jppoc)
     &     +trn(ji,jj,jk,jpgoc)+trn(ji,jj,jk,jpcal)+
     &      trn(ji,jj,jk,jpdsi))*1E6

         xscave(ji,jj,jk)=xfeequi*xlam1b*zstep
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif
C
C  Increased scavenging for very high iron concentrations
C  found near the coasts due to increased lithogenic particles
C  and let's say it unknown processes (precipitation, ...)
C  -----------------------------------------------------------
C
         zlamfac=max(0.,(gphit(ji,jj)+55.)/30.)
         zlamfac=min(1.,zlamfac)
         xlam1b=(80.*(trn(ji,jj,jk,jpdoc)+40E-6)+698.
     &    *trn(ji,jj,jk,jppoc)+1.05E4*trn(ji,jj,jk,jpgoc))
     &    *zdiss(ji,jj,jk)+1E-5*(1.-zlamfac)+xlam1*max(0.,
     &    (trn(ji,jj,jk,jpfer)*1E9-1.))

         xaggdfe(ji,jj,jk)=xlam1b*zstep*0.76*(trn(ji,jj,jk,jpfer)
     &     -xfeequi)
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif

C
             END DO
           END DO
         END DO
C
#endif
      RETURN
      END
