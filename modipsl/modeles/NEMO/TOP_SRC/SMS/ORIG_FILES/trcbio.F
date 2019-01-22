
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/trcbio.F,v 1.7 2006/04/10 15:39:56 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
      SUBROUTINE trcbio(kt)
#if defined key_passivetrc && defined key_trc_lobster1
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcbio
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     compute the now trend due to biogeochemical processes
CCC     and add it to the general trend of passive tracers equations.
CCC
CCC	Three options:
CCC  	  Default option  : no biological trend
CCC       IF 'key_trc_lobster1' : LOBSTER1 model
CCC
CC   METHOD :
CC   -------
CC      each now biological flux is calculated  in FUNCTION of now
CC      concentrations of tracers.
CC      depending on the tracer, these fluxes are sources or sinks.
CC      the total of the sources and sinks for each tracer
CC      is added to the general trend.
CC
CC		tra = tra + zf...tra - zftra...
CC                             |         |
CC                             |         |
CC                          source      sink
CC
CC
CC      IF 'key_trc_diabio' key is activated, the biogeochemical
CC	trends for passive tracers are saved for futher diagnostics.
CC
CC      multitasked on vertical slab (jj-loop)
CC
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC                                depths
CC            /cottrp/          : present and next fields for passive
CC                              : tracers
CC            /comtsk/          : multitasking
CC            /comtke/          : emin, en()
CC            /cotbio/          : biological parameters
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cottrp/ tra      : general tracer trend increased by the
CC                                now horizontal tracer advection trend
CC            /cottbd/ trbio    : now horizontal tracer advection trend
CC                                (IF 'key_trc_diabio' is activated)
CC
CC   WORKSPACE :
CC   ---------
CC      local
CC               zdet,zzoo,zphy,znh4,zno3,zdom    : now concentrations
CC               zlt,zlno3,zlnh4,zle              : limitation terms for phyto
CC               zfno3phy and so on..             : fluxes between bio boxes
CC               zphya,zzooa,zdeta, ...           : after bio trends
CC               zppz, zpdz, zpppz, zppdz, zfood  : preferences terms
CC               zfilpz, zfilpd                   : filtration terms
CC      COMMON
CC
CC   EXTERNAL :                   no
CC   --------
CC
CC   REFERENCES :                 no
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC       original : 99-07 (M. Levy)
CC                  00-12 (E. Kestenare): assign a parameter 
CC                                        to name individual tracers
CC                  01-03 (M. Levy) LNO3 + dia2d
CC----------------------------------------------------------------------
CC----------------------------------------------------------------------
      USE oce_trc
      USE trp_trc
      USE sms
      USE lbclnk
      IMPLICIT NONE
CC local declarations
CC ==================
      INTEGER kt
      INTEGER ji,jj,jk,jn
      REAL ztot(jpi), ze3t(jpk)
      REAL zdet,zzoo,zphy,zno3,znh4,zdom,zlno3,zlnh4,zle,zlt
      REAL zno3phy, znh4phy, zphynh4, zphydom, zphydet, zphyzoo, zdetzoo
     $    ,zzoonh4, zzoodom, zzoodet, zdetnh4, zdetdom, znh4no3, zdomnh4
     $    ,zppz,zpdz,zpppz,zppdz,zfood,zfilpz,zfildz,zphya,zzooa,zno3a
     $    ,znh4a,zdeta,zdoma, ztra, zzoobod, zboddet, zdomaju

CC----------------------------------------------------------------------
CC statement functions
CC ===================
! CDIR$ NOLIST
#include "domzgr_substitute.h90"
! CDIR$ LIST
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (07/99)
CCC---------------------------------------------------------------------
C   | --------------|
C   | LOBSTER1 MODEL| 
C   | --------------|

#if defined key_trc_diaadd
C convert fluxes in per day
      DO jk=1,jpkbm1
        ze3t(jk)=e3t(jk)*86400.
      END DO 
      DO jk=jpkb,jpk
        ze3t(jk)=0.
      END DO 
#endif
C
C vertical slab
C =============
C
      DO 1000 jj = 2,jpjm1
C
C 1. biological level
C ===================
C
        DO ji = 2,jpim1
          fbod(ji,jj)=0.
#if defined key_trc_diaadd
          DO jn=1,jpdia2d
            trc2d(ji,jj,jn)=0.          
          END DO 
#endif
        END DO 

        DO jk=1,jpkbm1
          DO ji = 2,jpim1
C
C
C 1.1 trophic variables( det, zoo, phy, no3, nh4, dom)
C ---------------------------------------------------
C
C negative trophic variables DO not contribute to the fluxes
C
            zdet = max(0.,trn(ji,jj,jk,jpdet))
            zzoo = max(0.,trn(ji,jj,jk,jpzoo))
            zphy = max(0.,trn(ji,jj,jk,jpphy))
            zno3 = max(0.,trn(ji,jj,jk,jpno3))
            znh4 = max(0.,trn(ji,jj,jk,jpnh4))
            zdom = max(0.,trn(ji,jj,jk,jpdom))
C
C
C 1.2  Limitations
C ----------------
C
            zlt = 1.
            zle = 1. - exp( -xpar(ji,jj,jk)/aki/zlt)
C psinut,akno3,aknh4 added by asklod AS Kremeur 2005-03
            zlno3 = zno3* exp(-psinut*znh4) / (akno3+zno3)
            zlnh4 = znh4 / (znh4+aknh4) 

C
C
C 1.3 sinks and sources
C ---------------------
C
C
C 1. phytoplankton production and exsudation
C
            zno3phy = tmumax * zle * zlt * zlno3 * zphy
            znh4phy = tmumax * zle * zlt * zlnh4 * zphy

C fphylab added by asklod AS Kremeur 2005-03
            zphydom = rgamma * (1 - fphylab) * (zno3phy + znh4phy)
            zphynh4 = rgamma * fphylab * (zno3phy + znh4phy)

C
C 2. zooplankton production
C
C preferences
C
            zppz = rppz
            zpdz = 1. - rppz
            zpppz = ( zppz * zphy ) /
     $          ( ( zppz * zphy + zpdz * zdet ) + 1.e-13 )
            zppdz = ( zpdz * zdet ) /
     $          ( ( zppz * zphy + zpdz * zdet ) + 1.e-13 )
            zfood = zpppz * zphy + zppdz * zdet
C
C filtration
C
            zfilpz = taus * zpppz / (aks + zfood)
            zfildz = taus * zppdz / (aks + zfood)
C
C grazing
C
            zphyzoo = zfilpz * zphy * zzoo
            zdetzoo = zfildz * zdet * zzoo
C
C 3. fecal pellets production
C
            zzoodet = rpnaz * zphyzoo + rdnaz * zdetzoo
C
C 4. zooplankton liquide excretion
C
            zzoonh4 = tauzn * fzoolab * zzoo 
            zzoodom = tauzn * (1 - fzoolab) * zzoo
C
C 5. mortality
C
C phytoplankton mortality 
C
            zphydet = tmminp * zphy
C
C
C zooplankton mortality
c closure : flux fbod is redistributed below level jpkbio
C
            zzoobod = tmminz * zzoo * zzoo
            fbod(ji,jj) = fbod(ji,jj) 
     $                 + (1-fdbod) * zzoobod * fse3t(ji,jj,jk)
            zboddet = fdbod * zzoobod
C
C
C 6. detritus and dom breakdown
C
C
            zdetnh4 = taudn * fdetlab * zdet
            zdetdom = taudn * (1 - fdetlab) * zdet 

            zdomnh4 = taudomn * zdom
C
C flux added to express how the excess of nitrogen from
C PHY, ZOO and DET to DOM goes directly to NH4 (flux of ajustment)
            zdomaju = (1 - redf/reddom) * (zphydom + zzoodom + zdetdom)
C
C 7. Nitrification
C
            znh4no3 = taunn * znh4
C
C
C
C 1.4 determination of trends
C ---------------------------
C
C total trend for each biological tracer
C
            zphya =   zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo
     $          - zphydet
            zzooa =   zphyzoo + zdetzoo - zzoodet - zzoodom - zzoonh4
     $          - zzoobod
            zno3a = - zno3phy + znh4no3
            znh4a = - znh4phy - znh4no3 + zphynh4 + zzoonh4 + zdomnh4
     $          + zdetnh4 + zdomaju
            zdeta = zphydet + zzoodet  - zdetzoo - zdetnh4 - zdetdom +
     $          zboddet
            zdoma = zphydom + zzoodom + zdetdom - zdomnh4 - zdomaju
C
#if defined key_trc_diabio
            trbio(ji,jj,jk,1) = zno3phy
            trbio(ji,jj,jk,2) = znh4phy
            trbio(ji,jj,jk,3) = zphynh4
            trbio(ji,jj,jk,4) = zphydom
            trbio(ji,jj,jk,5) = zphyzoo
            trbio(ji,jj,jk,6) = zphydet
            trbio(ji,jj,jk,7) = zdetzoo
            trbio(ji,jj,jk,9) = zzoodet
            trbio(ji,jj,jk,10) = zzoobod
            trbio(ji,jj,jk,11) = zzoonh4
            trbio(ji,jj,jk,12) = zzoodom
            trbio(ji,jj,jk,13) = znh4no3
            trbio(ji,jj,jk,14) = zdomnh4
            trbio(ji,jj,jk,15) = zdetnh4
#endif
#if defined key_trc_diaadd
            trc2d(ji,jj,1)=trc2d(ji,jj,1)+zno3phy*ze3t(jk)          
            trc2d(ji,jj,2)=trc2d(ji,jj,2)+znh4phy*ze3t(jk)
            trc2d(ji,jj,3)=trc2d(ji,jj,3)+zphydom*ze3t(jk)
            trc2d(ji,jj,4)=trc2d(ji,jj,4)+zphynh4*ze3t(jk)
            trc2d(ji,jj,5)=trc2d(ji,jj,5)+zphyzoo*ze3t(jk)
            trc2d(ji,jj,6)=trc2d(ji,jj,6)+zphydet*ze3t(jk)
            trc2d(ji,jj,7)=trc2d(ji,jj,7)+zdetzoo*ze3t(jk)
c trend number 8 is in trcsed.F            
            trc2d(ji,jj,9)=trc2d(ji,jj,9)+zzoodet*ze3t(jk)
            trc2d(ji,jj,10)=trc2d(ji,jj,10)+zzoobod*ze3t(jk)
            trc2d(ji,jj,11)=trc2d(ji,jj,11)+zzoonh4*ze3t(jk)
            trc2d(ji,jj,12)=trc2d(ji,jj,12)+zzoodom*ze3t(jk)
            trc2d(ji,jj,13)=trc2d(ji,jj,13)+znh4no3*ze3t(jk)
            trc2d(ji,jj,14)=trc2d(ji,jj,14)+zdomnh4*ze3t(jk)
            trc2d(ji,jj,15)=trc2d(ji,jj,15)+zdetnh4*ze3t(jk)
             
            trc2d(ji,jj,16)=trc2d(ji,jj,16)+(zno3phy+znh4phy-zphynh4
     $          -zphydom-zphyzoo-zphydet)*ze3t(jk)
            trc2d(ji,jj,17)=trc2d(ji,jj,17)+(zphyzoo+zdetzoo-zzoodet
     $          -zzoobod-zzoonh4-zzoodom) *ze3t(jk)
            trc2d(ji,jj,18)=trc2d(ji,jj,18)+zdetdom*ze3t(jk)
c trend number 19 is in trcexp.F
            trc3d(ji,jj,jk,1)= zno3phy *86400     
            trc3d(ji,jj,jk,2)= znh4phy *86400     
            trc3d(ji,jj,jk,3)= znh4no3 *86400     
#endif
C
C tracer flux at totox-point added to the general trend
C
            tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + zdeta
            tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + zzooa
            tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zphya
            tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zno3a
            tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + znh4a
            tra(ji,jj,jk,jpdom) = tra(ji,jj,jk,jpdom) + zdoma
C
          END DO
        END DO
C
C 2. under biological level
C =========================
C
        DO jk = jpkb,jpk
C
C 2.1 compute the remineralisation of all quantities towards nitrate 
C ------------------------------------------------------------------
C
          DO ji = 2,jpim1
C
C 2.1.1 trophic variables( det, zoo, phy, no3, nh4, dom)
C -----------------------------------------------------
C
C negative trophic variables DO not contribute to the fluxes
C
            zdet = max(0.,trn(ji,jj,jk,jpdet))
            zzoo = max(0.,trn(ji,jj,jk,jpzoo))
            zphy = max(0.,trn(ji,jj,jk,jpphy))
            zno3 = max(0.,trn(ji,jj,jk,jpno3))
            znh4 = max(0.,trn(ji,jj,jk,jpnh4))
            zdom = max(0.,trn(ji,jj,jk,jpdom))
CC
CC 2.1.2  Limitations
CC ----------------
CC
            zlt = 0.
            zle = 0.
            zlno3 = 0.
            zlnh4 = 0.
CC
CC
CC 2.1.3 sinks and sources
CC ---------------------
CC
CC
CC 1. phytoplankton production and exsudation
CC
            zno3phy = 0.
            znh4phy = 0.
C
            zphydom = 0.
            zphynh4 = 0.
CC
CC 2. zooplankton production
CC
CC grazing
CC
            zphyzoo = 0. 
            zdetzoo = 0.
CC
CC 3. fecal pellets production
CC
            zzoodet = 0.
CC
CC 4. zooplankton liquide excretion
CC
            zzoonh4 = tauzn * fzoolab * zzoo 
            zzoodom = tauzn * (1 - fzoolab) * zzoo
CC
CC 5. mortality
CC
CC phytoplankton mortality 
CC
            zphydet = tmminp * zphy
CC
CC
CC zooplankton mortality
Cc closure : flux fbod is redistributed below level jpkbio
CC
            zzoobod = 0.
            zboddet = 0. 
CC
CC
CC 6. detritus and dom breakdown
CC
            zdetnh4 = taudn * fdetlab * zdet
            zdetdom = taudn * (1 - fdetlab) * zdet 
C
            zdomnh4 = taudomn * zdom
            zdomaju = (1 - redf/reddom) * (zphydom + zzoodom + zdetdom)
CC
CC 7. Nitrification
CC
            znh4no3 = taunn * znh4
CC
CC
CC 2.1.4 determination of trends
CC ---------------------------
CC
CC total trend for each biological tracer
CC
            zphya =   zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo
     $          - zphydet
            zzooa =   zphyzoo + zdetzoo - zzoodet - zzoodom - zzoonh4
     $          - zzoobod
            zno3a = - zno3phy + znh4no3
            znh4a = - znh4phy - znh4no3 + zphynh4 + zzoonh4 + zdomnh4
     $          + zdetnh4 + zdomaju
            zdeta = zphydet + zzoodet  - zdetzoo - zdetnh4 - zdetdom +
     $          zboddet
            zdoma = zphydom + zzoodom + zdetdom - zdomnh4 - zdomaju
CC
#if defined key_trc_diabio
            trbio(ji,jj,jk,1) = zno3phy
            trbio(ji,jj,jk,2) = znh4phy
            trbio(ji,jj,jk,3) = zphynh4
            trbio(ji,jj,jk,4) = zphydom
            trbio(ji,jj,jk,5) = zphyzoo
            trbio(ji,jj,jk,6) = zphydet
            trbio(ji,jj,jk,7) = zdetzoo
            trbio(ji,jj,jk,9) = zzoodet
            trbio(ji,jj,jk,10) = zzoobod
            trbio(ji,jj,jk,11) = zzoonh4
            trbio(ji,jj,jk,12) = zzoodom
            trbio(ji,jj,jk,13) = znh4no3
            trbio(ji,jj,jk,14) = zdomnh4
            trbio(ji,jj,jk,15) = zdetnh4
#endif
#if defined key_trc_diaadd
            trc2d(ji,jj,1)=trc2d(ji,jj,1)+zno3phy*ze3t(jk)          
            trc2d(ji,jj,2)=trc2d(ji,jj,2)+znh4phy*ze3t(jk)
            trc2d(ji,jj,3)=trc2d(ji,jj,3)+zphydom*ze3t(jk)
            trc2d(ji,jj,4)=trc2d(ji,jj,4)+zphynh4*ze3t(jk)
            trc2d(ji,jj,5)=trc2d(ji,jj,5)+zphyzoo*ze3t(jk)
            trc2d(ji,jj,6)=trc2d(ji,jj,6)+zphydet*ze3t(jk)
            trc2d(ji,jj,7)=trc2d(ji,jj,7)+zdetzoo*ze3t(jk)
Cc trend number 8 is in trcsed.F            
            trc2d(ji,jj,9)=trc2d(ji,jj,9)+zzoodet*ze3t(jk)
            trc2d(ji,jj,10)=trc2d(ji,jj,10)+zzoobod*ze3t(jk)
            trc2d(ji,jj,11)=trc2d(ji,jj,11)+zzoonh4*ze3t(jk)
            trc2d(ji,jj,12)=trc2d(ji,jj,12)+zzoodom*ze3t(jk)
            trc2d(ji,jj,13)=trc2d(ji,jj,13)+znh4no3*ze3t(jk)
            trc2d(ji,jj,14)=trc2d(ji,jj,14)+zdomnh4*ze3t(jk)
            trc2d(ji,jj,15)=trc2d(ji,jj,15)+zdetnh4*ze3t(jk)
             
            trc2d(ji,jj,16)=trc2d(ji,jj,16)+(zno3phy+znh4phy-zphynh4
     $          -zphydom-zphyzoo-zphydet)*ze3t(jk)
            trc2d(ji,jj,17)=trc2d(ji,jj,17)+(zphyzoo+zdetzoo-zzoodet
     $          -zzoobod-zzoonh4-zzoodom) *ze3t(jk)
            trc2d(ji,jj,18)=trc2d(ji,jj,18)+zdetdom*ze3t(jk)

            trc3d(ji,jj,jk,1)= zno3phy *86400     
            trc3d(ji,jj,jk,2)= znh4phy *86400     
            trc3d(ji,jj,jk,3)= znh4no3 *86400     
#endif
CC
CC tracer flux at totox-point added to the general trend
CC
            tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + zdeta
            tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + zzooa
            tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zphya
            tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zno3a
            tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + znh4a
            tra(ji,jj,jk,jpdom) = tra(ji,jj,jk,jpdom) + zdoma
CC
          END DO
        END DO




c$$$        DO jk = jpkb,jpk
c$$$C
c$$$C 2.1 Old way to compute the remineralisation : asklod AS Kremeur (before 2005-03)
c$$$C ------------------------------------------------------------------
c$$$C
c$$$          DO ji=2,jpim1
c$$$            ztot(ji) = 0.
c$$$          END DO 
c$$$          DO jn=1,jptra
c$$$            IF (ctrcnm(jn).NE.'NO3') THEN 
c$$$                DO ji=2,jpim1
c$$$                  ztra = remdmp(jk,jn) * trn(ji,jj,jk,jn) 
c$$$                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) - ztra
c$$$                  ztot(ji) = ztot(ji) + ztra
c$$$                END DO 
c$$$            ENDIF
c$$$          END DO 
c$$$          DO jn=1,jptra
c$$$            IF (ctrcnm(jn).EQ.'NO3') THEN 
c$$$                DO ji=2,jpim1
c$$$                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztot(ji)
c$$$                END DO
c$$$#if defined key_trc_diabio
c$$$                trbio(ji,jj,jk,1)=ztot(ji)
c$$$#endif 
c$$$            ENDIF
c$$$          END DO
c$$$        END DO 

C
C
C END of slab
C ===========
C
 1000 CONTINUE

#if defined key_trc_diaadd

C Lateral boundary conditions on trc2d
      DO jn=1,jpdia2d
          CALL lbc_lnk(trc2d(:,:,jn),'T',1. )
      END DO 

C Lateral boundary conditions on trc3d
      DO jn=1,jpdia3d
          CALL lbc_lnk(trc3d(:,:,1,jn),'T',1. )
      END DO 

#endif

#if defined key_trc_diabio
C Lateral boundary conditions on trcbio
      DO jn=1,jpdiabio
          CALL lbc_lnk(trbio(:,:,1,jn),'T',1. )
      END DO 
#endif

#  else
C
C    no biological model
C
#  endif

C
C
      RETURN
      END
