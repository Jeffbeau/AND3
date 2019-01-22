MODULE trcdtr
   !!=======================================================================================
   !!
   !!                       *** MODULE trcdtr ***
   !!
   !!  Computes or READ initial DATA for passive tracer
   !!
   !!=======================================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trcdtr.F90,v 1.5 2006/04/10 15:40:28 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !! * Modules used
      !! ==============
      USE oce_trc
      USE trc
      USE sms
      USE trcdta
      USE lib_mpp

      IMPLICIT NONE
      PRIVATE
  !! * Accessibility
      PUBLIC trc_dtr

CONTAINS

#if defined key_passivetrc

SUBROUTINE trc_dtr
!!---------------------------------------------------------------------
!!
!!                       ROUTINE trci_dtr
!!                     ******************
!!  PURPOSE :
!!  ---------
!!     computes or READ initial DATA for passive tracer
!!   -----
!!      COMMON
!!            /comdom/          : domain PARAMETER
!!            /comcoo/          : orthogonal curvilinear coordinates
!!                                and scale factors
!!            /comask/          : masks, bathymetry
!!   OUTPUT :
!!   ------
!!      COMMON
!!            /cottrc/          : passive tracer field now and before
!!
!!
!!   History:
!!   --------
!!      original  : 96-11
!!      additions : 99-9
!!                : 00-12 (O. Aumont, E. Kestenare) add for POC in sediments 
!!                         add for POC in sediments  
!!    03/05  O. Aumont and A. El Moussaoui  F90 
!!----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! local declarations
!! ================== 
      INTEGER :: ji,jj,jk,jn 
#if defined key_trc_pisces
      REAL(wp) :: alka0,oxyg0,calc0,bioma0,    &
                  silic1,po4,no3,caralk,bicarb
#endif
!!---------------------------------------------------------------------
!!  OPA.9 
!!---------------------------------------------------------------------
!! 0. initialisations
!! ------------------

      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) ' *** trcdtr initialisation for '
      IF(lwp) WRITE(numout,*) '     passive tracers'
      IF(lwp) WRITE(numout,*) ' '


#if defined key_cfc
      trn(:,:,:,:)=0.0
#elif defined key_trc_pisces

      sco2   = 2.3e-3
      alka0  = 2.39e-3
      oxyg0  = 1.8e-4
      po4    = 2.165e-6/po4r
      bioma0 = 1.e-8
      silic1 = 91.51e-6
      calc0  = 1.e-6
      no3    = 30.88E-6*7.6

      trn(:,:,:,jpdic) = sco2
      trn(:,:,:,jptal) = alka0
      trn(:,:,:,jpoxy) = oxyg0
      trn(:,:,:,jppo4) = po4
      trn(:,:,:,jppoc) = bioma0
      trn(:,:,:,jpsil) = silic1
      trn(:,:,:,jpcal) = calc0
      trn(:,:,:,jpphy) = bioma0
      trn(:,:,:,jpzoo) = bioma0
      trn(:,:,:,jpdoc) = bioma0
      trn(:,:,:,jpdia) = bioma0
      trn(:,:,:,jpmes) = bioma0
      trn(:,:,:,jpbsi) = bioma0*0.15
      trn(:,:,:,jpfer) = 0.6E-9
      trn(:,:,:,jpbfe) = bioma0*5E-6
      trn(:,:,:,jpgoc) = bioma0
      trn(:,:,:,jpsfe) = bioma0*5.E-6
      trn(:,:,:,jpdfe) = bioma0*5.E-6
      trn(:,:,:,jpnfe) = bioma0*5.E-6
      trn(:,:,:,jpdsi) = bioma0*5.E-6
      trn(:,:,:,jpnch) = bioma0*12./55.
      trn(:,:,:,jpdch) = bioma0*12./55.
      trn(:,:,:,jpno3) = no3 
      trn(:,:,:,jpnh4) = bioma0

!!  Initialization of chemical variables of the carbon cycle
!!  --------------------------------------------------------

      DO jk = 1,jpk
        DO jj = 1,jpj
          DO ji = 1,jpi
              caralk = trn(ji,jj,jk,jptal)-         & 
                      borat(ji,jj,jk)/(1.+1.E-8/(rtrn+akb3(ji,jj,jk)))
              co3(ji,jj,jk)=(caralk-trn(ji,jj,jk,jpdic))*tmask(ji,jj,jk)       & 
                     +(1.-tmask(ji,jj,jk))*.5e-3
              bicarb = (2.*trn(ji,jj,jk,jpdic)-caralk)
              hi(ji,jj,jk) = (ak23(ji,jj,jk)*bicarb/co3(ji,jj,jk))             &  
                *tmask(ji,jj,jk)+(1.-tmask(ji,jj,jk))*1.e-9
          ENDDO
        ENDDO
      ENDDO

      h2co3(:,:) = 1.e-5

!!  initialize the half saturation constant for silicate
!!  ----------------------------------------------------

      xksi(:,:)=2.E-6

!! initialize Silicate and Calcite and POC in sediments
!! ---------------------------------------------------
      sedcal(:,:) = 0.
      sedsil(:,:) = 0.
      sedpoc(:,:) = 0.

      IF(lwp) WRITE(numout,*) 'Initialization of PISCES tracers done'
      IF(lwp) WRITE(numout,*) ' '

#elif defined key_trc_lobster1 && ( defined key_eel_r6 || defined key_eel_r2 )
! analytical initialisation used in Levy et al. (2001)
      
      DO jk=1,7
        trn(:,:,jk,jpdet)=0.016*tmask(:,:,jk)
        trn(:,:,jk,jpzoo)=0.018*tmask(:,:,jk)
        trn(:,:,jk,jpphy)=0.036*tmask(:,:,jk)
        trn(:,:,jk,jpno3)=1.e-5*tmask(:,:,jk)
        trn(:,:,jk,jpnh4)=0.0005*tmask(:,:,jk)
        trn(:,:,jk,jpdom)=0.017*tmask(:,:,jk)
      END DO

      trn(:,:,8,jpdet)=0.020*tmask(:,:,1)
      trn(:,:,8,jpzoo)=0.027*tmask(:,:,1)
      trn(:,:,8,jpphy)=0.041*tmask(:,:,1)
      trn(:,:,8,jpno3)=0.00022*tmask(:,:,1)
      trn(:,:,8,jpnh4)=0.0033*tmask(:,:,1)
      trn(:,:,8,jpdom)=0.021*tmask(:,:,1)

      trn(:,:,9,jpdet)=0.0556*tmask(:,:,1)
      trn(:,:,9,jpzoo)=0.123*tmask(:,:,1)
      trn(:,:,9,jpphy)=0.122*tmask(:,:,1)
      trn(:,:,9,jpno3)=0.028*tmask(:,:,1)
      trn(:,:,9,jpnh4)=0.024*tmask(:,:,1)
      trn(:,:,9,jpdom)=0.06*tmask(:,:,1)

      trn(:,:,10,jpdet)=0.025*tmask(:,:,1)
      trn(:,:,10,jpzoo)=0.016*tmask(:,:,1)
      trn(:,:,10,jpphy)=0.029*tmask(:,:,1)
      trn(:,:,10,jpno3)=2.462*tmask(:,:,1)
      trn(:,:,10,jpnh4)=0.04*tmask(:,:,1)
      trn(:,:,10,jpdom)=0.022*tmask(:,:,1)

      trn(:,:,11,jpdet)=0.0057*tmask(:,:,1)
      trn(:,:,11,jpzoo)=0.0005*tmask(:,:,1)
      trn(:,:,11,jpphy)=0.0006*tmask(:,:,1)
      trn(:,:,11,jpno3)=3.336*tmask(:,:,1)
      trn(:,:,11,jpnh4)=0.005*tmask(:,:,1)
      trn(:,:,11,jpdom)=0.004*tmask(:,:,1)

      trn(:,:,12,jpdet)=0.002*tmask(:,:,1)
      trn(:,:,12,jpzoo)=1.e-6*tmask(:,:,1)
      trn(:,:,12,jpphy)=5.e-6*tmask(:,:,1)
      trn(:,:,12,jpno3)=4.24*tmask(:,:,1)
      trn(:,:,12,jpnh4)=0.001*tmask(:,:,1)
      trn(:,:,12,jpdom)=3.e-5*tmask(:,:,1)

      DO jk=13,jpk
        trn(:,:,jk,jpdet)=0.0
        trn(:,:,jk,jpzoo)=0.0
        trn(:,:,jk,jpphy)=0.0
        trn(:,:,jk,jpnh4)=0.0
        trn(:,:,jk,jpdom)=0.0
      END DO

      trn(:,:,13,jpno3)=5.31*tmask(:,:,13)
      trn(:,:,14,jpno3)=6.73*tmask(:,:,14)
      trn(:,:,15,jpno3)=8.32*tmask(:,:,15)
      trn(:,:,16,jpno3)=10.13*tmask(:,:,16)
      trn(:,:,17,jpno3)=11.95*tmask(:,:,17)
      trn(:,:,18,jpno3)=13.57*tmask(:,:,18)
      trn(:,:,19,jpno3)=15.08*tmask(:,:,19)
      trn(:,:,20,jpno3)=16.41*tmask(:,:,20)
      trn(:,:,21,jpno3)=17.47*tmask(:,:,21)
      trn(:,:,22,jpno3)=18.29*tmask(:,:,22)
      trn(:,:,23,jpno3)=18.88*tmask(:,:,23)
      trn(:,:,24,jpno3)=19.30*tmask(:,:,24)
      trn(:,:,25,jpno3)=19.68*tmask(:,:,25)
      trn(:,:,26,jpno3)=19.91*tmask(:,:,26)
      trn(:,:,27,jpno3)=19.99*tmask(:,:,27)
      trn(:,:,28,jpno3)=20.01*tmask(:,:,28)
      trn(:,:,29,jpno3)=20.01*tmask(:,:,29)
      trn(:,:,30,jpno3)=20.01*tmask(:,:,30)

#elif defined key_trc_lobster1 && defined key_gyre
! init NO3=f(density) by asklod AS Kremeur 2005-07
      trn(:,:,:,jpdet)=0.1*tmask(:,:,:)
      trn(:,:,:,jpzoo)=0.1*tmask(:,:,:)
      trn(:,:,:,jpnh4)=0.1*tmask(:,:,:)
      trn(:,:,:,jpphy)=0.1*tmask(:,:,:)
      trn(:,:,:,jpdom)=1.*tmask(:,:,:)
      DO  jk=1,jpk
         DO  jj=1,jpj
            DO  ji=1,jpi
               IF (rhd(ji,jj,jk).LE.24.5e-3) THEN
                  trn(ji,jj,jk,jpno3)=2.*tmask(ji,jj,jk)
               ELSE
                  trn(ji,jj,jk,jpno3)=(15.55*(rhd(ji,jj,jk)*1000)-380.11)*tmask(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO

#else
 
!! general case
      do jn = 1, jptra
         trn(:,:,:,jn)=0.1*tmask(:,:,:)
      enddo

#endif

#if defined key_trc_dta
!!   Initialization of tracer from a file
!!   that may also be used for damping
      CALL dta_trc( nit000 )
      DO  jk = 1, jptra
        IF( lutini(jk) ) THEN 
!! initialisation from file
           trn(:,:,:,jk) = trdta(:,:,:,jk)*tmask(:,:,:)
        ENDIF
      END DO
#endif

!! before field :
!! -------------
      trb(:,:,:,:) = trn(:,:,:,:)

#if defined key_trc_lobster1
!!  initialize the POC in sediments

      sedpocb(:,:) = 0.
      sedpocn(:,:) = 0.
      sedpoca(:,:) = 0.
#endif
      
 END SUBROUTINE trc_dtr 

#else

SUBROUTINE  trc_dtr 
!!======================
   !! no passive tracers
!!======================
END SUBROUTINE  trc_dtr
#endif

END MODULE trcdtr
