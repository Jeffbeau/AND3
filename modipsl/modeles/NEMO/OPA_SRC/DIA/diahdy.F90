MODULE diahdy
   !!======================================================================
   !!                       ***  MODULE  diahdy  ***
   !! Ocean diagnostics : computation the dynamical heigh
   !!======================================================================
#if   defined key_diahdy   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_diahdy' :                          dynamical heigh diagnostics
   !!----------------------------------------------------------------------
   !!   dia_hdy      : dynamical heigh computation
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dia_hdy     ! called in step.F90 module

   !! * Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_diahdy = .TRUE.   !: dynamical heigh flag

   !! * Module variables
   REAL(wp), DIMENSION(jpk) ::   &
      rhosp         ! ???

   REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
      hdy           ! dynamical heigh

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diahdy.F90,v 1.3 2005/03/27 18:34:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE dia_hdy ( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hdy  ***
      !!
      !! ** Purpose :   Computes the dynamical heigh
      !!
      !! ** Method  : Millero + Poisson
      !!
      !! References : 
      !!	A. E. Gill, atmosphere-ocean dynamics 7.7 pp 215
      !!
      !! History :
      !!        !  9x-xx (P. Delecluse, C. Perigaud)  Original code
      !!        !  93-10  (C. Perigaud)  a trapezoidal vertical integration 
      !!                                 consistent WITH the code
      !!        !  93-12  (G. Madec M. Imbard)
      !!        !  96-03  (N. Ferry)  integration at t-points
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER :: ji, jj, jk
      INTEGER :: ihdsup, ik

      REAL(wp) :: zgdsup, za, zb, zciint, zfacto, zhd
      REAL(wp) :: zp, zh, zt, zs, zxk, zq, zsr, zr1, zr2, zr3, zr4
      REAL(wp) :: ze, zbw, zc, zd, zaw, zb1, za1, zkw, zk0
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zsva
      REAL(wp), DIMENSION(jpk)         :: zwkx, zwky, zwkz
      REAL(wp) :: fsatg
      REAL(wp) :: pfps, pfpt, pfphp  

      ! Adiabatic laspse rate fsatg, defined as the change of temperature
      ! per unit pressure for adiabatic change of pressure of an element
      ! of seawater (bryden,h.,1973,deep-sea res.,20,401-408).
      ! units:
      !      pressure        pfphp    decibars
      !      temperature     pfpt     deg celsius (ipts-68)
      !      salinity        pfps     (ipss-78)
      !      adiabatic       fsatg    deg. c/decibar
      ! checkvalue: atg=3.255976e-4 c/dbar for pfps=40 (ipss-78),
      ! pfpt=40 deg c, pfphp=10000 decibars
      
      fsatg(pfps,pfpt,pfphp)   &
         = (((-2.1687e-16*pfpt+1.8676e-14)*pfpt-4.6206e-13)*pfphp    &
         +((2.7759e-12*pfpt-1.1351e-10)*(pfps-35.)+((-5.4481e-14*pfpt    &
         +8.733e-12)*pfpt-6.7795e-10)*pfpt+1.8741e-8))*pfphp    &
         +(-4.2393e-8*pfpt+1.8932e-6)*(pfps-35.)    &
         +((6.6228e-10*pfpt-6.836e-8)*pfpt+8.5258e-6)*pfpt+3.5803e-5
      !!----------------------------------------------------------------------

      ! 1. height dynamic
      ! -----------------
      ! depth for reference

      zgdsup = 1500.
      
      ! below for hdyn levitus
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dia_hdy : computation of dynamical heigh'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
# if defined key_s_coord || defined key_partial_steps
         ! Dynamic height diagnostics  not yet implemented
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          key_s_coord or key_partial_steps used'
         IF(lwp) WRITE(numout,*) '          Dynamical height diagnostics not yet implemented'
         nstop = nstop + 1
# endif

         DO jk = 1, jpk
            IF( fsdepw(1,1,jk) > zgdsup ) GOTO 110
         END DO
         IF(lwp) WRITE(numout,*)'problem zgdsup greater than gdepw(jpk)'
         STOP 'dia_hdy'
110      CONTINUE
         ihdsup = jk - 1
         IF(lwp) WRITE(numout,*)' ihdsup = ', ihdsup

         ! Interpolation coefficients for zgdsup-gdepw(ihdsup) layer

         za = fsdepw(1,1,ihdsup  )
         zb = fsdepw(1,1,ihdsup+1)
         IF( za > zgdsup .OR. zb < zgdsup ) THEN
            IF(lwp) WRITE(numout,*) za, zb, ihdsup, zgdsup
            IF(lwp) WRITE(numout,*) ' bad ihdsup'
            STOP
         ENDIF
         
         zciint = (zgdsup - za) / (zb - za)

         ! Computes the specific volume reference in situ temperature
         
         DO jk = 1, jpk
            zp = 0.e0
            zh = fsdept(1,1,jk)
            zt = 0.e0
            zs = 35.
            zxk= zh * fsatg( zs, zt, zp )
            zt = zt + 0.5 * zxk
            zq = zxk
            zp = zp + 0.5 * zh
            zxk= zh*fsatg( zs, zt, zp )
            zt = zt + 0.29289322 * ( zxk - zq )
            zq = 0.58578644 * zxk + 0.121320344 * zq
            zxk= zh * fsatg( zs, zt, zp )
            zt = zt + 1.707106781 * ( zxk - zq )
            zq = 3.414213562 * zxk - 4.121320344 * zq
            zp = zp + 0.5 * zh
            zxk= zh * fsatg( zs, zt, zp )
            zwkx(jk) = zt + ( zxk - 2.0 * zq ) / 6.0
         END DO

         ! In situ density (add the compression terms)

         DO jk = 1, jpk
            zt = zwkx(jk)
            zs = 35.
            ! square root salinity
            zsr = sqrt( abs( zs ) )
            zwky(jk) = zsr
            ! compute density pure water at atm pressure
            zr1= ((((6.536332e-9*zt-1.120083e-6)*zt+1.001685e-4)*zt   &
               -9.095290e-3)*zt+6.793952e-2)*zt+999.842594
            ! seawater density atm pressure
            zr2= (((5.3875e-9*zt-8.2467e-7)*zt+7.6438e-5)*zt   &
               -4.0899e-3)*zt+8.24493e-1
            zr3= (-1.6546e-6*zt+1.0227e-4)*zt-5.72466e-3
            zr4= 4.8314e-4
            zwkz(jk)= (zr4*zs + zr3*zsr + zr2)*zs + zr1
         END DO

         DO jk = 1, jpk
            zt = zwkx(jk)
            zs = 35.
            zsr= zwky(jk)
            zh = fsdept(1,1,jk)

            ze = ( 9.1697e-11*zt+2.0816e-9 ) *zt-9.9348e-8
            zbw= ( 5.2787e-9*zt-6.12293e-7 ) * zt+8.50935e-6
            zb = zbw + ze * zs

            zd = 1.91075e-4
            zc = (-1.6078e-6*zt-1.0981e-5)*zt+2.2838e-3
            zaw= ((-5.77905e-7*zt+1.16092e-4)*zt+1.43713e-3)*zt+3.239908
            za = ( zd*zsr + zc)*zs + zaw

            zb1= (-5.3009e-3*zt+1.6483e-1)*zt+7.944e-1
            za1= ((-6.1670e-4*zt+1.09987e-1)*zt-6.03459)*zt+546.746
            zkw= (((-5.155288e-4*zt+1.360477e-1)*zt-23.27105)*zt   &
                +1484.206)*zt+196522.1
            zk0= (zb1*zsr + za1)*zs + zkw
            ! evaluate pressure polynomial
            zwkz(jk) = zwkz(jk) / ( 1.0 - zh / ( zk0+zh*(za+zb*zh) ) )
         END DO

         DO jk = 1, jpk
            rhosp(jk) = zwkz(jk)
         END DO
      ENDIF

      ! Computes the specific volume anomaly

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) /= 0. ) THEN
                  zsva(ji,jj,jk) = ( rau0*rhd(ji,jj,jk)+rau0 -rhosp(jk) ) / rhosp(jk)
               ELSE
                  zsva(ji,jj,jk)=0.
               ENDIF
            END DO
         END DO
      END DO

      ! zfacto coefficient to cmg
      
      ! zfacto= 1.  e+2
      !           mg->cmg
      zfacto = 1.0 * 1.e2
      
      ! Fisrt compute at depth ik=ihdsup
      
      ik = ihdsup
      DO jj = 1, jpj
         DO ji = 1, jpi
            zhd = zfacto * zciint * fse3t(ji,jj,ik) * zsva(ji,jj,ik)
            hdy(ji,jj,ik) = zhd * tmask(ji,jj,ik) * tmask(ji,jj,ik-1)
         END DO
      END DO
      
      ! Then compute other terms except level jk=1
      
      DO jk = ihdsup-1, 2, -1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zhd = hdy(ji,jj,jk+1) + zfacto * fse3t(ji,jj,jk) * zsva(ji,jj,jk)
               hdy(ji,jj,jk) = zhd * tmask(ji,jj,jk) * tmask(ji,jj,jk-1)
            END DO
         END DO
      END DO
      
      ! Then compute other the last layer term jk=1
      
      ik = 1
      DO jj = 1, jpj
         DO ji = 1, jpi
            zhd = hdy(ji,jj,ik+1) + zfacto * fse3t(ji,jj,ik) * zsva(ji,jj,ik)
            hdy(ji,jj,ik) = zhd * tmask(ji,jj,ik)
         END DO
      END DO

   END SUBROUTINE dia_hdy

#else
   !!----------------------------------------------------------------------
   !!   Default option :                       NO dynamic heigh diagnostics
   !!----------------------------------------------------------------------
   USE in_out_manager
   LOGICAL, PUBLIC, PARAMETER ::   lk_diahdy = .FALSE.   !: dynamical heigh flag
CONTAINS
   SUBROUTINE dia_hdy( kt )               ! Empty routine
      if(lwp) WRITE(numout,*) 'diahdy: You should not have seen this print! error?', kt
   END SUBROUTINE dia_hdy
#endif

   !!======================================================================
END MODULE diahdy
