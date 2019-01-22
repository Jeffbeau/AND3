MODULE trcfreons
   !!==============================================================
   !!                  ***  MODULE trcfreons  ***
   !!  Passive tracer : CFC main model
   !!==============================================================
#if defined key_cfc
   !!--------------------------------------------------------------
   !!   'key_cfc'                                         CFC model
   !!--------------------------------------------------------------
   !! * Modules used   
   USE daymod
   USE sms
   USE oce_trc
   USE trc


   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_freons        

   !! * Module variables
   REAL(wp), DIMENSION(jptra) :: & ! coefficient for solubility of CFC11 in  mol/l/atm
      soa1, soa2, soa3, soa4, &
      sob1, sob2, sob3

   REAL(wp), DIMENSION(jptra) :: & ! coefficients for schmidt number in degre Celcius
      sca1, sca2, sca3, sca4

   REAL(wp) ::              & ! coefficients for conversion
      xconv1 = 1.0       ,  & ! conversion from to 
      xconv2 = 0.01/3600.,  & ! conversion from cm/h to m/s: 
      xconv3 = 1.0e+3    ,  & ! conversion from mol/l/atm to mol/m3/atm
      xconv4 = 1.0e-12        ! conversion from mol/m3/atm to mol/m3/pptv 

   !! * Substitutions
#  include "passivetrc_substitute.h90"

   !!----------------------------------------------------------------------
   !!  TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/trcfreons.F90,v 1.2 2005/11/14 16:42:43 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_freons( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_freons  ***
      !!
      !! ** Purpose :   Compute the surface boundary contition on freon 11 
      !!      passive tracer associated with air-mer fluxes and add it to 
      !!      the general trend of tracers equations.
      !!
      !! ** Method :
      !!          - get the atmospheric partial pressure - given in pico -
      !!          - computation of solubility ( in 1.e-12 mol/l then in 1.e-9 mol/m3)
      !!          - computation of transfert speed ( given in cm/hour ----> cm/s )
      !!          - the input function is given by : 
      !!            speed * ( concentration at equilibrium - concemtration at surface )
      !!          - the input function is in pico-mol/m3/s and the
      !!            freons concentration in pico-mol/m3
      !!
      !! History :
      !!   8.1  !  99-10  (JC. Dutay)  original code
      !!   9.0  !  04-03  (C. Ethe)  free form + modularity
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step index

      !! * Local declarations
      INTEGER ::  &
         ji, jj, jn, jm

      INTEGER ::   &
         iyear_beg, iyear_end, &
         imonth, im1, im2

      REAL(wp) :: &
         ztap, zdtap, &
         zt1, zt2, zt3, zv2

      REAL(wp), DIMENSION(jphem,jptra)   ::  &   
         zpatm       ! atmospheric function

      REAL(wp) ::  & 
         zsol,     & ! solubility
         zsch        ! schmidt number 

      
      REAL(wp), DIMENSION(jpi,jpj,jptra)   ::  & 
         zca_cfc,  & ! concentration
         zak_cfc     ! transfert coefficients

      !!----------------------------------------------------------------------


      IF( kt == nittrc000 )   CALL trc_freons_cst

      ! Temporal interpolation
      ! ----------------------
      iyear_beg = nyear + ( nyear_res - 1900 - nyear_beg  )
      imonth    = nmonth

      IF ( imonth .LE. 6 ) THEN
         iyear_beg = iyear_beg - 2 + nyear_beg
         im1       = 6 - imonth + 1
         im2       = 6 + imonth - 1
      ELSE
         iyear_beg = iyear_beg - 1 + nyear_beg
         im1       = 12 - imonth + 7
         im2       =      imonth - 7
      ENDIF

      iyear_end = iyear_beg + 1




      !  Temporal and spatial interpolation at time k
      ! --------------------------------------------------
      DO jn = 1, jptra
         DO  jm = 1, jphem
            zpatm(jm,jn) = (  p_cfc(iyear_beg, jm, jn) * FLOAT (im1)  &
               &           +  p_cfc(iyear_end, jm, jn) * FLOAT (im2) ) / 12.
         ENDDO
      END DO

      DO jn = 1, jptra
         DO jj = 1, jpj 
            DO ji = 1, jpi
               pp_cfc(ji,jj,jn) =     xphem(ji,jj)   * zpatm(1,jn)  &
                  &           + ( 1.- xphem(ji,jj) ) * zpatm(2,jn)
            END DO
         END DO
      ENDDO


      !------------------------------------------------------------
      ! Computation of concentration at equilibrium : in picomol/l
      ! -----------------------------------------------------------

      DO jn = 1, jptra
         DO jj = 1 , jpj
            DO ji = 1 , jpi
               ! coefficient for solubility for CFC-11/12 in  mol/l/atm
               IF( tmask(ji,jj,1) .GE. 0.5 ) THEN
                  ztap  = ( tn(ji,jj,1) + 273.16 )* 0.01
                  zdtap = ( sob3(jn) * ztap + sob2(jn))* ztap + sob1(jn) 
                  zsol  =  EXP ( soa1(jn) + soa2(jn) / ztap + soa3(jn) * LOG ( ztap )   &
                     &                   + soa4(jn) * ztap * ztap + sn(ji,jj,1) * zdtap ) 
               ELSE
                  zsol  = 0.
               ENDIF
               ! conversion from mol/l/atm to mol/m3/atm and from mol/m3/atm to mol/m3/pptv    
               zsol = xconv4 * xconv3 * zsol * tmask(ji,jj,1)  
               ! concentration at equilibrium
               zca_cfc(ji,jj,jn) = xconv1 * pp_cfc(ji,jj,jn) * zsol * tmask(ji,jj,1)             
            END DO
         END DO
      ENDDO


      !-------------------------------
      ! Computation of speed transfert
      ! ------------------------------

      DO jn = 1, jptra
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Schmidt number
               zt1  = tn(ji,jj,1)
               zt2  = zt1 * zt1 
               zt3  = zt1 * zt2
               zsch = sca1(jn) + sca2(jn) * zt1 + sca3(jn) * zt2 + sca4(jn) * zt3
               ! speed transfert : formulae of wanninkhof 1992
               zv2 = vatm(ji,jj) * vatm(ji,jj)
               zsch = zsch / 660.
               zak_cfc(ji,jj,jn) = ( 0.39 * xconv2 * zv2 / SQRT(zsch) ) * tmask(ji,jj,1)
            END DO
         END DO
      ENDDO

      !----------------------------------------------------------------
      ! Input function  : speed *( conc. at equil - concen at surface )
      ! trn in pico-mol/l idem qtr; ak in en m/s
      !-----------------------------------------------------------------

      DO jn = 1, jptra
         DO jj = 1, jpj
            DO ji = 1, jpi
               qtr(ji,jj,jn) = -zak_cfc(ji,jj,jn) * ( trn(ji,jj,1,jn) - zca_cfc(ji,jj,jn) )   &
                  &                               * tmask(ji,jj,1) * ( 1. - freeze(ji,jj) )
            END DO
         END DO
      ENDDO

      ! ---------------------
      ! Add the trend
      ! ---------------------

      DO jn = 1, jptra
         DO  jj = 1, jpj
            DO  ji = 1, jpi
               tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + qtr(ji,jj,jn) / fse3t(ji,jj,1) 
            END DO
         END DO
      ENDDO

      ! --------------------------------------------
      ! cumulation of tracer flux at each time step
      ! --------------------------------------------
      DO jn = 1, jptra
         DO jj = 1, jpj
            DO ji = 1, jpi
               qint(ji,jj,jn) = qint (ji,jj,jn) + qtr(ji,jj,jn) * rdt
            END DO
         END DO
      ENDDO


   END SUBROUTINE trc_freons

   SUBROUTINE trc_freons_cst
      !!---------------------------------------------------------------------
      !!                     ***  trc_freons_cst  ***  
      !!
      !!   Purpose : sets constants for CFC model
      !!  ---------
      !!
      !!
      !! History :
      !!   8.2  !  04-06  (JC. Dutay)  original code
      !!   9.0  !  05-10  (C. Ethe) Modularity 
      !!---------------------------------------------------------------------
      !!  TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/trcfreons.F90,v 1.2 2005/11/14 16:42:43 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
      !!----------------------------------------------------------------
      !! Local declarations
      INTEGER :: jn

      DO jn = 1, jptra
         IF ( jn == jp11 ) THEN
            ! coefficient for solubility of CFC11 in  mol/l/atm
            soa1(jn) = -229.9261 
            soa2(jn) =  319.6552
            soa3(jn) =  119.4471
            soa4(jn) =  -1.39165
            sob1(jn) =  -0.142382
            sob2(jn) =   0.091459
            sob3(jn) =  -0.0157274
            
            ! coefficients for schmidt number in degre Celcius
            sca1(jn) = 3501.8
            sca2(jn) = -210.31
            sca3(jn) = 6.1851
            sca4(jn) = -0.07513

         ELSE IF( jn == jp12 ) THEN

            ! coefficient for solubility of CFC12 in  mol/l/atm
            soa1(jn) = -218.0971
            soa2(jn) =  298.9702
            soa3(jn) =  113.8049
            soa4(jn) =  -1.39165
            sob1(jn) =  -0.143566
            sob2(jn) =   0.091015
            sob3(jn) =  -0.0153924
                        
            ! coefficients for schmidt number in degre Celcius
            sca1(jn) =  3845.4 
            sca2(jn) = -228.95
            sca3(jn) = 6.1908 
            sca4(jn) = -0.067430
         ENDIF
      ENDDO

   END SUBROUTINE trc_freons_cst
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Dummy module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_freons( kt )       ! Empty routine
      WRITE(*,*) 'trc_freons: You should not have seen this print! error?', kt
   END SUBROUTINE trc_freons
#endif

   !!======================================================================
END MODULE trcfreons
