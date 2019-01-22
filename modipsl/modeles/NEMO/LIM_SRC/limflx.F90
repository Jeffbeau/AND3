MODULE limflx
   !!======================================================================
   !!                       ***  MODULE limflx   ***
   !!           computation of the flux at the sea ice/ocean interface
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim'                                     LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_flx  : flux at the ice / ocean interface
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce          ! ocean parameters
   USE phycst           ! physical constants
   USE ocfzpt           ! surface ocean freezing point
   USE ice_oce          ! sea-ice variable
   USE flx_oce          ! sea-ice/ocean forcings variables
   USE ice              ! LIM sea-ice variables
   USE flxblk           ! bulk formulea
   USE lbclnk           ! ocean lateral boundary condition
   USE in_out_manager   ! I/O manager
   USE albedo           ! albedo parameters
   USE prtctl           ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC lim_flx       ! called by lim_step

   !! * Module variables
   REAL(wp)  ::           &  ! constant values
      epsi16 = 1.e-16  ,  &
      rzero  = 0.e0    ,  &
      rone   = 1.e0

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limflx.F90,v 1.8 2005/12/21 10:46:27 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_flx
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE lim_flx ***
      !!  
      !! ** Purpose : Computes the mass and heat fluxes to the ocean
      !!         
      !! ** Action  : - Initialisation of some variables
      !!              - comput. of the fluxes at the sea ice/ocean interface
      !!     
      !! ** Outputs : - fsolar  : solar heat flux at sea ice/ocean interface
      !!              - fnsolar : non solar heat flux 
      !!              - fsalt   : salt flux at sea ice/ocean interface
      !!              - fmass   : freshwater flux at sea ice/ocean interface
      !!
      !! ** References :
      !!       H. Goosse et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90
      !!         original    : 00-01 (LIM)
      !!         addition    : 02-07 (C. Ethe, G. Madec)
      !!---------------------------------------------------------------------
      !! * Local variables
      INTEGER ::   ji, jj         ! dummy loop indices

      INTEGER ::   &
         ifvt, i1mfr, idfr ,   &  ! some switches
         iflt, ial, iadv, ifral, ifrdv
      
      REAL(wp) ::   &
         zinda  ,              &  ! switch for testing the values of ice concentration
         z1mthcm                  ! 1 - thcm
!!         zfcm1  ,              &  ! solar  heat fluxes
!!         zfcm2  ,              &  !  non solar heat fluxes
#if defined key_lim_fdd   
      REAL(wp) ::   &
         zfons,                &  ! salt exchanges at the ice/ocean interface
         zpme                     ! freshwater exchanges at the ice/ocean interface
#else
      REAL(wp) ::   &
         zprs  , zfons,        &  ! salt exchanges at the ice/ocean interface
         zpmess                   ! freshwater exchanges at the ice/ocean interface
#endif
      REAL(wp), DIMENSION(jpi,jpj) ::  &
         zfcm1  ,              &  ! solar  heat fluxes
         zfcm2                    !  non solar heat fluxes      
#if defined key_coupled    
      REAL(wp), DIMENSION(jpi,jpj) ::  &
         zalb  ,               &  ! albedo of ice under overcast sky
         zalcn ,               &  ! albedo of ocean under overcast sky
         zalbp ,               &  ! albedo of ice under clear sky
         zaldum                   ! albedo of ocean under clear sky
#endif
      !!---------------------------------------------------------------------
     
      !---------------------------------!
      !      Sea ice/ocean interface    !
      !---------------------------------!
       
       
      !      heat flux at the ocean surface
      !-------------------------------------------------------
       
      DO jj = 1, jpj
         DO ji = 1, jpi
            zinda   = 1.0 - MAX( rzero , SIGN( rone , - ( 1.0 - pfrld(ji,jj) ) ) )
            ifvt    = zinda  *  MAX( rzero , SIGN( rone, -phicif(ji,jj) ) )
            i1mfr   = 1.0 - MAX( rzero , SIGN( rone ,  - ( 1.0 - frld(ji,jj) ) ) )
            idfr    = 1.0 - MAX( rzero , SIGN( rone , frld(ji,jj) - pfrld(ji,jj) ) )
            iflt    = zinda  * (1 - i1mfr) * (1 - ifvt )
            ial     = ifvt   * i1mfr + ( 1 - ifvt ) * idfr
            iadv    = ( 1  - i1mfr ) * zinda
            ifral   = ( 1  - i1mfr * ( 1 - ial ) )   
            ifrdv   = ( 1  - ifral * ( 1 - ial ) ) * iadv 
            z1mthcm =   1. - thcm(ji,jj)       
            !   computation the solar flux at ocean surface
            zfcm1(ji,jj)   = pfrld(ji,jj) * qsr_oce(ji,jj)  + ( 1. - pfrld(ji,jj) ) * fstric(ji,jj)
            !  computation the non solar heat flux at ocean surface
            zfcm2(ji,jj) =  - z1mthcm * zfcm1(ji,jj)   &
               &           + iflt    * ( fscmbq(ji,jj) + ffltbif(ji,jj) )                            &
               &           + ifral   * ( ial * qcmif(ji,jj) + (1 - ial) * qldif(ji,jj) ) / rdt_ice   &
               &           + ifrdv   * ( qfvbq(ji,jj) + qdtcn(ji,jj) ) / rdt_ice

            fsbbq(ji,jj) = ( 1.0 - ( ifvt + iflt ) ) * fscmbq(ji,jj)     ! ???
            
            fsolar (ji,jj) = zfcm1(ji,jj)                                       ! solar heat flux 

            fnsolar(ji,jj) = zfcm2(ji,jj) - fdtcn(ji,jj)                        ! non solar heat flux
         END DO
      END DO
  
       
      !      mass flux at the ocean surface
      !-------------------------------------------------------
       
      DO jj = 1, jpj
         DO ji = 1, jpi
#if defined key_lim_fdd
            !  case of realistic freshwater flux (Tartinville et al., 2001)
            
            !  computing freshwater exchanges at the ice/ocean interface
            zpme = - evap(ji,jj) * frld(ji,jj)            &   !  evaporation over oceanic fraction
               &   + tprecip(ji,jj)                            &   !  total precipitation
               &   - sprecip(ji,jj) * ( 1. - pfrld(ji,jj) )  &   !  remov. snow precip over ice
               &   - rdmsnif(ji,jj) / rdt_ice                   !  freshwaterflux due to snow melting 
            
            !  computing salt exchanges at the ice/ocean interface
            zfons =  ( soce - sice ) * ( rdmicif(ji,jj) / rdt_ice ) 
            
            !  converting the salt flux from ice to a freshwater flux from ocean
            fsalt(ji,jj) = zfons / ( sss_io(ji,jj) + epsi16 )
            
            !  freshwater masses
            fmass(ji,jj) = - zpme 
#else
            !  case of freshwater flux equivalent as salt flux
            !  dilution effect due to evaporation and precipitation
            zprs  = ( tprecip(ji,jj) - sprecip(ji,jj) * ( 1. - pfrld(ji,jj) ) ) * soce  
            !  freshwater flux
            zfons = rdmicif(ji,jj) * ( soce - sice )  &  !  fwf : ice formation and melting
               &   -  dmgwi(ji,jj) * sice             &  !  fwf : salt flx needed to bring the fresh snow to sea/ice salinity
               &   + rdmsnif(ji,jj) * soce               !  fwf to ocean due to snow melting 
            !  salt exchanges at the ice/ocean interface
            zpmess         =  zprs - zfons / rdt_ice - evap(ji,jj) * soce * frld(ji,jj)
            fsalt(ji,jj) =  - zpmess
#endif
         END DO
      END DO


      !-------------------------------------------------------------------!
      !  computation of others transmitting variables from ice to ocean   !
      !------------------------------------------ ------------------------!

      !-----------------------------------------------!
      !   Storing the transmitted variables           !
      !-----------------------------------------------!

      ftaux (:,:) = - tio_u(:,:) * rau0   ! taux ( ice: N/m2/rau0, ocean: N/m2 )
      ftauy (:,:) = - tio_v(:,:) * rau0   ! tauy ( ice: N/m2/rau0, ocean: N/m2 )                
      freeze(:,:) = 1.0 - frld(:,:)       ! Sea ice cover            
      tn_ice(:,:) = sist(:,:)             ! Ice surface temperature                      

#if defined key_coupled            
      zalb  (:,:) = 0.e0
      zalcn (:,:) = 0.e0
      zalbp (:,:) = 0.e0
      zaldum(:,:) = 0.e0

      !------------------------------------------------!
      !  2) Computation of snow/ice and ocean albedo   !
      !------------------------------------------------!
      CALL flx_blk_albedo( zalb, zalcn, zalbp, zaldum )

      alb_ice(:,:) =  0.5 * zalbp(:,:) + 0.5 * zalb (:,:)   ! Ice albedo                       
#endif

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=fsolar, clinfo1=' lim_flx: fsolar : ', tab2d_2=fnsolar, clinfo2=' fnsolar : ')
         CALL prt_ctl(tab2d_1=fmass , clinfo1=' lim_flx: fmass  : ', tab2d_2=fsalt  , clinfo2=' fsalt   : ')
         CALL prt_ctl(tab2d_1=ftaux , clinfo1=' lim_flx: ftaux  : ', tab2d_2=ftauy  , clinfo2=' ftauy   : ')
         CALL prt_ctl(tab2d_1=freeze, clinfo1=' lim_flx: freeze : ', tab2d_2=tn_ice , clinfo2=' tn_ice  : ')
      ENDIF 
   
    END SUBROUTINE lim_flx

#else
   !!----------------------------------------------------------------------
   !!   Default option :        Dummy module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_flx         ! Dummy routine
   END SUBROUTINE lim_flx
#endif 

   !!======================================================================
END MODULE limflx
