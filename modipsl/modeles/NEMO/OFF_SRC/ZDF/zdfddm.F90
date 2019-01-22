MODULE zdfddm
   !!======================================================================
   !!                       ***  MODULE  zdfddm  ***
   !! Ocean physics : double diffusion mixing parameterization
   !!======================================================================
#if defined key_zdfddm   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_zdfddm' :                                     double diffusion
   !!----------------------------------------------------------------------
   !!   zdf_ddm       : compute the Ks for salinity
   !!   zdf_ddm_init  : read namelist and control the parameters
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/ZDF/zdfddm.F90,v 1.2 2005/11/16 16:16:03 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zdf_ddm     ! called by step.F90

   !! * Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfddm = .TRUE.    !: double diffusive mixing flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      avs ,               &  !: salinity vertical diffusivity coeff. at w-point
      rrau                   !: heat/salt buoyancy flux ratio

   !! * Module variables
   REAL(wp) ::            & !!! * double diffusive mixing namelist *
      avts  = 1.e-4_wp      ! maximum value of avs for salt fingering

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"

CONTAINS

   SUBROUTINE zdf_ddm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_ddm  ***
      !!                    
      !! ** Purpose :   Add to the vertical eddy diffusivity coefficient the 
      !!      effect of salt fingering and diffusive convection. 
      !!
      !! ** Method  :   Diapycnal mixing is increased in case of double
      !!      diffusive mixing (i.e. salt fingering and diffusive layering)
      !!      following Merryfield et al. (1999). The rate of double diffusive 
      !!      mixing depend on the buoyancy ratio: Rrau=alpha/beta dk[T]/dk[S]
      !!      which is computed in rn2.F
      !!         * salt fingering (Schmitt 1981):
      !!      for Rrau > 1 and rn2 > 0 : zavfs = avts / ( 1 + (Rrau/hsbfr)^6 )
      !!      for Rrau > 1 and rn2 > 0 : zavfs = O
      !!      otherwise                : zavft = 0.7 zavs / Rrau
      !!         * diffusive layering (Federov 1988):
      !!      for 0< Rrau < 1 and rn2 > 0 : zavdt = 1.3635e-6  
      !!                                 * exp( 4.6 exp(-0.54 (1/Rrau-1) ) )
      !!      otherwise                   : zavdt = 0 
      !!      for .5 < Rrau < 1 and rn2 > 0 : zavds = zavdt (1.885 Rrau -0.85)
      !!      for  0 < Rrau <.5 and rn2 > 0 : zavds = zavdt 0.15 Rrau      
      !!      otherwise                     : zavds = 0 
      !!         * update the eddy diffusivity:
      !!      avt = avt + zavft + zavdt
      !!      avs = avs + zavfs + zavds
      !!      avmu, avmv are required to remain at least above avt and avs.
      !!      
      !! ** Action  :   avt, avs : update vertical eddy diffusivity coef.
      !!                           for temperature and salinity
      !!
      !! References :
      !!      Merryfield et al., JPO, 29, 1124-1142, 1999.
      !! History :
      !!        !  00-08  (G. Madec)  double diffusive mixing
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step indexocean time step

      !!----------------------------------------------------------------------


      IF ( kt == nit000 )   CALL zdf_ddm_init          ! Initialization (first time-step only)


   END SUBROUTINE zdf_ddm
   
   
   SUBROUTINE zdf_ddm_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_ddm_init  ***
      !!
      !! ** Purpose :   Initialization of double diffusion mixing scheme
      !!
      !! ** Method  :   Read the nammbf namelist and check the parameter values
      !!      called by zdf_ddm at the first timestep (nit000)
      !!
      !! History :
      !!   8.5  !  02-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      NAMELIST/namddm/ avts
      !!----------------------------------------------------------------------

      ! Read Namelist namddm : double diffusion mixing scheme
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, namddm )


      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_ddm : double diffusive mixing'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '          Namelist namddm : set dd mixing parameter'
         WRITE(numout,*) '             maximum avs for dd mixing      avts   = ', avts
         WRITE(numout,*)
      ENDIF

   END SUBROUTINE zdf_ddm_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :          Dummy module          No double diffusion
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfddm = .FALSE.   !: double diffusion flag
CONTAINS
   SUBROUTINE zdf_ddm( kt )           ! Dummy routine
      WRITE(*,*) 'zdf_ddm: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_ddm
#endif

   !!======================================================================
END MODULE zdfddm
