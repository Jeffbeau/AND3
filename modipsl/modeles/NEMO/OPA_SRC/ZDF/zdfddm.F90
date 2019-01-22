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
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control

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
      avts  = 1.e-4_wp ,  &  ! maximum value of avs for salt fingering
      hsbfr = 1.6_wp         ! heat/salt buoyancy flux ratio

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfddm.F90,v 1.6 2005/09/02 15:45:43 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

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

      !! * Local declarations
      INTEGER ::   ji, jj , jk              ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zmsks, zmskf,                    & ! temporary workspace 
         zmskd1, zmskd2, zmskd3             !    "           "
      REAL(wp) ::   &
         zinr, zrr,                       & ! temporary scalars
         zavft, zavfs,                    & !    "         "
         zavdt, zavds                       !    "         "
      !!----------------------------------------------------------------------


      IF ( kt == nit000 )   CALL zdf_ddm_init          ! Initialization (first time-step only)


      ! Compute avs
      ! -----------
      !                                                ! ===============
      DO jk = 2, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Define the mask 
         ! ---------------
         ! only retains positive value of rrau
         rrau(:,:,jk) = MAX( 1.e-20, rrau(:,:,jk) )

         ! indicators:
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! stability indicator: msks=1 if rn2>0; 0 elsewhere
               IF( rn2(ji,jj,jk) + 1.e-12  <= 0. ) THEN
                  zmsks(ji,jj) = 0.e0
               ELSE
                  zmsks(ji,jj) = 1.e0
               ENDIF
               ! salt fingering indicator: msksf=1 if rrau>1; 0 elsewhere            
               IF( rrau(ji,jj,jk) <= 1. ) THEN
                  zmskf(ji,jj) = 0.e0
               ELSE
                  zmskf(ji,jj) = 1.e0
               ENDIF
               ! diffusive layering indicators: 
               !   mskdl1=1 if 0<rrau<1; 0 elsewhere
               IF( rrau(ji,jj,jk) >= 1. ) THEN
                  zmskd1(ji,jj) = 0.e0
               ELSE
                  zmskd1(ji,jj) = 1.e0
               ENDIF
               !   mskdl2=1 if 0<rrau<0.5; 0 elsewhere
               IF( rrau(ji,jj,jk) >= 0.5 ) THEN
                  zmskd2(ji,jj) = 0.e0
               ELSE
                  zmskd2(ji,jj) = 1.e0
               ENDIF
               !   mskdl3=1 if 0.5<rrau<1; 0 elsewhere
               IF( rrau(ji,jj,jk) <= 0.5 .OR. rrau(ji,jj,jk) >= 1. ) THEN
                  zmskd3(ji,jj) = 0.e0
               ELSE
                  zmskd3(ji,jj) = 1.e0
               ENDIF
            END DO
         END DO
         ! mask zmsk in order to have avt and avs masked
         zmsks(:,:) = zmsks(:,:) * tmask(:,:,jk)


         ! Update avt and avs
         ! ------------------
         ! Constant eddy coefficient: reset to the background value
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zinr = 1./rrau(ji,jj,jk)
               ! salt fingering
               zrr = rrau(ji,jj,jk)/hsbfr
               zrr = zrr * zrr
               zavfs = avts / ( 1 + zrr*zrr*zrr ) * zmsks(ji,jj) *zmskf(ji,jj)
               zavft = 0.7 * zavfs / rrau(ji,jj,jk)
               ! diffusive layering
               zavdt = 1.3635e-6 * EXP(4.6*EXP(-0.54*(zinr-1.) ) )   &
                                 * zmsks(ji,jj) * zmskd1(ji,jj)
               zavds = zavdt * zmsks(ji,jj)   &
                     * ( (1.85 * rrau(ji,jj,jk) - 0.85 ) * zmskd3(ji,jj)   &
                         + zavdt * 0.15 * rrau(ji,jj,jk) * zmskd2(ji,jj)  )
               ! add to the eddy viscosity coef. previously computed
               avs (ji,jj,jk) = avt(ji,jj,jk) + zavfs + zavds
               avt (ji,jj,jk) = avt(ji,jj,jk) + zavft + zavdt
            END DO
         END DO


         ! Increase avmu, avmv if necessary
         ! --------------------------------
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               avmu(ji,jj,jk) = MAX( avmu(ji,jj,jk),    &
                                     avt(ji,jj,jk), avt(ji+1,jj,jk),   &
                                     avs(ji,jj,jk), avs(ji+1,jj,jk) )   &
                              * umask(ji,jj,jk)
               avmv(ji,jj,jk) = MAX( avmv(ji,jj,jk),    &
                                     avt(ji,jj,jk), avt(ji,jj+1,jk),   &
                                     avs(ji,jj,jk), avs(ji,jj+1,jk) )   &
                              * vmask(ji,jj,jk)
            END DO
         END DO
         !                                                ! ===============
      END DO                                              !   End of slab
      !                                                   ! ===============
      
      ! Lateral boundary conditions on ( avt, avs, avmu, avmv )   (unchanged sign)
      ! -------------------------------========================
      CALL lbc_lnk( avt , 'W', 1. )
      CALL lbc_lnk( avs , 'W', 1. )
      CALL lbc_lnk( avmu, 'U', 1. ) 
      CALL lbc_lnk( avmv, 'V', 1. )

      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=avt , clinfo1=' ddm  - t: ', tab3d_2=avs , clinfo2=' s: ', ovlap=1, kdim=jpk)
         CALL prt_ctl(tab3d_1=avmu, clinfo1=' ddm  - u: ', tab3d_2=avmv, clinfo2=' v: ', ovlap=1, kdim=jpk)
      ENDIF
      
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
      NAMELIST/namddm/ avts, hsbfr
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
         WRITE(numout,*) '             heat/salt buoyancy flux ratio  hsbfr  = ', hsbfr
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
!      WRITE(*,*) 'zdf_ddm: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_ddm
#endif

   !!======================================================================
END MODULE zdfddm
