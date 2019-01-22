MODULE trcbbl
   !!==============================================================================
   !!                       ***  MODULE  trcbbl  ***
   !! Ocean passive tracers physics :  advective and/or diffusive bottom boundary 
   !!                                  layer scheme
   !!==============================================================================
#if  defined key_passivetrc && ( defined key_trcbbl_dif   ||   defined key_trcbbl_adv ) && ! defined key_cfg_1d
   !!----------------------------------------------------------------------
   !!   'key_trcbbl_dif'   or            diffusive bottom boundary layer
   !!   'key_trcbbl_adv'                 advective bottom boundary layer
   !!----------------------------------------------------------------------
   !!   trc_bbl_dif  : update the passive tracer trends due to the bottom
   !!                  boundary layer (diffusive only)
   !!   trc_bbl_adv  : update the passive tracer trends due to the bottom
   !!                  boundary layer (advective and/or diffusive)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables
   USE prtctl_trc          ! Print control for debbuging
   USE eosbn2
   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_bbl_dif    ! routine called by step.F90
   PUBLIC trc_bbl_adv    ! routine called by step.F90

   !! * Shared module variables
# if defined key_trcbbl_dif
   LOGICAL, PUBLIC, PARAMETER ::    &  !:
      lk_trcbbl_dif = .TRUE.   !: advective bottom boundary layer flag

# else
   LOGICAL, PUBLIC, PARAMETER ::    &  !:
      lk_trcbbl_dif = .FALSE.  !: advective bottom boundary layer flag
# endif

# if defined key_trcbbl_adv
   LOGICAL, PUBLIC, PARAMETER ::    &  !:
      lk_trcbbl_adv = .TRUE.   !: advective bottom boundary layer flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
       u_trc_bbl, v_trc_bbl, &  !: velocity involved in exhanges in the advective BBL
       w_trc_bbl                !: vertical increment of velocity due to advective BBL
       !                        !  only affect tracer vertical advection
# else
   LOGICAL, PUBLIC, PARAMETER ::    &  !:
      lk_trcbbl_adv = .FALSE.  !: advective bottom boundary layer flag
# endif

   !! * Module variables
   INTEGER, DIMENSION(jpi,jpj) ::   &  !:
      mbkt, mbku, mbkv                 ! ???

   REAL(wp) ::        &  !!! * trcbbl namelist *
      atrcbbl = 1.e+3      ! lateral coeff. for bottom boundary layer scheme (m2/s)

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcbbl.F90,v 1.11 2006/04/11 13:49:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_bbl_dif( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_bbl_dif  ***
      !!
      !! ** Purpose :   Compute the before tracer trend associated 
      !!      with the bottom boundary layer and add it to the general trend 
      !!      of tracer equations. The bottom boundary layer is supposed to be
      !!      a purely diffusive bottom boundary layer.
      !!
      !! ** Method  :   When the product grad( rho) * grad(h) < 0 (where grad 
      !!      is an along bottom slope gradient) an additional lateral diffu-
      !!      sive trend along the bottom slope is added to the general tracer
      !!      trend, otherwise nothing is done.
      !!      Second order operator (laplacian type) with variable coefficient
      !!      computed as follow for temperature (idem on s): 
      !!         difft = 1/(e1t*e2t*e3t) { di-1[ ahbt e2u*e3u/e1u di[ztb] ]
      !!                                 + dj-1[ ahbt e1v*e3v/e2v dj[ztb] ] }
      !!      where ztb is a 2D array: the bottom ocean temperature and ahtb
      !!      is a time and space varying diffusive coefficient defined by:
      !!         ahbt = zahbp    if grad(rho).grad(h) < 0
      !!              = 0.       otherwise.
      !!      Note that grad(.) is the along bottom slope gradient. grad(rho)
      !!      is evaluated using the local density (i.e. referenced at the
      !!      local depth). Typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!      Add this before trend to the general trend tra of the 
      !!      botton ocean tracer point:
      !!         tra = tra + difft
      !!
      !! ** Action  : - update tra at the bottom level with the bottom
      !!                boundary layer trend
      !!
      !! References :
      !!     Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!
      !! History :
      !!   8.0  !  96-06  (L. Mortier)  Original code
      !!   8.0  !  97-11  (G. Madec)  Optimization
      !!   8.5  !  02-08  (G. Madec)  free form + modules
      !!   9.0  !  04-03  (C. Ethe)   Adaptation for passive tracers
      !!----------------------------------------------------------------------
      !! * Arguments 
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      INTEGER ::   ji, jj,jn                ! dummy loop indices
      INTEGER ::   ik
      INTEGER ::   ii0, ii1, ij0, ij1       ! temporary integers
#  if defined key_partial_steps
      INTEGER  ::   iku1, iku2, ikv1,ikv2   ! temporary intergers
      REAL(wp) ::   ze3u, ze3v              ! temporary scalars
#  else
      INTEGER ::   iku, ikv
#  endif
      REAL(wp) ::   &
         zsign, zt, zs, zh, zalbet,      &  ! temporary scalars
         zgdrho, zbtr, ztra
      REAL(wp), DIMENSION(jpi,jpj) ::    &
        zki, zkj, zkx, zky,    &  ! temporary workspace arrays
        ztnb, zsnb, zdep,                &
        ztrb, zahu, zahv
      CHARACTER (len=22) :: charout
      REAL(wp) ::   &
         fsalbt, pft, pfs, pfh              ! statement function
      !!----------------------------------------------------------------------
      ! ratio alpha/beta
      ! ================
      !  fsalbt: ratio of thermal over saline expension coefficients
      !       pft :  potential temperature in degrees celcius
      !       pfs :  salinity anomaly (s-35) in psu
      !       pfh :  depth in meters

      fsalbt( pft, pfs, pfh ) =                                              &
         ( ( ( -0.255019e-07 * pft + 0.298357e-05 ) * pft                    &
                                   - 0.203814e-03 ) * pft                    &
                                   + 0.170907e-01 ) * pft                    &
                                   + 0.665157e-01                            &
         +(-0.678662e-05 * pfs - 0.846960e-04 * pft + 0.378110e-02 ) * pfs   &
         +  ( ( - 0.302285e-13 * pfh                                         &
                - 0.251520e-11 * pfs                                         &
                + 0.512857e-12 * pft * pft          ) * pfh                  &
                                     - 0.164759e-06   * pfs                  &
             +(   0.791325e-08 * pft - 0.933746e-06 ) * pft                  &
                                     + 0.380374e-04 ) * pfh   
      !!----------------------------------------------------------------------


      IF( kt == nittrc000 )   CALL trc_bbl_init


      ! 0. 2D fields of bottom temperature and salinity, and bottom slope
      ! -----------------------------------------------------------------
      ! mbathy= number of w-level, minimum value=1 (cf dommsk.F)

#  if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij   ! vector opt. (forced unrolling)
#  else
      DO jj = 1, jpj
         DO ji = 1, jpi
#  endif
            ik = mbkt(ji,jj)                              ! index of the bottom ocean T-level
            ztnb(ji,jj) = tn(ji,jj,ik) * tmask(ji,jj,1)   ! masked now T and S at ocean bottom 
            zsnb(ji,jj) = sn(ji,jj,ik) * tmask(ji,jj,1)
            zdep(ji,jj) = fsdept(ji,jj,ik)                ! depth of the ocean bottom T-level
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#  endif
      END DO

#  if defined key_partial_steps
      ! partial steps correction 
#   if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#   else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
#   endif
            iku1 = MAX( mbathy(ji+1,jj  )-1, 1 )
            iku2 = MAX( mbathy(ji  ,jj  )-1, 1 )
            ikv1 = MAX( mbathy(ji  ,jj+1)-1, 1 )
            ikv2 = MAX( mbathy(ji  ,jj  )-1, 1 )
            ze3u = MIN( fse3u(ji,jj,iku1), fse3u(ji,jj,iku2) ) 
            ze3v = MIN( fse3v(ji,jj,ikv1), fse3v(ji,jj,ikv2) ) 
            zahu(ji,jj) = atrcbbl * e2u(ji,jj) * ze3u / e1u(ji,jj) * umask(ji,jj,1)
            zahv(ji,jj) = atrcbbl * e1v(ji,jj) * ze3v / e2v(ji,jj) * vmask(ji,jj,1)
#   if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#   endif
      END DO
#  else
#   if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#   else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
#   endif
            iku = mbku(ji,jj)
            ikv = mbkv(ji,jj)
            zahu(ji,jj) = atrcbbl * e2u(ji,jj) * fse3u(ji,jj,iku) / e1u(ji,jj) * umask(ji,jj,1)
            zahv(ji,jj) = atrcbbl * e1v(ji,jj) * fse3v(ji,jj,ikv) / e2v(ji,jj) * vmask(ji,jj,1)
#   if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#   endif
      END DO
#  endif

!!
!!     OFFLINE VERSION OF DIFFUSIVE BBL
!!
#if defined key_off_tra

      ! 2. Additional second order diffusive trends
      ! -------------------------------------------

      DO jn = 1, jptra
         ! first derivative (gradient)
         
#  if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = 1, jpij   ! vector opt. (forced unrolling)
#  else
         DO jj = 1, jpj
            DO ji = 1, jpi
#  endif
               ik = mbkt(ji,jj) 
               ztrb(ji,jj) = trb(ji,jj,ik,jn) * tmask(ji,jj,1)
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
#  endif
         END DO

#  if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#  else
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
#  endif
               zkx(ji,jj) = bblx(ji,jj) * zahu(ji,jj) * ( ztrb(ji+1,jj) - ztrb(ji,jj) )
               zky(ji,jj) = bbly(ji,jj) * zahv(ji,jj) * ( ztrb(ji,jj+1) - ztrb(ji,jj) )
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
#  endif
         END DO
!!
!!  ONLINE VERSION OF DIFFUSIVE BBL
!!
#else
      ! 1. Criteria of additional bottom diffusivity: grad(rho).grad(h)<0
      ! --------------------------------------------
      ! Sign of the local density gradient along the i- and j-slopes
      ! multiplied by the slope of the ocean bottom
	SELECT CASE ( neos )

      	CASE ( 0 )               ! Jackett and McDougall (1994) formulation

#  if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#  else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
#  endif
            ! temperature, salinity anomalie and depth
            zt = 0.5 * ( ztnb(ji,jj) + ztnb(ji+1,jj) )
            zs = 0.5 * ( zsnb(ji,jj) + zsnb(ji+1,jj) ) - 35.0
            zh = 0.5 * ( zdep(ji,jj) + zdep(ji+1,jj) )
            ! masked ratio alpha/beta
            zalbet = fsalbt( zt, zs, zh )*umask(ji,jj,1)
            ! local density gradient along i-bathymetric slope
            zgdrho = zalbet * ( ztnb(ji+1,jj) - ztnb(ji,jj) )   &
                   -          ( zsnb(ji+1,jj) - zsnb(ji,jj) )
            ! sign of local i-gradient of density multiplied by the i-slope
            zsign = SIGN( 0.5, - zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
            zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#  endif
      END DO

#  if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#  else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
#  endif
            ! temperature, salinity anomalie and depth
            zt = 0.5 * ( ztnb(ji,jj+1) + ztnb(ji,jj) )
            zs = 0.5 * ( zsnb(ji,jj+1) + zsnb(ji,jj) ) - 35.0
            zh = 0.5 * ( zdep(ji,jj+1) + zdep(ji,jj) )
            ! masked ratio alpha/beta
            zalbet = fsalbt( zt, zs, zh )*vmask(ji,jj,1)
            ! local density gradient along j-bathymetric slope
            zgdrho = zalbet * ( ztnb(ji,jj+1) - ztnb(ji,jj) )   &
                   -          ( zsnb(ji,jj+1) - zsnb(ji,jj) )
            ! sign of local j-gradient of density multiplied by the j-slope
            zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
            zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#  endif
      END DO

	CASE ( 1 )               ! Linear formulation function of temperature only

#  if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#  else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
#  endif
            ! local density gradient along i-bathymetric slope
            zgdrho =  ( ztnb(ji+1,jj) - ztnb(ji,jj) )
            ! sign of local i-gradient of density multiplied by the i-slope
            zsign = SIGN( 0.5, - zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
            zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#  endif
      END DO

#  if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#  else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
#  endif
            ! local density gradient along j-bathymetric slope
            zgdrho =  ( ztnb(ji,jj+1) - ztnb(ji,jj) )
            ! sign of local j-gradient of density multiplied by the j-slope
            zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
            zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)

#  if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#  endif
      END DO

      CASE ( 2 )               ! Linear formulation function of temperature and salinity

      DO jj = 1, jpjm1
        DO ji = 1, fs_jpim1   ! vector opt.
            ! local density gradient along i-bathymetric slope
            zgdrho = - ( rbeta*( zsnb(ji+1,jj) - zsnb(ji,jj) )   &
                     -  ralpha*( ztnb(ji+1,jj) - ztnb(ji,jj) ) )
            ! sign of local i-gradient of density multiplied by the i-slope
            zsign = SIGN( 0.5, - zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
 	    zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
        END DO
      END DO

      DO jj = 1, jpjm1
        DO ji = 1, fs_jpim1   ! vector opt.
            ! local density gradient along j-bathymetric slope
            zgdrho = - ( rbeta*( zsnb(ji,jj+1) - zsnb(ji,jj) )   &
                     -  ralpha*( ztnb(ji,jj+1) - ztnb(ji,jj) ) )
            ! sign of local j-gradient of density multiplied by the j-slope
            zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
            zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
         END DO
      END DO


      CASE DEFAULT

         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
         nstop = nstop + 1

      END SELECT
      
      ! 2. Additional second order diffusive trends
      ! -------------------------------------------

      DO jn = 1, jptra
         ! first derivative (gradient)

#  if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = 1, jpij   ! vector opt. (forced unrolling)
#  else
         DO jj = 1, jpj
            DO ji = 1, jpi
#  endif
               ik = mbkt(ji,jj)
               ztrb(ji,jj) = trb(ji,jj,ik,jn) * tmask(ji,jj,1)
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
#  endif
         END DO
#  if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#  else
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
#  endif
               zkx(ji,jj) = zki(ji,jj) * ( ztrb(ji+1,jj) - ztrb(ji,jj) )
               zky(ji,jj) = zkj(ji,jj) * ( ztrb(ji,jj+1) - ztrb(ji,jj) )
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
#  endif
         END DO
#endif

         IF( cp_cfg == "orca" ) THEN
            
            SELECT CASE ( jp_cfg )
               !                                           ! =======================
            CASE ( 2 )                                  !  ORCA_R2 configuration
               !                                        ! =======================
               ! Gibraltar enhancement of BBL
               ij0 = 102   ;   ij1 = 102
               ii0 = 139   ;   ii1 = 140  
               zkx( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 4.e0 * zkx( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) )
               zky( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 4.e0 * zky( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) )
               
               ! Red Sea enhancement of BBL
               ij0 =  88   ;   ij1 =  88
               ii0 = 161   ;   ii1 = 162
               zkx( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 10.e0 * zkx( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) )
               zky( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 10.e0 * zky( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) )
               
               !                                        ! =======================
            CASE ( 4 )                                  !  ORCA_R4 configuration
               !                                        ! =======================
               ! Gibraltar enhancement of BBL
               ij0 =  52   ;   ij1 =  52
               ii0 =  70   ;   ii1 =  71  
               zkx( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 4.e0 * zkx( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) )
               zky( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 4.e0 * zky( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) )
               
            END SELECT
            
         ENDIF
         
         ! second derivative (divergence) and add to the general tracer trend
#  if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
#  else
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
#  endif
               ik = MAX( mbathy(ji,jj)-1, 1 )
               zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,ik) )
               ztra = (  zkx(ji,jj) - zkx(ji-1,jj  )    &
                  &    + zky(ji,jj) - zky(ji  ,jj-1)  ) * zbtr
               tra(ji,jj,ik,jn) = tra(ji,jj,ik,jn) + ztra
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
#  endif
         END DO

      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bbl - dif')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_bbl_dif

# if defined key_trcbbl_adv
   !!----------------------------------------------------------------------
   !!   'key_trcbbl_adv'                    advective bottom boundary layer
   !!----------------------------------------------------------------------
#  include "trcbbl_adv.h90"
# else
   !!----------------------------------------------------------------------
   !!   Default option :                 NO advective bottom boundary layer
   !!----------------------------------------------------------------------
   SUBROUTINE trc_bbl_adv (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_bbl_adv: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bbl_adv
# endif

   SUBROUTINE trc_bbl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_bbl_init  ***
      !!
      !! ** Purpose :   Initialization for the bottom boundary layer scheme.
      !!
      !! ** Method  :   Read the namtrcbbl namelist and check the parameters
      !!      called by tra_bbl at the first timestep (nittrc000)
      !!
      !! History :
      !!    8.5  !  02-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj      ! dummy loop indices
      INTEGER :: numnat=80
      NAMELIST/namtrcbbl/ atrcbbl

      !!----------------------------------------------------------------------
      ! Read Namelist namtrcbbl : bottom boundary layer scheme
      ! --------------------

      OPEN(numnat,FILE='namelist.trp.cfc')
      REWIND ( numnat )
      READ   ( numnat, namtrcbbl )
      CLOSE(numnat)


      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_bbl_init : * Diffusive Bottom Boundary Layer'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) ' bottom boundary layer coef.    atrcbbl = ', atrcbbl
# if defined key_trcbbl_adv
            WRITE(numout,*) '               * Advective Bottom Boundary Layer'
# endif
         WRITE(numout,*)
      ENDIF
 
      DO jj = 1, jpj
         DO ji = 1, jpi
            mbkt(ji,jj) = MAX( mbathy(ji,jj) - 1, 1 )   ! vertical index of the bottom ocean T-level
         END DO
      END DO
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
            mbku(ji,jj) = MAX( MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) ) - 1, 1 )
            mbkv(ji,jj) = MAX( MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) ) - 1, 1 )
         END DO
      END DO
!!bug ???
!!bug Caution : define the vakue of mbku & mbkv everywhere!!! but lbc mpp lnk : pb when closed (0)

# if defined key_trcbbl_adv
      w_trc_bbl(:,:,:) = 0.e0    ! initialisation of w_trc_bbl to zero
# endif

   END SUBROUTINE trc_bbl_init

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      No bottom boundary layer scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbl_dif = .FALSE.   !: diff bbl flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbl_adv = .FALSE.   !: adv  bbl flag
CONTAINS
   SUBROUTINE trc_bbl_dif (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_bbl_dif: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bbl_dif
   SUBROUTINE trc_bbl_adv (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_bbl_adv: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bbl_adv
#endif

   !!======================================================================
END MODULE trcbbl
