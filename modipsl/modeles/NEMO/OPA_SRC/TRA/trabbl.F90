MODULE trabbl
   !!==============================================================================
   !!                       ***  MODULE  trabbl  ***
   !! Ocean physics :  advective and/or diffusive bottom boundary layer scheme
   !!==============================================================================
#if   defined key_trabbl_dif   ||   defined key_trabbl_adv   || defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_trabbl_dif'   or            diffusive bottom boundary layer
   !!   'key_trabbl_adv'                 advective bottom boundary layer
   !!----------------------------------------------------------------------
   !!   tra_bbl_dif  : update the active tracer trends due to the bottom
   !!                  boundary layer (diffusive only)
   !!   tra_bbl_adv  : update the active tracer trends due to the bottom
   !!                  boundary layer (advective and/or diffusive)
   !!   tra_bbl_init : initialization, namlist read, parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce                  ! ocean dynamics and active tracers
   USE dom_oce              ! ocean space and time domain
   USE trdmod_oce           ! ocean variables trends
   USE in_out_manager       ! I/O manager
   USE prtctl               ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_bbl_dif    ! routine called by step.F90
   PUBLIC tra_bbl_adv    ! routine called by step.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::            &  !!: * bbl namelist *
      atrbbl = 1.e+3                  !: lateral coeff. for BBL scheme (m2/s)
#if defined key_trabbl_dif
   LOGICAL, PUBLIC, PARAMETER ::   &  !:
      lk_trabbl_dif = .TRUE.          !: diffusive bottom boundary layer flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   &  !:
      lk_trabbl_dif = .FALSE.         !: diffusive bottom boundary layer flag
#endif

# if defined key_trabbl_adv
   LOGICAL, PUBLIC, PARAMETER ::    &  !:
      lk_trabbl_adv = .TRUE.   !: advective bottom boundary layer flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
       u_bbl, v_bbl,  &  !: velocity involved in exhanges in the advective BBL
       w_bbl             !: vertical increment of velocity due to advective BBL
       !                 !  only affect tracer vertical advection
# else
   LOGICAL, PUBLIC, PARAMETER ::    &  !:
      lk_trabbl_adv = .FALSE.  !: advective bottom boundary layer flag
# endif

   !! * Module variables
   INTEGER, DIMENSION(jpi,jpj) ::   &  !:
      mbkt,           &   ! vertical index of the bottom ocean T-level
      mbku, mbkv          ! vertical index of the bottom ocean U/V-level

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trabbl.F90,v 1.10 2006/03/27 09:55:49 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_bbl_dif( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_bbl_dif  ***
      !!
      !! ** Purpose :   Compute the before tracer (t & s) trend associated 
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
      !!      Add this before trend to the general trend (ta,sa) of the 
      !!      botton ocean tracer point:
      !!         ta = ta + difft
      !!
      !! ** Action  : - update (ta,sa) at the bottom level with the bottom
      !!                boundary layer trend
      !!              - save the trends in tldfbbl/sldfbbl ('key_trdtra')
      !!
      !! References :
      !!     Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!
      !! History :
      !!   8.0  !  96-06  (L. Mortier)  Original code
      !!   8.0  !  97-11  (G. Madec)  Optimization
      !!   8.5  !  02-08  (G. Madec)  free form + modules
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdta => ua,     & ! use ua as 3D workspace   
                         ztdsa => va        ! use va as 3D workspace   
      USE eosbn2, ONLY : neos ! type of equation of state

      !! * Arguments 
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      INTEGER ::   ji, jj                   ! dummy loop indices
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
         zgdrho, zbtr, zta, zsa
      REAL(wp), DIMENSION(jpi,jpj) ::    &
        zki, zkj, zkw, zkx, zky, zkz,    &  ! temporary workspace arrays
        ztnb, zsnb, zdep,                &
        ztbb, zsbb, zahu, zahv
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

      IF( kt == nit000 )   CALL tra_bbl_init

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

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
            ztbb(ji,jj) = tb(ji,jj,ik) * tmask(ji,jj,1)   ! masked before T and S at ocean bottom 
            zsbb(ji,jj) = sb(ji,jj,ik) * tmask(ji,jj,1)
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
            zahu(ji,jj) = atrbbl * e2u(ji,jj) * ze3u / e1u(ji,jj) * umask(ji,jj,1)
            zahv(ji,jj) = atrbbl * e1v(ji,jj) * ze3v / e2v(ji,jj) * vmask(ji,jj,1)
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
            zahu(ji,jj) = atrbbl * e2u(ji,jj) * fse3u(ji,jj,iku) / e1u(ji,jj) * umask(ji,jj,1)
            zahv(ji,jj) = atrbbl * e1v(ji,jj) * fse3v(ji,jj,ikv) / e2v(ji,jj) * vmask(ji,jj,1)
#   if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#   endif
      END DO
#  endif

      ! 1. Criteria of additional bottom diffusivity: grad(rho).grad(h)<0
      ! --------------------------------------------
      ! Sign of the local density gradient along the i- and j-slopes
      ! multiplied by the slope of the ocean bottom

      SELECT CASE ( neos )

      CASE ( 0   )               ! 0 :Jackett and McDougall (1994) formulation

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
                               ! 
#  if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#  else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
#  endif
            ! local 'density/temperature' gradient along i-bathymetric slope
            zgdrho =  ztnb(ji+1,jj) - ztnb(ji,jj) 
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
            zgdrho =  ztnb(ji,jj+1) - ztnb(ji,jj) 
            ! sign of local j-gradient of density multiplied by the j-slope
            zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
            zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#  endif
      END DO

      CASE ( 2 )               ! Linear formulation function of temperature and salinity

         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          use of linear eos rho(T,S) = rau0 * ( rbeta * S - ralpha * T )'
         IF(lwp) WRITE(numout,*) '          bbl not implented: easy to do it '
         nstop = nstop + 1

      CASE DEFAULT

         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
         nstop = nstop + 1

      END SELECT

      ! 2. Additional second order diffusive trends
      ! -------------------------------------------

      ! first derivative (gradient)
#  if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
#  else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
#  endif
            zkx(ji,jj) = zki(ji,jj) * ( ztbb(ji+1,jj) - ztbb(ji,jj) )
            zkz(ji,jj) = zki(ji,jj) * ( zsbb(ji+1,jj) - zsbb(ji,jj) )

            zky(ji,jj) = zkj(ji,jj) * ( ztbb(ji,jj+1) - ztbb(ji,jj) )
            zkw(ji,jj) = zkj(ji,jj) * ( zsbb(ji,jj+1) - zsbb(ji,jj) )
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#  endif
      END DO

!!DB: deleted ORCA
!      IF( cp_cfg == "orca" ) THEN
      ! second derivative (divergence) and add to the general tracer trend


#  if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
#  else
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
#  endif
            ik = max( mbathy(ji,jj)-1, 1 )
            zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,ik) )
            zta = (  zkx(ji,jj) - zkx(ji-1,jj  )    &
                   + zky(ji,jj) - zky(ji  ,jj-1)  ) * zbtr
            zsa = (  zkz(ji,jj) - zkz(ji-1,jj  )    &
                   + zkw(ji,jj) - zkw(ji  ,jj-1)  ) * zbtr
            ta(ji,jj,ik) = ta(ji,jj,ik) + zta
            sa(ji,jj,ik) = sa(ji,jj,ik) + zsa
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
#  endif
      END DO

      ! save the trends for diagnostic
      ! BBL lateral diffusion tracers trends
      IF( l_trdtra )   THEN
#  if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
#  else
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
#  endif
            ik = max( mbathy(ji,jj)-1, 1 )
            tldfbbl(ji,jj) = ta(ji,jj,ik) - ztdta(ji,jj,ik)
            sldfbbl(ji,jj) = sa(ji,jj,ik) - ztdsa(ji,jj,ik)
#  if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
#  endif
         END DO

      ENDIF

      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=ta, clinfo1=' bbl  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

   END SUBROUTINE tra_bbl_dif

# if defined key_trabbl_adv
   !!----------------------------------------------------------------------
   !!   'key_trabbl_adv'                    advective bottom boundary layer
   !!----------------------------------------------------------------------
#  include "trabbl_adv.h90"
# else
   !!----------------------------------------------------------------------
   !!   Default option :                 NO advective bottom boundary layer
   !!----------------------------------------------------------------------
   SUBROUTINE tra_bbl_adv (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'tra_bbl_adv: You should not have seen this print! error?', kt
   END SUBROUTINE tra_bbl_adv
# endif

   SUBROUTINE tra_bbl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_init  ***
      !!
      !! ** Purpose :   Initialization for the bottom boundary layer scheme.
      !!
      !! ** Method  :   Read the nambbl namelist and check the parameters
      !!      called by tra_bbl at the first timestep (nit000)
      !!
      !! History :
      !!    8.5  !  02-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj      ! dummy loop indices
      NAMELIST/nambbl/ atrbbl
      !!----------------------------------------------------------------------

      ! Read Namelist nambbl : bottom boundary layer scheme
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, nambbl )


      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_bbl_init : '
         WRITE(numout,*) '~~~~~~~~~~~~'
         IF (lk_trabbl_dif ) THEN
            WRITE(numout,*) '               * Diffusive Bottom Boundary Layer'
         ENDIF 
         IF( lk_trabbl_adv ) THEN
            WRITE(numout,*) '               * Advective Bottom Boundary Layer'
         ENDIF
         WRITE(numout,*) '          Namelist nambbl : set bbl parameters'
         WRITE(numout,*)
         WRITE(numout,*) '          bottom boundary layer coef.    atrbbl = ', atrbbl
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

# if defined key_trabbl_adv
      ! initialisation of w_bbl to zero
      w_bbl(:,:,:) = 0.e0    
# endif

   END SUBROUTINE tra_bbl_init

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      No bottom boundary layer scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl_dif = .FALSE.   !: diff bbl flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl_adv = .FALSE.   !: adv  bbl flag
CONTAINS
   SUBROUTINE tra_bbl_dif (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'tra_bbl_dif: You should not have seen this print! error?', kt
   END SUBROUTINE tra_bbl_dif
   SUBROUTINE tra_bbl_adv (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'tra_bbl_adv: You should not have seen this print! error?', kt
   END SUBROUTINE tra_bbl_adv
#endif

   !!======================================================================
END MODULE trabbl
