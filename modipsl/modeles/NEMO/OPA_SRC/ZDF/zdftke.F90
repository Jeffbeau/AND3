MODULE zdftke
   !!======================================================================
   !!                       ***  MODULE  zdftke  ***
   !! Ocean physics:  vertical mixing coefficient compute from the tke 
   !!                 turbulent closure parameterization
   !!=====================================================================
#if defined key_zdftke   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_zdftke'                                             TKE scheme
   !!----------------------------------------------------------------------
   !!   zdf_tke      : update momentum and tracer Kz from a tke scheme
   !!   zdf_tke_init : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers 
   USE dom_oce         ! ocean space and time domain
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE taumod          ! surface stress
   USE prtctl          ! Print control
   USE traadv_ctl      ! advection scheme control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zdf_tke   ! routine called by step.F90

   !! * Share Module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdftke = .TRUE.    !: TKE vertical mixing flag
   LOGICAL, PUBLIC ::         & !!: ** tke namelist (namtke) **
     ln_rstke = .FALSE.          !: =T restart with tke from a run without tke with 
     !                           !  a none zero initial value for en
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      en                         !: now turbulent kinetic energy

   !! * Module variables
   INTEGER ::                 & !!! ** tke namelist (namtke) **
      nitke = 50 ,            &  ! number of restart iterative loops
      nmxl  =  2 ,            &  ! = 0/1/2/3 flag for the type of mixing length used
      npdl  =  1 ,            &  ! = 0/1/2 flag on prandtl number on vert. eddy coeff.
      nave  =  1 ,            &  ! = 0/1 flag for horizontal average on avt, avmu, avmv
      navb  =  0                 ! = 0/1 flag for constant or profile background avt
   REAL(wp) ::                & !!! ** tke namlist (namtke) **
      ediff = 0.1_wp       ,  &  ! coeff. for vertical eddy coef.; avt=ediff*mxl*sqrt(e)
      ediss = 0.7_wp       ,  &  ! coef. of the Kolmogoroff dissipation 
      ebb   = 3.75_wp      ,  &  ! coef. of the surface input of tke
      efave = 1._wp        ,  &  ! coef. for the tke vert. diff. coeff.; avtke=efave*avm
      emin  = 0.7071e-6_wp ,  &  ! minimum value of tke (m2/s2)
      emin0 = 1.e-4_wp     ,  &  ! surface minimum value of tke (m2/s2)
      ri_c  = 2._wp / 9._wp      ! critic Richardson number
   REAL(wp) ::   &  
      eboost        ! multiplicative coeff of the shear product.

   !! caution vectopt_memory change the solution (last digit of the solver stat)
# if defined key_vectopt_memory
   REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
      etmean,    &  ! coefficient used for horizontal smoothing
      eumean,    &  ! at t-, u- and v-points
      evmean        !
# endif

# if defined key_cfg_1d
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   
      e_dis,    &   ! dissipation turbulent lengh scale
      e_mix,    &   ! mixing turbulent lengh scale
      e_pdl,    &   ! prandl number
      e_ric         ! local Richardson number
#endif

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdftke.F90,v 1.9 2006/03/21 09:16:38 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

# if defined key_autotasking
   !!----------------------------------------------------------------------
   !!   'key_autotasking' :                             j-k-i loop (j-slab)
   !!----------------------------------------------------------------------
#  include "zdftke_atsk.h90"

# else
   !!----------------------------------------------------------------------
   !!   default option :                                         k-j-i loop
   !!----------------------------------------------------------------------

   SUBROUTINE zdf_tke ( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_tke  ***
      !!
      !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
      !!      coefficients using a 1.5 turbulent closure scheme.
      !!
      !! ** Method  :   The time evolution of the turbulent kinetic energy 
      !!      (tke) is computed from a prognostic equation :
      !!         d(en)/dt = eboost eav (d(u)/dz)**2       ! shear production
      !!                  + d( efave eav d(en)/dz )/dz    ! diffusion of tke
      !!                  + grav/rau0 pdl eav d(rau)/dz   ! stratif. destruc.
      !!                  - ediss / emxl en**(2/3)        ! dissipation
      !!      with the boundary conditions:
      !!         surface: en = max( emin0,ebb sqrt(taux^2 + tauy^2) )
      !!         bottom : en = emin
      !!      -1- The dissipation and mixing turbulent lengh scales are computed
      !!      from the usual diagnostic buoyancy length scale:  
      !!         mxl= 1/(sqrt(en)/N)  WHERE N is the brunt-vaisala frequency
      !!      Four cases : 
      !!         nmxl=0 : mxl bounded by the distance to surface and bottom.
      !!                  zmxld = zmxlm = mxl
      !!         nmxl=1 : mxl bounded by the vertical scale factor.
      !!                  zmxld = zmxlm = mxl
      !!         nmxl=2 : mxl bounded such that the vertical derivative of mxl 
      !!                  is less than 1 (|d/dz(xml)|<1). 
      !!                  zmxld = zmxlm = mxl
      !!         nmxl=3 : lup = mxl bounded using |d/dz(xml)|<1 from the surface 
      !!                        to the bottom
      !!                  ldown = mxl bounded using |d/dz(xml)|<1 from the bottom 
      !!                        to the surface
      !!                  zmxld = sqrt (lup*ldown) ; zmxlm = min(lup,ldown)
      !!      -2- Compute the now Turbulent kinetic energy. The time differencing
      !!      is implicit for vertical diffusion term, linearized for kolmo-
      !!      goroff dissipation term, and explicit forward for both buoyancy
      !!      and dynamic production terms. Thus a tridiagonal linear system is
      !!      solved.
      !!         Note that - the shear production is multiplied by eboost in order
      !!      to set the critic richardson number to ri_c (namelist parameter)
      !!                   - the destruction by stratification term is multiplied
      !!      by the Prandtl number (defined by an empirical funtion of the local 
      !!      Richardson number) if npdl=1 (namelist parameter)
      !!      coefficient (zesh2):
      !!      -3- Compute the now vertical eddy vicosity and diffusivity
      !!      coefficients from en (before the time stepping) and zmxlm:
      !!              avm = max( avtb, ediff*zmxlm*en^1/2 )
      !!              avt = max( avmb, pdl*avm )  (pdl=1 if npdl=0)
      !!              eav = max( avmb, avm )
      !!      avt and avm are horizontally averaged to avoid numerical insta-
      !!      bilities.
      !!        N.B. The computation is done from jk=2 to jpkm1 except for
      !!      en. Surface value of avt avmu avmv are set once a time to
      !!      their background value in routine zdf_tke_init.
      !!
      !! ** Action  :   compute en (now turbulent kinetic energy)
      !!                update avt, avmu, avmv (before vertical eddy coef.)
      !!
      !! References :
      !!      Gaspar et al., jgr, 95, 1990,
      !!      Blanke and Delecluse, jpo, 1991
      !! History :
      !!   6.0  !  91-03  (b. blanke)  Original code
      !!   7.0  !  91-11  (G. Madec)   bug fix
      !!   7.1  !  92-10  (G. Madec)   new mixing length and eav
      !!   7.2  !  93-03  (M. Guyon)   symetrical conditions
      !!   7.3  !  94-08  (G. Madec, M. Imbard)   npdl flag
      !!   7.5  !  96-01  (G. Madec)   s-coordinates
      !!   8.0  !  97-07  (G. Madec)   lbc
      !!   8.1  !  99-01  (E. Stretta) new option for the mixing length
      !!   8.5  !  02-08  (G. Madec)  ri_c and Free form, F90
      !!   9.0  !  04-10  (C. Ethe )  1D configuration
      !!----------------------------------------------------------------------
      !! * Modules used
      USE oce     , zwd   => ua,  &  ! use ua as workspace
         &          zmxlm => ta,  &  ! use ta as workspace
         &          zmxld => sa      ! use sa as workspace

      !! * arguments
      INTEGER, INTENT( in  ) ::   kt ! ocean time step

      !! * local declarations
      INTEGER ::   ji, jj, jk        ! dummy loop arguments
      REAL(wp) ::   &
         zmlmin, zbbrau,          &  ! temporary scalars
         zfact1, zfact2, zfact3,  &  !
         zrn2, zesurf,            &  !
         ztx2, zty2, zav,         &  !
         zcoef, zcof, zsh2,       &  !
         zdku, zdkv, zpdl, zri,   &  !
         zsqen, zesh2,            &  !
         zemxl, zemlm, zemlp
      !!--------------------------------------------------------------------

      ! Initialization (first time-step only)
      ! --------------
      IF( kt == nit000  )   CALL zdf_tke_init

      ! Local constant initialization
      zmlmin = 1.e-8
      zbbrau =  .5 * ebb / rau0
      zfact1 = -.5 * rdt * efave
      zfact2 = 1.5 * rdt * ediss
      zfact3 = 0.5 * rdt * ediss


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! I.  Mixing length
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


      ! Buoyancy length scale: l=sqrt(2*e/n**2)
      ! ---------------------
      zmxlm(:,:, 1 ) = zmlmin   ! surface set to the minimum value
      zmxlm(:,:,jpk) = zmlmin   ! bottom  set to the minimum value
!CDIR NOVERRCHK
      DO jk = 2, jpkm1
!CDIR NOVERRCHK
         DO jj = 2, jpjm1
!CDIR NOVERRCHK
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zrn2 = MAX( rn2(ji,jj,jk), rsmall )
               zmxlm(ji,jj,jk) = MAX( SQRT( 2. * en(ji,jj,jk) / zrn2 ), zmlmin  )
            END DO
         END DO
      END DO


      ! Physical limits for the mixing length
      ! -------------------------------------
      zmxld(:,:, 1 ) = zmlmin   ! surface set to the minimum value
      zmxld(:,:,jpk) = zmlmin   ! bottom  set to the minimum value

      SELECT CASE ( nmxl )

      CASE ( 0 )           ! bounded by the distance to surface and bottom
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zemxl = MIN( fsdepw(ji,jj,jk), zmxlm(ji,jj,jk),   &
                  &            fsdepw(ji,jj,mbathy(ji,jj)) - fsdepw(ji,jj,jk) )
                  zmxlm(ji,jj,jk) = zemxl
                  zmxld(ji,jj,jk) = zemxl
               END DO
            END DO
         END DO

      CASE ( 1 )           ! bounded by the vertical scale factor
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zemxl = MIN( fse3w(ji,jj,jk), zmxlm(ji,jj,jk) )
                  zmxlm(ji,jj,jk) = zemxl
                  zmxld(ji,jj,jk) = zemxl
               END DO
            END DO
         END DO

      CASE ( 2 )           ! |dk[xml]| bounded by e3t :
         DO jk = 2, jpkm1         ! from the surface to the bottom :
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zmxlm(ji,jj,jk) = MIN( zmxlm(ji,jj,jk-1) + fse3t(ji,jj,jk-1), zmxlm(ji,jj,jk) )
               END DO
            END DO
         END DO
         DO jk = jpkm1, 2, -1     ! from the bottom to the surface :
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zemxl = MIN( zmxlm(ji,jj,jk+1) + fse3t(ji,jj,jk+1), zmxlm(ji,jj,jk) )
                  zmxlm(ji,jj,jk) = zemxl
                  zmxld(ji,jj,jk) = zemxl
               END DO
            END DO
         END DO

      CASE ( 3 )           ! lup and ldown, |dk[xml]| bounded by e3t :
         DO jk = 2, jpkm1         ! from the surface to the bottom : lup
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zmxld(ji,jj,jk) = MIN( zmxld(ji,jj,jk-1) + fse3t(ji,jj,jk-1), zmxlm(ji,jj,jk) )
               END DO
            END DO
         END DO
         DO jk = jpkm1, 2, -1     ! from the bottom to the surface : ldown
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zmxlm(ji,jj,jk) = MIN( zmxlm(ji,jj,jk+1) + fse3t(ji,jj,jk+1), zmxlm(ji,jj,jk) )
               END DO
            END DO
         END DO
!CDIR NOVERRCHK
         DO jk = 2, jpkm1
!CDIR NOVERRCHK
            DO jj = 2, jpjm1
!CDIR NOVERRCHK
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zemlm = MIN ( zmxld(ji,jj,jk),  zmxlm(ji,jj,jk) )
                  zemlp = SQRT( zmxld(ji,jj,jk) * zmxlm(ji,jj,jk) )
                  zmxlm(ji,jj,jk) = zemlm
                  zmxld(ji,jj,jk) = zemlp
               END DO
            END DO
         END DO

      END SELECT

# if defined key_cfg_1d
      ! save mixing and dissipation turbulent length scales
      e_dis(:,:,:) = zmxld(:,:,:)
      e_mix(:,:,:) = zmxlm(:,:,:)
# endif


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! II  Tubulent kinetic energy time stepping
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


      ! 1. Vertical eddy viscosity on tke (put in zmxlm) and first estimate of avt
      ! ---------------------------------------------------------------------
!CDIR NOVERRCHK
      DO jk = 2, jpkm1
!CDIR NOVERRCHK
         DO jj = 2, jpjm1
!CDIR NOVERRCHK
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zsqen = SQRT( en(ji,jj,jk) )
               zav   = ediff * zmxlm(ji,jj,jk) * zsqen
!DL, october 2013, modifications to increase mixing in upper SL estuary
!DL, find local indices for global i=1:29, j=143:182
               if ((mig(ji).ge.12.AND.mig(ji).lt.30).AND.(mjg(jj).gt.150.AND. mjg(jj).lt.173)) then
                  avt  (ji,jj,jk) = MAX( zav, avtb(jk) )*6.0 * tmask(ji,jj,jk)
                  !avt  (ji,jj,jk) =( zav+ avtb(jk)*50.0 ) * tmask(ji,jj,jk)
               else
                  avt  (ji,jj,jk) = MAX( zav, avtb(jk) ) * tmask(ji,jj,jk) !original line
               endif
!DL End modif
               !avt  (ji,jj,jk) = MAX( zav, avtb(jk) ) * tmask(ji,jj,jk)
               zmxlm(ji,jj,jk) = MAX( zav, avmb(jk) ) * tmask(ji,jj,jk)
               zmxld(ji,jj,jk) = zsqen / zmxld(ji,jj,jk)
            END DO
         END DO
      END DO

      ! 2. Surface boundary condition on tke and its eddy viscosity (zmxlm)
      ! -------------------------------------------------
      ! en(1)   = ebb sqrt(taux^2+tauy^2) / rau0  (min value emin0)
      ! zmxlm(1) = avmb(1) and zmxlm(jpk) = 0.
!CDIR NOVERRCHK
      DO jj = 2, jpjm1
!CDIR NOVERRCHK
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ztx2 = taux(ji-1,jj  ) + taux(ji,jj)
            zty2 = tauy(ji  ,jj-1) + tauy(ji,jj)
            zesurf = zbbrau * SQRT( ztx2 * ztx2 + zty2 * zty2 )
            en (ji,jj,1) = MAX( zesurf, emin0 ) * tmask(ji,jj,1)
            zmxlm(ji,jj,1  ) = avmb(1) * tmask(ji,jj,1)
            zmxlm(ji,jj,jpk) = 0.e0
         END DO
      END DO

      ! 3. Now Turbulent kinetic energy (output in en)
      ! -------------------------------
      ! Resolution of a tridiagonal linear system by a "methode de chasse"
      ! computation from level 2 to jpkm1  (e(1) already computed and
      ! e(jpk)=0 ). 

      SELECT CASE ( npdl )

      CASE ( 0 )           ! No Prandtl number
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! zesh2 = eboost * (du/dz)^2 - N^2
                  zcoef = 0.5 / fse3w(ji,jj,jk)
                  ! shear
                  zdku = zcoef * (   ub(ji-1, jj ,jk-1) + ub(ji,jj,jk-1)   &
                  &                - ub(ji-1, jj ,jk  ) - ub(ji,jj,jk  )  )
                  zdkv = zcoef * (   vb( ji ,jj-1,jk-1) + vb(ji,jj,jk-1)   &
                  &                - vb( ji ,jj-1,jk  ) - vb(ji,jj,jk  )  )
                  ! coefficient (zesh2)
                  zesh2 =  eboost * ( zdku*zdku + zdkv*zdkv ) - rn2(ji,jj,jk)

                  ! Matrix
                  zcof = zfact1 * tmask(ji,jj,jk)
                  ! lower diagonal
                  avmv(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk  ) + zmxlm(ji,jj,jk-1) )   &
                  &                    / ( fse3t(ji,jj,jk-1) * fse3w(ji,jj,jk  ) )
                  ! upper diagonal
                  avmu(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk+1) + zmxlm(ji,jj,jk  ) )   &
                  &                    / ( fse3t(ji,jj,jk  ) * fse3w(ji,jj,jk) )
                  ! diagonal
                  zwd(ji,jj,jk) = 1. - avmv(ji,jj,jk) - avmu(ji,jj,jk) + zfact2 * zmxld(ji,jj,jk)
                  ! right hand side in en 
                  en(ji,jj,jk) = en(ji,jj,jk) + zfact3 * zmxld(ji,jj,jk) * en   (ji,jj,jk)   &
                  &                           +   rdt  * zmxlm(ji,jj,jk) * zesh2
               END DO
            END DO
         END DO

      CASE ( 1 )           ! Prandtl number
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! zesh2 =  eboost * (du/dz)^2 - pdl * N^2
                  zcoef = 0.5 / fse3w(ji,jj,jk)
                  ! shear
                  zdku = zcoef * (   ub(ji-1,jj  ,jk-1) + ub(ji,jj,jk-1) &
                  &                - ub(ji-1,jj  ,jk  ) - ub(ji,jj,jk  )   )
                  zdkv = zcoef * (   vb(ji  ,jj-1,jk-1) + vb(ji,jj,jk-1) &
                  &                - vb(ji  ,jj-1,jk  ) - vb(ji,jj,jk  )   )
                  ! square of vertical shear
                  zsh2 = zdku * zdku + zdkv * zdkv
                  ! local Richardson number
                  zri  = MAX( rn2(ji,jj,jk), 0. ) / ( zsh2 + 1.e-20 )
# if defined key_cfg_1d
                  ! save masked local Richardson number in zmxlm array
                  e_ric(ji,jj,jk) = zri * tmask(ji,jj,jk)
# endif
                  ! Prandtl number
                  zpdl = 1.0
                  IF( zri >= 0.2 ) zpdl = 0.2 / zri
                  zpdl = MAX( 0.1, zpdl )
                  ! coefficient (esh2)
                  zesh2 = eboost * zsh2 - zpdl * rn2(ji,jj,jk)

                  ! Matrix
                  zcof = zfact1 * tmask(ji,jj,jk)
                  ! lower diagonal
                  avmv(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk  ) + zmxlm(ji,jj,jk-1) )   &
                  &                     / ( fse3t(ji,jj,jk-1) * fse3w(ji,jj,jk  ) )
                  ! upper diagonal
                  avmu(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk+1) + zmxlm(ji,jj,jk  ) )   &
                  &                     / ( fse3t(ji,jj,jk  ) * fse3w(ji,jj,jk) )
                  ! diagonal
                  zwd(ji,jj,jk) = 1. - avmv(ji,jj,jk) - avmu(ji,jj,jk) + zfact2 * zmxld(ji,jj,jk)
                  ! right hand side in en 
                  en(ji,jj,jk) = en(ji,jj,jk) + zfact3 * zmxld(ji,jj,jk) * en   (ji,jj,jk)   &
                  &                           +   rdt  * zmxlm(ji,jj,jk) * zesh2
                  ! save masked Prandlt number in zmxlm array
                  zmxld(ji,jj,jk) = zpdl * tmask(ji,jj,jk)
               END DO
            END DO
         END DO

      END SELECT

# if defined key_cfg_1d
      !  save masked Prandlt number
      e_pdl(:,:,2:jpkm1) = zmxld(:,:,2:jpkm1)
      e_pdl(:,:,      1) = e_pdl(:,:,      2)
      e_pdl(:,:,    jpk) = e_pdl(:,:,  jpkm1)      
# endif

      ! 4. Matrix inversion from level 2 (tke prescribed at level 1)
      !!--------------------------------
      ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1
      DO jk = 3, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - avmv(ji,jj,jk) * avmu(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! Second recurrence : Lk = RHSk - Lk / Dk-1 * Lk-1
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1    ! vector opt.
            avmv(ji,jj,2) = en(ji,jj,2) - avmv(ji,jj,2) * en(ji,jj,1)    ! Surface boudary conditions on tke
         END DO
      END DO
      DO jk = 3, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               avmv(ji,jj,jk) = en(ji,jj,jk) - avmv(ji,jj,jk) / zwd(ji,jj,jk-1) *avmv(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! thrid recurrence : Ek = ( Lk - Uk * Ek+1 ) / Dk
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1    ! vector opt.
            en(ji,jj,jpkm1) = avmv(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 2, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1    ! vector opt.
               en(ji,jj,jk) = ( avmv(ji,jj,jk) - avmu(ji,jj,jk) * en(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Save the result in en and set minimum value of tke : emin
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               en(ji,jj,jk) = MAX( en(ji,jj,jk), emin ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Lateral boundary conditions on ( avt, en )   (sign unchanged)
      CALL lbc_lnk( en , 'W', 1. )   ;   CALL lbc_lnk( avt, 'W', 1. )


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! III.  Before vertical eddy vicosity and diffusivity coefficients
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      SELECT CASE ( nave )
         
      CASE ( 0 )                ! no horizontal average

         ! Vertical eddy viscosity

         DO jk = 2, jpkm1                                 ! Horizontal slab
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  avmu(ji,jj,jk) = ( avt  (ji,jj,jk) + avt  (ji+1,jj  ,jk) ) * umask(ji,jj,jk)   &
                  &     / MAX( 1.,   tmask(ji,jj,jk) + tmask(ji+1,jj  ,jk) )
                  avmv(ji,jj,jk) = ( avt  (ji,jj,jk) + avt  (ji  ,jj+1,jk) ) * vmask(ji,jj,jk)   &
                  &     / MAX( 1.,   tmask(ji,jj,jk) + tmask(ji  ,jj+1,jk) )
               END DO
            END DO
         END DO

!DB: jack up viscosity in shallow water 
!         do jj = 2, jpjm1
!            do ji = fs_2, fs_jpim1   ! vector opt.
!               if(mbathy(ji,jj) <= 5) then
!                  avmu(ji,jj,:) = 10. * avmu(ji,jj,:)
!                  avmv(ji,jj,:) = 10. * avmv(ji,jj,:)
!               endif
!            enddo
!         enddo
!!DB: END

         ! Lateral boundary conditions (avmu,avmv) (U- and V- points, sign unchanged)
         CALL lbc_lnk( avmu, 'U', 1. )   ;    CALL lbc_lnk( avmv, 'V', 1. )
         
      CASE ( 1 )                ! horizontal average

         !                                                ( 1/2  1/2 )
         ! Eddy viscosity: horizontal average: avmu = 1/4 ( 1    1   )
         !                      ( 1/2  1 1/2 )            ( 1/2  1/2 )
         !           avmv = 1/4 ( 1/2  1 1/2 )      
         
!! caution vectopt_memory change the solution (last digit of the solver stat)
#  if defined key_vectopt_memory
         DO jk = 2, jpkm1                                 ! Horizontal slab
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  avmu(ji,jj,jk) = (      avt(ji,jj  ,jk) + avt(ji+1,jj  ,jk)   &
                  &                 +.5*( avt(ji,jj-1,jk) + avt(ji+1,jj-1,jk)   &
                  &                      +avt(ji,jj+1,jk) + avt(ji+1,jj+1,jk) ) ) * eumean(ji,jj,jk)

                  avmv(ji,jj,jk) = (      avt(ji  ,jj,jk) + avt(ji  ,jj+1,jk)   &
                  &                 +.5*( avt(ji-1,jj,jk) + avt(ji-1,jj+1,jk)   &
                  &                      +avt(ji+1,jj,jk) + avt(ji+1,jj+1,jk) ) ) * evmean(ji,jj,jk)
               END DO
            END DO
         END DO
!DB: jack up viscosity in shallow water
!         do jj = 2, jpjm1
!            do ji = fs_2, fs_jpim1   ! vector opt.
!               if(mbathy(ji,jj) <= 5) then
!                  avmu(ji,jj,:) = 10. * avmu(ji,jj,:)
!                  avmv(ji,jj,:) = 10. * avmv(ji,jj,:)
!               endif
!            enddo
!         enddo
!!DB: END

#  else
         DO jk = 2, jpkm1                                 ! Horizontal slab
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  avmu(ji,jj,jk) = (   avt  (ji,jj  ,jk) + avt  (ji+1,jj  ,jk)   &
                  &              +.5*( avt  (ji,jj-1,jk) + avt  (ji+1,jj-1,jk)   &
                  &                   +avt  (ji,jj+1,jk) + avt  (ji+1,jj+1,jk) ) ) * umask(ji,jj,jk)  &
                  &       / MAX( 1.,   tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk)   &
                  &              +.5*( tmask(ji,jj-1,jk) + tmask(ji+1,jj-1,jk)   &
                  &                   +tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) )  )

                  avmv(ji,jj,jk) = (   avt  (ji  ,jj,jk) + avt  (ji  ,jj+1,jk)   &
                  &              +.5*( avt  (ji-1,jj,jk) + avt  (ji-1,jj+1,jk)   &
                  &                   +avt  (ji+1,jj,jk) + avt  (ji+1,jj+1,jk) ) ) * vmask(ji,jj,jk)  &
                  &      /  MAX( 1.,   tmask(ji  ,jj,jk) + tmask(ji  ,jj+1,jk)   &
                  &              +.5*( tmask(ji-1,jj,jk) + tmask(ji-1,jj+1,jk)   &
                  &                   +tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk) )  )
               END DO
            END DO
         END DO

!DB: jack up viscosity in shallow water
!         do jj = 2, jpjm1
!            do ji = fs_2, fs_jpim1   ! vector opt.
!               if(mbathy(ji,jj) <= 5) then
!                  avmu(ji,jj,:) = 10. * avmu(ji,jj,:)
!                  avmv(ji,jj,:) = 10. * avmv(ji,jj,:)
!               endif
!            enddo
!         enddo
!!DB: END

#  endif

         ! Lateral boundary conditions (avmu,avmv) (sign unchanged)
         CALL lbc_lnk( avmu, 'U', 1. )   ;    CALL lbc_lnk( avmv, 'V', 1. )

         ! Vertical eddy diffusivity
         ! ------------------------------
         !                                (1 2 1)
         ! horizontal average  avt = 1/16 (2 4 2)
         !                                (1 2 1)
         DO jk = 2, jpkm1                                 ! Horizontal slab
#  if defined key_vectopt_memory
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  avt(ji,jj,jk) = ( avmu(ji,jj,jk) + avmu(ji-1,jj  ,jk)    &
                  &               + avmv(ji,jj,jk) + avmv(ji  ,jj-1,jk)  ) * etmean(ji,jj,jk)
               END DO
            END DO
#  else
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  avt(ji,jj,jk) = ( avmu (ji,jj,jk) + avmu (ji-1,jj  ,jk)   &
                  &               + avmv (ji,jj,jk) + avmv (ji  ,jj-1,jk)  ) * tmask(ji,jj,jk)   &
                  &     / MAX( 1.,  umask(ji,jj,jk) + umask(ji-1,jj  ,jk)   &
                  &               + vmask(ji,jj,jk) + vmask(ji  ,jj-1,jk)  )
               END DO
            END DO
#  endif
         END DO

      END SELECT

      ! multiplied by the Prandtl number (npdl>1)
      ! ----------------------------------------
      IF( npdl == 1 ) THEN
         DO jk = 2, jpkm1                                 ! Horizontal slab
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zpdl = zmxld(ji,jj,jk)
                  avt(ji,jj,jk) = MAX( zpdl * avt(ji,jj,jk), avtb(jk) ) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! Minimum value on the eddy viscosity
      ! ----------------------------------------
      DO jk = 2, jpkm1                                 ! Horizontal slab
         DO jj = 1, jpj
            DO ji = 1, jpi
               avmu(ji,jj,jk) = MAX( avmu(ji,jj,jk), avmb(jk) ) * umask(ji,jj,jk)
               avmv(ji,jj,jk) = MAX( avmv(ji,jj,jk), avmb(jk) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Lateral boundary conditions on avt  (sign unchanged)
      ! ------------------------------=====
      CALL lbc_lnk( avt, 'W', 1. )

      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=en  , clinfo1=' tke  - e: ', tab3d_2=avt , clinfo2=' t: ', ovlap=1, kdim=jpk)
         CALL prt_ctl(tab3d_1=avmu, clinfo1=' tke  - u: ', tab3d_2=avmv, clinfo2=' v: ', ovlap=1, kdim=jpk)
      ENDIF

   END SUBROUTINE zdf_tke

# endif

   SUBROUTINE zdf_tke_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_tke_init  ***
      !!                     
      !! ** Purpose :   Initialization of the vertical eddy diffivity and 
      !!      viscosity when using a tke turbulent closure scheme
      !!
      !! ** Method  :   Read the namtke namelist and check the parameters
      !!      called at the first timestep (nit000)
      !!
      !! ** input   :   Namlist namtke
      !!
      !! ** Action  :   Increase by 1 the nstop flag is setting problem encounter
      !!
      !! history :
      !!  8.5  ! 02-06 (G. Madec) original code
      !!----------------------------------------------------------------------
      !! * Module used
      USE dynzdf_exp
      USE trazdf_exp

      !! * local declarations
      !! caution vectopt_memory change the solution (last digit of the solver stat)
# if defined key_vectopt_memory
      INTEGER ::   ji, jj, jk, jit   ! dummy loop indices
# else
      INTEGER ::           jk, jit   ! dummy loop indices
# endif

      NAMELIST/namtke/ ln_rstke, ediff, ediss, ebb, efave, emin, emin0,   &
         ri_c, nitke, nmxl, npdl, nave, navb
      !!----------------------------------------------------------------------

      ! Read Namelist namtke : Turbulente Kinetic Energy
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, namtke )

      ! Compute boost associated with the Richardson critic
      !     (control values: ri_c = 0.3   ==> eboost=1.25 for npdl=1 or 2)
      !     (                ri_c = 0.222 ==> eboost=1.                  )
      eboost = ri_c * ( 2. + ediss / ediff ) / 2.


      ! Parameter control and print
      ! ---------------------------
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_tke_init : tke turbulent closure scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist namtke : set tke mixing parameters'
         WRITE(numout,*) '             restart with tke from no tke ln_rstke = ', ln_rstke
         WRITE(numout,*) '             coef. to compute avt           ediff  = ', ediff
         WRITE(numout,*) '             Kolmogoroff dissipation coef.  ediss  = ', ediss
         WRITE(numout,*) '             tke surface input coef.        ebb    = ', ebb
         WRITE(numout,*) '             tke diffusion coef.            efave  = ', efave
         WRITE(numout,*) '             minimum value of tke           emin   = ', emin
         WRITE(numout,*) '             surface minimum value of tke   emin0  = ', emin0
         WRITE(numout,*) '             number of restart iter loops   nitke  = ', nitke
         WRITE(numout,*) '             mixing length type             nmxl   = ', nmxl
         WRITE(numout,*) '             prandl number flag             npdl   = ', npdl
         WRITE(numout,*) '             horizontal average flag        nave   = ', nave
         WRITE(numout,*) '             critic Richardson nb           ri_c   = ', ri_c
         WRITE(numout,*) '                and its associated coeff.   eboost = ', eboost
         WRITE(numout,*) '             constant background or profile navb   = ', navb
         WRITE(numout,*)
      ENDIF

      ! Check nmxl and npdl values
      IF( nmxl < 0 .OR. nmxl > 3 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag: nmxl is < 0 or > 3 '
         nstop = nstop + 1
      ENDIF
      IF ( npdl < 0 .OR. npdl > 1 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag: npdl is < 0 or > 1 '
         nstop = nstop + 1
      ENDIF


      ! Horizontal average : initialization of weighting arrays 
      ! -------------------
      
      SELECT CASE ( nave )

      CASE ( 0 )                ! no horizontal average
         IF(lwp) WRITE(numout,*) '          no horizontal average on avt, avmu, avmv'
         IF(lwp) WRITE(numout,*) '          only in very high horizontal resolution !'
!! caution vectopt_memory change the solution (last digit of the solver stat)
# if defined key_vectopt_memory
         ! weighting mean arrays etmean, eumean and evmean
         !           ( 1  1 )                                          ( 1 )
         ! avt = 1/4 ( 1  1 )     avmu = 1/2 ( 1  1 )       avmv=  1/2 ( 1 )
         !                          
         etmean(:,:,:) = 0.e0
         eumean(:,:,:) = 0.e0
         evmean(:,:,:) = 0.e0
         
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  etmean(ji,jj,jk) = tmask(ji,jj,jk)                     &
                  &  / MAX( 1.,  umask(ji-1,jj  ,jk) + umask(ji,jj,jk)   &
                  &            + vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)  )
                  
                  eumean(ji,jj,jk) = umask(ji,jj,jk)                     &
                  &  / MAX( 1.,  tmask(ji,jj,jk) + tmask(ji+1,jj  ,jk)  )

                  evmean(ji,jj,jk) = vmask(ji,jj,jk)                     &
                  &  / MAX( 1.,  tmask(ji,jj,jk) + tmask(ji  ,jj+1,jk)  )
               END DO
            END DO
         END DO
# endif

      CASE ( 1 )                ! horizontal average 
         IF(lwp) WRITE(numout,*) '          horizontal average on avt, avmu, avmv'
!! caution vectopt_memory change the solution (last digit of the solver stat)
# if defined key_vectopt_memory
         ! weighting mean arrays etmean, eumean and evmean
         !           ( 1  1 )              ( 1/2  1/2 )             ( 1/2  1  1/2 )
         ! avt = 1/4 ( 1  1 )   avmu = 1/4 ( 1    1   )   avmv= 1/4 ( 1/2  1  1/2 )
         !                                 ( 1/2  1/2 )
         etmean(:,:,:) = 0.e0
         eumean(:,:,:) = 0.e0
         evmean(:,:,:) = 0.e0
         
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  etmean(ji,jj,jk) = tmask(ji,jj,jk)                     &
                  &  / MAX( 1.,  umask(ji-1,jj  ,jk) + umask(ji,jj,jk)   &
                  &            + vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)  )
                  
                  eumean(ji,jj,jk) = umask(ji,jj,jk)                        &
                  &  / MAX( 1.,   tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk)   &
                  &       +.5 * ( tmask(ji,jj-1,jk) + tmask(ji+1,jj-1,jk)   &
                  &              +tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) )  )

                  evmean(ji,jj,jk) = vmask(ji,jj,jk)                        &
                  &  / MAX( 1.,   tmask(ji  ,jj,jk) + tmask(ji  ,jj+1,jk)   &
                  &       +.5 * ( tmask(ji-1,jj,jk) + tmask(ji-1,jj+1,jk)   &
                  &              +tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk) )  )
               END DO
            END DO
         END DO
# endif

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for nave = ', nave
         nstop = nstop + 1

      END SELECT


      ! Background eddy viscosity and diffusivity profil
      ! ------------------------------------------------
      IF( navb == 0 ) THEN
         !   Define avmb, avtb from namelist parameter
         avmb(:) = avm0
         avtb(:) = avt0
      ELSE
         !   Background profile of avt (fit a theoretical/observational profile (Krauss 1990) 
         avmb(:) = avm0
         avtb(:) = avt0 + ( 3.0e-4 - 2 * avt0 ) * 1.0e-4 * gdepw(:)   ! m2/s
      ENDIF

!!DB: delete ORCA
!      IF( cp_cfg == "orca"  .AND. jp_cfg == 2 .AND. ln_traadv_cen2 )  THEN

      ! Initialization of vertical eddy coef. to the background value
      ! -------------------------------------------------------------
      DO jk = 1, jpk
         avt (:,:,jk) = avtb(jk) * tmask(:,:,jk)
         avmu(:,:,jk) = avmb(jk) * umask(:,:,jk)
         avmv(:,:,jk) = avmb(jk) * vmask(:,:,jk)
      END DO


      ! Initialization of turbulent kinetic energy ( en )
      ! -------------------------------------------------
      IF( ln_rstart ) THEN
         ! no en field in the restart file, en set by iterative loop
         IF( ln_rstke ) THEN
            en (:,:,:) = emin * tmask(:,:,:)
            DO jit = 2, nitke+1
               CALL zdf_tke( jit )
            END DO
         ENDIF
         ! otherwise en is already read in dtrlec called by inidtr
      ELSE
         ! no restart: en set to emin
         en(:,:,:) = emin * tmask(:,:,:)
      ENDIF

   END SUBROUTINE zdf_tke_init

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        NO TKE scheme
   !!----------------------------------------------------------------------
   USE in_out_manager
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdftke = .FALSE.   !: TKE flag
CONTAINS
   SUBROUTINE zdf_tke( kt )          ! Empty routine
      if(lwp) WRITE(numout,*) 'zdf_tke: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_tke
#endif

   !!======================================================================
END MODULE zdftke
