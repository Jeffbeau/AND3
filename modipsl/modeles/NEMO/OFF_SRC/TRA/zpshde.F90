MODULE zpshde
   !!==============================================================================
   !!                       ***  MODULE zpshde   ***
   !! Ocean active tracers: 
   !!==============================================================================
#if defined key_partial_steps || defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_partial_steps' :               z-coordinate with partial steps
   !!----------------------------------------------------------------------
   !!   zps_hde      :  Horizontal DErivative of T, S and rd at the last
   !!                   ocean level (Z-coord. with Partial Steps)
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/TRA/zpshde.F90,v 1.2 2005/11/16 16:15:21 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space domain variables
   USE oce             ! ocean dynamics and tracers variables
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE eosbn2          ! ocean equation of state
   USE lbclnk          ! lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zps_hde          ! routine called by step.F90

   !! * module variables
   INTEGER, DIMENSION(jpi,jpj) ::   &
      mbatu, mbatv      ! bottom ocean level index at U- and V-points

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE zps_hde ( kt, ptem, psal, prd ,   &
                            pgtu, pgsu, pgru,   &
                            pgtv, pgsv, pgrv  )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE zps_hde  ***
      !!                    
      !! ** Purpose :   Compute the horizontal derivative of T, S and rd
      !!      at u- and v-points with a linear interpolation for z-coordinate
      !!      with partial steps.
      !!
      !! ** Method  :   In z-coord with partial steps, scale factors on last 
      !!      levels are different for each grid point, so that T, S and rd 
      !!      points are not at the same depth as in z-coord. To have horizontal
      !!      gradients again, we interpolate T and S at the good depth : 
      !!      Linear interpolation of T, S   
      !!         Computation of di(tb) and dj(tb) by vertical interpolation:
      !!          di(t) = t~ - t(i,j,k) or t(i+1,j,k) - t~
      !!          dj(t) = t~ - t(i,j,k) or t(i,j+1,k) - t~
      !!         This formulation computes the two cases:
      !!                 CASE 1                   CASE 2  
      !!         k-1  ___ ___________   k-1   ___ ___________
      !!                    Ti  T~                  T~  Ti+1
      !!                  _____                        _____
      !!         k        |   |Ti+1     k           Ti |   |
      !!                  |   |____                ____|   |
      !!              ___ |   |   |           ___  |   |   |
      !!                  
      !!      case 1->   e3w(i+1) >= e3w(i) ( and e3w(j+1) >= e3w(j) ) then
      !!          t~ = t(i+1,j  ,k) + (e3w(i+1) - e3w(i)) * dk(Ti+1)/e3w(i+1)
      !!        ( t~ = t(i  ,j+1,k) + (e3w(j+1) - e3w(j)) * dk(Tj+1)/e3w(j+1)  )
      !!          or
      !!      case 2->   e3w(i+1) <= e3w(i) ( and e3w(j+1) <= e3w(j) ) then
      !!          t~ = t(i,j,k) + (e3w(i) - e3w(i+1)) * dk(Ti)/e3w(i )
      !!        ( t~ = t(i,j,k) + (e3w(j) - e3w(j+1)) * dk(Tj)/e3w(j ) )
      !!          Idem for di(s) and dj(s)          
      !!
      !!      For rho, we call eos_insitu_2d which will compute rd~(t~,s~) at 
      !!      the good depth zh from interpolated T and S for the different
      !!      formulation of the equation of state (eos).
      !!      Gradient formulation for rho :
      !!          di(rho) = rd~ - rd(i,j,k) or rd (i+1,j,k) - rd~
      !!
      !! ** Action  : - pgtu, pgsu, pgru: horizontal gradient of T, S
      !!                and rd at U-points 
      !!              - pgtv, pgsv, pgrv: horizontal gradient of T, S
      !!                and rd at V-points 
      !!
      !! History :
      !!   8.5  !  02-04  (A. Bozec)  Original code
      !!   8.5  !  02-08  (G. Madec E. Durand)  Optimization and Free form
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
         ptem, psal, prd            ! 3D T, S and rd fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT( out ) ::   &
         pgtu, pgsu, pgru,       &  ! horizontal grad. of T, S and rd at u- 
         pgtv, pgsv, pgrv           ! and v-points of the partial step level

      !! * Local declarations
      INTEGER ::   ji, jj,       &  ! Dummy loop indices
                   iku,ikv          ! partial step level at u- and v-points
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zti, ztj, zsi, zsj,     &  ! interpolated value of T, S 
         zri, zrj,               &  ! and rd
         zhgi, zhgj                 ! depth of interpolation for eos2d
      REAL(wp) ::   &
         ze3wu, ze3wv,           &  ! temporary scalars
         zmaxu1, zmaxu2,         &  !    "         "
         zmaxv1, zmaxv2             !    "         "

      ! Initialization (first time-step only): compute mbatu and mbatv
      IF( kt == nit000 ) THEN
         mbatu(:,:) = 0
         mbatv(:,:) = 0
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               mbatu(ji,jj) = MAX( MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1, 2 )
               mbatv(ji,jj) = MAX( MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1, 2 )
            END DO
         END DO
         zti(:,:) = FLOAT( mbatu(:,:) )
         ztj(:,:) = FLOAT( mbatv(:,:) )
         ! lateral boundary conditions: T-point, sign unchanged
         CALL lbc_lnk( zti , 'U', 1. )
         CALL lbc_lnk( ztj , 'V', 1. )
         mbatu(:,:) = MAX( INT( zti(:,:) ), 2 )
         mbatv(:,:) = MAX( INT( ztj(:,:) ), 2 )
      ENDIF
      

      ! Interpolation of T and S at the last ocean level
# if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = 1, jpij-jpi   ! vector opt. (forced unrolled)
# else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
# endif
            ! last level
            iku = mbatu(ji,jj)
            ikv = mbatv(ji,jj)

            ze3wu  = fse3w(ji+1,jj  ,iku) - fse3w(ji,jj,iku)
            ze3wv  = fse3w(ji  ,jj+1,ikv) - fse3w(ji,jj,ikv)
            zmaxu1 =  ze3wu / fse3w(ji+1,jj  ,iku)
            zmaxu2 = -ze3wu / fse3w(ji  ,jj  ,iku)
            zmaxv1 =  ze3wv / fse3w(ji  ,jj+1,ikv)
            zmaxv2 = -ze3wv / fse3w(ji  ,jj  ,ikv)

            ! i- direction

            IF( ze3wu >= 0. ) THEN      ! case 1
               ! interpolated values of T and S
               zti(ji,jj) = ptem(ji+1,jj,iku) + zmaxu1 * ( ptem(ji+1,jj,iku-1) - ptem(ji+1,jj,iku) )
               zsi(ji,jj) = psal(ji+1,jj,iku) + zmaxu1 * ( psal(ji+1,jj,iku-1) - psal(ji+1,jj,iku) )
               ! depth of the partial step level
               zhgi(ji,jj) = fsdept(ji,jj,iku)
               ! gradient of T and S
               pgtu(ji,jj) = umask(ji,jj,1) * ( zti(ji,jj) - ptem(ji,jj,iku) )
               pgsu(ji,jj) = umask(ji,jj,1) * ( zsi(ji,jj) - psal(ji,jj,iku) )

            ELSE                        ! case 2
               ! interpolated values of T and S
               zti(ji,jj) = ptem(ji,jj,iku) + zmaxu2 * ( ptem(ji,jj,iku-1) - ptem(ji,jj,iku) )
               zsi(ji,jj) = psal(ji,jj,iku) + zmaxu2 * ( psal(ji,jj,iku-1) - psal(ji,jj,iku) )
               ! depth of the partial step level
               zhgi(ji,jj) = fsdept(ji+1,jj,iku)
               ! gradient of T and S 
               pgtu(ji,jj) = umask(ji,jj,1) * ( ptem(ji+1,jj,iku) - zti (ji,jj) )
               pgsu(ji,jj) = umask(ji,jj,1) * ( psal(ji+1,jj,iku) - zsi (ji,jj) )
            ENDIF

            ! j- direction

            IF( ze3wv >= 0. ) THEN      ! case 1
               ! interpolated values of T and S
               ztj(ji,jj) = ptem(ji,jj+1,ikv) + zmaxv1 * ( ptem(ji,jj+1,ikv-1) - ptem(ji,jj+1,ikv) )
               zsj(ji,jj) = psal(ji,jj+1,ikv) + zmaxv1 * ( psal(ji,jj+1,ikv-1) - psal(ji,jj+1,ikv) )
               ! depth of the partial step level
               zhgj(ji,jj) = fsdept(ji,jj,ikv) 
               ! gradient of T and S
               pgtv(ji,jj) = vmask(ji,jj,1) * ( ztj(ji,jj) - ptem(ji,jj,ikv) )
               pgsv(ji,jj) = vmask(ji,jj,1) * ( zsj(ji,jj) - psal(ji,jj,ikv) )

            ELSE                        ! case 2
               ! interpolated values of T and S
               ztj(ji,jj) = ptem(ji,jj,ikv) + zmaxv2 * ( ptem(ji,jj,ikv-1) - ptem(ji,jj,ikv) )
               zsj(ji,jj) = psal(ji,jj,ikv) + zmaxv2 * ( psal(ji,jj,ikv-1) - psal(ji,jj,ikv) ) 
               ! depth of the partial step level
               zhgj(ji,jj) = fsdept(ji,jj+1,ikv) 
               ! gradient of T and S
               pgtv(ji,jj) = vmask(ji,jj,1) * ( ptem(ji,jj+1,ikv) - ztj(ji,jj) )
               pgsv(ji,jj) = vmask(ji,jj,1) * ( psal(ji,jj+1,ikv) - zsj(ji,jj) )
            ENDIF
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO

      ! Compute interpolated rd from zti, zsi, ztj, zsj for the 2 cases at the depth of the partial
      ! step and store it in  zri, zrj for each  case
      CALL eos( zti, zsi, zhgi, zri )
      CALL eos( ztj, zsj, zhgj, zrj )


      ! Gradient of density at the last level 
# if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = 1, jpij-jpi   ! vector opt. (forced unrolled)
# else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
# endif
            iku = mbatu(ji,jj)
            ikv = mbatv(ji,jj)
            ze3wu  = fse3w(ji+1,jj  ,iku) - fse3w(ji,jj,iku)
            ze3wv  = fse3w(ji  ,jj+1,ikv) - fse3w(ji,jj,ikv)
            IF( ze3wu >= 0. ) THEN    ! i-direction: case 1
               pgru(ji,jj) = umask(ji,jj,1) * ( zri(ji,jj) - prd(ji,jj,iku) )
            ELSE                      ! i-direction: case 2
               pgru(ji,jj) = umask(ji,jj,1) * ( prd(ji+1,jj,iku) - zri(ji,jj) )
            ENDIF
            IF( ze3wv >= 0. ) THEN    ! j-direction: case 1
               pgrv(ji,jj) = vmask(ji,jj,1) * ( zrj(ji,jj) - prd(ji,jj,ikv) )  
            ELSE                      ! j-direction: case 2
               pgrv(ji,jj) = vmask(ji,jj,1) * ( prd(ji,jj+1,ikv) - zrj(ji,jj) )
            ENDIF
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO

      ! Lateral boundary conditions on each gradient
      CALL lbc_lnk( pgtu , 'U', -1. )   ;   CALL lbc_lnk( pgtv , 'V', -1. )
      CALL lbc_lnk( pgsu , 'U', -1. )   ;   CALL lbc_lnk( pgsv , 'V', -1. )
      CALL lbc_lnk( pgru , 'U', -1. )   ;   CALL lbc_lnk( pgrv , 'V', -1. )

   END SUBROUTINE zps_hde

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
   USE par_kind
CONTAINS
   SUBROUTINE zps_hde ( kt, ptem, psal, prd ,   &      ! Empty routine
                            pgtu, pgsu, pgru,   &
                            pgtv, pgsv, pgrv  )
      REAL(wp), DIMENSION(:,:,:) :: ptem, psal, prd
      REAL(wp) :: pgtu, pgsu, pgru, pgtv, pgsv, pgrv
      WRITE(*,*) 'zps_hde: You should not have seen this print! error?',   &
         kt, ptem, psal, prd, pgtu, pgsu, pgru, pgtv, pgsv, pgrv
   END SUBROUTINE zps_hde
#endif

   !!======================================================================
END MODULE zpshde
