MODULE zpshde_trc
   !!==============================================================================
   !!                       ***  MODULE zpshde_trc   ***
   !! Ocean passive tracers: 
   !!==============================================================================
#if defined key_passivetrc && ( defined key_partial_steps || defined key_esopa )
   !!----------------------------------------------------------------------
   !!   'key_partial_steps' :               z-coordinate with partial steps
   !!----------------------------------------------------------------------
   !!   zps_hde_trc  :  Horizontal DErivative of passive tracers at the last
   !!                   ocean level (Z-coord. with Partial Steps)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc         ! ocean dynamics and tracers variables
   USE trc             ! ocean passive tracers variables
   USE lbclnk          ! lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zps_hde_trc          ! routine called by step.F90

   !! * module variables
   INTEGER, DIMENSION(jpi,jpj) ::   &
      mbatu, mbatv      ! bottom ocean level index at U- and V-points

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/zpshde_trc.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE zps_hde_trc ( kt, ptra, pgtru, pgtrv )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE zps_hde_trc  ***
      !!                    
      !! ** Purpose :   Compute the horizontal derivative of passive tracers
      !!      TRA at u- and v-points with a linear interpolation for z-coordinate
      !!      with partial steps.
      !!
      !! ** Method  :   In z-coord with partial steps, scale factors on last 
      !!      levels are different for each grid point, so that TRA points 
      !!      are not at the same depth as in z-coord. To have horizontal
      !!      gradients again, we interpolate TRA at the good depth : 
      !!      Linear interpolation of TRA  
      !!         Computation of di(trb) and dj(trb) by vertical interpolation:
      !!          di(tra) = tra~ - tra(i,j,k) or tra(i+1,j,k) - tra~
      !!          dj(tra) = tra~ - tra(i,j,k) or tra(i,j+1,k) - tra~
      !!         This formulation computes the two cases:
      !!                 CASE 1                   CASE 2  
      !!         k-1  ___ ___________   k-1   ___ ___________
      !!                  TRAi  TRA~             TRA~  TRAi+1
      !!                  _____                        _____
      !!         k        |   |TRAi+1   k         TRAi |   |
      !!                  |   |____                ____|   |
      !!              ___ |   |   |           ___  |   |   |
      !!                  
      !!      case 1->   e3w(i+1) >= e3w(i) ( and e3w(j+1) >= e3w(j) ) then
      !!      tra~ = tra(i+1,j  ,k) + (e3w(i+1) - e3w(i)) * dk(TRAi+1)/e3w(i+1)
      !!    ( tra~ = tra(i  ,j+1,k) + (e3w(j+1) - e3w(j)) * dk(TRAj+1)/e3w(j+1))
      !!          or
      !!      case 2->   e3w(i+1) <= e3w(i) ( and e3w(j+1) <= e3w(j) ) then
      !!       tra~ = tra(i,j,k) + (e3w(i) - e3w(i+1)) * dk(TRAi)/e3w(i )
      !!     ( tra~ = tra(i,j,k) + (e3w(j) - e3w(j+1)) * dk(TRAj)/e3w(j ) )
      !!      
      !!
      !! ** Action  : - pgtru : horizontal gradient of TRA at U-points 
      !!              - pgtrv : horizontal gradient of TRA at V-points 
      !!
      !! History :
      !!   8.5  !  02-04  (A. Bozec)  Original code
      !!   8.5  !  02-08  (G. Madec E. Durand)  Optimization and Free form
      !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT( in ) ::   &
         ptra                       ! 4D tracers fields
      REAL(wp), DIMENSION(jpi,jpj,jptra), INTENT( out ) ::   &
         pgtru,                 &  ! horizontal grad. of TRA u- and v-points 
         pgtrv                     ! of the partial step level

      !! * Local declarations
      INTEGER ::   ji, jj,jn,     &  ! Dummy loop indices
                   iku,ikv          ! partial step level at u- and v-points
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zti, ztj                   ! tempory arrays

      REAL(wp), DIMENSION(jpi,jpj,jptra) ::   &
         ztrai, ztraj               ! interpolated value of TRA

      REAL(wp) ::   &
         ze3wu, ze3wv,           &  ! temporary scalars
         zmaxu1, zmaxu2,         &  !    "         "
         zmaxv1, zmaxv2             !    "         "
      !!----------------------------------------------------------------------

      ! Initialization (first time-step only): compute mbatu and mbatv
      IF( kt == nittrc000 ) THEN
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
      

      DO jn = 1, jptra
         ! Interpolation of passive tracers at the last ocean level
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
                  ! interpolated values of passive tracers
                  ztrai(ji,jj,jn) = ptra(ji+1,jj,iku,jn) + zmaxu1 * ( ptra(ji+1,jj,iku-1,jn) - ptra(ji+1,jj,iku,jn) )
                  ! gradient of passive tracers
                  pgtru(ji,jj,jn) = umask(ji,jj,1) * ( ztrai(ji,jj,jn) - ptra(ji,jj,iku,jn) )
               ELSE                        ! case 2
                  ! interpolated values of passive tracers
                  ztrai(ji,jj,jn) = ptra(ji,jj,iku,jn) + zmaxu2 * ( ptra(ji,jj,iku-1,jn) - ptra(ji,jj,iku,jn) )
                  ! gradient of passive tracers
                  pgtru(ji,jj,jn) = umask(ji,jj,1) * ( ptra(ji+1,jj,iku,jn) - ztrai (ji,jj,jn) )
               ENDIF

               ! j- direction

               IF( ze3wv >= 0. ) THEN      ! case 1
                  ! interpolated values of passive tracers
                  ztraj(ji,jj,jn) = ptra(ji,jj+1,ikv,jn) + zmaxv1 * ( ptra(ji,jj+1,ikv-1,jn) - ptra(ji,jj+1,ikv,jn) )
                  ! gradient of passive tracers
                  pgtrv(ji,jj,jn) = vmask(ji,jj,1) * ( ztraj(ji,jj,jn) - ptra(ji,jj,ikv,jn) )
               ELSE                        ! case 2
                  ! interpolated values of passive tracers
                  ztraj(ji,jj,jn) = ptra(ji,jj,ikv,jn) + zmaxv2 * ( ptra(ji,jj,ikv-1,jn) - ptra(ji,jj,ikv,jn) )
                  ! gradient of passive tracers
                  pgtrv(ji,jj,jn) = vmask(ji,jj,1) * ( ptra(ji,jj+1,ikv,jn) - ztraj(ji,jj,jn) )
               ENDIF
# if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
# endif
         END DO

         ! Lateral boundary conditions on each gradient
         CALL lbc_lnk( pgtru(:,:,jn) , 'U', -1. ) 
         CALL lbc_lnk( pgtrv(:,:,jn) , 'V', -1. )

      END DO

   END SUBROUTINE zps_hde_trc

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
   USE par_kind
CONTAINS
   SUBROUTINE zps_hde_trc ( kt, ptra, pgtru, pgtrv ) ! Empty routine
      INTEGER, INTENT( in) :: kt
      REAL(wp), DIMENSION(:,:,:,:) :: ptra
      REAL(wp) :: pgtru, pgtrv
!      WRITE(*,*) 'zps_hde_trc: You should not have seen this print! error?',   &
!         kt, ptra, pgtru, pgtrv
   END SUBROUTINE zps_hde_trc
#endif

   !!======================================================================
END MODULE zpshde_trc
