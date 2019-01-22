MODULE lbclnk
   !!======================================================================
   !!                       ***  MODULE  lbclnk  ***
   !! Ocean        : lateral boundary conditions
   !!=====================================================================
#if   defined key_mpp_mpi   ||   defined key_mpp_shmem
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'     OR      MPI massively parallel processing library
   !!   'key_mpp_shmem'         SHMEM massively parallel processing library
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   lbc_lnk      : generic interface for mpp_lnk_3d and mpp_lnk_2d
   !!                  routines defined in lib_mpp
   !!   lbc_lnk_e    : generic interface for mpp_lnk_2d_e
   !!                   routinee defined in lib_mpp
   !!----------------------------------------------------------------------
   !! * Modules used
   USE lib_mpp          ! distributed memory computing library

   INTERFACE lbc_lnk
      MODULE PROCEDURE mpp_lnk_3d, mpp_lnk_2d
   END INTERFACE

   INTERFACE lbc_lnk_e
      MODULE PROCEDURE mpp_lnk_2d_e
   END INTERFACE

   PUBLIC lbc_lnk       ! ocean lateral boundary conditions
   PUBLIC lbc_lnk_e
   !!----------------------------------------------------------------------

#else
   !!----------------------------------------------------------------------
   !!   Default option                              shared memory computing
   !!----------------------------------------------------------------------
   !!   lbc_lnk      : generic interface for lbc_lnk_3d and lbc_lnk_2d
   !!   lbc_lnk_3d   : set the lateral boundary condition on a 3D variable
   !!                  on OPA ocean mesh
   !!   lbc_lnk_2d   : set the lateral boundary condition on a 2D variable
   !!                  on OPA ocean mesh
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_lnk
      MODULE PROCEDURE lbc_lnk_3d, lbc_lnk_2d
   END INTERFACE

   INTERFACE lbc_lnk_e
      MODULE PROCEDURE lbc_lnk_2d
   END INTERFACE

   PUBLIC lbc_lnk       ! ocean/ice  lateral boundary conditions
   PUBLIC  lbc_lnk_e 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lbc_lnk_3d( pt3d, cd_type, psgn )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_lnk_3d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions (non mpp case)
      !!
      !! ** Method  :
      !!
      !! History :
      !!        !  97-06  (G. Madec)  Original code
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1), INTENT( in ) ::   &
         cd_type       ! nature of pt3d grid-points
         !             !   = T ,  U , V , F or W  gridpoints
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
         pt3d          ! 3D array on which the boundary condition is applied
      REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
         !             !   =-1 , the sign is changed if north fold boundary
         !             !   = 1 , no sign change
         !             !   = 0 , no sign change and > 0 required (use the inner
         !             !         row/column if closed boundary)

      !! * Local declarations
      INTEGER  ::   ji, jk
      INTEGER  ::   ijt, iju
      !!----------------------------------------------------------------------
      !!  OPA 9.0 , LOCEAN-IPSL (2005) 
      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/lbclnk.F90,v 1.2 2005/11/16 16:19:34 opalod Exp $ 
      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
      !!----------------------------------------------------------------------
      
      !                                                      ! ===============
      DO jk = 1, jpk                                         ! Horizontal slab
         !                                                   ! ===============

         !                                     ! East-West boundaries
         !                                     ! ====================
         SELECT CASE ( nperio )

         CASE ( 1 , 4 , 6 )                    ! * cyclic east-west
            pt3d( 1 ,:,jk) = pt3d(jpim1,:,jk)          ! all points
            pt3d(jpi,:,jk) = pt3d(  2  ,:,jk)

         CASE DEFAULT                          ! * closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt3d( 1 ,:,jk) = 0.e0
               pt3d(jpi,:,jk) = 0.e0
            CASE ( 'F' )                               ! F-point
               pt3d(jpi,:,jk) = 0.e0
            END SELECT

         END SELECT

         !                                     ! North-South boundaries
         !                                     ! ======================
         SELECT CASE ( nperio )

         CASE ( 2 )                            ! *  south symmetric

            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'W' )                   ! T-, U-, W-points
               pt3d(:, 1 ,jk) = pt3d(:,3,jk)
               pt3d(:,jpj,jk) = 0.e0
            CASE ( 'V' , 'F' )                         ! V-, F-points
               pt3d(:, 1 ,jk) = psgn * pt3d(:,2,jk)
               pt3d(:,jpj,jk) = 0.e0
            END SELECT

         CASE ( 3 , 4 )                        ! *  North fold  T-point pivot

            pt3d( 1 ,jpj,jk) = 0.e0
            pt3d(jpi,jpj,jk) = 0.e0

            SELECT CASE ( cd_type )
            CASE ( 'T' , 'W' )                         ! T-, W-point
               DO ji = 2, jpi
                  ijt = jpi-ji+2
                  pt3d(ji, 1 ,jk) = 0.e0
                  pt3d(ji,jpj,jk) = psgn * pt3d(ijt,jpj-2,jk)
               END DO
               DO ji = jpi/2+1, jpi
                  ijt = jpi-ji+2
                  pt3d(ji,jpjm1,jk) = psgn * pt3d(ijt,jpjm1,jk)
               END DO
            CASE ( 'U' )                               ! U-point
               DO ji = 1, jpi-1
                  iju = jpi-ji+1
                  pt3d(ji, 1 ,jk) = 0.e0
                  pt3d(ji,jpj,jk) = psgn * pt3d(iju,jpj-2,jk)
               END DO
               DO ji = jpi/2, jpi-1
                  iju = jpi-ji+1
                  pt3d(ji,jpjm1,jk) = psgn * pt3d(iju,jpjm1,jk)
               END DO
            CASE ( 'V' )                               ! V-point
                  DO ji = 2, jpi
                     ijt = jpi-ji+2
                     pt3d(ji,  1  ,jk) = 0.e0
                     pt3d(ji,jpj-1,jk) = psgn * pt3d(ijt,jpj-2,jk)
                     pt3d(ji,jpj  ,jk) = psgn * pt3d(ijt,jpj-3,jk)
                  END DO
            CASE ( 'F' )                               ! F-point
                  DO ji = 1, jpi-1
                     iju = jpi-ji+1
                     pt3d(ji,jpj-1,jk) = psgn * pt3d(iju,jpj-2,jk)
                     pt3d(ji,jpj  ,jk) = psgn * pt3d(iju,jpj-3,jk)
                  END DO
            END SELECT

         CASE ( 5 , 6 )                        ! *  North fold  F-point pivot

            pt3d( 1 ,jpj,jk) = 0.e0
            pt3d(jpi,jpj,jk) = 0.e0

            SELECT CASE ( cd_type )
            CASE ( 'T' , 'W' )                         ! T-, W-point
               DO ji = 1, jpi
                  ijt = jpi-ji+1
                  pt3d(ji, 1 ,jk) = 0.e0
                  pt3d(ji,jpj,jk) = psgn * pt3d(ijt,jpj-1,jk)
               END DO
            CASE ( 'U' )                               ! U-point
                  DO ji = 1, jpi-1
                     iju = jpi-ji
                     pt3d(ji, 1 ,jk) = 0.e0
                     pt3d(ji,jpj,jk) = psgn * pt3d(iju,jpj-1,jk)
                  END DO
            CASE ( 'V' )                               ! V-point
                  DO ji = 1, jpi
                     ijt = jpi-ji+1
                     pt3d(ji, 1 ,jk) = 0.e0
                     pt3d(ji,jpj,jk) = psgn * pt3d(ijt,jpj-2,jk)
                  END DO
                  DO ji = jpi/2+1, jpi
                     ijt = jpi-ji+1
                     pt3d(ji,jpjm1,jk) = psgn * pt3d(ijt,jpjm1,jk)
                  END DO
            CASE ( 'F' )                               ! F-point
                  DO ji = 1, jpi-1
                     iju = jpi-ji
                     pt3d(ji,jpj  ,jk) = psgn * pt3d(iju,jpj-2,jk)
                  END DO
                  DO ji = jpi/2+1, jpi-1
                     iju = jpi-ji
                     pt3d(ji,jpjm1,jk) = psgn * pt3d(iju,jpjm1,jk)
                  END DO
            END SELECT

         CASE DEFAULT                          ! *  closed

            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt3d(:, 1 ,jk) = 0.e0
               pt3d(:,jpj,jk) = 0.e0
            CASE ( 'F' )                               ! F-point
               pt3d(:,jpj,jk) = 0.e0
            END SELECT

         END SELECT
         !                                                   ! ===============
      END DO                                                 !   End of slab
      !                                                      ! ===============
   END SUBROUTINE lbc_lnk_3d


   SUBROUTINE lbc_lnk_2d( pt2d, cd_type, psgn )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk_2d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions (non mpp case)
      !!
      !! ** Method  :
      !!
      !! History :
      !!        !  97-06  (G. Madec)  Original code
      !!        !  01-05  (E. Durand)  correction
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1), INTENT( in ) ::   &
         cd_type       ! nature of pt2d grid-point
         !             !   = T , U , V , F or W  gridpoints
         !             !   = I sea-ice U-V gridpoint (= F ocean grid point with indice shift)
      REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
         !             !   =-1 , the sign is modified following the type of b.c. used
         !             !   = 1 , no sign change
      REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
         pt2d          ! 2D array on which the boundary condition is applied

      !! * Local declarations
      INTEGER  ::   ji
      INTEGER  ::   ijt, iju
            
      !                                        ! East-West boundaries
      !                                        ! ====================
      SELECT CASE ( nperio )

      CASE ( 1 , 4 , 6 )                       ! * cyclic east-west
         pt2d( 1 ,:) = pt2d(jpim1,:)
         pt2d(jpi,:) = pt2d(  2  ,:)

      CASE DEFAULT                             ! * closed 
         SELECT CASE ( cd_type )
         CASE ( 'T' , 'U' , 'V' , 'W' )                ! T-, U-, V-, W-points
            pt2d( 1 ,:) = 0.e0
            pt2d(jpi,:) = 0.e0
         CASE ( 'F' )                                  ! F-point, ice U-V point
            pt2d(jpi,:) = 0.e0 
         CASE ( 'I' )                                  ! F-point, ice U-V point
            pt2d( 1 ,:) = 0.e0 
            pt2d(jpi,:) = 0.e0 
         END SELECT

      END SELECT

      !                                        ! North-South boundaries
      !                                        ! ======================
      SELECT CASE ( nperio )

      CASE ( 2 )                               ! * South symmetric

         SELECT CASE ( cd_type )
         CASE ( 'T' , 'U' , 'W' )                      ! T-, U-, W-points
            pt2d(:, 1 ) = pt2d(:,3)
            pt2d(:,jpj) = 0.e0
         CASE ( 'V' , 'F' , 'I' )                      ! V-, F-points, ice U-V point
            pt2d(:, 1 ) = psgn * pt2d(:,2)
            pt2d(:,jpj) = 0.e0
         END SELECT

      CASE ( 3 , 4 )                           ! * North fold  T-point pivot

         pt2d( 1 , 1 ) = 0.e0        !!!!!  bug gm ??? !Edmee
         pt2d( 1 ,jpj) = 0.e0
         pt2d(jpi,jpj) = 0.e0

         SELECT CASE ( cd_type )

         CASE ( 'T' , 'W' )                            ! T-, W-point
            DO ji = 2, jpi
               ijt = jpi-ji+2
               pt2d(ji, 1 ) = 0.e0
               pt2d(ji,jpj) = psgn * pt2d(ijt,jpj-2)
            END DO
            DO ji = jpi/2+1, jpi
               ijt = jpi-ji+2
               pt2d(ji,jpjm1) = psgn * pt2d(ijt,jpjm1)
            END DO

         CASE ( 'U' )                                  ! U-point
            DO ji = 1, jpi-1
               iju = jpi-ji+1
               pt2d(ji, 1 ) = 0.e0
               pt2d(ji,jpj) = psgn * pt2d(iju,jpj-2)
            END DO
            DO ji = jpi/2, jpi-1
               iju = jpi-ji+1
               pt2d(ji,jpjm1) = psgn * pt2d(iju,jpjm1)
            END DO

         CASE ( 'V' )                                  ! V-point
            DO ji = 2, jpi
               ijt = jpi-ji+2
               pt2d(ji, 1   ) = 0.e0
               pt2d(ji,jpj-1) = psgn * pt2d(ijt,jpj-2)
               pt2d(ji,jpj  ) = psgn * pt2d(ijt,jpj-3)
            END DO

         CASE ( 'F' )                                  ! F-point
            DO ji = 1, jpi-1
               iju = jpi - ji + 1
               pt2d(ji,jpj-1) = psgn * pt2d(iju,jpj-2)
               pt2d(ji,jpj  ) = psgn * pt2d(iju,jpj-3)
            END DO

         CASE ( 'I' )                                  ! ice U-V point
            pt2d(:, 1 ) = 0.e0
            pt2d(2,jpj) = psgn * pt2d(3,jpj-1)
            DO ji = 3, jpi
               iju = jpi - ji + 3
               pt2d(ji,jpj) = psgn * pt2d(iju,jpj-1)
            END DO

         END SELECT

      CASE ( 5 , 6 )                           ! * North fold  F-point pivot

         pt2d( 1 , 1 ) = 0.e0           !!bug  ???
         pt2d( 1 ,jpj) = 0.e0
         pt2d(jpi,jpj) = 0.e0

         SELECT CASE ( cd_type )

         CASE ( 'T' , 'W' )                            ! T-, W-point
            DO ji = 1, jpi
               ijt = jpi-ji+1
               pt2d(ji, 1 ) = 0.e0
               pt2d(ji,jpj) = psgn * pt2d(ijt,jpj-1)
            END DO

         CASE ( 'U' )                                  ! U-point
            DO ji = 1, jpi-1
               iju = jpi-ji
               pt2d(ji, 1 ) = 0.e0
               pt2d(ji,jpj) = psgn * pt2d(iju,jpj-1)
            END DO

         CASE ( 'V' )                                  ! V-point
            DO ji = 1, jpi
               ijt = jpi-ji+1
               pt2d(ji, 1 ) = 0.e0
               pt2d(ji,jpj) = psgn * pt2d(ijt,jpj-2)
            END DO
            DO ji = jpi/2+1, jpi
               ijt = jpi-ji+1
               pt2d(ji,jpjm1) = psgn * pt2d(ijt,jpjm1)
            END DO

         CASE ( 'F' )                                  ! F-point
            DO ji = 1, jpi-1
               iju = jpi-ji
               pt2d(ji,jpj  ) = psgn * pt2d(iju,jpj-2)
            END DO
            DO ji = jpi/2+1, jpi-1
               iju = jpi-ji
               pt2d(ji,jpjm1) = psgn * pt2d(iju,jpjm1)
            END DO

         CASE ( 'I' )                                  ! ice U-V point
            pt2d( : , 1 ) = 0.e0
            pt2d( 2 ,jpj) = 0.e0
            DO ji = 2 , jpim1
               ijt = jpi - ji + 2
               pt2d(ji,jpj)= 0.5 * ( pt2d(ji,jpjm1) + psgn * pt2d(ijt,jpjm1) )
            END DO

         END SELECT

      CASE DEFAULT                             ! * closed

         SELECT CASE ( cd_type )
         CASE ( 'T' , 'U' , 'V' , 'W' )                ! T-, U-, V-, W-points
            pt2d(:, 1 ) = 0.e0
            pt2d(:,jpj) = 0.e0
         CASE ( 'F' )                                  ! F-point
            pt2d(:,jpj) = 0.e0
         CASE ( 'I' )                                  ! ice U-V point
            pt2d(:, 1 ) = 0.e0
            pt2d(:,jpj) = 0.e0
         END SELECT

      END SELECT

   END SUBROUTINE lbc_lnk_2d

#endif

   !!======================================================================
END MODULE lbclnk
