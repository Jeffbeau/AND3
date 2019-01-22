MODULE limrhg
   !!======================================================================
   !!                     ***  MODULE  limrhg  ***
   !!   Ice rheology :  performs sea ice rheology
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim'                                     LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_rhg   : computes ice velocities
   !!----------------------------------------------------------------------
   !! * Modules used
   USE phycst
   USE par_oce
   USE ice_oce         ! ice variables
   USE dom_ice
   USE ice
   USE lbclnk
   USE lib_mpp
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC lim_rhg  ! routine called by lim_dyn

   !! * Module variables
   REAL(wp)  ::           &  ! constant values
      rzero   = 0.e0   ,  &
      rone    = 1.e0
   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limrhg.F90,v 1.7 2006/03/10 10:35:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_rhg( k_j1, k_jpj )
      !!-------------------------------------------------------------------
      !!                 ***  SUBROUTINR lim_rhg  ***
      !!
      !! ** purpose :   determines the velocity field of sea ice by using
      !!  atmospheric (wind stress) and oceanic (water stress and surface
      !!  tilt) forcings. Ice-ice interaction is described by a non-linear
      !!  viscous-plastic law including shear strength and a bulk rheology.
      !!
      !! ** Action  : - compute u_ice, v_ice the sea-ice velocity
      !!
      !! History :
      !!   0.0  !  93-12  (M.A. Morales Maqueda.)  Original code
      !!   1.0  !  94-12  (H. Goosse) 
      !!   2.0  !  03-12  (C. Ethe, G. Madec)  F90, mpp
      !!-------------------------------------------------------------------
      ! * Arguments
      INTEGER, INTENT(in) ::   k_j1 ,  &  ! southern j-index for ice computation
         &                     k_jpj      ! northern j-index for ice computation

      ! * Local variables
      INTEGER ::   ji, jj              ! dummy loop indices

      INTEGER  :: &
         iim1, ijm1, iip1 , ijp1   , & ! temporary integers
         iter, jter                    !    "          "

      CHARACTER (len=50) :: charout

      REAL(wp) :: &
         ze11  , ze12  , ze22  , ze21  ,   &  ! temporary scalars
         zt11  , zt12  , zt21  , zt22  ,   &  !    "         "
         zvis11, zvis21, zvis12, zvis22,   &  !    "         "
         zgphsx, ztagnx, zusw  ,           &  !    "         "
         zgphsy, ztagny                       !    "         "
      REAL(wp) :: &
         zresm, zunw, zvnw, zur, zvr, zmod, za, zac, &
         zmpzas, zstms, zindu, zindu1, zusdtp, zmassdt, zcorlal,  &
         ztrace2, zdeter, zdelta, zsang, zmask, zdgp, zdgi, zdiag
      REAL(wp),DIMENSION(jpi,jpj) ::   &
         zpresh, zfrld, zmass, zcorl,     &
         zu0, zv0, zviszeta, zviseta,     &
         zc1u, zc1v, zc2u, zc2v, za1ct, za2ct, za1, za2, zb1, zb2,  &
         zc1, zc2, zd1, zd2, zden, zu_ice, zv_ice, zresr
      REAL(wp),DIMENSION(jpi,jpj,2,2) :: &
         zs11, zs12, zs22, zs21
      !!-------------------------------------------------------------------
      
      !  Store initial velocities
      !  ------------------------
      zu0(:,:) = u_ice(:,:)
      zv0(:,:) = v_ice(:,:)

      ! Ice mass, ice strength, and wind stress at the center            |
      ! of the grid squares.                                             |
      !-------------------------------------------------------------------

      DO jj = k_j1 , k_jpj-1
         DO ji = 1 , jpi
            za1(ji,jj)    = tms(ji,jj) * ( rhosn * hsnm(ji,jj) + rhoic * hicm(ji,jj) )
            zpresh(ji,jj) = tms(ji,jj) *  pstarh * hicm(ji,jj) * EXP( -c_rhg * frld(ji,jj) )
#if defined key_lim_cp1 && defined key_coupled
            zb1(ji,jj)    = tms(ji,jj) * gtaux(ji,jj) * ( 1.0 - frld(ji,jj) )
            zb2(ji,jj)    = tms(ji,jj) * gtauy(ji,jj) * ( 1.0 - frld(ji,jj) )
#else
            zb1(ji,jj)    = tms(ji,jj) * ( 1.0 - frld(ji,jj) )
            zb2(ji,jj)    = tms(ji,jj) * ( 1.0 - frld(ji,jj) )
#endif
         END DO
      END DO


      !---------------------------------------------------------------------------
      !  Wind stress, coriolis and mass terms at the corners of the grid squares |
      !  Gradient of ice strenght.                                               |
      !---------------------------------------------------------------------------
         
      DO jj = k_j1+1, k_jpj-1
         DO ji = 2, jpi
            zstms = tms(ji,jj  ) * wght(ji,jj,2,2) + tms(ji-1,jj  ) * wght(ji,jj,1,2)   &
               &  + tms(ji,jj-1) * wght(ji,jj,2,1) + tms(ji-1,jj-1) * wght(ji,jj,1,1)
            zusw  = 1.0 / MAX( zstms, epsd )

            zt11 = tms(ji  ,jj  ) * frld(ji  ,jj  ) 
            zt12 = tms(ji-1,jj  ) * frld(ji-1,jj  ) 
            zt21 = tms(ji  ,jj-1) * frld(ji  ,jj-1) 
            zt22 = tms(ji-1,jj-1) * frld(ji-1,jj-1)

            ! Leads area.
            zfrld(ji,jj) =  (  zt11 * wght(ji,jj,2,2) + zt12 * wght(ji,jj,1,2)   &
               &             + zt21 * wght(ji,jj,2,1) + zt22 * wght(ji,jj,1,1) ) * zusw

            ! Mass and coriolis coeff.
            zmass(ji,jj) = ( za1(ji,jj  ) * wght(ji,jj,2,2) + za1(ji-1,jj  ) * wght(ji,jj,1,2)   &
               &           + za1(ji,jj-1) * wght(ji,jj,2,1) + za1(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw
            zcorl(ji,jj) = zmass(ji,jj) * fcor(ji,jj)

            ! Wind stress.
#if defined key_lim_cp1 && defined key_coupled
            ztagnx = ( zb1(ji,jj  ) * wght(ji,jj,2,2) + zb1(ji-1,jj  ) * wght(ji,jj,1,2)   &
               &     + zb1(ji,jj-1) * wght(ji,jj,2,1) + zb1(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw
            ztagny = ( zb2(ji,jj  ) * wght(ji,jj,2,2) + zb2(ji-1,jj  ) * wght(ji,jj,1,2)   &
               &     + zb2(ji,jj-1) * wght(ji,jj,2,1) + zb2(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw
#else
            ztagnx = ( zb1(ji,jj  ) * wght(ji,jj,2,2) + zb1(ji-1,jj  ) * wght(ji,jj,1,2)   &
               &     + zb1(ji,jj-1) * wght(ji,jj,2,1) + zb1(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw * gtaux(ji,jj)
            ztagny = ( zb2(ji,jj  ) * wght(ji,jj,2,2) + zb2(ji-1,jj  ) * wght(ji,jj,1,2)   &
               &     + zb2(ji,jj-1) * wght(ji,jj,2,1) + zb2(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw * gtauy(ji,jj)
#endif

            ! Gradient of ice strength
            zgphsx =   ( alambd(ji,jj,2,2,2,1) - alambd(ji,jj,2,1,2,1) ) * zpresh(ji  ,jj-1)   &
               &     + ( alambd(ji,jj,2,2,2,2) - alambd(ji,jj,2,1,2,2) ) * zpresh(ji  ,jj  )   &
               &     - ( alambd(ji,jj,2,2,1,1) + alambd(ji,jj,2,1,1,1) ) * zpresh(ji-1,jj-1)   &
               &     - ( alambd(ji,jj,2,2,1,2) + alambd(ji,jj,2,1,1,2) ) * zpresh(ji-1,jj  )

            zgphsy = - ( alambd(ji,jj,1,1,2,1) + alambd(ji,jj,1,2,2,1) ) * zpresh(ji  ,jj-1)   &
               &     - ( alambd(ji,jj,1,1,1,1) + alambd(ji,jj,1,2,1,1) ) * zpresh(ji-1,jj-1)   &
               &     + ( alambd(ji,jj,1,1,2,2) - alambd(ji,jj,1,2,2,2) ) * zpresh(ji  ,jj  )   &
               &     + ( alambd(ji,jj,1,1,1,2) - alambd(ji,jj,1,2,1,2) ) * zpresh(ji-1,jj  )

            ! Computation of the velocity field taking into account the ice-ice interaction.                                 
            ! Terms that are independent of the velocity field.
            za1ct(ji,jj) = ztagnx - zcorl(ji,jj) * v_oce(ji,jj) - zgphsx
            za2ct(ji,jj) = ztagny + zcorl(ji,jj) * u_oce(ji,jj) - zgphsy
         END DO
      END DO

!! inutile!!
!!??    CALL lbc_lnk( za1ct, 'I', -1. )
!!??    CALL lbc_lnk( za2ct, 'I', -1. )


      ! SOLUTION OF THE MOMENTUM EQUATION.
      !------------------------------------------
      !                                                   ! ==================== !
      DO iter = 1 , 2 * nbiter                            !    loop over iter    !
         !                                                ! ==================== !        
         zindu = MOD( iter , 2 )
         zusdtp = ( zindu * 2.0 + ( 1.0 - zindu ) * 1.0 )  * REAL( nbiter ) / rdt_ice

         ! Computation of free drift field for free slip boundary conditions.

           DO jj = k_j1, k_jpj-1
              DO ji = 1, jpim1
                 !- Rate of strain tensor.
                 zt11 =   akappa(ji,jj,1,1) * ( u_ice(ji+1,jj) + u_ice(ji+1,jj+1) - u_ice(ji,jj  ) - u_ice(ji  ,jj+1) )  &
                    &   + akappa(ji,jj,1,2) * ( v_ice(ji+1,jj) + v_ice(ji+1,jj+1) + v_ice(ji,jj  ) + v_ice(ji  ,jj+1) )
                 zt12 = - akappa(ji,jj,2,2) * ( u_ice(ji  ,jj) + u_ice(ji+1,jj  ) - u_ice(ji,jj+1) - u_ice(ji+1,jj+1) )  &
                    &   - akappa(ji,jj,2,1) * ( v_ice(ji  ,jj) + v_ice(ji+1,jj  ) + v_ice(ji,jj+1) + v_ice(ji+1,jj+1) )
                 zt22 = - akappa(ji,jj,2,2) * ( v_ice(ji  ,jj) + v_ice(ji+1,jj  ) - v_ice(ji,jj+1) - v_ice(ji+1,jj+1) )  &
                    &   + akappa(ji,jj,2,1) * ( u_ice(ji  ,jj) + u_ice(ji+1,jj  ) + u_ice(ji,jj+1) + u_ice(ji+1,jj+1) )
                 zt21 =   akappa(ji,jj,1,1) * ( v_ice(ji+1,jj) + v_ice(ji+1,jj+1) - v_ice(ji,jj  ) - v_ice(ji  ,jj+1) )  &
                    &   - akappa(ji,jj,1,2) * ( u_ice(ji+1,jj) + u_ice(ji+1,jj+1) + u_ice(ji,jj  ) + u_ice(ji  ,jj+1) )

                 !- Rate of strain tensor. 
                 zdgp = zt11 + zt22
                 zdgi = zt12 + zt21
                 ztrace2 = zdgp * zdgp 
                 zdeter  = zt11 * zt22 - 0.25 * zdgi * zdgi

                 !  Creep limit depends on the size of the grid.
                 zdelta = MAX( SQRT( ztrace2 + ( ztrace2 - 4.0 * zdeter ) * usecc2),  creepl)

                 !-  Computation of viscosities.
                 zviszeta(ji,jj) = MAX( zpresh(ji,jj) / zdelta, etamn )
                 zviseta (ji,jj) = zviszeta(ji,jj) * usecc2
              END DO
           END DO
!!??       CALL lbc_lnk( zviszeta, 'I', -1. )  ! or T point???   semble reellement inutile
!!??       CALL lbc_lnk( zviseta , 'I', -1. )


           !-  Determination of zc1u, zc2u, zc1v and zc2v.
           DO jj = k_j1+1, k_jpj-1
              DO ji = 2, jpim1
                 ze11   =  akappa(ji-1,jj-1,1,1)
                 ze12   = +akappa(ji-1,jj-1,2,2)
                 ze22   =  akappa(ji-1,jj-1,2,1)
                 ze21   = -akappa(ji-1,jj-1,1,2)
                 zvis11 = 2.0 * zviseta (ji-1,jj-1) + dm
                 zvis22 =       zviszeta(ji-1,jj-1) - zviseta(ji-1,jj-1)
                 zvis12 =       zviseta (ji-1,jj-1) + dm
                 zvis21 =       zviseta (ji-1,jj-1)

                 zdiag = zvis22 * ( ze11 + ze22 )
                 zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
                 zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
                 zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag
                 zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12

                 ze11   = -akappa(ji,jj-1,1,1)
                 ze12   = +akappa(ji,jj-1,2,2)
                 ze22   =  akappa(ji,jj-1,2,1)
                 ze21   = -akappa(ji,jj-1,1,2)
                 zvis11 = 2.0 * zviseta (ji,jj-1) + dm
                 zvis22 =       zviszeta(ji,jj-1) - zviseta(ji,jj-1)
                 zvis12 =       zviseta (ji,jj-1) + dm
                 zvis21 =       zviseta (ji,jj-1)

                 zdiag = zvis22 * ( ze11 + ze22 )
                 zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
                 zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
                 zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
                 zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12

                 ze11   =  akappa(ji-1,jj,1,1)
                 ze12   = -akappa(ji-1,jj,2,2)
                 ze22   =  akappa(ji-1,jj,2,1)
                 ze21   = -akappa(ji-1,jj,1,2)
                 zvis11 = 2.0 * zviseta (ji-1,jj) + dm
                 zvis22 =       zviszeta(ji-1,jj) - zviseta(ji-1,jj)
                 zvis12 =       zviseta (ji-1,jj) + dm
                 zvis21 =       zviseta (ji-1,jj)

                 zdiag = zvis22 * ( ze11 + ze22 ) 
                 zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
                 zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
                 zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
                 zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12

                 ze11   = -akappa(ji,jj,1,1)
                 ze12   = -akappa(ji,jj,2,2)
                 ze22   =  akappa(ji,jj,2,1)
                 ze21   = -akappa(ji,jj,1,2)
                 zvis11 = 2.0 * zviseta (ji,jj) + dm
                 zvis22 =       zviszeta(ji,jj) - zviseta(ji,jj)
                 zvis12 =       zviseta (ji,jj) + dm
                 zvis21 =       zviseta (ji,jj)

                 zdiag = zvis22 * ( ze11 + ze22 )
                 zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
                 zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
                 zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag 
                 zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12
              END DO
           END DO

           DO jj = k_j1+1, k_jpj-1
              DO ji = 2, jpim1
                 zc1u(ji,jj) =   &
                    + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)   &
                    - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)   &
                    - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)   &
                    + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)   &
                    + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)   &
                    + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)   &
                    - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)   &
                    - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
                 
                 zc2u(ji,jj) =   &
                    + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)   &
                    - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)   &
                    - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)   &
                    + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)   &
                    - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)   &
                    - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)   &
                    + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)   &
                    + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
             END DO
           END DO

           DO jj = k_j1+1, k_jpj-1
              DO ji = 2, jpim1
                 !  zc1v , zc2v.
                 ze11   =  akappa(ji-1,jj-1,1,2)
                 ze12   = -akappa(ji-1,jj-1,2,1)
                 ze22   = +akappa(ji-1,jj-1,2,2)
                 ze21   =  akappa(ji-1,jj-1,1,1)
                 zvis11 = 2.0 * zviseta (ji-1,jj-1) + dm
                 zvis22 =       zviszeta(ji-1,jj-1) - zviseta(ji-1,jj-1)
                 zvis12 =       zviseta (ji-1,jj-1) + dm
                 zvis21 =       zviseta (ji-1,jj-1)

                 zdiag = zvis22 * ( ze11 + ze22 )
                 zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
                 zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
                 zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag
                 zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12
 
                 ze11   =  akappa(ji,jj-1,1,2)
                 ze12   = -akappa(ji,jj-1,2,1)
                 ze22   = +akappa(ji,jj-1,2,2)
                 ze21   = -akappa(ji,jj-1,1,1)
                 zvis11 = 2.0 * zviseta (ji,jj-1) + dm
                 zvis22 =       zviszeta(ji,jj-1) - zviseta(ji,jj-1)
                 zvis12 =       zviseta (ji,jj-1) + dm
                 zvis21 =       zviseta (ji,jj-1)

                 zdiag = zvis22 * ( ze11 + ze22 )
                 zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
                 zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
                 zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
                 zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12

                 ze11   =  akappa(ji-1,jj,1,2)
                 ze12   = -akappa(ji-1,jj,2,1)
                 ze22   = -akappa(ji-1,jj,2,2)
                 ze21   =  akappa(ji-1,jj,1,1)
                 zvis11 = 2.0 * zviseta (ji-1,jj) + dm
                 zvis22 =       zviszeta(ji-1,jj) - zviseta(ji-1,jj)
                 zvis12 =       zviseta (ji-1,jj) + dm
                 zvis21 =       zviseta (ji-1,jj)

                 zdiag = zvis22 * ( ze11 + ze22 )
                 zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
                 zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
                 zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
                 zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12

                 ze11   =  akappa(ji,jj,1,2)
                 ze12   = -akappa(ji,jj,2,1)
                 ze22   = -akappa(ji,jj,2,2)
                 ze21   = -akappa(ji,jj,1,1)
                 zvis11 = 2.0 * zviseta (ji,jj) + dm
                 zvis22 =       zviszeta(ji,jj) - zviseta(ji,jj)
                 zvis12 =       zviseta (ji,jj) + dm
                 zvis21 =       zviseta (ji,jj)

                 zdiag = zvis22 * ( ze11 + ze22 )
                 zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
                 zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
                 zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag
                 zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12

              END DO
           END DO

           DO jj = k_j1+1, k_jpj-1
              DO ji = 2, jpim1
                 zc1v(ji,jj) =   &
                    + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)   &
                    - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)   &
                    - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)   &
                    + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)   &
                    + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)   &
                    + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)   &
                    - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)   &
                    - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
                 zc2v(ji,jj) =   &
                    + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)   &
                    - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)   &
                    - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)   &
                    + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)   &
                    - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)   &
                    - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)   &
                    + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)   &
                    + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
              END DO
           END DO

         ! Relaxation.
           
iflag:   DO jter = 1 , nbitdr

            !  Store previous drift field.   
            DO jj = k_j1, k_jpj-1
               zu_ice(:,jj) = u_ice(:,jj)
               zv_ice(:,jj) = v_ice(:,jj)
            END DO

            DO jj = k_j1+1, k_jpj-1
               zsang  = SIGN( 1.e0, fcor(1,jj) ) * sangvg   ! only the sinus changes its sign with the hemisphere
               DO ji = 2, jpim1
                 zur     = u_ice(ji,jj) - u_oce(ji,jj)
                 zvr     = v_ice(ji,jj) - v_oce(ji,jj)
                 zmod    = SQRT( zur * zur + zvr * zvr) * ( 1.0 - zfrld(ji,jj) )
                 za      = rhoco * zmod
                 zac     = za * cangvg
                  zmpzas  = alpha * zcorl(ji,jj) + za * zsang
                  zmassdt = zusdtp * zmass(ji,jj)
                  zcorlal = ( 1.0 - alpha ) * zcorl(ji,jj)

                  za1(ji,jj) =  zmassdt * zu0(ji,jj) + zcorlal * zv0(ji,jj) + za1ct(ji,jj)   &
                     &        + za * ( cangvg * u_oce(ji,jj) - zsang * v_oce(ji,jj) )

                  za2(ji,jj) =  zmassdt * zv0(ji,jj) - zcorlal * zu0(ji,jj) + za2ct(ji,jj)   &
                     &        + za * ( cangvg * v_oce(ji,jj) + zsang * u_oce(ji,jj) )

                  zb1(ji,jj)  = zmassdt + zac - zc1u(ji,jj)
                  zb2(ji,jj)  = zmpzas  - zc2u(ji,jj)
                  zc1(ji,jj)  = zmpzas  + zc1v(ji,jj)
                  zc2(ji,jj)  = zmassdt + zac  - zc2v(ji,jj) 
                  zdeter      = zc1(ji,jj) * zb2(ji,jj) + zc2(ji,jj) * zb1(ji,jj)
                  zden(ji,jj) = SIGN( rone, zdeter) / MAX( epsd , ABS( zdeter ) )
               END DO
            END DO

            ! The computation of ice interaction term is splitted into two parts
            !-------------------------------------------------------------------------

            ! Terms that do not involve already up-dated velocities.
         
            DO jj = k_j1+1, k_jpj-1
               DO ji = 2, jpim1
                  iim1 = ji
                  ijm1 = jj - 1
                  iip1 = ji + 1
                  ijp1 = jj
                  ze11 =   akappa(iim1,ijm1,1,1) * u_ice(iip1,ijp1) + akappa(iim1,ijm1,1,2) * v_ice(iip1,ijp1)
                  ze12 = + akappa(iim1,ijm1,2,2) * u_ice(iip1,ijp1) - akappa(iim1,ijm1,2,1) * v_ice(iip1,ijp1)
                  ze22 = + akappa(iim1,ijm1,2,2) * v_ice(iip1,ijp1) + akappa(iim1,ijm1,2,1) * u_ice(iip1,ijp1)
                  ze21 =   akappa(iim1,ijm1,1,1) * v_ice(iip1,ijp1) - akappa(iim1,ijm1,1,2) * u_ice(iip1,ijp1)
                  zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
                  zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
                  zvis12 =       zviseta (iim1,ijm1) + dm
                  zvis21 =       zviseta (iim1,ijm1)
                  zdiag = zvis22 * ( ze11 + ze22 )
                  zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
                  zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
                  zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
                  zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12


                  iim1 = ji - 1
                  ijm1 = jj
                  iip1 = ji
                  ijp1 = jj + 1                   
                  ze11 =   akappa(iim1,ijm1,1,1) * ( u_ice(iip1,ijp1) - u_ice(iim1,ijp1) )   &
                     &   + akappa(iim1,ijm1,1,2) * ( v_ice(iip1,ijp1) + v_ice(iim1,ijp1) )
                  ze12 = + akappa(iim1,ijm1,2,2) * ( u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )   &
                     &   - akappa(iim1,ijm1,2,1) * ( v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )
                  ze22 = + akappa(iim1,ijm1,2,2) * ( v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )   &
                     &   + akappa(iim1,ijm1,2,1) * ( u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )
                  ze21 =   akappa(iim1,ijm1,1,1) * ( v_ice(iip1,ijp1) - v_ice(iim1,ijp1) )   &
                     &   - akappa(iim1,ijm1,1,2) * ( u_ice(iip1,ijp1) + u_ice(iim1,ijp1) )
                  zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
                  zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
                  zvis12 =       zviseta (iim1,ijm1) + dm
                  zvis21 =       zviseta (iim1,ijm1)
                  zdiag = zvis22 * ( ze11 + ze22 )
                  zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
                  zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
                  zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
                  zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12

                  iim1 = ji
                  ijm1 = jj
                  iip1 = ji + 1
                  ijp1 = jj + 1
                  ze11 =   akappa(iim1,ijm1,1,1) * ( u_ice(iip1,ijm1) + u_ice(iip1,ijp1) - u_ice(iim1,ijp1) )   &
                     &   + akappa(iim1,ijm1,1,2) * ( v_ice(iip1,ijm1) + v_ice(iip1,ijp1) + v_ice(iim1,ijp1) )
                  ze12 = - akappa(iim1,ijm1,2,2) * ( u_ice(iip1,ijm1) - u_ice(iim1,ijp1) - u_ice(iip1,ijp1) )   &
                     &   - akappa(iim1,ijm1,2,1) * ( v_ice(iip1,ijm1) + v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )
                  ze22 = - akappa(iim1,ijm1,2,2) * ( v_ice(iip1,ijm1) - v_ice(iim1,ijp1) - v_ice(iip1,ijp1) )   &
                     &   + akappa(iim1,ijm1,2,1) * ( u_ice(iip1,ijm1) + u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )
                  ze21 =   akappa(iim1,ijm1,1,1) * ( v_ice(iip1,ijm1) + v_ice(iip1,ijp1) - v_ice(iim1,ijp1) )   &
                     &   - akappa(iim1,ijm1,1,2) * ( u_ice(iip1,ijm1) + u_ice(iip1,ijp1) + u_ice(iim1,ijp1) ) 
                  zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
                  zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
                  zvis12 =       zviseta (iim1,ijm1) + dm
                  zvis21 =       zviseta (iim1,ijm1)

                  zdiag = zvis22 * ( ze11 + ze22 )
                  zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
                  zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
                  zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag
                  zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12

               END DO
            END DO

            ! Terms involving already up-dated velocities.
            !-Using the arrays zu_ice and zv_ice in the computation of the terms ze leads to JACOBI's method; 
            ! Using arrays u and v in the computation of the terms ze leads to GAUSS-SEIDEL method.
             
            DO jj = k_j1+1, k_jpj-1
               DO ji = 2, jpim1
                  iim1 = ji - 1
                  ijm1 = jj - 1
                  iip1 = ji
                  ijp1 = jj
                  ze11 =   akappa(iim1,ijm1,1,1) * ( zu_ice(iip1,ijm1) - zu_ice(iim1,ijm1) - zu_ice(iim1,ijp1) )   &
                     &   + akappa(iim1,ijm1,1,2) * ( zv_ice(iip1,ijm1) + zv_ice(iim1,ijm1) + zv_ice(iim1,ijp1) )
                  ze12 = - akappa(iim1,ijm1,2,2) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) - zu_ice(iim1,ijp1) )   &
                     &   - akappa(iim1,ijm1,2,1) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) + zv_ice(iim1,ijp1) )
                  ze22 = - akappa(iim1,ijm1,2,2) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) - zv_ice(iim1,ijp1) )   &
                     &   + akappa(iim1,ijm1,2,1) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) + zu_ice(iim1,ijp1) )
                  ze21 =   akappa(iim1,ijm1,1,1) * ( zv_ice(iip1,ijm1) - zv_ice(iim1,ijm1) - zv_ice(iim1,ijp1) )   &
                     &  -  akappa(iim1,ijm1,1,2) * ( zu_ice(iip1,ijm1) + zu_ice(iim1,ijm1) + zu_ice(iim1,ijp1) )
                  zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
                  zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
                  zvis12 =       zviseta (iim1,ijm1) + dm
                  zvis21 =       zviseta (iim1,ijm1)

                  zdiag = zvis22 * ( ze11 + ze22 )
                  zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
                  zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
                  zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag
                  zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12

#if defined key_agrif
             END DO
          END DO

          DO jj = k_j1+1, k_jpj-1
             DO ji = 2, jpim1
#endif

                  iim1 = ji
                  ijm1 = jj - 1
                  iip1 = ji + 1
                  ze11 =   akappa(iim1,ijm1,1,1) * ( zu_ice(iip1,ijm1) - zu_ice(iim1,ijm1) )   &
                     &   + akappa(iim1,ijm1,1,2) * ( zv_ice(iip1,ijm1) + zv_ice(iim1,ijm1) )
                  ze12 = - akappa(iim1,ijm1,2,2) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) )   &
                     &   - akappa(iim1,ijm1,2,1) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) )
                  ze22 = - akappa(iim1,ijm1,2,2) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) )   &
                     &   + akappa(iim1,ijm1,2,1) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) )
                  ze21 =   akappa(iim1,ijm1,1,1) * ( zv_ice(iip1,ijm1) - zv_ice(iim1,ijm1) )   &
                     &   - akappa(iim1,ijm1,1,2) * ( zu_ice(iip1,ijm1) + zu_ice(iim1,ijm1) )
                  zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
                  zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
                  zvis12 =       zviseta (iim1,ijm1) + dm
                  zvis21 =       zviseta (iim1,ijm1)

                  zdiag = zvis22 * ( ze11 + ze22 )
                  zs11(ji,jj,2,1) =  zs11(ji,jj,2,1) + zvis11 * ze11 + zdiag
                  zs12(ji,jj,2,1) =  zs12(ji,jj,2,1) + zvis12 * ze12 + zvis21 * ze21
                  zs22(ji,jj,2,1) =  zs22(ji,jj,2,1) + zvis11 * ze22 + zdiag
                  zs21(ji,jj,2,1) =  zs21(ji,jj,2,1) + zvis12 * ze21 + zvis21 * ze12


                  iim1 = ji - 1
                  ijm1 = jj 
                  ze11 = - akappa(iim1,ijm1,1,1) * zu_ice(iim1,ijm1) + akappa(iim1,ijm1,1,2) * zv_ice(iim1,ijm1)
                  ze12 = - akappa(iim1,ijm1,2,2) * zu_ice(iim1,ijm1) - akappa(iim1,ijm1,2,1) * zv_ice(iim1,ijm1)
                  ze22 = - akappa(iim1,ijm1,2,2) * zv_ice(iim1,ijm1) + akappa(iim1,ijm1,2,1) * zu_ice(iim1,ijm1)
                  ze21 = - akappa(iim1,ijm1,1,1) * zv_ice(iim1,ijm1) - akappa(iim1,ijm1,1,2) * zu_ice(iim1,ijm1)
                  zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
                  zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
                  zvis12 =       zviseta (iim1,ijm1) + dm
                  zvis21 =       zviseta (iim1,ijm1)

                  zdiag = zvis22 * ( ze11 + ze22 )
                  zs11(ji,jj,1,2) =  zs11(ji,jj,1,2) + zvis11 * ze11 + zdiag 
                  zs12(ji,jj,1,2) =  zs12(ji,jj,1,2) + zvis12 * ze12 + zvis21 * ze21
                  zs22(ji,jj,1,2) =  zs22(ji,jj,1,2) + zvis11 * ze22 + zdiag
                  zs21(ji,jj,1,2) =  zs21(ji,jj,1,2) + zvis12 * ze21 + zvis21 * ze12

#if defined key_agrif
             END DO
          END DO

          DO jj = k_j1+1, k_jpj-1
             DO ji = 2, jpim1
#endif
                  zd1(ji,jj) =   &
                     + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)  &
                     - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)  &
                     - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)  &
                     + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)  &
                     + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)  &
                     + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)  &
                     - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)  &
                     - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
                  zd2(ji,jj) =   &
                     + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)  &
                     - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)  &
                     - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)  &
                     + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)  &
                     - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)  &
                     - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)  &
                     + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)  &
                     + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
               END DO
            END DO

            DO jj = k_j1+1, k_jpj-1
               DO ji = 2, jpim1
                  zunw = (  ( za1(ji,jj) + zd1(ji,jj) ) * zc2(ji,jj)        &
                     &    + ( za2(ji,jj) + zd2(ji,jj) ) * zc1(ji,jj) ) * zden(ji,jj)

                  zvnw = (  ( za2(ji,jj) + zd2(ji,jj) ) * zb1(ji,jj)        &
                     &    - ( za1(ji,jj) + zd1(ji,jj) ) * zb2(ji,jj) ) * zden(ji,jj)

                  zmask = ( 1.0 - MAX( rzero, SIGN( rone , 1.0 - zmass(ji,jj) ) ) ) * tmu(ji,jj)

                  u_ice(ji,jj) = ( u_ice(ji,jj) + om * ( zunw - u_ice(ji,jj) ) * tmu(ji,jj) ) * zmask
                  v_ice(ji,jj) = ( v_ice(ji,jj) + om * ( zvnw - v_ice(ji,jj) ) * tmu(ji,jj) ) * zmask
               END DO
            END DO

            CALL lbc_lnk( u_ice, 'I', -1. )
            CALL lbc_lnk( v_ice, 'I', -1. )

            !---  5.2.5.4. Convergence test.
            DO jj = k_j1+1 , k_jpj-1
               zresr(:,jj) = MAX( ABS( u_ice(:,jj) - zu_ice(:,jj) ) , ABS( v_ice(:,jj) - zv_ice(:,jj) ) )
            END DO
            zresm = MAXVAL( zresr( 1:jpi , k_j1+1:k_jpj-1 ) )
            IF( lk_mpp )   CALL mpp_max( zresm )   ! max over the global domain

            IF ( zresm <= resl) EXIT iflag

         END DO iflag

         zindu1 = 1.0 - zindu
         DO jj = k_j1 , k_jpj-1
            zu0(:,jj) = zindu * zu0(:,jj) + zindu1 * u_ice(:,jj)
            zv0(:,jj) = zindu * zv0(:,jj) + zindu1 * v_ice(:,jj)
         END DO
      !                                                   ! ==================== !
      END DO                                              !  end loop over iter  !
      !                                                   ! ==================== !

      IF(ln_ctl) THEN
         WRITE(charout,FMT="('lim_rhg  : res =',D23.16, ' iter =',I4)") zresm, jter
         CALL prt_ctl_info(charout)
         CALL prt_ctl(tab2d_1=u_ice, clinfo1=' lim_rhg  : u_ice :', tab2d_2=v_ice, clinfo2=' v_ice :')
      ENDIF

   END SUBROUTINE lim_rhg

#else
   !!----------------------------------------------------------------------
   !!   Default option          Dummy module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE lim_rhg( k1 , k2 )         ! Dummy routine
      if(lwp) WRITE(numout,*) 'lim_rhg: You should not have seen this print! error?', k1, k2
   END SUBROUTINE lim_rhg
#endif

   !!==============================================================================
END MODULE limrhg
