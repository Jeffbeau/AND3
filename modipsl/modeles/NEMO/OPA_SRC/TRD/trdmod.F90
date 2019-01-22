MODULE trdmod
   !!======================================================================
   !!                       ***  MODULE  trdmod  ***
   !! Ocean diagnostics:  ocean tracers and dynamic trends
   !!=====================================================================
#if  defined key_trdtra || defined key_trddyn || defined key_trdmld || defined key_trdvor || defined key_esopa
   !!----------------------------------------------------------------------
   !!   trd_mod          : Call the trend to be computed
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce                     ! ocean dynamics and tracers variables
   USE dom_oce                 ! ocean space and time domain variables
   USE trdmod_oce              ! ocean variables trends
   USE trdvor                  ! ocean vorticity trends 
   USE trdicp                  ! ocean bassin integral constraints properties
   USE trdmld                  ! ocean active mixed layer tracers trends 
   USE trabbl                  ! bottom boundary layer variables
   USE in_out_manager          ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trd_mod        ! called by all dynXX or traXX modules

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmod.F90,v 1.2 2005/03/27 18:35:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trd_mod(ptrdx, ptrdy, ktrd, ctype, kt)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mod  ***
      !! 
      !! ** Purpose : Dispatch all trends computation, e.g. vorticity, mld or 
      !!              integral constrains
      !!
      !! ** Method :
      !!
      !! History :
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used
#if defined key_trabbl_adv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &  ! temporary arrays
         &         zun, zvn
#else
      USE oce                , zun => un,  &  ! When no bbl, zun == un
         &                     zvn => vn      ! When no bbl, zvn == vn
#endif

      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
         ptrdx,                      &   ! Temperature or U trend 
         ptrdy                           ! Salinity    or V trend

      INTEGER, INTENT( in ) ::   &
         kt  ,                   & ! time step
         ktrd                      ! tracer trend index

      CHARACTER(len=3), INTENT( in ) ::   &
         ctype                             ! momentum or tracers trends type
         !                                 ! 'DYN' or 'TRA'

      !! * Local save
      REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
         zbtr2

      !! * Local declarations
      INTEGER ::   ji, jj, jk    ! loop indices
      REAL(wp) ::   &
         zbtr,            &  ! temporary scalars
         zfui, zfvj,           &  !    "         "
         zfui1, zfvj1             !    "         "
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         z2dx, z2dy                        ! workspace arrays
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         z3dx, z3dy                            ! workspace arrays
      !!----------------------------------------------------------------------

      ! Initialization of workspace arrays
      z3dx(:,:,:) = 0.e0
      z3dy(:,:,:) = 0.e0
      z2dx(:,:) = 0.e0
      z2dy(:,:) = 0.e0

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! I. Bassin averaged properties for momentum and/or tracers trends
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF( ( mod(kt,ntrd) == 0 .OR. kt == nit000 .OR. kt == nitend) )   THEN

         ! Active tracers trends 
         IF( lk_trdtra .AND. ctype == 'TRA' )   THEN

            IF( ktrd == jpttdnsr )   THEN
               ! 2D array tracers surface forcing
               z2dx(:,:) = ptrdx(:,:,1)
               z2dy(:,:) = ptrdy(:,:,1)

               CALL trd(z2dx, z2dy, ktrd, ctype)
            ELSE
               ! 3D array
               CALL trd(ptrdx, ptrdy, ktrd, ctype)
            ENDIF

         ENDIF

         ! Momentum trends 
         IF( lk_trddyn .AND. ctype == 'DYN' )   THEN

            IF( ktrd == jpdtdswf .OR. ktrd == jpdtdbfr )   THEN
               ! momentum surface forcing/bottom friction  2D array
               z2dx(:,:) = ptrdx(:,:,1)
               z2dy(:,:) = ptrdy(:,:,1)

               CALL trd(z2dx, z2dy, ktrd, ctype)
            ELSE
               ! 3D array
               CALL trd(ptrdx, ptrdy, ktrd, ctype)
            ENDIF

         ENDIF

      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! II. Vorticity trends
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF( lk_trdvor .AND. ctype == 'DYN' )   THEN

         SELECT CASE ( ktrd )

         ! Pressure Gradient trend
         CASE ( jpdtdhpg )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorprg)

         ! KE Gradient trend
         CASE ( jpdtdkeg )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorkeg)

         ! Relative Vorticity trend
         CASE ( jpdtdrvo )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorrvo)

         ! Planetary Vorticity Term trend
         CASE ( jpdtdpvo )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorpvo)

         ! Horizontal Diffusion trend
         CASE ( jpdtdldf )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorldf)

         ! Vertical Advection trend
         CASE ( jpdtdzad )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorzad)

         ! Vertical Diffusion trend
         CASE ( jpdtdzdf )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorzdf)

         ! Surface Pressure Grad. trend
         CASE ( jpdtdspg )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorspg)

         ! Beta V trend 
         CASE ( jpdtddat )      
            CALL trd_vor_zint(ptrdx, ptrdy, jpvorbev)

         ! Wind stress forcing term
         CASE ( jpdtdswf )      
            z2dx(:,:) = ptrdx(:,:,1)
            z2dy(:,:) = ptrdy(:,:,1)

            CALL trd_vor_zint(z2dx, z2dy, jpvorswf)

         ! Bottom friction term
         CASE ( jpdtdbfr )      
            z2dx(:,:) = ptrdx(:,:,1)
            z2dy(:,:) = ptrdy(:,:,1)

            CALL trd_vor_zint(z2dx, z2dy, jpvorbfr)

         END SELECT

      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! III. Mixed layer trends
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF( lk_trdmld .AND. ctype == 'TRA' )   THEN
         
         SELECT CASE ( ktrd )

         ! horizontal advection trends
         CASE ( jpttdlad )      

#if defined key_trabbl_adv
            ! Advective bottom boundary layer 
            ! -------------------------------
            zun(:,:,:) = un(:,:,:) - u_bbl(:,:,:)
            zvn(:,:,:) = vn(:,:,:) - v_bbl(:,:,:)
#endif
            IF( kt == nit000 )   zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )

            SELECT CASE ( l_adv )

            CASE ( 'ce2' )

               ! Split horizontal trends into i- and j- compnents for trdmld case 
               ! ----------------------------------------------------------------

               ! i- advective trend computed as Uh gradh(T)
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1   ! vector opt.
# if defined key_s_coord || defined key_partial_steps
                        zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)

                        zfui = 0.5 * e2u(ji  ,jj) * fse3u(ji,  jj,jk) * zun(ji,  jj,jk)
                        zfui1= 0.5 * e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * zun(ji-1,jj,jk)
# else         
                        zbtr = zbtr2(ji,jj)

                        zfui = 0.5 * e2u(ji  ,jj) * zun(ji,  jj,jk)
                        zfui1= 0.5 * e2u(ji-1,jj) * zun(ji-1,jj,jk)
# endif
                        ! save i- advective trend 
                        z3dx(ji,jj,jk) = - zbtr * ( zfui  * ( tn(ji+1,jj,jk) - tn(ji  ,jj,jk) )    &
                            &                     + zfui1 * ( tn(ji  ,jj,jk) - tn(ji-1,jj,jk) ) )
                        z3dy(ji,jj,jk) = - zbtr * ( zfui  * ( sn(ji+1,jj,jk) - sn(ji  ,jj,jk) )    &
                            &                     + zfui1 * ( sn(ji  ,jj,jk) - sn(ji-1,jj,jk) ) )
                     END DO
                  END DO
               END DO

               ! save the i- horizontal trends for diagnostic
               CALL trd_mld_zint(z3dx, z3dy, jpmldxad, '3D')

               ! j- advective trend computed as Uh gradh(T)
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1   ! vector opt.
# if defined key_s_coord || defined key_partial_steps
                        zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)

                        zfvj = 0.5 * e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * zvn(ji,jj  ,jk)
                        zfvj1= 0.5 * e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * zvn(ji,jj-1,jk)
# else         
                        zbtr = zbtr2(ji,jj)

                        zfvj = 0.5 * e1v(ji,jj  ) * zvn(ji,jj  ,jk)
                        zfvj1= 0.5 * e1v(ji,jj-1) * zvn(ji,jj-1,jk)
# endif
                        ! save j- advective trend 
                        z3dx(ji,jj,jk) = - zbtr * ( zfvj  * ( tn(ji,jj+1,jk) - tn(ji,jj  ,jk) )   &
                            &                     + zfvj1 * ( tn(ji,jj  ,jk) - tn(ji,jj-1,jk) ) )
                        z3dy(ji,jj,jk) = - zbtr * ( zfvj  * ( sn(ji,jj+1,jk) - sn(ji,jj  ,jk) )   &
                            &                     + zfvj1 * ( sn(ji,jj  ,jk) - sn(ji,jj-1,jk) ) )
                     END DO
                  END DO
               END DO

               ! save the j- horizontal trend for diagnostic
               CALL trd_mld_zint(z3dx, z3dy, jpmldyad, '3D')

            CASE ( 'tvd' )

               ! Recompute the horizontal advection term Div(Uh.T) term 
               z3dx(:,:,:) = ptrdx(:,:,:) - tn(:,:,:) * hdivn(:,:,:)
               z3dy(:,:,:) = ptrdy(:,:,:) - sn(:,:,:) * hdivn(:,:,:)

               ! Deduce the i- horizontal advection in substracting the j- one.
               ! tladj()/sladj() are computed in traadv_tvd.F90 module
               z3dx(:,:,:) = z3dx(:,:,:) - tladj(:,:,:)
               z3dy(:,:,:) = z3dy(:,:,:) - sladj(:,:,:)

               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1
                        zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)

                        ! Compute the zonal et meridional divergence
                        zfui = e2u(ji  ,jj) * fse3u(ji  ,jj,jk) * zun(ji  ,jj,jk)  &
                             - e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * zun(ji-1,jj,jk)
                        zfvj = e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * zvn(ji,jj  ,jk)  &
                             - e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * zvn(ji,jj-1,jk)

                        ! i- advective trend computed as U gradx(T/S)
                        z3dx(ji,jj,jk) = z3dx(ji,jj,jk) + tn(ji,jj,jk) * zfui * zbtr
                        z3dy(ji,jj,jk) = z3dy(ji,jj,jk) + sn(ji,jj,jk) * zfui * zbtr

                        ! j- advective trend computed as V grady(T/S)
                        tladj(ji,jj,jk) = tladj(ji,jj,jk) + tn(ji,jj,jk) * zfvj * zbtr
                        sladj(ji,jj,jk) = sladj(ji,jj,jk) + sn(ji,jj,jk) * zfvj * zbtr

                     END DO
                  END DO
               END DO

               ! save the i- horizontal trend for diagnostic
               CALL trd_mld_zint(z3dx, z3dy, jpmldxad, '3D')

               ! save the j- horizontal trend for diagnostic
               CALL trd_mld_zint(tladj, sladi, jpmldyad, '3D')

            CASE ( 'mus', 'mu2' )

               !  Split horizontal trends in i- and j- direction for trdmld case 
               ! ----------------------------------------------------------------

               ! i- advective trend computed as U gradx(T/S)
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1      
                     DO ji = fs_2, fs_jpim1   ! vector opt.
# if defined key_s_coord || defined key_partial_steps
                        zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
                        zfui =  e2u(ji  ,jj) * fse3u(ji,  jj,jk) * zun(ji,  jj,jk)   &
                           & -  e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * zun(ji-1,jj,jk)
# else      
                        zbtr = zbtr2(ji,jj)
                        zfui = e2u(ji  ,jj) * zun(ji,  jj,jk)   &
                           & - e2u(ji-1,jj) * zun(ji-1,jj,jk)
# endif
                        ! save i- advective trend 
                        z3dx(ji,jj,jk) = - zbtr * ( tladi(ji,jj,jk) - tladi(ji-1,jj,jk) )   &
                            &                      + tn(ji,jj,jk) * zfui * zbtr
                        z3dy(ji,jj,jk) = - zbtr * ( sladi(ji,jj,jk) - sladi(ji-1,jj,jk) )  &
                            &                      + sn(ji,jj,jk) * zfui * zbtr
                     END DO
                  END DO
               END DO        

               ! save the i- horizontal trends for diagnostic
               CALL trd_mld_zint(z3dx, z3dy, jpmldxad, '3D')

               ! j- advective trend computed as V grady(T/S)
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1      
                     DO ji = fs_2, fs_jpim1   ! vector opt.
# if defined key_s_coord || defined key_partial_steps
                        zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
                        zfvj =  e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * zvn(ji,jj  ,jk)   &
                           & -  e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * zvn(ji,jj-1,jk)
# else      
                        zbtr = zbtr2(ji,jj)
                        zfvj = e1v(ji,jj  ) * zvn(ji,jj  ,jk)   &
                           & - e1v(ji,jj-1) * zvn(ji,jj-1,jk)
# endif
                        ! save j- advective trend 
                        z3dx(ji,jj,jk) =  - zbtr * ( tladj(ji,jj,jk) - tladj(ji,jj-1,jk) )   &
                            &                       + tn(ji,jj,jk) * zfvj * zbtr
                        z3dy(ji,jj,jk) =  - zbtr * ( sladj(ji,jj,jk) - sladj(ji,jj-1,jk) )   &
                            &                       + sn(ji,jj,jk) * zfvj * zbtr
                     END DO
                  END DO
               END DO        

               ! save the j- horizontal trends for diagnostic
               CALL trd_mld_zint(z3dx, z3dy, jpmldyad, '3D')

            END SELECT

         ! vertical advection trends
         CASE ( jpttdzad )      
            CALL trd_mld_zint(ptrdx, ptrdy, jpmldzad, '3D')

         ! lateral diffusion trends
         CASE ( jpttdldf )      
            CALL trd_mld_zint(ptrdx, ptrdy, jpmldldf, '3D')
# if defined key_traldf_eiv
            ! Save the i- and j- eddy induce velocity trends
            CALL trd_mld_zint(tladi, sladi, jpmldxei, '3D')
            CALL trd_mld_zint(tladj, sladj, jpmldyei, '3D')
# endif
            IF( lk_trabbl_dif )   THEN
               z3dx(:,:,:) = 0.e0
               z3dy(:,:,:) = 0.e0
               z3dx(:,:,1) = tldfbbl(:,:)
               z3dy(:,:,1) = sldfbbl(:,:)
               CALL trd_mld_zint(z3dx, z3dy, jpmldldf, '2D')
            ENDIF

         ! vertical diffusion trends
         CASE ( jpttdzdf )      
            CALL trd_mld_zint(ptrdx, ptrdy, jpmldzdf, '3D')

         ! vertical diffusion trends
         CASE ( jpttddoe )      
            CALL trd_mld_zint(ptrdx, ptrdy, jpmldzei, '3D')

         ! penetrative solar radiation trends
         CASE ( jpttdqsr )      
            CALL trd_mld_zint(ptrdx, ptrdy, jpmldfor, '3D')

         ! non penetrative solar radiation trends
         CASE ( jpttdnsr )
            ptrdx(:,:,2:jpk) = 0.e0
            ptrdy(:,:,2:jpk) = 0.e0
            CALL trd_mld_zint(ptrdx, ptrdy, jpmldfor, '2D')

         END SELECT   

      ENDIF


   END SUBROUTINE trd_mod

#   else
   !!----------------------------------------------------------------------
   !!   Default case :                                         Empty module
   !!----------------------------------------------------------------------
   USE trdmod_oce      ! ocean variables trends

CONTAINS
   SUBROUTINE trd_mod(ptrd3dx, ptrd3dy, ktrd , ctype, kt)       ! Empty routine
      REAL, DIMENSION(:,:,:), INTENT( in ) ::   &
          ptrd3dx,                     &   ! Temperature or U trend 
          ptrd3dy                          ! Salinity    or V trend
      INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
      INTEGER, INTENT( in ) ::   kt        ! Time step
      CHARACTER(len=3), INTENT( in ) ::   &
         ctype                             ! momentum or tracers trends type
!      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd3dx(1,1,1)
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd3dy(1,1,1)
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', kt
   END SUBROUTINE trd_mod
#   endif

   !!======================================================================
END MODULE trdmod
