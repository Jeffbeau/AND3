MODULE limdyn
   !!======================================================================
   !!                     ***  MODULE  limdyn  ***
   !!   Sea-Ice dynamics :  
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !!    lim_dyn      : computes ice velocities
   !!    lim_dyn_init : initialization and namelist read
   !!----------------------------------------------------------------------
   !! * Modules used
   USE phycst
   USE in_out_manager  ! I/O manager
   USE dom_ice
   USE dom_oce         ! ocean space and time domain
   USE ice
   USE ice_oce
   USE iceini
   USE limistate
   USE limrhg          ! ice rheology
   USE lbclnk
   USE lib_mpp
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC lim_dyn  ! routine called by ice_step

   !! * Module variables
   REAL(wp)  ::  rone    = 1.e0   ! constant value

   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdyn.F90,v 1.7 2005/09/22 10:43:39 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_dyn
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE lim_dyn  ***
      !!               
      !! ** Purpose :   compute ice velocity and ocean-ice stress
      !!                
      !! ** Method  : 
      !!
      !! ** Action  : - Initialisation
      !!              - Call of the dynamic routine for each hemisphere
      !!              - computation of the stress at the ocean surface         
      !!              - treatment of the case if no ice dynamic
      !! History :
      !!   1.0  !  01-04  (LIM)  Original code
      !!   2.0  !  02-08  (C. Ethe, G. Madec)  F90, mpp
      !!---------------------------------------------------------------------
      !! * Loal variables
      INTEGER ::   ji, jj             ! dummy loop indices
      INTEGER ::   i_j1, i_jpj        ! Starting/ending j-indices for rheology
      REAL(wp) ::   &
         ztairx, ztairy,           &  ! tempory scalars
         zsang , zmod,             &
         ztglx , ztgly ,           &
         zt11, zt12, zt21, zt22 ,  &
         zustm, zsfrld, zsfrldm4,  &
         zu_ice, zv_ice, ztair2
      REAL(wp),DIMENSION(jpj) ::   &
         zind,                     &  ! i-averaged indicator of sea-ice
         zmsk                         ! i-averaged of tmask
      !!---------------------------------------------------------------------

      IF( numit == nstart  )   CALL lim_dyn_init   ! Initialization (first time-step only)
      
      IF ( ln_limdyn ) THEN

         ! Mean ice and snow thicknesses.          
         hsnm(:,:)  = ( 1.0 - frld(:,:) ) * hsnif(:,:)
         hicm(:,:)  = ( 1.0 - frld(:,:) ) * hicif(:,:)

         u_oce(:,:)  = u_io(:,:) * tmu(:,:)
         v_oce(:,:)  = v_io(:,:) * tmu(:,:)
       
         !                                         ! Rheology (ice dynamics)
         !                                         ! ========
         
         !  Define the j-limits where ice rheology is computed
         ! ---------------------------------------------------
         
         IF( lk_mpp ) THEN                    ! mpp: compute over the whole domain
            i_j1 = 1   
            i_jpj = jpj
            IF(ln_ctl)    THEN
               CALL prt_ctl_info('lim_dyn  :    i_j1 = ', ivar1=i_j1, clinfo2=' ij_jpj = ', ivar2=i_jpj)
            ENDIF
            CALL lim_rhg( i_j1, i_jpj )

         ELSE                                 ! optimization of the computational area

            DO jj = 1, jpj
               zind(jj) = SUM( frld (:,jj  ) )   ! = FLOAT(jpj) if ocean everywhere on a j-line
               zmsk(jj) = SUM( tmask(:,jj,1) )   ! = 0          if land  everywhere on a j-line
            END DO

            IF( l_jeq ) THEN                     ! local domain include both hemisphere
               !                                 ! Rheology is computed in each hemisphere
               !                                 ! only over the ice cover latitude strip
               ! Northern hemisphere
               i_j1  = njeq
               i_jpj = jpj
               DO WHILE ( i_j1 <= jpj .AND. zind(i_j1) == FLOAT(jpi) .AND. zmsk(i_j1) /=0 )
                  i_j1 = i_j1 + 1
               END DO
               i_j1 = MAX( 1, i_j1-1 )
               IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : NH i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
    
               CALL lim_rhg( i_j1, i_jpj )
    
               ! Southern hemisphere
               i_j1  =  1 
               i_jpj = njeq
               DO WHILE ( i_jpj >= 1 .AND. zind(i_jpj) == FLOAT(jpi) .AND. zmsk(i_jpj) /=0 )
                  i_jpj = i_jpj - 1
               END DO
               i_jpj = MIN( jpj, i_jpj+2 )
               IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : SH i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
    
               CALL lim_rhg( i_j1, i_jpj )
    
            ELSE                                 ! local domain extends over one hemisphere only
               !                                 ! Rheology is computed only over the ice cover
               !                                 ! latitude strip
               i_j1  = 1
               DO WHILE ( i_j1 <= jpj .AND. zind(i_j1) == FLOAT(jpi) .AND. zmsk(i_j1) /=0 )
                  i_j1 = i_j1 + 1
               END DO
               i_j1 = MAX( 1, i_j1-1 )
    
               i_jpj  = jpj
               DO WHILE ( i_jpj >= 1  .AND. zind(i_jpj) == FLOAT(jpi) .AND. zmsk(i_jpj) /=0 )
                  i_jpj = i_jpj - 1
               END DO
               i_jpj = MIN( jpj, i_jpj+2)
    
               IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : one hemisphere: i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
    
               CALL lim_rhg( i_j1, i_jpj )

            ENDIF

         ENDIF

         IF(ln_ctl)   THEN 
            CALL prt_ctl(tab2d_1=u_oce , clinfo1=' lim_dyn  : u_oce :', tab2d_2=v_oce , clinfo2=' v_oce :')
            CALL prt_ctl(tab2d_1=u_ice , clinfo1=' lim_dyn  : u_ice :', tab2d_2=v_ice , clinfo2=' v_ice :')
         ENDIF
         
         !                                         ! Ice-Ocean stress
         !                                         ! ================
         DO jj = 2, jpjm1
            zsang  = SIGN(1.e0, gphif(1,jj-1) ) * sangvg
            DO ji = 2, jpim1
               ! computation of wind stress over ocean in X and Y direction
#if defined key_coupled && defined key_lim_cp1
               ztairx =  frld(ji-1,jj  ) * gtaux(ji-1,jj  ) + frld(ji,jj  ) * gtaux(ji,jj  )      &
                  &    + frld(ji-1,jj-1) * gtaux(ji-1,jj-1) + frld(ji,jj-1) * gtaux(ji,jj-1)

               ztairy =  frld(ji-1,jj  ) * gtauy(ji-1,jj  ) + frld(ji,jj  ) * gtauy(ji,jj  )      &
                  &    + frld(ji-1,jj-1) * gtauy(ji-1,jj-1) + frld(ji,jj-1) * gtauy(ji,jj-1)
#else
               zsfrld  = frld(ji,jj) + frld(ji-1,jj) + frld(ji-1,jj-1) + frld(ji,jj-1)
               ztairx  = zsfrld * gtaux(ji,jj)
               ztairy  = zsfrld * gtauy(ji,jj)
#endif
               zsfrldm4 = 4 - frld(ji,jj) - frld(ji-1,jj) - frld(ji-1,jj-1) - frld(ji,jj-1)
               zu_ice   = u_ice(ji,jj) - u_oce(ji,jj)
               zv_ice   = v_ice(ji,jj) - v_oce(ji,jj)
               zmod     = SQRT( zu_ice * zu_ice + zv_ice * zv_ice ) 
               ztglx   = zsfrldm4 * rhoco * zmod * ( cangvg * zu_ice - zsang * zv_ice ) 
               ztgly   = zsfrldm4 * rhoco * zmod * ( cangvg * zv_ice + zsang * zu_ice ) 

               tio_u(ji,jj) = - ( ztairx + 1.0 * ztglx ) / ( 4 * rau0 )
               tio_v(ji,jj) = - ( ztairy + 1.0 * ztgly ) / ( 4 * rau0 )
            END DO
         END DO
         
         ! computation of friction velocity
         DO jj = 2, jpjm1
            DO ji = 2, jpim1

               zu_ice   = u_ice(ji-1,jj-1) - u_oce(ji-1,jj-1)
               zv_ice   = v_ice(ji-1,jj-1) - v_oce(ji-1,jj-1)
               zt11  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice )

               zu_ice   = u_ice(ji-1,jj) - u_oce(ji-1,jj)
               zv_ice   = v_ice(ji-1,jj) - v_oce(ji-1,jj)
               zt12  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 

               zu_ice   = u_ice(ji,jj-1) - u_oce(ji,jj-1)
               zv_ice   = v_ice(ji,jj-1) - v_oce(ji,jj-1)
               zt21  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 

               zu_ice   = u_ice(ji,jj) - u_oce(ji,jj)
               zv_ice   = v_ice(ji,jj) - v_oce(ji,jj)
               zt22  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 

               ztair2 = gtaux(ji,jj) * gtaux(ji,jj) + gtauy(ji,jj) * gtauy(ji,jj)

               zustm =  ( 1 - frld(ji,jj) ) * 0.25 * ( zt11 + zt12 + zt21 + zt22 )        &
                  &  +        frld(ji,jj)   * SQRT( ztair2 )

               ust2s(ji,jj) = ( zustm / rau0 ) * ( rone + sdvt(ji,jj) ) * tms(ji,jj)
            END DO
         END DO

       ELSE      ! no ice dynamics : transmit directly the atmospheric stress to the ocean
                    
          DO jj = 2, jpjm1
             DO ji = 2, jpim1
#if defined key_coupled && defined key_lim_cp1
                tio_u(ji,jj) = - (  gtaux(ji  ,jj  ) + gtaux(ji-1,jj  )       &
                   &              + gtaux(ji-1,jj-1) + gtaux(ji  ,jj-1) ) / ( 4 * rau0 )

                tio_v(ji,jj) = - (  gtauy(ji  ,jj )  + gtauy(ji-1,jj  )       &
                   &              + gtauy(ji-1,jj-1) + gtauy(ji  ,jj-1) ) / ( 4 * rau0 )
#else
                tio_u(ji,jj) = - gtaux(ji,jj) / rau0
                tio_v(ji,jj) = - gtauy(ji,jj) / rau0 
#endif
                ztair2       = gtaux(ji,jj) * gtaux(ji,jj) + gtauy(ji,jj) * gtauy(ji,jj)
                zustm        = SQRT( ztair2  )

                ust2s(ji,jj) = ( zustm / rau0 ) * ( rone + sdvt(ji,jj) ) * tms(ji,jj)
            END DO
         END DO

      ENDIF

      CALL lbc_lnk( ust2s, 'T',  1. )   ! T-point
      CALL lbc_lnk( tio_u, 'I', -1. )   ! I-point (i.e. ice U-V point)
      CALL lbc_lnk( tio_v, 'I', -1. )   ! I-point (i.e. ice U-V point)

      IF(ln_ctl) THEN 
            CALL prt_ctl(tab2d_1=tio_u , clinfo1=' lim_dyn  : tio_u :', tab2d_2=tio_v , clinfo2=' tio_v :')
            CALL prt_ctl(tab2d_1=ust2s , clinfo1=' lim_dyn  : ust2s :')
      ENDIF

   END SUBROUTINE lim_dyn


   SUBROUTINE lim_dyn_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_dyn_init  ***
      !!
      !! ** Purpose : Physical constants and parameters linked to the ice
      !!      dynamics
      !!
      !! ** Method  :  Read the namicedyn namelist and check the ice-dynamic
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicedyn
      !!
      !! history :
      !!  8.5  ! 03-08 (C. Ethe) original code
      !!-------------------------------------------------------------------
      NAMELIST/namicedyn/ epsd, alpha,     &
         &                dm, nbiter, nbitdr, om, resl, cw, angvg, pstar,   &
         &                c_rhg, etamn, creepl, ecc, ahi0
      !!-------------------------------------------------------------------

      ! Define the initial parameters
      ! -------------------------

      ! Read Namelist namicedyn
      REWIND ( numnam_ice )
      READ   ( numnam_ice  , namicedyn )
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_dyn_init : ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '       tolerance parameter                              epsd   = ', epsd
         WRITE(numout,*) '       coefficient for semi-implicit coriolis           alpha  = ', alpha
         WRITE(numout,*) '       diffusion constant for dynamics                  dm     = ', dm
         WRITE(numout,*) '       number of sub-time steps for relaxation          nbiter = ', nbiter
         WRITE(numout,*) '       maximum number of iterations for relaxation      nbitdr = ', nbitdr
         WRITE(numout,*) '       relaxation constant                              om     = ', om
         WRITE(numout,*) '       maximum value for the residual of relaxation     resl   = ', resl
         WRITE(numout,*) '       drag coefficient for oceanic stress              cw     = ', cw
         WRITE(numout,*) '       turning angle for oceanic stress                 angvg  = ', angvg
         WRITE(numout,*) '       first bulk-rheology parameter                    pstar  = ', pstar
         WRITE(numout,*) '       second bulk-rhelogy parameter                    c_rhg  = ', c_rhg
         WRITE(numout,*) '       minimun value for viscosity                      etamn  = ', etamn
         WRITE(numout,*) '       creep limit                                      creepl = ', creepl
         WRITE(numout,*) '       eccentricity of the elliptical yield curve       ecc    = ', ecc
         WRITE(numout,*) '       horizontal diffusivity coeff. for sea-ice        ahi0   = ', ahi0
      ENDIF

      usecc2 = 1.0 / ( ecc * ecc )
      rhoco  = rau0 * cw
      angvg  = angvg * rad
      sangvg = SIN( angvg )
      cangvg = COS( angvg )
      pstarh = pstar / 2.0

      !  Diffusion coefficients.
      ahiu(:,:) = ahi0 * umask(:,:,1)
      ahiv(:,:) = ahi0 * vmask(:,:,1)

   END SUBROUTINE lim_dyn_init

#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_dyn         ! Empty routine
   END SUBROUTINE lim_dyn
#endif 

   !!======================================================================
END MODULE limdyn
