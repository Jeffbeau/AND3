MODULE limtrp
   !!======================================================================
   !!                       ***  MODULE limtrp   ***
   !! LIM transport ice model : sea-ice advection/diffusion
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_trp      : advection/diffusion process of sea ice
   !!   lim_trp_init : initialization and namelist read
   !!----------------------------------------------------------------------
   !! * Modules used
   USE phycst
   USE dom_oce
   USE daymod
   USE in_out_manager  ! I/O manager
   USE ice_oce         ! ice variables
   USE dom_ice
   USE ice
   USE iceini
   USE limistate
   USE limadv
   USE limhdf
   USE lbclnk
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC lim_trp       ! called by ice_step

   !! * Shared module variables
   REAL(wp), PUBLIC  ::   &  !:
      bound  = 0.e0          !: boundary condit. (0.0 no-slip, 1.0 free-slip)

   !! * Module variables
   REAL(wp)  ::           &  ! constant values
      epsi06 = 1.e-06  ,  &
      epsi03 = 1.e-03  ,  &
      epsi16 = 1.e-16  ,  &
      rzero  = 0.e0    ,  &
      rone   = 1.e0

   !! * Substitution
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limtrp.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_trp
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_trp ***
      !!                    
      !! ** purpose : advection/diffusion process of sea ice
      !!
      !! ** method  : variables included in the process are scalar,   
      !!     other values are considered as second order. 
      !!     For advection, a second order Prather scheme is used.  
      !!
      !! ** action :
      !!
      !! History :
      !!   1.0  !  00-01 (LIM)  Original code
      !!        !  01-05 (G. Madec, R. Hordoir) opa norm
      !!   2.0  !  04-01 (G. Madec, C. Ethe)  F90, mpp
      !!---------------------------------------------------------------------
      !! * Local Variables
      INTEGER  ::   ji, jj, jk,   &  ! dummy loop indices
         &          initad           ! number of sub-timestep for the advection

      REAL(wp) ::  &                              
         zindb  ,  &
         zacrith, &
         zindsn , &
         zindic , &
         zusvosn, &
         zusvoic, &
         zignm  , &
         zindhe , &
         zvbord , &
         zcfl   , &
         zusnit , &
         zrtt, ztsn, ztic1, ztic2

      REAL(wp), DIMENSION(jpi,jpj)  ::   &  ! temporary workspace
         zui_u , zvi_v , zsm   ,         &
         zs0ice, zs0sn , zs0a  ,         &
         zs0c0 , zs0c1 , zs0c2 ,         &
         zs0st
      !---------------------------------------------------------------------

      IF( numit == nstart  )   CALL lim_trp_init      ! Initialization (first time-step only)

      zsm(:,:) = area(:,:)
      
      IF( ln_limdyn ) THEN
         !-------------------------------------!
         !   Advection of sea ice properties   !
         !-------------------------------------!

         ! ice velocities at ocean U- and V-points (zui_u,zvi_v)
         ! ---------------------------------------
         ! zvbord factor between 1 and 2 to take into account slip or no-slip boundary conditions.        
         zvbord = 1.0 + ( 1.0 - bound )
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zui_u(ji,jj) = ( u_ice(ji+1,jj  ) + u_ice(ji+1,jj+1) ) / ( MAX( tmu(ji+1,jj  ) + tmu(ji+1,jj+1), zvbord ) )
               zvi_v(ji,jj) = ( v_ice(ji  ,jj+1) + v_ice(ji+1,jj+1) ) / ( MAX( tmu(ji  ,jj+1) + tmu(ji+1,jj+1), zvbord ) )
            END DO
         END DO
         ! Lateral boundary conditions on zui_u, zvi_v
         CALL lbc_lnk( zui_u, 'U', -1. )
         CALL lbc_lnk( zvi_v, 'V', -1. )

         ! CFL test for stability
         ! ----------------------
         zcfl  = 0.e0
         zcfl  = MAX( zcfl, MAXVAL( ABS( zui_u(1:jpim1, :     ) ) * rdt_ice / e1u(1:jpim1, :     ) ) )
         zcfl  = MAX( zcfl, MAXVAL( ABS( zvi_v( :     ,1:jpjm1) ) * rdt_ice / e2v( :     ,1:jpjm1) ) )

         IF (lk_mpp ) CALL mpp_max(zcfl)

         IF ( zcfl > 0.5 .AND. lwp )   WRITE(numout,*) 'lim_trp : violation of cfl criterion the ',nday,'th day, cfl = ',zcfl

         ! content of properties
         ! ---------------------
         zs0sn (:,:) =  hsnm(:,:) * area(:,:)                 ! Snow volume.
         zs0ice(:,:) =  hicm (:,:) * area(:,:)                ! Ice volume.
         zs0a  (:,:) =  ( 1.0 - frld(:,:) ) * area(:,:)       ! Surface covered by ice.
         zs0c0 (:,:) =  tbif(:,:,1) / rt0_snow * zs0sn(:,:)   ! Heat content of the snow layer.
         zs0c1 (:,:) =  tbif(:,:,2) / rt0_ice  * zs0ice(:,:)  ! Heat content of the first ice layer.
         zs0c2 (:,:) =  tbif(:,:,3) / rt0_ice  * zs0ice(:,:)  ! Heat content of the second ice layer.
         zs0st (:,:) =  qstoif(:,:) / xlic     * zs0a(:,:)    ! Heat reservoir for brine pockets.
         
 
         ! Advection 
         ! ---------
         ! If ice drift field is too fast, use an appropriate time step for advection.         
         initad = 1 + INT( MAX( rzero, SIGN( rone, zcfl-0.5 ) ) )
         zusnit = 1.0 / REAL( initad ) 
         
         IF ( MOD( nday , 2 ) == 0) THEN
            DO jk = 1,initad
               CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
               CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
               CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
               CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
               CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
               CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
               CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
               CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
               CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
               CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
               CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
               CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
               CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
               CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
            END DO
         ELSE
            DO jk = 1, initad
               CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
               CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
               CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
               CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
               CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
               CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
               CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
               CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
               CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
               CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
               CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
               CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
               CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
               CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
            END DO
         ENDIF
                        
         ! recover the properties from their contents
         ! ------------------------------------------
         zs0ice(:,:) = zs0ice(:,:) / area(:,:)
         zs0sn (:,:) = zs0sn (:,:) / area(:,:)
         zs0a  (:,:) = zs0a  (:,:) / area(:,:)
         zs0c0 (:,:) = zs0c0 (:,:) / area(:,:)
         zs0c1 (:,:) = zs0c1 (:,:) / area(:,:)
         zs0c2 (:,:) = zs0c2 (:,:) / area(:,:)
         zs0st (:,:) = zs0st (:,:) / area(:,:)


         !-------------------------------------!
         !   Diffusion of sea ice properties   !
         !-------------------------------------!

         ! Masked eddy diffusivity coefficient at ocean U- and V-points
         ! ------------------------------------------------------------
         DO jj = 1, jpjm1          ! NB: has not to be defined on jpj line and jpi row
            DO ji = 1 , fs_jpim1   ! vector opt.
               pahu(ji,jj) = ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji  ,jj) ) ) )   &
                  &        * ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji+1,jj) ) ) ) * ahiu(ji,jj)
               pahv(ji,jj) = ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji,jj  ) ) ) )   &
                  &        * ( 1.0 - MAX( rzero, SIGN( rone,- zs0a(ji,jj+1) ) ) ) * ahiv(ji,jj)
            END DO
         END DO

         ! diffusion
         ! ---------
         CALL lim_hdf( zs0ice )
         CALL lim_hdf( zs0sn  )
         CALL lim_hdf( zs0a   )
         CALL lim_hdf( zs0c0  )
         CALL lim_hdf( zs0c1  )
         CALL lim_hdf( zs0c2  )
         CALL lim_hdf( zs0st  )

         zs0ice(:,:) = MAX( rzero, zs0ice(:,:) * area(:,:) )    !!bug:  est-ce utile
         zs0sn (:,:) = MAX( rzero, zs0sn (:,:) * area(:,:) )    !!bug:  cf /area  juste apres
         zs0a  (:,:) = MAX( rzero, zs0a  (:,:) * area(:,:) )    !! suppression des 2 change le resultat...
         zs0c0 (:,:) = MAX( rzero, zs0c0 (:,:) * area(:,:) )
         zs0c1 (:,:) = MAX( rzero, zs0c1 (:,:) * area(:,:) )
         zs0c2 (:,:) = MAX( rzero, zs0c2 (:,:) * area(:,:) )
         zs0st (:,:) = MAX( rzero, zs0st (:,:) * area(:,:) )


         ! -------------------------------------------------------------------!
         !   Up-dating and limitation of sea ice properties after transport   !
         ! -------------------------------------------------------------------!

         ! Up-dating and limitation of sea ice properties after transport.
         DO jj = 1, jpj
!!!iii      zindhe = REAL( MAX( 0, isign(1, jj - njeq ) ) )              !ibug mpp  !!bugmpp  njeq!
            zindhe = MAX( 0.e0, SIGN( 1.e0, fcor(1,jj) ) )              ! = 0 for SH, =1 for NH
            DO ji = 1, jpi

               ! Recover mean values over the grid squares.
               zs0sn (ji,jj) = MAX( rzero, zs0sn (ji,jj)/area(ji,jj) )
               zs0ice(ji,jj) = MAX( rzero, zs0ice(ji,jj)/area(ji,jj) )
               zs0a  (ji,jj) = MAX( rzero, zs0a  (ji,jj)/area(ji,jj) )
               zs0c0 (ji,jj) = MAX( rzero, zs0c0 (ji,jj)/area(ji,jj) )
               zs0c1 (ji,jj) = MAX( rzero, zs0c1 (ji,jj)/area(ji,jj) )
               zs0c2 (ji,jj) = MAX( rzero, zs0c2 (ji,jj)/area(ji,jj) )
               zs0st (ji,jj) = MAX( rzero, zs0st (ji,jj)/area(ji,jj) )

               ! Recover in situ values.
               zindb         = MAX( rzero, SIGN( rone, zs0a(ji,jj) - epsi06 ) )
               zacrith       = 1.0 - ( zindhe * acrit(1) + ( 1.0 - zindhe ) * acrit(2) )
               zs0a (ji,jj)  = zindb * MIN( zs0a(ji,jj), zacrith )
               hsnif(ji,jj)  = zindb * ( zs0sn(ji,jj) /MAX( zs0a(ji,jj), epsi16 ) )
               hicif(ji,jj)  = zindb * ( zs0ice(ji,jj)/MAX( zs0a(ji,jj), epsi16 ) )
               zindsn        = MAX( rzero, SIGN( rone, hsnif(ji,jj) - epsi06 ) )
               zindic        = MAX( rzero, SIGN( rone, hicif(ji,jj) - epsi03 ) )
               zindb         = MAX( zindsn, zindic )
               zs0a (ji,jj)  = zindb * zs0a(ji,jj)
               frld (ji,jj)  = 1.0 - zs0a(ji,jj)
               hsnif(ji,jj)  = zindsn * hsnif(ji,jj)
               hicif(ji,jj)  = zindic * hicif(ji,jj)
               zusvosn       = 1.0/MAX( hsnif(ji,jj) * zs0a(ji,jj), epsi16 )
               zusvoic       = 1.0/MAX( hicif(ji,jj) * zs0a(ji,jj), epsi16 )
               zignm         = MAX( rzero,  SIGN( rone, hsndif - hsnif(ji,jj) ) )
               zrtt          = 173.15 * rone 
               ztsn          =          zignm   * tbif(ji,jj,1)  &
                              + ( 1.0 - zignm ) * MIN( MAX( zrtt, rt0_snow * zusvosn * zs0c0(ji,jj)) , tfu(ji,jj) ) 
               ztic1          = MIN( MAX( zrtt, rt0_ice * zusvoic * zs0c1(ji,jj) ) , tfu(ji,jj) )
               ztic2          = MIN( MAX( zrtt, rt0_ice * zusvoic * zs0c2(ji,jj) ) , tfu(ji,jj) )
 
               tbif(ji,jj,1) = zindsn * ztsn  + ( 1.0 - zindsn ) * tfu(ji,jj)               
               tbif(ji,jj,2) = zindic * ztic1 + ( 1.0 - zindic ) * tfu(ji,jj)
               tbif(ji,jj,3) = zindic * ztic2 + ( 1.0 - zindic ) * tfu(ji,jj)
               qstoif(ji,jj) = zindb  * xlic * zs0st(ji,jj) /  MAX( zs0a(ji,jj), epsi16 )
            END DO
         END DO
         
      ENDIF
      
   END SUBROUTINE lim_trp


   SUBROUTINE lim_trp_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_trp_init  ***
      !!
      !! ** Purpose :   initialization of ice advection parameters
      !!
      !! ** Method  : Read the namicetrp namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicetrp
      !!
      !! history :
      !!   2.0  !  03-08 (C. Ethe)  Original code
      !!-------------------------------------------------------------------
      NAMELIST/namicetrp/ bound
      !!-------------------------------------------------------------------

      ! Read Namelist namicetrp
      REWIND ( numnam_ice )
      READ   ( numnam_ice  , namicetrp )
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_trp_init : Ice parameters for advection '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   boundary conditions (0. no-slip, 1. free-slip) bound  = ', bound
      ENDIF
            
   END SUBROUTINE lim_trp_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty Module                No sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_trp        ! Empty routine
   END SUBROUTINE lim_trp
#endif

   !!======================================================================
END MODULE limtrp
