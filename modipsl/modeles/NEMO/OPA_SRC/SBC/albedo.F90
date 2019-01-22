MODULE albedo
   !!======================================================================
   !!                       ***  MODULE  albedo  ***
   !! Ocean forcing:  bulk thermohaline forcing of the ocean (or ice)
   !!=====================================================================
   !!----------------------------------------------------------------------
   !!   flx_blk_albedo : albedo for ocean and ice (clear and overcast skies)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE cpl_oce         ! ???
   USE phycst          ! physical constants
   USE daymod
   USE blk_oce         ! bulk variables
   USE flx_oce         ! forcings variables
   USE ocfzpt          ! ???
   USE in_out_manager
   USE lbclnk

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC flx_blk_albedo ! routine called by limflx.F90 in coupled
                         ! and in flxblk.F90 in forced
   !! * Module variables
   INTEGER  ::             &  !: nameos : ocean physical parameters
      albd_init = 0           !: control flag for initialization

   REAL(wp)  ::            &  ! constant values
      zzero   = 0.e0    ,  &
      zone    = 1.0

   !! * constants for albedo computation (flx_blk_albedo)
   REAL(wp) ::   &
      c1     = 0.05  ,     &   ! constants values
      c2     = 0.10  ,     &
      albice = 0.50  ,     &   !  albedo of melting ice in the arctic and antarctic (Shine & Hendersson-Sellers)
      cgren  = 0.06  ,     &   !  correction of the snow or ice albedo to take into account
                               !  effects of cloudiness (Grenfell & Perovich, 1984)
      alphd  = 0.80  ,     &   !  coefficients for linear interpolation used to compute
      alphdi = 0.72  ,     &   !  albedo between two extremes values (Pyane, 1972)
      alphc  = 0.65  ,     &
      zmue   = 0.40            !  cosine of local solar altitude

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/albedo.F90,v 1.3 2005/03/27 18:35:12 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim'                                         LIM ice model
   !!----------------------------------------------------------------------

   SUBROUTINE flx_blk_albedo( palb , palcn , palbp , palcnp )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE flx_blk_albedo  ***
      !!          
      !! ** Purpose :   Computation of the albedo of the snow/ice system 
      !!      as well as the ocean one
      !!       
      !! ** Method  : - Computation of the albedo of snow or ice (choose the 
      !!      rignt one by a large number of tests
      !!              - Computation of the albedo of the ocean
      !!
      !! References :
      !!      Shine and Hendersson-Sellers 1985, JGR, 90(D1), 2243-2250.
      !!
      !! History :
      !!  8.0   !  01-04  (LIM 1.0)
      !!  8.5   !  03-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ice                   ! ???

      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::  &
         palb         ,     &    !  albedo of ice under overcast sky
         palcn        ,     &    !  albedo of ocean under overcast sky
         palbp        ,     &    !  albedo of ice under clear sky 
         palcnp                  !  albedo of ocean under clear sky

      !! * Local variables
      INTEGER ::    &
         ji, jj                   ! dummy loop indices
      REAL(wp) ::   & 
         zmue14         ,     &   !  zmue**1.4
         zalbpsnm       ,     &   !  albedo of ice under clear sky when snow is melting
         zalbpsnf       ,     &   !  albedo of ice under clear sky when snow is freezing
         zalbpsn        ,     &   !  albedo of snow/ice system when ice is coverd by snow
         zalbpic        ,     &   !  albedo of snow/ice system when ice is free of snow
         zithsn         ,     &   !  = 1 for hsn >= 0 ( ice is cov. by snow ) ; = 0 otherwise (ice is free of snow)
         zitmlsn        ,     &   !  = 1 freezinz snow (sist >=rt0_snow) ; = 0 melting snow (sist<rt0_snow)
         zihsc1         ,     &   !  = 1 hsn <= c1 ; = 0 hsn > c1
         zihsc2                   !  = 1 hsn >= c2 ; = 0 hsn < c2
      REAL(wp), DIMENSION(jpi,jpj) ::  &
         zalbfz         ,     &   !  ( = alphdi for freezing ice ; = albice for melting ice )
         zficeth                  !  function of ice thickness
      LOGICAL , DIMENSION(jpi,jpj) ::  &
         llmask
      !!---------------------------------------------------------------------
      
      ! initialization 
      IF( albd_init == 0 )   CALL albedo_init

      !-------------------------                                                             
      !  Computation of  zficeth
      !-------------------------- 
      
      llmask = (hsnif == 0.e0) .AND. ( sist >= rt0_ice )
      WHERE ( llmask )   !  ice free of snow and melts
         zalbfz = albice
      ELSEWHERE                   
         zalbfz = alphdi
      END WHERE
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( hicif(ji,jj) > 1.5 ) THEN
               zficeth(ji,jj) = zalbfz(ji,jj)
            ELSEIF( hicif(ji,jj) > 1.0  .AND. hicif(ji,jj) <= 1.5 ) THEN
               zficeth(ji,jj) = 0.472 + 2.0 * ( zalbfz(ji,jj) - 0.472 ) * ( hicif(ji,jj) - 1.0 )
            ELSEIF( hicif(ji,jj) > 0.05 .AND. hicif(ji,jj) <= 1.0 ) THEN
               zficeth(ji,jj) = 0.2467 + 0.7049 * hicif(ji,jj)                                &
                  &                    - 0.8608 * hicif(ji,jj) * hicif(ji,jj)                 &
                  &                    + 0.3812 * hicif(ji,jj) * hicif(ji,jj) * hicif (ji,jj)
            ELSE
               zficeth(ji,jj) = 0.1 + 3.6 * hicif(ji,jj) 
            ENDIF
         END DO
      END DO
      
      !----------------------------------------------- 
      !    Computation of the snow/ice albedo system 
      !-------------------------- ---------------------
      
      !    Albedo of snow-ice for clear sky.
      !-----------------------------------------------    
      DO jj = 1, jpj
         DO ji = 1, jpi
            !  Case of ice covered by snow.             
            
            !  melting snow        
            zihsc1       = 1.0 - MAX ( zzero , SIGN ( zone , - ( hsnif(ji,jj) - c1 ) ) )
            zalbpsnm     = ( 1.0 - zihsc1 ) * ( zficeth(ji,jj) + hsnif(ji,jj) * ( alphd - zficeth(ji,jj) ) / c1 ) &
               &                 + zihsc1   * alphd  
            !  freezing snow                
            zihsc2       = MAX ( zzero , SIGN ( zone , hsnif(ji,jj) - c2 ) )
            zalbpsnf     = ( 1.0 - zihsc2 ) * ( albice + hsnif(ji,jj) * ( alphc - albice ) / c2 )                 &
               &                 + zihsc2   * alphc 
            
            zitmlsn      =  MAX ( zzero , SIGN ( zone , sist(ji,jj) - rt0_snow ) )   
            zalbpsn      =  zitmlsn * zalbpsnf + ( 1.0 - zitmlsn ) * zalbpsnm 
            
            !  Case of ice free of snow.
            zalbpic      = zficeth(ji,jj) 
            
            ! albedo of the system   
            zithsn       = 1.0 - MAX ( zzero , SIGN ( zone , - hsnif(ji,jj) ) )
            palbp(ji,jj) =  zithsn * zalbpsn + ( 1.0 - zithsn ) *  zalbpic
         END DO
      END DO
      
      !    Albedo of snow-ice for overcast sky.
      !----------------------------------------------  
      palb(:,:)   = palbp(:,:) + cgren                                           
      
      !--------------------------------------------
      !    Computation of the albedo of the ocean 
      !-------------------------- -----------------                                                          
      
      !  Parameterization of Briegled and Ramanathan, 1982 
      zmue14      = zmue**1.4                                       
      palcnp(:,:) = 0.05 / ( 1.1 * zmue14 + 0.15 )                
!DL test jan 2013, mais est recalcul/ dans flx_blk, modifi/ cette routine ...
      palcnp(:,:)=palcnp(:,:)*1.2
      
      !  Parameterization of Kondratyev, 1969 and Payne, 1972
!DL      palcn(:,:)  = 0.06
      palcn(:,:)  = 0.10
      
   END SUBROUTINE flx_blk_albedo

# else
   !!----------------------------------------------------------------------
   !!   Default option :                                   NO sea-ice model
   !!----------------------------------------------------------------------

   SUBROUTINE flx_blk_albedo( palb , palcn , palbp , palcnp )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE flx_blk_albedo  ***
      !! 
      !! ** Purpose :   Computation of the albedo of the snow/ice system
      !!      as well as the ocean one
      !!
      !! ** Method  :   Computation of the albedo of snow or ice (choose the
      !!      wright one by a large number of tests Computation of the albedo
      !!      of the ocean
      !!
      !! History :
      !!  8.0   !  01-04  (LIM 1.0)
      !!  8.5   !  03-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::  &
         palb         ,     &    !  albedo of ice under overcast sky
         palcn        ,     &    !  albedo of ocean under overcast sky
         palbp        ,     &    !  albedo of ice under clear sky
         palcnp                  !  albedo of ocean under clear sky

      REAL(wp) ::   &
         zmue14                 !  zmue**1.4
      !!----------------------------------------------------------------------

      !--------------------------------------------
      !    Computation of the albedo of the ocean
      !-------------------------- -----------------

      !  Parameterization of Briegled and Ramanathan, 1982
      zmue14      = zmue**1.4
      palcnp(:,:) = 0.05 / ( 1.1 * zmue14 + 0.15 )

!DL test jan 2013
      palcnp(:,:)=palcnp(:,:)*1.25

      !  Parameterization of Kondratyev, 1969 and Payne, 1972
!DL      palcn(:,:)  = 0.06
      palcn(:,:)  = 0.10

      palb (:,:)  = palcn(:,:)
      palbp(:,:)  = palcnp(:,:)

   END SUBROUTINE flx_blk_albedo

#endif

   SUBROUTINE albedo_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE albedo_init  ***
      !!
      !! ** Purpose :   initializations for the albedo parameters
      !!
      !! ** Method  :   Read the namelist namalb
      !!
      !! ** Action  :  
      !!
      !!
      !! History :
      !!   9.0  !  04-11  (C. Talandier)  Original code
      !!----------------------------------------------------------------------
      NAMELIST/namalb/ cgren, albice, alphd, alphdi, alphc
      !!----------------------------------------------------------------------
      !!  OPA 9.0, LODYC-IPSL (2004)
      !!----------------------------------------------------------------------

      ! set the initialization flag to 1
      albd_init = 1           ! indicate that the initialization has been done

      ! Read Namelist namalb : albedo parameters
      REWIND( numnam )
      READ  ( numnam, namalb )

      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'albedo_init : albedo '
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist namalb : set albedo parameters'
         WRITE(numout,*)
         WRITE(numout,*) '             correction of the snow or ice albedo to take into account cgren = ', cgren
         WRITE(numout,*) '             albedo of melting ice in the arctic and antarctic        albice = ', albice
         WRITE(numout,*) '             coefficients for linear                                   alphd = ', alphd
         WRITE(numout,*) '             interpolation used to compute albedo                     alphdi = ', alphdi
         WRITE(numout,*) '             between two extremes values (Pyane, 1972)                 alphc = ', alphc
         WRITE(numout,*)
      ENDIF

   END SUBROUTINE albedo_init
   !!======================================================================
END MODULE albedo
