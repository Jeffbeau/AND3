MODULE trdmld
   !!======================================================================
   !!                       ***  MODULE  trdmld  ***
   !! Ocean diagnostics:  mixed layer T-S trends 
   !!=====================================================================
#if   defined key_trdmld   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_trdmld'                          mixed layer trend diagnostics
   !!----------------------------------------------------------------------
   !!   trd_mld          : T and S cumulated trends averaged over the mixed layer
   !!   trd_mld_zint     : T and S trends vertical integration
   !!   trd_mld_init     : initialization step
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE trdmod_oce      ! ocean variables trends
   USE ldftra_oce      ! ocean active tracers lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE phycst          ! Define parameters for the routines
   USE daymod          ! calendar
   USE dianam          ! build the name of file (routine)
   USE ldfslp          ! iso-neutral slopes 
   USE zdfmxl          ! mixed layer depth
   USE zdfddm          ! ocean vertical physics: double diffusion
   USE ioipsl          ! NetCDF library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
!   USE diadimg         ! dimg direct access file format output

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trd_mld        ! routine called by step.F90
   PUBLIC trd_mld_init   ! routine called by opa.F90
   PUBLIC trd_mld_zint   ! routine called by tracers routines

   !! * Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_trdmld = .TRUE.    !: momentum trend flag

   !! * Module variables
   INTEGER ::   &
      nh_t, nmoymltrd,             &  ! ???
      nidtrd,                      &
      ndextrd1(jpi*jpj),           &
      ndimtrd1
   INTEGER, SAVE ::   &
      ionce, icount,               &
      idebug                          ! (0/1) set it to 1 in case of problem to have more print

   INTEGER, DIMENSION(jpi,jpj) ::   &
      nmld,                         & ! mixed layer depth
      nbol                

   REAL(wp), DIMENSION(jpi,jpj) ::   &
      rmld   ,          &  ! mld depth (m) corresponding to nmld
      tml    , sml  ,   &  ! average T and S over mixed layer
      tmlb   , smlb ,   &  ! before tml and sml (kt-1)
      tmlbb  , smlbb,   &  ! tml and sml at begining of the nwrite-1 
      !                    ! timestep averaging period
      tmlbn  , smlbn,   &  ! after tml and sml at time step after the
      !                    ! begining of the NWRITE-1 timesteps
      tmltrdm, smltrdm     !

   REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
      tmltrd ,          &  ! total cumulative trends of temperature and 
      smltrd ,          &  ! salinity over nwrite-1 time steps
      wkx

   CHARACTER(LEN=80) :: clname

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmld.F90,v 1.7 2005/12/12 14:18:08 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

SUBROUTINE trd_mld_zint( pttrdmld, pstrdmld, ktrd, ctype )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mld_zint  ***
      !! 
      !! ** Purpose :   computation of vertically integrated T and S budgets
      !!                from ocean surface down to control surface 
      !!
      !! ** Method/usage :
      !!      integration done over nwrite-1 time steps 
      !!      Control surface can be either a mixed layer depth (time varying)
      !!      or a fixed surface (jk level or bowl). 
      !!      Choose control surface with nctls in namelist NAMDIA.
      !!      nctls = 0  : use mixed layer with density criterion 
      !!      nctls = 1  : read index from file 'ctlsurf_idx'
      !!      nctls > 1  : use fixed level surface jk = nctls
      !!      Note: in the remainder of the routine, the volume between the 
      !!            surface and the control surface is called "mixed-layer"
      !!      Method check : if the control surface is fixed, the residual dh/dt
      !!                     entrainment should be zero
      !!
      !! ** Action :
      !!            /commld/   : rmld         mld depth corresponding to nmld
      !!                         tml          average T over mixed layer
      !!                         tmlb         tml at kt-1
      !!                         tmlbb        tml at begining of the NWRITE-1 
      !!                                      time steps averaging period
      !!                         tmlbn        tml at time step after the 
      !!                                      begining of the NWRITE-1 time
      !!                                      steps averaging period
      !!
      !!                  mixed layer trends :
      !!
      !!                  tmltrd (,,1) = zonal advection
      !!                  tmltrd (,,2) = meridional advection
      !!                  tmltrd (,,3) = vertical advection
      !!                  tmltrd (,,4) = lateral diffusion (horiz. component+Beckman)
      !!                  tmltrd (,,5) = forcing
      !!                  tmltrd (,,6) = entrainment due to vertical diffusion (TKE)
      !!          if iso  tmltrd (,,7) = lateral diffusion (vertical component)
      !!                  tmltrd (,,8) = eddy induced zonal advection
      !!                  tmltrd (,,9) = eddy induced meridional advection
      !!                  tmltrd (,,10) = eddy induced vertical advection
      !!
      !!                  tmltrdm(,) : total cumulative trends over nwrite-1 time steps
      !!                  ztmltot(,) : dT/dt over the NWRITE-1 time steps 
      !!                               averaging period (including Asselin 
      !!                               terms)
      !!                  ztmlres(,) : residual = dh/dt entrainment
      !!
      !!      trends output in netCDF format using ioipsl
      !!
      !! History :
      !!        !  95-04  (J. Vialard)  Original code
      !!        !  97-02  (E. Guilyardi)  Adaptation global + base cmo
      !!        !  99-09  (E. Guilyardi)  Re-writing + netCDF output
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index

      CHARACTER(len=2), INTENT( in ) ::   &
         ctype                                ! surface/bottom (2D arrays) or
                                              ! interior (3D arrays) physics

      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
         pttrdmld,                         &  ! Temperature trend 
         pstrdmld                             ! Salinity    trend

      !! * Local declarations
      INTEGER ::   ji, jj, jk, isum
# if defined key_trabbl_dif
      INTEGER ::   ikb
# endif

      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zvlmsk
      !!----------------------------------------------------------------------

      IF( icount == 1 ) THEN        

         zvlmsk(:,:)   = 0.e0
         tmltrd(:,:,:) = 0.e0
         smltrd(:,:,:) = 0.e0
         
         ! This computation should be done only once per time step

         !  ========================================================
         !   I. definition of control surface and associated fields
         !  ========================================================

         !    I.1 set nmld(ji,jj) = index of first T point below control surface
         !    -------------------                       or outside mixed-layer

         IF( nctls == 0 ) THEN
            ! control surface = mixed-layer with density criterion 
            ! (array nmln computed in zdfmxl.F90)
            nmld(:,:) = nmln(:,:)
         ELSE IF( nctls == 1 ) THEN
            ! control surface = read index from file 
            nmld(:,:) = nbol(:,:)
         ELSE IF( nctls >= 2 ) THEN
            ! control surface = model level
            nctls = MIN( nctls, jpktrd - 1 )
            nmld(:,:) = nctls + 1
         ENDIF

         IF( ionce == 1 ) THEN  ! compute ndextrd1 and ndimtrd1 only once
            ! Check of validity : nmld(ji,jj) =< jpktrd
            isum = 0

            IF( jpktrd < jpk ) THEN 
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     IF( nmld(ji,jj) <= jpktrd ) THEN
                        zvlmsk(ji,jj) = tmask(ji,jj,1)
                     ELSE
                        isum = isum + 1
                        zvlmsk(ji,jj) = 0.
                     ENDIF
                  END DO
               END DO
            ENDIF

            ! Index of ocean points (2D only)
            IF( isum > 0 ) THEN
               if(lwp) WRITE(numout,*)' Number of invalid points nmld > jpktrd', isum 
               CALL wheneq( jpi*jpj, zvlmsk(:,:) , 1, 1., ndextrd1, ndimtrd1 )    ! surface
            ELSE 
               CALL wheneq( jpi*jpj, tmask(:,:,1), 1, 1., ndextrd1, ndimtrd1 )    ! surface
            ENDIF                                

            ! no more pass here
            ionce = 0

         ENDIF
         
         IF( idebug /= 0 ) THEN
            ! CALL prihre (zvlmsk,jpi,jpj,1,jpi,2,1,jpj,2,3,numout)
            WRITE(numout,*) ' debuging trd_mld_zint: I.1 done '  
            CALL FLUSH(numout)
         ENDIF


         ! I.2 probability density function of presence in mixed-layer
         ! --------------------------------
         ! (i.e. weight of each grid point in vertical integration : wkx(ji,jj,jk)


         ! initialize wkx with vertical scale factor in mixed-layer

         wkx(:,:,:) = 0.e0
         DO jk = 1, jpktrd
            DO jj = 1,jpj
               DO ji = 1,jpi
                  IF( jk - nmld(ji,jj) < 0. )   wkx(ji,jj,jk) = fse3t(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         
         ! compute mixed-layer depth : rmld
         
         rmld(:,:) = 0.
         DO jk = 1, jpktrd
            rmld(:,:) = rmld(:,:) + wkx(:,:,jk)
         END DO
         
         ! compute PDF

         DO jk = 1, jpktrd
            wkx(:,:,jk) = wkx(:,:,jk) / MAX( 1., rmld(:,:) )
         END DO

         IF( idebug /= 0 ) THEN
             if(lwp) WRITE(numout,*) ' debuging trd_mld_zint: I.2 done '  
            CALL FLUSH(numout)
         ENDIF

         ! Set counter icount to 0 to avoid this part at each time step
         icount = 0

      ENDIF


      !  ====================================================
      !   II. vertical integration of trends in mixed-layer
      !  ====================================================

      ! II.1 vertical integration of 3D and 2D trends
      ! ---------------------------------------------

      SELECT CASE (ctype)

      CASE ('3D')       ! 3D treatment

         ! trends terms in the mixed-layer
         DO jk = 1, jpktrd
            ! Temperature
            tmltrd(:,:,ktrd) = tmltrd(:,:,ktrd) + pttrdmld(:,:,jk) * wkx(:,:,jk)   

            ! Salinity
            smltrd(:,:,ktrd) = smltrd(:,:,ktrd) + pstrdmld(:,:,jk) * wkx(:,:,jk)   
         ENDDO

      CASE ('2D')       ! 2D treatment

         SELECT CASE (ktrd) 

         CASE (jpmldldf)

# if defined key_trabbl_dif
               ! trends terms from Beckman over-flow parameterization
               DO jj = 1,jpj
                  DO ji = 1,jpi
                     ikb = MAX( mbathy(ji,jj)-1, 1 )
                     ! beckmann component -> horiz. part of lateral diffusion
                     tmltrd(ji,jj,ktrd) = tmltrd(ji,jj,ktrd) + pttrdmld(ji,jj,1) * wkx(ji,jj,ikb)
                     smltrd(ji,jj,ktrd) = smltrd(ji,jj,ktrd) + pstrdmld(ji,jj,1) * wkx(ji,jj,ikb)
                  END DO
               END DO
# endif

         CASE DEFAULT

            ! trends terms at upper boundary of mixed-layer

            ! forcing term (non penetrative)
            ! Temperature
            tmltrd(:,:,ktrd) = tmltrd(:,:,ktrd) + pttrdmld(:,:,1) * wkx(:,:,1)   

            ! forcing term
            ! Salinity
            smltrd(:,:,ktrd) = smltrd(:,:,ktrd) + pstrdmld(:,:,1) * wkx(:,:,1)   

         END SELECT

      END SELECT

      IF( idebug /= 0 ) THEN
         IF(lwp) WRITE(numout,*) ' debuging trd_mld_zint: II.1 done'  
         CALL FLUSH(numout)
      ENDIF

   END SUBROUTINE trd_mld_zint



   SUBROUTINE trd_mld( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mld  ***
      !! 
      !! ** Purpose :  computation of cumulated trends over analysis period
      !!               and make outputs (NetCDF or DIMG format)
      !!
      !! ** Method/usage :
      !!
      !! History :
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      !! * Local declarations
      INTEGER :: ji, jj, jk, jl, ik, it

      REAL(wp) :: zmean, zavt

      REAL(wp) ,DIMENSION(jpi,jpj) ::   &
         ztmltot, ztmlres,              &
         zsmltot, zsmlres,              & 
         z2d

!#if defined key_dimgout
!      INTEGER ::  iyear,imon,iday
!      CHARACTER(LEN=80) :: cltext, clmode
!#endif
      !!----------------------------------------------------------------------

      ! I. trends terms at lower boundary of mixed-layer
      ! ------------------------------------------------

      DO jj = 1,jpj
         DO ji = 1,jpi
            
            ik = nmld(ji,jj)
            
            ! Temperature
            ! entrainment due to vertical diffusion
            !       - due to vertical mixing scheme (TKE)
            zavt = avt(ji,jj,ik)
            tmltrd(ji,jj,jpmldevd) = - 1. * zavt / fse3w(ji,jj,ik) * tmask(ji,jj,ik)   &
               &                   * ( tn(ji,jj,ik-1) - tn(ji,jj,ik) )   &
               &                   / MAX( 1., rmld(ji,jj) ) * tmask(ji,jj,1)
            ! Salinity
            ! entrainment due to vertical diffusion
            !       - due to vertical mixing scheme (TKE)
            zavt = fsavs(ji,jj,ik)
            smltrd(ji,jj,jpmldevd) = -1. * zavt / fse3w(ji,jj,ik) * tmask(ji,jj,ik)   &
               &                  * ( sn(ji,jj,ik-1) - sn(ji,jj,ik) )   &
               &                  / MAX( 1., rmld(ji,jj) ) * tmask(ji,jj,1)
         END DO
      END DO

      IF( ln_traldf_iso ) THEN
         ! We substract to the TOTAL vertical diffusion tmltrd(:,:,jpmldzdf) 
         ! computed in subroutines trazdf_iso.F90 or trazdf_imp.F90
         ! the vertical part du to the Kz in order to keep only the vertical
         ! isopycnal diffusion (i.e the isopycnal diffusion componant on the vertical):
         tmltrd(:,:,jpmldzdf) = tmltrd(:,:,jpmldzdf) - tmltrd(:,:,jpmldevd)   ! - due to isopycnal mixing scheme (implicit part)
         smltrd(:,:,jpmldzdf) = smltrd(:,:,jpmldzdf) - smltrd(:,:,jpmldevd)   ! - due to isopycnal mixing scheme (implicit part)
      ENDIF

      ! Boundary conditions
      CALL lbc_lnk( tmltrd, 'T', 1. )
      CALL lbc_lnk( smltrd, 'T', 1. )

      IF( idebug /= 0 ) THEN
          if(lwp) WRITE(numout,*) ' debuging trd_mld: I. done'  
         CALL FLUSH(numout)
      ENDIF

      !  =================================
      !   II. Cumulated trends
      !  =================================

      ! II.1 set before values of vertically average T and S 
      ! ---------------------------------------------------

      IF( kt > nit000 ) THEN
         tmlb(:,:) = tml(:,:)
         smlb(:,:) = sml(:,:)
      ENDIF

      ! II.2 vertically integrated T and S
      ! ---------------------------------

      tml(:,:) = 0.
      sml(:,:) = 0.

      DO jk = 1, jpktrd - 1
         tml(:,:) = tml(:,:) + wkx(:,:,jk) * tn(:,:,jk)
         sml(:,:) = sml(:,:) + wkx(:,:,jk) * sn(:,:,jk) 
      END DO

      IF(idebug /= 0) THEN
         if(lwp) WRITE(numout,*) ' debuging trd_mld: II.2 done'  
         CALL FLUSH(numout)
      ENDIF

      ! II.3 set `before' mixed layer values for kt = nit000+1
      ! --------------------------------------------------------

      IF( kt == nit000+1 ) THEN
         tmlbb(:,:) = tmlb(:,:)
         tmlbn(:,:) = tml (:,:)
         smlbb(:,:) = smlb(:,:)
         smlbn(:,:) = sml (:,:)
      ENDIF

      IF( idebug /= 0 ) THEN
         if(lwp) WRITE(numout,*) ' debuging trd_mld: II.3 done'  
         CALL FLUSH(numout)
      ENDIF

      ! II.4 cumulated trends over analysis period (kt=2 to nwrite)
      ! -----------------------------------------------------------

      ! trends cumulated over nwrite-2 time steps

      IF( kt >= nit000+2 ) THEN
         nmoymltrd = nmoymltrd + 1
         DO jl = 1, jpltrd
            tmltrdm(:,:) = tmltrdm(:,:) + tmltrd(:,:,jl)
            smltrdm(:,:) = smltrdm(:,:) + smltrd(:,:,jl)
         END DO
      ENDIF

      IF( idebug /= 0 ) THEN
         if(lwp) WRITE(numout,*) ' debuging trd_mld: II.4 done'  
         CALL FLUSH(numout)
      ENDIF

      !  =============================================
      !   III. Output in netCDF + residual computation
      !  =============================================

      ztmltot(:,:) = 0.
      zsmltot(:,:) = 0.
      ztmlres(:,:) = 0.
      zsmlres(:,:) = 0.

      IF( MOD( kt - nit000+1, nwrite ) == 0 ) THEN

         ! III.1 compute total trend 
         ! ------------------------

         zmean = float(nmoymltrd)
         
         ztmltot(:,:) = ( tml(:,:) - tmlbn(:,:) + tmlb(:,:) - tmlbb(:,:) ) /  (zmean * 2. * rdt)
         zsmltot(:,:) = ( sml(:,:) - smlbn(:,:) + smlb(:,:) - smlbb(:,:) ) /  (zmean * 2. * rdt)

         IF(idebug /= 0) THEN
            if(lwp) WRITE(numout,*) ' zmean = ',zmean  
            if(lwp) WRITE(numout,*) ' debuging trd_mld: III.1 done'  
            CALL FLUSH(numout)
         ENDIF
          

         ! III.2 compute residual 
         ! ---------------------

         ztmlres(:,:) = ztmltot(:,:) - tmltrdm(:,:) / zmean
         zsmlres(:,:) = zsmltot(:,:) - smltrdm(:,:) / zmean


         ! Boundary conditions

         CALL lbc_lnk( ztmltot, 'T', 1. )
         CALL lbc_lnk( ztmlres, 'T', 1. )
         CALL lbc_lnk( zsmltot, 'T', 1. )
         CALL lbc_lnk( zsmlres, 'T', 1. )

         IF( idebug /= 0 ) THEN
            if(lwp) WRITE(numout,*) ' debuging trd_mld: III.2 done'  
            CALL FLUSH(numout)
         ENDIF


         ! III.3 time evolution array swap
         ! ------------------------------

         tmlbb(:,:) = tmlb(:,:)
         tmlbn(:,:) = tml (:,:)
         smlbb(:,:) = smlb(:,:)
         smlbn(:,:) = sml (:,:)

         IF( idebug /= 0 ) THEN
            if(lwp) WRITE(numout,*) ' debuging trd_mld: III.3 done'  
            CALL FLUSH(numout)
         ENDIF


         ! III.4 zero cumulative array
         ! ---------------------------

          nmoymltrd = 0

          tmltrdm(:,:) = 0.
          smltrdm(:,:) = 0.

          IF(idebug /= 0) THEN
             if(lwp) WRITE(numout,*) ' debuging trd_mld: III.4 done'  
             CALL FLUSH(numout)
          ENDIF
          
      ENDIF

      ! III.5 write trends to output
      ! ---------------------------

      IF( kt >=  nit000+1 ) THEN

         ! define time axis
         it= kt-nit000+1
         IF( lwp .AND. MOD( kt, nwrite ) == 0 ) THEN
            WRITE(numout,*) '     trd_mld : write NetCDF fields'
         ENDIF
         
         CALL histwrite( nidtrd,"somlttml",it,rmld          ,ndimtrd1,ndextrd1) ! Mixed-layer depth
         
         ! Temperature trends
         ! ------------------
         CALL histwrite( nidtrd,"somltemp",it,tml           ,ndimtrd1,ndextrd1) ! Mixed-layer temperature
         CALL histwrite( nidtrd,"somlttto",it,ztmltot       ,ndimtrd1,ndextrd1) ! total 
         CALL histwrite( nidtrd,"somlttax",it,tmltrd(:,:, 1),ndimtrd1,ndextrd1) ! i- adv.
         CALL histwrite( nidtrd,"somlttay",it,tmltrd(:,:, 2),ndimtrd1,ndextrd1) ! j- adv.
         CALL histwrite( nidtrd,"somlttaz",it,tmltrd(:,:, 3),ndimtrd1,ndextrd1) ! vertical adv.
         CALL histwrite( nidtrd,"somlttdh",it,tmltrd(:,:, 4),ndimtrd1,ndextrd1) ! hor. lateral diff.
         CALL histwrite( nidtrd,"somlttfo",it,tmltrd(:,:, 5),ndimtrd1,ndextrd1) ! forcing

         CALL histwrite( nidtrd,"somlbtdz",it,tmltrd(:,:, 6),ndimtrd1,ndextrd1) ! vert. diffusion 
         CALL histwrite( nidtrd,"somlbtdt",it,ztmlres       ,ndimtrd1,ndextrd1) ! dh/dt entrainment (residual)
         IF( ln_traldf_iso ) THEN
            CALL histwrite( nidtrd,"somlbtdv",it,tmltrd(:,:, 7),ndimtrd1,ndextrd1) ! vert. lateral diff.
         ENDIF
#if defined key_traldf_eiv
         CALL histwrite( nidtrd,"somlgtax",it,tmltrd(:,:, 8),ndimtrd1,ndextrd1) ! i- adv. (eiv)
         CALL histwrite( nidtrd,"somlgtay",it,tmltrd(:,:, 9),ndimtrd1,ndextrd1) ! j- adv. (eiv)
         CALL histwrite( nidtrd,"somlgtaz",it,tmltrd(:,:,10),ndimtrd1,ndextrd1) ! vert. adv. (eiv)
         z2d(:,:) = tmltrd(:,:,8) + tmltrd(:,:,9) + tmltrd(:,:,10)
         CALL histwrite( nidtrd,"somlgtat",it,z2d           ,ndimtrd1,ndextrd1) ! total adv. (eiv)
#endif   

         ! Salinity trends
         ! ---------------
         CALL histwrite( nidtrd,"somlsalt",it,sml           ,ndimtrd1,ndextrd1) ! Mixed-layer salinity
         CALL histwrite( nidtrd,"somltsto",it,zsmltot       ,ndimtrd1,ndextrd1) ! total 
         CALL histwrite( nidtrd,"somltsax",it,smltrd(:,:, 1),ndimtrd1,ndextrd1) ! i- adv.
         CALL histwrite( nidtrd,"somltsay",it,smltrd(:,:, 2),ndimtrd1,ndextrd1) ! j- adv.
         CALL histwrite( nidtrd,"somltsaz",it,smltrd(:,:, 3),ndimtrd1,ndextrd1) ! vert. adv.
         CALL histwrite( nidtrd,"somltsdh",it,smltrd(:,:, 4),ndimtrd1,ndextrd1) ! hor. lateral diff.
         CALL histwrite( nidtrd,"somltsfo",it,smltrd(:,:, 5),ndimtrd1,ndextrd1) ! forcing
         CALL histwrite( nidtrd,"somlbsdz",it,smltrd(:,:, 6),ndimtrd1,ndextrd1) ! vert. diff.
         CALL histwrite( nidtrd,"somlbsdt",it,zsmlres       ,ndimtrd1,ndextrd1) ! dh/dt entrainment (residual)
         IF( ln_traldf_iso ) THEN
            CALL histwrite( nidtrd,"somlbsdv",it,smltrd(:,:, 7),ndimtrd1,ndextrd1) ! vert. lateral diff;
         ENDIF
#if defined key_traldf_eiv
         CALL histwrite( nidtrd,"somlgsax",it,smltrd(:,:, 8),ndimtrd1,ndextrd1) ! i-adv. (eiv)
         CALL histwrite( nidtrd,"somlgsay",it,smltrd(:,:, 9),ndimtrd1,ndextrd1) ! j-adv. (eiv)
         CALL histwrite( nidtrd,"somlgsaz",it,smltrd(:,:,10),ndimtrd1,ndextrd1) ! vert. adv. (eiv)
         z2d(:,:) = smltrd(:,:,8) + smltrd(:,:,9) + smltrd(:,:,10)
         CALL histwrite( nidtrd,"somlgsat",it,z2d           ,ndimtrd1,ndextrd1) ! total adv. (eiv)
#endif

         IF( idebug /= 0 ) THEN
            if(lwp) WRITE(numout,*) ' debuging trd_mld: III.5 done'  
            CALL FLUSH(numout)
         ENDIF

         ! set counter icount to one to allow the calculation
         ! of the surface control in the next time step in the trd_mld_zint subroutine
         icount = 1

      ENDIF

      ! At the end of the 1st time step, set icount to 1 to be
      ! able to compute the surface control at the beginning of
      ! the second time step
      IF( kt == nit000 )   icount = 1

      IF( kt == nitend )   CALL histclo( nidtrd )

   END SUBROUTINE trd_mld



   SUBROUTINE trd_mld_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mld_init  ***
      !! 
      !! ** Purpose :   computation of vertically integrated T and S budgets
      !!      from ocean surface down to control surface (NetCDF output)
      !!
      !! ** Method/usage :
      !!
      !! History :
      !!        !  95-04  (J. Vialard)  Original code
      !!        !  97-02  (E. Guilyardi)  Adaptation global + base cmo
      !!        !  99-09  (E. Guilyardi)  Re-writing + netCDF output
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER :: ilseq

      REAL(wp) ::   zjulian, zsto, zout

      CHARACTER (LEN=21) ::   &
         clold ='OLD'        , & ! open specifier (direct access files)
         clunf ='UNFORMATTED', & ! open specifier (direct access files)
         clseq ='SEQUENTIAL'     ! open specifier (direct access files)
      CHARACTER (LEN=40) ::   clhstnam
      CHARACTER (LEN=40) ::   clop
      CHARACTER (LEN=12) ::   clmxl

      NAMELIST/namtrd/ ntrd, nctls
      !!----------------------------------------------------------------------

      !  ===================
      !   I. initialization
      !  ===================

      ! Open specifier
      ilseq  = 1
      idebug = 0      ! set it to 1 in case of problem to have more print
      icount = 1      
      ionce  = 1

      ! namelist namtrd : trend diagnostic
      REWIND( numnam )
      READ  ( numnam, namtrd )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) 'trd_mld_init: mixed layer heat & freshwater budget trends'
         WRITE(numout,*) '~~~~~~~~~~~~~'
         WRITE(numout,*) ' '
         WRITE(numout,*) '          Namelist namtrd : '
         WRITE(numout,*) '             control surface for trends      nctls = ',nctls
         WRITE(numout,*) ' '
      ENDIF

      ! cumulated trends array init
      nmoymltrd = 0
      tmltrdm(:,:) = 0.e0
      smltrdm(:,:) = 0.e0

      !  read control surface from file ctlsurf_idx

      IF( nctls == 1 ) THEN
         clname ='ctlsurf_idx'
         CALL ctlopn(numbol,clname,clold,clunf,clseq,   &
              ilseq,numout,lwp,1)
         REWIND (numbol)
         READ(numbol) nbol
      ENDIF


      IF( idebug /= 0 ) THEN
         if(lwp) WRITE(numout,*) ' debuging trd_mld_init: 0. done '  
         CALL FLUSH(numout)
      ENDIF

      !  ===================================
      !   II. netCDF output initialization
      !  ===================================

      !     clmxl = legend root for netCDF output
      IF( nctls == 0 ) THEN
         ! control surface = mixed-layer with density criterion 
         ! (array nmln computed in zdfmxl.F90)
         clmxl = 'Mixed Layer '
      ELSE IF( nctls == 1 ) THEN
         ! control surface = read index from file 
         clmxl = '      Bowl '
      ELSE IF( nctls >= 2 ) THEN
         ! control surface = model level
         WRITE(clmxl,'(A9,I2,1X)') 'Levels 1-', nctls
      ENDIF

      !-----------------------------------------
      ! II.1 Define frequency of output and means
      ! -----------------------------------------

#if defined key_diainstant
      zsto = nwrite*rdt
      clop ="inst(x)"
#else
      zsto = rdt
      clop ="ave(x)"
#endif
      zout = nwrite*rdt

      IF(lwp) WRITE (numout,*) ' trdmld_ncinit: netCDF initialization'

      ! II.2 Compute julian date from starting date of the run
      ! ------------------------

      CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
      IF (lwp) WRITE(numout,*)' '  
      IF (lwp) WRITE(numout,*)' Date 0 used :',nit000   &
           ,' YEAR ', nyear,' MONTH ', nmonth,' DAY ', nday   &
           ,'Julian day : ', zjulian


      ! II.3 Define the T grid trend file (nidtrd)
      ! ---------------------------------

      CALL dia_nam( clhstnam, nwrite, 'trends' )                  ! filename
      IF(lwp) WRITE(numout,*) ' Name of NETCDF file ', clhstnam
      CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,1, jpi,   &  ! Horizontal grid : glamt and gphit
         &          1, jpj, 0, zjulian, rdt, nh_t, nidtrd, domain_id=nidom )

      ! Declare output fields as netCDF variables

      ! Mixed layer Depth
      CALL histdef( nidtrd, "somlttml", clmxl//"Depth"              , "m"   ,   &  ! hmlp
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )

      ! Temperature
      CALL histdef( nidtrd, "somltemp", clmxl//"Temperature"        , "C"   ,   &  ! ???
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      ! Temperature trends
      CALL histdef( nidtrd, "somlttto", clmxl//"T Total"             , "C/s",   &  ! total
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zout, zout )
      CALL histdef( nidtrd, "somlttax", clmxl//"T Zonal Advection", "C/s",       & ! i-adv.
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlttay", clmxl//"T Meridional Advection", "C/s",   & ! j-adv.
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlttaz", clmxl//"T Vertical Advection", "C/s",   & ! vert. adv.
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlttdh", clmxl//"T Horizontal Diffusion ", "C/s",   & ! hor. lateral diff.
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlttfo", clmxl//"T Forcing", "C/s",   & ! forcing
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlbtdz", clmxl//"T Vertical Diffusion", "C/s",   & ! vert. diff.
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlbtdt", clmxl//"T dh/dt Entrainment (Residual)", "C/s",   & ! T * dh/dt 
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zout, zout )
      IF( ln_traldf_iso ) THEN
      CALL histdef( nidtrd, "somlbtdv", clmxl//"T Vert. lateral Diffusion","C/s",   & ! vertical diffusion entrainment (ISO)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      ENDIF
#if defined key_traldf_eiv
      CALL histdef( nidtrd, "somlgtax", clmxl//"T Zonal EIV Advection", "C/s",   & ! i-adv. (eiv)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlgtay", clmxl//"T Meridional EIV Advection", "C/s",   & ! j-adv. (eiv)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlgtaz", clmxl//"T Vertical EIV Advection", "C/s",   & ! vert. adv. (eiv)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlgtat", clmxl//"T Total EIV Advection", "C/s",   & ! total advection (eiv)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif
      ! Salinity
      CALL histdef( nidtrd, "somlsalt", clmxl//"Salinity", "PSU",   & ! Mixed-layer salinity
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      ! Salinity trends
      CALL histdef( nidtrd, "somltsto", clmxl//"S Total", "PSU/s",   & ! total 
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somltsax", clmxl//"S Zonal Advection", "PSU/s",   & ! i-advection
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somltsay", clmxl//"S Meridional Advection", "PSU/s",   & ! j-advection
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somltsaz", clmxl//"S Vertical Advection", "PSU/s",   & ! vertical advection
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somltsdh", clmxl//"S Horizontal Diffusion ", "PSU/s",   & ! hor. lat. diff.
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somltsfo", clmxl//"S Forcing", "PSU/s",   & ! forcing
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )

      CALL histdef( nidtrd, "somlbsdz", clmxl//"S Vertical Diffusion", "PSU/s",   & ! vert. diff.
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlbsdt", clmxl//"S dh/dt Entrainment (Residual)", "PSU/s",   & ! S * dh/dt 
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      IF( ln_traldf_iso ) THEN
      ! vertical diffusion entrainment (ISO)
      CALL histdef( nidtrd, "somlbsdv", clmxl//"S Vertical lateral Diffusion", "PSU/s",   & ! vert. lat. diff.
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      ENDIF
#if defined key_traldf_eiv
      CALL histdef( nidtrd, "somlgsax", clmxl//"S Zonal EIV Advection", "PSU/s",   & ! i-advection (eiv)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlgsay", clmxl//"S Meridional EIV Advection", "PSU/s",   & ! j-advection (eiv)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlgsaz", clmxl//"S Vertical EIV Advection", "PSU/s",   & ! vert. adv. (eiv)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( nidtrd, "somlgsat", clmxl//"S Total EIV Advection", "PSU/s",   & ! total adv. (eiv)
         &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif
      CALL histend( nidtrd )

      IF( idebug /= 0 ) THEN
         if(lwp) WRITE(numout,*) ' debuging trd_mld_init: II. done'  
         CALL FLUSH(numout)
      ENDIF


      END SUBROUTINE trd_mld_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_trdmld = .FALSE.   !: momentum trend flag
CONTAINS
   SUBROUTINE trd_mld( kt )             ! Empty routine
      INTEGER, INTENT( in) ::   kt
!      WRITE(*,*) 'trd_mld: You should not have seen this print! error?', kt
   END SUBROUTINE trd_mld
   SUBROUTINE trd_mld_zint( pttrdmld, pstrdmld, ktrd, ctype )
      REAL, DIMENSION(:,:,:), INTENT( in ) ::   &
         pttrdmld, pstrdmld                   ! Temperature and Salinity trends
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
      CHARACTER(len=2), INTENT( in ) ::   &  
         ctype                                ! surface/bottom (2D arrays) or
         !                                    ! interior (3D arrays) physics
!      WRITE(*,*) 'trd_mld_zint: You should not have seen this print! error?', pttrdmld(1,1,1)
!      WRITE(*,*) '  "      "  : You should not have seen this print! error?', pstrdmld(1,1,1)
!      WRITE(*,*) '  "      "  : You should not have seen this print! error?', ctype
!      WRITE(*,*) '  "      "  : You should not have seen this print! error?', ktrd
   END SUBROUTINE trd_mld_zint
   SUBROUTINE trd_mld_init              ! Empty routine
!      WRITE(*,*) 'trd_mld_init: You should not have seen this print! error?'
   END SUBROUTINE trd_mld_init
#endif

   !!======================================================================
END MODULE trdmld
