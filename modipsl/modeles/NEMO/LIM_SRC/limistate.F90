MODULE limistate
   !!======================================================================
   !!                     ***  MODULE  limistate  ***
   !!              Initialisation of diagnostics ice variables
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_istate      :  Initialisation of diagnostics ice variables
   !!   lim_istate_init :  initialization of ice state and namelist read
   !!----------------------------------------------------------------------
   !! * Modules used
   USE phycst
   USE ocfzpt
   USE oce             ! dynamics and tracers variables
   USE dom_oce
   USE par_ice         ! ice parameters
   USE ice_oce         ! ice variables
   USE in_out_manager
   USE dom_ice
   USE ice
   USE lbclnk

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC lim_istate      ! routine called by lim_init.F90

   !! * Module variables
   REAL(wp) ::           & !!! ** init namelist (namiceini) **
      ttest  = 2.0  ,    &  ! threshold water temperature for initial sea ice
      hninn  = 0.5  ,    &  ! initial snow thickness in the north
      hginn  = 3.0  ,    &  ! initial ice thickness in the north
      alinn  = 0.05 ,    &  ! initial leads area in the north
      hnins  = 0.1  ,    &  ! initial snow thickness in the south
      hgins  = 1.0  ,    &  ! initial ice thickness in the south
      alins  = 0.1          ! initial leads area in the south

   REAL(wp)  ::          &  ! constant values
      zzero   = 0.e0  ,  &
      zone    = 1.e0
   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limistate.F90,v 1.5 2006/03/21 08:38:38 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_istate
      !!-------------------------------------------------------------------
      !!                    ***  ROUTINE lim_istate  ***
      !!
      !! ** Purpose :   defined the sea-ice initial state
      !!
      !! ** Method  :   restart from a state defined in a binary file
      !!                or from arbitrary sea-ice conditions
      !!
      !! History :
      !!   2.0  !  01-04  (C. Ethe, G. Madec)  Original code
      !!        !  04-04  (S. Theetten) initialization from a file
      !!--------------------------------------------------------------------
      !! * Local variables
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zidto,    &  ! temporary scalar
         zs0, ztf, zbin
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         ztn
      !--------------------------------------------------------------------

 
      CALL lim_istate_init     !  reading the initials parameters of the ice

      !-- Initialisation of sst,sss,u,v do i=1,jpi
      u_io(:,:)  = 0.e0       ! ice velocity in x direction
      v_io(:,:)  = 0.e0       ! ice velocity in y direction

      IF( ln_limini ) THEN    ! 
        
         ! Initialisation at tn if no ice or sst_ini if ice
         ! Idem for salinity

      !--- Criterion for presence (zidto=1.) or absence (zidto=0.) of ice
         DO jj = 1 , jpj
            DO ji = 1 , jpi
               
               zidto = MAX(zzero, - SIGN(1.,frld(ji,jj) - 1.))
               
               sst_io(ji,jj) = ( nfice - 1 ) * (zidto * sst_ini(ji,jj)  + &   ! use the ocean initial values
                    &          (1.0 - zidto ) * ( tn(ji,jj,1) + rt0 ))        ! tricky trick *(nfice-1) !
               sss_io(ji,jj) = ( nfice - 1 ) * (zidto * sss_ini(ji,jj) + &
                    &          (1.0 - zidto ) *  sn(ji,jj,1) )

               ! to avoid the the melting of ice, several layers (mixed layer) should be
               ! set to sst_ini (sss_ini) if there is ice
               ! example for one layer 
               ! tn(ji,jj,1) = zidto * ( sst_ini(ji,jj) - rt0 )  + (1.0 - zidto ) *  tn(ji,jj,1)
               ! sn(ji,jj,1) = zidto * sss_ini(ji,jj)  + (1.0 - zidto ) *  sn(ji,jj,1)
               ! tb(ji,jj,1) = tn(ji,jj,1)
               ! sb(ji,jj,1) = sn(ji,jj,1)
            END DO
         END DO
         
         
         !  tfu: Melting point of sea water
         tfu(:,:)  = ztf   
         
         tfu(:,:)  = ABS ( rt0 - 0.0575       * sss_ini(:,:)                               &
              &                    + 1.710523e-03 * sss_ini(:,:) * SQRT( sss_ini(:,:) )    &
              &                    - 2.154996e-04 * sss_ini(:,:) * sss_ini(:,:) )
      ELSE                     !

         
         ! Initialisation at tn or -2 if ice
         DO jj = 1, jpj
            DO ji = 1, jpi
               zbin = MAX( 0., SIGN( 1., fzptn(ji,jj) - tn(ji,jj,1) ) )
               ztn(ji,jj) = ( (1.-zbin) * tn(ji,jj,1) - 2. * zbin + rt0 ) * tmask(ji,jj,1)
            END DO
         END DO
         
         u_io  (:,:) = 0.e0
         v_io  (:,:) = 0.e0
         sst_io(:,:) = ( nfice - 1 ) * ( tn(:,:,1) + rt0 )   ! use the ocean initial values
         sss_io(:,:) = ( nfice - 1 ) *   sn(:,:,1)           ! tricky trick *(nfice-1) !
         
         ! reference salinity 34psu
         zs0 = 34.e0
         ztf = ABS ( rt0 - 0.0575       * zs0                           &
              &                    + 1.710523e-03 * zs0 * SQRT( zs0 )   &
              &                    - 2.154996e-04 * zs0 *zs0          )
         
         !  tfu: Melting point of sea water
         tfu(:,:)  = ztf   
         
         DO jj = 1, jpj
            DO ji = 1, jpi
               !--- Criterion for presence (zidto=1) or absence (zidto=0) of ice
               zidto  = tms(ji,jj) * ( 1.0 - MAX(zzero, SIGN( zone, ztn(ji,jj) - tfu(ji,jj) - ttest) ) )
               
               IF( fcor(ji,jj) >= 0.e0 ) THEN     !--  Northern hemisphere.
                  hicif(ji,jj)   = zidto * hginn
                  frld(ji,jj)    = zidto * alinn + ( 1.0 - zidto ) * 1.0
                  hsnif(ji,jj)   = zidto * hninn
               ELSE                               !---  Southern hemisphere.
                  hicif(ji,jj)   = zidto * hgins
                  frld(ji,jj)    = zidto * alins + ( 1.0 - zidto ) * 1.0
                  hsnif(ji,jj)   = zidto * hnins
               ENDIF
            END DO
         END DO
         
         sist  (:,:)   = tfu(:,:)
         tbif  (:,:,1) = tfu(:,:)
         tbif  (:,:,2) = tfu(:,:)
         tbif  (:,:,3) = tfu(:,:)
      
      ENDIF
      fsbbq (:,:)   = 0.e0
      qstoif(:,:)   = 0.e0
      u_ice (:,:)   = 0.e0
      v_ice (:,:)   = 0.e0
# if defined key_coupled
      albege(:,:)   = 0.8 * tms(:,:)
# endif

      !---  Moments for advection.             

      sxice (:,:)  = 0.e0   ;   sxsn (:,:)  = 0.e0   ;   sxa  (:,:)  = 0.e0
      syice (:,:)  = 0.e0   ;   sysn (:,:)  = 0.e0   ;   sya  (:,:)  = 0.e0
      sxxice(:,:)  = 0.e0   ;   sxxsn(:,:)  = 0.e0   ;   sxxa (:,:)  = 0.e0
      syyice(:,:)  = 0.e0   ;   syysn(:,:)  = 0.e0   ;   syya (:,:)  = 0.e0
      sxyice(:,:)  = 0.e0   ;   sxysn(:,:)  = 0.e0   ;   sxya (:,:)  = 0.e0

      sxc0  (:,:)  = 0.e0   ;   sxc1 (:,:)  = 0.e0   ;   sxc2 (:,:)  = 0.e0
      syc0  (:,:)  = 0.e0   ;   syc1 (:,:)  = 0.e0   ;   syc2 (:,:)  = 0.e0
      sxxc0 (:,:)  = 0.e0   ;   sxxc1(:,:)  = 0.e0   ;   sxxc2(:,:)  = 0.e0
      syyc0 (:,:)  = 0.e0   ;   syyc1(:,:)  = 0.e0   ;   syyc2(:,:)  = 0.e0
      sxyc0 (:,:)  = 0.e0   ;   sxyc1(:,:)  = 0.e0   ;   sxyc2(:,:)  = 0.e0

      sxst  (:,:)  = 0.e0
      syst  (:,:)  = 0.e0
      sxxst (:,:)  = 0.e0
      syyst (:,:)  = 0.e0
      sxyst (:,:)  = 0.e0

      !-- lateral boundary conditions
      CALL lbc_lnk( hicif, 'T', 1. )
      CALL lbc_lnk( frld , 'T', 1. )

      ! C A U T I O N  frld = 1 over land and lbc_lnk put zero along 
      ! *************  closed boundaries herefore we force to one over land
      frld(:,:) = tms(:,:) * frld(:,:) + ( 1. - tms(:,:) )   

      CALL lbc_lnk( hsnif, 'T', 1. )
      CALL lbc_lnk( sist , 'T', 1. )
      DO jk = 1, jplayersp1
         CALL lbc_lnk(tbif(:,:,jk), 'T', 1. )
      END DO
      CALL lbc_lnk( fsbbq  , 'T', 1. )
      CALL lbc_lnk( qstoif , 'T', 1. )
      CALL lbc_lnk( sss_io , 'T', 1. )

   END SUBROUTINE lim_istate

   
   SUBROUTINE lim_istate_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_istate_init  ***
      !!        
      !! ** Purpose :   Definition of initial state of the ice 
      !!
      !! ** Method  :   Read the namiceini namelist and check the parameter 
      !!                values called at the first timestep (nit000)
      !!                or
      !!                Read 7 variables from a previous restart file
      !!                sst, sst, hicif, hsnif, frld, ts & tbif
      !!
      !! ** input   :   Namelist namiceini
      !!
      !! history
      !!  8.5  ! 03-08 (C. Ethe) original code
      !!  9.0  ! 04-04 (S. Theetten) read a file
      !!-------------------------------------------------------------------
      !! * Modules used
      USE ice
      USE ioipsl

      NAMELIST/namiceini/ ln_limini, ln_limdmp, ttest, hninn, hginn, alinn, &
         &                hnins, hgins, alins
      !!-------------------------------------------------------------------
      !! local declaration
      INTEGER, PARAMETER ::   jpmois=1
      
      INTEGER ::                   &
           itime, ipi, ipj, ipk  , & ! temporary integers
           inum_ice
      
      INTEGER ::  istep(jpmois)
      
      REAL(wp) ::   zdate0, zdt
      REAL(wp), DIMENSION(jpi,jpj) ::   zlon, zlat
      REAL(wp), DIMENSION(3) ::   zlev
      
      CHARACTER (len=32) :: cl_icedata
      
      LOGICAL :: llbon
      !!-------------------------------------------------------------------
      
      ! Read Namelist namiceini 

      REWIND ( numnam_ice )
      READ   ( numnam_ice , namiceini )
      
      IF(.NOT. ln_limini) THEN 
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'lim_istate_init : ice parameters inititialisation '
            WRITE(numout,*) '~~~~~~~~~~~~~~~'
            WRITE(numout,*) '         threshold water temp. for initial sea-ice    ttest      = ', ttest
            WRITE(numout,*) '         initial snow thickness in the north          hninn      = ', hninn
            WRITE(numout,*) '         initial ice thickness in the north           hginn      = ', hginn 
            WRITE(numout,*) '         initial leads area in the north              alinn      = ', alinn            
            WRITE(numout,*) '         initial snow thickness in the south          hnins      = ', hnins 
            WRITE(numout,*) '         initial ice thickness in the south           hgins      = ', hgins
            WRITE(numout,*) '         initial leads area in the south              alins      = ', alins
         ENDIF
      ENDIF

      IF( ln_limini ) THEN                      ! Ice initialization using input file

         cl_icedata = 'Ice_initialization.nc'
         INQUIRE( FILE=cl_icedata, EXIST=llbon )
         IF( llbon ) THEN
            IF(lwp) THEN
               WRITE(numout,*) ' '
               WRITE(numout,*) 'lim_istate_init : ice state initialization with : ',cl_icedata
               WRITE(numout,*) '~~~~~~~~~~~~~~~'
               WRITE(numout,*) '         Ice state initialization using input file    ln_limini  = ', ln_limini
               WRITE(numout,*) '         Ice damping                                  ln_limdmp  = ', ln_limdmp
               WRITE(numout,*) ' '
            ENDIF
            
            itime = 1
            ipi=jpiglo
            ipj=jpjglo
            ipk=1
            zdt=rdt
            
            CALL flinopen( TRIM(cl_icedata), mig(1), nlci, mjg(1), nlcj, .FALSE., &
               &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, inum_ice )
            
            CALL flinget( inum_ice, 'sst', jpidta, jpjdta, 1,  &
               &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sst_ini(1:nlci,1:nlcj) )
            
            CALL flinget( inum_ice, 'sss', jpidta, jpjdta, 1,  &
               &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sss_ini(1:nlci,1:nlcj) )
            
            CALL flinget( inum_ice, 'hicif', jpidta, jpjdta, 1,  &
               &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, hicif(1:nlci,1:nlcj) )
            
            CALL flinget( inum_ice, 'hsnif', jpidta, jpjdta, 1,  &
               &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, hsnif(1:nlci,1:nlcj) )
            
            CALL flinget( inum_ice, 'frld', jpidta, jpjdta, 1,  &
               &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, frld(1:nlci,1:nlcj) )
            
            CALL flinget( inum_ice, 'ts', jpidta, jpjdta, 1,  &
               &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sist(1:nlci,1:nlcj) )
            
            CALL flinclo( inum_ice)
            
            itime = 1
            ipi=jpiglo
            ipj=jpjglo
            ipk=jplayersp1
            
            CALL flinopen( TRIM(cl_icedata), mig(1), nlci, mjg(1), nlcj, .FALSE.,  &
               &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, inum_ice )
            
            CALL flinget( inum_ice, 'tbif', jpidta, jpjdta, ipk,  &
               &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, tbif(1:nlci,1:nlcj,1:ipk) )
            
            CALL flinclo( inum_ice)
            
         ELSE
            IF(lwp) WRITE(numout,cform_err) 
            IF(lwp) WRITE(numout,*) '            ',cl_icedata, ' not found !'
            nstop = nstop + 1
         ENDIF
      ENDIF

   END SUBROUTINE lim_istate_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module          NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_istate          ! Empty routine
   END SUBROUTINE lim_istate
#endif

   !!======================================================================
END MODULE limistate
