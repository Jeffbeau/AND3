MODULE limthd
   !!======================================================================
   !!                  ***  MODULE limthd   ***
   !!              LIM thermo ice model : ice thermodynamic
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_thd      : thermodynamic of sea ice
   !!   lim_thd_init : initialisation of sea-ice thermodynamic
   !!----------------------------------------------------------------------
   !! * Modules used
   USE phycst          ! physical constants
   USE dom_oce         ! ocean space and time domain variables
   USE lbclnk
   USE in_out_manager  ! I/O manager
   USE ice             ! LIM sea-ice variables
   USE ice_oce         ! sea-ice/ocean variables
   USE flx_oce         ! sea-ice/ocean forcings variables 
   USE thd_ice         ! LIM thermodynamic sea-ice variables
   USE dom_ice         ! LIM sea-ice domain
   USE iceini
   USE limthd_zdf
   USE limthd_lac
   USE limtab
   USE prtctl          ! Print control
      
   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC lim_thd       ! called by lim_step

   !! * Module variables
   REAL(wp)  ::            &  ! constant values
      epsi20 = 1.e-20   ,  &
      epsi16 = 1.e-16   ,  &
      epsi04 = 1.e-04   ,  &
      zzero  = 0.e0     ,  &
      zone   = 1.e0

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!-------- -------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limthd.F90,v 1.8 2006/03/21 08:44:32 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_thd
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE lim_thd  ***       
      !!  
      !! ** Purpose : This routine manages the ice thermodynamic.
      !!         
      !! ** Action : - Initialisation of some variables
      !!             - Some preliminary computation (oceanic heat flux
      !!               at the ice base, snow acc.,heat budget of the leads)
      !!             - selection of the icy points and put them in an array
      !!             - call lim_vert_ther for vert ice thermodynamic
      !!             - back to the geographic grid
      !!             - selection of points for lateral accretion
      !!             - call lim_lat_acc  for the ice accretion
      !!             - back to the geographic grid
      !!
      !! ** References :
      !!       H. Goosse et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90
      !!
      !! History :
      !!   1.0  !  00-01 (LIM)
      !!   2.0  !  02-07 (C. Ethe, G. Madec) F90
      !!---------------------------------------------------------------------
      !! * Local variables
      INTEGER  ::   ji, jj,    &   ! dummy loop indices
         nbpb  ,               &   ! nb of icy pts for thermo. cal.
         nbpac                     ! nb of pts for lateral accretion 
      CHARACTER (len=22) :: charout
      REAL(wp) ::  &
         zfric_umin = 5e-03 ,  &   ! lower bound for the friction velocity
         zfric_umax = 2e-02        ! upper bound for the friction velocity
      REAL(wp) ::   &
         zinda              ,  &   ! switch for test. the val. of concen.
         zindb, zindg       ,  &   ! switches for test. the val of arg
         za , zh, zthsnice  ,  &
         zfric_u            ,  &   ! friction velocity 
         zfnsol             ,  &   ! total non solar heat
         zfontn             ,  &   ! heat flux from snow thickness
         zfntlat, zpareff          ! test. the val. of lead heat budget
      REAL(wp), DIMENSION(jpi,jpj) :: &
         zhicifp            ,  &   ! ice thickness for outputs
         zqlbsbq                   ! link with lead energy budget qldif
      REAL(wp), DIMENSION(jpi,jpj,2) :: &
         zmsk                      ! working array
      !!-------------------------------------------------------------------

      IF( numit == nstart  )   CALL lim_thd_init  ! Initialization (first time-step only)
   
      !-------------------------------------------!
      !   Initilization of diagnostic variables   !
      !-------------------------------------------!
      
!i est-ce utile?  oui au moins en partie
      rdvosif(:,:) = 0.e0   ! variation of ice volume at surface
      rdvobif(:,:) = 0.e0   ! variation of ice volume at bottom
      fdvolif(:,:) = 0.e0   ! total variation of ice volume
      rdvonif(:,:) = 0.e0   ! lateral variation of ice volume
      fstric (:,:) = 0.e0   ! part of solar radiation absorbing inside the ice
      fscmbq (:,:) = 0.e0   ! linked with fstric
      ffltbif(:,:) = 0.e0   ! linked with fstric
      qfvbq  (:,:) = 0.e0   ! linked with fstric
      rdmsnif(:,:) = 0.e0   ! variation of snow mass per unit area
      rdmicif(:,:) = 0.e0   ! variation of ice mass per unit area
      hicifp (:,:) = 0.e0   ! daily thermodynamic ice production. 
      zmsk (:,:,:) = 0.e0

      DO jj = 1, jpj
         DO ji = 1, jpi
            hsnif(ji,jj)  = hsnif(ji,jj) *  MAX( zzero, SIGN( zone , hsnif(ji,jj) - epsi04 ) )
         END DO
      END DO

      IF(ln_ctl)   CALL prt_ctl(tab2d_1=hsnif     , clinfo1=' lim_thd: hsnif   : ')
      
      !-----------------------------------!
      !   Treatment of particular cases   !
      !-----------------------------------!
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            !  snow is transformed into ice if the original ice cover disappears.
            zindg         = tms(ji,jj) *  MAX( zzero , SIGN( zone , -hicif(ji,jj) ) )
            hicif(ji,jj)  = hicif(ji,jj) + zindg * rhosn * hsnif(ji,jj) / rau0
            hsnif(ji,jj)  = ( zone - zindg ) * hsnif(ji,jj) + zindg * hicif(ji,jj) * ( rau0 - rhoic ) / rhosn
            dmgwi(ji,jj)  = zindg * (1.0 - frld(ji,jj)) * rhoic * hicif(ji,jj)   ! snow/ice mass
            
            !  the lead fraction, frld, must be little than or equal to amax (ice ridging).
            zthsnice      = hsnif(ji,jj) + hicif(ji,jj)
            zindb         = tms(ji,jj) * ( 1.0 - MAX( zzero , SIGN( zone , - zthsnice ) ) ) 
            za            = zindb * MIN( zone, ( 1.0 - frld(ji,jj) ) * uscomi )
            hsnif (ji,jj) = hsnif(ji,jj)  * za
            hicif (ji,jj) = hicif(ji,jj)  * za
            qstoif(ji,jj) = qstoif(ji,jj) * za
            frld  (ji,jj) = 1.0 - zindb * ( 1.0 - frld(ji,jj) ) / MAX( za , epsi20 )
            
            !  the in situ ice thickness, hicif, must be equal to or greater than hiclim.
            zh            = MAX( zone , zindb * hiclim  / MAX( hicif(ji,jj) , epsi20 ) )
            hsnif (ji,jj) = hsnif(ji,jj)  * zh
            hicif (ji,jj) = hicif(ji,jj)  * zh
            qstoif(ji,jj) = qstoif(ji,jj) * zh
            frld  (ji,jj) = ( frld(ji,jj) + ( zh - 1.0 ) ) / zh
         END DO
      END DO

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=hicif  , clinfo1=' lim_thd: hicif   : ')
         CALL prt_ctl(tab2d_1=hsnif  , clinfo1=' lim_thd: hsnif   : ')
         CALL prt_ctl(tab2d_1=dmgwi  , clinfo1=' lim_thd: dmgwi   : ')
         CALL prt_ctl(tab2d_1=qstoif , clinfo1=' lim_thd: qstoif  : ')
         CALL prt_ctl(tab2d_1=frld   , clinfo1=' lim_thd: frld    : ')
      ENDIF

      
      !-------------------------------!
      !   Thermodynamics of sea ice   !
      !-------------------------------!
      
      !      Partial computation of forcing for the thermodynamic sea ice model.
      !--------------------------------------------------------------------------

      !CDIR NOVERRCHK
      DO jj = 1, jpj
         !CDIR NOVERRCHK
         DO ji = 1, jpi
            zthsnice       = hsnif(ji,jj) + hicif(ji,jj)
            zindb          = tms(ji,jj) * ( 1.0 - MAX( zzero , SIGN( zone , - zthsnice ) ) ) 
            pfrld(ji,jj)   = frld(ji,jj)
            zinda          = 1.0 - MAX( zzero , SIGN( zone , - ( 1.0 - pfrld(ji,jj) ) ) )
            
            !  solar irradiance transmission at the mixed layer bottom and used in the lead heat budget
            thcm(ji,jj)    = 0.e0 
            
            !  net downward heat flux from the ice to the ocean, expressed as a function of ocean 
            !  temperature and turbulent mixing (McPhee, 1992)
            zfric_u        = MAX ( MIN( SQRT( ust2s(ji,jj) ) , zfric_umax ) , zfric_umin )  ! friction velocity
            fdtcn(ji,jj)  = zindb * rau0 * rcp * 0.006  * zfric_u * ( sst_io(ji,jj) - tfu(ji,jj) ) 
            qdtcn(ji,jj)  = zindb * fdtcn(ji,jj) * frld(ji,jj) * rdt_ice
                        
            !  partial computation of the lead energy budget (qldif)
            zfontn         = ( sprecip(ji,jj) / rhosn ) * xlsn  !   energy for melting
            zfnsol         = qnsr_oce(ji,jj)  !  total non solar flux
            qldif(ji,jj)   = tms(ji,jj) * ( qsr_oce(ji,jj) * ( 1.0 - thcm(ji,jj) )   &
               &                               + zfnsol + fdtcn(ji,jj) - zfontn     &
               &                               + ( 1.0 - zindb ) * fsbbq(ji,jj) )   &
               &                               * frld(ji,jj) * rdt_ice    
            !  parlat : percentage of energy used for lateral ablation (0.0) 
            zfntlat        = 1.0 - MAX( zzero , SIGN( zone ,  - qldif(ji,jj) ) )
            zpareff        = 1.0 + ( parlat - 1.0 ) * zinda * zfntlat
            zqlbsbq(ji,jj) = qldif(ji,jj) * ( 1.0 - zpareff ) / MAX( (1.0 - frld(ji,jj)) * rdt_ice , epsi16 )
            qldif  (ji,jj) = zpareff *  qldif(ji,jj)
            qdtcn  (ji,jj) = zpareff * qdtcn(ji,jj)
            
            !  energy needed to bring ocean surface layer until its freezing
            qcmif  (ji,jj) =  rau0 * rcp * fse3t(ji,jj,1) * ( tfu(ji,jj) - sst_io(ji,jj) ) * ( 1 - zinda )
            
            !  calculate oceanic heat flux.
            fbif   (ji,jj) = zindb * (  fsbbq(ji,jj) / MAX( (1.0 - frld(ji,jj)) , epsi20 ) + fdtcn(ji,jj) )
            
            ! computation of the daily thermodynamic ice production (only needed for output)
            zhicifp(ji,jj) = hicif(ji,jj) * ( 1.0 - frld(ji,jj) )
         END DO
      END DO
      
      
      !         Select icy points and fulfill arrays for the vectorial grid.
      !----------------------------------------------------------------------
      nbpb = 0
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF ( frld(ji,jj) < 1.0 ) THEN     
               nbpb      = nbpb + 1
               npb(nbpb) = (jj - 1) * jpi + ji
            ENDIF
         END DO
      END DO

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=pfrld, clinfo1=' lim_thd: pfrld   : ', tab2d_2=thcm   , clinfo2='  thcm    : ')
         CALL prt_ctl(tab2d_1=fdtcn, clinfo1=' lim_thd: fdtcn   : ', tab2d_2=qdtcn  , clinfo2='  qdtcn   : ')
         CALL prt_ctl(tab2d_1=qldif, clinfo1=' lim_thd: qldif   : ', tab2d_2=zqlbsbq, clinfo2='  zqlbsbq : ')
         CALL prt_ctl(tab2d_1=qcmif, clinfo1=' lim_thd: qcmif   : ', tab2d_2=fbif   , clinfo2='  fbif    : ')
         zmsk(:,:,1) = tms(:,:)
         CALL prt_ctl(tab2d_1=qcmif  , clinfo1=' lim_thd: qcmif   : ', mask1=zmsk)
         CALL prt_ctl(tab2d_1=zhicifp, clinfo1=' lim_thd: zhicifp : ')
         WRITE(charout, FMT="('lim_thd: nbpb = ',I4)") nbpb
         CALL prt_ctl_info(charout)
      ENDIF
      
      
      ! If there is no ice, do nothing. Otherwise, compute Top and Bottom accretion/ablation 
      !------------------------------------------------------------------------------------ 

      IF ( nbpb > 0) THEN
         
         !  put the variable in a 1-D array for thermodynamics process
         CALL tab_2d_1d( nbpb, frld_1d    (1:nbpb)     , frld       , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, h_ice_1d   (1:nbpb)     , hicif      , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, h_snow_1d  (1:nbpb)     , hsnif      , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, sist_1d    (1:nbpb)     , sist       , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 1 ), tbif(:,:,1), jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 2 ), tbif(:,:,2), jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 3 ), tbif(:,:,3), jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, qsr_ice_1d (1:nbpb)     , qsr_ice    , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, fr1_i0_1d  (1:nbpb)     , fr1_i0     , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, fr2_i0_1d  (1:nbpb)     , fr2_i0     , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, qnsr_ice_1d(1:nbpb)     , qnsr_ice   , jpi, jpj, npb(1:nbpb) )
#if ! defined key_coupled
         CALL tab_2d_1d( nbpb, qla_ice_1d (1:nbpb)     , qla_ice    , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, dqla_ice_1d(1:nbpb)     , dqla_ice   , jpi, jpj, npb(1:nbpb) )
#endif
         CALL tab_2d_1d( nbpb, dqns_ice_1d(1:nbpb)     , dqns_ice   , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, tfu_1d     (1:nbpb)     , tfu        , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, sprecip_1d (1:nbpb)     , sprecip    , jpi, jpj, npb(1:nbpb) ) 
         CALL tab_2d_1d( nbpb, fbif_1d    (1:nbpb)     , fbif       , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, thcm_1d    (1:nbpb)     , thcm       , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, qldif_1d   (1:nbpb)     , qldif      , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, qstbif_1d  (1:nbpb)     , qstoif     , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, rdmicif_1d (1:nbpb)     , rdmicif    , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, dmgwi_1d   (1:nbpb)     , dmgwi      , jpi, jpj, npb(1:nbpb) )
         CALL tab_2d_1d( nbpb, qlbbq_1d   (1:nbpb)     , zqlbsbq    , jpi, jpj, npb(1:nbpb) )
 
         CALL lim_thd_zdf( 1, nbpb )       !  compute ice growth
         
         !  back to the geographic grid.
         CALL tab_1d_2d( nbpb, frld       , npb, frld_1d   (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, hicif      , npb, h_ice_1d  (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, hsnif      , npb, h_snow_1d (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, sist       , npb, sist_1d   (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, tbif(:,:,1), npb, tbif_1d   (1:nbpb , 1 ), jpi, jpj )   
         CALL tab_1d_2d( nbpb, tbif(:,:,2), npb, tbif_1d   (1:nbpb , 2 ), jpi, jpj )   
         CALL tab_1d_2d( nbpb, tbif(:,:,3), npb, tbif_1d   (1:nbpb , 3 ), jpi, jpj )   
         CALL tab_1d_2d( nbpb, fscmbq     , npb, fscbq_1d  (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, ffltbif    , npb, fltbif_1d (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, fstric     , npb, fstbif_1d (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, qldif      , npb, qldif_1d  (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, qfvbq      , npb, qfvbq_1d  (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, qstoif     , npb, qstbif_1d (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, rdmicif    , npb, rdmicif_1d(1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, dmgwi      , npb, dmgwi_1d  (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, rdmsnif    , npb, rdmsnif_1d(1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, rdvosif    , npb, dvsbq_1d  (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, rdvobif    , npb, dvbbq_1d  (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, fdvolif    , npb, dvlbq_1d  (1:nbpb)     , jpi, jpj )
         CALL tab_1d_2d( nbpb, rdvonif    , npb, dvnbq_1d  (1:nbpb)     , jpi, jpj ) 

 
      ENDIF

      
      !      Up-date sea ice thickness.
      !---------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            phicif(ji,jj) = hicif(ji,jj)  
            hicif(ji,jj)  = hicif(ji,jj) *  ( zone -  MAX( zzero, SIGN( zone, - ( 1.0 - frld(ji,jj) ) ) ) )
         END DO
      END DO

      
      !      Tricky trick : add 2 to frld in the Southern Hemisphere.
      !----------------------------------------------------------
      IF( fcor(1,1) < 0.e0 ) THEN
         DO jj = 1, njeqm1
            DO ji = 1, jpi
               frld(ji,jj) = frld(ji,jj) + 2.0
            END DO
         END DO
      ENDIF
      
      
      !     Select points for lateral accretion (this occurs when heat exchange
      !     between ice and ocean is negative; ocean losing heat) 
      !-----------------------------------------------------------------
      nbpac = 0
      DO jj = 1, jpj
         DO ji = 1, jpi
!i yes!     IF ( ( qcmif(ji,jj) - qldif(ji,jj) ) > 0.e0 ) THEN
            IF ( tms(ji,jj) * ( qcmif(ji,jj) - qldif(ji,jj) ) > 0.e0 ) THEN
               nbpac = nbpac + 1
               npac( nbpac ) = (jj - 1) * jpi + ji
            ENDIF
         END DO
      END DO
      
      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=phicif, clinfo1=' lim_thd: phicif  : ', tab2d_2=hicif, clinfo2=' hicif : ')
         WRITE(charout, FMT="('lim_thd: nbpac = ',I4)") nbpac
         CALL prt_ctl_info(charout)
      ENDIF

      
      !
      !     If ocean gains heat do nothing ; otherwise, one performs lateral accretion
      !--------------------------------------------------------------------------------

      IF( nbpac > 0 ) THEN
         
         !...Put the variable in a 1-D array for lateral accretion
         CALL tab_2d_1d( nbpac, frld_1d   (1:nbpac)     , frld       , jpi, jpj, npac(1:nbpac) )
         CALL tab_2d_1d( nbpac, h_snow_1d (1:nbpac)     , hsnif      , jpi, jpj, npac(1:nbpac) )
         CALL tab_2d_1d( nbpac, h_ice_1d  (1:nbpac)     , hicif      , jpi, jpj, npac(1:nbpac) )
         CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 1 ), tbif(:,:,1), jpi, jpj, npac(1:nbpac) )   
         CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 2 ), tbif(:,:,2), jpi, jpj, npac(1:nbpac) )   
         CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 3 ), tbif(:,:,3), jpi, jpj, npac(1:nbpac) )   
         CALL tab_2d_1d( nbpac, qldif_1d  (1:nbpac)     , qldif      , jpi, jpj, npac(1:nbpac) )
         CALL tab_2d_1d( nbpac, qcmif_1d  (1:nbpac)     , qcmif      , jpi, jpj, npac(1:nbpac) )
         CALL tab_2d_1d( nbpac, qstbif_1d (1:nbpac)     , qstoif     , jpi, jpj, npac(1:nbpac) )
         CALL tab_2d_1d( nbpac, rdmicif_1d(1:nbpac)     , rdmicif    , jpi, jpj, npac(1:nbpac) )
         CALL tab_2d_1d( nbpac, dvlbq_1d  (1:nbpac)     , fdvolif    , jpi, jpj, npac(1:nbpac) )
         CALL tab_2d_1d( nbpac, tfu_1d    (1:nbpac)     , tfu        , jpi, jpj, npac(1:nbpac) )
        
         !  call lateral accretion routine.
         CALL lim_thd_lac( 1 , nbpac )
         
         !   back to the geographic grid
         CALL tab_1d_2d( nbpac, frld       , npac(1:nbpac), frld_1d   (1:nbpac)     , jpi, jpj )
         CALL tab_1d_2d( nbpac, hsnif      , npac(1:nbpac), h_snow_1d (1:nbpac)     , jpi, jpj )
         CALL tab_1d_2d( nbpac, hicif      , npac(1:nbpac), h_ice_1d  (1:nbpac)     , jpi, jpj )
         CALL tab_1d_2d( nbpac, tbif(:,:,1), npac(1:nbpac), tbif_1d   (1:nbpac , 1 ), jpi, jpj )
         CALL tab_1d_2d( nbpac, tbif(:,:,2), npac(1:nbpac), tbif_1d   (1:nbpac , 2 ), jpi, jpj )
         CALL tab_1d_2d( nbpac, tbif(:,:,3), npac(1:nbpac), tbif_1d   (1:nbpac , 3 ), jpi, jpj )
         CALL tab_1d_2d( nbpac, qstoif     , npac(1:nbpac), qstbif_1d (1:nbpac)     , jpi, jpj )
         CALL tab_1d_2d( nbpac, rdmicif    , npac(1:nbpac), rdmicif_1d(1:nbpac)     , jpi, jpj )
         CALL tab_1d_2d( nbpac, fdvolif    , npac(1:nbpac), dvlbq_1d  (1:nbpac)     , jpi, jpj )
        
      ENDIF
       
       
      !      Recover frld values between 0 and 1 in the Southern Hemisphere (tricky trick)
      !      Update daily thermodynamic ice production.    
      !------------------------------------------------------------------------------
       
      DO jj = 1, jpj
         DO ji = 1, jpi
            frld  (ji,jj) = MIN( frld(ji,jj), ABS( frld(ji,jj) - 2.0 ) )
            hicifp(ji,jj) =  hicif(ji,jj) * ( 1.0 - frld(ji,jj) ) - zhicifp(ji,jj) + hicifp(ji,jj)
         END DO
      END DO

      IF(ln_ctl) THEN
         CALL prt_ctl_info(' lim_thd  end  ')
         CALL prt_ctl(tab2d_1=hicif , clinfo1=' lim_thd: hicif   : ', tab2d_2=hsnif , clinfo2=' hsnif  : ')
         CALL prt_ctl(tab2d_1=frld  , clinfo1=' lim_thd: frld    : ', tab2d_2=hicifp, clinfo2=' hicifp : ')
         CALL prt_ctl(tab2d_1=phicif, clinfo1=' lim_thd: phicif  : ', tab2d_2=pfrld , clinfo2=' pfrld  : ')
         CALL prt_ctl(tab2d_1=sist  , clinfo1=' lim_thd: sist    : ')
         CALL prt_ctl(tab2d_1=tbif(:,:,1), clinfo1=' lim_thd: tbif 1  : ')
         CALL prt_ctl(tab2d_1=tbif(:,:,2), clinfo1=' lim_thd: tbif 2  : ')
         CALL prt_ctl(tab2d_1=tbif(:,:,3), clinfo1=' lim_thd: tbif 3  : ')
         CALL prt_ctl(tab2d_1=fdtcn , clinfo1=' lim_thd: fdtcn   : ', tab2d_2=qdtcn , clinfo2=' qdtcn  : ')
         CALL prt_ctl(tab2d_1=qstoif, clinfo1=' lim_thd: qstoif  : ', tab2d_2=fsbbq , clinfo2=' fsbbq  : ')
      ENDIF

    END SUBROUTINE lim_thd


    SUBROUTINE lim_thd_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_thd_init *** 
      !!                 
      !! ** Purpose :   Physical constants and parameters linked to the ice 
      !!      thermodynamics
      !!
      !! ** Method  :   Read the namicethd namelist and check the ice-thermo
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicether
      !!
      !! history :
      !!  8.5  ! 03-08 (C. Ethe) original code
      !!-------------------------------------------------------------------
      NAMELIST/namicethd/ hmelt , hiccrit, hicmin, hiclim, amax  ,        &
         &                swiqst, sbeta  , parlat, hakspl, hibspl, exld,  &
         &                hakdif, hnzst  , thth  , parsub, alphs
      !!-------------------------------------------------------------------
      

      ! Define the initial parameters
      ! -------------------------
      REWIND( numnam_ice )
      READ  ( numnam_ice , namicethd )
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*)'lim_thd_init : ice parameters for ice thermodynamic computation '
         WRITE(numout,*)'~~~~~~~~~~~~'
         WRITE(numout,*)'       maximum melting at the bottom                           hmelt        = ', hmelt
         WRITE(numout,*)'       ice thick. for lateral accretion in NH (SH)             hiccrit(1/2) = ', hiccrit
         WRITE(numout,*)'       ice thick. corr. to max. energy stored in brine pocket  hicmin       = ', hicmin  
         WRITE(numout,*)'       minimum ice thickness                                   hiclim       = ', hiclim 
         WRITE(numout,*)'       maximum lead fraction                                   amax         = ', amax 
         WRITE(numout,*)'       energy stored in brine pocket (=1) or not (=0)	  swiqst       = ', swiqst 
         WRITE(numout,*)'       numerical carac. of the scheme for diffusion in ice '
         WRITE(numout,*)'       Cranck-Nicholson (=0.5), implicit (=1), explicit (=0)   sbeta        = ', sbeta
         WRITE(numout,*)'       percentage of energy used for lateral ablation          parlat       = ', parlat
         WRITE(numout,*)'       slope of distr. for Hakkinen-Mellor lateral melting     hakspl       = ', hakspl  
         WRITE(numout,*)'       slope of distribution for Hibler lateral melting        hibspl       = ', hibspl
         WRITE(numout,*)'       exponent for leads-closure rate                         exld         = ', exld
         WRITE(numout,*)'       coefficient for diffusions of ice and snow              hakdif       = ', hakdif
         WRITE(numout,*)'       threshold thick. for comp. of eq. thermal conductivity  zhth         = ', thth 
         WRITE(numout,*)'       thickness of the surf. layer in temp. computation       hnzst        = ', hnzst
         WRITE(numout,*)'       switch for snow sublimation  (=1) or not (=0)           parsub       = ', parsub  
         WRITE(numout,*)'       coefficient for snow density when snow ice formation    alphs        = ', alphs
      ENDIF
            
      uscomi = 1.0 / ( 1.0 - amax )   ! inverse of minimum lead fraction
      rcdsn = hakdif * rcdsn 
      rcdic = hakdif * rcdic
      
      IF ( ( hsndif > 100.e0 ) .OR. ( hicdif > 100.e0 ) ) THEN
         cnscg = 0.e0
      ELSE
         cnscg = rcpsn / rcpic   ! ratio  rcpsn/rcpic
      ENDIF
 
   END SUBROUTINE lim_thd_init

#else
   !!----------------------------------------------------------------------
   !!   Default option          Dummy module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_thd         ! Dummy routine
   END SUBROUTINE lim_thd
#endif

   !!======================================================================
END MODULE limthd
