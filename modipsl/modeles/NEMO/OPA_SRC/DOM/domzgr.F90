MODULE domzgr
   !!==============================================================================
   !!                       ***  MODULE domzgr   ***
   !! Ocean initialization : domain initialization
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dom_zgr     : defined the ocean vertical coordinate system
   !!       zgr_bat      : bathymetry fields (levels and meters)
   !!       zgr_bat_zoom : modify the bathymetry field if zoom domain
   !!       zgr_bat_ctl  : check the bathymetry files
   !!       zgr_z        : reference z-coordinate 
   !!       zgr_zps      : z-coordinate with partial steps
   !!       zgr_s        : s-coordinate
   !!---------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE closea
   USE solisl
   USE ini1d           ! initialization of the 1D configuration

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dom_zgr        ! called by dom_init.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr.F90,v 1.14 2006/04/28 12:24:19 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS       

   SUBROUTINE dom_zgr
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_zgr  ***
      !!                   
      !! ** Purpose :  set the depth of model levels and the resulting 
      !!      vertical scale factors.
      !!
      !! ** Method  reference vertical coordinate
      !!        Z-coordinates : The depth of model levels is defined
      !!      from an analytical function the derivative of which gives
      !!      the vertical scale factors.
      !!      both depth and scale factors only depend on k (1d arrays).
      !!              w-level: gdepw  = fsdep(k)
      !!                       e3w(k) = dk(fsdep)(k)     = fse3(k)
      !!              t-level: gdept  = fsdep(k+0.5)
      !!                       e3t(k) = dk(fsdep)(k+0.5) = fse3(k+0.5)
      !!
      !! ** Action : - gdept, gdepw : depth of T- and W-point (m)
      !!             -  e3t, e3w    : scale factors at T- and W-levels (m)
      !!
      !! Reference :
      !!      Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!
      !! History :
      !!   9.0  !  03-08  (G. Madec)  original code
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio = 0      ! temporary integer
      !!----------------------------------------------------------------------

      ! Check Vertical coordinate options
      ! ---------------------------------
      ioptio = 0
      IF( lk_sco ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_zgr : s-coordinate'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         ioptio = ioptio + 1
      ENDIF
      IF( lk_zps ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_zgr : z-coordinate with partial steps'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         ioptio = ioptio + 1
      ENDIF
      IF( ioptio == 0 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_zgr : z-coordinate'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      IF ( ioptio > 1 ) THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) ' several vertical coordinate options used'
          nstop = nstop + 1
      ENDIF

      ! Build the vertical coordinate system
      ! ------------------------------------

      CALL zgr_z                       ! Reference z-coordinate system

      CALL zgr_bat                     ! Bathymetry fields (levels and meters)

      CALL zgr_zps                     ! Partial step z-coordinate

      CALL zgr_s                       ! s-coordinate

   END SUBROUTINE dom_zgr


   SUBROUTINE zgr_z
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!                   
      !! ** Purpose :   set the depth of model levels and the resulting 
      !!      vertical scale factors.
      !!
      !! ** Method  :   z-coordinate system (use in all type of coordinate)
      !!        The depth of model levels is defined from an analytical
      !!      function the derivative of which gives the scale factors.
      !!        both depth and scale factors only depend on k (1d arrays).
      !!              w-level: gdepw  = fsdep(k)
      !!                       e3w(k) = dk(fsdep)(k)     = fse3(k)
      !!              t-level: gdept  = fsdep(k+0.5)
      !!                       e3t(k) = dk(fsdep)(k+0.5) = fse3(k+0.5)
      !!
      !! ** Action  : - gdept, gdepw : depth of T- and W-point (m)
      !!              -  e3t, e3w    : scale factors at T- and W-levels (m)
      !!
      !! Reference :
      !!      Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!
      !! History :
      !!   9.0  !  03-08  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER  ::   jk                     ! dummy loop indices
      REAL(wp) ::   zt, zw                 ! temporary scalars
      REAL(wp) ::   &
      zsur , za0, za1, zkth, zacr,      &  ! Values set from parameters in
      zdzmin, zhmax                        ! par_CONFIG_Rxx.h90
      !!----------------------------------------------------------------------

      ! Set variables from parameters
      ! ------------------------------
       zkth = ppkth       ;   zacr = ppacr
       zdzmin = ppdzmin   ;   zhmax = pphmax

      ! If ppa1 and ppa0 and ppsur are et to pp_to_be_computed
      !  za0, za1, zsur are computed from ppdzmin , pphmax, ppkth, ppacr
      !
       IF(  ppa1  == pp_to_be_computed  .AND.  &
         &  ppa0  == pp_to_be_computed  .AND.  &
         &  ppsur == pp_to_be_computed           ) THEN
         za1 = ( ppdzmin - pphmax / FLOAT(jpk-1) )          &
             / ( TANH((1-ppkth)/ppacr) - ppacr/FLOAT(jpk-1) &
             &                         *  (  LOG( COSH( (jpk - ppkth) / ppacr) )      &
             &                             - LOG( COSH( ( 1  - ppkth) / ppacr) )  )  )

         za0  = ppdzmin - za1 * TANH( (1-ppkth) / ppacr )
         zsur = - za0 - za1 * ppacr * LOG( COSH( (1-ppkth) / ppacr )  )

       ELSE
         za1 = ppa1 ;       za0 = ppa0 ;          zsur = ppsur
       ENDIF


      ! Parameter print
      ! ---------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates'
         WRITE(numout,*) '    ~~~~~~~'
         IF (  ppkth == 0. ) THEN              
              WRITE(numout,*) '            Uniform grid with ',jpk-1,' layers'
              WRITE(numout,*) '            Total depth    :', zhmax
              WRITE(numout,*) '            Layer thickness:', zhmax/(jpk-1)
         ELSE
            IF ( ppa1 == 0. .AND. ppa0 == 0. .AND. ppsur == 0. ) THEN
               WRITE(numout,*) '         zsur, za0, za1 computed from '
               WRITE(numout,*) '                 zdzmin = ', zdzmin
               WRITE(numout,*) '                 zhmax  = ', zhmax
            ENDIF
            WRITE(numout,*) '           Value of coefficients for vertical mesh:'
            WRITE(numout,*) '                 zsur = ', zsur
            WRITE(numout,*) '                 za0  = ', za0
            WRITE(numout,*) '                 za1  = ', za1
            WRITE(numout,*) '                 zkth = ', zkth
            WRITE(numout,*) '                 zacr = ', zacr
         ENDIF
      ENDIF


      ! Reference z-coordinate (depth - scale factor at T- and W-points)
      ! ======================
      IF (  ppkth == 0. ) THEN            !  uniform vertical grid       

         za1 = zhmax/FLOAT(jpk-1) 
         DO jk = 1, jpk
            zw = FLOAT( jk )
            zt = FLOAT( jk ) + 0.5
            gdepw(jk) = ( zw - 1 ) * za1
            gdept(jk) = ( zt - 1 ) * za1
            e3w  (jk) =  za1
            e3t  (jk) =  za1
         END DO

      ELSE

         DO jk = 1, jpk
            zw = FLOAT( jk )
            zt = FLOAT( jk ) + 0.5
            gdepw(jk) = ( zsur + za0 * zw + za1 * zacr * LOG( COSH( (zw-zkth)/zacr ) )  )
            gdept(jk) = ( zsur + za0 * zt + za1 * zacr * LOG( COSH( (zt-zkth)/zacr ) )  )
            e3w  (jk) =          za0      + za1        * TANH(      (zw-zkth)/zacr   )
            e3t  (jk) =          za0      + za1        * TANH(      (zt-zkth)/zacr   )
         END DO
         gdepw(1) = 0.e0   ! force first w-level to be exactly at zero

!!DBG: 2008.11.17: the params in par_SS008.h90 yield for the first 4 levels
!            level   gdept    gdepw     e3t      e3w  
!             1     3.05     0.00     6.19     6.00
!             2     9.45     6.19     6.64     6.40
!             3    16.36    12.84     7.20     6.90
!             4    23.90    20.04     7.89     7.53
!!Since the min #-levels = 2, this means that the min btm depth should be 12.84m
!!which we want to keep the same to start.
!!REM: e3t=thickness centered on T-pt; e3w=thickness centered on w-pt; 
!!e3w(1)=??? Pretend that e3w(1) ~ 2*gdept(1)
!!e.g. of change to a 3m thick top layer, while keeping gdepw(3) unchanged
!         gdepw(2) = 3.0 
!         gdept(1) = gdepw(2)/2.0
!         e3w(1) = 2.0*gdept(1)
!         e3t(1) = gdepw(2)-gdepw(1)
!         e3t(2) = gdepw(3)-gdepw(2) 
!         gdept(2) = gdepw(2) + e3t(2)/2.0 
!         e3w(2) = gdept(2)-gdept(1)
!         e3w(3) = gdept(3)-gdept(2)


      ENDIF

      ! Control and  print
      ! ==================

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '              Reference z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level   gdept    gdepw     e3t      e3w  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, gdept(jk), gdepw(jk), e3t(jk), e3w(jk), jk = 1, jpk )
      ENDIF

      DO jk = 1, jpk
         IF( e3w(jk) <= 0. .OR. e3t(jk) <= 0. ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' e3w or e3t =< 0 '
            nstop = nstop + 1
         ENDIF
         IF( gdepw(jk) < 0. .OR. gdept(jk) < 0.) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' gdepw or gdept < 0 '
            nstop = nstop + 1
         ENDIF
      END DO

   END SUBROUTINE zgr_z


   SUBROUTINE zgr_bat
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_bat  ***
      !! 
      !! ** Purpose :   set bathymetry both in levels and meters
      !!
      !! ** Method  :   read or define mbathy and bathy arrays
      !!       * level bathymetry:
      !!      The ocean basin geometry is given by a two-dimensional array,
      !!      mbathy, which is defined as follow :
      !!            mbathy(ji,jj) = 1, ..., jpk-1, the number of ocean level
      !!                              at t-point (ji,jj).
      !!                            = 0  over the continental t-point.
      !!                            = -n over the nth island t-point.
      !!      The array mbathy is checked to verified its consistency with
      !!      model option. in particular:
      !!            mbathy must have at least 1 land grid-points (mbathy<=0)
      !!                  along closed boundary.
      !!            mbathy must be cyclic IF jperio=1.
      !!            mbathy must be lower or equal to jpk-1.
      !!            isolated ocean grid points are suppressed from mbathy
      !!                  since they are only connected to remaining
      !!                  ocean through vertical diffusion.
      !!      ntopo=-1 :   rectangular channel or bassin with a bump 
      !!      ntopo= 0 :   flat rectangular channel or basin 
      !!      ntopo= 1 :   mbathy is read in 'bathy_level.nc' NetCDF file
      !!                   bathy  is read in 'bathy_meter.nc' NetCDF file
      !!      C A U T I O N : mbathy will be modified during the initializa-
      !!      tion phase to become the number of non-zero w-levels of a water
      !!      column, with a minimum value of 1.
      !!
      !! ** Action  : - mbathy: level bathymetry (in level index)
      !!              - bathy : meter bathymetry (in meters)
      !!
      !! History :
      !!   9.0  !  03-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Local declarations
      CHARACTER (len=18) ::   clname    ! temporary characters
      LOGICAL ::   llbon                ! check the existence of bathy files
      INTEGER ::   ji, jj, jl, jk       ! dummy loop indices
      INTEGER ::   inum = 11            ! temporary logical unit
      INTEGER  ::   &
         ipi, ipj, ipk,              &  ! temporary integers
         itime,                      &  !    "          "
         ii_bump, ij_bump               ! bump center position
      INTEGER, DIMENSION (1) ::   istep
      INTEGER , DIMENSION(jpidta,jpjdta) ::   &
         idta                           ! global domain integer data
      REAL(wp) ::   &
         r_bump, h_bump, h_oce,      &  ! bump characteristics 
         zi, zj, zdate0, zdt            ! temporary scalars
      REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
         zlamt, zphit,               &  ! temporary workspace (NetCDF read)
         zdta                           ! global domain scalar data
      REAL(wp), DIMENSION(jpk) ::   &
         zdept                          ! temporary workspace (NetCDF read)
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_bat : defines level and meter bathymetry'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~'

      ! ========================================
      ! global domain level and meter bathymetry (idta,zdta)
      ! ========================================
      !                                               ! =============== ! 
      IF( ntopo == 0 .OR. ntopo == -1 ) THEN          ! defined by hand !
         !                                            ! =============== !

         IF( ntopo == 0 ) THEN                        ! flat basin

            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '         bathymetry field: flat basin'

            idta(:,:) = jpkm1                            ! flat basin 
            zdta(:,:) = gdepw(jpk)

         ELSE                                         ! bump
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '         bathymetry field: flat basin with a bump'

            ii_bump = jpidta / 3 + 3       ! i-index of the bump center
            ij_bump = jpjdta / 2           ! j-index of the bump center
            r_bump  =    6              ! bump radius (index)       
            h_bump  =  240.e0           ! bump height (meters)
            h_oce   = gdepw(jpk)        ! background ocean depth (meters)
            IF(lwp) WRITE(numout,*) '            bump characteristics: '
            IF(lwp) WRITE(numout,*) '               bump center (i,j)   = ', ii_bump, ii_bump
            IF(lwp) WRITE(numout,*) '               bump height         = ', h_bump , ' meters'
            IF(lwp) WRITE(numout,*) '               bump radius         = ', r_bump , ' index'
            IF(lwp) WRITE(numout,*) '            background ocean depth = ', h_oce  , ' meters'
            ! zdta :
            DO jj = 1, jpjdta
               DO ji = 1, jpidta
                  zi = FLOAT( ji - ii_bump ) / r_bump      
                  zj = FLOAT( jj - ij_bump ) / r_bump       
                  zdta(ji,jj) = h_oce - h_bump * EXP( -( zi*zi + zj*zj ) )
               END DO
            END DO
            ! idta :
            idta(:,:) = jpkm1
            DO jk = 1, jpkm1
               DO jj = 1, jpjdta
                  DO ji = 1, jpidta
                     IF( gdept(jk) < zdta(ji,jj) .AND. zdta(ji,jj) <= gdept(jk+1) )   idta(ji,jj) = jk
                  END DO
               END DO
            END DO
         ENDIF

         ! set boundary conditions (caution, idta on the global domain: use of jperio, not nperio)
         IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) THEN
            idta( :    , 1    ) = -1                ;      zdta( :    , 1    ) = -1.e0
            idta( :    ,jpjdta) =  0                ;      zdta( :    ,jpjdta) =  0.e0
         ELSEIF( jperio == 2 ) THEN
            idta( :    , 1    ) = idta( : ,  3  )   ;      zdta( :    , 1    ) = zdta( : ,  3  )
            idta( :    ,jpjdta) = 0                 ;      zdta( :    ,jpjdta) =  0.e0
            idta( 1    , :    ) = 0                 ;      zdta( 1    , :    ) =  0.e0
            idta(jpidta, :    ) = 0                 ;      zdta(jpidta, :    ) =  0.e0
         ELSE
            idta( :    , 1    ) = 0                 ;      zdta( :    , 1    ) =  0.e0
            idta( :    ,jpjdta) = 0                 ;      zdta( :    ,jpjdta) =  0.e0
            idta( 1    , :    ) = 0                 ;      zdta( 1    , :    ) =  0.e0
            idta(jpidta, :    ) = 0                 ;      zdta(jpidta, :    ) =  0.e0
         ENDIF

         !  EEL R5 configuration with east and west open boundaries.
         !  Two rows of zeroes are needed at the south and north for OBCs
         !  This is for compatibility with the rigid lid option. 
          
!         IF( cp_cfg == "eel" .AND. jp_cfg == 5 ) THEN
!            idta( : , 2      ) = 0                 ;      zdta( : , 2      ) =  0.e0
!            idta( : ,jpjdta-1) = 0                 ;      zdta( : ,jpjdta-1) =  0.e0
!         ENDIF

         !                                            ! =============== !
      ELSEIF( ntopo == 1 ) THEN                       !   read in file  !
         !                                            ! =============== !
         IF( lk_zco ) THEN
            clname = 'bathy_level.nc'                       ! Level bathymetry
#if defined key_agrif
            IF( .NOT. Agrif_Root() ) THEN
               clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
            ENDIF
#endif         
            INQUIRE( FILE=clname, EXIST=llbon )
            IF( llbon ) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '         read level bathymetry in ', clname
               IF(lwp) WRITE(numout,*)
               IF(lwp) print*
               IF(lwp) print*, '         read level bathymetry in ', clname
               itime = 1
               ipi = jpidta
               ipj = jpjdta
               ipk = 1
               zdt = rdt
               CALL flinopen( clname, 1, jpidta, 1, jpjdta, .FALSE.,   &
                              ipi, ipj, ipk, zlamt, zphit, zdept, itime, istep, zdate0, zdt, inum )
               CALL flinget( inum, 'Bathy_level', jpidta, jpjdta, 1,   &
                             itime, 1, 1, 1, jpidta, 1, jpjdta, zdta(:,:) )
               idta(:,:) = zdta(:,:)
               CALL flinclo( inum )

            ELSE
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*)'    zgr_bat : unable to read the file', clname
               nstop = nstop + 1
            ENDIF   
  
         ELSEIF( lk_zps ) THEN
            clname = 'bathy_meter.nc'                       ! meter bathymetry
#if defined key_agrif
            IF( .NOT. Agrif_Root() ) THEN
               clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
            ENDIF
#endif   	 
            INQUIRE( FILE=clname, EXIST=llbon )
            IF( llbon ) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '         read meter bathymetry in ', clname
               IF(lwp) WRITE(numout,*)
               itime = 1
               ipi = jpidta
               ipj = jpjdta
               ipk = 1
               zdt = rdt
               CALL flinopen( clname, 1, jpidta, 1, jpjdta, .FALSE.,   &    
                              ipi, ipj, ipk, zlamt, zphit, zdept, itime, istep, zdate0, zdt, inum )
               CALL flinget( inum, 'Bathymetry', jpidta, jpjdta, 1,   &
                             itime, 1, 1, 1, jpidta, 1, jpjdta, zdta(:,:) ) 
               CALL flinclo( inum )
               idta(:,:) = jpkm1      ! initialisation
            ELSE
               IF(lwp) WRITE(numout,cform_err)       
               IF(lwp) WRITE(numout,*)'    zgr_bat : unable to read the file', clname
               nstop = nstop + 1
            ENDIF
         ENDIF
         !                                            ! =============== !
      ELSE                                            !      error      !
         !                                            ! =============== !
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          parameter , ntopo = ', ntopo
         nstop = nstop + 1
      ENDIF


      ! =======================================
      ! local domain level and meter bathymetry (mbathy,bathy)
      ! =======================================

      mbathy(:,:) = 0                                 ! set to zero extra halo points
      bathy (:,:) = 0.e0                              ! (require for mpp case)

      DO jj = 1, nlcj                                 ! interior values
         DO ji = 1, nlci
            mbathy(ji,jj) = idta( mig(ji), mjg(jj) )
            bathy (ji,jj) = zdta( mig(ji), mjg(jj) )
         END DO
      END DO

      ! =======================
      ! NO closed seas or lakes
      ! =======================

      IF( nclosea == 0 ) THEN
         DO jl = 1, jpncs
            DO jj = ncsj1(jl), ncsj2(jl)
               DO ji = ncsi1(jl), ncsi2(jl)
                  mbathy(ji,jj) = 0                   ! suppress closed seas
                  bathy (ji,jj) = 0.e0                ! and lakes
               END DO
            END DO
         END DO
      ENDIF

      ! ===========
      ! Zoom domain 
      ! ===========

      IF( lzoom )   CALL zgr_bat_zoom

      ! ================
      ! Bathymetry check
      ! ================

      CALL zgr_bat_ctl

   END SUBROUTINE zgr_bat


   SUBROUTINE zgr_bat_zoom
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_bat_zoom  ***
      !!
      !! ** Purpose : - Close zoom domain boundary if necessary
      !!              - Suppress Med Sea from ORCA R2 and R05 arctic zoom
      !!
      !! ** Method  : 
      !!
      !! ** Action  : - update mbathy: level bathymetry (in level index)
      !!
      !! History :
      !!   9.0  !  03-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local variables
      INTEGER ::   ii0, ii1, ij0, ij1   ! temporary integers
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_bat_zoom : modify the level bathymetry for zoom domain'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~'

      ! Zoom domain
      ! ===========

      ! Forced closed boundary if required
      IF( lzoom_w )   mbathy( mi0(jpizoom):mi1(jpizoom) , :  ) = 0
      IF( lzoom_s )   mbathy(  : , mj0(jpjzoom):mj1(jpjzoom) ) = 0
      IF( lzoom_e )   mbathy( mi0(jpiglo+jpizoom-1):mi1(jpiglo+jpizoom-1) , :  ) = 0
      IF( lzoom_n )   mbathy(  : , mj0(jpjglo+jpjzoom-1):mj1(jpjglo+jpjzoom-1) ) = 0

!!DB: delete ORCA
      ! Configuration specific domain modifications
      ! (here, ORCA arctic configuration: suppress Med Sea)
!      IF( cp_cfg == "orca" .AND. lzoom_arct ) THEN

   END SUBROUTINE zgr_bat_zoom


   SUBROUTINE zgr_bat_ctl
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_bat_ctl  ***
      !!
      !! ** Purpose :   check the bathymetry in levels
      !!
      !! ** Method  :   The array mbathy is checked to verified its consistency
      !!      with the model options. in particular:
      !!            mbathy must have at least 1 land grid-points (mbathy<=0)
      !!                  along closed boundary.
      !!            mbathy must be cyclic IF jperio=1.
      !!            mbathy must be lower or equal to jpk-1.
      !!            isolated ocean grid points are suppressed from mbathy
      !!                  since they are only connected to remaining
      !!                  ocean through vertical diffusion.
      !!      C A U T I O N : mbathy will be modified during the initializa-
      !!      tion phase to become the number of non-zero w-levels of a water
      !!      column, with a minimum value of 1.
      !!
      !! ** Action  : - update mbathy: level bathymetry (in level index)
      !!              - update bathy : meter bathymetry (in meters)
      !!
      !! History :
      !!   9.0  !  03-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj, jl           ! dummy loop indices
      INTEGER ::   &
         icompt, ibtest, ikmax          ! temporary integers
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zbathy                         ! temporary workspace
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_bat_ctl : check the bathymetry'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'

      ! ================
      ! Bathymetry check
      ! ================

      IF( .NOT. lk_cfg_1d )   THEN

         ! Suppress isolated ocean grid points

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'                   suppress isolated ocean grid points'
         IF(lwp) WRITE(numout,*)'                   -----------------------------------'

         icompt = 0

         DO jl = 1, 2

            IF( nperio == 1 .OR. nperio  ==  4 .OR. nperio  ==  6 ) THEN
               mbathy( 1 ,:) = mbathy(jpim1,:)
               mbathy(jpi,:) = mbathy(  2  ,:)
            ENDIF
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ibtest = MAX( mbathy(ji-1,jj), mbathy(ji+1,jj),   &
                     mbathy(ji,jj-1),mbathy(ji,jj+1) )
                  IF( ibtest < mbathy(ji,jj) ) THEN
                     IF(lwp) WRITE(numout,*) ' the number of ocean level at ',   &
                        'grid-point (i,j) =  ',ji,jj,' is changed from ',   &
                        mbathy(ji,jj),' to ', ibtest
                     mbathy(ji,jj) = ibtest
                     icompt = icompt + 1
                  ENDIF
               END DO
            END DO

         END DO
         IF( icompt == 0 ) THEN
            IF(lwp) WRITE(numout,*)'     no isolated ocean grid points'
         ELSE
            IF(lwp) WRITE(numout,*)'    ',icompt,' ocean grid points suppressed'
         ENDIF
         IF( lk_mpp ) THEN
            zbathy(:,:) = FLOAT( mbathy(:,:) )
            CALL lbc_lnk( zbathy, 'T', 1. )
            mbathy(:,:) = INT( zbathy(:,:) )
         ENDIF

         ! 3.2 East-west cyclic boundary conditions

         IF( nperio == 0 ) THEN
            IF(lwp) WRITE(numout,*) ' mbathy set to 0 along east and west',   &
               ' boundary: nperio = ', nperio
            IF( lk_mpp ) THEN
               IF( nbondi == -1 .OR. nbondi == 2 ) THEN
                  IF( jperio /= 1 )   mbathy(1,:) = 0
               ENDIF
               IF( nbondi == 1 .OR. nbondi == 2 ) THEN
                  IF( jperio /= 1 )   mbathy(nlci,:) = 0
               ENDIF
            ELSE
               mbathy( 1 ,:) = 0
               mbathy(jpi,:) = 0
            ENDIF
         ELSEIF( nperio == 1 .OR. nperio == 4 .OR. nperio ==  6 ) THEN
            IF(lwp) WRITE(numout,*)' east-west cyclic boundary conditions',   &
               ' on mbathy: nperio = ', nperio
            mbathy( 1 ,:) = mbathy(jpim1,:)
            mbathy(jpi,:) = mbathy(  2  ,:)
         ELSEIF( nperio == 2 ) THEN
            IF(lwp) WRITE(numout,*) '   equatorial boundary conditions',   &
               ' on mbathy: nperio = ', nperio
         ELSE
            IF(lwp) WRITE(numout,*) '    e r r o r'
            IF(lwp) WRITE(numout,*) '    parameter , nperio = ', nperio
            !         STOP 'dom_mba'
         ENDIF

         ! Set to zero mbathy over islands if necessary  (lk_isl=F)
         IF( .NOT. lk_isl ) THEN    ! No island
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '         mbathy set to 0 over islands'
            IF(lwp) WRITE(numout,*) '         ----------------------------'

            mbathy(:,:) = MAX( 0, mbathy(:,:) )

            !  Boundary condition on mbathy
            IF( .NOT.lk_mpp ) THEN 
               
               !!bug ???  y reflechir!
               !   ... mono- or macro-tasking: T-point, >0, 2D array, no slab
               
               zbathy(:,:) = FLOAT( mbathy(:,:) )
               CALL lbc_lnk( zbathy, 'T', 1. )
               mbathy(:,:) = INT( zbathy(:,:) )
            ENDIF

         ENDIF

      ENDIF

      ! Number of ocean level inferior or equal to jpkm1

      ikmax = 0
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikmax = MAX( ikmax, mbathy(ji,jj) )
         END DO
      END DO
      !!! test a faire:   ikmax = MAX( mbathy(:,:) )   ???

      IF( ikmax > jpkm1 ) THEN
         IF(lwp) WRITE(numout,*) ' maximum number of ocean level = ', ikmax,' >  jpk-1'
         IF(lwp) WRITE(numout,*) ' change jpk to ',ikmax+1,' to use the exact ead bathymetry'
      ELSE IF( ikmax < jpkm1 ) THEN
         IF(lwp) WRITE(numout,*) ' maximum number of ocean level = ', ikmax,' < jpk-1' 
         IF(lwp) WRITE(numout,*) ' you can decrease jpk to ', ikmax+1
      ENDIF

      IF( lwp .AND. nprint == 1 ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' bathymetric field '
         WRITE(numout,*) ' ------------------'
         WRITE(numout,*) ' number of non-zero T-levels '
         CALL prihin( mbathy, jpi, jpj, 1, jpi,   &
                      1     , 1  , jpj, 1, 3  ,   &
                      numout )
         WRITE(numout,*)
      ENDIF
   END SUBROUTINE zgr_bat_ctl


#  include "domzgr_zps.h90"


#  include "domzgr_s.h90"


   !!======================================================================
END MODULE domzgr
