MODULE trdvor
   !!======================================================================
   !!                       ***  MODULE  trdvor  ***
   !! Ocean diagnostics:  momentum trends
   !!=====================================================================
   
#if defined key_trdvor   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_trdvor'   : momentum trend diagnostics
   !!----------------------------------------------------------------------
   !!   trd_vor      : momentum trends averaged over the depth
   !!   trd_vor_zint : vorticity vertical integration
   !!   trd_vor_init : initialization step
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE trdmod_oce      ! ocean variables trends
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE phycst          ! Define parameters for the routines
   USE ldfdyn_oce      ! ocean active tracers: lateral physics
   USE daymod          ! calandar
   USE dianam          ! build the name of file (routine)
   USE ldfslp          ! iso-neutral slopes
   USE zdfmxl          ! mixed layer depth
   USE ioipsl          ! NetCDF library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)


   IMPLICIT NONE
   PRIVATE

   !! * Interfaces
   INTERFACE trd_vor_zint
      MODULE PROCEDURE trd_vor_zint_2d, trd_vor_zint_3d
   END INTERFACE

   !! * Accessibility
   PUBLIC trd_vor        ! routine called by step.F90
   PUBLIC trd_vor_zint   ! routine called by dynamics routines
   PUBLIC trd_vor_init   ! routine called by opa.F90

   !! * Shared module variables
   LOGICAL, PUBLIC ::   lk_trdvor = .TRUE.   ! momentum trend flag

   !! * Module variables
   INTEGER ::                &
      nh_t, nmoydpvor  ,     &
      nidvor, nhoridvor,     &
      ndexvor1(jpi*jpj),     &
      ndimvor1, icount,      &
      idebug                    ! (0/1) set it to 1 in case of problem to have more print

   REAL(wp), DIMENSION(jpi,jpj) ::  &
     vor_avr    ,     &  ! average
     vor_avrb   ,     &  ! before vorticity (kt-1)
     vor_avrbb  ,     &  ! vorticity at begining of the nwrite-1 timestep averaging period
     vor_avrbn  ,     &  ! after vorticity at time step after the
     rotot      ,     &  ! begining of the NWRITE-1 timesteps
     vor_avrtot ,     &
     vor_avrres

   REAL(wp), DIMENSION(jpi,jpj,jplvor)::   &  !: curl of trends
      vortrd   

   CHARACTER(len=12) ::   cvort

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldfdyn_substitute.h90"
#  include "vectopt_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdvor.F90,v 1.5 2005/12/12 14:18:08 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
  
CONTAINS

   SUBROUTINE trd_vor_zint_2d( putrdvor, pvtrdvor, ktrd )
      !!----------------------------------------------------------------------------
      !!                  ***  ROUTINE trd_vor_zint  ***
      !!
      !! ** Purpose :   computation of vertically integrated vorticity budgets
      !!      from ocean surface down to control surface (NetCDF output)
      !!
      !! ** Method/usage :
      !!      integration done over nwrite-1 time steps
      !!
      !!
      !! ** Action :
      !!            /comvor/   :
      !!                         vor_avr          average
      !!                         vor_avrb         vorticity at kt-1
      !!                         vor_avrbb        vorticity at begining of the NWRITE-1
      !!                                          time steps averaging period
      !!                         vor_avrbn         vorticity at time step after the
      !!                                          begining of the NWRITE-1 time
      !!                                          steps averaging period
      !!
      !!                 trends :
      !!
      !!                  vortrd (,,1) = Pressure Gradient Trend
      !!                  vortrd (,,2) = KE Gradient Trend
      !!                  vortrd (,,3) = Relative Vorticity Trend
      !!                  vortrd (,,4) = Coriolis Term Trend
      !!                  vortrd (,,5) = Horizontal Diffusion Trend
      !!                  vortrd (,,6) = Vertical Advection Trend
      !!                  vortrd (,,7) = Vertical Diffusion Trend
      !!                  vortrd (,,8) = Surface Pressure Grad. Trend
      !!                  vortrd (,,9) = Beta V
      !!                  vortrd (,,10) = forcing term
      !!		  vortrd (,,11) = bottom friction term
      !!                  rotot(,) : total cumulative trends over nwrite-1 time steps
      !!                  vor_avrtot(,) : first membre of vrticity equation
      !!                  vor_avrres(,) : residual = dh/dt entrainment
      !!
      !!      trends output in netCDF format using ioipsl
      !!
      !! History :
      !!   9.0  !  04-06  (L. Brunier, A-M. Treguier) Original code 
      !!        !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index

      REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
         putrdvor,                         &  ! u vorticity trend 
         pvtrdvor                             ! v vorticity trend

      !! * Local declarations
      INTEGER ::   ji, jj
      INTEGER ::   ikbu, ikbum1, ikbv, ikbvm1
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zudpvor,                       &  ! total cmulative trends
         zvdpvor                           !   "      "        "
      !!----------------------------------------------------------------------

      ! Initialization
      zudpvor(:,:) = 0.e0
      zvdpvor(:,:) = 0.e0

      CALL lbc_lnk( putrdvor,  'U' , -1. )
      CALL lbc_lnk( pvtrdvor,  'V' , -1. )

      !  =====================================
      !  I vertical integration of 2D trends
      !  =====================================

      SELECT CASE (ktrd) 

      CASE (jpvorbfr)        ! bottom friction

         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               ikbu   = min( mbathy(ji+1,jj), mbathy(ji,jj) )
               ikbum1 = max( ikbu-1, 1 )
               ikbv   = min( mbathy(ji,jj+1), mbathy(ji,jj) )
               ikbvm1 = max( ikbv-1, 1 )
            
               zudpvor(ji,jj) = putrdvor(ji,jj) * fse3u(ji,jj,ikbum1) * e1u(ji,jj) * umask(ji,jj,ikbum1)
               zvdpvor(ji,jj) = pvtrdvor(ji,jj) * fse3v(ji,jj,ikbvm1) * e2v(ji,jj) * vmask(ji,jj,ikbvm1)
            END DO
         END DO

      CASE (jpvorswf)        ! wind stress

         zudpvor(:,:) = putrdvor(:,:) * fse3u(:,:,1) * e1u(:,:) * umask(:,:,1)
         zvdpvor(:,:) = pvtrdvor(:,:) * fse3v(:,:,1) * e2v(:,:) * vmask(:,:,1)

      END SELECT

      ! Average except for Beta.V
      zudpvor(:,:) = zudpvor(:,:) * hur(:,:)
      zvdpvor(:,:) = zvdpvor(:,:) * hvr(:,:)
   
      ! Curl
      DO ji=1,jpim1
         DO jj=1,jpjm1
            vortrd(ji,jj,ktrd) = ( zvdpvor(ji+1,jj) - zvdpvor(ji,jj)        &
                 &                - ( zudpvor(ji,jj+1) - zudpvor(ji,jj) ) ) &
                 &               / ( e1f(ji,jj) * e2f(ji,jj) )
         END DO
      END DO

      ! Surface mask
      vortrd(:,:,ktrd) = vortrd(:,:,ktrd) * fmask(:,:,1)

      IF( idebug /= 0 ) THEN
         IF(lwp) WRITE(numout,*) ' debuging trd_vor_zint: I done'
         CALL FLUSH(numout)
      ENDIF

   END SUBROUTINE trd_vor_zint_2d



   SUBROUTINE trd_vor_zint_3d( putrdvor, pvtrdvor, ktrd )
      !!----------------------------------------------------------------------------
      !!                  ***  ROUTINE trd_vor_zint  ***
      !!
      !! ** Purpose :   computation of vertically integrated vorticity budgets
      !!      from ocean surface down to control surface (NetCDF output)
      !!
      !! ** Method/usage :
      !!      integration done over nwrite-1 time steps
      !!
      !!
      !! ** Action :
      !!            /comvor/   :
      !!                         vor_avr          average
      !!                         vor_avrb         vorticity at kt-1
      !!                         vor_avrbb        vorticity at begining of the NWRITE-1
      !!                                          time steps averaging period
      !!                         vor_avrbn         vorticity at time step after the
      !!                                          begining of the NWRITE-1 time
      !!                                          steps averaging period
      !!
      !!                 trends :
      !!
      !!                  vortrd (,,1) = Pressure Gradient Trend
      !!                  vortrd (,,2) = KE Gradient Trend
      !!                  vortrd (,,3) = Relative Vorticity Trend
      !!                  vortrd (,,4) = Coriolis Term Trend
      !!                  vortrd (,,5) = Horizontal Diffusion Trend
      !!                  vortrd (,,6) = Vertical Advection Trend
      !!                  vortrd (,,7) = Vertical Diffusion Trend
      !!                  vortrd (,,8) = Surface Pressure Grad. Trend
      !!                  vortrd (,,9) = Beta V
      !!                  vortrd (,,10) = forcing term
      !!		  vortrd (,,11) = bottom friction term
      !!                  rotot(,) : total cumulative trends over nwrite-1 time steps
      !!                  vor_avrtot(,) : first membre of vrticity equation
      !!                  vor_avrres(,) : residual = dh/dt entrainment
      !!
      !!      trends output in netCDF format using ioipsl
      !!
      !! History :
      !!   9.0  !  04-06  (L. Brunier, A-M. Treguier) Original code 
      !!        !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index

      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
         putrdvor,                         &  ! u vorticity trend 
         pvtrdvor                             ! v vorticity trend

      !! * Local declarations
      INTEGER ::   ji, jj, jk

      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zubet,                         &  ! u Beta.V case
         zvbet,                         &  ! v Beta.V case
         zudpvor,                       &  ! total cmulative trends
         zvdpvor                           !   "      "        "
      !!----------------------------------------------------------------------
     
      ! Initialization
      zubet(:,:) = 0.e0
      zvbet(:,:) = 0.e0
      zudpvor(:,:) = 0.e0
      zvdpvor(:,:) = 0.e0

      !  =====================================
      !  I vertical integration of 3D trends
      !  =====================================

      CALL lbc_lnk( putrdvor, 'U' , -1. )
      CALL lbc_lnk( pvtrdvor, 'V' , -1. )

      ! putrdvor and pvtrdvor terms
      DO jk = 1,jpk
        zudpvor(:,:) = zudpvor(:,:) + putrdvor(:,:,jk) * fse3u(:,:,jk) * e1u(:,:) * umask(:,:,jk)
        zvdpvor(:,:) = zvdpvor(:,:) + pvtrdvor(:,:,jk) * fse3v(:,:,jk) * e2v(:,:) * vmask(:,:,jk)
      END DO

      ! Save Beta.V term to avoid average before Curl
      ! Beta.V : intergration, no average
      IF( ktrd == jpvorbev ) THEN 
         zubet(:,:) = zudpvor(:,:)
         zvbet(:,:) = zvdpvor(:,:)
      ENDIF

      ! Average except for Beta.V
      zudpvor(:,:) = zudpvor(:,:) * hur(:,:)
      zvdpvor(:,:) = zvdpvor(:,:) * hvr(:,:)
   
      ! Curl
      DO ji=1,jpim1
         DO jj=1,jpjm1
            vortrd(ji,jj,ktrd) = (  zvdpvor(ji+1,jj) - zvdpvor(ji,jj) -   &
                 &                ( zudpvor(ji,jj+1) - zudpvor(ji,jj) ) ) &
                 &               / ( e1f(ji,jj) * e2f(ji,jj) )
         END DO
      END DO

      ! Surface mask
      vortrd(:,:,ktrd) = vortrd(:,:,ktrd) * fmask(:,:,1)

      ! Special treatement for the Beta.V term
      ! Compute the Curl of the Beta.V term which is not averaged
      IF( ktrd == jpvorbev ) THEN
         DO ji=1,jpim1
            DO jj=1,jpjm1
               vortrd(ji,jj,jpvorbev) = (  zvbet(ji+1,jj) - zvbet(ji,jj) -   &
                    &                    ( zubet(ji,jj+1) - zubet(ji,jj) ) ) &
                    &                   / ( e1f(ji,jj) * e2f(ji,jj) )
            END DO
         END DO

         ! Average on the Curl
         vortrd(:,:,jpvorbev) = vortrd(:,:,jpvorbev) * hur(:,:)

         ! Surface mask
         vortrd(:,:,jpvorbev) = vortrd(:,:,jpvorbev) * fmask(:,:,1)
      ENDIF
   
      IF( idebug /= 0 ) THEN
         IF(lwp) WRITE(numout,*) ' debuging trd_vor_zint: I done'
         CALL FLUSH(numout)
      ENDIF

   END SUBROUTINE trd_vor_zint_3d



   SUBROUTINE trd_vor( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_vor  ***
      !! 
      !! ** Purpose :  computation of cumulated trends over analysis period
      !!               and make outputs (NetCDF or DIMG format)
      !!
      !! ** Method/usage :
      !!
      !! History :
      !!   9.0  !  04-06  (L. Brunier, A-M. Treguier) Original code 
      !!        !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      !! * Local declarations
      INTEGER :: ji, jj, jk, jl, it

      REAL(wp) :: zmean

      REAL(wp) ,DIMENSION(jpi,jpj) ::   &
         zun, zvn
      !!----------------------------------------------------------------------

      !  =================
      !  I. Initialization
      !  =================
     
     
      ! I.1 set before values of vertically average u and v
      ! ---------------------------------------------------

      IF( kt > nit000 ) THEN
         vor_avrb(:,:) = vor_avr(:,:)
      ENDIF

       IF( idebug /= 0 ) THEN
          if(lwp)WRITE(numout,*) ' debuging trd_vor: I.1 done '
          CALL FLUSH(numout)
      ENDIF

      ! I.2 vertically integrated vorticity
      !  ----------------------------------

      vor_avr(:,:) = 0.
      zun(:,:)=0
      zvn(:,:)=0
      vor_avrtot(:,:)=0
      vor_avrres(:,:)=0
      
      ! Vertically averaged velocity
      DO jk = 1, jpk - 1
         zun(:,:)=zun(:,:) + e1u(:,:)*un(:,:,jk)*fse3u(:,:,jk)
         zvn(:,:)=zvn(:,:) + e2v(:,:)*vn(:,:,jk)*fse3v(:,:,jk)
      END DO
 
      zun(:,:)=zun(:,:)*hur(:,:)
      zvn(:,:)=zvn(:,:)*hvr(:,:)

      ! Curl
      DO ji=1,jpim1
         DO jj=1,jpjm1
            vor_avr(ji,jj) = ((zvn(ji+1,jj)-zvn(ji,jj))-   &
                              (zun(ji,jj+1)-zun(ji,jj)))   &
                             /( e1f(ji,jj) * e2f(ji,jj) )
            vor_avr(ji,jj) = vor_avr(ji,jj)*fmask(ji,jj,1)
         END DO
      END DO
      
      IF(idebug /= 0) THEN
         if(lwp) WRITE(numout,*) ' debuging trd_vor: I.2 done'
         CALL FLUSH(numout)
      ENDIF

      !  =================================
      !   II. Cumulated trends
      !  =================================

      ! II.1 set `before' mixed layer values for kt = nit000+1
      ! ------------------------------------------------------
      IF( kt == nit000+1 ) THEN
         vor_avrbb(:,:) = vor_avrb(:,:)
         vor_avrbn(:,:) = vor_avr (:,:)
      ENDIF

      IF( idebug /= 0 ) THEN
         if(lwp)  WRITE(numout,*) ' debuging trd_vor: I1.1 done'
         CALL FLUSH(numout)
      ENDIF

      ! II.2 cumulated trends over analysis period (kt=2 to nwrite)
      ! ----------------------
      ! trends cumulated over nwrite-2 time steps

      IF( kt >= nit000+2 ) THEN
         nmoydpvor = nmoydpvor + 1
         DO jl = 1, jplvor
            IF( jl /= 9 ) THEN
               rotot(:,:) = rotot(:,:) + vortrd(:,:,jl)
            ENDIF
         END DO
      ENDIF

      IF( idebug /= 0 ) THEN
         if(lwp)  WRITE(numout,*) ' debuging trd_vor: II.2 done'
         CALL FLUSH(numout)
      ENDIF

      !  =============================================
      !   III. Output in netCDF + residual computation
      !  =============================================

      IF( MOD( kt - nit000+1, ntrd ) == 0 ) THEN

         ! III.1 compute total trend
         ! ------------------------
         zmean = float(nmoydpvor)

         vor_avrtot(:,:) = ( vor_avr(:,:) - vor_avrbn(:,:) + vor_avrb(:,:) - &
                             vor_avrbb(:,:) ) /  (zmean * 2. * rdt)

         IF( idebug /= 0 ) THEN
             if(lwp)  WRITE(numout,*) ' zmean = ',zmean
             if(lwp)  WRITE(numout,*) ' debuging trd_vor: III.1 done'
             CALL FLUSH(numout)
         ENDIF

         ! III.2 compute residual
         ! ---------------------
         vor_avrres(:,:) = vor_avrtot(:,:) - rotot(:,:) / zmean

         ! Boundary conditions
         CALL lbc_lnk( vor_avrtot, 'F', 1. )
         CALL lbc_lnk( vor_avrres, 'F', 1. )

         IF( idebug /= 0 ) THEN
            if(lwp)  WRITE(numout,*) ' debuging trd_vor: III.2 done'
            CALL FLUSH(numout)
         ENDIF

         ! III.3 time evolution array swap
         ! ------------------------------
         vor_avrbb(:,:) = vor_avrb(:,:)
         vor_avrbn(:,:) = vor_avr(:,:)

         IF( idebug /= 0 ) THEN
            if(lwp)  WRITE(numout,*) ' debuging trd_vor: III.3 done'
            CALL FLUSH(numout)
         ENDIF

         nmoydpvor=0

      ENDIF

      ! III.4 write trends to output
      ! ---------------------------

      IF( kt >=  nit000+1 ) THEN

         ! define time axis
         it= kt-nit000+1
         IF( lwp .AND. MOD( kt, ntrd ) == 0 ) THEN
            WRITE(numout,*) '     trdvor_ncwrite : write NetCDF fields'
         ENDIF
 
         CALL histwrite( nidvor,"sovortPh",it,vortrd(:,:,1),ndimvor1,ndexvor1)  ! grad Ph
         CALL histwrite( nidvor,"sovortEk",it,vortrd(:,:,2),ndimvor1,ndexvor1)  ! Energy
         CALL histwrite( nidvor,"sovozeta",it,vortrd(:,:,3),ndimvor1,ndexvor1)  ! rel vorticity
         CALL histwrite( nidvor,"sovortif",it,vortrd(:,:,4),ndimvor1,ndexvor1)  ! coriolis
         CALL histwrite( nidvor,"sovodifl",it,vortrd(:,:,5),ndimvor1,ndexvor1)  ! lat diff
         CALL histwrite( nidvor,"sovoadvv",it,vortrd(:,:,6),ndimvor1,ndexvor1)  ! vert adv
         CALL histwrite( nidvor,"sovodifv",it,vortrd(:,:,7),ndimvor1,ndexvor1)  ! vert diff
         CALL histwrite( nidvor,"sovortPs",it,vortrd(:,:,8),ndimvor1,ndexvor1)  ! grad Ps
         CALL histwrite( nidvor,"sovortbv",it,vortrd(:,:,9),ndimvor1,ndexvor1)  ! beta.V
         CALL histwrite( nidvor,"sovowind",it,vortrd(:,:,10),ndimvor1,ndexvor1) ! wind stress
         CALL histwrite( nidvor,"sovobfri",it,vortrd(:,:,11),ndimvor1,ndexvor1) ! bottom friction
         CALL histwrite( nidvor,"1st_mbre",it,vor_avrtot    ,ndimvor1,ndexvor1) ! First membre
         CALL histwrite( nidvor,"sovorgap",it,vor_avrres    ,ndimvor1,ndexvor1) ! gap between 1st and 2 nd mbre

         IF( idebug /= 0 ) THEN
            if(lwp)  WRITE(numout,*) ' debuging trd_vor: III.4 done'
            CALL FLUSH(numout)
         ENDIF

      ENDIF

      IF( MOD( kt - nit000+1, ntrd ) == 0 ) rotot(:,:)=0

      IF( kt == nitend )   CALL histclo( nidvor )

   END SUBROUTINE trd_vor



   SUBROUTINE trd_vor_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_vor_init  ***
      !! 
      !! ** Purpose :   computation of vertically integrated T and S budgets
      !!      from ocean surface down to control surface (NetCDF output)
      !!
      !! ** Method/usage :
      !!
      !! History :
      !!   9.0  !  04-06  (L. Brunier, A-M. Treguier) Original code 
      !!        !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Local declarations
      REAL(wp) :: zjulian, zsto, zout

      CHARACTER (len=40) ::   clhstnam
      CHARACTER (len=40) ::   clop

      NAMELIST/namtrd/ ntrd,nctls
      !!----------------------------------------------------------------------

      !  ===================
      !   I. initialization
      !  ===================

      cvort='averaged-vor'

      ! Open specifier
      idebug = 0      ! set it to 1 in case of problem to have more Print

      ! namelist namtrd : trend diagnostic
      REWIND( numnam )
      READ  ( numnam, namtrd )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) 'trd_vor_init: vorticity trends'
         WRITE(numout,*) '~~~~~~~~~~~~~'
         WRITE(numout,*) ' '
         WRITE(numout,*) '          Namelist namtrd : '
         WRITE(numout,*) '             time step frequency trend       ntrd  = ',ntrd
         WRITE(numout,*) ' '
         WRITE(numout,*) '##########################################################################'
         WRITE(numout,*) ' CAUTION: The interpretation of the vorticity trends is'
         WRITE(numout,*) ' not obvious, please contact Anne-Marie TREGUIER at: treguier@ifremer.fr '
         WRITE(numout,*) '##########################################################################'
         WRITE(numout,*) ' '
      ENDIF

      ! cumulated trends array init
      nmoydpvor = 0
      rotot(:,:)=0
      vor_avrtot(:,:)=0
      vor_avrres(:,:)=0

      IF( idebug /= 0 ) THEN
         if(lwp)  WRITE(numout,*) ' debuging trd_vor_init: I. done'
         CALL FLUSH(numout)
      ENDIF

      !  =================================
      !   II. netCDF output initialization
      !  =================================

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
      zout = ntrd*rdt

      IF(lwp) WRITE (numout,*) ' trdvor_ncinit: netCDF initialization'

      ! II.2 Compute julian date from starting date of the run
      ! ------------------------
      CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
      IF (lwp) WRITE(numout,*)' '  
      IF (lwp) WRITE(numout,*)' Date 0 used :',nit000         &
           ,' YEAR ', nyear,' MONTH ', nmonth,' DAY ', nday   &
           ,'Julian day : ', zjulian

      ! II.3 Define the T grid trend file (nidvor)
      ! ---------------------------------
      CALL dia_nam( clhstnam, ntrd, 'vort' )                  ! filename
      IF(lwp) WRITE(numout,*) ' Name of NETCDF file ', clhstnam
      CALL histbeg( clhstnam, jpi, glamf, jpj, gphif,1, jpi,   &  ! Horizontal grid : glamt and gphit
         &          1, jpj, 0, zjulian, rdt, nh_t, nidvor, domain_id=nidom )
      CALL wheneq( jpi*jpj, fmask, 1, 1., ndexvor1, ndimvor1 )    ! surface

      ! Declare output fields as netCDF variables
      CALL histdef( nidvor, "sovortPh", cvort//"grad Ph" , "s-2",        & ! grad Ph
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovortEk", cvort//"Energy", "s-2",          & ! Energy
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovozeta", cvort//"rel vorticity", "s-2",   & ! rel vorticity
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovortif", cvort//"coriolis", "s-2",        & ! coriolis
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovodifl", cvort//"lat diff ", "s-2",       & ! lat diff
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovoadvv", cvort//"vert adv", "s-2",        & ! vert adv
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovodifv", cvort//"vert diff" , "s-2",      & ! vert diff
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovortPs", cvort//"grad Ps", "s-2",         & ! grad Ps
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovortbv", cvort//"Beta V", "s-2",          & ! beta.V
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovowind", cvort//"wind stress", "s-2",     & ! wind stress
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovobfri", cvort//"bottom friction", "s-2", & ! bottom friction
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "1st_mbre", cvort//"1st mbre", "s-2",        & ! First membre
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histdef( nidvor, "sovorgap", cvort//"gap", "s-2",             & ! gap between 1st and 2 nd mbre
         &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
      CALL histend( nidvor )

      IF( idebug /= 0 ) THEN
         if(lwp)  WRITE(numout,*) ' debuging trd_vor_init: II. done'
         CALL FLUSH(numout)
      ENDIF

   END SUBROUTINE trd_vor_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC ::   lk_trdvor = .FALSE.   ! momentum trend flag

   !! * Interfaces
   INTERFACE trd_vor_zint
      MODULE PROCEDURE trd_vor_zint_2d, trd_vor_zint_3d
   END INTERFACE

CONTAINS
   SUBROUTINE trd_vor( kt )        ! Empty routine
!      WRITE(*,*) 'trd_vor: You should not have seen this print! error?', kt
   END SUBROUTINE trd_vor
   SUBROUTINE trd_vor_zint_2d( putrdvor, pvtrdvor, ktrd )
      REAL, DIMENSION(:,:), INTENT( inout ) ::   &
         putrdvor, pvtrdvor                  ! U and V momentum trends
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
!      WRITE(*,*) 'trd_vor_zint_2d: You should not have seen this print! error?', putrdvor(1,1)
!      WRITE(*,*) '  "      "     : You should not have seen this print! error?', pvtrdvor(1,1)
!      WRITE(*,*) '  "      "     : You should not have seen this print! error?', ktrd
   END SUBROUTINE trd_vor_zint_2d
   SUBROUTINE trd_vor_zint_3d( putrdvor, pvtrdvor, ktrd )
      REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
         putrdvor, pvtrdvor                  ! U and V momentum trends
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
!      WRITE(*,*) 'trd_vor_zint_3d: You should not have seen this print! error?', putrdvor(1,1,1)
!      WRITE(*,*) '  "      "     : You should not have seen this print! error?', pvtrdvor(1,1,1)
!      WRITE(*,*) '  "      "     : You should not have seen this print! error?', ktrd
   END SUBROUTINE trd_vor_zint_3d
   SUBROUTINE trd_vor_init              ! Empty routine
!      WRITE(*,*) 'trd_vor_init: You should not have seen this print! error?'
   END SUBROUTINE trd_vor_init
#endif
   !!======================================================================
END MODULE trdvor
