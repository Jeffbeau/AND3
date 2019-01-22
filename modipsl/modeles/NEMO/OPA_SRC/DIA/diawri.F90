!!DB -- 2009.09.04 -- key_diadimg eliminated
!!DB - 2007.12.04
!!This version has eliminated the IOIPSL code


MODULE diawri
   !!======================================================================
   !!                     ***  MODULE  diawri  ***
   !! Ocean diagnostics :  write ocean output files
   !!=====================================================================

   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE zdf_oce         ! ocean vertical physics
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE sol_oce         ! solver variables
   USE ice_oce         ! ice variables
   USE phycst          ! physical constants
   USE ocfzpt          ! ocean freezing point
   USE ocesbc          ! surface thermohaline fluxes
   USE taumod          ! surface stress
   USE flxrnf          ! ocean runoffs
   USE zdfmxl          ! mixed layer
   USE daymod          ! calendar
   USE dianam          ! build name of file (routine)
   USE zdfddm          ! vertical  physics: double diffusion
   USE diahth          ! thermocline diagnostics
   USE diaspr          ! surface pressure diagnostics (rigid lid case)
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE flx_oce         ! sea-ice/ocean forcings variables

   USE lib_ncdf        ! netCDF I/O library

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dia_wri                 ! routines called by step.F90
   PUBLIC dia_wri_state
!!DB: tried to move to lib_ncdf but did not work (???)
   PUBLIC output_special             ! special routines called by step.F90 if ave flag is on
   PUBLIC output_aveTSUV             ! special routines called by step.F90 if ave flag is on

   !! * Module variables
!   INTEGER ::   &
!      nid_T, nz_T, nh_T, ndim_T, ndim_hT,      &   ! grid_T file
!      nid_U, nz_U, nh_U, ndim_U, ndim_hU,      &   ! grid_U file
!      nid_V, nz_V, nh_V, ndim_V, ndim_hV,      &   ! grid_V file
!      nid_W, nz_W, nh_W,                       &   ! grid_W file
!      ndex(1)                                      ! ???
!   INTEGER, DIMENSION(jpi*jpj) ::   &
!      ndex_hT, ndex_hU, ndex_hV
!   INTEGER, DIMENSION(jpi*jpj*jpk) ::   &
!      ndex_T, ndex_U, ndex_V

   !! * Substitutions
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diawri.F90,v 1.10 2006/03/09 17:21:54 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !!   Default option                                   NetCDF output file
   !!----------------------------------------------------------------------
   !!   dia_wri       : create the standart NetCDF output files
   !!   dia_wri_state : create an output NetCDF file for a single
   !!                   instantaeous ocean state and forcing fields
   !!----------------------------------------------------------------------

   SUBROUTINE dia_wri( kt, kindic )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!                   
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!
      !! ** Method  :   At the beginning of the first time step (nit000), 
      !!      define all the NETCDF files and fields
      !!      At each time step call histdef to compute the mean if ncessary
      !!      Each nwrite time step, output the instantaneous or mean fields
      !!      IF kindic <0, output of fields before the model interruption.
      !!      IF kindic =0, time step loop
      !!      IF kindic >0, output of fields before the time step loop
      !!
      !! History :
      !!        !  91-03  (M.-A. Foujols)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)  correction restart file
      !!        !  92-07  (M. Imbard)  split into diawri and rstwri
      !!        !  93-03  (M. Imbard)  suppress writibm
      !!        !  98-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
      !!        !  99-02  (E. Guilyardi)  name of netCDF files + variables
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
!      USE ioipsl

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   kindic  ! 

      !! * Local declarations
      LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
      CHARACTER (len=40) ::           &
         clhstnam, clop, clmx            ! temporary names
      INTEGER ::   inum = 11             ! temporary logical unit
      INTEGER ::   &
         iimi, iima, ipk, it,         &  ! temporary integers
         ijmi, ijma                      !    "          "

      INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls
      LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.    ! Use IOIPSL subroutines for output
                                                 ! Note that all W-grid calls are commented out!
      LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE.   ! Use lib_ncdf subroutines for output

      REAL(wp) :: test_array
      INTEGER :: a, b, c

      REAL(wp) ::   &
         zsto, zout, zmax,            &  ! temporary scalars
         zjulian, zdt                    !    "         "
      REAL(wp), DIMENSION(jpi,jpj) :: &
         zw2d                            ! temporary workspace
      CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT
!DB
      INTEGER :: tindex


      !!----------------------------------------------------------------------
      
      ! 0. Initialisation
      ! -----------------
      fnameU = trim(cexper)//'_grid_U.nc'
      fnameV = trim(cexper)//'_grid_V.nc'
      fnameT = trim(cexper)//'_grid_T.nc'

!!DB deleted ...
!      IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL for output?

   ! LIB_NCDF CALLS

!DB default to this code
!   IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Are we using lib_ncdf output routines?

      ! On the first timestep, create output files and write depth variables (since they don't change)
      IF(kt == nit000) THEN
!!DBG
         write(numout2,*)'DBG: diawri: cexper = ', cexper
         write(numout2,*)'DBG: diawri: ',fnameU, fnameV, fnameT

         CALL ncdf_create_file_u(fnameU, 'inst(x)', lncdf_stat)
         CALL ncdf_create_file_v(fnameV, 'inst(x)', lncdf_stat)
         CALL ncdf_create_file_t(fnameT, 'inst(x)', lncdf_stat)
         CALL ncdf_write(fnameU, 'depthu', gdept, lncdf_stat) 
         CALL ncdf_write(fnameV, 'depthv', gdept, lncdf_stat) 
         CALL ncdf_write(fnameT, 'deptht', gdept, lncdf_stat) 
      END IF
      
      ! On first write step only, write nav_lat and nav_lon grids (routine will abort if not on a write step)
      IF(kt == (nit000+ nwrite)) THEN
         CALL ncdf_write(fnameU, 'nav_lat', gphiu, kt, lncdf_stat)
         CALL ncdf_write(fnameU, 'nav_lon', glamu, kt, lncdf_stat)
         CALL ncdf_write(fnameV, 'nav_lat', gphiv, kt, lncdf_stat)
         CALL ncdf_write(fnameV, 'nav_lon', glamv, kt, lncdf_stat)
         CALL ncdf_write(fnameT, 'nav_lat', gphit, kt, lncdf_stat)
         CALL ncdf_write(fnameT, 'nav_lon', glamt, kt, lncdf_stat)
      END IF
      
      ! On write steps, output time counter
      IF(MOD(kt, nwrite) == 0) THEN
         tindex = kt/nwrite - nit000/nwrite
!         CALL ncdf_write(fnameU, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
!         CALL ncdf_write(fnameV, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
!         CALL ncdf_write(fnameT, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
         CALL ncdf_write(fnameU, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
         CALL ncdf_write(fnameV, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
         CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
      END IF
      
      ! Write other fields (calls will abort if this isn't a write step)
      CALL ncdf_write(fnameU, 'vozocrtx', un, kt, lncdf_stat)
      CALL ncdf_write(fnameU, 'sozotaux', taux, kt, lncdf_stat)
      CALL ncdf_write(fnameV, 'vomecrty', vn, kt, lncdf_stat)
      CALL ncdf_write(fnameV, 'sometauy', tauy, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'votemper', tn, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'vosaline', sn, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sosstsst', tn(:,:,1), kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sosaline', sn(:,:,1), kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sossheig', sshn, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sowaflup', emp, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sorunoff', runoff, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sowaflcd', emps, kt, lncdf_stat)
      zw2d(:,:) = emps(:,:) * sn(:,:,1) * tmask(:,:,1)
      CALL ncdf_write(fnameT, 'sosalflx', zw2d, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sohefldo', qt, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'soshfldo', qsr, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'somxl010', hmlp, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'somixhgt', hmld, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'soicecov', freeze, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sohefldp', qrp, kt, lncdf_stat)
      CALL ncdf_write(fnameT, 'sowafldp', erp, kt, lncdf_stat)
      zw2d(:,:) = erp(:,:) * sn(:,:,1) * tmask(:,:,1)
      CALL ncdf_write(fnameT, 'sosafldp', zw2d, kt, lncdf_stat)
      zw2d(:,:) = FLOAT( nmln(:,:) ) * tmask(:,:,1)
      CALL ncdf_write(fnameT, 'sobowlin', zw2d, kt, lncdf_stat)
      
   ! END OF LIB_NCDF CALLS

!DB code deleted
!   IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL output subroutines?

   END SUBROUTINE dia_wri


   SUBROUTINE dia_wri_state( cdfile_name )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.nc' is created in case of abnormal job end
      !!
      !! History :
      !!   8.2  !  00-06  (M. Imbard)  Original code (diabort.F)
      !!   8.5  !  02-06  (A.Bozec, E. Durand)  Original code (diainit.F)
      !!   9.0  !  02-12  (G. Madec)  merge of diabort and diainit, F90
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Arguments
      CHARACTER (len=* ), INTENT( in ) ::   &
         cdfile_name      ! name of the file created

      !! * Local declarations
      CHARACTER (len=40) :: clop
      INTEGER  ::   &
         id_i , nz_i, nh_i       
      INTEGER, DIMENSION(1) ::   &
         idex             ! temprary workspace
      REAL(wp) ::   &
         zsto, zout, zmax,   &
         zjulian, zdt
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
      IF(lwp) WRITE(numout,*) '                and named :', cdfile_name, '.nc'
      
      ! 0. Initialisation
      ! -----------------
      
      ! Define frequency of output and means
      zdt  = rdt
      zsto = rdt
      clop = "inst(x)"           ! no use of the mask value (require less cpu time)
      zout = rdt
      zmax = ( nitend - nit000 + 1 ) * zdt

      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      ! Compute julian date from starting date of the run
      CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )         ! time axis 
      CALL histbeg( cdfile_name, jpi, glamt, jpj, gphit,   &
          1, jpi, 1, jpj, 0, zjulian, zdt, nh_i, id_i, domain_id=nidom )          ! Horizontal grid : glamt and gphit
      CALL histvert( id_i, "deptht", "Vertical T levels",   &    ! Vertical grid : gdept
          "m", jpk, gdept, nz_i)

      ! Declare all the output fields as NetCDF variables

      CALL histdef( id_i, "vosaline", "Salinity"              , "PSU"    ,   &   ! salinity
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "votemper", "Temperature"           , "C"      ,   &   ! temperature
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
#if defined key_dynspg_rl
      CALL histdef( id_i, "sobarstf","Barotropic StreamFunction", "m3/s2"  ,   &  ! bsf
         &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
#else
      CALL histdef( id_i, "sossheig", "Sea Surface Height"    , "m"      ,   &  ! ssh
         &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
#endif
      CALL histdef( id_i, "vozocrtx", "Zonal Current"         , "m/s"    ,   &   ! zonal current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "vomecrty", "Meridional Current"    , "m/s"    ,   &   ! meridonal current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
      CALL histdef( id_i, "vovecrtz", "Vertical Velocity"     , "m/s"    ,   &   ! vertical current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
      CALL histdef( id_i, "sowaflup", "Net Upward Water Flux" , "Kg/m2/S",   &   ! net freshwater 
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sohefldo", "Net Downward Heat Flux", "W/m2"   ,   &   ! net heat flux
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "soshfldo", "Shortwave Radiation"   , "W/m2"   ,   &   ! solar flux
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "soicecov", "Ice fraction"          , "[0,1]"  ,   &   ! freeze
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sozotaux", "Zonal Wind Stress"     , "N/m2"   ,   &   ! i-wind stress
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sometauy", "Meridional Wind Stress", "N/m2"   ,   &   ! j-wind stress
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )

      CALL histend( id_i )

      ! 2. Start writing data
      ! ---------------------
      ! idex(1) est utilise ssi l'avant dernier argument est diffferent de 
      ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
      ! donne le nombre d'elements, et idex la liste des indices a sortir
      idex(1) = 1   ! init to avoid compil warning
      
      ! Write all fields on T grid
      CALL histwrite( id_i, "votemper", 1, tn    , jpi*jpj*jpk, idex )    ! now temperature
      CALL histwrite( id_i, "vosaline", 1, sn    , jpi*jpj*jpk, idex )    ! now salinity
#if defined key_dynspg_rl
      CALL histwrite( id_i, "sobarstf", 1, bsfn  , jpi*jpj    , idex )    ! barotropic streamfunction
#else
      CALL histwrite( id_i, "sossheig", 1, sshn  , jpi*jpj    , idex )    ! sea surface height
#endif
      CALL histwrite( id_i, "vozocrtx", 1, un    , jpi*jpj*jpk, idex )    ! now i-velocity
      CALL histwrite( id_i, "vomecrty", 1, vn    , jpi*jpj*jpk, idex )    ! now j-velocity
      CALL histwrite( id_i, "vovecrtz", 1, wn    , jpi*jpj*jpk, idex )    ! now k-velocity
      CALL histwrite( id_i, "sowaflup", 1, emp   , jpi*jpj    , idex )    ! freshwater budget
      CALL histwrite( id_i, "sohefldo", 1, qt    , jpi*jpj    , idex )    ! total heat flux
      CALL histwrite( id_i, "soshfldo", 1, qsr   , jpi*jpj    , idex )    ! total heat flux
      CALL histwrite( id_i, "soicecov", 1, freeze, jpi*jpj    , idex )    ! ice cover
      CALL histwrite( id_i, "sozotaux", 1, taux  , jpi*jpj    , idex )    ! i-wind stress
      CALL histwrite( id_i, "sometauy", 1, tauy  , jpi*jpj    , idex )    ! j-wind stress

      ! 3. Close the file
      ! -----------------
      CALL histclo( id_i )

   END SUBROUTINE dia_wri_state

   !!======================================================================

!!DB
   SUBROUTINE output_special( kt, kindic )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri_ave  ***
      !!                   
      !! ** Purpose :   Average and output special opa variable
      !!      NETCDF format is used by default 
      !!
      !! DB: 2007.12.04 ... 12.10 ... tested OK
      !!----------------------------------------------------------------------
!!DB -- notes
!1 - due to the way that lib_ncdf routine handles 'when to write' this routine will
! have to fool the call to ncdf_write...() to write at a non-nwrite timestep:
! fixed by inputting a -ve rec_num value
!2 - ALLOCATE arrays when routine is called ... DONE ...
!3 - Uses an improved M2-averaging routine (see dbrick:..MATLAB/test_M2ave.m)
! that interpolates to the correct (float)time (rM2dt). Routine also gives correct
! answer if rM2dt = integer

      !! * Modules used
!      USE ioipsl
!     USE oce
!     USE taumod

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   kindic  ! 

      !! * Local declarations
      LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
      CHARACTER (len=40) ::           &
         clhstnam, clop, clmx            ! temporary names
      INTEGER ::   inum = 11             ! temporary logical unit
      INTEGER ::   &
         iimi, iima, ipk, it,         &  ! temporary integers
         ijmi, ijma                      !    "          "

      INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls

      REAL(wp) ::   &
         zsto, zout, zmax,            &  ! temporary scalars
         zjulian, zdt                    !    "         "
      REAL(wp), DIMENSION(jpi,jpj) :: &
         zw2d                            ! temporary workspace
      CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT

!!DB
      INTEGER :: kkk, M2dt, ndt, ddt, M2dt_lo, M2dt_hi
      REAL(wp) :: M2_hr, rM2dt, zxy
      REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: aveU, aveV, aveT
      REAL(wp),DIMENSION(:,:),ALLOCATABLE, SAVE :: ave_ssh
      INTEGER, SAVE :: rec_num

      ! 0. Initialisation
!!DB
      M2_hr = 12.42
      ndt = 86400./rdt
      rM2dt = M2_hr*3600./rdt
      M2dt_lo = floor(M2_hr*3600./rdt)
!      M2dt_hi = ceil(M2_hr*3600./rdt) !no ceil() but variable is not needed
      M2dt = M2dt_lo;
      zxy = rM2dt-float(M2dt)

      fnameU = trim(cexper)//'_M2ave_U.nc'
      fnameV = trim(cexper)//'_M2ave_V.nc'
      fnameT = trim(cexper)//'_M2ave_T.nc'
!      fnameT = trim(cexper)//'_ave_TSUV.nc'

! On the first timestep, create output files and write depth variables (since they don't change)
      IF(kt == nit000) THEN
         rec_num = 0

         call ncdf_create_file_ave(fnameU, 'aveU', 'ave(x)', lncdf_stat)
         call ncdf_create_file_ave(fnameV, 'aveV', 'ave(x)', lncdf_stat)
!         call ncdf_create_file_ave(fnameT, 'aveT', 'ave(x)', lncdf_stat)

         ALLOCATE(aveU(jpi,jpj,jpk))
         ALLOCATE(aveV(jpi,jpj,jpk))
!         ALLOCATE(aveT(jpi,jpj,jpk))
         aveU = 0.0
         aveV = 0.0

!!DBG
         if(lwp) then
            write(numout2,*)'DBG: diawri: cexper = ', cexper
            write(numout2,*)'DBG: diawri: ',fnameU, fnameV
         endif


         CALL ncdf_write(fnameU, 'depthu', gdept, lncdf_stat) 
         CALL ncdf_write(fnameV, 'depthu', gdept, lncdf_stat) 
!         CALL ncdf_write(fnameT, 'deptht', gdept, lncdf_stat) 
      
! On first write step only, write nav_lat and nav_lon grids (routine will abort if not on a write step)
         CALL ncdf_write(fnameU, 'nav_lat', gphiu, -1, lncdf_stat)
         CALL ncdf_write(fnameU, 'nav_lon', glamu, -1, lncdf_stat)
         CALL ncdf_write(fnameV, 'nav_lat', gphiv, -1, lncdf_stat)
         CALL ncdf_write(fnameV, 'nav_lon', glamv, -1, lncdf_stat)
!         CALL ncdf_write(fnameT, 'nav_lat', gphit, -1, lncdf_stat)
!         CALL ncdf_write(fnameT, 'nav_lon', glamt, -1, lncdf_stat)

      END IF  ! first timestep

!!DB: 2007.12.04 -- Code to compute M2 ave fields, once per day
      ddt = mod(kt-nit000,ndt) + 1  !! = timestep-# within a day
!DB 12.10 -- the below actually gives 1 more dt than I wanted i.e.
! M2dt + 1, but now I need this extra dt
      if(ddt >= ndt-M2dt) then  !!accumulate sums
         if(ddt /= ndt) then    !! but not on the last dt of the window
            aveU = aveU + un
            aveV = aveV + vn
         endif
      endif

!!Output time if:
!!DB 12.10 -- the code below time-interpolates to the exact rM2dt, i.e.
!! val = (1-zxy)*sum_lo + zxy*sum_hi, where sum_hi = sum_lo + latest_field
!!     = (1-zxy)*sum_lo + zxy*(sum_lo+lf) = sum_lo + zxy*lf
!! val = val/rM2dt
!! This method avoids having to store any temporary 3D fields
      if(ddt == ndt) then
         rec_num = rec_num + 1

         aveU = (aveU + zxy*un)/rM2dt
         aveV = (aveV + zxy*vn)/rM2dt
!         aveU = aveU/float(M2dt)
!         aveV = aveV/float(M2dt)
         
         ! On write steps, output time counter
         CALL ncdf_write(fnameU, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)
         CALL ncdf_write(fnameV, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)
!         CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)

         CALL ncdf_write(fnameU, 'aveU', aveU, -rec_num, lncdf_stat)
         CALL ncdf_write(fnameV, 'aveV', aveV, -rec_num, lncdf_stat)
         
         aveU = 0.0
         aveV = 0.0
!!AD/DB 2009.09.30
         CALL ncdf_write(fnameU, 'ndastp',REAL(ndastp), rec_num, lncdf_stat)
         CALL ncdf_write(fnameU, 'model_time_step',REAL(kt), rec_num, lncdf_stat)
         CALL ncdf_write(fnameU, 'model_time',model_time, rec_num, lncdf_stat)
         CALL ncdf_write(fnameV, 'ndastp',REAL(ndastp), rec_num, lncdf_stat)
         CALL ncdf_write(fnameV, 'model_time_step',REAL(kt), rec_num, lncdf_stat)
         CALL ncdf_write(fnameV, 'model_time',model_time, rec_num, lncdf_stat)

      endif
      


   END SUBROUTINE output_special

!!DB
!!write time-averaged variables if the ioutput_ave flag > 0
!! averaging period = ioutput_ave; output is every mod(kt-nit000,ioutput_ave)
!!NB: More variables are created in the output file than are currently on the 
!!average list. 
!!2008.06.13 -- add w (but do not change name)

   SUBROUTINE output_aveTSUV( kt, kindic )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri_ave  ***
      !!                   
      !! ** Purpose :   Average and output special opa variable
      !!      NETCDF format is used by default 
      !!
      !! DB: 2007.12.11 
      !!----------------------------------------------------------------------

      !! * Modules used
!      USE ioipsl
!     USE oce
!     USE taumod, ONLY:  taux, tauy

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   kindic  ! 

      !! * Local declarations
      LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
      CHARACTER (len=40) ::           &
         clhstnam, clop, clmx            ! temporary names
      INTEGER ::   inum = 11             ! temporary logical unit
      INTEGER ::   &
         iimi, iima, ipk, it,         &  ! temporary integers
         ijmi, ijma                      !    "          "

      INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls

      REAL(wp) ::   &
         zsto, zout, zmax,            &  ! temporary scalars
         zjulian, zdt                    !    "         "
      REAL(wp), DIMENSION(jpi,jpj) :: &
         zw2d                            ! temporary workspace
      CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT

!!DB  
      REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: aveU2, aveV2, aveT2, aveS2, aveW2
      REAL(wp),DIMENSION(:,:),ALLOCATABLE, SAVE :: ave_ssh2, ave_taux, ave_tauy
      INTEGER, SAVE :: rec_num2

!DB  This routine should not be called if ioutput_ave = 0, but to be safe
      if(ioutput_ave == 0) then
         return
      endif

      ! 0. Initialisation
      fnameT = trim(cexper)//'_ave_TSUV.nc'

! On the first timestep, create output files and write depth variables (since they don't change)
      IF(kt == nit000) THEN
         rec_num2 = 0
         call ncdf_create_file_aveTSUV(fnameT, 'ave(x)', lncdf_stat)

! Allocate arrays, and to be safe set them to zero
         ALLOCATE(aveU2(jpi,jpj,jpk))
         ALLOCATE(aveV2(jpi,jpj,jpk))
         ALLOCATE(aveW2(jpi,jpj,jpk))
         ALLOCATE(aveT2(jpi,jpj,jpk))
         ALLOCATE(aveS2(jpi,jpj,jpk))
         ALLOCATE(ave_ssh2(jpi,jpj))
         ALLOCATE(ave_taux(jpi,jpj))
         ALLOCATE(ave_tauy(jpi,jpj))

         aveU2 = 0.0
         aveV2 = 0.0
         aveW2 = 0.0
         aveT2 = 0.0
         aveS2 = 0.0
         ave_ssh2 = 0.0
         ave_taux = 0.0
         ave_tauy = 0.0

!!DBG
         if(lwp) then 
            write(numout2,*)'DBG: diawri(...aveTSUV()) opens: ',fnameT
         endif
      
! On first write step only
         CALL ncdf_write(fnameT, 'depthu', gdept, lncdf_stat) 
         CALL ncdf_write(fnameT, 'depthw', gdepw, lncdf_stat) 
         CALL ncdf_write(fnameT, 'nav_lat', gphit, -1, lncdf_stat)
         CALL ncdf_write(fnameT, 'nav_lon', glamt, -1, lncdf_stat)

      END IF  ! first timestep

!!DB: 2007.12.11 -- Code to compute ave fields
      aveU2 = aveU2 + un
      aveV2 = aveV2 + vn
      aveW2 = aveW2 + wn
      aveT2 = aveT2 + tn 
      aveS2 = aveS2 + sn 
      ave_ssh2 = ave_ssh2 + sshn 
      ave_taux = ave_taux + taux
      ave_tauy = ave_tauy + tauy

!!Output time if:
      if(mod(kt-nit000+1,ioutput_ave) ==  0) then
         rec_num2 = rec_num2 + 1
         aveU2 = aveU2/float(ioutput_ave)
         aveV2 = aveV2/float(ioutput_ave)
         aveW2 = aveW2/float(ioutput_ave)
         aveT2 = aveT2/float(ioutput_ave)
         aveS2 = aveS2/float(ioutput_ave)
         ave_ssh2 = ave_ssh2/float(ioutput_ave)
         ave_taux = ave_taux/float(ioutput_ave)
         ave_tauy = ave_tauy/float(ioutput_ave)

         CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'vozocrtx', aveU2, -rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'vomecrty', aveV2, -rec_num2, lncdf_stat)
!!W
         CALL ncdf_write(fnameT, 'vovecrtz', aveW2, -rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'votemper', aveT2, -rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'vosaline', aveS2, -rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'sossheig', ave_ssh2, -rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'sozotaux', ave_taux, -rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'sometauy', ave_tauy, -rec_num2, lncdf_stat)
!!AD/DB 2009.09.30
         CALL ncdf_write(fnameT, 'ndastp',REAL(ndastp), rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'model_time_step',REAL(kt), rec_num2, lncdf_stat)
         CALL ncdf_write(fnameT, 'model_time',model_time, rec_num2, lncdf_stat)

         aveU2 = 0.0
         aveV2 = 0.0
         aveW2 = 0.0
         aveT2 = 0.0
         aveS2 = 0.0
         ave_ssh2 = 0.0
         ave_taux = 0.0
         ave_tauy = 0.0

      endif
      


   END SUBROUTINE output_aveTSUV


END MODULE diawri
