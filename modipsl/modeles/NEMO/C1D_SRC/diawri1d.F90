MODULE diawri1d
   !!======================================================================
   !!                     ***  MODULE  diawri1d  ***
   !! Ocean diagnostics :  write ocean output files
   !!=====================================================================
#if defined key_cfg_1d
   !!----------------------------------------------------------------------
   !!   'key_cfg_1d'               1D Configuration
   !!----------------------------------------------------------------------  
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE zdf_oce         ! ocean vertical physics
   USE zdftke          ! TKE vertical mixing
   USE zdfkpp          ! KPP vertical mixing
   USE sol_oce         ! solver variables
   USE ice_oce         ! ice variables
   USE phycst          ! physical constants
   USE ocfzpt          ! ???
   USE ocesbc          ! surface thermohaline fluxes
   USE taumod          ! surface stress
   USE flxrnf          ! ???
   USE zdfmxl          ! mixed layer
   USE daymod          ! calendar
   USE dianam          ! build name of file (routine)
   USE diawri
   USE zdfddm          ! vertical  physics: double diffusion
   USE diahth          ! thermocline diagnostics
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dia_wri_1d                 ! routines called by step.F90
   !! * Module variables
   INTEGER ::   &
      nid_T, nz_T, nh_T, ndim_T, ndim_hT,      &   ! grid_T file
      ndex(1)                                      ! ???
   INTEGER, DIMENSION(jpi*jpj) ::   &
      ndex_hT
   INTEGER, DIMENSION(jpi*jpj*jpk) ::   &
      ndex_T

   !! * Substitutions
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/diawri1d.F90,v 1.5 2005/12/21 10:46:26 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS
   !!----------------------------------------------------------------------
   !!   Default option                                   NetCDF output file
   !!----------------------------------------------------------------------
   !!   dia_wri_1d       : create the standart NetCDF output files
   !!   dia_wri_state_1d : create an output NetCDF file for a single
   !!                      instantaeous ocean state and forcing fields
   !!----------------------------------------------------------------------

   SUBROUTINE dia_wri_1d( kt, kindic )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri_1d  ***
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
      !!   9.0  !  04-10  (C. Ethe)   1D Configuration
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   kindic  ! 

      !! * Local declarations
      LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
      CHARACTER (len=40) ::           &
         clhstnam, clop, clmx            ! temporary names
      INTEGER ::   inum = 11             ! temporary logical unit
      INTEGER ::   &
         ji, jj, ik                      ! dummy loop indices
      INTEGER ::   &
         iimi, iima, ipk, it,         &  ! temporary integers
         ijmi, ijma                      !    "          "
      REAL(wp) ::   &
         zsto, zout, zmax,            &  ! temporary scalars
         zjulian, zdt                    !    "         "
      REAL(wp), DIMENSION(jpi,jpj) :: &
         zw2d                            ! temporary workspace
      !!----------------------------------------------------------------------
      
      ! 0. Initialisation
      ! -----------------
      
      ! local variable for debugging
      ll_print = .FALSE.
      ll_print = ll_print .AND. lwp

      ! Define frequency of output and means
      zdt = rdt
      IF( nacc == 1 ) zdt = rdtmin
#if defined key_diainstant
      zsto = nwrite * zdt
      clop = "inst(x)"           ! no use of the mask value (require less cpu time)
      !!! clop="inst(only(x))"   ! put 1.e+20 on land (very expensive!!)
#else
      zsto=zdt
      clop="ave(x)"              ! no use of the mask value (require less cpu time)
      !!! clop="ave(only(x))"    ! put 1.e+20 on land (very expensive!!)
#endif
      zout = nwrite * zdt
      zmax = ( nitend - nit000 + 1 ) * zdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      it = kt - nit000 + 1


      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      IF(ll_print) WRITE(numout,*) 'dia_wri_1d kt = ', kt, ' kindic ', kindic

      IF( kt == nit000 ) THEN

         ! Define the NETCDF files (one per grid)
         
         ! Compute julian date from starting date of the run
         CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'Date 0 used :', nit000, ' YEAR ', nyear,   &
            &                    ' MONTH ', nmonth, ' DAY ', nday, 'Julian day : ', zjulian
         IF(lwp)WRITE(numout,*) ' indexes of zoom = ', iimi, iima, ijmi, ijma,   &
                                 ' limit storage in depth = ', ipk

         ! WRITE root name in date.file for use by postpro
         CALL dia_nam( clhstnam, nwrite,' ' )
         CALL ctlopn( inum, 'date.file', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
         WRITE(inum,*) clhstnam
         CLOSE(inum)
         
         ! Define the T grid FILE ( nid_T )
         
         CALL dia_nam( clhstnam, nwrite, 'grid_T' )
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          0, zjulian, zdt, nh_T, nid_T, domain_id=nidom )
         CALL histvert( nid_T, "deptht", "Vertical T levels",      &  ! Vertical grid: gdept
            &           "m", ipk, gdept, nz_T )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, tmask, 1, 1., ndex_T , ndim_T  )      ! volume
         CALL wheneq( jpi*jpj    , tmask, 1, 1., ndex_hT, ndim_hT )      ! surface


         ! Declare all the output fields as NETCDF variables

         !                                                                                      !!! nid_T : 3D
         CALL histdef( nid_T, "votemper", "Temperature"                        , "C"      ,   &  ! tn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         CALL histdef( nid_T, "vosaline", "Salinity"                           , "PSU"    ,   &  ! sn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !                                                                                      !!! nid_T : 2D
         CALL histdef( nid_T, "sosstsst", "Sea Surface temperature"            , "C"      ,   &  ! sst
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosaline", "Sea Surface Salinity"               , "PSU"    ,   &  ! sss
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )

#if ! defined key_dynspg_rl && defined key_ice_lim
         ! sowaflup = sowaflep + sorunoff + sowafldp + a term associated to
         !    internal damping to Levitus that can be diagnosed from others
         ! sowaflcd = sowaflep + sorunoff + sowafldp + iowaflup
         CALL histdef( nid_T, "iowaflup", "Ice=>ocean net freshwater"          , "kg/m2/s",   &  ! fsalt
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowaflep", "atmos=>ocean net freshwater"        , "kg/m2/s",   &  ! fmass
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif
         CALL histdef( nid_T, "sowaflup", "Net Upward Water Flux"              , "Kg/m2/s",   &  ! emp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sorunoff", "Runoffs"                            , "Kg/m2/s",   &  ! runoffs
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowaflcd", "concentration/dilution water flux"  , "kg/m2/s",   &  ! emps
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosalflx", "Surface Salt Flux"                  , "Kg/m2/s",   &  ! emps * sn
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sohefldo", "Net Downward Heat Flux"             , "W/m2"   ,   &  ! qt
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "soshfldo", "Shortwave Radiation"                , "W/m2"   ,   &  ! qsr
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "somxl010", "Mixed Layer Depth 0.01"             , "m"      ,   &  ! hmlp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#if defined key_zdfkpp
         CALL histdef( nid_T, "sokppekd", "Ekman depth                     "   , "m"      ,   &  ! sokppekd
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sokppbld", "Boundary Layer Depth            "   , "m"      ,   &  ! sokppbld
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif
         CALL histdef( nid_T, "somxlavt", "AVT : bottom of the mixed layer    ", "m"      ,   &  ! avt_mxl
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "somixhgt", "Turbocline Depth"                   , "m"      ,   &  ! hmld
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "soicecov", "Ice Cover"                          , "[0,1]"  ,   &  ! freeze
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#if ! defined key_coupled 
         CALL histdef( nid_T, "sohefldp", "Surface Heat Flux: Damping"         , "W/m2"   ,   &  ! qrp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowafldp", "Surface Water Flux: Damping"        , "Kg/m2/s",   &  ! erp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosafldp", "Surface salt flux: damping"         , "Kg/m2/s",   &  ! erp * sn
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif

#if ( defined key_coupled && ! defined key_ice_lim ) 
         CALL histdef( nid_T, "sohefldp", "Surface Heat Flux: Damping"         , "W/m2"   ,   &  ! qrp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowafldp", "Surface Water Flux: Damping"        , "Kg/m2/s",   &  ! erp
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosafldp", "Surface salt flux: Damping"         , "Kg/m2/s",   &  ! erp * sn
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif
         clmx ="l_max(only(x))"    ! max index on a period
         CALL histdef( nid_T, "sobowlin", "Bowl Index"                         , "W-point",   &  ! bowl INDEX 
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clmx, zsto, zout )
#if defined key_diahth
         CALL histdef( nid_T, "sothedep", "Thermocline Depth"                  , "m"      ,   & ! hth
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "so20chgt", "Depth of 20C isotherm"              , "m"      ,   & ! hd20
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "so28chgt", "Depth of 28C isotherm"              , "m"      ,   & ! hd28
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sohtc300", "Heat content 300 m"                 , "W"      ,   & ! htc3
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif

#if defined key_ice_lim && defined key_coupled
         CALL histdef( nid_T,"soicetem" , "Ice Surface Temperature"            , "K"      ,   &  ! tn_ice
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T,"soicealb" , "Ice Albedo"                         , "[0,1]"  ,   &  ! alb_ice
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif 

         !                                                                                      !!! nid_U : 3D
         CALL histdef( nid_T, "vozocrtx", "Zonal Current"                      , "m/s"    ,   &  ! un
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
#if defined key_diaeiv
         CALL histdef( nid_T, "vozoeivu", "Zonal EIV Current"                  , "m/s"    ,   &  ! u_eiv
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
#endif
         !                                                                                      !!! nid_U : 2D
         CALL histdef( nid_T, "sozotaux", "Wind Stress along i-axis"           , "N/m2"   ,   &  ! taux
            &          jpi, jpj, nh_T, 1  , 1, 1  , - 99, 32, clop, zsto, zout )

         !                                                                                      !!! nid_V : 3D
         CALL histdef( nid_T, "vomecrty", "Meridional Current"                 , "m/s"    ,   &  ! vn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
#if defined key_diaeiv
         CALL histdef( nid_T, "vomeeivv", "Meridional EIV Current"             , "m/s"    ,   &  ! v_eiv
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
#endif
         !                                                                                      !!! nid_V : 2D
         CALL histdef( nid_T, "sometauy", "Wind Stress along j-axis"           , "N/m2"   ,   &  ! tauy
            &          jpi, jpj, nh_T, 1  , 1, 1  , - 99, 32, clop, zsto, zout )
#if defined key_zdftke
         !                                                                                      !!! nid_W : 3D
         CALL histdef( nid_T, "votlsdis", " Dissipation Turbulent Lenght Scale", "m"      ,   &  ! e_dis
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !
         CALL histdef( nid_T, "votlsmix", " Mixing Turbulent Lenght Scale"     , "m"      ,   &  ! e_mix
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !
         CALL histdef( nid_T, "votlspdl", " Prandl Number",                      "-"       ,   &  ! e_pdl
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !
         CALL histdef( nid_T, "votlsric", " Local Richardson Number",            "-"       ,   &  ! e_ric
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !
         CALL histdef( nid_T, "votkeend", "TKE: Turbulent kinetic energy"       , "m2/s"   ,   &  ! TKE
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
#endif
#if defined key_zdfkpp
         !                                                                                      !!! nid_W : 3D
         CALL histdef( nid_T, "vokpprig", " Gradient Richardson Number"        ,  "-"      ,   &  ! rig
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !
         CALL histdef( nid_T, "vokpprib", " Bulk Richardson Number    "        ,  "-"      ,   &   ! rib
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !
         CALL histdef( nid_T, "vokppbsf", " Buoyancy forcing          "        , "N/m2"    ,   &  ! sokppbsf
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !
         CALL histdef( nid_T, "vokppmol", "Moning Obukhov length scale     "   , "m"       ,   &  ! sokppmol
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
#endif
         !
         CALL histdef( nid_T, "voeosbn2", "Brunt-Vaisala Frequency"             , "m2/s2"  ,   &  ! rn2
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )

         CALL histdef( nid_T, "votkeavt", "Vertical Eddy Diffusivity"          , "m2/s"   ,   &  ! avt
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )

         CALL histdef( nid_T, "votkeevd", "Enhanced Vertical Diffusivity",       "m2/s"   ,   &  ! avt_evd
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )

         CALL histdef( nid_T, "votkeavm", "Vertical Eddy Viscosity",             "m2/s"   ,   &  ! avmu
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         !
         CALL histdef( nid_T, "votkeevm", "Enhanced Vertical Viscosity",         "m2/s"   ,   &  ! avmu_evd
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )

         IF( lk_zdfddm ) THEN
            CALL histdef( nid_T,"voddmavs","Salt Vertical Eddy Diffusivity"    , "m2/s"   ,   &  ! avs
               &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         ENDIF

         CALL histend( nid_T )

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization'
         IF(ll_print) CALL FLUSH(numout )

      ENDIF

      ! 2. Start writing data
      ! ---------------------

      ! ndex(1) est utilise ssi l'avant dernier argument est diffferent de 
      ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
      ! donne le nombre d'elements, et ndex la liste des indices a sortir

      IF( lwp .AND. MOD( kt, nwrite ) == 0 ) THEN 
         WRITE(numout,*) 'dia_wri : write model outputs in NetCDF files at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF

      ! Write fields on T grid
      CALL histwrite( nid_T, "votemper", it, tn            , ndim_T , ndex_T  )   ! temperature
      CALL histwrite( nid_T, "vosaline", it, sn            , ndim_T , ndex_T  )   ! salinity
      CALL histwrite( nid_T, "sosstsst", it, tn(:,:,1)     , ndim_hT, ndex_hT )   ! sea surface temperature
      CALL histwrite( nid_T, "sosaline", it, sn(:,:,1)     , ndim_hT, ndex_hT )   ! sea surface salinity
#if ! defined key_dynspg_rl && defined key_ice_lim
      CALL histwrite( nid_T, "iowaflup", it, fsalt(:,:)    , ndim_hT, ndex_hT )   ! ice=>ocean water flux
      CALL histwrite( nid_T, "sowaflep", it, fmass(:,:)    , ndim_hT, ndex_hT )   ! atmos=>ocean water flux
#endif
      CALL histwrite( nid_T, "sowaflup", it, emp           , ndim_hT, ndex_hT )   ! upward water flux
      CALL histwrite( nid_T, "sorunoff", it, runoff        , ndim_hT, ndex_hT )   ! runoff
      CALL histwrite( nid_T, "sowaflcd", it, emps          , ndim_hT, ndex_hT )   ! c/d water flux
      zw2d(:,:) = emps(:,:) * sn(:,:,1) * tmask(:,:,1)
      CALL histwrite( nid_T, "sosalflx", it, zw2d          , ndim_hT, ndex_hT )   ! c/d salt flux
      CALL histwrite( nid_T, "sohefldo", it, qt            , ndim_hT, ndex_hT )   ! total heat flux
      CALL histwrite( nid_T, "soshfldo", it, qsr           , ndim_hT, ndex_hT )   ! solar heat flux
      CALL histwrite( nid_T, "somxl010", it, hmlp          , ndim_hT, ndex_hT )   ! mixed layer depth
#if defined key_zdfkpp
      CALL histwrite( nid_T, "sokppekd", it, ekdp          , ndim_hT, ndex_hT )   ! Ekman depht
      CALL histwrite( nid_T, "sokppbld", it, hkpp          , ndim_hT, ndex_hT )   ! boundary layer depth 
#endif  
      ! store the vertical eddy diffusivity coef. at the bottom of the mixed layer
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = nmln(ji,jj)
            zw2d(ji,jj) = avt(ji,jj,ik) * tmask(ji,jj,1)
         END DO
      END DO
      CALL histwrite( nid_T, "somxlavt", it, zw2d          , ndim_hT, ndex_hT )   ! Kz at bottom of mixed layer 
      CALL histwrite( nid_T, "somixhgt", it, hmld          , ndim_hT, ndex_hT )   ! turbocline depth
      CALL histwrite( nid_T, "soicecov", it, freeze        , ndim_hT, ndex_hT )   ! ice cover 
#if ! defined key_coupled
      CALL histwrite( nid_T, "sohefldp", it, qrp           , ndim_hT, ndex_hT )   ! heat flux damping
      CALL histwrite( nid_T, "sowafldp", it, erp           , ndim_hT, ndex_hT )   ! freshwater flux damping
      zw2d(:,:) = erp(:,:) * sn(:,:,1) * tmask(:,:,1)
      CALL histwrite( nid_T, "sosafldp", it, zw2d          , ndim_hT, ndex_hT )   ! salt flux damping
#endif
#if ( defined key_coupled && ! defined key_ice_lim ) 
      CALL histwrite( nid_T, "sohefldp", it, qrp           , ndim_hT, ndex_hT )   ! heat flux damping
      CALL histwrite( nid_T, "sowafldp", it, erp           , ndim_hT, ndex_hT )   ! freshwater flux damping
         zw2d(:,:) = erp(:,:) * sn(:,:,1) * tmask(:,:,1)
      CALL histwrite( nid_T, "sosafldp", it, zw2d          , ndim_hT, ndex_hT )   ! salt flux damping
#endif
         zw2d(:,:) = FLOAT( nmln(:,:) ) * tmask(:,:,1)
      CALL histwrite( nid_T, "sobowlin", it, zw2d          , ndim_hT, ndex_hT )   ! ???

#if defined key_diahth
      CALL histwrite( nid_T, "sothedep", it, hth           , ndim_hT, ndex_hT )   ! depth of the thermocline
      CALL histwrite( nid_T, "so20chgt", it, hd20          , ndim_hT, ndex_hT )   ! depth of the 20 isotherm
      CALL histwrite( nid_T, "so28chgt", it, hd28          , ndim_hT, ndex_hT )   ! depth of the 28 isotherm
      CALL histwrite( nid_T, "sohtc300", it, htc3          , ndim_hT, ndex_hT )   ! first 300m heaat content
#endif
#if defined key_ice_lim &&  defined key_coupled 
      CALL histwrite( nid_T, "soicetem", it, tn_ice        , ndim_hT, ndex_hT )   ! surf. ice temperature
      CALL histwrite( nid_T, "soicealb", it, alb_ice       , ndim_hT, ndex_hT )   ! ice albedo
#endif

      CALL histwrite( nid_T, "vozocrtx", it, un            , ndim_T , ndex_T )    ! i-current
      CALL histwrite( nid_T, "sozotaux", it, taux          , ndim_hT, ndex_hT )   ! i-wind stress
      CALL histwrite( nid_T, "vomecrty", it, vn            , ndim_T , ndex_T  )   ! j-current
      CALL histwrite( nid_T, "sometauy", it, tauy          , ndim_hT, ndex_hT )   ! j-wind stress
#if defined key_zdftke
      CALL histwrite( nid_T, "votlsdis", it, e_dis         , ndim_T , ndex_T )    ! Diss. Turb. lenght scale
      CALL histwrite( nid_T, "votlsmix", it, e_mix         , ndim_T , ndex_T )    ! Mixing Turb. lenght scale
      CALL histwrite( nid_T, "votlspdl", it, e_pdl         , ndim_T , ndex_T )    ! Prandl number
      CALL histwrite( nid_T, "votlsric", it, e_ric         , ndim_T , ndex_T )    ! local Richardson number
      CALL histwrite( nid_T, "votkeend", it, en            , ndim_T , ndex_T )    ! TKE
#endif
#if defined key_zdfkpp
      CALL histwrite( nid_T, "vokpprig", it, rig           , ndim_T , ndex_T )    ! gradient Richardson number
      CALL histwrite( nid_T, "vokpprib", it, rib           , ndim_T , ndex_T )    ! bulk Richardson number
      CALL histwrite( nid_T, "vokppbsf", it, buof          , ndim_T , ndex_T )    ! buoyancy forcing
      CALL histwrite( nid_T, "vokppmol", it, mols          , ndim_T , ndex_T )    ! Moning-Obukov length scale
#endif
      CALL histwrite( nid_T, "voeosbn2", it, rn2           , ndim_T , ndex_T )    ! Brunt-Vaisala Frequency
      CALL histwrite( nid_T, "votkeavt", it, avt           , ndim_T , ndex_T )    ! T vert. eddy diff. coef.
      CALL histwrite( nid_T, "votkeevd", it, avt_evd       , ndim_T , ndex_T )    ! T enhan. vert. eddy diff. coef.
      CALL histwrite( nid_T, "votkeavm", it, avmu          , ndim_T , ndex_T )    ! T vert. eddy visc. coef.
      CALL histwrite( nid_T, "votkeevm", it, avmu_evd      , ndim_T , ndex_T )    ! T enhan. vert. eddy visc. coef.
      IF( lk_zdfddm ) THEN
         CALL histwrite( nid_T, "voddmavs", it, fsavs(:,:,:), ndim_T, ndex_T )    ! S vert. eddy diff. coef.
      ENDIF

      ! 3. Synchronise and close all files
      ! ---------------------------------------
      IF( MOD( kt, nwrite ) == 0 .OR. kindic < 0 ) THEN
         CALL histsync( nid_T )
      ENDIF

      !  Create an output files (output.abort.nc) if S < 0 or u > 20 m/s
      IF( kindic < 0 )   CALL dia_wri_state( 'output.abort' )

      IF( kt == nitend .OR. kindic < 0 ) THEN
         CALL histclo( nid_T )
      ENDIF

   END SUBROUTINE dia_wri_1d
#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO 1D Config
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dia_wri_1d ( kt, kindic )
      if(lwp) WRITE(numout,*) 'dia_wri_1d: You should not have seen this print! error?', kt, kindic
   END SUBROUTINE dia_wri_1d
#endif

   !!======================================================================
END MODULE diawri1d
