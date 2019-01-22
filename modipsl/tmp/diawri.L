


PGF90 (Version     10.5)          09/23/2010  09:49:30      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir /usr/local/include
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /usr/local/include
          -idir /usr/include/mpich2-x86_64
          -idir /usr/include/mpich2-x86_64
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: diawri.F90

(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) !!DB - 2007.12.04
(    3) !!This version has eliminated the IOIPSL code
(    4) 
(    5) 
(    6) MODULE diawri
(    7)    !!======================================================================
(    8)    !!                     ***  MODULE  diawri  ***
(    9)    !! Ocean diagnostics :  write ocean output files
(   10)    !!=====================================================================
(   11) 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers 
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE zdf_oce         ! ocean vertical physics
(   17)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   18)    USE ldfdyn_oce      ! ocean dynamics: lateral physics
(   19)    USE sol_oce         ! solver variables
(   20)    USE ice_oce         ! ice variables
(   21)    USE phycst          ! physical constants
(   22)    USE ocfzpt          ! ocean freezing point
(   23)    USE ocesbc          ! surface thermohaline fluxes
(   24)    USE taumod          ! surface stress
(   25)    USE flxrnf          ! ocean runoffs
(   26)    USE zdfmxl          ! mixed layer
(   27)    USE daymod          ! calendar
(   28)    USE dianam          ! build name of file (routine)
(   29)    USE zdfddm          ! vertical  physics: double diffusion
(   30)    USE diahth          ! thermocline diagnostics
(   31)    USE diaspr          ! surface pressure diagnostics (rigid lid case)
(   32)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   33)    USE in_out_manager  ! I/O manager
(   34)    USE flx_oce         ! sea-ice/ocean forcings variables
(   35) 
(   36)    USE lib_ncdf        ! netCDF I/O library
(   37) 
(   38)    IMPLICIT NONE
(   39)    PRIVATE
(   40) 
(   41)    !! * Accessibility
(   42)    PUBLIC dia_wri                 ! routines called by step.F90
(   43)    PUBLIC dia_wri_state
(   44) !!DB: tried to move to lib_ncdf but did not work (???)
(   45)    PUBLIC output_special             ! special routines called by step.F90 if ave flag is on
(   46)    PUBLIC output_aveTSUV             ! special routines called by step.F90 if ave flag is on






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 2

(   47) 
(   48)    !! * Module variables
(   49) !   INTEGER ::   &
(   50) !      nid_T, nz_T, nh_T, ndim_T, ndim_hT,      &   ! grid_T file
(   51) !      nid_U, nz_U, nh_U, ndim_U, ndim_hU,      &   ! grid_U file
(   52) !      nid_V, nz_V, nh_V, ndim_V, ndim_hV,      &   ! grid_V file
(   53) !      nid_W, nz_W, nh_W,                       &   ! grid_W file
(   54) !      ndex(1)                                      ! ???
(   55) !   INTEGER, DIMENSION(jpi*jpj) ::   &
(   56) !      ndex_hT, ndex_hU, ndex_hV
(   57) !   INTEGER, DIMENSION(jpi*jpj*jpk) ::   &
(   58) !      ndex_T, ndex_U, ndex_V
(   59) 
(   60)    !! * Substitutions
(   61) #  include "zdfddm_substitute.h90"
(   62)    !!----------------------------------------------------------------------
(   63)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   64)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diawri.F90,v 1.10 2006/03/09 17:21:54 opalod Exp $ 
(   65)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   66)    !!----------------------------------------------------------------------
(   67) 
(   68) CONTAINS
(   69) 
(   70)    !!----------------------------------------------------------------------
(   71)    !!   Default option                                   NetCDF output file
(   72)    !!----------------------------------------------------------------------
(   73)    !!   dia_wri       : create the standart NetCDF output files
(   74)    !!   dia_wri_state : create an output NetCDF file for a single
(   75)    !!                   instantaeous ocean state and forcing fields
(   76)    !!----------------------------------------------------------------------
(   77) 
(   78)    SUBROUTINE dia_wri( kt, kindic )
(   79)       !!---------------------------------------------------------------------
(   80)       !!                  ***  ROUTINE dia_wri  ***
(   81)       !!                   
(   82)       !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
(   83)       !!      NETCDF format is used by default 
(   84)       !!
(   85)       !! ** Method  :   At the beginning of the first time step (nit000), 
(   86)       !!      define all the NETCDF files and fields
(   87)       !!      At each time step call histdef to compute the mean if ncessary
(   88)       !!      Each nwrite time step, output the instantaneous or mean fields
(   89)       !!      IF kindic <0, output of fields before the model interruption.
(   90)       !!      IF kindic =0, time step loop
(   91)       !!      IF kindic >0, output of fields before the time step loop
(   92)       !!
(   93)       !! History :
(   94)       !!        !  91-03  (M.-A. Foujols)  Original code
(   95)       !!        !  91-11  (G. Madec)
(   96)       !!        !  92-06  (M. Imbard)  correction restart file
(   97)       !!        !  92-07  (M. Imbard)  split into diawri and rstwri
(   98)       !!        !  93-03  (M. Imbard)  suppress writibm
(   99)       !!        !  98-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
(  100)       !!        !  99-02  (E. Guilyardi)  name of netCDF files + variables
(  101)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  102)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  103)       !!----------------------------------------------------------------------
(  104)       !! * Modules used






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 3

(  105) !      USE ioipsl
(  106) 
(  107)       !! * Arguments
(  108)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  109)       INTEGER, INTENT( in ) ::   kindic  ! 
(  110) 
(  111)       !! * Local declarations
(  112)       LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
(  113)       CHARACTER (len=40) ::           &
(  114)          clhstnam, clop, clmx            ! temporary names
(  115)       INTEGER ::   inum = 11             ! temporary logical unit
(  116)       INTEGER ::   &
(  117)          iimi, iima, ipk, it,         &  ! temporary integers
(  118)          ijmi, ijma                      !    "          "
(  119) 
(  120)       INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls
(  121)       LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.    ! Use IOIPSL subroutines for output
(  122)                                                  ! Note that all W-grid calls are commented out!
(  123)       LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE.   ! Use lib_ncdf subroutines for output
(  124) 
(  125)       REAL(wp) :: test_array
(  126)       INTEGER :: a, b, c
(  127) 
(  128)       REAL(wp) ::   &
(  129)          zsto, zout, zmax,            &  ! temporary scalars
(  130)          zjulian, zdt                    !    "         "
(  131)       REAL(wp), DIMENSION(jpi,jpj) :: &
(  132)          zw2d                            ! temporary workspace
(  133)       CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT
(  134) !DB
(  135)       INTEGER :: tindex
(  136) 
(  137) 
(  138)       !!----------------------------------------------------------------------
(  139)       
(  140)       ! 0. Initialisation
(  141)       ! -----------------
(  142)       fnameU = trim(cexper)//'_grid_U.nc'
(  143)       fnameV = trim(cexper)//'_grid_V.nc'
(  144)       fnameT = trim(cexper)//'_grid_T.nc'
(  145) 
(  146) !!DB deleted ...
(  147) !      IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL for output?
(  148) 
(  149)    ! LIB_NCDF CALLS
(  150) 
(  151) !DB default to this code
(  152) !   IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Are we using lib_ncdf output routines?
(  153) 
(  154)       ! On the first timestep, create output files and write depth variables (since they don't change)
(  155)       IF(kt == nit000) THEN
(  156) !!DBG
(  157)          write(numout2,*)'DBG: diawri: cexper = ', cexper
(  158)          write(numout2,*)'DBG: diawri: ',fnameU, fnameV, fnameT
(  159) 
(  160)          CALL ncdf_create_file_u(fnameU, 'inst(x)', lncdf_stat)
(  161)          CALL ncdf_create_file_v(fnameV, 'inst(x)', lncdf_stat)
(  162)          CALL ncdf_create_file_t(fnameT, 'inst(x)', lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 4

(  163)          CALL ncdf_write(fnameU, 'depthu', gdept, lncdf_stat) 
(  164)          CALL ncdf_write(fnameV, 'depthv', gdept, lncdf_stat) 
(  165)          CALL ncdf_write(fnameT, 'deptht', gdept, lncdf_stat) 
(  166)       END IF
(  167)       
(  168)       ! On first write step only, write nav_lat and nav_lon grids (routine will abort if not on a write step)
(  169)       IF(kt == (nit000+ nwrite)) THEN
(  170)          CALL ncdf_write(fnameU, 'nav_lat', gphiu, kt, lncdf_stat)
(  171)          CALL ncdf_write(fnameU, 'nav_lon', glamu, kt, lncdf_stat)
(  172)          CALL ncdf_write(fnameV, 'nav_lat', gphiv, kt, lncdf_stat)
(  173)          CALL ncdf_write(fnameV, 'nav_lon', glamv, kt, lncdf_stat)
(  174)          CALL ncdf_write(fnameT, 'nav_lat', gphit, kt, lncdf_stat)
(  175)          CALL ncdf_write(fnameT, 'nav_lon', glamt, kt, lncdf_stat)
(  176)       END IF
(  177)       
(  178)       ! On write steps, output time counter
(  179)       IF(MOD(kt, nwrite) == 0) THEN
(  180)          tindex = kt/nwrite - nit000/nwrite
(  181) !         CALL ncdf_write(fnameU, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
(  182) !         CALL ncdf_write(fnameV, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
(  183) !         CALL ncdf_write(fnameT, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
(  184)          CALL ncdf_write(fnameU, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
(  185)          CALL ncdf_write(fnameV, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
(  186)          CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
(  187)       END IF
(  188)       
(  189)       ! Write other fields (calls will abort if this isn't a write step)
(  190)       CALL ncdf_write(fnameU, 'vozocrtx', un, kt, lncdf_stat)
(  191)       CALL ncdf_write(fnameU, 'sozotaux', taux, kt, lncdf_stat)
(  192)       CALL ncdf_write(fnameV, 'vomecrty', vn, kt, lncdf_stat)
(  193)       CALL ncdf_write(fnameV, 'sometauy', tauy, kt, lncdf_stat)
(  194)       CALL ncdf_write(fnameT, 'votemper', tn, kt, lncdf_stat)
(  195)       CALL ncdf_write(fnameT, 'vosaline', sn, kt, lncdf_stat)
(  196)       CALL ncdf_write(fnameT, 'sosstsst', tn(:,:,1), kt, lncdf_stat)
(  197)       CALL ncdf_write(fnameT, 'sosaline', sn(:,:,1), kt, lncdf_stat)
(  198)       CALL ncdf_write(fnameT, 'sossheig', sshn, kt, lncdf_stat)
(  199)       CALL ncdf_write(fnameT, 'sowaflup', emp, kt, lncdf_stat)
(  200)       CALL ncdf_write(fnameT, 'sorunoff', runoff, kt, lncdf_stat)
(  201)       CALL ncdf_write(fnameT, 'sowaflcd', emps, kt, lncdf_stat)
(  202)       zw2d(:,:) = emps(:,:) * sn(:,:,1) * tmask(:,:,1)
(  203)       CALL ncdf_write(fnameT, 'sosalflx', zw2d, kt, lncdf_stat)
(  204)       CALL ncdf_write(fnameT, 'sohefldo', qt, kt, lncdf_stat)
(  205)       CALL ncdf_write(fnameT, 'soshfldo', qsr, kt, lncdf_stat)
(  206)       CALL ncdf_write(fnameT, 'somxl010', hmlp, kt, lncdf_stat)
(  207)       CALL ncdf_write(fnameT, 'somixhgt', hmld, kt, lncdf_stat)
(  208)       CALL ncdf_write(fnameT, 'soicecov', freeze, kt, lncdf_stat)
(  209)       CALL ncdf_write(fnameT, 'sohefldp', qrp, kt, lncdf_stat)
(  210)       CALL ncdf_write(fnameT, 'sowafldp', erp, kt, lncdf_stat)
(  211)       zw2d(:,:) = erp(:,:) * sn(:,:,1) * tmask(:,:,1)
(  212)       CALL ncdf_write(fnameT, 'sosafldp', zw2d, kt, lncdf_stat)
(  213)       zw2d(:,:) = FLOAT( nmln(:,:) ) * tmask(:,:,1)
(  214)       CALL ncdf_write(fnameT, 'sobowlin', zw2d, kt, lncdf_stat)
(  215)       
(  216)    ! END OF LIB_NCDF CALLS
(  217) 
(  218) !DB code deleted
(  219) !   IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL output subroutines?
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 5

(  221)    END SUBROUTINE dia_wri
(  222) 
(  223) 
(  224)    SUBROUTINE dia_wri_state( cdfile_name )
(  225)       !!---------------------------------------------------------------------
(  226)       !!                 ***  ROUTINE dia_wri_state  ***
(  227)       !!        
(  228)       !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
(  229)       !!      the instantaneous ocean state and forcing fields.
(  230)       !!        Used to find errors in the initial state or save the last
(  231)       !!      ocean state in case of abnormal end of a simulation
(  232)       !!
(  233)       !! ** Method  :   NetCDF files using ioipsl
(  234)       !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
(  235)       !!      File 'output.abort.nc' is created in case of abnormal job end
(  236)       !!
(  237)       !! History :
(  238)       !!   8.2  !  00-06  (M. Imbard)  Original code (diabort.F)
(  239)       !!   8.5  !  02-06  (A.Bozec, E. Durand)  Original code (diainit.F)
(  240)       !!   9.0  !  02-12  (G. Madec)  merge of diabort and diainit, F90
(  241)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  242)       !!----------------------------------------------------------------------
(  243)       !! * Modules used
(  244)       USE ioipsl
(  245) 
(  246)       !! * Arguments
(  247)       CHARACTER (len=* ), INTENT( in ) ::   &
(  248)          cdfile_name      ! name of the file created
(  249) 
(  250)       !! * Local declarations
(  251)       CHARACTER (len=40) :: clop
(  252)       INTEGER  ::   &
(  253)          id_i , nz_i, nh_i       
(  254)       INTEGER, DIMENSION(1) ::   &
(  255)          idex             ! temprary workspace
(  256)       REAL(wp) ::   &
(  257)          zsto, zout, zmax,   &
(  258)          zjulian, zdt
(  259)       !!----------------------------------------------------------------------
(  260) 
(  261)       IF(lwp) WRITE(numout,*)
(  262)       IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
(  263)       IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
(  264)       IF(lwp) WRITE(numout,*) '                and named :', cdfile_name, '.nc'
(  265)       
(  266)       ! 0. Initialisation
(  267)       ! -----------------
(  268)       
(  269)       ! Define frequency of output and means
(  270)       zdt  = rdt
(  271)       zsto = rdt
(  272)       clop = "inst(x)"           ! no use of the mask value (require less cpu time)
(  273)       zout = rdt
(  274)       zmax = ( nitend - nit000 + 1 ) * zdt
(  275) 
(  276)       ! 1. Define NETCDF files and fields at beginning of first time step
(  277)       ! -----------------------------------------------------------------
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 6

(  279)       ! Compute julian date from starting date of the run
(  280)       CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )         ! time axis 
(  281)       CALL histbeg( cdfile_name, jpi, glamt, jpj, gphit,   &
(  282)           1, jpi, 1, jpj, 0, zjulian, zdt, nh_i, id_i, domain_id=nidom )          ! Horizontal grid : glamt and gphit
(  283)       CALL histvert( id_i, "deptht", "Vertical T levels",   &    ! Vertical grid : gdept
(  284)           "m", jpk, gdept, nz_i)
(  285) 
(  286)       ! Declare all the output fields as NetCDF variables
(  287) 
(  288)       CALL histdef( id_i, "vosaline", "Salinity"              , "PSU"    ,   &   ! salinity
(  289)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
(  290)       CALL histdef( id_i, "votemper", "Temperature"           , "C"      ,   &   ! temperature
(  291)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
(  292) #if defined key_dynspg_rl
(  293)       CALL histdef( id_i, "sobarstf","Barotropic StreamFunction", "m3/s2"  ,   &  ! bsf
(  294)          &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
(  295) #else
(  296)       CALL histdef( id_i, "sossheig", "Sea Surface Height"    , "m"      ,   &  ! ssh
(  297)          &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
(  298) #endif
(  299)       CALL histdef( id_i, "vozocrtx", "Zonal Current"         , "m/s"    ,   &   ! zonal current
(  300)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
(  301)       CALL histdef( id_i, "vomecrty", "Meridional Current"    , "m/s"    ,   &   ! meridonal current
(  302)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
(  303)       CALL histdef( id_i, "vovecrtz", "Vertical Velocity"     , "m/s"    ,   &   ! vertical current
(  304)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
(  305)       CALL histdef( id_i, "sowaflup", "Net Upward Water Flux" , "Kg/m2/S",   &   ! net freshwater 
(  306)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  307)       CALL histdef( id_i, "sohefldo", "Net Downward Heat Flux", "W/m2"   ,   &   ! net heat flux
(  308)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  309)       CALL histdef( id_i, "soshfldo", "Shortwave Radiation"   , "W/m2"   ,   &   ! solar flux
(  310)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  311)       CALL histdef( id_i, "soicecov", "Ice fraction"          , "[0,1]"  ,   &   ! freeze
(  312)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  313)       CALL histdef( id_i, "sozotaux", "Zonal Wind Stress"     , "N/m2"   ,   &   ! i-wind stress
(  314)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  315)       CALL histdef( id_i, "sometauy", "Meridional Wind Stress", "N/m2"   ,   &   ! j-wind stress
(  316)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  317) 
(  318)       CALL histend( id_i )
(  319) 
(  320)       ! 2. Start writing data
(  321)       ! ---------------------
(  322)       ! idex(1) est utilise ssi l'avant dernier argument est diffferent de 
(  323)       ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
(  324)       ! donne le nombre d'elements, et idex la liste des indices a sortir
(  325)       idex(1) = 1   ! init to avoid compil warning
(  326)       
(  327)       ! Write all fields on T grid
(  328)       CALL histwrite( id_i, "votemper", 1, tn    , jpi*jpj*jpk, idex )    ! now temperature
(  329)       CALL histwrite( id_i, "vosaline", 1, sn    , jpi*jpj*jpk, idex )    ! now salinity
(  330) #if defined key_dynspg_rl
(  331)       CALL histwrite( id_i, "sobarstf", 1, bsfn  , jpi*jpj    , idex )    ! barotropic streamfunction
(  332) #else
(  333)       CALL histwrite( id_i, "sossheig", 1, sshn  , jpi*jpj    , idex )    ! sea surface height
(  334) #endif
(  335)       CALL histwrite( id_i, "vozocrtx", 1, un    , jpi*jpj*jpk, idex )    ! now i-velocity
(  336)       CALL histwrite( id_i, "vomecrty", 1, vn    , jpi*jpj*jpk, idex )    ! now j-velocity






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 7

(  337)       CALL histwrite( id_i, "vovecrtz", 1, wn    , jpi*jpj*jpk, idex )    ! now k-velocity
(  338)       CALL histwrite( id_i, "sowaflup", 1, emp   , jpi*jpj    , idex )    ! freshwater budget
(  339)       CALL histwrite( id_i, "sohefldo", 1, qt    , jpi*jpj    , idex )    ! total heat flux
(  340)       CALL histwrite( id_i, "soshfldo", 1, qsr   , jpi*jpj    , idex )    ! total heat flux
(  341)       CALL histwrite( id_i, "soicecov", 1, freeze, jpi*jpj    , idex )    ! ice cover
(  342)       CALL histwrite( id_i, "sozotaux", 1, taux  , jpi*jpj    , idex )    ! i-wind stress
(  343)       CALL histwrite( id_i, "sometauy", 1, tauy  , jpi*jpj    , idex )    ! j-wind stress
(  344) 
(  345)       ! 3. Close the file
(  346)       ! -----------------
(  347)       CALL histclo( id_i )
(  348) 
(  349)    END SUBROUTINE dia_wri_state
(  350) 
(  351)    !!======================================================================
(  352) 
(  353) !!DB
(  354)    SUBROUTINE output_special( kt, kindic )
(  355)       !!---------------------------------------------------------------------
(  356)       !!                  ***  ROUTINE dia_wri_ave  ***
(  357)       !!                   
(  358)       !! ** Purpose :   Average and output special opa variable
(  359)       !!      NETCDF format is used by default 
(  360)       !!
(  361)       !! DB: 2007.12.04 ... 12.10 ... tested OK
(  362)       !!----------------------------------------------------------------------
(  363) !!DB -- notes
(  364) !1 - due to the way that lib_ncdf routine handles 'when to write' this routine will
(  365) ! have to fool the call to ncdf_write...() to write at a non-nwrite timestep:
(  366) ! fixed by inputting a -ve rec_num value
(  367) !2 - ALLOCATE arrays when routine is called ... DONE ...
(  368) !3 - Uses an improved M2-averaging routine (see dbrick:..MATLAB/test_M2ave.m)
(  369) ! that interpolates to the correct (float)time (rM2dt). Routine also gives correct
(  370) ! answer if rM2dt = integer
(  371) 
(  372)       !! * Modules used
(  373) !      USE ioipsl
(  374) !     USE oce
(  375) !     USE taumod
(  376) 
(  377)       !! * Arguments
(  378)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  379)       INTEGER, INTENT( in ) ::   kindic  ! 
(  380) 
(  381)       !! * Local declarations
(  382)       LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
(  383)       CHARACTER (len=40) ::           &
(  384)          clhstnam, clop, clmx            ! temporary names
(  385)       INTEGER ::   inum = 11             ! temporary logical unit
(  386)       INTEGER ::   &
(  387)          iimi, iima, ipk, it,         &  ! temporary integers
(  388)          ijmi, ijma                      !    "          "
(  389) 
(  390)       INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls
(  391) 
(  392)       REAL(wp) ::   &
(  393)          zsto, zout, zmax,            &  ! temporary scalars
(  394)          zjulian, zdt                    !    "         "






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 8

(  395)       REAL(wp), DIMENSION(jpi,jpj) :: &
(  396)          zw2d                            ! temporary workspace
(  397)       CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT
(  398) 
(  399) !!DB
(  400)       INTEGER :: kkk, M2dt, ndt, ddt, M2dt_lo, M2dt_hi
(  401)       REAL(wp) :: M2_hr, rM2dt, zxy
(  402)       REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: aveU, aveV, aveT
(  403)       REAL(wp),DIMENSION(:,:),ALLOCATABLE, SAVE :: ave_ssh
(  404)       INTEGER, SAVE :: rec_num
(  405) 
(  406)       ! 0. Initialisation
(  407) !!DB
(  408)       M2_hr = 12.42
(  409)       ndt = 86400./rdt
(  410)       rM2dt = M2_hr*3600./rdt
(  411)       M2dt_lo = floor(M2_hr*3600./rdt)
(  412) !      M2dt_hi = ceil(M2_hr*3600./rdt) !no ceil() but variable is not needed
(  413)       M2dt = M2dt_lo;
(  414)       zxy = rM2dt-float(M2dt)
(  415) 
(  416)       fnameU = trim(cexper)//'_M2ave_U.nc'
(  417)       fnameV = trim(cexper)//'_M2ave_V.nc'
(  418)       fnameT = trim(cexper)//'_M2ave_T.nc'
(  419) !      fnameT = trim(cexper)//'_ave_TSUV.nc'
(  420) 
(  421) ! On the first timestep, create output files and write depth variables (since they don't change)
(  422)       IF(kt == nit000) THEN
(  423)          rec_num = 0
(  424) 
(  425)          call ncdf_create_file_ave(fnameU, 'aveU', 'ave(x)', lncdf_stat)
(  426)          call ncdf_create_file_ave(fnameV, 'aveV', 'ave(x)', lncdf_stat)
(  427) !         call ncdf_create_file_ave(fnameT, 'aveT', 'ave(x)', lncdf_stat)
(  428) 
(  429)          ALLOCATE(aveU(jpi,jpj,jpk))
(  430)          ALLOCATE(aveV(jpi,jpj,jpk))
(  431) !         ALLOCATE(aveT(jpi,jpj,jpk))
(  432)          aveU = 0.0
(  433)          aveV = 0.0
(  434) 
(  435) !!DBG
(  436)          if(lwp) then
(  437)             write(numout2,*)'DBG: diawri: cexper = ', cexper
(  438)             write(numout2,*)'DBG: diawri: ',fnameU, fnameV
(  439)          endif
(  440) 
(  441) 
(  442)          CALL ncdf_write(fnameU, 'depthu', gdept, lncdf_stat) 
(  443)          CALL ncdf_write(fnameV, 'depthu', gdept, lncdf_stat) 
(  444) !         CALL ncdf_write(fnameT, 'deptht', gdept, lncdf_stat) 
(  445)       
(  446) ! On first write step only, write nav_lat and nav_lon grids (routine will abort if not on a write step)
(  447)          CALL ncdf_write(fnameU, 'nav_lat', gphiu, -1, lncdf_stat)
(  448)          CALL ncdf_write(fnameU, 'nav_lon', glamu, -1, lncdf_stat)
(  449)          CALL ncdf_write(fnameV, 'nav_lat', gphiv, -1, lncdf_stat)
(  450)          CALL ncdf_write(fnameV, 'nav_lon', glamv, -1, lncdf_stat)
(  451) !         CALL ncdf_write(fnameT, 'nav_lat', gphit, -1, lncdf_stat)
(  452) !         CALL ncdf_write(fnameT, 'nav_lon', glamt, -1, lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 9

(  453) 
(  454)       END IF  ! first timestep
(  455) 
(  456) !!DB: 2007.12.04 -- Code to compute M2 ave fields, once per day
(  457)       ddt = mod(kt-nit000,ndt) + 1  !! = timestep-# within a day
(  458) !DB 12.10 -- the below actually gives 1 more dt than I wanted i.e.
(  459) ! M2dt + 1, but now I need this extra dt
(  460)       if(ddt >= ndt-M2dt) then  !!accumulate sums
(  461)          if(ddt /= ndt) then    !! but not on the last dt of the window
(  462)             aveU = aveU + un
(  463)             aveV = aveV + vn
(  464)          endif
(  465)       endif
(  466) 
(  467) !!Output time if:
(  468) !!DB 12.10 -- the code below time-interpolates to the exact rM2dt, i.e.
(  469) !! val = (1-zxy)*sum_lo + zxy*sum_hi, where sum_hi = sum_lo + latest_field
(  470) !!     = (1-zxy)*sum_lo + zxy*(sum_lo+lf) = sum_lo + zxy*lf
(  471) !! val = val/rM2dt
(  472) !! This method avoids having to store any temporary 3D fields
(  473)       if(ddt == ndt) then
(  474)          rec_num = rec_num + 1
(  475) 
(  476)          aveU = (aveU + zxy*un)/rM2dt
(  477)          aveV = (aveV + zxy*vn)/rM2dt
(  478) !         aveU = aveU/float(M2dt)
(  479) !         aveV = aveV/float(M2dt)
(  480)          
(  481)          ! On write steps, output time counter
(  482)          CALL ncdf_write(fnameU, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)
(  483)          CALL ncdf_write(fnameV, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)
(  484) !         CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)
(  485) 
(  486)          CALL ncdf_write(fnameU, 'aveU', aveU, -rec_num, lncdf_stat)
(  487)          CALL ncdf_write(fnameV, 'aveV', aveV, -rec_num, lncdf_stat)
(  488)          
(  489)          aveU = 0.0
(  490)          aveV = 0.0
(  491) !!AD/DB 2009.09.30
(  492)          CALL ncdf_write(fnameU, 'ndastp',REAL(ndastp), rec_num, lncdf_stat)
(  493)          CALL ncdf_write(fnameU, 'model_time_step',REAL(kt), rec_num, lncdf_stat)
(  494)          CALL ncdf_write(fnameU, 'model_time',model_time, rec_num, lncdf_stat)
(  495)          CALL ncdf_write(fnameV, 'ndastp',REAL(ndastp), rec_num, lncdf_stat)
(  496)          CALL ncdf_write(fnameV, 'model_time_step',REAL(kt), rec_num, lncdf_stat)
(  497)          CALL ncdf_write(fnameV, 'model_time',model_time, rec_num, lncdf_stat)
(  498) 
(  499)       endif
(  500)       
(  501) 
(  502) 
(  503)    END SUBROUTINE output_special
(  504) 
(  505) !!DB
(  506) !!write time-averaged variables if the ioutput_ave flag > 0
(  507) !! averaging period = ioutput_ave; output is every mod(kt-nit000,ioutput_ave)
(  508) !!NB: More variables are created in the output file than are currently on the 
(  509) !!average list. 
(  510) !!2008.06.13 -- add w (but do not change name)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 10

(  511) 
(  512)    SUBROUTINE output_aveTSUV( kt, kindic )
(  513)       !!---------------------------------------------------------------------
(  514)       !!                  ***  ROUTINE dia_wri_ave  ***
(  515)       !!                   
(  516)       !! ** Purpose :   Average and output special opa variable
(  517)       !!      NETCDF format is used by default 
(  518)       !!
(  519)       !! DB: 2007.12.11 
(  520)       !!----------------------------------------------------------------------
(  521) 
(  522)       !! * Modules used
(  523) !      USE ioipsl
(  524) !     USE oce
(  525) !     USE taumod, ONLY:  taux, tauy
(  526) 
(  527)       !! * Arguments
(  528)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  529)       INTEGER, INTENT( in ) ::   kindic  ! 
(  530) 
(  531)       !! * Local declarations
(  532)       LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
(  533)       CHARACTER (len=40) ::           &
(  534)          clhstnam, clop, clmx            ! temporary names
(  535)       INTEGER ::   inum = 11             ! temporary logical unit
(  536)       INTEGER ::   &
(  537)          iimi, iima, ipk, it,         &  ! temporary integers
(  538)          ijmi, ijma                      !    "          "
(  539) 
(  540)       INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls
(  541) 
(  542)       REAL(wp) ::   &
(  543)          zsto, zout, zmax,            &  ! temporary scalars
(  544)          zjulian, zdt                    !    "         "
(  545)       REAL(wp), DIMENSION(jpi,jpj) :: &
(  546)          zw2d                            ! temporary workspace
(  547)       CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT
(  548) 
(  549) !!DB  
(  550)       REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: aveU2, aveV2, aveT2, aveS2, aveW2
(  551)       REAL(wp),DIMENSION(:,:),ALLOCATABLE, SAVE :: ave_ssh2, ave_taux, ave_tauy
(  552)       INTEGER, SAVE :: rec_num2
(  553) 
(  554) !DB  This routine should not be called if ioutput_ave = 0, but to be safe
(  555)       if(ioutput_ave == 0) then
(  556)          return
(  557)       endif
(  558) 
(  559)       ! 0. Initialisation
(  560)       fnameT = trim(cexper)//'_ave_TSUV.nc'
(  561) 
(  562) ! On the first timestep, create output files and write depth variables (since they don't change)
(  563)       IF(kt == nit000) THEN
(  564)          rec_num2 = 0
(  565)          call ncdf_create_file_aveTSUV(fnameT, 'ave(x)', lncdf_stat)
(  566) 
(  567) ! Allocate arrays, and to be safe set them to zero
(  568)          ALLOCATE(aveU2(jpi,jpj,jpk))






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 11

(  569)          ALLOCATE(aveV2(jpi,jpj,jpk))
(  570)          ALLOCATE(aveW2(jpi,jpj,jpk))
(  571)          ALLOCATE(aveT2(jpi,jpj,jpk))
(  572)          ALLOCATE(aveS2(jpi,jpj,jpk))
(  573)          ALLOCATE(ave_ssh2(jpi,jpj))
(  574)          ALLOCATE(ave_taux(jpi,jpj))
(  575)          ALLOCATE(ave_tauy(jpi,jpj))
(  576) 
(  577)          aveU2 = 0.0
(  578)          aveV2 = 0.0
(  579)          aveW2 = 0.0
(  580)          aveT2 = 0.0
(  581)          aveS2 = 0.0
(  582)          ave_ssh2 = 0.0
(  583)          ave_taux = 0.0
(  584)          ave_tauy = 0.0
(  585) 
(  586) !!DBG
(  587)          if(lwp) then 
(  588)             write(numout2,*)'DBG: diawri(...aveTSUV()) opens: ',fnameT
(  589)          endif
(  590)       
(  591) ! On first write step only
(  592)          CALL ncdf_write(fnameT, 'depthu', gdept, lncdf_stat) 
(  593)          CALL ncdf_write(fnameT, 'depthw', gdepw, lncdf_stat) 
(  594)          CALL ncdf_write(fnameT, 'nav_lat', gphit, -1, lncdf_stat)
(  595)          CALL ncdf_write(fnameT, 'nav_lon', glamt, -1, lncdf_stat)
(  596) 
(  597)       END IF  ! first timestep
(  598) 
(  599) !!DB: 2007.12.11 -- Code to compute ave fields
(  600)       aveU2 = aveU2 + un
(  601)       aveV2 = aveV2 + vn
(  602)       aveW2 = aveW2 + wn
(  603)       aveT2 = aveT2 + tn 
(  604)       aveS2 = aveS2 + sn 
(  605)       ave_ssh2 = ave_ssh2 + sshn 
(  606)       ave_taux = ave_taux + taux
(  607)       ave_tauy = ave_tauy + tauy
(  608) 
(  609) !!Output time if:
(  610)       if(mod(kt-nit000+1,ioutput_ave) ==  0) then
(  611)          rec_num2 = rec_num2 + 1
(  612)          aveU2 = aveU2/float(ioutput_ave)
(  613)          aveV2 = aveV2/float(ioutput_ave)
(  614)          aveW2 = aveW2/float(ioutput_ave)
(  615)          aveT2 = aveT2/float(ioutput_ave)
(  616)          aveS2 = aveS2/float(ioutput_ave)
(  617)          ave_ssh2 = ave_ssh2/float(ioutput_ave)
(  618)          ave_taux = ave_taux/float(ioutput_ave)
(  619)          ave_tauy = ave_tauy/float(ioutput_ave)
(  620) 
(  621)          CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), rec_num2, lncdf_stat)
(  622)          CALL ncdf_write(fnameT, 'vozocrtx', aveU2, -rec_num2, lncdf_stat)
(  623)          CALL ncdf_write(fnameT, 'vomecrty', aveV2, -rec_num2, lncdf_stat)
(  624) !!W
(  625)          CALL ncdf_write(fnameT, 'vovecrtz', aveW2, -rec_num2, lncdf_stat)
(  626)          CALL ncdf_write(fnameT, 'votemper', aveT2, -rec_num2, lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 12

(  627)          CALL ncdf_write(fnameT, 'vosaline', aveS2, -rec_num2, lncdf_stat)
(  628)          CALL ncdf_write(fnameT, 'sossheig', ave_ssh2, -rec_num2, lncdf_stat)
(  629)          CALL ncdf_write(fnameT, 'sozotaux', ave_taux, -rec_num2, lncdf_stat)
(  630)          CALL ncdf_write(fnameT, 'sometauy', ave_tauy, -rec_num2, lncdf_stat)
(  631) !!AD/DB 2009.09.30
(  632)          CALL ncdf_write(fnameT, 'ndastp',REAL(ndastp), rec_num2, lncdf_stat)
(  633)          CALL ncdf_write(fnameT, 'model_time_step',REAL(kt), rec_num2, lncdf_stat)
(  634)          CALL ncdf_write(fnameT, 'model_time',model_time, rec_num2, lncdf_stat)
(  635) 
(  636)          aveU2 = 0.0
(  637)          aveV2 = 0.0
(  638)          aveW2 = 0.0
(  639)          aveT2 = 0.0
(  640)          aveS2 = 0.0
(  641)          ave_ssh2 = 0.0
(  642)          ave_taux = 0.0
(  643)          ave_tauy = 0.0
(  644) 
(  645)       endif
(  646)       
(  647) 
(  648) 
(  649)    END SUBROUTINE output_aveTSUV
(  650) 
(  651) 
(  652) END MODULE diawri






































PGF90 (Version     10.5)          09/23/2010  09:49:30      page 13

(    1) # 1 "diawri.F90"
(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) !!DB - 2007.12.04
(    3) !!This version has eliminated the IOIPSL code
(    4) 
(    5) 
(    6) MODULE diawri
(    7)    !!======================================================================
(    8)    !!                     ***  MODULE  diawri  ***
(    9)    !! Ocean diagnostics :  write ocean output files
(   10)    !!=====================================================================
(   11) 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers 
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE zdf_oce         ! ocean vertical physics
(   17)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   18)    USE ldfdyn_oce      ! ocean dynamics: lateral physics
(   19)    USE sol_oce         ! solver variables
(   20)    USE ice_oce         ! ice variables
(   21)    USE phycst          ! physical constants
(   22)    USE ocfzpt          ! ocean freezing point
(   23)    USE ocesbc          ! surface thermohaline fluxes
(   24)    USE taumod          ! surface stress
(   25)    USE flxrnf          ! ocean runoffs
(   26)    USE zdfmxl          ! mixed layer
(   27)    USE daymod          ! calendar
(   28)    USE dianam          ! build name of file (routine)
(   29)    USE zdfddm          ! vertical  physics: double diffusion
(   30)    USE diahth          ! thermocline diagnostics
(   31)    USE diaspr          ! surface pressure diagnostics (rigid lid case)
(   32)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   33)    USE in_out_manager  ! I/O manager
(   34)    USE flx_oce         ! sea-ice/ocean forcings variables
(   35) 
(   36)    USE lib_ncdf        ! netCDF I/O library
(   37) 
(   38)    IMPLICIT NONE
(   39)    PRIVATE
(   40) 
(   41)    !! * Accessibility
(   42)    PUBLIC dia_wri                 ! routines called by step.F90
(   43)    PUBLIC dia_wri_state
(   44) !!DB: tried to move to lib_ncdf but did not work (???)
(   45)    PUBLIC output_special             ! special routines called by step.F90 if ave flag is on
(   46)    PUBLIC output_aveTSUV             ! special routines called by step.F90 if ave flag is on
(   47) 
(   48)    !! * Module variables
(   49) !   INTEGER ::   &
(   50) !      nid_T, nz_T, nh_T, ndim_T, ndim_hT,      &   ! grid_T file
(   51) !      nid_U, nz_U, nh_U, ndim_U, ndim_hU,      &   ! grid_U file
(   52) !      nid_V, nz_V, nh_V, ndim_V, ndim_hV,      &   ! grid_V file
(   53) !      nid_W, nz_W, nh_W,                       &   ! grid_W file
(   54) !      ndex(1)                                      ! ???
(   55) !   INTEGER, DIMENSION(jpi*jpj) ::   &
(   56) !      ndex_hT, ndex_hU, ndex_hV
(   57) !   INTEGER, DIMENSION(jpi*jpj*jpk) ::   &






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 14

(   58) !      ndex_T, ndex_U, ndex_V
(   59) 
(   60)    !! * Substitutions
(   61) # 1 "./zdfddm_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** zdfddm_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfddm_substitute.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 16
(   16) !   Defautl option :                     avs = avt
(   17) # 62 "diawri.F90"
(   62) # 62
(   62)    !!----------------------------------------------------------------------
(   63)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   64)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diawri.F90,v 1.10 2006/03/09 17:21:54 opalod Exp $ 
(   65)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   66)    !!----------------------------------------------------------------------
(   67) 
(   68) CONTAINS
(   69) 
(   70)    !!----------------------------------------------------------------------
(   71)    !!   Default option                                   NetCDF output file
(   72)    !!----------------------------------------------------------------------
(   73)    !!   dia_wri       : create the standart NetCDF output files
(   74)    !!   dia_wri_state : create an output NetCDF file for a single
(   75)    !!                   instantaeous ocean state and forcing fields
(   76)    !!----------------------------------------------------------------------
(   77) 
(   78)    SUBROUTINE dia_wri( kt, kindic )
(   79)       !!---------------------------------------------------------------------
(   80)       !!                  ***  ROUTINE dia_wri  ***
(   81)       !!                   
(   82)       !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
(   83)       !!      NETCDF format is used by default 
(   84)       !!
(   85)       !! ** Method  :   At the beginning of the first time step (nit000), 
(   86)       !!      define all the NETCDF files and fields
(   87)       !!      At each time step call histdef to compute the mean if ncessary
(   88)       !!      Each nwrite time step, output the instantaneous or mean fields
(   89)       !!      IF kindic <0, output of fields before the model interruption.
(   90)       !!      IF kindic =0, time step loop
(   91)       !!      IF kindic >0, output of fields before the time step loop
(   92)       !!
(   93)       !! History :
(   94)       !!        !  91-03  (M.-A. Foujols)  Original code
(   95)       !!        !  91-11  (G. Madec)
(   96)       !!        !  92-06  (M. Imbard)  correction restart file
(   97)       !!        !  92-07  (M. Imbard)  split into diawri and rstwri
(   98)       !!        !  93-03  (M. Imbard)  suppress writibm
(   99)       !!        !  98-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
(  100)       !!        !  99-02  (E. Guilyardi)  name of netCDF files + variables






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 15

(  101)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  102)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  103)       !!----------------------------------------------------------------------
(  104)       !! * Modules used
(  105) !      USE ioipsl
(  106) 
(  107)       !! * Arguments
(  108)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  109)       INTEGER, INTENT( in ) ::   kindic  ! 
(  110) 
(  111)       !! * Local declarations
(  112)       LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
(  113)       CHARACTER (len=40) ::           &
(  114)          clhstnam, clop, clmx            ! temporary names
(  115)       INTEGER ::   inum = 11             ! temporary logical unit
(  116)       INTEGER ::   &
(  117)          iimi, iima, ipk, it,         &  ! temporary integers
(  118)          ijmi, ijma                      !    "          "
(  119) 
(  120)       INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls
(  121)       LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.    ! Use IOIPSL subroutines for output
(  122)                                                  ! Note that all W-grid calls are commented out!
(  123)       LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE.   ! Use lib_ncdf subroutines for output
(  124) 
(  125)       REAL(wp) :: test_array
(  126)       INTEGER :: a, b, c
(  127) 
(  128)       REAL(wp) ::   &
(  129)          zsto, zout, zmax,            &  ! temporary scalars
(  130)          zjulian, zdt                    !    "         "
(  131)       REAL(wp), DIMENSION(jpi,jpj) :: &
(  132)          zw2d                            ! temporary workspace
(  133)       CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT
(  134) !DB
(  135)       INTEGER :: tindex
(  136) 
(  137) 
(  138)       !!----------------------------------------------------------------------
(  139)       
(  140)       ! 0. Initialisation
(  141)       ! -----------------
(  142)       fnameU = trim(cexper)//'_grid_U.nc'
(  143)       fnameV = trim(cexper)//'_grid_V.nc'
(  144)       fnameT = trim(cexper)//'_grid_T.nc'
(  145) 
(  146) !!DB deleted ...
(  147) !      IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL for output?
(  148) 
(  149)    ! LIB_NCDF CALLS
(  150) 
(  151) !DB default to this code
(  152) !   IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Are we using lib_ncdf output routines?
(  153) 
(  154)       ! On the first timestep, create output files and write depth variables (since they don't change)
(  155)       IF(kt == nit000) THEN
(  156) !!DBG
(  157)          write(numout2,*)'DBG: diawri: cexper = ', cexper
(  158)          write(numout2,*)'DBG: diawri: ',fnameU, fnameV, fnameT






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 16

(  159) 
(  160)          CALL ncdf_create_file_u(fnameU, 'inst(x)', lncdf_stat)
(  161)          CALL ncdf_create_file_v(fnameV, 'inst(x)', lncdf_stat)
(  162)          CALL ncdf_create_file_t(fnameT, 'inst(x)', lncdf_stat)
(  163)          CALL ncdf_write(fnameU, 'depthu', gdept, lncdf_stat) 
(  164)          CALL ncdf_write(fnameV, 'depthv', gdept, lncdf_stat) 
(  165)          CALL ncdf_write(fnameT, 'deptht', gdept, lncdf_stat) 
(  166)       END IF
(  167)       
(  168)       ! On first write step only, write nav_lat and nav_lon grids (routine will abort if not on a write step)
(  169)       IF(kt == (nit000+ nwrite)) THEN
(  170)          CALL ncdf_write(fnameU, 'nav_lat', gphiu, kt, lncdf_stat)
(  171)          CALL ncdf_write(fnameU, 'nav_lon', glamu, kt, lncdf_stat)
(  172)          CALL ncdf_write(fnameV, 'nav_lat', gphiv, kt, lncdf_stat)
(  173)          CALL ncdf_write(fnameV, 'nav_lon', glamv, kt, lncdf_stat)
(  174)          CALL ncdf_write(fnameT, 'nav_lat', gphit, kt, lncdf_stat)
(  175)          CALL ncdf_write(fnameT, 'nav_lon', glamt, kt, lncdf_stat)
(  176)       END IF
(  177)       
(  178)       ! On write steps, output time counter
(  179)       IF(MOD(kt, nwrite) == 0) THEN
(  180)          tindex = kt/nwrite - nit000/nwrite
(  181) !         CALL ncdf_write(fnameU, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
(  182) !         CALL ncdf_write(fnameV, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
(  183) !         CALL ncdf_write(fnameT, 'time_counter', REAL(kt * 480), NINT(REAL((kt - (nit000-1)) / nwrite)), lncdf_stat)
(  184)          CALL ncdf_write(fnameU, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
(  185)          CALL ncdf_write(fnameV, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
(  186)          CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), tindex, lncdf_stat)
(  187)       END IF
(  188)       
(  189)       ! Write other fields (calls will abort if this isn't a write step)
(  190)       CALL ncdf_write(fnameU, 'vozocrtx', un, kt, lncdf_stat)
(  191)       CALL ncdf_write(fnameU, 'sozotaux', taux, kt, lncdf_stat)
(  192)       CALL ncdf_write(fnameV, 'vomecrty', vn, kt, lncdf_stat)
(  193)       CALL ncdf_write(fnameV, 'sometauy', tauy, kt, lncdf_stat)
(  194)       CALL ncdf_write(fnameT, 'votemper', tn, kt, lncdf_stat)
(  195)       CALL ncdf_write(fnameT, 'vosaline', sn, kt, lncdf_stat)
(  196)       CALL ncdf_write(fnameT, 'sosstsst', tn(:,:,1), kt, lncdf_stat)
(  197)       CALL ncdf_write(fnameT, 'sosaline', sn(:,:,1), kt, lncdf_stat)
(  198)       CALL ncdf_write(fnameT, 'sossheig', sshn, kt, lncdf_stat)
(  199)       CALL ncdf_write(fnameT, 'sowaflup', emp, kt, lncdf_stat)
(  200)       CALL ncdf_write(fnameT, 'sorunoff', runoff, kt, lncdf_stat)
(  201)       CALL ncdf_write(fnameT, 'sowaflcd', emps, kt, lncdf_stat)
(  202)       zw2d(:,:) = emps(:,:) * sn(:,:,1) * tmask(:,:,1)
(  203)       CALL ncdf_write(fnameT, 'sosalflx', zw2d, kt, lncdf_stat)
(  204)       CALL ncdf_write(fnameT, 'sohefldo', qt, kt, lncdf_stat)
(  205)       CALL ncdf_write(fnameT, 'soshfldo', qsr, kt, lncdf_stat)
(  206)       CALL ncdf_write(fnameT, 'somxl010', hmlp, kt, lncdf_stat)
(  207)       CALL ncdf_write(fnameT, 'somixhgt', hmld, kt, lncdf_stat)
(  208)       CALL ncdf_write(fnameT, 'soicecov', freeze, kt, lncdf_stat)
(  209)       CALL ncdf_write(fnameT, 'sohefldp', qrp, kt, lncdf_stat)
(  210)       CALL ncdf_write(fnameT, 'sowafldp', erp, kt, lncdf_stat)
(  211)       zw2d(:,:) = erp(:,:) * sn(:,:,1) * tmask(:,:,1)
(  212)       CALL ncdf_write(fnameT, 'sosafldp', zw2d, kt, lncdf_stat)
(  213)       zw2d(:,:) = FLOAT( nmln(:,:) ) * tmask(:,:,1)
(  214)       CALL ncdf_write(fnameT, 'sobowlin', zw2d, kt, lncdf_stat)
(  215)       
(  216)    ! END OF LIB_NCDF CALLS






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 17

(  217) 
(  218) !DB code deleted
(  219) !   IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL output subroutines?
(  220) 
(  221)    END SUBROUTINE dia_wri
(  222) 
(  223) 
(  224)    SUBROUTINE dia_wri_state( cdfile_name )
(  225)       !!---------------------------------------------------------------------
(  226)       !!                 ***  ROUTINE dia_wri_state  ***
(  227)       !!        
(  228)       !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
(  229)       !!      the instantaneous ocean state and forcing fields.
(  230)       !!        Used to find errors in the initial state or save the last
(  231)       !!      ocean state in case of abnormal end of a simulation
(  232)       !!
(  233)       !! ** Method  :   NetCDF files using ioipsl
(  234)       !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
(  235)       !!      File 'output.abort.nc' is created in case of abnormal job end
(  236)       !!
(  237)       !! History :
(  238)       !!   8.2  !  00-06  (M. Imbard)  Original code (diabort.F)
(  239)       !!   8.5  !  02-06  (A.Bozec, E. Durand)  Original code (diainit.F)
(  240)       !!   9.0  !  02-12  (G. Madec)  merge of diabort and diainit, F90
(  241)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  242)       !!----------------------------------------------------------------------
(  243)       !! * Modules used
(  244)       USE ioipsl
(  245) 
(  246)       !! * Arguments
(  247)       CHARACTER (len=* ), INTENT( in ) ::   &
(  248)          cdfile_name      ! name of the file created
(  249) 
(  250)       !! * Local declarations
(  251)       CHARACTER (len=40) :: clop
(  252)       INTEGER  ::   &
(  253)          id_i , nz_i, nh_i       
(  254)       INTEGER, DIMENSION(1) ::   &
(  255)          idex             ! temprary workspace
(  256)       REAL(wp) ::   &
(  257)          zsto, zout, zmax,   &
(  258)          zjulian, zdt
(  259)       !!----------------------------------------------------------------------
(  260) 
(  261)       IF(lwp) WRITE(numout,*)
(  262)       IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
(  263)       IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
(  264)       IF(lwp) WRITE(numout,*) '                and named :', cdfile_name, '.nc'
(  265)       
(  266)       ! 0. Initialisation
(  267)       ! -----------------
(  268)       
(  269)       ! Define frequency of output and means
(  270)       zdt  = rdt
(  271)       zsto = rdt
(  272)       clop = "inst(x)"           ! no use of the mask value (require less cpu time)
(  273)       zout = rdt
(  274)       zmax = ( nitend - nit000 + 1 ) * zdt






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 18

(  275) 
(  276)       ! 1. Define NETCDF files and fields at beginning of first time step
(  277)       ! -----------------------------------------------------------------
(  278) 
(  279)       ! Compute julian date from starting date of the run
(  280)       CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )         ! time axis 
(  281)       CALL histbeg( cdfile_name, jpi, glamt, jpj, gphit,   &
(  282)           1, jpi, 1, jpj, 0, zjulian, zdt, nh_i, id_i, domain_id=nidom )          ! Horizontal grid : glamt and gphit
(  283)       CALL histvert( id_i, "deptht", "Vertical T levels",   &    ! Vertical grid : gdept
(  284)           "m", jpk, gdept, nz_i)
(  285) 
(  286)       ! Declare all the output fields as NetCDF variables
(  287) 
(  288)       CALL histdef( id_i, "vosaline", "Salinity"              , "PSU"    ,   &   ! salinity
(  289)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
(  290)       CALL histdef( id_i, "votemper", "Temperature"           , "C"      ,   &   ! temperature
(  291)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
(  292) # 296
(  296)       CALL histdef( id_i, "sossheig", "Sea Surface Height"    , "m"      ,   &  ! ssh
(  297)          &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
(  298) # 299
(  299)       CALL histdef( id_i, "vozocrtx", "Zonal Current"         , "m/s"    ,   &   ! zonal current
(  300)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
(  301)       CALL histdef( id_i, "vomecrty", "Meridional Current"    , "m/s"    ,   &   ! meridonal current
(  302)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
(  303)       CALL histdef( id_i, "vovecrtz", "Vertical Velocity"     , "m/s"    ,   &   ! vertical current
(  304)          &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
(  305)       CALL histdef( id_i, "sowaflup", "Net Upward Water Flux" , "Kg/m2/S",   &   ! net freshwater 
(  306)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  307)       CALL histdef( id_i, "sohefldo", "Net Downward Heat Flux", "W/m2"   ,   &   ! net heat flux
(  308)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  309)       CALL histdef( id_i, "soshfldo", "Shortwave Radiation"   , "W/m2"   ,   &   ! solar flux
(  310)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  311)       CALL histdef( id_i, "soicecov", "Ice fraction"          , "[0,1]"  ,   &   ! freeze
(  312)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  313)       CALL histdef( id_i, "sozotaux", "Zonal Wind Stress"     , "N/m2"   ,   &   ! i-wind stress
(  314)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  315)       CALL histdef( id_i, "sometauy", "Meridional Wind Stress", "N/m2"   ,   &   ! j-wind stress
(  316)          &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  317) 
(  318)       CALL histend( id_i )
(  319) 
(  320)       ! 2. Start writing data
(  321)       ! ---------------------
(  322)       ! idex(1) est utilise ssi l'avant dernier argument est diffferent de 
(  323)       ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
(  324)       ! donne le nombre d'elements, et idex la liste des indices a sortir
(  325)       idex(1) = 1   ! init to avoid compil warning
(  326)       
(  327)       ! Write all fields on T grid
(  328)       CALL histwrite( id_i, "votemper", 1, tn    , jpi*jpj*jpk, idex )    ! now temperature
(  329)       CALL histwrite( id_i, "vosaline", 1, sn    , jpi*jpj*jpk, idex )    ! now salinity
(  330) # 333
(  333)       CALL histwrite( id_i, "sossheig", 1, sshn  , jpi*jpj    , idex )    ! sea surface height
(  334) # 335
(  335)       CALL histwrite( id_i, "vozocrtx", 1, un    , jpi*jpj*jpk, idex )    ! now i-velocity
(  336)       CALL histwrite( id_i, "vomecrty", 1, vn    , jpi*jpj*jpk, idex )    ! now j-velocity
(  337)       CALL histwrite( id_i, "vovecrtz", 1, wn    , jpi*jpj*jpk, idex )    ! now k-velocity






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 19

(  338)       CALL histwrite( id_i, "sowaflup", 1, emp   , jpi*jpj    , idex )    ! freshwater budget
(  339)       CALL histwrite( id_i, "sohefldo", 1, qt    , jpi*jpj    , idex )    ! total heat flux
(  340)       CALL histwrite( id_i, "soshfldo", 1, qsr   , jpi*jpj    , idex )    ! total heat flux
(  341)       CALL histwrite( id_i, "soicecov", 1, freeze, jpi*jpj    , idex )    ! ice cover
(  342)       CALL histwrite( id_i, "sozotaux", 1, taux  , jpi*jpj    , idex )    ! i-wind stress
(  343)       CALL histwrite( id_i, "sometauy", 1, tauy  , jpi*jpj    , idex )    ! j-wind stress
(  344) 
(  345)       ! 3. Close the file
(  346)       ! -----------------
(  347)       CALL histclo( id_i )
(  348) 
(  349)    END SUBROUTINE dia_wri_state
(  350) 
(  351)    !!======================================================================
(  352) 
(  353) !!DB
(  354)    SUBROUTINE output_special( kt, kindic )
(  355)       !!---------------------------------------------------------------------
(  356)       !!                  ***  ROUTINE dia_wri_ave  ***
(  357)       !!                   
(  358)       !! ** Purpose :   Average and output special opa variable
(  359)       !!      NETCDF format is used by default 
(  360)       !!
(  361)       !! DB: 2007.12.04 ... 12.10 ... tested OK
(  362)       !!----------------------------------------------------------------------
(  363) !!DB -- notes
(  364) !1 - due to the way that lib_ncdf routine handles 'when to write' this routine will
(  365) ! have to fool the call to ncdf_write...() to write at a non-nwrite timestep:
(  366) ! fixed by inputting a -ve rec_num value
(  367) !2 - ALLOCATE arrays when routine is called ... DONE ...
(  368) !3 - Uses an improved M2-averaging routine (see dbrick:..MATLAB/test_M2ave.m)
(  369) ! that interpolates to the correct (float)time (rM2dt). Routine also gives correct
(  370) ! answer if rM2dt = integer
(  371) 
(  372)       !! * Modules used
(  373) !      USE ioipsl
(  374) !     USE oce
(  375) !     USE taumod
(  376) 
(  377)       !! * Arguments
(  378)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  379)       INTEGER, INTENT( in ) ::   kindic  ! 
(  380) 
(  381)       !! * Local declarations
(  382)       LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
(  383)       CHARACTER (len=40) ::           &
(  384)          clhstnam, clop, clmx            ! temporary names
(  385)       INTEGER ::   inum = 11             ! temporary logical unit
(  386)       INTEGER ::   &
(  387)          iimi, iima, ipk, it,         &  ! temporary integers
(  388)          ijmi, ijma                      !    "          "
(  389) 
(  390)       INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls
(  391) 
(  392)       REAL(wp) ::   &
(  393)          zsto, zout, zmax,            &  ! temporary scalars
(  394)          zjulian, zdt                    !    "         "
(  395)       REAL(wp), DIMENSION(jpi,jpj) :: &






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 20

(  396)          zw2d                            ! temporary workspace
(  397)       CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT
(  398) 
(  399) !!DB
(  400)       INTEGER :: kkk, M2dt, ndt, ddt, M2dt_lo, M2dt_hi
(  401)       REAL(wp) :: M2_hr, rM2dt, zxy
(  402)       REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: aveU, aveV, aveT
(  403)       REAL(wp),DIMENSION(:,:),ALLOCATABLE, SAVE :: ave_ssh
(  404)       INTEGER, SAVE :: rec_num
(  405) 
(  406)       ! 0. Initialisation
(  407) !!DB
(  408)       M2_hr = 12.42
(  409)       ndt = 86400./rdt
(  410)       rM2dt = M2_hr*3600./rdt
(  411)       M2dt_lo = floor(M2_hr*3600./rdt)
(  412) !      M2dt_hi = ceil(M2_hr*3600./rdt) !no ceil() but variable is not needed
(  413)       M2dt = M2dt_lo;
(  414)       zxy = rM2dt-float(M2dt)
(  415) 
(  416)       fnameU = trim(cexper)//'_M2ave_U.nc'
(  417)       fnameV = trim(cexper)//'_M2ave_V.nc'
(  418)       fnameT = trim(cexper)//'_M2ave_T.nc'
(  419) !      fnameT = trim(cexper)//'_ave_TSUV.nc'
(  420) 
(  421) ! On the first timestep, create output files and write depth variables (since they don't change)
(  422)       IF(kt == nit000) THEN
(  423)          rec_num = 0
(  424) 
(  425)          call ncdf_create_file_ave(fnameU, 'aveU', 'ave(x)', lncdf_stat)
(  426)          call ncdf_create_file_ave(fnameV, 'aveV', 'ave(x)', lncdf_stat)
(  427) !         call ncdf_create_file_ave(fnameT, 'aveT', 'ave(x)', lncdf_stat)
(  428) 
(  429)          ALLOCATE(aveU(jpi,jpj,jpk))
(  430)          ALLOCATE(aveV(jpi,jpj,jpk))
(  431) !         ALLOCATE(aveT(jpi,jpj,jpk))
(  432)          aveU = 0.0
(  433)          aveV = 0.0
(  434) 
(  435) !!DBG
(  436)          if(lwp) then
(  437)             write(numout2,*)'DBG: diawri: cexper = ', cexper
(  438)             write(numout2,*)'DBG: diawri: ',fnameU, fnameV
(  439)          endif
(  440) 
(  441) 
(  442)          CALL ncdf_write(fnameU, 'depthu', gdept, lncdf_stat) 
(  443)          CALL ncdf_write(fnameV, 'depthu', gdept, lncdf_stat) 
(  444) !         CALL ncdf_write(fnameT, 'deptht', gdept, lncdf_stat) 
(  445)       
(  446) ! On first write step only, write nav_lat and nav_lon grids (routine will abort if not on a write step)
(  447)          CALL ncdf_write(fnameU, 'nav_lat', gphiu, -1, lncdf_stat)
(  448)          CALL ncdf_write(fnameU, 'nav_lon', glamu, -1, lncdf_stat)
(  449)          CALL ncdf_write(fnameV, 'nav_lat', gphiv, -1, lncdf_stat)
(  450)          CALL ncdf_write(fnameV, 'nav_lon', glamv, -1, lncdf_stat)
(  451) !         CALL ncdf_write(fnameT, 'nav_lat', gphit, -1, lncdf_stat)
(  452) !         CALL ncdf_write(fnameT, 'nav_lon', glamt, -1, lncdf_stat)
(  453) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 21

(  454)       END IF  ! first timestep
(  455) 
(  456) !!DB: 2007.12.04 -- Code to compute M2 ave fields, once per day
(  457)       ddt = mod(kt-nit000,ndt) + 1  !! = timestep-# within a day
(  458) !DB 12.10 -- the below actually gives 1 more dt than I wanted i.e.
(  459) ! M2dt + 1, but now I need this extra dt
(  460)       if(ddt >= ndt-M2dt) then  !!accumulate sums
(  461)          if(ddt /= ndt) then    !! but not on the last dt of the window
(  462)             aveU = aveU + un
(  463)             aveV = aveV + vn
(  464)          endif
(  465)       endif
(  466) 
(  467) !!Output time if:
(  468) !!DB 12.10 -- the code below time-interpolates to the exact rM2dt, i.e.
(  469) !! val = (1-zxy)*sum_lo + zxy*sum_hi, where sum_hi = sum_lo + latest_field
(  470) !!     = (1-zxy)*sum_lo + zxy*(sum_lo+lf) = sum_lo + zxy*lf
(  471) !! val = val/rM2dt
(  472) !! This method avoids having to store any temporary 3D fields
(  473)       if(ddt == ndt) then
(  474)          rec_num = rec_num + 1
(  475) 
(  476)          aveU = (aveU + zxy*un)/rM2dt
(  477)          aveV = (aveV + zxy*vn)/rM2dt
(  478) !         aveU = aveU/float(M2dt)
(  479) !         aveV = aveV/float(M2dt)
(  480)          
(  481)          ! On write steps, output time counter
(  482)          CALL ncdf_write(fnameU, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)
(  483)          CALL ncdf_write(fnameV, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)
(  484) !         CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), rec_num, lncdf_stat)
(  485) 
(  486)          CALL ncdf_write(fnameU, 'aveU', aveU, -rec_num, lncdf_stat)
(  487)          CALL ncdf_write(fnameV, 'aveV', aveV, -rec_num, lncdf_stat)
(  488)          
(  489)          aveU = 0.0
(  490)          aveV = 0.0
(  491) !!AD/DB 2009.09.30
(  492)          CALL ncdf_write(fnameU, 'ndastp',REAL(ndastp), rec_num, lncdf_stat)
(  493)          CALL ncdf_write(fnameU, 'model_time_step',REAL(kt), rec_num, lncdf_stat)
(  494)          CALL ncdf_write(fnameU, 'model_time',model_time, rec_num, lncdf_stat)
(  495)          CALL ncdf_write(fnameV, 'ndastp',REAL(ndastp), rec_num, lncdf_stat)
(  496)          CALL ncdf_write(fnameV, 'model_time_step',REAL(kt), rec_num, lncdf_stat)
(  497)          CALL ncdf_write(fnameV, 'model_time',model_time, rec_num, lncdf_stat)
(  498) 
(  499)       endif
(  500)       
(  501) 
(  502) 
(  503)    END SUBROUTINE output_special
(  504) 
(  505) !!DB
(  506) !!write time-averaged variables if the ioutput_ave flag > 0
(  507) !! averaging period = ioutput_ave; output is every mod(kt-nit000,ioutput_ave)
(  508) !!NB: More variables are created in the output file than are currently on the 
(  509) !!average list. 
(  510) !!2008.06.13 -- add w (but do not change name)
(  511) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 22

(  512)    SUBROUTINE output_aveTSUV( kt, kindic )
(  513)       !!---------------------------------------------------------------------
(  514)       !!                  ***  ROUTINE dia_wri_ave  ***
(  515)       !!                   
(  516)       !! ** Purpose :   Average and output special opa variable
(  517)       !!      NETCDF format is used by default 
(  518)       !!
(  519)       !! DB: 2007.12.11 
(  520)       !!----------------------------------------------------------------------
(  521) 
(  522)       !! * Modules used
(  523) !      USE ioipsl
(  524) !     USE oce
(  525) !     USE taumod, ONLY:  taux, tauy
(  526) 
(  527)       !! * Arguments
(  528)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  529)       INTEGER, INTENT( in ) ::   kindic  ! 
(  530) 
(  531)       !! * Local declarations
(  532)       LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
(  533)       CHARACTER (len=40) ::           &
(  534)          clhstnam, clop, clmx            ! temporary names
(  535)       INTEGER ::   inum = 11             ! temporary logical unit
(  536)       INTEGER ::   &
(  537)          iimi, iima, ipk, it,         &  ! temporary integers
(  538)          ijmi, ijma                      !    "          "
(  539) 
(  540)       INTEGER :: lncdf_stat               ! Return status for lib_ncdf calls
(  541) 
(  542)       REAL(wp) ::   &
(  543)          zsto, zout, zmax,            &  ! temporary scalars
(  544)          zjulian, zdt                    !    "         "
(  545)       REAL(wp), DIMENSION(jpi,jpj) :: &
(  546)          zw2d                            ! temporary workspace
(  547)       CHARACTER (len=80) :: clname, fnameU, fnameV, fnameT
(  548) 
(  549) !!DB  
(  550)       REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: aveU2, aveV2, aveT2, aveS2, aveW2
(  551)       REAL(wp),DIMENSION(:,:),ALLOCATABLE, SAVE :: ave_ssh2, ave_taux, ave_tauy
(  552)       INTEGER, SAVE :: rec_num2
(  553) 
(  554) !DB  This routine should not be called if ioutput_ave = 0, but to be safe
(  555)       if(ioutput_ave == 0) then
(  556)          return
(  557)       endif
(  558) 
(  559)       ! 0. Initialisation
(  560)       fnameT = trim(cexper)//'_ave_TSUV.nc'
(  561) 
(  562) ! On the first timestep, create output files and write depth variables (since they don't change)
(  563)       IF(kt == nit000) THEN
(  564)          rec_num2 = 0
(  565)          call ncdf_create_file_aveTSUV(fnameT, 'ave(x)', lncdf_stat)
(  566) 
(  567) ! Allocate arrays, and to be safe set them to zero
(  568)          ALLOCATE(aveU2(jpi,jpj,jpk))
(  569)          ALLOCATE(aveV2(jpi,jpj,jpk))






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 23

(  570)          ALLOCATE(aveW2(jpi,jpj,jpk))
(  571)          ALLOCATE(aveT2(jpi,jpj,jpk))
(  572)          ALLOCATE(aveS2(jpi,jpj,jpk))
(  573)          ALLOCATE(ave_ssh2(jpi,jpj))
(  574)          ALLOCATE(ave_taux(jpi,jpj))
(  575)          ALLOCATE(ave_tauy(jpi,jpj))
(  576) 
(  577)          aveU2 = 0.0
(  578)          aveV2 = 0.0
(  579)          aveW2 = 0.0
(  580)          aveT2 = 0.0
(  581)          aveS2 = 0.0
(  582)          ave_ssh2 = 0.0
(  583)          ave_taux = 0.0
(  584)          ave_tauy = 0.0
(  585) 
(  586) !!DBG
(  587)          if(lwp) then 
(  588)             write(numout2,*)'DBG: diawri(...aveTSUV()) opens: ',fnameT
(  589)          endif
(  590)       
(  591) ! On first write step only
(  592)          CALL ncdf_write(fnameT, 'depthu', gdept, lncdf_stat) 
(  593)          CALL ncdf_write(fnameT, 'depthw', gdepw, lncdf_stat) 
(  594)          CALL ncdf_write(fnameT, 'nav_lat', gphit, -1, lncdf_stat)
(  595)          CALL ncdf_write(fnameT, 'nav_lon', glamt, -1, lncdf_stat)
(  596) 
(  597)       END IF  ! first timestep
(  598) 
(  599) !!DB: 2007.12.11 -- Code to compute ave fields
(  600)       aveU2 = aveU2 + un
(  601)       aveV2 = aveV2 + vn
(  602)       aveW2 = aveW2 + wn
(  603)       aveT2 = aveT2 + tn 
(  604)       aveS2 = aveS2 + sn 
(  605)       ave_ssh2 = ave_ssh2 + sshn 
(  606)       ave_taux = ave_taux + taux
(  607)       ave_tauy = ave_tauy + tauy
(  608) 
(  609) !!Output time if:
(  610)       if(mod(kt-nit000+1,ioutput_ave) ==  0) then
(  611)          rec_num2 = rec_num2 + 1
(  612)          aveU2 = aveU2/float(ioutput_ave)
(  613)          aveV2 = aveV2/float(ioutput_ave)
(  614)          aveW2 = aveW2/float(ioutput_ave)
(  615)          aveT2 = aveT2/float(ioutput_ave)
(  616)          aveS2 = aveS2/float(ioutput_ave)
(  617)          ave_ssh2 = ave_ssh2/float(ioutput_ave)
(  618)          ave_taux = ave_taux/float(ioutput_ave)
(  619)          ave_tauy = ave_tauy/float(ioutput_ave)
(  620) 
(  621)          CALL ncdf_write(fnameT, 'time_counter', REAL(kt * int(rdt)), rec_num2, lncdf_stat)
(  622)          CALL ncdf_write(fnameT, 'vozocrtx', aveU2, -rec_num2, lncdf_stat)
(  623)          CALL ncdf_write(fnameT, 'vomecrty', aveV2, -rec_num2, lncdf_stat)
(  624) !!W
(  625)          CALL ncdf_write(fnameT, 'vovecrtz', aveW2, -rec_num2, lncdf_stat)
(  626)          CALL ncdf_write(fnameT, 'votemper', aveT2, -rec_num2, lncdf_stat)
(  627)          CALL ncdf_write(fnameT, 'vosaline', aveS2, -rec_num2, lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 24

(  628)          CALL ncdf_write(fnameT, 'sossheig', ave_ssh2, -rec_num2, lncdf_stat)
(  629)          CALL ncdf_write(fnameT, 'sozotaux', ave_taux, -rec_num2, lncdf_stat)
(  630)          CALL ncdf_write(fnameT, 'sometauy', ave_tauy, -rec_num2, lncdf_stat)
(  631) !!AD/DB 2009.09.30
(  632)          CALL ncdf_write(fnameT, 'ndastp',REAL(ndastp), rec_num2, lncdf_stat)
(  633)          CALL ncdf_write(fnameT, 'model_time_step',REAL(kt), rec_num2, lncdf_stat)
(  634)          CALL ncdf_write(fnameT, 'model_time',model_time, rec_num2, lncdf_stat)
(  635) 
(  636)          aveU2 = 0.0
(  637)          aveV2 = 0.0
(  638)          aveW2 = 0.0
(  639)          aveT2 = 0.0
(  640)          aveS2 = 0.0
(  641)          ave_ssh2 = 0.0
(  642)          ave_taux = 0.0
(  643)          ave_tauy = 0.0
(  644) 
(  645)       endif
(  646)       
(  647) 
(  648) 
(  649)    END SUBROUTINE output_aveTSUV
(  650) 
(  651) 
(  652) END MODULE diawri
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (diawri.F90: 213)
