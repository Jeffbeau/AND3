


PGF90 (Version     10.5)          09/23/2010  09:49:31      page 1

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

Filename: diawri1d.F90

(    1) MODULE diawri1d
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  diawri1d  ***
(    4)    !! Ocean diagnostics :  write ocean output files
(    5)    !!=====================================================================
(    6) #if defined key_cfg_1d
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_cfg_1d'               1D Configuration
(    9)    !!----------------------------------------------------------------------  
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers 
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE zdf_oce         ! ocean vertical physics
(   15)    USE zdftke          ! TKE vertical mixing
(   16)    USE zdfkpp          ! KPP vertical mixing
(   17)    USE sol_oce         ! solver variables
(   18)    USE ice_oce         ! ice variables
(   19)    USE phycst          ! physical constants
(   20)    USE ocfzpt          ! ???
(   21)    USE ocesbc          ! surface thermohaline fluxes
(   22)    USE taumod          ! surface stress
(   23)    USE flxrnf          ! ???
(   24)    USE zdfmxl          ! mixed layer
(   25)    USE daymod          ! calendar
(   26)    USE dianam          ! build name of file (routine)
(   27)    USE diawri
(   28)    USE zdfddm          ! vertical  physics: double diffusion
(   29)    USE diahth          ! thermocline diagnostics
(   30)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   31)    USE in_out_manager  ! I/O manager
(   32) 
(   33)    IMPLICIT NONE
(   34)    PRIVATE
(   35) 
(   36)    !! * Accessibility
(   37)    PUBLIC dia_wri_1d                 ! routines called by step.F90
(   38)    !! * Module variables
(   39)    INTEGER ::   &
(   40)       nid_T, nz_T, nh_T, ndim_T, ndim_hT,      &   ! grid_T file
(   41)       ndex(1)                                      ! ???
(   42)    INTEGER, DIMENSION(jpi*jpj) ::   &
(   43)       ndex_hT
(   44)    INTEGER, DIMENSION(jpi*jpj*jpk) ::   &
(   45)       ndex_T
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 2

(   47)    !! * Substitutions
(   48) #  include "zdfddm_substitute.h90"
(   49)    !!----------------------------------------------------------------------
(   50)    !!   OPA 9.0 , LOCEAN-IPSL  (2005)
(   51)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/diawri1d.F90,v 1.5 2005/12/21 10:46:26 opalod Exp $ 
(   52)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   53)    !!----------------------------------------------------------------------
(   54) 
(   55) CONTAINS
(   56)    !!----------------------------------------------------------------------
(   57)    !!   Default option                                   NetCDF output file
(   58)    !!----------------------------------------------------------------------
(   59)    !!   dia_wri_1d       : create the standart NetCDF output files
(   60)    !!   dia_wri_state_1d : create an output NetCDF file for a single
(   61)    !!                      instantaeous ocean state and forcing fields
(   62)    !!----------------------------------------------------------------------
(   63) 
(   64)    SUBROUTINE dia_wri_1d( kt, kindic )
(   65)       !!---------------------------------------------------------------------
(   66)       !!                  ***  ROUTINE dia_wri_1d  ***
(   67)       !!                   
(   68)       !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
(   69)       !!      NETCDF format is used by default 
(   70)       !!
(   71)       !! ** Method  :   At the beginning of the first time step (nit000), 
(   72)       !!      define all the NETCDF files and fields
(   73)       !!      At each time step call histdef to compute the mean if ncessary
(   74)       !!      Each nwrite time step, output the instantaneous or mean fields
(   75)       !!      IF kindic <0, output of fields before the model interruption.
(   76)       !!      IF kindic =0, time step loop
(   77)       !!      IF kindic >0, output of fields before the time step loop
(   78)       !!
(   79)       !! History :
(   80)       !!        !  91-03  (M.-A. Foujols)  Original code
(   81)       !!        !  91-11  (G. Madec)
(   82)       !!        !  92-06  (M. Imbard)  correction restart file
(   83)       !!        !  92-07  (M. Imbard)  split into diawri and rstwri
(   84)       !!        !  93-03  (M. Imbard)  suppress writibm
(   85)       !!        !  98-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
(   86)       !!        !  99-02  (E. Guilyardi)  name of netCDF files + variables
(   87)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   88)       !!   9.0  !  04-10  (C. Ethe)   1D Configuration
(   89)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   90)       !!----------------------------------------------------------------------
(   91)       !! * Modules used
(   92)       USE ioipsl
(   93) 
(   94)       !! * Arguments
(   95)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   96)       INTEGER, INTENT( in ) ::   kindic  ! 
(   97) 
(   98)       !! * Local declarations
(   99)       LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
(  100)       CHARACTER (len=40) ::           &
(  101)          clhstnam, clop, clmx            ! temporary names
(  102)       INTEGER ::   inum = 11             ! temporary logical unit
(  103)       INTEGER ::   &
(  104)          ji, jj, ik                      ! dummy loop indices






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 3

(  105)       INTEGER ::   &
(  106)          iimi, iima, ipk, it,         &  ! temporary integers
(  107)          ijmi, ijma                      !    "          "
(  108)       REAL(wp) ::   &
(  109)          zsto, zout, zmax,            &  ! temporary scalars
(  110)          zjulian, zdt                    !    "         "
(  111)       REAL(wp), DIMENSION(jpi,jpj) :: &
(  112)          zw2d                            ! temporary workspace
(  113)       !!----------------------------------------------------------------------
(  114)       
(  115)       ! 0. Initialisation
(  116)       ! -----------------
(  117)       
(  118)       ! local variable for debugging
(  119)       ll_print = .FALSE.
(  120)       ll_print = ll_print .AND. lwp
(  121) 
(  122)       ! Define frequency of output and means
(  123)       zdt = rdt
(  124)       IF( nacc == 1 ) zdt = rdtmin
(  125) #if defined key_diainstant
(  126)       zsto = nwrite * zdt
(  127)       clop = "inst(x)"           ! no use of the mask value (require less cpu time)
(  128)       !!! clop="inst(only(x))"   ! put 1.e+20 on land (very expensive!!)
(  129) #else
(  130)       zsto=zdt
(  131)       clop="ave(x)"              ! no use of the mask value (require less cpu time)
(  132)       !!! clop="ave(only(x))"    ! put 1.e+20 on land (very expensive!!)
(  133) #endif
(  134)       zout = nwrite * zdt
(  135)       zmax = ( nitend - nit000 + 1 ) * zdt
(  136) 
(  137)       ! Define indices of the horizontal output zoom and vertical limit storage
(  138)       iimi = 1      ;      iima = jpi
(  139)       ijmi = 1      ;      ijma = jpj
(  140)       ipk = jpk
(  141) 
(  142)       ! define time axis
(  143)       it = kt - nit000 + 1
(  144) 
(  145) 
(  146)       ! 1. Define NETCDF files and fields at beginning of first time step
(  147)       ! -----------------------------------------------------------------
(  148) 
(  149)       IF(ll_print) WRITE(numout,*) 'dia_wri_1d kt = ', kt, ' kindic ', kindic
(  150) 
(  151)       IF( kt == nit000 ) THEN
(  152) 
(  153)          ! Define the NETCDF files (one per grid)
(  154)          
(  155)          ! Compute julian date from starting date of the run
(  156)          CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
(  157)          IF(lwp)WRITE(numout,*)
(  158)          IF(lwp)WRITE(numout,*) 'Date 0 used :', nit000, ' YEAR ', nyear,   &
(  159)             &                    ' MONTH ', nmonth, ' DAY ', nday, 'Julian day : ', zjulian
(  160)          IF(lwp)WRITE(numout,*) ' indexes of zoom = ', iimi, iima, ijmi, ijma,   &
(  161)                                  ' limit storage in depth = ', ipk
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 4

(  163)          ! WRITE root name in date.file for use by postpro
(  164)          CALL dia_nam( clhstnam, nwrite,' ' )
(  165)          CALL ctlopn( inum, 'date.file', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
(  166)          WRITE(inum,*) clhstnam
(  167)          CLOSE(inum)
(  168)          
(  169)          ! Define the T grid FILE ( nid_T )
(  170)          
(  171)          CALL dia_nam( clhstnam, nwrite, 'grid_T' )
(  172)          IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
(  173)          CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
(  174)             &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
(  175)             &          0, zjulian, zdt, nh_T, nid_T, domain_id=nidom )
(  176)          CALL histvert( nid_T, "deptht", "Vertical T levels",      &  ! Vertical grid: gdept
(  177)             &           "m", ipk, gdept, nz_T )
(  178)          !                                                            ! Index of ocean points
(  179)          CALL wheneq( jpi*jpj*ipk, tmask, 1, 1., ndex_T , ndim_T  )      ! volume
(  180)          CALL wheneq( jpi*jpj    , tmask, 1, 1., ndex_hT, ndim_hT )      ! surface
(  181) 
(  182) 
(  183)          ! Declare all the output fields as NETCDF variables
(  184) 
(  185)          !                                                                                      !!! nid_T : 3D
(  186)          CALL histdef( nid_T, "votemper", "Temperature"                        , "C"      ,   &  ! tn
(  187)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  188)          CALL histdef( nid_T, "vosaline", "Salinity"                           , "PSU"    ,   &  ! sn
(  189)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  190)          !                                                                                      !!! nid_T : 2D
(  191)          CALL histdef( nid_T, "sosstsst", "Sea Surface temperature"            , "C"      ,   &  ! sst
(  192)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  193)          CALL histdef( nid_T, "sosaline", "Sea Surface Salinity"               , "PSU"    ,   &  ! sss
(  194)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  195) 
(  196) #if ! defined key_dynspg_rl && defined key_ice_lim
(  197)          ! sowaflup = sowaflep + sorunoff + sowafldp + a term associated to
(  198)          !    internal damping to Levitus that can be diagnosed from others
(  199)          ! sowaflcd = sowaflep + sorunoff + sowafldp + iowaflup
(  200)          CALL histdef( nid_T, "iowaflup", "Ice=>ocean net freshwater"          , "kg/m2/s",   &  ! fsalt
(  201)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  202)          CALL histdef( nid_T, "sowaflep", "atmos=>ocean net freshwater"        , "kg/m2/s",   &  ! fmass
(  203)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  204) #endif
(  205)          CALL histdef( nid_T, "sowaflup", "Net Upward Water Flux"              , "Kg/m2/s",   &  ! emp
(  206)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  207)          CALL histdef( nid_T, "sorunoff", "Runoffs"                            , "Kg/m2/s",   &  ! runoffs
(  208)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  209)          CALL histdef( nid_T, "sowaflcd", "concentration/dilution water flux"  , "kg/m2/s",   &  ! emps
(  210)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  211)          CALL histdef( nid_T, "sosalflx", "Surface Salt Flux"                  , "Kg/m2/s",   &  ! emps * sn
(  212)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  213)          CALL histdef( nid_T, "sohefldo", "Net Downward Heat Flux"             , "W/m2"   ,   &  ! qt
(  214)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  215)          CALL histdef( nid_T, "soshfldo", "Shortwave Radiation"                , "W/m2"   ,   &  ! qsr
(  216)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  217)          CALL histdef( nid_T, "somxl010", "Mixed Layer Depth 0.01"             , "m"      ,   &  ! hmlp
(  218)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  219) #if defined key_zdfkpp
(  220)          CALL histdef( nid_T, "sokppekd", "Ekman depth                     "   , "m"      ,   &  ! sokppekd






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 5

(  221)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  222)          CALL histdef( nid_T, "sokppbld", "Boundary Layer Depth            "   , "m"      ,   &  ! sokppbld
(  223)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  224) #endif
(  225)          CALL histdef( nid_T, "somxlavt", "AVT : bottom of the mixed layer    ", "m"      ,   &  ! avt_mxl
(  226)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  227)          CALL histdef( nid_T, "somixhgt", "Turbocline Depth"                   , "m"      ,   &  ! hmld
(  228)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  229)          CALL histdef( nid_T, "soicecov", "Ice Cover"                          , "[0,1]"  ,   &  ! freeze
(  230)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  231) #if ! defined key_coupled 
(  232)          CALL histdef( nid_T, "sohefldp", "Surface Heat Flux: Damping"         , "W/m2"   ,   &  ! qrp
(  233)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  234)          CALL histdef( nid_T, "sowafldp", "Surface Water Flux: Damping"        , "Kg/m2/s",   &  ! erp
(  235)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  236)          CALL histdef( nid_T, "sosafldp", "Surface salt flux: damping"         , "Kg/m2/s",   &  ! erp * sn
(  237)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  238) #endif
(  239) 
(  240) #if ( defined key_coupled && ! defined key_ice_lim ) 
(  241)          CALL histdef( nid_T, "sohefldp", "Surface Heat Flux: Damping"         , "W/m2"   ,   &  ! qrp
(  242)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  243)          CALL histdef( nid_T, "sowafldp", "Surface Water Flux: Damping"        , "Kg/m2/s",   &  ! erp
(  244)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  245)          CALL histdef( nid_T, "sosafldp", "Surface salt flux: Damping"         , "Kg/m2/s",   &  ! erp * sn
(  246)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  247) #endif
(  248)          clmx ="l_max(only(x))"    ! max index on a period
(  249)          CALL histdef( nid_T, "sobowlin", "Bowl Index"                         , "W-point",   &  ! bowl INDEX 
(  250)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clmx, zsto, zout )
(  251) #if defined key_diahth
(  252)          CALL histdef( nid_T, "sothedep", "Thermocline Depth"                  , "m"      ,   & ! hth
(  253)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  254)          CALL histdef( nid_T, "so20chgt", "Depth of 20C isotherm"              , "m"      ,   & ! hd20
(  255)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  256)          CALL histdef( nid_T, "so28chgt", "Depth of 28C isotherm"              , "m"      ,   & ! hd28
(  257)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  258)          CALL histdef( nid_T, "sohtc300", "Heat content 300 m"                 , "W"      ,   & ! htc3
(  259)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  260) #endif
(  261) 
(  262) #if defined key_ice_lim && defined key_coupled
(  263)          CALL histdef( nid_T,"soicetem" , "Ice Surface Temperature"            , "K"      ,   &  ! tn_ice
(  264)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  265)          CALL histdef( nid_T,"soicealb" , "Ice Albedo"                         , "[0,1]"  ,   &  ! alb_ice
(  266)             &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  267) #endif 
(  268) 
(  269)          !                                                                                      !!! nid_U : 3D
(  270)          CALL histdef( nid_T, "vozocrtx", "Zonal Current"                      , "m/s"    ,   &  ! un
(  271)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  272) #if defined key_diaeiv
(  273)          CALL histdef( nid_T, "vozoeivu", "Zonal EIV Current"                  , "m/s"    ,   &  ! u_eiv
(  274)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  275) #endif
(  276)          !                                                                                      !!! nid_U : 2D
(  277)          CALL histdef( nid_T, "sozotaux", "Wind Stress along i-axis"           , "N/m2"   ,   &  ! taux
(  278)             &          jpi, jpj, nh_T, 1  , 1, 1  , - 99, 32, clop, zsto, zout )






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 6

(  279) 
(  280)          !                                                                                      !!! nid_V : 3D
(  281)          CALL histdef( nid_T, "vomecrty", "Meridional Current"                 , "m/s"    ,   &  ! vn
(  282)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  283) #if defined key_diaeiv
(  284)          CALL histdef( nid_T, "vomeeivv", "Meridional EIV Current"             , "m/s"    ,   &  ! v_eiv
(  285)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  286) #endif
(  287)          !                                                                                      !!! nid_V : 2D
(  288)          CALL histdef( nid_T, "sometauy", "Wind Stress along j-axis"           , "N/m2"   ,   &  ! tauy
(  289)             &          jpi, jpj, nh_T, 1  , 1, 1  , - 99, 32, clop, zsto, zout )
(  290) #if defined key_zdftke
(  291)          !                                                                                      !!! nid_W : 3D
(  292)          CALL histdef( nid_T, "votlsdis", " Dissipation Turbulent Lenght Scale", "m"      ,   &  ! e_dis
(  293)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  294)          !
(  295)          CALL histdef( nid_T, "votlsmix", " Mixing Turbulent Lenght Scale"     , "m"      ,   &  ! e_mix
(  296)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  297)          !
(  298)          CALL histdef( nid_T, "votlspdl", " Prandl Number",                      "-"       ,   &  ! e_pdl
(  299)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  300)          !
(  301)          CALL histdef( nid_T, "votlsric", " Local Richardson Number",            "-"       ,   &  ! e_ric
(  302)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  303)          !
(  304)          CALL histdef( nid_T, "votkeend", "TKE: Turbulent kinetic energy"       , "m2/s"   ,   &  ! TKE
(  305)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  306) #endif
(  307) #if defined key_zdfkpp
(  308)          !                                                                                      !!! nid_W : 3D
(  309)          CALL histdef( nid_T, "vokpprig", " Gradient Richardson Number"        ,  "-"      ,   &  ! rig
(  310)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  311)          !
(  312)          CALL histdef( nid_T, "vokpprib", " Bulk Richardson Number    "        ,  "-"      ,   &   ! rib
(  313)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  314)          !
(  315)          CALL histdef( nid_T, "vokppbsf", " Buoyancy forcing          "        , "N/m2"    ,   &  ! sokppbsf
(  316)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  317)          !
(  318)          CALL histdef( nid_T, "vokppmol", "Moning Obukhov length scale     "   , "m"       ,   &  ! sokppmol
(  319)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  320) #endif
(  321)          !
(  322)          CALL histdef( nid_T, "voeosbn2", "Brunt-Vaisala Frequency"             , "m2/s2"  ,   &  ! rn2
(  323)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  324) 
(  325)          CALL histdef( nid_T, "votkeavt", "Vertical Eddy Diffusivity"          , "m2/s"   ,   &  ! avt
(  326)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  327) 
(  328)          CALL histdef( nid_T, "votkeevd", "Enhanced Vertical Diffusivity",       "m2/s"   ,   &  ! avt_evd
(  329)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  330) 
(  331)          CALL histdef( nid_T, "votkeavm", "Vertical Eddy Viscosity",             "m2/s"   ,   &  ! avmu
(  332)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  333)          !
(  334)          CALL histdef( nid_T, "votkeevm", "Enhanced Vertical Viscosity",         "m2/s"   ,   &  ! avmu_evd
(  335)             &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 7

(  337)          IF( lk_zdfddm ) THEN
(  338)             CALL histdef( nid_T,"voddmavs","Salt Vertical Eddy Diffusivity"    , "m2/s"   ,   &  ! avs
(  339)                &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
(  340)          ENDIF
(  341) 
(  342)          CALL histend( nid_T )
(  343) 
(  344)          IF(lwp) WRITE(numout,*)
(  345)          IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization'
(  346)          IF(ll_print) CALL FLUSH(numout )
(  347) 
(  348)       ENDIF
(  349) 
(  350)       ! 2. Start writing data
(  351)       ! ---------------------
(  352) 
(  353)       ! ndex(1) est utilise ssi l'avant dernier argument est diffferent de 
(  354)       ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
(  355)       ! donne le nombre d'elements, et ndex la liste des indices a sortir
(  356) 
(  357)       IF( lwp .AND. MOD( kt, nwrite ) == 0 ) THEN 
(  358)          WRITE(numout,*) 'dia_wri : write model outputs in NetCDF files at ', kt, 'time-step'
(  359)          WRITE(numout,*) '~~~~~~ '
(  360)       ENDIF
(  361) 
(  362)       ! Write fields on T grid
(  363)       CALL histwrite( nid_T, "votemper", it, tn            , ndim_T , ndex_T  )   ! temperature
(  364)       CALL histwrite( nid_T, "vosaline", it, sn            , ndim_T , ndex_T  )   ! salinity
(  365)       CALL histwrite( nid_T, "sosstsst", it, tn(:,:,1)     , ndim_hT, ndex_hT )   ! sea surface temperature
(  366)       CALL histwrite( nid_T, "sosaline", it, sn(:,:,1)     , ndim_hT, ndex_hT )   ! sea surface salinity
(  367) #if ! defined key_dynspg_rl && defined key_ice_lim
(  368)       CALL histwrite( nid_T, "iowaflup", it, fsalt(:,:)    , ndim_hT, ndex_hT )   ! ice=>ocean water flux
(  369)       CALL histwrite( nid_T, "sowaflep", it, fmass(:,:)    , ndim_hT, ndex_hT )   ! atmos=>ocean water flux
(  370) #endif
(  371)       CALL histwrite( nid_T, "sowaflup", it, emp           , ndim_hT, ndex_hT )   ! upward water flux
(  372)       CALL histwrite( nid_T, "sorunoff", it, runoff        , ndim_hT, ndex_hT )   ! runoff
(  373)       CALL histwrite( nid_T, "sowaflcd", it, emps          , ndim_hT, ndex_hT )   ! c/d water flux
(  374)       zw2d(:,:) = emps(:,:) * sn(:,:,1) * tmask(:,:,1)
(  375)       CALL histwrite( nid_T, "sosalflx", it, zw2d          , ndim_hT, ndex_hT )   ! c/d salt flux
(  376)       CALL histwrite( nid_T, "sohefldo", it, qt            , ndim_hT, ndex_hT )   ! total heat flux
(  377)       CALL histwrite( nid_T, "soshfldo", it, qsr           , ndim_hT, ndex_hT )   ! solar heat flux
(  378)       CALL histwrite( nid_T, "somxl010", it, hmlp          , ndim_hT, ndex_hT )   ! mixed layer depth
(  379) #if defined key_zdfkpp
(  380)       CALL histwrite( nid_T, "sokppekd", it, ekdp          , ndim_hT, ndex_hT )   ! Ekman depht
(  381)       CALL histwrite( nid_T, "sokppbld", it, hkpp          , ndim_hT, ndex_hT )   ! boundary layer depth 
(  382) #endif  
(  383)       ! store the vertical eddy diffusivity coef. at the bottom of the mixed layer
(  384)       DO jj = 1, jpj
(  385)          DO ji = 1, jpi
(  386)             ik = nmln(ji,jj)
(  387)             zw2d(ji,jj) = avt(ji,jj,ik) * tmask(ji,jj,1)
(  388)          END DO
(  389)       END DO
(  390)       CALL histwrite( nid_T, "somxlavt", it, zw2d          , ndim_hT, ndex_hT )   ! Kz at bottom of mixed layer 
(  391)       CALL histwrite( nid_T, "somixhgt", it, hmld          , ndim_hT, ndex_hT )   ! turbocline depth
(  392)       CALL histwrite( nid_T, "soicecov", it, freeze        , ndim_hT, ndex_hT )   ! ice cover 
(  393) #if ! defined key_coupled
(  394)       CALL histwrite( nid_T, "sohefldp", it, qrp           , ndim_hT, ndex_hT )   ! heat flux damping






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 8

(  395)       CALL histwrite( nid_T, "sowafldp", it, erp           , ndim_hT, ndex_hT )   ! freshwater flux damping
(  396)       zw2d(:,:) = erp(:,:) * sn(:,:,1) * tmask(:,:,1)
(  397)       CALL histwrite( nid_T, "sosafldp", it, zw2d          , ndim_hT, ndex_hT )   ! salt flux damping
(  398) #endif
(  399) #if ( defined key_coupled && ! defined key_ice_lim ) 
(  400)       CALL histwrite( nid_T, "sohefldp", it, qrp           , ndim_hT, ndex_hT )   ! heat flux damping
(  401)       CALL histwrite( nid_T, "sowafldp", it, erp           , ndim_hT, ndex_hT )   ! freshwater flux damping
(  402)          zw2d(:,:) = erp(:,:) * sn(:,:,1) * tmask(:,:,1)
(  403)       CALL histwrite( nid_T, "sosafldp", it, zw2d          , ndim_hT, ndex_hT )   ! salt flux damping
(  404) #endif
(  405)          zw2d(:,:) = FLOAT( nmln(:,:) ) * tmask(:,:,1)
(  406)       CALL histwrite( nid_T, "sobowlin", it, zw2d          , ndim_hT, ndex_hT )   ! ???
(  407) 
(  408) #if defined key_diahth
(  409)       CALL histwrite( nid_T, "sothedep", it, hth           , ndim_hT, ndex_hT )   ! depth of the thermocline
(  410)       CALL histwrite( nid_T, "so20chgt", it, hd20          , ndim_hT, ndex_hT )   ! depth of the 20 isotherm
(  411)       CALL histwrite( nid_T, "so28chgt", it, hd28          , ndim_hT, ndex_hT )   ! depth of the 28 isotherm
(  412)       CALL histwrite( nid_T, "sohtc300", it, htc3          , ndim_hT, ndex_hT )   ! first 300m heaat content
(  413) #endif
(  414) #if defined key_ice_lim &&  defined key_coupled 
(  415)       CALL histwrite( nid_T, "soicetem", it, tn_ice        , ndim_hT, ndex_hT )   ! surf. ice temperature
(  416)       CALL histwrite( nid_T, "soicealb", it, alb_ice       , ndim_hT, ndex_hT )   ! ice albedo
(  417) #endif
(  418) 
(  419)       CALL histwrite( nid_T, "vozocrtx", it, un            , ndim_T , ndex_T )    ! i-current
(  420)       CALL histwrite( nid_T, "sozotaux", it, taux          , ndim_hT, ndex_hT )   ! i-wind stress
(  421)       CALL histwrite( nid_T, "vomecrty", it, vn            , ndim_T , ndex_T  )   ! j-current
(  422)       CALL histwrite( nid_T, "sometauy", it, tauy          , ndim_hT, ndex_hT )   ! j-wind stress
(  423) #if defined key_zdftke
(  424)       CALL histwrite( nid_T, "votlsdis", it, e_dis         , ndim_T , ndex_T )    ! Diss. Turb. lenght scale
(  425)       CALL histwrite( nid_T, "votlsmix", it, e_mix         , ndim_T , ndex_T )    ! Mixing Turb. lenght scale
(  426)       CALL histwrite( nid_T, "votlspdl", it, e_pdl         , ndim_T , ndex_T )    ! Prandl number
(  427)       CALL histwrite( nid_T, "votlsric", it, e_ric         , ndim_T , ndex_T )    ! local Richardson number
(  428)       CALL histwrite( nid_T, "votkeend", it, en            , ndim_T , ndex_T )    ! TKE
(  429) #endif
(  430) #if defined key_zdfkpp
(  431)       CALL histwrite( nid_T, "vokpprig", it, rig           , ndim_T , ndex_T )    ! gradient Richardson number
(  432)       CALL histwrite( nid_T, "vokpprib", it, rib           , ndim_T , ndex_T )    ! bulk Richardson number
(  433)       CALL histwrite( nid_T, "vokppbsf", it, buof          , ndim_T , ndex_T )    ! buoyancy forcing
(  434)       CALL histwrite( nid_T, "vokppmol", it, mols          , ndim_T , ndex_T )    ! Moning-Obukov length scale
(  435) #endif
(  436)       CALL histwrite( nid_T, "voeosbn2", it, rn2           , ndim_T , ndex_T )    ! Brunt-Vaisala Frequency
(  437)       CALL histwrite( nid_T, "votkeavt", it, avt           , ndim_T , ndex_T )    ! T vert. eddy diff. coef.
(  438)       CALL histwrite( nid_T, "votkeevd", it, avt_evd       , ndim_T , ndex_T )    ! T enhan. vert. eddy diff. coef.
(  439)       CALL histwrite( nid_T, "votkeavm", it, avmu          , ndim_T , ndex_T )    ! T vert. eddy visc. coef.
(  440)       CALL histwrite( nid_T, "votkeevm", it, avmu_evd      , ndim_T , ndex_T )    ! T enhan. vert. eddy visc. coef.
(  441)       IF( lk_zdfddm ) THEN
(  442)          CALL histwrite( nid_T, "voddmavs", it, fsavs(:,:,:), ndim_T, ndex_T )    ! S vert. eddy diff. coef.
(  443)       ENDIF
(  444) 
(  445)       ! 3. Synchronise and close all files
(  446)       ! ---------------------------------------
(  447)       IF( MOD( kt, nwrite ) == 0 .OR. kindic < 0 ) THEN
(  448)          CALL histsync( nid_T )
(  449)       ENDIF
(  450) 
(  451)       !  Create an output files (output.abort.nc) if S < 0 or u > 20 m/s
(  452)       IF( kindic < 0 )   CALL dia_wri_state( 'output.abort' )






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 9

(  453) 
(  454)       IF( kt == nitend .OR. kindic < 0 ) THEN
(  455)          CALL histclo( nid_T )
(  456)       ENDIF
(  457) 
(  458)    END SUBROUTINE dia_wri_1d
(  459) #else
(  460)    !!----------------------------------------------------------------------
(  461)    !!   Default key                                     NO 1D Config
(  462)    !!----------------------------------------------------------------------
(  463)    USE in_out_manager
(  464) CONTAINS
(  465)    SUBROUTINE dia_wri_1d ( kt, kindic )
(  466)       if(lwp) WRITE(numout,*) 'dia_wri_1d: You should not have seen this print! error?', kt, kindic
(  467)    END SUBROUTINE dia_wri_1d
(  468) #endif
(  469) 
(  470)    !!======================================================================
(  471) END MODULE diawri1d













































PGF90 (Version     10.5)          09/23/2010  09:49:31      page 10

(    1) # 1 "diawri1d.F90"
(    1) MODULE diawri1d
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  diawri1d  ***
(    4)    !! Ocean diagnostics :  write ocean output files
(    5)    !!=====================================================================
(    6) # 460
(  460)    !!----------------------------------------------------------------------
(  461)    !!   Default key                                     NO 1D Config
(  462)    !!----------------------------------------------------------------------
(  463)    USE in_out_manager
(  464) CONTAINS
(  465)    SUBROUTINE dia_wri_1d ( kt, kindic )
(  466)       if(lwp) WRITE(numout,*) 'dia_wri_1d: You should not have seen this print! error?', kt, kindic
(  467)    END SUBROUTINE dia_wri_1d
(  468) 
(  469) # 470
(  470)    !!======================================================================
(  471) END MODULE diawri1d
