


PGF90 (Version     10.5)          09/23/2010  09:49:05      page 1

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

Filename: lib_ncdf.F90

(    1) !-------------------------------------------------------------------------------
(    2) !      OPA NETCDF I/O LIBRARY
(    3) !     ========================
(    4) !
(    5) ! This module is intended to provide a workaround for some of the difficulties,
(    6) ! eccentricities or general insanities of the IOIPSL library that is integrated
(    7) ! with OPA. All output routines in this module write to a single file, rather
(    8) ! than a separate file for each processor (thus avoiding the necessity of
(    9) ! merging files offline). There is no measurable time overhead incurred by
(   10) ! doing this, so we're not sacrificing model speed to avoid merging files.
(   11) ! Similarly, all read routines need only a single file. This library also allows
(   12) ! for easy creation and structuring of new output files. Unlike IOIPSL, it
(   13) ! performs no averaging on the data - these routines simply write or read
(   14) ! the requested data, and nothing else. You also don't need to worry about
(   15) ! opening and closing files, ID numbers, or anything else - that's all handled
(   16) ! for you. This library should work for any number of processors in any layout.
(   17) !
(   18) ! The 2D and 3D variable writing routines attempt to deal with time axes
(   19) ! somewhat 'intelligently'. When a 2D array is passed to ncdf_write, the
(   20) ! corresponding subroutine checks the specified variable for a time dimension.
(   21) ! If one exists, the variable in the netCDF file is, obviously, 3D, even though
(   22) ! only a 2D array was passed. The subroutine determines the correct index on
(   23) ! the time axis to write to, then puts the data there. If no time axis is found,
(   24) ! the array is simply written to the variable as you'd expect. The same sort
(   25) ! of thing happens when a 3D array is passed. I hope that made some sense.
(   26) !
(   27) ! NOTE - I've had issues with the FLUSH subroutine causing segmentation faults
(   28) ! on some systems. If this happens, just comment out any calls to FLUSH and
(   29) ! recompiler. The problem should go away.
(   30) !
(   31) ! Chris Nickerson
(   32) ! October & November 2007
(   33) ! nickersonc@mar.dfo-mpo.gc.ca
(   34) !-------------------------------------------------------------------------------
(   35) 
(   36) MODULE lib_ncdf
(   37)   USE netcdf
(   38)   USE par_oce
(   39)   USE dom_oce
(   40)   USE in_out_manager
(   41)   USE calendar
(   42)   USE daymod
(   43) 
(   44)   IMPLICIT NONE
(   45) 
(   46) !!DB: PRIVATE statement omitted as the include below becomes hidden to USErs






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 2

(   47) # include <mpif.h>
(   48) 
(   49) 
(   50)   ! Error constant definitions
(   51)   INTEGER,PARAMETER :: NCDF_NOERR = 0, &  ! No error, normal return status
(   52)                        NCDF_NFERR = 1, &  ! netCDF-related error occurred
(   53)                        NCDF_MPERR = 2, &  ! MPI-related error occurred
(   54)                        NCDF_ARERR = 3, &  ! Invalid arguments were given
(   55)                        NCDF_OTHER = 4     ! Some other error happened
(   56)   
(   57)   ! Datatype definitions used for creating variables
(   58)   INTEGER,PARAMETER :: NCDF_FLOAT = nf90_float, &
(   59)                        NCDF_DOUBLE = nf90_double
(   60) 
(   61)   ! Turn debugging output on/off (goes to fort.100 - it's VERY verbose)
(   62)   LOGICAL,PARAMETER :: DEBUG_OUT = .TRUE. ! True if you want debugging output,
(   63)                                            ! false otherwise. Debug output is
(   64)                                            ! written to OPA.out in the run dir
(   65)                                            ! Note that debug output is extremely
(   66)                                            ! verbose!
(   67)   
(   68) !!DB 2008.05.22 -- added ice restart
(   69) !!DB 2008.06.26 -- added get_dim_size
(   70)   ! Interface definitions
(   71)   PUBLIC :: ncdf_create_file_u, ncdf_create_file_v, ncdf_create_file_t,  ncdf_write, &
(   72)        ncdf_errstr, ncdf_create_file, ncdf_create_dim, ncdf_create_var, ncdf_put_att, &
(   73)        ncdf_create_restart,  ncdf_create_ice_restart, ncdf_get_dim_size, ncdf_readdate, &
(   74)        ncdf_read_global, ncdf_create_file_aveTSUV, ncdf_create_file_ave
(   75) 
(   76) !  PUBLIC :: output_special             ! special routines called by step.F90 if ave flag is on
(   77) !  PUBLIC :: output_aveTSUV             ! special routines called by step.F90 if ave flag is on
(   78) 
(   79)   
(   80)   INTERFACE ncdf_write
(   81)      MODULE PROCEDURE ncdf_writesv, ncdf_write1d, ncdf_write2d, ncdf_write3d, ncdf_write4d
(   82)   END INTERFACE
(   83) 
(   84)   INTERFACE ncdf_read
(   85)      MODULE PROCEDURE ncdf_readsv, ncdf_read1d, ncdf_read2d, ncdf_read3d, ncdf_read4d
(   86)   END INTERFACE
(   87) 
(   88)   INTERFACE ncdf_read_global
(   89)      MODULE PROCEDURE ncdf_read2d_global, ncdf_read3d_global, ncdf_read4d_global
(   90)   END INTERFACE
(   91) 
(   92)   INTERFACE ncdf_put_att
(   93)      MODULE PROCEDURE ncdf_put_att_int, ncdf_put_att_real, ncdf_put_att_char
(   94)   END INTERFACE
(   95) 
(   96) CONTAINS
(   97) 
(   98)   ! ncdf_create_file_u builds a standard u-grid OPA output file with all the default
(   99)   ! dimensions, variables and attributes
(  100)   SUBROUTINE ncdf_create_file_u(filename, op_type, status)
(  101)     IMPLICIT NONE
(  102)     ! Subroutine argument declarations
(  103)     CHARACTER(LEN=*),INTENT(IN) :: filename
(  104)     CHARACTER(LEN=*),INTENT(IN) :: op_type






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 3

(  105)     INTEGER,INTENT(OUT) :: status
(  106) 
(  107)     ! Local declarations
(  108)     INTEGER :: ncid,    &  ! netCDF file ID
(  109)                varid,   &  ! ID of netCDF variable to be written to
(  110)                nfstat,  &  ! netCDF library call return status
(  111)                mpistat     ! MPI library call return status
(  112)     INTEGER,DIMENSION(1:4) :: dimids
(  113)     INTEGER,DIMENSION(1:6) :: varids
(  114)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  115)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  116)     CHARACTER(LEN=100) :: sec_since    
(  117)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  118)     INTEGER :: int_opp, &              ! Operation interval
(  119)                int_wri                 ! Write interval
(  120)     CHARACTER(LEN=3),PARAMETER :: &
(  121)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  122)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  123)     
(  124)     ! Initializations
(  125)     status = NCDF_NOERR
(  126)     CALL ioget_calendar(cal_type)
(  127)     CALL ioget_timestamp(timestamp)
(  128)     WRITE (UNIT=sec_since, &
(  129)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  130)          &  nyear,nmonth,nday,0, 0, 0
(  131)     WRITE(t_origin, &
(  132)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  133)          &   nyear,months(nmonth),nday,0,0,0
(  134)     int_opp = nwrite * rdt
(  135)     int_wri = nwrite * rdt
(  136) 
(  137)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  138)        WRITE(100,*) 'NCDF DEBUG: Creating default U output file:', filename
(  139)        CALL FLUSH
(  140)     END IF
(  141) 
(  142)     ! Only processor 0 does anything
(  143)     IF(nproc == 0) THEN
(  144)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  145)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Creating file:', filename
(  146)           CALL FLUSH
(  147)        END IF
(  148)        ! Create the file
(  149)        nfstat = nf90_create(filename, nf90_clobber, ncid)
(  150)        IF(nfstat /= nf90_noerr) THEN
(  151)           status = NCDF_NFERR
(  152)           RETURN
(  153)        END IF
(  154)        
(  155)        ! Define dimensions
(  156)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  157)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining dimensions in file:', filename
(  158)           CALL FLUSH
(  159)        END IF
(  160)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  161)        nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
(  162)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 4

(  163)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
(  164)        
(  165)        ! Define variables
(  166)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  167)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining variables in file:', filename
(  168)           CALL FLUSH
(  169)        END IF
(  170)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  171)             (/ dimids(4), dimids(3) /), &
(  172)             varids(1))
(  173)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  174)             (/ dimids(4), dimids(3) /), &
(  175)             varids(2))
(  176)        nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
(  177)             (/ dimids(2) /), &
(  178)             varids(3))
(  179)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  180)             (/ dimids(1) /), &
(  181)             varids(4))
(  182)        nfstat = nf90_def_var(ncid, 'vozocrtx', nf90_float, &
(  183)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  184)             varids(5))
(  185)        nfstat = nf90_def_var(ncid, 'sozotaux', nf90_float, &
(  186)             (/ dimids(4), dimids(3), dimids(1) /), &
(  187)             varids(6))
(  188)        
(  189)        ! Add attributes
(  190)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  191)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Writing attributes in file:', filename
(  192)           CALL FLUSH
(  193)        END IF
(  194)        ! nav_lon
(  195)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  196)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  197)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  198)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  199)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  200) 
(  201)        ! nav_lat
(  202)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  203)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  204)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  205)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  206)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  207) 
(  208)        ! depthu
(  209)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  210)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  211)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
(  212)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
(  213)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
(  214)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')
(  215) 
(  216)        ! time_counter
(  217)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  218)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  219)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  220)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 5

(  221)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
(  222) 
(  223)        ! vozocrtx
(  224)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
(  225)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
(  226)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
(  227)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
(  228)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Zonal Current')
(  229)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vozocrtx')
(  230)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
(  231)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
(  232)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
(  233)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
(  234)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')
(  235) 
(  236)        ! sozotaux
(  237)        nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
(  238)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
(  239)        nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
(  240)        nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
(  241)        nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along i-axis')
(  242)        nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sozotaux')
(  243)        nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
(  244)        nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')
(  245)        nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
(  246)        nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
(  247)        nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')
(  248) 
(  249)        ! global
(  250)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  251)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
(  252)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
(  253)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  254)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  255)       
(  256)        
(  257)        ! Close file
(  258)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  259)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Closing file:', filename
(  260)           CALL FLUSH
(  261)        END IF
(  262)        nfstat = nf90_close(ncid)
(  263)        IF(nfstat /= nf90_noerr) THEN
(  264)           status = NCDF_NFERR
(  265)           RETURN
(  266)        END IF
(  267)     END IF
(  268) 
(  269)     ! Sync up processors before returning from subroutine
(  270)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
(  271)     IF(mpistat /= 0) THEN
(  272)        status = NCDF_MPERR
(  273)        RETURN
(  274)     END IF
(  275)     
(  276)   END SUBROUTINE ncdf_create_file_u
(  277) 
(  278)   ! ncdf_create_file_v builds a standard v-grid OPA output file with all the default






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 6

(  279)   ! dimensions, variables and attributes
(  280)   SUBROUTINE ncdf_create_file_v(filename, op_type, status)
(  281)     IMPLICIT NONE
(  282)     ! Subroutine argument declarations
(  283)     CHARACTER(LEN=*),INTENT(IN) :: filename
(  284)     CHARACTER(LEN=*),INTENT(IN) :: op_type
(  285)     INTEGER,INTENT(OUT) :: status
(  286) 
(  287)     ! Local declarations
(  288)     INTEGER :: ncid,    &  ! netCDF file ID
(  289)                varid,   &  ! ID of netCDF variable to be written to
(  290)                nfstat,  &  ! netCDF library call return status
(  291)                mpistat     ! MPI library call return status
(  292)     INTEGER,DIMENSION(1:4) :: dimids
(  293)     INTEGER,DIMENSION(1:6) :: varids
(  294)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  295)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  296)     CHARACTER(LEN=100) :: sec_since
(  297)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  298)     INTEGER :: int_opp, &              ! Operation interval
(  299)                int_wri                 ! Write interval
(  300)     CHARACTER(LEN=3),PARAMETER :: &
(  301)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  302)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  303)     
(  304)     ! Initializations
(  305)     status = NCDF_NOERR
(  306)     CALL ioget_calendar(cal_type)
(  307)     CALL ioget_timestamp(timestamp)
(  308)     WRITE (UNIT=sec_since, &
(  309)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  310)          &  nyear,nmonth,nday,0, 0, 0
(  311)     WRITE(t_origin, &
(  312)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  313)          &   nyear,months(nmonth),nday,0,0,0
(  314)     int_opp = nwrite * rdt
(  315)     int_wri = nwrite * rdt
(  316) 
(  317)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  318)        WRITE(100,*) 'NCDF DEBUG: Creating default V output file:', filename
(  319)        CALL FLUSH
(  320)     END IF
(  321) 
(  322)     ! Only processor 0 does anything
(  323)     IF(nproc == 0) THEN
(  324)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  325)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Creating file:', filename
(  326)           CALL FLUSH
(  327)        END IF
(  328)        ! Create the file
(  329)        nfstat = nf90_create(filename, nf90_clobber, ncid)
(  330)        IF(nfstat /= nf90_noerr) THEN
(  331)           status = NCDF_NFERR
(  332)           RETURN
(  333)        END IF
(  334)        
(  335)        ! Define dimensions
(  336)        IF(DEBUG_OUT .EQV. .TRUE.) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 7

(  337)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Defining dimensions in file:', filename
(  338)           CALL FLUSH
(  339)        END IF
(  340)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  341)        nfstat = nf90_def_dim(ncid, 'depthv', jpkdta, dimids(2))
(  342)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
(  343)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
(  344)        
(  345)        ! Define variables
(  346)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  347)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Defining variables in file:', filename
(  348)           CALL FLUSH
(  349)        END IF
(  350)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  351)             (/ dimids(4), dimids(3) /), &
(  352)             varids(1))
(  353)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  354)             (/ dimids(4), dimids(3) /), &
(  355)             varids(2))
(  356)        nfstat = nf90_def_var(ncid, 'depthv', nf90_float, &
(  357)             (/ dimids(2) /), &
(  358)             varids(3))
(  359)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  360)             (/ dimids(1) /), &
(  361)             varids(4))
(  362)        nfstat = nf90_def_var(ncid, 'vomecrty', nf90_float, &
(  363)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  364)             varids(5))
(  365)        nfstat = nf90_def_var(ncid, 'sometauy', nf90_float, &
(  366)             (/ dimids(4), dimids(3), dimids(1) /), &
(  367)             varids(6))
(  368)        
(  369)        ! Add attributes
(  370)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  371)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Writing attributes in file:', filename
(  372)           CALL FLUSH
(  373)        END IF
(  374)        ! nav_lon
(  375)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  376)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  377)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  378)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  379)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  380) 
(  381)        ! nav_lat
(  382)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  383)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  384)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  385)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  386)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  387) 
(  388)        ! depthu
(  389)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  390) 
(  391)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  392)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
(  393)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
(  394)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthv')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 8

(  395)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical V levels')
(  396) 
(  397)        ! time_counter
(  398)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  399)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  400)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  401)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
(  402)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
(  403) 
(  404)        ! vomecrty
(  405)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
(  406)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
(  407)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
(  408)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
(  409)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Meridional Current')
(  410)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vomecrty')
(  411)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
(  412)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
(  413)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
(  414)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
(  415)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthv nav_lat nav_lon')
(  416) 
(  417)        ! sometauy
(  418)        nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
(  419)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
(  420)        nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
(  421)        nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
(  422)        nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along j-axis')
(  423)        nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sometauy')
(  424)        nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
(  425)        nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')
(  426)        nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
(  427)        nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
(  428)        nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')
(  429) 
(  430)        ! global
(  431)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  432)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
(  433)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
(  434)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  435)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  436)       
(  437)        
(  438)        ! Close file
(  439)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  440)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Closing file:', filename
(  441)           CALL FLUSH
(  442)        END IF
(  443)        nfstat = nf90_close(ncid)
(  444)        IF(nfstat /= nf90_noerr) THEN
(  445)           status = NCDF_NFERR
(  446)           RETURN
(  447)        END IF
(  448)     END IF
(  449) 
(  450)     ! Sync up processors before returning from subroutine
(  451)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
(  452)     IF(mpistat /= 0) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 9

(  453)        status = NCDF_MPERR
(  454)        RETURN
(  455)     END IF
(  456)     
(  457)   END SUBROUTINE ncdf_create_file_v
(  458) 
(  459)   ! ncdf_create_file_t builds a standard t-grid OPA output file with all the default
(  460)   ! dimensions, variables and attributes
(  461)   SUBROUTINE ncdf_create_file_t(filename, op_type, status)
(  462)     IMPLICIT NONE
(  463)     ! Subroutine argument declarations
(  464)     CHARACTER(LEN=*),INTENT(IN) :: filename
(  465)     CHARACTER(LEN=*),INTENT(IN) :: op_type
(  466)     INTEGER,INTENT(OUT) :: status
(  467) 
(  468)     ! Local declarations
(  469)     INTEGER :: ncid,    &  ! netCDF file ID
(  470)                varid,   &  ! ID of netCDF variable to be written to
(  471)                nfstat,  &  ! netCDF library call return status
(  472)                mpistat     ! MPI library call return status
(  473)     INTEGER,DIMENSION(1:4) :: dimids
(  474)     INTEGER,DIMENSION(1:22) :: varids
(  475)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  476)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  477)     CHARACTER(LEN=100) :: sec_since
(  478)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  479)     INTEGER :: int_opp, &              ! Operation interval
(  480)                int_wri                 ! Write interval
(  481)     CHARACTER(LEN=3),PARAMETER :: &
(  482)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  483)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  484)     
(  485)     ! Initializations
(  486)     status = NCDF_NOERR
(  487)     CALL ioget_calendar(cal_type)
(  488)     CALL ioget_timestamp(timestamp)
(  489)     WRITE (UNIT=sec_since, &
(  490)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  491)          &  nyear,nmonth,nday,0, 0, 0
(  492)     WRITE(t_origin, &
(  493)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  494)          &   nyear,months(nmonth),nday,0,0,0
(  495)     int_opp = nwrite * rdt
(  496)     int_wri = nwrite * rdt
(  497) 
(  498)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  499)        WRITE(100,*) 'NCDF DEBUG: Creating default T output file:', filename
(  500)        CALL FLUSH
(  501)     END IF
(  502) 
(  503)     ! Only processor 0 does anything
(  504)     IF(nproc == 0) THEN
(  505)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  506)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Creating file:', filename
(  507)           CALL FLUSH
(  508)        END IF
(  509)        ! Create the file
(  510)        nfstat = nf90_create(filename, nf90_clobber, ncid)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 10

(  511)        
(  512)        ! Define dimensions
(  513)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  514)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Defining dimensions in file:', filename
(  515)           CALL FLUSH
(  516)        END IF
(  517)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  518)        nfstat = nf90_def_dim(ncid, 'jpkdta', jpkdta, dimids(2))
(  519)        nfstat = nf90_def_dim(ncid, 'jpjdta', jpjdta, dimids(3))
(  520)        nfstat = nf90_def_dim(ncid, 'jpidta', jpidta, dimids(4))
(  521)        
(  522)        ! Define variables
(  523)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  524)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Defining variables in file:', filename
(  525)           CALL FLUSH
(  526)        END IF
(  527)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  528)             (/ dimids(4), dimids(3) /), &
(  529)             varids(1))
(  530)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  531)             (/ dimids(4), dimids(3) /), &
(  532)             varids(2))
(  533)        nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &
(  534)             (/ dimids(2) /), &
(  535)             varids(3))
(  536)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  537)             (/ dimids(1) /), &
(  538)             varids(4))
(  539)        nfstat = nf90_def_var(ncid, 'votemper', nf90_float, &
(  540)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  541)             varids(5))
(  542)        nfstat = nf90_def_var(ncid, 'vosaline', nf90_float, &
(  543)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  544)             varids(6))
(  545)        nfstat = nf90_def_var(ncid, 'sosstsst', nf90_float, &
(  546)             (/ dimids(4), dimids(3), dimids(1) /), &
(  547)             varids(7))
(  548)        nfstat = nf90_def_var(ncid, 'sosaline', nf90_float, &
(  549)             (/ dimids(4), dimids(3), dimids(1) /), &
(  550)             varids(8))
(  551)        nfstat = nf90_def_var(ncid, 'sossheig', nf90_float, &
(  552)             (/ dimids(4), dimids(3), dimids(1) /), &
(  553)             varids(9))
(  554)        nfstat = nf90_def_var(ncid, 'sowaflup', nf90_float, &
(  555)             (/ dimids(4), dimids(3), dimids(1) /), &
(  556)             varids(10))
(  557)        nfstat = nf90_def_var(ncid, 'sorunoff', nf90_float, &
(  558)             (/ dimids(4), dimids(3), dimids(1) /), &
(  559)             varids(11))
(  560)        nfstat = nf90_def_var(ncid, 'sowaflcd', nf90_float, &
(  561)             (/ dimids(4), dimids(3), dimids(1) /), &
(  562)             varids(12))
(  563)        nfstat = nf90_def_var(ncid, 'sosalflx', nf90_float, &
(  564)             (/ dimids(4), dimids(3), dimids(1) /), &
(  565)             varids(13))
(  566)        nfstat = nf90_def_var(ncid, 'sohefldo', nf90_float, &
(  567)             (/ dimids(4), dimids(3), dimids(1) /), &
(  568)             varids(14))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 11

(  569)        nfstat = nf90_def_var(ncid, 'soshfldo', nf90_float, &
(  570)             (/ dimids(4), dimids(3), dimids(1) /), &
(  571)             varids(15))
(  572)        nfstat = nf90_def_var(ncid, 'somxl010', nf90_float, &
(  573)             (/ dimids(4), dimids(3), dimids(1) /), &
(  574)             varids(16))
(  575)        nfstat = nf90_def_var(ncid, 'somixhgt', nf90_float, &
(  576)             (/ dimids(4), dimids(3), dimids(1) /), &
(  577)             varids(17))
(  578)        nfstat = nf90_def_var(ncid, 'soicecov', nf90_float, &
(  579)             (/ dimids(4), dimids(3), dimids(1) /), &
(  580)             varids(18))
(  581)        nfstat = nf90_def_var(ncid, 'sohefldp', nf90_float, &
(  582)             (/ dimids(4), dimids(3), dimids(1) /), &
(  583)             varids(19))
(  584)        nfstat = nf90_def_var(ncid, 'sowafldp', nf90_float, &
(  585)             (/ dimids(4), dimids(3), dimids(1) /), &
(  586)             varids(20))
(  587)        nfstat = nf90_def_var(ncid, 'sosafldp', nf90_float, &
(  588)             (/ dimids(4), dimids(3), dimids(1) /), &
(  589)             varids(21))
(  590)        nfstat = nf90_def_var(ncid, 'sobowlin', nf90_float, &
(  591)             (/ dimids(4), dimids(3) /), &
(  592)             varids(22))
(  593)        
(  594)        ! Add attributes
(  595)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  596)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Writing attributes in file:', filename
(  597)           CALL FLUSH
(  598)        END IF
(  599)        ! nav_lon
(  600)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  601)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  602)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  603)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  604)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  605) 
(  606)        ! nav_lat
(  607)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  608)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  609)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  610)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  611)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  612) 
(  613)        ! deptht
(  614)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  615)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  616)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
(  617)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
(  618)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
(  619)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical T levels')
(  620) 
(  621)        ! time_counter
(  622)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  623)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  624)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  625)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
(  626)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 12

(  627) 
(  628)        ! votemper
(  629)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'C')
(  630)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
(  631)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
(  632)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
(  633)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Temperature')
(  634)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'votemper')
(  635)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
(  636)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
(  637)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
(  638)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
(  639)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter deptht nav_lat nav_lon')
(  640) 
(  641)        ! vosaline
(  642)        nfstat = nf90_put_att(ncid, varids(6), 'units', 'PSU')
(  643)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
(  644)        nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
(  645)        nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Salinity')
(  646)        nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'vosaline')
(  647)        nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
(  648)        nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TZYX')
(  649)        nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
(  650)        nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
(  651)        nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter deptht nav_lat nav_lon')
(  652) 
(  653)        ! sosstsst
(  654)        nfstat = nf90_put_att(ncid, varids(7), 'units', 'C')
(  655)        nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
(  656)        nfstat = nf90_put_att(ncid, varids(7), 'valid_min', 1.000000E20 )
(  657)        nfstat = nf90_put_att(ncid, varids(7), 'valid_max', -1.000000E20)
(  658)        nfstat = nf90_put_att(ncid, varids(7), 'long_name', 'Sea Surface temperature')
(  659)        nfstat = nf90_put_att(ncid, varids(7), 'short_name', 'sosstsst')
(  660)        nfstat = nf90_put_att(ncid, varids(7), 'online_operation', TRIM(op_type))
(  661)        nfstat = nf90_put_att(ncid, varids(7), 'axis', 'TYX')
(  662)        nfstat = nf90_put_att(ncid, varids(7), 'interval_operation', int_opp)
(  663)        nfstat = nf90_put_att(ncid, varids(7), 'interval_write', int_wri)
(  664)        nfstat = nf90_put_att(ncid, varids(7), 'associate', 'time_counter nav_lat nav_lon')
(  665) 
(  666)        ! sosaline
(  667)        nfstat = nf90_put_att(ncid, varids(8), 'units', 'PSU')
(  668)        nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
(  669)        nfstat = nf90_put_att(ncid, varids(8), 'valid_min', 1.000000E20 )
(  670)        nfstat = nf90_put_att(ncid, varids(8), 'valid_max', -1.000000E20)
(  671)        nfstat = nf90_put_att(ncid, varids(8), 'long_name', 'Sea Surface Salinity')
(  672)        nfstat = nf90_put_att(ncid, varids(8), 'short_name', 'sosaline')
(  673)        nfstat = nf90_put_att(ncid, varids(8), 'online_operation', TRIM(op_type))
(  674)        nfstat = nf90_put_att(ncid, varids(8), 'axis', 'TYX')
(  675)        nfstat = nf90_put_att(ncid, varids(8), 'interval_operation', int_opp)
(  676)        nfstat = nf90_put_att(ncid, varids(8), 'interval_write', int_wri)
(  677)        nfstat = nf90_put_att(ncid, varids(8), 'associate', 'time_counter nav_lat nav_lon')
(  678) 
(  679)        ! sossheig
(  680)        nfstat = nf90_put_att(ncid, varids(9), 'units', 'm')
(  681)        nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)
(  682)        nfstat = nf90_put_att(ncid, varids(9), 'valid_min', 1.000000E20 )
(  683)        nfstat = nf90_put_att(ncid, varids(9), 'valid_max', -1.000000E20)
(  684)        nfstat = nf90_put_att(ncid, varids(9), 'long_name', 'Sea Surface Height')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 13

(  685)        nfstat = nf90_put_att(ncid, varids(9), 'short_name', 'sossheig')
(  686)        nfstat = nf90_put_att(ncid, varids(9), 'online_operation', TRIM(op_type))
(  687)        nfstat = nf90_put_att(ncid, varids(9), 'axis', 'TYX')
(  688)        nfstat = nf90_put_att(ncid, varids(9), 'interval_operation', int_opp)
(  689)        nfstat = nf90_put_att(ncid, varids(9), 'interval_write', int_wri)
(  690)        nfstat = nf90_put_att(ncid, varids(9), 'associate', 'time_counter nav_lat nav_lon')
(  691)        
(  692)        ! sowaflup
(  693)        nfstat = nf90_put_att(ncid, varids(10), 'units', 'Kg/m2/s')
(  694)        nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)
(  695)        nfstat = nf90_put_att(ncid, varids(10), 'valid_min', 1.000000E20 )
(  696)        nfstat = nf90_put_att(ncid, varids(10), 'valid_max', -1.000000E20)
(  697)        nfstat = nf90_put_att(ncid, varids(10), 'long_name', 'Net Upward Water Flux')
(  698)        nfstat = nf90_put_att(ncid, varids(10), 'short_name', 'sowaflup')
(  699)        nfstat = nf90_put_att(ncid, varids(10), 'online_operation', TRIM(op_type))
(  700)        nfstat = nf90_put_att(ncid, varids(10), 'axis', 'TYX')
(  701)        nfstat = nf90_put_att(ncid, varids(10), 'interval_operation', int_opp)
(  702)        nfstat = nf90_put_att(ncid, varids(10), 'interval_write', int_wri)
(  703)        nfstat = nf90_put_att(ncid, varids(10), 'associate', 'time_counter nav_lat nav_lon')
(  704)        
(  705)        ! sorunoff
(  706)        nfstat = nf90_put_att(ncid, varids(11), 'units', 'Kg/m2/s')
(  707)        nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)
(  708)        nfstat = nf90_put_att(ncid, varids(11), 'valid_min', 1.000000E20 )
(  709)        nfstat = nf90_put_att(ncid, varids(11), 'valid_max', -1.000000E20)
(  710)        nfstat = nf90_put_att(ncid, varids(11), 'long_name', 'Runoffs')
(  711)        nfstat = nf90_put_att(ncid, varids(11), 'short_name', 'sorunoff')
(  712)        nfstat = nf90_put_att(ncid, varids(11), 'online_operation', TRIM(op_type))
(  713)        nfstat = nf90_put_att(ncid, varids(11), 'axis', 'TYX')
(  714)        nfstat = nf90_put_att(ncid, varids(11), 'interval_operation', int_opp)
(  715)        nfstat = nf90_put_att(ncid, varids(11), 'interval_write', int_wri)
(  716)        nfstat = nf90_put_att(ncid, varids(11), 'associate', 'time_counter nav_lat nav_lon')
(  717)        
(  718)        ! sowaflcd
(  719)        nfstat = nf90_put_att(ncid, varids(12), 'units', 'kg/m2/s')
(  720)        nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)
(  721)        nfstat = nf90_put_att(ncid, varids(12), 'valid_min', 1.000000E20 )
(  722)        nfstat = nf90_put_att(ncid, varids(12), 'valid_max', -1.000000E20)
(  723)        nfstat = nf90_put_att(ncid, varids(12), 'long_name', 'concentration/dilution water flux')
(  724)        nfstat = nf90_put_att(ncid, varids(12), 'short_name', 'sowaflcd')
(  725)        nfstat = nf90_put_att(ncid, varids(12), 'online_operation', TRIM(op_type))
(  726)        nfstat = nf90_put_att(ncid, varids(12), 'axis', 'TYX')
(  727)        nfstat = nf90_put_att(ncid, varids(12), 'interval_operation', int_opp)
(  728)        nfstat = nf90_put_att(ncid, varids(12), 'interval_write', int_wri)
(  729)        nfstat = nf90_put_att(ncid, varids(12), 'associate', 'time_counter nav_lat nav_lon')
(  730)        
(  731)        ! sosalflx
(  732)        nfstat = nf90_put_att(ncid, varids(13), 'units', 'Kg/m2/s')
(  733)        nfstat = nf90_put_att(ncid, varids(13), 'missing_value', 1.000000E20)
(  734)        nfstat = nf90_put_att(ncid, varids(13), 'valid_min', 1.000000E20 )
(  735)        nfstat = nf90_put_att(ncid, varids(13), 'valid_max', -1.000000E20)
(  736)        nfstat = nf90_put_att(ncid, varids(13), 'long_name', 'Surface Salt Flux')
(  737)        nfstat = nf90_put_att(ncid, varids(13), 'short_name', 'sosalflx')
(  738)        nfstat = nf90_put_att(ncid, varids(13), 'online_operation', TRIM(op_type))
(  739)        nfstat = nf90_put_att(ncid, varids(13), 'axis', 'TYX')
(  740)        nfstat = nf90_put_att(ncid, varids(13), 'interval_operation', int_opp)
(  741)        nfstat = nf90_put_att(ncid, varids(13), 'interval_write', int_wri)
(  742)        nfstat = nf90_put_att(ncid, varids(13), 'associate', 'time_counter nav_lat nav_lon')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 14

(  743) 
(  744)        ! sohefldo
(  745)        nfstat = nf90_put_att(ncid, varids(14), 'units', 'W/m2')
(  746)        nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)
(  747)        nfstat = nf90_put_att(ncid, varids(14), 'valid_min', 1.000000E20 )
(  748)        nfstat = nf90_put_att(ncid, varids(14), 'valid_max', -1.000000E20)
(  749)        nfstat = nf90_put_att(ncid, varids(14), 'long_name', 'Net Downward Heat Flux')
(  750)        nfstat = nf90_put_att(ncid, varids(14), 'short_name', 'sohefldo')
(  751)        nfstat = nf90_put_att(ncid, varids(14), 'online_operation', TRIM(op_type))
(  752)        nfstat = nf90_put_att(ncid, varids(14), 'axis', 'TYX')
(  753)        nfstat = nf90_put_att(ncid, varids(14), 'interval_operation', int_opp)
(  754)        nfstat = nf90_put_att(ncid, varids(14), 'interval_write', int_wri)
(  755)        nfstat = nf90_put_att(ncid, varids(14), 'associate', 'time_counter nav_lat nav_lon')
(  756) 
(  757)        ! soshfldo
(  758)        nfstat = nf90_put_att(ncid, varids(15), 'units', 'W/m2')
(  759)        nfstat = nf90_put_att(ncid, varids(15), 'missing_value', 1.000000E20)
(  760)        nfstat = nf90_put_att(ncid, varids(15), 'valid_min', 1.000000E20 )
(  761)        nfstat = nf90_put_att(ncid, varids(15), 'valid_max', -1.000000E20)
(  762)        nfstat = nf90_put_att(ncid, varids(15), 'long_name', 'Shortwave Radiation')
(  763)        nfstat = nf90_put_att(ncid, varids(15), 'short_name', 'soshfldo')
(  764)        nfstat = nf90_put_att(ncid, varids(15), 'online_operation', TRIM(op_type))
(  765)        nfstat = nf90_put_att(ncid, varids(15), 'axis', 'TYX')
(  766)        nfstat = nf90_put_att(ncid, varids(15), 'interval_operation', int_opp)
(  767)        nfstat = nf90_put_att(ncid, varids(15), 'interval_write', int_wri)
(  768)        nfstat = nf90_put_att(ncid, varids(15), 'associate', 'time_counter nav_lat nav_lon')
(  769) 
(  770)        ! soml010
(  771)        nfstat = nf90_put_att(ncid, varids(16), 'units', 'm')
(  772)        nfstat = nf90_put_att(ncid, varids(16), 'missing_value', 1.000000E20)
(  773)        nfstat = nf90_put_att(ncid, varids(16), 'valid_min', 1.000000E20 )
(  774)        nfstat = nf90_put_att(ncid, varids(16), 'valid_max', -1.000000E20)
(  775)        nfstat = nf90_put_att(ncid, varids(16), 'long_name', 'Mixed Layer Depth 0.01')
(  776)        nfstat = nf90_put_att(ncid, varids(16), 'short_name', 'soml010')
(  777)        nfstat = nf90_put_att(ncid, varids(16), 'online_operation', TRIM(op_type))
(  778)        nfstat = nf90_put_att(ncid, varids(16), 'axis', 'TYX')
(  779)        nfstat = nf90_put_att(ncid, varids(16), 'interval_operation', int_opp)
(  780)        nfstat = nf90_put_att(ncid, varids(16), 'interval_write', int_wri)
(  781)        nfstat = nf90_put_att(ncid, varids(16), 'associate', 'time_counter nav_lat nav_lon')
(  782) 
(  783)        ! somixhgt
(  784)        nfstat = nf90_put_att(ncid, varids(17), 'units', 'm')
(  785)        nfstat = nf90_put_att(ncid, varids(17), 'missing_value', 1.000000E20)
(  786)        nfstat = nf90_put_att(ncid, varids(17), 'valid_min', 1.000000E20 )
(  787)        nfstat = nf90_put_att(ncid, varids(17), 'valid_max', -1.000000E20)
(  788)        nfstat = nf90_put_att(ncid, varids(17), 'long_name', 'Turbocline Depth')
(  789)        nfstat = nf90_put_att(ncid, varids(17), 'short_name', 'somixhgt')
(  790)        nfstat = nf90_put_att(ncid, varids(17), 'online_operation', TRIM(op_type))
(  791)        nfstat = nf90_put_att(ncid, varids(17), 'axis', 'TYX')
(  792)        nfstat = nf90_put_att(ncid, varids(17), 'interval_operation', int_opp)
(  793)        nfstat = nf90_put_att(ncid, varids(17), 'interval_write', int_wri)
(  794)        nfstat = nf90_put_att(ncid, varids(17), 'associate', 'time_counter nav_lat nav_lon')
(  795) 
(  796)        ! soicecov
(  797)        nfstat = nf90_put_att(ncid, varids(18), 'units', '[0,1]')
(  798)        nfstat = nf90_put_att(ncid, varids(18), 'missing_value', 1.000000E20)
(  799)        nfstat = nf90_put_att(ncid, varids(18), 'valid_min', 1.000000E20 )
(  800)        nfstat = nf90_put_att(ncid, varids(18), 'valid_max', -1.000000E20)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 15

(  801)        nfstat = nf90_put_att(ncid, varids(18), 'long_name', 'Ice Cover')
(  802)        nfstat = nf90_put_att(ncid, varids(18), 'short_name', 'soicecov')
(  803)        nfstat = nf90_put_att(ncid, varids(18), 'online_operation', TRIM(op_type))
(  804)        nfstat = nf90_put_att(ncid, varids(18), 'axis', 'TYX')
(  805)        nfstat = nf90_put_att(ncid, varids(18), 'interval_operation', int_opp)
(  806)        nfstat = nf90_put_att(ncid, varids(18), 'interval_write', int_wri)
(  807)        nfstat = nf90_put_att(ncid, varids(18), 'associate', 'time_counter nav_lat nav_lon')
(  808) 
(  809)        ! sohefldp
(  810)        nfstat = nf90_put_att(ncid, varids(19), 'units', 'W/m2')
(  811)        nfstat = nf90_put_att(ncid, varids(19), 'missing_value', 1.000000E20)
(  812)        nfstat = nf90_put_att(ncid, varids(19), 'valid_min', 1.000000E20 )
(  813)        nfstat = nf90_put_att(ncid, varids(19), 'valid_max', -1.000000E20)
(  814)        nfstat = nf90_put_att(ncid, varids(19), 'long_name', 'Surface Heat Flux: Damping')
(  815)        nfstat = nf90_put_att(ncid, varids(19), 'short_name', 'sohefldp')
(  816)        nfstat = nf90_put_att(ncid, varids(19), 'online_operation', TRIM(op_type))
(  817)        nfstat = nf90_put_att(ncid, varids(19), 'axis', 'TYX')
(  818)        nfstat = nf90_put_att(ncid, varids(19), 'interval_operation', int_opp)
(  819)        nfstat = nf90_put_att(ncid, varids(19), 'interval_write', int_wri)
(  820)        nfstat = nf90_put_att(ncid, varids(19), 'associate', 'time_counter nav_lat nav_lon')
(  821) 
(  822)        ! sowafldp
(  823)        nfstat = nf90_put_att(ncid, varids(20), 'units', 'Kg/m2/s')
(  824)        nfstat = nf90_put_att(ncid, varids(20), 'missing_value', 1.000000E20)
(  825)        nfstat = nf90_put_att(ncid, varids(20), 'valid_min', 1.000000E20 )
(  826)        nfstat = nf90_put_att(ncid, varids(20), 'valid_max', -1.000000E20)
(  827)        nfstat = nf90_put_att(ncid, varids(20), 'long_name', 'Surface Water Flux: Damping')
(  828)        nfstat = nf90_put_att(ncid, varids(20), 'short_name', 'sowafldp')
(  829)        nfstat = nf90_put_att(ncid, varids(20), 'online_operation', TRIM(op_type))
(  830)        nfstat = nf90_put_att(ncid, varids(20), 'axis', 'TYX')
(  831)        nfstat = nf90_put_att(ncid, varids(20), 'interval_operation', int_opp)
(  832)        nfstat = nf90_put_att(ncid, varids(20), 'interval_write', int_wri)
(  833)        nfstat = nf90_put_att(ncid, varids(20), 'associate', 'time_counter nav_lat nav_lon')
(  834) 
(  835)        ! sosafldp
(  836)        nfstat = nf90_put_att(ncid, varids(21), 'units', 'Kg/m2/s')
(  837)        nfstat = nf90_put_att(ncid, varids(21), 'missing_value', 1.000000E20)
(  838)        nfstat = nf90_put_att(ncid, varids(21), 'valid_min', 1.000000E20 )
(  839)        nfstat = nf90_put_att(ncid, varids(21), 'valid_max', -1.000000E20)
(  840)        nfstat = nf90_put_att(ncid, varids(21), 'long_name', 'Surface salt flux: damping')
(  841)        nfstat = nf90_put_att(ncid, varids(21), 'short_name', 'sosafldp')
(  842)        nfstat = nf90_put_att(ncid, varids(21), 'online_operation', TRIM(op_type))
(  843)        nfstat = nf90_put_att(ncid, varids(21), 'axis', 'TYX')
(  844)        nfstat = nf90_put_att(ncid, varids(21), 'interval_operation', int_opp)
(  845)        nfstat = nf90_put_att(ncid, varids(21), 'interval_write', int_wri)
(  846)        nfstat = nf90_put_att(ncid, varids(21), 'associate', 'time_counter nav_lat nav_lon')
(  847) 
(  848)        ! sobowlin
(  849)        nfstat = nf90_put_att(ncid, varids(22), 'units', 'W-point')
(  850)        nfstat = nf90_put_att(ncid, varids(22), 'missing_value', 1.000000E20)
(  851)        nfstat = nf90_put_att(ncid, varids(22), 'valid_min', 1.000000E20 )
(  852)        nfstat = nf90_put_att(ncid, varids(22), 'valid_max', -1.000000E20)
(  853)        nfstat = nf90_put_att(ncid, varids(22), 'long_name', 'Bowl Index')
(  854)        nfstat = nf90_put_att(ncid, varids(22), 'short_name', 'sobowlin')
(  855)        nfstat = nf90_put_att(ncid, varids(22), 'online_operation', 'l_max(only(x))')
(  856)        nfstat = nf90_put_att(ncid, varids(22), 'axis', 'TYX')
(  857)        nfstat = nf90_put_att(ncid, varids(22), 'associate', 'time_counter nav_lat nav_lon')
(  858) 






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 16

(  859)        ! global
(  860)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  861)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
(  862)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
(  863)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  864)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  865)       
(  866)        
(  867)        ! Close file
(  868)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  869)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_t - Closing file:', filename
(  870)           CALL FLUSH
(  871)        END IF
(  872)        nfstat = nf90_close(ncid)
(  873)        IF(nfstat /= nf90_noerr) THEN
(  874)           status = NCDF_NFERR
(  875)           RETURN
(  876)        END IF
(  877)     END IF
(  878) 
(  879)     ! Sync up processors before returning from subroutine
(  880)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
(  881)     IF(mpistat /= 0) THEN
(  882)        status = NCDF_MPERR
(  883)        RETURN
(  884)     END IF
(  885)     
(  886)   END SUBROUTINE ncdf_create_file_t
(  887) 
(  888)   ! ncdf_create_file_ave builds an OPA output file with a single variable
(  889)   ! whose name is specified by VARNAME. Other than that, it is similar to a
(  890)   ! U- or V-grid file
(  891)  SUBROUTINE ncdf_create_file_ave(filename, varname, op_type, status)
(  892)     IMPLICIT NONE
(  893)     ! Subroutine argument declarations
(  894)     CHARACTER(LEN=*),INTENT(IN) :: filename
(  895)     CHARACTER(LEN=*),INTENT(IN) :: op_type, varname
(  896)     INTEGER,INTENT(OUT) :: status
(  897) 
(  898)     ! Local declarations
(  899)     INTEGER :: ncid,    &  ! netCDF file ID
(  900)                varid,   &  ! ID of netCDF variable to be written to
(  901)                nfstat,  &  ! netCDF library call return status
(  902)                mpistat     ! MPI library call return status
(  903)     INTEGER,DIMENSION(1:4) :: dimids
(  904)     INTEGER,DIMENSION(1:5) :: varids
(  905)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  906)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  907)     CHARACTER(LEN=100) :: sec_since    
(  908)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  909)     INTEGER :: int_opp, &              ! Operation interval
(  910)                int_wri                 ! Write interval
(  911)     INTEGER :: varnum
(  912)     CHARACTER(LEN=3),PARAMETER :: &
(  913)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  914)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  915)     
(  916)     ! Initializations






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 17

(  917)     status = NCDF_NOERR
(  918)     CALL ioget_calendar(cal_type)
(  919)     CALL ioget_timestamp(timestamp)
(  920)     WRITE (UNIT=sec_since, &
(  921)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  922)          &  nyear,nmonth,nday,0, 0, 0
(  923)     WRITE(t_origin, &
(  924)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  925)          &   nyear,months(nmonth),nday,0,0,0
(  926)     int_opp = nwrite * rdt
(  927)     int_wri = nwrite * rdt
(  928) 
(  929)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  930)        WRITE(100,*) 'NCDF DEBUG: Creating output file:', filename
(  931)        CALL FLUSH
(  932)     END IF
(  933) 
(  934)     ! Only processor 0 does anything
(  935)     IF(nproc == 0) THEN
(  936)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  937)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Creating file:', filename
(  938)           CALL FLUSH
(  939)        END IF
(  940)        ! Create the file
(  941)        nfstat = nf90_create(filename, nf90_clobber, ncid)
(  942)        IF(nfstat /= nf90_noerr) THEN
(  943)           status = NCDF_NFERR
(  944)           RETURN
(  945)        END IF
(  946)        
(  947)        ! Define dimensions
(  948)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  949)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Defining dimensions in file:', filename
(  950)           CALL FLUSH
(  951)        END IF
(  952)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  953)        nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
(  954)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
(  955)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
(  956)        
(  957)        ! Define variables
(  958)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  959)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Defining variables in file:', filename
(  960)           CALL FLUSH
(  961)        END IF
(  962)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  963)             (/ dimids(4), dimids(3) /), &
(  964)             varids(1))
(  965)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  966)             (/ dimids(4), dimids(3) /), &
(  967)             varids(2))
(  968)        nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
(  969)             (/ dimids(2) /), &
(  970)             varids(3))
(  971)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  972)             (/ dimids(1) /), &
(  973)             varids(4))
(  974)        nfstat = nf90_def_var(ncid, varname, nf90_float, &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 18

(  975)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  976)             varids(5))
(  977)        
(  978)        ! Add attributes
(  979)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  980)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Writing attributes in file:', filename
(  981)           CALL FLUSH
(  982)        END IF
(  983)        ! nav_lon
(  984)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  985)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  986)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  987)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  988)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  989) 
(  990)        ! nav_lat
(  991)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  992)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  993)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  994)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  995)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  996) 
(  997)        ! depthu
(  998)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  999)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
( 1000)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
( 1001)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
( 1002)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
( 1003)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')
( 1004) 
( 1005)        ! time_counter
( 1006)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
( 1007)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
( 1008)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
( 1009)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
( 1010)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
( 1011) 
( 1012)        ! custom variable
( 1013)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
( 1014)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
( 1015)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
( 1016)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
( 1017)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
( 1018)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', varname)
( 1019)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
( 1020)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
( 1021)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
( 1022)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')
( 1023) 
( 1024)        ! global
( 1025)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
( 1026)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
( 1027)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
( 1028)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
( 1029)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
( 1030) 
( 1031)        varnum = 5
( 1032) !AD/DB: add new time-related variables






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 19

( 1033)        varnum = varnum + 1
( 1034)        ! ndate (ndastp)
( 1035)        nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
( 1036)             (/ dimids(1) /), &
( 1037)             varids(varnum))
( 1038)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
( 1039)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'time step date in year/month/day aammjj')
( 1040) 
( 1041)        varnum = varnum + 1
( 1042)        ! ndate (model_time)
( 1043)        nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
( 1044)             (/ dimids(1) /), &
( 1045)             varids(varnum))
( 1046)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', &
( 1047)             'time step date (when output is writen) in year/month/day aammjj (decimal day)')
( 1048)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
( 1049)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula1', 'nyear  =   model_time / 10000')       
( 1050)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula2', & 
( 1051)             'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
( 1052)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula3', & 
( 1053)             'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           
( 1054) 
( 1055)        varnum = varnum + 1
( 1056)        ! kt 
( 1057)        nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
( 1058)             (/ dimids(1) /), &
( 1059)             varids(varnum))
( 1060) 
( 1061)        
( 1062)        ! Close file
( 1063)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1064)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Closing file:', filename
( 1065)           CALL FLUSH
( 1066)        END IF
( 1067)        nfstat = nf90_close(ncid)
( 1068)        IF(nfstat /= nf90_noerr) THEN
( 1069)           status = NCDF_NFERR
( 1070)           RETURN
( 1071)        END IF
( 1072)     END IF
( 1073) 
( 1074)     ! Sync up processors before returning from subroutine
( 1075)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1076)     IF(mpistat /= 0) THEN
( 1077)        status = NCDF_MPERR
( 1078)        RETURN
( 1079)     END IF
( 1080)     
( 1081)   END SUBROUTINE ncdf_create_file_ave
( 1082)   
( 1083)   ! ncdf_create_restart builds a standard OPA restart file with all the default
( 1084)   ! dimensions, variables and attributes. Note that this is a single restart
( 1085)   ! file, and therefore incompatible with the default IOIPSL-based restart
( 1086)   ! routines (though it holds the same data)
( 1087)   ! NOTE: Some keys we haven't been using are untested, so I don't know if those
( 1088)   ! fields will be created correctly.
( 1089)   SUBROUTINE ncdf_create_restart(filename, status)
( 1090)     IMPLICIT NONE






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 20

( 1091)     ! Subroutine argument declarations
( 1092)     CHARACTER(LEN=*),INTENT(IN) :: filename
( 1093)     INTEGER,INTENT(OUT) :: status
( 1094)     
( 1095)     ! Local declarations
( 1096)     INTEGER :: ncid,    &  ! netCDF file ID
( 1097)          varid,   &  ! ID of netCDF variable to be written to
( 1098)          nfstat,  &  ! netCDF library call return status
( 1099)          mpistat     ! MPI library call return status
( 1100)     INTEGER,DIMENSION(1:8) :: dimids
( 1101)     INTEGER,DIMENSION(1:36) :: varids
( 1102)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
( 1103)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
( 1104)     CHARACTER(LEN=100) :: sec_since
( 1105)     CHARACTER(LEN=100) :: tstp_since
( 1106)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
( 1107)     CHARACTER(LEN=3),PARAMETER :: &
( 1108)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
( 1109)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
( 1110)     
( 1111)     ! Initializations
( 1112)     status = NCDF_NOERR
( 1113)     
( 1114)     ! Initializations
( 1115)     status = NCDF_NOERR
( 1116)     CALL ioget_calendar(cal_type)
( 1117)     CALL ioget_timestamp(timestamp)
( 1118)     WRITE (UNIT=sec_since, &
( 1119)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 1120)          &  nyear,nmonth,nday,0, 0, 0
( 1121)     WRITE (UNIT=tstp_since, &
( 1122)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 1123)          &  nyear,nmonth,nday,0, 0, 0
( 1124)     WRITE(t_origin, &
( 1125)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
( 1126)          &   nyear,months(nmonth),nday,0,0,0
( 1127) 
( 1128)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1129)        WRITE(100,*) 'NCDF DEBUG: Creating default restart file:', filename
( 1130)        CALL FLUSH
( 1131)     END IF
( 1132)     
( 1133)     ! Only processor 0 does anything
( 1134)     IF(nproc == 0) THEN
( 1135)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1136)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Creating file:', filename
( 1137)           CALL FLUSH
( 1138)        END IF
( 1139)        ! Create the file
( 1140)        nfstat = nf90_create(filename, nf90_clobber, ncid)
( 1141)        
( 1142)        ! Define dimensions
( 1143)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1144)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining dimensions in file:', filename
( 1145)           CALL FLUSH
( 1146)        END IF
( 1147)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(1))
( 1148)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(2))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 21

( 1149)        nfstat = nf90_def_dim(ncid, 'z', jpkdta, dimids(3))
( 1150)        nfstat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimids(4))
( 1151)        nfstat = nf90_def_dim(ncid, 'x_a', 1, dimids(5))
( 1152)        nfstat = nf90_def_dim(ncid, 'y_a', 1, dimids(6))
( 1153)        nfstat = nf90_def_dim(ncid, 'z_a', 10, dimids(7))
( 1154)        nfstat = nf90_def_dim(ncid, 'z_b', 1, dimids(8))
( 1155)        
( 1156)        ! Define variables
( 1157)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1158)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining variables in file:', filename
( 1159)           CALL FLUSH
( 1160)        END IF
( 1161)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
( 1162)             (/ dimids(1), dimids(2) /), &
( 1163)             varids(1))
( 1164)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
( 1165)             (/ dimids(1), dimids(2) /), &
( 1166)             varids(2))
( 1167)        nfstat = nf90_def_var(ncid, 'nav_lev', nf90_float, &
( 1168)             (/ dimids(3) /), &
( 1169)             varids(3))
( 1170)        nfstat = nf90_def_var(ncid, 'time', nf90_float, &
( 1171)             (/ dimids(4) /), &
( 1172)             varids(4))
( 1173)        nfstat = nf90_def_var(ncid, 'time_steps', nf90_int, &
( 1174)             (/ dimids(4) /), &
( 1175)             varids(5))
( 1176)        nfstat = nf90_def_var(ncid, 'info', nf90_double, &
( 1177)             (/ dimids(7) /), &
( 1178)             varids(6))
( 1179)        nfstat = nf90_def_var(ncid, 'ub', nf90_double, &
( 1180)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1181)             varids(7))
( 1182)        nfstat = nf90_def_var(ncid, 'vb', nf90_double, &
( 1183)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1184)             varids(8))
( 1185)        nfstat = nf90_def_var(ncid, 'tb', nf90_double, &
( 1186)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1187)             varids(9))
( 1188)        nfstat = nf90_def_var(ncid, 'sb', nf90_double, &
( 1189)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1190)             varids(10))
( 1191)        nfstat = nf90_def_var(ncid, 'rotb', nf90_double, &
( 1192)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1193)             varids(11))
( 1194)        nfstat = nf90_def_var(ncid, 'hdivb', nf90_double, &
( 1195)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1196)             varids(12))
( 1197)        nfstat = nf90_def_var(ncid, 'un', nf90_double, &
( 1198)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1199)             varids(13))
( 1200)        nfstat = nf90_def_var(ncid, 'vn', nf90_double, &
( 1201)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1202)             varids(14))
( 1203)        nfstat = nf90_def_var(ncid, 'tn', nf90_double, &
( 1204)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1205)             varids(15))
( 1206)        nfstat = nf90_def_var(ncid, 'sn', nf90_double, &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 22

( 1207)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1208)             varids(16))
( 1209)        nfstat = nf90_def_var(ncid, 'rotn', nf90_double, &
( 1210)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1211)             varids(17))
( 1212)        nfstat = nf90_def_var(ncid, 'hdivn', nf90_double, &
( 1213)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1214)             varids(18))
( 1215)        nfstat = nf90_def_var(ncid, 'gcx', nf90_double, &
( 1216)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1217)             varids(19))
( 1218)        nfstat = nf90_def_var(ncid, 'gcxb', nf90_double, &
( 1219)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1220)             varids(20))
( 1221) # if defined key_dynspg_rl
( 1222)        nfstat = nf90_def_var(ncid, 'bsfb', nf90_double, &
( 1223)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1224)             varids(21))
( 1225)        nfstat = nf90_def_var(ncid, 'bsfn', nf90_double, &
( 1226)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1227)             varids(22))
( 1228)        nfstat = nf90_def_var(ncid, 'bsfd', nf90_double, &
( 1229)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1230)             varids(27))
( 1231) # else 
( 1232)        nfstat = nf90_def_var(ncid, 'sshb', nf90_double, &
( 1233)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1234)             varids(21))
( 1235)        nfstat = nf90_def_var(ncid, 'sshn', nf90_double, &
( 1236)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1237)             varids(22))
( 1238) #  if defined key_dynspg_ts
( 1239)        nfstat = nf90_def_var(ncid, 'sshb_b', nf90_double, &
( 1240)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1241)             varids(23))
( 1242)        nfstat = nf90_def_var(ncid, 'sshn_b', nf90_double, &
( 1243)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1244)             varids(24))
( 1245)        nfstat = nf90_def_var(ncid, 'un_b', nf90_double, &
( 1246)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1247)             varids(25))
( 1248)        nfstat = nf90_def_var(ncid, 'vn_b', nf90_double, &
( 1249)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1250)             varids(26))
( 1251) #endif
( 1252) #endif
( 1253) 
( 1254)        ! Fields that are only defined if specific keys are set
( 1255) # if defined key_zdftke   ||   defined key_esopa
( 1256)        nfstat = nf90_def_var(ncid, 'en', nf90_double, &
( 1257)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1258)             varids(28))
( 1259) # endif
( 1260) # if defined key_ice_lim
( 1261)        nfstat = nf90_def_var(ncid, 'nfice', nf90_double, &
( 1262)             (/dimids(8) /), &
( 1263)             varids(29))
( 1264)        nfstat = nf90_def_var(ncid, 'sst_io', nf90_double, &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 23

( 1265)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1266)             varids(30))
( 1267)        nfstat = nf90_def_var(ncid, 'sss_io', nf90_double, &
( 1268)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1269)             varids(31))
( 1270)        nfstat = nf90_def_var(ncid, 'u_io', nf90_double, &
( 1271)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1272)             varids(32))
( 1273)        nfstat = nf90_def_var(ncid, 'v_io', nf90_double, &
( 1274)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1275)             varids(33))
( 1276) # if defined key_coupled
( 1277)        nfstat = nf90_def_var(ncid, 'alb_ice', nf90_double, &
( 1278)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1279)             varids(34))
( 1280) # endif
( 1281) # endif
( 1282) # if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
( 1283)        nfstat = nf90_def_var(ncid, 'nfbulk', nf90_double, &
( 1284)             (/ dimids(8) /), &
( 1285)             varids(35))
( 1286) 
( 1287)        nfstat = nf90_def_var(ncid, 'gsst', nf90_double, &
( 1288)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1289)             varids(36))
( 1290) # endif
( 1291)        
( 1292)        ! Add attributes
( 1293)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1294)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Writing attributes in file:', filename
( 1295)           CALL FLUSH
( 1296)        END IF
( 1297)        ! nav_lon
( 1298)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
( 1299)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -1.800000E2)
( 1300)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', 1.800000E2)
( 1301)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
( 1302) 
( 1303)        ! nav_lat
( 1304)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
( 1305)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', -9.000000E1)
( 1306)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 9.000000E1)
( 1307)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
( 1308) 
( 1309)        ! nav_lev
( 1310)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'model_levels')
( 1311)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E0)
( 1312)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E3)
( 1313)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Model levels')
( 1314) 
( 1315)        ! time
( 1316)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
( 1317)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
( 1318)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
( 1319)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'Time axis')
( 1320)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', '0001-JUL-01 00:00:00')
( 1321) 
( 1322)        ! time_steps






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 24

( 1323)        nfstat = nf90_put_att(ncid, varids(5), 'units', TRIM(tstp_since))
( 1324)        nfstat = nf90_put_att(ncid, varids(5), 'title', 'Time steps')
( 1325)        nfstat = nf90_put_att(ncid, varids(5), 'tstep_sec', 1.877760E7)
( 1326)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Time step axis')
( 1327)        nfstat = nf90_put_att(ncid, varids(5), 'time_origin', TRIM(t_origin))
( 1328) 
( 1329)        ! info
( 1330)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
( 1331) 
( 1332)        ! ub
( 1333)        nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
( 1334) 
( 1335)        ! vb
( 1336)        nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
( 1337) 
( 1338)        ! tb
( 1339)        nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)
( 1340) 
( 1341)        ! sb
( 1342)        nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)
( 1343) 
( 1344)        ! rotb
( 1345)        nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)
( 1346) 
( 1347)        ! hdivb
( 1348)        nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)
( 1349) 
( 1350)        ! un
( 1351)        nfstat = nf90_put_att(ncid, varids(13), 'missing_value', 1.000000E20)
( 1352) 
( 1353)        ! vn
( 1354)        nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)
( 1355) 
( 1356)        ! tn
( 1357)        nfstat = nf90_put_att(ncid, varids(15), 'missing_value', 1.000000E20)
( 1358) 
( 1359)        ! sn
( 1360)        nfstat = nf90_put_att(ncid, varids(16), 'missing_value', 1.000000E20)
( 1361) 
( 1362)        ! rotn
( 1363)        nfstat = nf90_put_att(ncid, varids(17), 'missing_value', 1.000000E20)
( 1364) 
( 1365)        ! hdivn
( 1366)        nfstat = nf90_put_att(ncid, varids(18), 'missing_value', 1.000000E20)
( 1367) 
( 1368)        ! gcx
( 1369)        nfstat = nf90_put_att(ncid, varids(19), 'missing_value', 1.000000E20)
( 1370) 
( 1371)        ! gcxb
( 1372)        nfstat = nf90_put_att(ncid, varids(20), 'missing_value', 1.000000E20)
( 1373) # if defined key_dynspg_rl
( 1374)        ! bsfb
( 1375)        nfstat = nf90_put_att(ncid, varids(21), 'missing_value', 1.000000E20)
( 1376) 
( 1377)        ! bsfn
( 1378)        nfstat = nf90_put_att(ncid, varids(22), 'missing_value', 1.000000E20)
( 1379) 
( 1380)        !bsfd






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 25

( 1381)        nfstat = nf90_put_att(ncid, varids(27), 'missing_value', 1.000000E20)
( 1382) # else
( 1383)        ! sshb
( 1384)        nfstat = nf90_put_att(ncid, varids(21), 'missing_value', 1.000000E20)
( 1385) 
( 1386)        ! sshn
( 1387)        nfstat = nf90_put_att(ncid, varids(22), 'missing_value', 1.000000E20)
( 1388) 
( 1389) #  if defined key_dynspg_ts
( 1390)        ! sshb_b
( 1391)        nfstat = nf90_put_att(ncid, varids(23), 'missing_value', 1.000000E20)
( 1392) 
( 1393)        ! sshn_b
( 1394)        nfstat = nf90_put_att(ncid, varids(24), 'missing_value', 1.000000E20)
( 1395) 
( 1396)        ! un_b
( 1397)        nfstat = nf90_put_att(ncid, varids(25), 'missing_value', 1.000000E20)
( 1398) 
( 1399)        ! vn_b
( 1400)        nfstat = nf90_put_att(ncid, varids(26), 'missing_value', 1.000000E20)
( 1401) # endif
( 1402) # endif
( 1403) # if defined key_zdftke   ||   defined key_esopa
( 1404)        ! en
( 1405)        nfstat = nf90_put_att(ncid, varids(28), 'missing_value', 1.000000E20)
( 1406) # endif
( 1407) # if defined key_ice_lim
( 1408)        ! nfice
( 1409)        nfstat = nf90_put_att(ncid, varids(29), 'missing_value', 1.000000E20)
( 1410) 
( 1411)        ! sst_io
( 1412)        nfstat = nf90_put_att(ncid, varids(30), 'missing_value', 1.000000E20)
( 1413) 
( 1414)        ! sss_io
( 1415)        nfstat = nf90_put_att(ncid, varids(31), 'missing_value', 1.000000E20)
( 1416) 
( 1417)        ! u_io
( 1418)        nfstat = nf90_put_att(ncid, varids(32), 'missing_value', 1.000000E20)
( 1419) 
( 1420)        ! v_io
( 1421)        nfstat = nf90_put_att(ncid, varids(33), 'missing_value', 1.000000E20)
( 1422) 
( 1423) # if defined key_coupled
( 1424)        ! alb_ice
( 1425)        nfstat = nf90_put_att(ncid, varids(34), 'missing_value', 1.000000E20)
( 1426) 
( 1427) # endif
( 1428) # endif
( 1429) # if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
( 1430)        ! nfbulk
( 1431)        nfstat = nf90_put_att(ncid, varids(35), 'missing_value', 1.000000E20)
( 1432)  
( 1433)        ! gsst
( 1434)        nfstat = nf90_put_att(ncid, varids(36), 'missing_value', 1.000000E20)
( 1435)  
( 1436) # endif
( 1437) 
( 1438)        ! global






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 26

( 1439)        nfstat = nf90_put_att(ncid, NF90_GLOBAL , 'Conventions', 'GDT 1.2')
( 1440)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
( 1441)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
( 1442)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number_total', 0)
( 1443)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number', 0)
( 1444)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids', 0)
( 1445)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_global', 0)
( 1446)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_local', 0)
( 1447)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_first', 0)
( 1448)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_last', 0)
( 1449)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', 0)
( 1450)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', 0)
( 1451)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_type','box')
( 1452) 
( 1453)        ! Close file
( 1454)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1455)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Closing file:', filename
( 1456)           CALL FLUSH
( 1457)        END IF
( 1458)        nfstat = nf90_close(ncid)
( 1459)        IF(nfstat /= nf90_noerr) THEN
( 1460)           status = NCDF_NFERR
( 1461)           RETURN
( 1462)        END IF
( 1463)     END IF
( 1464) 
( 1465)     ! Sync up processors before returning from subroutine
( 1466)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1467)     IF(mpistat /= 0) THEN
( 1468)        status = NCDF_MPERR
( 1469)        RETURN
( 1470)     END IF
( 1471)        
( 1472)   END SUBROUTINE ncdf_create_restart
( 1473) 
( 1474)   ! ncdf_writesv (single value) writes a single, scalar value to a specified
( 1475)   ! index in a netCDF variable. The netCDF variable is assumed to be
( 1476)   ! 1-dimensional.
( 1477)   ! filename - file to write to
( 1478)   ! varname - variable to write to
( 1479)   ! data - scalar value to write
( 1480)   ! index - where in the variable to put the data
( 1481)   ! status - return status of the subroutine
( 1482)   SUBROUTINE ncdf_writesv(filename, varname, data, index, status)
( 1483)      IMPLICIT NONE
( 1484)     ! Subroutine argument declarations
( 1485)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 1486)     REAL(wp),INTENT(IN) :: data
( 1487)     INTEGER,INTENT(IN) :: index   ! Where in the variable to write the value
( 1488)     INTEGER,INTENT(OUT) :: status
( 1489) 
( 1490)     ! Local declarations
( 1491)     INTEGER :: ncid,    &  ! netCDF file ID
( 1492)                varid,   &  ! ID of netCDF variable to be written to
( 1493)                nfstat,  &  ! netCDF library call return status
( 1494)                mpistat     ! MPI library call return status
( 1495) 
( 1496)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 27

( 1497)        WRITE(100,*) 'NCDF DEBUG: Writing single value to file:', filename
( 1498)        CALL FLUSH
( 1499)     END IF
( 1500) 
( 1501)     ! Initializations
( 1502)     status = NCDF_NOERR
( 1503) 
( 1504)     ! Open netCDF file and get info
( 1505)     IF(nproc == 0) THEN
( 1506)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1507)           WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Opening file:', filename
( 1508)           CALL FLUSH
( 1509)        END IF
( 1510)        nfstat = nf90_open(filename, nf90_write, ncid)
( 1511)        IF(nfstat /= nf90_noerr) THEN
( 1512)           status = NCDF_NFERR
( 1513)           RETURN
( 1514)        END IF
( 1515)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1516)           WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Getting info from file:', filename
( 1517)           CALL FLUSH
( 1518)        END IF
( 1519)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 1520)        IF(nfstat /= nf90_noerr) THEN
( 1521)           status = NCDF_NFERR
( 1522)           RETURN
( 1523)        END IF
( 1524)        
( 1525)        ! Write data to netCDF file
( 1526)        ! This subroutine assumes all processors have the same data in this
( 1527)        ! variable, no attempt to sync data or merge data in the file is made
( 1528)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1529)           WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Writing data to file:', filename
( 1530)           CALL FLUSH
( 1531)        END IF
( 1532)        nfstat = nf90_put_var(ncid, varid, data, &
( 1533)            (/ index /))
( 1534)        IF(nfstat /= nf90_noerr) THEN
( 1535)           status = NCDF_NFERR
( 1536)           RETURN
( 1537)        END IF
( 1538)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1539)           WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Closing file:', filename
( 1540)           CALL FLUSH
( 1541)        END IF
( 1542)        nfstat = nf90_close(ncid)
( 1543)        IF(nfstat /= nf90_noerr) THEN
( 1544)           status = NCDF_NFERR
( 1545)           RETURN
( 1546)        END IF
( 1547)     END IF
( 1548) 
( 1549)     ! Sync up processors before returning from subroutine
( 1550)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1551)     IF(mpistat /= 0) THEN
( 1552)        status = NCDF_MPERR
( 1553)        RETURN
( 1554)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 28

( 1555)     
( 1556)   END SUBROUTINE ncdf_writesv
( 1557) 
( 1558)   ! ncdf_write1d writes a 1-dimensional array to a 1-D variable in a netCDF
( 1559)   ! file. It assumes that all processors have the same data in the array, 
( 1560)   ! so no attempt is made to sync data between processors or merge data
( 1561)   ! together in the file, unlike ncdf_write2d and ncdf_write3d
( 1562)   ! filename - file to write to
( 1563)   ! varname - variable to write to
( 1564)   ! data - 1-D array to write
( 1565)   ! status - return status of the subroutine
( 1566)   SUBROUTINE ncdf_write1d(filename, varname, data, status)
( 1567)     IMPLICIT NONE
( 1568)     ! Subroutine argument declarations
( 1569)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 1570)     REAL(wp),DIMENSION(:),INTENT(IN) :: data
( 1571)     INTEGER,INTENT(OUT) :: status
( 1572) 
( 1573)     ! Local declarations
( 1574)     INTEGER :: ncid,    &  ! netCDF file ID
( 1575)                varid,   &  ! ID of netCDF variable to be written to
( 1576)                nfstat,  &  ! netCDF library call return status
( 1577)                mpistat     ! MPI library call return status
( 1578) 
( 1579)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1580)        WRITE(100,*) 'NCDF DEBUG: Writing 1D array to file:', filename
( 1581)        CALL FLUSH
( 1582)     END IF
( 1583) 
( 1584)     ! Initializations
( 1585)     status = NCDF_NOERR
( 1586) 
( 1587)     ! Open netCDF file and get info
( 1588)     IF(nproc == 0) THEN
( 1589)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1590)           WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Opening file:', filename
( 1591)           CALL FLUSH
( 1592)        END IF
( 1593)        nfstat = nf90_open(filename, nf90_write, ncid)
( 1594)        IF(nfstat /= nf90_noerr) THEN
( 1595)           status = NCDF_NFERR
( 1596)           RETURN
( 1597)        END IF
( 1598)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1599)           WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Getting info from file:', filename
( 1600)           CALL FLUSH
( 1601)        END IF
( 1602)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 1603)        IF(nfstat /= nf90_noerr) THEN
( 1604)           status = NCDF_NFERR
( 1605)           RETURN
( 1606)        END IF
( 1607)        
( 1608)        ! Write data to netCDF file
( 1609)        ! This subroutine assumes all processors have the same data in this
( 1610)        ! array, no attempt to sync data or merge data in the file is made
( 1611)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1612)           WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Writing data to file:', filename






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 29

( 1613)           CALL FLUSH
( 1614)        END IF
( 1615)        nfstat = nf90_put_var(ncid, varid, data)
( 1616)        IF(nfstat /= nf90_noerr) THEN
( 1617)           status = NCDF_NFERR
( 1618)           WRITE(125,*) 'ERROR: ', nf90_strerror(nfstat)
( 1619)           CALL FLUSH
( 1620)           RETURN
( 1621)        END IF
( 1622)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1623)           WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Closing file:', filename
( 1624)           CALL FLUSH
( 1625)        END IF
( 1626)        nfstat = nf90_close(ncid)
( 1627)        IF(nfstat /= nf90_noerr) THEN
( 1628)           status = NCDF_NFERR
( 1629)           RETURN
( 1630)        END IF
( 1631)     END IF
( 1632) 
( 1633)     ! Sync up processors before returning from subroutine
( 1634)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1635)     IF(mpistat /= 0) THEN
( 1636)        status = NCDF_MPERR
( 1637)        RETURN
( 1638)     END IF
( 1639) 
( 1640)   END SUBROUTINE ncdf_write1d
( 1641) 
( 1642)   ! ncdf_write2d writes a 2-dimensional array to a 2-D variable in a netCDF
( 1643)   ! file. It assumes the data array being passed is a subdomain - that is,
( 1644)   ! each processor's array contains different values, and the data needs to
( 1645)   ! be written to the location in the netCDF file which corresponds to the
( 1646)   ! processor writing the subdomain. This subroutine will correctly merge
( 1647)   ! all subroutines so that the variable in the netCDF file contains
( 1648)   ! the entire global domain of values
( 1649)   ! filename - file to write to
( 1650)   ! varname - variable to write to
( 1651)   ! data - 2-D array to write
( 1652)   ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
( 1653)   !                                      will force an immediate write, with
( 1654)   !                                      data being written at a time index in
( 1655)   !                                      the file equal to (tstep * -1). Use this
( 1656)   !                                      to output on timesteps which are not
( 1657)   !                                      multiples of nwrite
( 1658)   ! status - return status of the subroutine
( 1659)   ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
( 1660)   ! a tstep which is <0
( 1661)   SUBROUTINE ncdf_write2d(filename, varname, data, tstep, status)
( 1662)     IMPLICIT NONE
( 1663)     ! Subroutine argument declarations
( 1664)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 1665)     INTEGER,INTENT(IN) :: tstep
( 1666)     REAL(wp),DIMENSION(:,:),INTENT(IN) :: data
( 1667)     INTEGER,INTENT(OUT) :: status
( 1668) 
( 1669)     ! Local declarations
( 1670)     INTEGER :: ncid,    &  ! netCDF file ID






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 30

( 1671)                varid,   &  ! ID of netCDF variable to be written to
( 1672)                tindex,  &  ! Index to write to on the time axis
( 1673)                mpistat, &  ! MPI library call return status
( 1674)                nfstat,  &  ! netCDF library call return status
( 1675)                is,      &  ! i-axis start index for writing to file
( 1676)                js,      &  ! j-axis   "     "    "     "     "  "
( 1677)                i,k,       &  ! Loop counter
( 1678)                ndims,   &  ! Number of dimensions in this variabl
( 1679)                iz,      &  ! I-size of array to write to file
( 1680)                jz          ! J-size of array to write to file
( 1681)     
( 1682)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 1683)     INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
( 1684)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 1685)     REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: buf   ! Send/receive buffer for data array
( 1686)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 1687) 
( 1688)     ! Exit subroutine if no write should take place on current timestep
( 1689)     IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >= 0)) THEN
( 1690)        status = NCDF_NOERR
( 1691)        RETURN
( 1692)     END IF
( 1693) 
( 1694)     ! Initializations
( 1695)     is = nimpp                          ! Get global i-index for this subdomain
( 1696)     js = njmpp                          ! Get global j-index  "   "      "
( 1697)     iz = nlei                           ! Get i-size of subdomain section to be written
( 1698)     jz = nlej                           ! Get j-size of subdomain section to be written
( 1699)     ALLOCATE(buf(1:iz,1:jz))            ! Allocate send/receive buffer
( 1700)     buf = data(1:iz,1:jz)               ! Copy section of subdomain to buffer
( 1701)     IF(tstep >= 0) THEN
( 1702)        tindex = tstep/nwrite - nit000/nwrite
( 1703)     ELSE
( 1704)        tindex = tstep * (-1)
( 1705)     END IF
( 1706)     status = NCDF_NOERR
( 1707)     hastime = .FALSE.
( 1708) 
( 1709)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1710)        WRITE(100,*) 'NCDF DEBUG: Writing 2D array to file:', filename
( 1711)        CALL FLUSH
( 1712)     END IF
( 1713) 
( 1714)     ! Open netCDF file and get info
( 1715)     IF(nproc == 0) THEN
( 1716)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1717)           WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Opening file:', filename
( 1718)           CALL FLUSH
( 1719)        END IF
( 1720)        nfstat = nf90_open(filename, nf90_write, ncid)
( 1721)        IF(nfstat /= nf90_noerr) THEN
( 1722)           status = NCDF_NFERR
( 1723)           RETURN
( 1724)        END IF
( 1725)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1726)           WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Getting info from file:', filename
( 1727)           CALL FLUSH
( 1728)        END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 31

( 1729)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 1730)        IF(nfstat /= nf90_noerr) THEN
( 1731)           status = NCDF_NFERR
( 1732)           RETURN
( 1733)        END IF
( 1734)        
( 1735)        ! Determine if this variable contains a time axis
( 1736)        nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 1737)             ndims=ndims, dimids=var_dimids)
( 1738)        IF(nfstat /= nf90_noerr) THEN
( 1739)           status = NCDF_NFERR
( 1740)           RETURN
( 1741)        END IF
( 1742)        DO i = 1, ndims
( 1743)           nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 1744)                name=dname, len=k)
( 1745) !          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 1746) !               name=dname)
( 1747)           dname=TRIM(dname)
( 1748)           IF(nfstat /= nf90_noerr) THEN
( 1749)              status = NCDF_NFERR
( 1750)              RETURN
( 1751)           END IF
( 1752)           IF(dname=='time_counter') THEN
( 1753)              hastime = .TRUE.
( 1754)           END IF
( 1755)        END DO
( 1756)     END IF
( 1757) 
( 1758)     ! Get data from each processor in turn, processor 0 writes to file
( 1759)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1760)        WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Writing data to file:', filename
( 1761)        CALL FLUSH
( 1762)     END IF
( 1763) 
( 1764)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1765)     IF(mpistat /= 0) THEN
( 1766)        status = NCDF_MPERR
( 1767)        RETURN
( 1768)     END IF
( 1769)     DO i = 0, (jpnij - 1)
( 1770)        ! Since processor 0 doesn't need to send data to itself, only
( 1771)        ! processors with IDs > 0 make send calls
( 1772)        IF(i > 0) THEN
( 1773)           ! If it's the local processor's turn to send, go ahead and
( 1774)           ! send data and write indices to processor 0
( 1775)           IF(nproc == i) THEN
( 1776)              CALL MPI_SEND(is, 1, MPI_INTEGER, &
( 1777)                   0, i, MPI_COMM_WORLD, mpistat)
( 1778)              IF(mpistat /= 0) THEN
( 1779)                 status = NCDF_MPERR
( 1780)                 RETURN
( 1781)              END IF
( 1782)              CALL MPI_SEND(js, 1, MPI_INTEGER, &
( 1783)                   0, i+1, MPI_COMM_WORLD, mpistat)
( 1784)              IF(mpistat /= 0) THEN
( 1785)                 status = NCDF_MPERR
( 1786)                 RETURN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 32

( 1787)              END IF
( 1788)              CALL MPI_SEND(iz, 1, MPI_INTEGER, &
( 1789)                   0, i+2, MPI_COMM_WORLD, mpistat)
( 1790)              IF(mpistat /= 0) THEN
( 1791)                 status = NCDF_MPERR
( 1792)                 RETURN
( 1793)              END IF
( 1794)              CALL MPI_SEND(jz, 1, MPI_INTEGER, &
( 1795)                   0, i+3, MPI_COMM_WORLD, mpistat)
( 1796)              IF(mpistat /= 0) THEN
( 1797)                 status = NCDF_MPERR
( 1798)                 RETURN
( 1799)              END IF
( 1800)              CALL MPI_SEND(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 1801)                   0, i+4, MPI_COMM_WORLD, mpistat)
( 1802)              IF(mpistat /= 0) THEN
( 1803)                 status = NCDF_MPERR
( 1804)                 RETURN
( 1805)              END IF
( 1806)           END IF
( 1807)        END IF
( 1808) 
( 1809)        ! Processor 0 receives data and writes it to file
( 1810)        IF(nproc == 0) THEN
( 1811)           ! Since processor 0 doesn't need to receive data from itself, 
( 1812)           ! we only make receive calls if it's another processor's turn
( 1813)           ! to send data
( 1814)           IF(i > 0) THEN
( 1815)              CALL MPI_RECV(is, 1, MPI_INTEGER, &
( 1816)                   i, i, MPI_COMM_WORLD, rcvstat, mpistat)
( 1817)              IF(mpistat /= 0) THEN
( 1818)                 status = NCDF_MPERR
( 1819)                 RETURN
( 1820)              END IF
( 1821)              CALL MPI_RECV(js, 1, MPI_INTEGER, &
( 1822)                   i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)
( 1823)              IF(mpistat /= 0) THEN
( 1824)                 status = NCDF_MPERR
( 1825)                 RETURN
( 1826)              END IF
( 1827)              CALL MPI_RECV(iz, 1, MPI_INTEGER, &
( 1828)                   i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
( 1829)              IF(mpistat /= 0) THEN
( 1830)                 status = NCDF_MPERR
( 1831)                 RETURN
( 1832)              END IF
( 1833)              CALL MPI_RECV(jz, 1, MPI_INTEGER, &
( 1834)                   i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
( 1835)              IF(mpistat /= 0) THEN
( 1836)                 status = NCDF_MPERR
( 1837)                 RETURN
( 1838)              END IF
( 1839)              DEALLOCATE(buf)
( 1840)              ALLOCATE(buf(1:iz,1:jz))
( 1841)              CALL MPI_RECV(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 1842)                   i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
( 1843)              IF(mpistat /= 0) THEN
( 1844)                 status = NCDF_MPERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 33

( 1845)                 RETURN
( 1846)              END IF
( 1847)           END IF
( 1848)           ! Got the data, write it to the proper location in the netCDF file
( 1849)           IF(hastime .EQV. .TRUE.) THEN
( 1850)              nfstat = nf90_put_var(ncid, varid, buf, &
( 1851)                   (/ is, js, tindex/), &
( 1852)                   (/ iz, jz, 1 /))
( 1853)           ELSE
( 1854)              nfstat = nf90_put_var(ncid, varid, buf, &
( 1855)                   (/ is, js /), &
( 1856)                   (/ iz, jz /))
( 1857)           END IF
( 1858)           IF(nfstat /= nf90_noerr) THEN
( 1859)              status = NCDF_NFERR
( 1860)              RETURN
( 1861)           END IF
( 1862)        END IF
( 1863) 
( 1864)        ! Sync all processors at the end of each loop iteration to prevent
( 1865)        ! concurrency issues/race conditions/various MPI-related ugliness
( 1866)        CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1867)        IF(mpistat /= 0) THEN
( 1868)           status = NCDF_MPERR
( 1869)           RETURN
( 1870)        END IF
( 1871)     END DO
( 1872) 
( 1873)     ! All done, close up the netCDF file
( 1874)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1875)        WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Closing file:', filename
( 1876)        CALL FLUSH
( 1877)     END IF
( 1878)     IF(nproc == 0) THEN
( 1879)        nfstat = nf90_close(ncid)
( 1880)        IF(nfstat /= nf90_noerr) THEN
( 1881)           status = NCDF_NFERR
( 1882)           RETURN
( 1883)        END IF
( 1884)     END IF
( 1885)     DEALLOCATE(buf)
( 1886) 
( 1887)     ! Sync up processors before returning from subroutine
( 1888)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1889)     IF(mpistat /= 0) THEN
( 1890)        status = NCDF_MPERR
( 1891)        RETURN
( 1892)     END IF
( 1893) 
( 1894)   END SUBROUTINE ncdf_write2d
( 1895) 
( 1896)   ! ncdf_write3d writes a 3-dimensional array to a 3-D variable in a netCDF
( 1897)   ! file. It assumes the data array being passed is a subdomain - that is,
( 1898)   ! each processor's array contains different values, and the data needs to
( 1899)   ! be written to the location in the netCDF file which corresponds to the
( 1900)   ! processor writing the subdomain. This subroutine will correctly merge
( 1901)   ! all subroutines so that the variable in the netCDF file contains
( 1902)   ! the entire global domain of values






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 34

( 1903)   ! filename - file to write to
( 1904)   ! varname - variable to write to
( 1905)   ! data - 3-D array to write
( 1906)   ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
( 1907)   !                                      will force an immediate write, with
( 1908)   !                                      data being written at a time index in
( 1909)   !                                      the file equal to (tstep * -1). Use this
( 1910)   !                                      to output on timesteps which are not
( 1911)   !                                      multiples of nwrite
( 1912)   ! status - return status of the subroutine
( 1913)   ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
( 1914)   ! a tstep which is <0
( 1915)   SUBROUTINE ncdf_write3d(filename, varname, data, tstep, status)
( 1916)     IMPLICIT NONE
( 1917)     ! Subroutine argument declarations
( 1918)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 1919)     INTEGER,INTENT(IN) :: tstep
( 1920)     REAL(wp),DIMENSION(:,:,:),INTENT(IN) :: data
( 1921)     INTEGER,INTENT(OUT) :: status
( 1922) 
( 1923)     ! Local declarations
( 1924)     INTEGER :: ncid,    &  ! netCDF file ID
( 1925)                varid,   &  ! ID of netCDF variable to be written to
( 1926)                tindex,  &  ! Index to write to on the time axis
( 1927)                mpistat, &  ! MPI library call return status
( 1928)                nfstat,  &  ! netCDF library call return status
( 1929)                !rcvstat, &  ! MPI receive status
( 1930)                is,      &  ! i-axis start index for writing to file
( 1931)                js,      &  ! j-axis   "     "    "     "     "  "
( 1932)                ks,      &  ! k-axis   "     "    "     "     "  "
( 1933)                i,k,     &  ! Loop counter
( 1934)                ndims,   &  ! Number of dimensions in this variable
( 1935)                iz,      &  ! I-size of array to write to file
( 1936)                jz          ! J-size of array to write to file
( 1937) 
( 1938)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 1939)     INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
( 1940)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: buf   ! Send/receive buffer for data array
( 1941)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 1942)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 1943) 
( 1944)     ! Exit subroutine if no write should take place on current timestep
( 1945)     IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >=0)) THEN
( 1946)        status = NCDF_NOERR
( 1947)        RETURN
( 1948)     END IF
( 1949) 
( 1950)     ! Initializations
( 1951)     is = nimpp                          ! Get global i-index for this subdomain
( 1952)     js = njmpp                          ! Get global j-index  "   "      "
( 1953)     ks = 1                              ! Depth writes always start at k=1
( 1954)     iz = nlei                           ! Get i-size of subdomain section to be written
( 1955)     jz = nlej                           ! Get j-size of subdomain section to be written
( 1956)     ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,DIM=3)))      ! Allocate send/recieve buffer
( 1957)     buf = data(1:iz,1:jz,1:SIZE(data,DIM=3))         ! Copy section of subdomain to buffer
( 1958)     IF(tstep >= 0) THEN
( 1959)        tindex = tstep/nwrite - nit000/nwrite
( 1960)     ELSE






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 35

( 1961)        tindex = tstep * (-1)
( 1962)     END IF
( 1963)     status = NCDF_NOERR
( 1964)     hastime = .FALSE.
( 1965) 
( 1966)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1967)        WRITE(100,*) 'NCDF DEBUG: Writing 3D array to file:', filename
( 1968)        CALL FLUSH
( 1969)     END IF
( 1970) 
( 1971)     ! Open netCDF file and get info
( 1972)     IF(nproc == 0) THEN
( 1973)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1974)           WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Opening file:', filename
( 1975)           CALL FLUSH
( 1976)        END IF
( 1977)        nfstat = nf90_open(filename, nf90_write, ncid)
( 1978)        IF(nfstat /= nf90_noerr) THEN
( 1979)           status = NCDF_NFERR
( 1980)           RETURN
( 1981)        END IF
( 1982)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1983)           WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Getting info from file:', filename
( 1984)           CALL FLUSH
( 1985)        END IF
( 1986)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 1987)        IF(nfstat /= nf90_noerr) THEN
( 1988)           status = NCDF_NFERR
( 1989)           RETURN
( 1990)        END IF
( 1991) 
( 1992)        ! Determine if this variable contains a time axis
( 1993)        nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 1994)             ndims=ndims, dimids=var_dimids)
( 1995)        IF(nfstat /= nf90_noerr) THEN
( 1996)           status = NCDF_NFERR
( 1997)           RETURN
( 1998)        END IF
( 1999)        DO i = 1, ndims
( 2000)           nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2001)                name=dname, len=k)
( 2002) !          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2003) !               name=dname)
( 2004)           dname=TRIM(dname)
( 2005)           IF(nfstat /= nf90_noerr) THEN
( 2006)              status = NCDF_NFERR
( 2007)              RETURN
( 2008)           END IF
( 2009)           IF(dname=='time_counter') THEN
( 2010)              hastime = .TRUE.
( 2011)              EXIT
( 2012)           END IF
( 2013)        END DO
( 2014)     END IF
( 2015) 
( 2016)     ! Get data from each processor in turn, processor 0 writes to file
( 2017)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2018)        WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Writing data to file:', filename






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 36

( 2019)        CALL FLUSH
( 2020)     END IF
( 2021)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2022)     IF(mpistat /= 0) THEN
( 2023)        status = NCDF_MPERR
( 2024)        RETURN
( 2025)     END IF
( 2026)     DO i=0, (jpnij - 1)
( 2027)        ! Since processor 0 doesn't need to send data to itself, only
( 2028)        ! processors with IDs > 0 make send calls
( 2029)        IF(i > 0) THEN
( 2030)           ! If it's the local processor's turn to send, go ahead and
( 2031)           ! send data and write indices to processor 0
( 2032)           IF(nproc == i) THEN
( 2033)              CALL MPI_SEND(is, 1, MPI_INTEGER, &
( 2034)                   0, i, MPI_COMM_WORLD, mpistat)
( 2035)              IF(mpistat /= 0) THEN
( 2036)                 status = NCDF_MPERR
( 2037)                 RETURN
( 2038)              END IF
( 2039)              CALL MPI_SEND(js, 1, MPI_INTEGER, &
( 2040)                   0, i+1, MPI_COMM_WORLD, mpistat)
( 2041)              IF(mpistat /= 0) THEN
( 2042)                 status = NCDF_MPERR
( 2043)                 RETURN
( 2044)              END IF
( 2045)              CALL MPI_SEND(ks, 1, MPI_INTEGER, &
( 2046)                   0, i+2, MPI_COMM_WORLD, mpistat)
( 2047)              IF(mpistat /= 0) THEN
( 2048)                 status = NCDF_MPERR
( 2049)                 RETURN
( 2050)              END IF
( 2051)              CALL MPI_SEND(iz, 1, MPI_INTEGER, &
( 2052)                   0, i+3, MPI_COMM_WORLD, mpistat)
( 2053)              IF(mpistat /= 0) THEN
( 2054)                 status = NCDF_MPERR
( 2055)                 RETURN
( 2056)              END IF
( 2057)              CALL MPI_SEND(jz, 1, MPI_INTEGER, &
( 2058)                   0, i+4, MPI_COMM_WORLD, mpistat)
( 2059)              IF(mpistat /= 0) THEN
( 2060)                 status = NCDF_MPERR
( 2061)                 RETURN
( 2062)              END IF
( 2063)              CALL MPI_SEND(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 2064)                   0, i+5, MPI_COMM_WORLD, mpistat)
( 2065)              IF(mpistat /= 0) THEN
( 2066)                 status = NCDF_MPERR
( 2067)                 RETURN
( 2068)              END IF
( 2069)           END IF
( 2070)        END IF
( 2071) 
( 2072)        ! Processor 0 receives data and writes it to file
( 2073)        IF(nproc == 0) THEN
( 2074)           ! Since processor 0 doesn't need to receive data from itself, 
( 2075)           ! we only make receive calls if it's another processor's turn
( 2076)           ! to send data






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 37

( 2077)           IF(i > 0) THEN
( 2078)              CALL MPI_RECV(is, 1, MPI_INTEGER, &
( 2079)                   i, i, MPI_COMM_WORLD, rcvstat, mpistat)
( 2080)              IF(mpistat /= 0) THEN
( 2081)                 status = NCDF_MPERR
( 2082)                 RETURN
( 2083)              END IF
( 2084)              CALL MPI_RECV(js, 1, MPI_INTEGER, &
( 2085)                   i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)
( 2086)              IF(mpistat /= 0) THEN
( 2087)                 status = NCDF_MPERR
( 2088)                 RETURN
( 2089)              END IF
( 2090)              CALL MPI_RECV(ks, 1, MPI_INTEGER, &
( 2091)                   i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
( 2092)              IF(mpistat /= 0) THEN
( 2093)                 status = NCDF_MPERR
( 2094)                 RETURN
( 2095)              END IF
( 2096)              CALL MPI_RECV(iz, 1, MPI_INTEGER, &
( 2097)                   i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
( 2098)              IF(mpistat /= 0) THEN
( 2099)                 status = NCDF_MPERR
( 2100)                 RETURN
( 2101)              END IF
( 2102)              CALL MPI_RECV(jz, 1, MPI_INTEGER, &
( 2103)                   i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
( 2104)              IF(mpistat /= 0) THEN
( 2105)                 status = NCDF_MPERR
( 2106)                 RETURN
( 2107)              END IF
( 2108)              DEALLOCATE(buf)
( 2109)              ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,DIM=3)))
( 2110)              CALL MPI_RECV(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 2111)                   i, i+5, MPI_COMM_WORLD, rcvstat, mpistat)
( 2112)              IF(mpistat /= 0) THEN
( 2113)                 status = NCDF_MPERR
( 2114)                 RETURN
( 2115)              END IF
( 2116)           END IF
( 2117)           
( 2118)           ! Got the data, write it to the proper location in the netCDF file
( 2119)           IF(hastime .EQV. .TRUE.) THEN
( 2120)              nfstat = nf90_put_var(ncid, varid, buf, &
( 2121)                   (/ is, js, ks, tindex/), &
( 2122)                   (/ iz, jz, SIZE(data,DIM=3), 1 /))
( 2123)           ELSE
( 2124)              nfstat = nf90_put_var(ncid, varid, buf, &
( 2125)                   (/ is, js, ks /), &
( 2126)                   (/ iz, jz, SIZE(data,DIM=3) /))
( 2127)           END IF
( 2128)           IF(nfstat /= nf90_noerr) THEN
( 2129)              status = NCDF_NFERR
( 2130)              RETURN
( 2131)           END IF
( 2132)        END IF
( 2133)        
( 2134)        ! Sync all processors at the end of each loop iteration to prevent






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 38

( 2135)        ! concurrency issues/race conditions/various MPI-related ugliness
( 2136)        CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2137)        IF(mpistat /= 0) THEN
( 2138)           status = NCDF_MPERR
( 2139)           RETURN
( 2140)        END IF
( 2141)     END DO
( 2142) 
( 2143)     ! All done, close up the netCDF file
( 2144)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2145)        WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Closing file:', filename
( 2146)        CALL FLUSH
( 2147)     END IF
( 2148)     IF(nproc == 0) THEN
( 2149)        nfstat = nf90_close(ncid)
( 2150)        IF(nfstat /= nf90_noerr) THEN
( 2151)           status = NCDF_NFERR
( 2152)           RETURN
( 2153)        END IF
( 2154)     END IF
( 2155)     DEALLOCATE(buf)
( 2156) 
( 2157)     ! Sync up processors before returning from subroutine
( 2158)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2159)     IF(mpistat /= 0) THEN
( 2160)        status = NCDF_MPERR
( 2161)        RETURN
( 2162)     END IF
( 2163)     
( 2164)   END SUBROUTINE ncdf_write3d
( 2165)   
( 2166)   ! ncdf_write4d writes a 4-dimensional array to a 4-D variable in a netCDF
( 2167)   ! file. It assumes the data array being passed is a subdomain - that is,
( 2168)   ! each processor's array contains different values, and the data needs to
( 2169)   ! be written to the location in the netCDF file which corresponds to the
( 2170)   ! processor writing the subdomain. This subroutine will correctly merge
( 2171)   ! all subroutines so that the variable in the netCDF file contains
( 2172)   ! the entire global domain of values
( 2173)   ! filename - file to write to
( 2174)   ! varname - variable to write to
( 2175)   ! data - 3-D array to write
( 2176)   ! dsz - Dimension SiZe - size of the 4th dimension in the variable
( 2177)   ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
( 2178)   !                                      will force an immediate write, with
( 2179)   !                                      data being written at a time index in
( 2180)   !                                      the file equal to (tstep * -1). Use this
( 2181)   !                                      to output on timesteps which are not
( 2182)   !                                      multiples of nwrite
( 2183)   ! status - return status of the subroutine
( 2184)   ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
( 2185)   ! a tstep which is < 0
( 2186)   ! NOT EXTENSIVELY TESTED - MAY STILL BE BUGGY
( 2187)   SUBROUTINE ncdf_write4d(filename, varname, data, dsz, tstep, status)
( 2188)     IMPLICIT NONE
( 2189)     ! Subroutine argument declarations
( 2190)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2191)     INTEGER,INTENT(IN) :: tstep, dsz
( 2192)     REAL,DIMENSION(:,:,:,:),INTENT(IN) :: data






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 39

( 2193)     INTEGER,INTENT(OUT) :: status
( 2194) 
( 2195)     ! Local declarations
( 2196)     INTEGER :: ncid,    &  ! netCDF file ID
( 2197)                varid,   &  ! ID of netCDF variable to be written to
( 2198)                tindex,  &  ! Index to write to on the time axis
( 2199)                mpistat, &  ! MPI library call return status
( 2200)                nfstat,  &  ! netCDF library call return status
( 2201)                is,      &  ! i-axis start index for writing to file
( 2202)                js,      &  ! j-axis   "     "    "     "     "  "
( 2203)                ks,      &  ! k-axis   "     "    "     "     "  "
( 2204)                i,k,       &  ! Loop counter
( 2205)                ndims,   &  ! Number of dimensions in this variable
( 2206)                iz,      &  ! I-size of array to write to file
( 2207)                jz          ! J-size of array to write to file
( 2208) 
( 2209)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 2210)     INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
( 2211)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:,:) :: buf   ! Send/receive buffer for data array
( 2212)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 2213)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 2214) 
( 2215)     ! Exit subroutine if no write should take place on current timestep
( 2216)     IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >= 0)) THEN
( 2217)        status = NCDF_NOERR
( 2218)        RETURN
( 2219)     END IF
( 2220) 
( 2221)     ! Initializations
( 2222)     is = nimpp                          ! Get global i-index for this subdomain
( 2223)     js = njmpp                          ! Get global j-index  "   "      "
( 2224)     ks = 1                              ! Depth writes always start at k=1
( 2225)     is = nimpp                          ! Get global i-index for this subdomain
( 2226)     js = njmpp                          ! Get global j-index  "   "      "
( 2227)     iz = nlei                           ! Get i-size of subdomain section to be written
( 2228)     jz = nlej                           ! Get j-size of subdomain section to be written
( 2229)     ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,3),1:SIZE(data,4)))      ! Allocate send/recieve buffer
( 2230)     buf = data(1:iz,1:jz,1:SIZE(data,3), 1:SIZE(data,4))         ! Copy section of subdomain to buffer
( 2231)     IF(tstep >= 0) THEN
( 2232)        tindex = tstep/nwrite - nit000/nwrite
( 2233)     ELSE
( 2234)        tindex = tstep * (-1)
( 2235)     END IF
( 2236)     status = NCDF_NOERR
( 2237)     hastime = .FALSE.
( 2238) 
( 2239)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2240)        WRITE(100,*) 'NCDF DEBUG: Writing 4D array to file:', filename
( 2241)     END IF
( 2242) 
( 2243)     ! Open netCDF file and get info
( 2244)     IF(nproc == 0) THEN
( 2245)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 2246)           WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Opening file:', filename
( 2247)           CALL FLUSH
( 2248)        END IF
( 2249)        nfstat = nf90_open(filename, nf90_write, ncid)
( 2250)        IF(nfstat /= nf90_noerr) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 40

( 2251)           status = NCDF_NFERR
( 2252)           RETURN
( 2253)        END IF
( 2254)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 2255)           WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Getting info from file:', filename
( 2256)           CALL FLUSH
( 2257)        END IF
( 2258)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 2259)        IF(nfstat /= nf90_noerr) THEN
( 2260)           status = NCDF_NFERR
( 2261)           RETURN
( 2262)        END IF
( 2263) 
( 2264)        ! Determine if this variable contains a time axis
( 2265)        nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 2266)             ndims=ndims, dimids=var_dimids)
( 2267)        IF(nfstat /= nf90_noerr) THEN
( 2268)           status = NCDF_NFERR
( 2269)           RETURN
( 2270)        END IF
( 2271)        DO i = 1, ndims
( 2272)           nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2273)                name=dname, len=k)
( 2274) !!DB: next line not correct 
( 2275) !          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2276) !               name=dname)
( 2277)           dname=TRIM(dname)
( 2278)           IF(nfstat /= nf90_noerr) THEN
( 2279)              status = NCDF_NFERR
( 2280)              RETURN
( 2281)           END IF
( 2282)           IF(dname=='time_counter') THEN
( 2283)              hastime = .TRUE.
( 2284)              EXIT
( 2285)           END IF
( 2286)        END DO
( 2287)     END IF
( 2288) 
( 2289)     ! Get data from each processor in turn, processor 0 writes to file
( 2290)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2291)        WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Writing data to file:', filename
( 2292)        CALL FLUSH
( 2293)     END IF
( 2294)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2295)     IF(mpistat /= 0) THEN
( 2296)        status = NCDF_MPERR
( 2297)        RETURN
( 2298)     END IF
( 2299)    
( 2300)     DO i=0, (jpnij - 1)
( 2301)        ! Since processor 0 doesn't need to send data to itself, only
( 2302)        ! processors with IDs > 0 make send calls
( 2303)        IF(i > 0) THEN
( 2304)           ! If it's the local processor's turn to send, go ahead and
( 2305)           ! send data and write indices to processor 0
( 2306)           IF(nproc == i) THEN
( 2307)              CALL MPI_SEND(is, 1, MPI_INTEGER, &
( 2308)                   0, i, MPI_COMM_WORLD, mpistat)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 41

( 2309)              IF(mpistat /= 0) THEN
( 2310)                 status = NCDF_MPERR
( 2311)                 RETURN
( 2312)              END IF
( 2313)              CALL MPI_SEND(js, 1, MPI_INTEGER, &
( 2314)                   0, i+1, MPI_COMM_WORLD, mpistat)
( 2315)              IF(mpistat /= 0) THEN
( 2316)                 status = NCDF_MPERR
( 2317)                 RETURN
( 2318)              END IF
( 2319)              CALL MPI_SEND(ks, 1, MPI_INTEGER, &
( 2320)                   0, i+2, MPI_COMM_WORLD, mpistat)
( 2321)              IF(mpistat /= 0) THEN
( 2322)                 status = NCDF_MPERR
( 2323)                 RETURN
( 2324)              END IF
( 2325)              CALL MPI_SEND(iz, 1, MPI_INTEGER, &
( 2326)                   0, i+3, MPI_COMM_WORLD, mpistat)
( 2327)              IF(mpistat /= 0) THEN
( 2328)                 status = NCDF_MPERR
( 2329)                 RETURN
( 2330)              END IF
( 2331)              CALL MPI_SEND(jz, 1, MPI_INTEGER, &
( 2332)                   0, i+4, MPI_COMM_WORLD, mpistat)
( 2333)              IF(mpistat /= 0) THEN
( 2334)                 status = NCDF_MPERR
( 2335)                 RETURN
( 2336)              END IF
( 2337)              CALL MPI_SEND(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 2338)                   0, i+5, MPI_COMM_WORLD, mpistat)
( 2339)              IF(mpistat /= 0) THEN
( 2340)                 status = NCDF_MPERR
( 2341)                 RETURN
( 2342)              END IF
( 2343)           END IF
( 2344)        END IF
( 2345) 
( 2346)        ! Processor 0 receives data and writes it to file
( 2347)        IF(nproc == 0) THEN
( 2348)           ! Since processor 0 doesn't need to receive data from itself, 
( 2349)           ! we only make receive calls if it's another processor's turn
( 2350)           ! to send data
( 2351)           IF(i > 0) THEN
( 2352)              CALL MPI_RECV(is, 1, MPI_INTEGER, &
( 2353)                   i, i, MPI_COMM_WORLD, rcvstat, mpistat)
( 2354)              IF(mpistat /= 0) THEN
( 2355)                 status = NCDF_MPERR
( 2356)                 RETURN
( 2357)              END IF
( 2358)              CALL MPI_RECV(js, 1, MPI_INTEGER, &
( 2359)                   i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)
( 2360)              IF(mpistat /= 0) THEN
( 2361)                 status = NCDF_MPERR
( 2362)                 RETURN
( 2363)              END IF
( 2364)              CALL MPI_RECV(ks, 1, MPI_INTEGER, &
( 2365)                   i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
( 2366)              IF(mpistat /= 0) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 42

( 2367)                 status = NCDF_MPERR
( 2368)                 RETURN
( 2369)              END IF
( 2370)              CALL MPI_RECV(iz, 1, MPI_INTEGER, &
( 2371)                   i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
( 2372)              IF(mpistat /= 0) THEN
( 2373)                 status = NCDF_MPERR
( 2374)                 RETURN
( 2375)              END IF
( 2376)              CALL MPI_RECV(jz, 1, MPI_INTEGER, &
( 2377)                   i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
( 2378)              IF(mpistat /= 0) THEN
( 2379)                 status = NCDF_MPERR
( 2380)                 RETURN
( 2381)              END IF
( 2382)              DEALLOCATE(buf)
( 2383)              ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,3),1:SIZE(data,4)))
( 2384)              CALL MPI_RECV(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 2385)                   i, i+5, MPI_COMM_WORLD, rcvstat, mpistat)
( 2386)              IF(mpistat /= 0) THEN
( 2387)                 status = NCDF_MPERR
( 2388)                 RETURN
( 2389)              END IF
( 2390)           END IF
( 2391)           ! Got the data, write it to the proper location in the netCDF file
( 2392)           IF(hastime .EQV. .TRUE.) THEN
( 2393)              nfstat = nf90_put_var(ncid, varid, buf, &
( 2394)                   (/ is, js, ks, 1, tindex/), &
( 2395)                   (/ iz, jz, SIZE(data,3), SIZE(data,4), 1 /))
( 2396)           ELSE
( 2397)              nfstat = nf90_put_var(ncid, varid, buf, &
( 2398)                   (/ is, js, ks, 1 /), &
( 2399)                   (/ iz, jz, SIZE(data,3), SIZE(data,4) /))
( 2400)           END IF
( 2401)           IF(nfstat /= nf90_noerr) THEN
( 2402)              status = NCDF_NFERR
( 2403)              RETURN
( 2404)           END IF
( 2405)        END IF
( 2406) 
( 2407)        ! Sync all processors at the end of each loop iteration to prevent
( 2408)        ! concurrency issues/race conditions/various MPI-related ugliness
( 2409)        CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2410)        IF(mpistat /= 0) THEN
( 2411)           status = NCDF_MPERR
( 2412)           RETURN
( 2413)        END IF
( 2414)     END DO
( 2415) 
( 2416)     ! All done, close up the netCDF file
( 2417)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2418)        WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Closing file:', filename
( 2419)        CALL FLUSH
( 2420)     END IF
( 2421)     IF(nproc == 0) THEN
( 2422)        nfstat = nf90_close(ncid)
( 2423)        IF(nfstat /= nf90_noerr) THEN
( 2424)           status = NCDF_NFERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 43

( 2425)           RETURN
( 2426)        END IF
( 2427)     END IF
( 2428)     DEALLOCATE(buf)
( 2429) 
( 2430)     ! Sync up processors before returning from subroutine
( 2431)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2432)     IF(mpistat /= 0) THEN
( 2433)        status = NCDF_MPERR
( 2434)        RETURN
( 2435)     END IF
( 2436) 
( 2437)   END SUBROUTINE ncdf_write4d
( 2438) 
( 2439)   ! The inverse of ncdf_writesv. Reads a scalar value from a 1-D array in a
( 2440)   ! netCDF file. All processors will get the same value.
( 2441)   ! filename - file to read from
( 2442)   ! varname - variable to read from
( 2443)   ! data - scalar to put the data into (must be a REAL)
( 2444)   ! index -the index to read from
( 2445)   ! status - return status of the subroutine
( 2446)   SUBROUTINE ncdf_readsv(filename, varname, data, index, status)
( 2447)     IMPLICIT NONE
( 2448)     ! Subroutine argument declarations
( 2449)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2450)     INTEGER,INTENT(IN) :: index
( 2451)     REAL,INTENT(OUT) :: data
( 2452)     INTEGER,INTENT(OUT) :: status
( 2453) 
( 2454)     ! Local declarations
( 2455)     INTEGER :: ncid,    &  ! netCDF file ID
( 2456)                varid,   &  ! ID of netCDF variable to be written to
( 2457)                nfstat,  &  ! netCDF library call return status
( 2458)                mpistat     ! MPI library call return status
( 2459)     INTEGER,DIMENSION(1:5) :: count
( 2460) 
( 2461)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2462)        WRITE(100,*) 'NCDF DEBUG: Reading single value from file:', filename
( 2463)        CALL FLUSH
( 2464)     END IF
( 2465) 
( 2466)     ! Initializations
( 2467)     status = NCDF_NOERR
( 2468)     count(1) = 1
( 2469) 
( 2470)     ! Open netCDF file and get info
( 2471) !!DB: think that there should be no if(nproc==0) then ...
( 2472) !    IF(nproc == 0) THEN
( 2473)        IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2474)           WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Opening file:', filename
( 2475)           CALL FLUSH
( 2476)        END IF
( 2477) !DB
( 2478) !       nfstat = nf90_open(filename, nf90_write, ncid)
( 2479)        nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2480)        IF(nfstat /= nf90_noerr) THEN
( 2481)           status = NCDF_NFERR
( 2482)           RETURN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 44

( 2483)        END IF
( 2484)        IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2485)           WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Getting info from file:', filename
( 2486)           CALL FLUSH
( 2487)        END IF
( 2488)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 2489)        IF(nfstat /= nf90_noerr) THEN
( 2490)           status = NCDF_NFERR
( 2491)           RETURN
( 2492)        END IF
( 2493)        
( 2494)        ! Read data from netCDF file
( 2495)        ! This subroutine assumes all processors will get the same data from this
( 2496)        ! variable
( 2497)        IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2498)           WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Reading data from file:', filename
( 2499)           CALL FLUSH
( 2500)        END IF
( 2501)        nfstat = nf90_get_var(ncid, varid, data, (/ index /))
( 2502)        IF(nfstat /= nf90_noerr) THEN
( 2503)           status = NCDF_NFERR
( 2504)           RETURN
( 2505)        END IF
( 2506)        IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2507)           WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Closing file:', filename
( 2508)           CALL FLUSH
( 2509)        END IF
( 2510)        nfstat = nf90_close(ncid)
( 2511)        IF(nfstat /= nf90_noerr) THEN
( 2512)           status = NCDF_NFERR
( 2513)           RETURN
( 2514)        END IF
( 2515) !    END IF
( 2516) 
( 2517)     ! Sync up processors before returning from subroutine
( 2518)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2519)     IF(mpistat /= 0) THEN
( 2520)        status = NCDF_MPERR
( 2521)        RETURN
( 2522)     END IF
( 2523) 
( 2524)   END SUBROUTINE ncdf_readsv
( 2525) 
( 2526)   ! The inverse of ncdf_write1d. Reads a 1-D array from a netCDF file
( 2527)   ! All processors will get the same values.
( 2528)   ! filename - file to read from
( 2529)   ! varname - variable to read from
( 2530)   ! data - 1-D array to put the data into (must be a REAL)
( 2531)   ! status - return status of the subroutine
( 2532)   SUBROUTINE ncdf_read1d(filename, varname, data, status, tide)
( 2533)     IMPLICIT NONE
( 2534)     ! Subroutine argument declarations
( 2535)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2536)     REAL,DIMENSION(:),INTENT(OUT) :: data
( 2537)     INTEGER,INTENT(OUT) :: status
( 2538)     INTEGER,INTENT(IN),OPTIONAL :: tide
( 2539) 
( 2540)     ! Local declarations






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 45

( 2541)     INTEGER :: ncid,    &  ! netCDF file ID
( 2542)                varid,   &  ! ID of netCDF variable to be written to
( 2543)                nfstat,  &  ! netCDF library call return status
( 2544)                mpistat     ! MPI library call return status
( 2545)     
( 2546) 
( 2547) 
( 2548)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2549)        WRITE(100,*) 'NCDF DEBUG: Reading 1D array from file:', filename
( 2550)        CALL FLUSH
( 2551)     END IF
( 2552) 
( 2553)     ! Initializations
( 2554)     status = NCDF_NOERR
( 2555) 
( 2556)     ! Open netCDF file and get info
( 2557)     
( 2558)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2559)        WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Opening file:', filename
( 2560)        CALL FLUSH
( 2561)     END IF
( 2562)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2563)     IF(nfstat /= nf90_noerr) THEN
( 2564)        status = NCDF_NFERR
( 2565)        RETURN
( 2566)     END IF
( 2567)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2568)        WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Getting info from file:', filename
( 2569)        CALL FLUSH
( 2570)     END IF
( 2571)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 2572)     IF(nfstat /= nf90_noerr) THEN
( 2573)        status = NCDF_NFERR
( 2574)        RETURN
( 2575)     END IF
( 2576)     
( 2577)     ! Read data from netCDF file
( 2578)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2579)        WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Reading data from file:', filename
( 2580)        CALL FLUSH
( 2581)     END IF
( 2582)     nfstat = nf90_get_var(ncid, varid, data, (/ 1 /), (/ SIZE(data) /))
( 2583)     IF(nfstat /= nf90_noerr) THEN
( 2584)        status = NCDF_NFERR
( 2585)        CALL FLUSH
( 2586)        RETURN
( 2587)     END IF
( 2588)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2589)        WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Closing file:', filename
( 2590)        CALL FLUSH
( 2591)     END IF
( 2592)     nfstat = nf90_close(ncid)
( 2593)     IF(nfstat /= nf90_noerr) THEN
( 2594)        status = NCDF_NFERR
( 2595)        RETURN
( 2596)     END IF
( 2597)  
( 2598)     ! Sync up processors before returning from subroutine






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 46

( 2599)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2600)     IF(mpistat /= 0) THEN
( 2601)        status = NCDF_MPERR
( 2602)        RETURN
( 2603)     END IF
( 2604)     
( 2605)   END SUBROUTINE ncdf_read1d
( 2606) 
( 2607)   ! The inverse of ncdf_write2d. Reads a 2-D array from a netCDF file.
( 2608)   ! Each processor only reads its local subdomain.
( 2609)   ! filename - file to read from
( 2610)   ! varname - variable to read from
( 2611)   ! data - 2-D array to put the data into (must be a REAL)
( 2612)   ! time - time index to read (If the target variable in the netCDF file is
( 2613)   !        actually a 3-D variable with a time axis, ignored otherwise)
( 2614)   ! status - return status of the subroutine
( 2615)   SUBROUTINE ncdf_read2d(filename, varname, data, time, status)
( 2616)     IMPLICIT NONE
( 2617)     ! Subroutine argument declarations
( 2618)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2619)     INTEGER,INTENT(IN) :: time
( 2620)     REAL,DIMENSION(:,:),INTENT(OUT) :: data
( 2621)     INTEGER,INTENT(OUT) :: status
( 2622) 
( 2623)     ! Local declarations
( 2624)     INTEGER :: ncid,    &  ! netCDF file ID
( 2625)                varid,   &  ! ID of netCDF variable to be written to
( 2626)                mpistat, &  ! MPI library call return status
( 2627)                i,       &
( 2628)                j,       &
( 2629)                k,       &
( 2630)                nfstat,  &  ! netCDF library call return status
( 2631)                ndims       ! Number of dimensions in this variable
( 2632)     
( 2633)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 2634)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 2635)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 2636)     REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: databuf
( 2637)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 2638)     INTEGER :: tindex
( 2639) 
( 2640)     ! Initializations
( 2641)     status = NCDF_NOERR
( 2642)     hastime = .FALSE.
( 2643)     j = 0
( 2644)     var_dimlens = 0
( 2645)     if(time > 0) then
( 2646)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 2647)     else
( 2648)        tindex = -time
( 2649)     endif
( 2650) 
( 2651) 
( 2652)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2653)        WRITE(100,*) 'NCDF DEBUG: Reading 2D array from file:', filename
( 2654)        CALL FLUSH
( 2655)     END IF
( 2656) 






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 47

( 2657)     ! Open netCDF file and get info
( 2658)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2659)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Opening file:', filename
( 2660)        CALL FLUSH
( 2661)     END IF
( 2662) 
( 2663)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2664)     IF(nfstat /= nf90_noerr) THEN
( 2665)        status = NCDF_NFERR
( 2666)        RETURN
( 2667)     END IF
( 2668)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2669)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Getting info from file:', filename
( 2670)        CALL FLUSH
( 2671)     END IF
( 2672)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 2673)     IF(nfstat /= nf90_noerr) THEN
( 2674)        status = NCDF_NFERR
( 2675)        RETURN
( 2676)     END IF
( 2677)     
( 2678)     ! Determine if this variable contains a time axis
( 2679)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 2680)          ndims=ndims, dimids=var_dimids)
( 2681)     IF(nfstat /= nf90_noerr) THEN
( 2682)        status = NCDF_NFERR
( 2683)        RETURN
( 2684)     END IF
( 2685)     DO i = 1, ndims
( 2686) !!DB 2008.07.07
( 2687) !       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2688)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2689)             name=dname, len=k)
( 2690)        dname=TRIM(dname)
( 2691)        var_dimlens(i - j) = k
( 2692)        IF(nfstat /= nf90_noerr) THEN
( 2693)           status = NCDF_NFERR
( 2694)           RETURN
( 2695)        END IF
( 2696)        IF((dname=='time_counter') .OR. (dname=='time')) THEN
( 2697)           hastime = .TRUE.
( 2698)           j = 1
( 2699)        END IF
( 2700)     END DO
( 2701)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2)))
( 2702)     
( 2703)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2704)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Reading record : ', tindex, 'from file:', filename
( 2705)        CALL FLUSH
( 2706)     END IF
( 2707)     
( 2708)     ! Read the appropriate subdomain for each processor
( 2709)     IF(hastime .EQV. .TRUE.) THEN
( 2710)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 2711)             (/ 1, 1, tindex/), &
( 2712)             (/ var_dimlens(1), var_dimlens(2), 1 /))
( 2713)     ELSE
( 2714)        nfstat = nf90_get_var(ncid, varid, databuf, &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 48

( 2715)             (/ 1, 1 /), &
( 2716)             (/ var_dimlens(1), var_dimlens(2) /))
( 2717)     END IF
( 2718) 
( 2719)     DO j = 1, nlcj
( 2720)        DO i = 1, nlci
( 2721)           data(i, j) = databuf(mig(i), mjg(j))
( 2722)        END DO
( 2723)     END DO
( 2724)     
( 2725) 
( 2726)     ! All done, close up the netCDF file
( 2727)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2728)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Closing file:', filename
( 2729)        CALL FLUSH
( 2730)     END IF
( 2731)     nfstat = nf90_close(ncid)
( 2732)     IF(nfstat /= nf90_noerr) THEN
( 2733)        status = NCDF_NFERR
( 2734)        RETURN
( 2735)     END IF
( 2736)     DEALLOCATE(databuf)
( 2737)     
( 2738)     ! Sync up processors before returning from subroutine
( 2739)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2740)     IF(mpistat /= 0) THEN
( 2741)        status = NCDF_MPERR
( 2742)        RETURN
( 2743)     END IF
( 2744) 
( 2745)   END SUBROUTINE ncdf_read2d
( 2746) 
( 2747)   ! The inverse of ncdf_write3d. Reads a 3-D array from a netCDF file.
( 2748)   ! Each processor only reads its local subdomain.
( 2749)   ! filename - file to read from
( 2750)   ! varname - variable to read from
( 2751)   ! data - 3-D array to put the data into (must be a REAL)
( 2752)   ! time - time index to read (If the target variable in the netCDF file is
( 2753)   !        actually a 4-D variable with a time axis, ignored otherwise)
( 2754)   ! status - return status of the subroutine
( 2755)   SUBROUTINE ncdf_read3d(filename, varname, data, time, status)
( 2756)     IMPLICIT NONE
( 2757)     ! Subroutine argument declarations
( 2758)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2759)     INTEGER,INTENT(IN) :: time
( 2760)     REAL,DIMENSION(:,:,:),INTENT(OUT) :: data
( 2761)     INTEGER,INTENT(OUT) :: status
( 2762) 
( 2763)     ! Local declarations
( 2764)     INTEGER :: ncid,    &  ! netCDF file ID
( 2765)                varid,   &  ! ID of netCDF variable to be written to
( 2766)                mpistat, &  ! MPI library call return status
( 2767)                nfstat,  &  ! netCDF library call return status
( 2768)                i,       &  ! Loop counter
( 2769)                j,       &  !     "
( 2770)                k,       &  !     "
( 2771)                ndims       ! Number of dimensions in this variable
( 2772) 






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 49

( 2773)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 2774)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 2775)     INTEGER :: tindex
( 2776)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 2777)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 2778)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: databuf
( 2779) 
( 2780)     ! Initializations
( 2781)     status = NCDF_NOERR
( 2782)     hastime = .FALSE.
( 2783)     var_dimlens = 0
( 2784)     j = 0
( 2785)     if(time > 0) then
( 2786)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 2787)     else
( 2788)        tindex = -time
( 2789)     endif
( 2790) 
( 2791)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2792)        WRITE(100,*) 'NCDF DEBUG: Reading 3D array from file:', filename
( 2793)        CALL FLUSH
( 2794)     END IF
( 2795) 
( 2796)     ! Open netCDF file and get info
( 2797)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2798)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Opening file:', filename
( 2799)        CALL FLUSH
( 2800)     END IF
( 2801) 
( 2802)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2803)     IF(nfstat /= nf90_noerr) THEN
( 2804)        status = NCDF_NFERR
( 2805)        RETURN
( 2806)     END IF
( 2807)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2808)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Getting info from file:', filename
( 2809)        CALL FLUSH
( 2810)     END IF
( 2811)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 2812)     IF(nfstat /= nf90_noerr) THEN
( 2813)        status = NCDF_NFERR
( 2814)        RETURN
( 2815)     END IF
( 2816)     
( 2817)     ! Determine if this variable contains a time axis
( 2818)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 2819)          ndims=ndims, dimids=var_dimids)
( 2820)     IF(nfstat /= nf90_noerr) THEN
( 2821)        status = NCDF_NFERR
( 2822)        RETURN
( 2823)     END IF
( 2824)     DO i = 1, ndims
( 2825)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2826)             name=dname, len=k)
( 2827) !       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2828) !            name=dname, len=k)
( 2829)        dname=TRIM(dname)
( 2830)        var_dimlens(i-j) = k






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 50

( 2831)        IF(nfstat /= nf90_noerr) THEN
( 2832)           status = NCDF_NFERR
( 2833)           RETURN
( 2834)        END IF
( 2835)        IF((dname=='time_counter') .OR. (dname=='time')) THEN
( 2836)           hastime = .TRUE.
( 2837)           j = 1
( 2838)           EXIT
( 2839)        END IF
( 2840)     END DO
( 2841)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3)))
( 2842) 
( 2843)     ! Read the appropriate subdomain for each processor
( 2844)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2845)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Reading data from file:', filename
( 2846)        CALL FLUSH
( 2847)     END IF
( 2848)     IF(hastime .EQV. .TRUE.) THEN
( 2849)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 2850)             (/ 1, 1, 1, tindex/), &
( 2851)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), 1 /))
( 2852)     ELSE
( 2853)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 2854)             (/ 1, 1, 1 /), &
( 2855)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3) /))
( 2856)     END IF
( 2857)     IF(nfstat /= nf90_noerr) THEN
( 2858)        status = NCDF_NFERR
( 2859)        RETURN
( 2860)     END IF
( 2861) 
( 2862)     DO k = 1, var_dimlens(3)
( 2863)        DO j = 1, nlcj
( 2864)           DO i = 1, nlci
( 2865)              data(i, j, k) = databuf(mig(i), mjg(j), k)
( 2866)           END DO
( 2867)        END DO
( 2868)     END DO
( 2869) 
( 2870) 
( 2871)     ! All done, close up the netCDF file
( 2872)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2873)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Closing file:', filename
( 2874)        CALL FLUSH
( 2875)     END IF
( 2876)     nfstat = nf90_close(ncid)
( 2877)     IF(nfstat /= nf90_noerr) THEN
( 2878)        status = NCDF_NFERR
( 2879)        RETURN
( 2880)     END IF
( 2881) 
( 2882) 
( 2883)     DEALLOCATE(databuf)
( 2884) 
( 2885)     ! Sync up processors before returning from subroutine
( 2886)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2887)     IF(mpistat /= 0) THEN
( 2888)        status = NCDF_MPERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 51

( 2889)        RETURN
( 2890)     END IF
( 2891) 
( 2892)   END SUBROUTINE ncdf_read3d
( 2893) 
( 2894)   ! The inverse of ncdf_write3d. Reads a 3-D array from a netCDF file.
( 2895)   ! Each processor only reads its local subdomain.
( 2896)   ! filename - file to read from
( 2897)   ! varname - variable to read from
( 2898)   ! data - 3-D array to put the data into (must be a REAL)
( 2899)   ! dsz - Dimension SiZe - size of the 4th dimension in the variable
( 2900)   ! time - time index to read (If the target variable in the netCDF file is
( 2901)   !        actually a 4-D variable with a time axis, ignored otherwise)
( 2902)   ! status - return status of the subroutine
( 2903)   ! NOT EXTENSIVELY TESTED - MAY STILL BE BUGGY
( 2904)   SUBROUTINE ncdf_read4d(filename, varname, data, dsz, time, status)
( 2905)     IMPLICIT NONE
( 2906)     ! Subroutine argument declarations
( 2907)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2908)     INTEGER,INTENT(IN) :: time, dsz
( 2909)     REAL,DIMENSION(:,:,:,:),INTENT(OUT) :: data
( 2910)     INTEGER,INTENT(OUT) :: status
( 2911)     
( 2912)     ! Local declarations
( 2913)     INTEGER :: ncid,    &  ! netCDF file ID
( 2914)                varid,   &  ! ID of netCDF variable to be written to
( 2915)                mpistat, &  ! MPI library call return status
( 2916)                nfstat,  &  ! netCDF library call return status
( 2917)                i,k,       &  ! Loop counter
( 2918)                ndims       ! Number of dimensions in this variable
( 2919) 
( 2920)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 2921)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 2922)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 2923)     INTEGER :: tindex
( 2924)     
( 2925)     ! Initializations
( 2926)     status = NCDF_NOERR
( 2927)     hastime = .FALSE.
( 2928) 
( 2929)     if(time > 0) then
( 2930)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 2931)     else
( 2932)        tindex = -time
( 2933)     endif
( 2934) !    tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 2935)     
( 2936)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2937)        WRITE(100,*) 'NCDF DEBUG: Reading 4D array from file:', filename
( 2938)        CALL FLUSH
( 2939)     END IF
( 2940)     
( 2941)     ! Open netCDF file and get info
( 2942)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2943)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Opening file:', filename
( 2944)        CALL FLUSH
( 2945)     END IF
( 2946) 






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 52

( 2947)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2948)     IF(nfstat /= nf90_noerr) THEN
( 2949)        status = NCDF_NFERR
( 2950)        RETURN
( 2951)     END IF
( 2952)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2953)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Getting info from file:', filename
( 2954)        CALL FLUSH
( 2955)     END IF
( 2956)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 2957)     IF(nfstat /= nf90_noerr) THEN
( 2958)        status = NCDF_NFERR
( 2959)        RETURN
( 2960)     END IF
( 2961)     
( 2962)     ! Determine if this variable contains a time axis
( 2963)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 2964)          ndims=ndims, dimids=var_dimids)
( 2965)     IF(nfstat /= nf90_noerr) THEN
( 2966)        status = NCDF_NFERR
( 2967)        RETURN
( 2968)     END IF
( 2969)     DO i = 1, ndims
( 2970)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2971)             name=dname, len=k)
( 2972) !       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2973) !            name=dname)
( 2974)        dname=TRIM(dname)
( 2975)        IF(nfstat /= nf90_noerr) THEN
( 2976)           status = NCDF_NFERR
( 2977)           RETURN
( 2978)        END IF
( 2979)        IF(dname=='time_counter') THEN
( 2980)           hastime = .TRUE.
( 2981)           EXIT
( 2982)        END IF
( 2983)     END DO
( 2984)     
( 2985)     ! Got the data, write it to the proper location in the netCDF file
( 2986)     IF(hastime .EQV. .TRUE.) THEN
( 2987)        nfstat = nf90_get_var(ncid, varid, data, &
( 2988)             (/ nimpp, njmpp, 1, 1, tindex/), &
( 2989)             (/ jpi, jpj, SIZE(data,3)-1, SIZE(data,4)-1, 1 /))
( 2990)     ELSE
( 2991)        nfstat = nf90_get_var(ncid, varid, data, &
( 2992)             (/ nimpp, njmpp, 1, 1 /), &
( 2993)             (/ jpj, jpi, SIZE(data,3)-1, SIZE(data,4)-1 /))
( 2994)     END IF
( 2995)     IF(nfstat /= nf90_noerr) THEN
( 2996)        status = NCDF_NFERR
( 2997)        RETURN
( 2998)     END IF
( 2999)     
( 3000)     ! All done, close up the netCDF file
( 3001)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3002)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Closing file:', filename
( 3003)        CALL FLUSH
( 3004)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 53

( 3005)     nfstat = nf90_close(ncid)
( 3006)     IF(nfstat /= nf90_noerr) THEN
( 3007)        status = NCDF_NFERR
( 3008)        RETURN
( 3009)     END IF
( 3010)     
( 3011)     ! Sync up processors before returning from subroutine
( 3012)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3013)     IF(mpistat /= 0) THEN
( 3014)        status = NCDF_MPERR
( 3015)        RETURN
( 3016)     END IF
( 3017) 
( 3018)   END SUBROUTINE ncdf_read4d
( 3019) 
( 3020)   ! ncdf_create_file creates an empty netCDF dataset with the specified filename
( 3021)   ! filename - name of the file to be created
( 3022)   ! status - return status
( 3023)   SUBROUTINE ncdf_create_file(filename, status)
( 3024)     IMPLICIT NONE
( 3025)     
( 3026)     ! Subroutine argument declarations
( 3027)     CHARACTER(LEN=*),INTENT(IN) :: filename
( 3028)     INTEGER,INTENT(OUT) :: status
( 3029) 
( 3030)     ! Local declarations
( 3031)     INTEGER :: ncid,   &  ! file ID
( 3032)                nfstat, &  ! netCDF call return status
( 3033)                mpistat    ! MPI call return status
( 3034) 
( 3035)     ! Initialization
( 3036)     status = NCDF_NOERR
( 3037) 
( 3038)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3039)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_file - Creating:', filename
( 3040)        CALL FLUSH
( 3041)     END IF
( 3042) 
( 3043)     IF(nproc == 0) THEN
( 3044)        nfstat = nf90_create(filename, NF90_CLOBBER, ncid)
( 3045)        IF(nfstat /= nf90_noerr) THEN
( 3046)           status = NCDF_NFERR
( 3047)           RETURN
( 3048)        END IF
( 3049)        nfstat = nf90_close(ncid)
( 3050)        IF(nfstat /= nf90_noerr) THEN
( 3051)           status = NCDF_NFERR
( 3052)           RETURN
( 3053)        END IF
( 3054)     END IF
( 3055) 
( 3056)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3057)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_file - Done creating:', filename
( 3058)        CALL FLUSH
( 3059)     END IF
( 3060) 
( 3061)     ! Sync up processors before returning from subroutine
( 3062)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 54

( 3063)     IF(mpistat /= 0) THEN
( 3064)        status = NCDF_MPERR
( 3065)        RETURN
( 3066)     END IF
( 3067) 
( 3068)   END SUBROUTINE ncdf_create_file
( 3069) 
( 3070)   ! ncdf_create_dim creates a new dimension in an existing dataset
( 3071)   ! filename - name of the file to modified
( 3072)   ! dimname - name of the dimension to create
( 3073)   ! dim_len - length of the new dimension
( 3074)   ! status - return status
( 3075)   SUBROUTINE ncdf_create_dim(filename, dimname, dim_len, status)
( 3076)     IMPLICIT NONE
( 3077) 
( 3078)     ! Subroutine argument declarations
( 3079)     CHARACTER(LEN=*),INTENT(IN) :: filename, dimname
( 3080)     INTEGER,INTENT(IN) :: dim_len
( 3081)     INTEGER,INTENT(OUT) :: status
( 3082) 
( 3083)     ! Local declarations
( 3084)     INTEGER :: ncid,    &  ! file ID
( 3085)                nfstat,  &  ! netCDF call return status
( 3086)                mpistat, &  ! MPI call return status
( 3087)                dimid,   &  ! New dimension ID 
( 3088)                dimlen
( 3089) 
( 3090)     ! Initialization
( 3091)     status = NCDF_NOERR
( 3092)     dimlen = dim_len
( 3093) 
( 3094)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3095)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_dim - Creating dimension in file:', filename
( 3096)        CALL FLUSH
( 3097)     END IF
( 3098) 
( 3099)     IF(nproc == 0) THEN
( 3100)        IF(dimlen == -1) THEN
( 3101)           dimlen = NF90_UNLIMITED
( 3102)        END IF
( 3103)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3104)        IF(nfstat /= nf90_noerr) THEN
( 3105)           status = NCDF_NFERR
( 3106)           RETURN
( 3107)        END IF
( 3108)        nfstat = nf90_redef(ncid)
( 3109)        IF(nfstat /= nf90_noerr) THEN
( 3110)           status = NCDF_NFERR
( 3111)           RETURN
( 3112)        END IF
( 3113)        nfstat = nf90_def_dim(ncid, dimname, dimlen, dimid)
( 3114)        IF(nfstat /= nf90_noerr) THEN
( 3115)           status = NCDF_NFERR
( 3116)           RETURN
( 3117)        END IF
( 3118)        nfstat = nf90_close(ncid)
( 3119)        IF(nfstat /= nf90_noerr) THEN
( 3120)           status = NCDF_NFERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 55

( 3121)           RETURN
( 3122)        END IF
( 3123)     END IF
( 3124) 
( 3125)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3126)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_dim - Done creating dimension in file:', filename
( 3127)        CALL FLUSH
( 3128)     END IF
( 3129)     
( 3130)     ! Sync up processors before returning from subroutine
( 3131)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3132)     IF(mpistat /= 0) THEN
( 3133)        status = NCDF_MPERR
( 3134)        RETURN
( 3135)     END IF
( 3136) 
( 3137)   END SUBROUTINE ncdf_create_dim
( 3138) 
( 3139)   ! ncdf_create_var creates a new variable in an existing dataset
( 3140)   ! filename - name of the file to modified
( 3141)   ! varname - name of the variable to create
( 3142)   ! dimname - an array containing the names (in order) of the dimensions
( 3143)   !           the new variable will use
( 3144)   ! datatype - type of the variable (NCDF_FLOAT or NCDF_DOUBLE)
( 3145)   ! status - return status
( 3146)   SUBROUTINE ncdf_create_var(filename, varname, dimnames, datatype, status)
( 3147)     IMPLICIT NONE
( 3148)     
( 3149)     ! Subroutine argument declarations
( 3150)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 3151)     CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: dimnames
( 3152)     INTEGER,INTENT(IN) :: datatype
( 3153)     INTEGER,INTENT(OUT) :: status
( 3154) 
( 3155)     ! Local declarations
( 3156)     INTEGER :: ncid,    &  ! file ID
( 3157)                nfstat,  &  ! netCDF call return status
( 3158)                mpistat, &  ! MPI call return status
( 3159)                varid,   &  ! New variable ID
( 3160)                i
( 3161)     INTEGER,DIMENSION(1:NF90_MAX_VAR_DIMS) :: dimids
( 3162) 
( 3163)     ! Initialization
( 3164)     status = NCDF_NOERR
( 3165)     dimids = -1
( 3166) 
( 3167)     IF((datatype /= NCDF_FLOAT) .AND. (datatype /= NCDF_DOUBLE)) THEN
( 3168)        status = NCDF_ARERR
( 3169)        RETURN
( 3170)     END IF
( 3171) 
( 3172)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3173)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_var - Creating variable in file:', filename
( 3174)        CALL FLUSH
( 3175)     END IF
( 3176) 
( 3177)     IF(nproc == 0) THEN
( 3178)        nfstat = nf90_open(filename, nf90_write, ncid)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 56

( 3179)        IF(nfstat /= nf90_noerr) THEN
( 3180)           status = NCDF_NFERR
( 3181)           RETURN
( 3182)        END IF
( 3183)        ! Get dimension IDs
( 3184)        DO i = 1, SIZE(dimnames)
( 3185)           nfstat = nf90_inq_dimid(ncid, dimnames(i), dimids(i))
( 3186)           IF(nfstat /= nf90_noerr) THEN
( 3187)              status = NCDF_NFERR
( 3188)              RETURN
( 3189)           END IF
( 3190)        END DO
( 3191)        ! Create variable
( 3192)        nfstat = nf90_redef(ncid)
( 3193)        IF(nfstat /= nf90_noerr) THEN
( 3194)           status = NCDF_NFERR
( 3195)           RETURN
( 3196)        END IF
( 3197)        nfstat = nf90_def_var(ncid, varname, datatype, &
( 3198)             dimids(1:SIZE(dimnames)), varid)
( 3199)        IF(nfstat /= nf90_noerr) THEN
( 3200)           status = NCDF_NFERR
( 3201)           RETURN
( 3202)        END IF
( 3203)        nfstat = nf90_close(ncid)
( 3204)        IF(nfstat /= nf90_noerr) THEN
( 3205)           status = NCDF_NFERR
( 3206)           RETURN
( 3207)        END IF
( 3208)     END IF
( 3209) 
( 3210)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3211)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_var - Done creating variable in file:', filename
( 3212)        CALL FLUSH
( 3213)     END IF
( 3214) 
( 3215)     ! Sync up processors before returning from subroutine
( 3216)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3217)     IF(mpistat /= 0) THEN
( 3218)        status = NCDF_MPERR
( 3219)        RETURN
( 3220)     END IF
( 3221) 
( 3222)   END SUBROUTINE ncdf_create_var
( 3223) 
( 3224)   ! ncdf_put_att_int creates a new attribute in an existing dataset
( 3225)   ! filename - name of the file to modified
( 3226)   ! varname - name of the variable to attach the attribute to (or "GLOBAL")
( 3227)   ! attname - name of the new attribute
( 3228)   ! attval - value the attribute should have
( 3229)   ! status - return status
( 3230)   SUBROUTINE ncdf_put_att_int(filename, varname, attname, attval, status)
( 3231)     IMPLICIT NONE
( 3232) 
( 3233)     ! Subroutine argument declarations
( 3234)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname
( 3235)     INTEGER,INTENT(IN) :: attval
( 3236)     INTEGER,INTENT(OUT) :: status






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 57

( 3237) 
( 3238)     ! Local declarations
( 3239)     INTEGER :: ncid,    &  ! file ID
( 3240)                nfstat,  &  ! netCDF call return status
( 3241)                mpistat, &  ! MPI call return status
( 3242)                varid       ! ID of variable to attach attribute to
( 3243) 
( 3244)     ! Initialization
( 3245)     status = NCDF_NOERR
( 3246) 
( 3247)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3248)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_int - Creating attribute in file:', filename
( 3249)        CALL FLUSH
( 3250)     END IF
( 3251) 
( 3252)     IF(nproc == 0) THEN
( 3253)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3254)        IF(nfstat /= nf90_noerr) THEN
( 3255)           status = NCDF_NFERR
( 3256)           RETURN
( 3257)        END IF
( 3258)        ! Get variable ID
( 3259)        IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN
( 3260)           varid = NF90_GLOBAL
( 3261)        ELSE
( 3262)           nfstat = nf90_inq_varid(ncid, varname, varid)
( 3263)           IF(nfstat /= nf90_noerr) THEN
( 3264)              status = NCDF_NFERR
( 3265)              RETURN
( 3266)           END IF
( 3267)        END IF
( 3268)        ! Write attribute
( 3269)        nfstat = nf90_redef(ncid)
( 3270)        IF(nfstat /= nf90_noerr) THEN
( 3271)           status = NCDF_NFERR
( 3272)           RETURN
( 3273)        END IF
( 3274)        nfstat = nf90_put_att(ncid, varid, attname, attval)
( 3275)        IF(nfstat /= nf90_noerr) THEN
( 3276)           status = NCDF_NFERR
( 3277)           RETURN
( 3278)        END IF
( 3279)        nfstat = nf90_close(ncid)
( 3280)        IF(nfstat /= nf90_noerr) THEN
( 3281)           status = NCDF_NFERR
( 3282)           RETURN
( 3283)        END IF
( 3284)     END IF
( 3285) 
( 3286)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3287)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_int - Done creating attribute in file:', filename
( 3288)        CALL FLUSH
( 3289)     END IF
( 3290) 
( 3291)     ! Sync up processors before returning from subroutine
( 3292)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3293)     IF(mpistat /= 0) THEN
( 3294)        status = NCDF_MPERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 58

( 3295)        RETURN
( 3296)     END IF
( 3297) 
( 3298)   END SUBROUTINE ncdf_put_att_int
( 3299) 
( 3300)   ! ncdf_put_att_real creates a new attribute in an existing dataset
( 3301)   ! filename - name of the file to modified
( 3302)   ! varname - name of the variable to attach the attribute to (or "GLOBAL")
( 3303)   ! attname - name of the new attribute
( 3304)   ! attval - value the attribute should have
( 3305)   ! status - return status
( 3306)   SUBROUTINE ncdf_put_att_real(filename, varname, attname, attval, status)
( 3307)     IMPLICIT NONE
( 3308) 
( 3309)     ! Subroutine argument declarations
( 3310)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname
( 3311)     REAL(wp),INTENT(IN) :: attval
( 3312)     INTEGER,INTENT(OUT) :: status
( 3313) 
( 3314)     ! Local declarations
( 3315)     INTEGER :: ncid,    &  ! file ID
( 3316)                nfstat,  &  ! netCDF call return status
( 3317)                mpistat, &  ! MPI call return status
( 3318)                varid       ! ID of variable to attach attribute to
( 3319) 
( 3320)     ! Initialization
( 3321)     status = NCDF_NOERR
( 3322) 
( 3323)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3324)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_real - Creating attribute in file:', filename
( 3325)        CALL FLUSH
( 3326)     END IF
( 3327) 
( 3328)     IF(nproc == 0) THEN
( 3329)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3330)        IF(nfstat /= nf90_noerr) THEN
( 3331)           status = NCDF_NFERR
( 3332)           RETURN
( 3333)        END IF
( 3334)        ! Get variable ID
( 3335)        IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN
( 3336)           varid = NF90_GLOBAL
( 3337)        ELSE
( 3338)           nfstat = nf90_inq_varid(ncid, varname, varid)
( 3339)           IF(nfstat /= nf90_noerr) THEN
( 3340)              status = NCDF_NFERR
( 3341)              RETURN
( 3342)           END IF
( 3343)        END IF
( 3344)        ! Write attribute
( 3345)        nfstat = nf90_redef(ncid)
( 3346)        IF(nfstat /= nf90_noerr) THEN
( 3347)           status = NCDF_NFERR
( 3348)           RETURN
( 3349)        END IF
( 3350)        nfstat = nf90_put_att(ncid, varid, attname, attval)
( 3351)        IF(nfstat /= nf90_noerr) THEN
( 3352)           status = NCDF_NFERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 59

( 3353)           RETURN
( 3354)        END IF
( 3355)        nfstat = nf90_close(ncid)
( 3356)        IF(nfstat /= nf90_noerr) THEN
( 3357)           status = NCDF_NFERR
( 3358)           RETURN
( 3359)        END IF
( 3360)     END IF
( 3361) 
( 3362)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3363)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_real - Done creating attribute in file:', filename
( 3364)        CALL FLUSH
( 3365)     END IF
( 3366) 
( 3367)     ! Sync up processors before returning from subroutine
( 3368)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3369)     IF(mpistat /= 0) THEN
( 3370)        status = NCDF_MPERR
( 3371)        RETURN
( 3372)     END IF
( 3373) 
( 3374)   END SUBROUTINE ncdf_put_att_real
( 3375) 
( 3376)   ! ncdf_put_att_char creates a new attribute in an existing dataset
( 3377)   ! filename - name of the file to modified
( 3378)   ! varname - name of the variable to attach the attribute to (or "GLOBAL")
( 3379)   ! attname - name of the new attribute
( 3380)   ! attval - value the attribute should have
( 3381)   ! status - return status
( 3382)   SUBROUTINE ncdf_put_att_char(filename, varname, attname, attval, status)
( 3383)     IMPLICIT NONE
( 3384) 
( 3385)     ! Subroutine argument declarations
( 3386)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname, attval
( 3387)     INTEGER,INTENT(OUT) :: status
( 3388) 
( 3389)     ! Local declarations
( 3390)     INTEGER :: ncid,    &  ! file ID
( 3391)                nfstat,  &  ! netCDF call return status
( 3392)                mpistat, &  ! MPI call return status
( 3393)                varid       ! ID of variable to attach attribute to
( 3394) 
( 3395)     ! Initialization
( 3396)     status = NCDF_NOERR
( 3397) 
( 3398)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3399)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_char - Creating attribute in file:', filename
( 3400)        CALL FLUSH
( 3401)     END IF
( 3402) 
( 3403)     IF(nproc == 0) THEN
( 3404)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3405)        IF(nfstat /= nf90_noerr) THEN
( 3406)           status = NCDF_NFERR
( 3407)           RETURN
( 3408)        END IF
( 3409)        ! Get variable ID
( 3410)        IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 60

( 3411)           varid = NF90_GLOBAL
( 3412)        ELSE
( 3413)           nfstat = nf90_inq_varid(ncid, varname, varid)
( 3414)           IF(nfstat /= nf90_noerr) THEN
( 3415)              status = NCDF_NFERR
( 3416)              RETURN
( 3417)           END IF
( 3418)        END IF
( 3419)        ! Write attribute
( 3420)        nfstat = nf90_redef(ncid)
( 3421)        IF(nfstat /= nf90_noerr) THEN
( 3422)           status = NCDF_NFERR
( 3423)           RETURN
( 3424)        END IF
( 3425)        nfstat = nf90_put_att(ncid, varid, attname, attval)
( 3426)        IF(nfstat /= nf90_noerr) THEN
( 3427)           status = NCDF_NFERR
( 3428)           RETURN
( 3429)        END IF
( 3430)        nfstat = nf90_close(ncid)
( 3431)        IF(nfstat /= nf90_noerr) THEN
( 3432)           status = NCDF_NFERR
( 3433)           RETURN
( 3434)        END IF
( 3435)     END IF
( 3436) 
( 3437)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3438)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_char - Done creating attribute in file:', filename
( 3439)        CALL FLUSH
( 3440)     END IF
( 3441) 
( 3442)     ! Sync up processors before returning from subroutine
( 3443)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3444)     IF(mpistat /= 0) THEN
( 3445)        status = NCDF_MPERR
( 3446)        RETURN
( 3447)     END IF
( 3448) 
( 3449)   END SUBROUTINE ncdf_put_att_char
( 3450) 
( 3451)   ! ncdf_errstr takes an error code returned by a subroutine in this module
( 3452)   ! and returns a string indicating the meaning of the error code
( 3453)   SUBROUTINE ncdf_errstr(errcode, errstr)
( 3454)     INTEGER,INTENT(IN) :: errcode
( 3455)     CHARACTER(LEN=80),INTENT(OUT) :: errstr
( 3456)     
( 3457)     errstr = 'No error'
( 3458)     
( 3459)     SELECT CASE (errcode)
( 3460)     CASE (0)
( 3461)        errstr = 'No error'
( 3462)     CASE(1)
( 3463)        errstr = 'netCDF error'
( 3464)     CASE(2)
( 3465)        errstr = 'MPI error'
( 3466)     CASE(3)
( 3467)        errstr = 'Invalid arguments to subroutine'
( 3468)     CASE(4)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 61

( 3469)        errstr = 'Other/unrecognized error'
( 3470)     CASE DEFAULT
( 3471)        errstr = 'Unknown error code'
( 3472)     END SELECT
( 3473) 
( 3474)   END SUBROUTINE ncdf_errstr
( 3475) 
( 3476) 
( 3477) !!DB: 2007.12.11
( 3478)   ! ncdf_create_file_aveTSUV builds a special OPA output file with most of the default
( 3479)   ! dimensions, variables and attributes for TSUV fields
( 3480)   ! It is designed for average outputs of useful variables:
( 3481)   ! u, v, T, S, taux, tauy, ssh, ...
( 3482) !!2008.06.13 w added ----> 1 more dimid; 2 more varids
( 3483) !float vovecrtz (time_counter, depthw, y, x)
( 3484) !                char units = "m/s"
( 3485) !                float missing_value = 1.000000e+20
( 3486) !                float valid_min = 1.000000e+20
( 3487) !                float valid_max = -1.000000e+20
( 3488) !                char long_name = "Vertical Velocity"
( 3489) !                char short_name = "vovecrtz"
( 3490) !                char online_operation = "ave(x)"
( 3491) !                char axis = "TZYX"
( 3492) !                float interval_operation = 4.800000e+02
( 3493) !                float interval_write = 8.640000e+04
( 3494) !                char associate = "time_counter depthw nav_lat nav_lon"
( 3495) 
( 3496) 
( 3497)   SUBROUTINE ncdf_create_file_aveTSUV(filename, op_type, status)
( 3498)     IMPLICIT NONE
( 3499)     ! Subroutine argument declarations
( 3500)     CHARACTER(LEN=*),INTENT(IN) :: filename
( 3501)     CHARACTER(LEN=*),INTENT(IN) :: op_type
( 3502)     INTEGER,INTENT(OUT) :: status
( 3503) 
( 3504)     ! Local declarations
( 3505)     INTEGER :: ncid,    &  ! netCDF file ID
( 3506)                varid,   &  ! ID of netCDF variable to be written to
( 3507)                nfstat,  &  ! netCDF library call return status
( 3508)                mpistat     ! MPI library call return status
( 3509) !    INTEGER,DIMENSION(1:4) :: dimids
( 3510)     INTEGER,DIMENSION(1:5) :: dimids
( 3511) !    INTEGER,DIMENSION(1:12) :: varids
( 3512) !AD:    INTEGER,DIMENSION(1:14) :: varids
( 3513)     INTEGER,DIMENSION(1:17) :: varids    
( 3514)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
( 3515)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
( 3516)     CHARACTER(LEN=100) :: sec_since    
( 3517)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
( 3518)     INTEGER :: int_opp, &              ! Operation interval
( 3519)                int_wri                 ! Write interval
( 3520)     CHARACTER(LEN=3),PARAMETER :: &
( 3521)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
( 3522)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
( 3523)     
( 3524)     ! Initializations
( 3525)     status = NCDF_NOERR
( 3526)     CALL ioget_calendar(cal_type)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 62

( 3527)     CALL ioget_timestamp(timestamp)
( 3528)     WRITE (UNIT=sec_since, &
( 3529)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 3530)          &  nyear,nmonth,nday,0, 0, 0
( 3531)     WRITE(t_origin, &
( 3532)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
( 3533)          &   nyear,months(nmonth),nday,0,0,0
( 3534)     int_opp = nwrite * rdt
( 3535)     int_wri = nwrite * rdt
( 3536) 
( 3537)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3538)        WRITE(100,*) 'NCDF DEBUG: Creating default U output file:', filename
( 3539)        CALL FLUSH
( 3540)     END IF
( 3541) 
( 3542)     ! Only processor 0 does anything
( 3543)     IF(nproc == 0) THEN
( 3544)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3545)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Creating file:', filename
( 3546)           CALL FLUSH
( 3547)        END IF
( 3548)        ! Create the file
( 3549)        nfstat = nf90_create(filename, nf90_clobber, ncid)
( 3550)        IF(nfstat /= nf90_noerr) THEN
( 3551)           status = NCDF_NFERR
( 3552)           RETURN
( 3553)        END IF
( 3554)        
( 3555)        ! Define dimensions
( 3556)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3557)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining dimensions in file:', filename
( 3558)           CALL FLUSH
( 3559)        END IF
( 3560)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
( 3561)        nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
( 3562)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
( 3563)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
( 3564)        nfstat = nf90_def_dim(ncid, 'depthw', jpkdta, dimids(5))
( 3565)        
( 3566)        ! Define variables
( 3567)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3568)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining variables in file:', filename
( 3569)           CALL FLUSH
( 3570)        END IF
( 3571)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
( 3572)             (/ dimids(4), dimids(3) /), &
( 3573)             varids(1))
( 3574)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
( 3575)             (/ dimids(4), dimids(3) /), &
( 3576)             varids(2))
( 3577)        nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
( 3578)             (/ dimids(2) /), &
( 3579)             varids(3))
( 3580)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
( 3581)             (/ dimids(1) /), &
( 3582)             varids(4))
( 3583)        nfstat = nf90_def_var(ncid, 'vozocrtx', nf90_float, &
( 3584)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 63

( 3585)             varids(5))
( 3586)        nfstat = nf90_def_var(ncid, 'sozotaux', nf90_float, &
( 3587)             (/ dimids(4), dimids(3), dimids(1) /), &
( 3588)             varids(6))
( 3589) 
( 3590)        nfstat = nf90_def_var(ncid, 'vomecrty', nf90_float, &
( 3591)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
( 3592)             varids(7))
( 3593)        nfstat = nf90_def_var(ncid, 'sometauy', nf90_float, &
( 3594)             (/ dimids(4), dimids(3), dimids(1) /), &
( 3595)             varids(8))
( 3596) 
( 3597)        nfstat = nf90_def_var(ncid, 'votemper', nf90_float, &
( 3598)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
( 3599)             varids(9))
( 3600)        nfstat = nf90_def_var(ncid, 'vosaline', nf90_float, &
( 3601)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
( 3602)             varids(10))
( 3603)        nfstat = nf90_def_var(ncid, 'sossheig', nf90_float, &
( 3604)             (/ dimids(4), dimids(3), dimids(1) /), &
( 3605)             varids(11))
( 3606)         nfstat = nf90_def_var(ncid, 'soicecov', nf90_float, &
( 3607)             (/ dimids(4), dimids(3), dimids(1) /), &
( 3608)             varids(12))
( 3609) 
( 3610)        nfstat = nf90_def_var(ncid, 'depthw', nf90_float, &
( 3611)             (/ dimids(5) /), &
( 3612)             varids(13))
( 3613)        nfstat = nf90_def_var(ncid, 'vovecrtz', nf90_float, &
( 3614)             (/ dimids(4), dimids(3), dimids(5), dimids(1) /), &
( 3615)             varids(14))
( 3616) !AD: add 
( 3617)        ! ndate (ndastp)
( 3618)        nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
( 3619)             (/ dimids(1) /), &
( 3620)             varids(15))
( 3621) 
( 3622)        ! ndate (model_time)
( 3623)        nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
( 3624)             (/ dimids(1) /), &
( 3625)             varids(16))
( 3626) 
( 3627)        ! kt 
( 3628)        nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
( 3629)             (/ dimids(1) /), &
( 3630)             varids(17))
( 3631) 
( 3632) 
( 3633)        
( 3634)        ! Add attributes
( 3635)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3636)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_aveTSUV - Writing attributes in file:', filename
( 3637)           CALL FLUSH
( 3638)        END IF
( 3639)        ! nav_lon
( 3640)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
( 3641)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
( 3642)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 64

( 3643)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
( 3644)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
( 3645) 
( 3646)        ! nav_lat
( 3647)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
( 3648)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
( 3649)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
( 3650)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
( 3651)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
( 3652) 
( 3653)        ! depthu
( 3654)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
( 3655)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
( 3656)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
( 3657)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
( 3658)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
( 3659)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')
( 3660) 
( 3661)        ! depthw
( 3662)        nfstat = nf90_put_att(ncid, varids(13), 'units', 'm')
( 3663)        nfstat = nf90_put_att(ncid, varids(13), 'positive', 'unknown')
( 3664)        nfstat = nf90_put_att(ncid, varids(13), 'valid_min', 0.00E00)
( 3665)        nfstat = nf90_put_att(ncid, varids(13), 'valid_max', 5.875141E03)
( 3666)        nfstat = nf90_put_att(ncid, varids(13), 'title', 'depthw')
( 3667)        nfstat = nf90_put_att(ncid, varids(13), 'long_name', 'Vertical W levels')
( 3668) 
( 3669)        ! time_counter
( 3670)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
( 3671)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
( 3672)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
( 3673)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
( 3674)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
( 3675) 
( 3676)        ! vozocrtx
( 3677)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
( 3678)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
( 3679)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
( 3680)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
( 3681)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Zonal Current')
( 3682)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vozocrtx')
( 3683)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
( 3684)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
( 3685)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
( 3686)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
( 3687)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')
( 3688) 
( 3689)        ! sozotaux
( 3690)        nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
( 3691)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
( 3692)        nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
( 3693)        nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
( 3694)        nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along i-axis')
( 3695)        nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sozotaux')
( 3696)        nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
( 3697)        nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')
( 3698)        nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
( 3699)        nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
( 3700)        nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 65

( 3701) 
( 3702)        ! vomecrty
( 3703)        nfstat = nf90_put_att(ncid, varids(7), 'units', 'm/s')
( 3704)        nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
( 3705)        nfstat = nf90_put_att(ncid, varids(7), 'valid_min', 1.000000E20 )
( 3706)        nfstat = nf90_put_att(ncid, varids(7), 'valid_max', -1.000000E20)
( 3707)        nfstat = nf90_put_att(ncid, varids(7), 'long_name', 'Meridional Current')
( 3708)        nfstat = nf90_put_att(ncid, varids(7), 'short_name', 'vomecrty')
( 3709)        nfstat = nf90_put_att(ncid, varids(7), 'online_operation', TRIM(op_type))
( 3710)        nfstat = nf90_put_att(ncid, varids(7), 'axis', 'TZYX')
( 3711)        nfstat = nf90_put_att(ncid, varids(7), 'interval_operation', int_opp)
( 3712)        nfstat = nf90_put_att(ncid, varids(7), 'interval_write', int_wri)
( 3713)        nfstat = nf90_put_att(ncid, varids(7), 'associate', 'time_counter depthv nav_lat nav_lon')
( 3714) 
( 3715)        ! sometauy
( 3716)        nfstat = nf90_put_att(ncid, varids(8), 'units', 'N/m2')
( 3717)        nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
( 3718)        nfstat = nf90_put_att(ncid, varids(8), 'valid_min', 1.000000E20 )
( 3719)        nfstat = nf90_put_att(ncid, varids(8), 'valid_max', -1.000000E20)
( 3720)        nfstat = nf90_put_att(ncid, varids(8), 'long_name', 'Wind Stress along j-axis')
( 3721)        nfstat = nf90_put_att(ncid, varids(8), 'short_name', 'sometauy')
( 3722)        nfstat = nf90_put_att(ncid, varids(8), 'online_operation', TRIM(op_type))
( 3723)        nfstat = nf90_put_att(ncid, varids(8), 'axis', 'TYX')
( 3724)        nfstat = nf90_put_att(ncid, varids(8), 'interval_operation', int_opp)
( 3725)        nfstat = nf90_put_att(ncid, varids(8), 'interval_write', int_wri)
( 3726)        nfstat = nf90_put_att(ncid, varids(8), 'associate', 'time_counter nav_lat nav_lon')
( 3727) 
( 3728)        ! votemper
( 3729)        nfstat = nf90_put_att(ncid, varids(9), 'units', 'C')
( 3730)        nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)
( 3731)        nfstat = nf90_put_att(ncid, varids(9), 'valid_min', 1.000000E20 )
( 3732)        nfstat = nf90_put_att(ncid, varids(9), 'valid_max', -1.000000E20)
( 3733)        nfstat = nf90_put_att(ncid, varids(9), 'long_name', 'Temperature')
( 3734)        nfstat = nf90_put_att(ncid, varids(9), 'short_name', 'votemper')
( 3735)        nfstat = nf90_put_att(ncid, varids(9), 'online_operation', TRIM(op_type))
( 3736)        nfstat = nf90_put_att(ncid, varids(9), 'axis', 'TZYX')
( 3737)        nfstat = nf90_put_att(ncid, varids(9), 'interval_operation', int_opp)
( 3738)        nfstat = nf90_put_att(ncid, varids(9), 'interval_write', int_wri)
( 3739)        nfstat = nf90_put_att(ncid, varids(9), 'associate', 'time_counter deptht nav_lat nav_lon')
( 3740) 
( 3741)        ! vosaline
( 3742)        nfstat = nf90_put_att(ncid, varids(10), 'units', 'PSU')
( 3743)        nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)
( 3744)        nfstat = nf90_put_att(ncid, varids(10), 'valid_min', 1.000000E20 )
( 3745)        nfstat = nf90_put_att(ncid, varids(10), 'long_name', 'Salinity')
( 3746)        nfstat = nf90_put_att(ncid, varids(10), 'short_name', 'vosaline')
( 3747)        nfstat = nf90_put_att(ncid, varids(10), 'online_operation', TRIM(op_type))
( 3748)        nfstat = nf90_put_att(ncid, varids(10), 'axis', 'TZYX')
( 3749)        nfstat = nf90_put_att(ncid, varids(10), 'interval_operation', int_opp)
( 3750)        nfstat = nf90_put_att(ncid, varids(10), 'interval_write', int_wri)
( 3751)        nfstat = nf90_put_att(ncid, varids(10), 'associate', 'time_counter deptht nav_lat nav_lon')
( 3752) 
( 3753)        ! sossheig
( 3754)        nfstat = nf90_put_att(ncid, varids(11), 'units', 'm')
( 3755)        nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)
( 3756)        nfstat = nf90_put_att(ncid, varids(11), 'valid_min', 1.000000E20 )
( 3757)        nfstat = nf90_put_att(ncid, varids(11), 'valid_max', -1.000000E20)
( 3758)        nfstat = nf90_put_att(ncid, varids(11), 'long_name', 'Sea Surface Height')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 66

( 3759)        nfstat = nf90_put_att(ncid, varids(11), 'short_name', 'sossheig')
( 3760)        nfstat = nf90_put_att(ncid, varids(11), 'online_operation', TRIM(op_type))
( 3761)        nfstat = nf90_put_att(ncid, varids(11), 'axis', 'TYX')
( 3762)        nfstat = nf90_put_att(ncid, varids(11), 'interval_operation', int_opp)
( 3763)        nfstat = nf90_put_att(ncid, varids(11), 'interval_write', int_wri)
( 3764)        nfstat = nf90_put_att(ncid, varids(11), 'associate', 'time_counter nav_lat nav_lon')
( 3765) 
( 3766)        ! soicecov
( 3767)        nfstat = nf90_put_att(ncid, varids(12), 'units', '[0,1]')
( 3768)        nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)
( 3769)        nfstat = nf90_put_att(ncid, varids(12), 'valid_min', 1.000000E20 )
( 3770)        nfstat = nf90_put_att(ncid, varids(12), 'valid_max', -1.000000E20)
( 3771)        nfstat = nf90_put_att(ncid, varids(12), 'long_name', 'Ice Cover')
( 3772)        nfstat = nf90_put_att(ncid, varids(12), 'short_name', 'soicecov')
( 3773)        nfstat = nf90_put_att(ncid, varids(12), 'online_operation', TRIM(op_type))
( 3774)        nfstat = nf90_put_att(ncid, varids(12), 'axis', 'TYX')
( 3775)        nfstat = nf90_put_att(ncid, varids(12), 'interval_operation', int_opp)
( 3776)        nfstat = nf90_put_att(ncid, varids(12), 'interval_write', int_wri)
( 3777)        nfstat = nf90_put_att(ncid, varids(12), 'associate', 'time_counter nav_lat nav_lon')
( 3778) 
( 3779)        ! vovecrtz
( 3780)        nfstat = nf90_put_att(ncid, varids(14), 'units', 'm/s')
( 3781)        nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)
( 3782)        nfstat = nf90_put_att(ncid, varids(14), 'valid_min', 1.000000E20 )
( 3783)        nfstat = nf90_put_att(ncid, varids(14), 'valid_max', -1.000000E20)
( 3784)        nfstat = nf90_put_att(ncid, varids(14), 'long_name', 'Vertical Velocity')
( 3785)        nfstat = nf90_put_att(ncid, varids(14), 'short_name', 'vovecrtz')
( 3786)        nfstat = nf90_put_att(ncid, varids(14), 'online_operation', TRIM(op_type))
( 3787)        nfstat = nf90_put_att(ncid, varids(14), 'axis', 'TZYX')
( 3788)        nfstat = nf90_put_att(ncid, varids(14), 'interval_operation', int_opp)
( 3789)        nfstat = nf90_put_att(ncid, varids(14), 'interval_write', int_wri)
( 3790)        nfstat = nf90_put_att(ncid, varids(14), 'associate', 'time_counter depthw nav_lat nav_lon')
( 3791) 
( 3792)        ! ndate (ndastp)
( 3793)        nfstat = nf90_put_att(ncid, varids(15), 'units', '=nyear*10000+nmonth*100+nday')
( 3794)        nfstat = nf90_put_att(ncid, varids(15), 'long_name', 'time step date in year/month/day aammjj')
( 3795) 
( 3796)        ! ndate (model_time)
( 3797)        nfstat = nf90_put_att(ncid, varids(16), 'long_name', &
( 3798)             'time step date (when output is writen) in year/month/day aammjj (decimal day)')
( 3799)        nfstat = nf90_put_att(ncid, varids(16), 'units', '=nyear*10000+nmonth*100+nday')
( 3800)        nfstat = nf90_put_att(ncid, varids(16), 'formula1', 'nyear  =   model_time / 10000')       
( 3801)        nfstat = nf90_put_att(ncid, varids(16), 'formula2', & 
( 3802)             'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
( 3803)        nfstat = nf90_put_att(ncid, varids(16), 'formula3', & 
( 3804)             'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           
( 3805) 
( 3806) 
( 3807)        ! global
( 3808)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
( 3809)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
( 3810)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
( 3811)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
( 3812)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
( 3813)       
( 3814)        
( 3815)        ! Close file
( 3816)        IF(DEBUG_OUT .EQV. .TRUE.) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 67

( 3817)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_TSUV - Closing file:', filename
( 3818)           CALL FLUSH
( 3819)        END IF
( 3820)        nfstat = nf90_close(ncid)
( 3821)        IF(nfstat /= nf90_noerr) THEN
( 3822)           status = NCDF_NFERR
( 3823)           RETURN
( 3824)        END IF
( 3825)     END IF
( 3826) 
( 3827)     ! Sync up processors before returning from subroutine
( 3828)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3829)     IF(mpistat /= 0) THEN
( 3830)        status = NCDF_MPERR
( 3831)        RETURN
( 3832)     END IF
( 3833)     
( 3834)   END SUBROUTINE ncdf_create_file_aveTSUV
( 3835) 
( 3836) !!DB 2008.05.22
( 3837)   ! ncdf_create_ice_restart builds a standard OPA lim ice restart file with all the default
( 3838)   ! dimensions, variables and attributes. Note that this is a single restart
( 3839)   ! file, and therefore incompatible with the default IOIPSL-based restart
( 3840)   ! routines (though it holds the same data)
( 3841)   SUBROUTINE ncdf_create_ice_restart(filename, status)
( 3842)     IMPLICIT NONE
( 3843)     ! Subroutine argument declarations
( 3844)     CHARACTER(LEN=*),INTENT(IN) :: filename
( 3845)     INTEGER,INTENT(OUT) :: status
( 3846)     
( 3847)     ! Local declarations
( 3848)     INTEGER :: ncid,    &  ! netCDF file ID
( 3849)          varid,   &  ! ID of netCDF variable to be written to
( 3850)          nfstat,  &  ! netCDF library call return status
( 3851)          mpistat     ! MPI library call return status
( 3852)     INTEGER,DIMENSION(1:9) :: dimids
( 3853)     INTEGER,DIMENSION(1:36) :: varids
( 3854)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
( 3855)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
( 3856)     CHARACTER(LEN=100) :: sec_since
( 3857)     CHARACTER(LEN=100) :: tstp_since
( 3858)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
( 3859)     CHARACTER(LEN=3),PARAMETER :: &
( 3860)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
( 3861)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
( 3862)     INTEGER :: ji
( 3863) 
( 3864)     
( 3865)     ! Initializations
( 3866)     status = NCDF_NOERR
( 3867)     
( 3868)     ! Initializations
( 3869)     status = NCDF_NOERR
( 3870)     CALL ioget_calendar(cal_type)
( 3871)     CALL ioget_timestamp(timestamp)
( 3872)     WRITE (UNIT=sec_since, &
( 3873)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 3874)          &  nyear,nmonth,nday,0, 0, 0






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 68

( 3875)     WRITE (UNIT=tstp_since, &
( 3876)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 3877)          &  nyear,nmonth,nday,0, 0, 0
( 3878)     WRITE(t_origin, &
( 3879)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
( 3880)          &   nyear,months(nmonth),nday,0,0,0
( 3881) 
( 3882)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3883)        WRITE(100,*) 'NCDF DEBUG: Creating default restart file:', filename
( 3884)        CALL FLUSH
( 3885)     END IF
( 3886)     
( 3887)     ! Only processor 0 does anything
( 3888)     IF(nproc == 0) THEN
( 3889)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3890)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Creating file:', filename
( 3891)           CALL FLUSH
( 3892)        END IF
( 3893)        ! Create the file
( 3894)        nfstat = nf90_create(filename, nf90_clobber, ncid)
( 3895)        
( 3896)        ! Define dimensions
( 3897)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3898)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining dimensions in file:', filename
( 3899)           CALL FLUSH
( 3900)        END IF
( 3901)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(1))
( 3902)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(2))
( 3903)        nfstat = nf90_def_dim(ncid, 'z', 1, dimids(3))
( 3904)        nfstat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimids(4))
( 3905)        nfstat = nf90_def_dim(ncid, 'x_a', 1, dimids(5))
( 3906)        nfstat = nf90_def_dim(ncid, 'y_a', 1, dimids(6))
( 3907)        nfstat = nf90_def_dim(ncid, 'z_a', 2, dimids(7))
( 3908)        nfstat = nf90_def_dim(ncid, 'z_b', 3, dimids(8))
( 3909)        nfstat = nf90_def_dim(ncid, 'z_c', 35, dimids(9))
( 3910)        
( 3911)        ! Define variables
( 3912)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3913)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining variables in file:', filename
( 3914)           CALL FLUSH
( 3915)        END IF
( 3916)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
( 3917)             (/ dimids(1), dimids(2) /), &
( 3918)             varids(1))
( 3919)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
( 3920)             (/ dimids(1), dimids(2) /), &
( 3921)             varids(2))
( 3922)        nfstat = nf90_def_var(ncid, 'nav_lev', nf90_float, &
( 3923)             (/ dimids(3) /), &
( 3924)             varids(3))
( 3925)        nfstat = nf90_def_var(ncid, 'time', nf90_float, &
( 3926)             (/ dimids(4) /), &
( 3927)             varids(4))
( 3928)        nfstat = nf90_def_var(ncid, 'time_steps', nf90_int, &
( 3929)             (/ dimids(4) /), &
( 3930)             varids(5))
( 3931)        nfstat = nf90_def_var(ncid, 'info', nf90_double, &
( 3932)             (/ dimids(7) /), &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 69

( 3933)             varids(6))
( 3934) 
( 3935)        nfstat = nf90_def_var(ncid, 'hicif'  , nf90_double, &
( 3936)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3937)             varids(7))
( 3938)        nfstat = nf90_def_var(ncid, 'hsnif'   , nf90_double, &
( 3939)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3940)             varids(8))
( 3941)        nfstat = nf90_def_var(ncid, 'frld'    , nf90_double, &
( 3942)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3943)             varids(9))
( 3944)        nfstat = nf90_def_var(ncid, 'sist'    , nf90_double, &
( 3945)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3946)             varids(10))
( 3947)        nfstat = nf90_def_var(ncid, 'tbif'   ,  nf90_double, &
( 3948)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 3949)             varids(11))
( 3950)        nfstat = nf90_def_var(ncid, 'u_ice'  ,  nf90_double, &
( 3951)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3952)             varids(12))
( 3953)        nfstat = nf90_def_var(ncid, 'v_ice'  ,  nf90_double, &
( 3954)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3955)             varids(13))
( 3956)        nfstat = nf90_def_var(ncid, 'gtaux'  ,  nf90_double, &
( 3957)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3958)             varids(14))
( 3959)        nfstat = nf90_def_var(ncid, 'gtauy'  ,  nf90_double, &
( 3960)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3961)             varids(15))
( 3962)        nfstat = nf90_def_var(ncid, 'qstoif' ,  nf90_double, &
( 3963)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3964)             varids(16))
( 3965)        nfstat = nf90_def_var(ncid, 'fsbbq'  ,  nf90_double, &
( 3966)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3967)             varids(17))
( 3968)        nfstat = nf90_def_var(ncid, 'moment' ,  nf90_double, &
( 3969)             (/ dimids(1), dimids(2), dimids(9), dimids(4) /), &
( 3970)             varids(18))
( 3971) # if defined key_coupled
( 3972)       nfstat = nf90_def_var(ncid, 'albege' , , nf90_double, &
( 3973)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3974)             varids(19))
( 3975) # endif
( 3976) 
( 3977)        
( 3978)        ! Add attributes
( 3979)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3980)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Writing attributes in file:', filename
( 3981)           CALL FLUSH
( 3982)        END IF
( 3983)        ! nav_lon
( 3984)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
( 3985)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -1.800000E2)
( 3986)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', 1.800000E2)
( 3987)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
( 3988) 
( 3989)        ! nav_lat
( 3990)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 70

( 3991)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', -9.000000E1)
( 3992)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 9.000000E1)
( 3993)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
( 3994) 
( 3995)        ! nav_lev
( 3996)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'model_levels')
( 3997)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 0.0)
( 3998)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 0.0)
( 3999)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Model levels')
( 4000) 
( 4001)        ! time
( 4002)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
( 4003)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
( 4004)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
( 4005)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'Time axis')
( 4006)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', '0001-JUL-01 00:00:00')
( 4007) 
( 4008)        ! time_steps
( 4009)        nfstat = nf90_put_att(ncid, varids(5), 'units', TRIM(tstp_since))
( 4010)        nfstat = nf90_put_att(ncid, varids(5), 'title', 'Time steps')
( 4011)        nfstat = nf90_put_att(ncid, varids(5), 'tstep_sec', 1.877760E7)
( 4012)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Time step axis')
( 4013)        nfstat = nf90_put_att(ncid, varids(5), 'time_origin', TRIM(t_origin))
( 4014) 
( 4015)        ! info
( 4016)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
( 4017) !DB as per OLD code
( 4018)        do ji = 7, 18
( 4019)           nfstat = nf90_put_att(ncid, varids(ji), 'missing_value', 1.000000E20)
( 4020)        enddo
( 4021) #if defined key_coupled
( 4022)           nfstat = nf90_put_att(ncid, varids(19), 'missing_value', 1.000000E20)
( 4023) #endif
( 4024) 
( 4025)        ! global
( 4026)        nfstat = nf90_put_att(ncid, NF90_GLOBAL , 'Conventions', 'GDT 1.2')
( 4027)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
( 4028)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
( 4029)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number_total', 0)
( 4030)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number', 0)
( 4031)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids', 0)
( 4032)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_global', 0)
( 4033)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_local', 0)
( 4034)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_first', 0)
( 4035)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_last', 0)
( 4036)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', 0)
( 4037)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', 0)
( 4038)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_type','box')
( 4039) 
( 4040)        ! Close file
( 4041)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 4042)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Closing file:', filename
( 4043)           CALL FLUSH
( 4044)        END IF
( 4045)        nfstat = nf90_close(ncid)
( 4046)        IF(nfstat /= nf90_noerr) THEN
( 4047)           status = NCDF_NFERR
( 4048)           RETURN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 71

( 4049)        END IF
( 4050)     END IF
( 4051) 
( 4052)     ! Sync up processors before returning from subroutine
( 4053)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4054)     IF(mpistat /= 0) THEN
( 4055)        status = NCDF_MPERR
( 4056)        RETURN
( 4057)     END IF
( 4058)        
( 4059)   END SUBROUTINE ncdf_create_ice_restart
( 4060) 
( 4061) !! CN 10/14/2008 - finds number of dimensions in the specified variable
( 4062)   SUBROUTINE ncdf_get_num_dims(filename, varname, ndims, status)
( 4063)     IMPLICIT NONE
( 4064)     ! Subroutine argument declarations
( 4065)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4066)     INTEGER,INTENT(OUT) :: ndims, status
( 4067) 
( 4068)     ! Local declarations
( 4069)     INTEGER :: ncid,    &  ! netCDF file ID
( 4070)                id,      &  ! ID of netCDF variable of interest
( 4071)                nfstat,  &  ! netCDF library call return status
( 4072)                mpistat     ! MPI library call return status
( 4073) 
( 4074) 
( 4075)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4076)        WRITE(100,*) 'NCDF DEBUG: Getting number of dimensions in variable from file:', filename
( 4077)        CALL FLUSH
( 4078)     END IF
( 4079) 
( 4080)     ! Initializations
( 4081)     status = NCDF_NOERR
( 4082) 
( 4083)     ! Open netCDF file and get info
( 4084)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4085)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Opening file:', filename
( 4086)        CALL FLUSH
( 4087)     END IF
( 4088)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4089)     IF(nfstat /= nf90_noerr) THEN
( 4090)        status = NCDF_NFERR
( 4091)        RETURN
( 4092)     END IF
( 4093)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4094)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Getting varid from file:', filename
( 4095)        CALL FLUSH
( 4096)     END IF
( 4097)     nfstat = nf90_inq_varid(ncid, varname, id)
( 4098)     IF(nfstat /= nf90_noerr) THEN
( 4099)        status = NCDF_NFERR
( 4100)        RETURN
( 4101)     END IF
( 4102) 
( 4103)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4104)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Getting number of dims in variable from file:', filename
( 4105)        CALL FLUSH
( 4106)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 72

( 4107)     nfstat = nf90_inquire_variable(ncid, id, ndims=ndims)
( 4108)     IF(nfstat /= nf90_noerr) THEN
( 4109)        status = NCDF_NFERR
( 4110)        RETURN
( 4111)     END IF
( 4112) 
( 4113)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4114)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Closing file:', filename
( 4115)        CALL FLUSH
( 4116)     END IF
( 4117)     nfstat = nf90_close(ncid)
( 4118)     IF(nfstat /= nf90_noerr) THEN
( 4119)        status = NCDF_NFERR
( 4120)        RETURN
( 4121)     END IF
( 4122)  
( 4123)     ! Sync up processors before returning from subroutine
( 4124)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4125)     IF(mpistat /= 0) THEN
( 4126)        status = NCDF_MPERR
( 4127)        RETURN
( 4128)     END IF
( 4129) 
( 4130)   END SUBROUTINE ncdf_get_num_dims
( 4131) 
( 4132) !!DB
( 4133) !!Get size of dimension -- common use is to get size of unlimited dimension
( 4134) !!(i.e. # records)
( 4135) !!Note that all processors open the file and get the dimension size
( 4136)   SUBROUTINE ncdf_get_dim_size(filename, dimname, len, status)
( 4137)     IMPLICIT NONE
( 4138)     ! Subroutine argument declarations
( 4139)     CHARACTER(LEN=*),INTENT(IN) :: filename, dimname
( 4140)     INTEGER,INTENT(OUT) :: len, status
( 4141) 
( 4142)     ! Local declarations
( 4143)     INTEGER :: ncid,    &  ! netCDF file ID
( 4144)                id,      &  ! ID of netCDF dimension of interest
( 4145)                nfstat,  &  ! netCDF library call return status
( 4146)                mpistat     ! MPI library call return status
( 4147)     CHARACTER(LEN=nf90_max_name) dname
( 4148) 
( 4149) 
( 4150)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4151)        WRITE(100,*) 'NCDF DEBUG: Getting dim size from file:', filename
( 4152)        CALL FLUSH
( 4153)     END IF
( 4154) 
( 4155)     ! Initializations
( 4156)     status = NCDF_NOERR
( 4157) 
( 4158)     ! Open netCDF file and get info
( 4159)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4160)        WRITE(100,*) 'NCDF DEBUG: ncdf_get dim_size - Opening file:', filename
( 4161)        CALL FLUSH
( 4162)     END IF
( 4163)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4164)     IF(nfstat /= nf90_noerr) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 73

( 4165)        status = NCDF_NFERR
( 4166)        RETURN
( 4167)     END IF
( 4168)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4169)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Getting dimid from file:', filename
( 4170)        CALL FLUSH
( 4171)     END IF
( 4172)     nfstat = NF90_INQ_DIMID (ncid, dimname, id)
( 4173)     IF(nfstat /= nf90_noerr) THEN
( 4174)        status = NCDF_NFERR
( 4175)        RETURN
( 4176)     END IF
( 4177) 
( 4178)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4179)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Getting dimsize from file:', filename
( 4180)        CALL FLUSH
( 4181)     END IF
( 4182)     nfstat = NF90_INQUIRE_DIMENSION (ncid, id, len=len)
( 4183) 
( 4184)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4185)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Closing file:', filename
( 4186)        CALL FLUSH
( 4187)     END IF
( 4188)     nfstat = nf90_close(ncid)
( 4189)     IF(nfstat /= nf90_noerr) THEN
( 4190)        status = NCDF_NFERR
( 4191)        RETURN
( 4192)     END IF
( 4193)  
( 4194)     ! Sync up processors before returning from subroutine
( 4195)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4196)     IF(mpistat /= 0) THEN
( 4197)        status = NCDF_MPERR
( 4198)        RETURN
( 4199)     END IF
( 4200) 
( 4201)   END SUBROUTINE ncdf_get_dim_size
( 4202) 
( 4203)   
( 4204) !>>> A.D:  2008.09.29
( 4205)   
( 4206)     SUBROUTINE ncdf_readdate(filename, varname, cdata, index, status)
( 4207)     IMPLICIT NONE
( 4208)     ! Subroutine argument declarations
( 4209)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4210)     INTEGER,INTENT(IN) :: index
( 4211)     CHARACTER(LEN=19),INTENT(OUT) :: cdata
( 4212)     INTEGER,INTENT(OUT) :: status
( 4213) 
( 4214)     ! Local declarations
( 4215)     INTEGER :: ncid,    &  ! netCDF file ID
( 4216)                varid,   &  ! ID of netCDF variable to be written to
( 4217)                nfstat,  &  ! netCDF library call return status
( 4218)                mpistat     ! MPI library call return status
( 4219)     INTEGER,DIMENSION(1:5) :: count
( 4220) 
( 4221)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4222)        WRITE(100,*) 'NCDF DEBUG: Reading single date from file:', filename






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 74

( 4223)        CALL FLUSH
( 4224)     END IF
( 4225) 
( 4226)     ! Initializations
( 4227)     status = NCDF_NOERR
( 4228)     count(1) = 1
( 4229) 
( 4230)     ! Open netCDF file and get info
( 4231)        
( 4232)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4233)        WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Opening file:', filename
( 4234)        CALL FLUSH
( 4235)     END IF
( 4236)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4237)     IF(nfstat /= nf90_noerr) THEN
( 4238)        status = NCDF_NFERR
( 4239)        RETURN
( 4240)     END IF
( 4241) 
( 4242)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4243)        WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Getting info from file:', filename
( 4244)        CALL FLUSH
( 4245)     END IF
( 4246)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 4247)     IF(nfstat /= nf90_noerr) THEN
( 4248)        status = NCDF_NFERR
( 4249)        RETURN
( 4250)     END IF
( 4251)     
( 4252)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4253)        WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Reading data from file:', filename
( 4254)        CALL FLUSH
( 4255)     END IF
( 4256)     !       nfstat = nf90_get_var(ncid, varid, data, (/ index /))
( 4257)     nfstat = nf90_get_var(ncid, varid, cdata, (/ 1,index /), (/ 19,1 /))           
( 4258)     IF(nfstat /= nf90_noerr) THEN
( 4259)        status = NCDF_NFERR
( 4260)        RETURN
( 4261)     END IF
( 4262) 
( 4263)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4264)        WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Closing file:', filename
( 4265)        CALL FLUSH
( 4266)     END IF
( 4267)     nfstat = nf90_close(ncid)
( 4268)     IF(nfstat /= nf90_noerr) THEN
( 4269)        status = NCDF_NFERR
( 4270)        RETURN
( 4271)     END IF
( 4272)     
( 4273)     ! Sync up processors before returning from subroutine
( 4274)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4275)     IF(mpistat /= 0) THEN
( 4276)        status = NCDF_MPERR
( 4277)        RETURN
( 4278)     END IF
( 4279)     
( 4280)   END SUBROUTINE ncdf_readdate






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 75

( 4281) 
( 4282)   
( 4283) !<<< A.D:  
( 4284)   
( 4285) SUBROUTINE ncdf_read2d_global(filename, varname, data, time, status)
( 4286)     IMPLICIT NONE
( 4287)     ! Subroutine argument declarations
( 4288)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4289)     INTEGER,INTENT(IN) :: time
( 4290)     REAL,DIMENSION(:,:),INTENT(OUT) :: data
( 4291)     INTEGER,INTENT(OUT) :: status
( 4292) 
( 4293)     ! Local declarations
( 4294)     INTEGER :: ncid,    &  ! netCDF file ID
( 4295)                varid,   &  ! ID of netCDF variable to be written to
( 4296)                mpistat, &  ! MPI library call return status
( 4297)                nfstat,  &  ! netCDF library call return status
( 4298)                i,       &  ! Loop counter
( 4299)                j,       &  !     "
( 4300)                k,       &  !     "
( 4301)                ndims       ! Number of dimensions in this variable
( 4302) 
( 4303)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 4304)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 4305)     INTEGER :: tindex
( 4306)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 4307)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 4308)     REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: databuf
( 4309) 
( 4310)     ! Initializations
( 4311)     status = NCDF_NOERR
( 4312)     hastime = .FALSE.
( 4313)     var_dimlens = 0
( 4314)     j = 0
( 4315)     if(time > 0) then
( 4316)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 4317)     else
( 4318)        tindex = -time
( 4319)     endif
( 4320) 
( 4321) 
( 4322)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4323)        WRITE(100,*) 'NCDF DEBUG: Reading 2D array from file:', filename
( 4324)        CALL FLUSH
( 4325)     END IF
( 4326) 
( 4327)     ! Open netCDF file and get info
( 4328)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4329)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Opening file:', filename
( 4330)        CALL FLUSH
( 4331)     END IF
( 4332) 
( 4333)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4334)     IF(nfstat /= nf90_noerr) THEN
( 4335)        status = NCDF_NFERR
( 4336)        RETURN
( 4337)     END IF
( 4338)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 76

( 4339)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Getting info from file:', filename
( 4340)        CALL FLUSH
( 4341)     END IF
( 4342)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 4343)     IF(nfstat /= nf90_noerr) THEN
( 4344)        status = NCDF_NFERR
( 4345)        RETURN
( 4346)     END IF
( 4347)     
( 4348)     ! Determine if this variable contains a time axis
( 4349)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 4350)          ndims=ndims, dimids=var_dimids)
( 4351)     IF(nfstat /= nf90_noerr) THEN
( 4352)        status = NCDF_NFERR
( 4353)        RETURN
( 4354)     END IF
( 4355)     DO i = 1, ndims
( 4356)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 4357)             name=dname, len=k)
( 4358)        dname=TRIM(dname)
( 4359)        var_dimlens(i-j) = k
( 4360)        IF(nfstat /= nf90_noerr) THEN
( 4361)           status = NCDF_NFERR
( 4362)           RETURN
( 4363)        END IF
( 4364)        IF((dname=='time_counter') .OR. (dname=='time')) THEN
( 4365)           hastime = .TRUE.
( 4366)           j = 1
( 4367)           EXIT
( 4368)        END IF
( 4369)     END DO
( 4370)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2)))
( 4371) 
( 4372)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4373)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Reading data from file:', filename
( 4374)        CALL FLUSH
( 4375)     END IF
( 4376)     IF(hastime .EQV. .TRUE.) THEN
( 4377)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4378)             (/ 1, 1, tindex/), &
( 4379)             (/ var_dimlens(1), var_dimlens(2), 1 /))
( 4380)     ELSE
( 4381)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4382)             (/ 1, 1 /), &
( 4383)             (/ var_dimlens(1), var_dimlens(2) /))
( 4384)     END IF
( 4385)     IF(nfstat /= nf90_noerr) THEN
( 4386)        status = NCDF_NFERR
( 4387)        RETURN
( 4388)     END IF
( 4389) 
( 4390)     data = databuf
( 4391) 
( 4392)     ! All done, close up the netCDF file
( 4393)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4394)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Closing file:', filename
( 4395)        CALL FLUSH
( 4396)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 77

( 4397)     nfstat = nf90_close(ncid)
( 4398)     IF(nfstat /= nf90_noerr) THEN
( 4399)        status = NCDF_NFERR
( 4400)        RETURN
( 4401)     END IF
( 4402) 
( 4403)     DEALLOCATE(databuf)
( 4404) 
( 4405)     ! Sync up processors before returning from subroutine
( 4406)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4407)     IF(mpistat /= 0) THEN
( 4408)        status = NCDF_MPERR
( 4409)        RETURN
( 4410)     END IF
( 4411) 
( 4412)   END SUBROUTINE ncdf_read2d_global
( 4413) 
( 4414) 
( 4415) !!DB
( 4416)   ! Reads a 3-D array from a netCDF file.
( 4417)   ! Returns global array to wach processor
( 4418)   ! filename - file to read from
( 4419)   ! varname - variable to read from
( 4420)   ! data - 3-D array to put the data into (must be a REAL)
( 4421)   ! time - time index to read (If the target variable in the netCDF file is
( 4422)   !        actually a 4-D variable with a time axis, ignored otherwise)
( 4423)   ! status - return status of the subroutine
( 4424)   SUBROUTINE ncdf_read3d_global(filename, varname, data, time, status)
( 4425)     IMPLICIT NONE
( 4426)     ! Subroutine argument declarations
( 4427)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4428)     INTEGER,INTENT(IN) :: time
( 4429)     REAL,DIMENSION(:,:,:),INTENT(OUT) :: data
( 4430)     INTEGER,INTENT(OUT) :: status
( 4431) 
( 4432)     ! Local declarations
( 4433)     INTEGER :: ncid,    &  ! netCDF file ID
( 4434)                varid,   &  ! ID of netCDF variable to be written to
( 4435)                mpistat, &  ! MPI library call return status
( 4436)                nfstat,  &  ! netCDF library call return status
( 4437)                i,       &  ! Loop counter
( 4438)                j,       &  !     "
( 4439)                k,       &  !     "
( 4440)                ndims       ! Number of dimensions in this variable
( 4441) 
( 4442)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 4443)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 4444)     INTEGER :: tindex
( 4445)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 4446)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 4447)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: databuf
( 4448) 
( 4449)     ! Initializations
( 4450)     status = NCDF_NOERR
( 4451)     hastime = .FALSE.
( 4452)     var_dimlens = 0
( 4453)     j = 0
( 4454)     if(time > 0) then






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 78

( 4455)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 4456)     else
( 4457)        tindex = -time
( 4458)     endif
( 4459) 
( 4460)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4461)        WRITE(100,*) 'NCDF DEBUG: Reading 3D array from file:', filename
( 4462)        CALL FLUSH
( 4463)     END IF
( 4464) 
( 4465)     ! Open netCDF file and get info
( 4466)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4467)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Opening file:', filename
( 4468)        CALL FLUSH
( 4469)     END IF
( 4470) 
( 4471)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4472)     IF(nfstat /= nf90_noerr) THEN
( 4473)        status = NCDF_NFERR
( 4474)        RETURN
( 4475)     END IF
( 4476)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4477)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Getting info from file:', filename
( 4478)        CALL FLUSH
( 4479)     END IF
( 4480)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 4481)     IF(nfstat /= nf90_noerr) THEN
( 4482)        status = NCDF_NFERR
( 4483)        RETURN
( 4484)     END IF
( 4485)     
( 4486)     ! Determine if this variable contains a time axis
( 4487)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 4488)          ndims=ndims, dimids=var_dimids)
( 4489)     IF(nfstat /= nf90_noerr) THEN
( 4490)        status = NCDF_NFERR
( 4491)        RETURN
( 4492)     END IF
( 4493)     DO i = 1, ndims
( 4494)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 4495)             name=dname, len=k)
( 4496)        dname=TRIM(dname)
( 4497)        var_dimlens(i-j) = k
( 4498)        IF(nfstat /= nf90_noerr) THEN
( 4499)           status = NCDF_NFERR
( 4500)           RETURN
( 4501)        END IF
( 4502)        IF((dname=='time_counter') .OR. (dname=='time')) THEN
( 4503)           hastime = .TRUE.
( 4504)           j = 1
( 4505)           EXIT
( 4506)        END IF
( 4507)     END DO
( 4508)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3)))
( 4509) 
( 4510)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4511)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Reading data from file:', filename
( 4512)        CALL FLUSH






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 79

( 4513)     END IF
( 4514)     IF(hastime .EQV. .TRUE.) THEN
( 4515)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4516)             (/ 1, 1, 1, tindex/), &
( 4517)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), 1 /))
( 4518)     ELSE
( 4519)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4520)             (/ 1, 1, 1 /), &
( 4521)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3) /))
( 4522)     END IF
( 4523)     IF(nfstat /= nf90_noerr) THEN
( 4524)        status = NCDF_NFERR
( 4525)        RETURN
( 4526)     END IF
( 4527) 
( 4528)     data = databuf
( 4529) 
( 4530)     ! All done, close up the netCDF file
( 4531)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4532)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Closing file:', filename
( 4533)        CALL FLUSH
( 4534)     END IF
( 4535)     nfstat = nf90_close(ncid)
( 4536)     IF(nfstat /= nf90_noerr) THEN
( 4537)        status = NCDF_NFERR
( 4538)        RETURN
( 4539)     END IF
( 4540) 
( 4541)     DEALLOCATE(databuf)
( 4542) 
( 4543)     ! Sync up processors before returning from subroutine
( 4544)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4545)     IF(mpistat /= 0) THEN
( 4546)        status = NCDF_MPERR
( 4547)        RETURN
( 4548)     END IF
( 4549) 
( 4550)   END SUBROUTINE ncdf_read3d_global
( 4551) 
( 4552)   SUBROUTINE ncdf_read4d_global(filename, varname, data, dsz, time, status)
( 4553)     IMPLICIT NONE
( 4554)     ! Subroutine argument declarations
( 4555)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4556)     INTEGER,INTENT(IN) :: time, dsz
( 4557)     REAL,DIMENSION(:,:,:,:),INTENT(OUT) :: data
( 4558)     INTEGER,INTENT(OUT) :: status
( 4559)     
( 4560)     ! Local declarations
( 4561)     INTEGER :: ncid,    &  ! netCDF file ID
( 4562)                varid,   &  ! ID of netCDF variable to be written to
( 4563)                mpistat, &  ! MPI library call return status
( 4564)                nfstat,  &  ! netCDF library call return status
( 4565)                i,j,k,       &  ! Loop counter
( 4566)                ndims       ! Number of dimensions in this variable
( 4567) 
( 4568)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 4569)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 4570)     LOGICAL :: hastime   ! Whether this variable has a time axis






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 80

( 4571)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 4572)     INTEGER :: tindex
( 4573)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:,:) :: databuf
( 4574)     
( 4575)     ! Initializations
( 4576)     status = NCDF_NOERR
( 4577)     hastime = .FALSE.
( 4578)     j=0
( 4579) 
( 4580)     if(time > 0) then
( 4581)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 4582)     else
( 4583)        tindex = -time
( 4584)     endif
( 4585) !    tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 4586)     
( 4587)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4588)        WRITE(100,*) 'NCDF DEBUG: Reading 4D array from file:', filename
( 4589)        CALL FLUSH
( 4590)     END IF
( 4591)     
( 4592)     ! Open netCDF file and get info
( 4593)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4594)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Opening file:', filename
( 4595)        CALL FLUSH
( 4596)     END IF
( 4597) 
( 4598)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4599)     IF(nfstat /= nf90_noerr) THEN
( 4600)        status = NCDF_NFERR
( 4601)        RETURN
( 4602)     END IF
( 4603)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4604)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Getting info from file:', filename
( 4605)        CALL FLUSH
( 4606)     END IF
( 4607)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 4608)     IF(nfstat /= nf90_noerr) THEN
( 4609)        status = NCDF_NFERR
( 4610)        RETURN
( 4611)     END IF
( 4612)     
( 4613)     ! Determine if this variable contains a time axis
( 4614)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 4615)          ndims=ndims, dimids=var_dimids)
( 4616)     IF(nfstat /= nf90_noerr) THEN
( 4617)        status = NCDF_NFERR
( 4618)        RETURN
( 4619)     END IF
( 4620)     DO i = 1, ndims
( 4621)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 4622)             name=dname, len=k)
( 4623) !       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 4624) !            name=dname)
( 4625)        dname=TRIM(dname)
( 4626)        var_dimlens(i-j) = k
( 4627)        IF(nfstat /= nf90_noerr) THEN
( 4628)           status = NCDF_NFERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 81

( 4629)           RETURN
( 4630)        END IF
( 4631)        IF(dname=='time_counter') THEN
( 4632)           hastime = .TRUE.
( 4633)           j = 1
( 4634)           EXIT
( 4635)        END IF
( 4636)     END DO
( 4637)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3),1:var_dimlens(4)))
( 4638)     
( 4639)     ! Read data
( 4640)     IF(hastime .EQV. .TRUE.) THEN
( 4641)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4642)             (/ 1, 1, 1, 1, tindex/), &
( 4643)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), var_dimlens(4), 1 /))
( 4644)     ELSE
( 4645)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4646)             (/ 1, 1, 1, 1 /), &
( 4647)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), var_dimlens(4) /))
( 4648)     END IF
( 4649)     IF(nfstat /= nf90_noerr) THEN
( 4650)        status = NCDF_NFERR
( 4651)        RETURN
( 4652)     END IF
( 4653) 
( 4654)     data = databuf
( 4655)     DEALLOCATE(databuf)
( 4656)     
( 4657)     ! All done, close up the netCDF file
( 4658)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4659)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Closing file:', filename
( 4660)        CALL FLUSH
( 4661)     END IF
( 4662)     nfstat = nf90_close(ncid)
( 4663)     IF(nfstat /= nf90_noerr) THEN
( 4664)        status = NCDF_NFERR
( 4665)        RETURN
( 4666)     END IF
( 4667)     
( 4668)     ! Sync up processors before returning from subroutine
( 4669)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4670)     IF(mpistat /= 0) THEN
( 4671)        status = NCDF_MPERR
( 4672)        RETURN
( 4673)     END IF
( 4674) 
( 4675)   END SUBROUTINE ncdf_read4d_global
( 4676) 
( 4677) 
( 4678) 
( 4679) END MODULE lib_ncdf













PGF90 (Version     10.5)          09/23/2010  09:49:05      page 82

(    1) # 1 "lib_ncdf.F90"
(    1) !-------------------------------------------------------------------------------
(    2) !      OPA NETCDF I/O LIBRARY
(    3) !     ========================
(    4) !
(    5) ! This module is intended to provide a workaround for some of the difficulties,
(    6) ! eccentricities or general insanities of the IOIPSL library that is integrated
(    7) ! with OPA. All output routines in this module write to a single file, rather
(    8) ! than a separate file for each processor (thus avoiding the necessity of
(    9) ! merging files offline). There is no measurable time overhead incurred by
(   10) ! doing this, so we're not sacrificing model speed to avoid merging files.
(   11) ! Similarly, all read routines need only a single file. This library also allows
(   12) ! for easy creation and structuring of new output files. Unlike IOIPSL, it
(   13) ! performs no averaging on the data - these routines simply write or read
(   14) ! the requested data, and nothing else. You also don't need to worry about
(   15) ! opening and closing files, ID numbers, or anything else - that's all handled
(   16) ! for you. This library should work for any number of processors in any layout.
(   17) !
(   18) ! The 2D and 3D variable writing routines attempt to deal with time axes
(   19) ! somewhat 'intelligently'. When a 2D array is passed to ncdf_write, the
(   20) ! corresponding subroutine checks the specified variable for a time dimension.
(   21) ! If one exists, the variable in the netCDF file is, obviously, 3D, even though
(   22) ! only a 2D array was passed. The subroutine determines the correct index on
(   23) ! the time axis to write to, then puts the data there. If no time axis is found,
(   24) ! the array is simply written to the variable as you'd expect. The same sort
(   25) ! of thing happens when a 3D array is passed. I hope that made some sense.
(   26) !
(   27) ! NOTE - I've had issues with the FLUSH subroutine causing segmentation faults
(   28) ! on some systems. If this happens, just comment out any calls to FLUSH and
(   29) ! recompiler. The problem should go away.
(   30) !
(   31) ! Chris Nickerson
(   32) ! October & November 2007
(   33) ! nickersonc@mar.dfo-mpo.gc.ca
(   34) !-------------------------------------------------------------------------------
(   35) 
(   36) MODULE lib_ncdf
(   37)   USE netcdf
(   38)   USE par_oce
(   39)   USE dom_oce
(   40)   USE in_out_manager
(   41)   USE calendar
(   42)   USE daymod
(   43) 
(   44)   IMPLICIT NONE
(   45) 
(   46) !!DB: PRIVATE statement omitted as the include below becomes hidden to USErs
(   47) # 1 "/usr/include/mpich2-x86_64/mpif.h"
(    1) !      /* -*- Mode: Fortran; -*- */
(    2) !      
(    3) !      (C) 2001 by Argonne National Laboratory.
(    4) !      See COPYRIGHT in top-level directory.
(    5) !      
(    6) !      DO NOT EDIT
(    7) !      This file created by buildiface 
(    8) !      
(    9)        INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
(   10)        PARAMETER (MPI_SOURCE=3,MPI_TAG=4,MPI_ERROR=5)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 83

(   11)        INTEGER MPI_STATUS_SIZE
(   12)        PARAMETER (MPI_STATUS_SIZE=5)
(   13)        INTEGER MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
(   14)        INTEGER MPI_STATUSES_IGNORE(MPI_STATUS_SIZE,1)
(   15)        INTEGER MPI_ERRCODES_IGNORE(1)
(   16)        CHARACTER*1 MPI_ARGVS_NULL(1,1)
(   17)        CHARACTER*1 MPI_ARGV_NULL(1)
(   18)        INTEGER MPI_SUCCESS
(   19)        PARAMETER (MPI_SUCCESS=0)
(   20)        INTEGER MPI_ERR_OTHER
(   21)        PARAMETER (MPI_ERR_OTHER=15)
(   22)        INTEGER MPI_ERR_WIN
(   23)        PARAMETER (MPI_ERR_WIN=45)
(   24)        INTEGER MPI_ERR_FILE
(   25)        PARAMETER (MPI_ERR_FILE=27)
(   26)        INTEGER MPI_ERR_COUNT
(   27)        PARAMETER (MPI_ERR_COUNT=2)
(   28)        INTEGER MPI_ERR_SPAWN
(   29)        PARAMETER (MPI_ERR_SPAWN=42)
(   30)        INTEGER MPI_ERR_BASE
(   31)        PARAMETER (MPI_ERR_BASE=46)
(   32)        INTEGER MPI_ERR_RMA_CONFLICT
(   33)        PARAMETER (MPI_ERR_RMA_CONFLICT=49)
(   34)        INTEGER MPI_ERR_IN_STATUS
(   35)        PARAMETER (MPI_ERR_IN_STATUS=17)
(   36)        INTEGER MPI_ERR_INFO_KEY
(   37)        PARAMETER (MPI_ERR_INFO_KEY=29)
(   38)        INTEGER MPI_ERR_LOCKTYPE
(   39)        PARAMETER (MPI_ERR_LOCKTYPE=47)
(   40)        INTEGER MPI_ERR_OP
(   41)        PARAMETER (MPI_ERR_OP=9)
(   42)        INTEGER MPI_ERR_ARG
(   43)        PARAMETER (MPI_ERR_ARG=12)
(   44)        INTEGER MPI_ERR_READ_ONLY
(   45)        PARAMETER (MPI_ERR_READ_ONLY=40)
(   46)        INTEGER MPI_ERR_SIZE
(   47)        PARAMETER (MPI_ERR_SIZE=51)
(   48)        INTEGER MPI_ERR_BUFFER
(   49)        PARAMETER (MPI_ERR_BUFFER=1)
(   50)        INTEGER MPI_ERR_DUP_DATAREP
(   51)        PARAMETER (MPI_ERR_DUP_DATAREP=24)
(   52)        INTEGER MPI_ERR_UNSUPPORTED_DATAREP
(   53)        PARAMETER (MPI_ERR_UNSUPPORTED_DATAREP=43)
(   54)        INTEGER MPI_ERR_LASTCODE
(   55)        PARAMETER (MPI_ERR_LASTCODE=1073741823)
(   56)        INTEGER MPI_ERR_TRUNCATE
(   57)        PARAMETER (MPI_ERR_TRUNCATE=14)
(   58)        INTEGER MPI_ERR_DISP
(   59)        PARAMETER (MPI_ERR_DISP=52)
(   60)        INTEGER MPI_ERR_PORT
(   61)        PARAMETER (MPI_ERR_PORT=38)
(   62)        INTEGER MPI_ERR_INFO_NOKEY
(   63)        PARAMETER (MPI_ERR_INFO_NOKEY=31)
(   64)        INTEGER MPI_ERR_ASSERT
(   65)        PARAMETER (MPI_ERR_ASSERT=53)
(   66)        INTEGER MPI_ERR_FILE_EXISTS
(   67)        PARAMETER (MPI_ERR_FILE_EXISTS=25)
(   68)        INTEGER MPI_ERR_PENDING






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 84

(   69)        PARAMETER (MPI_ERR_PENDING=18)
(   70)        INTEGER MPI_ERR_COMM
(   71)        PARAMETER (MPI_ERR_COMM=5)
(   72)        INTEGER MPI_ERR_KEYVAL
(   73)        PARAMETER (MPI_ERR_KEYVAL=48)
(   74)        INTEGER MPI_ERR_NAME
(   75)        PARAMETER (MPI_ERR_NAME=33)
(   76)        INTEGER MPI_ERR_REQUEST
(   77)        PARAMETER (MPI_ERR_REQUEST=19)
(   78)        INTEGER MPI_ERR_GROUP
(   79)        PARAMETER (MPI_ERR_GROUP=8)
(   80)        INTEGER MPI_ERR_TOPOLOGY
(   81)        PARAMETER (MPI_ERR_TOPOLOGY=10)
(   82)        INTEGER MPI_ERR_TYPE
(   83)        PARAMETER (MPI_ERR_TYPE=3)
(   84)        INTEGER MPI_ERR_TAG
(   85)        PARAMETER (MPI_ERR_TAG=4)
(   86)        INTEGER MPI_ERR_INFO_VALUE
(   87)        PARAMETER (MPI_ERR_INFO_VALUE=30)
(   88)        INTEGER MPI_ERR_NOT_SAME
(   89)        PARAMETER (MPI_ERR_NOT_SAME=35)
(   90)        INTEGER MPI_ERR_RMA_SYNC
(   91)        PARAMETER (MPI_ERR_RMA_SYNC=50)
(   92)        INTEGER MPI_ERR_INFO
(   93)        PARAMETER (MPI_ERR_INFO=28)
(   94)        INTEGER MPI_ERR_NO_MEM
(   95)        PARAMETER (MPI_ERR_NO_MEM=34)
(   96)        INTEGER MPI_ERR_BAD_FILE
(   97)        PARAMETER (MPI_ERR_BAD_FILE=22)
(   98)        INTEGER MPI_ERR_FILE_IN_USE
(   99)        PARAMETER (MPI_ERR_FILE_IN_USE=26)
(  100)        INTEGER MPI_ERR_UNKNOWN
(  101)        PARAMETER (MPI_ERR_UNKNOWN=13)
(  102)        INTEGER MPI_ERR_UNSUPPORTED_OPERATION
(  103)        PARAMETER (MPI_ERR_UNSUPPORTED_OPERATION=44)
(  104)        INTEGER MPI_ERR_QUOTA
(  105)        PARAMETER (MPI_ERR_QUOTA=39)
(  106)        INTEGER MPI_ERR_AMODE
(  107)        PARAMETER (MPI_ERR_AMODE=21)
(  108)        INTEGER MPI_ERR_ROOT
(  109)        PARAMETER (MPI_ERR_ROOT=7)
(  110)        INTEGER MPI_ERR_RANK
(  111)        PARAMETER (MPI_ERR_RANK=6)
(  112)        INTEGER MPI_ERR_DIMS
(  113)        PARAMETER (MPI_ERR_DIMS=11)
(  114)        INTEGER MPI_ERR_NO_SUCH_FILE
(  115)        PARAMETER (MPI_ERR_NO_SUCH_FILE=37)
(  116)        INTEGER MPI_ERR_SERVICE
(  117)        PARAMETER (MPI_ERR_SERVICE=41)
(  118)        INTEGER MPI_ERR_INTERN
(  119)        PARAMETER (MPI_ERR_INTERN=16)
(  120)        INTEGER MPI_ERR_IO
(  121)        PARAMETER (MPI_ERR_IO=32)
(  122)        INTEGER MPI_ERR_ACCESS
(  123)        PARAMETER (MPI_ERR_ACCESS=20)
(  124)        INTEGER MPI_ERR_NO_SPACE
(  125)        PARAMETER (MPI_ERR_NO_SPACE=36)
(  126)        INTEGER MPI_ERR_CONVERSION






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 85

(  127)        PARAMETER (MPI_ERR_CONVERSION=23)
(  128)        INTEGER MPI_ERRORS_ARE_FATAL
(  129)        PARAMETER (MPI_ERRORS_ARE_FATAL=1409286144)
(  130)        INTEGER MPI_ERRORS_RETURN
(  131)        PARAMETER (MPI_ERRORS_RETURN=1409286145)
(  132)        INTEGER MPI_IDENT
(  133)        PARAMETER (MPI_IDENT=0)
(  134)        INTEGER MPI_CONGRUENT
(  135)        PARAMETER (MPI_CONGRUENT=1)
(  136)        INTEGER MPI_SIMILAR
(  137)        PARAMETER (MPI_SIMILAR=2)
(  138)        INTEGER MPI_UNEQUAL
(  139)        PARAMETER (MPI_UNEQUAL=3)
(  140)        INTEGER MPI_MAX
(  141)        PARAMETER (MPI_MAX=1476395009)
(  142)        INTEGER MPI_MIN
(  143)        PARAMETER (MPI_MIN=1476395010)
(  144)        INTEGER MPI_SUM
(  145)        PARAMETER (MPI_SUM=1476395011)
(  146)        INTEGER MPI_PROD
(  147)        PARAMETER (MPI_PROD=1476395012)
(  148)        INTEGER MPI_LAND
(  149)        PARAMETER (MPI_LAND=1476395013)
(  150)        INTEGER MPI_BAND
(  151)        PARAMETER (MPI_BAND=1476395014)
(  152)        INTEGER MPI_LOR
(  153)        PARAMETER (MPI_LOR=1476395015)
(  154)        INTEGER MPI_BOR
(  155)        PARAMETER (MPI_BOR=1476395016)
(  156)        INTEGER MPI_LXOR
(  157)        PARAMETER (MPI_LXOR=1476395017)
(  158)        INTEGER MPI_BXOR
(  159)        PARAMETER (MPI_BXOR=1476395018)
(  160)        INTEGER MPI_MINLOC
(  161)        PARAMETER (MPI_MINLOC=1476395019)
(  162)        INTEGER MPI_MAXLOC
(  163)        PARAMETER (MPI_MAXLOC=1476395020)
(  164)        INTEGER MPI_REPLACE
(  165)        PARAMETER (MPI_REPLACE=1476395021)
(  166)        INTEGER MPI_COMM_WORLD
(  167)        PARAMETER (MPI_COMM_WORLD=1140850688)
(  168)        INTEGER MPI_COMM_SELF
(  169)        PARAMETER (MPI_COMM_SELF=1140850689)
(  170)        INTEGER MPI_GROUP_EMPTY
(  171)        PARAMETER (MPI_GROUP_EMPTY=1207959552)
(  172)        INTEGER MPI_COMM_NULL
(  173)        PARAMETER (MPI_COMM_NULL=67108864)
(  174)        INTEGER MPI_WIN_NULL
(  175)        PARAMETER (MPI_WIN_NULL=536870912)
(  176)        INTEGER MPI_FILE_NULL
(  177)        PARAMETER (MPI_FILE_NULL=0)
(  178)        INTEGER MPI_GROUP_NULL
(  179)        PARAMETER (MPI_GROUP_NULL=134217728)
(  180)        INTEGER MPI_OP_NULL
(  181)        PARAMETER (MPI_OP_NULL=402653184)
(  182)        INTEGER MPI_DATATYPE_NULL
(  183)        PARAMETER (MPI_DATATYPE_NULL=201326592)
(  184)        INTEGER MPI_REQUEST_NULL






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 86

(  185)        PARAMETER (MPI_REQUEST_NULL=738197504)
(  186)        INTEGER MPI_ERRHANDLER_NULL
(  187)        PARAMETER (MPI_ERRHANDLER_NULL=335544320)
(  188)        INTEGER MPI_INFO_NULL
(  189)        PARAMETER (MPI_INFO_NULL=469762048)
(  190)        INTEGER MPI_TAG_UB
(  191)        PARAMETER (MPI_TAG_UB=1681915906)
(  192)        INTEGER MPI_HOST
(  193)        PARAMETER (MPI_HOST=1681915908)
(  194)        INTEGER MPI_IO
(  195)        PARAMETER (MPI_IO=1681915910)
(  196)        INTEGER MPI_WTIME_IS_GLOBAL
(  197)        PARAMETER (MPI_WTIME_IS_GLOBAL=1681915912)
(  198)        INTEGER MPI_UNIVERSE_SIZE
(  199)        PARAMETER (MPI_UNIVERSE_SIZE=1681915914)
(  200)        INTEGER MPI_LASTUSEDCODE
(  201)        PARAMETER (MPI_LASTUSEDCODE=1681915916)
(  202)        INTEGER MPI_APPNUM
(  203)        PARAMETER (MPI_APPNUM=1681915918)
(  204)        INTEGER MPI_WIN_BASE
(  205)        PARAMETER (MPI_WIN_BASE=1711276034)
(  206)        INTEGER MPI_WIN_SIZE
(  207)        PARAMETER (MPI_WIN_SIZE=1711276036)
(  208)        INTEGER MPI_WIN_DISP_UNIT
(  209)        PARAMETER (MPI_WIN_DISP_UNIT=1711276038)
(  210)        INTEGER MPI_MAX_ERROR_STRING
(  211)        PARAMETER (MPI_MAX_ERROR_STRING=1023)
(  212)        INTEGER MPI_MAX_PORT_NAME
(  213)        PARAMETER (MPI_MAX_PORT_NAME=255)
(  214)        INTEGER MPI_MAX_OBJECT_NAME
(  215)        PARAMETER (MPI_MAX_OBJECT_NAME=127)
(  216)        INTEGER MPI_MAX_INFO_KEY
(  217)        PARAMETER (MPI_MAX_INFO_KEY=254)
(  218)        INTEGER MPI_MAX_INFO_VAL
(  219)        PARAMETER (MPI_MAX_INFO_VAL=1023)
(  220)        INTEGER MPI_MAX_PROCESSOR_NAME
(  221)        PARAMETER (MPI_MAX_PROCESSOR_NAME=128-1)
(  222)        INTEGER MPI_MAX_DATAREP_STRING
(  223)        PARAMETER (MPI_MAX_DATAREP_STRING=127)
(  224)        INTEGER MPI_UNDEFINED
(  225)        PARAMETER (MPI_UNDEFINED=(-32766))
(  226)        INTEGER MPI_KEYVAL_INVALID
(  227)        PARAMETER (MPI_KEYVAL_INVALID=603979776)
(  228)        INTEGER MPI_BSEND_OVERHEAD
(  229)        PARAMETER (MPI_BSEND_OVERHEAD=88)
(  230)        INTEGER MPI_PROC_NULL
(  231)        PARAMETER (MPI_PROC_NULL=-1)
(  232)        INTEGER MPI_ANY_SOURCE
(  233)        PARAMETER (MPI_ANY_SOURCE=-2)
(  234)        INTEGER MPI_ANY_TAG
(  235)        PARAMETER (MPI_ANY_TAG=-1)
(  236)        INTEGER MPI_ROOT
(  237)        PARAMETER (MPI_ROOT=-3)
(  238)        INTEGER MPI_GRAPH
(  239)        PARAMETER (MPI_GRAPH=1)
(  240)        INTEGER MPI_CART
(  241)        PARAMETER (MPI_CART=2)
(  242)        INTEGER MPI_VERSION






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 87

(  243)        PARAMETER (MPI_VERSION=2)
(  244)        INTEGER MPI_SUBVERSION
(  245)        PARAMETER (MPI_SUBVERSION=2)
(  246)        INTEGER MPI_LOCK_EXCLUSIVE
(  247)        PARAMETER (MPI_LOCK_EXCLUSIVE=234)
(  248)        INTEGER MPI_LOCK_SHARED
(  249)        PARAMETER (MPI_LOCK_SHARED=235)
(  250)        INTEGER MPI_COMPLEX
(  251)        PARAMETER (MPI_COMPLEX=1275070494)
(  252)        INTEGER MPI_DOUBLE_COMPLEX
(  253)        PARAMETER (MPI_DOUBLE_COMPLEX=1275072546)
(  254)        INTEGER MPI_LOGICAL
(  255)        PARAMETER (MPI_LOGICAL=1275069469)
(  256)        INTEGER MPI_REAL
(  257)        PARAMETER (MPI_REAL=1275069468)
(  258)        INTEGER MPI_DOUBLE_PRECISION
(  259)        PARAMETER (MPI_DOUBLE_PRECISION=1275070495)
(  260)        INTEGER MPI_INTEGER
(  261)        PARAMETER (MPI_INTEGER=1275069467)
(  262)        INTEGER MPI_2INTEGER
(  263)        PARAMETER (MPI_2INTEGER=1275070496)
(  264)        INTEGER MPI_2COMPLEX
(  265)        PARAMETER (MPI_2COMPLEX=1275072548)
(  266)        INTEGER MPI_2DOUBLE_PRECISION
(  267)        PARAMETER (MPI_2DOUBLE_PRECISION=1275072547)
(  268)        INTEGER MPI_2REAL
(  269)        PARAMETER (MPI_2REAL=1275070497)
(  270)        INTEGER MPI_2DOUBLE_COMPLEX
(  271)        PARAMETER (MPI_2DOUBLE_COMPLEX=1275076645)
(  272)        INTEGER MPI_CHARACTER
(  273)        PARAMETER (MPI_CHARACTER=1275068698)
(  274)        INTEGER MPI_BYTE
(  275)        PARAMETER (MPI_BYTE=1275068685)
(  276)        INTEGER MPI_UB
(  277)        PARAMETER (MPI_UB=1275068433)
(  278)        INTEGER MPI_LB
(  279)        PARAMETER (MPI_LB=1275068432)
(  280)        INTEGER MPI_PACKED
(  281)        PARAMETER (MPI_PACKED=1275068687)
(  282)        INTEGER MPI_INTEGER1
(  283)        PARAMETER (MPI_INTEGER1=1275068717)
(  284)        INTEGER MPI_INTEGER2
(  285)        PARAMETER (MPI_INTEGER2=1275068975)
(  286)        INTEGER MPI_INTEGER4
(  287)        PARAMETER (MPI_INTEGER4=1275069488)
(  288)        INTEGER MPI_INTEGER8
(  289)        PARAMETER (MPI_INTEGER8=1275070513)
(  290)        INTEGER MPI_INTEGER16
(  291)        PARAMETER (MPI_INTEGER16=MPI_DATATYPE_NULL)
(  292)        INTEGER MPI_REAL4
(  293)        PARAMETER (MPI_REAL4=1275069479)
(  294)        INTEGER MPI_REAL8
(  295)        PARAMETER (MPI_REAL8=1275070505)
(  296)        INTEGER MPI_REAL16
(  297)        PARAMETER (MPI_REAL16=1275072555)
(  298)        INTEGER MPI_COMPLEX8
(  299)        PARAMETER (MPI_COMPLEX8=1275070504)
(  300)        INTEGER MPI_COMPLEX16






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 88

(  301)        PARAMETER (MPI_COMPLEX16=1275072554)
(  302)        INTEGER MPI_COMPLEX32
(  303)        PARAMETER (MPI_COMPLEX32=1275076652)
(  304)        INTEGER MPI_ADDRESS_KIND, MPI_OFFSET_KIND
(  305)        PARAMETER (MPI_ADDRESS_KIND=8)
(  306)        PARAMETER (MPI_OFFSET_KIND=8)
(  307)        INTEGER MPI_CHAR
(  308)        PARAMETER (MPI_CHAR=1275068673)
(  309)        INTEGER MPI_SIGNED_CHAR
(  310)        PARAMETER (MPI_SIGNED_CHAR=1275068696)
(  311)        INTEGER MPI_UNSIGNED_CHAR
(  312)        PARAMETER (MPI_UNSIGNED_CHAR=1275068674)
(  313)        INTEGER MPI_WCHAR
(  314)        PARAMETER (MPI_WCHAR=1275069454)
(  315)        INTEGER MPI_SHORT
(  316)        PARAMETER (MPI_SHORT=1275068931)
(  317)        INTEGER MPI_UNSIGNED_SHORT
(  318)        PARAMETER (MPI_UNSIGNED_SHORT=1275068932)
(  319)        INTEGER MPI_INT
(  320)        PARAMETER (MPI_INT=1275069445)
(  321)        INTEGER MPI_UNSIGNED
(  322)        PARAMETER (MPI_UNSIGNED=1275069446)
(  323)        INTEGER MPI_LONG
(  324)        PARAMETER (MPI_LONG=1275070471)
(  325)        INTEGER MPI_UNSIGNED_LONG
(  326)        PARAMETER (MPI_UNSIGNED_LONG=1275070472)
(  327)        INTEGER MPI_FLOAT
(  328)        PARAMETER (MPI_FLOAT=1275069450)
(  329)        INTEGER MPI_DOUBLE
(  330)        PARAMETER (MPI_DOUBLE=1275070475)
(  331)        INTEGER MPI_LONG_DOUBLE
(  332)        PARAMETER (MPI_LONG_DOUBLE=1275072524)
(  333)        INTEGER MPI_LONG_LONG_INT
(  334)        PARAMETER (MPI_LONG_LONG_INT=1275070473)
(  335)        INTEGER MPI_UNSIGNED_LONG_LONG
(  336)        PARAMETER (MPI_UNSIGNED_LONG_LONG=1275070489)
(  337)        INTEGER MPI_LONG_LONG
(  338)        PARAMETER (MPI_LONG_LONG=1275070473)
(  339)        INTEGER MPI_FLOAT_INT
(  340)        PARAMETER (MPI_FLOAT_INT=-1946157056)
(  341)        INTEGER MPI_DOUBLE_INT
(  342)        PARAMETER (MPI_DOUBLE_INT=-1946157055)
(  343)        INTEGER MPI_LONG_INT
(  344)        PARAMETER (MPI_LONG_INT=-1946157054)
(  345)        INTEGER MPI_SHORT_INT
(  346)        PARAMETER (MPI_SHORT_INT=-1946157053)
(  347)        INTEGER MPI_2INT
(  348)        PARAMETER (MPI_2INT=1275070486)
(  349)        INTEGER MPI_LONG_DOUBLE_INT
(  350)        PARAMETER (MPI_LONG_DOUBLE_INT=-1946157052)
(  351)        INTEGER MPI_INT8_T
(  352)        PARAMETER (MPI_INT8_T=1275068727)
(  353)        INTEGER MPI_INT16_T
(  354)        PARAMETER (MPI_INT16_T=1275068984)
(  355)        INTEGER MPI_INT32_T
(  356)        PARAMETER (MPI_INT32_T=1275069497)
(  357)        INTEGER MPI_INT64_T
(  358)        PARAMETER (MPI_INT64_T=1275070522)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 89

(  359)        INTEGER MPI_UINT8_T
(  360)        PARAMETER (MPI_UINT8_T=1275068731)
(  361)        INTEGER MPI_UINT16_T
(  362)        PARAMETER (MPI_UINT16_T=1275068988)
(  363)        INTEGER MPI_UINT32_T
(  364)        PARAMETER (MPI_UINT32_T=1275069501)
(  365)        INTEGER MPI_UINT64_T
(  366)        PARAMETER (MPI_UINT64_T=1275070526)
(  367)        INTEGER MPI_C_BOOL
(  368)        PARAMETER (MPI_C_BOOL=1275068735)
(  369)        INTEGER MPI_C_FLOAT_COMPLEX
(  370)        PARAMETER (MPI_C_FLOAT_COMPLEX=1275070528)
(  371)        INTEGER MPI_C_COMPLEX
(  372)        PARAMETER (MPI_C_COMPLEX=1275070528)
(  373)        INTEGER MPI_C_DOUBLE_COMPLEX
(  374)        PARAMETER (MPI_C_DOUBLE_COMPLEX=1275072577)
(  375)        INTEGER MPI_C_LONG_DOUBLE_COMPLEX
(  376)        PARAMETER (MPI_C_LONG_DOUBLE_COMPLEX=1275076674)
(  377)        INTEGER MPI_AINT
(  378)        PARAMETER (MPI_AINT=1275070531)
(  379)        INTEGER MPI_OFFSET
(  380)        PARAMETER (MPI_OFFSET=1275070532)
(  381)        INTEGER MPI_COMBINER_NAMED
(  382)        PARAMETER (MPI_COMBINER_NAMED=1)
(  383)        INTEGER MPI_COMBINER_DUP
(  384)        PARAMETER (MPI_COMBINER_DUP=2)
(  385)        INTEGER MPI_COMBINER_CONTIGUOUS
(  386)        PARAMETER (MPI_COMBINER_CONTIGUOUS=3)
(  387)        INTEGER MPI_COMBINER_VECTOR
(  388)        PARAMETER (MPI_COMBINER_VECTOR=4)
(  389)        INTEGER MPI_COMBINER_HVECTOR_INTEGER
(  390)        PARAMETER (MPI_COMBINER_HVECTOR_INTEGER=5)
(  391)        INTEGER MPI_COMBINER_HVECTOR
(  392)        PARAMETER (MPI_COMBINER_HVECTOR=6)
(  393)        INTEGER MPI_COMBINER_INDEXED
(  394)        PARAMETER (MPI_COMBINER_INDEXED=7)
(  395)        INTEGER MPI_COMBINER_HINDEXED_INTEGER
(  396)        PARAMETER (MPI_COMBINER_HINDEXED_INTEGER=8)
(  397)        INTEGER MPI_COMBINER_HINDEXED
(  398)        PARAMETER (MPI_COMBINER_HINDEXED=9)
(  399)        INTEGER MPI_COMBINER_INDEXED_BLOCK
(  400)        PARAMETER (MPI_COMBINER_INDEXED_BLOCK=10)
(  401)        INTEGER MPI_COMBINER_STRUCT_INTEGER
(  402)        PARAMETER (MPI_COMBINER_STRUCT_INTEGER=11)
(  403)        INTEGER MPI_COMBINER_STRUCT
(  404)        PARAMETER (MPI_COMBINER_STRUCT=12)
(  405)        INTEGER MPI_COMBINER_SUBARRAY
(  406)        PARAMETER (MPI_COMBINER_SUBARRAY=13)
(  407)        INTEGER MPI_COMBINER_DARRAY
(  408)        PARAMETER (MPI_COMBINER_DARRAY=14)
(  409)        INTEGER MPI_COMBINER_F90_REAL
(  410)        PARAMETER (MPI_COMBINER_F90_REAL=15)
(  411)        INTEGER MPI_COMBINER_F90_COMPLEX
(  412)        PARAMETER (MPI_COMBINER_F90_COMPLEX=16)
(  413)        INTEGER MPI_COMBINER_F90_INTEGER
(  414)        PARAMETER (MPI_COMBINER_F90_INTEGER=17)
(  415)        INTEGER MPI_COMBINER_RESIZED
(  416)        PARAMETER (MPI_COMBINER_RESIZED=18)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 90

(  417)        INTEGER MPI_TYPECLASS_REAL
(  418)        PARAMETER (MPI_TYPECLASS_REAL=1)
(  419)        INTEGER MPI_TYPECLASS_INTEGER
(  420)        PARAMETER (MPI_TYPECLASS_INTEGER=2)
(  421)        INTEGER MPI_TYPECLASS_COMPLEX
(  422)        PARAMETER (MPI_TYPECLASS_COMPLEX=3)
(  423)        INTEGER MPI_MODE_NOCHECK
(  424)        PARAMETER (MPI_MODE_NOCHECK=1024)
(  425)        INTEGER MPI_MODE_NOSTORE
(  426)        PARAMETER (MPI_MODE_NOSTORE=2048)
(  427)        INTEGER MPI_MODE_NOPUT
(  428)        PARAMETER (MPI_MODE_NOPUT=4096)
(  429)        INTEGER MPI_MODE_NOPRECEDE
(  430)        PARAMETER (MPI_MODE_NOPRECEDE=8192)
(  431)        INTEGER MPI_MODE_NOSUCCEED
(  432)        PARAMETER (MPI_MODE_NOSUCCEED=16384)
(  433)        INTEGER MPI_THREAD_SINGLE
(  434)        PARAMETER (MPI_THREAD_SINGLE=0)
(  435)        INTEGER MPI_THREAD_FUNNELED
(  436)        PARAMETER (MPI_THREAD_FUNNELED=1)
(  437)        INTEGER MPI_THREAD_SERIALIZED
(  438)        PARAMETER (MPI_THREAD_SERIALIZED=2)
(  439)        INTEGER MPI_THREAD_MULTIPLE
(  440)        PARAMETER (MPI_THREAD_MULTIPLE=3)
(  441)        INTEGER MPI_MODE_RDONLY
(  442)        PARAMETER (MPI_MODE_RDONLY=2)
(  443)        INTEGER MPI_MODE_RDWR
(  444)        PARAMETER (MPI_MODE_RDWR=8)
(  445)        INTEGER MPI_MODE_WRONLY
(  446)        PARAMETER (MPI_MODE_WRONLY=4)
(  447)        INTEGER MPI_MODE_DELETE_ON_CLOSE
(  448)        PARAMETER (MPI_MODE_DELETE_ON_CLOSE=16)
(  449)        INTEGER MPI_MODE_UNIQUE_OPEN
(  450)        PARAMETER (MPI_MODE_UNIQUE_OPEN=32)
(  451)        INTEGER MPI_MODE_CREATE
(  452)        PARAMETER (MPI_MODE_CREATE=1)
(  453)        INTEGER MPI_MODE_EXCL
(  454)        PARAMETER (MPI_MODE_EXCL=64)
(  455)        INTEGER MPI_MODE_APPEND
(  456)        PARAMETER (MPI_MODE_APPEND=128)
(  457)        INTEGER MPI_MODE_SEQUENTIAL
(  458)        PARAMETER (MPI_MODE_SEQUENTIAL=256)
(  459)        INTEGER MPI_SEEK_SET
(  460)        PARAMETER (MPI_SEEK_SET=600)
(  461)        INTEGER MPI_SEEK_CUR
(  462)        PARAMETER (MPI_SEEK_CUR=602)
(  463)        INTEGER MPI_SEEK_END
(  464)        PARAMETER (MPI_SEEK_END=604)
(  465)        INTEGER MPI_ORDER_C
(  466)        PARAMETER (MPI_ORDER_C=56)
(  467)        INTEGER MPI_ORDER_FORTRAN
(  468)        PARAMETER (MPI_ORDER_FORTRAN=57)
(  469)        INTEGER MPI_DISTRIBUTE_BLOCK
(  470)        PARAMETER (MPI_DISTRIBUTE_BLOCK=121)
(  471)        INTEGER MPI_DISTRIBUTE_CYCLIC
(  472)        PARAMETER (MPI_DISTRIBUTE_CYCLIC=122)
(  473)        INTEGER MPI_DISTRIBUTE_NONE
(  474)        PARAMETER (MPI_DISTRIBUTE_NONE=123)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 91

(  475)        INTEGER MPI_DISTRIBUTE_DFLT_DARG
(  476)        PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
(  477)        integer*8 MPI_DISPLACEMENT_CURRENT
(  478)        PARAMETER (MPI_DISPLACEMENT_CURRENT=-54278278)
(  479)        INTEGER MPI_BOTTOM, MPI_IN_PLACE
(  480)        EXTERNAL MPI_DUP_FN, MPI_NULL_DELETE_FN, MPI_NULL_COPY_FN
(  481)        EXTERNAL MPI_WTIME, MPI_WTICK
(  482)        EXTERNAL PMPI_WTIME, PMPI_WTICK
(  483)        EXTERNAL MPI_COMM_DUP_FN, MPI_COMM_NULL_DELETE_FN
(  484)        EXTERNAL MPI_COMM_NULL_COPY_FN
(  485)        EXTERNAL MPI_WIN_DUP_FN, MPI_WIN_NULL_DELETE_FN
(  486)        EXTERNAL MPI_WIN_NULL_COPY_FN
(  487)        EXTERNAL MPI_TYPE_DUP_FN, MPI_TYPE_NULL_DELETE_FN
(  488)        EXTERNAL MPI_TYPE_NULL_COPY_FN
(  489)        EXTERNAL MPI_CONVERSION_FN_NULL
(  490)        REAL*8 MPI_WTIME, MPI_WTICK
(  491)        REAL*8 PMPI_WTIME, PMPI_WTICK
(  492) 
(  493) 
(  494)        COMMON /MPIPRIV1/ MPI_BOTTOM, MPI_IN_PLACE, MPI_STATUS_IGNORE
(  495) 
(  496)        COMMON /MPIPRIV2/ MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE
(  497)        SAVE /MPIPRIV1/,/MPIPRIV2/
(  498) 
(  499)        COMMON /MPIPRIVC/ MPI_ARGVS_NULL, MPI_ARGV_NULL
(  500)        SAVE   /MPIPRIVC/
(  501) # 48 "lib_ncdf.F90"
(   48) 
(   49) 
(   50) # 50
(   50)   ! Error constant definitions
(   51)   INTEGER,PARAMETER :: NCDF_NOERR = 0, &  ! No error, normal return status
(   52)                        NCDF_NFERR = 1, &  ! netCDF-related error occurred
(   53)                        NCDF_MPERR = 2, &  ! MPI-related error occurred
(   54)                        NCDF_ARERR = 3, &  ! Invalid arguments were given
(   55)                        NCDF_OTHER = 4     ! Some other error happened
(   56)   
(   57)   ! Datatype definitions used for creating variables
(   58)   INTEGER,PARAMETER :: NCDF_FLOAT = nf90_float, &
(   59)                        NCDF_DOUBLE = nf90_double
(   60) 
(   61)   ! Turn debugging output on/off (goes to fort.100 - it's VERY verbose)
(   62)   LOGICAL,PARAMETER :: DEBUG_OUT = .TRUE. ! True if you want debugging output,
(   63)                                            ! false otherwise. Debug output is
(   64)                                            ! written to OPA.out in the run dir
(   65)                                            ! Note that debug output is extremely
(   66)                                            ! verbose!
(   67)   
(   68) !!DB 2008.05.22 -- added ice restart
(   69) !!DB 2008.06.26 -- added get_dim_size
(   70)   ! Interface definitions
(   71)   PUBLIC :: ncdf_create_file_u, ncdf_create_file_v, ncdf_create_file_t,  ncdf_write, &
(   72)        ncdf_errstr, ncdf_create_file, ncdf_create_dim, ncdf_create_var, ncdf_put_att, &
(   73)        ncdf_create_restart,  ncdf_create_ice_restart, ncdf_get_dim_size, ncdf_readdate, &
(   74)        ncdf_read_global, ncdf_create_file_aveTSUV, ncdf_create_file_ave
(   75) 
(   76) !  PUBLIC :: output_special             ! special routines called by step.F90 if ave flag is on
(   77) !  PUBLIC :: output_aveTSUV             ! special routines called by step.F90 if ave flag is on






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 92

(   78) 
(   79)   
(   80)   INTERFACE ncdf_write
(   81)      MODULE PROCEDURE ncdf_writesv, ncdf_write1d, ncdf_write2d, ncdf_write3d, ncdf_write4d
(   82)   END INTERFACE
(   83) 
(   84)   INTERFACE ncdf_read
(   85)      MODULE PROCEDURE ncdf_readsv, ncdf_read1d, ncdf_read2d, ncdf_read3d, ncdf_read4d
(   86)   END INTERFACE
(   87) 
(   88)   INTERFACE ncdf_read_global
(   89)      MODULE PROCEDURE ncdf_read2d_global, ncdf_read3d_global, ncdf_read4d_global
(   90)   END INTERFACE
(   91) 
(   92)   INTERFACE ncdf_put_att
(   93)      MODULE PROCEDURE ncdf_put_att_int, ncdf_put_att_real, ncdf_put_att_char
(   94)   END INTERFACE
(   95) 
(   96) CONTAINS
(   97) 
(   98)   ! ncdf_create_file_u builds a standard u-grid OPA output file with all the default
(   99)   ! dimensions, variables and attributes
(  100)   SUBROUTINE ncdf_create_file_u(filename, op_type, status)
(  101)     IMPLICIT NONE
(  102)     ! Subroutine argument declarations
(  103)     CHARACTER(LEN=*),INTENT(IN) :: filename
(  104)     CHARACTER(LEN=*),INTENT(IN) :: op_type
(  105)     INTEGER,INTENT(OUT) :: status
(  106) 
(  107)     ! Local declarations
(  108)     INTEGER :: ncid,    &  ! netCDF file ID
(  109)                varid,   &  ! ID of netCDF variable to be written to
(  110)                nfstat,  &  ! netCDF library call return status
(  111)                mpistat     ! MPI library call return status
(  112)     INTEGER,DIMENSION(1:4) :: dimids
(  113)     INTEGER,DIMENSION(1:6) :: varids
(  114)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  115)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  116)     CHARACTER(LEN=100) :: sec_since    
(  117)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  118)     INTEGER :: int_opp, &              ! Operation interval
(  119)                int_wri                 ! Write interval
(  120)     CHARACTER(LEN=3),PARAMETER :: &
(  121)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  122)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  123)     
(  124)     ! Initializations
(  125)     status = NCDF_NOERR
(  126)     CALL ioget_calendar(cal_type)
(  127)     CALL ioget_timestamp(timestamp)
(  128)     WRITE (UNIT=sec_since, &
(  129)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  130)          &  nyear,nmonth,nday,0, 0, 0
(  131)     WRITE(t_origin, &
(  132)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  133)          &   nyear,months(nmonth),nday,0,0,0
(  134)     int_opp = nwrite * rdt
(  135)     int_wri = nwrite * rdt






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 93

(  136) 
(  137)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  138)        WRITE(100,*) 'NCDF DEBUG: Creating default U output file:', filename
(  139)        CALL FLUSH
(  140)     END IF
(  141) 
(  142)     ! Only processor 0 does anything
(  143)     IF(nproc == 0) THEN
(  144)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  145)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Creating file:', filename
(  146)           CALL FLUSH
(  147)        END IF
(  148)        ! Create the file
(  149)        nfstat = nf90_create(filename, nf90_clobber, ncid)
(  150)        IF(nfstat /= nf90_noerr) THEN
(  151)           status = NCDF_NFERR
(  152)           RETURN
(  153)        END IF
(  154)        
(  155)        ! Define dimensions
(  156)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  157)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining dimensions in file:', filename
(  158)           CALL FLUSH
(  159)        END IF
(  160)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  161)        nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
(  162)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
(  163)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
(  164)        
(  165)        ! Define variables
(  166)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  167)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining variables in file:', filename
(  168)           CALL FLUSH
(  169)        END IF
(  170)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  171)             (/ dimids(4), dimids(3) /), &
(  172)             varids(1))
(  173)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  174)             (/ dimids(4), dimids(3) /), &
(  175)             varids(2))
(  176)        nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
(  177)             (/ dimids(2) /), &
(  178)             varids(3))
(  179)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  180)             (/ dimids(1) /), &
(  181)             varids(4))
(  182)        nfstat = nf90_def_var(ncid, 'vozocrtx', nf90_float, &
(  183)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  184)             varids(5))
(  185)        nfstat = nf90_def_var(ncid, 'sozotaux', nf90_float, &
(  186)             (/ dimids(4), dimids(3), dimids(1) /), &
(  187)             varids(6))
(  188)        
(  189)        ! Add attributes
(  190)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  191)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Writing attributes in file:', filename
(  192)           CALL FLUSH
(  193)        END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 94

(  194)        ! nav_lon
(  195)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  196)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  197)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  198)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  199)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  200) 
(  201)        ! nav_lat
(  202)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  203)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  204)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  205)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  206)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  207) 
(  208)        ! depthu
(  209)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  210)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  211)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
(  212)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
(  213)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
(  214)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')
(  215) 
(  216)        ! time_counter
(  217)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  218)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  219)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  220)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
(  221)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
(  222) 
(  223)        ! vozocrtx
(  224)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
(  225)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
(  226)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
(  227)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
(  228)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Zonal Current')
(  229)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vozocrtx')
(  230)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
(  231)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
(  232)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
(  233)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
(  234)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')
(  235) 
(  236)        ! sozotaux
(  237)        nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
(  238)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
(  239)        nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
(  240)        nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
(  241)        nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along i-axis')
(  242)        nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sozotaux')
(  243)        nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
(  244)        nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')
(  245)        nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
(  246)        nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
(  247)        nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')
(  248) 
(  249)        ! global
(  250)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  251)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 95

(  252)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
(  253)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  254)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  255)       
(  256)        
(  257)        ! Close file
(  258)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  259)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Closing file:', filename
(  260)           CALL FLUSH
(  261)        END IF
(  262)        nfstat = nf90_close(ncid)
(  263)        IF(nfstat /= nf90_noerr) THEN
(  264)           status = NCDF_NFERR
(  265)           RETURN
(  266)        END IF
(  267)     END IF
(  268) 
(  269)     ! Sync up processors before returning from subroutine
(  270)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
(  271)     IF(mpistat /= 0) THEN
(  272)        status = NCDF_MPERR
(  273)        RETURN
(  274)     END IF
(  275)     
(  276)   END SUBROUTINE ncdf_create_file_u
(  277) 
(  278)   ! ncdf_create_file_v builds a standard v-grid OPA output file with all the default
(  279)   ! dimensions, variables and attributes
(  280)   SUBROUTINE ncdf_create_file_v(filename, op_type, status)
(  281)     IMPLICIT NONE
(  282)     ! Subroutine argument declarations
(  283)     CHARACTER(LEN=*),INTENT(IN) :: filename
(  284)     CHARACTER(LEN=*),INTENT(IN) :: op_type
(  285)     INTEGER,INTENT(OUT) :: status
(  286) 
(  287)     ! Local declarations
(  288)     INTEGER :: ncid,    &  ! netCDF file ID
(  289)                varid,   &  ! ID of netCDF variable to be written to
(  290)                nfstat,  &  ! netCDF library call return status
(  291)                mpistat     ! MPI library call return status
(  292)     INTEGER,DIMENSION(1:4) :: dimids
(  293)     INTEGER,DIMENSION(1:6) :: varids
(  294)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  295)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  296)     CHARACTER(LEN=100) :: sec_since
(  297)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  298)     INTEGER :: int_opp, &              ! Operation interval
(  299)                int_wri                 ! Write interval
(  300)     CHARACTER(LEN=3),PARAMETER :: &
(  301)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  302)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  303)     
(  304)     ! Initializations
(  305)     status = NCDF_NOERR
(  306)     CALL ioget_calendar(cal_type)
(  307)     CALL ioget_timestamp(timestamp)
(  308)     WRITE (UNIT=sec_since, &
(  309)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 96

(  310)          &  nyear,nmonth,nday,0, 0, 0
(  311)     WRITE(t_origin, &
(  312)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  313)          &   nyear,months(nmonth),nday,0,0,0
(  314)     int_opp = nwrite * rdt
(  315)     int_wri = nwrite * rdt
(  316) 
(  317)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  318)        WRITE(100,*) 'NCDF DEBUG: Creating default V output file:', filename
(  319)        CALL FLUSH
(  320)     END IF
(  321) 
(  322)     ! Only processor 0 does anything
(  323)     IF(nproc == 0) THEN
(  324)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  325)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Creating file:', filename
(  326)           CALL FLUSH
(  327)        END IF
(  328)        ! Create the file
(  329)        nfstat = nf90_create(filename, nf90_clobber, ncid)
(  330)        IF(nfstat /= nf90_noerr) THEN
(  331)           status = NCDF_NFERR
(  332)           RETURN
(  333)        END IF
(  334)        
(  335)        ! Define dimensions
(  336)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  337)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Defining dimensions in file:', filename
(  338)           CALL FLUSH
(  339)        END IF
(  340)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  341)        nfstat = nf90_def_dim(ncid, 'depthv', jpkdta, dimids(2))
(  342)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
(  343)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
(  344)        
(  345)        ! Define variables
(  346)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  347)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Defining variables in file:', filename
(  348)           CALL FLUSH
(  349)        END IF
(  350)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  351)             (/ dimids(4), dimids(3) /), &
(  352)             varids(1))
(  353)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  354)             (/ dimids(4), dimids(3) /), &
(  355)             varids(2))
(  356)        nfstat = nf90_def_var(ncid, 'depthv', nf90_float, &
(  357)             (/ dimids(2) /), &
(  358)             varids(3))
(  359)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  360)             (/ dimids(1) /), &
(  361)             varids(4))
(  362)        nfstat = nf90_def_var(ncid, 'vomecrty', nf90_float, &
(  363)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  364)             varids(5))
(  365)        nfstat = nf90_def_var(ncid, 'sometauy', nf90_float, &
(  366)             (/ dimids(4), dimids(3), dimids(1) /), &
(  367)             varids(6))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 97

(  368)        
(  369)        ! Add attributes
(  370)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  371)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Writing attributes in file:', filename
(  372)           CALL FLUSH
(  373)        END IF
(  374)        ! nav_lon
(  375)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  376)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  377)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  378)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  379)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  380) 
(  381)        ! nav_lat
(  382)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  383)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  384)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  385)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  386)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  387) 
(  388)        ! depthu
(  389)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  390) 
(  391)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  392)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
(  393)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
(  394)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthv')
(  395)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical V levels')
(  396) 
(  397)        ! time_counter
(  398)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  399)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  400)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  401)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
(  402)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
(  403) 
(  404)        ! vomecrty
(  405)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
(  406)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
(  407)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
(  408)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
(  409)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Meridional Current')
(  410)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vomecrty')
(  411)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
(  412)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
(  413)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
(  414)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
(  415)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthv nav_lat nav_lon')
(  416) 
(  417)        ! sometauy
(  418)        nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
(  419)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
(  420)        nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
(  421)        nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
(  422)        nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along j-axis')
(  423)        nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sometauy')
(  424)        nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
(  425)        nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 98

(  426)        nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
(  427)        nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
(  428)        nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')
(  429) 
(  430)        ! global
(  431)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  432)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
(  433)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
(  434)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  435)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  436)       
(  437)        
(  438)        ! Close file
(  439)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  440)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Closing file:', filename
(  441)           CALL FLUSH
(  442)        END IF
(  443)        nfstat = nf90_close(ncid)
(  444)        IF(nfstat /= nf90_noerr) THEN
(  445)           status = NCDF_NFERR
(  446)           RETURN
(  447)        END IF
(  448)     END IF
(  449) 
(  450)     ! Sync up processors before returning from subroutine
(  451)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
(  452)     IF(mpistat /= 0) THEN
(  453)        status = NCDF_MPERR
(  454)        RETURN
(  455)     END IF
(  456)     
(  457)   END SUBROUTINE ncdf_create_file_v
(  458) 
(  459)   ! ncdf_create_file_t builds a standard t-grid OPA output file with all the default
(  460)   ! dimensions, variables and attributes
(  461)   SUBROUTINE ncdf_create_file_t(filename, op_type, status)
(  462)     IMPLICIT NONE
(  463)     ! Subroutine argument declarations
(  464)     CHARACTER(LEN=*),INTENT(IN) :: filename
(  465)     CHARACTER(LEN=*),INTENT(IN) :: op_type
(  466)     INTEGER,INTENT(OUT) :: status
(  467) 
(  468)     ! Local declarations
(  469)     INTEGER :: ncid,    &  ! netCDF file ID
(  470)                varid,   &  ! ID of netCDF variable to be written to
(  471)                nfstat,  &  ! netCDF library call return status
(  472)                mpistat     ! MPI library call return status
(  473)     INTEGER,DIMENSION(1:4) :: dimids
(  474)     INTEGER,DIMENSION(1:22) :: varids
(  475)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  476)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  477)     CHARACTER(LEN=100) :: sec_since
(  478)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  479)     INTEGER :: int_opp, &              ! Operation interval
(  480)                int_wri                 ! Write interval
(  481)     CHARACTER(LEN=3),PARAMETER :: &
(  482)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  483)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 99

(  484)     
(  485)     ! Initializations
(  486)     status = NCDF_NOERR
(  487)     CALL ioget_calendar(cal_type)
(  488)     CALL ioget_timestamp(timestamp)
(  489)     WRITE (UNIT=sec_since, &
(  490)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  491)          &  nyear,nmonth,nday,0, 0, 0
(  492)     WRITE(t_origin, &
(  493)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  494)          &   nyear,months(nmonth),nday,0,0,0
(  495)     int_opp = nwrite * rdt
(  496)     int_wri = nwrite * rdt
(  497) 
(  498)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  499)        WRITE(100,*) 'NCDF DEBUG: Creating default T output file:', filename
(  500)        CALL FLUSH
(  501)     END IF
(  502) 
(  503)     ! Only processor 0 does anything
(  504)     IF(nproc == 0) THEN
(  505)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  506)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Creating file:', filename
(  507)           CALL FLUSH
(  508)        END IF
(  509)        ! Create the file
(  510)        nfstat = nf90_create(filename, nf90_clobber, ncid)
(  511)        
(  512)        ! Define dimensions
(  513)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  514)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Defining dimensions in file:', filename
(  515)           CALL FLUSH
(  516)        END IF
(  517)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  518)        nfstat = nf90_def_dim(ncid, 'jpkdta', jpkdta, dimids(2))
(  519)        nfstat = nf90_def_dim(ncid, 'jpjdta', jpjdta, dimids(3))
(  520)        nfstat = nf90_def_dim(ncid, 'jpidta', jpidta, dimids(4))
(  521)        
(  522)        ! Define variables
(  523)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  524)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Defining variables in file:', filename
(  525)           CALL FLUSH
(  526)        END IF
(  527)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  528)             (/ dimids(4), dimids(3) /), &
(  529)             varids(1))
(  530)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  531)             (/ dimids(4), dimids(3) /), &
(  532)             varids(2))
(  533)        nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &
(  534)             (/ dimids(2) /), &
(  535)             varids(3))
(  536)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  537)             (/ dimids(1) /), &
(  538)             varids(4))
(  539)        nfstat = nf90_def_var(ncid, 'votemper', nf90_float, &
(  540)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  541)             varids(5))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 100

(  542)        nfstat = nf90_def_var(ncid, 'vosaline', nf90_float, &
(  543)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  544)             varids(6))
(  545)        nfstat = nf90_def_var(ncid, 'sosstsst', nf90_float, &
(  546)             (/ dimids(4), dimids(3), dimids(1) /), &
(  547)             varids(7))
(  548)        nfstat = nf90_def_var(ncid, 'sosaline', nf90_float, &
(  549)             (/ dimids(4), dimids(3), dimids(1) /), &
(  550)             varids(8))
(  551)        nfstat = nf90_def_var(ncid, 'sossheig', nf90_float, &
(  552)             (/ dimids(4), dimids(3), dimids(1) /), &
(  553)             varids(9))
(  554)        nfstat = nf90_def_var(ncid, 'sowaflup', nf90_float, &
(  555)             (/ dimids(4), dimids(3), dimids(1) /), &
(  556)             varids(10))
(  557)        nfstat = nf90_def_var(ncid, 'sorunoff', nf90_float, &
(  558)             (/ dimids(4), dimids(3), dimids(1) /), &
(  559)             varids(11))
(  560)        nfstat = nf90_def_var(ncid, 'sowaflcd', nf90_float, &
(  561)             (/ dimids(4), dimids(3), dimids(1) /), &
(  562)             varids(12))
(  563)        nfstat = nf90_def_var(ncid, 'sosalflx', nf90_float, &
(  564)             (/ dimids(4), dimids(3), dimids(1) /), &
(  565)             varids(13))
(  566)        nfstat = nf90_def_var(ncid, 'sohefldo', nf90_float, &
(  567)             (/ dimids(4), dimids(3), dimids(1) /), &
(  568)             varids(14))
(  569)        nfstat = nf90_def_var(ncid, 'soshfldo', nf90_float, &
(  570)             (/ dimids(4), dimids(3), dimids(1) /), &
(  571)             varids(15))
(  572)        nfstat = nf90_def_var(ncid, 'somxl010', nf90_float, &
(  573)             (/ dimids(4), dimids(3), dimids(1) /), &
(  574)             varids(16))
(  575)        nfstat = nf90_def_var(ncid, 'somixhgt', nf90_float, &
(  576)             (/ dimids(4), dimids(3), dimids(1) /), &
(  577)             varids(17))
(  578)        nfstat = nf90_def_var(ncid, 'soicecov', nf90_float, &
(  579)             (/ dimids(4), dimids(3), dimids(1) /), &
(  580)             varids(18))
(  581)        nfstat = nf90_def_var(ncid, 'sohefldp', nf90_float, &
(  582)             (/ dimids(4), dimids(3), dimids(1) /), &
(  583)             varids(19))
(  584)        nfstat = nf90_def_var(ncid, 'sowafldp', nf90_float, &
(  585)             (/ dimids(4), dimids(3), dimids(1) /), &
(  586)             varids(20))
(  587)        nfstat = nf90_def_var(ncid, 'sosafldp', nf90_float, &
(  588)             (/ dimids(4), dimids(3), dimids(1) /), &
(  589)             varids(21))
(  590)        nfstat = nf90_def_var(ncid, 'sobowlin', nf90_float, &
(  591)             (/ dimids(4), dimids(3) /), &
(  592)             varids(22))
(  593)        
(  594)        ! Add attributes
(  595)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  596)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Writing attributes in file:', filename
(  597)           CALL FLUSH
(  598)        END IF
(  599)        ! nav_lon






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 101

(  600)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  601)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  602)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  603)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  604)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  605) 
(  606)        ! nav_lat
(  607)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  608)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  609)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  610)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  611)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  612) 
(  613)        ! deptht
(  614)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  615)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  616)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
(  617)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
(  618)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
(  619)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical T levels')
(  620) 
(  621)        ! time_counter
(  622)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  623)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  624)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  625)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
(  626)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
(  627) 
(  628)        ! votemper
(  629)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'C')
(  630)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
(  631)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
(  632)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
(  633)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Temperature')
(  634)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'votemper')
(  635)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
(  636)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
(  637)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
(  638)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
(  639)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter deptht nav_lat nav_lon')
(  640) 
(  641)        ! vosaline
(  642)        nfstat = nf90_put_att(ncid, varids(6), 'units', 'PSU')
(  643)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
(  644)        nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
(  645)        nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Salinity')
(  646)        nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'vosaline')
(  647)        nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
(  648)        nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TZYX')
(  649)        nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
(  650)        nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
(  651)        nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter deptht nav_lat nav_lon')
(  652) 
(  653)        ! sosstsst
(  654)        nfstat = nf90_put_att(ncid, varids(7), 'units', 'C')
(  655)        nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
(  656)        nfstat = nf90_put_att(ncid, varids(7), 'valid_min', 1.000000E20 )
(  657)        nfstat = nf90_put_att(ncid, varids(7), 'valid_max', -1.000000E20)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 102

(  658)        nfstat = nf90_put_att(ncid, varids(7), 'long_name', 'Sea Surface temperature')
(  659)        nfstat = nf90_put_att(ncid, varids(7), 'short_name', 'sosstsst')
(  660)        nfstat = nf90_put_att(ncid, varids(7), 'online_operation', TRIM(op_type))
(  661)        nfstat = nf90_put_att(ncid, varids(7), 'axis', 'TYX')
(  662)        nfstat = nf90_put_att(ncid, varids(7), 'interval_operation', int_opp)
(  663)        nfstat = nf90_put_att(ncid, varids(7), 'interval_write', int_wri)
(  664)        nfstat = nf90_put_att(ncid, varids(7), 'associate', 'time_counter nav_lat nav_lon')
(  665) 
(  666)        ! sosaline
(  667)        nfstat = nf90_put_att(ncid, varids(8), 'units', 'PSU')
(  668)        nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
(  669)        nfstat = nf90_put_att(ncid, varids(8), 'valid_min', 1.000000E20 )
(  670)        nfstat = nf90_put_att(ncid, varids(8), 'valid_max', -1.000000E20)
(  671)        nfstat = nf90_put_att(ncid, varids(8), 'long_name', 'Sea Surface Salinity')
(  672)        nfstat = nf90_put_att(ncid, varids(8), 'short_name', 'sosaline')
(  673)        nfstat = nf90_put_att(ncid, varids(8), 'online_operation', TRIM(op_type))
(  674)        nfstat = nf90_put_att(ncid, varids(8), 'axis', 'TYX')
(  675)        nfstat = nf90_put_att(ncid, varids(8), 'interval_operation', int_opp)
(  676)        nfstat = nf90_put_att(ncid, varids(8), 'interval_write', int_wri)
(  677)        nfstat = nf90_put_att(ncid, varids(8), 'associate', 'time_counter nav_lat nav_lon')
(  678) 
(  679)        ! sossheig
(  680)        nfstat = nf90_put_att(ncid, varids(9), 'units', 'm')
(  681)        nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)
(  682)        nfstat = nf90_put_att(ncid, varids(9), 'valid_min', 1.000000E20 )
(  683)        nfstat = nf90_put_att(ncid, varids(9), 'valid_max', -1.000000E20)
(  684)        nfstat = nf90_put_att(ncid, varids(9), 'long_name', 'Sea Surface Height')
(  685)        nfstat = nf90_put_att(ncid, varids(9), 'short_name', 'sossheig')
(  686)        nfstat = nf90_put_att(ncid, varids(9), 'online_operation', TRIM(op_type))
(  687)        nfstat = nf90_put_att(ncid, varids(9), 'axis', 'TYX')
(  688)        nfstat = nf90_put_att(ncid, varids(9), 'interval_operation', int_opp)
(  689)        nfstat = nf90_put_att(ncid, varids(9), 'interval_write', int_wri)
(  690)        nfstat = nf90_put_att(ncid, varids(9), 'associate', 'time_counter nav_lat nav_lon')
(  691)        
(  692)        ! sowaflup
(  693)        nfstat = nf90_put_att(ncid, varids(10), 'units', 'Kg/m2/s')
(  694)        nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)
(  695)        nfstat = nf90_put_att(ncid, varids(10), 'valid_min', 1.000000E20 )
(  696)        nfstat = nf90_put_att(ncid, varids(10), 'valid_max', -1.000000E20)
(  697)        nfstat = nf90_put_att(ncid, varids(10), 'long_name', 'Net Upward Water Flux')
(  698)        nfstat = nf90_put_att(ncid, varids(10), 'short_name', 'sowaflup')
(  699)        nfstat = nf90_put_att(ncid, varids(10), 'online_operation', TRIM(op_type))
(  700)        nfstat = nf90_put_att(ncid, varids(10), 'axis', 'TYX')
(  701)        nfstat = nf90_put_att(ncid, varids(10), 'interval_operation', int_opp)
(  702)        nfstat = nf90_put_att(ncid, varids(10), 'interval_write', int_wri)
(  703)        nfstat = nf90_put_att(ncid, varids(10), 'associate', 'time_counter nav_lat nav_lon')
(  704)        
(  705)        ! sorunoff
(  706)        nfstat = nf90_put_att(ncid, varids(11), 'units', 'Kg/m2/s')
(  707)        nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)
(  708)        nfstat = nf90_put_att(ncid, varids(11), 'valid_min', 1.000000E20 )
(  709)        nfstat = nf90_put_att(ncid, varids(11), 'valid_max', -1.000000E20)
(  710)        nfstat = nf90_put_att(ncid, varids(11), 'long_name', 'Runoffs')
(  711)        nfstat = nf90_put_att(ncid, varids(11), 'short_name', 'sorunoff')
(  712)        nfstat = nf90_put_att(ncid, varids(11), 'online_operation', TRIM(op_type))
(  713)        nfstat = nf90_put_att(ncid, varids(11), 'axis', 'TYX')
(  714)        nfstat = nf90_put_att(ncid, varids(11), 'interval_operation', int_opp)
(  715)        nfstat = nf90_put_att(ncid, varids(11), 'interval_write', int_wri)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 103

(  716)        nfstat = nf90_put_att(ncid, varids(11), 'associate', 'time_counter nav_lat nav_lon')
(  717)        
(  718)        ! sowaflcd
(  719)        nfstat = nf90_put_att(ncid, varids(12), 'units', 'kg/m2/s')
(  720)        nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)
(  721)        nfstat = nf90_put_att(ncid, varids(12), 'valid_min', 1.000000E20 )
(  722)        nfstat = nf90_put_att(ncid, varids(12), 'valid_max', -1.000000E20)
(  723)        nfstat = nf90_put_att(ncid, varids(12), 'long_name', 'concentration/dilution water flux')
(  724)        nfstat = nf90_put_att(ncid, varids(12), 'short_name', 'sowaflcd')
(  725)        nfstat = nf90_put_att(ncid, varids(12), 'online_operation', TRIM(op_type))
(  726)        nfstat = nf90_put_att(ncid, varids(12), 'axis', 'TYX')
(  727)        nfstat = nf90_put_att(ncid, varids(12), 'interval_operation', int_opp)
(  728)        nfstat = nf90_put_att(ncid, varids(12), 'interval_write', int_wri)
(  729)        nfstat = nf90_put_att(ncid, varids(12), 'associate', 'time_counter nav_lat nav_lon')
(  730)        
(  731)        ! sosalflx
(  732)        nfstat = nf90_put_att(ncid, varids(13), 'units', 'Kg/m2/s')
(  733)        nfstat = nf90_put_att(ncid, varids(13), 'missing_value', 1.000000E20)
(  734)        nfstat = nf90_put_att(ncid, varids(13), 'valid_min', 1.000000E20 )
(  735)        nfstat = nf90_put_att(ncid, varids(13), 'valid_max', -1.000000E20)
(  736)        nfstat = nf90_put_att(ncid, varids(13), 'long_name', 'Surface Salt Flux')
(  737)        nfstat = nf90_put_att(ncid, varids(13), 'short_name', 'sosalflx')
(  738)        nfstat = nf90_put_att(ncid, varids(13), 'online_operation', TRIM(op_type))
(  739)        nfstat = nf90_put_att(ncid, varids(13), 'axis', 'TYX')
(  740)        nfstat = nf90_put_att(ncid, varids(13), 'interval_operation', int_opp)
(  741)        nfstat = nf90_put_att(ncid, varids(13), 'interval_write', int_wri)
(  742)        nfstat = nf90_put_att(ncid, varids(13), 'associate', 'time_counter nav_lat nav_lon')
(  743) 
(  744)        ! sohefldo
(  745)        nfstat = nf90_put_att(ncid, varids(14), 'units', 'W/m2')
(  746)        nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)
(  747)        nfstat = nf90_put_att(ncid, varids(14), 'valid_min', 1.000000E20 )
(  748)        nfstat = nf90_put_att(ncid, varids(14), 'valid_max', -1.000000E20)
(  749)        nfstat = nf90_put_att(ncid, varids(14), 'long_name', 'Net Downward Heat Flux')
(  750)        nfstat = nf90_put_att(ncid, varids(14), 'short_name', 'sohefldo')
(  751)        nfstat = nf90_put_att(ncid, varids(14), 'online_operation', TRIM(op_type))
(  752)        nfstat = nf90_put_att(ncid, varids(14), 'axis', 'TYX')
(  753)        nfstat = nf90_put_att(ncid, varids(14), 'interval_operation', int_opp)
(  754)        nfstat = nf90_put_att(ncid, varids(14), 'interval_write', int_wri)
(  755)        nfstat = nf90_put_att(ncid, varids(14), 'associate', 'time_counter nav_lat nav_lon')
(  756) 
(  757)        ! soshfldo
(  758)        nfstat = nf90_put_att(ncid, varids(15), 'units', 'W/m2')
(  759)        nfstat = nf90_put_att(ncid, varids(15), 'missing_value', 1.000000E20)
(  760)        nfstat = nf90_put_att(ncid, varids(15), 'valid_min', 1.000000E20 )
(  761)        nfstat = nf90_put_att(ncid, varids(15), 'valid_max', -1.000000E20)
(  762)        nfstat = nf90_put_att(ncid, varids(15), 'long_name', 'Shortwave Radiation')
(  763)        nfstat = nf90_put_att(ncid, varids(15), 'short_name', 'soshfldo')
(  764)        nfstat = nf90_put_att(ncid, varids(15), 'online_operation', TRIM(op_type))
(  765)        nfstat = nf90_put_att(ncid, varids(15), 'axis', 'TYX')
(  766)        nfstat = nf90_put_att(ncid, varids(15), 'interval_operation', int_opp)
(  767)        nfstat = nf90_put_att(ncid, varids(15), 'interval_write', int_wri)
(  768)        nfstat = nf90_put_att(ncid, varids(15), 'associate', 'time_counter nav_lat nav_lon')
(  769) 
(  770)        ! soml010
(  771)        nfstat = nf90_put_att(ncid, varids(16), 'units', 'm')
(  772)        nfstat = nf90_put_att(ncid, varids(16), 'missing_value', 1.000000E20)
(  773)        nfstat = nf90_put_att(ncid, varids(16), 'valid_min', 1.000000E20 )






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 104

(  774)        nfstat = nf90_put_att(ncid, varids(16), 'valid_max', -1.000000E20)
(  775)        nfstat = nf90_put_att(ncid, varids(16), 'long_name', 'Mixed Layer Depth 0.01')
(  776)        nfstat = nf90_put_att(ncid, varids(16), 'short_name', 'soml010')
(  777)        nfstat = nf90_put_att(ncid, varids(16), 'online_operation', TRIM(op_type))
(  778)        nfstat = nf90_put_att(ncid, varids(16), 'axis', 'TYX')
(  779)        nfstat = nf90_put_att(ncid, varids(16), 'interval_operation', int_opp)
(  780)        nfstat = nf90_put_att(ncid, varids(16), 'interval_write', int_wri)
(  781)        nfstat = nf90_put_att(ncid, varids(16), 'associate', 'time_counter nav_lat nav_lon')
(  782) 
(  783)        ! somixhgt
(  784)        nfstat = nf90_put_att(ncid, varids(17), 'units', 'm')
(  785)        nfstat = nf90_put_att(ncid, varids(17), 'missing_value', 1.000000E20)
(  786)        nfstat = nf90_put_att(ncid, varids(17), 'valid_min', 1.000000E20 )
(  787)        nfstat = nf90_put_att(ncid, varids(17), 'valid_max', -1.000000E20)
(  788)        nfstat = nf90_put_att(ncid, varids(17), 'long_name', 'Turbocline Depth')
(  789)        nfstat = nf90_put_att(ncid, varids(17), 'short_name', 'somixhgt')
(  790)        nfstat = nf90_put_att(ncid, varids(17), 'online_operation', TRIM(op_type))
(  791)        nfstat = nf90_put_att(ncid, varids(17), 'axis', 'TYX')
(  792)        nfstat = nf90_put_att(ncid, varids(17), 'interval_operation', int_opp)
(  793)        nfstat = nf90_put_att(ncid, varids(17), 'interval_write', int_wri)
(  794)        nfstat = nf90_put_att(ncid, varids(17), 'associate', 'time_counter nav_lat nav_lon')
(  795) 
(  796)        ! soicecov
(  797)        nfstat = nf90_put_att(ncid, varids(18), 'units', '[0,1]')
(  798)        nfstat = nf90_put_att(ncid, varids(18), 'missing_value', 1.000000E20)
(  799)        nfstat = nf90_put_att(ncid, varids(18), 'valid_min', 1.000000E20 )
(  800)        nfstat = nf90_put_att(ncid, varids(18), 'valid_max', -1.000000E20)
(  801)        nfstat = nf90_put_att(ncid, varids(18), 'long_name', 'Ice Cover')
(  802)        nfstat = nf90_put_att(ncid, varids(18), 'short_name', 'soicecov')
(  803)        nfstat = nf90_put_att(ncid, varids(18), 'online_operation', TRIM(op_type))
(  804)        nfstat = nf90_put_att(ncid, varids(18), 'axis', 'TYX')
(  805)        nfstat = nf90_put_att(ncid, varids(18), 'interval_operation', int_opp)
(  806)        nfstat = nf90_put_att(ncid, varids(18), 'interval_write', int_wri)
(  807)        nfstat = nf90_put_att(ncid, varids(18), 'associate', 'time_counter nav_lat nav_lon')
(  808) 
(  809)        ! sohefldp
(  810)        nfstat = nf90_put_att(ncid, varids(19), 'units', 'W/m2')
(  811)        nfstat = nf90_put_att(ncid, varids(19), 'missing_value', 1.000000E20)
(  812)        nfstat = nf90_put_att(ncid, varids(19), 'valid_min', 1.000000E20 )
(  813)        nfstat = nf90_put_att(ncid, varids(19), 'valid_max', -1.000000E20)
(  814)        nfstat = nf90_put_att(ncid, varids(19), 'long_name', 'Surface Heat Flux: Damping')
(  815)        nfstat = nf90_put_att(ncid, varids(19), 'short_name', 'sohefldp')
(  816)        nfstat = nf90_put_att(ncid, varids(19), 'online_operation', TRIM(op_type))
(  817)        nfstat = nf90_put_att(ncid, varids(19), 'axis', 'TYX')
(  818)        nfstat = nf90_put_att(ncid, varids(19), 'interval_operation', int_opp)
(  819)        nfstat = nf90_put_att(ncid, varids(19), 'interval_write', int_wri)
(  820)        nfstat = nf90_put_att(ncid, varids(19), 'associate', 'time_counter nav_lat nav_lon')
(  821) 
(  822)        ! sowafldp
(  823)        nfstat = nf90_put_att(ncid, varids(20), 'units', 'Kg/m2/s')
(  824)        nfstat = nf90_put_att(ncid, varids(20), 'missing_value', 1.000000E20)
(  825)        nfstat = nf90_put_att(ncid, varids(20), 'valid_min', 1.000000E20 )
(  826)        nfstat = nf90_put_att(ncid, varids(20), 'valid_max', -1.000000E20)
(  827)        nfstat = nf90_put_att(ncid, varids(20), 'long_name', 'Surface Water Flux: Damping')
(  828)        nfstat = nf90_put_att(ncid, varids(20), 'short_name', 'sowafldp')
(  829)        nfstat = nf90_put_att(ncid, varids(20), 'online_operation', TRIM(op_type))
(  830)        nfstat = nf90_put_att(ncid, varids(20), 'axis', 'TYX')
(  831)        nfstat = nf90_put_att(ncid, varids(20), 'interval_operation', int_opp)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 105

(  832)        nfstat = nf90_put_att(ncid, varids(20), 'interval_write', int_wri)
(  833)        nfstat = nf90_put_att(ncid, varids(20), 'associate', 'time_counter nav_lat nav_lon')
(  834) 
(  835)        ! sosafldp
(  836)        nfstat = nf90_put_att(ncid, varids(21), 'units', 'Kg/m2/s')
(  837)        nfstat = nf90_put_att(ncid, varids(21), 'missing_value', 1.000000E20)
(  838)        nfstat = nf90_put_att(ncid, varids(21), 'valid_min', 1.000000E20 )
(  839)        nfstat = nf90_put_att(ncid, varids(21), 'valid_max', -1.000000E20)
(  840)        nfstat = nf90_put_att(ncid, varids(21), 'long_name', 'Surface salt flux: damping')
(  841)        nfstat = nf90_put_att(ncid, varids(21), 'short_name', 'sosafldp')
(  842)        nfstat = nf90_put_att(ncid, varids(21), 'online_operation', TRIM(op_type))
(  843)        nfstat = nf90_put_att(ncid, varids(21), 'axis', 'TYX')
(  844)        nfstat = nf90_put_att(ncid, varids(21), 'interval_operation', int_opp)
(  845)        nfstat = nf90_put_att(ncid, varids(21), 'interval_write', int_wri)
(  846)        nfstat = nf90_put_att(ncid, varids(21), 'associate', 'time_counter nav_lat nav_lon')
(  847) 
(  848)        ! sobowlin
(  849)        nfstat = nf90_put_att(ncid, varids(22), 'units', 'W-point')
(  850)        nfstat = nf90_put_att(ncid, varids(22), 'missing_value', 1.000000E20)
(  851)        nfstat = nf90_put_att(ncid, varids(22), 'valid_min', 1.000000E20 )
(  852)        nfstat = nf90_put_att(ncid, varids(22), 'valid_max', -1.000000E20)
(  853)        nfstat = nf90_put_att(ncid, varids(22), 'long_name', 'Bowl Index')
(  854)        nfstat = nf90_put_att(ncid, varids(22), 'short_name', 'sobowlin')
(  855)        nfstat = nf90_put_att(ncid, varids(22), 'online_operation', 'l_max(only(x))')
(  856)        nfstat = nf90_put_att(ncid, varids(22), 'axis', 'TYX')
(  857)        nfstat = nf90_put_att(ncid, varids(22), 'associate', 'time_counter nav_lat nav_lon')
(  858) 
(  859)        ! global
(  860)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  861)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
(  862)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
(  863)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  864)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  865)       
(  866)        
(  867)        ! Close file
(  868)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  869)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_t - Closing file:', filename
(  870)           CALL FLUSH
(  871)        END IF
(  872)        nfstat = nf90_close(ncid)
(  873)        IF(nfstat /= nf90_noerr) THEN
(  874)           status = NCDF_NFERR
(  875)           RETURN
(  876)        END IF
(  877)     END IF
(  878) 
(  879)     ! Sync up processors before returning from subroutine
(  880)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
(  881)     IF(mpistat /= 0) THEN
(  882)        status = NCDF_MPERR
(  883)        RETURN
(  884)     END IF
(  885)     
(  886)   END SUBROUTINE ncdf_create_file_t
(  887) 
(  888)   ! ncdf_create_file_ave builds an OPA output file with a single variable
(  889)   ! whose name is specified by VARNAME. Other than that, it is similar to a






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 106

(  890)   ! U- or V-grid file
(  891)  SUBROUTINE ncdf_create_file_ave(filename, varname, op_type, status)
(  892)     IMPLICIT NONE
(  893)     ! Subroutine argument declarations
(  894)     CHARACTER(LEN=*),INTENT(IN) :: filename
(  895)     CHARACTER(LEN=*),INTENT(IN) :: op_type, varname
(  896)     INTEGER,INTENT(OUT) :: status
(  897) 
(  898)     ! Local declarations
(  899)     INTEGER :: ncid,    &  ! netCDF file ID
(  900)                varid,   &  ! ID of netCDF variable to be written to
(  901)                nfstat,  &  ! netCDF library call return status
(  902)                mpistat     ! MPI library call return status
(  903)     INTEGER,DIMENSION(1:4) :: dimids
(  904)     INTEGER,DIMENSION(1:5) :: varids
(  905)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  906)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  907)     CHARACTER(LEN=100) :: sec_since    
(  908)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  909)     INTEGER :: int_opp, &              ! Operation interval
(  910)                int_wri                 ! Write interval
(  911)     INTEGER :: varnum
(  912)     CHARACTER(LEN=3),PARAMETER :: &
(  913)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  914)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  915)     
(  916)     ! Initializations
(  917)     status = NCDF_NOERR
(  918)     CALL ioget_calendar(cal_type)
(  919)     CALL ioget_timestamp(timestamp)
(  920)     WRITE (UNIT=sec_since, &
(  921)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  922)          &  nyear,nmonth,nday,0, 0, 0
(  923)     WRITE(t_origin, &
(  924)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  925)          &   nyear,months(nmonth),nday,0,0,0
(  926)     int_opp = nwrite * rdt
(  927)     int_wri = nwrite * rdt
(  928) 
(  929)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  930)        WRITE(100,*) 'NCDF DEBUG: Creating output file:', filename
(  931)        CALL FLUSH
(  932)     END IF
(  933) 
(  934)     ! Only processor 0 does anything
(  935)     IF(nproc == 0) THEN
(  936)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  937)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Creating file:', filename
(  938)           CALL FLUSH
(  939)        END IF
(  940)        ! Create the file
(  941)        nfstat = nf90_create(filename, nf90_clobber, ncid)
(  942)        IF(nfstat /= nf90_noerr) THEN
(  943)           status = NCDF_NFERR
(  944)           RETURN
(  945)        END IF
(  946)        
(  947)        ! Define dimensions






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 107

(  948)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  949)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Defining dimensions in file:', filename
(  950)           CALL FLUSH
(  951)        END IF
(  952)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  953)        nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
(  954)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
(  955)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
(  956)        
(  957)        ! Define variables
(  958)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  959)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Defining variables in file:', filename
(  960)           CALL FLUSH
(  961)        END IF
(  962)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  963)             (/ dimids(4), dimids(3) /), &
(  964)             varids(1))
(  965)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  966)             (/ dimids(4), dimids(3) /), &
(  967)             varids(2))
(  968)        nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
(  969)             (/ dimids(2) /), &
(  970)             varids(3))
(  971)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  972)             (/ dimids(1) /), &
(  973)             varids(4))
(  974)        nfstat = nf90_def_var(ncid, varname, nf90_float, &
(  975)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  976)             varids(5))
(  977)        
(  978)        ! Add attributes
(  979)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  980)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Writing attributes in file:', filename
(  981)           CALL FLUSH
(  982)        END IF
(  983)        ! nav_lon
(  984)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  985)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  986)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  987)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  988)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  989) 
(  990)        ! nav_lat
(  991)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  992)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  993)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  994)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  995)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  996) 
(  997)        ! depthu
(  998)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  999)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
( 1000)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
( 1001)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
( 1002)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
( 1003)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')
( 1004) 
( 1005)        ! time_counter






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 108

( 1006)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
( 1007)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
( 1008)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
( 1009)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
( 1010)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
( 1011) 
( 1012)        ! custom variable
( 1013)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
( 1014)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
( 1015)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
( 1016)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
( 1017)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
( 1018)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', varname)
( 1019)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
( 1020)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
( 1021)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
( 1022)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')
( 1023) 
( 1024)        ! global
( 1025)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
( 1026)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
( 1027)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
( 1028)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
( 1029)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
( 1030) 
( 1031)        varnum = 5
( 1032) !AD/DB: add new time-related variables
( 1033)        varnum = varnum + 1
( 1034)        ! ndate (ndastp)
( 1035)        nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
( 1036)             (/ dimids(1) /), &
( 1037)             varids(varnum))
( 1038)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
( 1039)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'time step date in year/month/day aammjj')
( 1040) 
( 1041)        varnum = varnum + 1
( 1042)        ! ndate (model_time)
( 1043)        nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
( 1044)             (/ dimids(1) /), &
( 1045)             varids(varnum))
( 1046)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', &
( 1047)             'time step date (when output is writen) in year/month/day aammjj (decimal day)')
( 1048)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
( 1049)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula1', 'nyear  =   model_time / 10000')       
( 1050)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula2', & 
( 1051)             'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
( 1052)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula3', & 
( 1053)             'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           
( 1054) 
( 1055)        varnum = varnum + 1
( 1056)        ! kt 
( 1057)        nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
( 1058)             (/ dimids(1) /), &
( 1059)             varids(varnum))
( 1060) 
( 1061)        
( 1062)        ! Close file
( 1063)        IF(DEBUG_OUT .EQV. .TRUE.) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 109

( 1064)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Closing file:', filename
( 1065)           CALL FLUSH
( 1066)        END IF
( 1067)        nfstat = nf90_close(ncid)
( 1068)        IF(nfstat /= nf90_noerr) THEN
( 1069)           status = NCDF_NFERR
( 1070)           RETURN
( 1071)        END IF
( 1072)     END IF
( 1073) 
( 1074)     ! Sync up processors before returning from subroutine
( 1075)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1076)     IF(mpistat /= 0) THEN
( 1077)        status = NCDF_MPERR
( 1078)        RETURN
( 1079)     END IF
( 1080)     
( 1081)   END SUBROUTINE ncdf_create_file_ave
( 1082)   
( 1083)   ! ncdf_create_restart builds a standard OPA restart file with all the default
( 1084)   ! dimensions, variables and attributes. Note that this is a single restart
( 1085)   ! file, and therefore incompatible with the default IOIPSL-based restart
( 1086)   ! routines (though it holds the same data)
( 1087)   ! NOTE: Some keys we haven't been using are untested, so I don't know if those
( 1088)   ! fields will be created correctly.
( 1089)   SUBROUTINE ncdf_create_restart(filename, status)
( 1090)     IMPLICIT NONE
( 1091)     ! Subroutine argument declarations
( 1092)     CHARACTER(LEN=*),INTENT(IN) :: filename
( 1093)     INTEGER,INTENT(OUT) :: status
( 1094)     
( 1095)     ! Local declarations
( 1096)     INTEGER :: ncid,    &  ! netCDF file ID
( 1097)          varid,   &  ! ID of netCDF variable to be written to
( 1098)          nfstat,  &  ! netCDF library call return status
( 1099)          mpistat     ! MPI library call return status
( 1100)     INTEGER,DIMENSION(1:8) :: dimids
( 1101)     INTEGER,DIMENSION(1:36) :: varids
( 1102)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
( 1103)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
( 1104)     CHARACTER(LEN=100) :: sec_since
( 1105)     CHARACTER(LEN=100) :: tstp_since
( 1106)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
( 1107)     CHARACTER(LEN=3),PARAMETER :: &
( 1108)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
( 1109)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
( 1110)     
( 1111)     ! Initializations
( 1112)     status = NCDF_NOERR
( 1113)     
( 1114)     ! Initializations
( 1115)     status = NCDF_NOERR
( 1116)     CALL ioget_calendar(cal_type)
( 1117)     CALL ioget_timestamp(timestamp)
( 1118)     WRITE (UNIT=sec_since, &
( 1119)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 1120)          &  nyear,nmonth,nday,0, 0, 0
( 1121)     WRITE (UNIT=tstp_since, &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 110

( 1122)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 1123)          &  nyear,nmonth,nday,0, 0, 0
( 1124)     WRITE(t_origin, &
( 1125)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
( 1126)          &   nyear,months(nmonth),nday,0,0,0
( 1127) 
( 1128)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1129)        WRITE(100,*) 'NCDF DEBUG: Creating default restart file:', filename
( 1130)        CALL FLUSH
( 1131)     END IF
( 1132)     
( 1133)     ! Only processor 0 does anything
( 1134)     IF(nproc == 0) THEN
( 1135)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1136)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Creating file:', filename
( 1137)           CALL FLUSH
( 1138)        END IF
( 1139)        ! Create the file
( 1140)        nfstat = nf90_create(filename, nf90_clobber, ncid)
( 1141)        
( 1142)        ! Define dimensions
( 1143)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1144)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining dimensions in file:', filename
( 1145)           CALL FLUSH
( 1146)        END IF
( 1147)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(1))
( 1148)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(2))
( 1149)        nfstat = nf90_def_dim(ncid, 'z', jpkdta, dimids(3))
( 1150)        nfstat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimids(4))
( 1151)        nfstat = nf90_def_dim(ncid, 'x_a', 1, dimids(5))
( 1152)        nfstat = nf90_def_dim(ncid, 'y_a', 1, dimids(6))
( 1153)        nfstat = nf90_def_dim(ncid, 'z_a', 10, dimids(7))
( 1154)        nfstat = nf90_def_dim(ncid, 'z_b', 1, dimids(8))
( 1155)        
( 1156)        ! Define variables
( 1157)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1158)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining variables in file:', filename
( 1159)           CALL FLUSH
( 1160)        END IF
( 1161)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
( 1162)             (/ dimids(1), dimids(2) /), &
( 1163)             varids(1))
( 1164)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
( 1165)             (/ dimids(1), dimids(2) /), &
( 1166)             varids(2))
( 1167)        nfstat = nf90_def_var(ncid, 'nav_lev', nf90_float, &
( 1168)             (/ dimids(3) /), &
( 1169)             varids(3))
( 1170)        nfstat = nf90_def_var(ncid, 'time', nf90_float, &
( 1171)             (/ dimids(4) /), &
( 1172)             varids(4))
( 1173)        nfstat = nf90_def_var(ncid, 'time_steps', nf90_int, &
( 1174)             (/ dimids(4) /), &
( 1175)             varids(5))
( 1176)        nfstat = nf90_def_var(ncid, 'info', nf90_double, &
( 1177)             (/ dimids(7) /), &
( 1178)             varids(6))
( 1179)        nfstat = nf90_def_var(ncid, 'ub', nf90_double, &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 111

( 1180)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1181)             varids(7))
( 1182)        nfstat = nf90_def_var(ncid, 'vb', nf90_double, &
( 1183)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1184)             varids(8))
( 1185)        nfstat = nf90_def_var(ncid, 'tb', nf90_double, &
( 1186)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1187)             varids(9))
( 1188)        nfstat = nf90_def_var(ncid, 'sb', nf90_double, &
( 1189)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1190)             varids(10))
( 1191)        nfstat = nf90_def_var(ncid, 'rotb', nf90_double, &
( 1192)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1193)             varids(11))
( 1194)        nfstat = nf90_def_var(ncid, 'hdivb', nf90_double, &
( 1195)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1196)             varids(12))
( 1197)        nfstat = nf90_def_var(ncid, 'un', nf90_double, &
( 1198)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1199)             varids(13))
( 1200)        nfstat = nf90_def_var(ncid, 'vn', nf90_double, &
( 1201)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1202)             varids(14))
( 1203)        nfstat = nf90_def_var(ncid, 'tn', nf90_double, &
( 1204)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1205)             varids(15))
( 1206)        nfstat = nf90_def_var(ncid, 'sn', nf90_double, &
( 1207)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1208)             varids(16))
( 1209)        nfstat = nf90_def_var(ncid, 'rotn', nf90_double, &
( 1210)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1211)             varids(17))
( 1212)        nfstat = nf90_def_var(ncid, 'hdivn', nf90_double, &
( 1213)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1214)             varids(18))
( 1215)        nfstat = nf90_def_var(ncid, 'gcx', nf90_double, &
( 1216)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1217)             varids(19))
( 1218)        nfstat = nf90_def_var(ncid, 'gcxb', nf90_double, &
( 1219)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1220)             varids(20))
( 1221) # 1232
( 1232)        nfstat = nf90_def_var(ncid, 'sshb', nf90_double, &
( 1233)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1234)             varids(21))
( 1235)        nfstat = nf90_def_var(ncid, 'sshn', nf90_double, &
( 1236)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1237)             varids(22))
( 1238) # 1239
( 1239)        nfstat = nf90_def_var(ncid, 'sshb_b', nf90_double, &
( 1240)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1241)             varids(23))
( 1242)        nfstat = nf90_def_var(ncid, 'sshn_b', nf90_double, &
( 1243)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1244)             varids(24))
( 1245)        nfstat = nf90_def_var(ncid, 'un_b', nf90_double, &
( 1246)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1247)             varids(25))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 112

( 1248)        nfstat = nf90_def_var(ncid, 'vn_b', nf90_double, &
( 1249)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1250)             varids(26))
( 1251) 
( 1252) # 1254
( 1254)        ! Fields that are only defined if specific keys are set
( 1255) # 1256
( 1256)        nfstat = nf90_def_var(ncid, 'en', nf90_double, &
( 1257)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 1258)             varids(28))
( 1259) # 1261
( 1261)        nfstat = nf90_def_var(ncid, 'nfice', nf90_double, &
( 1262)             (/dimids(8) /), &
( 1263)             varids(29))
( 1264)        nfstat = nf90_def_var(ncid, 'sst_io', nf90_double, &
( 1265)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1266)             varids(30))
( 1267)        nfstat = nf90_def_var(ncid, 'sss_io', nf90_double, &
( 1268)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1269)             varids(31))
( 1270)        nfstat = nf90_def_var(ncid, 'u_io', nf90_double, &
( 1271)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1272)             varids(32))
( 1273)        nfstat = nf90_def_var(ncid, 'v_io', nf90_double, &
( 1274)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1275)             varids(33))
( 1276) # 1283
( 1283)        nfstat = nf90_def_var(ncid, 'nfbulk', nf90_double, &
( 1284)             (/ dimids(8) /), &
( 1285)             varids(35))
( 1286) 
( 1287)        nfstat = nf90_def_var(ncid, 'gsst', nf90_double, &
( 1288)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 1289)             varids(36))
( 1290) # 1291
( 1291)        
( 1292)        ! Add attributes
( 1293)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1294)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Writing attributes in file:', filename
( 1295)           CALL FLUSH
( 1296)        END IF
( 1297)        ! nav_lon
( 1298)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
( 1299)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -1.800000E2)
( 1300)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', 1.800000E2)
( 1301)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
( 1302) 
( 1303)        ! nav_lat
( 1304)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
( 1305)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', -9.000000E1)
( 1306)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 9.000000E1)
( 1307)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
( 1308) 
( 1309)        ! nav_lev
( 1310)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'model_levels')
( 1311)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E0)
( 1312)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E3)
( 1313)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Model levels')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 113

( 1314) 
( 1315)        ! time
( 1316)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
( 1317)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
( 1318)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
( 1319)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'Time axis')
( 1320)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', '0001-JUL-01 00:00:00')
( 1321) 
( 1322)        ! time_steps
( 1323)        nfstat = nf90_put_att(ncid, varids(5), 'units', TRIM(tstp_since))
( 1324)        nfstat = nf90_put_att(ncid, varids(5), 'title', 'Time steps')
( 1325)        nfstat = nf90_put_att(ncid, varids(5), 'tstep_sec', 1.877760E7)
( 1326)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Time step axis')
( 1327)        nfstat = nf90_put_att(ncid, varids(5), 'time_origin', TRIM(t_origin))
( 1328) 
( 1329)        ! info
( 1330)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
( 1331) 
( 1332)        ! ub
( 1333)        nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
( 1334) 
( 1335)        ! vb
( 1336)        nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
( 1337) 
( 1338)        ! tb
( 1339)        nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)
( 1340) 
( 1341)        ! sb
( 1342)        nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)
( 1343) 
( 1344)        ! rotb
( 1345)        nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)
( 1346) 
( 1347)        ! hdivb
( 1348)        nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)
( 1349) 
( 1350)        ! un
( 1351)        nfstat = nf90_put_att(ncid, varids(13), 'missing_value', 1.000000E20)
( 1352) 
( 1353)        ! vn
( 1354)        nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)
( 1355) 
( 1356)        ! tn
( 1357)        nfstat = nf90_put_att(ncid, varids(15), 'missing_value', 1.000000E20)
( 1358) 
( 1359)        ! sn
( 1360)        nfstat = nf90_put_att(ncid, varids(16), 'missing_value', 1.000000E20)
( 1361) 
( 1362)        ! rotn
( 1363)        nfstat = nf90_put_att(ncid, varids(17), 'missing_value', 1.000000E20)
( 1364) 
( 1365)        ! hdivn
( 1366)        nfstat = nf90_put_att(ncid, varids(18), 'missing_value', 1.000000E20)
( 1367) 
( 1368)        ! gcx
( 1369)        nfstat = nf90_put_att(ncid, varids(19), 'missing_value', 1.000000E20)
( 1370) 
( 1371)        ! gcxb






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 114

( 1372)        nfstat = nf90_put_att(ncid, varids(20), 'missing_value', 1.000000E20)
( 1373) # 1383
( 1383)        ! sshb
( 1384)        nfstat = nf90_put_att(ncid, varids(21), 'missing_value', 1.000000E20)
( 1385) 
( 1386)        ! sshn
( 1387)        nfstat = nf90_put_att(ncid, varids(22), 'missing_value', 1.000000E20)
( 1388) 
( 1389) # 1390
( 1390)        ! sshb_b
( 1391)        nfstat = nf90_put_att(ncid, varids(23), 'missing_value', 1.000000E20)
( 1392) 
( 1393)        ! sshn_b
( 1394)        nfstat = nf90_put_att(ncid, varids(24), 'missing_value', 1.000000E20)
( 1395) 
( 1396)        ! un_b
( 1397)        nfstat = nf90_put_att(ncid, varids(25), 'missing_value', 1.000000E20)
( 1398) 
( 1399)        ! vn_b
( 1400)        nfstat = nf90_put_att(ncid, varids(26), 'missing_value', 1.000000E20)
( 1401) # 1404
( 1404)        ! en
( 1405)        nfstat = nf90_put_att(ncid, varids(28), 'missing_value', 1.000000E20)
( 1406) # 1408
( 1408)        ! nfice
( 1409)        nfstat = nf90_put_att(ncid, varids(29), 'missing_value', 1.000000E20)
( 1410) 
( 1411)        ! sst_io
( 1412)        nfstat = nf90_put_att(ncid, varids(30), 'missing_value', 1.000000E20)
( 1413) 
( 1414)        ! sss_io
( 1415)        nfstat = nf90_put_att(ncid, varids(31), 'missing_value', 1.000000E20)
( 1416) 
( 1417)        ! u_io
( 1418)        nfstat = nf90_put_att(ncid, varids(32), 'missing_value', 1.000000E20)
( 1419) 
( 1420)        ! v_io
( 1421)        nfstat = nf90_put_att(ncid, varids(33), 'missing_value', 1.000000E20)
( 1422) 
( 1423) # 1430
( 1430)        ! nfbulk
( 1431)        nfstat = nf90_put_att(ncid, varids(35), 'missing_value', 1.000000E20)
( 1432)  
( 1433)        ! gsst
( 1434)        nfstat = nf90_put_att(ncid, varids(36), 'missing_value', 1.000000E20)
( 1435)  
( 1436) 
( 1437) # 1438
( 1438)        ! global
( 1439)        nfstat = nf90_put_att(ncid, NF90_GLOBAL , 'Conventions', 'GDT 1.2')
( 1440)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
( 1441)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
( 1442)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number_total', 0)
( 1443)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number', 0)
( 1444)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids', 0)
( 1445)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_global', 0)
( 1446)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_local', 0)
( 1447)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_first', 0)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 115

( 1448)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_last', 0)
( 1449)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', 0)
( 1450)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', 0)
( 1451)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_type','box')
( 1452) 
( 1453)        ! Close file
( 1454)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1455)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Closing file:', filename
( 1456)           CALL FLUSH
( 1457)        END IF
( 1458)        nfstat = nf90_close(ncid)
( 1459)        IF(nfstat /= nf90_noerr) THEN
( 1460)           status = NCDF_NFERR
( 1461)           RETURN
( 1462)        END IF
( 1463)     END IF
( 1464) 
( 1465)     ! Sync up processors before returning from subroutine
( 1466)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1467)     IF(mpistat /= 0) THEN
( 1468)        status = NCDF_MPERR
( 1469)        RETURN
( 1470)     END IF
( 1471)        
( 1472)   END SUBROUTINE ncdf_create_restart
( 1473) 
( 1474)   ! ncdf_writesv (single value) writes a single, scalar value to a specified
( 1475)   ! index in a netCDF variable. The netCDF variable is assumed to be
( 1476)   ! 1-dimensional.
( 1477)   ! filename - file to write to
( 1478)   ! varname - variable to write to
( 1479)   ! data - scalar value to write
( 1480)   ! index - where in the variable to put the data
( 1481)   ! status - return status of the subroutine
( 1482)   SUBROUTINE ncdf_writesv(filename, varname, data, index, status)
( 1483)      IMPLICIT NONE
( 1484)     ! Subroutine argument declarations
( 1485)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 1486)     REAL(wp),INTENT(IN) :: data
( 1487)     INTEGER,INTENT(IN) :: index   ! Where in the variable to write the value
( 1488)     INTEGER,INTENT(OUT) :: status
( 1489) 
( 1490)     ! Local declarations
( 1491)     INTEGER :: ncid,    &  ! netCDF file ID
( 1492)                varid,   &  ! ID of netCDF variable to be written to
( 1493)                nfstat,  &  ! netCDF library call return status
( 1494)                mpistat     ! MPI library call return status
( 1495) 
( 1496)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1497)        WRITE(100,*) 'NCDF DEBUG: Writing single value to file:', filename
( 1498)        CALL FLUSH
( 1499)     END IF
( 1500) 
( 1501)     ! Initializations
( 1502)     status = NCDF_NOERR
( 1503) 
( 1504)     ! Open netCDF file and get info
( 1505)     IF(nproc == 0) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 116

( 1506)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1507)           WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Opening file:', filename
( 1508)           CALL FLUSH
( 1509)        END IF
( 1510)        nfstat = nf90_open(filename, nf90_write, ncid)
( 1511)        IF(nfstat /= nf90_noerr) THEN
( 1512)           status = NCDF_NFERR
( 1513)           RETURN
( 1514)        END IF
( 1515)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1516)           WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Getting info from file:', filename
( 1517)           CALL FLUSH
( 1518)        END IF
( 1519)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 1520)        IF(nfstat /= nf90_noerr) THEN
( 1521)           status = NCDF_NFERR
( 1522)           RETURN
( 1523)        END IF
( 1524)        
( 1525)        ! Write data to netCDF file
( 1526)        ! This subroutine assumes all processors have the same data in this
( 1527)        ! variable, no attempt to sync data or merge data in the file is made
( 1528)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1529)           WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Writing data to file:', filename
( 1530)           CALL FLUSH
( 1531)        END IF
( 1532)        nfstat = nf90_put_var(ncid, varid, data, &
( 1533)            (/ index /))
( 1534)        IF(nfstat /= nf90_noerr) THEN
( 1535)           status = NCDF_NFERR
( 1536)           RETURN
( 1537)        END IF
( 1538)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1539)           WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Closing file:', filename
( 1540)           CALL FLUSH
( 1541)        END IF
( 1542)        nfstat = nf90_close(ncid)
( 1543)        IF(nfstat /= nf90_noerr) THEN
( 1544)           status = NCDF_NFERR
( 1545)           RETURN
( 1546)        END IF
( 1547)     END IF
( 1548) 
( 1549)     ! Sync up processors before returning from subroutine
( 1550)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1551)     IF(mpistat /= 0) THEN
( 1552)        status = NCDF_MPERR
( 1553)        RETURN
( 1554)     END IF
( 1555)     
( 1556)   END SUBROUTINE ncdf_writesv
( 1557) 
( 1558)   ! ncdf_write1d writes a 1-dimensional array to a 1-D variable in a netCDF
( 1559)   ! file. It assumes that all processors have the same data in the array, 
( 1560)   ! so no attempt is made to sync data between processors or merge data
( 1561)   ! together in the file, unlike ncdf_write2d and ncdf_write3d
( 1562)   ! filename - file to write to
( 1563)   ! varname - variable to write to






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 117

( 1564)   ! data - 1-D array to write
( 1565)   ! status - return status of the subroutine
( 1566)   SUBROUTINE ncdf_write1d(filename, varname, data, status)
( 1567)     IMPLICIT NONE
( 1568)     ! Subroutine argument declarations
( 1569)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 1570)     REAL(wp),DIMENSION(:),INTENT(IN) :: data
( 1571)     INTEGER,INTENT(OUT) :: status
( 1572) 
( 1573)     ! Local declarations
( 1574)     INTEGER :: ncid,    &  ! netCDF file ID
( 1575)                varid,   &  ! ID of netCDF variable to be written to
( 1576)                nfstat,  &  ! netCDF library call return status
( 1577)                mpistat     ! MPI library call return status
( 1578) 
( 1579)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1580)        WRITE(100,*) 'NCDF DEBUG: Writing 1D array to file:', filename
( 1581)        CALL FLUSH
( 1582)     END IF
( 1583) 
( 1584)     ! Initializations
( 1585)     status = NCDF_NOERR
( 1586) 
( 1587)     ! Open netCDF file and get info
( 1588)     IF(nproc == 0) THEN
( 1589)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1590)           WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Opening file:', filename
( 1591)           CALL FLUSH
( 1592)        END IF
( 1593)        nfstat = nf90_open(filename, nf90_write, ncid)
( 1594)        IF(nfstat /= nf90_noerr) THEN
( 1595)           status = NCDF_NFERR
( 1596)           RETURN
( 1597)        END IF
( 1598)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1599)           WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Getting info from file:', filename
( 1600)           CALL FLUSH
( 1601)        END IF
( 1602)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 1603)        IF(nfstat /= nf90_noerr) THEN
( 1604)           status = NCDF_NFERR
( 1605)           RETURN
( 1606)        END IF
( 1607)        
( 1608)        ! Write data to netCDF file
( 1609)        ! This subroutine assumes all processors have the same data in this
( 1610)        ! array, no attempt to sync data or merge data in the file is made
( 1611)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1612)           WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Writing data to file:', filename
( 1613)           CALL FLUSH
( 1614)        END IF
( 1615)        nfstat = nf90_put_var(ncid, varid, data)
( 1616)        IF(nfstat /= nf90_noerr) THEN
( 1617)           status = NCDF_NFERR
( 1618)           WRITE(125,*) 'ERROR: ', nf90_strerror(nfstat)
( 1619)           CALL FLUSH
( 1620)           RETURN
( 1621)        END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 118

( 1622)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1623)           WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Closing file:', filename
( 1624)           CALL FLUSH
( 1625)        END IF
( 1626)        nfstat = nf90_close(ncid)
( 1627)        IF(nfstat /= nf90_noerr) THEN
( 1628)           status = NCDF_NFERR
( 1629)           RETURN
( 1630)        END IF
( 1631)     END IF
( 1632) 
( 1633)     ! Sync up processors before returning from subroutine
( 1634)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1635)     IF(mpistat /= 0) THEN
( 1636)        status = NCDF_MPERR
( 1637)        RETURN
( 1638)     END IF
( 1639) 
( 1640)   END SUBROUTINE ncdf_write1d
( 1641) 
( 1642)   ! ncdf_write2d writes a 2-dimensional array to a 2-D variable in a netCDF
( 1643)   ! file. It assumes the data array being passed is a subdomain - that is,
( 1644)   ! each processor's array contains different values, and the data needs to
( 1645)   ! be written to the location in the netCDF file which corresponds to the
( 1646)   ! processor writing the subdomain. This subroutine will correctly merge
( 1647)   ! all subroutines so that the variable in the netCDF file contains
( 1648)   ! the entire global domain of values
( 1649)   ! filename - file to write to
( 1650)   ! varname - variable to write to
( 1651)   ! data - 2-D array to write
( 1652)   ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
( 1653)   !                                      will force an immediate write, with
( 1654)   !                                      data being written at a time index in
( 1655)   !                                      the file equal to (tstep * -1). Use this
( 1656)   !                                      to output on timesteps which are not
( 1657)   !                                      multiples of nwrite
( 1658)   ! status - return status of the subroutine
( 1659)   ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
( 1660)   ! a tstep which is <0
( 1661)   SUBROUTINE ncdf_write2d(filename, varname, data, tstep, status)
( 1662)     IMPLICIT NONE
( 1663)     ! Subroutine argument declarations
( 1664)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 1665)     INTEGER,INTENT(IN) :: tstep
( 1666)     REAL(wp),DIMENSION(:,:),INTENT(IN) :: data
( 1667)     INTEGER,INTENT(OUT) :: status
( 1668) 
( 1669)     ! Local declarations
( 1670)     INTEGER :: ncid,    &  ! netCDF file ID
( 1671)                varid,   &  ! ID of netCDF variable to be written to
( 1672)                tindex,  &  ! Index to write to on the time axis
( 1673)                mpistat, &  ! MPI library call return status
( 1674)                nfstat,  &  ! netCDF library call return status
( 1675)                is,      &  ! i-axis start index for writing to file
( 1676)                js,      &  ! j-axis   "     "    "     "     "  "
( 1677)                i,k,       &  ! Loop counter
( 1678)                ndims,   &  ! Number of dimensions in this variabl
( 1679)                iz,      &  ! I-size of array to write to file






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 119

( 1680)                jz          ! J-size of array to write to file
( 1681)     
( 1682)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 1683)     INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
( 1684)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 1685)     REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: buf   ! Send/receive buffer for data array
( 1686)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 1687) 
( 1688)     ! Exit subroutine if no write should take place on current timestep
( 1689)     IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >= 0)) THEN
( 1690)        status = NCDF_NOERR
( 1691)        RETURN
( 1692)     END IF
( 1693) 
( 1694)     ! Initializations
( 1695)     is = nimpp                          ! Get global i-index for this subdomain
( 1696)     js = njmpp                          ! Get global j-index  "   "      "
( 1697)     iz = nlei                           ! Get i-size of subdomain section to be written
( 1698)     jz = nlej                           ! Get j-size of subdomain section to be written
( 1699)     ALLOCATE(buf(1:iz,1:jz))            ! Allocate send/receive buffer
( 1700)     buf = data(1:iz,1:jz)               ! Copy section of subdomain to buffer
( 1701)     IF(tstep >= 0) THEN
( 1702)        tindex = tstep/nwrite - nit000/nwrite
( 1703)     ELSE
( 1704)        tindex = tstep * (-1)
( 1705)     END IF
( 1706)     status = NCDF_NOERR
( 1707)     hastime = .FALSE.
( 1708) 
( 1709)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1710)        WRITE(100,*) 'NCDF DEBUG: Writing 2D array to file:', filename
( 1711)        CALL FLUSH
( 1712)     END IF
( 1713) 
( 1714)     ! Open netCDF file and get info
( 1715)     IF(nproc == 0) THEN
( 1716)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1717)           WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Opening file:', filename
( 1718)           CALL FLUSH
( 1719)        END IF
( 1720)        nfstat = nf90_open(filename, nf90_write, ncid)
( 1721)        IF(nfstat /= nf90_noerr) THEN
( 1722)           status = NCDF_NFERR
( 1723)           RETURN
( 1724)        END IF
( 1725)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1726)           WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Getting info from file:', filename
( 1727)           CALL FLUSH
( 1728)        END IF
( 1729)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 1730)        IF(nfstat /= nf90_noerr) THEN
( 1731)           status = NCDF_NFERR
( 1732)           RETURN
( 1733)        END IF
( 1734)        
( 1735)        ! Determine if this variable contains a time axis
( 1736)        nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 1737)             ndims=ndims, dimids=var_dimids)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 120

( 1738)        IF(nfstat /= nf90_noerr) THEN
( 1739)           status = NCDF_NFERR
( 1740)           RETURN
( 1741)        END IF
( 1742)        DO i = 1, ndims
( 1743)           nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 1744)                name=dname, len=k)
( 1745) !          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 1746) !               name=dname)
( 1747)           dname=TRIM(dname)
( 1748)           IF(nfstat /= nf90_noerr) THEN
( 1749)              status = NCDF_NFERR
( 1750)              RETURN
( 1751)           END IF
( 1752)           IF(dname=='time_counter') THEN
( 1753)              hastime = .TRUE.
( 1754)           END IF
( 1755)        END DO
( 1756)     END IF
( 1757) 
( 1758)     ! Get data from each processor in turn, processor 0 writes to file
( 1759)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1760)        WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Writing data to file:', filename
( 1761)        CALL FLUSH
( 1762)     END IF
( 1763) 
( 1764)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1765)     IF(mpistat /= 0) THEN
( 1766)        status = NCDF_MPERR
( 1767)        RETURN
( 1768)     END IF
( 1769)     DO i = 0, (jpnij - 1)
( 1770)        ! Since processor 0 doesn't need to send data to itself, only
( 1771)        ! processors with IDs > 0 make send calls
( 1772)        IF(i > 0) THEN
( 1773)           ! If it's the local processor's turn to send, go ahead and
( 1774)           ! send data and write indices to processor 0
( 1775)           IF(nproc == i) THEN
( 1776)              CALL MPI_SEND(is, 1, MPI_INTEGER, &
( 1777)                   0, i, MPI_COMM_WORLD, mpistat)
( 1778)              IF(mpistat /= 0) THEN
( 1779)                 status = NCDF_MPERR
( 1780)                 RETURN
( 1781)              END IF
( 1782)              CALL MPI_SEND(js, 1, MPI_INTEGER, &
( 1783)                   0, i+1, MPI_COMM_WORLD, mpistat)
( 1784)              IF(mpistat /= 0) THEN
( 1785)                 status = NCDF_MPERR
( 1786)                 RETURN
( 1787)              END IF
( 1788)              CALL MPI_SEND(iz, 1, MPI_INTEGER, &
( 1789)                   0, i+2, MPI_COMM_WORLD, mpistat)
( 1790)              IF(mpistat /= 0) THEN
( 1791)                 status = NCDF_MPERR
( 1792)                 RETURN
( 1793)              END IF
( 1794)              CALL MPI_SEND(jz, 1, MPI_INTEGER, &
( 1795)                   0, i+3, MPI_COMM_WORLD, mpistat)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 121

( 1796)              IF(mpistat /= 0) THEN
( 1797)                 status = NCDF_MPERR
( 1798)                 RETURN
( 1799)              END IF
( 1800)              CALL MPI_SEND(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 1801)                   0, i+4, MPI_COMM_WORLD, mpistat)
( 1802)              IF(mpistat /= 0) THEN
( 1803)                 status = NCDF_MPERR
( 1804)                 RETURN
( 1805)              END IF
( 1806)           END IF
( 1807)        END IF
( 1808) 
( 1809)        ! Processor 0 receives data and writes it to file
( 1810)        IF(nproc == 0) THEN
( 1811)           ! Since processor 0 doesn't need to receive data from itself, 
( 1812)           ! we only make receive calls if it's another processor's turn
( 1813)           ! to send data
( 1814)           IF(i > 0) THEN
( 1815)              CALL MPI_RECV(is, 1, MPI_INTEGER, &
( 1816)                   i, i, MPI_COMM_WORLD, rcvstat, mpistat)
( 1817)              IF(mpistat /= 0) THEN
( 1818)                 status = NCDF_MPERR
( 1819)                 RETURN
( 1820)              END IF
( 1821)              CALL MPI_RECV(js, 1, MPI_INTEGER, &
( 1822)                   i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)
( 1823)              IF(mpistat /= 0) THEN
( 1824)                 status = NCDF_MPERR
( 1825)                 RETURN
( 1826)              END IF
( 1827)              CALL MPI_RECV(iz, 1, MPI_INTEGER, &
( 1828)                   i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
( 1829)              IF(mpistat /= 0) THEN
( 1830)                 status = NCDF_MPERR
( 1831)                 RETURN
( 1832)              END IF
( 1833)              CALL MPI_RECV(jz, 1, MPI_INTEGER, &
( 1834)                   i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
( 1835)              IF(mpistat /= 0) THEN
( 1836)                 status = NCDF_MPERR
( 1837)                 RETURN
( 1838)              END IF
( 1839)              DEALLOCATE(buf)
( 1840)              ALLOCATE(buf(1:iz,1:jz))
( 1841)              CALL MPI_RECV(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 1842)                   i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
( 1843)              IF(mpistat /= 0) THEN
( 1844)                 status = NCDF_MPERR
( 1845)                 RETURN
( 1846)              END IF
( 1847)           END IF
( 1848)           ! Got the data, write it to the proper location in the netCDF file
( 1849)           IF(hastime .EQV. .TRUE.) THEN
( 1850)              nfstat = nf90_put_var(ncid, varid, buf, &
( 1851)                   (/ is, js, tindex/), &
( 1852)                   (/ iz, jz, 1 /))
( 1853)           ELSE






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 122

( 1854)              nfstat = nf90_put_var(ncid, varid, buf, &
( 1855)                   (/ is, js /), &
( 1856)                   (/ iz, jz /))
( 1857)           END IF
( 1858)           IF(nfstat /= nf90_noerr) THEN
( 1859)              status = NCDF_NFERR
( 1860)              RETURN
( 1861)           END IF
( 1862)        END IF
( 1863) 
( 1864)        ! Sync all processors at the end of each loop iteration to prevent
( 1865)        ! concurrency issues/race conditions/various MPI-related ugliness
( 1866)        CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1867)        IF(mpistat /= 0) THEN
( 1868)           status = NCDF_MPERR
( 1869)           RETURN
( 1870)        END IF
( 1871)     END DO
( 1872) 
( 1873)     ! All done, close up the netCDF file
( 1874)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1875)        WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Closing file:', filename
( 1876)        CALL FLUSH
( 1877)     END IF
( 1878)     IF(nproc == 0) THEN
( 1879)        nfstat = nf90_close(ncid)
( 1880)        IF(nfstat /= nf90_noerr) THEN
( 1881)           status = NCDF_NFERR
( 1882)           RETURN
( 1883)        END IF
( 1884)     END IF
( 1885)     DEALLOCATE(buf)
( 1886) 
( 1887)     ! Sync up processors before returning from subroutine
( 1888)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1889)     IF(mpistat /= 0) THEN
( 1890)        status = NCDF_MPERR
( 1891)        RETURN
( 1892)     END IF
( 1893) 
( 1894)   END SUBROUTINE ncdf_write2d
( 1895) 
( 1896)   ! ncdf_write3d writes a 3-dimensional array to a 3-D variable in a netCDF
( 1897)   ! file. It assumes the data array being passed is a subdomain - that is,
( 1898)   ! each processor's array contains different values, and the data needs to
( 1899)   ! be written to the location in the netCDF file which corresponds to the
( 1900)   ! processor writing the subdomain. This subroutine will correctly merge
( 1901)   ! all subroutines so that the variable in the netCDF file contains
( 1902)   ! the entire global domain of values
( 1903)   ! filename - file to write to
( 1904)   ! varname - variable to write to
( 1905)   ! data - 3-D array to write
( 1906)   ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
( 1907)   !                                      will force an immediate write, with
( 1908)   !                                      data being written at a time index in
( 1909)   !                                      the file equal to (tstep * -1). Use this
( 1910)   !                                      to output on timesteps which are not
( 1911)   !                                      multiples of nwrite






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 123

( 1912)   ! status - return status of the subroutine
( 1913)   ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
( 1914)   ! a tstep which is <0
( 1915)   SUBROUTINE ncdf_write3d(filename, varname, data, tstep, status)
( 1916)     IMPLICIT NONE
( 1917)     ! Subroutine argument declarations
( 1918)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 1919)     INTEGER,INTENT(IN) :: tstep
( 1920)     REAL(wp),DIMENSION(:,:,:),INTENT(IN) :: data
( 1921)     INTEGER,INTENT(OUT) :: status
( 1922) 
( 1923)     ! Local declarations
( 1924)     INTEGER :: ncid,    &  ! netCDF file ID
( 1925)                varid,   &  ! ID of netCDF variable to be written to
( 1926)                tindex,  &  ! Index to write to on the time axis
( 1927)                mpistat, &  ! MPI library call return status
( 1928)                nfstat,  &  ! netCDF library call return status
( 1929)                !rcvstat, &  ! MPI receive status
( 1930)                is,      &  ! i-axis start index for writing to file
( 1931)                js,      &  ! j-axis   "     "    "     "     "  "
( 1932)                ks,      &  ! k-axis   "     "    "     "     "  "
( 1933)                i,k,     &  ! Loop counter
( 1934)                ndims,   &  ! Number of dimensions in this variable
( 1935)                iz,      &  ! I-size of array to write to file
( 1936)                jz          ! J-size of array to write to file
( 1937) 
( 1938)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 1939)     INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
( 1940)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: buf   ! Send/receive buffer for data array
( 1941)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 1942)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 1943) 
( 1944)     ! Exit subroutine if no write should take place on current timestep
( 1945)     IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >=0)) THEN
( 1946)        status = NCDF_NOERR
( 1947)        RETURN
( 1948)     END IF
( 1949) 
( 1950)     ! Initializations
( 1951)     is = nimpp                          ! Get global i-index for this subdomain
( 1952)     js = njmpp                          ! Get global j-index  "   "      "
( 1953)     ks = 1                              ! Depth writes always start at k=1
( 1954)     iz = nlei                           ! Get i-size of subdomain section to be written
( 1955)     jz = nlej                           ! Get j-size of subdomain section to be written
( 1956)     ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,DIM=3)))      ! Allocate send/recieve buffer
( 1957)     buf = data(1:iz,1:jz,1:SIZE(data,DIM=3))         ! Copy section of subdomain to buffer
( 1958)     IF(tstep >= 0) THEN
( 1959)        tindex = tstep/nwrite - nit000/nwrite
( 1960)     ELSE
( 1961)        tindex = tstep * (-1)
( 1962)     END IF
( 1963)     status = NCDF_NOERR
( 1964)     hastime = .FALSE.
( 1965) 
( 1966)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 1967)        WRITE(100,*) 'NCDF DEBUG: Writing 3D array to file:', filename
( 1968)        CALL FLUSH
( 1969)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 124

( 1970) 
( 1971)     ! Open netCDF file and get info
( 1972)     IF(nproc == 0) THEN
( 1973)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1974)           WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Opening file:', filename
( 1975)           CALL FLUSH
( 1976)        END IF
( 1977)        nfstat = nf90_open(filename, nf90_write, ncid)
( 1978)        IF(nfstat /= nf90_noerr) THEN
( 1979)           status = NCDF_NFERR
( 1980)           RETURN
( 1981)        END IF
( 1982)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 1983)           WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Getting info from file:', filename
( 1984)           CALL FLUSH
( 1985)        END IF
( 1986)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 1987)        IF(nfstat /= nf90_noerr) THEN
( 1988)           status = NCDF_NFERR
( 1989)           RETURN
( 1990)        END IF
( 1991) 
( 1992)        ! Determine if this variable contains a time axis
( 1993)        nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 1994)             ndims=ndims, dimids=var_dimids)
( 1995)        IF(nfstat /= nf90_noerr) THEN
( 1996)           status = NCDF_NFERR
( 1997)           RETURN
( 1998)        END IF
( 1999)        DO i = 1, ndims
( 2000)           nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2001)                name=dname, len=k)
( 2002) !          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2003) !               name=dname)
( 2004)           dname=TRIM(dname)
( 2005)           IF(nfstat /= nf90_noerr) THEN
( 2006)              status = NCDF_NFERR
( 2007)              RETURN
( 2008)           END IF
( 2009)           IF(dname=='time_counter') THEN
( 2010)              hastime = .TRUE.
( 2011)              EXIT
( 2012)           END IF
( 2013)        END DO
( 2014)     END IF
( 2015) 
( 2016)     ! Get data from each processor in turn, processor 0 writes to file
( 2017)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2018)        WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Writing data to file:', filename
( 2019)        CALL FLUSH
( 2020)     END IF
( 2021)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2022)     IF(mpistat /= 0) THEN
( 2023)        status = NCDF_MPERR
( 2024)        RETURN
( 2025)     END IF
( 2026)     DO i=0, (jpnij - 1)
( 2027)        ! Since processor 0 doesn't need to send data to itself, only






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 125

( 2028)        ! processors with IDs > 0 make send calls
( 2029)        IF(i > 0) THEN
( 2030)           ! If it's the local processor's turn to send, go ahead and
( 2031)           ! send data and write indices to processor 0
( 2032)           IF(nproc == i) THEN
( 2033)              CALL MPI_SEND(is, 1, MPI_INTEGER, &
( 2034)                   0, i, MPI_COMM_WORLD, mpistat)
( 2035)              IF(mpistat /= 0) THEN
( 2036)                 status = NCDF_MPERR
( 2037)                 RETURN
( 2038)              END IF
( 2039)              CALL MPI_SEND(js, 1, MPI_INTEGER, &
( 2040)                   0, i+1, MPI_COMM_WORLD, mpistat)
( 2041)              IF(mpistat /= 0) THEN
( 2042)                 status = NCDF_MPERR
( 2043)                 RETURN
( 2044)              END IF
( 2045)              CALL MPI_SEND(ks, 1, MPI_INTEGER, &
( 2046)                   0, i+2, MPI_COMM_WORLD, mpistat)
( 2047)              IF(mpistat /= 0) THEN
( 2048)                 status = NCDF_MPERR
( 2049)                 RETURN
( 2050)              END IF
( 2051)              CALL MPI_SEND(iz, 1, MPI_INTEGER, &
( 2052)                   0, i+3, MPI_COMM_WORLD, mpistat)
( 2053)              IF(mpistat /= 0) THEN
( 2054)                 status = NCDF_MPERR
( 2055)                 RETURN
( 2056)              END IF
( 2057)              CALL MPI_SEND(jz, 1, MPI_INTEGER, &
( 2058)                   0, i+4, MPI_COMM_WORLD, mpistat)
( 2059)              IF(mpistat /= 0) THEN
( 2060)                 status = NCDF_MPERR
( 2061)                 RETURN
( 2062)              END IF
( 2063)              CALL MPI_SEND(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 2064)                   0, i+5, MPI_COMM_WORLD, mpistat)
( 2065)              IF(mpistat /= 0) THEN
( 2066)                 status = NCDF_MPERR
( 2067)                 RETURN
( 2068)              END IF
( 2069)           END IF
( 2070)        END IF
( 2071) 
( 2072)        ! Processor 0 receives data and writes it to file
( 2073)        IF(nproc == 0) THEN
( 2074)           ! Since processor 0 doesn't need to receive data from itself, 
( 2075)           ! we only make receive calls if it's another processor's turn
( 2076)           ! to send data
( 2077)           IF(i > 0) THEN
( 2078)              CALL MPI_RECV(is, 1, MPI_INTEGER, &
( 2079)                   i, i, MPI_COMM_WORLD, rcvstat, mpistat)
( 2080)              IF(mpistat /= 0) THEN
( 2081)                 status = NCDF_MPERR
( 2082)                 RETURN
( 2083)              END IF
( 2084)              CALL MPI_RECV(js, 1, MPI_INTEGER, &
( 2085)                   i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 126

( 2086)              IF(mpistat /= 0) THEN
( 2087)                 status = NCDF_MPERR
( 2088)                 RETURN
( 2089)              END IF
( 2090)              CALL MPI_RECV(ks, 1, MPI_INTEGER, &
( 2091)                   i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
( 2092)              IF(mpistat /= 0) THEN
( 2093)                 status = NCDF_MPERR
( 2094)                 RETURN
( 2095)              END IF
( 2096)              CALL MPI_RECV(iz, 1, MPI_INTEGER, &
( 2097)                   i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
( 2098)              IF(mpistat /= 0) THEN
( 2099)                 status = NCDF_MPERR
( 2100)                 RETURN
( 2101)              END IF
( 2102)              CALL MPI_RECV(jz, 1, MPI_INTEGER, &
( 2103)                   i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
( 2104)              IF(mpistat /= 0) THEN
( 2105)                 status = NCDF_MPERR
( 2106)                 RETURN
( 2107)              END IF
( 2108)              DEALLOCATE(buf)
( 2109)              ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,DIM=3)))
( 2110)              CALL MPI_RECV(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 2111)                   i, i+5, MPI_COMM_WORLD, rcvstat, mpistat)
( 2112)              IF(mpistat /= 0) THEN
( 2113)                 status = NCDF_MPERR
( 2114)                 RETURN
( 2115)              END IF
( 2116)           END IF
( 2117)           
( 2118)           ! Got the data, write it to the proper location in the netCDF file
( 2119)           IF(hastime .EQV. .TRUE.) THEN
( 2120)              nfstat = nf90_put_var(ncid, varid, buf, &
( 2121)                   (/ is, js, ks, tindex/), &
( 2122)                   (/ iz, jz, SIZE(data,DIM=3), 1 /))
( 2123)           ELSE
( 2124)              nfstat = nf90_put_var(ncid, varid, buf, &
( 2125)                   (/ is, js, ks /), &
( 2126)                   (/ iz, jz, SIZE(data,DIM=3) /))
( 2127)           END IF
( 2128)           IF(nfstat /= nf90_noerr) THEN
( 2129)              status = NCDF_NFERR
( 2130)              RETURN
( 2131)           END IF
( 2132)        END IF
( 2133)        
( 2134)        ! Sync all processors at the end of each loop iteration to prevent
( 2135)        ! concurrency issues/race conditions/various MPI-related ugliness
( 2136)        CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2137)        IF(mpistat /= 0) THEN
( 2138)           status = NCDF_MPERR
( 2139)           RETURN
( 2140)        END IF
( 2141)     END DO
( 2142) 
( 2143)     ! All done, close up the netCDF file






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 127

( 2144)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2145)        WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Closing file:', filename
( 2146)        CALL FLUSH
( 2147)     END IF
( 2148)     IF(nproc == 0) THEN
( 2149)        nfstat = nf90_close(ncid)
( 2150)        IF(nfstat /= nf90_noerr) THEN
( 2151)           status = NCDF_NFERR
( 2152)           RETURN
( 2153)        END IF
( 2154)     END IF
( 2155)     DEALLOCATE(buf)
( 2156) 
( 2157)     ! Sync up processors before returning from subroutine
( 2158)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2159)     IF(mpistat /= 0) THEN
( 2160)        status = NCDF_MPERR
( 2161)        RETURN
( 2162)     END IF
( 2163)     
( 2164)   END SUBROUTINE ncdf_write3d
( 2165)   
( 2166)   ! ncdf_write4d writes a 4-dimensional array to a 4-D variable in a netCDF
( 2167)   ! file. It assumes the data array being passed is a subdomain - that is,
( 2168)   ! each processor's array contains different values, and the data needs to
( 2169)   ! be written to the location in the netCDF file which corresponds to the
( 2170)   ! processor writing the subdomain. This subroutine will correctly merge
( 2171)   ! all subroutines so that the variable in the netCDF file contains
( 2172)   ! the entire global domain of values
( 2173)   ! filename - file to write to
( 2174)   ! varname - variable to write to
( 2175)   ! data - 3-D array to write
( 2176)   ! dsz - Dimension SiZe - size of the 4th dimension in the variable
( 2177)   ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
( 2178)   !                                      will force an immediate write, with
( 2179)   !                                      data being written at a time index in
( 2180)   !                                      the file equal to (tstep * -1). Use this
( 2181)   !                                      to output on timesteps which are not
( 2182)   !                                      multiples of nwrite
( 2183)   ! status - return status of the subroutine
( 2184)   ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
( 2185)   ! a tstep which is < 0
( 2186)   ! NOT EXTENSIVELY TESTED - MAY STILL BE BUGGY
( 2187)   SUBROUTINE ncdf_write4d(filename, varname, data, dsz, tstep, status)
( 2188)     IMPLICIT NONE
( 2189)     ! Subroutine argument declarations
( 2190)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2191)     INTEGER,INTENT(IN) :: tstep, dsz
( 2192)     REAL,DIMENSION(:,:,:,:),INTENT(IN) :: data
( 2193)     INTEGER,INTENT(OUT) :: status
( 2194) 
( 2195)     ! Local declarations
( 2196)     INTEGER :: ncid,    &  ! netCDF file ID
( 2197)                varid,   &  ! ID of netCDF variable to be written to
( 2198)                tindex,  &  ! Index to write to on the time axis
( 2199)                mpistat, &  ! MPI library call return status
( 2200)                nfstat,  &  ! netCDF library call return status
( 2201)                is,      &  ! i-axis start index for writing to file






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 128

( 2202)                js,      &  ! j-axis   "     "    "     "     "  "
( 2203)                ks,      &  ! k-axis   "     "    "     "     "  "
( 2204)                i,k,       &  ! Loop counter
( 2205)                ndims,   &  ! Number of dimensions in this variable
( 2206)                iz,      &  ! I-size of array to write to file
( 2207)                jz          ! J-size of array to write to file
( 2208) 
( 2209)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 2210)     INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
( 2211)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:,:) :: buf   ! Send/receive buffer for data array
( 2212)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 2213)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 2214) 
( 2215)     ! Exit subroutine if no write should take place on current timestep
( 2216)     IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >= 0)) THEN
( 2217)        status = NCDF_NOERR
( 2218)        RETURN
( 2219)     END IF
( 2220) 
( 2221)     ! Initializations
( 2222)     is = nimpp                          ! Get global i-index for this subdomain
( 2223)     js = njmpp                          ! Get global j-index  "   "      "
( 2224)     ks = 1                              ! Depth writes always start at k=1
( 2225)     is = nimpp                          ! Get global i-index for this subdomain
( 2226)     js = njmpp                          ! Get global j-index  "   "      "
( 2227)     iz = nlei                           ! Get i-size of subdomain section to be written
( 2228)     jz = nlej                           ! Get j-size of subdomain section to be written
( 2229)     ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,3),1:SIZE(data,4)))      ! Allocate send/recieve buffer
( 2230)     buf = data(1:iz,1:jz,1:SIZE(data,3), 1:SIZE(data,4))         ! Copy section of subdomain to buffer
( 2231)     IF(tstep >= 0) THEN
( 2232)        tindex = tstep/nwrite - nit000/nwrite
( 2233)     ELSE
( 2234)        tindex = tstep * (-1)
( 2235)     END IF
( 2236)     status = NCDF_NOERR
( 2237)     hastime = .FALSE.
( 2238) 
( 2239)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2240)        WRITE(100,*) 'NCDF DEBUG: Writing 4D array to file:', filename
( 2241)     END IF
( 2242) 
( 2243)     ! Open netCDF file and get info
( 2244)     IF(nproc == 0) THEN
( 2245)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 2246)           WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Opening file:', filename
( 2247)           CALL FLUSH
( 2248)        END IF
( 2249)        nfstat = nf90_open(filename, nf90_write, ncid)
( 2250)        IF(nfstat /= nf90_noerr) THEN
( 2251)           status = NCDF_NFERR
( 2252)           RETURN
( 2253)        END IF
( 2254)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 2255)           WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Getting info from file:', filename
( 2256)           CALL FLUSH
( 2257)        END IF
( 2258)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 2259)        IF(nfstat /= nf90_noerr) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 129

( 2260)           status = NCDF_NFERR
( 2261)           RETURN
( 2262)        END IF
( 2263) 
( 2264)        ! Determine if this variable contains a time axis
( 2265)        nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 2266)             ndims=ndims, dimids=var_dimids)
( 2267)        IF(nfstat /= nf90_noerr) THEN
( 2268)           status = NCDF_NFERR
( 2269)           RETURN
( 2270)        END IF
( 2271)        DO i = 1, ndims
( 2272)           nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2273)                name=dname, len=k)
( 2274) !!DB: next line not correct 
( 2275) !          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2276) !               name=dname)
( 2277)           dname=TRIM(dname)
( 2278)           IF(nfstat /= nf90_noerr) THEN
( 2279)              status = NCDF_NFERR
( 2280)              RETURN
( 2281)           END IF
( 2282)           IF(dname=='time_counter') THEN
( 2283)              hastime = .TRUE.
( 2284)              EXIT
( 2285)           END IF
( 2286)        END DO
( 2287)     END IF
( 2288) 
( 2289)     ! Get data from each processor in turn, processor 0 writes to file
( 2290)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2291)        WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Writing data to file:', filename
( 2292)        CALL FLUSH
( 2293)     END IF
( 2294)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2295)     IF(mpistat /= 0) THEN
( 2296)        status = NCDF_MPERR
( 2297)        RETURN
( 2298)     END IF
( 2299)    
( 2300)     DO i=0, (jpnij - 1)
( 2301)        ! Since processor 0 doesn't need to send data to itself, only
( 2302)        ! processors with IDs > 0 make send calls
( 2303)        IF(i > 0) THEN
( 2304)           ! If it's the local processor's turn to send, go ahead and
( 2305)           ! send data and write indices to processor 0
( 2306)           IF(nproc == i) THEN
( 2307)              CALL MPI_SEND(is, 1, MPI_INTEGER, &
( 2308)                   0, i, MPI_COMM_WORLD, mpistat)
( 2309)              IF(mpistat /= 0) THEN
( 2310)                 status = NCDF_MPERR
( 2311)                 RETURN
( 2312)              END IF
( 2313)              CALL MPI_SEND(js, 1, MPI_INTEGER, &
( 2314)                   0, i+1, MPI_COMM_WORLD, mpistat)
( 2315)              IF(mpistat /= 0) THEN
( 2316)                 status = NCDF_MPERR
( 2317)                 RETURN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 130

( 2318)              END IF
( 2319)              CALL MPI_SEND(ks, 1, MPI_INTEGER, &
( 2320)                   0, i+2, MPI_COMM_WORLD, mpistat)
( 2321)              IF(mpistat /= 0) THEN
( 2322)                 status = NCDF_MPERR
( 2323)                 RETURN
( 2324)              END IF
( 2325)              CALL MPI_SEND(iz, 1, MPI_INTEGER, &
( 2326)                   0, i+3, MPI_COMM_WORLD, mpistat)
( 2327)              IF(mpistat /= 0) THEN
( 2328)                 status = NCDF_MPERR
( 2329)                 RETURN
( 2330)              END IF
( 2331)              CALL MPI_SEND(jz, 1, MPI_INTEGER, &
( 2332)                   0, i+4, MPI_COMM_WORLD, mpistat)
( 2333)              IF(mpistat /= 0) THEN
( 2334)                 status = NCDF_MPERR
( 2335)                 RETURN
( 2336)              END IF
( 2337)              CALL MPI_SEND(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 2338)                   0, i+5, MPI_COMM_WORLD, mpistat)
( 2339)              IF(mpistat /= 0) THEN
( 2340)                 status = NCDF_MPERR
( 2341)                 RETURN
( 2342)              END IF
( 2343)           END IF
( 2344)        END IF
( 2345) 
( 2346)        ! Processor 0 receives data and writes it to file
( 2347)        IF(nproc == 0) THEN
( 2348)           ! Since processor 0 doesn't need to receive data from itself, 
( 2349)           ! we only make receive calls if it's another processor's turn
( 2350)           ! to send data
( 2351)           IF(i > 0) THEN
( 2352)              CALL MPI_RECV(is, 1, MPI_INTEGER, &
( 2353)                   i, i, MPI_COMM_WORLD, rcvstat, mpistat)
( 2354)              IF(mpistat /= 0) THEN
( 2355)                 status = NCDF_MPERR
( 2356)                 RETURN
( 2357)              END IF
( 2358)              CALL MPI_RECV(js, 1, MPI_INTEGER, &
( 2359)                   i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)
( 2360)              IF(mpistat /= 0) THEN
( 2361)                 status = NCDF_MPERR
( 2362)                 RETURN
( 2363)              END IF
( 2364)              CALL MPI_RECV(ks, 1, MPI_INTEGER, &
( 2365)                   i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
( 2366)              IF(mpistat /= 0) THEN
( 2367)                 status = NCDF_MPERR
( 2368)                 RETURN
( 2369)              END IF
( 2370)              CALL MPI_RECV(iz, 1, MPI_INTEGER, &
( 2371)                   i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
( 2372)              IF(mpistat /= 0) THEN
( 2373)                 status = NCDF_MPERR
( 2374)                 RETURN
( 2375)              END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 131

( 2376)              CALL MPI_RECV(jz, 1, MPI_INTEGER, &
( 2377)                   i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
( 2378)              IF(mpistat /= 0) THEN
( 2379)                 status = NCDF_MPERR
( 2380)                 RETURN
( 2381)              END IF
( 2382)              DEALLOCATE(buf)
( 2383)              ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,3),1:SIZE(data,4)))
( 2384)              CALL MPI_RECV(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
( 2385)                   i, i+5, MPI_COMM_WORLD, rcvstat, mpistat)
( 2386)              IF(mpistat /= 0) THEN
( 2387)                 status = NCDF_MPERR
( 2388)                 RETURN
( 2389)              END IF
( 2390)           END IF
( 2391)           ! Got the data, write it to the proper location in the netCDF file
( 2392)           IF(hastime .EQV. .TRUE.) THEN
( 2393)              nfstat = nf90_put_var(ncid, varid, buf, &
( 2394)                   (/ is, js, ks, 1, tindex/), &
( 2395)                   (/ iz, jz, SIZE(data,3), SIZE(data,4), 1 /))
( 2396)           ELSE
( 2397)              nfstat = nf90_put_var(ncid, varid, buf, &
( 2398)                   (/ is, js, ks, 1 /), &
( 2399)                   (/ iz, jz, SIZE(data,3), SIZE(data,4) /))
( 2400)           END IF
( 2401)           IF(nfstat /= nf90_noerr) THEN
( 2402)              status = NCDF_NFERR
( 2403)              RETURN
( 2404)           END IF
( 2405)        END IF
( 2406) 
( 2407)        ! Sync all processors at the end of each loop iteration to prevent
( 2408)        ! concurrency issues/race conditions/various MPI-related ugliness
( 2409)        CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2410)        IF(mpistat /= 0) THEN
( 2411)           status = NCDF_MPERR
( 2412)           RETURN
( 2413)        END IF
( 2414)     END DO
( 2415) 
( 2416)     ! All done, close up the netCDF file
( 2417)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2418)        WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Closing file:', filename
( 2419)        CALL FLUSH
( 2420)     END IF
( 2421)     IF(nproc == 0) THEN
( 2422)        nfstat = nf90_close(ncid)
( 2423)        IF(nfstat /= nf90_noerr) THEN
( 2424)           status = NCDF_NFERR
( 2425)           RETURN
( 2426)        END IF
( 2427)     END IF
( 2428)     DEALLOCATE(buf)
( 2429) 
( 2430)     ! Sync up processors before returning from subroutine
( 2431)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2432)     IF(mpistat /= 0) THEN
( 2433)        status = NCDF_MPERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 132

( 2434)        RETURN
( 2435)     END IF
( 2436) 
( 2437)   END SUBROUTINE ncdf_write4d
( 2438) 
( 2439)   ! The inverse of ncdf_writesv. Reads a scalar value from a 1-D array in a
( 2440)   ! netCDF file. All processors will get the same value.
( 2441)   ! filename - file to read from
( 2442)   ! varname - variable to read from
( 2443)   ! data - scalar to put the data into (must be a REAL)
( 2444)   ! index -the index to read from
( 2445)   ! status - return status of the subroutine
( 2446)   SUBROUTINE ncdf_readsv(filename, varname, data, index, status)
( 2447)     IMPLICIT NONE
( 2448)     ! Subroutine argument declarations
( 2449)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2450)     INTEGER,INTENT(IN) :: index
( 2451)     REAL,INTENT(OUT) :: data
( 2452)     INTEGER,INTENT(OUT) :: status
( 2453) 
( 2454)     ! Local declarations
( 2455)     INTEGER :: ncid,    &  ! netCDF file ID
( 2456)                varid,   &  ! ID of netCDF variable to be written to
( 2457)                nfstat,  &  ! netCDF library call return status
( 2458)                mpistat     ! MPI library call return status
( 2459)     INTEGER,DIMENSION(1:5) :: count
( 2460) 
( 2461)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2462)        WRITE(100,*) 'NCDF DEBUG: Reading single value from file:', filename
( 2463)        CALL FLUSH
( 2464)     END IF
( 2465) 
( 2466)     ! Initializations
( 2467)     status = NCDF_NOERR
( 2468)     count(1) = 1
( 2469) 
( 2470)     ! Open netCDF file and get info
( 2471) !!DB: think that there should be no if(nproc==0) then ...
( 2472) !    IF(nproc == 0) THEN
( 2473)        IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2474)           WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Opening file:', filename
( 2475)           CALL FLUSH
( 2476)        END IF
( 2477) !DB
( 2478) !       nfstat = nf90_open(filename, nf90_write, ncid)
( 2479)        nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2480)        IF(nfstat /= nf90_noerr) THEN
( 2481)           status = NCDF_NFERR
( 2482)           RETURN
( 2483)        END IF
( 2484)        IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2485)           WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Getting info from file:', filename
( 2486)           CALL FLUSH
( 2487)        END IF
( 2488)        nfstat = nf90_inq_varid(ncid, varname, varid)
( 2489)        IF(nfstat /= nf90_noerr) THEN
( 2490)           status = NCDF_NFERR
( 2491)           RETURN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 133

( 2492)        END IF
( 2493)        
( 2494)        ! Read data from netCDF file
( 2495)        ! This subroutine assumes all processors will get the same data from this
( 2496)        ! variable
( 2497)        IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2498)           WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Reading data from file:', filename
( 2499)           CALL FLUSH
( 2500)        END IF
( 2501)        nfstat = nf90_get_var(ncid, varid, data, (/ index /))
( 2502)        IF(nfstat /= nf90_noerr) THEN
( 2503)           status = NCDF_NFERR
( 2504)           RETURN
( 2505)        END IF
( 2506)        IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2507)           WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Closing file:', filename
( 2508)           CALL FLUSH
( 2509)        END IF
( 2510)        nfstat = nf90_close(ncid)
( 2511)        IF(nfstat /= nf90_noerr) THEN
( 2512)           status = NCDF_NFERR
( 2513)           RETURN
( 2514)        END IF
( 2515) !    END IF
( 2516) 
( 2517)     ! Sync up processors before returning from subroutine
( 2518)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2519)     IF(mpistat /= 0) THEN
( 2520)        status = NCDF_MPERR
( 2521)        RETURN
( 2522)     END IF
( 2523) 
( 2524)   END SUBROUTINE ncdf_readsv
( 2525) 
( 2526)   ! The inverse of ncdf_write1d. Reads a 1-D array from a netCDF file
( 2527)   ! All processors will get the same values.
( 2528)   ! filename - file to read from
( 2529)   ! varname - variable to read from
( 2530)   ! data - 1-D array to put the data into (must be a REAL)
( 2531)   ! status - return status of the subroutine
( 2532)   SUBROUTINE ncdf_read1d(filename, varname, data, status, tide)
( 2533)     IMPLICIT NONE
( 2534)     ! Subroutine argument declarations
( 2535)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2536)     REAL,DIMENSION(:),INTENT(OUT) :: data
( 2537)     INTEGER,INTENT(OUT) :: status
( 2538)     INTEGER,INTENT(IN),OPTIONAL :: tide
( 2539) 
( 2540)     ! Local declarations
( 2541)     INTEGER :: ncid,    &  ! netCDF file ID
( 2542)                varid,   &  ! ID of netCDF variable to be written to
( 2543)                nfstat,  &  ! netCDF library call return status
( 2544)                mpistat     ! MPI library call return status
( 2545)     
( 2546) 
( 2547) 
( 2548)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2549)        WRITE(100,*) 'NCDF DEBUG: Reading 1D array from file:', filename






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 134

( 2550)        CALL FLUSH
( 2551)     END IF
( 2552) 
( 2553)     ! Initializations
( 2554)     status = NCDF_NOERR
( 2555) 
( 2556)     ! Open netCDF file and get info
( 2557)     
( 2558)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2559)        WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Opening file:', filename
( 2560)        CALL FLUSH
( 2561)     END IF
( 2562)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2563)     IF(nfstat /= nf90_noerr) THEN
( 2564)        status = NCDF_NFERR
( 2565)        RETURN
( 2566)     END IF
( 2567)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2568)        WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Getting info from file:', filename
( 2569)        CALL FLUSH
( 2570)     END IF
( 2571)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 2572)     IF(nfstat /= nf90_noerr) THEN
( 2573)        status = NCDF_NFERR
( 2574)        RETURN
( 2575)     END IF
( 2576)     
( 2577)     ! Read data from netCDF file
( 2578)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2579)        WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Reading data from file:', filename
( 2580)        CALL FLUSH
( 2581)     END IF
( 2582)     nfstat = nf90_get_var(ncid, varid, data, (/ 1 /), (/ SIZE(data) /))
( 2583)     IF(nfstat /= nf90_noerr) THEN
( 2584)        status = NCDF_NFERR
( 2585)        CALL FLUSH
( 2586)        RETURN
( 2587)     END IF
( 2588)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2589)        WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Closing file:', filename
( 2590)        CALL FLUSH
( 2591)     END IF
( 2592)     nfstat = nf90_close(ncid)
( 2593)     IF(nfstat /= nf90_noerr) THEN
( 2594)        status = NCDF_NFERR
( 2595)        RETURN
( 2596)     END IF
( 2597)  
( 2598)     ! Sync up processors before returning from subroutine
( 2599)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2600)     IF(mpistat /= 0) THEN
( 2601)        status = NCDF_MPERR
( 2602)        RETURN
( 2603)     END IF
( 2604)     
( 2605)   END SUBROUTINE ncdf_read1d
( 2606) 
( 2607)   ! The inverse of ncdf_write2d. Reads a 2-D array from a netCDF file.






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 135

( 2608)   ! Each processor only reads its local subdomain.
( 2609)   ! filename - file to read from
( 2610)   ! varname - variable to read from
( 2611)   ! data - 2-D array to put the data into (must be a REAL)
( 2612)   ! time - time index to read (If the target variable in the netCDF file is
( 2613)   !        actually a 3-D variable with a time axis, ignored otherwise)
( 2614)   ! status - return status of the subroutine
( 2615)   SUBROUTINE ncdf_read2d(filename, varname, data, time, status)
( 2616)     IMPLICIT NONE
( 2617)     ! Subroutine argument declarations
( 2618)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2619)     INTEGER,INTENT(IN) :: time
( 2620)     REAL,DIMENSION(:,:),INTENT(OUT) :: data
( 2621)     INTEGER,INTENT(OUT) :: status
( 2622) 
( 2623)     ! Local declarations
( 2624)     INTEGER :: ncid,    &  ! netCDF file ID
( 2625)                varid,   &  ! ID of netCDF variable to be written to
( 2626)                mpistat, &  ! MPI library call return status
( 2627)                i,       &
( 2628)                j,       &
( 2629)                k,       &
( 2630)                nfstat,  &  ! netCDF library call return status
( 2631)                ndims       ! Number of dimensions in this variable
( 2632)     
( 2633)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 2634)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 2635)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 2636)     REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: databuf
( 2637)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 2638)     INTEGER :: tindex
( 2639) 
( 2640)     ! Initializations
( 2641)     status = NCDF_NOERR
( 2642)     hastime = .FALSE.
( 2643)     j = 0
( 2644)     var_dimlens = 0
( 2645)     if(time > 0) then
( 2646)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 2647)     else
( 2648)        tindex = -time
( 2649)     endif
( 2650) 
( 2651) 
( 2652)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2653)        WRITE(100,*) 'NCDF DEBUG: Reading 2D array from file:', filename
( 2654)        CALL FLUSH
( 2655)     END IF
( 2656) 
( 2657)     ! Open netCDF file and get info
( 2658)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2659)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Opening file:', filename
( 2660)        CALL FLUSH
( 2661)     END IF
( 2662) 
( 2663)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2664)     IF(nfstat /= nf90_noerr) THEN
( 2665)        status = NCDF_NFERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 136

( 2666)        RETURN
( 2667)     END IF
( 2668)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2669)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Getting info from file:', filename
( 2670)        CALL FLUSH
( 2671)     END IF
( 2672)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 2673)     IF(nfstat /= nf90_noerr) THEN
( 2674)        status = NCDF_NFERR
( 2675)        RETURN
( 2676)     END IF
( 2677)     
( 2678)     ! Determine if this variable contains a time axis
( 2679)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 2680)          ndims=ndims, dimids=var_dimids)
( 2681)     IF(nfstat /= nf90_noerr) THEN
( 2682)        status = NCDF_NFERR
( 2683)        RETURN
( 2684)     END IF
( 2685)     DO i = 1, ndims
( 2686) !!DB 2008.07.07
( 2687) !       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2688)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2689)             name=dname, len=k)
( 2690)        dname=TRIM(dname)
( 2691)        var_dimlens(i - j) = k
( 2692)        IF(nfstat /= nf90_noerr) THEN
( 2693)           status = NCDF_NFERR
( 2694)           RETURN
( 2695)        END IF
( 2696)        IF((dname=='time_counter') .OR. (dname=='time')) THEN
( 2697)           hastime = .TRUE.
( 2698)           j = 1
( 2699)        END IF
( 2700)     END DO
( 2701)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2)))
( 2702)     
( 2703)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2704)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Reading record : ', tindex, 'from file:', filename
( 2705)        CALL FLUSH
( 2706)     END IF
( 2707)     
( 2708)     ! Read the appropriate subdomain for each processor
( 2709)     IF(hastime .EQV. .TRUE.) THEN
( 2710)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 2711)             (/ 1, 1, tindex/), &
( 2712)             (/ var_dimlens(1), var_dimlens(2), 1 /))
( 2713)     ELSE
( 2714)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 2715)             (/ 1, 1 /), &
( 2716)             (/ var_dimlens(1), var_dimlens(2) /))
( 2717)     END IF
( 2718) 
( 2719)     DO j = 1, nlcj
( 2720)        DO i = 1, nlci
( 2721)           data(i, j) = databuf(mig(i), mjg(j))
( 2722)        END DO
( 2723)     END DO






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 137

( 2724)     
( 2725) 
( 2726)     ! All done, close up the netCDF file
( 2727)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2728)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Closing file:', filename
( 2729)        CALL FLUSH
( 2730)     END IF
( 2731)     nfstat = nf90_close(ncid)
( 2732)     IF(nfstat /= nf90_noerr) THEN
( 2733)        status = NCDF_NFERR
( 2734)        RETURN
( 2735)     END IF
( 2736)     DEALLOCATE(databuf)
( 2737)     
( 2738)     ! Sync up processors before returning from subroutine
( 2739)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2740)     IF(mpistat /= 0) THEN
( 2741)        status = NCDF_MPERR
( 2742)        RETURN
( 2743)     END IF
( 2744) 
( 2745)   END SUBROUTINE ncdf_read2d
( 2746) 
( 2747)   ! The inverse of ncdf_write3d. Reads a 3-D array from a netCDF file.
( 2748)   ! Each processor only reads its local subdomain.
( 2749)   ! filename - file to read from
( 2750)   ! varname - variable to read from
( 2751)   ! data - 3-D array to put the data into (must be a REAL)
( 2752)   ! time - time index to read (If the target variable in the netCDF file is
( 2753)   !        actually a 4-D variable with a time axis, ignored otherwise)
( 2754)   ! status - return status of the subroutine
( 2755)   SUBROUTINE ncdf_read3d(filename, varname, data, time, status)
( 2756)     IMPLICIT NONE
( 2757)     ! Subroutine argument declarations
( 2758)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2759)     INTEGER,INTENT(IN) :: time
( 2760)     REAL,DIMENSION(:,:,:),INTENT(OUT) :: data
( 2761)     INTEGER,INTENT(OUT) :: status
( 2762) 
( 2763)     ! Local declarations
( 2764)     INTEGER :: ncid,    &  ! netCDF file ID
( 2765)                varid,   &  ! ID of netCDF variable to be written to
( 2766)                mpistat, &  ! MPI library call return status
( 2767)                nfstat,  &  ! netCDF library call return status
( 2768)                i,       &  ! Loop counter
( 2769)                j,       &  !     "
( 2770)                k,       &  !     "
( 2771)                ndims       ! Number of dimensions in this variable
( 2772) 
( 2773)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 2774)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 2775)     INTEGER :: tindex
( 2776)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 2777)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 2778)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: databuf
( 2779) 
( 2780)     ! Initializations
( 2781)     status = NCDF_NOERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 138

( 2782)     hastime = .FALSE.
( 2783)     var_dimlens = 0
( 2784)     j = 0
( 2785)     if(time > 0) then
( 2786)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 2787)     else
( 2788)        tindex = -time
( 2789)     endif
( 2790) 
( 2791)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2792)        WRITE(100,*) 'NCDF DEBUG: Reading 3D array from file:', filename
( 2793)        CALL FLUSH
( 2794)     END IF
( 2795) 
( 2796)     ! Open netCDF file and get info
( 2797)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2798)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Opening file:', filename
( 2799)        CALL FLUSH
( 2800)     END IF
( 2801) 
( 2802)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2803)     IF(nfstat /= nf90_noerr) THEN
( 2804)        status = NCDF_NFERR
( 2805)        RETURN
( 2806)     END IF
( 2807)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2808)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Getting info from file:', filename
( 2809)        CALL FLUSH
( 2810)     END IF
( 2811)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 2812)     IF(nfstat /= nf90_noerr) THEN
( 2813)        status = NCDF_NFERR
( 2814)        RETURN
( 2815)     END IF
( 2816)     
( 2817)     ! Determine if this variable contains a time axis
( 2818)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 2819)          ndims=ndims, dimids=var_dimids)
( 2820)     IF(nfstat /= nf90_noerr) THEN
( 2821)        status = NCDF_NFERR
( 2822)        RETURN
( 2823)     END IF
( 2824)     DO i = 1, ndims
( 2825)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2826)             name=dname, len=k)
( 2827) !       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2828) !            name=dname, len=k)
( 2829)        dname=TRIM(dname)
( 2830)        var_dimlens(i-j) = k
( 2831)        IF(nfstat /= nf90_noerr) THEN
( 2832)           status = NCDF_NFERR
( 2833)           RETURN
( 2834)        END IF
( 2835)        IF((dname=='time_counter') .OR. (dname=='time')) THEN
( 2836)           hastime = .TRUE.
( 2837)           j = 1
( 2838)           EXIT
( 2839)        END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 139

( 2840)     END DO
( 2841)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3)))
( 2842) 
( 2843)     ! Read the appropriate subdomain for each processor
( 2844)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2845)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Reading data from file:', filename
( 2846)        CALL FLUSH
( 2847)     END IF
( 2848)     IF(hastime .EQV. .TRUE.) THEN
( 2849)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 2850)             (/ 1, 1, 1, tindex/), &
( 2851)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), 1 /))
( 2852)     ELSE
( 2853)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 2854)             (/ 1, 1, 1 /), &
( 2855)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3) /))
( 2856)     END IF
( 2857)     IF(nfstat /= nf90_noerr) THEN
( 2858)        status = NCDF_NFERR
( 2859)        RETURN
( 2860)     END IF
( 2861) 
( 2862)     DO k = 1, var_dimlens(3)
( 2863)        DO j = 1, nlcj
( 2864)           DO i = 1, nlci
( 2865)              data(i, j, k) = databuf(mig(i), mjg(j), k)
( 2866)           END DO
( 2867)        END DO
( 2868)     END DO
( 2869) 
( 2870) 
( 2871)     ! All done, close up the netCDF file
( 2872)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2873)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Closing file:', filename
( 2874)        CALL FLUSH
( 2875)     END IF
( 2876)     nfstat = nf90_close(ncid)
( 2877)     IF(nfstat /= nf90_noerr) THEN
( 2878)        status = NCDF_NFERR
( 2879)        RETURN
( 2880)     END IF
( 2881) 
( 2882) 
( 2883)     DEALLOCATE(databuf)
( 2884) 
( 2885)     ! Sync up processors before returning from subroutine
( 2886)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 2887)     IF(mpistat /= 0) THEN
( 2888)        status = NCDF_MPERR
( 2889)        RETURN
( 2890)     END IF
( 2891) 
( 2892)   END SUBROUTINE ncdf_read3d
( 2893) 
( 2894)   ! The inverse of ncdf_write3d. Reads a 3-D array from a netCDF file.
( 2895)   ! Each processor only reads its local subdomain.
( 2896)   ! filename - file to read from
( 2897)   ! varname - variable to read from






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 140

( 2898)   ! data - 3-D array to put the data into (must be a REAL)
( 2899)   ! dsz - Dimension SiZe - size of the 4th dimension in the variable
( 2900)   ! time - time index to read (If the target variable in the netCDF file is
( 2901)   !        actually a 4-D variable with a time axis, ignored otherwise)
( 2902)   ! status - return status of the subroutine
( 2903)   ! NOT EXTENSIVELY TESTED - MAY STILL BE BUGGY
( 2904)   SUBROUTINE ncdf_read4d(filename, varname, data, dsz, time, status)
( 2905)     IMPLICIT NONE
( 2906)     ! Subroutine argument declarations
( 2907)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 2908)     INTEGER,INTENT(IN) :: time, dsz
( 2909)     REAL,DIMENSION(:,:,:,:),INTENT(OUT) :: data
( 2910)     INTEGER,INTENT(OUT) :: status
( 2911)     
( 2912)     ! Local declarations
( 2913)     INTEGER :: ncid,    &  ! netCDF file ID
( 2914)                varid,   &  ! ID of netCDF variable to be written to
( 2915)                mpistat, &  ! MPI library call return status
( 2916)                nfstat,  &  ! netCDF library call return status
( 2917)                i,k,       &  ! Loop counter
( 2918)                ndims       ! Number of dimensions in this variable
( 2919) 
( 2920)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 2921)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 2922)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 2923)     INTEGER :: tindex
( 2924)     
( 2925)     ! Initializations
( 2926)     status = NCDF_NOERR
( 2927)     hastime = .FALSE.
( 2928) 
( 2929)     if(time > 0) then
( 2930)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 2931)     else
( 2932)        tindex = -time
( 2933)     endif
( 2934) !    tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 2935)     
( 2936)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 2937)        WRITE(100,*) 'NCDF DEBUG: Reading 4D array from file:', filename
( 2938)        CALL FLUSH
( 2939)     END IF
( 2940)     
( 2941)     ! Open netCDF file and get info
( 2942)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2943)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Opening file:', filename
( 2944)        CALL FLUSH
( 2945)     END IF
( 2946) 
( 2947)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 2948)     IF(nfstat /= nf90_noerr) THEN
( 2949)        status = NCDF_NFERR
( 2950)        RETURN
( 2951)     END IF
( 2952)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 2953)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Getting info from file:', filename
( 2954)        CALL FLUSH
( 2955)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 141

( 2956)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 2957)     IF(nfstat /= nf90_noerr) THEN
( 2958)        status = NCDF_NFERR
( 2959)        RETURN
( 2960)     END IF
( 2961)     
( 2962)     ! Determine if this variable contains a time axis
( 2963)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 2964)          ndims=ndims, dimids=var_dimids)
( 2965)     IF(nfstat /= nf90_noerr) THEN
( 2966)        status = NCDF_NFERR
( 2967)        RETURN
( 2968)     END IF
( 2969)     DO i = 1, ndims
( 2970)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 2971)             name=dname, len=k)
( 2972) !       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 2973) !            name=dname)
( 2974)        dname=TRIM(dname)
( 2975)        IF(nfstat /= nf90_noerr) THEN
( 2976)           status = NCDF_NFERR
( 2977)           RETURN
( 2978)        END IF
( 2979)        IF(dname=='time_counter') THEN
( 2980)           hastime = .TRUE.
( 2981)           EXIT
( 2982)        END IF
( 2983)     END DO
( 2984)     
( 2985)     ! Got the data, write it to the proper location in the netCDF file
( 2986)     IF(hastime .EQV. .TRUE.) THEN
( 2987)        nfstat = nf90_get_var(ncid, varid, data, &
( 2988)             (/ nimpp, njmpp, 1, 1, tindex/), &
( 2989)             (/ jpi, jpj, SIZE(data,3)-1, SIZE(data,4)-1, 1 /))
( 2990)     ELSE
( 2991)        nfstat = nf90_get_var(ncid, varid, data, &
( 2992)             (/ nimpp, njmpp, 1, 1 /), &
( 2993)             (/ jpj, jpi, SIZE(data,3)-1, SIZE(data,4)-1 /))
( 2994)     END IF
( 2995)     IF(nfstat /= nf90_noerr) THEN
( 2996)        status = NCDF_NFERR
( 2997)        RETURN
( 2998)     END IF
( 2999)     
( 3000)     ! All done, close up the netCDF file
( 3001)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3002)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Closing file:', filename
( 3003)        CALL FLUSH
( 3004)     END IF
( 3005)     nfstat = nf90_close(ncid)
( 3006)     IF(nfstat /= nf90_noerr) THEN
( 3007)        status = NCDF_NFERR
( 3008)        RETURN
( 3009)     END IF
( 3010)     
( 3011)     ! Sync up processors before returning from subroutine
( 3012)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3013)     IF(mpistat /= 0) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 142

( 3014)        status = NCDF_MPERR
( 3015)        RETURN
( 3016)     END IF
( 3017) 
( 3018)   END SUBROUTINE ncdf_read4d
( 3019) 
( 3020)   ! ncdf_create_file creates an empty netCDF dataset with the specified filename
( 3021)   ! filename - name of the file to be created
( 3022)   ! status - return status
( 3023)   SUBROUTINE ncdf_create_file(filename, status)
( 3024)     IMPLICIT NONE
( 3025)     
( 3026)     ! Subroutine argument declarations
( 3027)     CHARACTER(LEN=*),INTENT(IN) :: filename
( 3028)     INTEGER,INTENT(OUT) :: status
( 3029) 
( 3030)     ! Local declarations
( 3031)     INTEGER :: ncid,   &  ! file ID
( 3032)                nfstat, &  ! netCDF call return status
( 3033)                mpistat    ! MPI call return status
( 3034) 
( 3035)     ! Initialization
( 3036)     status = NCDF_NOERR
( 3037) 
( 3038)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3039)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_file - Creating:', filename
( 3040)        CALL FLUSH
( 3041)     END IF
( 3042) 
( 3043)     IF(nproc == 0) THEN
( 3044)        nfstat = nf90_create(filename, NF90_CLOBBER, ncid)
( 3045)        IF(nfstat /= nf90_noerr) THEN
( 3046)           status = NCDF_NFERR
( 3047)           RETURN
( 3048)        END IF
( 3049)        nfstat = nf90_close(ncid)
( 3050)        IF(nfstat /= nf90_noerr) THEN
( 3051)           status = NCDF_NFERR
( 3052)           RETURN
( 3053)        END IF
( 3054)     END IF
( 3055) 
( 3056)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3057)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_file - Done creating:', filename
( 3058)        CALL FLUSH
( 3059)     END IF
( 3060) 
( 3061)     ! Sync up processors before returning from subroutine
( 3062)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3063)     IF(mpistat /= 0) THEN
( 3064)        status = NCDF_MPERR
( 3065)        RETURN
( 3066)     END IF
( 3067) 
( 3068)   END SUBROUTINE ncdf_create_file
( 3069) 
( 3070)   ! ncdf_create_dim creates a new dimension in an existing dataset
( 3071)   ! filename - name of the file to modified






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 143

( 3072)   ! dimname - name of the dimension to create
( 3073)   ! dim_len - length of the new dimension
( 3074)   ! status - return status
( 3075)   SUBROUTINE ncdf_create_dim(filename, dimname, dim_len, status)
( 3076)     IMPLICIT NONE
( 3077) 
( 3078)     ! Subroutine argument declarations
( 3079)     CHARACTER(LEN=*),INTENT(IN) :: filename, dimname
( 3080)     INTEGER,INTENT(IN) :: dim_len
( 3081)     INTEGER,INTENT(OUT) :: status
( 3082) 
( 3083)     ! Local declarations
( 3084)     INTEGER :: ncid,    &  ! file ID
( 3085)                nfstat,  &  ! netCDF call return status
( 3086)                mpistat, &  ! MPI call return status
( 3087)                dimid,   &  ! New dimension ID 
( 3088)                dimlen
( 3089) 
( 3090)     ! Initialization
( 3091)     status = NCDF_NOERR
( 3092)     dimlen = dim_len
( 3093) 
( 3094)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3095)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_dim - Creating dimension in file:', filename
( 3096)        CALL FLUSH
( 3097)     END IF
( 3098) 
( 3099)     IF(nproc == 0) THEN
( 3100)        IF(dimlen == -1) THEN
( 3101)           dimlen = NF90_UNLIMITED
( 3102)        END IF
( 3103)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3104)        IF(nfstat /= nf90_noerr) THEN
( 3105)           status = NCDF_NFERR
( 3106)           RETURN
( 3107)        END IF
( 3108)        nfstat = nf90_redef(ncid)
( 3109)        IF(nfstat /= nf90_noerr) THEN
( 3110)           status = NCDF_NFERR
( 3111)           RETURN
( 3112)        END IF
( 3113)        nfstat = nf90_def_dim(ncid, dimname, dimlen, dimid)
( 3114)        IF(nfstat /= nf90_noerr) THEN
( 3115)           status = NCDF_NFERR
( 3116)           RETURN
( 3117)        END IF
( 3118)        nfstat = nf90_close(ncid)
( 3119)        IF(nfstat /= nf90_noerr) THEN
( 3120)           status = NCDF_NFERR
( 3121)           RETURN
( 3122)        END IF
( 3123)     END IF
( 3124) 
( 3125)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3126)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_dim - Done creating dimension in file:', filename
( 3127)        CALL FLUSH
( 3128)     END IF
( 3129)     






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 144

( 3130)     ! Sync up processors before returning from subroutine
( 3131)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3132)     IF(mpistat /= 0) THEN
( 3133)        status = NCDF_MPERR
( 3134)        RETURN
( 3135)     END IF
( 3136) 
( 3137)   END SUBROUTINE ncdf_create_dim
( 3138) 
( 3139)   ! ncdf_create_var creates a new variable in an existing dataset
( 3140)   ! filename - name of the file to modified
( 3141)   ! varname - name of the variable to create
( 3142)   ! dimname - an array containing the names (in order) of the dimensions
( 3143)   !           the new variable will use
( 3144)   ! datatype - type of the variable (NCDF_FLOAT or NCDF_DOUBLE)
( 3145)   ! status - return status
( 3146)   SUBROUTINE ncdf_create_var(filename, varname, dimnames, datatype, status)
( 3147)     IMPLICIT NONE
( 3148)     
( 3149)     ! Subroutine argument declarations
( 3150)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 3151)     CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: dimnames
( 3152)     INTEGER,INTENT(IN) :: datatype
( 3153)     INTEGER,INTENT(OUT) :: status
( 3154) 
( 3155)     ! Local declarations
( 3156)     INTEGER :: ncid,    &  ! file ID
( 3157)                nfstat,  &  ! netCDF call return status
( 3158)                mpistat, &  ! MPI call return status
( 3159)                varid,   &  ! New variable ID
( 3160)                i
( 3161)     INTEGER,DIMENSION(1:NF90_MAX_VAR_DIMS) :: dimids
( 3162) 
( 3163)     ! Initialization
( 3164)     status = NCDF_NOERR
( 3165)     dimids = -1
( 3166) 
( 3167)     IF((datatype /= NCDF_FLOAT) .AND. (datatype /= NCDF_DOUBLE)) THEN
( 3168)        status = NCDF_ARERR
( 3169)        RETURN
( 3170)     END IF
( 3171) 
( 3172)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3173)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_var - Creating variable in file:', filename
( 3174)        CALL FLUSH
( 3175)     END IF
( 3176) 
( 3177)     IF(nproc == 0) THEN
( 3178)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3179)        IF(nfstat /= nf90_noerr) THEN
( 3180)           status = NCDF_NFERR
( 3181)           RETURN
( 3182)        END IF
( 3183)        ! Get dimension IDs
( 3184)        DO i = 1, SIZE(dimnames)
( 3185)           nfstat = nf90_inq_dimid(ncid, dimnames(i), dimids(i))
( 3186)           IF(nfstat /= nf90_noerr) THEN
( 3187)              status = NCDF_NFERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 145

( 3188)              RETURN
( 3189)           END IF
( 3190)        END DO
( 3191)        ! Create variable
( 3192)        nfstat = nf90_redef(ncid)
( 3193)        IF(nfstat /= nf90_noerr) THEN
( 3194)           status = NCDF_NFERR
( 3195)           RETURN
( 3196)        END IF
( 3197)        nfstat = nf90_def_var(ncid, varname, datatype, &
( 3198)             dimids(1:SIZE(dimnames)), varid)
( 3199)        IF(nfstat /= nf90_noerr) THEN
( 3200)           status = NCDF_NFERR
( 3201)           RETURN
( 3202)        END IF
( 3203)        nfstat = nf90_close(ncid)
( 3204)        IF(nfstat /= nf90_noerr) THEN
( 3205)           status = NCDF_NFERR
( 3206)           RETURN
( 3207)        END IF
( 3208)     END IF
( 3209) 
( 3210)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3211)        WRITE(100,*) 'NCDF DEBUG: ncdf_create_var - Done creating variable in file:', filename
( 3212)        CALL FLUSH
( 3213)     END IF
( 3214) 
( 3215)     ! Sync up processors before returning from subroutine
( 3216)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3217)     IF(mpistat /= 0) THEN
( 3218)        status = NCDF_MPERR
( 3219)        RETURN
( 3220)     END IF
( 3221) 
( 3222)   END SUBROUTINE ncdf_create_var
( 3223) 
( 3224)   ! ncdf_put_att_int creates a new attribute in an existing dataset
( 3225)   ! filename - name of the file to modified
( 3226)   ! varname - name of the variable to attach the attribute to (or "GLOBAL")
( 3227)   ! attname - name of the new attribute
( 3228)   ! attval - value the attribute should have
( 3229)   ! status - return status
( 3230)   SUBROUTINE ncdf_put_att_int(filename, varname, attname, attval, status)
( 3231)     IMPLICIT NONE
( 3232) 
( 3233)     ! Subroutine argument declarations
( 3234)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname
( 3235)     INTEGER,INTENT(IN) :: attval
( 3236)     INTEGER,INTENT(OUT) :: status
( 3237) 
( 3238)     ! Local declarations
( 3239)     INTEGER :: ncid,    &  ! file ID
( 3240)                nfstat,  &  ! netCDF call return status
( 3241)                mpistat, &  ! MPI call return status
( 3242)                varid       ! ID of variable to attach attribute to
( 3243) 
( 3244)     ! Initialization
( 3245)     status = NCDF_NOERR






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 146

( 3246) 
( 3247)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3248)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_int - Creating attribute in file:', filename
( 3249)        CALL FLUSH
( 3250)     END IF
( 3251) 
( 3252)     IF(nproc == 0) THEN
( 3253)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3254)        IF(nfstat /= nf90_noerr) THEN
( 3255)           status = NCDF_NFERR
( 3256)           RETURN
( 3257)        END IF
( 3258)        ! Get variable ID
( 3259)        IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN
( 3260)           varid = NF90_GLOBAL
( 3261)        ELSE
( 3262)           nfstat = nf90_inq_varid(ncid, varname, varid)
( 3263)           IF(nfstat /= nf90_noerr) THEN
( 3264)              status = NCDF_NFERR
( 3265)              RETURN
( 3266)           END IF
( 3267)        END IF
( 3268)        ! Write attribute
( 3269)        nfstat = nf90_redef(ncid)
( 3270)        IF(nfstat /= nf90_noerr) THEN
( 3271)           status = NCDF_NFERR
( 3272)           RETURN
( 3273)        END IF
( 3274)        nfstat = nf90_put_att(ncid, varid, attname, attval)
( 3275)        IF(nfstat /= nf90_noerr) THEN
( 3276)           status = NCDF_NFERR
( 3277)           RETURN
( 3278)        END IF
( 3279)        nfstat = nf90_close(ncid)
( 3280)        IF(nfstat /= nf90_noerr) THEN
( 3281)           status = NCDF_NFERR
( 3282)           RETURN
( 3283)        END IF
( 3284)     END IF
( 3285) 
( 3286)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3287)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_int - Done creating attribute in file:', filename
( 3288)        CALL FLUSH
( 3289)     END IF
( 3290) 
( 3291)     ! Sync up processors before returning from subroutine
( 3292)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3293)     IF(mpistat /= 0) THEN
( 3294)        status = NCDF_MPERR
( 3295)        RETURN
( 3296)     END IF
( 3297) 
( 3298)   END SUBROUTINE ncdf_put_att_int
( 3299) 
( 3300)   ! ncdf_put_att_real creates a new attribute in an existing dataset
( 3301)   ! filename - name of the file to modified
( 3302)   ! varname - name of the variable to attach the attribute to (or "GLOBAL")
( 3303)   ! attname - name of the new attribute






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 147

( 3304)   ! attval - value the attribute should have
( 3305)   ! status - return status
( 3306)   SUBROUTINE ncdf_put_att_real(filename, varname, attname, attval, status)
( 3307)     IMPLICIT NONE
( 3308) 
( 3309)     ! Subroutine argument declarations
( 3310)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname
( 3311)     REAL(wp),INTENT(IN) :: attval
( 3312)     INTEGER,INTENT(OUT) :: status
( 3313) 
( 3314)     ! Local declarations
( 3315)     INTEGER :: ncid,    &  ! file ID
( 3316)                nfstat,  &  ! netCDF call return status
( 3317)                mpistat, &  ! MPI call return status
( 3318)                varid       ! ID of variable to attach attribute to
( 3319) 
( 3320)     ! Initialization
( 3321)     status = NCDF_NOERR
( 3322) 
( 3323)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3324)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_real - Creating attribute in file:', filename
( 3325)        CALL FLUSH
( 3326)     END IF
( 3327) 
( 3328)     IF(nproc == 0) THEN
( 3329)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3330)        IF(nfstat /= nf90_noerr) THEN
( 3331)           status = NCDF_NFERR
( 3332)           RETURN
( 3333)        END IF
( 3334)        ! Get variable ID
( 3335)        IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN
( 3336)           varid = NF90_GLOBAL
( 3337)        ELSE
( 3338)           nfstat = nf90_inq_varid(ncid, varname, varid)
( 3339)           IF(nfstat /= nf90_noerr) THEN
( 3340)              status = NCDF_NFERR
( 3341)              RETURN
( 3342)           END IF
( 3343)        END IF
( 3344)        ! Write attribute
( 3345)        nfstat = nf90_redef(ncid)
( 3346)        IF(nfstat /= nf90_noerr) THEN
( 3347)           status = NCDF_NFERR
( 3348)           RETURN
( 3349)        END IF
( 3350)        nfstat = nf90_put_att(ncid, varid, attname, attval)
( 3351)        IF(nfstat /= nf90_noerr) THEN
( 3352)           status = NCDF_NFERR
( 3353)           RETURN
( 3354)        END IF
( 3355)        nfstat = nf90_close(ncid)
( 3356)        IF(nfstat /= nf90_noerr) THEN
( 3357)           status = NCDF_NFERR
( 3358)           RETURN
( 3359)        END IF
( 3360)     END IF
( 3361) 






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 148

( 3362)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3363)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_real - Done creating attribute in file:', filename
( 3364)        CALL FLUSH
( 3365)     END IF
( 3366) 
( 3367)     ! Sync up processors before returning from subroutine
( 3368)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3369)     IF(mpistat /= 0) THEN
( 3370)        status = NCDF_MPERR
( 3371)        RETURN
( 3372)     END IF
( 3373) 
( 3374)   END SUBROUTINE ncdf_put_att_real
( 3375) 
( 3376)   ! ncdf_put_att_char creates a new attribute in an existing dataset
( 3377)   ! filename - name of the file to modified
( 3378)   ! varname - name of the variable to attach the attribute to (or "GLOBAL")
( 3379)   ! attname - name of the new attribute
( 3380)   ! attval - value the attribute should have
( 3381)   ! status - return status
( 3382)   SUBROUTINE ncdf_put_att_char(filename, varname, attname, attval, status)
( 3383)     IMPLICIT NONE
( 3384) 
( 3385)     ! Subroutine argument declarations
( 3386)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname, attval
( 3387)     INTEGER,INTENT(OUT) :: status
( 3388) 
( 3389)     ! Local declarations
( 3390)     INTEGER :: ncid,    &  ! file ID
( 3391)                nfstat,  &  ! netCDF call return status
( 3392)                mpistat, &  ! MPI call return status
( 3393)                varid       ! ID of variable to attach attribute to
( 3394) 
( 3395)     ! Initialization
( 3396)     status = NCDF_NOERR
( 3397) 
( 3398)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3399)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_char - Creating attribute in file:', filename
( 3400)        CALL FLUSH
( 3401)     END IF
( 3402) 
( 3403)     IF(nproc == 0) THEN
( 3404)        nfstat = nf90_open(filename, nf90_write, ncid)
( 3405)        IF(nfstat /= nf90_noerr) THEN
( 3406)           status = NCDF_NFERR
( 3407)           RETURN
( 3408)        END IF
( 3409)        ! Get variable ID
( 3410)        IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN
( 3411)           varid = NF90_GLOBAL
( 3412)        ELSE
( 3413)           nfstat = nf90_inq_varid(ncid, varname, varid)
( 3414)           IF(nfstat /= nf90_noerr) THEN
( 3415)              status = NCDF_NFERR
( 3416)              RETURN
( 3417)           END IF
( 3418)        END IF
( 3419)        ! Write attribute






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 149

( 3420)        nfstat = nf90_redef(ncid)
( 3421)        IF(nfstat /= nf90_noerr) THEN
( 3422)           status = NCDF_NFERR
( 3423)           RETURN
( 3424)        END IF
( 3425)        nfstat = nf90_put_att(ncid, varid, attname, attval)
( 3426)        IF(nfstat /= nf90_noerr) THEN
( 3427)           status = NCDF_NFERR
( 3428)           RETURN
( 3429)        END IF
( 3430)        nfstat = nf90_close(ncid)
( 3431)        IF(nfstat /= nf90_noerr) THEN
( 3432)           status = NCDF_NFERR
( 3433)           RETURN
( 3434)        END IF
( 3435)     END IF
( 3436) 
( 3437)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3438)        WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_char - Done creating attribute in file:', filename
( 3439)        CALL FLUSH
( 3440)     END IF
( 3441) 
( 3442)     ! Sync up processors before returning from subroutine
( 3443)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3444)     IF(mpistat /= 0) THEN
( 3445)        status = NCDF_MPERR
( 3446)        RETURN
( 3447)     END IF
( 3448) 
( 3449)   END SUBROUTINE ncdf_put_att_char
( 3450) 
( 3451)   ! ncdf_errstr takes an error code returned by a subroutine in this module
( 3452)   ! and returns a string indicating the meaning of the error code
( 3453)   SUBROUTINE ncdf_errstr(errcode, errstr)
( 3454)     INTEGER,INTENT(IN) :: errcode
( 3455)     CHARACTER(LEN=80),INTENT(OUT) :: errstr
( 3456)     
( 3457)     errstr = 'No error'
( 3458)     
( 3459)     SELECT CASE (errcode)
( 3460)     CASE (0)
( 3461)        errstr = 'No error'
( 3462)     CASE(1)
( 3463)        errstr = 'netCDF error'
( 3464)     CASE(2)
( 3465)        errstr = 'MPI error'
( 3466)     CASE(3)
( 3467)        errstr = 'Invalid arguments to subroutine'
( 3468)     CASE(4)
( 3469)        errstr = 'Other/unrecognized error'
( 3470)     CASE DEFAULT
( 3471)        errstr = 'Unknown error code'
( 3472)     END SELECT
( 3473) 
( 3474)   END SUBROUTINE ncdf_errstr
( 3475) 
( 3476) 
( 3477) !!DB: 2007.12.11






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 150

( 3478)   ! ncdf_create_file_aveTSUV builds a special OPA output file with most of the default
( 3479)   ! dimensions, variables and attributes for TSUV fields
( 3480)   ! It is designed for average outputs of useful variables:
( 3481)   ! u, v, T, S, taux, tauy, ssh, ...
( 3482) !!2008.06.13 w added ----> 1 more dimid; 2 more varids
( 3483) !float vovecrtz (time_counter, depthw, y, x)
( 3484) !                char units = "m/s"
( 3485) !                float missing_value = 1.000000e+20
( 3486) !                float valid_min = 1.000000e+20
( 3487) !                float valid_max = -1.000000e+20
( 3488) !                char long_name = "Vertical Velocity"
( 3489) !                char short_name = "vovecrtz"
( 3490) !                char online_operation = "ave(x)"
( 3491) !                char axis = "TZYX"
( 3492) !                float interval_operation = 4.800000e+02
( 3493) !                float interval_write = 8.640000e+04
( 3494) !                char associate = "time_counter depthw nav_lat nav_lon"
( 3495) 
( 3496) 
( 3497)   SUBROUTINE ncdf_create_file_aveTSUV(filename, op_type, status)
( 3498)     IMPLICIT NONE
( 3499)     ! Subroutine argument declarations
( 3500)     CHARACTER(LEN=*),INTENT(IN) :: filename
( 3501)     CHARACTER(LEN=*),INTENT(IN) :: op_type
( 3502)     INTEGER,INTENT(OUT) :: status
( 3503) 
( 3504)     ! Local declarations
( 3505)     INTEGER :: ncid,    &  ! netCDF file ID
( 3506)                varid,   &  ! ID of netCDF variable to be written to
( 3507)                nfstat,  &  ! netCDF library call return status
( 3508)                mpistat     ! MPI library call return status
( 3509) !    INTEGER,DIMENSION(1:4) :: dimids
( 3510)     INTEGER,DIMENSION(1:5) :: dimids
( 3511) !    INTEGER,DIMENSION(1:12) :: varids
( 3512) !AD:    INTEGER,DIMENSION(1:14) :: varids
( 3513)     INTEGER,DIMENSION(1:17) :: varids    
( 3514)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
( 3515)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
( 3516)     CHARACTER(LEN=100) :: sec_since    
( 3517)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
( 3518)     INTEGER :: int_opp, &              ! Operation interval
( 3519)                int_wri                 ! Write interval
( 3520)     CHARACTER(LEN=3),PARAMETER :: &
( 3521)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
( 3522)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
( 3523)     
( 3524)     ! Initializations
( 3525)     status = NCDF_NOERR
( 3526)     CALL ioget_calendar(cal_type)
( 3527)     CALL ioget_timestamp(timestamp)
( 3528)     WRITE (UNIT=sec_since, &
( 3529)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 3530)          &  nyear,nmonth,nday,0, 0, 0
( 3531)     WRITE(t_origin, &
( 3532)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
( 3533)          &   nyear,months(nmonth),nday,0,0,0
( 3534)     int_opp = nwrite * rdt
( 3535)     int_wri = nwrite * rdt






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 151

( 3536) 
( 3537)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3538)        WRITE(100,*) 'NCDF DEBUG: Creating default U output file:', filename
( 3539)        CALL FLUSH
( 3540)     END IF
( 3541) 
( 3542)     ! Only processor 0 does anything
( 3543)     IF(nproc == 0) THEN
( 3544)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3545)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Creating file:', filename
( 3546)           CALL FLUSH
( 3547)        END IF
( 3548)        ! Create the file
( 3549)        nfstat = nf90_create(filename, nf90_clobber, ncid)
( 3550)        IF(nfstat /= nf90_noerr) THEN
( 3551)           status = NCDF_NFERR
( 3552)           RETURN
( 3553)        END IF
( 3554)        
( 3555)        ! Define dimensions
( 3556)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3557)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining dimensions in file:', filename
( 3558)           CALL FLUSH
( 3559)        END IF
( 3560)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
( 3561)        nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
( 3562)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
( 3563)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
( 3564)        nfstat = nf90_def_dim(ncid, 'depthw', jpkdta, dimids(5))
( 3565)        
( 3566)        ! Define variables
( 3567)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3568)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining variables in file:', filename
( 3569)           CALL FLUSH
( 3570)        END IF
( 3571)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
( 3572)             (/ dimids(4), dimids(3) /), &
( 3573)             varids(1))
( 3574)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
( 3575)             (/ dimids(4), dimids(3) /), &
( 3576)             varids(2))
( 3577)        nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
( 3578)             (/ dimids(2) /), &
( 3579)             varids(3))
( 3580)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
( 3581)             (/ dimids(1) /), &
( 3582)             varids(4))
( 3583)        nfstat = nf90_def_var(ncid, 'vozocrtx', nf90_float, &
( 3584)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
( 3585)             varids(5))
( 3586)        nfstat = nf90_def_var(ncid, 'sozotaux', nf90_float, &
( 3587)             (/ dimids(4), dimids(3), dimids(1) /), &
( 3588)             varids(6))
( 3589) 
( 3590)        nfstat = nf90_def_var(ncid, 'vomecrty', nf90_float, &
( 3591)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
( 3592)             varids(7))
( 3593)        nfstat = nf90_def_var(ncid, 'sometauy', nf90_float, &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 152

( 3594)             (/ dimids(4), dimids(3), dimids(1) /), &
( 3595)             varids(8))
( 3596) 
( 3597)        nfstat = nf90_def_var(ncid, 'votemper', nf90_float, &
( 3598)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
( 3599)             varids(9))
( 3600)        nfstat = nf90_def_var(ncid, 'vosaline', nf90_float, &
( 3601)             (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
( 3602)             varids(10))
( 3603)        nfstat = nf90_def_var(ncid, 'sossheig', nf90_float, &
( 3604)             (/ dimids(4), dimids(3), dimids(1) /), &
( 3605)             varids(11))
( 3606)         nfstat = nf90_def_var(ncid, 'soicecov', nf90_float, &
( 3607)             (/ dimids(4), dimids(3), dimids(1) /), &
( 3608)             varids(12))
( 3609) 
( 3610)        nfstat = nf90_def_var(ncid, 'depthw', nf90_float, &
( 3611)             (/ dimids(5) /), &
( 3612)             varids(13))
( 3613)        nfstat = nf90_def_var(ncid, 'vovecrtz', nf90_float, &
( 3614)             (/ dimids(4), dimids(3), dimids(5), dimids(1) /), &
( 3615)             varids(14))
( 3616) !AD: add 
( 3617)        ! ndate (ndastp)
( 3618)        nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
( 3619)             (/ dimids(1) /), &
( 3620)             varids(15))
( 3621) 
( 3622)        ! ndate (model_time)
( 3623)        nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
( 3624)             (/ dimids(1) /), &
( 3625)             varids(16))
( 3626) 
( 3627)        ! kt 
( 3628)        nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
( 3629)             (/ dimids(1) /), &
( 3630)             varids(17))
( 3631) 
( 3632) 
( 3633)        
( 3634)        ! Add attributes
( 3635)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3636)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_aveTSUV - Writing attributes in file:', filename
( 3637)           CALL FLUSH
( 3638)        END IF
( 3639)        ! nav_lon
( 3640)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
( 3641)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
( 3642)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
( 3643)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
( 3644)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
( 3645) 
( 3646)        ! nav_lat
( 3647)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
( 3648)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
( 3649)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
( 3650)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
( 3651)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 153

( 3652) 
( 3653)        ! depthu
( 3654)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
( 3655)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
( 3656)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
( 3657)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
( 3658)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
( 3659)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')
( 3660) 
( 3661)        ! depthw
( 3662)        nfstat = nf90_put_att(ncid, varids(13), 'units', 'm')
( 3663)        nfstat = nf90_put_att(ncid, varids(13), 'positive', 'unknown')
( 3664)        nfstat = nf90_put_att(ncid, varids(13), 'valid_min', 0.00E00)
( 3665)        nfstat = nf90_put_att(ncid, varids(13), 'valid_max', 5.875141E03)
( 3666)        nfstat = nf90_put_att(ncid, varids(13), 'title', 'depthw')
( 3667)        nfstat = nf90_put_att(ncid, varids(13), 'long_name', 'Vertical W levels')
( 3668) 
( 3669)        ! time_counter
( 3670)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
( 3671)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
( 3672)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
( 3673)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
( 3674)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
( 3675) 
( 3676)        ! vozocrtx
( 3677)        nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
( 3678)        nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
( 3679)        nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
( 3680)        nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
( 3681)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Zonal Current')
( 3682)        nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vozocrtx')
( 3683)        nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
( 3684)        nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
( 3685)        nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
( 3686)        nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
( 3687)        nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')
( 3688) 
( 3689)        ! sozotaux
( 3690)        nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
( 3691)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
( 3692)        nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
( 3693)        nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
( 3694)        nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along i-axis')
( 3695)        nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sozotaux')
( 3696)        nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
( 3697)        nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')
( 3698)        nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
( 3699)        nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
( 3700)        nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')
( 3701) 
( 3702)        ! vomecrty
( 3703)        nfstat = nf90_put_att(ncid, varids(7), 'units', 'm/s')
( 3704)        nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
( 3705)        nfstat = nf90_put_att(ncid, varids(7), 'valid_min', 1.000000E20 )
( 3706)        nfstat = nf90_put_att(ncid, varids(7), 'valid_max', -1.000000E20)
( 3707)        nfstat = nf90_put_att(ncid, varids(7), 'long_name', 'Meridional Current')
( 3708)        nfstat = nf90_put_att(ncid, varids(7), 'short_name', 'vomecrty')
( 3709)        nfstat = nf90_put_att(ncid, varids(7), 'online_operation', TRIM(op_type))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 154

( 3710)        nfstat = nf90_put_att(ncid, varids(7), 'axis', 'TZYX')
( 3711)        nfstat = nf90_put_att(ncid, varids(7), 'interval_operation', int_opp)
( 3712)        nfstat = nf90_put_att(ncid, varids(7), 'interval_write', int_wri)
( 3713)        nfstat = nf90_put_att(ncid, varids(7), 'associate', 'time_counter depthv nav_lat nav_lon')
( 3714) 
( 3715)        ! sometauy
( 3716)        nfstat = nf90_put_att(ncid, varids(8), 'units', 'N/m2')
( 3717)        nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
( 3718)        nfstat = nf90_put_att(ncid, varids(8), 'valid_min', 1.000000E20 )
( 3719)        nfstat = nf90_put_att(ncid, varids(8), 'valid_max', -1.000000E20)
( 3720)        nfstat = nf90_put_att(ncid, varids(8), 'long_name', 'Wind Stress along j-axis')
( 3721)        nfstat = nf90_put_att(ncid, varids(8), 'short_name', 'sometauy')
( 3722)        nfstat = nf90_put_att(ncid, varids(8), 'online_operation', TRIM(op_type))
( 3723)        nfstat = nf90_put_att(ncid, varids(8), 'axis', 'TYX')
( 3724)        nfstat = nf90_put_att(ncid, varids(8), 'interval_operation', int_opp)
( 3725)        nfstat = nf90_put_att(ncid, varids(8), 'interval_write', int_wri)
( 3726)        nfstat = nf90_put_att(ncid, varids(8), 'associate', 'time_counter nav_lat nav_lon')
( 3727) 
( 3728)        ! votemper
( 3729)        nfstat = nf90_put_att(ncid, varids(9), 'units', 'C')
( 3730)        nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)
( 3731)        nfstat = nf90_put_att(ncid, varids(9), 'valid_min', 1.000000E20 )
( 3732)        nfstat = nf90_put_att(ncid, varids(9), 'valid_max', -1.000000E20)
( 3733)        nfstat = nf90_put_att(ncid, varids(9), 'long_name', 'Temperature')
( 3734)        nfstat = nf90_put_att(ncid, varids(9), 'short_name', 'votemper')
( 3735)        nfstat = nf90_put_att(ncid, varids(9), 'online_operation', TRIM(op_type))
( 3736)        nfstat = nf90_put_att(ncid, varids(9), 'axis', 'TZYX')
( 3737)        nfstat = nf90_put_att(ncid, varids(9), 'interval_operation', int_opp)
( 3738)        nfstat = nf90_put_att(ncid, varids(9), 'interval_write', int_wri)
( 3739)        nfstat = nf90_put_att(ncid, varids(9), 'associate', 'time_counter deptht nav_lat nav_lon')
( 3740) 
( 3741)        ! vosaline
( 3742)        nfstat = nf90_put_att(ncid, varids(10), 'units', 'PSU')
( 3743)        nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)
( 3744)        nfstat = nf90_put_att(ncid, varids(10), 'valid_min', 1.000000E20 )
( 3745)        nfstat = nf90_put_att(ncid, varids(10), 'long_name', 'Salinity')
( 3746)        nfstat = nf90_put_att(ncid, varids(10), 'short_name', 'vosaline')
( 3747)        nfstat = nf90_put_att(ncid, varids(10), 'online_operation', TRIM(op_type))
( 3748)        nfstat = nf90_put_att(ncid, varids(10), 'axis', 'TZYX')
( 3749)        nfstat = nf90_put_att(ncid, varids(10), 'interval_operation', int_opp)
( 3750)        nfstat = nf90_put_att(ncid, varids(10), 'interval_write', int_wri)
( 3751)        nfstat = nf90_put_att(ncid, varids(10), 'associate', 'time_counter deptht nav_lat nav_lon')
( 3752) 
( 3753)        ! sossheig
( 3754)        nfstat = nf90_put_att(ncid, varids(11), 'units', 'm')
( 3755)        nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)
( 3756)        nfstat = nf90_put_att(ncid, varids(11), 'valid_min', 1.000000E20 )
( 3757)        nfstat = nf90_put_att(ncid, varids(11), 'valid_max', -1.000000E20)
( 3758)        nfstat = nf90_put_att(ncid, varids(11), 'long_name', 'Sea Surface Height')
( 3759)        nfstat = nf90_put_att(ncid, varids(11), 'short_name', 'sossheig')
( 3760)        nfstat = nf90_put_att(ncid, varids(11), 'online_operation', TRIM(op_type))
( 3761)        nfstat = nf90_put_att(ncid, varids(11), 'axis', 'TYX')
( 3762)        nfstat = nf90_put_att(ncid, varids(11), 'interval_operation', int_opp)
( 3763)        nfstat = nf90_put_att(ncid, varids(11), 'interval_write', int_wri)
( 3764)        nfstat = nf90_put_att(ncid, varids(11), 'associate', 'time_counter nav_lat nav_lon')
( 3765) 
( 3766)        ! soicecov
( 3767)        nfstat = nf90_put_att(ncid, varids(12), 'units', '[0,1]')






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 155

( 3768)        nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)
( 3769)        nfstat = nf90_put_att(ncid, varids(12), 'valid_min', 1.000000E20 )
( 3770)        nfstat = nf90_put_att(ncid, varids(12), 'valid_max', -1.000000E20)
( 3771)        nfstat = nf90_put_att(ncid, varids(12), 'long_name', 'Ice Cover')
( 3772)        nfstat = nf90_put_att(ncid, varids(12), 'short_name', 'soicecov')
( 3773)        nfstat = nf90_put_att(ncid, varids(12), 'online_operation', TRIM(op_type))
( 3774)        nfstat = nf90_put_att(ncid, varids(12), 'axis', 'TYX')
( 3775)        nfstat = nf90_put_att(ncid, varids(12), 'interval_operation', int_opp)
( 3776)        nfstat = nf90_put_att(ncid, varids(12), 'interval_write', int_wri)
( 3777)        nfstat = nf90_put_att(ncid, varids(12), 'associate', 'time_counter nav_lat nav_lon')
( 3778) 
( 3779)        ! vovecrtz
( 3780)        nfstat = nf90_put_att(ncid, varids(14), 'units', 'm/s')
( 3781)        nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)
( 3782)        nfstat = nf90_put_att(ncid, varids(14), 'valid_min', 1.000000E20 )
( 3783)        nfstat = nf90_put_att(ncid, varids(14), 'valid_max', -1.000000E20)
( 3784)        nfstat = nf90_put_att(ncid, varids(14), 'long_name', 'Vertical Velocity')
( 3785)        nfstat = nf90_put_att(ncid, varids(14), 'short_name', 'vovecrtz')
( 3786)        nfstat = nf90_put_att(ncid, varids(14), 'online_operation', TRIM(op_type))
( 3787)        nfstat = nf90_put_att(ncid, varids(14), 'axis', 'TZYX')
( 3788)        nfstat = nf90_put_att(ncid, varids(14), 'interval_operation', int_opp)
( 3789)        nfstat = nf90_put_att(ncid, varids(14), 'interval_write', int_wri)
( 3790)        nfstat = nf90_put_att(ncid, varids(14), 'associate', 'time_counter depthw nav_lat nav_lon')
( 3791) 
( 3792)        ! ndate (ndastp)
( 3793)        nfstat = nf90_put_att(ncid, varids(15), 'units', '=nyear*10000+nmonth*100+nday')
( 3794)        nfstat = nf90_put_att(ncid, varids(15), 'long_name', 'time step date in year/month/day aammjj')
( 3795) 
( 3796)        ! ndate (model_time)
( 3797)        nfstat = nf90_put_att(ncid, varids(16), 'long_name', &
( 3798)             'time step date (when output is writen) in year/month/day aammjj (decimal day)')
( 3799)        nfstat = nf90_put_att(ncid, varids(16), 'units', '=nyear*10000+nmonth*100+nday')
( 3800)        nfstat = nf90_put_att(ncid, varids(16), 'formula1', 'nyear  =   model_time / 10000')       
( 3801)        nfstat = nf90_put_att(ncid, varids(16), 'formula2', & 
( 3802)             'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
( 3803)        nfstat = nf90_put_att(ncid, varids(16), 'formula3', & 
( 3804)             'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           
( 3805) 
( 3806) 
( 3807)        ! global
( 3808)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
( 3809)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
( 3810)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
( 3811)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
( 3812)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
( 3813)       
( 3814)        
( 3815)        ! Close file
( 3816)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3817)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_TSUV - Closing file:', filename
( 3818)           CALL FLUSH
( 3819)        END IF
( 3820)        nfstat = nf90_close(ncid)
( 3821)        IF(nfstat /= nf90_noerr) THEN
( 3822)           status = NCDF_NFERR
( 3823)           RETURN
( 3824)        END IF
( 3825)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 156

( 3826) 
( 3827)     ! Sync up processors before returning from subroutine
( 3828)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 3829)     IF(mpistat /= 0) THEN
( 3830)        status = NCDF_MPERR
( 3831)        RETURN
( 3832)     END IF
( 3833)     
( 3834)   END SUBROUTINE ncdf_create_file_aveTSUV
( 3835) 
( 3836) !!DB 2008.05.22
( 3837)   ! ncdf_create_ice_restart builds a standard OPA lim ice restart file with all the default
( 3838)   ! dimensions, variables and attributes. Note that this is a single restart
( 3839)   ! file, and therefore incompatible with the default IOIPSL-based restart
( 3840)   ! routines (though it holds the same data)
( 3841)   SUBROUTINE ncdf_create_ice_restart(filename, status)
( 3842)     IMPLICIT NONE
( 3843)     ! Subroutine argument declarations
( 3844)     CHARACTER(LEN=*),INTENT(IN) :: filename
( 3845)     INTEGER,INTENT(OUT) :: status
( 3846)     
( 3847)     ! Local declarations
( 3848)     INTEGER :: ncid,    &  ! netCDF file ID
( 3849)          varid,   &  ! ID of netCDF variable to be written to
( 3850)          nfstat,  &  ! netCDF library call return status
( 3851)          mpistat     ! MPI library call return status
( 3852)     INTEGER,DIMENSION(1:9) :: dimids
( 3853)     INTEGER,DIMENSION(1:36) :: varids
( 3854)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
( 3855)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
( 3856)     CHARACTER(LEN=100) :: sec_since
( 3857)     CHARACTER(LEN=100) :: tstp_since
( 3858)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
( 3859)     CHARACTER(LEN=3),PARAMETER :: &
( 3860)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
( 3861)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
( 3862)     INTEGER :: ji
( 3863) 
( 3864)     
( 3865)     ! Initializations
( 3866)     status = NCDF_NOERR
( 3867)     
( 3868)     ! Initializations
( 3869)     status = NCDF_NOERR
( 3870)     CALL ioget_calendar(cal_type)
( 3871)     CALL ioget_timestamp(timestamp)
( 3872)     WRITE (UNIT=sec_since, &
( 3873)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 3874)          &  nyear,nmonth,nday,0, 0, 0
( 3875)     WRITE (UNIT=tstp_since, &
( 3876)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
( 3877)          &  nyear,nmonth,nday,0, 0, 0
( 3878)     WRITE(t_origin, &
( 3879)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
( 3880)          &   nyear,months(nmonth),nday,0,0,0
( 3881) 
( 3882)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 3883)        WRITE(100,*) 'NCDF DEBUG: Creating default restart file:', filename






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 157

( 3884)        CALL FLUSH
( 3885)     END IF
( 3886)     
( 3887)     ! Only processor 0 does anything
( 3888)     IF(nproc == 0) THEN
( 3889)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3890)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Creating file:', filename
( 3891)           CALL FLUSH
( 3892)        END IF
( 3893)        ! Create the file
( 3894)        nfstat = nf90_create(filename, nf90_clobber, ncid)
( 3895)        
( 3896)        ! Define dimensions
( 3897)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3898)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining dimensions in file:', filename
( 3899)           CALL FLUSH
( 3900)        END IF
( 3901)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(1))
( 3902)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(2))
( 3903)        nfstat = nf90_def_dim(ncid, 'z', 1, dimids(3))
( 3904)        nfstat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimids(4))
( 3905)        nfstat = nf90_def_dim(ncid, 'x_a', 1, dimids(5))
( 3906)        nfstat = nf90_def_dim(ncid, 'y_a', 1, dimids(6))
( 3907)        nfstat = nf90_def_dim(ncid, 'z_a', 2, dimids(7))
( 3908)        nfstat = nf90_def_dim(ncid, 'z_b', 3, dimids(8))
( 3909)        nfstat = nf90_def_dim(ncid, 'z_c', 35, dimids(9))
( 3910)        
( 3911)        ! Define variables
( 3912)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3913)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining variables in file:', filename
( 3914)           CALL FLUSH
( 3915)        END IF
( 3916)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
( 3917)             (/ dimids(1), dimids(2) /), &
( 3918)             varids(1))
( 3919)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
( 3920)             (/ dimids(1), dimids(2) /), &
( 3921)             varids(2))
( 3922)        nfstat = nf90_def_var(ncid, 'nav_lev', nf90_float, &
( 3923)             (/ dimids(3) /), &
( 3924)             varids(3))
( 3925)        nfstat = nf90_def_var(ncid, 'time', nf90_float, &
( 3926)             (/ dimids(4) /), &
( 3927)             varids(4))
( 3928)        nfstat = nf90_def_var(ncid, 'time_steps', nf90_int, &
( 3929)             (/ dimids(4) /), &
( 3930)             varids(5))
( 3931)        nfstat = nf90_def_var(ncid, 'info', nf90_double, &
( 3932)             (/ dimids(7) /), &
( 3933)             varids(6))
( 3934) 
( 3935)        nfstat = nf90_def_var(ncid, 'hicif'  , nf90_double, &
( 3936)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3937)             varids(7))
( 3938)        nfstat = nf90_def_var(ncid, 'hsnif'   , nf90_double, &
( 3939)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3940)             varids(8))
( 3941)        nfstat = nf90_def_var(ncid, 'frld'    , nf90_double, &






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 158

( 3942)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3943)             varids(9))
( 3944)        nfstat = nf90_def_var(ncid, 'sist'    , nf90_double, &
( 3945)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3946)             varids(10))
( 3947)        nfstat = nf90_def_var(ncid, 'tbif'   ,  nf90_double, &
( 3948)             (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
( 3949)             varids(11))
( 3950)        nfstat = nf90_def_var(ncid, 'u_ice'  ,  nf90_double, &
( 3951)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3952)             varids(12))
( 3953)        nfstat = nf90_def_var(ncid, 'v_ice'  ,  nf90_double, &
( 3954)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3955)             varids(13))
( 3956)        nfstat = nf90_def_var(ncid, 'gtaux'  ,  nf90_double, &
( 3957)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3958)             varids(14))
( 3959)        nfstat = nf90_def_var(ncid, 'gtauy'  ,  nf90_double, &
( 3960)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3961)             varids(15))
( 3962)        nfstat = nf90_def_var(ncid, 'qstoif' ,  nf90_double, &
( 3963)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3964)             varids(16))
( 3965)        nfstat = nf90_def_var(ncid, 'fsbbq'  ,  nf90_double, &
( 3966)             (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
( 3967)             varids(17))
( 3968)        nfstat = nf90_def_var(ncid, 'moment' ,  nf90_double, &
( 3969)             (/ dimids(1), dimids(2), dimids(9), dimids(4) /), &
( 3970)             varids(18))
( 3971) 
( 3972) # 3977
( 3977)        
( 3978)        ! Add attributes
( 3979)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 3980)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Writing attributes in file:', filename
( 3981)           CALL FLUSH
( 3982)        END IF
( 3983)        ! nav_lon
( 3984)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
( 3985)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -1.800000E2)
( 3986)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', 1.800000E2)
( 3987)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
( 3988) 
( 3989)        ! nav_lat
( 3990)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
( 3991)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', -9.000000E1)
( 3992)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 9.000000E1)
( 3993)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
( 3994) 
( 3995)        ! nav_lev
( 3996)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'model_levels')
( 3997)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 0.0)
( 3998)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 0.0)
( 3999)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Model levels')
( 4000) 
( 4001)        ! time
( 4002)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
( 4003)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 159

( 4004)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
( 4005)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'Time axis')
( 4006)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', '0001-JUL-01 00:00:00')
( 4007) 
( 4008)        ! time_steps
( 4009)        nfstat = nf90_put_att(ncid, varids(5), 'units', TRIM(tstp_since))
( 4010)        nfstat = nf90_put_att(ncid, varids(5), 'title', 'Time steps')
( 4011)        nfstat = nf90_put_att(ncid, varids(5), 'tstep_sec', 1.877760E7)
( 4012)        nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Time step axis')
( 4013)        nfstat = nf90_put_att(ncid, varids(5), 'time_origin', TRIM(t_origin))
( 4014) 
( 4015)        ! info
( 4016)        nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
( 4017) !DB as per OLD code
( 4018)        do ji = 7, 18
( 4019)           nfstat = nf90_put_att(ncid, varids(ji), 'missing_value', 1.000000E20)
( 4020)        enddo
( 4021) 
( 4022) # 4025
( 4025)        ! global
( 4026)        nfstat = nf90_put_att(ncid, NF90_GLOBAL , 'Conventions', 'GDT 1.2')
( 4027)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
( 4028)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
( 4029)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number_total', 0)
( 4030)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number', 0)
( 4031)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids', 0)
( 4032)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_global', 0)
( 4033)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_local', 0)
( 4034)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_first', 0)
( 4035)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_last', 0)
( 4036)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', 0)
( 4037)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', 0)
( 4038)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_type','box')
( 4039) 
( 4040)        ! Close file
( 4041)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
( 4042)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Closing file:', filename
( 4043)           CALL FLUSH
( 4044)        END IF
( 4045)        nfstat = nf90_close(ncid)
( 4046)        IF(nfstat /= nf90_noerr) THEN
( 4047)           status = NCDF_NFERR
( 4048)           RETURN
( 4049)        END IF
( 4050)     END IF
( 4051) 
( 4052)     ! Sync up processors before returning from subroutine
( 4053)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4054)     IF(mpistat /= 0) THEN
( 4055)        status = NCDF_MPERR
( 4056)        RETURN
( 4057)     END IF
( 4058)        
( 4059)   END SUBROUTINE ncdf_create_ice_restart
( 4060) 
( 4061) !! CN 10/14/2008 - finds number of dimensions in the specified variable
( 4062)   SUBROUTINE ncdf_get_num_dims(filename, varname, ndims, status)
( 4063)     IMPLICIT NONE






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 160

( 4064)     ! Subroutine argument declarations
( 4065)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4066)     INTEGER,INTENT(OUT) :: ndims, status
( 4067) 
( 4068)     ! Local declarations
( 4069)     INTEGER :: ncid,    &  ! netCDF file ID
( 4070)                id,      &  ! ID of netCDF variable of interest
( 4071)                nfstat,  &  ! netCDF library call return status
( 4072)                mpistat     ! MPI library call return status
( 4073) 
( 4074) 
( 4075)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4076)        WRITE(100,*) 'NCDF DEBUG: Getting number of dimensions in variable from file:', filename
( 4077)        CALL FLUSH
( 4078)     END IF
( 4079) 
( 4080)     ! Initializations
( 4081)     status = NCDF_NOERR
( 4082) 
( 4083)     ! Open netCDF file and get info
( 4084)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4085)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Opening file:', filename
( 4086)        CALL FLUSH
( 4087)     END IF
( 4088)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4089)     IF(nfstat /= nf90_noerr) THEN
( 4090)        status = NCDF_NFERR
( 4091)        RETURN
( 4092)     END IF
( 4093)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4094)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Getting varid from file:', filename
( 4095)        CALL FLUSH
( 4096)     END IF
( 4097)     nfstat = nf90_inq_varid(ncid, varname, id)
( 4098)     IF(nfstat /= nf90_noerr) THEN
( 4099)        status = NCDF_NFERR
( 4100)        RETURN
( 4101)     END IF
( 4102) 
( 4103)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4104)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Getting number of dims in variable from file:', filename
( 4105)        CALL FLUSH
( 4106)     END IF
( 4107)     nfstat = nf90_inquire_variable(ncid, id, ndims=ndims)
( 4108)     IF(nfstat /= nf90_noerr) THEN
( 4109)        status = NCDF_NFERR
( 4110)        RETURN
( 4111)     END IF
( 4112) 
( 4113)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4114)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Closing file:', filename
( 4115)        CALL FLUSH
( 4116)     END IF
( 4117)     nfstat = nf90_close(ncid)
( 4118)     IF(nfstat /= nf90_noerr) THEN
( 4119)        status = NCDF_NFERR
( 4120)        RETURN
( 4121)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 161

( 4122)  
( 4123)     ! Sync up processors before returning from subroutine
( 4124)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4125)     IF(mpistat /= 0) THEN
( 4126)        status = NCDF_MPERR
( 4127)        RETURN
( 4128)     END IF
( 4129) 
( 4130)   END SUBROUTINE ncdf_get_num_dims
( 4131) 
( 4132) !!DB
( 4133) !!Get size of dimension -- common use is to get size of unlimited dimension
( 4134) !!(i.e. # records)
( 4135) !!Note that all processors open the file and get the dimension size
( 4136)   SUBROUTINE ncdf_get_dim_size(filename, dimname, len, status)
( 4137)     IMPLICIT NONE
( 4138)     ! Subroutine argument declarations
( 4139)     CHARACTER(LEN=*),INTENT(IN) :: filename, dimname
( 4140)     INTEGER,INTENT(OUT) :: len, status
( 4141) 
( 4142)     ! Local declarations
( 4143)     INTEGER :: ncid,    &  ! netCDF file ID
( 4144)                id,      &  ! ID of netCDF dimension of interest
( 4145)                nfstat,  &  ! netCDF library call return status
( 4146)                mpistat     ! MPI library call return status
( 4147)     CHARACTER(LEN=nf90_max_name) dname
( 4148) 
( 4149) 
( 4150)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4151)        WRITE(100,*) 'NCDF DEBUG: Getting dim size from file:', filename
( 4152)        CALL FLUSH
( 4153)     END IF
( 4154) 
( 4155)     ! Initializations
( 4156)     status = NCDF_NOERR
( 4157) 
( 4158)     ! Open netCDF file and get info
( 4159)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4160)        WRITE(100,*) 'NCDF DEBUG: ncdf_get dim_size - Opening file:', filename
( 4161)        CALL FLUSH
( 4162)     END IF
( 4163)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4164)     IF(nfstat /= nf90_noerr) THEN
( 4165)        status = NCDF_NFERR
( 4166)        RETURN
( 4167)     END IF
( 4168)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4169)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Getting dimid from file:', filename
( 4170)        CALL FLUSH
( 4171)     END IF
( 4172)     nfstat = NF90_INQ_DIMID (ncid, dimname, id)
( 4173)     IF(nfstat /= nf90_noerr) THEN
( 4174)        status = NCDF_NFERR
( 4175)        RETURN
( 4176)     END IF
( 4177) 
( 4178)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4179)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Getting dimsize from file:', filename






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 162

( 4180)        CALL FLUSH
( 4181)     END IF
( 4182)     nfstat = NF90_INQUIRE_DIMENSION (ncid, id, len=len)
( 4183) 
( 4184)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4185)        WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Closing file:', filename
( 4186)        CALL FLUSH
( 4187)     END IF
( 4188)     nfstat = nf90_close(ncid)
( 4189)     IF(nfstat /= nf90_noerr) THEN
( 4190)        status = NCDF_NFERR
( 4191)        RETURN
( 4192)     END IF
( 4193)  
( 4194)     ! Sync up processors before returning from subroutine
( 4195)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4196)     IF(mpistat /= 0) THEN
( 4197)        status = NCDF_MPERR
( 4198)        RETURN
( 4199)     END IF
( 4200) 
( 4201)   END SUBROUTINE ncdf_get_dim_size
( 4202) 
( 4203)   
( 4204) !>>> A.D:  2008.09.29
( 4205)   
( 4206)     SUBROUTINE ncdf_readdate(filename, varname, cdata, index, status)
( 4207)     IMPLICIT NONE
( 4208)     ! Subroutine argument declarations
( 4209)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4210)     INTEGER,INTENT(IN) :: index
( 4211)     CHARACTER(LEN=19),INTENT(OUT) :: cdata
( 4212)     INTEGER,INTENT(OUT) :: status
( 4213) 
( 4214)     ! Local declarations
( 4215)     INTEGER :: ncid,    &  ! netCDF file ID
( 4216)                varid,   &  ! ID of netCDF variable to be written to
( 4217)                nfstat,  &  ! netCDF library call return status
( 4218)                mpistat     ! MPI library call return status
( 4219)     INTEGER,DIMENSION(1:5) :: count
( 4220) 
( 4221)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4222)        WRITE(100,*) 'NCDF DEBUG: Reading single date from file:', filename
( 4223)        CALL FLUSH
( 4224)     END IF
( 4225) 
( 4226)     ! Initializations
( 4227)     status = NCDF_NOERR
( 4228)     count(1) = 1
( 4229) 
( 4230)     ! Open netCDF file and get info
( 4231)        
( 4232)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4233)        WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Opening file:', filename
( 4234)        CALL FLUSH
( 4235)     END IF
( 4236)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4237)     IF(nfstat /= nf90_noerr) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 163

( 4238)        status = NCDF_NFERR
( 4239)        RETURN
( 4240)     END IF
( 4241) 
( 4242)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4243)        WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Getting info from file:', filename
( 4244)        CALL FLUSH
( 4245)     END IF
( 4246)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 4247)     IF(nfstat /= nf90_noerr) THEN
( 4248)        status = NCDF_NFERR
( 4249)        RETURN
( 4250)     END IF
( 4251)     
( 4252)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4253)        WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Reading data from file:', filename
( 4254)        CALL FLUSH
( 4255)     END IF
( 4256)     !       nfstat = nf90_get_var(ncid, varid, data, (/ index /))
( 4257)     nfstat = nf90_get_var(ncid, varid, cdata, (/ 1,index /), (/ 19,1 /))           
( 4258)     IF(nfstat /= nf90_noerr) THEN
( 4259)        status = NCDF_NFERR
( 4260)        RETURN
( 4261)     END IF
( 4262) 
( 4263)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4264)        WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Closing file:', filename
( 4265)        CALL FLUSH
( 4266)     END IF
( 4267)     nfstat = nf90_close(ncid)
( 4268)     IF(nfstat /= nf90_noerr) THEN
( 4269)        status = NCDF_NFERR
( 4270)        RETURN
( 4271)     END IF
( 4272)     
( 4273)     ! Sync up processors before returning from subroutine
( 4274)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4275)     IF(mpistat /= 0) THEN
( 4276)        status = NCDF_MPERR
( 4277)        RETURN
( 4278)     END IF
( 4279)     
( 4280)   END SUBROUTINE ncdf_readdate
( 4281) 
( 4282)   
( 4283) !<<< A.D:  
( 4284)   
( 4285) SUBROUTINE ncdf_read2d_global(filename, varname, data, time, status)
( 4286)     IMPLICIT NONE
( 4287)     ! Subroutine argument declarations
( 4288)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4289)     INTEGER,INTENT(IN) :: time
( 4290)     REAL,DIMENSION(:,:),INTENT(OUT) :: data
( 4291)     INTEGER,INTENT(OUT) :: status
( 4292) 
( 4293)     ! Local declarations
( 4294)     INTEGER :: ncid,    &  ! netCDF file ID
( 4295)                varid,   &  ! ID of netCDF variable to be written to






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 164

( 4296)                mpistat, &  ! MPI library call return status
( 4297)                nfstat,  &  ! netCDF library call return status
( 4298)                i,       &  ! Loop counter
( 4299)                j,       &  !     "
( 4300)                k,       &  !     "
( 4301)                ndims       ! Number of dimensions in this variable
( 4302) 
( 4303)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 4304)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 4305)     INTEGER :: tindex
( 4306)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 4307)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 4308)     REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: databuf
( 4309) 
( 4310)     ! Initializations
( 4311)     status = NCDF_NOERR
( 4312)     hastime = .FALSE.
( 4313)     var_dimlens = 0
( 4314)     j = 0
( 4315)     if(time > 0) then
( 4316)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 4317)     else
( 4318)        tindex = -time
( 4319)     endif
( 4320) 
( 4321) 
( 4322)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4323)        WRITE(100,*) 'NCDF DEBUG: Reading 2D array from file:', filename
( 4324)        CALL FLUSH
( 4325)     END IF
( 4326) 
( 4327)     ! Open netCDF file and get info
( 4328)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4329)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Opening file:', filename
( 4330)        CALL FLUSH
( 4331)     END IF
( 4332) 
( 4333)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4334)     IF(nfstat /= nf90_noerr) THEN
( 4335)        status = NCDF_NFERR
( 4336)        RETURN
( 4337)     END IF
( 4338)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4339)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Getting info from file:', filename
( 4340)        CALL FLUSH
( 4341)     END IF
( 4342)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 4343)     IF(nfstat /= nf90_noerr) THEN
( 4344)        status = NCDF_NFERR
( 4345)        RETURN
( 4346)     END IF
( 4347)     
( 4348)     ! Determine if this variable contains a time axis
( 4349)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 4350)          ndims=ndims, dimids=var_dimids)
( 4351)     IF(nfstat /= nf90_noerr) THEN
( 4352)        status = NCDF_NFERR
( 4353)        RETURN






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 165

( 4354)     END IF
( 4355)     DO i = 1, ndims
( 4356)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 4357)             name=dname, len=k)
( 4358)        dname=TRIM(dname)
( 4359)        var_dimlens(i-j) = k
( 4360)        IF(nfstat /= nf90_noerr) THEN
( 4361)           status = NCDF_NFERR
( 4362)           RETURN
( 4363)        END IF
( 4364)        IF((dname=='time_counter') .OR. (dname=='time')) THEN
( 4365)           hastime = .TRUE.
( 4366)           j = 1
( 4367)           EXIT
( 4368)        END IF
( 4369)     END DO
( 4370)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2)))
( 4371) 
( 4372)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4373)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Reading data from file:', filename
( 4374)        CALL FLUSH
( 4375)     END IF
( 4376)     IF(hastime .EQV. .TRUE.) THEN
( 4377)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4378)             (/ 1, 1, tindex/), &
( 4379)             (/ var_dimlens(1), var_dimlens(2), 1 /))
( 4380)     ELSE
( 4381)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4382)             (/ 1, 1 /), &
( 4383)             (/ var_dimlens(1), var_dimlens(2) /))
( 4384)     END IF
( 4385)     IF(nfstat /= nf90_noerr) THEN
( 4386)        status = NCDF_NFERR
( 4387)        RETURN
( 4388)     END IF
( 4389) 
( 4390)     data = databuf
( 4391) 
( 4392)     ! All done, close up the netCDF file
( 4393)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4394)        WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Closing file:', filename
( 4395)        CALL FLUSH
( 4396)     END IF
( 4397)     nfstat = nf90_close(ncid)
( 4398)     IF(nfstat /= nf90_noerr) THEN
( 4399)        status = NCDF_NFERR
( 4400)        RETURN
( 4401)     END IF
( 4402) 
( 4403)     DEALLOCATE(databuf)
( 4404) 
( 4405)     ! Sync up processors before returning from subroutine
( 4406)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4407)     IF(mpistat /= 0) THEN
( 4408)        status = NCDF_MPERR
( 4409)        RETURN
( 4410)     END IF
( 4411) 






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 166

( 4412)   END SUBROUTINE ncdf_read2d_global
( 4413) 
( 4414) 
( 4415) !!DB
( 4416)   ! Reads a 3-D array from a netCDF file.
( 4417)   ! Returns global array to wach processor
( 4418)   ! filename - file to read from
( 4419)   ! varname - variable to read from
( 4420)   ! data - 3-D array to put the data into (must be a REAL)
( 4421)   ! time - time index to read (If the target variable in the netCDF file is
( 4422)   !        actually a 4-D variable with a time axis, ignored otherwise)
( 4423)   ! status - return status of the subroutine
( 4424)   SUBROUTINE ncdf_read3d_global(filename, varname, data, time, status)
( 4425)     IMPLICIT NONE
( 4426)     ! Subroutine argument declarations
( 4427)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4428)     INTEGER,INTENT(IN) :: time
( 4429)     REAL,DIMENSION(:,:,:),INTENT(OUT) :: data
( 4430)     INTEGER,INTENT(OUT) :: status
( 4431) 
( 4432)     ! Local declarations
( 4433)     INTEGER :: ncid,    &  ! netCDF file ID
( 4434)                varid,   &  ! ID of netCDF variable to be written to
( 4435)                mpistat, &  ! MPI library call return status
( 4436)                nfstat,  &  ! netCDF library call return status
( 4437)                i,       &  ! Loop counter
( 4438)                j,       &  !     "
( 4439)                k,       &  !     "
( 4440)                ndims       ! Number of dimensions in this variable
( 4441) 
( 4442)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 4443)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 4444)     INTEGER :: tindex
( 4445)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 4446)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 4447)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: databuf
( 4448) 
( 4449)     ! Initializations
( 4450)     status = NCDF_NOERR
( 4451)     hastime = .FALSE.
( 4452)     var_dimlens = 0
( 4453)     j = 0
( 4454)     if(time > 0) then
( 4455)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 4456)     else
( 4457)        tindex = -time
( 4458)     endif
( 4459) 
( 4460)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4461)        WRITE(100,*) 'NCDF DEBUG: Reading 3D array from file:', filename
( 4462)        CALL FLUSH
( 4463)     END IF
( 4464) 
( 4465)     ! Open netCDF file and get info
( 4466)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4467)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Opening file:', filename
( 4468)        CALL FLUSH
( 4469)     END IF






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 167

( 4470) 
( 4471)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4472)     IF(nfstat /= nf90_noerr) THEN
( 4473)        status = NCDF_NFERR
( 4474)        RETURN
( 4475)     END IF
( 4476)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4477)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Getting info from file:', filename
( 4478)        CALL FLUSH
( 4479)     END IF
( 4480)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 4481)     IF(nfstat /= nf90_noerr) THEN
( 4482)        status = NCDF_NFERR
( 4483)        RETURN
( 4484)     END IF
( 4485)     
( 4486)     ! Determine if this variable contains a time axis
( 4487)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 4488)          ndims=ndims, dimids=var_dimids)
( 4489)     IF(nfstat /= nf90_noerr) THEN
( 4490)        status = NCDF_NFERR
( 4491)        RETURN
( 4492)     END IF
( 4493)     DO i = 1, ndims
( 4494)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 4495)             name=dname, len=k)
( 4496)        dname=TRIM(dname)
( 4497)        var_dimlens(i-j) = k
( 4498)        IF(nfstat /= nf90_noerr) THEN
( 4499)           status = NCDF_NFERR
( 4500)           RETURN
( 4501)        END IF
( 4502)        IF((dname=='time_counter') .OR. (dname=='time')) THEN
( 4503)           hastime = .TRUE.
( 4504)           j = 1
( 4505)           EXIT
( 4506)        END IF
( 4507)     END DO
( 4508)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3)))
( 4509) 
( 4510)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4511)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Reading data from file:', filename
( 4512)        CALL FLUSH
( 4513)     END IF
( 4514)     IF(hastime .EQV. .TRUE.) THEN
( 4515)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4516)             (/ 1, 1, 1, tindex/), &
( 4517)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), 1 /))
( 4518)     ELSE
( 4519)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4520)             (/ 1, 1, 1 /), &
( 4521)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3) /))
( 4522)     END IF
( 4523)     IF(nfstat /= nf90_noerr) THEN
( 4524)        status = NCDF_NFERR
( 4525)        RETURN
( 4526)     END IF
( 4527) 






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 168

( 4528)     data = databuf
( 4529) 
( 4530)     ! All done, close up the netCDF file
( 4531)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4532)        WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Closing file:', filename
( 4533)        CALL FLUSH
( 4534)     END IF
( 4535)     nfstat = nf90_close(ncid)
( 4536)     IF(nfstat /= nf90_noerr) THEN
( 4537)        status = NCDF_NFERR
( 4538)        RETURN
( 4539)     END IF
( 4540) 
( 4541)     DEALLOCATE(databuf)
( 4542) 
( 4543)     ! Sync up processors before returning from subroutine
( 4544)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4545)     IF(mpistat /= 0) THEN
( 4546)        status = NCDF_MPERR
( 4547)        RETURN
( 4548)     END IF
( 4549) 
( 4550)   END SUBROUTINE ncdf_read3d_global
( 4551) 
( 4552)   SUBROUTINE ncdf_read4d_global(filename, varname, data, dsz, time, status)
( 4553)     IMPLICIT NONE
( 4554)     ! Subroutine argument declarations
( 4555)     CHARACTER(LEN=*),INTENT(IN) :: filename, varname
( 4556)     INTEGER,INTENT(IN) :: time, dsz
( 4557)     REAL,DIMENSION(:,:,:,:),INTENT(OUT) :: data
( 4558)     INTEGER,INTENT(OUT) :: status
( 4559)     
( 4560)     ! Local declarations
( 4561)     INTEGER :: ncid,    &  ! netCDF file ID
( 4562)                varid,   &  ! ID of netCDF variable to be written to
( 4563)                mpistat, &  ! MPI library call return status
( 4564)                nfstat,  &  ! netCDF library call return status
( 4565)                i,j,k,       &  ! Loop counter
( 4566)                ndims       ! Number of dimensions in this variable
( 4567) 
( 4568)     INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
( 4569)     INTEGER,DIMENSION(1:5) :: var_dimlens
( 4570)     LOGICAL :: hastime   ! Whether this variable has a time axis
( 4571)     CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
( 4572)     INTEGER :: tindex
( 4573)     REAL(wp),ALLOCATABLE,DIMENSION(:,:,:,:) :: databuf
( 4574)     
( 4575)     ! Initializations
( 4576)     status = NCDF_NOERR
( 4577)     hastime = .FALSE.
( 4578)     j=0
( 4579) 
( 4580)     if(time > 0) then
( 4581)        tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
( 4582)     else
( 4583)        tindex = -time
( 4584)     endif
( 4585) !    tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 169

( 4586)     
( 4587)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4588)        WRITE(100,*) 'NCDF DEBUG: Reading 4D array from file:', filename
( 4589)        CALL FLUSH
( 4590)     END IF
( 4591)     
( 4592)     ! Open netCDF file and get info
( 4593)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4594)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Opening file:', filename
( 4595)        CALL FLUSH
( 4596)     END IF
( 4597) 
( 4598)     nfstat = nf90_open(filename, nf90_nowrite, ncid)
( 4599)     IF(nfstat /= nf90_noerr) THEN
( 4600)        status = NCDF_NFERR
( 4601)        RETURN
( 4602)     END IF
( 4603)     IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
( 4604)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Getting info from file:', filename
( 4605)        CALL FLUSH
( 4606)     END IF
( 4607)     nfstat = nf90_inq_varid(ncid, varname, varid)
( 4608)     IF(nfstat /= nf90_noerr) THEN
( 4609)        status = NCDF_NFERR
( 4610)        RETURN
( 4611)     END IF
( 4612)     
( 4613)     ! Determine if this variable contains a time axis
( 4614)     nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
( 4615)          ndims=ndims, dimids=var_dimids)
( 4616)     IF(nfstat /= nf90_noerr) THEN
( 4617)        status = NCDF_NFERR
( 4618)        RETURN
( 4619)     END IF
( 4620)     DO i = 1, ndims
( 4621)        nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
( 4622)             name=dname, len=k)
( 4623) !       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
( 4624) !            name=dname)
( 4625)        dname=TRIM(dname)
( 4626)        var_dimlens(i-j) = k
( 4627)        IF(nfstat /= nf90_noerr) THEN
( 4628)           status = NCDF_NFERR
( 4629)           RETURN
( 4630)        END IF
( 4631)        IF(dname=='time_counter') THEN
( 4632)           hastime = .TRUE.
( 4633)           j = 1
( 4634)           EXIT
( 4635)        END IF
( 4636)     END DO
( 4637)     ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3),1:var_dimlens(4)))
( 4638)     
( 4639)     ! Read data
( 4640)     IF(hastime .EQV. .TRUE.) THEN
( 4641)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4642)             (/ 1, 1, 1, 1, tindex/), &
( 4643)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), var_dimlens(4), 1 /))






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 170

( 4644)     ELSE
( 4645)        nfstat = nf90_get_var(ncid, varid, databuf, &
( 4646)             (/ 1, 1, 1, 1 /), &
( 4647)             (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), var_dimlens(4) /))
( 4648)     END IF
( 4649)     IF(nfstat /= nf90_noerr) THEN
( 4650)        status = NCDF_NFERR
( 4651)        RETURN
( 4652)     END IF
( 4653) 
( 4654)     data = databuf
( 4655)     DEALLOCATE(databuf)
( 4656)     
( 4657)     ! All done, close up the netCDF file
( 4658)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
( 4659)        WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Closing file:', filename
( 4660)        CALL FLUSH
( 4661)     END IF
( 4662)     nfstat = nf90_close(ncid)
( 4663)     IF(nfstat /= nf90_noerr) THEN
( 4664)        status = NCDF_NFERR
( 4665)        RETURN
( 4666)     END IF
( 4667)     
( 4668)     ! Sync up processors before returning from subroutine
( 4669)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 4670)     IF(mpistat /= 0) THEN
( 4671)        status = NCDF_MPERR
( 4672)        RETURN
( 4673)     END IF
( 4674) 
( 4675)   END SUBROUTINE ncdf_read4d_global
( 4676) 
( 4677) 
( 4678) 
( 4679) END MODULE lib_ncdf
