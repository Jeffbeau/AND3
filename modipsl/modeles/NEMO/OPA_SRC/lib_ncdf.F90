!-------------------------------------------------------------------------------
!      OPA NETCDF I/O LIBRARY
!     ========================
!
! This module is intended to provide a workaround for some of the difficulties,
! eccentricities or general insanities of the IOIPSL library that is integrated
! with OPA. All output routines in this module write to a single file, rather
! than a separate file for each processor (thus avoiding the necessity of
! merging files offline). There is no measurable time overhead incurred by
! doing this, so we're not sacrificing model speed to avoid merging files.
! Similarly, all read routines need only a single file. This library also allows
! for easy creation and structuring of new output files. Unlike IOIPSL, it
! performs no averaging on the data - these routines simply write or read
! the requested data, and nothing else. You also don't need to worry about
! opening and closing files, ID numbers, or anything else - that's all handled
! for you. This library should work for any number of processors in any layout.
!
! The 2D and 3D variable writing routines attempt to deal with time axes
! somewhat 'intelligently'. When a 2D array is passed to ncdf_write, the
! corresponding subroutine checks the specified variable for a time dimension.
! If one exists, the variable in the netCDF file is, obviously, 3D, even though
! only a 2D array was passed. The subroutine determines the correct index on
! the time axis to write to, then puts the data there. If no time axis is found,
! the array is simply written to the variable as you'd expect. The same sort
! of thing happens when a 3D array is passed. I hope that made some sense.
!
! NOTE - I've had issues with the FLUSH subroutine causing segmentation faults
! on some systems. If this happens, just comment out any calls to FLUSH and
! recompiler. The problem should go away.
!
! Chris Nickerson
! October & November 2007
! nickersonc@mar.dfo-mpo.gc.ca
!-------------------------------------------------------------------------------

MODULE lib_ncdf
  USE netcdf
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE calendar
  USE daymod

  IMPLICIT NONE

!!DB: PRIVATE statement omitted as the include below becomes hidden to USErs
# include <mpif.h>


  ! Error constant definitions
  INTEGER,PARAMETER :: NCDF_NOERR = 0, &  ! No error, normal return status
                       NCDF_NFERR = 1, &  ! netCDF-related error occurred
                       NCDF_MPERR = 2, &  ! MPI-related error occurred
                       NCDF_ARERR = 3, &  ! Invalid arguments were given
                       NCDF_OTHER = 4     ! Some other error happened
  
  ! Datatype definitions used for creating variables
  INTEGER,PARAMETER :: NCDF_FLOAT = nf90_float, &
                       NCDF_DOUBLE = nf90_double

  ! Turn debugging output on/off (goes to fort.100 - it's VERY verbose)
  LOGICAL,PARAMETER :: DEBUG_OUT = .FALSE. ! True if you want debugging output,
                                           ! false otherwise. Debug output is
                                           ! written to OPA.out in the run dir
                                           ! Note that debug output is extremely
                                           ! verbose!
  
!!DB 2008.05.22 -- added ice restart
!!DB 2008.06.26 -- added get_dim_size
  ! Interface definitions
  PUBLIC :: ncdf_create_file_u, ncdf_create_file_v, ncdf_create_file_t,  ncdf_write, &
       ncdf_errstr, ncdf_create_file, ncdf_create_dim, ncdf_create_var, ncdf_put_att, &
       ncdf_create_restart,  ncdf_create_ice_restart, ncdf_get_dim_size, ncdf_readdate, &
       ncdf_read_global, ncdf_create_file_aveTSUV, ncdf_create_file_ave

!  PUBLIC :: output_special             ! special routines called by step.F90 if ave flag is on
!  PUBLIC :: output_aveTSUV             ! special routines called by step.F90 if ave flag is on

  
  INTERFACE ncdf_write
     MODULE PROCEDURE ncdf_writesv, ncdf_write1d, ncdf_write2d, ncdf_write3d, ncdf_write4d
  END INTERFACE

  INTERFACE ncdf_read
     MODULE PROCEDURE ncdf_readsv, ncdf_read1d, ncdf_read2d, ncdf_read3d, ncdf_read4d
  END INTERFACE

  INTERFACE ncdf_read_global
     MODULE PROCEDURE ncdf_read2d_global, ncdf_read3d_global, ncdf_read4d_global
  END INTERFACE

  INTERFACE ncdf_put_att
     MODULE PROCEDURE ncdf_put_att_int, ncdf_put_att_real, ncdf_put_att_char
  END INTERFACE

CONTAINS

  ! ncdf_create_file_u builds a standard u-grid OPA output file with all the default
  ! dimensions, variables and attributes
  SUBROUTINE ncdf_create_file_u(filename, op_type, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    CHARACTER(LEN=*),INTENT(IN) :: op_type
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:4) :: dimids
    INTEGER,DIMENSION(1:6) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since    
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    INTEGER :: int_opp, &              ! Operation interval
               int_wri                 ! Write interval
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    
    ! Initializations
    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0
    int_opp = nwrite * rdt
    int_wri = nwrite * rdt

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating default U output file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Creating file'
!JC          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(filename, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining dimensions in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining variables in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(2))
       nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
            (/ dimids(2) /), &
            varids(3))
       nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
            (/ dimids(1) /), &
            varids(4))
       nfstat = nf90_def_var(ncid, 'vozocrtx', nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(5))
       nfstat = nf90_def_var(ncid, 'sozotaux', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(6))
       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Writing attributes in file'
!JC          CALL FLUSH
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! depthu
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')

       ! time_counter
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! vozocrtx
       nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
       nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Zonal Current')
       nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vozocrtx')
       nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')

       ! sozotaux
       nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
       nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along i-axis')
       nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sozotaux')
       nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
      
       
       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat) ,trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_file_u

  ! ncdf_create_file_v builds a standard v-grid OPA output file with all the default
  ! dimensions, variables and attributes
  SUBROUTINE ncdf_create_file_v(filename, op_type, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    CHARACTER(LEN=*),INTENT(IN) :: op_type
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:4) :: dimids
    INTEGER,DIMENSION(1:6) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    INTEGER :: int_opp, &              ! Operation interval
               int_wri                 ! Write interval
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    
    ! Initializations
    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0
    int_opp = nwrite * rdt
    int_wri = nwrite * rdt

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating default V output file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Creating file'
!JC          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(filename, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Defining dimensions in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'depthv', jpkdta, dimids(2))
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Defining variables in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(2))
       nfstat = nf90_def_var(ncid, 'depthv', nf90_float, &
            (/ dimids(2) /), &
            varids(3))
       nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
            (/ dimids(1) /), &
            varids(4))
       nfstat = nf90_def_var(ncid, 'vomecrty', nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(5))
       nfstat = nf90_def_var(ncid, 'sometauy', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(6))
       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Writing attributes in file'
!JC          CALL FLUSH
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! depthu
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')

       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthv')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical V levels')

       ! time_counter
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! vomecrty
       nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
       nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Meridional Current')
       nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vomecrty')
       nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthv nav_lat nav_lon')

       ! sometauy
       nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
       nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along j-axis')
       nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sometauy')
       nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
      
       
       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_v - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_file_v

  ! ncdf_create_file_t builds a standard t-grid OPA output file with all the default
  ! dimensions, variables and attributes
  SUBROUTINE ncdf_create_file_t(filename, op_type, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    CHARACTER(LEN=*),INTENT(IN) :: op_type
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:4) :: dimids
    INTEGER,DIMENSION(1:22) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    INTEGER :: int_opp, &              ! Operation interval
               int_wri                 ! Write interval
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    
    ! Initializations
    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0
    int_opp = nwrite * rdt
    int_wri = nwrite * rdt

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating default T output file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Creating file'
!JC          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(filename, nf90_clobber, ncid)
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Defining dimensions in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'jpkdta', jpkdta, dimids(2))
       nfstat = nf90_def_dim(ncid, 'jpjdta', jpjdta, dimids(3))
       nfstat = nf90_def_dim(ncid, 'jpidta', jpidta, dimids(4))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Defining variables in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(2))
       nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &
            (/ dimids(2) /), &
            varids(3))
       nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
            (/ dimids(1) /), &
            varids(4))
       nfstat = nf90_def_var(ncid, 'votemper', nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(5))
       nfstat = nf90_def_var(ncid, 'vosaline', nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(6))
       nfstat = nf90_def_var(ncid, 'sosstsst', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(7))
       nfstat = nf90_def_var(ncid, 'sosaline', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(8))
       nfstat = nf90_def_var(ncid, 'sossheig', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(9))
       nfstat = nf90_def_var(ncid, 'sowaflup', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(10))
       nfstat = nf90_def_var(ncid, 'sorunoff', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(11))
       nfstat = nf90_def_var(ncid, 'sowaflcd', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(12))
       nfstat = nf90_def_var(ncid, 'sosalflx', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(13))
       nfstat = nf90_def_var(ncid, 'sohefldo', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(14))
       nfstat = nf90_def_var(ncid, 'soshfldo', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(15))
       nfstat = nf90_def_var(ncid, 'somxl010', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(16))
       nfstat = nf90_def_var(ncid, 'somixhgt', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(17))
       nfstat = nf90_def_var(ncid, 'soicecov', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(18))
       nfstat = nf90_def_var(ncid, 'sohefldp', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(19))
       nfstat = nf90_def_var(ncid, 'sowafldp', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(20))
       nfstat = nf90_def_var(ncid, 'sosafldp', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(21))
       nfstat = nf90_def_var(ncid, 'sobowlin', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(22))
       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_T - Writing attributes in file'
!JC          CALL FLUSH
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! deptht
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical T levels')

       ! time_counter
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! votemper
       nfstat = nf90_put_att(ncid, varids(5), 'units', 'C')
       nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Temperature')
       nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'votemper')
       nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter deptht nav_lat nav_lon')

       ! vosaline
       nfstat = nf90_put_att(ncid, varids(6), 'units', 'PSU')
       nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Salinity')
       nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'vosaline')
       nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter deptht nav_lat nav_lon')

       ! sosstsst
       nfstat = nf90_put_att(ncid, varids(7), 'units', 'C')
       nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(7), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(7), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(7), 'long_name', 'Sea Surface temperature')
       nfstat = nf90_put_att(ncid, varids(7), 'short_name', 'sosstsst')
       nfstat = nf90_put_att(ncid, varids(7), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(7), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(7), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(7), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(7), 'associate', 'time_counter nav_lat nav_lon')

       ! sosaline
       nfstat = nf90_put_att(ncid, varids(8), 'units', 'PSU')
       nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(8), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(8), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(8), 'long_name', 'Sea Surface Salinity')
       nfstat = nf90_put_att(ncid, varids(8), 'short_name', 'sosaline')
       nfstat = nf90_put_att(ncid, varids(8), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(8), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(8), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(8), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(8), 'associate', 'time_counter nav_lat nav_lon')

       ! sossheig
       nfstat = nf90_put_att(ncid, varids(9), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(9), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(9), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(9), 'long_name', 'Sea Surface Height')
       nfstat = nf90_put_att(ncid, varids(9), 'short_name', 'sossheig')
       nfstat = nf90_put_att(ncid, varids(9), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(9), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(9), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(9), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(9), 'associate', 'time_counter nav_lat nav_lon')
       
       ! sowaflup
       nfstat = nf90_put_att(ncid, varids(10), 'units', 'Kg/m2/s')
       nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(10), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(10), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(10), 'long_name', 'Net Upward Water Flux')
       nfstat = nf90_put_att(ncid, varids(10), 'short_name', 'sowaflup')
       nfstat = nf90_put_att(ncid, varids(10), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(10), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(10), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(10), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(10), 'associate', 'time_counter nav_lat nav_lon')
       
       ! sorunoff
       nfstat = nf90_put_att(ncid, varids(11), 'units', 'Kg/m2/s')
       nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(11), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(11), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(11), 'long_name', 'Runoffs')
       nfstat = nf90_put_att(ncid, varids(11), 'short_name', 'sorunoff')
       nfstat = nf90_put_att(ncid, varids(11), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(11), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(11), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(11), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(11), 'associate', 'time_counter nav_lat nav_lon')
       
       ! sowaflcd
       nfstat = nf90_put_att(ncid, varids(12), 'units', 'kg/m2/s')
       nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(12), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(12), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(12), 'long_name', 'concentration/dilution water flux')
       nfstat = nf90_put_att(ncid, varids(12), 'short_name', 'sowaflcd')
       nfstat = nf90_put_att(ncid, varids(12), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(12), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(12), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(12), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(12), 'associate', 'time_counter nav_lat nav_lon')
       
       ! sosalflx
       nfstat = nf90_put_att(ncid, varids(13), 'units', 'Kg/m2/s')
       nfstat = nf90_put_att(ncid, varids(13), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(13), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(13), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(13), 'long_name', 'Surface Salt Flux')
       nfstat = nf90_put_att(ncid, varids(13), 'short_name', 'sosalflx')
       nfstat = nf90_put_att(ncid, varids(13), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(13), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(13), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(13), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(13), 'associate', 'time_counter nav_lat nav_lon')

       ! sohefldo
       nfstat = nf90_put_att(ncid, varids(14), 'units', 'W/m2')
       nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(14), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(14), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(14), 'long_name', 'Net Downward Heat Flux')
       nfstat = nf90_put_att(ncid, varids(14), 'short_name', 'sohefldo')
       nfstat = nf90_put_att(ncid, varids(14), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(14), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(14), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(14), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(14), 'associate', 'time_counter nav_lat nav_lon')

       ! soshfldo
       nfstat = nf90_put_att(ncid, varids(15), 'units', 'W/m2')
       nfstat = nf90_put_att(ncid, varids(15), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(15), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(15), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(15), 'long_name', 'Shortwave Radiation')
       nfstat = nf90_put_att(ncid, varids(15), 'short_name', 'soshfldo')
       nfstat = nf90_put_att(ncid, varids(15), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(15), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(15), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(15), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(15), 'associate', 'time_counter nav_lat nav_lon')

       ! soml010
       nfstat = nf90_put_att(ncid, varids(16), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(16), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(16), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(16), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(16), 'long_name', 'Mixed Layer Depth 0.01')
       nfstat = nf90_put_att(ncid, varids(16), 'short_name', 'soml010')
       nfstat = nf90_put_att(ncid, varids(16), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(16), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(16), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(16), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(16), 'associate', 'time_counter nav_lat nav_lon')

       ! somixhgt
       nfstat = nf90_put_att(ncid, varids(17), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(17), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(17), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(17), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(17), 'long_name', 'Turbocline Depth')
       nfstat = nf90_put_att(ncid, varids(17), 'short_name', 'somixhgt')
       nfstat = nf90_put_att(ncid, varids(17), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(17), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(17), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(17), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(17), 'associate', 'time_counter nav_lat nav_lon')

       ! soicecov
       nfstat = nf90_put_att(ncid, varids(18), 'units', '[0,1]')
       nfstat = nf90_put_att(ncid, varids(18), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(18), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(18), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(18), 'long_name', 'Ice Cover')
       nfstat = nf90_put_att(ncid, varids(18), 'short_name', 'soicecov')
       nfstat = nf90_put_att(ncid, varids(18), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(18), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(18), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(18), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(18), 'associate', 'time_counter nav_lat nav_lon')

       ! sohefldp
       nfstat = nf90_put_att(ncid, varids(19), 'units', 'W/m2')
       nfstat = nf90_put_att(ncid, varids(19), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(19), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(19), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(19), 'long_name', 'Surface Heat Flux: Damping')
       nfstat = nf90_put_att(ncid, varids(19), 'short_name', 'sohefldp')
       nfstat = nf90_put_att(ncid, varids(19), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(19), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(19), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(19), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(19), 'associate', 'time_counter nav_lat nav_lon')

       ! sowafldp
       nfstat = nf90_put_att(ncid, varids(20), 'units', 'Kg/m2/s')
       nfstat = nf90_put_att(ncid, varids(20), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(20), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(20), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(20), 'long_name', 'Surface Water Flux: Damping')
       nfstat = nf90_put_att(ncid, varids(20), 'short_name', 'sowafldp')
       nfstat = nf90_put_att(ncid, varids(20), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(20), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(20), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(20), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(20), 'associate', 'time_counter nav_lat nav_lon')

       ! sosafldp
       nfstat = nf90_put_att(ncid, varids(21), 'units', 'Kg/m2/s')
       nfstat = nf90_put_att(ncid, varids(21), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(21), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(21), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(21), 'long_name', 'Surface salt flux: damping')
       nfstat = nf90_put_att(ncid, varids(21), 'short_name', 'sosafldp')
       nfstat = nf90_put_att(ncid, varids(21), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(21), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(21), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(21), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(21), 'associate', 'time_counter nav_lat nav_lon')

       ! sobowlin
       nfstat = nf90_put_att(ncid, varids(22), 'units', 'W-point')
       nfstat = nf90_put_att(ncid, varids(22), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(22), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(22), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(22), 'long_name', 'Bowl Index')
       nfstat = nf90_put_att(ncid, varids(22), 'short_name', 'sobowlin')
       nfstat = nf90_put_att(ncid, varids(22), 'online_operation', 'l_max(only(x))')
       nfstat = nf90_put_att(ncid, varids(22), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(22), 'associate', 'time_counter nav_lat nav_lon')

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
      
       
       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_t - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_file_t

  ! ncdf_create_file_ave builds an OPA output file with a single variable
  ! whose name is specified by VARNAME. Other than that, it is similar to a
  ! U- or V-grid file
 SUBROUTINE ncdf_create_file_ave(filename, varname, op_type, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    CHARACTER(LEN=*),INTENT(IN) :: op_type, varname
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:4) :: dimids
    INTEGER,DIMENSION(1:5) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since    
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    INTEGER :: int_opp, &              ! Operation interval
               int_wri                 ! Write interval
    INTEGER :: varnum
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    
    ! Initializations
    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0
    int_opp = nwrite * rdt
    int_wri = nwrite * rdt

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating output file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Creating file'
!JC          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(filename, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename, trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Defining dimensions in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Defining variables in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(2))
       nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
            (/ dimids(2) /), &
            varids(3))
       nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
            (/ dimids(1) /), &
            varids(4))
       nfstat = nf90_def_var(ncid, varname, nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(5))
       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Writing attributes in file'
!JC          CALL FLUSH
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! depthu
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')

       ! time_counter
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! custom variable
       nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
       nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(5), 'short_name', varname)
       nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')

       varnum = 5
!AD/DB: add new time-related variables
       varnum = varnum + 1
       ! ndate (ndastp)
       nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
            (/ dimids(1) /), &
            varids(varnum))
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'time step date in year/month/day aammjj')

       varnum = varnum + 1
       ! ndate (model_time)
       nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
            (/ dimids(1) /), &
            varids(varnum))
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', &
            'time step date (when output is writen) in year/month/day aammjj (decimal day)')
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula1', 'nyear  =   model_time / 10000')       
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula2', & 
            'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula3', & 
            'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           

       varnum = varnum + 1
       ! kt 
       nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
            (/ dimids(1) /), &
            varids(varnum))

       
       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_ave - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_file_ave
  
  ! ncdf_create_restart builds a standard OPA restart file with all the default
  ! dimensions, variables and attributes. Note that this is a single restart
  ! file, and therefore incompatible with the default IOIPSL-based restart
  ! routines (though it holds the same data)
  ! NOTE: Some keys we haven't been using are untested, so I don't know if those
  ! fields will be created correctly.
  SUBROUTINE ncdf_create_restart(filename, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER,INTENT(OUT) :: status
    
    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
         varid,   &  ! ID of netCDF variable to be written to
         nfstat,  &  ! netCDF library call return status
         mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:8) :: dimids
    INTEGER,DIMENSION(1:36) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since
    CHARACTER(LEN=100) :: tstp_since
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    
    ! Initializations
    status = NCDF_NOERR
    
    ! Initializations
    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE (UNIT=tstp_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating default restart file:', trim(filename)
!JC       CALL FLUSH
    END IF
    
    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Creating file'
!JC          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(filename, nf90_clobber, ncid)
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining dimensions in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(1))
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(2))
       nfstat = nf90_def_dim(ncid, 'z', jpkdta, dimids(3))
       nfstat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimids(4))
       nfstat = nf90_def_dim(ncid, 'x_a', 1, dimids(5))
       nfstat = nf90_def_dim(ncid, 'y_a', 1, dimids(6))
       nfstat = nf90_def_dim(ncid, 'z_a', 10, dimids(7))
       nfstat = nf90_def_dim(ncid, 'z_b', 1, dimids(8))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining variables in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(1), dimids(2) /), &
            varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(1), dimids(2) /), &
            varids(2))
       nfstat = nf90_def_var(ncid, 'nav_lev', nf90_float, &
            (/ dimids(3) /), &
            varids(3))
       nfstat = nf90_def_var(ncid, 'time', nf90_float, &
            (/ dimids(4) /), &
            varids(4))
       nfstat = nf90_def_var(ncid, 'time_steps', nf90_float, &
            (/ dimids(4) /), &
            varids(5))
       nfstat = nf90_def_var(ncid, 'info', nf90_float, &
            (/ dimids(7) /), &
            varids(6))
       nfstat = nf90_def_var(ncid, 'ub', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(7))
       nfstat = nf90_def_var(ncid, 'vb', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(8))
       nfstat = nf90_def_var(ncid, 'tb', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(9))
       nfstat = nf90_def_var(ncid, 'sb', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(10))
       nfstat = nf90_def_var(ncid, 'rotb', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(11))
       nfstat = nf90_def_var(ncid, 'hdivb', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(12))
       nfstat = nf90_def_var(ncid, 'un', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(13))
       nfstat = nf90_def_var(ncid, 'vn', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(14))
       nfstat = nf90_def_var(ncid, 'tn', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(15))
       nfstat = nf90_def_var(ncid, 'sn', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(16))
       nfstat = nf90_def_var(ncid, 'rotn', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(17))
       nfstat = nf90_def_var(ncid, 'hdivn', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(18))
       nfstat = nf90_def_var(ncid, 'gcx', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(19))
       nfstat = nf90_def_var(ncid, 'gcxb', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(20))
# if defined key_dynspg_rl
       nfstat = nf90_def_var(ncid, 'bsfb', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(21))
       nfstat = nf90_def_var(ncid, 'bsfn', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(22))
       nfstat = nf90_def_var(ncid, 'bsfd', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(27))
# else 
       nfstat = nf90_def_var(ncid, 'sshb', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(21))
       nfstat = nf90_def_var(ncid, 'sshn', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(22))
#  if defined key_dynspg_ts
       nfstat = nf90_def_var(ncid, 'sshb_b', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(23))
       nfstat = nf90_def_var(ncid, 'sshn_b', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(24))
       nfstat = nf90_def_var(ncid, 'un_b', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(25))
       nfstat = nf90_def_var(ncid, 'vn_b', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(26))
#endif
#endif

       ! Fields that are only defined if specific keys are set
# if defined key_zdftke   ||   defined key_esopa
       nfstat = nf90_def_var(ncid, 'en', nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(28))
# endif
# if defined key_ice_lim
       nfstat = nf90_def_var(ncid, 'nfice', nf90_float, &
            (/dimids(8) /), &
            varids(29))
       nfstat = nf90_def_var(ncid, 'sst_io', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(30))
       nfstat = nf90_def_var(ncid, 'sss_io', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(31))
       nfstat = nf90_def_var(ncid, 'u_io', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(32))
       nfstat = nf90_def_var(ncid, 'v_io', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(33))
# if defined key_coupled
       nfstat = nf90_def_var(ncid, 'alb_ice', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(34))
# endif
# endif
# if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
       nfstat = nf90_def_var(ncid, 'nfbulk', nf90_float, &
            (/ dimids(8) /), &
            varids(35))

       nfstat = nf90_def_var(ncid, 'gsst', nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(36))
# endif
       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Writing attributes in file'
!JC          CALL FLUSH
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -1.800000E2)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', 1.800000E2)
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', -9.000000E1)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 9.000000E1)
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')

       ! nav_lev
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'model_levels')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E0)
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E3)
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Model levels')

       ! time
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'Time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', '0001-JUL-01 00:00:00')

       ! time_steps
       nfstat = nf90_put_att(ncid, varids(5), 'units', TRIM(tstp_since))
       nfstat = nf90_put_att(ncid, varids(5), 'title', 'Time steps')
       nfstat = nf90_put_att(ncid, varids(5), 'tstep_sec', 1.877760E7)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Time step axis')
       nfstat = nf90_put_att(ncid, varids(5), 'time_origin', TRIM(t_origin))

       ! info
       nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)

       ! ub
       nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)

       ! vb
       nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)

       ! tb
       nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)

       ! sb
       nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)

       ! rotb
       nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)

       ! hdivb
       nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)

       ! un
       nfstat = nf90_put_att(ncid, varids(13), 'missing_value', 1.000000E20)

       ! vn
       nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)

       ! tn
       nfstat = nf90_put_att(ncid, varids(15), 'missing_value', 1.000000E20)

       ! sn
       nfstat = nf90_put_att(ncid, varids(16), 'missing_value', 1.000000E20)

       ! rotn
       nfstat = nf90_put_att(ncid, varids(17), 'missing_value', 1.000000E20)

       ! hdivn
       nfstat = nf90_put_att(ncid, varids(18), 'missing_value', 1.000000E20)

       ! gcx
       nfstat = nf90_put_att(ncid, varids(19), 'missing_value', 1.000000E20)

       ! gcxb
       nfstat = nf90_put_att(ncid, varids(20), 'missing_value', 1.000000E20)
# if defined key_dynspg_rl
       ! bsfb
       nfstat = nf90_put_att(ncid, varids(21), 'missing_value', 1.000000E20)

       ! bsfn
       nfstat = nf90_put_att(ncid, varids(22), 'missing_value', 1.000000E20)

       !bsfd
       nfstat = nf90_put_att(ncid, varids(27), 'missing_value', 1.000000E20)
# else
       ! sshb
       nfstat = nf90_put_att(ncid, varids(21), 'missing_value', 1.000000E20)

       ! sshn
       nfstat = nf90_put_att(ncid, varids(22), 'missing_value', 1.000000E20)

#  if defined key_dynspg_ts
       ! sshb_b
       nfstat = nf90_put_att(ncid, varids(23), 'missing_value', 1.000000E20)

       ! sshn_b
       nfstat = nf90_put_att(ncid, varids(24), 'missing_value', 1.000000E20)

       ! un_b
       nfstat = nf90_put_att(ncid, varids(25), 'missing_value', 1.000000E20)

       ! vn_b
       nfstat = nf90_put_att(ncid, varids(26), 'missing_value', 1.000000E20)
# endif
# endif
# if defined key_zdftke   ||   defined key_esopa
       ! en
       nfstat = nf90_put_att(ncid, varids(28), 'missing_value', 1.000000E20)
# endif
# if defined key_ice_lim
       ! nfice
       nfstat = nf90_put_att(ncid, varids(29), 'missing_value', 1.000000E20)

       ! sst_io
       nfstat = nf90_put_att(ncid, varids(30), 'missing_value', 1.000000E20)

       ! sss_io
       nfstat = nf90_put_att(ncid, varids(31), 'missing_value', 1.000000E20)

       ! u_io
       nfstat = nf90_put_att(ncid, varids(32), 'missing_value', 1.000000E20)

       ! v_io
       nfstat = nf90_put_att(ncid, varids(33), 'missing_value', 1.000000E20)

# if defined key_coupled
       ! alb_ice
       nfstat = nf90_put_att(ncid, varids(34), 'missing_value', 1.000000E20)

# endif
# endif
# if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
       ! nfbulk
       nfstat = nf90_put_att(ncid, varids(35), 'missing_value', 1.000000E20)
 
       ! gsst
       nfstat = nf90_put_att(ncid, varids(36), 'missing_value', 1.000000E20)
 
# endif

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL , 'Conventions', 'GDT 1.2')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number_total', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_global', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_local', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_first', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_last', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_type','box')

       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
       
  END SUBROUTINE ncdf_create_restart

  ! ncdf_writesv (single value) writes a single, scalar value to a specified
  ! index in a netCDF variable. The netCDF variable is assumed to be
  ! 1-dimensional.
  ! filename - file to write to
  ! varname - variable to write to
  ! data - scalar value to write
  ! index - where in the variable to put the data
  ! status - return status of the subroutine
  SUBROUTINE ncdf_writesv(filename, varname, data, index, status)
     IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    REAL(wp),INTENT(IN) :: data
    INTEGER,INTENT(IN) :: index   ! Where in the variable to write the value
    INTEGER,INTENT(OUT) :: status
    real*4 buf

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Writing single value to file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Initializations
    status = NCDF_NOERR

    ! Open netCDF file and get info
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Opening file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename, trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Getting info from file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_inq_varid(ncid, varname, varid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       
       ! Write data to netCDF file
       ! This subroutine assumes all processors have the same data in this
       ! variable, no attempt to sync data or merge data in the file is made
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Writing data to file'
!JC          CALL FLUSH
       END IF
       buf = real(data,4)
       nfstat = nf90_put_var(ncid, varid, buf, &
           (/ index /))
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_writesv - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_writesv

  ! ncdf_write1d writes a 1-dimensional array to a 1-D variable in a netCDF
  ! file. It assumes that all processors have the same data in the array, 
  ! so no attempt is made to sync data between processors or merge data
  ! together in the file, unlike ncdf_write2d and ncdf_write3d
  ! filename - file to write to
  ! varname - variable to write to
  ! data - 1-D array to write
  ! status - return status of the subroutine
  SUBROUTINE ncdf_write1d(filename, varname, data, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    REAL(wp),DIMENSION(:),INTENT(IN) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    REAL*4,ALLOCATABLE,DIMENSION(:) :: buf   ! Send/receive buffer for data array

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Writing 1D array to file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Initializations
    status = NCDF_NOERR

    ! Open netCDF file and get info
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Opening file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Getting info from file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_inq_varid(ncid, varname, varid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       
       ! Write data to netCDF file
       ! This subroutine assumes all processors have the same data in this
       ! array, no attempt to sync data or merge data in the file is made
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Writing data to file'
!JC          CALL FLUSH
       END IF
       ALLOCATE(buf(size(data,1)))            ! Allocate send/receive buffer
       buf = real(data,4)               ! Copy section of subdomain to buffer
       nfstat = nf90_put_var(ncid, varid, buf )
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
!JC          CALL FLUSH
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write1d - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
      deallocate(buf)
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_write1d

  ! ncdf_write2d writes a 2-dimensional array to a 2-D variable in a netCDF
  ! file. It assumes the data array being passed is a subdomain - that is,
  ! each processor's array contains different values, and the data needs to
  ! be written to the location in the netCDF file which corresponds to the
  ! processor writing the subdomain. This subroutine will correctly merge
  ! all subroutines so that the variable in the netCDF file contains
  ! the entire global domain of values
  ! filename - file to write to
  ! varname - variable to write to
  ! data - 2-D array to write
  ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
  !                                      will force an immediate write, with
  !                                      data being written at a time index in
  !                                      the file equal to (tstep * -1). Use this
  !                                      to output on timesteps which are not
  !                                      multiples of nwrite
  ! status - return status of the subroutine
  ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
  ! a tstep which is <0
  SUBROUTINE ncdf_write2d(filename, varname, data, tstep, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: tstep
    REAL(wp),DIMENSION(:,:),INTENT(IN) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               tindex,  &  ! Index to write to on the time axis
               mpistat, &  ! MPI library call return status
               nfstat,  &  ! netCDF library call return status
               is,      &  ! i-axis start index for writing to file
               js,      &  ! j-axis   "     "    "     "     "  "
               i,k,       &  ! Loop counter
               ndims,   &  ! Number of dimensions in this variabl
               iz,      &  ! I-size of array to write to file
               jz          ! J-size of array to write to file
    
    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
    REAL*4,ALLOCATABLE,DIMENSION(:,:) :: buf   ! Send/receive buffer for data array
    LOGICAL :: hastime   ! Whether this variable has a time axis

    ! Exit subroutine if no write should take place on current timestep
    IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >= 0)) THEN
       status = NCDF_NOERR
       RETURN
    END IF

    ! Initializations
    is = nimpp                          ! Get global i-index for this subdomain
    js = njmpp                          ! Get global j-index  "   "      "
    iz = nlei                           ! Get i-size of subdomain section to be written
    jz = nlej                           ! Get j-size of subdomain section to be written
    ALLOCATE(buf(1:iz,1:jz))            ! Allocate send/receive buffer
    buf = real(data(1:iz,1:jz),4)               ! Copy section of subdomain to buffer
    IF(tstep >= 0) THEN
       tindex = tstep/nwrite - nit000/nwrite
    ELSE
       tindex = tstep * (-1)
    END IF
    status = NCDF_NOERR
    hastime = .FALSE.

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Writing 2D array to file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Open netCDF file and get info
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Opening file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname)
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Getting info from file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_inq_varid(ncid, varname, varid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname)
          RETURN
       END IF
       
       ! Determine if this variable contains a time axis
       nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
            ndims=ndims, dimids=var_dimids)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
          RETURN
       END IF
       DO i = 1, ndims
          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
               name=dname, len=k)
!          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
!               name=dname)
          dname=TRIM(dname)
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
             RETURN
          END IF
          IF(dname=='time_counter') THEN
             hastime = .TRUE.
          END IF
       END DO
    END IF

    ! Get data from each processor in turn, processor 0 writes to file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Writing data to file'
!JC       CALL FLUSH
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    DO i = 0, (jpnij - 1)
       ! Since processor 0 doesn't need to send data to itself, only
       ! processors with IDs > 0 make send calls
       IF(i > 0) THEN
          ! If it's the local processor's turn to send, go ahead and
          ! send data and write indices to processor 0
          IF(nproc == i) THEN
             CALL MPI_SEND(is, 1, MPI_INTEGER, &
                  0, i, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(js, 1, MPI_INTEGER, &
                  0, i+1, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(iz, 1, MPI_INTEGER, &
                  0, i+2, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(jz, 1, MPI_INTEGER, &
                  0, i+3, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(buf, SIZE(buf), MPI_REAL4, &
                  0, i+4, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
          END IF
       END IF

       ! Processor 0 receives data and writes it to file
       IF(nproc == 0) THEN
          ! Since processor 0 doesn't need to receive data from itself, 
          ! we only make receive calls if it's another processor's turn
          ! to send data
          IF(i > 0) THEN
             CALL MPI_RECV(is, 1, MPI_INTEGER, &
                  i, i, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(js, 1, MPI_INTEGER, &
                  i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(iz, 1, MPI_INTEGER, &
                  i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(jz, 1, MPI_INTEGER, &
                  i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             DEALLOCATE(buf)
             ALLOCATE(buf(1:iz,1:jz))
             CALL MPI_RECV(buf, SIZE(buf), MPI_REAL4, &
                  i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
          END IF
          ! Got the data, write it to the proper location in the netCDF file
          IF(hastime .EQV. .TRUE.) THEN
             nfstat = nf90_put_var(ncid, varid, buf, &
                  (/ is, js, tindex/), &
                  (/ iz, jz, 1 /))
          ELSE
             nfstat = nf90_put_var(ncid, varid, buf, &
                  (/ is, js /), &
                  (/ iz, jz /))
          END IF
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
             RETURN
          END IF
       END IF

       ! Sync all processors at the end of each loop iteration to prevent
       ! concurrency issues/race conditions/various MPI-related ugliness
       CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
       IF(mpistat /= 0) THEN
          status = NCDF_MPERR
          RETURN
       END IF
    END DO

    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_write2d - Closing file'
!JC       CALL FLUSH
    END IF
    IF(nproc == 0) THEN
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname)
          RETURN
       END IF
    END IF
    DEALLOCATE(buf)

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_write2d

  ! ncdf_write3d writes a 3-dimensional array to a 3-D variable in a netCDF
  ! file. It assumes the data array being passed is a subdomain - that is,
  ! each processor's array contains different values, and the data needs to
  ! be written to the location in the netCDF file which corresponds to the
  ! processor writing the subdomain. This subroutine will correctly merge
  ! all subroutines so that the variable in the netCDF file contains
  ! the entire global domain of values
  ! filename - file to write to
  ! varname - variable to write to
  ! data - 3-D array to write
  ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
  !                                      will force an immediate write, with
  !                                      data being written at a time index in
  !                                      the file equal to (tstep * -1). Use this
  !                                      to output on timesteps which are not
  !                                      multiples of nwrite
  ! status - return status of the subroutine
  ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
  ! a tstep which is <0
  SUBROUTINE ncdf_write3d(filename, varname, data, tstep, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: tstep
    REAL(wp),DIMENSION(:,:,:),INTENT(IN) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               tindex,  &  ! Index to write to on the time axis
               mpistat, &  ! MPI library call return status
               nfstat,  &  ! netCDF library call return status
               !rcvstat, &  ! MPI receive status
               is,      &  ! i-axis start index for writing to file
               js,      &  ! j-axis   "     "    "     "     "  "
               ks,      &  ! k-axis   "     "    "     "     "  "
               i,k,     &  ! Loop counter
               ndims,   &  ! Number of dimensions in this variable
               iz,      &  ! I-size of array to write to file
               jz          ! J-size of array to write to file

    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
    REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: buf   ! Send/receive buffer for data array
    LOGICAL :: hastime   ! Whether this variable has a time axis
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 

    ! Exit subroutine if no write should take place on current timestep
    IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >=0)) THEN
       status = NCDF_NOERR
       RETURN
    END IF

    ! Initializations
    is = nimpp                          ! Get global i-index for this subdomain
    js = njmpp                          ! Get global j-index  "   "      "
    ks = 1                              ! Depth writes always start at k=1
    iz = nlei                           ! Get i-size of subdomain section to be written
    jz = nlej                           ! Get j-size of subdomain section to be written
    ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,DIM=3)))      ! Allocate send/recieve buffer
    buf = real(data(1:iz,1:jz,1:SIZE(data,DIM=3)),4)         ! Copy section of subdomain to buffer
    IF(tstep >= 0) THEN
       tindex = tstep/nwrite - nit000/nwrite
    ELSE
       tindex = tstep * (-1)
    END IF
    status = NCDF_NOERR
    hastime = .FALSE.

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Writing 3D array to file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Open netCDF file and get info
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Opening file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname)
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Getting info from file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_inq_varid(ncid, varname, varid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname)
          RETURN
       END IF

       ! Determine if this variable contains a time axis
       nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
            ndims=ndims, dimids=var_dimids)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
          RETURN
       END IF
       DO i = 1, ndims
          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
               name=dname, len=k)
!          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
!               name=dname)
          dname=TRIM(dname)
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
             RETURN
          END IF
          IF(dname=='time_counter') THEN
             hastime = .TRUE.
             EXIT
          END IF
       END DO
    END IF

    ! Get data from each processor in turn, processor 0 writes to file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Writing data to file'
!JC       CALL FLUSH
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    DO i=0, (jpnij - 1)
       ! Since processor 0 doesn't need to send data to itself, only
       ! processors with IDs > 0 make send calls
       IF(i > 0) THEN
          ! If it's the local processor's turn to send, go ahead and
          ! send data and write indices to processor 0
          IF(nproc == i) THEN
             CALL MPI_SEND(is, 1, MPI_INTEGER, &
                  0, i, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(js, 1, MPI_INTEGER, &
                  0, i+1, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(ks, 1, MPI_INTEGER, &
                  0, i+2, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(iz, 1, MPI_INTEGER, &
                  0, i+3, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(jz, 1, MPI_INTEGER, &
                  0, i+4, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(buf, SIZE(buf), MPI_REAL4, &
                  0, i+5, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
          END IF
       END IF

       ! Processor 0 receives data and writes it to file
       IF(nproc == 0) THEN
          ! Since processor 0 doesn't need to receive data from itself, 
          ! we only make receive calls if it's another processor's turn
          ! to send data
          IF(i > 0) THEN
             CALL MPI_RECV(is, 1, MPI_INTEGER, &
                  i, i, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(js, 1, MPI_INTEGER, &
                  i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(ks, 1, MPI_INTEGER, &
                  i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(iz, 1, MPI_INTEGER, &
                  i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(jz, 1, MPI_INTEGER, &
                  i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             DEALLOCATE(buf)
             ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,DIM=3)))
             CALL MPI_RECV(buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
                  i, i+5, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
          END IF
          
          ! Got the data, write it to the proper location in the netCDF file
          IF(hastime .EQV. .TRUE.) THEN
             nfstat = nf90_put_var(ncid, varid, buf,&
                  (/ is, js, ks, tindex/), &
                  (/ iz, jz, SIZE(data,DIM=3), 1 /))
          ELSE
             nfstat = nf90_put_var(ncid, varid, buf,&
                  (/ is, js, ks /), &
                  (/ iz, jz, SIZE(data,DIM=3) /))
          END IF
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
             RETURN
          END IF
       END IF
       
       ! Sync all processors at the end of each loop iteration to prevent
       ! concurrency issues/race conditions/various MPI-related ugliness
       CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
       IF(mpistat /= 0) THEN
          status = NCDF_MPERR
          RETURN
       END IF
    END DO

    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_write3d - Closing file'
!JC       CALL FLUSH
    END IF
    IF(nproc == 0) THEN
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
          RETURN
       END IF
    END IF
    DEALLOCATE(buf)

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_write3d
  
  ! ncdf_write4d writes a 4-dimensional array to a 4-D variable in a netCDF
  ! file. It assumes the data array being passed is a subdomain - that is,
  ! each processor's array contains different values, and the data needs to
  ! be written to the location in the netCDF file which corresponds to the
  ! processor writing the subdomain. This subroutine will correctly merge
  ! all subroutines so that the variable in the netCDF file contains
  ! the entire global domain of values
  ! filename - file to write to
  ! varname - variable to write to
  ! data - 3-D array to write
  ! dsz - Dimension SiZe - size of the 4th dimension in the variable
  ! tstep - the current timestep - NOTE: giving a timestep which is less than 0
  !                                      will force an immediate write, with
  !                                      data being written at a time index in
  !                                      the file equal to (tstep * -1). Use this
  !                                      to output on timesteps which are not
  !                                      multiples of nwrite
  ! status - return status of the subroutine
  ! NOTE: Writes will only happen if MOD(tstep, nwrite) == 0 unless you provide
  ! a tstep which is < 0
  ! NOT EXTENSIVELY TESTED - MAY STILL BE BUGGY
  SUBROUTINE ncdf_write4d(filename, varname, data, dsz, tstep, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: tstep, dsz
    REAL(wp),DIMENSION(:,:,:,:),INTENT(IN) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               tindex,  &  ! Index to write to on the time axis
               mpistat, &  ! MPI library call return status
               nfstat,  &  ! netCDF library call return status
               is,      &  ! i-axis start index for writing to file
               js,      &  ! j-axis   "     "    "     "     "  "
               ks,      &  ! k-axis   "     "    "     "     "  "
               i,k,jn,       &  ! Loop counter  !NL##
               ndims,   &  ! Number of dimensions in this variable
               iz,      &  ! I-size of array to write to file
               jz          ! J-size of array to write to file

    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    INTEGER,DIMENSION(1:MPI_STATUS_SIZE) :: rcvstat
    REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: buf   ! Send/receive buffer for data array
    REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: buf3D   ! Send/receive buffer for data array  !NL##
    LOGICAL :: hastime   ! Whether this variable has a time axis
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 

    ! Exit subroutine if no write should take place on current timestep
    IF((MOD(tstep, nwrite) /= 0) .AND. (tstep >= 0)) THEN
       status = NCDF_NOERR
       RETURN
    END IF

    ! Initializations
    is = nimpp                          ! Get global i-index for this subdomain
    js = njmpp                          ! Get global j-index  "   "      "
    ks = 1                              ! Depth writes always start at k=1
    is = nimpp                          ! Get global i-index for this subdomain
    js = njmpp                          ! Get global j-index  "   "      "
    iz = nlei                           ! Get i-size of subdomain section to be written
    jz = nlej                           ! Get j-size of subdomain section to be written
    ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,3),1:SIZE(data,4)))      ! Allocate send/recieve buffer
    ALLOCATE(buf3D(1:iz,1:jz,1:SIZE(data,3)))      ! Allocate send/recieve buffer
    buf = real(data(1:iz,1:jz,1:SIZE(data,3), 1:SIZE(data,4)),4)         ! Copy section of subdomain to buffer
    IF(tstep >= 0) THEN
       tindex = tstep/nwrite - nit000/nwrite
    ELSE
       tindex = tstep * (-1)
    END IF
    status = NCDF_NOERR
    hastime = .FALSE.

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Writing 4D array to file:', trim(filename)
    END IF

    ! Open netCDF file and get info
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Opening file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Getting info from file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_inq_varid(ncid, varname, varid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
          RETURN
       END IF

       ! Determine if this variable contains a time axis
       nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
            ndims=ndims, dimids=var_dimids)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
          RETURN
       END IF
       DO i = 1, ndims
          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
               name=dname, len=k)
!!DB: next line not correct 
!          nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
!               name=dname)
          dname=TRIM(dname)
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
             RETURN
          END IF
          IF(dname=='time_counter') THEN
             hastime = .TRUE.
             EXIT
          END IF
       END DO
    END IF

    ! Get data from each processor in turn, processor 0 writes to file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Writing data to file'
!JC       CALL FLUSH
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
   
    DO i=0, (jpnij - 1)
       ! Since processor 0 doesn't need to send data to itself, only
       ! processors with IDs > 0 make send calls
       IF(i > 0) THEN
          ! If it's the local processor's turn to send, go ahead and
          ! send data and write indices to processor 0
          IF(nproc == i) THEN
             CALL MPI_SEND(is, 1, MPI_INTEGER, &
                  0, i, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(js, 1, MPI_INTEGER, &
                  0, i+1, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(ks, 1, MPI_INTEGER, &
                  0, i+2, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(iz, 1, MPI_INTEGER, &
                  0, i+3, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_SEND(jz, 1, MPI_INTEGER, &
                  0, i+4, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             ! HERE !!
             do jn=1,size(buf,4)
               buf3D=buf(:,:,:,jn)
               CALL MPI_SEND(buf3D, SIZE(buf3D), MPI_REAL4, &
                  0, i+5, MPI_COMM_WORLD, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
               buf(:,:,:,jn)=buf3D
             enddo 
          END IF
       END IF

       ! Processor 0 receives data and writes it to file
       IF(nproc == 0) THEN
          ! Since processor 0 doesn't need to receive data from itself, 
          ! we only make receive calls if it's another processor's turn
          ! to send data
          IF(i > 0) THEN
             CALL MPI_RECV(is, 1, MPI_INTEGER, &
                  i, i, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(js, 1, MPI_INTEGER, &
                  i, i+1, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(ks, 1, MPI_INTEGER, &
                  i, i+2, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(iz, 1, MPI_INTEGER, &
                  i, i+3, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             CALL MPI_RECV(jz, 1, MPI_INTEGER, &
                  i, i+4, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
             DEALLOCATE(buf,buf3D)
             ALLOCATE(buf3D(1:iz,1:jz,1:SIZE(data,3)))
             ALLOCATE(buf(1:iz,1:jz,1:SIZE(data,3),1:SIZE(data,4)))
              ! HERE !!
             do jn =1, size(buf,4)
               buf3D=buf(:,:,:,jn) 
             CALL MPI_RECV(buf3D, SIZE(buf3D), MPI_REAL4, &
                  i, i+5, MPI_COMM_WORLD, rcvstat, mpistat)
             IF(mpistat /= 0) THEN
                status = NCDF_MPERR
                RETURN
             END IF
               buf(:,:,:,jn)=buf3D
             enddo
          END IF
          ! Got the data, write it to the proper location in the netCDF file
          IF(hastime .EQV. .TRUE.) THEN
             nfstat = nf90_put_var(ncid, varid, buf,&
                  (/ is, js, ks, 1, tindex/), &
                  (/ iz, jz, SIZE(data,3), SIZE(data,4), 1 /))
          ELSE
             nfstat = nf90_put_var(ncid, varid, buf,&
                  (/ is, js, ks, 1 /), &
                  (/ iz, jz, SIZE(data,3), SIZE(data,4) /))
          END IF
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
             RETURN
          END IF
       END IF

       ! Sync all processors at the end of each loop iteration to prevent
       ! concurrency issues/race conditions/various MPI-related ugliness
       CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
       IF(mpistat /= 0) THEN
          status = NCDF_MPERR
          RETURN
       END IF
    END DO

    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_write4d - Closing file'
!JC       CALL FLUSH
    END IF
    IF(nproc == 0) THEN
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename),trim(varname), trim(filename)
          RETURN
       END IF
    END IF
    DEALLOCATE(buf,buf3D)

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_write4d

  ! The inverse of ncdf_writesv. Reads a scalar value from a 1-D array in a
  ! netCDF file. All processors will get the same value.
  ! filename - file to read from
  ! varname - variable to read from
  ! data - scalar to put the data into (must be a REAL)
  ! index -the index to read from
  ! status - return status of the subroutine
  SUBROUTINE ncdf_readsv(filename, varname, data, index, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: index
    REAL(wp),INTENT(OUT) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:5) :: count
    REAL*4 buf

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading single value from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Initializations
    status = NCDF_NOERR
    count(1) = 1

    ! Open netCDF file and get info
!!DB: think that there should be no if(nproc==0) then ...
!    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Opening file'
!JC          CALL FLUSH
       END IF
!DB
!       nfstat = nf90_open(filename, nf90_write, ncid)
       nfstat = nf90_open(filename, nf90_nowrite, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename, trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Getting info from file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_inq_varid(ncid, varname, varid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       
       ! Read data from netCDF file
       ! This subroutine assumes all processors will get the same data from this
       ! variable
       IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Reading data from file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_get_var(ncid, varid, buf, (/ index /))
       data=real(buf,8)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_readsv - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
!    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_readsv

  ! The inverse of ncdf_write1d. Reads a 1-D array from a netCDF file
  ! All processors will get the same values.
  ! filename - file to read from
  ! varname - variable to read from
  ! data - 1-D array to put the data into (must be a REAL)
  ! status - return status of the subroutine
  SUBROUTINE ncdf_read1d(filename, varname, data, status, tide)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    REAL(wp),DIMENSION(:),INTENT(OUT) :: data
    INTEGER,INTENT(OUT) :: status
    INTEGER,INTENT(IN),OPTIONAL :: tide

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    REAL*4, ALLOCATABLE,DIMENSION(:) :: buf 


    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading 1D array from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Initializations
    status = NCDF_NOERR

    ! Open netCDF file and get info
    
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Opening file'
    END IF
    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
          WRITE(100,*) 'JC: DEBUG: ncdf_read1d - problem opening file.'
       END IF
    nfstat = nf90_open(filename, nf90_nowrite, ncid)
       RETURN
    END IF
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Getting info from file'
       WRITE(100,*) 'JC: DEBUG: ncdf_read1d - reading file.'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, varid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Read data from netCDF file
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Reading data from file'
!JC       CALL FLUSH
    END IF
    allocate(buf(size(data,1)))
    nfstat = nf90_get_var(ncid, varid, buf, (/ 1 /), (/ SIZE(data) /))
    data=real(buf,8)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
!JC       CALL FLUSH
       RETURN
    END IF
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read1d - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
 
    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_read1d

  ! The inverse of ncdf_write2d. Reads a 2-D array from a netCDF file.
  ! Each processor only reads its local subdomain.
  ! filename - file to read from
  ! varname - variable to read from
  ! data - 2-D array to put the data into (must be a REAL)
  ! time - time index to read (If the target variable in the netCDF file is
  !        actually a 3-D variable with a time axis, ignored otherwise)
  ! status - return status of the subroutine
  SUBROUTINE ncdf_read2d(filename, varname, data, time, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: time
    REAL(wp),DIMENSION(:,:),INTENT(OUT) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               mpistat, &  ! MPI library call return status
               i,       &
               j,       &
               k,       &
               nfstat,  &  ! netCDF library call return status
               ndims       ! Number of dimensions in this variable
    
    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
    LOGICAL :: hastime   ! Whether this variable has a time axis
    REAL*4,ALLOCATABLE,DIMENSION(:,:) :: databuf
    INTEGER,DIMENSION(1:5) :: var_dimlens
    INTEGER :: tindex

    ! Initializations
    status = NCDF_NOERR
    hastime = .FALSE.
    j = 0
    var_dimlens = 0
    if(time > 0) then
       tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
    else
       tindex = -time
    endif


    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading 2D array from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Open netCDF file and get info
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Opening file'
!JC       CALL FLUSH
    END IF

    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Getting info from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, varid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Determine if this variable contains a time axis
    nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
         ndims=ndims, dimids=var_dimids)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    DO i = 1, ndims
!!DB 2008.07.07
!       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
            name=dname, len=k)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       dname=TRIM(dname)
       var_dimlens(i - j) = k
       IF((dname=='time_counter') .OR. (dname=='time')) THEN
          hastime = .TRUE.
          j = 1
       END IF
    END DO
    ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2)))
    
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Reading record : ', tindex, 'from file'
!JC       CALL FLUSH
    END IF
    
    ! Read the appropriate subdomain for each processor
    IF(hastime .EQV. .TRUE.) THEN
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, tindex/), &
            (/ var_dimlens(1), var_dimlens(2), 1 /))
    ELSE
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1 /), &
            (/ var_dimlens(1), var_dimlens(2) /))
    END IF
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF

    DO j = 1, nlcj
       DO i = 1, nlci
          data(i, j) = real(databuf(mig(i), mjg(j)),8)
       END DO
    END DO
    

    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read2d - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    DEALLOCATE(databuf)
    
    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_read2d

  ! The inverse of ncdf_write3d. Reads a 3-D array from a netCDF file.
  ! Each processor only reads its local subdomain.
  ! filename - file to read from
  ! varname - variable to read from
  ! data - 3-D array to put the data into (must be a REAL)
  ! time - time index to read (If the target variable in the netCDF file is
  !        actually a 4-D variable with a time axis, ignored otherwise)
  ! status - return status of the subroutine
  SUBROUTINE ncdf_read3d(filename, varname, data, time, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: time
    REAL(wp),DIMENSION(:,:,:),INTENT(OUT) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               mpistat, &  ! MPI library call return status
               nfstat,  &  ! netCDF library call return status
               i,       &  ! Loop counter
               j,       &  !     "
               k,       &  !     "
               ndims       ! Number of dimensions in this variable

    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    INTEGER,DIMENSION(1:5) :: var_dimlens
    INTEGER :: tindex
    LOGICAL :: hastime   ! Whether this variable has a time axis
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
    REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: databuf

!DL test, 1 octobre 2012
    !if (lwp) write(100,*), filename, varname, time

    ! Initializations
    status = NCDF_NOERR
    hastime = .FALSE.
    var_dimlens = 0
    j = 0
    if(time > 0) then
       tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
    else
       tindex = -time
    endif

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading 3D array from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Open netCDF file and get info
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Opening file'
!JC       CALL FLUSH
    END IF

    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Getting info from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, varid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Determine if this variable contains a time axis
    nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
         ndims=ndims, dimids=var_dimids)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    DO i = 1, ndims
       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
            name=dname, len=k)
!       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
!            name=dname, len=k)
       dname=TRIM(dname)
       var_dimlens(i-j) = k
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
!DL added test on Time, october 1sr, 2012
       IF((dname=='time_counter') .OR. (dname=='time') .OR. (dname=='Time')) THEN
          hastime = .TRUE.
          j = 1
          EXIT
       END IF
    END DO
    ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3)))

    ! Read the appropriate subdomain for each processor
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Reading data from file'
!JC       CALL FLUSH
    END IF

    IF(hastime .EQV. .TRUE.) THEN
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, 1, tindex/), &
            (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), 1 /))
    ELSE
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, 1 /), &
            (/ var_dimlens(1), var_dimlens(2), var_dimlens(3) /))
    END IF
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    DO k = 1, var_dimlens(3)
       DO j = 1, nlcj
          DO i = 1, nlci
             data(i, j, k) = real(databuf(mig(i), mjg(j), k),8)
          END DO
       END DO
    END DO


    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read3d - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF


    DEALLOCATE(databuf)

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_read3d

  ! The inverse of ncdf_write4d. Reads a 4-D array from a netCDF file.
  ! Each processor only reads its local subdomain.
  ! filename - file to read from
  ! varname - variable to read from
  ! data - 4-D array to put the data into (must be a REAL)
  ! dsz - Dimension SiZe - size of the 4th dimension in the variable
  ! time - time index to read (If the target variable in the netCDF file is
  !        actually a 4-D variable with a time axis, ignored otherwise)
  ! status - return status of the subroutine
  ! NOT EXTENSIVELY TESTED - MAY STILL BE BUGGY

! =================================
! V V V V V       NL#2     V V V V V V
  ! Update 2013.08.23 : Nicolas Lambert
  !      => Add of a buffer (databuf) in the reading of the data to solve a bug about the reading 
  !         of the data in local domain. Now, the code look like the ncdf_read3d but with 1 more 
  !         dimension. ncdf_write4d maybe need to be corect too... 


  SUBROUTINE ncdf_read4d(filename, varname, data, dsz, time, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: time, dsz
    REAL(wp),DIMENSION(:,:,:,:),INTENT(OUT) :: data
    INTEGER,INTENT(OUT) :: status
    
    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               mpistat, &  ! MPI library call return status
               nfstat,  &  ! netCDF library call return status
               i,k,j,jn,     &  ! Loop counter
               ndims       ! Number of dimensions in this variable

    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    INTEGER,DIMENSION(1:5) :: var_dimlens
    LOGICAL :: hastime   ! Whether this variable has a time axis
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
    REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: databuf
    INTEGER :: tindex
    
    ! Initializations
    status = NCDF_NOERR
    hastime = .FALSE.
    var_dimlens = 0
    j = 0

    if(time > 0) then
       tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
    else
       tindex = -time
    endif
!    tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
    
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading 4D array from file:', trim(filename)
!JC       CALL FLUSH
    END IF
    
    ! Open netCDF file and get info
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Opening file'
!JC       CALL FLUSH
    END IF

    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Getting info from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, varid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Determine if this variable contains a time axis
    nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
         ndims=ndims, dimids=var_dimids)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    DO i = 1, ndims
       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
            name=dname, len=k)
!       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
!            name=dname, len=k)
       dname=TRIM(dname)
       var_dimlens(i-j) = k
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
!DL added test on Time, october 1sr, 2012
       IF((dname=='time_counter') .OR. (dname=='time') .OR. (dname=='Time')) THEN
          hastime = .TRUE.
          j = 1
          EXIT
       END IF
    END DO
    ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3),1:var_dimlens(4)))
    
    ! Got the data, write it to the proper location in the netCDF file
    IF(hastime .EQV. .TRUE.) THEN
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, 1, 1, tindex/), &
            (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), var_dimlens(4), 1 /)) 
    ELSE
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, 1, 1 /), &
            (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), var_dimlens(4) /)) 
    END IF
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    DO jn = 1,var_dimlens(4)
    DO k = 1, var_dimlens(3)
       DO j = 1, nlcj
          DO i = 1, nlci
             data(i, j, k, jn) = real(databuf(mig(i), mjg(j), k ,jn),8)
          END DO
       END DO
    END DO
    END DO  
    
    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read4d - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_read4d
! ^ ^ ^ ^ ^       NL#2   ^ ^ ^ ^ ^   

  ! ncdf_create_file creates an empty netCDF dataset with the specified filename
  ! filename - name of the file to be created
  ! status - return status
  SUBROUTINE ncdf_create_file(filename, status)
    IMPLICIT NONE
    
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,   &  ! file ID
               nfstat, &  ! netCDF call return status
               mpistat    ! MPI call return status

    ! Initialization
    status = NCDF_NOERR

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_create_file - Creating:', trim(filename)
!JC       CALL FLUSH
    END IF

    IF(nproc == 0) THEN
       nfstat = nf90_create(filename, NF90_CLOBBER, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_create_file - Done creating'
!JC       CALL FLUSH
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_create_file

  ! ncdf_create_dim creates a new dimension in an existing dataset
  ! filename - name of the file to modified
  ! dimname - name of the dimension to create
  ! dim_len - length of the new dimension
  ! status - return status
  SUBROUTINE ncdf_create_dim(filename, dimname, dim_len, status)
    IMPLICIT NONE

    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, dimname
    INTEGER,INTENT(IN) :: dim_len
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! file ID
               nfstat,  &  ! netCDF call return status
               mpistat, &  ! MPI call return status
               dimid,   &  ! New dimension ID 
               dimlen

    ! Initialization
    status = NCDF_NOERR
    dimlen = dim_len

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_create_dim - Creating dimension in file:', trim(filename)
!JC       CALL FLUSH
    END IF

    IF(nproc == 0) THEN
       IF(dimlen == -1) THEN
          dimlen = NF90_UNLIMITED
       END IF
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_redef(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_def_dim(ncid, dimname, dimlen, dimid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_create_dim - Done creating dimension in file'
!JC       CALL FLUSH
    END IF
    
    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_create_dim

  ! ncdf_create_var creates a new variable in an existing dataset
  ! filename - name of the file to modified
  ! varname - name of the variable to create
  ! dimname - an array containing the names (in order) of the dimensions
  !           the new variable will use
  ! datatype - type of the variable (NCDF_FLOAT or NCDF_DOUBLE)
  ! status - return status
  SUBROUTINE ncdf_create_var(filename, varname, dimnames, datatype, status)
    IMPLICIT NONE
    
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: dimnames
    INTEGER,INTENT(IN) :: datatype
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! file ID
               nfstat,  &  ! netCDF call return status
               mpistat, &  ! MPI call return status
               varid,   &  ! New variable ID
               i
    INTEGER,DIMENSION(1:NF90_MAX_VAR_DIMS) :: dimids

    ! Initialization
    status = NCDF_NOERR
    dimids = -1

    IF((datatype /= NCDF_FLOAT) .AND. (datatype /= NCDF_DOUBLE)) THEN
       status = NCDF_ARERR
       RETURN
    END IF

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_create_var - Creating variable in file:', trim(filename)
!JC       CALL FLUSH
    END IF

    IF(nproc == 0) THEN
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       ! Get dimension IDs
       DO i = 1, SIZE(dimnames)
          nfstat = nf90_inq_dimid(ncid, dimnames(i), dimids(i))
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
             RETURN
          END IF
       END DO
       ! Create variable
       nfstat = nf90_redef(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_def_var(ncid, varname, datatype, &
            dimids(1:SIZE(dimnames)), varid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_create_var - Done creating variable in file'
!JC       CALL FLUSH
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_create_var

  ! ncdf_put_att_int creates a new attribute in an existing dataset
  ! filename - name of the file to modified
  ! varname - name of the variable to attach the attribute to (or "GLOBAL")
  ! attname - name of the new attribute
  ! attval - value the attribute should have
  ! status - return status
  SUBROUTINE ncdf_put_att_int(filename, varname, attname, attval, status)
    IMPLICIT NONE

    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname
    INTEGER,INTENT(IN) :: attval
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! file ID
               nfstat,  &  ! netCDF call return status
               mpistat, &  ! MPI call return status
               varid       ! ID of variable to attach attribute to

    ! Initialization
    status = NCDF_NOERR

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_int - Creating attribute in file:', trim(filename)
!JC       CALL FLUSH
    END IF

    IF(nproc == 0) THEN
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       ! Get variable ID
       IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN
          varid = NF90_GLOBAL
       ELSE
          nfstat = nf90_inq_varid(ncid, varname, varid)
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
             RETURN
          END IF
       END IF
       ! Write attribute
       nfstat = nf90_redef(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_put_att(ncid, varid, attname, attval)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_int - Done creating attribute in file'
!JC       CALL FLUSH
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_put_att_int

  ! ncdf_put_att_real creates a new attribute in an existing dataset
  ! filename - name of the file to modified
  ! varname - name of the variable to attach the attribute to (or "GLOBAL")
  ! attname - name of the new attribute
  ! attval - value the attribute should have
  ! status - return status
  SUBROUTINE ncdf_put_att_real(filename, varname, attname, attval, status)
    IMPLICIT NONE

    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname
    REAL(wp),INTENT(IN) :: attval
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! file ID
               nfstat,  &  ! netCDF call return status
               mpistat, &  ! MPI call return status
               varid       ! ID of variable to attach attribute to

    ! Initialization
    status = NCDF_NOERR

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_real - Creating attribute in file:', trim(filename)
!JC       CALL FLUSH
    END IF

    IF(nproc == 0) THEN
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       ! Get variable ID
       IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN
          varid = NF90_GLOBAL
       ELSE
          nfstat = nf90_inq_varid(ncid, varname, varid)
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
             RETURN
          END IF
       END IF
       ! Write attribute
       nfstat = nf90_redef(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_put_att(ncid, varid, attname, attval)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_real - Done creating attribute in file'
!JC       CALL FLUSH
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_put_att_real

  ! ncdf_put_att_char creates a new attribute in an existing dataset
  ! filename - name of the file to modified
  ! varname - name of the variable to attach the attribute to (or "GLOBAL")
  ! attname - name of the new attribute
  ! attval - value the attribute should have
  ! status - return status
  SUBROUTINE ncdf_put_att_char(filename, varname, attname, attval, status)
    IMPLICIT NONE

    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname, attname, attval
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! file ID
               nfstat,  &  ! netCDF call return status
               mpistat, &  ! MPI call return status
               varid       ! ID of variable to attach attribute to

    ! Initialization
    status = NCDF_NOERR

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_char - Creating attribute in file:', trim(filename)
!JC       CALL FLUSH
    END IF

    IF(nproc == 0) THEN
       nfstat = nf90_open(filename, nf90_write, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       ! Get variable ID
       IF((varname == 'global') .OR. (varname == 'GLOBAL')) THEN
          varid = NF90_GLOBAL
       ELSE
          nfstat = nf90_inq_varid(ncid, varname, varid)
          IF(nfstat /= nf90_noerr) THEN
             status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
             RETURN
          END IF
       END IF
       ! Write attribute
       nfstat = nf90_redef(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_put_att(ncid, varid, attname, attval)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_put_att_char - Done creating attribute in file'
!JC       CALL FLUSH
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_put_att_char

  ! ncdf_errstr takes an error code returned by a subroutine in this module
  ! and returns a string indicating the meaning of the error code
  SUBROUTINE ncdf_errstr(errcode, errstr)
    INTEGER,INTENT(IN) :: errcode
    CHARACTER(LEN=80),INTENT(OUT) :: errstr
    
    errstr = 'No error'
    
    SELECT CASE (errcode)
    CASE (0)
       errstr = 'No error'
    CASE(1)
       errstr = 'netCDF error'
    CASE(2)
       errstr = 'MPI error'
    CASE(3)
       errstr = 'Invalid arguments to subroutine'
    CASE(4)
       errstr = 'Other/unrecognized error'
    CASE DEFAULT
       errstr = 'Unknown error code'
    END SELECT

  END SUBROUTINE ncdf_errstr


!!DB: 2007.12.11
  ! ncdf_create_file_aveTSUV builds a special OPA output file with most of the default
  ! dimensions, variables and attributes for TSUV fields
  ! It is designed for average outputs of useful variables:
  ! u, v, T, S, taux, tauy, ssh, ...
!!2008.06.13 w added ----> 1 more dimid; 2 more varids
!float vovecrtz (time_counter, depthw, y, x)
!                char units = "m/s"
!                float missing_value = 1.000000e+20
!                float valid_min = 1.000000e+20
!                float valid_max = -1.000000e+20
!                char long_name = "Vertical Velocity"
!                char short_name = "vovecrtz"
!                char online_operation = "ave(x)"
!                char axis = "TZYX"
!                float interval_operation = 4.800000e+02
!                float interval_write = 8.640000e+04
!                char associate = "time_counter depthw nav_lat nav_lon"


  SUBROUTINE ncdf_create_file_aveTSUV(filename, op_type, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    CHARACTER(LEN=*),INTENT(IN) :: op_type
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
!    INTEGER,DIMENSION(1:4) :: dimids
    INTEGER,DIMENSION(1:5) :: dimids
!    INTEGER,DIMENSION(1:12) :: varids
!AD:    INTEGER,DIMENSION(1:14) :: varids
    INTEGER,DIMENSION(1:17) :: varids    
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since    
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    INTEGER :: int_opp, &              ! Operation interval
               int_wri                 ! Write interval
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    
    ! Initializations
    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0
    int_opp = nwrite * rdt
    int_wri = nwrite * rdt

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating default U output file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Creating file'
!JC          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(filename, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining dimensions in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'depthu', jpkdta, dimids(2))
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
       nfstat = nf90_def_dim(ncid, 'depthw', jpkdta, dimids(5))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_u - Defining variables in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(4), dimids(3) /), &
            varids(2))
       nfstat = nf90_def_var(ncid, 'depthu', nf90_float, &
            (/ dimids(2) /), &
            varids(3))
       nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
            (/ dimids(1) /), &
            varids(4))
       nfstat = nf90_def_var(ncid, 'vozocrtx', nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(5))
       nfstat = nf90_def_var(ncid, 'sozotaux', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(6))

       nfstat = nf90_def_var(ncid, 'vomecrty', nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(7))
       nfstat = nf90_def_var(ncid, 'sometauy', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(8))

       nfstat = nf90_def_var(ncid, 'votemper', nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(9))
       nfstat = nf90_def_var(ncid, 'vosaline', nf90_float, &
            (/ dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(10))
       nfstat = nf90_def_var(ncid, 'sossheig', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(11))
        nfstat = nf90_def_var(ncid, 'soicecov', nf90_float, &
            (/ dimids(4), dimids(3), dimids(1) /), &
            varids(12))

       nfstat = nf90_def_var(ncid, 'depthw', nf90_float, &
            (/ dimids(5) /), &
            varids(13))
       nfstat = nf90_def_var(ncid, 'vovecrtz', nf90_float, &
            (/ dimids(4), dimids(3), dimids(5), dimids(1) /), &
            varids(14))
!AD: add 
       ! ndate (ndastp)
       nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
            (/ dimids(1) /), &
            varids(15))

       ! ndate (model_time)
       nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
            (/ dimids(1) /), &
            varids(16))

       ! kt 
       nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
            (/ dimids(1) /), &
            varids(17))


       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_aveTSUV - Writing attributes in file'
!JC          CALL FLUSH
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! depthu
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 3.046773E00)
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 5.875141E03)
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'depthu')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical U levels')

       ! depthw
       nfstat = nf90_put_att(ncid, varids(13), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(13), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(13), 'valid_min', 0.00E00)
       nfstat = nf90_put_att(ncid, varids(13), 'valid_max', 5.875141E03)
       nfstat = nf90_put_att(ncid, varids(13), 'title', 'depthw')
       nfstat = nf90_put_att(ncid, varids(13), 'long_name', 'Vertical W levels')

       ! time_counter
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! vozocrtx
       nfstat = nf90_put_att(ncid, varids(5), 'units', 'm/s')
       nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Zonal Current')
       nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'vozocrtx')
       nfstat = nf90_put_att(ncid, varids(5), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(5), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter depthu nav_lat nav_lon')

       ! sozotaux
       nfstat = nf90_put_att(ncid, varids(6), 'units', 'N/m2')
       nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'Wind Stress along i-axis')
       nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'sozotaux')
       nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(6), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter nav_lat nav_lon')

       ! vomecrty
       nfstat = nf90_put_att(ncid, varids(7), 'units', 'm/s')
       nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(7), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(7), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(7), 'long_name', 'Meridional Current')
       nfstat = nf90_put_att(ncid, varids(7), 'short_name', 'vomecrty')
       nfstat = nf90_put_att(ncid, varids(7), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(7), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(7), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(7), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(7), 'associate', 'time_counter depthv nav_lat nav_lon')

       ! sometauy
       nfstat = nf90_put_att(ncid, varids(8), 'units', 'N/m2')
       nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(8), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(8), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(8), 'long_name', 'Wind Stress along j-axis')
       nfstat = nf90_put_att(ncid, varids(8), 'short_name', 'sometauy')
       nfstat = nf90_put_att(ncid, varids(8), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(8), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(8), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(8), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(8), 'associate', 'time_counter nav_lat nav_lon')

       ! votemper
       nfstat = nf90_put_att(ncid, varids(9), 'units', 'C')
       nfstat = nf90_put_att(ncid, varids(9), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(9), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(9), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(9), 'long_name', 'Temperature')
       nfstat = nf90_put_att(ncid, varids(9), 'short_name', 'votemper')
       nfstat = nf90_put_att(ncid, varids(9), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(9), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(9), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(9), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(9), 'associate', 'time_counter deptht nav_lat nav_lon')

       ! vosaline
       nfstat = nf90_put_att(ncid, varids(10), 'units', 'PSU')
       nfstat = nf90_put_att(ncid, varids(10), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(10), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(10), 'long_name', 'Salinity')
       nfstat = nf90_put_att(ncid, varids(10), 'short_name', 'vosaline')
       nfstat = nf90_put_att(ncid, varids(10), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(10), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(10), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(10), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(10), 'associate', 'time_counter deptht nav_lat nav_lon')

       ! sossheig
       nfstat = nf90_put_att(ncid, varids(11), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(11), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(11), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(11), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(11), 'long_name', 'Sea Surface Height')
       nfstat = nf90_put_att(ncid, varids(11), 'short_name', 'sossheig')
       nfstat = nf90_put_att(ncid, varids(11), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(11), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(11), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(11), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(11), 'associate', 'time_counter nav_lat nav_lon')

       ! soicecov
       nfstat = nf90_put_att(ncid, varids(12), 'units', '[0,1]')
       nfstat = nf90_put_att(ncid, varids(12), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(12), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(12), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(12), 'long_name', 'Ice Cover')
       nfstat = nf90_put_att(ncid, varids(12), 'short_name', 'soicecov')
       nfstat = nf90_put_att(ncid, varids(12), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(12), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(12), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(12), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(12), 'associate', 'time_counter nav_lat nav_lon')

       ! vovecrtz
       nfstat = nf90_put_att(ncid, varids(14), 'units', 'm/s')
       nfstat = nf90_put_att(ncid, varids(14), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(14), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(14), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(14), 'long_name', 'Vertical Velocity')
       nfstat = nf90_put_att(ncid, varids(14), 'short_name', 'vovecrtz')
       nfstat = nf90_put_att(ncid, varids(14), 'online_operation', TRIM(op_type))
       nfstat = nf90_put_att(ncid, varids(14), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(14), 'interval_operation', int_opp)
       nfstat = nf90_put_att(ncid, varids(14), 'interval_write', int_wri)
       nfstat = nf90_put_att(ncid, varids(14), 'associate', 'time_counter depthw nav_lat nav_lon')

       ! ndate (ndastp)
       nfstat = nf90_put_att(ncid, varids(15), 'units', '=nyear*10000+nmonth*100+nday')
       nfstat = nf90_put_att(ncid, varids(15), 'long_name', 'time step date in year/month/day aammjj')

       ! ndate (model_time)
       nfstat = nf90_put_att(ncid, varids(16), 'long_name', &
            'time step date (when output is writen) in year/month/day aammjj (decimal day)')
       nfstat = nf90_put_att(ncid, varids(16), 'units', '=nyear*10000+nmonth*100+nday')
       nfstat = nf90_put_att(ncid, varids(16), 'formula1', 'nyear  =   model_time / 10000')       
       nfstat = nf90_put_att(ncid, varids(16), 'formula2', & 
            'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
       nfstat = nf90_put_att(ncid, varids(16), 'formula3', & 
            'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           


       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
      
       
       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_file_TSUV - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_file_aveTSUV

!!DB 2008.05.22
  ! ncdf_create_ice_restart builds a standard OPA lim ice restart file with all the default
  ! dimensions, variables and attributes. Note that this is a single restart
  ! file, and therefore incompatible with the default IOIPSL-based restart
  ! routines (though it holds the same data)
  SUBROUTINE ncdf_create_ice_restart(filename, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER,INTENT(OUT) :: status
    
    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
         varid,   &  ! ID of netCDF variable to be written to
         nfstat,  &  ! netCDF library call return status
         mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:9) :: dimids
    INTEGER,DIMENSION(1:36) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since
    CHARACTER(LEN=100) :: tstp_since
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    INTEGER :: ji

    
    ! Initializations
    status = NCDF_NOERR
    
    ! Initializations
    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE (UNIT=tstp_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating default restart file:', trim(filename)
!JC       CALL FLUSH
    END IF
    
    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Creating file'
!JC          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(filename, nf90_clobber, ncid)
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining dimensions in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(1))
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(2))
       nfstat = nf90_def_dim(ncid, 'z', 1, dimids(3))
       nfstat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimids(4))
       nfstat = nf90_def_dim(ncid, 'x_a', 1, dimids(5))
       nfstat = nf90_def_dim(ncid, 'y_a', 1, dimids(6))
       nfstat = nf90_def_dim(ncid, 'z_a', 2, dimids(7))
       nfstat = nf90_def_dim(ncid, 'z_b', 3, dimids(8))
       nfstat = nf90_def_dim(ncid, 'z_c', 35, dimids(9))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining variables in file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(1), dimids(2) /), &
            varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(1), dimids(2) /), &
            varids(2))
       nfstat = nf90_def_var(ncid, 'nav_lev', nf90_float, &
            (/ dimids(3) /), &
            varids(3))
       nfstat = nf90_def_var(ncid, 'time', nf90_float, &
            (/ dimids(4) /), &
            varids(4))
       nfstat = nf90_def_var(ncid, 'time_steps', nf90_float, &
            (/ dimids(4) /), &
            varids(5))
       nfstat = nf90_def_var(ncid, 'info', nf90_float, &
            (/ dimids(7) /), &
            varids(6))

       nfstat = nf90_def_var(ncid, 'hicif'  , nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(7))
       nfstat = nf90_def_var(ncid, 'hsnif'   , nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(8))
       nfstat = nf90_def_var(ncid, 'frld'    , nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(9))
       nfstat = nf90_def_var(ncid, 'sist'    , nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(10))
       nfstat = nf90_def_var(ncid, 'tbif'   ,  nf90_float, &
            (/ dimids(1), dimids(2), dimids(8), dimids(4) /), &
            varids(11))
       nfstat = nf90_def_var(ncid, 'u_ice'  ,  nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(12))
       nfstat = nf90_def_var(ncid, 'v_ice'  ,  nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(13))
       nfstat = nf90_def_var(ncid, 'gtaux'  ,  nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(14))
       nfstat = nf90_def_var(ncid, 'gtauy'  ,  nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(15))
       nfstat = nf90_def_var(ncid, 'qstoif' ,  nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(16))
       nfstat = nf90_def_var(ncid, 'fsbbq'  ,  nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(17))
       nfstat = nf90_def_var(ncid, 'moment' ,  nf90_float, &
            (/ dimids(1), dimids(2), dimids(9), dimids(4) /), &
            varids(18))
# if defined key_coupled
      nfstat = nf90_def_var(ncid, 'albege' , , nf90_float, &
            (/ dimids(1), dimids(2), dimids(3), dimids(4) /), &
            varids(19))
# endif

       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Writing attributes in file'
!JC          CALL FLUSH
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -1.800000E2)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', 1.800000E2)
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', -9.000000E1)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 9.000000E1)
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')

       ! nav_lev
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'model_levels')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 0.0)
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 0.0)
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Model levels')

       ! time
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'Time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', '0001-JUL-01 00:00:00')

       ! time_steps
       nfstat = nf90_put_att(ncid, varids(5), 'units', TRIM(tstp_since))
       nfstat = nf90_put_att(ncid, varids(5), 'title', 'Time steps')
       nfstat = nf90_put_att(ncid, varids(5), 'tstep_sec', 1.877760E7)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Time step axis')
       nfstat = nf90_put_att(ncid, varids(5), 'time_origin', TRIM(t_origin))

       ! info
       nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
!DB as per OLD code
       do ji = 7, 18
          nfstat = nf90_put_att(ncid, varids(ji), 'missing_value', 1.000000E20)
       enddo
#if defined key_coupled
          nfstat = nf90_put_att(ncid, varids(19), 'missing_value', 1.000000E20)
#endif

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL , 'Conventions', 'GDT 1.2')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number_total', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_global', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_local', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_first', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_last', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_type','box')

       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Closing file'
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
       
  END SUBROUTINE ncdf_create_ice_restart

!! CN 10/14/2008 - finds number of dimensions in the specified variable
  SUBROUTINE ncdf_get_num_dims(filename, varname, ndims, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(OUT) :: ndims, status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               id,      &  ! ID of netCDF variable of interest
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status


    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Getting number of dimensions in variable from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Initializations
    status = NCDF_NOERR

    ! Open netCDF file and get info
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Opening file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Getting varid from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, id)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Getting number of dims in variable from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inquire_variable(ncid, id, ndims=ndims)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_get_num_dims - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
 
    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_get_num_dims

!!DB
!!Get size of dimension -- common use is to get size of unlimited dimension
!!(i.e. # records)
!!Note that all processors open the file and get the dimension size
  SUBROUTINE ncdf_get_dim_size(filename, dimname, len, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, dimname
    INTEGER,INTENT(OUT) :: len, status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               id,      &  ! ID of netCDF dimension of interest
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    CHARACTER(LEN=nf90_max_name) dname


    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Getting dim size from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Initializations
    status = NCDF_NOERR

    ! Open netCDF file and get info
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_get dim_size - Opening file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Getting dimid from file'
!JC       CALL FLUSH
    END IF
    nfstat = NF90_INQ_DIMID (ncid, dimname, id)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Getting dimsize from file'
!JC       CALL FLUSH
    END IF
    nfstat = NF90_INQUIRE_DIMENSION (ncid, id, len=len)

    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_get_dim_size - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
 
    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_get_dim_size

  
!>>> A.D:  2008.09.29
  
    SUBROUTINE ncdf_readdate(filename, varname, cdata, index, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: index
    CHARACTER(LEN=19),INTENT(OUT) :: cdata
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:5) :: count

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading single date from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Initializations
    status = NCDF_NOERR
    count(1) = 1

    ! Open netCDF file and get info
       
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Opening file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Getting info from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, varid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Reading data from file'
!JC       CALL FLUSH
    END IF
    !       nfstat = nf90_get_var(ncid, varid, data, (/ index /))
    nfstat = nf90_get_var(ncid, varid, cdata, (/ 1,index /), (/ 19,1 /))           
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_readdate - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_readdate

  
!<<< A.D:  
  
SUBROUTINE ncdf_read2d_global(filename, varname, data, time, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: time
    REAL(wp),DIMENSION(:,:),INTENT(OUT) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               mpistat, &  ! MPI library call return status
               nfstat,  &  ! netCDF library call return status
               i,       &  ! Loop counter
               j,       &  !     "
               k,       &  !     "
               ndims       ! Number of dimensions in this variable

    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    INTEGER,DIMENSION(1:5) :: var_dimlens
    INTEGER :: tindex
    LOGICAL :: hastime   ! Whether this variable has a time axis
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
    REAL*4,ALLOCATABLE,DIMENSION(:,:) :: databuf

    ! Initializations
    status = NCDF_NOERR
    hastime = .FALSE.
    var_dimlens = 0
    j = 0
    if(time > 0) then
       tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
    else
       tindex = -time
    endif


    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading 2D array from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Open netCDF file and get info
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Opening file'
!JC       CALL FLUSH
    END IF

    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
      WRITE(100,*) 'NCDF ERROR open:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Getting info from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, varid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
      WRITE(100,*) 'NCDF ERROR inq_varID:',nf90_strerror(nfstat)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Determine if this variable contains a time axis
    nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
         ndims=ndims, dimids=var_dimids)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
      WRITE(100,*) 'NCDF ERROR inquire_variable:',nf90_strerror(nfstat)
       RETURN
    END IF
    DO i = 1, ndims
       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
            name=dname, len=k)
       dname=TRIM(dname)
       var_dimlens(i-j) = k
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          WRITE(100,*) 'NCDF ERROR inquire_dimension:',nf90_strerror(nfstat), trim(filename)

          RETURN
       END IF
       IF((dname=='time_counter') .OR. (dname=='time')) THEN
          hastime = .TRUE.
          j = 1
          EXIT
       END IF
    END DO
    ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2)))

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Reading data from file'
!JC       CALL FLUSH
    END IF
    IF(hastime .EQV. .TRUE.) THEN
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, tindex/), &
            (/ var_dimlens(1), var_dimlens(2), 1 /))
    ELSE
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1 /), &
            (/ var_dimlens(1), var_dimlens(2) /))
    END IF
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;
      WRITE(100,*) 'NCDF ERROR get_var:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    data = real(databuf,8)

    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read2d_global - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    DEALLOCATE(databuf)

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_read2d_global


!!DB
  ! Reads a 3-D array from a netCDF file.
  ! Returns global array to wach processor
  ! filename - file to read from
  ! varname - variable to read from
  ! data - 3-D array to put the data into (must be a REAL)
  ! time - time index to read (If the target variable in the netCDF file is
  !        actually a 4-D variable with a time axis, ignored otherwise)
  ! status - return status of the subroutine
  SUBROUTINE ncdf_read3d_global(filename, varname, data, time, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: time
    REAL(wp),DIMENSION(:,:,:),INTENT(OUT) :: data
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               mpistat, &  ! MPI library call return status
               nfstat,  &  ! netCDF library call return status
               i,       &  ! Loop counter
               j,       &  !     "
               k,       &  !     "
               ndims       ! Number of dimensions in this variable

    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    INTEGER,DIMENSION(1:5) :: var_dimlens
    INTEGER :: tindex
    LOGICAL :: hastime   ! Whether this variable has a time axis
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
    REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: databuf

    ! Initializations
    status = NCDF_NOERR
    hastime = .FALSE.
    var_dimlens = 0
    j = 0
    if(time > 0) then
       tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
    else
       tindex = -time
    endif

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading 3D array from file:', trim(filename)
!JC       CALL FLUSH
    END IF

    ! Open netCDF file and get info
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Opening file'
!JC       CALL FLUSH
    END IF

    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Getting info from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, varid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Determine if this variable contains a time axis
    nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
         ndims=ndims, dimids=var_dimids)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    DO i = 1, ndims
       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
            name=dname, len=k)
       dname=TRIM(dname)
       var_dimlens(i-j) = k
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       IF((dname=='time_counter') .OR. (dname=='time')) THEN
          hastime = .TRUE.
          j = 1
          EXIT
       END IF
    END DO
    ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3)))

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Reading data from file'
!JC       CALL FLUSH
    END IF
    IF(hastime .EQV. .TRUE.) THEN
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, 1, tindex/), &
            (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), 1 /))
    ELSE
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, 1 /), &
            (/ var_dimlens(1), var_dimlens(2), var_dimlens(3) /))
    END IF
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    data = real(databuf,8)

    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read3d_global - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    DEALLOCATE(databuf)

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_read3d_global

  SUBROUTINE ncdf_read4d_global(filename, varname, data, dsz, time, status)
    IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename, varname
    INTEGER,INTENT(IN) :: time, dsz
    REAL(wp),DIMENSION(:,:,:,:),INTENT(OUT) :: data
    INTEGER,INTENT(OUT) :: status
    
    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               mpistat, &  ! MPI library call return status
               nfstat,  &  ! netCDF library call return status
               i,j,k,       &  ! Loop counter
               ndims       ! Number of dimensions in this variable

    INTEGER,DIMENSION(1:5) :: var_dimids ! Array to hold dimension ids
    INTEGER,DIMENSION(1:5) :: var_dimlens
    LOGICAL :: hastime   ! Whether this variable has a time axis
    CHARACTER(LEN=NF90_MAX_NAME) :: dname ! dimension name 
    INTEGER :: tindex
    REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: databuf
    
    ! Initializations
    status = NCDF_NOERR
    hastime = .FALSE.
    j=0

    if(time > 0) then
       tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
    else
       tindex = -time
    endif
!    tindex = NINT(REAL(time / nwrite)) ! Calculate time index for read
    
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Reading 4D array from file:', trim(filename)
!JC       CALL FLUSH
    END IF
    
    ! Open netCDF file and get info
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Opening file'
!JC       CALL FLUSH
    END IF

    nfstat = nf90_open(filename, nf90_nowrite, ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat),filename
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    IF(DEBUG_OUT .EQV. .TRUE. .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Getting info from file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_inq_varid(ncid, varname, varid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Determine if this variable contains a time axis
    nfstat = nf90_inquire_variable(ncid=ncid, varid=varid, &
         ndims=ndims, dimids=var_dimids)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    DO i = 1, ndims
       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=var_dimids(i), &
            name=dname, len=k)
!       nfstat = nf90_inquire_dimension(ncid=ncid, dimid=i, &
!            name=dname)
       dname=TRIM(dname)
       var_dimlens(i-j) = k
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
          RETURN
       END IF
       IF(dname=='time_counter') THEN
          hastime = .TRUE.
          j = 1
          EXIT
       END IF
    END DO
    ALLOCATE(databuf(1:var_dimlens(1),1:var_dimlens(2),1:var_dimlens(3),1:var_dimlens(4)))
    
    ! Read data
    IF(hastime .EQV. .TRUE.) THEN
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, 1, 1, tindex/), &
            (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), var_dimlens(4), 1 /))
    ELSE
       nfstat = nf90_get_var(ncid, varid, databuf, &
            (/ 1, 1, 1, 1 /), &
            (/ var_dimlens(1), var_dimlens(2), var_dimlens(3), var_dimlens(4) /))
    END IF
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF

    data = real(databuf,8)
    DEALLOCATE(databuf)
    
    ! All done, close up the netCDF file
    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: ncdf_read4d_global - Closing file'
!JC       CALL FLUSH
    END IF
    nfstat = nf90_close(ncid)
    IF(nfstat /= nf90_noerr) THEN
       status = NCDF_NFERR;write(100,*),'ERROR:',nf90_strerror(nfstat), trim(filename)
          print*,'ERROR:',nf90_strerror(nfstat),trim(filename)
       RETURN
    END IF
    
    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_read4d_global



END MODULE lib_ncdf
