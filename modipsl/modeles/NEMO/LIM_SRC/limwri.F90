!!DB -- 2009.09.04 -- key_diadimg eliminated

!!DB 2008.05.21
!!Modified to create and output 1 .nc file
!!Tested vs IOIPSL output = identical
!!05.27 -- Added averaging of arrays 

MODULE limwri
   !!======================================================================
   !!                     ***  MODULE  limwri  ***
   !!         Ice diagnostics :  write ice output files
   !!======================================================================
   !!----------------------------------------------------------------------
   !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limwri.F90,v 1.6 2005/12/12 14:18:00 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim'                                     LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_wri      : write of the diagnostics variables in ouput file 
   !!   lim_wri_init : initialization and namelist read
   !!----------------------------------------------------------------------
   !! * Modules used
   USE ioipsl
   USE dianam    ! build name of file (routine)
   USE phycst
   USE dom_oce
   USE daymod
   USE in_out_manager
   USE ice_oce         ! ice variables
   USE flx_oce
   USE dom_ice
   USE ice
   USE iceini
   USE lbclnk

!!DB
   USE lib_ncdf        ! netCDF I/O library

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC lim_wri        ! routine called by lim_step.F90

   !! * Module variables
   INTEGER, PARAMETER ::   &  !:
      jpnoumax = 40             !: maximum number of variable for ice output
   INTEGER  ::                                &
      noumef                                     ! number of fields
   REAL(wp)           , DIMENSION(jpnoumax) ::  &
      cmulti ,                                &  ! multiplicative constant
      cadd                                       ! additive constant
   CHARACTER(len = 35), DIMENSION(jpnoumax) ::  &
      titn                                       ! title of the field
   CHARACTER(len = 8 ), DIMENSION(jpnoumax) ::  &
      nam                                        ! name of the field
   CHARACTER(len = 8 ), DIMENSION(jpnoumax) ::  &
      uni                                        ! unit of the field
   INTEGER            , DIMENSION(jpnoumax) ::  &
      nc                                         ! switch for saving field ( = 1 ) or not ( = 0 )

   REAL(wp)  ::            &  ! constant values
      epsi16 = 1.e-16   ,  &
      zzero  = 0.e0     ,  &
      zone   = 1.e0
   !!-------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_wri
      !!-------------------------------------------------------------------
      !!  This routine computes the average of some variables and write it
      !!  on the ouput files.
      !!  ATTENTION cette routine n'est valable que si le pas de temps est
      !!  egale a une fraction entiere de 1 jours.
      !!  Diff 1-D 3-D : suppress common also included in etat
      !!                 suppress cmoymo 11-18
      !!  modif : 03/06/98
      !!-------------------------------------------------------------------

!!DB
     

      !! * Local variables
      REAL(wp),DIMENSION(1) ::   zdept
      
      REAL(wp) :: &
         zsto, zsec, zjulian,zout, &
         zindh,zinda,zindb,  &
         ztmu
      REAL(wp), DIMENSION(jpi,jpj,jpnoumax) :: &
         zcmo
      REAL(wp), DIMENSION(jpi,jpj) ::  &
         zfield
      INTEGER ::  ji, jj, jf   ! dummy loop indices
!!DB
      INTEGER :: status, rec_num, kt
      CHARACTER(len = 80) :: fname
!!DB -- related to time averaging of ice fields
!!NB: ice_field_index could be allocated, but instead just give it a v.big dimension
      INTEGER, SAVE ::  num_ice_fields = 0
      INTEGER, DIMENSION(100), SAVE :: ice_field_index = 0
      REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: ave_ice_field
      LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.    ! Use IOIPSL subroutines for output

      CHARACTER(len = 40)  :: &
         clhstnam, clop

      INTEGER , SAVE ::      &
         nice, nhorid, ndim, niter, ndepid
      INTEGER , DIMENSION( jpij ) , SAVE ::  &
         ndex51  
      !!-------------------------------------------------------------------
      

!!DB
      fname = trim(cexper)//'_icemod.nc'

      IF ( numit == nstart ) THEN 
!!!DB
         call ncdf_create_ice_file(status)

         CALL lim_wri_init                 

!!DB -- determine # of ice output variables and allocate memory
!!Eventually only do if clop="ave(x)"
         do jf = 1, noumef
            if ( nc(jf) == 1 ) then 
               num_ice_fields = num_ice_fields + 1
               ice_field_index(jf) = num_ice_fields
            endif
         enddo
         ALLOCATE(ave_ice_field(jpi,jpj,num_ice_fields))
         ave_ice_field = 0.0



         
         !---5----|----5----|----5----|----5----|----5----|----5----|----5----|72
         !  1) INITIALIZATIONS.                                                 |
         !-----------------------------------------------------------------------
         
         !-- essai NetCDF
         
         zsto= rdt_ice
!!Chris         clop     = "ave(only(x))"      !ibug  namelist parameter a ajouter
         clop = "ave(x)"
         zout = nwrite * rdt_ice / nfice
         zsec     = 0.
         niter    = 0
         zdept(1) = 0.

         if(USE_IOIPSL) then
            CALL ymds2ju ( nyear, nmonth, nday, zsec, zjulian )
            CALL dia_nam ( clhstnam, nwrite, 'icemod' )
            CALL histbeg ( clhstnam, jpi, glamt, jpj, gphit, 1, jpi, 1, jpj, 0, zjulian, rdt_ice, nhorid, nice , domain_id=nidom)
            CALL histvert( nice, "deptht", "Vertical T levels", "m", 1, zdept, ndepid)
            CALL wheneq  ( jpij , tmask(:,:,1), 1, 1., ndex51, ndim)
            
            DO jf = 1, noumef
               IF ( nc(jf) == 1 ) THEN
                  CALL histdef( nice, nam(jf), titn(jf), uni(jf), jpi, jpj   &
                       , nhorid, 1, 1, 1, -99, 32, clop, zsto, zout )
               ENDIF
            END DO
            CALL histend(nice)
         endif

      ENDIF
      
      !---5----|----5----|----5----|----5----|----5----|----5----|----5----|72
      !--2. Computation of instantaneous values                                         |
      !-----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_wri : write ice outputs in NetCDF files at time : ', nyear, nmonth, nday, numit
         WRITE(numout,*) '~~~~~~~ '
      ENDIF

      !-- calculs des valeurs instantanees
      
      zcmo(:,:, 1:jpnoumax ) = 0.e0 
      DO jj = 2 , jpjm1
         DO ji = 2 , jpim1
            zindh  = MAX( zzero , SIGN( zone , hicif(ji,jj) * (1.0 - frld(ji,jj) ) - 0.10 ) )
            zinda  = MAX( zzero , SIGN( zone , ( 1.0 - frld(ji,jj) ) - 0.10 ) )
            zindb  = zindh * zinda
            ztmu   = MAX( 0.5 * zone , ( tmu(ji,jj) + tmu(ji+1,jj) + tmu(ji,jj+1) + tmu(ji+1,jj+1) ) ) 
            zcmo(ji,jj,1)  = hsnif (ji,jj)
            zcmo(ji,jj,2)  = hicif (ji,jj)
            zcmo(ji,jj,3)  = hicifp(ji,jj)
            zcmo(ji,jj,4)  = frld  (ji,jj)
            zcmo(ji,jj,5)  = sist  (ji,jj)
            zcmo(ji,jj,6)  = fbif  (ji,jj)
            zcmo(ji,jj,7)  = zindb * (  u_ice(ji,jj  ) * tmu(ji,jj  ) + u_ice(ji+1,jj  ) * tmu(ji+1,jj  )   &
                                      + u_ice(ji,jj+1) * tmu(ji,jj+1) + u_ice(ji+1,jj+1) * tmu(ji+1,jj+1) ) &
                                  / ztmu 

            zcmo(ji,jj,8)  = zindb * (  v_ice(ji,jj  ) * tmu(ji,jj  ) + v_ice(ji+1,jj  ) * tmu(ji+1,jj  )   &
                                      + v_ice(ji,jj+1) * tmu(ji,jj+1) + v_ice(ji+1,jj+1) * tmu(ji+1,jj+1) ) &
                                  / ztmu
            zcmo(ji,jj,9)  = sst_io(ji,jj)
            zcmo(ji,jj,10) = sss_io(ji,jj)

            zcmo(ji,jj,11) = fnsolar(ji,jj) + fsolar(ji,jj)
            zcmo(ji,jj,12) = fsolar (ji,jj)
            zcmo(ji,jj,13) = fnsolar(ji,jj)
            ! See thersf for the coefficient
            zcmo(ji,jj,14) = - fsalt(ji,jj) * rday * ( sss_io(ji,jj) + epsi16 ) / soce
            zcmo(ji,jj,15) = gtaux(ji,jj)
            zcmo(ji,jj,16) = gtauy(ji,jj)
            zcmo(ji,jj,17) = ( 1.0 - frld(ji,jj) ) * qsr_ice (ji,jj) + frld(ji,jj) * qsr_oce (ji,jj)
            zcmo(ji,jj,18) = ( 1.0 - frld(ji,jj) ) * qnsr_ice(ji,jj) + frld(ji,jj) * qnsr_oce(ji,jj)
            zcmo(ji,jj,19) = sprecip(ji,jj)
         END DO
      END DO
                
      !
      ! ecriture d'un fichier netcdf
      !
      niter = niter + 1
      DO jf = 1 , noumef
         DO jj = 1 , jpj
            DO ji = 1 , jpi
               zfield(ji,jj) = zcmo(ji,jj,jf) * cmulti(jf) + cadd(jf)
            END DO
         END DO
         
         IF ( jf == 7  .OR. jf == 8  .OR. jf == 11 .OR. jf == 12 .OR. jf == 15 .OR.   &
            jf == 23 .OR. jf == 24 .OR. jf == 16 ) THEN 
            CALL lbc_lnk( zfield, 'T', -1. )
         ELSE 
            CALL lbc_lnk( zfield, 'T',  1. )
         ENDIF
         
         if(USE_IOIPSL) then
            IF ( nc(jf) == 1 ) CALL histwrite( nice, nam(jf), niter, zfield, ndim, ndex51 )
         endif

!!DB output routines; 
         if(nc(jf) == 1) then  !!accumulate field
                  ave_ice_field(:,:,ice_field_index(jf)) =  & 
                       ave_ice_field(:,:,ice_field_index(jf)) + zfield(:,:) 
         endif

         if(mod(niter,int(nwrite/nfice)) ==  0) then

            if(nc(jf) == 1) then
               kt = nfice*niter + nit000 - 1
               rec_num = niter/int(nwrite/nfice)
               zfield(:,:) = ave_ice_field(:,:,ice_field_index(jf))/(float(nwrite)/float(nfice))
               CALL ncdf_write(fname, 'time_counter', REAL(kt * rdt), rec_num, status)
               CALL ncdf_write(fname, nam(jf), zfield, -rec_num, status)
               ave_ice_field(:,:,ice_field_index(jf)) = 0.0

!!AD/DB 2009.09.30
               CALL ncdf_write(fname, 'ndastp',REAL(ndastp), rec_num, status)
               CALL ncdf_write(fname, 'model_time_step',REAL(kt), rec_num, status)
               CALL ncdf_write(fname, 'model_time',model_time, rec_num, status)

            endif

         endif
         
      END DO
      
      if(USE_IOIPSL) then
         IF ( ( nfice * niter + nit000 - 1 ) >= nitend ) THEN
            CALL histclo( nice ) 
         ENDIF
      endif

   END SUBROUTINE lim_wri

   
   SUBROUTINE lim_wri_init
      !!-------------------------------------------------------------------
      !!                    ***   ROUTINE lim_wri_init  ***
      !!                
      !! ** Purpose :   ???
      !!
      !! ** Method  : Read the namicewri namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicewri
      !!
      !! history :
      !!  8.5  ! 03-08 (C. Ethe) original code
      !!-------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   nf      ! ???

      TYPE FIELD 
         CHARACTER(len = 35) :: ztitle 
         CHARACTER(len = 8 ) :: zname          
         CHARACTER(len = 8 ) :: zunit
         INTEGER             :: znc   
         REAL                :: zcmulti 
         REAL                :: zcadd        
      END TYPE FIELD

      TYPE(FIELD) ::  &
         field_1 , field_2 , field_3 , field_4 , field_5 , field_6 ,   &
         field_7 , field_8 , field_9 , field_10, field_11, field_12,   &
         field_13, field_14, field_15, field_16, field_17, field_18,   &
         field_19

      TYPE(FIELD) , DIMENSION(jpnoumax) :: zfield

      NAMELIST/namiceout/ noumef, &
         field_1 , field_2 , field_3 , field_4 , field_5 , field_6 ,   &
         field_7 , field_8 , field_9 , field_10, field_11, field_12,   &
         field_13, field_14, field_15, field_16, field_17, field_18,   &
         field_19
      !!-------------------------------------------------------------------


      ! Read Namelist namicewri
      REWIND ( numnam_ice )
      READ   ( numnam_ice  , namiceout )
      zfield(1)  = field_1
      zfield(2)  = field_2
      zfield(3)  = field_3
      zfield(4)  = field_4
      zfield(5)  = field_5
      zfield(6)  = field_6
      zfield(7)  = field_7
      zfield(8)  = field_8
      zfield(9)  = field_9
      zfield(10) = field_10
      zfield(11) = field_11
      zfield(12) = field_12
      zfield(13) = field_13
      zfield(14) = field_14
      zfield(15) = field_15
      zfield(16) = field_16
      zfield(17) = field_17
      zfield(18) = field_18
      zfield(19) = field_19
      
      DO nf = 1, noumef
         titn  (nf) = zfield(nf)%ztitle
         nam   (nf) = zfield(nf)%zname
         uni   (nf) = zfield(nf)%zunit
         nc    (nf) = zfield(nf)%znc
         cmulti(nf) = zfield(nf)%zcmulti
         cadd  (nf) = zfield(nf)%zcadd
      END DO

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_wri_init : Ice parameters for outputs'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '    number of fields to be stored         noumef = ', noumef
         WRITE(numout,*) '           title                            name     unit      Saving (1/0) ',   &
            &            '    multiplicative constant       additive constant '
         DO nf = 1 , noumef         
            WRITE(numout,*) '   ', titn(nf), '   ', nam(nf),'      ', uni(nf),'  ', nc(nf),'        ', cmulti(nf),   &
               '        ', cadd(nf)
         END DO
      ENDIF
            
   END SUBROUTINE lim_wri_init


!!DB
!!Create 1 ice file equivalent to the one created in lim_wri

 SUBROUTINE ncdf_create_ice_file(status)
    IMPLICIT NONE
    ! Subroutine argument declarations

    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:4) :: dimids
!!DB -- unsure how to dimension varids, so choose a large number 
    INTEGER,DIMENSION(1:50) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since    
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    CHARACTER (len=80) :: fname
    CHARACTER(LEN=20) :: op_type, varname
    INTEGER :: int_opp, &              ! Operation interval
               int_wri                 ! Write interval
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    INTEGER :: jf, varnum

! Get required info from namelist_ice
    if ( numit == nstart ) THEN 
       CALL lim_wri_init 
    endif
    
    
    ! Initializations
    op_type = "ave(x)"       !!default, dangerous, requires routine to do this  
    fname = trim(cexper)//'_icemod.nc'

    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0

!!DB see above
    int_opp = rdt_ice
    int_wri = nwrite * rdt_ice/nfice

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating output file:', fname
!JC       CALL FLUSH
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Creating file:', fname
!JC          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(fname, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Defining dimensions in file:', fname
!JC          CALL FLUSH
       END IF
       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'z', 1, dimids(2))          !!ice model
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Defining variables in file:', fname
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
       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Writing attributes in file:', fname
!JC          CALL FLUSH
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
!       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
!       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', minval(glamt))
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', maxval(glamt))
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
!       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
!       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', minval(gphit))
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', maxval(gphit))
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! deptht
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 0.0)      !for ice model
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 0.0)      !for ice model
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical T levels')

       ! time_counter
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(fname))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'MISC', 'DB-created file')

!!Loop over variables in namelist_ice as done above
       varnum = 4      !start at # of vars defined above
       do jf = 1, noumef
          if ( nc(jf) == 1 ) THEN

             varnum = varnum + 1
             nfstat = nf90_def_var(ncid, nam(jf), nf90_float, &
                  (/ dimids(4), dimids(3), dimids(1) /), varids(varnum))

             nfstat = nf90_put_att(ncid, varids(varnum), 'units', uni(jf))
             nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
             nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
             nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
             nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', TRIM(op_type))
             nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', titn(jf))
             nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', nam(jf))
             nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TYX')
             nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
             nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
             nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

          endif
       enddo



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
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Closing file:', fname
!JC          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF
    END IF



! On first write step only
    CALL ncdf_write(fname, 'nav_lat', gphit, -1, status)
    CALL ncdf_write(fname, 'nav_lon', glamt, -1, status)


    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_ice_file




#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module          NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_wri          ! Empty routine
   END SUBROUTINE lim_wri
#endif

   !!======================================================================
END MODULE limwri
