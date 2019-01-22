


PGF90 (Version     10.5)          09/23/2010  09:50:00      page 1

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

Filename: limwri.F90

(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) 
(    3) !!DB 2008.05.21
(    4) !!Modified to create and output 1 .nc file
(    5) !!Tested vs IOIPSL output = identical
(    6) !!05.27 -- Added averaging of arrays 
(    7) 
(    8) MODULE limwri
(    9)    !!======================================================================
(   10)    !!                     ***  MODULE  limwri  ***
(   11)    !!         Ice diagnostics :  write ice output files
(   12)    !!======================================================================
(   13)    !!----------------------------------------------------------------------
(   14)    !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
(   15)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limwri.F90,v 1.6 2005/12/12 14:18:00 opalod Exp $
(   16)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   17)    !!----------------------------------------------------------------------
(   18) #if defined key_ice_lim
(   19)    !!----------------------------------------------------------------------
(   20)    !!   'key_ice_lim'                                     LIM sea-ice model
(   21)    !!----------------------------------------------------------------------
(   22)    !!   lim_wri      : write of the diagnostics variables in ouput file 
(   23)    !!   lim_wri_init : initialization and namelist read
(   24)    !!----------------------------------------------------------------------
(   25)    !! * Modules used
(   26)    USE ioipsl
(   27)    USE dianam    ! build name of file (routine)
(   28)    USE phycst
(   29)    USE dom_oce
(   30)    USE daymod
(   31)    USE in_out_manager
(   32)    USE ice_oce         ! ice variables
(   33)    USE flx_oce
(   34)    USE dom_ice
(   35)    USE ice
(   36)    USE iceini
(   37)    USE lbclnk
(   38) 
(   39) !!DB
(   40)    USE lib_ncdf        ! netCDF I/O library
(   41) 
(   42)    IMPLICIT NONE
(   43)    PRIVATE
(   44) 
(   45)    !! * Accessibility
(   46)    PUBLIC lim_wri        ! routine called by lim_step.F90






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 2

(   47) 
(   48)    !! * Module variables
(   49)    INTEGER, PARAMETER ::   &  !:
(   50)       jpnoumax = 40             !: maximum number of variable for ice output
(   51)    INTEGER  ::                                &
(   52)       noumef                                     ! number of fields
(   53)    REAL(wp)           , DIMENSION(jpnoumax) ::  &
(   54)       cmulti ,                                &  ! multiplicative constant
(   55)       cadd                                       ! additive constant
(   56)    CHARACTER(len = 35), DIMENSION(jpnoumax) ::  &
(   57)       titn                                       ! title of the field
(   58)    CHARACTER(len = 8 ), DIMENSION(jpnoumax) ::  &
(   59)       nam                                        ! name of the field
(   60)    CHARACTER(len = 8 ), DIMENSION(jpnoumax) ::  &
(   61)       uni                                        ! unit of the field
(   62)    INTEGER            , DIMENSION(jpnoumax) ::  &
(   63)       nc                                         ! switch for saving field ( = 1 ) or not ( = 0 )
(   64) 
(   65)    REAL(wp)  ::            &  ! constant values
(   66)       epsi16 = 1.e-16   ,  &
(   67)       zzero  = 0.e0     ,  &
(   68)       zone   = 1.e0
(   69)    !!-------------------------------------------------------------------
(   70) 
(   71) CONTAINS
(   72) 
(   73)    SUBROUTINE lim_wri
(   74)       !!-------------------------------------------------------------------
(   75)       !!  This routine computes the average of some variables and write it
(   76)       !!  on the ouput files.
(   77)       !!  ATTENTION cette routine n'est valable que si le pas de temps est
(   78)       !!  egale a une fraction entiere de 1 jours.
(   79)       !!  Diff 1-D 3-D : suppress common also included in etat
(   80)       !!                 suppress cmoymo 11-18
(   81)       !!  modif : 03/06/98
(   82)       !!-------------------------------------------------------------------
(   83) 
(   84) !!DB
(   85)      
(   86) 
(   87)       !! * Local variables
(   88)       REAL(wp),DIMENSION(1) ::   zdept
(   89)       
(   90)       REAL(wp) :: &
(   91)          zsto, zsec, zjulian,zout, &
(   92)          zindh,zinda,zindb,  &
(   93)          ztmu
(   94)       REAL(wp), DIMENSION(jpi,jpj,jpnoumax) :: &
(   95)          zcmo
(   96)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   97)          zfield
(   98)       INTEGER ::  ji, jj, jf   ! dummy loop indices
(   99) !!DB
(  100)       INTEGER :: status, rec_num, kt
(  101)       CHARACTER(len = 80) :: fname
(  102) !!DB -- related to time averaging of ice fields
(  103) !!NB: ice_field_index could be allocated, but instead just give it a v.big dimension
(  104)       INTEGER, SAVE ::  num_ice_fields = 0






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 3

(  105)       INTEGER, DIMENSION(100), SAVE :: ice_field_index = 0
(  106)       REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: ave_ice_field
(  107)       LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.    ! Use IOIPSL subroutines for output
(  108) 
(  109)       CHARACTER(len = 40)  :: &
(  110)          clhstnam, clop
(  111) 
(  112)       INTEGER , SAVE ::      &
(  113)          nice, nhorid, ndim, niter, ndepid
(  114)       INTEGER , DIMENSION( jpij ) , SAVE ::  &
(  115)          ndex51  
(  116)       !!-------------------------------------------------------------------
(  117)       
(  118) 
(  119) !!DB
(  120)       fname = trim(cexper)//'_icemod.nc'
(  121) 
(  122)       IF ( numit == nstart ) THEN 
(  123) !!!DB
(  124)          call ncdf_create_ice_file(status)
(  125) 
(  126)          CALL lim_wri_init                 
(  127) 
(  128) !!DB -- determine # of ice output variables and allocate memory
(  129) !!Eventually only do if clop="ave(x)"
(  130)          do jf = 1, noumef
(  131)             if ( nc(jf) == 1 ) then 
(  132)                num_ice_fields = num_ice_fields + 1
(  133)                ice_field_index(jf) = num_ice_fields
(  134)             endif
(  135)          enddo
(  136)          ALLOCATE(ave_ice_field(jpi,jpj,num_ice_fields))
(  137)          ave_ice_field = 0.0
(  138) 
(  139) 
(  140) 
(  141)          
(  142)          !---5----|----5----|----5----|----5----|----5----|----5----|----5----|72
(  143)          !  1) INITIALIZATIONS.                                                 |
(  144)          !-----------------------------------------------------------------------
(  145)          
(  146)          !-- essai NetCDF
(  147)          
(  148)          zsto= rdt_ice
(  149) !!Chris         clop     = "ave(only(x))"      !ibug  namelist parameter a ajouter
(  150)          clop = "ave(x)"
(  151)          zout = nwrite * rdt_ice / nfice
(  152)          zsec     = 0.
(  153)          niter    = 0
(  154)          zdept(1) = 0.
(  155) 
(  156)          if(USE_IOIPSL) then
(  157)             CALL ymds2ju ( nyear, nmonth, nday, zsec, zjulian )
(  158)             CALL dia_nam ( clhstnam, nwrite, 'icemod' )
(  159)             CALL histbeg ( clhstnam, jpi, glamt, jpj, gphit, 1, jpi, 1, jpj, 0, zjulian, rdt_ice, nhorid, nice , domain_id=n
(  160)             CALL histvert( nice, "deptht", "Vertical T levels", "m", 1, zdept, ndepid)
(  161)             CALL wheneq  ( jpij , tmask(:,:,1), 1, 1., ndex51, ndim)
(  162)             






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 4

(  163)             DO jf = 1, noumef
(  164)                IF ( nc(jf) == 1 ) THEN
(  165)                   CALL histdef( nice, nam(jf), titn(jf), uni(jf), jpi, jpj   &
(  166)                        , nhorid, 1, 1, 1, -99, 32, clop, zsto, zout )
(  167)                ENDIF
(  168)             END DO
(  169)             CALL histend(nice)
(  170)          endif
(  171) 
(  172)       ENDIF
(  173)       
(  174)       !---5----|----5----|----5----|----5----|----5----|----5----|----5----|72
(  175)       !--2. Computation of instantaneous values                                         |
(  176)       !-----------------------------------------------------------------------
(  177) 
(  178)       IF(lwp) THEN
(  179)          WRITE(numout,*)
(  180)          WRITE(numout,*) 'lim_wri : write ice outputs in NetCDF files at time : ', nyear, nmonth, nday, numit
(  181)          WRITE(numout,*) '~~~~~~~ '
(  182)       ENDIF
(  183) 
(  184)       !-- calculs des valeurs instantanees
(  185)       
(  186)       zcmo(:,:, 1:jpnoumax ) = 0.e0 
(  187)       DO jj = 2 , jpjm1
(  188)          DO ji = 2 , jpim1
(  189)             zindh  = MAX( zzero , SIGN( zone , hicif(ji,jj) * (1.0 - frld(ji,jj) ) - 0.10 ) )
(  190)             zinda  = MAX( zzero , SIGN( zone , ( 1.0 - frld(ji,jj) ) - 0.10 ) )
(  191)             zindb  = zindh * zinda
(  192)             ztmu   = MAX( 0.5 * zone , ( tmu(ji,jj) + tmu(ji+1,jj) + tmu(ji,jj+1) + tmu(ji+1,jj+1) ) ) 
(  193)             zcmo(ji,jj,1)  = hsnif (ji,jj)
(  194)             zcmo(ji,jj,2)  = hicif (ji,jj)
(  195)             zcmo(ji,jj,3)  = hicifp(ji,jj)
(  196)             zcmo(ji,jj,4)  = frld  (ji,jj)
(  197)             zcmo(ji,jj,5)  = sist  (ji,jj)
(  198)             zcmo(ji,jj,6)  = fbif  (ji,jj)
(  199)             zcmo(ji,jj,7)  = zindb * (  u_ice(ji,jj  ) * tmu(ji,jj  ) + u_ice(ji+1,jj  ) * tmu(ji+1,jj  )   &
(  200)                                       + u_ice(ji,jj+1) * tmu(ji,jj+1) + u_ice(ji+1,jj+1) * tmu(ji+1,jj+1) ) &
(  201)                                   / ztmu 
(  202) 
(  203)             zcmo(ji,jj,8)  = zindb * (  v_ice(ji,jj  ) * tmu(ji,jj  ) + v_ice(ji+1,jj  ) * tmu(ji+1,jj  )   &
(  204)                                       + v_ice(ji,jj+1) * tmu(ji,jj+1) + v_ice(ji+1,jj+1) * tmu(ji+1,jj+1) ) &
(  205)                                   / ztmu
(  206)             zcmo(ji,jj,9)  = sst_io(ji,jj)
(  207)             zcmo(ji,jj,10) = sss_io(ji,jj)
(  208) 
(  209)             zcmo(ji,jj,11) = fnsolar(ji,jj) + fsolar(ji,jj)
(  210)             zcmo(ji,jj,12) = fsolar (ji,jj)
(  211)             zcmo(ji,jj,13) = fnsolar(ji,jj)
(  212)             ! See thersf for the coefficient
(  213)             zcmo(ji,jj,14) = - fsalt(ji,jj) * rday * ( sss_io(ji,jj) + epsi16 ) / soce
(  214)             zcmo(ji,jj,15) = gtaux(ji,jj)
(  215)             zcmo(ji,jj,16) = gtauy(ji,jj)
(  216)             zcmo(ji,jj,17) = ( 1.0 - frld(ji,jj) ) * qsr_ice (ji,jj) + frld(ji,jj) * qsr_oce (ji,jj)
(  217)             zcmo(ji,jj,18) = ( 1.0 - frld(ji,jj) ) * qnsr_ice(ji,jj) + frld(ji,jj) * qnsr_oce(ji,jj)
(  218)             zcmo(ji,jj,19) = sprecip(ji,jj)
(  219)          END DO
(  220)       END DO






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 5

(  221)                 
(  222)       !
(  223)       ! ecriture d'un fichier netcdf
(  224)       !
(  225)       niter = niter + 1
(  226)       DO jf = 1 , noumef
(  227)          DO jj = 1 , jpj
(  228)             DO ji = 1 , jpi
(  229)                zfield(ji,jj) = zcmo(ji,jj,jf) * cmulti(jf) + cadd(jf)
(  230)             END DO
(  231)          END DO
(  232)          
(  233)          IF ( jf == 7  .OR. jf == 8  .OR. jf == 11 .OR. jf == 12 .OR. jf == 15 .OR.   &
(  234)             jf == 23 .OR. jf == 24 .OR. jf == 16 ) THEN 
(  235)             CALL lbc_lnk( zfield, 'T', -1. )
(  236)          ELSE 
(  237)             CALL lbc_lnk( zfield, 'T',  1. )
(  238)          ENDIF
(  239)          
(  240)          if(USE_IOIPSL) then
(  241)             IF ( nc(jf) == 1 ) CALL histwrite( nice, nam(jf), niter, zfield, ndim, ndex51 )
(  242)          endif
(  243) 
(  244) !!DB output routines; 
(  245)          if(nc(jf) == 1) then  !!accumulate field
(  246)                   ave_ice_field(:,:,ice_field_index(jf)) =  & 
(  247)                        ave_ice_field(:,:,ice_field_index(jf)) + zfield(:,:) 
(  248)          endif
(  249) 
(  250)          if(mod(niter,int(nwrite/nfice)) ==  0) then
(  251) 
(  252)             if(nc(jf) == 1) then
(  253)                kt = nfice*niter + nit000 - 1
(  254)                rec_num = niter/int(nwrite/nfice)
(  255)                zfield(:,:) = ave_ice_field(:,:,ice_field_index(jf))/(float(nwrite)/float(nfice))
(  256)                CALL ncdf_write(fname, 'time_counter', REAL(kt * rdt), rec_num, status)
(  257)                CALL ncdf_write(fname, nam(jf), zfield, -rec_num, status)
(  258)                ave_ice_field(:,:,ice_field_index(jf)) = 0.0
(  259) 
(  260) !!AD/DB 2009.09.30
(  261)                CALL ncdf_write(fname, 'ndastp',REAL(ndastp), rec_num, status)
(  262)                CALL ncdf_write(fname, 'model_time_step',REAL(kt), rec_num, status)
(  263)                CALL ncdf_write(fname, 'model_time',model_time, rec_num, status)
(  264) 
(  265)             endif
(  266) 
(  267)          endif
(  268)          
(  269)       END DO
(  270)       
(  271)       if(USE_IOIPSL) then
(  272)          IF ( ( nfice * niter + nit000 - 1 ) >= nitend ) THEN
(  273)             CALL histclo( nice ) 
(  274)          ENDIF
(  275)       endif
(  276) 
(  277)    END SUBROUTINE lim_wri
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 6

(  279)    
(  280)    SUBROUTINE lim_wri_init
(  281)       !!-------------------------------------------------------------------
(  282)       !!                    ***   ROUTINE lim_wri_init  ***
(  283)       !!                
(  284)       !! ** Purpose :   ???
(  285)       !!
(  286)       !! ** Method  : Read the namicewri namelist and check the parameter 
(  287)       !!       values called at the first timestep (nit000)
(  288)       !!
(  289)       !! ** input   :   Namelist namicewri
(  290)       !!
(  291)       !! history :
(  292)       !!  8.5  ! 03-08 (C. Ethe) original code
(  293)       !!-------------------------------------------------------------------
(  294)       !! * Local declarations
(  295)       INTEGER ::   nf      ! ???
(  296) 
(  297)       TYPE FIELD 
(  298)          CHARACTER(len = 35) :: ztitle 
(  299)          CHARACTER(len = 8 ) :: zname          
(  300)          CHARACTER(len = 8 ) :: zunit
(  301)          INTEGER             :: znc   
(  302)          REAL                :: zcmulti 
(  303)          REAL                :: zcadd        
(  304)       END TYPE FIELD
(  305) 
(  306)       TYPE(FIELD) ::  &
(  307)          field_1 , field_2 , field_3 , field_4 , field_5 , field_6 ,   &
(  308)          field_7 , field_8 , field_9 , field_10, field_11, field_12,   &
(  309)          field_13, field_14, field_15, field_16, field_17, field_18,   &
(  310)          field_19
(  311) 
(  312)       TYPE(FIELD) , DIMENSION(jpnoumax) :: zfield
(  313) 
(  314)       NAMELIST/namiceout/ noumef, &
(  315)          field_1 , field_2 , field_3 , field_4 , field_5 , field_6 ,   &
(  316)          field_7 , field_8 , field_9 , field_10, field_11, field_12,   &
(  317)          field_13, field_14, field_15, field_16, field_17, field_18,   &
(  318)          field_19
(  319)       !!-------------------------------------------------------------------
(  320) 
(  321) 
(  322)       ! Read Namelist namicewri
(  323)       REWIND ( numnam_ice )
(  324)       READ   ( numnam_ice  , namiceout )
(  325)       zfield(1)  = field_1
(  326)       zfield(2)  = field_2
(  327)       zfield(3)  = field_3
(  328)       zfield(4)  = field_4
(  329)       zfield(5)  = field_5
(  330)       zfield(6)  = field_6
(  331)       zfield(7)  = field_7
(  332)       zfield(8)  = field_8
(  333)       zfield(9)  = field_9
(  334)       zfield(10) = field_10
(  335)       zfield(11) = field_11
(  336)       zfield(12) = field_12






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 7

(  337)       zfield(13) = field_13
(  338)       zfield(14) = field_14
(  339)       zfield(15) = field_15
(  340)       zfield(16) = field_16
(  341)       zfield(17) = field_17
(  342)       zfield(18) = field_18
(  343)       zfield(19) = field_19
(  344)       
(  345)       DO nf = 1, noumef
(  346)          titn  (nf) = zfield(nf)%ztitle
(  347)          nam   (nf) = zfield(nf)%zname
(  348)          uni   (nf) = zfield(nf)%zunit
(  349)          nc    (nf) = zfield(nf)%znc
(  350)          cmulti(nf) = zfield(nf)%zcmulti
(  351)          cadd  (nf) = zfield(nf)%zcadd
(  352)       END DO
(  353) 
(  354)       IF(lwp) THEN
(  355)          WRITE(numout,*)
(  356)          WRITE(numout,*) 'lim_wri_init : Ice parameters for outputs'
(  357)          WRITE(numout,*) '~~~~~~~~~~~~'
(  358)          WRITE(numout,*) '    number of fields to be stored         noumef = ', noumef
(  359)          WRITE(numout,*) '           title                            name     unit      Saving (1/0) ',   &
(  360)             &            '    multiplicative constant       additive constant '
(  361)          DO nf = 1 , noumef         
(  362)             WRITE(numout,*) '   ', titn(nf), '   ', nam(nf),'      ', uni(nf),'  ', nc(nf),'        ', cmulti(nf),   &
(  363)                '        ', cadd(nf)
(  364)          END DO
(  365)       ENDIF
(  366)             
(  367)    END SUBROUTINE lim_wri_init
(  368) 
(  369) 
(  370) !!DB
(  371) !!Create 1 ice file equivalent to the one created in lim_wri
(  372) 
(  373)  SUBROUTINE ncdf_create_ice_file(status)
(  374)     IMPLICIT NONE
(  375)     ! Subroutine argument declarations
(  376) 
(  377)     INTEGER,INTENT(OUT) :: status
(  378) 
(  379)     ! Local declarations
(  380)     INTEGER :: ncid,    &  ! netCDF file ID
(  381)                varid,   &  ! ID of netCDF variable to be written to
(  382)                nfstat,  &  ! netCDF library call return status
(  383)                mpistat     ! MPI library call return status
(  384)     INTEGER,DIMENSION(1:4) :: dimids
(  385) !!DB -- unsure how to dimension varids, so choose a large number 
(  386)     INTEGER,DIMENSION(1:50) :: varids
(  387)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  388)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  389)     CHARACTER(LEN=100) :: sec_since    
(  390)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  391)     CHARACTER (len=80) :: fname
(  392)     CHARACTER(LEN=20) :: op_type, varname
(  393)     INTEGER :: int_opp, &              ! Operation interval
(  394)                int_wri                 ! Write interval






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 8

(  395)     CHARACTER(LEN=3),PARAMETER :: &
(  396)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  397)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  398)     INTEGER :: jf, varnum
(  399) 
(  400) ! Get required info from namelist_ice
(  401)     if ( numit == nstart ) THEN 
(  402)        CALL lim_wri_init 
(  403)     endif
(  404)     
(  405)     
(  406)     ! Initializations
(  407)     op_type = "ave(x)"       !!default, dangerous, requires routine to do this  
(  408)     fname = trim(cexper)//'_icemod.nc'
(  409) 
(  410)     status = NCDF_NOERR
(  411)     CALL ioget_calendar(cal_type)
(  412)     CALL ioget_timestamp(timestamp)
(  413)     WRITE (UNIT=sec_since, &
(  414)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  415)          &  nyear,nmonth,nday,0, 0, 0
(  416)     WRITE(t_origin, &
(  417)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  418)          &   nyear,months(nmonth),nday,0,0,0
(  419) 
(  420) !!DB see above
(  421)     int_opp = rdt_ice
(  422)     int_wri = nwrite * rdt_ice/nfice
(  423) 
(  424)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  425)        WRITE(100,*) 'NCDF DEBUG: Creating output file:', fname
(  426)        CALL FLUSH
(  427)     END IF
(  428) 
(  429)     ! Only processor 0 does anything
(  430)     IF(nproc == 0) THEN
(  431)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  432)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Creating file:', fname
(  433)           CALL FLUSH
(  434)        END IF
(  435)        ! Create the file
(  436)        nfstat = nf90_create(fname, nf90_clobber, ncid)
(  437)        IF(nfstat /= nf90_noerr) THEN
(  438)           status = NCDF_NFERR
(  439)           RETURN
(  440)        END IF
(  441)        
(  442)        ! Define dimensions
(  443)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  444)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Defining dimensions in file:', fname
(  445)           CALL FLUSH
(  446)        END IF
(  447)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  448)        nfstat = nf90_def_dim(ncid, 'z', 1, dimids(2))          !!ice model
(  449)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
(  450)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
(  451)        
(  452)        ! Define variables






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 9

(  453)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  454)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Defining variables in file:', fname
(  455)           CALL FLUSH
(  456)        END IF
(  457)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  458)             (/ dimids(4), dimids(3) /), &
(  459)             varids(1))
(  460)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  461)             (/ dimids(4), dimids(3) /), &
(  462)             varids(2))
(  463)        nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &
(  464)             (/ dimids(2) /), &
(  465)             varids(3))
(  466)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  467)             (/ dimids(1) /), &
(  468)             varids(4))
(  469)        
(  470)        ! Add attributes
(  471)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  472)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Writing attributes in file:', fname
(  473)           CALL FLUSH
(  474)        END IF
(  475)        ! nav_lon
(  476)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  477) !       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  478) !       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  479)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', minval(glamt))
(  480)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', maxval(glamt))
(  481)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  482)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  483) 
(  484)        ! nav_lat
(  485)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  486) !       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  487) !       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  488)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', minval(gphit))
(  489)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', maxval(gphit))
(  490)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  491)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  492) 
(  493)        ! deptht
(  494)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  495)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  496)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 0.0)      !for ice model
(  497)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 0.0)      !for ice model
(  498)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
(  499)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical T levels')
(  500) 
(  501)        ! time_counter
(  502)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  503)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  504)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  505)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
(  506)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
(  507) 
(  508)        ! global
(  509)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  510)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(fname))






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 10

(  511)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
(  512)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  513)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  514)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'MISC', 'DB-created file')
(  515) 
(  516) !!Loop over variables in namelist_ice as done above
(  517)        varnum = 4      !start at # of vars defined above
(  518)        do jf = 1, noumef
(  519)           if ( nc(jf) == 1 ) THEN
(  520) 
(  521)              varnum = varnum + 1
(  522)              nfstat = nf90_def_var(ncid, nam(jf), nf90_float, &
(  523)                   (/ dimids(4), dimids(3), dimids(1) /), varids(varnum))
(  524) 
(  525)              nfstat = nf90_put_att(ncid, varids(varnum), 'units', uni(jf))
(  526)              nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
(  527)              nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
(  528)              nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
(  529)              nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', TRIM(op_type))
(  530)              nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', titn(jf))
(  531)              nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', nam(jf))
(  532)              nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TYX')
(  533)              nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
(  534)              nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
(  535)              nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')
(  536) 
(  537)           endif
(  538)        enddo
(  539) 
(  540) 
(  541) 
(  542) !AD/DB: add new time-related variables
(  543)        varnum = varnum + 1
(  544)        ! ndate (ndastp)
(  545)        nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
(  546)             (/ dimids(1) /), &
(  547)             varids(varnum))
(  548)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
(  549)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'time step date in year/month/day aammjj')
(  550) 
(  551)        varnum = varnum + 1
(  552)        ! ndate (model_time)
(  553)        nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
(  554)             (/ dimids(1) /), &
(  555)             varids(varnum))
(  556)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', &
(  557)             'time step date (when output is writen) in year/month/day aammjj (decimal day)')
(  558)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
(  559)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula1', 'nyear  =   model_time / 10000')       
(  560)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula2', & 
(  561)             'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
(  562)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula3', & 
(  563)             'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           
(  564) 
(  565)        varnum = varnum + 1
(  566)        ! kt 
(  567)        nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
(  568)             (/ dimids(1) /), &






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 11

(  569)             varids(varnum))
(  570) 
(  571) 
(  572)        ! Close file
(  573)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  574)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Closing file:', fname
(  575)           CALL FLUSH
(  576)        END IF
(  577)        nfstat = nf90_close(ncid)
(  578)        IF(nfstat /= nf90_noerr) THEN
(  579)           status = NCDF_NFERR
(  580)           RETURN
(  581)        END IF
(  582)     END IF
(  583) 
(  584) 
(  585) 
(  586) ! On first write step only
(  587)     CALL ncdf_write(fname, 'nav_lat', gphit, -1, status)
(  588)     CALL ncdf_write(fname, 'nav_lon', glamt, -1, status)
(  589) 
(  590) 
(  591)     ! Sync up processors before returning from subroutine
(  592)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
(  593)     IF(mpistat /= 0) THEN
(  594)        status = NCDF_MPERR
(  595)        RETURN
(  596)     END IF
(  597)     
(  598)   END SUBROUTINE ncdf_create_ice_file
(  599) 
(  600) 
(  601) 
(  602) 
(  603) #else
(  604)    !!----------------------------------------------------------------------
(  605)    !!   Default option :         Empty module          NO LIM sea-ice model
(  606)    !!----------------------------------------------------------------------
(  607) CONTAINS
(  608)    SUBROUTINE lim_wri          ! Empty routine
(  609)    END SUBROUTINE lim_wri
(  610) #endif
(  611) 
(  612)    !!======================================================================
(  613) END MODULE limwri



















PGF90 (Version     10.5)          09/23/2010  09:50:00      page 12

(    1) # 1 "limwri.F90"
(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) 
(    3) !!DB 2008.05.21
(    4) !!Modified to create and output 1 .nc file
(    5) !!Tested vs IOIPSL output = identical
(    6) !!05.27 -- Added averaging of arrays 
(    7) 
(    8) MODULE limwri
(    9)    !!======================================================================
(   10)    !!                     ***  MODULE  limwri  ***
(   11)    !!         Ice diagnostics :  write ice output files
(   12)    !!======================================================================
(   13)    !!----------------------------------------------------------------------
(   14)    !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
(   15)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limwri.F90,v 1.6 2005/12/12 14:18:00 opalod Exp $
(   16)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   17)    !!----------------------------------------------------------------------
(   18) # 19
(   19)    !!----------------------------------------------------------------------
(   20)    !!   'key_ice_lim'                                     LIM sea-ice model
(   21)    !!----------------------------------------------------------------------
(   22)    !!   lim_wri      : write of the diagnostics variables in ouput file 
(   23)    !!   lim_wri_init : initialization and namelist read
(   24)    !!----------------------------------------------------------------------
(   25)    !! * Modules used
(   26)    USE ioipsl
(   27)    USE dianam    ! build name of file (routine)
(   28)    USE phycst
(   29)    USE dom_oce
(   30)    USE daymod
(   31)    USE in_out_manager
(   32)    USE ice_oce         ! ice variables
(   33)    USE flx_oce
(   34)    USE dom_ice
(   35)    USE ice
(   36)    USE iceini
(   37)    USE lbclnk
(   38) 
(   39) !!DB
(   40)    USE lib_ncdf        ! netCDF I/O library
(   41) 
(   42)    IMPLICIT NONE
(   43)    PRIVATE
(   44) 
(   45)    !! * Accessibility
(   46)    PUBLIC lim_wri        ! routine called by lim_step.F90
(   47) 
(   48)    !! * Module variables
(   49)    INTEGER, PARAMETER ::   &  !:
(   50)       jpnoumax = 40             !: maximum number of variable for ice output
(   51)    INTEGER  ::                                &
(   52)       noumef                                     ! number of fields
(   53)    REAL(wp)           , DIMENSION(jpnoumax) ::  &
(   54)       cmulti ,                                &  ! multiplicative constant
(   55)       cadd                                       ! additive constant
(   56)    CHARACTER(len = 35), DIMENSION(jpnoumax) ::  &
(   57)       titn                                       ! title of the field






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 13

(   58)    CHARACTER(len = 8 ), DIMENSION(jpnoumax) ::  &
(   59)       nam                                        ! name of the field
(   60)    CHARACTER(len = 8 ), DIMENSION(jpnoumax) ::  &
(   61)       uni                                        ! unit of the field
(   62)    INTEGER            , DIMENSION(jpnoumax) ::  &
(   63)       nc                                         ! switch for saving field ( = 1 ) or not ( = 0 )
(   64) 
(   65)    REAL(wp)  ::            &  ! constant values
(   66)       epsi16 = 1.e-16   ,  &
(   67)       zzero  = 0.e0     ,  &
(   68)       zone   = 1.e0
(   69)    !!-------------------------------------------------------------------
(   70) 
(   71) CONTAINS
(   72) 
(   73)    SUBROUTINE lim_wri
(   74)       !!-------------------------------------------------------------------
(   75)       !!  This routine computes the average of some variables and write it
(   76)       !!  on the ouput files.
(   77)       !!  ATTENTION cette routine n'est valable que si le pas de temps est
(   78)       !!  egale a une fraction entiere de 1 jours.
(   79)       !!  Diff 1-D 3-D : suppress common also included in etat
(   80)       !!                 suppress cmoymo 11-18
(   81)       !!  modif : 03/06/98
(   82)       !!-------------------------------------------------------------------
(   83) 
(   84) !!DB
(   85)      
(   86) 
(   87)       !! * Local variables
(   88)       REAL(wp),DIMENSION(1) ::   zdept
(   89)       
(   90)       REAL(wp) :: &
(   91)          zsto, zsec, zjulian,zout, &
(   92)          zindh,zinda,zindb,  &
(   93)          ztmu
(   94)       REAL(wp), DIMENSION(jpi,jpj,jpnoumax) :: &
(   95)          zcmo
(   96)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   97)          zfield
(   98)       INTEGER ::  ji, jj, jf   ! dummy loop indices
(   99) !!DB
(  100)       INTEGER :: status, rec_num, kt
(  101)       CHARACTER(len = 80) :: fname
(  102) !!DB -- related to time averaging of ice fields
(  103) !!NB: ice_field_index could be allocated, but instead just give it a v.big dimension
(  104)       INTEGER, SAVE ::  num_ice_fields = 0
(  105)       INTEGER, DIMENSION(100), SAVE :: ice_field_index = 0
(  106)       REAL(wp),DIMENSION(:,:,:),ALLOCATABLE, SAVE :: ave_ice_field
(  107)       LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.    ! Use IOIPSL subroutines for output
(  108) 
(  109)       CHARACTER(len = 40)  :: &
(  110)          clhstnam, clop
(  111) 
(  112)       INTEGER , SAVE ::      &
(  113)          nice, nhorid, ndim, niter, ndepid
(  114)       INTEGER , DIMENSION( jpij ) , SAVE ::  &
(  115)          ndex51  






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 14

(  116)       !!-------------------------------------------------------------------
(  117)       
(  118) 
(  119) !!DB
(  120)       fname = trim(cexper)//'_icemod.nc'
(  121) 
(  122)       IF ( numit == nstart ) THEN 
(  123) !!!DB
(  124)          call ncdf_create_ice_file(status)
(  125) 
(  126)          CALL lim_wri_init                 
(  127) 
(  128) !!DB -- determine # of ice output variables and allocate memory
(  129) !!Eventually only do if clop="ave(x)"
(  130)          do jf = 1, noumef
(  131)             if ( nc(jf) == 1 ) then 
(  132)                num_ice_fields = num_ice_fields + 1
(  133)                ice_field_index(jf) = num_ice_fields
(  134)             endif
(  135)          enddo
(  136)          ALLOCATE(ave_ice_field(jpi,jpj,num_ice_fields))
(  137)          ave_ice_field = 0.0
(  138) 
(  139) 
(  140) 
(  141)          
(  142)          !---5----|----5----|----5----|----5----|----5----|----5----|----5----|72
(  143)          !  1) INITIALIZATIONS.                                                 |
(  144)          !-----------------------------------------------------------------------
(  145)          
(  146)          !-- essai NetCDF
(  147)          
(  148)          zsto= rdt_ice
(  149) !!Chris         clop     = "ave(only(x))"      !ibug  namelist parameter a ajouter
(  150)          clop = "ave(x)"
(  151)          zout = nwrite * rdt_ice / nfice
(  152)          zsec     = 0.
(  153)          niter    = 0
(  154)          zdept(1) = 0.
(  155) 
(  156)          if(USE_IOIPSL) then
(  157)             CALL ymds2ju ( nyear, nmonth, nday, zsec, zjulian )
(  158)             CALL dia_nam ( clhstnam, nwrite, 'icemod' )
(  159)             CALL histbeg ( clhstnam, jpi, glamt, jpj, gphit, 1, jpi, 1, jpj, 0, zjulian, rdt_ice, nhorid, nice , domain_id=nidom)
(  160)             CALL histvert( nice, "deptht", "Vertical T levels", "m", 1, zdept, ndepid)
(  161)             CALL wheneq  ( jpij , tmask(:,:,1), 1, 1., ndex51, ndim)
(  162)             
(  163)             DO jf = 1, noumef
(  164)                IF ( nc(jf) == 1 ) THEN
(  165)                   CALL histdef( nice, nam(jf), titn(jf), uni(jf), jpi, jpj   &
(  166)                        , nhorid, 1, 1, 1, -99, 32, clop, zsto, zout )
(  167)                ENDIF
(  168)             END DO
(  169)             CALL histend(nice)
(  170)          endif
(  171) 
(  172)       ENDIF
(  173)       






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 15

(  174)       !---5----|----5----|----5----|----5----|----5----|----5----|----5----|72
(  175)       !--2. Computation of instantaneous values                                         |
(  176)       !-----------------------------------------------------------------------
(  177) 
(  178)       IF(lwp) THEN
(  179)          WRITE(numout,*)
(  180)          WRITE(numout,*) 'lim_wri : write ice outputs in NetCDF files at time : ', nyear, nmonth, nday, numit
(  181)          WRITE(numout,*) '~~~~~~~ '
(  182)       ENDIF
(  183) 
(  184)       !-- calculs des valeurs instantanees
(  185)       
(  186)       zcmo(:,:, 1:jpnoumax ) = 0.e0 
(  187)       DO jj = 2 , jpjm1
(  188)          DO ji = 2 , jpim1
(  189)             zindh  = MAX( zzero , SIGN( zone , hicif(ji,jj) * (1.0 - frld(ji,jj) ) - 0.10 ) )
(  190)             zinda  = MAX( zzero , SIGN( zone , ( 1.0 - frld(ji,jj) ) - 0.10 ) )
(  191)             zindb  = zindh * zinda
(  192)             ztmu   = MAX( 0.5 * zone , ( tmu(ji,jj) + tmu(ji+1,jj) + tmu(ji,jj+1) + tmu(ji+1,jj+1) ) ) 
(  193)             zcmo(ji,jj,1)  = hsnif (ji,jj)
(  194)             zcmo(ji,jj,2)  = hicif (ji,jj)
(  195)             zcmo(ji,jj,3)  = hicifp(ji,jj)
(  196)             zcmo(ji,jj,4)  = frld  (ji,jj)
(  197)             zcmo(ji,jj,5)  = sist  (ji,jj)
(  198)             zcmo(ji,jj,6)  = fbif  (ji,jj)
(  199)             zcmo(ji,jj,7)  = zindb * (  u_ice(ji,jj  ) * tmu(ji,jj  ) + u_ice(ji+1,jj  ) * tmu(ji+1,jj  )   &
(  200)                                       + u_ice(ji,jj+1) * tmu(ji,jj+1) + u_ice(ji+1,jj+1) * tmu(ji+1,jj+1) ) &
(  201)                                   / ztmu 
(  202) 
(  203)             zcmo(ji,jj,8)  = zindb * (  v_ice(ji,jj  ) * tmu(ji,jj  ) + v_ice(ji+1,jj  ) * tmu(ji+1,jj  )   &
(  204)                                       + v_ice(ji,jj+1) * tmu(ji,jj+1) + v_ice(ji+1,jj+1) * tmu(ji+1,jj+1) ) &
(  205)                                   / ztmu
(  206)             zcmo(ji,jj,9)  = sst_io(ji,jj)
(  207)             zcmo(ji,jj,10) = sss_io(ji,jj)
(  208) 
(  209)             zcmo(ji,jj,11) = fnsolar(ji,jj) + fsolar(ji,jj)
(  210)             zcmo(ji,jj,12) = fsolar (ji,jj)
(  211)             zcmo(ji,jj,13) = fnsolar(ji,jj)
(  212)             ! See thersf for the coefficient
(  213)             zcmo(ji,jj,14) = - fsalt(ji,jj) * rday * ( sss_io(ji,jj) + epsi16 ) / soce
(  214)             zcmo(ji,jj,15) = gtaux(ji,jj)
(  215)             zcmo(ji,jj,16) = gtauy(ji,jj)
(  216)             zcmo(ji,jj,17) = ( 1.0 - frld(ji,jj) ) * qsr_ice (ji,jj) + frld(ji,jj) * qsr_oce (ji,jj)
(  217)             zcmo(ji,jj,18) = ( 1.0 - frld(ji,jj) ) * qnsr_ice(ji,jj) + frld(ji,jj) * qnsr_oce(ji,jj)
(  218)             zcmo(ji,jj,19) = sprecip(ji,jj)
(  219)          END DO
(  220)       END DO
(  221)                 
(  222)       !
(  223)       ! ecriture d'un fichier netcdf
(  224)       !
(  225)       niter = niter + 1
(  226)       DO jf = 1 , noumef
(  227)          DO jj = 1 , jpj
(  228)             DO ji = 1 , jpi
(  229)                zfield(ji,jj) = zcmo(ji,jj,jf) * cmulti(jf) + cadd(jf)
(  230)             END DO
(  231)          END DO






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 16

(  232)          
(  233)          IF ( jf == 7  .OR. jf == 8  .OR. jf == 11 .OR. jf == 12 .OR. jf == 15 .OR.   &
(  234)             jf == 23 .OR. jf == 24 .OR. jf == 16 ) THEN 
(  235)             CALL lbc_lnk( zfield, 'T', -1. )
(  236)          ELSE 
(  237)             CALL lbc_lnk( zfield, 'T',  1. )
(  238)          ENDIF
(  239)          
(  240)          if(USE_IOIPSL) then
(  241)             IF ( nc(jf) == 1 ) CALL histwrite( nice, nam(jf), niter, zfield, ndim, ndex51 )
(  242)          endif
(  243) 
(  244) !!DB output routines; 
(  245)          if(nc(jf) == 1) then  !!accumulate field
(  246)                   ave_ice_field(:,:,ice_field_index(jf)) =  & 
(  247)                        ave_ice_field(:,:,ice_field_index(jf)) + zfield(:,:) 
(  248)          endif
(  249) 
(  250)          if(mod(niter,int(nwrite/nfice)) ==  0) then
(  251) 
(  252)             if(nc(jf) == 1) then
(  253)                kt = nfice*niter + nit000 - 1
(  254)                rec_num = niter/int(nwrite/nfice)
(  255)                zfield(:,:) = ave_ice_field(:,:,ice_field_index(jf))/(float(nwrite)/float(nfice))
(  256)                CALL ncdf_write(fname, 'time_counter', REAL(kt * rdt), rec_num, status)
(  257)                CALL ncdf_write(fname, nam(jf), zfield, -rec_num, status)
(  258)                ave_ice_field(:,:,ice_field_index(jf)) = 0.0
(  259) 
(  260) !!AD/DB 2009.09.30
(  261)                CALL ncdf_write(fname, 'ndastp',REAL(ndastp), rec_num, status)
(  262)                CALL ncdf_write(fname, 'model_time_step',REAL(kt), rec_num, status)
(  263)                CALL ncdf_write(fname, 'model_time',model_time, rec_num, status)
(  264) 
(  265)             endif
(  266) 
(  267)          endif
(  268)          
(  269)       END DO
(  270)       
(  271)       if(USE_IOIPSL) then
(  272)          IF ( ( nfice * niter + nit000 - 1 ) >= nitend ) THEN
(  273)             CALL histclo( nice ) 
(  274)          ENDIF
(  275)       endif
(  276) 
(  277)    END SUBROUTINE lim_wri
(  278) 
(  279)    
(  280)    SUBROUTINE lim_wri_init
(  281)       !!-------------------------------------------------------------------
(  282)       !!                    ***   ROUTINE lim_wri_init  ***
(  283)       !!                
(  284)       !! ** Purpose :   ???
(  285)       !!
(  286)       !! ** Method  : Read the namicewri namelist and check the parameter 
(  287)       !!       values called at the first timestep (nit000)
(  288)       !!
(  289)       !! ** input   :   Namelist namicewri






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 17

(  290)       !!
(  291)       !! history :
(  292)       !!  8.5  ! 03-08 (C. Ethe) original code
(  293)       !!-------------------------------------------------------------------
(  294)       !! * Local declarations
(  295)       INTEGER ::   nf      ! ???
(  296) 
(  297)       TYPE FIELD 
(  298)          CHARACTER(len = 35) :: ztitle 
(  299)          CHARACTER(len = 8 ) :: zname          
(  300)          CHARACTER(len = 8 ) :: zunit
(  301)          INTEGER             :: znc   
(  302)          REAL                :: zcmulti 
(  303)          REAL                :: zcadd        
(  304)       END TYPE FIELD
(  305) 
(  306)       TYPE(FIELD) ::  &
(  307)          field_1 , field_2 , field_3 , field_4 , field_5 , field_6 ,   &
(  308)          field_7 , field_8 , field_9 , field_10, field_11, field_12,   &
(  309)          field_13, field_14, field_15, field_16, field_17, field_18,   &
(  310)          field_19
(  311) 
(  312)       TYPE(FIELD) , DIMENSION(jpnoumax) :: zfield
(  313) 
(  314)       NAMELIST/namiceout/ noumef, &
(  315)          field_1 , field_2 , field_3 , field_4 , field_5 , field_6 ,   &
(  316)          field_7 , field_8 , field_9 , field_10, field_11, field_12,   &
(  317)          field_13, field_14, field_15, field_16, field_17, field_18,   &
(  318)          field_19
(  319)       !!-------------------------------------------------------------------
(  320) 
(  321) 
(  322)       ! Read Namelist namicewri
(  323)       REWIND ( numnam_ice )
(  324)       READ   ( numnam_ice  , namiceout )
(  325)       zfield(1)  = field_1
(  326)       zfield(2)  = field_2
(  327)       zfield(3)  = field_3
(  328)       zfield(4)  = field_4
(  329)       zfield(5)  = field_5
(  330)       zfield(6)  = field_6
(  331)       zfield(7)  = field_7
(  332)       zfield(8)  = field_8
(  333)       zfield(9)  = field_9
(  334)       zfield(10) = field_10
(  335)       zfield(11) = field_11
(  336)       zfield(12) = field_12
(  337)       zfield(13) = field_13
(  338)       zfield(14) = field_14
(  339)       zfield(15) = field_15
(  340)       zfield(16) = field_16
(  341)       zfield(17) = field_17
(  342)       zfield(18) = field_18
(  343)       zfield(19) = field_19
(  344)       
(  345)       DO nf = 1, noumef
(  346)          titn  (nf) = zfield(nf)%ztitle
(  347)          nam   (nf) = zfield(nf)%zname






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 18

(  348)          uni   (nf) = zfield(nf)%zunit
(  349)          nc    (nf) = zfield(nf)%znc
(  350)          cmulti(nf) = zfield(nf)%zcmulti
(  351)          cadd  (nf) = zfield(nf)%zcadd
(  352)       END DO
(  353) 
(  354)       IF(lwp) THEN
(  355)          WRITE(numout,*)
(  356)          WRITE(numout,*) 'lim_wri_init : Ice parameters for outputs'
(  357)          WRITE(numout,*) '~~~~~~~~~~~~'
(  358)          WRITE(numout,*) '    number of fields to be stored         noumef = ', noumef
(  359)          WRITE(numout,*) '           title                            name     unit      Saving (1/0) ',   &
(  360)             &            '    multiplicative constant       additive constant '
(  361)          DO nf = 1 , noumef         
(  362)             WRITE(numout,*) '   ', titn(nf), '   ', nam(nf),'      ', uni(nf),'  ', nc(nf),'        ', cmulti(nf),   &
(  363)                '        ', cadd(nf)
(  364)          END DO
(  365)       ENDIF
(  366)             
(  367)    END SUBROUTINE lim_wri_init
(  368) 
(  369) 
(  370) !!DB
(  371) !!Create 1 ice file equivalent to the one created in lim_wri
(  372) 
(  373)  SUBROUTINE ncdf_create_ice_file(status)
(  374)     IMPLICIT NONE
(  375)     ! Subroutine argument declarations
(  376) 
(  377)     INTEGER,INTENT(OUT) :: status
(  378) 
(  379)     ! Local declarations
(  380)     INTEGER :: ncid,    &  ! netCDF file ID
(  381)                varid,   &  ! ID of netCDF variable to be written to
(  382)                nfstat,  &  ! netCDF library call return status
(  383)                mpistat     ! MPI library call return status
(  384)     INTEGER,DIMENSION(1:4) :: dimids
(  385) !!DB -- unsure how to dimension varids, so choose a large number 
(  386)     INTEGER,DIMENSION(1:50) :: varids
(  387)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  388)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  389)     CHARACTER(LEN=100) :: sec_since    
(  390)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  391)     CHARACTER (len=80) :: fname
(  392)     CHARACTER(LEN=20) :: op_type, varname
(  393)     INTEGER :: int_opp, &              ! Operation interval
(  394)                int_wri                 ! Write interval
(  395)     CHARACTER(LEN=3),PARAMETER :: &
(  396)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  397)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  398)     INTEGER :: jf, varnum
(  399) 
(  400) ! Get required info from namelist_ice
(  401)     if ( numit == nstart ) THEN 
(  402)        CALL lim_wri_init 
(  403)     endif
(  404)     
(  405)     






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 19

(  406)     ! Initializations
(  407)     op_type = "ave(x)"       !!default, dangerous, requires routine to do this  
(  408)     fname = trim(cexper)//'_icemod.nc'
(  409) 
(  410)     status = NCDF_NOERR
(  411)     CALL ioget_calendar(cal_type)
(  412)     CALL ioget_timestamp(timestamp)
(  413)     WRITE (UNIT=sec_since, &
(  414)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  415)          &  nyear,nmonth,nday,0, 0, 0
(  416)     WRITE(t_origin, &
(  417)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  418)          &   nyear,months(nmonth),nday,0,0,0
(  419) 
(  420) !!DB see above
(  421)     int_opp = rdt_ice
(  422)     int_wri = nwrite * rdt_ice/nfice
(  423) 
(  424)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  425)        WRITE(100,*) 'NCDF DEBUG: Creating output file:', fname
(  426)        CALL FLUSH
(  427)     END IF
(  428) 
(  429)     ! Only processor 0 does anything
(  430)     IF(nproc == 0) THEN
(  431)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  432)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Creating file:', fname
(  433)           CALL FLUSH
(  434)        END IF
(  435)        ! Create the file
(  436)        nfstat = nf90_create(fname, nf90_clobber, ncid)
(  437)        IF(nfstat /= nf90_noerr) THEN
(  438)           status = NCDF_NFERR
(  439)           RETURN
(  440)        END IF
(  441)        
(  442)        ! Define dimensions
(  443)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  444)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Defining dimensions in file:', fname
(  445)           CALL FLUSH
(  446)        END IF
(  447)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  448)        nfstat = nf90_def_dim(ncid, 'z', 1, dimids(2))          !!ice model
(  449)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(3))
(  450)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(4))
(  451)        
(  452)        ! Define variables
(  453)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  454)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Defining variables in file:', fname
(  455)           CALL FLUSH
(  456)        END IF
(  457)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  458)             (/ dimids(4), dimids(3) /), &
(  459)             varids(1))
(  460)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  461)             (/ dimids(4), dimids(3) /), &
(  462)             varids(2))
(  463)        nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 20

(  464)             (/ dimids(2) /), &
(  465)             varids(3))
(  466)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  467)             (/ dimids(1) /), &
(  468)             varids(4))
(  469)        
(  470)        ! Add attributes
(  471)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  472)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Writing attributes in file:', fname
(  473)           CALL FLUSH
(  474)        END IF
(  475)        ! nav_lon
(  476)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  477) !       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -7.138476E01)
(  478) !       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', -6.299908E01)
(  479)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', minval(glamt))
(  480)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', maxval(glamt))
(  481)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  482)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  483) 
(  484)        ! nav_lat
(  485)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  486) !       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', 3.852767E01)
(  487) !       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 4.223735e+01)
(  488)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', minval(gphit))
(  489)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', maxval(gphit))
(  490)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  491)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  492) 
(  493)        ! deptht
(  494)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  495)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  496)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 0.0)      !for ice model
(  497)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 0.0)      !for ice model
(  498)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
(  499)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical T levels')
(  500) 
(  501)        ! time_counter
(  502)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  503)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  504)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  505)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
(  506)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
(  507) 
(  508)        ! global
(  509)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  510)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(fname))
(  511)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'production', 'An IPSL model')
(  512)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  513)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  514)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'MISC', 'DB-created file')
(  515) 
(  516) !!Loop over variables in namelist_ice as done above
(  517)        varnum = 4      !start at # of vars defined above
(  518)        do jf = 1, noumef
(  519)           if ( nc(jf) == 1 ) THEN
(  520) 
(  521)              varnum = varnum + 1






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 21

(  522)              nfstat = nf90_def_var(ncid, nam(jf), nf90_float, &
(  523)                   (/ dimids(4), dimids(3), dimids(1) /), varids(varnum))
(  524) 
(  525)              nfstat = nf90_put_att(ncid, varids(varnum), 'units', uni(jf))
(  526)              nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
(  527)              nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
(  528)              nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
(  529)              nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', TRIM(op_type))
(  530)              nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', titn(jf))
(  531)              nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', nam(jf))
(  532)              nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TYX')
(  533)              nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
(  534)              nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
(  535)              nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')
(  536) 
(  537)           endif
(  538)        enddo
(  539) 
(  540) 
(  541) 
(  542) !AD/DB: add new time-related variables
(  543)        varnum = varnum + 1
(  544)        ! ndate (ndastp)
(  545)        nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
(  546)             (/ dimids(1) /), &
(  547)             varids(varnum))
(  548)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
(  549)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'time step date in year/month/day aammjj')
(  550) 
(  551)        varnum = varnum + 1
(  552)        ! ndate (model_time)
(  553)        nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
(  554)             (/ dimids(1) /), &
(  555)             varids(varnum))
(  556)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', &
(  557)             'time step date (when output is writen) in year/month/day aammjj (decimal day)')
(  558)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
(  559)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula1', 'nyear  =   model_time / 10000')       
(  560)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula2', & 
(  561)             'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
(  562)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula3', & 
(  563)             'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           
(  564) 
(  565)        varnum = varnum + 1
(  566)        ! kt 
(  567)        nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
(  568)             (/ dimids(1) /), &
(  569)             varids(varnum))
(  570) 
(  571) 
(  572)        ! Close file
(  573)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  574)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_ice_file - Closing file:', fname
(  575)           CALL FLUSH
(  576)        END IF
(  577)        nfstat = nf90_close(ncid)
(  578)        IF(nfstat /= nf90_noerr) THEN
(  579)           status = NCDF_NFERR






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 22

(  580)           RETURN
(  581)        END IF
(  582)     END IF
(  583) 
(  584) 
(  585) 
(  586) ! On first write step only
(  587)     CALL ncdf_write(fname, 'nav_lat', gphit, -1, status)
(  588)     CALL ncdf_write(fname, 'nav_lon', glamt, -1, status)
(  589) 
(  590) 
(  591)     ! Sync up processors before returning from subroutine
(  592)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
(  593)     IF(mpistat /= 0) THEN
(  594)        status = NCDF_MPERR
(  595)        RETURN
(  596)     END IF
(  597)     
(  598)   END SUBROUTINE ncdf_create_ice_file
(  599) 
(  600) 
(  601) 
(  602) 
(  603) 
(  604) # 612
(  612)    !!======================================================================
(  613) END MODULE limwri
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (limwri.F90: 255)
