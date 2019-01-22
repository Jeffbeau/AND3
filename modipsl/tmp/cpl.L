


PGF90 (Version     10.5)          09/23/2010  09:49:27      page 1

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

Filename: cpl.F90

(    1) MODULE cpl
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE cpl  ***
(    4)    !! Coupled O/A : coupled ocean-atmosphere case using OASIS
(    5)    !!=====================================================================
(    6) #if defined key_coupled
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_coupled'                              coupled Ocean/Atmosphere
(    9)    !!----------------------------------------------------------------------
(   10)    !!   cpl_init     : initialization of coupled mode communication
(   11)    !!   cpl_read     : read the coupled namelist
(   12)    !!   cpl_stp      : exchange fields in coupled mode
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE cpl_oce         ! coupled exchange variables (???)
(   19)    USE flx_oce         ! in case of ice 
(   20)    USE ocfzpt          ! ocean freezing point
(   21)    USE daymod          ! calendar
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! Routine accessibility
(   27)    PUBLIC cpl_init     ! routine called in opa module
(   28)    PUBLIC cpl_stp      ! routine called in step module
(   29)    !!----------------------------------------------------------------------
(   30)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/cpl.F90,v 1.6 2005/12/12 14:18:01 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35) CONTAINS
(   36) 
(   37)    SUBROUTINE cpl_init( kastp, kexch, kstep )
(   38)       !!---------------------------------------------------------------------
(   39)       !!                  ***  ROUTINE cpl_init  ***
(   40)       !!
(   41)       !! ** Purpose :   Initialize coupled mode communication for ocean
(   42)       !!    exchange process identifiers and timestep information
(   43)       !!    between AGCM, OGCM and COUPLER. (OASIS software)
(   44)       !!
(   45)       !! ** Method  :  3 types :
(   46)       !!      a) Use named pipes(FIFO) to exchange process identifiers






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 2

(   47)       !!          between the 3 programs
(   48)       !!      b) USE a messag passing method with PVM language (CLIM)
(   49)       !!      c) USE SVIPC method
(   50)       !!
(   51)       !! ** Input   :   npiat     : agcm process identifier
(   52)       !!                npicp     : coupler process identifier
(   53)       !!                npioc     : ogcm process identifier
(   54)       !!
(   55)       !! Reference  :   see Epicoa 0803 (1992)
(   56)       !!
(   57)       !! History :
(   58)       !!        !  92-09  (L. Terray)  Original code
(   59)       !!        !  96-07  (L. Terray)  OASIS version 2
(   60)       !!        !  96-11  (E. Guilyardi)  run-off + Taux,Tauy x2
(   61)       !!        !  98-04  (M.A Foujols, S. Valcke, M. Imbard)  OASIS2.2
(   62)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   63)       !!----------------------------------------------------------------------
(   64)       !! * Arguments
(   65)       INTEGER, INTENT(in ) ::   &
(   66)          kastp,      &  ! total number of timesteps in oceanic model
(   67)          kexch,      &  ! frequency of exchange for the fluxes (in time steps)
(   68)          kstep          ! timestep value (in seconds)
(   69) 
(   70)       !! * Local declarations
(   71)       INTEGER,DIMENSION(3)  :: iparal
(   72)       INTEGER               :: ifcpl, idt, info, imxtag, istep
(   73) 
(   74)       CHARACTER (len=9) ::   clpoolnam 
(   75)       INTEGER           :: ipoolhandle, ipoolsize, jf
(   76)       CHARACTER (len=3) ::   cljobnam      ! experiment name
(   77)       INTEGER           :: ierror
(   78) !      INTEGER,DIMENSION(4) ::  imess
(   79)       INTEGER,DIMENSION(4) ::  imesso
(   80)       !!----------------------------------------------------------------------
(   81)      
(   82)       IF(lwp) WRITE(numout,*)
(   83)       IF(lwp) WRITE(numout,*) 'cpl_init : initialization in coupled ocean/atmosphere case'
(   84)       IF(lwp) WRITE(numout,*) '~~~~~~~~'
(   85)       IF(lwp) WRITE(numout,*)
(   86)      
(   87) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily || defined key_flx_forced_daily
(   88)       IF(lwp)WRITE(numout,cform_err)
(   89)       IF(lwp)WRITE(numout,*) ' key_coupled and key_flx_bulk_... are incompatible'
(   90)       nstop = nstop + 1
(   91) #endif
(   92)  
(   93)       IF(lwp)WRITE(numout,*)'     coupled simulation'
(   94)       IF(lwp)WRITE(numout,*)'        unit ',numlhf,' receive atm fluxes'
(   95)       IF(lwp)WRITE(numout,*)'        unit ',numlws,' receive windstress'
(   96)       IF(lwp)WRITE(numout,*)'        unit ',numlts,' transfer sst'
(   97)       IF(lwp)WRITE(numout,*)'        unit ',numlic,' transfer ice cover'
(   98) 
(   99) 
(  100)       CALL cpl_read           ! read the coupled mode information in namelist
(  101) 
(  102)       CALL flush(numout)
(  103) 
(  104)       ! I- PIPE






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 3

(  105)       ! --------
(  106)       ! W A R N I N G : PIPE technique is temporary diseable (nov. 2000)
(  107) 
(  108)       IF( cchan == 'PIPE' ) THEN 
(  109) 
(  110)          IF(lwp) WRITE(numout,*)
(  111)          IF(lwp) WRITE(numout,*) 'Making pipes for fields to receive from CPL'
(  112)          IF(lwp) WRITE(numout,*)
(  113)          CALL flush(numout)
(  114)          ierror = 0
(  115) 
(  116)          ! loop to define pipes (CPL=atmos to ocean)
(  117) 
(  118)          DO jf = 1, nflxc2o  
(  119)             ! CALL PIPE_Model_Define( numout, cpl_readflx(jf), jpread, info )
(  120)             IF( info /= 0 ) ierror = ierror + 1
(  121)          END DO
(  122)          DO jf = 1, ntauc2o
(  123)             ! CALL PIPE_Model_Define( numout, cpl_readtau(jf), jpread, info )
(  124)             IF( info /= 0 ) ierror = ierror + 1
(  125)          END DO
(  126)          
(  127)          IF(lwp) WRITE(numout,*) ' '
(  128)          IF(lwp) WRITE(numout,*) 'Making pipes for fields to send to CPL'
(  129)          IF(lwp) WRITE(numout,*) ' '
(  130)          
(  131)          ! loop to define pipes (ocean to atmos=CPL)
(  132) 
(  133)          DO jf = 1, nfldo2c
(  134)             ! CALL PIPE_Model_Define( numout, cpl_writ(jf), jpwrit, info )
(  135)             IF( info /= 0 ) ierror = ierror + 1
(  136)          END DO
(  137) 
(  138)          IF( ierror /= 0 ) THEN
(  139)             IF(lwp) WRITE(numout,*)
(  140)             IF(lwp) WRITE(numout,*) 'cpl_init: end of job due to error in pipes definitions'
(  141)             CALL abort('')
(  142)          END IF
(  143)          
(  144)          IF(lwp) WRITE(numout,*)
(  145)          IF(lwp) WRITE(numout,*) 'All pipes have been made'
(  146)          
(  147)          IF(lwp) WRITE(numout,*)
(  148)          IF(lwp) WRITE(numout,*) 'Communication test between OCE and CPL'
(  149)          CALL flush(numout)
(  150)          
(  151)          ! CALL PIPE_Model_Stepi(numout, imess, nmodcpl, imesso, ierror)
(  152)          
(  153)          IF( ierror /= 0 ) THEN
(  154)             IF(lwp) WRITE(numout,*)
(  155)             IF(lwp) WRITE(numout,*) 'cpl_init: end of job due to error in exchange first informations with Oasis'
(  156)             CALL abort('')
(  157)          END IF
(  158)          
(  159)          IF(lwp) WRITE(numout,*)
(  160)          IF(lwp) WRITE(numout,*) 'Communication test between OCE and CPL is OK'
(  161)          IF(lwp) WRITE(numout,*) ' total simulation time in oasis = ',imesso(1)
(  162)          IF(lwp) WRITE(numout,*) ' total number of iterations is  = ',imesso(2)






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 4

(  163)          IF(lwp) WRITE(numout,*) ' value of oasis timestep  is    = ',imesso(3)
(  164)          IF(lwp) WRITE(numout,*) ' process id for oasis  is       = ',imesso(4)
(  165)          CALL flush(numout)
(  166)          
(  167)          ! II SVIPC
(  168)          ! ---------
(  169)          ! W A R N I N G : SVIPC technique is temporary diseable (nov. 2000)
(  170)          
(  171)          
(  172)       ELSE IF( cchan == 'SIPC' ) THEN
(  173) 
(  174)          ! debug for more information
(  175) 
(  176)          ! CALL SVIPC_debug(1)
(  177) 
(  178)          ! Define the experiment name :
(  179) 
(  180)           cljobnam = 'IPC'      ! as $JOBNAM in namcouple
(  181) 
(  182)           ! Attach to shared memory pool used to exchange initial infos 
(  183) 
(  184)           info = 0
(  185)           ! CALL SIPC_Init_Model (cljobnam, cplmodnam, 1, info)
(  186)           IF( info /= 0 ) THEN 
(  187)              IF(lwp) WRITE(numout,*)
(  188)              IF(lwp) WRITE(numout,*)'WARNING: Problem with attachement to',info 
(  189)              IF(lwp) WRITE(numout,*)'         initial memory pool(s) in ocean'
(  190)              IF(lwp) WRITE(numout,*)
(  191)              CALL abort('STOP in ocean')
(  192)           ENDIF
(  193) 
(  194)           ! Attach to pools used to exchange fields from ocean to coupler
(  195) 
(  196)           DO jf = 1, nfldo2c
(  197)              ! Pool name:
(  198)              clpoolnam = 'P'//cpl_writ(jf)
(  199)              ! CALL SIPC_Attach(clpoolnam, ipoolhandle)
(  200)              ! Resulting pool handle:
(  201)              mpoolwrit(jf) = ipoolhandle  
(  202)           END DO
(  203) 
(  204)           ! Attach to pools used to exchange fields from coupler to ocean
(  205)           
(  206)           DO jf = 1, nflxc2o
(  207)              ! Pool name:
(  208)              clpoolnam = 'P'//cpl_readflx(jf)
(  209)              ! CALL SIPC_Attach(clpoolnam, ipoolhandle)
(  210)              ! Resulting pool handle:
(  211)              mpoolread(jf) = ipoolhandle  
(  212)           END DO  
(  213) 
(  214)           DO jf = 1, ntauc2o
(  215)              ! Pool name:
(  216)              clpoolnam = 'P'//cpl_readtau(jf)
(  217)              ! CALL SIPC_Attach(clpoolnam, ipoolhandle)
(  218)              ! Resulting pool handle:
(  219)              mpoolread(jf+nflxc2o) = ipoolhandle  
(  220)           END DO 






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 5

(  221) 
(  222)           ! Exchange of initial infos
(  223) 
(  224)           ! Write data array isend to pool READ by Oasis
(  225) 
(  226)           info = 0
(  227)           ipoolsize = 4*jpbyteint
(  228)           ! CALL SVIPC_write(mpoolinitr, imess, ipoolsize, info)
(  229) 
(  230)           ! Find error if any
(  231) 
(  232)           IF( info < 0 ) THEN
(  233)              IF(lwp) WRITE(numout,*)
(  234)              IF(lwp) WRITE(numout,*) 'Problem in ocean in writing initial' 
(  235)              IF(lwp) WRITE(numout,*) 'infos to the shared memory segment(s)'
(  236)              IF(lwp) WRITE(numout,*)
(  237)           ELSE
(  238)              IF(lwp) WRITE(numout,*)
(  239)              IF(lwp) WRITE(numout,*) 'Initial infos written in ocean'            
(  240)              IF(lwp) WRITE(numout,*) 'to the shared memory segment(s)'
(  241)              IF(lwp) WRITE(numout,*)
(  242)           ENDIF
(  243) 
(  244)           ! Read data array irecv from pool written by Oasis
(  245) 
(  246)           info = 0
(  247)           ipoolsize = 4 * jpbyteint
(  248)           CALL svipc_read(mpoolinitw, imesso, ipoolsize, info)
(  249) 
(  250)           ! Find error if any
(  251) 
(  252)           IF( info < 0 ) THEN
(  253)              IF(lwp) WRITE(numout,*) '   '
(  254)              IF(lwp) WRITE(numout,*) 'Problem in ocean in reading initial' 
(  255)              IF(lwp) WRITE(numout,*) 'infos from the shared memory segment(s)'
(  256)              IF(lwp) WRITE(numout,*) '   '
(  257)           ELSE
(  258)              IF(lwp) WRITE(numout,*) '   '
(  259)              IF(lwp) WRITE(numout,*) 'Initial infos read by ocean'                
(  260)              IF(lwp) WRITE(numout,*) 'from the shared memory segment(s)'
(  261)              IF(lwp) WRITE(numout,*) '   '
(  262)              IF(lwp) WRITE(numout,*) ' ntime, niter, nstep, Oasis pid:'
(  263)              IF(lwp) WRITE(numout,*) imesso(1), imesso(2), imesso(3), imesso(4) 
(  264)           ENDIF
(  265) 
(  266)           ! Detach from shared memory segment(s)
(  267) 
(  268)           info = 0
(  269)           ! CALL SVIPC_close(mpoolinitw, 0, info)
(  270)           
(  271)           ! Find error if any
(  272) 
(  273)           IF( info < 0 ) THEN
(  274)              IF(lwp) WRITE(numout,*) 'Problem in detaching from shared memory segment(s)'
(  275)              IF(lwp) WRITE(numout,*) 'used by ocean to read initial infos' 
(  276)           ENDIF
(  277) 
(  278)           ! III CLIM






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 6

(  279)           ! --------
(  280) 
(  281)       ELSE IF( cchan == 'CLIM' ) THEN
(  282) 
(  283)          ! Define the number of processors involved in the coupling for
(  284)          ! Oasis (=1) and each model (as last two INTEGER on $CHATYPE line
(  285)          ! in the namcouple); they will be stored in a COMMON in mpiclim.h
(  286)          ! (used for CLIM/MPI2 only)
(  287)          mpi_nproc(0)=1
(  288)          mpi_nproc(1)=1
(  289)          mpi_nproc(2)=1 
(  290) 
(  291)          ! Define the number of processors used by each model as in
(  292)          ! $CHATYPE line of namcouple (used for CLIM/MPI2 only)
(  293)          mpi_totproc(1)=1
(  294)          mpi_totproc(2)=1
(  295)          
(  296)          ! Define names of each model as in $NBMODEL line of namcouple
(  297)          ! (used for CLIM/MPI2 only)
(  298)          cmpi_modnam(1)='lmdz.x'
(  299)          cmpi_modnam(2)=cplmodnam
(  300)          !  
(  301)          ! 1.1-Define the experiment name :
(  302)          
(  303)          cljobnam = 'CLI'      ! as $JOBNAM in namcouple
(  304)          
(  305)          OPEN ( UNIT = 7, FILE = 'trace', STATUS = 'unknown', FORM = 'formatted')
(  306)          CALL clim_init ( cljobnam, cplmodnam, 3, 7,   &
(  307)                           kastp, kexch, kstep,   &
(  308)                           5, 3600, 3600, info )
(  309) 
(  310)          IF( info /= clim_ok ) THEN
(  311)             IF(lwp) WRITE( numout, *) 'cpl_init : pb init clim, error code is = ', info
(  312)             CALL abort( 'STOP in cpl_init' )
(  313)          ELSE
(  314)             IF(lwp) WRITE(numout,*) 'cpl_init : init clim ok '
(  315)          ENDIF
(  316)          
(  317)          iparal ( clim_strategy ) = clim_serial
(  318)          iparal ( clim_length   ) = jpiglo*jpjglo
(  319)          iparal ( clim_offset   ) = 0
(  320)          
(  321)          ! loop to define messages (CPL=atmos to ocean)
(  322)          DO jf = 1, nflxc2o
(  323)             CALL CLIM_Define ( cpl_readflx(jf), clim_in, clim_double, iparal, info )  
(  324)          END DO
(  325)          DO jf = 1, ntauc2o
(  326)             CALL CLIM_Define ( cpl_readtau(jf), clim_in, clim_double, iparal, info )  
(  327)          END DO
(  328)          
(  329)          ! loop to define messages (ocean to CPL=atmos)
(  330)          DO jf = 1, nfldo2c
(  331)             CALL CLIM_Define ( cpl_writ(jf), clim_out, clim_double, iparal, info )   
(  332)          END DO
(  333)          
(  334)          IF(lwp) WRITE(numout,*) 'cpl_init : clim_define ok '
(  335)          
(  336)          CALL CLIM_Start( imxtag, info )






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 7

(  337)          
(  338)          IF( info /= clim_ok ) THEN
(  339)             IF(lwp) WRITE(numout,*) 'cpl_init : pb start clim, error code is = ', info
(  340)             CALL abort( 'stop in cpl_init' )
(  341)          ELSE
(  342)             IF(lwp) WRITE(numout,*) 'cpl_init : start clim ok '
(  343)          ENDIF
(  344)          
(  345)          CALL CLIM_Stepi ( cploasis, istep, ifcpl, idt, info )
(  346) 
(  347)          IF( info /= clim_ok ) THEN
(  348)             IF(lwp) WRITE(numout,*) ' warning : problem in getting step info from oasis '
(  349)             IF(lwp) WRITE(numout,*) ' =======   error code number = ', info
(  350)          ELSE
(  351)             IF(lwp) WRITE(numout,*) ' got step information from oasis '
(  352)          ENDIF
(  353)          IF(lwp) WRITE(numout,*) ' number of tstep in oasis ', istep
(  354)          IF(lwp) WRITE(numout,*) ' exchange frequency in oasis ', ifcpl
(  355)          IF(lwp) WRITE(numout,*) ' length of tstep in oasis ', idt
(  356)       ENDIF
(  357) 
(  358)    END SUBROUTINE cpl_init
(  359) 
(  360) 
(  361)    SUBROUTINE cpl_read
(  362)       !!---------------------------------------------------------------------
(  363)       !!                  ***  ROUTINE cpl_read  ***
(  364)       !!                    
(  365)       !! ** Purpose :   Read and print options for the coupled run (namelist)
(  366)       !!
(  367)       !! ** Method  :   ???
(  368)       !!
(  369)       !! History :
(  370)       !!   8.5  !  02-12  (G. Madec)  F90: Free form and module
(  371)       !!----------------------------------------------------------------------
(  372)       !! * Local declarations
(  373)       INTEGER :: jf
(  374) 
(  375)       NAMELIST/namcpl/ nexco, cchan, nmodcpl, cplmodnam, cploasis   &
(  376)           , nfldo2c, nflxc2o, ntauc2o, cpl_f_readflx, cpl_f_readtau   &
(  377)           , cpl_f_writ, cpl_readflx, cpl_readtau, cpl_writ
(  378)       !!----------------------------------------------------------------------
(  379)       
(  380)       IF(lwp) WRITE(numout,*)
(  381)       IF(lwp) WRITE(numout,*) ' cpl_read : read the coupled parameters in namelist'
(  382)       IF(lwp) WRITE(numout,*) ' ~~~~~~~~'
(  383)       IF(lwp) WRITE(numout,*)
(  384) 
(  385)       ! Default values
(  386)       
(  387)       nexco = 24
(  388)       cchan='PIPE'              ! echange TYPE
(  389)       nmodcpl = 2
(  390)       cplmodnam = 'opa.xx'      ! model name : as $NBMODEL in namcouple
(  391)       cploasis = 'Oasis'        ! coupler name : as in coupler
(  392) 
(  393)       ! -Define symbolic name for fields exchanged from ocean to coupler,
(  394)       ! must be the same as (1) of the field  definition in namcouple:






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 8

(  395)       nfldo2c=2
(  396)       cpl_writ(1)='SOSSTSST'
(  397)       cpl_writ(2)='SOICECOV'
(  398) 
(  399)       ! -Define files name for fields exchanged from ocean to coupler,
(  400)       ! must be the same as (6) of the field  definition in namcouple:
(  401)       nflxc2o=6
(  402)       cpl_readflx(1)='SONSFLDO' ! non solar heat flux (positif vers l'ocean)
(  403)       cpl_readflx(2)='SOSHFLDO' ! solar flux
(  404)       cpl_readflx(3)='SOTOPRSU' ! precip
(  405)       cpl_readflx(4)='SOTFSHSU' ! evaporation
(  406)       cpl_readflx(5)='SORUNCOA' ! coastal run-off
(  407)       cpl_readflx(6)='SORIVFLU' ! river run-off
(  408)       ntauc2o=2
(  409)       cpl_readtau(1)='SOZOTAUX' ! first zonal wind stress
(  410)       cpl_readtau(2)='SOZOTAU2' ! second zonal wind stress
(  411)       cpl_readtau(3)='SOMETAUY' ! first meridien wind stress
(  412)       cpl_readtau(4)='SOMETAU2' ! second meridien wind stress
(  413) 
(  414)       ! -Define files name for fields exchanged from ocean to coupler,
(  415)       ! must be the same as (6) of the field  definition in namcouple:
(  416)       cpl_f_writ(1)='ocesst'
(  417)       cpl_f_writ(2)='oceice'
(  418) 
(  419)       ! -Define files name for fields exchanged from coupler to ocean
(  420)       ! must be the same as (7) of the field  definition in namcouple:
(  421)       cpl_f_readflx(1)='oceflx'
(  422)       cpl_f_readflx(2)='oceflx'
(  423)       cpl_f_readflx(3)='oceflx'
(  424)       cpl_f_readflx(4)='oceflx'
(  425)       cpl_f_readflx(5)='oceflx'
(  426)       cpl_f_readflx(6)='oceflx'
(  427)       cpl_f_readtau(1)='ocetau'
(  428)       cpl_f_readtau(2)='ocetau'
(  429)       cpl_f_readtau(3)='ocetau'
(  430)       cpl_f_readtau(4)='ocetau'
(  431) 
(  432)       ! Namelist namcpl : coupling mode and frequency
(  433)       REWIND( numnam )
(  434)       READ  ( numnam, namcpl )
(  435) 
(  436)       IF(lwp) THEN
(  437)          WRITE(numout,*) 'namcpl'
(  438)          WRITE(numout,*) 
(  439)          WRITE(numout,*) ' Coupling exchange frequency    nexco  = ',nexco
(  440)          WRITE(numout,*) ' Coupling exchange technique    cchan  = ',cchan
(  441)          WRITE(numout,*) ' Mode Coupling technique      nmodcpl  = ',nmodcpl
(  442)          WRITE(numout,*) ' Define the model name      cplmodnam  = ',cplmodnam
(  443)          WRITE(numout,*) ' Define the coupler name      cploasis = ',cploasis
(  444)          WRITE(numout,*) ' Fields number ocean to coupler nfldo2c= ',nfldo2c
(  445)          WRITE(numout,*) ' Flux fields coupler to ocean nflxc2o  = ',nflxc2o
(  446)          WRITE(numout,*) ' Stress fields coupler to ocean ntauc2o= ',ntauc2o
(  447)          IF ( cchan == 'PIPE' .OR.  cchan == 'pipe' ) THEN
(  448)             cchan='PIPE'
(  449)             WRITE(numout,*)'(communication between models made by pipes)'
(  450)          ELSEIF( cchan == 'CLIM' .OR. cchan == 'clim' ) THEN
(  451)             cchan='CLIM'
(  452)             WRITE(numout,*)'(communication between models made by CLIM-PVM)'






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 9

(  453)          ELSEIF( cchan == 'SIPC' .OR. cchan == 'sipc' ) THEN
(  454)             cchan='SIPC'
(  455)             WRITE(numout,*)'(communication between models made by the',  &
(  456)                'Share Memory Segment and Semaphore)'
(  457)          ELSE
(  458)             WRITE(numout,*) 'technic not yet implemented ',cchan
(  459)             STOP 'in cpl_read'
(  460)          ENDIF
(  461)          DO jf = 1, nflxc2o
(  462)             WRITE(numout,*) 'file to receive field number = ',jf,'  ',cpl_f_readflx(jf) 
(  463)          END DO
(  464)          DO jf = 1, ntauc2o
(  465)             WRITE(numout,*) 'file to receive field number = ',jf,'  ',cpl_f_readtau(jf) 
(  466)          END DO
(  467)          DO jf = 1, nfldo2c
(  468)             WRITE(numout,*) 'file to send field number = ',jf,'  ',cpl_f_writ(jf)
(  469)          END DO
(  470)          WRITE(numout,*) ' fields received from coupler'
(  471)          DO jf = 1, nflxc2o
(  472)             WRITE(numout,*) 'symbolic name for field number = ',jf,'  ',cpl_readflx(jf) 
(  473)          END DO
(  474)          DO jf = 1, ntauc2o
(  475)             WRITE(numout,*) 'symbolic name for field number = ',jf,'  ',cpl_readtau(jf) 
(  476)          END DO
(  477)          WRITE(numout,*) ' fields send to coupler'
(  478)          DO jf = 1, nfldo2c
(  479)             WRITE(numout,*) 'symbolic name for field number = ',jf,'  ',cpl_writ(jf)
(  480)          END DO
(  481)       ENDIF
(  482) 
(  483)    END SUBROUTINE cpl_read
(  484) 
(  485) 
(  486)    SUBROUTINE cpl_stp( kt )
(  487)       !!---------------------------------------------------------------------
(  488)       !!                  ***  ROUTINE cpl_stp  ***
(  489)       !!                      *****************
(  490)       !!                      * OASIS routine *
(  491)       !!                      *****************
(  492)       !! ** Purpose : - At each coupling time-step,this routine sends fields
(  493)       !!      like sst or ice cover to the coupler.
(  494)       !!              - At each time-step computes time average values 
(  495)       !!              - Specific treatment for the last time-step
(  496)       !!
(  497)       !! ** Method  :   3 types available:
(  498)       !!      a) Use named pipes(FIFO) to exchange process identifiers
(  499)       !!         between the 3 programs
(  500)       !!      b) USE a messag passing method with PVM language (CLIM)
(  501)       !!      c) USE SVIPC method
(  502)       !!
(  503)       !! History :
(  504)       !!        !  92-09 (L. Terray)  Original code
(  505)       !!        !  96-07 (L. Terray)  OASIS version 2
(  506)       !!        !  96-11 (E. Guilyardi)  run-off + Taux,Tauy x2
(  507)       !!        !  98-04 (M.A Foujols, S. Valcke, M. Imbard)  OASIS2.2
(  508)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  509)       !!----------------------------------------------------------------------
(  510)       !! * modules used






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 10

(  511)       USE ioipsl
(  512)       USE phycst          ! physical constants
(  513) 
(  514)       !! * Arguments
(  515)       INTEGER, INTENT( in ) ::   kt    ! ocean time-step
(  516) 
(  517)       !! * Local declarations
(  518)       INTEGER :: ji, jj, jn, jf           ! dummy loop indexes
(  519)       INTEGER :: icstep, info, ierror, isize
(  520)       INTEGER :: iflmax, iunmax
(  521)       INTEGER :: ifile(jpmaxfld), ifield(jpmaxfld)
(  522)       CHARACTER (len=8) ::  clfile(jpmaxfld) 
(  523)       LOGICAL :: llfind
(  524)       REAL(wp), DIMENSION(jpiglo,jpjglo) ::    &
(  525)          zstoc, zieoc, zalboc, zticoc
(  526)       
(  527)       ! netcdf outputs
(  528)        
(  529)       CHARACTER (len=80) ::  clcplsnam
(  530)       INTEGER, SAVE ::  nhoridcs, nidcs, ndexcs(1)
(  531)       LOGICAL, SAVE :: lfirsts = .TRUE.
(  532)       REAL(wp) ::    zjulian
(  533)       
(  534)       ! Additions for SVIPC
(  535)       
(  536)       INTEGER  :: index
(  537) !      INTEGER, DIMENSION(3) :: infos
(  538)       CHARACTER (len=3) ::  clmodinf       ! Header or not
(  539) !      CHARACTER (len=3) ::  cljobnam      ! experiment name
(  540)       !!----------------------------------------------------------------------
(  541) 
(  542)       ! coupled mode Ocean/Atmosphere
(  543) 
(  544)       ! 0. Initialization
(  545)       ! -----------------
(  546) 
(  547)       isize = jpiglo * jpjglo
(  548) 
(  549)       ! First time step: ocean sst and ice sea-ice extend set to zero
(  550)       IF( kt == nit000 ) THEN
(  551)          sstoc(:,:) = 0.e0
(  552)          sieoc(:,:) = 0.e0
(  553)          alboc(:,:) = 0.e0
(  554)          ticoc(:,:) = 0.e0
(  555) 
(  556)          ! initialisation for netcdf outputs
(  557)          ! 
(  558)          ndexcs(:) = 0
(  559)          clcplsnam = "cpl_oce_sst"
(  560) 
(  561)          ! Compute julian date from starting date of the run
(  562)          CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
(  563)          CALL histbeg(clcplsnam, jpiglo, glamt, jpjglo, gphit,   &
(  564)             1, jpiglo, 1, jpjglo, 0,   &
(  565)             zjulian, rdt, nhoridcs, nidcs, domain_id=nidom)
(  566)          ! no vertical axis
(  567)          DO jf = 1, nfldo2c
(  568)             CALL histdef(nidcs, cpl_writ(jf),cpl_writ(jf),"-",jpi,    &






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 11

(  569)                jpj, nhoridcs, 1, 1, 1, -99, 32, "inst", rdt, rdt)
(  570)          END DO
(  571)          CALL histend(nidcs)
(  572)       ENDIF
(  573) 
(  574)       ! 1. Cumulated sst and sea-ice extend
(  575)       !------------------------------------
(  576) 
(  577)       sstoc(:,:) = sstoc(:,:) + ( 1.0 - freeze(:,:) ) * ( tn(:,:,1) + rt0 )
(  578)       sieoc(:,:) = sieoc(:,:) + freeze(:,:)
(  579) 
(  580) #if defined key_ice_lim
(  581)       alboc(:,:) = alboc(:,:) + freeze(:,:) * alb_ice(:,:)
(  582)       ticoc(:,:) = ticoc(:,:) + freeze(:,:) * tn_ice(:,:) 
(  583) #else
(  584)       alboc(:,:) = alboc(:,:) + freeze(:,:) * 0.8
(  585)       ticoc(:,:) = ticoc(:,:) + freeze(:,:) * ( -10.e0 + rt0 )
(  586) #endif
(  587) 
(  588) 
(  589)       ! 2. Send coupling fields to OASIS
(  590)       !---------------------------------
(  591) 
(  592)       IF( MOD( kt, nexco ) == 0 ) THEN
(  593) 
(  594)          ! 2.1 Average : mean coupling fields
(  595)          zstoc (:,:) = 0.e0
(  596)          zieoc (:,:) = 0.e0
(  597)          zalboc(:,:) = 0.e0
(  598)          zticoc(:,:) = 0.e0
(  599)          DO jj = 1, nlcj
(  600)             DO ji = 1, nlci
(  601)                zstoc (mig(ji),mjg(jj)) = sstoc(ji,jj) / FLOAT( nexco )
(  602)                zieoc (mig(ji),mjg(jj)) = sieoc(ji,jj) / FLOAT( nexco )
(  603)                zalboc(mig(ji),mjg(jj)) = alboc(ji,jj) / FLOAT( nexco )
(  604)                zticoc(mig(ji),mjg(jj)) = ticoc(ji,jj) / FLOAT( nexco )
(  605)             END DO
(  606)          END DO
(  607)          icstep = kt - nit000 + 1
(  608) 
(  609)          if(lwp) then
(  610)             WRITE(numout,*)
(  611)             WRITE(numout,*) 'STEP: Send fields to CPL with kt= ', kt
(  612)             WRITE(numout,*)
(  613)          endif
(  614) 
(  615)          ! outputs
(  616) 
(  617)          CALL histwrite( nidcs, cpl_writ(1), kt, zstoc , jpi*jpj, ndexcs )
(  618)          CALL histwrite( nidcs, cpl_writ(2), kt, zieoc , jpi*jpj, ndexcs )
(  619)          CALL histwrite( nidcs, cpl_writ(3), kt, zalboc, jpi*jpj, ndexcs )
(  620)          CALL histwrite( nidcs, cpl_writ(4), kt, zticoc, jpi*jpj, ndexcs )
(  621)          CALL histsync ( nidcs )
(  622) 
(  623)          ! 2.2 Last time step (clim or pipe) or pipe mode
(  624)          ! 
(  625)          IF( kt == nitend .OR. cchan == 'PIPE' ) THEN 
(  626) 






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 12

(  627)             ! finalize outputs
(  628) 
(  629)             CALL histclo( nidcs )
(  630) 
(  631)             ! WRITE fields for coupler with pipe technique or for last time step
(  632) 
(  633)             ! initialisation
(  634) 
(  635)             iflmax =  1
(  636)             iunmax = 99
(  637)             
(  638)             clfile(iflmax) = cpl_f_writ(iflmax)     ! keeps first file name
(  639)             ifile(iflmax) = iunmax                  ! keeps first file unit
(  640)             iunmax = iunmax - 1                     ! decrements file unit maximum
(  641)             ifield(1) = ifile(iflmax)               ! keeps file unit for field
(  642) 
(  643)             ! different files names counter
(  644)             DO jf = 2, nfldo2c
(  645)                llfind = .FALSE.
(  646)                DO jn = 1, iflmax
(  647)                   IF( .NOT. llfind ) THEN
(  648)                      IF( cpl_f_writ(jf) == clfile(jn) ) THEN
(  649)                         ifield(jf) = ifile(jn)      ! keep file unit for field
(  650)                         llfind = .TRUE.
(  651)                      ENDIF
(  652)                   END IF
(  653)                END DO
(  654)                IF( .NOT. llfind) THEN
(  655)                   iflmax = iflmax + 1               ! increment the number of different files
(  656)                   clfile(iflmax) = cpl_f_writ(jf)   ! keep file name
(  657)                   ifile (iflmax) = iunmax           ! keep file unit for file
(  658)                   ifield(jf) = ifile(iflmax)        ! keep file unit for field
(  659)                   iunmax = iunmax-1                 ! decrement unit maximum number from 99 to 98...
(  660)                ENDIF
(  661)             END DO
(  662)             !          
(  663)             DO jn = 1, iflmax 
(  664)                OPEN (ifile(jn), FILE=clfile(jn), FORM='UNFORMATTED')
(  665)             END DO
(  666)             !              
(  667)             DO jf = 1, nfldo2c
(  668)                IF( jf == 1 ) CALL locwrite(cpl_writ(jf),zstoc , isize, ifield(jf), ierror, numout) 
(  669)                IF( jf == 2 ) CALL locwrite(cpl_writ(jf),zieoc , isize, ifield(jf), ierror, numout) 
(  670)                IF( jf == 3 ) CALL locwrite(cpl_writ(jf),zalboc, isize, ifield(jf), ierror, numout) 
(  671)                IF( jf == 4 ) CALL locwrite(cpl_writ(jf),zticoc, isize, ifield(jf), ierror, numout) 
(  672)             END DO
(  673) 
(  674)             ! simulate a FLUSH
(  675) 
(  676)             DO jn = 1, iflmax 
(  677)                CLOSE( ifile(jn) )
(  678)             END DO
(  679) 
(  680)             ! Clim mode
(  681)             IF( cchan == 'CLIM' ) THEN  ! inform PVM daemon, I have finished
(  682)                CALL CLIM_Quit( CLIM_ContPvm, info )
(  683)                IF( info /= CLIM_Ok ) THEN
(  684)                   WRITE (6, *) 'An error occured while leaving CLIM. Error = ',info






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 13

(  685)                ENDIF
(  686)             ENDIF
(  687) 
(  688)          ENDIF
(  689) 
(  690)          ! IF last we have finished if not pass info to the atmosphere
(  691) 
(  692)          IF ( kt /= nitend ) THEN
(  693) 
(  694)             ! 2.3 normal exchange
(  695) 
(  696)             ! PIPE mode      
(  697)             IF( cchan == 'PIPE' ) THEN 
(  698) 
(  699)                ! Send message to pipes for CPL=atmosphere
(  700) 
(  701)                ! DO jf = 1, nfldo2c
(  702)                !    CALL PIPE_Model_Send(cpl_writ(jf), icstep, numout)
(  703)                ! END DO 
(  704) 
(  705)                ! SIPC mode
(  706)             ELSE IF( cchan == 'SIPC' ) THEN
(  707) 
(  708)                ! Define IF a header must be encapsulated within the field brick :
(  709)                clmodinf = 'NOT'  ! as $MODINFO in namcouple  
(  710) 
(  711)                ! IF clmodinf = 'YES', define encapsulated infos to be exchanged
(  712)                !    infos(1) = initial date
(  713)                !    infos(2) = timestep
(  714)                !    infos(3) = actual time
(  715)                !
(  716)                ! Writing of output field SST SOSSTSST
(  717)                !
(  718)                ! Index of SST in total number of fields jpfldo2a: 
(  719)                index = 1
(  720)                !
(  721)                ! CALL SIPC_Write_Model(index, isize, clmodinf, cljobnam, infos, zstoc)
(  722)                !
(  723)                ! Writing of output field Sea-Ice SOICECOV 
(  724)                !
(  725)                ! Index of sea-ice in total number of fields jpfldo2a: 
(  726)                index = 2
(  727)                !
(  728)                ! CALL SIPC_Write_Model(index, isize, clmodinf, cljobnam, infos, zieoc)
(  729)    
(  730)                ! CLIM mode
(  731)             ELSE IF( cchan == 'CLIM' ) THEN
(  732)    
(  733)                DO jn = 1, nfldo2c
(  734)    
(  735)                   IF (jn == 1) CALL CLIM_Export(cpl_writ(jn), icstep, zstoc , info)
(  736)                   IF (jn == 2) CALL CLIM_Export(cpl_writ(jn), icstep, zieoc , info)
(  737)                   IF (jn == 3) CALL CLIM_Export(cpl_writ(jn), icstep, zalboc, info)
(  738)                   IF (jn == 4) CALL CLIM_Export(cpl_writ(jn), icstep, zticoc, info)
(  739) 
(  740)                   IF (info /= CLIM_Ok) THEN
(  741)                      WRITE (numout,*) 'STEP : Pb giving', cpl_writ(jn), ':', jn
(  742)                      WRITE (numout,*) ' at timestep = ', icstep, 'kt = ', kt






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 14

(  743)                      WRITE (numout,*) 'Clim error code is = ',info
(  744)                      WRITE (numout,*) 'STOP in stpcpl '
(  745)                      CALL abort(' stpcpl ')
(  746)                   ENDIF
(  747)                END DO
(  748)             ENDIF
(  749) 
(  750)             ! reset cumulative sst and sea-ice extend to zero
(  751)             sstoc(:,:) = 0.e0
(  752)             sieoc(:,:) = 0.e0
(  753)             alboc(:,:) = 0.e0
(  754)             ticoc(:,:) = 0.e0
(  755)          ENDIF
(  756)       ENDIF
(  757) 
(  758)    END SUBROUTINE cpl_stp
(  759) 
(  760) #else
(  761)    !!----------------------------------------------------------------------
(  762)    !!   Default case           Dummy module         forced Ocean/Atmosphere
(  763)    !!----------------------------------------------------------------------
(  764)    USE in_out_manager
(  765) CONTAINS
(  766)    SUBROUTINE cpl_init            ! Dummy routine
(  767)       if(lwp) WRITE(numout,*) 'cpl_init: You should have not see this print! error?'
(  768)    END SUBROUTINE cpl_init
(  769)    SUBROUTINE cpl_stp( kt )       ! Dummy routine
(  770)       if(lwp) WRITE(numout,*) 'cpl_stp: You should have not see this print! error?', kt
(  771)    END SUBROUTINE cpl_stp
(  772) #endif
(  773) 
(  774)    !!======================================================================
(  775) END MODULE cpl































PGF90 (Version     10.5)          09/23/2010  09:49:27      page 15

(    1) # 1 "cpl.F90"
(    1) MODULE cpl
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE cpl  ***
(    4)    !! Coupled O/A : coupled ocean-atmosphere case using OASIS
(    5)    !!=====================================================================
(    6) # 761
(  761)    !!----------------------------------------------------------------------
(  762)    !!   Default case           Dummy module         forced Ocean/Atmosphere
(  763)    !!----------------------------------------------------------------------
(  764)    USE in_out_manager
(  765) CONTAINS
(  766)    SUBROUTINE cpl_init            ! Dummy routine
(  767)       if(lwp) WRITE(numout,*) 'cpl_init: You should have not see this print! error?'
(  768)    END SUBROUTINE cpl_init
(  769)    SUBROUTINE cpl_stp( kt )       ! Dummy routine
(  770)       if(lwp) WRITE(numout,*) 'cpl_stp: You should have not see this print! error?', kt
(  771)    END SUBROUTINE cpl_stp
(  772) 
(  773) # 774
(  774)    !!======================================================================
(  775) END MODULE cpl
