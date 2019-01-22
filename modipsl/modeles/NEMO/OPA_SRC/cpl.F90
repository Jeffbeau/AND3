MODULE cpl
   !!======================================================================
   !!                       ***  MODULE cpl  ***
   !! Coupled O/A : coupled ocean-atmosphere case using OASIS
   !!=====================================================================
#if defined key_coupled
   !!----------------------------------------------------------------------
   !!   'key_coupled'                              coupled Ocean/Atmosphere
   !!----------------------------------------------------------------------
   !!   cpl_init     : initialization of coupled mode communication
   !!   cpl_read     : read the coupled namelist
   !!   cpl_stp      : exchange fields in coupled mode
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE cpl_oce         ! coupled exchange variables (???)
   USE flx_oce         ! in case of ice 
   USE ocfzpt          ! ocean freezing point
   USE daymod          ! calendar

   IMPLICIT NONE
   PRIVATE

   !! Routine accessibility
   PUBLIC cpl_init     ! routine called in opa module
   PUBLIC cpl_stp      ! routine called in step module
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/cpl.F90,v 1.6 2005/12/12 14:18:01 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE cpl_init( kastp, kexch, kstep )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE cpl_init  ***
      !!
      !! ** Purpose :   Initialize coupled mode communication for ocean
      !!    exchange process identifiers and timestep information
      !!    between AGCM, OGCM and COUPLER. (OASIS software)
      !!
      !! ** Method  :  3 types :
      !!      a) Use named pipes(FIFO) to exchange process identifiers
      !!          between the 3 programs
      !!      b) USE a messag passing method with PVM language (CLIM)
      !!      c) USE SVIPC method
      !!
      !! ** Input   :   npiat     : agcm process identifier
      !!                npicp     : coupler process identifier
      !!                npioc     : ogcm process identifier
      !!
      !! Reference  :   see Epicoa 0803 (1992)
      !!
      !! History :
      !!        !  92-09  (L. Terray)  Original code
      !!        !  96-07  (L. Terray)  OASIS version 2
      !!        !  96-11  (E. Guilyardi)  run-off + Taux,Tauy x2
      !!        !  98-04  (M.A Foujols, S. Valcke, M. Imbard)  OASIS2.2
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(in ) ::   &
         kastp,      &  ! total number of timesteps in oceanic model
         kexch,      &  ! frequency of exchange for the fluxes (in time steps)
         kstep          ! timestep value (in seconds)

      !! * Local declarations
      INTEGER,DIMENSION(3)  :: iparal
      INTEGER               :: ifcpl, idt, info, imxtag, istep

      CHARACTER (len=9) ::   clpoolnam 
      INTEGER           :: ipoolhandle, ipoolsize, jf
      CHARACTER (len=3) ::   cljobnam      ! experiment name
      INTEGER           :: ierror
!      INTEGER,DIMENSION(4) ::  imess
      INTEGER,DIMENSION(4) ::  imesso
      !!----------------------------------------------------------------------
     
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cpl_init : initialization in coupled ocean/atmosphere case'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      IF(lwp) WRITE(numout,*)
     
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily || defined key_flx_forced_daily
      IF(lwp)WRITE(numout,cform_err)
      IF(lwp)WRITE(numout,*) ' key_coupled and key_flx_bulk_... are incompatible'
      nstop = nstop + 1
#endif
 
      IF(lwp)WRITE(numout,*)'     coupled simulation'
      IF(lwp)WRITE(numout,*)'        unit ',numlhf,' receive atm fluxes'
      IF(lwp)WRITE(numout,*)'        unit ',numlws,' receive windstress'
      IF(lwp)WRITE(numout,*)'        unit ',numlts,' transfer sst'
      IF(lwp)WRITE(numout,*)'        unit ',numlic,' transfer ice cover'


      CALL cpl_read           ! read the coupled mode information in namelist

      CALL flush(numout)

      ! I- PIPE
      ! --------
      ! W A R N I N G : PIPE technique is temporary diseable (nov. 2000)

      IF( cchan == 'PIPE' ) THEN 

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'Making pipes for fields to receive from CPL'
         IF(lwp) WRITE(numout,*)
         CALL flush(numout)
         ierror = 0

         ! loop to define pipes (CPL=atmos to ocean)

         DO jf = 1, nflxc2o  
            ! CALL PIPE_Model_Define( numout, cpl_readflx(jf), jpread, info )
            IF( info /= 0 ) ierror = ierror + 1
         END DO
         DO jf = 1, ntauc2o
            ! CALL PIPE_Model_Define( numout, cpl_readtau(jf), jpread, info )
            IF( info /= 0 ) ierror = ierror + 1
         END DO
         
         IF(lwp) WRITE(numout,*) ' '
         IF(lwp) WRITE(numout,*) 'Making pipes for fields to send to CPL'
         IF(lwp) WRITE(numout,*) ' '
         
         ! loop to define pipes (ocean to atmos=CPL)

         DO jf = 1, nfldo2c
            ! CALL PIPE_Model_Define( numout, cpl_writ(jf), jpwrit, info )
            IF( info /= 0 ) ierror = ierror + 1
         END DO

         IF( ierror /= 0 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'cpl_init: end of job due to error in pipes definitions'
            CALL abort('')
         END IF
         
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'All pipes have been made'
         
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'Communication test between OCE and CPL'
         CALL flush(numout)
         
         ! CALL PIPE_Model_Stepi(numout, imess, nmodcpl, imesso, ierror)
         
         IF( ierror /= 0 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'cpl_init: end of job due to error in exchange first informations with Oasis'
            CALL abort('')
         END IF
         
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'Communication test between OCE and CPL is OK'
         IF(lwp) WRITE(numout,*) ' total simulation time in oasis = ',imesso(1)
         IF(lwp) WRITE(numout,*) ' total number of iterations is  = ',imesso(2)
         IF(lwp) WRITE(numout,*) ' value of oasis timestep  is    = ',imesso(3)
         IF(lwp) WRITE(numout,*) ' process id for oasis  is       = ',imesso(4)
         CALL flush(numout)
         
         ! II SVIPC
         ! ---------
         ! W A R N I N G : SVIPC technique is temporary diseable (nov. 2000)
         
         
      ELSE IF( cchan == 'SIPC' ) THEN

         ! debug for more information

         ! CALL SVIPC_debug(1)

         ! Define the experiment name :

          cljobnam = 'IPC'      ! as $JOBNAM in namcouple

          ! Attach to shared memory pool used to exchange initial infos 

          info = 0
          ! CALL SIPC_Init_Model (cljobnam, cplmodnam, 1, info)
          IF( info /= 0 ) THEN 
             IF(lwp) WRITE(numout,*)
             IF(lwp) WRITE(numout,*)'WARNING: Problem with attachement to',info 
             IF(lwp) WRITE(numout,*)'         initial memory pool(s) in ocean'
             IF(lwp) WRITE(numout,*)
             CALL abort('STOP in ocean')
          ENDIF

          ! Attach to pools used to exchange fields from ocean to coupler

          DO jf = 1, nfldo2c
             ! Pool name:
             clpoolnam = 'P'//cpl_writ(jf)
             ! CALL SIPC_Attach(clpoolnam, ipoolhandle)
             ! Resulting pool handle:
             mpoolwrit(jf) = ipoolhandle  
          END DO

          ! Attach to pools used to exchange fields from coupler to ocean
          
          DO jf = 1, nflxc2o
             ! Pool name:
             clpoolnam = 'P'//cpl_readflx(jf)
             ! CALL SIPC_Attach(clpoolnam, ipoolhandle)
             ! Resulting pool handle:
             mpoolread(jf) = ipoolhandle  
          END DO  

          DO jf = 1, ntauc2o
             ! Pool name:
             clpoolnam = 'P'//cpl_readtau(jf)
             ! CALL SIPC_Attach(clpoolnam, ipoolhandle)
             ! Resulting pool handle:
             mpoolread(jf+nflxc2o) = ipoolhandle  
          END DO 

          ! Exchange of initial infos

          ! Write data array isend to pool READ by Oasis

          info = 0
          ipoolsize = 4*jpbyteint
          ! CALL SVIPC_write(mpoolinitr, imess, ipoolsize, info)

          ! Find error if any

          IF( info < 0 ) THEN
             IF(lwp) WRITE(numout,*)
             IF(lwp) WRITE(numout,*) 'Problem in ocean in writing initial' 
             IF(lwp) WRITE(numout,*) 'infos to the shared memory segment(s)'
             IF(lwp) WRITE(numout,*)
          ELSE
             IF(lwp) WRITE(numout,*)
             IF(lwp) WRITE(numout,*) 'Initial infos written in ocean'            
             IF(lwp) WRITE(numout,*) 'to the shared memory segment(s)'
             IF(lwp) WRITE(numout,*)
          ENDIF

          ! Read data array irecv from pool written by Oasis

          info = 0
          ipoolsize = 4 * jpbyteint
          CALL svipc_read(mpoolinitw, imesso, ipoolsize, info)

          ! Find error if any

          IF( info < 0 ) THEN
             IF(lwp) WRITE(numout,*) '   '
             IF(lwp) WRITE(numout,*) 'Problem in ocean in reading initial' 
             IF(lwp) WRITE(numout,*) 'infos from the shared memory segment(s)'
             IF(lwp) WRITE(numout,*) '   '
          ELSE
             IF(lwp) WRITE(numout,*) '   '
             IF(lwp) WRITE(numout,*) 'Initial infos read by ocean'                
             IF(lwp) WRITE(numout,*) 'from the shared memory segment(s)'
             IF(lwp) WRITE(numout,*) '   '
             IF(lwp) WRITE(numout,*) ' ntime, niter, nstep, Oasis pid:'
             IF(lwp) WRITE(numout,*) imesso(1), imesso(2), imesso(3), imesso(4) 
          ENDIF

          ! Detach from shared memory segment(s)

          info = 0
          ! CALL SVIPC_close(mpoolinitw, 0, info)
          
          ! Find error if any

          IF( info < 0 ) THEN
             IF(lwp) WRITE(numout,*) 'Problem in detaching from shared memory segment(s)'
             IF(lwp) WRITE(numout,*) 'used by ocean to read initial infos' 
          ENDIF

          ! III CLIM
          ! --------

      ELSE IF( cchan == 'CLIM' ) THEN

         ! Define the number of processors involved in the coupling for
         ! Oasis (=1) and each model (as last two INTEGER on $CHATYPE line
         ! in the namcouple); they will be stored in a COMMON in mpiclim.h
         ! (used for CLIM/MPI2 only)
         mpi_nproc(0)=1
         mpi_nproc(1)=1
         mpi_nproc(2)=1 

         ! Define the number of processors used by each model as in
         ! $CHATYPE line of namcouple (used for CLIM/MPI2 only)
         mpi_totproc(1)=1
         mpi_totproc(2)=1
         
         ! Define names of each model as in $NBMODEL line of namcouple
         ! (used for CLIM/MPI2 only)
         cmpi_modnam(1)='lmdz.x'
         cmpi_modnam(2)=cplmodnam
         !  
         ! 1.1-Define the experiment name :
         
         cljobnam = 'CLI'      ! as $JOBNAM in namcouple
         
         OPEN ( UNIT = 7, FILE = 'trace', STATUS = 'unknown', FORM = 'formatted')
         CALL clim_init ( cljobnam, cplmodnam, 3, 7,   &
                          kastp, kexch, kstep,   &
                          5, 3600, 3600, info )

         IF( info /= clim_ok ) THEN
            IF(lwp) WRITE( numout, *) 'cpl_init : pb init clim, error code is = ', info
            CALL abort( 'STOP in cpl_init' )
         ELSE
            IF(lwp) WRITE(numout,*) 'cpl_init : init clim ok '
         ENDIF
         
         iparal ( clim_strategy ) = clim_serial
         iparal ( clim_length   ) = jpiglo*jpjglo
         iparal ( clim_offset   ) = 0
         
         ! loop to define messages (CPL=atmos to ocean)
         DO jf = 1, nflxc2o
            CALL CLIM_Define ( cpl_readflx(jf), clim_in, clim_double, iparal, info )  
         END DO
         DO jf = 1, ntauc2o
            CALL CLIM_Define ( cpl_readtau(jf), clim_in, clim_double, iparal, info )  
         END DO
         
         ! loop to define messages (ocean to CPL=atmos)
         DO jf = 1, nfldo2c
            CALL CLIM_Define ( cpl_writ(jf), clim_out, clim_double, iparal, info )   
         END DO
         
         IF(lwp) WRITE(numout,*) 'cpl_init : clim_define ok '
         
         CALL CLIM_Start( imxtag, info )
         
         IF( info /= clim_ok ) THEN
            IF(lwp) WRITE(numout,*) 'cpl_init : pb start clim, error code is = ', info
            CALL abort( 'stop in cpl_init' )
         ELSE
            IF(lwp) WRITE(numout,*) 'cpl_init : start clim ok '
         ENDIF
         
         CALL CLIM_Stepi ( cploasis, istep, ifcpl, idt, info )

         IF( info /= clim_ok ) THEN
            IF(lwp) WRITE(numout,*) ' warning : problem in getting step info from oasis '
            IF(lwp) WRITE(numout,*) ' =======   error code number = ', info
         ELSE
            IF(lwp) WRITE(numout,*) ' got step information from oasis '
         ENDIF
         IF(lwp) WRITE(numout,*) ' number of tstep in oasis ', istep
         IF(lwp) WRITE(numout,*) ' exchange frequency in oasis ', ifcpl
         IF(lwp) WRITE(numout,*) ' length of tstep in oasis ', idt
      ENDIF

   END SUBROUTINE cpl_init


   SUBROUTINE cpl_read
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE cpl_read  ***
      !!                    
      !! ** Purpose :   Read and print options for the coupled run (namelist)
      !!
      !! ** Method  :   ???
      !!
      !! History :
      !!   8.5  !  02-12  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER :: jf

      NAMELIST/namcpl/ nexco, cchan, nmodcpl, cplmodnam, cploasis   &
          , nfldo2c, nflxc2o, ntauc2o, cpl_f_readflx, cpl_f_readtau   &
          , cpl_f_writ, cpl_readflx, cpl_readtau, cpl_writ
      !!----------------------------------------------------------------------
      
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' cpl_read : read the coupled parameters in namelist'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~'
      IF(lwp) WRITE(numout,*)

      ! Default values
      
      nexco = 24
      cchan='PIPE'              ! echange TYPE
      nmodcpl = 2
      cplmodnam = 'opa.xx'      ! model name : as $NBMODEL in namcouple
      cploasis = 'Oasis'        ! coupler name : as in coupler

      ! -Define symbolic name for fields exchanged from ocean to coupler,
      ! must be the same as (1) of the field  definition in namcouple:
      nfldo2c=2
      cpl_writ(1)='SOSSTSST'
      cpl_writ(2)='SOICECOV'

      ! -Define files name for fields exchanged from ocean to coupler,
      ! must be the same as (6) of the field  definition in namcouple:
      nflxc2o=6
      cpl_readflx(1)='SONSFLDO' ! non solar heat flux (positif vers l'ocean)
      cpl_readflx(2)='SOSHFLDO' ! solar flux
      cpl_readflx(3)='SOTOPRSU' ! precip
      cpl_readflx(4)='SOTFSHSU' ! evaporation
      cpl_readflx(5)='SORUNCOA' ! coastal run-off
      cpl_readflx(6)='SORIVFLU' ! river run-off
      ntauc2o=2
      cpl_readtau(1)='SOZOTAUX' ! first zonal wind stress
      cpl_readtau(2)='SOZOTAU2' ! second zonal wind stress
      cpl_readtau(3)='SOMETAUY' ! first meridien wind stress
      cpl_readtau(4)='SOMETAU2' ! second meridien wind stress

      ! -Define files name for fields exchanged from ocean to coupler,
      ! must be the same as (6) of the field  definition in namcouple:
      cpl_f_writ(1)='ocesst'
      cpl_f_writ(2)='oceice'

      ! -Define files name for fields exchanged from coupler to ocean
      ! must be the same as (7) of the field  definition in namcouple:
      cpl_f_readflx(1)='oceflx'
      cpl_f_readflx(2)='oceflx'
      cpl_f_readflx(3)='oceflx'
      cpl_f_readflx(4)='oceflx'
      cpl_f_readflx(5)='oceflx'
      cpl_f_readflx(6)='oceflx'
      cpl_f_readtau(1)='ocetau'
      cpl_f_readtau(2)='ocetau'
      cpl_f_readtau(3)='ocetau'
      cpl_f_readtau(4)='ocetau'

      ! Namelist namcpl : coupling mode and frequency
      REWIND( numnam )
      READ  ( numnam, namcpl )

      IF(lwp) THEN
         WRITE(numout,*) 'namcpl'
         WRITE(numout,*) 
         WRITE(numout,*) ' Coupling exchange frequency    nexco  = ',nexco
         WRITE(numout,*) ' Coupling exchange technique    cchan  = ',cchan
         WRITE(numout,*) ' Mode Coupling technique      nmodcpl  = ',nmodcpl
         WRITE(numout,*) ' Define the model name      cplmodnam  = ',cplmodnam
         WRITE(numout,*) ' Define the coupler name      cploasis = ',cploasis
         WRITE(numout,*) ' Fields number ocean to coupler nfldo2c= ',nfldo2c
         WRITE(numout,*) ' Flux fields coupler to ocean nflxc2o  = ',nflxc2o
         WRITE(numout,*) ' Stress fields coupler to ocean ntauc2o= ',ntauc2o
         IF ( cchan == 'PIPE' .OR.  cchan == 'pipe' ) THEN
            cchan='PIPE'
            WRITE(numout,*)'(communication between models made by pipes)'
         ELSEIF( cchan == 'CLIM' .OR. cchan == 'clim' ) THEN
            cchan='CLIM'
            WRITE(numout,*)'(communication between models made by CLIM-PVM)'
         ELSEIF( cchan == 'SIPC' .OR. cchan == 'sipc' ) THEN
            cchan='SIPC'
            WRITE(numout,*)'(communication between models made by the',  &
               'Share Memory Segment and Semaphore)'
         ELSE
            WRITE(numout,*) 'technic not yet implemented ',cchan
            STOP 'in cpl_read'
         ENDIF
         DO jf = 1, nflxc2o
            WRITE(numout,*) 'file to receive field number = ',jf,'  ',cpl_f_readflx(jf) 
         END DO
         DO jf = 1, ntauc2o
            WRITE(numout,*) 'file to receive field number = ',jf,'  ',cpl_f_readtau(jf) 
         END DO
         DO jf = 1, nfldo2c
            WRITE(numout,*) 'file to send field number = ',jf,'  ',cpl_f_writ(jf)
         END DO
         WRITE(numout,*) ' fields received from coupler'
         DO jf = 1, nflxc2o
            WRITE(numout,*) 'symbolic name for field number = ',jf,'  ',cpl_readflx(jf) 
         END DO
         DO jf = 1, ntauc2o
            WRITE(numout,*) 'symbolic name for field number = ',jf,'  ',cpl_readtau(jf) 
         END DO
         WRITE(numout,*) ' fields send to coupler'
         DO jf = 1, nfldo2c
            WRITE(numout,*) 'symbolic name for field number = ',jf,'  ',cpl_writ(jf)
         END DO
      ENDIF

   END SUBROUTINE cpl_read


   SUBROUTINE cpl_stp( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE cpl_stp  ***
      !!                      *****************
      !!                      * OASIS routine *
      !!                      *****************
      !! ** Purpose : - At each coupling time-step,this routine sends fields
      !!      like sst or ice cover to the coupler.
      !!              - At each time-step computes time average values 
      !!              - Specific treatment for the last time-step
      !!
      !! ** Method  :   3 types available:
      !!      a) Use named pipes(FIFO) to exchange process identifiers
      !!         between the 3 programs
      !!      b) USE a messag passing method with PVM language (CLIM)
      !!      c) USE SVIPC method
      !!
      !! History :
      !!        !  92-09 (L. Terray)  Original code
      !!        !  96-07 (L. Terray)  OASIS version 2
      !!        !  96-11 (E. Guilyardi)  run-off + Taux,Tauy x2
      !!        !  98-04 (M.A Foujols, S. Valcke, M. Imbard)  OASIS2.2
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * modules used
      USE ioipsl
      USE phycst          ! physical constants

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step

      !! * Local declarations
      INTEGER :: ji, jj, jn, jf           ! dummy loop indexes
      INTEGER :: icstep, info, ierror, isize
      INTEGER :: iflmax, iunmax
      INTEGER :: ifile(jpmaxfld), ifield(jpmaxfld)
      CHARACTER (len=8) ::  clfile(jpmaxfld) 
      LOGICAL :: llfind
      REAL(wp), DIMENSION(jpiglo,jpjglo) ::    &
         zstoc, zieoc, zalboc, zticoc
      
      ! netcdf outputs
       
      CHARACTER (len=80) ::  clcplsnam
      INTEGER, SAVE ::  nhoridcs, nidcs, ndexcs(1)
      LOGICAL, SAVE :: lfirsts = .TRUE.
      REAL(wp) ::    zjulian
      
      ! Additions for SVIPC
      
      INTEGER  :: index
!      INTEGER, DIMENSION(3) :: infos
      CHARACTER (len=3) ::  clmodinf       ! Header or not
!      CHARACTER (len=3) ::  cljobnam      ! experiment name
      !!----------------------------------------------------------------------

      ! coupled mode Ocean/Atmosphere

      ! 0. Initialization
      ! -----------------

      isize = jpiglo * jpjglo

      ! First time step: ocean sst and ice sea-ice extend set to zero
      IF( kt == nit000 ) THEN
         sstoc(:,:) = 0.e0
         sieoc(:,:) = 0.e0
         alboc(:,:) = 0.e0
         ticoc(:,:) = 0.e0

         ! initialisation for netcdf outputs
         ! 
         ndexcs(:) = 0
         clcplsnam = "cpl_oce_sst"

         ! Compute julian date from starting date of the run
         CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
         CALL histbeg(clcplsnam, jpiglo, glamt, jpjglo, gphit,   &
            1, jpiglo, 1, jpjglo, 0,   &
            zjulian, rdt, nhoridcs, nidcs, domain_id=nidom)
         ! no vertical axis
         DO jf = 1, nfldo2c
            CALL histdef(nidcs, cpl_writ(jf),cpl_writ(jf),"-",jpi,    &
               jpj, nhoridcs, 1, 1, 1, -99, 32, "inst", rdt, rdt)
         END DO
         CALL histend(nidcs)
      ENDIF

      ! 1. Cumulated sst and sea-ice extend
      !------------------------------------

      sstoc(:,:) = sstoc(:,:) + ( 1.0 - freeze(:,:) ) * ( tn(:,:,1) + rt0 )
      sieoc(:,:) = sieoc(:,:) + freeze(:,:)

#if defined key_ice_lim
      alboc(:,:) = alboc(:,:) + freeze(:,:) * alb_ice(:,:)
      ticoc(:,:) = ticoc(:,:) + freeze(:,:) * tn_ice(:,:) 
#else
      alboc(:,:) = alboc(:,:) + freeze(:,:) * 0.8
      ticoc(:,:) = ticoc(:,:) + freeze(:,:) * ( -10.e0 + rt0 )
#endif


      ! 2. Send coupling fields to OASIS
      !---------------------------------

      IF( MOD( kt, nexco ) == 0 ) THEN

         ! 2.1 Average : mean coupling fields
         zstoc (:,:) = 0.e0
         zieoc (:,:) = 0.e0
         zalboc(:,:) = 0.e0
         zticoc(:,:) = 0.e0
         DO jj = 1, nlcj
            DO ji = 1, nlci
               zstoc (mig(ji),mjg(jj)) = sstoc(ji,jj) / FLOAT( nexco )
               zieoc (mig(ji),mjg(jj)) = sieoc(ji,jj) / FLOAT( nexco )
               zalboc(mig(ji),mjg(jj)) = alboc(ji,jj) / FLOAT( nexco )
               zticoc(mig(ji),mjg(jj)) = ticoc(ji,jj) / FLOAT( nexco )
            END DO
         END DO
         icstep = kt - nit000 + 1

         if(lwp) then
            WRITE(numout,*)
            WRITE(numout,*) 'STEP: Send fields to CPL with kt= ', kt
            WRITE(numout,*)
         endif

         ! outputs

         CALL histwrite( nidcs, cpl_writ(1), kt, zstoc , jpi*jpj, ndexcs )
         CALL histwrite( nidcs, cpl_writ(2), kt, zieoc , jpi*jpj, ndexcs )
         CALL histwrite( nidcs, cpl_writ(3), kt, zalboc, jpi*jpj, ndexcs )
         CALL histwrite( nidcs, cpl_writ(4), kt, zticoc, jpi*jpj, ndexcs )
         CALL histsync ( nidcs )

         ! 2.2 Last time step (clim or pipe) or pipe mode
         ! 
         IF( kt == nitend .OR. cchan == 'PIPE' ) THEN 

            ! finalize outputs

            CALL histclo( nidcs )

            ! WRITE fields for coupler with pipe technique or for last time step

            ! initialisation

            iflmax =  1
            iunmax = 99
            
            clfile(iflmax) = cpl_f_writ(iflmax)     ! keeps first file name
            ifile(iflmax) = iunmax                  ! keeps first file unit
            iunmax = iunmax - 1                     ! decrements file unit maximum
            ifield(1) = ifile(iflmax)               ! keeps file unit for field

            ! different files names counter
            DO jf = 2, nfldo2c
               llfind = .FALSE.
               DO jn = 1, iflmax
                  IF( .NOT. llfind ) THEN
                     IF( cpl_f_writ(jf) == clfile(jn) ) THEN
                        ifield(jf) = ifile(jn)      ! keep file unit for field
                        llfind = .TRUE.
                     ENDIF
                  END IF
               END DO
               IF( .NOT. llfind) THEN
                  iflmax = iflmax + 1               ! increment the number of different files
                  clfile(iflmax) = cpl_f_writ(jf)   ! keep file name
                  ifile (iflmax) = iunmax           ! keep file unit for file
                  ifield(jf) = ifile(iflmax)        ! keep file unit for field
                  iunmax = iunmax-1                 ! decrement unit maximum number from 99 to 98...
               ENDIF
            END DO
            !          
            DO jn = 1, iflmax 
               OPEN (ifile(jn), FILE=clfile(jn), FORM='UNFORMATTED')
            END DO
            !              
            DO jf = 1, nfldo2c
               IF( jf == 1 ) CALL locwrite(cpl_writ(jf),zstoc , isize, ifield(jf), ierror, numout) 
               IF( jf == 2 ) CALL locwrite(cpl_writ(jf),zieoc , isize, ifield(jf), ierror, numout) 
               IF( jf == 3 ) CALL locwrite(cpl_writ(jf),zalboc, isize, ifield(jf), ierror, numout) 
               IF( jf == 4 ) CALL locwrite(cpl_writ(jf),zticoc, isize, ifield(jf), ierror, numout) 
            END DO

            ! simulate a FLUSH

            DO jn = 1, iflmax 
               CLOSE( ifile(jn) )
            END DO

            ! Clim mode
            IF( cchan == 'CLIM' ) THEN  ! inform PVM daemon, I have finished
               CALL CLIM_Quit( CLIM_ContPvm, info )
               IF( info /= CLIM_Ok ) THEN
                  WRITE (6, *) 'An error occured while leaving CLIM. Error = ',info
               ENDIF
            ENDIF

         ENDIF

         ! IF last we have finished if not pass info to the atmosphere

         IF ( kt /= nitend ) THEN

            ! 2.3 normal exchange

            ! PIPE mode      
            IF( cchan == 'PIPE' ) THEN 

               ! Send message to pipes for CPL=atmosphere

               ! DO jf = 1, nfldo2c
               !    CALL PIPE_Model_Send(cpl_writ(jf), icstep, numout)
               ! END DO 

               ! SIPC mode
            ELSE IF( cchan == 'SIPC' ) THEN

               ! Define IF a header must be encapsulated within the field brick :
               clmodinf = 'NOT'  ! as $MODINFO in namcouple  

               ! IF clmodinf = 'YES', define encapsulated infos to be exchanged
               !    infos(1) = initial date
               !    infos(2) = timestep
               !    infos(3) = actual time
               !
               ! Writing of output field SST SOSSTSST
               !
               ! Index of SST in total number of fields jpfldo2a: 
               index = 1
               !
               ! CALL SIPC_Write_Model(index, isize, clmodinf, cljobnam, infos, zstoc)
               !
               ! Writing of output field Sea-Ice SOICECOV 
               !
               ! Index of sea-ice in total number of fields jpfldo2a: 
               index = 2
               !
               ! CALL SIPC_Write_Model(index, isize, clmodinf, cljobnam, infos, zieoc)
   
               ! CLIM mode
            ELSE IF( cchan == 'CLIM' ) THEN
   
               DO jn = 1, nfldo2c
   
                  IF (jn == 1) CALL CLIM_Export(cpl_writ(jn), icstep, zstoc , info)
                  IF (jn == 2) CALL CLIM_Export(cpl_writ(jn), icstep, zieoc , info)
                  IF (jn == 3) CALL CLIM_Export(cpl_writ(jn), icstep, zalboc, info)
                  IF (jn == 4) CALL CLIM_Export(cpl_writ(jn), icstep, zticoc, info)

                  IF (info /= CLIM_Ok) THEN
                     WRITE (numout,*) 'STEP : Pb giving', cpl_writ(jn), ':', jn
                     WRITE (numout,*) ' at timestep = ', icstep, 'kt = ', kt
                     WRITE (numout,*) 'Clim error code is = ',info
                     WRITE (numout,*) 'STOP in stpcpl '
                     CALL abort(' stpcpl ')
                  ENDIF
               END DO
            ENDIF

            ! reset cumulative sst and sea-ice extend to zero
            sstoc(:,:) = 0.e0
            sieoc(:,:) = 0.e0
            alboc(:,:) = 0.e0
            ticoc(:,:) = 0.e0
         ENDIF
      ENDIF

   END SUBROUTINE cpl_stp

#else
   !!----------------------------------------------------------------------
   !!   Default case           Dummy module         forced Ocean/Atmosphere
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE cpl_init            ! Dummy routine
      if(lwp) WRITE(numout,*) 'cpl_init: You should have not see this print! error?'
   END SUBROUTINE cpl_init
   SUBROUTINE cpl_stp( kt )       ! Dummy routine
      if(lwp) WRITE(numout,*) 'cpl_stp: You should have not see this print! error?', kt
   END SUBROUTINE cpl_stp
#endif

   !!======================================================================
END MODULE cpl
