MODULE iceini
   !!======================================================================
   !!                       ***  MODULE iceini   ***
   !!   Sea-ice model : LIM Sea ice model Initialization
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_init       : sea-ice model initialization
   !!----------------------------------------------------------------------
   USE dom_oce
   USE in_out_manager
   USE ice_oce         ! ice variables
   USE flx_oce
   USE phycst          ! Define parameters for the routines
   USE ocfzpt
   USE ice
   USE limmsh
   USE limistate
   USE limrst
   USE ini1d           ! initialization of the 1D configuration
!!DB
   USE restart


   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC ice_init                 ! called by opa.F90

   !! * Share Module variables
   LOGICAL , PUBLIC  ::   & !!! ** init namelist (namicerun) **
      ln_limdyn   = .TRUE.   !: flag for ice dynamics (T) or not (F)
   INTEGER , PUBLIC  ::   &  !:
      nstart ,            &  !: iteration number of the begining of the run 
      nlast  ,            &  !: iteration number of the end of the run 
      nitrun ,            &  !: number of iteration
      numit                  !: iteration number
   REAL(wp), PUBLIC  ::   &  !:
      hsndif = 0.e0 ,     &  !: computation of temp. in snow (0) or not (9999)
      hicdif = 0.e0 ,     &  !: computation of temp. in ice (0) or not (9999)
      tpstot                 !: time of the run in seconds
   REAL(wp), PUBLIC, DIMENSION(2)  ::  &  !:
      acrit  = (/ 1.e-06 , 1.e-06 /)    !: minimum fraction for leads in 
      !                                   !  north and south hemisphere
   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/iceini.F90,v 1.6 2006/03/10 10:35:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE ice_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_init  ***
      !!
      !! ** purpose :   
      !!
      !! History :
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and modules
      !!----------------------------------------------------------------------
       CHARACTER(len=80) :: namelist_icename
       
      ! Open the namelist file 
      namelist_icename = 'namelist_ice'
      CALL ctlopn(numnam_ice,namelist_icename,'OLD', 'FORMATTED', 'SEQUENTIAL',   &
                     1,numout,.FALSE.,1)      

      CALL ice_run                    !  read in namelist some run parameters
                 
      ! Louvain la Neuve Ice model
      IF( nacc == 1 ) THEN
          dtsd2   = nfice * rdtmin * 0.5
          rdt_ice = nfice * rdtmin
      ELSE
          dtsd2   = nfice * rdt * 0.5
          rdt_ice = nfice * rdt
      ENDIF

      CALL lim_msh                    ! ice mesh initialization
     
      ! Initial sea-ice state
      IF( .NOT.ln_rstart ) THEN
         numit = 0
         CALL lim_istate              ! start from rest: sea-ice deduced from sst
      ELSE
!!DB
         call rst_ice_read(numit)
      ENDIF
      
      tn_ice(:,:) = sist(:,:)         ! initialisation of ice temperature   
      freeze(:,:) = 1.0 - frld(:,:)   ! initialisation of sea/ice cover    
# if defined key_coupled
      alb_ice(:,:) = albege(:,:)      ! sea-ice albedo
# endif
      
      nstart = numit  + nfice      
      nitrun = nitend - nit000 + 1 
      nlast  = numit  + nitrun 

      IF( nstock == 0  )  nstock = nlast + 1

   END SUBROUTINE ice_init


   SUBROUTINE ice_run
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_run ***
      !!                 
      !! ** Purpose :   Definition some run parameter for ice model
      !!
      !! ** Method  :   Read the namicerun namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicerun
      !!
      !! history :
      !!   2.0  !  03-08 (C. Ethe)  Original code
      !!-------------------------------------------------------------------

      NAMELIST/namicerun/ ln_limdyn, acrit, hsndif, hicdif
      !!-------------------------------------------------------------------

      !                                           ! Read Namelist namicerun 
      REWIND ( numnam_ice )
      READ   ( numnam_ice , namicerun )

      IF( lk_cfg_1d  )  ln_limdyn = .FALSE.       ! No ice transport in 1D configuration

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_run : ice share parameters for dynamics/advection/thermo of sea-ice'
         WRITE(numout,*) ' ~~~~~~'
         WRITE(numout,*) '   switch for ice dynamics (1) or not (0)      ln_limdyn   = ', ln_limdyn
         WRITE(numout,*) '   minimum fraction for leads in the NH (SH)  acrit(1/2)   = ', acrit(:)
         WRITE(numout,*) '   computation of temp. in snow (=0) or not (=9999) hsndif = ', hsndif
         WRITE(numout,*) '   computation of temp. in ice  (=0) or not (=9999) hicdif = ', hicdif
      ENDIF
   END SUBROUTINE ice_run

#else
   !!----------------------------------------------------------------------
   !!   Default option :        Empty module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ice_init        ! Empty routine
   END SUBROUTINE ice_init
#endif

   !!======================================================================
END MODULE iceini
