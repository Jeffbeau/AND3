MODULE cpl_oce
   !!======================================================================
   !!                   ***  MODULE  cpl_oce  ***
   !! Ocean coupling:  ocean-atmosphere-sea ice coupled exchanges 
   !!=====================================================================
#if defined key_coupled
   !!----------------------------------------------------------------------
   !!   'key_coupled'                              Coupled Ocean/Atmosphere
   !!----------------------------------------------------------------------
   !! ** Purpose :   Atmosphere/Ice/Ocean Coupling
   !!      GASTON TEAM (CERFACS, Meteo-France, LMD, LSCE, IPSL, LODYC)
   !!
   !! history :
   !!  8.0   ! 08-98  (M.A. Foujols, M. Imbard)  Original code
   !!  8.5   ! 06/02  (G. Madec)  modules
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/cpl_oce.F90,v 1.5 2005/03/27 18:34:46 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce          ! ocean parameters

   IMPLICIT NONE

   !! ---------------------------------------------------------------------
   !! Ocean/Ice/Atmosphere Coupling
   !! -----------------------------
   LOGICAL, PUBLIC, PARAMETER :: lk_cpl = .TRUE.   !: coupled flag

   INTEGER, PARAMETER ::   &  !:  
      jpmaxfld = 40           !: Number of maximum fields exchange betwwen
      !                       ! the ocean and the coupler

   !!---------------------------------------------------------------------
   !! SIPC Method   (L. Terray, S. Valcke, CERFACS)
   !! -----------

   INTEGER, PARAMETER ::   &  !:
      jpbyteint = 4,       &  !: number of bytes per integer
      jpbyterea = 8,       &  !: number of bytes per real
      jpbytecha = 1           !: number of bytes per character  

   INTEGER, PARAMETER ::   &  !:
      jptest = 100            !: The models will test during 2*jptest 
      !                       ! seconds if the file DUMMY_SIPC has been
      !                       ! created by OASIS, signaling that the
      !                       ! SHM pools are opened. After, it aborts.

   !!---------------------------------------------------------------------
   !! PIPE Method   (L. Terray, CERFACS)
   !! -----------

   INTEGER, PARAMETER ::   &  !:
      jpread = 0,          &  !: 
      jpwrit = 1              !: 

   !!---------------------------------------------------------------------
   !! Messag Passing Method (CLIM)
   !! ----------------------------
!!!INCLUDE '../../CPL/include/clim.h90'
!!
!! -- clim.h   18-08-95   Version 2.0   Author: Laurent Terray
!!    ******
!!             26-10-99   Version 2.4   Jean Latour (F.S.E.) MPI-2 support
!!
!!   clim.h90  13-08-04  Change to F90 C. Levy
!!@
!!@  Contents : variables related to the CLIM library
!!@  --------
!!@ For complete definition, see the CLIM manual
!!@
      INTEGER (kind=4)	CLIM_MaxMod,    CLIM_MaxPort,  CLIM_MaxSegments, &
               CLIM_MaxTag, &
               CLIM_MaxLink, &
               CLIM_ParSize, & 
               CLIM_Clength, &
               CLIM_MaxCodes
!!
      INTEGER (kind=4) CLIM_Void
!!
      INTEGER (kind=4) CLIM_In,	CLIM_Out,	CLIM_InOut
!!
      INTEGER (kind=4) CLIM_Strategy,  CLIM_Segments,  &
               CLIM_Serial,    CLIM_Length,    CLIM_Orange, &
               CLIM_Apple,     CLIM_Offset, &
               CLIM_Box,	CLIM_SizeX,	CLIM_SizeY, &
               CLIM_LdX
!!
      INTEGER  (kind=4)CLIM_Integer,	CLIM_Real,	CLIM_Double
!!
      INTEGER  (kind=4)CLIM_StopPvm,   CLIM_ContPvm
!!
      INTEGER (kind=4)	CLIM_Ok
      INTEGER (kind=4) CLIM_FastExit, 	CLIM_BadName, 	CLIM_BadPort, &
               CLIM_BadType, 	CLIM_DoubleDef, CLIM_NotStep, &
               CLIM_IncStep, 	CLIM_IncSize, 	CLIM_NotClim, &
               CLIM_TimeOut, &
               CLIM_Pvm, 	CLIM_FirstCall, CLIM_PbRoute, &
               CLIM_Group, 	CLIM_BadTaskId, CLIM_NoTask, &
               CLIM_InitBuff, 	CLIM_Pack, 	CLIM_Unpack, &
               CLIM_Down, 	CLIM_PvmExit
!!
      INTEGER (kind=4) CLIM_jpmax, 	CLIM_jpmx8, 	CLIM_Mpi
!!
!!-----Parameter sizes
!!
      PARAMETER ( CLIM_Void    = 0  )
      PARAMETER ( CLIM_MaxMod  = 8 )
      PARAMETER ( CLIM_MaxPort = 50 )
      PARAMETER ( CLIM_MaxSegments = 50 )
      PARAMETER ( CLIM_MaxLink = CLIM_MaxMod * CLIM_MaxPort )
      PARAMETER ( CLIM_ParSize = 2*CLIM_MaxSegments+2 )
      PARAMETER ( CLIM_MaxTag  = 16777215 )
      PARAMETER ( CLIM_Clength = 32 )
!!
!!-----Dimension of buffer for packing / unpacking messages with MPI
!!     (must be equal to jpmax of Oasis)
!!
      PARAMETER ( CLIM_jpmax = 400000 )
      PARAMETER ( CLIM_jpmx8 = CLIM_jpmax*8 )
!!
!!-----Ports status
!!
      PARAMETER ( CLIM_In      = 1 )
      PARAMETER ( CLIM_Out     = 0 )
      PARAMETER ( CLIM_InOut   = 2 )
!!
!!-----Parallel distribution
!!
      PARAMETER ( CLIM_Strategy = 1 )
      PARAMETER ( CLIM_Segments = 2 )
      PARAMETER ( CLIM_Serial   = 0 )
      PARAMETER ( CLIM_Apple    = 1 )
      PARAMETER ( CLIM_Box      = 2 )
      PARAMETER ( CLIM_Orange   = 3 )
      PARAMETER ( CLIM_Offset   = 2 )
      PARAMETER ( CLIM_Length   = 3 )
      PARAMETER ( CLIM_SizeX    = 3 )
      PARAMETER ( CLIM_SizeY    = 4 )
      PARAMETER ( CLIM_LdX      = 5 )
!!
!!-----Datatypes
!!
      PARAMETER ( CLIM_Integer = 1 )
      PARAMETER ( CLIM_Real    = 4 ) 
      PARAMETER ( CLIM_Double  = 8 )
!!
!!-----Quit parameters
!!
      PARAMETER ( CLIM_ContPvm = 0 )
      PARAMETER ( CLIM_StopPvm = 1 )
!!
!!-----Error Codes
!!
      PARAMETER ( CLIM_MaxCodes  = -22 )
!!
      PARAMETER ( CLIM_Ok	 = 0 )
      PARAMETER ( CLIM_FastExit  = -1 )
      PARAMETER ( CLIM_BadName   = -2 )
      PARAMETER ( CLIM_BadPort   = -3 )
      PARAMETER ( CLIM_BadType   = -4 )
      PARAMETER ( CLIM_DoubleDef = -5 )
      PARAMETER ( CLIM_NotStep   = -6 )
      PARAMETER ( CLIM_IncStep   = -7 )
      PARAMETER ( CLIM_IncSize   = -8 )
      PARAMETER ( CLIM_NotClim   = -9 )
      PARAMETER ( CLIM_TimeOut   = -10 )
      PARAMETER ( CLIM_Pvm       = -11 )
      PARAMETER ( CLIM_FirstCall = -12 )
      PARAMETER ( CLIM_PbRoute   = -13 )
      PARAMETER	( CLIM_Group     = -14 )
      PARAMETER ( CLIM_BadTaskId = -15 )
      PARAMETER ( CLIM_NoTask    = -16 )
      PARAMETER ( CLIM_InitBuff  = -17 )
      PARAMETER ( CLIM_Pack      = -18 )
      PARAMETER ( CLIM_Unpack    = -19 )
      PARAMETER ( CLIM_Down      = -20 )
      PARAMETER ( CLIM_PvmExit   = -21 )
      PARAMETER ( CLIM_Mpi       = -22 )

!!
!     --- end of clim.h90
!!!END-----------------------------------------------------------------

!!!INCLUDE '../../CPL/include/mpiclim.h90'
!!
!! -- mpiclim.h  26-10-99   Version 2.4   Author: Jean Latour (F.S.E.)
!!    *********
!!    mpiclim.h90 13-08-04 change to F90 C. Levy
!!@
!!@  Contents : variables related to MPI-2 message passing
!!@  --------
!!@
!!@ -- mpi_totproc: number of processors on which to launch each model
!!@
!!@ -- mpi_nproc: number of processors involved in the coupling for
!!@               each model
!!@ -- cmpi_modnam: models name
!!     -----------------------------------------------------------------
!!
      INTEGER (kind=4) mpi_totproc(1:CLIM_MaxMod-1),mpi_nproc(0:CLIM_MaxMod-1)
!!
      CHARACTER (len=6) cmpi_modnam(1:CLIM_MaxMod-1)
!!
      common/CLIM_mpiclim/mpi_totproc, mpi_nproc, cmpi_modnam 
!!
!!!END-----------------------------------------------------------------


   !!----------------------------------------------------------------------
   !!  Atmosphere/Ice/Ocean Coupling
   !!---------------------------------

   REAL(wp), DIMENSION(jpi,jpj) ::   &   !: data from an atmospheric model
      qc  ,       &  !: total surf. total heat flux (wm-2)
      ec  ,       &  !: surface water flux (kg m-2s-1)
      qsrc           !: solar radiation (w m-2)

#  if defined key_ice_lim
   REAL(wp), DIMENSION(jpi,jpj) ::   &  !:
      watm        ,    &  !:
      tatm        ,    &  !:
      hatm        ,    &  !:
      vatm        ,    &  !:
      catm                !:
#  endif

   !! Coupling

   INTEGER ::          &  !:
      npioc       ,    &  !: process-id of ocean PROGRAM
      nexco       ,    &  !: exchange frequency for fluxes
      nmodcpl     ,    &  !: coupling mode
      nflxc2o     ,    &  !: fluxes field number coupler to ocean
      ntauc2o     ,    &  !: stress field number coupler to ocean
      nfldo2c             !: surface field number ocean to coupler

   CHARACTER(len=4) ::   cchan       !: type of message passing (pipe or clim) 
   CHARACTER(len=6) ::   cplmodnam   !: model name
   CHARACTER(len=5) ::   cploasis    !: coupler name

   CHARACTER(len=8), DIMENSION(jpmaxfld) ::   &  !:
      cpl_f_readflx,   &  !: coupler to ocean file name for flx.coupled
      cpl_f_readtau,   &  !: coupler to ocean file name for tau.coupled
      cpl_f_writ   ,   &  !: ocean to coupler file name for cpl_stp
      cpl_readflx  ,   &  !: coupler to ocean field name for flx.coupled
      cpl_readtau  ,   &  !: coupler to ocean field name for tau.coupled
      cpl_writ            !: ocean to coupler field name for cpl_stp

  REAL(wp), DIMENSION(jpi,jpj) ::   &  !:
      sstoc,     &  !: work array to average sst
      sieoc,     &  !: work array to average Ice index
      alboc,     &  !: work array to average Ice Albedo
      ticoc         !: work array to average Ice temperature


   !! -- inc_sipc.h   97-08-11   Version 2.0   Author: S&A
   !!    **********
   !! variables describing pools formed of shared memory segments 

   INTEGER ::   &  !:
      mpoolinitr,  &  !: handles associated to model pools for passing
      mpoolinitw      !: initial info (r=read, w=write)

   INTEGER, DIMENSION(jpmaxfld) ::   &  !:
      mpoolwrit,   &  !: handles associated to pools used to pass fields
      !               !  exchanged from model to coupler 
      !               !  (see libsipc/SIPC_Write_Model.f)
      mpoolread       !: handles associated to pools used to pass fields
      !               !  exchanged from model to coupler
      !               !  (see libsipc/SIPC_Read_Model.f)

#else
   !!----------------------------------------------------------------------
   !!   Default case                                Forced Ocean/Atmosphere
   !!----------------------------------------------------------------------
   !!   Empty module
   LOGICAL, PUBLIC, PARAMETER :: lk_cpl = .FALSE.   !: coupled flag
#endif

   !!----------------------------------------------------------------------
END MODULE cpl_oce
