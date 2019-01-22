!!DB 2009.08.24 -- eliminate non BGCM code options 
MODULE par_trc_trp
   !!======================================================================
   !!                        ***  par_trc_trp  ***
   !! passive tracers :   set the number of passive tracers
   !!======================================================================
   !! History :
   !!   9.0  !  04-03  (C. Ethe)  Orignal
   !!----------------------------------------------------------------------
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/par_trc_trp.F90,v 1.7 2006/04/10 15:40:28 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_passivetrc
   !!---------------------------------------------------------------------
   !!   'key_passivetrc'   :                               Passive tracer
   !!---------------------------------------------------------------------

   IMPLICIT NONE
   PUBLIC
   
   !! jptra   : number of passive tracers
   !! jpdia2d : additional 2d output
   !! jpdia3d : additional 3d output

#if defined key_BGCM_01
   !!---------------------------------------------------------------------
   !!   'key_BGCM_01     simple NPZ model 
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER :: jptra   = 3  !! at this time
#if defined key_trc_diaadd    
   INTEGER, PUBLIC, PARAMETER :: jpdia2d = 1
   INTEGER, PUBLIC, PARAMETER :: jpdia3d = 1
#endif

#elif defined key_PASSIVE_TRACER_01
   !!---------------------------------------------------------------------
   !!   multiple passive tracer module
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER :: jptra   = 5  !! arbitrary test value
!   INTEGER, PUBLIC, PARAMETER :: jptra   = 15  !! large test value (OK)
#if defined key_trc_diaadd    
   INTEGER, PUBLIC, PARAMETER :: jpdia2d = 1
   INTEGER, PUBLIC, PARAMETER :: jpdia3d = 1
#endif

#elif defined key_BGCM_02
   !!---------------------------------------------------------------------
   !!  IML model
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER :: jptra   = 10  !! guess value
#if defined key_trc_diaadd    
   INTEGER, PUBLIC, PARAMETER :: jpdia2d = 1
   INTEGER, PUBLIC, PARAMETER :: jpdia3d = 1
#endif


#else
   !!---------------------------------------------------------------------
   !!   'default'   :          temperature and salinity as passive tracers
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER :: jptra   = 2
#if defined key_trc_diaadd
   INTEGER, PUBLIC, PARAMETER :: jpdia2d = 1
   INTEGER, PUBLIC, PARAMETER :: jpdia3d = 1
#endif
#endif

#else
   !!======================================================================
   !!  Empty module : No passive tracer 
   !!======================================================================
#endif

END MODULE par_trc_trp
