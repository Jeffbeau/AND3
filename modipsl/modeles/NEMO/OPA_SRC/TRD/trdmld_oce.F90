MODULE trdmld_oce
   !!======================================================================
   !!                   ***  MODULE trdmld_oce  ***
   !! Ocean trends :   set tracer and momentum trend variables
   !!======================================================================
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmld_oce.F90,v 1.2 2005/03/27 18:35:23 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! ocean parameters

   IMPLICIT NONE
   PUBLIC

   INTEGER, PARAMETER ::            &  !: mixed layer trends index
      jpmldxad = 1,   &  !: zonal advection
      jpmldyad = 2,   &  !: meridionnal advection
      jpmldzad = 3,   &  !: vertical advection
      jpmldldf = 4,   &  !: lateral diffusion (horiz. component+Beckman)
      jpmldfor = 5,   &  !: forcing 
      jpmldevd = 6,   &  !: entrainment due to vertical diffusion (TKE)
      jpmldzdf = 7,   &  !: explicit vertical part if isopycnal diffusion
      jpmldxei = 8,   &  !: eddy induced zonal advection
      jpmldyei = 9,   &  !: eddy induced meridional advection
      jpmldzei =10       !: eddy induced vertical advection

#if   defined  key_trdmld   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_trdmld'                         mixed layer trends diagnostics
   !!----------------------------------------------------------------------

   !! Trends diagnostics parameters
   !!---------------------------------------------------------------------
   INTEGER, PARAMETER ::            &  !:
# if defined key_traldf_eiv
      jpltrd = 10,  &  !: number of mixed-layer trends arrays
      jpktrd = jpk     !: max level for mixed-layer trends diag.
# else
      jpltrd = 7,   &  !: number of mixed-layer trends arrays
      jpktrd = jpk     !: max level for mixed-layer trends diag.
# endif

#endif
  !!======================================================================
END MODULE trdmld_oce
