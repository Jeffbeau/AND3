MODULE trdmod_oce
   !!======================================================================
   !!                   ***  MODULE trdmod_oce  ***
   !! Ocean trends :   set tracer and momentum trend variables
   !!======================================================================
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmod_oce.F90,v 1.2 2005/03/27 18:35:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE trdicp_oce              ! ocean momentum/tracers bassin properties trends variables
   USE trdmld_oce              ! ocean active mixed layer tracers trends variables
   USE trdvor_oce              ! ocean vorticity trends variables

   !! Control parameters
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC ::   l_trdtra = .FALSE.    !: tracers  trend flag
   LOGICAL, PUBLIC ::   l_trddyn = .FALSE.    !: momentum trend flag

  !!======================================================================
END MODULE trdmod_oce
