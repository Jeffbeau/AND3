MODULE par_ice
   !!======================================================================
   !!                       ***  MODULE par_ice   ***
   !! Sea-Ice model : definition of the parameters
   !!======================================================================
   !!----------------------------------------------------------------------
   !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/par_ice.F90,v 1.4 2005/03/27 18:34:42 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce

   IMPLICIT NONE
   PUBLIC               ! allows par_oce and par_kind to be known in ice modules

   INTEGER, PUBLIC, PARAMETER ::   &  !:
      jpkmax =  1    ,      &  !: ???
      jpsmax =  2              !: ???

   INTEGER, PUBLIC, PARAMETER ::   &  !: 
      jplayers   = 2 ,           &  !: number of vertical ice layers
      jplayersp1 = jplayers + 1     !: ???

   !!======================================================================
END MODULE par_ice
