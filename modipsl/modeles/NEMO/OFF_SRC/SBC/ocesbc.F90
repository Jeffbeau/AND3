MODULE ocesbc
   !!======================================================================
   !!                     ***  MODULE  ocesbc  ***
   !!                     Ocean surface boundary conditions
   !!======================================================================
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/SBC/ocesbc.F90,v 1.2 2005/11/16 16:14:39 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

   !! * Modules used
   USE oce
   USE lib_mpp
   USE in_out_manager ! I/O manager


   IMPLICIT NONE
   PRIVATE


   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      qt  ,         &  !: total surface heat flux (w/m2)
      qsr ,         &  !: solar radiation (w/m2)
      emp ,         &  !: evaporation minus precipitation (kg/m2/s = mm/s)
      emps             !: evaporation - precipitation (free surface)


   !!======================================================================
END MODULE ocesbc
