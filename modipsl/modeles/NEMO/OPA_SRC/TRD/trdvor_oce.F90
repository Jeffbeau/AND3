MODULE trdvor_oce
   !!======================================================================
   !!                   ***  MODULE trdvor_oce  ***
   !! Ocean trends :   set vorticity trend variables
   !!======================================================================
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdvor_oce.F90,v 1.2 2005/03/27 18:35:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce      ! ocean parameters

   IMPLICIT NONE
   PUBLIC

   INTEGER,PARAMETER :: jplvor = 11     ! Number of vorticity trend terms

   INTEGER, PARAMETER ::            &  !: vorticity trends index
      jpvorprg = 1,   &  !: Pressure Gradient Trend
      jpvorkeg = 2,   &  !: KE Gradient Trend
      jpvorrvo = 3,   &  !: Relative Vorticity Trend
      jpvorpvo = 4,   &  !: Planetary Vorticity Term Trend
      jpvorldf = 5,   &  !: Horizontal Diffusion Trend
      jpvorzad = 6,   &  !: Vertical Advection Trend
      jpvorzdf = 7,   &  !: Vertical Diffusion Trend
      jpvorspg = 8,   &  !: Surface Pressure Grad. Trend
      jpvorbev = 9,   &  !: Beta V
      jpvorswf =10,   &  !: wind stress forcing term
      jpvorbfr =11       !: bottom friction term

  !!======================================================================
END MODULE trdvor_oce
