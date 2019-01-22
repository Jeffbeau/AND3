PROGRAM model
   !!----------------------------------------------------------------------
   !!                     ***  PROGRAM model  ***
   !!
   !! ** Purpose :   encapsulate the opa model so that opa can also be 
   !!      called together with the adjoint and linear tangent models
   !!
   !! History :
   !!   8.0  !  01-02  (M. Imbard, A. Weaver)  Original code
   !!   9.0  !  03-10  (G. Madec) F90
   !!----------------------------------------------------------------------
   USE opa                  ! OPA system     (opa_model routine)
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/model.f90,v 1.2 2005/03/27 18:34:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
 
   CALL opa_model           ! OPA system
 
END PROGRAM model
