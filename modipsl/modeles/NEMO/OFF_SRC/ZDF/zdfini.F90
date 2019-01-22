MODULE zdfini
   !!======================================================================
   !!              ***  MODULE  zdfini  ***
   !! Ocean physics : define vertical mixing variables
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   zdf_init    : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/ZDF/zdfini.F90,v 1.2 2005/11/16 16:16:03 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! mesh and scale factors
   USE zdf_oce         ! TKE vertical mixing          
   USE ldfslp          ! ???

   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC zdf_init          ! routine called by opa.F90
   
CONTAINS

   SUBROUTINE zdf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_init  ***
      !! 
      !! ** Purpose :   initializations of the vertical ocean physics
      !!
      !! ** Method  :   Read namelist namzdf, control cpp keys
      !!
      !! History :
      !!        !  97-06  (G. Madec)  Original code from inimix
      !!   8.5  !  02-08  (G. Madec)  F90 : free form
      !!----------------------------------------------------------------------
      !! * Local declarations

      !! * Namelist
      NAMELIST/namzdf/ avt0, ln_zdfnpc

      ! Read namzdf namelist : vertical mixing parameters
      ! --------------------
      REWIND( numnam )
      READ  ( numnam, namzdf )

      ! Parameter print
      ! ---------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_init: vertical physics'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '          Namelist namzdf : set vertical mixing mixing parameters'
         WRITE(numout,*) '             non-penetrative convection       ln_zdfnpc = ', ln_zdfnpc
         WRITE(numout,*) '             vertical eddy diffusivity           avt0   = ', avt0
      ENDIF

   END SUBROUTINE zdf_init

   !!======================================================================
END MODULE zdfini
