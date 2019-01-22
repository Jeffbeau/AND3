MODULE traqsr
   !!======================================================================
   !!                       ***  MODULE  traqsr  ***
   !! Ocean physics: solar radiation penetration in the top ocean levels
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   tra_qsr      : trend due to the solar radiation penetration
   !!   tra_qsr_init : solar radiation penetration initialization
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/TRA/traqsr.F90,v 1.2 2005/11/16 16:15:21 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_qsr_init ! routine called by opa.F90

   !! * Shared module variables
   LOGICAL, PUBLIC ::   ln_traqsr = .TRUE.   !: qsr flag (Default=T)

   !! * Module variables
   REAL(wp), PUBLIC ::  & !!! * penetrative solar radiation namelist *
      xsi1 = 0.35_wp     ! first depth of extinction 

CONTAINS

   SUBROUTINE tra_qsr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr_init  ***
      !!
      !! ** Purpose :   Initialization for the penetrative solar radiation
      !!
      !! ** Method  :   The profile of solar radiation within the ocean is set
      !!      from two length scale of penetration (xsr1,xsr2) and a ratio
      !!      (rabs). These parameters are read in the namqsr namelist. The
      !!      default values correspond to clear water (type I in Jerlov' 
      !!      (1968) classification.
      !!         called by tra_qsr at the first timestep (nit000)
      !!
      !! ** Action  : - initialize xsr1, xsr2 and rabs
      !!
      !! Reference :
      !!   Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!
      !! History :
      !!   8.5  !  02-06  (G. Madec) Original code
      !!----------------------------------------------------------------------
      !! * Local declarations

      NAMELIST/namqsr/ ln_traqsr, xsi1
      !!----------------------------------------------------------------------

      ! Read Namelist namqsr : ratio and length of penetration
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, namqsr )

      ! Parameter control and print
      ! ---------------------------
      IF( ln_traqsr  ) THEN
        IF ( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_qsr_init : penetration of the surface solar radiation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '    Namelist namqsr : set the parameter of penetration'
         WRITE(numout,*) '        first depth of extinction        xsi1        = ',xsi1
         WRITE(numout,*) ' '
        END IF
      ELSE
        IF ( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_qsr_init : NO solar flux penetration'
         WRITE(numout,*) '~~~~~~~~~~~~'
        END IF
      ENDIF

      IF( xsi1 < 0.e0 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '             0<xsi1 not satisfied'
         nstop = nstop + 1
      ENDIF


   END SUBROUTINE tra_qsr_init

   !!======================================================================
END MODULE traqsr
