MODULE trcini
   !!==========================================================================
   !!                       *** MODULE trcini ***  
   !! Ocean passive tracers:  Manage the passive tracer initialization 
   !!=========================================================================   
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_ini : Initialization for passive tracer
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trcini.F90,v 1.5 2006/04/10 15:40:28 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE ioipsl
   USE oce_trc
   USE trc
   USE sms
   USE lib_mpp
   USE lbclnk

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_ini

#if defined key_trc_lobster1
   !!----------------------------------------------------------------------
   !!   'key_trc_lobster1'                        LOBSTER1 biological model  
   !!----------------------------------------------------------------------
#  include "trcini.lobster1.h90"

#elif defined key_trc_pisces
   !!----------------------------------------------------------------------
   !!   'key_trc_pisces'                            PISCES biological model                  
   !!----------------------------------------------------------------------
#  include "trcini.pisces.h90"

#elif defined key_cfc
   !!----------------------------------------------------------------------
   !!   'key_cfc  '                                          CFC model                  
   !!----------------------------------------------------------------------
#  include "trcini.cfc.h90"

#else
   !!----------------------------------------------------------------------
   !!   Default option                               
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE trc_ini  ***
      !!              
      !! ** Purpose : Initialization for passive tracer
      !!              for restart or not
      !!
      !! History :
      !!        !  00-04  O. Aumont, M.A. Foujols HAMOCC3 and P3ZD
      !!   8.5  !  05-03  O.Aumont and A.El Moussaoui  F90
      !!   9.0  !  05-10  C. Ethe  Modularity
      !!----------------------------------------------------------------------
      !! * local declarations
      INTEGER ::                   & 
         ji ,jj ,jk ,jn, jl        ! dummy loop indices  
      !!---------------------------------------------------------------------


      !! 1. initialization of passives tracers field
      !! -------------------------------------------
      DO jn = 1, jptra
         trn(:,:,:,jn)=0.e0
         tra(:,:,:,jn)=0.e0
      END DO

#if defined key_trc_diaadd
      !! initialization of output 2d and 3d arrays

      DO jn = 1, jpdia2d
         trc2d(:,:,jn)=0.e0
      END DO

      DO jn = 1, jpdia3d
         trc3d(:,:,:,jn)=0.e0
      END DO
#endif

#if defined key_trc_diabio
      !! initialization of biological trends
      DO jn=1,jpdiabio
         trbio(:,:,:,jn) = 0.e0
      END DO
#endif

#if defined key_trc_diatrd
      !! initialization of tracer trends
      DO jl = 1, jpdiatrc
         DO jn = 1, jptra
            IF (luttrd(jn)) trtrd(:,:,:,ikeep(jn),jl) = 0.e0
         END DO
      END DO
#endif      

      IF( lwp ) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' trcini: generic initialisation done '
         WRITE(numout,*) ' '
      ENDIF

   END SUBROUTINE trc_ini

#endif

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini              ! Empty routine

   END SUBROUTINE trc_ini
#endif

   !!======================================================================
END MODULE trcini
