MODULE domstp
   !!==============================================================================
   !!                       ***  MODULE domstp   ***
   !! Ocean initialization : time domain
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dom_stp        : ocean time domain initialization
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * routine accessibility
   PUBLIC dom_stp        ! routine called by inidom.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/DOM/domstp.F90,v 1.2 2005/11/16 16:12:12 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dom_stp
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dom_stp  ***
      !!          
      !! ** Purpose :   Intialize ocean time step for the run
      !!
      !! ** Method  : - Initialization of a coef. use in the Asselin time
      !!      filter:  atfp1 = 1 - 2 * atfp  where atfp is the Asselin time
      !!      filter parameter read in namelist
      !!              - Model time step:
      !!      nacc = 0 : synchronous time intergration. 
      !!      There is one time step only, defined by: rdt, rdttra(k)=rdt
      !!      nacc = 1 : accelerating the convergence. There is 2 different
      !!      time steps for dynamics and tracers:
      !!        rdt      : dynamical part
      !!        rdttra(k): temperature and salinity
      !!      The tracer time step is a function of vertical level. the model
      !!      reference time step ( i.e. for wind stress, surface heat and
      !!      salt fluxes) is the surface tracer time step is rdttra(1).
      !!         N.B. depth dependent acceleration of convergence is not im-
      !!      plemented for s-coordinate.
      !!
      !! ** Action  : - rdttra   : vertical profile of tracer time step
      !!              - atfp1    : = 1 - 2*atfp
      !!
      !! References :
      !!      Bryan, K., 1984, J. Phys. Oceanogr., 14, 666-673.
      !!
      !! History :
      !!        !  90-10  (O. Marti)  Original code
      !!        !  96-01  (G. Madec)  terrain following coordinates
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   jk              ! dummy loop indice
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_stp : time stepping setting'
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! 0. Asselin Time filter
      ! ----------------------
      
      atfp1 = 1. - 2. * atfp


      SELECT CASE ( nacc )

         CASE ( 0 )                ! Synchronous time stepping
            IF(lwp) WRITE(numout,*)'               synchronous time stepping'
            IF(lwp) WRITE(numout,*)'               dynamics and tracer time step = ', rdt/3600., ' hours'

            rdttra(:) = rdt

         CASE ( 1 )                ! Accelerating the convergence
            IF(lwp) WRITE(numout,*) '              no tracer damping in the turbocline'
            IF(lwp) WRITE(numout,*)'               accelerating the convergence'
            IF(lwp) WRITE(numout,*)'               dynamics time step = ', rdt/3600., ' hours'
#if defined key_s_coord
            IF( rdtmin /= rdtmax ) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*)' depth dependent acceleration of &
                                      &convergence not implemented in s-coordinates'
               nstop = nstop + 1
            ENDIF
#endif
#if defined key_partial_steps
            IF( rdtmin /= rdtmax ) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*)' depth dependent acceleration of &
                                      &convergence not implemented for partial steps case'
               nstop = nstop + 1
            ENDIF
#endif
            IF(lwp) WRITE(numout,*)'         tracers   time step :  dt (hours)  level'

            DO jk = 1, jpk
               IF( fsdept(1,1,jk) <= rdth ) rdttra(jk) = rdtmin
               IF( fsdept(1,1,jk) >  rdth ) THEN
                  rdttra(jk) = rdtmin + ( rdtmax - rdtmin )   &
                                      * ( EXP( ( fsdept(1,1,jk ) - rdth ) / rdth ) - 1. )   &
                                      / ( EXP( ( fsdept(1,1,jpk) - rdth ) / rdth ) - 1. )
               ENDIF
               IF(lwp) WRITE(numout,9200) rdttra(jk)/3600., jk
            END DO  
 9200       FORMAT(36x,f5.2,'     ',i3)

         CASE DEFAULT              ! E R R O R 

            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' nacc value e r r o r, nacc= ',nacc
            IF(lwp) WRITE(numout,*) ' we stop'
            nstop = nstop + 1

      END SELECT

   END SUBROUTINE dom_stp

   !!======================================================================
END MODULE domstp
