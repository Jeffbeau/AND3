!!DB: 2008.10.01 -- Routine modified 
!!-1- deleted reference to orca
!!-2- replaced IOIPSL with lib_ncdf 


MODULE dtatem
   !!======================================================================
   !!                     ***  MODULE  dtatem  ***
   !! Ocean data  :  read ocean temperature data from monthly atlas data
   !!=====================================================================
#if defined key_dtatem   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dtatem'                              3D temperature data field
   !!----------------------------------------------------------------------
   !!   dta_tem      : read ocean temperature data
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE daymod          ! calendar

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dta_tem   ! called by step.F90 and inidta.F90

   !! * Shared module variables
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtatem = .TRUE.   !: temperature data flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
      t_dta             !: temperature data at given time-step

   !! * Module variables
   CHARACTER (len=45) ::   &
      cl_tdata
   INTEGER ::   &
      nlecte =  0,   &  ! switch for the first read
      ntem1      ,   &  ! first record used
      ntem2             ! second record used
   REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
      temdta            ! temperature data at two consecutive times

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DTA/dtatem.F90,v 1.12 2006/04/19 14:43:15 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !!   Default case                                            NetCDF file
   !!----------------------------------------------------------------------

   SUBROUTINE dta_tem( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_tem  ***
      !!                    
      !! ** Purpose :   Reads monthly temperature data 
      !! 
      !! ** Method  :   Read on unit numtdt the interpolated temperature 
      !!      onto the model grid.
      !!      Data begin at january. 
      !!      The value is centered at the middle of month.
      !!      In the opa model, kt=1 agree with january 1.
      !!      At each time step, a linear interpolation is applied between 
      !!      two monthly values.
      !!      Read on unit numtdt
      !!
      !! ** Action  :   define t_dta array at time-step kt
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!        !  92-07  (M. Imbard)
      !!        !  99-10  (M.A. Foujols, M. Imbard)  NetCDF FORMAT 
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
!!DB
      USE lib_ncdf

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step

      !! * Local declarations
      INTEGER, PARAMETER ::   &
         jpmois = 12                    ! number of month
      INTEGER ::   ji, jj, jl           ! dummy loop indicies
      INTEGER ::   &
         imois, iman, itime, ik ,    &  ! temporary integers
         i15, ipi, ipj, ipk             !    "          "
#  if defined key_tradmp
      INTEGER ::   &
         il0, il1, ii0, ii1, ij0, ij1   ! temporary integers
# endif

      INTEGER, DIMENSION(jpmois) ::   istep
      REAL(wp) ::   zxy, zl, zdate0
      REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
      REAL(wp), DIMENSION(jpk) ::   zlev
!!DB
      INTEGER :: len, status 
      !!----------------------------------------------------------------------

      ! 0. Initialization
      ! -----------------

      iman  = jpmois
      i15   = nday / 16
      imois = nmonth + i15 - 1
      IF( imois == 0 )   imois = iman

      itime = jpmois
      ipi = jpiglo
      ipj = jpjglo
      ipk = jpk

      ! 1. First call kt=nit000
      ! -----------------------

      IF( kt == nit000 .AND. nlecte == 0 ) THEN
         ntem1 = 0
         cl_tdata = 'data_1m_potential_temperature_nomask.nc'
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' dtatem : Monthly Temperature fields'
         IF(lwp) WRITE(numout,*) ' ~~~~~~'
         IF(lwp) WRITE(numout,*) '             NetCDF File'
         IF(lwp) WRITE(numout,*)cl_tdata
         IF(lwp) WRITE(numout,*)
         
#if defined key_agrif
         if ( .NOT. Agrif_Root() ) then
            cl_tdata = TRIM(Agrif_CFixed())//'_'//TRIM(cl_tdata)
         endif
#endif         
!!DB: 
         call ncdf_get_dim_size(cl_tdata, 'time_counter', itime, status)
         call ncdf_get_dim_size(cl_tdata, 'x', ipi, status)
         call ncdf_get_dim_size(cl_tdata, 'y', ipj, status)
         call ncdf_get_dim_size(cl_tdata, 'z', ipk, status)
         if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with dimensions'
               write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               write(numout,*) ' ipk ',ipk,' jpk ',jpk
            endif
            stop 'dtatem'
         endif
!!DB: 
         if( itime /= jpmois ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with time coordinates'
               write(numout,*) ' itime ',itime,' jpmois ',jpmois
            endif
            stop 'dta_sal'
         endif

      ENDIF


      ! 2. Read monthly file
      ! -------------------

      IF( ( kt == nit000 .AND. nlecte == 0 ) .OR. imois /= ntem1 ) THEN
         nlecte = 1

         ! Calendar computation
         
         ntem1 = imois        ! first file record used 
         ntem2 = ntem1 + 1    ! last  file record used
         ntem1 = MOD( ntem1, iman )
         IF( ntem1 == 0 )   ntem1 = iman
         ntem2 = MOD( ntem2, iman )
         IF( ntem2 == 0 )   ntem2 = iman
         if(lwp) write(numout,*) 'dtatem reading records ',ntem1, ntem2
         
         ! Read monthly temperature data
!!DB
         call ncdf_read(cl_tdata,'votemper',temdta(:,:,:,1),-ntem1,status)
         call ncdf_read(cl_tdata,'votemper',temdta(:,:,:,2),-ntem2,status)


         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' Monthly temperature input OK'
         IF(lwp) WRITE(numout,*)
         
#if defined key_tradmp
!!DB -- orca-related code deleted
#endif

         !                                  ! Mask
         DO jl = 1, 2
            temdta(:,:,:,jl) = temdta(:,:,:,jl) * tmask(:,:,:)
            temdta(:,:,jpk,jl) = 0.
            IF( lk_zps ) THEN                ! z-coord. with partial steps
               DO jj = 1, jpj                  ! interpolation of temperature at the last level
                  DO ji = 1, jpi
                     ik = mbathy(ji,jj) - 1
                     IF( ik > 2 ) THEN
                        zl = ( gdept(ik) - fsdept(ji,jj,ik) ) / ( gdept(ik) - gdept(ik-1) )
                        temdta(ji,jj,ik,jl) = (1.-zl) * temdta(ji,jj,ik,jl) + zl * temdta(ji,jj,ik-1,jl) 
                     ENDIF
                  END DO
               END DO
            ENDIF
         END DO

!         IF(lwp) THEN
!            WRITE(numout,*) ' temperature month ', ntem1, ntem2
!            WRITE(numout,*)
!            WRITE(numout,*) '  month = ', ntem1, '  level = 1'
!            CALL prihre( temdta(:,:,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
!            WRITE(numout,*) '  month = ', ntem1, '  level = ', jpk/2
!            CALL prihre( temdta(:,:,jpk/2,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
!            WRITE(numout,*) '  month = ',ntem1,'  level = ', jpkm1
!            CALL prihre( temdta(:,:,jpkm1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
!         ENDIF
      ENDIF

 
      ! 2. At every time step compute temperature data
      ! ----------------------------------------------
      zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!!DB: Note that for monthly data it is not necessary to do this every kt.
!!    So in the future use a code fragment like: 
!      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
!         t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)
!      endif

      t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)


   END SUBROUTINE dta_tem

#else
   !!----------------------------------------------------------------------
   !!   Default case                           NO 3D temperature data field
   !!----------------------------------------------------------------------
   USE in_out_manager
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtatem = .FALSE.   !: temperature data flag
CONTAINS
   SUBROUTINE dta_tem( kt )        ! Empty routine
      if(lwp) WRITE(numout,*) 'dta_tem: You should not have seen this print! error?', kt
   END SUBROUTINE dta_tem
#endif
   !!======================================================================
END MODULE dtatem
