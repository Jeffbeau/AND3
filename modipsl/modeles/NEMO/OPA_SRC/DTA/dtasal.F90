!!DB: 2008.10.01 -- Routine modified 
!!-1- deleted reference to orca
!!-2- replaced IOIPSL with lib_ncdf 

MODULE dtasal
   !!======================================================================
   !!                     ***  MODULE  dtasal  ***
   !! Ocean data  :  read ocean salinity data from monthly atlas data
   !!=====================================================================
#if defined key_dtasal   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dtasal'                                          salinity data
   !!----------------------------------------------------------------------
   !!   dta_sal      : read ocean salinity data
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE daymod          ! calendar

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dta_sal   ! called by step.F90 and inidta.F90
   
   !! * Shared module variables
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtasal = .TRUE.    !: salinity data flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
      s_dta       !: salinity data at given time-step

   !! * Module variables
   CHARACTER (len=32) ::   clname
   INTEGER ::   &
      nlecsa = 0,   &  ! switch for the first read
      nsal1     ,   &  ! first record used
      nsal2            ! second record used
   REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
      saldta    ! salinity data at two consecutive times

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DTA/dtasal.F90,v 1.12 2006/04/19 14:43:14 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !!   Default option:                                         NetCDF file
   !!----------------------------------------------------------------------

   SUBROUTINE dta_sal( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_sal  ***
      !!        
      !! ** Purpose :   Reads monthly salinity data
      !!             
      !! ** Method  : - Read on unit numsdt the monthly salinity data interpo-
      !!     lated onto the model grid.
      !!              - At each time step, a linear interpolation is applied
      !!     between two monthly values.
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!        !  92-07  (M. Imbard)
      !!   9.0  !  02-06  (G. Madec)  F90: Free form and module 
      !!----------------------------------------------------------------------
      !! * Modules used
!!DB
      USE lib_ncdf

      !! * Arguments
      INTEGER, INTENT(in) ::   kt             ! ocean time step

      !! * Local declarations


      INTEGER, PARAMETER ::   jpmois = 12, jpf = 1
      INTEGER ::   ji, jj, jl           ! dummy loop indicies
      INTEGER ::   &
         imois, iman, ik, i15,       &  ! temporary integers
         ipi, ipj, ipk, itime           !    "          "
#if defined key_tradmp
      INTEGER ::   &
         jk, il0, il1,               &  ! temporary integers
         ii0, ii1, ij0, ij1             !    "          "
#endif

      INTEGER, DIMENSION(jpmois) ::   istep
      REAL(wp) ::   &
         zxy, zl, zdate0
      REAL(wp), DIMENSION(jpi,jpj) ::   zlon, zlat
      REAL(wp), DIMENSION(jpk) ::   zlev
!!DB
      INTEGER :: len, status 

      !!----------------------------------------------------------------------

      ! 0. Initialization
      ! -----------------

      iman  = jpmois
      i15   = nday / 16

      imois = nmonth + i15 - 1
      IF( imois == 0 ) imois = iman

      itime = jpmois
      ipi=jpiglo
      ipj=jpjglo
      ipk = jpk

      ! 1. First call kt=nit000
      ! -----------------------

      IF( kt == nit000 .AND. nlecsa == 0 ) THEN
         nsal1 = 0
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' dta_sal : monthly salinity data in NetCDF file'
            WRITE(numout,*) ' ~~~~~~~'
            WRITE(numout,*)
         ENDIF

         ! open file
         clname = 'data_1m_salinity_nomask.nc'
#if defined key_agrif
         if ( .NOT. Agrif_Root() ) then
            clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
         endif
#endif          
!!DB: 
         call ncdf_get_dim_size(clname, 'time_counter', itime, status)
         call ncdf_get_dim_size(clname, 'x', ipi, status)
         call ncdf_get_dim_size(clname, 'y', ipj, status)
         call ncdf_get_dim_size(clname, 'z', ipk, status)
         if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with dimensions'
               write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               write(numout,*) ' ipk ',ipk,' jpk ',jpk
            endif
            stop 'dtasal'
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

      IF( ( kt == nit000 .AND. nlecsa == 0) .OR. imois /= nsal1 ) THEN
         nlecsa = 1
         
         ! 2.1 Calendar computation
         
         nsal1 = imois        ! first file record used 
         nsal2 = nsal1 + 1    ! last  file record used
         nsal1 = MOD( nsal1, iman )
         IF( nsal1 == 0 ) nsal1 = iman
         nsal2 = MOD( nsal2, iman )
         IF( nsal2 == 0 ) nsal2 = iman
         if(lwp) write(numout,*) 'dtasal reading records ',nsal1, nsal2
         
         ! 2.3 Read monthly salinity data

!!DB
         call ncdf_read(clname,'vosaline',saldta(:,:,:,1),-nsal1,status)
         call ncdf_read(clname,'vosaline',saldta(:,:,:,2),-nsal2,status)

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' Monthly salinity input OK'
         IF(lwp) WRITE(numout,*)
         
#if defined key_tradmp
!!DB -- orca code deleted
#endif
         
         !                                     ! Mask
         DO jl = 1, 2
            saldta(:,:,:,jl) = saldta(:,:,:,jl)*tmask(:,:,:)
            saldta(:,:,jpk,jl) = 0.
            IF( lk_zps ) THEN          ! z-coord. partial steps
               DO jj = 1, jpj          ! interpolation of salinity at the last ocean level (i.e. the partial step)
                  DO ji = 1, jpi
                     ik = mbathy(ji,jj) - 1
                     IF( ik > 2 ) THEN
                        zl = ( gdept(ik) - fsdept(ji,jj,ik) ) / ( gdept(ik) - gdept(ik-1) )
                        saldta(ji,jj,ik,jl) = (1.-zl) * saldta(ji,jj,ik,jl) +zl * saldta(ji,jj,ik-1,jl) 
                     ENDIF
                  END DO
               END DO
            ENDIF
         END DO
         

!         IF(lwp) THEN
!            WRITE(numout,*)' salinity month ',nsal1,nsal2
!            WRITE(numout,*)
!            WRITE(numout,*) ' month = ',nsal1,'  level = 1'
!            CALL prihre(saldta(:,:,1,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
!            WRITE(numout,*) ' month = ',nsal1,'  level = ',jpk/2
!            CALL prihre(saldta(:,:,jpk/2,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
!            WRITE(numout,*) ' month = ',nsal1,'  level = ',jpkm1
!            CALL prihre(saldta(:,:,jpkm1,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
!         ENDIF
      ENDIF
      
 
      ! 3. At every time step compute salinity data
      ! -------------------------------------------
      zxy = FLOAT(nday + 15 - 30*i15)/30.
!!DB: Note that for monthly data it is not necessary to do this every kt.
!!    So in the future use a code fragment like: 
!      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
!         s_dta(:,:,:) = (1.-zxy) * saldta(:,:,:,1) + zxy * saldta(:,:,:,2)
!      endif

      s_dta(:,:,:) = ( 1.- zxy ) * saldta(:,:,:,1) + zxy * saldta(:,:,:,2)

   END SUBROUTINE dta_sal

#else
   !!----------------------------------------------------------------------
   !!   Default option:                                    NO salinity data
   !!----------------------------------------------------------------------
   USE in_out_manager
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtasal = .FALSE.   !: salinity data flag
CONTAINS
   SUBROUTINE dta_sal( kt )        ! Empty routine
      if(lwp) WRITE(numout,*) 'dta_sal: You should not have seen this print! error?', kt
   END SUBROUTINE dta_sal
#endif
   !!======================================================================
END MODULE dtasal
