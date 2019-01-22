MODULE limdmp
   !!======================================================================
   !!                       ***  MODULE limdmp   ***
   !!  Ice model : restoring Ice thickness and Fraction leads
   !!======================================================================
!byoung #if defined key_ice_lim
#if defined key_ice_lim && defined key_tradmp
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_dmp      : ice model damping
   !!----------------------------------------------------------------------
   !! * Modules used
   USE in_out_manager  ! I/O manager
   USE ice
   USE ice_oce
   USE tradmp
   USE dom_oce
   USE oce
   USE daymod          ! calendar
   
   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC lim_dmp     ! called by ice_step
   
   !! * Shared module variables
   CHARACTER (len=38) ::   &
      cl_icedata = 'ice_damping_ATL4.nc'
   INTEGER ::   &
        nice1      ,   &  ! first record used
        nice2             ! second record used
   
    REAL(wp), DIMENSION(jpi,jpj,2) ::   &
         hicif_data ,   & ! ice thickness data at two consecutive times
         frld_data        ! fraction lead data at two consecutive times

    REAL(wp), DIMENSION(jpi,jpj) ::   &
         hicif_dta ,   &  ! ice thickness at a given time
         frld_dta         ! fraction lead at a given time

   !! * Substitution
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   LIM 2.0 , UCL-LOCEAN-IPSL  (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdmp.F90,v 1.1 2006/03/21 08:42:21 opalod Exp $
   !! This software is governed by the CeCILL licence see !modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_dmp(kt)
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_dmp ***
      !!
      !! ** purpose : ice model damping : restoring ice thickness and 
      !!              fraction leads
      !!
      !! ** method  : the key_tradmp must be used to compute resto(:,:) coef.
      !!     
      !! ** action :
      !!
      !! History :
      !!
      !!   2.0  !  04-04 (S. Theetten) Original
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step

      !! * Local Variables
      INTEGER  ::   ji, jj         ! dummy loop indices
      !!---------------------------------------------------------------------
    
      CALL dta_lim(kt)

      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.

            hicif(ji,jj) = hicif(ji,jj) - rdt_ice * resto(ji,jj,1) * ( hicif(ji,jj) -  hicif_dta(ji,jj))
            frld(ji,jj)  = frld(ji,jj)  - rdt_ice * resto(ji,jj,1) * ( frld(ji,jj)  - frld_dta(ji,jj))  

         ENDDO
      ENDDO

   END SUBROUTINE lim_dmp



   SUBROUTINE dta_lim( kt ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_lim  ***
      !!
      !! ** Purpose :   Reads monthly ice thickness and fraction lead  data
      !!
      !! ** Method  :   Read on unit numicedt the interpolated ice variable
      !!      onto the model grid.
      !!      Data begin at january.
      !!      The value is centered at the middle of month.
      !!      In the opa model, kt=1 agree with january 1.
      !!      At each time step, a linear interpolation is applied between
      !!      two monthly values.
      !!      
      !!
      !! ** Action  :   define hicif_dta and frld_dta arrays at time-step kt
      !!
      !! History :
      !!   2.0   !   04-04 (S. Theetten) Original
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step

      !! * Local declarations
      INTEGER, PARAMETER ::   jpmois = 12       ! number of month
     
      INTEGER ::   &
         imois, iman, itime ,    &  ! temporary integers
         i15, ipi, ipj, ipk         !    "          "

      INTEGER, DIMENSION(jpmois) ::   istep
      REAL(wp) ::   zxy, zdate0, zdt
      REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
      REAL(wp), DIMENSION(jpk) ::   zlev
      !!----------------------------------------------------------------------

      ! 0. Initialization
      ! -----------------
      iman  = jpmois
      i15   = nday / 16
      imois = nmonth + i15 - 1
      IF( imois == 0 )   imois = iman

      itime = jpmois
      ipi=jpiglo
      ipj=jpjglo
      ipk=1
      zdt=rdt

      ! 1. First call kt=nit000
      ! -----------------------

      IF( kt == nit000 ) THEN
         nice1 = 0
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dtalim : Ice thickness and lead fraction  monthly fields'
         IF(lwp) WRITE(numout,*) '~~~~~~'
         IF(lwp) WRITE(numout,*) '             NetCDF FORMAT'
         IF(lwp) WRITE(numout,*)
         
         ! open file
         
         CALL flinopen( TRIM(cl_icedata), mig(1), nlci , mjg(1),  nlcj, .FALSE.,  &
            &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, numice_dmp )

          ! title, dimensions and tests
         IF( itime /= jpmois ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'problem with time coordinates'
               WRITE(numout,*) ' itime ',itime,' jpmois ',jpmois
            ENDIF
            STOP 'dta_lim'
         ENDIF
         IF( ipi /= jpidta .OR. ipj /= jpjdta ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'problem with dimensions'
               WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
            ENDIF
            STOP 'dta_lim'
         ENDIF
         IF(lwp) WRITE(numout,*) itime,istep,zdate0,zdt,numice_dmp

      ENDIF


      ! 2. Read monthly file
      ! -------------------

      IF( ( kt == nit000 ) .OR. imois /= nice1 ) THEN

         ! Calendar computation
         
         nice1 = imois        ! first file record used 
         nice2 = nice1 + 1    ! last  file record used
         nice1 = MOD( nice1, iman )
         IF( nice1 == 0 )   nice1 = iman
         nice2 = MOD( nice2, iman )
         IF( nice2 == 0 )   nice2 = iman
         IF(lwp) WRITE(numout,*) 'first record file used nice1 ', nice1
         IF(lwp) WRITE(numout,*) 'last  record file used nice2 ', nice2
         
         ! Read monthly ice thickness Levitus 
         
         CALL flinget( numice_dmp, 'icethic', jpidta, jpjdta, jpk,  &
            &          jpmois, nice1, nice1, mig(1), nlci, mjg(1), nlcj, hicif_data(1:nlci,1:nlcj,1) )
         CALL flinget( numice_dmp, 'icethic', jpidta, jpjdta, jpk,  &
            &          jpmois, nice2, nice2, mig(1), nlci, mjg(1), nlcj, hicif_data(1:nlci,1:nlcj,2) )
         
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' read ice thickness ok'
         IF(lwp) WRITE(numout,*)

         ! Read monthly ice thickness Levitus 
         
         CALL flinget( numice_dmp, 'ileadfra', jpidta, jpjdta, jpk,  &
            &          jpmois, nice1, nice1, mig(1), nlci, mjg(1), nlcj, frld_data(1:nlci,1:nlcj,1) )
         CALL flinget( numice_dmp, 'ileadfra', jpidta, jpjdta, jpk,  &
            &          jpmois, nice2, nice2, mig(1), nlci, mjg(1), nlcj, frld_data(1:nlci,1:nlcj,2) )
         
         ! The fraction lead read in the file is in fact the 
         ! ice concentration which is 1 - the fraction lead
         frld_data = 1 - frld_data          
         
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' read fraction lead ok'
         IF(lwp) WRITE(numout,*)


         IF(lwp) THEN
            WRITE(numout,*) ' Ice thickness month ', nice1,' and ', nice2
            WRITE(numout,*)
            WRITE(numout,*) ' Ice thickness month = ', nice1
            CALL prihre( hicif_data(1,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
            WRITE(numout,*)
            WRITE(numout,*) ' Fraction lead months ', nice1,' and ', nice2
            WRITE(numout,*)
            WRITE(numout,*) ' Fraction lead month = ', nice1
            CALL prihre( frld_data(1,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         ENDIF
         
         ! 2. At every time step compute ice thickness and fraction lead data
         ! ------------------------------------------------------------------
         
         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
         hicif_dta(:,:) = (1.-zxy) * hicif_data(:,:,1) + zxy * hicif_data(:,:,2)
         frld_dta(:,:) = (1.-zxy) * frld_data(:,:,1) + zxy * frld_data(:,:,2)

      ENDIF


   END SUBROUTINE dta_lim

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty Module                  No ice damping
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_dmp(kt)        ! Dummy routine
      WRITE(*,*) 'lim_dmp: You should not see this print! error? ', kt
   END SUBROUTINE lim_dmp
#endif

   !!======================================================================

END MODULE limdmp
