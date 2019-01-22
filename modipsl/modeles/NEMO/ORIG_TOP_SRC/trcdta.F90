MODULE trcdta
   !!======================================================================
   !!                     ***  MODULE  trcdta  ***
   !! Ocean data  :  reads passive tracer data 
   !!=====================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trcdta.F90,v 1.6 2006/04/10 15:40:28 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

#if  defined key_passivetrc && defined key_trc_dta
   !!----------------------------------------------------------------------
   !!   'key_trc_dta'                           3D tracer data field
   !!----------------------------------------------------------------------
   !!   dta_trc      : read ocean passive tracer data
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc
   USE trc
   USE par_sms
   USE lib_print

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dta_trc   ! called by trcdtr.F90 and trcdmp.F90

   !! * Shared module variables
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk,jptra) ::   &  !:
      trdta             !: tracer data at given time-step

   !! * Module variables
   REAL(wp), DIMENSION(jpi,jpj,jpk,jptra,2) ::   &
      tracdta            ! tracer data at two consecutive times
   INTEGER , DIMENSION(jptra) :: &
      nlectr  ,   &    !!: switch for reading once
      ntrc1   ,   &    !!: number of first month when reading 12 monthly value
      ntrc2            !!: number of second month when reading 12 monthly value

   !! * Substitutions
#  include "passivetrc_substitute.h90"

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LODYC-IPSL  (2003)
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !!   Default case                                            NetCDF file
   !!----------------------------------------------------------------------
   
   SUBROUTINE dta_trc( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_trc  ***
      !!
      !! ** Purpose :   Reads passive tracer data (Levitus monthly data)
      !!
      !! ** Method  :   Read on unit numtr the interpolated tracer concentra-
      !!      tion onto the global grid. Data begin at january. 
      !!      The value is centered at the middle of month. 
      !!      In the opa model, kt=1 agree with january 1. 
      !!      At each time step, a linear interpolation is applied between 
      !!      two monthly values.
      !!
      !! History :
      !!   8.2  !  02-04  (O. Aumont)  Original code
      !!   9.0  !  04-03  (C. Ethe)    
      !!   9.0  !  05-03  (O. Aumont and A. El Moussaoui) F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Arguments
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step

      !! * Local declarations
      INTEGER :: ji, jj, jn, jl 
      INTEGER, PARAMETER ::  &
         jpmois  = 12        ! number of months

      INTEGER ::   &
         imois, iman, i15, itime, ik, &  ! temporary integers 
         ipi, ipj, ipk                              !    "        "
      INTEGER :: istep(jpmois)
      CHARACTER (len=39) :: clname(jptra)
      REAL(wp), DIMENSION (jpi,jpj) ::  zlon, zlat
      REAL(wp), DIMENSION (jpk) ::  zlev
      REAL(wp) :: zdate0, zxy, zl
      !!----------------------------------------------------------------------

      DO jn = 1, jptra

         IF( lutini(jn) ) THEN 

            IF ( kt == nit000 ) THEN
               !! 3D tracer data
               IF(lwp)WRITE(numout,*)
               IF(lwp)WRITE(numout,*) ' dta_trc: reading tracer' 
               IF(lwp)WRITE(numout,*) ' data file ', jn, ctrcnm(jn)
               IF(lwp)WRITE(numout,*)
               nlectr(jn) = 0
            ENDIF
            ! Initialization
            iman = jpmois
            i15  = nday/16
            imois = nmonth + i15 -1
            IF( imois == 0 ) imois = iman
            itime = jpmois
            ipi = jpiglo
            ipj = jpjglo
            ipk = jpk

            ! First call kt=nit000
            ! --------------------

            IF ( kt == nit000 .AND. nlectr(jn) == 0 ) THEN
               ntrc1(jn) = 0
               IF(lwp) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) ' Tracer data  fields' 
                  WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~'
                  WRITE(numout,*) ' NetCDF FORMAT'
                  WRITE(numout,*)
               ENDIF

               ! open file 
#if defined key_trc_pisces
               clname(jn) = 'LEVITUS_'//ctrcnm(jn)
#else
               itime=1
               clname(jn) = ctrcnm(jn)
#endif
               CALL flinopen(TRIM(clname(jn)),mig(1),nlci,mjg(1),nlcj,    &
                  .FALSE.,ipi,ipj,ipk,zlon,zlat,zlev,itime,    &
                  istep,zdate0,rdt,numtr(jn)               )

#if defined key_trc_pisces
               ! title, dimensions and tests
               IF( itime /= jpmois ) THEN
                  IF(lwp) THEN
                     WRITE(numout,*) ' '
                     WRITE(numout,*) 'problem with time coordinates'
                     WRITE(numout,*) ' itime ',itime,' jpmois ',jpmois
                  ENDIF
                  STOP 'dta_trc'
               ENDIF
#endif

               IF( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) THEN
                  IF(lwp) THEN
                     WRITE(numout,*) ' '
                     WRITE(numout,*) 'problem with dimensions'
                     WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
                     WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
                     WRITE(numout,*) ' ipk ',ipk,' jpk ',jpk
                  ENDIF
                  STOP 'dta_trc'
               ENDIF
               IF(lwp)WRITE(numout,*) itime,istep(1),zdate0,rdt,numtr(jn)
               trdta(:,:,:,jn) = 0.

            ENDIF

#if defined key_trc_pisces
            ! Read montly file
            IF( ( kt == nit000 .AND. nlectr(jn) == 0)   & 
               .OR. imois /= ntrc1(jn) ) THEN
               nlectr(jn) = 1

               ! Calendar computation

               ! ntrc1 number of the first file record used in the simulation
               ! ntrc2 number of the last  file record

               ntrc1(jn) = imois
               ntrc2(jn) = ntrc1(jn) + 1
               ntrc1(jn) = MOD( ntrc1(jn), iman )
               IF ( ntrc1(jn) == 0 ) ntrc1(jn) = iman
               ntrc2(jn) = MOD( ntrc2(jn), iman )
               IF ( ntrc2(jn) == 0 ) ntrc2(jn) = iman
               IF(lwp) WRITE(numout,*) 'first record file used ntrc1 ', ntrc1(jn) 
               IF(lwp) WRITE(numout,*) 'last  record file used ntrc2 ', ntrc2(jn)

               ! Read montly passive tracer data Levitus 

               CALL flinget( numtr(jn),ctrcnm(jn),jpidta,jpjdta,jpk,    &
                  jpmois,ntrc1(jn),ntrc1(jn),mig(1),nlci,mjg(1),nlcj,  &
                  tracdta(1:nlci,1:nlcj,1:jpk,jn,1)                  )

               CALL flinget( numtr(jn),ctrcnm(jn),jpidta,jpjdta,jpk,     &
                  jpmois,ntrc2(jn),ntrc2(jn),mig(1),nlci,mjg(1),nlcj,   &
                  tracdta(1:nlci,1:nlcj,1:jpk,jn,2)                  )

               IF(lwp) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) ' read tracer data ', ctrcnm(jn),' ok'
                  WRITE(numout,*)
               ENDIF

               ! Apply Mask
               DO jl = 1, 2
                  tracdta(:,:,:  ,jn,jl) = tracdta(:,:,:,jn,jl) * tmask(:,:,:) 
                  tracdta(:,:,jpk,jn,jl) = 0.
                  IF( lk_zps ) THEN                ! z-coord. with partial steps
                     DO jj = 1, jpj                ! interpolation of temperature at the last level
                        DO ji = 1, jpi
                           ik = mbathy(ji,jj) - 1
                           IF( ik > 2 ) THEN
                              zl = ( gdept(ik) - fsdept(ji,jj,ik) ) / ( gdept(ik) - gdept(ik-1) )
                              tracdta(ji,jj,ik,jn,jl) = (1.-zl) * tracdta(ji,jj,ik,jn,jl) + zl * tracdta(ji,jj,ik-1,jn,jl)
                           ENDIF
                        END DO
                     END DO
                  ENDIF

               END DO

            ENDIF

            IF(lwp) THEN
               WRITE(numout,*) ctrcnm(jn), 'Levitus month ', ntrc1(jn),   &
                  ntrc2(jn)
               WRITE(numout,*)
               WRITE(numout,*) ' Levitus month = ', ntrc1(jn),   &
                  '  level = 1'
               CALL prihre( tracdta(1,1,1,jn,1), jpi, jpj, 1, jpi, 20, 1   &
                  ,jpj, 20, 1., numout )
               WRITE(numout,*) ' Levitus month = ', ntrc1(jn),    &
                  '  level = ',jpk/2
               CALL prihre( tracdta(1,1,jpk/2,jn,1), jpi, jpj, 1, jpi,    &
                  20, 1, jpj, 20, 1., numout )
               WRITE(numout,*) ' Levitus month = ',ntrc1(jn)     &
                  ,'  level = ',jpkm1
               CALL prihre( tracdta(1,1,jpkm1,jn,1), jpi, jpj, 1, jpi,     &
                  20, 1, jpj, 20, 1., numout )
            ENDIF

            ! At every time step compute temperature data

            zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
            trdta(:,:,:,jn)=  ( 1. - zxy ) * tracdta(:,:,:,jn,1)    &
               +       zxy   * tracdta(:,:,:,jn,2) 

            IF( jn == jpno3) trdta(:,:,:,jn) = trdta(:,:,:,jn) * 7.6E-6
            IF( jn == jpdic) trdta(:,:,:,jn) = trdta(:,:,:,jn) * 1.E-6
            IF( jn == jptal) trdta(:,:,:,jn) = trdta(:,:,:,jn) * 1.E-6
            IF( jn == jpoxy) trdta(:,:,:,jn) = trdta(:,:,:,jn) * 44.6E-6
            IF( jn == jpsil) trdta(:,:,:,jn) = trdta(:,:,:,jn) * 1.E-6
            IF( jn == jppo4) trdta(:,:,:,jn) = trdta(:,:,:,jn) * 122.E-6
#else
            ! Read init file only
            IF( kt == nit000  ) THEN
               CALL flinget( numtr(jn),ctrcnm(jn),jpidta,jpjdta,jpk,    &
               1,1,1,mig(1),nlci,mjg(1),nlcj,  &
               trdta(1:nlci,1:nlcj,1:jpk,jn) )
               trdta(:,:,:,jn)=trdta(:,:,:,jn)*tmask(:,:,:)
            ENDIF 
#endif

        ENDIF

       END DO

   END SUBROUTINE dta_trc

#else

   !!----------------------------------------------------------------------
   !!   Default case                        NO 3D passive tracer data field
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dta_trc( kt )        ! Empty routine
      WRITE(*,*) 'dta_trc: You should not have seen this print! error?', kt
   END SUBROUTINE dta_trc

#endif

END MODULE trcdta
