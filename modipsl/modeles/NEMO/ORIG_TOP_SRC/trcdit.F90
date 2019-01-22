MODULE trcdit
   !!----------------------------------------------------------------------
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trcdit.F90,v 1.6 2006/04/10 15:40:28 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   USE dianam    ! build name of file (routine)
   USE in_out_manager  ! I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trcdit_wr
   PUBLIC trcdid_wr
   PUBLIC trcdii_wr
   PUBLIC trcdib_wr

   !! * Module variables
   INTEGER            ::  &
      nit5     ,  &   !!: id for tracer output file
      ndepit5  ,  &   !!: id for depth mesh
      nhorit5  ,  &   !!: id for horizontal mesh
      ndimt50  ,  &   !!: number of ocean points in index array
      ndimt51         !!: number of ocean points in index array
   REAL(wp) :: zjulian
   INTEGER , DIMENSION (jpij*jpk) ::  ndext50 !!: integer arrays for ocean 3D index
   INTEGER , DIMENSION (jpij)     ::  ndext51 !!: integer arrays for ocean surface index
#    if defined key_passivetrc && defined key_trc_diaadd
   INTEGER            :: &
      nitd     ,  &   !!: id for additional array output file
      ndepitd  ,  &   !!: id for depth mesh
      nhoritd         !!: id for horizontal mesh
#    endif
#    if defined key_passivetrc && defined key_trc_diatrd
   INTEGER , DIMENSION (jptra)  :: &
      nit6    ,   &   !!: id for additional array output file
      ndepit6 ,   &   !!: id for depth mesh
      nhorit6         !!: id for horizontal mesh
#    endif
#    if defined key_passivetrc && defined key_trc_diabio
   INTEGER            :: &
      nitb     ,   &  !!:  id for additional array output FILE
      ndepitb  ,   &  !!:  id for depth mesh
      nhoritb         !!:  id for horizontal mesh

#    endif


   !! * Substitutions
#  include "passivetrc_substitute.h90"

CONTAINS

#    if defined key_passivetrc

      SUBROUTINE trcdit_wr(kt,kindic)
   !!===========================================================================================
   !!
   !!                       ROUTINE trcdit_wr
   !!===========================================================================================
   !!
   !! Purpose :
   !!---------
   !!          Standard output of passive tracer : concentration fields
   !!
   !!
   !! Method :
   !! -------
   !!
   !!        At the beginning of the first time step (nit000), define all
   !!        the NETCDF files and fields for concentration of passive tracer
   !!
   !!        At each time step call histdef to compute the mean if necessary
   !!        Each nwritetrc time step, output the instantaneous or mean fields
   !!
   !!        IF kindic <0, output of fields before the model interruption.
   !!        IF kindic =0, time step loop
   !!        IF kindic >0, output of fields before the time step loop
   !!
   !! Input :
   !! -----
   !!   argument
   !!           kt              : time step
   !!           kindic          : indicator of abnormal termination
   !!
   !! EXTERNAL :
   !! --------
   !! prihre, hist..., dianam
   !!
   !! History:
   !! --------
   !!   original  : 95-01  passive tracers  (M. Levy)
   !!   additions : 98-01 (C. Levy) NETCDF format using ioipsl interface
   !!   additions : 99-01 (M.A. Foujols) adapted for passive tracer
   !!   additions : 99-09 (M.A. Foujols) split into three parts
   !!   05-03 (O. Aumont and A. El Moussaoui) F90
   !!==================================================================================================!

      !! Modules used
      USE ioipsl


      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,kindic         ! ocean time-step

      !! * Local declarations
      INTEGER :: jn
      LOGICAL :: ll_print = .FALSE.

      CHARACTER (len=40) :: clhstnam, clop
      CHARACTER (len=20) :: cltra, cltrau
      CHARACTER (len=80) :: cltral

      REAL(wp) :: zsto, zout, zdt
      INTEGER  :: iimi, iima, ijmi, ijma, ipk, it
!
! 0. Initialisation
! -----------------

! local variable for debugging
      ll_print = .FALSE.
      ll_print = ll_print .AND. lwp

! Define frequency of output and means

      zdt = rdt
#        if defined key_diainstant
      zsto=nwritetrc*rdt
      clop='inst(only(x))'
#        else
      zsto=zdt
      clop='ave(only(x))'
#        endif
      zout=nwritetrc*zdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      it = kt - nit000 + 1

! 1. Define NETCDF files and fields at beginning of first time step
! -----------------------------------------------------------------

      IF(ll_print)WRITE(numout,*)'trcdit_wr kt=',kt,' kindic ',kindic
      IF(kt == nit000) THEN

! Compute julian date from starting date of the run

         CALL ymds2ju(nyear,nmonth,nday,0.0,zjulian)
         IF(lwp)WRITE(numout,*)' '  
         IF(lwp)WRITE(numout,*)' Date 0 used :',nit000     &
       &     ,' YEAR ',nyear,' MONTH ',nmonth,' DAY ',nday   &
       &     ,'Julian day : ',zjulian    
         IF(lwp)WRITE(numout,*) ' indexes of zoom = ', iimi, iima, ijmi, ijma,  &
                                 ' limit storage in depth = ', ipk


! Define the NETCDF files for passive tracer concentration

         CALL dia_nam(clhstnam,nwritetrc,'ptrc_T')

         IF(lwp)WRITE(numout,*)" Name of NETCDF file ", clhstnam
! Horizontal grid : glamt and gphit
 
         CALL histbeg(clhstnam, jpi, glamt, jpj, gphit,     &
         &    iimi, iima-iimi+1, ijmi, ijma-ijmi+1,         & 
         &    0, zjulian, zdt, nhorit5, nit5 , domain_id=nidom)
! Vertical grid for tracer : gdept
         CALL histvert(nit5, 'deptht', 'Vertical T levels', &
         &    'm', ipk, gdept, ndepit5)

! Index of ocean points in 3D and 2D (surface)
         CALL wheneq(jpi*jpj*ipk,tmask,1,1.,ndext50,ndimt50)
         CALL wheneq(jpi*jpj,tmask,1,1.,ndext51,ndimt51)

! Declare all the output fields as NETCDF variables

! tracer concentrations

         DO jn=1,jptra
           cltra=ctrcnm(jn)    ! short title for tracer
           cltral=ctrcnl(jn)   ! long title for tracer
           cltrau=ctrcun(jn)   ! UNIT for tracer
           CALL histdef(nit5, cltra, cltral, cltrau, jpi, jpj, nhorit5,  &
         &          ipk, 1, ipk,  ndepit5, 32, clop, zsto, zout) 
         END DO           

! CLOSE netcdf Files
          
         CALL histend(nit5)

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization in trcdit_wr'
         IF(ll_print) CALL FLUSH(numout )

      ENDIF

! 2. Start writing data
! ---------------------

! tracer concentrations

      IF( lwp .AND. MOD( kt, nwritetrc ) == 0 ) THEN
         WRITE(numout,*) 'trcdit_wr : write NetCDF passive tracer concentrations at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF

      DO jn=1,jptra
         cltra=ctrcnm(jn) ! short title for tracer
         CALL histwrite(nit5, cltra, it, trn(:,:,:,jn), ndimt50,   &
      &          ndext50)
      END DO 

! synchronise FILE

      IF( MOD( kt, nwritetrc ) == 0 .OR. kindic < 0 ) THEN
              CALL histsync(nit5)
      ENDIF

! 3. Closing all files
! --------------------
      IF( kt == nitend .OR. kindic < 0 ) THEN
          CALL histclo(nit5)
      ENDIF

END SUBROUTINE trcdit_wr

#    else

! no passive tracers

SUBROUTINE trcdit_wr(kt,kindic)
     !!! no passive tracers
     INTEGER, INTENT ( in ) :: kt, kindic
     WRITE(*,*) 'trcdit_wr: You should not have seen this print! error?', kt, kindic
END SUBROUTINE trcdit_wr

#    endif

#    if defined key_passivetrc && defined key_trc_diatrd

      SUBROUTINE trcdid_wr(kt,kindic)
 !!===========================================================================================
   !!
   !!                       ROUTINE trcdid_wr
   !!===========================================================================================
   !!
   !! Purpose :
   !!---------
   !!          output of opa: passive tracer dynamical trends
   !!
   !!
   !! Method :
   !! -------
   !!
   !!        At the beginning of the first time step (nit000), define all
   !!        the NETCDF files and fields for dynamical trends of tracers
   !!
   !!        At each time step call histdef to compute the mean if necessary
   !!        Each nwritetrd time step, output the instantaneous or mean fields
   !!
   !!        IF kindic <0, output of fields before the model interruption.
   !!        IF kindic =0, time step loop
   !!        IF kindic >0, output of fields before the time step loop
   !!
   !! Input :
   !! -----
   !!   argument
   !!           kt              : time step
   !!           kindic          : indicator of abnormal termination
   !!
   !! Output :
   !! ------
   !!   file
   !!           "clhstnam" files : one for concentration
   !!
   !! History:
   !! --------
   !!   original  : 95-01  passive tracers  (M. Levy)
   !!   additions : 98-01 (C. Levy) NETCDF format using ioipsl interface
   !!   additions : 99-01 (M.A. Foujols) adapted for passive tracer
   !!   additions : 99-09 (M.A. Foujols) split into three parts
   !!   additions : 01-06 (Mehdi B, Elodie K): suppress initialization
   !!                                          of nit6,nhorit6,ndepit6
   !!   05-03 (O. Aumont and A. El Moussaoui) F90
   !!==================================================================================================!

      !! Modules used
      USE ioipsl

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,kindic         ! ocean time-step

      INTEGER :: jn, jl
      LOGICAL :: ll_print = .FALSE.

      CHARACTER (len=40) :: clhstnam, clop
      CHARACTER (len=20) :: cltra, cltrau
      CHARACTER (len=80) :: cltral
      CHARACTER (len=10) :: csuff

      REAL(wp) :: zsto, zout, zdt
      INTEGER :: iimi, iima, ijmi, ijma, ipk, it

!
! 0. Initialisation
! -----------------

! local variable for debugging
      ll_print = .FALSE.
      ll_print = ll_print .AND. lwp
!
! Define frequency of output and means
!
      zdt = rdt
#        if defined key_diainstant
      zsto=nwritetrd*rdt
      clop='inst(only(x))'
#        else
      zsto=zdt
      clop='ave(only(x))'
#        endif
      zout=nwritetrd*zdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      it = kt - nit000 + 1

! Define the NETCDF files (one per tracer)
!
      IF(ll_print)WRITE(numout,*)'trcdid kt=',kt,' kindic ',kindic
      IF(kt == nit000) THEN

          DO jn=1,jptra

            IF (luttrd(jn)) THEN

! Define the file for dynamical trends - one per each tracer IF required

         IF(lwp)WRITE(numout,*) ' indexes of zoom = ', iimi, iima, ijmi, ijma,  &
                                 ' limit storage in depth = ', ipk
                csuff='DY_'//ctrcnm(jn)
                CALL dia_nam(clhstnam,nwritetrd,csuff)
                IF(lwp)WRITE(numout,*)     &
                &      " Name of NETCDF file for dynamical trends",   &
                &      " of tracer number : ",clhstnam

                CALL histbeg(clhstnam, jpi, glamt, jpj, gphit,   &
                &    iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
                &    0, zjulian, rdt, nhorit6(jn),               &
                &    nit6(jn) , domain_id=nidom)

! Vertical grid for tracer trend - one per each tracer IF needed
                CALL histvert(nit6(jn), 'deptht', 'Vertical T levels',  &
                &    'm', ipk, gdept, ndepit6(jn)) 


            END IF
          END DO

! Declare all the output fields as NETCDF variables


! trends for tracer concentrations
          DO jn=1,jptra
            IF (luttrd(jn)) THEN
                DO jl=1,jpdiatrc
                  IF (jl.eq.1) THEN
! short and long title for x advection for tracer
                      WRITE (cltra,'("XAD_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("X advective trend for ",58a)')  &
                      &      ctrcnl(jn)(1:58)
                  END IF
                  IF (jl.eq.2)  THEN
! short and long title for y advection for tracer
                      WRITE (cltra,'("YAD_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("Y advective trend for ",58a)')  &
                      &      ctrcnl(jn)(1:58)
                  END IF
                  IF (jl.eq.3)  THEN
! short and long title for Z advection for tracer
                      WRITE (cltra,'("ZAD_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("Z advective trend for ",58a)')  &
                      &      ctrcnl(jn)(1:58)
                  END IF
                  IF (jl.eq.4)  THEN
! short and long title for X diffusion for tracer
                      WRITE (cltra,'("XDF_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("X diffusion trend for ",58a)')  &
                      &      ctrcnl(jn)(1:58)
                  END IF
                  IF (jl.eq.5)  THEN
! short and long title for Y diffusion for tracer
                      WRITE (cltra,'("YDF_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("Y diffusion trend for ",58a)')  &
                      &      ctrcnl(jn)(1:58)
                  END IF
                  IF (jl.eq.6)  THEN
! short and long title for Z diffusion for tracer
                      WRITE (cltra,'("ZDF_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("Z diffusion trend for ",58a)')  &
                      &      ctrcnl(jn)(1:58)
                  END IF
# if defined key_trc_ldfeiv
                  IF (jl.eq.7) THEN
! short and long title for x gent velocity for tracer
                      WRITE (cltra,'("XGV",16a)') ctrcnm(jn)
                      WRITE (cltral,'("X gent velocity trend for ",53a)')  &
                      &      ctrcnl(jn)(1:53)
                  END IF
                  IF (jl.eq.8)  THEN
! short and long title for y gent velocity for tracer
                      WRITE (cltra,'("YGV_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("Y gent velocity trend for ",53a)')  &
                      &      ctrcnl(jn)(1:53)
                  END IF
                  IF (jl.eq.9)  THEN
! short and long title for Z gent velocity for tracer
                      WRITE (cltra,'("ZGV_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("Z gent velocity trend for ",53a)')  &
                      &      ctrcnl(jn)(1:53)
                  END IF
# endif
# if defined key_trcdmp
                  IF (jl.eq.jpdiatrc)  THEN
! last trends for tracer damping : short and long title
                      WRITE (cltra,'("TDM_",16a)') ctrcnm(jn)
                      WRITE (cltral,'("Tracer damping trend for ",55a)')  &
                      &      ctrcnl(jn)(1:55)
                  END IF
# endif
                  call flush(numout)
                  cltrau=ctrcun(jn) ! UNIT for tracer /trends
                  CALL histdef(nit6(jn), cltra, cltral, cltrau, jpi,jpj,  &
                  &   nhorit6(jn), ipk, 1, ipk,  ndepit6(jn), 32, clop ,  &
                  &   zsto,zout)
                END DO
            END IF
          END DO

! CLOSE netcdf Files

          DO jn=1,jptra
             IF (luttrd(jn)) CALL histend(nit6(jn))
          END DO

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization in trcdid'
         IF(ll_print) CALL FLUSH(numout )

      ENDIF

! SOME diagnostics to DO first time

! 2. Start writing data
! ---------------------

! trends for tracer concentrations

      IF( lwp .AND. MOD( kt, nwritetrd ) == 0 ) THEN
         WRITE(numout,*) 'trcdid_wr : write NetCDF dynamical trends at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF

          DO jn=1,jptra
            IF (luttrd(jn)) THEN
                DO jl=1,jpdiatrc
                  IF (jl.eq.1) THEN
! short title for x advection for tracer
                      WRITE (cltra,'("XAD_",16a)') ctrcnm(jn)
                  END IF
                  IF (jl.eq.2)  THEN
! short title for y advection for tracer
                      WRITE (cltra,'("YAD_",16a)') ctrcnm(jn)
                  END IF
                  IF (jl.eq.3)  THEN
! short title for z advection for tracer
                      WRITE (cltra,'("ZAD_",16a)') ctrcnm(jn)
                  END IF
                  IF (jl.eq.4)  THEN
! short title for x diffusion for tracer
                      WRITE (cltra,'("XDF_",16a)') ctrcnm(jn)
                  END IF
                  IF (jl.eq.5)  THEN
! short title for y diffusion for tracer
                      WRITE (cltra,'("YDF_",16a)') ctrcnm(jn)
                  END IF
                  IF (jl.eq.6)  THEN
! short title for z diffusion for tracer
                      WRITE (cltra,'("ZDF_",16a)') ctrcnm(jn)
                  END IF
# if defined key_trc_ldfeiv
                  IF (jl.eq.7) THEN
! short for x gent velocity for tracer
                      WRITE (cltra,'("XGV_",16a)') ctrcnm(jn)
                  END IF
                  IF (jl.eq.8)  THEN
! short for y gent velocity for tracer
                      WRITE (cltra,'("YGV_",16a)') ctrcnm(jn)
                  END IF
                  IF (jl.eq.9)  THEN
! short title for Z gent velocity for tracer
                      WRITE (cltra,'("ZGV_",16a)') ctrcnm(jn)
                  END IF
# endif
# if defined key_trcdmp
                  IF (jl.eq.jpdiatrc) THEN
! short for x gent velocity for tracer
                      WRITE (cltra,'("TDM_",16a)') ctrcnm(jn)
                  END IF
# endif

                  CALL histwrite(nit6(jn), cltra, it, trtrd(:,:,:,ikeep(jn),jl)  &
                  &    ,ndimt50, ndext50)
                END DO
            END IF
          END DO

! synchronise FILE

      IF( MOD( kt, nwritetrd ) == 0 .OR. kindic < 0 ) THEN
          DO jn=1,jptra
             IF (luttrd(jn)) CALL histsync(nit6(jn))
          END DO
      ENDIF

! 3. Closing all files
! --------------------

      IF( kt == nitend .OR. kindic < 0 ) THEN
          DO jn=1,jptra
             IF (luttrd(jn)) CALL histclo(nit6(jn))
          END DO
      ENDIF

END SUBROUTINE trcdid_wr

#    else

SUBROUTINE trcdid_wr(kt,kindic)
     !!! no passive tracers
     INTEGER, INTENT ( in ) :: kt, kindic
     WRITE(*,*) 'trcdid_wr: You should not have seen this print! error?', kt, kindic
END SUBROUTINE trcdid_wr

#    endif

#    if defined key_passivetrc && defined key_trc_diaadd

      SUBROUTINE trcdii_wr(kt,kindic)
   !!===========================================================================================
   !!
   !!                       ROUTINE trcdii_wr
   !!===========================================================================================
   !!
   !! Purpose :
   !!---------
   !!          output of passive tracer : additional 2D and 3D arrays
   !!
   !!
   !! Method :
   !! -------
   !!
   !!        At the beginning of the first time step (nit000), define all
   !!        the NETCDF files and fields for additional arrays
   !!
   !!        At each time step call histdef to compute the mean if necessary
   !!        Each nwritetrc time step, output the instantaneous or mean fields
   !!
   !!
   !!        IF kindic <0, output of fields before the model interruption.
   !!        IF kindic =0, time step loop
   !!        IF kindic >0, output of fields before the time step loop
   !!
   !! Input :
   !! -----
   !!   argument
   !!           kt              : time step
   !!           kindic          : indicator of abnormal termination
   !!
   !! EXTERNAL :
   !! --------
   !! prihre, hist..., dianam
   !!
   !! History:
   !! --------
   !!   original  : 95-01  passive tracers  (M. Levy)
   !!   additions : 98-01 (C. Levy) NETCDF format using ioipsl interface
   !!   additions : 99-01 (M.A. Foujols) adapted for passive tracer
   !!   additions : 99-09 (M.A. Foujols) split into three parts
   !!   05-03 (O. Aumont and A. El Moussaoui) F90
   !!==================================================================================================!

      !! Modules used
      USE ioipsl

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,kindic         ! ocean time-step

      INTEGER :: jn
      LOGICAL :: ll_print = .FALSE.

      CHARACTER (len=40) :: clhstnam, clop
      CHARACTER (len=20) :: cltra, cltrau
      CHARACTER (len=80) :: cltral

      REAL(wp) :: zsto, zout, zdt
      INTEGER :: iimi, iima, ijmi, ijma, ipk, it

!
! 0. Initialisation
! -----------------

! local variable for debugging
      ll_print = .FALSE.
      ll_print = ll_print .AND. lwp
!
! Define frequency of output and means
!
      zdt = rdt
#        if defined key_diainstant
      zsto=nwriteadd*zdt
      clop='inst(only(x))'
#        else
      zsto=zdt
      clop='ave(only(x))'
#        endif
      zout=nwriteadd*zdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      it = kt - nit000 + 1

! 1. Define NETCDF files and fields at beginning of first time step
! -----------------------------------------------------------------

      IF(ll_print)WRITE(numout,*)'trcdii_wr kt=',kt,' kindic ',kindic
      IF(kt == nit000) THEN

! Define the NETCDF files for additional arrays : 2D or 3D

! Define the T grid file for tracer auxiliary files

          CALL dia_nam(clhstnam,nwrite,'diad_T')
          IF(lwp)WRITE(numout,*)" Name of NETCDF file ", clhstnam

! Define a netcdf FILE for 2d and 3d arrays

          CALL histbeg(clhstnam, jpi, glamt, jpj, gphit,     &
          &    iimi, iima-iimi+1, ijmi, ijma-ijmi+1,         &
          &    0, zjulian, zdt, nhoritd, nitd , domain_id=nidom)

! Vertical grid for 2d and 3d arrays

          CALL histvert(nitd, 'deptht', 'Vertical T levels', &
          &    'm', ipk, gdept, ndepitd)


! Declare all the output fields as NETCDF variables

! more 3D horizontal arrays

          DO jn=1,jpdia3d
            cltra=ctrc3d(jn)    ! short title for 3D diagnostic
            cltral=ctrc3l(jn)   ! long title for 3D diagnostic
            cltrau=ctrc3u(jn)   ! UNIT for 3D diagnostic
            CALL histdef(nitd, cltra, cltral, cltrau, jpi, jpj, nhoritd,  &
            &    ipk, 1, ipk,  ndepitd, 32, clop, zsto, zout)
          END DO


! more 2D horizontal arrays

          DO jn=1,jpdia2d
            cltra=ctrc2d(jn)    ! short title for 2D diagnostic
            cltral=ctrc2l(jn)   ! long title for 2D diagnostic
            cltrau=ctrc2u(jn)   ! UNIT for 2D diagnostic
            CALL histdef(nitd, cltra, cltral, cltrau, jpi, jpj, nhoritd,  &
            &    1, 1, 1,  -99, 32, clop, zsto, zout)
          END DO

! TODO: more 2D vertical sections arrays : I or J indice fixed

! CLOSE netcdf Files

          CALL histend(nitd)

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization in trcdii_wr'
         IF(ll_print) CALL FLUSH(numout )

      ENDIF

! 2. Start writing data
! ---------------------

      IF( lwp .AND. MOD( kt, nwriteadd ) == 0 ) THEN
         WRITE(numout,*) 'trcdii_wr : write NetCDF additional arrays at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF

! more 3D horizontal arrays

          DO jn=1,jpdia3d
            cltra=ctrc3d(jn) ! short title for 3D diagnostic
            CALL histwrite(nitd, cltra, it, trc3d(:,:,:,jn), ndimt50  &
            &   ,ndext50)
          END DO

! more 2D horizontal arrays

          DO jn=1,jpdia2d
            cltra=ctrc2d(jn) ! short title for 2D diagnostic
            CALL histwrite(nitd, cltra, kt, trc2d(:,:,jn), ndimt51    &
            &   ,ndext51)
          END DO

! synchronise FILE

      IF( MOD( kt, nwriteadd ) == 0 .OR. kindic < 0 ) THEN
              CALL histsync(nitd)
      ENDIF

! 3. Closing all files
! --------------------

      IF( kt == nitend .OR. kindic < 0 ) THEN
          CALL histclo(nitd)
      ENDIF

END SUBROUTINE trcdii_wr

#    else

SUBROUTINE trcdii_wr(kt,kindic)
     !!! no passive tracers
     INTEGER, INTENT ( in ) :: kt, kindic
     WRITE(*,*) 'trcdii_wr: You should not have seen this print! error?', kt, kindic
END SUBROUTINE trcdii_wr

#    endif

#    if defined key_passivetrc && defined key_trc_diabio

      SUBROUTINE trcdib_wr(kt,kindic)
 !!===========================================================================================
   !!
   !!                       ROUTINE trcdib_wr
   !!===========================================================================================
   !!
   !! Purpose :
   !!---------
   !!          Specific output of opa: biological fields
   !!
   !!
   !! Method :
   !! -------
   !!
   !!        At the beginning of the first time step (nit000), define all
   !!        the NETCDF files and fields for biological fields
   !!
   !!        At each time step call histdef to compute the mean if necessary
   !!        Each nwritetrd time step, output the instantaneous or mean fields
   !!
   !!        IF kindic <0, output of fields before the model interruption.
   !!        IF kindic =0, time step loop
   !!        IF kindic >0, output of fields before the time step loop
   !!
   !! Input :
   !! -----
   !!   argument
   !!           kt              : time step
   !!           kindic          : indicator of abnormal termination
   !!
   !! Output :
   !! ------
   !!   file
   !!           "histname" files : at least one file for each grid
   !!
   !! History:
   !! --------
   !!   original  : 95-01  passive tracers  (M. Levy)
   !!   additions : 98-01 (C. Levy) NETCDF format using ioipsl interface
   !!   additions : 99-01 (M.A. Foujols) adapted for passive tracer
   !!   additions : 99-09 (M.A. Foujols) split into three parts
   !!   additions : 01-06 (E Kestenare) assign a parameter to name
   !!                                          individual tracers
   !!   additions : 05-03 (O. Aumont and A El Moussaoui) F90
   !!==================================================================================================!

      !! Modules used
      USE ioipsl
      USE sms

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,kindic         ! ocean time-step

      INTEGER :: ji, jj, jk, jn
      LOGICAL :: ll_print = .FALSE.

      CHARACTER (len=40) :: clhstnam, clop
      CHARACTER (len=20) :: cltra, cltrau
      CHARACTER (len=80) :: cltral

      REAL(wp) :: zsto, zout, zdt
      INTEGER  :: iimi, iima, ijmi, ijma, ipk, it

!
! 0. Initialisation
! -----------------

! local variable for debugging
      ll_print = .FALSE.
      ll_print = ll_print .AND. lwp
!
! Define frequency of output and means
!
      zdt = rdt
#        if defined key_diainstant
      zsto=nwritebio*zdt
      clop='inst(only(x))'
#        else
      zsto=zdt
      clop='ave(only(x))'
#        endif
      zout=nwritebio*zdt

      ! Define indices of the horizontal output zoom and vertical limit storage      iimi = 1      ;      iima = jpi
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      it = kt - nit000 + 1

! 1. Define NETCDF files and fields at beginning of first time step
! -----------------------------------------------------------------

      IF(ll_print)WRITE(numout,*)'trcdib_wr kt=',kt,' kindic ',kindic
      IF(kt == nit000) THEN

! Define the NETCDF files for biological trends

          CALL dia_nam(clhstnam,nwrite,'biolog')
          IF(lwp)WRITE(numout,*)        &
          &      " Name of NETCDF file for biological trends ",clhstnam
! Horizontal grid : glamt and gphit
          CALL histbeg(clhstnam, jpi, glamt, jpj, gphit,      &
          &    iimi, iima-iimi+1, ijmi, ijma-ijmi+1,          &
          &    0, zjulian, rdt, nhoritb, nitb , domain_id=nidom)
! Vertical grid for biological trends
          CALL histvert(nitb, 'deptht', 'Vertical T levels',  &
          &    'm', ipk, gdept, ndepitb)

! Declare all the output fields as NETCDF variables

! biological trends

          DO jn=1,jpdiabio
            cltra=ctrbio(jn)    ! short title for biological diagnostic
            cltral=ctrbil(jn)   ! long title for biological diagnostic
            cltrau=ctrbiu(jn)   ! UNIT for biological diagnostic
            CALL histdef(nitb, cltra, cltral, cltrau, jpi, jpj, nhoritb,  &
            &    ipk, 1, ipk,  ndepitb, 32, clop, zsto, zout)
          END DO

! CLOSE netcdf Files

          CALL histend(nitb)

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization in trcdib_wr'
         IF(ll_print) CALL FLUSH(numout )

     ENDIF

! 2. Start writing data
! ---------------------

! biological trends

      IF( lwp .AND. MOD( kt, nwritebio ) == 0 ) THEN
         WRITE(numout,*) 'trcdit_wr : write NetCDF biological trends at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF


      DO jn=1,jpdiabio
         cltra=ctrbio(jn)  ! short title for biological diagnostic
         CALL histwrite(nitb, cltra, kt, trbio(:,:,:,jn), ndimt50,ndext50)
      END DO

! synchronise FILE

      IF( MOD( kt, nwritebio ) == 0 .OR. kindic < 0 ) THEN
              CALL histsync(nitb)
      ENDIF

! 3. Closing all files
! --------------------
      IF( kt == nitend .OR. kindic < 0 ) THEN
          CALL histclo(nitb)
      ENDIF

END SUBROUTINE trcdib_wr

#    else

SUBROUTINE trcdib_wr(kt,kindic)
     !!! no passive tracers
     INTEGER, INTENT ( in ) :: kt, kindic
     WRITE(*,*) 'trcdib_wr: You should not have seen this print! error?', kt, kindic
END SUBROUTINE trcdib_wr

#    endif

END MODULE trcdit
