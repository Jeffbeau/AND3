MODULE lib_mpp
   !!======================================================================
   !!                       ***  MODULE  lib_mpp  ***
   !! Ocean numerics:  massively parallel processing librairy
   !!=====================================================================
#if   defined key_mpp_mpi   ||   defined key_mpp_shmem
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'     OR      MPI massively parallel processing library
   !!   'key_mpp_shmem'         SHMEM massively parallel processing library
   !!----------------------------------------------------------------------
   !!   mynode
   !!   mpparent
   !!   mppshmem
   !!   mpp_lnk     : generic interface (defined in lbclnk) for :
   !!                 mpp_lnk_2d, mpp_lnk_3d
   !!   mpp_lnk_3d_gather :  Message passing manadgement for two 3D arrays
   !!   mpp_lnk_e   : interface defined in lbclnk
   !!   mpplnks
   !!   mpprecv
   !!   mppsend
   !!   mppscatter
   !!   mppgather
   !!   mpp_isl    : generic inteface  for :
   !!                mppisl_int , mppisl_a_int , mppisl_real, mppisl_a_real
   !!   mpp_min    : generic interface for : 
   !!                mppmin_int , mppmin_a_int , mppmin_real, mppmin_a_real
   !!   mpp_max    : generic interface for :
   !!                mppmax_int , mppmax_a_int , mppmax_real, mppmax_a_real
   !!   mpp_sum    : generic interface for :
   !!                mppsum_int , mppsum_a_int , mppsum_real, mppsum_a_real
   !!                mppsum_realdd, mppsum_a_realdd
   !!   mpp_minloc
   !!   mpp_maxloc
   !!   mppsync
   !!   mppstop
   !!   mppobc     : variant of mpp_lnk for open boundaries
   !!   mpp_ini_north
   !!   mpp_lbc_north
   !!   mpp_lbc_north_e : variant of mpp_lbc_north for extra outer halo (nsolv=4)
   !!----------------------------------------------------------------------
   !! History :
   !!        !  94 (M. Guyon, J. Escobar, M. Imbard)  Original code
   !!        !  97  (A.M. Treguier)  SHMEM additions
   !!        !  98  (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI
   !!   9.0  !  03  (J.-M. Molines, G. Madec)  F90, free form
   !!        !  04  (R. Bourdalle Badie)  isend option in mpi
   !!        !  05  (G. Madec, S. Masson)  npolj=5,6 F-point & ice cases
   !!        !  05  (R. Redler) Replacement of MPI_COMM_WORLD except for MPI_Abort
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_mpp.F90,v 1.21 2007/06/05 10:28:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!---------------------------------------------------------------------
   !! * Modules used
   USE dom_oce                    ! ocean space and time domain 
   USE in_out_manager             ! I/O manager

   IMPLICIT NONE

   PRIVATE
   PUBLIC  mynode, mpparent, mpp_isl, mpp_min, mpp_max, mpp_sum,  mpp_lbc_north
   PUBLIC  mpp_lbc_north_e, mpp_minloc, mpp_maxloc, mpp_lnk_3d, mpp_lnk_2d, mpp_lnk_3d_gather, mpp_lnk_2d_e, mpplnks
   PUBLIC  mpprecv, mppsend, mppscatter, mppgather, mppobc, mpp_ini_north, mppstop, mppsync
#if defined key_oasis3 || defined key_oasis4
   PUBLIC  mppsize, mpprank
#endif
!!DB
   PUBLIC ctl_stop

   !! * Interfaces
   !! define generic interface for these routine as they are called sometimes
   !!        with scalar arguments instead of array arguments, which causes problems
   !!        for the compilation on AIX system as well as NEC and SGI. Ok on COMPACQ

   INTERFACE mpp_isl
      MODULE PROCEDURE mppisl_a_int, mppisl_int, mppisl_a_real, mppisl_real
   END INTERFACE
   INTERFACE mpp_min
      MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
   END INTERFACE
   INTERFACE mpp_max
      MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
   END INTERFACE
   INTERFACE mpp_sum
      MODULE PROCEDURE mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real, &
                       mppsum_realdd, mppsum_a_realdd
   END INTERFACE
   INTERFACE mpp_lbc_north
      MODULE PROCEDURE mpp_lbc_north_3d, mpp_lbc_north_2d 
   END INTERFACE
  INTERFACE mpp_minloc
     MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
  END INTERFACE
  INTERFACE mpp_maxloc
     MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
  END INTERFACE


   !! * Share module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .TRUE.       !: mpp flag

   !! The processor number is a required power of two : 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,...
   INTEGER, PARAMETER ::   &
      nprocmax = 2**10     ! maximun dimension

#if defined key_mpp_mpi
   !! ========================= !!
   !!  MPI  variable definition !!
   !! ========================= !!
!$AGRIF_DO_NOT_TREAT
#  include <mpif.h>
!$AGRIF_END_DO_NOT_TREAT

   INTEGER ::   &
      mppsize,  &  ! number of process
      mpprank,  &  ! process number  [ 0 - size-1 ]
      mpi_comm_opa ! opa local communicator

   ! MD: needed to define MPI_SUMDD as a reduction operator
   INTEGER, PUBLIC :: MPI_SUMDD

   ! variables used in case of north fold condition in mpp_mpi with jpni > 1
   INTEGER ::      &       !
      ngrp_world,  &       ! group ID for the world processors
      ngrp_north,  &       ! group ID for the northern processors (to be fold)
      ncomm_north, &       ! communicator made by the processors belonging to ngrp_north
      ndim_rank_north, &   ! number of 'sea' processor in the northern line (can be /= jpni !)
      njmppmax             ! value of njmpp for the processors of the northern line
   INTEGER ::      &       !
      north_root           ! number (in the comm_opa) of proc 0 in the northern comm
   INTEGER, DIMENSION(:), ALLOCATABLE ::   &
      nrank_north          ! dimension ndim_rank_north, number of the procs belonging to ncomm_north
   CHARACTER (len=1) ::  &
      c_mpi_send = 'S'     ! type od mpi send/recieve (S=standard, B=bsend, I=isend)
   LOGICAL  ::           &
      l_isend = .FALSE.    ! isend use indicator (T if c_mpi_send='I')
!!DB: 
   CHARACTER(len=200) ::   ctmp1, ctmp2, ctmp3      !: temporary character


#elif defined key_mpp_shmem
   !! ========================= !!
   !! SHMEM variable definition !!
   !! ========================= !!
#  include  <fpvm3.h>
#  include <mpp/shmem.fh>

   CHARACTER (len=80), PARAMETER ::   simfile    = 'pvm3_ndim'   ! file name
   CHARACTER (len=47), PARAMETER ::   executable = 'opa'         ! executable name
   CHARACTER, PARAMETER ::            opaall     = ""            ! group name (old def opaall*(*))

   INTEGER, PARAMETER ::   & !! SHMEM control print
      mynode_print   = 0,  &  ! flag for print, mynode   routine
      mpprecv_print  = 0,  &  ! flag for print, mpprecv  routine
      mppsend_print  = 0,  &  ! flag for print, mppsend  routine
      mppsync_print  = 0,  &  ! flag for print, mppsync  routine
      mppsum_print   = 0,  &  ! flag for print, mpp_sum  routine
      mppisl_print   = 0,  &  ! flag for print, mpp_isl  routine
      mppmin_print   = 0,  &  ! flag for print, mpp_min  routine
      mppmax_print   = 0,  &  ! flag for print, mpp_max  routine
      mpparent_print = 0      ! flag for print, mpparent routine

   INTEGER, PARAMETER ::   & !! Variable definition
      jpvmint = 21            ! ???

   INTEGER, PARAMETER ::   & !! Maximum  dimension of array to sum on the processors
      jpmsec   = 50000,    &  ! ???
      jpmpplat =    30,    &  ! ???
      jpmppsum = MAX( jpisl*jpisl, jpmpplat*jpk, jpmsec )   ! ???

   INTEGER ::   &
      npvm_ipas ,  &  ! pvm initialization flag
      npvm_mytid,  &  ! pvm tid
      npvm_me   ,  &  ! node number [ 0 - nproc-1 ]
      npvm_nproc,  &  ! real number of nodes
      npvm_inum       ! ???
   INTEGER, DIMENSION(0:nprocmax-1) ::   &
      npvm_tids       ! tids array [ 0 - nproc-1 ]

   INTEGER ::   &
      nt3d_ipas ,  &  ! pvm initialization flag
      nt3d_mytid,  &  ! pvm tid
      nt3d_me   ,  &  ! node number [ 0 - nproc-1 ]
      nt3d_nproc      ! real number of nodes
   INTEGER, DIMENSION(0:nprocmax-1) ::   &
      nt3d_tids       ! tids array [ 0 - nproc-1 ]

   !! real sum reduction
   INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
       nrs1sync_shmem,   &  ! 
       nrs2sync_shmem
   REAL(wp), DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
       wrs1wrk_shmem,    &  !
       wrs2wrk_shmem        !
   REAL(wp), DIMENSION(jpmppsum) ::   &
       wrstab_shmem         !

   !! minimum and maximum reduction
   INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
       ni1sync_shmem,    &  ! 
       ni2sync_shmem        ! 
   REAL(wp), DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
       wi1wrk_shmem,     &  !
       wi2wrk_shmem
   REAL(wp), DIMENSION(jpmppsum) ::   &
       wintab_shmem,     &  ! 
       wi1tab_shmem,     &  ! 
       wi2tab_shmem         ! 
       
       !! value not equal zero for barotropic stream function around islands
   INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
       ni11sync_shmem,   &  !
       ni12sync_shmem,   &  !
       ni21sync_shmem,   &  !
       ni22sync_shmem       !
   REAL(wp), DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
       wi11wrk_shmem,    &  ! 
       wi12wrk_shmem,    &  !
       wi21wrk_shmem,    &  !
       wi22wrk_shmem        !
   REAL(wp), DIMENSION(jpmppsum) ::   &
       wiltab_shmem ,    &  !
       wi11tab_shmem,    &  !
       wi12tab_shmem,    &  ! 
       wi21tab_shmem,    &  ! 
       wi22tab_shmem

   INTEGER, DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
       ni11wrk_shmem,    &  !
       ni12wrk_shmem,    &  !
       ni21wrk_shmem,    &  !
       ni22wrk_shmem        !
   INTEGER, DIMENSION(jpmppsum) ::   &
       niitab_shmem ,    &  !
       ni11tab_shmem,    &  !
       ni12tab_shmem        !
   INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
       nis1sync_shmem,   &  !
       nis2sync_shmem       !
   INTEGER, DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
       nis1wrk_shmem,    &  ! 
       nis2wrk_shmem        !
   INTEGER, DIMENSION(jpmppsum) ::   &
       nistab_shmem

   !! integer sum reduction
   INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
       nil1sync_shmem,   &  !
       nil2sync_shmem       !
   INTEGER, DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
       nil1wrk_shmem,    &  !
       nil2wrk_shmem        !
   INTEGER, DIMENSION(jpmppsum) ::   &
       niltab_shmem
#endif

   REAL(wp), DIMENSION(jpi,jprecj,jpk,2,2) ::   &
       t4ns, t4sn  ! 3d message passing arrays north-south & south-north
   REAL(wp), DIMENSION(jpj,jpreci,jpk,2,2) ::   &
       t4ew, t4we  ! 3d message passing arrays east-west & west-east
   REAL(wp), DIMENSION(jpi,jprecj,jpk,2,2) ::   &
       t4p1, t4p2  ! 3d message passing arrays north fold
   REAL(wp), DIMENSION(jpi,jprecj,jpk,2) ::   &
       t3ns, t3sn  ! 3d message passing arrays north-south & south-north
   REAL(wp), DIMENSION(jpj,jpreci,jpk,2) ::   &
       t3ew, t3we  ! 3d message passing arrays east-west & west-east
   REAL(wp), DIMENSION(jpi,jprecj,jpk,2) ::   &
       t3p1, t3p2  ! 3d message passing arrays north fold
   REAL(wp), DIMENSION(jpi,jprecj,2) ::   &
       t2ns, t2sn  ! 2d message passing arrays north-south & south-north
   REAL(wp), DIMENSION(jpj,jpreci,2) ::   &
       t2ew, t2we  ! 2d message passing arrays east-west & west-east
   REAL(wp), DIMENSION(jpi,jprecj,2) ::   &
       t2p1, t2p2  ! 2d message passing arrays north fold
   REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,jprecj+jpr2dj,2) ::   &
       tr2ns, tr2sn  ! 2d message passing arrays north-south & south-north including extra outer halo
   REAL(wp), DIMENSION(1-jpr2dj:jpj+jpr2dj,jpreci+jpr2di,2) ::   &
       tr2ew, tr2we  ! 2d message passing arrays east-west & west-east including extra outer halo
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_mpp.F90,v 1.20 2007/03/07 13:26:44 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!---------------------------------------------------------------------

CONTAINS
   !!---------------------------------------------------------------------
   !!   Routine DDPDD_MPI: used by reduction operator MPI_SUMDD
   !!
   !!   Modification of original codes written by David H. Bailey
   !!   This subroutine computes ddb(i) = dda(i)+ddb(i)
   !!---------------------------------------------------------------------
      subroutine DDPDD_MPI (dda, ddb, len, itype)
      implicit none

    !! * Arguments
      INTEGER, INTENT(in)                         :: len, itype
      COMPLEX, DIMENSION(len), INTENT(in)     :: dda
      COMPLEX, DIMENSION(len), INTENT(inout)  :: ddb

    !! * Local variables
      REAL :: e, t1, t2  ! local work variables
      INTEGER :: i           ! local scalar

      DO i=1,len
   ! Compute dda + ddb using Knuth's trick.
        t1 = real(dda(i)) + real(ddb(i))
        e = t1 - real(dda(i))
        t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e))) &
             +imag(dda(i)) + imag(ddb(i))

   ! The result is t1 + t2, after normalization.
        ddb(i) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
      END DO
      END SUBROUTINE DDPDD_MPI
   !!---------------------------------------------------------------------


   FUNCTION mynode(localComm)
      !!----------------------------------------------------------------------
      !!                  ***  routine mynode  ***
      !!                    
      !! ** Purpose :   Find processor unit
      !!
      !!----------------------------------------------------------------------
#if defined key_mpp_mpi
      !! * Local variables   (MPI version)
      INTEGER ::   mynode, ierr, code
      LOGICAL ::   mpi_was_called
      INTEGER,OPTIONAL ::   localComm
      NAMELIST/nam_mpp/ c_mpi_send
      !!----------------------------------------------------------------------

!      WRITE(numout,*)
!      WRITE(numout,*) 'mynode : mpi initialisation'
!      WRITE(numout,*) '~~~~~~ '
!      WRITE(numout,*)

      ! Namelist namrun : parameters of the run
!      REWIND( numnam )
!      READ  ( numnam, nam_mpp )

!      WRITE(numout,*) '        Namelist nam_mpp'
!      WRITE(numout,*) '           mpi send type            c_mpi_send = ', c_mpi_send

! MD
! June 2006
! MPI_INIT must be called before any files are opened!
! So it's not possible to determine c_mpi_send from namelist.
! Hence this is hacked to work with a hardcoded mpi send type:

c_mpi_send = 'I'
!c_mpi_send = 'I'
!c_mpi_send = 'B'

#if defined key_agrif
      IF( Agrif_Root() ) THEN
#endif
!!bug RB : should be clean to use Agrif in coupled mode
#if ! defined key_agrif
         CALL mpi_initialized ( mpi_was_called, code )
         IF( code /= MPI_SUCCESS ) THEN
            CALL ctl_stop( ' lib_mpp: Error in routine mpi_initialized' )
            CALL mpi_abort( mpi_comm_world, code, ierr )
         ENDIF

         IF( PRESENT(localComm) .and. mpi_was_called ) THEN
            mpi_comm_opa = localComm
            SELECT CASE ( c_mpi_send )
            CASE ( 'S' )                ! Standard mpi send (blocking)
!               WRITE(numout,*) '           Standard blocking mpi send (send)'
            CASE ( 'B' )                ! Buffer mpi send (blocking)
!               WRITE(numout,*) '           Buffer blocking mpi send (bsend)'
               CALL mpi_init_opa( ierr )
            CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
!               WRITE(numout,*) '           Immediate non-blocking send (isend)'
               l_isend = .TRUE.
            CASE DEFAULT
!               WRITE(numout,cform_err)
!               WRITE(numout,*) '           bad value for c_mpi_send = ', c_mpi_send
               nstop = nstop + 1
            END SELECT
         ELSE IF ( PRESENT(localComm) .and. .not. mpi_was_called ) THEN
!            WRITE(numout,*) ' lib_mpp: You cannot provide a local communicator '
!            WRITE(numout,*) '          without calling MPI_Init before ! '
         ELSE
#endif
            SELECT CASE ( c_mpi_send )
            CASE ( 'S' )                ! Standard mpi send (blocking)
!               WRITE(numout,*) '           Standard blocking mpi send (send)'
               CALL mpi_init( ierr )
            CASE ( 'B' )                ! Buffer mpi send (blocking)
!               WRITE(numout,*) '           Buffer blocking mpi send (bsend)'
               CALL mpi_init_opa( ierr )
            CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
!               WRITE(numout,*) '           Immediate non-blocking send (isend)'
               l_isend = .TRUE.
               CALL mpi_init( ierr )
            CASE DEFAULT
!               WRITE(ctmp1,*) '           bad value for c_mpi_send = ', c_mpi_send
               CALL ctl_stop( ctmp1 )
            END SELECT

#if ! defined key_agrif
            CALL mpi_comm_dup( mpi_comm_world, mpi_comm_opa, code)
            IF( code /= MPI_SUCCESS ) THEN
               CALL ctl_stop( ' lib_mpp: Error in routine mpi_comm_dup' )
               CALL mpi_abort( mpi_comm_world, code, ierr )
            ENDIF
            !
         ENDIF
#endif
#if defined key_agrif
      ELSE
         SELECT CASE ( c_mpi_send )
         CASE ( 'S' )                ! Standard mpi send (blocking)
!            WRITE(numout,*) '           Standard blocking mpi send (send)'
         CASE ( 'B' )                ! Buffer mpi send (blocking)
!            WRITE(numout,*) '           Buffer blocking mpi send (bsend)'
         CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
!            WRITE(numout,*) '           Immediate non-blocking send (isend)'
            l_isend = .TRUE.
         CASE DEFAULT
!            WRITE(numout,cform_err)
!            WRITE(numout,*) '           bad value for c_mpi_send = ', c_mpi_send
            nstop = nstop + 1
         END SELECT
      ENDIF

      mpi_comm_opa = mpi_comm_world
#endif

      CALL mpi_comm_rank( mpi_comm_opa, mpprank, ierr )
      CALL mpi_comm_size( mpi_comm_opa, mppsize, ierr )
      mynode = mpprank

!MD: Create MPISUMDD
      call MPI_OP_CREATE(DDPDD_MPI, .TRUE., MPI_SUMDD, ierr)

#else
      !! * Local variables   (SHMEM version)
      INTEGER ::   mynode
      INTEGER ::   &
           imypid, imyhost, ji, info, iparent_tid
      !!----------------------------------------------------------------------

      IF( npvm_ipas /= nprocmax ) THEN
         !         ---   first passage in mynode
         !         -------------
         !         enroll in pvm
         !         -------------
         CALL pvmfmytid( npvm_mytid )
         IF( mynode_print /= 0 ) THEN
            WRITE(numout,*) 'mynode, npvm_ipas =', npvm_ipas, ' nprocmax=', nprocmax
            WRITE(numout,*) 'mynode, npvm_mytid=', npvm_mytid, ' after pvmfmytid'
         ENDIF

         !         ---------------------------------------------------------------
         !         find out IF i am parent or child spawned processes have parents
         !         ---------------------------------------------------------------
         CALL mpparent( iparent_tid )
         IF( mynode_print /= 0 ) THEN
            WRITE(numout,*) 'mynode, npvm_mytid=', npvm_mytid,   &
               &            ' after mpparent, npvm_tids(0) = ',   &
               &            npvm_tids(0), ' iparent_tid=', iparent_tid
         ENDIF
         IF( iparent_tid < 0 )  THEN
            WRITE(numout,*) 'mynode, npvm_mytid=', npvm_mytid,   &
               &            ' after mpparent, npvm_tids(0) = ',   &
               &            npvm_tids(0), ' iparent_tid=', iparent_tid
            npvm_tids(0) = npvm_mytid
            npvm_me = 0
            IF( jpnij > nprocmax ) THEN
               WRITE(ctmp1,*) 'npvm_mytid=', npvm_mytid, ' too great'
               CALL ctl_stop( ctmp1 )

            ELSE
               npvm_nproc = jpnij
            ENDIF

            ! -------------------------
            ! start up copies of myself
            ! -------------------------
            IF( npvm_nproc > 1 ) THEN
               DO ji = 1, npvm_nproc-1
                  npvm_tids(ji) = nt3d_tids(ji)
               END DO
               info=npvm_nproc-1
  
               IF( mynode_print /= 0 ) THEN
                  WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid,   &
                     &            ' maitre=',executable,' info=', info   &
                     &            ,' npvm_nproc=',npvm_nproc
                  WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid,   &
                     &            ' npvm_tids ',(npvm_tids(ji),ji=0,npvm_nproc-1)
               ENDIF

               ! ---------------------------
               ! multicast tids array to children
               ! ---------------------------
               CALL pvmfinitsend( pvmdefault, info )
               CALL pvmfpack ( jpvmint, npvm_nproc, 1         , 1, info )
               CALL pvmfpack ( jpvmint, npvm_tids , npvm_nproc, 1, info )
               CALL pvmfmcast( npvm_nproc-1, npvm_tids(1), 10, info )
            ENDIF
         ELSE

            ! ---------------------------------
            ! receive the tids array and set me
            ! ---------------------------------
            IF( mynode_print /= 0 )   WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid, ' pvmfrecv'
            CALL pvmfrecv( iparent_tid, 10, info )
            IF( mynode_print /= 0 )   WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid, " fin pvmfrecv"
            CALL pvmfunpack( jpvmint, npvm_nproc, 1         , 1, info )
            CALL pvmfunpack( jpvmint, npvm_tids , npvm_nproc, 1, info )
            IF( mynode_print /= 0 ) THEN
               WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid,   &
                  &            ' esclave=', executable,' info=', info,' npvm_nproc=',npvm_nproc
               WRITE(numout,*) 'mynode, npvm_mytid=', npvm_mytid,   &
                  &            'npvm_tids', ( npvm_tids(ji), ji = 0, npvm_nproc-1 )
            ENDIF
            DO ji = 0, npvm_nproc-1
               IF( npvm_mytid == npvm_tids(ji) ) npvm_me = ji
            END DO
         ENDIF

         ! ------------------------------------------------------------
         ! all nproc tasks are equal now
         ! and can address each other by tids(0) thru tids(nproc-1)
         ! for each process me => process number [0-(nproc-1)]
         ! ------------------------------------------------------------
         CALL pvmfjoingroup ( "bidon", info )
         CALL pvmfbarrier   ( "bidon", npvm_nproc, info )
         DO ji = 0, npvm_nproc-1
            IF( ji == npvm_me ) THEN
               CALL pvmfjoingroup ( opaall, npvm_inum )
               IF( npvm_inum /= npvm_me )   WRITE(numout,*) 'mynode not arrived in the good order for opaall'
            ENDIF
            CALL pvmfbarrier( "bidon", npvm_nproc, info )
         END DO
         CALL pvmfbarrier( opaall, npvm_nproc, info )
  
      ELSE
         ! ---   other passage in mynode
      ENDIF
 
      npvm_ipas = nprocmax
      mynode    = npvm_me
      imypid    = npvm_mytid
      imyhost   = npvm_tids(0)
      IF( mynode_print /= 0 ) THEN
         WRITE(numout,*)'mynode: npvm_mytid=', npvm_mytid, ' npvm_me=', npvm_me,   &
            &           ' npvm_nproc=', npvm_nproc , ' npvm_ipas=', npvm_ipas
      ENDIF
#endif
   END FUNCTION mynode


   SUBROUTINE mpparent( kparent_tid )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpparent  ***
      !!
      !! ** Purpose :   use an pvmfparent routine for T3E (key_mpp_shmem)
      !!              or  only return -1 (key_mpp_mpi)
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(inout) ::   kparent_tid      ! ???
  
#if defined key_mpp_mpi
      ! MPI version : retour -1

      kparent_tid = -1

#else
      !! * Local variables   (SHMEN onto T3E version)
      INTEGER ::   &
           it3d_my_pe, LEADZ, ji, info
  
      CALL pvmfmytid( nt3d_mytid )
      CALL pvmfgetpe( nt3d_mytid, it3d_my_pe )
      IF( mpparent_print /= 0 ) THEN
         WRITE(numout,*) 'mpparent: nt3d_mytid= ', nt3d_mytid ,' it3d_my_pe=',it3d_my_pe
      ENDIF
      IF( it3d_my_pe == 0 ) THEN
         !-----------------------------------------------------------------!
         !     process = 0 => receive other tids                           !
         !-----------------------------------------------------------------!
         kparent_tid = -1
         IF(mpparent_print /= 0 ) THEN
            WRITE(numout,*) 'mpparent, nt3d_mytid=',nt3d_mytid ,' kparent_tid=',kparent_tid
         ENDIF
         !          --- END receive dimension ---
         IF( jpnij > nprocmax ) THEN
            WRITE(ctmp1,*) 'mytid=',nt3d_mytid,' too great'
            CALL ctl_stop( ctmp1 )
         ELSE
            nt3d_nproc =  jpnij
         ENDIF
         IF( mpparent_print /= 0 ) THEN
            WRITE(numout,*) 'mpparent, nt3d_mytid=', nt3d_mytid , ' nt3d_nproc=', nt3d_nproc
         ENDIF
         !-------- receive tids from others process --------
         DO ji = 1, nt3d_nproc-1
            CALL pvmfrecv( ji , 100, info )
            CALL pvmfunpack( jpvmint, nt3d_tids(ji), 1, 1, info )
            IF( mpparent_print /= 0 ) THEN
               WRITE(numout,*) 'mpparent, nt3d_mytid=', nt3d_mytid , ' receive=', nt3d_tids(ji), ' from = ', ji
            ENDIF
         END DO
         nt3d_tids(0) = nt3d_mytid
         IF( mpparent_print /= 0 ) THEN
            WRITE(numout,*) 'mpparent, nt3d_mytid=', nt3d_mytid , ' nt3d_tids(ji) =', (nt3d_tids(ji),   &
                 ji = 0, nt3d_nproc-1 )
            WRITE(numout,*) 'mpparent, nt3d_mytid=', nt3d_mytid , ' kparent_tid=', kparent_tid
         ENDIF

      ELSE
         !!----------------------------------------------------------------!
         !     process <> 0 => send  other tids                            !
         !!----------------------------------------------------------------!
         kparent_tid = 0
         CALL pvmfinitsend( pvmdataraw, info )
         CALL pvmfpack( jpvmint, nt3d_mytid, 1, 1, info )
         CALL pvmfsend( kparent_tid, 100, info )
      ENDIF
#endif

   END SUBROUTINE mpparent

#if defined key_mpp_shmem

   SUBROUTINE mppshmem
      !!----------------------------------------------------------------------
      !!                  ***  routine mppshmem  ***
      !!
      !! ** Purpose :   SHMEM ROUTINE
      !!
      !!----------------------------------------------------------------------
      nrs1sync_shmem = SHMEM_SYNC_VALUE
      nrs2sync_shmem = SHMEM_SYNC_VALUE
      nis1sync_shmem = SHMEM_SYNC_VALUE
      nis2sync_shmem = SHMEM_SYNC_VALUE
      nil1sync_shmem = SHMEM_SYNC_VALUE
      nil2sync_shmem = SHMEM_SYNC_VALUE
      ni11sync_shmem = SHMEM_SYNC_VALUE
      ni12sync_shmem = SHMEM_SYNC_VALUE
      ni21sync_shmem = SHMEM_SYNC_VALUE
      ni22sync_shmem = SHMEM_SYNC_VALUE
      CALL barrier()
  
   END SUBROUTINE mppshmem

#endif

   SUBROUTINE mpp_lnk_3d( ptab, cd_type, psgn, cd_mpp )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_3d  ***
      !!
      !! ** Purpose :   Message passing manadgement
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask 
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors 
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!
      !! ** Action  :   ptab with update value at its periphery
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1) , INTENT( in ) ::   &
         cd_type       ! define the nature of ptab array grid-points
         !             ! = T , U , V , F , W points
         !             ! = S : T-point, north fold treatment ???
         !             ! = G : F-point, north fold treatment ???
      REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
         !             !   = -1. , the sign is changed if north fold boundary
         !             !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
         ptab          ! 3D array on which the boundary condition is applied
      CHARACTER(len=3), INTENT( in ), OPTIONAL ::    &
         cd_mpp        ! fill the overlap area only 

      !! * Local variables
      INTEGER ::   ji, jj, jk, jl                        ! dummy loop indices
      INTEGER ::   imigr, iihom, ijhom, iloc, ijt, iju   ! temporary integers
      INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
      INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
      !!----------------------------------------------------------------------

      ! 1. standard boundary treatment
      ! ------------------------------

      IF( PRESENT( cd_mpp ) ) THEN
         DO jj = nlcj+1, jpj   ! only fill extra allows last line
            ptab(1:nlci, jj, :) = ptab(1:nlci, nlej, :)
         END DO
         DO ji = nlci+1, jpi   ! only fill extra allows last column
            ptab(ji    , : , :) = ptab(nlei  , :   , :)
         END DO
      ELSE      

         !                                        ! East-West boundaries
         !                                        ! ====================
         IF( nbondi == 2 .AND.   &      ! Cyclic east-west
            &   (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
            ptab( 1 ,:,:) = ptab(jpim1,:,:)
            ptab(jpi,:,:) = ptab(  2  ,:,:)

         ELSE                           ! closed
            SELECT CASE ( cd_type )
            CASE ( 'T', 'U', 'V', 'W' )
               ptab(     1       :jpreci,:,:) = 0.e0
               ptab(nlci-jpreci+1:jpi   ,:,:) = 0.e0
            CASE ( 'F' )
               ptab(nlci-jpreci+1:jpi   ,:,:) = 0.e0
            END SELECT 
         ENDIF

         !                                        ! North-South boundaries
         !                                        ! ======================
         SELECT CASE ( cd_type )
         CASE ( 'T', 'U', 'V', 'W' )
            ptab(:,     1       :jprecj,:) = 0.e0
            ptab(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
         CASE ( 'F' )
            ptab(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
         END SELECT
     
      ENDIF

      ! 2. East and west directions exchange
      ! ------------------------------------

      ! 2.1 Read Dirichlet lateral conditions

      SELECT CASE ( nbondi )
      CASE ( -1, 0, 1 )    ! all exept 2 
         iihom = nlci-nreci
         DO jl = 1, jpreci
            t3ew(:,jl,:,1) = ptab(jpreci+jl,:,:)
            t3we(:,jl,:,1) = ptab(iihom +jl,:,:)
         END DO
      END SELECT

      ! 2.2 Migrations

#if defined key_mpp_shmem
      !! * SHMEM version

      imigr = jpreci * jpj * jpk

      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL shmem_put( t3we(1,1,1,2), t3we(1,1,1,1), imigr, noea )
      CASE ( 0 )
         CALL shmem_put( t3ew(1,1,1,2), t3ew(1,1,1,1), imigr, nowe )
         CALL shmem_put( t3we(1,1,1,2), t3we(1,1,1,1), imigr, noea )
      CASE ( 1 )
         CALL shmem_put( t3ew(1,1,1,2), t3ew(1,1,1,1), imigr, nowe )
      END SELECT

      CALL barrier()
      CALL shmem_udcflush()

#elif defined key_mpp_mpi
      !! * Local variables   (MPI version)

      imigr = jpreci * jpj * jpk

      SELECT CASE ( nbondi ) 
      CASE ( -1 )
         CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req1 )
         CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
         CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req2 )
         CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
         CALL mpprecv( 2, t3we(1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 )
         CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, t3we(1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
#endif

      ! 2.3 Write Dirichlet lateral conditions

      iihom = nlci-jpreci

      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, jpreci
            ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
         END DO
      CASE ( 0 ) 
         DO jl = 1, jpreci
            ptab(jl      ,:,:) = t3we(:,jl,:,2)
            ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
         END DO
      CASE ( 1 )
         DO jl = 1, jpreci
            ptab(jl      ,:,:) = t3we(:,jl,:,2)
         END DO
      END SELECT


      ! 3. North and south directions
      ! -----------------------------

      ! 3.1 Read Dirichlet lateral conditions

      IF( nbondj /= 2 ) THEN
         ijhom = nlcj-nrecj
         DO jl = 1, jprecj
            t3sn(:,jl,:,1) = ptab(:,ijhom +jl,:)
            t3ns(:,jl,:,1) = ptab(:,jprecj+jl,:)
         END DO
      ENDIF

      ! 3.2 Migrations

#if defined key_mpp_shmem
      !! * SHMEM version

      imigr = jprecj * jpi * jpk

      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL shmem_put( t3sn(1,1,1,2), t3sn(1,1,1,1), imigr, nono )
      CASE ( 0 )
         CALL shmem_put( t3ns(1,1,1,2), t3ns(1,1,1,1), imigr, noso )
         CALL shmem_put( t3sn(1,1,1,2), t3sn(1,1,1,1), imigr, nono )
      CASE ( 1 )
         CALL shmem_put( t3ns(1,1,1,2), t3ns(1,1,1,1), imigr, noso )
      END SELECT

      CALL barrier()
      CALL shmem_udcflush()

#elif defined key_mpp_mpi
      !! * Local variables   (MPI version)
  
      imigr=jprecj*jpi*jpk

      SELECT CASE ( nbondj )     
      CASE ( -1 )
         CALL mppsend( 4, t3sn(1,1,1,1), imigr, nono, ml_req1 )
         CALL mpprecv( 3, t3ns(1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mppsend( 3, t3ns(1,1,1,1), imigr, noso, ml_req1 )
         CALL mppsend( 4, t3sn(1,1,1,1), imigr, nono, ml_req2 )
         CALL mpprecv( 3, t3ns(1,1,1,2), imigr )
         CALL mpprecv( 4, t3sn(1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 ) 
         CALL mppsend( 3, t3ns(1,1,1,1), imigr, noso, ml_req1 )
         CALL mpprecv( 4, t3sn(1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT

#endif

      ! 3.3 Write Dirichlet lateral conditions

      ijhom = nlcj-jprecj

      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, jprecj
            ptab(:,ijhom+jl,:) = t3ns(:,jl,:,2)
         END DO
      CASE ( 0 ) 
         DO jl = 1, jprecj
            ptab(:,jl      ,:) = t3sn(:,jl,:,2)
            ptab(:,ijhom+jl,:) = t3ns(:,jl,:,2)
         END DO
      CASE ( 1 )
         DO jl = 1, jprecj
            ptab(:,jl,:) = t3sn(:,jl,:,2)
         END DO
      END SELECT


      ! 4. north fold treatment
      ! -----------------------

      IF (PRESENT(cd_mpp)) THEN
         ! No north fold treatment (it is assumed to be already OK)
     
      ELSE      

      ! 4.1 treatment without exchange (jpni odd)
      !     T-point pivot  

      SELECT CASE ( jpni )

      CASE ( 1 )  ! only one proc along I, no mpp exchange

         SELECT CASE ( npolj )
  
         CASE ( 3 , 4 )    ! T pivot
            iloc = jpiglo - 2 * ( nimpp - 1 )

            SELECT CASE ( cd_type )

            CASE ( 'T' , 'S', 'W' )
               DO jk = 1, jpk
                  DO ji = 2, nlci
                     ijt=iloc-ji+2
                     ptab(ji,nlcj,jk) = psgn * ptab(ijt,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci
                     ijt=iloc-ji+2
                     ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-1,jk)
                  END DO
               END DO
          
            CASE ( 'U' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji+1
                     ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2, nlci-1
                     iju=iloc-ji+1
                     ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'V' )
               DO jk = 1, jpk
                  DO ji = 2, nlci
                     ijt=iloc-ji+2
                     ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-2,jk)
                     ptab(ji,nlcj  ,jk) = psgn * ptab(ijt,nlcj-3,jk)
                  END DO
               END DO

            CASE ( 'F', 'G' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji+1
                     ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-2,jk)
                     ptab(ji,nlcj  ,jk) = psgn * ptab(iju,nlcj-3,jk)
                  END DO
               END DO
  
          END SELECT
       
         CASE ( 5 , 6 ) ! F pivot
            iloc=jpiglo-2*(nimpp-1)
  
            SELECT CASE ( cd_type )

            CASE ( 'T' , 'S', 'W' )
               DO jk = 1, jpk
                  DO ji = 1, nlci
                     ijt=iloc-ji+1
                     ptab(ji,nlcj,jk) = psgn * ptab(ijt,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'U' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji
                     ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'V' )
               DO jk = 1, jpk
                  DO ji = 1, nlci
                     ijt=iloc-ji+1
                     ptab(ji,nlcj  ,jk) = psgn * ptab(ijt,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci
                     ijt=iloc-ji+1
                     ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'F', 'G' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji
                     ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci-1
                     iju=iloc-ji
                     ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-1,jk)
                  END DO
               END DO
            END SELECT  ! cd_type

         END SELECT     !  npolj
  
      CASE DEFAULT ! more than 1 proc along I
         IF ( npolj /= 0 ) CALL mpp_lbc_north (ptab, cd_type, psgn)  ! only for northern procs.

      END SELECT ! jpni 

      ENDIF
      

      ! 5. East and west directions exchange
      ! ------------------------------------

      SELECT CASE ( npolj )

      CASE ( 3, 4, 5, 6 )

         ! 5.1 Read Dirichlet lateral conditions

         SELECT CASE ( nbondi )

         CASE ( -1, 0, 1 )
            iihom = nlci-nreci
            DO jl = 1, jpreci
               t3ew(:,jl,:,1) = ptab(jpreci+jl,:,:)
               t3we(:,jl,:,1) = ptab(iihom +jl,:,:)
            END DO

         END SELECT

         ! 5.2 Migrations

#if defined key_mpp_shmem
         !! SHMEM version

         imigr = jpreci * jpj * jpk

         SELECT CASE ( nbondi )
         CASE ( -1 )
            CALL shmem_put( t3we(1,1,1,2), t3we(1,1,1,1), imigr, noea )
         CASE ( 0 )
            CALL shmem_put( t3ew(1,1,1,2), t3ew(1,1,1,1), imigr, nowe )
            CALL shmem_put( t3we(1,1,1,2), t3we(1,1,1,1), imigr, noea )
         CASE ( 1 )
            CALL shmem_put( t3ew(1,1,1,2), t3ew(1,1,1,1), imigr, nowe )
         END SELECT

         CALL barrier()
         CALL shmem_udcflush()

#elif defined key_mpp_mpi
         !! MPI version

         imigr=jpreci*jpj*jpk
  
         SELECT CASE ( nbondi )
         CASE ( -1 )
            CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req1 )
            CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         CASE ( 0 )
            CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
            CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req2 )
            CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
            CALL mpprecv( 2, t3we(1,1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
            IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
         CASE ( 1 )
            CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
            CALL mpprecv( 2, t3we(1,1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         END SELECT
#endif

         ! 5.3 Write Dirichlet lateral conditions

         iihom = nlci-jpreci

         SELECT CASE ( nbondi)
         CASE ( -1 )
            DO jl = 1, jpreci
               ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
            END DO
         CASE ( 0 ) 
            DO jl = 1, jpreci
               ptab(jl      ,:,:) = t3we(:,jl,:,2)
               ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
            END DO
         CASE ( 1 )
            DO jl = 1, jpreci
               ptab(jl      ,:,:) = t3we(:,jl,:,2)
            END DO
         END SELECT

      END SELECT    ! npolj 

   END SUBROUTINE mpp_lnk_3d


   SUBROUTINE mpp_lnk_2d( pt2d, cd_type, psgn, cd_mpp )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d  ***
      !!                  
      !! ** Purpose :   Message passing manadgement for 2d array
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask 
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors 
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1) , INTENT( in ) ::   &
         cd_type       ! define the nature of pt2d array grid-points
         !             !  = T , U , V , F , W 
         !             !  = S : T-point, north fold treatment
         !             !  = G : F-point, north fold treatment
         !             !  = I : sea-ice velocity at F-point with index shift
      REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
         !             !   = -1. , the sign is changed if north fold boundary
         !             !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
         pt2d          ! 2D array on which the boundary condition is applied
      CHARACTER(len=3), INTENT( in ), OPTIONAL ::    &
         cd_mpp        ! fill the overlap area only 

      !! * Local variables
      INTEGER  ::   ji, jj, jl      ! dummy loop indices
      INTEGER  ::   &
         imigr, iihom, ijhom,    &  ! temporary integers
         iloc, ijt, iju             !    "          "
      INTEGER  ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
      INTEGER  ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
      !!----------------------------------------------------------------------

      ! 1. standard boundary treatment
      ! ------------------------------
      IF (PRESENT(cd_mpp)) THEN
         DO jj = nlcj+1, jpj   ! only fill extra allows last line
            pt2d(1:nlci, jj) = pt2d(1:nlci, nlej)
         END DO
         DO ji = nlci+1, jpi   ! only fill extra allows last column
            pt2d(ji    , : ) = pt2d(nlei  , :   )
         END DO     
      ELSE      

         !                                        ! East-West boundaries
         !                                        ! ====================
         IF( nbondi == 2 .AND.   &      ! Cyclic east-west
            &    (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
            pt2d( 1 ,:) = pt2d(jpim1,:)
            pt2d(jpi,:) = pt2d(  2  ,:)

         ELSE                           ! ... closed
            SELECT CASE ( cd_type )
            CASE ( 'T', 'U', 'V', 'W' , 'I' )
               pt2d(     1       :jpreci,:) = 0.e0
               pt2d(nlci-jpreci+1:jpi   ,:) = 0.e0
            CASE ( 'F' )
               pt2d(nlci-jpreci+1:jpi   ,:) = 0.e0
            END SELECT
         ENDIF

         !                                        ! North-South boundaries
         !                                        ! ======================
         SELECT CASE ( cd_type )
         CASE ( 'T', 'U', 'V', 'W' , 'I' )
            pt2d(:,     1       :jprecj) = 0.e0
            pt2d(:,nlcj-jprecj+1:jpj   ) = 0.e0
         CASE ( 'F' )
            pt2d(:,nlcj-jprecj+1:jpj   ) = 0.e0
         END SELECT

      ENDIF


      ! 2. East and west directions
      ! ---------------------------

      ! 2.1 Read Dirichlet lateral conditions

      SELECT CASE ( nbondi )
      CASE ( -1, 0, 1 )    ! all except 2
         iihom = nlci-nreci
         DO jl = 1, jpreci
            t2ew(:,jl,1) = pt2d(jpreci+jl,:)
            t2we(:,jl,1) = pt2d(iihom +jl,:)
         END DO
      END SELECT

      ! 2.2 Migrations

#if defined key_mpp_shmem
      !! * SHMEM version

      imigr = jpreci * jpj

      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr, noea )
      CASE ( 0 )
         CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr, nowe )
         CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr, noea )
      CASE ( 1 )
         CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr, nowe )
      END SELECT

      CALL barrier()
      CALL shmem_udcflush()

#elif defined key_mpp_mpi
      !! * MPI version

      imigr = jpreci * jpj

      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req1 )
         CALL mpprecv( 1, t2ew(1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
         CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req2 )
         CALL mpprecv( 1, t2ew(1,1,2), imigr )
         CALL mpprecv( 2, t2we(1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, t2we(1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT

#endif

      ! 2.3 Write Dirichlet lateral conditions

      iihom = nlci - jpreci
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, jpreci
            pt2d(iihom+jl,:) = t2ew(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, jpreci
            pt2d(jl      ,:) = t2we(:,jl,2)
            pt2d(iihom+jl,:) = t2ew(:,jl,2)
         END DO
      CASE ( 1 )
         DO jl = 1, jpreci
            pt2d(jl      ,:) = t2we(:,jl,2)
         END DO
      END SELECT


      ! 3. North and south directions
      ! -----------------------------

      ! 3.1 Read Dirichlet lateral conditions

      IF( nbondj /= 2 ) THEN
         ijhom = nlcj-nrecj
         DO jl = 1, jprecj
            t2sn(:,jl,1) = pt2d(:,ijhom +jl)
            t2ns(:,jl,1) = pt2d(:,jprecj+jl)
         END DO
      ENDIF

      ! 3.2 Migrations

#if defined key_mpp_shmem
      !! * SHMEM version

      imigr = jprecj * jpi

      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL shmem_put( t2sn(1,1,2), t2sn(1,1,1), imigr, nono )
      CASE ( 0 )
         CALL shmem_put( t2ns(1,1,2), t2ns(1,1,1), imigr, noso )
         CALL shmem_put( t2sn(1,1,2), t2sn(1,1,1), imigr, nono )
      CASE ( 1 )
         CALL shmem_put( t2ns(1,1,2), t2ns(1,1,1), imigr, noso )
      END SELECT 
      CALL barrier()
      CALL shmem_udcflush()

#elif defined key_mpp_mpi
      !! * MPI version

      imigr = jprecj * jpi

      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, t2sn(1,1,1), imigr, nono, ml_req1 )
         CALL mpprecv( 3, t2ns(1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 3, t2ns(1,1,1), imigr, noso, ml_req1 )
         CALL mppsend( 4, t2sn(1,1,1), imigr, nono, ml_req2 )
         CALL mpprecv( 3, t2ns(1,1,2), imigr )
         CALL mpprecv( 4, t2sn(1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 3, t2ns(1,1,1), imigr, noso, ml_req1 )
         CALL mpprecv( 4, t2sn(1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
  
#endif

      ! 3.3 Write Dirichlet lateral conditions

      ijhom = nlcj - jprecj

      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, jprecj
            pt2d(:,ijhom+jl) = t2ns(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, jprecj
            pt2d(:,jl      ) = t2sn(:,jl,2)
            pt2d(:,ijhom+jl) = t2ns(:,jl,2)
         END DO
      CASE ( 1 ) 
         DO jl = 1, jprecj
            pt2d(:,jl      ) = t2sn(:,jl,2)
         END DO
      END SELECT 
  

      ! 4. north fold treatment
      ! -----------------------
  
      IF (PRESENT(cd_mpp)) THEN
         ! No north fold treatment (it is assumed to be already OK)
     
      ELSE      

      ! 4.1 treatment without exchange (jpni odd)
      
      SELECT CASE ( jpni )
  
      CASE ( 1 ) ! only one proc along I, no mpp exchange
  
         SELECT CASE ( npolj )
  
         CASE ( 3 , 4 )   !  T pivot
            iloc = jpiglo - 2 * ( nimpp - 1 )
  
            SELECT CASE ( cd_type )
  
            CASE ( 'T' , 'S', 'W' )
               DO ji = 2, nlci
                  ijt=iloc-ji+2
                  pt2d(ji,nlcj) = psgn * pt2d(ijt,nlcj-2)
               END DO
               DO ji = nlci/2+1, nlci
                  ijt=iloc-ji+2
                  pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
               END DO
  
            CASE ( 'U' )
               DO ji = 1, nlci-1
                  iju=iloc-ji+1
                  pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-2)
               END DO
               DO ji = nlci/2, nlci-1
                  iju=iloc-ji+1
                  pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
               END DO
  
            CASE ( 'V' )
               DO ji = 2, nlci
                  ijt=iloc-ji+2
                  pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-2)
                  pt2d(ji,nlcj  ) = psgn * pt2d(ijt,nlcj-3)
               END DO
  
            CASE ( 'F', 'G' )
               DO ji = 1, nlci-1
                  iju=iloc-ji+1
                  pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-2)
                  pt2d(ji,nlcj  ) = psgn * pt2d(iju,nlcj-3)
               END DO
  
            CASE ( 'I' )                                  ! ice U-V point
               pt2d(2,nlcj) = psgn * pt2d(3,nlcj-1)
               DO ji = 3, nlci
                  iju = iloc - ji + 3
                  pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-1)
               END DO
  
            END SELECT
  
         CASE ( 5 , 6 )                 ! F pivot
            iloc=jpiglo-2*(nimpp-1)
  
            SELECT CASE (cd_type )
  
            CASE ( 'T', 'S', 'W' )
               DO ji = 1, nlci
                  ijt=iloc-ji+1
                  pt2d(ji,nlcj) = psgn * pt2d(ijt,nlcj-1)
               END DO
  
            CASE ( 'U' )
               DO ji = 1, nlci-1
                  iju=iloc-ji
                  pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-1)
               END DO

            CASE ( 'V' )
               DO ji = 1, nlci
                  ijt=iloc-ji+1
                  pt2d(ji,nlcj  ) = psgn * pt2d(ijt,nlcj-2)
               END DO
               DO ji = nlci/2+1, nlci
                  ijt=iloc-ji+1
                  pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
               END DO
  
            CASE ( 'F', 'G' )
               DO ji = 1, nlci-1
                  iju=iloc-ji
                  pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-2)
               END DO
               DO ji = nlci/2+1, nlci-1
                  iju=iloc-ji
                  pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
               END DO
  
            CASE ( 'I' )                                  ! ice U-V point
               pt2d( 2 ,nlcj) = 0.e0
               DO ji = 2 , nlci-1
                  ijt = iloc - ji + 2
                  pt2d(ji,nlcj)= 0.5 * ( pt2d(ji,nlcj-1) + psgn * pt2d(ijt,nlcj-1) )
               END DO
  
            END SELECT   ! cd_type
  
         END SELECT   ! npolj

      CASE DEFAULT   ! more than 1 proc along I
         IF( npolj /= 0 )   CALL mpp_lbc_north( pt2d, cd_type, psgn )   ! only for northern procs.

      END SELECT   ! jpni

      ENDIF

      ! 5. East and west directions
      ! ---------------------------

      SELECT CASE ( npolj )

      CASE ( 3, 4, 5, 6 )

         ! 5.1 Read Dirichlet lateral conditions

         SELECT CASE ( nbondi )
         CASE ( -1, 0, 1 )
            iihom = nlci-nreci
            DO jl = 1, jpreci
               DO jj = 1, jpj
                  t2ew(jj,jl,1) = pt2d(jpreci+jl,jj)
                  t2we(jj,jl,1) = pt2d(iihom +jl,jj)
               END DO
            END DO
         END SELECT

         ! 5.2 Migrations

#if defined key_mpp_shmem
         !! * SHMEM version

         imigr=jpreci*jpj

         SELECT CASE ( nbondi )
         CASE ( -1 )
            CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr, noea )
         CASE ( 0 )
            CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr, nowe )
            CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr, noea )
         CASE ( 1 )
            CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr, nowe )
         END SELECT

         CALL barrier()
         CALL shmem_udcflush()
  
#elif defined key_mpp_mpi
         !! * MPI version
  
         imigr=jpreci*jpj
  
         SELECT CASE ( nbondi )
         CASE ( -1 )
            CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req1 )
            CALL mpprecv( 1, t2ew(1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         CASE ( 0 )
            CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
            CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req2 )
            CALL mpprecv( 1, t2ew(1,1,2), imigr )
            CALL mpprecv( 2, t2we(1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
            IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
         CASE ( 1 )
            CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
            CALL mpprecv( 2, t2we(1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         END SELECT 
#endif

         ! 5.3 Write Dirichlet lateral conditions
  
         iihom = nlci - jpreci
  
         SELECT CASE ( nbondi )
         CASE ( -1 )
            DO jl = 1, jpreci
               pt2d(iihom+jl,:) = t2ew(:,jl,2)
            END DO
         CASE ( 0 )
            DO jl = 1, jpreci
               pt2d(jl      ,:) = t2we(:,jl,2)
               pt2d(iihom+jl,:) = t2ew(:,jl,2)
            END DO
         CASE ( 1 )
            DO jl = 1, jpreci
               pt2d(jl,:) = t2we(:,jl,2)
            END DO
         END SELECT 
  
      END SELECT   ! npolj
  
   END SUBROUTINE mpp_lnk_2d


   SUBROUTINE mpp_lnk_3d_gather( ptab1, cd_type1, ptab2, cd_type2, psgn )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_3d_gather  ***
      !!
      !! ** Purpose :   Message passing manadgement for two 3D arrays
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask 
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors 
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!
      !! ** Action  :   ptab1 and ptab2  with update value at its periphery
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1) , INTENT( in ) ::   &
         cd_type1, cd_type2       ! define the nature of ptab array grid-points
         !                        ! = T , U , V , F , W points
         !                        ! = S : T-point, north fold treatment ???
         !                        ! = G : F-point, north fold treatment ???
      REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
         !             !   = -1. , the sign is changed if north fold boundary
         !             !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
         ptab1, ptab2             ! 3D array on which the boundary condition is applied

      !! * Local variables
      INTEGER ::   ji, jk, jl   ! dummy loop indices
      INTEGER ::   imigr, iihom, ijhom, iloc, ijt, iju   ! temporary integers
      INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
      INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
      !!----------------------------------------------------------------------

      ! 1. standard boundary treatment
      ! ------------------------------
      !                                        ! East-West boundaries
      !                                        ! ====================
      IF( nbondi == 2 .AND.   &      ! Cyclic east-west
         &   (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
         ptab1( 1 ,:,:) = ptab1(jpim1,:,:)
         ptab1(jpi,:,:) = ptab1(  2  ,:,:)
         ptab2( 1 ,:,:) = ptab2(jpim1,:,:)
         ptab2(jpi,:,:) = ptab2(  2  ,:,:)

      ELSE                           ! closed
         SELECT CASE ( cd_type1 )
         CASE ( 'T', 'U', 'V', 'W' )
            ptab1(     1       :jpreci,:,:) = 0.e0
            ptab1(nlci-jpreci+1:jpi   ,:,:) = 0.e0
         CASE ( 'F' )
            ptab1(nlci-jpreci+1:jpi   ,:,:) = 0.e0
         END SELECT 
         SELECT CASE ( cd_type2 )
         CASE ( 'T', 'U', 'V', 'W' )
            ptab2(     1       :jpreci,:,:) = 0.e0
            ptab2(nlci-jpreci+1:jpi   ,:,:) = 0.e0
         CASE ( 'F' )
            ptab2(nlci-jpreci+1:jpi   ,:,:) = 0.e0
         END SELECT 
      ENDIF

      !                                        ! North-South boundaries
      !                                        ! ======================
      SELECT CASE ( cd_type1 )
      CASE ( 'T', 'U', 'V', 'W' )
         ptab1(:,     1       :jprecj,:) = 0.e0
         ptab1(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
      CASE ( 'F' )
         ptab1(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
      END SELECT

      SELECT CASE ( cd_type2 )
      CASE ( 'T', 'U', 'V', 'W' )
         ptab2(:,     1       :jprecj,:) = 0.e0
         ptab2(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
      CASE ( 'F' )
         ptab2(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
      END SELECT


      ! 2. East and west directions exchange
      ! ------------------------------------

      ! 2.1 Read Dirichlet lateral conditions

      SELECT CASE ( nbondi )
      CASE ( -1, 0, 1 )    ! all exept 2 
         iihom = nlci-nreci
         DO jl = 1, jpreci
            t4ew(:,jl,:,1,1) = ptab1(jpreci+jl,:,:)
            t4we(:,jl,:,1,1) = ptab1(iihom +jl,:,:)
            t4ew(:,jl,:,2,1) = ptab2(jpreci+jl,:,:)
            t4we(:,jl,:,2,1) = ptab2(iihom +jl,:,:)
         END DO
      END SELECT

      ! 2.2 Migrations

#if defined key_mpp_shmem
      !! * SHMEM version

      imigr = jpreci * jpj * jpk *2

      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL shmem_put( t4we(1,1,1,1,2), t4we(1,1,1,1,1), imigr, noea )
      CASE ( 0 )
         CALL shmem_put( t4ew(1,1,1,1,2), t4ew(1,1,1,1,1), imigr, nowe )
         CALL shmem_put( t4we(1,1,1,1,2), t4we(1,1,1,1,1), imigr, noea )
      CASE ( 1 )
         CALL shmem_put( t4ew(1,1,1,1,2), t4ew(1,1,1,1,1), imigr, nowe )
      END SELECT

      CALL barrier()
      CALL shmem_udcflush()

#elif defined key_mpp_mpi
      !! * Local variables   (MPI version)

      imigr = jpreci * jpj * jpk *2

      SELECT CASE ( nbondi ) 
      CASE ( -1 )
         CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req1 )
         CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
         CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req2 )
         CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
         CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 )
         CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
#endif

      ! 2.3 Write Dirichlet lateral conditions

      iihom = nlci-jpreci

      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, jpreci
            ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
            ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
         END DO
      CASE ( 0 ) 
         DO jl = 1, jpreci
            ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
            ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
            ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
            ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
         END DO
      CASE ( 1 )
         DO jl = 1, jpreci
            ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
            ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
         END DO
      END SELECT


      ! 3. North and south directions
      ! -----------------------------

      ! 3.1 Read Dirichlet lateral conditions

      IF( nbondj /= 2 ) THEN
         ijhom = nlcj-nrecj
         DO jl = 1, jprecj
            t4sn(:,jl,:,1,1) = ptab1(:,ijhom +jl,:)
            t4ns(:,jl,:,1,1) = ptab1(:,jprecj+jl,:)
            t4sn(:,jl,:,2,1) = ptab2(:,ijhom +jl,:)
            t4ns(:,jl,:,2,1) = ptab2(:,jprecj+jl,:)
         END DO
      ENDIF

      ! 3.2 Migrations

#if defined key_mpp_shmem
      !! * SHMEM version

      imigr = jprecj * jpi * jpk * 2

      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL shmem_put( t4sn(1,1,1,1,2), t4sn(1,1,1,1,1), imigr, nono )
      CASE ( 0 )
         CALL shmem_put( t4ns(1,1,1,1,2), t4ns(1,1,1,1,1), imigr, noso )
         CALL shmem_put( t4sn(1,1,1,1,2), t4sn(1,1,1,1,1), imigr, nono )
      CASE ( 1 )
         CALL shmem_put( t4ns(1,1,1,1,2), t4ns(1,1,1,1;,1), imigr, noso )
      END SELECT

      CALL barrier()
      CALL shmem_udcflush()

#elif defined key_mpp_mpi
      !! * Local variables   (MPI version)
  
      imigr=jprecj * jpi * jpk * 2

      SELECT CASE ( nbondj )     
      CASE ( -1 )
         CALL mppsend( 4, t4sn(1,1,1,1,1), imigr, nono, ml_req1 )
         CALL mpprecv( 3, t4ns(1,1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mppsend( 3, t4ns(1,1,1,1,1), imigr, noso, ml_req1 )
         CALL mppsend( 4, t4sn(1,1,1,1,1), imigr, nono, ml_req2 )
         CALL mpprecv( 3, t4ns(1,1,1,1,2), imigr )
         CALL mpprecv( 4, t4sn(1,1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 ) 
         CALL mppsend( 3, t4ns(1,1,1,1,1), imigr, noso, ml_req1 )
         CALL mpprecv( 4, t4sn(1,1,1,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT

#endif

      ! 3.3 Write Dirichlet lateral conditions

      ijhom = nlcj-jprecj

      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, jprecj
            ptab1(:,ijhom+jl,:) = t4ns(:,jl,:,1,2)
            ptab2(:,ijhom+jl,:) = t4ns(:,jl,:,2,2)
         END DO
      CASE ( 0 ) 
         DO jl = 1, jprecj
            ptab1(:,jl      ,:) = t4sn(:,jl,:,1,2)
            ptab1(:,ijhom+jl,:) = t4ns(:,jl,:,1,2)
            ptab2(:,jl      ,:) = t4sn(:,jl,:,2,2)
            ptab2(:,ijhom+jl,:) = t4ns(:,jl,:,2,2)
         END DO
      CASE ( 1 )
         DO jl = 1, jprecj
            ptab1(:,jl,:) = t4sn(:,jl,:,1,2)
            ptab2(:,jl,:) = t4sn(:,jl,:,2,2)
         END DO
      END SELECT


      ! 4. north fold treatment
      ! -----------------------

      ! 4.1 treatment without exchange (jpni odd)
      !     T-point pivot  

      SELECT CASE ( jpni )

      CASE ( 1 )  ! only one proc along I, no mpp exchange

      SELECT CASE ( npolj )
  
         CASE ( 3 , 4 )    ! T pivot
            iloc = jpiglo - 2 * ( nimpp - 1 )

            SELECT CASE ( cd_type1 )

            CASE ( 'T' , 'S', 'W' )
               DO jk = 1, jpk
                  DO ji = 2, nlci
                     ijt=iloc-ji+2
                     ptab1(ji,nlcj,jk) = psgn * ptab1(ijt,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci
                     ijt=iloc-ji+2
                     ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-1,jk)
                  END DO
               END DO
          
            CASE ( 'U' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji+1
                     ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2, nlci-1
                     iju=iloc-ji+1
                     ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'V' )
               DO jk = 1, jpk
                  DO ji = 2, nlci
                     ijt=iloc-ji+2
                     ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-2,jk)
                     ptab1(ji,nlcj  ,jk) = psgn * ptab1(ijt,nlcj-3,jk)
                  END DO
               END DO

            CASE ( 'F', 'G' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji+1
                     ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-2,jk)
                     ptab1(ji,nlcj  ,jk) = psgn * ptab1(iju,nlcj-3,jk)
                  END DO
               END DO
  
            END SELECT
            
            SELECT CASE ( cd_type2 )

            CASE ( 'T' , 'S', 'W' )
               DO jk = 1, jpk
                  DO ji = 2, nlci
                     ijt=iloc-ji+2
                     ptab2(ji,nlcj,jk) = psgn * ptab2(ijt,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci
                     ijt=iloc-ji+2
                     ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-1,jk)
                  END DO
               END DO
          
            CASE ( 'U' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji+1
                     ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2, nlci-1
                     iju=iloc-ji+1
                     ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'V' )
               DO jk = 1, jpk
                  DO ji = 2, nlci
                     ijt=iloc-ji+2
                     ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-2,jk)
                     ptab2(ji,nlcj  ,jk) = psgn * ptab2(ijt,nlcj-3,jk)
                  END DO
               END DO

            CASE ( 'F', 'G' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji+1
                     ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-2,jk)
                     ptab2(ji,nlcj  ,jk) = psgn * ptab2(iju,nlcj-3,jk)
                  END DO
               END DO
  
          END SELECT
       
         CASE ( 5 , 6 ) ! F pivot
            iloc=jpiglo-2*(nimpp-1)
  
            SELECT CASE ( cd_type1 )

            CASE ( 'T' , 'S', 'W' )
               DO jk = 1, jpk
                  DO ji = 1, nlci
                     ijt=iloc-ji+1
                     ptab1(ji,nlcj,jk) = psgn * ptab1(ijt,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'U' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji
                     ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'V' )
               DO jk = 1, jpk
                  DO ji = 1, nlci
                     ijt=iloc-ji+1
                     ptab1(ji,nlcj  ,jk) = psgn * ptab1(ijt,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci
                     ijt=iloc-ji+1
                     ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'F', 'G' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji
                     ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci-1
                     iju=iloc-ji
                     ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-1,jk)
                  END DO
               END DO
            END SELECT  ! cd_type1

            SELECT CASE ( cd_type2 )

            CASE ( 'T' , 'S', 'W' )
               DO jk = 1, jpk
                  DO ji = 1, nlci
                     ijt=iloc-ji+1
                     ptab2(ji,nlcj,jk) = psgn * ptab2(ijt,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'U' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji
                     ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'V' )
               DO jk = 1, jpk
                  DO ji = 1, nlci
                     ijt=iloc-ji+1
                     ptab2(ji,nlcj  ,jk) = psgn * ptab2(ijt,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci
                     ijt=iloc-ji+1
                     ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-1,jk)
                  END DO
               END DO

            CASE ( 'F', 'G' )
               DO jk = 1, jpk
                  DO ji = 1, nlci-1
                     iju=iloc-ji
                     ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-2,jk)
                  END DO
                  DO ji = nlci/2+1, nlci-1
                     iju=iloc-ji
                     ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-1,jk)
                  END DO
               END DO

            END SELECT  ! cd_type2

         END SELECT     !  npolj
  
      CASE DEFAULT ! more than 1 proc along I
         IF ( npolj /= 0 ) THEN
            CALL mpp_lbc_north (ptab1, cd_type1, psgn)  ! only for northern procs.
            CALL mpp_lbc_north (ptab2, cd_type2, psgn)  ! only for northern procs.
         ENDIF

      END SELECT ! jpni 


      ! 5. East and west directions exchange
      ! ------------------------------------

      SELECT CASE ( npolj )

      CASE ( 3, 4, 5, 6 )

         ! 5.1 Read Dirichlet lateral conditions

         SELECT CASE ( nbondi )

         CASE ( -1, 0, 1 )
            iihom = nlci-nreci
            DO jl = 1, jpreci
               t4ew(:,jl,:,1,1) = ptab1(jpreci+jl,:,:)
               t4we(:,jl,:,1,1) = ptab1(iihom +jl,:,:)
               t4ew(:,jl,:,2,1) = ptab2(jpreci+jl,:,:)
               t4we(:,jl,:,2,1) = ptab2(iihom +jl,:,:)
            END DO

         END SELECT

         ! 5.2 Migrations

#if defined key_mpp_shmem
         !! SHMEM version

         imigr = jpreci * jpj * jpk * 2

         SELECT CASE ( nbondi )
         CASE ( -1 )
            CALL shmem_put( t4we(1,1,1,1,2), t4we(1,1,1,1,1), imigr, noea )
         CASE ( 0 )
            CALL shmem_put( t4ew(1,1,1,1,2), t4ew(1,1,1,1,1), imigr, nowe )
            CALL shmem_put( t4we(1,1,1,1,2), t4we(1,1,1,1,1), imigr, noea )
         CASE ( 1 )
            CALL shmem_put( t4ew(1,1,1,1,2), t4ew(1,1,1,1,1), imigr, nowe )
         END SELECT

         CALL barrier()
         CALL shmem_udcflush()

#elif defined key_mpp_mpi
         !! MPI version

         imigr = jpreci * jpj * jpk * 2
  
         SELECT CASE ( nbondi )
         CASE ( -1 )
            CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req1 )
            CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         CASE ( 0 )
            CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
            CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req2 )
            CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
            CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
            IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
         CASE ( 1 )
            CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
            CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         END SELECT
#endif

         ! 5.3 Write Dirichlet lateral conditions

         iihom = nlci-jpreci

         SELECT CASE ( nbondi)
         CASE ( -1 )
            DO jl = 1, jpreci
               ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
               ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
            END DO
         CASE ( 0 ) 
            DO jl = 1, jpreci
               ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
               ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
               ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
               ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
            END DO
         CASE ( 1 )
            DO jl = 1, jpreci
               ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
               ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
            END DO
         END SELECT

      END SELECT    ! npolj 

   END SUBROUTINE mpp_lnk_3d_gather


   SUBROUTINE mpp_lnk_2d_e( pt2d, cd_type, psgn )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_e  ***
      !!                  
      !! ** Purpose :   Message passing manadgement for 2d array (with halo)
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask 
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    jpr2di : number of rows for extra outer halo
      !!                    jpr2dj : number of columns for extra outer halo
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors 
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!   
      !! History :
      !!       
      !!   9.0  !  05-09  (R. Benshila, G. Madec)  original code
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1) , INTENT( in ) ::   &
         cd_type       ! define the nature of pt2d array grid-points
         !             !  = T , U , V , F , W 
         !             !  = S : T-point, north fold treatment
         !             !  = G : F-point, north fold treatment
         !             !  = I : sea-ice velocity at F-point with index shift
      REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
         !             !   = -1. , the sign is changed if north fold boundary
         !             !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj), INTENT( inout ) ::   &
         pt2d          ! 2D array on which the boundary condition is applied

      !! * Local variables
      INTEGER  ::   ji, jl      ! dummy loop indices
      INTEGER  ::   &
         imigr, iihom, ijhom,    &  ! temporary integers
         iloc, ijt, iju             !    "          "
      INTEGER  ::   &
         ipreci, iprecj             ! temporary integers
      INTEGER  ::   ml_req1, ml_req2, ml_err     ! for isend
      INTEGER  ::   ml_stat(MPI_STATUS_SIZE)     ! for isend
     !!---------------------------------------------------------------------

      ! take into account outer extra 2D overlap area
      ipreci = jpreci + jpr2di
      iprecj = jprecj + jpr2dj


      ! 1. standard boundary treatment
      ! ------------------------------

      !                                        ! East-West boundaries
      !                                        ! ====================
      IF( nbondi == 2 .AND.   &      ! Cyclic east-west
         &    (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
         pt2d(1-jpr2di:     1    ,:) = pt2d(jpim1-jpr2di:  jpim1 ,:)
         pt2d(   jpi  :jpi+jpr2di,:) = pt2d(     2      :2+jpr2di,:)

      ELSE                           ! ... closed
         SELECT CASE ( cd_type )
         CASE ( 'T', 'U', 'V', 'W' , 'I' )
            pt2d(  1-jpr2di   :jpreci    ,:) = 0.e0
            pt2d(nlci-jpreci+1:jpi+jpr2di,:) = 0.e0
         CASE ( 'F' )
            pt2d(nlci-jpreci+1:jpi+jpr2di,:) = 0.e0
         END SELECT
      ENDIF

      !                                        ! North-South boundaries
      !                                        ! ======================
      SELECT CASE ( cd_type )
      CASE ( 'T', 'U', 'V', 'W' , 'I' )
         pt2d(:,  1-jpr2dj   :  jprecj  ) = 0.e0
         pt2d(:,nlcj-jprecj+1:jpj+jpr2dj) = 0.e0
      CASE ( 'F' )
         pt2d(:,nlcj-jprecj+1:jpj+jpr2dj) = 0.e0
      END SELECT


      ! 2. East and west directions
      ! ---------------------------

      ! 2.1 Read Dirichlet lateral conditions

      SELECT CASE ( nbondi )
      CASE ( -1, 0, 1 )    ! all except 2
         iihom = nlci-nreci-jpr2di
         DO jl = 1, ipreci
            tr2ew(:,jl,1) = pt2d(jpreci+jl,:)
            tr2we(:,jl,1) = pt2d(iihom +jl,:)
         END DO
      END SELECT

      ! 2.2 Migrations

#if defined key_mpp_shmem
      !! * SHMEM version

      imigr = ipreci * ( jpj + 2*jpr2dj)

      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL shmem_put( tr2we(1-jpr2dj,1,2), tr2we(1,1,1), imigr, noea )
      CASE ( 0 )
         CALL shmem_put( tr2ew(1-jpr2dj,1,2), tr2ew(1,1,1), imigr, nowe )
         CALL shmem_put( tr2we(1-jpr2dj,1,2), tr2we(1,1,1), imigr, noea )
      CASE ( 1 )
         CALL shmem_put( tr2ew(1-jpr2dj,1,2), tr2ew(1,1,1), imigr, nowe )
      END SELECT

      CALL barrier()
      CALL shmem_udcflush()

#elif defined key_mpp_mpi
      !! * MPI version

      imigr = ipreci * ( jpj + 2*jpr2dj)

      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req1 )
         CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
         CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req2 )
         CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
         CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT

#endif

      ! 2.3 Write Dirichlet lateral conditions

      iihom = nlci - jpreci

      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, ipreci
            pt2d(iihom+jl,:) = tr2ew(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, ipreci
            pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
            pt2d( iihom+jl,:) = tr2ew(:,jl,2)
         END DO
      CASE ( 1 )
         DO jl = 1, ipreci
            pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
         END DO
      END SELECT


      ! 3. North and south directions
      ! -----------------------------

      ! 3.1 Read Dirichlet lateral conditions

      IF( nbondj /= 2 ) THEN
         ijhom = nlcj-nrecj-jpr2dj
         DO jl = 1, iprecj
            tr2sn(:,jl,1) = pt2d(:,ijhom +jl)
            tr2ns(:,jl,1) = pt2d(:,jprecj+jl)
         END DO
      ENDIF

      ! 3.2 Migrations

#if defined key_mpp_shmem
      !! * SHMEM version

      imigr = iprecj * ( jpi + 2*jpr2di )

      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL shmem_put( tr2sn(1-jpr2di,1,2), tr2sn(1,1,1), imigr, nono )
      CASE ( 0 )
         CALL shmem_put( tr2ns(1-jpr2di,1,2), tr2ns(1,1,1), imigr, noso )
         CALL shmem_put( tr2sn(1-jpr2di,1,2), tr2sn(1,1,1), imigr, nono )
      CASE ( 1 )
         CALL shmem_put( tr2ns(1-jpr2di,1,2), tr2ns(1,1,1), imigr, noso )
      END SELECT 
      CALL barrier()
      CALL shmem_udcflush()

#elif defined key_mpp_mpi
      !! * MPI version

      imigr = iprecj * ( jpi + 2*jpr2di )

      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, tr2sn(1-jpr2di,1,1), imigr, nono, ml_req1 )
         CALL mpprecv( 3, tr2ns(1-jpr2di,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 3, tr2ns(1-jpr2di,1,1), imigr, noso, ml_req1 )
         CALL mppsend( 4, tr2sn(1-jpr2di,1,1), imigr, nono, ml_req2 )
         CALL mpprecv( 3, tr2ns(1-jpr2di,1,2), imigr )
         CALL mpprecv( 4, tr2sn(1-jpr2di,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 3, tr2ns(1-jpr2di,1,1), imigr, noso, ml_req1 )
         CALL mpprecv( 4, tr2sn(1-jpr2di,1,2), imigr )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
  
#endif

      ! 3.3 Write Dirichlet lateral conditions

      ijhom = nlcj - jprecj  

      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, iprecj
            pt2d(:,ijhom+jl) = tr2ns(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, iprecj
            pt2d(:,jl-jpr2dj) = tr2sn(:,jl,2)
            pt2d(:,ijhom+jl ) = tr2ns(:,jl,2)
         END DO
      CASE ( 1 ) 
         DO jl = 1, iprecj
            pt2d(:,jl-jpr2dj) = tr2sn(:,jl,2)
         END DO
      END SELECT 
  

      ! 4. north fold treatment
      ! -----------------------
  
      ! 4.1 treatment without exchange (jpni odd)
      
      SELECT CASE ( jpni )
  
      CASE ( 1 ) ! only one proc along I, no mpp exchange
  
         SELECT CASE ( npolj )
  
         CASE ( 3 , 4 )   !  T pivot
            iloc = jpiglo - 2 * ( nimpp - 1 )
  
            SELECT CASE ( cd_type )
  
            CASE ( 'T', 'S', 'W' )
               DO jl = 0, iprecj-1
                  DO ji = 2-jpr2di, nlci+jpr2di
                     ijt=iloc-ji+2
                     pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-2-jl)
                  END DO
               END DO
               DO ji = nlci/2+1, nlci+jpr2di
                  ijt=iloc-ji+2
                  pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
               END DO
 
            CASE ( 'U' )
               DO jl =0, iprecj-1
                  DO ji = 1-jpr2di, nlci-1-jpr2di
                     iju=iloc-ji+1
                     pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-2-jl)
                  END DO
               END DO
               DO ji = nlci/2, nlci-1+jpr2di
                  iju=iloc-ji+1
                  pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
               END DO
  
            CASE ( 'V' )
               DO jl = -1, iprecj-1
                  DO ji = 2-jpr2di, nlci+jpr2di
                     ijt=iloc-ji+2
                     pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-3-jl)
                  END DO
               END DO
  
            CASE ( 'F', 'G' )
               DO jl = -1, iprecj-1
                  DO ji = 1-jpr2di, nlci-1+jpr2di
                     iju=iloc-ji+1
                     pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-3-jl)
                  END DO
               END DO
  
            CASE ( 'I' )                                  ! ice U-V point
               DO jl = 0, iprecj-1
                  pt2d(2,nlcj+jl) = psgn * pt2d(3,nlcj-1-jl)
                  DO ji = 3, nlci+jpr2di
                     iju = iloc - ji + 3
                     pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-1-jl)
                  END DO
               END DO
  
            END SELECT
  
         CASE ( 5 , 6 )                 ! F pivot
            iloc=jpiglo-2*(nimpp-1)
  
            SELECT CASE (cd_type )
  
            CASE ( 'T', 'S', 'W' )
               DO jl = 0, iprecj-1
                  DO ji = 1-jpr2di, nlci+jpr2di
                     ijt=iloc-ji+1
                     pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-1-jl)
                  END DO
               END DO
  
            CASE ( 'U' )
               DO jl = 0, iprecj-1
                  DO ji = 1-jpr2di, nlci-1+jpr2di
                     iju=iloc-ji
                     pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-1-jl)
                  END DO
               END DO
 
            CASE ( 'V' )
               DO jl = 0, iprecj-1
                  DO ji = 1-jpr2di, nlci+jpr2di
                     ijt=iloc-ji+1
                     pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-2-jl)
                  END DO
               END DO 
               DO ji = nlci/2+1, nlci+jpr2di
                  ijt=iloc-ji+1
                  pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
               END DO
  
            CASE ( 'F', 'G' )
               DO jl = 0, iprecj-1
                  DO ji = 1-jpr2di, nlci-1+jpr2di
                     iju=iloc-ji
                     pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-2-jl)
                  END DO
               END DO
               DO ji = nlci/2+1, nlci-1+jpr2di
                  iju=iloc-ji
                  pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
               END DO
  
            CASE ( 'I' )                                  ! ice U-V point
               pt2d( 2 ,nlcj) = 0.e0
               DO jl = 0, iprecj-1
                  DO ji = 2 , nlci-1+jpr2di
                     ijt = iloc - ji + 2
                     pt2d(ji,nlcj+jl)= 0.5 * ( pt2d(ji,nlcj-1-jl) + psgn * pt2d(ijt,nlcj-1-jl) )
                  END DO
               END DO
  
            END SELECT   ! cd_type
  
         END SELECT   ! npolj

      CASE DEFAULT   ! more than 1 proc along I
         IF( npolj /= 0 )   CALL mpp_lbc_north_e( pt2d, cd_type, psgn )   ! only for northern procs
         
      END SELECT   ! jpni


      ! 5. East and west directions
      ! ---------------------------

      SELECT CASE ( npolj )

      CASE ( 3, 4, 5, 6 )

         ! 5.1 Read Dirichlet lateral conditions

         SELECT CASE ( nbondi )
         CASE ( -1, 0, 1 )
            iihom = nlci-nreci-jpr2di
            DO jl = 1, ipreci
               tr2ew(:,jl,1) = pt2d(jpreci+jl,:)
               tr2we(:,jl,1) = pt2d(iihom +jl,:)
            END DO
         END SELECT

         ! 5.2 Migrations

#if defined key_mpp_shmem
         !! * SHMEM version

         imigr = ipreci * ( jpj + 2*jpr2dj )

         SELECT CASE ( nbondi )
         CASE ( -1 )
            CALL shmem_put( tr2we(1-jpr2dj,1,2), tr2we(1,1,1), imigr, noea )
         CASE ( 0 )
            CALL shmem_put( tr2ew(1-jpr2dj,1,2), tr2ew(1,1,1), imigr, nowe )
            CALL shmem_put( tr2we(1-jpr2dj,1,2), tr2we(1,1,1), imigr, noea )
         CASE ( 1 )
            CALL shmem_put( tr2ew(1-jpr2dj,1,2), tr2ew(1,1,1), imigr, nowe )
         END SELECT

         CALL barrier()
         CALL shmem_udcflush()
  
#elif defined key_mpp_mpi
         !! * MPI version
  
         imigr=ipreci* ( jpj + 2*jpr2dj )
  
         SELECT CASE ( nbondi )
         CASE ( -1 )
            CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req1 )
            CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         CASE ( 0 )
            CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
            CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req2 )
            CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
            CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
            IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
         CASE ( 1 )
            CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
            CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
            IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         END SELECT 
#endif

         ! 5.3 Write Dirichlet lateral conditions
  
         iihom = nlci - jpreci
  
         SELECT CASE ( nbondi )
         CASE ( -1 )
            DO jl = 1, ipreci
               pt2d(iihom+jl,:) = tr2ew(:,jl,2)
            END DO
         CASE ( 0 )
            DO jl = 1, ipreci
               pt2d(jl- jpr2di,:) = tr2we(:,jl,2)
               pt2d(iihom+jl,:) = tr2ew(:,jl,2)
            END DO
         CASE ( 1 )
            DO jl = 1, ipreci
               pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
            END DO
         END SELECT 
  
      END SELECT   ! npolj
  
   END SUBROUTINE mpp_lnk_2d_e


   SUBROUTINE mpplnks( ptab )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpplnks  ***
      !!
      !! ** Purpose :   Message passing manadgement for add 2d array local boundary
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask between
      !!       processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors 
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   &
         ptab                     ! 2D array
  
      !! * Local variables
      INTEGER ::   ji, jl         ! dummy loop indices
      INTEGER ::   &
         imigr, iihom, ijhom      ! temporary integers
      INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
      INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
      !!----------------------------------------------------------------------


      ! 1. north fold treatment
      ! -----------------------

      ! 1.1 treatment without exchange (jpni odd)
  
      SELECT CASE ( npolj )
      CASE ( 4 )
         DO ji = 1, nlci
            ptab(ji,nlcj-2) = ptab(ji,nlcj-2) + t2p1(ji,1,1)
         END DO
      CASE ( 6 )
         DO ji = 1, nlci
            ptab(ji,nlcj-1) = ptab(ji,nlcj-1) + t2p1(ji,1,1)
         END DO

      ! 1.2 treatment with exchange (jpni greater than 1)
      ! 
      CASE ( 3 )
#if defined key_mpp_shmem
  
         !! * SHMEN version
  
         imigr=jprecj*jpi
  
         CALL shmem_put(t2p1(1,1,2),t2p1(1,1,1),imigr,nono)
         CALL barrier()
         CALL shmem_udcflush()

#  elif defined key_mpp_mpi
       !! * MPI version

       imigr=jprecj*jpi

       CALL mppsend(3,t2p1(1,1,1),imigr,nono, ml_req1)
       CALL mpprecv(3,t2p1(1,1,2),imigr)
       IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)

#endif      

       ! Write north fold conditions

       DO ji = 1, nlci
          ptab(ji,nlcj-2) = ptab(ji,nlcj-2)+t2p1(ji,1,2)
       END DO

    CASE ( 5 )

#if defined key_mpp_shmem

       !! * SHMEN version

       imigr=jprecj*jpi

       CALL shmem_put(t2p1(1,1,2),t2p1(1,1,1),imigr,nono)
       CALL barrier()
       CALL shmem_udcflush()

#  elif defined key_mpp_mpi
       !! * Local variables   (MPI version)

       imigr=jprecj*jpi

       CALL mppsend(3,t2p1(1,1,1),imigr,nono, ml_req1)
       CALL mpprecv(3,t2p1(1,1,2),imigr)
       IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)

#endif      

       ! Write north fold conditions

       DO ji = 1, nlci
          ptab(ji,nlcj-1) = ptab(ji,nlcj-1)+t2p1(ji,1,2)
       END DO

    END SELECT


    ! 2. East and west directions
    ! ---------------------------

    ! 2.1 Read Dirichlet lateral conditions

    iihom = nlci-jpreci

    SELECT CASE ( nbondi )

    CASE ( -1, 0, 1 )  ! all except 2
       DO jl = 1, jpreci
             t2ew(:,jl,1) = ptab(  jl    ,:)
             t2we(:,jl,1) = ptab(iihom+jl,:)
       END DO
    END SELECT

    ! 2.2 Migrations

#if defined key_mpp_shmem

    !! * SHMEN version

    imigr=jpreci*jpj

    SELECT CASE ( nbondi )

    CASE ( -1 )
       CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)

    CASE ( 0 )
       CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
       CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)

    CASE ( 1 )
       CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)

    END SELECT
    CALL  barrier()
    CALL  shmem_udcflush()

#  elif defined key_mpp_mpi
    !! * Local variables   (MPI version)

    imigr=jpreci*jpj

    SELECT CASE ( nbondi )

    CASE ( -1 )
       CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req1)
       CALL mpprecv(1,t2ew(1,1,2),imigr)
       IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
    CASE ( 0 )
       CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
       CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req2)
       CALL mpprecv(1,t2ew(1,1,2),imigr)
       CALL mpprecv(2,t2we(1,1,2),imigr)
       IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
       IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)

    CASE ( 1 )
       CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
       CALL mpprecv(2,t2we(1,1,2),imigr)
       IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)

    END SELECT

#endif

    ! 2.3 Write Dirichlet lateral conditions

       iihom = nlci-nreci

    SELECT CASE ( nbondi )

    CASE ( -1 )
       DO jl = 1, jpreci
             ptab(iihom +jl,:) = ptab(iihom +jl,:)+t2ew(:,jl,2)
       END DO

    CASE ( 0 )
       DO jl = 1, jpreci
             ptab(jpreci+jl,:) = ptab(jpreci+jl,:)+t2we(:,jl,2)
             ptab(iihom +jl,:) = ptab(iihom +jl,:)+t2ew(:,jl,2)
       END DO

    CASE ( 1 )
       DO jl = 1, jpreci
             ptab(jpreci+jl,:) = ptab(jpreci+jl,:)+t2we(:,jl,2)
       END DO
    END SELECT


    ! 3. North and south directions
    ! -----------------------------

    ! 3.1 Read Dirichlet lateral conditions

    ijhom = nlcj-jprecj

    SELECT CASE ( nbondj )

    CASE ( -1, 0, 1 )
       DO jl = 1, jprecj
             t2sn(:,jl,1) = ptab(:,ijhom+jl)
             t2ns(:,jl,1) = ptab(:,   jl   )
       END DO

    END SELECT 

    ! 3.2 Migrations

#if defined key_mpp_shmem

    !! * SHMEN version

    imigr=jprecj*jpi

    SELECT CASE ( nbondj )

    CASE ( -1 )
       CALL shmem_put(t2sn(1,1,2),t2sn(1,1,1),imigr,nono)

    CASE ( 0 )
       CALL shmem_put(t2ns(1,1,2),t2ns(1,1,1),imigr,noso)
       CALL shmem_put(t2sn(1,1,2),t2sn(1,1,1),imigr,nono)

    CASE ( 1 )
       CALL shmem_put(t2ns(1,1,2),t2ns(1,1,1),imigr,noso)

    END SELECT
    CALL  barrier()
    CALL  shmem_udcflush()

#  elif defined key_mpp_mpi
    !! * Local variables   (MPI version)

    imigr=jprecj*jpi

    SELECT CASE ( nbondj )

    CASE ( -1 )
       CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req1)
       CALL mpprecv(3,t2ns(1,1,2),imigr)
       IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)

    CASE ( 0 )
       CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
       CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req2)
       CALL mpprecv(3,t2ns(1,1,2),imigr)
       CALL mpprecv(4,t2sn(1,1,2),imigr)
       IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
       IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)

    CASE ( 1 )
       CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
       CALL mpprecv(4,t2sn(1,1,2),imigr)
       IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
    END SELECT

#endif

    ! 3.3 Write Dirichlet lateral conditions

       ijhom = nlcj-nrecj

    SELECT CASE ( nbondj )

    CASE ( -1 )
       DO jl = 1, jprecj
             ptab(:,ijhom +jl) = ptab(:,ijhom +jl)+t2ns(:,jl,2)
       END DO

    CASE ( 0 )
       DO jl = 1, jprecj
             ptab(:,jprecj+jl) = ptab(:,jprecj+jl)+t2sn(:,jl,2)
             ptab(:,ijhom +jl) = ptab(:,ijhom +jl)+t2ns(:,jl,2)
       END DO

    CASE ( 1 ) 
       DO jl = 1, jprecj
             ptab(:,jprecj+jl) = ptab(:,jprecj+jl)+t2sn(:,jl,2)
       END DO

    END SELECT

  END SUBROUTINE mpplnks


   SUBROUTINE mppsend( ktyp, pmess, kbytes, kdest, md_req)
      !!----------------------------------------------------------------------
      !!                  ***  routine mppsend  ***
      !!                   
      !! ** Purpose :   Send messag passing array
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), INTENT(inout) ::   pmess(*)       ! array of real
      INTEGER , INTENT( in  ) ::   kbytes,     &  ! size of the array pmess
         &                         kdest ,     &  ! receive process number
         &                         ktyp,       &  ! Tag of the message
         &                         md_req         ! Argument for isend
      !!----------------------------------------------------------------------
#if defined key_mpp_shmem
      !! * SHMEM version  :    routine not used

#elif defined key_mpp_mpi
      !! * MPI version
      INTEGER ::   iflag

      SELECT CASE ( c_mpi_send )
      CASE ( 'S' )                ! Standard mpi send (blocking)
         CALL mpi_send ( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
            &                          mpi_comm_opa, iflag )
      CASE ( 'B' )                ! Buffer mpi send (blocking)
         CALL mpi_bsend( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
            &                          mpi_comm_opa, iflag )
      CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
         ! Be carefull, one more argument here : the mpi request identifier..
         CALL mpi_isend( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
            &                          mpi_comm_opa, md_req, iflag )
      END SELECT
#endif

   END SUBROUTINE mppsend


   SUBROUTINE mpprecv( ktyp, pmess, kbytes )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpprecv  ***
      !!
      !! ** Purpose :   Receive messag passing array
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), INTENT(inout) ::   pmess(*)       ! array of real
      INTEGER , INTENT( in  ) ::   kbytes,     &  ! suze of the array pmess
         &                         ktyp           ! Tag of the recevied message
      !!----------------------------------------------------------------------
#if defined key_mpp_shmem
      !! * SHMEM version  :    routine not used

#  elif defined key_mpp_mpi
      !! * MPI version
      INTEGER :: istatus(mpi_status_size)
      INTEGER :: iflag

      CALL mpi_recv( pmess, kbytes, mpi_double_precision, mpi_any_source, ktyp,   &
         &                          mpi_comm_opa, istatus, iflag )
#endif

   END SUBROUTINE mpprecv


   SUBROUTINE mppgather( ptab, kp, pio )
      !!----------------------------------------------------------------------
      !!                   ***  routine mppgather  ***
      !!                   
      !! ** Purpose :   Transfert between a local subdomain array and a work 
      !!     array which is distributed following the vertical level.
      !!
      !! ** Method  :
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj),       INTENT( in  ) ::   ptab   ! subdomain input array
      INTEGER ,                           INTENT( in  ) ::   kp     ! record length
      REAL(wp), DIMENSION(jpi,jpj,jpnij), INTENT( out ) ::   pio    ! subdomain input array
      !!---------------------------------------------------------------------
#if defined key_mpp_shmem
      !! * SHMEM version

      CALL barrier()
      CALL shmem_put( pio(1,1,npvm_me+1), ptab, jpi*jpj, kp )
      CALL barrier()

#elif defined key_mpp_mpi
      !! * Local variables   (MPI version)
      INTEGER :: itaille,ierror
  
      itaille=jpi*jpj
      CALL mpi_gather( ptab, itaille, mpi_double_precision, pio, itaille,   &
         &                            mpi_double_precision, kp , mpi_comm_opa, ierror ) 
#endif

   END SUBROUTINE mppgather


   SUBROUTINE mppscatter( pio, kp, ptab )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppscatter  ***
      !!
      !! ** Purpose :   Transfert between awork array which is distributed 
      !!      following the vertical level and the local subdomain array.
      !!
      !! ** Method :
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpnij)  ::  pio        ! output array
      INTEGER                             ::   kp        ! Tag (not used with MPI
      REAL(wp), DIMENSION(jpi,jpj)        ::  ptab       ! subdomain array input
      !!---------------------------------------------------------------------
#if defined key_mpp_shmem
      !! * SHMEM version

      CALL barrier()
      CALL shmem_get( ptab, pio(1,1,npvm_me+1), jpi*jpj, kp )
      CALL barrier()

#  elif defined key_mpp_mpi
      !! * Local variables   (MPI version)
      INTEGER :: itaille, ierror
  
      itaille=jpi*jpj

      CALL mpi_scatter( pio, itaille, mpi_double_precision, ptab, itaille,   &
         &                            mpi_double_precision, kp, mpi_comm_opa, ierror )
#endif

   END SUBROUTINE mppscatter


   SUBROUTINE mppisl_a_int( ktab, kdim )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppisl_a_int  ***
      !!                   
      !! ** Purpose :   Massively parallel processors
      !!                Find the  non zero value
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in  )                  ::   kdim       ! ???
      INTEGER, INTENT(inout), DIMENSION(kdim) ::   ktab       ! ???
  
#if defined key_mpp_shmem
      !! * Local variables   (SHMEM version)
      INTEGER :: ji
      INTEGER, SAVE :: ibool=0

      IF( kdim > jpmppsum ) CALL ctl_stop( 'mppisl_a_int routine : kdim is too big', &
           &                               'change jpmppsum dimension in mpp.h' )

      DO ji = 1, kdim
         niitab_shmem(ji) = ktab(ji)
      END DO
      CALL  barrier()
      IF(ibool == 0 ) THEN 
         CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,kdim,0   &
              ,0,N$PES,ni11wrk_shmem,ni11sync_shmem)
         CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,kdim,0   &
              ,0,N$PES,ni12wrk_shmem,ni12sync_shmem)
      ELSE
         CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,kdim,0   &
              ,0,N$PES,ni21wrk_shmem,ni21sync_shmem)
         CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,kdim,0   &
              ,0,N$PES,ni22wrk_shmem,ni22sync_shmem)
      ENDIF
      CALL  barrier()
      ibool=ibool+1
      ibool=MOD( ibool,2)
      DO ji = 1, kdim
         IF( ni11tab_shmem(ji) /= 0. ) THEN
            ktab(ji) = ni11tab_shmem(ji)
         ELSE
            ktab(ji) = ni12tab_shmem(ji)
         ENDIF
      END DO
  
#  elif defined key_mpp_mpi
      !! * Local variables   (MPI version)
      LOGICAL  :: lcommute
      INTEGER, DIMENSION(kdim) ::   iwork
      INTEGER  :: mpi_isl,ierror
  
      lcommute = .TRUE.
      CALL mpi_op_create( lc_isl, lcommute, mpi_isl, ierror )
      CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer   &
           , mpi_isl, mpi_comm_opa, ierror )
      ktab(:) = iwork(:)
#endif

   END SUBROUTINE mppisl_a_int


   SUBROUTINE mppisl_int( ktab )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppisl_int  ***
      !!                   
      !! ** Purpose :   Massively parallel processors
      !!                Find the non zero value
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER , INTENT( inout ) ::   ktab        ! 

#if defined key_mpp_shmem
      !! * Local variables   (SHMEM version)
      INTEGER, SAVE :: ibool=0

      niitab_shmem(1) = ktab
      CALL  barrier()
      IF(ibool == 0 ) THEN 
         CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,1,0   &
              ,0,N$PES,ni11wrk_shmem,ni11sync_shmem)
         CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,1,0   &
              ,0,N$PES,ni12wrk_shmem,ni12sync_shmem)
      ELSE
         CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,1,0   &
              ,0,N$PES,ni21wrk_shmem,ni21sync_shmem)
         CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,1,0   &
              ,0,N$PES,ni22wrk_shmem,ni22sync_shmem)
      ENDIF
      CALL  barrier()
      ibool=ibool+1
      ibool=MOD( ibool,2)
      IF( ni11tab_shmem(1) /= 0. ) THEN
         ktab = ni11tab_shmem(1)
      ELSE
         ktab = ni12tab_shmem(1)
      ENDIF
  
#  elif defined key_mpp_mpi
  
      !! * Local variables   (MPI version)
      LOGICAL :: lcommute
      INTEGER :: mpi_isl,ierror
      INTEGER ::   iwork
  
      lcommute = .TRUE.
      CALL mpi_op_create(lc_isl,lcommute,mpi_isl,ierror)
      CALL mpi_allreduce(ktab, iwork, 1,mpi_integer   &
           ,mpi_isl,mpi_comm_opa,ierror)
      ktab = iwork
#endif

   END SUBROUTINE mppisl_int


   SUBROUTINE mppmax_a_int( ktab, kdim )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppmax_a_int  ***
      !! 
      !! ** Purpose :   Find maximum value in an integer layout array
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER , INTENT( in  )                  ::   kdim        ! size of array
      INTEGER , INTENT(inout), DIMENSION(kdim) ::   ktab        ! input array
  
#if defined key_mpp_shmem
      !! * Local declarations    (SHMEM version)
      INTEGER :: ji
      INTEGER, SAVE :: ibool=0
  
      IF( kdim > jpmppsum ) CALL ctl_stop( 'mppmax_a_int routine : kdim is too big', &
           &                               'change jpmppsum dimension in mpp.h' )
  
      DO ji = 1, kdim
         niltab_shmem(ji) = ktab(ji)
      END DO
      CALL  barrier()
      IF(ibool == 0 ) THEN 
         CALL shmem_int8_max_to_all (niltab_shmem,niltab_shmem,kdim,0,0   &
              ,N$PES,nil1wrk_shmem,nil1sync_shmem )
      ELSE
         CALL shmem_int8_max_to_all (niltab_shmem,niltab_shmem,kdim,0,0   &
              ,N$PES,nil2wrk_shmem,nil2sync_shmem )
      ENDIF
      CALL  barrier()
      ibool=ibool+1
      ibool=MOD( ibool,2)
      DO ji = 1, kdim
         ktab(ji) = niltab_shmem(ji)
      END DO
  
#  elif defined key_mpp_mpi
  
      !! * Local variables   (MPI version)
      INTEGER :: ierror
      INTEGER, DIMENSION(kdim) ::   iwork
  
      CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer,   &
           &                mpi_max, mpi_comm_opa, ierror )
  
      ktab(:) = iwork(:)
#endif

   END SUBROUTINE mppmax_a_int


   SUBROUTINE mppmax_int( ktab )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppmax_int  ***
      !!
      !! ** Purpose :
      !!     Massively parallel processors
      !!     Find maximum value in an integer layout array
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(inout) ::   ktab      ! ???
  
      !! * Local declarations

#if defined key_mpp_shmem

      !! * Local variables   (SHMEM version)
      INTEGER :: ji
      INTEGER, SAVE :: ibool=0
  
      niltab_shmem(1) = ktab
      CALL  barrier()
      IF(ibool == 0 ) THEN 
         CALL shmem_int8_max_to_all (niltab_shmem,niltab_shmem, 1,0,0   &
              ,N$PES,nil1wrk_shmem,nil1sync_shmem )
      ELSE
         CALL shmem_int8_max_to_all (niltab_shmem,niltab_shmem, 1,0,0   &
              ,N$PES,nil2wrk_shmem,nil2sync_shmem )
      ENDIF
      CALL  barrier()
      ibool=ibool+1
      ibool=MOD( ibool,2)
      ktab = niltab_shmem(1)
  
#  elif defined key_mpp_mpi

      !! * Local variables   (MPI version)
      INTEGER ::  ierror, iwork
  
      CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
           &              ,mpi_max,mpi_comm_opa,ierror)
  
      ktab = iwork
#endif

   END SUBROUTINE mppmax_int


   SUBROUTINE mppmin_a_int( ktab, kdim )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppmin_a_int  ***
      !! 
      !! ** Purpose :   Find minimum value in an integer layout array
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER , INTENT( in  )                  ::   kdim        ! size of array
      INTEGER , INTENT(inout), DIMENSION(kdim) ::   ktab        ! input array
  
#if defined key_mpp_shmem
      !! * Local declarations    (SHMEM version)
      INTEGER :: ji
      INTEGER, SAVE :: ibool=0
  
      IF( kdim > jpmppsum ) CALL ctl_stop( 'mppmin_a_int routine : kdim is too big', &
           &                               'change jpmppsum dimension in mpp.h' )
  
      DO ji = 1, kdim
         niltab_shmem(ji) = ktab(ji)
      END DO
      CALL  barrier()
      IF(ibool == 0 ) THEN 
         CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem,kdim,0,0   &
              ,N$PES,nil1wrk_shmem,nil1sync_shmem )
      ELSE
         CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem,kdim,0,0   &
              ,N$PES,nil2wrk_shmem,nil2sync_shmem )
      ENDIF
      CALL  barrier()
      ibool=ibool+1
      ibool=MOD( ibool,2)
      DO ji = 1, kdim
         ktab(ji) = niltab_shmem(ji)
      END DO
  
#  elif defined key_mpp_mpi
  
      !! * Local variables   (MPI version)
      INTEGER :: ierror
      INTEGER, DIMENSION(kdim) ::   iwork
  
      CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer,   &
           &                mpi_min, mpi_comm_opa, ierror )
  
      ktab(:) = iwork(:)
#endif

   END SUBROUTINE mppmin_a_int


   SUBROUTINE mppmin_int( ktab )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppmin_int  ***
      !!
      !! ** Purpose :
      !!     Massively parallel processors
      !!     Find minimum value in an integer layout array
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(inout) ::   ktab      ! ???
  
      !! * Local declarations

#if defined key_mpp_shmem

      !! * Local variables   (SHMEM version)
      INTEGER :: ji
      INTEGER, SAVE :: ibool=0
  
      niltab_shmem(1) = ktab
      CALL  barrier()
      IF(ibool == 0 ) THEN 
         CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem, 1,0,0   &
              ,N$PES,nil1wrk_shmem,nil1sync_shmem )
      ELSE
         CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem, 1,0,0   &
              ,N$PES,nil2wrk_shmem,nil2sync_shmem )
      ENDIF
      CALL  barrier()
      ibool=ibool+1
      ibool=MOD( ibool,2)
      ktab = niltab_shmem(1)
  
#  elif defined key_mpp_mpi

      !! * Local variables   (MPI version)
      INTEGER ::  ierror, iwork
  
      CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
           &              ,mpi_min,mpi_comm_opa,ierror)
  
      ktab = iwork
#endif

   END SUBROUTINE mppmin_int


   SUBROUTINE mppsum_a_int( ktab, kdim )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppsum_a_int  ***
      !!                    
      !! ** Purpose :   Massively parallel processors
      !!                Global integer sum
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in  )                   ::   kdim      ! ???
      INTEGER, INTENT(inout), DIMENSION (kdim) ::   ktab      ! ???
  
#if defined key_mpp_shmem

      !! * Local variables   (SHMEM version)
      INTEGER :: ji
      INTEGER, SAVE :: ibool=0

      IF( kdim > jpmppsum ) CALL ctl_stop( 'mppsum_a_int routine : kdim is too big', &
           &                               'change jpmppsum dimension in mpp.h' )

      DO ji = 1, kdim
         nistab_shmem(ji) = ktab(ji)
      END DO
      CALL  barrier()
      IF(ibool == 0 ) THEN 
         CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem,kdim,0,0,   &
              N$PES,nis1wrk_shmem,nis1sync_shmem)
      ELSE
         CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem,kdim,0,0,   &
              N$PES,nis2wrk_shmem,nis2sync_shmem)
      ENDIF
      CALL  barrier()
      ibool = ibool + 1
      ibool = MOD( ibool, 2 )
      DO ji = 1, kdim
         ktab(ji) = nistab_shmem(ji)
      END DO
  
#  elif defined key_mpp_mpi

      !! * Local variables   (MPI version)
      INTEGER :: ierror
      INTEGER, DIMENSION (kdim) ::  iwork
  
      CALL mpi_allreduce(ktab, iwork,kdim,mpi_integer   &
           ,mpi_sum,mpi_comm_opa,ierror)
  
      ktab(:) = iwork(:)
#endif

   END SUBROUTINE mppsum_a_int


  SUBROUTINE mppsum_int( ktab )
    !!----------------------------------------------------------------------
    !!                 ***  routine mppsum_int  ***
    !!                  
    !! ** Purpose :   Global integer sum
    !!
    !!----------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT(inout) ::   ktab

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    INTEGER, SAVE :: ibool=0

    nistab_shmem(1) = ktab
    CALL  barrier()
    IF(ibool == 0 ) THEN 
       CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem, 1,0,0,   &
            N$PES,nis1wrk_shmem,nis1sync_shmem)
    ELSE
       CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem, 1,0,0,   &
            N$PES,nis2wrk_shmem,nis2sync_shmem)
    ENDIF
    CALL  barrier()
    ibool=ibool+1
    ibool=MOD( ibool,2)
    ktab = nistab_shmem(1)

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER :: ierror, iwork

    CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
         ,mpi_sum,mpi_comm_opa,ierror)

    ktab = iwork

#endif

  END SUBROUTINE mppsum_int


  SUBROUTINE mppisl_a_real( ptab, kdim )
    !!----------------------------------------------------------------------
    !!                 ***  routine mppisl_a_real  ***
    !!         
    !! ** Purpose :   Massively parallel processors
    !!           Find the non zero island barotropic stream function value
    !!
    !!   Modifications:
    !!        !  93-09 (M. Imbard)
    !!        !  96-05 (j. Escobar)
    !!        !  98-05 (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI 
    !!----------------------------------------------------------------------
    INTEGER , INTENT( in  )                  ::   kdim      ! ???
    REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab      ! ???

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    INTEGER :: ji
    INTEGER, SAVE :: ibool=0

    IF( kdim > jpmppsum ) CALL ctl_stop( 'mppisl_a_real routine : kdim is too big', &
         &                               'change jpmppsum dimension in mpp.h' )

    DO ji = 1, kdim
       wiltab_shmem(ji) = ptab(ji)
    END DO
    CALL  barrier()
    IF(ibool == 0 ) THEN 
       CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem,kdim,0   &
            ,0,N$PES,wi11wrk_shmem,ni11sync_shmem)
       CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem,kdim,0   &
            ,0,N$PES,wi12wrk_shmem,ni12sync_shmem)
    ELSE
       CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem,kdim,0   &
            ,0,N$PES,wi21wrk_shmem,ni21sync_shmem)
       CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem,kdim,0   &
            ,0,N$PES,wi22wrk_shmem,ni22sync_shmem)
    ENDIF
    CALL  barrier()
    ibool=ibool+1
    ibool=MOD( ibool,2)
    DO ji = 1, kdim
       IF(wi1tab_shmem(ji) /= 0. ) THEN
          ptab(ji) = wi1tab_shmem(ji)
       ELSE
          ptab(ji) = wi2tab_shmem(ji)
       ENDIF
    END DO

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    LOGICAL ::   lcommute = .TRUE.
    INTEGER ::   mpi_isl, ierror
    REAL(wp), DIMENSION(kdim) ::  zwork

    CALL mpi_op_create(lc_isl,lcommute,mpi_isl,ierror)
    CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
         ,mpi_isl,mpi_comm_opa,ierror)
    ptab(:) = zwork(:)

#endif

  END SUBROUTINE mppisl_a_real


   SUBROUTINE mppisl_real( ptab )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppisl_real  ***
      !!                  
      !! ** Purpose :   Massively parallel processors
      !!       Find the  non zero island barotropic stream function value
      !!
      !!     Modifications:
      !!        !  93-09 (M. Imbard)
      !!        !  96-05 (j. Escobar)
      !!        !  98-05 (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI 
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout) ::   ptab

#if defined key_mpp_shmem

      !! * Local variables   (SHMEM version)
      INTEGER, SAVE :: ibool=0

      wiltab_shmem(1) = ptab
      CALL  barrier()
      IF(ibool == 0 ) THEN 
         CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem, 1,0   &
            ,0,N$PES,wi11wrk_shmem,ni11sync_shmem)
         CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem, 1,0   &
            ,0,N$PES,wi12wrk_shmem,ni12sync_shmem)
      ELSE
         CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem, 1,0   &
            ,0,N$PES,wi21wrk_shmem,ni21sync_shmem)
         CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem, 1,0   &
            ,0,N$PES,wi22wrk_shmem,ni22sync_shmem)
      ENDIF
      CALL barrier()
      ibool = ibool + 1
      ibool = MOD( ibool, 2 )
      IF( wi1tab_shmem(1) /= 0. ) THEN
         ptab = wi1tab_shmem(1)
      ELSE
         ptab = wi2tab_shmem(1)
      ENDIF

#  elif defined key_mpp_mpi

      !! * Local variables   (MPI version)
      LOGICAL  ::   lcommute = .TRUE.
      INTEGER  ::   mpi_isl, ierror
      REAL(wp) ::   zwork

      CALL mpi_op_create( lc_isl, lcommute, mpi_isl, ierror )
      CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision,   &
         &                                mpi_isl  , mpi_comm_opa, ierror )
      ptab = zwork

#endif

   END SUBROUTINE mppisl_real


  FUNCTION lc_isl( py, px, kdim, kdtatyp )
    INTEGER :: kdim
    REAL(wp), DIMENSION(kdim) ::  px, py
    INTEGER :: kdtatyp, ji
    INTEGER :: lc_isl
    DO ji = 1, kdim
       IF( py(ji) /= 0. )   px(ji) = py(ji)
    END DO
    lc_isl=0

  END FUNCTION lc_isl


  SUBROUTINE mppmax_a_real( ptab, kdim )
    !!----------------------------------------------------------------------
    !!                 ***  routine mppmax_a_real  ***
    !!                  
    !! ** Purpose :   Maximum
    !!
    !!----------------------------------------------------------------------
    !! * Arguments
    INTEGER , INTENT( in  )                  ::   kdim
    REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    INTEGER :: ji
    INTEGER, SAVE :: ibool=0

    IF( kdim > jpmppsum ) CALL ctl_stop( 'mppmax_a_real routine : kdim is too big', &
         &                               'change jpmppsum dimension in mpp.h' )

    DO ji = 1, kdim
       wintab_shmem(ji) = ptab(ji)
    END DO
    CALL  barrier()
    IF(ibool == 0 ) THEN 
       CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem,kdim,0   &
            ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
    ELSE
       CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem,kdim,0   &
            ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
    ENDIF
    CALL  barrier()
    ibool=ibool+1
    ibool=MOD( ibool,2)
    DO ji = 1, kdim
       ptab(ji) = wintab_shmem(ji)
    END DO

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER :: ierror
    REAL(wp), DIMENSION(kdim) ::  zwork

    CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
         ,mpi_max,mpi_comm_opa,ierror)
    ptab(:) = zwork(:)

#endif

  END SUBROUTINE mppmax_a_real


  SUBROUTINE mppmax_real( ptab )
    !!----------------------------------------------------------------------
    !!                  ***  routine mppmax_real  ***
    !!                    
    !! ** Purpose :   Maximum
    !!
    !!----------------------------------------------------------------------
    !! * Arguments
    REAL(wp), INTENT(inout) ::   ptab      ! ???

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    INTEGER, SAVE :: ibool=0

    wintab_shmem(1) = ptab
    CALL  barrier()
    IF(ibool == 0 ) THEN 
       CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem, 1,0   &
            ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
    ELSE
       CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem, 1,0   &
            ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
    ENDIF
    CALL  barrier()
    ibool=ibool+1
    ibool=MOD( ibool,2)
    ptab = wintab_shmem(1)

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER  ::   ierror
    REAL(wp) ::   zwork

    CALL mpi_allreduce( ptab, zwork  , 1             , mpi_double_precision,   &
       &                      mpi_max, mpi_comm_opa, ierror     )
    ptab = zwork

#endif

  END SUBROUTINE mppmax_real


  SUBROUTINE mppmin_a_real( ptab, kdim )
    !!----------------------------------------------------------------------
    !!                 ***  routine mppmin_a_real  ***
    !!                  
    !! ** Purpose :   Minimum
    !!
    !!-----------------------------------------------------------------------
    !! * Arguments
    INTEGER , INTENT( in  )                  ::   kdim
    REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    INTEGER :: ji
    INTEGER, SAVE :: ibool=0

    IF( kdim > jpmppsum ) CALL ctl_stop( 'mpprmin routine : kdim is too big', &
         &                               'change jpmppsum dimension in mpp.h' )

    DO ji = 1, kdim
       wintab_shmem(ji) = ptab(ji)
    END DO
    CALL  barrier()
    IF(ibool == 0 ) THEN 
       CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem,kdim,0   &
            ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
    ELSE
       CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem,kdim,0   &
            ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
    ENDIF
    CALL  barrier()
    ibool=ibool+1
    ibool=MOD( ibool,2)
    DO ji = 1, kdim
       ptab(ji) = wintab_shmem(ji)
    END DO

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER :: ierror
    REAL(wp), DIMENSION(kdim) ::   zwork

    CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
         ,mpi_min,mpi_comm_opa,ierror)
    ptab(:) = zwork(:)

#endif

  END SUBROUTINE mppmin_a_real


  SUBROUTINE mppmin_real( ptab )
    !!----------------------------------------------------------------------
    !!                  ***  routine mppmin_real  ***
    !! 
    !! ** Purpose :   minimum in Massively Parallel Processing
    !!                REAL scalar case
    !!
    !!-----------------------------------------------------------------------
    !! * Arguments
    REAL(wp), INTENT( inout ) ::   ptab        ! 

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    INTEGER, SAVE :: ibool=0

    wintab_shmem(1) = ptab
    CALL  barrier()
    IF(ibool == 0 ) THEN 
       CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem, 1,0   &
            ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
    ELSE
       CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem, 1,0   &
            ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
    ENDIF
    CALL  barrier()
    ibool=ibool+1
    ibool=MOD( ibool,2)
    ptab = wintab_shmem(1)

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER  ::   ierror
    REAL(wp) ::   zwork

    CALL mpi_allreduce( ptab, zwork, 1,mpi_double_precision   &
         &               ,mpi_min,mpi_comm_opa,ierror)
    ptab = zwork

#endif

  END SUBROUTINE mppmin_real


  SUBROUTINE mppsum_a_real( ptab, kdim )
    !!----------------------------------------------------------------------
    !!                  ***  routine mppsum_a_real  ***
    !! 
    !! ** Purpose :   global sum in Massively Parallel Processing
    !!                REAL ARRAY argument case
    !!
    !!-----------------------------------------------------------------------
    INTEGER , INTENT( in )                     ::   kdim      ! size of ptab
    REAL(wp), DIMENSION(kdim), INTENT( inout ) ::   ptab      ! input array

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    INTEGER :: ji
    INTEGER, SAVE :: ibool=0

    IF( kdim > jpmppsum ) CALL ctl_stop( 'mppsum_a_real routine : kdim is too big', &
         &                               'change jpmppsum dimension in mpp.h' )

    DO ji = 1, kdim
       wrstab_shmem(ji) = ptab(ji)
    END DO
    CALL  barrier()
    IF(ibool == 0 ) THEN 
       CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem,kdim,0   &
            ,0,N$PES,wrs1wrk_shmem,nrs1sync_shmem )
    ELSE
       CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem,kdim,0   &
            ,0,N$PES,wrs2wrk_shmem,nrs2sync_shmem )
    ENDIF
    CALL  barrier()
    ibool=ibool+1
    ibool=MOD( ibool,2)
    DO ji = 1, kdim
       ptab(ji) = wrstab_shmem(ji)
    END DO

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER                   ::   ierror    ! temporary integer
    REAL(wp), DIMENSION(kdim) ::   zwork     ! temporary workspace 

    CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
         &              ,mpi_sum,mpi_comm_opa,ierror)
    ptab(:) = zwork(:)

#endif

  END SUBROUTINE mppsum_a_real


  SUBROUTINE mppsum_real( ptab )
    !!----------------------------------------------------------------------
    !!                  ***  routine mppsum_real  ***
    !!              
    !! ** Purpose :   global sum in Massively Parallel Processing
    !!                SCALAR argument case
    !!
    !!-----------------------------------------------------------------------
    REAL(wp), INTENT(inout) ::   ptab        ! input scalar

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    INTEGER, SAVE :: ibool=0

    wrstab_shmem(1) = ptab
    CALL  barrier()
    IF(ibool == 0 ) THEN 
       CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem, 1,0   &
            ,0,N$PES,wrs1wrk_shmem,nrs1sync_shmem )
    ELSE
       CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem, 1,0   &
            ,0,N$PES,wrs2wrk_shmem,nrs2sync_shmem )
    ENDIF
    CALL  barrier()
    ibool = ibool + 1
    ibool = MOD( ibool, 2 )
    ptab = wrstab_shmem(1)

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER  ::   ierror
    REAL(wp) ::   zwork

    CALL mpi_allreduce(ptab, zwork, 1,mpi_double_precision   &
         &              ,mpi_sum,mpi_comm_opa,ierror)
    ptab = zwork

#endif

  END SUBROUTINE mppsum_real

  SUBROUTINE mppsum_realdd( ptab )
    !!----------------------------------------------------------------------
    !!                  ***  routine mppsum_realdd ***
    !!
    !! ** Purpose :   global sum in Massively Parallel Processing
    !!                SCALAR argument case for double-double precision
    !!
    !!-----------------------------------------------------------------------
    COMPLEX(wp), INTENT(inout)         :: ptab    ! input scalar

#if defined key_mpp_shmem

COMPILATION PROBLEM: MPPSUM_REALDD IMPLEMENTATION DOES NOT EXIST FOR SHMEM CASE

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER  ::    ierror
    COMPLEX(wp) :: zwork

! reduce local sums into global sum
    call MPI_ALLREDUCE (ptab, zwork, 1, MPI_DOUBLE_COMPLEX, &
                     MPI_SUMDD, MPI_COMM_WORLD, ierror)
    ptab = zwork
#endif

  END SUBROUTINE mppsum_realdd

  SUBROUTINE mppsum_a_realdd( ptab, kdim )
    !!----------------------------------------------------------------------
    !!                  ***  routine mppsum_a_realdd  ***
    !!
    !! ** Purpose :   global sum in Massively Parallel Processing
    !!                COMPLEX ARRAY case for double-double precision
    !!
    !!-----------------------------------------------------------------------
    INTEGER , INTENT( in )                        ::   kdim      ! size of ptab
    COMPLEX(wp), DIMENSION(kdim), INTENT( inout ) ::   ptab      ! input array

#if defined key_mpp_shmem

COMPILATION PROBLEM: MPPSUM_A_REALDD IMPLEMENTATION DOES NOT EXIST FOR SHMEM CASE

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER                      :: ierror    ! temporary integer
    COMPLEX(wp), DIMENSION(kdim) :: zwork     ! temporary workspace


    call MPI_ALLREDUCE (ptab, zwork, kdim, MPI_DOUBLE_COMPLEX, &
                     MPI_SUMDD, MPI_COMM_WORLD, ierror)
    ptab(:) = zwork(:)

#endif

  END SUBROUTINE mppsum_a_realdd



  SUBROUTINE mpp_minloc2d(ptab, pmask, pmin, ki,kj )
    !!------------------------------------------------------------------------
    !!             ***  routine mpp_minloc  ***
    !!
    !! ** Purpose :  Compute the global minimum of an array ptab
    !!              and also give its global position
    !!
    !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
    !!
    !! ** Arguments : I : ptab =local 2D array
    !!                O : pmin = global minimum
    !!                O : ki,kj = global position of minimum
    !!
    !! ** Author : J.M. Molines 10/10/2004
    !!--------------------------------------------------------------------------
#ifdef key_mpp_shmem
    CALL ctl_stop( ' mpp_minloc not yet available in SHMEM' )
# elif key_mpp_mpi
    !! * Arguments
    REAL(wp), DIMENSION (jpi,jpj), INTENT (in)  :: ptab ,& ! Local 2D array
         &                                         pmask   ! Local mask
    REAL(wp)                     , INTENT (out) :: pmin    ! Global minimum of ptab
    INTEGER                      , INTENT (out) :: ki,kj   ! index of minimum in global frame

    !! * Local variables
    REAL(wp) :: zmin   ! local minimum
    REAL(wp) ,DIMENSION(2,1) :: zain, zaout
    INTEGER, DIMENSION (2)  :: ilocs
    INTEGER :: ierror


    zmin  = MINVAL( ptab(:,:) , mask= pmask == 1.e0 )
    ilocs = MINLOC( ptab(:,:) , mask= pmask == 1.e0 )

    ki = ilocs(1) + nimpp - 1
    kj = ilocs(2) + njmpp - 1

    zain(1,:)=zmin
    zain(2,:)=ki+10000.*kj

    CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_OPA,ierror)

    pmin=zaout(1,1)
    kj= INT(zaout(2,1)/10000.)
    ki= INT(zaout(2,1) - 10000.*kj )
#endif

  END SUBROUTINE mpp_minloc2d


  SUBROUTINE mpp_minloc3d(ptab, pmask, pmin, ki,kj ,kk)
    !!------------------------------------------------------------------------
    !!             ***  routine mpp_minloc  ***
    !!
    !! ** Purpose :  Compute the global minimum of an array ptab
    !!              and also give its global position
    !!
    !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
    !!
    !! ** Arguments : I : ptab =local 2D array
    !!                O : pmin = global minimum
    !!                O : ki,kj = global position of minimum
    !!
    !! ** Author : J.M. Molines 10/10/2004
    !!--------------------------------------------------------------------------
#ifdef key_mpp_shmem
    CALL ctl_stop( ' mpp_minloc not yet available in SHMEM' )
# elif key_mpp_mpi
    !! * Arguments
    REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT (in)  :: ptab ,& ! Local 2D array
         &                                         pmask   ! Local mask
    REAL(wp)                     , INTENT (out) :: pmin    ! Global minimum of ptab
    INTEGER                      , INTENT (out) :: ki,kj,kk ! index of minimum in global frame

    !! * Local variables
    REAL(wp) :: zmin   ! local minimum
    REAL(wp) ,DIMENSION(2,1) :: zain, zaout
    INTEGER, DIMENSION (3)  :: ilocs
    INTEGER :: ierror


    zmin  = MINVAL( ptab(:,:,:) , mask= pmask == 1.e0 )
    ilocs = MINLOC( ptab(:,:,:) , mask= pmask == 1.e0 )

    ki = ilocs(1) + nimpp - 1
    kj = ilocs(2) + njmpp - 1
    kk = ilocs(3)

    zain(1,:)=zmin
    zain(2,:)=ki+10000.*kj+100000000.*kk

    CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_OPA,ierror)

    pmin=zaout(1,1)
    kk= INT(zaout(2,1)/100000000.)
    kj= INT(zaout(2,1) - kk * 100000000. )/10000
    ki= INT(zaout(2,1) - kk * 100000000. -kj * 10000. )
#endif

  END SUBROUTINE mpp_minloc3d


  SUBROUTINE mpp_maxloc2d(ptab, pmask, pmax, ki,kj )
    !!------------------------------------------------------------------------
    !!             ***  routine mpp_maxloc  ***
    !!
    !! ** Purpose :  Compute the global maximum of an array ptab
    !!              and also give its global position
    !!
    !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
    !!
    !! ** Arguments : I : ptab =local 2D array
    !!                O : pmax = global maximum
    !!                O : ki,kj = global position of maximum
    !!
    !! ** Author : J.M. Molines 10/10/2004
    !!--------------------------------------------------------------------------
#ifdef key_mpp_shmem
    CALL ctl_stop( ' mpp_maxloc not yet available in SHMEM' )
# elif key_mpp_mpi
    !! * Arguments
    REAL(wp), DIMENSION (jpi,jpj), INTENT (in)  :: ptab ,& ! Local 2D array
         &                                         pmask   ! Local mask
    REAL(wp)                     , INTENT (out) :: pmax    ! Global maximum of ptab
    INTEGER                      , INTENT (out) :: ki,kj   ! index of maximum in global frame

    !! * Local variables
    REAL(wp) :: zmax   ! local maximum
    REAL(wp) ,DIMENSION(2,1) :: zain, zaout
    INTEGER, DIMENSION (2)  :: ilocs
    INTEGER :: ierror


    zmax  = MAXVAL( ptab(:,:) , mask= pmask == 1.e0 )
    ilocs = MAXLOC( ptab(:,:) , mask= pmask == 1.e0 )

    ki = ilocs(1) + nimpp - 1
    kj = ilocs(2) + njmpp - 1

    zain(1,:)=zmax
    zain(2,:)=ki+10000.*kj

    CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_OPA,ierror)

    pmax=zaout(1,1)
    kj= INT(zaout(2,1)/10000.)
    ki= INT(zaout(2,1) - 10000.*kj )
#endif

  END SUBROUTINE mpp_maxloc2d

  SUBROUTINE mpp_maxloc3d(ptab, pmask, pmax, ki,kj,kk )
    !!------------------------------------------------------------------------
    !!             ***  routine mpp_maxloc  ***
    !!
    !! ** Purpose :  Compute the global maximum of an array ptab
    !!              and also give its global position
    !!
    !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
    !!
    !! ** Arguments : I : ptab =local 2D array
    !!                O : pmax = global maximum
    !!                O : ki,kj = global position of maximum
    !!
    !! ** Author : J.M. Molines 10/10/2004
    !!--------------------------------------------------------------------------
#ifdef key_mpp_shmem
    CALL ctl_stop( ' mpp_maxloc not yet available in SHMEM' )
# elif key_mpp_mpi
    !! * Arguments
    REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT (in)  :: ptab ,& ! Local 2D array
         &                                         pmask   ! Local mask
    REAL(wp)                     , INTENT (out) :: pmax    ! Global maximum of ptab
    INTEGER                      , INTENT (out) :: ki,kj,kk   ! index of maximum in global frame

    !! * Local variables
    REAL(wp) :: zmax   ! local maximum
    REAL(wp) ,DIMENSION(2,1) :: zain, zaout
    INTEGER, DIMENSION (3)  :: ilocs
    INTEGER :: ierror


    zmax  = MAXVAL( ptab(:,:,:) , mask= pmask == 1.e0 )
    ilocs = MAXLOC( ptab(:,:,:) , mask= pmask == 1.e0 )

    ki = ilocs(1) + nimpp - 1
    kj = ilocs(2) + njmpp - 1
    kk = ilocs(3)

    zain(1,:)=zmax
    zain(2,:)=ki+10000.*kj+100000000.*kk

    CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_OPA,ierror)

    pmax=zaout(1,1)
    kk= INT(zaout(2,1)/100000000.)
    kj= INT(zaout(2,1) - kk * 100000000. )/10000
    ki= INT(zaout(2,1) - kk * 100000000. -kj * 10000. )
#endif

  END SUBROUTINE mpp_maxloc3d

  SUBROUTINE mppsync()
    !!----------------------------------------------------------------------
    !!                  ***  routine mppsync  ***
    !!                   
    !! ** Purpose :   Massively parallel processors, synchroneous
    !!
    !!-----------------------------------------------------------------------

#if defined key_mpp_shmem

    !! * Local variables   (SHMEM version)
    CALL barrier()

#  elif defined key_mpp_mpi

    !! * Local variables   (MPI version)
    INTEGER :: ierror

    CALL mpi_barrier(mpi_comm_opa,ierror)

#endif

  END SUBROUTINE mppsync


  SUBROUTINE mppstop
    !!----------------------------------------------------------------------
    !!                  ***  routine mppstop  ***
    !!                   
    !! ** purpose :   Stop massilively parallel processors method
    !!
    !!----------------------------------------------------------------------
    !! * Local declarations
    INTEGER ::   info
    !!----------------------------------------------------------------------

    ! 1. Mpp synchroneus
    ! ------------------

    CALL mppsync
#if defined key_mpp_mpi
    CALL mpi_finalize( info )
#endif

  END SUBROUTINE mppstop


  SUBROUTINE mppobc( ptab, kd1, kd2, kl, kk, ktype, kij )
    !!----------------------------------------------------------------------
    !!                  ***  routine mppobc  ***
    !! 
    !! ** Purpose :   Message passing manadgement for open boundary
    !!     conditions array
    !!
    !! ** Method  :   Use mppsend and mpprecv function for passing mask
    !!       between processors following neighboring subdomains.
    !!       domain parameters
    !!                    nlci   : first dimension of the local subdomain
    !!                    nlcj   : second dimension of the local subdomain
    !!                    nbondi : mark for "east-west local boundary"
    !!                    nbondj : mark for "north-south local boundary"
    !!                    noea   : number for local neighboring processors 
    !!                    nowe   : number for local neighboring processors
    !!                    noso   : number for local neighboring processors
    !!                    nono   : number for local neighboring processors
    !!
    !! History :
    !!        !  98-07 (J.M. Molines) Open boundary conditions
    !!----------------------------------------------------------------------
    !! * Arguments
    INTEGER , INTENT( in ) ::   &
         kd1, kd2,   &  ! starting and ending indices
         kl ,        &  ! index of open boundary
         kk,         &  ! vertical dimension
         ktype,      &  ! define north/south or east/west cdt
         !              !  = 1  north/south  ;  = 2  east/west
         kij            ! horizontal dimension
    REAL(wp), DIMENSION(kij,kk), INTENT( inout )  ::   &
         ptab           ! variable array

    !! * Local variables
    INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
    INTEGER  ::   &
         iipt0, iipt1, ilpt1,     &  ! temporary integers
         ijpt0, ijpt1,            &  !    "          "
         imigr, iihom, ijhom         !    "          "
    INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
    INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
    REAL(wp), DIMENSION(jpi,jpj) ::   &
         ztab                        ! temporary workspace
    !!----------------------------------------------------------------------


    ! boundary condition initialization
    ! ---------------------------------

    ztab(:,:) = 0.e0

    IF( ktype==1 ) THEN                                  ! north/south boundaries
       iipt0 = MAX( 1, MIN(kd1 - nimpp+1, nlci     ) )
       iipt1 = MAX( 0, MIN(kd2 - nimpp+1, nlci - 1 ) )
       ilpt1 = MAX( 1, MIN(kd2 - nimpp+1, nlci     ) )
       ijpt0 = MAX( 1, MIN(kl  - njmpp+1, nlcj     ) )
       ijpt1 = MAX( 0, MIN(kl  - njmpp+1, nlcj - 1 ) )
    ELSEIF( ktype==2 ) THEN                              ! east/west boundaries
       iipt0 = MAX( 1, MIN(kl  - nimpp+1, nlci     ) )
       iipt1 = MAX( 0, MIN(kl  - nimpp+1, nlci - 1 ) )
       ijpt0 = MAX( 1, MIN(kd1 - njmpp+1, nlcj     ) )
       ijpt1 = MAX( 0, MIN(kd2 - njmpp+1, nlcj - 1 ) )
       ilpt1 = MAX( 1, MIN(kd2 - njmpp+1, nlcj     ) )
    ELSE
       CALL ctl_stop( 'mppobc: bad ktype' )
    ENDIF

    DO jk = 1, kk
       IF( ktype==1 ) THEN                               ! north/south boundaries
          DO jj = ijpt0, ijpt1
             DO ji = iipt0, iipt1
                ztab(ji,jj) = ptab(ji,jk)
             END DO
          END DO
       ELSEIF( ktype==2 ) THEN                           ! east/west boundaries
          DO jj = ijpt0, ijpt1
             DO ji = iipt0, iipt1
                ztab(ji,jj) = ptab(jj,jk)
             END DO
          END DO
       ENDIF


       ! 1. East and west directions
       ! ---------------------------

       ! 1.1 Read Dirichlet lateral conditions

       IF( nbondi /= 2 ) THEN
          iihom = nlci-nreci

          DO jl = 1, jpreci
             t2ew(:,jl,1) = ztab(jpreci+jl,:)
             t2we(:,jl,1) = ztab(iihom +jl,:)
          END DO
       ENDIF

       ! 1.2 Migrations

#if defined key_mpp_shmem
       !! *  (SHMEM version)
       imigr=jpreci*jpj*jpbyt

       IF( nbondi == -1 ) THEN
          CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr/jpbyt, noea )
       ELSEIF( nbondi == 0 ) THEN
          CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr/jpbyt, nowe )
          CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr/jpbyt, noea )
       ELSEIF( nbondi == 1 ) THEN
          CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr/jpbyt, nowe )
       ENDIF
       CALL barrier()
       CALL shmem_udcflush()

#  elif key_mpp_mpi
       !! * (MPI version)

       imigr=jpreci*jpj

       IF( nbondi == -1 ) THEN
          CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req1)
          CALL mpprecv(1,t2ew(1,1,2),imigr)
          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
       ELSEIF( nbondi == 0 ) THEN
          CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
          CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req2)
          CALL mpprecv(1,t2ew(1,1,2),imigr)
          CALL mpprecv(2,t2we(1,1,2),imigr)
          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
       ELSEIF( nbondi == 1 ) THEN
          CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
          CALL mpprecv(2,t2we(1,1,2),imigr)
          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
       ENDIF
#endif


       ! 1.3 Write Dirichlet lateral conditions

       iihom = nlci-jpreci
       IF( nbondi == 0 .OR. nbondi == 1 ) THEN
          DO jl = 1, jpreci
             ztab(jl,:) = t2we(:,jl,2)
          END DO
       ENDIF

       IF( nbondi == -1 .OR. nbondi == 0 ) THEN
          DO jl = 1, jpreci
             ztab(iihom+jl,:) = t2ew(:,jl,2)
          END DO
       ENDIF


       ! 2. North and south directions
       ! -----------------------------

       ! 2.1 Read Dirichlet lateral conditions

       IF( nbondj /= 2 ) THEN
          ijhom = nlcj-nrecj
          DO jl = 1, jprecj
             t2sn(:,jl,1) = ztab(:,ijhom +jl)
             t2ns(:,jl,1) = ztab(:,jprecj+jl)
          END DO
       ENDIF

       ! 2.2 Migrations

#if defined key_mpp_shmem
       !! * SHMEM version

       imigr=jprecj*jpi*jpbyt

       IF( nbondj == -1 ) THEN
          CALL shmem_put( t2sn(1,1,2), t2sn(1,1,1), imigr/jpbyt, nono )
       ELSEIF( nbondj == 0 ) THEN
          CALL shmem_put( t2ns(1,1,2), t2ns(1,1,1), imigr/jpbyt, noso )
          CALL shmem_put( t2sn(1,1,2), t2sn(1,1,1), imigr/jpbyt, nono )
       ELSEIF( nbondj == 1 ) THEN
          CALL shmem_put( t2ns(1,1,2), t2ns(1,1,1), imigr/jpbyt, noso )
       ENDIF
       CALL barrier()
       CALL shmem_udcflush()

#  elif key_mpp_mpi
       !! * Local variables   (MPI version)

       imigr=jprecj*jpi

       IF( nbondj == -1 ) THEN
          CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req1)
          CALL mpprecv(3,t2ns(1,1,2),imigr)
          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
       ELSEIF( nbondj == 0 ) THEN
          CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
          CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req2)
          CALL mpprecv(3,t2ns(1,1,2),imigr)
          CALL mpprecv(4,t2sn(1,1,2),imigr)
          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
       ELSEIF( nbondj == 1 ) THEN
          CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
          CALL mpprecv(4,t2sn(1,1,2),imigr)
          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
       ENDIF

#endif

       ! 2.3 Write Dirichlet lateral conditions

       ijhom = nlcj - jprecj
       IF( nbondj == 0 .OR. nbondj == 1 ) THEN
          DO jl = 1, jprecj
             ztab(:,jl) = t2sn(:,jl,2)
          END DO
       ENDIF

       IF( nbondj == 0 .OR. nbondj == -1 ) THEN
          DO jl = 1, jprecj
             ztab(:,ijhom+jl) = t2ns(:,jl,2)
          END DO
       ENDIF

       IF( ktype==1 .AND. kd1 <= jpi+nimpp-1 .AND. nimpp <= kd2 ) THEN
          ! north/south boundaries
          DO jj = ijpt0,ijpt1
             DO ji = iipt0,ilpt1
                ptab(ji,jk) = ztab(ji,jj)  
             END DO
          END DO
       ELSEIF( ktype==2 .AND. kd1 <= jpj+njmpp-1 .AND. njmpp <= kd2 ) THEN
          ! east/west boundaries
          DO jj = ijpt0,ilpt1
             DO ji = iipt0,iipt1
                ptab(jj,jk) = ztab(ji,jj) 
             END DO
          END DO
       ENDIF

    END DO

  END SUBROUTINE mppobc


  SUBROUTINE mpp_ini_north
    !!----------------------------------------------------------------------
    !!               ***  routine mpp_ini_north  ***
    !!
    !! ** Purpose :   Initialize special communicator for north folding 
    !!      condition together with global variables needed in the mpp folding
    !!
    !! ** Method  : - Look for northern processors
    !!              - Put their number in nrank_north
    !!              - Create groups for the world processors and the north processors
    !!              - Create a communicator for northern processors
    !!
    !! ** output
    !!      njmppmax = njmpp for northern procs
    !!      ndim_rank_north = number of processors in the northern line
    !!      nrank_north (ndim_rank_north) = number  of the northern procs.
    !!      ngrp_world = group ID for the world processors
    !!      ngrp_north = group ID for the northern processors
    !!      ncomm_north = communicator for the northern procs.
    !!      north_root = number (in the world) of proc 0 in the northern comm.
    !!
    !! History :
    !!        !  03-09 (J.M. Molines, MPI only )
    !!----------------------------------------------------------------------
#ifdef key_mpp_shmem
    CALL ctl_stop( ' mpp_ini_north not available in SHMEM' )
# elif key_mpp_mpi
    INTEGER :: ierr
    INTEGER :: jproc
    INTEGER :: ii,ji
    !!----------------------------------------------------------------------

    njmppmax=MAXVAL(njmppt)

    ! Look for how many procs on the northern boundary
    !
    ndim_rank_north=0
    DO jproc=1,jpnij
       IF ( njmppt(jproc) == njmppmax ) THEN
          ndim_rank_north = ndim_rank_north + 1
       END IF
    END DO


    ! Allocate the right size to nrank_north
    !
    ALLOCATE(nrank_north(ndim_rank_north))

    ! Fill the nrank_north array with proc. number of northern procs.
    ! Note : the rank start at 0 in MPI
    !
    ii=0
    DO ji = 1, jpnij
       IF ( njmppt(ji) == njmppmax   ) THEN
          ii=ii+1
          nrank_north(ii)=ji-1
       END IF
    END DO
    ! create the world group
    !
    CALL MPI_COMM_GROUP(mpi_comm_opa,ngrp_world,ierr)
    !
    ! Create the North group from the world group
    CALL MPI_GROUP_INCL(ngrp_world,ndim_rank_north,nrank_north,ngrp_north,ierr)

    ! Create the North communicator , ie the pool of procs in the north group
    !
    CALL MPI_COMM_CREATE(mpi_comm_opa,ngrp_north,ncomm_north,ierr)


    ! find proc number in the world of proc 0 in the north
    CALL MPI_GROUP_TRANSLATE_RANKS(ngrp_north,1,0,ngrp_world,north_root,ierr)
#endif

  END SUBROUTINE mpp_ini_north


   SUBROUTINE mpp_lbc_north_3d ( pt3d, cd_type, psgn )
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_3d  ***
      !!
      !! ** Purpose :
      !!      Ensure proper north fold horizontal bondary condition in mpp configuration
      !!      in case of jpn1 > 1
      !!
      !! ** Method :
      !!      Gather the 4 northern lines of the global domain on 1 processor and 
      !!      apply lbc north-fold on this sub array. Then scatter the fold array 
      !!      back to the processors.
      !!
      !! History :
      !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
      !!                                  from lbc routine
      !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1), INTENT( in ) ::   &
         cd_type       ! nature of pt3d grid-points
         !             !   = T ,  U , V , F or W  gridpoints
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
         pt3d          ! 3D array on which the boundary condition is applied
      REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
         !             !   = -1. , the sign is changed if north fold boundary
         !             !   =  1. , the sign is kept  if north fold boundary

      !! * Local declarations
      INTEGER :: ji, jj, jk, jr, jproc
      INTEGER :: ierr
      INTEGER :: ildi,ilei,iilb
      INTEGER :: ijpj,ijpjm1,ij,ijt,iju
      INTEGER :: itaille
      REAL(wp), DIMENSION(jpiglo,4,jpk) :: ztab
      REAL(wp), DIMENSION(jpi,4,jpk,jpni) :: znorthgloio
      REAL(wp), DIMENSION(jpi,4,jpk) :: znorthloc
      !!----------------------------------------------------------------------

    ! If we get in this routine it s because : North fold condition and mpp with more
    !   than one proc across i : we deal only with the North condition

    ! 0. Sign setting
    ! ---------------

    ijpj=4
    ijpjm1=3

    ! put in znorthloc the last 4 jlines of pt3d
    DO jk = 1, jpk 
       DO jj = nlcj - ijpj +1, nlcj
          ij = jj - nlcj + ijpj
          znorthloc(:,ij,jk) = pt3d(:,jj,jk)
       END DO
    END DO


    IF (npolj /= 0 ) THEN
       ! Build in proc 0 of ncomm_north the znorthgloio
       znorthgloio(:,:,:,:) = 0_wp

#ifdef key_mpp_shmem
       not done : compiler error
#elif defined key_mpp_mpi
       itaille=jpi*jpk*ijpj
       CALL MPI_GATHER(znorthloc,itaille,MPI_DOUBLE_PRECISION,znorthgloio,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
#endif

    ENDIF

    IF (narea == north_root+1 ) THEN
       ! recover the global north array
       ztab(:,:,:) = 0_wp

       DO jr = 1, ndim_rank_north
          jproc = nrank_north(jr) + 1
          ildi  = nldit (jproc)
          ilei  = nleit (jproc)
          iilb  = nimppt(jproc)
          DO jk = 1, jpk 
             DO jj = 1, 4
                DO ji = ildi, ilei
                   ztab(ji+iilb-1,jj,jk) = znorthgloio(ji,jj,jk,jr)
                END DO
             END DO
          END DO
       END DO


       ! Horizontal slab
       ! ===============

       DO jk = 1, jpk 


          ! 2. North-Fold boundary conditions
          ! ----------------------------------

          SELECT CASE ( npolj )

          CASE ( 3, 4 )                       ! *  North fold  T-point pivot

             ztab( 1    ,ijpj,jk) = 0.e0
             ztab(jpiglo,ijpj,jk) = 0.e0

             SELECT CASE ( cd_type )

             CASE ( 'T' , 'S' , 'W' )                   ! T-, W-point
                DO ji = 2, jpiglo
                   ijt = jpiglo-ji+2
                   ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-2,jk)
                END DO
                DO ji = jpiglo/2+1, jpiglo
                   ijt = jpiglo-ji+2
                   ztab(ji,ijpjm1,jk) = psgn * ztab(ijt,ijpjm1,jk)
                END DO

             CASE ( 'U' )                               ! U-point
                DO ji = 1, jpiglo-1
                   iju = jpiglo-ji+1
                   ztab(ji,ijpj,jk) = psgn * ztab(iju,ijpj-2,jk)
                END DO
                DO ji = jpiglo/2, jpiglo-1
                   iju = jpiglo-ji+1
                   ztab(ji,ijpjm1,jk) = psgn * ztab(iju,ijpjm1,jk)
                END DO

             CASE ( 'V' )                               ! V-point
                DO ji = 2, jpiglo
                   ijt = jpiglo-ji+2
                   ztab(ji,ijpj-1,jk) = psgn * ztab(ijt,ijpj-2,jk)
                   ztab(ji,ijpj  ,jk) = psgn * ztab(ijt,ijpj-3,jk)
                END DO

             CASE ( 'F' , 'G' )                         ! F-point
                DO ji = 1, jpiglo-1
                   iju = jpiglo-ji+1
                   ztab(ji,ijpj-1,jk) = psgn * ztab(iju,ijpj-2,jk)
                   ztab(ji,ijpj  ,jk) = psgn * ztab(iju,ijpj-3,jk)
                END DO

             END SELECT

          CASE ( 5, 6 )                        ! *  North fold  F-point pivot

             ztab( 1    ,ijpj,jk) = 0.e0
             ztab(jpiglo,ijpj,jk) = 0.e0

             SELECT CASE ( cd_type )

             CASE ( 'T' , 'S' , 'W' )                   ! T-, W-point
                DO ji = 1, jpiglo
                   ijt = jpiglo-ji+1
                   ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-1,jk)
                END DO

             CASE ( 'U' )                               ! U-point
                DO ji = 1, jpiglo-1
                   iju = jpiglo-ji
                   ztab(ji,ijpj,jk) = psgn * ztab(iju,ijpj-1,jk)
                END DO

             CASE ( 'V' )                               ! V-point
                DO ji = 1, jpiglo
                   ijt = jpiglo-ji+1
                   ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-2,jk)
                END DO
                DO ji = jpiglo/2+1, jpiglo
                   ijt = jpiglo-ji+1
                   ztab(ji,ijpjm1,jk) = psgn * ztab(ijt,ijpjm1,jk)
                END DO

             CASE ( 'F' , 'G' )                         ! F-point
                DO ji = 1, jpiglo-1
                   iju = jpiglo-ji
                   ztab(ji,ijpj  ,jk) = psgn * ztab(iju,ijpj-2,jk)
                END DO
                DO ji = jpiglo/2+1, jpiglo-1
                   iju = jpiglo-ji
                   ztab(ji,ijpjm1,jk) = psgn * ztab(iju,ijpjm1,jk)
                END DO

             END SELECT

          CASE DEFAULT                           ! *  closed

             SELECT CASE ( cd_type) 

             CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
                ztab(:, 1  ,jk) = 0.e0
                ztab(:,ijpj,jk) = 0.e0

             CASE ( 'F' )                               ! F-point
                ztab(:,ijpj,jk) = 0.e0

             END SELECT

          END SELECT

          !     End of slab
          !     ===========

       END DO

       !! Scatter back to pt3d
       DO jr = 1, ndim_rank_north
          jproc=nrank_north(jr)+1
          ildi=nldit (jproc)
          ilei=nleit (jproc)
          iilb=nimppt(jproc)
          DO jk=  1, jpk
             DO jj=1,ijpj
                DO ji=ildi,ilei
                   znorthgloio(ji,jj,jk,jr)=ztab(ji+iilb-1,jj,jk)
                END DO
             END DO
          END DO
       END DO

    ENDIF      ! only done on proc 0 of ncomm_north

#ifdef key_mpp_shmem
    not done yet in shmem : compiler error
#elif key_mpp_mpi
    IF ( npolj /= 0 ) THEN
       itaille=jpi*jpk*ijpj
       CALL MPI_SCATTER(znorthgloio,itaille,MPI_DOUBLE_PRECISION,znorthloc,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
    ENDIF
#endif

    ! put in the last ijpj jlines of pt3d znorthloc
    DO jk = 1 , jpk 
       DO jj = nlcj - ijpj + 1 , nlcj
          ij = jj - nlcj + ijpj
          pt3d(:,jj,jk)= znorthloc(:,ij,jk)
       END DO
    END DO

  END SUBROUTINE mpp_lbc_north_3d


  SUBROUTINE mpp_lbc_north_2d ( pt2d, cd_type, psgn)
    !!---------------------------------------------------------------------
    !!                   ***  routine mpp_lbc_north_2d  ***
    !!
    !! ** Purpose :
    !!      Ensure proper north fold horizontal bondary condition in mpp configuration
    !!      in case of jpn1 > 1 (for 2d array )
    !!
    !! ** Method :
    !!      Gather the 4 northern lines of the global domain on 1 processor and 
    !!      apply lbc north-fold on this sub array. Then scatter the fold array 
    !!      back to the processors.
    !!
    !! History :
    !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
    !!                                  from lbc routine
    !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
    !!----------------------------------------------------------------------

    !! * Arguments
    CHARACTER(len=1), INTENT( in ) ::   &
         cd_type       ! nature of pt2d grid-points
    !             !   = T ,  U , V , F or W  gridpoints
    REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
         pt2d          ! 2D array on which the boundary condition is applied
    REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
    !             !   = -1. , the sign is changed if north fold boundary
    !             !   =  1. , the sign is kept  if north fold boundary


    !! * Local declarations

    INTEGER :: ji, jj,  jr, jproc
    INTEGER :: ierr
    INTEGER :: ildi,ilei,iilb
    INTEGER :: ijpj,ijpjm1,ij,ijt,iju
    INTEGER :: itaille

    REAL(wp), DIMENSION(jpiglo,4) :: ztab
    REAL(wp), DIMENSION(jpi,4,jpni) :: znorthgloio
    REAL(wp), DIMENSION(jpi,4) :: znorthloc
    !!----------------------------------------------------------------------
    !!  OPA 8.5, LODYC-IPSL (2002)
    !!----------------------------------------------------------------------
    ! If we get in this routine it s because : North fold condition and mpp with more
    !   than one proc across i : we deal only with the North condition

    ! 0. Sign setting
    ! ---------------

    ijpj=4
    ijpjm1=3


    ! put in znorthloc the last 4 jlines of pt2d
    DO jj = nlcj - ijpj +1, nlcj
       ij = jj - nlcj + ijpj
       znorthloc(:,ij)=pt2d(:,jj)
    END DO

    IF (npolj /= 0 ) THEN
       ! Build in proc 0 of ncomm_north the znorthgloio
       znorthgloio(:,:,:) = 0_wp
#ifdef key_mpp_shmem
       not done : compiler error
#elif defined key_mpp_mpi
       itaille=jpi*ijpj
       CALL MPI_GATHER(znorthloc,itaille,MPI_DOUBLE_PRECISION,znorthgloio,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
#endif
    ENDIF

    IF (narea == north_root+1 ) THEN
       ! recover the global north array
       ztab(:,:) = 0_wp

       DO jr = 1, ndim_rank_north
          jproc=nrank_north(jr)+1
          ildi=nldit (jproc)
          ilei=nleit (jproc)
          iilb=nimppt(jproc)
          DO jj=1,4
             DO ji=ildi,ilei
                ztab(ji+iilb-1,jj)=znorthgloio(ji,jj,jr)
             END DO
          END DO
       END DO


       ! 2. North-Fold boundary conditions
       ! ----------------------------------

       SELECT CASE ( npolj )

       CASE ( 3, 4 )                       ! *  North fold  T-point pivot

          ztab( 1    ,ijpj) = 0.e0
          ztab(jpiglo,ijpj) = 0.e0

          SELECT CASE ( cd_type )

          CASE ( 'T' , 'W' , 'S' )                         ! T-, W-point
             DO ji = 2, jpiglo
                ijt = jpiglo-ji+2
                ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-2)
             END DO
             DO ji = jpiglo/2+1, jpiglo
                ijt = jpiglo-ji+2
                ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
             END DO

          CASE ( 'U' )                                     ! U-point
             DO ji = 1, jpiglo-1
                iju = jpiglo-ji+1
                ztab(ji,ijpj) = psgn * ztab(iju,ijpj-2)
             END DO
             DO ji = jpiglo/2, jpiglo-1
                iju = jpiglo-ji+1
                ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
             END DO

          CASE ( 'V' )                                     ! V-point
             DO ji = 2, jpiglo
                ijt = jpiglo-ji+2
                ztab(ji,ijpj-1) = psgn * ztab(ijt,ijpj-2)
                ztab(ji,ijpj  ) = psgn * ztab(ijt,ijpj-3)
             END DO

          CASE ( 'F' , 'G' )                               ! F-point
             DO ji = 1, jpiglo-1
                iju = jpiglo-ji+1
                ztab(ji,ijpj-1) = psgn * ztab(iju,ijpj-2)
                ztab(ji,ijpj  ) = psgn * ztab(iju,ijpj-3)
             END DO

          CASE ( 'I' )                                     ! ice U-V point
             ztab(2,ijpj) = psgn * ztab(3,ijpj-1)
             DO ji = 3, jpiglo
                iju = jpiglo - ji + 3
                ztab(ji,ijpj) = psgn * ztab(iju,ijpj-1)
             END DO

          END SELECT

       CASE ( 5, 6 )                        ! *  North fold  F-point pivot

          ztab( 1 ,ijpj) = 0.e0
          ztab(jpiglo,ijpj) = 0.e0

          SELECT CASE ( cd_type )

          CASE ( 'T' , 'W' ,'S' )                          ! T-, W-point
             DO ji = 1, jpiglo
                ijt = jpiglo-ji+1
                ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-1)
             END DO

          CASE ( 'U' )                                     ! U-point
             DO ji = 1, jpiglo-1
                iju = jpiglo-ji
                ztab(ji,ijpj) = psgn * ztab(iju,ijpj-1)
             END DO

          CASE ( 'V' )                                     ! V-point
             DO ji = 1, jpiglo
                ijt = jpiglo-ji+1
                ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-2)
             END DO
             DO ji = jpiglo/2+1, jpiglo
                ijt = jpiglo-ji+1
                ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
             END DO

          CASE ( 'F' , 'G' )                               ! F-point
             DO ji = 1, jpiglo-1
                iju = jpiglo-ji
                ztab(ji,ijpj  ) = psgn * ztab(iju,ijpj-2)
             END DO
             DO ji = jpiglo/2+1, jpiglo-1
                iju = jpiglo-ji
                ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
             END DO

             CASE ( 'I' )                                  ! ice U-V point
                ztab( 2 ,ijpj) = 0.e0
                DO ji = 2 , jpiglo-1
                   ijt = jpiglo - ji + 2
                   ztab(ji,ijpj)= 0.5 * ( ztab(ji,ijpj-1) + psgn * ztab(ijt,ijpj-1) )
                END DO

          END SELECT

       CASE DEFAULT                           ! *  closed : the code probably never go through

            SELECT CASE ( cd_type) 
  
            CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
               ztab(:, 1 ) = 0.e0
               ztab(:,ijpj) = 0.e0

            CASE ( 'F' )                                   ! F-point
               ztab(:,ijpj) = 0.e0

            CASE ( 'I' )                                   ! ice U-V point
               ztab(:, 1 ) = 0.e0
               ztab(:,ijpj) = 0.e0

            END SELECT

         END SELECT

         !     End of slab
         !     ===========

         !! Scatter back to pt2d
         DO jr = 1, ndim_rank_north
            jproc=nrank_north(jr)+1
            ildi=nldit (jproc)
            ilei=nleit (jproc)
            iilb=nimppt(jproc)
            DO jj=1,ijpj
               DO ji=ildi,ilei
                  znorthgloio(ji,jj,jr)=ztab(ji+iilb-1,jj)
               END DO
            END DO
         END DO

      ENDIF      ! only done on proc 0 of ncomm_north

#ifdef key_mpp_shmem
      not done yet in shmem : compiler error
#elif key_mpp_mpi
      IF ( npolj /= 0 ) THEN
         itaille=jpi*ijpj
         CALL MPI_SCATTER(znorthgloio,itaille,MPI_DOUBLE_PRECISION,znorthloc,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
      ENDIF
#endif

      ! put in the last ijpj jlines of pt2d znorthloc
      DO jj = nlcj - ijpj + 1 , nlcj
         ij = jj - nlcj + ijpj
         pt2d(:,jj)= znorthloc(:,ij)
      END DO

   END SUBROUTINE mpp_lbc_north_2d


   SUBROUTINE mpp_lbc_north_e ( pt2d, cd_type, psgn)
    !!---------------------------------------------------------------------
    !!                   ***  routine mpp_lbc_north_2d  ***
    !!
    !! ** Purpose :
    !!      Ensure proper north fold horizontal bondary condition in mpp configuration
    !!      in case of jpn1 > 1 (for 2d array with outer extra halo)
    !!
    !! ** Method :
    !!      Gather the 4+2*jpr2dj northern lines of the global domain on 1 processor and 
    !!      apply lbc north-fold on this sub array. Then scatter the fold array 
    !!      back to the processors.
    !!
    !! History :
    !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
    !!                                  from lbc routine
    !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
    !!   9.0  !  05-09  (R. Benshila )   adapt mpp_lbc_north_2d 
    !!----------------------------------------------------------------------

    !! * Arguments
    CHARACTER(len=1), INTENT( in ) ::   &
         cd_type       ! nature of pt2d grid-points
    !             !   = T ,  U , V , F or W  gridpoints
    REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj), INTENT( inout ) ::   &
         pt2d          ! 2D array on which the boundary condition is applied
    REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
    !             !   = -1. , the sign is changed if north fold boundary
    !             !   =  1. , the sign is kept  if north fold boundary


    !! * Local declarations

    INTEGER :: ji, jj,  jr, jproc, jl
    INTEGER :: ierr
    INTEGER :: ildi,ilei,iilb
    INTEGER :: ijpj,ijpjm1,ij,ijt,iju, iprecj
    INTEGER :: itaille

    REAL(wp), DIMENSION(jpiglo,1-jpr2dj:4+jpr2dj) :: ztab
    REAL(wp), DIMENSION(jpi,1-jpr2dj:4+jpr2dj,jpni) :: znorthgloio
    REAL(wp), DIMENSION(jpi,1-jpr2dj:4+jpr2dj) :: znorthloc

    ! If we get in this routine it s because : North fold condition and mpp with more
    !   than one proc across i : we deal only with the North condition

    ! 0. Sign setting
    ! ---------------

    ijpj=4
    ijpjm1=3
    iprecj = jpr2dj+jprecj

    ! put in znorthloc the last 4 jlines of pt2d
    DO jj = nlcj - ijpj + 1 - jpr2dj, nlcj +jpr2dj
       ij = jj - nlcj + ijpj
       znorthloc(:,ij)=pt2d(1:jpi,jj)
    END DO

    IF (npolj /= 0 ) THEN
       ! Build in proc 0 of ncomm_north the znorthgloio
       znorthgloio(:,:,:) = 0_wp
#ifdef key_mpp_shmem
       not done : compiler error
#elif defined key_mpp_mpi
       itaille=jpi*(ijpj+2*jpr2dj)
       CALL MPI_GATHER(znorthloc(1,1-jpr2dj),itaille,MPI_DOUBLE_PRECISION, &
                     & znorthgloio(1,1-jpr2dj,1),itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
#endif
    ENDIF

    IF (narea == north_root+1 ) THEN
       ! recover the global north array
       ztab(:,:) = 0_wp

       DO jr = 1, ndim_rank_north
          jproc=nrank_north(jr)+1
          ildi=nldit (jproc)
          ilei=nleit (jproc)
          iilb=nimppt(jproc)
          DO jj=1-jpr2dj,ijpj+jpr2dj
             DO ji=ildi,ilei
                ztab(ji+iilb-1,jj)=znorthgloio(ji,jj,jr)
             END DO
          END DO
       END DO


       ! 2. North-Fold boundary conditions
       ! ----------------------------------

       SELECT CASE ( npolj )

       CASE ( 3, 4 )                       ! *  North fold  T-point pivot

          ztab( 1    ,ijpj:ijpj+jpr2dj) = 0.e0
          ztab(jpiglo,ijpj:ijpj+jpr2dj) = 0.e0

          SELECT CASE ( cd_type )

          CASE ( 'T' , 'W' , 'S' )                         ! T-, W-point
             DO jl =0, iprecj-1
                DO ji = 2, jpiglo
                   ijt = jpiglo-ji+2
                   ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-2-jl)
                END DO
             END DO
             DO ji = jpiglo/2+1, jpiglo
                ijt = jpiglo-ji+2
                ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
             END DO

          CASE ( 'U' )                                     ! U-point
             DO jl =0, iprecj-1
                DO ji = 1, jpiglo-1
                   iju = jpiglo-ji+1
                   ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-2-jl)
                END DO
             END DO
             DO ji = jpiglo/2, jpiglo-1
                iju = jpiglo-ji+1
                ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
             END DO

          CASE ( 'V' )                                     ! V-point
            DO jl =-1, iprecj-1
               DO ji = 2, jpiglo
                  ijt = jpiglo-ji+2
                  ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-3-jl)
               END DO
            END DO

          CASE ( 'F' , 'G' )                               ! F-point
            DO jl =-1, iprecj-1
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-3-jl)
               END DO
             END DO

          CASE ( 'I' )                                     ! ice U-V point
             DO jl =0, iprecj-1
                ztab(2,ijpj+jl) = psgn * ztab(3,ijpj-1+jl)
                DO ji = 3, jpiglo
                   iju = jpiglo - ji + 3
                   ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-1-jl)
                END DO
             END DO

          END SELECT

       CASE ( 5, 6 )                        ! *  North fold  F-point pivot

          ztab( 1 ,ijpj:ijpj+jpr2dj) = 0.e0
          ztab(jpiglo,ijpj:ijpj+jpr2dj) = 0.e0

          SELECT CASE ( cd_type )

          CASE ( 'T' , 'W' ,'S' )                          ! T-, W-point
             DO jl = 0, iprecj-1
                DO ji = 1, jpiglo
                   ijt = jpiglo-ji+1
                   ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-1-jl)
                END DO
             END DO

          CASE ( 'U' )                                     ! U-point
             DO jl = 0, iprecj-1
                DO ji = 1, jpiglo-1
                   iju = jpiglo-ji
                   ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-1-jl)
                END DO
             END DO

          CASE ( 'V' )                                     ! V-point
             DO jl = 0, iprecj-1
                DO ji = 1, jpiglo
                   ijt = jpiglo-ji+1
                   ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-2-jl)
                END DO
             END DO
             DO ji = jpiglo/2+1, jpiglo
                ijt = jpiglo-ji+1
                ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
             END DO

          CASE ( 'F' , 'G' )                               ! F-point
             DO jl = 0, iprecj-1
                DO ji = 1, jpiglo-1
                   iju = jpiglo-ji
                   ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-2-jl)
                END DO
             END DO
             DO ji = jpiglo/2+1, jpiglo-1
                iju = jpiglo-ji
                ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
             END DO

             CASE ( 'I' )                                  ! ice U-V point
                ztab( 2 ,ijpj:ijpj+jpr2dj) = 0.e0
                DO jl = 0, jpr2dj
                   DO ji = 2 , jpiglo-1
                      ijt = jpiglo - ji + 2
                      ztab(ji,ijpj+jl)= 0.5 * ( ztab(ji,ijpj-1-jl) + psgn * ztab(ijt,ijpj-1-jl) )
                   END DO
                END DO

          END SELECT

       CASE DEFAULT                           ! *  closed : the code probably never go through

            SELECT CASE ( cd_type) 
  
            CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
               ztab(:, 1:1-jpr2dj     ) = 0.e0
               ztab(:,ijpj:ijpj+jpr2dj) = 0.e0

            CASE ( 'F' )                                   ! F-point
               ztab(:,ijpj:ijpj+jpr2dj) = 0.e0

            CASE ( 'I' )                                   ! ice U-V point
               ztab(:, 1:1-jpr2dj     ) = 0.e0
               ztab(:,ijpj:ijpj+jpr2dj) = 0.e0

            END SELECT

         END SELECT

         !     End of slab
         !     ===========

         !! Scatter back to pt2d
         DO jr = 1, ndim_rank_north
            jproc=nrank_north(jr)+1
            ildi=nldit (jproc)
            ilei=nleit (jproc)
            iilb=nimppt(jproc)
            DO jj=1-jpr2dj,ijpj+jpr2dj
               DO ji=ildi,ilei
                  znorthgloio(ji,jj,jr)=ztab(ji+iilb-1,jj)
               END DO
            END DO
         END DO

      ENDIF      ! only done on proc 0 of ncomm_north

#ifdef key_mpp_shmem
      not done yet in shmem : compiler error
#elif key_mpp_mpi
      IF ( npolj /= 0 ) THEN
         itaille=jpi*(ijpj+2*jpr2dj)
         CALL MPI_SCATTER(znorthgloio(1,1-jpr2dj,1),itaille,MPI_DOUBLE_PRECISION, &
                        & znorthloc(1,1-jpr2dj),itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
      ENDIF
#endif

      ! put in the last ijpj jlines of pt2d znorthloc
      DO jj = nlcj - ijpj  -jpr2dj + 1 , nlcj +jpr2dj
         ij = jj - nlcj + ijpj 
         pt2d(1:jpi,jj)= znorthloc(:,ij)
      END DO

   END SUBROUTINE mpp_lbc_north_e


   !!!!!


   !! 
   !!    This is valid on IBM machine ONLY. 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! mpi_init_opa.f90 : Redefinition du point d'entree MPI_INIT de la bibliotheque
   !!                MPI afin de faire, en plus de l'initialisation de
   !!                l'environnement MPI, l'allocation d'une zone tampon
   !!                qui sera ulterieurement utilisee automatiquement lors
   !!                de tous les envois de messages par MPI_BSEND
   !!
   !! Auteur : CNRS/IDRIS
   !! Date   : Tue Nov 13 12:02:14 2001
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE mpi_init_opa(code)
      IMPLICIT NONE

!$AGRIF_DO_NOT_TREAT
#     include <mpif.h>
!$AGRIF_END_DO_NOT_TREAT

      INTEGER                                 :: code,rang,ierr
      LOGICAL                                 :: mpi_was_called
 
      ! La valeur suivante doit etre au moins egale a la taille
      ! du plus grand message qui sera transfere dans le programme
      ! (de toute facon, il y aura un message d'erreur si cette
      ! valeur s'avere trop petite)
      INTEGER                                 :: taille_tampon
      CHARACTER(len=9)                        :: taille_tampon_alphanum
      REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: tampon
 
      ! Le point d'entree dans la bibliotheque MPI elle-meme
      CALL mpi_initialized(mpi_was_called, code)
      IF ( code /= MPI_SUCCESS ) THEN
        CALL ctl_stop( ' lib_mpp: Error in routine mpi_initialized' )
        CALL mpi_abort( mpi_comm_world, code, ierr )
      ENDIF

      IF ( .NOT. mpi_was_called ) THEN
         CALL mpi_init(code)
         CALL mpi_comm_dup( mpi_comm_world, mpi_comm_opa, code)
         IF ( code /= MPI_SUCCESS ) THEN
            CALL ctl_stop( ' lib_mpp: Error in routine mpi_comm_dup' )
            CALL mpi_abort( mpi_comm_world, code, ierr )
         ENDIF
      ENDIF
      ! La definition de la zone tampon pour les futurs envois
      ! par MPI_BSEND (on alloue une fois pour toute cette zone
      ! tampon, qui sera automatiquement utilisee lors de chaque
      ! appel  a MPI_BSEND).
      ! La desallocation sera implicite quand on sortira de
      ! l'environnement MPI.

      ! Recuperation de la valeur de la variable d'environnement
      ! BUFFER_LENGTH
      ! qui, si elle est definie, doit contenir une valeur superieure
      ! a  la taille en octets du plus gros message
      CALL getenv('BUFFER_LENGTH',taille_tampon_alphanum)
  
      ! Si la variable BUFFER_LENGTH n'est pas positionnee, on lui met par
      ! defaut la plus grande valeur de la variable MP_EAGER_LIMIT, soit
      ! 65 536 octets
      IF (taille_tampon_alphanum == ' ') THEN
         taille_tampon = 65536
      ELSE
         READ(taille_tampon_alphanum,'(i9)') taille_tampon
      END IF

      ! On est limite en mode d'adressage 32 bits a  1750 Mo pour la zone
      ! "data" soit 7 segments, c.-a -d. 1750/8 = 210 Mo
      IF (taille_tampon > 210000000) THEN
         CALL ctl_stop( ' lib_mpp: Attention la valeur BUFFER_LENGTH doit etre <= 210000000' )
         CALL mpi_abort(MPI_COMM_WORLD,2,code)
      END IF

      CALL mpi_comm_rank(MPI_COMM_OPA,rang,code)
      IF (rang == 0 ) PRINT *,'Taille du buffer alloue : ',taille_tampon

      ! Allocation du tampon et attachement
      ALLOCATE(tampon(taille_tampon))
      CALL mpi_buffer_attach(tampon,taille_tampon,code)

   END SUBROUTINE mpi_init_opa

!!DB 2008.04.24
!!Routine copied from  /export/home/wangz/NEMO_FOR_DAVE/.../OPA_SRC/IOM/in_out_manager.F90 
   SUBROUTINE ctl_stop( cd1, cd2, cd3, cd4, cd5,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_opa  ***
      !!
      !! ** Purpose : ??? blah blah....
      !!-----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!-----------------------------------------------------------------------
      !
      nstop = nstop + 1 
      IF(lwp) THEN
         WRITE(numout,"(/,' ===>>> : E R R O R',     /,'         ===========',/)") 
         IF( PRESENT(cd1 ) ) WRITE(numout,*) cd1
         IF( PRESENT(cd2 ) ) WRITE(numout,*) cd2
         IF( PRESENT(cd3 ) ) WRITE(numout,*) cd3
         IF( PRESENT(cd4 ) ) WRITE(numout,*) cd4
         IF( PRESENT(cd5 ) ) WRITE(numout,*) cd5
         IF( PRESENT(cd6 ) ) WRITE(numout,*) cd6
         IF( PRESENT(cd7 ) ) WRITE(numout,*) cd7
         IF( PRESENT(cd8 ) ) WRITE(numout,*) cd8
         IF( PRESENT(cd9 ) ) WRITE(numout,*) cd9
         IF( PRESENT(cd10) ) WRITE(numout,*) cd10
      ENDIF
      CALL FLUSH(numout)
      !
   END SUBROUTINE ctl_stop




#else
   !!----------------------------------------------------------------------
   !!   Default case:            Dummy module        share memory computing
   !!----------------------------------------------------------------------
   INTERFACE mpp_sum
      MODULE PROCEDURE mpp_sum_a2s, mpp_sum_as, mpp_sum_ai, mpp_sum_s, mpp_sum_i, mpp_sum_c
   END INTERFACE
   INTERFACE mpp_max
      MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
   END INTERFACE
   INTERFACE mpp_min
      MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
   END INTERFACE
   INTERFACE mpp_isl
      MODULE PROCEDURE mppisl_a_int, mppisl_int, mppisl_a_real, mppisl_real
   END INTERFACE
   INTERFACE mppobc
      MODULE PROCEDURE mppobc_1d, mppobc_2d, mppobc_3d, mppobc_4d
   END INTERFACE
  INTERFACE mpp_minloc
     MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
  END INTERFACE
  INTERFACE mpp_maxloc
     MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
  END INTERFACE


   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .FALSE.      !: mpp flag

CONTAINS

   FUNCTION mynode(localComm) RESULT (function_value)
      INTEGER, OPTIONAL :: localComm
      function_value = 0
   END FUNCTION mynode

   SUBROUTINE mppsync                       ! Dummy routine
   END SUBROUTINE mppsync

   SUBROUTINE mpp_sum_as( parr, kdim )      ! Dummy routine
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      WRITE(*,*) 'mpp_sum_as: You should not have seen this print! error?', kdim, parr(1)
   END SUBROUTINE mpp_sum_as

   SUBROUTINE mpp_sum_a2s( parr, kdim )      ! Dummy routine
      REAL   , DIMENSION(:,:) :: parr
      INTEGER               :: kdim
      WRITE(*,*) 'mpp_sum_a2s: You should not have seen this print! error?', kdim, parr(1,1)
   END SUBROUTINE mpp_sum_a2s

   SUBROUTINE mpp_sum_ai( karr, kdim )      ! Dummy routine
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      WRITE(*,*) 'mpp_sum_ai: You should not have seen this print! error?', kdim, karr(1)
   END SUBROUTINE mpp_sum_ai

   SUBROUTINE mpp_sum_s( psca )            ! Dummy routine
      REAL                  :: psca
      WRITE(*,*) 'mpp_sum_s: You should not have seen this print! error?', psca
   END SUBROUTINE mpp_sum_s

   SUBROUTINE mpp_sum_i( kint )            ! Dummy routine
      integer               :: kint
      WRITE(*,*) 'mpp_sum_i: You should not have seen this print! error?', kint
   END SUBROUTINE mpp_sum_i

   SUBROUTINE mpp_sum_c( kcmplx )            ! Dummy routine
      complex               :: kcmplx
      WRITE(*,*) 'mpp_sum_c: You should not have seen this print! error?', kint
   END SUBROUTINE mpp_sum_c

   SUBROUTINE mppmax_a_real( parr, kdim )
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      WRITE(*,*) 'mppmax_a_real: You should not have seen this print! error?', kdim, parr(1)
   END SUBROUTINE mppmax_a_real

   SUBROUTINE mppmax_real( psca )
      REAL                  :: psca
      WRITE(*,*) 'mppmax_real: You should not have seen this print! error?', psca
   END SUBROUTINE mppmax_real

   SUBROUTINE mppmin_a_real( parr, kdim )
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      WRITE(*,*) 'mppmin_a_real: You should not have seen this print! error?', kdim, parr(1)
   END SUBROUTINE mppmin_a_real

   SUBROUTINE mppmin_real( psca )
      REAL                  :: psca
      WRITE(*,*) 'mppmin_real: You should not have seen this print! error?', psca
   END SUBROUTINE mppmin_real

   SUBROUTINE mppmax_a_int( karr, kdim )
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      WRITE(*,*) 'mppmax_a_int: You should not have seen this print! error?', kdim, karr(1)
   END SUBROUTINE mppmax_a_int

   SUBROUTINE mppmax_int( kint )
      INTEGER               :: kint
      WRITE(*,*) 'mppmax_int: You should not have seen this print! error?', kint
   END SUBROUTINE mppmax_int

   SUBROUTINE mppmin_a_int( karr, kdim )
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      WRITE(*,*) 'mppmin_a_int: You should not have seen this print! error?', kdim, karr(1)
   END SUBROUTINE mppmin_a_int

   SUBROUTINE mppmin_int( kint )
      INTEGER               :: kint
      WRITE(*,*) 'mppmin_int: You should not have seen this print! error?', kint
   END SUBROUTINE mppmin_int

   SUBROUTINE mppobc_1d( parr, kd1, kd2, kl, kk, ktype, kij )
    INTEGER  ::   kd1, kd2, kl , kk, ktype, kij
    REAL, DIMENSION(:) ::   parr           ! variable array
      WRITE(*,*) 'mppobc: You should not have seen this print! error?',   &
         &        parr(1), kd1, kd2, kl, kk, ktype, kij
   END SUBROUTINE mppobc_1d

   SUBROUTINE mppobc_2d( parr, kd1, kd2, kl, kk, ktype, kij )
    INTEGER  ::   kd1, kd2, kl , kk, ktype, kij
    REAL, DIMENSION(:,:) ::   parr           ! variable array
      WRITE(*,*) 'mppobc: You should not have seen this print! error?',   &
         &        parr(1,1), kd1, kd2, kl, kk, ktype, kij
   END SUBROUTINE mppobc_2d

   SUBROUTINE mppobc_3d( parr, kd1, kd2, kl, kk, ktype, kij )
    INTEGER  ::   kd1, kd2, kl , kk, ktype, kij
    REAL, DIMENSION(:,:,:) ::   parr           ! variable array
      WRITE(*,*) 'mppobc: You should not have seen this print! error?',   &
         &        parr(1,1,1), kd1, kd2, kl, kk, ktype, kij
   END SUBROUTINE mppobc_3d

   SUBROUTINE mppobc_4d( parr, kd1, kd2, kl, kk, ktype, kij )
    INTEGER  ::   kd1, kd2, kl , kk, ktype, kij
    REAL, DIMENSION(:,:,:,:) ::   parr           ! variable array
      WRITE(*,*) 'mppobc: You should not have seen this print! error?',   &
         &        parr(1,1,1,1), kd1, kd2, kl, kk, ktype, kij
   END SUBROUTINE mppobc_4d


   SUBROUTINE mpplnks( parr )            ! Dummy routine
      REAL, DIMENSION(:,:) :: parr
      WRITE(*,*) 'mpplnks: You should not have seen this print! error?', parr(1,1)
   END SUBROUTINE mpplnks

   SUBROUTINE mppisl_a_int( karr, kdim )
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      WRITE(*,*) 'mppisl_a_int: You should not have seen this print! error?', kdim, karr(1)
   END SUBROUTINE mppisl_a_int

   SUBROUTINE mppisl_int( kint )
      INTEGER               :: kint
      WRITE(*,*) 'mppisl_int: You should not have seen this print! error?', kint
   END SUBROUTINE mppisl_int

   SUBROUTINE mppisl_a_real( parr, kdim )
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      WRITE(*,*) 'mppisl_a_real: You should not have seen this print! error?', kdim, parr(1)
   END SUBROUTINE mppisl_a_real

   SUBROUTINE mppisl_real( psca )
      REAL                  :: psca
      WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', psca
   END SUBROUTINE mppisl_real

   SUBROUTINE mpp_minloc2d ( ptab, pmask, pmin, ki, kj )
      REAL                   :: pmin
      REAL , DIMENSION (:,:) :: ptab, pmask
      INTEGER :: ki, kj
      WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', pmin, ki, kj
      WRITE(*,*) '   "      ":             "                 "            ', ptab(1,1), pmask(1,1)
   END SUBROUTINE mpp_minloc2d

   SUBROUTINE mpp_minloc3d ( ptab, pmask, pmin, ki, kj, kk )
      REAL                     :: pmin
      REAL , DIMENSION (:,:,:) :: ptab, pmask
      INTEGER :: ki, kj, kk
      WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', pmin, ki, kj, kk
      WRITE(*,*) '   "      ":             "                 "            ', ptab(1,1,1), pmask(1,1,1)
   END SUBROUTINE mpp_minloc3d

   SUBROUTINE mpp_maxloc2d ( ptab, pmask, pmax, ki, kj )
      REAL                   :: pmax
      REAL , DIMENSION (:,:) :: ptab, pmask
      INTEGER :: ki, kj
      WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', pmax, ki, kj
      WRITE(*,*) '   "      ":             "                 "            ', ptab(1,1), pmask(1,1)
   END SUBROUTINE mpp_maxloc2d

   SUBROUTINE mpp_maxloc3d ( ptab, pmask, pmax, ki, kj, kk )
      REAL                     :: pmax
      REAL , DIMENSION (:,:,:) :: ptab, pmask
      INTEGER :: ki, kj, kk
      WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', pmax, ki, kj, kk
      WRITE(*,*) '   "      ":             "                 "            ', ptab(1,1,1), pmask(1,1,1)
   END SUBROUTINE mpp_maxloc3d

   SUBROUTINE mppstop
      WRITE(*,*) 'mppstop: You should not have seen this print! error?'
   END SUBROUTINE mppstop

#endif
   !!----------------------------------------------------------------------
END MODULE lib_mpp
