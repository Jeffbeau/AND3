


PGF90 (Version     10.5)          09/23/2010  09:48:58      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir /usr/local/include
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /usr/local/include
          -idir /usr/include/mpich2-x86_64
          -idir /usr/include/mpich2-x86_64
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: lib_mpp.F90

(    1) MODULE lib_mpp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  lib_mpp  ***
(    4)    !! Ocean numerics:  massively parallel processing librairy
(    5)    !!=====================================================================
(    6) #if   defined key_mpp_mpi   ||   defined key_mpp_shmem
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_mpp_mpi'     OR      MPI massively parallel processing library
(    9)    !!   'key_mpp_shmem'         SHMEM massively parallel processing library
(   10)    !!----------------------------------------------------------------------
(   11)    !!   mynode
(   12)    !!   mpparent
(   13)    !!   mppshmem
(   14)    !!   mpp_lnk     : generic interface (defined in lbclnk) for :
(   15)    !!                 mpp_lnk_2d, mpp_lnk_3d
(   16)    !!   mpp_lnk_3d_gather :  Message passing manadgement for two 3D arrays
(   17)    !!   mpp_lnk_e   : interface defined in lbclnk
(   18)    !!   mpplnks
(   19)    !!   mpprecv
(   20)    !!   mppsend
(   21)    !!   mppscatter
(   22)    !!   mppgather
(   23)    !!   mpp_isl    : generic inteface  for :
(   24)    !!                mppisl_int , mppisl_a_int , mppisl_real, mppisl_a_real
(   25)    !!   mpp_min    : generic interface for : 
(   26)    !!                mppmin_int , mppmin_a_int , mppmin_real, mppmin_a_real
(   27)    !!   mpp_max    : generic interface for :
(   28)    !!                mppmax_int , mppmax_a_int , mppmax_real, mppmax_a_real
(   29)    !!   mpp_sum    : generic interface for :
(   30)    !!                mppsum_int , mppsum_a_int , mppsum_real, mppsum_a_real
(   31)    !!                mppsum_realdd, mppsum_a_realdd
(   32)    !!   mpp_minloc
(   33)    !!   mpp_maxloc
(   34)    !!   mppsync
(   35)    !!   mppstop
(   36)    !!   mppobc     : variant of mpp_lnk for open boundaries
(   37)    !!   mpp_ini_north
(   38)    !!   mpp_lbc_north
(   39)    !!   mpp_lbc_north_e : variant of mpp_lbc_north for extra outer halo (nsolv=4)
(   40)    !!----------------------------------------------------------------------
(   41)    !! History :
(   42)    !!        !  94 (M. Guyon, J. Escobar, M. Imbard)  Original code
(   43)    !!        !  97  (A.M. Treguier)  SHMEM additions
(   44)    !!        !  98  (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI
(   45)    !!   9.0  !  03  (J.-M. Molines, G. Madec)  F90, free form
(   46)    !!        !  04  (R. Bourdalle Badie)  isend option in mpi






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 2

(   47)    !!        !  05  (G. Madec, S. Masson)  npolj=5,6 F-point & ice cases
(   48)    !!        !  05  (R. Redler) Replacement of MPI_COMM_WORLD except for MPI_Abort
(   49)    !!----------------------------------------------------------------------
(   50)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   51)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_mpp.F90,v 1.21 2007/06/05 10:28:55 opalod Exp $ 
(   52)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   53)    !!---------------------------------------------------------------------
(   54)    !! * Modules used
(   55)    USE dom_oce                    ! ocean space and time domain 
(   56)    USE in_out_manager             ! I/O manager
(   57) 
(   58)    IMPLICIT NONE
(   59) 
(   60)    PRIVATE
(   61)    PUBLIC  mynode, mpparent, mpp_isl, mpp_min, mpp_max, mpp_sum,  mpp_lbc_north
(   62)    PUBLIC  mpp_lbc_north_e, mpp_minloc, mpp_maxloc, mpp_lnk_3d, mpp_lnk_2d, mpp_lnk_3d_gather, mpp_lnk_2d_e, mpplnks
(   63)    PUBLIC  mpprecv, mppsend, mppscatter, mppgather, mppobc, mpp_ini_north, mppstop, mppsync
(   64) #if defined key_oasis3 || defined key_oasis4
(   65)    PUBLIC  mppsize, mpprank
(   66) #endif
(   67) !!DB
(   68)    PUBLIC ctl_stop
(   69) 
(   70)    !! * Interfaces
(   71)    !! define generic interface for these routine as they are called sometimes
(   72)    !!        with scalar arguments instead of array arguments, which causes problems
(   73)    !!        for the compilation on AIX system as well as NEC and SGI. Ok on COMPACQ
(   74) 
(   75)    INTERFACE mpp_isl
(   76)       MODULE PROCEDURE mppisl_a_int, mppisl_int, mppisl_a_real, mppisl_real
(   77)    END INTERFACE
(   78)    INTERFACE mpp_min
(   79)       MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
(   80)    END INTERFACE
(   81)    INTERFACE mpp_max
(   82)       MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
(   83)    END INTERFACE
(   84)    INTERFACE mpp_sum
(   85)       MODULE PROCEDURE mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real, &
(   86)                        mppsum_realdd, mppsum_a_realdd
(   87)    END INTERFACE
(   88)    INTERFACE mpp_lbc_north
(   89)       MODULE PROCEDURE mpp_lbc_north_3d, mpp_lbc_north_2d 
(   90)    END INTERFACE
(   91)   INTERFACE mpp_minloc
(   92)      MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
(   93)   END INTERFACE
(   94)   INTERFACE mpp_maxloc
(   95)      MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
(   96)   END INTERFACE
(   97) 
(   98) 
(   99)    !! * Share module variables
(  100)    LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .TRUE.       !: mpp flag
(  101) 
(  102)    !! The processor number is a required power of two : 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,...
(  103)    INTEGER, PARAMETER ::   &
(  104)       nprocmax = 2**10     ! maximun dimension






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 3

(  105) 
(  106) #if defined key_mpp_mpi
(  107)    !! ========================= !!
(  108)    !!  MPI  variable definition !!
(  109)    !! ========================= !!
(  110) !$AGRIF_DO_NOT_TREAT
(  111) #  include <mpif.h>
(  112) !$AGRIF_END_DO_NOT_TREAT
(  113) 
(  114)    INTEGER ::   &
(  115)       mppsize,  &  ! number of process
(  116)       mpprank,  &  ! process number  [ 0 - size-1 ]
(  117)       mpi_comm_opa ! opa local communicator
(  118) 
(  119)    ! MD: needed to define MPI_SUMDD as a reduction operator
(  120)    INTEGER, PUBLIC :: MPI_SUMDD
(  121) 
(  122)    ! variables used in case of north fold condition in mpp_mpi with jpni > 1
(  123)    INTEGER ::      &       !
(  124)       ngrp_world,  &       ! group ID for the world processors
(  125)       ngrp_north,  &       ! group ID for the northern processors (to be fold)
(  126)       ncomm_north, &       ! communicator made by the processors belonging to ngrp_north
(  127)       ndim_rank_north, &   ! number of 'sea' processor in the northern line (can be /= jpni !)
(  128)       njmppmax             ! value of njmpp for the processors of the northern line
(  129)    INTEGER ::      &       !
(  130)       north_root           ! number (in the comm_opa) of proc 0 in the northern comm
(  131)    INTEGER, DIMENSION(:), ALLOCATABLE ::   &
(  132)       nrank_north          ! dimension ndim_rank_north, number of the procs belonging to ncomm_north
(  133)    CHARACTER (len=1) ::  &
(  134)       c_mpi_send = 'S'     ! type od mpi send/recieve (S=standard, B=bsend, I=isend)
(  135)    LOGICAL  ::           &
(  136)       l_isend = .FALSE.    ! isend use indicator (T if c_mpi_send='I')
(  137) !!DB: 
(  138)    CHARACTER(len=200) ::   ctmp1, ctmp2, ctmp3      !: temporary character
(  139) 
(  140) 
(  141) #elif defined key_mpp_shmem
(  142)    !! ========================= !!
(  143)    !! SHMEM variable definition !!
(  144)    !! ========================= !!
(  145) #  include  <fpvm3.h>
(  146) #  include <mpp/shmem.fh>
(  147) 
(  148)    CHARACTER (len=80), PARAMETER ::   simfile    = 'pvm3_ndim'   ! file name
(  149)    CHARACTER (len=47), PARAMETER ::   executable = 'opa'         ! executable name
(  150)    CHARACTER, PARAMETER ::            opaall     = ""            ! group name (old def opaall*(*))
(  151) 
(  152)    INTEGER, PARAMETER ::   & !! SHMEM control print
(  153)       mynode_print   = 0,  &  ! flag for print, mynode   routine
(  154)       mpprecv_print  = 0,  &  ! flag for print, mpprecv  routine
(  155)       mppsend_print  = 0,  &  ! flag for print, mppsend  routine
(  156)       mppsync_print  = 0,  &  ! flag for print, mppsync  routine
(  157)       mppsum_print   = 0,  &  ! flag for print, mpp_sum  routine
(  158)       mppisl_print   = 0,  &  ! flag for print, mpp_isl  routine
(  159)       mppmin_print   = 0,  &  ! flag for print, mpp_min  routine
(  160)       mppmax_print   = 0,  &  ! flag for print, mpp_max  routine
(  161)       mpparent_print = 0      ! flag for print, mpparent routine
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 4

(  163)    INTEGER, PARAMETER ::   & !! Variable definition
(  164)       jpvmint = 21            ! ???
(  165) 
(  166)    INTEGER, PARAMETER ::   & !! Maximum  dimension of array to sum on the processors
(  167)       jpmsec   = 50000,    &  ! ???
(  168)       jpmpplat =    30,    &  ! ???
(  169)       jpmppsum = MAX( jpisl*jpisl, jpmpplat*jpk, jpmsec )   ! ???
(  170) 
(  171)    INTEGER ::   &
(  172)       npvm_ipas ,  &  ! pvm initialization flag
(  173)       npvm_mytid,  &  ! pvm tid
(  174)       npvm_me   ,  &  ! node number [ 0 - nproc-1 ]
(  175)       npvm_nproc,  &  ! real number of nodes
(  176)       npvm_inum       ! ???
(  177)    INTEGER, DIMENSION(0:nprocmax-1) ::   &
(  178)       npvm_tids       ! tids array [ 0 - nproc-1 ]
(  179) 
(  180)    INTEGER ::   &
(  181)       nt3d_ipas ,  &  ! pvm initialization flag
(  182)       nt3d_mytid,  &  ! pvm tid
(  183)       nt3d_me   ,  &  ! node number [ 0 - nproc-1 ]
(  184)       nt3d_nproc      ! real number of nodes
(  185)    INTEGER, DIMENSION(0:nprocmax-1) ::   &
(  186)       nt3d_tids       ! tids array [ 0 - nproc-1 ]
(  187) 
(  188)    !! real sum reduction
(  189)    INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
(  190)        nrs1sync_shmem,   &  ! 
(  191)        nrs2sync_shmem
(  192)    REAL(wp), DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
(  193)        wrs1wrk_shmem,    &  !
(  194)        wrs2wrk_shmem        !
(  195)    REAL(wp), DIMENSION(jpmppsum) ::   &
(  196)        wrstab_shmem         !
(  197) 
(  198)    !! minimum and maximum reduction
(  199)    INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
(  200)        ni1sync_shmem,    &  ! 
(  201)        ni2sync_shmem        ! 
(  202)    REAL(wp), DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
(  203)        wi1wrk_shmem,     &  !
(  204)        wi2wrk_shmem
(  205)    REAL(wp), DIMENSION(jpmppsum) ::   &
(  206)        wintab_shmem,     &  ! 
(  207)        wi1tab_shmem,     &  ! 
(  208)        wi2tab_shmem         ! 
(  209)        
(  210)        !! value not equal zero for barotropic stream function around islands
(  211)    INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
(  212)        ni11sync_shmem,   &  !
(  213)        ni12sync_shmem,   &  !
(  214)        ni21sync_shmem,   &  !
(  215)        ni22sync_shmem       !
(  216)    REAL(wp), DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
(  217)        wi11wrk_shmem,    &  ! 
(  218)        wi12wrk_shmem,    &  !
(  219)        wi21wrk_shmem,    &  !
(  220)        wi22wrk_shmem        !






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 5

(  221)    REAL(wp), DIMENSION(jpmppsum) ::   &
(  222)        wiltab_shmem ,    &  !
(  223)        wi11tab_shmem,    &  !
(  224)        wi12tab_shmem,    &  ! 
(  225)        wi21tab_shmem,    &  ! 
(  226)        wi22tab_shmem
(  227) 
(  228)    INTEGER, DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
(  229)        ni11wrk_shmem,    &  !
(  230)        ni12wrk_shmem,    &  !
(  231)        ni21wrk_shmem,    &  !
(  232)        ni22wrk_shmem        !
(  233)    INTEGER, DIMENSION(jpmppsum) ::   &
(  234)        niitab_shmem ,    &  !
(  235)        ni11tab_shmem,    &  !
(  236)        ni12tab_shmem        !
(  237)    INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
(  238)        nis1sync_shmem,   &  !
(  239)        nis2sync_shmem       !
(  240)    INTEGER, DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
(  241)        nis1wrk_shmem,    &  ! 
(  242)        nis2wrk_shmem        !
(  243)    INTEGER, DIMENSION(jpmppsum) ::   &
(  244)        nistab_shmem
(  245) 
(  246)    !! integer sum reduction
(  247)    INTEGER, DIMENSION(SHMEM_REDUCE_SYNC_SIZE) ::   &
(  248)        nil1sync_shmem,   &  !
(  249)        nil2sync_shmem       !
(  250)    INTEGER, DIMENSION( MAX( SHMEM_REDUCE_MIN_WRKDATA_SIZE, jpmppsum/2+1 ) ) ::   &
(  251)        nil1wrk_shmem,    &  !
(  252)        nil2wrk_shmem        !
(  253)    INTEGER, DIMENSION(jpmppsum) ::   &
(  254)        niltab_shmem
(  255) #endif
(  256) 
(  257)    REAL(wp), DIMENSION(jpi,jprecj,jpk,2,2) ::   &
(  258)        t4ns, t4sn  ! 3d message passing arrays north-south & south-north
(  259)    REAL(wp), DIMENSION(jpj,jpreci,jpk,2,2) ::   &
(  260)        t4ew, t4we  ! 3d message passing arrays east-west & west-east
(  261)    REAL(wp), DIMENSION(jpi,jprecj,jpk,2,2) ::   &
(  262)        t4p1, t4p2  ! 3d message passing arrays north fold
(  263)    REAL(wp), DIMENSION(jpi,jprecj,jpk,2) ::   &
(  264)        t3ns, t3sn  ! 3d message passing arrays north-south & south-north
(  265)    REAL(wp), DIMENSION(jpj,jpreci,jpk,2) ::   &
(  266)        t3ew, t3we  ! 3d message passing arrays east-west & west-east
(  267)    REAL(wp), DIMENSION(jpi,jprecj,jpk,2) ::   &
(  268)        t3p1, t3p2  ! 3d message passing arrays north fold
(  269)    REAL(wp), DIMENSION(jpi,jprecj,2) ::   &
(  270)        t2ns, t2sn  ! 2d message passing arrays north-south & south-north
(  271)    REAL(wp), DIMENSION(jpj,jpreci,2) ::   &
(  272)        t2ew, t2we  ! 2d message passing arrays east-west & west-east
(  273)    REAL(wp), DIMENSION(jpi,jprecj,2) ::   &
(  274)        t2p1, t2p2  ! 2d message passing arrays north fold
(  275)    REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,jprecj+jpr2dj,2) ::   &
(  276)        tr2ns, tr2sn  ! 2d message passing arrays north-south & south-north including extra outer halo
(  277)    REAL(wp), DIMENSION(1-jpr2dj:jpj+jpr2dj,jpreci+jpr2di,2) ::   &
(  278)        tr2ew, tr2we  ! 2d message passing arrays east-west & west-east including extra outer halo






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 6

(  279)    !!----------------------------------------------------------------------
(  280)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(  281)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_mpp.F90,v 1.20 2007/03/07 13:26:44 opalod Exp $ 
(  282)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(  283)    !!---------------------------------------------------------------------
(  284) 
(  285) CONTAINS
(  286)    !!---------------------------------------------------------------------
(  287)    !!   Routine DDPDD_MPI: used by reduction operator MPI_SUMDD
(  288)    !!
(  289)    !!   Modification of original codes written by David H. Bailey
(  290)    !!   This subroutine computes ddb(i) = dda(i)+ddb(i)
(  291)    !!---------------------------------------------------------------------
(  292)       subroutine DDPDD_MPI (dda, ddb, len, itype)
(  293)       implicit none
(  294) 
(  295)     !! * Arguments
(  296)       INTEGER, INTENT(in)                         :: len, itype
(  297)       COMPLEX, DIMENSION(len), INTENT(in)     :: dda
(  298)       COMPLEX, DIMENSION(len), INTENT(inout)  :: ddb
(  299) 
(  300)     !! * Local variables
(  301)       REAL :: e, t1, t2  ! local work variables
(  302)       INTEGER :: i           ! local scalar
(  303) 
(  304)       DO i=1,len
(  305)    ! Compute dda + ddb using Knuth's trick.
(  306)         t1 = real(dda(i)) + real(ddb(i))
(  307)         e = t1 - real(dda(i))
(  308)         t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e))) &
(  309)              +imag(dda(i)) + imag(ddb(i))
(  310) 
(  311)    ! The result is t1 + t2, after normalization.
(  312)         ddb(i) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
(  313)       END DO
(  314)       END SUBROUTINE DDPDD_MPI
(  315)    !!---------------------------------------------------------------------
(  316) 
(  317) 
(  318)    FUNCTION mynode(localComm)
(  319)       !!----------------------------------------------------------------------
(  320)       !!                  ***  routine mynode  ***
(  321)       !!                    
(  322)       !! ** Purpose :   Find processor unit
(  323)       !!
(  324)       !!----------------------------------------------------------------------
(  325) #if defined key_mpp_mpi
(  326)       !! * Local variables   (MPI version)
(  327)       INTEGER ::   mynode, ierr, code
(  328)       LOGICAL ::   mpi_was_called
(  329)       INTEGER,OPTIONAL ::   localComm
(  330)       NAMELIST/nam_mpp/ c_mpi_send
(  331)       !!----------------------------------------------------------------------
(  332) 
(  333) !      WRITE(numout,*)
(  334) !      WRITE(numout,*) 'mynode : mpi initialisation'
(  335) !      WRITE(numout,*) '~~~~~~ '
(  336) !      WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 7

(  337) 
(  338)       ! Namelist namrun : parameters of the run
(  339) !      REWIND( numnam )
(  340) !      READ  ( numnam, nam_mpp )
(  341) 
(  342) !      WRITE(numout,*) '        Namelist nam_mpp'
(  343) !      WRITE(numout,*) '           mpi send type            c_mpi_send = ', c_mpi_send
(  344) 
(  345) ! MD
(  346) ! June 2006
(  347) ! MPI_INIT must be called before any files are opened!
(  348) ! So it's not possible to determine c_mpi_send from namelist.
(  349) ! Hence this is hacked to work with a hardcoded mpi send type:
(  350) 
(  351) c_mpi_send = 'I'
(  352) !c_mpi_send = 'I'
(  353) !c_mpi_send = 'B'
(  354) 
(  355) #if defined key_agrif
(  356)       IF( Agrif_Root() ) THEN
(  357) #endif
(  358) !!bug RB : should be clean to use Agrif in coupled mode
(  359) #if ! defined key_agrif
(  360)          CALL mpi_initialized ( mpi_was_called, code )
(  361)          IF( code /= MPI_SUCCESS ) THEN
(  362)             CALL ctl_stop( ' lib_mpp: Error in routine mpi_initialized' )
(  363)             CALL mpi_abort( mpi_comm_world, code, ierr )
(  364)          ENDIF
(  365) 
(  366)          IF( PRESENT(localComm) .and. mpi_was_called ) THEN
(  367)             mpi_comm_opa = localComm
(  368)             SELECT CASE ( c_mpi_send )
(  369)             CASE ( 'S' )                ! Standard mpi send (blocking)
(  370) !               WRITE(numout,*) '           Standard blocking mpi send (send)'
(  371)             CASE ( 'B' )                ! Buffer mpi send (blocking)
(  372) !               WRITE(numout,*) '           Buffer blocking mpi send (bsend)'
(  373)                CALL mpi_init_opa( ierr )
(  374)             CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
(  375) !               WRITE(numout,*) '           Immediate non-blocking send (isend)'
(  376)                l_isend = .TRUE.
(  377)             CASE DEFAULT
(  378) !               WRITE(numout,cform_err)
(  379) !               WRITE(numout,*) '           bad value for c_mpi_send = ', c_mpi_send
(  380)                nstop = nstop + 1
(  381)             END SELECT
(  382)          ELSE IF ( PRESENT(localComm) .and. .not. mpi_was_called ) THEN
(  383) !            WRITE(numout,*) ' lib_mpp: You cannot provide a local communicator '
(  384) !            WRITE(numout,*) '          without calling MPI_Init before ! '
(  385)          ELSE
(  386) #endif
(  387)             SELECT CASE ( c_mpi_send )
(  388)             CASE ( 'S' )                ! Standard mpi send (blocking)
(  389) !               WRITE(numout,*) '           Standard blocking mpi send (send)'
(  390)                CALL mpi_init( ierr )
(  391)             CASE ( 'B' )                ! Buffer mpi send (blocking)
(  392) !               WRITE(numout,*) '           Buffer blocking mpi send (bsend)'
(  393)                CALL mpi_init_opa( ierr )
(  394)             CASE ( 'I' )                ! Immediate mpi send (non-blocking send)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 8

(  395) !               WRITE(numout,*) '           Immediate non-blocking send (isend)'
(  396)                l_isend = .TRUE.
(  397)                CALL mpi_init( ierr )
(  398)             CASE DEFAULT
(  399) !               WRITE(ctmp1,*) '           bad value for c_mpi_send = ', c_mpi_send
(  400)                CALL ctl_stop( ctmp1 )
(  401)             END SELECT
(  402) 
(  403) #if ! defined key_agrif
(  404)             CALL mpi_comm_dup( mpi_comm_world, mpi_comm_opa, code)
(  405)             IF( code /= MPI_SUCCESS ) THEN
(  406)                CALL ctl_stop( ' lib_mpp: Error in routine mpi_comm_dup' )
(  407)                CALL mpi_abort( mpi_comm_world, code, ierr )
(  408)             ENDIF
(  409)             !
(  410)          ENDIF
(  411) #endif
(  412) #if defined key_agrif
(  413)       ELSE
(  414)          SELECT CASE ( c_mpi_send )
(  415)          CASE ( 'S' )                ! Standard mpi send (blocking)
(  416) !            WRITE(numout,*) '           Standard blocking mpi send (send)'
(  417)          CASE ( 'B' )                ! Buffer mpi send (blocking)
(  418) !            WRITE(numout,*) '           Buffer blocking mpi send (bsend)'
(  419)          CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
(  420) !            WRITE(numout,*) '           Immediate non-blocking send (isend)'
(  421)             l_isend = .TRUE.
(  422)          CASE DEFAULT
(  423) !            WRITE(numout,cform_err)
(  424) !            WRITE(numout,*) '           bad value for c_mpi_send = ', c_mpi_send
(  425)             nstop = nstop + 1
(  426)          END SELECT
(  427)       ENDIF
(  428) 
(  429)       mpi_comm_opa = mpi_comm_world
(  430) #endif
(  431) 
(  432)       CALL mpi_comm_rank( mpi_comm_opa, mpprank, ierr )
(  433)       CALL mpi_comm_size( mpi_comm_opa, mppsize, ierr )
(  434)       mynode = mpprank
(  435) 
(  436) !MD: Create MPISUMDD
(  437)       call MPI_OP_CREATE(DDPDD_MPI, .TRUE., MPI_SUMDD, ierr)
(  438) 
(  439) #else
(  440)       !! * Local variables   (SHMEM version)
(  441)       INTEGER ::   mynode
(  442)       INTEGER ::   &
(  443)            imypid, imyhost, ji, info, iparent_tid
(  444)       !!----------------------------------------------------------------------
(  445) 
(  446)       IF( npvm_ipas /= nprocmax ) THEN
(  447)          !         ---   first passage in mynode
(  448)          !         -------------
(  449)          !         enroll in pvm
(  450)          !         -------------
(  451)          CALL pvmfmytid( npvm_mytid )
(  452)          IF( mynode_print /= 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 9

(  453)             WRITE(numout,*) 'mynode, npvm_ipas =', npvm_ipas, ' nprocmax=', nprocmax
(  454)             WRITE(numout,*) 'mynode, npvm_mytid=', npvm_mytid, ' after pvmfmytid'
(  455)          ENDIF
(  456) 
(  457)          !         ---------------------------------------------------------------
(  458)          !         find out IF i am parent or child spawned processes have parents
(  459)          !         ---------------------------------------------------------------
(  460)          CALL mpparent( iparent_tid )
(  461)          IF( mynode_print /= 0 ) THEN
(  462)             WRITE(numout,*) 'mynode, npvm_mytid=', npvm_mytid,   &
(  463)                &            ' after mpparent, npvm_tids(0) = ',   &
(  464)                &            npvm_tids(0), ' iparent_tid=', iparent_tid
(  465)          ENDIF
(  466)          IF( iparent_tid < 0 )  THEN
(  467)             WRITE(numout,*) 'mynode, npvm_mytid=', npvm_mytid,   &
(  468)                &            ' after mpparent, npvm_tids(0) = ',   &
(  469)                &            npvm_tids(0), ' iparent_tid=', iparent_tid
(  470)             npvm_tids(0) = npvm_mytid
(  471)             npvm_me = 0
(  472)             IF( jpnij > nprocmax ) THEN
(  473)                WRITE(ctmp1,*) 'npvm_mytid=', npvm_mytid, ' too great'
(  474)                CALL ctl_stop( ctmp1 )
(  475) 
(  476)             ELSE
(  477)                npvm_nproc = jpnij
(  478)             ENDIF
(  479) 
(  480)             ! -------------------------
(  481)             ! start up copies of myself
(  482)             ! -------------------------
(  483)             IF( npvm_nproc > 1 ) THEN
(  484)                DO ji = 1, npvm_nproc-1
(  485)                   npvm_tids(ji) = nt3d_tids(ji)
(  486)                END DO
(  487)                info=npvm_nproc-1
(  488)   
(  489)                IF( mynode_print /= 0 ) THEN
(  490)                   WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid,   &
(  491)                      &            ' maitre=',executable,' info=', info   &
(  492)                      &            ,' npvm_nproc=',npvm_nproc
(  493)                   WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid,   &
(  494)                      &            ' npvm_tids ',(npvm_tids(ji),ji=0,npvm_nproc-1)
(  495)                ENDIF
(  496) 
(  497)                ! ---------------------------
(  498)                ! multicast tids array to children
(  499)                ! ---------------------------
(  500)                CALL pvmfinitsend( pvmdefault, info )
(  501)                CALL pvmfpack ( jpvmint, npvm_nproc, 1         , 1, info )
(  502)                CALL pvmfpack ( jpvmint, npvm_tids , npvm_nproc, 1, info )
(  503)                CALL pvmfmcast( npvm_nproc-1, npvm_tids(1), 10, info )
(  504)             ENDIF
(  505)          ELSE
(  506) 
(  507)             ! ---------------------------------
(  508)             ! receive the tids array and set me
(  509)             ! ---------------------------------
(  510)             IF( mynode_print /= 0 )   WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid, ' pvmfrecv'






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 10

(  511)             CALL pvmfrecv( iparent_tid, 10, info )
(  512)             IF( mynode_print /= 0 )   WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid, " fin pvmfrecv"
(  513)             CALL pvmfunpack( jpvmint, npvm_nproc, 1         , 1, info )
(  514)             CALL pvmfunpack( jpvmint, npvm_tids , npvm_nproc, 1, info )
(  515)             IF( mynode_print /= 0 ) THEN
(  516)                WRITE(numout,*) 'mynode, npvm_mytid=',npvm_mytid,   &
(  517)                   &            ' esclave=', executable,' info=', info,' npvm_nproc=',npvm_nproc
(  518)                WRITE(numout,*) 'mynode, npvm_mytid=', npvm_mytid,   &
(  519)                   &            'npvm_tids', ( npvm_tids(ji), ji = 0, npvm_nproc-1 )
(  520)             ENDIF
(  521)             DO ji = 0, npvm_nproc-1
(  522)                IF( npvm_mytid == npvm_tids(ji) ) npvm_me = ji
(  523)             END DO
(  524)          ENDIF
(  525) 
(  526)          ! ------------------------------------------------------------
(  527)          ! all nproc tasks are equal now
(  528)          ! and can address each other by tids(0) thru tids(nproc-1)
(  529)          ! for each process me => process number [0-(nproc-1)]
(  530)          ! ------------------------------------------------------------
(  531)          CALL pvmfjoingroup ( "bidon", info )
(  532)          CALL pvmfbarrier   ( "bidon", npvm_nproc, info )
(  533)          DO ji = 0, npvm_nproc-1
(  534)             IF( ji == npvm_me ) THEN
(  535)                CALL pvmfjoingroup ( opaall, npvm_inum )
(  536)                IF( npvm_inum /= npvm_me )   WRITE(numout,*) 'mynode not arrived in the good order for opaall'
(  537)             ENDIF
(  538)             CALL pvmfbarrier( "bidon", npvm_nproc, info )
(  539)          END DO
(  540)          CALL pvmfbarrier( opaall, npvm_nproc, info )
(  541)   
(  542)       ELSE
(  543)          ! ---   other passage in mynode
(  544)       ENDIF
(  545)  
(  546)       npvm_ipas = nprocmax
(  547)       mynode    = npvm_me
(  548)       imypid    = npvm_mytid
(  549)       imyhost   = npvm_tids(0)
(  550)       IF( mynode_print /= 0 ) THEN
(  551)          WRITE(numout,*)'mynode: npvm_mytid=', npvm_mytid, ' npvm_me=', npvm_me,   &
(  552)             &           ' npvm_nproc=', npvm_nproc , ' npvm_ipas=', npvm_ipas
(  553)       ENDIF
(  554) #endif
(  555)    END FUNCTION mynode
(  556) 
(  557) 
(  558)    SUBROUTINE mpparent( kparent_tid )
(  559)       !!----------------------------------------------------------------------
(  560)       !!                  ***  routine mpparent  ***
(  561)       !!
(  562)       !! ** Purpose :   use an pvmfparent routine for T3E (key_mpp_shmem)
(  563)       !!              or  only return -1 (key_mpp_mpi)
(  564)       !!----------------------------------------------------------------------
(  565)       !! * Arguments
(  566)       INTEGER, INTENT(inout) ::   kparent_tid      ! ???
(  567)   
(  568) #if defined key_mpp_mpi






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 11

(  569)       ! MPI version : retour -1
(  570) 
(  571)       kparent_tid = -1
(  572) 
(  573) #else
(  574)       !! * Local variables   (SHMEN onto T3E version)
(  575)       INTEGER ::   &
(  576)            it3d_my_pe, LEADZ, ji, info
(  577)   
(  578)       CALL pvmfmytid( nt3d_mytid )
(  579)       CALL pvmfgetpe( nt3d_mytid, it3d_my_pe )
(  580)       IF( mpparent_print /= 0 ) THEN
(  581)          WRITE(numout,*) 'mpparent: nt3d_mytid= ', nt3d_mytid ,' it3d_my_pe=',it3d_my_pe
(  582)       ENDIF
(  583)       IF( it3d_my_pe == 0 ) THEN
(  584)          !-----------------------------------------------------------------!
(  585)          !     process = 0 => receive other tids                           !
(  586)          !-----------------------------------------------------------------!
(  587)          kparent_tid = -1
(  588)          IF(mpparent_print /= 0 ) THEN
(  589)             WRITE(numout,*) 'mpparent, nt3d_mytid=',nt3d_mytid ,' kparent_tid=',kparent_tid
(  590)          ENDIF
(  591)          !          --- END receive dimension ---
(  592)          IF( jpnij > nprocmax ) THEN
(  593)             WRITE(ctmp1,*) 'mytid=',nt3d_mytid,' too great'
(  594)             CALL ctl_stop( ctmp1 )
(  595)          ELSE
(  596)             nt3d_nproc =  jpnij
(  597)          ENDIF
(  598)          IF( mpparent_print /= 0 ) THEN
(  599)             WRITE(numout,*) 'mpparent, nt3d_mytid=', nt3d_mytid , ' nt3d_nproc=', nt3d_nproc
(  600)          ENDIF
(  601)          !-------- receive tids from others process --------
(  602)          DO ji = 1, nt3d_nproc-1
(  603)             CALL pvmfrecv( ji , 100, info )
(  604)             CALL pvmfunpack( jpvmint, nt3d_tids(ji), 1, 1, info )
(  605)             IF( mpparent_print /= 0 ) THEN
(  606)                WRITE(numout,*) 'mpparent, nt3d_mytid=', nt3d_mytid , ' receive=', nt3d_tids(ji), ' from = ', ji
(  607)             ENDIF
(  608)          END DO
(  609)          nt3d_tids(0) = nt3d_mytid
(  610)          IF( mpparent_print /= 0 ) THEN
(  611)             WRITE(numout,*) 'mpparent, nt3d_mytid=', nt3d_mytid , ' nt3d_tids(ji) =', (nt3d_tids(ji),   &
(  612)                  ji = 0, nt3d_nproc-1 )
(  613)             WRITE(numout,*) 'mpparent, nt3d_mytid=', nt3d_mytid , ' kparent_tid=', kparent_tid
(  614)          ENDIF
(  615) 
(  616)       ELSE
(  617)          !!----------------------------------------------------------------!
(  618)          !     process <> 0 => send  other tids                            !
(  619)          !!----------------------------------------------------------------!
(  620)          kparent_tid = 0
(  621)          CALL pvmfinitsend( pvmdataraw, info )
(  622)          CALL pvmfpack( jpvmint, nt3d_mytid, 1, 1, info )
(  623)          CALL pvmfsend( kparent_tid, 100, info )
(  624)       ENDIF
(  625) #endif
(  626) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 12

(  627)    END SUBROUTINE mpparent
(  628) 
(  629) #if defined key_mpp_shmem
(  630) 
(  631)    SUBROUTINE mppshmem
(  632)       !!----------------------------------------------------------------------
(  633)       !!                  ***  routine mppshmem  ***
(  634)       !!
(  635)       !! ** Purpose :   SHMEM ROUTINE
(  636)       !!
(  637)       !!----------------------------------------------------------------------
(  638)       nrs1sync_shmem = SHMEM_SYNC_VALUE
(  639)       nrs2sync_shmem = SHMEM_SYNC_VALUE
(  640)       nis1sync_shmem = SHMEM_SYNC_VALUE
(  641)       nis2sync_shmem = SHMEM_SYNC_VALUE
(  642)       nil1sync_shmem = SHMEM_SYNC_VALUE
(  643)       nil2sync_shmem = SHMEM_SYNC_VALUE
(  644)       ni11sync_shmem = SHMEM_SYNC_VALUE
(  645)       ni12sync_shmem = SHMEM_SYNC_VALUE
(  646)       ni21sync_shmem = SHMEM_SYNC_VALUE
(  647)       ni22sync_shmem = SHMEM_SYNC_VALUE
(  648)       CALL barrier()
(  649)   
(  650)    END SUBROUTINE mppshmem
(  651) 
(  652) #endif
(  653) 
(  654)    SUBROUTINE mpp_lnk_3d( ptab, cd_type, psgn, cd_mpp )
(  655)       !!----------------------------------------------------------------------
(  656)       !!                  ***  routine mpp_lnk_3d  ***
(  657)       !!
(  658)       !! ** Purpose :   Message passing manadgement
(  659)       !!
(  660)       !! ** Method  :   Use mppsend and mpprecv function for passing mask 
(  661)       !!      between processors following neighboring subdomains.
(  662)       !!            domain parameters
(  663)       !!                    nlci   : first dimension of the local subdomain
(  664)       !!                    nlcj   : second dimension of the local subdomain
(  665)       !!                    nbondi : mark for "east-west local boundary"
(  666)       !!                    nbondj : mark for "north-south local boundary"
(  667)       !!                    noea   : number for local neighboring processors 
(  668)       !!                    nowe   : number for local neighboring processors
(  669)       !!                    noso   : number for local neighboring processors
(  670)       !!                    nono   : number for local neighboring processors
(  671)       !!
(  672)       !! ** Action  :   ptab with update value at its periphery
(  673)       !!
(  674)       !!----------------------------------------------------------------------
(  675)       !! * Arguments
(  676)       CHARACTER(len=1) , INTENT( in ) ::   &
(  677)          cd_type       ! define the nature of ptab array grid-points
(  678)          !             ! = T , U , V , F , W points
(  679)          !             ! = S : T-point, north fold treatment ???
(  680)          !             ! = G : F-point, north fold treatment ???
(  681)       REAL(wp), INTENT( in ) ::   &
(  682)          psgn          ! control of the sign change
(  683)          !             !   = -1. , the sign is changed if north fold boundary
(  684)          !             !   =  1. , the sign is kept  if north fold boundary






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 13

(  685)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
(  686)          ptab          ! 3D array on which the boundary condition is applied
(  687)       CHARACTER(len=3), INTENT( in ), OPTIONAL ::    &
(  688)          cd_mpp        ! fill the overlap area only 
(  689) 
(  690)       !! * Local variables
(  691)       INTEGER ::   ji, jj, jk, jl                        ! dummy loop indices
(  692)       INTEGER ::   imigr, iihom, ijhom, iloc, ijt, iju   ! temporary integers
(  693)       INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
(  694)       INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
(  695)       !!----------------------------------------------------------------------
(  696) 
(  697)       ! 1. standard boundary treatment
(  698)       ! ------------------------------
(  699) 
(  700)       IF( PRESENT( cd_mpp ) ) THEN
(  701)          DO jj = nlcj+1, jpj   ! only fill extra allows last line
(  702)             ptab(1:nlci, jj, :) = ptab(1:nlci, nlej, :)
(  703)          END DO
(  704)          DO ji = nlci+1, jpi   ! only fill extra allows last column
(  705)             ptab(ji    , : , :) = ptab(nlei  , :   , :)
(  706)          END DO
(  707)       ELSE      
(  708) 
(  709)          !                                        ! East-West boundaries
(  710)          !                                        ! ====================
(  711)          IF( nbondi == 2 .AND.   &      ! Cyclic east-west
(  712)             &   (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
(  713)             ptab( 1 ,:,:) = ptab(jpim1,:,:)
(  714)             ptab(jpi,:,:) = ptab(  2  ,:,:)
(  715) 
(  716)          ELSE                           ! closed
(  717)             SELECT CASE ( cd_type )
(  718)             CASE ( 'T', 'U', 'V', 'W' )
(  719)                ptab(     1       :jpreci,:,:) = 0.e0
(  720)                ptab(nlci-jpreci+1:jpi   ,:,:) = 0.e0
(  721)             CASE ( 'F' )
(  722)                ptab(nlci-jpreci+1:jpi   ,:,:) = 0.e0
(  723)             END SELECT 
(  724)          ENDIF
(  725) 
(  726)          !                                        ! North-South boundaries
(  727)          !                                        ! ======================
(  728)          SELECT CASE ( cd_type )
(  729)          CASE ( 'T', 'U', 'V', 'W' )
(  730)             ptab(:,     1       :jprecj,:) = 0.e0
(  731)             ptab(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
(  732)          CASE ( 'F' )
(  733)             ptab(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
(  734)          END SELECT
(  735)      
(  736)       ENDIF
(  737) 
(  738)       ! 2. East and west directions exchange
(  739)       ! ------------------------------------
(  740) 
(  741)       ! 2.1 Read Dirichlet lateral conditions
(  742) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 14

(  743)       SELECT CASE ( nbondi )
(  744)       CASE ( -1, 0, 1 )    ! all exept 2 
(  745)          iihom = nlci-nreci
(  746)          DO jl = 1, jpreci
(  747)             t3ew(:,jl,:,1) = ptab(jpreci+jl,:,:)
(  748)             t3we(:,jl,:,1) = ptab(iihom +jl,:,:)
(  749)          END DO
(  750)       END SELECT
(  751) 
(  752)       ! 2.2 Migrations
(  753) 
(  754) #if defined key_mpp_shmem
(  755)       !! * SHMEM version
(  756) 
(  757)       imigr = jpreci * jpj * jpk
(  758) 
(  759)       SELECT CASE ( nbondi )
(  760)       CASE ( -1 )
(  761)          CALL shmem_put( t3we(1,1,1,2), t3we(1,1,1,1), imigr, noea )
(  762)       CASE ( 0 )
(  763)          CALL shmem_put( t3ew(1,1,1,2), t3ew(1,1,1,1), imigr, nowe )
(  764)          CALL shmem_put( t3we(1,1,1,2), t3we(1,1,1,1), imigr, noea )
(  765)       CASE ( 1 )
(  766)          CALL shmem_put( t3ew(1,1,1,2), t3ew(1,1,1,1), imigr, nowe )
(  767)       END SELECT
(  768) 
(  769)       CALL barrier()
(  770)       CALL shmem_udcflush()
(  771) 
(  772) #elif defined key_mpp_mpi
(  773)       !! * Local variables   (MPI version)
(  774) 
(  775)       imigr = jpreci * jpj * jpk
(  776) 
(  777)       SELECT CASE ( nbondi ) 
(  778)       CASE ( -1 )
(  779)          CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req1 )
(  780)          CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
(  781)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  782)       CASE ( 0 )
(  783)          CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
(  784)          CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req2 )
(  785)          CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
(  786)          CALL mpprecv( 2, t3we(1,1,1,2), imigr )
(  787)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  788)          IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
(  789)       CASE ( 1 )
(  790)          CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
(  791)          CALL mpprecv( 2, t3we(1,1,1,2), imigr )
(  792)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  793)       END SELECT
(  794) #endif
(  795) 
(  796)       ! 2.3 Write Dirichlet lateral conditions
(  797) 
(  798)       iihom = nlci-jpreci
(  799) 
(  800)       SELECT CASE ( nbondi )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 15

(  801)       CASE ( -1 )
(  802)          DO jl = 1, jpreci
(  803)             ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
(  804)          END DO
(  805)       CASE ( 0 ) 
(  806)          DO jl = 1, jpreci
(  807)             ptab(jl      ,:,:) = t3we(:,jl,:,2)
(  808)             ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
(  809)          END DO
(  810)       CASE ( 1 )
(  811)          DO jl = 1, jpreci
(  812)             ptab(jl      ,:,:) = t3we(:,jl,:,2)
(  813)          END DO
(  814)       END SELECT
(  815) 
(  816) 
(  817)       ! 3. North and south directions
(  818)       ! -----------------------------
(  819) 
(  820)       ! 3.1 Read Dirichlet lateral conditions
(  821) 
(  822)       IF( nbondj /= 2 ) THEN
(  823)          ijhom = nlcj-nrecj
(  824)          DO jl = 1, jprecj
(  825)             t3sn(:,jl,:,1) = ptab(:,ijhom +jl,:)
(  826)             t3ns(:,jl,:,1) = ptab(:,jprecj+jl,:)
(  827)          END DO
(  828)       ENDIF
(  829) 
(  830)       ! 3.2 Migrations
(  831) 
(  832) #if defined key_mpp_shmem
(  833)       !! * SHMEM version
(  834) 
(  835)       imigr = jprecj * jpi * jpk
(  836) 
(  837)       SELECT CASE ( nbondj )
(  838)       CASE ( -1 )
(  839)          CALL shmem_put( t3sn(1,1,1,2), t3sn(1,1,1,1), imigr, nono )
(  840)       CASE ( 0 )
(  841)          CALL shmem_put( t3ns(1,1,1,2), t3ns(1,1,1,1), imigr, noso )
(  842)          CALL shmem_put( t3sn(1,1,1,2), t3sn(1,1,1,1), imigr, nono )
(  843)       CASE ( 1 )
(  844)          CALL shmem_put( t3ns(1,1,1,2), t3ns(1,1,1,1), imigr, noso )
(  845)       END SELECT
(  846) 
(  847)       CALL barrier()
(  848)       CALL shmem_udcflush()
(  849) 
(  850) #elif defined key_mpp_mpi
(  851)       !! * Local variables   (MPI version)
(  852)   
(  853)       imigr=jprecj*jpi*jpk
(  854) 
(  855)       SELECT CASE ( nbondj )     
(  856)       CASE ( -1 )
(  857)          CALL mppsend( 4, t3sn(1,1,1,1), imigr, nono, ml_req1 )
(  858)          CALL mpprecv( 3, t3ns(1,1,1,2), imigr )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 16

(  859)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  860)       CASE ( 0 )
(  861)          CALL mppsend( 3, t3ns(1,1,1,1), imigr, noso, ml_req1 )
(  862)          CALL mppsend( 4, t3sn(1,1,1,1), imigr, nono, ml_req2 )
(  863)          CALL mpprecv( 3, t3ns(1,1,1,2), imigr )
(  864)          CALL mpprecv( 4, t3sn(1,1,1,2), imigr )
(  865)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  866)          IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
(  867)       CASE ( 1 ) 
(  868)          CALL mppsend( 3, t3ns(1,1,1,1), imigr, noso, ml_req1 )
(  869)          CALL mpprecv( 4, t3sn(1,1,1,2), imigr )
(  870)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  871)       END SELECT
(  872) 
(  873) #endif
(  874) 
(  875)       ! 3.3 Write Dirichlet lateral conditions
(  876) 
(  877)       ijhom = nlcj-jprecj
(  878) 
(  879)       SELECT CASE ( nbondj )
(  880)       CASE ( -1 )
(  881)          DO jl = 1, jprecj
(  882)             ptab(:,ijhom+jl,:) = t3ns(:,jl,:,2)
(  883)          END DO
(  884)       CASE ( 0 ) 
(  885)          DO jl = 1, jprecj
(  886)             ptab(:,jl      ,:) = t3sn(:,jl,:,2)
(  887)             ptab(:,ijhom+jl,:) = t3ns(:,jl,:,2)
(  888)          END DO
(  889)       CASE ( 1 )
(  890)          DO jl = 1, jprecj
(  891)             ptab(:,jl,:) = t3sn(:,jl,:,2)
(  892)          END DO
(  893)       END SELECT
(  894) 
(  895) 
(  896)       ! 4. north fold treatment
(  897)       ! -----------------------
(  898) 
(  899)       IF (PRESENT(cd_mpp)) THEN
(  900)          ! No north fold treatment (it is assumed to be already OK)
(  901)      
(  902)       ELSE      
(  903) 
(  904)       ! 4.1 treatment without exchange (jpni odd)
(  905)       !     T-point pivot  
(  906) 
(  907)       SELECT CASE ( jpni )
(  908) 
(  909)       CASE ( 1 )  ! only one proc along I, no mpp exchange
(  910) 
(  911)          SELECT CASE ( npolj )
(  912)   
(  913)          CASE ( 3 , 4 )    ! T pivot
(  914)             iloc = jpiglo - 2 * ( nimpp - 1 )
(  915) 
(  916)             SELECT CASE ( cd_type )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 17

(  917) 
(  918)             CASE ( 'T' , 'S', 'W' )
(  919)                DO jk = 1, jpk
(  920)                   DO ji = 2, nlci
(  921)                      ijt=iloc-ji+2
(  922)                      ptab(ji,nlcj,jk) = psgn * ptab(ijt,nlcj-2,jk)
(  923)                   END DO
(  924)                   DO ji = nlci/2+1, nlci
(  925)                      ijt=iloc-ji+2
(  926)                      ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-1,jk)
(  927)                   END DO
(  928)                END DO
(  929)           
(  930)             CASE ( 'U' )
(  931)                DO jk = 1, jpk
(  932)                   DO ji = 1, nlci-1
(  933)                      iju=iloc-ji+1
(  934)                      ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-2,jk)
(  935)                   END DO
(  936)                   DO ji = nlci/2, nlci-1
(  937)                      iju=iloc-ji+1
(  938)                      ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-1,jk)
(  939)                   END DO
(  940)                END DO
(  941) 
(  942)             CASE ( 'V' )
(  943)                DO jk = 1, jpk
(  944)                   DO ji = 2, nlci
(  945)                      ijt=iloc-ji+2
(  946)                      ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-2,jk)
(  947)                      ptab(ji,nlcj  ,jk) = psgn * ptab(ijt,nlcj-3,jk)
(  948)                   END DO
(  949)                END DO
(  950) 
(  951)             CASE ( 'F', 'G' )
(  952)                DO jk = 1, jpk
(  953)                   DO ji = 1, nlci-1
(  954)                      iju=iloc-ji+1
(  955)                      ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-2,jk)
(  956)                      ptab(ji,nlcj  ,jk) = psgn * ptab(iju,nlcj-3,jk)
(  957)                   END DO
(  958)                END DO
(  959)   
(  960)           END SELECT
(  961)        
(  962)          CASE ( 5 , 6 ) ! F pivot
(  963)             iloc=jpiglo-2*(nimpp-1)
(  964)   
(  965)             SELECT CASE ( cd_type )
(  966) 
(  967)             CASE ( 'T' , 'S', 'W' )
(  968)                DO jk = 1, jpk
(  969)                   DO ji = 1, nlci
(  970)                      ijt=iloc-ji+1
(  971)                      ptab(ji,nlcj,jk) = psgn * ptab(ijt,nlcj-1,jk)
(  972)                   END DO
(  973)                END DO
(  974) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 18

(  975)             CASE ( 'U' )
(  976)                DO jk = 1, jpk
(  977)                   DO ji = 1, nlci-1
(  978)                      iju=iloc-ji
(  979)                      ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-1,jk)
(  980)                   END DO
(  981)                END DO
(  982) 
(  983)             CASE ( 'V' )
(  984)                DO jk = 1, jpk
(  985)                   DO ji = 1, nlci
(  986)                      ijt=iloc-ji+1
(  987)                      ptab(ji,nlcj  ,jk) = psgn * ptab(ijt,nlcj-2,jk)
(  988)                   END DO
(  989)                   DO ji = nlci/2+1, nlci
(  990)                      ijt=iloc-ji+1
(  991)                      ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-1,jk)
(  992)                   END DO
(  993)                END DO
(  994) 
(  995)             CASE ( 'F', 'G' )
(  996)                DO jk = 1, jpk
(  997)                   DO ji = 1, nlci-1
(  998)                      iju=iloc-ji
(  999)                      ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-2,jk)
( 1000)                   END DO
( 1001)                   DO ji = nlci/2+1, nlci-1
( 1002)                      iju=iloc-ji
( 1003)                      ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-1,jk)
( 1004)                   END DO
( 1005)                END DO
( 1006)             END SELECT  ! cd_type
( 1007) 
( 1008)          END SELECT     !  npolj
( 1009)   
( 1010)       CASE DEFAULT ! more than 1 proc along I
( 1011)          IF ( npolj /= 0 ) CALL mpp_lbc_north (ptab, cd_type, psgn)  ! only for northern procs.
( 1012) 
( 1013)       END SELECT ! jpni 
( 1014) 
( 1015)       ENDIF
( 1016)       
( 1017) 
( 1018)       ! 5. East and west directions exchange
( 1019)       ! ------------------------------------
( 1020) 
( 1021)       SELECT CASE ( npolj )
( 1022) 
( 1023)       CASE ( 3, 4, 5, 6 )
( 1024) 
( 1025)          ! 5.1 Read Dirichlet lateral conditions
( 1026) 
( 1027)          SELECT CASE ( nbondi )
( 1028) 
( 1029)          CASE ( -1, 0, 1 )
( 1030)             iihom = nlci-nreci
( 1031)             DO jl = 1, jpreci
( 1032)                t3ew(:,jl,:,1) = ptab(jpreci+jl,:,:)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 19

( 1033)                t3we(:,jl,:,1) = ptab(iihom +jl,:,:)
( 1034)             END DO
( 1035) 
( 1036)          END SELECT
( 1037) 
( 1038)          ! 5.2 Migrations
( 1039) 
( 1040) #if defined key_mpp_shmem
( 1041)          !! SHMEM version
( 1042) 
( 1043)          imigr = jpreci * jpj * jpk
( 1044) 
( 1045)          SELECT CASE ( nbondi )
( 1046)          CASE ( -1 )
( 1047)             CALL shmem_put( t3we(1,1,1,2), t3we(1,1,1,1), imigr, noea )
( 1048)          CASE ( 0 )
( 1049)             CALL shmem_put( t3ew(1,1,1,2), t3ew(1,1,1,1), imigr, nowe )
( 1050)             CALL shmem_put( t3we(1,1,1,2), t3we(1,1,1,1), imigr, noea )
( 1051)          CASE ( 1 )
( 1052)             CALL shmem_put( t3ew(1,1,1,2), t3ew(1,1,1,1), imigr, nowe )
( 1053)          END SELECT
( 1054) 
( 1055)          CALL barrier()
( 1056)          CALL shmem_udcflush()
( 1057) 
( 1058) #elif defined key_mpp_mpi
( 1059)          !! MPI version
( 1060) 
( 1061)          imigr=jpreci*jpj*jpk
( 1062)   
( 1063)          SELECT CASE ( nbondi )
( 1064)          CASE ( -1 )
( 1065)             CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req1 )
( 1066)             CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
( 1067)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1068)          CASE ( 0 )
( 1069)             CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
( 1070)             CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req2 )
( 1071)             CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
( 1072)             CALL mpprecv( 2, t3we(1,1,1,2), imigr )
( 1073)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1074)             IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 1075)          CASE ( 1 )
( 1076)             CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
( 1077)             CALL mpprecv( 2, t3we(1,1,1,2), imigr )
( 1078)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1079)          END SELECT
( 1080) #endif
( 1081) 
( 1082)          ! 5.3 Write Dirichlet lateral conditions
( 1083) 
( 1084)          iihom = nlci-jpreci
( 1085) 
( 1086)          SELECT CASE ( nbondi)
( 1087)          CASE ( -1 )
( 1088)             DO jl = 1, jpreci
( 1089)                ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
( 1090)             END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 20

( 1091)          CASE ( 0 ) 
( 1092)             DO jl = 1, jpreci
( 1093)                ptab(jl      ,:,:) = t3we(:,jl,:,2)
( 1094)                ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
( 1095)             END DO
( 1096)          CASE ( 1 )
( 1097)             DO jl = 1, jpreci
( 1098)                ptab(jl      ,:,:) = t3we(:,jl,:,2)
( 1099)             END DO
( 1100)          END SELECT
( 1101) 
( 1102)       END SELECT    ! npolj 
( 1103) 
( 1104)    END SUBROUTINE mpp_lnk_3d
( 1105) 
( 1106) 
( 1107)    SUBROUTINE mpp_lnk_2d( pt2d, cd_type, psgn, cd_mpp )
( 1108)       !!----------------------------------------------------------------------
( 1109)       !!                  ***  routine mpp_lnk_2d  ***
( 1110)       !!                  
( 1111)       !! ** Purpose :   Message passing manadgement for 2d array
( 1112)       !!
( 1113)       !! ** Method  :   Use mppsend and mpprecv function for passing mask 
( 1114)       !!      between processors following neighboring subdomains.
( 1115)       !!            domain parameters
( 1116)       !!                    nlci   : first dimension of the local subdomain
( 1117)       !!                    nlcj   : second dimension of the local subdomain
( 1118)       !!                    nbondi : mark for "east-west local boundary"
( 1119)       !!                    nbondj : mark for "north-south local boundary"
( 1120)       !!                    noea   : number for local neighboring processors 
( 1121)       !!                    nowe   : number for local neighboring processors
( 1122)       !!                    noso   : number for local neighboring processors
( 1123)       !!                    nono   : number for local neighboring processors
( 1124)       !!
( 1125)       !!----------------------------------------------------------------------
( 1126)       !! * Arguments
( 1127)       CHARACTER(len=1) , INTENT( in ) ::   &
( 1128)          cd_type       ! define the nature of pt2d array grid-points
( 1129)          !             !  = T , U , V , F , W 
( 1130)          !             !  = S : T-point, north fold treatment
( 1131)          !             !  = G : F-point, north fold treatment
( 1132)          !             !  = I : sea-ice velocity at F-point with index shift
( 1133)       REAL(wp), INTENT( in ) ::   &
( 1134)          psgn          ! control of the sign change
( 1135)          !             !   = -1. , the sign is changed if north fold boundary
( 1136)          !             !   =  1. , the sign is kept  if north fold boundary
( 1137)       REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
( 1138)          pt2d          ! 2D array on which the boundary condition is applied
( 1139)       CHARACTER(len=3), INTENT( in ), OPTIONAL ::    &
( 1140)          cd_mpp        ! fill the overlap area only 
( 1141) 
( 1142)       !! * Local variables
( 1143)       INTEGER  ::   ji, jj, jl      ! dummy loop indices
( 1144)       INTEGER  ::   &
( 1145)          imigr, iihom, ijhom,    &  ! temporary integers
( 1146)          iloc, ijt, iju             !    "          "
( 1147)       INTEGER  ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
( 1148)       INTEGER  ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 21

( 1149)       !!----------------------------------------------------------------------
( 1150) 
( 1151)       ! 1. standard boundary treatment
( 1152)       ! ------------------------------
( 1153)       IF (PRESENT(cd_mpp)) THEN
( 1154)          DO jj = nlcj+1, jpj   ! only fill extra allows last line
( 1155)             pt2d(1:nlci, jj) = pt2d(1:nlci, nlej)
( 1156)          END DO
( 1157)          DO ji = nlci+1, jpi   ! only fill extra allows last column
( 1158)             pt2d(ji    , : ) = pt2d(nlei  , :   )
( 1159)          END DO     
( 1160)       ELSE      
( 1161) 
( 1162)          !                                        ! East-West boundaries
( 1163)          !                                        ! ====================
( 1164)          IF( nbondi == 2 .AND.   &      ! Cyclic east-west
( 1165)             &    (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
( 1166)             pt2d( 1 ,:) = pt2d(jpim1,:)
( 1167)             pt2d(jpi,:) = pt2d(  2  ,:)
( 1168) 
( 1169)          ELSE                           ! ... closed
( 1170)             SELECT CASE ( cd_type )
( 1171)             CASE ( 'T', 'U', 'V', 'W' , 'I' )
( 1172)                pt2d(     1       :jpreci,:) = 0.e0
( 1173)                pt2d(nlci-jpreci+1:jpi   ,:) = 0.e0
( 1174)             CASE ( 'F' )
( 1175)                pt2d(nlci-jpreci+1:jpi   ,:) = 0.e0
( 1176)             END SELECT
( 1177)          ENDIF
( 1178) 
( 1179)          !                                        ! North-South boundaries
( 1180)          !                                        ! ======================
( 1181)          SELECT CASE ( cd_type )
( 1182)          CASE ( 'T', 'U', 'V', 'W' , 'I' )
( 1183)             pt2d(:,     1       :jprecj) = 0.e0
( 1184)             pt2d(:,nlcj-jprecj+1:jpj   ) = 0.e0
( 1185)          CASE ( 'F' )
( 1186)             pt2d(:,nlcj-jprecj+1:jpj   ) = 0.e0
( 1187)          END SELECT
( 1188) 
( 1189)       ENDIF
( 1190) 
( 1191) 
( 1192)       ! 2. East and west directions
( 1193)       ! ---------------------------
( 1194) 
( 1195)       ! 2.1 Read Dirichlet lateral conditions
( 1196) 
( 1197)       SELECT CASE ( nbondi )
( 1198)       CASE ( -1, 0, 1 )    ! all except 2
( 1199)          iihom = nlci-nreci
( 1200)          DO jl = 1, jpreci
( 1201)             t2ew(:,jl,1) = pt2d(jpreci+jl,:)
( 1202)             t2we(:,jl,1) = pt2d(iihom +jl,:)
( 1203)          END DO
( 1204)       END SELECT
( 1205) 
( 1206)       ! 2.2 Migrations






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 22

( 1207) 
( 1208) #if defined key_mpp_shmem
( 1209)       !! * SHMEM version
( 1210) 
( 1211)       imigr = jpreci * jpj
( 1212) 
( 1213)       SELECT CASE ( nbondi )
( 1214)       CASE ( -1 )
( 1215)          CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr, noea )
( 1216)       CASE ( 0 )
( 1217)          CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr, nowe )
( 1218)          CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr, noea )
( 1219)       CASE ( 1 )
( 1220)          CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr, nowe )
( 1221)       END SELECT
( 1222) 
( 1223)       CALL barrier()
( 1224)       CALL shmem_udcflush()
( 1225) 
( 1226) #elif defined key_mpp_mpi
( 1227)       !! * MPI version
( 1228) 
( 1229)       imigr = jpreci * jpj
( 1230) 
( 1231)       SELECT CASE ( nbondi )
( 1232)       CASE ( -1 )
( 1233)          CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req1 )
( 1234)          CALL mpprecv( 1, t2ew(1,1,2), imigr )
( 1235)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1236)       CASE ( 0 )
( 1237)          CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
( 1238)          CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req2 )
( 1239)          CALL mpprecv( 1, t2ew(1,1,2), imigr )
( 1240)          CALL mpprecv( 2, t2we(1,1,2), imigr )
( 1241)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1242)          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 1243)       CASE ( 1 )
( 1244)          CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
( 1245)          CALL mpprecv( 2, t2we(1,1,2), imigr )
( 1246)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1247)       END SELECT
( 1248) 
( 1249) #endif
( 1250) 
( 1251)       ! 2.3 Write Dirichlet lateral conditions
( 1252) 
( 1253)       iihom = nlci - jpreci
( 1254)       SELECT CASE ( nbondi )
( 1255)       CASE ( -1 )
( 1256)          DO jl = 1, jpreci
( 1257)             pt2d(iihom+jl,:) = t2ew(:,jl,2)
( 1258)          END DO
( 1259)       CASE ( 0 )
( 1260)          DO jl = 1, jpreci
( 1261)             pt2d(jl      ,:) = t2we(:,jl,2)
( 1262)             pt2d(iihom+jl,:) = t2ew(:,jl,2)
( 1263)          END DO
( 1264)       CASE ( 1 )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 23

( 1265)          DO jl = 1, jpreci
( 1266)             pt2d(jl      ,:) = t2we(:,jl,2)
( 1267)          END DO
( 1268)       END SELECT
( 1269) 
( 1270) 
( 1271)       ! 3. North and south directions
( 1272)       ! -----------------------------
( 1273) 
( 1274)       ! 3.1 Read Dirichlet lateral conditions
( 1275) 
( 1276)       IF( nbondj /= 2 ) THEN
( 1277)          ijhom = nlcj-nrecj
( 1278)          DO jl = 1, jprecj
( 1279)             t2sn(:,jl,1) = pt2d(:,ijhom +jl)
( 1280)             t2ns(:,jl,1) = pt2d(:,jprecj+jl)
( 1281)          END DO
( 1282)       ENDIF
( 1283) 
( 1284)       ! 3.2 Migrations
( 1285) 
( 1286) #if defined key_mpp_shmem
( 1287)       !! * SHMEM version
( 1288) 
( 1289)       imigr = jprecj * jpi
( 1290) 
( 1291)       SELECT CASE ( nbondj )
( 1292)       CASE ( -1 )
( 1293)          CALL shmem_put( t2sn(1,1,2), t2sn(1,1,1), imigr, nono )
( 1294)       CASE ( 0 )
( 1295)          CALL shmem_put( t2ns(1,1,2), t2ns(1,1,1), imigr, noso )
( 1296)          CALL shmem_put( t2sn(1,1,2), t2sn(1,1,1), imigr, nono )
( 1297)       CASE ( 1 )
( 1298)          CALL shmem_put( t2ns(1,1,2), t2ns(1,1,1), imigr, noso )
( 1299)       END SELECT 
( 1300)       CALL barrier()
( 1301)       CALL shmem_udcflush()
( 1302) 
( 1303) #elif defined key_mpp_mpi
( 1304)       !! * MPI version
( 1305) 
( 1306)       imigr = jprecj * jpi
( 1307) 
( 1308)       SELECT CASE ( nbondj )
( 1309)       CASE ( -1 )
( 1310)          CALL mppsend( 4, t2sn(1,1,1), imigr, nono, ml_req1 )
( 1311)          CALL mpprecv( 3, t2ns(1,1,2), imigr )
( 1312)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1313)       CASE ( 0 )
( 1314)          CALL mppsend( 3, t2ns(1,1,1), imigr, noso, ml_req1 )
( 1315)          CALL mppsend( 4, t2sn(1,1,1), imigr, nono, ml_req2 )
( 1316)          CALL mpprecv( 3, t2ns(1,1,2), imigr )
( 1317)          CALL mpprecv( 4, t2sn(1,1,2), imigr )
( 1318)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1319)          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 1320)       CASE ( 1 )
( 1321)          CALL mppsend( 3, t2ns(1,1,1), imigr, noso, ml_req1 )
( 1322)          CALL mpprecv( 4, t2sn(1,1,2), imigr )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 24

( 1323)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1324)       END SELECT
( 1325)   
( 1326) #endif
( 1327) 
( 1328)       ! 3.3 Write Dirichlet lateral conditions
( 1329) 
( 1330)       ijhom = nlcj - jprecj
( 1331) 
( 1332)       SELECT CASE ( nbondj )
( 1333)       CASE ( -1 )
( 1334)          DO jl = 1, jprecj
( 1335)             pt2d(:,ijhom+jl) = t2ns(:,jl,2)
( 1336)          END DO
( 1337)       CASE ( 0 )
( 1338)          DO jl = 1, jprecj
( 1339)             pt2d(:,jl      ) = t2sn(:,jl,2)
( 1340)             pt2d(:,ijhom+jl) = t2ns(:,jl,2)
( 1341)          END DO
( 1342)       CASE ( 1 ) 
( 1343)          DO jl = 1, jprecj
( 1344)             pt2d(:,jl      ) = t2sn(:,jl,2)
( 1345)          END DO
( 1346)       END SELECT 
( 1347)   
( 1348) 
( 1349)       ! 4. north fold treatment
( 1350)       ! -----------------------
( 1351)   
( 1352)       IF (PRESENT(cd_mpp)) THEN
( 1353)          ! No north fold treatment (it is assumed to be already OK)
( 1354)      
( 1355)       ELSE      
( 1356) 
( 1357)       ! 4.1 treatment without exchange (jpni odd)
( 1358)       
( 1359)       SELECT CASE ( jpni )
( 1360)   
( 1361)       CASE ( 1 ) ! only one proc along I, no mpp exchange
( 1362)   
( 1363)          SELECT CASE ( npolj )
( 1364)   
( 1365)          CASE ( 3 , 4 )   !  T pivot
( 1366)             iloc = jpiglo - 2 * ( nimpp - 1 )
( 1367)   
( 1368)             SELECT CASE ( cd_type )
( 1369)   
( 1370)             CASE ( 'T' , 'S', 'W' )
( 1371)                DO ji = 2, nlci
( 1372)                   ijt=iloc-ji+2
( 1373)                   pt2d(ji,nlcj) = psgn * pt2d(ijt,nlcj-2)
( 1374)                END DO
( 1375)                DO ji = nlci/2+1, nlci
( 1376)                   ijt=iloc-ji+2
( 1377)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
( 1378)                END DO
( 1379)   
( 1380)             CASE ( 'U' )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 25

( 1381)                DO ji = 1, nlci-1
( 1382)                   iju=iloc-ji+1
( 1383)                   pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-2)
( 1384)                END DO
( 1385)                DO ji = nlci/2, nlci-1
( 1386)                   iju=iloc-ji+1
( 1387)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
( 1388)                END DO
( 1389)   
( 1390)             CASE ( 'V' )
( 1391)                DO ji = 2, nlci
( 1392)                   ijt=iloc-ji+2
( 1393)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-2)
( 1394)                   pt2d(ji,nlcj  ) = psgn * pt2d(ijt,nlcj-3)
( 1395)                END DO
( 1396)   
( 1397)             CASE ( 'F', 'G' )
( 1398)                DO ji = 1, nlci-1
( 1399)                   iju=iloc-ji+1
( 1400)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-2)
( 1401)                   pt2d(ji,nlcj  ) = psgn * pt2d(iju,nlcj-3)
( 1402)                END DO
( 1403)   
( 1404)             CASE ( 'I' )                                  ! ice U-V point
( 1405)                pt2d(2,nlcj) = psgn * pt2d(3,nlcj-1)
( 1406)                DO ji = 3, nlci
( 1407)                   iju = iloc - ji + 3
( 1408)                   pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-1)
( 1409)                END DO
( 1410)   
( 1411)             END SELECT
( 1412)   
( 1413)          CASE ( 5 , 6 )                 ! F pivot
( 1414)             iloc=jpiglo-2*(nimpp-1)
( 1415)   
( 1416)             SELECT CASE (cd_type )
( 1417)   
( 1418)             CASE ( 'T', 'S', 'W' )
( 1419)                DO ji = 1, nlci
( 1420)                   ijt=iloc-ji+1
( 1421)                   pt2d(ji,nlcj) = psgn * pt2d(ijt,nlcj-1)
( 1422)                END DO
( 1423)   
( 1424)             CASE ( 'U' )
( 1425)                DO ji = 1, nlci-1
( 1426)                   iju=iloc-ji
( 1427)                   pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-1)
( 1428)                END DO
( 1429) 
( 1430)             CASE ( 'V' )
( 1431)                DO ji = 1, nlci
( 1432)                   ijt=iloc-ji+1
( 1433)                   pt2d(ji,nlcj  ) = psgn * pt2d(ijt,nlcj-2)
( 1434)                END DO
( 1435)                DO ji = nlci/2+1, nlci
( 1436)                   ijt=iloc-ji+1
( 1437)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
( 1438)                END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 26

( 1439)   
( 1440)             CASE ( 'F', 'G' )
( 1441)                DO ji = 1, nlci-1
( 1442)                   iju=iloc-ji
( 1443)                   pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-2)
( 1444)                END DO
( 1445)                DO ji = nlci/2+1, nlci-1
( 1446)                   iju=iloc-ji
( 1447)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
( 1448)                END DO
( 1449)   
( 1450)             CASE ( 'I' )                                  ! ice U-V point
( 1451)                pt2d( 2 ,nlcj) = 0.e0
( 1452)                DO ji = 2 , nlci-1
( 1453)                   ijt = iloc - ji + 2
( 1454)                   pt2d(ji,nlcj)= 0.5 * ( pt2d(ji,nlcj-1) + psgn * pt2d(ijt,nlcj-1) )
( 1455)                END DO
( 1456)   
( 1457)             END SELECT   ! cd_type
( 1458)   
( 1459)          END SELECT   ! npolj
( 1460) 
( 1461)       CASE DEFAULT   ! more than 1 proc along I
( 1462)          IF( npolj /= 0 )   CALL mpp_lbc_north( pt2d, cd_type, psgn )   ! only for northern procs.
( 1463) 
( 1464)       END SELECT   ! jpni
( 1465) 
( 1466)       ENDIF
( 1467) 
( 1468)       ! 5. East and west directions
( 1469)       ! ---------------------------
( 1470) 
( 1471)       SELECT CASE ( npolj )
( 1472) 
( 1473)       CASE ( 3, 4, 5, 6 )
( 1474) 
( 1475)          ! 5.1 Read Dirichlet lateral conditions
( 1476) 
( 1477)          SELECT CASE ( nbondi )
( 1478)          CASE ( -1, 0, 1 )
( 1479)             iihom = nlci-nreci
( 1480)             DO jl = 1, jpreci
( 1481)                DO jj = 1, jpj
( 1482)                   t2ew(jj,jl,1) = pt2d(jpreci+jl,jj)
( 1483)                   t2we(jj,jl,1) = pt2d(iihom +jl,jj)
( 1484)                END DO
( 1485)             END DO
( 1486)          END SELECT
( 1487) 
( 1488)          ! 5.2 Migrations
( 1489) 
( 1490) #if defined key_mpp_shmem
( 1491)          !! * SHMEM version
( 1492) 
( 1493)          imigr=jpreci*jpj
( 1494) 
( 1495)          SELECT CASE ( nbondi )
( 1496)          CASE ( -1 )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 27

( 1497)             CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr, noea )
( 1498)          CASE ( 0 )
( 1499)             CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr, nowe )
( 1500)             CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr, noea )
( 1501)          CASE ( 1 )
( 1502)             CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr, nowe )
( 1503)          END SELECT
( 1504) 
( 1505)          CALL barrier()
( 1506)          CALL shmem_udcflush()
( 1507)   
( 1508) #elif defined key_mpp_mpi
( 1509)          !! * MPI version
( 1510)   
( 1511)          imigr=jpreci*jpj
( 1512)   
( 1513)          SELECT CASE ( nbondi )
( 1514)          CASE ( -1 )
( 1515)             CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req1 )
( 1516)             CALL mpprecv( 1, t2ew(1,1,2), imigr )
( 1517)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1518)          CASE ( 0 )
( 1519)             CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
( 1520)             CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req2 )
( 1521)             CALL mpprecv( 1, t2ew(1,1,2), imigr )
( 1522)             CALL mpprecv( 2, t2we(1,1,2), imigr )
( 1523)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1524)             IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 1525)          CASE ( 1 )
( 1526)             CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
( 1527)             CALL mpprecv( 2, t2we(1,1,2), imigr )
( 1528)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1529)          END SELECT 
( 1530) #endif
( 1531) 
( 1532)          ! 5.3 Write Dirichlet lateral conditions
( 1533)   
( 1534)          iihom = nlci - jpreci
( 1535)   
( 1536)          SELECT CASE ( nbondi )
( 1537)          CASE ( -1 )
( 1538)             DO jl = 1, jpreci
( 1539)                pt2d(iihom+jl,:) = t2ew(:,jl,2)
( 1540)             END DO
( 1541)          CASE ( 0 )
( 1542)             DO jl = 1, jpreci
( 1543)                pt2d(jl      ,:) = t2we(:,jl,2)
( 1544)                pt2d(iihom+jl,:) = t2ew(:,jl,2)
( 1545)             END DO
( 1546)          CASE ( 1 )
( 1547)             DO jl = 1, jpreci
( 1548)                pt2d(jl,:) = t2we(:,jl,2)
( 1549)             END DO
( 1550)          END SELECT 
( 1551)   
( 1552)       END SELECT   ! npolj
( 1553)   
( 1554)    END SUBROUTINE mpp_lnk_2d






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 28

( 1555) 
( 1556) 
( 1557)    SUBROUTINE mpp_lnk_3d_gather( ptab1, cd_type1, ptab2, cd_type2, psgn )
( 1558)       !!----------------------------------------------------------------------
( 1559)       !!                  ***  routine mpp_lnk_3d_gather  ***
( 1560)       !!
( 1561)       !! ** Purpose :   Message passing manadgement for two 3D arrays
( 1562)       !!
( 1563)       !! ** Method  :   Use mppsend and mpprecv function for passing mask 
( 1564)       !!      between processors following neighboring subdomains.
( 1565)       !!            domain parameters
( 1566)       !!                    nlci   : first dimension of the local subdomain
( 1567)       !!                    nlcj   : second dimension of the local subdomain
( 1568)       !!                    nbondi : mark for "east-west local boundary"
( 1569)       !!                    nbondj : mark for "north-south local boundary"
( 1570)       !!                    noea   : number for local neighboring processors 
( 1571)       !!                    nowe   : number for local neighboring processors
( 1572)       !!                    noso   : number for local neighboring processors
( 1573)       !!                    nono   : number for local neighboring processors
( 1574)       !!
( 1575)       !! ** Action  :   ptab1 and ptab2  with update value at its periphery
( 1576)       !!
( 1577)       !!----------------------------------------------------------------------
( 1578)       !! * Arguments
( 1579)       CHARACTER(len=1) , INTENT( in ) ::   &
( 1580)          cd_type1, cd_type2       ! define the nature of ptab array grid-points
( 1581)          !                        ! = T , U , V , F , W points
( 1582)          !                        ! = S : T-point, north fold treatment ???
( 1583)          !                        ! = G : F-point, north fold treatment ???
( 1584)       REAL(wp), INTENT( in ) ::   &
( 1585)          psgn          ! control of the sign change
( 1586)          !             !   = -1. , the sign is changed if north fold boundary
( 1587)          !             !   =  1. , the sign is kept  if north fold boundary
( 1588)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
( 1589)          ptab1, ptab2             ! 3D array on which the boundary condition is applied
( 1590) 
( 1591)       !! * Local variables
( 1592)       INTEGER ::   ji, jk, jl   ! dummy loop indices
( 1593)       INTEGER ::   imigr, iihom, ijhom, iloc, ijt, iju   ! temporary integers
( 1594)       INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
( 1595)       INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
( 1596)       !!----------------------------------------------------------------------
( 1597) 
( 1598)       ! 1. standard boundary treatment
( 1599)       ! ------------------------------
( 1600)       !                                        ! East-West boundaries
( 1601)       !                                        ! ====================
( 1602)       IF( nbondi == 2 .AND.   &      ! Cyclic east-west
( 1603)          &   (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
( 1604)          ptab1( 1 ,:,:) = ptab1(jpim1,:,:)
( 1605)          ptab1(jpi,:,:) = ptab1(  2  ,:,:)
( 1606)          ptab2( 1 ,:,:) = ptab2(jpim1,:,:)
( 1607)          ptab2(jpi,:,:) = ptab2(  2  ,:,:)
( 1608) 
( 1609)       ELSE                           ! closed
( 1610)          SELECT CASE ( cd_type1 )
( 1611)          CASE ( 'T', 'U', 'V', 'W' )
( 1612)             ptab1(     1       :jpreci,:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 29

( 1613)             ptab1(nlci-jpreci+1:jpi   ,:,:) = 0.e0
( 1614)          CASE ( 'F' )
( 1615)             ptab1(nlci-jpreci+1:jpi   ,:,:) = 0.e0
( 1616)          END SELECT 
( 1617)          SELECT CASE ( cd_type2 )
( 1618)          CASE ( 'T', 'U', 'V', 'W' )
( 1619)             ptab2(     1       :jpreci,:,:) = 0.e0
( 1620)             ptab2(nlci-jpreci+1:jpi   ,:,:) = 0.e0
( 1621)          CASE ( 'F' )
( 1622)             ptab2(nlci-jpreci+1:jpi   ,:,:) = 0.e0
( 1623)          END SELECT 
( 1624)       ENDIF
( 1625) 
( 1626)       !                                        ! North-South boundaries
( 1627)       !                                        ! ======================
( 1628)       SELECT CASE ( cd_type1 )
( 1629)       CASE ( 'T', 'U', 'V', 'W' )
( 1630)          ptab1(:,     1       :jprecj,:) = 0.e0
( 1631)          ptab1(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
( 1632)       CASE ( 'F' )
( 1633)          ptab1(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
( 1634)       END SELECT
( 1635) 
( 1636)       SELECT CASE ( cd_type2 )
( 1637)       CASE ( 'T', 'U', 'V', 'W' )
( 1638)          ptab2(:,     1       :jprecj,:) = 0.e0
( 1639)          ptab2(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
( 1640)       CASE ( 'F' )
( 1641)          ptab2(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
( 1642)       END SELECT
( 1643) 
( 1644) 
( 1645)       ! 2. East and west directions exchange
( 1646)       ! ------------------------------------
( 1647) 
( 1648)       ! 2.1 Read Dirichlet lateral conditions
( 1649) 
( 1650)       SELECT CASE ( nbondi )
( 1651)       CASE ( -1, 0, 1 )    ! all exept 2 
( 1652)          iihom = nlci-nreci
( 1653)          DO jl = 1, jpreci
( 1654)             t4ew(:,jl,:,1,1) = ptab1(jpreci+jl,:,:)
( 1655)             t4we(:,jl,:,1,1) = ptab1(iihom +jl,:,:)
( 1656)             t4ew(:,jl,:,2,1) = ptab2(jpreci+jl,:,:)
( 1657)             t4we(:,jl,:,2,1) = ptab2(iihom +jl,:,:)
( 1658)          END DO
( 1659)       END SELECT
( 1660) 
( 1661)       ! 2.2 Migrations
( 1662) 
( 1663) #if defined key_mpp_shmem
( 1664)       !! * SHMEM version
( 1665) 
( 1666)       imigr = jpreci * jpj * jpk *2
( 1667) 
( 1668)       SELECT CASE ( nbondi )
( 1669)       CASE ( -1 )
( 1670)          CALL shmem_put( t4we(1,1,1,1,2), t4we(1,1,1,1,1), imigr, noea )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 30

( 1671)       CASE ( 0 )
( 1672)          CALL shmem_put( t4ew(1,1,1,1,2), t4ew(1,1,1,1,1), imigr, nowe )
( 1673)          CALL shmem_put( t4we(1,1,1,1,2), t4we(1,1,1,1,1), imigr, noea )
( 1674)       CASE ( 1 )
( 1675)          CALL shmem_put( t4ew(1,1,1,1,2), t4ew(1,1,1,1,1), imigr, nowe )
( 1676)       END SELECT
( 1677) 
( 1678)       CALL barrier()
( 1679)       CALL shmem_udcflush()
( 1680) 
( 1681) #elif defined key_mpp_mpi
( 1682)       !! * Local variables   (MPI version)
( 1683) 
( 1684)       imigr = jpreci * jpj * jpk *2
( 1685) 
( 1686)       SELECT CASE ( nbondi ) 
( 1687)       CASE ( -1 )
( 1688)          CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req1 )
( 1689)          CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
( 1690)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1691)       CASE ( 0 )
( 1692)          CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
( 1693)          CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req2 )
( 1694)          CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
( 1695)          CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
( 1696)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1697)          IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
( 1698)       CASE ( 1 )
( 1699)          CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
( 1700)          CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
( 1701)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1702)       END SELECT
( 1703) #endif
( 1704) 
( 1705)       ! 2.3 Write Dirichlet lateral conditions
( 1706) 
( 1707)       iihom = nlci-jpreci
( 1708) 
( 1709)       SELECT CASE ( nbondi )
( 1710)       CASE ( -1 )
( 1711)          DO jl = 1, jpreci
( 1712)             ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
( 1713)             ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
( 1714)          END DO
( 1715)       CASE ( 0 ) 
( 1716)          DO jl = 1, jpreci
( 1717)             ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
( 1718)             ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
( 1719)             ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
( 1720)             ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
( 1721)          END DO
( 1722)       CASE ( 1 )
( 1723)          DO jl = 1, jpreci
( 1724)             ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
( 1725)             ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
( 1726)          END DO
( 1727)       END SELECT
( 1728) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 31

( 1729) 
( 1730)       ! 3. North and south directions
( 1731)       ! -----------------------------
( 1732) 
( 1733)       ! 3.1 Read Dirichlet lateral conditions
( 1734) 
( 1735)       IF( nbondj /= 2 ) THEN
( 1736)          ijhom = nlcj-nrecj
( 1737)          DO jl = 1, jprecj
( 1738)             t4sn(:,jl,:,1,1) = ptab1(:,ijhom +jl,:)
( 1739)             t4ns(:,jl,:,1,1) = ptab1(:,jprecj+jl,:)
( 1740)             t4sn(:,jl,:,2,1) = ptab2(:,ijhom +jl,:)
( 1741)             t4ns(:,jl,:,2,1) = ptab2(:,jprecj+jl,:)
( 1742)          END DO
( 1743)       ENDIF
( 1744) 
( 1745)       ! 3.2 Migrations
( 1746) 
( 1747) #if defined key_mpp_shmem
( 1748)       !! * SHMEM version
( 1749) 
( 1750)       imigr = jprecj * jpi * jpk * 2
( 1751) 
( 1752)       SELECT CASE ( nbondj )
( 1753)       CASE ( -1 )
( 1754)          CALL shmem_put( t4sn(1,1,1,1,2), t4sn(1,1,1,1,1), imigr, nono )
( 1755)       CASE ( 0 )
( 1756)          CALL shmem_put( t4ns(1,1,1,1,2), t4ns(1,1,1,1,1), imigr, noso )
( 1757)          CALL shmem_put( t4sn(1,1,1,1,2), t4sn(1,1,1,1,1), imigr, nono )
( 1758)       CASE ( 1 )
( 1759)          CALL shmem_put( t4ns(1,1,1,1,2), t4ns(1,1,1,1;,1), imigr, noso )
( 1760)       END SELECT
( 1761) 
( 1762)       CALL barrier()
( 1763)       CALL shmem_udcflush()
( 1764) 
( 1765) #elif defined key_mpp_mpi
( 1766)       !! * Local variables   (MPI version)
( 1767)   
( 1768)       imigr=jprecj * jpi * jpk * 2
( 1769) 
( 1770)       SELECT CASE ( nbondj )     
( 1771)       CASE ( -1 )
( 1772)          CALL mppsend( 4, t4sn(1,1,1,1,1), imigr, nono, ml_req1 )
( 1773)          CALL mpprecv( 3, t4ns(1,1,1,1,2), imigr )
( 1774)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1775)       CASE ( 0 )
( 1776)          CALL mppsend( 3, t4ns(1,1,1,1,1), imigr, noso, ml_req1 )
( 1777)          CALL mppsend( 4, t4sn(1,1,1,1,1), imigr, nono, ml_req2 )
( 1778)          CALL mpprecv( 3, t4ns(1,1,1,1,2), imigr )
( 1779)          CALL mpprecv( 4, t4sn(1,1,1,1,2), imigr )
( 1780)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1781)          IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
( 1782)       CASE ( 1 ) 
( 1783)          CALL mppsend( 3, t4ns(1,1,1,1,1), imigr, noso, ml_req1 )
( 1784)          CALL mpprecv( 4, t4sn(1,1,1,1,2), imigr )
( 1785)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1786)       END SELECT






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 32

( 1787) 
( 1788) #endif
( 1789) 
( 1790)       ! 3.3 Write Dirichlet lateral conditions
( 1791) 
( 1792)       ijhom = nlcj-jprecj
( 1793) 
( 1794)       SELECT CASE ( nbondj )
( 1795)       CASE ( -1 )
( 1796)          DO jl = 1, jprecj
( 1797)             ptab1(:,ijhom+jl,:) = t4ns(:,jl,:,1,2)
( 1798)             ptab2(:,ijhom+jl,:) = t4ns(:,jl,:,2,2)
( 1799)          END DO
( 1800)       CASE ( 0 ) 
( 1801)          DO jl = 1, jprecj
( 1802)             ptab1(:,jl      ,:) = t4sn(:,jl,:,1,2)
( 1803)             ptab1(:,ijhom+jl,:) = t4ns(:,jl,:,1,2)
( 1804)             ptab2(:,jl      ,:) = t4sn(:,jl,:,2,2)
( 1805)             ptab2(:,ijhom+jl,:) = t4ns(:,jl,:,2,2)
( 1806)          END DO
( 1807)       CASE ( 1 )
( 1808)          DO jl = 1, jprecj
( 1809)             ptab1(:,jl,:) = t4sn(:,jl,:,1,2)
( 1810)             ptab2(:,jl,:) = t4sn(:,jl,:,2,2)
( 1811)          END DO
( 1812)       END SELECT
( 1813) 
( 1814) 
( 1815)       ! 4. north fold treatment
( 1816)       ! -----------------------
( 1817) 
( 1818)       ! 4.1 treatment without exchange (jpni odd)
( 1819)       !     T-point pivot  
( 1820) 
( 1821)       SELECT CASE ( jpni )
( 1822) 
( 1823)       CASE ( 1 )  ! only one proc along I, no mpp exchange
( 1824) 
( 1825)       SELECT CASE ( npolj )
( 1826)   
( 1827)          CASE ( 3 , 4 )    ! T pivot
( 1828)             iloc = jpiglo - 2 * ( nimpp - 1 )
( 1829) 
( 1830)             SELECT CASE ( cd_type1 )
( 1831) 
( 1832)             CASE ( 'T' , 'S', 'W' )
( 1833)                DO jk = 1, jpk
( 1834)                   DO ji = 2, nlci
( 1835)                      ijt=iloc-ji+2
( 1836)                      ptab1(ji,nlcj,jk) = psgn * ptab1(ijt,nlcj-2,jk)
( 1837)                   END DO
( 1838)                   DO ji = nlci/2+1, nlci
( 1839)                      ijt=iloc-ji+2
( 1840)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-1,jk)
( 1841)                   END DO
( 1842)                END DO
( 1843)           
( 1844)             CASE ( 'U' )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 33

( 1845)                DO jk = 1, jpk
( 1846)                   DO ji = 1, nlci-1
( 1847)                      iju=iloc-ji+1
( 1848)                      ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-2,jk)
( 1849)                   END DO
( 1850)                   DO ji = nlci/2, nlci-1
( 1851)                      iju=iloc-ji+1
( 1852)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-1,jk)
( 1853)                   END DO
( 1854)                END DO
( 1855) 
( 1856)             CASE ( 'V' )
( 1857)                DO jk = 1, jpk
( 1858)                   DO ji = 2, nlci
( 1859)                      ijt=iloc-ji+2
( 1860)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-2,jk)
( 1861)                      ptab1(ji,nlcj  ,jk) = psgn * ptab1(ijt,nlcj-3,jk)
( 1862)                   END DO
( 1863)                END DO
( 1864) 
( 1865)             CASE ( 'F', 'G' )
( 1866)                DO jk = 1, jpk
( 1867)                   DO ji = 1, nlci-1
( 1868)                      iju=iloc-ji+1
( 1869)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-2,jk)
( 1870)                      ptab1(ji,nlcj  ,jk) = psgn * ptab1(iju,nlcj-3,jk)
( 1871)                   END DO
( 1872)                END DO
( 1873)   
( 1874)             END SELECT
( 1875)             
( 1876)             SELECT CASE ( cd_type2 )
( 1877) 
( 1878)             CASE ( 'T' , 'S', 'W' )
( 1879)                DO jk = 1, jpk
( 1880)                   DO ji = 2, nlci
( 1881)                      ijt=iloc-ji+2
( 1882)                      ptab2(ji,nlcj,jk) = psgn * ptab2(ijt,nlcj-2,jk)
( 1883)                   END DO
( 1884)                   DO ji = nlci/2+1, nlci
( 1885)                      ijt=iloc-ji+2
( 1886)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-1,jk)
( 1887)                   END DO
( 1888)                END DO
( 1889)           
( 1890)             CASE ( 'U' )
( 1891)                DO jk = 1, jpk
( 1892)                   DO ji = 1, nlci-1
( 1893)                      iju=iloc-ji+1
( 1894)                      ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-2,jk)
( 1895)                   END DO
( 1896)                   DO ji = nlci/2, nlci-1
( 1897)                      iju=iloc-ji+1
( 1898)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-1,jk)
( 1899)                   END DO
( 1900)                END DO
( 1901) 
( 1902)             CASE ( 'V' )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 34

( 1903)                DO jk = 1, jpk
( 1904)                   DO ji = 2, nlci
( 1905)                      ijt=iloc-ji+2
( 1906)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-2,jk)
( 1907)                      ptab2(ji,nlcj  ,jk) = psgn * ptab2(ijt,nlcj-3,jk)
( 1908)                   END DO
( 1909)                END DO
( 1910) 
( 1911)             CASE ( 'F', 'G' )
( 1912)                DO jk = 1, jpk
( 1913)                   DO ji = 1, nlci-1
( 1914)                      iju=iloc-ji+1
( 1915)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-2,jk)
( 1916)                      ptab2(ji,nlcj  ,jk) = psgn * ptab2(iju,nlcj-3,jk)
( 1917)                   END DO
( 1918)                END DO
( 1919)   
( 1920)           END SELECT
( 1921)        
( 1922)          CASE ( 5 , 6 ) ! F pivot
( 1923)             iloc=jpiglo-2*(nimpp-1)
( 1924)   
( 1925)             SELECT CASE ( cd_type1 )
( 1926) 
( 1927)             CASE ( 'T' , 'S', 'W' )
( 1928)                DO jk = 1, jpk
( 1929)                   DO ji = 1, nlci
( 1930)                      ijt=iloc-ji+1
( 1931)                      ptab1(ji,nlcj,jk) = psgn * ptab1(ijt,nlcj-1,jk)
( 1932)                   END DO
( 1933)                END DO
( 1934) 
( 1935)             CASE ( 'U' )
( 1936)                DO jk = 1, jpk
( 1937)                   DO ji = 1, nlci-1
( 1938)                      iju=iloc-ji
( 1939)                      ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-1,jk)
( 1940)                   END DO
( 1941)                END DO
( 1942) 
( 1943)             CASE ( 'V' )
( 1944)                DO jk = 1, jpk
( 1945)                   DO ji = 1, nlci
( 1946)                      ijt=iloc-ji+1
( 1947)                      ptab1(ji,nlcj  ,jk) = psgn * ptab1(ijt,nlcj-2,jk)
( 1948)                   END DO
( 1949)                   DO ji = nlci/2+1, nlci
( 1950)                      ijt=iloc-ji+1
( 1951)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-1,jk)
( 1952)                   END DO
( 1953)                END DO
( 1954) 
( 1955)             CASE ( 'F', 'G' )
( 1956)                DO jk = 1, jpk
( 1957)                   DO ji = 1, nlci-1
( 1958)                      iju=iloc-ji
( 1959)                      ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-2,jk)
( 1960)                   END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 35

( 1961)                   DO ji = nlci/2+1, nlci-1
( 1962)                      iju=iloc-ji
( 1963)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-1,jk)
( 1964)                   END DO
( 1965)                END DO
( 1966)             END SELECT  ! cd_type1
( 1967) 
( 1968)             SELECT CASE ( cd_type2 )
( 1969) 
( 1970)             CASE ( 'T' , 'S', 'W' )
( 1971)                DO jk = 1, jpk
( 1972)                   DO ji = 1, nlci
( 1973)                      ijt=iloc-ji+1
( 1974)                      ptab2(ji,nlcj,jk) = psgn * ptab2(ijt,nlcj-1,jk)
( 1975)                   END DO
( 1976)                END DO
( 1977) 
( 1978)             CASE ( 'U' )
( 1979)                DO jk = 1, jpk
( 1980)                   DO ji = 1, nlci-1
( 1981)                      iju=iloc-ji
( 1982)                      ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-1,jk)
( 1983)                   END DO
( 1984)                END DO
( 1985) 
( 1986)             CASE ( 'V' )
( 1987)                DO jk = 1, jpk
( 1988)                   DO ji = 1, nlci
( 1989)                      ijt=iloc-ji+1
( 1990)                      ptab2(ji,nlcj  ,jk) = psgn * ptab2(ijt,nlcj-2,jk)
( 1991)                   END DO
( 1992)                   DO ji = nlci/2+1, nlci
( 1993)                      ijt=iloc-ji+1
( 1994)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-1,jk)
( 1995)                   END DO
( 1996)                END DO
( 1997) 
( 1998)             CASE ( 'F', 'G' )
( 1999)                DO jk = 1, jpk
( 2000)                   DO ji = 1, nlci-1
( 2001)                      iju=iloc-ji
( 2002)                      ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-2,jk)
( 2003)                   END DO
( 2004)                   DO ji = nlci/2+1, nlci-1
( 2005)                      iju=iloc-ji
( 2006)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-1,jk)
( 2007)                   END DO
( 2008)                END DO
( 2009) 
( 2010)             END SELECT  ! cd_type2
( 2011) 
( 2012)          END SELECT     !  npolj
( 2013)   
( 2014)       CASE DEFAULT ! more than 1 proc along I
( 2015)          IF ( npolj /= 0 ) THEN
( 2016)             CALL mpp_lbc_north (ptab1, cd_type1, psgn)  ! only for northern procs.
( 2017)             CALL mpp_lbc_north (ptab2, cd_type2, psgn)  ! only for northern procs.
( 2018)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 36

( 2019) 
( 2020)       END SELECT ! jpni 
( 2021) 
( 2022) 
( 2023)       ! 5. East and west directions exchange
( 2024)       ! ------------------------------------
( 2025) 
( 2026)       SELECT CASE ( npolj )
( 2027) 
( 2028)       CASE ( 3, 4, 5, 6 )
( 2029) 
( 2030)          ! 5.1 Read Dirichlet lateral conditions
( 2031) 
( 2032)          SELECT CASE ( nbondi )
( 2033) 
( 2034)          CASE ( -1, 0, 1 )
( 2035)             iihom = nlci-nreci
( 2036)             DO jl = 1, jpreci
( 2037)                t4ew(:,jl,:,1,1) = ptab1(jpreci+jl,:,:)
( 2038)                t4we(:,jl,:,1,1) = ptab1(iihom +jl,:,:)
( 2039)                t4ew(:,jl,:,2,1) = ptab2(jpreci+jl,:,:)
( 2040)                t4we(:,jl,:,2,1) = ptab2(iihom +jl,:,:)
( 2041)             END DO
( 2042) 
( 2043)          END SELECT
( 2044) 
( 2045)          ! 5.2 Migrations
( 2046) 
( 2047) #if defined key_mpp_shmem
( 2048)          !! SHMEM version
( 2049) 
( 2050)          imigr = jpreci * jpj * jpk * 2
( 2051) 
( 2052)          SELECT CASE ( nbondi )
( 2053)          CASE ( -1 )
( 2054)             CALL shmem_put( t4we(1,1,1,1,2), t4we(1,1,1,1,1), imigr, noea )
( 2055)          CASE ( 0 )
( 2056)             CALL shmem_put( t4ew(1,1,1,1,2), t4ew(1,1,1,1,1), imigr, nowe )
( 2057)             CALL shmem_put( t4we(1,1,1,1,2), t4we(1,1,1,1,1), imigr, noea )
( 2058)          CASE ( 1 )
( 2059)             CALL shmem_put( t4ew(1,1,1,1,2), t4ew(1,1,1,1,1), imigr, nowe )
( 2060)          END SELECT
( 2061) 
( 2062)          CALL barrier()
( 2063)          CALL shmem_udcflush()
( 2064) 
( 2065) #elif defined key_mpp_mpi
( 2066)          !! MPI version
( 2067) 
( 2068)          imigr = jpreci * jpj * jpk * 2
( 2069)   
( 2070)          SELECT CASE ( nbondi )
( 2071)          CASE ( -1 )
( 2072)             CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req1 )
( 2073)             CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
( 2074)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2075)          CASE ( 0 )
( 2076)             CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 37

( 2077)             CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req2 )
( 2078)             CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
( 2079)             CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
( 2080)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2081)             IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2082)          CASE ( 1 )
( 2083)             CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
( 2084)             CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
( 2085)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2086)          END SELECT
( 2087) #endif
( 2088) 
( 2089)          ! 5.3 Write Dirichlet lateral conditions
( 2090) 
( 2091)          iihom = nlci-jpreci
( 2092) 
( 2093)          SELECT CASE ( nbondi)
( 2094)          CASE ( -1 )
( 2095)             DO jl = 1, jpreci
( 2096)                ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
( 2097)                ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
( 2098)             END DO
( 2099)          CASE ( 0 ) 
( 2100)             DO jl = 1, jpreci
( 2101)                ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
( 2102)                ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
( 2103)                ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
( 2104)                ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
( 2105)             END DO
( 2106)          CASE ( 1 )
( 2107)             DO jl = 1, jpreci
( 2108)                ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
( 2109)                ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
( 2110)             END DO
( 2111)          END SELECT
( 2112) 
( 2113)       END SELECT    ! npolj 
( 2114) 
( 2115)    END SUBROUTINE mpp_lnk_3d_gather
( 2116) 
( 2117) 
( 2118)    SUBROUTINE mpp_lnk_2d_e( pt2d, cd_type, psgn )
( 2119)       !!----------------------------------------------------------------------
( 2120)       !!                  ***  routine mpp_lnk_2d_e  ***
( 2121)       !!                  
( 2122)       !! ** Purpose :   Message passing manadgement for 2d array (with halo)
( 2123)       !!
( 2124)       !! ** Method  :   Use mppsend and mpprecv function for passing mask 
( 2125)       !!      between processors following neighboring subdomains.
( 2126)       !!            domain parameters
( 2127)       !!                    nlci   : first dimension of the local subdomain
( 2128)       !!                    nlcj   : second dimension of the local subdomain
( 2129)       !!                    jpr2di : number of rows for extra outer halo
( 2130)       !!                    jpr2dj : number of columns for extra outer halo
( 2131)       !!                    nbondi : mark for "east-west local boundary"
( 2132)       !!                    nbondj : mark for "north-south local boundary"
( 2133)       !!                    noea   : number for local neighboring processors 
( 2134)       !!                    nowe   : number for local neighboring processors






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 38

( 2135)       !!                    noso   : number for local neighboring processors
( 2136)       !!                    nono   : number for local neighboring processors
( 2137)       !!   
( 2138)       !! History :
( 2139)       !!       
( 2140)       !!   9.0  !  05-09  (R. Benshila, G. Madec)  original code
( 2141)       !!
( 2142)       !!----------------------------------------------------------------------
( 2143)       !! * Arguments
( 2144)       CHARACTER(len=1) , INTENT( in ) ::   &
( 2145)          cd_type       ! define the nature of pt2d array grid-points
( 2146)          !             !  = T , U , V , F , W 
( 2147)          !             !  = S : T-point, north fold treatment
( 2148)          !             !  = G : F-point, north fold treatment
( 2149)          !             !  = I : sea-ice velocity at F-point with index shift
( 2150)       REAL(wp), INTENT( in ) ::   &
( 2151)          psgn          ! control of the sign change
( 2152)          !             !   = -1. , the sign is changed if north fold boundary
( 2153)          !             !   =  1. , the sign is kept  if north fold boundary
( 2154)       REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj), INTENT( inout ) ::   &
( 2155)          pt2d          ! 2D array on which the boundary condition is applied
( 2156) 
( 2157)       !! * Local variables
( 2158)       INTEGER  ::   ji, jl      ! dummy loop indices
( 2159)       INTEGER  ::   &
( 2160)          imigr, iihom, ijhom,    &  ! temporary integers
( 2161)          iloc, ijt, iju             !    "          "
( 2162)       INTEGER  ::   &
( 2163)          ipreci, iprecj             ! temporary integers
( 2164)       INTEGER  ::   ml_req1, ml_req2, ml_err     ! for isend
( 2165)       INTEGER  ::   ml_stat(MPI_STATUS_SIZE)     ! for isend
( 2166)      !!---------------------------------------------------------------------
( 2167) 
( 2168)       ! take into account outer extra 2D overlap area
( 2169)       ipreci = jpreci + jpr2di
( 2170)       iprecj = jprecj + jpr2dj
( 2171) 
( 2172) 
( 2173)       ! 1. standard boundary treatment
( 2174)       ! ------------------------------
( 2175) 
( 2176)       !                                        ! East-West boundaries
( 2177)       !                                        ! ====================
( 2178)       IF( nbondi == 2 .AND.   &      ! Cyclic east-west
( 2179)          &    (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
( 2180)          pt2d(1-jpr2di:     1    ,:) = pt2d(jpim1-jpr2di:  jpim1 ,:)
( 2181)          pt2d(   jpi  :jpi+jpr2di,:) = pt2d(     2      :2+jpr2di,:)
( 2182) 
( 2183)       ELSE                           ! ... closed
( 2184)          SELECT CASE ( cd_type )
( 2185)          CASE ( 'T', 'U', 'V', 'W' , 'I' )
( 2186)             pt2d(  1-jpr2di   :jpreci    ,:) = 0.e0
( 2187)             pt2d(nlci-jpreci+1:jpi+jpr2di,:) = 0.e0
( 2188)          CASE ( 'F' )
( 2189)             pt2d(nlci-jpreci+1:jpi+jpr2di,:) = 0.e0
( 2190)          END SELECT
( 2191)       ENDIF
( 2192) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 39

( 2193)       !                                        ! North-South boundaries
( 2194)       !                                        ! ======================
( 2195)       SELECT CASE ( cd_type )
( 2196)       CASE ( 'T', 'U', 'V', 'W' , 'I' )
( 2197)          pt2d(:,  1-jpr2dj   :  jprecj  ) = 0.e0
( 2198)          pt2d(:,nlcj-jprecj+1:jpj+jpr2dj) = 0.e0
( 2199)       CASE ( 'F' )
( 2200)          pt2d(:,nlcj-jprecj+1:jpj+jpr2dj) = 0.e0
( 2201)       END SELECT
( 2202) 
( 2203) 
( 2204)       ! 2. East and west directions
( 2205)       ! ---------------------------
( 2206) 
( 2207)       ! 2.1 Read Dirichlet lateral conditions
( 2208) 
( 2209)       SELECT CASE ( nbondi )
( 2210)       CASE ( -1, 0, 1 )    ! all except 2
( 2211)          iihom = nlci-nreci-jpr2di
( 2212)          DO jl = 1, ipreci
( 2213)             tr2ew(:,jl,1) = pt2d(jpreci+jl,:)
( 2214)             tr2we(:,jl,1) = pt2d(iihom +jl,:)
( 2215)          END DO
( 2216)       END SELECT
( 2217) 
( 2218)       ! 2.2 Migrations
( 2219) 
( 2220) #if defined key_mpp_shmem
( 2221)       !! * SHMEM version
( 2222) 
( 2223)       imigr = ipreci * ( jpj + 2*jpr2dj)
( 2224) 
( 2225)       SELECT CASE ( nbondi )
( 2226)       CASE ( -1 )
( 2227)          CALL shmem_put( tr2we(1-jpr2dj,1,2), tr2we(1,1,1), imigr, noea )
( 2228)       CASE ( 0 )
( 2229)          CALL shmem_put( tr2ew(1-jpr2dj,1,2), tr2ew(1,1,1), imigr, nowe )
( 2230)          CALL shmem_put( tr2we(1-jpr2dj,1,2), tr2we(1,1,1), imigr, noea )
( 2231)       CASE ( 1 )
( 2232)          CALL shmem_put( tr2ew(1-jpr2dj,1,2), tr2ew(1,1,1), imigr, nowe )
( 2233)       END SELECT
( 2234) 
( 2235)       CALL barrier()
( 2236)       CALL shmem_udcflush()
( 2237) 
( 2238) #elif defined key_mpp_mpi
( 2239)       !! * MPI version
( 2240) 
( 2241)       imigr = ipreci * ( jpj + 2*jpr2dj)
( 2242) 
( 2243)       SELECT CASE ( nbondi )
( 2244)       CASE ( -1 )
( 2245)          CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req1 )
( 2246)          CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
( 2247)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2248)       CASE ( 0 )
( 2249)          CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
( 2250)          CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req2 )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 40

( 2251)          CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
( 2252)          CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
( 2253)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2254)          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2255)       CASE ( 1 )
( 2256)          CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
( 2257)          CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
( 2258)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2259)       END SELECT
( 2260) 
( 2261) #endif
( 2262) 
( 2263)       ! 2.3 Write Dirichlet lateral conditions
( 2264) 
( 2265)       iihom = nlci - jpreci
( 2266) 
( 2267)       SELECT CASE ( nbondi )
( 2268)       CASE ( -1 )
( 2269)          DO jl = 1, ipreci
( 2270)             pt2d(iihom+jl,:) = tr2ew(:,jl,2)
( 2271)          END DO
( 2272)       CASE ( 0 )
( 2273)          DO jl = 1, ipreci
( 2274)             pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
( 2275)             pt2d( iihom+jl,:) = tr2ew(:,jl,2)
( 2276)          END DO
( 2277)       CASE ( 1 )
( 2278)          DO jl = 1, ipreci
( 2279)             pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
( 2280)          END DO
( 2281)       END SELECT
( 2282) 
( 2283) 
( 2284)       ! 3. North and south directions
( 2285)       ! -----------------------------
( 2286) 
( 2287)       ! 3.1 Read Dirichlet lateral conditions
( 2288) 
( 2289)       IF( nbondj /= 2 ) THEN
( 2290)          ijhom = nlcj-nrecj-jpr2dj
( 2291)          DO jl = 1, iprecj
( 2292)             tr2sn(:,jl,1) = pt2d(:,ijhom +jl)
( 2293)             tr2ns(:,jl,1) = pt2d(:,jprecj+jl)
( 2294)          END DO
( 2295)       ENDIF
( 2296) 
( 2297)       ! 3.2 Migrations
( 2298) 
( 2299) #if defined key_mpp_shmem
( 2300)       !! * SHMEM version
( 2301) 
( 2302)       imigr = iprecj * ( jpi + 2*jpr2di )
( 2303) 
( 2304)       SELECT CASE ( nbondj )
( 2305)       CASE ( -1 )
( 2306)          CALL shmem_put( tr2sn(1-jpr2di,1,2), tr2sn(1,1,1), imigr, nono )
( 2307)       CASE ( 0 )
( 2308)          CALL shmem_put( tr2ns(1-jpr2di,1,2), tr2ns(1,1,1), imigr, noso )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 41

( 2309)          CALL shmem_put( tr2sn(1-jpr2di,1,2), tr2sn(1,1,1), imigr, nono )
( 2310)       CASE ( 1 )
( 2311)          CALL shmem_put( tr2ns(1-jpr2di,1,2), tr2ns(1,1,1), imigr, noso )
( 2312)       END SELECT 
( 2313)       CALL barrier()
( 2314)       CALL shmem_udcflush()
( 2315) 
( 2316) #elif defined key_mpp_mpi
( 2317)       !! * MPI version
( 2318) 
( 2319)       imigr = iprecj * ( jpi + 2*jpr2di )
( 2320) 
( 2321)       SELECT CASE ( nbondj )
( 2322)       CASE ( -1 )
( 2323)          CALL mppsend( 4, tr2sn(1-jpr2di,1,1), imigr, nono, ml_req1 )
( 2324)          CALL mpprecv( 3, tr2ns(1-jpr2di,1,2), imigr )
( 2325)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2326)       CASE ( 0 )
( 2327)          CALL mppsend( 3, tr2ns(1-jpr2di,1,1), imigr, noso, ml_req1 )
( 2328)          CALL mppsend( 4, tr2sn(1-jpr2di,1,1), imigr, nono, ml_req2 )
( 2329)          CALL mpprecv( 3, tr2ns(1-jpr2di,1,2), imigr )
( 2330)          CALL mpprecv( 4, tr2sn(1-jpr2di,1,2), imigr )
( 2331)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2332)          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2333)       CASE ( 1 )
( 2334)          CALL mppsend( 3, tr2ns(1-jpr2di,1,1), imigr, noso, ml_req1 )
( 2335)          CALL mpprecv( 4, tr2sn(1-jpr2di,1,2), imigr )
( 2336)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2337)       END SELECT
( 2338)   
( 2339) #endif
( 2340) 
( 2341)       ! 3.3 Write Dirichlet lateral conditions
( 2342) 
( 2343)       ijhom = nlcj - jprecj  
( 2344) 
( 2345)       SELECT CASE ( nbondj )
( 2346)       CASE ( -1 )
( 2347)          DO jl = 1, iprecj
( 2348)             pt2d(:,ijhom+jl) = tr2ns(:,jl,2)
( 2349)          END DO
( 2350)       CASE ( 0 )
( 2351)          DO jl = 1, iprecj
( 2352)             pt2d(:,jl-jpr2dj) = tr2sn(:,jl,2)
( 2353)             pt2d(:,ijhom+jl ) = tr2ns(:,jl,2)
( 2354)          END DO
( 2355)       CASE ( 1 ) 
( 2356)          DO jl = 1, iprecj
( 2357)             pt2d(:,jl-jpr2dj) = tr2sn(:,jl,2)
( 2358)          END DO
( 2359)       END SELECT 
( 2360)   
( 2361) 
( 2362)       ! 4. north fold treatment
( 2363)       ! -----------------------
( 2364)   
( 2365)       ! 4.1 treatment without exchange (jpni odd)
( 2366)       






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 42

( 2367)       SELECT CASE ( jpni )
( 2368)   
( 2369)       CASE ( 1 ) ! only one proc along I, no mpp exchange
( 2370)   
( 2371)          SELECT CASE ( npolj )
( 2372)   
( 2373)          CASE ( 3 , 4 )   !  T pivot
( 2374)             iloc = jpiglo - 2 * ( nimpp - 1 )
( 2375)   
( 2376)             SELECT CASE ( cd_type )
( 2377)   
( 2378)             CASE ( 'T', 'S', 'W' )
( 2379)                DO jl = 0, iprecj-1
( 2380)                   DO ji = 2-jpr2di, nlci+jpr2di
( 2381)                      ijt=iloc-ji+2
( 2382)                      pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-2-jl)
( 2383)                   END DO
( 2384)                END DO
( 2385)                DO ji = nlci/2+1, nlci+jpr2di
( 2386)                   ijt=iloc-ji+2
( 2387)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
( 2388)                END DO
( 2389)  
( 2390)             CASE ( 'U' )
( 2391)                DO jl =0, iprecj-1
( 2392)                   DO ji = 1-jpr2di, nlci-1-jpr2di
( 2393)                      iju=iloc-ji+1
( 2394)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-2-jl)
( 2395)                   END DO
( 2396)                END DO
( 2397)                DO ji = nlci/2, nlci-1+jpr2di
( 2398)                   iju=iloc-ji+1
( 2399)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
( 2400)                END DO
( 2401)   
( 2402)             CASE ( 'V' )
( 2403)                DO jl = -1, iprecj-1
( 2404)                   DO ji = 2-jpr2di, nlci+jpr2di
( 2405)                      ijt=iloc-ji+2
( 2406)                      pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-3-jl)
( 2407)                   END DO
( 2408)                END DO
( 2409)   
( 2410)             CASE ( 'F', 'G' )
( 2411)                DO jl = -1, iprecj-1
( 2412)                   DO ji = 1-jpr2di, nlci-1+jpr2di
( 2413)                      iju=iloc-ji+1
( 2414)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-3-jl)
( 2415)                   END DO
( 2416)                END DO
( 2417)   
( 2418)             CASE ( 'I' )                                  ! ice U-V point
( 2419)                DO jl = 0, iprecj-1
( 2420)                   pt2d(2,nlcj+jl) = psgn * pt2d(3,nlcj-1-jl)
( 2421)                   DO ji = 3, nlci+jpr2di
( 2422)                      iju = iloc - ji + 3
( 2423)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-1-jl)
( 2424)                   END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 43

( 2425)                END DO
( 2426)   
( 2427)             END SELECT
( 2428)   
( 2429)          CASE ( 5 , 6 )                 ! F pivot
( 2430)             iloc=jpiglo-2*(nimpp-1)
( 2431)   
( 2432)             SELECT CASE (cd_type )
( 2433)   
( 2434)             CASE ( 'T', 'S', 'W' )
( 2435)                DO jl = 0, iprecj-1
( 2436)                   DO ji = 1-jpr2di, nlci+jpr2di
( 2437)                      ijt=iloc-ji+1
( 2438)                      pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-1-jl)
( 2439)                   END DO
( 2440)                END DO
( 2441)   
( 2442)             CASE ( 'U' )
( 2443)                DO jl = 0, iprecj-1
( 2444)                   DO ji = 1-jpr2di, nlci-1+jpr2di
( 2445)                      iju=iloc-ji
( 2446)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-1-jl)
( 2447)                   END DO
( 2448)                END DO
( 2449)  
( 2450)             CASE ( 'V' )
( 2451)                DO jl = 0, iprecj-1
( 2452)                   DO ji = 1-jpr2di, nlci+jpr2di
( 2453)                      ijt=iloc-ji+1
( 2454)                      pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-2-jl)
( 2455)                   END DO
( 2456)                END DO 
( 2457)                DO ji = nlci/2+1, nlci+jpr2di
( 2458)                   ijt=iloc-ji+1
( 2459)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
( 2460)                END DO
( 2461)   
( 2462)             CASE ( 'F', 'G' )
( 2463)                DO jl = 0, iprecj-1
( 2464)                   DO ji = 1-jpr2di, nlci-1+jpr2di
( 2465)                      iju=iloc-ji
( 2466)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-2-jl)
( 2467)                   END DO
( 2468)                END DO
( 2469)                DO ji = nlci/2+1, nlci-1+jpr2di
( 2470)                   iju=iloc-ji
( 2471)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
( 2472)                END DO
( 2473)   
( 2474)             CASE ( 'I' )                                  ! ice U-V point
( 2475)                pt2d( 2 ,nlcj) = 0.e0
( 2476)                DO jl = 0, iprecj-1
( 2477)                   DO ji = 2 , nlci-1+jpr2di
( 2478)                      ijt = iloc - ji + 2
( 2479)                      pt2d(ji,nlcj+jl)= 0.5 * ( pt2d(ji,nlcj-1-jl) + psgn * pt2d(ijt,nlcj-1-jl) )
( 2480)                   END DO
( 2481)                END DO
( 2482)   






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 44

( 2483)             END SELECT   ! cd_type
( 2484)   
( 2485)          END SELECT   ! npolj
( 2486) 
( 2487)       CASE DEFAULT   ! more than 1 proc along I
( 2488)          IF( npolj /= 0 )   CALL mpp_lbc_north_e( pt2d, cd_type, psgn )   ! only for northern procs
( 2489)          
( 2490)       END SELECT   ! jpni
( 2491) 
( 2492) 
( 2493)       ! 5. East and west directions
( 2494)       ! ---------------------------
( 2495) 
( 2496)       SELECT CASE ( npolj )
( 2497) 
( 2498)       CASE ( 3, 4, 5, 6 )
( 2499) 
( 2500)          ! 5.1 Read Dirichlet lateral conditions
( 2501) 
( 2502)          SELECT CASE ( nbondi )
( 2503)          CASE ( -1, 0, 1 )
( 2504)             iihom = nlci-nreci-jpr2di
( 2505)             DO jl = 1, ipreci
( 2506)                tr2ew(:,jl,1) = pt2d(jpreci+jl,:)
( 2507)                tr2we(:,jl,1) = pt2d(iihom +jl,:)
( 2508)             END DO
( 2509)          END SELECT
( 2510) 
( 2511)          ! 5.2 Migrations
( 2512) 
( 2513) #if defined key_mpp_shmem
( 2514)          !! * SHMEM version
( 2515) 
( 2516)          imigr = ipreci * ( jpj + 2*jpr2dj )
( 2517) 
( 2518)          SELECT CASE ( nbondi )
( 2519)          CASE ( -1 )
( 2520)             CALL shmem_put( tr2we(1-jpr2dj,1,2), tr2we(1,1,1), imigr, noea )
( 2521)          CASE ( 0 )
( 2522)             CALL shmem_put( tr2ew(1-jpr2dj,1,2), tr2ew(1,1,1), imigr, nowe )
( 2523)             CALL shmem_put( tr2we(1-jpr2dj,1,2), tr2we(1,1,1), imigr, noea )
( 2524)          CASE ( 1 )
( 2525)             CALL shmem_put( tr2ew(1-jpr2dj,1,2), tr2ew(1,1,1), imigr, nowe )
( 2526)          END SELECT
( 2527) 
( 2528)          CALL barrier()
( 2529)          CALL shmem_udcflush()
( 2530)   
( 2531) #elif defined key_mpp_mpi
( 2532)          !! * MPI version
( 2533)   
( 2534)          imigr=ipreci* ( jpj + 2*jpr2dj )
( 2535)   
( 2536)          SELECT CASE ( nbondi )
( 2537)          CASE ( -1 )
( 2538)             CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req1 )
( 2539)             CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
( 2540)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 45

( 2541)          CASE ( 0 )
( 2542)             CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
( 2543)             CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req2 )
( 2544)             CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
( 2545)             CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
( 2546)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2547)             IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2548)          CASE ( 1 )
( 2549)             CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
( 2550)             CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
( 2551)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2552)          END SELECT 
( 2553) #endif
( 2554) 
( 2555)          ! 5.3 Write Dirichlet lateral conditions
( 2556)   
( 2557)          iihom = nlci - jpreci
( 2558)   
( 2559)          SELECT CASE ( nbondi )
( 2560)          CASE ( -1 )
( 2561)             DO jl = 1, ipreci
( 2562)                pt2d(iihom+jl,:) = tr2ew(:,jl,2)
( 2563)             END DO
( 2564)          CASE ( 0 )
( 2565)             DO jl = 1, ipreci
( 2566)                pt2d(jl- jpr2di,:) = tr2we(:,jl,2)
( 2567)                pt2d(iihom+jl,:) = tr2ew(:,jl,2)
( 2568)             END DO
( 2569)          CASE ( 1 )
( 2570)             DO jl = 1, ipreci
( 2571)                pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
( 2572)             END DO
( 2573)          END SELECT 
( 2574)   
( 2575)       END SELECT   ! npolj
( 2576)   
( 2577)    END SUBROUTINE mpp_lnk_2d_e
( 2578) 
( 2579) 
( 2580)    SUBROUTINE mpplnks( ptab )
( 2581)       !!----------------------------------------------------------------------
( 2582)       !!                  ***  routine mpplnks  ***
( 2583)       !!
( 2584)       !! ** Purpose :   Message passing manadgement for add 2d array local boundary
( 2585)       !!
( 2586)       !! ** Method  :   Use mppsend and mpprecv function for passing mask between
( 2587)       !!       processors following neighboring subdomains.
( 2588)       !!            domain parameters
( 2589)       !!                    nlci   : first dimension of the local subdomain
( 2590)       !!                    nlcj   : second dimension of the local subdomain
( 2591)       !!                    nbondi : mark for "east-west local boundary"
( 2592)       !!                    nbondj : mark for "north-south local boundary"
( 2593)       !!                    noea   : number for local neighboring processors 
( 2594)       !!                    nowe   : number for local neighboring processors
( 2595)       !!                    noso   : number for local neighboring processors
( 2596)       !!                    nono   : number for local neighboring processors
( 2597)       !!
( 2598)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 46

( 2599)       !! * Arguments
( 2600)       REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   &
( 2601)          ptab                     ! 2D array
( 2602)   
( 2603)       !! * Local variables
( 2604)       INTEGER ::   ji, jl         ! dummy loop indices
( 2605)       INTEGER ::   &
( 2606)          imigr, iihom, ijhom      ! temporary integers
( 2607)       INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
( 2608)       INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
( 2609)       !!----------------------------------------------------------------------
( 2610) 
( 2611) 
( 2612)       ! 1. north fold treatment
( 2613)       ! -----------------------
( 2614) 
( 2615)       ! 1.1 treatment without exchange (jpni odd)
( 2616)   
( 2617)       SELECT CASE ( npolj )
( 2618)       CASE ( 4 )
( 2619)          DO ji = 1, nlci
( 2620)             ptab(ji,nlcj-2) = ptab(ji,nlcj-2) + t2p1(ji,1,1)
( 2621)          END DO
( 2622)       CASE ( 6 )
( 2623)          DO ji = 1, nlci
( 2624)             ptab(ji,nlcj-1) = ptab(ji,nlcj-1) + t2p1(ji,1,1)
( 2625)          END DO
( 2626) 
( 2627)       ! 1.2 treatment with exchange (jpni greater than 1)
( 2628)       ! 
( 2629)       CASE ( 3 )
( 2630) #if defined key_mpp_shmem
( 2631)   
( 2632)          !! * SHMEN version
( 2633)   
( 2634)          imigr=jprecj*jpi
( 2635)   
( 2636)          CALL shmem_put(t2p1(1,1,2),t2p1(1,1,1),imigr,nono)
( 2637)          CALL barrier()
( 2638)          CALL shmem_udcflush()
( 2639) 
( 2640) #  elif defined key_mpp_mpi
( 2641)        !! * MPI version
( 2642) 
( 2643)        imigr=jprecj*jpi
( 2644) 
( 2645)        CALL mppsend(3,t2p1(1,1,1),imigr,nono, ml_req1)
( 2646)        CALL mpprecv(3,t2p1(1,1,2),imigr)
( 2647)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2648) 
( 2649) #endif      
( 2650) 
( 2651)        ! Write north fold conditions
( 2652) 
( 2653)        DO ji = 1, nlci
( 2654)           ptab(ji,nlcj-2) = ptab(ji,nlcj-2)+t2p1(ji,1,2)
( 2655)        END DO
( 2656) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 47

( 2657)     CASE ( 5 )
( 2658) 
( 2659) #if defined key_mpp_shmem
( 2660) 
( 2661)        !! * SHMEN version
( 2662) 
( 2663)        imigr=jprecj*jpi
( 2664) 
( 2665)        CALL shmem_put(t2p1(1,1,2),t2p1(1,1,1),imigr,nono)
( 2666)        CALL barrier()
( 2667)        CALL shmem_udcflush()
( 2668) 
( 2669) #  elif defined key_mpp_mpi
( 2670)        !! * Local variables   (MPI version)
( 2671) 
( 2672)        imigr=jprecj*jpi
( 2673) 
( 2674)        CALL mppsend(3,t2p1(1,1,1),imigr,nono, ml_req1)
( 2675)        CALL mpprecv(3,t2p1(1,1,2),imigr)
( 2676)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2677) 
( 2678) #endif      
( 2679) 
( 2680)        ! Write north fold conditions
( 2681) 
( 2682)        DO ji = 1, nlci
( 2683)           ptab(ji,nlcj-1) = ptab(ji,nlcj-1)+t2p1(ji,1,2)
( 2684)        END DO
( 2685) 
( 2686)     END SELECT
( 2687) 
( 2688) 
( 2689)     ! 2. East and west directions
( 2690)     ! ---------------------------
( 2691) 
( 2692)     ! 2.1 Read Dirichlet lateral conditions
( 2693) 
( 2694)     iihom = nlci-jpreci
( 2695) 
( 2696)     SELECT CASE ( nbondi )
( 2697) 
( 2698)     CASE ( -1, 0, 1 )  ! all except 2
( 2699)        DO jl = 1, jpreci
( 2700)              t2ew(:,jl,1) = ptab(  jl    ,:)
( 2701)              t2we(:,jl,1) = ptab(iihom+jl,:)
( 2702)        END DO
( 2703)     END SELECT
( 2704) 
( 2705)     ! 2.2 Migrations
( 2706) 
( 2707) #if defined key_mpp_shmem
( 2708) 
( 2709)     !! * SHMEN version
( 2710) 
( 2711)     imigr=jpreci*jpj
( 2712) 
( 2713)     SELECT CASE ( nbondi )
( 2714) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 48

( 2715)     CASE ( -1 )
( 2716)        CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)
( 2717) 
( 2718)     CASE ( 0 )
( 2719)        CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
( 2720)        CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)
( 2721) 
( 2722)     CASE ( 1 )
( 2723)        CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
( 2724) 
( 2725)     END SELECT
( 2726)     CALL  barrier()
( 2727)     CALL  shmem_udcflush()
( 2728) 
( 2729) #  elif defined key_mpp_mpi
( 2730)     !! * Local variables   (MPI version)
( 2731) 
( 2732)     imigr=jpreci*jpj
( 2733) 
( 2734)     SELECT CASE ( nbondi )
( 2735) 
( 2736)     CASE ( -1 )
( 2737)        CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req1)
( 2738)        CALL mpprecv(1,t2ew(1,1,2),imigr)
( 2739)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2740)     CASE ( 0 )
( 2741)        CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
( 2742)        CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req2)
( 2743)        CALL mpprecv(1,t2ew(1,1,2),imigr)
( 2744)        CALL mpprecv(2,t2we(1,1,2),imigr)
( 2745)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2746)        IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2747) 
( 2748)     CASE ( 1 )
( 2749)        CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
( 2750)        CALL mpprecv(2,t2we(1,1,2),imigr)
( 2751)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2752) 
( 2753)     END SELECT
( 2754) 
( 2755) #endif
( 2756) 
( 2757)     ! 2.3 Write Dirichlet lateral conditions
( 2758) 
( 2759)        iihom = nlci-nreci
( 2760) 
( 2761)     SELECT CASE ( nbondi )
( 2762) 
( 2763)     CASE ( -1 )
( 2764)        DO jl = 1, jpreci
( 2765)              ptab(iihom +jl,:) = ptab(iihom +jl,:)+t2ew(:,jl,2)
( 2766)        END DO
( 2767) 
( 2768)     CASE ( 0 )
( 2769)        DO jl = 1, jpreci
( 2770)              ptab(jpreci+jl,:) = ptab(jpreci+jl,:)+t2we(:,jl,2)
( 2771)              ptab(iihom +jl,:) = ptab(iihom +jl,:)+t2ew(:,jl,2)
( 2772)        END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 49

( 2773) 
( 2774)     CASE ( 1 )
( 2775)        DO jl = 1, jpreci
( 2776)              ptab(jpreci+jl,:) = ptab(jpreci+jl,:)+t2we(:,jl,2)
( 2777)        END DO
( 2778)     END SELECT
( 2779) 
( 2780) 
( 2781)     ! 3. North and south directions
( 2782)     ! -----------------------------
( 2783) 
( 2784)     ! 3.1 Read Dirichlet lateral conditions
( 2785) 
( 2786)     ijhom = nlcj-jprecj
( 2787) 
( 2788)     SELECT CASE ( nbondj )
( 2789) 
( 2790)     CASE ( -1, 0, 1 )
( 2791)        DO jl = 1, jprecj
( 2792)              t2sn(:,jl,1) = ptab(:,ijhom+jl)
( 2793)              t2ns(:,jl,1) = ptab(:,   jl   )
( 2794)        END DO
( 2795) 
( 2796)     END SELECT 
( 2797) 
( 2798)     ! 3.2 Migrations
( 2799) 
( 2800) #if defined key_mpp_shmem
( 2801) 
( 2802)     !! * SHMEN version
( 2803) 
( 2804)     imigr=jprecj*jpi
( 2805) 
( 2806)     SELECT CASE ( nbondj )
( 2807) 
( 2808)     CASE ( -1 )
( 2809)        CALL shmem_put(t2sn(1,1,2),t2sn(1,1,1),imigr,nono)
( 2810) 
( 2811)     CASE ( 0 )
( 2812)        CALL shmem_put(t2ns(1,1,2),t2ns(1,1,1),imigr,noso)
( 2813)        CALL shmem_put(t2sn(1,1,2),t2sn(1,1,1),imigr,nono)
( 2814) 
( 2815)     CASE ( 1 )
( 2816)        CALL shmem_put(t2ns(1,1,2),t2ns(1,1,1),imigr,noso)
( 2817) 
( 2818)     END SELECT
( 2819)     CALL  barrier()
( 2820)     CALL  shmem_udcflush()
( 2821) 
( 2822) #  elif defined key_mpp_mpi
( 2823)     !! * Local variables   (MPI version)
( 2824) 
( 2825)     imigr=jprecj*jpi
( 2826) 
( 2827)     SELECT CASE ( nbondj )
( 2828) 
( 2829)     CASE ( -1 )
( 2830)        CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req1)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 50

( 2831)        CALL mpprecv(3,t2ns(1,1,2),imigr)
( 2832)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2833) 
( 2834)     CASE ( 0 )
( 2835)        CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
( 2836)        CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req2)
( 2837)        CALL mpprecv(3,t2ns(1,1,2),imigr)
( 2838)        CALL mpprecv(4,t2sn(1,1,2),imigr)
( 2839)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2840)        IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2841) 
( 2842)     CASE ( 1 )
( 2843)        CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
( 2844)        CALL mpprecv(4,t2sn(1,1,2),imigr)
( 2845)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2846)     END SELECT
( 2847) 
( 2848) #endif
( 2849) 
( 2850)     ! 3.3 Write Dirichlet lateral conditions
( 2851) 
( 2852)        ijhom = nlcj-nrecj
( 2853) 
( 2854)     SELECT CASE ( nbondj )
( 2855) 
( 2856)     CASE ( -1 )
( 2857)        DO jl = 1, jprecj
( 2858)              ptab(:,ijhom +jl) = ptab(:,ijhom +jl)+t2ns(:,jl,2)
( 2859)        END DO
( 2860) 
( 2861)     CASE ( 0 )
( 2862)        DO jl = 1, jprecj
( 2863)              ptab(:,jprecj+jl) = ptab(:,jprecj+jl)+t2sn(:,jl,2)
( 2864)              ptab(:,ijhom +jl) = ptab(:,ijhom +jl)+t2ns(:,jl,2)
( 2865)        END DO
( 2866) 
( 2867)     CASE ( 1 ) 
( 2868)        DO jl = 1, jprecj
( 2869)              ptab(:,jprecj+jl) = ptab(:,jprecj+jl)+t2sn(:,jl,2)
( 2870)        END DO
( 2871) 
( 2872)     END SELECT
( 2873) 
( 2874)   END SUBROUTINE mpplnks
( 2875) 
( 2876) 
( 2877)    SUBROUTINE mppsend( ktyp, pmess, kbytes, kdest, md_req)
( 2878)       !!----------------------------------------------------------------------
( 2879)       !!                  ***  routine mppsend  ***
( 2880)       !!                   
( 2881)       !! ** Purpose :   Send messag passing array
( 2882)       !!
( 2883)       !!----------------------------------------------------------------------
( 2884)       !! * Arguments
( 2885)       REAL(wp), INTENT(inout) ::   pmess(*)       ! array of real
( 2886)       INTEGER , INTENT( in  ) ::   kbytes,     &  ! size of the array pmess
( 2887)          &                         kdest ,     &  ! receive process number
( 2888)          &                         ktyp,       &  ! Tag of the message






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 51

( 2889)          &                         md_req         ! Argument for isend
( 2890)       !!----------------------------------------------------------------------
( 2891) #if defined key_mpp_shmem
( 2892)       !! * SHMEM version  :    routine not used
( 2893) 
( 2894) #elif defined key_mpp_mpi
( 2895)       !! * MPI version
( 2896)       INTEGER ::   iflag
( 2897) 
( 2898)       SELECT CASE ( c_mpi_send )
( 2899)       CASE ( 'S' )                ! Standard mpi send (blocking)
( 2900)          CALL mpi_send ( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
( 2901)             &                          mpi_comm_opa, iflag )
( 2902)       CASE ( 'B' )                ! Buffer mpi send (blocking)
( 2903)          CALL mpi_bsend( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
( 2904)             &                          mpi_comm_opa, iflag )
( 2905)       CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
( 2906)          ! Be carefull, one more argument here : the mpi request identifier..
( 2907)          CALL mpi_isend( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
( 2908)             &                          mpi_comm_opa, md_req, iflag )
( 2909)       END SELECT
( 2910) #endif
( 2911) 
( 2912)    END SUBROUTINE mppsend
( 2913) 
( 2914) 
( 2915)    SUBROUTINE mpprecv( ktyp, pmess, kbytes )
( 2916)       !!----------------------------------------------------------------------
( 2917)       !!                  ***  routine mpprecv  ***
( 2918)       !!
( 2919)       !! ** Purpose :   Receive messag passing array
( 2920)       !!
( 2921)       !!----------------------------------------------------------------------
( 2922)       !! * Arguments
( 2923)       REAL(wp), INTENT(inout) ::   pmess(*)       ! array of real
( 2924)       INTEGER , INTENT( in  ) ::   kbytes,     &  ! suze of the array pmess
( 2925)          &                         ktyp           ! Tag of the recevied message
( 2926)       !!----------------------------------------------------------------------
( 2927) #if defined key_mpp_shmem
( 2928)       !! * SHMEM version  :    routine not used
( 2929) 
( 2930) #  elif defined key_mpp_mpi
( 2931)       !! * MPI version
( 2932)       INTEGER :: istatus(mpi_status_size)
( 2933)       INTEGER :: iflag
( 2934) 
( 2935)       CALL mpi_recv( pmess, kbytes, mpi_double_precision, mpi_any_source, ktyp,   &
( 2936)          &                          mpi_comm_opa, istatus, iflag )
( 2937) #endif
( 2938) 
( 2939)    END SUBROUTINE mpprecv
( 2940) 
( 2941) 
( 2942)    SUBROUTINE mppgather( ptab, kp, pio )
( 2943)       !!----------------------------------------------------------------------
( 2944)       !!                   ***  routine mppgather  ***
( 2945)       !!                   
( 2946)       !! ** Purpose :   Transfert between a local subdomain array and a work 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 52

( 2947)       !!     array which is distributed following the vertical level.
( 2948)       !!
( 2949)       !! ** Method  :
( 2950)       !!
( 2951)       !!----------------------------------------------------------------------
( 2952)       !! * Arguments
( 2953)       REAL(wp), DIMENSION(jpi,jpj),       INTENT( in  ) ::   ptab   ! subdomain input array
( 2954)       INTEGER ,                           INTENT( in  ) ::   kp     ! record length
( 2955)       REAL(wp), DIMENSION(jpi,jpj,jpnij), INTENT( out ) ::   pio    ! subdomain input array
( 2956)       !!---------------------------------------------------------------------
( 2957) #if defined key_mpp_shmem
( 2958)       !! * SHMEM version
( 2959) 
( 2960)       CALL barrier()
( 2961)       CALL shmem_put( pio(1,1,npvm_me+1), ptab, jpi*jpj, kp )
( 2962)       CALL barrier()
( 2963) 
( 2964) #elif defined key_mpp_mpi
( 2965)       !! * Local variables   (MPI version)
( 2966)       INTEGER :: itaille,ierror
( 2967)   
( 2968)       itaille=jpi*jpj
( 2969)       CALL mpi_gather( ptab, itaille, mpi_double_precision, pio, itaille,   &
( 2970)          &                            mpi_double_precision, kp , mpi_comm_opa, ierror ) 
( 2971) #endif
( 2972) 
( 2973)    END SUBROUTINE mppgather
( 2974) 
( 2975) 
( 2976)    SUBROUTINE mppscatter( pio, kp, ptab )
( 2977)       !!----------------------------------------------------------------------
( 2978)       !!                  ***  routine mppscatter  ***
( 2979)       !!
( 2980)       !! ** Purpose :   Transfert between awork array which is distributed 
( 2981)       !!      following the vertical level and the local subdomain array.
( 2982)       !!
( 2983)       !! ** Method :
( 2984)       !!
( 2985)       !!----------------------------------------------------------------------
( 2986)       REAL(wp), DIMENSION(jpi,jpj,jpnij)  ::  pio        ! output array
( 2987)       INTEGER                             ::   kp        ! Tag (not used with MPI
( 2988)       REAL(wp), DIMENSION(jpi,jpj)        ::  ptab       ! subdomain array input
( 2989)       !!---------------------------------------------------------------------
( 2990) #if defined key_mpp_shmem
( 2991)       !! * SHMEM version
( 2992) 
( 2993)       CALL barrier()
( 2994)       CALL shmem_get( ptab, pio(1,1,npvm_me+1), jpi*jpj, kp )
( 2995)       CALL barrier()
( 2996) 
( 2997) #  elif defined key_mpp_mpi
( 2998)       !! * Local variables   (MPI version)
( 2999)       INTEGER :: itaille, ierror
( 3000)   
( 3001)       itaille=jpi*jpj
( 3002) 
( 3003)       CALL mpi_scatter( pio, itaille, mpi_double_precision, ptab, itaille,   &
( 3004)          &                            mpi_double_precision, kp, mpi_comm_opa, ierror )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 53

( 3005) #endif
( 3006) 
( 3007)    END SUBROUTINE mppscatter
( 3008) 
( 3009) 
( 3010)    SUBROUTINE mppisl_a_int( ktab, kdim )
( 3011)       !!----------------------------------------------------------------------
( 3012)       !!                  ***  routine mppisl_a_int  ***
( 3013)       !!                   
( 3014)       !! ** Purpose :   Massively parallel processors
( 3015)       !!                Find the  non zero value
( 3016)       !!
( 3017)       !!----------------------------------------------------------------------
( 3018)       !! * Arguments
( 3019)       INTEGER, INTENT( in  )                  ::   kdim       ! ???
( 3020)       INTEGER, INTENT(inout), DIMENSION(kdim) ::   ktab       ! ???
( 3021)   
( 3022) #if defined key_mpp_shmem
( 3023)       !! * Local variables   (SHMEM version)
( 3024)       INTEGER :: ji
( 3025)       INTEGER, SAVE :: ibool=0
( 3026) 
( 3027)       IF( kdim > jpmppsum ) CALL ctl_stop( 'mppisl_a_int routine : kdim is too big', &
( 3028)            &                               'change jpmppsum dimension in mpp.h' )
( 3029) 
( 3030)       DO ji = 1, kdim
( 3031)          niitab_shmem(ji) = ktab(ji)
( 3032)       END DO
( 3033)       CALL  barrier()
( 3034)       IF(ibool == 0 ) THEN 
( 3035)          CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,kdim,0   &
( 3036)               ,0,N$PES,ni11wrk_shmem,ni11sync_shmem)
( 3037)          CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,kdim,0   &
( 3038)               ,0,N$PES,ni12wrk_shmem,ni12sync_shmem)
( 3039)       ELSE
( 3040)          CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,kdim,0   &
( 3041)               ,0,N$PES,ni21wrk_shmem,ni21sync_shmem)
( 3042)          CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,kdim,0   &
( 3043)               ,0,N$PES,ni22wrk_shmem,ni22sync_shmem)
( 3044)       ENDIF
( 3045)       CALL  barrier()
( 3046)       ibool=ibool+1
( 3047)       ibool=MOD( ibool,2)
( 3048)       DO ji = 1, kdim
( 3049)          IF( ni11tab_shmem(ji) /= 0. ) THEN
( 3050)             ktab(ji) = ni11tab_shmem(ji)
( 3051)          ELSE
( 3052)             ktab(ji) = ni12tab_shmem(ji)
( 3053)          ENDIF
( 3054)       END DO
( 3055)   
( 3056) #  elif defined key_mpp_mpi
( 3057)       !! * Local variables   (MPI version)
( 3058)       LOGICAL  :: lcommute
( 3059)       INTEGER, DIMENSION(kdim) ::   iwork
( 3060)       INTEGER  :: mpi_isl,ierror
( 3061)   
( 3062)       lcommute = .TRUE.






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 54

( 3063)       CALL mpi_op_create( lc_isl, lcommute, mpi_isl, ierror )
( 3064)       CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer   &
( 3065)            , mpi_isl, mpi_comm_opa, ierror )
( 3066)       ktab(:) = iwork(:)
( 3067) #endif
( 3068) 
( 3069)    END SUBROUTINE mppisl_a_int
( 3070) 
( 3071) 
( 3072)    SUBROUTINE mppisl_int( ktab )
( 3073)       !!----------------------------------------------------------------------
( 3074)       !!                  ***  routine mppisl_int  ***
( 3075)       !!                   
( 3076)       !! ** Purpose :   Massively parallel processors
( 3077)       !!                Find the non zero value
( 3078)       !!
( 3079)       !!----------------------------------------------------------------------
( 3080)       !! * Arguments
( 3081)       INTEGER , INTENT( inout ) ::   ktab        ! 
( 3082) 
( 3083) #if defined key_mpp_shmem
( 3084)       !! * Local variables   (SHMEM version)
( 3085)       INTEGER, SAVE :: ibool=0
( 3086) 
( 3087)       niitab_shmem(1) = ktab
( 3088)       CALL  barrier()
( 3089)       IF(ibool == 0 ) THEN 
( 3090)          CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,1,0   &
( 3091)               ,0,N$PES,ni11wrk_shmem,ni11sync_shmem)
( 3092)          CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,1,0   &
( 3093)               ,0,N$PES,ni12wrk_shmem,ni12sync_shmem)
( 3094)       ELSE
( 3095)          CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,1,0   &
( 3096)               ,0,N$PES,ni21wrk_shmem,ni21sync_shmem)
( 3097)          CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,1,0   &
( 3098)               ,0,N$PES,ni22wrk_shmem,ni22sync_shmem)
( 3099)       ENDIF
( 3100)       CALL  barrier()
( 3101)       ibool=ibool+1
( 3102)       ibool=MOD( ibool,2)
( 3103)       IF( ni11tab_shmem(1) /= 0. ) THEN
( 3104)          ktab = ni11tab_shmem(1)
( 3105)       ELSE
( 3106)          ktab = ni12tab_shmem(1)
( 3107)       ENDIF
( 3108)   
( 3109) #  elif defined key_mpp_mpi
( 3110)   
( 3111)       !! * Local variables   (MPI version)
( 3112)       LOGICAL :: lcommute
( 3113)       INTEGER :: mpi_isl,ierror
( 3114)       INTEGER ::   iwork
( 3115)   
( 3116)       lcommute = .TRUE.
( 3117)       CALL mpi_op_create(lc_isl,lcommute,mpi_isl,ierror)
( 3118)       CALL mpi_allreduce(ktab, iwork, 1,mpi_integer   &
( 3119)            ,mpi_isl,mpi_comm_opa,ierror)
( 3120)       ktab = iwork






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 55

( 3121) #endif
( 3122) 
( 3123)    END SUBROUTINE mppisl_int
( 3124) 
( 3125) 
( 3126)    SUBROUTINE mppmax_a_int( ktab, kdim )
( 3127)       !!----------------------------------------------------------------------
( 3128)       !!                  ***  routine mppmax_a_int  ***
( 3129)       !! 
( 3130)       !! ** Purpose :   Find maximum value in an integer layout array
( 3131)       !!
( 3132)       !!----------------------------------------------------------------------
( 3133)       !! * Arguments
( 3134)       INTEGER , INTENT( in  )                  ::   kdim        ! size of array
( 3135)       INTEGER , INTENT(inout), DIMENSION(kdim) ::   ktab        ! input array
( 3136)   
( 3137) #if defined key_mpp_shmem
( 3138)       !! * Local declarations    (SHMEM version)
( 3139)       INTEGER :: ji
( 3140)       INTEGER, SAVE :: ibool=0
( 3141)   
( 3142)       IF( kdim > jpmppsum ) CALL ctl_stop( 'mppmax_a_int routine : kdim is too big', &
( 3143)            &                               'change jpmppsum dimension in mpp.h' )
( 3144)   
( 3145)       DO ji = 1, kdim
( 3146)          niltab_shmem(ji) = ktab(ji)
( 3147)       END DO
( 3148)       CALL  barrier()
( 3149)       IF(ibool == 0 ) THEN 
( 3150)          CALL shmem_int8_max_to_all (niltab_shmem,niltab_shmem,kdim,0,0   &
( 3151)               ,N$PES,nil1wrk_shmem,nil1sync_shmem )
( 3152)       ELSE
( 3153)          CALL shmem_int8_max_to_all (niltab_shmem,niltab_shmem,kdim,0,0   &
( 3154)               ,N$PES,nil2wrk_shmem,nil2sync_shmem )
( 3155)       ENDIF
( 3156)       CALL  barrier()
( 3157)       ibool=ibool+1
( 3158)       ibool=MOD( ibool,2)
( 3159)       DO ji = 1, kdim
( 3160)          ktab(ji) = niltab_shmem(ji)
( 3161)       END DO
( 3162)   
( 3163) #  elif defined key_mpp_mpi
( 3164)   
( 3165)       !! * Local variables   (MPI version)
( 3166)       INTEGER :: ierror
( 3167)       INTEGER, DIMENSION(kdim) ::   iwork
( 3168)   
( 3169)       CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer,   &
( 3170)            &                mpi_max, mpi_comm_opa, ierror )
( 3171)   
( 3172)       ktab(:) = iwork(:)
( 3173) #endif
( 3174) 
( 3175)    END SUBROUTINE mppmax_a_int
( 3176) 
( 3177) 
( 3178)    SUBROUTINE mppmax_int( ktab )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 56

( 3179)       !!----------------------------------------------------------------------
( 3180)       !!                  ***  routine mppmax_int  ***
( 3181)       !!
( 3182)       !! ** Purpose :
( 3183)       !!     Massively parallel processors
( 3184)       !!     Find maximum value in an integer layout array
( 3185)       !!
( 3186)       !!----------------------------------------------------------------------
( 3187)       !! * Arguments
( 3188)       INTEGER, INTENT(inout) ::   ktab      ! ???
( 3189)   
( 3190)       !! * Local declarations
( 3191) 
( 3192) #if defined key_mpp_shmem
( 3193) 
( 3194)       !! * Local variables   (SHMEM version)
( 3195)       INTEGER :: ji
( 3196)       INTEGER, SAVE :: ibool=0
( 3197)   
( 3198)       niltab_shmem(1) = ktab
( 3199)       CALL  barrier()
( 3200)       IF(ibool == 0 ) THEN 
( 3201)          CALL shmem_int8_max_to_all (niltab_shmem,niltab_shmem, 1,0,0   &
( 3202)               ,N$PES,nil1wrk_shmem,nil1sync_shmem )
( 3203)       ELSE
( 3204)          CALL shmem_int8_max_to_all (niltab_shmem,niltab_shmem, 1,0,0   &
( 3205)               ,N$PES,nil2wrk_shmem,nil2sync_shmem )
( 3206)       ENDIF
( 3207)       CALL  barrier()
( 3208)       ibool=ibool+1
( 3209)       ibool=MOD( ibool,2)
( 3210)       ktab = niltab_shmem(1)
( 3211)   
( 3212) #  elif defined key_mpp_mpi
( 3213) 
( 3214)       !! * Local variables   (MPI version)
( 3215)       INTEGER ::  ierror, iwork
( 3216)   
( 3217)       CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
( 3218)            &              ,mpi_max,mpi_comm_opa,ierror)
( 3219)   
( 3220)       ktab = iwork
( 3221) #endif
( 3222) 
( 3223)    END SUBROUTINE mppmax_int
( 3224) 
( 3225) 
( 3226)    SUBROUTINE mppmin_a_int( ktab, kdim )
( 3227)       !!----------------------------------------------------------------------
( 3228)       !!                  ***  routine mppmin_a_int  ***
( 3229)       !! 
( 3230)       !! ** Purpose :   Find minimum value in an integer layout array
( 3231)       !!
( 3232)       !!----------------------------------------------------------------------
( 3233)       !! * Arguments
( 3234)       INTEGER , INTENT( in  )                  ::   kdim        ! size of array
( 3235)       INTEGER , INTENT(inout), DIMENSION(kdim) ::   ktab        ! input array
( 3236)   






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 57

( 3237) #if defined key_mpp_shmem
( 3238)       !! * Local declarations    (SHMEM version)
( 3239)       INTEGER :: ji
( 3240)       INTEGER, SAVE :: ibool=0
( 3241)   
( 3242)       IF( kdim > jpmppsum ) CALL ctl_stop( 'mppmin_a_int routine : kdim is too big', &
( 3243)            &                               'change jpmppsum dimension in mpp.h' )
( 3244)   
( 3245)       DO ji = 1, kdim
( 3246)          niltab_shmem(ji) = ktab(ji)
( 3247)       END DO
( 3248)       CALL  barrier()
( 3249)       IF(ibool == 0 ) THEN 
( 3250)          CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem,kdim,0,0   &
( 3251)               ,N$PES,nil1wrk_shmem,nil1sync_shmem )
( 3252)       ELSE
( 3253)          CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem,kdim,0,0   &
( 3254)               ,N$PES,nil2wrk_shmem,nil2sync_shmem )
( 3255)       ENDIF
( 3256)       CALL  barrier()
( 3257)       ibool=ibool+1
( 3258)       ibool=MOD( ibool,2)
( 3259)       DO ji = 1, kdim
( 3260)          ktab(ji) = niltab_shmem(ji)
( 3261)       END DO
( 3262)   
( 3263) #  elif defined key_mpp_mpi
( 3264)   
( 3265)       !! * Local variables   (MPI version)
( 3266)       INTEGER :: ierror
( 3267)       INTEGER, DIMENSION(kdim) ::   iwork
( 3268)   
( 3269)       CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer,   &
( 3270)            &                mpi_min, mpi_comm_opa, ierror )
( 3271)   
( 3272)       ktab(:) = iwork(:)
( 3273) #endif
( 3274) 
( 3275)    END SUBROUTINE mppmin_a_int
( 3276) 
( 3277) 
( 3278)    SUBROUTINE mppmin_int( ktab )
( 3279)       !!----------------------------------------------------------------------
( 3280)       !!                  ***  routine mppmin_int  ***
( 3281)       !!
( 3282)       !! ** Purpose :
( 3283)       !!     Massively parallel processors
( 3284)       !!     Find minimum value in an integer layout array
( 3285)       !!
( 3286)       !!----------------------------------------------------------------------
( 3287)       !! * Arguments
( 3288)       INTEGER, INTENT(inout) ::   ktab      ! ???
( 3289)   
( 3290)       !! * Local declarations
( 3291) 
( 3292) #if defined key_mpp_shmem
( 3293) 
( 3294)       !! * Local variables   (SHMEM version)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 58

( 3295)       INTEGER :: ji
( 3296)       INTEGER, SAVE :: ibool=0
( 3297)   
( 3298)       niltab_shmem(1) = ktab
( 3299)       CALL  barrier()
( 3300)       IF(ibool == 0 ) THEN 
( 3301)          CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem, 1,0,0   &
( 3302)               ,N$PES,nil1wrk_shmem,nil1sync_shmem )
( 3303)       ELSE
( 3304)          CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem, 1,0,0   &
( 3305)               ,N$PES,nil2wrk_shmem,nil2sync_shmem )
( 3306)       ENDIF
( 3307)       CALL  barrier()
( 3308)       ibool=ibool+1
( 3309)       ibool=MOD( ibool,2)
( 3310)       ktab = niltab_shmem(1)
( 3311)   
( 3312) #  elif defined key_mpp_mpi
( 3313) 
( 3314)       !! * Local variables   (MPI version)
( 3315)       INTEGER ::  ierror, iwork
( 3316)   
( 3317)       CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
( 3318)            &              ,mpi_min,mpi_comm_opa,ierror)
( 3319)   
( 3320)       ktab = iwork
( 3321) #endif
( 3322) 
( 3323)    END SUBROUTINE mppmin_int
( 3324) 
( 3325) 
( 3326)    SUBROUTINE mppsum_a_int( ktab, kdim )
( 3327)       !!----------------------------------------------------------------------
( 3328)       !!                  ***  routine mppsum_a_int  ***
( 3329)       !!                    
( 3330)       !! ** Purpose :   Massively parallel processors
( 3331)       !!                Global integer sum
( 3332)       !!
( 3333)       !!----------------------------------------------------------------------
( 3334)       !! * Arguments
( 3335)       INTEGER, INTENT( in  )                   ::   kdim      ! ???
( 3336)       INTEGER, INTENT(inout), DIMENSION (kdim) ::   ktab      ! ???
( 3337)   
( 3338) #if defined key_mpp_shmem
( 3339) 
( 3340)       !! * Local variables   (SHMEM version)
( 3341)       INTEGER :: ji
( 3342)       INTEGER, SAVE :: ibool=0
( 3343) 
( 3344)       IF( kdim > jpmppsum ) CALL ctl_stop( 'mppsum_a_int routine : kdim is too big', &
( 3345)            &                               'change jpmppsum dimension in mpp.h' )
( 3346) 
( 3347)       DO ji = 1, kdim
( 3348)          nistab_shmem(ji) = ktab(ji)
( 3349)       END DO
( 3350)       CALL  barrier()
( 3351)       IF(ibool == 0 ) THEN 
( 3352)          CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem,kdim,0,0,   &






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 59

( 3353)               N$PES,nis1wrk_shmem,nis1sync_shmem)
( 3354)       ELSE
( 3355)          CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem,kdim,0,0,   &
( 3356)               N$PES,nis2wrk_shmem,nis2sync_shmem)
( 3357)       ENDIF
( 3358)       CALL  barrier()
( 3359)       ibool = ibool + 1
( 3360)       ibool = MOD( ibool, 2 )
( 3361)       DO ji = 1, kdim
( 3362)          ktab(ji) = nistab_shmem(ji)
( 3363)       END DO
( 3364)   
( 3365) #  elif defined key_mpp_mpi
( 3366) 
( 3367)       !! * Local variables   (MPI version)
( 3368)       INTEGER :: ierror
( 3369)       INTEGER, DIMENSION (kdim) ::  iwork
( 3370)   
( 3371)       CALL mpi_allreduce(ktab, iwork,kdim,mpi_integer   &
( 3372)            ,mpi_sum,mpi_comm_opa,ierror)
( 3373)   
( 3374)       ktab(:) = iwork(:)
( 3375) #endif
( 3376) 
( 3377)    END SUBROUTINE mppsum_a_int
( 3378) 
( 3379) 
( 3380)   SUBROUTINE mppsum_int( ktab )
( 3381)     !!----------------------------------------------------------------------
( 3382)     !!                 ***  routine mppsum_int  ***
( 3383)     !!                  
( 3384)     !! ** Purpose :   Global integer sum
( 3385)     !!
( 3386)     !!----------------------------------------------------------------------
( 3387)     !! * Arguments
( 3388)     INTEGER, INTENT(inout) ::   ktab
( 3389) 
( 3390) #if defined key_mpp_shmem
( 3391) 
( 3392)     !! * Local variables   (SHMEM version)
( 3393)     INTEGER, SAVE :: ibool=0
( 3394) 
( 3395)     nistab_shmem(1) = ktab
( 3396)     CALL  barrier()
( 3397)     IF(ibool == 0 ) THEN 
( 3398)        CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem, 1,0,0,   &
( 3399)             N$PES,nis1wrk_shmem,nis1sync_shmem)
( 3400)     ELSE
( 3401)        CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem, 1,0,0,   &
( 3402)             N$PES,nis2wrk_shmem,nis2sync_shmem)
( 3403)     ENDIF
( 3404)     CALL  barrier()
( 3405)     ibool=ibool+1
( 3406)     ibool=MOD( ibool,2)
( 3407)     ktab = nistab_shmem(1)
( 3408) 
( 3409) #  elif defined key_mpp_mpi
( 3410) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 60

( 3411)     !! * Local variables   (MPI version)
( 3412)     INTEGER :: ierror, iwork
( 3413) 
( 3414)     CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
( 3415)          ,mpi_sum,mpi_comm_opa,ierror)
( 3416) 
( 3417)     ktab = iwork
( 3418) 
( 3419) #endif
( 3420) 
( 3421)   END SUBROUTINE mppsum_int
( 3422) 
( 3423) 
( 3424)   SUBROUTINE mppisl_a_real( ptab, kdim )
( 3425)     !!----------------------------------------------------------------------
( 3426)     !!                 ***  routine mppisl_a_real  ***
( 3427)     !!         
( 3428)     !! ** Purpose :   Massively parallel processors
( 3429)     !!           Find the non zero island barotropic stream function value
( 3430)     !!
( 3431)     !!   Modifications:
( 3432)     !!        !  93-09 (M. Imbard)
( 3433)     !!        !  96-05 (j. Escobar)
( 3434)     !!        !  98-05 (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI 
( 3435)     !!----------------------------------------------------------------------
( 3436)     INTEGER , INTENT( in  )                  ::   kdim      ! ???
( 3437)     REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab      ! ???
( 3438) 
( 3439) #if defined key_mpp_shmem
( 3440) 
( 3441)     !! * Local variables   (SHMEM version)
( 3442)     INTEGER :: ji
( 3443)     INTEGER, SAVE :: ibool=0
( 3444) 
( 3445)     IF( kdim > jpmppsum ) CALL ctl_stop( 'mppisl_a_real routine : kdim is too big', &
( 3446)          &                               'change jpmppsum dimension in mpp.h' )
( 3447) 
( 3448)     DO ji = 1, kdim
( 3449)        wiltab_shmem(ji) = ptab(ji)
( 3450)     END DO
( 3451)     CALL  barrier()
( 3452)     IF(ibool == 0 ) THEN 
( 3453)        CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem,kdim,0   &
( 3454)             ,0,N$PES,wi11wrk_shmem,ni11sync_shmem)
( 3455)        CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem,kdim,0   &
( 3456)             ,0,N$PES,wi12wrk_shmem,ni12sync_shmem)
( 3457)     ELSE
( 3458)        CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem,kdim,0   &
( 3459)             ,0,N$PES,wi21wrk_shmem,ni21sync_shmem)
( 3460)        CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem,kdim,0   &
( 3461)             ,0,N$PES,wi22wrk_shmem,ni22sync_shmem)
( 3462)     ENDIF
( 3463)     CALL  barrier()
( 3464)     ibool=ibool+1
( 3465)     ibool=MOD( ibool,2)
( 3466)     DO ji = 1, kdim
( 3467)        IF(wi1tab_shmem(ji) /= 0. ) THEN
( 3468)           ptab(ji) = wi1tab_shmem(ji)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 61

( 3469)        ELSE
( 3470)           ptab(ji) = wi2tab_shmem(ji)
( 3471)        ENDIF
( 3472)     END DO
( 3473) 
( 3474) #  elif defined key_mpp_mpi
( 3475) 
( 3476)     !! * Local variables   (MPI version)
( 3477)     LOGICAL ::   lcommute = .TRUE.
( 3478)     INTEGER ::   mpi_isl, ierror
( 3479)     REAL(wp), DIMENSION(kdim) ::  zwork
( 3480) 
( 3481)     CALL mpi_op_create(lc_isl,lcommute,mpi_isl,ierror)
( 3482)     CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
( 3483)          ,mpi_isl,mpi_comm_opa,ierror)
( 3484)     ptab(:) = zwork(:)
( 3485) 
( 3486) #endif
( 3487) 
( 3488)   END SUBROUTINE mppisl_a_real
( 3489) 
( 3490) 
( 3491)    SUBROUTINE mppisl_real( ptab )
( 3492)       !!----------------------------------------------------------------------
( 3493)       !!                  ***  routine mppisl_real  ***
( 3494)       !!                  
( 3495)       !! ** Purpose :   Massively parallel processors
( 3496)       !!       Find the  non zero island barotropic stream function value
( 3497)       !!
( 3498)       !!     Modifications:
( 3499)       !!        !  93-09 (M. Imbard)
( 3500)       !!        !  96-05 (j. Escobar)
( 3501)       !!        !  98-05 (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI 
( 3502)       !!----------------------------------------------------------------------
( 3503)       REAL(wp), INTENT(inout) ::   ptab
( 3504) 
( 3505) #if defined key_mpp_shmem
( 3506) 
( 3507)       !! * Local variables   (SHMEM version)
( 3508)       INTEGER, SAVE :: ibool=0
( 3509) 
( 3510)       wiltab_shmem(1) = ptab
( 3511)       CALL  barrier()
( 3512)       IF(ibool == 0 ) THEN 
( 3513)          CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem, 1,0   &
( 3514)             ,0,N$PES,wi11wrk_shmem,ni11sync_shmem)
( 3515)          CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem, 1,0   &
( 3516)             ,0,N$PES,wi12wrk_shmem,ni12sync_shmem)
( 3517)       ELSE
( 3518)          CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem, 1,0   &
( 3519)             ,0,N$PES,wi21wrk_shmem,ni21sync_shmem)
( 3520)          CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem, 1,0   &
( 3521)             ,0,N$PES,wi22wrk_shmem,ni22sync_shmem)
( 3522)       ENDIF
( 3523)       CALL barrier()
( 3524)       ibool = ibool + 1
( 3525)       ibool = MOD( ibool, 2 )
( 3526)       IF( wi1tab_shmem(1) /= 0. ) THEN






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 62

( 3527)          ptab = wi1tab_shmem(1)
( 3528)       ELSE
( 3529)          ptab = wi2tab_shmem(1)
( 3530)       ENDIF
( 3531) 
( 3532) #  elif defined key_mpp_mpi
( 3533) 
( 3534)       !! * Local variables   (MPI version)
( 3535)       LOGICAL  ::   lcommute = .TRUE.
( 3536)       INTEGER  ::   mpi_isl, ierror
( 3537)       REAL(wp) ::   zwork
( 3538) 
( 3539)       CALL mpi_op_create( lc_isl, lcommute, mpi_isl, ierror )
( 3540)       CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision,   &
( 3541)          &                                mpi_isl  , mpi_comm_opa, ierror )
( 3542)       ptab = zwork
( 3543) 
( 3544) #endif
( 3545) 
( 3546)    END SUBROUTINE mppisl_real
( 3547) 
( 3548) 
( 3549)   FUNCTION lc_isl( py, px, kdim, kdtatyp )
( 3550)     INTEGER :: kdim
( 3551)     REAL(wp), DIMENSION(kdim) ::  px, py
( 3552)     INTEGER :: kdtatyp, ji
( 3553)     INTEGER :: lc_isl
( 3554)     DO ji = 1, kdim
( 3555)        IF( py(ji) /= 0. )   px(ji) = py(ji)
( 3556)     END DO
( 3557)     lc_isl=0
( 3558) 
( 3559)   END FUNCTION lc_isl
( 3560) 
( 3561) 
( 3562)   SUBROUTINE mppmax_a_real( ptab, kdim )
( 3563)     !!----------------------------------------------------------------------
( 3564)     !!                 ***  routine mppmax_a_real  ***
( 3565)     !!                  
( 3566)     !! ** Purpose :   Maximum
( 3567)     !!
( 3568)     !!----------------------------------------------------------------------
( 3569)     !! * Arguments
( 3570)     INTEGER , INTENT( in  )                  ::   kdim
( 3571)     REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab
( 3572) 
( 3573) #if defined key_mpp_shmem
( 3574) 
( 3575)     !! * Local variables   (SHMEM version)
( 3576)     INTEGER :: ji
( 3577)     INTEGER, SAVE :: ibool=0
( 3578) 
( 3579)     IF( kdim > jpmppsum ) CALL ctl_stop( 'mppmax_a_real routine : kdim is too big', &
( 3580)          &                               'change jpmppsum dimension in mpp.h' )
( 3581) 
( 3582)     DO ji = 1, kdim
( 3583)        wintab_shmem(ji) = ptab(ji)
( 3584)     END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 63

( 3585)     CALL  barrier()
( 3586)     IF(ibool == 0 ) THEN 
( 3587)        CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem,kdim,0   &
( 3588)             ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
( 3589)     ELSE
( 3590)        CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem,kdim,0   &
( 3591)             ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
( 3592)     ENDIF
( 3593)     CALL  barrier()
( 3594)     ibool=ibool+1
( 3595)     ibool=MOD( ibool,2)
( 3596)     DO ji = 1, kdim
( 3597)        ptab(ji) = wintab_shmem(ji)
( 3598)     END DO
( 3599) 
( 3600) #  elif defined key_mpp_mpi
( 3601) 
( 3602)     !! * Local variables   (MPI version)
( 3603)     INTEGER :: ierror
( 3604)     REAL(wp), DIMENSION(kdim) ::  zwork
( 3605) 
( 3606)     CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
( 3607)          ,mpi_max,mpi_comm_opa,ierror)
( 3608)     ptab(:) = zwork(:)
( 3609) 
( 3610) #endif
( 3611) 
( 3612)   END SUBROUTINE mppmax_a_real
( 3613) 
( 3614) 
( 3615)   SUBROUTINE mppmax_real( ptab )
( 3616)     !!----------------------------------------------------------------------
( 3617)     !!                  ***  routine mppmax_real  ***
( 3618)     !!                    
( 3619)     !! ** Purpose :   Maximum
( 3620)     !!
( 3621)     !!----------------------------------------------------------------------
( 3622)     !! * Arguments
( 3623)     REAL(wp), INTENT(inout) ::   ptab      ! ???
( 3624) 
( 3625) #if defined key_mpp_shmem
( 3626) 
( 3627)     !! * Local variables   (SHMEM version)
( 3628)     INTEGER, SAVE :: ibool=0
( 3629) 
( 3630)     wintab_shmem(1) = ptab
( 3631)     CALL  barrier()
( 3632)     IF(ibool == 0 ) THEN 
( 3633)        CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem, 1,0   &
( 3634)             ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
( 3635)     ELSE
( 3636)        CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem, 1,0   &
( 3637)             ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
( 3638)     ENDIF
( 3639)     CALL  barrier()
( 3640)     ibool=ibool+1
( 3641)     ibool=MOD( ibool,2)
( 3642)     ptab = wintab_shmem(1)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 64

( 3643) 
( 3644) #  elif defined key_mpp_mpi
( 3645) 
( 3646)     !! * Local variables   (MPI version)
( 3647)     INTEGER  ::   ierror
( 3648)     REAL(wp) ::   zwork
( 3649) 
( 3650)     CALL mpi_allreduce( ptab, zwork  , 1             , mpi_double_precision,   &
( 3651)        &                      mpi_max, mpi_comm_opa, ierror     )
( 3652)     ptab = zwork
( 3653) 
( 3654) #endif
( 3655) 
( 3656)   END SUBROUTINE mppmax_real
( 3657) 
( 3658) 
( 3659)   SUBROUTINE mppmin_a_real( ptab, kdim )
( 3660)     !!----------------------------------------------------------------------
( 3661)     !!                 ***  routine mppmin_a_real  ***
( 3662)     !!                  
( 3663)     !! ** Purpose :   Minimum
( 3664)     !!
( 3665)     !!-----------------------------------------------------------------------
( 3666)     !! * Arguments
( 3667)     INTEGER , INTENT( in  )                  ::   kdim
( 3668)     REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab
( 3669) 
( 3670) #if defined key_mpp_shmem
( 3671) 
( 3672)     !! * Local variables   (SHMEM version)
( 3673)     INTEGER :: ji
( 3674)     INTEGER, SAVE :: ibool=0
( 3675) 
( 3676)     IF( kdim > jpmppsum ) CALL ctl_stop( 'mpprmin routine : kdim is too big', &
( 3677)          &                               'change jpmppsum dimension in mpp.h' )
( 3678) 
( 3679)     DO ji = 1, kdim
( 3680)        wintab_shmem(ji) = ptab(ji)
( 3681)     END DO
( 3682)     CALL  barrier()
( 3683)     IF(ibool == 0 ) THEN 
( 3684)        CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem,kdim,0   &
( 3685)             ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
( 3686)     ELSE
( 3687)        CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem,kdim,0   &
( 3688)             ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
( 3689)     ENDIF
( 3690)     CALL  barrier()
( 3691)     ibool=ibool+1
( 3692)     ibool=MOD( ibool,2)
( 3693)     DO ji = 1, kdim
( 3694)        ptab(ji) = wintab_shmem(ji)
( 3695)     END DO
( 3696) 
( 3697) #  elif defined key_mpp_mpi
( 3698) 
( 3699)     !! * Local variables   (MPI version)
( 3700)     INTEGER :: ierror






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 65

( 3701)     REAL(wp), DIMENSION(kdim) ::   zwork
( 3702) 
( 3703)     CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
( 3704)          ,mpi_min,mpi_comm_opa,ierror)
( 3705)     ptab(:) = zwork(:)
( 3706) 
( 3707) #endif
( 3708) 
( 3709)   END SUBROUTINE mppmin_a_real
( 3710) 
( 3711) 
( 3712)   SUBROUTINE mppmin_real( ptab )
( 3713)     !!----------------------------------------------------------------------
( 3714)     !!                  ***  routine mppmin_real  ***
( 3715)     !! 
( 3716)     !! ** Purpose :   minimum in Massively Parallel Processing
( 3717)     !!                REAL scalar case
( 3718)     !!
( 3719)     !!-----------------------------------------------------------------------
( 3720)     !! * Arguments
( 3721)     REAL(wp), INTENT( inout ) ::   ptab        ! 
( 3722) 
( 3723) #if defined key_mpp_shmem
( 3724) 
( 3725)     !! * Local variables   (SHMEM version)
( 3726)     INTEGER, SAVE :: ibool=0
( 3727) 
( 3728)     wintab_shmem(1) = ptab
( 3729)     CALL  barrier()
( 3730)     IF(ibool == 0 ) THEN 
( 3731)        CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem, 1,0   &
( 3732)             ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
( 3733)     ELSE
( 3734)        CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem, 1,0   &
( 3735)             ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
( 3736)     ENDIF
( 3737)     CALL  barrier()
( 3738)     ibool=ibool+1
( 3739)     ibool=MOD( ibool,2)
( 3740)     ptab = wintab_shmem(1)
( 3741) 
( 3742) #  elif defined key_mpp_mpi
( 3743) 
( 3744)     !! * Local variables   (MPI version)
( 3745)     INTEGER  ::   ierror
( 3746)     REAL(wp) ::   zwork
( 3747) 
( 3748)     CALL mpi_allreduce( ptab, zwork, 1,mpi_double_precision   &
( 3749)          &               ,mpi_min,mpi_comm_opa,ierror)
( 3750)     ptab = zwork
( 3751) 
( 3752) #endif
( 3753) 
( 3754)   END SUBROUTINE mppmin_real
( 3755) 
( 3756) 
( 3757)   SUBROUTINE mppsum_a_real( ptab, kdim )
( 3758)     !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 66

( 3759)     !!                  ***  routine mppsum_a_real  ***
( 3760)     !! 
( 3761)     !! ** Purpose :   global sum in Massively Parallel Processing
( 3762)     !!                REAL ARRAY argument case
( 3763)     !!
( 3764)     !!-----------------------------------------------------------------------
( 3765)     INTEGER , INTENT( in )                     ::   kdim      ! size of ptab
( 3766)     REAL(wp), DIMENSION(kdim), INTENT( inout ) ::   ptab      ! input array
( 3767) 
( 3768) #if defined key_mpp_shmem
( 3769) 
( 3770)     !! * Local variables   (SHMEM version)
( 3771)     INTEGER :: ji
( 3772)     INTEGER, SAVE :: ibool=0
( 3773) 
( 3774)     IF( kdim > jpmppsum ) CALL ctl_stop( 'mppsum_a_real routine : kdim is too big', &
( 3775)          &                               'change jpmppsum dimension in mpp.h' )
( 3776) 
( 3777)     DO ji = 1, kdim
( 3778)        wrstab_shmem(ji) = ptab(ji)
( 3779)     END DO
( 3780)     CALL  barrier()
( 3781)     IF(ibool == 0 ) THEN 
( 3782)        CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem,kdim,0   &
( 3783)             ,0,N$PES,wrs1wrk_shmem,nrs1sync_shmem )
( 3784)     ELSE
( 3785)        CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem,kdim,0   &
( 3786)             ,0,N$PES,wrs2wrk_shmem,nrs2sync_shmem )
( 3787)     ENDIF
( 3788)     CALL  barrier()
( 3789)     ibool=ibool+1
( 3790)     ibool=MOD( ibool,2)
( 3791)     DO ji = 1, kdim
( 3792)        ptab(ji) = wrstab_shmem(ji)
( 3793)     END DO
( 3794) 
( 3795) #  elif defined key_mpp_mpi
( 3796) 
( 3797)     !! * Local variables   (MPI version)
( 3798)     INTEGER                   ::   ierror    ! temporary integer
( 3799)     REAL(wp), DIMENSION(kdim) ::   zwork     ! temporary workspace 
( 3800) 
( 3801)     CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
( 3802)          &              ,mpi_sum,mpi_comm_opa,ierror)
( 3803)     ptab(:) = zwork(:)
( 3804) 
( 3805) #endif
( 3806) 
( 3807)   END SUBROUTINE mppsum_a_real
( 3808) 
( 3809) 
( 3810)   SUBROUTINE mppsum_real( ptab )
( 3811)     !!----------------------------------------------------------------------
( 3812)     !!                  ***  routine mppsum_real  ***
( 3813)     !!              
( 3814)     !! ** Purpose :   global sum in Massively Parallel Processing
( 3815)     !!                SCALAR argument case
( 3816)     !!






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 67

( 3817)     !!-----------------------------------------------------------------------
( 3818)     REAL(wp), INTENT(inout) ::   ptab        ! input scalar
( 3819) 
( 3820) #if defined key_mpp_shmem
( 3821) 
( 3822)     !! * Local variables   (SHMEM version)
( 3823)     INTEGER, SAVE :: ibool=0
( 3824) 
( 3825)     wrstab_shmem(1) = ptab
( 3826)     CALL  barrier()
( 3827)     IF(ibool == 0 ) THEN 
( 3828)        CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem, 1,0   &
( 3829)             ,0,N$PES,wrs1wrk_shmem,nrs1sync_shmem )
( 3830)     ELSE
( 3831)        CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem, 1,0   &
( 3832)             ,0,N$PES,wrs2wrk_shmem,nrs2sync_shmem )
( 3833)     ENDIF
( 3834)     CALL  barrier()
( 3835)     ibool = ibool + 1
( 3836)     ibool = MOD( ibool, 2 )
( 3837)     ptab = wrstab_shmem(1)
( 3838) 
( 3839) #  elif defined key_mpp_mpi
( 3840) 
( 3841)     !! * Local variables   (MPI version)
( 3842)     INTEGER  ::   ierror
( 3843)     REAL(wp) ::   zwork
( 3844) 
( 3845)     CALL mpi_allreduce(ptab, zwork, 1,mpi_double_precision   &
( 3846)          &              ,mpi_sum,mpi_comm_opa,ierror)
( 3847)     ptab = zwork
( 3848) 
( 3849) #endif
( 3850) 
( 3851)   END SUBROUTINE mppsum_real
( 3852) 
( 3853)   SUBROUTINE mppsum_realdd( ptab )
( 3854)     !!----------------------------------------------------------------------
( 3855)     !!                  ***  routine mppsum_realdd ***
( 3856)     !!
( 3857)     !! ** Purpose :   global sum in Massively Parallel Processing
( 3858)     !!                SCALAR argument case for double-double precision
( 3859)     !!
( 3860)     !!-----------------------------------------------------------------------
( 3861)     COMPLEX(wp), INTENT(inout)         :: ptab    ! input scalar
( 3862) 
( 3863) #if defined key_mpp_shmem
( 3864) 
( 3865) COMPILATION PROBLEM: MPPSUM_REALDD IMPLEMENTATION DOES NOT EXIST FOR SHMEM CASE
( 3866) 
( 3867) #  elif defined key_mpp_mpi
( 3868) 
( 3869)     !! * Local variables   (MPI version)
( 3870)     INTEGER  ::    ierror
( 3871)     COMPLEX(wp) :: zwork
( 3872) 
( 3873) ! reduce local sums into global sum
( 3874)     call MPI_ALLREDUCE (ptab, zwork, 1, MPI_DOUBLE_COMPLEX, &






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 68

( 3875)                      MPI_SUMDD, MPI_COMM_WORLD, ierror)
( 3876)     ptab = zwork
( 3877) #endif
( 3878) 
( 3879)   END SUBROUTINE mppsum_realdd
( 3880) 
( 3881)   SUBROUTINE mppsum_a_realdd( ptab, kdim )
( 3882)     !!----------------------------------------------------------------------
( 3883)     !!                  ***  routine mppsum_a_realdd  ***
( 3884)     !!
( 3885)     !! ** Purpose :   global sum in Massively Parallel Processing
( 3886)     !!                COMPLEX ARRAY case for double-double precision
( 3887)     !!
( 3888)     !!-----------------------------------------------------------------------
( 3889)     INTEGER , INTENT( in )                        ::   kdim      ! size of ptab
( 3890)     COMPLEX(wp), DIMENSION(kdim), INTENT( inout ) ::   ptab      ! input array
( 3891) 
( 3892) #if defined key_mpp_shmem
( 3893) 
( 3894) COMPILATION PROBLEM: MPPSUM_A_REALDD IMPLEMENTATION DOES NOT EXIST FOR SHMEM CASE
( 3895) 
( 3896) #  elif defined key_mpp_mpi
( 3897) 
( 3898)     !! * Local variables   (MPI version)
( 3899)     INTEGER                      :: ierror    ! temporary integer
( 3900)     COMPLEX(wp), DIMENSION(kdim) :: zwork     ! temporary workspace
( 3901) 
( 3902) 
( 3903)     call MPI_ALLREDUCE (ptab, zwork, kdim, MPI_DOUBLE_COMPLEX, &
( 3904)                      MPI_SUMDD, MPI_COMM_WORLD, ierror)
( 3905)     ptab(:) = zwork(:)
( 3906) 
( 3907) #endif
( 3908) 
( 3909)   END SUBROUTINE mppsum_a_realdd
( 3910) 
( 3911) 
( 3912) 
( 3913)   SUBROUTINE mpp_minloc2d(ptab, pmask, pmin, ki,kj )
( 3914)     !!------------------------------------------------------------------------
( 3915)     !!             ***  routine mpp_minloc  ***
( 3916)     !!
( 3917)     !! ** Purpose :  Compute the global minimum of an array ptab
( 3918)     !!              and also give its global position
( 3919)     !!
( 3920)     !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
( 3921)     !!
( 3922)     !! ** Arguments : I : ptab =local 2D array
( 3923)     !!                O : pmin = global minimum
( 3924)     !!                O : ki,kj = global position of minimum
( 3925)     !!
( 3926)     !! ** Author : J.M. Molines 10/10/2004
( 3927)     !!--------------------------------------------------------------------------
( 3928) #ifdef key_mpp_shmem
( 3929)     CALL ctl_stop( ' mpp_minloc not yet available in SHMEM' )
( 3930) # elif key_mpp_mpi
( 3931)     !! * Arguments
( 3932)     REAL(wp), DIMENSION (jpi,jpj), INTENT (in)  :: ptab ,& ! Local 2D array






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 69

( 3933)          &                                         pmask   ! Local mask
( 3934)     REAL(wp)                     , INTENT (out) :: pmin    ! Global minimum of ptab
( 3935)     INTEGER                      , INTENT (out) :: ki,kj   ! index of minimum in global frame
( 3936) 
( 3937)     !! * Local variables
( 3938)     REAL(wp) :: zmin   ! local minimum
( 3939)     REAL(wp) ,DIMENSION(2,1) :: zain, zaout
( 3940)     INTEGER, DIMENSION (2)  :: ilocs
( 3941)     INTEGER :: ierror
( 3942) 
( 3943) 
( 3944)     zmin  = MINVAL( ptab(:,:) , mask= pmask == 1.e0 )
( 3945)     ilocs = MINLOC( ptab(:,:) , mask= pmask == 1.e0 )
( 3946) 
( 3947)     ki = ilocs(1) + nimpp - 1
( 3948)     kj = ilocs(2) + njmpp - 1
( 3949) 
( 3950)     zain(1,:)=zmin
( 3951)     zain(2,:)=ki+10000.*kj
( 3952) 
( 3953)     CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_OPA,ierror)
( 3954) 
( 3955)     pmin=zaout(1,1)
( 3956)     kj= INT(zaout(2,1)/10000.)
( 3957)     ki= INT(zaout(2,1) - 10000.*kj )
( 3958) #endif
( 3959) 
( 3960)   END SUBROUTINE mpp_minloc2d
( 3961) 
( 3962) 
( 3963)   SUBROUTINE mpp_minloc3d(ptab, pmask, pmin, ki,kj ,kk)
( 3964)     !!------------------------------------------------------------------------
( 3965)     !!             ***  routine mpp_minloc  ***
( 3966)     !!
( 3967)     !! ** Purpose :  Compute the global minimum of an array ptab
( 3968)     !!              and also give its global position
( 3969)     !!
( 3970)     !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
( 3971)     !!
( 3972)     !! ** Arguments : I : ptab =local 2D array
( 3973)     !!                O : pmin = global minimum
( 3974)     !!                O : ki,kj = global position of minimum
( 3975)     !!
( 3976)     !! ** Author : J.M. Molines 10/10/2004
( 3977)     !!--------------------------------------------------------------------------
( 3978) #ifdef key_mpp_shmem
( 3979)     CALL ctl_stop( ' mpp_minloc not yet available in SHMEM' )
( 3980) # elif key_mpp_mpi
( 3981)     !! * Arguments
( 3982)     REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT (in)  :: ptab ,& ! Local 2D array
( 3983)          &                                         pmask   ! Local mask
( 3984)     REAL(wp)                     , INTENT (out) :: pmin    ! Global minimum of ptab
( 3985)     INTEGER                      , INTENT (out) :: ki,kj,kk ! index of minimum in global frame
( 3986) 
( 3987)     !! * Local variables
( 3988)     REAL(wp) :: zmin   ! local minimum
( 3989)     REAL(wp) ,DIMENSION(2,1) :: zain, zaout
( 3990)     INTEGER, DIMENSION (3)  :: ilocs






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 70

( 3991)     INTEGER :: ierror
( 3992) 
( 3993) 
( 3994)     zmin  = MINVAL( ptab(:,:,:) , mask= pmask == 1.e0 )
( 3995)     ilocs = MINLOC( ptab(:,:,:) , mask= pmask == 1.e0 )
( 3996) 
( 3997)     ki = ilocs(1) + nimpp - 1
( 3998)     kj = ilocs(2) + njmpp - 1
( 3999)     kk = ilocs(3)
( 4000) 
( 4001)     zain(1,:)=zmin
( 4002)     zain(2,:)=ki+10000.*kj+100000000.*kk
( 4003) 
( 4004)     CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_OPA,ierror)
( 4005) 
( 4006)     pmin=zaout(1,1)
( 4007)     kk= INT(zaout(2,1)/100000000.)
( 4008)     kj= INT(zaout(2,1) - kk * 100000000. )/10000
( 4009)     ki= INT(zaout(2,1) - kk * 100000000. -kj * 10000. )
( 4010) #endif
( 4011) 
( 4012)   END SUBROUTINE mpp_minloc3d
( 4013) 
( 4014) 
( 4015)   SUBROUTINE mpp_maxloc2d(ptab, pmask, pmax, ki,kj )
( 4016)     !!------------------------------------------------------------------------
( 4017)     !!             ***  routine mpp_maxloc  ***
( 4018)     !!
( 4019)     !! ** Purpose :  Compute the global maximum of an array ptab
( 4020)     !!              and also give its global position
( 4021)     !!
( 4022)     !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
( 4023)     !!
( 4024)     !! ** Arguments : I : ptab =local 2D array
( 4025)     !!                O : pmax = global maximum
( 4026)     !!                O : ki,kj = global position of maximum
( 4027)     !!
( 4028)     !! ** Author : J.M. Molines 10/10/2004
( 4029)     !!--------------------------------------------------------------------------
( 4030) #ifdef key_mpp_shmem
( 4031)     CALL ctl_stop( ' mpp_maxloc not yet available in SHMEM' )
( 4032) # elif key_mpp_mpi
( 4033)     !! * Arguments
( 4034)     REAL(wp), DIMENSION (jpi,jpj), INTENT (in)  :: ptab ,& ! Local 2D array
( 4035)          &                                         pmask   ! Local mask
( 4036)     REAL(wp)                     , INTENT (out) :: pmax    ! Global maximum of ptab
( 4037)     INTEGER                      , INTENT (out) :: ki,kj   ! index of maximum in global frame
( 4038) 
( 4039)     !! * Local variables
( 4040)     REAL(wp) :: zmax   ! local maximum
( 4041)     REAL(wp) ,DIMENSION(2,1) :: zain, zaout
( 4042)     INTEGER, DIMENSION (2)  :: ilocs
( 4043)     INTEGER :: ierror
( 4044) 
( 4045) 
( 4046)     zmax  = MAXVAL( ptab(:,:) , mask= pmask == 1.e0 )
( 4047)     ilocs = MAXLOC( ptab(:,:) , mask= pmask == 1.e0 )
( 4048) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 71

( 4049)     ki = ilocs(1) + nimpp - 1
( 4050)     kj = ilocs(2) + njmpp - 1
( 4051) 
( 4052)     zain(1,:)=zmax
( 4053)     zain(2,:)=ki+10000.*kj
( 4054) 
( 4055)     CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_OPA,ierror)
( 4056) 
( 4057)     pmax=zaout(1,1)
( 4058)     kj= INT(zaout(2,1)/10000.)
( 4059)     ki= INT(zaout(2,1) - 10000.*kj )
( 4060) #endif
( 4061) 
( 4062)   END SUBROUTINE mpp_maxloc2d
( 4063) 
( 4064)   SUBROUTINE mpp_maxloc3d(ptab, pmask, pmax, ki,kj,kk )
( 4065)     !!------------------------------------------------------------------------
( 4066)     !!             ***  routine mpp_maxloc  ***
( 4067)     !!
( 4068)     !! ** Purpose :  Compute the global maximum of an array ptab
( 4069)     !!              and also give its global position
( 4070)     !!
( 4071)     !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
( 4072)     !!
( 4073)     !! ** Arguments : I : ptab =local 2D array
( 4074)     !!                O : pmax = global maximum
( 4075)     !!                O : ki,kj = global position of maximum
( 4076)     !!
( 4077)     !! ** Author : J.M. Molines 10/10/2004
( 4078)     !!--------------------------------------------------------------------------
( 4079) #ifdef key_mpp_shmem
( 4080)     CALL ctl_stop( ' mpp_maxloc not yet available in SHMEM' )
( 4081) # elif key_mpp_mpi
( 4082)     !! * Arguments
( 4083)     REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT (in)  :: ptab ,& ! Local 2D array
( 4084)          &                                         pmask   ! Local mask
( 4085)     REAL(wp)                     , INTENT (out) :: pmax    ! Global maximum of ptab
( 4086)     INTEGER                      , INTENT (out) :: ki,kj,kk   ! index of maximum in global frame
( 4087) 
( 4088)     !! * Local variables
( 4089)     REAL(wp) :: zmax   ! local maximum
( 4090)     REAL(wp) ,DIMENSION(2,1) :: zain, zaout
( 4091)     INTEGER, DIMENSION (3)  :: ilocs
( 4092)     INTEGER :: ierror
( 4093) 
( 4094) 
( 4095)     zmax  = MAXVAL( ptab(:,:,:) , mask= pmask == 1.e0 )
( 4096)     ilocs = MAXLOC( ptab(:,:,:) , mask= pmask == 1.e0 )
( 4097) 
( 4098)     ki = ilocs(1) + nimpp - 1
( 4099)     kj = ilocs(2) + njmpp - 1
( 4100)     kk = ilocs(3)
( 4101) 
( 4102)     zain(1,:)=zmax
( 4103)     zain(2,:)=ki+10000.*kj+100000000.*kk
( 4104) 
( 4105)     CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_OPA,ierror)
( 4106) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 72

( 4107)     pmax=zaout(1,1)
( 4108)     kk= INT(zaout(2,1)/100000000.)
( 4109)     kj= INT(zaout(2,1) - kk * 100000000. )/10000
( 4110)     ki= INT(zaout(2,1) - kk * 100000000. -kj * 10000. )
( 4111) #endif
( 4112) 
( 4113)   END SUBROUTINE mpp_maxloc3d
( 4114) 
( 4115)   SUBROUTINE mppsync()
( 4116)     !!----------------------------------------------------------------------
( 4117)     !!                  ***  routine mppsync  ***
( 4118)     !!                   
( 4119)     !! ** Purpose :   Massively parallel processors, synchroneous
( 4120)     !!
( 4121)     !!-----------------------------------------------------------------------
( 4122) 
( 4123) #if defined key_mpp_shmem
( 4124) 
( 4125)     !! * Local variables   (SHMEM version)
( 4126)     CALL barrier()
( 4127) 
( 4128) #  elif defined key_mpp_mpi
( 4129) 
( 4130)     !! * Local variables   (MPI version)
( 4131)     INTEGER :: ierror
( 4132) 
( 4133)     CALL mpi_barrier(mpi_comm_opa,ierror)
( 4134) 
( 4135) #endif
( 4136) 
( 4137)   END SUBROUTINE mppsync
( 4138) 
( 4139) 
( 4140)   SUBROUTINE mppstop
( 4141)     !!----------------------------------------------------------------------
( 4142)     !!                  ***  routine mppstop  ***
( 4143)     !!                   
( 4144)     !! ** purpose :   Stop massilively parallel processors method
( 4145)     !!
( 4146)     !!----------------------------------------------------------------------
( 4147)     !! * Local declarations
( 4148)     INTEGER ::   info
( 4149)     !!----------------------------------------------------------------------
( 4150) 
( 4151)     ! 1. Mpp synchroneus
( 4152)     ! ------------------
( 4153) 
( 4154)     CALL mppsync
( 4155) #if defined key_mpp_mpi
( 4156)     CALL mpi_finalize( info )
( 4157) #endif
( 4158) 
( 4159)   END SUBROUTINE mppstop
( 4160) 
( 4161) 
( 4162)   SUBROUTINE mppobc( ptab, kd1, kd2, kl, kk, ktype, kij )
( 4163)     !!----------------------------------------------------------------------
( 4164)     !!                  ***  routine mppobc  ***






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 73

( 4165)     !! 
( 4166)     !! ** Purpose :   Message passing manadgement for open boundary
( 4167)     !!     conditions array
( 4168)     !!
( 4169)     !! ** Method  :   Use mppsend and mpprecv function for passing mask
( 4170)     !!       between processors following neighboring subdomains.
( 4171)     !!       domain parameters
( 4172)     !!                    nlci   : first dimension of the local subdomain
( 4173)     !!                    nlcj   : second dimension of the local subdomain
( 4174)     !!                    nbondi : mark for "east-west local boundary"
( 4175)     !!                    nbondj : mark for "north-south local boundary"
( 4176)     !!                    noea   : number for local neighboring processors 
( 4177)     !!                    nowe   : number for local neighboring processors
( 4178)     !!                    noso   : number for local neighboring processors
( 4179)     !!                    nono   : number for local neighboring processors
( 4180)     !!
( 4181)     !! History :
( 4182)     !!        !  98-07 (J.M. Molines) Open boundary conditions
( 4183)     !!----------------------------------------------------------------------
( 4184)     !! * Arguments
( 4185)     INTEGER , INTENT( in ) ::   &
( 4186)          kd1, kd2,   &  ! starting and ending indices
( 4187)          kl ,        &  ! index of open boundary
( 4188)          kk,         &  ! vertical dimension
( 4189)          ktype,      &  ! define north/south or east/west cdt
( 4190)          !              !  = 1  north/south  ;  = 2  east/west
( 4191)          kij            ! horizontal dimension
( 4192)     REAL(wp), DIMENSION(kij,kk), INTENT( inout )  ::   &
( 4193)          ptab           ! variable array
( 4194) 
( 4195)     !! * Local variables
( 4196)     INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
( 4197)     INTEGER  ::   &
( 4198)          iipt0, iipt1, ilpt1,     &  ! temporary integers
( 4199)          ijpt0, ijpt1,            &  !    "          "
( 4200)          imigr, iihom, ijhom         !    "          "
( 4201)     INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
( 4202)     INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
( 4203)     REAL(wp), DIMENSION(jpi,jpj) ::   &
( 4204)          ztab                        ! temporary workspace
( 4205)     !!----------------------------------------------------------------------
( 4206) 
( 4207) 
( 4208)     ! boundary condition initialization
( 4209)     ! ---------------------------------
( 4210) 
( 4211)     ztab(:,:) = 0.e0
( 4212) 
( 4213)     IF( ktype==1 ) THEN                                  ! north/south boundaries
( 4214)        iipt0 = MAX( 1, MIN(kd1 - nimpp+1, nlci     ) )
( 4215)        iipt1 = MAX( 0, MIN(kd2 - nimpp+1, nlci - 1 ) )
( 4216)        ilpt1 = MAX( 1, MIN(kd2 - nimpp+1, nlci     ) )
( 4217)        ijpt0 = MAX( 1, MIN(kl  - njmpp+1, nlcj     ) )
( 4218)        ijpt1 = MAX( 0, MIN(kl  - njmpp+1, nlcj - 1 ) )
( 4219)     ELSEIF( ktype==2 ) THEN                              ! east/west boundaries
( 4220)        iipt0 = MAX( 1, MIN(kl  - nimpp+1, nlci     ) )
( 4221)        iipt1 = MAX( 0, MIN(kl  - nimpp+1, nlci - 1 ) )
( 4222)        ijpt0 = MAX( 1, MIN(kd1 - njmpp+1, nlcj     ) )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 74

( 4223)        ijpt1 = MAX( 0, MIN(kd2 - njmpp+1, nlcj - 1 ) )
( 4224)        ilpt1 = MAX( 1, MIN(kd2 - njmpp+1, nlcj     ) )
( 4225)     ELSE
( 4226)        CALL ctl_stop( 'mppobc: bad ktype' )
( 4227)     ENDIF
( 4228) 
( 4229)     DO jk = 1, kk
( 4230)        IF( ktype==1 ) THEN                               ! north/south boundaries
( 4231)           DO jj = ijpt0, ijpt1
( 4232)              DO ji = iipt0, iipt1
( 4233)                 ztab(ji,jj) = ptab(ji,jk)
( 4234)              END DO
( 4235)           END DO
( 4236)        ELSEIF( ktype==2 ) THEN                           ! east/west boundaries
( 4237)           DO jj = ijpt0, ijpt1
( 4238)              DO ji = iipt0, iipt1
( 4239)                 ztab(ji,jj) = ptab(jj,jk)
( 4240)              END DO
( 4241)           END DO
( 4242)        ENDIF
( 4243) 
( 4244) 
( 4245)        ! 1. East and west directions
( 4246)        ! ---------------------------
( 4247) 
( 4248)        ! 1.1 Read Dirichlet lateral conditions
( 4249) 
( 4250)        IF( nbondi /= 2 ) THEN
( 4251)           iihom = nlci-nreci
( 4252) 
( 4253)           DO jl = 1, jpreci
( 4254)              t2ew(:,jl,1) = ztab(jpreci+jl,:)
( 4255)              t2we(:,jl,1) = ztab(iihom +jl,:)
( 4256)           END DO
( 4257)        ENDIF
( 4258) 
( 4259)        ! 1.2 Migrations
( 4260) 
( 4261) #if defined key_mpp_shmem
( 4262)        !! *  (SHMEM version)
( 4263)        imigr=jpreci*jpj*jpbyt
( 4264) 
( 4265)        IF( nbondi == -1 ) THEN
( 4266)           CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr/jpbyt, noea )
( 4267)        ELSEIF( nbondi == 0 ) THEN
( 4268)           CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr/jpbyt, nowe )
( 4269)           CALL shmem_put( t2we(1,1,2), t2we(1,1,1), imigr/jpbyt, noea )
( 4270)        ELSEIF( nbondi == 1 ) THEN
( 4271)           CALL shmem_put( t2ew(1,1,2), t2ew(1,1,1), imigr/jpbyt, nowe )
( 4272)        ENDIF
( 4273)        CALL barrier()
( 4274)        CALL shmem_udcflush()
( 4275) 
( 4276) #  elif key_mpp_mpi
( 4277)        !! * (MPI version)
( 4278) 
( 4279)        imigr=jpreci*jpj
( 4280) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 75

( 4281)        IF( nbondi == -1 ) THEN
( 4282)           CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req1)
( 4283)           CALL mpprecv(1,t2ew(1,1,2),imigr)
( 4284)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4285)        ELSEIF( nbondi == 0 ) THEN
( 4286)           CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
( 4287)           CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req2)
( 4288)           CALL mpprecv(1,t2ew(1,1,2),imigr)
( 4289)           CALL mpprecv(2,t2we(1,1,2),imigr)
( 4290)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4291)           IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 4292)        ELSEIF( nbondi == 1 ) THEN
( 4293)           CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
( 4294)           CALL mpprecv(2,t2we(1,1,2),imigr)
( 4295)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4296)        ENDIF
( 4297) #endif
( 4298) 
( 4299) 
( 4300)        ! 1.3 Write Dirichlet lateral conditions
( 4301) 
( 4302)        iihom = nlci-jpreci
( 4303)        IF( nbondi == 0 .OR. nbondi == 1 ) THEN
( 4304)           DO jl = 1, jpreci
( 4305)              ztab(jl,:) = t2we(:,jl,2)
( 4306)           END DO
( 4307)        ENDIF
( 4308) 
( 4309)        IF( nbondi == -1 .OR. nbondi == 0 ) THEN
( 4310)           DO jl = 1, jpreci
( 4311)              ztab(iihom+jl,:) = t2ew(:,jl,2)
( 4312)           END DO
( 4313)        ENDIF
( 4314) 
( 4315) 
( 4316)        ! 2. North and south directions
( 4317)        ! -----------------------------
( 4318) 
( 4319)        ! 2.1 Read Dirichlet lateral conditions
( 4320) 
( 4321)        IF( nbondj /= 2 ) THEN
( 4322)           ijhom = nlcj-nrecj
( 4323)           DO jl = 1, jprecj
( 4324)              t2sn(:,jl,1) = ztab(:,ijhom +jl)
( 4325)              t2ns(:,jl,1) = ztab(:,jprecj+jl)
( 4326)           END DO
( 4327)        ENDIF
( 4328) 
( 4329)        ! 2.2 Migrations
( 4330) 
( 4331) #if defined key_mpp_shmem
( 4332)        !! * SHMEM version
( 4333) 
( 4334)        imigr=jprecj*jpi*jpbyt
( 4335) 
( 4336)        IF( nbondj == -1 ) THEN
( 4337)           CALL shmem_put( t2sn(1,1,2), t2sn(1,1,1), imigr/jpbyt, nono )
( 4338)        ELSEIF( nbondj == 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 76

( 4339)           CALL shmem_put( t2ns(1,1,2), t2ns(1,1,1), imigr/jpbyt, noso )
( 4340)           CALL shmem_put( t2sn(1,1,2), t2sn(1,1,1), imigr/jpbyt, nono )
( 4341)        ELSEIF( nbondj == 1 ) THEN
( 4342)           CALL shmem_put( t2ns(1,1,2), t2ns(1,1,1), imigr/jpbyt, noso )
( 4343)        ENDIF
( 4344)        CALL barrier()
( 4345)        CALL shmem_udcflush()
( 4346) 
( 4347) #  elif key_mpp_mpi
( 4348)        !! * Local variables   (MPI version)
( 4349) 
( 4350)        imigr=jprecj*jpi
( 4351) 
( 4352)        IF( nbondj == -1 ) THEN
( 4353)           CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req1)
( 4354)           CALL mpprecv(3,t2ns(1,1,2),imigr)
( 4355)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4356)        ELSEIF( nbondj == 0 ) THEN
( 4357)           CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
( 4358)           CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req2)
( 4359)           CALL mpprecv(3,t2ns(1,1,2),imigr)
( 4360)           CALL mpprecv(4,t2sn(1,1,2),imigr)
( 4361)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4362)           IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 4363)        ELSEIF( nbondj == 1 ) THEN
( 4364)           CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
( 4365)           CALL mpprecv(4,t2sn(1,1,2),imigr)
( 4366)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4367)        ENDIF
( 4368) 
( 4369) #endif
( 4370) 
( 4371)        ! 2.3 Write Dirichlet lateral conditions
( 4372) 
( 4373)        ijhom = nlcj - jprecj
( 4374)        IF( nbondj == 0 .OR. nbondj == 1 ) THEN
( 4375)           DO jl = 1, jprecj
( 4376)              ztab(:,jl) = t2sn(:,jl,2)
( 4377)           END DO
( 4378)        ENDIF
( 4379) 
( 4380)        IF( nbondj == 0 .OR. nbondj == -1 ) THEN
( 4381)           DO jl = 1, jprecj
( 4382)              ztab(:,ijhom+jl) = t2ns(:,jl,2)
( 4383)           END DO
( 4384)        ENDIF
( 4385) 
( 4386)        IF( ktype==1 .AND. kd1 <= jpi+nimpp-1 .AND. nimpp <= kd2 ) THEN
( 4387)           ! north/south boundaries
( 4388)           DO jj = ijpt0,ijpt1
( 4389)              DO ji = iipt0,ilpt1
( 4390)                 ptab(ji,jk) = ztab(ji,jj)  
( 4391)              END DO
( 4392)           END DO
( 4393)        ELSEIF( ktype==2 .AND. kd1 <= jpj+njmpp-1 .AND. njmpp <= kd2 ) THEN
( 4394)           ! east/west boundaries
( 4395)           DO jj = ijpt0,ilpt1
( 4396)              DO ji = iipt0,iipt1






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 77

( 4397)                 ptab(jj,jk) = ztab(ji,jj) 
( 4398)              END DO
( 4399)           END DO
( 4400)        ENDIF
( 4401) 
( 4402)     END DO
( 4403) 
( 4404)   END SUBROUTINE mppobc
( 4405) 
( 4406) 
( 4407)   SUBROUTINE mpp_ini_north
( 4408)     !!----------------------------------------------------------------------
( 4409)     !!               ***  routine mpp_ini_north  ***
( 4410)     !!
( 4411)     !! ** Purpose :   Initialize special communicator for north folding 
( 4412)     !!      condition together with global variables needed in the mpp folding
( 4413)     !!
( 4414)     !! ** Method  : - Look for northern processors
( 4415)     !!              - Put their number in nrank_north
( 4416)     !!              - Create groups for the world processors and the north processors
( 4417)     !!              - Create a communicator for northern processors
( 4418)     !!
( 4419)     !! ** output
( 4420)     !!      njmppmax = njmpp for northern procs
( 4421)     !!      ndim_rank_north = number of processors in the northern line
( 4422)     !!      nrank_north (ndim_rank_north) = number  of the northern procs.
( 4423)     !!      ngrp_world = group ID for the world processors
( 4424)     !!      ngrp_north = group ID for the northern processors
( 4425)     !!      ncomm_north = communicator for the northern procs.
( 4426)     !!      north_root = number (in the world) of proc 0 in the northern comm.
( 4427)     !!
( 4428)     !! History :
( 4429)     !!        !  03-09 (J.M. Molines, MPI only )
( 4430)     !!----------------------------------------------------------------------
( 4431) #ifdef key_mpp_shmem
( 4432)     CALL ctl_stop( ' mpp_ini_north not available in SHMEM' )
( 4433) # elif key_mpp_mpi
( 4434)     INTEGER :: ierr
( 4435)     INTEGER :: jproc
( 4436)     INTEGER :: ii,ji
( 4437)     !!----------------------------------------------------------------------
( 4438) 
( 4439)     njmppmax=MAXVAL(njmppt)
( 4440) 
( 4441)     ! Look for how many procs on the northern boundary
( 4442)     !
( 4443)     ndim_rank_north=0
( 4444)     DO jproc=1,jpnij
( 4445)        IF ( njmppt(jproc) == njmppmax ) THEN
( 4446)           ndim_rank_north = ndim_rank_north + 1
( 4447)        END IF
( 4448)     END DO
( 4449) 
( 4450) 
( 4451)     ! Allocate the right size to nrank_north
( 4452)     !
( 4453)     ALLOCATE(nrank_north(ndim_rank_north))
( 4454) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 78

( 4455)     ! Fill the nrank_north array with proc. number of northern procs.
( 4456)     ! Note : the rank start at 0 in MPI
( 4457)     !
( 4458)     ii=0
( 4459)     DO ji = 1, jpnij
( 4460)        IF ( njmppt(ji) == njmppmax   ) THEN
( 4461)           ii=ii+1
( 4462)           nrank_north(ii)=ji-1
( 4463)        END IF
( 4464)     END DO
( 4465)     ! create the world group
( 4466)     !
( 4467)     CALL MPI_COMM_GROUP(mpi_comm_opa,ngrp_world,ierr)
( 4468)     !
( 4469)     ! Create the North group from the world group
( 4470)     CALL MPI_GROUP_INCL(ngrp_world,ndim_rank_north,nrank_north,ngrp_north,ierr)
( 4471) 
( 4472)     ! Create the North communicator , ie the pool of procs in the north group
( 4473)     !
( 4474)     CALL MPI_COMM_CREATE(mpi_comm_opa,ngrp_north,ncomm_north,ierr)
( 4475) 
( 4476) 
( 4477)     ! find proc number in the world of proc 0 in the north
( 4478)     CALL MPI_GROUP_TRANSLATE_RANKS(ngrp_north,1,0,ngrp_world,north_root,ierr)
( 4479) #endif
( 4480) 
( 4481)   END SUBROUTINE mpp_ini_north
( 4482) 
( 4483) 
( 4484)    SUBROUTINE mpp_lbc_north_3d ( pt3d, cd_type, psgn )
( 4485)       !!---------------------------------------------------------------------
( 4486)       !!                   ***  routine mpp_lbc_north_3d  ***
( 4487)       !!
( 4488)       !! ** Purpose :
( 4489)       !!      Ensure proper north fold horizontal bondary condition in mpp configuration
( 4490)       !!      in case of jpn1 > 1
( 4491)       !!
( 4492)       !! ** Method :
( 4493)       !!      Gather the 4 northern lines of the global domain on 1 processor and 
( 4494)       !!      apply lbc north-fold on this sub array. Then scatter the fold array 
( 4495)       !!      back to the processors.
( 4496)       !!
( 4497)       !! History :
( 4498)       !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
( 4499)       !!                                  from lbc routine
( 4500)       !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
( 4501)       !!----------------------------------------------------------------------
( 4502)       !! * Arguments
( 4503)       CHARACTER(len=1), INTENT( in ) ::   &
( 4504)          cd_type       ! nature of pt3d grid-points
( 4505)          !             !   = T ,  U , V , F or W  gridpoints
( 4506)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
( 4507)          pt3d          ! 3D array on which the boundary condition is applied
( 4508)       REAL(wp), INTENT( in ) ::   &
( 4509)          psgn          ! control of the sign change
( 4510)          !             !   = -1. , the sign is changed if north fold boundary
( 4511)          !             !   =  1. , the sign is kept  if north fold boundary
( 4512) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 79

( 4513)       !! * Local declarations
( 4514)       INTEGER :: ji, jj, jk, jr, jproc
( 4515)       INTEGER :: ierr
( 4516)       INTEGER :: ildi,ilei,iilb
( 4517)       INTEGER :: ijpj,ijpjm1,ij,ijt,iju
( 4518)       INTEGER :: itaille
( 4519)       REAL(wp), DIMENSION(jpiglo,4,jpk) :: ztab
( 4520)       REAL(wp), DIMENSION(jpi,4,jpk,jpni) :: znorthgloio
( 4521)       REAL(wp), DIMENSION(jpi,4,jpk) :: znorthloc
( 4522)       !!----------------------------------------------------------------------
( 4523) 
( 4524)     ! If we get in this routine it s because : North fold condition and mpp with more
( 4525)     !   than one proc across i : we deal only with the North condition
( 4526) 
( 4527)     ! 0. Sign setting
( 4528)     ! ---------------
( 4529) 
( 4530)     ijpj=4
( 4531)     ijpjm1=3
( 4532) 
( 4533)     ! put in znorthloc the last 4 jlines of pt3d
( 4534)     DO jk = 1, jpk 
( 4535)        DO jj = nlcj - ijpj +1, nlcj
( 4536)           ij = jj - nlcj + ijpj
( 4537)           znorthloc(:,ij,jk) = pt3d(:,jj,jk)
( 4538)        END DO
( 4539)     END DO
( 4540) 
( 4541) 
( 4542)     IF (npolj /= 0 ) THEN
( 4543)        ! Build in proc 0 of ncomm_north the znorthgloio
( 4544)        znorthgloio(:,:,:,:) = 0_wp
( 4545) 
( 4546) #ifdef key_mpp_shmem
( 4547)        not done : compiler error
( 4548) #elif defined key_mpp_mpi
( 4549)        itaille=jpi*jpk*ijpj
( 4550)        CALL MPI_GATHER(znorthloc,itaille,MPI_DOUBLE_PRECISION,znorthgloio,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 4551) #endif
( 4552) 
( 4553)     ENDIF
( 4554) 
( 4555)     IF (narea == north_root+1 ) THEN
( 4556)        ! recover the global north array
( 4557)        ztab(:,:,:) = 0_wp
( 4558) 
( 4559)        DO jr = 1, ndim_rank_north
( 4560)           jproc = nrank_north(jr) + 1
( 4561)           ildi  = nldit (jproc)
( 4562)           ilei  = nleit (jproc)
( 4563)           iilb  = nimppt(jproc)
( 4564)           DO jk = 1, jpk 
( 4565)              DO jj = 1, 4
( 4566)                 DO ji = ildi, ilei
( 4567)                    ztab(ji+iilb-1,jj,jk) = znorthgloio(ji,jj,jk,jr)
( 4568)                 END DO
( 4569)              END DO
( 4570)           END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 80

( 4571)        END DO
( 4572) 
( 4573) 
( 4574)        ! Horizontal slab
( 4575)        ! ===============
( 4576) 
( 4577)        DO jk = 1, jpk 
( 4578) 
( 4579) 
( 4580)           ! 2. North-Fold boundary conditions
( 4581)           ! ----------------------------------
( 4582) 
( 4583)           SELECT CASE ( npolj )
( 4584) 
( 4585)           CASE ( 3, 4 )                       ! *  North fold  T-point pivot
( 4586) 
( 4587)              ztab( 1    ,ijpj,jk) = 0.e0
( 4588)              ztab(jpiglo,ijpj,jk) = 0.e0
( 4589) 
( 4590)              SELECT CASE ( cd_type )
( 4591) 
( 4592)              CASE ( 'T' , 'S' , 'W' )                   ! T-, W-point
( 4593)                 DO ji = 2, jpiglo
( 4594)                    ijt = jpiglo-ji+2
( 4595)                    ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-2,jk)
( 4596)                 END DO
( 4597)                 DO ji = jpiglo/2+1, jpiglo
( 4598)                    ijt = jpiglo-ji+2
( 4599)                    ztab(ji,ijpjm1,jk) = psgn * ztab(ijt,ijpjm1,jk)
( 4600)                 END DO
( 4601) 
( 4602)              CASE ( 'U' )                               ! U-point
( 4603)                 DO ji = 1, jpiglo-1
( 4604)                    iju = jpiglo-ji+1
( 4605)                    ztab(ji,ijpj,jk) = psgn * ztab(iju,ijpj-2,jk)
( 4606)                 END DO
( 4607)                 DO ji = jpiglo/2, jpiglo-1
( 4608)                    iju = jpiglo-ji+1
( 4609)                    ztab(ji,ijpjm1,jk) = psgn * ztab(iju,ijpjm1,jk)
( 4610)                 END DO
( 4611) 
( 4612)              CASE ( 'V' )                               ! V-point
( 4613)                 DO ji = 2, jpiglo
( 4614)                    ijt = jpiglo-ji+2
( 4615)                    ztab(ji,ijpj-1,jk) = psgn * ztab(ijt,ijpj-2,jk)
( 4616)                    ztab(ji,ijpj  ,jk) = psgn * ztab(ijt,ijpj-3,jk)
( 4617)                 END DO
( 4618) 
( 4619)              CASE ( 'F' , 'G' )                         ! F-point
( 4620)                 DO ji = 1, jpiglo-1
( 4621)                    iju = jpiglo-ji+1
( 4622)                    ztab(ji,ijpj-1,jk) = psgn * ztab(iju,ijpj-2,jk)
( 4623)                    ztab(ji,ijpj  ,jk) = psgn * ztab(iju,ijpj-3,jk)
( 4624)                 END DO
( 4625) 
( 4626)              END SELECT
( 4627) 
( 4628)           CASE ( 5, 6 )                        ! *  North fold  F-point pivot






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 81

( 4629) 
( 4630)              ztab( 1    ,ijpj,jk) = 0.e0
( 4631)              ztab(jpiglo,ijpj,jk) = 0.e0
( 4632) 
( 4633)              SELECT CASE ( cd_type )
( 4634) 
( 4635)              CASE ( 'T' , 'S' , 'W' )                   ! T-, W-point
( 4636)                 DO ji = 1, jpiglo
( 4637)                    ijt = jpiglo-ji+1
( 4638)                    ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-1,jk)
( 4639)                 END DO
( 4640) 
( 4641)              CASE ( 'U' )                               ! U-point
( 4642)                 DO ji = 1, jpiglo-1
( 4643)                    iju = jpiglo-ji
( 4644)                    ztab(ji,ijpj,jk) = psgn * ztab(iju,ijpj-1,jk)
( 4645)                 END DO
( 4646) 
( 4647)              CASE ( 'V' )                               ! V-point
( 4648)                 DO ji = 1, jpiglo
( 4649)                    ijt = jpiglo-ji+1
( 4650)                    ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-2,jk)
( 4651)                 END DO
( 4652)                 DO ji = jpiglo/2+1, jpiglo
( 4653)                    ijt = jpiglo-ji+1
( 4654)                    ztab(ji,ijpjm1,jk) = psgn * ztab(ijt,ijpjm1,jk)
( 4655)                 END DO
( 4656) 
( 4657)              CASE ( 'F' , 'G' )                         ! F-point
( 4658)                 DO ji = 1, jpiglo-1
( 4659)                    iju = jpiglo-ji
( 4660)                    ztab(ji,ijpj  ,jk) = psgn * ztab(iju,ijpj-2,jk)
( 4661)                 END DO
( 4662)                 DO ji = jpiglo/2+1, jpiglo-1
( 4663)                    iju = jpiglo-ji
( 4664)                    ztab(ji,ijpjm1,jk) = psgn * ztab(iju,ijpjm1,jk)
( 4665)                 END DO
( 4666) 
( 4667)              END SELECT
( 4668) 
( 4669)           CASE DEFAULT                           ! *  closed
( 4670) 
( 4671)              SELECT CASE ( cd_type) 
( 4672) 
( 4673)              CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
( 4674)                 ztab(:, 1  ,jk) = 0.e0
( 4675)                 ztab(:,ijpj,jk) = 0.e0
( 4676) 
( 4677)              CASE ( 'F' )                               ! F-point
( 4678)                 ztab(:,ijpj,jk) = 0.e0
( 4679) 
( 4680)              END SELECT
( 4681) 
( 4682)           END SELECT
( 4683) 
( 4684)           !     End of slab
( 4685)           !     ===========
( 4686) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 82

( 4687)        END DO
( 4688) 
( 4689)        !! Scatter back to pt3d
( 4690)        DO jr = 1, ndim_rank_north
( 4691)           jproc=nrank_north(jr)+1
( 4692)           ildi=nldit (jproc)
( 4693)           ilei=nleit (jproc)
( 4694)           iilb=nimppt(jproc)
( 4695)           DO jk=  1, jpk
( 4696)              DO jj=1,ijpj
( 4697)                 DO ji=ildi,ilei
( 4698)                    znorthgloio(ji,jj,jk,jr)=ztab(ji+iilb-1,jj,jk)
( 4699)                 END DO
( 4700)              END DO
( 4701)           END DO
( 4702)        END DO
( 4703) 
( 4704)     ENDIF      ! only done on proc 0 of ncomm_north
( 4705) 
( 4706) #ifdef key_mpp_shmem
( 4707)     not done yet in shmem : compiler error
( 4708) #elif key_mpp_mpi
( 4709)     IF ( npolj /= 0 ) THEN
( 4710)        itaille=jpi*jpk*ijpj
( 4711)        CALL MPI_SCATTER(znorthgloio,itaille,MPI_DOUBLE_PRECISION,znorthloc,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 4712)     ENDIF
( 4713) #endif
( 4714) 
( 4715)     ! put in the last ijpj jlines of pt3d znorthloc
( 4716)     DO jk = 1 , jpk 
( 4717)        DO jj = nlcj - ijpj + 1 , nlcj
( 4718)           ij = jj - nlcj + ijpj
( 4719)           pt3d(:,jj,jk)= znorthloc(:,ij,jk)
( 4720)        END DO
( 4721)     END DO
( 4722) 
( 4723)   END SUBROUTINE mpp_lbc_north_3d
( 4724) 
( 4725) 
( 4726)   SUBROUTINE mpp_lbc_north_2d ( pt2d, cd_type, psgn)
( 4727)     !!---------------------------------------------------------------------
( 4728)     !!                   ***  routine mpp_lbc_north_2d  ***
( 4729)     !!
( 4730)     !! ** Purpose :
( 4731)     !!      Ensure proper north fold horizontal bondary condition in mpp configuration
( 4732)     !!      in case of jpn1 > 1 (for 2d array )
( 4733)     !!
( 4734)     !! ** Method :
( 4735)     !!      Gather the 4 northern lines of the global domain on 1 processor and 
( 4736)     !!      apply lbc north-fold on this sub array. Then scatter the fold array 
( 4737)     !!      back to the processors.
( 4738)     !!
( 4739)     !! History :
( 4740)     !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
( 4741)     !!                                  from lbc routine
( 4742)     !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
( 4743)     !!----------------------------------------------------------------------
( 4744) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 83

( 4745)     !! * Arguments
( 4746)     CHARACTER(len=1), INTENT( in ) ::   &
( 4747)          cd_type       ! nature of pt2d grid-points
( 4748)     !             !   = T ,  U , V , F or W  gridpoints
( 4749)     REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
( 4750)          pt2d          ! 2D array on which the boundary condition is applied
( 4751)     REAL(wp), INTENT( in ) ::   &
( 4752)          psgn          ! control of the sign change
( 4753)     !             !   = -1. , the sign is changed if north fold boundary
( 4754)     !             !   =  1. , the sign is kept  if north fold boundary
( 4755) 
( 4756) 
( 4757)     !! * Local declarations
( 4758) 
( 4759)     INTEGER :: ji, jj,  jr, jproc
( 4760)     INTEGER :: ierr
( 4761)     INTEGER :: ildi,ilei,iilb
( 4762)     INTEGER :: ijpj,ijpjm1,ij,ijt,iju
( 4763)     INTEGER :: itaille
( 4764) 
( 4765)     REAL(wp), DIMENSION(jpiglo,4) :: ztab
( 4766)     REAL(wp), DIMENSION(jpi,4,jpni) :: znorthgloio
( 4767)     REAL(wp), DIMENSION(jpi,4) :: znorthloc
( 4768)     !!----------------------------------------------------------------------
( 4769)     !!  OPA 8.5, LODYC-IPSL (2002)
( 4770)     !!----------------------------------------------------------------------
( 4771)     ! If we get in this routine it s because : North fold condition and mpp with more
( 4772)     !   than one proc across i : we deal only with the North condition
( 4773) 
( 4774)     ! 0. Sign setting
( 4775)     ! ---------------
( 4776) 
( 4777)     ijpj=4
( 4778)     ijpjm1=3
( 4779) 
( 4780) 
( 4781)     ! put in znorthloc the last 4 jlines of pt2d
( 4782)     DO jj = nlcj - ijpj +1, nlcj
( 4783)        ij = jj - nlcj + ijpj
( 4784)        znorthloc(:,ij)=pt2d(:,jj)
( 4785)     END DO
( 4786) 
( 4787)     IF (npolj /= 0 ) THEN
( 4788)        ! Build in proc 0 of ncomm_north the znorthgloio
( 4789)        znorthgloio(:,:,:) = 0_wp
( 4790) #ifdef key_mpp_shmem
( 4791)        not done : compiler error
( 4792) #elif defined key_mpp_mpi
( 4793)        itaille=jpi*ijpj
( 4794)        CALL MPI_GATHER(znorthloc,itaille,MPI_DOUBLE_PRECISION,znorthgloio,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 4795) #endif
( 4796)     ENDIF
( 4797) 
( 4798)     IF (narea == north_root+1 ) THEN
( 4799)        ! recover the global north array
( 4800)        ztab(:,:) = 0_wp
( 4801) 
( 4802)        DO jr = 1, ndim_rank_north






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 84

( 4803)           jproc=nrank_north(jr)+1
( 4804)           ildi=nldit (jproc)
( 4805)           ilei=nleit (jproc)
( 4806)           iilb=nimppt(jproc)
( 4807)           DO jj=1,4
( 4808)              DO ji=ildi,ilei
( 4809)                 ztab(ji+iilb-1,jj)=znorthgloio(ji,jj,jr)
( 4810)              END DO
( 4811)           END DO
( 4812)        END DO
( 4813) 
( 4814) 
( 4815)        ! 2. North-Fold boundary conditions
( 4816)        ! ----------------------------------
( 4817) 
( 4818)        SELECT CASE ( npolj )
( 4819) 
( 4820)        CASE ( 3, 4 )                       ! *  North fold  T-point pivot
( 4821) 
( 4822)           ztab( 1    ,ijpj) = 0.e0
( 4823)           ztab(jpiglo,ijpj) = 0.e0
( 4824) 
( 4825)           SELECT CASE ( cd_type )
( 4826) 
( 4827)           CASE ( 'T' , 'W' , 'S' )                         ! T-, W-point
( 4828)              DO ji = 2, jpiglo
( 4829)                 ijt = jpiglo-ji+2
( 4830)                 ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-2)
( 4831)              END DO
( 4832)              DO ji = jpiglo/2+1, jpiglo
( 4833)                 ijt = jpiglo-ji+2
( 4834)                 ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
( 4835)              END DO
( 4836) 
( 4837)           CASE ( 'U' )                                     ! U-point
( 4838)              DO ji = 1, jpiglo-1
( 4839)                 iju = jpiglo-ji+1
( 4840)                 ztab(ji,ijpj) = psgn * ztab(iju,ijpj-2)
( 4841)              END DO
( 4842)              DO ji = jpiglo/2, jpiglo-1
( 4843)                 iju = jpiglo-ji+1
( 4844)                 ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
( 4845)              END DO
( 4846) 
( 4847)           CASE ( 'V' )                                     ! V-point
( 4848)              DO ji = 2, jpiglo
( 4849)                 ijt = jpiglo-ji+2
( 4850)                 ztab(ji,ijpj-1) = psgn * ztab(ijt,ijpj-2)
( 4851)                 ztab(ji,ijpj  ) = psgn * ztab(ijt,ijpj-3)
( 4852)              END DO
( 4853) 
( 4854)           CASE ( 'F' , 'G' )                               ! F-point
( 4855)              DO ji = 1, jpiglo-1
( 4856)                 iju = jpiglo-ji+1
( 4857)                 ztab(ji,ijpj-1) = psgn * ztab(iju,ijpj-2)
( 4858)                 ztab(ji,ijpj  ) = psgn * ztab(iju,ijpj-3)
( 4859)              END DO
( 4860) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 85

( 4861)           CASE ( 'I' )                                     ! ice U-V point
( 4862)              ztab(2,ijpj) = psgn * ztab(3,ijpj-1)
( 4863)              DO ji = 3, jpiglo
( 4864)                 iju = jpiglo - ji + 3
( 4865)                 ztab(ji,ijpj) = psgn * ztab(iju,ijpj-1)
( 4866)              END DO
( 4867) 
( 4868)           END SELECT
( 4869) 
( 4870)        CASE ( 5, 6 )                        ! *  North fold  F-point pivot
( 4871) 
( 4872)           ztab( 1 ,ijpj) = 0.e0
( 4873)           ztab(jpiglo,ijpj) = 0.e0
( 4874) 
( 4875)           SELECT CASE ( cd_type )
( 4876) 
( 4877)           CASE ( 'T' , 'W' ,'S' )                          ! T-, W-point
( 4878)              DO ji = 1, jpiglo
( 4879)                 ijt = jpiglo-ji+1
( 4880)                 ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-1)
( 4881)              END DO
( 4882) 
( 4883)           CASE ( 'U' )                                     ! U-point
( 4884)              DO ji = 1, jpiglo-1
( 4885)                 iju = jpiglo-ji
( 4886)                 ztab(ji,ijpj) = psgn * ztab(iju,ijpj-1)
( 4887)              END DO
( 4888) 
( 4889)           CASE ( 'V' )                                     ! V-point
( 4890)              DO ji = 1, jpiglo
( 4891)                 ijt = jpiglo-ji+1
( 4892)                 ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-2)
( 4893)              END DO
( 4894)              DO ji = jpiglo/2+1, jpiglo
( 4895)                 ijt = jpiglo-ji+1
( 4896)                 ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
( 4897)              END DO
( 4898) 
( 4899)           CASE ( 'F' , 'G' )                               ! F-point
( 4900)              DO ji = 1, jpiglo-1
( 4901)                 iju = jpiglo-ji
( 4902)                 ztab(ji,ijpj  ) = psgn * ztab(iju,ijpj-2)
( 4903)              END DO
( 4904)              DO ji = jpiglo/2+1, jpiglo-1
( 4905)                 iju = jpiglo-ji
( 4906)                 ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
( 4907)              END DO
( 4908) 
( 4909)              CASE ( 'I' )                                  ! ice U-V point
( 4910)                 ztab( 2 ,ijpj) = 0.e0
( 4911)                 DO ji = 2 , jpiglo-1
( 4912)                    ijt = jpiglo - ji + 2
( 4913)                    ztab(ji,ijpj)= 0.5 * ( ztab(ji,ijpj-1) + psgn * ztab(ijt,ijpj-1) )
( 4914)                 END DO
( 4915) 
( 4916)           END SELECT
( 4917) 
( 4918)        CASE DEFAULT                           ! *  closed : the code probably never go through






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 86

( 4919) 
( 4920)             SELECT CASE ( cd_type) 
( 4921)   
( 4922)             CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
( 4923)                ztab(:, 1 ) = 0.e0
( 4924)                ztab(:,ijpj) = 0.e0
( 4925) 
( 4926)             CASE ( 'F' )                                   ! F-point
( 4927)                ztab(:,ijpj) = 0.e0
( 4928) 
( 4929)             CASE ( 'I' )                                   ! ice U-V point
( 4930)                ztab(:, 1 ) = 0.e0
( 4931)                ztab(:,ijpj) = 0.e0
( 4932) 
( 4933)             END SELECT
( 4934) 
( 4935)          END SELECT
( 4936) 
( 4937)          !     End of slab
( 4938)          !     ===========
( 4939) 
( 4940)          !! Scatter back to pt2d
( 4941)          DO jr = 1, ndim_rank_north
( 4942)             jproc=nrank_north(jr)+1
( 4943)             ildi=nldit (jproc)
( 4944)             ilei=nleit (jproc)
( 4945)             iilb=nimppt(jproc)
( 4946)             DO jj=1,ijpj
( 4947)                DO ji=ildi,ilei
( 4948)                   znorthgloio(ji,jj,jr)=ztab(ji+iilb-1,jj)
( 4949)                END DO
( 4950)             END DO
( 4951)          END DO
( 4952) 
( 4953)       ENDIF      ! only done on proc 0 of ncomm_north
( 4954) 
( 4955) #ifdef key_mpp_shmem
( 4956)       not done yet in shmem : compiler error
( 4957) #elif key_mpp_mpi
( 4958)       IF ( npolj /= 0 ) THEN
( 4959)          itaille=jpi*ijpj
( 4960)          CALL MPI_SCATTER(znorthgloio,itaille,MPI_DOUBLE_PRECISION,znorthloc,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr
( 4961)       ENDIF
( 4962) #endif
( 4963) 
( 4964)       ! put in the last ijpj jlines of pt2d znorthloc
( 4965)       DO jj = nlcj - ijpj + 1 , nlcj
( 4966)          ij = jj - nlcj + ijpj
( 4967)          pt2d(:,jj)= znorthloc(:,ij)
( 4968)       END DO
( 4969) 
( 4970)    END SUBROUTINE mpp_lbc_north_2d
( 4971) 
( 4972) 
( 4973)    SUBROUTINE mpp_lbc_north_e ( pt2d, cd_type, psgn)
( 4974)     !!---------------------------------------------------------------------
( 4975)     !!                   ***  routine mpp_lbc_north_2d  ***
( 4976)     !!






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 87

( 4977)     !! ** Purpose :
( 4978)     !!      Ensure proper north fold horizontal bondary condition in mpp configuration
( 4979)     !!      in case of jpn1 > 1 (for 2d array with outer extra halo)
( 4980)     !!
( 4981)     !! ** Method :
( 4982)     !!      Gather the 4+2*jpr2dj northern lines of the global domain on 1 processor and 
( 4983)     !!      apply lbc north-fold on this sub array. Then scatter the fold array 
( 4984)     !!      back to the processors.
( 4985)     !!
( 4986)     !! History :
( 4987)     !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
( 4988)     !!                                  from lbc routine
( 4989)     !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
( 4990)     !!   9.0  !  05-09  (R. Benshila )   adapt mpp_lbc_north_2d 
( 4991)     !!----------------------------------------------------------------------
( 4992) 
( 4993)     !! * Arguments
( 4994)     CHARACTER(len=1), INTENT( in ) ::   &
( 4995)          cd_type       ! nature of pt2d grid-points
( 4996)     !             !   = T ,  U , V , F or W  gridpoints
( 4997)     REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj), INTENT( inout ) ::   &
( 4998)          pt2d          ! 2D array on which the boundary condition is applied
( 4999)     REAL(wp), INTENT( in ) ::   &
( 5000)          psgn          ! control of the sign change
( 5001)     !             !   = -1. , the sign is changed if north fold boundary
( 5002)     !             !   =  1. , the sign is kept  if north fold boundary
( 5003) 
( 5004) 
( 5005)     !! * Local declarations
( 5006) 
( 5007)     INTEGER :: ji, jj,  jr, jproc, jl
( 5008)     INTEGER :: ierr
( 5009)     INTEGER :: ildi,ilei,iilb
( 5010)     INTEGER :: ijpj,ijpjm1,ij,ijt,iju, iprecj
( 5011)     INTEGER :: itaille
( 5012) 
( 5013)     REAL(wp), DIMENSION(jpiglo,1-jpr2dj:4+jpr2dj) :: ztab
( 5014)     REAL(wp), DIMENSION(jpi,1-jpr2dj:4+jpr2dj,jpni) :: znorthgloio
( 5015)     REAL(wp), DIMENSION(jpi,1-jpr2dj:4+jpr2dj) :: znorthloc
( 5016) 
( 5017)     ! If we get in this routine it s because : North fold condition and mpp with more
( 5018)     !   than one proc across i : we deal only with the North condition
( 5019) 
( 5020)     ! 0. Sign setting
( 5021)     ! ---------------
( 5022) 
( 5023)     ijpj=4
( 5024)     ijpjm1=3
( 5025)     iprecj = jpr2dj+jprecj
( 5026) 
( 5027)     ! put in znorthloc the last 4 jlines of pt2d
( 5028)     DO jj = nlcj - ijpj + 1 - jpr2dj, nlcj +jpr2dj
( 5029)        ij = jj - nlcj + ijpj
( 5030)        znorthloc(:,ij)=pt2d(1:jpi,jj)
( 5031)     END DO
( 5032) 
( 5033)     IF (npolj /= 0 ) THEN
( 5034)        ! Build in proc 0 of ncomm_north the znorthgloio






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 88

( 5035)        znorthgloio(:,:,:) = 0_wp
( 5036) #ifdef key_mpp_shmem
( 5037)        not done : compiler error
( 5038) #elif defined key_mpp_mpi
( 5039)        itaille=jpi*(ijpj+2*jpr2dj)
( 5040)        CALL MPI_GATHER(znorthloc(1,1-jpr2dj),itaille,MPI_DOUBLE_PRECISION, &
( 5041)                      & znorthgloio(1,1-jpr2dj,1),itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 5042) #endif
( 5043)     ENDIF
( 5044) 
( 5045)     IF (narea == north_root+1 ) THEN
( 5046)        ! recover the global north array
( 5047)        ztab(:,:) = 0_wp
( 5048) 
( 5049)        DO jr = 1, ndim_rank_north
( 5050)           jproc=nrank_north(jr)+1
( 5051)           ildi=nldit (jproc)
( 5052)           ilei=nleit (jproc)
( 5053)           iilb=nimppt(jproc)
( 5054)           DO jj=1-jpr2dj,ijpj+jpr2dj
( 5055)              DO ji=ildi,ilei
( 5056)                 ztab(ji+iilb-1,jj)=znorthgloio(ji,jj,jr)
( 5057)              END DO
( 5058)           END DO
( 5059)        END DO
( 5060) 
( 5061) 
( 5062)        ! 2. North-Fold boundary conditions
( 5063)        ! ----------------------------------
( 5064) 
( 5065)        SELECT CASE ( npolj )
( 5066) 
( 5067)        CASE ( 3, 4 )                       ! *  North fold  T-point pivot
( 5068) 
( 5069)           ztab( 1    ,ijpj:ijpj+jpr2dj) = 0.e0
( 5070)           ztab(jpiglo,ijpj:ijpj+jpr2dj) = 0.e0
( 5071) 
( 5072)           SELECT CASE ( cd_type )
( 5073) 
( 5074)           CASE ( 'T' , 'W' , 'S' )                         ! T-, W-point
( 5075)              DO jl =0, iprecj-1
( 5076)                 DO ji = 2, jpiglo
( 5077)                    ijt = jpiglo-ji+2
( 5078)                    ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-2-jl)
( 5079)                 END DO
( 5080)              END DO
( 5081)              DO ji = jpiglo/2+1, jpiglo
( 5082)                 ijt = jpiglo-ji+2
( 5083)                 ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
( 5084)              END DO
( 5085) 
( 5086)           CASE ( 'U' )                                     ! U-point
( 5087)              DO jl =0, iprecj-1
( 5088)                 DO ji = 1, jpiglo-1
( 5089)                    iju = jpiglo-ji+1
( 5090)                    ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-2-jl)
( 5091)                 END DO
( 5092)              END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 89

( 5093)              DO ji = jpiglo/2, jpiglo-1
( 5094)                 iju = jpiglo-ji+1
( 5095)                 ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
( 5096)              END DO
( 5097) 
( 5098)           CASE ( 'V' )                                     ! V-point
( 5099)             DO jl =-1, iprecj-1
( 5100)                DO ji = 2, jpiglo
( 5101)                   ijt = jpiglo-ji+2
( 5102)                   ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-3-jl)
( 5103)                END DO
( 5104)             END DO
( 5105) 
( 5106)           CASE ( 'F' , 'G' )                               ! F-point
( 5107)             DO jl =-1, iprecj-1
( 5108)                DO ji = 1, jpiglo-1
( 5109)                   iju = jpiglo-ji+1
( 5110)                   ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-3-jl)
( 5111)                END DO
( 5112)              END DO
( 5113) 
( 5114)           CASE ( 'I' )                                     ! ice U-V point
( 5115)              DO jl =0, iprecj-1
( 5116)                 ztab(2,ijpj+jl) = psgn * ztab(3,ijpj-1+jl)
( 5117)                 DO ji = 3, jpiglo
( 5118)                    iju = jpiglo - ji + 3
( 5119)                    ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-1-jl)
( 5120)                 END DO
( 5121)              END DO
( 5122) 
( 5123)           END SELECT
( 5124) 
( 5125)        CASE ( 5, 6 )                        ! *  North fold  F-point pivot
( 5126) 
( 5127)           ztab( 1 ,ijpj:ijpj+jpr2dj) = 0.e0
( 5128)           ztab(jpiglo,ijpj:ijpj+jpr2dj) = 0.e0
( 5129) 
( 5130)           SELECT CASE ( cd_type )
( 5131) 
( 5132)           CASE ( 'T' , 'W' ,'S' )                          ! T-, W-point
( 5133)              DO jl = 0, iprecj-1
( 5134)                 DO ji = 1, jpiglo
( 5135)                    ijt = jpiglo-ji+1
( 5136)                    ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-1-jl)
( 5137)                 END DO
( 5138)              END DO
( 5139) 
( 5140)           CASE ( 'U' )                                     ! U-point
( 5141)              DO jl = 0, iprecj-1
( 5142)                 DO ji = 1, jpiglo-1
( 5143)                    iju = jpiglo-ji
( 5144)                    ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-1-jl)
( 5145)                 END DO
( 5146)              END DO
( 5147) 
( 5148)           CASE ( 'V' )                                     ! V-point
( 5149)              DO jl = 0, iprecj-1
( 5150)                 DO ji = 1, jpiglo






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 90

( 5151)                    ijt = jpiglo-ji+1
( 5152)                    ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-2-jl)
( 5153)                 END DO
( 5154)              END DO
( 5155)              DO ji = jpiglo/2+1, jpiglo
( 5156)                 ijt = jpiglo-ji+1
( 5157)                 ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
( 5158)              END DO
( 5159) 
( 5160)           CASE ( 'F' , 'G' )                               ! F-point
( 5161)              DO jl = 0, iprecj-1
( 5162)                 DO ji = 1, jpiglo-1
( 5163)                    iju = jpiglo-ji
( 5164)                    ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-2-jl)
( 5165)                 END DO
( 5166)              END DO
( 5167)              DO ji = jpiglo/2+1, jpiglo-1
( 5168)                 iju = jpiglo-ji
( 5169)                 ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
( 5170)              END DO
( 5171) 
( 5172)              CASE ( 'I' )                                  ! ice U-V point
( 5173)                 ztab( 2 ,ijpj:ijpj+jpr2dj) = 0.e0
( 5174)                 DO jl = 0, jpr2dj
( 5175)                    DO ji = 2 , jpiglo-1
( 5176)                       ijt = jpiglo - ji + 2
( 5177)                       ztab(ji,ijpj+jl)= 0.5 * ( ztab(ji,ijpj-1-jl) + psgn * ztab(ijt,ijpj-1-jl) )
( 5178)                    END DO
( 5179)                 END DO
( 5180) 
( 5181)           END SELECT
( 5182) 
( 5183)        CASE DEFAULT                           ! *  closed : the code probably never go through
( 5184) 
( 5185)             SELECT CASE ( cd_type) 
( 5186)   
( 5187)             CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
( 5188)                ztab(:, 1:1-jpr2dj     ) = 0.e0
( 5189)                ztab(:,ijpj:ijpj+jpr2dj) = 0.e0
( 5190) 
( 5191)             CASE ( 'F' )                                   ! F-point
( 5192)                ztab(:,ijpj:ijpj+jpr2dj) = 0.e0
( 5193) 
( 5194)             CASE ( 'I' )                                   ! ice U-V point
( 5195)                ztab(:, 1:1-jpr2dj     ) = 0.e0
( 5196)                ztab(:,ijpj:ijpj+jpr2dj) = 0.e0
( 5197) 
( 5198)             END SELECT
( 5199) 
( 5200)          END SELECT
( 5201) 
( 5202)          !     End of slab
( 5203)          !     ===========
( 5204) 
( 5205)          !! Scatter back to pt2d
( 5206)          DO jr = 1, ndim_rank_north
( 5207)             jproc=nrank_north(jr)+1
( 5208)             ildi=nldit (jproc)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 91

( 5209)             ilei=nleit (jproc)
( 5210)             iilb=nimppt(jproc)
( 5211)             DO jj=1-jpr2dj,ijpj+jpr2dj
( 5212)                DO ji=ildi,ilei
( 5213)                   znorthgloio(ji,jj,jr)=ztab(ji+iilb-1,jj)
( 5214)                END DO
( 5215)             END DO
( 5216)          END DO
( 5217) 
( 5218)       ENDIF      ! only done on proc 0 of ncomm_north
( 5219) 
( 5220) #ifdef key_mpp_shmem
( 5221)       not done yet in shmem : compiler error
( 5222) #elif key_mpp_mpi
( 5223)       IF ( npolj /= 0 ) THEN
( 5224)          itaille=jpi*(ijpj+2*jpr2dj)
( 5225)          CALL MPI_SCATTER(znorthgloio(1,1-jpr2dj,1),itaille,MPI_DOUBLE_PRECISION, &
( 5226)                         & znorthloc(1,1-jpr2dj),itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 5227)       ENDIF
( 5228) #endif
( 5229) 
( 5230)       ! put in the last ijpj jlines of pt2d znorthloc
( 5231)       DO jj = nlcj - ijpj  -jpr2dj + 1 , nlcj +jpr2dj
( 5232)          ij = jj - nlcj + ijpj 
( 5233)          pt2d(1:jpi,jj)= znorthloc(:,ij)
( 5234)       END DO
( 5235) 
( 5236)    END SUBROUTINE mpp_lbc_north_e
( 5237) 
( 5238) 
( 5239)    !!!!!
( 5240) 
( 5241) 
( 5242)    !! 
( 5243)    !!    This is valid on IBM machine ONLY. 
( 5244)    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
( 5245)    !! mpi_init_opa.f90 : Redefinition du point d'entree MPI_INIT de la bibliotheque
( 5246)    !!                MPI afin de faire, en plus de l'initialisation de
( 5247)    !!                l'environnement MPI, l'allocation d'une zone tampon
( 5248)    !!                qui sera ulterieurement utilisee automatiquement lors
( 5249)    !!                de tous les envois de messages par MPI_BSEND
( 5250)    !!
( 5251)    !! Auteur : CNRS/IDRIS
( 5252)    !! Date   : Tue Nov 13 12:02:14 2001
( 5253)    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
( 5254) 
( 5255)    SUBROUTINE mpi_init_opa(code)
( 5256)       IMPLICIT NONE
( 5257) 
( 5258) !$AGRIF_DO_NOT_TREAT
( 5259) #     include <mpif.h>
( 5260) !$AGRIF_END_DO_NOT_TREAT
( 5261) 
( 5262)       INTEGER                                 :: code,rang,ierr
( 5263)       LOGICAL                                 :: mpi_was_called
( 5264)  
( 5265)       ! La valeur suivante doit etre au moins egale a la taille
( 5266)       ! du plus grand message qui sera transfere dans le programme






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 92

( 5267)       ! (de toute facon, il y aura un message d'erreur si cette
( 5268)       ! valeur s'avere trop petite)
( 5269)       INTEGER                                 :: taille_tampon
( 5270)       CHARACTER(len=9)                        :: taille_tampon_alphanum
( 5271)       REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: tampon
( 5272)  
( 5273)       ! Le point d'entree dans la bibliotheque MPI elle-meme
( 5274)       CALL mpi_initialized(mpi_was_called, code)
( 5275)       IF ( code /= MPI_SUCCESS ) THEN
( 5276)         CALL ctl_stop( ' lib_mpp: Error in routine mpi_initialized' )
( 5277)         CALL mpi_abort( mpi_comm_world, code, ierr )
( 5278)       ENDIF
( 5279) 
( 5280)       IF ( .NOT. mpi_was_called ) THEN
( 5281)          CALL mpi_init(code)
( 5282)          CALL mpi_comm_dup( mpi_comm_world, mpi_comm_opa, code)
( 5283)          IF ( code /= MPI_SUCCESS ) THEN
( 5284)             CALL ctl_stop( ' lib_mpp: Error in routine mpi_comm_dup' )
( 5285)             CALL mpi_abort( mpi_comm_world, code, ierr )
( 5286)          ENDIF
( 5287)       ENDIF
( 5288)       ! La definition de la zone tampon pour les futurs envois
( 5289)       ! par MPI_BSEND (on alloue une fois pour toute cette zone
( 5290)       ! tampon, qui sera automatiquement utilisee lors de chaque
( 5291)       ! appel  a MPI_BSEND).
( 5292)       ! La desallocation sera implicite quand on sortira de
( 5293)       ! l'environnement MPI.
( 5294) 
( 5295)       ! Recuperation de la valeur de la variable d'environnement
( 5296)       ! BUFFER_LENGTH
( 5297)       ! qui, si elle est definie, doit contenir une valeur superieure
( 5298)       ! a  la taille en octets du plus gros message
( 5299)       CALL getenv('BUFFER_LENGTH',taille_tampon_alphanum)
( 5300)   
( 5301)       ! Si la variable BUFFER_LENGTH n'est pas positionnee, on lui met par
( 5302)       ! defaut la plus grande valeur de la variable MP_EAGER_LIMIT, soit
( 5303)       ! 65 536 octets
( 5304)       IF (taille_tampon_alphanum == ' ') THEN
( 5305)          taille_tampon = 65536
( 5306)       ELSE
( 5307)          READ(taille_tampon_alphanum,'(i9)') taille_tampon
( 5308)       END IF
( 5309) 
( 5310)       ! On est limite en mode d'adressage 32 bits a  1750 Mo pour la zone
( 5311)       ! "data" soit 7 segments, c.-a -d. 1750/8 = 210 Mo
( 5312)       IF (taille_tampon > 210000000) THEN
( 5313)          CALL ctl_stop( ' lib_mpp: Attention la valeur BUFFER_LENGTH doit etre <= 210000000' )
( 5314)          CALL mpi_abort(MPI_COMM_WORLD,2,code)
( 5315)       END IF
( 5316) 
( 5317)       CALL mpi_comm_rank(MPI_COMM_OPA,rang,code)
( 5318)       IF (rang == 0 ) PRINT *,'Taille du buffer alloue : ',taille_tampon
( 5319) 
( 5320)       ! Allocation du tampon et attachement
( 5321)       ALLOCATE(tampon(taille_tampon))
( 5322)       CALL mpi_buffer_attach(tampon,taille_tampon,code)
( 5323) 
( 5324)    END SUBROUTINE mpi_init_opa






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 93

( 5325) 
( 5326) !!DB 2008.04.24
( 5327) !!Routine copied from  /export/home/wangz/NEMO_FOR_DAVE/.../OPA_SRC/IOM/in_out_manager.F90 
( 5328)    SUBROUTINE ctl_stop( cd1, cd2, cd3, cd4, cd5,   &
( 5329)       &                 cd6, cd7, cd8, cd9, cd10 )
( 5330)       !!-----------------------------------------------------------------------
( 5331)       !!                  ***  ROUTINE  stop_opa  ***
( 5332)       !!
( 5333)       !! ** Purpose : ??? blah blah....
( 5334)       !!-----------------------------------------------------------------------
( 5335)       CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
( 5336)       CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
( 5337)       !!-----------------------------------------------------------------------
( 5338)       !
( 5339)       nstop = nstop + 1 
( 5340)       IF(lwp) THEN
( 5341)          WRITE(numout,"(/,' ===>>> : E R R O R',     /,'         ===========',/)") 
( 5342)          IF( PRESENT(cd1 ) ) WRITE(numout,*) cd1
( 5343)          IF( PRESENT(cd2 ) ) WRITE(numout,*) cd2
( 5344)          IF( PRESENT(cd3 ) ) WRITE(numout,*) cd3
( 5345)          IF( PRESENT(cd4 ) ) WRITE(numout,*) cd4
( 5346)          IF( PRESENT(cd5 ) ) WRITE(numout,*) cd5
( 5347)          IF( PRESENT(cd6 ) ) WRITE(numout,*) cd6
( 5348)          IF( PRESENT(cd7 ) ) WRITE(numout,*) cd7
( 5349)          IF( PRESENT(cd8 ) ) WRITE(numout,*) cd8
( 5350)          IF( PRESENT(cd9 ) ) WRITE(numout,*) cd9
( 5351)          IF( PRESENT(cd10) ) WRITE(numout,*) cd10
( 5352)       ENDIF
( 5353)       CALL FLUSH(numout)
( 5354)       !
( 5355)    END SUBROUTINE ctl_stop
( 5356) 
( 5357) 
( 5358) 
( 5359) 
( 5360) #else
( 5361)    !!----------------------------------------------------------------------
( 5362)    !!   Default case:            Dummy module        share memory computing
( 5363)    !!----------------------------------------------------------------------
( 5364)    INTERFACE mpp_sum
( 5365)       MODULE PROCEDURE mpp_sum_a2s, mpp_sum_as, mpp_sum_ai, mpp_sum_s, mpp_sum_i, mpp_sum_c
( 5366)    END INTERFACE
( 5367)    INTERFACE mpp_max
( 5368)       MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
( 5369)    END INTERFACE
( 5370)    INTERFACE mpp_min
( 5371)       MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
( 5372)    END INTERFACE
( 5373)    INTERFACE mpp_isl
( 5374)       MODULE PROCEDURE mppisl_a_int, mppisl_int, mppisl_a_real, mppisl_real
( 5375)    END INTERFACE
( 5376)    INTERFACE mppobc
( 5377)       MODULE PROCEDURE mppobc_1d, mppobc_2d, mppobc_3d, mppobc_4d
( 5378)    END INTERFACE
( 5379)   INTERFACE mpp_minloc
( 5380)      MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
( 5381)   END INTERFACE
( 5382)   INTERFACE mpp_maxloc






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 94

( 5383)      MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
( 5384)   END INTERFACE
( 5385) 
( 5386) 
( 5387)    LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .FALSE.      !: mpp flag
( 5388) 
( 5389) CONTAINS
( 5390) 
( 5391)    FUNCTION mynode(localComm) RESULT (function_value)
( 5392)       INTEGER, OPTIONAL :: localComm
( 5393)       function_value = 0
( 5394)    END FUNCTION mynode
( 5395) 
( 5396)    SUBROUTINE mppsync                       ! Dummy routine
( 5397)    END SUBROUTINE mppsync
( 5398) 
( 5399)    SUBROUTINE mpp_sum_as( parr, kdim )      ! Dummy routine
( 5400)       REAL   , DIMENSION(:) :: parr
( 5401)       INTEGER               :: kdim
( 5402)       WRITE(*,*) 'mpp_sum_as: You should not have seen this print! error?', kdim, parr(1)
( 5403)    END SUBROUTINE mpp_sum_as
( 5404) 
( 5405)    SUBROUTINE mpp_sum_a2s( parr, kdim )      ! Dummy routine
( 5406)       REAL   , DIMENSION(:,:) :: parr
( 5407)       INTEGER               :: kdim
( 5408)       WRITE(*,*) 'mpp_sum_a2s: You should not have seen this print! error?', kdim, parr(1,1)
( 5409)    END SUBROUTINE mpp_sum_a2s
( 5410) 
( 5411)    SUBROUTINE mpp_sum_ai( karr, kdim )      ! Dummy routine
( 5412)       INTEGER, DIMENSION(:) :: karr
( 5413)       INTEGER               :: kdim
( 5414)       WRITE(*,*) 'mpp_sum_ai: You should not have seen this print! error?', kdim, karr(1)
( 5415)    END SUBROUTINE mpp_sum_ai
( 5416) 
( 5417)    SUBROUTINE mpp_sum_s( psca )            ! Dummy routine
( 5418)       REAL                  :: psca
( 5419)       WRITE(*,*) 'mpp_sum_s: You should not have seen this print! error?', psca
( 5420)    END SUBROUTINE mpp_sum_s
( 5421) 
( 5422)    SUBROUTINE mpp_sum_i( kint )            ! Dummy routine
( 5423)       integer               :: kint
( 5424)       WRITE(*,*) 'mpp_sum_i: You should not have seen this print! error?', kint
( 5425)    END SUBROUTINE mpp_sum_i
( 5426) 
( 5427)    SUBROUTINE mpp_sum_c( kcmplx )            ! Dummy routine
( 5428)       complex               :: kcmplx
( 5429)       WRITE(*,*) 'mpp_sum_c: You should not have seen this print! error?', kint
( 5430)    END SUBROUTINE mpp_sum_c
( 5431) 
( 5432)    SUBROUTINE mppmax_a_real( parr, kdim )
( 5433)       REAL   , DIMENSION(:) :: parr
( 5434)       INTEGER               :: kdim
( 5435)       WRITE(*,*) 'mppmax_a_real: You should not have seen this print! error?', kdim, parr(1)
( 5436)    END SUBROUTINE mppmax_a_real
( 5437) 
( 5438)    SUBROUTINE mppmax_real( psca )
( 5439)       REAL                  :: psca
( 5440)       WRITE(*,*) 'mppmax_real: You should not have seen this print! error?', psca






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 95

( 5441)    END SUBROUTINE mppmax_real
( 5442) 
( 5443)    SUBROUTINE mppmin_a_real( parr, kdim )
( 5444)       REAL   , DIMENSION(:) :: parr
( 5445)       INTEGER               :: kdim
( 5446)       WRITE(*,*) 'mppmin_a_real: You should not have seen this print! error?', kdim, parr(1)
( 5447)    END SUBROUTINE mppmin_a_real
( 5448) 
( 5449)    SUBROUTINE mppmin_real( psca )
( 5450)       REAL                  :: psca
( 5451)       WRITE(*,*) 'mppmin_real: You should not have seen this print! error?', psca
( 5452)    END SUBROUTINE mppmin_real
( 5453) 
( 5454)    SUBROUTINE mppmax_a_int( karr, kdim )
( 5455)       INTEGER, DIMENSION(:) :: karr
( 5456)       INTEGER               :: kdim
( 5457)       WRITE(*,*) 'mppmax_a_int: You should not have seen this print! error?', kdim, karr(1)
( 5458)    END SUBROUTINE mppmax_a_int
( 5459) 
( 5460)    SUBROUTINE mppmax_int( kint )
( 5461)       INTEGER               :: kint
( 5462)       WRITE(*,*) 'mppmax_int: You should not have seen this print! error?', kint
( 5463)    END SUBROUTINE mppmax_int
( 5464) 
( 5465)    SUBROUTINE mppmin_a_int( karr, kdim )
( 5466)       INTEGER, DIMENSION(:) :: karr
( 5467)       INTEGER               :: kdim
( 5468)       WRITE(*,*) 'mppmin_a_int: You should not have seen this print! error?', kdim, karr(1)
( 5469)    END SUBROUTINE mppmin_a_int
( 5470) 
( 5471)    SUBROUTINE mppmin_int( kint )
( 5472)       INTEGER               :: kint
( 5473)       WRITE(*,*) 'mppmin_int: You should not have seen this print! error?', kint
( 5474)    END SUBROUTINE mppmin_int
( 5475) 
( 5476)    SUBROUTINE mppobc_1d( parr, kd1, kd2, kl, kk, ktype, kij )
( 5477)     INTEGER  ::   kd1, kd2, kl , kk, ktype, kij
( 5478)     REAL, DIMENSION(:) ::   parr           ! variable array
( 5479)       WRITE(*,*) 'mppobc: You should not have seen this print! error?',   &
( 5480)          &        parr(1), kd1, kd2, kl, kk, ktype, kij
( 5481)    END SUBROUTINE mppobc_1d
( 5482) 
( 5483)    SUBROUTINE mppobc_2d( parr, kd1, kd2, kl, kk, ktype, kij )
( 5484)     INTEGER  ::   kd1, kd2, kl , kk, ktype, kij
( 5485)     REAL, DIMENSION(:,:) ::   parr           ! variable array
( 5486)       WRITE(*,*) 'mppobc: You should not have seen this print! error?',   &
( 5487)          &        parr(1,1), kd1, kd2, kl, kk, ktype, kij
( 5488)    END SUBROUTINE mppobc_2d
( 5489) 
( 5490)    SUBROUTINE mppobc_3d( parr, kd1, kd2, kl, kk, ktype, kij )
( 5491)     INTEGER  ::   kd1, kd2, kl , kk, ktype, kij
( 5492)     REAL, DIMENSION(:,:,:) ::   parr           ! variable array
( 5493)       WRITE(*,*) 'mppobc: You should not have seen this print! error?',   &
( 5494)          &        parr(1,1,1), kd1, kd2, kl, kk, ktype, kij
( 5495)    END SUBROUTINE mppobc_3d
( 5496) 
( 5497)    SUBROUTINE mppobc_4d( parr, kd1, kd2, kl, kk, ktype, kij )
( 5498)     INTEGER  ::   kd1, kd2, kl , kk, ktype, kij






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 96

( 5499)     REAL, DIMENSION(:,:,:,:) ::   parr           ! variable array
( 5500)       WRITE(*,*) 'mppobc: You should not have seen this print! error?',   &
( 5501)          &        parr(1,1,1,1), kd1, kd2, kl, kk, ktype, kij
( 5502)    END SUBROUTINE mppobc_4d
( 5503) 
( 5504) 
( 5505)    SUBROUTINE mpplnks( parr )            ! Dummy routine
( 5506)       REAL, DIMENSION(:,:) :: parr
( 5507)       WRITE(*,*) 'mpplnks: You should not have seen this print! error?', parr(1,1)
( 5508)    END SUBROUTINE mpplnks
( 5509) 
( 5510)    SUBROUTINE mppisl_a_int( karr, kdim )
( 5511)       INTEGER, DIMENSION(:) :: karr
( 5512)       INTEGER               :: kdim
( 5513)       WRITE(*,*) 'mppisl_a_int: You should not have seen this print! error?', kdim, karr(1)
( 5514)    END SUBROUTINE mppisl_a_int
( 5515) 
( 5516)    SUBROUTINE mppisl_int( kint )
( 5517)       INTEGER               :: kint
( 5518)       WRITE(*,*) 'mppisl_int: You should not have seen this print! error?', kint
( 5519)    END SUBROUTINE mppisl_int
( 5520) 
( 5521)    SUBROUTINE mppisl_a_real( parr, kdim )
( 5522)       REAL   , DIMENSION(:) :: parr
( 5523)       INTEGER               :: kdim
( 5524)       WRITE(*,*) 'mppisl_a_real: You should not have seen this print! error?', kdim, parr(1)
( 5525)    END SUBROUTINE mppisl_a_real
( 5526) 
( 5527)    SUBROUTINE mppisl_real( psca )
( 5528)       REAL                  :: psca
( 5529)       WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', psca
( 5530)    END SUBROUTINE mppisl_real
( 5531) 
( 5532)    SUBROUTINE mpp_minloc2d ( ptab, pmask, pmin, ki, kj )
( 5533)       REAL                   :: pmin
( 5534)       REAL , DIMENSION (:,:) :: ptab, pmask
( 5535)       INTEGER :: ki, kj
( 5536)       WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', pmin, ki, kj
( 5537)       WRITE(*,*) '   "      ":             "                 "            ', ptab(1,1), pmask(1,1)
( 5538)    END SUBROUTINE mpp_minloc2d
( 5539) 
( 5540)    SUBROUTINE mpp_minloc3d ( ptab, pmask, pmin, ki, kj, kk )
( 5541)       REAL                     :: pmin
( 5542)       REAL , DIMENSION (:,:,:) :: ptab, pmask
( 5543)       INTEGER :: ki, kj, kk
( 5544)       WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', pmin, ki, kj, kk
( 5545)       WRITE(*,*) '   "      ":             "                 "            ', ptab(1,1,1), pmask(1,1,1)
( 5546)    END SUBROUTINE mpp_minloc3d
( 5547) 
( 5548)    SUBROUTINE mpp_maxloc2d ( ptab, pmask, pmax, ki, kj )
( 5549)       REAL                   :: pmax
( 5550)       REAL , DIMENSION (:,:) :: ptab, pmask
( 5551)       INTEGER :: ki, kj
( 5552)       WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', pmax, ki, kj
( 5553)       WRITE(*,*) '   "      ":             "                 "            ', ptab(1,1), pmask(1,1)
( 5554)    END SUBROUTINE mpp_maxloc2d
( 5555) 
( 5556)    SUBROUTINE mpp_maxloc3d ( ptab, pmask, pmax, ki, kj, kk )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 97

( 5557)       REAL                     :: pmax
( 5558)       REAL , DIMENSION (:,:,:) :: ptab, pmask
( 5559)       INTEGER :: ki, kj, kk
( 5560)       WRITE(*,*) 'mppisl_real: You should not have seen this print! error?', pmax, ki, kj, kk
( 5561)       WRITE(*,*) '   "      ":             "                 "            ', ptab(1,1,1), pmask(1,1,1)
( 5562)    END SUBROUTINE mpp_maxloc3d
( 5563) 
( 5564)    SUBROUTINE mppstop
( 5565)       WRITE(*,*) 'mppstop: You should not have seen this print! error?'
( 5566)    END SUBROUTINE mppstop
( 5567) 
( 5568) #endif
( 5569)    !!----------------------------------------------------------------------
( 5570) END MODULE lib_mpp


















































PGF90 (Version     10.5)          09/23/2010  09:48:58      page 98

(    1) # 1 "lib_mpp.F90"
(    1) MODULE lib_mpp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  lib_mpp  ***
(    4)    !! Ocean numerics:  massively parallel processing librairy
(    5)    !!=====================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_mpp_mpi'     OR      MPI massively parallel processing library
(    9)    !!   'key_mpp_shmem'         SHMEM massively parallel processing library
(   10)    !!----------------------------------------------------------------------
(   11)    !!   mynode
(   12)    !!   mpparent
(   13)    !!   mppshmem
(   14)    !!   mpp_lnk     : generic interface (defined in lbclnk) for :
(   15)    !!                 mpp_lnk_2d, mpp_lnk_3d
(   16)    !!   mpp_lnk_3d_gather :  Message passing manadgement for two 3D arrays
(   17)    !!   mpp_lnk_e   : interface defined in lbclnk
(   18)    !!   mpplnks
(   19)    !!   mpprecv
(   20)    !!   mppsend
(   21)    !!   mppscatter
(   22)    !!   mppgather
(   23)    !!   mpp_isl    : generic inteface  for :
(   24)    !!                mppisl_int , mppisl_a_int , mppisl_real, mppisl_a_real
(   25)    !!   mpp_min    : generic interface for : 
(   26)    !!                mppmin_int , mppmin_a_int , mppmin_real, mppmin_a_real
(   27)    !!   mpp_max    : generic interface for :
(   28)    !!                mppmax_int , mppmax_a_int , mppmax_real, mppmax_a_real
(   29)    !!   mpp_sum    : generic interface for :
(   30)    !!                mppsum_int , mppsum_a_int , mppsum_real, mppsum_a_real
(   31)    !!                mppsum_realdd, mppsum_a_realdd
(   32)    !!   mpp_minloc
(   33)    !!   mpp_maxloc
(   34)    !!   mppsync
(   35)    !!   mppstop
(   36)    !!   mppobc     : variant of mpp_lnk for open boundaries
(   37)    !!   mpp_ini_north
(   38)    !!   mpp_lbc_north
(   39)    !!   mpp_lbc_north_e : variant of mpp_lbc_north for extra outer halo (nsolv=4)
(   40)    !!----------------------------------------------------------------------
(   41)    !! History :
(   42)    !!        !  94 (M. Guyon, J. Escobar, M. Imbard)  Original code
(   43)    !!        !  97  (A.M. Treguier)  SHMEM additions
(   44)    !!        !  98  (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI
(   45)    !!   9.0  !  03  (J.-M. Molines, G. Madec)  F90, free form
(   46)    !!        !  04  (R. Bourdalle Badie)  isend option in mpi
(   47)    !!        !  05  (G. Madec, S. Masson)  npolj=5,6 F-point & ice cases
(   48)    !!        !  05  (R. Redler) Replacement of MPI_COMM_WORLD except for MPI_Abort
(   49)    !!----------------------------------------------------------------------
(   50)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   51)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_mpp.F90,v 1.21 2007/06/05 10:28:55 opalod Exp $ 
(   52)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   53)    !!---------------------------------------------------------------------
(   54)    !! * Modules used
(   55)    USE dom_oce                    ! ocean space and time domain 
(   56)    USE in_out_manager             ! I/O manager
(   57) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 99

(   58)    IMPLICIT NONE
(   59) 
(   60)    PRIVATE
(   61)    PUBLIC  mynode, mpparent, mpp_isl, mpp_min, mpp_max, mpp_sum,  mpp_lbc_north
(   62)    PUBLIC  mpp_lbc_north_e, mpp_minloc, mpp_maxloc, mpp_lnk_3d, mpp_lnk_2d, mpp_lnk_3d_gather, mpp_lnk_2d_e, mpplnks
(   63)    PUBLIC  mpprecv, mppsend, mppscatter, mppgather, mppobc, mpp_ini_north, mppstop, mppsync
(   64) # 67
(   67) !!DB
(   68)    PUBLIC ctl_stop
(   69) 
(   70)    !! * Interfaces
(   71)    !! define generic interface for these routine as they are called sometimes
(   72)    !!        with scalar arguments instead of array arguments, which causes problems
(   73)    !!        for the compilation on AIX system as well as NEC and SGI. Ok on COMPACQ
(   74) 
(   75)    INTERFACE mpp_isl
(   76)       MODULE PROCEDURE mppisl_a_int, mppisl_int, mppisl_a_real, mppisl_real
(   77)    END INTERFACE
(   78)    INTERFACE mpp_min
(   79)       MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
(   80)    END INTERFACE
(   81)    INTERFACE mpp_max
(   82)       MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
(   83)    END INTERFACE
(   84)    INTERFACE mpp_sum
(   85)       MODULE PROCEDURE mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real, &
(   86)                        mppsum_realdd, mppsum_a_realdd
(   87)    END INTERFACE
(   88)    INTERFACE mpp_lbc_north
(   89)       MODULE PROCEDURE mpp_lbc_north_3d, mpp_lbc_north_2d 
(   90)    END INTERFACE
(   91)   INTERFACE mpp_minloc
(   92)      MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
(   93)   END INTERFACE
(   94)   INTERFACE mpp_maxloc
(   95)      MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
(   96)   END INTERFACE
(   97) 
(   98) 
(   99)    !! * Share module variables
(  100)    LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .TRUE.       !: mpp flag
(  101) 
(  102)    !! The processor number is a required power of two : 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,...
(  103)    INTEGER, PARAMETER ::   &
(  104)       nprocmax = 2**10     ! maximun dimension
(  105) 
(  106) # 107
(  107)    !! ========================= !!
(  108)    !!  MPI  variable definition !!
(  109)    !! ========================= !!
(  110) !$AGRIF_DO_NOT_TREAT
(  111) # 1 "/usr/include/mpich2-x86_64/mpif.h"
(    1) !      /* -*- Mode: Fortran; -*- */
(    2) !      
(    3) !      (C) 2001 by Argonne National Laboratory.
(    4) !      See COPYRIGHT in top-level directory.
(    5) !      
(    6) !      DO NOT EDIT






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 100

(    7) !      This file created by buildiface 
(    8) !      
(    9)        INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
(   10)        PARAMETER (MPI_SOURCE=3,MPI_TAG=4,MPI_ERROR=5)
(   11)        INTEGER MPI_STATUS_SIZE
(   12)        PARAMETER (MPI_STATUS_SIZE=5)
(   13)        INTEGER MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
(   14)        INTEGER MPI_STATUSES_IGNORE(MPI_STATUS_SIZE,1)
(   15)        INTEGER MPI_ERRCODES_IGNORE(1)
(   16)        CHARACTER*1 MPI_ARGVS_NULL(1,1)
(   17)        CHARACTER*1 MPI_ARGV_NULL(1)
(   18)        INTEGER MPI_SUCCESS
(   19)        PARAMETER (MPI_SUCCESS=0)
(   20)        INTEGER MPI_ERR_OTHER
(   21)        PARAMETER (MPI_ERR_OTHER=15)
(   22)        INTEGER MPI_ERR_WIN
(   23)        PARAMETER (MPI_ERR_WIN=45)
(   24)        INTEGER MPI_ERR_FILE
(   25)        PARAMETER (MPI_ERR_FILE=27)
(   26)        INTEGER MPI_ERR_COUNT
(   27)        PARAMETER (MPI_ERR_COUNT=2)
(   28)        INTEGER MPI_ERR_SPAWN
(   29)        PARAMETER (MPI_ERR_SPAWN=42)
(   30)        INTEGER MPI_ERR_BASE
(   31)        PARAMETER (MPI_ERR_BASE=46)
(   32)        INTEGER MPI_ERR_RMA_CONFLICT
(   33)        PARAMETER (MPI_ERR_RMA_CONFLICT=49)
(   34)        INTEGER MPI_ERR_IN_STATUS
(   35)        PARAMETER (MPI_ERR_IN_STATUS=17)
(   36)        INTEGER MPI_ERR_INFO_KEY
(   37)        PARAMETER (MPI_ERR_INFO_KEY=29)
(   38)        INTEGER MPI_ERR_LOCKTYPE
(   39)        PARAMETER (MPI_ERR_LOCKTYPE=47)
(   40)        INTEGER MPI_ERR_OP
(   41)        PARAMETER (MPI_ERR_OP=9)
(   42)        INTEGER MPI_ERR_ARG
(   43)        PARAMETER (MPI_ERR_ARG=12)
(   44)        INTEGER MPI_ERR_READ_ONLY
(   45)        PARAMETER (MPI_ERR_READ_ONLY=40)
(   46)        INTEGER MPI_ERR_SIZE
(   47)        PARAMETER (MPI_ERR_SIZE=51)
(   48)        INTEGER MPI_ERR_BUFFER
(   49)        PARAMETER (MPI_ERR_BUFFER=1)
(   50)        INTEGER MPI_ERR_DUP_DATAREP
(   51)        PARAMETER (MPI_ERR_DUP_DATAREP=24)
(   52)        INTEGER MPI_ERR_UNSUPPORTED_DATAREP
(   53)        PARAMETER (MPI_ERR_UNSUPPORTED_DATAREP=43)
(   54)        INTEGER MPI_ERR_LASTCODE
(   55)        PARAMETER (MPI_ERR_LASTCODE=1073741823)
(   56)        INTEGER MPI_ERR_TRUNCATE
(   57)        PARAMETER (MPI_ERR_TRUNCATE=14)
(   58)        INTEGER MPI_ERR_DISP
(   59)        PARAMETER (MPI_ERR_DISP=52)
(   60)        INTEGER MPI_ERR_PORT
(   61)        PARAMETER (MPI_ERR_PORT=38)
(   62)        INTEGER MPI_ERR_INFO_NOKEY
(   63)        PARAMETER (MPI_ERR_INFO_NOKEY=31)
(   64)        INTEGER MPI_ERR_ASSERT






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 101

(   65)        PARAMETER (MPI_ERR_ASSERT=53)
(   66)        INTEGER MPI_ERR_FILE_EXISTS
(   67)        PARAMETER (MPI_ERR_FILE_EXISTS=25)
(   68)        INTEGER MPI_ERR_PENDING
(   69)        PARAMETER (MPI_ERR_PENDING=18)
(   70)        INTEGER MPI_ERR_COMM
(   71)        PARAMETER (MPI_ERR_COMM=5)
(   72)        INTEGER MPI_ERR_KEYVAL
(   73)        PARAMETER (MPI_ERR_KEYVAL=48)
(   74)        INTEGER MPI_ERR_NAME
(   75)        PARAMETER (MPI_ERR_NAME=33)
(   76)        INTEGER MPI_ERR_REQUEST
(   77)        PARAMETER (MPI_ERR_REQUEST=19)
(   78)        INTEGER MPI_ERR_GROUP
(   79)        PARAMETER (MPI_ERR_GROUP=8)
(   80)        INTEGER MPI_ERR_TOPOLOGY
(   81)        PARAMETER (MPI_ERR_TOPOLOGY=10)
(   82)        INTEGER MPI_ERR_TYPE
(   83)        PARAMETER (MPI_ERR_TYPE=3)
(   84)        INTEGER MPI_ERR_TAG
(   85)        PARAMETER (MPI_ERR_TAG=4)
(   86)        INTEGER MPI_ERR_INFO_VALUE
(   87)        PARAMETER (MPI_ERR_INFO_VALUE=30)
(   88)        INTEGER MPI_ERR_NOT_SAME
(   89)        PARAMETER (MPI_ERR_NOT_SAME=35)
(   90)        INTEGER MPI_ERR_RMA_SYNC
(   91)        PARAMETER (MPI_ERR_RMA_SYNC=50)
(   92)        INTEGER MPI_ERR_INFO
(   93)        PARAMETER (MPI_ERR_INFO=28)
(   94)        INTEGER MPI_ERR_NO_MEM
(   95)        PARAMETER (MPI_ERR_NO_MEM=34)
(   96)        INTEGER MPI_ERR_BAD_FILE
(   97)        PARAMETER (MPI_ERR_BAD_FILE=22)
(   98)        INTEGER MPI_ERR_FILE_IN_USE
(   99)        PARAMETER (MPI_ERR_FILE_IN_USE=26)
(  100)        INTEGER MPI_ERR_UNKNOWN
(  101)        PARAMETER (MPI_ERR_UNKNOWN=13)
(  102)        INTEGER MPI_ERR_UNSUPPORTED_OPERATION
(  103)        PARAMETER (MPI_ERR_UNSUPPORTED_OPERATION=44)
(  104)        INTEGER MPI_ERR_QUOTA
(  105)        PARAMETER (MPI_ERR_QUOTA=39)
(  106)        INTEGER MPI_ERR_AMODE
(  107)        PARAMETER (MPI_ERR_AMODE=21)
(  108)        INTEGER MPI_ERR_ROOT
(  109)        PARAMETER (MPI_ERR_ROOT=7)
(  110)        INTEGER MPI_ERR_RANK
(  111)        PARAMETER (MPI_ERR_RANK=6)
(  112)        INTEGER MPI_ERR_DIMS
(  113)        PARAMETER (MPI_ERR_DIMS=11)
(  114)        INTEGER MPI_ERR_NO_SUCH_FILE
(  115)        PARAMETER (MPI_ERR_NO_SUCH_FILE=37)
(  116)        INTEGER MPI_ERR_SERVICE
(  117)        PARAMETER (MPI_ERR_SERVICE=41)
(  118)        INTEGER MPI_ERR_INTERN
(  119)        PARAMETER (MPI_ERR_INTERN=16)
(  120)        INTEGER MPI_ERR_IO
(  121)        PARAMETER (MPI_ERR_IO=32)
(  122)        INTEGER MPI_ERR_ACCESS






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 102

(  123)        PARAMETER (MPI_ERR_ACCESS=20)
(  124)        INTEGER MPI_ERR_NO_SPACE
(  125)        PARAMETER (MPI_ERR_NO_SPACE=36)
(  126)        INTEGER MPI_ERR_CONVERSION
(  127)        PARAMETER (MPI_ERR_CONVERSION=23)
(  128)        INTEGER MPI_ERRORS_ARE_FATAL
(  129)        PARAMETER (MPI_ERRORS_ARE_FATAL=1409286144)
(  130)        INTEGER MPI_ERRORS_RETURN
(  131)        PARAMETER (MPI_ERRORS_RETURN=1409286145)
(  132)        INTEGER MPI_IDENT
(  133)        PARAMETER (MPI_IDENT=0)
(  134)        INTEGER MPI_CONGRUENT
(  135)        PARAMETER (MPI_CONGRUENT=1)
(  136)        INTEGER MPI_SIMILAR
(  137)        PARAMETER (MPI_SIMILAR=2)
(  138)        INTEGER MPI_UNEQUAL
(  139)        PARAMETER (MPI_UNEQUAL=3)
(  140)        INTEGER MPI_MAX
(  141)        PARAMETER (MPI_MAX=1476395009)
(  142)        INTEGER MPI_MIN
(  143)        PARAMETER (MPI_MIN=1476395010)
(  144)        INTEGER MPI_SUM
(  145)        PARAMETER (MPI_SUM=1476395011)
(  146)        INTEGER MPI_PROD
(  147)        PARAMETER (MPI_PROD=1476395012)
(  148)        INTEGER MPI_LAND
(  149)        PARAMETER (MPI_LAND=1476395013)
(  150)        INTEGER MPI_BAND
(  151)        PARAMETER (MPI_BAND=1476395014)
(  152)        INTEGER MPI_LOR
(  153)        PARAMETER (MPI_LOR=1476395015)
(  154)        INTEGER MPI_BOR
(  155)        PARAMETER (MPI_BOR=1476395016)
(  156)        INTEGER MPI_LXOR
(  157)        PARAMETER (MPI_LXOR=1476395017)
(  158)        INTEGER MPI_BXOR
(  159)        PARAMETER (MPI_BXOR=1476395018)
(  160)        INTEGER MPI_MINLOC
(  161)        PARAMETER (MPI_MINLOC=1476395019)
(  162)        INTEGER MPI_MAXLOC
(  163)        PARAMETER (MPI_MAXLOC=1476395020)
(  164)        INTEGER MPI_REPLACE
(  165)        PARAMETER (MPI_REPLACE=1476395021)
(  166)        INTEGER MPI_COMM_WORLD
(  167)        PARAMETER (MPI_COMM_WORLD=1140850688)
(  168)        INTEGER MPI_COMM_SELF
(  169)        PARAMETER (MPI_COMM_SELF=1140850689)
(  170)        INTEGER MPI_GROUP_EMPTY
(  171)        PARAMETER (MPI_GROUP_EMPTY=1207959552)
(  172)        INTEGER MPI_COMM_NULL
(  173)        PARAMETER (MPI_COMM_NULL=67108864)
(  174)        INTEGER MPI_WIN_NULL
(  175)        PARAMETER (MPI_WIN_NULL=536870912)
(  176)        INTEGER MPI_FILE_NULL
(  177)        PARAMETER (MPI_FILE_NULL=0)
(  178)        INTEGER MPI_GROUP_NULL
(  179)        PARAMETER (MPI_GROUP_NULL=134217728)
(  180)        INTEGER MPI_OP_NULL






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 103

(  181)        PARAMETER (MPI_OP_NULL=402653184)
(  182)        INTEGER MPI_DATATYPE_NULL
(  183)        PARAMETER (MPI_DATATYPE_NULL=201326592)
(  184)        INTEGER MPI_REQUEST_NULL
(  185)        PARAMETER (MPI_REQUEST_NULL=738197504)
(  186)        INTEGER MPI_ERRHANDLER_NULL
(  187)        PARAMETER (MPI_ERRHANDLER_NULL=335544320)
(  188)        INTEGER MPI_INFO_NULL
(  189)        PARAMETER (MPI_INFO_NULL=469762048)
(  190)        INTEGER MPI_TAG_UB
(  191)        PARAMETER (MPI_TAG_UB=1681915906)
(  192)        INTEGER MPI_HOST
(  193)        PARAMETER (MPI_HOST=1681915908)
(  194)        INTEGER MPI_IO
(  195)        PARAMETER (MPI_IO=1681915910)
(  196)        INTEGER MPI_WTIME_IS_GLOBAL
(  197)        PARAMETER (MPI_WTIME_IS_GLOBAL=1681915912)
(  198)        INTEGER MPI_UNIVERSE_SIZE
(  199)        PARAMETER (MPI_UNIVERSE_SIZE=1681915914)
(  200)        INTEGER MPI_LASTUSEDCODE
(  201)        PARAMETER (MPI_LASTUSEDCODE=1681915916)
(  202)        INTEGER MPI_APPNUM
(  203)        PARAMETER (MPI_APPNUM=1681915918)
(  204)        INTEGER MPI_WIN_BASE
(  205)        PARAMETER (MPI_WIN_BASE=1711276034)
(  206)        INTEGER MPI_WIN_SIZE
(  207)        PARAMETER (MPI_WIN_SIZE=1711276036)
(  208)        INTEGER MPI_WIN_DISP_UNIT
(  209)        PARAMETER (MPI_WIN_DISP_UNIT=1711276038)
(  210)        INTEGER MPI_MAX_ERROR_STRING
(  211)        PARAMETER (MPI_MAX_ERROR_STRING=1023)
(  212)        INTEGER MPI_MAX_PORT_NAME
(  213)        PARAMETER (MPI_MAX_PORT_NAME=255)
(  214)        INTEGER MPI_MAX_OBJECT_NAME
(  215)        PARAMETER (MPI_MAX_OBJECT_NAME=127)
(  216)        INTEGER MPI_MAX_INFO_KEY
(  217)        PARAMETER (MPI_MAX_INFO_KEY=254)
(  218)        INTEGER MPI_MAX_INFO_VAL
(  219)        PARAMETER (MPI_MAX_INFO_VAL=1023)
(  220)        INTEGER MPI_MAX_PROCESSOR_NAME
(  221)        PARAMETER (MPI_MAX_PROCESSOR_NAME=128-1)
(  222)        INTEGER MPI_MAX_DATAREP_STRING
(  223)        PARAMETER (MPI_MAX_DATAREP_STRING=127)
(  224)        INTEGER MPI_UNDEFINED
(  225)        PARAMETER (MPI_UNDEFINED=(-32766))
(  226)        INTEGER MPI_KEYVAL_INVALID
(  227)        PARAMETER (MPI_KEYVAL_INVALID=603979776)
(  228)        INTEGER MPI_BSEND_OVERHEAD
(  229)        PARAMETER (MPI_BSEND_OVERHEAD=88)
(  230)        INTEGER MPI_PROC_NULL
(  231)        PARAMETER (MPI_PROC_NULL=-1)
(  232)        INTEGER MPI_ANY_SOURCE
(  233)        PARAMETER (MPI_ANY_SOURCE=-2)
(  234)        INTEGER MPI_ANY_TAG
(  235)        PARAMETER (MPI_ANY_TAG=-1)
(  236)        INTEGER MPI_ROOT
(  237)        PARAMETER (MPI_ROOT=-3)
(  238)        INTEGER MPI_GRAPH






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 104

(  239)        PARAMETER (MPI_GRAPH=1)
(  240)        INTEGER MPI_CART
(  241)        PARAMETER (MPI_CART=2)
(  242)        INTEGER MPI_VERSION
(  243)        PARAMETER (MPI_VERSION=2)
(  244)        INTEGER MPI_SUBVERSION
(  245)        PARAMETER (MPI_SUBVERSION=2)
(  246)        INTEGER MPI_LOCK_EXCLUSIVE
(  247)        PARAMETER (MPI_LOCK_EXCLUSIVE=234)
(  248)        INTEGER MPI_LOCK_SHARED
(  249)        PARAMETER (MPI_LOCK_SHARED=235)
(  250)        INTEGER MPI_COMPLEX
(  251)        PARAMETER (MPI_COMPLEX=1275070494)
(  252)        INTEGER MPI_DOUBLE_COMPLEX
(  253)        PARAMETER (MPI_DOUBLE_COMPLEX=1275072546)
(  254)        INTEGER MPI_LOGICAL
(  255)        PARAMETER (MPI_LOGICAL=1275069469)
(  256)        INTEGER MPI_REAL
(  257)        PARAMETER (MPI_REAL=1275069468)
(  258)        INTEGER MPI_DOUBLE_PRECISION
(  259)        PARAMETER (MPI_DOUBLE_PRECISION=1275070495)
(  260)        INTEGER MPI_INTEGER
(  261)        PARAMETER (MPI_INTEGER=1275069467)
(  262)        INTEGER MPI_2INTEGER
(  263)        PARAMETER (MPI_2INTEGER=1275070496)
(  264)        INTEGER MPI_2COMPLEX
(  265)        PARAMETER (MPI_2COMPLEX=1275072548)
(  266)        INTEGER MPI_2DOUBLE_PRECISION
(  267)        PARAMETER (MPI_2DOUBLE_PRECISION=1275072547)
(  268)        INTEGER MPI_2REAL
(  269)        PARAMETER (MPI_2REAL=1275070497)
(  270)        INTEGER MPI_2DOUBLE_COMPLEX
(  271)        PARAMETER (MPI_2DOUBLE_COMPLEX=1275076645)
(  272)        INTEGER MPI_CHARACTER
(  273)        PARAMETER (MPI_CHARACTER=1275068698)
(  274)        INTEGER MPI_BYTE
(  275)        PARAMETER (MPI_BYTE=1275068685)
(  276)        INTEGER MPI_UB
(  277)        PARAMETER (MPI_UB=1275068433)
(  278)        INTEGER MPI_LB
(  279)        PARAMETER (MPI_LB=1275068432)
(  280)        INTEGER MPI_PACKED
(  281)        PARAMETER (MPI_PACKED=1275068687)
(  282)        INTEGER MPI_INTEGER1
(  283)        PARAMETER (MPI_INTEGER1=1275068717)
(  284)        INTEGER MPI_INTEGER2
(  285)        PARAMETER (MPI_INTEGER2=1275068975)
(  286)        INTEGER MPI_INTEGER4
(  287)        PARAMETER (MPI_INTEGER4=1275069488)
(  288)        INTEGER MPI_INTEGER8
(  289)        PARAMETER (MPI_INTEGER8=1275070513)
(  290)        INTEGER MPI_INTEGER16
(  291)        PARAMETER (MPI_INTEGER16=MPI_DATATYPE_NULL)
(  292)        INTEGER MPI_REAL4
(  293)        PARAMETER (MPI_REAL4=1275069479)
(  294)        INTEGER MPI_REAL8
(  295)        PARAMETER (MPI_REAL8=1275070505)
(  296)        INTEGER MPI_REAL16






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 105

(  297)        PARAMETER (MPI_REAL16=1275072555)
(  298)        INTEGER MPI_COMPLEX8
(  299)        PARAMETER (MPI_COMPLEX8=1275070504)
(  300)        INTEGER MPI_COMPLEX16
(  301)        PARAMETER (MPI_COMPLEX16=1275072554)
(  302)        INTEGER MPI_COMPLEX32
(  303)        PARAMETER (MPI_COMPLEX32=1275076652)
(  304)        INTEGER MPI_ADDRESS_KIND, MPI_OFFSET_KIND
(  305)        PARAMETER (MPI_ADDRESS_KIND=8)
(  306)        PARAMETER (MPI_OFFSET_KIND=8)
(  307)        INTEGER MPI_CHAR
(  308)        PARAMETER (MPI_CHAR=1275068673)
(  309)        INTEGER MPI_SIGNED_CHAR
(  310)        PARAMETER (MPI_SIGNED_CHAR=1275068696)
(  311)        INTEGER MPI_UNSIGNED_CHAR
(  312)        PARAMETER (MPI_UNSIGNED_CHAR=1275068674)
(  313)        INTEGER MPI_WCHAR
(  314)        PARAMETER (MPI_WCHAR=1275069454)
(  315)        INTEGER MPI_SHORT
(  316)        PARAMETER (MPI_SHORT=1275068931)
(  317)        INTEGER MPI_UNSIGNED_SHORT
(  318)        PARAMETER (MPI_UNSIGNED_SHORT=1275068932)
(  319)        INTEGER MPI_INT
(  320)        PARAMETER (MPI_INT=1275069445)
(  321)        INTEGER MPI_UNSIGNED
(  322)        PARAMETER (MPI_UNSIGNED=1275069446)
(  323)        INTEGER MPI_LONG
(  324)        PARAMETER (MPI_LONG=1275070471)
(  325)        INTEGER MPI_UNSIGNED_LONG
(  326)        PARAMETER (MPI_UNSIGNED_LONG=1275070472)
(  327)        INTEGER MPI_FLOAT
(  328)        PARAMETER (MPI_FLOAT=1275069450)
(  329)        INTEGER MPI_DOUBLE
(  330)        PARAMETER (MPI_DOUBLE=1275070475)
(  331)        INTEGER MPI_LONG_DOUBLE
(  332)        PARAMETER (MPI_LONG_DOUBLE=1275072524)
(  333)        INTEGER MPI_LONG_LONG_INT
(  334)        PARAMETER (MPI_LONG_LONG_INT=1275070473)
(  335)        INTEGER MPI_UNSIGNED_LONG_LONG
(  336)        PARAMETER (MPI_UNSIGNED_LONG_LONG=1275070489)
(  337)        INTEGER MPI_LONG_LONG
(  338)        PARAMETER (MPI_LONG_LONG=1275070473)
(  339)        INTEGER MPI_FLOAT_INT
(  340)        PARAMETER (MPI_FLOAT_INT=-1946157056)
(  341)        INTEGER MPI_DOUBLE_INT
(  342)        PARAMETER (MPI_DOUBLE_INT=-1946157055)
(  343)        INTEGER MPI_LONG_INT
(  344)        PARAMETER (MPI_LONG_INT=-1946157054)
(  345)        INTEGER MPI_SHORT_INT
(  346)        PARAMETER (MPI_SHORT_INT=-1946157053)
(  347)        INTEGER MPI_2INT
(  348)        PARAMETER (MPI_2INT=1275070486)
(  349)        INTEGER MPI_LONG_DOUBLE_INT
(  350)        PARAMETER (MPI_LONG_DOUBLE_INT=-1946157052)
(  351)        INTEGER MPI_INT8_T
(  352)        PARAMETER (MPI_INT8_T=1275068727)
(  353)        INTEGER MPI_INT16_T
(  354)        PARAMETER (MPI_INT16_T=1275068984)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 106

(  355)        INTEGER MPI_INT32_T
(  356)        PARAMETER (MPI_INT32_T=1275069497)
(  357)        INTEGER MPI_INT64_T
(  358)        PARAMETER (MPI_INT64_T=1275070522)
(  359)        INTEGER MPI_UINT8_T
(  360)        PARAMETER (MPI_UINT8_T=1275068731)
(  361)        INTEGER MPI_UINT16_T
(  362)        PARAMETER (MPI_UINT16_T=1275068988)
(  363)        INTEGER MPI_UINT32_T
(  364)        PARAMETER (MPI_UINT32_T=1275069501)
(  365)        INTEGER MPI_UINT64_T
(  366)        PARAMETER (MPI_UINT64_T=1275070526)
(  367)        INTEGER MPI_C_BOOL
(  368)        PARAMETER (MPI_C_BOOL=1275068735)
(  369)        INTEGER MPI_C_FLOAT_COMPLEX
(  370)        PARAMETER (MPI_C_FLOAT_COMPLEX=1275070528)
(  371)        INTEGER MPI_C_COMPLEX
(  372)        PARAMETER (MPI_C_COMPLEX=1275070528)
(  373)        INTEGER MPI_C_DOUBLE_COMPLEX
(  374)        PARAMETER (MPI_C_DOUBLE_COMPLEX=1275072577)
(  375)        INTEGER MPI_C_LONG_DOUBLE_COMPLEX
(  376)        PARAMETER (MPI_C_LONG_DOUBLE_COMPLEX=1275076674)
(  377)        INTEGER MPI_AINT
(  378)        PARAMETER (MPI_AINT=1275070531)
(  379)        INTEGER MPI_OFFSET
(  380)        PARAMETER (MPI_OFFSET=1275070532)
(  381)        INTEGER MPI_COMBINER_NAMED
(  382)        PARAMETER (MPI_COMBINER_NAMED=1)
(  383)        INTEGER MPI_COMBINER_DUP
(  384)        PARAMETER (MPI_COMBINER_DUP=2)
(  385)        INTEGER MPI_COMBINER_CONTIGUOUS
(  386)        PARAMETER (MPI_COMBINER_CONTIGUOUS=3)
(  387)        INTEGER MPI_COMBINER_VECTOR
(  388)        PARAMETER (MPI_COMBINER_VECTOR=4)
(  389)        INTEGER MPI_COMBINER_HVECTOR_INTEGER
(  390)        PARAMETER (MPI_COMBINER_HVECTOR_INTEGER=5)
(  391)        INTEGER MPI_COMBINER_HVECTOR
(  392)        PARAMETER (MPI_COMBINER_HVECTOR=6)
(  393)        INTEGER MPI_COMBINER_INDEXED
(  394)        PARAMETER (MPI_COMBINER_INDEXED=7)
(  395)        INTEGER MPI_COMBINER_HINDEXED_INTEGER
(  396)        PARAMETER (MPI_COMBINER_HINDEXED_INTEGER=8)
(  397)        INTEGER MPI_COMBINER_HINDEXED
(  398)        PARAMETER (MPI_COMBINER_HINDEXED=9)
(  399)        INTEGER MPI_COMBINER_INDEXED_BLOCK
(  400)        PARAMETER (MPI_COMBINER_INDEXED_BLOCK=10)
(  401)        INTEGER MPI_COMBINER_STRUCT_INTEGER
(  402)        PARAMETER (MPI_COMBINER_STRUCT_INTEGER=11)
(  403)        INTEGER MPI_COMBINER_STRUCT
(  404)        PARAMETER (MPI_COMBINER_STRUCT=12)
(  405)        INTEGER MPI_COMBINER_SUBARRAY
(  406)        PARAMETER (MPI_COMBINER_SUBARRAY=13)
(  407)        INTEGER MPI_COMBINER_DARRAY
(  408)        PARAMETER (MPI_COMBINER_DARRAY=14)
(  409)        INTEGER MPI_COMBINER_F90_REAL
(  410)        PARAMETER (MPI_COMBINER_F90_REAL=15)
(  411)        INTEGER MPI_COMBINER_F90_COMPLEX
(  412)        PARAMETER (MPI_COMBINER_F90_COMPLEX=16)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 107

(  413)        INTEGER MPI_COMBINER_F90_INTEGER
(  414)        PARAMETER (MPI_COMBINER_F90_INTEGER=17)
(  415)        INTEGER MPI_COMBINER_RESIZED
(  416)        PARAMETER (MPI_COMBINER_RESIZED=18)
(  417)        INTEGER MPI_TYPECLASS_REAL
(  418)        PARAMETER (MPI_TYPECLASS_REAL=1)
(  419)        INTEGER MPI_TYPECLASS_INTEGER
(  420)        PARAMETER (MPI_TYPECLASS_INTEGER=2)
(  421)        INTEGER MPI_TYPECLASS_COMPLEX
(  422)        PARAMETER (MPI_TYPECLASS_COMPLEX=3)
(  423)        INTEGER MPI_MODE_NOCHECK
(  424)        PARAMETER (MPI_MODE_NOCHECK=1024)
(  425)        INTEGER MPI_MODE_NOSTORE
(  426)        PARAMETER (MPI_MODE_NOSTORE=2048)
(  427)        INTEGER MPI_MODE_NOPUT
(  428)        PARAMETER (MPI_MODE_NOPUT=4096)
(  429)        INTEGER MPI_MODE_NOPRECEDE
(  430)        PARAMETER (MPI_MODE_NOPRECEDE=8192)
(  431)        INTEGER MPI_MODE_NOSUCCEED
(  432)        PARAMETER (MPI_MODE_NOSUCCEED=16384)
(  433)        INTEGER MPI_THREAD_SINGLE
(  434)        PARAMETER (MPI_THREAD_SINGLE=0)
(  435)        INTEGER MPI_THREAD_FUNNELED
(  436)        PARAMETER (MPI_THREAD_FUNNELED=1)
(  437)        INTEGER MPI_THREAD_SERIALIZED
(  438)        PARAMETER (MPI_THREAD_SERIALIZED=2)
(  439)        INTEGER MPI_THREAD_MULTIPLE
(  440)        PARAMETER (MPI_THREAD_MULTIPLE=3)
(  441)        INTEGER MPI_MODE_RDONLY
(  442)        PARAMETER (MPI_MODE_RDONLY=2)
(  443)        INTEGER MPI_MODE_RDWR
(  444)        PARAMETER (MPI_MODE_RDWR=8)
(  445)        INTEGER MPI_MODE_WRONLY
(  446)        PARAMETER (MPI_MODE_WRONLY=4)
(  447)        INTEGER MPI_MODE_DELETE_ON_CLOSE
(  448)        PARAMETER (MPI_MODE_DELETE_ON_CLOSE=16)
(  449)        INTEGER MPI_MODE_UNIQUE_OPEN
(  450)        PARAMETER (MPI_MODE_UNIQUE_OPEN=32)
(  451)        INTEGER MPI_MODE_CREATE
(  452)        PARAMETER (MPI_MODE_CREATE=1)
(  453)        INTEGER MPI_MODE_EXCL
(  454)        PARAMETER (MPI_MODE_EXCL=64)
(  455)        INTEGER MPI_MODE_APPEND
(  456)        PARAMETER (MPI_MODE_APPEND=128)
(  457)        INTEGER MPI_MODE_SEQUENTIAL
(  458)        PARAMETER (MPI_MODE_SEQUENTIAL=256)
(  459)        INTEGER MPI_SEEK_SET
(  460)        PARAMETER (MPI_SEEK_SET=600)
(  461)        INTEGER MPI_SEEK_CUR
(  462)        PARAMETER (MPI_SEEK_CUR=602)
(  463)        INTEGER MPI_SEEK_END
(  464)        PARAMETER (MPI_SEEK_END=604)
(  465)        INTEGER MPI_ORDER_C
(  466)        PARAMETER (MPI_ORDER_C=56)
(  467)        INTEGER MPI_ORDER_FORTRAN
(  468)        PARAMETER (MPI_ORDER_FORTRAN=57)
(  469)        INTEGER MPI_DISTRIBUTE_BLOCK
(  470)        PARAMETER (MPI_DISTRIBUTE_BLOCK=121)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 108

(  471)        INTEGER MPI_DISTRIBUTE_CYCLIC
(  472)        PARAMETER (MPI_DISTRIBUTE_CYCLIC=122)
(  473)        INTEGER MPI_DISTRIBUTE_NONE
(  474)        PARAMETER (MPI_DISTRIBUTE_NONE=123)
(  475)        INTEGER MPI_DISTRIBUTE_DFLT_DARG
(  476)        PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
(  477)        integer*8 MPI_DISPLACEMENT_CURRENT
(  478)        PARAMETER (MPI_DISPLACEMENT_CURRENT=-54278278)
(  479)        INTEGER MPI_BOTTOM, MPI_IN_PLACE
(  480)        EXTERNAL MPI_DUP_FN, MPI_NULL_DELETE_FN, MPI_NULL_COPY_FN
(  481)        EXTERNAL MPI_WTIME, MPI_WTICK
(  482)        EXTERNAL PMPI_WTIME, PMPI_WTICK
(  483)        EXTERNAL MPI_COMM_DUP_FN, MPI_COMM_NULL_DELETE_FN
(  484)        EXTERNAL MPI_COMM_NULL_COPY_FN
(  485)        EXTERNAL MPI_WIN_DUP_FN, MPI_WIN_NULL_DELETE_FN
(  486)        EXTERNAL MPI_WIN_NULL_COPY_FN
(  487)        EXTERNAL MPI_TYPE_DUP_FN, MPI_TYPE_NULL_DELETE_FN
(  488)        EXTERNAL MPI_TYPE_NULL_COPY_FN
(  489)        EXTERNAL MPI_CONVERSION_FN_NULL
(  490)        REAL*8 MPI_WTIME, MPI_WTICK
(  491)        REAL*8 PMPI_WTIME, PMPI_WTICK
(  492) 
(  493) 
(  494)        COMMON /MPIPRIV1/ MPI_BOTTOM, MPI_IN_PLACE, MPI_STATUS_IGNORE
(  495) 
(  496)        COMMON /MPIPRIV2/ MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE
(  497)        SAVE /MPIPRIV1/,/MPIPRIV2/
(  498) 
(  499)        COMMON /MPIPRIVC/ MPI_ARGVS_NULL, MPI_ARGV_NULL
(  500)        SAVE   /MPIPRIVC/
(  501) # 112 "lib_mpp.F90"
(  112) # 112
(  112) !$AGRIF_END_DO_NOT_TREAT
(  113) 
(  114)    INTEGER ::   &
(  115)       mppsize,  &  ! number of process
(  116)       mpprank,  &  ! process number  [ 0 - size-1 ]
(  117)       mpi_comm_opa ! opa local communicator
(  118) 
(  119)    ! MD: needed to define MPI_SUMDD as a reduction operator
(  120)    INTEGER, PUBLIC :: MPI_SUMDD
(  121) 
(  122)    ! variables used in case of north fold condition in mpp_mpi with jpni > 1
(  123)    INTEGER ::      &       !
(  124)       ngrp_world,  &       ! group ID for the world processors
(  125)       ngrp_north,  &       ! group ID for the northern processors (to be fold)
(  126)       ncomm_north, &       ! communicator made by the processors belonging to ngrp_north
(  127)       ndim_rank_north, &   ! number of 'sea' processor in the northern line (can be /= jpni !)
(  128)       njmppmax             ! value of njmpp for the processors of the northern line
(  129)    INTEGER ::      &       !
(  130)       north_root           ! number (in the comm_opa) of proc 0 in the northern comm
(  131)    INTEGER, DIMENSION(:), ALLOCATABLE ::   &
(  132)       nrank_north          ! dimension ndim_rank_north, number of the procs belonging to ncomm_north
(  133)    CHARACTER (len=1) ::  &
(  134)       c_mpi_send = 'S'     ! type od mpi send/recieve (S=standard, B=bsend, I=isend)
(  135)    LOGICAL  ::           &
(  136)       l_isend = .FALSE.    ! isend use indicator (T if c_mpi_send='I')
(  137) !!DB: 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 109

(  138)    CHARACTER(len=200) ::   ctmp1, ctmp2, ctmp3      !: temporary character
(  139) 
(  140) 
(  141) 
(  142) # 257
(  257)    REAL(wp), DIMENSION(jpi,jprecj,jpk,2,2) ::   &
(  258)        t4ns, t4sn  ! 3d message passing arrays north-south & south-north
(  259)    REAL(wp), DIMENSION(jpj,jpreci,jpk,2,2) ::   &
(  260)        t4ew, t4we  ! 3d message passing arrays east-west & west-east
(  261)    REAL(wp), DIMENSION(jpi,jprecj,jpk,2,2) ::   &
(  262)        t4p1, t4p2  ! 3d message passing arrays north fold
(  263)    REAL(wp), DIMENSION(jpi,jprecj,jpk,2) ::   &
(  264)        t3ns, t3sn  ! 3d message passing arrays north-south & south-north
(  265)    REAL(wp), DIMENSION(jpj,jpreci,jpk,2) ::   &
(  266)        t3ew, t3we  ! 3d message passing arrays east-west & west-east
(  267)    REAL(wp), DIMENSION(jpi,jprecj,jpk,2) ::   &
(  268)        t3p1, t3p2  ! 3d message passing arrays north fold
(  269)    REAL(wp), DIMENSION(jpi,jprecj,2) ::   &
(  270)        t2ns, t2sn  ! 2d message passing arrays north-south & south-north
(  271)    REAL(wp), DIMENSION(jpj,jpreci,2) ::   &
(  272)        t2ew, t2we  ! 2d message passing arrays east-west & west-east
(  273)    REAL(wp), DIMENSION(jpi,jprecj,2) ::   &
(  274)        t2p1, t2p2  ! 2d message passing arrays north fold
(  275)    REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,jprecj+jpr2dj,2) ::   &
(  276)        tr2ns, tr2sn  ! 2d message passing arrays north-south & south-north including extra outer halo
(  277)    REAL(wp), DIMENSION(1-jpr2dj:jpj+jpr2dj,jpreci+jpr2di,2) ::   &
(  278)        tr2ew, tr2we  ! 2d message passing arrays east-west & west-east including extra outer halo
(  279)    !!----------------------------------------------------------------------
(  280)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(  281)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_mpp.F90,v 1.20 2007/03/07 13:26:44 opalod Exp $ 
(  282)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(  283)    !!---------------------------------------------------------------------
(  284) 
(  285) CONTAINS
(  286)    !!---------------------------------------------------------------------
(  287)    !!   Routine DDPDD_MPI: used by reduction operator MPI_SUMDD
(  288)    !!
(  289)    !!   Modification of original codes written by David H. Bailey
(  290)    !!   This subroutine computes ddb(i) = dda(i)+ddb(i)
(  291)    !!---------------------------------------------------------------------
(  292)       subroutine DDPDD_MPI (dda, ddb, len, itype)
(  293)       implicit none
(  294) 
(  295)     !! * Arguments
(  296)       INTEGER, INTENT(in)                         :: len, itype
(  297)       COMPLEX, DIMENSION(len), INTENT(in)     :: dda
(  298)       COMPLEX, DIMENSION(len), INTENT(inout)  :: ddb
(  299) 
(  300)     !! * Local variables
(  301)       REAL :: e, t1, t2  ! local work variables
(  302)       INTEGER :: i           ! local scalar
(  303) 
(  304)       DO i=1,len
(  305)    ! Compute dda + ddb using Knuth's trick.
(  306)         t1 = real(dda(i)) + real(ddb(i))
(  307)         e = t1 - real(dda(i))
(  308)         t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e))) &
(  309)              +imag(dda(i)) + imag(ddb(i))






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 110

(  310) 
(  311)    ! The result is t1 + t2, after normalization.
(  312)         ddb(i) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
(  313)       END DO
(  314)       END SUBROUTINE DDPDD_MPI
(  315)    !!---------------------------------------------------------------------
(  316) 
(  317) 
(  318)    FUNCTION mynode(localComm)
(  319)       !!----------------------------------------------------------------------
(  320)       !!                  ***  routine mynode  ***
(  321)       !!                    
(  322)       !! ** Purpose :   Find processor unit
(  323)       !!
(  324)       !!----------------------------------------------------------------------
(  325) # 326
(  326)       !! * Local variables   (MPI version)
(  327)       INTEGER ::   mynode, ierr, code
(  328)       LOGICAL ::   mpi_was_called
(  329)       INTEGER,OPTIONAL ::   localComm
(  330)       NAMELIST/nam_mpp/ c_mpi_send
(  331)       !!----------------------------------------------------------------------
(  332) 
(  333) !      WRITE(numout,*)
(  334) !      WRITE(numout,*) 'mynode : mpi initialisation'
(  335) !      WRITE(numout,*) '~~~~~~ '
(  336) !      WRITE(numout,*)
(  337) 
(  338)       ! Namelist namrun : parameters of the run
(  339) !      REWIND( numnam )
(  340) !      READ  ( numnam, nam_mpp )
(  341) 
(  342) !      WRITE(numout,*) '        Namelist nam_mpp'
(  343) !      WRITE(numout,*) '           mpi send type            c_mpi_send = ', c_mpi_send
(  344) 
(  345) ! MD
(  346) ! June 2006
(  347) ! MPI_INIT must be called before any files are opened!
(  348) ! So it's not possible to determine c_mpi_send from namelist.
(  349) ! Hence this is hacked to work with a hardcoded mpi send type:
(  350) 
(  351) c_mpi_send = 'I'
(  352) !c_mpi_send = 'I'
(  353) !c_mpi_send = 'B'
(  354) 
(  355) # 358
(  358) !!bug RB : should be clean to use Agrif in coupled mode
(  359) # 360
(  360)          CALL mpi_initialized ( mpi_was_called, code )
(  361)          IF( code /= MPI_SUCCESS ) THEN
(  362)             CALL ctl_stop( ' lib_mpp: Error in routine mpi_initialized' )
(  363)             CALL mpi_abort( mpi_comm_world, code, ierr )
(  364)          ENDIF
(  365) 
(  366)          IF( PRESENT(localComm) .and. mpi_was_called ) THEN
(  367)             mpi_comm_opa = localComm
(  368)             SELECT CASE ( c_mpi_send )
(  369)             CASE ( 'S' )                ! Standard mpi send (blocking)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 111

(  370) !               WRITE(numout,*) '           Standard blocking mpi send (send)'
(  371)             CASE ( 'B' )                ! Buffer mpi send (blocking)
(  372) !               WRITE(numout,*) '           Buffer blocking mpi send (bsend)'
(  373)                CALL mpi_init_opa( ierr )
(  374)             CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
(  375) !               WRITE(numout,*) '           Immediate non-blocking send (isend)'
(  376)                l_isend = .TRUE.
(  377)             CASE DEFAULT
(  378) !               WRITE(numout,cform_err)
(  379) !               WRITE(numout,*) '           bad value for c_mpi_send = ', c_mpi_send
(  380)                nstop = nstop + 1
(  381)             END SELECT
(  382)          ELSE IF ( PRESENT(localComm) .and. .not. mpi_was_called ) THEN
(  383) !            WRITE(numout,*) ' lib_mpp: You cannot provide a local communicator '
(  384) !            WRITE(numout,*) '          without calling MPI_Init before ! '
(  385)          ELSE
(  386) # 387
(  387)             SELECT CASE ( c_mpi_send )
(  388)             CASE ( 'S' )                ! Standard mpi send (blocking)
(  389) !               WRITE(numout,*) '           Standard blocking mpi send (send)'
(  390)                CALL mpi_init( ierr )
(  391)             CASE ( 'B' )                ! Buffer mpi send (blocking)
(  392) !               WRITE(numout,*) '           Buffer blocking mpi send (bsend)'
(  393)                CALL mpi_init_opa( ierr )
(  394)             CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
(  395) !               WRITE(numout,*) '           Immediate non-blocking send (isend)'
(  396)                l_isend = .TRUE.
(  397)                CALL mpi_init( ierr )
(  398)             CASE DEFAULT
(  399) !               WRITE(ctmp1,*) '           bad value for c_mpi_send = ', c_mpi_send
(  400)                CALL ctl_stop( ctmp1 )
(  401)             END SELECT
(  402) 
(  403) # 404
(  404)             CALL mpi_comm_dup( mpi_comm_world, mpi_comm_opa, code)
(  405)             IF( code /= MPI_SUCCESS ) THEN
(  406)                CALL ctl_stop( ' lib_mpp: Error in routine mpi_comm_dup' )
(  407)                CALL mpi_abort( mpi_comm_world, code, ierr )
(  408)             ENDIF
(  409)             !
(  410)          ENDIF
(  411) 
(  412) # 432
(  432)       CALL mpi_comm_rank( mpi_comm_opa, mpprank, ierr )
(  433)       CALL mpi_comm_size( mpi_comm_opa, mppsize, ierr )
(  434)       mynode = mpprank
(  435) 
(  436) !MD: Create MPISUMDD
(  437)       call MPI_OP_CREATE(DDPDD_MPI, .TRUE., MPI_SUMDD, ierr)
(  438) 
(  439) # 555
(  555)    END FUNCTION mynode
(  556) 
(  557) 
(  558)    SUBROUTINE mpparent( kparent_tid )
(  559)       !!----------------------------------------------------------------------
(  560)       !!                  ***  routine mpparent  ***
(  561)       !!






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 112

(  562)       !! ** Purpose :   use an pvmfparent routine for T3E (key_mpp_shmem)
(  563)       !!              or  only return -1 (1)
(  564)       !!----------------------------------------------------------------------
(  565)       !! * Arguments
(  566)       INTEGER, INTENT(inout) ::   kparent_tid      ! ???
(  567)   
(  568) # 569
(  569)       ! MPI version : retour -1
(  570) 
(  571)       kparent_tid = -1
(  572) 
(  573) 
(  574) # 627
(  627)    END SUBROUTINE mpparent
(  628) 
(  629) 
(  630) # 654
(  654)    SUBROUTINE mpp_lnk_3d( ptab, cd_type, psgn, cd_mpp )
(  655)       !!----------------------------------------------------------------------
(  656)       !!                  ***  routine mpp_lnk_3d  ***
(  657)       !!
(  658)       !! ** Purpose :   Message passing manadgement
(  659)       !!
(  660)       !! ** Method  :   Use mppsend and mpprecv function for passing mask 
(  661)       !!      between processors following neighboring subdomains.
(  662)       !!            domain parameters
(  663)       !!                    nlci   : first dimension of the local subdomain
(  664)       !!                    nlcj   : second dimension of the local subdomain
(  665)       !!                    nbondi : mark for "east-west local boundary"
(  666)       !!                    nbondj : mark for "north-south local boundary"
(  667)       !!                    noea   : number for local neighboring processors 
(  668)       !!                    nowe   : number for local neighboring processors
(  669)       !!                    noso   : number for local neighboring processors
(  670)       !!                    nono   : number for local neighboring processors
(  671)       !!
(  672)       !! ** Action  :   ptab with update value at its periphery
(  673)       !!
(  674)       !!----------------------------------------------------------------------
(  675)       !! * Arguments
(  676)       CHARACTER(len=1) , INTENT( in ) ::   &
(  677)          cd_type       ! define the nature of ptab array grid-points
(  678)          !             ! = T , U , V , F , W points
(  679)          !             ! = S : T-point, north fold treatment ???
(  680)          !             ! = G : F-point, north fold treatment ???
(  681)       REAL(wp), INTENT( in ) ::   &
(  682)          psgn          ! control of the sign change
(  683)          !             !   = -1. , the sign is changed if north fold boundary
(  684)          !             !   =  1. , the sign is kept  if north fold boundary
(  685)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
(  686)          ptab          ! 3D array on which the boundary condition is applied
(  687)       CHARACTER(len=3), INTENT( in ), OPTIONAL ::    &
(  688)          cd_mpp        ! fill the overlap area only 
(  689) 
(  690)       !! * Local variables
(  691)       INTEGER ::   ji, jj, jk, jl                        ! dummy loop indices
(  692)       INTEGER ::   imigr, iihom, ijhom, iloc, ijt, iju   ! temporary integers
(  693)       INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
(  694)       INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 113

(  695)       !!----------------------------------------------------------------------
(  696) 
(  697)       ! 1. standard boundary treatment
(  698)       ! ------------------------------
(  699) 
(  700)       IF( PRESENT( cd_mpp ) ) THEN
(  701)          DO jj = nlcj+1, jpj   ! only fill extra allows last line
(  702)             ptab(1:nlci, jj, :) = ptab(1:nlci, nlej, :)
(  703)          END DO
(  704)          DO ji = nlci+1, jpi   ! only fill extra allows last column
(  705)             ptab(ji    , : , :) = ptab(nlei  , :   , :)
(  706)          END DO
(  707)       ELSE      
(  708) 
(  709)          !                                        ! East-West boundaries
(  710)          !                                        ! ====================
(  711)          IF( nbondi == 2 .AND.   &      ! Cyclic east-west
(  712)             &   (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
(  713)             ptab( 1 ,:,:) = ptab(jpim1,:,:)
(  714)             ptab(jpi,:,:) = ptab(  2  ,:,:)
(  715) 
(  716)          ELSE                           ! closed
(  717)             SELECT CASE ( cd_type )
(  718)             CASE ( 'T', 'U', 'V', 'W' )
(  719)                ptab(     1       :jpreci,:,:) = 0.e0
(  720)                ptab(nlci-jpreci+1:jpi   ,:,:) = 0.e0
(  721)             CASE ( 'F' )
(  722)                ptab(nlci-jpreci+1:jpi   ,:,:) = 0.e0
(  723)             END SELECT 
(  724)          ENDIF
(  725) 
(  726)          !                                        ! North-South boundaries
(  727)          !                                        ! ======================
(  728)          SELECT CASE ( cd_type )
(  729)          CASE ( 'T', 'U', 'V', 'W' )
(  730)             ptab(:,     1       :jprecj,:) = 0.e0
(  731)             ptab(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
(  732)          CASE ( 'F' )
(  733)             ptab(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
(  734)          END SELECT
(  735)      
(  736)       ENDIF
(  737) 
(  738)       ! 2. East and west directions exchange
(  739)       ! ------------------------------------
(  740) 
(  741)       ! 2.1 Read Dirichlet lateral conditions
(  742) 
(  743)       SELECT CASE ( nbondi )
(  744)       CASE ( -1, 0, 1 )    ! all exept 2 
(  745)          iihom = nlci-nreci
(  746)          DO jl = 1, jpreci
(  747)             t3ew(:,jl,:,1) = ptab(jpreci+jl,:,:)
(  748)             t3we(:,jl,:,1) = ptab(iihom +jl,:,:)
(  749)          END DO
(  750)       END SELECT
(  751) 
(  752)       ! 2.2 Migrations






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 114

(  753) 
(  754) # 773
(  773)       !! * Local variables   (MPI version)
(  774) 
(  775)       imigr = jpreci * jpj * jpk
(  776) 
(  777)       SELECT CASE ( nbondi ) 
(  778)       CASE ( -1 )
(  779)          CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req1 )
(  780)          CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
(  781)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  782)       CASE ( 0 )
(  783)          CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
(  784)          CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req2 )
(  785)          CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
(  786)          CALL mpprecv( 2, t3we(1,1,1,2), imigr )
(  787)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  788)          IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
(  789)       CASE ( 1 )
(  790)          CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
(  791)          CALL mpprecv( 2, t3we(1,1,1,2), imigr )
(  792)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  793)       END SELECT
(  794) 
(  795) # 796
(  796)       ! 2.3 Write Dirichlet lateral conditions
(  797) 
(  798)       iihom = nlci-jpreci
(  799) 
(  800)       SELECT CASE ( nbondi )
(  801)       CASE ( -1 )
(  802)          DO jl = 1, jpreci
(  803)             ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
(  804)          END DO
(  805)       CASE ( 0 ) 
(  806)          DO jl = 1, jpreci
(  807)             ptab(jl      ,:,:) = t3we(:,jl,:,2)
(  808)             ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
(  809)          END DO
(  810)       CASE ( 1 )
(  811)          DO jl = 1, jpreci
(  812)             ptab(jl      ,:,:) = t3we(:,jl,:,2)
(  813)          END DO
(  814)       END SELECT
(  815) 
(  816) 
(  817)       ! 3. North and south directions
(  818)       ! -----------------------------
(  819) 
(  820)       ! 3.1 Read Dirichlet lateral conditions
(  821) 
(  822)       IF( nbondj /= 2 ) THEN
(  823)          ijhom = nlcj-nrecj
(  824)          DO jl = 1, jprecj
(  825)             t3sn(:,jl,:,1) = ptab(:,ijhom +jl,:)
(  826)             t3ns(:,jl,:,1) = ptab(:,jprecj+jl,:)
(  827)          END DO
(  828)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 115

(  829) 
(  830)       ! 3.2 Migrations
(  831) 
(  832) # 851
(  851)       !! * Local variables   (MPI version)
(  852)   
(  853)       imigr=jprecj*jpi*jpk
(  854) 
(  855)       SELECT CASE ( nbondj )     
(  856)       CASE ( -1 )
(  857)          CALL mppsend( 4, t3sn(1,1,1,1), imigr, nono, ml_req1 )
(  858)          CALL mpprecv( 3, t3ns(1,1,1,2), imigr )
(  859)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  860)       CASE ( 0 )
(  861)          CALL mppsend( 3, t3ns(1,1,1,1), imigr, noso, ml_req1 )
(  862)          CALL mppsend( 4, t3sn(1,1,1,1), imigr, nono, ml_req2 )
(  863)          CALL mpprecv( 3, t3ns(1,1,1,2), imigr )
(  864)          CALL mpprecv( 4, t3sn(1,1,1,2), imigr )
(  865)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  866)          IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
(  867)       CASE ( 1 ) 
(  868)          CALL mppsend( 3, t3ns(1,1,1,1), imigr, noso, ml_req1 )
(  869)          CALL mpprecv( 4, t3sn(1,1,1,2), imigr )
(  870)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
(  871)       END SELECT
(  872) 
(  873) 
(  874) # 875
(  875)       ! 3.3 Write Dirichlet lateral conditions
(  876) 
(  877)       ijhom = nlcj-jprecj
(  878) 
(  879)       SELECT CASE ( nbondj )
(  880)       CASE ( -1 )
(  881)          DO jl = 1, jprecj
(  882)             ptab(:,ijhom+jl,:) = t3ns(:,jl,:,2)
(  883)          END DO
(  884)       CASE ( 0 ) 
(  885)          DO jl = 1, jprecj
(  886)             ptab(:,jl      ,:) = t3sn(:,jl,:,2)
(  887)             ptab(:,ijhom+jl,:) = t3ns(:,jl,:,2)
(  888)          END DO
(  889)       CASE ( 1 )
(  890)          DO jl = 1, jprecj
(  891)             ptab(:,jl,:) = t3sn(:,jl,:,2)
(  892)          END DO
(  893)       END SELECT
(  894) 
(  895) 
(  896)       ! 4. north fold treatment
(  897)       ! -----------------------
(  898) 
(  899)       IF (PRESENT(cd_mpp)) THEN
(  900)          ! No north fold treatment (it is assumed to be already OK)
(  901)      
(  902)       ELSE      
(  903) 
(  904)       ! 4.1 treatment without exchange (jpni odd)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 116

(  905)       !     T-point pivot  
(  906) 
(  907)       SELECT CASE ( jpni )
(  908) 
(  909)       CASE ( 1 )  ! only one proc along I, no mpp exchange
(  910) 
(  911)          SELECT CASE ( npolj )
(  912)   
(  913)          CASE ( 3 , 4 )    ! T pivot
(  914)             iloc = jpiglo - 2 * ( nimpp - 1 )
(  915) 
(  916)             SELECT CASE ( cd_type )
(  917) 
(  918)             CASE ( 'T' , 'S', 'W' )
(  919)                DO jk = 1, jpk
(  920)                   DO ji = 2, nlci
(  921)                      ijt=iloc-ji+2
(  922)                      ptab(ji,nlcj,jk) = psgn * ptab(ijt,nlcj-2,jk)
(  923)                   END DO
(  924)                   DO ji = nlci/2+1, nlci
(  925)                      ijt=iloc-ji+2
(  926)                      ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-1,jk)
(  927)                   END DO
(  928)                END DO
(  929)           
(  930)             CASE ( 'U' )
(  931)                DO jk = 1, jpk
(  932)                   DO ji = 1, nlci-1
(  933)                      iju=iloc-ji+1
(  934)                      ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-2,jk)
(  935)                   END DO
(  936)                   DO ji = nlci/2, nlci-1
(  937)                      iju=iloc-ji+1
(  938)                      ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-1,jk)
(  939)                   END DO
(  940)                END DO
(  941) 
(  942)             CASE ( 'V' )
(  943)                DO jk = 1, jpk
(  944)                   DO ji = 2, nlci
(  945)                      ijt=iloc-ji+2
(  946)                      ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-2,jk)
(  947)                      ptab(ji,nlcj  ,jk) = psgn * ptab(ijt,nlcj-3,jk)
(  948)                   END DO
(  949)                END DO
(  950) 
(  951)             CASE ( 'F', 'G' )
(  952)                DO jk = 1, jpk
(  953)                   DO ji = 1, nlci-1
(  954)                      iju=iloc-ji+1
(  955)                      ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-2,jk)
(  956)                      ptab(ji,nlcj  ,jk) = psgn * ptab(iju,nlcj-3,jk)
(  957)                   END DO
(  958)                END DO
(  959)   
(  960)           END SELECT
(  961)        
(  962)          CASE ( 5 , 6 ) ! F pivot






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 117

(  963)             iloc=jpiglo-2*(nimpp-1)
(  964)   
(  965)             SELECT CASE ( cd_type )
(  966) 
(  967)             CASE ( 'T' , 'S', 'W' )
(  968)                DO jk = 1, jpk
(  969)                   DO ji = 1, nlci
(  970)                      ijt=iloc-ji+1
(  971)                      ptab(ji,nlcj,jk) = psgn * ptab(ijt,nlcj-1,jk)
(  972)                   END DO
(  973)                END DO
(  974) 
(  975)             CASE ( 'U' )
(  976)                DO jk = 1, jpk
(  977)                   DO ji = 1, nlci-1
(  978)                      iju=iloc-ji
(  979)                      ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-1,jk)
(  980)                   END DO
(  981)                END DO
(  982) 
(  983)             CASE ( 'V' )
(  984)                DO jk = 1, jpk
(  985)                   DO ji = 1, nlci
(  986)                      ijt=iloc-ji+1
(  987)                      ptab(ji,nlcj  ,jk) = psgn * ptab(ijt,nlcj-2,jk)
(  988)                   END DO
(  989)                   DO ji = nlci/2+1, nlci
(  990)                      ijt=iloc-ji+1
(  991)                      ptab(ji,nlcj-1,jk) = psgn * ptab(ijt,nlcj-1,jk)
(  992)                   END DO
(  993)                END DO
(  994) 
(  995)             CASE ( 'F', 'G' )
(  996)                DO jk = 1, jpk
(  997)                   DO ji = 1, nlci-1
(  998)                      iju=iloc-ji
(  999)                      ptab(ji,nlcj,jk) = psgn * ptab(iju,nlcj-2,jk)
( 1000)                   END DO
( 1001)                   DO ji = nlci/2+1, nlci-1
( 1002)                      iju=iloc-ji
( 1003)                      ptab(ji,nlcj-1,jk) = psgn * ptab(iju,nlcj-1,jk)
( 1004)                   END DO
( 1005)                END DO
( 1006)             END SELECT  ! cd_type
( 1007) 
( 1008)          END SELECT     !  npolj
( 1009)   
( 1010)       CASE DEFAULT ! more than 1 proc along I
( 1011)          IF ( npolj /= 0 ) CALL mpp_lbc_north (ptab, cd_type, psgn)  ! only for northern procs.
( 1012) 
( 1013)       END SELECT ! jpni 
( 1014) 
( 1015)       ENDIF
( 1016)       
( 1017) 
( 1018)       ! 5. East and west directions exchange
( 1019)       ! ------------------------------------
( 1020) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 118

( 1021)       SELECT CASE ( npolj )
( 1022) 
( 1023)       CASE ( 3, 4, 5, 6 )
( 1024) 
( 1025)          ! 5.1 Read Dirichlet lateral conditions
( 1026) 
( 1027)          SELECT CASE ( nbondi )
( 1028) 
( 1029)          CASE ( -1, 0, 1 )
( 1030)             iihom = nlci-nreci
( 1031)             DO jl = 1, jpreci
( 1032)                t3ew(:,jl,:,1) = ptab(jpreci+jl,:,:)
( 1033)                t3we(:,jl,:,1) = ptab(iihom +jl,:,:)
( 1034)             END DO
( 1035) 
( 1036)          END SELECT
( 1037) 
( 1038)          ! 5.2 Migrations
( 1039) 
( 1040) # 1059
( 1059)          !! MPI version
( 1060) 
( 1061)          imigr=jpreci*jpj*jpk
( 1062)   
( 1063)          SELECT CASE ( nbondi )
( 1064)          CASE ( -1 )
( 1065)             CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req1 )
( 1066)             CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
( 1067)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1068)          CASE ( 0 )
( 1069)             CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
( 1070)             CALL mppsend( 2, t3we(1,1,1,1), imigr, noea, ml_req2 )
( 1071)             CALL mpprecv( 1, t3ew(1,1,1,2), imigr )
( 1072)             CALL mpprecv( 2, t3we(1,1,1,2), imigr )
( 1073)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1074)             IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 1075)          CASE ( 1 )
( 1076)             CALL mppsend( 1, t3ew(1,1,1,1), imigr, nowe, ml_req1 )
( 1077)             CALL mpprecv( 2, t3we(1,1,1,2), imigr )
( 1078)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1079)          END SELECT
( 1080) 
( 1081) # 1082
( 1082)          ! 5.3 Write Dirichlet lateral conditions
( 1083) 
( 1084)          iihom = nlci-jpreci
( 1085) 
( 1086)          SELECT CASE ( nbondi)
( 1087)          CASE ( -1 )
( 1088)             DO jl = 1, jpreci
( 1089)                ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
( 1090)             END DO
( 1091)          CASE ( 0 ) 
( 1092)             DO jl = 1, jpreci
( 1093)                ptab(jl      ,:,:) = t3we(:,jl,:,2)
( 1094)                ptab(iihom+jl,:,:) = t3ew(:,jl,:,2)
( 1095)             END DO
( 1096)          CASE ( 1 )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 119

( 1097)             DO jl = 1, jpreci
( 1098)                ptab(jl      ,:,:) = t3we(:,jl,:,2)
( 1099)             END DO
( 1100)          END SELECT
( 1101) 
( 1102)       END SELECT    ! npolj 
( 1103) 
( 1104)    END SUBROUTINE mpp_lnk_3d
( 1105) 
( 1106) 
( 1107)    SUBROUTINE mpp_lnk_2d( pt2d, cd_type, psgn, cd_mpp )
( 1108)       !!----------------------------------------------------------------------
( 1109)       !!                  ***  routine mpp_lnk_2d  ***
( 1110)       !!                  
( 1111)       !! ** Purpose :   Message passing manadgement for 2d array
( 1112)       !!
( 1113)       !! ** Method  :   Use mppsend and mpprecv function for passing mask 
( 1114)       !!      between processors following neighboring subdomains.
( 1115)       !!            domain parameters
( 1116)       !!                    nlci   : first dimension of the local subdomain
( 1117)       !!                    nlcj   : second dimension of the local subdomain
( 1118)       !!                    nbondi : mark for "east-west local boundary"
( 1119)       !!                    nbondj : mark for "north-south local boundary"
( 1120)       !!                    noea   : number for local neighboring processors 
( 1121)       !!                    nowe   : number for local neighboring processors
( 1122)       !!                    noso   : number for local neighboring processors
( 1123)       !!                    nono   : number for local neighboring processors
( 1124)       !!
( 1125)       !!----------------------------------------------------------------------
( 1126)       !! * Arguments
( 1127)       CHARACTER(len=1) , INTENT( in ) ::   &
( 1128)          cd_type       ! define the nature of pt2d array grid-points
( 1129)          !             !  = T , U , V , F , W 
( 1130)          !             !  = S : T-point, north fold treatment
( 1131)          !             !  = G : F-point, north fold treatment
( 1132)          !             !  = I : sea-ice velocity at F-point with index shift
( 1133)       REAL(wp), INTENT( in ) ::   &
( 1134)          psgn          ! control of the sign change
( 1135)          !             !   = -1. , the sign is changed if north fold boundary
( 1136)          !             !   =  1. , the sign is kept  if north fold boundary
( 1137)       REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
( 1138)          pt2d          ! 2D array on which the boundary condition is applied
( 1139)       CHARACTER(len=3), INTENT( in ), OPTIONAL ::    &
( 1140)          cd_mpp        ! fill the overlap area only 
( 1141) 
( 1142)       !! * Local variables
( 1143)       INTEGER  ::   ji, jj, jl      ! dummy loop indices
( 1144)       INTEGER  ::   &
( 1145)          imigr, iihom, ijhom,    &  ! temporary integers
( 1146)          iloc, ijt, iju             !    "          "
( 1147)       INTEGER  ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
( 1148)       INTEGER  ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
( 1149)       !!----------------------------------------------------------------------
( 1150) 
( 1151)       ! 1. standard boundary treatment
( 1152)       ! ------------------------------
( 1153)       IF (PRESENT(cd_mpp)) THEN
( 1154)          DO jj = nlcj+1, jpj   ! only fill extra allows last line






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 120

( 1155)             pt2d(1:nlci, jj) = pt2d(1:nlci, nlej)
( 1156)          END DO
( 1157)          DO ji = nlci+1, jpi   ! only fill extra allows last column
( 1158)             pt2d(ji    , : ) = pt2d(nlei  , :   )
( 1159)          END DO     
( 1160)       ELSE      
( 1161) 
( 1162)          !                                        ! East-West boundaries
( 1163)          !                                        ! ====================
( 1164)          IF( nbondi == 2 .AND.   &      ! Cyclic east-west
( 1165)             &    (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
( 1166)             pt2d( 1 ,:) = pt2d(jpim1,:)
( 1167)             pt2d(jpi,:) = pt2d(  2  ,:)
( 1168) 
( 1169)          ELSE                           ! ... closed
( 1170)             SELECT CASE ( cd_type )
( 1171)             CASE ( 'T', 'U', 'V', 'W' , 'I' )
( 1172)                pt2d(     1       :jpreci,:) = 0.e0
( 1173)                pt2d(nlci-jpreci+1:jpi   ,:) = 0.e0
( 1174)             CASE ( 'F' )
( 1175)                pt2d(nlci-jpreci+1:jpi   ,:) = 0.e0
( 1176)             END SELECT
( 1177)          ENDIF
( 1178) 
( 1179)          !                                        ! North-South boundaries
( 1180)          !                                        ! ======================
( 1181)          SELECT CASE ( cd_type )
( 1182)          CASE ( 'T', 'U', 'V', 'W' , 'I' )
( 1183)             pt2d(:,     1       :jprecj) = 0.e0
( 1184)             pt2d(:,nlcj-jprecj+1:jpj   ) = 0.e0
( 1185)          CASE ( 'F' )
( 1186)             pt2d(:,nlcj-jprecj+1:jpj   ) = 0.e0
( 1187)          END SELECT
( 1188) 
( 1189)       ENDIF
( 1190) 
( 1191) 
( 1192)       ! 2. East and west directions
( 1193)       ! ---------------------------
( 1194) 
( 1195)       ! 2.1 Read Dirichlet lateral conditions
( 1196) 
( 1197)       SELECT CASE ( nbondi )
( 1198)       CASE ( -1, 0, 1 )    ! all except 2
( 1199)          iihom = nlci-nreci
( 1200)          DO jl = 1, jpreci
( 1201)             t2ew(:,jl,1) = pt2d(jpreci+jl,:)
( 1202)             t2we(:,jl,1) = pt2d(iihom +jl,:)
( 1203)          END DO
( 1204)       END SELECT
( 1205) 
( 1206)       ! 2.2 Migrations
( 1207) 
( 1208) # 1227
( 1227)       !! * MPI version
( 1228) 
( 1229)       imigr = jpreci * jpj
( 1230) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 121

( 1231)       SELECT CASE ( nbondi )
( 1232)       CASE ( -1 )
( 1233)          CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req1 )
( 1234)          CALL mpprecv( 1, t2ew(1,1,2), imigr )
( 1235)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1236)       CASE ( 0 )
( 1237)          CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
( 1238)          CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req2 )
( 1239)          CALL mpprecv( 1, t2ew(1,1,2), imigr )
( 1240)          CALL mpprecv( 2, t2we(1,1,2), imigr )
( 1241)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1242)          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 1243)       CASE ( 1 )
( 1244)          CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
( 1245)          CALL mpprecv( 2, t2we(1,1,2), imigr )
( 1246)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1247)       END SELECT
( 1248) 
( 1249) 
( 1250) # 1251
( 1251)       ! 2.3 Write Dirichlet lateral conditions
( 1252) 
( 1253)       iihom = nlci - jpreci
( 1254)       SELECT CASE ( nbondi )
( 1255)       CASE ( -1 )
( 1256)          DO jl = 1, jpreci
( 1257)             pt2d(iihom+jl,:) = t2ew(:,jl,2)
( 1258)          END DO
( 1259)       CASE ( 0 )
( 1260)          DO jl = 1, jpreci
( 1261)             pt2d(jl      ,:) = t2we(:,jl,2)
( 1262)             pt2d(iihom+jl,:) = t2ew(:,jl,2)
( 1263)          END DO
( 1264)       CASE ( 1 )
( 1265)          DO jl = 1, jpreci
( 1266)             pt2d(jl      ,:) = t2we(:,jl,2)
( 1267)          END DO
( 1268)       END SELECT
( 1269) 
( 1270) 
( 1271)       ! 3. North and south directions
( 1272)       ! -----------------------------
( 1273) 
( 1274)       ! 3.1 Read Dirichlet lateral conditions
( 1275) 
( 1276)       IF( nbondj /= 2 ) THEN
( 1277)          ijhom = nlcj-nrecj
( 1278)          DO jl = 1, jprecj
( 1279)             t2sn(:,jl,1) = pt2d(:,ijhom +jl)
( 1280)             t2ns(:,jl,1) = pt2d(:,jprecj+jl)
( 1281)          END DO
( 1282)       ENDIF
( 1283) 
( 1284)       ! 3.2 Migrations
( 1285) 
( 1286) # 1304
( 1304)       !! * MPI version
( 1305) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 122

( 1306)       imigr = jprecj * jpi
( 1307) 
( 1308)       SELECT CASE ( nbondj )
( 1309)       CASE ( -1 )
( 1310)          CALL mppsend( 4, t2sn(1,1,1), imigr, nono, ml_req1 )
( 1311)          CALL mpprecv( 3, t2ns(1,1,2), imigr )
( 1312)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1313)       CASE ( 0 )
( 1314)          CALL mppsend( 3, t2ns(1,1,1), imigr, noso, ml_req1 )
( 1315)          CALL mppsend( 4, t2sn(1,1,1), imigr, nono, ml_req2 )
( 1316)          CALL mpprecv( 3, t2ns(1,1,2), imigr )
( 1317)          CALL mpprecv( 4, t2sn(1,1,2), imigr )
( 1318)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1319)          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 1320)       CASE ( 1 )
( 1321)          CALL mppsend( 3, t2ns(1,1,1), imigr, noso, ml_req1 )
( 1322)          CALL mpprecv( 4, t2sn(1,1,2), imigr )
( 1323)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1324)       END SELECT
( 1325)   
( 1326) 
( 1327) # 1328
( 1328)       ! 3.3 Write Dirichlet lateral conditions
( 1329) 
( 1330)       ijhom = nlcj - jprecj
( 1331) 
( 1332)       SELECT CASE ( nbondj )
( 1333)       CASE ( -1 )
( 1334)          DO jl = 1, jprecj
( 1335)             pt2d(:,ijhom+jl) = t2ns(:,jl,2)
( 1336)          END DO
( 1337)       CASE ( 0 )
( 1338)          DO jl = 1, jprecj
( 1339)             pt2d(:,jl      ) = t2sn(:,jl,2)
( 1340)             pt2d(:,ijhom+jl) = t2ns(:,jl,2)
( 1341)          END DO
( 1342)       CASE ( 1 ) 
( 1343)          DO jl = 1, jprecj
( 1344)             pt2d(:,jl      ) = t2sn(:,jl,2)
( 1345)          END DO
( 1346)       END SELECT 
( 1347)   
( 1348) 
( 1349)       ! 4. north fold treatment
( 1350)       ! -----------------------
( 1351)   
( 1352)       IF (PRESENT(cd_mpp)) THEN
( 1353)          ! No north fold treatment (it is assumed to be already OK)
( 1354)      
( 1355)       ELSE      
( 1356) 
( 1357)       ! 4.1 treatment without exchange (jpni odd)
( 1358)       
( 1359)       SELECT CASE ( jpni )
( 1360)   
( 1361)       CASE ( 1 ) ! only one proc along I, no mpp exchange
( 1362)   
( 1363)          SELECT CASE ( npolj )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 123

( 1364)   
( 1365)          CASE ( 3 , 4 )   !  T pivot
( 1366)             iloc = jpiglo - 2 * ( nimpp - 1 )
( 1367)   
( 1368)             SELECT CASE ( cd_type )
( 1369)   
( 1370)             CASE ( 'T' , 'S', 'W' )
( 1371)                DO ji = 2, nlci
( 1372)                   ijt=iloc-ji+2
( 1373)                   pt2d(ji,nlcj) = psgn * pt2d(ijt,nlcj-2)
( 1374)                END DO
( 1375)                DO ji = nlci/2+1, nlci
( 1376)                   ijt=iloc-ji+2
( 1377)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
( 1378)                END DO
( 1379)   
( 1380)             CASE ( 'U' )
( 1381)                DO ji = 1, nlci-1
( 1382)                   iju=iloc-ji+1
( 1383)                   pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-2)
( 1384)                END DO
( 1385)                DO ji = nlci/2, nlci-1
( 1386)                   iju=iloc-ji+1
( 1387)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
( 1388)                END DO
( 1389)   
( 1390)             CASE ( 'V' )
( 1391)                DO ji = 2, nlci
( 1392)                   ijt=iloc-ji+2
( 1393)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-2)
( 1394)                   pt2d(ji,nlcj  ) = psgn * pt2d(ijt,nlcj-3)
( 1395)                END DO
( 1396)   
( 1397)             CASE ( 'F', 'G' )
( 1398)                DO ji = 1, nlci-1
( 1399)                   iju=iloc-ji+1
( 1400)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-2)
( 1401)                   pt2d(ji,nlcj  ) = psgn * pt2d(iju,nlcj-3)
( 1402)                END DO
( 1403)   
( 1404)             CASE ( 'I' )                                  ! ice U-V point
( 1405)                pt2d(2,nlcj) = psgn * pt2d(3,nlcj-1)
( 1406)                DO ji = 3, nlci
( 1407)                   iju = iloc - ji + 3
( 1408)                   pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-1)
( 1409)                END DO
( 1410)   
( 1411)             END SELECT
( 1412)   
( 1413)          CASE ( 5 , 6 )                 ! F pivot
( 1414)             iloc=jpiglo-2*(nimpp-1)
( 1415)   
( 1416)             SELECT CASE (cd_type )
( 1417)   
( 1418)             CASE ( 'T', 'S', 'W' )
( 1419)                DO ji = 1, nlci
( 1420)                   ijt=iloc-ji+1
( 1421)                   pt2d(ji,nlcj) = psgn * pt2d(ijt,nlcj-1)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 124

( 1422)                END DO
( 1423)   
( 1424)             CASE ( 'U' )
( 1425)                DO ji = 1, nlci-1
( 1426)                   iju=iloc-ji
( 1427)                   pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-1)
( 1428)                END DO
( 1429) 
( 1430)             CASE ( 'V' )
( 1431)                DO ji = 1, nlci
( 1432)                   ijt=iloc-ji+1
( 1433)                   pt2d(ji,nlcj  ) = psgn * pt2d(ijt,nlcj-2)
( 1434)                END DO
( 1435)                DO ji = nlci/2+1, nlci
( 1436)                   ijt=iloc-ji+1
( 1437)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
( 1438)                END DO
( 1439)   
( 1440)             CASE ( 'F', 'G' )
( 1441)                DO ji = 1, nlci-1
( 1442)                   iju=iloc-ji
( 1443)                   pt2d(ji,nlcj) = psgn * pt2d(iju,nlcj-2)
( 1444)                END DO
( 1445)                DO ji = nlci/2+1, nlci-1
( 1446)                   iju=iloc-ji
( 1447)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
( 1448)                END DO
( 1449)   
( 1450)             CASE ( 'I' )                                  ! ice U-V point
( 1451)                pt2d( 2 ,nlcj) = 0.e0
( 1452)                DO ji = 2 , nlci-1
( 1453)                   ijt = iloc - ji + 2
( 1454)                   pt2d(ji,nlcj)= 0.5 * ( pt2d(ji,nlcj-1) + psgn * pt2d(ijt,nlcj-1) )
( 1455)                END DO
( 1456)   
( 1457)             END SELECT   ! cd_type
( 1458)   
( 1459)          END SELECT   ! npolj
( 1460) 
( 1461)       CASE DEFAULT   ! more than 1 proc along I
( 1462)          IF( npolj /= 0 )   CALL mpp_lbc_north( pt2d, cd_type, psgn )   ! only for northern procs.
( 1463) 
( 1464)       END SELECT   ! jpni
( 1465) 
( 1466)       ENDIF
( 1467) 
( 1468)       ! 5. East and west directions
( 1469)       ! ---------------------------
( 1470) 
( 1471)       SELECT CASE ( npolj )
( 1472) 
( 1473)       CASE ( 3, 4, 5, 6 )
( 1474) 
( 1475)          ! 5.1 Read Dirichlet lateral conditions
( 1476) 
( 1477)          SELECT CASE ( nbondi )
( 1478)          CASE ( -1, 0, 1 )
( 1479)             iihom = nlci-nreci






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 125

( 1480)             DO jl = 1, jpreci
( 1481)                DO jj = 1, jpj
( 1482)                   t2ew(jj,jl,1) = pt2d(jpreci+jl,jj)
( 1483)                   t2we(jj,jl,1) = pt2d(iihom +jl,jj)
( 1484)                END DO
( 1485)             END DO
( 1486)          END SELECT
( 1487) 
( 1488)          ! 5.2 Migrations
( 1489) 
( 1490) # 1509
( 1509)          !! * MPI version
( 1510)   
( 1511)          imigr=jpreci*jpj
( 1512)   
( 1513)          SELECT CASE ( nbondi )
( 1514)          CASE ( -1 )
( 1515)             CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req1 )
( 1516)             CALL mpprecv( 1, t2ew(1,1,2), imigr )
( 1517)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1518)          CASE ( 0 )
( 1519)             CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
( 1520)             CALL mppsend( 2, t2we(1,1,1), imigr, noea, ml_req2 )
( 1521)             CALL mpprecv( 1, t2ew(1,1,2), imigr )
( 1522)             CALL mpprecv( 2, t2we(1,1,2), imigr )
( 1523)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1524)             IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 1525)          CASE ( 1 )
( 1526)             CALL mppsend( 1, t2ew(1,1,1), imigr, nowe, ml_req1 )
( 1527)             CALL mpprecv( 2, t2we(1,1,2), imigr )
( 1528)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 1529)          END SELECT 
( 1530) 
( 1531) # 1532
( 1532)          ! 5.3 Write Dirichlet lateral conditions
( 1533)   
( 1534)          iihom = nlci - jpreci
( 1535)   
( 1536)          SELECT CASE ( nbondi )
( 1537)          CASE ( -1 )
( 1538)             DO jl = 1, jpreci
( 1539)                pt2d(iihom+jl,:) = t2ew(:,jl,2)
( 1540)             END DO
( 1541)          CASE ( 0 )
( 1542)             DO jl = 1, jpreci
( 1543)                pt2d(jl      ,:) = t2we(:,jl,2)
( 1544)                pt2d(iihom+jl,:) = t2ew(:,jl,2)
( 1545)             END DO
( 1546)          CASE ( 1 )
( 1547)             DO jl = 1, jpreci
( 1548)                pt2d(jl,:) = t2we(:,jl,2)
( 1549)             END DO
( 1550)          END SELECT 
( 1551)   
( 1552)       END SELECT   ! npolj
( 1553)   
( 1554)    END SUBROUTINE mpp_lnk_2d
( 1555) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 126

( 1556) 
( 1557)    SUBROUTINE mpp_lnk_3d_gather( ptab1, cd_type1, ptab2, cd_type2, psgn )
( 1558)       !!----------------------------------------------------------------------
( 1559)       !!                  ***  routine mpp_lnk_3d_gather  ***
( 1560)       !!
( 1561)       !! ** Purpose :   Message passing manadgement for two 3D arrays
( 1562)       !!
( 1563)       !! ** Method  :   Use mppsend and mpprecv function for passing mask 
( 1564)       !!      between processors following neighboring subdomains.
( 1565)       !!            domain parameters
( 1566)       !!                    nlci   : first dimension of the local subdomain
( 1567)       !!                    nlcj   : second dimension of the local subdomain
( 1568)       !!                    nbondi : mark for "east-west local boundary"
( 1569)       !!                    nbondj : mark for "north-south local boundary"
( 1570)       !!                    noea   : number for local neighboring processors 
( 1571)       !!                    nowe   : number for local neighboring processors
( 1572)       !!                    noso   : number for local neighboring processors
( 1573)       !!                    nono   : number for local neighboring processors
( 1574)       !!
( 1575)       !! ** Action  :   ptab1 and ptab2  with update value at its periphery
( 1576)       !!
( 1577)       !!----------------------------------------------------------------------
( 1578)       !! * Arguments
( 1579)       CHARACTER(len=1) , INTENT( in ) ::   &
( 1580)          cd_type1, cd_type2       ! define the nature of ptab array grid-points
( 1581)          !                        ! = T , U , V , F , W points
( 1582)          !                        ! = S : T-point, north fold treatment ???
( 1583)          !                        ! = G : F-point, north fold treatment ???
( 1584)       REAL(wp), INTENT( in ) ::   &
( 1585)          psgn          ! control of the sign change
( 1586)          !             !   = -1. , the sign is changed if north fold boundary
( 1587)          !             !   =  1. , the sign is kept  if north fold boundary
( 1588)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
( 1589)          ptab1, ptab2             ! 3D array on which the boundary condition is applied
( 1590) 
( 1591)       !! * Local variables
( 1592)       INTEGER ::   ji, jk, jl   ! dummy loop indices
( 1593)       INTEGER ::   imigr, iihom, ijhom, iloc, ijt, iju   ! temporary integers
( 1594)       INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
( 1595)       INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
( 1596)       !!----------------------------------------------------------------------
( 1597) 
( 1598)       ! 1. standard boundary treatment
( 1599)       ! ------------------------------
( 1600)       !                                        ! East-West boundaries
( 1601)       !                                        ! ====================
( 1602)       IF( nbondi == 2 .AND.   &      ! Cyclic east-west
( 1603)          &   (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
( 1604)          ptab1( 1 ,:,:) = ptab1(jpim1,:,:)
( 1605)          ptab1(jpi,:,:) = ptab1(  2  ,:,:)
( 1606)          ptab2( 1 ,:,:) = ptab2(jpim1,:,:)
( 1607)          ptab2(jpi,:,:) = ptab2(  2  ,:,:)
( 1608) 
( 1609)       ELSE                           ! closed
( 1610)          SELECT CASE ( cd_type1 )
( 1611)          CASE ( 'T', 'U', 'V', 'W' )
( 1612)             ptab1(     1       :jpreci,:,:) = 0.e0
( 1613)             ptab1(nlci-jpreci+1:jpi   ,:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 127

( 1614)          CASE ( 'F' )
( 1615)             ptab1(nlci-jpreci+1:jpi   ,:,:) = 0.e0
( 1616)          END SELECT 
( 1617)          SELECT CASE ( cd_type2 )
( 1618)          CASE ( 'T', 'U', 'V', 'W' )
( 1619)             ptab2(     1       :jpreci,:,:) = 0.e0
( 1620)             ptab2(nlci-jpreci+1:jpi   ,:,:) = 0.e0
( 1621)          CASE ( 'F' )
( 1622)             ptab2(nlci-jpreci+1:jpi   ,:,:) = 0.e0
( 1623)          END SELECT 
( 1624)       ENDIF
( 1625) 
( 1626)       !                                        ! North-South boundaries
( 1627)       !                                        ! ======================
( 1628)       SELECT CASE ( cd_type1 )
( 1629)       CASE ( 'T', 'U', 'V', 'W' )
( 1630)          ptab1(:,     1       :jprecj,:) = 0.e0
( 1631)          ptab1(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
( 1632)       CASE ( 'F' )
( 1633)          ptab1(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
( 1634)       END SELECT
( 1635) 
( 1636)       SELECT CASE ( cd_type2 )
( 1637)       CASE ( 'T', 'U', 'V', 'W' )
( 1638)          ptab2(:,     1       :jprecj,:) = 0.e0
( 1639)          ptab2(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
( 1640)       CASE ( 'F' )
( 1641)          ptab2(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
( 1642)       END SELECT
( 1643) 
( 1644) 
( 1645)       ! 2. East and west directions exchange
( 1646)       ! ------------------------------------
( 1647) 
( 1648)       ! 2.1 Read Dirichlet lateral conditions
( 1649) 
( 1650)       SELECT CASE ( nbondi )
( 1651)       CASE ( -1, 0, 1 )    ! all exept 2 
( 1652)          iihom = nlci-nreci
( 1653)          DO jl = 1, jpreci
( 1654)             t4ew(:,jl,:,1,1) = ptab1(jpreci+jl,:,:)
( 1655)             t4we(:,jl,:,1,1) = ptab1(iihom +jl,:,:)
( 1656)             t4ew(:,jl,:,2,1) = ptab2(jpreci+jl,:,:)
( 1657)             t4we(:,jl,:,2,1) = ptab2(iihom +jl,:,:)
( 1658)          END DO
( 1659)       END SELECT
( 1660) 
( 1661)       ! 2.2 Migrations
( 1662) 
( 1663) # 1682
( 1682)       !! * Local variables   (MPI version)
( 1683) 
( 1684)       imigr = jpreci * jpj * jpk *2
( 1685) 
( 1686)       SELECT CASE ( nbondi ) 
( 1687)       CASE ( -1 )
( 1688)          CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req1 )
( 1689)          CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 128

( 1690)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1691)       CASE ( 0 )
( 1692)          CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
( 1693)          CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req2 )
( 1694)          CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
( 1695)          CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
( 1696)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1697)          IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
( 1698)       CASE ( 1 )
( 1699)          CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
( 1700)          CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
( 1701)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1702)       END SELECT
( 1703) 
( 1704) # 1705
( 1705)       ! 2.3 Write Dirichlet lateral conditions
( 1706) 
( 1707)       iihom = nlci-jpreci
( 1708) 
( 1709)       SELECT CASE ( nbondi )
( 1710)       CASE ( -1 )
( 1711)          DO jl = 1, jpreci
( 1712)             ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
( 1713)             ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
( 1714)          END DO
( 1715)       CASE ( 0 ) 
( 1716)          DO jl = 1, jpreci
( 1717)             ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
( 1718)             ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
( 1719)             ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
( 1720)             ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
( 1721)          END DO
( 1722)       CASE ( 1 )
( 1723)          DO jl = 1, jpreci
( 1724)             ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
( 1725)             ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
( 1726)          END DO
( 1727)       END SELECT
( 1728) 
( 1729) 
( 1730)       ! 3. North and south directions
( 1731)       ! -----------------------------
( 1732) 
( 1733)       ! 3.1 Read Dirichlet lateral conditions
( 1734) 
( 1735)       IF( nbondj /= 2 ) THEN
( 1736)          ijhom = nlcj-nrecj
( 1737)          DO jl = 1, jprecj
( 1738)             t4sn(:,jl,:,1,1) = ptab1(:,ijhom +jl,:)
( 1739)             t4ns(:,jl,:,1,1) = ptab1(:,jprecj+jl,:)
( 1740)             t4sn(:,jl,:,2,1) = ptab2(:,ijhom +jl,:)
( 1741)             t4ns(:,jl,:,2,1) = ptab2(:,jprecj+jl,:)
( 1742)          END DO
( 1743)       ENDIF
( 1744) 
( 1745)       ! 3.2 Migrations
( 1746) 
( 1747) # 1766






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 129

( 1766)       !! * Local variables   (MPI version)
( 1767)   
( 1768)       imigr=jprecj * jpi * jpk * 2
( 1769) 
( 1770)       SELECT CASE ( nbondj )     
( 1771)       CASE ( -1 )
( 1772)          CALL mppsend( 4, t4sn(1,1,1,1,1), imigr, nono, ml_req1 )
( 1773)          CALL mpprecv( 3, t4ns(1,1,1,1,2), imigr )
( 1774)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1775)       CASE ( 0 )
( 1776)          CALL mppsend( 3, t4ns(1,1,1,1,1), imigr, noso, ml_req1 )
( 1777)          CALL mppsend( 4, t4sn(1,1,1,1,1), imigr, nono, ml_req2 )
( 1778)          CALL mpprecv( 3, t4ns(1,1,1,1,2), imigr )
( 1779)          CALL mpprecv( 4, t4sn(1,1,1,1,2), imigr )
( 1780)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1781)          IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
( 1782)       CASE ( 1 ) 
( 1783)          CALL mppsend( 3, t4ns(1,1,1,1,1), imigr, noso, ml_req1 )
( 1784)          CALL mpprecv( 4, t4sn(1,1,1,1,2), imigr )
( 1785)          IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
( 1786)       END SELECT
( 1787) 
( 1788) 
( 1789) # 1790
( 1790)       ! 3.3 Write Dirichlet lateral conditions
( 1791) 
( 1792)       ijhom = nlcj-jprecj
( 1793) 
( 1794)       SELECT CASE ( nbondj )
( 1795)       CASE ( -1 )
( 1796)          DO jl = 1, jprecj
( 1797)             ptab1(:,ijhom+jl,:) = t4ns(:,jl,:,1,2)
( 1798)             ptab2(:,ijhom+jl,:) = t4ns(:,jl,:,2,2)
( 1799)          END DO
( 1800)       CASE ( 0 ) 
( 1801)          DO jl = 1, jprecj
( 1802)             ptab1(:,jl      ,:) = t4sn(:,jl,:,1,2)
( 1803)             ptab1(:,ijhom+jl,:) = t4ns(:,jl,:,1,2)
( 1804)             ptab2(:,jl      ,:) = t4sn(:,jl,:,2,2)
( 1805)             ptab2(:,ijhom+jl,:) = t4ns(:,jl,:,2,2)
( 1806)          END DO
( 1807)       CASE ( 1 )
( 1808)          DO jl = 1, jprecj
( 1809)             ptab1(:,jl,:) = t4sn(:,jl,:,1,2)
( 1810)             ptab2(:,jl,:) = t4sn(:,jl,:,2,2)
( 1811)          END DO
( 1812)       END SELECT
( 1813) 
( 1814) 
( 1815)       ! 4. north fold treatment
( 1816)       ! -----------------------
( 1817) 
( 1818)       ! 4.1 treatment without exchange (jpni odd)
( 1819)       !     T-point pivot  
( 1820) 
( 1821)       SELECT CASE ( jpni )
( 1822) 
( 1823)       CASE ( 1 )  ! only one proc along I, no mpp exchange






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 130

( 1824) 
( 1825)       SELECT CASE ( npolj )
( 1826)   
( 1827)          CASE ( 3 , 4 )    ! T pivot
( 1828)             iloc = jpiglo - 2 * ( nimpp - 1 )
( 1829) 
( 1830)             SELECT CASE ( cd_type1 )
( 1831) 
( 1832)             CASE ( 'T' , 'S', 'W' )
( 1833)                DO jk = 1, jpk
( 1834)                   DO ji = 2, nlci
( 1835)                      ijt=iloc-ji+2
( 1836)                      ptab1(ji,nlcj,jk) = psgn * ptab1(ijt,nlcj-2,jk)
( 1837)                   END DO
( 1838)                   DO ji = nlci/2+1, nlci
( 1839)                      ijt=iloc-ji+2
( 1840)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-1,jk)
( 1841)                   END DO
( 1842)                END DO
( 1843)           
( 1844)             CASE ( 'U' )
( 1845)                DO jk = 1, jpk
( 1846)                   DO ji = 1, nlci-1
( 1847)                      iju=iloc-ji+1
( 1848)                      ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-2,jk)
( 1849)                   END DO
( 1850)                   DO ji = nlci/2, nlci-1
( 1851)                      iju=iloc-ji+1
( 1852)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-1,jk)
( 1853)                   END DO
( 1854)                END DO
( 1855) 
( 1856)             CASE ( 'V' )
( 1857)                DO jk = 1, jpk
( 1858)                   DO ji = 2, nlci
( 1859)                      ijt=iloc-ji+2
( 1860)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-2,jk)
( 1861)                      ptab1(ji,nlcj  ,jk) = psgn * ptab1(ijt,nlcj-3,jk)
( 1862)                   END DO
( 1863)                END DO
( 1864) 
( 1865)             CASE ( 'F', 'G' )
( 1866)                DO jk = 1, jpk
( 1867)                   DO ji = 1, nlci-1
( 1868)                      iju=iloc-ji+1
( 1869)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-2,jk)
( 1870)                      ptab1(ji,nlcj  ,jk) = psgn * ptab1(iju,nlcj-3,jk)
( 1871)                   END DO
( 1872)                END DO
( 1873)   
( 1874)             END SELECT
( 1875)             
( 1876)             SELECT CASE ( cd_type2 )
( 1877) 
( 1878)             CASE ( 'T' , 'S', 'W' )
( 1879)                DO jk = 1, jpk
( 1880)                   DO ji = 2, nlci
( 1881)                      ijt=iloc-ji+2






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 131

( 1882)                      ptab2(ji,nlcj,jk) = psgn * ptab2(ijt,nlcj-2,jk)
( 1883)                   END DO
( 1884)                   DO ji = nlci/2+1, nlci
( 1885)                      ijt=iloc-ji+2
( 1886)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-1,jk)
( 1887)                   END DO
( 1888)                END DO
( 1889)           
( 1890)             CASE ( 'U' )
( 1891)                DO jk = 1, jpk
( 1892)                   DO ji = 1, nlci-1
( 1893)                      iju=iloc-ji+1
( 1894)                      ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-2,jk)
( 1895)                   END DO
( 1896)                   DO ji = nlci/2, nlci-1
( 1897)                      iju=iloc-ji+1
( 1898)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-1,jk)
( 1899)                   END DO
( 1900)                END DO
( 1901) 
( 1902)             CASE ( 'V' )
( 1903)                DO jk = 1, jpk
( 1904)                   DO ji = 2, nlci
( 1905)                      ijt=iloc-ji+2
( 1906)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-2,jk)
( 1907)                      ptab2(ji,nlcj  ,jk) = psgn * ptab2(ijt,nlcj-3,jk)
( 1908)                   END DO
( 1909)                END DO
( 1910) 
( 1911)             CASE ( 'F', 'G' )
( 1912)                DO jk = 1, jpk
( 1913)                   DO ji = 1, nlci-1
( 1914)                      iju=iloc-ji+1
( 1915)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-2,jk)
( 1916)                      ptab2(ji,nlcj  ,jk) = psgn * ptab2(iju,nlcj-3,jk)
( 1917)                   END DO
( 1918)                END DO
( 1919)   
( 1920)           END SELECT
( 1921)        
( 1922)          CASE ( 5 , 6 ) ! F pivot
( 1923)             iloc=jpiglo-2*(nimpp-1)
( 1924)   
( 1925)             SELECT CASE ( cd_type1 )
( 1926) 
( 1927)             CASE ( 'T' , 'S', 'W' )
( 1928)                DO jk = 1, jpk
( 1929)                   DO ji = 1, nlci
( 1930)                      ijt=iloc-ji+1
( 1931)                      ptab1(ji,nlcj,jk) = psgn * ptab1(ijt,nlcj-1,jk)
( 1932)                   END DO
( 1933)                END DO
( 1934) 
( 1935)             CASE ( 'U' )
( 1936)                DO jk = 1, jpk
( 1937)                   DO ji = 1, nlci-1
( 1938)                      iju=iloc-ji
( 1939)                      ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-1,jk)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 132

( 1940)                   END DO
( 1941)                END DO
( 1942) 
( 1943)             CASE ( 'V' )
( 1944)                DO jk = 1, jpk
( 1945)                   DO ji = 1, nlci
( 1946)                      ijt=iloc-ji+1
( 1947)                      ptab1(ji,nlcj  ,jk) = psgn * ptab1(ijt,nlcj-2,jk)
( 1948)                   END DO
( 1949)                   DO ji = nlci/2+1, nlci
( 1950)                      ijt=iloc-ji+1
( 1951)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(ijt,nlcj-1,jk)
( 1952)                   END DO
( 1953)                END DO
( 1954) 
( 1955)             CASE ( 'F', 'G' )
( 1956)                DO jk = 1, jpk
( 1957)                   DO ji = 1, nlci-1
( 1958)                      iju=iloc-ji
( 1959)                      ptab1(ji,nlcj,jk) = psgn * ptab1(iju,nlcj-2,jk)
( 1960)                   END DO
( 1961)                   DO ji = nlci/2+1, nlci-1
( 1962)                      iju=iloc-ji
( 1963)                      ptab1(ji,nlcj-1,jk) = psgn * ptab1(iju,nlcj-1,jk)
( 1964)                   END DO
( 1965)                END DO
( 1966)             END SELECT  ! cd_type1
( 1967) 
( 1968)             SELECT CASE ( cd_type2 )
( 1969) 
( 1970)             CASE ( 'T' , 'S', 'W' )
( 1971)                DO jk = 1, jpk
( 1972)                   DO ji = 1, nlci
( 1973)                      ijt=iloc-ji+1
( 1974)                      ptab2(ji,nlcj,jk) = psgn * ptab2(ijt,nlcj-1,jk)
( 1975)                   END DO
( 1976)                END DO
( 1977) 
( 1978)             CASE ( 'U' )
( 1979)                DO jk = 1, jpk
( 1980)                   DO ji = 1, nlci-1
( 1981)                      iju=iloc-ji
( 1982)                      ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-1,jk)
( 1983)                   END DO
( 1984)                END DO
( 1985) 
( 1986)             CASE ( 'V' )
( 1987)                DO jk = 1, jpk
( 1988)                   DO ji = 1, nlci
( 1989)                      ijt=iloc-ji+1
( 1990)                      ptab2(ji,nlcj  ,jk) = psgn * ptab2(ijt,nlcj-2,jk)
( 1991)                   END DO
( 1992)                   DO ji = nlci/2+1, nlci
( 1993)                      ijt=iloc-ji+1
( 1994)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(ijt,nlcj-1,jk)
( 1995)                   END DO
( 1996)                END DO
( 1997) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 133

( 1998)             CASE ( 'F', 'G' )
( 1999)                DO jk = 1, jpk
( 2000)                   DO ji = 1, nlci-1
( 2001)                      iju=iloc-ji
( 2002)                      ptab2(ji,nlcj,jk) = psgn * ptab2(iju,nlcj-2,jk)
( 2003)                   END DO
( 2004)                   DO ji = nlci/2+1, nlci-1
( 2005)                      iju=iloc-ji
( 2006)                      ptab2(ji,nlcj-1,jk) = psgn * ptab2(iju,nlcj-1,jk)
( 2007)                   END DO
( 2008)                END DO
( 2009) 
( 2010)             END SELECT  ! cd_type2
( 2011) 
( 2012)          END SELECT     !  npolj
( 2013)   
( 2014)       CASE DEFAULT ! more than 1 proc along I
( 2015)          IF ( npolj /= 0 ) THEN
( 2016)             CALL mpp_lbc_north (ptab1, cd_type1, psgn)  ! only for northern procs.
( 2017)             CALL mpp_lbc_north (ptab2, cd_type2, psgn)  ! only for northern procs.
( 2018)          ENDIF
( 2019) 
( 2020)       END SELECT ! jpni 
( 2021) 
( 2022) 
( 2023)       ! 5. East and west directions exchange
( 2024)       ! ------------------------------------
( 2025) 
( 2026)       SELECT CASE ( npolj )
( 2027) 
( 2028)       CASE ( 3, 4, 5, 6 )
( 2029) 
( 2030)          ! 5.1 Read Dirichlet lateral conditions
( 2031) 
( 2032)          SELECT CASE ( nbondi )
( 2033) 
( 2034)          CASE ( -1, 0, 1 )
( 2035)             iihom = nlci-nreci
( 2036)             DO jl = 1, jpreci
( 2037)                t4ew(:,jl,:,1,1) = ptab1(jpreci+jl,:,:)
( 2038)                t4we(:,jl,:,1,1) = ptab1(iihom +jl,:,:)
( 2039)                t4ew(:,jl,:,2,1) = ptab2(jpreci+jl,:,:)
( 2040)                t4we(:,jl,:,2,1) = ptab2(iihom +jl,:,:)
( 2041)             END DO
( 2042) 
( 2043)          END SELECT
( 2044) 
( 2045)          ! 5.2 Migrations
( 2046) 
( 2047) # 2066
( 2066)          !! MPI version
( 2067) 
( 2068)          imigr = jpreci * jpj * jpk * 2
( 2069)   
( 2070)          SELECT CASE ( nbondi )
( 2071)          CASE ( -1 )
( 2072)             CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req1 )
( 2073)             CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 134

( 2074)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2075)          CASE ( 0 )
( 2076)             CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
( 2077)             CALL mppsend( 2, t4we(1,1,1,1,1), imigr, noea, ml_req2 )
( 2078)             CALL mpprecv( 1, t4ew(1,1,1,1,2), imigr )
( 2079)             CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
( 2080)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2081)             IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2082)          CASE ( 1 )
( 2083)             CALL mppsend( 1, t4ew(1,1,1,1,1), imigr, nowe, ml_req1 )
( 2084)             CALL mpprecv( 2, t4we(1,1,1,1,2), imigr )
( 2085)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2086)          END SELECT
( 2087) 
( 2088) # 2089
( 2089)          ! 5.3 Write Dirichlet lateral conditions
( 2090) 
( 2091)          iihom = nlci-jpreci
( 2092) 
( 2093)          SELECT CASE ( nbondi)
( 2094)          CASE ( -1 )
( 2095)             DO jl = 1, jpreci
( 2096)                ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
( 2097)                ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
( 2098)             END DO
( 2099)          CASE ( 0 ) 
( 2100)             DO jl = 1, jpreci
( 2101)                ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
( 2102)                ptab1(iihom+jl,:,:) = t4ew(:,jl,:,1,2)
( 2103)                ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
( 2104)                ptab2(iihom+jl,:,:) = t4ew(:,jl,:,2,2)
( 2105)             END DO
( 2106)          CASE ( 1 )
( 2107)             DO jl = 1, jpreci
( 2108)                ptab1(jl      ,:,:) = t4we(:,jl,:,1,2)
( 2109)                ptab2(jl      ,:,:) = t4we(:,jl,:,2,2)
( 2110)             END DO
( 2111)          END SELECT
( 2112) 
( 2113)       END SELECT    ! npolj 
( 2114) 
( 2115)    END SUBROUTINE mpp_lnk_3d_gather
( 2116) 
( 2117) 
( 2118)    SUBROUTINE mpp_lnk_2d_e( pt2d, cd_type, psgn )
( 2119)       !!----------------------------------------------------------------------
( 2120)       !!                  ***  routine mpp_lnk_2d_e  ***
( 2121)       !!                  
( 2122)       !! ** Purpose :   Message passing manadgement for 2d array (with halo)
( 2123)       !!
( 2124)       !! ** Method  :   Use mppsend and mpprecv function for passing mask 
( 2125)       !!      between processors following neighboring subdomains.
( 2126)       !!            domain parameters
( 2127)       !!                    nlci   : first dimension of the local subdomain
( 2128)       !!                    nlcj   : second dimension of the local subdomain
( 2129)       !!                    jpr2di : number of rows for extra outer halo
( 2130)       !!                    jpr2dj : number of columns for extra outer halo
( 2131)       !!                    nbondi : mark for "east-west local boundary"






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 135

( 2132)       !!                    nbondj : mark for "north-south local boundary"
( 2133)       !!                    noea   : number for local neighboring processors 
( 2134)       !!                    nowe   : number for local neighboring processors
( 2135)       !!                    noso   : number for local neighboring processors
( 2136)       !!                    nono   : number for local neighboring processors
( 2137)       !!   
( 2138)       !! History :
( 2139)       !!       
( 2140)       !!   9.0  !  05-09  (R. Benshila, G. Madec)  original code
( 2141)       !!
( 2142)       !!----------------------------------------------------------------------
( 2143)       !! * Arguments
( 2144)       CHARACTER(len=1) , INTENT( in ) ::   &
( 2145)          cd_type       ! define the nature of pt2d array grid-points
( 2146)          !             !  = T , U , V , F , W 
( 2147)          !             !  = S : T-point, north fold treatment
( 2148)          !             !  = G : F-point, north fold treatment
( 2149)          !             !  = I : sea-ice velocity at F-point with index shift
( 2150)       REAL(wp), INTENT( in ) ::   &
( 2151)          psgn          ! control of the sign change
( 2152)          !             !   = -1. , the sign is changed if north fold boundary
( 2153)          !             !   =  1. , the sign is kept  if north fold boundary
( 2154)       REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj), INTENT( inout ) ::   &
( 2155)          pt2d          ! 2D array on which the boundary condition is applied
( 2156) 
( 2157)       !! * Local variables
( 2158)       INTEGER  ::   ji, jl      ! dummy loop indices
( 2159)       INTEGER  ::   &
( 2160)          imigr, iihom, ijhom,    &  ! temporary integers
( 2161)          iloc, ijt, iju             !    "          "
( 2162)       INTEGER  ::   &
( 2163)          ipreci, iprecj             ! temporary integers
( 2164)       INTEGER  ::   ml_req1, ml_req2, ml_err     ! for isend
( 2165)       INTEGER  ::   ml_stat(MPI_STATUS_SIZE)     ! for isend
( 2166)      !!---------------------------------------------------------------------
( 2167) 
( 2168)       ! take into account outer extra 2D overlap area
( 2169)       ipreci = jpreci + jpr2di
( 2170)       iprecj = jprecj + jpr2dj
( 2171) 
( 2172) 
( 2173)       ! 1. standard boundary treatment
( 2174)       ! ------------------------------
( 2175) 
( 2176)       !                                        ! East-West boundaries
( 2177)       !                                        ! ====================
( 2178)       IF( nbondi == 2 .AND.   &      ! Cyclic east-west
( 2179)          &    (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
( 2180)          pt2d(1-jpr2di:     1    ,:) = pt2d(jpim1-jpr2di:  jpim1 ,:)
( 2181)          pt2d(   jpi  :jpi+jpr2di,:) = pt2d(     2      :2+jpr2di,:)
( 2182) 
( 2183)       ELSE                           ! ... closed
( 2184)          SELECT CASE ( cd_type )
( 2185)          CASE ( 'T', 'U', 'V', 'W' , 'I' )
( 2186)             pt2d(  1-jpr2di   :jpreci    ,:) = 0.e0
( 2187)             pt2d(nlci-jpreci+1:jpi+jpr2di,:) = 0.e0
( 2188)          CASE ( 'F' )
( 2189)             pt2d(nlci-jpreci+1:jpi+jpr2di,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 136

( 2190)          END SELECT
( 2191)       ENDIF
( 2192) 
( 2193)       !                                        ! North-South boundaries
( 2194)       !                                        ! ======================
( 2195)       SELECT CASE ( cd_type )
( 2196)       CASE ( 'T', 'U', 'V', 'W' , 'I' )
( 2197)          pt2d(:,  1-jpr2dj   :  jprecj  ) = 0.e0
( 2198)          pt2d(:,nlcj-jprecj+1:jpj+jpr2dj) = 0.e0
( 2199)       CASE ( 'F' )
( 2200)          pt2d(:,nlcj-jprecj+1:jpj+jpr2dj) = 0.e0
( 2201)       END SELECT
( 2202) 
( 2203) 
( 2204)       ! 2. East and west directions
( 2205)       ! ---------------------------
( 2206) 
( 2207)       ! 2.1 Read Dirichlet lateral conditions
( 2208) 
( 2209)       SELECT CASE ( nbondi )
( 2210)       CASE ( -1, 0, 1 )    ! all except 2
( 2211)          iihom = nlci-nreci-jpr2di
( 2212)          DO jl = 1, ipreci
( 2213)             tr2ew(:,jl,1) = pt2d(jpreci+jl,:)
( 2214)             tr2we(:,jl,1) = pt2d(iihom +jl,:)
( 2215)          END DO
( 2216)       END SELECT
( 2217) 
( 2218)       ! 2.2 Migrations
( 2219) 
( 2220) # 2239
( 2239)       !! * MPI version
( 2240) 
( 2241)       imigr = ipreci * ( jpj + 2*jpr2dj)
( 2242) 
( 2243)       SELECT CASE ( nbondi )
( 2244)       CASE ( -1 )
( 2245)          CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req1 )
( 2246)          CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
( 2247)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2248)       CASE ( 0 )
( 2249)          CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
( 2250)          CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req2 )
( 2251)          CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
( 2252)          CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
( 2253)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2254)          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2255)       CASE ( 1 )
( 2256)          CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
( 2257)          CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
( 2258)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2259)       END SELECT
( 2260) 
( 2261) 
( 2262) # 2263
( 2263)       ! 2.3 Write Dirichlet lateral conditions
( 2264) 
( 2265)       iihom = nlci - jpreci






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 137

( 2266) 
( 2267)       SELECT CASE ( nbondi )
( 2268)       CASE ( -1 )
( 2269)          DO jl = 1, ipreci
( 2270)             pt2d(iihom+jl,:) = tr2ew(:,jl,2)
( 2271)          END DO
( 2272)       CASE ( 0 )
( 2273)          DO jl = 1, ipreci
( 2274)             pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
( 2275)             pt2d( iihom+jl,:) = tr2ew(:,jl,2)
( 2276)          END DO
( 2277)       CASE ( 1 )
( 2278)          DO jl = 1, ipreci
( 2279)             pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
( 2280)          END DO
( 2281)       END SELECT
( 2282) 
( 2283) 
( 2284)       ! 3. North and south directions
( 2285)       ! -----------------------------
( 2286) 
( 2287)       ! 3.1 Read Dirichlet lateral conditions
( 2288) 
( 2289)       IF( nbondj /= 2 ) THEN
( 2290)          ijhom = nlcj-nrecj-jpr2dj
( 2291)          DO jl = 1, iprecj
( 2292)             tr2sn(:,jl,1) = pt2d(:,ijhom +jl)
( 2293)             tr2ns(:,jl,1) = pt2d(:,jprecj+jl)
( 2294)          END DO
( 2295)       ENDIF
( 2296) 
( 2297)       ! 3.2 Migrations
( 2298) 
( 2299) # 2317
( 2317)       !! * MPI version
( 2318) 
( 2319)       imigr = iprecj * ( jpi + 2*jpr2di )
( 2320) 
( 2321)       SELECT CASE ( nbondj )
( 2322)       CASE ( -1 )
( 2323)          CALL mppsend( 4, tr2sn(1-jpr2di,1,1), imigr, nono, ml_req1 )
( 2324)          CALL mpprecv( 3, tr2ns(1-jpr2di,1,2), imigr )
( 2325)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2326)       CASE ( 0 )
( 2327)          CALL mppsend( 3, tr2ns(1-jpr2di,1,1), imigr, noso, ml_req1 )
( 2328)          CALL mppsend( 4, tr2sn(1-jpr2di,1,1), imigr, nono, ml_req2 )
( 2329)          CALL mpprecv( 3, tr2ns(1-jpr2di,1,2), imigr )
( 2330)          CALL mpprecv( 4, tr2sn(1-jpr2di,1,2), imigr )
( 2331)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2332)          IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2333)       CASE ( 1 )
( 2334)          CALL mppsend( 3, tr2ns(1-jpr2di,1,1), imigr, noso, ml_req1 )
( 2335)          CALL mpprecv( 4, tr2sn(1-jpr2di,1,2), imigr )
( 2336)          IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2337)       END SELECT
( 2338)   
( 2339) 
( 2340) # 2341






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 138

( 2341)       ! 3.3 Write Dirichlet lateral conditions
( 2342) 
( 2343)       ijhom = nlcj - jprecj  
( 2344) 
( 2345)       SELECT CASE ( nbondj )
( 2346)       CASE ( -1 )
( 2347)          DO jl = 1, iprecj
( 2348)             pt2d(:,ijhom+jl) = tr2ns(:,jl,2)
( 2349)          END DO
( 2350)       CASE ( 0 )
( 2351)          DO jl = 1, iprecj
( 2352)             pt2d(:,jl-jpr2dj) = tr2sn(:,jl,2)
( 2353)             pt2d(:,ijhom+jl ) = tr2ns(:,jl,2)
( 2354)          END DO
( 2355)       CASE ( 1 ) 
( 2356)          DO jl = 1, iprecj
( 2357)             pt2d(:,jl-jpr2dj) = tr2sn(:,jl,2)
( 2358)          END DO
( 2359)       END SELECT 
( 2360)   
( 2361) 
( 2362)       ! 4. north fold treatment
( 2363)       ! -----------------------
( 2364)   
( 2365)       ! 4.1 treatment without exchange (jpni odd)
( 2366)       
( 2367)       SELECT CASE ( jpni )
( 2368)   
( 2369)       CASE ( 1 ) ! only one proc along I, no mpp exchange
( 2370)   
( 2371)          SELECT CASE ( npolj )
( 2372)   
( 2373)          CASE ( 3 , 4 )   !  T pivot
( 2374)             iloc = jpiglo - 2 * ( nimpp - 1 )
( 2375)   
( 2376)             SELECT CASE ( cd_type )
( 2377)   
( 2378)             CASE ( 'T', 'S', 'W' )
( 2379)                DO jl = 0, iprecj-1
( 2380)                   DO ji = 2-jpr2di, nlci+jpr2di
( 2381)                      ijt=iloc-ji+2
( 2382)                      pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-2-jl)
( 2383)                   END DO
( 2384)                END DO
( 2385)                DO ji = nlci/2+1, nlci+jpr2di
( 2386)                   ijt=iloc-ji+2
( 2387)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
( 2388)                END DO
( 2389)  
( 2390)             CASE ( 'U' )
( 2391)                DO jl =0, iprecj-1
( 2392)                   DO ji = 1-jpr2di, nlci-1-jpr2di
( 2393)                      iju=iloc-ji+1
( 2394)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-2-jl)
( 2395)                   END DO
( 2396)                END DO
( 2397)                DO ji = nlci/2, nlci-1+jpr2di
( 2398)                   iju=iloc-ji+1






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 139

( 2399)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
( 2400)                END DO
( 2401)   
( 2402)             CASE ( 'V' )
( 2403)                DO jl = -1, iprecj-1
( 2404)                   DO ji = 2-jpr2di, nlci+jpr2di
( 2405)                      ijt=iloc-ji+2
( 2406)                      pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-3-jl)
( 2407)                   END DO
( 2408)                END DO
( 2409)   
( 2410)             CASE ( 'F', 'G' )
( 2411)                DO jl = -1, iprecj-1
( 2412)                   DO ji = 1-jpr2di, nlci-1+jpr2di
( 2413)                      iju=iloc-ji+1
( 2414)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-3-jl)
( 2415)                   END DO
( 2416)                END DO
( 2417)   
( 2418)             CASE ( 'I' )                                  ! ice U-V point
( 2419)                DO jl = 0, iprecj-1
( 2420)                   pt2d(2,nlcj+jl) = psgn * pt2d(3,nlcj-1-jl)
( 2421)                   DO ji = 3, nlci+jpr2di
( 2422)                      iju = iloc - ji + 3
( 2423)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-1-jl)
( 2424)                   END DO
( 2425)                END DO
( 2426)   
( 2427)             END SELECT
( 2428)   
( 2429)          CASE ( 5 , 6 )                 ! F pivot
( 2430)             iloc=jpiglo-2*(nimpp-1)
( 2431)   
( 2432)             SELECT CASE (cd_type )
( 2433)   
( 2434)             CASE ( 'T', 'S', 'W' )
( 2435)                DO jl = 0, iprecj-1
( 2436)                   DO ji = 1-jpr2di, nlci+jpr2di
( 2437)                      ijt=iloc-ji+1
( 2438)                      pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-1-jl)
( 2439)                   END DO
( 2440)                END DO
( 2441)   
( 2442)             CASE ( 'U' )
( 2443)                DO jl = 0, iprecj-1
( 2444)                   DO ji = 1-jpr2di, nlci-1+jpr2di
( 2445)                      iju=iloc-ji
( 2446)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-1-jl)
( 2447)                   END DO
( 2448)                END DO
( 2449)  
( 2450)             CASE ( 'V' )
( 2451)                DO jl = 0, iprecj-1
( 2452)                   DO ji = 1-jpr2di, nlci+jpr2di
( 2453)                      ijt=iloc-ji+1
( 2454)                      pt2d(ji,nlcj+jl) = psgn * pt2d(ijt,nlcj-2-jl)
( 2455)                   END DO
( 2456)                END DO 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 140

( 2457)                DO ji = nlci/2+1, nlci+jpr2di
( 2458)                   ijt=iloc-ji+1
( 2459)                   pt2d(ji,nlcj-1) = psgn * pt2d(ijt,nlcj-1)
( 2460)                END DO
( 2461)   
( 2462)             CASE ( 'F', 'G' )
( 2463)                DO jl = 0, iprecj-1
( 2464)                   DO ji = 1-jpr2di, nlci-1+jpr2di
( 2465)                      iju=iloc-ji
( 2466)                      pt2d(ji,nlcj+jl) = psgn * pt2d(iju,nlcj-2-jl)
( 2467)                   END DO
( 2468)                END DO
( 2469)                DO ji = nlci/2+1, nlci-1+jpr2di
( 2470)                   iju=iloc-ji
( 2471)                   pt2d(ji,nlcj-1) = psgn * pt2d(iju,nlcj-1)
( 2472)                END DO
( 2473)   
( 2474)             CASE ( 'I' )                                  ! ice U-V point
( 2475)                pt2d( 2 ,nlcj) = 0.e0
( 2476)                DO jl = 0, iprecj-1
( 2477)                   DO ji = 2 , nlci-1+jpr2di
( 2478)                      ijt = iloc - ji + 2
( 2479)                      pt2d(ji,nlcj+jl)= 0.5 * ( pt2d(ji,nlcj-1-jl) + psgn * pt2d(ijt,nlcj-1-jl) )
( 2480)                   END DO
( 2481)                END DO
( 2482)   
( 2483)             END SELECT   ! cd_type
( 2484)   
( 2485)          END SELECT   ! npolj
( 2486) 
( 2487)       CASE DEFAULT   ! more than 1 proc along I
( 2488)          IF( npolj /= 0 )   CALL mpp_lbc_north_e( pt2d, cd_type, psgn )   ! only for northern procs
( 2489)          
( 2490)       END SELECT   ! jpni
( 2491) 
( 2492) 
( 2493)       ! 5. East and west directions
( 2494)       ! ---------------------------
( 2495) 
( 2496)       SELECT CASE ( npolj )
( 2497) 
( 2498)       CASE ( 3, 4, 5, 6 )
( 2499) 
( 2500)          ! 5.1 Read Dirichlet lateral conditions
( 2501) 
( 2502)          SELECT CASE ( nbondi )
( 2503)          CASE ( -1, 0, 1 )
( 2504)             iihom = nlci-nreci-jpr2di
( 2505)             DO jl = 1, ipreci
( 2506)                tr2ew(:,jl,1) = pt2d(jpreci+jl,:)
( 2507)                tr2we(:,jl,1) = pt2d(iihom +jl,:)
( 2508)             END DO
( 2509)          END SELECT
( 2510) 
( 2511)          ! 5.2 Migrations
( 2512) 
( 2513) # 2532
( 2532)          !! * MPI version






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 141

( 2533)   
( 2534)          imigr=ipreci* ( jpj + 2*jpr2dj )
( 2535)   
( 2536)          SELECT CASE ( nbondi )
( 2537)          CASE ( -1 )
( 2538)             CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req1 )
( 2539)             CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
( 2540)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2541)          CASE ( 0 )
( 2542)             CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
( 2543)             CALL mppsend( 2, tr2we(1-jpr2dj,1,1), imigr, noea, ml_req2 )
( 2544)             CALL mpprecv( 1, tr2ew(1-jpr2dj,1,2), imigr )
( 2545)             CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
( 2546)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2547)             IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2548)          CASE ( 1 )
( 2549)             CALL mppsend( 1, tr2ew(1-jpr2dj,1,1), imigr, nowe, ml_req1 )
( 2550)             CALL mpprecv( 2, tr2we(1-jpr2dj,1,2), imigr )
( 2551)             IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2552)          END SELECT 
( 2553) 
( 2554) # 2555
( 2555)          ! 5.3 Write Dirichlet lateral conditions
( 2556)   
( 2557)          iihom = nlci - jpreci
( 2558)   
( 2559)          SELECT CASE ( nbondi )
( 2560)          CASE ( -1 )
( 2561)             DO jl = 1, ipreci
( 2562)                pt2d(iihom+jl,:) = tr2ew(:,jl,2)
( 2563)             END DO
( 2564)          CASE ( 0 )
( 2565)             DO jl = 1, ipreci
( 2566)                pt2d(jl- jpr2di,:) = tr2we(:,jl,2)
( 2567)                pt2d(iihom+jl,:) = tr2ew(:,jl,2)
( 2568)             END DO
( 2569)          CASE ( 1 )
( 2570)             DO jl = 1, ipreci
( 2571)                pt2d(jl-jpr2di,:) = tr2we(:,jl,2)
( 2572)             END DO
( 2573)          END SELECT 
( 2574)   
( 2575)       END SELECT   ! npolj
( 2576)   
( 2577)    END SUBROUTINE mpp_lnk_2d_e
( 2578) 
( 2579) 
( 2580)    SUBROUTINE mpplnks( ptab )
( 2581)       !!----------------------------------------------------------------------
( 2582)       !!                  ***  routine mpplnks  ***
( 2583)       !!
( 2584)       !! ** Purpose :   Message passing manadgement for add 2d array local boundary
( 2585)       !!
( 2586)       !! ** Method  :   Use mppsend and mpprecv function for passing mask between
( 2587)       !!       processors following neighboring subdomains.
( 2588)       !!            domain parameters
( 2589)       !!                    nlci   : first dimension of the local subdomain
( 2590)       !!                    nlcj   : second dimension of the local subdomain






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 142

( 2591)       !!                    nbondi : mark for "east-west local boundary"
( 2592)       !!                    nbondj : mark for "north-south local boundary"
( 2593)       !!                    noea   : number for local neighboring processors 
( 2594)       !!                    nowe   : number for local neighboring processors
( 2595)       !!                    noso   : number for local neighboring processors
( 2596)       !!                    nono   : number for local neighboring processors
( 2597)       !!
( 2598)       !!----------------------------------------------------------------------
( 2599)       !! * Arguments
( 2600)       REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   &
( 2601)          ptab                     ! 2D array
( 2602)   
( 2603)       !! * Local variables
( 2604)       INTEGER ::   ji, jl         ! dummy loop indices
( 2605)       INTEGER ::   &
( 2606)          imigr, iihom, ijhom      ! temporary integers
( 2607)       INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
( 2608)       INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
( 2609)       !!----------------------------------------------------------------------
( 2610) 
( 2611) 
( 2612)       ! 1. north fold treatment
( 2613)       ! -----------------------
( 2614) 
( 2615)       ! 1.1 treatment without exchange (jpni odd)
( 2616)   
( 2617)       SELECT CASE ( npolj )
( 2618)       CASE ( 4 )
( 2619)          DO ji = 1, nlci
( 2620)             ptab(ji,nlcj-2) = ptab(ji,nlcj-2) + t2p1(ji,1,1)
( 2621)          END DO
( 2622)       CASE ( 6 )
( 2623)          DO ji = 1, nlci
( 2624)             ptab(ji,nlcj-1) = ptab(ji,nlcj-1) + t2p1(ji,1,1)
( 2625)          END DO
( 2626) 
( 2627)       ! 1.2 treatment with exchange (jpni greater than 1)
( 2628)       ! 
( 2629)       CASE ( 3 )
( 2630) # 2641
( 2641)        !! * MPI version
( 2642) 
( 2643)        imigr=jprecj*jpi
( 2644) 
( 2645)        CALL mppsend(3,t2p1(1,1,1),imigr,nono, ml_req1)
( 2646)        CALL mpprecv(3,t2p1(1,1,2),imigr)
( 2647)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2648) 
( 2649) 
( 2650) # 2651
( 2651)        ! Write north fold conditions
( 2652) 
( 2653)        DO ji = 1, nlci
( 2654)           ptab(ji,nlcj-2) = ptab(ji,nlcj-2)+t2p1(ji,1,2)
( 2655)        END DO
( 2656) 
( 2657)     CASE ( 5 )
( 2658) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 143

( 2659) # 2670
( 2670)        !! * Local variables   (MPI version)
( 2671) 
( 2672)        imigr=jprecj*jpi
( 2673) 
( 2674)        CALL mppsend(3,t2p1(1,1,1),imigr,nono, ml_req1)
( 2675)        CALL mpprecv(3,t2p1(1,1,2),imigr)
( 2676)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2677) 
( 2678) 
( 2679) # 2680
( 2680)        ! Write north fold conditions
( 2681) 
( 2682)        DO ji = 1, nlci
( 2683)           ptab(ji,nlcj-1) = ptab(ji,nlcj-1)+t2p1(ji,1,2)
( 2684)        END DO
( 2685) 
( 2686)     END SELECT
( 2687) 
( 2688) 
( 2689)     ! 2. East and west directions
( 2690)     ! ---------------------------
( 2691) 
( 2692)     ! 2.1 Read Dirichlet lateral conditions
( 2693) 
( 2694)     iihom = nlci-jpreci
( 2695) 
( 2696)     SELECT CASE ( nbondi )
( 2697) 
( 2698)     CASE ( -1, 0, 1 )  ! all except 2
( 2699)        DO jl = 1, jpreci
( 2700)              t2ew(:,jl,1) = ptab(  jl    ,:)
( 2701)              t2we(:,jl,1) = ptab(iihom+jl,:)
( 2702)        END DO
( 2703)     END SELECT
( 2704) 
( 2705)     ! 2.2 Migrations
( 2706) 
( 2707) # 2730
( 2730)     !! * Local variables   (MPI version)
( 2731) 
( 2732)     imigr=jpreci*jpj
( 2733) 
( 2734)     SELECT CASE ( nbondi )
( 2735) 
( 2736)     CASE ( -1 )
( 2737)        CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req1)
( 2738)        CALL mpprecv(1,t2ew(1,1,2),imigr)
( 2739)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2740)     CASE ( 0 )
( 2741)        CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
( 2742)        CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req2)
( 2743)        CALL mpprecv(1,t2ew(1,1,2),imigr)
( 2744)        CALL mpprecv(2,t2we(1,1,2),imigr)
( 2745)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2746)        IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2747) 
( 2748)     CASE ( 1 )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 144

( 2749)        CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
( 2750)        CALL mpprecv(2,t2we(1,1,2),imigr)
( 2751)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2752) 
( 2753)     END SELECT
( 2754) 
( 2755) 
( 2756) # 2757
( 2757)     ! 2.3 Write Dirichlet lateral conditions
( 2758) 
( 2759)        iihom = nlci-nreci
( 2760) 
( 2761)     SELECT CASE ( nbondi )
( 2762) 
( 2763)     CASE ( -1 )
( 2764)        DO jl = 1, jpreci
( 2765)              ptab(iihom +jl,:) = ptab(iihom +jl,:)+t2ew(:,jl,2)
( 2766)        END DO
( 2767) 
( 2768)     CASE ( 0 )
( 2769)        DO jl = 1, jpreci
( 2770)              ptab(jpreci+jl,:) = ptab(jpreci+jl,:)+t2we(:,jl,2)
( 2771)              ptab(iihom +jl,:) = ptab(iihom +jl,:)+t2ew(:,jl,2)
( 2772)        END DO
( 2773) 
( 2774)     CASE ( 1 )
( 2775)        DO jl = 1, jpreci
( 2776)              ptab(jpreci+jl,:) = ptab(jpreci+jl,:)+t2we(:,jl,2)
( 2777)        END DO
( 2778)     END SELECT
( 2779) 
( 2780) 
( 2781)     ! 3. North and south directions
( 2782)     ! -----------------------------
( 2783) 
( 2784)     ! 3.1 Read Dirichlet lateral conditions
( 2785) 
( 2786)     ijhom = nlcj-jprecj
( 2787) 
( 2788)     SELECT CASE ( nbondj )
( 2789) 
( 2790)     CASE ( -1, 0, 1 )
( 2791)        DO jl = 1, jprecj
( 2792)              t2sn(:,jl,1) = ptab(:,ijhom+jl)
( 2793)              t2ns(:,jl,1) = ptab(:,   jl   )
( 2794)        END DO
( 2795) 
( 2796)     END SELECT 
( 2797) 
( 2798)     ! 3.2 Migrations
( 2799) 
( 2800) # 2823
( 2823)     !! * Local variables   (MPI version)
( 2824) 
( 2825)     imigr=jprecj*jpi
( 2826) 
( 2827)     SELECT CASE ( nbondj )
( 2828) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 145

( 2829)     CASE ( -1 )
( 2830)        CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req1)
( 2831)        CALL mpprecv(3,t2ns(1,1,2),imigr)
( 2832)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2833) 
( 2834)     CASE ( 0 )
( 2835)        CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
( 2836)        CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req2)
( 2837)        CALL mpprecv(3,t2ns(1,1,2),imigr)
( 2838)        CALL mpprecv(4,t2sn(1,1,2),imigr)
( 2839)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2840)        IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 2841) 
( 2842)     CASE ( 1 )
( 2843)        CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
( 2844)        CALL mpprecv(4,t2sn(1,1,2),imigr)
( 2845)        IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 2846)     END SELECT
( 2847) 
( 2848) 
( 2849) # 2850
( 2850)     ! 3.3 Write Dirichlet lateral conditions
( 2851) 
( 2852)        ijhom = nlcj-nrecj
( 2853) 
( 2854)     SELECT CASE ( nbondj )
( 2855) 
( 2856)     CASE ( -1 )
( 2857)        DO jl = 1, jprecj
( 2858)              ptab(:,ijhom +jl) = ptab(:,ijhom +jl)+t2ns(:,jl,2)
( 2859)        END DO
( 2860) 
( 2861)     CASE ( 0 )
( 2862)        DO jl = 1, jprecj
( 2863)              ptab(:,jprecj+jl) = ptab(:,jprecj+jl)+t2sn(:,jl,2)
( 2864)              ptab(:,ijhom +jl) = ptab(:,ijhom +jl)+t2ns(:,jl,2)
( 2865)        END DO
( 2866) 
( 2867)     CASE ( 1 ) 
( 2868)        DO jl = 1, jprecj
( 2869)              ptab(:,jprecj+jl) = ptab(:,jprecj+jl)+t2sn(:,jl,2)
( 2870)        END DO
( 2871) 
( 2872)     END SELECT
( 2873) 
( 2874)   END SUBROUTINE mpplnks
( 2875) 
( 2876) 
( 2877)    SUBROUTINE mppsend( ktyp, pmess, kbytes, kdest, md_req)
( 2878)       !!----------------------------------------------------------------------
( 2879)       !!                  ***  routine mppsend  ***
( 2880)       !!                   
( 2881)       !! ** Purpose :   Send messag passing array
( 2882)       !!
( 2883)       !!----------------------------------------------------------------------
( 2884)       !! * Arguments
( 2885)       REAL(wp), INTENT(inout) ::   pmess(*)       ! array of real
( 2886)       INTEGER , INTENT( in  ) ::   kbytes,     &  ! size of the array pmess






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 146

( 2887)          &                         kdest ,     &  ! receive process number
( 2888)          &                         ktyp,       &  ! Tag of the message
( 2889)          &                         md_req         ! Argument for isend
( 2890)       !!----------------------------------------------------------------------
( 2891) # 2895
( 2895)       !! * MPI version
( 2896)       INTEGER ::   iflag
( 2897) 
( 2898)       SELECT CASE ( c_mpi_send )
( 2899)       CASE ( 'S' )                ! Standard mpi send (blocking)
( 2900)          CALL mpi_send ( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
( 2901)             &                          mpi_comm_opa, iflag )
( 2902)       CASE ( 'B' )                ! Buffer mpi send (blocking)
( 2903)          CALL mpi_bsend( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
( 2904)             &                          mpi_comm_opa, iflag )
( 2905)       CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
( 2906)          ! Be carefull, one more argument here : the mpi request identifier..
( 2907)          CALL mpi_isend( pmess, kbytes, mpi_double_precision, kdest, ktyp,   &
( 2908)             &                          mpi_comm_opa, md_req, iflag )
( 2909)       END SELECT
( 2910) 
( 2911) # 2912
( 2912)    END SUBROUTINE mppsend
( 2913) 
( 2914) 
( 2915)    SUBROUTINE mpprecv( ktyp, pmess, kbytes )
( 2916)       !!----------------------------------------------------------------------
( 2917)       !!                  ***  routine mpprecv  ***
( 2918)       !!
( 2919)       !! ** Purpose :   Receive messag passing array
( 2920)       !!
( 2921)       !!----------------------------------------------------------------------
( 2922)       !! * Arguments
( 2923)       REAL(wp), INTENT(inout) ::   pmess(*)       ! array of real
( 2924)       INTEGER , INTENT( in  ) ::   kbytes,     &  ! suze of the array pmess
( 2925)          &                         ktyp           ! Tag of the recevied message
( 2926)       !!----------------------------------------------------------------------
( 2927) # 2931
( 2931)       !! * MPI version
( 2932)       INTEGER :: istatus(mpi_status_size)
( 2933)       INTEGER :: iflag
( 2934) 
( 2935)       CALL mpi_recv( pmess, kbytes, mpi_double_precision, mpi_any_source, ktyp,   &
( 2936)          &                          mpi_comm_opa, istatus, iflag )
( 2937) 
( 2938) # 2939
( 2939)    END SUBROUTINE mpprecv
( 2940) 
( 2941) 
( 2942)    SUBROUTINE mppgather( ptab, kp, pio )
( 2943)       !!----------------------------------------------------------------------
( 2944)       !!                   ***  routine mppgather  ***
( 2945)       !!                   
( 2946)       !! ** Purpose :   Transfert between a local subdomain array and a work 
( 2947)       !!     array which is distributed following the vertical level.
( 2948)       !!
( 2949)       !! ** Method  :
( 2950)       !!






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 147

( 2951)       !!----------------------------------------------------------------------
( 2952)       !! * Arguments
( 2953)       REAL(wp), DIMENSION(jpi,jpj),       INTENT( in  ) ::   ptab   ! subdomain input array
( 2954)       INTEGER ,                           INTENT( in  ) ::   kp     ! record length
( 2955)       REAL(wp), DIMENSION(jpi,jpj,jpnij), INTENT( out ) ::   pio    ! subdomain input array
( 2956)       !!---------------------------------------------------------------------
( 2957) # 2965
( 2965)       !! * Local variables   (MPI version)
( 2966)       INTEGER :: itaille,ierror
( 2967)   
( 2968)       itaille=jpi*jpj
( 2969)       CALL mpi_gather( ptab, itaille, mpi_double_precision, pio, itaille,   &
( 2970)          &                            mpi_double_precision, kp , mpi_comm_opa, ierror ) 
( 2971) 
( 2972) # 2973
( 2973)    END SUBROUTINE mppgather
( 2974) 
( 2975) 
( 2976)    SUBROUTINE mppscatter( pio, kp, ptab )
( 2977)       !!----------------------------------------------------------------------
( 2978)       !!                  ***  routine mppscatter  ***
( 2979)       !!
( 2980)       !! ** Purpose :   Transfert between awork array which is distributed 
( 2981)       !!      following the vertical level and the local subdomain array.
( 2982)       !!
( 2983)       !! ** Method :
( 2984)       !!
( 2985)       !!----------------------------------------------------------------------
( 2986)       REAL(wp), DIMENSION(jpi,jpj,jpnij)  ::  pio        ! output array
( 2987)       INTEGER                             ::   kp        ! Tag (not used with MPI
( 2988)       REAL(wp), DIMENSION(jpi,jpj)        ::  ptab       ! subdomain array input
( 2989)       !!---------------------------------------------------------------------
( 2990) # 2998
( 2998)       !! * Local variables   (MPI version)
( 2999)       INTEGER :: itaille, ierror
( 3000)   
( 3001)       itaille=jpi*jpj
( 3002) 
( 3003)       CALL mpi_scatter( pio, itaille, mpi_double_precision, ptab, itaille,   &
( 3004)          &                            mpi_double_precision, kp, mpi_comm_opa, ierror )
( 3005) 
( 3006) # 3007
( 3007)    END SUBROUTINE mppscatter
( 3008) 
( 3009) 
( 3010)    SUBROUTINE mppisl_a_int( ktab, kdim )
( 3011)       !!----------------------------------------------------------------------
( 3012)       !!                  ***  routine mppisl_a_int  ***
( 3013)       !!                   
( 3014)       !! ** Purpose :   Massively parallel processors
( 3015)       !!                Find the  non zero value
( 3016)       !!
( 3017)       !!----------------------------------------------------------------------
( 3018)       !! * Arguments
( 3019)       INTEGER, INTENT( in  )                  ::   kdim       ! ???
( 3020)       INTEGER, INTENT(inout), DIMENSION(kdim) ::   ktab       ! ???
( 3021)   
( 3022) # 3057






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 148

( 3057)       !! * Local variables   (MPI version)
( 3058)       LOGICAL  :: lcommute
( 3059)       INTEGER, DIMENSION(kdim) ::   iwork
( 3060)       INTEGER  :: mpi_isl,ierror
( 3061)   
( 3062)       lcommute = .TRUE.
( 3063)       CALL mpi_op_create( lc_isl, lcommute, mpi_isl, ierror )
( 3064)       CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer   &
( 3065)            , mpi_isl, mpi_comm_opa, ierror )
( 3066)       ktab(:) = iwork(:)
( 3067) 
( 3068) # 3069
( 3069)    END SUBROUTINE mppisl_a_int
( 3070) 
( 3071) 
( 3072)    SUBROUTINE mppisl_int( ktab )
( 3073)       !!----------------------------------------------------------------------
( 3074)       !!                  ***  routine mppisl_int  ***
( 3075)       !!                   
( 3076)       !! ** Purpose :   Massively parallel processors
( 3077)       !!                Find the non zero value
( 3078)       !!
( 3079)       !!----------------------------------------------------------------------
( 3080)       !! * Arguments
( 3081)       INTEGER , INTENT( inout ) ::   ktab        ! 
( 3082) 
( 3083) # 3110
( 3110)   
( 3111)       !! * Local variables   (MPI version)
( 3112)       LOGICAL :: lcommute
( 3113)       INTEGER :: mpi_isl,ierror
( 3114)       INTEGER ::   iwork
( 3115)   
( 3116)       lcommute = .TRUE.
( 3117)       CALL mpi_op_create(lc_isl,lcommute,mpi_isl,ierror)
( 3118)       CALL mpi_allreduce(ktab, iwork, 1,mpi_integer   &
( 3119)            ,mpi_isl,mpi_comm_opa,ierror)
( 3120)       ktab = iwork
( 3121) 
( 3122) # 3123
( 3123)    END SUBROUTINE mppisl_int
( 3124) 
( 3125) 
( 3126)    SUBROUTINE mppmax_a_int( ktab, kdim )
( 3127)       !!----------------------------------------------------------------------
( 3128)       !!                  ***  routine mppmax_a_int  ***
( 3129)       !! 
( 3130)       !! ** Purpose :   Find maximum value in an integer layout array
( 3131)       !!
( 3132)       !!----------------------------------------------------------------------
( 3133)       !! * Arguments
( 3134)       INTEGER , INTENT( in  )                  ::   kdim        ! size of array
( 3135)       INTEGER , INTENT(inout), DIMENSION(kdim) ::   ktab        ! input array
( 3136)   
( 3137) # 3164
( 3164)   
( 3165)       !! * Local variables   (MPI version)
( 3166)       INTEGER :: ierror






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 149

( 3167)       INTEGER, DIMENSION(kdim) ::   iwork
( 3168)   
( 3169)       CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer,   &
( 3170)            &                mpi_max, mpi_comm_opa, ierror )
( 3171)   
( 3172)       ktab(:) = iwork(:)
( 3173) 
( 3174) # 3175
( 3175)    END SUBROUTINE mppmax_a_int
( 3176) 
( 3177) 
( 3178)    SUBROUTINE mppmax_int( ktab )
( 3179)       !!----------------------------------------------------------------------
( 3180)       !!                  ***  routine mppmax_int  ***
( 3181)       !!
( 3182)       !! ** Purpose :
( 3183)       !!     Massively parallel processors
( 3184)       !!     Find maximum value in an integer layout array
( 3185)       !!
( 3186)       !!----------------------------------------------------------------------
( 3187)       !! * Arguments
( 3188)       INTEGER, INTENT(inout) ::   ktab      ! ???
( 3189)   
( 3190)       !! * Local declarations
( 3191) 
( 3192) 
( 3193) # 3214
( 3214)       !! * Local variables   (MPI version)
( 3215)       INTEGER ::  ierror, iwork
( 3216)   
( 3217)       CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
( 3218)            &              ,mpi_max,mpi_comm_opa,ierror)
( 3219)   
( 3220)       ktab = iwork
( 3221) 
( 3222) # 3223
( 3223)    END SUBROUTINE mppmax_int
( 3224) 
( 3225) 
( 3226)    SUBROUTINE mppmin_a_int( ktab, kdim )
( 3227)       !!----------------------------------------------------------------------
( 3228)       !!                  ***  routine mppmin_a_int  ***
( 3229)       !! 
( 3230)       !! ** Purpose :   Find minimum value in an integer layout array
( 3231)       !!
( 3232)       !!----------------------------------------------------------------------
( 3233)       !! * Arguments
( 3234)       INTEGER , INTENT( in  )                  ::   kdim        ! size of array
( 3235)       INTEGER , INTENT(inout), DIMENSION(kdim) ::   ktab        ! input array
( 3236)   
( 3237) # 3264
( 3264)   
( 3265)       !! * Local variables   (MPI version)
( 3266)       INTEGER :: ierror
( 3267)       INTEGER, DIMENSION(kdim) ::   iwork
( 3268)   
( 3269)       CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer,   &
( 3270)            &                mpi_min, mpi_comm_opa, ierror )






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 150

( 3271)   
( 3272)       ktab(:) = iwork(:)
( 3273) 
( 3274) # 3275
( 3275)    END SUBROUTINE mppmin_a_int
( 3276) 
( 3277) 
( 3278)    SUBROUTINE mppmin_int( ktab )
( 3279)       !!----------------------------------------------------------------------
( 3280)       !!                  ***  routine mppmin_int  ***
( 3281)       !!
( 3282)       !! ** Purpose :
( 3283)       !!     Massively parallel processors
( 3284)       !!     Find minimum value in an integer layout array
( 3285)       !!
( 3286)       !!----------------------------------------------------------------------
( 3287)       !! * Arguments
( 3288)       INTEGER, INTENT(inout) ::   ktab      ! ???
( 3289)   
( 3290)       !! * Local declarations
( 3291) 
( 3292) 
( 3293) # 3314
( 3314)       !! * Local variables   (MPI version)
( 3315)       INTEGER ::  ierror, iwork
( 3316)   
( 3317)       CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
( 3318)            &              ,mpi_min,mpi_comm_opa,ierror)
( 3319)   
( 3320)       ktab = iwork
( 3321) 
( 3322) # 3323
( 3323)    END SUBROUTINE mppmin_int
( 3324) 
( 3325) 
( 3326)    SUBROUTINE mppsum_a_int( ktab, kdim )
( 3327)       !!----------------------------------------------------------------------
( 3328)       !!                  ***  routine mppsum_a_int  ***
( 3329)       !!                    
( 3330)       !! ** Purpose :   Massively parallel processors
( 3331)       !!                Global integer sum
( 3332)       !!
( 3333)       !!----------------------------------------------------------------------
( 3334)       !! * Arguments
( 3335)       INTEGER, INTENT( in  )                   ::   kdim      ! ???
( 3336)       INTEGER, INTENT(inout), DIMENSION (kdim) ::   ktab      ! ???
( 3337)   
( 3338) 
( 3339) # 3367
( 3367)       !! * Local variables   (MPI version)
( 3368)       INTEGER :: ierror
( 3369)       INTEGER, DIMENSION (kdim) ::  iwork
( 3370)   
( 3371)       CALL mpi_allreduce(ktab, iwork,kdim,mpi_integer   &
( 3372)            ,mpi_sum,mpi_comm_opa,ierror)
( 3373)   
( 3374)       ktab(:) = iwork(:)
( 3375) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 151

( 3376) # 3377
( 3377)    END SUBROUTINE mppsum_a_int
( 3378) 
( 3379) 
( 3380)   SUBROUTINE mppsum_int( ktab )
( 3381)     !!----------------------------------------------------------------------
( 3382)     !!                 ***  routine mppsum_int  ***
( 3383)     !!                  
( 3384)     !! ** Purpose :   Global integer sum
( 3385)     !!
( 3386)     !!----------------------------------------------------------------------
( 3387)     !! * Arguments
( 3388)     INTEGER, INTENT(inout) ::   ktab
( 3389) 
( 3390) 
( 3391) # 3411
( 3411)     !! * Local variables   (MPI version)
( 3412)     INTEGER :: ierror, iwork
( 3413) 
( 3414)     CALL mpi_allreduce(ktab,iwork, 1,mpi_integer   &
( 3415)          ,mpi_sum,mpi_comm_opa,ierror)
( 3416) 
( 3417)     ktab = iwork
( 3418) 
( 3419) 
( 3420) # 3421
( 3421)   END SUBROUTINE mppsum_int
( 3422) 
( 3423) 
( 3424)   SUBROUTINE mppisl_a_real( ptab, kdim )
( 3425)     !!----------------------------------------------------------------------
( 3426)     !!                 ***  routine mppisl_a_real  ***
( 3427)     !!         
( 3428)     !! ** Purpose :   Massively parallel processors
( 3429)     !!           Find the non zero island barotropic stream function value
( 3430)     !!
( 3431)     !!   Modifications:
( 3432)     !!        !  93-09 (M. Imbard)
( 3433)     !!        !  96-05 (j. Escobar)
( 3434)     !!        !  98-05 (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI 
( 3435)     !!----------------------------------------------------------------------
( 3436)     INTEGER , INTENT( in  )                  ::   kdim      ! ???
( 3437)     REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab      ! ???
( 3438) 
( 3439) 
( 3440) # 3476
( 3476)     !! * Local variables   (MPI version)
( 3477)     LOGICAL ::   lcommute = .TRUE.
( 3478)     INTEGER ::   mpi_isl, ierror
( 3479)     REAL(wp), DIMENSION(kdim) ::  zwork
( 3480) 
( 3481)     CALL mpi_op_create(lc_isl,lcommute,mpi_isl,ierror)
( 3482)     CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
( 3483)          ,mpi_isl,mpi_comm_opa,ierror)
( 3484)     ptab(:) = zwork(:)
( 3485) 
( 3486) 
( 3487) # 3488






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 152

( 3488)   END SUBROUTINE mppisl_a_real
( 3489) 
( 3490) 
( 3491)    SUBROUTINE mppisl_real( ptab )
( 3492)       !!----------------------------------------------------------------------
( 3493)       !!                  ***  routine mppisl_real  ***
( 3494)       !!                  
( 3495)       !! ** Purpose :   Massively parallel processors
( 3496)       !!       Find the  non zero island barotropic stream function value
( 3497)       !!
( 3498)       !!     Modifications:
( 3499)       !!        !  93-09 (M. Imbard)
( 3500)       !!        !  96-05 (j. Escobar)
( 3501)       !!        !  98-05 (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI 
( 3502)       !!----------------------------------------------------------------------
( 3503)       REAL(wp), INTENT(inout) ::   ptab
( 3504) 
( 3505) 
( 3506) # 3534
( 3534)       !! * Local variables   (MPI version)
( 3535)       LOGICAL  ::   lcommute = .TRUE.
( 3536)       INTEGER  ::   mpi_isl, ierror
( 3537)       REAL(wp) ::   zwork
( 3538) 
( 3539)       CALL mpi_op_create( lc_isl, lcommute, mpi_isl, ierror )
( 3540)       CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision,   &
( 3541)          &                                mpi_isl  , mpi_comm_opa, ierror )
( 3542)       ptab = zwork
( 3543) 
( 3544) 
( 3545) # 3546
( 3546)    END SUBROUTINE mppisl_real
( 3547) 
( 3548) 
( 3549)   FUNCTION lc_isl( py, px, kdim, kdtatyp )
( 3550)     INTEGER :: kdim
( 3551)     REAL(wp), DIMENSION(kdim) ::  px, py
( 3552)     INTEGER :: kdtatyp, ji
( 3553)     INTEGER :: lc_isl
( 3554)     DO ji = 1, kdim
( 3555)        IF( py(ji) /= 0. )   px(ji) = py(ji)
( 3556)     END DO
( 3557)     lc_isl=0
( 3558) 
( 3559)   END FUNCTION lc_isl
( 3560) 
( 3561) 
( 3562)   SUBROUTINE mppmax_a_real( ptab, kdim )
( 3563)     !!----------------------------------------------------------------------
( 3564)     !!                 ***  routine mppmax_a_real  ***
( 3565)     !!                  
( 3566)     !! ** Purpose :   Maximum
( 3567)     !!
( 3568)     !!----------------------------------------------------------------------
( 3569)     !! * Arguments
( 3570)     INTEGER , INTENT( in  )                  ::   kdim
( 3571)     REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab
( 3572) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 153

( 3573) 
( 3574) # 3602
( 3602)     !! * Local variables   (MPI version)
( 3603)     INTEGER :: ierror
( 3604)     REAL(wp), DIMENSION(kdim) ::  zwork
( 3605) 
( 3606)     CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
( 3607)          ,mpi_max,mpi_comm_opa,ierror)
( 3608)     ptab(:) = zwork(:)
( 3609) 
( 3610) 
( 3611) # 3612
( 3612)   END SUBROUTINE mppmax_a_real
( 3613) 
( 3614) 
( 3615)   SUBROUTINE mppmax_real( ptab )
( 3616)     !!----------------------------------------------------------------------
( 3617)     !!                  ***  routine mppmax_real  ***
( 3618)     !!                    
( 3619)     !! ** Purpose :   Maximum
( 3620)     !!
( 3621)     !!----------------------------------------------------------------------
( 3622)     !! * Arguments
( 3623)     REAL(wp), INTENT(inout) ::   ptab      ! ???
( 3624) 
( 3625) 
( 3626) # 3646
( 3646)     !! * Local variables   (MPI version)
( 3647)     INTEGER  ::   ierror
( 3648)     REAL(wp) ::   zwork
( 3649) 
( 3650)     CALL mpi_allreduce( ptab, zwork  , 1             , mpi_double_precision,   &
( 3651)        &                      mpi_max, mpi_comm_opa, ierror     )
( 3652)     ptab = zwork
( 3653) 
( 3654) 
( 3655) # 3656
( 3656)   END SUBROUTINE mppmax_real
( 3657) 
( 3658) 
( 3659)   SUBROUTINE mppmin_a_real( ptab, kdim )
( 3660)     !!----------------------------------------------------------------------
( 3661)     !!                 ***  routine mppmin_a_real  ***
( 3662)     !!                  
( 3663)     !! ** Purpose :   Minimum
( 3664)     !!
( 3665)     !!-----------------------------------------------------------------------
( 3666)     !! * Arguments
( 3667)     INTEGER , INTENT( in  )                  ::   kdim
( 3668)     REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab
( 3669) 
( 3670) 
( 3671) # 3699
( 3699)     !! * Local variables   (MPI version)
( 3700)     INTEGER :: ierror
( 3701)     REAL(wp), DIMENSION(kdim) ::   zwork
( 3702) 
( 3703)     CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 154

( 3704)          ,mpi_min,mpi_comm_opa,ierror)
( 3705)     ptab(:) = zwork(:)
( 3706) 
( 3707) 
( 3708) # 3709
( 3709)   END SUBROUTINE mppmin_a_real
( 3710) 
( 3711) 
( 3712)   SUBROUTINE mppmin_real( ptab )
( 3713)     !!----------------------------------------------------------------------
( 3714)     !!                  ***  routine mppmin_real  ***
( 3715)     !! 
( 3716)     !! ** Purpose :   minimum in Massively Parallel Processing
( 3717)     !!                REAL scalar case
( 3718)     !!
( 3719)     !!-----------------------------------------------------------------------
( 3720)     !! * Arguments
( 3721)     REAL(wp), INTENT( inout ) ::   ptab        ! 
( 3722) 
( 3723) 
( 3724) # 3744
( 3744)     !! * Local variables   (MPI version)
( 3745)     INTEGER  ::   ierror
( 3746)     REAL(wp) ::   zwork
( 3747) 
( 3748)     CALL mpi_allreduce( ptab, zwork, 1,mpi_double_precision   &
( 3749)          &               ,mpi_min,mpi_comm_opa,ierror)
( 3750)     ptab = zwork
( 3751) 
( 3752) 
( 3753) # 3754
( 3754)   END SUBROUTINE mppmin_real
( 3755) 
( 3756) 
( 3757)   SUBROUTINE mppsum_a_real( ptab, kdim )
( 3758)     !!----------------------------------------------------------------------
( 3759)     !!                  ***  routine mppsum_a_real  ***
( 3760)     !! 
( 3761)     !! ** Purpose :   global sum in Massively Parallel Processing
( 3762)     !!                REAL ARRAY argument case
( 3763)     !!
( 3764)     !!-----------------------------------------------------------------------
( 3765)     INTEGER , INTENT( in )                     ::   kdim      ! size of ptab
( 3766)     REAL(wp), DIMENSION(kdim), INTENT( inout ) ::   ptab      ! input array
( 3767) 
( 3768) 
( 3769) # 3797
( 3797)     !! * Local variables   (MPI version)
( 3798)     INTEGER                   ::   ierror    ! temporary integer
( 3799)     REAL(wp), DIMENSION(kdim) ::   zwork     ! temporary workspace 
( 3800) 
( 3801)     CALL mpi_allreduce(ptab, zwork,kdim,mpi_double_precision   &
( 3802)          &              ,mpi_sum,mpi_comm_opa,ierror)
( 3803)     ptab(:) = zwork(:)
( 3804) 
( 3805) 
( 3806) # 3807
( 3807)   END SUBROUTINE mppsum_a_real






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 155

( 3808) 
( 3809) 
( 3810)   SUBROUTINE mppsum_real( ptab )
( 3811)     !!----------------------------------------------------------------------
( 3812)     !!                  ***  routine mppsum_real  ***
( 3813)     !!              
( 3814)     !! ** Purpose :   global sum in Massively Parallel Processing
( 3815)     !!                SCALAR argument case
( 3816)     !!
( 3817)     !!-----------------------------------------------------------------------
( 3818)     REAL(wp), INTENT(inout) ::   ptab        ! input scalar
( 3819) 
( 3820) 
( 3821) # 3841
( 3841)     !! * Local variables   (MPI version)
( 3842)     INTEGER  ::   ierror
( 3843)     REAL(wp) ::   zwork
( 3844) 
( 3845)     CALL mpi_allreduce(ptab, zwork, 1,mpi_double_precision   &
( 3846)          &              ,mpi_sum,mpi_comm_opa,ierror)
( 3847)     ptab = zwork
( 3848) 
( 3849) 
( 3850) # 3851
( 3851)   END SUBROUTINE mppsum_real
( 3852) 
( 3853)   SUBROUTINE mppsum_realdd( ptab )
( 3854)     !!----------------------------------------------------------------------
( 3855)     !!                  ***  routine mppsum_realdd ***
( 3856)     !!
( 3857)     !! ** Purpose :   global sum in Massively Parallel Processing
( 3858)     !!                SCALAR argument case for double-double precision
( 3859)     !!
( 3860)     !!-----------------------------------------------------------------------
( 3861)     COMPLEX(wp), INTENT(inout)         :: ptab    ! input scalar
( 3862) 
( 3863) 
( 3864) # 3869
( 3869)     !! * Local variables   (MPI version)
( 3870)     INTEGER  ::    ierror
( 3871)     COMPLEX(wp) :: zwork
( 3872) 
( 3873) ! reduce local sums into global sum
( 3874)     call MPI_ALLREDUCE (ptab, zwork, 1, MPI_DOUBLE_COMPLEX, &
( 3875)                      MPI_SUMDD, MPI_COMM_WORLD, ierror)
( 3876)     ptab = zwork
( 3877) 
( 3878) # 3879
( 3879)   END SUBROUTINE mppsum_realdd
( 3880) 
( 3881)   SUBROUTINE mppsum_a_realdd( ptab, kdim )
( 3882)     !!----------------------------------------------------------------------
( 3883)     !!                  ***  routine mppsum_a_realdd  ***
( 3884)     !!
( 3885)     !! ** Purpose :   global sum in Massively Parallel Processing
( 3886)     !!                COMPLEX ARRAY case for double-double precision
( 3887)     !!
( 3888)     !!-----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 156

( 3889)     INTEGER , INTENT( in )                        ::   kdim      ! size of ptab
( 3890)     COMPLEX(wp), DIMENSION(kdim), INTENT( inout ) ::   ptab      ! input array
( 3891) 
( 3892) 
( 3893) # 3898
( 3898)     !! * Local variables   (MPI version)
( 3899)     INTEGER                      :: ierror    ! temporary integer
( 3900)     COMPLEX(wp), DIMENSION(kdim) :: zwork     ! temporary workspace
( 3901) 
( 3902) 
( 3903)     call MPI_ALLREDUCE (ptab, zwork, kdim, MPI_DOUBLE_COMPLEX, &
( 3904)                      MPI_SUMDD, MPI_COMM_WORLD, ierror)
( 3905)     ptab(:) = zwork(:)
( 3906) 
( 3907) 
( 3908) # 3909
( 3909)   END SUBROUTINE mppsum_a_realdd
( 3910) 
( 3911) 
( 3912) 
( 3913)   SUBROUTINE mpp_minloc2d(ptab, pmask, pmin, ki,kj )
( 3914)     !!------------------------------------------------------------------------
( 3915)     !!             ***  routine mpp_minloc  ***
( 3916)     !!
( 3917)     !! ** Purpose :  Compute the global minimum of an array ptab
( 3918)     !!              and also give its global position
( 3919)     !!
( 3920)     !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
( 3921)     !!
( 3922)     !! ** Arguments : I : ptab =local 2D array
( 3923)     !!                O : pmin = global minimum
( 3924)     !!                O : ki,kj = global position of minimum
( 3925)     !!
( 3926)     !! ** Author : J.M. Molines 10/10/2004
( 3927)     !!--------------------------------------------------------------------------
( 3928) # 3931
( 3931)     !! * Arguments
( 3932)     REAL(wp), DIMENSION (jpi,jpj), INTENT (in)  :: ptab ,& ! Local 2D array
( 3933)          &                                         pmask   ! Local mask
( 3934)     REAL(wp)                     , INTENT (out) :: pmin    ! Global minimum of ptab
( 3935)     INTEGER                      , INTENT (out) :: ki,kj   ! index of minimum in global frame
( 3936) 
( 3937)     !! * Local variables
( 3938)     REAL(wp) :: zmin   ! local minimum
( 3939)     REAL(wp) ,DIMENSION(2,1) :: zain, zaout
( 3940)     INTEGER, DIMENSION (2)  :: ilocs
( 3941)     INTEGER :: ierror
( 3942) 
( 3943) 
( 3944)     zmin  = MINVAL( ptab(:,:) , mask= pmask == 1.e0 )
( 3945)     ilocs = MINLOC( ptab(:,:) , mask= pmask == 1.e0 )
( 3946) 
( 3947)     ki = ilocs(1) + nimpp - 1
( 3948)     kj = ilocs(2) + njmpp - 1
( 3949) 
( 3950)     zain(1,:)=zmin
( 3951)     zain(2,:)=ki+10000.*kj
( 3952) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 157

( 3953)     CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_OPA,ierror)
( 3954) 
( 3955)     pmin=zaout(1,1)
( 3956)     kj= INT(zaout(2,1)/10000.)
( 3957)     ki= INT(zaout(2,1) - 10000.*kj )
( 3958) 
( 3959) # 3960
( 3960)   END SUBROUTINE mpp_minloc2d
( 3961) 
( 3962) 
( 3963)   SUBROUTINE mpp_minloc3d(ptab, pmask, pmin, ki,kj ,kk)
( 3964)     !!------------------------------------------------------------------------
( 3965)     !!             ***  routine mpp_minloc  ***
( 3966)     !!
( 3967)     !! ** Purpose :  Compute the global minimum of an array ptab
( 3968)     !!              and also give its global position
( 3969)     !!
( 3970)     !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
( 3971)     !!
( 3972)     !! ** Arguments : I : ptab =local 2D array
( 3973)     !!                O : pmin = global minimum
( 3974)     !!                O : ki,kj = global position of minimum
( 3975)     !!
( 3976)     !! ** Author : J.M. Molines 10/10/2004
( 3977)     !!--------------------------------------------------------------------------
( 3978) # 3981
( 3981)     !! * Arguments
( 3982)     REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT (in)  :: ptab ,& ! Local 2D array
( 3983)          &                                         pmask   ! Local mask
( 3984)     REAL(wp)                     , INTENT (out) :: pmin    ! Global minimum of ptab
( 3985)     INTEGER                      , INTENT (out) :: ki,kj,kk ! index of minimum in global frame
( 3986) 
( 3987)     !! * Local variables
( 3988)     REAL(wp) :: zmin   ! local minimum
( 3989)     REAL(wp) ,DIMENSION(2,1) :: zain, zaout
( 3990)     INTEGER, DIMENSION (3)  :: ilocs
( 3991)     INTEGER :: ierror
( 3992) 
( 3993) 
( 3994)     zmin  = MINVAL( ptab(:,:,:) , mask= pmask == 1.e0 )
( 3995)     ilocs = MINLOC( ptab(:,:,:) , mask= pmask == 1.e0 )
( 3996) 
( 3997)     ki = ilocs(1) + nimpp - 1
( 3998)     kj = ilocs(2) + njmpp - 1
( 3999)     kk = ilocs(3)
( 4000) 
( 4001)     zain(1,:)=zmin
( 4002)     zain(2,:)=ki+10000.*kj+100000000.*kk
( 4003) 
( 4004)     CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_OPA,ierror)
( 4005) 
( 4006)     pmin=zaout(1,1)
( 4007)     kk= INT(zaout(2,1)/100000000.)
( 4008)     kj= INT(zaout(2,1) - kk * 100000000. )/10000
( 4009)     ki= INT(zaout(2,1) - kk * 100000000. -kj * 10000. )
( 4010) 
( 4011) # 4012
( 4012)   END SUBROUTINE mpp_minloc3d






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 158

( 4013) 
( 4014) 
( 4015)   SUBROUTINE mpp_maxloc2d(ptab, pmask, pmax, ki,kj )
( 4016)     !!------------------------------------------------------------------------
( 4017)     !!             ***  routine mpp_maxloc  ***
( 4018)     !!
( 4019)     !! ** Purpose :  Compute the global maximum of an array ptab
( 4020)     !!              and also give its global position
( 4021)     !!
( 4022)     !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
( 4023)     !!
( 4024)     !! ** Arguments : I : ptab =local 2D array
( 4025)     !!                O : pmax = global maximum
( 4026)     !!                O : ki,kj = global position of maximum
( 4027)     !!
( 4028)     !! ** Author : J.M. Molines 10/10/2004
( 4029)     !!--------------------------------------------------------------------------
( 4030) # 4033
( 4033)     !! * Arguments
( 4034)     REAL(wp), DIMENSION (jpi,jpj), INTENT (in)  :: ptab ,& ! Local 2D array
( 4035)          &                                         pmask   ! Local mask
( 4036)     REAL(wp)                     , INTENT (out) :: pmax    ! Global maximum of ptab
( 4037)     INTEGER                      , INTENT (out) :: ki,kj   ! index of maximum in global frame
( 4038) 
( 4039)     !! * Local variables
( 4040)     REAL(wp) :: zmax   ! local maximum
( 4041)     REAL(wp) ,DIMENSION(2,1) :: zain, zaout
( 4042)     INTEGER, DIMENSION (2)  :: ilocs
( 4043)     INTEGER :: ierror
( 4044) 
( 4045) 
( 4046)     zmax  = MAXVAL( ptab(:,:) , mask= pmask == 1.e0 )
( 4047)     ilocs = MAXLOC( ptab(:,:) , mask= pmask == 1.e0 )
( 4048) 
( 4049)     ki = ilocs(1) + nimpp - 1
( 4050)     kj = ilocs(2) + njmpp - 1
( 4051) 
( 4052)     zain(1,:)=zmax
( 4053)     zain(2,:)=ki+10000.*kj
( 4054) 
( 4055)     CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_OPA,ierror)
( 4056) 
( 4057)     pmax=zaout(1,1)
( 4058)     kj= INT(zaout(2,1)/10000.)
( 4059)     ki= INT(zaout(2,1) - 10000.*kj )
( 4060) 
( 4061) # 4062
( 4062)   END SUBROUTINE mpp_maxloc2d
( 4063) 
( 4064)   SUBROUTINE mpp_maxloc3d(ptab, pmask, pmax, ki,kj,kk )
( 4065)     !!------------------------------------------------------------------------
( 4066)     !!             ***  routine mpp_maxloc  ***
( 4067)     !!
( 4068)     !! ** Purpose :  Compute the global maximum of an array ptab
( 4069)     !!              and also give its global position
( 4070)     !!
( 4071)     !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
( 4072)     !!






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 159

( 4073)     !! ** Arguments : I : ptab =local 2D array
( 4074)     !!                O : pmax = global maximum
( 4075)     !!                O : ki,kj = global position of maximum
( 4076)     !!
( 4077)     !! ** Author : J.M. Molines 10/10/2004
( 4078)     !!--------------------------------------------------------------------------
( 4079) # 4082
( 4082)     !! * Arguments
( 4083)     REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT (in)  :: ptab ,& ! Local 2D array
( 4084)          &                                         pmask   ! Local mask
( 4085)     REAL(wp)                     , INTENT (out) :: pmax    ! Global maximum of ptab
( 4086)     INTEGER                      , INTENT (out) :: ki,kj,kk   ! index of maximum in global frame
( 4087) 
( 4088)     !! * Local variables
( 4089)     REAL(wp) :: zmax   ! local maximum
( 4090)     REAL(wp) ,DIMENSION(2,1) :: zain, zaout
( 4091)     INTEGER, DIMENSION (3)  :: ilocs
( 4092)     INTEGER :: ierror
( 4093) 
( 4094) 
( 4095)     zmax  = MAXVAL( ptab(:,:,:) , mask= pmask == 1.e0 )
( 4096)     ilocs = MAXLOC( ptab(:,:,:) , mask= pmask == 1.e0 )
( 4097) 
( 4098)     ki = ilocs(1) + nimpp - 1
( 4099)     kj = ilocs(2) + njmpp - 1
( 4100)     kk = ilocs(3)
( 4101) 
( 4102)     zain(1,:)=zmax
( 4103)     zain(2,:)=ki+10000.*kj+100000000.*kk
( 4104) 
( 4105)     CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_OPA,ierror)
( 4106) 
( 4107)     pmax=zaout(1,1)
( 4108)     kk= INT(zaout(2,1)/100000000.)
( 4109)     kj= INT(zaout(2,1) - kk * 100000000. )/10000
( 4110)     ki= INT(zaout(2,1) - kk * 100000000. -kj * 10000. )
( 4111) 
( 4112) # 4113
( 4113)   END SUBROUTINE mpp_maxloc3d
( 4114) 
( 4115)   SUBROUTINE mppsync()
( 4116)     !!----------------------------------------------------------------------
( 4117)     !!                  ***  routine mppsync  ***
( 4118)     !!                   
( 4119)     !! ** Purpose :   Massively parallel processors, synchroneous
( 4120)     !!
( 4121)     !!-----------------------------------------------------------------------
( 4122) 
( 4123) 
( 4124) # 4130
( 4130)     !! * Local variables   (MPI version)
( 4131)     INTEGER :: ierror
( 4132) 
( 4133)     CALL mpi_barrier(mpi_comm_opa,ierror)
( 4134) 
( 4135) 
( 4136) # 4137
( 4137)   END SUBROUTINE mppsync






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 160

( 4138) 
( 4139) 
( 4140)   SUBROUTINE mppstop
( 4141)     !!----------------------------------------------------------------------
( 4142)     !!                  ***  routine mppstop  ***
( 4143)     !!                   
( 4144)     !! ** purpose :   Stop massilively parallel processors method
( 4145)     !!
( 4146)     !!----------------------------------------------------------------------
( 4147)     !! * Local declarations
( 4148)     INTEGER ::   info
( 4149)     !!----------------------------------------------------------------------
( 4150) 
( 4151)     ! 1. Mpp synchroneus
( 4152)     ! ------------------
( 4153) 
( 4154)     CALL mppsync
( 4155) # 4156
( 4156)     CALL mpi_finalize( info )
( 4157) 
( 4158) # 4159
( 4159)   END SUBROUTINE mppstop
( 4160) 
( 4161) 
( 4162)   SUBROUTINE mppobc( ptab, kd1, kd2, kl, kk, ktype, kij )
( 4163)     !!----------------------------------------------------------------------
( 4164)     !!                  ***  routine mppobc  ***
( 4165)     !! 
( 4166)     !! ** Purpose :   Message passing manadgement for open boundary
( 4167)     !!     conditions array
( 4168)     !!
( 4169)     !! ** Method  :   Use mppsend and mpprecv function for passing mask
( 4170)     !!       between processors following neighboring subdomains.
( 4171)     !!       domain parameters
( 4172)     !!                    nlci   : first dimension of the local subdomain
( 4173)     !!                    nlcj   : second dimension of the local subdomain
( 4174)     !!                    nbondi : mark for "east-west local boundary"
( 4175)     !!                    nbondj : mark for "north-south local boundary"
( 4176)     !!                    noea   : number for local neighboring processors 
( 4177)     !!                    nowe   : number for local neighboring processors
( 4178)     !!                    noso   : number for local neighboring processors
( 4179)     !!                    nono   : number for local neighboring processors
( 4180)     !!
( 4181)     !! History :
( 4182)     !!        !  98-07 (J.M. Molines) Open boundary conditions
( 4183)     !!----------------------------------------------------------------------
( 4184)     !! * Arguments
( 4185)     INTEGER , INTENT( in ) ::   &
( 4186)          kd1, kd2,   &  ! starting and ending indices
( 4187)          kl ,        &  ! index of open boundary
( 4188)          kk,         &  ! vertical dimension
( 4189)          ktype,      &  ! define north/south or east/west cdt
( 4190)          !              !  = 1  north/south  ;  = 2  east/west
( 4191)          kij            ! horizontal dimension
( 4192)     REAL(wp), DIMENSION(kij,kk), INTENT( inout )  ::   &
( 4193)          ptab           ! variable array
( 4194) 
( 4195)     !! * Local variables






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 161

( 4196)     INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
( 4197)     INTEGER  ::   &
( 4198)          iipt0, iipt1, ilpt1,     &  ! temporary integers
( 4199)          ijpt0, ijpt1,            &  !    "          "
( 4200)          imigr, iihom, ijhom         !    "          "
( 4201)     INTEGER ::   ml_req1, ml_req2, ml_err     ! for key_mpi_isend
( 4202)     INTEGER ::   ml_stat(MPI_STATUS_SIZE)     ! for key_mpi_isend
( 4203)     REAL(wp), DIMENSION(jpi,jpj) ::   &
( 4204)          ztab                        ! temporary workspace
( 4205)     !!----------------------------------------------------------------------
( 4206) 
( 4207) 
( 4208)     ! boundary condition initialization
( 4209)     ! ---------------------------------
( 4210) 
( 4211)     ztab(:,:) = 0.e0
( 4212) 
( 4213)     IF( ktype==1 ) THEN                                  ! north/south boundaries
( 4214)        iipt0 = MAX( 1, MIN(kd1 - nimpp+1, nlci     ) )
( 4215)        iipt1 = MAX( 0, MIN(kd2 - nimpp+1, nlci - 1 ) )
( 4216)        ilpt1 = MAX( 1, MIN(kd2 - nimpp+1, nlci     ) )
( 4217)        ijpt0 = MAX( 1, MIN(kl  - njmpp+1, nlcj     ) )
( 4218)        ijpt1 = MAX( 0, MIN(kl  - njmpp+1, nlcj - 1 ) )
( 4219)     ELSEIF( ktype==2 ) THEN                              ! east/west boundaries
( 4220)        iipt0 = MAX( 1, MIN(kl  - nimpp+1, nlci     ) )
( 4221)        iipt1 = MAX( 0, MIN(kl  - nimpp+1, nlci - 1 ) )
( 4222)        ijpt0 = MAX( 1, MIN(kd1 - njmpp+1, nlcj     ) )
( 4223)        ijpt1 = MAX( 0, MIN(kd2 - njmpp+1, nlcj - 1 ) )
( 4224)        ilpt1 = MAX( 1, MIN(kd2 - njmpp+1, nlcj     ) )
( 4225)     ELSE
( 4226)        CALL ctl_stop( 'mppobc: bad ktype' )
( 4227)     ENDIF
( 4228) 
( 4229)     DO jk = 1, kk
( 4230)        IF( ktype==1 ) THEN                               ! north/south boundaries
( 4231)           DO jj = ijpt0, ijpt1
( 4232)              DO ji = iipt0, iipt1
( 4233)                 ztab(ji,jj) = ptab(ji,jk)
( 4234)              END DO
( 4235)           END DO
( 4236)        ELSEIF( ktype==2 ) THEN                           ! east/west boundaries
( 4237)           DO jj = ijpt0, ijpt1
( 4238)              DO ji = iipt0, iipt1
( 4239)                 ztab(ji,jj) = ptab(jj,jk)
( 4240)              END DO
( 4241)           END DO
( 4242)        ENDIF
( 4243) 
( 4244) 
( 4245)        ! 1. East and west directions
( 4246)        ! ---------------------------
( 4247) 
( 4248)        ! 1.1 Read Dirichlet lateral conditions
( 4249) 
( 4250)        IF( nbondi /= 2 ) THEN
( 4251)           iihom = nlci-nreci
( 4252) 
( 4253)           DO jl = 1, jpreci






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 162

( 4254)              t2ew(:,jl,1) = ztab(jpreci+jl,:)
( 4255)              t2we(:,jl,1) = ztab(iihom +jl,:)
( 4256)           END DO
( 4257)        ENDIF
( 4258) 
( 4259)        ! 1.2 Migrations
( 4260) 
( 4261) # 4277
( 4277)        !! * (MPI version)
( 4278) 
( 4279)        imigr=jpreci*jpj
( 4280) 
( 4281)        IF( nbondi == -1 ) THEN
( 4282)           CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req1)
( 4283)           CALL mpprecv(1,t2ew(1,1,2),imigr)
( 4284)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4285)        ELSEIF( nbondi == 0 ) THEN
( 4286)           CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
( 4287)           CALL mppsend(2,t2we(1,1,1),imigr,noea, ml_req2)
( 4288)           CALL mpprecv(1,t2ew(1,1,2),imigr)
( 4289)           CALL mpprecv(2,t2we(1,1,2),imigr)
( 4290)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4291)           IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 4292)        ELSEIF( nbondi == 1 ) THEN
( 4293)           CALL mppsend(1,t2ew(1,1,1),imigr,nowe, ml_req1)
( 4294)           CALL mpprecv(2,t2we(1,1,2),imigr)
( 4295)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4296)        ENDIF
( 4297) 
( 4298) 
( 4299) # 4300
( 4300)        ! 1.3 Write Dirichlet lateral conditions
( 4301) 
( 4302)        iihom = nlci-jpreci
( 4303)        IF( nbondi == 0 .OR. nbondi == 1 ) THEN
( 4304)           DO jl = 1, jpreci
( 4305)              ztab(jl,:) = t2we(:,jl,2)
( 4306)           END DO
( 4307)        ENDIF
( 4308) 
( 4309)        IF( nbondi == -1 .OR. nbondi == 0 ) THEN
( 4310)           DO jl = 1, jpreci
( 4311)              ztab(iihom+jl,:) = t2ew(:,jl,2)
( 4312)           END DO
( 4313)        ENDIF
( 4314) 
( 4315) 
( 4316)        ! 2. North and south directions
( 4317)        ! -----------------------------
( 4318) 
( 4319)        ! 2.1 Read Dirichlet lateral conditions
( 4320) 
( 4321)        IF( nbondj /= 2 ) THEN
( 4322)           ijhom = nlcj-nrecj
( 4323)           DO jl = 1, jprecj
( 4324)              t2sn(:,jl,1) = ztab(:,ijhom +jl)
( 4325)              t2ns(:,jl,1) = ztab(:,jprecj+jl)
( 4326)           END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 163

( 4327)        ENDIF
( 4328) 
( 4329)        ! 2.2 Migrations
( 4330) 
( 4331) # 4348
( 4348)        !! * Local variables   (MPI version)
( 4349) 
( 4350)        imigr=jprecj*jpi
( 4351) 
( 4352)        IF( nbondj == -1 ) THEN
( 4353)           CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req1)
( 4354)           CALL mpprecv(3,t2ns(1,1,2),imigr)
( 4355)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4356)        ELSEIF( nbondj == 0 ) THEN
( 4357)           CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
( 4358)           CALL mppsend(4,t2sn(1,1,1),imigr,nono, ml_req2)
( 4359)           CALL mpprecv(3,t2ns(1,1,2),imigr)
( 4360)           CALL mpprecv(4,t2sn(1,1,2),imigr)
( 4361)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4362)           IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
( 4363)        ELSEIF( nbondj == 1 ) THEN
( 4364)           CALL mppsend(3,t2ns(1,1,1),imigr,noso, ml_req1)
( 4365)           CALL mpprecv(4,t2sn(1,1,2),imigr)
( 4366)           IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
( 4367)        ENDIF
( 4368) 
( 4369) 
( 4370) # 4371
( 4371)        ! 2.3 Write Dirichlet lateral conditions
( 4372) 
( 4373)        ijhom = nlcj - jprecj
( 4374)        IF( nbondj == 0 .OR. nbondj == 1 ) THEN
( 4375)           DO jl = 1, jprecj
( 4376)              ztab(:,jl) = t2sn(:,jl,2)
( 4377)           END DO
( 4378)        ENDIF
( 4379) 
( 4380)        IF( nbondj == 0 .OR. nbondj == -1 ) THEN
( 4381)           DO jl = 1, jprecj
( 4382)              ztab(:,ijhom+jl) = t2ns(:,jl,2)
( 4383)           END DO
( 4384)        ENDIF
( 4385) 
( 4386)        IF( ktype==1 .AND. kd1 <= jpi+nimpp-1 .AND. nimpp <= kd2 ) THEN
( 4387)           ! north/south boundaries
( 4388)           DO jj = ijpt0,ijpt1
( 4389)              DO ji = iipt0,ilpt1
( 4390)                 ptab(ji,jk) = ztab(ji,jj)  
( 4391)              END DO
( 4392)           END DO
( 4393)        ELSEIF( ktype==2 .AND. kd1 <= jpj+njmpp-1 .AND. njmpp <= kd2 ) THEN
( 4394)           ! east/west boundaries
( 4395)           DO jj = ijpt0,ilpt1
( 4396)              DO ji = iipt0,iipt1
( 4397)                 ptab(jj,jk) = ztab(ji,jj) 
( 4398)              END DO
( 4399)           END DO
( 4400)        ENDIF






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 164

( 4401) 
( 4402)     END DO
( 4403) 
( 4404)   END SUBROUTINE mppobc
( 4405) 
( 4406) 
( 4407)   SUBROUTINE mpp_ini_north
( 4408)     !!----------------------------------------------------------------------
( 4409)     !!               ***  routine mpp_ini_north  ***
( 4410)     !!
( 4411)     !! ** Purpose :   Initialize special communicator for north folding 
( 4412)     !!      condition together with global variables needed in the mpp folding
( 4413)     !!
( 4414)     !! ** Method  : - Look for northern processors
( 4415)     !!              - Put their number in nrank_north
( 4416)     !!              - Create groups for the world processors and the north processors
( 4417)     !!              - Create a communicator for northern processors
( 4418)     !!
( 4419)     !! ** output
( 4420)     !!      njmppmax = njmpp for northern procs
( 4421)     !!      ndim_rank_north = number of processors in the northern line
( 4422)     !!      nrank_north (ndim_rank_north) = number  of the northern procs.
( 4423)     !!      ngrp_world = group ID for the world processors
( 4424)     !!      ngrp_north = group ID for the northern processors
( 4425)     !!      ncomm_north = communicator for the northern procs.
( 4426)     !!      north_root = number (in the world) of proc 0 in the northern comm.
( 4427)     !!
( 4428)     !! History :
( 4429)     !!        !  03-09 (J.M. Molines, MPI only )
( 4430)     !!----------------------------------------------------------------------
( 4431) # 4434
( 4434)     INTEGER :: ierr
( 4435)     INTEGER :: jproc
( 4436)     INTEGER :: ii,ji
( 4437)     !!----------------------------------------------------------------------
( 4438) 
( 4439)     njmppmax=MAXVAL(njmppt)
( 4440) 
( 4441)     ! Look for how many procs on the northern boundary
( 4442)     !
( 4443)     ndim_rank_north=0
( 4444)     DO jproc=1,jpnij
( 4445)        IF ( njmppt(jproc) == njmppmax ) THEN
( 4446)           ndim_rank_north = ndim_rank_north + 1
( 4447)        END IF
( 4448)     END DO
( 4449) 
( 4450) 
( 4451)     ! Allocate the right size to nrank_north
( 4452)     !
( 4453)     ALLOCATE(nrank_north(ndim_rank_north))
( 4454) 
( 4455)     ! Fill the nrank_north array with proc. number of northern procs.
( 4456)     ! Note : the rank start at 0 in MPI
( 4457)     !
( 4458)     ii=0
( 4459)     DO ji = 1, jpnij
( 4460)        IF ( njmppt(ji) == njmppmax   ) THEN






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 165

( 4461)           ii=ii+1
( 4462)           nrank_north(ii)=ji-1
( 4463)        END IF
( 4464)     END DO
( 4465)     ! create the world group
( 4466)     !
( 4467)     CALL MPI_COMM_GROUP(mpi_comm_opa,ngrp_world,ierr)
( 4468)     !
( 4469)     ! Create the North group from the world group
( 4470)     CALL MPI_GROUP_INCL(ngrp_world,ndim_rank_north,nrank_north,ngrp_north,ierr)
( 4471) 
( 4472)     ! Create the North communicator , ie the pool of procs in the north group
( 4473)     !
( 4474)     CALL MPI_COMM_CREATE(mpi_comm_opa,ngrp_north,ncomm_north,ierr)
( 4475) 
( 4476) 
( 4477)     ! find proc number in the world of proc 0 in the north
( 4478)     CALL MPI_GROUP_TRANSLATE_RANKS(ngrp_north,1,0,ngrp_world,north_root,ierr)
( 4479) 
( 4480) # 4481
( 4481)   END SUBROUTINE mpp_ini_north
( 4482) 
( 4483) 
( 4484)    SUBROUTINE mpp_lbc_north_3d ( pt3d, cd_type, psgn )
( 4485)       !!---------------------------------------------------------------------
( 4486)       !!                   ***  routine mpp_lbc_north_3d  ***
( 4487)       !!
( 4488)       !! ** Purpose :
( 4489)       !!      Ensure proper north fold horizontal bondary condition in mpp configuration
( 4490)       !!      in case of jpn1 > 1
( 4491)       !!
( 4492)       !! ** Method :
( 4493)       !!      Gather the 4 northern lines of the global domain on 1 processor and 
( 4494)       !!      apply lbc north-fold on this sub array. Then scatter the fold array 
( 4495)       !!      back to the processors.
( 4496)       !!
( 4497)       !! History :
( 4498)       !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
( 4499)       !!                                  from lbc routine
( 4500)       !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
( 4501)       !!----------------------------------------------------------------------
( 4502)       !! * Arguments
( 4503)       CHARACTER(len=1), INTENT( in ) ::   &
( 4504)          cd_type       ! nature of pt3d grid-points
( 4505)          !             !   = T ,  U , V , F or W  gridpoints
( 4506)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
( 4507)          pt3d          ! 3D array on which the boundary condition is applied
( 4508)       REAL(wp), INTENT( in ) ::   &
( 4509)          psgn          ! control of the sign change
( 4510)          !             !   = -1. , the sign is changed if north fold boundary
( 4511)          !             !   =  1. , the sign is kept  if north fold boundary
( 4512) 
( 4513)       !! * Local declarations
( 4514)       INTEGER :: ji, jj, jk, jr, jproc
( 4515)       INTEGER :: ierr
( 4516)       INTEGER :: ildi,ilei,iilb
( 4517)       INTEGER :: ijpj,ijpjm1,ij,ijt,iju
( 4518)       INTEGER :: itaille






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 166

( 4519)       REAL(wp), DIMENSION(jpiglo,4,jpk) :: ztab
( 4520)       REAL(wp), DIMENSION(jpi,4,jpk,jpni) :: znorthgloio
( 4521)       REAL(wp), DIMENSION(jpi,4,jpk) :: znorthloc
( 4522)       !!----------------------------------------------------------------------
( 4523) 
( 4524)     ! If we get in this routine it s because : North fold condition and mpp with more
( 4525)     !   than one proc across i : we deal only with the North condition
( 4526) 
( 4527)     ! 0. Sign setting
( 4528)     ! ---------------
( 4529) 
( 4530)     ijpj=4
( 4531)     ijpjm1=3
( 4532) 
( 4533)     ! put in znorthloc the last 4 jlines of pt3d
( 4534)     DO jk = 1, jpk 
( 4535)        DO jj = nlcj - ijpj +1, nlcj
( 4536)           ij = jj - nlcj + ijpj
( 4537)           znorthloc(:,ij,jk) = pt3d(:,jj,jk)
( 4538)        END DO
( 4539)     END DO
( 4540) 
( 4541) 
( 4542)     IF (npolj /= 0 ) THEN
( 4543)        ! Build in proc 0 of ncomm_north the znorthgloio
( 4544)        znorthgloio(:,:,:,:) = 0_wp
( 4545) 
( 4546) # 4549
( 4549)        itaille=jpi*jpk*ijpj
( 4550)        CALL MPI_GATHER(znorthloc,itaille,MPI_DOUBLE_PRECISION,znorthgloio,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 4551) 
( 4552) # 4553
( 4553)     ENDIF
( 4554) 
( 4555)     IF (narea == north_root+1 ) THEN
( 4556)        ! recover the global north array
( 4557)        ztab(:,:,:) = 0_wp
( 4558) 
( 4559)        DO jr = 1, ndim_rank_north
( 4560)           jproc = nrank_north(jr) + 1
( 4561)           ildi  = nldit (jproc)
( 4562)           ilei  = nleit (jproc)
( 4563)           iilb  = nimppt(jproc)
( 4564)           DO jk = 1, jpk 
( 4565)              DO jj = 1, 4
( 4566)                 DO ji = ildi, ilei
( 4567)                    ztab(ji+iilb-1,jj,jk) = znorthgloio(ji,jj,jk,jr)
( 4568)                 END DO
( 4569)              END DO
( 4570)           END DO
( 4571)        END DO
( 4572) 
( 4573) 
( 4574)        ! Horizontal slab
( 4575)        ! ===============
( 4576) 
( 4577)        DO jk = 1, jpk 
( 4578) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 167

( 4579) 
( 4580)           ! 2. North-Fold boundary conditions
( 4581)           ! ----------------------------------
( 4582) 
( 4583)           SELECT CASE ( npolj )
( 4584) 
( 4585)           CASE ( 3, 4 )                       ! *  North fold  T-point pivot
( 4586) 
( 4587)              ztab( 1    ,ijpj,jk) = 0.e0
( 4588)              ztab(jpiglo,ijpj,jk) = 0.e0
( 4589) 
( 4590)              SELECT CASE ( cd_type )
( 4591) 
( 4592)              CASE ( 'T' , 'S' , 'W' )                   ! T-, W-point
( 4593)                 DO ji = 2, jpiglo
( 4594)                    ijt = jpiglo-ji+2
( 4595)                    ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-2,jk)
( 4596)                 END DO
( 4597)                 DO ji = jpiglo/2+1, jpiglo
( 4598)                    ijt = jpiglo-ji+2
( 4599)                    ztab(ji,ijpjm1,jk) = psgn * ztab(ijt,ijpjm1,jk)
( 4600)                 END DO
( 4601) 
( 4602)              CASE ( 'U' )                               ! U-point
( 4603)                 DO ji = 1, jpiglo-1
( 4604)                    iju = jpiglo-ji+1
( 4605)                    ztab(ji,ijpj,jk) = psgn * ztab(iju,ijpj-2,jk)
( 4606)                 END DO
( 4607)                 DO ji = jpiglo/2, jpiglo-1
( 4608)                    iju = jpiglo-ji+1
( 4609)                    ztab(ji,ijpjm1,jk) = psgn * ztab(iju,ijpjm1,jk)
( 4610)                 END DO
( 4611) 
( 4612)              CASE ( 'V' )                               ! V-point
( 4613)                 DO ji = 2, jpiglo
( 4614)                    ijt = jpiglo-ji+2
( 4615)                    ztab(ji,ijpj-1,jk) = psgn * ztab(ijt,ijpj-2,jk)
( 4616)                    ztab(ji,ijpj  ,jk) = psgn * ztab(ijt,ijpj-3,jk)
( 4617)                 END DO
( 4618) 
( 4619)              CASE ( 'F' , 'G' )                         ! F-point
( 4620)                 DO ji = 1, jpiglo-1
( 4621)                    iju = jpiglo-ji+1
( 4622)                    ztab(ji,ijpj-1,jk) = psgn * ztab(iju,ijpj-2,jk)
( 4623)                    ztab(ji,ijpj  ,jk) = psgn * ztab(iju,ijpj-3,jk)
( 4624)                 END DO
( 4625) 
( 4626)              END SELECT
( 4627) 
( 4628)           CASE ( 5, 6 )                        ! *  North fold  F-point pivot
( 4629) 
( 4630)              ztab( 1    ,ijpj,jk) = 0.e0
( 4631)              ztab(jpiglo,ijpj,jk) = 0.e0
( 4632) 
( 4633)              SELECT CASE ( cd_type )
( 4634) 
( 4635)              CASE ( 'T' , 'S' , 'W' )                   ! T-, W-point
( 4636)                 DO ji = 1, jpiglo






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 168

( 4637)                    ijt = jpiglo-ji+1
( 4638)                    ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-1,jk)
( 4639)                 END DO
( 4640) 
( 4641)              CASE ( 'U' )                               ! U-point
( 4642)                 DO ji = 1, jpiglo-1
( 4643)                    iju = jpiglo-ji
( 4644)                    ztab(ji,ijpj,jk) = psgn * ztab(iju,ijpj-1,jk)
( 4645)                 END DO
( 4646) 
( 4647)              CASE ( 'V' )                               ! V-point
( 4648)                 DO ji = 1, jpiglo
( 4649)                    ijt = jpiglo-ji+1
( 4650)                    ztab(ji,ijpj,jk) = psgn * ztab(ijt,ijpj-2,jk)
( 4651)                 END DO
( 4652)                 DO ji = jpiglo/2+1, jpiglo
( 4653)                    ijt = jpiglo-ji+1
( 4654)                    ztab(ji,ijpjm1,jk) = psgn * ztab(ijt,ijpjm1,jk)
( 4655)                 END DO
( 4656) 
( 4657)              CASE ( 'F' , 'G' )                         ! F-point
( 4658)                 DO ji = 1, jpiglo-1
( 4659)                    iju = jpiglo-ji
( 4660)                    ztab(ji,ijpj  ,jk) = psgn * ztab(iju,ijpj-2,jk)
( 4661)                 END DO
( 4662)                 DO ji = jpiglo/2+1, jpiglo-1
( 4663)                    iju = jpiglo-ji
( 4664)                    ztab(ji,ijpjm1,jk) = psgn * ztab(iju,ijpjm1,jk)
( 4665)                 END DO
( 4666) 
( 4667)              END SELECT
( 4668) 
( 4669)           CASE DEFAULT                           ! *  closed
( 4670) 
( 4671)              SELECT CASE ( cd_type) 
( 4672) 
( 4673)              CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
( 4674)                 ztab(:, 1  ,jk) = 0.e0
( 4675)                 ztab(:,ijpj,jk) = 0.e0
( 4676) 
( 4677)              CASE ( 'F' )                               ! F-point
( 4678)                 ztab(:,ijpj,jk) = 0.e0
( 4679) 
( 4680)              END SELECT
( 4681) 
( 4682)           END SELECT
( 4683) 
( 4684)           !     End of slab
( 4685)           !     ===========
( 4686) 
( 4687)        END DO
( 4688) 
( 4689)        !! Scatter back to pt3d
( 4690)        DO jr = 1, ndim_rank_north
( 4691)           jproc=nrank_north(jr)+1
( 4692)           ildi=nldit (jproc)
( 4693)           ilei=nleit (jproc)
( 4694)           iilb=nimppt(jproc)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 169

( 4695)           DO jk=  1, jpk
( 4696)              DO jj=1,ijpj
( 4697)                 DO ji=ildi,ilei
( 4698)                    znorthgloio(ji,jj,jk,jr)=ztab(ji+iilb-1,jj,jk)
( 4699)                 END DO
( 4700)              END DO
( 4701)           END DO
( 4702)        END DO
( 4703) 
( 4704)     ENDIF      ! only done on proc 0 of ncomm_north
( 4705) 
( 4706) # 4709
( 4709)     IF ( npolj /= 0 ) THEN
( 4710)        itaille=jpi*jpk*ijpj
( 4711)        CALL MPI_SCATTER(znorthgloio,itaille,MPI_DOUBLE_PRECISION,znorthloc,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 4712)     ENDIF
( 4713) 
( 4714) # 4715
( 4715)     ! put in the last ijpj jlines of pt3d znorthloc
( 4716)     DO jk = 1 , jpk 
( 4717)        DO jj = nlcj - ijpj + 1 , nlcj
( 4718)           ij = jj - nlcj + ijpj
( 4719)           pt3d(:,jj,jk)= znorthloc(:,ij,jk)
( 4720)        END DO
( 4721)     END DO
( 4722) 
( 4723)   END SUBROUTINE mpp_lbc_north_3d
( 4724) 
( 4725) 
( 4726)   SUBROUTINE mpp_lbc_north_2d ( pt2d, cd_type, psgn)
( 4727)     !!---------------------------------------------------------------------
( 4728)     !!                   ***  routine mpp_lbc_north_2d  ***
( 4729)     !!
( 4730)     !! ** Purpose :
( 4731)     !!      Ensure proper north fold horizontal bondary condition in mpp configuration
( 4732)     !!      in case of jpn1 > 1 (for 2d array )
( 4733)     !!
( 4734)     !! ** Method :
( 4735)     !!      Gather the 4 northern lines of the global domain on 1 processor and 
( 4736)     !!      apply lbc north-fold on this sub array. Then scatter the fold array 
( 4737)     !!      back to the processors.
( 4738)     !!
( 4739)     !! History :
( 4740)     !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
( 4741)     !!                                  from lbc routine
( 4742)     !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
( 4743)     !!----------------------------------------------------------------------
( 4744) 
( 4745)     !! * Arguments
( 4746)     CHARACTER(len=1), INTENT( in ) ::   &
( 4747)          cd_type       ! nature of pt2d grid-points
( 4748)     !             !   = T ,  U , V , F or W  gridpoints
( 4749)     REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
( 4750)          pt2d          ! 2D array on which the boundary condition is applied
( 4751)     REAL(wp), INTENT( in ) ::   &
( 4752)          psgn          ! control of the sign change
( 4753)     !             !   = -1. , the sign is changed if north fold boundary
( 4754)     !             !   =  1. , the sign is kept  if north fold boundary






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 170

( 4755) 
( 4756) 
( 4757)     !! * Local declarations
( 4758) 
( 4759)     INTEGER :: ji, jj,  jr, jproc
( 4760)     INTEGER :: ierr
( 4761)     INTEGER :: ildi,ilei,iilb
( 4762)     INTEGER :: ijpj,ijpjm1,ij,ijt,iju
( 4763)     INTEGER :: itaille
( 4764) 
( 4765)     REAL(wp), DIMENSION(jpiglo,4) :: ztab
( 4766)     REAL(wp), DIMENSION(jpi,4,jpni) :: znorthgloio
( 4767)     REAL(wp), DIMENSION(jpi,4) :: znorthloc
( 4768)     !!----------------------------------------------------------------------
( 4769)     !!  OPA 8.5, LODYC-IPSL (2002)
( 4770)     !!----------------------------------------------------------------------
( 4771)     ! If we get in this routine it s because : North fold condition and mpp with more
( 4772)     !   than one proc across i : we deal only with the North condition
( 4773) 
( 4774)     ! 0. Sign setting
( 4775)     ! ---------------
( 4776) 
( 4777)     ijpj=4
( 4778)     ijpjm1=3
( 4779) 
( 4780) 
( 4781)     ! put in znorthloc the last 4 jlines of pt2d
( 4782)     DO jj = nlcj - ijpj +1, nlcj
( 4783)        ij = jj - nlcj + ijpj
( 4784)        znorthloc(:,ij)=pt2d(:,jj)
( 4785)     END DO
( 4786) 
( 4787)     IF (npolj /= 0 ) THEN
( 4788)        ! Build in proc 0 of ncomm_north the znorthgloio
( 4789)        znorthgloio(:,:,:) = 0_wp
( 4790) # 4793
( 4793)        itaille=jpi*ijpj
( 4794)        CALL MPI_GATHER(znorthloc,itaille,MPI_DOUBLE_PRECISION,znorthgloio,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 4795) # 4796
( 4796)     ENDIF
( 4797) 
( 4798)     IF (narea == north_root+1 ) THEN
( 4799)        ! recover the global north array
( 4800)        ztab(:,:) = 0_wp
( 4801) 
( 4802)        DO jr = 1, ndim_rank_north
( 4803)           jproc=nrank_north(jr)+1
( 4804)           ildi=nldit (jproc)
( 4805)           ilei=nleit (jproc)
( 4806)           iilb=nimppt(jproc)
( 4807)           DO jj=1,4
( 4808)              DO ji=ildi,ilei
( 4809)                 ztab(ji+iilb-1,jj)=znorthgloio(ji,jj,jr)
( 4810)              END DO
( 4811)           END DO
( 4812)        END DO
( 4813) 
( 4814) 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 171

( 4815)        ! 2. North-Fold boundary conditions
( 4816)        ! ----------------------------------
( 4817) 
( 4818)        SELECT CASE ( npolj )
( 4819) 
( 4820)        CASE ( 3, 4 )                       ! *  North fold  T-point pivot
( 4821) 
( 4822)           ztab( 1    ,ijpj) = 0.e0
( 4823)           ztab(jpiglo,ijpj) = 0.e0
( 4824) 
( 4825)           SELECT CASE ( cd_type )
( 4826) 
( 4827)           CASE ( 'T' , 'W' , 'S' )                         ! T-, W-point
( 4828)              DO ji = 2, jpiglo
( 4829)                 ijt = jpiglo-ji+2
( 4830)                 ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-2)
( 4831)              END DO
( 4832)              DO ji = jpiglo/2+1, jpiglo
( 4833)                 ijt = jpiglo-ji+2
( 4834)                 ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
( 4835)              END DO
( 4836) 
( 4837)           CASE ( 'U' )                                     ! U-point
( 4838)              DO ji = 1, jpiglo-1
( 4839)                 iju = jpiglo-ji+1
( 4840)                 ztab(ji,ijpj) = psgn * ztab(iju,ijpj-2)
( 4841)              END DO
( 4842)              DO ji = jpiglo/2, jpiglo-1
( 4843)                 iju = jpiglo-ji+1
( 4844)                 ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
( 4845)              END DO
( 4846) 
( 4847)           CASE ( 'V' )                                     ! V-point
( 4848)              DO ji = 2, jpiglo
( 4849)                 ijt = jpiglo-ji+2
( 4850)                 ztab(ji,ijpj-1) = psgn * ztab(ijt,ijpj-2)
( 4851)                 ztab(ji,ijpj  ) = psgn * ztab(ijt,ijpj-3)
( 4852)              END DO
( 4853) 
( 4854)           CASE ( 'F' , 'G' )                               ! F-point
( 4855)              DO ji = 1, jpiglo-1
( 4856)                 iju = jpiglo-ji+1
( 4857)                 ztab(ji,ijpj-1) = psgn * ztab(iju,ijpj-2)
( 4858)                 ztab(ji,ijpj  ) = psgn * ztab(iju,ijpj-3)
( 4859)              END DO
( 4860) 
( 4861)           CASE ( 'I' )                                     ! ice U-V point
( 4862)              ztab(2,ijpj) = psgn * ztab(3,ijpj-1)
( 4863)              DO ji = 3, jpiglo
( 4864)                 iju = jpiglo - ji + 3
( 4865)                 ztab(ji,ijpj) = psgn * ztab(iju,ijpj-1)
( 4866)              END DO
( 4867) 
( 4868)           END SELECT
( 4869) 
( 4870)        CASE ( 5, 6 )                        ! *  North fold  F-point pivot
( 4871) 
( 4872)           ztab( 1 ,ijpj) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 172

( 4873)           ztab(jpiglo,ijpj) = 0.e0
( 4874) 
( 4875)           SELECT CASE ( cd_type )
( 4876) 
( 4877)           CASE ( 'T' , 'W' ,'S' )                          ! T-, W-point
( 4878)              DO ji = 1, jpiglo
( 4879)                 ijt = jpiglo-ji+1
( 4880)                 ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-1)
( 4881)              END DO
( 4882) 
( 4883)           CASE ( 'U' )                                     ! U-point
( 4884)              DO ji = 1, jpiglo-1
( 4885)                 iju = jpiglo-ji
( 4886)                 ztab(ji,ijpj) = psgn * ztab(iju,ijpj-1)
( 4887)              END DO
( 4888) 
( 4889)           CASE ( 'V' )                                     ! V-point
( 4890)              DO ji = 1, jpiglo
( 4891)                 ijt = jpiglo-ji+1
( 4892)                 ztab(ji,ijpj) = psgn * ztab(ijt,ijpj-2)
( 4893)              END DO
( 4894)              DO ji = jpiglo/2+1, jpiglo
( 4895)                 ijt = jpiglo-ji+1
( 4896)                 ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
( 4897)              END DO
( 4898) 
( 4899)           CASE ( 'F' , 'G' )                               ! F-point
( 4900)              DO ji = 1, jpiglo-1
( 4901)                 iju = jpiglo-ji
( 4902)                 ztab(ji,ijpj  ) = psgn * ztab(iju,ijpj-2)
( 4903)              END DO
( 4904)              DO ji = jpiglo/2+1, jpiglo-1
( 4905)                 iju = jpiglo-ji
( 4906)                 ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
( 4907)              END DO
( 4908) 
( 4909)              CASE ( 'I' )                                  ! ice U-V point
( 4910)                 ztab( 2 ,ijpj) = 0.e0
( 4911)                 DO ji = 2 , jpiglo-1
( 4912)                    ijt = jpiglo - ji + 2
( 4913)                    ztab(ji,ijpj)= 0.5 * ( ztab(ji,ijpj-1) + psgn * ztab(ijt,ijpj-1) )
( 4914)                 END DO
( 4915) 
( 4916)           END SELECT
( 4917) 
( 4918)        CASE DEFAULT                           ! *  closed : the code probably never go through
( 4919) 
( 4920)             SELECT CASE ( cd_type) 
( 4921)   
( 4922)             CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
( 4923)                ztab(:, 1 ) = 0.e0
( 4924)                ztab(:,ijpj) = 0.e0
( 4925) 
( 4926)             CASE ( 'F' )                                   ! F-point
( 4927)                ztab(:,ijpj) = 0.e0
( 4928) 
( 4929)             CASE ( 'I' )                                   ! ice U-V point
( 4930)                ztab(:, 1 ) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 173

( 4931)                ztab(:,ijpj) = 0.e0
( 4932) 
( 4933)             END SELECT
( 4934) 
( 4935)          END SELECT
( 4936) 
( 4937)          !     End of slab
( 4938)          !     ===========
( 4939) 
( 4940)          !! Scatter back to pt2d
( 4941)          DO jr = 1, ndim_rank_north
( 4942)             jproc=nrank_north(jr)+1
( 4943)             ildi=nldit (jproc)
( 4944)             ilei=nleit (jproc)
( 4945)             iilb=nimppt(jproc)
( 4946)             DO jj=1,ijpj
( 4947)                DO ji=ildi,ilei
( 4948)                   znorthgloio(ji,jj,jr)=ztab(ji+iilb-1,jj)
( 4949)                END DO
( 4950)             END DO
( 4951)          END DO
( 4952) 
( 4953)       ENDIF      ! only done on proc 0 of ncomm_north
( 4954) 
( 4955) # 4958
( 4958)       IF ( npolj /= 0 ) THEN
( 4959)          itaille=jpi*ijpj
( 4960)          CALL MPI_SCATTER(znorthgloio,itaille,MPI_DOUBLE_PRECISION,znorthloc,itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 4961)       ENDIF
( 4962) 
( 4963) # 4964
( 4964)       ! put in the last ijpj jlines of pt2d znorthloc
( 4965)       DO jj = nlcj - ijpj + 1 , nlcj
( 4966)          ij = jj - nlcj + ijpj
( 4967)          pt2d(:,jj)= znorthloc(:,ij)
( 4968)       END DO
( 4969) 
( 4970)    END SUBROUTINE mpp_lbc_north_2d
( 4971) 
( 4972) 
( 4973)    SUBROUTINE mpp_lbc_north_e ( pt2d, cd_type, psgn)
( 4974)     !!---------------------------------------------------------------------
( 4975)     !!                   ***  routine mpp_lbc_north_2d  ***
( 4976)     !!
( 4977)     !! ** Purpose :
( 4978)     !!      Ensure proper north fold horizontal bondary condition in mpp configuration
( 4979)     !!      in case of jpn1 > 1 (for 2d array with outer extra halo)
( 4980)     !!
( 4981)     !! ** Method :
( 4982)     !!      Gather the 4+2*jpr2dj northern lines of the global domain on 1 processor and 
( 4983)     !!      apply lbc north-fold on this sub array. Then scatter the fold array 
( 4984)     !!      back to the processors.
( 4985)     !!
( 4986)     !! History :
( 4987)     !!   8.5  !  03-09  (J.M. Molines ) For mpp folding condition at north
( 4988)     !!                                  from lbc routine
( 4989)     !!   9.0  !  03-12  (J.M. Molines ) encapsulation into lib_mpp, coding rules of lbc_lnk
( 4990)     !!   9.0  !  05-09  (R. Benshila )   adapt mpp_lbc_north_2d 






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 174

( 4991)     !!----------------------------------------------------------------------
( 4992) 
( 4993)     !! * Arguments
( 4994)     CHARACTER(len=1), INTENT( in ) ::   &
( 4995)          cd_type       ! nature of pt2d grid-points
( 4996)     !             !   = T ,  U , V , F or W  gridpoints
( 4997)     REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj), INTENT( inout ) ::   &
( 4998)          pt2d          ! 2D array on which the boundary condition is applied
( 4999)     REAL(wp), INTENT( in ) ::   &
( 5000)          psgn          ! control of the sign change
( 5001)     !             !   = -1. , the sign is changed if north fold boundary
( 5002)     !             !   =  1. , the sign is kept  if north fold boundary
( 5003) 
( 5004) 
( 5005)     !! * Local declarations
( 5006) 
( 5007)     INTEGER :: ji, jj,  jr, jproc, jl
( 5008)     INTEGER :: ierr
( 5009)     INTEGER :: ildi,ilei,iilb
( 5010)     INTEGER :: ijpj,ijpjm1,ij,ijt,iju, iprecj
( 5011)     INTEGER :: itaille
( 5012) 
( 5013)     REAL(wp), DIMENSION(jpiglo,1-jpr2dj:4+jpr2dj) :: ztab
( 5014)     REAL(wp), DIMENSION(jpi,1-jpr2dj:4+jpr2dj,jpni) :: znorthgloio
( 5015)     REAL(wp), DIMENSION(jpi,1-jpr2dj:4+jpr2dj) :: znorthloc
( 5016) 
( 5017)     ! If we get in this routine it s because : North fold condition and mpp with more
( 5018)     !   than one proc across i : we deal only with the North condition
( 5019) 
( 5020)     ! 0. Sign setting
( 5021)     ! ---------------
( 5022) 
( 5023)     ijpj=4
( 5024)     ijpjm1=3
( 5025)     iprecj = jpr2dj+jprecj
( 5026) 
( 5027)     ! put in znorthloc the last 4 jlines of pt2d
( 5028)     DO jj = nlcj - ijpj + 1 - jpr2dj, nlcj +jpr2dj
( 5029)        ij = jj - nlcj + ijpj
( 5030)        znorthloc(:,ij)=pt2d(1:jpi,jj)
( 5031)     END DO
( 5032) 
( 5033)     IF (npolj /= 0 ) THEN
( 5034)        ! Build in proc 0 of ncomm_north the znorthgloio
( 5035)        znorthgloio(:,:,:) = 0_wp
( 5036) # 5039
( 5039)        itaille=jpi*(ijpj+2*jpr2dj)
( 5040)        CALL MPI_GATHER(znorthloc(1,1-jpr2dj),itaille,MPI_DOUBLE_PRECISION, &
( 5041)                      & znorthgloio(1,1-jpr2dj,1),itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)
( 5042) # 5043
( 5043)     ENDIF
( 5044) 
( 5045)     IF (narea == north_root+1 ) THEN
( 5046)        ! recover the global north array
( 5047)        ztab(:,:) = 0_wp
( 5048) 
( 5049)        DO jr = 1, ndim_rank_north
( 5050)           jproc=nrank_north(jr)+1






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 175

( 5051)           ildi=nldit (jproc)
( 5052)           ilei=nleit (jproc)
( 5053)           iilb=nimppt(jproc)
( 5054)           DO jj=1-jpr2dj,ijpj+jpr2dj
( 5055)              DO ji=ildi,ilei
( 5056)                 ztab(ji+iilb-1,jj)=znorthgloio(ji,jj,jr)
( 5057)              END DO
( 5058)           END DO
( 5059)        END DO
( 5060) 
( 5061) 
( 5062)        ! 2. North-Fold boundary conditions
( 5063)        ! ----------------------------------
( 5064) 
( 5065)        SELECT CASE ( npolj )
( 5066) 
( 5067)        CASE ( 3, 4 )                       ! *  North fold  T-point pivot
( 5068) 
( 5069)           ztab( 1    ,ijpj:ijpj+jpr2dj) = 0.e0
( 5070)           ztab(jpiglo,ijpj:ijpj+jpr2dj) = 0.e0
( 5071) 
( 5072)           SELECT CASE ( cd_type )
( 5073) 
( 5074)           CASE ( 'T' , 'W' , 'S' )                         ! T-, W-point
( 5075)              DO jl =0, iprecj-1
( 5076)                 DO ji = 2, jpiglo
( 5077)                    ijt = jpiglo-ji+2
( 5078)                    ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-2-jl)
( 5079)                 END DO
( 5080)              END DO
( 5081)              DO ji = jpiglo/2+1, jpiglo
( 5082)                 ijt = jpiglo-ji+2
( 5083)                 ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
( 5084)              END DO
( 5085) 
( 5086)           CASE ( 'U' )                                     ! U-point
( 5087)              DO jl =0, iprecj-1
( 5088)                 DO ji = 1, jpiglo-1
( 5089)                    iju = jpiglo-ji+1
( 5090)                    ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-2-jl)
( 5091)                 END DO
( 5092)              END DO
( 5093)              DO ji = jpiglo/2, jpiglo-1
( 5094)                 iju = jpiglo-ji+1
( 5095)                 ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
( 5096)              END DO
( 5097) 
( 5098)           CASE ( 'V' )                                     ! V-point
( 5099)             DO jl =-1, iprecj-1
( 5100)                DO ji = 2, jpiglo
( 5101)                   ijt = jpiglo-ji+2
( 5102)                   ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-3-jl)
( 5103)                END DO
( 5104)             END DO
( 5105) 
( 5106)           CASE ( 'F' , 'G' )                               ! F-point
( 5107)             DO jl =-1, iprecj-1
( 5108)                DO ji = 1, jpiglo-1






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 176

( 5109)                   iju = jpiglo-ji+1
( 5110)                   ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-3-jl)
( 5111)                END DO
( 5112)              END DO
( 5113) 
( 5114)           CASE ( 'I' )                                     ! ice U-V point
( 5115)              DO jl =0, iprecj-1
( 5116)                 ztab(2,ijpj+jl) = psgn * ztab(3,ijpj-1+jl)
( 5117)                 DO ji = 3, jpiglo
( 5118)                    iju = jpiglo - ji + 3
( 5119)                    ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-1-jl)
( 5120)                 END DO
( 5121)              END DO
( 5122) 
( 5123)           END SELECT
( 5124) 
( 5125)        CASE ( 5, 6 )                        ! *  North fold  F-point pivot
( 5126) 
( 5127)           ztab( 1 ,ijpj:ijpj+jpr2dj) = 0.e0
( 5128)           ztab(jpiglo,ijpj:ijpj+jpr2dj) = 0.e0
( 5129) 
( 5130)           SELECT CASE ( cd_type )
( 5131) 
( 5132)           CASE ( 'T' , 'W' ,'S' )                          ! T-, W-point
( 5133)              DO jl = 0, iprecj-1
( 5134)                 DO ji = 1, jpiglo
( 5135)                    ijt = jpiglo-ji+1
( 5136)                    ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-1-jl)
( 5137)                 END DO
( 5138)              END DO
( 5139) 
( 5140)           CASE ( 'U' )                                     ! U-point
( 5141)              DO jl = 0, iprecj-1
( 5142)                 DO ji = 1, jpiglo-1
( 5143)                    iju = jpiglo-ji
( 5144)                    ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-1-jl)
( 5145)                 END DO
( 5146)              END DO
( 5147) 
( 5148)           CASE ( 'V' )                                     ! V-point
( 5149)              DO jl = 0, iprecj-1
( 5150)                 DO ji = 1, jpiglo
( 5151)                    ijt = jpiglo-ji+1
( 5152)                    ztab(ji,ijpj+jl) = psgn * ztab(ijt,ijpj-2-jl)
( 5153)                 END DO
( 5154)              END DO
( 5155)              DO ji = jpiglo/2+1, jpiglo
( 5156)                 ijt = jpiglo-ji+1
( 5157)                 ztab(ji,ijpjm1) = psgn * ztab(ijt,ijpjm1)
( 5158)              END DO
( 5159) 
( 5160)           CASE ( 'F' , 'G' )                               ! F-point
( 5161)              DO jl = 0, iprecj-1
( 5162)                 DO ji = 1, jpiglo-1
( 5163)                    iju = jpiglo-ji
( 5164)                    ztab(ji,ijpj+jl) = psgn * ztab(iju,ijpj-2-jl)
( 5165)                 END DO
( 5166)              END DO






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 177

( 5167)              DO ji = jpiglo/2+1, jpiglo-1
( 5168)                 iju = jpiglo-ji
( 5169)                 ztab(ji,ijpjm1) = psgn * ztab(iju,ijpjm1)
( 5170)              END DO
( 5171) 
( 5172)              CASE ( 'I' )                                  ! ice U-V point
( 5173)                 ztab( 2 ,ijpj:ijpj+jpr2dj) = 0.e0
( 5174)                 DO jl = 0, jpr2dj
( 5175)                    DO ji = 2 , jpiglo-1
( 5176)                       ijt = jpiglo - ji + 2
( 5177)                       ztab(ji,ijpj+jl)= 0.5 * ( ztab(ji,ijpj-1-jl) + psgn * ztab(ijt,ijpj-1-jl) )
( 5178)                    END DO
( 5179)                 END DO
( 5180) 
( 5181)           END SELECT
( 5182) 
( 5183)        CASE DEFAULT                           ! *  closed : the code probably never go through
( 5184) 
( 5185)             SELECT CASE ( cd_type) 
( 5186)   
( 5187)             CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
( 5188)                ztab(:, 1:1-jpr2dj     ) = 0.e0
( 5189)                ztab(:,ijpj:ijpj+jpr2dj) = 0.e0
( 5190) 
( 5191)             CASE ( 'F' )                                   ! F-point
( 5192)                ztab(:,ijpj:ijpj+jpr2dj) = 0.e0
( 5193) 
( 5194)             CASE ( 'I' )                                   ! ice U-V point
( 5195)                ztab(:, 1:1-jpr2dj     ) = 0.e0
( 5196)                ztab(:,ijpj:ijpj+jpr2dj) = 0.e0
( 5197) 
( 5198)             END SELECT
( 5199) 
( 5200)          END SELECT
( 5201) 
( 5202)          !     End of slab
( 5203)          !     ===========
( 5204) 
( 5205)          !! Scatter back to pt2d
( 5206)          DO jr = 1, ndim_rank_north
( 5207)             jproc=nrank_north(jr)+1
( 5208)             ildi=nldit (jproc)
( 5209)             ilei=nleit (jproc)
( 5210)             iilb=nimppt(jproc)
( 5211)             DO jj=1-jpr2dj,ijpj+jpr2dj
( 5212)                DO ji=ildi,ilei
( 5213)                   znorthgloio(ji,jj,jr)=ztab(ji+iilb-1,jj)
( 5214)                END DO
( 5215)             END DO
( 5216)          END DO
( 5217) 
( 5218)       ENDIF      ! only done on proc 0 of ncomm_north
( 5219) 
( 5220) # 5223
( 5223)       IF ( npolj /= 0 ) THEN
( 5224)          itaille=jpi*(ijpj+2*jpr2dj)
( 5225)          CALL MPI_SCATTER(znorthgloio(1,1-jpr2dj,1),itaille,MPI_DOUBLE_PRECISION, &
( 5226)                         & znorthloc(1,1-jpr2dj),itaille,MPI_DOUBLE_PRECISION,0,ncomm_north,ierr)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 178

( 5227)       ENDIF
( 5228) 
( 5229) # 5230
( 5230)       ! put in the last ijpj jlines of pt2d znorthloc
( 5231)       DO jj = nlcj - ijpj  -jpr2dj + 1 , nlcj +jpr2dj
( 5232)          ij = jj - nlcj + ijpj 
( 5233)          pt2d(1:jpi,jj)= znorthloc(:,ij)
( 5234)       END DO
( 5235) 
( 5236)    END SUBROUTINE mpp_lbc_north_e
( 5237) 
( 5238) 
( 5239)    !!!!!
( 5240) 
( 5241) 
( 5242)    !! 
( 5243)    !!    This is valid on IBM machine ONLY. 
( 5244)    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
( 5245)    !! mpi_init_opa.f90 : Redefinition du point d'entree MPI_INIT de la bibliotheque
( 5246)    !!                MPI afin de faire, en plus de l'initialisation de
( 5247)    !!                l'environnement MPI, l'allocation d'une zone tampon
( 5248)    !!                qui sera ulterieurement utilisee automatiquement lors
( 5249)    !!                de tous les envois de messages par MPI_BSEND
( 5250)    !!
( 5251)    !! Auteur : CNRS/IDRIS
( 5252)    !! Date   : Tue Nov 13 12:02:14 2001
( 5253)    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
( 5254) 
( 5255)    SUBROUTINE mpi_init_opa(code)
( 5256)       IMPLICIT NONE
( 5257) 
( 5258) !$AGRIF_DO_NOT_TREAT
( 5259) # 1 "/usr/include/mpich2-x86_64/mpif.h"
(    1) !      /* -*- Mode: Fortran; -*- */
(    2) !      
(    3) !      (C) 2001 by Argonne National Laboratory.
(    4) !      See COPYRIGHT in top-level directory.
(    5) !      
(    6) !      DO NOT EDIT
(    7) !      This file created by buildiface 
(    8) !      
(    9)        INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
(   10)        PARAMETER (MPI_SOURCE=3,MPI_TAG=4,MPI_ERROR=5)
(   11)        INTEGER MPI_STATUS_SIZE
(   12)        PARAMETER (MPI_STATUS_SIZE=5)
(   13)        INTEGER MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
(   14)        INTEGER MPI_STATUSES_IGNORE(MPI_STATUS_SIZE,1)
(   15)        INTEGER MPI_ERRCODES_IGNORE(1)
(   16)        CHARACTER*1 MPI_ARGVS_NULL(1,1)
(   17)        CHARACTER*1 MPI_ARGV_NULL(1)
(   18)        INTEGER MPI_SUCCESS
(   19)        PARAMETER (MPI_SUCCESS=0)
(   20)        INTEGER MPI_ERR_OTHER
(   21)        PARAMETER (MPI_ERR_OTHER=15)
(   22)        INTEGER MPI_ERR_WIN
(   23)        PARAMETER (MPI_ERR_WIN=45)
(   24)        INTEGER MPI_ERR_FILE
(   25)        PARAMETER (MPI_ERR_FILE=27)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 179

(   26)        INTEGER MPI_ERR_COUNT
(   27)        PARAMETER (MPI_ERR_COUNT=2)
(   28)        INTEGER MPI_ERR_SPAWN
(   29)        PARAMETER (MPI_ERR_SPAWN=42)
(   30)        INTEGER MPI_ERR_BASE
(   31)        PARAMETER (MPI_ERR_BASE=46)
(   32)        INTEGER MPI_ERR_RMA_CONFLICT
(   33)        PARAMETER (MPI_ERR_RMA_CONFLICT=49)
(   34)        INTEGER MPI_ERR_IN_STATUS
(   35)        PARAMETER (MPI_ERR_IN_STATUS=17)
(   36)        INTEGER MPI_ERR_INFO_KEY
(   37)        PARAMETER (MPI_ERR_INFO_KEY=29)
(   38)        INTEGER MPI_ERR_LOCKTYPE
(   39)        PARAMETER (MPI_ERR_LOCKTYPE=47)
(   40)        INTEGER MPI_ERR_OP
(   41)        PARAMETER (MPI_ERR_OP=9)
(   42)        INTEGER MPI_ERR_ARG
(   43)        PARAMETER (MPI_ERR_ARG=12)
(   44)        INTEGER MPI_ERR_READ_ONLY
(   45)        PARAMETER (MPI_ERR_READ_ONLY=40)
(   46)        INTEGER MPI_ERR_SIZE
(   47)        PARAMETER (MPI_ERR_SIZE=51)
(   48)        INTEGER MPI_ERR_BUFFER
(   49)        PARAMETER (MPI_ERR_BUFFER=1)
(   50)        INTEGER MPI_ERR_DUP_DATAREP
(   51)        PARAMETER (MPI_ERR_DUP_DATAREP=24)
(   52)        INTEGER MPI_ERR_UNSUPPORTED_DATAREP
(   53)        PARAMETER (MPI_ERR_UNSUPPORTED_DATAREP=43)
(   54)        INTEGER MPI_ERR_LASTCODE
(   55)        PARAMETER (MPI_ERR_LASTCODE=1073741823)
(   56)        INTEGER MPI_ERR_TRUNCATE
(   57)        PARAMETER (MPI_ERR_TRUNCATE=14)
(   58)        INTEGER MPI_ERR_DISP
(   59)        PARAMETER (MPI_ERR_DISP=52)
(   60)        INTEGER MPI_ERR_PORT
(   61)        PARAMETER (MPI_ERR_PORT=38)
(   62)        INTEGER MPI_ERR_INFO_NOKEY
(   63)        PARAMETER (MPI_ERR_INFO_NOKEY=31)
(   64)        INTEGER MPI_ERR_ASSERT
(   65)        PARAMETER (MPI_ERR_ASSERT=53)
(   66)        INTEGER MPI_ERR_FILE_EXISTS
(   67)        PARAMETER (MPI_ERR_FILE_EXISTS=25)
(   68)        INTEGER MPI_ERR_PENDING
(   69)        PARAMETER (MPI_ERR_PENDING=18)
(   70)        INTEGER MPI_ERR_COMM
(   71)        PARAMETER (MPI_ERR_COMM=5)
(   72)        INTEGER MPI_ERR_KEYVAL
(   73)        PARAMETER (MPI_ERR_KEYVAL=48)
(   74)        INTEGER MPI_ERR_NAME
(   75)        PARAMETER (MPI_ERR_NAME=33)
(   76)        INTEGER MPI_ERR_REQUEST
(   77)        PARAMETER (MPI_ERR_REQUEST=19)
(   78)        INTEGER MPI_ERR_GROUP
(   79)        PARAMETER (MPI_ERR_GROUP=8)
(   80)        INTEGER MPI_ERR_TOPOLOGY
(   81)        PARAMETER (MPI_ERR_TOPOLOGY=10)
(   82)        INTEGER MPI_ERR_TYPE
(   83)        PARAMETER (MPI_ERR_TYPE=3)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 180

(   84)        INTEGER MPI_ERR_TAG
(   85)        PARAMETER (MPI_ERR_TAG=4)
(   86)        INTEGER MPI_ERR_INFO_VALUE
(   87)        PARAMETER (MPI_ERR_INFO_VALUE=30)
(   88)        INTEGER MPI_ERR_NOT_SAME
(   89)        PARAMETER (MPI_ERR_NOT_SAME=35)
(   90)        INTEGER MPI_ERR_RMA_SYNC
(   91)        PARAMETER (MPI_ERR_RMA_SYNC=50)
(   92)        INTEGER MPI_ERR_INFO
(   93)        PARAMETER (MPI_ERR_INFO=28)
(   94)        INTEGER MPI_ERR_NO_MEM
(   95)        PARAMETER (MPI_ERR_NO_MEM=34)
(   96)        INTEGER MPI_ERR_BAD_FILE
(   97)        PARAMETER (MPI_ERR_BAD_FILE=22)
(   98)        INTEGER MPI_ERR_FILE_IN_USE
(   99)        PARAMETER (MPI_ERR_FILE_IN_USE=26)
(  100)        INTEGER MPI_ERR_UNKNOWN
(  101)        PARAMETER (MPI_ERR_UNKNOWN=13)
(  102)        INTEGER MPI_ERR_UNSUPPORTED_OPERATION
(  103)        PARAMETER (MPI_ERR_UNSUPPORTED_OPERATION=44)
(  104)        INTEGER MPI_ERR_QUOTA
(  105)        PARAMETER (MPI_ERR_QUOTA=39)
(  106)        INTEGER MPI_ERR_AMODE
(  107)        PARAMETER (MPI_ERR_AMODE=21)
(  108)        INTEGER MPI_ERR_ROOT
(  109)        PARAMETER (MPI_ERR_ROOT=7)
(  110)        INTEGER MPI_ERR_RANK
(  111)        PARAMETER (MPI_ERR_RANK=6)
(  112)        INTEGER MPI_ERR_DIMS
(  113)        PARAMETER (MPI_ERR_DIMS=11)
(  114)        INTEGER MPI_ERR_NO_SUCH_FILE
(  115)        PARAMETER (MPI_ERR_NO_SUCH_FILE=37)
(  116)        INTEGER MPI_ERR_SERVICE
(  117)        PARAMETER (MPI_ERR_SERVICE=41)
(  118)        INTEGER MPI_ERR_INTERN
(  119)        PARAMETER (MPI_ERR_INTERN=16)
(  120)        INTEGER MPI_ERR_IO
(  121)        PARAMETER (MPI_ERR_IO=32)
(  122)        INTEGER MPI_ERR_ACCESS
(  123)        PARAMETER (MPI_ERR_ACCESS=20)
(  124)        INTEGER MPI_ERR_NO_SPACE
(  125)        PARAMETER (MPI_ERR_NO_SPACE=36)
(  126)        INTEGER MPI_ERR_CONVERSION
(  127)        PARAMETER (MPI_ERR_CONVERSION=23)
(  128)        INTEGER MPI_ERRORS_ARE_FATAL
(  129)        PARAMETER (MPI_ERRORS_ARE_FATAL=1409286144)
(  130)        INTEGER MPI_ERRORS_RETURN
(  131)        PARAMETER (MPI_ERRORS_RETURN=1409286145)
(  132)        INTEGER MPI_IDENT
(  133)        PARAMETER (MPI_IDENT=0)
(  134)        INTEGER MPI_CONGRUENT
(  135)        PARAMETER (MPI_CONGRUENT=1)
(  136)        INTEGER MPI_SIMILAR
(  137)        PARAMETER (MPI_SIMILAR=2)
(  138)        INTEGER MPI_UNEQUAL
(  139)        PARAMETER (MPI_UNEQUAL=3)
(  140)        INTEGER MPI_MAX
(  141)        PARAMETER (MPI_MAX=1476395009)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 181

(  142)        INTEGER MPI_MIN
(  143)        PARAMETER (MPI_MIN=1476395010)
(  144)        INTEGER MPI_SUM
(  145)        PARAMETER (MPI_SUM=1476395011)
(  146)        INTEGER MPI_PROD
(  147)        PARAMETER (MPI_PROD=1476395012)
(  148)        INTEGER MPI_LAND
(  149)        PARAMETER (MPI_LAND=1476395013)
(  150)        INTEGER MPI_BAND
(  151)        PARAMETER (MPI_BAND=1476395014)
(  152)        INTEGER MPI_LOR
(  153)        PARAMETER (MPI_LOR=1476395015)
(  154)        INTEGER MPI_BOR
(  155)        PARAMETER (MPI_BOR=1476395016)
(  156)        INTEGER MPI_LXOR
(  157)        PARAMETER (MPI_LXOR=1476395017)
(  158)        INTEGER MPI_BXOR
(  159)        PARAMETER (MPI_BXOR=1476395018)
(  160)        INTEGER MPI_MINLOC
(  161)        PARAMETER (MPI_MINLOC=1476395019)
(  162)        INTEGER MPI_MAXLOC
(  163)        PARAMETER (MPI_MAXLOC=1476395020)
(  164)        INTEGER MPI_REPLACE
(  165)        PARAMETER (MPI_REPLACE=1476395021)
(  166)        INTEGER MPI_COMM_WORLD
(  167)        PARAMETER (MPI_COMM_WORLD=1140850688)
(  168)        INTEGER MPI_COMM_SELF
(  169)        PARAMETER (MPI_COMM_SELF=1140850689)
(  170)        INTEGER MPI_GROUP_EMPTY
(  171)        PARAMETER (MPI_GROUP_EMPTY=1207959552)
(  172)        INTEGER MPI_COMM_NULL
(  173)        PARAMETER (MPI_COMM_NULL=67108864)
(  174)        INTEGER MPI_WIN_NULL
(  175)        PARAMETER (MPI_WIN_NULL=536870912)
(  176)        INTEGER MPI_FILE_NULL
(  177)        PARAMETER (MPI_FILE_NULL=0)
(  178)        INTEGER MPI_GROUP_NULL
(  179)        PARAMETER (MPI_GROUP_NULL=134217728)
(  180)        INTEGER MPI_OP_NULL
(  181)        PARAMETER (MPI_OP_NULL=402653184)
(  182)        INTEGER MPI_DATATYPE_NULL
(  183)        PARAMETER (MPI_DATATYPE_NULL=201326592)
(  184)        INTEGER MPI_REQUEST_NULL
(  185)        PARAMETER (MPI_REQUEST_NULL=738197504)
(  186)        INTEGER MPI_ERRHANDLER_NULL
(  187)        PARAMETER (MPI_ERRHANDLER_NULL=335544320)
(  188)        INTEGER MPI_INFO_NULL
(  189)        PARAMETER (MPI_INFO_NULL=469762048)
(  190)        INTEGER MPI_TAG_UB
(  191)        PARAMETER (MPI_TAG_UB=1681915906)
(  192)        INTEGER MPI_HOST
(  193)        PARAMETER (MPI_HOST=1681915908)
(  194)        INTEGER MPI_IO
(  195)        PARAMETER (MPI_IO=1681915910)
(  196)        INTEGER MPI_WTIME_IS_GLOBAL
(  197)        PARAMETER (MPI_WTIME_IS_GLOBAL=1681915912)
(  198)        INTEGER MPI_UNIVERSE_SIZE
(  199)        PARAMETER (MPI_UNIVERSE_SIZE=1681915914)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 182

(  200)        INTEGER MPI_LASTUSEDCODE
(  201)        PARAMETER (MPI_LASTUSEDCODE=1681915916)
(  202)        INTEGER MPI_APPNUM
(  203)        PARAMETER (MPI_APPNUM=1681915918)
(  204)        INTEGER MPI_WIN_BASE
(  205)        PARAMETER (MPI_WIN_BASE=1711276034)
(  206)        INTEGER MPI_WIN_SIZE
(  207)        PARAMETER (MPI_WIN_SIZE=1711276036)
(  208)        INTEGER MPI_WIN_DISP_UNIT
(  209)        PARAMETER (MPI_WIN_DISP_UNIT=1711276038)
(  210)        INTEGER MPI_MAX_ERROR_STRING
(  211)        PARAMETER (MPI_MAX_ERROR_STRING=1023)
(  212)        INTEGER MPI_MAX_PORT_NAME
(  213)        PARAMETER (MPI_MAX_PORT_NAME=255)
(  214)        INTEGER MPI_MAX_OBJECT_NAME
(  215)        PARAMETER (MPI_MAX_OBJECT_NAME=127)
(  216)        INTEGER MPI_MAX_INFO_KEY
(  217)        PARAMETER (MPI_MAX_INFO_KEY=254)
(  218)        INTEGER MPI_MAX_INFO_VAL
(  219)        PARAMETER (MPI_MAX_INFO_VAL=1023)
(  220)        INTEGER MPI_MAX_PROCESSOR_NAME
(  221)        PARAMETER (MPI_MAX_PROCESSOR_NAME=128-1)
(  222)        INTEGER MPI_MAX_DATAREP_STRING
(  223)        PARAMETER (MPI_MAX_DATAREP_STRING=127)
(  224)        INTEGER MPI_UNDEFINED
(  225)        PARAMETER (MPI_UNDEFINED=(-32766))
(  226)        INTEGER MPI_KEYVAL_INVALID
(  227)        PARAMETER (MPI_KEYVAL_INVALID=603979776)
(  228)        INTEGER MPI_BSEND_OVERHEAD
(  229)        PARAMETER (MPI_BSEND_OVERHEAD=88)
(  230)        INTEGER MPI_PROC_NULL
(  231)        PARAMETER (MPI_PROC_NULL=-1)
(  232)        INTEGER MPI_ANY_SOURCE
(  233)        PARAMETER (MPI_ANY_SOURCE=-2)
(  234)        INTEGER MPI_ANY_TAG
(  235)        PARAMETER (MPI_ANY_TAG=-1)
(  236)        INTEGER MPI_ROOT
(  237)        PARAMETER (MPI_ROOT=-3)
(  238)        INTEGER MPI_GRAPH
(  239)        PARAMETER (MPI_GRAPH=1)
(  240)        INTEGER MPI_CART
(  241)        PARAMETER (MPI_CART=2)
(  242)        INTEGER MPI_VERSION
(  243)        PARAMETER (MPI_VERSION=2)
(  244)        INTEGER MPI_SUBVERSION
(  245)        PARAMETER (MPI_SUBVERSION=2)
(  246)        INTEGER MPI_LOCK_EXCLUSIVE
(  247)        PARAMETER (MPI_LOCK_EXCLUSIVE=234)
(  248)        INTEGER MPI_LOCK_SHARED
(  249)        PARAMETER (MPI_LOCK_SHARED=235)
(  250)        INTEGER MPI_COMPLEX
(  251)        PARAMETER (MPI_COMPLEX=1275070494)
(  252)        INTEGER MPI_DOUBLE_COMPLEX
(  253)        PARAMETER (MPI_DOUBLE_COMPLEX=1275072546)
(  254)        INTEGER MPI_LOGICAL
(  255)        PARAMETER (MPI_LOGICAL=1275069469)
(  256)        INTEGER MPI_REAL
(  257)        PARAMETER (MPI_REAL=1275069468)






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 183

(  258)        INTEGER MPI_DOUBLE_PRECISION
(  259)        PARAMETER (MPI_DOUBLE_PRECISION=1275070495)
(  260)        INTEGER MPI_INTEGER
(  261)        PARAMETER (MPI_INTEGER=1275069467)
(  262)        INTEGER MPI_2INTEGER
(  263)        PARAMETER (MPI_2INTEGER=1275070496)
(  264)        INTEGER MPI_2COMPLEX
(  265)        PARAMETER (MPI_2COMPLEX=1275072548)
(  266)        INTEGER MPI_2DOUBLE_PRECISION
(  267)        PARAMETER (MPI_2DOUBLE_PRECISION=1275072547)
(  268)        INTEGER MPI_2REAL
(  269)        PARAMETER (MPI_2REAL=1275070497)
(  270)        INTEGER MPI_2DOUBLE_COMPLEX
(  271)        PARAMETER (MPI_2DOUBLE_COMPLEX=1275076645)
(  272)        INTEGER MPI_CHARACTER
(  273)        PARAMETER (MPI_CHARACTER=1275068698)
(  274)        INTEGER MPI_BYTE
(  275)        PARAMETER (MPI_BYTE=1275068685)
(  276)        INTEGER MPI_UB
(  277)        PARAMETER (MPI_UB=1275068433)
(  278)        INTEGER MPI_LB
(  279)        PARAMETER (MPI_LB=1275068432)
(  280)        INTEGER MPI_PACKED
(  281)        PARAMETER (MPI_PACKED=1275068687)
(  282)        INTEGER MPI_INTEGER1
(  283)        PARAMETER (MPI_INTEGER1=1275068717)
(  284)        INTEGER MPI_INTEGER2
(  285)        PARAMETER (MPI_INTEGER2=1275068975)
(  286)        INTEGER MPI_INTEGER4
(  287)        PARAMETER (MPI_INTEGER4=1275069488)
(  288)        INTEGER MPI_INTEGER8
(  289)        PARAMETER (MPI_INTEGER8=1275070513)
(  290)        INTEGER MPI_INTEGER16
(  291)        PARAMETER (MPI_INTEGER16=MPI_DATATYPE_NULL)
(  292)        INTEGER MPI_REAL4
(  293)        PARAMETER (MPI_REAL4=1275069479)
(  294)        INTEGER MPI_REAL8
(  295)        PARAMETER (MPI_REAL8=1275070505)
(  296)        INTEGER MPI_REAL16
(  297)        PARAMETER (MPI_REAL16=1275072555)
(  298)        INTEGER MPI_COMPLEX8
(  299)        PARAMETER (MPI_COMPLEX8=1275070504)
(  300)        INTEGER MPI_COMPLEX16
(  301)        PARAMETER (MPI_COMPLEX16=1275072554)
(  302)        INTEGER MPI_COMPLEX32
(  303)        PARAMETER (MPI_COMPLEX32=1275076652)
(  304)        INTEGER MPI_ADDRESS_KIND, MPI_OFFSET_KIND
(  305)        PARAMETER (MPI_ADDRESS_KIND=8)
(  306)        PARAMETER (MPI_OFFSET_KIND=8)
(  307)        INTEGER MPI_CHAR
(  308)        PARAMETER (MPI_CHAR=1275068673)
(  309)        INTEGER MPI_SIGNED_CHAR
(  310)        PARAMETER (MPI_SIGNED_CHAR=1275068696)
(  311)        INTEGER MPI_UNSIGNED_CHAR
(  312)        PARAMETER (MPI_UNSIGNED_CHAR=1275068674)
(  313)        INTEGER MPI_WCHAR
(  314)        PARAMETER (MPI_WCHAR=1275069454)
(  315)        INTEGER MPI_SHORT






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 184

(  316)        PARAMETER (MPI_SHORT=1275068931)
(  317)        INTEGER MPI_UNSIGNED_SHORT
(  318)        PARAMETER (MPI_UNSIGNED_SHORT=1275068932)
(  319)        INTEGER MPI_INT
(  320)        PARAMETER (MPI_INT=1275069445)
(  321)        INTEGER MPI_UNSIGNED
(  322)        PARAMETER (MPI_UNSIGNED=1275069446)
(  323)        INTEGER MPI_LONG
(  324)        PARAMETER (MPI_LONG=1275070471)
(  325)        INTEGER MPI_UNSIGNED_LONG
(  326)        PARAMETER (MPI_UNSIGNED_LONG=1275070472)
(  327)        INTEGER MPI_FLOAT
(  328)        PARAMETER (MPI_FLOAT=1275069450)
(  329)        INTEGER MPI_DOUBLE
(  330)        PARAMETER (MPI_DOUBLE=1275070475)
(  331)        INTEGER MPI_LONG_DOUBLE
(  332)        PARAMETER (MPI_LONG_DOUBLE=1275072524)
(  333)        INTEGER MPI_LONG_LONG_INT
(  334)        PARAMETER (MPI_LONG_LONG_INT=1275070473)
(  335)        INTEGER MPI_UNSIGNED_LONG_LONG
(  336)        PARAMETER (MPI_UNSIGNED_LONG_LONG=1275070489)
(  337)        INTEGER MPI_LONG_LONG
(  338)        PARAMETER (MPI_LONG_LONG=1275070473)
(  339)        INTEGER MPI_FLOAT_INT
(  340)        PARAMETER (MPI_FLOAT_INT=-1946157056)
(  341)        INTEGER MPI_DOUBLE_INT
(  342)        PARAMETER (MPI_DOUBLE_INT=-1946157055)
(  343)        INTEGER MPI_LONG_INT
(  344)        PARAMETER (MPI_LONG_INT=-1946157054)
(  345)        INTEGER MPI_SHORT_INT
(  346)        PARAMETER (MPI_SHORT_INT=-1946157053)
(  347)        INTEGER MPI_2INT
(  348)        PARAMETER (MPI_2INT=1275070486)
(  349)        INTEGER MPI_LONG_DOUBLE_INT
(  350)        PARAMETER (MPI_LONG_DOUBLE_INT=-1946157052)
(  351)        INTEGER MPI_INT8_T
(  352)        PARAMETER (MPI_INT8_T=1275068727)
(  353)        INTEGER MPI_INT16_T
(  354)        PARAMETER (MPI_INT16_T=1275068984)
(  355)        INTEGER MPI_INT32_T
(  356)        PARAMETER (MPI_INT32_T=1275069497)
(  357)        INTEGER MPI_INT64_T
(  358)        PARAMETER (MPI_INT64_T=1275070522)
(  359)        INTEGER MPI_UINT8_T
(  360)        PARAMETER (MPI_UINT8_T=1275068731)
(  361)        INTEGER MPI_UINT16_T
(  362)        PARAMETER (MPI_UINT16_T=1275068988)
(  363)        INTEGER MPI_UINT32_T
(  364)        PARAMETER (MPI_UINT32_T=1275069501)
(  365)        INTEGER MPI_UINT64_T
(  366)        PARAMETER (MPI_UINT64_T=1275070526)
(  367)        INTEGER MPI_C_BOOL
(  368)        PARAMETER (MPI_C_BOOL=1275068735)
(  369)        INTEGER MPI_C_FLOAT_COMPLEX
(  370)        PARAMETER (MPI_C_FLOAT_COMPLEX=1275070528)
(  371)        INTEGER MPI_C_COMPLEX
(  372)        PARAMETER (MPI_C_COMPLEX=1275070528)
(  373)        INTEGER MPI_C_DOUBLE_COMPLEX






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 185

(  374)        PARAMETER (MPI_C_DOUBLE_COMPLEX=1275072577)
(  375)        INTEGER MPI_C_LONG_DOUBLE_COMPLEX
(  376)        PARAMETER (MPI_C_LONG_DOUBLE_COMPLEX=1275076674)
(  377)        INTEGER MPI_AINT
(  378)        PARAMETER (MPI_AINT=1275070531)
(  379)        INTEGER MPI_OFFSET
(  380)        PARAMETER (MPI_OFFSET=1275070532)
(  381)        INTEGER MPI_COMBINER_NAMED
(  382)        PARAMETER (MPI_COMBINER_NAMED=1)
(  383)        INTEGER MPI_COMBINER_DUP
(  384)        PARAMETER (MPI_COMBINER_DUP=2)
(  385)        INTEGER MPI_COMBINER_CONTIGUOUS
(  386)        PARAMETER (MPI_COMBINER_CONTIGUOUS=3)
(  387)        INTEGER MPI_COMBINER_VECTOR
(  388)        PARAMETER (MPI_COMBINER_VECTOR=4)
(  389)        INTEGER MPI_COMBINER_HVECTOR_INTEGER
(  390)        PARAMETER (MPI_COMBINER_HVECTOR_INTEGER=5)
(  391)        INTEGER MPI_COMBINER_HVECTOR
(  392)        PARAMETER (MPI_COMBINER_HVECTOR=6)
(  393)        INTEGER MPI_COMBINER_INDEXED
(  394)        PARAMETER (MPI_COMBINER_INDEXED=7)
(  395)        INTEGER MPI_COMBINER_HINDEXED_INTEGER
(  396)        PARAMETER (MPI_COMBINER_HINDEXED_INTEGER=8)
(  397)        INTEGER MPI_COMBINER_HINDEXED
(  398)        PARAMETER (MPI_COMBINER_HINDEXED=9)
(  399)        INTEGER MPI_COMBINER_INDEXED_BLOCK
(  400)        PARAMETER (MPI_COMBINER_INDEXED_BLOCK=10)
(  401)        INTEGER MPI_COMBINER_STRUCT_INTEGER
(  402)        PARAMETER (MPI_COMBINER_STRUCT_INTEGER=11)
(  403)        INTEGER MPI_COMBINER_STRUCT
(  404)        PARAMETER (MPI_COMBINER_STRUCT=12)
(  405)        INTEGER MPI_COMBINER_SUBARRAY
(  406)        PARAMETER (MPI_COMBINER_SUBARRAY=13)
(  407)        INTEGER MPI_COMBINER_DARRAY
(  408)        PARAMETER (MPI_COMBINER_DARRAY=14)
(  409)        INTEGER MPI_COMBINER_F90_REAL
(  410)        PARAMETER (MPI_COMBINER_F90_REAL=15)
(  411)        INTEGER MPI_COMBINER_F90_COMPLEX
(  412)        PARAMETER (MPI_COMBINER_F90_COMPLEX=16)
(  413)        INTEGER MPI_COMBINER_F90_INTEGER
(  414)        PARAMETER (MPI_COMBINER_F90_INTEGER=17)
(  415)        INTEGER MPI_COMBINER_RESIZED
(  416)        PARAMETER (MPI_COMBINER_RESIZED=18)
(  417)        INTEGER MPI_TYPECLASS_REAL
(  418)        PARAMETER (MPI_TYPECLASS_REAL=1)
(  419)        INTEGER MPI_TYPECLASS_INTEGER
(  420)        PARAMETER (MPI_TYPECLASS_INTEGER=2)
(  421)        INTEGER MPI_TYPECLASS_COMPLEX
(  422)        PARAMETER (MPI_TYPECLASS_COMPLEX=3)
(  423)        INTEGER MPI_MODE_NOCHECK
(  424)        PARAMETER (MPI_MODE_NOCHECK=1024)
(  425)        INTEGER MPI_MODE_NOSTORE
(  426)        PARAMETER (MPI_MODE_NOSTORE=2048)
(  427)        INTEGER MPI_MODE_NOPUT
(  428)        PARAMETER (MPI_MODE_NOPUT=4096)
(  429)        INTEGER MPI_MODE_NOPRECEDE
(  430)        PARAMETER (MPI_MODE_NOPRECEDE=8192)
(  431)        INTEGER MPI_MODE_NOSUCCEED






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 186

(  432)        PARAMETER (MPI_MODE_NOSUCCEED=16384)
(  433)        INTEGER MPI_THREAD_SINGLE
(  434)        PARAMETER (MPI_THREAD_SINGLE=0)
(  435)        INTEGER MPI_THREAD_FUNNELED
(  436)        PARAMETER (MPI_THREAD_FUNNELED=1)
(  437)        INTEGER MPI_THREAD_SERIALIZED
(  438)        PARAMETER (MPI_THREAD_SERIALIZED=2)
(  439)        INTEGER MPI_THREAD_MULTIPLE
(  440)        PARAMETER (MPI_THREAD_MULTIPLE=3)
(  441)        INTEGER MPI_MODE_RDONLY
(  442)        PARAMETER (MPI_MODE_RDONLY=2)
(  443)        INTEGER MPI_MODE_RDWR
(  444)        PARAMETER (MPI_MODE_RDWR=8)
(  445)        INTEGER MPI_MODE_WRONLY
(  446)        PARAMETER (MPI_MODE_WRONLY=4)
(  447)        INTEGER MPI_MODE_DELETE_ON_CLOSE
(  448)        PARAMETER (MPI_MODE_DELETE_ON_CLOSE=16)
(  449)        INTEGER MPI_MODE_UNIQUE_OPEN
(  450)        PARAMETER (MPI_MODE_UNIQUE_OPEN=32)
(  451)        INTEGER MPI_MODE_CREATE
(  452)        PARAMETER (MPI_MODE_CREATE=1)
(  453)        INTEGER MPI_MODE_EXCL
(  454)        PARAMETER (MPI_MODE_EXCL=64)
(  455)        INTEGER MPI_MODE_APPEND
(  456)        PARAMETER (MPI_MODE_APPEND=128)
(  457)        INTEGER MPI_MODE_SEQUENTIAL
(  458)        PARAMETER (MPI_MODE_SEQUENTIAL=256)
(  459)        INTEGER MPI_SEEK_SET
(  460)        PARAMETER (MPI_SEEK_SET=600)
(  461)        INTEGER MPI_SEEK_CUR
(  462)        PARAMETER (MPI_SEEK_CUR=602)
(  463)        INTEGER MPI_SEEK_END
(  464)        PARAMETER (MPI_SEEK_END=604)
(  465)        INTEGER MPI_ORDER_C
(  466)        PARAMETER (MPI_ORDER_C=56)
(  467)        INTEGER MPI_ORDER_FORTRAN
(  468)        PARAMETER (MPI_ORDER_FORTRAN=57)
(  469)        INTEGER MPI_DISTRIBUTE_BLOCK
(  470)        PARAMETER (MPI_DISTRIBUTE_BLOCK=121)
(  471)        INTEGER MPI_DISTRIBUTE_CYCLIC
(  472)        PARAMETER (MPI_DISTRIBUTE_CYCLIC=122)
(  473)        INTEGER MPI_DISTRIBUTE_NONE
(  474)        PARAMETER (MPI_DISTRIBUTE_NONE=123)
(  475)        INTEGER MPI_DISTRIBUTE_DFLT_DARG
(  476)        PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
(  477)        integer*8 MPI_DISPLACEMENT_CURRENT
(  478)        PARAMETER (MPI_DISPLACEMENT_CURRENT=-54278278)
(  479)        INTEGER MPI_BOTTOM, MPI_IN_PLACE
(  480)        EXTERNAL MPI_DUP_FN, MPI_NULL_DELETE_FN, MPI_NULL_COPY_FN
(  481)        EXTERNAL MPI_WTIME, MPI_WTICK
(  482)        EXTERNAL PMPI_WTIME, PMPI_WTICK
(  483)        EXTERNAL MPI_COMM_DUP_FN, MPI_COMM_NULL_DELETE_FN
(  484)        EXTERNAL MPI_COMM_NULL_COPY_FN
(  485)        EXTERNAL MPI_WIN_DUP_FN, MPI_WIN_NULL_DELETE_FN
(  486)        EXTERNAL MPI_WIN_NULL_COPY_FN
(  487)        EXTERNAL MPI_TYPE_DUP_FN, MPI_TYPE_NULL_DELETE_FN
(  488)        EXTERNAL MPI_TYPE_NULL_COPY_FN
(  489)        EXTERNAL MPI_CONVERSION_FN_NULL






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 187

(  490)        REAL*8 MPI_WTIME, MPI_WTICK
(  491)        REAL*8 PMPI_WTIME, PMPI_WTICK
(  492) 
(  493) 
(  494)        COMMON /MPIPRIV1/ MPI_BOTTOM, MPI_IN_PLACE, MPI_STATUS_IGNORE
(  495) 
(  496)        COMMON /MPIPRIV2/ MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE
(  497)        SAVE /MPIPRIV1/,/MPIPRIV2/
(  498) 
(  499)        COMMON /MPIPRIVC/ MPI_ARGVS_NULL, MPI_ARGV_NULL
(  500)        SAVE   /MPIPRIVC/
(  501) # 5260 "lib_mpp.F90"
( 5260) # 5260
( 5260) !$AGRIF_END_DO_NOT_TREAT
( 5261) 
( 5262)       INTEGER                                 :: code,rang,ierr
( 5263)       LOGICAL                                 :: mpi_was_called
( 5264)  
( 5265)       ! La valeur suivante doit etre au moins egale a la taille
( 5266)       ! du plus grand message qui sera transfere dans le programme
( 5267)       ! (de toute facon, il y aura un message d'erreur si cette
( 5268)       ! valeur s'avere trop petite)
( 5269)       INTEGER                                 :: taille_tampon
( 5270)       CHARACTER(len=9)                        :: taille_tampon_alphanum
( 5271)       REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: tampon
( 5272)  
( 5273)       ! Le point d'entree dans la bibliotheque MPI elle-meme
( 5274)       CALL mpi_initialized(mpi_was_called, code)
( 5275)       IF ( code /= MPI_SUCCESS ) THEN
( 5276)         CALL ctl_stop( ' lib_mpp: Error in routine mpi_initialized' )
( 5277)         CALL mpi_abort( mpi_comm_world, code, ierr )
( 5278)       ENDIF
( 5279) 
( 5280)       IF ( .NOT. mpi_was_called ) THEN
( 5281)          CALL mpi_init(code)
( 5282)          CALL mpi_comm_dup( mpi_comm_world, mpi_comm_opa, code)
( 5283)          IF ( code /= MPI_SUCCESS ) THEN
( 5284)             CALL ctl_stop( ' lib_mpp: Error in routine mpi_comm_dup' )
( 5285)             CALL mpi_abort( mpi_comm_world, code, ierr )
( 5286)          ENDIF
( 5287)       ENDIF
( 5288)       ! La definition de la zone tampon pour les futurs envois
( 5289)       ! par MPI_BSEND (on alloue une fois pour toute cette zone
( 5290)       ! tampon, qui sera automatiquement utilisee lors de chaque
( 5291)       ! appel  a MPI_BSEND).
( 5292)       ! La desallocation sera implicite quand on sortira de
( 5293)       ! l'environnement MPI.
( 5294) 
( 5295)       ! Recuperation de la valeur de la variable d'environnement
( 5296)       ! BUFFER_LENGTH
( 5297)       ! qui, si elle est definie, doit contenir une valeur superieure
( 5298)       ! a  la taille en octets du plus gros message
( 5299)       CALL getenv('BUFFER_LENGTH',taille_tampon_alphanum)
( 5300)   
( 5301)       ! Si la variable BUFFER_LENGTH n'est pas positionnee, on lui met par
( 5302)       ! defaut la plus grande valeur de la variable MP_EAGER_LIMIT, soit
( 5303)       ! 65 536 octets
( 5304)       IF (taille_tampon_alphanum == ' ') THEN






PGF90 (Version     10.5)          09/23/2010  09:48:58      page 188

( 5305)          taille_tampon = 65536
( 5306)       ELSE
( 5307)          READ(taille_tampon_alphanum,'(i9)') taille_tampon
( 5308)       END IF
( 5309) 
( 5310)       ! On est limite en mode d'adressage 32 bits a  1750 Mo pour la zone
( 5311)       ! "data" soit 7 segments, c.-a -d. 1750/8 = 210 Mo
( 5312)       IF (taille_tampon > 210000000) THEN
( 5313)          CALL ctl_stop( ' lib_mpp: Attention la valeur BUFFER_LENGTH doit etre <= 210000000' )
( 5314)          CALL mpi_abort(MPI_COMM_WORLD,2,code)
( 5315)       END IF
( 5316) 
( 5317)       CALL mpi_comm_rank(MPI_COMM_OPA,rang,code)
( 5318)       IF (rang == 0 ) PRINT *,'Taille du buffer alloue : ',taille_tampon
( 5319) 
( 5320)       ! Allocation du tampon et attachement
( 5321)       ALLOCATE(tampon(taille_tampon))
( 5322)       CALL mpi_buffer_attach(tampon,taille_tampon,code)
( 5323) 
( 5324)    END SUBROUTINE mpi_init_opa
( 5325) 
( 5326) !!DB 2008.04.24
( 5327) !!Routine copied from  /export/home/wangz/NEMO_FOR_DAVE/.../OPA_SRC/IOM/in_out_manager.F90 
( 5328)    SUBROUTINE ctl_stop( cd1, cd2, cd3, cd4, cd5,   &
( 5329)       &                 cd6, cd7, cd8, cd9, cd10 )
( 5330)       !!-----------------------------------------------------------------------
( 5331)       !!                  ***  ROUTINE  stop_opa  ***
( 5332)       !!
( 5333)       !! ** Purpose : ??? blah blah....
( 5334)       !!-----------------------------------------------------------------------
( 5335)       CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
( 5336)       CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
( 5337)       !!-----------------------------------------------------------------------
( 5338)       !
( 5339)       nstop = nstop + 1 
( 5340)       IF(lwp) THEN
( 5341)          WRITE(numout,"(/,' ===>>> : E R R O R',     /,'         ===========',/)") 
( 5342)          IF( PRESENT(cd1 ) ) WRITE(numout,*) cd1
( 5343)          IF( PRESENT(cd2 ) ) WRITE(numout,*) cd2
( 5344)          IF( PRESENT(cd3 ) ) WRITE(numout,*) cd3
( 5345)          IF( PRESENT(cd4 ) ) WRITE(numout,*) cd4
( 5346)          IF( PRESENT(cd5 ) ) WRITE(numout,*) cd5
( 5347)          IF( PRESENT(cd6 ) ) WRITE(numout,*) cd6
( 5348)          IF( PRESENT(cd7 ) ) WRITE(numout,*) cd7
( 5349)          IF( PRESENT(cd8 ) ) WRITE(numout,*) cd8
( 5350)          IF( PRESENT(cd9 ) ) WRITE(numout,*) cd9
( 5351)          IF( PRESENT(cd10) ) WRITE(numout,*) cd10
( 5352)       ENDIF
( 5353)       CALL FLUSH(numout)
( 5354)       !
( 5355)    END SUBROUTINE ctl_stop
( 5356) 
( 5357) 
( 5358) 
( 5359) 
( 5360) # 5569
( 5569)    !!----------------------------------------------------------------------
( 5570) END MODULE lib_mpp



