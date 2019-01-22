!!DB -- 2009.09.04 -- key_diadimg eliminated
MODULE limrst
   !!======================================================================
   !!                     ***  MODULE  limrst  ***
   !! Ice restart :  write the ice restart file
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_rst_write   : write of the restart file 
   !!   lim_rst_read    : read  the restart file 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE in_out_manager
   USE ice
   USE ioipsl
   USE dom_oce
   USE ice_oce         ! ice variables
   USE daymod

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC lim_rst_write  ! routine called by lim_step.F90
   PUBLIC lim_rst_read   ! routine called by lim_init.F90

   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limrst.F90,v 1.8 2006/03/10 10:35:43 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !!   Default option                                          NetCDF file
   !!----------------------------------------------------------------------

   SUBROUTINE lim_rst_write( niter )
      !!----------------------------------------------------------------------
      !!                    ***  lim_rst_write  ***
      !!
      !! ** purpose  :   output of sea-ice variable in a netcdf file
      !!
      !!----------------------------------------------------------------------
      ! Arguments
      INTEGER  ::    niter        ! number of iteration

      !- dummy variables :
      LOGICAL :: &
         llbon
      INTEGER :: &
         ji, jj
      INTEGER :: &
         inumwrs, it0, itime
      REAL(wp), DIMENSION(1) :: &
         zdept
      REAL(wp), DIMENSION(2) :: &
         zinfo
      REAL(wp),DIMENSION(jpi,jpj,35) :: &
         zmoment
      REAL(wp) :: &
         zsec, zdate0, zdt

      CHARACTER(len=45)  ::  ccfile

!sujie add-----------------------------------
      CHARACTER (len=8) ::   cln
     ! Name of the new restart file
       WRITE(cln,'(i4.4,i2.2,i2.2)') nyear, nmonth,nday
       ccfile = 'restart_ice_'//cln//'_out.nc'
!      ccfile = 'restart_ice_out.nc'
!--------------------------------------------

#if defined key_agrif
      if ( .NOT. Agrif_Root() ) then
         ccfile= TRIM(Agrif_CFixed())//'_'//TRIM(ccfile)
      endif
#endif

      inumwrs  = 61
      INQUIRE ( FILE = ccfile, EXIST = llbon )
      IF( llbon ) THEN
         OPEN ( UNIT = inumwrs , FILE = ccfile, STATUS = 'old' )
         CLOSE( inumwrs , STATUS = 'delete' )
      ENDIF


      it0      = niter
      zinfo(1) = FLOAT( nfice  )  ! coupling frequency OPA ICELLN  nfice
      zinfo(2) = FLOAT( it0   )   ! iteration number

      zsec     = 0.e0
      itime    = 0
      zdept(1) = 0.e0
      zdt      = rdt_ice * nstock

      ! Write in inumwrs

      DO jj = 1, jpj              ! 3D array: 10 time faster than 35 restput
         DO ji = 1, jpi
            zmoment(ji,jj,1)  = sxice(ji,jj)
            zmoment(ji,jj,2)  = syice(ji,jj)
            zmoment(ji,jj,3)  = sxxice(ji,jj)
            zmoment(ji,jj,4)  = syyice(ji,jj)
            zmoment(ji,jj,5)  = sxyice(ji,jj)
            zmoment(ji,jj,6)  = sxsn(ji,jj)
            zmoment(ji,jj,7)  = sysn(ji,jj)
            zmoment(ji,jj,8)  = sxxsn(ji,jj)
            zmoment(ji,jj,9)  = syysn(ji,jj)
            zmoment(ji,jj,10) = sxysn(ji,jj)
            zmoment(ji,jj,11) = sxa(ji,jj)
            zmoment(ji,jj,12) = sya(ji,jj)
            zmoment(ji,jj,13) = sxxa(ji,jj)
            zmoment(ji,jj,14) = syya(ji,jj)
            zmoment(ji,jj,15) = sxya(ji,jj)
            zmoment(ji,jj,16) = sxc0(ji,jj)
            zmoment(ji,jj,17) = syc0(ji,jj)
            zmoment(ji,jj,18) = sxxc0(ji,jj)
            zmoment(ji,jj,19) = syyc0(ji,jj)
            zmoment(ji,jj,20) = sxyc0(ji,jj)
            zmoment(ji,jj,21) = sxc1(ji,jj)
            zmoment(ji,jj,22) = syc1(ji,jj)
            zmoment(ji,jj,23) = sxxc1(ji,jj)
            zmoment(ji,jj,24) = syyc1(ji,jj)
            zmoment(ji,jj,25) = sxyc1(ji,jj)
            zmoment(ji,jj,26) = sxc2(ji,jj)
            zmoment(ji,jj,27) = syc2(ji,jj)
            zmoment(ji,jj,28) = sxxc2(ji,jj)
            zmoment(ji,jj,29) = syyc2(ji,jj)
            zmoment(ji,jj,30) = sxyc2(ji,jj)
            zmoment(ji,jj,31) = sxst(ji,jj)
            zmoment(ji,jj,32) = syst(ji,jj)
            zmoment(ji,jj,33) = sxxst(ji,jj)
            zmoment(ji,jj,34) = syyst(ji,jj)
            zmoment(ji,jj,35) = sxyst(ji,jj)
         END DO
      END DO

      CALL ymds2ju( nyear, nmonth, nday, zsec, zdate0 )
      CALL restini( 'NONE', jpi, jpj, glamt, gphit, 1 , zdept, ccfile, itime, zdate0, zdt, &
         &         inumwrs, domain_id=nidom )
      
      CALL restput( inumwrs, 'info'   ,   1,   1, 2 , 0, zinfo   )  ! restart informations
       
      CALL restput( inumwrs, 'hicif'  , jpi, jpj, 1 , 0, hicif   )  ! prognostic variables 
      CALL restput( inumwrs, 'hsnif'  , jpi, jpj, 1 , 0, hsnif   )
      CALL restput( inumwrs, 'frld'   , jpi, jpj, 1 , 0, frld    )
      CALL restput( inumwrs, 'sist'   , jpi, jpj, 1 , 0, sist    )
# if defined key_coupled
      CALL restput( inumwrs, 'albege' , jpi, jpj, 1 , 0, albege  )
# endif
      CALL restput( inumwrs, 'tbif'   , jpi, jpj, 3 , 0, tbif    )
      CALL restput( inumwrs, 'u_ice'  , jpi, jpj, 1 , 0, u_ice   )
      CALL restput( inumwrs, 'v_ice'  , jpi, jpj, 1 , 0, v_ice   )
      CALL restput( inumwrs, 'gtaux'  , jpi, jpj, 1 , 0, gtaux  )
      CALL restput( inumwrs, 'gtauy'  , jpi, jpj, 1 , 0, gtauy  )
      CALL restput( inumwrs, 'qstoif' , jpi, jpj, 1 , 0, qstoif  )
      CALL restput( inumwrs, 'fsbbq'  , jpi, jpj, 1 , 0, fsbbq   )
      CALL restput( inumwrs, 'moment' , jpi, jpj, 35, 0, zmoment )

      
      CALL restclo( inumwrs )

   END SUBROUTINE lim_rst_write


   SUBROUTINE lim_rst_read( niter )
      !-----------------------------------------------------------------------
      !  restart from a state defined in a binary file
      !-----------------------------------------------------------------------
      ! Arguments
      INTEGER  ::   niter        ! number of iteration

      !- dummy variables :
      CHARACTER(len=45)  ::  ccfile
      INTEGER :: &
        ji, jj
      INTEGER :: &
         inumrst, it0, it1, itime, ibvar, ifice
      LOGICAL :: &
         llog
      REAL(wp),DIMENSION(jpi,jpj) :: &
         zlamt, zphit
      REAL(wp),DIMENSION(jpi,jpj,35) :: &
         zmoment
      REAL(wp),DIMENSION(1) :: &
         zdept
      REAL(wp),DIMENSION(2) :: &
         zinfo
      REAL(wp) :: &
         zdate0, zdt
      CHARACTER ( len = 10 ) ::  &
         clvnames(60)       


       ccfile = 'restart_ice_in.nc'
#if defined key_agrif
      if ( .NOT. Agrif_Root() ) then
         ccfile= TRIM(Agrif_CFixed())//'_'//TRIM(ccfile)
      endif
#endif

      !Initialisations
      inumrst    = 71
      it0        = nit000
      itime      = 0
      llog       = .FALSE.
      zlamt(:,:) = 0.
      zphit(:,:) = 0.
      zdept(1)   = 0.

      CALL restini(ccfile , jpi, jpj, zlamt, zphit, 1 , zdept, 'NONE', itime, zdate0, zdt, inumrst, &
         &         domain_id=nidom )      
      CALL ioget_vname( inumrst, ibvar, clvnames )

      CALL restget    ( inumrst,'info', 1, 1 , 2, 0, llog, zinfo )
 
      ifice   = INT( zinfo(1) )
      it1     = INT( zinfo(2) )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_rst_read : READ restart file name ', ccfile,' at time step : ', it1
         WRITE(numout,*) '~~~~~~~~~~~~   number of variables   : ', ibvar
         WRITE(numout,*) '               NetCDF variables      : ', clvnames(1:ibvar)
      ENDIF

      
      !Control of date
      
      IF( ( it0 - it1 ) /= 1 .AND. ABS( nrstdt ) == 1 ) THEN
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) 'lim_rst_read ===>>>> : problem with nit000 for the restart'
            WRITE(numout,*) '   we stop. verify the file or rerun with the value  0 for the'
            WRITE(numout,*) '   control of time parameter  nrstdt'
            nstop = nstop + 1
         ENDIF
      ENDIF

      CALL restget( inumrst, 'hicif'  , jpi, jpj, 1 , 0, llog, hicif   )
      CALL restget( inumrst, 'hsnif'  , jpi, jpj, 1 , 0, llog, hsnif   )
      CALL restget( inumrst, 'frld'   , jpi, jpj, 1 , 0, llog, frld    )
      CALL restget( inumrst, 'sist'   , jpi, jpj, 1 , 0, llog, sist    )
# if defined key_coupled 
      CALL restget( inumrst, 'albege' , jpi, jpj, 1 , 0, llog, albege  )
# endif
      CALL restget( inumrst, 'tbif'   , jpi, jpj, 3 , 0, llog, tbif    )
      CALL restget( inumrst, 'u_ice'  , jpi, jpj, 1 , 0, llog, u_ice   )
      CALL restget( inumrst, 'v_ice'  , jpi, jpj, 1 , 0, llog, v_ice   )
      CALL restget( inumrst, 'gtaux'  , jpi, jpj, 1 , 0, llog, gtaux  )
      CALL restget( inumrst, 'gtauy'  , jpi, jpj, 1 , 0, llog, gtauy  )
      CALL restget( inumrst, 'qstoif' , jpi, jpj, 1 , 0, llog, qstoif  )
      CALL restget( inumrst, 'fsbbq'  , jpi, jpj, 1 , 0, llog, fsbbq   )
      CALL restget( inumrst, 'moment' , jpi, jpj, 35, 0, llog, zmoment )

      CALL restclo( inumrst )

      niter = it1
      DO jj = 1, jpj
         DO ji = 1, jpi
            sxice(ji,jj)  = zmoment(ji,jj,1)
            syice(ji,jj)  = zmoment(ji,jj,2)
            sxxice(ji,jj) = zmoment(ji,jj,3)
            syyice(ji,jj) = zmoment(ji,jj,4)
            sxyice(ji,jj) = zmoment(ji,jj,5)
            sxsn(ji,jj)   = zmoment(ji,jj,6)
            sysn(ji,jj)   = zmoment(ji,jj,7)
            sxxsn(ji,jj)  = zmoment(ji,jj,8)
            syysn(ji,jj)  = zmoment(ji,jj,9)
            sxysn(ji,jj)  = zmoment(ji,jj,10)
            sxa(ji,jj)    = zmoment(ji,jj,11)
            sya(ji,jj)    = zmoment(ji,jj,12)
            sxxa(ji,jj)   = zmoment(ji,jj,13)
            syya(ji,jj)   = zmoment(ji,jj,14)
            sxya(ji,jj)   = zmoment(ji,jj,15)
            sxc0(ji,jj)   = zmoment(ji,jj,16)
            syc0(ji,jj)   = zmoment(ji,jj,17)
            sxxc0(ji,jj)  = zmoment(ji,jj,18)
            syyc0(ji,jj)  = zmoment(ji,jj,19)
            sxyc0(ji,jj)  = zmoment(ji,jj,20)
            sxc1(ji,jj)   = zmoment(ji,jj,21)
            syc1(ji,jj)   = zmoment(ji,jj,22)
            sxxc1(ji,jj)  = zmoment(ji,jj,23)
            syyc1(ji,jj)  = zmoment(ji,jj,24)
            sxyc1(ji,jj)  = zmoment(ji,jj,25)
            sxc2(ji,jj)   = zmoment(ji,jj,26)
            syc2(ji,jj)   = zmoment(ji,jj,27)
            sxxc2(ji,jj)  = zmoment(ji,jj,28)
            syyc2(ji,jj)  = zmoment(ji,jj,29)
            sxyc2(ji,jj)  = zmoment(ji,jj,30)
            sxst(ji,jj)   = zmoment(ji,jj,31)
            syst(ji,jj)   = zmoment(ji,jj,32)
            sxxst(ji,jj)  = zmoment(ji,jj,33)
            syyst(ji,jj)  = zmoment(ji,jj,34)
            sxyst(ji,jj)  = zmoment(ji,jj,35)
         END DO
      END DO

      
   END SUBROUTINE lim_rst_read


#else
   !!----------------------------------------------------------------------
   !!   Default option :       Empty module            NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_rst_read             ! Empty routine
   END SUBROUTINE lim_rst_read
   SUBROUTINE lim_rst_write            ! Empty routine
   END SUBROUTINE lim_rst_write
#endif

   !!======================================================================
END MODULE limrst
