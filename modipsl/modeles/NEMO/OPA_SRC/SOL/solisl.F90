MODULE solisl
   !!==============================================================================
   !!                       ***  MODULE  solisl  ***
   !! Ocean island : specific treatment of island in rigid-lid case
   !!==============================================================================
#if defined key_islands
   !!----------------------------------------------------------------------
   !!   'key_islands' :                           islands in rigid-lid case
   !!----------------------------------------------------------------------
   !!   isl_dom     : locate islands in the domain
   !!   isl_pri     : control print of island gridpoint position
   !!   isl_pth     : Compute coeff. associated with path round each island
   !!   isl_mat     : Compute the matrix associated with islands
   !!   isl_bsf     : Compute the barotropic streamfunction of each island
   !!   isl_dyn_spg : Update the barotropic streamfunction trend with the
   !!                 Island contribution (call by dyn_spg)
   !!   isl_stp_ctl : print island information (call by stp_ctl)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE in_out_manager  ! I/O manager
   USE sol_oce         ! ocean solver 
   USE obc_oce         ! ocean open boundary condition
   USE lib_mpp         ! distributed memory computing

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC isl_dom        ! routine called by solver_init
   PUBLIC isl_mat        ! routine called by solver_init
   PUBLIC isl_bsf        ! routine called by solver_init
   PUBLIC isl_dyn_spg    ! routine called by dyn_spg
   PUBLIC isl_stp_ctl    ! routine called by stp_ctl

   !! * Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_isl = .TRUE.    !: 'key_islands' flag

   !! * module variable
   INTEGER :: numisl = 11      ! logical unit for island file only used
   !                           ! here during the initialization phase
   INTEGER ::   &
      nimlu, njmlu, nkmlu,  &  ! i-j-k-dimensions read
      nlmlu, nmmlu, nnmlu      ! read islands number
   INTEGER, DIMENSION(jpnisl,0:4,jpisl) ::   &
      miisl, mjisl             ! position of island grid-points
   INTEGER, DIMENSION(0:4,jpisl) ::   &
      mnisl                    ! number of grid-points for each island
 
   REAL(wp) ::    &
      replu                    ! read absolute precision
   REAL(wp), DIMENSION(jpisl,jpisl) ::   &
      aisl, aislm1             ! island matrix and its inverse
   REAL(wp), DIMENSION(jpi,jpj,jpisl) ::   &
      bsfisl                   ! barotropic streamfunction of island
   REAL(wp), DIMENSION(jpi,jpj,2,jpisl) ::   &
      acisl1, acisl2           ! coef. to compute circulations round islands
   REAL(wp), DIMENSION(jpisl) ::   &
      bisl,                 &  ! second member of island linear system
      visl                     ! trend of island stream function
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solisl.F90,v 1.6 2005/12/12 14:18:06 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE isl_dom
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE isl_dom  ***
      !!                    
      !! ** Purpose :   Locate island grid-points from mbathy array
      !!
      !! ** Method  :   The coordinates of ocean grid-points round an island 
      !!      are found from mbathy array read or computed in dommba.
      !!      we first compute zwb, an ocean/land mask defined as follows:
      !!              zwb(i,j)  =  0. over the main land and the ocean
      !!                        = -n  over the nth island
      !!      With the proper boundary conditions (defined by nperio)
      !!        Note: IF i,j are the coordinates of an ocean grid-point west
      !!      or east of a island, the corresponding coordinates miisl(ip,4,n)
      !!      or mjisl(ip,2,n) are those of the western or southern side of
      !!      the island (i.e. i+1 ou j+1, respectively)
      !!
      !! ** Action  :   compute mnisl, miisl, mjisl arrays defined as:
      !!      mnisl(i,n) : nb of grid-points along (i=0), north (i=1), north
      !!      (i=1), south (i=2), east (i=3), or west (i=4) of the nth island
      !!      miisl(ip,i,n), mjisl(ip,i,n) : (i,j) index of ipth u- or v-point
      !!      along (i=0), north (i=1),south (i=2), east (i=3), or west (i=4)
      !!      of the nth island
      !!
      !! History :
      !!    4.0  !  88-03  (G. Madec)  Original code
      !!    6.0  !  96-01  (G. Madec)  Suppress common workspace, use of bmask
      !!    8.5  !  96-01  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj, jn, jnil   ! dummy loop indices
      INTEGER ::   ind, inilt, ip, ipn, ips, ipe, ipw, ii, ij, iju
      INTEGER ::   iista, iiend, ijsta, ijend, ijstm1, ijenm1
      INTEGER ::   isrchne
      INTEGER, DIMENSION(jpj) ::   indil
      INTEGER, DIMENSION(jpi,jpj) ::   idil

      REAL(wp) znil
      REAL(wp), DIMENSION(jpi,jpj) ::   zwb
      !!----------------------------------------------------------------------


      ! 0. Islands index computed from mbathy
      ! -------------------------------------
      ! zwb=0 over the continent and ocean, =-n over the nth island
      
      ! Computation
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
            zwb(ji,jj) = MIN( 0 , mbathy(ji,jj  ), mbathy(ji+1,jj  ),   &
                                  mbathy(ji,jj+1), mbathy(ji+1,jj+1)  )
         END DO
      END DO
      
      ! Lateral boundary conditions on zwb

      !  mono- or macro-tasking environnement:
      IF( nperio == 2 ) THEN
         zwb(:, 1 ) = zwb(:, 2 )
      ELSEIF( nperio == 3 .OR. nperio == 4) THEN
! om 
!$$$          DO ji = 1, jpim1
! On ne peut pas partir de ji=1, car alors jiu=jpi, et cette valeur
! n'est pas encore initialiee. 
! Le cas ji=1 est de toute facon traite a partir de la ligne 135
         DO ji = 2, jpim1
            iju = jpi-ji+1
            zwb(ji,jpj  ) = zwb(iju,jpj-3)
            zwb(ji,jpjm1) = zwb(iju,jpj-2)
         END DO
      ELSEIF( nperio == 5 .OR. nperio == 6) THEN
         DO ji = 1, jpim1
            iju = jpi-ji
            zwb(ji,jpj  ) = zwb(iju,jpj-2)
         END DO
         DO ji = jpi/2+1, jpi-1
            iju = jpi-ji
            zwb(ji,jpjm1) = zwb(iju,jpjm1)
         END DO
      ELSE
         zwb(:,jpj) = 0.e0
      ENDIF

      IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
         zwb( 1 ,:) = zwb(jpim1,:)
         zwb(jpi,:) = zwb(  2  ,:)
      ELSE
         zwb( 1 ,:) = 0.e0
         zwb(jpi,:) = 0.e0
      ENDIF
      IF( lk_mpp )   CALL lbc_lnk( zwb, 'G', 1. )


      ! 1. Initialization for the search of island grid-points
      ! ------------------------------------------------------

      IF( lk_mpp ) THEN

         ! Mpp : The overlap region are not taken into account
         ! (islands bondaries are searched over subdomain only)
         iista =  1   + jpreci
         iiend = nlci - jpreci
         ijsta =  1   + jprecj
         ijend = nlcj - jprecj
         ijstm1=  1   + jprecj
         ijenm1= nlcj - jprecj
         IF( nbondi == -1 .OR. nbondi == 2 ) THEN
            iista  = 1
         ENDIF
         IF( nbondi ==  1 .OR. nbondi == 2 ) THEN
            iiend  = nlci
         ENDIF
         IF( nbondj == -1 .OR. nbondj == 2 ) THEN
            ijsta  = 1
            ijstm1 = 2
         ENDIF
         IF( nbondj ==  1 .OR. nbondj == 2 ) THEN
            ijend  = nlcj
            ijenm1 = nlcj-1
         ENDIF
         IF( npolj == 3 .OR. npolj == 4 ) THEN
            ijend  = nlcj-2
            ijenm1 = nlcj-2
         ENDIF 
      ELSE
         ! mono- or macro-tasking environnement: full domain scan
         iista  = 1
         iiend  = jpi
         ijsta  = 1
         ijstm1 = 2
         IF( nperio == 3 .OR. nperio == 4 ) THEN
            ijend  = jpj-2
            ijenm1 = jpj-2
         ELSEIF( nperio == 5 .OR. nperio == 6 ) THEN
            ijend  = jpj-1
            ijenm1 = jpj-1
         ELSE
            ijend  = jpj
            ijenm1 = jpj-1
         ENDIF
      ENDIF


      ! 2. Loop over island
      ! -------------------

      DO jnil = 1, jpisl

         ! 2.1 Initialization to zero of miisl, mjisl of the jnil island

         miisl(:,:,jnil) = 0
         mjisl(:,:,jnil) = 0
         
         ! 2.2 Search grid-points of island jnil
 
         indil(:)   = 0
         idil (:,:) = 0
         
         znil = - FLOAT( jnil )
         DO jj = ijsta, ijend
            indil(jj) = 0
            ind = 0
            DO ji = iista, iiend
               IF( zwb(ji,jj) == znil ) THEN
                  ind = ind + 1
                  idil(ind,jj) = ji
                  indil(jj) = indil(jj) + 1
               ENDIF
            END DO
         END DO
         
         ! 2.3 Check the number of island
         
         inilt = 0
         DO jj = ijsta, ijend
            inilt = inilt + indil(jj)
         END DO
         IF( lk_mpp )   CALL mpp_sum( inilt )   ! sum over the global domain

         IF( inilt == 0 ) THEN
            IF(lwp) THEN
               WRITE(numout,*) ' isldom: there is not island number: ', jnil,' while jpisl= ', jpisl
               WRITE(numout,*) ' change parameter.h'
            ENDIF
            STOP 'isldom'      !cr replace by nstop
         ENDIF
         
         ! 2.4 Coastal island grid-points (miisl,mjisl)
         
         ip  = 0
         ipn = 0
         ips = 0
         ipe = 0
         ipw = 0

         ! South line (ij=1)
         ij = 1
         DO jn = 1, indil(ij)
            ii = idil(jn,ij)
            IF( (ij+njmpp-1) == 1 .AND. ii > jpreci .AND. ii < (nlci-jpreci+1) ) THEN
               IF(  zwb(ii-1,ij) * zwb(ii+1,ij) * zwb(ii,ij+1)  == 0. ) THEN
                  ip = ip + 1
                  miisl(ip,0,jnil) = ii
                  mjisl(ip,0,jnil) = ij
                  IF( zwb(ii-1,ij) == 0. ) THEN
                     ipw = ipw + 1
                     miisl(ipw,4,jnil) = ii
                     mjisl(ipw,4,jnil) = ij
                  ENDIF
                  IF( zwb(ii+1,ij) == 0. ) THEN
                     ipe = ipe+1
                     IF( (nperio == 1 .OR. nperio == 4.OR. nperio == 6)  &
                        .AND. ii == (nlci-jpreci) ) THEN
                        miisl(ipe,3,jnil) = 1+jpreci
                     ELSE
                        miisl(ipe,3,jnil) = ii + 1
                     ENDIF
                     mjisl(ipe,3,jnil) = ij
                  ENDIF
                  IF( zwb(ii,ij+1) == 0. ) THEN
                     ipn = ipn+1
                     miisl(ipn,1,jnil) = ii
                     mjisl(ipn,1,jnil) = ij + 1
                  ENDIF
               ENDIF
            ENDIF
         END DO
         
         ! Middle lines (2=<jj=<jpjm1 or jpj-2 if north fold b.c.)
         DO jj = ijstm1, ijenm1
            DO jn = 1, indil(jj)
               ii = idil(jn,jj)
               IF( ii > jpreci .AND. ii < (nlci-jpreci+1) ) THEN
                  IF( (zwb(ii-1, jj )*zwb(ii+1, jj )*   &
                       zwb( ii ,jj-1)*zwb( ii ,jj+1)) == 0. ) THEN
                     ip = ip + 1
                     miisl(ip,0,jnil) = ii
                     mjisl(ip,0,jnil) = jj
                     IF( zwb(ii-1,jj) == 0. ) THEN
                        ipw = ipw + 1
                        miisl(ipw,4,jnil) = ii
                        mjisl(ipw,4,jnil) = jj
                     ENDIF
                     IF( zwb(ii+1,jj) == 0. ) THEN
                        ipe = ipe + 1
                        IF((nperio == 1.OR.nperio == 4.OR.nperio == 6)   &
                           .AND.ii == (nlci-jpreci) ) THEN
                           miisl(ipe,3,jnil) = 1 + jpreci
                        ELSE
                           miisl(ipe,3,jnil) = ii + 1
                        ENDIF
                        mjisl(ipe,3,jnil) = jj
                     ENDIF
                     IF( zwb(ii,jj-1) == 0. ) THEN
                        ips = ips + 1
                        miisl(ips,2,jnil) = ii
                        mjisl(ips,2,jnil) = jj
                     ENDIF
                     IF( zwb(ii,jj+1) == 0. ) THEN
                        ipn = ipn + 1
                        miisl(ipn,1,jnil) = ii
                        mjisl(ipn,1,jnil) = jj + 1
                     ENDIF
                  ENDIF
               ENDIF
            END DO
         END DO
         
         ! North line (jj=jpj) only if not north fold b.c.
         IF( nperio /= 3 .AND. nperio /= 4 .AND. nperio /= 5 .AND. nperio /= 6 ) THEN
            ij = jpj
            DO jn = 1, indil(ij)
               ii = idil(jn,ij)
               IF( (ij+njmpp-1) == jpjglo .AND. ii > jpreci .AND. ii < (nlci-jpreci+1) ) THEN
                  IF( (zwb(ii-1, ij )*zwb(ii+1, ij )* zwb( ii ,ij-1) ) == 0. ) THEN
                     ip = ip+1
                     miisl(ip,0,jnil) = ii
                     mjisl(ip,0,jnil) = ij
                     IF( zwb(ii-1,ij) == 0. ) THEN
                        ipw = ipw+1
                        miisl(ipw,4,jnil) = ii
                        mjisl(ipw,4,jnil) = ij
                     ENDIF
                     IF( zwb(ii+1,ij) == 0. ) THEN
                        ipe = ipe+1
                        IF( (nperio == 1) .AND. ii == (nlci-jpreci) ) THEN
                           miisl(ipe,3,jnil) = 1+jpreci
                        ELSE
                           miisl(ipe,3,jnil) = ii+1
                        ENDIF
                        mjisl(ipe,3,jnil) = ij
                     ENDIF
                     IF( zwb(ii,ij-1) == 0. ) THEN
                        ips = ips+1
                        miisl(ips,2,jnil) = ii
                        mjisl(ips,2,jnil) = ij
                     ENDIF
                  ENDIF
               ENDIF
            END DO
         ENDIF
         
         mnisl(0,jnil) = ip
         mnisl(1,jnil) = ipn
         mnisl(2,jnil) = ips
         mnisl(3,jnil) = ipe
         mnisl(4,jnil) = ipw
         
         ! Take account of redundant points
         
         IF( lk_mpp )   CALL mpp_sum( ip )   ! sum over the global domain
         
         IF( ip > jpnisl ) THEN
            IF(lwp) THEN
               WRITE(numout,*) ' isldom: the island ',jnil,' has ',   &
                  mnisl(0,jnil),' grid-points, while jpnisl= ', jpnisl,ip
               WRITE(numout,*) ' change parameter.h'
            ENDIF
            STOP 'isldom'    !cr => nstop
         ENDIF
         
         ! 2.5 Set to zero the grid-points of the jnil island in x
         
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zwb(ji,jj)+FLOAT(jnil) == 0. ) zwb(ji,jj) = 0. 
            END DO
         END DO
         
      END DO
      

      ! 3. Check the number of island
      ! -----------------------------

      inilt = isrchne( jpij, zwb(1,1), 1, 0. )
      IF( lk_mpp )   CALL mpp_min( inilt )   ! min over the global domain

      IF( inilt /= jpij+1 ) THEN
         IF(lwp) THEN
            WRITE(numout,*) ' isldom: there is at least one more ',   &
                  'island in the domain and jpisl=', jpisl
            WRITE(numout,*) ' change parameter.h'
         ENDIF
         STOP 'isldom'
      ENDIF


      ! 4. Print of island parametres and arrays
      ! ----------------------------------------

      CALL isl_pri


      ! 5. Array for computation of circulation arround islands
      ! -------------------------------------------------------

      CALL isl_pth

   END SUBROUTINE isl_dom


   SUBROUTINE isl_pri
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE isl_pri  ***
      !!              
      !! ** Purpose :   Print islands variables and islands arrays
      !!
      !! ** Method  :
      !!
      !! History :
      !!   1.0  !  88-03  (G. Madec)  Original code
      !!   8.5  !  02-08  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jni, jnp      ! dummy loop variables
      INTEGER ::   &
         ip, ipn, ips, ipe, ipw      ! temporary integers
      !!----------------------------------------------------------------------
      
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '*** islpri number of islands : jpisl =',jpisl
      ENDIF

      DO jni = 1, jpisl
         ip  = mnisl(0,jni)
         ipn = mnisl(1,jni)
         ips = mnisl(2,jni)
         ipe = mnisl(3,jni)
         ipw = mnisl(4,jni)
         IF( lk_mpp ) THEN
            CALL mpp_sum( ip  )   ! sums over the global domain
            CALL mpp_sum( ipn )
            CALL mpp_sum( ips )
            CALL mpp_sum( ipe )
            CALL mpp_sum( ipw )
         ENDIF
         IF(lwp) THEN
            WRITE(numout,9000) jni
            WRITE(numout,9010) ip, ipn, ips, ipe, ipw
            WRITE(numout,9020)
            DO jnp = 1, mnisl(0,jni)
               WRITE(numout,9030) jnp, ( miisl(jnp,ji,jni)+nimpp-1,   &
                                         mjisl(jnp,ji,jni)+njmpp-1, ji = 0, 4 )
            END DO
         ENDIF
      END DO

      ! FORMAT   !!cr => no more format
 9000 FORMAT(/, /, 'island number= ', i2 )
 9010 FORMAT(/, 'npil=',i4,' npn=',i3,' nps=',i3,' npe=',i3,' npw=',i3 )
 9020 FORMAT(/,'     * isl point *  point n  *  point s  *  point e ', '*  point w  *')
 9030 FORMAT(i4,' * (',i3,',',i3,') * (',i3,',',i3,') * (',i3,',',i3,   &
          ') * (',i3,',',i3,') * (',i3,',',i3,') *' )

   END SUBROUTINE isl_pri


   SUBROUTINE isl_pth
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE isl_pth  ***
      !!             
      !! ** Purpose :   intialize arrays for the computation of streamfunction
      !!      around islands
      !!
      !! ** Method  :
      !!
      !! ** Action  : - acisl1(i,j,ii,n): coefficient n-s (ii=1) and e-w (ii=2)
      !!                for the calculation of mu and mv around the island n      
      !!              - acisl2(i,j,ii,n): coefficient n-s (ii=1) and e-w (ii=2)
      !!                for the calculation of bsfd around the island n
      !!
      !! History :
      !!   1.0  !  88-03  (G. Madec)  Original code
      !!   8.5  !  02-08  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   jni, jii, jnp    ! dummy loop indices
      INTEGER ::   ii, ij           ! temporary integers
      !!----------------------------------------------------------------------
      
      ! 1. Initialisation
      ! -----------------
      acisl1(:,:,:,:) = 0.e0
      acisl2(:,:,:,:) = 0.e0
      bsfisl(:,:,:)   = 0.e0

      ! 2. Coefficient arrays
      ! ---------------------
      DO jni = 1, jpisl
         DO jii = 1, 4
            DO jnp = 1,mnisl(jii,jni)
               ii = miisl(jnp,jii,jni)
               ij = mjisl(jnp,jii,jni)
               IF( jii <= 2 ) THEN
                  ! north and south points
                  acisl1(ii,ij,1,jni) = float( 2*jii-3) * e1u(ii,ij)
                  acisl2(ii,ij,1,jni) = float( 2*jii-3) * e1u(ii,ij) * hur(ii,ij) / e2u(ii,ij)
               ELSE
                  ! east and west points
                  acisl1(ii,ij,2,jni) = float(-2*jii+7) * e2v(ii,ij)
                  acisl2(ii,ij,2,jni) = float(-2*jii+7) * e2v(ii,ij) * hvr(ii,ij) / e1v(ii,ij)
               ENDIF
            END DO
         END DO
      END DO
      
   END SUBROUTINE isl_pth

   !!----------------------------------------------------------------------
   !!   Default option :                                        NetCDF file
   !!----------------------------------------------------------------------

   SUBROUTINE isl_mat
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE isl_mat  ***
      !!                
      !! ** Purpose :   Compute and invert the island matrix 
      !!
      !! ** Method  :   aisl(jni,jnj) matrix constituted by bsfisl circulation
      !!
      !! ** Action  : - aisl     : island matrix
      !!              - aislm1   : invert of the island matrix
      !!      file
      !!              - islands  : contain bsfisl, aisl and aislm1
      !!
      !! History :
      !!   1.0  !  88-10  (G. Madec)  Original code
      !!   7.0  !  96-01  (G. Madec)  suppression of common work arrays
      !!   8.5  !  02-08  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl
      
      !! * Local declarations
      INTEGER ::   ji, jj, jni, jnj, jn, jl   ! dummy loop indices
      INTEGER ::   itime, ibvar, ios          ! temporary integers
      LOGICAL ::   llog
      CHARACTER (len=32) ::   clname
      CHARACTER (len=8 ) ::   clvnames(100)
      REAL(wp), DIMENSION(1) ::   zdept
      REAL(wp), DIMENSION(jpi,jpj) ::   zlamt, zphit
      REAL(wp), DIMENSION(jpi,jpj,2) ::   zwx
      REAL(wp), DIMENSION(jpisl*jpisl) ::   ztab
      !!----------------------------------------------------------------------


      ! I. Island matrix lecture in numisl (if it exists)
      ! ==================================

      ! Lecture
      zlamt(:,:) = 0.
      zphit(:,:) = 0.
      zdept(1)   = 0.
      itime = 0
      clvnames="        "
      clname = 'islands'
      CALL ioget_vname(numisl, ibvar, clvnames)
      IF(lwp) WRITE(numout,*) clvnames
      ios=0
      DO jn=1,100
        IF(clvnames(jn) == 'aisl') ios=1
      END DO
      IF( ios == 0 ) go to 110 

      CALL restget( numisl, 'aisl'  , jpisl, jpisl, 1, 0, llog, aisl   )
      CALL restget( numisl, 'aislm1', jpisl, jpisl, 1, 0, llog, aislm1 )
      CALL restclo( numisl )
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*)' islmat: lecture aisl/aislm1 in numisl done'
         WRITE(numout,*)' ~~~~~~'
         WRITE(numout,*)
         WRITE(numout,*) '        island matrix : '
         WRITE(numout,*)
         
         DO jnj = 1, jpisl
            WRITE(numout,'(8e12.4)') ( aisl(jni,jnj), jni = 1, jpisl )
         END DO

         WRITE(numout,*)
         WRITE(numout,*) '       inverse of the island matrix'
         WRITE(numout,*)

         DO jnj = 1, jpisl
            WRITE(numout,'(12e11.3)') ( aislm1(jni,jnj), jni=1,jpisl )
         END DO
      ENDIF
      
      RETURN

 110  CONTINUE


      ! II. Island matrix computation
      ! =============================

      DO jnj = 1, jpisl

         ! Circulation of bsf(jnj) around island jni
         
         DO jj = 2, jpj
               zwx(:,jj,1) = -( bsfisl(:,jj,jnj) - bsfisl(:,jj-1,jnj) )
         END DO
         zwx(:,1,1) = 0.e0
         
         DO jj = 1, jpj
            DO ji = 2, jpi
               zwx(ji,jj,2) = ( bsfisl(ji,jj,jnj) - bsfisl(ji-1,jj,jnj) )
            END DO
         END DO
         zwx(1,:,2) = 0.e0
         
         ! Island matrix
         
         DO jni = 1, jpisl
            aisl(jni,jnj) = 0.e0
            DO jl = 1, 2
               DO jj=1,jpj
                  DO ji=1,jpi
                     aisl(jni,jnj) = aisl(jni,jnj) + acisl2(ji,jj,jl,jni)*zwx(ji,jj,jl)
                  END DO
               END DO
            END DO
         END DO
         
      END DO
      IF( lk_mpp ) THEN
         DO jnj = 1, jpisl
            DO jni = 1, jpisl
               ztab(jni+(jnj-1)*jpisl) = aisl(jni,jnj)
            END DO
         END DO
         CALL mpp_sum( ztab, jpisl*jpisl )   ! sum over the global domain
      ENDIF

      ! 1.3 Control print

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'islmat : island matrix'
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*)
         
         DO jnj = 1, jpisl
            WRITE(numout,'(8e12.4)') ( aisl(jni,jnj), jni = 1, jpisl )
         END DO
      ENDIF
      

      ! 2. Invertion of the island matrix
      ! ---------------------------------
      
      ! 2.1 Call of an imsl routine for the matrix invertion

      CALL linrg( jpisl, aisl, jpisl, aislm1, jpisl )

      ! 2.2 Control print
      
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'islmat : inverse of the island matrix'
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*)
         
         DO jnj = 1, jpisl
            WRITE(numout, '(12e11.3)')  '        ', ( aislm1(jni,jnj), jni=1, jpisl )
         END DO
      ENDIF
      

      ! 3. Output of aisl and aislm1 in numisl
      ! --------------------------------------

      CALL restput( numisl, 'aisl'  , jpisl, jpisl, 1, 0, aisl   )
      CALL restput( numisl, 'aislm1', jpisl, jpisl, 1, 0, aislm1 )
      CALL restclo( numisl )

   END SUBROUTINE isl_mat


   SUBROUTINE isl_bsf
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE isl_bsf  ***
      !!          
      !! ** Purpose :
      !!         Compute the barotropic stream function associated with each
      !!      island using FETI or preconditioned conjugate gradient method
      !!      or read them in numisl.
      !!
      !! ** Method :
      !!
      !! ** input/output file :
      !!            numisl            : barotropic stream function associated
      !!                          with each island of the domain
      !!
      !! ** Action :
      !!       bsfisl, the streamfunction which takes the value 1 over island ni,
      !!    and 0 over the others islands
      !!       file 'numisl' barotropic stream function associated
      !!                              with each island of the domain
      !!
      !! History :
      !!        !  87-10  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  93-03  (G. Madec)  release 7.1
      !!        !  93-04  (M. Guyon)  loops and suppress pointers
      !!        !  96-11  (A. Weaver)  correction to preconditioning
      !!        !  98-02  (M. Guyon)  FETI method
      !!        !  99-11  (M. Imbard)  NetCDF FORMAT with IOIPSL
      !!   8.5  !  02-08  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl
      USE solpcg
      USE solfet
      USE solsor

      !! * Local declarations
      LOGICAL  ::   llog, llbon
      CHARACTER (len=10) ::  clisl
      CHARACTER (len=32) ::  clname, clname2
      INTEGER  ::   ji, jj, jni, jii, jnp, je   ! dummy loop indices
      INTEGER  ::   iimlu, ijmlu, inmlu, iju
      INTEGER  ::   ii, ij, icile, icut, inmax, indic
      INTEGER  ::   itime, ie
      REAL(wp) ::   zepsr, zeplu, zgwgt
      REAL(wp) ::   zep(jpisl), zlamt(jpi,jpj), zphit(jpi,jpj), zdept(1), zprec(4)
      REAL(wp) ::   zdate0, zdt
      REAL(wp) ::   t2p1(jpi,1,1)
      INTEGER  ::   iloc
      !!----------------------------------------------------------------------


      ! 0. Initializations
      ! ==================

      icile         = 0      ! set to zero the convergence indicator
      bsfisl(:,:,:) = 0.e0   ! set to zero of bsfisl


      ! I. Lecture of bsfisl in numisl (if it exists)
      ! =============================================

      icut  = 0
      iimlu = 0
      ijmlu = 0
      inmlu = 0
      zeplu = 0.
      zlamt(:,:) = 0.
      zphit(:,:) = 0.
      zdept(1)   = 0.
      itime = 0
      clname = 'islands'
      ie=1
      DO je = 1, 32
        IF( clname(je:je) /= ' ' ) ie = je
      END DO
      clname2 = clname(1:ie)//".nc"
      INQUIRE( FILE=clname2, EXIST=llbon )
! islands FILE does not EXIST : icut=999
      IF( llbon ) THEN 
         ! island FILE is present 
         CALL restini(clname,jpi,jpj,zlamt,zphit,1,zdept,  &
            &         'NONE',itime,zdate0,zdt,numisl,domain_id=nidom)
         CALL restget(numisl,'PRECISION',1,1,4,0,llog,zprec)
         iimlu = NINT( zprec(1) )
         ijmlu = NINT( zprec(2) )
         inmlu = NINT( zprec(3) )
         zeplu = zprec(4)
         ! the read domain does not correspond to the model one : icut=999
         IF( iimlu /= jpi .OR. ijmlu /= jpj .OR. inmlu /= jpisl ) THEN
            icut = 999
            CALL restclo(numisl)
         ELSE 
            DO jni = 1, jpisl
               IF( jni < 10 ) THEN
                  WRITE(clisl,'("island",I1)') jni
               ELSEIF( jni < 100 ) THEN
                  WRITE(clisl,'("island",I2)') jni
               ELSE 
                  WRITE(clisl,'("island",I3)') jni
               ENDIF
               CALL restget(numisl,clisl,jpi,jpj,1,0,llog, bsfisl(:,:,jni))
            END DO
         ENDIF
      ELSE 
         ! islands FILE does not EXIST : icut=999 
         icut = 999
         CALL restclo(numisl)
      ENDIF
      
      ! the read precision is not the required one : icut=888
      IF( zeplu > epsisl ) THEN 
         icut = 888
         CALL restclo(numisl)
      ENDIF

      ! Control print
      IF( icut == 999 ) THEN
          IF(lwp) THEN
              WRITE(numout,*)
              WRITE(numout,*) 'islbsf : lecture bsfisl in numisl failed'
              WRITE(numout,*) '~~~~~~'
              WRITE(numout,*) '         icut= ', icut
              WRITE(numout,*) '         imlu= ', iimlu, ' jmlu= ', ijmlu,   &
                                      ' nilu= ', inmlu, ' epsisl lu= ', zeplu
              WRITE(numout,*) '         the bsfisl are computed from zero'
          ENDIF
        ELSEIF( icut == 888 ) THEN
          IF(lwp) THEN
              WRITE(numout,*)
              WRITE(numout,*) 'islbsf : lecture bsfisl in numisl done'
              WRITE(numout,*) '~~~~~~'
              WRITE(numout,*) '         the required accuracy is not reached'
              WRITE(numout,*) '         epsisl lu= ', zeplu,' epsisl required ', epsisl
              WRITE(numout,*) '         the bsfisl are computed from the read values'
          ENDIF
        ELSE
          IF(lwp) THEN
              WRITE(numout,*)
              WRITE(numout,*) 'islbsf : lecture bsfisl in numisl done'
              WRITE(numout,*) '~~~~~~'
          ENDIF
          RETURN
      ENDIF
      
      
      ! II. Compute the bsfisl (if icut=888 or 999)
      ! ============================================

      ! save nmax
      inmax = nmax

      ! set the number of iteration of island computation
      nmax = nmisl

      ! Loop over islands
      ! -----------------
      
      DO jni = 1, jpisl


         ! 1. Initalizations of island computation
         ! ---------------------------------------
         
         ! Set the pcg solution gcb to zero
         gcb(:,:) = 0.e0

         ! Set first guess gcx either to zero or to the read bsfisl
         IF( icut == 999 ) THEN
            gcx(:,:) = 0.e0
         ELSEIF( icut == 888 ) THEN
            IF(lwp) WRITE(numout,*) ' islbsf: bsfisl read are used as first guess'
            ! c a u t i o n: bsfisl masked because it contains 1 along island seaside
            gcx(:,:) = bsfisl(:,:,jni) * bmask(:,:)
         ENDIF
         
         ! Right hand side of the streamfunction equation
         
         IF( lk_mpp ) THEN

            ! north fold treatment
            IF( npolj == 3 .OR. npolj == 5)   iloc=jpiglo-(nimpp-1+nimppt(nono+1)-1)
            IF( npolj == 4 .OR. npolj == 6)   iloc=jpiglo-2*(nimpp-1)
            t2p1(:,1,1) = 0.e0
            ! north and south grid-points
            DO jii = 1, 2
               DO jnp = 1, mnisl(jii,jni)
                  ii = miisl(jnp,jii,jni)
                  ij = mjisl(jnp,jii,jni)
                  IF( ( npolj == 3 .OR. npolj == 4 ) .AND.   &
                     ( ij == nlcj-1 .AND. jii == 1) ) THEN 
                     iju=iloc-ii+1
                     t2p1(iju,1,1) =  t2p1(iju,1,1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
                  ELSEIF( ( npolj == 5 .OR. npolj == 6 ) .AND.   &
                     ( ij == nlcj-1 .AND. jii == 1) ) THEN 
                     iju=iloc-ii
                     gcb(ii,ij) =  gcb(ii,ij) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
                     t2p1(iju,1,1) =  t2p1(iju,1,1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
                  ELSE  
                     gcb(ii,ij-jii+1) =  gcb(ii,ij-jii+1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
                  ENDIF
               END DO
            END DO
         
            ! east and west grid-points
         
            DO jii = 3, 4
               DO jnp = 1, mnisl(jii,jni)
                  ii = miisl(jnp,jii,jni)
                  ij = mjisl(jnp,jii,jni)
                  gcb(ii-jii+3,ij) = gcb(ii-jii+3,ij) + hvr(ii,ij) * e2v(ii,ij) / e1v(ii,ij)
               END DO
            END DO

            IF( lk_mpp )   CALL mpplnks( gcb )   !!bug ? should use an lbclnk ? is it possible?

         ELSE

            ! north and south grid-points
            DO jii = 1, 2
               DO jnp = 1, mnisl(jii,jni)
                  ii = miisl(jnp,jii,jni)
                  ij = mjisl(jnp,jii,jni)
                  IF( ( nperio == 3 .OR. nperio == 4 ) .AND.   &
                     ( ij == jpj-1 .AND. jii == 1) ) THEN 
                     gcb(jpi-ii+1,ij-1) = gcb(jpi-ii+1,ij-1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
                  ELSEIF( ( nperio == 5 .OR. nperio == 6 ) .AND.   &
                     ( ij == jpj-1 .AND. jii == 1) ) THEN 
                     gcb(ii,ij) =  gcb(ii,ij) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij)
                     gcb(jpi-ii,ij) = gcb(jpi-ii,ij) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
                  ELSE  
                     gcb(ii,ij-jii+1) =  gcb(ii,ij-jii+1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij)
                  ENDIF
               END DO
            END DO

            ! east and west grid-points
            DO jii = 3, 4
               DO jnp = 1, mnisl(jii,jni)
                  ii = miisl(jnp,jii,jni)
                  ij = mjisl(jnp,jii,jni)
                  IF( bmask(ii-jii+3,ij) /= 0. ) THEN
                     gcb(ii-jii+3,ij) = gcb(ii-jii+3,ij) + hvr(ii,ij) * e2v(ii,ij) / e1v(ii,ij)
                  ELSE
                     ! east-west cyclic boundary conditions
                     IF( ii-jii+3 == 1 ) THEN
                        gcb(jpim1,ij) = gcb(jpim1,ij) + hvr(ii,ij) * e2v(ii,ij) / e1v(ii,ij)
                     ENDIF
                  ENDIF
               END DO
            END DO
         ENDIF

         ! Preconditioned right hand side and absolute precision

         IF( nsolv == 3 ) THEN 
            ! FETI method
            ncut    = 0
            rnorme  = 0.e0
            gcb(:,:) = bmask(:,:) * gcb(:,:)
            DO jj = 1, jpj
               DO ji = 1, jpi
                  rnorme = rnorme + gcb(ji,jj) * gcb(ji,jj)
               END DO
            END DO
            
            IF( lk_mpp )   CALL mpp_sum( rnorme )

            IF(lwp) WRITE(numout,*) 'rnorme ', rnorme
            epsr  = epsisl * epsisl * rnorme
            indic = 0
         ELSE 
            ncut   = 0
            rnorme = 0.e0
            DO jj = 1, jpj
               DO ji = 1, jpi
                  gcb  (ji,jj) = gcdprc(ji,jj) * gcb(ji,jj)
                  zgwgt  = gcdmat(ji,jj) * gcb(ji,jj)
                  rnorme = rnorme + gcb(ji,jj) * zgwgt
               END DO
            END DO
            IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain

            IF(lwp) WRITE(numout,*) 'rnorme ', rnorme
            epsr = epsisl * epsisl * rnorme
            indic = 0
         ENDIF


         ! 3. PCG solver for gcp.gcx=gcb in monotask
         ! -------------------------------------------

         IF( nsolv == 3 ) THEN 
            epsilo = epsisl       ! precision to compute Islands matrix A
            CALL sol_fet( indic )  ! FETI method
            epsilo = eps          ! precision to compute grad PS
         ELSE 
            CALL sol_pcg( indic )  ! pcg method
         ENDIF


         ! 4. Save the solution in bsfisl
         ! ------------------------------

         bsfisl(:,:,jni) = gcx(:,:)


         ! 5. Boundary conditions
         ! ----------------------

         ! set to 1. coastal gridpoints of the island
         DO jnp = 1, mnisl(0,jni)
            ii = miisl(jnp,0,jni)
            ij = mjisl(jnp,0,jni)
            bsfisl(ii,ij,jni) = 1.
         END DO

         ! cyclic boundary conditions
         IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
            bsfisl( 1 ,:,jni) = bsfisl(jpim1,:,jni)
            bsfisl(jpi,:,jni) = bsfisl(  2  ,:,jni)
         ENDIF
         IF( nperio == 3 .OR. nperio == 4 ) THEN
            DO ji = 1, jpim1
               iju = jpi-ji+1
               bsfisl(ji,jpj-1,jni) = bsfisl(iju,jpj-2,jni)
               bsfisl(ji, jpj ,jni) = bsfisl(iju,jpj-3,jni)
            END DO
         ENDIF
         IF( nperio == 5 .OR. nperio == 6 ) THEN
            DO ji = 1, jpi-1
               iju=jpi-ji
               bsfisl(ji,jpj,jni) = bsfisl(iju,jpj-2,jni)
            END DO
            DO ji = jpi/2+1, jpi-1
               iju=jpi-ji
               bsfisl (ji,jpjm1,jni) = bsfisl(iju,jpjm1,jni)
            END DO
         ENDIF
         IF( lk_mpp )   CALL lbc_lnk( bsfisl(:,:,jni), 'G', 1. )   ! link at G-point


         ! 6. Control print
         ! ----------------

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' islbsf: island number: ', jni
         IF(lwp) WRITE (numout,9290) niter, res, SQRT(epsr)/epsisl
         zep(jni) =  MAX(epsisl, res/(SQRT(epsr)/epsisl))
         IF( indic <  0 ) THEN
            icile = icile-1
            IF(lwp) WRITE(numout,*) '  pcg do not converge for island: ', jni
            IF(lwp) WRITE(numout,*) '      Precision reached: ',zep(jni)
         ENDIF

9290     FORMAT('        niter :',i4,' , res :',e20.10,' , gcb :',e20.10)

         !                                          !====================
      END DO                                        !  End Loop islands
      !                                             !====================

      ! 7. Reset PCG
      ! ------------

      ! reset the number of iteration for pcg
      nmax = inmax

      ! reset to zero pcg arrays
      gcx  (:,:) = 0.e0
      gcxb (:,:) = 0.e0
      gcb  (:,:) = 0.e0
      gcr  (:,:) = 0.e0
      gcdes(:,:) = 0.e0
      gccd (:,:) = 0.e0


      ! III. Output of bsfisl in numisl
      ! ===============================

      CALL ymds2ju( 0, 1, 1, 0.e0, zdate0 )
      zprec(1) = FLOAT(jpi)
      zprec(2) = FLOAT(jpj)
      zprec(3) = FLOAT(jpisl)
      IF(lwp) WRITE(numout,*) clname
      CALL restini( 'NONE', jpi, jpj, glamt, gphit, 1, zdept,  &
         &          clname, itime, zdate0, rdt, numisl, domain_id=nidom )
      IF( icile == 0 .AND. icut /= 0 ) THEN
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*)' islbsf: write bsfisl in numisl ', numisl
            WRITE(numout,*)' -------------------------------------'
         ENDIF
         zprec(4) = epsisl
         CALL restput(numisl,'PRECISION',1,1,4,0,zprec)
         DO jni = 1, jpisl
            IF(jni < 10) THEN
               WRITE(clisl,'("island",I1)') jni
            ELSE IF(jni < 100) THEN
               WRITE(clisl,'("island",I2)') jni
            ELSE 
               WRITE(clisl,'("island",I3)') jni
            ENDIF
            CALL restput( numisl, clisl, jpi, jpj, 1, 0, bsfisl(:,:,jni) )
         END DO
      ENDIF

      IF( icile < 0 ) THEN
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' islbsf: number of island without convergence : ',ABS(icile)
            WRITE(numout,*) ' ---------------------------------------------'
         ENDIF
         zepsr = epsisl
         DO jni = 1, jpisl
            IF(lwp) WRITE(numout,*) '    isl ',jni,' precision reached ', zep(jni)
            zepsr = MAX( zep(jni), zepsr )
         END DO
         IF( zepsr == 0. ) zepsr = epsisl
         IF(lwp) THEN
            WRITE(numout,*) ' save value of precision reached: ',zepsr
            WRITE(numout,*)
            WRITE(numout,*)' islbsf: save bsfisl in numisl ',numisl
            WRITE(numout,*)' -------------------------------------'
         ENDIF

         zprec(4) = zepsr
         CALL restput( numisl, 'PRECISION', 1, 1, 1, 0, zprec )
         DO jni = 1, jpisl
            IF( jni < 10 ) THEN
               WRITE(clisl,'("island",I1)') jni
            ELSE IF( jni < 100 ) THEN
               WRITE(clisl,'("island",I2)') jni
            ELSE 
               WRITE(clisl,'("island",I3)') jni
            ENDIF
            CALL restput( numisl, clisl, jpi, jpj, 1, 0, bsfisl(:,:,jni) )
         END DO
         CALL restclo(numisl)
         nstop = nstop + 1
      ENDIF

   END SUBROUTINE isl_bsf


   SUBROUTINE isl_dyn_spg
      !!----------------------------------------------------------------------
      !!                  ***  routine isl_dyn_spg  ***
      !!
      !! ** Purpose :   Compute and add the island contribution to the 
      !!      barotropic stream function trend.
      !!
      !!
      !! ** Method  :   Rigid-lid appromimation: ...????
      !!
      !! ** Action : - Update bsfd with the island contribution
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  isolate island computation
      !!---------------------------------------------------------------------

      !! * Local declarations
      INTEGER ::   ji, jj, jni, jnj    ! dummy loop indices
      !!----------------------------------------------------------------------
      

      ! compute the island potential
      ! ----------------------------
      DO jni = 1, jpisl                      ! second member
         bisl(jni) =  0.e0
         DO jj = 2, jpj
            DO ji = 2, jpi
               bisl(jni) =  bisl(jni) + acisl1(ji,jj,1,jni) *  spgu(ji,jj)                  &
                  &                   + acisl1(ji,jj,2,jni) *  spgv(ji,jj)                  &
                  &                   + acisl2(ji,jj,1,jni) * ( gcx(ji,jj)-gcx(ji,jj-1) )   &
                  &                   - acisl2(ji,jj,2,jni) * ( gcx(ji,jj)-gcx(ji-1,jj) )
            END DO
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( bisl, jpisl )   ! sum over the global domain

      DO jni = 1, jpisl                     ! Island stream function trend
         visl(jni) = 0.e0
         DO jnj = 1, jpisl
            visl(jni) = visl(jni) + aislm1(jni,jnj) * bisl(jnj)
         END DO
      END DO

      ! update the bsf trend ( caution : bsfd is not zero along island coastlines, dont mask it ! )
      ! --------------------
      DO jj = 1, jpj
         DO jni = 1, jpisl
            DO ji = 1, jpi
               bsfd(ji,jj) = bsfd(ji,jj) + visl(jni) * bsfisl(ji,jj,jni)
            END DO
         END DO
      END DO

   END SUBROUTINE isl_dyn_spg


   SUBROUTINE isl_stp_ctl( kt, kindic )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE isl_stp_ctl  ***
      !!                     
      !! ** Purpose :   ???
      !!
      !! ** Method  : - print island potential
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  isolated from stp_ctl
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(in   ) ::   kt        ! ocean time-step index
      INTEGER, INTENT(inout) ::   kindic    ! indicator of solver convergence

      !! * local declarations
      INTEGER  ::   jni                     ! dummy loop indice
      REAL(wp) ::   zfact                   ! temporary scalar
      !!----------------------------------------------------------------------
      !!  OPA 8.5, LODYC-IPSL (2002)
      !!----------------------------------------------------------------------

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'isl_stp_ctl : time-stepping control'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      ! Island trends
      DO jni = 1, jpisl
         zfact = 0.
         IF( miisl(1,0,jni) /= 0 .AND. mjisl(1,0,jni) /= 0 ) THEN
            zfact = 1.e-6 * bsfn(miisl(1,0,jni),mjisl(1,0,jni))
         ENDIF
         IF( lk_mpp )   CALL mpp_isl( zfact )

         IF(lwp) WRITE(numisp,9300) kt, jni, zfact, visl(jni)
         IF( MOD( kt, nwrite ) == 0 .OR. kindic < 0     &
            .OR. ( kt == nit000 .AND. kindic > 0 )      &
            .OR.   kt == nitend                         ) THEN
            IF( jni == 1 .AND. lwp ) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'isl_stp_ctl : island bsf'
               WRITE(numout,*) '~~~~~~~~~~~'
            ENDIF
            IF(lwp) WRITE(numout,9300) kt, jni, zfact, visl(jni)
         ENDIF
      END DO
9300  FORMAT(' it : ',i8,' island  :',i4,'   BSF (Sverdrup) : ',f7.2, '   visl : ',e15.6)

   END SUBROUTINE isl_stp_ctl

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_isl = .FALSE.    !: 'key_islands' flag
CONTAINS
   SUBROUTINE isl_dom                        ! Empty routine
   END SUBROUTINE isl_dom
   SUBROUTINE isl_bsf                        ! Empty routine
   END SUBROUTINE isl_bsf
   SUBROUTINE isl_mat                        ! Empty routine
   END SUBROUTINE isl_mat
   SUBROUTINE isl_dyn_spg                    ! Empty routine
   END SUBROUTINE isl_dyn_spg
   SUBROUTINE isl_stp_ctl( kt, kindic )      ! Empty routine
!      WRITE(*,*) 'isl_stp_ctl: You should not have seen this print! error?', kt, kindic
   END SUBROUTINE isl_stp_ctl
#endif

   !!======================================================================
END MODULE solisl
