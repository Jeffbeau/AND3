MODULE obcdom
   !!=================================================================================
   !!                       ***  MODULE  obcdom  ***
   !! Space domain  :  get all the isolated coastline points needed to resolve the 
   !!                  barotropic streamfunction elliptic equation associated with 
   !!                  the open boundaries.
   !!=================================================================================
#if defined key_obc && defined key_dynspg_rl
   !!---------------------------------------------------------------------------------
   !!   'key_obc'           AND                                Open Boundary Condition
   !!   'key_dynspg_rl'                                          Rigid-Lid formulation
   !!---------------------------------------------------------------------------------
   !!   obc_dom        : domain initialization in rid-lid formulation
   !!---------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC obc_dom        ! routine called by iniobc.F90
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_dom
      !!------------------------------------------------------------------------------
      !!                       SUBROUTINE obc_dom
      !!                      ********************
      !! ** Purpose :   Initialize the array used for the computation of the part of 
      !!        the right hand side of the barotropic streamfunction elliptic equation
      !!        associated with the open boundaries
      !!
      !! **  Method :
      !!      + The (i,j) indices of ocean grid-points round isolated coastlines
      !!        are found (isolated coastlines = coast lines separated by an
      !!        open boundary) from icoast array read in coastlines file.
      !!
      !!      + read 'coastline' file  initialize icoast()
      !!        modify icoast() depending on the number of open boundaries 
      !!        specified through key_obc
      !!
      !!      + compute zwb, an ocean/land mask defined as follows:
      !!             zwb(i,j)  =  0. over the one isolated coastline
      !!                       = -1, -2, -3 over the orthers
      !!      + for example, when 4 open boundaries are specified:
      !!
      !!                       //|           |// 
      !!        North          //|   North   |// -1 -1   North
      !!         West     0  0 //| - - - - - |// -1 -1    East
      !!                       //| open bnd  |//
      !!              ///////////|           |/////////
      !!              ------------            ----------
      !!                                               
      !!             west   |                  |    east
      !!           open bnd                       open bnd    
      !!                    |                  |         
      !!              ___________             _________ 
      !!              ///////////|           |/////////
      !!                       //|   south   |// 
      !!        South    -3 -3 //| - - - - - |// -2 -2   South
      !!         West    -3 -3 //| open bnd  |// -2 -2    East
      !!                       //|           |// 
      !!
      !!        With the proper boundary conditions (defined by nperio)
      !!
      !!        C a u t i o n :  no check, the user must enter a well defined
      !!        coastline file. Further more, he must verify that isolated
      !!        coastlines have been well located dans that the right potential
      !!        is affected to the right coastline in obc.F
      !!
      !! History :
      !!   8.1  !  09-97  (J.M. Molines, G. Madec)  Original code
      !!   8.2  !  06-99  (J.M. Molines) suppress zwb(,) for ATL6 (memory saving)
      !!        !  02-02  (A.M. Treguier) icoast in 2 dimension
      !!   8.5  !  02-08  (G. Madec)  F90 : free form
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj, jn, jnic, jnp, jii    ! dummy loop indices
      INTEGER ::   inum = 11         ! temporary logical unit
      INTEGER ::   ifreq, il1, il2, ii, ij, icheck
      INTEGER ::   ip, ipn, ips, ipe, ipw
      INTEGER ::   iim, ijm, iii, ijj
      INTEGER, DIMENSION(jpidta,jpjdta) ::   icoast
      CHARACTER (len=15) ::   clexp
      REAL(wp) ::   zzic, zland
      REAL(wp) ::   zwb, zwbn, zwbs, zwbe, zwbw
      REAL(wp) ::   zglo(jpiglo,jpjglo)
      !!---------------------------------------------------------------------
      !!  OPA 9.0 , LOCEAN-IPSL (2005) 
      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdom.F90,v 1.4 2005/03/27 18:35:10 opalod Exp $ 
      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
      !!---------------------------------------------------------------------
      
      ! 0. initialization of gcfobc to zero
      ! -----------------------------------
      
      DO jn = 1, 3
         gcfobc(:,:,jn) = 0.e0
      END DO
      
      ! 1. Only 1 open boundary : gcfobc is zero, return
      ! ------------------------------------------------
      
      IF( nbobc == 1 .OR. nbic == 0 ) THEN 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' obc_dom: No isolated coastlines gcfobc is set to zero'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
         nstop = nstop + 1
      END IF

      ! 2. Lecture of 'coastlines' file
      ! -------------------------------
      
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'obc_dom: define isolated coastlines from "coastlines" file'
      IF(lwp) WRITE(numout,*) '~~~~~~~'
      IF(lwp) WRITE(numout,*)

      ! open coastlines file'
      CALL ctlopn( inum, 'coastlines', 'OLD', 'FORMATTED', 'SEQUENTIAL',   &
                   1 , numout, lwp, 1 )

      ! lecture of coastlines, set icoast array
      ! Note that this is coded for jpjdta > 1000
      REWIND(inum)
      READ(inum,9101) clexp, iim, ijm
      READ(inum,'(/)')
      ifreq = 40
      il1 = 1
      IF( jpjglo < 1000 ) THEN
         DO jn = 1, jpidta/ifreq+1
            READ(inum,'(/)')
            il2 = min0( jpidta, il1+ifreq-1 )
            READ(inum,9201) ( ii, ji = il1, il2, 5 )
            READ(inum,'(/)')
            DO jj = jpjdta, 1, -1
               READ(inum,9202) ij, ( icoast(ji,jj), ji = il1, il2 )
            END DO
            il1 = il1 + ifreq
         END DO
      ELSE
         DO jn = 1, jpidta/ifreq+1
            READ(inum,'(/)')
            il2 = min0( jpidta, il1+ifreq-1 )
            READ(inum,9221) ( ii, ji = il1, il2, 5 )
            READ(inum,'(/)')
            DO jj = jpjdta, 1, -1
               READ(inum,9222) ij, ( icoast(ji,jj), ji = il1, il2 )
            END DO
            il1 = il1 + ifreq
         END DO
      END IF
      CLOSE(inum)

   ! in case of zoom, icoast must be set to 0 on the domain border
   ! it must be the same for the bathymetry
   IF (lzoom_w) icoast(jpiglo            ,:) = 0 
   IF (lzoom_e) icoast(jpiglo +jpizoom -1,:) = 0 
   IF (lzoom_s) icoast(:,jpjzoom           ) = 0 
   IF (lzoom_n) icoast(:,jpjglo+jpjzoom -1 ) = 0 

      DO jj = 1, jpjglo
         DO ji = 1, jpiglo
            zglo(ji,jj) = icoast( ji+jpizoom-1, jj+jpjzoom-1)
         END DO
      END DO

 9101 FORMAT(1x,a15,2i8)
 9201 FORMAT(3x,13(i3,12x))
 9202 FORMAT(i3,41i3)
 9221 FORMAT(4x,13(i3,12x))
 9222 FORMAT(i4,41i3)

      ! check consistency between tmask and icoast
      
      icheck = 0
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
            icheck = icheck +  INT( tmask(ji,jj,1) ) - MAX(  0, icoast( mig(ji), mjg(jj) )  )
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum(icheck)   ! sum over the global domain

      IF( icheck /= 0 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) 'obc_dom : tmask and isolated coastlines mask are not equal', icheck
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         nstop = nstop + 1
      END IF

      ! 3. transfer the coastline information from T- to f-points
      !    (i.e. from icoast to zwb with  zwb=0 over the continent
      !     and ocean, =-n over the nth isolated coastline)
      ! -----------------------------------------------------------
      
      ! east open boundary
      IF( lp_obc_east .AND. ( jpieob /= 0 ) ) THEN
         IF(lwp) WRITE(numout,*) '         East open boundary: from coastline S.E : ', &
                                 INT(zglo(jpieob,jpjed)),' to N.E : ',               &
                                 INT(zglo(jpieob,jpjef))
      END IF
      ! west open boundary
      IF( lp_obc_west .AND. ( jpiwob /= 0 ) ) THEN
         IF(lwp) WRITE(numout,*) '         West open boundary: from coastline S.W : ', &
                                 INT(zglo(jpiwob,jpjwd)),' to N.W : ',               &
                                 INT(zglo(jpiwob,jpjwf))
      END IF
      ! north open boundary
      IF( lp_obc_north .AND. ( jpjnob /= 0 ) ) THEN
         IF(lwp) WRITE(numout,*) '         North open boundary: from coastline N.W : ', &
                                 INT(zglo(jpind,jpjnob)),' to N.E : ',                &
                                 INT(zglo(jpinf,jpjnob))
      END IF
      ! south open boundary
      IF( lp_obc_south .AND. ( jpjsob /= 0 ) ) THEN
         IF(lwp) WRITE(numout,*) '         South open boundary: from coastline S.W : ', &
                                 INT(zglo(jpisd,jpjsob)),' to S.E : ',                &
                                 INT(zglo(jpisf,jpjsob))
      END IF

      ! 4. Identify the isolated coastline grid point position
      ! ------------------------------------------------------

      ! Loop over isolated coastlines

      DO jnic = 1, nbobc-1
         ! set to zero of miic, mjic of the jnic isolated coastline
         DO jn = 0, 4
            DO ji = 1, jpnic
               miic(ji,jn,jnic) = 0
               mjic(ji,jn,jnic) = 0
            END DO
         END DO
         
         ! Coastal isolated coastline grid-points (miic,mjic)
         ip  = 0
         ipn = 0
         ips = 0
         ipe = 0
         ipw = 0
         
         ! Middle lines (1=<jj=<jpjm1)

         !       jj+1  --zwb--v--ZWB--v--zwb--
         !                |       |       |
         !          jj+1  u   T   u   T   u
         !                |       |       |
         !       jj    --ZWB--v--ZWB--v--ZWB--
         !                |       |       |
         !              jj    u   T   u   T   u
         !                |       |       |
         !       jj-1  --zwb--v--ZWB--v--zwb--
         !                |       |       | 
         !                |   ii  | ii+1  |
         !                |       |       | 
         !               ii-1    ii      ii+1
         
         DO jj = 1, jpjglo-1
            DO ji = 1, jpiglo-1
               ii = ji
               zwb = MIN( 0., zglo(ji,jj), zglo(ji+1,jj), zglo(ji,jj+1), zglo(ji+1,jj+1) )
               IF( jj == jpjglo -1 ) THEN
                  zwbn = zwb
               ELSE
                  zwbn= MIN( 0., zglo(ji,jj+1), zglo(ji+1,jj+1), zglo(ji,jj+2), zglo(ji+1,jj+2) )
               END IF
               IF( jj == 1 ) THEN
                  zwbs = zwb
               ELSE
                  zwbs= MIN( 0., zglo(ji,jj-1), zglo(ji+1,jj-1), zglo(ji,jj), zglo(ji+1,jj) )
               END IF
               IF( ji == jpiglo -1 ) THEN
                  zwbe = zwb
               ELSE
                  zwbe= MIN( 0., zglo(ji+1,jj), zglo(ji+2,jj), zglo(ji+1,jj+1), zglo(ji+2,jj+1) )
               END IF
               IF( ji == 1 ) THEN
                  zwbw = zwb
               ELSE
                  zwbw= MIN( 0., zglo(ji-1,jj), zglo(ji,jj), zglo(ji-1,jj+1), zglo(ji,jj+1) )
               END IF
               
               ! inside coastlines indicator
               zzic  =                zwbn                 &
                     * zwbw                       * zwbe   &
                                    * zwbs
               ! inside land indicator
               zland = MAX( 0., zglo(ji,jj+1) ) + MAX( 0., zglo(ji+1,jj+1) )   &
                     + MAX( 0., zglo(ji,jj  ) ) + MAX( 0., zglo(ji+1,jj  ) )
               ! if isolated coastline grid-point 
               IF( zwb == float( -jnic ) .AND.   &
                  ! not inside the isolated coastline
                             zzic == 0.  .AND.   &
                  ! not inside the land
                             zland >= 2.         ) THEN
                  ! coastal point of the isolated coastline jnic
                  ip = ip + 1
                  miic(ip,0,jnic) = ii
                  mjic(ip,0,jnic) = jj
                  ! which has a west ocean grid point 
                  IF( zwbw == 0. ) THEN
                     ipw = ipw + 1
                     miic(ipw,4,jnic) = ii
                     mjic(ipw,4,jnic) = jj
                  END IF
                  ! which has a east ocean grid point
                  IF( zwbe == 0. ) THEN
                     ipe = ipe + 1
                     IF( nperio == 1 .AND. ii == jpiglo-1 ) THEN
                        miic(ipe,3,jnic) = 2
                     ELSE
                        miic(ipe,3,jnic) = ii + 1
                     END IF
                     mjic(ipe,3,jnic) = jj
                  END IF
                  ! which has a south ocean grid point 
                  IF( zwbs == 0. ) THEN
                     ips = ips + 1
                     miic(ips,2,jnic) = ii
                     mjic(ips,2,jnic) = jj
                  END IF
                  ! which has a north ocean grid point not out of north open b.
                  IF( zwbn == 0. ) THEN
                     ipn = ipn + 1
                     miic(ipn,1,jnic) = ii
                     mjic(ipn,1,jnic) = jj + 1
                  END IF
               END IF
            END DO
         END DO
         
         mnic(0,jnic) = ip
         mnic(1,jnic) = ipn
         mnic(2,jnic) = ips
         mnic(3,jnic) = ipe
         mnic(4,jnic) = ipw
         
      END DO

      ! 5. Check the number of isolated coastline
      ! -----------------------------------------

      DO jnic = 1, nbobc-1
         IF( mnic(0,jnic) > jpnic ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) 'obc_dom: isolated coastline ',jnic,   &
               ' has ',ip,' grid-points > ',jpnic 
            IF(lwp) WRITE(numout,*) '~~~~~~~'
            IF(lwp) WRITE(numout,*) ' modify this dimension in obc_dom'
            nstop = nstop + 1
         END IF
         IF( mnic(0,jnic) == 0 ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) 'obc_dom: isolated coastline ',jnic,   &
               ' has 0  grid-points verify coastlines file'
            IF(lwp) WRITE(numout,*) '~~~~~~~'
            nstop = nstop + 1
         END IF
      END DO

      ! 6. Print of isolated coastline parametres and arrays
      ! -----------------------------------------------------

      IF(lwp) WRITE(numout,*) '           '
      IF(lwp) WRITE(numout,*) '         isolated coastlines found:', nbobc - 1

      DO jnic = 1, nbobc-1
         ip  = mnic(0,jnic)
         ipn = mnic(1,jnic)
         ips = mnic(2,jnic)
         ipe = mnic(3,jnic)
         ipw = mnic(4,jnic)
         IF(lwp) THEN
            WRITE(numout,9000) jnic
            WRITE(numout,9010) ip, ipn, ips, ipe, ipw
            WRITE(numout,9020)
            DO jnp = 1, mnic(0,jnic)
               WRITE(numout,9030) jnp,( miic(jnp,ji,jnic)+nimpp-1, mjic(jnp,ji,jnic)+njmpp-1, ji=0,4 )
            END DO
         END IF
         
         ! format

 9000 FORMAT(/,'          isolated coastline number= ',i2)
 9010 FORMAT(/,'          npic=',i4,' npn=',i4,' nps=',i4,' npe=',i4,' npw=',i4)
 9020 FORMAT(/,'              * ic  point *  point n  *  point s  *  point e ','*  point w  *')
 9030 FORMAT('         ',i4,' * (',i4,',',i4,') * (',i4,',',i4,') * (',i4,',',i4,') * (',i4,',',i4,') * (',i4,',',i4,') *')

      END DO

      ! 7. Construct the gcfobc array associated with each isolated coastline
      ! ----------------------------------------------------------------------

      DO jnic = 1, nbobc-1

         ! north and south grid-points
         DO jii = 1, 2
            DO jnp = 1, mnic(jii,jnic)
               ii = miic(jnp,jii,jnic)
               ij = mjic(jnp,jii,jnic)
               ! take only into account gridpoint of the model domain
               IF( ii >= nldi+nimpp-1 .AND. ii <= nlci+nimpp-1 .AND.   &
                   ij >= nldj+njmpp-1 .AND. ij <= nlcj+njmpp-1       ) THEN
                  iii=ii-nimpp+1
                  ijj=ij-njmpp+1
                  gcfobc(iii,ijj-jii+1,jnic) = gcfobc(iii,ijj-jii+1,jnic)   &
                                             - hur(iii,ijj) * e1u(iii,ijj) / e2u(iii,ijj)

               END IF
            END DO
         END DO

         ! east and west grid-points
         DO jii = 3, 4
            DO jnp = 1, mnic(jii,jnic)
               ii = miic(jnp,jii,jnic)
               ij = mjic(jnp,jii,jnic)
               ! take only into account gridpoint of the model domain
               IF( ii >= nldi+nimpp-1 .AND. ii <= nlci+nimpp-1 .AND.    &
                   ij >= nldj+njmpp-1 .AND. ij <= nlcj+njmpp-1 ) THEN
                  iii=ii-nimpp+1
                  ijj=ij-njmpp+1
                  IF( iii-jii+3 == 1 ) THEN
                     ! cyclic east-west boundary
                     gcfobc(jpim1    ,ijj,jnic) = gcfobc(jpim1    ,ijj,jnic)   &
                                                - hvr(iii,ijj) * e2v(iii,ijj) / e1v(iii,ijj)
                  ELSE
                     ! interior points
                     gcfobc(iii-jii+3,ijj,jnic) = gcfobc(iii-jii+3,ijj,jnic)   &
                                                - hvr(iii,ijj) * e2v(iii,ijj) / e1v(iii,ijj)
                  END IF
               END IF
            END DO
         END DO

         ! applied bmask to suppress coastal open boundary influence
         DO jj = 1, jpj
            DO ji = 1, jpi
               gcfobc(ji,jj,jnic) = gcfobc(ji,jj,jnic) * bmask(ji,jj)
            END DO
         END DO
         
      END DO


      ! 8. check the grid point which value controls the isolated coastline potential
      !    Note: in order to activate those tests you need to make zwb a global array,
      !    which is not done usually to spare memory.
      !    n.b. here at least 2 open boundaries
      ! ------------------------------------------------------------------------------
!
! east open boundary:
!     IF( nieob /= 0 ) THEN
!   east open & south open  :                              Ed === Sf
!         IF( njsob /= 0 ) THEN
!             IF( zwb(nieob,jped) /= zwb(jpsf,njsob) ) THEN
!                 IF(lwp)WRITE(numout,*) ' E R R O R : east d # south f'
!             END IF
!   east open, south closed & west open :                  Ed === Wd
!           ELSEIF( niwob /= 0 ) THEN
!             IF( zwb(nieob,jped) /= zwb(niwob,jpwd) ) THEN
!                 IF(lwp)WRITE(numout,*) ' E R R O R : east d # west d'
!             END IF
!   east open, south closed, west closed & north open :    Ed === Nd
!           ELSEIF( njnob /= 0 ) THEN
!             IF( zwb(nieob,jped) /= zwb(jpnd,njnob) ) THEN
!                 IF(lwp)WRITE(numout,*) ' E R R O R : east d # north d'
!             END IF
!         END IF
!   east open & north open :                               Ef === Nf
!         IF( njnob /= 0 ) THEN
!             IF( zwb(nieob,jpef) /= zwb(jpnf,njnob) ) THEN
!                 IF(lwp)WRITE(numout,*) ' E R R O R : east f # north f'
!             END IF
!   east open, north closed & west open :                  Ef === Wf
!           ELSEIF( niwob /= 0 ) THEN
!             IF( zwb(nieob,jpef) /= zwb(niwob,jpwf) ) THEN
!                 IF(lwp)WRITE(numout,*) ' E R R O R : east f # west f'
!             END IF
!   east open, north closed, west closed & south open :    Ef === Sd
!           ELSEIF( njsob /= 0 ) THEN
!             IF( zwb(nieob,jpef) /= zwb(jpsd,njnob) ) THEN
!                 IF(lwp)WRITE(numout,*) ' E R R O R : east f # south d'
!             END IF
!         END IF
!
! east closed 
!       ELSE
! east closed, south open
!         IF( njsob /= 0 ) THEN
!   east closed, south open & west open :                  Sd === Wd
!             IF( niwob /= 0 ) THEN
!                 IF( zwb(jpsd,njsob) /= zwb(niwob,jpwd) ) THEN
!                     IF(lwp)WRITE(numout,*) ' E R R O R :',
!    $                                       ' south d # west d'
!                 END IF
!   east closed, south open, west closed & north open :    Sd === Nd
!               ELSEIF( njnob /= 0 ) THEN
!                 IF( zwb(jpsd,njsob) /= zwb(jpnd,njnob) ) THEN
!                     IF(lwp)WRITE(numout,*) ' E R R O R : ',
!    $                                       ' south d # north d'
!                 END IF
!             END IF
!   south open, east closed & north open :                 Sf === Nf
!             IF( njnob /= 0 ) THEN
!                 IF( zwb(jpsf,njsob) /= zwb(jpnf,njnob) ) THEN
!                     IF(lwp)WRITE(numout,*) ' E R R O R : ',
!    $                                       ' south f # north f'
!                 END IF
!   south open, east closed, north closed & west open :    Sf === Wf
!               ELSEIF( niwob /= 0 ) THEN
!                 IF( zwb(jpsf,njsob) /= zwb(niwob,jpwf) ) THEN
!                     IF(lwp)WRITE(numout,*) ' E R R O R : ',
!    $                                       ' south f # west f'
!                 END IF
!             END IF
!
! east & south closed ==> north & west open :              Nd === Wf
!                                                              Nf === Wd
!           ELSE
!             IF( zwb(jpnd,njnob) /= zwb(niwob,jpwf) ) THEN
!                 IF(lwp)WRITE(numout,*) ' E R R O R : north d # west f'
!             END IF
!             IF( zwb(jpnf,njnob) /= zwb(niwob,jpwd) ) THEN
!                 IF(lwp)WRITE(numout,*) ' E R R O R : north f # west d'
!             END IF
!         END IF
!
!     END IF
!
!

   END SUBROUTINE obc_dom
#else
   !!=================================================================================
   !!                       ***  MODULE  obcdom  ***
   !! Space domain :  get all the isolated coastline points needed to resolve the 
   !!                 barotropic streamfunction elliptic equation associated with 
   !!                 the open boundaries.
   !!=================================================================================
CONTAINS

   SUBROUTINE obc_dom                 

      ! No isolated coastline OR No Open Boundaries ==> empty routine

   END SUBROUTINE obc_dom
#endif

END MODULE obcdom
