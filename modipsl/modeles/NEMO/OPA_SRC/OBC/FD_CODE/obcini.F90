 MODULE obcini
   !!=================================================================================
   !!                       ***  MODULE  obcini  ***
   !! OBC initial state :  Open boundary initial state
   !!=================================================================================
#if defined key_obc
   !!---------------------------------------------------------------------------------
   !!   'key_obc'                                             Open Boundary Conditions
   !!---------------------------------------------------------------------------------
   !!   obc_init       : initialization for the open boundary condition
   !!---------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lib_mpp         ! for mpp_sum
   USE in_out_manager  ! I/O units

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC obc_init        ! routine called by opa.F90

   !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcini.F90,v 1.9 2006/05/11 15:15:31 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!---------------------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE obc_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE obc_init  ***
      !!         
      !! ** Purpose :   Initialization of the dynamics and tracer fields at 
      !!      the open boundaries.
      !!
      !! ** Method  :   initialization of open boundary variables
      !!      (u, v, bsf) over 3 time step and 3 rows
      !!      (t, s)      over 2 time step and 2 rows
      !!      if ln_rstart = .FALSE. : no restart, fields set to zero
      !!      if ln_rstart = .TRUE.  : restart, fields are read in a file 
      !!      if rdpxxx = 0 then lfbc is set true for this boundary.
      !!
      !! ** Input   :   restart.obc file, restart file for open boundaries 
      !!
      !! History :
      !!   8.0  !  97-07  (G. Madec)  Original code
      !!        !  97-11  (J.M. Molines)
      !!   8.5  !  02-11  (C. Talandier, A-M. Treguier) Free surface, F90
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE obcrst,   ONLY :   obc_rst_lec   ! Make obc_rst_lec routine available
      USE obcdom,   ONLY :   obc_dom       ! Make obc_dom routine available

      !! * Local declarations
      INTEGER  ::   ji, jj, istop , inumfbc
      INTEGER, DIMENSION(4) ::   icorner
      REAL(wp) ::   zbsic1, zbsic2, zbsic3
      REAL(wp), DIMENSION(2) ::   ztestmask
! FD needed for MD double double precision in summation
      COMPLEX(wp) :: ctmp
      REAL(wp) :: tmp1
      INTEGER :: itype                !MD: tmp var for DDPDD

      NAMELIST/namobc/ rdpein, rdpwin, rdpnin, rdpsin,   &
         &             rdpeob, rdpwob, rdpnob, rdpsob,   &
         &             zbsic1, zbsic2, zbsic3,           &
         &             nbic, volemp, nobc_dta,           &
         &             ln_obc_clim, ln_vol_cst, ln_obc_fla
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'obc_init : initialization of open boundaries'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'


      ! 0. read namelist parameters
      ! ---------------------------
      ! default values already set except:
      zbsic1 = 0.e0
      zbsic2 = 0.e0
      zbsic3 = 0.e0

      ! Namelist namobc : open boundaries
      REWIND( numnam )
      READ  ( numnam, namobc )

      bsfic0(1) = zbsic1
      bsfic (2) = zbsic2
      bsfic (3) = zbsic3

      ! By security we set rdpxin and rdpxob respectively
      ! to 1. and 15. if the corresponding OBC is not activated
      IF( .NOT.lp_obc_east ) THEN
         rdpein = 1.
         rdpeob = 15.
      END IF
      IF( .NOT.lp_obc_west ) THEN
         rdpwin = 1.
         rdpwob = 15.
      END IF
      IF( .NOT.lp_obc_north ) THEN
         rdpnin = 1.
         rdpnob = 15.
      END IF
      IF( .NOT.lp_obc_south ) THEN
         rdpsin = 1.
         rdpsob = 15.
      END IF

      ! number of open boudaries and open boundary indicators
      nbobc = 0
      IF( lp_obc_east  )   nbobc = nbobc + 1
      IF( lp_obc_west  )   nbobc = nbobc + 1
      IF( lp_obc_north )   nbobc = nbobc + 1
      IF( lp_obc_south )   nbobc = nbobc + 1

      IF(lwp) WRITE(numout,*) '         Number of open boundaries    nbobc = ',nbobc
      IF(lwp) WRITE(numout,*)
      IF( nbobc /= 0 .AND. jperio /= 0 ) &
           &   CALL ctl_stop( ' Cyclic or symmetric, and open boundary condition are not compatible' )

      ! control prints
      IF(lwp) WRITE(numout,*) '         namobc'
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '         data in file (=1) or     nobc_dta = ', nobc_dta
      IF(lwp) WRITE(numout,*) '         initial state used (=0)             '
      IF(lwp) WRITE(numout,*) '         climatology (true) or not:', ln_obc_clim
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '                                 WARNING                                                  '
      IF(lwp) WRITE(numout,*) '         Flather"s algorithm is applied with explicit free surface scheme                 '
      IF(lwp) WRITE(numout,*) '         or with free surface time-splitting scheme                                       '
      IF(lwp) WRITE(numout,*) '         Nor radiation neither relaxation is allowed with explicit free surface scheme:   '
      IF(lwp) WRITE(numout,*) '         Radiation and/or relaxation is allowed with free surface time-splitting scheme '
      IF(lwp) WRITE(numout,*) '         depending of the choice of rdpXin = rdpXob  = 0. for open boundaries             '
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '         For the rigid-lid case or the filtered free surface case,                        '
      IF(lwp) WRITE(numout,*) '         radiation, relaxation or presciption of data can be applied                      '
      IF( lwp.AND.lp_obc_east ) THEN
         WRITE(numout,*) '         East open boundary :'
         WRITE(numout,*) '              i index                    jpieob = ', jpieob
         WRITE(numout,*) '              damping time scale (days)  rdpeob = ', rdpeob
         WRITE(numout,*) '              damping time scale (days)  rdpein = ', rdpein
      ENDIF

      IF( lwp.AND.lp_obc_west ) THEN
         WRITE(numout,*) '         West open boundary :'
         WRITE(numout,*) '              i index                    jpiwob = ', jpiwob
         WRITE(numout,*) '              damping time scale (days)  rdpwob = ', rdpwob
         WRITE(numout,*) '              damping time scale (days)  rdpwin = ', rdpwin
      ENDIF

      IF( lwp.AND.lp_obc_north ) THEN
         WRITE(numout,*) '         North open boundary :'
         WRITE(numout,*) '               j index                    jpjnob = ', jpjnob
         WRITE(numout,*) '               damping time scale (days)  rdpnob = ', rdpnob
         WRITE(numout,*) '               damping time scale (days)  rdpnin = ', rdpnin
      ENDIF

      IF( lwp.AND.lp_obc_south ) THEN
         WRITE(numout,*) '         South open boundary :'
         WRITE(numout,*) '               j index                    jpjsob = ', jpjsob
         WRITE(numout,*) '               damping time scale (days)  rdpsob = ', rdpsob
         WRITE(numout,*) '               damping time scale (days)  rdpsin = ', rdpsin
         WRITE(numout,*) ' '
      ENDIF

      ! 1. Initialisation of constants 
      ! ------------------------------

      ! ... convert rdp$ob in seconds
      rdpein = rdpein * rday
      rdpwin = rdpwin * rday
      rdpnin = rdpnin * rday
      rdpsin = rdpsin * rday
      rdpeob = rdpeob * rday
      rdpwob = rdpwob * rday
      rdpnob = rdpnob * rday
      rdpsob = rdpsob * rday
      lfbceast  = .FALSE.
      lfbcwest  = .FALSE.
      lfbcnorth = .FALSE.
      lfbcsouth = .FALSE.
      inumfbc = 0
      ! ... look for Fixed Boundaries (rdp = 0 )
      ! ... When specified, lbcxxx flags are set to TRUE and rdpxxx are set to
      ! ...  a small arbitrary value, (to avoid division by zero further on). 
      ! ...  rdpxxx is not used anymore.
      IF( lp_obc_east )  THEN
         IF( (rdpein+rdpeob) == 0 )  THEN
            lfbceast = .TRUE.
            rdpein = 1e-3
            rdpeob = 1e-3
            inumfbc = inumfbc+1
         ELSEIF ( (rdpein*rdpeob) == 0 )  THEN
            CALL ctl_stop( 'obc_init : rdpein & rdpeob must be both zero or non zero' )
         END IF
      END IF
      IF( lp_obc_west )  THEN
         IF( (rdpwin + rdpwob) == 0 )  THEN
            lfbcwest = .TRUE.
            rdpwin = 1e-3
            rdpwob = 1e-3
            inumfbc = inumfbc+1
         ELSEIF ( (rdpwin*rdpwob) == 0 )  THEN
            CALL ctl_stop( 'obc_init : rdpwin & rdpwob must be both zero or non zero' )
         END IF
      END IF
      IF( lp_obc_north )  THEN
         IF( (rdpnin + rdpnob) == 0 )  THEN
            lfbcnorth = .TRUE.
            rdpnin = 1e-3
            rdpnob = 1e-3
            inumfbc = inumfbc+1
         ELSEIF ( (rdpnin*rdpnob) == 0 )  THEN
            CALL ctl_stop( 'obc_init : rdpnin & rdpnob must be both zero or non zero' )
         END IF
      END IF
      IF( lp_obc_south )  THEN
         IF( (rdpsin + rdpsob) == 0 )  THEN
            lfbcsouth = .TRUE.
            rdpsin = 1e-3
            rdpsob = 1e-3
            inumfbc = inumfbc+1
         ELSEIF ( (rdpsin*rdpsob) == 0 )  THEN
            CALL ctl_stop( 'obc_init : rdpsin & rdpsob must be both zero or non zero' )
         END IF
      END IF

      ! 2.  Clever mpp indices for loops on the open boundaries. 
      !     The loops will be performed only on the processors 
      !     that contain a given open boundary.
      ! --------------------------------------------------------

      IF( lp_obc_east ) THEN
         ! ...   mpp initialization
         nie0   = max( 1, min(jpieob   - nimpp+1, jpi     ) )
         nie1   = max( 0, min(jpieob   - nimpp+1, jpi - 1 ) )
         nie0p1 = max( 1, min(jpieob+1 - nimpp+1, jpi     ) )
         nie1p1 = max( 0, min(jpieob+1 - nimpp+1, jpi - 1 ) )
         nie0m1 = max( 1, min(jpieob-1 - nimpp+1, jpi     ) )
         nie1m1 = max( 0, min(jpieob-1 - nimpp+1, jpi - 1 ) )
         nje0   = max( 2, min(jpjed    - njmpp+1, jpj     ) )
         nje1   = max( 0, min(jpjef    - njmpp+1, jpj - 1 ) )
         nje0p1 = max( 1, min(jpjedp1  - njmpp+1, jpj     ) )
         nje0m1 = max( 1, min(jpjed    - njmpp+1, jpj     ) )
         nje1m1 = max( 0, min(jpjefm1  - njmpp+1, jpj - 1 ) )
         nje1m2 = max( 0, min(jpjefm1-1- njmpp+1, jpj - 1 ) )
         IF(lwp) THEN
            IF( lfbceast ) THEN
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Specified East Open Boundary'
            ELSE
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Radiative East Open Boundary'
            END IF
         END IF
      END IF

      IF( lp_obc_west ) THEN
         ! ...   mpp initialization
         niw0   = max( 1, min(jpiwob   - nimpp+1, jpi     ) )
         niw1   = max( 0, min(jpiwob   - nimpp+1, jpi - 1 ) )
         niw0p1 = max( 1, min(jpiwob+1 - nimpp+1, jpi     ) )
         niw1p1 = max( 0, min(jpiwob+1 - nimpp+1, jpi - 1 ) )
         njw0   = max( 2, min(jpjwd    - njmpp+1, jpj     ) )
         njw1   = max( 0, min(jpjwf    - njmpp+1, jpj - 1 ) )
         njw0p1 = max( 1, min(jpjwdp1  - njmpp+1, jpj     ) )
         njw0m1 = max( 1, min(jpjwd    - njmpp+1, jpj     ) )
         njw1m1 = max( 0, min(jpjwfm1  - njmpp+1, jpj - 1 ) )
         njw1m2 = max( 0, min(jpjwfm1-1- njmpp+1, jpj - 1 ) )
         IF(lwp) THEN
            IF( lfbcwest ) THEN
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Specified West Open Boundary'
            ELSE
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Radiative West Open Boundary'
            END IF
         END IF
      END IF
 
      IF( lp_obc_north ) THEN
         ! ...   mpp initialization
         nin0   = max( 2, min(jpind    - nimpp+1, jpi     ) )
         nin1   = max( 0, min(jpinf    - nimpp+1, jpi - 1 ) )
         nin0p1 = max( 1, min(jpindp1  - nimpp+1, jpi     ) )
         nin0m1 = max( 1, min(jpind    - nimpp+1, jpi     ) )
         nin1m1 = max( 0, min(jpinfm1  - nimpp+1, jpi - 1 ) )
         nin1m2 = max( 0, min(jpinfm1-1- nimpp+1, jpi - 1 ) )
         njn0   = max( 1, min(jpjnob   - njmpp+1, jpj     ) )
         njn1   = max( 0, min(jpjnob   - njmpp+1, jpj - 1 ) )
         njn0p1 = max( 1, min(jpjnob+1 - njmpp+1, jpj     ) )
         njn1p1 = max( 0, min(jpjnob+1 - njmpp+1, jpj - 1 ) )
         njn0m1 = max( 1, min(jpjnob-1 - njmpp+1, jpj     ) )
         njn1m1 = max( 0, min(jpjnob-1 - njmpp+1, jpj - 1 ) )
         IF(lwp) THEN
            IF( lfbcnorth ) THEN
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Specified North Open Boundary'
            ELSE
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Radiative North Open Boundary'
            END IF
         END IF
      END IF

      IF( lp_obc_south ) THEN
         ! ...   mpp initialization
         nis0   = max( 2, min(jpisd    - nimpp+1, jpi     ) )
         nis1   = max( 0, min(jpisf    - nimpp+1, jpi - 1 ) )
         nis0p1 = max( 1, min(jpisdp1  - nimpp+1, jpi     ) )
         nis0m1 = max( 1, min(jpisd    - nimpp+1, jpi     ) )
         nis1m1 = max( 0, min(jpisfm1  - nimpp+1, jpi - 1 ) )
         nis1m2 = max( 0, min(jpisfm1-1- nimpp+1, jpi - 1 ) )
         njs0   = max( 1, min(jpjsob   - njmpp+1, jpj     ) )
         njs1   = max( 0, min(jpjsob   - njmpp+1, jpj - 1 ) )
         njs0p1 = max( 1, min(jpjsob+1 - njmpp+1, jpj     ) )
         njs1p1 = max( 0, min(jpjsob+1 - njmpp+1, jpj - 1 ) )
         IF(lwp) THEN
            IF( lfbcsouth ) THEN
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Specified South Open Boundary'
            ELSE
               WRITE(numout,*)'     '
               WRITE(numout,*)'         Radiative South Open Boundary'
            END IF
         END IF
      END IF

! FD: the following test on corners should be done before
! FD: the OB masks and the total lateral surface are set or computed
! FD: so that, any adjustment will be viewed by the rest of the code.
      ! 5. Control print on mask 
      !    The extremities of the open boundaries must be in land
      !    or else correspond to an "ocean corner" between two open boundaries. 
      !    corner 1 is southwest, 2 is south east, 3 is northeast, 4 is northwest. 
      ! --------------------------------------------------------------------------

      icorner(:)=0

      ! ... control of the west boundary
      IF( lp_obc_west ) THEN
         IF( jpiwob < 2 .OR.  jpiwob >= jpiglo-2 ) THEN
            WRITE(ctmp1,*) ' jpiwob exceed ', jpiglo-2, 'or less than 2'
            CALL ctl_stop( ctmp1 )
         END IF
         ztestmask(:)=0.
         DO ji=niw0,niw1
            IF( (njw0 + njmpp - 1) == jpjwd ) ztestmask(1)=ztestmask(1)+ tmask(ji,njw0,1)
            IF( (njw1 + njmpp - 1) == jpjwf ) ztestmask(2)=ztestmask(2)+ tmask(ji,njw1,1)
         END DO
         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain

         IF( ztestmask(1) /= 0. ) icorner(1)=icorner(1)+1
         IF( ztestmask(2) /= 0. ) icorner(4)=icorner(4)+1
      END IF

      ! ... control of the east boundary
      IF( lp_obc_east ) THEN
         IF( jpieob < 4 .OR.  jpieob >= jpiglo ) THEN
            WRITE(ctmp1,*) ' jpieob exceed ', jpiglo, ' or less than 4'
            CALL ctl_stop( ctmp1 )
         END IF
         ztestmask(:)=0.
         DO ji=nie0p1,nie1p1
            IF( (nje0 + njmpp - 1) == jpjed ) ztestmask(1)=ztestmask(1)+ tmask(ji,nje0,1)
            IF( (nje1 + njmpp - 1) == jpjef ) ztestmask(2)=ztestmask(2)+ tmask(ji,nje1,1)
         END DO
         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain

        IF( ztestmask(1) /= 0. ) icorner(2)=icorner(2)+1
        IF( ztestmask(2) /= 0. ) icorner(3)=icorner(3)+1
      END IF

      ! ... control of the north boundary
      IF( lp_obc_north ) THEN
         IF( jpjnob < 4 .OR.  jpjnob >= jpjglo ) THEN
            WRITE(ctmp1,*) 'jpjnob exceed ', jpjglo, ' or less than 4'
            CALL ctl_stop( ctmp1 )
         END IF
         ztestmask(:)=0.
         DO jj=njn0p1,njn1p1
            IF( (nin0 + nimpp - 1) == jpind ) ztestmask(1)=ztestmask(1)+ tmask(nin0,jj,1)
            IF( (nin1 + nimpp - 1) == jpinf ) ztestmask(2)=ztestmask(2)+ tmask(nin1,jj,1)
         END DO
         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain

         IF( ztestmask(1) /= 0. ) icorner(4)=icorner(4)+1
         IF( ztestmask(2) /= 0. ) icorner(3)=icorner(3)+1
      END IF

      ! ... control of the south boundary
      IF( lp_obc_south ) THEN
         IF( jpjsob < 2 .OR.  jpjsob >= jpjglo-2 ) THEN
            WRITE(ctmp1,*) ' jpjsob exceed ', jpjglo-2, ' or less than 2'
            CALL ctl_stop( ctmp1 )
         END IF
         ztestmask(:)=0.
         DO jj=njs0,njs1
            IF( (nis0 + nimpp - 1) == jpisd ) ztestmask(1)=ztestmask(1)+ tmask(nis0,jj,1)
            IF( (nis1 + nimpp - 1) == jpisf ) ztestmask(2)=ztestmask(2)+ tmask(nis1,jj,1)
         END DO
         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain

         IF( ztestmask(1) /= 0. ) icorner(1)=icorner(1)+1
         IF( ztestmask(2) /= 0. ) icorner(2)=icorner(2)+1
      END IF

      IF( icorner(1) == 2 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' South West ocean corner, two open boudaries'
         IF(lwp) WRITE(numout,*) ' ========== '
         IF(lwp) WRITE(numout,*)
         IF( jpisd /= jpiwob.OR.jpjsob /= jpjwd ) &
              &   CALL ctl_stop( ' Open boundaries do not fit, we stop' )

! FD mask the T-point around the corner, otherwise the flux will be counting an extra useless point in obcvol
! WEST:
       IF( lp_obc_west ) THEN
         DO ji=niw0,niw1
            IF( (njw0 + njmpp - 1) == jpjwd ) THEN
               tmask(ji,njw0,:) = 0.e0
               umask(ji,njw0,:) = 0.e0
            ENDIF
         END DO
       END IF
! SOUTH
       IF( lp_obc_south ) THEN
         DO jj=njs0,njs1
            IF( (nis0 + nimpp - 1) == jpisd ) THEN
               tmask(nis0,jj,:) = 0.e0
               vmask(nis0,jj,:) = 0.e0
            ENDIF
         END DO
       END IF

!FD      ELSE IF( icorner(1) == 1 ) THEN
!FD         CALL ctl_stop( ' Open boundaries do not fit at SW corner, we stop' )
      END IF 

      IF( icorner(2) == 2 ) THEN
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) ' South East ocean corner, two open boudaries'
          IF(lwp) WRITE(numout,*) ' ========== '
          IF(lwp) WRITE(numout,*)
          IF( jpisf /= jpieob+1.OR.jpjsob /= jpjed ) &
               &   CALL ctl_stop( ' Open boundaries do not fit, we stop' )

! FD mask the T-point around the corner, otherwise the flux will be counting an extra useless point in obcvol
! EAST:
       IF( lp_obc_east ) THEN
         DO ji=nie0,nie1
            IF( (nje0 + njmpp - 1) == jpjed ) THEN
               tmask(ji+1,nje0,:) = 0.e0
               umask(ji,nje0  ,:) = 0.e0
            ENDIF
         END DO
       END IF
! SOUTH
       IF( lp_obc_south ) THEN
         DO jj=njs0,njs1
            IF( (nis1 + nimpp - 1) == jpisf ) THEN
               tmask(nis1,jj,:) = 0.e0
               vmask(nis1,jj,:) = 0.e0
            ENDIF
         END DO
       END IF

!FD      ELSE IF( icorner(2) == 1 ) THEN
!FD         CALL ctl_stop( ' Open boundaries do not fit at SE corner, we stop' )
      END IF 

      IF( icorner(3) == 2 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' North East ocean corner, two open boudaries'
         IF(lwp) WRITE(numout,*) ' ========== '
         IF(lwp) WRITE(numout,*)
         IF( jpinf /= jpieob+1 .OR. jpjnob+1 /= jpjef ) &
              &   CALL ctl_stop( ' Open boundaries do not fit, we stop' )

! FD mask the T-point around the corner, otherwise the flux will be counting an extra useless point in obcvol
! EAST:
       IF( lp_obc_east ) THEN
         DO ji=nie0,nie1
            IF( (nje1 + njmpp - 1) == jpjef ) THEN
               tmask(ji+1,nje1,:) = 0.e0
               umask(ji,nje1  ,:) = 0.e0
            ENDIF
         END DO
       END IF
! NORTH
       IF( lp_obc_north ) THEN
         DO jj=njn0,njn1
            IF( (nin1 + nimpp - 1) == jpinf ) THEN
              tmask(nin1,jj+1,:) = 0.e0
              vmask(nin1,jj  ,:) = 0.e0
            ENDIF
         END DO
       END IF

!FD       ELSE IF( icorner(3) == 1 ) THEN
!FD          CALL ctl_stop( ' Open boundaries do not fit at NE corner, we stop' )
       END IF 

      IF( icorner(4) == 2 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' North West ocean corner, two open boudaries'
         IF(lwp) WRITE(numout,*) ' ========== '
         IF(lwp) WRITE(numout,*)
         IF( jpind /= jpiwob.OR.jpjnob+1 /= jpjwf ) &
              &   CALL ctl_stop( ' Open boundaries do not fit, we stop' )

! FD mask the T-point around the corner, otherwise the flux will be counting an extra useless point in obcvol
! WEST:
       IF( lp_obc_west ) THEN
         DO ji=niw0,niw1
            IF( (njw1 + njmpp - 1) == jpjwf ) THEN
               tmask(ji,njw1,:) = 0.e0
               umask(ji,njw1,:) = 0.e0
            ENDIF
         END DO
       END IF
! NORTH
       IF( lp_obc_north ) THEN
         DO jj=njn0,njn1
            IF( (nin0 + nimpp - 1) == jpind ) THEN
              tmask(nin0,jj+1,1) = 0.e0
              vmask(nin0,jj  ,1) = 0.e0
            ENDIF
         END DO
       END IF

!FD       ELSE IF( icorner(4) == 1 ) THEN
!FD          CALL ctl_stop( ' Open boundaries do not fit at NW corner, we stop' )
       END IF 

      ! 3. mask correction for OBCs
      ! ---------------------------

      IF( lp_obc_east ) THEN
         !... (jpjed,jpjefm1),jpieob
! FD        DO jj = nje0, nje1m1
         DO jj = nje0, nje1
# if defined key_dynspg_rl
            DO ji = nie0, nie1
# else
            DO ji = nie0p1, nie1p1
# endif
               bmask(ji,jj) = 0.e0
            END DO
         END DO

         ! ... initilization to zero
         uemsk(:,:) = 0.e0
         vemsk(:,:) = 0.e0
         temsk(:,:) = 0.e0

         ! ... set 2D mask on East OBC,  Vopt
         DO ji = fs_nie0, fs_nie1
            DO jj = nje0, nje1
               uemsk(jj,:) = umask(ji,  jj,:)
               vemsk(jj,:) = vmask(ji+1,jj,:)
               temsk(jj,:) = tmask(ji+1,jj,:)
            END DO
         END DO

      END IF

      IF( lp_obc_west ) THEN
         ! ... (jpjwd,jpjwfm1),jpiwob
! FD        DO jj = njw0, njw1m1
         DO jj = njw0, njw1
            DO ji = niw0, niw1
               bmask(ji,jj) = 0.e0
            END DO
         END DO

         ! ... initilization to zero
         uwmsk(:,:) = 0.e0
         vwmsk(:,:) = 0.e0
         twmsk(:,:) = 0.e0

         ! ... set 2D mask on West OBC,  Vopt
         DO ji = fs_niw0, fs_niw1
            DO jj = njw0, njw1
               uwmsk(jj,:) = umask(ji,jj,:)
               vwmsk(jj,:) = vmask(ji,jj,:)
               twmsk(jj,:) = tmask(ji,jj,:)
            END DO
         END DO

      END IF

      IF( lp_obc_north ) THEN
         ! ... jpjnob,(jpind,jpisfm1)
# if defined key_dynspg_rl
         DO jj = njn0, njn1
# else
         DO jj = njn0p1, njn1p1
# endif
! FD           DO ji = nin0, nin1m1
            DO ji = nin0, nin1
               bmask(ji,jj) = 0.e0
            END DO
         END DO

         ! ... initilization to zero
         unmsk(:,:) = 0.e0
         vnmsk(:,:) = 0.e0
         tnmsk(:,:) = 0.e0

         ! ... set 2D mask on North OBC,  Vopt
         DO jj = fs_njn0, fs_njn1
            DO ji = nin0, nin1
               unmsk(ji,:) = umask(ji,jj+1,:)
               vnmsk(ji,:) = vmask(ji,jj  ,:)
               tnmsk(ji,:) = tmask(ji,jj+1,:)
            END DO
         END DO

      END IF

      IF( lp_obc_south ) THEN 
         ! ... jpjsob,(jpisd,jpisfm1)
         DO jj = njs0, njs1
! FD           DO ji = nis0, nis1m1
            DO ji = nis0, nis1
               bmask(ji,jj) = 0.e0
            END DO
         END DO

         ! ... initilization to zero
         usmsk(:,:) = 0.e0
         vsmsk(:,:) = 0.e0
         tsmsk(:,:) = 0.e0

         ! ... set 2D mask on South OBC,  Vopt
         DO jj = njs0, njs1
            DO ji = nis0, nis1
               usmsk(ji,:) = umask(ji,jj,:)
               vsmsk(ji,:) = vmask(ji,jj,:)
               tsmsk(ji,:) = tmask(ji,jj,:)
            END DO
         END DO

      END IF

# if defined key_dynspg_flt

      ! ... Initialize obcumask and obcvmask for the Force filtering 
      !     boundary condition in dynspg_flt
      obcumask(:,:) = umask(:,:,1)
      obcvmask(:,:) = vmask(:,:,1)

      ! ... Initialize obctmsk on overlap region and obcs. This mask
      !     is used in obcvol.F90 to calculate cumulate flux E-P. 
      !     - no flux E-P on obcs and overlap region (jpereci = jprecj = 1)
      obctmsk(:,:) = tmask(:,:,1)     
      obctmsk(1  ,:) = 0.e0           
      obctmsk(jpi,:) = 0.e0           
      obctmsk(:  ,1) = 0.e0           
      obctmsk(:,jpj) = 0.e0           

      IF( lp_obc_east ) THEN
         ! ... East obc Force filtering mask for the grad D
         DO ji = nie0, nie1
! FD           DO jj = nje0p1, nje1m1
            DO jj = nje0, nje1
               obcumask(ji  ,jj)=0.e0
               obcvmask(ji+1,jj)=0.e0
            END DO
         END DO

         ! ... set to 0 on East OBC
! FD        DO jj = nje0p1, nje1m1
         DO jj = nje0, nje1
            DO ji = nie0p1, nie1p1
               obctmsk(ji,jj) = 0.e0
            END DO
         END DO
      END IF

      IF( lp_obc_west ) THEN
         ! ... West obc Force filtering mask for the grad D
         DO ji = niw0, niw1
! FD           DO jj = njw0p1, njw1m1
            DO jj = njw0, njw1
               obcumask(ji,jj)=0.e0
               obcvmask(ji,jj)=0.e0
            END DO
         END DO

         ! ... set to 0 on West OBC
! FD        DO jj = njw0p1, njw1m1
         DO jj = njw0, njw1
            DO ji = niw0, niw1
               obctmsk(ji,jj) = 0.e0
            END DO
         END DO
      END IF

      IF( lp_obc_north ) THEN
         ! ... North obc Force filtering mask for the grad D
         DO jj = njn0, njn1
! FD           DO ji = nin0p1, nin1m1
            DO ji = nin0, nin1
               obcvmask(ji,jj  )=0.e0
               obcumask(ji,jj+1)=0.e0
            END DO
         END DO

         ! ... set to 0 on North OBC
         DO jj = njn0p1, njn1p1
! FD           DO ji = nin0p1, nin1m1
            DO ji = nin0, nin1
               obctmsk(ji,jj) = 0.e0
            END DO
         END DO
      END IF

      IF( lp_obc_south ) THEN
         ! ... South obc Force filtering mask for the grad D
         DO jj = njs0, njs1 
! FD           DO ji = nis0p1, nis1m1
            DO ji = nis0, nis1
               obcumask(ji,jj)=0.e0
               obcvmask(ji,jj)=0.e0
            END DO
         END DO

         ! ... set to 0 on South OBC
         DO jj = njs0, njs1
! FD           DO ji = nis0p1, nis1m1
            DO ji = nis0, nis1
               obctmsk(ji,jj) = 0.e0
            END DO
         END DO
      END IF

# endif

# if ! defined key_dynspg_rl

! FD: better compute this in any case
! FD  IF ( ln_vol_cst ) THEN

         ! 3.1 Total lateral surface for each open boundary
         ! ------------------------------------------------
!MD START
         ctmp = cmplx(0.e0,0.e0,wp)

         ! ... West open boundary surface
         IF( lp_obc_west ) THEN
            DO ji = niw0, niw1
               DO jj = 1, jpj 
! FD                  obcsurftot = obcsurftot+hu(ji,jj)*e2u(ji,jj)*uwmsk(jj,1)*tmask_i(ji,jj)
                  tmp1 = hu(ji,jj)*e2u(ji,jj)*uwmsk(jj,1)*tmask_i(ji,jj)
                  CALL DDPDD( cmplx(tmp1,0.e0,wp), ctmp, 1, itype)
               END DO
            END DO
         END IF

         ! ... East open boundary surface
         IF( lp_obc_east ) THEN
            DO ji = nie0, nie1
               DO jj = 1, jpj 
! FD                  obcsurftot = obcsurftot+hu(ji,jj)*e2u(ji,jj)*uemsk(jj,1)*tmask_i(ji+1,jj)
                  tmp1 = hu(ji,jj)*e2u(ji,jj)*uemsk(jj,1)*tmask_i(ji+1,jj)
                  CALL DDPDD( cmplx(tmp1,0.e0,wp), ctmp, 1, itype)
               END DO
            END DO
         END IF

         ! ... North open boundary vertical surface
         IF( lp_obc_north ) THEN
            DO jj = njn0, njn1
               DO ji = 1, jpi
! FD                  obcsurftot = obcsurftot+hv(ji,jj)*e1v(ji,jj)*vnmsk(ji,1)*tmask_i(ji,jj+1)
                  tmp1 = hv(ji,jj)*e1v(ji,jj)*vnmsk(ji,1)*tmask_i(ji,jj+1)
                  CALL DDPDD( cmplx(tmp1,0.e0,wp), ctmp, 1, itype)
               END DO
            END DO
         END IF

         ! ... South open boundary vertical surface
         IF( lp_obc_south ) THEN
            DO jj = njs0, njs1
               DO ji = 1, jpi
! FD                  obcsurftot = obcsurftot+hv(ji,jj)*e1v(ji,jj)*vsmsk(ji,1)*tmask_i(ji,jj)
                  tmp1 = hv(ji,jj)*e1v(ji,jj)*vsmsk(ji,1)*tmask_i(ji,jj)
                  CALL DDPDD( cmplx(tmp1,0.e0,wp), ctmp, 1, itype)
               END DO
            END DO
         END IF
! FD         IF( lk_mpp )   CALL mpp_sum( obcsurftot )   ! sum over the global domain
!MD END
         IF( lk_mpp )   CALL mpp_sum( ctmp )   ! sum over the global domain
         obcsurftot = real( ctmp,wp )
! FD  ENDIF
# endif
      ! 6. Initialization of open boundary variables (u, v, bsf, t, s)
      ! --------------------------------------------------------------
      !   only if at least one boundary is  radiative 

      ! ... Restart from restart.obc
      IF ( inumfbc < nbobc .AND.  ln_rstart ) THEN
         CALL obc_rst_lec
      ELSE

          ! ... Initialization to zero of radiation arrays.
          !     Those have dimensions of local subdomains

          bebnd(:,:,:)   = 0.e0   ;   bnbnd(:,:,:)   = 0.e0
          uebnd(:,:,:,:) = 0.e0   ;   unbnd(:,:,:,:) = 0.e0
          vebnd(:,:,:,:) = 0.e0   ;   vnbnd(:,:,:,:) = 0.e0
          tebnd(:,:,:,:) = 0.e0   ;   tnbnd(:,:,:,:) = 0.e0 
          sebnd(:,:,:,:) = 0.e0   ;   snbnd(:,:,:,:) = 0.e0

          bwbnd(:,:,:)   = 0.e0   ;   bsbnd(:,:,:)   = 0.e0
          uwbnd(:,:,:,:) = 0.e0   ;   usbnd(:,:,:,:) = 0.e0
          vwbnd(:,:,:,:) = 0.e0   ;   vsbnd(:,:,:,:) = 0.e0
          twbnd(:,:,:,:) = 0.e0   ;   tsbnd(:,:,:,:) = 0.e0 
          swbnd(:,:,:,:) = 0.e0   ;   ssbnd(:,:,:,:) = 0.e0

      END IF

# if defined key_dynspg_rl
      ! 7. Isolated coastline arrays initialization (rigid lid case only)
      ! -----------------------------------------------------------------
      CALL obc_dom
# endif

      ! 8. Control print
      ! ... control of the east boundary
      IF( lp_obc_east ) THEN
         istop = 0
         IF( jpieob < 4 .OR.  jpieob >= jpiglo ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '            jpieob exceed ', jpim1, ' or less than 4'
            istop = istop + 1
         END IF

         IF( lk_mpp ) THEN
            ! ... 
            IF( nimpp > jpieob-5) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the East OBC'
               IF(lwp) WRITE(numout,*) '        nimpp must be < jpieob-5'
               istop = istop + 1
            ENDIF
         ELSE

            ! ... stop if  e r r o r (s)   detected
            IF( istop /= 0 ) THEN
               WRITE(ctmp1,*) istop,' obcini : E R R O R (S) detected : stop'
               CALL ctl_stop( ctmp1 )
            ENDIF
         ENDIF
      ENDIF

      ! ... control of the west boundary
      IF( lp_obc_west ) THEN
         istop = 0
         IF( jpiwob < 2 .OR.  jpiwob >= jpiglo ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '            jpiwob exceed ', jpim1, ' or less than 2'
            istop = istop + 1
         END IF

         IF( lk_mpp ) THEN
            IF( (nimpp < jpiwob+5) .AND. (nimpp > 1) ) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the West OBC'
               IF(lwp) WRITE(numout,*) '        nimpp must be > jpiwob-5 or =1'
               istop = istop + 1
            ENDIF
         ELSE
   
            ! ... stop if  e r r o r (s)   detected
            IF( istop /= 0 ) THEN
               WRITE(ctmp1,*) istop,' obcini : E R R O R (S) detected : stop'
               CALL ctl_stop( ctmp1 )
            ENDIF
         ENDIF
      ENDIF

      ! control of the north boundary
      IF( lp_obc_north ) THEN
         istop = 0
         IF( jpjnob < 4 .OR.  jpjnob >= jpjglo ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          jpjnob exceed ', jpjm1,' or less than 4'
            istop = istop + 1
         END IF

         IF( lk_mpp ) THEN
            IF( njmpp > jpjnob-5) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the North OBC'
               IF(lwp) WRITE(numout,*) '        njmpp must be < jpjnob-5'
               istop = istop + 1
            ENDIF
         ELSE
   
            ! ... stop if  e r r o r (s)   detected
            IF( istop /= 0 ) THEN
                WRITE(ctmp1,*) istop,' obcini : E R R O R (S) detected : stop'
               CALL ctl_stop( ctmp1 )
           ENDIF
         ENDIF
      ENDIF

      ! control of the south boundary
      IF( lp_obc_south ) THEN
         istop = 0
         IF( jpjsob < 2 .OR. jpjsob >= jpjglo ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          jpjsob exceed ', jpjm1,' or less than 2'
            istop = istop + 1
         END IF

         IF( lk_mpp ) THEN
            IF( (njmpp < jpjsob+5) .AND. (njmpp > 1) ) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the South OBC'
               IF(lwp) WRITE(numout,*) '        njmpp must be > jpjsob+5 or =1'
               istop = istop + 1
            ENDIF
         ELSE
   
            ! ... stop if  e r r o r (s)   detected
            IF( istop /= 0 ) THEN
               WRITE(ctmp1,*) istop,' obcini : E R R O R (S) detected : stop'
               CALL ctl_stop( ctmp1 )
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE obc_init

#else
   !!---------------------------------------------------------------------------------
   !!   Dummy module                                                NO open boundaries
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE obc_init      ! Dummy routine
   END SUBROUTINE obc_init
#endif

   !!=================================================================================
END MODULE obcini
