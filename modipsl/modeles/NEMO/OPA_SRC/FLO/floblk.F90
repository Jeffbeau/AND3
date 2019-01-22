!!DB 2008.03.19 ...
!!Added a Random Displacement Model (RDM) = particle tracking plus turbulence
!!Routine is called flo_RDM. It replaces the flo_4rk() routine if(ln_flork4)
!!NB: the Rk4 routine did not work for me.
!!To use RDM scheme: in namelist set ln_flork4 = .true. 
!!Requires a init_float file with float positions which I have also modified. 
!!See flodom.F90 and look for DB
!!NB: Routine works but it is not really complete. For example:
!!   (a) it ignores w-direction, i.e. floats remain at initial level
!!   (b) following from (a): there is no vertical turbulence
!!   (c) Add more notes ... TO DO ...


MODULE floblk
   !!======================================================================
   !!                     ***  MODULE  floblk  ***
   !! Ocean floats :   trajectory computation
   !!======================================================================
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                     float trajectories
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!    flotblk     : compute float trajectories with Blanke algorithme
   !!----------------------------------------------------------------------
   !! * Modules used
   USE flo_oce         ! ocean drifting floats
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
!!DB
   USE obc_oce         ! ocean open boundary conditions and other definitions
   USE ldfdyn_oce          ! to get momentum diffusivities ??? 
   USE ldftra_oce          ! to get tracer diffusivities ??? 

   IMPLICIT NONE
   PRIVATE

! CN: For Mersenne twister RNG
! Default seed
    integer, parameter :: defaultsd = 4357
! Period parameters
    integer, parameter :: N = 624, N1 = N + 1

! the array for the state vector
    integer, save, dimension(0:N-1) :: mt
    integer, save                   :: mti = N1

   !! * Accessibility
!!DB
   PRIVATE gasdev      ! gauss function 
   PUBLIC flo_blk      ! routine called by floats.F90
   PUBLIC flo_RDM      ! routine called by floats.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/floblk.F90,v 1.4 2005/03/27 18:35:05 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE flo_blk( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_blk  ***
      !!           
      !! ** Purpose :   Compute the geographical position,latitude, longitude
      !!      and depth of each float at each time step.
      !! 
      !! ** Method  :   The position of a float is computed with Bruno Blanke
      !!      algorithm. We need to know the velocity field, the old positions
      !!      of the floats and the grid defined on the domain.
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt ! ocean time step

      !! * Local declarations
      INTEGER :: jfl              ! dummy loop arguments
      INTEGER :: ind, ifin, iloop
!!DB
!      INTEGER , DIMENSION ( jpnfl )  ::   &
      INTEGER , DIMENSION (:),ALLOCATABLE, SAVE  ::   &
         iil, ijl, ikl,             &     ! index of nearest mesh
         iiloc , ijloc,             &
         iiinfl, ijinfl, ikinfl,    &     ! index of input mesh of the float.
         iioutfl, ijoutfl, ikoutfl        ! index of output mesh of the float.
!!DB
!      REAL(wp) , DIMENSION ( jpnfl )  ::    &
      REAL(wp) , DIMENSION (:),ALLOCATABLE, SAVE  ::    &
         zgifl, zgjfl, zgkfl,       &     ! position of floats, index on 
                                          ! velocity mesh.
         ztxfl, ztyfl, ztzfl,       &     ! time for a float to quit the mesh
                                          ! across one of the face x,y and z 
         zttfl,                     &     ! time for a float to quit the mesh 
         zagefl,                    &     ! time during which, trajectorie of 
                                          ! the float has been computed
         zagenewfl,                 &     ! new age of float after calculation 
                                          ! of new position
         zufl, zvfl, zwfl,          &     ! interpolated vel. at float position
         zudfl, zvdfl, zwdfl,       &     ! velocity diff input/output of mesh
         zgidfl, zgjdfl, zgkdfl           ! direction index of float 
      REAL(wp)   ::       &
         zuinfl,zvinfl,zwinfl,      &     ! transport across the input face
         zuoutfl,zvoutfl,zwoutfl,   &     ! transport across the ouput face
         zvol,                      &     ! volume of the mesh
         zsurfz,                    &     ! surface of the face of the mesh 
         zind
      REAL(wp), DIMENSION ( 2 )  ::   &
         zsurfx, zsurfy                   ! surface of the face of the mesh


      !!---------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'flo_blk : compute Blanke trajectories for floats '
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
!!DB Allocate arrays
         ALLOCATE(iil(jpnfl));ALLOCATE(ijl(jpnfl));ALLOCATE(ikl(jpnfl))
         ALLOCATE(iiloc(jpnfl));ALLOCATE(ijloc(jpnfl))
         ALLOCATE(iiinfl(jpnfl));ALLOCATE(ijinfl(jpnfl));ALLOCATE(ikinfl(jpnfl))
         ALLOCATE(iioutfl(jpnfl));ALLOCATE(ijoutfl(jpnfl));ALLOCATE(ikoutfl(jpnfl))

         ALLOCATE(zgifl(jpnfl));ALLOCATE(zgjfl(jpnfl));ALLOCATE(zgkfl(jpnfl))
         ALLOCATE(zgidfl(jpnfl));ALLOCATE(zgjdfl(jpnfl));ALLOCATE(zgkdfl(jpnfl))
         ALLOCATE(ztxfl(jpnfl));ALLOCATE(ztyfl(jpnfl));ALLOCATE(ztzfl(jpnfl))
         ALLOCATE(zudfl(jpnfl));ALLOCATE(zvdfl(jpnfl));ALLOCATE(zwdfl(jpnfl))
         ALLOCATE(zufl(jpnfl));ALLOCATE(zvfl(jpnfl));ALLOCATE(zwfl(jpnfl))
         ALLOCATE(zttfl(jpnfl));ALLOCATE(zagefl(jpnfl));ALLOCATE(zagenewfl(jpnfl))
      ENDIF


      ! Initialisation of parameters
      
      DO jfl = 1, jpnfl
         ! ages of floats are put at zero
         zagefl(jfl) = 0.
         ! index on the velocity grid 
         ! We considere k coordinate negative, with this transformation 
         ! the computation in the 3 direction is the same. 
         zgifl(jfl) = tpifl(jfl) - 0.5
         zgjfl(jfl) = tpjfl(jfl) - 0.5
         zgkfl(jfl) = MIN(-1.,-(tpkfl(jfl)))
         ! surface drift every 10 days 
         IF( ln_argo ) THEN
            IF( MOD(kt,150) >= 146 .OR. MOD(kt,150) == 0 )  zgkfl(jfl) = -1.
         ENDIF
         ! index of T mesh
         iil(jfl) = 1 + INT(zgifl(jfl))
         ijl(jfl) = 1 + INT(zgjfl(jfl))
         ikl(jfl) =     INT(zgkfl(jfl))
      END DO
       
      iloop = 0
222   DO jfl = 1, jpnfl
# if   defined key_mpp_mpi   ||   defined key_mpp_shmem
         IF( (iil(jfl) >= (mig(nldi)-jpizoom+1)) .AND. (iil(jfl) <= (mig(nlei)-jpizoom+1)) .AND.   &
             (ijl(jfl) >= (mjg(nldj)-jpjzoom+1)) .AND. (ijl(jfl) <= (mjg(nlej)-jpjzoom+1)) ) THEN
            iiloc(jfl) = iil(jfl) - (mig(1)-jpizoom+1) + 1
            ijloc(jfl) = ijl(jfl) - (mjg(1)-jpjzoom+1) + 1
# else 
            iiloc(jfl) = iil(jfl)
            ijloc(jfl) = ijl(jfl)
# endif
            
            ! compute the transport across the mesh where the float is.            
            zsurfx(1) = e2u(iiloc(jfl)-1,ijloc(jfl)  ) * e3t(-ikl(jfl))
            zsurfx(2) = e2u(iiloc(jfl)  ,ijloc(jfl)  ) * e3t(-ikl(jfl))
            zsurfy(1) = e1v(iiloc(jfl)  ,ijloc(jfl)-1) * e3t(-ikl(jfl))
            zsurfy(2) = e1v(iiloc(jfl)  ,ijloc(jfl)  ) * e3t(-ikl(jfl))

            ! for a isobar float zsurfz is put to zero. The vertical velocity will be zero too.
            zsurfz=  e1t(iiloc(jfl),ijloc(jfl)) * e2t(iiloc(jfl),ijloc(jfl))
            zvol  =( e1t(iiloc(jfl),ijloc(jfl)) * e2t(iiloc(jfl),ijloc(jfl)) * e3t(-ikl(jfl)) )

            !
            zuinfl =( ub(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) )/2.*zsurfx(1)
            zuoutfl=( ub(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) )/2.*zsurfx(2)
            zvinfl =( vb(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) )/2.*zsurfy(1)
            zvoutfl=( vb(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) )/2.*zsurfy(2)
            zwinfl =-(wb(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1))    &
               &   +  wn(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1)) )/2. *  zsurfz*nisobfl(jfl)
            zwoutfl=-(wb(iiloc(jfl),ijloc(jfl),- ikl(jfl)   )   &
               &   +  wn(iiloc(jfl),ijloc(jfl),- ikl(jfl)   ) )/2. *  zsurfz*nisobfl(jfl)
            
            ! interpolation of velocity field on the float initial position            
            zufl(jfl)=  zuinfl  + ( zgifl(jfl) - float(iil(jfl)-1) ) * ( zuoutfl - zuinfl)
            zvfl(jfl)=  zvinfl  + ( zgjfl(jfl) - float(ijl(jfl)-1) ) * ( zvoutfl - zvinfl)
            zwfl(jfl)=  zwinfl  + ( zgkfl(jfl) - float(ikl(jfl)-1) ) * ( zwoutfl - zwinfl)
            
            ! faces of input and output
            ! u-direction
            IF( zufl(jfl) < 0. ) THEN
               iioutfl(jfl) = iil(jfl) - 1.
               iiinfl (jfl) = iil(jfl)
               zind   = zuinfl
               zuinfl = zuoutfl
               zuoutfl= zind
            ELSE
               iioutfl(jfl) = iil(jfl)
               iiinfl (jfl) = iil(jfl) - 1
            ENDIF
            ! v-direction       
            IF( zvfl(jfl) < 0. ) THEN
               ijoutfl(jfl) = ijl(jfl) - 1.
               ijinfl (jfl) = ijl(jfl)
               zind    = zvinfl
               zvinfl  = zvoutfl
               zvoutfl = zind
            ELSE
               ijoutfl(jfl) = ijl(jfl)
               ijinfl (jfl) = ijl(jfl) - 1.
            ENDIF
            ! w-direction
            IF( zwfl(jfl) < 0. ) THEN
               ikoutfl(jfl) = ikl(jfl) - 1.
               ikinfl (jfl) = ikl(jfl)
               zind    = zwinfl
               zwinfl  = zwoutfl
               zwoutfl = zind
            ELSE
               ikoutfl(jfl) = ikl(jfl)
               ikinfl (jfl) = ikl(jfl) - 1.
            ENDIF
            
            ! compute the time to go out the mesh across a face
            ! u-direction
            zudfl (jfl) = zuoutfl - zuinfl
            zgidfl(jfl) = float(iioutfl(jfl) - iiinfl(jfl))
            IF( zufl(jfl)*zuoutfl <= 0. ) THEN
               ztxfl(jfl) = 1.E99
            ELSE
               IF( ABS(zudfl(jfl)) >= 1.E-5 ) THEN
                  ztxfl(jfl)= zgidfl(jfl)/zudfl(jfl) * LOG(zuoutfl/zufl (jfl))
               ELSE
                  ztxfl(jfl)=(float(iioutfl(jfl))-zgifl(jfl))/zufl(jfl)
               ENDIF
               IF( (ABS(zgifl(jfl)-float(iiinfl (jfl))) <=  1.E-7) .OR.   &
                   (ABS(zgifl(jfl)-float(iioutfl(jfl))) <=  1.E-7) ) THEN
                  ztxfl(jfl)=(zgidfl(jfl))/zufl(jfl)
               ENDIF
            ENDIF
            ! v-direction
            zvdfl (jfl) = zvoutfl - zvinfl
            zgjdfl(jfl) = float(ijoutfl(jfl)-ijinfl(jfl))
            IF( zvfl(jfl)*zvoutfl <= 0. ) THEN
               ztyfl(jfl) = 1.E99
            ELSE
               IF( ABS(zvdfl(jfl)) >= 1.E-5 ) THEN
                  ztyfl(jfl) = zgjdfl(jfl)/zvdfl(jfl) * LOG(zvoutfl/zvfl (jfl))
               ELSE
                  ztyfl(jfl) = (float(ijoutfl(jfl)) - zgjfl(jfl))/zvfl(jfl)
               ENDIF
               IF( (ABS(zgjfl(jfl)-float(ijinfl (jfl))) <= 1.E-7) .OR.   &
                   (ABS(zgjfl(jfl)-float(ijoutfl(jfl))) <=  1.E-7) ) THEN
                  ztyfl(jfl) = (zgjdfl(jfl)) / zvfl(jfl)
               ENDIF
            ENDIF
            ! w-direction        
            IF( nisobfl(jfl) == 1. ) THEN 
               zwdfl (jfl) = zwoutfl - zwinfl
               zgkdfl(jfl) = float(ikoutfl(jfl) - ikinfl(jfl))
               IF( zwfl(jfl)*zwoutfl <= 0. ) THEN
                  ztzfl(jfl) = 1.E99
               ELSE
                  IF( ABS(zwdfl(jfl)) >= 1.E-5 ) THEN
                     ztzfl(jfl) = zgkdfl(jfl)/zwdfl(jfl) * LOG(zwoutfl/zwfl (jfl))
                  ELSE
                     ztzfl(jfl) = (float(ikoutfl(jfl)) - zgkfl(jfl))/zwfl(jfl)
                  ENDIF
                  IF( (ABS(zgkfl(jfl)-float(ikinfl (jfl))) <=  1.E-7) .OR.   &
                      (ABS(zgkfl(jfl)-float(ikoutfl(jfl))) <= 1.E-7) ) THEN
                     ztzfl(jfl) = (zgkdfl(jfl)) / zwfl(jfl)
                  ENDIF
               ENDIF
            ENDIF
            
            ! the time to go leave the mesh is the smallest time
                   
            IF( nisobfl(jfl) == 1. ) THEN 
               zttfl(jfl) = MIN(ztxfl(jfl),ztyfl(jfl),ztzfl(jfl))
            ELSE
               zttfl(jfl) = MIN(ztxfl(jfl),ztyfl(jfl))
            ENDIF
            ! new age of the FLOAT
            zagenewfl(jfl) = zagefl(jfl) + zttfl(jfl)*zvol
            ! test to know if the "age" of the float is not bigger than the 
            ! time step
            IF( zagenewfl(jfl) > rdt ) THEN
               zttfl(jfl) = (rdt-zagefl(jfl)) / zvol
               zagenewfl(jfl) = rdt
            ENDIF
            
            ! In the "minimal" direction we compute the index of new mesh
            ! on i-direction
            IF( ztxfl(jfl) <=  zttfl(jfl) ) THEN
               zgifl(jfl) = float(iioutfl(jfl))
               ind = iioutfl(jfl)
               IF( iioutfl(jfl) >= iiinfl(jfl) ) THEN
                  iioutfl(jfl) = iioutfl(jfl) + 1
               ELSE
                  iioutfl(jfl) = iioutfl(jfl) - 1
               ENDIF
               iiinfl(jfl) = ind
            ELSE
               IF( ABS(zudfl(jfl)) >= 1.E-5 ) THEN 
                  zgifl(jfl) = zgifl(jfl) + zgidfl(jfl)*zufl(jfl)    &
                     &       * ( EXP( zudfl(jfl)/zgidfl(jfl)*zttfl(jfl) ) - 1. ) /  zudfl(jfl)
               ELSE
                  zgifl(jfl) = zgifl(jfl) + zufl(jfl) * zttfl(jfl)
               ENDIF
            ENDIF
            ! on j-direction
            IF( ztyfl(jfl) <= zttfl(jfl) ) THEN
               zgjfl(jfl) = float(ijoutfl(jfl))
               ind = ijoutfl(jfl)
               IF( ijoutfl(jfl) >= ijinfl(jfl) ) THEN
                  ijoutfl(jfl) = ijoutfl(jfl) + 1
               ELSE
                  ijoutfl(jfl) = ijoutfl(jfl) - 1
               ENDIF
               ijinfl(jfl) = ind
            ELSE
               IF( ABS(zvdfl(jfl)) >= 1.E-5 ) THEN 
                  zgjfl(jfl) = zgjfl(jfl)+zgjdfl(jfl)*zvfl(jfl)   &
                     &       * ( EXP(zvdfl(jfl)/zgjdfl(jfl)*zttfl(jfl)) - 1. ) /  zvdfl(jfl)
               ELSE
                  zgjfl(jfl) = zgjfl(jfl)+zvfl(jfl)*zttfl(jfl)
               ENDIF
            ENDIF
            ! on k-direction
            IF( nisobfl(jfl) == 1. ) THEN
               IF( ztzfl(jfl) <= zttfl(jfl) ) THEN
                  zgkfl(jfl) = float(ikoutfl(jfl))
                  ind = ikoutfl(jfl)
                  IF( ikoutfl(jfl) >= ikinfl(jfl) ) THEN
                     ikoutfl(jfl) = ikoutfl(jfl)+1
                  ELSE
                     ikoutfl(jfl) = ikoutfl(jfl)-1
                  ENDIF
                  ikinfl(jfl) = ind
               ELSE
                  IF( ABS(zwdfl(jfl)) >= 1.E-5 ) THEN 
                     zgkfl(jfl) = zgkfl(jfl)+zgkdfl(jfl)*zwfl(jfl)    &
                        &       * ( EXP(zwdfl(jfl)/zgkdfl(jfl)*zttfl(jfl)) - 1. ) /  zwdfl(jfl)
                  ELSE
                     zgkfl(jfl) = zgkfl(jfl)+zwfl(jfl)*zttfl(jfl)
                  ENDIF
               ENDIF
            ENDIF
            
            ! coordinate of the new point on the temperature grid
            
            iil(jfl) = MAX(iiinfl(jfl),iioutfl(jfl))
            ijl(jfl) = MAX(ijinfl(jfl),ijoutfl(jfl))
            IF( nisobfl(jfl) ==  1 ) ikl(jfl) = MAX(ikinfl(jfl),ikoutfl(jfl))
            ! reinitialisation of the age of FLOAT
            zagefl(jfl) = zagenewfl(jfl)
# if   defined key_mpp_mpi   ||   defined key_mpp_shmem
         ELSE
            ! we put zgifl, zgjfl, zgkfl, zagefl
            zgifl (jfl) = 0.
            zgjfl (jfl) = 0.
            zgkfl (jfl) = 0.
            zagefl(jfl) = 0.
            iil(jfl) = 0
            ijl(jfl) = 0
         ENDIF
# endif
      END DO
      
      ! synchronisation
      IF( lk_mpp )   CALL mpp_sum( zgifl , jpnfl )   ! sums over the global domain
      IF( lk_mpp )   CALL mpp_sum( zgjfl , jpnfl )
      IF( lk_mpp )   CALL mpp_sum( zgkfl , jpnfl )
      IF( lk_mpp )   CALL mpp_sum( zagefl, jpnfl )
      IF( lk_mpp )   CALL mpp_sum( iil   , jpnfl )
      IF( lk_mpp )   CALL mpp_sum( ijl   , jpnfl )
      
      ! in the case of open boundaries we need to test if the floats don't
      ! go out of the domain. If it goes out, the float is put at the 
      ! middle of the mesh in the domain but the trajectory isn't compute 
      ! more time.      
# if defined key_obc
!!DB - 2008.03.10 -- errors in variable names if free-surface ----> fix 
      DO jfl = 1, jpnfl
         IF( lp_obc_east ) THEN
!            IF( jped <=  zgjfl(jfl) .AND. zgjfl(jfl) <= jpef .AND. nieob-1 <=  zgifl(jfl) ) THEN
            IF( jpjed <=  zgjfl(jfl) .AND. zgjfl(jfl) <= jpjef .AND. jpieob-1 <=  zgifl(jfl) ) THEN
               zgifl (jfl) = INT(zgifl(jfl)) + 0.5
               zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
               zagefl(jfl) = rdt
            END IF
         END IF
         IF( lp_obc_west ) THEN
!            IF( jpwd <= zgjfl(jfl) .AND. zgjfl(jfl) <= jpwf .AND. niwob >=  zgifl(jfl) ) THEN
            IF( jpjwd <= zgjfl(jfl) .AND. zgjfl(jfl) <= jpjwf .AND. jpiwob >=  zgifl(jfl) ) THEN
               zgifl (jfl) = INT(zgifl(jfl)) + 0.5
               zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
               zagefl(jfl) = rdt
            END IF
         END IF
         IF( lp_obc_north ) THEN
!            IF( jpnd <=  zgifl(jfl) .AND. zgifl(jfl) <= jpnf .AND. njnob-1 >=  zgjfl(jfl) ) THEN
            IF( jpind <=  zgifl(jfl) .AND. zgifl(jfl) <= jpinf .AND. jpjnob-1 >=  zgjfl(jfl) ) THEN
               zgifl (jfl) = INT(zgifl(jfl)) + 0.5
               zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
               zagefl(jfl) = rdt
            END IF
         END IF
         IF( lp_obc_south ) THEN
!            IF( jpsd <=  zgifl(jfl) .AND. zgifl(jfl) <= jpsf .AND.  njsob >= zgjfl(jfl) ) THEN
            IF( jpisd <=  zgifl(jfl) .AND. zgifl(jfl) <= jpisf .AND.  jpjsob >= zgjfl(jfl) ) THEN
               zgifl (jfl) = INT(zgifl(jfl)) + 0.5
               zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
               zagefl(jfl) = rdt
            END IF
         END IF
      END DO
#endif

      ! Test to know if a  float hasn't integrated enought time
      IF( ln_argo ) THEN
         ifin = 1
         DO jfl = 1, jpnfl
            IF( zagefl(jfl) < rdt )   ifin = 0
            tpifl(jfl) = zgifl(jfl) + 0.5
            tpjfl(jfl) = zgjfl(jfl) + 0.5
         END DO
      ELSE
         ifin = 1
         DO jfl = 1, jpnfl
            IF( zagefl(jfl) < rdt )   ifin = 0
            tpifl(jfl) = zgifl(jfl) + 0.5
            tpjfl(jfl) = zgjfl(jfl) + 0.5
            IF( nisobfl(jfl) == 1 ) tpkfl(jfl) = -(zgkfl(jfl))
         END DO
      ENDIF
      IF( ifin == 0 ) THEN
         iloop = iloop + 1 
         GO TO 222
      ENDIF

   END SUBROUTINE flo_blk

!!DB
   SUBROUTINE flo_RDM( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_RDM  ***
      !!           
      !! ** Purpose :   Compute the geographical position,latitude, longitude
      !!      and depth of each float at each time step.
      !! 
      !! ** Method  :   The position of a float is updated using a simple Euler 
      !!      algorithm. We need to know the velocity field, the old positions
      !!      of the floats and the grid defined on the domain.
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt ! ocean time step

      !! * Local declarations
      INTEGER :: jfl              ! dummy loop arguments
      INTEGER :: ind, ifin, iloop
!!DB
!      INTEGER , DIMENSION ( jpnfl )  ::   &
      INTEGER , DIMENSION (:),ALLOCATABLE, SAVE  ::   &
         iil, ijl, ikl,             &     ! index of nearest mesh
         iiloc , ijloc 
!!DB
!      REAL(wp) , DIMENSION ( jpnfl )  ::    &
      REAL(wp) , DIMENSION (:),ALLOCATABLE, SAVE  ::    &
         zgifl, zgjfl, zgkfl,       &     ! position of floats, index on velocity mesh.
         zufl, zvfl, zwfl
      REAL(wp)   ::       &
         zuinfl,zvinfl,zwinfl,      &     ! transport across the input face
         zuoutfl,zvoutfl,zwoutfl
!!DB
      REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: dxfl, dyfl, dzfl, dxifl, dyjfl, dzkfl
      REAL(wp) ::   &
           val, hor_diff, dx_rdm, dy_rdm, Ah_lo, Ah_hi, dAh_dx, dAh_dy, &
           Ahx_top, Ahx_btm, Ahy_left, Ahy_right, Ah_bl, Ah_br, Ah_tl, Ah_tr
      integer :: iseed, i1, i2, i3, j1, j2, k3
      LOGICAL :: RDM = .TRUE.

      !!---------------------------------------------------------------------
      !!---------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'DBG: flo_RDM : RDM trajectories for floats '
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
!!DB
!!NB: iseed = 5 means that the same random numbers will be generated every run
!!    For a different series must "randomize" the assignment of iseed
         iseed = 5
         call sgrnd(iseed)   !!initialize Mersenne Twister
         val = gasdev(iseed)
         if(lwp) write(numout,*)'DBG: kt, gasdev = ', kt, val
!!DB Allocate arrays
         ALLOCATE(iil(jpnfl));ALLOCATE(ijl(jpnfl));ALLOCATE(ikl(jpnfl))
         ALLOCATE(iiloc(jpnfl));ALLOCATE(ijloc(jpnfl))

         ALLOCATE(zgifl(jpnfl));ALLOCATE(zgjfl(jpnfl));ALLOCATE(zgkfl(jpnfl))
         ALLOCATE(zufl(jpnfl));ALLOCATE(zvfl(jpnfl));ALLOCATE(zwfl(jpnfl))

         ALLOCATE(dxfl(jpnfl)); ALLOCATE(dyfl(jpnfl)); ALLOCATE(dzfl(jpnfl))
         ALLOCATE(dxifl(jpnfl)); ALLOCATE(dyjfl(jpnfl)); ALLOCATE(dzkfl(jpnfl))

      ENDIF

      ! Initialisation of parameters
      
      do jfl = 1, jpnfl
         ! index on the velocity grid 
         ! We considere k coordinate negative, with this transformation 
         ! the computation in the 3 direction is the same. 
         zgifl(jfl) = tpifl(jfl) - 0.5
         zgjfl(jfl) = tpjfl(jfl) - 0.5
         zgkfl(jfl) = MIN(-1.,-(tpkfl(jfl)))
         ! index of T mesh
         iil(jfl) = 1 + INT(zgifl(jfl))
         ijl(jfl) = 1 + INT(zgjfl(jfl))
         ikl(jfl) =     INT(zgkfl(jfl))
      enddo
       
      do jfl = 1, jpnfl
# if   defined key_mpp_mpi   ||   defined key_mpp_shmem
         IF( (iil(jfl) >= (mig(nldi)-jpizoom+1)) .AND. (iil(jfl) <= (mig(nlei)-jpizoom+1)) .AND.   &
             (ijl(jfl) >= (mjg(nldj)-jpjzoom+1)) .AND. (ijl(jfl) <= (mjg(nlej)-jpjzoom+1)) ) THEN
            iiloc(jfl) = iil(jfl) - (mig(1)-jpizoom+1) + 1
            ijloc(jfl) = ijl(jfl) - (mjg(1)-jpjzoom+1) + 1
# else 
            iiloc(jfl) = iil(jfl)
            ijloc(jfl) = ijl(jfl)
# endif
            
            ! compute the transport across the mesh where the float is.            
!            zsurfx(1) = e2u(iiloc(jfl)-1,ijloc(jfl)  ) * e3t(-ikl(jfl))
!            zsurfx(2) = e2u(iiloc(jfl)  ,ijloc(jfl)  ) * e3t(-ikl(jfl))
!            zsurfy(1) = e1v(iiloc(jfl)  ,ijloc(jfl)-1) * e3t(-ikl(jfl))
!            zsurfy(2) = e1v(iiloc(jfl)  ,ijloc(jfl)  ) * e3t(-ikl(jfl))

            ! for a isobar float zsurfz is put to zero. The vertical velocity will be zero too.
!            zsurfz=  e1t(iiloc(jfl),ijloc(jfl)) * e2t(iiloc(jfl),ijloc(jfl))
!            zvol  =( e1t(iiloc(jfl),ijloc(jfl)) * e2t(iiloc(jfl),ijloc(jfl)) * e3t(-ikl(jfl)) )


!!DB  Modify to compute velocity only
            !
            zuinfl =( ub(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) )/2.
            zuoutfl=( ub(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) )/2.
            zvinfl =( vb(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) )/2.
            zvoutfl=( vb(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) )/2.
            zwinfl =-(wb(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1))    &
               &   +  wn(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1)) )/2. * nisobfl(jfl)
            zwoutfl=-(wb(iiloc(jfl),ijloc(jfl),- ikl(jfl)   )   &
               &   +  wn(iiloc(jfl),ijloc(jfl),- ikl(jfl)   ) )/2. * nisobfl(jfl)
            
            ! interpolation of velocity field on the float initial position            
            zufl(jfl)=  zuinfl  + ( zgifl(jfl) - float(iil(jfl)-1) ) * ( zuoutfl - zuinfl)
            zvfl(jfl)=  zvinfl  + ( zgjfl(jfl) - float(ijl(jfl)-1) ) * ( zvoutfl - zvinfl)
            zwfl(jfl)=  zwinfl  + ( zgkfl(jfl) - float(ikl(jfl)-1) ) * ( zwoutfl - zwinfl)

!!DB update position
!!the zgifl index holds the global U index == fractional position on the vel-grid
!!tpifl holds the global T index == fractional position on the T-grid
!!These are displacements in meters (NB: could be -ve)
            dxfl(jfl) = zufl(jfl)*rdt
            dyfl(jfl) = zvfl(jfl)*rdt
! ignore z-direction for now
!            dzfl(jfl) = zwfl(jfl)*rdt


!!DB: correct Random Displacement Model
!! ie it uses the local Ah plus its spatial derivative
!!Algorithm: 
!!Part 1: Ah interpolated to float location; use F-pt diffusivity and bilinear interp
!! Displacement ~ N(0,Ah_local)
!!Part 2: Grad of Ah
!! Add correction to displacement to account for horizontal gradients in diffusivity
!!NB: z-direction ignored
            if(RDM) then

!!Code for tracer diffusivity
#ifdef  key_traldf_smag
               i1 = max(iiloc(jfl)-1,1)
               i2 = iiloc(jfl)
               i3 = min(iiloc(jfl)+1,jpi)
               j1 = max(ijloc(jfl)-1,1)
               j2 = ijloc(jfl)
               k3 = -ikl(jfl)
               Ah_tl = 0.5*(ahtv(i1,j2,k3)+ahtv(i2,j2,k3))
               Ah_tr = 0.5*(ahtv(i2,j2,k3)+ahtv(i3,j2,k3))
               Ah_br = 0.5*(ahtv(i2,j1,k3)+ahtv(i3,j1,k3))
               Ah_bl = 0.5*(ahtv(i1,j1,k3)+ahtv(i2,j1,k3))
#else
!Use momentum diffusivity 
!-0- Ah at 4 corners of cell: tr=top-right; bl=btm-left, etc
               Ah_tl = ahm2(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl))
               Ah_tr = ahm2(iiloc(jfl),ijloc(jfl),-ikl(jfl))
               Ah_br = ahm2(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl))
               Ah_bl = ahm2(iiloc(jfl)-1,ijloc(jfl)-1,-ikl(jfl))
#endif

!!Part 1:
!-1- Ah at x-position at top of cell
               Ahx_top = Ah_tl + (zgifl(jfl) - float(iil(jfl)-1))*(Ah_tr-Ah_tl) 
!-2- Ah at x-position at btm of cell
               Ahx_btm = Ah_bl + (zgifl(jfl) - float(iil(jfl)-1))*(Ah_br-Ah_bl) 
!-3- Interp to y-position
               hor_diff = Ahx_btm + (zgjfl(jfl) - float(ijl(jfl)-1))*(Ahx_top-Ahx_btm)
               dx_rdm = sqrt(2.0*hor_diff*rdt)*gasdev(iseed)
               dy_rdm = sqrt(2.0*hor_diff*rdt)*gasdev(iseed)
               dxfl(jfl) = dxfl(jfl) + dx_rdm
               dyfl(jfl) = dyfl(jfl) + dy_rdm

!!Part 2: Compute drift correction term ~ d(Ah)/dx * rdt
!!DB: need local gradients of Ah, which is more complicated as local Ah indices +/- 2 are needed
!!which can overindex local domain. A solution is to compute grad(Ah) everywhere, call lbc_lnk,
!!and then interpolate the halo-updated field. This is too computationally expensive, so I will
!!use a simpler method. 
!!Note that the bilinear interpolation used above is equivalent to assuming constant gradients
!!within the cell ====>
               dAh_dy = (Ahx_top-Ahx_btm)/e2v(iiloc(jfl),ijloc(jfl)-1)
               Ahy_left = Ah_bl + (zgjfl(jfl) - float(ijl(jfl)-1))*(Ah_tl-Ah_bl)
               Ahy_right = Ah_br + (zgjfl(jfl) - float(ijl(jfl)-1))*(Ah_tr-Ah_br)
               dAh_dx = (Ahy_right-Ahy_left)/e1u(iiloc(jfl)-1,ijloc(jfl))

               hor_diff = dAh_dx * rdt
               dxfl(jfl) = dxfl(jfl) + hor_diff
               hor_diff = dAh_dy * rdt
               dyfl(jfl) = dyfl(jfl) + hor_diff

            endif   !!RDM true

!displacement relative to cell size
            dxifl(jfl) = dxfl(jfl)/e1u(iiloc(jfl)-1,ijloc(jfl))
            dyjfl(jfl) = dyfl(jfl)/e2v(iiloc(jfl),ijloc(jfl)-1)



!!DB: 04.04: I'll ignore 2 possibilities (for now):
!-1- If particle leaves cell then relative displacement should be adjusted to account
!for part of its travel in the new cell where e1u may be different. As this is likely a
!small difference, I'll ignore it for now and fix it at some future time (likely never!)
!-2- displacement is > cell size; ie dx = U*rdt > e1u: This is like a CFL condition
!which I'll assume is satisfied as dt = rdt (ie routine called every kt). If the routine
!changes so that it is called every N*rdt seconds then we would have to check if displacement
!is greater than 1 cell size and adjust as necessary.

            zgifl(jfl) = zgifl(jfl) + dxifl(jfl)
            zgjfl(jfl) = zgjfl(jfl) + dyjfl(jfl)
! ignore w-direction for now


# if   defined key_mpp_mpi   ||   defined key_mpp_shmem
         ELSE
            ! we put zgifl, zgjfl, zgkfl
            zgifl (jfl) = 0.
            zgjfl (jfl) = 0.
            zgkfl (jfl) = 0.
            iil(jfl) = 0
            ijl(jfl) = 0
         ENDIF
# endif
      END DO
      
      ! synchronisation
      IF( lk_mpp )   CALL mpp_sum( zgifl , jpnfl )   ! sums over the global domain
      IF( lk_mpp )   CALL mpp_sum( zgjfl , jpnfl )
      IF( lk_mpp )   CALL mpp_sum( zgkfl , jpnfl )
      IF( lk_mpp )   CALL mpp_sum( iil   , jpnfl )
      IF( lk_mpp )   CALL mpp_sum( ijl   , jpnfl )
      
!!DB: Make sure float is inside domain. 
!!At this point all processors hold the same zg*fl variables, in global indices 
!!Cavalier at this stage but DBG 
      do jfl = 1, jpnfl      
         if(zgifl(jfl) > jpidta-1) zgifl(jfl) = jpidta-1
         if(zgjfl(jfl) > jpjdta-1) zgjfl(jfl) = jpjdta-1
         if(zgifl(jfl) < 1) zgifl(jfl) = 1
         if(zgjfl(jfl) < 1) zgjfl(jfl) = 1
      enddo


!!DB: the below (comment and code) is old and I do not think it is necessary 
!!I think what I did above is OK
      ! in the case of open boundaries we need to test if the floats don't
      ! go out of the domain. If it goes out, the float is put at the 
      ! middle of the mesh in the domain 
# if defined key_obc
      DO jfl = 1, jpnfl
         IF( lp_obc_east ) THEN
            IF( jpjed <=  zgjfl(jfl) .AND. zgjfl(jfl) <= jpjef .AND. jpieob-1 <=  zgifl(jfl) ) THEN
               zgifl (jfl) = INT(zgifl(jfl)) + 0.5
               zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
            END IF
         END IF
         IF( lp_obc_west ) THEN
            IF( jpjwd <= zgjfl(jfl) .AND. zgjfl(jfl) <= jpjwf .AND. jpiwob >=  zgifl(jfl) ) THEN
               zgifl (jfl) = INT(zgifl(jfl)) + 0.5
               zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
            END IF
         END IF
         IF( lp_obc_north ) THEN
            IF( jpind <=  zgifl(jfl) .AND. zgifl(jfl) <= jpinf .AND. jpjnob-1 >=  zgjfl(jfl) ) THEN
               zgifl (jfl) = INT(zgifl(jfl)) + 0.5
               zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
            END IF
         END IF
         IF( lp_obc_south ) THEN
            IF( jpisd <=  zgifl(jfl) .AND. zgifl(jfl) <= jpisf .AND.  jpjsob >= zgjfl(jfl) ) THEN
               zgifl (jfl) = INT(zgifl(jfl)) + 0.5
               zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
            END IF
         END IF
      END DO
#endif

!!DB
      do jfl = 1, jpnfl
         tpifl(jfl) = zgifl(jfl) + 0.5
         tpjfl(jfl) = zgjfl(jfl) + 0.5
!!Not sure what to do here
!         IF( nisobfl(jfl) == 1 ) tpkfl(jfl) = -(zgkfl(jfl))
      enddo



   END SUBROUTINE flo_RDM


!-----------------------------------------------------------------------
!
      FUNCTION GASDEV(IDUM)
!
! GASDEV - GASDEV RETURNS A NORMALLY DISTRIBUTED DEVIATE
!          WITH ZERO MEAN AND UNIT VARIAN!E, USING RAN1(IDUM) AS
!          THE SOURCE FOR THE UNIFORM DEVIATES.
!
      INTEGER :: IDUM
      REAL :: GASDEV
      REAL(8) :: n1, n2
      INTEGER :: ISET
      REAL :: FAC,GSET,RSQ,V1,V2,RAN1,RAN2
      !CN: no longer needed since grnd is a subroutine
      !double precision :: grnd
      SAVE ISET,GSET
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
        !CN: Changed grnd to subroutine to fix linker errors
        CALL grnd(n1)
        CALL grnd(n2)
1       V1=2.*n1-1.
        V2=2.*n2-1.
        !V1=2.*grnd()-1.
        !V2=2.*grnd()-1. 
        RSQ=V1**2+V2**2
        IF(RSQ.GE.1..OR.RSQ.EQ.0.)GOTO 1
        FAC=SQRT(-2.*LOG(RSQ)/RSQ)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN

      end function GASDEV



!  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE &OL`.
!

!     MERSENNE TWISTER CODE
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999

!Initialization subroutine
  subroutine sgrnd(seed)
    implicit none
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
    integer, intent(in) :: seed

    mt(0) = iand(seed,-1)
    do mti=1,N-1
      mt(mti) = iand(69069 * mt(mti-1),-1)
    enddo
!
    return
  end subroutine sgrnd

!Random number generator
!CN: Changed grnd to a subroutine to fix linker errors
  subroutine grnd(rnd_num)
    implicit integer(a-z)

    REAL(8),INTENT(OUT) :: rnd_num

! Period parameters
    integer, parameter :: M = 397, MATA  = -1727483681
!                                    constant vector a
    integer, parameter :: LMASK =  2147483647
!                                    least significant r bits
    integer, parameter :: UMASK = -LMASK - 1
!                                    most significant w-r bits
! Tempering parameters
    integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544

    dimension mag01(0:1)
    data mag01/0, MATA/
    save mag01
!                        mag01(x) = x * MATA for x=0,1

    TSHFTU(y)=ishft(y,-11)
    TSHFTS(y)=ishft(y,7)
    TSHFTT(y)=ishft(y,15)
    TSHFTL(y)=ishft(y,-18)

    if(mti.ge.N) then
!                       generate N words at one time
      if(mti.eq.N+1) then
!                            if sgrnd() has not been called,
        call sgrnd( defaultsd )
!                              a default initial seed is used
      endif

      do kk=0,N-M-1
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      do kk=N-M,N-2
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
      mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
      mti = 0
    endif

    y=mt(mti)
    mti = mti + 1 
    y=ieor(y,TSHFTU(y))
    y=ieor(y,iand(TSHFTS(y),TMASKB))
    y=ieor(y,iand(TSHFTT(y),TMASKC))
    y=ieor(y,TSHFTL(y))

    if(y .lt. 0) then
      rnd_num=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
    else
      rnd_num=dble(y)/(2.0d0**32-1.0d0)
    endif

    return
  end subroutine grnd

#  else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flo_blk                  ! Empty routine
   END SUBROUTINE flo_blk 
#endif
   
   !!======================================================================
END MODULE floblk 
