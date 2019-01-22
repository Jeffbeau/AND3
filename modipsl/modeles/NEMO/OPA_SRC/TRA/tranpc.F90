MODULE tranpc
   !!==============================================================================
   !!                       ***  MODULE  tranpc  ***
   !! Ocean active tracers:  non penetrative convection scheme
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_npc      : apply the non penetrative convection scheme
   !!   tra_npc_init : initialization and control of the scheme
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers 
   USE dom_oce         ! ocean space and time domain
   USE trdmod          ! ocean active tracer trends
   USE trdmod_oce      ! ocean variables trends
   USE eosbn2          ! equation of state (eos routine) 
   USE lbclnk          ! lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_npc      ! routine called by step.F90

   !! * Module variable
   INTEGER ::       &
      nnpc1 =   1,  &  ! nnpc1   non penetrative convective scheme frequency
      nnpc2 =  15      ! nnpc2   non penetrative convective scheme print frequency

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tranpc.F90,v 1.3 2005/03/27 18:35:20 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_npc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tranpc  ***
      !!
      !! ** Purpose :   Non penetrative convective adjustment scheme. solve 
      !!      the static instability of the water column (now, after the swap) 
      !!      while conserving heat and salt contents.
      !!
      !! ** Method  :   The algorithm used converges in a maximium of jpk 
      !!      iterations. instabilities are treated when the vertical density
      !!      gradient is less than 1.e-5.
      !!
      !!      'key_trdtra' defined: the trend associated with this
      !!                               algorithm is saved.
      !!
      !!      macro-tasked on vertical slab (jj-loop)
      !!
      !! ** Action  : - (tn,sn) after the application od the npc scheme
      !!              - save the associated trends (ttrd,strd) ('key_trdtra')
      !!
      !! References :
      !!      Madec, et al., 1991, JPO, 21, 9, 1349-1371.
      !!
      !! History :
      !!   1.0  !  90-09  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)  periodic conditions on t and s
      !!        !  93-03  (M. Guyon)  symetrical conditions 
      !!        !  96-01  (G. Madec)  statement function for e3
      !!                                  suppression of common work arrays
      !!   8.5  !  02-06  (G. Madec)  free form F90
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
                         ztdsa => va      ! use va as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk             ! dummy loop indices
      INTEGER ::   &
         inpcc ,                        & ! number of statically instable water column
         inpci ,                        & ! number of iteration for npc scheme
         jiter, jkdown, jkp,            & ! ???
         ikbot, ik, ikup, ikdown          ! ???
      REAL(wp) ::   &                     ! temporary arrays
         ze3tot, zta, zsa, zraua, ze3dwn
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         zwx, zwy, zwz                    ! temporary arrays
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zrhop                            ! temporary arrays
      !!----------------------------------------------------------------------

      IF( kt == nit000  )   CALL tra_npc_init


      IF( MOD( kt, nnpc1 ) == 0 ) THEN

         inpcc = 0
         inpci = 0

         ! 0. Potential density
         ! --------------------

         CALL eos( tn, sn, rhd, zrhop )

         ! Save tn and sn trends
         IF( l_trdtra )   THEN
            ztdta(:,:,:) = tn(:,:,:) 
            ztdsa(:,:,:) = sn(:,:,:) 
         ENDIF

         !                                                ! ===============
         DO jj = 1, jpj                                   !  Vertical slab
            !                                             ! ===============

            ! 1. Static instability pointer 
            ! -----------------------------

            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  zwx(ji,jk) = ( zrhop(ji,jj,jk) - zrhop(ji,jj,jk+1) ) * tmask(ji,jj,jk+1)
               END DO
            END DO

            ! 1.1 do not consider the boundary points

            ! even if east-west cyclic b. c. do not considere ji=1 or jpi
            DO jk = 1, jpkm1
               zwx( 1 ,jk) = 0.e0
               zwx(jpi,jk) = 0.e0
            END DO
            ! even if south-symmetric b. c. used, do not considere jj=1
            IF( jj == 1 ) zwx(:,:) = 0.e0

            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  zwx(ji,jk) = 1.
                  IF( zwx(ji,jk) < 1.e-5 ) zwx(ji,jk)=0.
               END DO
            END DO

            zwy(:,1) = 0.
            DO ji = 1, jpi
               DO jk = 1, jpkm1
                  zwy(ji,1) = zwy(ji,1) + zwx(ji,jk)
               END DO
            END DO

            zwz(1,1) = 0.
            DO ji = 1, jpi
               zwz(1,1) = zwz(1,1) + zwy(ji,1)
            END DO

            inpcc = inpcc + NINT( zwz(1,1) )


            ! 2. Vertical mixing for each instable portion of the density profil
            ! ------------------------------------------------------------------

            IF (zwz(1,1) /= 0.) THEN

               ! -->> the density profil is statically instable :

               DO ji = 1, jpi
                  IF( zwy(ji,1) /= 0. ) THEN

                     ! ikbot: ocean bottom level

                     ikbot = mbathy(ji,jj)

                     ! vertical iteration

                     DO jiter = 1, jpk

                        ! search of ikup : the first static instability from the sea surface

                        ik = 0
220                     CONTINUE
                        ik = ik + 1
                        IF( ik >= ikbot-1 ) GO TO 200
                        zwx(ji,ik) = zrhop(ji,jj,ik) - zrhop(ji,jj,ik+1)
                        IF( zwx(ji,ik) <= 0. ) GO TO 220
                        ikup = ik
                        ! the density profil is instable below ikup
   
                        ! ikdown : bottom of the instable portion of the density profil

                        ! search of ikdown and vertical mixing from ikup to ikdown

                        ze3tot= fse3t(ji,jj,ikup)
                        zta   = tn   (ji,jj,ikup)
                        zsa   = sn   (ji,jj,ikup)
                        zraua = zrhop(ji,jj,ikup)

                        DO jkdown = ikup+1, ikbot-1
                           IF( zraua <= zrhop(ji,jj,jkdown) ) THEN
                              ikdown = jkdown
                              GO TO 240
                           ENDIF
                           ze3dwn =  fse3t(ji,jj,jkdown)
                           ze3tot =  ze3tot + ze3dwn
                           zta   = ( zta*(ze3tot-ze3dwn) + tn(ji,jj,jkdown)*ze3dwn )/ze3tot
                           zsa   = ( zsa*(ze3tot-ze3dwn) + sn(ji,jj,jkdown)*ze3dwn )/ze3tot
                           zraua = ( zraua*(ze3tot-ze3dwn) + zrhop(ji,jj,jkdown)*ze3dwn )/ze3tot
                           inpci = inpci+1
                        END DO
                        ikdown = ikbot-1
240                     CONTINUE

                        DO jkp = ikup, ikdown-1
                           tn(ji,jj,jkp) = zta
                           sn(ji,jj,jkp) = zsa
                           zrhop(ji,jj,jkp) = zraua
                        END DO
                        IF (ikdown == ikbot-1 .AND. zraua >= zrhop(ji,jj,ikdown) ) THEN
                           tn(ji,jj,ikdown) = zta
                           sn(ji,jj,ikdown) = zsa
                           zrhop(ji,jj,ikdown) = zraua
                        ENDIF

                     END DO
                  ENDIF
200               CONTINUE
               END DO

               ! <<-- no more static instability on slab jj

            ENDIF
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


         ! save the trends for diagnostic
         ! Non penetrative mixing trends
         IF( l_trdtra )   THEN
            ztdta(:,:,:) = tn(:,:,:) - ztdta(:,:,:)
            ztdsa(:,:,:) = sn(:,:,:) - ztdsa(:,:,:)

            CALL trd_mod(ztdta, ztdsa, jpttdnpc, 'TRA', kt)
         ENDIF
      
         ! Lateral boundary conditions on ( tn, sn )   ( Unchanged sign)
         ! ------------------------------============
         CALL lbc_lnk( tn, 'T', 1. )
         CALL lbc_lnk( sn, 'T', 1. )
      

         !  2. non penetrative convective scheme statistics
         !  -----------------------------------------------

         IF( nnpc2 /= 0 .AND. MOD( kt, nnpc2 ) == 0 ) THEN
            IF(lwp) WRITE(numout,*)' kt=',kt, ' number of statically instable',   &
               ' water column : ',inpcc, ' number of iteration : ',inpci
         ENDIF

      ENDIF
      
   END SUBROUTINE tra_npc


   SUBROUTINE tra_npc_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_npc_init  ***
      !!                   
      !! ** Purpose :   initializations of the non-penetrative adjustment scheme
      !!
      !! History :
      !!   8.5  !  02-12  (G. Madec)  F90 : free form
      !!----------------------------------------------------------------------
      !! * Namelist
      NAMELIST/namnpc/ nnpc1, nnpc2
      !!----------------------------------------------------------------------

      ! Namelist namzdf : vertical diffusion
      REWIND( numnam )
      READ  ( numnam, namnpc )

      ! Parameter print
      ! ---------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_npc_init : Non Penetrative Convection (npc) scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist namnpc : set npc scheme parameters'
         WRITE(numout,*)
         WRITE(numout,*) '             npc scheme frequency           nnpc1  = ', nnpc1
         WRITE(numout,*) '             npc scheme print frequency     nnpc2  = ', nnpc2
         WRITE(numout,*)
      ENDIF


      ! Parameter controls
      ! ------------------
      IF ( nnpc1 == 0 ) THEN
          IF(lwp) WRITE(numout,cform_war)
          IF(lwp) WRITE(numout,*) '             nnpc1 = ', nnpc1, ' is forced to 1'
          nnpc1 = 1
          nwarn = nwarn + 1
      ENDIF
      
   END SUBROUTINE tra_npc_init

   !!======================================================================
END MODULE tranpc
