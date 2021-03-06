MODULE divcur
   !!==============================================================================
   !!                       ***  MODULE  divcur  ***
   !! Ocean diagnostic variable : horizontal divergence and relative vorticity
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   div_cur    : Compute the horizontal divergence and relative
   !!                vorticity fields
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE obc_oce         ! ocean lateral open boundary condition
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC div_cur    ! routine called by step.F90 and istate.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/divcur.F90,v 1.6 2006/03/10 10:55:41 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

#if defined key_noslip_accurate
   !!----------------------------------------------------------------------
   !!   'key_noslip_accurate'                     2nd order centered scheme
   !!                                                4th order at the coast
   !!----------------------------------------------------------------------

   SUBROUTINE div_cur( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE div_cur  ***
      !!
      !! ** Purpose :   compute the horizontal divergence and the relative
      !!      vorticity at before and now time-step
      !!
      !! ** Method  : 
      !!      I.  divergence :
      !!         - save the divergence computed at the previous time-step
      !!      (note that the Asselin filter has not been applied on hdivb)
      !!         - compute the now divergence given by :
      !!            * s-coordinate ('key_s_coord' defined)
      !!         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!         * z-coordinate (default key)
      !!         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
      !!         - apply lateral boundary conditions on hdivn 
      !!      II. vorticity :
      !!         - save the curl computed at the previous time-step
      !!            rotb = rotn
      !!      (note that the Asselin time filter has not been applied to rotb)
      !!         - compute the now curl in tensorial formalism:
      !!            rotn = 1/(e1f*e2f) ( di[e2v vn] - dj[e1u un] )
      !!         - apply lateral boundary conditions on rotn through a call
      !!      of lbc_lnk routine.
      !!         - Coastal boundary condition: 'key_noslip_accurate' defined,
      !!      the no-slip boundary condition is computed using Schchepetkin
      !!      and O'Brien (1996) scheme (i.e. 4th order at the coast).
      !!      For example, along east coast, the one-sided finite difference
      !!      approximation used for di[v] is:
      !!         di[e2v vn] =  1/(e1f*e2f)
      !!                    * ( (e2v vn)(i) + (e2v vn)(i-1) + (e2v vn)(i-2) )
      !!
      !! ** Action  : - update hdivb, hdivn, the before & now hor. divergence
      !!              - update rotb , rotn , the before & now rel. vorticity
      !!
      !! History :
      !!   8.2  !  00-03  (G. Madec)  no slip accurate
      !!   9.0  !  03-08  (G. Madec)  merged of cur and div, free form, F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk     ! dummy loop indices
      INTEGER ::   ii, ij, jl     ! temporary integer
      INTEGER ::   ijt, iju       ! temporary integer
      REAL(wp) ::   zdiv, zdju
      REAL(wp), DIMENSION(   jpi  ,1:jpj+2) ::   zwu   ! workspace
      REAL(wp), DIMENSION(-1:jpi+2,  jpj  ) ::   zwv   ! workspace
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_cur : horizontal velocity divergence and relative vorticity'
         IF(lwp) WRITE(numout,*) '~~~~~~~   NOT optimal for auto-tasking case'
      ENDIF

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         hdivb(:,:,jk) = hdivn(:,:,jk)    ! time swap of div arrays
         rotb (:,:,jk) = rotn (:,:,jk)    ! time swap of rot arrays

         !                                             ! --------
         ! Horizontal divergence                       !   div
         !                                             ! --------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               hdivn(ji,jj,jk) =   &
                  (  e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk) - e2u(ji-1,jj  )*fse3u(ji-1,jj  ,jk)  * un(ji-1,jj  ,jk)       &
                   + e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk) - e1v(ji  ,jj-1)*fse3v(ji  ,jj-1,jk)  * vn(ji  ,jj-1,jk)  )    &
                  / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
#else
               hdivn(ji,jj,jk) = (  e2u(ji,jj) * un(ji,jj,jk) - e2u(ji-1,jj  ) * un(ji-1,jj  ,jk)      &
                  &               + e1v(ji,jj) * vn(ji,jj,jk) - e1v(ji  ,jj-1) * vn(ji  ,jj-1,jk)  )   &
     &            / ( e1t(ji,jj) * e2t(ji,jj) )
#endif
            END DO
         END DO

#if defined key_obc
         ! open boundaries (div must be zero behind the open boundary)
         !  mpp remark: The zeroing of hdivn can probably be extended to 1->jpi/jpj for the correct row/column
         IF( lp_obc_east  )   hdivn(nie0p1:nie1p1,nje0  :nje1  ,jk) = 0.e0      ! east
         IF( lp_obc_west  )   hdivn(niw0  :niw1  ,njw0  :njw1  ,jk) = 0.e0      ! west
         IF( lp_obc_north )   hdivn(nin0  :nin1  ,njn0p1:njn1p1,jk) = 0.e0      ! north
         IF( lp_obc_south )   hdivn(nis0  :nis1  ,njs0  :njs1  ,jk) = 0.e0      ! south
#endif         
#if defined key_agrif
         if ( .NOT. AGRIF_Root() ) then
            IF ((nbondi ==  1).OR.(nbondi == 2)) hdivn(nlci-1 , :     ,jk) = 0.e0      ! east
            IF ((nbondi == -1).OR.(nbondi == 2)) hdivn(2      , :     ,jk) = 0.e0      ! west
            IF ((nbondj ==  1).OR.(nbondj == 2)) hdivn(:      ,nlcj-1 ,jk) = 0.e0      ! north
            IF ((nbondj == -1).OR.(nbondj == 2)) hdivn(:      ,2      ,jk) = 0.e0      ! south
         endif
#endif       

         !                                             ! --------
         ! relative vorticity                          !   rot 
         !                                             ! --------
         ! contravariant velocity (extended for lateral b.c.)
         ! inside the model domain
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwu(ji,jj) = e1u(ji,jj) * un(ji,jj,jk)
               zwv(ji,jj) = e2v(ji,jj) * vn(ji,jj,jk)
            END DO  
         END DO  
 
         ! East-West boundary conditions
         IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6) THEN
            zwv(  0  ,:) = zwv(jpi-2,:)
            zwv( -1  ,:) = zwv(jpi-3,:)
            zwv(jpi+1,:) = zwv(  3  ,:)
            zwv(jpi+2,:) = zwv(  4  ,:)
         ELSE
            zwv(  0  ,:) = 0.e0
            zwv( -1  ,:) = 0.e0
            zwv(jpi+1,:) = 0.e0
            zwv(jpi+2,:) = 0.e0
         ENDIF

         ! North-South boundary conditions
         IF( nperio == 3 .OR. nperio == 4 ) THEN
            ! north fold ( Grid defined with a T-point pivot) ORCA 2 degre
            zwu(jpi,jpj+1) = 0.e0
            zwu(jpi,jpj+2) = 0.e0
            DO ji = 1, jpi-1
               iju = jpi - ji + 1
               zwu(ji,jpj+1) = - zwu(iju,jpj-3)
               zwu(ji,jpj+2) = - zwu(iju,jpj-4)
            END DO
         ELSEIF( nperio == 5 .OR. nperio == 6 ) THEN
            ! north fold ( Grid defined with a F-point pivot) ORCA 0.5 degre\
            zwu(jpi,jpj+1) = 0.e0
            zwu(jpi,jpj+2) = 0.e0
            DO ji = 1, jpi-1
               iju = jpi - ji
               zwu(ji,jpj  ) = - zwu(iju,jpj-1)
               zwu(ji,jpj+1) = - zwu(iju,jpj-2)
               zwu(ji,jpj+2) = - zwu(iju,jpj-3)
            END DO
            DO ji = -1, jpi+2
               ijt = jpi - ji + 1
               zwv(ji,jpj) = - zwv(ijt,jpj-2)
            END DO
            DO ji = jpi/2+1, jpi+2
               ijt = jpi - ji + 1
               zwv(ji,jpjm1) = - zwv(ijt,jpjm1)
            END DO
         ELSE
            ! closed
            zwu(:,jpj+1) = 0.e0
            zwu(:,jpj+2) = 0.e0
         ENDIF

         ! relative vorticity (vertical component of the velocity curl) 
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               rotn(ji,jj,jk) = (  zwv(ji+1,jj  ) - zwv(ji,jj)      &
                                 - zwu(ji  ,jj+1) + zwu(ji,jj)  )   &
                              * fmask(ji,jj,jk) / ( e1f(ji,jj)*e2f(ji,jj) )
            END DO
         END DO

         ! second order accurate scheme along straight coast
         DO jl = 1, npcoa(1,jk)
            ii = nicoa(jl,1,jk)
            ij = njcoa(jl,1,jk)
            rotn(ii,ij,jk) = 1. / ( e1f(ii,ij) * e2f(ii,ij) )   &
                           * ( + 4. * zwv(ii+1,ij) - zwv(ii+2,ij) + 0.2 * zwv(ii+3,ij) )
         END DO
         DO jl = 1, npcoa(2,jk)
            ii = nicoa(jl,2,jk)
            ij = njcoa(jl,2,jk)
            rotn(ii,ij,jk) = 1./(e1f(ii,ij)*e2f(ii,ij))   &
               *(-4.*zwv(ii,ij)+zwv(ii-1,ij)-0.2*zwv(ii-2,ij))
         END DO
         DO jl = 1, npcoa(3,jk)
            ii = nicoa(jl,3,jk)
            ij = njcoa(jl,3,jk)
            rotn(ii,ij,jk) = -1. / ( e1f(ii,ij)*e2f(ii,ij) )   &
               * ( +4. * zwu(ii,ij+1) - zwu(ii,ij+2) + 0.2 * zwu(ii,ij+3) )
         END DO
         DO jl = 1, npcoa(4,jk)
            ii = nicoa(jl,4,jk)
            ij = njcoa(jl,4,jk)
            rotn(ii,ij,jk) = -1. / ( e1f(ii,ij)*e2f(ii,ij) )   &
               * ( -4. * zwu(ii,ij) + zwu(ii,ij-1) - 0.2 * zwu(ii,ij-2) )
         END DO

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      
      ! 4. Lateral boundary conditions on hdivn and rotn
      ! ---------------------------------=======---======
      CALL lbc_lnk( hdivn, 'T', 1. )     ! T-point, no sign change
      CALL lbc_lnk( rotn , 'F', 1. )     ! F-point, no sign change

   END SUBROUTINE div_cur
   
#else
   !!----------------------------------------------------------------------
   !!   Default option                           2nd order centered schemes
   !!----------------------------------------------------------------------

   SUBROUTINE div_cur( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE div_cur  ***
      !!                    
      !! ** Purpose :   compute the horizontal divergence and the relative
      !!      vorticity at before and now time-step
      !!
      !! ** Method  : - Divergence:
      !!      - save the divergence computed at the previous time-step
      !!      (note that the Asselin filter has not been applied on hdivb)
      !!      - compute the now divergence given by :
      !!         * s-coordinate ('key_s_coord' defined)
      !!         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!         * z-coordinate (default key)
      !!         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
      !!      - apply lateral boundary conditions on hdivn 
      !!              - Relavtive Vorticity :
      !!      - save the curl computed at the previous time-step (rotb = rotn)
      !!      (note that the Asselin time filter has not been applied to rotb)
      !!      - compute the now curl in tensorial formalism:
      !!            rotn = 1/(e1f*e2f) ( di[e2v vn] - dj[e1u un] )
      !!      - apply lateral boundary conditions on rotn through a call to
      !!      routine lbc_lnk routine.
      !!      Note: Coastal boundary condition: lateral friction set through
      !!      the value of fmask along the coast (see dommsk.F90) and shlat
      !!      (namelist parameter)
      !!
      !! ** Action  : - update hdivb, hdivn, the before & now hor. divergence
      !!              - update rotb , rotn , the before & now rel. vorticity
      !!
      !! History :
      !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
      !!   4.0  !  91-11  (G. Madec)
      !!   6.0  !  93-03  (M. Guyon)  symetrical conditions
      !!   7.0  !  96-01  (G. Madec)  s-coordinates
      !!   8.0  !  97-06  (G. Madec)  lateral boundary cond., lbc
      !!   8.1  !  97-08  (J.M. Molines)  Open boundaries
      !!   9.0  !  02-09  (G. Madec, E. Durand)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
      
      !! * Local declarations
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_cur : horizontal velocity divergence and'
         IF(lwp) WRITE(numout,*) '~~~~~~~   relative vorticity'
      ENDIF

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         hdivb(:,:,jk) = hdivn(:,:,jk)    ! time swap of div arrays
         rotb (:,:,jk) = rotn (:,:,jk)    ! time swap of rot arrays

         !                                             ! --------
         ! Horizontal divergence                       !   div 
         !                                             ! --------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               hdivn(ji,jj,jk) =   &
                  (  e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk) - e2u(ji-1,jj  )*fse3u(ji-1,jj  ,jk)  * un(ji-1,jj  ,jk)       &
                   + e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk) - e1v(ji  ,jj-1)*fse3v(ji  ,jj-1,jk)  * vn(ji  ,jj-1,jk)  )    &
                  / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
#else
               hdivn(ji,jj,jk) = (  e2u(ji,jj) * un(ji,jj,jk) - e2u(ji-1,jj  ) * un(ji-1,jj  ,jk)      &
                  &               + e1v(ji,jj) * vn(ji,jj,jk) - e1v(ji  ,jj-1) * vn(ji  ,jj-1,jk)  )   & 
                  / ( e1t(ji,jj) * e2t(ji,jj) )
#endif
            END DO  
         END DO  

#if defined key_obc
         ! open boundaries (div must be zero behind the open boundary)
         !  mpp remark: The zeroing of hdivn can probably be extended to 1->jpi/jpj for the correct row/column
         IF( lp_obc_east  )   hdivn(nie0p1:nie1p1,nje0  :nje1  ,jk) = 0.e0      ! east
         IF( lp_obc_west  )   hdivn(niw0  :niw1  ,njw0  :njw1  ,jk) = 0.e0      ! west
         IF( lp_obc_north )   hdivn(nin0  :nin1  ,njn0p1:njn1p1,jk) = 0.e0      ! north
         IF( lp_obc_south )   hdivn(nis0  :nis1  ,njs0  :njs1  ,jk) = 0.e0      ! south
#endif         
#if defined key_agrif
         if ( .NOT. AGRIF_Root() ) then
            IF ((nbondi ==  1).OR.(nbondi == 2)) hdivn(nlci-1 , :     ,jk) = 0.e0      ! east
            IF ((nbondi == -1).OR.(nbondi == 2)) hdivn(2      , :     ,jk) = 0.e0      ! west
            IF ((nbondj ==  1).OR.(nbondj == 2)) hdivn(:      ,nlcj-1 ,jk) = 0.e0      ! north
            IF ((nbondj == -1).OR.(nbondj == 2)) hdivn(:      ,2      ,jk) = 0.e0      ! south
         endif
#endif       
         !                                             ! --------
         ! relative vorticity                          !   rot 
         !                                             ! --------
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               rotn(ji,jj,jk) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ,jk) - e2v(ji,jj) * vn(ji,jj,jk)    &
                  &              - e1u(ji  ,jj+1) * un(ji  ,jj+1,jk) + e1u(ji,jj) * un(ji,jj,jk)  ) &
                  &           * fmask(ji,jj,jk) / ( e1f(ji,jj) * e2f(ji,jj) )
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      
      ! 4. Lateral boundary conditions on hdivn and rotn
      ! ---------------------------------=======---======
      CALL lbc_lnk( hdivn, 'T', 1. )       ! T-point, no sign change
      CALL lbc_lnk( rotn , 'F', 1. )       ! F-point, no sign change

   END SUBROUTINE div_cur
   
#endif
   !!======================================================================
END MODULE divcur
