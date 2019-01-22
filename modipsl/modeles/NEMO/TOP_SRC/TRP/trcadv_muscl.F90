MODULE trcadv_muscl
   !!==============================================================================
   !!                       ***  MODULE  trcadv_muscl  ***
   !! Ocean passive tracers:  horizontal & vertical advective trend
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_adv_muscl : update the tracer trend with the horizontal
   !!                   and vertical advection trends using MUSCL scheme
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean passive tracers variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE trcbbl          ! advective passive tracers in the BBL
   USE lib_mpp
   USE prtctl_trc      ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_adv_muscl  ! routine called by trcstp.F90

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_muscl.F90,v 1.11 2006/04/10 15:38:54 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_adv_muscl( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE trc_adv_muscl  ***
      !!          
      !! ** Purpose :   Compute the now trend due to total advection of any pas-
      !!      sive tracer using a MUSCL scheme (Monotone Upstream-centered Scheme
      !!      for Conservation Laws) and add it to the general tracer trend.
      !!
      !! ** Method  :
      !!
      !! ** Action  : - update tra with the now advective tracer trends
      !!              - save trends in trtrd ('key_trc_diatrd')
      !!
      !! References :                
      !!      Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
      !!	IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
      !!
      !! History :
      !!        !  06-00  (A.Estublier)  for passive tracers
      !!   9.0  !  03-04  (C. Ethe, G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * modules used
#if defined key_trcbbl_adv
      USE oce_trc            , zun => ua,  &  ! use ua as workspace
         &                     zvn => va      ! use va as workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
#else
      USE oce_trc            , zun => un,  &  ! When no bbl, zun == un
                               zvn => vn,  &  !              zvn == vn
                               zwn => wn      !              zwn == wn
#endif

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      INTEGER ::   ji, jj, jk,jn            ! dummy loop indices
      REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
         zt1, zt2, ztp1, ztp2

      REAL(wp) ::   zu, zv, zw, zeu, zev, zew, zbtr, ztra
      REAL(wp) ::   z0u, z0v, z0w
      REAL(wp) ::   zzt1, zzt2, zalpha, z2dtt
#if defined key_trc_diatrd
      REAL(wp) ::   ztai, ztaj
      REAL(wp) ::   zfui, zfvj
#endif
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------


      IF( kt == nittrc000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_adv : MUSCL advection scheme'
         WRITE(numout,*) '~~~~~~~'
      ENDIF

 

#if defined key_trcbbl_adv
      ! Advective bottom boundary layer
      ! -------------------------------
      zun(:,:,:) = un (:,:,:) - u_trc_bbl(:,:,:)
      zvn(:,:,:) = vn (:,:,:) - v_trc_bbl(:,:,:)
      zwn(:,:,:) = wn (:,:,:) + w_trc_bbl(:,:,:)
#endif

 

      DO jn = 1, jptra

         ! I. Horizontal advective fluxes
         ! ------------------------------

         ! first guess of the slopes
         ! interior values
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1      
               DO ji = 1, fs_jpim1   ! vector opt.
                  zt1(ji,jj,jk) = umask(ji,jj,jk) * ( trb(ji+1,jj,jk,jn) - trb(ji,jj,jk,jn) )
                  zt2(ji,jj,jk) = vmask(ji,jj,jk) * ( trb(ji,jj+1,jk,jn) - trb(ji,jj,jk,jn) )
               END DO
            END DO
         END DO
         ! bottom values
         zt1(:,:,jpk) = 0.e0
         zt2(:,:,jpk) = 0.e0

         ! lateral boundary conditions on zt1, zt2
         CALL lbc_lnk( zt1, 'U', -1. )   
         CALL lbc_lnk( zt2, 'V', -1. ) 


         ! Slopes
         ! interior values
         DO jk = 1, jpkm1
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji-1,jj  ,jk) )   &
                     &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji-1,jj  ,jk) ) )
                  ztp2(ji,jj,jk) =                    ( zt2(ji,jj,jk) + zt2(ji  ,jj-1,jk) )   &
                     &           * ( 0.25 + SIGN( 0.25, zt2(ji,jj,jk) * zt2(ji  ,jj-1,jk) ) )
               END DO
            END DO
         END DO
         ! bottom values
         ztp1(:,:,jpk) = 0.e0
         ztp2(:,:,jpk) = 0.e0

         ! Slopes limitation
         DO jk = 1, jpkm1
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
                     &           * MIN(    ABS( ztp1(ji  ,jj,jk) ),   &
                     &                  2.*ABS( zt1 (ji-1,jj,jk) ),   &
                     &                  2.*ABS( zt1 (ji  ,jj,jk) ) )

                  ztp2(ji,jj,jk) = SIGN( 1., ztp2(ji,jj,jk) )   &
                     &           * MIN(    ABS( ztp2(ji,jj  ,jk) ),   &
                     &                  2.*ABS( zt2 (ji,jj-1,jk) ),   &
                     &                  2.*ABS( zt2 (ji,jj  ,jk) ) )

               END DO
            END DO
         END DO

         ! Advection terms
         ! interior values
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! volume fluxes
#if defined key_s_coord || defined key_partial_steps
                  zeu = e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
                  zev = e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
#else
                  zeu = e2u(ji,jj) * zun(ji,jj,jk)
                  zev = e1v(ji,jj) * zvn(ji,jj,jk)
#endif
                  ! MUSCL fluxes
                  z2dtt = rdttra(jk) * FLOAT(ndttrc)
                  z0u = SIGN( 0.5, zun(ji,jj,jk) )            
                  zalpha = 0.5 - z0u
                  zu  = z0u - 0.5 * zun(ji,jj,jk) * z2dtt / e1u(ji,jj)
                  zzt1 = trb(ji+1,jj,jk,jn) + zu*ztp1(ji+1,jj,jk)
                  zzt2 = trb(ji  ,jj,jk,jn) + zu*ztp1(ji  ,jj,jk)
                  zt1(ji,jj,jk) = zeu * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
                  z0v = SIGN( 0.5, zvn(ji,jj,jk) )            
                  zalpha = 0.5 - z0v
                  zv  = z0v - 0.5 * zvn(ji,jj,jk) * z2dtt / e2v(ji,jj)
                  zzt1 = trb(ji,jj+1,jk,jn) + zv*ztp2(ji,jj+1,jk)
                  zzt2 = trb(ji,jj  ,jk,jn) + zv*ztp2(ji,jj  ,jk)
                  zt2(ji,jj,jk) = zev * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
               END DO
            END DO
         END DO

         ! lateral boundary conditions on zt1, zt2 (changed sign)
         CALL lbc_lnk( zt1, 'U', -1. ) 
         CALL lbc_lnk( zt2, 'V', -1. ) 

         ! Compute and add the horizontal advective trend

         DO jk = 1, jpkm1
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
                  zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
#else
                  zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
#endif
                  ! horizontal advective trends
                  ztra = - zbtr * ( zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk  )   &
                     &            + zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk  ) )
                  ! add it to the general tracer trends
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
#if defined key_trc_diatrd
                  ! recompute the trends in i- and j-direction as Uh gradh(T)
#   if defined key_s_coord || defined key_partial_steps
                  zfui =  e2u(ji  ,jj) * fse3u(ji,  jj,jk) * un(ji,  jj,jk)   &
                     & -  e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * un(ji-1,jj,jk)
                  zfvj =  e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * vn(ji,jj  ,jk)   &
                     & -  e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * vn(ji,jj-1,jk)
#   else
                  zfui = e2u(ji  ,jj) * un(ji,  jj,jk)   &
                     & - e2u(ji-1,jj) * un(ji-1,jj,jk)
                  zfvj = e1v(ji,jj  ) * vn(ji,jj  ,jk)   &
                     & - e1v(ji,jj-1) * vn(ji,jj-1,jk)
#   endif
                  ztai =-zbtr * (  zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk) - trn(ji,jj,jk,jn) * zfui  )
                  ztaj =-zbtr * (  zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk) - trn(ji,jj,jk,jn) * zfvj  )
                  ! save i- and j- advective trends computed as Uh gradh(T)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = ztai
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = ztaj
#endif
               END DO
            END DO
         END DO
      ENDDO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('muscl - had')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

         ! II. Vertical advective fluxes
         ! -----------------------------

      DO jn = 1, jptra
         ! First guess of the slope
         ! interior values
         DO jk = 2, jpkm1
            zt1(:,:,jk) = tmask(:,:,jk) * ( trb(:,:,jk-1,jn) - trb(:,:,jk,jn) )
         END DO
         ! surface and bottom boundary conditions
         zt1 (:,:, 1 ) = 0.e0 
         zt1 (:,:,jpk) = 0.e0
         ! Slopes
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji,jj,jk+1) )   &
                     &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji,jj,jk+1) ) )
               END DO
            END DO
         END DO

         ! Slopes limitation
         ! interior values
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
                     &           * MIN(    ABS( ztp1(ji,jj,jk  ) ),   &
                     &                  2.*ABS( zt1 (ji,jj,jk+1) ),   &
                     &                  2.*ABS( zt1 (ji,jj,jk  ) ) )
               END DO
            END DO
         END DO
         ! surface values
         ztp1(:,:,1) = 0. 
         ! vertical advective flux
         ! interior values
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z2dtt = rdttra(jk) * FLOAT(ndttrc)
                  zew = zwn(ji,jj,jk+1)
                  z0w = SIGN( 0.5, zwn(ji,jj,jk+1) )
                  zalpha = 0.5 + z0w
                  zw  = z0w - 0.5 * zwn(ji,jj,jk+1)*z2dtt / fse3w(ji,jj,jk+1)
                  zzt1 = trb(ji,jj,jk+1,jn) + zw*ztp1(ji,jj,jk+1)
                  zzt2 = trb(ji,jj,jk  ,jn) + zw*ztp1(ji,jj,jk  )
                  zt1(ji,jj,jk+1) = zew * ( zalpha * zzt1 + (1.-zalpha)*zzt2 )
               END DO
            END DO
         END DO
         ! surface values
         IF( lk_dynspg_rl ) THEN        ! rigid lid : flux set to zero
            zt1(:,:, 1 ) = 0.e0
         ELSE                           ! free surface
            zt1(:,:, 1 ) = zwn(:,:,1) * trb(:,:,1,jn)
         ENDIF

         ! bottom values
         zt1(:,:,jpk) = 0.e0

         ! Compute & add the vertical advective trend

         DO jk = 1, jpkm1
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1. / fse3t(ji,jj,jk)
                  ! horizontal advective trends
                  ztra = - zbtr * ( zt1(ji,jj,jk) - zt1(ji,jj,jk+1) )
                  ! add it to the general tracer trends
                  tra(ji,jj,jk,jn) =  tra(ji,jj,jk,jn) + ztra
#if defined key_trc_diatrd
                  ! save the vertical advective trends computed as w gradz(T)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = ztra - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
#endif
               END DO
            END DO
         END DO

      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('muscl - zad')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

END SUBROUTINE trc_adv_muscl

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_adv_muscl( kt )  
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_adv_muscl: You should not have seen this print! error?', kt
   END SUBROUTINE trc_adv_muscl
#endif

   !!======================================================================
END MODULE trcadv_muscl
