MODULE trcadv_smolar
   !!==============================================================================
   !!                       ***  MODULE  trcadv_smolar  ***
   !! Ocean passive tracers:  horizontal & vertical advective trend
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_adv_smolar : update the passive tracer trend with the horizontal
   !!                  and vertical advection trends using a Smolarkiewicz 
   !!                  FCT scheme
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables
   USE lbclnk              ! ocean lateral boundary conditions (or mpp link)
   USE trcbbl              ! advective passive tracers in the BBL
   USE prtctl_trc      ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_adv_smolar    ! routine called by trcstp.F90

   !! * Module variable
   REAL(wp), DIMENSION(jpk) ::   &
      rdttrc                     ! vertical profile of tracer time-step
 
   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_smolar.F90,v 1.11 2006/04/10 15:38:54 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_adv_smolar( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_adv_smolar  ***
      !!
      !! ** Purpose :   Compute the now trend due to total advection of passi-
      !!      ve tracer using a Smolarkiewicz FCT (Flux Corrected Transport ) 
      !!      scheme and add it to the general tracer trend.
      !!
      !! ** Method : Computation of not exactly the advection but the
      !!             transport term, i.e.  div(u*tra). 
      !!             Computes the now horizontal and vertical advection with
      !!             the complete 3d method.
      !! 
      !!       note: - sc is an empirical factor to be used with care
      !!             - this advection scheme needs an euler-forward time scheme
      !!
      !! ** Action : - update tra with the now advective tracer trends
      !!             - save trends in trtrd ('key_trc_diatrd')
      !!
      !! References :                
      !!     Piotr K. Smolarkiewicz, 1983,
      !!       "A simple positive definit advection
      !!        scheme with small IMPLICIT diffusion"
      !!        Monthly Weather Review, pp 479-486
      !!
      !! History :
      !!         !  87-06 (pa-dl) Original
      !!         !  91-11 (G. Madec)
      !!         !  94-08 (A. Czaja)
      !!         !  95-09 (M. Levy) passive tracers
      !!         !  98-03 (M.A. Foujols) lateral boundary conditions
      !!         !  99-02 (M.A. Foujols) lbc in conjonction with ORCA
      !!         !  00-05 (MA Foujols) add lbc for tracer trends
      !!         !  00-10 (MA Foujols and E.Kestenare) INCLUDE instead of routine
      !!         !  01-05 (E.Kestenare) fix bug in trtrd indexes
      !!         !  02-05 (M-A Filiberti, and M.Levy) correction in trtrd computation
      !!   9.0   !  03-04  (C. Ethe)  F90: Free form and module
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
      INTEGER :: ji, jj, jk,jt, jn            ! dummy loop indices

      REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &      
         zti, ztj,              &
         zaa, zbb, zcc,         &
         zx , zy , zz ,         &
         zkx, zky, zkz,         &
         zbuf

#if defined key_trc_diatrd
      REAL(wp) :: zgm, zgz
#endif

      REAL(wp) :: zbtr, ztra
      REAL(wp) :: zfp_ui, zfp_vj, zfm_ui, zfm_vj, zfp_w, zfm_w
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------


      IF( kt == nittrc000  .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_adv_smolar : SMOLARKIEWICZ advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         rdttrc(:) = rdttra(:) * FLOAT(ndttrc)
      ENDIF


#if defined key_trcbbl_adv        
      ! Advective bottom boundary layer
      ! -------------------------------
      zun(:,:,:) = un (:,:,:) - u_trc_bbl(:,:,:)
      zvn(:,:,:) = vn (:,:,:) - v_trc_bbl(:,:,:)
      zwn(:,:,:) = wn (:,:,:) + w_trc_bbl( :,:,:)
#endif

      ! tracer loop parallelized (macrotasking)
      ! =======================================
      
      DO jn = 1, jptra
         
         ! 1. tracer flux in the 3 directions
         ! ----------------------------------
         
         ! 1.1 mass flux at u v and t-points and initialization

        DO jk = 1,jpk

           DO jj = 1,jpj
              DO ji = 1,jpi
                 zaa(ji,jj,jk) = e2u(ji,jj)*fse3u(ji,jj,jk) * zun(ji,jj,jk)
                 zbb(ji,jj,jk) = e1v(ji,jj)*fse3v(ji,jj,jk) * zvn(ji,jj,jk)
                 zcc(ji,jj,jk) = e1t(ji,jj)*e2t(ji,jj)      * zwn(ji,jj,jk)
                 zbuf(ji,jj,jk) = 0.
                 ztj(ji,jj,jk) = 0.
                 zx(ji,jj,jk) = 0.
                 zy(ji,jj,jk) = 0.
                 zz(ji,jj,jk) = 0.
                 zti(ji,jj,jk) = trn(ji,jj,jk,jn)
#if defined key_trc_diatrd
                 IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = 0.
                 IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = 0.
                 IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = 0.
#endif
              END DO
           END DO
           
           ! 1.2 calcul of intermediate field with an upstream advection scheme
           !     and mass fluxes calculated above
           
           ! calcul of tracer flux in the i and j direction
           
           DO jj=1,jpj
              zkx(  1,jj,jk)=0.
              zkx(jpi,jj,jk)=0.
           END DO
           
           DO ji=1,jpi
              zky(ji,  1,jk)=0.
              zky(ji,jpj,jk)=0.
           END DO
           
           DO jj = 2,jpjm1
              DO ji = 2,jpim1
                 zfp_ui = 0.5 * ( zaa(ji,jj,jk) + ABS( zaa(ji,jj,jk) ) )
                 zfp_vj = 0.5 * ( zbb(ji,jj,jk) + ABS( zbb(ji,jj,jk) ) )
                 zfm_ui = 0.5 * ( zaa(ji,jj,jk) - ABS( zaa(ji,jj,jk) ) )
                 zfm_vj = 0.5 * ( zbb(ji,jj,jk) - ABS( zbb(ji,jj,jk) ) )            
                 zkx(ji,jj,jk) = zfp_ui * zti(ji,jj,jk) + zfm_ui * zti(ji+1,jj  ,jk) 
                 zky(ji,jj,jk) = zfp_vj * zti(ji,jj,jk) + zfm_vj * zti(ji    ,jj+1,jk)             
              END DO
           END DO

        END DO

         ! II. Vertical advection
         ! ----------------------

         ! Surface value
         IF( lk_dynspg_rl ) THEN        ! rigid lid : flux set to zero
            zkz(:,:, 1 ) = 0.e0  
         ELSE                           ! free surface
            zkz(:,:, 1 ) = zwn(:,:,1) * trn(:,:,1,jn) * tmask(ji,jj,1)
         ENDIF

        DO jk = 2,jpk
          DO jj = 1,jpj
            DO ji = 1,jpi
               zfp_w = 0.5 * ( zcc(ji,jj,jk) + ABS( zcc(ji,jj,jk) ) )
               zfm_w = 0.5 * ( zcc(ji,jj,jk) - ABS( zcc(ji,jj,jk) ) )        
               zkz(ji,jj,jk) = zfp_w * zti(ji,jj,jk) + zfm_w * zti(ji,jj,jk-1)     
            END DO
          END DO
        END DO

! ... Lateral boundary conditions on zk[xy]
      CALL lbc_lnk( zkx, 'U', -1. )
      CALL lbc_lnk( zky, 'V', -1. )


! 2. calcul of after field using an upstream advection scheme
! -----------------------------------------------------------

        DO jk = 1,jpkm1
          DO jj = 2,jpjm1
            DO ji = 2,jpim1
              zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
              ztj(ji,jj,jk) = -zbtr*    &
     &            ( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk)  &
     &            + zky(ji,jj,jk) - zky(ji,jj - 1,jk)  &
     &            + zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
#if defined key_trc_diatrd
              IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) -  &
     &                       zbtr*( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk) )

              IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) -  &
     &            zbtr*( zky(ji,jj,jk) - zky(ji,jj - 1,jk) )

              IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3) -  &
     &            zbtr*( zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
#endif
            END DO
          END DO
        END DO

! 2.1 start of antidiffusive correction loop

        DO jt = 1,ncortrc

! 2.2 calcul of intermediary field zti

          DO jk = 1,jpkm1
            DO jj = 2,jpjm1
              DO ji = 2,jpim1
                zti(ji,jj,jk) = zti(ji,jj,jk)+rdttrc(jk)*ztj(ji,jj,jk)
                zbuf(ji,jj,jk) = zbuf(ji,jj,jk) + ztj(ji,jj,jk)
              END DO
            END DO
          END DO

! ... Lateral boundary conditions on zti
      CALL lbc_lnk( zti, 'T', 1. )


! 2.3 calcul of the antidiffusive flux

          DO jk = 1,jpkm1
            DO jj = 2,jpjm1
              DO ji = 2,jpim1
                zx(ji,jj,jk) = ( abs(zaa(ji,jj,jk)) - rdttrc(jk)       &
     &              *zaa(ji,jj,jk)**2/                          &
     &              (e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) ) )  &
     &              *(zti(ji + 1,jj,jk) - zti( ji ,jj,jk))      &
     &              /(zti( ji ,jj,jk) + zti(ji + 1,jj,jk) + rtrn)    &
     &              * rsc

                zy(ji,jj,jk) = ( abs(zbb(ji,jj,jk)) - rdttrc(jk)       &
     &              *zbb(ji,jj,jk)**2/                          &
     &              (e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) ) )  &
     &              *(zti(ji,jj + 1,jk) - zti(ji, jj ,jk))      &
     &              /(zti(ji, jj ,jk) + zti(ji,jj + 1,jk) + rtrn)    &
     &              * rsc
              END DO
            END DO
          END DO

          DO jk = 2,jpkm1
            DO jj = 2,jpjm1
              DO ji = 2,jpim1
                zz(ji,jj,jk) = ( abs(zcc(ji,jj,jk)) - rdttrc(jk)*zcc(ji,jj,jk)**2  &  
     &              /( e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk) ) ) &
     &               *( zti(ji,jj,jk) - zti(ji,jj,jk - 1) )/  &
     &                ( zti(ji,jj,jk) + zti(ji,jj,jk - 1) + rtrn )* rsc*( -1.)
              END DO
            END DO
          END DO

! 2.4 cross terms

          IF (crosster) THEN
              DO jk = 2,jpkm1
                DO jj = 2,jpjm1
                  DO ji = 2,jpim1
                    zx(ji,jj,jk) = zx(ji,jj,jk) &
     &                  - 0.5*rdttrc(jk)*rsc*zaa(ji,jj,jk)*0.25* &
     &                  (    (zbb(ji  ,jj - 1,jk  ) + zbb(ji + 1,jj - 1 &
     &                  ,jk  ) + zbb(ji + 1,jj  ,jk  ) + zbb(ji  ,jj  &
     &                  ,jk))* (zti(ji  ,jj + 1,jk  ) + zti(ji + 1,jj + &
     &                  1,jk  ) - zti(ji + 1,jj - 1,jk  ) - zti(ji  ,jj &
     &                  - 1,jk  ))/ (zti(ji  ,jj + 1,jk  ) + zti(ji + 1 &
     &                  ,jj + 1,jk  ) + zti(ji + 1,jj - 1,jk  ) + zti(ji &
     &                  ,jj - 1,jk  ) + rtrn) + (zcc(ji  ,jj  ,jk  ) + &
     &                  zcc(ji + 1,jj  ,jk  ) + zcc(ji  ,jj  ,jk + 1) + &
     &                  zcc(ji + 1,jj  ,jk + 1))* (zti(ji  ,jj  ,jk - 1) &
     &                  + zti(ji + 1,jj  ,jk - 1) - zti(ji  ,jj  ,jk + 1 &
     &                  )- zti(ji + 1,jj  ,jk + 1))/ (zti(ji  ,jj  ,jk - &
     &                  1) + zti(ji + 1,jj  ,jk - 1) + zti(ji  ,jj  ,jk &
     &                  +1) + zti(ji + 1,jj  ,jk + 1) + rtrn))/(e1u(ji &
     &                  ,jj)*e2u(ji,jj)*fse3u(ji,jj,jk))*vmask(ji  ,jj - &
     &                  1,jk  )*vmask(ji + 1,jj - 1,jk  )*vmask(ji + 1 &
     &                  ,jj,jk)*vmask(ji  ,jj  ,jk  )*tmask(ji  ,jj  ,jk &
     &                  )*tmask(ji + 1,jj  ,jk  )*tmask(ji  ,jj  ,jk + 1 &
     &                  )*tmask(ji + 1,jj  ,jk + 1)

                    zy(ji,jj,jk) = zy(ji,jj,jk)    &   
     &                  - 0.5*rdttrc(jk)*rsc*zbb(ji,jj,jk)*0.25*    &   
     &                  (    (zaa(ji - 1,jj  ,jk  ) + zaa(ji - 1,jj + 1    &   
     &                  ,jk  ) + zaa(ji  ,jj  ,jk  ) + zaa(ji  ,jj + 1    &   
     &                  ,jk))* (zti(ji + 1,jj + 1,jk  ) + zti(ji + 1,jj    &   
     &                  ,jk  ) - zti(ji - 1,jj + 1,jk  ) - zti(ji - 1,jj    &   
     &                  ,jk  ))/ (zti(ji + 1,jj + 1,jk  ) + zti(ji + 1    &   
     &                  ,jj  ,jk  ) + zti(ji - 1,jj + 1,jk  ) + zti(ji    &   
     &                  - 1,jj  ,jk  ) + rtrn) + (zcc(ji  ,jj  ,jk  )    &   
     &                  + zcc(ji  ,jj  ,jk + 1) + zcc(ji  ,jj + 1,jk  )    &   
     &                  + zcc(ji  ,jj + 1,jk + 1))* (zti(ji  ,jj  ,jk -    &   
     &                  1) + zti(ji  ,jj + 1,jk - 1) - zti(ji  ,jj  ,jk    &   
     &                  +1) - zti(ji  ,jj + 1,jk + 1))/ (zti(ji  ,jj    &   
     &                  ,jk- 1) + zti(ji  ,jj + 1,jk - 1) + zti(ji  ,jj    &   
     &                  ,jk+ 1) + zti(ji  ,jj + 1,jk + 1) + rtrn))    &   
     &                  /(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk))    &   
     &                  *umask(ji - 1,jj,jk  )*umask(ji - 1,jj + 1,jk  )    &   
     &                  *umask(ji  ,jj,jk  )*umask(ji  ,jj + 1,jk  )    &   
     &                  *tmask(ji  ,jj,jk)*tmask(ji  ,jj  ,jk + 1)    &   
     &                  *tmask(ji  ,jj + 1,jk)*tmask(ji  ,jj + 1,jk + 1)      

                    zz(ji,jj,jk) = zz(ji,jj,jk)    &   
     &                  - 0.5*rdttrc(jk)*rsc*zcc(ji,jj,jk)*0.25*    &   
     &                  (    (zaa(ji - 1,jj  ,jk  ) + zaa(ji  ,jj  ,jk    &   
     &                  ) + zaa(ji  ,jj  ,jk - 1) + zaa(ji - 1,jj  ,jk -    &   
     &                  1))*(zti(ji + 1,jj  ,jk - 1) + zti(ji + 1,jj    &   
     &                  ,jk  ) - zti(ji - 1,jj  ,jk  ) - zti(ji - 1,jj    &   
     &                  ,jk - 1))/(zti(ji + 1,jj  ,jk - 1) + zti(ji + 1    &   
     &                  ,jj,jk  ) + zti(ji - 1,jj  ,jk  ) + zti(ji - 1    &   
     &                  ,jj,jk - 1) + rtrn) + (zbb(ji  ,jj - 1,jk  )    &   
     &                  + zbb(ji  ,jj  ,jk  ) + zbb(ji  ,jj  ,jk - 1)    &   
     &                  + zbb(ji  ,jj - 1,jk - 1))*(zti(ji  ,jj + 1,jk -    &   
     &                  1) + zti(ji  ,jj + 1,jk  ) - zti(ji  ,jj - 1,jk    &   
     &                  ) - zti(ji  ,jj - 1,jk - 1))/(zti(ji  ,jj + 1,jk    &   
     &                  - 1) + zti(ji  ,jj + 1,jk  ) + zti(ji  ,jj - 1    &   
     &                  ,jk  ) + zti(ji  ,jj - 1,jk - 1) + rtrn))    &   
     &                  /(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk))    &   
     &                  *umask(ji - 1,jj,jk  )*umask(ji  ,jj  ,jk  )    &   
     &                  *umask(ji  ,jj,jk- 1)*umask(ji - 1,jj  ,jk - 1)    &   
     &                  *vmask(ji  ,jj- 1,jk)*vmask(ji  ,jj  ,jk  )    &   
     &                  *vmask(ji  ,jj  ,jk-1)*vmask(ji  ,jj - 1,jk - 1)       
                  END DO
                END DO
              END DO

              DO jj = 2,jpjm1
                DO ji = 2,jpim1
                  zx(ji,jj,1) = zx(ji,jj,1)    &   
     &                - 0.5*rdttrc(jk)*rsc*zaa(ji,jj,1)*0.25*    &   
     &                ( (zbb(ji  ,jj - 1,1  ) + zbb(ji + 1,jj - 1,1  )    &   
     &                + zbb(ji + 1,jj  ,1  ) + zbb(ji  ,jj  ,1  ))    &   
     &                *(zti(ji  ,jj + 1,1  ) + zti(ji + 1,jj + 1,1  )    &   
     &                - zti(ji + 1,jj - 1,1  ) - zti(ji  ,jj - 1,1  ))    &   
     &                /(zti(ji  ,jj + 1,1  ) + zti(ji + 1,jj + 1,1  )    &   
     &                + zti(ji + 1,jj - 1,1  ) + zti(ji  ,jj - 1,1  ) +    &   
     &                rtrn))/(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,1))    &   
     &                *vmask(ji  ,jj - 1,1  )*vmask(ji + 1,jj - 1,1  )    &   
     &                *vmask(ji + 1,jj  ,1  )*vmask(ji  ,jj  ,1  )    

                 zy(ji,jj,1) = zy(ji,jj,1)    &   
     &                - 0.5*rdttrc(jk)*rsc*zbb(ji,jj,1)*0.25*    &   
     &                ( (zaa(ji-1  ,jj ,1  ) + zaa(ji - 1,jj + 1,1  )    &   
     &                + zaa(ji ,jj  ,1  ) + zaa(ji  ,jj + 1  ,1  ))    &   
     &                *(zti(ji + 1,jj + 1,1  ) + zti(ji + 1,jj ,1  )    &   
     &                - zti(ji - 1,jj + 1,1  ) - zti(ji - 1,jj ,1  ))    &   
     &                /(zti(ji + 1,jj + 1,1  ) + zti(ji + 1,jj ,1  )    &   
     &                + zti(ji - 1,jj + 1,1  ) + zti(ji - 1,jj ,1  ) +    &   
     &                rtrn))/(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,1))    &   
     &                *umask(ji - 1,jj,1  )*umask(ji - 1,jj + 1,1  )    &   
     &                *umask(ji    ,jj,1  )*umask(ji  ,jj + 1 ,1  )    

                END DO
              END DO
          ENDIF

          ! ... Lateral boundary conditions on z[xyz]
          CALL lbc_lnk( zx, 'U', -1. )
          CALL lbc_lnk( zy, 'V', -1. )
          CALL lbc_lnk( zz, 'W',  1. )

! 2.4 reinitialization

          DO jk = 1,jpk
            DO jj = 1,jpj
              DO ji = 1,jpi
                zaa(ji,jj,jk) = zx(ji,jj,jk)
                zbb(ji,jj,jk) = zy(ji,jj,jk)
                zcc(ji,jj,jk) = zz(ji,jj,jk)
              END DO
            END DO
          END DO

! 2.5 calcul of the final field:
!    advection by antidiffusive mass fluxes and an upstream scheme

          DO jk = 1,jpk
             DO jj = 2,jpjm1
                DO ji = 2,jpim1
                   zfp_ui = 0.5 * ( zaa(ji,jj,jk) + ABS( zaa(ji,jj,jk) ) )
                   zfp_vj = 0.5 * ( zbb(ji,jj,jk) + ABS( zbb(ji,jj,jk) ) )
                   zfm_ui = 0.5 * ( zaa(ji,jj,jk) - ABS( zaa(ji,jj,jk) ) )
                   zfm_vj = 0.5 * ( zbb(ji,jj,jk) - ABS( zbb(ji,jj,jk) ) )            
                   zkx(ji,jj,jk) = zfp_ui * zti(ji,jj,jk) + zfm_ui * zti(ji+1,jj  ,jk) 
                   zky(ji,jj,jk) = zfp_vj * zti(ji,jj,jk) + zfm_vj * zti(ji    ,jj+1,jk)             
                END DO
             END DO
          END DO

          DO jk = 2,jpk
             DO jj = 1,jpj
                DO ji = 1,jpi
                   zfp_w = 0.5 * ( zcc(ji,jj,jk) + ABS( zcc(ji,jj,jk) ) )
                   zfm_w = 0.5 * ( zcc(ji,jj,jk) - ABS( zcc(ji,jj,jk) ) )       
                   zkz(ji,jj,jk) = zfp_w * zti(ji,jj,jk) + zfm_w * zti(ji,jj,jk-1)     
                END DO
             END DO
          END DO


! ... Lateral boundary conditions on zk[xy]
      CALL lbc_lnk( zkx, 'U', -1. )
      CALL lbc_lnk( zky, 'V', -1. )


! 2.6. calcul of after field using an upstream advection scheme

          DO jk = 1,jpkm1
            DO jj = 2,jpjm1
              DO ji = 2,jpim1
                zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
                ztj(ji,jj,jk) = -zbtr*     &  
     &              ( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk)    &  
     &              + zky(ji,jj,jk) - zky(ji,jj - 1,jk)    &  
     &              + zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
#if defined key_trc_diatrd
                IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) -    &  
     &              zbtr*( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk) )   

                IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) -    &  
     &              zbtr*( zky(ji,jj,jk) - zky(ji,jj - 1,jk) )

                IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3) -    &  
     &              zbtr*( zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
#endif
              END DO
            END DO
          END DO

! 2.6 END of antidiffusive correction loop

        END DO

! 3. trend due to horizontal and vertical advection of tracer jn
! --------------------------------------------------------------

        DO jk = 1,jpk
          DO jj = 2,jpjm1
            DO ji = 2,jpim1
              ztra = ( zbuf(ji,jj,jk) + ztj(ji,jj,jk) ) * tmask(ji,jj,jk)
              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
            END DO
          END DO
        END DO

! 4.0 convert the transport trend into advection trend
! ----------------------------------------------------

#if defined key_trc_diatrd
        DO jk = 1,jpk
          DO jj = 2,jpjm1
            DO  ji = 2,jpim1
              zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
              zgm = zbtr * trn(ji,jj,jk,jn) *     &  
     &            ( zun(ji  ,jj,jk) * e2u(ji  ,jj) * fse3u(ji  ,jj,jk)    &  
     &            -zun(ji-1,jj,jk) * e2u(ji-1,jj) * fse3u(ji-1,jj,jk))

              zgz = zbtr * trn(ji,jj,jk,jn) *     &  
     &            ( zvn(ji,jj  ,jk) * e1v(ji,jj  ) * fse3v(ji,jj  ,jk)    &  
     &            -zvn(ji,jj-1,jk) * e1v(ji,jj-1) * fse3v(ji,jj-1,jk))

              IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) + zgm
              IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) + zgz
              IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3)    &  
     &            - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
            END DO
          END DO
        END DO

        ! Lateral boundary conditions on trtrd:

        IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),1), 'T', 1. )
        IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),2), 'T', 1. )
        IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),3), 'T', 1. )
#endif

 
        ! END of tracer loop
        ! ==================
     ENDDO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('smolar - adv')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF
     
  END SUBROUTINE trc_adv_smolar

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_adv_smolar( kt ) 
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_adv_smolar: You should not have seen this print! error?', kt
   END SUBROUTINE trc_adv_smolar
#endif

   !!======================================================================
END MODULE trcadv_smolar
