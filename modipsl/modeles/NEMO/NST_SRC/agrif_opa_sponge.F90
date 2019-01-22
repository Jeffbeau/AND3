#define SPONGE

      Module agrif_opa_sponge
#if defined key_agrif
      USE par_oce
      USE oce
      USE dom_oce
      
      Contains


      Subroutine Agrif_Sponge_Tra( kt )

      implicit none

      INTEGER :: kt
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: tabtemp, tbdiff, sbdiff
      INTEGER :: ji,jj,jk
      REAL(wp) :: viscsponge
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: umasktemp,vmasktemp
      INTEGER :: spongearea
      integer ipt,jpt
      real,dimension(:,:),pointer :: e1tparent,e2tparent
      REAL(wp), DIMENSION(jpi,jpj) :: localviscsponge
      real(wp) :: timecoeff
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztu ,ztv, zsu ,zsv
      REAL(wp) :: zta, zsa, zabe1, zabe2, zbtr

#include "domzgr_substitute.h90"

        
#if defined SPONGE

      timecoeff = real(Agrif_NbStepint())/Agrif_rhot()

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .TRUE.
      tabtemp = 0.
      Call Agrif_Bc_Variable(tabtemp, ta,calledweight=timecoeff,procname=interptn)
      Agrif_UseSpecialValue = .FALSE.

      tbdiff(:,:,:) = tb(:,:,:) - tabtemp(:,:,:)

      tabtemp = 0.
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .TRUE.
      Call Agrif_Bc_Variable(tabtemp, sa,calledweight=timecoeff,procname=interpsn)
      Agrif_UseSpecialValue = .FALSE.

      sbdiff(:,:,:) = sb(:,:,:) - tabtemp(:,:,:)
       
      viscsponge = rdt

      spongearea = 2 + 2 * Agrif_irhox()

      localviscsponge = 0.
      umasktemp = 0.
      vmasktemp = 0.

      IF ((nbondi == -1).OR.(nbondi == 2)) THEN

        DO ji = 2, spongearea
          localviscsponge(ji,:) = viscsponge * (spongearea-ji)/real(spongearea-2)
        ENDDO

        DO jk = 1, jpkm1
          umasktemp(2:spongearea-1,:,jk) = umask(2:spongearea-1,:,jk) &
           * 0.5 * (localviscsponge(2:spongearea-1,:) + localviscsponge(3:spongearea,:))
        ENDDO

       DO jk = 1, jpkm1
         vmasktemp(2:spongearea,1:jpjm1,jk) = vmask(2:spongearea,1:jpjm1,jk) &
          * 0.5 * (localviscsponge(2:spongearea,1:jpjm1) + localviscsponge(2:spongearea,2:jpj))
       ENDDO

      ENDIF

      IF ((nbondi == 1).OR.(nbondi == 2)) THEN

        DO ji = nlci-spongearea + 1,nlci-1
          localviscsponge(ji,:) = viscsponge * (ji - (nlci-spongearea+1))/real(spongearea-2)
        ENDDO

       DO jk = 1, jpkm1
        umasktemp(nlci-spongearea + 1:nlci-2,:,jk) = umask(nlci-spongearea + 1:nlci-2,:,jk) &
         * 0.5 * (localviscsponge(nlci-spongearea + 1:nlci-2,:) + localviscsponge(nlci-spongearea + 2:nlci-1,:))
       ENDDO

       DO jk = 1, jpkm1
        vmasktemp(nlci-spongearea + 1:nlci-1,1:jpjm1,jk) = vmask(nlci-spongearea + 1:nlci-1,1:jpjm1,jk) &
          * 0.5 * (localviscsponge(nlci-spongearea + 1:nlci-1,1:jpjm1) + localviscsponge(nlci-spongearea + 1:nlci-1,2:jpj))
       ENDDO

      ENDIF



      IF ((nbondj == -1).OR.(nbondj == 2)) THEN

        DO jj = 2, spongearea
        localviscsponge(:,jj) = viscsponge * (spongearea-jj)/real(spongearea-2)
        ENDDO

      DO jk = 1, jpkm1
      vmasktemp(:,2:spongearea-1,jk) = vmask(:,2:spongearea-1,jk) &
         * 0.5 * (localviscsponge(:,2:spongearea-1) + localviscsponge(:,3:spongearea))
        ENDDO

      DO jk = 1, jpkm1
       umasktemp(1:jpim1,2:spongearea,jk) = umask(1:jpim1,2:spongearea,jk) &
         * 0.5 * (localviscsponge(1:jpim1,2:spongearea) + localviscsponge(2:jpi,2:spongearea))
      ENDDO

      ENDIF

      IF ((nbondj == 1).OR.(nbondj == 2)) THEN

        DO jj = nlcj-spongearea + 1,nlcj-1
       localviscsponge(:,jj) = viscsponge * (jj - (nlcj-spongearea+1))/real(spongearea-2)
       ENDDO

      DO jk = 1, jpkm1
       vmasktemp(:,nlcj-spongearea + 1:nlcj-2,jk) = vmask(:,nlcj-spongearea + 1:nlcj-2,jk) &
          * 0.5 * (localviscsponge(:,nlcj-spongearea + 1:nlcj-2) + localviscsponge(:,nlcj-spongearea + 2:nlcj-1))
       ENDDO

      DO jk = 1, jpkm1
        umasktemp(1:jpim1,nlcj-spongearea + 1:nlcj-1,jk) = umask(1:jpim1,nlcj-spongearea + 1:nlcj-1,jk) &
         * 0.5 * (localviscsponge(1:jpim1,nlcj-spongearea + 1:nlcj-1) + localviscsponge(2:jpi,nlcj-spongearea + 1:nlcj-1))
      ENDDO

      ENDIF

      IF (.Not. spongedoneT) THEN
         zspe1ur(:,:) = e2u(:,:) / e1u(:,:)
         zspe2vr(:,:) = e1v(:,:) / e2v(:,:)
         zspbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:))
         
         spongedoneT = .TRUE.
      ENDIF

        DO jk = 1, jpkm1
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
#if defined key_s_coord || defined key_partial_steps
               zabe1 = umasktemp(ji,jj,jk) * zspe1ur(ji,jj) * fse3u(ji,jj,jk)
               zabe2 = vmasktemp(ji,jj,jk) * zspe2vr(ji,jj) * fse3v(ji,jj,jk)
#else
               zabe1 = umasktemp(ji,jj,jk) * zspe1ur(ji,jj)
               zabe2 = vmasktemp(ji,jj,jk) * zspe2vr(ji,jj)
#endif
               ztu(ji,jj,jk) = zabe1 * ( tbdiff(ji+1,jj  ,jk) - tbdiff(ji,jj,jk) )
               zsu(ji,jj,jk) = zabe1 * ( sbdiff(ji+1,jj  ,jk) - sbdiff(ji,jj,jk) )
               ztv(ji,jj,jk) = zabe2 * ( tbdiff(ji  ,jj+1,jk) - tbdiff(ji,jj,jk) )
               zsv(ji,jj,jk) = zabe2 * ( sbdiff(ji  ,jj+1,jk) - sbdiff(ji,jj,jk) )
            ENDDO
          ENDDO

         DO jj = 2,jpjm1
            DO ji = 2,jpim1
#if defined key_s_coord || defined key_partial_steps
               zbtr = zspbtr2(ji,jj) / fse3t(ji,jj,jk)
#else
               zbtr = zspbtr2(ji,jj)
#endif
               ! horizontal diffusive trends
               zta = zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
                  &          + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
               zsa = zbtr * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk)   &
                  &          + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
               ! add it to the general tracer trends
               ta(ji,jj,jk) = (ta(ji,jj,jk) + zta)
               sa(ji,jj,jk) = (sa(ji,jj,jk) + zsa)
            END DO
         END DO

        ENDDO

#endif

      Return
      End Subroutine Agrif_Sponge_Tra
      
      Subroutine Agrif_Sponge_dyn( kt )

      implicit none

      INTEGER :: kt
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: tabtemp, ubdiff, vbdiff,rotdiff,hdivdiff
      INTEGER :: ji,jj,jk
      REAL(wp) :: viscsponge
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: umasktemp,vmasktemp
      INTEGER :: spongearea
      integer ipt,jpt
      real,dimension(:,:),pointer :: e1tparent,e2tparent
      REAL(wp), DIMENSION(jpi,jpj) :: localviscsponge
      real(wp) :: timecoeff
      REAL(wp):: ze2u, ze1v, zua, zva

#include "domzgr_substitute.h90"

#if defined SPONGE

      timecoeff = real(Agrif_NbStepint())/Agrif_rhot()

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .TRUE.
      tabtemp = 0.
      Call Agrif_Bc_Variable(tabtemp, ua,calledweight=timecoeff,procname=interpun)
      Agrif_UseSpecialValue = .FALSE.

      ubdiff(:,:,:) = ub(:,:,:) - tabtemp(:,:,:)

      tabtemp = 0.
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .TRUE.
      Call Agrif_Bc_Variable(tabtemp, va,calledweight=timecoeff,procname=interpvn)
      Agrif_UseSpecialValue = .FALSE.

      vbdiff(:,:,:) = vb(:,:,:) - tabtemp(:,:,:)
       
      viscsponge = rdt

      spongearea = 2 + 2 * Agrif_irhox()

      localviscsponge = 0.
      umasktemp = 0.
      vmasktemp = 0.

      IF ((nbondi == -1).OR.(nbondi == 2)) THEN

        DO ji = 2, spongearea
          localviscsponge(ji,:) = viscsponge * (spongearea-ji)/real(spongearea-2)
        ENDDO

        DO jk = 1, jpkm1
          umasktemp(2:spongearea-1,:,jk) = umask(2:spongearea-1,:,jk) &
           * 0.5 * (localviscsponge(2:spongearea-1,:) + localviscsponge(3:spongearea,:))
        ENDDO

       DO jk = 1, jpkm1
         vmasktemp(2:spongearea,1:jpjm1,jk) = vmask(2:spongearea,1:jpjm1,jk) &
          * 0.5 * (localviscsponge(2:spongearea,1:jpjm1) + localviscsponge(2:spongearea,2:jpj))
       ENDDO

      ENDIF

      IF ((nbondi == 1).OR.(nbondi == 2)) THEN

        DO ji = nlci-spongearea + 1,nlci-1
          localviscsponge(ji,:) = viscsponge * (ji - (nlci-spongearea+1))/real(spongearea-2)
        ENDDO

       DO jk = 1, jpkm1
        umasktemp(nlci-spongearea + 1:nlci-2,:,jk) = umask(nlci-spongearea + 1:nlci-2,:,jk) &
         * 0.5 * (localviscsponge(nlci-spongearea + 1:nlci-2,:) + localviscsponge(nlci-spongearea + 2:nlci-1,:))
       ENDDO

       DO jk = 1, jpkm1
        vmasktemp(nlci-spongearea + 1:nlci-1,1:jpjm1,jk) = vmask(nlci-spongearea + 1:nlci-1,1:jpjm1,jk) &
          * 0.5 * (localviscsponge(nlci-spongearea + 1:nlci-1,1:jpjm1) + localviscsponge(nlci-spongearea + 1:nlci-1,2:jpj))
       ENDDO

      ENDIF



      IF ((nbondj == -1).OR.(nbondj == 2)) THEN

        DO jj = 2, spongearea
        localviscsponge(:,jj) = viscsponge * (spongearea-jj)/real(spongearea-2)
        ENDDO

      DO jk = 1, jpkm1
      vmasktemp(:,2:spongearea-1,jk) = vmask(:,2:spongearea-1,jk) &
         * 0.5 * (localviscsponge(:,2:spongearea-1) + localviscsponge(:,3:spongearea))
        ENDDO

      DO jk = 1, jpkm1
       umasktemp(1:jpim1,2:spongearea,jk) = umask(1:jpim1,2:spongearea,jk) &
         * 0.5 * (localviscsponge(1:jpim1,2:spongearea) + localviscsponge(2:jpi,2:spongearea))
      ENDDO

      ENDIF

      IF ((nbondj == 1).OR.(nbondj == 2)) THEN

        DO jj = nlcj-spongearea + 1,nlcj-1
       localviscsponge(:,jj) = viscsponge * (jj - (nlcj-spongearea+1))/real(spongearea-2)
       ENDDO

      DO jk = 1, jpkm1
       vmasktemp(:,nlcj-spongearea + 1:nlcj-2,jk) = vmask(:,nlcj-spongearea + 1:nlcj-2,jk) &
          * 0.5 * (localviscsponge(:,nlcj-spongearea + 1:nlcj-2) + localviscsponge(:,nlcj-spongearea + 2:nlcj-1))
       ENDDO

      DO jk = 1, jpkm1
        umasktemp(1:jpim1,nlcj-spongearea + 1:nlcj-1,jk) = umask(1:jpim1,nlcj-spongearea + 1:nlcj-1,jk) &
         * 0.5 * (localviscsponge(1:jpim1,nlcj-spongearea + 1:nlcj-1) + localviscsponge(2:jpi,nlcj-spongearea + 1:nlcj-1))
      ENDDO

      ENDIF
      
      ubdiff = ubdiff * umasktemp
      vbdiff = vbdiff * vmasktemp
      
      hdivdiff = 0.
      rotdiff = 0.

      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         !                                             ! --------
         ! Horizontal divergence                       !   div
         !                                             ! --------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               hdivdiff(ji,jj,jk) =   &
                  (  e2u(ji,jj)*fse3u(ji,jj,jk) * & 
                  ubdiff(ji,jj,jk) - e2u(ji-1,jj  )* &
                  fse3u(ji-1,jj  ,jk)  * ubdiff(ji-1,jj  ,jk)       &
                   + e1v(ji,jj)*fse3v(ji,jj,jk) * &
                  vbdiff(ji,jj,jk) - e1v(ji  ,jj-1)* &
                  fse3v(ji  ,jj-1,jk)  * vbdiff(ji  ,jj-1,jk)  )    &
                  / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
#else
               hdivdiff(ji,jj,jk) = (  e2u(ji,jj) * ubdiff(ji,jj,jk) &
               - e2u(ji-1,jj  ) * ubdiff(ji-1,jj  ,jk)      &
     &               + e1v(ji,jj) * vbdiff(ji,jj,jk) - & 
     &              e1v(ji  ,jj-1) * vbdiff(ji  ,jj-1,jk)  )   &
     &            / ( e1t(ji,jj) * e2t(ji,jj) )
#endif
            END DO
         END DO
         
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               rotdiff(ji,jj,jk) = (  e2v(ji+1,jj  ) * vbdiff(ji+1,jj  ,jk) - e2v(ji,jj) * vbdiff(ji,jj,jk)    &
                  &              - e1u(ji  ,jj+1) * ubdiff(ji  ,jj+1,jk) + e1u(ji,jj) * ubdiff(ji,jj,jk)  ) &
                  &           * fmask(ji,jj,jk) / ( e1f(ji,jj) * e2f(ji,jj) )
            END DO
         END DO
                  
         ENDDO
         
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               ze2u = rotdiff (ji,jj,jk)*fse3f(ji,jj,jk)
               ze1v = hdivdiff(ji,jj,jk)
               ! horizontal diffusive trends
               zua = - ( ze2u - rotdiff (ji,jj-1,jk)* &
               fse3f(ji,jj-1,jk) ) / ( e2u(ji,jj) * fse3u(ji,jj,jk) )   &
                     + ( hdivdiff(ji+1,jj,jk) - ze1v      &
               ) / e1u(ji,jj)

               zva = + ( ze2u - rotdiff (ji-1,jj,jk)* &
               fse3f(ji-1,jj,jk) ) / ( e1v(ji,jj) * fse3v(ji,jj,jk) )   &
                     + ( hdivdiff(ji,jj+1,jk) - ze1v    &
                    ) / e2v(ji,jj)
#else
               ! horizontal diffusive trends
               ze2u = rotdiff (ji,jj,jk)
               ze1v = hdivdiff(ji,jj,jk)
               zua = - (                ze2u                  - &
               rotdiff (ji,jj-1,jk) ) / e2u(ji,jj)   &
                     + ( hdivdiff(ji+1,jj,jk) -     &
                ze1v                  ) / e1u(ji,jj)

               zva = + (                ze2u                  - &
               rotdiff (ji-1,jj,jk) ) / e1v(ji,jj)   &
                     + ( hdivdiff(ji,jj+1,jk) -       &
                ze1v                  ) / e2v(ji,jj)
#endif

               ! add it to the general momentum trends

               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

#endif

      Return
      End Subroutine Agrif_Sponge_dyn

       subroutine interptn(tabres,i1,i2,j1,j2,k1,k2)
       Implicit none
#  include "domzgr_substitute.h90"       
       integer i1,i2,j1,j2,k1,k2
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres

       tabres(i1:i2,j1:j2,k1:k2) = tn(i1:i2,j1:j2,k1:k2)

       end subroutine interptn    
       
       subroutine interpsn(tabres,i1,i2,j1,j2,k1,k2)
       Implicit none
#  include "domzgr_substitute.h90"       
       integer i1,i2,j1,j2,k1,k2
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres

       tabres(i1:i2,j1:j2,k1:k2) = sn(i1:i2,j1:j2,k1:k2)

       end subroutine interpsn  
                 
 
       subroutine interpun(tabres,i1,i2,j1,j2,k1,k2)
       Implicit none
#  include "domzgr_substitute.h90"       
       integer i1,i2,j1,j2,k1,k2
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres

       tabres(i1:i2,j1:j2,k1:k2) = un(i1:i2,j1:j2,k1:k2)

       end subroutine interpun 
       
       subroutine interpvn(tabres,i1,i2,j1,j2,k1,k2)
       Implicit none
#  include "domzgr_substitute.h90"       
       integer i1,i2,j1,j2,k1,k2
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres

       tabres(i1:i2,j1:j2,k1:k2) = vn(i1:i2,j1:j2,k1:k2)

       end subroutine interpvn 

#else
       CONTAINS
       subroutine agrif_opa_sponge_empty
       end subroutine agrif_opa_sponge_empty
#endif
                    
       End Module agrif_opa_sponge
