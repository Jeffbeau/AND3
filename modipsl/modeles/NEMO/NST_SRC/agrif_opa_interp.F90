!
      Module agrif_opa_interp
#if defined key_agrif
      USE par_oce
      USE oce
      USE dom_oce      
      USE sol_oce

      CONTAINS
      SUBROUTINE Agrif_tra( kt )

      Implicit none
      
   !! * Substitutions
#  include "domzgr_substitute.h90"  
#  include "vectopt_loop_substitute.h90"
!
      INTEGER :: kt
      REAL(wp) tatemp(jpi,jpj,jpk) , satemp(jpi,jpj,jpk)
      INTEGER :: ji,jj,jk
      REAL(wp) :: rhox
      REAL(wp) :: alpha1, alpha2, alpha3, alpha4
      REAL(wp) :: alpha5, alpha6, alpha7
!
        IF (Agrif_Root()) RETURN

           Agrif_SpecialValue=0.
           Agrif_UseSpecialValue = .TRUE.
           tatemp = 0.
           satemp = 0.

           Call Agrif_Bc_variable(tatemp,tn)
           Call Agrif_Bc_variable(satemp,sn)
           Agrif_UseSpecialValue = .FALSE.
       
           rhox = Agrif_Rhox()
   
           alpha1 = (rhox-1.)/2.
           alpha2 = 1.-alpha1
   
           alpha3 = (rhox-1)/(rhox+1)
           alpha4 = 1.-alpha3
   
           alpha6 = 2.*(rhox-1.)/(rhox+1.)
           alpha7 = -(rhox-1)/(rhox+3)
           alpha5 = 1. - alpha6 - alpha7
   
!
      If ((nbondi == 1).OR.(nbondi == 2)) THEN
      
      ta(nlci,:,:) = alpha1 * tatemp(nlci,:,:) + alpha2 * tatemp(nlci-1,:,:)
      sa(nlci,:,:) = alpha1 * satemp(nlci,:,:) + alpha2 * satemp(nlci-1,:,:)
      
      Do jk=1,jpk      
      Do jj=1,jpj
        IF (umask(nlci-2,jj,jk).EQ.0.) THEN
        ta(nlci-1,jj,jk) = ta(nlci,jj,jk) * tmask(nlci-1,jj,jk)
        sa(nlci-1,jj,jk) = sa(nlci,jj,jk) * tmask(nlci-1,jj,jk)
        ELSE
        ta(nlci-1,jj,jk)=(alpha4*ta(nlci,jj,jk)+alpha3*ta(nlci-2,jj,jk))*tmask(nlci-1,jj,jk)
        sa(nlci-1,jj,jk)=(alpha4*sa(nlci,jj,jk)+alpha3*sa(nlci-2,jj,jk))*tmask(nlci-1,jj,jk)
         IF (un(nlci-2,jj,jk).GT.0.) THEN
          ta(nlci-1,jj,jk)=(alpha6*ta(nlci-2,jj,jk)+alpha5*ta(nlci,jj,jk)+alpha7*ta(nlci-3,jj,jk))*tmask(nlci-1,jj,jk)
          sa(nlci-1,jj,jk)=(alpha6*sa(nlci-2,jj,jk)+alpha5*sa(nlci,jj,jk)+alpha7*sa(nlci-3,jj,jk))*tmask(nlci-1,jj,jk)
         ENDIF
        ENDIF
      End Do
      enddo 
      ENDIF        

      If ((nbondj == 1).OR.(nbondj == 2)) THEN
      
      ta(:,nlcj,:) = alpha1 * tatemp(:,nlcj,:) + alpha2 * tatemp(:,nlcj-1,:)
      sa(:,nlcj,:) = alpha1 * satemp(:,nlcj,:) + alpha2 * satemp(:,nlcj-1,:)
              
      Do jk=1,jpk      
      Do ji=1,jpi
        IF (vmask(ji,nlcj-2,jk).EQ.0.) THEN
        ta(ji,nlcj-1,jk) = ta(ji,nlcj,jk) * tmask(ji,nlcj-1,jk)
        sa(ji,nlcj-1,jk) = sa(ji,nlcj,jk) * tmask(ji,nlcj-1,jk)
        ELSE
        ta(ji,nlcj-1,jk)=(alpha4*ta(ji,nlcj,jk)+alpha3*ta(ji,nlcj-2,jk))*tmask(ji,nlcj-1,jk)        
        sa(ji,nlcj-1,jk)=(alpha4*sa(ji,nlcj,jk)+alpha3*sa(ji,nlcj-2,jk))*tmask(ji,nlcj-1,jk)
          IF (vn(ji,nlcj-2,jk) .GT. 0.) THEN
           ta(ji,nlcj-1,jk)=(alpha6*ta(ji,nlcj-2,jk)+alpha5*ta(ji,nlcj,jk)+alpha7*ta(ji,nlcj-3,jk))*tmask(ji,nlcj-1,jk)
           sa(ji,nlcj-1,jk)=(alpha6*sa(ji,nlcj-2,jk)+alpha5*sa(ji,nlcj,jk)+alpha7*sa(ji,nlcj-3,jk))*tmask(ji,nlcj-1,jk)
          ENDIF
        ENDIF
      End Do
      enddo
      ENDIF

      IF ((nbondi == -1).OR.(nbondi == 2)) THEN
      
      ta(1,:,:) = alpha1 * tatemp(1,:,:) + alpha2 * tatemp(2,:,:)
      sa(1,:,:) = alpha1 * satemp(1,:,:) + alpha2 * satemp(2,:,:)      
      
      Do jk=1,jpk      
      Do jj=1,jpj
        IF (umask(2,jj,jk).EQ.0.) THEN
        ta(2,jj,jk) = ta(1,jj,jk) * tmask(2,jj,jk)
        sa(2,jj,jk) = sa(1,jj,jk) * tmask(2,jj,jk)
        ELSE
        ta(2,jj,jk)=(alpha4*ta(1,jj,jk)+alpha3*ta(3,jj,jk))*tmask(2,jj,jk)        
        sa(2,jj,jk)=(alpha4*sa(1,jj,jk)+alpha3*sa(3,jj,jk))*tmask(2,jj,jk)
         IF (un(2,jj,jk).LT.0.) THEN
           ta(2,jj,jk)=(alpha6*ta(3,jj,jk)+alpha5*ta(1,jj,jk)+alpha7*ta(4,jj,jk))*tmask(2,jj,jk)
           sa(2,jj,jk)=(alpha6*sa(3,jj,jk)+alpha5*sa(1,jj,jk)+alpha7*sa(4,jj,jk))*tmask(2,jj,jk)
         ENDIF
        ENDIF
      End Do
      enddo
      ENDIF

      IF ((nbondj == -1).OR.(nbondj == 2)) THEN
      
      ta(:,1,:) = alpha1 * tatemp(:,1,:) + alpha2 * tatemp(:,2,:)
      sa(:,1,:) = alpha1 * satemp(:,1,:) + alpha2 * satemp(:,2,:)
            
      Do jk=1,jpk      
      Do ji=1,jpi
        IF (vmask(ji,2,jk).EQ.0.) THEN
        ta(ji,2,jk)=ta(ji,1,jk) * tmask(ji,2,jk)
        sa(ji,2,jk)=sa(ji,1,jk) * tmask(ji,2,jk)
        ELSE
        ta(ji,2,jk)=(alpha4*ta(ji,1,jk)+alpha3*ta(ji,3,jk))*tmask(ji,2,jk)
        sa(ji,2,jk)=(alpha4*sa(ji,1,jk)+alpha3*sa(ji,3,jk))*tmask(ji,2,jk) 
          IF (vn(ji,2,jk) .LT. 0.) THEN
            ta(ji,2,jk)=(alpha6*ta(ji,3,jk)+alpha5*ta(ji,1,jk)+alpha7*ta(ji,4,jk))*tmask(ji,2,jk)
            sa(ji,2,jk)=(alpha6*sa(ji,3,jk)+alpha5*sa(ji,1,jk)+alpha7*sa(ji,4,jk))*tmask(ji,2,jk)
          ENDIF
        ENDIF
      End Do
      enddo 
      ENDIF

      Return
      End Subroutine Agrif_tra
!
!
       SUBROUTINE Agrif_dyn(kt)
!
      USE phycst
      USE sol_oce
      USE in_out_manager

      implicit none
#  include "domzgr_substitute.h90"
!
      REAL(wp) uatemp(jpi,jpj,jpk) , vatemp(jpi,jpj,jpk)
      INTEGER :: ji,jj,jk
      INTEGER kt
      REAL(wp) :: z2dt, znugdt
      REAL(wp), DIMENSION(jpi,jpj) :: uatemp2D, vatemp2D
      REAL(wp) :: timeref
      REAL(wp), DIMENSION(jpi,jpj) :: spgu1,spgv1
      REAL(wp) :: rhox, rhoy

      IF (Agrif_Root()) RETURN

      rhox = Agrif_Rhox()
      rhoy = Agrif_Rhoy()

      timeref = 1.

      ! time step: leap-frog
      z2dt = 2. * rdt
      ! time step: Euler if restart from rest
      IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt
      ! coefficients
      znugdt =  rnu * grav * z2dt    

        Agrif_SpecialValue=0.
        Agrif_UseSpecialValue = .TRUE.
        uatemp = 0.
        vatemp = 0.
        Call Agrif_Bc_variable(uatemp,un,procname=interpu)
        Call Agrif_Bc_variable(vatemp,vn,procname=interpv)
        uatemp2d = 0.
        vatemp2d = 0.

          Agrif_SpecialValue=0.
        Agrif_UseSpecialValue = .TRUE.
       Call Agrif_Bc_variable(uatemp2d,e1u,calledweight=1.,procname=interpu2d)
       Call Agrif_Bc_variable(vatemp2d,e2v,calledweight=1.,procname=interpv2d)
        Agrif_UseSpecialValue = .FALSE.


        If ((nbondi == -1).OR.(nbondi == 2)) THEN

        DO jj=1,jpj
          laplacu(2,jj) = timeref * (uatemp2d(2,jj)/(rhoy*e2u(2,jj)))*umask(2,jj,1)
        ENDDO
        
        Do jk=1,jpkm1
        DO jj=1,jpj
          ua(1:2,jj,jk) = (uatemp(1:2,jj,jk)/(rhoy*e2u(1:2,jj)))
#if defined key_partial_steps
           ua(1:2,jj,jk) = ua(1:2,jj,jk) / fse3u(1:2,jj,jk)
#endif
        ENDDO
        ENDDO

        Do jk=1,jpkm1
        DO jj=1,jpj
          ua(2,jj,jk) = (ua(2,jj,jk) - z2dt * znugdt * laplacu(2,jj))*umask(2,jj,jk)
        ENDDO
        ENDDO

        spgu(2,:)=0.

        do jk=1,jpkm1
        do jj=1,jpj
        spgu(2,jj)=spgu(2,jj)+fse3u(2,jj,jk)*ua(2,jj,jk)
        enddo
        enddo

        DO jj=1,jpj
        IF (umask(2,jj,1).NE.0.) THEN
         spgu(2,jj)=spgu(2,jj)/hu(2,jj)
        ENDIF
        enddo

        Do jk=1,jpkm1
        DO jj=1,jpj
          ua(2,jj,jk) = 0.25*(ua(1,jj,jk)+2.*ua(2,jj,jk)+ua(3,jj,jk))
          ua(2,jj,jk) = ua(2,jj,jk) * umask(2,jj,jk)
        ENDDO
        ENDDO

        spgu1(2,:)=0.

        do jk=1,jpkm1
        do jj=1,jpj
        spgu1(2,jj)=spgu1(2,jj)+fse3u(2,jj,jk)*ua(2,jj,jk)
        enddo
        enddo

        DO jj=1,jpj
        IF (umask(2,jj,1).NE.0.) THEN
         spgu1(2,jj)=spgu1(2,jj)/hu(2,jj)
        ENDIF
        enddo

        DO jk=1,jpkm1
        DO jj=1,jpj
         ua(2,jj,jk) = (ua(2,jj,jk)+spgu(2,jj)-spgu1(2,jj))*umask(2,jj,jk)
        ENDDO
        ENDDO

        Do jk=1,jpkm1
        Do jj=1,jpj
           va(2,jj,jk) = (vatemp(2,jj,jk)/(rhox*e1v(2,jj)))*vmask(2,jj,jk)
#if defined key_partial_steps
           va(2,jj,jk) = va(2,jj,jk) / fse3v(2,jj,jk)
#endif           
        End Do
        End Do

        sshn(2,:)=sshn(3,:)
        sshb(2,:)=sshb(3,:)
                                
        ENDIF

        If ((nbondi == 1).OR.(nbondi == 2)) THEN

        DO jj=1,jpj
          laplacu(nlci-2,jj) = timeref * (uatemp2d(nlci-2,jj)/(rhoy*e2u(nlci-2,jj)))
        ENDDO

        Do jk=1,jpkm1
        DO jj=1,jpj
          ua(nlci-2:nlci-1,jj,jk) = (uatemp(nlci-2:nlci-1,jj,jk)/(rhoy*e2u(nlci-2:nlci-1,jj)))

#if defined key_partial_steps
           ua(nlci-2:nlci-1,jj,jk) = ua(nlci-2:nlci-1,jj,jk) / fse3u(nlci-2:nlci-1,jj,jk)
#endif

        ENDDO
        ENDDO

        Do jk=1,jpkm1
        DO jj=1,jpj
          ua(nlci-2,jj,jk) = (ua(nlci-2,jj,jk)- z2dt * znugdt * laplacu(nlci-2,jj))*umask(nlci-2,jj,jk)
        ENDDO
        ENDDO


        spgu(nlci-2,:)=0.

        do jk=1,jpkm1
        do jj=1,jpj
        spgu(nlci-2,jj)=spgu(nlci-2,jj)+fse3u(nlci-2,jj,jk)*ua(nlci-2,jj,jk)
        enddo
        enddo

        DO jj=1,jpj
        IF (umask(nlci-2,jj,1).NE.0.) THEN
         spgu(nlci-2,jj)=spgu(nlci-2,jj)/hu(nlci-2,jj)
        ENDIF
        enddo

        Do jk=1,jpkm1
        DO jj=1,jpj
         ua(nlci-2,jj,jk) = 0.25*(ua(nlci-3,jj,jk)+2.*ua(nlci-2,jj,jk)+ua(nlci-1,jj,jk))

          ua(nlci-2,jj,jk) = ua(nlci-2,jj,jk) * umask(nlci-2,jj,jk)

        ENDDO
        ENDDO

        spgu1(nlci-2,:)=0.

        do jk=1,jpkm1
        do jj=1,jpj
        spgu1(nlci-2,jj)=spgu1(nlci-2,jj)+fse3u(nlci-2,jj,jk)*ua(nlci-2,jj,jk)*umask(nlci-2,jj,jk)
        enddo
        enddo

        DO jj=1,jpj
        IF (umask(nlci-2,jj,1).NE.0.) THEN
         spgu1(nlci-2,jj)=spgu1(nlci-2,jj)/hu(nlci-2,jj)
        ENDIF
        enddo

        DO jk=1,jpkm1
        DO jj=1,jpj
         ua(nlci-2,jj,jk) = (ua(nlci-2,jj,jk)+spgu(nlci-2,jj)-spgu1(nlci-2,jj))*umask(nlci-2,jj,jk)
        ENDDO
        ENDDO

        Do jk=1,jpkm1
        Do jj=1,jpj-1
           va(nlci-1,jj,jk) = (vatemp(nlci-1,jj,jk)/(rhox*e1v(nlci-1,jj)))*vmask(nlci-1,jj,jk)
#if defined key_partial_steps
           va(nlci-1,jj,jk) = va(nlci-1,jj,jk) / fse3v(nlci-1,jj,jk)
#endif
        End Do
        End Do

        sshn(nlci-1,:)=sshn(nlci-2,:)
        sshb(nlci-1,:)=sshb(nlci-2,:)        
        ENDIF

        If ((nbondj == -1).OR.(nbondj == 2)) THEN

        DO ji=1,jpi
          laplacv(ji,2) = timeref * (vatemp2d(ji,2)/(rhox*e1v(ji,2)))
        ENDDO

        DO jk=1,jpkm1
        DO ji=1,jpi
          va(ji,1:2,jk) = (vatemp(ji,1:2,jk)/(rhox*e1v(ji,1:2)))
#if defined key_partial_steps
           va(ji,1:2,jk) = va(ji,1:2,jk) / fse3v(ji,1:2,jk)
#endif
        ENDDO
        ENDDO

        DO jk=1,jpkm1
        DO ji=1,jpi
          va(ji,2,jk) = (va(ji,2,jk) - z2dt * znugdt * laplacv(ji,2))*vmask(ji,2,jk)
        ENDDO
        ENDDO

        spgv(:,2)=0.

        do jk=1,jpkm1
        do ji=1,jpi
        spgv(ji,2)=spgv(ji,2)+fse3v(ji,2,jk)*va(ji,2,jk)
        enddo
        enddo

        DO ji=1,jpi
        IF (vmask(ji,2,1).NE.0.) THEN
         spgv(ji,2)=spgv(ji,2)/hv(ji,2)
        ENDIF
        enddo

        DO jk=1,jpkm1
        DO ji=1,jpi
          va(ji,2,jk)=0.25*(va(ji,1,jk)+2.*va(ji,2,jk)+va(ji,3,jk))
           va(ji,2,jk)=va(ji,2,jk)*vmask(ji,2,jk)
        ENDDO
        ENDDO

        spgv1(:,2)=0.

        do jk=1,jpkm1
        do ji=1,jpi
        spgv1(ji,2)=spgv1(ji,2)+fse3v(ji,2,jk)*va(ji,2,jk)*vmask(ji,2,jk)
        enddo
        enddo

        DO ji=1,jpi
        IF (vmask(ji,2,1).NE.0.) THEN
         spgv1(ji,2)=spgv1(ji,2)/hv(ji,2)
        ENDIF
        enddo

        DO jk=1,jpkm1
        DO ji=1,jpi
         va(ji,2,jk) = (va(ji,2,jk)+spgv(ji,2)-spgv1(ji,2))*vmask(ji,2,jk)
        ENDDO
        ENDDO

        DO jk=1,jpkm1
        DO ji=1,jpi
        ua(ji,2,jk) = (uatemp(ji,2,jk)/(rhoy*e2u(ji,2)))*umask(ji,2,jk) 
#if defined key_partial_steps
           ua(ji,2,jk) = ua(ji,2,jk) / fse3u(ji,2,jk)
#endif                
        ENDDO
        ENDDO

        sshn(:,2)=sshn(:,3)
        sshb(:,2)=sshb(:,3)
        ENDIF

        If ((nbondj == 1).OR.(nbondj == 2)) THEN

        DO ji=1,jpi
          laplacv(ji,nlcj-2) = timeref * (vatemp2d(ji,nlcj-2)/(rhox*e1v(ji,nlcj-2)))
        ENDDO

        DO jk=1,jpkm1
        DO ji=1,jpi
          va(ji,nlcj-2:nlcj-1,jk) = (vatemp(ji,nlcj-2:nlcj-1,jk)/(rhox*e1v(ji,nlcj-2:nlcj-1)))
#if defined key_partial_steps
           va(ji,nlcj-2:nlcj-1,jk) = va(ji,nlcj-2:nlcj-1,jk) / fse3v(ji,nlcj-2:nlcj-1,jk)
#endif
        ENDDO
        ENDDO

        DO jk=1,jpkm1
        DO ji=1,jpi
          va(ji,nlcj-2,jk) = (va(ji,nlcj-2,jk)-z2dt * znugdt * laplacv(ji,nlcj-2))*vmask(ji,nlcj-2,jk)
        ENDDO
        ENDDO


        spgv(:,nlcj-2)=0.

        do jk=1,jpkm1
        do ji=1,jpi
        spgv(ji,nlcj-2)=spgv(ji,nlcj-2)+fse3v(ji,nlcj-2,jk)*va(ji,nlcj-2,jk)
        enddo
        enddo

        DO ji=1,jpi
        IF (vmask(ji,nlcj-2,1).NE.0.) THEN
         spgv(ji,nlcj-2)=spgv(ji,nlcj-2)/hv(ji,nlcj-2)
        ENDIF
        enddo

        DO jk=1,jpkm1
        DO ji=1,jpi
           va(ji,nlcj-2,jk)=0.25*(va(ji,nlcj-3,jk)+2.*va(ji,nlcj-2,jk)+va(ji,nlcj-1,jk))
           va(ji,nlcj-2,jk) = va(ji,nlcj-2,jk) * vmask(ji,nlcj-2,jk)
        ENDDO
        ENDDO

        spgv1(:,nlcj-2)=0.

        do jk=1,jpkm1
        do ji=1,jpi
        spgv1(ji,nlcj-2)=spgv1(ji,nlcj-2)+fse3v(ji,nlcj-2,jk)*va(ji,nlcj-2,jk)
        enddo
        enddo

        DO ji=1,jpi
        IF (vmask(ji,nlcj-2,1).NE.0.) THEN
         spgv1(ji,nlcj-2)=spgv1(ji,nlcj-2)/hv(ji,nlcj-2)
        ENDIF
        enddo

        DO jk=1,jpkm1
        DO ji=1,jpi
        va(ji,nlcj-2,jk) = (va(ji,nlcj-2,jk)+spgv(ji,nlcj-2)-spgv1(ji,nlcj-2))*vmask(ji,nlcj-2,jk)
        ENDDO
        ENDDO

        DO jk=1,jpkm1
        DO ji=1,jpi
          ua(ji,nlcj-1,jk) = (uatemp(ji,nlcj-1,jk)/(rhoy*e2u(ji,nlcj-1)))*umask(ji,nlcj-1,jk)
#if defined key_partial_steps
           ua(ji,nlcj-1,jk) = ua(ji,nlcj-1,jk) / fse3u(ji,nlcj-1,jk)
#endif          
        ENDDO
        ENDDO
        
        sshn(:,nlcj-1)=sshn(:,nlcj-2)
        sshb(:,nlcj-1)=sshb(:,nlcj-2)                
        ENDIF
            
!
      Return
      End Subroutine Agrif_dyn


       subroutine interpu(tabres,i1,i2,j1,j2,k1,k2)
       Implicit none
#  include "domzgr_substitute.h90"       
       integer i1,i2,j1,j2,k1,k2
       integer ji,jj,jk
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres

       do jk=k1,k2
       DO jj=j1,j2
       DO ji=i1,i2
         tabres(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
#if defined key_partial_steps
          tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3u(ji,jj,jk)
#endif
       ENDDO
       ENDDO
       ENDDO
       end subroutine interpu

       subroutine interpu2d(tabres,i1,i2,j1,j2)
       Implicit none
       integer i1,i2,j1,j2
       integer ji,jj
       real,dimension(i1:i2,j1:j2) :: tabres

       DO jj=j1,j2
       DO ji=i1,i2
         tabres(ji,jj) = e2u(ji,jj) * ((gcx(ji+1,jj) - gcx(ji,jj))/e1u(ji,jj)) &
                                       *umask(ji,jj,1)
       ENDDO
       ENDDO
       end subroutine interpu2d

       subroutine interpv(tabres,i1,i2,j1,j2,k1,k2)
       Implicit none
#  include "domzgr_substitute.h90"       
       integer i1,i2,j1,j2,k1,k2
       integer ji,jj,jk
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres

       do jk=k1,k2
       DO jj=j1,j2
       DO ji=i1,i2
         tabres(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
#if defined key_partial_steps
          tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3v(ji,jj,jk)
#endif           
       ENDDO
       ENDDO
       ENDDO
       end subroutine interpv

       subroutine interpv2d(tabres,i1,i2,j1,j2)
       Implicit none
       integer i1,i2,j1,j2
       integer ji,jj
       real,dimension(i1:i2,j1:j2) :: tabres

       DO jj=j1,j2
       DO ji=i1,i2
         tabres(ji,jj) = e1v(ji,jj) * ((gcx(ji,jj+1) - gcx(ji,jj))/e2v(ji,jj)) &
                                       * vmask(ji,jj,1)
       ENDDO
       ENDDO
       end subroutine interpv2d

#else
      CONTAINS
      subroutine Agrif_OPA_Interp_empty

      end subroutine Agrif_OPA_Interp_empty
#endif
      End Module agrif_opa_interp

