#define TWO_WAY

      Module agrif_opa_update
#if defined key_agrif
      USE par_oce
      USE oce
      USE dom_oce
      
      Integer, Parameter :: nbclineupdate = 3
      Integer :: nbcline

      Contains

      Subroutine Agrif_Update_Tra( kt )
!
!     Modules used:
!

      implicit none
!
!     Declarations:
      INTEGER :: kt
!
!
!     Variables
!
      Real :: tabtemp(jpi,jpj,jpk)
!
!     Begin
!

      IF ((Agrif_NbStepint() .NE. (Agrif_irhot()-1)).AND.(kt /= 0)) Return
#if defined TWO_WAY
      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid = 0.
      IF (mod(nbcline,nbclineupdate) == 0) THEN
      Call Agrif_Update_Variable(tabtemp,tn, procname=updateT)
      Call Agrif_Update_Variable(tabtemp,sn, procname=updateS)
      ELSE
      Call Agrif_Update_Variable(tabtemp,tn,locupdate=(/0,2/), procname=updateT)
      Call Agrif_Update_Variable(tabtemp,sn,locupdate=(/0,2/), procname=updateS)
      ENDIF


      Agrif_UseSpecialValueInUpdate = .FALSE.
#endif

      Return
      End subroutine Agrif_Update_Tra

      Subroutine Agrif_Update_Dyn( kt )
!
!     Modules used:
!
!
!     Declarations:
!
      INTEGER :: kt
!
!     Variables
!
      Real :: tabtemp(jpi,jpj,jpk)
      Real :: tabtemp2d(jpi,jpj)
!
!     Begin
!
!
       return
       
      IF ((Agrif_NbStepint() .NE. (Agrif_irhot()-1)).AND.(kt /= 0)) Return
#if defined TWO_WAY

      IF (mod(nbcline,nbclineupdate) == 0) THEN
      Call Agrif_Update_Variable(tabtemp,un,procname = updateU)
      Call Agrif_Update_Variable(tabtemp,vn,procname = updateV)
      ELSE
      Call Agrif_Update_Variable(tabtemp,un,locupdate=(/0,1/),procname = updateU)
      Call Agrif_Update_Variable(tabtemp,vn,locupdate=(/0,1/),procname = updateV)         
      ENDIF

      Call Agrif_Update_Variable(tabtemp2d,e1u,procname = updateU2d)
      Call Agrif_Update_Variable(tabtemp2d,e2v,procname = updateV2d)  
      
      nbcline = nbcline + 1

       Agrif_UseSpecialValueInUpdate = .TRUE.
       Agrif_SpecialValueFineGrid = 0.
       Call Agrif_Update_Variable(tabtemp2d,sshn,procname = updateSSH)
       Agrif_UseSpecialValueInUpdate = .FALSE.


      Call Agrif_ChildGrid_To_ParentGrid()
      Call recompute_diags( kt )
      Call Agrif_ParentGrid_To_ChildGrid()

#endif
!
      Return
      End subroutine Agrif_Update_Dyn

      Subroutine recompute_diags(kt)
      Use divcur
      Use wzvmod
      Use cla_div
      Use  ocfzpt
      Implicit None
      INTEGER kt
      
      ta = hdivb
      sa = rotb
      CALL oc_fz_pt
      Call div_cur(kt)

      hdivb = ta
      rotb  = sa

      IF( n_cla == 1 ) CALL div_cla( kt )
      Call wzv( kt )
      
      End Subroutine recompute_diags

       subroutine updateT(tabres,i1,i2,j1,j2,k1,k2,before)
       Implicit none
#  include "domzgr_substitute.h90"
       integer i1,i2,j1,j2,k1,k2
       integer ji,jj,jk
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres
       LOGICAL :: before

       IF (before) THEN
       
         DO jk=k1,k2
           DO jj=j1,j2
             DO ji=i1,i2
               tabres(ji,jj,jk) = tn(ji,jj,jk)
             ENDDO
           ENDDO
         ENDDO
         
       ELSE

         DO jk=k1,k2
           DO jj=j1,j2
             DO ji=i1,i2
               IF (tabres(ji,jj,jk).NE.0.) THEN
               tn(ji,jj,jk) = tabres(ji,jj,jk) * tmask(ji,jj,jk)
               ENDIF
             ENDDO
            ENDDO
          ENDDO
       ENDIF

       end subroutine updateT

       
       subroutine updateS(tabres,i1,i2,j1,j2,k1,k2,before)
       Implicit none
#  include "domzgr_substitute.h90"
       integer i1,i2,j1,j2,k1,k2
       integer ji,jj,jk
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres
       LOGICAL :: before


       IF (before) THEN
       
         DO jk=k1,k2
           DO jj=j1,j2
             DO ji=i1,i2
               tabres(ji,jj,jk) = sn(ji,jj,jk)
             ENDDO
           ENDDO
         ENDDO
         
       ELSE

         DO jk=k1,k2
           DO jj=j1,j2
             DO ji=i1,i2
               IF (tabres(ji,jj,jk).NE.0.) THEN
               sn(ji,jj,jk) = tabres(ji,jj,jk) * tmask(ji,jj,jk)
               ENDIF
             ENDDO
           ENDDO
         ENDDO
       ENDIF

       end subroutine updateS

       subroutine updateu(tabres,i1,i2,j1,j2,k1,k2,before)
       Implicit none
#  include "domzgr_substitute.h90"
       integer i1,i2,j1,j2,k1,k2
       integer ji,jj,jk
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres
       LOGICAL :: before
       REAL(wp) :: rhoy


       IF (before) THEN
       
       rhoy = Agrif_Rhoy()
       
         DO jk=k1,k2
           DO jj=j1,j2
             DO ji=i1,i2
               tabres(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
#if defined key_partial_steps
               tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3u(ji,jj,jk)
#endif
             ENDDO
           ENDDO
         ENDDO
 
         tabres = rhoy * tabres
 
       ELSE

         DO jk=k1,k2
           DO jj=j1,j2
             DO ji=i1,i2
               un(ji,jj,jk) = tabres(ji,jj,jk) / (e2u(ji,jj))
               un(ji,jj,jk) = un(ji,jj,jk) * umask(ji,jj,jk)
#if defined key_partial_steps
               un(ji,jj,jk) = un(ji,jj,jk) / fse3u(ji,jj,jk)
#endif
       ENDDO
       ENDDO
       ENDDO
       ENDIF

       end subroutine updateu

       subroutine updatev(tabres,i1,i2,j1,j2,k1,k2,before)
       Implicit none
#  include "domzgr_substitute.h90"
       integer i1,i2,j1,j2,k1,k2
       integer ji,jj,jk
       real,dimension(i1:i2,j1:j2,k1:k2) :: tabres
       LOGICAL :: before
       REAL(wp) :: rhox


       IF (before) THEN
       
       rhox = Agrif_Rhox()
       
         DO jk=k1,k2
           DO jj=j1,j2
             DO ji=i1,i2
               tabres(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
#if defined key_partial_steps
               tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3v(ji,jj,jk)
#endif
             ENDDO
           ENDDO
         ENDDO
 
        tabres = rhox * tabres
 
       ELSE

         DO jk=k1,k2
           DO jj=j1,j2
             DO ji=i1,i2
               vn(ji,jj,jk) = tabres(ji,jj,jk) / (e1v(ji,jj))
               vn(ji,jj,jk) = vn(ji,jj,jk) * vmask(ji,jj,jk)
#if defined key_partial_steps
               vn(ji,jj,jk) = vn(ji,jj,jk) / fse3v(ji,jj,jk)
#endif
       ENDDO
       ENDDO
       ENDDO
       ENDIF

       end subroutine updatev

       subroutine updateu2d(tabres,i1,i2,j1,j2,before)
       Implicit none
#  include "domzgr_substitute.h90"
       integer i1,i2,j1,j2
       integer ji,jj,jk
       real,dimension(i1:i2,j1:j2) :: tabres
       LOGICAL :: before
       REAL(wp) :: rhoy
       REAL(wp) :: hinv


       IF (before) THEN
       
       rhoy = Agrif_Rhoy()
       
           DO jk = 1,jpkm1
             DO jj=j1,j2
             DO ji=i1,i2
                tabres(ji,jj) = tabres(ji,jj) + fse3u(ji,jj,jk) * un(ji,jj,jk)
             ENDDO
             ENDDO
           ENDDO
           
           DO jj=j1,j2
           DO ji=i1,i2
             tabres(ji,jj) = tabres(ji,jj) * e2u(ji,jj)
           ENDDO
           ENDDO
   
          tabres = rhoy * tabres
   
       ELSE

           DO jj=j1,j2
             DO ji=i1,i2
               IF (umask(ji,jj,1) .NE. 0.) THEN             
               spgu(ji,jj) = 0.
               Do jk=1,jpk
                spgu(ji,jj) = spgu(ji,jj) + fse3u(ji,jj,jk) * un(ji,jj,jk)
               EndDo
               spgu(ji,jj) = spgu(ji,jj) * e2u(ji,jj)
               hinv = (tabres(ji,jj)-spgu(ji,jj))/(hu(ji,jj)*e2u(ji,jj))
               Do jk=1,jpk              
               un(ji,jj,jk) = un(ji,jj,jk) + hinv
               un(ji,jj,jk) = un(ji,jj,jk) * umask(ji,jj,jk)            
               EndDo
               ENDIF
             ENDDO
           ENDDO
       ENDIF

       end subroutine updateu2d

       subroutine updatev2d(tabres,i1,i2,j1,j2,before)
       Implicit none
       integer i1,i2,j1,j2
       integer ji,jj,jk
       real,dimension(i1:i2,j1:j2) :: tabres
       LOGICAL :: before
       REAL(wp) :: rhox
       REAL(wp) :: hinv


       IF (before) THEN
       
       rhox = Agrif_Rhox()
       
           tabres = 0.
           
           DO jk = 1,jpkm1
             DO jj=j1,j2
             DO ji=i1,i2
                tabres(ji,jj) = tabres(ji,jj) + fse3v(ji,jj,jk) * vn(ji,jj,jk)
             ENDDO
             ENDDO
           ENDDO
           
           DO jj=j1,j2
           DO ji=i1,i2
              tabres(ji,jj) = tabres(ji,jj) * e1v(ji,jj)
           ENDDO
           ENDDO
   
         tabres = rhox * tabres
   
       ELSE

           DO jj=j1,j2
             DO ji=i1,i2
               IF (vmask(ji,jj,1) .NE. 0.) THEN             
               spgv(ji,jj) = 0.
               Do jk=1,jpk
                spgv(ji,jj) = spgv(ji,jj) + fse3v(ji,jj,jk) * vn(ji,jj,jk)
               EndDo
               spgv(ji,jj) = spgv(ji,jj) * e1v(ji,jj)
               hinv = (tabres(ji,jj)-spgv(ji,jj))/(hv(ji,jj)*e1v(ji,jj))

               Do jk=1,jpk             
               vn(ji,jj,jk) = vn(ji,jj,jk) + hinv
               vn(ji,jj,jk) = vn(ji,jj,jk) * vmask(ji,jj,jk)
               EndDo
               ENDIF
           ENDDO
           ENDDO
           
       ENDIF

       end subroutine updatev2d

       subroutine updateSSH(tabres,i1,i2,j1,j2,before)
       Implicit none
#  include "domzgr_substitute.h90"
       integer i1,i2,j1,j2
       integer ji,jj
       real,dimension(i1:i2,j1:j2) :: tabres
       LOGICAL :: before
       REAL(wp) :: rhox, rhoy


       IF (before) THEN
       rhox = Agrif_Rhox()
       rhoy = Agrif_Rhoy()
       
           DO jj=j1,j2
             DO ji=i1,i2
               tabres(ji,jj) = e1t(ji,jj) * e2t(ji,jj) * sshn(ji,jj)
             ENDDO
           ENDDO
   
         tabres = rhox * rhoy * tabres
 
       ELSE
           DO jj=j1,j2
             DO ji=i1,i2
               sshn(ji,jj) = tabres(ji,jj) / (e1t(ji,jj) * e2t(ji,jj))
               sshn(ji,jj) = sshn(ji,jj) * tmask(ji,jj,1)
       ENDDO
       ENDDO
       ENDIF

       end subroutine updateSSH
       
#else
       CONTAINS
       subroutine agrif_opa_update_empty
       end subroutine agrif_opa_update_empty
#endif
       End Module agrif_opa_update
