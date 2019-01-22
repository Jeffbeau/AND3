MODULE trcrst
   !!======================================================================
   !!
   !!                       *** MODULE trcrst ***
   !!
   !!   Read the restart files for passive tracers
   !!
   !!======================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trcrst.F90,v 1.7 2006/04/10 15:40:28 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_passivetrc   
   !!----------------------------------------------------------------------
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   USE sms
   USE trctrp_lec   
   USE lib_mpp
   
   IMPLICIT NONE
   PRIVATE
   
   !! * Accessibility
   PUBLIC trc_rst
   PUBLIC trc_wri
   
   !! * Module variables
   CHARACTER (len=48) ::   &
      trestart = 'initial.trc.nc'   ! restart file name

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   
CONTAINS

#if defined key_fdir
   !!----------------------------------------------------------------------
   !!   'key_fdir'                                       direct access file
   !!----------------------------------------------------------------------
#include "trcrst_fdir.h90"
   
#else 

   SUBROUTINE trc_rst 
      !!===========================================================================================
      !!
      !!                       ROUTINE trc_rst
      !!                     *******************
      !!
      !!  PURPOSE :
      !!  ---------
      !!     READ files for restart for passive tracer
      !!
      !!   METHOD :
      !!   -------
      !!      READ the previous fields on the FILE nutrst
      !!      the first record indicates previous characterics
      !!      after control with the present run, we READ :
      !!      - prognostic variables on the second and more record
      !!
      !!   History:
      !!   --------
      !!  original  : 96-11
      !!  00-05 (A. Estublier) TVD Limiter Scheme key_trc_tvd 
      !!  00-12 (O. Aumont, E. Kestenare) read restart file for sediments 
      !!  01-05 (O. Aumont, E. Kestenare) read restart file for calcite and silicate sediments 
      !!  05-03 (O. Aumont and A. El Moussaoui) F90           
      !!------------------------------------------------------------------------
      !! * Modules used
      USE ioipsl


      !! local declarations
      !! ==================
      LOGICAL ::  llog       !!!
      CHARACTER (len=32) :: clname1,clname2
      CHARACTER (len=32) :: clname = 'restart.trc'
      CHARACTER (len=12) :: clvnames(80)  

      INTEGER :: ino1,jn,iarak0,iarak1,          &
         ji, jj, jk,                   &
         itime, ibvar
      REAL(wp) :: caralk,bicarb,zdt,        &     
         zdate0
      REAL(wp) ::   zdept(jpk), zlamt(jpi,jpj), zphit(jpi,jpj)

      REAL(wp), DIMENSION(3) :: zinfo

#if defined key_trc_pisces 
#if ! defined key_cfg_1d && ( defined key_orca_r4 || defined key_orca_r2 || defined key_orca_r05 || defined key_orca_r025 )
      REAL(wp) ::   zareatot, zpo4tot
#endif
#endif

      !!---------------------------------------------------------------------
      !!  OPA.9 03-2005  
      !!---------------------------------------------------------------------
      !! 0. initialisations
      !!------------------


      IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
         iarak0 = 1
      ELSE
         iarak0=0
      ENDIF


      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) ' *** trc_rst beginning of restart for'
      IF(lwp) WRITE(numout,*) ' passive tracer'
      IF(lwp) WRITE(numout,*) ' the present run :'
      IF(lwp) WRITE(numout,*) '   number job is  : ',no
      IF(lwp) WRITE(numout,*) '   with the time nit000 : ',nit000
      IF(lwp) THEN
         IF(iarak0.eq.1) then
            WRITE(numout,*) '   and before fields for Arakawa sheme '
         ENDIF
         WRITE(numout,*) ' '
      ENDIF

      ! Time domain : restart
      ! -------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' *** passive tracer restart option'
      SELECT CASE ( nrsttr )
      CASE ( 0 )
         IF(lwp) WRITE(numout,*) ' nrsttr = 0 no control of nit000'
      CASE ( 1 )
         IF(lwp) WRITE(numout,*) ' nrsttr = 1 we control the date of nit000'
      CASE ( 2 )
         IF(lwp) WRITE(numout,*) ' nrsttr = 2 the date adatrj is read in restart file'
      CASE DEFAULT
         IF(lwp) WRITE(numout,*) '  ===>>>> nrsttr not equal 0, 1 or 2 : no control of the date'
         IF(lwp) WRITE(numout,*) ' =======                   ========='
      END SELECT


      !! 1. READ nutrst
      !! --------------
      !! ... first information
      !! ---------------------
      itime=0
      llog=.false.           !!!
      zlamt(:,:) = 0.e0
      zphit(:,:) = 0.e0
      zdept(:)   = 0.e0
      CALL restini(clname,jpi,jpj,zlamt,zphit,jpk,zdept,clname         & 
         &           ,itime,zdate0,zdt,nutrst,domain_id=nidom)

      CALL ioget_vname(nutrst, ibvar, clvnames)
      CALL restget(nutrst,'info',1,1,3,0,llog,zinfo)
      ino1  = nint(zinfo(1))
      iarak1 = nint(zinfo(3))

      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) ' READ nutrst with '
      IF(lwp) WRITE(numout,*) '   number job is  : ',ino1
      IF(lwp) WRITE(numout,*) '   with the time it : ',nint(zinfo(2))
      IF(lwp) THEN
         IF(iarak1.eq.1) then
            WRITE(numout,*) '   and before fields for Arakawa sheme '
         ENDIF
      ENDIF
      IF(lwp) WRITE(numout,*) '   number of variables   : ', ibvar
      IF(lwp) WRITE(numout,*) '   NetCDF variables      : '
      IF(lwp) WRITE(numout,*) ' ',clvnames (:ibvar)
      IF(lwp) WRITE(numout,*) ' '

      !! 1.2 control of date
      !! -------------------

      IF( nit000- NINT( zinfo(2) ) /= 1 .AND. nrsttr /= 0 ) THEN
         IF(lwp) THEN
            WRITE(numout,*) ' ===>>>> : problem with nit000 for the',    &  
               ' passive tracer restart'
            WRITE(numout,*) ' =======                              ',    &    
               ' ======================'
            WRITE(numout,*) ' we stop. verify the FILE'
            WRITE(numout,*) ' or rerun with the value  0 for the'
            WRITE(numout,*) ' control of time PARAMETER   nrstdt'
            WRITE(numout,*) ' '
         ENDIF
         STOP 'trc_rst'       !!
      ENDIF

      !! 1.3 Control of the sheme
      !! ------------------------

      IF(iarak0.ne.iarak1) THEN
         IF(lwp) THEN
            WRITE(numout,*) ' ===>>>> : problem with the',       &   
               ' passive tracer restart file'
            WRITE(numout,*) ' =======                              ',        & 
               ' ==========================='
            WRITE(numout,*) ' we stop. verify the FILE'
            WRITE(numout,*) ' before field required IF 1=',iarak0
            WRITE(numout,*) ' before field present in file IF 1=',           & 
               iarak1
            WRITE(numout,*) ' '
         ENDIF
         STOP 'trc_rst'       !!!!!    AVERIFIER AU NIV F90'
      ENDIF


      !! ... READ prognostic variables and computes diagnostic variable
      !! ---------------------------------------------------------------

      DO jn=1,jptra
         clname='TRN'//ctrcnm(jn)
         CALL restget(nutrst,clname,jpi,jpj,jpk,0,llog,trn(:,:,:,jn))
      END DO

      DO jn=1,jptra
         clname='TRB'//ctrcnm(jn)
         CALL restget(nutrst,clname,jpi,jpj,jpk,0,llog,trb(:,:,:,jn))
      END DO


#if defined key_trc_lobster1
      clname='SEDB'//ctrcnm(jpdet)
      clname1='SEDN'//ctrcnm(jpdet)
      CALL restget(nutrst,clname,jpi,jpj,1,0,llog,sedpocb(:,:))
      CALL restget(nutrst,clname1,jpi,jpj,1,0,llog,sedpocn(:,:))
#elif defined key_trc_pisces
      clname='Silicalim'
      CALL restget(nutrst,clname,jpi,jpj,1,0,llog,xksi)
      xksimax=xksi

      clname='SED'//ctrcnm(jppoc)
      clname1='SED'//ctrcnm(jpcal)
      clname2='SED'//ctrcnm(jpsil)
      CALL restget(nutrst,clname1,jpi,jpj,1,0,llog,sedcal(:,:))
      CALL restget(nutrst,clname2,jpi,jpj,1,0,llog,sedsil(:,:))
      CALL restget(nutrst,clname,jpi,jpj,1,0,llog,sedpoc(:,:))

#elif defined key_cfc
      clname='qint'
      CALL restget(nutrst,clname,jpi,jpj,jptra,0,llog,qint(:,:,:))
      clname1='qtr'
      CALL restget(nutrst,clname1,jpi,jpj,jptra,0,llog,qtr(:,:,:))         
#endif

#if defined key_trc_pisces 

#if ! defined key_cfg_1d && ( defined key_orca_r4 || defined key_orca_r2 || defined key_orca_r05 || defined key_orca_r025 ) 

      zareatot = 0.
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zareatot = zareatot + tmask(ji,jj,jk) * tmask_i(ji,jj) * &
                  &                 e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) 
            END DO
         END DO
      END DO

      IF( lk_mpp ) THEN 
         CALL mpp_sum( zareatot )     ! sum over the global domain  
      END IF

      zpo4tot = 0.
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zpo4tot = zpo4tot + trn(ji,jj,jk,jptal) * tmask(ji,jj,jk) * tmask_i(ji,jj) *   &
                  &                e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( lk_mpp ) THEN 
         CALL mpp_sum( zpo4tot )     ! sum over the global domain  
      END IF

      WRITE(0,*) 'TALK moyen ', zpo4tot/zareatot*1E6
      zpo4tot = zpo4tot/zareatot*1E6
      trn(:,:,:,jptal) = trn(:,:,:,jptal)*2391./zpo4tot

      zpo4tot = 0.
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zpo4tot = zpo4tot + trn(ji,jj,jk,jppo4) * tmask(ji,jj,jk) * tmask_i(ji,jj) *   &
                  &                e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( lk_mpp ) THEN 
         CALL mpp_sum( zpo4tot )     ! sum over the global domain  
      END IF


      WRITE(0,*) 'PO4 moyen ', zpo4tot/zareatot*1E6/122.
      zpo4tot = zpo4tot/zareatot*1E6/122.
      trn(:,:,:,jppo4) = trn(:,:,:,jppo4)*2.165/zpo4tot

      zpo4tot = 0.
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zpo4tot = zpo4tot + trn(ji,jj,jk,jpno3) * tmask(ji,jj,jk) * tmask_i(ji,jj) *   &
                  &                e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( lk_mpp ) THEN 
         CALL mpp_sum( zpo4tot )     ! sum over the global domain  
      END IF


      WRITE(0,*) 'NO3 moyen ', zpo4tot/zareatot*1E6/7.6
      zpo4tot = zpo4tot/zareatot*1E6/7.6
      trn(:,:,:,jpno3) = trn(:,:,:,jpno3)*30.9/zpo4tot

      zpo4tot = 0.
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zpo4tot = zpo4tot + trn(ji,jj,jk,jpsil) * tmask(ji,jj,jk) * tmask_i(ji,jj) *   &
                  &                e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( lk_mpp ) THEN 
         CALL mpp_sum( zpo4tot )     ! sum over the global domain  
      END IF

      WRITE(0,*) 'SiO3 moyen ', zpo4tot/zareatot*1E6
      zpo4tot = zpo4tot/zareatot*1E6
      trn(:,:,:,jpsil) = MIN( 400E-6,trn(:,:,:,jpsil)*91.51/zpo4tot) 

#endif
      !!  Initialization of chemical variables of the carbon cycle
      !!  --------------------------------------------------------
      DO jk = 1,jpk
         DO jj = 1,jpj
            DO ji = 1,jpi
               caralk = trn(ji,jj,jk,jptal)-       &
                  &        borat(ji,jj,jk)/(1.+1.E-8/(rtrn+akb3(ji,jj,jk)))
               co3(ji,jj,jk)=(caralk-trn(ji,jj,jk,jpdic))*tmask(ji,jj,jk)   &
                  &        +(1.-tmask(ji,jj,jk))*.5e-3
               bicarb = (2.*trn(ji,jj,jk,jpdic)-caralk)
               hi(ji,jj,jk) = (ak23(ji,jj,jk)*bicarb/co3(ji,jj,jk))     &
                  &  *tmask(ji,jj,jk)+(1.-tmask(ji,jj,jk))*1.e-9
               h2co3(ji,jj) = 1.e-5
            ENDDO
         ENDDO
      ENDDO
#endif

   END SUBROUTINE trc_rst

   SUBROUTINE trc_wri(kt)
      !! ==================================================================================
      !!
      !!                       ROUTINE trc_wri
      !!                     ******************
      !!
      !!  PURPOSE :
      !!  ---------
      !!     WRITE restart fields in nutwrs
      !!   METHOD :
      !!   -------
      !!
      !!   nutwrs FILE:
      !!   each nstock time step , SAVE fields which are necessary for
      !!   passive tracer restart
      !!
      !!
      !!   INPUT :
      !!   -----
      !!      argument
      !!              kt              : time step
      !!      COMMON
      !!            /cottrc/          : passive tracers fields (before,now
      !!                                  ,after)
      !!
      !!   OUTPUT :
      !!   ------
      !!      FILE
      !!           nutwrs          : standard restart fields OUTPUT
      !!
      !!   WORKSPACE :
      !!   ---------
      !!      ji,jj,jk,jl,ino0,it0,iarak0
      !!
      !!   History:
      !!   --------
      !!      original : 96-12
      !!      addition : 99-12 (M.-A. Foujols) NetCDF FORMAT with ioipsl
      !!      additions : 00-05 (A. Estublier)
      !!                  TVD Limiter Scheme : key_trc_tvd
      !!      additions : 01-01 (M.A Foujols, E. Kestenare) bug fix: restclo
      !!      additions : 01-01 (O. Aumont, E. Kestenare)
      !!                  write restart file for sediments
      !!      additions : 01-05 (O. Aumont, E. Kestenare)
      !!                  write restart file for calcite and silicate sediments
      !!   05-03 (O. Aumont and A. El Moussaoui) F90
      !!========================================================================================!
      !! * Modules used
      USE ioipsl

      !! * Arguments
      !! -----------
      INTEGER, INTENT( in ) :: kt

      !! * local declarations
      !! ====================

      LOGICAL :: clbon         !!!
      CHARACTER (len=50) :: clname,clname1,clname2,cln

      INTEGER :: jn,   &
         ino0,it0,iarak0,     &
         ic,jc,ji,jj,jk,      &
         itime

      REAL(wp) :: zdate0, zinfo(3),zdiag_var,    &
         zdiag_varmin, zdiag_varmax, zdiag_tot, zder


      !! 1. OUTPUT of restart fields (nutwrs)
      !! ---------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_wri : write passive tracers restart.output NetCDF file'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
       

         areatot = 0.
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  areatot = areatot + tmask(ji,jj,jk)*tmask_i(ji,jj)*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) 
               END DO
            END DO
         END DO

         IF( lk_mpp ) THEN 
             CALL mpp_sum(areatot)     ! sum over the global domain  
         END IF

         trai = 0.
         DO jn = 1, jptra
            DO jk = 1,jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     trai=trai+tmask(ji,jj,jk)*trn(ji,jj,jk,jn)*     &
                        &    tmask_i(ji,jj)* e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) 
                  END DO
               END DO
            END DO
         END DO

         IF( lk_mpp ) THEN 
             CALL mpp_sum(trai)         ! sum over the global domain  
         END IF

         IF (lwp) WRITE(numout,*) 'Integral of all tracers over the full domain at NIT000 =',trai

      ENDIF


      IF( MOD(kt,nstock) == 0 .OR. kt == nitend ) THEN

         !! 0. initialisations
         !! ------------------

         IF(lwp) WRITE(numout,*) ' '
         IF(lwp) WRITE(numout,*) 'trc_wri : write the passive tracer restart file in NetCDF format ',   &
            'at it= ',kt,' date= ',ndastp
         IF(lwp) WRITE(numout,*) '~~~~~~~~~'


         ino0 =no
         it0  =kt
         IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
            iarak0 = 1
         ELSE
            iarak0=0
         ENDIF

         zinfo(1)=FLOAT(ino0)
         zinfo(2)=FLOAT(it0)
         zinfo(3)=FLOAT(iarak0)

         !! 1. WRITE in nutwrs
         !! ------------------
         !!... first information

         INQUIRE (FILE=trestart,EXIST=clbon)
         IF(clbon) THEN
            OPEN(UNIT=nutwrs,FILE=trestart,STATUS='old')
            CLOSE(nutwrs,STATUS='delete')
         ENDIF

         ic=1
         DO jc=1,16
            IF(cexper(jc:jc) /= ' ') ic = jc
         END DO
         WRITE(cln,'("_",i4.4,i2.2,i2.2,"_restart.trc")') nyear, nmonth, nday
         clname=cexper(1:ic)//cln
         ic=1
         DO jc=1,48
            IF(clname(jc:jc) /= ' ') ic=jc
         END DO
         trestart=clname(1:ic)//".nc"
         itime=0
         CALL ymds2ju(nyear,nmonth,nday,0.0,zdate0)
         CALL restini('NONE',jpi,jpj,glamt,gphit,jpk,gdept,clname           &
            &        ,itime,zdate0,rdt*nstock,nutwrs,domain_id=nidom)

         CALL restput(nutwrs,'info',1,1,3,0,zinfo)

         ! prognostic variables
         ! --------------------

         IF (lwp) WRITE(numout,*) '----TRACER STAT----'
         zdiag_tot=0.
         DO jn=1,jptra
            clname='TRN'//ctrcnm(jn)
            CALL restput(nutwrs,clname,jpi,jpj,jpk,0,trn(:,:,:,jn))

            zdiag_var=0.
            zdiag_varmin=0.
            zdiag_varmax=0.

            DO ji=1, jpi
               DO jj=1, jpj
                  DO jk=1,jpk
                    zdiag_var=zdiag_var+tmask(ji,jj,jk)*trn(ji,jj,jk,jn)*     &
                               tmask_i(ji,jj)* e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) 

                  END DO
               END DO
            END DO

            zdiag_varmin=MINVAL(trn(:,:,:,jn), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)))
            zdiag_varmax=MAXVAL(trn(:,:,:,jn), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)))

            IF( lk_mpp ) THEN 
               CALL mpp_min(zdiag_varmin)      ! min over the global domain  
               CALL mpp_max(zdiag_varmax)      ! max over the global domain  
               CALL mpp_sum(zdiag_var)         ! sum over the global domain  
            END IF

            zdiag_tot=zdiag_tot+zdiag_var
            zdiag_var=zdiag_var/areatot

            IF (lwp) WRITE(numout,*) 'MEAN NO ',jn,ctrcnm(jn),' =',zdiag_var,'MIN= '  &
               ,zdiag_varmin,'MAX= ',zdiag_varmax

         END DO

         zdiag_tot=zdiag_tot
         zder=((zdiag_tot-trai)/trai)*100._wp
         IF (lwp) WRITE(numout,*) 'Integral of all tracers over the full domain  =',zdiag_tot  
         IF (lwp) WRITE(numout,*) 'Drift of the sum of all tracers =',zder, '%'  

         DO jn=1,jptra
            clname='TRB'//ctrcnm(jn)
            CALL restput(nutwrs,clname,jpi,jpj,jpk,0,trb(:,:,:,jn))
         END DO

#if defined key_trc_lobster1
         clname='SEDB'//ctrcnm(jpdet)
         clname1='SEDN'//ctrcnm(jpdet)
         CALL restput(nutwrs,clname,jpi,jpj,1,0,sedpocb(:,:))
         CALL restput(nutwrs,clname1,jpi,jpj,1,0,sedpocn(:,:))
#elif defined key_trc_pisces
         clname='SED'//ctrcnm(jppoc)
         clname1='SED'//ctrcnm(jpcal)
         clname2='SED'//ctrcnm(jpsil)
         CALL restput(nutwrs,clname1,jpi,jpj,1,0,sedcal(:,:))
         CALL restput(nutwrs,clname2,jpi,jpj,1,0,sedsil(:,:))
         CALL restput(nutwrs,clname,jpi,jpj,1,0,sedpoc(:,:))

         clname='Silicalim'
         CALL restput(nutwrs,clname,jpi,jpj,1,0,xksi(:,:))
#elif defined key_cfc
         clname='qint'
         CALL restput(nutwrs,clname,jpi,jpj,jptra,0,qint(:,:,:))
         clname1='qtr'
         CALL restput(nutwrs,clname1,jpi,jpj,jptra,0,qtr(:,:,:))
#endif


         CALL restclo(nutwrs)

      ENDIF

   END SUBROUTINE trc_wri

#endif

#else
   !!======================================================================
   !!  Empty module : No passive tracer
   !!======================================================================
CONTAINS

   SUBROUTINE trc_rst
      !! no passive tracers
   END SUBROUTINE trc_rst

   SUBROUTINE trc_wri(kt)
      !! no passive tracers
      INTEGER, INTENT ( in ) :: kt
      WRITE(*,*) 'trc_wri: You should not have seen this print! error?', kt
   END SUBROUTINE trc_wri
   
#endif
   
END MODULE trcrst
