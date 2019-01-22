MODULE flowri
   !!======================================================================
   !!                       ***  MODULE  flowri  ***
   !! 
   !!======================================================================
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                     float trajectories
   !!----------------------------------------------------------------------
   !!    flowri     : write trajectories of floats in file 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE flo_oce         ! ocean drifting floats
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! distribued memory computing library
   USE daymod
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE

   !! * Accessibility
   PRIVATE
   PUBLIC flo_wri     ! routine called by floats.F90

   !! * Module variables
      INTEGER :: jfl              ! number of floats

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/flowri.F90,v 1.3 2005/03/27 18:35:05 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE flo_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_wri  ***
      !!             
      !! ** Purpose :   Write position of floats in "trajec_float" file
      !!      and the temperature and salinity at this position
      !!      
      !! ** Method  :   The frequency is nwritefl
      !!      
      !!  History :
      !!    8.0  !  99-09  (Y. Drillet)  Original code
      !!         !  00-06  (J.-M. Molines)  Profiling floats for CLS 
      !!    8.5  !  02-10  (A. Bozec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER  :: kt                               ! time step

      !! * Local declarations
      CHARACTER (len=21) ::  clname
      INTEGER ::   inum = 11       ! temporary logical unit for restart file
      INTEGER  ::   &
         iafl,ibfl,icfl,ia1fl,ib1fl,ic1fl,jfl,irecflo,   &
         iafloc,ibfloc,ia1floc,ib1floc,   &
         iafln, ibfln
      INTEGER  ::    ic, jc , jpn
      INTEGER, DIMENSION ( jpnij )  :: iproc

      REAL(wp) :: zafl,zbfl,zcfl,zdtj
      REAL(wp) :: zxxu, zxxu_01,zxxu_10, zxxu_11
      REAL(wp), DIMENSION (jpk,jpnfl) :: ztemp, zsal
      !!---------------------------------------------------------------------
      
      IF( kt == nit000 .OR. MOD( kt,nwritefl)== 0 ) THEN 

         ! header of output floats file
      
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'flo_wri : write in trajec_float file '
            WRITE(numout,*) '~~~~~~~    '
         ENDIF

         ! open the file numflo 
         clname='trajec_float'
         OPEN (numflo,FILE=clname,FORM='UNFORMATTED')
         !      REWIND numflo

         IF( kt == nit000 ) THEN
            irecflo = NINT( (nitend-nit000) / FLOAT(nwritefl) )
            IF(lwp) WRITE(numflo)cexper,no,irecflo,jpnfl,nwritefl
         ENDIF
         zdtj = rdt / 86400.      !!bug   use of 86400 instead of the phycst parameter

         ! translation of index position in geographical position

         IF( lk_mpp ) THEN
            DO jfl = 1, jpnfl
               iafl  = INT ( tpifl(jfl) )
               ibfl  = INT ( tpjfl(jfl) )
               icfl  = INT ( tpkfl(jfl) )
               iafln = NINT( tpifl(jfl) )
               ibfln = NINT( tpjfl(jfl) )
               ia1fl = iafl + 1
               ib1fl = ibfl + 1
               ic1fl = icfl + 1
               zafl  = tpifl(jfl) - FLOAT( iafl )
               zbfl  = tpjfl(jfl) - FLOAT( ibfl )
               zcfl  = tpkfl(jfl) - FLOAT( icfl )
               IF(   iafl >= mig(nldi)-jpizoom+1 .AND. iafl <= mig(nlei)-jpizoom+1 .AND.   &
                  &  ibfl >= mjg(nldj)-jpjzoom+1 .AND. ibfl <= mjg(nlej)-jpjzoom+1       ) THEN

                  ! local index

                  iafloc  = iafl -(mig(1)-jpizoom+1) + 1
                  ibfloc  = ibfl -(mjg(1)-jpjzoom+1) + 1
                  ia1floc = iafloc + 1
                  ib1floc = ibfloc + 1

                  flyy(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
                     &      +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)
                  flxx(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
                     &      +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
                  flzz(jfl) = (1.-zcfl)*fsdepw(iafloc,ibfloc,icfl ) + zcfl * fsdepw(iafloc,ibfloc,ic1fl)

                  ! Change  by Alexandra Bozec et Jean-Philippe Boulanger
                  ! We save  the instantaneous profile of T and S of the column     
                  ! ztemp(jfl)=tn(iafloc,ibfloc,icfl)
                  ! zsal(jfl)=sn(iafloc,ibfloc,icfl)
                  ztemp(1:jpk,jfl) = tn(iafloc,ibfloc,1:jpk)
                  zsal (1:jpk,jfl) = sn(iafloc,ibfloc,1:jpk)            
               ELSE
                  flxx(jfl) = 0.
                  flyy(jfl) = 0.
                  flzz(jfl) = 0.
                  ztemp(1:jpk,jfl) = 0.
                  zsal (1:jpk,jfl) = 0.
               ENDIF
            END DO

            CALL mpp_sum( flxx, jpnfl )   ! sums over the global domain
            CALL mpp_sum( flyy, jpnfl )
            CALL mpp_sum( flzz, jpnfl )
            ! these 2 lines have accendentaly been removed from ATL6-V8 run hence
            ! giving 0 salinity and temperature on the float trajectory
            CALL mpp_sum( ztemp, jpk*jpnfl )
            CALL mpp_sum( zsal , jpk*jpnfl )

         ELSE
            DO jfl = 1, jpnfl
               iafl  = INT (tpifl(jfl))
               ibfl  = INT (tpjfl(jfl))
               icfl  = INT (tpkfl(jfl))
               iafln = NINT(tpifl(jfl))
               ibfln = NINT(tpjfl(jfl))
               ia1fl = iafl+1
               ib1fl = ibfl+1
               ic1fl = icfl+1
               zafl  = tpifl(jfl) - FLOAT(iafl)
               zbfl  = tpjfl(jfl) - FLOAT(ibfl)
               zcfl  = tpkfl(jfl) - FLOAT(icfl)
               iafloc  = iafl
               ibfloc  = ibfl
               ia1floc = iafloc + 1
               ib1floc = ibfloc + 1
               !
               flyy(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
                         +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)
               flxx(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
                         +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
               flzz(jfl) = (1.-zcfl)*fsdepw(iafloc,ibfloc,icfl ) + zcfl * fsdepw(iafloc,ibfloc,ic1fl)
               !ALEX
               ! Astuce pour ne pas avoir des flotteurs qui se baladent sur IDL
               zxxu_11 = glamt(iafloc ,ibfloc )
               zxxu_10 = glamt(iafloc ,ib1floc)
               zxxu_01 = glamt(ia1floc,ibfloc )
               zxxu    = glamt(ia1floc,ib1floc)

               IF( iafloc == 52 )  zxxu_10 = -181
               IF( iafloc == 52 )  zxxu_11 = -181
               flxx(jfl)=(1.-zafl)*(1.-zbfl)* zxxu_11 + (1.-zafl)*    zbfl * zxxu_10   &
                        +    zafl *(1.-zbfl)* zxxu_01 +     zafl *    zbfl * zxxu
               !ALEX         
               ! Change  by Alexandra Bozec et Jean-Philippe Boulanger
               ! We save  the instantaneous profile of T and S of the column     
               !     ztemp(jfl)=tn(iafloc,ibfloc,icfl)
               !     zsal(jfl)=sn(iafloc,ibfloc,icfl)
               ztemp(1:jpk,jfl) = tn(iafloc,ibfloc,1:jpk)
               zsal (1:jpk,jfl) = sn(iafloc,ibfloc,1:jpk)
            END DO
         ENDIF

         !
         WRITE(numflo) flxx,flyy,flzz,nisobfl,ngrpfl,ztemp,zsal, FLOAT(ndastp)
      !!
      !! case when profiles are dumped. In order to save memory, dumps are
      !! done level by level.
      !      IF (mod(kt,nflclean) == 0.) THEN
      !!     IF ( nwflo == nwprofil ) THEN
      !        DO jk = 1,jpk
      !         DO jfl=1,jpnfl
      !         iafl= INT(tpifl(jfl))
      !         ibfl=INT(tpjfl(jfl))
      !         iafln=NINT(tpifl(jfl))
      !         ibfln=NINT(tpjfl(jfl))
      !# if defined key_mpp_mpi   ||   defined key_mpp_shmem
      !        IF ( (iafl >= (mig(nldi)-jpizoom+1)) .AND.
      !     $       (iafl <= (mig(nlei)-jpizoom+1)) .AND.
      !     $       (ibfl >= (mjg(nldj)-jpjzoom+1)) .AND.
      !     $       (ibfl <= (mjg(nlej)-jpjzoom+1)) ) THEN
      !!
      !! local index
      !!
      !         iafloc=iafln-(mig(1)-jpizoom+1)+1
      !         ibfloc=ibfln-(mjg(1)-jpjzoom+1)+1
      !!         IF (jk == 1 ) THEN
      !!      PRINT *,'<<<>>> ',jfl,narea, iafloc ,ibfloc, iafln, ibfln,adatrj
      !!         ENDIF
      !# else
      !         iafloc=iafln
      !         ibfloc=ibfln
      !# endif
      !         ztemp(jfl)=tn(iafloc,ibfloc,jk)
      !         zsal(jfl)=sn(iaflo!,ibfloc,jk)
      !# if defined key_mpp_mpi   ||   defined key_mpp_shmem
      !        ELSE
      !         ztemp(jfl) = 0.
      !         zsal(jfl) = 0.
      !        ENDIF
      !# endif
      !! ... next float
      !        END DO
      !      IF( lk_mpp )   CALL mpp_sum( ztemp, jpnfl )
      !      IF( lk_mpp )   CALL mpp_sum( zsal , jpnfl )
      !
      !      IF (lwp) THEN 
      !         WRITE(numflo) ztemp, zsal
      !      ENDIF
      !! ... next level jk
      !      END DO
      !! ... reset nwflo to 0 for ALL processors, if profile has been written
      !!       nwflo = 0
      !      ENDIF
      !!
      !      CALL flush (numflo)
      !! ... time of dumping floats
      !!      END IF
      ENDIF
      
      IF( (MOD(kt,nstockfl) == 0) .OR. ( kt == nitend ) ) THEN 
         ! Writing the restart file 
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'flo_wri : write in  restart_float file '
            WRITE(numout,*) '~~~~~~~    '
         ENDIF

         ! file is opened and closed every time it is used.

         clname = 'restart.float.'
         ic = 1
         DO jc = 1, 16
            IF( cexper(jc:jc) /= ' ' ) ic = jc
         END DO
         clname = clname(1:14)//cexper(1:ic)
         ic = 1
         DO jc = 1, 48
            IF( clname(jc:jc) /= ' ' ) ic = jc
         END DO

         OPEN (inum,FILE=clname,FORM='UNFORMATTED')
         REWIND inum
         !
         DO jpn = 1, jpnij
            iproc(jpn) = 0
         END DO
         !
         IF(lwp) THEN
            REWIND(inum)
            WRITE (inum) tpifl,tpjfl,tpkfl,nisobfl,ngrpfl
            CLOSE (inum) 
         ENDIF
         !
         ! Compute the number of trajectories for each processor
         !
         IF( lk_mpp ) THEN
            DO jfl = 1, jpnfl
               IF( (INT(tpifl(jfl)) >= (mig(nldi)-jpizoom+1)) .AND.   &
                  &(INT(tpifl(jfl)) <= (mig(nlei)-jpizoom+1)) .AND.   &
                  &(INT(tpjfl(jfl)) >= (mjg(nldj)-jpjzoom+1)) .AND.   &
                  &(INT(tpjfl(jfl)) <= (mjg(nlej)-jpjzoom+1)) ) THEN
                  iproc(narea) = iproc(narea)+1
               ENDIF
            END DO
            CALL mpp_sum( iproc, jpnij )
            !
            IF(lwp) THEN 
               WRITE(numout,*) 'DATE',adatrj
               DO jpn = 1, jpnij
                  IF( iproc(jpn) /= 0 ) THEN
                     WRITE(numout,*)'PROCESSOR',jpn-1,'compute',iproc(jpn), 'trajectories.'
                  ENDIF
               END DO
            ENDIF
         ENDIF
      ENDIF 

   END SUBROUTINE flo_wri

#  else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flo_wri                 ! Empty routine
   END SUBROUTINE flo_wri
#endif
   
   !!======================================================================
END MODULE flowri
