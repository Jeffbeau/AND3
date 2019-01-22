MODULE trclec
   !!==========================================================================
   !!
   !!                       *** MODULE trclec ***
   !! Read and print options for the passive tracer run (namelist)
   !! O.Aumont and A.El Moussaoui 03/05 F90
   !!=========================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trclec.F90,v 1.5 2006/04/10 15:40:28 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_passivetrc
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   USE trctrp_lec
   USE trclsm

   IMPLICIT NONE
   PRIVATE 

   !! * Accessibility
   PUBLIC trc_lec

#include "passivetrc_substitute.h90"

CONTAINS

   SUBROUTINE trc_lec
      !!---------------------------------------------------------------------
      !!                       ROUTINE trclec
      !!                     ******************
      !!  PURPOSE :
      !!  ---------
      !!     READ and PRINT options for the passive tracer run (namelist)
      !!
      !!   History:
      !!   --------
      !!      original  : 96-11 (M.A. Foujols, M. Levy) passive tracer
      !!      modification : 98-04 (M.A Foujols, L. Bopp) ahtrb0 for isopycnal
      !!                                                  diffusion
      !!      modification : 99-10(M.A. Foujols, M. Levy) separation of sms
      !!      additions : 00-05(A. Estublier) TVD Limiter Scheme : Tests 
      !!                                      on ndttrc
      !!      additions : 00-06(A. Estublier) MUSCL Scheme : Tests 
      !!                                      on ndttrc
      !!      additions : 00-07(A. Estublier) PPM Scheme : Tests on ndttrc
      !!      modification : 00-11 (M.A Foujols, E Kestenare) trcrat, ahtrc0 and aeivtr0
      !!      modification : 01-01 (E Kestenare) suppress ndttrc=1 
      !!                                         for Arakawa and TVD schemes
      !!     O.Aumont and A.El Moussaoui 03/05 F90
      !!----------------------------------------------------------------------

      !! local declarations
      !! ==================

      INTEGER ::  ji
      CHARACTER (len=32) :: clname

      !!---------------------------------------------------------------------
      !!  OPA.90   03/2005 
      !!---------------------------------------------------------------------

      !! 0. initializations
      !! ------------------

      namelist/nattrc/nwritetrc,lrsttr,nrsttr, ctrcnm,ctrcnl,ctrcun,lutini     !general   

      namelist/natnum/rsc,rtrn,ncortrc,ndttrc,crosster

#if defined key_trc_diatrd
      namelist/natrtd/luttrd,nwritetrd                      ! dynamical trends
#endif

#if defined key_trc_diaadd
      namelist/natadd/ctrc3d,ctrc3l,ctrc2d,ctrc2l, ctrc3u, ctrc2u,     &
         nwriteadd                             !additional diagnostics
#endif

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' ROUTINE trclec'
         WRITE(numout,*) ' **************'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist for passive tracer'
         WRITE(numout,*) ' ***************************'
         WRITE(numout,*) ' '
      ENDIF

      numnat=80
      REWIND (numnat)

      clname='namelist.passivetrc'
      OPEN( numnat, FILE= clname, FORM='formatted', STATUS = 'old')

      !! 1., 2. & 3. initialization with namelist files
      !! ----------------------------------------------
      !! 1.0 namelist nattrc :

      nwritetrc = 10
      lrsttr=.FALSE.
      nrsttr = 0

      DO ji=1,jptra
         WRITE (ctrcnm(ji),'("TR_",I1)') ji
         WRITE (ctrcnl(ji),'("TRACER NUMBER ",I1)') ji
         ctrcun(ji)='mmole/m3'
         lutini(ji)=.FALSE. 
      END DO


      REWIND(numnat)
      READ(numnat,nattrc)

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) 'nattrc'
         WRITE(numout,*) ' '
         WRITE(numout,*)          &
            ' frequency of outputs for passive tracers nwritetrc = '    &
            ,nwritetrc  
         WRITE(numout,*) ' restart LOGICAL for passive tr. lrsttr = ',   &
            &         lrsttr
         WRITE(numout,*) ' control of time step for p. tr. nrsttr = ',   & 
            &         nrsttr
         DO ji=1,jptra
            WRITE(numout,*) ' tracer nb: ',ji,' name = ',ctrcnm(ji)       & 
               &           ,ctrcnl(ji) 
            WRITE(numout,*) ' in unit = ',ctrcun(ji)
            WRITE(numout,*) ' initial value in FILE : ',lutini(ji) 
            WRITE(numout,*) ' '
         END DO
         WRITE(numout,*) ' '
      ENDIF

#if defined key_trc_diatrd

      !! 1.2 namelist nattrd : passive tracers dynamical trends

      nwritetrd=10

      !! default : no dynamical trend recording
      !! --------------------------------------
      DO ji=1,jptra
         luttrd(ji) = .FALSE.
      END DO

      REWIND(numnat)
      READ(numnat,natrtd)

      nkeep=0
      ikeep(:)=0
      DO ji=1,jptra
         IF (luttrd(ji)) THEN 
             nkeep=nkeep+1
             ikeep(ji)=nkeep
         END IF 
      END DO
      IF (nkeep.GT.0) THEN  
        IF (.NOT. ALLOCATED(trtrd)) ALLOCATE(trtrd(jpi,jpj,jpk,nkeep,jpdiatrc)) 
        trtrd(:,:,:,:,:)=0.0
      ENDIF 
      IF(lwp) THEN
         WRITE(numout,*) 'natrtd'
         WRITE(numout,*) ' '
         WRITE(numout,*)                        &
            ' frequency of outputs for dynamical trends nwritetrd = '   &
            ,nwritetrd
         DO ji=1,jptra
            WRITE(numout,*)                      &
               ' keep dynamical trends for tracer number :',ji          &
               ,luttrd(ji), ikeep(ji)
         END DO
         WRITE(numout,*) 'total = ',nkeep,' tracers dyn trends saved'
         WRITE(numout,*) 'size of trtrd = ',jpi*jpj*jpk*nkeep*jpdiatrc
      ENDIF
#endif

      !!1.3 namelist natadd : passive tracers diagnostics
      !!-------------------------------------------------

#if defined key_trc_diaadd

      nwriteadd = 10

      !! default value for 3D output arrays : short and long name, units

      DO ji=1,jpdia3d
         WRITE (ctrc3d(ji),'("3D_",I1)') ji
         WRITE (ctrc3l(ji),'("3D DIAGNOSTIC NUMBER ",I2)') ji
         ctrc3u(ji)=' '
      END DO


      !! default value for 2D output arrays : short and long name, units
      !! ---------------------------------------------------------------
      DO ji=1,jpdia2d
         WRITE (ctrc2d(ji),'("2D_",I1)') ji
         WRITE (ctrc2l(ji),'("2D DIAGNOSTIC NUMBER ",I2)') ji
         ctrc2u(ji)=' '
      END DO

      REWIND(numnat)
      READ(numnat,natadd)

      IF(lwp) THEN
         WRITE(numout,*) ' natadd'
         WRITE(numout,*) ' '
         WRITE(numout,*)                          &
            ' frequency of outputs for additional arrays nwriteadd = '   &
            ,nwriteadd
         DO ji=1,jpdia3d
            WRITE(numout,*)                     &
               'name of 3d output field number :',ji,' : ',ctrc3d(ji)  
            WRITE(numout,*) ctrc3l(ji)  
            WRITE(numout,*) ' in unit = ',ctrc3u(ji)
         END DO
         WRITE(numout,*) ' '
         DO ji=1,jpdia2d
            WRITE(numout,*)                    &
               'name of 2d output field number :',ji,' : ',ctrc2d(ji)  
            WRITE(numout,*) ctrc2l(ji)  
            WRITE(numout,*) ' in unit = ',ctrc2u(ji)
         END DO
         WRITE(numout,*) ' '
      ENDIF
#endif

      !! 1.1 namelist natnum :
      !! ---------------------
      rsc=1.
      rtrn=1.e-15
      ncortrc=1
      ndttrc=4
      crosster=.FALSE.

      REWIND(numnat)
      READ(numnat,natnum)

!!Chris  computes the first time step of tracer model
      nittrc000 = nit000 + ndttrc - 1

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) 'natnum'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' tuning coefficient              rsc     = ',    &
            rsc
         WRITE(numout,*) ' truncation value                rtrn    = ',    &
            rtrn
         WRITE(numout,*) ' number of corrective phase      ncortrc = ',    &
            ncortrc
         WRITE(numout,*) ' time step freq. for pass. trac. ndttrc  = ',    &
            ndttrc
         WRITE(numout,*) ' 1st time step for pass. trac. nittrc000 = ',    &
            nittrc000
         WRITE(numout,*) ' computes or not crossterms    crosster  = ',    &
            crosster
      ENDIF


      !! namelist of transport
      !! ---------------------
      CALL trc_trp_lec

      !! namelist of SMS
      !! ---------------      
      CALL trc_lsm

   END SUBROUTINE trc_lec

#else
   !!======================================================================
   !!  Empty module : No passive tracer
   !!======================================================================
CONTAINS

   SUBROUTINE trc_lec

   END SUBROUTINE trc_lec

#endif

END MODULE  trclec
