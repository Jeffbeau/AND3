CC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/trcopt.F,v 1.5 2005/11/14 12:30:28 opalod Exp $
CDIR$ LIST
      SUBROUTINE trcopt(kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcopt
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     computes the light propagation in the water column
CCC     and the euphotic layer depth
CCC
CCC
CC   METHOD :
CC   -------
CC
CC      multitasked on vertical slab (jj-loop)
CC      local par is computed in w layers using light propagation
CC      mean par in t layers are computed by integration
CC
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC                                depths
CC            /comzdf/          : avt vertical eddy diffusivity
CC            /comqsr/          : solar radiation
CC            /comtsk/          : multitasking
CC            /cotopt/          : optical parameters
CC            /cotbio/          : biological parameters
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cotopt/          : optical parameters
CC
CC   WORKSPACE :
CC   ---------
CC      local     zparr         : red compound of par
CC                zparg         : green compound of par
CC                zpar0m        : irradiance just below the surface
CC                zpar100       : irradiance at euphotic layer depth
CC                zkr           : total absorption coefficient in red
CC                zkg           : total absorption coefficient in green
CC                zpig          : total pigment
CC                imaske        : euphotic layer mask
CC                itabe         : euphotic layer last k index
CC
CC      COMMON
CC
CC   EXTERNAL :                   no
CC   --------
CC
CC   REFERENCES :                 no
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC       original : 95-05 (M. Levy)
CC                  99-09 (J-M Andre & M. Levy)
CC       modifications : 99-11 (C. Menkes M.A. Foujols) itabe initial. 
CC       modifications : 00-02 (M.A. Foujols) change x**y par exp(y*log(x))
CC----------------------------------------------------------------------
CDIR$ NOLIST

      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
CDIR$ LIST
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (11/96)
CCC---------------------------------------------------------------------
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kt

#if defined key_passivetrc && defined key_trc_lobster1
C
      INTEGER ji,jj,jk,jn,in

      REAL zpig,zkr,zkg

      REAL zparr(jpi,jpk),zparg(jpi,jpk)
      REAL zpar0m(jpi),zpar100(jpi)
      INTEGER itabe(jpi),imaske(jpi,jpk)
CC----------------------------------------------------------------------
CC statement functions
CC ===================
CDIR$ NOLIST
#include "domzgr_substitute.h90"
CDIR$ LIST
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
C
C find Phytoplancton index - test CTRCNM
C
      in=0
      DO jn = 1,jptra
        IF ((ctrcnm(jn) .EQ. 'PHY') .OR.
     $      (ctrcnm(jn) .EQ. 'PHYTO') ) THEN
            
            in = jn
        END IF
      END DO
      IF (in.eq.0) THEN 
          IF (lwp) THEN 
              WRITE (numout,*)
     $            ' Problem trcopt : PHY or PHYTO not found '
              CALL FLUSH(numout)
          ENDIF 
      ENDIF 
C
C vertical slab
C ===============
C
      DO 1000 jj = 1,jpj
C
C
C 1. determination of surface irradiance
C --------------------------------------
C
C
        DO ji = 1,jpi
          zpar0m(ji) = qsr(ji,jj)*0.43
          zpar100(ji) = zpar0m(ji)*0.01
          xpar(ji,jj,1) = zpar0m(ji)
          zparr(ji,1) = 0.5* zpar0m(ji)
          zparg(ji,1) = 0.5* zpar0m(ji)
        END DO

C
C 2. determination of xpar
C ------------------------
C
C determination of local par in w levels
        DO jk = 2,jpk
          DO ji = 1,jpi
            zpig = max(tiny(0.),trn(ji,jj,jk - 1,in))*12*redf/rcchl/rpig
            zkr = xkr0 + xkrp*exp(xlr*log(zpig))
            zkg = xkg0 + xkgp*exp(xlg*log(zpig))
            zparr(ji,jk) = zparr(ji,jk - 1)
     $          *exp( -zkr*fse3t(ji,jj,jk - 1) )
            zparg(ji,jk) = zparg(ji,jk - 1)
     $          *exp( -zkg*fse3t(ji,jj,jk - 1) )
          END DO
        END DO

C
C mean par in t levels
        DO jk = 1,jpkm1
          DO ji = 1,jpi
            zpig = max(tiny(0.),trn(ji,jj,jk  ,in))*12*redf/rcchl/rpig
            zkr = xkr0 + xkrp*exp(xlr*log(zpig))
            zkg = xkg0 + xkgp*exp(xlg*log(zpig))
            zparr(ji,jk) = zparr(ji,jk) / zkr / fse3t(ji,jj,jk)
     $          * ( 1 - exp( -zkr*fse3t(ji,jj,jk) ) )
            zparg(ji,jk) = zparg(ji,jk) / zkg / fse3t(ji,jj,jk)
     $          * ( 1 - exp( -zkg*fse3t(ji,jj,jk) ) )
            xpar(ji,jj,jk) = max(zparr(ji,jk)
     $          + zparg(ji,jk),1.e-15)
          END DO
        END DO
C
C
C 4. determination of euphotic layer depth 
C ----------------------------------------
C
C imaske equal 1 in the euphotic layer, and 0 without
C
        DO jk = 1,jpk
          DO ji = 1,jpi
            imaske(ji,jk) = 0
            IF (xpar(ji,jj,jk) .GE. zpar100(ji)) imaske(ji,jk) = 1
          END DO
        END DO
C
        DO ji = 1,jpi
          itabe(ji) = 0
        END DO
C
        DO jk = 1,jpk
          DO ji = 1,jpi
            itabe(ji) = itabe(ji) + imaske(ji,jk)
          END DO
        END DO
C
        DO ji = 1,jpi
          itabe(ji) = max(1,itabe(ji))
          xze(ji,jj) = fsdepw(ji,jj,itabe(ji) + 1)
        END DO 
C
C
C END of slab
C ===========
C
 1000 CONTINUE
C
#else
C
C no passive tracers
C
#endif
C
      RETURN
      END
