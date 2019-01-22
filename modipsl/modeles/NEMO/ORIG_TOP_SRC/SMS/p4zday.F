CCC$Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zday.F,v 1.4 2005/09/12 09:04:57 opalod Exp $
CCC  TOP 1.0 , LOCEAN-IPSL (2005)
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
C ---------------------------------------------------------------------------
CDIR$ LIST
      SUBROUTINE p4zday
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC          ROUTINE p4zday : PISCES MODEL
CCC          *****************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC        PISCES : compute the day length depending on latitude
CCC                 and the day
CCC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC      common
CC              all the common defined in opa
CC
CC
CC   OUTPUT :                   : no
CC   ------
CC
CC   EXTERNAL :
CC   --------
CC            None
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : E. Maier-Reimer (GBC 1993)
CC      additions : C. Le Quere (1999)
CC      modifications : O. Aumont (2004)
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj, iyy
      REAL rum, delta, codel, phi, argu
C
C Get year
C --------
C
      iyy = ndastp/10000

       IF(lwp) write(numout,*)
       IF(lwp) write(numout,*) 'p4zday - Julian day ', nday_year
       IF(lwp) write(numout,*)


      IF (nleapy.EQ.1 .AND. MOD(iyy,4).EQ.0) THEN
          rum = FLOAT(nday_year-80)/366.
      ELSE
          rum = FLOAT(nday_year-80)/365.
      ENDIF

C
      delta = SIN(rum*rpi*2.)*sin(rpi*23.5/180.)
      codel = asin(delta)

      DO jj = 1,jpj
        DO ji = 1,jpi
          phi = gphit(ji,jj)*rpi/180.
          argu = tan(codel)*tan(phi)
          strn(ji,jj) = 0.
          argu=min(1.,argu)
          argu=max(-1.,argu)
          strn(ji,jj)=24.-2.*acos(argu)*180./rpi/15.
          strn(ji,jj)=max(strn(ji,jj),0.)
        END DO
      END DO
C
#endif
      RETURN
      END
