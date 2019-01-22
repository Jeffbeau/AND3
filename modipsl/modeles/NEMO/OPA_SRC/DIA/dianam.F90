MODULE dianam
   !!======================================================================
   !!                       ***  MODULE  dianam  ***
   !! Ocean diagnostics:  Builds output file name
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   dia_nam       : Builds output file name
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE daymod          ! calendar

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dia_nam   ! routine called by step.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/dianam.F90,v 1.3 2005/03/27 18:34:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_nam( cdfnam, kfreq, cdsuff )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_nam  ***
      !!                   
      !! ** Purpose :   Builds output file name
      !!
      !! ** Method  :   File name is a function of date and output frequency
      !!      cdfnam=<cexper>_<clave>_<idtbeg>_<idtend>_grid_<cdsuff>
      !!      <clave> = averaging frequency (DA, MO, etc...)
      !!      <idtbeg>,<idtend> date of beginning and end of run
      !!
      !! History :
      !!        !  99-02  (E. Guilyardi)  Creation for 30 days/month
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER (len=*), INTENT( out ) ::   cdfnam   ! file name
      CHARACTER (len=*), INTENT( in  ) ::   cdsuff   ! ???
      INTEGER,           INTENT( in  ) ::   kfreq    ! ???

      !! * Local declarations
      CHARACTER (len=8) ::   clexper
      CHARACTER (len=2) ::   clave
      CHARACTER (len=5) ::   clout
      CHARACTER (len=6) ::   clsuff
      INTEGER :: jt, jc, jd, je            ! dummy loop indices
      INTEGER ::   &
         ic, id, ie, ig, ijjmm, iout,   &  ! temporary integers
         iyear1, imonth1, iday1,        &  !    "          "
         iyear2, imonth2, iday2            !    "          "
      REAL(wp) ::   &
         z5j, znbsec, zdate1, zdate2, zdrun, zdt   ! temporary scalars
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' dia_nam: building output file name'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~'
      IF(lwp) WRITE(numout,*)

      ! 0. Initialisation
      ! -----------------

      cdfnam = ''

      !    number of seconds of the run

      z5j = 5*rjjss
      zdt = rdt
      IF( nacc == 1 ) zdt = rdtmin
      zdrun = FLOAT( nitend - nit000 ) * zdt

      !  date of beginning of run

      iyear1  = ndastp/10000
      imonth1 = ndastp/100 - iyear1*100
      iday1   = ndastp - imonth1*100 - iyear1*10000
      IF( nleapy == 1) THEN 
         ijjmm=0
         IF( MOD( iyear1, 4 ) == 0 ) THEN
            DO jt = 1, imonth1-1
               ijjmm = ijjmm + nbiss(jt)
            END DO
         ELSE
            DO jt = 1, imonth1-1
               ijjmm = ijjmm + nobis(jt)
            END DO
         ENDIF
         ijjmm = ijjmm + (iyear1-1)/4
         zdate1 = ( (iyear1-1)*365 + ijjmm +iday1-1 ) * rjjss   
      ELSE IF( nleapy == 0 ) THEN
         ijjmm = 0
         DO jt = 1, imonth1-1
            ijjmm = ijjmm + nobis(jt)
         END DO
         zdate1 = ( (iyear1-1)*raajj + ijjmm + iday1-1)* rjjss
      ELSE 
         zdate1 = ( (iyear1-1)*nleapy*raamo + (imonth1-1)*nleapy + iday1-1)* rjjss
      ENDIF

      !  date of end of run (= date of beginning of next run)

      zdate2 = zdate1 + zdrun
      IF( nleapy == 1 ) THEN 
         iyear2 = zdate2/(365.25*rjjss)+1
         ijjmm = INT(zdate2/rjjss)-365*(iyear2-1)-(iyear2-1)/4
         IF( ijjmm < 0 ) THEN
            iyear2 = iyear2-1
            ijjmm = zdate2/rjjss-365.*(iyear2-1)-(iyear2-1)/4
         ENDIF
         IF( MOD( iyear2, 4 ) == 0 ) THEN
            DO jt = 1, 12
               ijjmm = ijjmm - nbiss(jt)
               IF( ijjmm <= 0 ) go to 10
            END DO
            jt = 12
10          CONTINUE
            imonth2 = jt
            ijjmm = 0
            DO jt = 1, jt-1
               ijjmm = ijjmm + nbiss(jt)
            END DO
         ELSE
            DO jt = 1, 12
               ijjmm = ijjmm - nobis(jt)
               IF( ijjmm <= 0 ) go to 15
            END DO
            jt = 12
15          CONTINUE
            imonth2 = jt
            ijjmm = 0
            DO jt = 1, jt-1
               ijjmm = ijjmm + nobis(jt)
            END DO
         ENDIF
         iday2 = zdate2/rjjss-365.*(iyear2-1)-ijjmm+1-(iyear2-1)/4     
      ELSE IF( nleapy == 0 ) THEN
         iyear2 = zdate2/raass+1
         ijjmm  = zdate2/rjjss-raajj*(iyear2-1)
         DO jt = 1, 12
            ijjmm = ijjmm - nobis(jt)
            IF(ijjmm <= 0) go to 20
         END DO
         jt = 12
20       CONTINUE
         imonth2 = jt
         ijjmm = 0
         DO jt = 1, jt-1
            ijjmm = ijjmm + nobis(jt)
         END DO
         iday2 = zdate2/rjjss-raajj*(iyear2-1)-ijjmm+1          
      ELSE 
         zdate2 = zdate2 / rjjss
         imonth2 = zdate2/FLOAT(nleapy)
         iday2 = zdate2 - imonth2*FLOAT(nleapy) + 1.
         iyear2 = imonth2/12
         imonth2 = imonth2 - iyear2*12
         imonth2 = imonth2 + 1
         iyear2 = iyear2 + 1
         IF( iday2 == 0 ) THEN
            iday2 = nleapy
            imonth2 = imonth2 - 1
            IF( imonth2 == 0 ) THEN
               imonth2 = 12
               iyear2 = iyear2 - 1
            ENDIF
         ENDIF
      ENDIF


      ! 1. Define time averaging period <nn><type>
      !    ---------------------------------------

      iout = 0
#if defined key_diainstant
      clave = 'IN'
      IF( iyear2 <= 99 ) THEN 
         WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
      ELSE IF( iyear2 <= 999 ) THEN 
         WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
      ELSE IF( iyear2 <= 9999 ) THEN 
         WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
      ELSE
         WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
      ENDIF
#else

      znbsec=kfreq*zdt
      ! daily output
      IF( znbsec == rjjss ) THEN
         clave = '1d'
         IF( iyear2 <= 99 ) THEN 
            WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 999 ) THEN 
            WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 9999 ) THEN 
            WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE
            WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ENDIF
         ! 5 day output 
      ELSE IF( znbsec == z5j ) THEN
         clave='5d'
         IF( iyear2 <= 99 ) THEN 
            WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 999 ) THEN 
            WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 9999 ) THEN 
            WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE
            WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ENDIF
         ! monthly ouput 
      ELSE IF( (znbsec == rmoss .AND. nleapy > 1) .OR.   &
               (znbsec >= 28*rjjss .AND. znbsec <= 31*rjjss .AND. nleapy <= 1) ) THEN
         clave = '1m'
         IF( iyear2 <= 99 ) THEN 
            WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 999 ) THEN 
            WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 9999 ) THEN 
            WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE
            WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ENDIF
         ! annual output
      ELSE IF( (znbsec == raass .AND. nleapy > 1) .OR.   &
               (znbsec >= 365*rjjss .AND. znbsec <= 366*rjjss .AND. nleapy <= 1) ) THEN
         clave = '1y'
         IF( iyear2 <= 99 ) THEN 
            WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 999 ) THEN 
            WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 9999 ) THEN 
            WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE
            WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ENDIF
      ELSE
         ! others
         iout = kfreq
         ig = 0
         clout = ''
         IF( iout <= 9 ) THEN 
            ig = 1
            WRITE(clout,'(i1.1)') iout
         ELSE IF( iout <= 99 ) THEN 
            ig = 2
            WRITE(clout,'(i2.2)') iout
         ELSE IF( iout <= 999 ) THEN 
            ig = 3
            WRITE(clout,'(i3.3)') iout
         ELSE IF( iout <= 9999 ) THEN 
            ig = 4
            WRITE(clout,'(i4.4)') iout
         ELSE
            ig = 5
            WRITE(clout,'(i5.5)') iout
         ENDIF
         clave = 'CU'
         IF( iyear2 <= 99 ) THEN 
            WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 999 ) THEN 
            WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE IF( iyear2 <= 9999 ) THEN 
            WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ELSE
            WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
         ENDIF
      ENDIF
#endif
      DO jc = 1, 8
         IF( cexper(jc:jc)==' ' ) go to 120
      END DO
120   CONTINUE
      ic = jc
      clexper = cexper
      IF( jc-1 == 0 ) THEN
         clexper = 'orcafile'
         ic = 9
      ENDIF
      DO jd = 1, 6
         IF( cdsuff(jd:jd) == ' ' ) go to 130
      END DO
130   CONTINUE
      id = jd
      clsuff = cdsuff
      IF( jd-1 == 0 ) THEN
          clsuff = 'output'
          id = 7
      ENDIF
      DO je = 1, 45
        IF( cdfnam(je:je) == ' ' ) go to 140
      END DO
140   CONTINUE
      ie = je
      IF( iout == 0 ) THEN 
         cdfnam=clexper(1:ic-1)//"_"//clave//cdfnam(1:ie-1)//clsuff(1:id-1)
      ELSE 
         cdfnam=clexper(1:ic-1)//"_"//clave//clout(1:ig)//cdfnam(1:ie-1)//clsuff(1:id-1)
      ENDIF
      IF(lwp) WRITE(numout,*) cdfnam     
      IF(lwp) WRITE(numout,*)          

      ! FORMATS

 9001 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
 9002 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
 9003 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
 9004 FORMAT("_",I6.6,2I2.2,"_",I6.6,2I2.2,"_")
 9011 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
 9012 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
 9013 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
 9014 FORMAT("_",I6.6,I2.2,"_",I6.6,I2.2,"_")
 9021 FORMAT("_",I4.4,"_",I4.4,"_")
 9022 FORMAT("_",I4.4,"_",I4.4,"_")
 9023 FORMAT("_",I4.4,"_",I4.4,"_")
 9024 FORMAT("_",I6.6,"_",I6.6,"_")

   END SUBROUTINE dia_nam

   !!======================================================================
END MODULE dianam
