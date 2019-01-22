   SUBROUTINE ctlopn ( knum, cdfile, cdstat, cdform, cdacce,   &
                       klengh, kout, ldwp, krequ )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctlopn  ***
      !!
      !! ** Purpose :   Open file and check if required file is available.
      !!
      !! ** Method  :   Fortan open
      !!
      !! History :
      !!        !  95-12  (G. Madec)  Original code
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT( in ) ::   &
         knum,     & ! logical unit to open
         krequ,    & ! =1 file required (stop if not exist)
         !           ! =0 file not required (create the file if does not exist)
         kout,     & ! number of logical units for write
         klengh      ! record length

      INTEGER ::   iost
      CHARACTER (len=* ), INTENT( in ) ::   &
         cdacce,   & ! access specifier
         cdform,   & ! formatting specifier
         cdstat      ! disposition specifier
      CHARACTER (len=* ), INTENT( in ) ::   &
         cdfile      ! file name to open

      LOGICAL ::  ldwp   ! boolean term for print
      !!----------------------------------------------------------------------
      !!  OPA 9.0 , LOCEAN-IPSL (2005) 
      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/ctlopn.f90,v 1.1.1.1 2005/11/14 10:41:07 opalod Exp $ 
      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
      !!----------------------------------------------------------------------


      ! 1. Required file
      ! ----------------

      IF( krequ == 1 ) THEN

         iost=0
         IF( cdacce(1:6) == 'DIRECT' )  THEN
            OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
               STATUS=cdstat, RECL=klengh, ERR=100, IOSTAT=iost )
         ELSE
            OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
               STATUS=cdstat, ERR=100, IOSTAT=iost)
         ENDIF
         IF( iost == 0 ) THEN
            IF(ldwp) THEN
               WRITE(kout,*) '     file   : ', cdfile,' open ok'
               WRITE(kout,*) '     unit   = ', knum
               WRITE(kout,*) '     status = ', cdstat
               WRITE(kout,*) '     form   = ', cdform
               WRITE(kout,*) '     access = ', cdacce
               WRITE(kout,*)
            ENDIF
         ENDIF
100      CONTINUE
         IF( iost /= 0 ) THEN
            IF(ldwp) THEN
               WRITE(kout,*)
               WRITE(kout,*) ' ===>>>> : bad opening file: ', cdfile
               WRITE(kout,*) ' =======   ===  '
               WRITE(kout,*) '           unit   = ', knum
               WRITE(kout,*) '           status = ', cdstat
               WRITE(kout,*) '           form   = ', cdform
               WRITE(kout,*) '           access = ', cdacce
               WRITE(kout,*) '           iostat = ', iost
               WRITE(kout,*) '           we stop. verify the file '
               WRITE(kout,*)
            ENDIF
            STOP 'ctlopn bad opening'
         ENDIF
         
         
         ! 2. Not required, file create if not exist
         ! -----------------------------------------
         
      ELSEIF( krequ == 0 ) THEN

         iost = 0
         IF( cdacce(1:6) == 'DIRECT' ) THEN
            OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
               STATUS=cdstat, RECL=klengh, ERR=200, IOSTAT=iost )
         ELSE
            OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
               STATUS=cdstat, ERR=200, IOSTAT=iost )
         ENDIF
         IF(iost == 0) THEN
            IF(ldwp) THEN
               WRITE(kout,*) '     file   : ', cdfile,' open ok'
               WRITE(kout,*) '     unit   = ', knum
               WRITE(kout,*) '     status = ', cdstat
               WRITE(kout,*) '     form   = ', cdform
               WRITE(kout,*) '     access = ', cdacce
               WRITE(kout,*)
            ENDIF
         ENDIF
200      CONTINUE
         IF( iost /= 0 ) THEN
            iost = 0
            IF(ldwp) THEN
               WRITE(kout,*)
               WRITE(kout,*) '     ===>>>> : file ', cdfile,   &
                  ' does not exist: it is created'
               WRITE(kout,*) ' =======   ===  '
            ENDIF
            IF( cdacce(1:6) == 'DIRECT' ) THEN
               OPEN( UNIT=knum, FILE=cdfile, FORM=cdform,   &
                  ACCESS=cdacce, STATUS=cdstat,   &
                  RECL=klengh, ERR=210, IOSTAT=iost )
            ELSE
               OPEN( UNIT=knum, FILE=cdfile, FORM=cdform,   &
                  ACCESS=cdacce, STATUS=cdstat, ERR=210,   &
                  IOSTAT=iost )
            ENDIF
            IF(ldwp) THEN
               WRITE(kout,*) '     file   : ', cdfile,' open ok'
               WRITE(kout,*) '     unit   = ', knum
               WRITE(kout,*) '     status = ', cdstat
               WRITE(kout,*) '     form   = ', cdform
               WRITE(kout,*) '     access = ', cdacce
               WRITE(kout,*)
            ENDIF
210         CONTINUE
            IF( iost /= 0 ) THEN
               IF(ldwp) THEN
                  WRITE(kout,*) ' logical unit ',knum,' iostat = ', iost
                  WRITE(kout,*) ' we stop. verify the file ', cdfile
                  WRITE(kout,*)
               ENDIF
               STOP '001'
            ENDIF
         ENDIF
         
      ELSE
         
         IF(ldwp) THEN
            WRITE(kout,*)
            WRITE(kout,*) ' ctlopn : invalid option, krequ = ', krequ
            WRITE(kout,*) ' ~~~~~~   call for file ', cdfile
            WRITE(kout,*)
         ENDIF


         STOP 'ctlopn invalid option'
      ENDIF
      
   END SUBROUTINE ctlopn
