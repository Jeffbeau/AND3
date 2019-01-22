


PGF90 (Version     10.5)          09/23/2010  09:49:27      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir /usr/local/include
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /usr/local/include
          -idir /usr/include/mpich2-x86_64
          -idir /usr/include/mpich2-x86_64
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: ctlopn.F90

(    1)    SUBROUTINE ctlopn ( knum, cdfile, cdstat, cdform, cdacce,   &
(    2)                        klengh, kout, ldwp, krequ )
(    3)       !!----------------------------------------------------------------------
(    4)       !!                  ***  ROUTINE ctlopn  ***
(    5)       !!
(    6)       !! ** Purpose :   Open file and check if required file is available.
(    7)       !!
(    8)       !! ** Method  :   Fortan open
(    9)       !!
(   10)       !! History :
(   11)       !!        !  95-12  (G. Madec)  Original code
(   12)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   13)       !!----------------------------------------------------------------------
(   14) #if defined key_agrif
(   15)       USE Agrif_Util
(   16) #endif
(   17)       IMPLICIT NONE
(   18) 
(   19)       INTEGER, INTENT( in ) ::   &
(   20) #if !defined key_agrif
(   21)          knum,     & ! logical unit to open
(   22) #endif
(   23)          krequ,    & ! =1 file required (stop if not exist)
(   24)          !           ! =0 file not required (create the file if does not exist)
(   25)          kout,     & ! number of logical units for write
(   26)          klengh      ! record length
(   27) 
(   28)       INTEGER ::   iost
(   29)       CHARACTER (len=* ), INTENT( in ) ::   &
(   30)          cdacce,   & ! access specifier
(   31)          cdform,   & ! formatting specifier
(   32)          cdstat      ! disposition specifier
(   33) #if !defined key_agrif
(   34)       CHARACTER (len=* ), INTENT( in ) ::   &
(   35)          cdfile      ! file name to open
(   36) #endif
(   37) 
(   38)       LOGICAL ::  ldwp   ! boolean term for print
(   39) 
(   40) #if defined key_agrif
(   41)       LOGICAL :: op
(   42)       INTEGER, INTENT( inout ) ::   &
(   43)          knum      ! logical unit to open
(   44)       CHARACTER (len=* ), INTENT( inout ) ::   &
(   45)          cdfile      ! file name to open
(   46) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 2

(   47)       !!----------------------------------------------------------------------
(   48)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   49)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ctlopn.F90,v 1.2 2006/03/10 10:55:33 opalod Exp $ 
(   50)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   51)       !!----------------------------------------------------------------------
(   52) 
(   53) #if defined key_agrif
(   54)       if ( .NOT. Agrif_Root() ) then
(   55)          cdfile= TRIM(Agrif_CFixed())//'_'//TRIM(cdfile)
(   56)       endif
(   57)       INQUIRE(UNIT=knum,OPENED=op)
(   58)       IF (op) THEN
(   59)         knum=Agrif_Get_Unit()
(   60)       ENDIF
(   61) #endif
(   62) 
(   63)       ! 1. Required file
(   64)       ! ----------------
(   65) 
(   66)       IF( krequ == 1 ) THEN
(   67) 
(   68)          iost=0
(   69)          IF( cdacce(1:6) == 'DIRECT' )  THEN
(   70)             OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
(   71)                STATUS=cdstat, RECL=klengh, ERR=100, IOSTAT=iost )
(   72)          ELSE
(   73)             OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
(   74)                STATUS=cdstat, ERR=100, IOSTAT=iost)
(   75)          ENDIF
(   76)          IF( iost == 0 ) THEN
(   77)             IF(ldwp) THEN
(   78)                WRITE(kout,*) '     file   : ', cdfile,' open ok'
(   79)                WRITE(kout,*) '     unit   = ', knum
(   80)                WRITE(kout,*) '     status = ', cdstat
(   81)                WRITE(kout,*) '     form   = ', cdform
(   82)                WRITE(kout,*) '     access = ', cdacce
(   83)                WRITE(kout,*)
(   84)             ENDIF
(   85)          ENDIF
(   86) 100      CONTINUE
(   87)          IF( iost /= 0 ) THEN
(   88)             IF(ldwp) THEN
(   89)                WRITE(kout,*)
(   90)                WRITE(kout,*) ' ===>>>> : bad opening file: ', cdfile
(   91)                WRITE(kout,*) ' =======   ===  '
(   92)                WRITE(kout,*) '           unit   = ', knum
(   93)                WRITE(kout,*) '           status = ', cdstat
(   94)                WRITE(kout,*) '           form   = ', cdform
(   95)                WRITE(kout,*) '           access = ', cdacce
(   96)                WRITE(kout,*) '           iostat = ', iost
(   97)                WRITE(kout,*) '           we stop. verify the file '
(   98)                WRITE(kout,*)
(   99)             ENDIF
(  100)             STOP 'ctlopn bad opening'
(  101)          ENDIF
(  102)          
(  103)          
(  104)          ! 2. Not required, file create if not exist






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 3

(  105)          ! -----------------------------------------
(  106)          
(  107)       ELSEIF( krequ == 0 ) THEN
(  108) 
(  109)          iost = 0
(  110)          IF( cdacce(1:6) == 'DIRECT' ) THEN
(  111)             OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
(  112)                STATUS=cdstat, RECL=klengh, ERR=200, IOSTAT=iost )
(  113)          ELSE
(  114)             OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
(  115)                STATUS=cdstat, ERR=200, IOSTAT=iost )
(  116)          ENDIF
(  117)          IF(iost == 0) THEN
(  118)             IF(ldwp) THEN
(  119)                WRITE(kout,*) '     file   : ', cdfile,' open ok'
(  120)                WRITE(kout,*) '     unit   = ', knum
(  121)                WRITE(kout,*) '     status = ', cdstat
(  122)                WRITE(kout,*) '     form   = ', cdform
(  123)                WRITE(kout,*) '     access = ', cdacce
(  124)                WRITE(kout,*)
(  125)             ENDIF
(  126)          ENDIF
(  127) 200      CONTINUE
(  128)          IF( iost /= 0 ) THEN
(  129)             iost = 0
(  130)             IF(ldwp) THEN
(  131)                WRITE(kout,*)
(  132)                WRITE(kout,*) '     ===>>>> : file ', cdfile,   &
(  133)                   ' does not exist: it is created'
(  134)                WRITE(kout,*) ' =======   ===  '
(  135)             ENDIF
(  136)             IF( cdacce(1:6) == 'DIRECT' ) THEN
(  137)                OPEN( UNIT=knum, FILE=cdfile, FORM=cdform,   &
(  138)                   ACCESS=cdacce, STATUS=cdstat,   &
(  139)                   RECL=klengh, ERR=210, IOSTAT=iost )
(  140)             ELSE
(  141)                OPEN( UNIT=knum, FILE=cdfile, FORM=cdform,   &
(  142)                   ACCESS=cdacce, STATUS=cdstat, ERR=210,   &
(  143)                   IOSTAT=iost )
(  144)             ENDIF
(  145)             IF(ldwp) THEN
(  146)                WRITE(kout,*) '     file   : ', cdfile,' open ok'
(  147)                WRITE(kout,*) '     unit   = ', knum
(  148)                WRITE(kout,*) '     status = ', cdstat
(  149)                WRITE(kout,*) '     form   = ', cdform
(  150)                WRITE(kout,*) '     access = ', cdacce
(  151)                WRITE(kout,*)
(  152)             ENDIF
(  153) 210         CONTINUE
(  154)             IF( iost /= 0 ) THEN
(  155)                IF(ldwp) THEN
(  156)                   WRITE(kout,*) ' logical unit ',knum,' iostat = ', iost
(  157)                   WRITE(kout,*) ' we stop. verify the file ', cdfile
(  158)                   WRITE(kout,*)
(  159)                ENDIF
(  160)                STOP '001'
(  161)             ENDIF
(  162)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 4

(  163)          
(  164)       ELSE
(  165)          
(  166)          IF(ldwp) THEN
(  167)             WRITE(kout,*)
(  168)             WRITE(kout,*) ' ctlopn : invalid option, krequ = ', krequ
(  169)             WRITE(kout,*) ' ~~~~~~   call for file ', cdfile
(  170)             WRITE(kout,*)
(  171)          ENDIF
(  172) 
(  173) 
(  174)          STOP 'ctlopn invalid option'
(  175)       ENDIF
(  176)       
(  177)    END SUBROUTINE ctlopn

















































PGF90 (Version     10.5)          09/23/2010  09:49:27      page 5

(    1) # 1 "ctlopn.F90"
(    1)    SUBROUTINE ctlopn ( knum, cdfile, cdstat, cdform, cdacce,   &
(    2)                        klengh, kout, ldwp, krequ )
(    3)       !!----------------------------------------------------------------------
(    4)       !!                  ***  ROUTINE ctlopn  ***
(    5)       !!
(    6)       !! ** Purpose :   Open file and check if required file is available.
(    7)       !!
(    8)       !! ** Method  :   Fortan open
(    9)       !!
(   10)       !! History :
(   11)       !!        !  95-12  (G. Madec)  Original code
(   12)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   13)       !!----------------------------------------------------------------------
(   14) # 17
(   17)       IMPLICIT NONE
(   18) 
(   19)       INTEGER, INTENT( in ) ::   &
(   21)          knum,     & ! logical unit to open
(   23)          krequ,    & ! =1 file required (stop if not exist)
(   24)          !           ! =0 file not required (create the file if does not exist)
(   25)          kout,     & ! number of logical units for write
(   26)          klengh      ! record length
(   27) 
(   28)       INTEGER ::   iost
(   29)       CHARACTER (len=* ), INTENT( in ) ::   &
(   30)          cdacce,   & ! access specifier
(   31)          cdform,   & ! formatting specifier
(   32)          cdstat      ! disposition specifier
(   33) # 34
(   34)       CHARACTER (len=* ), INTENT( in ) ::   &
(   35)          cdfile      ! file name to open
(   36) 
(   37) # 38
(   38)       LOGICAL ::  ldwp   ! boolean term for print
(   39) 
(   40) # 47
(   47)       !!----------------------------------------------------------------------
(   48)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   49)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ctlopn.F90,v 1.2 2006/03/10 10:55:33 opalod Exp $ 
(   50)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   51)       !!----------------------------------------------------------------------
(   52) 
(   53) 
(   54) # 63
(   63)       ! 1. Required file
(   64)       ! ----------------
(   65) 
(   66)       IF( krequ == 1 ) THEN
(   67) 
(   68)          iost=0
(   69)          IF( cdacce(1:6) == 'DIRECT' )  THEN
(   70)             OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
(   71)                STATUS=cdstat, RECL=klengh, ERR=100, IOSTAT=iost )
(   72)          ELSE
(   73)             OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
(   74)                STATUS=cdstat, ERR=100, IOSTAT=iost)
(   75)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 6

(   76)          IF( iost == 0 ) THEN
(   77)             IF(ldwp) THEN
(   78)                WRITE(kout,*) '     file   : ', cdfile,' open ok'
(   79)                WRITE(kout,*) '     unit   = ', knum
(   80)                WRITE(kout,*) '     status = ', cdstat
(   81)                WRITE(kout,*) '     form   = ', cdform
(   82)                WRITE(kout,*) '     access = ', cdacce
(   83)                WRITE(kout,*)
(   84)             ENDIF
(   85)          ENDIF
(   86) 100      CONTINUE
(   87)          IF( iost /= 0 ) THEN
(   88)             IF(ldwp) THEN
(   89)                WRITE(kout,*)
(   90)                WRITE(kout,*) ' ===>>>> : bad opening file: ', cdfile
(   91)                WRITE(kout,*) ' =======   ===  '
(   92)                WRITE(kout,*) '           unit   = ', knum
(   93)                WRITE(kout,*) '           status = ', cdstat
(   94)                WRITE(kout,*) '           form   = ', cdform
(   95)                WRITE(kout,*) '           access = ', cdacce
(   96)                WRITE(kout,*) '           iostat = ', iost
(   97)                WRITE(kout,*) '           we stop. verify the file '
(   98)                WRITE(kout,*)
(   99)             ENDIF
(  100)             STOP 'ctlopn bad opening'
(  101)          ENDIF
(  102)          
(  103)          
(  104)          ! 2. Not required, file create if not exist
(  105)          ! -----------------------------------------
(  106)          
(  107)       ELSEIF( krequ == 0 ) THEN
(  108) 
(  109)          iost = 0
(  110)          IF( cdacce(1:6) == 'DIRECT' ) THEN
(  111)             OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
(  112)                STATUS=cdstat, RECL=klengh, ERR=200, IOSTAT=iost )
(  113)          ELSE
(  114)             OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,   &
(  115)                STATUS=cdstat, ERR=200, IOSTAT=iost )
(  116)          ENDIF
(  117)          IF(iost == 0) THEN
(  118)             IF(ldwp) THEN
(  119)                WRITE(kout,*) '     file   : ', cdfile,' open ok'
(  120)                WRITE(kout,*) '     unit   = ', knum
(  121)                WRITE(kout,*) '     status = ', cdstat
(  122)                WRITE(kout,*) '     form   = ', cdform
(  123)                WRITE(kout,*) '     access = ', cdacce
(  124)                WRITE(kout,*)
(  125)             ENDIF
(  126)          ENDIF
(  127) 200      CONTINUE
(  128)          IF( iost /= 0 ) THEN
(  129)             iost = 0
(  130)             IF(ldwp) THEN
(  131)                WRITE(kout,*)
(  132)                WRITE(kout,*) '     ===>>>> : file ', cdfile,   &
(  133)                   ' does not exist: it is created'






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 7

(  134)                WRITE(kout,*) ' =======   ===  '
(  135)             ENDIF
(  136)             IF( cdacce(1:6) == 'DIRECT' ) THEN
(  137)                OPEN( UNIT=knum, FILE=cdfile, FORM=cdform,   &
(  138)                   ACCESS=cdacce, STATUS=cdstat,   &
(  139)                   RECL=klengh, ERR=210, IOSTAT=iost )
(  140)             ELSE
(  141)                OPEN( UNIT=knum, FILE=cdfile, FORM=cdform,   &
(  142)                   ACCESS=cdacce, STATUS=cdstat, ERR=210,   &
(  143)                   IOSTAT=iost )
(  144)             ENDIF
(  145)             IF(ldwp) THEN
(  146)                WRITE(kout,*) '     file   : ', cdfile,' open ok'
(  147)                WRITE(kout,*) '     unit   = ', knum
(  148)                WRITE(kout,*) '     status = ', cdstat
(  149)                WRITE(kout,*) '     form   = ', cdform
(  150)                WRITE(kout,*) '     access = ', cdacce
(  151)                WRITE(kout,*)
(  152)             ENDIF
(  153) 210         CONTINUE
(  154)             IF( iost /= 0 ) THEN
(  155)                IF(ldwp) THEN
(  156)                   WRITE(kout,*) ' logical unit ',knum,' iostat = ', iost
(  157)                   WRITE(kout,*) ' we stop. verify the file ', cdfile
(  158)                   WRITE(kout,*)
(  159)                ENDIF
(  160)                STOP '001'
(  161)             ENDIF
(  162)          ENDIF
(  163)          
(  164)       ELSE
(  165)          
(  166)          IF(ldwp) THEN
(  167)             WRITE(kout,*)
(  168)             WRITE(kout,*) ' ctlopn : invalid option, krequ = ', krequ
(  169)             WRITE(kout,*) ' ~~~~~~   call for file ', cdfile
(  170)             WRITE(kout,*)
(  171)          ENDIF
(  172) 
(  173) 
(  174)          STOP 'ctlopn invalid option'
(  175)       ENDIF
(  176)       
(  177)    END SUBROUTINE ctlopn
