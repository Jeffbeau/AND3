


PGF90 (Version     10.5)          09/23/2010  09:49:07      page 1

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

Filename: dianam.F90

(    1) MODULE dianam
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dianam  ***
(    4)    !! Ocean diagnostics:  Builds output file name
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dia_nam       : Builds output file name
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE phycst          ! physical constants
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE daymod          ! calendar
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC dia_nam   ! routine called by step.F90
(   21)    !!----------------------------------------------------------------------
(   22)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   23)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/dianam.F90,v 1.3 2005/03/27 18:34:55 opalod Exp $ 
(   24)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   25)    !!----------------------------------------------------------------------
(   26) 
(   27) CONTAINS
(   28) 
(   29)    SUBROUTINE dia_nam( cdfnam, kfreq, cdsuff )
(   30)       !!---------------------------------------------------------------------
(   31)       !!                  ***  ROUTINE dia_nam  ***
(   32)       !!                   
(   33)       !! ** Purpose :   Builds output file name
(   34)       !!
(   35)       !! ** Method  :   File name is a function of date and output frequency
(   36)       !!      cdfnam=<cexper>_<clave>_<idtbeg>_<idtend>_grid_<cdsuff>
(   37)       !!      <clave> = averaging frequency (DA, MO, etc...)
(   38)       !!      <idtbeg>,<idtend> date of beginning and end of run
(   39)       !!
(   40)       !! History :
(   41)       !!        !  99-02  (E. Guilyardi)  Creation for 30 days/month
(   42)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   43)       !!----------------------------------------------------------------------
(   44)       !! * Arguments
(   45)       CHARACTER (len=*), INTENT( out ) ::   cdfnam   ! file name
(   46)       CHARACTER (len=*), INTENT( in  ) ::   cdsuff   ! ???






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 2

(   47)       INTEGER,           INTENT( in  ) ::   kfreq    ! ???
(   48) 
(   49)       !! * Local declarations
(   50)       CHARACTER (len=8) ::   clexper
(   51)       CHARACTER (len=2) ::   clave
(   52)       CHARACTER (len=5) ::   clout
(   53)       CHARACTER (len=6) ::   clsuff
(   54)       INTEGER :: jt, jc, jd, je            ! dummy loop indices
(   55)       INTEGER ::   &
(   56)          ic, id, ie, ig, ijjmm, iout,   &  ! temporary integers
(   57)          iyear1, imonth1, iday1,        &  !    "          "
(   58)          iyear2, imonth2, iday2            !    "          "
(   59)       REAL(wp) ::   &
(   60)          z5j, znbsec, zdate1, zdate2, zdrun, zdt   ! temporary scalars
(   61)       !!----------------------------------------------------------------------
(   62) 
(   63)       IF(lwp) WRITE(numout,*)
(   64)       IF(lwp) WRITE(numout,*) ' dia_nam: building output file name'
(   65)       IF(lwp) WRITE(numout,*) ' ~~~~~~~'
(   66)       IF(lwp) WRITE(numout,*)
(   67) 
(   68)       ! 0. Initialisation
(   69)       ! -----------------
(   70) 
(   71)       cdfnam = ''
(   72) 
(   73)       !    number of seconds of the run
(   74) 
(   75)       z5j = 5*rjjss
(   76)       zdt = rdt
(   77)       IF( nacc == 1 ) zdt = rdtmin
(   78)       zdrun = FLOAT( nitend - nit000 ) * zdt
(   79) 
(   80)       !  date of beginning of run
(   81) 
(   82)       iyear1  = ndastp/10000
(   83)       imonth1 = ndastp/100 - iyear1*100
(   84)       iday1   = ndastp - imonth1*100 - iyear1*10000
(   85)       IF( nleapy == 1) THEN 
(   86)          ijjmm=0
(   87)          IF( MOD( iyear1, 4 ) == 0 ) THEN
(   88)             DO jt = 1, imonth1-1
(   89)                ijjmm = ijjmm + nbiss(jt)
(   90)             END DO
(   91)          ELSE
(   92)             DO jt = 1, imonth1-1
(   93)                ijjmm = ijjmm + nobis(jt)
(   94)             END DO
(   95)          ENDIF
(   96)          ijjmm = ijjmm + (iyear1-1)/4
(   97)          zdate1 = ( (iyear1-1)*365 + ijjmm +iday1-1 ) * rjjss   
(   98)       ELSE IF( nleapy == 0 ) THEN
(   99)          ijjmm = 0
(  100)          DO jt = 1, imonth1-1
(  101)             ijjmm = ijjmm + nobis(jt)
(  102)          END DO
(  103)          zdate1 = ( (iyear1-1)*raajj + ijjmm + iday1-1)* rjjss
(  104)       ELSE 






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 3

(  105)          zdate1 = ( (iyear1-1)*nleapy*raamo + (imonth1-1)*nleapy + iday1-1)* rjjss
(  106)       ENDIF
(  107) 
(  108)       !  date of end of run (= date of beginning of next run)
(  109) 
(  110)       zdate2 = zdate1 + zdrun
(  111)       IF( nleapy == 1 ) THEN 
(  112)          iyear2 = zdate2/(365.25*rjjss)+1
(  113)          ijjmm = INT(zdate2/rjjss)-365*(iyear2-1)-(iyear2-1)/4
(  114)          IF( ijjmm < 0 ) THEN
(  115)             iyear2 = iyear2-1
(  116)             ijjmm = zdate2/rjjss-365.*(iyear2-1)-(iyear2-1)/4
(  117)          ENDIF
(  118)          IF( MOD( iyear2, 4 ) == 0 ) THEN
(  119)             DO jt = 1, 12
(  120)                ijjmm = ijjmm - nbiss(jt)
(  121)                IF( ijjmm <= 0 ) go to 10
(  122)             END DO
(  123)             jt = 12
(  124) 10          CONTINUE
(  125)             imonth2 = jt
(  126)             ijjmm = 0
(  127)             DO jt = 1, jt-1
(  128)                ijjmm = ijjmm + nbiss(jt)
(  129)             END DO
(  130)          ELSE
(  131)             DO jt = 1, 12
(  132)                ijjmm = ijjmm - nobis(jt)
(  133)                IF( ijjmm <= 0 ) go to 15
(  134)             END DO
(  135)             jt = 12
(  136) 15          CONTINUE
(  137)             imonth2 = jt
(  138)             ijjmm = 0
(  139)             DO jt = 1, jt-1
(  140)                ijjmm = ijjmm + nobis(jt)
(  141)             END DO
(  142)          ENDIF
(  143)          iday2 = zdate2/rjjss-365.*(iyear2-1)-ijjmm+1-(iyear2-1)/4     
(  144)       ELSE IF( nleapy == 0 ) THEN
(  145)          iyear2 = zdate2/raass+1
(  146)          ijjmm  = zdate2/rjjss-raajj*(iyear2-1)
(  147)          DO jt = 1, 12
(  148)             ijjmm = ijjmm - nobis(jt)
(  149)             IF(ijjmm <= 0) go to 20
(  150)          END DO
(  151)          jt = 12
(  152) 20       CONTINUE
(  153)          imonth2 = jt
(  154)          ijjmm = 0
(  155)          DO jt = 1, jt-1
(  156)             ijjmm = ijjmm + nobis(jt)
(  157)          END DO
(  158)          iday2 = zdate2/rjjss-raajj*(iyear2-1)-ijjmm+1          
(  159)       ELSE 
(  160)          zdate2 = zdate2 / rjjss
(  161)          imonth2 = zdate2/FLOAT(nleapy)
(  162)          iday2 = zdate2 - imonth2*FLOAT(nleapy) + 1.






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 4

(  163)          iyear2 = imonth2/12
(  164)          imonth2 = imonth2 - iyear2*12
(  165)          imonth2 = imonth2 + 1
(  166)          iyear2 = iyear2 + 1
(  167)          IF( iday2 == 0 ) THEN
(  168)             iday2 = nleapy
(  169)             imonth2 = imonth2 - 1
(  170)             IF( imonth2 == 0 ) THEN
(  171)                imonth2 = 12
(  172)                iyear2 = iyear2 - 1
(  173)             ENDIF
(  174)          ENDIF
(  175)       ENDIF
(  176) 
(  177) 
(  178)       ! 1. Define time averaging period <nn><type>
(  179)       !    ---------------------------------------
(  180) 
(  181)       iout = 0
(  182) #if defined key_diainstant
(  183)       clave = 'IN'
(  184)       IF( iyear2 <= 99 ) THEN 
(  185)          WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  186)       ELSE IF( iyear2 <= 999 ) THEN 
(  187)          WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  188)       ELSE IF( iyear2 <= 9999 ) THEN 
(  189)          WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  190)       ELSE
(  191)          WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  192)       ENDIF
(  193) #else
(  194) 
(  195)       znbsec=kfreq*zdt
(  196)       ! daily output
(  197)       IF( znbsec == rjjss ) THEN
(  198)          clave = '1d'
(  199)          IF( iyear2 <= 99 ) THEN 
(  200)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  201)          ELSE IF( iyear2 <= 999 ) THEN 
(  202)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  203)          ELSE IF( iyear2 <= 9999 ) THEN 
(  204)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  205)          ELSE
(  206)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  207)          ENDIF
(  208)          ! 5 day output 
(  209)       ELSE IF( znbsec == z5j ) THEN
(  210)          clave='5d'
(  211)          IF( iyear2 <= 99 ) THEN 
(  212)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  213)          ELSE IF( iyear2 <= 999 ) THEN 
(  214)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  215)          ELSE IF( iyear2 <= 9999 ) THEN 
(  216)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  217)          ELSE
(  218)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  219)          ENDIF
(  220)          ! monthly ouput 






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 5

(  221)       ELSE IF( (znbsec == rmoss .AND. nleapy > 1) .OR.   &
(  222)                (znbsec >= 28*rjjss .AND. znbsec <= 31*rjjss .AND. nleapy <= 1) ) THEN
(  223)          clave = '1m'
(  224)          IF( iyear2 <= 99 ) THEN 
(  225)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  226)          ELSE IF( iyear2 <= 999 ) THEN 
(  227)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  228)          ELSE IF( iyear2 <= 9999 ) THEN 
(  229)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  230)          ELSE
(  231)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  232)          ENDIF
(  233)          ! annual output
(  234)       ELSE IF( (znbsec == raass .AND. nleapy > 1) .OR.   &
(  235)                (znbsec >= 365*rjjss .AND. znbsec <= 366*rjjss .AND. nleapy <= 1) ) THEN
(  236)          clave = '1y'
(  237)          IF( iyear2 <= 99 ) THEN 
(  238)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  239)          ELSE IF( iyear2 <= 999 ) THEN 
(  240)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  241)          ELSE IF( iyear2 <= 9999 ) THEN 
(  242)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  243)          ELSE
(  244)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  245)          ENDIF
(  246)       ELSE
(  247)          ! others
(  248)          iout = kfreq
(  249)          ig = 0
(  250)          clout = ''
(  251)          IF( iout <= 9 ) THEN 
(  252)             ig = 1
(  253)             WRITE(clout,'(i1.1)') iout
(  254)          ELSE IF( iout <= 99 ) THEN 
(  255)             ig = 2
(  256)             WRITE(clout,'(i2.2)') iout
(  257)          ELSE IF( iout <= 999 ) THEN 
(  258)             ig = 3
(  259)             WRITE(clout,'(i3.3)') iout
(  260)          ELSE IF( iout <= 9999 ) THEN 
(  261)             ig = 4
(  262)             WRITE(clout,'(i4.4)') iout
(  263)          ELSE
(  264)             ig = 5
(  265)             WRITE(clout,'(i5.5)') iout
(  266)          ENDIF
(  267)          clave = 'CU'
(  268)          IF( iyear2 <= 99 ) THEN 
(  269)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  270)          ELSE IF( iyear2 <= 999 ) THEN 
(  271)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  272)          ELSE IF( iyear2 <= 9999 ) THEN 
(  273)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  274)          ELSE
(  275)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  276)          ENDIF
(  277)       ENDIF
(  278) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 6

(  279)       DO jc = 1, 8
(  280)          IF( cexper(jc:jc)==' ' ) go to 120
(  281)       END DO
(  282) 120   CONTINUE
(  283)       ic = jc
(  284)       clexper = cexper
(  285)       IF( jc-1 == 0 ) THEN
(  286)          clexper = 'orcafile'
(  287)          ic = 9
(  288)       ENDIF
(  289)       DO jd = 1, 6
(  290)          IF( cdsuff(jd:jd) == ' ' ) go to 130
(  291)       END DO
(  292) 130   CONTINUE
(  293)       id = jd
(  294)       clsuff = cdsuff
(  295)       IF( jd-1 == 0 ) THEN
(  296)           clsuff = 'output'
(  297)           id = 7
(  298)       ENDIF
(  299)       DO je = 1, 45
(  300)         IF( cdfnam(je:je) == ' ' ) go to 140
(  301)       END DO
(  302) 140   CONTINUE
(  303)       ie = je
(  304)       IF( iout == 0 ) THEN 
(  305)          cdfnam=clexper(1:ic-1)//"_"//clave//cdfnam(1:ie-1)//clsuff(1:id-1)
(  306)       ELSE 
(  307)          cdfnam=clexper(1:ic-1)//"_"//clave//clout(1:ig)//cdfnam(1:ie-1)//clsuff(1:id-1)
(  308)       ENDIF
(  309)       IF(lwp) WRITE(numout,*) cdfnam     
(  310)       IF(lwp) WRITE(numout,*)          
(  311) 
(  312)       ! FORMATS
(  313) 
(  314)  9001 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
(  315)  9002 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
(  316)  9003 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
(  317)  9004 FORMAT("_",I6.6,2I2.2,"_",I6.6,2I2.2,"_")
(  318)  9011 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
(  319)  9012 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
(  320)  9013 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
(  321)  9014 FORMAT("_",I6.6,I2.2,"_",I6.6,I2.2,"_")
(  322)  9021 FORMAT("_",I4.4,"_",I4.4,"_")
(  323)  9022 FORMAT("_",I4.4,"_",I4.4,"_")
(  324)  9023 FORMAT("_",I4.4,"_",I4.4,"_")
(  325)  9024 FORMAT("_",I6.6,"_",I6.6,"_")
(  326) 
(  327)    END SUBROUTINE dia_nam
(  328) 
(  329)    !!======================================================================
(  330) END MODULE dianam












PGF90 (Version     10.5)          09/23/2010  09:49:07      page 7

(    1) # 1 "dianam.F90"
(    1) MODULE dianam
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dianam  ***
(    4)    !! Ocean diagnostics:  Builds output file name
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dia_nam       : Builds output file name
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE phycst          ! physical constants
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE daymod          ! calendar
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC dia_nam   ! routine called by step.F90
(   21)    !!----------------------------------------------------------------------
(   22)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   23)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/dianam.F90,v 1.3 2005/03/27 18:34:55 opalod Exp $ 
(   24)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   25)    !!----------------------------------------------------------------------
(   26) 
(   27) CONTAINS
(   28) 
(   29)    SUBROUTINE dia_nam( cdfnam, kfreq, cdsuff )
(   30)       !!---------------------------------------------------------------------
(   31)       !!                  ***  ROUTINE dia_nam  ***
(   32)       !!                   
(   33)       !! ** Purpose :   Builds output file name
(   34)       !!
(   35)       !! ** Method  :   File name is a function of date and output frequency
(   36)       !!      cdfnam=<cexper>_<clave>_<idtbeg>_<idtend>_grid_<cdsuff>
(   37)       !!      <clave> = averaging frequency (DA, MO, etc...)
(   38)       !!      <idtbeg>,<idtend> date of beginning and end of run
(   39)       !!
(   40)       !! History :
(   41)       !!        !  99-02  (E. Guilyardi)  Creation for 30 days/month
(   42)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   43)       !!----------------------------------------------------------------------
(   44)       !! * Arguments
(   45)       CHARACTER (len=*), INTENT( out ) ::   cdfnam   ! file name
(   46)       CHARACTER (len=*), INTENT( in  ) ::   cdsuff   ! ???
(   47)       INTEGER,           INTENT( in  ) ::   kfreq    ! ???
(   48) 
(   49)       !! * Local declarations
(   50)       CHARACTER (len=8) ::   clexper
(   51)       CHARACTER (len=2) ::   clave
(   52)       CHARACTER (len=5) ::   clout
(   53)       CHARACTER (len=6) ::   clsuff
(   54)       INTEGER :: jt, jc, jd, je            ! dummy loop indices
(   55)       INTEGER ::   &
(   56)          ic, id, ie, ig, ijjmm, iout,   &  ! temporary integers
(   57)          iyear1, imonth1, iday1,        &  !    "          "






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 8

(   58)          iyear2, imonth2, iday2            !    "          "
(   59)       REAL(wp) ::   &
(   60)          z5j, znbsec, zdate1, zdate2, zdrun, zdt   ! temporary scalars
(   61)       !!----------------------------------------------------------------------
(   62) 
(   63)       IF(lwp) WRITE(numout,*)
(   64)       IF(lwp) WRITE(numout,*) ' dia_nam: building output file name'
(   65)       IF(lwp) WRITE(numout,*) ' ~~~~~~~'
(   66)       IF(lwp) WRITE(numout,*)
(   67) 
(   68)       ! 0. Initialisation
(   69)       ! -----------------
(   70) 
(   71)       cdfnam = ''
(   72) 
(   73)       !    number of seconds of the run
(   74) 
(   75)       z5j = 5*rjjss
(   76)       zdt = rdt
(   77)       IF( nacc == 1 ) zdt = rdtmin
(   78)       zdrun = FLOAT( nitend - nit000 ) * zdt
(   79) 
(   80)       !  date of beginning of run
(   81) 
(   82)       iyear1  = ndastp/10000
(   83)       imonth1 = ndastp/100 - iyear1*100
(   84)       iday1   = ndastp - imonth1*100 - iyear1*10000
(   85)       IF( nleapy == 1) THEN 
(   86)          ijjmm=0
(   87)          IF( MOD( iyear1, 4 ) == 0 ) THEN
(   88)             DO jt = 1, imonth1-1
(   89)                ijjmm = ijjmm + nbiss(jt)
(   90)             END DO
(   91)          ELSE
(   92)             DO jt = 1, imonth1-1
(   93)                ijjmm = ijjmm + nobis(jt)
(   94)             END DO
(   95)          ENDIF
(   96)          ijjmm = ijjmm + (iyear1-1)/4
(   97)          zdate1 = ( (iyear1-1)*365 + ijjmm +iday1-1 ) * rjjss   
(   98)       ELSE IF( nleapy == 0 ) THEN
(   99)          ijjmm = 0
(  100)          DO jt = 1, imonth1-1
(  101)             ijjmm = ijjmm + nobis(jt)
(  102)          END DO
(  103)          zdate1 = ( (iyear1-1)*raajj + ijjmm + iday1-1)* rjjss
(  104)       ELSE 
(  105)          zdate1 = ( (iyear1-1)*nleapy*raamo + (imonth1-1)*nleapy + iday1-1)* rjjss
(  106)       ENDIF
(  107) 
(  108)       !  date of end of run (= date of beginning of next run)
(  109) 
(  110)       zdate2 = zdate1 + zdrun
(  111)       IF( nleapy == 1 ) THEN 
(  112)          iyear2 = zdate2/(365.25*rjjss)+1
(  113)          ijjmm = INT(zdate2/rjjss)-365*(iyear2-1)-(iyear2-1)/4
(  114)          IF( ijjmm < 0 ) THEN
(  115)             iyear2 = iyear2-1






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 9

(  116)             ijjmm = zdate2/rjjss-365.*(iyear2-1)-(iyear2-1)/4
(  117)          ENDIF
(  118)          IF( MOD( iyear2, 4 ) == 0 ) THEN
(  119)             DO jt = 1, 12
(  120)                ijjmm = ijjmm - nbiss(jt)
(  121)                IF( ijjmm <= 0 ) go to 10
(  122)             END DO
(  123)             jt = 12
(  124) 10          CONTINUE
(  125)             imonth2 = jt
(  126)             ijjmm = 0
(  127)             DO jt = 1, jt-1
(  128)                ijjmm = ijjmm + nbiss(jt)
(  129)             END DO
(  130)          ELSE
(  131)             DO jt = 1, 12
(  132)                ijjmm = ijjmm - nobis(jt)
(  133)                IF( ijjmm <= 0 ) go to 15
(  134)             END DO
(  135)             jt = 12
(  136) 15          CONTINUE
(  137)             imonth2 = jt
(  138)             ijjmm = 0
(  139)             DO jt = 1, jt-1
(  140)                ijjmm = ijjmm + nobis(jt)
(  141)             END DO
(  142)          ENDIF
(  143)          iday2 = zdate2/rjjss-365.*(iyear2-1)-ijjmm+1-(iyear2-1)/4     
(  144)       ELSE IF( nleapy == 0 ) THEN
(  145)          iyear2 = zdate2/raass+1
(  146)          ijjmm  = zdate2/rjjss-raajj*(iyear2-1)
(  147)          DO jt = 1, 12
(  148)             ijjmm = ijjmm - nobis(jt)
(  149)             IF(ijjmm <= 0) go to 20
(  150)          END DO
(  151)          jt = 12
(  152) 20       CONTINUE
(  153)          imonth2 = jt
(  154)          ijjmm = 0
(  155)          DO jt = 1, jt-1
(  156)             ijjmm = ijjmm + nobis(jt)
(  157)          END DO
(  158)          iday2 = zdate2/rjjss-raajj*(iyear2-1)-ijjmm+1          
(  159)       ELSE 
(  160)          zdate2 = zdate2 / rjjss
(  161)          imonth2 = zdate2/FLOAT(nleapy)
(  162)          iday2 = zdate2 - imonth2*FLOAT(nleapy) + 1.
(  163)          iyear2 = imonth2/12
(  164)          imonth2 = imonth2 - iyear2*12
(  165)          imonth2 = imonth2 + 1
(  166)          iyear2 = iyear2 + 1
(  167)          IF( iday2 == 0 ) THEN
(  168)             iday2 = nleapy
(  169)             imonth2 = imonth2 - 1
(  170)             IF( imonth2 == 0 ) THEN
(  171)                imonth2 = 12
(  172)                iyear2 = iyear2 - 1
(  173)             ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 10

(  174)          ENDIF
(  175)       ENDIF
(  176) 
(  177) 
(  178)       ! 1. Define time averaging period <nn><type>
(  179)       !    ---------------------------------------
(  180) 
(  181)       iout = 0
(  182) 
(  183) # 195
(  195)       znbsec=kfreq*zdt
(  196)       ! daily output
(  197)       IF( znbsec == rjjss ) THEN
(  198)          clave = '1d'
(  199)          IF( iyear2 <= 99 ) THEN 
(  200)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  201)          ELSE IF( iyear2 <= 999 ) THEN 
(  202)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  203)          ELSE IF( iyear2 <= 9999 ) THEN 
(  204)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  205)          ELSE
(  206)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  207)          ENDIF
(  208)          ! 5 day output 
(  209)       ELSE IF( znbsec == z5j ) THEN
(  210)          clave='5d'
(  211)          IF( iyear2 <= 99 ) THEN 
(  212)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  213)          ELSE IF( iyear2 <= 999 ) THEN 
(  214)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  215)          ELSE IF( iyear2 <= 9999 ) THEN 
(  216)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  217)          ELSE
(  218)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  219)          ENDIF
(  220)          ! monthly ouput 
(  221)       ELSE IF( (znbsec == rmoss .AND. nleapy > 1) .OR.   &
(  222)                (znbsec >= 28*rjjss .AND. znbsec <= 31*rjjss .AND. nleapy <= 1) ) THEN
(  223)          clave = '1m'
(  224)          IF( iyear2 <= 99 ) THEN 
(  225)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  226)          ELSE IF( iyear2 <= 999 ) THEN 
(  227)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  228)          ELSE IF( iyear2 <= 9999 ) THEN 
(  229)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  230)          ELSE
(  231)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  232)          ENDIF
(  233)          ! annual output
(  234)       ELSE IF( (znbsec == raass .AND. nleapy > 1) .OR.   &
(  235)                (znbsec >= 365*rjjss .AND. znbsec <= 366*rjjss .AND. nleapy <= 1) ) THEN
(  236)          clave = '1y'
(  237)          IF( iyear2 <= 99 ) THEN 
(  238)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  239)          ELSE IF( iyear2 <= 999 ) THEN 
(  240)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  241)          ELSE IF( iyear2 <= 9999 ) THEN 
(  242)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 11

(  243)          ELSE
(  244)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  245)          ENDIF
(  246)       ELSE
(  247)          ! others
(  248)          iout = kfreq
(  249)          ig = 0
(  250)          clout = ''
(  251)          IF( iout <= 9 ) THEN 
(  252)             ig = 1
(  253)             WRITE(clout,'(i1.1)') iout
(  254)          ELSE IF( iout <= 99 ) THEN 
(  255)             ig = 2
(  256)             WRITE(clout,'(i2.2)') iout
(  257)          ELSE IF( iout <= 999 ) THEN 
(  258)             ig = 3
(  259)             WRITE(clout,'(i3.3)') iout
(  260)          ELSE IF( iout <= 9999 ) THEN 
(  261)             ig = 4
(  262)             WRITE(clout,'(i4.4)') iout
(  263)          ELSE
(  264)             ig = 5
(  265)             WRITE(clout,'(i5.5)') iout
(  266)          ENDIF
(  267)          clave = 'CU'
(  268)          IF( iyear2 <= 99 ) THEN 
(  269)             WRITE(cdfnam,9001) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  270)          ELSE IF( iyear2 <= 999 ) THEN 
(  271)             WRITE(cdfnam,9002) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  272)          ELSE IF( iyear2 <= 9999 ) THEN 
(  273)             WRITE(cdfnam,9003) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  274)          ELSE
(  275)             WRITE(cdfnam,9004) iyear1,imonth1,iday1,iyear2,imonth2,iday2
(  276)          ENDIF
(  277)       ENDIF
(  278) # 279
(  279)       DO jc = 1, 8
(  280)          IF( cexper(jc:jc)==' ' ) go to 120
(  281)       END DO
(  282) 120   CONTINUE
(  283)       ic = jc
(  284)       clexper = cexper
(  285)       IF( jc-1 == 0 ) THEN
(  286)          clexper = 'orcafile'
(  287)          ic = 9
(  288)       ENDIF
(  289)       DO jd = 1, 6
(  290)          IF( cdsuff(jd:jd) == ' ' ) go to 130
(  291)       END DO
(  292) 130   CONTINUE
(  293)       id = jd
(  294)       clsuff = cdsuff
(  295)       IF( jd-1 == 0 ) THEN
(  296)           clsuff = 'output'
(  297)           id = 7
(  298)       ENDIF
(  299)       DO je = 1, 45
(  300)         IF( cdfnam(je:je) == ' ' ) go to 140






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 12

(  301)       END DO
(  302) 140   CONTINUE
(  303)       ie = je
(  304)       IF( iout == 0 ) THEN 
(  305)          cdfnam=clexper(1:ic-1)//"_"//clave//cdfnam(1:ie-1)//clsuff(1:id-1)
(  306)       ELSE 
(  307)          cdfnam=clexper(1:ic-1)//"_"//clave//clout(1:ig)//cdfnam(1:ie-1)//clsuff(1:id-1)
(  308)       ENDIF
(  309)       IF(lwp) WRITE(numout,*) cdfnam     
(  310)       IF(lwp) WRITE(numout,*)          
(  311) 
(  312)       ! FORMATS
(  313) 
(  314)  9001 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
(  315)  9002 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
(  316)  9003 FORMAT("_",I4.4,2I2.2,"_",I4.4,2I2.2,"_")
(  317)  9004 FORMAT("_",I6.6,2I2.2,"_",I6.6,2I2.2,"_")
(  318)  9011 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
(  319)  9012 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
(  320)  9013 FORMAT("_",I4.4,I2.2,"_",I4.4,I2.2,"_")
(  321)  9014 FORMAT("_",I6.6,I2.2,"_",I6.6,I2.2,"_")
(  322)  9021 FORMAT("_",I4.4,"_",I4.4,"_")
(  323)  9022 FORMAT("_",I4.4,"_",I4.4,"_")
(  324)  9023 FORMAT("_",I4.4,"_",I4.4,"_")
(  325)  9024 FORMAT("_",I6.6,"_",I6.6,"_")
(  326) 
(  327)    END SUBROUTINE dia_nam
(  328) 
(  329)    !!======================================================================
(  330) END MODULE dianam
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (dianam.F90: 78)
