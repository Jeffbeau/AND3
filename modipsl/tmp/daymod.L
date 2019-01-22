


PGF90 (Version     10.5)          09/23/2010  09:49:01      page 1

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

Filename: daymod.F90

(    1) MODULE daymod
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  daymod  ***
(    4)    !! Ocean        :  calendar 
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   day        : calendar
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE phycst          ! physical constants
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE prtctl          ! Print control
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC day        ! called by step.F90
(   21) 
(   22)    !! * Shared module variables
(   23)    INTEGER , PUBLIC ::   &  !:
(   24)       nyear     ,   &  !: current year
(   25)       nmonth    ,   &  !: current month
(   26)       nday      ,   &  !: current day of the month
(   27)       nday_year ,   &  !: curent day counted from jan 1st of the current year
(   28)       ndastp           !: time step date in year/month/day aammjj
(   29)    REAL(wp), PUBLIC ::   &  !:
(   30)        model_time, &   !AD: decimal version of ndastp date for output files
(   31)        adatrj   ,   &  !: number of elapsed days since the begining of the run
(   32)        adatrj0         !: value of adatrj at nit000-1 (before the present run).
(   33)        !               !  it is the accumulated duration of previous runs
(   34)        !               !  that may have been run with different time steps.
(   35)    !!----------------------------------------------------------------------
(   36)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/daymod.F90,v 1.5 2005/09/02 15:45:19 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE day( kt )
(   44)       !!----------------------------------------------------------------------
(   45)       !!                      ***  ROUTINE day  ***
(   46)       !! 






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 2

(   47)       !! ** Purpose :   Compute the date with a day iteration IF necessary.
(   48)       !!
(   49)       !! ** Method  : - ???
(   50)       !!
(   51)       !! ** Action  : - nyear     : current year
(   52)       !!              - nmonth    : current month of the year nyear
(   53)       !!              - nday      : current day of the month nmonth
(   54)       !!              - nday_year : current day of the year nyear
(   55)       !!              - ndastp    : =nyear*10000+nmonth*100+nday
(   56)       !!              - adatrj    : date in days since the beginning of the run
(   57)       !!
(   58)       !! History :
(   59)       !!        !  94-09  (M. Pontaud M. Imbard)  Original code
(   60)       !!        !  97-03  (O. Marti)
(   61)       !!        !  97-05  (G. Madec) 
(   62)       !!        !  97-08  (M. Imbard)
(   63)       !!   9.0  !  03-09  (G. Madec)  F90 + nyear, nmonth, nday
(   64)       !!        !  04-01  (A.M. Treguier) new calculation based on adatrj
(   65)       !!----------------------------------------------------------------------      
(   66)       !! * Arguments
(   67)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step indices
(   68) 
(   69)       !! * Local declarations
(   70)       INTEGER  ::   js                   ! dummy loop indice
(   71)       INTEGER  ::   iend, iday0, iday1   ! temporary integers
(   72)       REAL(wp) :: zadatrjn, zadatrjb     ! adatrj at timestep kt-1 and kt-2 
(   73)       CHARACTER (len=25) :: charout
(   74)       !!----------------------------------------------------------------------
(   75) 
(   76)       ! 0.  initialization of adatrj0 and nday, nmonth,nyear, nday_year.
(   77)       !     ndastp has been initialized in domain.F90 or restart.F90
(   78)       !-----------------------------------------------------------------
(   79) 
(   80)       IF( kt == nit000 ) THEN
(   81) 
(   82)          IF( .NOT.ln_rstart )   adatrj0 = 0.e0      ! adatrj0 initialized in rst_read when restart 
(   83) 
(   84)          adatrj  = adatrj0
(   85)          nyear   =   ndastp / 10000
(   86)          nmonth  = ( ndastp - (nyear * 10000) ) / 100
(   87)          nday    =   ndastp - (nyear * 10000) - ( nmonth * 100 ) 
(   88) 
(   89)          ! Calculates nday_year, day since january 1st (useful to read  daily forcing fields)
(   90)          nday_year =  nday
(   91)          !                               ! accumulates days of previous months of this year
(   92)          DO js = 1, nmonth-1
(   93)             IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
(   94)                nday_year = nday_year + nbiss(js)
(   95)             ELSE
(   96)                nday_year = nday_year + nobis(js)
(   97)             ENDIF
(   98)          END DO
(   99) 
(  100)       ENDIF
(  101) 
(  102)       ! I.  calculates adatrj, zadatrjn, zadatrjb.
(  103)       ! ------------------------------------------------------------------
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 3

(  105)       adatrj    = adatrj0 + ( kt - nit000 + 1 ) * rdttra(1) / rday
(  106)       zadatrjn  = adatrj0 + ( kt - nit000     ) * rdttra(1) / rday
(  107)       zadatrjb  = adatrj0 + ( kt - nit000 - 1 ) * rdttra(1) / rday
(  108) 
(  109) 
(  110)       ! II.  increment the date.  The date corresponds to 'now' variables (kt-1),
(  111)       !      which is the time step of forcing fields. 
(  112)       !      Do not do this at nit000  unless nrstdt= 2
(  113)       !      In that case ndastp (read in restart) was for step nit000-2
(  114)       ! -------------------------------------------------------------------
(  115) 
(  116)       iday0 = INT( zadatrjb )
(  117)       iday1 = INT( zadatrjn )
(  118) 
(  119)       IF( iday1 - iday0 >= 1 .AND. ( kt /= nit000 .OR. nrstdt == 2 ) ) THEN
(  120) 
(  121)          ! increase calendar
(  122)          nyear  =   ndastp / 10000
(  123)          nmonth = ( ndastp - (nyear * 10000) ) / 100
(  124)          nday   =   ndastp - (nyear * 10000) - ( nmonth * 100 ) 
(  125)          nday = nday + 1
(  126)          IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
(  127)             iend = nbiss(nmonth)
(  128)          ELSEIF( nleapy > 1 ) THEN 
(  129)             iend = nleapy
(  130)          ELSE 
(  131)             iend = nobis(nmonth)
(  132)          ENDIF
(  133)          IF( nday == iend + 1 ) THEN
(  134)             nday  = 1
(  135)             nmonth = nmonth + 1
(  136)             IF( nmonth == 13 ) THEN
(  137)                nmonth  = 1
(  138)                nyear = nyear + 1
(  139)             ENDIF
(  140)          ENDIF
(  141)          ndastp = nyear * 10000 + nmonth * 100 + nday
(  142) 
(  143)          ! Calculates nday_year, day since january 1st (useful to read  daily forcing fields)
(  144)          nday_year =  nday
(  145)          !                                ! accumulates days of previous months of this year
(  146)          DO js = 1, nmonth-1
(  147)             IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
(  148)                nday_year = nday_year + nbiss(js)
(  149)             ELSE
(  150)                nday_year = nday_year + nobis(js)
(  151)             ENDIF
(  152)          END DO
(  153) 
(  154)          IF(lwp) WRITE(numout,*)' ==============>> time-step =', kt, ' New day, DATE= ',   &
(  155)             &                   nyear, '/', nmonth, '/', nday, 'nday_year:', nday_year
(  156)       ENDIF
(  157) 
(  158)       !AD: NOTE: implemented to have decimal date in output file. Decimal day is relative to 
(  159)       !          ndate0( which is an integer date only  i.e starts at 00 ZZ). So for example
(  160)       !          20051225.5 is noon of december 25, 2005. In future a decimal date could be
(  161)       !          added to ndate0 (i.e convert it to a real). In this case the decimal date
(  162)       !          Here will have to modified to account for the offset (29 Sept 2009)






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 4

(  163) 
(  164)       !    NOTE: The date (ndastp) here is for the now model time (kt-1)*rdt should 
(  165)       !           techincally be for end of of step time (kt)*rdt. 
(  166)       !           This will give a adate of 20010101.99 for
(  167)       !           the first output (nwrite=180,rdt=480,nit000=1) at kt=180.     
(  168)       !           To avoid this we add at the time of one time step
(  169) 
(  170)       model_time=FLOAT(ndastp)+(zadatrjn-FLOAT(INT(zadatrjn))) + rdttra(1) / rday 
(  171) 
(  172)       IF(ln_ctl) THEN
(  173)          WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
(  174)          CALL prt_ctl_info(charout)
(  175)       ENDIF
(  176) 
(  177)    END SUBROUTINE day
(  178) 
(  179)    !!======================================================================
(  180) END MODULE daymod














































PGF90 (Version     10.5)          09/23/2010  09:49:01      page 5

(    1) # 1 "daymod.F90"
(    1) MODULE daymod
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  daymod  ***
(    4)    !! Ocean        :  calendar 
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   day        : calendar
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE phycst          ! physical constants
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE prtctl          ! Print control
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC day        ! called by step.F90
(   21) 
(   22)    !! * Shared module variables
(   23)    INTEGER , PUBLIC ::   &  !:
(   24)       nyear     ,   &  !: current year
(   25)       nmonth    ,   &  !: current month
(   26)       nday      ,   &  !: current day of the month
(   27)       nday_year ,   &  !: curent day counted from jan 1st of the current year
(   28)       ndastp           !: time step date in year/month/day aammjj
(   29)    REAL(wp), PUBLIC ::   &  !:
(   30)        model_time, &   !AD: decimal version of ndastp date for output files
(   31)        adatrj   ,   &  !: number of elapsed days since the begining of the run
(   32)        adatrj0         !: value of adatrj at nit000-1 (before the present run).
(   33)        !               !  it is the accumulated duration of previous runs
(   34)        !               !  that may have been run with different time steps.
(   35)    !!----------------------------------------------------------------------
(   36)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/daymod.F90,v 1.5 2005/09/02 15:45:19 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE day( kt )
(   44)       !!----------------------------------------------------------------------
(   45)       !!                      ***  ROUTINE day  ***
(   46)       !! 
(   47)       !! ** Purpose :   Compute the date with a day iteration IF necessary.
(   48)       !!
(   49)       !! ** Method  : - ???
(   50)       !!
(   51)       !! ** Action  : - nyear     : current year
(   52)       !!              - nmonth    : current month of the year nyear
(   53)       !!              - nday      : current day of the month nmonth
(   54)       !!              - nday_year : current day of the year nyear
(   55)       !!              - ndastp    : =nyear*10000+nmonth*100+nday
(   56)       !!              - adatrj    : date in days since the beginning of the run
(   57)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 6

(   58)       !! History :
(   59)       !!        !  94-09  (M. Pontaud M. Imbard)  Original code
(   60)       !!        !  97-03  (O. Marti)
(   61)       !!        !  97-05  (G. Madec) 
(   62)       !!        !  97-08  (M. Imbard)
(   63)       !!   9.0  !  03-09  (G. Madec)  F90 + nyear, nmonth, nday
(   64)       !!        !  04-01  (A.M. Treguier) new calculation based on adatrj
(   65)       !!----------------------------------------------------------------------      
(   66)       !! * Arguments
(   67)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step indices
(   68) 
(   69)       !! * Local declarations
(   70)       INTEGER  ::   js                   ! dummy loop indice
(   71)       INTEGER  ::   iend, iday0, iday1   ! temporary integers
(   72)       REAL(wp) :: zadatrjn, zadatrjb     ! adatrj at timestep kt-1 and kt-2 
(   73)       CHARACTER (len=25) :: charout
(   74)       !!----------------------------------------------------------------------
(   75) 
(   76)       ! 0.  initialization of adatrj0 and nday, nmonth,nyear, nday_year.
(   77)       !     ndastp has been initialized in domain.F90 or restart.F90
(   78)       !-----------------------------------------------------------------
(   79) 
(   80)       IF( kt == nit000 ) THEN
(   81) 
(   82)          IF( .NOT.ln_rstart )   adatrj0 = 0.e0      ! adatrj0 initialized in rst_read when restart 
(   83) 
(   84)          adatrj  = adatrj0
(   85)          nyear   =   ndastp / 10000
(   86)          nmonth  = ( ndastp - (nyear * 10000) ) / 100
(   87)          nday    =   ndastp - (nyear * 10000) - ( nmonth * 100 ) 
(   88) 
(   89)          ! Calculates nday_year, day since january 1st (useful to read  daily forcing fields)
(   90)          nday_year =  nday
(   91)          !                               ! accumulates days of previous months of this year
(   92)          DO js = 1, nmonth-1
(   93)             IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
(   94)                nday_year = nday_year + nbiss(js)
(   95)             ELSE
(   96)                nday_year = nday_year + nobis(js)
(   97)             ENDIF
(   98)          END DO
(   99) 
(  100)       ENDIF
(  101) 
(  102)       ! I.  calculates adatrj, zadatrjn, zadatrjb.
(  103)       ! ------------------------------------------------------------------
(  104) 
(  105)       adatrj    = adatrj0 + ( kt - nit000 + 1 ) * rdttra(1) / rday
(  106)       zadatrjn  = adatrj0 + ( kt - nit000     ) * rdttra(1) / rday
(  107)       zadatrjb  = adatrj0 + ( kt - nit000 - 1 ) * rdttra(1) / rday
(  108) 
(  109) 
(  110)       ! II.  increment the date.  The date corresponds to 'now' variables (kt-1),
(  111)       !      which is the time step of forcing fields. 
(  112)       !      Do not do this at nit000  unless nrstdt= 2
(  113)       !      In that case ndastp (read in restart) was for step nit000-2
(  114)       ! -------------------------------------------------------------------
(  115) 






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 7

(  116)       iday0 = INT( zadatrjb )
(  117)       iday1 = INT( zadatrjn )
(  118) 
(  119)       IF( iday1 - iday0 >= 1 .AND. ( kt /= nit000 .OR. nrstdt == 2 ) ) THEN
(  120) 
(  121)          ! increase calendar
(  122)          nyear  =   ndastp / 10000
(  123)          nmonth = ( ndastp - (nyear * 10000) ) / 100
(  124)          nday   =   ndastp - (nyear * 10000) - ( nmonth * 100 ) 
(  125)          nday = nday + 1
(  126)          IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
(  127)             iend = nbiss(nmonth)
(  128)          ELSEIF( nleapy > 1 ) THEN 
(  129)             iend = nleapy
(  130)          ELSE 
(  131)             iend = nobis(nmonth)
(  132)          ENDIF
(  133)          IF( nday == iend + 1 ) THEN
(  134)             nday  = 1
(  135)             nmonth = nmonth + 1
(  136)             IF( nmonth == 13 ) THEN
(  137)                nmonth  = 1
(  138)                nyear = nyear + 1
(  139)             ENDIF
(  140)          ENDIF
(  141)          ndastp = nyear * 10000 + nmonth * 100 + nday
(  142) 
(  143)          ! Calculates nday_year, day since january 1st (useful to read  daily forcing fields)
(  144)          nday_year =  nday
(  145)          !                                ! accumulates days of previous months of this year
(  146)          DO js = 1, nmonth-1
(  147)             IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
(  148)                nday_year = nday_year + nbiss(js)
(  149)             ELSE
(  150)                nday_year = nday_year + nobis(js)
(  151)             ENDIF
(  152)          END DO
(  153) 
(  154)          IF(lwp) WRITE(numout,*)' ==============>> time-step =', kt, ' New day, DATE= ',   &
(  155)             &                   nyear, '/', nmonth, '/', nday, 'nday_year:', nday_year
(  156)       ENDIF
(  157) 
(  158)       !AD: NOTE: implemented to have decimal date in output file. Decimal day is relative to 
(  159)       !          ndate0( which is an integer date only  i.e starts at 00 ZZ). So for example
(  160)       !          20051225.5 is noon of december 25, 2005. In future a decimal date could be
(  161)       !          added to ndate0 (i.e convert it to a real). In this case the decimal date
(  162)       !          Here will have to modified to account for the offset (29 Sept 2009)
(  163) 
(  164)       !    NOTE: The date (ndastp) here is for the now model time (kt-1)*rdt should 
(  165)       !           techincally be for end of of step time (kt)*rdt. 
(  166)       !           This will give a adate of 20010101.99 for
(  167)       !           the first output (nwrite=180,rdt=480,nit000=1) at kt=180.     
(  168)       !           To avoid this we add at the time of one time step
(  169) 
(  170)       model_time=FLOAT(ndastp)+(zadatrjn-FLOAT(INT(zadatrjn))) + rdttra(1) / rday 
(  171) 
(  172)       IF(ln_ctl) THEN
(  173)          WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 8

(  174)          CALL prt_ctl_info(charout)
(  175)       ENDIF
(  176) 
(  177)    END SUBROUTINE day
(  178) 
(  179)    !!======================================================================
(  180) END MODULE daymod
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (daymod.F90: 170)
