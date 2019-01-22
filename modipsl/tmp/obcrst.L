


PGF90 (Version     10.5)          09/23/2010  09:50:07      page 1

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

Filename: obcrst.F90

(    1) MODULE obcrst
(    2) #if defined key_obc
(    3)    !!=================================================================================
(    4)    !!                       ***  MODULE  obcrst  ***
(    5)    !! Ocean dynamic :  Input/Output files for restart on OBC
(    6)    !!=================================================================================
(    7) 
(    8)    !!---------------------------------------------------------------------------------
(    9)    !! * Modules used
(   10)    USE oce             ! ocean dynamics and tracers variables
(   11)    USE dom_oce         ! ocean space and time domain variables
(   12)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   13)    USE phycst          ! physical constants
(   14)    USE obc_oce         ! ocean open boundary conditions
(   15)    USE daymod          ! calendar
(   16)    USE lib_mpp         ! for mppobc
(   17)    USE in_out_manager  ! I/O manager
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Accessibility
(   23)    PUBLIC obc_rst_lec        ! routine called by iniobc.F90
(   24)    PUBLIC obc_rst_wri        ! routine called by step.F90
(   25) 
(   26)    !!---------------------------------------------------------------------------------
(   27)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcrst.F90,v 1.5 2005/12/28 09:25:08 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!---------------------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE obc_rst_wri ( kt )
(   35)       !!--------------------------------------------------------------------------------
(   36)       !!                  ***  SUBROUTINE obc_rst_wri  ***
(   37)       !!                
(   38)       !! ** Purpose :   Write open boundary restart fields in restart.obc.output file 
(   39)       !!
(   40)       !! ** Method  :   restart.obc.output file: Direct access non formatted file.
(   41)       !!      Each nstock time step , save fields which are necessary for restart.
(   42)       !!      - This routine is called if at least the key_obc is defined. It is called
(   43)       !!        at the same time step than rstwri.
(   44)       !!      - First record holds OBC parameters nbobc,jpieob,jpiwob,jpjnob,jpjsob and 
(   45)       !!        the OBC layout jpjed, jpjef ... for checking purposes.
(   46)       !!      - Following records hold the boundary arrays, in the order east west north






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 2

(   47)       !!        south, if they exist.
(   48)       !!      - The writing is realised by vertical slab across the boundary, for bsf, u,
(   49)       !!        v, t, and s boundary arrays. Each record hold a vertical slab.
(   50)       !!      - For mpp, this allows each processor to write only the correct informations
(   51)       !!        it hold. If a processor has no valid informations on boundary, it just 
(   52)       !!        skip the writing part (automatically).
(   53)       !!      - Special care is taken for dumping the starting point of a boundary (jpjed,
(   54)       !!        jpjwd, jpind, jpisd) because of the general definition of nje0 njw0,nin0,
(   55)       !!        nis0. This is done to avoid records to be written by 2 adjacent processors.
(   56)       !!
(   57)       !!  History :
(   58)       !!         ! 97-11 (J.M. Molines) Original code
(   59)       !!         ! 98-11 (J.M. Molines) Bug fix for adjacent processors
(   60)       !!   8.5   ! 02-10 (C. Talandier, A-M. Treguier) F90
(   61)       !!         ! 03-06 (J.M. Molines) Bug fix for adjacent processors
(   62)       !!   9.0   ! 04-02 (G. Madec)  suppression of numwob, use inum
(   63)       !!-----------------------------------------------------------------------------------
(   64)       !! * Arguments
(   65)       INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
(   66) 
(   67)       !! * Local declarations
(   68)       INTEGER ::   ji, jj, jk, ios
(   69)       INTEGER ::   inum = 11          ! temporary logical unit
(   70)       INTEGER ::   ibloc, nreclo, jrec, jt, jb 
(   71)       INTEGER ::   jfoe, jfow, ifon, ifos
(   72)       INTEGER ::   ino0, it0
(   73)       !!-----------------------------------------------------------------------------
(   74) 
(   75)       ! 1. Output of restart fields (inum)
(   76)       ! ------------------------------------
(   77)  
(   78)       IF( ( mod(kt,nstock) == 0 ) .OR. ( kt == nitend ) ) THEN
(   79) 
(   80)          ! 1.0 Initializations
(   81)          ! -------------------
(   82)          IF(lwp) THEN
(   83)               WRITE(numout,*) ' '
(   84)               WRITE(numout,*) 'obcrst: OBC output for restart with obc_rst_wri routine'
(   85)               WRITE(numout,*) '~~~~~~'
(   86)               WRITE(numout,*) '        output done in restart.obc.output file at it= ', kt, ' date= ', ndastp
(   87)          END IF
(   88) 
(   89)          ino0 = no
(   90)          it0  = kt
(   91)          ibloc  = 4096*4
(   92)          nreclo = ibloc*( ( ( 26 *jpk + 9 )*jpbyt -1)/ibloc + 1)
(   93)          IF(lwp) WRITE(numout,*) '             '
(   94)          IF(lwp) WRITE(numout,*) '        OBC restart file opened with nreclo = ',nreclo
(   95) 
(   96)          ! 1.1 Open file
(   97)          ! -------------
(   98)          OPEN( UNIT   = inum,                 &
(   99)                IOSTAT = ios,                  &
(  100)                FILE   = 'restart.obc.output', &
(  101)                ACCESS = 'DIRECT',             &
(  102)                RECL   =  nreclo,              &
(  103)                FORM   = 'UNFORMATTED' )
(  104)          IF( ios > 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 3

(  105)             IF(lwp) WRITE(numout,*) '      '
(  106)             IF(lwp) WRITE(numout,*) '        Pbm to OPEN the restart.obc.output file ' 
(  107)             IF(lwp) WRITE(numout,*) '      '
(  108)             nstop = nstop + 1   
(  109)          END IF
(  110) 
(  111)          ! 1.2 Write header
(  112)          ! ----------------
(  113)          WRITE (inum,REC=1) ino0,it0,nbobc,jpieob,jpiwob,jpjnob,jpjsob,     &
(  114)                               jpjed,jpjef,jpjwd,jpjwf,jpind,jpinf,jpisd,jpisf
(  115) 
(  116)          ! 1.3 Write east boundary array if any.
(  117)          ! -------------------------------------
(  118)          IF( lp_obc_east ) THEN
(  119)             IF( lfbceast ) THEN
(  120)                IF(lwp) THEN
(  121)                   WRITE(numout,*) ' '
(  122)                   WRITE(numout,*) '        No restart file for the fixed east OBC'
(  123)                END IF
(  124)             ELSE
(  125)                IF( jpieob /= 0 ) THEN
(  126)                   IF( nje0+njmpp-1  == jpjed .AND. nie1 >= nie0 ) THEN
(  127)             ! ... dump of jpjed if it is on this proc.
(  128)                      jrec = 2
(  129)                      jfoe = jpjed - njmpp + 1
(  130)                      PRINT *,'Narea =',narea,' write jrec =2 east'
(  131)                      WRITE(inum,REC=jrec)                                    &
(  132) # if defined key_dynspg_rl
(  133)                            ((  bebnd(jfoe,   jb,jt),          jb=1,3),jt=1,3), &
(  134) # endif
(  135)                            ((( uebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  136)                            ((( vebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  137)                            ((( tebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  138)                            ((( sebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  139)                   ENDIF
(  140)                   DO ji = nie0, nie1
(  141)                      DO jj = nje0, nje1
(  142)             ! ... only interested processors go through the following lines
(  143)             !           jfoe = jj + njmpp -1
(  144)                         jfoe = jj 
(  145)                         jrec = 2 + jj + njmpp -1 -jpjed
(  146)                         WRITE (inum,REC=jrec)                                   &
(  147) # if defined key_dynspg_rl
(  148)                               ((  bebnd(jfoe,   jb,jt),          jb=1,3),jt=1,3), &
(  149) # endif
(  150)                               ((( uebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  151)                               ((( vebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  152)                               ((( tebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  153)                               ((( sebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  154)                      END DO
(  155)                   END DO
(  156)                END IF
(  157)             END IF
(  158)          END IF
(  159)  
(  160)          ! 1.4 Write west boundary arrays if any
(  161)          ! -------------------------------------
(  162)          IF( lp_obc_west ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 4

(  163)             IF( lfbcwest ) THEN
(  164)                IF(lwp) THEN
(  165)                   WRITE(numout,*) ' '
(  166)                   WRITE(numout,*) '        No restart file for the fixed west OBC'
(  167)                END IF
(  168)             ELSE
(  169)                IF( jpiwob /= 0 ) THEN
(  170)                   IF( njw0+njmpp+1 == jpjwd .AND. niw1 >= niw0 ) THEN
(  171)             ! ... dump of jpjwd if it is on this proc.
(  172)                      jrec = 3 + jpjef - jpjed
(  173)             !        jfow = jpjwd
(  174)                      jfow = jpjwd -njmpp + 1
(  175)                      PRINT *,'Narea =',narea,' write jrec =',jrec,' west'
(  176)                      WRITE (inum,REC=jrec)                                   &
(  177) # if defined key_dynspg_rl
(  178)                            ((  bwbnd(jfow,   jb,jt),          jb=1,3),jt=1,3), &
(  179) # endif
(  180)                            ((( uwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  181)                            ((( vwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  182)                            ((( twbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  183)                            ((( swbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  184)                   END IF
(  185)                   DO ji = niw0, niw1
(  186)                      DO jj = njw0, njw1
(  187)             ! ... only interested processors go through the following lines
(  188)             !           jfow = jj + njmpp -1
(  189)                         jfow = jj 
(  190)                         jrec = 3 + jpjef -jpjed + jj + njmpp -1 -jpjwd
(  191)                         WRITE (inum,REC=jrec)                                   &
(  192) # if defined key_dynspg_rl
(  193)                               ((  bwbnd(jfow,   jb,jt),          jb=1,3),jt=1,3), &
(  194) # endif
(  195)                               ((( uwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  196)                               ((( vwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  197)                               ((( twbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  198)                               ((( swbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  199)                      END DO
(  200)                   END DO
(  201)                END IF
(  202)             END IF
(  203)          END IF
(  204)  
(  205)          ! 1.5 Write north boundary arrays if any
(  206)          ! --------------------------------------
(  207)          IF( lp_obc_north ) THEN
(  208)             IF( lfbcnorth ) THEN
(  209)                IF(lwp) THEN
(  210)                   WRITE(numout,*) ' '
(  211)                   WRITE(numout,*) '        No restart file for the fixed north OBC'
(  212)                END IF
(  213)             ELSE
(  214)                IF( jpjnob /= 0) THEN
(  215)                   IF( nin0+nimpp-1 == jpind .AND. njn1 >= njn0 ) THEN
(  216)             ! ... dump of jpind if it is on this proc.
(  217)                      jrec = 4 + jpjef -jpjed + jpjwf -jpjwd
(  218)             !        ifon = jpind
(  219)                      ifon = jpind -nimpp +1
(  220)                      WRITE (inum,REC=jrec)                                   &






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 5

(  221) # if defined key_dynspg_rl
(  222)                            ((  bnbnd(ifon,   jb,jt),          jb=1,3),jt=1,3), &
(  223) # endif
(  224)                            ((( unbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  225)                            ((( vnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  226)                            ((( tnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  227)                            ((( snbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  228)                   END IF
(  229)                   DO jj = njn0, njn1
(  230)                      DO ji = nin0, nin1
(  231)             ! ... only interested processors go through the following lines
(  232)             !           ifon = ji + nimpp -1
(  233)                         ifon = ji 
(  234)                         jrec = 4 + jpjef -jpjed + jpjwf -jpjwd +ji + nimpp -1  -jpind
(  235)                         WRITE (inum,REC=jrec)                                   &
(  236) # if defined key_dynspg_rl
(  237)                               ((  bnbnd(ifon,   jb,jt),          jb=1,3),jt=1,3), &
(  238) # endif
(  239)                               ((( unbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  240)                               ((( vnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  241)                               ((( tnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  242)                               ((( snbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  243)                      END DO
(  244)                   END DO
(  245)                END IF
(  246)             END IF
(  247)          END IF
(  248)  
(  249)          ! 1.6 Write south boundary arrays if any
(  250)          ! --------------------------------------
(  251)          IF( lp_obc_south ) THEN
(  252)             IF( lfbcsouth ) THEN
(  253)                IF(lwp) THEN
(  254)                   WRITE(numout,*) ' '
(  255)                   WRITE(numout,*) '        No restart file for the fixed south OBC'
(  256)                   WRITE(numout,*) ' '
(  257)                END IF
(  258)             ELSE
(  259)                IF( jpjsob /= 0 ) THEN
(  260)                   IF( nis0+nimpp-1 == jpisd .AND. njs1 >= njs0 ) THEN
(  261)             ! ... dump of jpisd if it is on this proc.
(  262)                      jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind
(  263)             !        ifos = jpisd
(  264)                      ifos = jpisd -nimpp + 1
(  265)                      WRITE (inum,REC=jrec)                                   &
(  266) # if defined key_dynspg_rl
(  267)                            ((  bsbnd(ifos,   jb,jt),          jb=1,3),jt=1,3), &
(  268) # endif
(  269)                            ((( usbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  270)                            ((( vsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  271)                            ((( tsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  272)                            ((( ssbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  273)                   END IF
(  274)                   DO jj = njs0, njs1
(  275)                      DO ji = nis0, nis1
(  276)             ! ... only interested processors go through the following lines
(  277)             !           ifos = ji + nimpp -1
(  278)                         ifos = ji 






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 6

(  279)                         jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind + &
(  280)                               ji + nimpp -1 -jpisd
(  281)                         WRITE (inum,REC=jrec) &
(  282) # if defined key_dynspg_rl
(  283)                               ((  bsbnd(ifos,   jb,jt),          jb=1,3),jt=1,3), &
(  284) # endif 
(  285)                               ((( usbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  286)                               ((( vsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  287)                               ((( tsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  288)                               ((( ssbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  289)                      END DO
(  290)                   END DO
(  291)                END IF
(  292)             END IF
(  293)          END IF
(  294)       END IF
(  295)       CLOSE(inum)
(  296) 
(  297)    END SUBROUTINE obc_rst_wri
(  298) 
(  299) 
(  300)    SUBROUTINE obc_rst_lec
(  301)       !!----------------------------------------------------------------------------
(  302)       !!                   ***  SUBROUTINE obc_rst_lec  ***
(  303)       !!                   
(  304)       !! ** Purpose :   Read files for restart at open boundaries
(  305)       !!
(  306)       !! ** Method  :   Read the previous boundary arrays on unit inum
(  307)       !!      The first record indicates previous characterics
(  308)       !!
(  309)       !! History :
(  310)       !!        ! 97-11 (J.M. Molines) Original code
(  311)       !!   8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(  312)       !!----------------------------------------------------------------------------
(  313)       !! * Local declarations
(  314)       INTEGER ::   inum = 11            ! temporary logical unit
(  315)       INTEGER ::   ji,jj,jk,ios
(  316)       INTEGER ::   ino0,it0,nbobc0,jpieob0,jpiwob0,jpjnob0,jpjsob0
(  317)       INTEGER ::   ino1,it1,nbobc1,jpieob1,jpiwob1,jpjnob1,jpjsob1
(  318)       INTEGER ::   ied0,ief0,iwd0,iwf0,ind0,inf0,isd0,isf0
(  319)       INTEGER ::   ied1,ief1,iwd1,iwf1,ind1,inf1,isd1,isf1
(  320)       INTEGER ::   ibloc, nreclo, jrec, jt, jb
(  321)       INTEGER ::   jfoe, jfow, ifon, ifos
(  322)       !!-----------------------------------------------------------------------------
(  323) 
(  324)       ! 0. Initialisations
(  325)       ! ------------------
(  326)  
(  327)       ino0    = no
(  328)       it0     = nit000
(  329)       nbobc0  = nbobc
(  330)       jpieob0 = jpieob
(  331)       jpiwob0 = jpiwob
(  332)       jpjnob0 = jpjnob
(  333)       jpjsob0 = jpjsob
(  334)  
(  335)       ied0   = jpjed
(  336)       ief0   = jpjef






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 7

(  337)       iwd0   = jpjwd
(  338)       iwf0   = jpjwf
(  339)       ind0   = jpind
(  340)       inf0   = jpinf
(  341)       isd0   = jpisd
(  342)       isf0   = jpisf
(  343)  
(  344)       ibloc  = 4096*4
(  345)       nreclo = ibloc *( ( ( 26 *jpk + 9 )*jpbyt -1)/ibloc + 1)
(  346)  
(  347)       IF(lwp) THEN
(  348)          WRITE(numout,*) 'obcrst: beginning of restart with obc_rst_lec routine'
(  349)          WRITE(numout,*) '~~~~~~'
(  350)          WRITE(numout,*) ' '
(  351)          WRITE(numout,*) '        The present run :'
(  352)          WRITE(numout,*) '        number job is  : ',no 
(  353)          WRITE(numout,*) '        with the time nit000 : ',nit000
(  354)          WRITE(numout,*) '        OBC restart file opened with nreclo = ',nreclo 
(  355)       END IF
(  356)  
(  357)       ! 0.1 Open files
(  358)       ! ---------------
(  359)       OPEN( UNIT   =  inum,       &
(  360)             IOSTAT =  ios,          &
(  361)             FILE   = 'restart.obc', &
(  362)             ACCESS = 'DIRECT',      &
(  363)             RECL   =  nreclo,       &
(  364)             FORM   = 'UNFORMATTED' )
(  365)       IF( ios > 0 ) THEN
(  366)          IF(lwp) WRITE(numout,*) '        Pbm to OPEN the restart.obc file ' 
(  367)          nstop = nstop + 1   
(  368)       END IF
(  369) 
(  370)       ! 1. Read
(  371)       ! -------
(  372)  
(  373)       ! 1.1 First record
(  374)       ! -----------------
(  375)       READ(inum,REC=1) ino1,it1,nbobc1,jpieob1,jpiwob1,jpjnob1,     &
(  376)                          jpjsob1,ied1,ief1,iwd1,iwf1,ind1,inf1,isd1,isf1
(  377)  
(  378)       IF(lwp) THEN
(  379)          WRITE(numout,*) ' '
(  380)          WRITE(numout,*) '        READ inum with number job : ',ino1,' with the time it: ',it1
(  381)          WRITE(numout,*) ' '
(  382)       END IF
(  383)  
(  384)       ! 1.2 Control of date
(  385)       ! --------------------
(  386)       IF( ( it0-it1 ) /= 1 .AND. abs(nrstdt) == 1 ) THEN
(  387)           IF(lwp) THEN
(  388)              WRITE(numout,*) '        ===>>>> : problem with nit000 for the restart'
(  389)              WRITE(numout,*) '        =============='
(  390)              WRITE(numout,*) '        we stop in obc_rst_lec routine. Verify the file or rerun with the value'
(  391)              WRITE(numout,*) '        0 for the control of time parameter nrstdt'
(  392)              WRITE(numout,*) ' ' 
(  393)           END IF
(  394)           nstop = nstop + 1






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 8

(  395)       END IF
(  396)  
(  397)       ! 1.3 Control of number of open boundaries
(  398)       ! ----------------------------------------
(  399)       IF( nbobc1 /= nbobc0 ) THEN
(  400)          IF(lwp) THEN
(  401)             WRITE(numout,*) '        ===> W A R N I N G: The number of OBC have changed:'
(  402)             WRITE(numout,*) '        Last run : ',nbobc0,' obcs'
(  403)             WRITE(numout,*) '        This run : ',nbobc1,' obcs'
(  404)          END IF
(  405)       END IF
(  406)  
(  407)       ! 1.4 Control of which boundary is open
(  408)       ! -------------------------------------
(  409)       IF( lp_obc_east .AND. ( jpieob1 /= 0 ) ) THEN
(  410)          IF(lwp) THEN
(  411)             WRITE(numout,*) '         '
(  412)             WRITE(numout,*) '        East open boundary'
(  413)             IF( jpieob0 /= jpieob1 ) THEN
(  414)                WRITE(numout,*) '         ==>>>> : Problem in obc_rst_lec, jpieob have changed'
(  415)                nstop = nstop + 1
(  416)             END IF
(  417)          END IF
(  418)       END IF
(  419)  
(  420)       IF( lp_obc_west .AND. ( jpiwob1 /= 0 ) ) THEN
(  421)          IF(lwp) THEN
(  422)             WRITE(numout,*) '         '
(  423)             WRITE(numout,*) '        West open boundary'
(  424)             IF( jpiwob0 /= jpiwob1 ) THEN
(  425)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpiwob has changed'
(  426)                nstop = nstop + 1
(  427)             END IF
(  428)          END IF
(  429)       END IF
(  430)  
(  431)       IF( lp_obc_north .AND. ( jpjnob1 /= 0 ) ) THEN
(  432)          IF(lwp) THEN
(  433)             WRITE(numout,*) '         '
(  434)             WRITE(numout,*) '        North open boundary'
(  435)             IF( jpjnob0 /= jpjnob1 ) THEN
(  436)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjnob has changed'
(  437)                nstop = nstop + 1
(  438)             END IF
(  439)          END IF
(  440)       END IF
(  441)  
(  442)       IF( lp_obc_south .AND. ( jpjsob1 /= 0 ) ) THEN
(  443)          IF(lwp) THEN
(  444)             WRITE(numout,*) '         '
(  445)             WRITE(numout,*) '        South open boundary'
(  446)             IF( jpjsob0 /= jpjsob1) THEN
(  447)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjsob has changed'
(  448)                nstop = nstop + 1
(  449)             END IF
(  450)          END IF
(  451)       END IF
(  452)  






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 9

(  453)  
(  454)       ! 1.5 Control of the limit of the boundaries
(  455)       ! ------------------------------------------
(  456)       IF( lp_obc_east .AND. ( jpieob1 /= 0 ) ) THEN
(  457)          IF(lwp) THEN
(  458)             IF( ied1 /= ied0 ) THEN
(  459)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjed has changed'
(  460)                nstop = nstop + 1
(  461)             END IF
(  462)             IF( ief1 /= ief0 ) THEN
(  463)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjef has changed'
(  464)                nstop = nstop + 1
(  465)             END IF
(  466)          END IF
(  467)       END IF
(  468) 
(  469)       IF( lp_obc_west .AND. ( jpiwob1 /= 0 ) ) THEN
(  470)          IF(lwp) THEN
(  471)             IF( iwd1 /= iwd0 ) THEN
(  472)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjwd has changed'
(  473)                nstop = nstop + 1
(  474)             END IF
(  475)             IF( iwf1 /= iwf0 ) THEN
(  476)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjwf has changed'
(  477)                nstop = nstop + 1
(  478)             END IF
(  479)          END IF
(  480)       END IF
(  481)  
(  482)       IF( lp_obc_north .AND. ( jpjnob1 /= 0 ) ) THEN
(  483)          IF(lwp) THEN
(  484)             IF( ind1 /= ind0 ) THEN
(  485)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpind has changed'
(  486)                nstop = nstop + 1
(  487)             END IF
(  488)             IF( inf1 /= inf0 ) THEN
(  489)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpinf has changed'
(  490)                nstop = nstop + 1
(  491)              END IF
(  492)           END IF
(  493)       END IF
(  494)  
(  495)       IF( lp_obc_south .AND. ( jpjsob1 /= 0 ) ) THEN
(  496)          IF(lwp) THEN
(  497)             IF( isd1 /= isd0 ) THEN
(  498)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpisd has changed'
(  499)                nstop = nstop + 1
(  500)             END IF
(  501)             IF( isf1 /= isf0 ) THEN
(  502)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpisf has changed'
(  503)                nstop = nstop + 1
(  504)             END IF
(  505)          END IF
(  506)       END IF
(  507)  
(  508)  
(  509)       ! 2. Now read the boundary arrays
(  510)       ! -------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 10

(  511)  
(  512)       ! 2.1 Read east boundary array if any.
(  513)       ! ------------------------------------
(  514)       IF( lp_obc_east ) THEN
(  515)          IF( jpieob1 /= 0) THEN
(  516)             IF( nje0+njmpp-1 == jpjed .AND. nie1 >= nie0 ) THEN
(  517)       ! ... read of jpjed if it is on this proc.
(  518)                jrec = 2
(  519)       !        jfoe = jpjed
(  520)                jfoe = jpjed -njmpp + 1
(  521)                READ (inum,REC=jrec)                                   &
(  522) # if defined key_dynspg_rl
(  523)                     ((  bebnd(jfoe,   jb,jt),          jb=1,3),jt=1,3), &
(  524) # endif 
(  525)                     ((( uebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  526)                     ((( vebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  527)                     ((( tebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  528)                     ((( sebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  529)             END IF
(  530)             DO ji = nie0, nie1
(  531)                DO jj = nje0, nje1
(  532)       ! ... only interested processors go through the following lines
(  533)       !           jfoe = jj + njmpp -1
(  534)                   jfoe = jj 
(  535)                   jrec = 2 + jj + njmpp -1 -jpjed
(  536)                   READ (inum,REC=jrec)                                   &
(  537) # if defined key_dynspg_rl
(  538)                        ((  bebnd(jfoe,   jb,jt),          jb=1,3),jt=1,3), &
(  539) # endif
(  540)                        ((( uebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  541)                        ((( vebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  542)                        ((( tebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  543)                        ((( sebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  544)                END DO
(  545)             END DO
(  546) 
(  547)          ELSE
(  548) 
(  549)             !  lp_obc_east was not TRUE previously
(  550)          END IF
(  551) 
(  552)       END IF
(  553)  
(  554)       ! 2.2 Read west boundary arrays if any.
(  555)       ! -------------------------------------
(  556)       IF( lp_obc_west ) THEN
(  557)          IF( jpiwob1 /= 0) THEN
(  558)             IF( njw0+njmpp-1 == jpjwd .AND. niw1 >= niw0 ) THEN
(  559)       ! ... read of jpjwd if it is on this proc.
(  560)                jrec = 3 + jpjef - jpjed
(  561)       !        jfow = jpjwd
(  562)                jfow = jpjwd -njmpp + 1
(  563)                READ (inum,REC=jrec)                                   &
(  564) # if defined key_dynspg_rl
(  565)                     ((  bwbnd(jfow,   jb,jt),          jb=1,3),jt=1,3), &
(  566) # endif
(  567)                     ((( uwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  568)                     ((( vwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 11

(  569)                     ((( twbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  570)                     ((( swbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  571)             END IF
(  572)             DO ji = niw0, niw1
(  573)                DO jj = njw0, njw1
(  574)       ! ... only interested processors go through the following lines
(  575)       !           jfow = jj + njmpp -1
(  576)                   jfow = jj 
(  577)                   jrec = 3 + jpjef -jpjed + jj + njmpp -1 -jpjwd
(  578)                   READ (inum,REC=jrec)                                   &
(  579) # if defined key_dynspg_rl
(  580)                        ((  bwbnd(jfow,   jb,jt),          jb=1,3),jt=1,3), &
(  581) # endif
(  582)                        ((( uwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  583)                        ((( vwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  584)                        ((( twbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  585)                        ((( swbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  586)                END DO
(  587)             END DO
(  588) 
(  589)          ELSE
(  590) 
(  591)             !  lp_obc_west was not TRUE previously
(  592)          END IF
(  593) 
(  594)       END IF
(  595)  
(  596)       ! 2.3 Read north boundary arrays if any.
(  597)       ! --------------------------------------
(  598)       IF( lp_obc_north ) THEN
(  599)          IF( jpjnob1 /= 0) THEN
(  600)             IF( nin0+nimpp-1 == jpind .AND. njn1 >= njn0 ) THEN
(  601)       ! ... read of jpind if it is on this proc.
(  602)                jrec = 4 + jpjef -jpjed + jpjwf -jpjwd
(  603)       !        ifon = jpind
(  604)                ifon = jpind -nimpp +1
(  605)                READ (inum,REC=jrec)                                   &
(  606) # if defined key_dynspg_rl
(  607)                     ((  bnbnd(ifon,   jb,jt),          jb=1,3),jt=1,3), &
(  608) # endif
(  609)                     ((( unbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  610)                     ((( vnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  611)                     ((( tnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), & 
(  612)                     ((( snbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  613)             END IF
(  614)             DO jj = njn0, njn1
(  615)                DO ji = nin0, nin1
(  616)       ! ... only interested processors go through the following lines
(  617)       !           ifon = ji + nimpp -1
(  618)                   ifon = ji 
(  619)                   jrec = 4 + jpjef -jpjed + jpjwf -jpjwd +ji + nimpp -1  -jpind
(  620)                   READ (inum,REC=jrec)                                   & 
(  621) # if defined key_dynspg_rl
(  622)                        ((  bnbnd(ifon,   jb,jt),          jb=1,3),jt=1,3), &
(  623) # endif  
(  624)                        ((( unbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  625)                        ((( vnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  626)                        ((( tnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 12

(  627)                        ((( snbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  628)                END DO
(  629)             END DO
(  630) 
(  631)          ELSE
(  632) 
(  633)            !  lp_obc_north was not TRUE previously
(  634)          END IF
(  635) 
(  636)       END IF
(  637)  
(  638)       ! 2.4 Read south boundary arrays if any.
(  639)       ! -------------------------------------
(  640)       IF( lp_obc_south ) THEN
(  641)          IF( jpjsob1 /= 0) THEN
(  642)             IF( nis0+nimpp-1 == jpisd .AND. njs1 >= njs0 ) THEN
(  643)       ! ... read of jpisd if it is on this proc.
(  644)                jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind
(  645)       !        ifos = jpisd
(  646)                ifos = jpisd -nimpp + 1
(  647)                READ (inum,REC=jrec)                                   &
(  648) # if defined key_dynspg_rl
(  649)                     ((  bsbnd(ifos,   jb,jt),          jb=1,3),jt=1,3), &
(  650) # endif
(  651)                     ((( usbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  652)                     ((( vsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  653)                     ((( tsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  654)                     ((( ssbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  655)             END IF
(  656)             DO jj = njs0, njs1
(  657)                DO ji = nis0, nis1
(  658)       ! ... only interested processors go through the following lines
(  659)       !           ifos = ji + nimpp -1
(  660)                   ifos = ji 
(  661)                   jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind +  &
(  662)                         ji + nimpp -1 -jpisd
(  663)                   READ (inum,REC=jrec)                                   & 
(  664) # if defined key_dynspg_rl
(  665)                        ((  bsbnd(ifos,   jb,jt),          jb=1,3),jt=1,3), &
(  666) # endif
(  667)                        ((( usbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  668)                        ((( vsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  669)                        ((( tsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  670)                        ((( ssbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  671)                END DO
(  672)             END DO
(  673)          ELSE
(  674)             !  lp_obc_south was not TRUE previously
(  675)          END IF
(  676) 
(  677)       END IF
(  678)       CLOSE(inum)
(  679) 
(  680)       IF( lk_mpp ) THEN
(  681)          IF( lp_obc_east ) THEN
(  682) # if defined key_dynspg_rl
(  683)             CALL mppobc(bebnd,jpjed,jpjef,jpieob,3*3,2,jpj)
(  684) # endif






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 13

(  685)             CALL mppobc(uebnd,jpjed,jpjef,jpieob,jpk*3*3,2,jpj)
(  686)             CALL mppobc(vebnd,jpjed,jpjef,jpieob+1,jpk*3*3,2,jpj)
(  687)             CALL mppobc(tebnd,jpjed,jpjef,jpieob+1,jpk*2*2,2,jpj)
(  688)             CALL mppobc(sebnd,jpjed,jpjef,jpieob+1,jpk*2*2,2,jpj)
(  689)          ENDIF
(  690)          IF( lp_obc_west ) THEN
(  691) # if defined key_dynspg_rl
(  692)             CALL mppobc(bwbnd,jpjwd,jpjwf,jpiwob,3*3,2,jpj)
(  693) # endif
(  694)             CALL mppobc(uwbnd,jpjwd,jpjwf,jpiwob,jpk*3*3,2,jpj)
(  695)             CALL mppobc(vwbnd,jpjwd,jpjwf,jpiwob,jpk*3*3,2,jpj)
(  696)             CALL mppobc(twbnd,jpjwd,jpjwf,jpiwob,jpk*2*2,2,jpj)
(  697)             CALL mppobc(swbnd,jpjwd,jpjwf,jpiwob,jpk*2*2,2,jpj)
(  698)          ENDIF
(  699)          IF( lp_obc_north ) THEN 
(  700) # if defined key_dynspg_rl
(  701)             CALL mppobc(bnbnd,jpind,jpinf,jpjnob  ,3*3    ,1,jpi)
(  702) # endif
(  703)             CALL mppobc(unbnd,jpind,jpinf,jpjnob+1,jpk*3*3,1,jpi)
(  704)             CALL mppobc(vnbnd,jpind,jpinf,jpjnob  ,jpk*3*3,1,jpi)
(  705)             CALL mppobc(tnbnd,jpind,jpinf,jpjnob+1,jpk*2*2,1,jpi)
(  706)             CALL mppobc(snbnd,jpind,jpinf,jpjnob+1,jpk*2*2,1,jpi)
(  707)          ENDIF
(  708)          IF( lp_obc_south ) THEN
(  709) # if defined key_dynspg_rl
(  710)             CALL mppobc(bsbnd,jpisd,jpisf,jpjsob,    3*3,1,jpi)
(  711) # endif
(  712)             CALL mppobc(usbnd,jpisd,jpisf,jpjsob,jpk*3*3,1,jpi)
(  713)             CALL mppobc(vsbnd,jpisd,jpisf,jpjsob,jpk*3*3,1,jpi)
(  714)             CALL mppobc(tsbnd,jpisd,jpisf,jpjsob,jpk*2*2,1,jpi)
(  715)             CALL mppobc(ssbnd,jpisd,jpisf,jpjsob,jpk*2*2,1,jpi)
(  716)          ENDIF
(  717)       ENDIF
(  718)  
(  719)    END SUBROUTINE obc_rst_lec
(  720) #else
(  721)    !!=================================================================================
(  722)    !!                       ***  MODULE  obcrst  ***
(  723)    !! Ocean dynamic :  Input/Output files for restart on OBC
(  724)    !!=================================================================================
(  725) CONTAINS
(  726)    SUBROUTINE obc_rst_wri( kt )           !  No Open boundary ==> empty routine
(  727)       INTEGER,INTENT(in) :: kt
(  728) !      WRITE(*,*) 'obc_rst_wri: You should not have seen this print! error?', kt
(  729)    END SUBROUTINE obc_rst_wri
(  730)    SUBROUTINE obc_rst_lec                 !  No Open boundary ==> empty routine
(  731)    END SUBROUTINE obc_rst_lec
(  732) #endif
(  733) 
(  734)    !!=================================================================================
(  735) END MODULE obcrst













PGF90 (Version     10.5)          09/23/2010  09:50:07      page 14

(    1) # 1 "obcrst.F90"
(    1) MODULE obcrst
(    2) # 3
(    3)    !!=================================================================================
(    4)    !!                       ***  MODULE  obcrst  ***
(    5)    !! Ocean dynamic :  Input/Output files for restart on OBC
(    6)    !!=================================================================================
(    7) 
(    8)    !!---------------------------------------------------------------------------------
(    9)    !! * Modules used
(   10)    USE oce             ! ocean dynamics and tracers variables
(   11)    USE dom_oce         ! ocean space and time domain variables
(   12)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   13)    USE phycst          ! physical constants
(   14)    USE obc_oce         ! ocean open boundary conditions
(   15)    USE daymod          ! calendar
(   16)    USE lib_mpp         ! for mppobc
(   17)    USE in_out_manager  ! I/O manager
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Accessibility
(   23)    PUBLIC obc_rst_lec        ! routine called by iniobc.F90
(   24)    PUBLIC obc_rst_wri        ! routine called by step.F90
(   25) 
(   26)    !!---------------------------------------------------------------------------------
(   27)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcrst.F90,v 1.5 2005/12/28 09:25:08 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!---------------------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE obc_rst_wri ( kt )
(   35)       !!--------------------------------------------------------------------------------
(   36)       !!                  ***  SUBROUTINE obc_rst_wri  ***
(   37)       !!                
(   38)       !! ** Purpose :   Write open boundary restart fields in restart.obc.output file 
(   39)       !!
(   40)       !! ** Method  :   restart.obc.output file: Direct access non formatted file.
(   41)       !!      Each nstock time step , save fields which are necessary for restart.
(   42)       !!      - This routine is called if at least the 1 is defined. It is called
(   43)       !!        at the same time step than rstwri.
(   44)       !!      - First record holds OBC parameters nbobc,jpieob,jpiwob,jpjnob,jpjsob and 
(   45)       !!        the OBC layout jpjed, jpjef ... for checking purposes.
(   46)       !!      - Following records hold the boundary arrays, in the order east west north
(   47)       !!        south, if they exist.
(   48)       !!      - The writing is realised by vertical slab across the boundary, for bsf, u,
(   49)       !!        v, t, and s boundary arrays. Each record hold a vertical slab.
(   50)       !!      - For mpp, this allows each processor to write only the correct informations
(   51)       !!        it hold. If a processor has no valid informations on boundary, it just 
(   52)       !!        skip the writing part (automatically).
(   53)       !!      - Special care is taken for dumping the starting point of a boundary (jpjed,
(   54)       !!        jpjwd, jpind, jpisd) because of the general definition of nje0 njw0,nin0,
(   55)       !!        nis0. This is done to avoid records to be written by 2 adjacent processors.
(   56)       !!
(   57)       !!  History :






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 15

(   58)       !!         ! 97-11 (J.M. Molines) Original code
(   59)       !!         ! 98-11 (J.M. Molines) Bug fix for adjacent processors
(   60)       !!   8.5   ! 02-10 (C. Talandier, A-M. Treguier) F90
(   61)       !!         ! 03-06 (J.M. Molines) Bug fix for adjacent processors
(   62)       !!   9.0   ! 04-02 (G. Madec)  suppression of numwob, use inum
(   63)       !!-----------------------------------------------------------------------------------
(   64)       !! * Arguments
(   65)       INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
(   66) 
(   67)       !! * Local declarations
(   68)       INTEGER ::   ji, jj, jk, ios
(   69)       INTEGER ::   inum = 11          ! temporary logical unit
(   70)       INTEGER ::   ibloc, nreclo, jrec, jt, jb 
(   71)       INTEGER ::   jfoe, jfow, ifon, ifos
(   72)       INTEGER ::   ino0, it0
(   73)       !!-----------------------------------------------------------------------------
(   74) 
(   75)       ! 1. Output of restart fields (inum)
(   76)       ! ------------------------------------
(   77)  
(   78)       IF( ( mod(kt,nstock) == 0 ) .OR. ( kt == nitend ) ) THEN
(   79) 
(   80)          ! 1.0 Initializations
(   81)          ! -------------------
(   82)          IF(lwp) THEN
(   83)               WRITE(numout,*) ' '
(   84)               WRITE(numout,*) 'obcrst: OBC output for restart with obc_rst_wri routine'
(   85)               WRITE(numout,*) '~~~~~~'
(   86)               WRITE(numout,*) '        output done in restart.obc.output file at it= ', kt, ' date= ', ndastp
(   87)          END IF
(   88) 
(   89)          ino0 = no
(   90)          it0  = kt
(   91)          ibloc  = 4096*4
(   92)          nreclo = ibloc*( ( ( 26 *jpk + 9 )*jpbyt -1)/ibloc + 1)
(   93)          IF(lwp) WRITE(numout,*) '             '
(   94)          IF(lwp) WRITE(numout,*) '        OBC restart file opened with nreclo = ',nreclo
(   95) 
(   96)          ! 1.1 Open file
(   97)          ! -------------
(   98)          OPEN( UNIT   = inum,                 &
(   99)                IOSTAT = ios,                  &
(  100)                FILE   = 'restart.obc.output', &
(  101)                ACCESS = 'DIRECT',             &
(  102)                RECL   =  nreclo,              &
(  103)                FORM   = 'UNFORMATTED' )
(  104)          IF( ios > 0 ) THEN
(  105)             IF(lwp) WRITE(numout,*) '      '
(  106)             IF(lwp) WRITE(numout,*) '        Pbm to OPEN the restart.obc.output file ' 
(  107)             IF(lwp) WRITE(numout,*) '      '
(  108)             nstop = nstop + 1   
(  109)          END IF
(  110) 
(  111)          ! 1.2 Write header
(  112)          ! ----------------
(  113)          WRITE (inum,REC=1) ino0,it0,nbobc,jpieob,jpiwob,jpjnob,jpjsob,     &
(  114)                               jpjed,jpjef,jpjwd,jpjwf,jpind,jpinf,jpisd,jpisf
(  115) 






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 16

(  116)          ! 1.3 Write east boundary array if any.
(  117)          ! -------------------------------------
(  118)          IF( lp_obc_east ) THEN
(  119)             IF( lfbceast ) THEN
(  120)                IF(lwp) THEN
(  121)                   WRITE(numout,*) ' '
(  122)                   WRITE(numout,*) '        No restart file for the fixed east OBC'
(  123)                END IF
(  124)             ELSE
(  125)                IF( jpieob /= 0 ) THEN
(  126)                   IF( nje0+njmpp-1  == jpjed .AND. nie1 >= nie0 ) THEN
(  127)             ! ... dump of jpjed if it is on this proc.
(  128)                      jrec = 2
(  129)                      jfoe = jpjed - njmpp + 1
(  130)                      PRINT *,'Narea =',narea,' write jrec =2 east'
(  131)                      WRITE(inum,REC=jrec)                                    &
(  135)                            ((( uebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  136)                            ((( vebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  137)                            ((( tebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  138)                            ((( sebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  139)                   ENDIF
(  140)                   DO ji = nie0, nie1
(  141)                      DO jj = nje0, nje1
(  142)             ! ... only interested processors go through the following lines
(  143)             !           jfoe = jj + njmpp -1
(  144)                         jfoe = jj 
(  145)                         jrec = 2 + jj + njmpp -1 -jpjed
(  146)                         WRITE (inum,REC=jrec)                                   &
(  150)                               ((( uebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  151)                               ((( vebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  152)                               ((( tebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  153)                               ((( sebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  154)                      END DO
(  155)                   END DO
(  156)                END IF
(  157)             END IF
(  158)          END IF
(  159)  
(  160)          ! 1.4 Write west boundary arrays if any
(  161)          ! -------------------------------------
(  162)          IF( lp_obc_west ) THEN
(  163)             IF( lfbcwest ) THEN
(  164)                IF(lwp) THEN
(  165)                   WRITE(numout,*) ' '
(  166)                   WRITE(numout,*) '        No restart file for the fixed west OBC'
(  167)                END IF
(  168)             ELSE
(  169)                IF( jpiwob /= 0 ) THEN
(  170)                   IF( njw0+njmpp+1 == jpjwd .AND. niw1 >= niw0 ) THEN
(  171)             ! ... dump of jpjwd if it is on this proc.
(  172)                      jrec = 3 + jpjef - jpjed
(  173)             !        jfow = jpjwd
(  174)                      jfow = jpjwd -njmpp + 1
(  175)                      PRINT *,'Narea =',narea,' write jrec =',jrec,' west'
(  176)                      WRITE (inum,REC=jrec)                                   &
(  180)                            ((( uwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  181)                            ((( vwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  182)                            ((( twbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 17

(  183)                            ((( swbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  184)                   END IF
(  185)                   DO ji = niw0, niw1
(  186)                      DO jj = njw0, njw1
(  187)             ! ... only interested processors go through the following lines
(  188)             !           jfow = jj + njmpp -1
(  189)                         jfow = jj 
(  190)                         jrec = 3 + jpjef -jpjed + jj + njmpp -1 -jpjwd
(  191)                         WRITE (inum,REC=jrec)                                   &
(  195)                               ((( uwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  196)                               ((( vwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  197)                               ((( twbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  198)                               ((( swbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  199)                      END DO
(  200)                   END DO
(  201)                END IF
(  202)             END IF
(  203)          END IF
(  204)  
(  205)          ! 1.5 Write north boundary arrays if any
(  206)          ! --------------------------------------
(  207)          IF( lp_obc_north ) THEN
(  208)             IF( lfbcnorth ) THEN
(  209)                IF(lwp) THEN
(  210)                   WRITE(numout,*) ' '
(  211)                   WRITE(numout,*) '        No restart file for the fixed north OBC'
(  212)                END IF
(  213)             ELSE
(  214)                IF( jpjnob /= 0) THEN
(  215)                   IF( nin0+nimpp-1 == jpind .AND. njn1 >= njn0 ) THEN
(  216)             ! ... dump of jpind if it is on this proc.
(  217)                      jrec = 4 + jpjef -jpjed + jpjwf -jpjwd
(  218)             !        ifon = jpind
(  219)                      ifon = jpind -nimpp +1
(  220)                      WRITE (inum,REC=jrec)                                   &
(  224)                            ((( unbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  225)                            ((( vnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  226)                            ((( tnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  227)                            ((( snbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  228)                   END IF
(  229)                   DO jj = njn0, njn1
(  230)                      DO ji = nin0, nin1
(  231)             ! ... only interested processors go through the following lines
(  232)             !           ifon = ji + nimpp -1
(  233)                         ifon = ji 
(  234)                         jrec = 4 + jpjef -jpjed + jpjwf -jpjwd +ji + nimpp -1  -jpind
(  235)                         WRITE (inum,REC=jrec)                                   &
(  239)                               ((( unbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  240)                               ((( vnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  241)                               ((( tnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  242)                               ((( snbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  243)                      END DO
(  244)                   END DO
(  245)                END IF
(  246)             END IF
(  247)          END IF
(  248)  
(  249)          ! 1.6 Write south boundary arrays if any






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 18

(  250)          ! --------------------------------------
(  251)          IF( lp_obc_south ) THEN
(  252)             IF( lfbcsouth ) THEN
(  253)                IF(lwp) THEN
(  254)                   WRITE(numout,*) ' '
(  255)                   WRITE(numout,*) '        No restart file for the fixed south OBC'
(  256)                   WRITE(numout,*) ' '
(  257)                END IF
(  258)             ELSE
(  259)                IF( jpjsob /= 0 ) THEN
(  260)                   IF( nis0+nimpp-1 == jpisd .AND. njs1 >= njs0 ) THEN
(  261)             ! ... dump of jpisd if it is on this proc.
(  262)                      jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind
(  263)             !        ifos = jpisd
(  264)                      ifos = jpisd -nimpp + 1
(  265)                      WRITE (inum,REC=jrec)                                   &
(  269)                            ((( usbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  270)                            ((( vsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  271)                            ((( tsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  272)                            ((( ssbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  273)                   END IF
(  274)                   DO jj = njs0, njs1
(  275)                      DO ji = nis0, nis1
(  276)             ! ... only interested processors go through the following lines
(  277)             !           ifos = ji + nimpp -1
(  278)                         ifos = ji 
(  279)                         jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind + &
(  280)                               ji + nimpp -1 -jpisd
(  281)                         WRITE (inum,REC=jrec) &
(  285)                               ((( usbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  286)                               ((( vsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  287)                               ((( tsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  288)                               ((( ssbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  289)                      END DO
(  290)                   END DO
(  291)                END IF
(  292)             END IF
(  293)          END IF
(  294)       END IF
(  295)       CLOSE(inum)
(  296) 
(  297)    END SUBROUTINE obc_rst_wri
(  298) 
(  299) 
(  300)    SUBROUTINE obc_rst_lec
(  301)       !!----------------------------------------------------------------------------
(  302)       !!                   ***  SUBROUTINE obc_rst_lec  ***
(  303)       !!                   
(  304)       !! ** Purpose :   Read files for restart at open boundaries
(  305)       !!
(  306)       !! ** Method  :   Read the previous boundary arrays on unit inum
(  307)       !!      The first record indicates previous characterics
(  308)       !!
(  309)       !! History :
(  310)       !!        ! 97-11 (J.M. Molines) Original code
(  311)       !!   8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(  312)       !!----------------------------------------------------------------------------
(  313)       !! * Local declarations






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 19

(  314)       INTEGER ::   inum = 11            ! temporary logical unit
(  315)       INTEGER ::   ji,jj,jk,ios
(  316)       INTEGER ::   ino0,it0,nbobc0,jpieob0,jpiwob0,jpjnob0,jpjsob0
(  317)       INTEGER ::   ino1,it1,nbobc1,jpieob1,jpiwob1,jpjnob1,jpjsob1
(  318)       INTEGER ::   ied0,ief0,iwd0,iwf0,ind0,inf0,isd0,isf0
(  319)       INTEGER ::   ied1,ief1,iwd1,iwf1,ind1,inf1,isd1,isf1
(  320)       INTEGER ::   ibloc, nreclo, jrec, jt, jb
(  321)       INTEGER ::   jfoe, jfow, ifon, ifos
(  322)       !!-----------------------------------------------------------------------------
(  323) 
(  324)       ! 0. Initialisations
(  325)       ! ------------------
(  326)  
(  327)       ino0    = no
(  328)       it0     = nit000
(  329)       nbobc0  = nbobc
(  330)       jpieob0 = jpieob
(  331)       jpiwob0 = jpiwob
(  332)       jpjnob0 = jpjnob
(  333)       jpjsob0 = jpjsob
(  334)  
(  335)       ied0   = jpjed
(  336)       ief0   = jpjef
(  337)       iwd0   = jpjwd
(  338)       iwf0   = jpjwf
(  339)       ind0   = jpind
(  340)       inf0   = jpinf
(  341)       isd0   = jpisd
(  342)       isf0   = jpisf
(  343)  
(  344)       ibloc  = 4096*4
(  345)       nreclo = ibloc *( ( ( 26 *jpk + 9 )*jpbyt -1)/ibloc + 1)
(  346)  
(  347)       IF(lwp) THEN
(  348)          WRITE(numout,*) 'obcrst: beginning of restart with obc_rst_lec routine'
(  349)          WRITE(numout,*) '~~~~~~'
(  350)          WRITE(numout,*) ' '
(  351)          WRITE(numout,*) '        The present run :'
(  352)          WRITE(numout,*) '        number job is  : ',no 
(  353)          WRITE(numout,*) '        with the time nit000 : ',nit000
(  354)          WRITE(numout,*) '        OBC restart file opened with nreclo = ',nreclo 
(  355)       END IF
(  356)  
(  357)       ! 0.1 Open files
(  358)       ! ---------------
(  359)       OPEN( UNIT   =  inum,       &
(  360)             IOSTAT =  ios,          &
(  361)             FILE   = 'restart.obc', &
(  362)             ACCESS = 'DIRECT',      &
(  363)             RECL   =  nreclo,       &
(  364)             FORM   = 'UNFORMATTED' )
(  365)       IF( ios > 0 ) THEN
(  366)          IF(lwp) WRITE(numout,*) '        Pbm to OPEN the restart.obc file ' 
(  367)          nstop = nstop + 1   
(  368)       END IF
(  369) 
(  370)       ! 1. Read
(  371)       ! -------






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 20

(  372)  
(  373)       ! 1.1 First record
(  374)       ! -----------------
(  375)       READ(inum,REC=1) ino1,it1,nbobc1,jpieob1,jpiwob1,jpjnob1,     &
(  376)                          jpjsob1,ied1,ief1,iwd1,iwf1,ind1,inf1,isd1,isf1
(  377)  
(  378)       IF(lwp) THEN
(  379)          WRITE(numout,*) ' '
(  380)          WRITE(numout,*) '        READ inum with number job : ',ino1,' with the time it: ',it1
(  381)          WRITE(numout,*) ' '
(  382)       END IF
(  383)  
(  384)       ! 1.2 Control of date
(  385)       ! --------------------
(  386)       IF( ( it0-it1 ) /= 1 .AND. abs(nrstdt) == 1 ) THEN
(  387)           IF(lwp) THEN
(  388)              WRITE(numout,*) '        ===>>>> : problem with nit000 for the restart'
(  389)              WRITE(numout,*) '        =============='
(  390)              WRITE(numout,*) '        we stop in obc_rst_lec routine. Verify the file or rerun with the value'
(  391)              WRITE(numout,*) '        0 for the control of time parameter nrstdt'
(  392)              WRITE(numout,*) ' ' 
(  393)           END IF
(  394)           nstop = nstop + 1
(  395)       END IF
(  396)  
(  397)       ! 1.3 Control of number of open boundaries
(  398)       ! ----------------------------------------
(  399)       IF( nbobc1 /= nbobc0 ) THEN
(  400)          IF(lwp) THEN
(  401)             WRITE(numout,*) '        ===> W A R N I N G: The number of OBC have changed:'
(  402)             WRITE(numout,*) '        Last run : ',nbobc0,' obcs'
(  403)             WRITE(numout,*) '        This run : ',nbobc1,' obcs'
(  404)          END IF
(  405)       END IF
(  406)  
(  407)       ! 1.4 Control of which boundary is open
(  408)       ! -------------------------------------
(  409)       IF( lp_obc_east .AND. ( jpieob1 /= 0 ) ) THEN
(  410)          IF(lwp) THEN
(  411)             WRITE(numout,*) '         '
(  412)             WRITE(numout,*) '        East open boundary'
(  413)             IF( jpieob0 /= jpieob1 ) THEN
(  414)                WRITE(numout,*) '         ==>>>> : Problem in obc_rst_lec, jpieob have changed'
(  415)                nstop = nstop + 1
(  416)             END IF
(  417)          END IF
(  418)       END IF
(  419)  
(  420)       IF( lp_obc_west .AND. ( jpiwob1 /= 0 ) ) THEN
(  421)          IF(lwp) THEN
(  422)             WRITE(numout,*) '         '
(  423)             WRITE(numout,*) '        West open boundary'
(  424)             IF( jpiwob0 /= jpiwob1 ) THEN
(  425)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpiwob has changed'
(  426)                nstop = nstop + 1
(  427)             END IF
(  428)          END IF
(  429)       END IF






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 21

(  430)  
(  431)       IF( lp_obc_north .AND. ( jpjnob1 /= 0 ) ) THEN
(  432)          IF(lwp) THEN
(  433)             WRITE(numout,*) '         '
(  434)             WRITE(numout,*) '        North open boundary'
(  435)             IF( jpjnob0 /= jpjnob1 ) THEN
(  436)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjnob has changed'
(  437)                nstop = nstop + 1
(  438)             END IF
(  439)          END IF
(  440)       END IF
(  441)  
(  442)       IF( lp_obc_south .AND. ( jpjsob1 /= 0 ) ) THEN
(  443)          IF(lwp) THEN
(  444)             WRITE(numout,*) '         '
(  445)             WRITE(numout,*) '        South open boundary'
(  446)             IF( jpjsob0 /= jpjsob1) THEN
(  447)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjsob has changed'
(  448)                nstop = nstop + 1
(  449)             END IF
(  450)          END IF
(  451)       END IF
(  452)  
(  453)  
(  454)       ! 1.5 Control of the limit of the boundaries
(  455)       ! ------------------------------------------
(  456)       IF( lp_obc_east .AND. ( jpieob1 /= 0 ) ) THEN
(  457)          IF(lwp) THEN
(  458)             IF( ied1 /= ied0 ) THEN
(  459)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjed has changed'
(  460)                nstop = nstop + 1
(  461)             END IF
(  462)             IF( ief1 /= ief0 ) THEN
(  463)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjef has changed'
(  464)                nstop = nstop + 1
(  465)             END IF
(  466)          END IF
(  467)       END IF
(  468) 
(  469)       IF( lp_obc_west .AND. ( jpiwob1 /= 0 ) ) THEN
(  470)          IF(lwp) THEN
(  471)             IF( iwd1 /= iwd0 ) THEN
(  472)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjwd has changed'
(  473)                nstop = nstop + 1
(  474)             END IF
(  475)             IF( iwf1 /= iwf0 ) THEN
(  476)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpjwf has changed'
(  477)                nstop = nstop + 1
(  478)             END IF
(  479)          END IF
(  480)       END IF
(  481)  
(  482)       IF( lp_obc_north .AND. ( jpjnob1 /= 0 ) ) THEN
(  483)          IF(lwp) THEN
(  484)             IF( ind1 /= ind0 ) THEN
(  485)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpind has changed'
(  486)                nstop = nstop + 1
(  487)             END IF






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 22

(  488)             IF( inf1 /= inf0 ) THEN
(  489)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpinf has changed'
(  490)                nstop = nstop + 1
(  491)              END IF
(  492)           END IF
(  493)       END IF
(  494)  
(  495)       IF( lp_obc_south .AND. ( jpjsob1 /= 0 ) ) THEN
(  496)          IF(lwp) THEN
(  497)             IF( isd1 /= isd0 ) THEN
(  498)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpisd has changed'
(  499)                nstop = nstop + 1
(  500)             END IF
(  501)             IF( isf1 /= isf0 ) THEN
(  502)                WRITE(numout,*) '        ==>>>> : Problem in obc_rst_lec, jpisf has changed'
(  503)                nstop = nstop + 1
(  504)             END IF
(  505)          END IF
(  506)       END IF
(  507)  
(  508)  
(  509)       ! 2. Now read the boundary arrays
(  510)       ! -------------------------------
(  511)  
(  512)       ! 2.1 Read east boundary array if any.
(  513)       ! ------------------------------------
(  514)       IF( lp_obc_east ) THEN
(  515)          IF( jpieob1 /= 0) THEN
(  516)             IF( nje0+njmpp-1 == jpjed .AND. nie1 >= nie0 ) THEN
(  517)       ! ... read of jpjed if it is on this proc.
(  518)                jrec = 2
(  519)       !        jfoe = jpjed
(  520)                jfoe = jpjed -njmpp + 1
(  521)                READ (inum,REC=jrec)                                   &
(  525)                     ((( uebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  526)                     ((( vebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  527)                     ((( tebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  528)                     ((( sebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  529)             END IF
(  530)             DO ji = nie0, nie1
(  531)                DO jj = nje0, nje1
(  532)       ! ... only interested processors go through the following lines
(  533)       !           jfoe = jj + njmpp -1
(  534)                   jfoe = jj 
(  535)                   jrec = 2 + jj + njmpp -1 -jpjed
(  536)                   READ (inum,REC=jrec)                                   &
(  540)                        ((( uebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  541)                        ((( vebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  542)                        ((( tebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  543)                        ((( sebnd(jfoe,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  544)                END DO
(  545)             END DO
(  546) 
(  547)          ELSE
(  548) 
(  549)             !  lp_obc_east was not TRUE previously
(  550)          END IF
(  551) 






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 23

(  552)       END IF
(  553)  
(  554)       ! 2.2 Read west boundary arrays if any.
(  555)       ! -------------------------------------
(  556)       IF( lp_obc_west ) THEN
(  557)          IF( jpiwob1 /= 0) THEN
(  558)             IF( njw0+njmpp-1 == jpjwd .AND. niw1 >= niw0 ) THEN
(  559)       ! ... read of jpjwd if it is on this proc.
(  560)                jrec = 3 + jpjef - jpjed
(  561)       !        jfow = jpjwd
(  562)                jfow = jpjwd -njmpp + 1
(  563)                READ (inum,REC=jrec)                                   &
(  567)                     ((( uwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  568)                     ((( vwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  569)                     ((( twbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  570)                     ((( swbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  571)             END IF
(  572)             DO ji = niw0, niw1
(  573)                DO jj = njw0, njw1
(  574)       ! ... only interested processors go through the following lines
(  575)       !           jfow = jj + njmpp -1
(  576)                   jfow = jj 
(  577)                   jrec = 3 + jpjef -jpjed + jj + njmpp -1 -jpjwd
(  578)                   READ (inum,REC=jrec)                                   &
(  582)                        ((( uwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  583)                        ((( vwbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  584)                        ((( twbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  585)                        ((( swbnd(jfow,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  586)                END DO
(  587)             END DO
(  588) 
(  589)          ELSE
(  590) 
(  591)             !  lp_obc_west was not TRUE previously
(  592)          END IF
(  593) 
(  594)       END IF
(  595)  
(  596)       ! 2.3 Read north boundary arrays if any.
(  597)       ! --------------------------------------
(  598)       IF( lp_obc_north ) THEN
(  599)          IF( jpjnob1 /= 0) THEN
(  600)             IF( nin0+nimpp-1 == jpind .AND. njn1 >= njn0 ) THEN
(  601)       ! ... read of jpind if it is on this proc.
(  602)                jrec = 4 + jpjef -jpjed + jpjwf -jpjwd
(  603)       !        ifon = jpind
(  604)                ifon = jpind -nimpp +1
(  605)                READ (inum,REC=jrec)                                   &
(  609)                     ((( unbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  610)                     ((( vnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  611)                     ((( tnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), & 
(  612)                     ((( snbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  613)             END IF
(  614)             DO jj = njn0, njn1
(  615)                DO ji = nin0, nin1
(  616)       ! ... only interested processors go through the following lines
(  617)       !           ifon = ji + nimpp -1
(  618)                   ifon = ji 






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 24

(  619)                   jrec = 4 + jpjef -jpjed + jpjwf -jpjwd +ji + nimpp -1  -jpind
(  620)                   READ (inum,REC=jrec)                                   & 
(  624)                        ((( unbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  625)                        ((( vnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  626)                        ((( tnbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  627)                        ((( snbnd(ifon,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  628)                END DO
(  629)             END DO
(  630) 
(  631)          ELSE
(  632) 
(  633)            !  lp_obc_north was not TRUE previously
(  634)          END IF
(  635) 
(  636)       END IF
(  637)  
(  638)       ! 2.4 Read south boundary arrays if any.
(  639)       ! -------------------------------------
(  640)       IF( lp_obc_south ) THEN
(  641)          IF( jpjsob1 /= 0) THEN
(  642)             IF( nis0+nimpp-1 == jpisd .AND. njs1 >= njs0 ) THEN
(  643)       ! ... read of jpisd if it is on this proc.
(  644)                jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind
(  645)       !        ifos = jpisd
(  646)                ifos = jpisd -nimpp + 1
(  647)                READ (inum,REC=jrec)                                   &
(  651)                     ((( usbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  652)                     ((( vsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  653)                     ((( tsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  654)                     ((( ssbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  655)             END IF
(  656)             DO jj = njs0, njs1
(  657)                DO ji = nis0, nis1
(  658)       ! ... only interested processors go through the following lines
(  659)       !           ifos = ji + nimpp -1
(  660)                   ifos = ji 
(  661)                   jrec = 5 + jpjef -jpjed + jpjwf -jpjwd +jpinf -jpind +  &
(  662)                         ji + nimpp -1 -jpisd
(  663)                   READ (inum,REC=jrec)                                   & 
(  667)                        ((( usbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  668)                        ((( vsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,3),jt=1,3), &
(  669)                        ((( tsbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2), &
(  670)                        ((( ssbnd(ifos,jk,jb,jt),jk=1,jpk),jb=1,2),jt=1,2)
(  671)                END DO
(  672)             END DO
(  673)          ELSE
(  674)             !  lp_obc_south was not TRUE previously
(  675)          END IF
(  676) 
(  677)       END IF
(  678)       CLOSE(inum)
(  679) 
(  680)       IF( lk_mpp ) THEN
(  681)          IF( lp_obc_east ) THEN
(  682) # 685
(  685)             CALL mppobc(uebnd,jpjed,jpjef,jpieob,jpk*3*3,2,jpj)
(  686)             CALL mppobc(vebnd,jpjed,jpjef,jpieob+1,jpk*3*3,2,jpj)
(  687)             CALL mppobc(tebnd,jpjed,jpjef,jpieob+1,jpk*2*2,2,jpj)






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 25

(  688)             CALL mppobc(sebnd,jpjed,jpjef,jpieob+1,jpk*2*2,2,jpj)
(  689)          ENDIF
(  690)          IF( lp_obc_west ) THEN
(  691) # 694
(  694)             CALL mppobc(uwbnd,jpjwd,jpjwf,jpiwob,jpk*3*3,2,jpj)
(  695)             CALL mppobc(vwbnd,jpjwd,jpjwf,jpiwob,jpk*3*3,2,jpj)
(  696)             CALL mppobc(twbnd,jpjwd,jpjwf,jpiwob,jpk*2*2,2,jpj)
(  697)             CALL mppobc(swbnd,jpjwd,jpjwf,jpiwob,jpk*2*2,2,jpj)
(  698)          ENDIF
(  699)          IF( lp_obc_north ) THEN 
(  700) # 703
(  703)             CALL mppobc(unbnd,jpind,jpinf,jpjnob+1,jpk*3*3,1,jpi)
(  704)             CALL mppobc(vnbnd,jpind,jpinf,jpjnob  ,jpk*3*3,1,jpi)
(  705)             CALL mppobc(tnbnd,jpind,jpinf,jpjnob+1,jpk*2*2,1,jpi)
(  706)             CALL mppobc(snbnd,jpind,jpinf,jpjnob+1,jpk*2*2,1,jpi)
(  707)          ENDIF
(  708)          IF( lp_obc_south ) THEN
(  709) # 712
(  712)             CALL mppobc(usbnd,jpisd,jpisf,jpjsob,jpk*3*3,1,jpi)
(  713)             CALL mppobc(vsbnd,jpisd,jpisf,jpjsob,jpk*3*3,1,jpi)
(  714)             CALL mppobc(tsbnd,jpisd,jpisf,jpjsob,jpk*2*2,1,jpi)
(  715)             CALL mppobc(ssbnd,jpisd,jpisf,jpjsob,jpk*2*2,1,jpi)
(  716)          ENDIF
(  717)       ENDIF
(  718)  
(  719)    END SUBROUTINE obc_rst_lec
(  720) 
(  721) # 734
(  734)    !!=================================================================================
(  735) END MODULE obcrst
