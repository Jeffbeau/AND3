


PGF90 (Version     10.5)          09/23/2010  09:49:49      page 1

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

Filename: flowri.F90

(    1) !!DB 2008.05.16
(    2) !!Changes:
(    3) !!Output float trajectory data to float_traject.dat, written by lwp
(    4) !!At every output time: Output a file that could be used as a restart file
(    5) 
(    6) MODULE flowri
(    7)    !!======================================================================
(    8)    !!                       ***  MODULE  flowri  ***
(    9)    !! 
(   10)    !!======================================================================
(   11) #if   defined key_floats   ||   defined key_esopa
(   12)    !!----------------------------------------------------------------------
(   13)    !!   'key_floats'                                     float trajectories
(   14)    !!----------------------------------------------------------------------
(   15)    !!    flowri     : write trajectories of floats in file 
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE flo_oce         ! ocean drifting floats
(   19)    USE oce             ! ocean dynamics and tracers
(   20)    USE dom_oce         ! ocean space and time domain
(   21)    USE lib_mpp         ! distribued memory computing library
(   22)    USE daymod
(   23)    USE in_out_manager  ! I/O manager
(   24) 
(   25)    IMPLICIT NONE
(   26) 
(   27)    !! * Accessibility
(   28)    PRIVATE
(   29)    PUBLIC flo_wri     ! routine called by floats.F90
(   30) 
(   31)    !! * Module variables
(   32)       INTEGER :: jfl              ! number of floats
(   33) 
(   34)    !! * Substitutions
(   35) #  include "domzgr_substitute.h90"
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/flowri.F90,v 1.3 2005/03/27 18:35:05 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE flo_wri( kt )
(   45)       !!---------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE flo_wri  ***






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 2

(   47)       !!             
(   48)       !! ** Purpose :   Write position of floats in "trajec_float" file
(   49)       !!      and the temperature and salinity at this position
(   50)       !!      
(   51)       !! ** Method  :   The frequency is nwritefl
(   52)       !!      
(   53)       !!  History :
(   54)       !!    8.0  !  99-09  (Y. Drillet)  Original code
(   55)       !!         !  00-06  (J.-M. Molines)  Profiling floats for CLS 
(   56)       !!    8.5  !  02-10  (A. Bozec)  F90: Free form and module
(   57)       !!----------------------------------------------------------------------
(   58)       !! * Arguments
(   59)       INTEGER  :: kt                               ! time step
(   60) 
(   61)       !! * Local declarations
(   62)       CHARACTER (len=25) ::  clname
(   63)       INTEGER ::   inum = 11       ! temporary logical unit for restart file
(   64)       INTEGER  ::   &
(   65)          iafl,ibfl,icfl,ia1fl,ib1fl,ic1fl,jfl,irecflo,   &
(   66)          iafloc,ibfloc,ia1floc,ib1floc,   &
(   67)          iafln, ibfln
(   68)       INTEGER  ::    ic, jc , jpn
(   69)       INTEGER, DIMENSION ( jpnij )  :: iproc
(   70) 
(   71)       REAL(wp) :: zafl,zbfl,zcfl,zdtj
(   72)       REAL(wp) :: zxxu, zxxu_01,zxxu_10, zxxu_11
(   73) !!DB
(   74)       REAL(wp), DIMENSION (:,:),ALLOCATABLE, SAVE :: ztemp, zsal
(   75)       REAL(wp), DIMENSION (:),ALLOCATABLE, SAVE :: fltemp, flsal, flu, flv, flw
(   76) 
(   77) 
(   78) !!---------------------------------------------------------------------
(   79)       
(   80) !!DB
(   81)       IF( kt == nit000 ) THEN
(   82) !!allocated but not used
(   83)          ALLOCATE(ztemp(jpk,jpnfl));ALLOCATE(zsal(jpk,jpnfl))
(   84) !!float T,S,U,V -- NB do not bother with w at this time
(   85)          ALLOCATE(fltemp(jpnfl)); ALLOCATE(flsal(jpnfl)); ALLOCATE(flu(jpnfl)); ALLOCATE(flv(jpnfl))
(   86)          clname='float_traject.dat'         
(   87)          OPEN (numflo,FILE=clname)
(   88)          if(lwp) then
(   89)             open(555,file='float_traject.info')
(   90)             write(555,*)'Run name  #-floats nwritefl: ', cexper,jpnfl,nwritefl
(   91)             write(555,*)'Columns in ',clname
(   92)             write(555,*)'float-# kt i j k lon lat z niso ngrp T S U V'
(   93)             close(555)
(   94)          endif
(   95)       ENDIF
(   96) 
(   97)       IF( kt == nit000 .OR. MOD( kt,nwritefl)== 0 ) THEN 
(   98) 
(   99)          ! header of output floats file
(  100)       
(  101)          IF(lwp) THEN
(  102)             WRITE(numout,*)
(  103)             WRITE(numout,*) 'flo_wri : write in trajec_float file '
(  104)             WRITE(numout,*) '~~~~~~~    '






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 3

(  105)          ENDIF
(  106) 
(  107) 
(  108)          zdtj = rdt / 86400.      !!bug   use of 86400 instead of the phycst parameter
(  109) 
(  110)          ! translation of index position in geographical position
(  111) 
(  112)          IF( lk_mpp ) THEN
(  113)             DO jfl = 1, jpnfl
(  114)                iafl  = INT ( tpifl(jfl) )
(  115)                ibfl  = INT ( tpjfl(jfl) )
(  116)                icfl  = INT ( tpkfl(jfl) )
(  117)                iafln = NINT( tpifl(jfl) )
(  118)                ibfln = NINT( tpjfl(jfl) )
(  119)                ia1fl = iafl + 1
(  120)                ib1fl = ibfl + 1
(  121)                ic1fl = icfl + 1
(  122)                zafl  = tpifl(jfl) - FLOAT( iafl )
(  123)                zbfl  = tpjfl(jfl) - FLOAT( ibfl )
(  124)                zcfl  = tpkfl(jfl) - FLOAT( icfl )
(  125)                IF(   iafl >= mig(nldi)-jpizoom+1 .AND. iafl <= mig(nlei)-jpizoom+1 .AND.   &
(  126)                   &  ibfl >= mjg(nldj)-jpjzoom+1 .AND. ibfl <= mjg(nlej)-jpjzoom+1       ) THEN
(  127) 
(  128)                   ! local index
(  129) 
(  130)                   iafloc  = iafl -(mig(1)-jpizoom+1) + 1
(  131)                   ibfloc  = ibfl -(mjg(1)-jpjzoom+1) + 1
(  132)                   ia1floc = iafloc + 1
(  133)                   ib1floc = ibfloc + 1
(  134) 
(  135)                   flyy(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
(  136)                      &      +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)
(  137)                   flxx(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
(  138)                      &      +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
(  139)                   flzz(jfl) = (1.-zcfl)*fsdepw(iafloc,ibfloc,icfl ) + zcfl * fsdepw(iafloc,ibfloc,ic1fl)
(  140) 
(  141)                   ! Change  by Alexandra Bozec et Jean-Philippe Boulanger
(  142)                   ! We save  the instantaneous profile of T and S of the column     
(  143)                   ! ztemp(jfl)=tn(iafloc,ibfloc,icfl)
(  144)                   ! zsal(jfl)=sn(iafloc,ibfloc,icfl)
(  145)                   ztemp(1:jpk,jfl) = tn(iafloc,ibfloc,1:jpk)
(  146)                   zsal (1:jpk,jfl) = sn(iafloc,ibfloc,1:jpk)            
(  147) !!DB: float T,S,U,V 
(  148) !!Note that I do not interpolate these to the actual float position
(  149) !!I just use the (i,j,k) values as done just above for ztemp, zsal
(  150)                   fltemp(jfl) =  tn(iafloc,ibfloc,icfl)
(  151)                   flsal(jfl) =  sn(iafloc,ibfloc,icfl)
(  152)                   flu(jfl) =  un(iafloc,ibfloc,icfl)
(  153)                   flv(jfl) =  vn(iafloc,ibfloc,icfl)
(  154)                ELSE
(  155)                   flxx(jfl) = 0.
(  156)                   flyy(jfl) = 0.
(  157)                   flzz(jfl) = 0.
(  158)                   fltemp(jfl) = 0.
(  159)                   flsal(jfl) = 0.
(  160)                   flu(jfl) = 0.
(  161)                   flv(jfl) = 0.
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 4

(  163)                   ztemp(1:jpk,jfl) = 0.
(  164)                   zsal (1:jpk,jfl) = 0.
(  165)                ENDIF
(  166)             END DO
(  167) 
(  168)             CALL mpp_sum( flxx, jpnfl )   ! sums over the global domain
(  169)             CALL mpp_sum( flyy, jpnfl )
(  170)             CALL mpp_sum( flzz, jpnfl )
(  171) !!DB
(  172)             CALL mpp_sum( fltemp, jpnfl )
(  173)             CALL mpp_sum( flsal, jpnfl )
(  174)             CALL mpp_sum( flu, jpnfl )
(  175)             CALL mpp_sum( flv, jpnfl )
(  176) 
(  177) !!DB 2008.03.10 -- does not compile 
(  178) !Error = There is no matching specific subroutine for this generic subroutine call.
(  179) !            CALL mpp_sum( ztemp, jpk*jpnfl )
(  180) !            CALL mpp_sum( zsal , jpk*jpnfl )
(  181) 
(  182) 
(  183) !DB
(  184)             if(lwp) then
(  185)                do jfl = 1, jpnfl
(  186)                   iafl  = INT ( tpifl(jfl) )
(  187)                   ibfl  = INT ( tpjfl(jfl) )
(  188)                   icfl  = INT ( tpkfl(jfl) )
(  189)                   write(numflo,'(5(i6,2x),3(f12.4,1x),2(i3,1x),4(f12.4,1x))')   &
(  190)                        jfl,kt,iafl,ibfl,icfl,flxx(jfl),flyy(jfl),flzz(jfl),nisobfl(jfl),ngrpfl(jfl), &
(  191)                        fltemp(jfl),flsal(jfl),flu(jfl),flv(jfl)
(  192)                enddo
(  193) !!Note that this file is in the input file format (init_float), for ready use
(  194)                open(555,file='restart_float.dat')
(  195)                do jfl = 1, jpnfl
(  196)                   write(555,'(3(f12.4,1x),2(i3,1x))')   &
(  197)                        flxx(jfl),flyy(jfl),flzz(jfl),nisobfl(jfl),ngrpfl(jfl)
(  198)                enddo
(  199)                close(555)
(  200)             endif
(  201) 
(  202)          ELSE
(  203)             DO jfl = 1, jpnfl
(  204)                iafl  = INT (tpifl(jfl))
(  205)                ibfl  = INT (tpjfl(jfl))
(  206)                icfl  = INT (tpkfl(jfl))
(  207)                iafln = NINT(tpifl(jfl))
(  208)                ibfln = NINT(tpjfl(jfl))
(  209)                ia1fl = iafl+1
(  210)                ib1fl = ibfl+1
(  211)                ic1fl = icfl+1
(  212)                zafl  = tpifl(jfl) - FLOAT(iafl)
(  213)                zbfl  = tpjfl(jfl) - FLOAT(ibfl)
(  214)                zcfl  = tpkfl(jfl) - FLOAT(icfl)
(  215)                iafloc  = iafl
(  216)                ibfloc  = ibfl
(  217)                ia1floc = iafloc + 1
(  218)                ib1floc = ibfloc + 1
(  219)                !
(  220)                flyy(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 5

(  221)                          +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)
(  222)                flxx(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
(  223)                          +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
(  224)                flzz(jfl) = (1.-zcfl)*fsdepw(iafloc,ibfloc,icfl ) + zcfl * fsdepw(iafloc,ibfloc,ic1fl)
(  225)                !ALEX
(  226)                ! Astuce pour ne pas avoir des flotteurs qui se baladent sur IDL
(  227)                zxxu_11 = glamt(iafloc ,ibfloc )
(  228)                zxxu_10 = glamt(iafloc ,ib1floc)
(  229)                zxxu_01 = glamt(ia1floc,ibfloc )
(  230)                zxxu    = glamt(ia1floc,ib1floc)
(  231) 
(  232)                IF( iafloc == 52 )  zxxu_10 = -181
(  233)                IF( iafloc == 52 )  zxxu_11 = -181
(  234)                flxx(jfl)=(1.-zafl)*(1.-zbfl)* zxxu_11 + (1.-zafl)*    zbfl * zxxu_10   &
(  235)                         +    zafl *(1.-zbfl)* zxxu_01 +     zafl *    zbfl * zxxu
(  236)                !ALEX         
(  237)                ! Change  by Alexandra Bozec et Jean-Philippe Boulanger
(  238)                ! We save  the instantaneous profile of T and S of the column     
(  239)                !     ztemp(jfl)=tn(iafloc,ibfloc,icfl)
(  240)                !     zsal(jfl)=sn(iafloc,ibfloc,icfl)
(  241)                ztemp(1:jpk,jfl) = tn(iafloc,ibfloc,1:jpk)
(  242)                zsal (1:jpk,jfl) = sn(iafloc,ibfloc,1:jpk)
(  243)             END DO
(  244) !!DBG -- no not cover single CPU case for now
(  245)             write(numflo,*)'DBG: 1 CPU ----> no float output'
(  246) 
(  247) !            WRITE(numflo) flxx,flyy,flzz,nisobfl,ngrpfl,ztemp,zsal, FLOAT(ndastp)
(  248) 
(  249) 
(  250)          ENDIF
(  251) 
(  252)          !
(  253) 
(  254) 
(  255)       !!
(  256)       !! case when profiles are dumped. In order to save memory, dumps are
(  257)       !! done level by level.
(  258)       !      IF (mod(kt,nflclean) == 0.) THEN
(  259)       !!     IF ( nwflo == nwprofil ) THEN
(  260)       !        DO jk = 1,jpk
(  261)       !         DO jfl=1,jpnfl
(  262)       !         iafl= INT(tpifl(jfl))
(  263)       !         ibfl=INT(tpjfl(jfl))
(  264)       !         iafln=NINT(tpifl(jfl))
(  265)       !         ibfln=NINT(tpjfl(jfl))
(  266)       !# if defined key_mpp_mpi   ||   defined key_mpp_shmem
(  267)       !        IF ( (iafl >= (mig(nldi)-jpizoom+1)) .AND.
(  268)       !     $       (iafl <= (mig(nlei)-jpizoom+1)) .AND.
(  269)       !     $       (ibfl >= (mjg(nldj)-jpjzoom+1)) .AND.
(  270)       !     $       (ibfl <= (mjg(nlej)-jpjzoom+1)) ) THEN
(  271)       !!
(  272)       !! local index
(  273)       !!
(  274)       !         iafloc=iafln-(mig(1)-jpizoom+1)+1
(  275)       !         ibfloc=ibfln-(mjg(1)-jpjzoom+1)+1
(  276)       !!         IF (jk == 1 ) THEN
(  277)       !!      PRINT *,'<<<>>> ',jfl,narea, iafloc ,ibfloc, iafln, ibfln,adatrj
(  278)       !!         ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 6

(  279)       !# else
(  280)       !         iafloc=iafln
(  281)       !         ibfloc=ibfln
(  282)       !# endif
(  283)       !         ztemp(jfl)=tn(iafloc,ibfloc,jk)
(  284)       !         zsal(jfl)=sn(iaflo!,ibfloc,jk)
(  285)       !# if defined key_mpp_mpi   ||   defined key_mpp_shmem
(  286)       !        ELSE
(  287)       !         ztemp(jfl) = 0.
(  288)       !         zsal(jfl) = 0.
(  289)       !        ENDIF
(  290)       !# endif
(  291)       !! ... next float
(  292)       !        END DO
(  293)       !      IF( lk_mpp )   CALL mpp_sum( ztemp, jpnfl )
(  294)       !      IF( lk_mpp )   CALL mpp_sum( zsal , jpnfl )
(  295)       !
(  296)       !      IF (lwp) THEN 
(  297)       !         WRITE(numflo) ztemp, zsal
(  298)       !      ENDIF
(  299)       !! ... next level jk
(  300)       !      END DO
(  301)       !! ... reset nwflo to 0 for ALL processors, if profile has been written
(  302)       !!       nwflo = 0
(  303)       !      ENDIF
(  304)       !!
(  305)       !      CALL flush (numflo)
(  306)       !! ... time of dumping floats
(  307)       !!      END IF
(  308)       ENDIF
(  309)       
(  310)       IF( (MOD(kt,nstockfl) == 0) .OR. ( kt == nitend ) ) THEN 
(  311)          ! Writing the restart file 
(  312)          IF(lwp) THEN
(  313)             WRITE(numout,*)
(  314)             WRITE(numout,*) 'flo_wri : write in  restart_float file '
(  315)             WRITE(numout,*) '~~~~~~~    '
(  316)          ENDIF
(  317) 
(  318)          ! file is opened and closed every time it is used.
(  319) 
(  320)          clname = 'restart.float.'
(  321)          ic = 1
(  322)          DO jc = 1, 16
(  323)             IF( cexper(jc:jc) /= ' ' ) ic = jc
(  324)          END DO
(  325)          clname = clname(1:14)//cexper(1:ic)
(  326)          ic = 1
(  327)          DO jc = 1, 48
(  328)             IF( clname(jc:jc) /= ' ' ) ic = jc
(  329)          END DO
(  330) 
(  331)          OPEN (inum,FILE=clname,FORM='UNFORMATTED')
(  332)          REWIND inum
(  333)          !
(  334)          DO jpn = 1, jpnij
(  335)             iproc(jpn) = 0
(  336)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 7

(  337)          !
(  338)          IF(lwp) THEN
(  339)             REWIND(inum)
(  340)             WRITE (inum) tpifl,tpjfl,tpkfl,nisobfl,ngrpfl
(  341)             CLOSE (inum) 
(  342)          ENDIF
(  343)          !
(  344)          ! Compute the number of trajectories for each processor
(  345)          !
(  346)          IF( lk_mpp ) THEN
(  347)             DO jfl = 1, jpnfl
(  348)                IF( (INT(tpifl(jfl)) >= (mig(nldi)-jpizoom+1)) .AND.   &
(  349)                   &(INT(tpifl(jfl)) <= (mig(nlei)-jpizoom+1)) .AND.   &
(  350)                   &(INT(tpjfl(jfl)) >= (mjg(nldj)-jpjzoom+1)) .AND.   &
(  351)                   &(INT(tpjfl(jfl)) <= (mjg(nlej)-jpjzoom+1)) ) THEN
(  352)                   iproc(narea) = iproc(narea)+1
(  353)                ENDIF
(  354)             END DO
(  355)             CALL mpp_sum( iproc, jpnij )
(  356)             !
(  357)             IF(lwp) THEN 
(  358)                WRITE(numout,*) 'DATE',adatrj
(  359)                DO jpn = 1, jpnij
(  360)                   IF( iproc(jpn) /= 0 ) THEN
(  361)                      WRITE(numout,*)'PROCESSOR',jpn-1,'compute',iproc(jpn), 'trajectories.'
(  362)                   ENDIF
(  363)                END DO
(  364)             ENDIF
(  365)          ENDIF
(  366)       ENDIF 
(  367) 
(  368)    END SUBROUTINE flo_wri
(  369) 
(  370) #  else
(  371)    !!----------------------------------------------------------------------
(  372)    !!   Default option                                         Empty module
(  373)    !!----------------------------------------------------------------------
(  374) CONTAINS
(  375)    SUBROUTINE flo_wri                 ! Empty routine
(  376)    END SUBROUTINE flo_wri
(  377) #endif
(  378)    
(  379)    !!======================================================================
(  380) END MODULE flowri




















PGF90 (Version     10.5)          09/23/2010  09:49:49      page 8

(    1) # 1 "flowri.F90"
(    1) !!DB 2008.05.16
(    2) !!Changes:
(    3) !!Output float trajectory data to float_traject.dat, written by lwp
(    4) !!At every output time: Output a file that could be used as a restart file
(    5) 
(    6) MODULE flowri
(    7)    !!======================================================================
(    8)    !!                       ***  MODULE  flowri  ***
(    9)    !! 
(   10)    !!======================================================================
(   11) # 371
(  371)    !!----------------------------------------------------------------------
(  372)    !!   Default option                                         Empty module
(  373)    !!----------------------------------------------------------------------
(  374) CONTAINS
(  375)    SUBROUTINE flo_wri                 ! Empty routine
(  376)    END SUBROUTINE flo_wri
(  377) # 378
(  378)    
(  379)    !!======================================================================
(  380) END MODULE flowri
