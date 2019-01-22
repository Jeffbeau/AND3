


PGF90 (Version     10.5)          09/23/2010  09:50:06      page 1

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

Filename: obcdom.F90

(    1) MODULE obcdom
(    2)    !!=================================================================================
(    3)    !!                       ***  MODULE  obcdom  ***
(    4)    !! Space domain  :  get all the isolated coastline points needed to resolve the 
(    5)    !!                  barotropic streamfunction elliptic equation associated with 
(    6)    !!                  the open boundaries.
(    7)    !!=================================================================================
(    8) #if defined key_obc && defined key_dynspg_rl
(    9)    !!---------------------------------------------------------------------------------
(   10)    !!   'key_obc'           AND                                Open Boundary Condition
(   11)    !!   'key_dynspg_rl'                                          Rigid-Lid formulation
(   12)    !!---------------------------------------------------------------------------------
(   13)    !!   obc_dom        : domain initialization in rid-lid formulation
(   14)    !!---------------------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers   
(   17)    USE dom_oce         ! ocean space and time domain 
(   18)    USE phycst          ! physical constants
(   19)    USE obc_oce         ! ocean open boundary conditions
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE lib_mpp         ! distributed memory computing library
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Accessibility
(   27)    PUBLIC obc_dom        ! routine called by iniobc.F90
(   28)    !!---------------------------------------------------------------------------------
(   29) 
(   30) CONTAINS
(   31) 
(   32)    SUBROUTINE obc_dom
(   33)       !!------------------------------------------------------------------------------
(   34)       !!                       SUBROUTINE obc_dom
(   35)       !!                      ********************
(   36)       !! ** Purpose :   Initialize the array used for the computation of the part of 
(   37)       !!        the right hand side of the barotropic streamfunction elliptic equation
(   38)       !!        associated with the open boundaries
(   39)       !!
(   40)       !! **  Method :
(   41)       !!      + The (i,j) indices of ocean grid-points round isolated coastlines
(   42)       !!        are found (isolated coastlines = coast lines separated by an
(   43)       !!        open boundary) from icoast array read in coastlines file.
(   44)       !!
(   45)       !!      + read 'coastline' file  initialize icoast()
(   46)       !!        modify icoast() depending on the number of open boundaries 






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 2

(   47)       !!        specified through key_obc
(   48)       !!
(   49)       !!      + compute zwb, an ocean/land mask defined as follows:
(   50)       !!             zwb(i,j)  =  0. over the one isolated coastline
(   51)       !!                       = -1, -2, -3 over the orthers
(   52)       !!      + for example, when 4 open boundaries are specified:
(   53)       !!
(   54)       !!                       //|           |// 
(   55)       !!        North          //|   North   |// -1 -1   North
(   56)       !!         West     0  0 //| - - - - - |// -1 -1    East
(   57)       !!                       //| open bnd  |//
(   58)       !!              ///////////|           |/////////
(   59)       !!              ------------            ----------
(   60)       !!                                               
(   61)       !!             west   |                  |    east
(   62)       !!           open bnd                       open bnd    
(   63)       !!                    |                  |         
(   64)       !!              ___________             _________ 
(   65)       !!              ///////////|           |/////////
(   66)       !!                       //|   south   |// 
(   67)       !!        South    -3 -3 //| - - - - - |// -2 -2   South
(   68)       !!         West    -3 -3 //| open bnd  |// -2 -2    East
(   69)       !!                       //|           |// 
(   70)       !!
(   71)       !!        With the proper boundary conditions (defined by nperio)
(   72)       !!
(   73)       !!        C a u t i o n :  no check, the user must enter a well defined
(   74)       !!        coastline file. Further more, he must verify that isolated
(   75)       !!        coastlines have been well located dans that the right potential
(   76)       !!        is affected to the right coastline in obc.F
(   77)       !!
(   78)       !! History :
(   79)       !!   8.1  !  09-97  (J.M. Molines, G. Madec)  Original code
(   80)       !!   8.2  !  06-99  (J.M. Molines) suppress zwb(,) for ATL6 (memory saving)
(   81)       !!        !  02-02  (A.M. Treguier) icoast in 2 dimension
(   82)       !!   8.5  !  02-08  (G. Madec)  F90 : free form
(   83)       !!----------------------------------------------------------------------
(   84)       !! * Local declarations
(   85)       INTEGER ::   ji, jj, jn, jnic, jnp, jii    ! dummy loop indices
(   86)       INTEGER ::   inum = 11         ! temporary logical unit
(   87)       INTEGER ::   ifreq, il1, il2, ii, ij, icheck
(   88)       INTEGER ::   ip, ipn, ips, ipe, ipw
(   89)       INTEGER ::   iim, ijm, iii, ijj
(   90)       INTEGER, DIMENSION(jpidta,jpjdta) ::   icoast
(   91)       CHARACTER (len=15) ::   clexp
(   92)       REAL(wp) ::   zzic, zland
(   93)       REAL(wp) ::   zwb, zwbn, zwbs, zwbe, zwbw
(   94)       REAL(wp) ::   zglo(jpiglo,jpjglo)
(   95)       !!---------------------------------------------------------------------
(   96)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   97)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdom.F90,v 1.4 2005/03/27 18:35:10 opalod Exp $ 
(   98)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   99)       !!---------------------------------------------------------------------
(  100)       
(  101)       ! 0. initialization of gcfobc to zero
(  102)       ! -----------------------------------
(  103)       
(  104)       DO jn = 1, 3






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 3

(  105)          gcfobc(:,:,jn) = 0.e0
(  106)       END DO
(  107)       
(  108)       ! 1. Only 1 open boundary : gcfobc is zero, return
(  109)       ! ------------------------------------------------
(  110)       
(  111)       IF( nbobc == 1 .OR. nbic == 0 ) THEN 
(  112)          IF(lwp) WRITE(numout,*)
(  113)          IF(lwp) WRITE(numout,*) ' obc_dom: No isolated coastlines gcfobc is set to zero'
(  114)          IF(lwp) WRITE(numout,*) ' ~~~~~~~'
(  115)          nstop = nstop + 1
(  116)       END IF
(  117) 
(  118)       ! 2. Lecture of 'coastlines' file
(  119)       ! -------------------------------
(  120)       
(  121)       IF(lwp) WRITE(numout,*)
(  122)       IF(lwp) WRITE(numout,*) 'obc_dom: define isolated coastlines from "coastlines" file'
(  123)       IF(lwp) WRITE(numout,*) '~~~~~~~'
(  124)       IF(lwp) WRITE(numout,*)
(  125) 
(  126)       ! open coastlines file'
(  127)       CALL ctlopn( inum, 'coastlines', 'OLD', 'FORMATTED', 'SEQUENTIAL',   &
(  128)                    1 , numout, lwp, 1 )
(  129) 
(  130)       ! lecture of coastlines, set icoast array
(  131)       ! Note that this is coded for jpjdta > 1000
(  132)       REWIND(inum)
(  133)       READ(inum,9101) clexp, iim, ijm
(  134)       READ(inum,'(/)')
(  135)       ifreq = 40
(  136)       il1 = 1
(  137)       IF( jpjglo < 1000 ) THEN
(  138)          DO jn = 1, jpidta/ifreq+1
(  139)             READ(inum,'(/)')
(  140)             il2 = min0( jpidta, il1+ifreq-1 )
(  141)             READ(inum,9201) ( ii, ji = il1, il2, 5 )
(  142)             READ(inum,'(/)')
(  143)             DO jj = jpjdta, 1, -1
(  144)                READ(inum,9202) ij, ( icoast(ji,jj), ji = il1, il2 )
(  145)             END DO
(  146)             il1 = il1 + ifreq
(  147)          END DO
(  148)       ELSE
(  149)          DO jn = 1, jpidta/ifreq+1
(  150)             READ(inum,'(/)')
(  151)             il2 = min0( jpidta, il1+ifreq-1 )
(  152)             READ(inum,9221) ( ii, ji = il1, il2, 5 )
(  153)             READ(inum,'(/)')
(  154)             DO jj = jpjdta, 1, -1
(  155)                READ(inum,9222) ij, ( icoast(ji,jj), ji = il1, il2 )
(  156)             END DO
(  157)             il1 = il1 + ifreq
(  158)          END DO
(  159)       END IF
(  160)       CLOSE(inum)
(  161) 
(  162)    ! in case of zoom, icoast must be set to 0 on the domain border






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 4

(  163)    ! it must be the same for the bathymetry
(  164)    IF (lzoom_w) icoast(jpiglo            ,:) = 0 
(  165)    IF (lzoom_e) icoast(jpiglo +jpizoom -1,:) = 0 
(  166)    IF (lzoom_s) icoast(:,jpjzoom           ) = 0 
(  167)    IF (lzoom_n) icoast(:,jpjglo+jpjzoom -1 ) = 0 
(  168) 
(  169)       DO jj = 1, jpjglo
(  170)          DO ji = 1, jpiglo
(  171)             zglo(ji,jj) = icoast( ji+jpizoom-1, jj+jpjzoom-1)
(  172)          END DO
(  173)       END DO
(  174) 
(  175)  9101 FORMAT(1x,a15,2i8)
(  176)  9201 FORMAT(3x,13(i3,12x))
(  177)  9202 FORMAT(i3,41i3)
(  178)  9221 FORMAT(4x,13(i3,12x))
(  179)  9222 FORMAT(i4,41i3)
(  180) 
(  181)       ! check consistency between tmask and icoast
(  182)       
(  183)       icheck = 0
(  184)       DO jj = 1, jpjm1
(  185)          DO ji = 1, jpim1
(  186)             icheck = icheck +  INT( tmask(ji,jj,1) ) - MAX(  0, icoast( mig(ji), mjg(jj) )  )
(  187)          END DO
(  188)       END DO
(  189)       IF( lk_mpp )   CALL mpp_sum(icheck)   ! sum over the global domain
(  190) 
(  191)       IF( icheck /= 0 ) THEN
(  192)          IF(lwp) WRITE(numout,cform_err)
(  193)          IF(lwp) WRITE(numout,*) 'obc_dom : tmask and isolated coastlines mask are not equal', icheck
(  194)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(  195)          nstop = nstop + 1
(  196)       END IF
(  197) 
(  198)       ! 3. transfer the coastline information from T- to f-points
(  199)       !    (i.e. from icoast to zwb with  zwb=0 over the continent
(  200)       !     and ocean, =-n over the nth isolated coastline)
(  201)       ! -----------------------------------------------------------
(  202)       
(  203)       ! east open boundary
(  204)       IF( lp_obc_east .AND. ( jpieob /= 0 ) ) THEN
(  205)          IF(lwp) WRITE(numout,*) '         East open boundary: from coastline S.E : ', &
(  206)                                  INT(zglo(jpieob,jpjed)),' to N.E : ',               &
(  207)                                  INT(zglo(jpieob,jpjef))
(  208)       END IF
(  209)       ! west open boundary
(  210)       IF( lp_obc_west .AND. ( jpiwob /= 0 ) ) THEN
(  211)          IF(lwp) WRITE(numout,*) '         West open boundary: from coastline S.W : ', &
(  212)                                  INT(zglo(jpiwob,jpjwd)),' to N.W : ',               &
(  213)                                  INT(zglo(jpiwob,jpjwf))
(  214)       END IF
(  215)       ! north open boundary
(  216)       IF( lp_obc_north .AND. ( jpjnob /= 0 ) ) THEN
(  217)          IF(lwp) WRITE(numout,*) '         North open boundary: from coastline N.W : ', &
(  218)                                  INT(zglo(jpind,jpjnob)),' to N.E : ',                &
(  219)                                  INT(zglo(jpinf,jpjnob))
(  220)       END IF






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 5

(  221)       ! south open boundary
(  222)       IF( lp_obc_south .AND. ( jpjsob /= 0 ) ) THEN
(  223)          IF(lwp) WRITE(numout,*) '         South open boundary: from coastline S.W : ', &
(  224)                                  INT(zglo(jpisd,jpjsob)),' to S.E : ',                &
(  225)                                  INT(zglo(jpisf,jpjsob))
(  226)       END IF
(  227) 
(  228)       ! 4. Identify the isolated coastline grid point position
(  229)       ! ------------------------------------------------------
(  230) 
(  231)       ! Loop over isolated coastlines
(  232) 
(  233)       DO jnic = 1, nbobc-1
(  234)          ! set to zero of miic, mjic of the jnic isolated coastline
(  235)          DO jn = 0, 4
(  236)             DO ji = 1, jpnic
(  237)                miic(ji,jn,jnic) = 0
(  238)                mjic(ji,jn,jnic) = 0
(  239)             END DO
(  240)          END DO
(  241)          
(  242)          ! Coastal isolated coastline grid-points (miic,mjic)
(  243)          ip  = 0
(  244)          ipn = 0
(  245)          ips = 0
(  246)          ipe = 0
(  247)          ipw = 0
(  248)          
(  249)          ! Middle lines (1=<jj=<jpjm1)
(  250) 
(  251)          !       jj+1  --zwb--v--ZWB--v--zwb--
(  252)          !                |       |       |
(  253)          !          jj+1  u   T   u   T   u
(  254)          !                |       |       |
(  255)          !       jj    --ZWB--v--ZWB--v--ZWB--
(  256)          !                |       |       |
(  257)          !              jj    u   T   u   T   u
(  258)          !                |       |       |
(  259)          !       jj-1  --zwb--v--ZWB--v--zwb--
(  260)          !                |       |       | 
(  261)          !                |   ii  | ii+1  |
(  262)          !                |       |       | 
(  263)          !               ii-1    ii      ii+1
(  264)          
(  265)          DO jj = 1, jpjglo-1
(  266)             DO ji = 1, jpiglo-1
(  267)                ii = ji
(  268)                zwb = MIN( 0., zglo(ji,jj), zglo(ji+1,jj), zglo(ji,jj+1), zglo(ji+1,jj+1) )
(  269)                IF( jj == jpjglo -1 ) THEN
(  270)                   zwbn = zwb
(  271)                ELSE
(  272)                   zwbn= MIN( 0., zglo(ji,jj+1), zglo(ji+1,jj+1), zglo(ji,jj+2), zglo(ji+1,jj+2) )
(  273)                END IF
(  274)                IF( jj == 1 ) THEN
(  275)                   zwbs = zwb
(  276)                ELSE
(  277)                   zwbs= MIN( 0., zglo(ji,jj-1), zglo(ji+1,jj-1), zglo(ji,jj), zglo(ji+1,jj) )
(  278)                END IF






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 6

(  279)                IF( ji == jpiglo -1 ) THEN
(  280)                   zwbe = zwb
(  281)                ELSE
(  282)                   zwbe= MIN( 0., zglo(ji+1,jj), zglo(ji+2,jj), zglo(ji+1,jj+1), zglo(ji+2,jj+1) )
(  283)                END IF
(  284)                IF( ji == 1 ) THEN
(  285)                   zwbw = zwb
(  286)                ELSE
(  287)                   zwbw= MIN( 0., zglo(ji-1,jj), zglo(ji,jj), zglo(ji-1,jj+1), zglo(ji,jj+1) )
(  288)                END IF
(  289)                
(  290)                ! inside coastlines indicator
(  291)                zzic  =                zwbn                 &
(  292)                      * zwbw                       * zwbe   &
(  293)                                     * zwbs
(  294)                ! inside land indicator
(  295)                zland = MAX( 0., zglo(ji,jj+1) ) + MAX( 0., zglo(ji+1,jj+1) )   &
(  296)                      + MAX( 0., zglo(ji,jj  ) ) + MAX( 0., zglo(ji+1,jj  ) )
(  297)                ! if isolated coastline grid-point 
(  298)                IF( zwb == float( -jnic ) .AND.   &
(  299)                   ! not inside the isolated coastline
(  300)                              zzic == 0.  .AND.   &
(  301)                   ! not inside the land
(  302)                              zland >= 2.         ) THEN
(  303)                   ! coastal point of the isolated coastline jnic
(  304)                   ip = ip + 1
(  305)                   miic(ip,0,jnic) = ii
(  306)                   mjic(ip,0,jnic) = jj
(  307)                   ! which has a west ocean grid point 
(  308)                   IF( zwbw == 0. ) THEN
(  309)                      ipw = ipw + 1
(  310)                      miic(ipw,4,jnic) = ii
(  311)                      mjic(ipw,4,jnic) = jj
(  312)                   END IF
(  313)                   ! which has a east ocean grid point
(  314)                   IF( zwbe == 0. ) THEN
(  315)                      ipe = ipe + 1
(  316)                      IF( nperio == 1 .AND. ii == jpiglo-1 ) THEN
(  317)                         miic(ipe,3,jnic) = 2
(  318)                      ELSE
(  319)                         miic(ipe,3,jnic) = ii + 1
(  320)                      END IF
(  321)                      mjic(ipe,3,jnic) = jj
(  322)                   END IF
(  323)                   ! which has a south ocean grid point 
(  324)                   IF( zwbs == 0. ) THEN
(  325)                      ips = ips + 1
(  326)                      miic(ips,2,jnic) = ii
(  327)                      mjic(ips,2,jnic) = jj
(  328)                   END IF
(  329)                   ! which has a north ocean grid point not out of north open b.
(  330)                   IF( zwbn == 0. ) THEN
(  331)                      ipn = ipn + 1
(  332)                      miic(ipn,1,jnic) = ii
(  333)                      mjic(ipn,1,jnic) = jj + 1
(  334)                   END IF
(  335)                END IF
(  336)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 7

(  337)          END DO
(  338)          
(  339)          mnic(0,jnic) = ip
(  340)          mnic(1,jnic) = ipn
(  341)          mnic(2,jnic) = ips
(  342)          mnic(3,jnic) = ipe
(  343)          mnic(4,jnic) = ipw
(  344)          
(  345)       END DO
(  346) 
(  347)       ! 5. Check the number of isolated coastline
(  348)       ! -----------------------------------------
(  349) 
(  350)       DO jnic = 1, nbobc-1
(  351)          IF( mnic(0,jnic) > jpnic ) THEN
(  352)             IF(lwp) WRITE(numout,cform_err)
(  353)             IF(lwp) WRITE(numout,*) 'obc_dom: isolated coastline ',jnic,   &
(  354)                ' has ',ip,' grid-points > ',jpnic 
(  355)             IF(lwp) WRITE(numout,*) '~~~~~~~'
(  356)             IF(lwp) WRITE(numout,*) ' modify this dimension in obc_dom'
(  357)             nstop = nstop + 1
(  358)          END IF
(  359)          IF( mnic(0,jnic) == 0 ) THEN
(  360)             IF(lwp) WRITE(numout,cform_err)
(  361)             IF(lwp) WRITE(numout,*) 'obc_dom: isolated coastline ',jnic,   &
(  362)                ' has 0  grid-points verify coastlines file'
(  363)             IF(lwp) WRITE(numout,*) '~~~~~~~'
(  364)             nstop = nstop + 1
(  365)          END IF
(  366)       END DO
(  367) 
(  368)       ! 6. Print of isolated coastline parametres and arrays
(  369)       ! -----------------------------------------------------
(  370) 
(  371)       IF(lwp) WRITE(numout,*) '           '
(  372)       IF(lwp) WRITE(numout,*) '         isolated coastlines found:', nbobc - 1
(  373) 
(  374)       DO jnic = 1, nbobc-1
(  375)          ip  = mnic(0,jnic)
(  376)          ipn = mnic(1,jnic)
(  377)          ips = mnic(2,jnic)
(  378)          ipe = mnic(3,jnic)
(  379)          ipw = mnic(4,jnic)
(  380)          IF(lwp) THEN
(  381)             WRITE(numout,9000) jnic
(  382)             WRITE(numout,9010) ip, ipn, ips, ipe, ipw
(  383)             WRITE(numout,9020)
(  384)             DO jnp = 1, mnic(0,jnic)
(  385)                WRITE(numout,9030) jnp,( miic(jnp,ji,jnic)+nimpp-1, mjic(jnp,ji,jnic)+njmpp-1, ji=0,4 )
(  386)             END DO
(  387)          END IF
(  388)          
(  389)          ! format
(  390) 
(  391)  9000 FORMAT(/,'          isolated coastline number= ',i2)
(  392)  9010 FORMAT(/,'          npic=',i4,' npn=',i4,' nps=',i4,' npe=',i4,' npw=',i4)
(  393)  9020 FORMAT(/,'              * ic  point *  point n  *  point s  *  point e ','*  point w  *')
(  394)  9030 FORMAT('         ',i4,' * (',i4,',',i4,') * (',i4,',',i4,') * (',i4,',',i4,') * (',i4,',',i4,') * (',i4,',',i4,') *')






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 8

(  395) 
(  396)       END DO
(  397) 
(  398)       ! 7. Construct the gcfobc array associated with each isolated coastline
(  399)       ! ----------------------------------------------------------------------
(  400) 
(  401)       DO jnic = 1, nbobc-1
(  402) 
(  403)          ! north and south grid-points
(  404)          DO jii = 1, 2
(  405)             DO jnp = 1, mnic(jii,jnic)
(  406)                ii = miic(jnp,jii,jnic)
(  407)                ij = mjic(jnp,jii,jnic)
(  408)                ! take only into account gridpoint of the model domain
(  409)                IF( ii >= nldi+nimpp-1 .AND. ii <= nlci+nimpp-1 .AND.   &
(  410)                    ij >= nldj+njmpp-1 .AND. ij <= nlcj+njmpp-1       ) THEN
(  411)                   iii=ii-nimpp+1
(  412)                   ijj=ij-njmpp+1
(  413)                   gcfobc(iii,ijj-jii+1,jnic) = gcfobc(iii,ijj-jii+1,jnic)   &
(  414)                                              - hur(iii,ijj) * e1u(iii,ijj) / e2u(iii,ijj)
(  415) 
(  416)                END IF
(  417)             END DO
(  418)          END DO
(  419) 
(  420)          ! east and west grid-points
(  421)          DO jii = 3, 4
(  422)             DO jnp = 1, mnic(jii,jnic)
(  423)                ii = miic(jnp,jii,jnic)
(  424)                ij = mjic(jnp,jii,jnic)
(  425)                ! take only into account gridpoint of the model domain
(  426)                IF( ii >= nldi+nimpp-1 .AND. ii <= nlci+nimpp-1 .AND.    &
(  427)                    ij >= nldj+njmpp-1 .AND. ij <= nlcj+njmpp-1 ) THEN
(  428)                   iii=ii-nimpp+1
(  429)                   ijj=ij-njmpp+1
(  430)                   IF( iii-jii+3 == 1 ) THEN
(  431)                      ! cyclic east-west boundary
(  432)                      gcfobc(jpim1    ,ijj,jnic) = gcfobc(jpim1    ,ijj,jnic)   &
(  433)                                                 - hvr(iii,ijj) * e2v(iii,ijj) / e1v(iii,ijj)
(  434)                   ELSE
(  435)                      ! interior points
(  436)                      gcfobc(iii-jii+3,ijj,jnic) = gcfobc(iii-jii+3,ijj,jnic)   &
(  437)                                                 - hvr(iii,ijj) * e2v(iii,ijj) / e1v(iii,ijj)
(  438)                   END IF
(  439)                END IF
(  440)             END DO
(  441)          END DO
(  442) 
(  443)          ! applied bmask to suppress coastal open boundary influence
(  444)          DO jj = 1, jpj
(  445)             DO ji = 1, jpi
(  446)                gcfobc(ji,jj,jnic) = gcfobc(ji,jj,jnic) * bmask(ji,jj)
(  447)             END DO
(  448)          END DO
(  449)          
(  450)       END DO
(  451) 
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 9

(  453)       ! 8. check the grid point which value controls the isolated coastline potential
(  454)       !    Note: in order to activate those tests you need to make zwb a global array,
(  455)       !    which is not done usually to spare memory.
(  456)       !    n.b. here at least 2 open boundaries
(  457)       ! ------------------------------------------------------------------------------
(  458) !
(  459) ! east open boundary:
(  460) !     IF( nieob /= 0 ) THEN
(  461) !   east open & south open  :                              Ed === Sf
(  462) !         IF( njsob /= 0 ) THEN
(  463) !             IF( zwb(nieob,jped) /= zwb(jpsf,njsob) ) THEN
(  464) !                 IF(lwp)WRITE(numout,*) ' E R R O R : east d # south f'
(  465) !             END IF
(  466) !   east open, south closed & west open :                  Ed === Wd
(  467) !           ELSEIF( niwob /= 0 ) THEN
(  468) !             IF( zwb(nieob,jped) /= zwb(niwob,jpwd) ) THEN
(  469) !                 IF(lwp)WRITE(numout,*) ' E R R O R : east d # west d'
(  470) !             END IF
(  471) !   east open, south closed, west closed & north open :    Ed === Nd
(  472) !           ELSEIF( njnob /= 0 ) THEN
(  473) !             IF( zwb(nieob,jped) /= zwb(jpnd,njnob) ) THEN
(  474) !                 IF(lwp)WRITE(numout,*) ' E R R O R : east d # north d'
(  475) !             END IF
(  476) !         END IF
(  477) !   east open & north open :                               Ef === Nf
(  478) !         IF( njnob /= 0 ) THEN
(  479) !             IF( zwb(nieob,jpef) /= zwb(jpnf,njnob) ) THEN
(  480) !                 IF(lwp)WRITE(numout,*) ' E R R O R : east f # north f'
(  481) !             END IF
(  482) !   east open, north closed & west open :                  Ef === Wf
(  483) !           ELSEIF( niwob /= 0 ) THEN
(  484) !             IF( zwb(nieob,jpef) /= zwb(niwob,jpwf) ) THEN
(  485) !                 IF(lwp)WRITE(numout,*) ' E R R O R : east f # west f'
(  486) !             END IF
(  487) !   east open, north closed, west closed & south open :    Ef === Sd
(  488) !           ELSEIF( njsob /= 0 ) THEN
(  489) !             IF( zwb(nieob,jpef) /= zwb(jpsd,njnob) ) THEN
(  490) !                 IF(lwp)WRITE(numout,*) ' E R R O R : east f # south d'
(  491) !             END IF
(  492) !         END IF
(  493) !
(  494) ! east closed 
(  495) !       ELSE
(  496) ! east closed, south open
(  497) !         IF( njsob /= 0 ) THEN
(  498) !   east closed, south open & west open :                  Sd === Wd
(  499) !             IF( niwob /= 0 ) THEN
(  500) !                 IF( zwb(jpsd,njsob) /= zwb(niwob,jpwd) ) THEN
(  501) !                     IF(lwp)WRITE(numout,*) ' E R R O R :',
(  502) !    $                                       ' south d # west d'
(  503) !                 END IF
(  504) !   east closed, south open, west closed & north open :    Sd === Nd
(  505) !               ELSEIF( njnob /= 0 ) THEN
(  506) !                 IF( zwb(jpsd,njsob) /= zwb(jpnd,njnob) ) THEN
(  507) !                     IF(lwp)WRITE(numout,*) ' E R R O R : ',
(  508) !    $                                       ' south d # north d'
(  509) !                 END IF
(  510) !             END IF






PGF90 (Version     10.5)          09/23/2010  09:50:06      page 10

(  511) !   south open, east closed & north open :                 Sf === Nf
(  512) !             IF( njnob /= 0 ) THEN
(  513) !                 IF( zwb(jpsf,njsob) /= zwb(jpnf,njnob) ) THEN
(  514) !                     IF(lwp)WRITE(numout,*) ' E R R O R : ',
(  515) !    $                                       ' south f # north f'
(  516) !                 END IF
(  517) !   south open, east closed, north closed & west open :    Sf === Wf
(  518) !               ELSEIF( niwob /= 0 ) THEN
(  519) !                 IF( zwb(jpsf,njsob) /= zwb(niwob,jpwf) ) THEN
(  520) !                     IF(lwp)WRITE(numout,*) ' E R R O R : ',
(  521) !    $                                       ' south f # west f'
(  522) !                 END IF
(  523) !             END IF
(  524) !
(  525) ! east & south closed ==> north & west open :              Nd === Wf
(  526) !                                                              Nf === Wd
(  527) !           ELSE
(  528) !             IF( zwb(jpnd,njnob) /= zwb(niwob,jpwf) ) THEN
(  529) !                 IF(lwp)WRITE(numout,*) ' E R R O R : north d # west f'
(  530) !             END IF
(  531) !             IF( zwb(jpnf,njnob) /= zwb(niwob,jpwd) ) THEN
(  532) !                 IF(lwp)WRITE(numout,*) ' E R R O R : north f # west d'
(  533) !             END IF
(  534) !         END IF
(  535) !
(  536) !     END IF
(  537) !
(  538) !
(  539) 
(  540)    END SUBROUTINE obc_dom
(  541) #else
(  542)    !!=================================================================================
(  543)    !!                       ***  MODULE  obcdom  ***
(  544)    !! Space domain :  get all the isolated coastline points needed to resolve the 
(  545)    !!                 barotropic streamfunction elliptic equation associated with 
(  546)    !!                 the open boundaries.
(  547)    !!=================================================================================
(  548) CONTAINS
(  549) 
(  550)    SUBROUTINE obc_dom                 
(  551) 
(  552)       ! No isolated coastline OR No Open Boundaries ==> empty routine
(  553) 
(  554)    END SUBROUTINE obc_dom
(  555) #endif
(  556) 
(  557) END MODULE obcdom

















PGF90 (Version     10.5)          09/23/2010  09:50:06      page 11

(    1) # 1 "obcdom.F90"
(    1) MODULE obcdom
(    2)    !!=================================================================================
(    3)    !!                       ***  MODULE  obcdom  ***
(    4)    !! Space domain  :  get all the isolated coastline points needed to resolve the 
(    5)    !!                  barotropic streamfunction elliptic equation associated with 
(    6)    !!                  the open boundaries.
(    7)    !!=================================================================================
(    8) # 542
(  542)    !!=================================================================================
(  543)    !!                       ***  MODULE  obcdom  ***
(  544)    !! Space domain :  get all the isolated coastline points needed to resolve the 
(  545)    !!                 barotropic streamfunction elliptic equation associated with 
(  546)    !!                 the open boundaries.
(  547)    !!=================================================================================
(  548) CONTAINS
(  549) 
(  550)    SUBROUTINE obc_dom                 
(  551) 
(  552)       ! No isolated coastline OR No Open Boundaries ==> empty routine
(  553) 
(  554)    END SUBROUTINE obc_dom
(  555) 
(  556) # 557
(  557) END MODULE obcdom
