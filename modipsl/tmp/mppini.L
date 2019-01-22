


PGF90 (Version     10.5)          09/23/2010  09:50:05      page 1

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

Filename: mppini.F90

(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) MODULE mppini
(    3)    !!==============================================================================
(    4)    !!                       ***  MODULE mppini   ***
(    5)    !! Ocean initialization : distributed memory computing initialization
(    6)    !!==============================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   mpp_init       : Lay out the global domain over processors
(   10)    !!   mpp_init2      : Lay out the global domain over processors 
(   11)    !!                    with land processor elimination
(   12)    !!   mpp_init_ioispl: IOIPSL initialization in mpp
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE dom_oce         ! ocean space and time domain 
(   16)    USE in_out_manager  ! I/O Manager
(   17)    USE sol_oce         ! ocean elliptic solver
(   18)    USE lib_mpp         ! distribued memory computing library
(   19)    USE ioipsl
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC mpp_init       ! called by opa.F90
(   26)    PUBLIC mpp_init2      ! called by opa.F90
(   27) 
(   28)    !! * Substitutions
(   29) #  include "domzgr_substitute.h90"
(   30)    !!----------------------------------------------------------------------
(   31)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/mppini.F90,v 1.9 2006/04/10 15:46:05 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36) CONTAINS
(   37) 
(   38) #if ! defined key_mpp_mpi   &&   ! defined key_mpp_shmem
(   39)    !!----------------------------------------------------------------------
(   40)    !!   Default option :                            shared memory computing
(   41)    !!----------------------------------------------------------------------
(   42) 
(   43)    SUBROUTINE mpp_init
(   44)       !!----------------------------------------------------------------------
(   45)       !!                  ***  ROUTINE mpp_init  ***
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 2

(   47)       !! ** Purpose :   Lay out the global domain over processors.
(   48)       !!
(   49)       !! ** Method  :   Shared memory computing, set the local processor
(   50)       !!      variables to the value of the global domain
(   51)       !!
(   52)       !! History :
(   53)       !!   9.0  !  04-01  (G. Madec, J.M. Molines)  F90 : free form, north fold jpni >1
(   54)       !!----------------------------------------------------------------------
(   55) 
(   56)       ! No mpp computation
(   57)       nimpp  = 1
(   58)       njmpp  = 1
(   59)       nlci   = jpi
(   60)       nlcj   = jpj
(   61)       nldi   = 1
(   62)       nldj   = 1
(   63)       nlei   = jpi
(   64)       nlej   = jpj
(   65)       nperio = jperio
(   66)       nbondi = 2
(   67)       nbondj = 2
(   68)       nidom  = FLIO_DOM_NONE
(   69) 
(   70)       IF(lwp) THEN
(   71)          WRITE(numout,*)
(   72)          WRITE(numout,*) 'mpp_init(2) : NO massively parallel processing'
(   73)          WRITE(numout,*) '~~~~~~~~~~~: '
(   74)          WRITE(numout,*) '         nperio = ', nperio
(   75)          WRITE(numout,*) '         npolj  = ', npolj
(   76)          WRITE(numout,*) '         nimpp  = ', nimpp
(   77)          WRITE(numout,*) '         njmpp  = ', njmpp
(   78)       ENDIF
(   79) 
(   80)       IF(  jpni /= 1 .OR. jpnj /= 1 .OR. jpnij /= 1 ) THEN
(   81)           IF(lwp)WRITE(numout,cform_err)
(   82)           IF(lwp)WRITE(numout,*) 'equality  jpni = jpnj = jpnij = 1 is not satisfied'
(   83)           IF(lwp)WRITE(numout,*) 'the domain is lay out for distributed memory computing! '
(   84)           nstop = nstop + 1
(   85)       ENDIF
(   86) 
(   87) 
(   88)    END SUBROUTINE mpp_init
(   89) 
(   90) 
(   91)    SUBROUTINE mpp_init2 
(   92)       CALL mpp_init                             ! same routine as mpp_init
(   93)    END SUBROUTINE mpp_init2
(   94) 
(   95) #else
(   96)    !!----------------------------------------------------------------------
(   97)    !!   'key_mpp_mpi'          OR         MPI massively parallel processing
(   98)    !!   'key_mpp_shmem'                 SHMEM massively parallel processing
(   99)    !!----------------------------------------------------------------------
(  100) 
(  101)    SUBROUTINE mpp_init
(  102)       !!----------------------------------------------------------------------
(  103)       !!                  ***  ROUTINE mpp_init  ***
(  104)       !!                    






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 3

(  105)       !! ** Purpose :   Lay out the global domain over processors.
(  106)       !!
(  107)       !! ** Method  :   Global domain is distributed in smaller local domains.
(  108)       !!      Periodic condition is a function of the local domain position
(  109)       !!      (global boundary or neighbouring domain) and of the global
(  110)       !!      periodic
(  111)       !!      Type :         jperio global periodic condition
(  112)       !!                     nperio local  periodic condition
(  113)       !!
(  114)       !! ** Action  : - set domain parameters
(  115)       !!                    nimpp     : longitudinal index 
(  116)       !!                    njmpp     : latitudinal  index
(  117)       !!                    nperio    : lateral condition type 
(  118)       !!                    narea     : number for local area
(  119)       !!                    nlci      : first dimension
(  120)       !!                    nlcj      : second dimension
(  121)       !!                    nbondi    : mark for "east-west local boundary"
(  122)       !!                    nbondj    : mark for "north-south local boundary"
(  123)       !!                    nproc     : number for local processor
(  124)       !!                    noea      : number for local neighboring processor
(  125)       !!                    nowe      : number for local neighboring processor
(  126)       !!                    noso      : number for local neighboring processor
(  127)       !!                    nono      : number for local neighboring processor
(  128)       !!
(  129)       !! History :
(  130)       !!        !  94-11  (M. Guyon)  Original code
(  131)       !!        !  95-04  (J. Escobar, M. Imbard)
(  132)       !!        !  98-02  (M. Guyon)  FETI method
(  133)       !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
(  134)       !!   8.5  !  02-08  (G. Madec)  F90 : free form
(  135)       !!----------------------------------------------------------------------
(  136)       !! * Local variables
(  137)       INTEGER ::   ji, jj, jn               ! dummy loop indices
(  138)       INTEGER ::   &
(  139)          ii, ij, ifreq, il1, il2,        &  ! temporary integers
(  140)          iresti, irestj, ijm1, imil,     &  !    "          "
(  141)          inum                               ! temporary logical unit
(  142) 
(  143)       INTEGER, DIMENSION(jpni,jpnj) ::   &
(  144)          iimppt, ijmppt, ilcit, ilcjt       ! temporary workspace
(  145)       REAL(wp) ::   zidom, zjdom            ! temporary scalars
(  146)       !!----------------------------------------------------------------------
(  147) 
(  148) #if defined key_mpp_shmem
(  149)       IF(lwp) WRITE(numout,*)
(  150)       IF(lwp) WRITE(numout,*) 'mpp_init : Message Passing PVM T3E + SHMEM'
(  151)       IF(lwp) WRITE(numout,*) '~~~~~~~~'
(  152) 
(  153)       CALL mppshmem                           ! Initialisation of shmem array
(  154) 
(  155) #endif
(  156) #if defined key_mpp_mpi
(  157)       IF(lwp) WRITE(numout,*)
(  158)       IF(lwp) WRITE(numout,*) 'mpp_init : Message Passing MPI'
(  159)       IF(lwp) WRITE(numout,*) '~~~~~~~~'
(  160) #endif
(  161) 
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 4

(  163)       !  1. Dimension arrays for subdomains
(  164)       ! -----------------------------------
(  165)       !  Computation of local domain sizes ilcit() ilcjt()
(  166)       !  These dimensions depend on global sizes jpni,jpnj and jpiglo,jpjglo
(  167)       !  The subdomains are squares leeser than or equal to the global
(  168)       !  dimensions divided by the number of processors minus the overlap
(  169)       !  array (cf. par_oce.F90).
(  170)       
(  171)       nreci  = 2 * jpreci
(  172)       nrecj  = 2 * jprecj
(  173)       iresti = MOD( jpiglo - nreci , jpni )
(  174)       irestj = MOD( jpjglo - nrecj , jpnj )
(  175) 
(  176)       IF(  iresti == 0 )   iresti = jpni
(  177)       DO jj = 1, jpnj
(  178)          DO ji = 1, iresti
(  179)             ilcit(ji,jj) = jpi
(  180)          END DO
(  181)          DO ji = iresti+1, jpni
(  182)             ilcit(ji,jj) = jpi -1
(  183)          END DO
(  184)       END DO
(  185)       
(  186)       IF( irestj == 0 )   irestj = jpnj
(  187)       DO ji = 1, jpni
(  188)          DO jj = 1, irestj
(  189)             ilcjt(ji,jj) = jpj
(  190)          END DO
(  191)          DO jj = irestj+1, jpnj
(  192)             ilcjt(ji,jj) = jpj -1
(  193)          END DO
(  194)       END DO
(  195)       
(  196)       IF(lwp) THEN
(  197)          WRITE(numout,*)
(  198)          WRITE(numout,*) '           defines mpp subdomains'
(  199)          WRITE(numout,*) '           ----------------------'
(  200)          WRITE(numout,*) '              iresti=',iresti,' irestj=',irestj
(  201)          WRITE(numout,*) '              jpni=',jpni,' jpnj=',jpnj
(  202)          ifreq = 4
(  203)          il1   = 1
(  204)          DO jn = 1, (jpni-1)/ifreq+1
(  205)             il2 = MIN( jpni, il1+ifreq-1 )
(  206)             WRITE(numout,*)
(  207)             WRITE(numout,9201) (ji,ji = il1,il2)
(  208)             WRITE(numout,9200) ('***',ji = il1,il2-1)
(  209)             DO jj = 1, jpnj
(  210)                WRITE(numout,9203) ('   ',ji = il1,il2-1)
(  211)                WRITE(numout,9202) jj, ( ilcit(ji,jj),ilcjt(ji,jj),ji = il1,il2 )
(  212)                WRITE(numout,9203) ('   ',ji = il1,il2-1)
(  213)                WRITE(numout,9200) ('***',ji = il1,il2-1)
(  214)             END DO
(  215)             il1 = il1+ifreq
(  216)          END DO
(  217)  9200    FORMAT('     ***',20('*************',a3))
(  218)  9203    FORMAT('     *     ',20('         *   ',a3))
(  219)  9201    FORMAT('        ',20('   ',i3,'          '))
(  220)  9202    FORMAT(' ',i3,' *  ',20(i3,'  x',i3,'   *   '))






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 5

(  221)       ENDIF
(  222) 
(  223)       zidom = nreci
(  224)       DO ji = 1, jpni
(  225)          zidom = zidom + ilcit(ji,1) - nreci
(  226)       END DO
(  227)       IF(lwp) WRITE(numout,*)
(  228)       IF(lwp) WRITE(numout,*)' sum ilcit(i,1) = ', zidom, ' jpiglo = ', jpiglo
(  229)       
(  230)       zjdom = nrecj
(  231)       DO jj = 1, jpnj
(  232)          zjdom = zjdom + ilcjt(1,jj) - nrecj
(  233)       END DO
(  234)       IF(lwp) WRITE(numout,*)' sum ilcit(1,j) = ', zjdom, ' jpjglo = ', jpjglo
(  235)       IF(lwp) WRITE(numout,*)
(  236)       
(  237) 
(  238)       !  2. Index arrays for subdomains
(  239)       ! -------------------------------
(  240)       
(  241)       iimppt(:,:) = 1
(  242)       ijmppt(:,:) = 1
(  243)       
(  244)       IF( jpni > 1 ) THEN
(  245)          DO jj = 1, jpnj
(  246)             DO ji = 2, jpni
(  247)                iimppt(ji,jj) = iimppt(ji-1,jj) + ilcit(ji-1,jj) - nreci
(  248)             END DO
(  249)          END DO
(  250)       ENDIF
(  251) 
(  252)       IF( jpnj > 1 ) THEN
(  253)          DO jj = 2, jpnj
(  254)             DO ji = 1, jpni
(  255)                ijmppt(ji,jj) = ijmppt(ji,jj-1)+ilcjt(ji,jj-1)-nrecj
(  256)             END DO
(  257)          END DO
(  258)       ENDIF
(  259)       
(  260)       ! 3. Subdomain description
(  261)       ! ------------------------
(  262) 
(  263)       DO jn = 1, jpnij
(  264)          ii = 1 + MOD( jn-1, jpni )
(  265)          ij = 1 + (jn-1) / jpni
(  266)          nimppt(jn) = iimppt(ii,ij)
(  267)          njmppt(jn) = ijmppt(ii,ij)
(  268)          nlcit (jn) = ilcit (ii,ij)     
(  269)          nlci       = nlcit (jn)     
(  270)          nlcjt (jn) = ilcjt (ii,ij)     
(  271)          nlcj       = nlcjt (jn)
(  272)          nbondj = -1                                   ! general case
(  273)          IF( jn   >  jpni          )   nbondj = 0      ! first row of processor
(  274)          IF( jn   >  (jpnj-1)*jpni )   nbondj = 1      ! last  row of processor
(  275)          IF( jpnj == 1             )   nbondj = 2      ! one processor only in j-direction
(  276)          ibonjt(jn) = nbondj
(  277)          
(  278)          nbondi = 0                                    ! 






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 6

(  279)          IF( MOD( jn, jpni ) == 1 )   nbondi = -1      !
(  280)          IF( MOD( jn, jpni ) == 0 )   nbondi =  1      !
(  281)          IF( jpni            == 1 )   nbondi =  2      ! one processor only in i-direction
(  282)          ibonit(jn) = nbondi
(  283)          
(  284)          nldi =  1   + jpreci
(  285)          nlei = nlci - jpreci
(  286)          IF( nbondi == -1 .OR. nbondi == 2 )   nldi = 1
(  287)          IF( nbondi ==  1 .OR. nbondi == 2 )   nlei = nlci
(  288)          nldj =  1   + jprecj
(  289)          nlej = nlcj - jprecj
(  290)          IF( nbondj == -1 .OR. nbondj == 2 )   nldj = 1
(  291)          IF( nbondj ==  1 .OR. nbondj == 2 )   nlej = nlcj
(  292)          nldit(jn) = nldi
(  293)          nleit(jn) = nlei
(  294)          nldjt(jn) = nldj
(  295)          nlejt(jn) = nlej
(  296)       END DO
(  297)       
(  298) 
(  299)       ! 4. From global to local
(  300)       ! -----------------------
(  301) 
(  302)       nperio = 0
(  303)       IF( jperio == 2 .AND. nbondj == -1 )   nperio = 2
(  304) 
(  305) 
(  306)       ! 5. Subdomain neighbours
(  307)       ! ----------------------
(  308) 
(  309)       nproc = narea - 1
(  310)       noso  = nproc - jpni
(  311)       nowe  = nproc - 1
(  312)       noea  = nproc + 1
(  313)       nono  = nproc + jpni
(  314)       ! great neighbours
(  315)       npnw = nono - 1
(  316)       npne = nono + 1
(  317)       npsw = noso - 1
(  318)       npse = noso + 1
(  319)       nbsw = 1
(  320)       nbnw = 1
(  321)       IF( MOD( nproc, jpni ) == 0 ) THEN
(  322)          nbsw = 0
(  323)          nbnw = 0
(  324)       ENDIF
(  325)       nbse = 1
(  326)       nbne = 1
(  327)       IF( MOD( nproc, jpni ) == jpni-1 ) THEN
(  328)          nbse = 0
(  329)          nbne = 0
(  330)       ENDIF
(  331)       IF(nproc < jpni) THEN
(  332)          nbsw = 0
(  333)          nbse = 0
(  334)       ENDIF
(  335)       IF( nproc >= (jpnj-1)*jpni ) THEN
(  336)          nbnw = 0






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 7

(  337)          nbne = 0
(  338)       ENDIF
(  339)       nlcj = nlcjt(narea)  
(  340)       nlci = nlcit(narea)  
(  341)       nldi = nldit(narea)
(  342)       nlei = nleit(narea)
(  343)       nldj = nldjt(narea)
(  344)       nlej = nlejt(narea)
(  345)       nbondi = ibonit(narea)
(  346)       nbondj = ibonjt(narea)
(  347)       nimpp  = nimppt(narea)  
(  348)       njmpp  = njmppt(narea)  
(  349) 
(  350)      ! Save processor layout in layout.dat file 
(  351)        IF (lwp) THEN
(  352)         inum = 11    
(  353) 
(  354)         OPEN(inum,FILE='layout.dat')
(  355)         WRITE(inum,'(a)') '   jpnij     jpi     jpj     jpk  jpiglo  jpjglo'
(  356)         WRITE(inum,'(6i8)') jpnij,jpi,jpj,jpk,jpiglo,jpjglo
(  357)         WRITE(inum,'(a)') 'NAREA nlci nlcj nldi nldj nlei nlej nimpp njmpp'
(  358) 
(  359)         DO  jn = 1, jpnij
(  360)          WRITE(inum,'(9i5)') jn, nlcit(jn), nlcjt(jn), &
(  361)                                       nldit(jn), nldjt(jn), &
(  362)                                       nleit(jn), nlejt(jn), &
(  363)                                       nimppt(jn), njmppt(jn)
(  364)         END DO
(  365)         CLOSE(inum)   
(  366)       END IF
(  367) 
(  368) 
(  369)       ! w a r n i n g  narea (zone) /= nproc (processors)!
(  370) 
(  371)       IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) THEN
(  372)          IF( jpni == 1 )THEN
(  373)             nbondi = 2
(  374)             nperio = 1
(  375)          ELSE
(  376)             nbondi = 0
(  377)          ENDIF
(  378)          IF( MOD( narea, jpni ) == 0 ) THEN
(  379)             noea = nproc-(jpni-1)
(  380)             npne = npne-jpni
(  381)             npse = npse-jpni
(  382)          ENDIF
(  383)          IF( MOD( narea, jpni ) == 1 ) THEN
(  384)             nowe = nproc+(jpni-1)
(  385)             npnw = npnw+jpni
(  386)             npsw = npsw+jpni
(  387)          ENDIF
(  388)          nbsw = 1
(  389)          nbnw = 1
(  390)          nbse = 1
(  391)          nbne = 1
(  392)          IF( nproc < jpni ) THEN
(  393)             nbsw = 0
(  394)             nbse = 0






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 8

(  395)          ENDIF
(  396)          IF( nproc >= (jpnj-1)*jpni ) THEN
(  397)             nbnw = 0
(  398)             nbne = 0
(  399)          ENDIF
(  400)       ENDIF
(  401)       npolj = 0
(  402)       IF( jperio == 3 .OR. jperio == 4 ) THEN
(  403)          ijm1 = jpni*(jpnj-1)
(  404)          imil = ijm1+(jpni+1)/2
(  405)          IF( narea > ijm1 ) npolj = 3
(  406)          IF( MOD(jpni,2) == 1 .AND. narea == imil ) npolj = 4
(  407)          IF( npolj == 3 ) nono = jpni*jpnj-narea+ijm1
(  408)       ENDIF
(  409)       IF( jperio == 5 .OR. jperio == 6 ) THEN
(  410)           ijm1 = jpni*(jpnj-1)
(  411)           imil = ijm1+(jpni+1)/2
(  412)           IF( narea > ijm1) npolj = 5
(  413)           IF( MOD(jpni,2) == 1 .AND. narea == imil ) npolj = 6
(  414)           IF( npolj == 5 ) nono = jpni*jpnj-narea+ijm1
(  415)       ENDIF
(  416) 
(  417)       ! FETI method
(  418) 
(  419)       IF( nperio == 1 .AND. nsolv == 3 ) THEN 
(  420)          
(  421)          ! general case : Earth == infinite tube
(  422)          
(  423)          nbnw = 1
(  424)          npnw = narea
(  425)          nbne = 1
(  426)          npne = narea
(  427)          nbsw = 1
(  428)          npsw = (narea-2)
(  429)          nbse = 1
(  430)          npse = (narea-2)
(  431)          
(  432)          ! REAL boundary condition
(  433)          
(  434)          IF( nbondj == -1 .OR. nbondj == 2 ) THEN 
(  435)             nbsw = 0
(  436)             nbse = 0
(  437)          ENDIF
(  438)          
(  439)          IF( nbondj == -1 .OR. nbondj == 2 ) THEN 
(  440)             nbsw = 0
(  441)             nbse = 0
(  442)          ENDIF
(  443)          
(  444)          IF( nbondj == 1  .OR. nbondj == 2 ) THEN 
(  445)             nbnw = 0
(  446)             nbne = 0
(  447)          ENDIF
(  448)       ENDIF
(  449) 
(  450)       ! Periodicity : no corner if nbondi = 2 and nperio != 1
(  451) 
(  452)       IF(lwp) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 9

(  453)          WRITE(numout,*) ' nproc  = ', nproc
(  454)          WRITE(numout,*) ' nowe   = ', nowe  , ' noea   =  ', noea
(  455)          WRITE(numout,*) ' nono   = ', nono  , ' noso   =  ', noso
(  456)          WRITE(numout,*) ' nbondi = ', nbondi
(  457)          WRITE(numout,*) ' nbondj = ', nbondj
(  458)          WRITE(numout,*) ' npolj  = ', npolj
(  459)          WRITE(numout,*) ' nperio = ', nperio
(  460)          WRITE(numout,*) ' nlci   = ', nlci
(  461)          WRITE(numout,*) ' nlcj   = ', nlcj
(  462)          WRITE(numout,*) ' nimpp  = ', nimpp
(  463)          WRITE(numout,*) ' njmpp  = ', njmpp
(  464)          WRITE(numout,*) ' nbse   = ', nbse  , ' npse   = ', npse
(  465)          WRITE(numout,*) ' nbsw   = ', nbsw  , ' npsw   = ', npsw
(  466)          WRITE(numout,*) ' nbne   = ', nbne  , ' npne   = ', npne
(  467)          WRITE(numout,*) ' nbnw   = ', nbnw  , ' npnw   = ', npnw
(  468)       ENDIF
(  469) 
(  470)       IF( nperio == 1 .AND. jpni /= 1 )THEN
(  471)          IF(lwp) WRITE(numout,cform_err)
(  472)          IF(lwp) WRITE(numout,*) ' mpp_init: error on cyclicity'
(  473)          nstop = nstop + 1
(  474)       ENDIF
(  475) 
(  476)       ! Prepare mpp north fold
(  477) 
(  478)       IF (jperio >= 3 .AND. jperio <= 6 .AND. jpni > 1 ) THEN
(  479)          CALL mpp_ini_north
(  480)       END IF
(  481) 
(  482)       ! Prepare NetCDF output file (if necessary)
(  483)       CALL mpp_init_ioipsl
(  484) 
(  485)    END SUBROUTINE mpp_init
(  486) 
(  487) #  include "mppini_2.h90"
(  488) 
(  489)    SUBROUTINE mpp_init_ioipsl
(  490)       !!----------------------------------------------------------------------
(  491)       !!                  ***  ROUTINE mpp_init_ioipsl  ***
(  492)       !!
(  493)       !! ** Purpose :   
(  494)       !!
(  495)       !! ** Method  :   
(  496)       !!
(  497)       !! History :
(  498)       !!   9.0  !  04-03  (G. Madec)  MPP-IOIPSL 
(  499)       !!----------------------------------------------------------------------
(  500)       !! Local declarations
(  501) 
(  502)       INTEGER, DIMENSION(2) ::   &
(  503)          iglo, iloc, iabsf, iabsl, ihals, ihale, idid
(  504)       INTEGER ::   &
(  505)          inbondi, inbondj
(  506)       !!----------------------------------------------------------------------
(  507) 
(  508)       ! The domain is splitted only horizontally along i- or/and j- direction
(  509)       ! So we need at the most only 1D arrays with 2 elements
(  510)       iglo(1) = jpiglo






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 10

(  511)       iglo(2) = jpjglo
(  512)       iloc(1) = nlci
(  513)       iloc(2) = nlcj
(  514)       iabsf(1) = nimppt(narea)
(  515)       iabsf(2) = njmppt(narea)
(  516)       iabsl(:) = iabsf(:) + iloc(:) - 1
(  517)       ihals(1) = jpreci
(  518)       ihals(2) = jprecj
(  519)       ihale(1) = jpreci
(  520)       ihale(2) = jprecj
(  521)       idid(1) = 1
(  522)       idid(2) = 2
(  523) 
(  524)       inbondj = -1                                   ! general case
(  525)       IF( narea   >  jpni          )   inbondj = 0      ! first row of processor
(  526)       IF( narea   >  (jpnj-1)*jpni )   inbondj = 1      ! last  row of processor
(  527)       IF( jpnj    == 1             )   inbondj = 2      ! one processor only in j-direction
(  528)       IF( inbondj == -1 .OR. inbondj == 2 )   ihals(2) = 0
(  529)       IF( inbondj ==  1 .OR. inbondj == 2 )   ihale(2) = 0
(  530)          
(  531)       inbondi = 0                                    ! 
(  532)       IF( MOD( narea, jpni ) == 1 )   inbondi = -1      !
(  533)       IF( MOD( narea, jpni ) == 0 )   inbondi =  1      !
(  534)       IF( jpni               == 1 )   inbondi =  2      ! one processor only in i-direction
(  535)       IF( inbondi == -1 .OR. inbondi == 2 )   ihals(1) = 0
(  536)       IF( inbondi ==  1 .OR. inbondi == 2 )   ihale(1) = 0
(  537) 
(  538) 
(  539)       IF(lwp) THEN
(  540)           WRITE(numout,*) 'mpp_init_ioipsl :   iloc  = ', iloc (1), iloc (2)
(  541)           WRITE(numout,*) '~~~~~~~~~~~~~~~     iabsf = ', iabsf(1), iabsf(2)
(  542)           WRITE(numout,*) '                    ihals = ', ihals(1), ihals(2)
(  543)           WRITE(numout,*) '                    ihale = ', ihale(1), ihale(2)
(  544)       ENDIF
(  545) 
(  546)       CALL flio_dom_set ( jpnij, nproc, idid, iglo, iloc, iabsf, iabsl, ihals, ihale, 'BOX', nidom)
(  547) 
(  548)    END SUBROUTINE mpp_init_ioipsl  
(  549) 
(  550) #endif
(  551) 
(  552)    !!======================================================================
(  553) END MODULE mppini





















PGF90 (Version     10.5)          09/23/2010  09:50:05      page 11

(    1) # 1 "mppini.F90"
(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) MODULE mppini
(    3)    !!==============================================================================
(    4)    !!                       ***  MODULE mppini   ***
(    5)    !! Ocean initialization : distributed memory computing initialization
(    6)    !!==============================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   mpp_init       : Lay out the global domain over processors
(   10)    !!   mpp_init2      : Lay out the global domain over processors 
(   11)    !!                    with land processor elimination
(   12)    !!   mpp_init_ioispl: IOIPSL initialization in mpp
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE dom_oce         ! ocean space and time domain 
(   16)    USE in_out_manager  ! I/O Manager
(   17)    USE sol_oce         ! ocean elliptic solver
(   18)    USE lib_mpp         ! distribued memory computing library
(   19)    USE ioipsl
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC mpp_init       ! called by opa.F90
(   26)    PUBLIC mpp_init2      ! called by opa.F90
(   27) 
(   28)    !! * Substitutions
(   29) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 30 "mppini.F90"
(   30) # 30
(   30)    !!----------------------------------------------------------------------
(   31)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/mppini.F90,v 1.9 2006/04/10 15:46:05 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36) CONTAINS
(   37) 
(   38) # 96
(   96)    !!----------------------------------------------------------------------
(   97)    !!   'key_mpp_mpi'          OR         MPI massively parallel processing
(   98)    !!   'key_mpp_shmem'                 SHMEM massively parallel processing






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 12

(   99)    !!----------------------------------------------------------------------
(  100) 
(  101)    SUBROUTINE mpp_init
(  102)       !!----------------------------------------------------------------------
(  103)       !!                  ***  ROUTINE mpp_init  ***
(  104)       !!                    
(  105)       !! ** Purpose :   Lay out the global domain over processors.
(  106)       !!
(  107)       !! ** Method  :   Global domain is distributed in smaller local domains.
(  108)       !!      Periodic condition is a function of the local domain position
(  109)       !!      (global boundary or neighbouring domain) and of the global
(  110)       !!      periodic
(  111)       !!      Type :         jperio global periodic condition
(  112)       !!                     nperio local  periodic condition
(  113)       !!
(  114)       !! ** Action  : - set domain parameters
(  115)       !!                    nimpp     : longitudinal index 
(  116)       !!                    njmpp     : latitudinal  index
(  117)       !!                    nperio    : lateral condition type 
(  118)       !!                    narea     : number for local area
(  119)       !!                    nlci      : first dimension
(  120)       !!                    nlcj      : second dimension
(  121)       !!                    nbondi    : mark for "east-west local boundary"
(  122)       !!                    nbondj    : mark for "north-south local boundary"
(  123)       !!                    nproc     : number for local processor
(  124)       !!                    noea      : number for local neighboring processor
(  125)       !!                    nowe      : number for local neighboring processor
(  126)       !!                    noso      : number for local neighboring processor
(  127)       !!                    nono      : number for local neighboring processor
(  128)       !!
(  129)       !! History :
(  130)       !!        !  94-11  (M. Guyon)  Original code
(  131)       !!        !  95-04  (J. Escobar, M. Imbard)
(  132)       !!        !  98-02  (M. Guyon)  FETI method
(  133)       !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
(  134)       !!   8.5  !  02-08  (G. Madec)  F90 : free form
(  135)       !!----------------------------------------------------------------------
(  136)       !! * Local variables
(  137)       INTEGER ::   ji, jj, jn               ! dummy loop indices
(  138)       INTEGER ::   &
(  139)          ii, ij, ifreq, il1, il2,        &  ! temporary integers
(  140)          iresti, irestj, ijm1, imil,     &  !    "          "
(  141)          inum                               ! temporary logical unit
(  142) 
(  143)       INTEGER, DIMENSION(jpni,jpnj) ::   &
(  144)          iimppt, ijmppt, ilcit, ilcjt       ! temporary workspace
(  145)       REAL(wp) ::   zidom, zjdom            ! temporary scalars
(  146)       !!----------------------------------------------------------------------
(  147) 
(  148) # 157
(  157)       IF(lwp) WRITE(numout,*)
(  158)       IF(lwp) WRITE(numout,*) 'mpp_init : Message Passing MPI'
(  159)       IF(lwp) WRITE(numout,*) '~~~~~~~~'
(  160) 
(  161) 
(  162) # 163
(  163)       !  1. Dimension arrays for subdomains
(  164)       ! -----------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 13

(  165)       !  Computation of local domain sizes ilcit() ilcjt()
(  166)       !  These dimensions depend on global sizes jpni,jpnj and jpiglo,jpjglo
(  167)       !  The subdomains are squares leeser than or equal to the global
(  168)       !  dimensions divided by the number of processors minus the overlap
(  169)       !  array (cf. par_oce.F90).
(  170)       
(  171)       nreci  = 2 * jpreci
(  172)       nrecj  = 2 * jprecj
(  173)       iresti = MOD( jpiglo - nreci , jpni )
(  174)       irestj = MOD( jpjglo - nrecj , jpnj )
(  175) 
(  176)       IF(  iresti == 0 )   iresti = jpni
(  177)       DO jj = 1, jpnj
(  178)          DO ji = 1, iresti
(  179)             ilcit(ji,jj) = jpi
(  180)          END DO
(  181)          DO ji = iresti+1, jpni
(  182)             ilcit(ji,jj) = jpi -1
(  183)          END DO
(  184)       END DO
(  185)       
(  186)       IF( irestj == 0 )   irestj = jpnj
(  187)       DO ji = 1, jpni
(  188)          DO jj = 1, irestj
(  189)             ilcjt(ji,jj) = jpj
(  190)          END DO
(  191)          DO jj = irestj+1, jpnj
(  192)             ilcjt(ji,jj) = jpj -1
(  193)          END DO
(  194)       END DO
(  195)       
(  196)       IF(lwp) THEN
(  197)          WRITE(numout,*)
(  198)          WRITE(numout,*) '           defines mpp subdomains'
(  199)          WRITE(numout,*) '           ----------------------'
(  200)          WRITE(numout,*) '              iresti=',iresti,' irestj=',irestj
(  201)          WRITE(numout,*) '              jpni=',jpni,' jpnj=',jpnj
(  202)          ifreq = 4
(  203)          il1   = 1
(  204)          DO jn = 1, (jpni-1)/ifreq+1
(  205)             il2 = MIN( jpni, il1+ifreq-1 )
(  206)             WRITE(numout,*)
(  207)             WRITE(numout,9201) (ji,ji = il1,il2)
(  208)             WRITE(numout,9200) ('***',ji = il1,il2-1)
(  209)             DO jj = 1, jpnj
(  210)                WRITE(numout,9203) ('   ',ji = il1,il2-1)
(  211)                WRITE(numout,9202) jj, ( ilcit(ji,jj),ilcjt(ji,jj),ji = il1,il2 )
(  212)                WRITE(numout,9203) ('   ',ji = il1,il2-1)
(  213)                WRITE(numout,9200) ('***',ji = il1,il2-1)
(  214)             END DO
(  215)             il1 = il1+ifreq
(  216)          END DO
(  217)  9200    FORMAT('     ***',20('*************',a3))
(  218)  9203    FORMAT('     *     ',20('         *   ',a3))
(  219)  9201    FORMAT('        ',20('   ',i3,'          '))
(  220)  9202    FORMAT(' ',i3,' *  ',20(i3,'  x',i3,'   *   '))
(  221)       ENDIF
(  222) 






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 14

(  223)       zidom = nreci
(  224)       DO ji = 1, jpni
(  225)          zidom = zidom + ilcit(ji,1) - nreci
(  226)       END DO
(  227)       IF(lwp) WRITE(numout,*)
(  228)       IF(lwp) WRITE(numout,*)' sum ilcit(i,1) = ', zidom, ' jpiglo = ', jpiglo
(  229)       
(  230)       zjdom = nrecj
(  231)       DO jj = 1, jpnj
(  232)          zjdom = zjdom + ilcjt(1,jj) - nrecj
(  233)       END DO
(  234)       IF(lwp) WRITE(numout,*)' sum ilcit(1,j) = ', zjdom, ' jpjglo = ', jpjglo
(  235)       IF(lwp) WRITE(numout,*)
(  236)       
(  237) 
(  238)       !  2. Index arrays for subdomains
(  239)       ! -------------------------------
(  240)       
(  241)       iimppt(:,:) = 1
(  242)       ijmppt(:,:) = 1
(  243)       
(  244)       IF( jpni > 1 ) THEN
(  245)          DO jj = 1, jpnj
(  246)             DO ji = 2, jpni
(  247)                iimppt(ji,jj) = iimppt(ji-1,jj) + ilcit(ji-1,jj) - nreci
(  248)             END DO
(  249)          END DO
(  250)       ENDIF
(  251) 
(  252)       IF( jpnj > 1 ) THEN
(  253)          DO jj = 2, jpnj
(  254)             DO ji = 1, jpni
(  255)                ijmppt(ji,jj) = ijmppt(ji,jj-1)+ilcjt(ji,jj-1)-nrecj
(  256)             END DO
(  257)          END DO
(  258)       ENDIF
(  259)       
(  260)       ! 3. Subdomain description
(  261)       ! ------------------------
(  262) 
(  263)       DO jn = 1, jpnij
(  264)          ii = 1 + MOD( jn-1, jpni )
(  265)          ij = 1 + (jn-1) / jpni
(  266)          nimppt(jn) = iimppt(ii,ij)
(  267)          njmppt(jn) = ijmppt(ii,ij)
(  268)          nlcit (jn) = ilcit (ii,ij)     
(  269)          nlci       = nlcit (jn)     
(  270)          nlcjt (jn) = ilcjt (ii,ij)     
(  271)          nlcj       = nlcjt (jn)
(  272)          nbondj = -1                                   ! general case
(  273)          IF( jn   >  jpni          )   nbondj = 0      ! first row of processor
(  274)          IF( jn   >  (jpnj-1)*jpni )   nbondj = 1      ! last  row of processor
(  275)          IF( jpnj == 1             )   nbondj = 2      ! one processor only in j-direction
(  276)          ibonjt(jn) = nbondj
(  277)          
(  278)          nbondi = 0                                    ! 
(  279)          IF( MOD( jn, jpni ) == 1 )   nbondi = -1      !
(  280)          IF( MOD( jn, jpni ) == 0 )   nbondi =  1      !






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 15

(  281)          IF( jpni            == 1 )   nbondi =  2      ! one processor only in i-direction
(  282)          ibonit(jn) = nbondi
(  283)          
(  284)          nldi =  1   + jpreci
(  285)          nlei = nlci - jpreci
(  286)          IF( nbondi == -1 .OR. nbondi == 2 )   nldi = 1
(  287)          IF( nbondi ==  1 .OR. nbondi == 2 )   nlei = nlci
(  288)          nldj =  1   + jprecj
(  289)          nlej = nlcj - jprecj
(  290)          IF( nbondj == -1 .OR. nbondj == 2 )   nldj = 1
(  291)          IF( nbondj ==  1 .OR. nbondj == 2 )   nlej = nlcj
(  292)          nldit(jn) = nldi
(  293)          nleit(jn) = nlei
(  294)          nldjt(jn) = nldj
(  295)          nlejt(jn) = nlej
(  296)       END DO
(  297)       
(  298) 
(  299)       ! 4. From global to local
(  300)       ! -----------------------
(  301) 
(  302)       nperio = 0
(  303)       IF( jperio == 2 .AND. nbondj == -1 )   nperio = 2
(  304) 
(  305) 
(  306)       ! 5. Subdomain neighbours
(  307)       ! ----------------------
(  308) 
(  309)       nproc = narea - 1
(  310)       noso  = nproc - jpni
(  311)       nowe  = nproc - 1
(  312)       noea  = nproc + 1
(  313)       nono  = nproc + jpni
(  314)       ! great neighbours
(  315)       npnw = nono - 1
(  316)       npne = nono + 1
(  317)       npsw = noso - 1
(  318)       npse = noso + 1
(  319)       nbsw = 1
(  320)       nbnw = 1
(  321)       IF( MOD( nproc, jpni ) == 0 ) THEN
(  322)          nbsw = 0
(  323)          nbnw = 0
(  324)       ENDIF
(  325)       nbse = 1
(  326)       nbne = 1
(  327)       IF( MOD( nproc, jpni ) == jpni-1 ) THEN
(  328)          nbse = 0
(  329)          nbne = 0
(  330)       ENDIF
(  331)       IF(nproc < jpni) THEN
(  332)          nbsw = 0
(  333)          nbse = 0
(  334)       ENDIF
(  335)       IF( nproc >= (jpnj-1)*jpni ) THEN
(  336)          nbnw = 0
(  337)          nbne = 0
(  338)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 16

(  339)       nlcj = nlcjt(narea)  
(  340)       nlci = nlcit(narea)  
(  341)       nldi = nldit(narea)
(  342)       nlei = nleit(narea)
(  343)       nldj = nldjt(narea)
(  344)       nlej = nlejt(narea)
(  345)       nbondi = ibonit(narea)
(  346)       nbondj = ibonjt(narea)
(  347)       nimpp  = nimppt(narea)  
(  348)       njmpp  = njmppt(narea)  
(  349) 
(  350)      ! Save processor layout in layout.dat file 
(  351)        IF (lwp) THEN
(  352)         inum = 11    
(  353) 
(  354)         OPEN(inum,FILE='layout.dat')
(  355)         WRITE(inum,'(a)') '   jpnij     jpi     jpj     jpk  jpiglo  jpjglo'
(  356)         WRITE(inum,'(6i8)') jpnij,jpi,jpj,jpk,jpiglo,jpjglo
(  357)         WRITE(inum,'(a)') 'NAREA nlci nlcj nldi nldj nlei nlej nimpp njmpp'
(  358) 
(  359)         DO  jn = 1, jpnij
(  360)          WRITE(inum,'(9i5)') jn, nlcit(jn), nlcjt(jn), &
(  361)                                       nldit(jn), nldjt(jn), &
(  362)                                       nleit(jn), nlejt(jn), &
(  363)                                       nimppt(jn), njmppt(jn)
(  364)         END DO
(  365)         CLOSE(inum)   
(  366)       END IF
(  367) 
(  368) 
(  369)       ! w a r n i n g  narea (zone) /= nproc (processors)!
(  370) 
(  371)       IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) THEN
(  372)          IF( jpni == 1 )THEN
(  373)             nbondi = 2
(  374)             nperio = 1
(  375)          ELSE
(  376)             nbondi = 0
(  377)          ENDIF
(  378)          IF( MOD( narea, jpni ) == 0 ) THEN
(  379)             noea = nproc-(jpni-1)
(  380)             npne = npne-jpni
(  381)             npse = npse-jpni
(  382)          ENDIF
(  383)          IF( MOD( narea, jpni ) == 1 ) THEN
(  384)             nowe = nproc+(jpni-1)
(  385)             npnw = npnw+jpni
(  386)             npsw = npsw+jpni
(  387)          ENDIF
(  388)          nbsw = 1
(  389)          nbnw = 1
(  390)          nbse = 1
(  391)          nbne = 1
(  392)          IF( nproc < jpni ) THEN
(  393)             nbsw = 0
(  394)             nbse = 0
(  395)          ENDIF
(  396)          IF( nproc >= (jpnj-1)*jpni ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 17

(  397)             nbnw = 0
(  398)             nbne = 0
(  399)          ENDIF
(  400)       ENDIF
(  401)       npolj = 0
(  402)       IF( jperio == 3 .OR. jperio == 4 ) THEN
(  403)          ijm1 = jpni*(jpnj-1)
(  404)          imil = ijm1+(jpni+1)/2
(  405)          IF( narea > ijm1 ) npolj = 3
(  406)          IF( MOD(jpni,2) == 1 .AND. narea == imil ) npolj = 4
(  407)          IF( npolj == 3 ) nono = jpni*jpnj-narea+ijm1
(  408)       ENDIF
(  409)       IF( jperio == 5 .OR. jperio == 6 ) THEN
(  410)           ijm1 = jpni*(jpnj-1)
(  411)           imil = ijm1+(jpni+1)/2
(  412)           IF( narea > ijm1) npolj = 5
(  413)           IF( MOD(jpni,2) == 1 .AND. narea == imil ) npolj = 6
(  414)           IF( npolj == 5 ) nono = jpni*jpnj-narea+ijm1
(  415)       ENDIF
(  416) 
(  417)       ! FETI method
(  418) 
(  419)       IF( nperio == 1 .AND. nsolv == 3 ) THEN 
(  420)          
(  421)          ! general case : Earth == infinite tube
(  422)          
(  423)          nbnw = 1
(  424)          npnw = narea
(  425)          nbne = 1
(  426)          npne = narea
(  427)          nbsw = 1
(  428)          npsw = (narea-2)
(  429)          nbse = 1
(  430)          npse = (narea-2)
(  431)          
(  432)          ! REAL boundary condition
(  433)          
(  434)          IF( nbondj == -1 .OR. nbondj == 2 ) THEN 
(  435)             nbsw = 0
(  436)             nbse = 0
(  437)          ENDIF
(  438)          
(  439)          IF( nbondj == -1 .OR. nbondj == 2 ) THEN 
(  440)             nbsw = 0
(  441)             nbse = 0
(  442)          ENDIF
(  443)          
(  444)          IF( nbondj == 1  .OR. nbondj == 2 ) THEN 
(  445)             nbnw = 0
(  446)             nbne = 0
(  447)          ENDIF
(  448)       ENDIF
(  449) 
(  450)       ! Periodicity : no corner if nbondi = 2 and nperio != 1
(  451) 
(  452)       IF(lwp) THEN
(  453)          WRITE(numout,*) ' nproc  = ', nproc
(  454)          WRITE(numout,*) ' nowe   = ', nowe  , ' noea   =  ', noea






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 18

(  455)          WRITE(numout,*) ' nono   = ', nono  , ' noso   =  ', noso
(  456)          WRITE(numout,*) ' nbondi = ', nbondi
(  457)          WRITE(numout,*) ' nbondj = ', nbondj
(  458)          WRITE(numout,*) ' npolj  = ', npolj
(  459)          WRITE(numout,*) ' nperio = ', nperio
(  460)          WRITE(numout,*) ' nlci   = ', nlci
(  461)          WRITE(numout,*) ' nlcj   = ', nlcj
(  462)          WRITE(numout,*) ' nimpp  = ', nimpp
(  463)          WRITE(numout,*) ' njmpp  = ', njmpp
(  464)          WRITE(numout,*) ' nbse   = ', nbse  , ' npse   = ', npse
(  465)          WRITE(numout,*) ' nbsw   = ', nbsw  , ' npsw   = ', npsw
(  466)          WRITE(numout,*) ' nbne   = ', nbne  , ' npne   = ', npne
(  467)          WRITE(numout,*) ' nbnw   = ', nbnw  , ' npnw   = ', npnw
(  468)       ENDIF
(  469) 
(  470)       IF( nperio == 1 .AND. jpni /= 1 )THEN
(  471)          IF(lwp) WRITE(numout,cform_err)
(  472)          IF(lwp) WRITE(numout,*) ' mpp_init: error on cyclicity'
(  473)          nstop = nstop + 1
(  474)       ENDIF
(  475) 
(  476)       ! Prepare mpp north fold
(  477) 
(  478)       IF (jperio >= 3 .AND. jperio <= 6 .AND. jpni > 1 ) THEN
(  479)          CALL mpp_ini_north
(  480)       END IF
(  481) 
(  482)       ! Prepare NetCDF output file (if necessary)
(  483)       CALL mpp_init_ioipsl
(  484) 
(  485)    END SUBROUTINE mpp_init
(  486) 
(  487) # 1 "./mppini_2.h90"
(    1)    SUBROUTINE mpp_init2
(    2)       !!----------------------------------------------------------------------
(    3)       !!                  ***  ROUTINE mpp_init2  ***
(    4)       !!
(    5)       !! * Purpose :   Lay out the global domain over processors.
(    6)       !!     FOR USING THIS VERSION, A PREPROCESSING TRAITMENT IS RECOMMENDED
(    7)       !!     FOR DEFINING BETTER CUTTING OUT.
(    8)       !!       This routine is used with a the bathymetry file.
(    9)       !!       In this version, the land processors are avoided and the adress
(   10)       !!     processor (nproc, narea,noea, ...) are calculated again.
(   11)       !!     The jpnij parameter can be lesser than jpni x jpnj
(   12)       !!     and this jpnij parameter must be calculated before with an
(   13)       !!     algoritmic preprocessing program.
(   14)       !!
(   15)       !! ** Method  :   Global domain is distributed in smaller local domains.
(   16)       !!      Periodic condition is a function of the local domain position
(   17)       !!      (global boundary or neighbouring domain) and of the global
(   18)       !!      periodic
(   19)       !!      Type :         jperio global periodic condition
(   20)       !!                     nperio local  periodic condition
(   21)       !!
(   22)       !! ** Action :        nimpp     : longitudinal index 
(   23)       !!                    njmpp     : latitudinal  index
(   24)       !!                    nperio    : lateral condition type 
(   25)       !!                    narea     : number for local area






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 19

(   26)       !!                    nlci      : first dimension
(   27)       !!                    nlcj      : second dimension
(   28)       !!                    nproc     : number for local processor
(   29)       !!                    noea      : number for local neighboring processor
(   30)       !!                    nowe      : number for local neighboring processor
(   31)       !!                    noso      : number for local neighboring processor
(   32)       !!                    nono      : number for local neighboring processor
(   33)       !!
(   34)       !! History :
(   35)       !!        !  94-11  (M. Guyon)  Original code
(   36)       !!        !  95-04  (J. Escobar, M. Imbard)
(   37)       !!        !  98-02  (M. Guyon)  FETI method
(   38)       !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
(   39)       !!   9.0  !  04-01  (G. Madec, J.M Molines)  F90 : free form , north fold jpni > 1
(   40)       !!----------------------------------------------------------------------
(   41)       !! * Modules used
(   42)       USE ioipsl
(   43) 
(   44)       !! Local variables
(   45)       CHARACTER (len=25) ::               &  ! temporary name
(   46)                 clname , clvar               ! filename and cdf variable name for bathy
(   47)       LOGICAL ::   llbon                      ! check the existence of bathy files
(   48)       INTEGER :: ji, jj, jn, jproc, jarea     ! dummy loop indices
(   49)       INTEGER ::   inum = 11                  ! temporary logical unit
(   50)       INTEGER ::   &
(   51)          ii, ij, ifreq, il1, il2,          &  ! temporary integers
(   52)          icont, ili, ilj,                  &  !    "          "
(   53)          isurf, ijm1, imil,                &  !    "          "
(   54)          iino, ijno, iiso, ijso,           &  !    "          " 
(   55)          iiea, ijea, iiwe, ijwe,           &  !    "          "
(   56)          iinw, ijnw, iine, ijne,           &  !    "          "
(   57)          iisw, ijsw, iise, ijse,           &  !    "          "
(   58)          iresti, irestj, iproc                !    "          "
(   59)       INTEGER, DIMENSION(jpnij) ::   &
(   60)          iin, ijn          
(   61)       INTEGER, DIMENSION(jpni,jpnj) ::   &
(   62)          iimppt, ijmppt, ilci  , ilcj  ,   &  ! temporary workspace
(   63)          ipproc, ibondj, ibondi, ipolj ,   &  !    "           "
(   64)          ilei  , ilej  , ildi  , ildj  ,   &  !    "           "
(   65)          ioea  , iowe  , ioso  , iono  ,   &  !    "           "
(   66)          ione  , ionw  , iose  , iosw  ,   &  !    "           "
(   67)          ibne  , ibnw  , ibse  , ibsw         !    "           "
(   68)       INTEGER  ::   &
(   69)          ipi, ipj, ipk,              &  ! temporary integers
(   70)          itime                          !    "          "
(   71)       INTEGER, DIMENSION (1) ::   istep
(   72) 
(   73)       INTEGER, DIMENSION(jpiglo,jpjglo) ::   &
(   74)          imask                                ! temporary global workspace
(   75) 
(   76)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(   77)          zlamt, zphit, zdta                   ! temporary data workspace
(   78)       REAL(wp), DIMENSION(jpk) ::   &   
(   79)          zdept                                ! temporary workspace (NetCDF read)
(   80)       REAL(wp) ::   zidom , zjdom,   &        ! temporary scalars
(   81)          zdt, zdate0
(   82) 
(   83)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 20

(   84)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   85)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/mppini_2.h90,v 1.11 2006/03/10 10:55:34 opalod Exp $ 
(   86)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   87)       !!----------------------------------------------------------------------
(   88) 
(   89) # 99
(   99)       IF(lwp)WRITE(numout,*)
(  100)       IF(lwp)WRITE(numout,*) 'mpp_init : Message Passing MPI'
(  101)       IF(lwp)WRITE(numout,*) '~~~~~~~~'
(  102)       IF(lwp)WRITE(numout,*) ' '
(  103) 
(  104) 
(  105) # 106
(  106)       IF( jpni*jpnj < jpnij ) THEN
(  107)          IF(lwp) WRITE(numout,cform_err)
(  108)          IF(lwp) WRITE(numout,*) ' jpnij > jpni x jpnj impossible'
(  109)          nstop = nstop + 1
(  110)       ENDIF
(  111) 
(  112) 
(  113)       ! 0. initialisation
(  114)       ! -----------------
(  115) 
(  116)       ! open the file
(  117)          IF ( lk_zps ) THEN 
(  118)             clname = 'bathy_meter.nc'         ! Meter bathy in case of partial steps
(  119)             clvar = 'Bathymetry'
(  120)          ELSE
(  121)             clname = 'bathy_level.nc'                       ! Level bathymetry
(  122)             clvar = 'Bathy_level'
(  123)          ENDIF
(  124) 
(  125) # 130
(  130)          INQUIRE( FILE=clname, EXIST=llbon )
(  131)       IF( llbon ) THEN
(  132)             IF(lwp) WRITE(numout,*)
(  133)             IF(lwp) WRITE(numout,*) '         read bathymetry in ', clname
(  134)             IF(lwp) WRITE(numout,*)
(  135)             itime = 1
(  136)             ipi = jpidta
(  137)             ipj = jpjdta
(  138)             ipk = 1
(  139)             zdt = rdt
(  140) 
(  141)             CALL flinopen( clname, 1, jpidta, 1, jpjdta, .FALSE.,   &
(  142)                            ipi, ipj, ipk, zlamt, zphit, zdept, itime, istep, zdate0, zdt, inum )
(  143)             CALL flinget( inum, clvar, jpidta, jpjdta, 1,   &
(  144)                           itime, 1, 1, 1, jpidta, 1, jpjdta, zdta(:,:) )
(  145)             CALL flinclo( inum )
(  146)       ELSE
(  147)          IF(lwp) WRITE(numout,cform_err)
(  148)          IF(lwp) WRITE(numout,*)'    mppini_2 : unable to read the file ', clname
(  149)          nstop = nstop + 1
(  150)       ENDIF
(  151) 
(  152)       ! land/sea mask over the global/zoom domain
(  153) 
(  154)       imask(:,:)=1






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 21

(  155)       WHERE ( zdta(jpizoom:(jpizoom+jpiglo-1),jpjzoom:(jpjglo+jpjzoom-1)) <= 0. ) imask = 0
(  156) 
(  157)       !  1. Dimension arrays for subdomains
(  158)       ! -----------------------------------
(  159) 
(  160)       !  Computation of local domain sizes ilci() ilcj()
(  161)       !  These dimensions depend on global sizes jpni,jpnj and jpiglo,jpjglo
(  162)       !  The subdomains are squares leeser than or equal to the global
(  163)       !  dimensions divided by the number of processors minus the overlap
(  164)       !  array.
(  165) 
(  166)       nreci=2*jpreci
(  167)       nrecj=2*jprecj
(  168)       iresti = 1 + MOD( jpiglo - nreci -1 , jpni )
(  169)       irestj = 1 + MOD( jpjglo - nrecj -1 , jpnj )
(  170) 
(  171)       ilci(1:iresti      ,:) = jpi
(  172)       ilci(iresti+1:jpni ,:) = jpi-1
(  173) 
(  174)       ilcj(:,      1:irestj) = jpj
(  175)       ilcj(:, irestj+1:jpnj) = jpj-1
(  176) 
(  177)       IF(lwp) WRITE(numout,*)
(  178)       IF(lwp) WRITE(numout,*) ' mpp_init2: defines mpp subdomains'
(  179)       IF(lwp) WRITE(numout,*) ' ~~~~~~  ----------------------'
(  180)       IF(lwp) WRITE(numout,*)
(  181)       IF(lwp) WRITE(numout,*) 'iresti=',iresti,' irestj=',irestj
(  182)       IF(lwp) WRITE(numout,*)
(  183)       IF(lwp) WRITE(numout,*) 'jpni=',jpni,' jpnj=',jpnj
(  184) 
(  185)       zidom = nreci + sum(ilci(:,1) - nreci ) 
(  186)       IF(lwp) WRITE(numout,*)
(  187)       IF(lwp) WRITE(numout,*)' sum ilci(i,1)=',zidom,' jpiglo=',jpiglo
(  188) 
(  189)       zjdom = nrecj + sum(ilcj(1,:) - nrecj ) 
(  190)       IF(lwp) WRITE(numout,*) ' sum ilcj(1,j)=',zjdom,' jpjglo=',jpjglo
(  191)       IF(lwp) WRITE(numout,*)
(  192) 
(  193) 
(  194)       !  2. Index arrays for subdomains
(  195)       ! -------------------------------
(  196) 
(  197)       iimppt(:,:) = 1
(  198)       ijmppt(:,:) = 1
(  199)       ipproc(:,:) = -1
(  200) 
(  201)       IF( jpni > 1 )THEN
(  202)          DO jj = 1, jpnj
(  203)             DO ji = 2, jpni
(  204)                iimppt(ji,jj) = iimppt(ji-1,jj) + ilci(ji-1,jj) - nreci
(  205)             END DO
(  206)          END DO
(  207)       ENDIF
(  208) 
(  209)       IF( jpnj > 1 )THEN
(  210)          DO jj = 2, jpnj
(  211)             DO ji = 1, jpni
(  212)                ijmppt(ji,jj) = ijmppt(ji,jj-1) + ilcj(ji,jj-1) - nrecj






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 22

(  213)             END DO
(  214)          END DO
(  215)       ENDIF
(  216) 
(  217) 
(  218)       ! 3. Subdomain description in the Regular Case
(  219)       ! --------------------------------------------
(  220) 
(  221)       nperio = 0
(  222)       icont = -1
(  223)       DO jarea = 1, jpni*jpnj
(  224)          ii = 1 + MOD(jarea-1,jpni)
(  225)          ij = 1 +    (jarea-1)/jpni
(  226)          ili = ilci(ii,ij)
(  227)          ilj = ilcj(ii,ij)
(  228) 
(  229)          ibondj(ii,ij) = -1
(  230)          IF( jarea >  jpni          )   ibondj(ii,ij) = 0
(  231)          IF( jarea >  (jpnj-1)*jpni )   ibondj(ii,ij) = 1
(  232)          IF( jpnj  == 1             )   ibondj(ii,ij) = 2
(  233) 
(  234)          ibondi(ii,ij) = 0
(  235)          IF( MOD(jarea,jpni) == 1 )   ibondi(ii,ij) = -1
(  236)          IF( MOD(jarea,jpni) == 0 )   ibondi(ii,ij) =  1
(  237)          IF( jpni            == 1 )   ibondi(ii,ij) =  2
(  238) 
(  239)          ! 2.4 Subdomain neighbors
(  240) 
(  241)          iproc = jarea - 1
(  242)          ioso(ii,ij) = iproc - jpni
(  243)          iowe(ii,ij) = iproc - 1
(  244)          ioea(ii,ij) = iproc + 1
(  245)          iono(ii,ij) = iproc + jpni
(  246)          ildi(ii,ij) = 1 + jpreci
(  247)          ilei(ii,ij) = ili -jpreci
(  248)          ionw(ii,ij) = iono(ii,ij) - 1
(  249)          ione(ii,ij) = iono(ii,ij) + 1
(  250)          iosw(ii,ij) = ioso(ii,ij) - 1
(  251)          iose(ii,ij) = ioso(ii,ij) + 1
(  252)          ibsw(ii,ij) = 1
(  253)          ibnw(ii,ij) = 1
(  254)          IF( MOD(iproc,jpni) == 0 ) THEN
(  255)             ibsw(ii,ij) = 0
(  256)             ibnw(ii,ij) = 0
(  257)          ENDIF
(  258)          ibse(ii,ij) = 1
(  259)          ibne(ii,ij) = 1
(  260)          IF( MOD(iproc,jpni) == jpni-1 ) THEN
(  261)             ibse(ii,ij) = 0
(  262)             ibne(ii,ij) = 0
(  263)          ENDIF
(  264)          IF( iproc < jpni ) THEN
(  265)             ibsw(ii,ij) = 0
(  266)             ibse(ii,ij) = 0
(  267)          ENDIF
(  268)          IF( iproc >= (jpnj-1)*jpni ) THEN
(  269)             ibnw(ii,ij) = 0
(  270)             ibne(ii,ij) = 0






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 23

(  271)          ENDIF
(  272)          IF( ibondi(ii,ij) == -1 .OR. ibondi(ii,ij) == 2 ) ildi(ii,ij) = 1
(  273)          IF( ibondi(ii,ij) ==  1 .OR. ibondi(ii,ij) == 2 ) ilei(ii,ij) = ili
(  274)          ildj(ii,ij) =  1  + jprecj
(  275)          ilej(ii,ij) = ilj - jprecj
(  276)          IF( ibondj(ii,ij) == -1 .OR. ibondj(ii,ij) == 2 ) ildj(ii,ij) = 1
(  277)          IF( ibondj(ii,ij) ==  1 .OR. ibondj(ii,ij) == 2 ) ilej(ii,ij) = ilj
(  278) 
(  279)          ! warning ii*ij (zone) /= nproc (processors)!
(  280) 
(  281)          IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) THEN
(  282)             IF( jpni == 1 )THEN
(  283)                ibondi(ii,ij) = 2
(  284)                nperio = 1
(  285)             ELSE
(  286)                ibondi(ii,ij) = 0
(  287)             ENDIF
(  288)             IF( MOD(jarea,jpni) == 0 ) THEN
(  289)                ioea(ii,ij) = iproc - (jpni-1)
(  290)                ione(ii,ij) = ione(ii,ij) - jpni
(  291)                iose(ii,ij) = iose(ii,ij) - jpni
(  292)             ENDIF
(  293)             IF( MOD(jarea,jpni) == 1 ) THEN
(  294)                iowe(ii,ij) = iproc + jpni - 1
(  295)                ionw(ii,ij) = ionw(ii,ij) + jpni
(  296)                iosw(ii,ij) = iosw(ii,ij) + jpni 
(  297)             ENDIF
(  298)             ibsw(ii,ij) = 1
(  299)             ibnw(ii,ij) = 1
(  300)             ibse(ii,ij) = 1
(  301)             ibne(ii,ij) = 1
(  302)             IF( iproc < jpni ) THEN
(  303)                ibsw(ii,ij) = 0
(  304)                ibse(ii,ij) = 0
(  305)             ENDIF
(  306)             IF( iproc >= (jpnj-1)*jpni ) THEN
(  307)                ibnw(ii,ij) = 0
(  308)                ibne(ii,ij) = 0
(  309)             ENDIF
(  310)          ENDIF
(  311)          ipolj(ii,ij) = 0
(  312)          IF( jperio == 3 .OR. jperio == 4 ) THEN
(  313)             ijm1 = jpni*(jpnj-1)
(  314)             imil = ijm1+(jpni+1)/2
(  315)             IF( jarea > ijm1 ) ipolj(ii,ij) = 3
(  316)             IF( MOD(jpni,2) == 1 .AND. jarea == imil ) ipolj(ii,ij) = 4
(  317)             IF( ipolj(ii,ij) == 3 ) iono(ii,ij) = jpni*jpnj-jarea+ijm1
(  318)          ENDIF
(  319)          IF( jperio == 5 .OR. jperio == 6 ) THEN
(  320)             ijm1 = jpni*(jpnj-1)
(  321)             imil = ijm1+(jpni+1)/2
(  322)             IF( jarea > ijm1) ipolj(ii,ij) = 5
(  323)             IF( MOD(jpni,2) == 1 .AND. jarea == imil ) ipolj(ii,ij) = 6
(  324)             IF( ipolj(ii,ij) == 5) iono(ii,ij) = jpni*jpnj-jarea+ijm1
(  325)          ENDIF
(  326) 
(  327)          isurf = 0
(  328)          DO jj = 1+jprecj, ilj-jprecj






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 24

(  329)             DO  ji = 1+jpreci, ili-jpreci
(  330)                IF( imask(ji+iimppt(ii,ij)-1, jj+ijmppt(ii,ij)-1) == 1) isurf = isurf+1
(  331)             END DO
(  332)          END DO
(  333)          IF(isurf /= 0) THEN
(  334)             icont = icont + 1
(  335)             ipproc(ii,ij) = icont
(  336)             iin(icont+1) = ii
(  337)             ijn(icont+1) = ij
(  338)          ENDIF
(  339)       END DO
(  340) 
(  341)       ! Control
(  342)       IF(icont+1 /= jpnij) THEN
(  343)          IF(lwp) THEN 
(  344)             WRITE(numout,*) ' Eliminate land processors algorithm'
(  345)             WRITE(numout,*)
(  346)             WRITE(numout,*) ' jpni =',jpni,' jpnj =',jpnj
(  347)             WRITE(numout,*) ' jpnij =',jpnij, '< jpni x jpnj' 
(  348)             WRITE(numout,*)
(  349)             WRITE(numout,*) ' E R R O R '
(  350)             WRITE(numout,*) ' ***********, mpp_init2 finds jpnij=',icont+1
(  351)             WRITE(numout,*) ' we stop'
(  352)          ENDIF
(  353)          STOP 'mpp_init2'
(  354)       ENDIF
(  355) 
(  356) 
(  357)       ! 4. Subdomain print
(  358)       ! ------------------
(  359) 
(  360)       IF(lwp) THEN
(  361)          ifreq = 4
(  362)          il1 = 1
(  363)          DO jn = 1,(jpni-1)/ifreq+1
(  364)             il2 = MIN(jpni,il1+ifreq-1)
(  365)             WRITE(numout,*)
(  366)             WRITE(numout,9401) (ji,ji=il1,il2)
(  367)             WRITE(numout,9400) ('***',ji=il1,il2-1)
(  368)             DO jj = 1, jpnj
(  369)                ! WRITE(numout,9400)
(  370)                WRITE(numout,9403) ('   ',ji=il1,il2-1)
(  371)                WRITE(numout,9402) jj, (ilci(ji,jj),ilcj(ji,jj),ji=il1,il2)
(  372)                WRITE(numout,9404) (ipproc(ji,jj),ji=il1,il2)
(  373)                WRITE(numout,9403) ('   ',ji=il1,il2-1)
(  374)                WRITE(numout,9400) ('***',ji=il1,il2-1)
(  375)                ! WRITE(numout,9400)
(  376)             END DO
(  377)             il1 = il1+ifreq
(  378)          END DO
(  379)  9400     FORMAT('     ***',20('*************',a3))
(  380)  9403     FORMAT('     *     ',20('         *   ',a3))
(  381)  9401     FORMAT('        ',20('   ',i3,'          '))
(  382)  9402     FORMAT(' ',i3,' *  ',20(i3,'  x',i3,'   *   '))
(  383)  9404     FORMAT('     *  ',20('      ',i3,'   *   '))
(  384)       ENDIF
(  385) 
(  386) 






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 25

(  387)       ! 5. neighbour treatment
(  388)       ! ----------------------
(  389) 
(  390)       DO jarea = 1, jpni*jpnj
(  391)          iproc = jarea-1
(  392)          ii = 1 + MOD(jarea-1,jpni)
(  393)          ij = 1 +    (jarea-1)/jpni
(  394)          IF( ipproc(ii,ij) == -1 .AND. iono(ii,ij) >= 0   &
(  395)             .AND. iono(ii,ij) <= jpni*jpnj-1 ) THEN
(  396)             iino = 1 + MOD(iono(ii,ij),jpni)
(  397)             ijno = 1 +    (iono(ii,ij))/jpni
(  398)             IF( ibondj(iino,ijno) == 1 ) ibondj(iino,ijno)=2
(  399)             IF( ibondj(iino,ijno) == 0 ) ibondj(iino,ijno) = -1
(  400)          ENDIF
(  401)          IF( ipproc(ii,ij) == -1 .AND. ioso(ii,ij) >= 0   &
(  402)             .AND. ioso(ii,ij) <= jpni*jpnj-1 ) THEN
(  403)             iiso = 1 + MOD(ioso(ii,ij),jpni)
(  404)             ijso = 1 +    (ioso(ii,ij))/jpni
(  405)             IF( ibondj(iiso,ijso) == -1 ) ibondj(iiso,ijso) = 2
(  406)             IF( ibondj(iiso,ijso) ==  0 ) ibondj(iiso,ijso) = 1
(  407)          ENDIF
(  408)          IF( ipproc(ii,ij) == -1 .AND. ioea(ii,ij) >= 0   &
(  409)             .AND. ioea(ii,ij) <= jpni*jpnj-1) THEN
(  410)             iiea = 1 + MOD(ioea(ii,ij),jpni)
(  411)             ijea = 1 +    (ioea(ii,ij))/jpni
(  412)             IF( ibondi(iiea,ijea) == 1 ) ibondi(iiea,ijea) = 2
(  413)             IF( ibondi(iiea,ijea) == 0 ) ibondi(iiea,ijea) = -1
(  414)          ENDIF
(  415)          IF( ipproc(ii,ij) == -1 .AND. iowe(ii,ij) >= 0   &
(  416)             .AND. iowe(ii,ij) <= jpni*jpnj-1) THEN
(  417)             iiwe = 1 + MOD(iowe(ii,ij),jpni)
(  418)             ijwe = 1 +    (iowe(ii,ij))/jpni
(  419)             IF( ibondi(iiwe,ijwe) == -1 ) ibondi(iiwe,ijwe) = 2
(  420)             IF( ibondi(iiwe,ijwe) ==  0 ) ibondi(iiwe,ijwe) = 1
(  421)          ENDIF
(  422)          IF( ipproc(ii,ij) == -1 .AND. ibne(ii,ij) == 1 ) THEN
(  423)             iine = 1 + MOD(ione(ii,ij),jpni)
(  424)             ijne = 1 +    (ione(ii,ij))/jpni
(  425)             IF( ibsw(iine,ijne) == 1 ) ibsw(iine,ijne) = 0
(  426)          ENDIF
(  427)          IF( ipproc(ii,ij) == -1 .AND. ibsw(ii,ij) == 1 ) THEN
(  428)             iisw = 1 + MOD(iosw(ii,ij),jpni)
(  429)             ijsw = 1 +    (iosw(ii,ij))/jpni
(  430)             IF( ibne(iisw,ijsw) == 1 ) ibne(iisw,ijsw) = 0
(  431)          ENDIF
(  432)          IF( ipproc(ii,ij) == -1 .AND. ibnw(ii,ij) == 1 ) THEN
(  433)             iinw = 1 + MOD(ionw(ii,ij),jpni)
(  434)             ijnw = 1 +    (ionw(ii,ij))/jpni
(  435)             IF( ibse(iinw,ijnw) == 1 ) ibse(iinw,ijnw)=0
(  436)          ENDIF
(  437)          IF( ipproc(ii,ij) == -1 .AND. ibse(ii,ij) == 1 ) THEN
(  438)             iise = 1 + MOD(iose(ii,ij),jpni)
(  439)             ijse = 1 +    (iose(ii,ij))/jpni
(  440)             IF( ibnw(iise,ijse) == 1 ) ibnw(iise,ijse) = 0
(  441)          ENDIF
(  442)       END DO
(  443) 
(  444) 






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 26

(  445)       ! 6. Change processor name
(  446)       ! ------------------------
(  447) 
(  448)       nproc = narea-1
(  449)       ii = iin(narea)
(  450)       ij = ijn(narea)
(  451)       IF( ioso(ii,ij) >= 0 .AND. ioso(ii,ij) <= (jpni*jpnj-1) ) THEN 
(  452)          iiso = 1 + MOD(ioso(ii,ij),jpni)
(  453)          ijso = 1 +    (ioso(ii,ij))/jpni
(  454)          noso = ipproc(iiso,ijso)
(  455)       ENDIF
(  456)       IF( iowe(ii,ij) >= 0 .AND. iowe(ii,ij) <= (jpni*jpnj-1) ) THEN 
(  457)          iiwe = 1 + MOD(iowe(ii,ij),jpni)
(  458)          ijwe = 1 +    (iowe(ii,ij))/jpni
(  459)          nowe = ipproc(iiwe,ijwe)
(  460)       ENDIF
(  461)       IF( ioea(ii,ij) >= 0 .AND. ioea(ii,ij) <= (jpni*jpnj-1) ) THEN 
(  462)          iiea = 1 + MOD(ioea(ii,ij),jpni)
(  463)          ijea = 1 +    (ioea(ii,ij))/jpni
(  464)          noea = ipproc(iiea,ijea)
(  465)       ENDIF
(  466)       IF( iono(ii,ij) >= 0 .AND. iono(ii,ij) <= (jpni*jpnj-1) ) THEN 
(  467)          iino = 1 + MOD(iono(ii,ij),jpni)
(  468)          ijno = 1 +    (iono(ii,ij))/jpni
(  469)          nono = ipproc(iino,ijno)
(  470)       ENDIF
(  471)       IF( iose(ii,ij) >= 0 .AND. iose(ii,ij) <= (jpni*jpnj-1) ) THEN 
(  472)          iise = 1 + MOD(iose(ii,ij),jpni)
(  473)          ijse = 1 +    (iose(ii,ij))/jpni
(  474)          npse = ipproc(iise,ijse)
(  475)       ENDIF
(  476)       IF( iosw(ii,ij) >= 0 .AND. iosw(ii,ij) <= (jpni*jpnj-1) ) THEN 
(  477)          iisw = 1 + MOD(iosw(ii,ij),jpni)
(  478)          ijsw = 1 +    (iosw(ii,ij))/jpni
(  479)          npsw = ipproc(iisw,ijsw)
(  480)       ENDIF
(  481)       IF( ione(ii,ij) >= 0 .AND. ione(ii,ij) <= (jpni*jpnj-1) ) THEN 
(  482)          iine = 1 + MOD(ione(ii,ij),jpni)
(  483)          ijne = 1 +    (ione(ii,ij))/jpni
(  484)          npne = ipproc(iine,ijne)
(  485)       ENDIF
(  486)       IF( ionw(ii,ij) >= 0 .AND. ionw(ii,ij) <= (jpni*jpnj-1) ) THEN 
(  487)          iinw = 1 + MOD(ionw(ii,ij),jpni)
(  488)          ijnw = 1 +    (ionw(ii,ij))/jpni
(  489)          npnw = ipproc(iinw,ijnw)
(  490)       ENDIF
(  491)       nbnw = ibnw(ii,ij)
(  492)       nbne = ibne(ii,ij)
(  493)       nbsw = ibsw(ii,ij)
(  494)       nbse = ibse(ii,ij)
(  495)       nlcj = ilcj(ii,ij)  
(  496)       nlci = ilci(ii,ij)  
(  497)       nldi = ildi(ii,ij)
(  498)       nlei = ilei(ii,ij)
(  499)       nldj = ildj(ii,ij)
(  500)       nlej = ilej(ii,ij)
(  501)       nbondi = ibondi(ii,ij)
(  502)       nbondj = ibondj(ii,ij)






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 27

(  503)       nimpp = iimppt(ii,ij)  
(  504)       njmpp = ijmppt(ii,ij)  
(  505)       DO jproc = 1, jpnij
(  506)          ii = iin(jproc)
(  507)          ij = ijn(jproc)
(  508)          nimppt(jproc) = iimppt(ii,ij)  
(  509)          njmppt(jproc) = ijmppt(ii,ij)  
(  510)          nlcjt(jproc) = ilcj(ii,ij)
(  511)          nlcit(jproc) = ilci(ii,ij)
(  512)          nldit(jproc) = ildi(ii,ij)
(  513)          nleit(jproc) = ilei(ii,ij)
(  514)          nldjt(jproc) = ildj(ii,ij)
(  515)          nlejt(jproc) = ilej(ii,ij)
(  516)       END DO
(  517) 
(  518)       ! Save processor layout in ascii file
(  519)       IF (lwp) THEN
(  520)         OPEN(inum,FILE='layout.dat')
(  521)         WRITE(inum,'(6i8)') jpnij,jpi,jpj,jpk,jpiglo,jpjglo
(  522)         WRITE(inum,'(a)') 'NAREA nlci nlcj nldi nldj nlei nlej nimpp njmpp'
(  523) 
(  524)         DO  jproc = 1, jpnij
(  525)          WRITE(inum,'(9i5)') jproc, nlcit(jproc), nlcjt(jproc), &
(  526)                                       nldit(jproc), nldjt(jproc), &
(  527)                                       nleit(jproc), nlejt(jproc), &
(  528)                                       nimppt(jproc), njmppt(jproc)
(  529)         END DO
(  530)         CLOSE(inum)   
(  531)       END IF
(  532) 
(  533) 
(  534)       ! FETI method
(  535) 
(  536)       IF( nperio == 1 .AND. nsolv == 3 ) THEN 
(  537) 
(  538)          ! general CASE : Earth == infinite tube
(  539) 
(  540)          nbnw = 1
(  541)          npnw = narea
(  542)          nbne = 1
(  543)          npne = narea
(  544)          nbsw = 1
(  545)          npsw = (narea-2)
(  546)          nbse = 1
(  547)          npse = (narea-2)
(  548) 
(  549)          ! REAL boundary condition
(  550) 
(  551)          IF( nbondj == -1 .OR. nbondj == 2 ) THEN 
(  552)             nbsw = 0
(  553)             nbse = 0
(  554)          ENDIF
(  555) 
(  556)          IF( nbondj == -1 .OR. nbondj == 2 ) THEN 
(  557)             nbsw = 0
(  558)             nbse = 0
(  559)          ENDIF
(  560) 






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 28

(  561)          IF( nbondj ==  1 .OR. nbondj == 2 ) THEN 
(  562)             nbnw = 0
(  563)             nbne = 0
(  564)          ENDIF
(  565)       ENDIF
(  566) 
(  567)       IF( nperio == 1 .AND.jpni /= 1 ) THEN
(  568)          IF(lwp) WRITE(numout,cform_err)
(  569)          IF(lwp) WRITE(numout,*) ' mpp_init2:  error on cyclicity'
(  570)          nstop = nstop + 1
(  571)       ENDIF
(  572) 
(  573)       ! Prepare mpp north fold
(  574) 
(  575)       IF( jperio >= 3 .AND. jperio <= 6 .AND. jpni > 1 ) THEN
(  576)          CALL mpp_ini_north
(  577)          IF(lwp) WRITE(numout,*) ' mpp_init2 : North fold boundary prepared for jpni >1'
(  578)       ENDIF
(  579) 
(  580)       ! Defined npolj, either 0, 3 , 4 , 5 , 6
(  581)       ! In this case the important thing is that npolj /= 0
(  582)       ! Because if we go through these line it is because jpni >1 and thus
(  583)       ! we must use lbcnorthmpp, which tests only npolj =0 or npolj /= 0
(  584) 
(  585)       npolj = 0
(  586)       ij = ijn(narea)
(  587) 
(  588)       IF( jperio == 3 .OR. jperio == 4 ) THEN
(  589)          IF( ij == jpnj ) npolj = 3
(  590)       ENDIF
(  591) 
(  592)       IF( jperio == 5 .OR. jperio == 6 ) THEN
(  593)          IF( ij == jpnj ) npolj = 5
(  594)       ENDIF
(  595) 
(  596)       ! Prepare NetCDF output file (if necessary)
(  597)       CALL mpp_init_ioipsl
(  598) 
(  599)       ! Periodicity : no corner if nbondi = 2 and nperio != 1
(  600) 
(  601)       IF(lwp) THEN
(  602)          WRITE(numout,*) ' nproc=  ',nproc
(  603)          WRITE(numout,*) ' nowe=   ',nowe
(  604)          WRITE(numout,*) ' noea=   ',noea
(  605)          WRITE(numout,*) ' nono=   ',nono
(  606)          WRITE(numout,*) ' noso=   ',noso
(  607)          WRITE(numout,*) ' nbondi= ',nbondi
(  608)          WRITE(numout,*) ' nbondj= ',nbondj
(  609)          WRITE(numout,*) ' npolj=  ',npolj
(  610)          WRITE(numout,*) ' nperio= ',nperio
(  611)          WRITE(numout,*) ' nlci=   ',nlci
(  612)          WRITE(numout,*) ' nlcj=   ',nlcj
(  613)          WRITE(numout,*) ' nimpp=  ',nimpp
(  614)          WRITE(numout,*) ' njmpp=  ',njmpp
(  615)          WRITE(numout,*) ' nbse=   ',nbse,' npse= ',npse
(  616)          WRITE(numout,*) ' nbsw=   ',nbsw,' npsw= ',npsw
(  617)          WRITE(numout,*) ' nbne=   ',nbne,' npne= ',npne
(  618)          WRITE(numout,*) ' nbnw=   ',nbnw,' npnw= ',npnw






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 29

(  619)       ENDIF
(  620) 
(  621)    END SUBROUTINE mpp_init2
(  622) # 488 "mppini.F90"
(  488) 
(  489) # 489
(  489)    SUBROUTINE mpp_init_ioipsl
(  490)       !!----------------------------------------------------------------------
(  491)       !!                  ***  ROUTINE mpp_init_ioipsl  ***
(  492)       !!
(  493)       !! ** Purpose :   
(  494)       !!
(  495)       !! ** Method  :   
(  496)       !!
(  497)       !! History :
(  498)       !!   9.0  !  04-03  (G. Madec)  MPP-IOIPSL 
(  499)       !!----------------------------------------------------------------------
(  500)       !! Local declarations
(  501) 
(  502)       INTEGER, DIMENSION(2) ::   &
(  503)          iglo, iloc, iabsf, iabsl, ihals, ihale, idid
(  504)       INTEGER ::   &
(  505)          inbondi, inbondj
(  506)       !!----------------------------------------------------------------------
(  507) 
(  508)       ! The domain is splitted only horizontally along i- or/and j- direction
(  509)       ! So we need at the most only 1D arrays with 2 elements
(  510)       iglo(1) = jpiglo
(  511)       iglo(2) = jpjglo
(  512)       iloc(1) = nlci
(  513)       iloc(2) = nlcj
(  514)       iabsf(1) = nimppt(narea)
(  515)       iabsf(2) = njmppt(narea)
(  516)       iabsl(:) = iabsf(:) + iloc(:) - 1
(  517)       ihals(1) = jpreci
(  518)       ihals(2) = jprecj
(  519)       ihale(1) = jpreci
(  520)       ihale(2) = jprecj
(  521)       idid(1) = 1
(  522)       idid(2) = 2
(  523) 
(  524)       inbondj = -1                                   ! general case
(  525)       IF( narea   >  jpni          )   inbondj = 0      ! first row of processor
(  526)       IF( narea   >  (jpnj-1)*jpni )   inbondj = 1      ! last  row of processor
(  527)       IF( jpnj    == 1             )   inbondj = 2      ! one processor only in j-direction
(  528)       IF( inbondj == -1 .OR. inbondj == 2 )   ihals(2) = 0
(  529)       IF( inbondj ==  1 .OR. inbondj == 2 )   ihale(2) = 0
(  530)          
(  531)       inbondi = 0                                    ! 
(  532)       IF( MOD( narea, jpni ) == 1 )   inbondi = -1      !
(  533)       IF( MOD( narea, jpni ) == 0 )   inbondi =  1      !
(  534)       IF( jpni               == 1 )   inbondi =  2      ! one processor only in i-direction
(  535)       IF( inbondi == -1 .OR. inbondi == 2 )   ihals(1) = 0
(  536)       IF( inbondi ==  1 .OR. inbondi == 2 )   ihale(1) = 0
(  537) 
(  538) 
(  539)       IF(lwp) THEN
(  540)           WRITE(numout,*) 'mpp_init_ioipsl :   iloc  = ', iloc (1), iloc (2)






PGF90 (Version     10.5)          09/23/2010  09:50:05      page 30

(  541)           WRITE(numout,*) '~~~~~~~~~~~~~~~     iabsf = ', iabsf(1), iabsf(2)
(  542)           WRITE(numout,*) '                    ihals = ', ihals(1), ihals(2)
(  543)           WRITE(numout,*) '                    ihale = ', ihale(1), ihale(2)
(  544)       ENDIF
(  545) 
(  546)       CALL flio_dom_set ( jpnij, nproc, idid, iglo, iloc, iabsf, iabsl, ihals, ihale, 'BOX', nidom)
(  547) 
(  548)    END SUBROUTINE mpp_init_ioipsl  
(  549) 
(  550) 
(  551) # 552
(  552)    !!======================================================================
(  553) END MODULE mppini
