


PGF90 (Version     10.5)          09/23/2010  09:49:17      page 1

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

Filename: prtctl_trc.F90

(    1) MODULE prtctl_trc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE prtctl   ***
(    4)    !! Ocean system   : print all SUM trends for each processor domain
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7) 
(    8)    USE par_trc_trp
(    9)    USE oce_trc          ! ocean space and time domain variables
(   10)    USE in_out_manager   ! I/O manager
(   11)    USE lib_mpp          ! distributed memory computing
(   12) 
(   13)    IMPLICIT NONE
(   14)    PRIVATE
(   15) 
(   16)    !! * Module declaration
(   17)    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::   &  !:
(   18)       nlditl , nldjtl ,   &  !: first, last indoor index for each i-domain
(   19)       nleitl , nlejtl ,   &  !: first, last indoor index for each j-domain
(   20)       nimpptl, njmpptl,   &  !: i-, j-indexes for each processor
(   21)       nlcitl , nlcjtl ,   &  !: dimensions of every subdomain
(   22)       ibonitl, ibonjtl
(   23) 
(   24)    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   &  !:
(   25)       tra_ctl                   !: previous trend values
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC prt_ctl_trc         ! called by all subroutines
(   29)    PUBLIC prt_ctl_trc_info    !
(   30)    PUBLIC prt_ctl_trc_init    ! called by opa.F90
(   31)    !!----------------------------------------------------------------------
(   32)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/prtctl_trc.F90,v 1.2 2006/03/21 15:53:52 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36) 
(   37) 
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE prt_ctl_trc (tab4d, mask, clinfo, ovlap, kdim, clinfo2)
(   41)       !!----------------------------------------------------------------------
(   42)       !!                     ***  ROUTINE prt_ctl  ***
(   43)       !!
(   44)       !! ** Purpose : - print sum control 3D arrays over the same area 
(   45)       !!                in mono and mpp case. This way can be usefull when
(   46)       !!                debugging a new parametrization in mono or mpp. 






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 2

(   47)       !!
(   48)       !! ** Method  : 2 possibilities exist when setting the ln_ctl parameter to
(   49)       !!                .true. in the ocean namelist:
(   50)       !!              - to debug a MPI run .vs. a mono-processor one; 
(   51)       !!                the control print will be done over each sub-domain.
(   52)       !!                The nictl[se] and njctl[se] parameters in the namelist must 
(   53)       !!                be set to zero and [ij]splt to the corresponding splitted
(   54)       !!                domain in MPI along respectively i-, j- directions.
(   55)       !!              - to debug a mono-processor run over the whole domain/a specific area; 
(   56)       !!                in the first case the nictl[se] and njctl[se] parameters must be set
(   57)       !!                to zero else to the indices of the area to be controled. In both cases
(   58)       !!                isplt and jsplt must be set to 1.
(   59)       !!              - All arguments of the above calling sequence are optional so their
(   60)       !!                name must be explicitly typed if used. For instance if the mask
(   61)       !!                array tmask(:,:,:) must be passed through the prt_ctl subroutine, 
(   62)       !!                it must looks like: CALL prt_ctl(mask=tmask).
(   63)       !!
(   64)       !!                    tab4d   : 4D array
(   65)       !!                    mask    : mask (3D) to apply to the tab4d array
(   66)       !!                    clinfo  : information about the tab3d array
(   67)       !!                    ovlap   : overlap value
(   68)       !!                    kdim    : k- direction for 4D arrays 
(   69)       !!
(   70)       !! History :
(   71)       !!   9.0  !  05-07  (C. Talandier) original code
(   72)       !!        !  05-10  (C. Ethe     ) adapted to passive tracer
(   73)       !!----------------------------------------------------------------------
(   74)       !! * Arguments
(   75)       REAL(wp), DIMENSION(:,:,:,:), INTENT(in), OPTIONAL :: tab4d
(   76)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
(   77)       CHARACTER (len=*), DIMENSION(:), INTENT(in), OPTIONAL :: clinfo
(   78)       CHARACTER (len=*), INTENT(in), OPTIONAL :: clinfo2
(   79)       INTEGER, INTENT(in), OPTIONAL :: ovlap
(   80)       INTEGER, INTENT(in), OPTIONAL :: kdim
(   81) 
(   82)       !! * Local declarations
(   83)       INTEGER  :: overlap, numid, jn, js, sind, eind, kdir
(   84)       REAL(wp) :: zsum, zvctl
(   85)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmask, ztab3d
(   86)       CHARACTER (len=20), DIMENSION(jptra) :: cl
(   87)       CHARACTER (len=10) :: cl2
(   88)       !!----------------------------------------------------------------------
(   89) 
(   90)       ! Arrays, scalars initialization 
(   91)       overlap       = 0
(   92)       kdir          = jpkm1
(   93)       zsum          = 0.e0
(   94)       zvctl         = 0.e0
(   95)       cl(:)         = ''
(   96)       cl2           = ''
(   97)       ztab3d(:,:,:) = 0.e0
(   98)       zmask (:,:,:) = 1.e0
(   99) 
(  100)       ! Control of optional arguments
(  101) 
(  102)       IF( PRESENT(ovlap)   )  overlap       = ovlap
(  103)       IF( PRESENT(kdim)    )  kdir          = kdim
(  104)       IF( PRESENT(clinfo ) )  cl(:)         = clinfo(:)






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 3

(  105)       IF( PRESENT(clinfo2) )  cl2           = clinfo2
(  106)       IF( PRESENT(mask)    )  zmask (:,:,:) = mask(:,:,:)
(  107) 
(  108)       IF( lk_mpp )   THEN
(  109)          ! processor number
(  110)          sind = narea
(  111)          eind = narea
(  112)       ELSE
(  113)          ! processors total number
(  114)          sind = 1
(  115)          eind = ijsplt
(  116)       ENDIF
(  117) 
(  118)       ! Loop over each sub-domain, i.e. the total number of processors ijsplt
(  119)       DO js = sind, eind
(  120) 
(  121)          numid = 90 + js
(  122) 
(  123)          ! Set indices for the SUM control
(  124)          IF( .NOT. lsp_area ) THEN
(  125)             IF (lk_mpp )   THEN
(  126)                nictls = MAX( 1, nlditl(js) - overlap )
(  127)                nictle = nleitl(js) + overlap * MIN( 1, nlcitl(js) - nleitl(js)) 
(  128)                njctls = MAX( 1, nldjtl(js) - overlap )
(  129)                njctle = nlejtl(js) + overlap * MIN( 1, nlcjtl(js) - nlejtl(js))
(  130)                ! Do not take into account the bound of the domain
(  131)                IF( ibonitl(js) == -1 .OR. ibonitl(js) == 2 ) nictls = MAX(2, nictls)
(  132)                IF( ibonitl(js) ==  1 .OR. ibonitl(js) == 2 ) nictle = MIN(nictle, nleitl(js) - 1)
(  133)                IF( ibonjtl(js) == -1 .OR. ibonjtl(js) == 2 ) njctls = MAX(2, njctls)
(  134)                IF( ibonjtl(js) ==  1 .OR. ibonjtl(js) == 2 ) njctle = MIN(njctle, nlejtl(js) - 1)
(  135)             ELSE
(  136)                nictls = MAX( 1, nimpptl(js) + nlditl(js) - 1 - overlap )
(  137)                nictle = nimpptl(js) + nleitl(js) - 1 + overlap * MIN( 1, nlcitl(js) - nleitl(js) ) 
(  138)                njctls = MAX( 1, njmpptl(js) + nldjtl(js) - 1 - overlap )
(  139)                njctle = njmpptl(js) + nlejtl(js) - 1 + overlap * MIN( 1, nlcjtl(js) - nlejtl(js) ) 
(  140)                ! Do not take into account the bound of the domain
(  141)                IF( ibonitl(js) == -1 .OR. ibonitl(js) == 2 ) nictls = MAX(2, nictls)
(  142)                IF( ibonjtl(js) == -1 .OR. ibonjtl(js) == 2 ) njctls = MAX(2, njctls)
(  143)                IF( ibonitl(js) ==  1 .OR. ibonitl(js) == 2 ) nictle = MIN(nictle, nimpptl(js) + nleitl(js) - 2)
(  144)                IF( ibonjtl(js) ==  1 .OR. ibonjtl(js) == 2 ) njctle = MIN(njctle, njmpptl(js) + nlejtl(js) - 2)
(  145)             ENDIF
(  146)          ENDIF
(  147)          
(  148)          IF( PRESENT(clinfo2) ) THEN
(  149)             DO jn = 1, jptra
(  150)                zvctl  = tra_ctl(jn,js)
(  151)                ztab3d(:,:,:) = tab4d(:,:,:,jn)
(  152)                zsum          = SUM( ztab3d(nictls:nictle,njctls:njctle,1:kdir) &
(  153)                   &                 *zmask(nictls:nictle,njctls:njctle,1:kdir) )
(  154)                WRITE(numid,FMT="(3x,a,' : ',D23.16)") cl(jn), zsum-zvctl
(  155)                tra_ctl(jn,js) = zsum
(  156)             ENDDO
(  157)          ELSE
(  158)             DO jn = 1, jptra
(  159)                ztab3d(:,:,:) = tab4d(:,:,:,jn)
(  160)                zsum          = SUM( ztab3d(nictls:nictle,njctls:njctle,1:kdir) &
(  161)                   &               * zmask(nictls:nictle,njctls:njctle,1:kdir) )
(  162)                WRITE(numid,FMT="(3x,a,' : ',D23.16)") cl(jn), zsum






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 4

(  163)             END DO
(  164)          ENDIF
(  165)          
(  166) 
(  167)       ENDDO
(  168) 
(  169)    END SUBROUTINE prt_ctl_trc
(  170) 
(  171)    SUBROUTINE prt_ctl_trc_info (clinfo)
(  172)       !!----------------------------------------------------------------------
(  173)       !!                     ***  ROUTINE prt_ctl_trc_info  ***
(  174)       !!
(  175)       !! ** Purpose : - print information without any computation
(  176)       !!
(  177)       !! ** Action  : - input arguments
(  178)       !!                    clinfo : information to print
(  179)       !!
(  180)       !! History :
(  181)       !!   9.0  !  05-07  (C. Talandier) original code
(  182)       !!----------------------------------------------------------------------
(  183)       !! * Arguments
(  184)       CHARACTER (len=*), INTENT(in) ::   clinfo
(  185) 
(  186)       !! * Local declarations
(  187)       INTEGER ::  numid, js, sind, eind
(  188)       !!----------------------------------------------------------------------
(  189) 
(  190)       IF( lk_mpp )   THEN
(  191)          ! processor number
(  192)          sind = narea
(  193)          eind = narea
(  194)       ELSE
(  195)          ! total number of processors
(  196)          sind = 1
(  197)          eind = ijsplt
(  198)       ENDIF
(  199) 
(  200)       ! Loop over each sub-domain, i.e. number of processors ijsplt
(  201)       DO js = sind, eind
(  202)          numid = 90 + js
(  203)          WRITE(numid,*)clinfo
(  204)       ENDDO
(  205) 
(  206) 
(  207)    END SUBROUTINE prt_ctl_trc_info
(  208) 
(  209)    SUBROUTINE prt_ctl_trc_init
(  210)       !!----------------------------------------------------------------------
(  211)       !!                     ***  ROUTINE prt_ctl_trc_init  ***
(  212)       !!
(  213)       !! ** Purpose :   open ASCII files & compute indices
(  214)       !!
(  215)       !! History :
(  216)       !!   9.0  !  05-07  (C. Talandier) original code
(  217)       !!        !  05-10  (C. Ethe     ) adapted to passive tracer
(  218)       !!----------------------------------------------------------------------
(  219)       !! * Local declarations
(  220)       INTEGER ::   js, numid, sind, eind






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 5

(  221)       CHARACTER (len=31) :: clfile_out
(  222)       CHARACTER (len=27) :: clb_name
(  223)       CHARACTER (len=19) :: cl_run
(  224)       !!----------------------------------------------------------------------
(  225) 
(  226)       ! Allocate arrays
(  227)       ALLOCATE(nlditl (ijsplt))
(  228)       ALLOCATE(nldjtl (ijsplt))
(  229)       ALLOCATE(nleitl (ijsplt))
(  230)       ALLOCATE(nlejtl (ijsplt))
(  231)       ALLOCATE(nimpptl(ijsplt))
(  232)       ALLOCATE(njmpptl(ijsplt))
(  233)       ALLOCATE(nlcitl (ijsplt))
(  234)       ALLOCATE(nlcjtl (ijsplt))
(  235)       ALLOCATE(tra_ctl(jptra,ijsplt))
(  236)       ALLOCATE(ibonitl(ijsplt))
(  237)       ALLOCATE(ibonjtl(ijsplt))
(  238) 
(  239)       ! Initialization 
(  240)       tra_ctl (:,:)=0.e0
(  241) 
(  242)       IF( lk_mpp ) THEN
(  243)          sind = narea
(  244)          eind = narea
(  245)          clb_name = "('mpp.top.output_',I3.3)"
(  246)          cl_run = 'MULTI processor run'
(  247)          ! use indices for each area computed by mpp_init subroutine
(  248)          nlditl(:) = nldit(:) 
(  249)          nleitl(:) = nleit(:) 
(  250)          nldjtl(:) = nldjt(:) 
(  251)          nlejtl(:) = nlejt(:) 
(  252)          !
(  253)          nimpptl(:) = nimppt(:)
(  254)          njmpptl(:) = njmppt(:)
(  255)          !
(  256)          nlcitl(:) = nlcit(:)
(  257)          nlcjtl(:) = nlcjt(:)
(  258)          !
(  259)          ibonitl(:) = ibonit(:)
(  260)          ibonjtl(:) = ibonjt(:)
(  261)       ELSE
(  262)          sind = 1
(  263)          eind = ijsplt
(  264)          clb_name = "('mono.top.output_',I3.3)"
(  265)          cl_run = 'MONO processor run '
(  266)          ! compute indices for each area as done in mpp_init subroutine
(  267)          CALL sub_dom
(  268)       ENDIF
(  269) 
(  270)       DO js = sind, eind
(  271)          numid = 90 + js
(  272)          WRITE(clfile_out,FMT=clb_name) js-1
(  273)          OPEN ( UNIT=numid, FILE=TRIM(clfile_out),FORM='FORMATTED' )
(  274)          WRITE(numid,*)
(  275)          WRITE(numid,*) '                 L O D Y C - I P S L'
(  276)          WRITE(numid,*) '                     O P A model'
(  277)          WRITE(numid,*) '            Ocean General Circulation Model'
(  278)          WRITE(numid,*) '               version OPA 9.0  (2005) '






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 6

(  279)          WRITE(numid,*)
(  280)          WRITE(numid,*) '                   PROC number: ', js
(  281)          WRITE(numid,*)
(  282)          WRITE(numid,FMT="(19x,a20)")cl_run
(  283) 
(  284)          ! Print the SUM control indices
(  285)          IF( .NOT. lsp_area )   THEN
(  286)             IF ( lk_mpp )   THEN
(  287)                nictls = nlditl(js) 
(  288)                nictle = nleitl(js)
(  289)                njctls = nldjtl(js)
(  290)                njctle = nlejtl(js)
(  291)             ELSE
(  292)                nictls = nimpptl(js) + nlditl(js) - 1
(  293)                nictle = nimpptl(js) + nleitl(js) - 1
(  294)                njctls = njmpptl(js) + nldjtl(js) - 1
(  295)                njctle = njmpptl(js) + nlejtl(js) - 1
(  296)             ENDIF
(  297)          ENDIF
(  298)          WRITE(numid,*) 
(  299)          WRITE(numid,*) 'prt_tra_ctl :  Sum control indices'
(  300)          WRITE(numid,*) '~~~~~~~'
(  301)          WRITE(numid,*)
(  302)          WRITE(numid,9000)'                                nlej   = ', nlejtl(js), '              '
(  303)          WRITE(numid,9000)'                  ------------- njctle = ', njctle, ' -------------'
(  304)          WRITE(numid,9001)'                  |                                       |'
(  305)          WRITE(numid,9001)'                  |                                       |'
(  306)          WRITE(numid,9001)'                  |                                       |'
(  307)          WRITE(numid,9002)'           nictls = ', nictls,  '                           nictle = ', nictle
(  308)          WRITE(numid,9002)'           nldi   = ', nlditl(js),  '                           nlei   = ', nleitl(js)
(  309)          WRITE(numid,9001)'                  |                                       |'
(  310)          WRITE(numid,9001)'                  |                                       |'
(  311)          WRITE(numid,9001)'                  |                                       |'
(  312)          WRITE(numid,9004)'  njmpp  = ',njmpptl(js),'   ------------- njctls = ', njctls, ' -------------'
(  313)          WRITE(numid,9003)'           nimpp  = ', nimpptl(js), '        nldj   = ', nldjtl(js), '              '
(  314)          WRITE(numid,*)
(  315)          WRITE(numid,*)
(  316) 
(  317) 9000     FORMAT(a41,i4.4,a14)
(  318) 9001     FORMAT(a59)
(  319) 9002     FORMAT(a20,i4.4,a36,i3.3)
(  320) 9003     FORMAT(a20,i4.4,a17,i4.4)
(  321) 9004     FORMAT(a11,i4.4,a26,i4.4,a14)
(  322)       ENDDO
(  323) 
(  324)    END SUBROUTINE prt_ctl_trc_init
(  325) 
(  326) 
(  327)    SUBROUTINE sub_dom
(  328)       !!----------------------------------------------------------------------
(  329)       !!                  ***  ROUTINE sub_dom  ***
(  330)       !!                    
(  331)       !! ** Purpose :   Lay out the global domain over processors. 
(  332)       !!                CAUTION: 
(  333)       !!                This part has been extracted from the mpp_init
(  334)       !!                subroutine and names of variables/arrays have been 
(  335)       !!                slightly changed to avoid confusion but the computation
(  336)       !!                is exactly the same. Any modification about indices of






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 7

(  337)       !!                each sub-domain in the mppini.F90 module should be reported 
(  338)       !!                here.
(  339)       !!
(  340)       !! ** Method  :   Global domain is distributed in smaller local domains.
(  341)       !!                Periodic condition is a function of the local domain position
(  342)       !!                (global boundary or neighbouring domain) and of the global
(  343)       !!                periodic
(  344)       !!                Type :         jperio global periodic condition
(  345)       !!                               nperio local  periodic condition
(  346)       !!
(  347)       !! ** Action  : - set domain parameters
(  348)       !!                    nimpp     : longitudinal index 
(  349)       !!                    njmpp     : latitudinal  index
(  350)       !!                    nperio    : lateral condition type 
(  351)       !!                    narea     : number for local area
(  352)       !!                    nlcil      : first dimension
(  353)       !!                    nlcjl      : second dimension
(  354)       !!                    nbondil    : mark for "east-west local boundary"
(  355)       !!                    nbondjl    : mark for "north-south local boundary"
(  356)       !!
(  357)       !! History :
(  358)       !!        !  94-11  (M. Guyon)  Original code
(  359)       !!        !  95-04  (J. Escobar, M. Imbard)
(  360)       !!        !  98-02  (M. Guyon)  FETI method
(  361)       !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
(  362)       !!   8.5  !  02-08  (G. Madec)  F90 : free form
(  363)       !!----------------------------------------------------------------------
(  364)       !! * Local variables
(  365)       INTEGER ::   ji, jj, js               ! dummy loop indices
(  366)       INTEGER ::   &
(  367)          ii, ij,                         &  ! temporary integers
(  368)          irestil, irestjl,               &  !    "          "
(  369)          ijpi  , ijpj, nlcil,            &  ! temporary logical unit
(  370)          nlcjl , nbondil, nbondjl,       &
(  371)          nrecil, nrecjl, nldil, nleil, nldjl, nlejl
(  372) 
(  373)       INTEGER, DIMENSION(:,:), ALLOCATABLE ::   &
(  374)          iimpptl, ijmpptl, ilcitl, ilcjtl       ! temporary workspace
(  375)       REAL(wp) ::   zidom, zjdom            ! temporary scalars
(  376)       !!----------------------------------------------------------------------
(  377) 
(  378)       !  1. Dimension arrays for subdomains
(  379)       ! -----------------------------------
(  380)       !  Computation of local domain sizes ilcitl() ilcjtl()
(  381)       !  These dimensions depend on global sizes isplt,jsplt and jpiglo,jpjglo
(  382)       !  The subdomains are squares leeser than or equal to the global
(  383)       !  dimensions divided by the number of processors minus the overlap
(  384)       !  array (cf. par_oce.F90).
(  385) 
(  386)       ijpi = ( jpiglo-2*jpreci + (isplt-1) ) / isplt + 2*jpreci
(  387)       ijpj = ( jpjglo-2*jprecj + (jsplt-1) ) / jsplt + 2*jprecj
(  388) 
(  389)       ALLOCATE(ilcitl (isplt,jsplt))
(  390)       ALLOCATE(ilcjtl (isplt,jsplt))
(  391) 
(  392)       nrecil  = 2 * jpreci
(  393)       nrecjl  = 2 * jprecj
(  394)       irestil = MOD( jpiglo - nrecil , isplt )






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 8

(  395)       irestjl = MOD( jpjglo - nrecjl , jsplt )
(  396) 
(  397)       IF(  irestil == 0 )   irestil = isplt
(  398)       DO jj = 1, jsplt
(  399)          DO ji = 1, irestil
(  400)             ilcitl(ji,jj) = ijpi
(  401)          END DO
(  402)          DO ji = irestil+1, isplt
(  403)             ilcitl(ji,jj) = ijpi -1
(  404)          END DO
(  405)       END DO
(  406)       
(  407)       IF( irestjl == 0 )   irestjl = jsplt
(  408)       DO ji = 1, isplt
(  409)          DO jj = 1, irestjl
(  410)             ilcjtl(ji,jj) = ijpj
(  411)          END DO
(  412)          DO jj = irestjl+1, jsplt
(  413)             ilcjtl(ji,jj) = ijpj -1
(  414)          END DO
(  415)       END DO
(  416)       
(  417)       zidom = nrecil
(  418)       DO ji = 1, isplt
(  419)          zidom = zidom + ilcitl(ji,1) - nrecil
(  420)       END DO
(  421)       
(  422)       zjdom = nrecjl
(  423)       DO jj = 1, jsplt
(  424)          zjdom = zjdom + ilcjtl(1,jj) - nrecjl
(  425)       END DO
(  426) 
(  427)       !  2. Index arrays for subdomains
(  428)       ! -------------------------------
(  429) 
(  430)       ALLOCATE(iimpptl(isplt,jsplt))
(  431)       ALLOCATE(ijmpptl(isplt,jsplt))
(  432)       
(  433)       iimpptl(:,:) = 1
(  434)       ijmpptl(:,:) = 1
(  435)       
(  436)       IF( isplt > 1 ) THEN
(  437)          DO jj = 1, jsplt
(  438)             DO ji = 2, isplt
(  439)                iimpptl(ji,jj) = iimpptl(ji-1,jj) + ilcitl(ji-1,jj) - nrecil
(  440)             END DO
(  441)          END DO
(  442)       ENDIF
(  443) 
(  444)       IF( jsplt > 1 ) THEN
(  445)          DO jj = 2, jsplt
(  446)             DO ji = 1, isplt
(  447)                ijmpptl(ji,jj) = ijmpptl(ji,jj-1)+ilcjtl(ji,jj-1)-nrecjl
(  448)             END DO
(  449)          END DO
(  450)       ENDIF
(  451)       
(  452)       ! 3. Subdomain description






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 9

(  453)       ! ------------------------
(  454) 
(  455)       DO js = 1, ijsplt
(  456)          ii = 1 + MOD( js-1, isplt )
(  457)          ij = 1 + (js-1) / isplt
(  458)          nimpptl(js) = iimpptl(ii,ij)
(  459)          njmpptl(js) = ijmpptl(ii,ij)
(  460)          nlcitl (js) = ilcitl (ii,ij)     
(  461)          nlcil       = nlcitl (js)     
(  462)          nlcjtl (js) = ilcjtl (ii,ij)     
(  463)          nlcjl       = nlcjtl (js)
(  464)          nbondjl = -1                                    ! general case
(  465)          IF( js   >  isplt          )   nbondjl = 0      ! first row of processor
(  466)          IF( js   >  (jsplt-1)*isplt )  nbondjl = 1     ! last  row of processor
(  467)          IF( jsplt == 1             )   nbondjl = 2      ! one processor only in j-direction
(  468)          ibonjtl(js) = nbondjl
(  469)          
(  470)          nbondil = 0                                     ! 
(  471)          IF( MOD( js, isplt ) == 1 )   nbondil = -1      !
(  472)          IF( MOD( js, isplt ) == 0 )   nbondil =  1      !
(  473)          IF( isplt            == 1 )   nbondil =  2      ! one processor only in i-direction
(  474)          ibonitl(js) = nbondil
(  475)          
(  476)          nldil =  1   + jpreci
(  477)          nleil = nlcil - jpreci
(  478)          IF( nbondil == -1 .OR. nbondil == 2 )   nldil = 1
(  479)          IF( nbondil ==  1 .OR. nbondil == 2 )   nleil = nlcil
(  480)          nldjl =  1   + jprecj
(  481)          nlejl = nlcjl - jprecj
(  482)          IF( nbondjl == -1 .OR. nbondjl == 2 )   nldjl = 1
(  483)          IF( nbondjl ==  1 .OR. nbondjl == 2 )   nlejl = nlcjl
(  484)          nlditl(js) = nldil
(  485)          nleitl(js) = nleil
(  486)          nldjtl(js) = nldjl
(  487)          nlejtl(js) = nlejl
(  488)       END DO
(  489) 
(  490)       DEALLOCATE(iimpptl)
(  491)       DEALLOCATE(ijmpptl)
(  492)       DEALLOCATE(ilcitl)
(  493)       DEALLOCATE(ilcjtl)
(  494) 
(  495)    END SUBROUTINE sub_dom
(  496)  
(  497) #else
(  498)    !!----------------------------------------------------------------------
(  499)    !!   Dummy module :                      NO passive tracer
(  500)    !!----------------------------------------------------------------------
(  501) #endif
(  502)     
(  503)    !!======================================================================
(  504) 
(  505) END MODULE prtctl_trc











PGF90 (Version     10.5)          09/23/2010  09:49:17      page 10

(    1) # 1 "prtctl_trc.F90"
(    1) MODULE prtctl_trc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE prtctl   ***
(    4)    !! Ocean system   : print all SUM trends for each processor domain
(    5)    !!==============================================================================
(    6) # 498
(  498)    !!----------------------------------------------------------------------
(  499)    !!   Dummy module :                      NO passive tracer
(  500)    !!----------------------------------------------------------------------
(  501) # 502
(  502)     
(  503)    !!======================================================================
(  504) 
(  505) END MODULE prtctl_trc
