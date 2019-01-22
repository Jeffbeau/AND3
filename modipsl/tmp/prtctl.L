


PGF90 (Version     10.5)          09/23/2010  09:49:00      page 1

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

Filename: prtctl.F90

(    1) MODULE prtctl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE prtctl   ***
(    4)    !! Ocean system   : print all SUM trends for each processor domain
(    5)    !!==============================================================================
(    6)    USE dom_oce          ! ocean space and time domain variables
(    7)    USE in_out_manager   ! I/O manager
(    8)    USE lib_mpp          ! distributed memory computing
(    9) 
(   10)    IMPLICIT NONE
(   11)    PRIVATE
(   12) 
(   13)    !! * Module declaration
(   14)    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::   &  !:
(   15)       nlditl , nldjtl ,   &  !: first, last indoor index for each i-domain
(   16)       nleitl , nlejtl ,   &  !: first, last indoor index for each j-domain
(   17)       nimpptl, njmpptl,   &  !: i-, j-indexes for each processor
(   18)       nlcitl , nlcjtl ,   &  !: dimensions of every subdomain
(   19)       ibonitl, ibonjtl
(   20) 
(   21)    REAL(wp), DIMENSION(:), ALLOCATABLE ::   &  !:
(   22)       t_ctll , s_ctll ,   &  !: previous trend values
(   23)       u_ctll , v_ctll
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC prt_ctl         ! called by all subroutines
(   27)    PUBLIC prt_ctl_info    ! called by all subroutines
(   28)    PUBLIC prt_ctl_init    ! called by opa.F90
(   29)    !!----------------------------------------------------------------------
(   30)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/prtctl.F90,v 1.2 2006/03/21 15:54:15 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35) 
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE prt_ctl (tab2d_1, tab3d_1, mask1, clinfo1, tab2d_2, tab3d_2, mask2, clinfo2, ovlap, kdim, clinfo3)
(   39)       !!----------------------------------------------------------------------
(   40)       !!                     ***  ROUTINE prt_ctl  ***
(   41)       !!
(   42)       !! ** Purpose : - print sum control of 2D or 3D arrays over the same area 
(   43)       !!                in mono and mpp case. This way can be usefull when
(   44)       !!                debugging a new parametrization in mono or mpp. 
(   45)       !!
(   46)       !! ** Method  : 2 possibilities exist when setting the ln_ctl parameter to






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 2

(   47)       !!                .true. in the ocean namelist:
(   48)       !!              - to debug a MPI run .vs. a mono-processor one; 
(   49)       !!                the control print will be done over each sub-domain.
(   50)       !!                The nictl[se] and njctl[se] parameters in the namelist must 
(   51)       !!                be set to zero and [ij]splt to the corresponding splitted
(   52)       !!                domain in MPI along respectively i-, j- directions.
(   53)       !!              - to debug a mono-processor run over the whole domain/a specific area; 
(   54)       !!                in the first case the nictl[se] and njctl[se] parameters must be set
(   55)       !!                to zero else to the indices of the area to be controled. In both cases
(   56)       !!                isplt and jsplt must be set to 1.
(   57)       !!              - All arguments of the above calling sequence are optional so their
(   58)       !!                name must be explicitly typed if used. For instance if the 3D
(   59)       !!                array tn(:,:,:) must be passed through the prt_ctl subroutine, 
(   60)       !!                it must looks like: CALL prt_ctl(tab3d_1=tn).
(   61)       !!
(   62)       !!                    tab2d_1 : first 2D array
(   63)       !!                    tab3d_1 : first 3D array
(   64)       !!                    mask1   : mask (3D) to apply to the tab[23]d_1 array
(   65)       !!                    clinfo1 : information about the tab[23]d_1 array
(   66)       !!                    tab2d_2 : second 2D array
(   67)       !!                    tab3d_2 : second 3D array
(   68)       !!                    mask2   : mask (3D) to apply to the tab[23]d_2 array
(   69)       !!                    clinfo2 : information about the tab[23]d_2 array
(   70)       !!                    ovlap   : overlap value
(   71)       !!                    kdim    : k- direction for 3D arrays 
(   72)       !!                    clinfo3 : additional information 
(   73)       !!
(   74)       !! History :
(   75)       !!   9.0  !  05-07  (C. Talandier) original code
(   76)       !!----------------------------------------------------------------------
(   77)       !! * Arguments
(   78)       REAL(wp), DIMENSION(:,:)  , INTENT(in), OPTIONAL :: tab2d_1
(   79)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: tab3d_1
(   80)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask1
(   81)       CHARACTER (len=*), INTENT(in), OPTIONAL :: clinfo1
(   82)       REAL(wp), DIMENSION(:,:)  , INTENT(in), OPTIONAL :: tab2d_2
(   83)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: tab3d_2
(   84)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask2
(   85)       CHARACTER (len=*), INTENT(in), OPTIONAL :: clinfo2
(   86)       INTEGER, INTENT(in), OPTIONAL :: ovlap
(   87)       INTEGER, INTENT(in), OPTIONAL :: kdim
(   88)       CHARACTER (len=*), INTENT(in), OPTIONAL :: clinfo3
(   89) 
(   90)       !! * Local declarations
(   91)       INTEGER :: overlap, numid, jn, sind, eind, kdir
(   92)       CHARACTER (len=15) :: cl2
(   93)       REAL(wp) :: zsum1, zsum2, zvctl1, zvctl2
(   94)       REAL(wp), DIMENSION(jpi,jpj) :: ztab2d_1, ztab2d_2
(   95)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmask1, zmask2, ztab3d_1, ztab3d_2
(   96)       !!----------------------------------------------------------------------
(   97) 
(   98)       ! Arrays, scalars initialization 
(   99)       overlap   = 0
(  100)       kdir      = jpkm1
(  101)       cl2       = ''
(  102)       zsum1     = 0.e0
(  103)       zsum2     = 0.e0
(  104)       zvctl1    = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 3

(  105)       zvctl2    = 0.e0
(  106)       ztab2d_1(:,:)   = 0.e0
(  107)       ztab2d_2(:,:)   = 0.e0
(  108)       ztab3d_1(:,:,:) = 0.e0
(  109)       ztab3d_2(:,:,:) = 0.e0
(  110)       zmask1  (:,:,:) = 1.e0
(  111)       zmask2  (:,:,:) = 1.e0
(  112) 
(  113)       ! Control of optional arguments
(  114)       IF( PRESENT(clinfo2) )  cl2            = clinfo2
(  115)       IF( PRESENT(ovlap)   )  overlap        = ovlap
(  116)       IF( PRESENT(kdim)    )  kdir           = kdim
(  117)       IF( PRESENT(tab2d_1) )  ztab2d_1(:,:)  = tab2d_1(:,:)
(  118)       IF( PRESENT(tab2d_2) )  ztab2d_2(:,:)  = tab2d_2(:,:)
(  119)       IF( PRESENT(tab3d_1) )  ztab3d_1(:,:,:)= tab3d_1(:,:,:)
(  120)       IF( PRESENT(tab3d_2) )  ztab3d_2(:,:,:)= tab3d_2(:,:,:)
(  121)       IF( PRESENT(mask1)   )  zmask1  (:,:,:)= mask1  (:,:,:)
(  122)       IF( PRESENT(mask2)   )  zmask2  (:,:,:)= mask2  (:,:,:)
(  123) 
(  124)       IF( lk_mpp )   THEN
(  125)          ! processor number
(  126)          sind = narea
(  127)          eind = narea
(  128)       ELSE
(  129)          ! processors total number
(  130)          sind = 1
(  131)          eind = ijsplt
(  132)       ENDIF
(  133) 
(  134)       ! Loop over each sub-domain, i.e. the total number of processors ijsplt
(  135)       DO jn = sind, eind
(  136) 
(  137)          numid = 80 + jn
(  138) 
(  139)          ! Set indices for the SUM control
(  140)          IF( .NOT. lsp_area ) THEN
(  141)             IF (lk_mpp )   THEN
(  142)                nictls = MAX( 1, nlditl(jn) - overlap )
(  143)                nictle = nleitl(jn) + overlap * MIN( 1, nlcitl(jn) - nleitl(jn)) 
(  144)                njctls = MAX( 1, nldjtl(jn) - overlap )
(  145)                njctle = nlejtl(jn) + overlap * MIN( 1, nlcjtl(jn) - nlejtl(jn))
(  146)                ! Do not take into account the bound of the domain
(  147)                IF( ibonitl(jn) == -1 .OR. ibonitl(jn) == 2 ) nictls = MAX(2, nictls)
(  148)                IF( ibonjtl(jn) == -1 .OR. ibonjtl(jn) == 2 ) njctls = MAX(2, njctls)
(  149)                IF( ibonitl(jn) ==  1 .OR. ibonitl(jn) == 2 ) nictle = MIN(nictle, nleitl(jn) - 1)
(  150)                IF( ibonjtl(jn) ==  1 .OR. ibonjtl(jn) == 2 ) njctle = MIN(njctle, nlejtl(jn) - 1)
(  151)             ELSE
(  152)                nictls = MAX( 1, nimpptl(jn) + nlditl(jn) - 1 - overlap )
(  153)                nictle = nimpptl(jn) + nleitl(jn) - 1 + overlap * MIN( 1, nlcitl(jn) - nleitl(jn) ) 
(  154)                njctls = MAX( 1, njmpptl(jn) + nldjtl(jn) - 1 - overlap )
(  155)                njctle = njmpptl(jn) + nlejtl(jn) - 1 + overlap * MIN( 1, nlcjtl(jn) - nlejtl(jn) ) 
(  156)                ! Do not take into account the bound of the domain
(  157)                IF( ibonitl(jn) == -1 .OR. ibonitl(jn) == 2 ) nictls = MAX(2, nictls)
(  158)                IF( ibonjtl(jn) == -1 .OR. ibonjtl(jn) == 2 ) njctls = MAX(2, njctls)
(  159)                IF( ibonitl(jn) ==  1 .OR. ibonitl(jn) == 2 ) nictle = MIN(nictle, nimpptl(jn) + nleitl(jn) - 2)
(  160)                IF( ibonjtl(jn) ==  1 .OR. ibonjtl(jn) == 2 ) njctle = MIN(njctle, njmpptl(jn) + nlejtl(jn) - 2)
(  161)             ENDIF
(  162)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 4

(  163) 
(  164)          IF ( clinfo3 == 'tra' )  THEN
(  165)              zvctl1 = t_ctll(jn)
(  166)              zvctl2 = s_ctll(jn)
(  167)          ELSEIF ( clinfo3 == 'dyn' )   THEN
(  168)              zvctl1 = u_ctll(jn)
(  169)              zvctl2 = v_ctll(jn)
(  170)          ENDIF
(  171) 
(  172)          ! Compute the sum control
(  173)          ! 2D arrays
(  174)          IF( PRESENT(tab2d_1) )   THEN
(  175)             zsum1 = SUM( ztab2d_1(nictls:nictle,njctls:njctle)*zmask1(nictls:nictle,njctls:njctle,1) )
(  176)             zsum2 = SUM( ztab2d_2(nictls:nictle,njctls:njctle)*zmask2(nictls:nictle,njctls:njctle,1) )
(  177)          ENDIF
(  178) 
(  179)          ! 3D arrays
(  180)          IF( PRESENT(tab3d_1) )   THEN
(  181)             zsum1 = SUM( ztab3d_1(nictls:nictle,njctls:njctle,1:kdir)*zmask1(nictls:nictle,njctls:njctle,1:kdir) )
(  182)             zsum2 = SUM( ztab3d_2(nictls:nictle,njctls:njctle,1:kdir)*zmask2(nictls:nictle,njctls:njctle,1:kdir) )
(  183)          ENDIF
(  184) 
(  185)          ! Print the result
(  186)          IF( PRESENT(clinfo2) .OR. PRESENT(tab2d_2) .OR. PRESENT(tab3d_2) )   THEN
(  187)             WRITE(numid,FMT='(a,D23.16,3x,a,D23.16)')clinfo1, zsum1-zvctl1, cl2, zsum2-zvctl2
(  188)              IF( clinfo3 == 'tra-ta' )  THEN
(  189)                 t_ctll(jn) = zsum1
(  190)              ENDIF 
(  191)              IF( clinfo3 == 'tra' )  THEN
(  192)                 t_ctll(jn) = zsum1
(  193)                 s_ctll(jn) = zsum2
(  194)              ENDIF 
(  195)              IF( clinfo3 == 'dyn' )   THEN
(  196)                 u_ctll(jn) = zsum1
(  197)                 v_ctll(jn) = zsum2 
(  198)              ENDIF
(  199)          ELSE
(  200)             WRITE(numid,FMT='(a,D23.16)')clinfo1, zsum1
(  201)          ENDIF
(  202) 
(  203)       ENDDO
(  204) 
(  205)    END SUBROUTINE prt_ctl
(  206) 
(  207) 
(  208)    SUBROUTINE prt_ctl_info (clinfo1, ivar1, clinfo2, ivar2)
(  209)       !!----------------------------------------------------------------------
(  210)       !!                     ***  ROUTINE prt_ctl_info  ***
(  211)       !!
(  212)       !! ** Purpose : - print information without any computation
(  213)       !!
(  214)       !! ** Action  : - input arguments
(  215)       !!                    clinfo1 : information about the ivar1
(  216)       !!                    ivar1   : value to print
(  217)       !!                    clinfo2 : information about the ivar2
(  218)       !!                    ivar2   : value to print
(  219)       !!
(  220)       !! History :






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 5

(  221)       !!   9.0  !  05-07  (C. Talandier) original code
(  222)       !!----------------------------------------------------------------------
(  223)       !! * Arguments
(  224)       CHARACTER (len=*), INTENT(in) ::   clinfo1
(  225)       INTEGER, INTENT(in), OPTIONAL ::   ivar1
(  226)       CHARACTER (len=*), INTENT(in), OPTIONAL ::   clinfo2
(  227)       INTEGER, INTENT(in), OPTIONAL ::   ivar2
(  228) 
(  229)       !! * Local declarations
(  230)       INTEGER ::  numid, jn, sind, eind
(  231)       !!----------------------------------------------------------------------
(  232) 
(  233)       IF( lk_mpp )   THEN
(  234)          ! processor number
(  235)          sind = narea
(  236)          eind = narea
(  237)       ELSE
(  238)          ! total number of processors
(  239)          sind = 1
(  240)          eind = ijsplt
(  241)       ENDIF
(  242) 
(  243)       ! Loop over each sub-domain, i.e. number of processors ijsplt
(  244)       DO jn = sind, eind
(  245) 
(  246)          numid = 80 + jn
(  247) 
(  248)          IF( PRESENT(ivar1) .AND. PRESENT(clinfo2) .AND. PRESENT(ivar2) )   THEN
(  249)             WRITE(numid,*)clinfo1, ivar1, clinfo2, ivar2
(  250)          ELSEIF ( PRESENT(ivar1) .AND. PRESENT(clinfo2) .AND. .NOT. PRESENT(ivar2) )   THEN
(  251)             WRITE(numid,*)clinfo1, ivar1, clinfo2
(  252)          ELSEIF ( PRESENT(ivar1) .AND. .NOT. PRESENT(clinfo2) .AND. PRESENT(ivar2) )   THEN
(  253)             WRITE(numid,*)clinfo1, ivar1, ivar2
(  254)          ELSEIF ( PRESENT(ivar1) .AND. .NOT. PRESENT(clinfo2) .AND. .NOT. PRESENT(ivar2) )   THEN
(  255)             WRITE(numid,*)clinfo1, ivar1
(  256)          ELSE
(  257)             WRITE(numid,*)clinfo1
(  258)          ENDIF
(  259) 
(  260)       ENDDO
(  261) 
(  262) 
(  263)       END SUBROUTINE prt_ctl_info
(  264) 
(  265)    SUBROUTINE prt_ctl_init
(  266)       !!----------------------------------------------------------------------
(  267)       !!                     ***  ROUTINE prt_ctl_init  ***
(  268)       !!
(  269)       !! ** Purpose :   open ASCII files & compute indices
(  270)       !!
(  271)       !! History :
(  272)       !!   9.0  !  05-07  (C. Talandier) original code
(  273)       !!----------------------------------------------------------------------
(  274)       !! * Local declarations
(  275)       INTEGER ::   jn, numid, sind, eind
(  276)       CHARACTER (len=28) :: clfile_out
(  277)       CHARACTER (len=23) :: clb_name
(  278)       CHARACTER (len=19) :: cl_run






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 6

(  279)       !!----------------------------------------------------------------------
(  280) 
(  281)       ! Allocate arrays
(  282)       ALLOCATE(nlditl (ijsplt))
(  283)       ALLOCATE(nldjtl (ijsplt))
(  284)       ALLOCATE(nleitl (ijsplt))
(  285)       ALLOCATE(nlejtl (ijsplt))
(  286)       ALLOCATE(nimpptl(ijsplt))
(  287)       ALLOCATE(njmpptl(ijsplt))
(  288)       ALLOCATE(nlcitl (ijsplt))
(  289)       ALLOCATE(nlcjtl (ijsplt))
(  290)       ALLOCATE(t_ctll (ijsplt))
(  291)       ALLOCATE(s_ctll (ijsplt))
(  292)       ALLOCATE(u_ctll (ijsplt))
(  293)       ALLOCATE(v_ctll (ijsplt))
(  294)       ALLOCATE(ibonitl(ijsplt))
(  295)       ALLOCATE(ibonjtl(ijsplt))
(  296) 
(  297)       ! Initialization 
(  298)       t_ctll (:)=0.e0
(  299)       s_ctll (:)=0.e0
(  300)       u_ctll (:)=0.e0
(  301)       v_ctll (:)=0.e0
(  302) 
(  303)       IF( lk_mpp ) THEN
(  304)          sind = narea
(  305)          eind = narea
(  306)          clb_name = "('mpp.output_',I4.4)"
(  307)          cl_run = 'MULTI processor run'
(  308)          ! use indices for each area computed by mpp_init subroutine
(  309)          nlditl(:) = nldit(:) 
(  310)          nleitl(:) = nleit(:) 
(  311)          nldjtl(:) = nldjt(:) 
(  312)          nlejtl(:) = nlejt(:) 
(  313)          !
(  314)          nimpptl(:) = nimppt(:)
(  315)          njmpptl(:) = njmppt(:)
(  316)          !
(  317)          nlcitl(:) = nlcit(:)
(  318)          nlcjtl(:) = nlcjt(:)
(  319)          !
(  320)          ibonitl(:) = ibonit(:)
(  321)          ibonjtl(:) = ibonjt(:)
(  322)       ELSE
(  323)          sind = 1
(  324)          eind = ijsplt
(  325)          clb_name = "('mono.output_',I4.4)"
(  326)          cl_run = 'MONO processor run '
(  327)          ! compute indices for each area as done in mpp_init subroutine
(  328)          CALL sub_dom
(  329)       ENDIF
(  330) 
(  331)       DO jn = sind, eind
(  332)          numid = 80 + jn
(  333)          WRITE(clfile_out,FMT=clb_name) jn-1
(  334)          OPEN ( UNIT=numid, FILE=TRIM(clfile_out),FORM='FORMATTED' )
(  335)          WRITE(numid,*)
(  336)          WRITE(numid,*) '                 L O D Y C - I P S L'






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 7

(  337)          WRITE(numid,*) '                     O P A model'
(  338)          WRITE(numid,*) '            Ocean General Circulation Model'
(  339)          WRITE(numid,*) '               version OPA 9.0  (2005) '
(  340)          WRITE(numid,*)
(  341)          WRITE(numid,*) '                   PROC number: ', jn
(  342)          WRITE(numid,*)
(  343)          WRITE(numid,FMT="(19x,a20)")cl_run
(  344) 
(  345)          ! Print the SUM control indices
(  346)          IF( .NOT. lsp_area )   THEN
(  347)             IF ( lk_mpp )   THEN
(  348)                nictls = nlditl(jn) 
(  349)                nictle = nleitl(jn)
(  350)                njctls = nldjtl(jn)
(  351)                njctle = nlejtl(jn)
(  352)             ELSE
(  353)                nictls = nimpptl(jn) + nlditl(jn) - 1
(  354)                nictle = nimpptl(jn) + nleitl(jn) - 1
(  355)                njctls = njmpptl(jn) + nldjtl(jn) - 1
(  356)                njctle = njmpptl(jn) + nlejtl(jn) - 1
(  357)             ENDIF
(  358)          ENDIF
(  359)          WRITE(numid,*) 
(  360)          WRITE(numid,*) 'prt_ctl :  Sum control indices'
(  361)          WRITE(numid,*) '~~~~~~~'
(  362)          WRITE(numid,*)
(  363)          WRITE(numid,9000)'                                nlej   = ', nlejtl(jn), '              '
(  364)          WRITE(numid,9000)'                  ------------- njctle = ', njctle, ' -------------'
(  365)          WRITE(numid,9001)'                  |                                       |'
(  366)          WRITE(numid,9001)'                  |                                       |'
(  367)          WRITE(numid,9001)'                  |                                       |'
(  368)          WRITE(numid,9002)'           nictls = ', nictls,  '                           nictle = ', nictle
(  369)          WRITE(numid,9002)'           nldi   = ', nlditl(jn),  '                           nlei   = ', nleitl(jn)
(  370)          WRITE(numid,9001)'                  |                                       |'
(  371)          WRITE(numid,9001)'                  |                                       |'
(  372)          WRITE(numid,9001)'                  |                                       |'
(  373)          WRITE(numid,9004)'  njmpp  = ',njmpptl(jn),'   ------------- njctls = ', njctls, ' -------------'
(  374)          WRITE(numid,9003)'           nimpp  = ', nimpptl(jn), '        nldj   = ', nldjtl(jn), '              '
(  375)          WRITE(numid,*)
(  376)          WRITE(numid,*)
(  377) 
(  378) 9000     FORMAT(a41,i4.4,a14)
(  379) 9001     FORMAT(a59)
(  380) 9002     FORMAT(a20,i4.4,a36,i3.3)
(  381) 9003     FORMAT(a20,i4.4,a17,i4.4)
(  382) 9004     FORMAT(a11,i4.4,a26,i4.4,a14)
(  383)       ENDDO
(  384) 
(  385)    END SUBROUTINE prt_ctl_init
(  386) 
(  387) 
(  388)    SUBROUTINE sub_dom
(  389)       !!----------------------------------------------------------------------
(  390)       !!                  ***  ROUTINE sub_dom  ***
(  391)       !!                    
(  392)       !! ** Purpose :   Lay out the global domain over processors. 
(  393)       !!                CAUTION: 
(  394)       !!                This part has been extracted from the mpp_init






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 8

(  395)       !!                subroutine and names of variables/arrays have been 
(  396)       !!                slightly changed to avoid confusion but the computation
(  397)       !!                is exactly the same. Any modification about indices of
(  398)       !!                each sub-domain in the mppini.F90 module should be reported 
(  399)       !!                here.
(  400)       !!
(  401)       !! ** Method  :   Global domain is distributed in smaller local domains.
(  402)       !!                Periodic condition is a function of the local domain position
(  403)       !!                (global boundary or neighbouring domain) and of the global
(  404)       !!                periodic
(  405)       !!                Type :         jperio global periodic condition
(  406)       !!                               nperio local  periodic condition
(  407)       !!
(  408)       !! ** Action  : - set domain parameters
(  409)       !!                    nimpp     : longitudinal index 
(  410)       !!                    njmpp     : latitudinal  index
(  411)       !!                    nperio    : lateral condition type 
(  412)       !!                    narea     : number for local area
(  413)       !!                    nlcil      : first dimension
(  414)       !!                    nlcjl      : second dimension
(  415)       !!                    nbondil    : mark for "east-west local boundary"
(  416)       !!                    nbondjl    : mark for "north-south local boundary"
(  417)       !!
(  418)       !! History :
(  419)       !!        !  94-11  (M. Guyon)  Original code
(  420)       !!        !  95-04  (J. Escobar, M. Imbard)
(  421)       !!        !  98-02  (M. Guyon)  FETI method
(  422)       !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
(  423)       !!   8.5  !  02-08  (G. Madec)  F90 : free form
(  424)       !!----------------------------------------------------------------------
(  425)       !! * Local variables
(  426)       INTEGER ::   ji, jj, jn               ! dummy loop indices
(  427)       INTEGER ::   &
(  428)          ii, ij,                         &  ! temporary integers
(  429)          irestil, irestjl,               &  !    "          "
(  430)          ijpi  , ijpj, nlcil,            &  ! temporary logical unit
(  431)          nlcjl , nbondil, nbondjl,       &
(  432)          nrecil, nrecjl, nldil, nleil, nldjl, nlejl
(  433) 
(  434)       INTEGER, DIMENSION(:,:), ALLOCATABLE ::   &
(  435)          iimpptl, ijmpptl, ilcitl, ilcjtl       ! temporary workspace
(  436)       REAL(wp) ::   zidom, zjdom            ! temporary scalars
(  437)       !!----------------------------------------------------------------------
(  438) 
(  439)       !  1. Dimension arrays for subdomains
(  440)       ! -----------------------------------
(  441)       !  Computation of local domain sizes ilcitl() ilcjtl()
(  442)       !  These dimensions depend on global sizes isplt,jsplt and jpiglo,jpjglo
(  443)       !  The subdomains are squares leeser than or equal to the global
(  444)       !  dimensions divided by the number of processors minus the overlap
(  445)       !  array (cf. par_oce.F90).
(  446) 
(  447)       ijpi = ( jpiglo-2*jpreci + (isplt-1) ) / isplt + 2*jpreci
(  448)       ijpj = ( jpjglo-2*jprecj + (jsplt-1) ) / jsplt + 2*jprecj
(  449) 
(  450)       ALLOCATE(ilcitl (isplt,jsplt))
(  451)       ALLOCATE(ilcjtl (isplt,jsplt))
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 9

(  453)       nrecil  = 2 * jpreci
(  454)       nrecjl  = 2 * jprecj
(  455)       irestil = MOD( jpiglo - nrecil , isplt )
(  456)       irestjl = MOD( jpjglo - nrecjl , jsplt )
(  457) 
(  458)       IF(  irestil == 0 )   irestil = isplt
(  459)       DO jj = 1, jsplt
(  460)          DO ji = 1, irestil
(  461)             ilcitl(ji,jj) = ijpi
(  462)          END DO
(  463)          DO ji = irestil+1, isplt
(  464)             ilcitl(ji,jj) = ijpi -1
(  465)          END DO
(  466)       END DO
(  467)       
(  468)       IF( irestjl == 0 )   irestjl = jsplt
(  469)       DO ji = 1, isplt
(  470)          DO jj = 1, irestjl
(  471)             ilcjtl(ji,jj) = ijpj
(  472)          END DO
(  473)          DO jj = irestjl+1, jsplt
(  474)             ilcjtl(ji,jj) = ijpj -1
(  475)          END DO
(  476)       END DO
(  477)       
(  478)       zidom = nrecil
(  479)       DO ji = 1, isplt
(  480)          zidom = zidom + ilcitl(ji,1) - nrecil
(  481)       END DO
(  482)       IF(lwp) WRITE(numout,*)
(  483)       IF(lwp) WRITE(numout,*)' sum ilcitl(i,1) = ', zidom, ' jpiglo = ', jpiglo
(  484)       
(  485)       zjdom = nrecjl
(  486)       DO jj = 1, jsplt
(  487)          zjdom = zjdom + ilcjtl(1,jj) - nrecjl
(  488)       END DO
(  489)       IF(lwp) WRITE(numout,*)' sum ilcitl(1,j) = ', zjdom, ' jpjglo = ', jpjglo
(  490)       IF(lwp) WRITE(numout,*)
(  491)       
(  492) 
(  493)       !  2. Index arrays for subdomains
(  494)       ! -------------------------------
(  495) 
(  496)       ALLOCATE(iimpptl(isplt,jsplt))
(  497)       ALLOCATE(ijmpptl(isplt,jsplt))
(  498)       
(  499)       iimpptl(:,:) = 1
(  500)       ijmpptl(:,:) = 1
(  501)       
(  502)       IF( isplt > 1 ) THEN
(  503)          DO jj = 1, jsplt
(  504)             DO ji = 2, isplt
(  505)                iimpptl(ji,jj) = iimpptl(ji-1,jj) + ilcitl(ji-1,jj) - nrecil
(  506)             END DO
(  507)          END DO
(  508)       ENDIF
(  509) 
(  510)       IF( jsplt > 1 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 10

(  511)          DO jj = 2, jsplt
(  512)             DO ji = 1, isplt
(  513)                ijmpptl(ji,jj) = ijmpptl(ji,jj-1)+ilcjtl(ji,jj-1)-nrecjl
(  514)             END DO
(  515)          END DO
(  516)       ENDIF
(  517)       
(  518)       ! 3. Subdomain description
(  519)       ! ------------------------
(  520) 
(  521)       DO jn = 1, ijsplt
(  522)          ii = 1 + MOD( jn-1, isplt )
(  523)          ij = 1 + (jn-1) / isplt
(  524)          nimpptl(jn) = iimpptl(ii,ij)
(  525)          njmpptl(jn) = ijmpptl(ii,ij)
(  526)          nlcitl (jn) = ilcitl (ii,ij)     
(  527)          nlcil       = nlcitl (jn)     
(  528)          nlcjtl (jn) = ilcjtl (ii,ij)     
(  529)          nlcjl       = nlcjtl (jn)
(  530)          nbondjl = -1                                    ! general case
(  531)          IF( jn   >  isplt          )   nbondjl = 0      ! first row of processor
(  532)          IF( jn   >  (jsplt-1)*isplt )  nbondjl = 1     ! last  row of processor
(  533)          IF( jsplt == 1             )   nbondjl = 2      ! one processor only in j-direction
(  534)          ibonjtl(jn) = nbondjl
(  535)          
(  536)          nbondil = 0                                     ! 
(  537)          IF( MOD( jn, isplt ) == 1 )   nbondil = -1      !
(  538)          IF( MOD( jn, isplt ) == 0 )   nbondil =  1      !
(  539)          IF( isplt            == 1 )   nbondil =  2      ! one processor only in i-direction
(  540)          ibonitl(jn) = nbondil
(  541)          
(  542)          nldil =  1   + jpreci
(  543)          nleil = nlcil - jpreci
(  544)          IF( nbondil == -1 .OR. nbondil == 2 )   nldil = 1
(  545)          IF( nbondil ==  1 .OR. nbondil == 2 )   nleil = nlcil
(  546)          nldjl =  1   + jprecj
(  547)          nlejl = nlcjl - jprecj
(  548)          IF( nbondjl == -1 .OR. nbondjl == 2 )   nldjl = 1
(  549)          IF( nbondjl ==  1 .OR. nbondjl == 2 )   nlejl = nlcjl
(  550)          nlditl(jn) = nldil
(  551)          nleitl(jn) = nleil
(  552)          nldjtl(jn) = nldjl
(  553)          nlejtl(jn) = nlejl
(  554)       END DO
(  555) 
(  556)       DEALLOCATE(iimpptl)
(  557)       DEALLOCATE(ijmpptl)
(  558)       DEALLOCATE(ilcitl)
(  559)       DEALLOCATE(ilcjtl)
(  560) 
(  561)    END SUBROUTINE sub_dom
(  562) 
(  563) END MODULE prtctl











PGF90 (Version     10.5)          09/23/2010  09:49:00      page 11

(    1) # 1 "prtctl.F90"
(    1) MODULE prtctl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE prtctl   ***
(    4)    !! Ocean system   : print all SUM trends for each processor domain
(    5)    !!==============================================================================
(    6)    USE dom_oce          ! ocean space and time domain variables
(    7)    USE in_out_manager   ! I/O manager
(    8)    USE lib_mpp          ! distributed memory computing
(    9) 
(   10)    IMPLICIT NONE
(   11)    PRIVATE
(   12) 
(   13)    !! * Module declaration
(   14)    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::   &  !:
(   15)       nlditl , nldjtl ,   &  !: first, last indoor index for each i-domain
(   16)       nleitl , nlejtl ,   &  !: first, last indoor index for each j-domain
(   17)       nimpptl, njmpptl,   &  !: i-, j-indexes for each processor
(   18)       nlcitl , nlcjtl ,   &  !: dimensions of every subdomain
(   19)       ibonitl, ibonjtl
(   20) 
(   21)    REAL(wp), DIMENSION(:), ALLOCATABLE ::   &  !:
(   22)       t_ctll , s_ctll ,   &  !: previous trend values
(   23)       u_ctll , v_ctll
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC prt_ctl         ! called by all subroutines
(   27)    PUBLIC prt_ctl_info    ! called by all subroutines
(   28)    PUBLIC prt_ctl_init    ! called by opa.F90
(   29)    !!----------------------------------------------------------------------
(   30)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/prtctl.F90,v 1.2 2006/03/21 15:54:15 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35) 
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE prt_ctl (tab2d_1, tab3d_1, mask1, clinfo1, tab2d_2, tab3d_2, mask2, clinfo2, ovlap, kdim, clinfo3)
(   39)       !!----------------------------------------------------------------------
(   40)       !!                     ***  ROUTINE prt_ctl  ***
(   41)       !!
(   42)       !! ** Purpose : - print sum control of 2D or 3D arrays over the same area 
(   43)       !!                in mono and mpp case. This way can be usefull when
(   44)       !!                debugging a new parametrization in mono or mpp. 
(   45)       !!
(   46)       !! ** Method  : 2 possibilities exist when setting the ln_ctl parameter to
(   47)       !!                .true. in the ocean namelist:
(   48)       !!              - to debug a MPI run .vs. a mono-processor one; 
(   49)       !!                the control print will be done over each sub-domain.
(   50)       !!                The nictl[se] and njctl[se] parameters in the namelist must 
(   51)       !!                be set to zero and [ij]splt to the corresponding splitted
(   52)       !!                domain in MPI along respectively i-, j- directions.
(   53)       !!              - to debug a mono-processor run over the whole domain/a specific area; 
(   54)       !!                in the first case the nictl[se] and njctl[se] parameters must be set
(   55)       !!                to zero else to the indices of the area to be controled. In both cases
(   56)       !!                isplt and jsplt must be set to 1.
(   57)       !!              - All arguments of the above calling sequence are optional so their






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 12

(   58)       !!                name must be explicitly typed if used. For instance if the 3D
(   59)       !!                array tn(:,:,:) must be passed through the prt_ctl subroutine, 
(   60)       !!                it must looks like: CALL prt_ctl(tab3d_1=tn).
(   61)       !!
(   62)       !!                    tab2d_1 : first 2D array
(   63)       !!                    tab3d_1 : first 3D array
(   64)       !!                    mask1   : mask (3D) to apply to the tab[23]d_1 array
(   65)       !!                    clinfo1 : information about the tab[23]d_1 array
(   66)       !!                    tab2d_2 : second 2D array
(   67)       !!                    tab3d_2 : second 3D array
(   68)       !!                    mask2   : mask (3D) to apply to the tab[23]d_2 array
(   69)       !!                    clinfo2 : information about the tab[23]d_2 array
(   70)       !!                    ovlap   : overlap value
(   71)       !!                    kdim    : k- direction for 3D arrays 
(   72)       !!                    clinfo3 : additional information 
(   73)       !!
(   74)       !! History :
(   75)       !!   9.0  !  05-07  (C. Talandier) original code
(   76)       !!----------------------------------------------------------------------
(   77)       !! * Arguments
(   78)       REAL(wp), DIMENSION(:,:)  , INTENT(in), OPTIONAL :: tab2d_1
(   79)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: tab3d_1
(   80)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask1
(   81)       CHARACTER (len=*), INTENT(in), OPTIONAL :: clinfo1
(   82)       REAL(wp), DIMENSION(:,:)  , INTENT(in), OPTIONAL :: tab2d_2
(   83)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: tab3d_2
(   84)       REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask2
(   85)       CHARACTER (len=*), INTENT(in), OPTIONAL :: clinfo2
(   86)       INTEGER, INTENT(in), OPTIONAL :: ovlap
(   87)       INTEGER, INTENT(in), OPTIONAL :: kdim
(   88)       CHARACTER (len=*), INTENT(in), OPTIONAL :: clinfo3
(   89) 
(   90)       !! * Local declarations
(   91)       INTEGER :: overlap, numid, jn, sind, eind, kdir
(   92)       CHARACTER (len=15) :: cl2
(   93)       REAL(wp) :: zsum1, zsum2, zvctl1, zvctl2
(   94)       REAL(wp), DIMENSION(jpi,jpj) :: ztab2d_1, ztab2d_2
(   95)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmask1, zmask2, ztab3d_1, ztab3d_2
(   96)       !!----------------------------------------------------------------------
(   97) 
(   98)       ! Arrays, scalars initialization 
(   99)       overlap   = 0
(  100)       kdir      = jpkm1
(  101)       cl2       = ''
(  102)       zsum1     = 0.e0
(  103)       zsum2     = 0.e0
(  104)       zvctl1    = 0.e0
(  105)       zvctl2    = 0.e0
(  106)       ztab2d_1(:,:)   = 0.e0
(  107)       ztab2d_2(:,:)   = 0.e0
(  108)       ztab3d_1(:,:,:) = 0.e0
(  109)       ztab3d_2(:,:,:) = 0.e0
(  110)       zmask1  (:,:,:) = 1.e0
(  111)       zmask2  (:,:,:) = 1.e0
(  112) 
(  113)       ! Control of optional arguments
(  114)       IF( PRESENT(clinfo2) )  cl2            = clinfo2
(  115)       IF( PRESENT(ovlap)   )  overlap        = ovlap






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 13

(  116)       IF( PRESENT(kdim)    )  kdir           = kdim
(  117)       IF( PRESENT(tab2d_1) )  ztab2d_1(:,:)  = tab2d_1(:,:)
(  118)       IF( PRESENT(tab2d_2) )  ztab2d_2(:,:)  = tab2d_2(:,:)
(  119)       IF( PRESENT(tab3d_1) )  ztab3d_1(:,:,:)= tab3d_1(:,:,:)
(  120)       IF( PRESENT(tab3d_2) )  ztab3d_2(:,:,:)= tab3d_2(:,:,:)
(  121)       IF( PRESENT(mask1)   )  zmask1  (:,:,:)= mask1  (:,:,:)
(  122)       IF( PRESENT(mask2)   )  zmask2  (:,:,:)= mask2  (:,:,:)
(  123) 
(  124)       IF( lk_mpp )   THEN
(  125)          ! processor number
(  126)          sind = narea
(  127)          eind = narea
(  128)       ELSE
(  129)          ! processors total number
(  130)          sind = 1
(  131)          eind = ijsplt
(  132)       ENDIF
(  133) 
(  134)       ! Loop over each sub-domain, i.e. the total number of processors ijsplt
(  135)       DO jn = sind, eind
(  136) 
(  137)          numid = 80 + jn
(  138) 
(  139)          ! Set indices for the SUM control
(  140)          IF( .NOT. lsp_area ) THEN
(  141)             IF (lk_mpp )   THEN
(  142)                nictls = MAX( 1, nlditl(jn) - overlap )
(  143)                nictle = nleitl(jn) + overlap * MIN( 1, nlcitl(jn) - nleitl(jn)) 
(  144)                njctls = MAX( 1, nldjtl(jn) - overlap )
(  145)                njctle = nlejtl(jn) + overlap * MIN( 1, nlcjtl(jn) - nlejtl(jn))
(  146)                ! Do not take into account the bound of the domain
(  147)                IF( ibonitl(jn) == -1 .OR. ibonitl(jn) == 2 ) nictls = MAX(2, nictls)
(  148)                IF( ibonjtl(jn) == -1 .OR. ibonjtl(jn) == 2 ) njctls = MAX(2, njctls)
(  149)                IF( ibonitl(jn) ==  1 .OR. ibonitl(jn) == 2 ) nictle = MIN(nictle, nleitl(jn) - 1)
(  150)                IF( ibonjtl(jn) ==  1 .OR. ibonjtl(jn) == 2 ) njctle = MIN(njctle, nlejtl(jn) - 1)
(  151)             ELSE
(  152)                nictls = MAX( 1, nimpptl(jn) + nlditl(jn) - 1 - overlap )
(  153)                nictle = nimpptl(jn) + nleitl(jn) - 1 + overlap * MIN( 1, nlcitl(jn) - nleitl(jn) ) 
(  154)                njctls = MAX( 1, njmpptl(jn) + nldjtl(jn) - 1 - overlap )
(  155)                njctle = njmpptl(jn) + nlejtl(jn) - 1 + overlap * MIN( 1, nlcjtl(jn) - nlejtl(jn) ) 
(  156)                ! Do not take into account the bound of the domain
(  157)                IF( ibonitl(jn) == -1 .OR. ibonitl(jn) == 2 ) nictls = MAX(2, nictls)
(  158)                IF( ibonjtl(jn) == -1 .OR. ibonjtl(jn) == 2 ) njctls = MAX(2, njctls)
(  159)                IF( ibonitl(jn) ==  1 .OR. ibonitl(jn) == 2 ) nictle = MIN(nictle, nimpptl(jn) + nleitl(jn) - 2)
(  160)                IF( ibonjtl(jn) ==  1 .OR. ibonjtl(jn) == 2 ) njctle = MIN(njctle, njmpptl(jn) + nlejtl(jn) - 2)
(  161)             ENDIF
(  162)          ENDIF
(  163) 
(  164)          IF ( clinfo3 == 'tra' )  THEN
(  165)              zvctl1 = t_ctll(jn)
(  166)              zvctl2 = s_ctll(jn)
(  167)          ELSEIF ( clinfo3 == 'dyn' )   THEN
(  168)              zvctl1 = u_ctll(jn)
(  169)              zvctl2 = v_ctll(jn)
(  170)          ENDIF
(  171) 
(  172)          ! Compute the sum control
(  173)          ! 2D arrays






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 14

(  174)          IF( PRESENT(tab2d_1) )   THEN
(  175)             zsum1 = SUM( ztab2d_1(nictls:nictle,njctls:njctle)*zmask1(nictls:nictle,njctls:njctle,1) )
(  176)             zsum2 = SUM( ztab2d_2(nictls:nictle,njctls:njctle)*zmask2(nictls:nictle,njctls:njctle,1) )
(  177)          ENDIF
(  178) 
(  179)          ! 3D arrays
(  180)          IF( PRESENT(tab3d_1) )   THEN
(  181)             zsum1 = SUM( ztab3d_1(nictls:nictle,njctls:njctle,1:kdir)*zmask1(nictls:nictle,njctls:njctle,1:kdir) )
(  182)             zsum2 = SUM( ztab3d_2(nictls:nictle,njctls:njctle,1:kdir)*zmask2(nictls:nictle,njctls:njctle,1:kdir) )
(  183)          ENDIF
(  184) 
(  185)          ! Print the result
(  186)          IF( PRESENT(clinfo2) .OR. PRESENT(tab2d_2) .OR. PRESENT(tab3d_2) )   THEN
(  187)             WRITE(numid,FMT='(a,D23.16,3x,a,D23.16)')clinfo1, zsum1-zvctl1, cl2, zsum2-zvctl2
(  188)              IF( clinfo3 == 'tra-ta' )  THEN
(  189)                 t_ctll(jn) = zsum1
(  190)              ENDIF 
(  191)              IF( clinfo3 == 'tra' )  THEN
(  192)                 t_ctll(jn) = zsum1
(  193)                 s_ctll(jn) = zsum2
(  194)              ENDIF 
(  195)              IF( clinfo3 == 'dyn' )   THEN
(  196)                 u_ctll(jn) = zsum1
(  197)                 v_ctll(jn) = zsum2 
(  198)              ENDIF
(  199)          ELSE
(  200)             WRITE(numid,FMT='(a,D23.16)')clinfo1, zsum1
(  201)          ENDIF
(  202) 
(  203)       ENDDO
(  204) 
(  205)    END SUBROUTINE prt_ctl
(  206) 
(  207) 
(  208)    SUBROUTINE prt_ctl_info (clinfo1, ivar1, clinfo2, ivar2)
(  209)       !!----------------------------------------------------------------------
(  210)       !!                     ***  ROUTINE prt_ctl_info  ***
(  211)       !!
(  212)       !! ** Purpose : - print information without any computation
(  213)       !!
(  214)       !! ** Action  : - input arguments
(  215)       !!                    clinfo1 : information about the ivar1
(  216)       !!                    ivar1   : value to print
(  217)       !!                    clinfo2 : information about the ivar2
(  218)       !!                    ivar2   : value to print
(  219)       !!
(  220)       !! History :
(  221)       !!   9.0  !  05-07  (C. Talandier) original code
(  222)       !!----------------------------------------------------------------------
(  223)       !! * Arguments
(  224)       CHARACTER (len=*), INTENT(in) ::   clinfo1
(  225)       INTEGER, INTENT(in), OPTIONAL ::   ivar1
(  226)       CHARACTER (len=*), INTENT(in), OPTIONAL ::   clinfo2
(  227)       INTEGER, INTENT(in), OPTIONAL ::   ivar2
(  228) 
(  229)       !! * Local declarations
(  230)       INTEGER ::  numid, jn, sind, eind
(  231)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 15

(  232) 
(  233)       IF( lk_mpp )   THEN
(  234)          ! processor number
(  235)          sind = narea
(  236)          eind = narea
(  237)       ELSE
(  238)          ! total number of processors
(  239)          sind = 1
(  240)          eind = ijsplt
(  241)       ENDIF
(  242) 
(  243)       ! Loop over each sub-domain, i.e. number of processors ijsplt
(  244)       DO jn = sind, eind
(  245) 
(  246)          numid = 80 + jn
(  247) 
(  248)          IF( PRESENT(ivar1) .AND. PRESENT(clinfo2) .AND. PRESENT(ivar2) )   THEN
(  249)             WRITE(numid,*)clinfo1, ivar1, clinfo2, ivar2
(  250)          ELSEIF ( PRESENT(ivar1) .AND. PRESENT(clinfo2) .AND. .NOT. PRESENT(ivar2) )   THEN
(  251)             WRITE(numid,*)clinfo1, ivar1, clinfo2
(  252)          ELSEIF ( PRESENT(ivar1) .AND. .NOT. PRESENT(clinfo2) .AND. PRESENT(ivar2) )   THEN
(  253)             WRITE(numid,*)clinfo1, ivar1, ivar2
(  254)          ELSEIF ( PRESENT(ivar1) .AND. .NOT. PRESENT(clinfo2) .AND. .NOT. PRESENT(ivar2) )   THEN
(  255)             WRITE(numid,*)clinfo1, ivar1
(  256)          ELSE
(  257)             WRITE(numid,*)clinfo1
(  258)          ENDIF
(  259) 
(  260)       ENDDO
(  261) 
(  262) 
(  263)       END SUBROUTINE prt_ctl_info
(  264) 
(  265)    SUBROUTINE prt_ctl_init
(  266)       !!----------------------------------------------------------------------
(  267)       !!                     ***  ROUTINE prt_ctl_init  ***
(  268)       !!
(  269)       !! ** Purpose :   open ASCII files & compute indices
(  270)       !!
(  271)       !! History :
(  272)       !!   9.0  !  05-07  (C. Talandier) original code
(  273)       !!----------------------------------------------------------------------
(  274)       !! * Local declarations
(  275)       INTEGER ::   jn, numid, sind, eind
(  276)       CHARACTER (len=28) :: clfile_out
(  277)       CHARACTER (len=23) :: clb_name
(  278)       CHARACTER (len=19) :: cl_run
(  279)       !!----------------------------------------------------------------------
(  280) 
(  281)       ! Allocate arrays
(  282)       ALLOCATE(nlditl (ijsplt))
(  283)       ALLOCATE(nldjtl (ijsplt))
(  284)       ALLOCATE(nleitl (ijsplt))
(  285)       ALLOCATE(nlejtl (ijsplt))
(  286)       ALLOCATE(nimpptl(ijsplt))
(  287)       ALLOCATE(njmpptl(ijsplt))
(  288)       ALLOCATE(nlcitl (ijsplt))
(  289)       ALLOCATE(nlcjtl (ijsplt))






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 16

(  290)       ALLOCATE(t_ctll (ijsplt))
(  291)       ALLOCATE(s_ctll (ijsplt))
(  292)       ALLOCATE(u_ctll (ijsplt))
(  293)       ALLOCATE(v_ctll (ijsplt))
(  294)       ALLOCATE(ibonitl(ijsplt))
(  295)       ALLOCATE(ibonjtl(ijsplt))
(  296) 
(  297)       ! Initialization 
(  298)       t_ctll (:)=0.e0
(  299)       s_ctll (:)=0.e0
(  300)       u_ctll (:)=0.e0
(  301)       v_ctll (:)=0.e0
(  302) 
(  303)       IF( lk_mpp ) THEN
(  304)          sind = narea
(  305)          eind = narea
(  306)          clb_name = "('mpp.output_',I4.4)"
(  307)          cl_run = 'MULTI processor run'
(  308)          ! use indices for each area computed by mpp_init subroutine
(  309)          nlditl(:) = nldit(:) 
(  310)          nleitl(:) = nleit(:) 
(  311)          nldjtl(:) = nldjt(:) 
(  312)          nlejtl(:) = nlejt(:) 
(  313)          !
(  314)          nimpptl(:) = nimppt(:)
(  315)          njmpptl(:) = njmppt(:)
(  316)          !
(  317)          nlcitl(:) = nlcit(:)
(  318)          nlcjtl(:) = nlcjt(:)
(  319)          !
(  320)          ibonitl(:) = ibonit(:)
(  321)          ibonjtl(:) = ibonjt(:)
(  322)       ELSE
(  323)          sind = 1
(  324)          eind = ijsplt
(  325)          clb_name = "('mono.output_',I4.4)"
(  326)          cl_run = 'MONO processor run '
(  327)          ! compute indices for each area as done in mpp_init subroutine
(  328)          CALL sub_dom
(  329)       ENDIF
(  330) 
(  331)       DO jn = sind, eind
(  332)          numid = 80 + jn
(  333)          WRITE(clfile_out,FMT=clb_name) jn-1
(  334)          OPEN ( UNIT=numid, FILE=TRIM(clfile_out),FORM='FORMATTED' )
(  335)          WRITE(numid,*)
(  336)          WRITE(numid,*) '                 L O D Y C - I P S L'
(  337)          WRITE(numid,*) '                     O P A model'
(  338)          WRITE(numid,*) '            Ocean General Circulation Model'
(  339)          WRITE(numid,*) '               version OPA 9.0  (2005) '
(  340)          WRITE(numid,*)
(  341)          WRITE(numid,*) '                   PROC number: ', jn
(  342)          WRITE(numid,*)
(  343)          WRITE(numid,FMT="(19x,a20)")cl_run
(  344) 
(  345)          ! Print the SUM control indices
(  346)          IF( .NOT. lsp_area )   THEN
(  347)             IF ( lk_mpp )   THEN






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 17

(  348)                nictls = nlditl(jn) 
(  349)                nictle = nleitl(jn)
(  350)                njctls = nldjtl(jn)
(  351)                njctle = nlejtl(jn)
(  352)             ELSE
(  353)                nictls = nimpptl(jn) + nlditl(jn) - 1
(  354)                nictle = nimpptl(jn) + nleitl(jn) - 1
(  355)                njctls = njmpptl(jn) + nldjtl(jn) - 1
(  356)                njctle = njmpptl(jn) + nlejtl(jn) - 1
(  357)             ENDIF
(  358)          ENDIF
(  359)          WRITE(numid,*) 
(  360)          WRITE(numid,*) 'prt_ctl :  Sum control indices'
(  361)          WRITE(numid,*) '~~~~~~~'
(  362)          WRITE(numid,*)
(  363)          WRITE(numid,9000)'                                nlej   = ', nlejtl(jn), '              '
(  364)          WRITE(numid,9000)'                  ------------- njctle = ', njctle, ' -------------'
(  365)          WRITE(numid,9001)'                  |                                       |'
(  366)          WRITE(numid,9001)'                  |                                       |'
(  367)          WRITE(numid,9001)'                  |                                       |'
(  368)          WRITE(numid,9002)'           nictls = ', nictls,  '                           nictle = ', nictle
(  369)          WRITE(numid,9002)'           nldi   = ', nlditl(jn),  '                           nlei   = ', nleitl(jn)
(  370)          WRITE(numid,9001)'                  |                                       |'
(  371)          WRITE(numid,9001)'                  |                                       |'
(  372)          WRITE(numid,9001)'                  |                                       |'
(  373)          WRITE(numid,9004)'  njmpp  = ',njmpptl(jn),'   ------------- njctls = ', njctls, ' -------------'
(  374)          WRITE(numid,9003)'           nimpp  = ', nimpptl(jn), '        nldj   = ', nldjtl(jn), '              '
(  375)          WRITE(numid,*)
(  376)          WRITE(numid,*)
(  377) 
(  378) 9000     FORMAT(a41,i4.4,a14)
(  379) 9001     FORMAT(a59)
(  380) 9002     FORMAT(a20,i4.4,a36,i3.3)
(  381) 9003     FORMAT(a20,i4.4,a17,i4.4)
(  382) 9004     FORMAT(a11,i4.4,a26,i4.4,a14)
(  383)       ENDDO
(  384) 
(  385)    END SUBROUTINE prt_ctl_init
(  386) 
(  387) 
(  388)    SUBROUTINE sub_dom
(  389)       !!----------------------------------------------------------------------
(  390)       !!                  ***  ROUTINE sub_dom  ***
(  391)       !!                    
(  392)       !! ** Purpose :   Lay out the global domain over processors. 
(  393)       !!                CAUTION: 
(  394)       !!                This part has been extracted from the mpp_init
(  395)       !!                subroutine and names of variables/arrays have been 
(  396)       !!                slightly changed to avoid confusion but the computation
(  397)       !!                is exactly the same. Any modification about indices of
(  398)       !!                each sub-domain in the mppini.F90 module should be reported 
(  399)       !!                here.
(  400)       !!
(  401)       !! ** Method  :   Global domain is distributed in smaller local domains.
(  402)       !!                Periodic condition is a function of the local domain position
(  403)       !!                (global boundary or neighbouring domain) and of the global
(  404)       !!                periodic
(  405)       !!                Type :         jperio global periodic condition






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 18

(  406)       !!                               nperio local  periodic condition
(  407)       !!
(  408)       !! ** Action  : - set domain parameters
(  409)       !!                    nimpp     : longitudinal index 
(  410)       !!                    njmpp     : latitudinal  index
(  411)       !!                    nperio    : lateral condition type 
(  412)       !!                    narea     : number for local area
(  413)       !!                    nlcil      : first dimension
(  414)       !!                    nlcjl      : second dimension
(  415)       !!                    nbondil    : mark for "east-west local boundary"
(  416)       !!                    nbondjl    : mark for "north-south local boundary"
(  417)       !!
(  418)       !! History :
(  419)       !!        !  94-11  (M. Guyon)  Original code
(  420)       !!        !  95-04  (J. Escobar, M. Imbard)
(  421)       !!        !  98-02  (M. Guyon)  FETI method
(  422)       !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
(  423)       !!   8.5  !  02-08  (G. Madec)  F90 : free form
(  424)       !!----------------------------------------------------------------------
(  425)       !! * Local variables
(  426)       INTEGER ::   ji, jj, jn               ! dummy loop indices
(  427)       INTEGER ::   &
(  428)          ii, ij,                         &  ! temporary integers
(  429)          irestil, irestjl,               &  !    "          "
(  430)          ijpi  , ijpj, nlcil,            &  ! temporary logical unit
(  431)          nlcjl , nbondil, nbondjl,       &
(  432)          nrecil, nrecjl, nldil, nleil, nldjl, nlejl
(  433) 
(  434)       INTEGER, DIMENSION(:,:), ALLOCATABLE ::   &
(  435)          iimpptl, ijmpptl, ilcitl, ilcjtl       ! temporary workspace
(  436)       REAL(wp) ::   zidom, zjdom            ! temporary scalars
(  437)       !!----------------------------------------------------------------------
(  438) 
(  439)       !  1. Dimension arrays for subdomains
(  440)       ! -----------------------------------
(  441)       !  Computation of local domain sizes ilcitl() ilcjtl()
(  442)       !  These dimensions depend on global sizes isplt,jsplt and jpiglo,jpjglo
(  443)       !  The subdomains are squares leeser than or equal to the global
(  444)       !  dimensions divided by the number of processors minus the overlap
(  445)       !  array (cf. par_oce.F90).
(  446) 
(  447)       ijpi = ( jpiglo-2*jpreci + (isplt-1) ) / isplt + 2*jpreci
(  448)       ijpj = ( jpjglo-2*jprecj + (jsplt-1) ) / jsplt + 2*jprecj
(  449) 
(  450)       ALLOCATE(ilcitl (isplt,jsplt))
(  451)       ALLOCATE(ilcjtl (isplt,jsplt))
(  452) 
(  453)       nrecil  = 2 * jpreci
(  454)       nrecjl  = 2 * jprecj
(  455)       irestil = MOD( jpiglo - nrecil , isplt )
(  456)       irestjl = MOD( jpjglo - nrecjl , jsplt )
(  457) 
(  458)       IF(  irestil == 0 )   irestil = isplt
(  459)       DO jj = 1, jsplt
(  460)          DO ji = 1, irestil
(  461)             ilcitl(ji,jj) = ijpi
(  462)          END DO
(  463)          DO ji = irestil+1, isplt






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 19

(  464)             ilcitl(ji,jj) = ijpi -1
(  465)          END DO
(  466)       END DO
(  467)       
(  468)       IF( irestjl == 0 )   irestjl = jsplt
(  469)       DO ji = 1, isplt
(  470)          DO jj = 1, irestjl
(  471)             ilcjtl(ji,jj) = ijpj
(  472)          END DO
(  473)          DO jj = irestjl+1, jsplt
(  474)             ilcjtl(ji,jj) = ijpj -1
(  475)          END DO
(  476)       END DO
(  477)       
(  478)       zidom = nrecil
(  479)       DO ji = 1, isplt
(  480)          zidom = zidom + ilcitl(ji,1) - nrecil
(  481)       END DO
(  482)       IF(lwp) WRITE(numout,*)
(  483)       IF(lwp) WRITE(numout,*)' sum ilcitl(i,1) = ', zidom, ' jpiglo = ', jpiglo
(  484)       
(  485)       zjdom = nrecjl
(  486)       DO jj = 1, jsplt
(  487)          zjdom = zjdom + ilcjtl(1,jj) - nrecjl
(  488)       END DO
(  489)       IF(lwp) WRITE(numout,*)' sum ilcitl(1,j) = ', zjdom, ' jpjglo = ', jpjglo
(  490)       IF(lwp) WRITE(numout,*)
(  491)       
(  492) 
(  493)       !  2. Index arrays for subdomains
(  494)       ! -------------------------------
(  495) 
(  496)       ALLOCATE(iimpptl(isplt,jsplt))
(  497)       ALLOCATE(ijmpptl(isplt,jsplt))
(  498)       
(  499)       iimpptl(:,:) = 1
(  500)       ijmpptl(:,:) = 1
(  501)       
(  502)       IF( isplt > 1 ) THEN
(  503)          DO jj = 1, jsplt
(  504)             DO ji = 2, isplt
(  505)                iimpptl(ji,jj) = iimpptl(ji-1,jj) + ilcitl(ji-1,jj) - nrecil
(  506)             END DO
(  507)          END DO
(  508)       ENDIF
(  509) 
(  510)       IF( jsplt > 1 ) THEN
(  511)          DO jj = 2, jsplt
(  512)             DO ji = 1, isplt
(  513)                ijmpptl(ji,jj) = ijmpptl(ji,jj-1)+ilcjtl(ji,jj-1)-nrecjl
(  514)             END DO
(  515)          END DO
(  516)       ENDIF
(  517)       
(  518)       ! 3. Subdomain description
(  519)       ! ------------------------
(  520) 
(  521)       DO jn = 1, ijsplt






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 20

(  522)          ii = 1 + MOD( jn-1, isplt )
(  523)          ij = 1 + (jn-1) / isplt
(  524)          nimpptl(jn) = iimpptl(ii,ij)
(  525)          njmpptl(jn) = ijmpptl(ii,ij)
(  526)          nlcitl (jn) = ilcitl (ii,ij)     
(  527)          nlcil       = nlcitl (jn)     
(  528)          nlcjtl (jn) = ilcjtl (ii,ij)     
(  529)          nlcjl       = nlcjtl (jn)
(  530)          nbondjl = -1                                    ! general case
(  531)          IF( jn   >  isplt          )   nbondjl = 0      ! first row of processor
(  532)          IF( jn   >  (jsplt-1)*isplt )  nbondjl = 1     ! last  row of processor
(  533)          IF( jsplt == 1             )   nbondjl = 2      ! one processor only in j-direction
(  534)          ibonjtl(jn) = nbondjl
(  535)          
(  536)          nbondil = 0                                     ! 
(  537)          IF( MOD( jn, isplt ) == 1 )   nbondil = -1      !
(  538)          IF( MOD( jn, isplt ) == 0 )   nbondil =  1      !
(  539)          IF( isplt            == 1 )   nbondil =  2      ! one processor only in i-direction
(  540)          ibonitl(jn) = nbondil
(  541)          
(  542)          nldil =  1   + jpreci
(  543)          nleil = nlcil - jpreci
(  544)          IF( nbondil == -1 .OR. nbondil == 2 )   nldil = 1
(  545)          IF( nbondil ==  1 .OR. nbondil == 2 )   nleil = nlcil
(  546)          nldjl =  1   + jprecj
(  547)          nlejl = nlcjl - jprecj
(  548)          IF( nbondjl == -1 .OR. nbondjl == 2 )   nldjl = 1
(  549)          IF( nbondjl ==  1 .OR. nbondjl == 2 )   nlejl = nlcjl
(  550)          nlditl(jn) = nldil
(  551)          nleitl(jn) = nleil
(  552)          nldjtl(jn) = nldjl
(  553)          nlejtl(jn) = nlejl
(  554)       END DO
(  555) 
(  556)       DEALLOCATE(iimpptl)
(  557)       DEALLOCATE(ijmpptl)
(  558)       DEALLOCATE(ilcitl)
(  559)       DEALLOCATE(ilcjtl)
(  560) 
(  561)    END SUBROUTINE sub_dom
(  562) 
(  563) END MODULE prtctl
