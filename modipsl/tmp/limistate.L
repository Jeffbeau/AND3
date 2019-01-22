


PGF90 (Version     10.5)          09/23/2010  09:49:50      page 1

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

Filename: limistate.F90

(    1) MODULE limistate
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  limistate  ***
(    4)    !!              Initialisation of diagnostics ice variables
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_istate      :  Initialisation of diagnostics ice variables
(   11)    !!   lim_istate_init :  initialization of ice state and namelist read
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst
(   15)    USE ocfzpt
(   16)    USE oce             ! dynamics and tracers variables
(   17)    USE dom_oce
(   18)    USE par_ice         ! ice parameters
(   19)    USE ice_oce         ! ice variables
(   20)    USE in_out_manager
(   21)    USE dom_ice
(   22)    USE ice
(   23)    USE lbclnk
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Accessibility
(   29)    PUBLIC lim_istate      ! routine called by lim_init.F90
(   30) 
(   31)    !! * Module variables
(   32)    REAL(wp) ::           & !!! ** init namelist (namiceini) **
(   33)       ttest  = 2.0  ,    &  ! threshold water temperature for initial sea ice
(   34)       hninn  = 0.5  ,    &  ! initial snow thickness in the north
(   35)       hginn  = 3.0  ,    &  ! initial ice thickness in the north
(   36)       alinn  = 0.05 ,    &  ! initial leads area in the north
(   37)       hnins  = 0.1  ,    &  ! initial snow thickness in the south
(   38)       hgins  = 1.0  ,    &  ! initial ice thickness in the south
(   39)       alins  = 0.1          ! initial leads area in the south
(   40) 
(   41)    REAL(wp)  ::          &  ! constant values
(   42)       zzero   = 0.e0  ,  &
(   43)       zone    = 1.e0
(   44)    !!----------------------------------------------------------------------
(   45)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   46)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limistate.F90,v 1.5 2006/03/21 08:38:38 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 2

(   47)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   48)    !!----------------------------------------------------------------------
(   49) 
(   50) CONTAINS
(   51) 
(   52)    SUBROUTINE lim_istate
(   53)       !!-------------------------------------------------------------------
(   54)       !!                    ***  ROUTINE lim_istate  ***
(   55)       !!
(   56)       !! ** Purpose :   defined the sea-ice initial state
(   57)       !!
(   58)       !! ** Method  :   restart from a state defined in a binary file
(   59)       !!                or from arbitrary sea-ice conditions
(   60)       !!
(   61)       !! History :
(   62)       !!   2.0  !  01-04  (C. Ethe, G. Madec)  Original code
(   63)       !!        !  04-04  (S. Theetten) initialization from a file
(   64)       !!--------------------------------------------------------------------
(   65)       !! * Local variables
(   66)       INTEGER  ::   ji, jj, jk   ! dummy loop indices
(   67)       REAL(wp) ::   zidto,    &  ! temporary scalar
(   68)          zs0, ztf, zbin
(   69)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   70)          ztn
(   71)       !--------------------------------------------------------------------
(   72) 
(   73)  
(   74)       CALL lim_istate_init     !  reading the initials parameters of the ice
(   75) 
(   76)       !-- Initialisation of sst,sss,u,v do i=1,jpi
(   77)       u_io(:,:)  = 0.e0       ! ice velocity in x direction
(   78)       v_io(:,:)  = 0.e0       ! ice velocity in y direction
(   79) 
(   80)       IF( ln_limini ) THEN    ! 
(   81)         
(   82)          ! Initialisation at tn if no ice or sst_ini if ice
(   83)          ! Idem for salinity
(   84) 
(   85)       !--- Criterion for presence (zidto=1.) or absence (zidto=0.) of ice
(   86)          DO jj = 1 , jpj
(   87)             DO ji = 1 , jpi
(   88)                
(   89)                zidto = MAX(zzero, - SIGN(1.,frld(ji,jj) - 1.))
(   90)                
(   91)                sst_io(ji,jj) = ( nfice - 1 ) * (zidto * sst_ini(ji,jj)  + &   ! use the ocean initial values
(   92)                     &          (1.0 - zidto ) * ( tn(ji,jj,1) + rt0 ))        ! tricky trick *(nfice-1) !
(   93)                sss_io(ji,jj) = ( nfice - 1 ) * (zidto * sss_ini(ji,jj) + &
(   94)                     &          (1.0 - zidto ) *  sn(ji,jj,1) )
(   95) 
(   96)                ! to avoid the the melting of ice, several layers (mixed layer) should be
(   97)                ! set to sst_ini (sss_ini) if there is ice
(   98)                ! example for one layer 
(   99)                ! tn(ji,jj,1) = zidto * ( sst_ini(ji,jj) - rt0 )  + (1.0 - zidto ) *  tn(ji,jj,1)
(  100)                ! sn(ji,jj,1) = zidto * sss_ini(ji,jj)  + (1.0 - zidto ) *  sn(ji,jj,1)
(  101)                ! tb(ji,jj,1) = tn(ji,jj,1)
(  102)                ! sb(ji,jj,1) = sn(ji,jj,1)
(  103)             END DO
(  104)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 3

(  105)          
(  106)          
(  107)          !  tfu: Melting point of sea water
(  108)          tfu(:,:)  = ztf   
(  109)          
(  110)          tfu(:,:)  = ABS ( rt0 - 0.0575       * sss_ini(:,:)                               &
(  111)               &                    + 1.710523e-03 * sss_ini(:,:) * SQRT( sss_ini(:,:) )    &
(  112)               &                    - 2.154996e-04 * sss_ini(:,:) * sss_ini(:,:) )
(  113)       ELSE                     !
(  114) 
(  115)          
(  116)          ! Initialisation at tn or -2 if ice
(  117)          DO jj = 1, jpj
(  118)             DO ji = 1, jpi
(  119)                zbin = MAX( 0., SIGN( 1., fzptn(ji,jj) - tn(ji,jj,1) ) )
(  120)                ztn(ji,jj) = ( (1.-zbin) * tn(ji,jj,1) - 2. * zbin + rt0 ) * tmask(ji,jj,1)
(  121)             END DO
(  122)          END DO
(  123)          
(  124)          u_io  (:,:) = 0.e0
(  125)          v_io  (:,:) = 0.e0
(  126)          sst_io(:,:) = ( nfice - 1 ) * ( tn(:,:,1) + rt0 )   ! use the ocean initial values
(  127)          sss_io(:,:) = ( nfice - 1 ) *   sn(:,:,1)           ! tricky trick *(nfice-1) !
(  128)          
(  129)          ! reference salinity 34psu
(  130)          zs0 = 34.e0
(  131)          ztf = ABS ( rt0 - 0.0575       * zs0                           &
(  132)               &                    + 1.710523e-03 * zs0 * SQRT( zs0 )   &
(  133)               &                    - 2.154996e-04 * zs0 *zs0          )
(  134)          
(  135)          !  tfu: Melting point of sea water
(  136)          tfu(:,:)  = ztf   
(  137)          
(  138)          DO jj = 1, jpj
(  139)             DO ji = 1, jpi
(  140)                !--- Criterion for presence (zidto=1) or absence (zidto=0) of ice
(  141)                zidto  = tms(ji,jj) * ( 1.0 - MAX(zzero, SIGN( zone, ztn(ji,jj) - tfu(ji,jj) - ttest) ) )
(  142)                
(  143)                IF( fcor(ji,jj) >= 0.e0 ) THEN     !--  Northern hemisphere.
(  144)                   hicif(ji,jj)   = zidto * hginn
(  145)                   frld(ji,jj)    = zidto * alinn + ( 1.0 - zidto ) * 1.0
(  146)                   hsnif(ji,jj)   = zidto * hninn
(  147)                ELSE                               !---  Southern hemisphere.
(  148)                   hicif(ji,jj)   = zidto * hgins
(  149)                   frld(ji,jj)    = zidto * alins + ( 1.0 - zidto ) * 1.0
(  150)                   hsnif(ji,jj)   = zidto * hnins
(  151)                ENDIF
(  152)             END DO
(  153)          END DO
(  154)          
(  155)          sist  (:,:)   = tfu(:,:)
(  156)          tbif  (:,:,1) = tfu(:,:)
(  157)          tbif  (:,:,2) = tfu(:,:)
(  158)          tbif  (:,:,3) = tfu(:,:)
(  159)       
(  160)       ENDIF
(  161)       fsbbq (:,:)   = 0.e0
(  162)       qstoif(:,:)   = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 4

(  163)       u_ice (:,:)   = 0.e0
(  164)       v_ice (:,:)   = 0.e0
(  165) # if defined key_coupled
(  166)       albege(:,:)   = 0.8 * tms(:,:)
(  167) # endif
(  168) 
(  169)       !---  Moments for advection.             
(  170) 
(  171)       sxice (:,:)  = 0.e0   ;   sxsn (:,:)  = 0.e0   ;   sxa  (:,:)  = 0.e0
(  172)       syice (:,:)  = 0.e0   ;   sysn (:,:)  = 0.e0   ;   sya  (:,:)  = 0.e0
(  173)       sxxice(:,:)  = 0.e0   ;   sxxsn(:,:)  = 0.e0   ;   sxxa (:,:)  = 0.e0
(  174)       syyice(:,:)  = 0.e0   ;   syysn(:,:)  = 0.e0   ;   syya (:,:)  = 0.e0
(  175)       sxyice(:,:)  = 0.e0   ;   sxysn(:,:)  = 0.e0   ;   sxya (:,:)  = 0.e0
(  176) 
(  177)       sxc0  (:,:)  = 0.e0   ;   sxc1 (:,:)  = 0.e0   ;   sxc2 (:,:)  = 0.e0
(  178)       syc0  (:,:)  = 0.e0   ;   syc1 (:,:)  = 0.e0   ;   syc2 (:,:)  = 0.e0
(  179)       sxxc0 (:,:)  = 0.e0   ;   sxxc1(:,:)  = 0.e0   ;   sxxc2(:,:)  = 0.e0
(  180)       syyc0 (:,:)  = 0.e0   ;   syyc1(:,:)  = 0.e0   ;   syyc2(:,:)  = 0.e0
(  181)       sxyc0 (:,:)  = 0.e0   ;   sxyc1(:,:)  = 0.e0   ;   sxyc2(:,:)  = 0.e0
(  182) 
(  183)       sxst  (:,:)  = 0.e0
(  184)       syst  (:,:)  = 0.e0
(  185)       sxxst (:,:)  = 0.e0
(  186)       syyst (:,:)  = 0.e0
(  187)       sxyst (:,:)  = 0.e0
(  188) 
(  189)       !-- lateral boundary conditions
(  190)       CALL lbc_lnk( hicif, 'T', 1. )
(  191)       CALL lbc_lnk( frld , 'T', 1. )
(  192) 
(  193)       ! C A U T I O N  frld = 1 over land and lbc_lnk put zero along 
(  194)       ! *************  closed boundaries herefore we force to one over land
(  195)       frld(:,:) = tms(:,:) * frld(:,:) + ( 1. - tms(:,:) )   
(  196) 
(  197)       CALL lbc_lnk( hsnif, 'T', 1. )
(  198)       CALL lbc_lnk( sist , 'T', 1. )
(  199)       DO jk = 1, jplayersp1
(  200)          CALL lbc_lnk(tbif(:,:,jk), 'T', 1. )
(  201)       END DO
(  202)       CALL lbc_lnk( fsbbq  , 'T', 1. )
(  203)       CALL lbc_lnk( qstoif , 'T', 1. )
(  204)       CALL lbc_lnk( sss_io , 'T', 1. )
(  205) 
(  206)    END SUBROUTINE lim_istate
(  207) 
(  208)    
(  209)    SUBROUTINE lim_istate_init
(  210)       !!-------------------------------------------------------------------
(  211)       !!                   ***  ROUTINE lim_istate_init  ***
(  212)       !!        
(  213)       !! ** Purpose :   Definition of initial state of the ice 
(  214)       !!
(  215)       !! ** Method  :   Read the namiceini namelist and check the parameter 
(  216)       !!                values called at the first timestep (nit000)
(  217)       !!                or
(  218)       !!                Read 7 variables from a previous restart file
(  219)       !!                sst, sst, hicif, hsnif, frld, ts & tbif
(  220)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 5

(  221)       !! ** input   :   Namelist namiceini
(  222)       !!
(  223)       !! history
(  224)       !!  8.5  ! 03-08 (C. Ethe) original code
(  225)       !!  9.0  ! 04-04 (S. Theetten) read a file
(  226)       !!-------------------------------------------------------------------
(  227)       !! * Modules used
(  228)       USE ice
(  229)       USE ioipsl
(  230) 
(  231)       NAMELIST/namiceini/ ln_limini, ln_limdmp, ttest, hninn, hginn, alinn, &
(  232)          &                hnins, hgins, alins
(  233)       !!-------------------------------------------------------------------
(  234)       !! local declaration
(  235)       INTEGER, PARAMETER ::   jpmois=1
(  236)       
(  237)       INTEGER ::                   &
(  238)            itime, ipi, ipj, ipk  , & ! temporary integers
(  239)            inum_ice
(  240)       
(  241)       INTEGER ::  istep(jpmois)
(  242)       
(  243)       REAL(wp) ::   zdate0, zdt
(  244)       REAL(wp), DIMENSION(jpi,jpj) ::   zlon, zlat
(  245)       REAL(wp), DIMENSION(3) ::   zlev
(  246)       
(  247)       CHARACTER (len=32) :: cl_icedata
(  248)       
(  249)       LOGICAL :: llbon
(  250)       !!-------------------------------------------------------------------
(  251)       
(  252)       ! Read Namelist namiceini 
(  253) 
(  254)       REWIND ( numnam_ice )
(  255)       READ   ( numnam_ice , namiceini )
(  256)       
(  257)       IF(.NOT. ln_limini) THEN 
(  258)          IF(lwp) THEN
(  259)             WRITE(numout,*)
(  260)             WRITE(numout,*) 'lim_istate_init : ice parameters inititialisation '
(  261)             WRITE(numout,*) '~~~~~~~~~~~~~~~'
(  262)             WRITE(numout,*) '         threshold water temp. for initial sea-ice    ttest      = ', ttest
(  263)             WRITE(numout,*) '         initial snow thickness in the north          hninn      = ', hninn
(  264)             WRITE(numout,*) '         initial ice thickness in the north           hginn      = ', hginn 
(  265)             WRITE(numout,*) '         initial leads area in the north              alinn      = ', alinn            
(  266)             WRITE(numout,*) '         initial snow thickness in the south          hnins      = ', hnins 
(  267)             WRITE(numout,*) '         initial ice thickness in the south           hgins      = ', hgins
(  268)             WRITE(numout,*) '         initial leads area in the south              alins      = ', alins
(  269)          ENDIF
(  270)       ENDIF
(  271) 
(  272)       IF( ln_limini ) THEN                      ! Ice initialization using input file
(  273) 
(  274)          cl_icedata = 'Ice_initialization.nc'
(  275)          INQUIRE( FILE=cl_icedata, EXIST=llbon )
(  276)          IF( llbon ) THEN
(  277)             IF(lwp) THEN
(  278)                WRITE(numout,*) ' '






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 6

(  279)                WRITE(numout,*) 'lim_istate_init : ice state initialization with : ',cl_icedata
(  280)                WRITE(numout,*) '~~~~~~~~~~~~~~~'
(  281)                WRITE(numout,*) '         Ice state initialization using input file    ln_limini  = ', ln_limini
(  282)                WRITE(numout,*) '         Ice damping                                  ln_limdmp  = ', ln_limdmp
(  283)                WRITE(numout,*) ' '
(  284)             ENDIF
(  285)             
(  286)             itime = 1
(  287)             ipi=jpiglo
(  288)             ipj=jpjglo
(  289)             ipk=1
(  290)             zdt=rdt
(  291)             
(  292)             CALL flinopen( TRIM(cl_icedata), mig(1), nlci, mjg(1), nlcj, .FALSE., &
(  293)                &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, inum_ice )
(  294)             
(  295)             CALL flinget( inum_ice, 'sst', jpidta, jpjdta, 1,  &
(  296)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sst_ini(1:nlci,1:nlcj) )
(  297)             
(  298)             CALL flinget( inum_ice, 'sss', jpidta, jpjdta, 1,  &
(  299)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sss_ini(1:nlci,1:nlcj) )
(  300)             
(  301)             CALL flinget( inum_ice, 'hicif', jpidta, jpjdta, 1,  &
(  302)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, hicif(1:nlci,1:nlcj) )
(  303)             
(  304)             CALL flinget( inum_ice, 'hsnif', jpidta, jpjdta, 1,  &
(  305)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, hsnif(1:nlci,1:nlcj) )
(  306)             
(  307)             CALL flinget( inum_ice, 'frld', jpidta, jpjdta, 1,  &
(  308)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, frld(1:nlci,1:nlcj) )
(  309)             
(  310)             CALL flinget( inum_ice, 'ts', jpidta, jpjdta, 1,  &
(  311)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sist(1:nlci,1:nlcj) )
(  312)             
(  313)             CALL flinclo( inum_ice)
(  314)             
(  315)             itime = 1
(  316)             ipi=jpiglo
(  317)             ipj=jpjglo
(  318)             ipk=jplayersp1
(  319)             
(  320)             CALL flinopen( TRIM(cl_icedata), mig(1), nlci, mjg(1), nlcj, .FALSE.,  &
(  321)                &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, inum_ice )
(  322)             
(  323)             CALL flinget( inum_ice, 'tbif', jpidta, jpjdta, ipk,  &
(  324)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, tbif(1:nlci,1:nlcj,1:ipk) )
(  325)             
(  326)             CALL flinclo( inum_ice)
(  327)             
(  328)          ELSE
(  329)             IF(lwp) WRITE(numout,cform_err) 
(  330)             IF(lwp) WRITE(numout,*) '            ',cl_icedata, ' not found !'
(  331)             nstop = nstop + 1
(  332)          ENDIF
(  333)       ENDIF
(  334) 
(  335)    END SUBROUTINE lim_istate_init
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 7

(  337) #else
(  338)    !!----------------------------------------------------------------------
(  339)    !!   Default option :         Empty module          NO LIM sea-ice model
(  340)    !!----------------------------------------------------------------------
(  341) CONTAINS
(  342)    SUBROUTINE lim_istate          ! Empty routine
(  343)    END SUBROUTINE lim_istate
(  344) #endif
(  345) 
(  346)    !!======================================================================
(  347) END MODULE limistate





















































PGF90 (Version     10.5)          09/23/2010  09:49:50      page 8

(    1) # 1 "limistate.F90"
(    1) MODULE limistate
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  limistate  ***
(    4)    !!              Initialisation of diagnostics ice variables
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_istate      :  Initialisation of diagnostics ice variables
(   11)    !!   lim_istate_init :  initialization of ice state and namelist read
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst
(   15)    USE ocfzpt
(   16)    USE oce             ! dynamics and tracers variables
(   17)    USE dom_oce
(   18)    USE par_ice         ! ice parameters
(   19)    USE ice_oce         ! ice variables
(   20)    USE in_out_manager
(   21)    USE dom_ice
(   22)    USE ice
(   23)    USE lbclnk
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Accessibility
(   29)    PUBLIC lim_istate      ! routine called by lim_init.F90
(   30) 
(   31)    !! * Module variables
(   32)    REAL(wp) ::           & !!! ** init namelist (namiceini) **
(   33)       ttest  = 2.0  ,    &  ! threshold water temperature for initial sea ice
(   34)       hninn  = 0.5  ,    &  ! initial snow thickness in the north
(   35)       hginn  = 3.0  ,    &  ! initial ice thickness in the north
(   36)       alinn  = 0.05 ,    &  ! initial leads area in the north
(   37)       hnins  = 0.1  ,    &  ! initial snow thickness in the south
(   38)       hgins  = 1.0  ,    &  ! initial ice thickness in the south
(   39)       alins  = 0.1          ! initial leads area in the south
(   40) 
(   41)    REAL(wp)  ::          &  ! constant values
(   42)       zzero   = 0.e0  ,  &
(   43)       zone    = 1.e0
(   44)    !!----------------------------------------------------------------------
(   45)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   46)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limistate.F90,v 1.5 2006/03/21 08:38:38 opalod Exp $ 
(   47)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   48)    !!----------------------------------------------------------------------
(   49) 
(   50) CONTAINS
(   51) 
(   52)    SUBROUTINE lim_istate
(   53)       !!-------------------------------------------------------------------
(   54)       !!                    ***  ROUTINE lim_istate  ***
(   55)       !!
(   56)       !! ** Purpose :   defined the sea-ice initial state
(   57)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 9

(   58)       !! ** Method  :   restart from a state defined in a binary file
(   59)       !!                or from arbitrary sea-ice conditions
(   60)       !!
(   61)       !! History :
(   62)       !!   2.0  !  01-04  (C. Ethe, G. Madec)  Original code
(   63)       !!        !  04-04  (S. Theetten) initialization from a file
(   64)       !!--------------------------------------------------------------------
(   65)       !! * Local variables
(   66)       INTEGER  ::   ji, jj, jk   ! dummy loop indices
(   67)       REAL(wp) ::   zidto,    &  ! temporary scalar
(   68)          zs0, ztf, zbin
(   69)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   70)          ztn
(   71)       !--------------------------------------------------------------------
(   72) 
(   73)  
(   74)       CALL lim_istate_init     !  reading the initials parameters of the ice
(   75) 
(   76)       !-- Initialisation of sst,sss,u,v do i=1,jpi
(   77)       u_io(:,:)  = 0.e0       ! ice velocity in x direction
(   78)       v_io(:,:)  = 0.e0       ! ice velocity in y direction
(   79) 
(   80)       IF( ln_limini ) THEN    ! 
(   81)         
(   82)          ! Initialisation at tn if no ice or sst_ini if ice
(   83)          ! Idem for salinity
(   84) 
(   85)       !--- Criterion for presence (zidto=1.) or absence (zidto=0.) of ice
(   86)          DO jj = 1 , jpj
(   87)             DO ji = 1 , jpi
(   88)                
(   89)                zidto = MAX(zzero, - SIGN(1.,frld(ji,jj) - 1.))
(   90)                
(   91)                sst_io(ji,jj) = ( nfice - 1 ) * (zidto * sst_ini(ji,jj)  + &   ! use the ocean initial values
(   92)                     &          (1.0 - zidto ) * ( tn(ji,jj,1) + rt0 ))        ! tricky trick *(nfice-1) !
(   93)                sss_io(ji,jj) = ( nfice - 1 ) * (zidto * sss_ini(ji,jj) + &
(   94)                     &          (1.0 - zidto ) *  sn(ji,jj,1) )
(   95) 
(   96)                ! to avoid the the melting of ice, several layers (mixed layer) should be
(   97)                ! set to sst_ini (sss_ini) if there is ice
(   98)                ! example for one layer 
(   99)                ! tn(ji,jj,1) = zidto * ( sst_ini(ji,jj) - rt0 )  + (1.0 - zidto ) *  tn(ji,jj,1)
(  100)                ! sn(ji,jj,1) = zidto * sss_ini(ji,jj)  + (1.0 - zidto ) *  sn(ji,jj,1)
(  101)                ! tb(ji,jj,1) = tn(ji,jj,1)
(  102)                ! sb(ji,jj,1) = sn(ji,jj,1)
(  103)             END DO
(  104)          END DO
(  105)          
(  106)          
(  107)          !  tfu: Melting point of sea water
(  108)          tfu(:,:)  = ztf   
(  109)          
(  110)          tfu(:,:)  = ABS ( rt0 - 0.0575       * sss_ini(:,:)                               &
(  111)               &                    + 1.710523e-03 * sss_ini(:,:) * SQRT( sss_ini(:,:) )    &
(  112)               &                    - 2.154996e-04 * sss_ini(:,:) * sss_ini(:,:) )
(  113)       ELSE                     !
(  114) 
(  115)          






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 10

(  116)          ! Initialisation at tn or -2 if ice
(  117)          DO jj = 1, jpj
(  118)             DO ji = 1, jpi
(  119)                zbin = MAX( 0., SIGN( 1., fzptn(ji,jj) - tn(ji,jj,1) ) )
(  120)                ztn(ji,jj) = ( (1.-zbin) * tn(ji,jj,1) - 2. * zbin + rt0 ) * tmask(ji,jj,1)
(  121)             END DO
(  122)          END DO
(  123)          
(  124)          u_io  (:,:) = 0.e0
(  125)          v_io  (:,:) = 0.e0
(  126)          sst_io(:,:) = ( nfice - 1 ) * ( tn(:,:,1) + rt0 )   ! use the ocean initial values
(  127)          sss_io(:,:) = ( nfice - 1 ) *   sn(:,:,1)           ! tricky trick *(nfice-1) !
(  128)          
(  129)          ! reference salinity 34psu
(  130)          zs0 = 34.e0
(  131)          ztf = ABS ( rt0 - 0.0575       * zs0                           &
(  132)               &                    + 1.710523e-03 * zs0 * SQRT( zs0 )   &
(  133)               &                    - 2.154996e-04 * zs0 *zs0          )
(  134)          
(  135)          !  tfu: Melting point of sea water
(  136)          tfu(:,:)  = ztf   
(  137)          
(  138)          DO jj = 1, jpj
(  139)             DO ji = 1, jpi
(  140)                !--- Criterion for presence (zidto=1) or absence (zidto=0) of ice
(  141)                zidto  = tms(ji,jj) * ( 1.0 - MAX(zzero, SIGN( zone, ztn(ji,jj) - tfu(ji,jj) - ttest) ) )
(  142)                
(  143)                IF( fcor(ji,jj) >= 0.e0 ) THEN     !--  Northern hemisphere.
(  144)                   hicif(ji,jj)   = zidto * hginn
(  145)                   frld(ji,jj)    = zidto * alinn + ( 1.0 - zidto ) * 1.0
(  146)                   hsnif(ji,jj)   = zidto * hninn
(  147)                ELSE                               !---  Southern hemisphere.
(  148)                   hicif(ji,jj)   = zidto * hgins
(  149)                   frld(ji,jj)    = zidto * alins + ( 1.0 - zidto ) * 1.0
(  150)                   hsnif(ji,jj)   = zidto * hnins
(  151)                ENDIF
(  152)             END DO
(  153)          END DO
(  154)          
(  155)          sist  (:,:)   = tfu(:,:)
(  156)          tbif  (:,:,1) = tfu(:,:)
(  157)          tbif  (:,:,2) = tfu(:,:)
(  158)          tbif  (:,:,3) = tfu(:,:)
(  159)       
(  160)       ENDIF
(  161)       fsbbq (:,:)   = 0.e0
(  162)       qstoif(:,:)   = 0.e0
(  163)       u_ice (:,:)   = 0.e0
(  164)       v_ice (:,:)   = 0.e0
(  165) 
(  166) # 169
(  169)       !---  Moments for advection.             
(  170) 
(  171)       sxice (:,:)  = 0.e0   ;   sxsn (:,:)  = 0.e0   ;   sxa  (:,:)  = 0.e0
(  172)       syice (:,:)  = 0.e0   ;   sysn (:,:)  = 0.e0   ;   sya  (:,:)  = 0.e0
(  173)       sxxice(:,:)  = 0.e0   ;   sxxsn(:,:)  = 0.e0   ;   sxxa (:,:)  = 0.e0
(  174)       syyice(:,:)  = 0.e0   ;   syysn(:,:)  = 0.e0   ;   syya (:,:)  = 0.e0
(  175)       sxyice(:,:)  = 0.e0   ;   sxysn(:,:)  = 0.e0   ;   sxya (:,:)  = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 11

(  176) 
(  177)       sxc0  (:,:)  = 0.e0   ;   sxc1 (:,:)  = 0.e0   ;   sxc2 (:,:)  = 0.e0
(  178)       syc0  (:,:)  = 0.e0   ;   syc1 (:,:)  = 0.e0   ;   syc2 (:,:)  = 0.e0
(  179)       sxxc0 (:,:)  = 0.e0   ;   sxxc1(:,:)  = 0.e0   ;   sxxc2(:,:)  = 0.e0
(  180)       syyc0 (:,:)  = 0.e0   ;   syyc1(:,:)  = 0.e0   ;   syyc2(:,:)  = 0.e0
(  181)       sxyc0 (:,:)  = 0.e0   ;   sxyc1(:,:)  = 0.e0   ;   sxyc2(:,:)  = 0.e0
(  182) 
(  183)       sxst  (:,:)  = 0.e0
(  184)       syst  (:,:)  = 0.e0
(  185)       sxxst (:,:)  = 0.e0
(  186)       syyst (:,:)  = 0.e0
(  187)       sxyst (:,:)  = 0.e0
(  188) 
(  189)       !-- lateral boundary conditions
(  190)       CALL lbc_lnk( hicif, 'T', 1. )
(  191)       CALL lbc_lnk( frld , 'T', 1. )
(  192) 
(  193)       ! C A U T I O N  frld = 1 over land and lbc_lnk put zero along 
(  194)       ! *************  closed boundaries herefore we force to one over land
(  195)       frld(:,:) = tms(:,:) * frld(:,:) + ( 1. - tms(:,:) )   
(  196) 
(  197)       CALL lbc_lnk( hsnif, 'T', 1. )
(  198)       CALL lbc_lnk( sist , 'T', 1. )
(  199)       DO jk = 1, jplayersp1
(  200)          CALL lbc_lnk(tbif(:,:,jk), 'T', 1. )
(  201)       END DO
(  202)       CALL lbc_lnk( fsbbq  , 'T', 1. )
(  203)       CALL lbc_lnk( qstoif , 'T', 1. )
(  204)       CALL lbc_lnk( sss_io , 'T', 1. )
(  205) 
(  206)    END SUBROUTINE lim_istate
(  207) 
(  208)    
(  209)    SUBROUTINE lim_istate_init
(  210)       !!-------------------------------------------------------------------
(  211)       !!                   ***  ROUTINE lim_istate_init  ***
(  212)       !!        
(  213)       !! ** Purpose :   Definition of initial state of the ice 
(  214)       !!
(  215)       !! ** Method  :   Read the namiceini namelist and check the parameter 
(  216)       !!                values called at the first timestep (nit000)
(  217)       !!                or
(  218)       !!                Read 7 variables from a previous restart file
(  219)       !!                sst, sst, hicif, hsnif, frld, ts & tbif
(  220)       !!
(  221)       !! ** input   :   Namelist namiceini
(  222)       !!
(  223)       !! history
(  224)       !!  8.5  ! 03-08 (C. Ethe) original code
(  225)       !!  9.0  ! 04-04 (S. Theetten) read a file
(  226)       !!-------------------------------------------------------------------
(  227)       !! * Modules used
(  228)       USE ice
(  229)       USE ioipsl
(  230) 
(  231)       NAMELIST/namiceini/ ln_limini, ln_limdmp, ttest, hninn, hginn, alinn, &
(  232)          &                hnins, hgins, alins
(  233)       !!-------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 12

(  234)       !! local declaration
(  235)       INTEGER, PARAMETER ::   jpmois=1
(  236)       
(  237)       INTEGER ::                   &
(  238)            itime, ipi, ipj, ipk  , & ! temporary integers
(  239)            inum_ice
(  240)       
(  241)       INTEGER ::  istep(jpmois)
(  242)       
(  243)       REAL(wp) ::   zdate0, zdt
(  244)       REAL(wp), DIMENSION(jpi,jpj) ::   zlon, zlat
(  245)       REAL(wp), DIMENSION(3) ::   zlev
(  246)       
(  247)       CHARACTER (len=32) :: cl_icedata
(  248)       
(  249)       LOGICAL :: llbon
(  250)       !!-------------------------------------------------------------------
(  251)       
(  252)       ! Read Namelist namiceini 
(  253) 
(  254)       REWIND ( numnam_ice )
(  255)       READ   ( numnam_ice , namiceini )
(  256)       
(  257)       IF(.NOT. ln_limini) THEN 
(  258)          IF(lwp) THEN
(  259)             WRITE(numout,*)
(  260)             WRITE(numout,*) 'lim_istate_init : ice parameters inititialisation '
(  261)             WRITE(numout,*) '~~~~~~~~~~~~~~~'
(  262)             WRITE(numout,*) '         threshold water temp. for initial sea-ice    ttest      = ', ttest
(  263)             WRITE(numout,*) '         initial snow thickness in the north          hninn      = ', hninn
(  264)             WRITE(numout,*) '         initial ice thickness in the north           hginn      = ', hginn 
(  265)             WRITE(numout,*) '         initial leads area in the north              alinn      = ', alinn            
(  266)             WRITE(numout,*) '         initial snow thickness in the south          hnins      = ', hnins 
(  267)             WRITE(numout,*) '         initial ice thickness in the south           hgins      = ', hgins
(  268)             WRITE(numout,*) '         initial leads area in the south              alins      = ', alins
(  269)          ENDIF
(  270)       ENDIF
(  271) 
(  272)       IF( ln_limini ) THEN                      ! Ice initialization using input file
(  273) 
(  274)          cl_icedata = 'Ice_initialization.nc'
(  275)          INQUIRE( FILE=cl_icedata, EXIST=llbon )
(  276)          IF( llbon ) THEN
(  277)             IF(lwp) THEN
(  278)                WRITE(numout,*) ' '
(  279)                WRITE(numout,*) 'lim_istate_init : ice state initialization with : ',cl_icedata
(  280)                WRITE(numout,*) '~~~~~~~~~~~~~~~'
(  281)                WRITE(numout,*) '         Ice state initialization using input file    ln_limini  = ', ln_limini
(  282)                WRITE(numout,*) '         Ice damping                                  ln_limdmp  = ', ln_limdmp
(  283)                WRITE(numout,*) ' '
(  284)             ENDIF
(  285)             
(  286)             itime = 1
(  287)             ipi=jpiglo
(  288)             ipj=jpjglo
(  289)             ipk=1
(  290)             zdt=rdt
(  291)             






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 13

(  292)             CALL flinopen( TRIM(cl_icedata), mig(1), nlci, mjg(1), nlcj, .FALSE., &
(  293)                &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, inum_ice )
(  294)             
(  295)             CALL flinget( inum_ice, 'sst', jpidta, jpjdta, 1,  &
(  296)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sst_ini(1:nlci,1:nlcj) )
(  297)             
(  298)             CALL flinget( inum_ice, 'sss', jpidta, jpjdta, 1,  &
(  299)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sss_ini(1:nlci,1:nlcj) )
(  300)             
(  301)             CALL flinget( inum_ice, 'hicif', jpidta, jpjdta, 1,  &
(  302)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, hicif(1:nlci,1:nlcj) )
(  303)             
(  304)             CALL flinget( inum_ice, 'hsnif', jpidta, jpjdta, 1,  &
(  305)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, hsnif(1:nlci,1:nlcj) )
(  306)             
(  307)             CALL flinget( inum_ice, 'frld', jpidta, jpjdta, 1,  &
(  308)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, frld(1:nlci,1:nlcj) )
(  309)             
(  310)             CALL flinget( inum_ice, 'ts', jpidta, jpjdta, 1,  &
(  311)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, sist(1:nlci,1:nlcj) )
(  312)             
(  313)             CALL flinclo( inum_ice)
(  314)             
(  315)             itime = 1
(  316)             ipi=jpiglo
(  317)             ipj=jpjglo
(  318)             ipk=jplayersp1
(  319)             
(  320)             CALL flinopen( TRIM(cl_icedata), mig(1), nlci, mjg(1), nlcj, .FALSE.,  &
(  321)                &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, inum_ice )
(  322)             
(  323)             CALL flinget( inum_ice, 'tbif', jpidta, jpjdta, ipk,  &
(  324)                &          jpmois, 1, 0, mig(1), nlci, mjg(1), nlcj, tbif(1:nlci,1:nlcj,1:ipk) )
(  325)             
(  326)             CALL flinclo( inum_ice)
(  327)             
(  328)          ELSE
(  329)             IF(lwp) WRITE(numout,cform_err) 
(  330)             IF(lwp) WRITE(numout,*) '            ',cl_icedata, ' not found !'
(  331)             nstop = nstop + 1
(  332)          ENDIF
(  333)       ENDIF
(  334) 
(  335)    END SUBROUTINE lim_istate_init
(  336) 
(  337) 
(  338) # 346
(  346)    !!======================================================================
(  347) END MODULE limistate
