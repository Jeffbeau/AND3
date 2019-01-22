


PGF90 (Version     10.5)          09/23/2010  09:49:16      page 1

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

Filename: obc_oce.F90

(    1) !!DB 2008.04.14
(    2) !!see obc_oce.Q90 for last version of this program
(    3) 
(    4) MODULE obc_oce
(    5)   !!==============================================================================
(    6)   !!                       ***  MODULE obc_oce   ***
(    7)   !! Open Boundary Cond. :   define related variables
(    8)   !!==============================================================================
(    9) #if defined key_obc
(   10)   !!----------------------------------------------------------------------
(   11)   !!   'key_obc' :                                 Open Boundary Condition
(   12)   !!----------------------------------------------------------------------
(   13)   !! history :
(   14)   !!  8.0   01/91   (CLIPPER)  Original code
(   15)   !!  8.5   06/02   (C. Talandier)  modules
(   16)   !!        06/04   (F. Durand) Dimensions of arrays vsdta, tsdta, ssdta,
(   17)   !!                vndta, tndta, sndta, uwdta, twdta, swdta, uedta, tedta, sedta
(   18)   !!                are defined to the actual dimensions of the OBs i.e.
(   19)   !!		     (jpisd:jpisf,jpk,jptobc) for the South OB
(   20)   !! 		     (jpind:jpinf,jpk,jptobc) for the North OB
(   21)   !!		     (jpjwd:jpjwf,jpk,jptobc) for the West OB
(   22)   !! 		     (jpjed:jpjef,jpk,jptobc) for the East OB
(   23)   !!
(   24)   !! 2007.12.07      D.Brickman 
(   25)   !!                 Dimensions as listed above are potentially incorrect
(   26)   !!                 as looping indices in the routines that use these arrays
(   27)   !!                 can start at a number less than the starting index
(   28)   !!                 e.g. jpjed = 2, but loop can start at 1 ----> ERROR
(   29)   !!                 Solution: set starting index = 1 
(   30)   !!		     (1:jpisf,jpk,jptobc) for the South OB
(   31)   !! 		     (1:jpinf,jpk,jptobc) for the North OB 
(   32)   !!		     (1:jpjwf,jpk,jptobc) for the West OB
(   33)   !! 		     (1:jpjef,jpk,jptobc) for the East OB 
(   34)   !!
(   35)   !!----------------------------------------------------------------------
(   36)   !! * Modules used
(   37)   USE par_oce         ! ocean parameters
(   38)   USE obc_par         ! open boundary condition parameters
(   39) 
(   40)   IMPLICIT NONE
(   41)   PUBLIC
(   42)   !!----------------------------------------------------------------------
(   43)   !!  OPA 9.0 , LOCEAN-IPSL (2005)
(   44)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_oce.F90,v 1.5 2005/12/28 09:25:07 opalod Exp $
(   45)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   46)   !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 2

(   47) 
(   48)   !!----------------------------------------------------------------------
(   49)   !! open boundary variables
(   50)   !!----------------------------------------------------------------------
(   51)   !!
(   52)   !!General variables for open boundaries:
(   53)   !!--------------------------------------
(   54)   INTEGER ::              & !: * namelist ??? *
(   55)      nbobc    = 3  ,      & !: number of open boundaries ( 1=< nbobc =< 4 ) !byoung
(   56)      nobc_dta = 0           !:  = 0 use the initial state as obc data
(   57)      !                      !   = 1 read obc data in obcxxx.dta files
(   58) 
(   59)   LOGICAL ::  ln_obc_clim = .true.  !:  obc data files are climatological
(   60)   LOGICAL ::  ln_obc_fla  = .false. !:  Flather open boundary condition not used
(   61)   LOGICAL ::  ln_vol_cst  = .true.  !:  Conservation of the whole volume
(   62) 
(   63)   REAL(wp) ::             & !!: open boundary namelist (namobc)
(   64)      rdpein =  1.  ,      &  !: damping time scale for inflow at East open boundary
(   65)      rdpwin =  1.  ,      &  !:    "                      "   at West open boundary
(   66)      rdpsin =  1.  ,      &  !:    "                      "   at South open boundary
(   67)      rdpnin =  1.  ,      &  !:    "                      "   at North open boundary
(   68)      rdpeob = 15.  ,      &  !: damping time scale for the climatology at East open boundary
(   69)      rdpwob = 15.  ,      &  !:    "                           "       at West open boundary
(   70)      rdpsob = 15.  ,      &  !:    "                           "       at South open boundary
(   71)      rdpnob = 15.  ,      &  !:    "                           "       at North open boundary
(   72)      volemp =  1.            !: = 0 the total volume will have the variability of the
(   73)                              !      surface Flux E-P else (volemp = 1) the volume will be constant
(   74)                              !  = 1 the volume will be constant during all the integration.
(   75) 
(   76)   LOGICAL ::              &  !:
(   77)      lfbceast, lfbcwest,  &  !: logical flag for a fixed East and West open boundaries
(   78)      lfbcnorth, lfbcsouth    !: logical flag for a fixed North and South open boundaries
(   79)      !                       !  These logical flags are set to 'true' if damping time
(   80)      !                       !  scale are set to 0 in the namelist, for both inflow and outflow).
(   81) 
(   82)   REAL(wp), DIMENSION(jpi,jpj) :: &  !:
(   83)      obctmsk                !: mask array identical to tmask, execpt along OBC where it is set to 0
(   84)      !                      !  it used to calculate the cumulate flux E-P in the obcvol.F90 routine
(   85) 
(   86)   !!----------------
(   87)   !! Rigid lid case:
(   88)   !!----------------
(   89)   INTEGER ::   nbic !: number of isolated coastlines ( 0 <= nbic <= 3 )
(   90) 
(   91)   INTEGER, DIMENSION(jpnic,0:4,3) ::   &  !:
(   92)      miic, mjic     !: position of isolated coastlines points
(   93) 
(   94)   INTEGER, DIMENSION(0:4,3) ::   &  !:
(   95)      mnic           !: number of points on isolated coastlines
(   96) 
(   97)   REAL(wp), DIMENSION(jpi,jpj) ::   &  !:
(   98)      gcbob          !: right hand side of the barotropic elliptic equation associated
(   99)      !              !  with the OBC
(  100) 
(  101)   REAL(wp), DIMENSION(jpi,jpj,3) ::   &  !:
(  102)      gcfobc         !: coef. associated with the contribution of isolated coastlines
(  103)      !              !  to the right hand side of the barotropic elliptic equation
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 3

(  105)   REAL(wp), DIMENSION(3) ::   &  !:
(  106)      gcbic          !: time variation of the barotropic stream function along the
(  107)      !              !  isolated coastlines
(  108) 
(  109)   REAL(wp), DIMENSION(1) ::   &  !:
(  110)      bsfic0         !: barotropic stream function on isolated coastline
(  111) 
(  112)   REAL(wp), DIMENSION(3) ::   &  !:
(  113)      bsfic          !: barotropic stream function on isolated coastline
(  114) 
(  115) !!DB 2008.05.06 global masking arrays for obcdta calcs
(  116) !!Arrays computed in obcini
(  117) !!del the first one later
(  118)   REAL(wp) :: emaskg(jpjdta)=0.0, wmaskg(jpjdta)=0.0, smaskg(jpidta)=0.0, nmaskg(jpidta)=0.0
(  119)   REAL(wp) :: emaskg2(jpjdta,jpk)=0.0, wmaskg2(jpjdta,jpk)=0.0, & 
(  120)        smaskg2(jpidta,jpk)=0.0, nmaskg2(jpidta,jpk)=0.0
(  121) 
(  122) 
(  123) 
(  124)   !!--------------------
(  125)   !! East open boundary:
(  126)   !!--------------------
(  127)   INTEGER ::   nie0  , nie1      !: do loop index in mpp case for jpieob
(  128)   INTEGER ::   nie0p1, nie1p1    !: do loop index in mpp case for jpieob+1
(  129)   INTEGER ::   nie0m1, nie1m1    !: do loop index in mpp case for jpieob-1
(  130)   INTEGER ::   nje0  , nje1      !: do loop index in mpp case for jpjed, jpjef
(  131)   INTEGER ::   nje0p1, nje1m1    !: do loop index in mpp case for jpjedp1,jpjefm1
(  132)   INTEGER ::   nje1m2, nje0m1    !: do loop index in mpp case for jpjefm1-1,jpjed
(  133) 
(  134)   REAL(wp), DIMENSION(jpj) ::    &  !:
(  135)      bsfeob              !: now barotropic stream fuction computed at the OBC. The corres-
(  136)      !                   ! ponding bsfn will be computed by the forward time step in dynspg.
(  137) 
(  138)   REAL(wp), DIMENSION(jpj,3,3) ::   &  !:
(  139)      bebnd               !: east boundary barotropic streamfunction over 3 rows
(  140)      !                   ! and 3 time step (now, before, and before before)
(  141) 
(  142) !DB
(  143)   REAL(wp), DIMENSION(1:jpjef) ::   &  !:
(  144)      bfoe,             & !: now climatology of the east boundary barotropic stream function
(  145)      sshfoe,           & !: now climatology of the east boundary sea surface height
(  146)      eta_e,      & !: AD3Aug07:  elev. geostropic elev computed from uedta
(  147)      ubtfoe,vbtfoe       !: now climatology of the east boundary barotropic transport
(  148) 
(  149)   REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
(  150)      ufoe, vfoe,       & !: now climatology of the east boundary velocities
(  151) !     tfoe, sfoe,       & !: now climatology of the east boundary temperature and salinity
(  152)      uclie               !: baroclinic componant of the zonal velocity after radiation
(  153)      !                   ! in the obcdyn.F90 routine
(  154) !byoung
(  155)   REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
(  156)      tfoe, sfoe        !: now climatology of the east boundary temperature and salinity
(  157)      
(  158) !DB jpjed <---- 1
(  159)   REAL(wp), DIMENSION(1:jpjef,jpj) ::   &  !:
(  160)      sshfoe_b            !: east boundary ssh correction averaged over the barotropic loop
(  161)                          !: (if Flather's algoritm applied at open boundary)
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 4

(  163) !DB jpjed <---- 1
(  164)   REAL(wp), DIMENSION(1:jpjef,0:jptobc+1) ::   &  !:
(  165)      sshedta, ubtedta    !: array used for interpolating monthly data on the east boundary
(  166) 
(  167) !DB jpjed <---- 1
(  168)   REAL(wp), DIMENSION(1:jpjef,5) ::    &    !:
(  169)      tidesshemag, tidesshepha, tidevbtemag, tidevbtepha    
(  170) 
(  171) !DB jpjed <---- 1
(  172)   REAL(wp), DIMENSION(1:jpjef,jpk,jptobc) ::   &  !:
(  173)      uedta 
(  174) 
(  175) !AD
(  176) !DB jpjed <---- 1
(  177)   REAL(wp), DIMENSION(1:jpjef,jpk) ::   &  !:
(  178)      uedta1 !: array used for interpolating monthly data on the east boundary
(  179) 
(  180) !DB jpjed <---- 1
(  181)   REAL(wp), DIMENSION(1:jpjef,jpk,5,jptobc) ::   &  !:
(  182)      tedta, sedta, tedta1, sedta1 !: array used for interpolating monthly data on the east boundary
(  183) 
(  184)   !!-------------------------------
(  185)   !! Arrays for radiative East OBC:
(  186)   !!-------------------------------
(  187)   REAL(wp), DIMENSION(jpj,jpk,3,3) ::   &  !:
(  188)      uebnd, vebnd                  !: baroclinic u & v component of the velocity over 3 rows
(  189)                                    ! and 3 time step (now, before, and before before)
(  190) 
(  191)   REAL(wp), DIMENSION(jpj,jpk,2,2) ::   &  !:
(  192)      tebnd, sebnd                  !: East boundary temperature and salinity over 2 rows
(  193)                                    ! and 2 time step (now and before)
(  194) 
(  195)   REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
(  196)      u_cxebnd, v_cxebnd            !: Zonal component of the phase speed ratio computed with
(  197)                                    ! radiation of u and v velocity (respectively) at the
(  198)                                    ! east open boundary (u_cxebnd = cx rdt )
(  199) 
(  200)   REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
(  201)      uemsk, vemsk, temsk           !: 2D mask for the East OB
(  202) 
(  203) !byoung
(  204)   REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
(  205)      temsk5
(  206) 
(  207)   ! Note that those arrays are optimized for mpp case
(  208)   ! (hence the dimension jpj is the size of one processor subdomain)
(  209) 
(  210)   !!--------------------
(  211)   !! West open boundary
(  212)   !!--------------------
(  213)   INTEGER ::   niw0  , niw1       !: do loop index in mpp case for jpiwob
(  214)   INTEGER ::   niw0p1, niw1p1     !: do loop index in mpp case for jpiwob+1
(  215)   INTEGER ::   njw0  , njw1       !: do loop index in mpp case for jpjwd, jpjwf
(  216)   INTEGER ::   njw0p1, njw1m1     !: do loop index in mpp case for jpjwdp1,jpjwfm1
(  217)   INTEGER ::   njw1m2, njw0m1     !: do loop index in mpp case for jpjwfm2,jpjwd
(  218) 
(  219)   REAL(wp), DIMENSION(jpj) ::   &  !:
(  220)      bsfwob              !: now barotropic stream fuction computed at the OBC. The corres-






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 5

(  221)      !                   !  ponding bsfn will be computed by the forward time step in dynspg.
(  222) 
(  223)   REAL(wp), DIMENSION(jpj,3,3) ::   &  !:
(  224)      bwbnd               !: West boundary barotropic streamfunction over
(  225)      !                   !  3 rows and 3 time step (now, before, and before before)
(  226) 
(  227) !DB
(  228)   REAL(wp), DIMENSION(1:jpjwf) ::   &  !:
(  229)      bfow,             & !: now climatology of the west boundary barotropic stream function
(  230)      sshfow,           & !: now climatology of the west boundary sea surface height
(  231)      eta_w,       & !: elevation computed from uwdta
(  232)      ubtfow,vbtfow       !: now climatology of the west boundary barotropic transport
(  233)      
(  234) 
(  235)   REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
(  236)      ufow, vfow,       & !: now climatology of the west velocities
(  237) !     tfow, sfow,       & !: now climatology of the west temperature and salinity
(  238)      ucliw               !: baroclinic componant of the zonal velocity after the radiation
(  239)      !                   !  in the obcdyn.F90 routine
(  240) 
(  241) !byoung
(  242)   REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
(  243)      tfow, sfow        !: now climatology of the west temperature and salinity
(  244)      
(  245)   REAL(wp), DIMENSION(1:jpjwf,jpj) ::   &  !:
(  246)      sshfow_b            !: west boundary ssh correction averaged over the barotropic loop
(  247)                          !: (if Flather's algoritm applied at open boundary)
(  248) 
(  249)   REAL(wp), DIMENSION(1:jpjwf,0:jptobc+1) ::   &  !:
(  250)      sshwdta, ubtwdta    !: array used for interpolating monthly data on the west boundary
(  251) 
(  252) !ylu
(  253) !  REAL(wp), DIMENSION(1:jpjwf,4) ::    &    !:
(  254)   REAL(wp), DIMENSION(1:jpjwf,5) ::    &    !:
(  255)      tidesshwmag, tidesshwpha, tidevbtwmag, tidevbtwpha    
(  256) 
(  257)   REAL(wp), DIMENSION(1:jpjwf,jpk,jptobc) ::   &  !:
(  258)      uwdta !, uwdta1 !: array used for interpolating monthly data on the west boundary
(  259) !     uwdta, twdta, swdta !: array used for interpolating monthly data on the west boundary
(  260) 
(  261) !DBG
(  262)   REAL(wp), DIMENSION(1:jpjwf,jpk) ::   &  !:
(  263)      uwdta1 !: array used for interpolating monthly data on the west boundary
(  264) 
(  265) 
(  266) !byoung
(  267)   REAL(wp), DIMENSION(1:jpjwf,jpk,5,jptobc) ::   &  !:
(  268)      twdta, swdta, twdta1, swdta1 !: array used for interpolating monthly data on the west boundary
(  269) 
(  270)   !!-------------------------------
(  271)   !! Arrays for radiative West OBC
(  272)   !!-------------------------------
(  273)   REAL(wp), DIMENSION(jpj,jpk,3,3) ::   &  !:
(  274)      uwbnd, vwbnd                  !: baroclinic u & v components of the velocity over 3 rows
(  275)      !                             !  and 3 time step (now, before, and before before)
(  276) 
(  277)   REAL(wp), DIMENSION(jpj,jpk,2,2) ::   &  !:
(  278)      twbnd, swbnd                  !: west boundary temperature and salinity over 2 rows and






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 6

(  279)      !                             !  2 time step (now and before)
(  280) 
(  281)   REAL(wp), DIMENSION(jpj,jpk) ::    &  !:
(  282)      u_cxwbnd, v_cxwbnd            !: Zonal component of the phase speed ratio computed with
(  283)      !                             !  radiation of zonal and meridional velocity (respectively)
(  284)      !                             !  at the west open boundary (u_cxwbnd = cx rdt )
(  285) 
(  286)   REAL(wp), DIMENSION(jpj,jpk) ::    &  !:
(  287)      uwmsk, vwmsk, twmsk           !: 2D mask for the West OB
(  288) 
(  289) !byoung
(  290)   REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
(  291)      twmsk5
(  292) 
(  293)   ! Note that those arrays are optimized for mpp case
(  294)   ! (hence the dimension jpj is the size of one processor subdomain)
(  295) 
(  296)   !!---------------------
(  297)   !! North open boundary
(  298)   !!---------------------
(  299)   INTEGER ::   nin0  , nin1       !: do loop index in mpp case for jpind, jpinf
(  300)   INTEGER ::   nin0p1, nin1m1     !: do loop index in mpp case for jpindp1, jpinfm1
(  301)   INTEGER ::   nin1m2, nin0m1     !: do loop index in mpp case for jpinfm1-1,jpind
(  302)   INTEGER ::   njn0  , njn1       !: do loop index in mpp case for jpnob
(  303)   INTEGER ::   njn0p1, njn1p1     !: do loop index in mpp case for jpnob+1
(  304)   INTEGER ::   njn0m1, njn1m1     !: do loop index in mpp case for jpnob-1
(  305) 
(  306)   REAL(wp), DIMENSION(jpi) ::   &  !:
(  307)      bsfnob              !: now barotropic stream fuction computed at the OBC. The corres-
(  308)      !                   !  ponding bsfn will be computed by the forward time step in dynspg.
(  309) 
(  310)   REAL(wp), DIMENSION(jpi,3,3) ::   &  !:
(  311)      bnbnd               !: north boundary barotropic streamfunction over
(  312)      !                   !  3 rows and 3 time step (now, before, and before before)
(  313) 
(  314)   REAL(wp), DIMENSION(1:jpinf) ::   &  !:
(  315)      bfon,             & !: now climatology of the north boundary barotropic stream function
(  316)      sshfon,           & !: now climatology of the north boundary sea surface height
(  317)      ubtfon,vbtfon       !: now climatology of the north boundary barotropic transport
(  318) 
(  319)   REAL(wp), DIMENSION(jpi,jpk) ::   &    !:
(  320)      ufon, vfon,       & !: now climatology of the north boundary velocities
(  321) !      tfon, sfon,       & !: now climatology of the north boundary temperature and salinity
(  322)      vclin               !: baroclinic componant of the meridian velocity after the radiation
(  323)      !                   !  in yhe obcdyn.F90 routine
(  324) 
(  325) !sujie
(  326)   REAL(wp), DIMENSION(jpi,jpk,5) ::   &    !:
(  327)      tfon, sfon        !: now climatology of the north boundary temperature and salinity
(  328) 
(  329)   REAL(wp), DIMENSION(1:jpinf,jpj) ::   &  !:
(  330)      sshfon_b            !: north boundary ssh correction averaged over the barotropic loop
(  331)                          !: (if Flather's algoritm applied at open boundary)
(  332) 
(  333)   REAL(wp), DIMENSION(1:jpinf,0:jptobc+1) ::   &  !:
(  334)      sshndta, vbtndta   !: array used for interpolating monthly data on the north boundary
(  335) 
(  336) !ylu






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 7

(  337)   REAL(wp), DIMENSION(1:jpinf,5) ::    &    !:
(  338)      tidesshnmag, tidesshnpha, tidevbtnmag, tidevbtnpha    
(  339) 
(  340)   REAL(wp), DIMENSION(1:jpinf,jpk,jptobc) ::   &  !:
(  341)      vndta !sujie, tndta, sndta !: array used for interpolating monthly data on the north boundary
(  342) !sujie
(  343)   REAL(wp), DIMENSION(1:jpinf,jpk,5,jptobc) ::   &  !:
(  344)      tndta, sndta !: array used for interpolating monthly data on the north boundary
(  345) 
(  346)   !!--------------------------------
(  347)   !! Arrays for radiative North OBC
(  348)   !!--------------------------------
(  349)   !!
(  350)   REAL(wp), DIMENSION(jpi,jpk,3,3) ::   &   !:
(  351)      unbnd, vnbnd                  !: baroclinic u & v components of the velocity over 3
(  352)      !                             !  rows and 3 time step (now, before, and before before)
(  353) 
(  354)   REAL(wp), DIMENSION(jpi,jpk,2,2) ::   &   !:
(  355)      tnbnd, snbnd                  !: north boundary temperature and salinity over
(  356)      !                             !  2 rows and 2 time step (now and before)
(  357) 
(  358)   REAL(wp), DIMENSION(jpi,jpk) ::   &     !:
(  359)      u_cynbnd, v_cynbnd            !: Meridional component of the phase speed ratio compu-
(  360)      !                             !  ted with radiation of zonal and meridional velocity
(  361)      !                             !  (respectively) at the north OB (u_cynbnd = cx rdt )
(  362) 
(  363)   REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
(  364)      unmsk, vnmsk, tnmsk           !: 2D mask for the North OB
(  365) 
(  366) !sujie
(  367)   REAL(wp), DIMENSION(jpi,jpk,5) ::   &  !:
(  368)      tnmsk5
(  369) 
(  370)   ! Note that those arrays are optimized for mpp case
(  371)   ! (hence the dimension jpj is the size of one processor subdomain)
(  372) 
(  373)   !!---------------------
(  374)   !! South open boundary
(  375)   !!---------------------
(  376)   INTEGER ::   nis0  , nis1       !: do loop index in mpp case for jpisd, jpisf
(  377)   INTEGER ::   nis0p1, nis1m1     !: do loop index in mpp case for jpisdp1, jpisfm1
(  378)   INTEGER ::   nis1m2, nis0m1     !: do loop index in mpp case for jpisfm1-1,jpisd
(  379)   INTEGER ::   njs0  , njs1       !: do loop index in mpp case for jpsob
(  380)   INTEGER ::   njs0p1, njs1p1     !: do loop index in mpp case for jpsob+1
(  381) 
(  382)   REAL(wp), DIMENSION(jpi) ::    &   !:
(  383)      bsfsob              !: now barotropic stream fuction computed at the OBC.The corres-
(  384)      !                   !  ponding bsfn will be computed by the forward time step in dynspg.
(  385)   REAL(wp), DIMENSION(jpi,3,3) ::   &   !:
(  386)      bsbnd               !: south boundary barotropic stream function over
(  387)      !                   !  3 rows and 3 time step (now, before, and before before)
(  388) 
(  389) !DB
(  390)   REAL(wp), DIMENSION(1:jpisf) ::    &   !:
(  391)      bfos,             & !: now climatology of the south boundary barotropic stream function
(  392)      sshfos,           & !: now climatology of the south boundary sea surface height
(  393)      eta_s,            & !: AD3Aug07:  elev. geostropic elev computed from vsdta
(  394)      ubtfos,vbtfos       !: now climatology of the south boundary barotropic transport






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 8

(  395) 
(  396)   REAL(wp), DIMENSION(jpi,jpk) ::    &   !:
(  397)      ufos, vfos ,       & !: now climatology of the south boundary velocities
(  398) !      tfos, sfos,       & !: now climatology of the south boundary temperature and salinity
(  399)      vclis               !: baroclinic componant of the meridian velocity after the radiation
(  400)      !                   !  in the obcdyn.F90 routine
(  401) 
(  402) !sujie
(  403)   REAL(wp), DIMENSION(jpi,jpk,5) ::    &   !:
(  404)      tfos, sfos        !: now climatology of the south boundary temperature and salinity
(  405) 
(  406) !DB
(  407)   REAL(wp), DIMENSION(1:jpisf,jpj) ::   &  !:
(  408)      sshfos_b            !: south boundary ssh correction averaged over the barotropic loop
(  409)                          !: (if Flather's algoritm applied at open boundary)
(  410) 
(  411)   REAL(wp), DIMENSION(1:jpisf,0:jptobc+1) ::    &    !:
(  412)      sshsdta, vbtsdta    !: array used for interpolating monthly data on the south boundary
(  413) 
(  414) !ylu
(  415)   REAL(wp), DIMENSION(1:jpisf,5) ::    &    !:
(  416)      tidesshsmag, tidesshspha, tidevbtsmag, tidevbtspha    
(  417) 
(  418)   REAL(wp), DIMENSION(1:jpisf,jpk,jptobc) ::    &    !:
(  419)      vsdta !sujie, tsdta, ssdta   !: array used for interpolating monthly data on the south boundary
(  420) 
(  421) !DB
(  422)   REAL(wp), DIMENSION(1:jpisf,jpk) ::    &    !:
(  423)      vsdta1  !AD: array used for interpolatated global obc array
(  424)      
(  425) !sujie
(  426)   REAL(wp), DIMENSION(1:jpisf,jpk,5,jptobc) ::    &    !:
(  427)      tsdta, ssdta, tsdta1, ssdta1   !: array used for interpolating monthly data on the south boundary
(  428) 
(  429)   !!--------------------------------
(  430)   !! Arrays for radiative South OBC
(  431)   !!--------------------------------
(  432)   !!                        computed by the forward time step in dynspg.
(  433)   REAL(wp), DIMENSION(jpi,jpk,3,3) ::   &   !:
(  434)      usbnd, vsbnd                  !: baroclinic u & v components of the velocity over 3
(  435)      !                             !  rows and 3 time step (now, before, and before before)
(  436) 
(  437)   REAL(wp), DIMENSION(jpi,jpk,2,2) ::   &  !:
(  438)      tsbnd, ssbnd                  !: south boundary temperature and salinity over
(  439)      !                             !  2 rows and 2 time step (now and before)
(  440) 
(  441)   REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
(  442)      u_cysbnd, v_cysbnd            !: Meridional component of the phase speed ratio compu-
(  443)      !                             !  ted with radiation of zonal and meridional velocity
(  444)      !                             !  (repsectively) at the south OB (u_cynbnd = cx rdt )
(  445) 
(  446)   REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
(  447)      usmsk, vsmsk, tsmsk           !: 2D mask for the South OB
(  448) 
(  449) !sujie
(  450)   REAL(wp), DIMENSION(jpi,jpk,5) ::   &  !:
(  451)      tsmsk5
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 9

(  453)   ! Note that those arrays are optimized for mpp case
(  454)   ! (hence the dimension jpj is the size of one processor subdomain)
(  455) 
(  456) #else
(  457)   !!----------------------------------------------------------------------
(  458)   !!   Default option :                                       Empty module
(  459)   !!----------------------------------------------------------------------
(  460) #endif
(  461) 
(  462)   !!======================================================================
(  463) END MODULE obc_oce





















































PGF90 (Version     10.5)          09/23/2010  09:49:16      page 10

(    1) # 1 "obc_oce.F90"
(    1) !!DB 2008.04.14
(    2) !!see obc_oce.Q90 for last version of this program
(    3) 
(    4) MODULE obc_oce
(    5)   !!==============================================================================
(    6)   !!                       ***  MODULE obc_oce   ***
(    7)   !! Open Boundary Cond. :   define related variables
(    8)   !!==============================================================================
(    9) # 10
(   10)   !!----------------------------------------------------------------------
(   11)   !!   'key_obc' :                                 Open Boundary Condition
(   12)   !!----------------------------------------------------------------------
(   13)   !! history :
(   14)   !!  8.0   01/91   (CLIPPER)  Original code
(   15)   !!  8.5   06/02   (C. Talandier)  modules
(   16)   !!        06/04   (F. Durand) Dimensions of arrays vsdta, tsdta, ssdta,
(   17)   !!                vndta, tndta, sndta, uwdta, twdta, swdta, uedta, tedta, sedta
(   18)   !!                are defined to the actual dimensions of the OBs i.e.
(   19)   !!		     (jpisd:jpisf,jpk,jptobc) for the South OB
(   20)   !! 		     (jpind:jpinf,jpk,jptobc) for the North OB
(   21)   !!		     (jpjwd:jpjwf,jpk,jptobc) for the West OB
(   22)   !! 		     (jpjed:jpjef,jpk,jptobc) for the East OB
(   23)   !!
(   24)   !! 2007.12.07      D.Brickman 
(   25)   !!                 Dimensions as listed above are potentially incorrect
(   26)   !!                 as looping indices in the routines that use these arrays
(   27)   !!                 can start at a number less than the starting index
(   28)   !!                 e.g. jpjed = 2, but loop can start at 1 ----> ERROR
(   29)   !!                 Solution: set starting index = 1 
(   30)   !!		     (1:jpisf,jpk,jptobc) for the South OB
(   31)   !! 		     (1:jpinf,jpk,jptobc) for the North OB 
(   32)   !!		     (1:jpjwf,jpk,jptobc) for the West OB
(   33)   !! 		     (1:jpjef,jpk,jptobc) for the East OB 
(   34)   !!
(   35)   !!----------------------------------------------------------------------
(   36)   !! * Modules used
(   37)   USE par_oce         ! ocean parameters
(   38)   USE obc_par         ! open boundary condition parameters
(   39) 
(   40)   IMPLICIT NONE
(   41)   PUBLIC
(   42)   !!----------------------------------------------------------------------
(   43)   !!  OPA 9.0 , LOCEAN-IPSL (2005)
(   44)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_oce.F90,v 1.5 2005/12/28 09:25:07 opalod Exp $
(   45)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   46)   !!----------------------------------------------------------------------
(   47) 
(   48)   !!----------------------------------------------------------------------
(   49)   !! open boundary variables
(   50)   !!----------------------------------------------------------------------
(   51)   !!
(   52)   !!General variables for open boundaries:
(   53)   !!--------------------------------------
(   54)   INTEGER ::              & !: * namelist ??? *
(   55)      nbobc    = 3  ,      & !: number of open boundaries ( 1=< nbobc =< 4 ) !byoung
(   56)      nobc_dta = 0           !:  = 0 use the initial state as obc data
(   57)      !                      !   = 1 read obc data in obcxxx.dta files






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 11

(   58) 
(   59)   LOGICAL ::  ln_obc_clim = .true.  !:  obc data files are climatological
(   60)   LOGICAL ::  ln_obc_fla  = .false. !:  Flather open boundary condition not used
(   61)   LOGICAL ::  ln_vol_cst  = .true.  !:  Conservation of the whole volume
(   62) 
(   63)   REAL(wp) ::             & !!: open boundary namelist (namobc)
(   64)      rdpein =  1.  ,      &  !: damping time scale for inflow at East open boundary
(   65)      rdpwin =  1.  ,      &  !:    "                      "   at West open boundary
(   66)      rdpsin =  1.  ,      &  !:    "                      "   at South open boundary
(   67)      rdpnin =  1.  ,      &  !:    "                      "   at North open boundary
(   68)      rdpeob = 15.  ,      &  !: damping time scale for the climatology at East open boundary
(   69)      rdpwob = 15.  ,      &  !:    "                           "       at West open boundary
(   70)      rdpsob = 15.  ,      &  !:    "                           "       at South open boundary
(   71)      rdpnob = 15.  ,      &  !:    "                           "       at North open boundary
(   72)      volemp =  1.            !: = 0 the total volume will have the variability of the
(   73)                              !      surface Flux E-P else (volemp = 1) the volume will be constant
(   74)                              !  = 1 the volume will be constant during all the integration.
(   75) 
(   76)   LOGICAL ::              &  !:
(   77)      lfbceast, lfbcwest,  &  !: logical flag for a fixed East and West open boundaries
(   78)      lfbcnorth, lfbcsouth    !: logical flag for a fixed North and South open boundaries
(   79)      !                       !  These logical flags are set to 'true' if damping time
(   80)      !                       !  scale are set to 0 in the namelist, for both inflow and outflow).
(   81) 
(   82)   REAL(wp), DIMENSION(jpi,jpj) :: &  !:
(   83)      obctmsk                !: mask array identical to tmask, execpt along OBC where it is set to 0
(   84)      !                      !  it used to calculate the cumulate flux E-P in the obcvol.F90 routine
(   85) 
(   86)   !!----------------
(   87)   !! Rigid lid case:
(   88)   !!----------------
(   89)   INTEGER ::   nbic !: number of isolated coastlines ( 0 <= nbic <= 3 )
(   90) 
(   91)   INTEGER, DIMENSION(jpnic,0:4,3) ::   &  !:
(   92)      miic, mjic     !: position of isolated coastlines points
(   93) 
(   94)   INTEGER, DIMENSION(0:4,3) ::   &  !:
(   95)      mnic           !: number of points on isolated coastlines
(   96) 
(   97)   REAL(wp), DIMENSION(jpi,jpj) ::   &  !:
(   98)      gcbob          !: right hand side of the barotropic elliptic equation associated
(   99)      !              !  with the OBC
(  100) 
(  101)   REAL(wp), DIMENSION(jpi,jpj,3) ::   &  !:
(  102)      gcfobc         !: coef. associated with the contribution of isolated coastlines
(  103)      !              !  to the right hand side of the barotropic elliptic equation
(  104) 
(  105)   REAL(wp), DIMENSION(3) ::   &  !:
(  106)      gcbic          !: time variation of the barotropic stream function along the
(  107)      !              !  isolated coastlines
(  108) 
(  109)   REAL(wp), DIMENSION(1) ::   &  !:
(  110)      bsfic0         !: barotropic stream function on isolated coastline
(  111) 
(  112)   REAL(wp), DIMENSION(3) ::   &  !:
(  113)      bsfic          !: barotropic stream function on isolated coastline
(  114) 
(  115) !!DB 2008.05.06 global masking arrays for obcdta calcs






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 12

(  116) !!Arrays computed in obcini
(  117) !!del the first one later
(  118)   REAL(wp) :: emaskg(jpjdta)=0.0, wmaskg(jpjdta)=0.0, smaskg(jpidta)=0.0, nmaskg(jpidta)=0.0
(  119)   REAL(wp) :: emaskg2(jpjdta,jpk)=0.0, wmaskg2(jpjdta,jpk)=0.0, & 
(  120)        smaskg2(jpidta,jpk)=0.0, nmaskg2(jpidta,jpk)=0.0
(  121) 
(  122) 
(  123) 
(  124)   !!--------------------
(  125)   !! East open boundary:
(  126)   !!--------------------
(  127)   INTEGER ::   nie0  , nie1      !: do loop index in mpp case for jpieob
(  128)   INTEGER ::   nie0p1, nie1p1    !: do loop index in mpp case for jpieob+1
(  129)   INTEGER ::   nie0m1, nie1m1    !: do loop index in mpp case for jpieob-1
(  130)   INTEGER ::   nje0  , nje1      !: do loop index in mpp case for jpjed, jpjef
(  131)   INTEGER ::   nje0p1, nje1m1    !: do loop index in mpp case for jpjedp1,jpjefm1
(  132)   INTEGER ::   nje1m2, nje0m1    !: do loop index in mpp case for jpjefm1-1,jpjed
(  133) 
(  134)   REAL(wp), DIMENSION(jpj) ::    &  !:
(  135)      bsfeob              !: now barotropic stream fuction computed at the OBC. The corres-
(  136)      !                   ! ponding bsfn will be computed by the forward time step in dynspg.
(  137) 
(  138)   REAL(wp), DIMENSION(jpj,3,3) ::   &  !:
(  139)      bebnd               !: east boundary barotropic streamfunction over 3 rows
(  140)      !                   ! and 3 time step (now, before, and before before)
(  141) 
(  142) !DB
(  143)   REAL(wp), DIMENSION(1:jpjef) ::   &  !:
(  144)      bfoe,             & !: now climatology of the east boundary barotropic stream function
(  145)      sshfoe,           & !: now climatology of the east boundary sea surface height
(  146)      eta_e,      & !: AD3Aug07:  elev. geostropic elev computed from uedta
(  147)      ubtfoe,vbtfoe       !: now climatology of the east boundary barotropic transport
(  148) 
(  149)   REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
(  150)      ufoe, vfoe,       & !: now climatology of the east boundary velocities
(  151) !     tfoe, sfoe,       & !: now climatology of the east boundary temperature and salinity
(  152)      uclie               !: baroclinic componant of the zonal velocity after radiation
(  153)      !                   ! in the obcdyn.F90 routine
(  154) !byoung
(  155)   REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
(  156)      tfoe, sfoe        !: now climatology of the east boundary temperature and salinity
(  157)      
(  158) !DB jpjed <---- 1
(  159)   REAL(wp), DIMENSION(1:jpjef,jpj) ::   &  !:
(  160)      sshfoe_b            !: east boundary ssh correction averaged over the barotropic loop
(  161)                          !: (if Flather's algoritm applied at open boundary)
(  162) 
(  163) !DB jpjed <---- 1
(  164)   REAL(wp), DIMENSION(1:jpjef,0:jptobc+1) ::   &  !:
(  165)      sshedta, ubtedta    !: array used for interpolating monthly data on the east boundary
(  166) 
(  167) !DB jpjed <---- 1
(  168)   REAL(wp), DIMENSION(1:jpjef,5) ::    &    !:
(  169)      tidesshemag, tidesshepha, tidevbtemag, tidevbtepha    
(  170) 
(  171) !DB jpjed <---- 1
(  172)   REAL(wp), DIMENSION(1:jpjef,jpk,jptobc) ::   &  !:
(  173)      uedta 






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 13

(  174) 
(  175) !AD
(  176) !DB jpjed <---- 1
(  177)   REAL(wp), DIMENSION(1:jpjef,jpk) ::   &  !:
(  178)      uedta1 !: array used for interpolating monthly data on the east boundary
(  179) 
(  180) !DB jpjed <---- 1
(  181)   REAL(wp), DIMENSION(1:jpjef,jpk,5,jptobc) ::   &  !:
(  182)      tedta, sedta, tedta1, sedta1 !: array used for interpolating monthly data on the east boundary
(  183) 
(  184)   !!-------------------------------
(  185)   !! Arrays for radiative East OBC:
(  186)   !!-------------------------------
(  187)   REAL(wp), DIMENSION(jpj,jpk,3,3) ::   &  !:
(  188)      uebnd, vebnd                  !: baroclinic u & v component of the velocity over 3 rows
(  189)                                    ! and 3 time step (now, before, and before before)
(  190) 
(  191)   REAL(wp), DIMENSION(jpj,jpk,2,2) ::   &  !:
(  192)      tebnd, sebnd                  !: East boundary temperature and salinity over 2 rows
(  193)                                    ! and 2 time step (now and before)
(  194) 
(  195)   REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
(  196)      u_cxebnd, v_cxebnd            !: Zonal component of the phase speed ratio computed with
(  197)                                    ! radiation of u and v velocity (respectively) at the
(  198)                                    ! east open boundary (u_cxebnd = cx rdt )
(  199) 
(  200)   REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
(  201)      uemsk, vemsk, temsk           !: 2D mask for the East OB
(  202) 
(  203) !byoung
(  204)   REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
(  205)      temsk5
(  206) 
(  207)   ! Note that those arrays are optimized for mpp case
(  208)   ! (hence the dimension jpj is the size of one processor subdomain)
(  209) 
(  210)   !!--------------------
(  211)   !! West open boundary
(  212)   !!--------------------
(  213)   INTEGER ::   niw0  , niw1       !: do loop index in mpp case for jpiwob
(  214)   INTEGER ::   niw0p1, niw1p1     !: do loop index in mpp case for jpiwob+1
(  215)   INTEGER ::   njw0  , njw1       !: do loop index in mpp case for jpjwd, jpjwf
(  216)   INTEGER ::   njw0p1, njw1m1     !: do loop index in mpp case for jpjwdp1,jpjwfm1
(  217)   INTEGER ::   njw1m2, njw0m1     !: do loop index in mpp case for jpjwfm2,jpjwd
(  218) 
(  219)   REAL(wp), DIMENSION(jpj) ::   &  !:
(  220)      bsfwob              !: now barotropic stream fuction computed at the OBC. The corres-
(  221)      !                   !  ponding bsfn will be computed by the forward time step in dynspg.
(  222) 
(  223)   REAL(wp), DIMENSION(jpj,3,3) ::   &  !:
(  224)      bwbnd               !: West boundary barotropic streamfunction over
(  225)      !                   !  3 rows and 3 time step (now, before, and before before)
(  226) 
(  227) !DB
(  228)   REAL(wp), DIMENSION(1:jpjwf) ::   &  !:
(  229)      bfow,             & !: now climatology of the west boundary barotropic stream function
(  230)      sshfow,           & !: now climatology of the west boundary sea surface height
(  231)      eta_w,       & !: elevation computed from uwdta






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 14

(  232)      ubtfow,vbtfow       !: now climatology of the west boundary barotropic transport
(  233)      
(  234) 
(  235)   REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
(  236)      ufow, vfow,       & !: now climatology of the west velocities
(  237) !     tfow, sfow,       & !: now climatology of the west temperature and salinity
(  238)      ucliw               !: baroclinic componant of the zonal velocity after the radiation
(  239)      !                   !  in the obcdyn.F90 routine
(  240) 
(  241) !byoung
(  242)   REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
(  243)      tfow, sfow        !: now climatology of the west temperature and salinity
(  244)      
(  245)   REAL(wp), DIMENSION(1:jpjwf,jpj) ::   &  !:
(  246)      sshfow_b            !: west boundary ssh correction averaged over the barotropic loop
(  247)                          !: (if Flather's algoritm applied at open boundary)
(  248) 
(  249)   REAL(wp), DIMENSION(1:jpjwf,0:jptobc+1) ::   &  !:
(  250)      sshwdta, ubtwdta    !: array used for interpolating monthly data on the west boundary
(  251) 
(  252) !ylu
(  253) !  REAL(wp), DIMENSION(1:jpjwf,4) ::    &    !:
(  254)   REAL(wp), DIMENSION(1:jpjwf,5) ::    &    !:
(  255)      tidesshwmag, tidesshwpha, tidevbtwmag, tidevbtwpha    
(  256) 
(  257)   REAL(wp), DIMENSION(1:jpjwf,jpk,jptobc) ::   &  !:
(  258)      uwdta !, uwdta1 !: array used for interpolating monthly data on the west boundary
(  259) !     uwdta, twdta, swdta !: array used for interpolating monthly data on the west boundary
(  260) 
(  261) !DBG
(  262)   REAL(wp), DIMENSION(1:jpjwf,jpk) ::   &  !:
(  263)      uwdta1 !: array used for interpolating monthly data on the west boundary
(  264) 
(  265) 
(  266) !byoung
(  267)   REAL(wp), DIMENSION(1:jpjwf,jpk,5,jptobc) ::   &  !:
(  268)      twdta, swdta, twdta1, swdta1 !: array used for interpolating monthly data on the west boundary
(  269) 
(  270)   !!-------------------------------
(  271)   !! Arrays for radiative West OBC
(  272)   !!-------------------------------
(  273)   REAL(wp), DIMENSION(jpj,jpk,3,3) ::   &  !:
(  274)      uwbnd, vwbnd                  !: baroclinic u & v components of the velocity over 3 rows
(  275)      !                             !  and 3 time step (now, before, and before before)
(  276) 
(  277)   REAL(wp), DIMENSION(jpj,jpk,2,2) ::   &  !:
(  278)      twbnd, swbnd                  !: west boundary temperature and salinity over 2 rows and
(  279)      !                             !  2 time step (now and before)
(  280) 
(  281)   REAL(wp), DIMENSION(jpj,jpk) ::    &  !:
(  282)      u_cxwbnd, v_cxwbnd            !: Zonal component of the phase speed ratio computed with
(  283)      !                             !  radiation of zonal and meridional velocity (respectively)
(  284)      !                             !  at the west open boundary (u_cxwbnd = cx rdt )
(  285) 
(  286)   REAL(wp), DIMENSION(jpj,jpk) ::    &  !:
(  287)      uwmsk, vwmsk, twmsk           !: 2D mask for the West OB
(  288) 
(  289) !byoung






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 15

(  290)   REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
(  291)      twmsk5
(  292) 
(  293)   ! Note that those arrays are optimized for mpp case
(  294)   ! (hence the dimension jpj is the size of one processor subdomain)
(  295) 
(  296)   !!---------------------
(  297)   !! North open boundary
(  298)   !!---------------------
(  299)   INTEGER ::   nin0  , nin1       !: do loop index in mpp case for jpind, jpinf
(  300)   INTEGER ::   nin0p1, nin1m1     !: do loop index in mpp case for jpindp1, jpinfm1
(  301)   INTEGER ::   nin1m2, nin0m1     !: do loop index in mpp case for jpinfm1-1,jpind
(  302)   INTEGER ::   njn0  , njn1       !: do loop index in mpp case for jpnob
(  303)   INTEGER ::   njn0p1, njn1p1     !: do loop index in mpp case for jpnob+1
(  304)   INTEGER ::   njn0m1, njn1m1     !: do loop index in mpp case for jpnob-1
(  305) 
(  306)   REAL(wp), DIMENSION(jpi) ::   &  !:
(  307)      bsfnob              !: now barotropic stream fuction computed at the OBC. The corres-
(  308)      !                   !  ponding bsfn will be computed by the forward time step in dynspg.
(  309) 
(  310)   REAL(wp), DIMENSION(jpi,3,3) ::   &  !:
(  311)      bnbnd               !: north boundary barotropic streamfunction over
(  312)      !                   !  3 rows and 3 time step (now, before, and before before)
(  313) 
(  314)   REAL(wp), DIMENSION(1:jpinf) ::   &  !:
(  315)      bfon,             & !: now climatology of the north boundary barotropic stream function
(  316)      sshfon,           & !: now climatology of the north boundary sea surface height
(  317)      ubtfon,vbtfon       !: now climatology of the north boundary barotropic transport
(  318) 
(  319)   REAL(wp), DIMENSION(jpi,jpk) ::   &    !:
(  320)      ufon, vfon,       & !: now climatology of the north boundary velocities
(  321) !      tfon, sfon,       & !: now climatology of the north boundary temperature and salinity
(  322)      vclin               !: baroclinic componant of the meridian velocity after the radiation
(  323)      !                   !  in yhe obcdyn.F90 routine
(  324) 
(  325) !sujie
(  326)   REAL(wp), DIMENSION(jpi,jpk,5) ::   &    !:
(  327)      tfon, sfon        !: now climatology of the north boundary temperature and salinity
(  328) 
(  329)   REAL(wp), DIMENSION(1:jpinf,jpj) ::   &  !:
(  330)      sshfon_b            !: north boundary ssh correction averaged over the barotropic loop
(  331)                          !: (if Flather's algoritm applied at open boundary)
(  332) 
(  333)   REAL(wp), DIMENSION(1:jpinf,0:jptobc+1) ::   &  !:
(  334)      sshndta, vbtndta   !: array used for interpolating monthly data on the north boundary
(  335) 
(  336) !ylu
(  337)   REAL(wp), DIMENSION(1:jpinf,5) ::    &    !:
(  338)      tidesshnmag, tidesshnpha, tidevbtnmag, tidevbtnpha    
(  339) 
(  340)   REAL(wp), DIMENSION(1:jpinf,jpk,jptobc) ::   &  !:
(  341)      vndta !sujie, tndta, sndta !: array used for interpolating monthly data on the north boundary
(  342) !sujie
(  343)   REAL(wp), DIMENSION(1:jpinf,jpk,5,jptobc) ::   &  !:
(  344)      tndta, sndta !: array used for interpolating monthly data on the north boundary
(  345) 
(  346)   !!--------------------------------
(  347)   !! Arrays for radiative North OBC






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 16

(  348)   !!--------------------------------
(  349)   !!
(  350)   REAL(wp), DIMENSION(jpi,jpk,3,3) ::   &   !:
(  351)      unbnd, vnbnd                  !: baroclinic u & v components of the velocity over 3
(  352)      !                             !  rows and 3 time step (now, before, and before before)
(  353) 
(  354)   REAL(wp), DIMENSION(jpi,jpk,2,2) ::   &   !:
(  355)      tnbnd, snbnd                  !: north boundary temperature and salinity over
(  356)      !                             !  2 rows and 2 time step (now and before)
(  357) 
(  358)   REAL(wp), DIMENSION(jpi,jpk) ::   &     !:
(  359)      u_cynbnd, v_cynbnd            !: Meridional component of the phase speed ratio compu-
(  360)      !                             !  ted with radiation of zonal and meridional velocity
(  361)      !                             !  (respectively) at the north OB (u_cynbnd = cx rdt )
(  362) 
(  363)   REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
(  364)      unmsk, vnmsk, tnmsk           !: 2D mask for the North OB
(  365) 
(  366) !sujie
(  367)   REAL(wp), DIMENSION(jpi,jpk,5) ::   &  !:
(  368)      tnmsk5
(  369) 
(  370)   ! Note that those arrays are optimized for mpp case
(  371)   ! (hence the dimension jpj is the size of one processor subdomain)
(  372) 
(  373)   !!---------------------
(  374)   !! South open boundary
(  375)   !!---------------------
(  376)   INTEGER ::   nis0  , nis1       !: do loop index in mpp case for jpisd, jpisf
(  377)   INTEGER ::   nis0p1, nis1m1     !: do loop index in mpp case for jpisdp1, jpisfm1
(  378)   INTEGER ::   nis1m2, nis0m1     !: do loop index in mpp case for jpisfm1-1,jpisd
(  379)   INTEGER ::   njs0  , njs1       !: do loop index in mpp case for jpsob
(  380)   INTEGER ::   njs0p1, njs1p1     !: do loop index in mpp case for jpsob+1
(  381) 
(  382)   REAL(wp), DIMENSION(jpi) ::    &   !:
(  383)      bsfsob              !: now barotropic stream fuction computed at the OBC.The corres-
(  384)      !                   !  ponding bsfn will be computed by the forward time step in dynspg.
(  385)   REAL(wp), DIMENSION(jpi,3,3) ::   &   !:
(  386)      bsbnd               !: south boundary barotropic stream function over
(  387)      !                   !  3 rows and 3 time step (now, before, and before before)
(  388) 
(  389) !DB
(  390)   REAL(wp), DIMENSION(1:jpisf) ::    &   !:
(  391)      bfos,             & !: now climatology of the south boundary barotropic stream function
(  392)      sshfos,           & !: now climatology of the south boundary sea surface height
(  393)      eta_s,            & !: AD3Aug07:  elev. geostropic elev computed from vsdta
(  394)      ubtfos,vbtfos       !: now climatology of the south boundary barotropic transport
(  395) 
(  396)   REAL(wp), DIMENSION(jpi,jpk) ::    &   !:
(  397)      ufos, vfos ,       & !: now climatology of the south boundary velocities
(  398) !      tfos, sfos,       & !: now climatology of the south boundary temperature and salinity
(  399)      vclis               !: baroclinic componant of the meridian velocity after the radiation
(  400)      !                   !  in the obcdyn.F90 routine
(  401) 
(  402) !sujie
(  403)   REAL(wp), DIMENSION(jpi,jpk,5) ::    &   !:
(  404)      tfos, sfos        !: now climatology of the south boundary temperature and salinity
(  405) 






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 17

(  406) !DB
(  407)   REAL(wp), DIMENSION(1:jpisf,jpj) ::   &  !:
(  408)      sshfos_b            !: south boundary ssh correction averaged over the barotropic loop
(  409)                          !: (if Flather's algoritm applied at open boundary)
(  410) 
(  411)   REAL(wp), DIMENSION(1:jpisf,0:jptobc+1) ::    &    !:
(  412)      sshsdta, vbtsdta    !: array used for interpolating monthly data on the south boundary
(  413) 
(  414) !ylu
(  415)   REAL(wp), DIMENSION(1:jpisf,5) ::    &    !:
(  416)      tidesshsmag, tidesshspha, tidevbtsmag, tidevbtspha    
(  417) 
(  418)   REAL(wp), DIMENSION(1:jpisf,jpk,jptobc) ::    &    !:
(  419)      vsdta !sujie, tsdta, ssdta   !: array used for interpolating monthly data on the south boundary
(  420) 
(  421) !DB
(  422)   REAL(wp), DIMENSION(1:jpisf,jpk) ::    &    !:
(  423)      vsdta1  !AD: array used for interpolatated global obc array
(  424)      
(  425) !sujie
(  426)   REAL(wp), DIMENSION(1:jpisf,jpk,5,jptobc) ::    &    !:
(  427)      tsdta, ssdta, tsdta1, ssdta1   !: array used for interpolating monthly data on the south boundary
(  428) 
(  429)   !!--------------------------------
(  430)   !! Arrays for radiative South OBC
(  431)   !!--------------------------------
(  432)   !!                        computed by the forward time step in dynspg.
(  433)   REAL(wp), DIMENSION(jpi,jpk,3,3) ::   &   !:
(  434)      usbnd, vsbnd                  !: baroclinic u & v components of the velocity over 3
(  435)      !                             !  rows and 3 time step (now, before, and before before)
(  436) 
(  437)   REAL(wp), DIMENSION(jpi,jpk,2,2) ::   &  !:
(  438)      tsbnd, ssbnd                  !: south boundary temperature and salinity over
(  439)      !                             !  2 rows and 2 time step (now and before)
(  440) 
(  441)   REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
(  442)      u_cysbnd, v_cysbnd            !: Meridional component of the phase speed ratio compu-
(  443)      !                             !  ted with radiation of zonal and meridional velocity
(  444)      !                             !  (repsectively) at the south OB (u_cynbnd = cx rdt )
(  445) 
(  446)   REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
(  447)      usmsk, vsmsk, tsmsk           !: 2D mask for the South OB
(  448) 
(  449) !sujie
(  450)   REAL(wp), DIMENSION(jpi,jpk,5) ::   &  !:
(  451)      tsmsk5
(  452) 
(  453)   ! Note that those arrays are optimized for mpp case
(  454)   ! (hence the dimension jpj is the size of one processor subdomain)
(  455) 
(  456) 
(  457) # 462
(  462)   !!======================================================================
(  463) END MODULE obc_oce
