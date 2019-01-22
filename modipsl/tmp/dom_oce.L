


PGF90 (Version     10.5)          09/23/2010  18:36:35      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /distrib/local/netcdf/pgf/include/
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: dom_oce.F90

(    1) MODULE dom_oce
(    2)    !!----------------------------------------------------------------------
(    3)    !!                       ***  MODULE dom_oce  ***
(    4)    !!       
(    5)    !! ** Purpose :   Define in memory all the ocean space domain variables
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/dom_oce.F90,v 1.10 2006/03/10 10:55:38 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE par_oce      ! ocean parameters
(   13) 
(   14)    IMPLICIT NONE
(   15)    PUBLIC           ! allows the acces to par_oce when dom_oce is used
(   16)    !                ! exception to coding rules... to be suppressed ???
(   17) 
(   18)    !!----------------------------------------------------------------------
(   19)    !! space domain parameters
(   20)    !! -----------------------
(   21)    LOGICAL, PUBLIC ::   &   !:
(   22)       lzoom      =  .FALSE. ,   &  !: zoom flag
(   23)       lzoom_e    =  .FALSE. ,   &  !: East  zoom type flag
(   24)       lzoom_w    =  .FALSE. ,   &  !: West  zoom type flag
(   25)       lzoom_s    =  .FALSE. ,   &  !: South zoom type flag
(   26)       lzoom_n    =  .FALSE. ,   &  !: North zoom type flag
(   27)       lzoom_arct =  .FALSE. ,   &  !: ORCA    arctic zoom flag
(   28)       lzoom_anta =  .FALSE.        !: ORCA antarctic zoom flag
(   29) 
(   30)    INTEGER, PUBLIC ::           & !!: namdom : space domain (bathymetry, mesh)
(   31)       ntopo   =  0 ,            &  !: = 0/1 ,compute/read the bathymetry file
(   32)       ngrid   =  0 ,            &  !: = 0/1, compute/read the horizontal mesh file
(   33)       nmsh    =  0                 !: = 1 create a mesh-mask file
(   34) 
(   35)    INTEGER, PUBLIC ::         &   !:
(   36)       ! domain parameters linked to mpp
(   37)       nperio,          &  !: type of lateral boundary condition
(   38)       nimpp, njmpp,    &  !: i- & j-indexes for mpp-subdomain left bottom
(   39)       nreci, nrecj,    &  !: overlap region in i and j
(   40)       nproc,           &  !: number for local processor
(   41)       narea,           &  !: number for local area
(   42)       nbondi, nbondj,  &  !: mark of i- and j-direction local boundaries
(   43)       npolj,           &  !: north fold mark (0, 3 or 4)
(   44)       nlci, nlcj,      &  !: i- & j-dimensions of the local subdomain
(   45)       nldi, nlei,      &  !: first and last indoor i- and j-indexes
(   46)       nldj, nlej,      &  !:
(   47)       noea, nowe,      &  !: index of the local neighboring processors in
(   48)       noso, nono,      &  !: east, west, south and north directions
(   49)       npne, npnw,      &  !: index of north east and north west processor






PGF90 (Version     10.5)          09/23/2010  18:36:35      page 2

(   50)       npse, npsw,      &  !: index of south east and south west processor
(   51)       nbne, nbnw,      &  !: logical of north east & north west processor
(   52)       nbse, nbsw,      &  !: logical of south east & south west processor
(   53)       nidom
(   54) 
(   55)    INTEGER, PUBLIC, DIMENSION(jpi) ::   &   !:
(   56)       mig                 !: local  ==> global  domain i-indice
(   57)    INTEGER, PUBLIC, DIMENSION(jpj) ::   &   !:
(   58)       mjg                 !: local  ==> global  domain j-indice
(   59)    INTEGER, PUBLIC, DIMENSION( jpidta ) ::   &  !:  !!bug ==> other solution?
(   60)       mi0, mi1            !: global ==> local domain i-indice
(   61)       !                   !  (mi0=1 and mi1=0 if the global indice is not in the local domain)
(   62)    INTEGER, PUBLIC, DIMENSION( jpjdta ) ::   &  !:
(   63)       mj0, mj1            !: global ==> local domain j-indice
(   64)       !                   ! (mi0=1 and mi1=0 if the global indice is not in the local domain)
(   65) 
(   66)    INTEGER, PUBLIC, DIMENSION(jpnij) ::   &  !:
(   67)       nimppt, njmppt,  &  !: i-, j-indexes for each processor
(   68)       ibonit, ibonjt,  &  !: i-, j- processor neighbour existence
(   69)       nlcit, nlcjt,    &  !: dimensions of every subdomain
(   70)       nldit, nldjt,    &  !: first, last indoor index for each i-domain
(   71)       nleit, nlejt        !: first, last indoor index for each j-domain
(   72) 
(   73)    !!----------------------------------------------------------------------
(   74)    !! horizontal curvilinear coordinate and scale factors
(   75)    !! ---------------------------------------------------------------------
(   76) 
(   77)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   78)       glamt, glamu,    &  !: longitude of t-, u-, v- and f-points (degre)
(   79)       glamv, glamf,    &  !:
(   80)       gphit, gphiu,    &  !: latitude  of t-, u-, v- and f-points (degre)
(   81)       gphiv, gphif,    &  !:
(   82)       e1t, e2t,        &  !: horizontal scale factors at t-point (m)
(   83)       e1u, e2u,        &  !: horizontal scale factors at u-point (m)
(   84)       e1v, e2v,        &  !: horizontal scale factors at v-point (m)
(   85)       e1f, e2f,        &  !: horizontal scale factors at f-point (m)
(   86)       ff                  !: coriolis factor (2.*omega*sin(yphi) ) (s-1)
(   87) 
(   88) !!DB 2007.12.11 -- special arrays needed in obcdta
(   89) !!The coriolis arrays added 2008.04.17
(   90) !!They are initialized in domhgr <--- hgr_read()
(   91) !!NB: the *_n (north) arrays are not used as this boundary is closed
(   92)       REAL(wp), dimension(jpidta) :: e2u_s, e2u_n, ff_s, ff_n
(   93)       REAL(wp), dimension(jpjdta) :: e2v_e, e2v_w, ff_e, ff_w
(   94) 
(   95) 
(   96)    !!----------------------------------------------------------------------
(   97)    !! vertical coordinate and scale factors
(   98)    !! --------------------------------------
(   99) 
(  100)    REAL(wp), PUBLIC ::   & !!: * namelist namdom *
(  101)       e3zps_min = 5.0,   &  !: miminum thickness for partial steps (meters)
(  102)       e3zps_rat = 0.1       !: minimum thickness ration for partial steps
(  103) 
(  104)    !! z-coordinate (default option) (also used in the other cases
(  105)    !! -----------------------------  as reference z-coordinate)
(  106)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
(  107)       gdept, gdepw,    &  !: reference depth of t- and w-points (m)






PGF90 (Version     10.5)          09/23/2010  18:36:35      page 3

(  108)       e3t, e3w            !: reference vertical scale factors at T- and W-pts (m)
(  109) 
(  110) #if defined key_partial_steps
(  111)    !! Partial steps ('key_partial_steps')
(  112)    !! -----------------------------------
(  113)    LOGICAL, PUBLIC, PARAMETER ::   lk_zps = .TRUE.   !: partial steps flag
(  114)    LOGICAL, PUBLIC, PARAMETER ::   lk_sco = .FALSE.  !: s-coordinate flag
(  115)    LOGICAL, PUBLIC, PARAMETER ::   lk_zco = .FALSE.  !: z-coordinate flag
(  116)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(  117)       gdep3w,                 &  !: ???
(  118)       gdept_ps, gdepw_ps,     &  !: depth of t- and w-points (m)
(  119)       e3t_ps, e3u_ps, e3v_ps, &  !: vertical scale factors at t-, u-, w-,
(  120)       e3w_ps, e3f_ps,         &  !: w- and f- points (m)
(  121)       e3uw_ps, e3vw_ps           !: uw- and vw- points (m)
(  122) 
(  123)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(  124)       hdept, hdepw, e3tp, e3wp   !: ???
(  125) 
(  126) #elif defined key_s_coord
(  127)    !! s-coordinate ('key_s_coord')
(  128)    !! ----------------------------
(  129)    LOGICAL, PUBLIC, PARAMETER ::   lk_zps = .FALSE.   !: partial steps flag
(  130)    LOGICAL, PUBLIC, PARAMETER ::   lk_sco = .TRUE.    !: s-coordinate flag
(  131)    LOGICAL, PUBLIC, PARAMETER ::   lk_zco = .FALSE.   !: z-coordinate flag
(  132)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(  133)       hbatt, hbatu,    &  !: ocean depth at the vertical of  t-, u-, v-
(  134)       hbatv, hbatf        !: and f-point (m)
(  135) 
(  136)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &   !:
(  137)       gsigt, gsigw ,   &  !: model level depth coefficient at t-, w-levels
(  138)       gsi3w,           &  !: model level depth coefficient at w-level
(  139)                           !  defined as the sum of e3w scale factors
(  140)       esigt, esigw        !: vertical scale factor coef. at t-, w-levels
(  141) 
(  142) #else
(  143)    !! z-coordinate (Default option)
(  144)    !! -----------------------------
(  145)    LOGICAL, PUBLIC, PARAMETER ::   lk_zps = .FALSE.   !: partial steps flag
(  146)    LOGICAL, PUBLIC, PARAMETER ::   lk_sco = .FALSE.   !: s-coordinate flag
(  147)    LOGICAL, PUBLIC, PARAMETER ::   lk_zco = .TRUE.    !: s-coordinate flag
(  148) #endif
(  149)    !!----------------------------------------------------------------------
(  150)    !! masks, bathymetry
(  151)    !! -----------------
(  152) 
(  153)    INTEGER , PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(  154)       mbathy     !: number of ocean level (=0, 1, ... , jpk-1)
(  155) 
(  156)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(  157)       bathy  ,         &  !: ocean depth (meters)
(  158)       tmask_i             !: interior domain T-point mask
(  159) 
(  160)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(  161)       tmask, umask,    &  !: land/ocean mask at T-, U-, V- and F-points
(  162)       vmask, fmask,nmask        !:
(  163) 
(  164)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(  165)       bmask               !: land/ocean mask of barotropic stream function






PGF90 (Version     10.5)          09/23/2010  18:36:35      page 4

(  166) 
(  167)    REAL(wp), PUBLIC, DIMENSION(jpiglo) ::   &   !:
(  168)       tpol, fpol          !: north fold mask (nperio= 3 or 4)
(  169) 
(  170) #if defined key_noslip_accurate
(  171)    INTEGER, PUBLIC, DIMENSION(4,jpk) ::   &   !:
(  172)       npcoa               !: ???
(  173)    INTEGER, PUBLIC, DIMENSION(2*(jpi+jpj),4,jpk) ::   &   !:
(  174)       nicoa,           &  !: ???
(  175)       njcoa               !: ???
(  176) 
(  177) #endif
(  178) 
(  179)    !!----------------------------------------------------------------------
(  180)    !! time domain
(  181)    !!----------------------------------------------------------------------
(  182)    INTEGER, PUBLIC ::    & !!: * Namelist * ???
(  183)       nacc   = 0 ,       &  !: = 0/1 use of the acceleration of convergence technique
(  184)       neuler                !: restart euler forward option (0=Euler)
(  185) 
(  186) 
(  187)    REAL(wp), PUBLIC ::       & !!: * Namelist ??? *
(  188)       rdt    = 3600._wp ,    &  !: time step for the dynamics (and tracer if nacc=0)
(  189)       rdtmin = 3600._wp ,    &  !: minimum time step on tracers
(  190)       rdtmax = 3600._wp ,    &  !: maximum time step on tracers
(  191)       rdth   =  800._wp ,    &  !: depth variation of tracer step
(  192)       rdtbt  =   90._wp ,    &  !: barotropic time step for the dynamics (lk_dynspg_ts=T)
(  193)       atfp   = 0.1_wp   ,    &  !: asselin time filter parameter
(  194)       atfp1                     !: asselin time filter coeff. (atfp1= 1-2*atfp)
(  195) 
(  196)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
(  197)       rdttra                    !: vertical profile of tracer time step
(  198) 
(  199)    !!----------------------------------------------------------------------
(  200)    !! cross land advection
(  201)    !!----------------------------------------------------------------------
(  202) 
(  203)    INTEGER, PUBLIC ::       & !!: namelist ???
(  204)       n_cla                    !: flag (0/1) for cross land advection to
(  205)       !                        ! parameterize exchanges through straits
(  206) 
(  207) #if defined key_agrif
(  208)    !!----------------------------------------------------------------------
(  209)    !! agrif sponge layer
(  210)    !!----------------------------------------------------------------------
(  211)       LOGICAL :: spongedoneT = .FALSE.
(  212)       REAL(wp), DIMENSION(jpi,jpj) :: zspe1ur, zspe2vr ,zspbtr2
(  213)    !!----------------------------------------------------------------------
(  214) #endif
(  215) 
(  216) END MODULE dom_oce













PGF90 (Version     10.5)          09/23/2010  18:36:35      page 5

(    1) # 1 "dom_oce.F90"
(    1) MODULE dom_oce
(    2)    !!----------------------------------------------------------------------
(    3)    !!                       ***  MODULE dom_oce  ***
(    4)    !!       
(    5)    !! ** Purpose :   Define in memory all the ocean space domain variables
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/dom_oce.F90,v 1.10 2006/03/10 10:55:38 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE par_oce      ! ocean parameters
(   13) 
(   14)    IMPLICIT NONE
(   15)    PUBLIC           ! allows the acces to par_oce when dom_oce is used
(   16)    !                ! exception to coding rules... to be suppressed ???
(   17) 
(   18)    !!----------------------------------------------------------------------
(   19)    !! space domain parameters
(   20)    !! -----------------------
(   21)    LOGICAL, PUBLIC ::   &   !:
(   22)       lzoom      =  .FALSE. ,   &  !: zoom flag
(   23)       lzoom_e    =  .FALSE. ,   &  !: East  zoom type flag
(   24)       lzoom_w    =  .FALSE. ,   &  !: West  zoom type flag
(   25)       lzoom_s    =  .FALSE. ,   &  !: South zoom type flag
(   26)       lzoom_n    =  .FALSE. ,   &  !: North zoom type flag
(   27)       lzoom_arct =  .FALSE. ,   &  !: ORCA    arctic zoom flag
(   28)       lzoom_anta =  .FALSE.        !: ORCA antarctic zoom flag
(   29) 
(   30)    INTEGER, PUBLIC ::           & !!: namdom : space domain (bathymetry, mesh)
(   31)       ntopo   =  0 ,            &  !: = 0/1 ,compute/read the bathymetry file
(   32)       ngrid   =  0 ,            &  !: = 0/1, compute/read the horizontal mesh file
(   33)       nmsh    =  0                 !: = 1 create a mesh-mask file
(   34) 
(   35)    INTEGER, PUBLIC ::         &   !:
(   36)       ! domain parameters linked to mpp
(   37)       nperio,          &  !: type of lateral boundary condition
(   38)       nimpp, njmpp,    &  !: i- & j-indexes for mpp-subdomain left bottom
(   39)       nreci, nrecj,    &  !: overlap region in i and j
(   40)       nproc,           &  !: number for local processor
(   41)       narea,           &  !: number for local area
(   42)       nbondi, nbondj,  &  !: mark of i- and j-direction local boundaries
(   43)       npolj,           &  !: north fold mark (0, 3 or 4)
(   44)       nlci, nlcj,      &  !: i- & j-dimensions of the local subdomain
(   45)       nldi, nlei,      &  !: first and last indoor i- and j-indexes
(   46)       nldj, nlej,      &  !:
(   47)       noea, nowe,      &  !: index of the local neighboring processors in
(   48)       noso, nono,      &  !: east, west, south and north directions
(   49)       npne, npnw,      &  !: index of north east and north west processor
(   50)       npse, npsw,      &  !: index of south east and south west processor
(   51)       nbne, nbnw,      &  !: logical of north east & north west processor
(   52)       nbse, nbsw,      &  !: logical of south east & south west processor
(   53)       nidom
(   54) 
(   55)    INTEGER, PUBLIC, DIMENSION(jpi) ::   &   !:
(   56)       mig                 !: local  ==> global  domain i-indice
(   57)    INTEGER, PUBLIC, DIMENSION(jpj) ::   &   !:






PGF90 (Version     10.5)          09/23/2010  18:36:35      page 6

(   58)       mjg                 !: local  ==> global  domain j-indice
(   59)    INTEGER, PUBLIC, DIMENSION( jpidta ) ::   &  !:  !!bug ==> other solution?
(   60)       mi0, mi1            !: global ==> local domain i-indice
(   61)       !                   !  (mi0=1 and mi1=0 if the global indice is not in the local domain)
(   62)    INTEGER, PUBLIC, DIMENSION( jpjdta ) ::   &  !:
(   63)       mj0, mj1            !: global ==> local domain j-indice
(   64)       !                   ! (mi0=1 and mi1=0 if the global indice is not in the local domain)
(   65) 
(   66)    INTEGER, PUBLIC, DIMENSION(jpnij) ::   &  !:
(   67)       nimppt, njmppt,  &  !: i-, j-indexes for each processor
(   68)       ibonit, ibonjt,  &  !: i-, j- processor neighbour existence
(   69)       nlcit, nlcjt,    &  !: dimensions of every subdomain
(   70)       nldit, nldjt,    &  !: first, last indoor index for each i-domain
(   71)       nleit, nlejt        !: first, last indoor index for each j-domain
(   72) 
(   73)    !!----------------------------------------------------------------------
(   74)    !! horizontal curvilinear coordinate and scale factors
(   75)    !! ---------------------------------------------------------------------
(   76) 
(   77)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   78)       glamt, glamu,    &  !: longitude of t-, u-, v- and f-points (degre)
(   79)       glamv, glamf,    &  !:
(   80)       gphit, gphiu,    &  !: latitude  of t-, u-, v- and f-points (degre)
(   81)       gphiv, gphif,    &  !:
(   82)       e1t, e2t,        &  !: horizontal scale factors at t-point (m)
(   83)       e1u, e2u,        &  !: horizontal scale factors at u-point (m)
(   84)       e1v, e2v,        &  !: horizontal scale factors at v-point (m)
(   85)       e1f, e2f,        &  !: horizontal scale factors at f-point (m)
(   86)       ff                  !: coriolis factor (2.*omega*sin(yphi) ) (s-1)
(   87) 
(   88) !!DB 2007.12.11 -- special arrays needed in obcdta
(   89) !!The coriolis arrays added 2008.04.17
(   90) !!They are initialized in domhgr <--- hgr_read()
(   91) !!NB: the *_n (north) arrays are not used as this boundary is closed
(   92)       REAL(wp), dimension(jpidta) :: e2u_s, e2u_n, ff_s, ff_n
(   93)       REAL(wp), dimension(jpjdta) :: e2v_e, e2v_w, ff_e, ff_w
(   94) 
(   95) 
(   96)    !!----------------------------------------------------------------------
(   97)    !! vertical coordinate and scale factors
(   98)    !! --------------------------------------
(   99) 
(  100)    REAL(wp), PUBLIC ::   & !!: * namelist namdom *
(  101)       e3zps_min = 5.0,   &  !: miminum thickness for partial steps (meters)
(  102)       e3zps_rat = 0.1       !: minimum thickness ration for partial steps
(  103) 
(  104)    !! z-coordinate (default option) (also used in the other cases
(  105)    !! -----------------------------  as reference z-coordinate)
(  106)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
(  107)       gdept, gdepw,    &  !: reference depth of t- and w-points (m)
(  108)       e3t, e3w            !: reference vertical scale factors at T- and W-pts (m)
(  109) 
(  110) # 143
(  143)    !! z-coordinate (Default option)
(  144)    !! -----------------------------
(  145)    LOGICAL, PUBLIC, PARAMETER ::   lk_zps = .FALSE.   !: partial steps flag
(  146)    LOGICAL, PUBLIC, PARAMETER ::   lk_sco = .FALSE.   !: s-coordinate flag
(  147)    LOGICAL, PUBLIC, PARAMETER ::   lk_zco = .TRUE.    !: s-coordinate flag






PGF90 (Version     10.5)          09/23/2010  18:36:35      page 7

(  148) # 149
(  149)    !!----------------------------------------------------------------------
(  150)    !! masks, bathymetry
(  151)    !! -----------------
(  152) 
(  153)    INTEGER , PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(  154)       mbathy     !: number of ocean level (=0, 1, ... , jpk-1)
(  155) 
(  156)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(  157)       bathy  ,         &  !: ocean depth (meters)
(  158)       tmask_i             !: interior domain T-point mask
(  159) 
(  160)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(  161)       tmask, umask,    &  !: land/ocean mask at T-, U-, V- and F-points
(  162)       vmask, fmask,nmask        !:
(  163) 
(  164)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(  165)       bmask               !: land/ocean mask of barotropic stream function
(  166) 
(  167)    REAL(wp), PUBLIC, DIMENSION(jpiglo) ::   &   !:
(  168)       tpol, fpol          !: north fold mask (nperio= 3 or 4)
(  169) 
(  170) 
(  171) # 179
(  179)    !!----------------------------------------------------------------------
(  180)    !! time domain
(  181)    !!----------------------------------------------------------------------
(  182)    INTEGER, PUBLIC ::    & !!: * Namelist * ???
(  183)       nacc   = 0 ,       &  !: = 0/1 use of the acceleration of convergence technique
(  184)       neuler                !: restart euler forward option (0=Euler)
(  185) 
(  186) 
(  187)    REAL(wp), PUBLIC ::       & !!: * Namelist ??? *
(  188)       rdt    = 3600._wp ,    &  !: time step for the dynamics (and tracer if nacc=0)
(  189)       rdtmin = 3600._wp ,    &  !: minimum time step on tracers
(  190)       rdtmax = 3600._wp ,    &  !: maximum time step on tracers
(  191)       rdth   =  800._wp ,    &  !: depth variation of tracer step
(  192)       rdtbt  =   90._wp ,    &  !: barotropic time step for the dynamics (lk_dynspg_ts=T)
(  193)       atfp   = 0.1_wp   ,    &  !: asselin time filter parameter
(  194)       atfp1                     !: asselin time filter coeff. (atfp1= 1-2*atfp)
(  195) 
(  196)    REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
(  197)       rdttra                    !: vertical profile of tracer time step
(  198) 
(  199)    !!----------------------------------------------------------------------
(  200)    !! cross land advection
(  201)    !!----------------------------------------------------------------------
(  202) 
(  203)    INTEGER, PUBLIC ::       & !!: namelist ???
(  204)       n_cla                    !: flag (0/1) for cross land advection to
(  205)       !                        ! parameterize exchanges through straits
(  206) 
(  207) 
(  208) # 216
(  216) END MODULE dom_oce
