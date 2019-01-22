MODULE dom_oce
   !!----------------------------------------------------------------------
   !!                       ***  MODULE dom_oce  ***
   !!       
   !! ** Purpose :   Define in memory all the ocean space domain variables
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/dom_oce.F90,v 1.10 2006/03/10 10:55:38 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce      ! ocean parameters

   IMPLICIT NONE
   PUBLIC           ! allows the acces to par_oce when dom_oce is used
   !                ! exception to coding rules... to be suppressed ???

   !!----------------------------------------------------------------------
   !! space domain parameters
   !! -----------------------
   LOGICAL, PUBLIC ::   &   !:
      lzoom      =  .FALSE. ,   &  !: zoom flag
      lzoom_e    =  .FALSE. ,   &  !: East  zoom type flag
      lzoom_w    =  .FALSE. ,   &  !: West  zoom type flag
      lzoom_s    =  .FALSE. ,   &  !: South zoom type flag
      lzoom_n    =  .FALSE. ,   &  !: North zoom type flag
      lzoom_arct =  .FALSE. ,   &  !: ORCA    arctic zoom flag
      lzoom_anta =  .FALSE.        !: ORCA antarctic zoom flag

   INTEGER, PUBLIC ::           & !!: namdom : space domain (bathymetry, mesh)
      ntopo   =  0 ,            &  !: = 0/1 ,compute/read the bathymetry file
      ngrid   =  0 ,            &  !: = 0/1, compute/read the horizontal mesh file
      nmsh    =  0                 !: = 1 create a mesh-mask file

   INTEGER, PUBLIC ::         &   !:
      ! domain parameters linked to mpp
      nperio,          &  !: type of lateral boundary condition
      nimpp, njmpp,    &  !: i- & j-indexes for mpp-subdomain left bottom
      nreci, nrecj,    &  !: overlap region in i and j
      nproc,           &  !: number for local processor
      narea,           &  !: number for local area
      nbondi, nbondj,  &  !: mark of i- and j-direction local boundaries
      npolj,           &  !: north fold mark (0, 3 or 4)
      nlci, nlcj,      &  !: i- & j-dimensions of the local subdomain
      nldi, nlei,      &  !: first and last indoor i- and j-indexes
      nldj, nlej,      &  !:
      noea, nowe,      &  !: index of the local neighboring processors in
      noso, nono,      &  !: east, west, south and north directions
      npne, npnw,      &  !: index of north east and north west processor
      npse, npsw,      &  !: index of south east and south west processor
      nbne, nbnw,      &  !: logical of north east & north west processor
      nbse, nbsw,      &  !: logical of south east & south west processor
      nidom

   INTEGER, PUBLIC, DIMENSION(jpi) ::   &   !:
      mig                 !: local  ==> global  domain i-indice
   INTEGER, PUBLIC, DIMENSION(jpj) ::   &   !:
      mjg                 !: local  ==> global  domain j-indice
   INTEGER, PUBLIC, DIMENSION( jpidta ) ::   &  !:  !!bug ==> other solution?
      mi0, mi1            !: global ==> local domain i-indice
      !                   !  (mi0=1 and mi1=0 if the global indice is not in the local domain)
   INTEGER, PUBLIC, DIMENSION( jpjdta ) ::   &  !:
      mj0, mj1            !: global ==> local domain j-indice
      !                   ! (mi0=1 and mi1=0 if the global indice is not in the local domain)

   INTEGER, PUBLIC, DIMENSION(jpnij) ::   &  !:
      nimppt, njmppt,  &  !: i-, j-indexes for each processor
      ibonit, ibonjt,  &  !: i-, j- processor neighbour existence
      nlcit, nlcjt,    &  !: dimensions of every subdomain
      nldit, nldjt,    &  !: first, last indoor index for each i-domain
      nleit, nlejt        !: first, last indoor index for each j-domain

   !!----------------------------------------------------------------------
   !! horizontal curvilinear coordinate and scale factors
   !! ---------------------------------------------------------------------

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      glamt, glamu,    &  !: longitude of t-, u-, v- and f-points (degre)
      glamv, glamf,    &  !:
      gphit, gphiu,    &  !: latitude  of t-, u-, v- and f-points (degre)
      gphiv, gphif,    &  !:
      e1t, e2t,        &  !: horizontal scale factors at t-point (m)
      e1u, e2u,        &  !: horizontal scale factors at u-point (m)
      e1v, e2v,        &  !: horizontal scale factors at v-point (m)
      e1f, e2f,        &  !: horizontal scale factors at f-point (m)
      ff                  !: coriolis factor (2.*omega*sin(yphi) ) (s-1)

!!DB 2007.12.11 -- special arrays needed in obcdta
!!The coriolis arrays added 2008.04.17
!!They are initialized in domhgr <--- hgr_read()
!!NB: the *_n (north) arrays are not used as this boundary is closed
      REAL(wp), dimension(jpidta) :: e2u_s, e2u_n, ff_s, ff_n
      REAL(wp), dimension(jpjdta) :: e2v_e, e2v_w, ff_e, ff_w


   !!----------------------------------------------------------------------
   !! vertical coordinate and scale factors
   !! --------------------------------------

   REAL(wp), PUBLIC ::   & !!: * namelist namdom *
      e3zps_min = 5.0,   &  !: miminum thickness for partial steps (meters)
      e3zps_rat = 0.1       !: minimum thickness ration for partial steps

   !! z-coordinate (default option) (also used in the other cases
   !! -----------------------------  as reference z-coordinate)
   REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
      gdept, gdepw,    &  !: reference depth of t- and w-points (m)
      e3t, e3w            !: reference vertical scale factors at T- and W-pts (m)

#if defined key_partial_steps
   !! Partial steps ('key_partial_steps')
   !! -----------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zps = .TRUE.   !: partial steps flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_sco = .FALSE.  !: s-coordinate flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_zco = .FALSE.  !: z-coordinate flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
      gdep3w,                 &  !: ???
      gdept_ps, gdepw_ps,     &  !: depth of t- and w-points (m)
      e3t_ps, e3u_ps, e3v_ps, &  !: vertical scale factors at t-, u-, w-,
      e3w_ps, e3f_ps,         &  !: w- and f- points (m)
      e3uw_ps, e3vw_ps           !: uw- and vw- points (m)

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      hdept, hdepw, e3tp, e3wp   !: ???

#elif defined key_s_coord
   !! s-coordinate ('key_s_coord')
   !! ----------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zps = .FALSE.   !: partial steps flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_sco = .TRUE.    !: s-coordinate flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_zco = .FALSE.   !: z-coordinate flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      hbatt, hbatu,    &  !: ocean depth at the vertical of  t-, u-, v-
      hbatv, hbatf        !: and f-point (m)

   REAL(wp), PUBLIC, DIMENSION(jpk) ::   &   !:
      gsigt, gsigw ,   &  !: model level depth coefficient at t-, w-levels
      gsi3w,           &  !: model level depth coefficient at w-level
                          !  defined as the sum of e3w scale factors
      esigt, esigw        !: vertical scale factor coef. at t-, w-levels

#else
   !! z-coordinate (Default option)
   !! -----------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zps = .FALSE.   !: partial steps flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_sco = .FALSE.   !: s-coordinate flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_zco = .TRUE.    !: s-coordinate flag
#endif
   !!----------------------------------------------------------------------
   !! masks, bathymetry
   !! -----------------

   INTEGER , PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      mbathy     !: number of ocean level (=0, 1, ... , jpk-1)

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      bathy  ,         &  !: ocean depth (meters)
      tmask_i             !: interior domain T-point mask

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      tmask, umask,    &  !: land/ocean mask at T-, U-, V- and F-points
      vmask, fmask,nmask        !:

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      bmask               !: land/ocean mask of barotropic stream function

   REAL(wp), PUBLIC, DIMENSION(jpiglo) ::   &   !:
      tpol, fpol          !: north fold mask (nperio= 3 or 4)

#if defined key_noslip_accurate
   INTEGER, PUBLIC, DIMENSION(4,jpk) ::   &   !:
      npcoa               !: ???
   INTEGER, PUBLIC, DIMENSION(2*(jpi+jpj),4,jpk) ::   &   !:
      nicoa,           &  !: ???
      njcoa               !: ???

#endif

   !!----------------------------------------------------------------------
   !! time domain
   !!----------------------------------------------------------------------
   INTEGER, PUBLIC ::    & !!: * Namelist * ???
      nacc   = 0 ,       &  !: = 0/1 use of the acceleration of convergence technique
      neuler                !: restart euler forward option (0=Euler)


   REAL(wp), PUBLIC ::       & !!: * Namelist ??? *
      rdt    = 3600._wp ,    &  !: time step for the dynamics (and tracer if nacc=0)
      rdtmin = 3600._wp ,    &  !: minimum time step on tracers
      rdtmax = 3600._wp ,    &  !: maximum time step on tracers
      rdth   =  800._wp ,    &  !: depth variation of tracer step
      rdtbt  =   90._wp ,    &  !: barotropic time step for the dynamics (lk_dynspg_ts=T)
      atfp   = 0.1_wp   ,    &  !: asselin time filter parameter
      atfp1                     !: asselin time filter coeff. (atfp1= 1-2*atfp)

   REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
      rdttra                    !: vertical profile of tracer time step

   !!----------------------------------------------------------------------
   !! cross land advection
   !!----------------------------------------------------------------------

   INTEGER, PUBLIC ::       & !!: namelist ???
      n_cla                    !: flag (0/1) for cross land advection to
      !                        ! parameterize exchanges through straits

#if defined key_agrif
   !!----------------------------------------------------------------------
   !! agrif sponge layer
   !!----------------------------------------------------------------------
      LOGICAL :: spongedoneT = .FALSE.
      REAL(wp), DIMENSION(jpi,jpj) :: zspe1ur, zspe2vr ,zspbtr2
   !!----------------------------------------------------------------------
#endif

END MODULE dom_oce
