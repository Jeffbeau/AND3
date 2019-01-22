!!DB 2008.04.14
!!see obc_oce.Q90 for last version of this program

MODULE obc_oce
  !!==============================================================================
  !!                       ***  MODULE obc_oce   ***
  !! Open Boundary Cond. :   define related variables
  !!==============================================================================
#if defined key_obc
  !!----------------------------------------------------------------------
  !!   'key_obc' :                                 Open Boundary Condition
  !!----------------------------------------------------------------------
  !! history :
  !!  8.0   01/91   (CLIPPER)  Original code
  !!  8.5   06/02   (C. Talandier)  modules
  !!        06/04   (F. Durand) ORCA2_ZIND config
  !!        06/04   (F. Durand) Dimensions of arrays vsdta, tsdta, ssdta,
  !!                vndta, tndta, sndta, uwdta, twdta, swdta, uedta, tedta, sedta
  !!                are defined to the actual dimensions of the OBs i.e.
  !!		     (jpisd:jpisf,jpk,jptobc) for the South OB
  !! 		     (jpind:jpinf,jpk,jptobc) for the North OB
  !!		     (jpjwd:jpjwf,jpk,jptobc) for the West OB
  !! 		     (jpjed:jpjef,jpk,jptobc) for the East OB
  !!
  !! 2007.12.07      D.Brickman 
  !!                 Dimensions as listed above are potentially incorrect
  !!                 as looping indices in the routines that use these arrays
  !!                 can start at a number less than the starting index
  !!                 e.g. jpjed = 2, but loop can start at 1 ----> ERROR
  !!                 Solution: set starting index = 1 
  !!		     (1:jpisf,jpk,jptobc) for the South OB
  !! 		     (1:jpinf,jpk,jptobc) for the North OB 
  !!		     (1:jpjwf,jpk,jptobc) for the West OB
  !! 		     (1:jpjef,jpk,jptobc) for the East OB 
  !!
  !!----------------------------------------------------------------------
  !! * Modules used
  USE par_oce         ! ocean parameters
  USE obc_par         ! open boundary condition parameters

  IMPLICIT NONE
  PUBLIC
  !!----------------------------------------------------------------------
  !!  OPA 9.0 , LOCEAN-IPSL (2005)
  !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_oce.F90,v 1.5 2005/12/28 09:25:07 opalod Exp $
  !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
  !!----------------------------------------------------------------------

  !!----------------------------------------------------------------------
  !! open boundary variables
  !!----------------------------------------------------------------------
  !!
  !!General variables for open boundaries:
  !!--------------------------------------
  INTEGER ::              & !: * namelist ??? *
     nbobc    = 3  ,      & !: number of open boundaries ( 1=< nbobc =< 4 ) !byoung
     nobc_dta = 0           !:  = 0 use the initial state as obc data
     !                      !   = 1 read obc data in obcxxx.dta files

  LOGICAL ::  ln_obc_clim = .true.  !:  obc data files are climatological
  LOGICAL ::  ln_obc_fla  = .false. !:  Flather open boundary condition not used
  LOGICAL ::  ln_vol_cst  = .true.  !:  Conservation of the whole volume

  REAL(wp) ::             & !!: open boundary namelist (namobc)
     rdpein =  1.  ,      &  !: damping time scale for inflow at East open boundary
     rdpwin =  1.  ,      &  !:    "                      "   at West open boundary
     rdpsin =  1.  ,      &  !:    "                      "   at South open boundary
     rdpnin =  1.  ,      &  !:    "                      "   at North open boundary
     rdpeob = 15.  ,      &  !: damping time scale for the climatology at East open boundary
     rdpwob = 15.  ,      &  !:    "                           "       at West open boundary
     rdpsob = 15.  ,      &  !:    "                           "       at South open boundary
     rdpnob = 15.  ,      &  !:    "                           "       at North open boundary
     volemp =  1.            !: = 0 the total volume will have the variability of the
                             !      surface Flux E-P else (volemp = 1) the volume will be constant
                             !  = 1 the volume will be constant during all the integration.

  LOGICAL ::              &  !:
     lfbceast, lfbcwest,  &  !: logical flag for a fixed East and West open boundaries
     lfbcnorth, lfbcsouth    !: logical flag for a fixed North and South open boundaries
     !                       !  These logical flags are set to 'true' if damping time
     !                       !  scale are set to 0 in the namelist, for both inflow and outflow).

  REAL(wp), DIMENSION(jpi,jpj) :: &  !:
     obctmsk                !: mask array identical to tmask, execpt along OBC where it is set to 0
     !                      !  it used to calculate the cumulate flux E-P in the obcvol.F90 routine

  !!----------------
  !! Rigid lid case:
  !!----------------
  INTEGER ::   nbic !: number of isolated coastlines ( 0 <= nbic <= 3 )

  INTEGER, DIMENSION(jpnic,0:4,3) ::   &  !:
     miic, mjic     !: position of isolated coastlines points

  INTEGER, DIMENSION(0:4,3) ::   &  !:
     mnic           !: number of points on isolated coastlines

  REAL(wp), DIMENSION(jpi,jpj) ::   &  !:
     gcbob          !: right hand side of the barotropic elliptic equation associated
     !              !  with the OBC

  REAL(wp), DIMENSION(jpi,jpj,3) ::   &  !:
     gcfobc         !: coef. associated with the contribution of isolated coastlines
     !              !  to the right hand side of the barotropic elliptic equation

  REAL(wp), DIMENSION(3) ::   &  !:
     gcbic          !: time variation of the barotropic stream function along the
     !              !  isolated coastlines

  REAL(wp), DIMENSION(1) ::   &  !:
     bsfic0         !: barotropic stream function on isolated coastline

  REAL(wp), DIMENSION(3) ::   &  !:
     bsfic          !: barotropic stream function on isolated coastline

!!DB 2008.05.06 global masking arrays for obcdta calcs
!!Arrays computed in obcini
!!del the first one later
  REAL(wp) :: emaskg(jpjdta)=0.0, wmaskg(jpjdta)=0.0, smaskg(jpidta)=0.0, nmaskg(jpidta)=0.0
  REAL(wp) :: emaskg2(jpjdta,jpk)=0.0, wmaskg2(jpjdta,jpk)=0.0, & 
       smaskg2(jpidta,jpk)=0.0, nmaskg2(jpidta,jpk)=0.0



  !!--------------------
  !! East open boundary:
  !!--------------------
  INTEGER ::   nie0  , nie1      !: do loop index in mpp case for jpieob
  INTEGER ::   nie0p1, nie1p1    !: do loop index in mpp case for jpieob+1
  INTEGER ::   nie0m1, nie1m1    !: do loop index in mpp case for jpieob-1
  INTEGER ::   nje0  , nje1      !: do loop index in mpp case for jpjed, jpjef
  INTEGER ::   nje0p1, nje1m1    !: do loop index in mpp case for jpjedp1,jpjefm1
  INTEGER ::   nje1m2, nje0m1    !: do loop index in mpp case for jpjefm1-1,jpjed

  REAL(wp), DIMENSION(jpj) ::    &  !:
     bsfeob              !: now barotropic stream fuction computed at the OBC. The corres-
     !                   ! ponding bsfn will be computed by the forward time step in dynspg.

  REAL(wp), DIMENSION(jpj,3,3) ::   &  !:
     bebnd               !: east boundary barotropic streamfunction over 3 rows
     !                   ! and 3 time step (now, before, and before before)

!DB
  REAL(wp), DIMENSION(1:jpjef) ::   &  !:
     bfoe,             & !: now climatology of the east boundary barotropic stream function
     sshfoe,           & !: now climatology of the east boundary sea surface height
     eta_e,      & !: AD3Aug07:  elev. geostropic elev computed from uedta
     ubtfoe,vbtfoe       !: now climatology of the east boundary barotropic transport

  REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
     ufoe, vfoe,       & !: now climatology of the east boundary velocities
!     tfoe, sfoe,       & !: now climatology of the east boundary temperature and salinity
     uclie               !: baroclinic componant of the zonal velocity after radiation
     !                   ! in the obcdyn.F90 routine
!byoung
  REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
     tfoe, sfoe        !: now climatology of the east boundary temperature and salinity
     
!DB jpjed <---- 1
  REAL(wp), DIMENSION(1:jpjef,jpj) ::   &  !:
     sshfoe_b            !: east boundary ssh correction averaged over the barotropic loop
                         !: (if Flather's algoritm applied at open boundary)

!DB jpjed <---- 1
  REAL(wp), DIMENSION(1:jpjef,0:jptobc+1) ::   &  !:
     sshedta, ubtedta    !: array used for interpolating monthly data on the east boundary

!DB jpjed <---- 1
  REAL(wp), DIMENSION(1:jpjef,5) ::    &    !:
     tidesshemag, tidesshepha, tidevbtemag, tidevbtepha    

!DB jpjed <---- 1
  REAL(wp), DIMENSION(1:jpjef,jpk,jptobc) ::   &  !:
     uedta 

!AD
!DB jpjed <---- 1
  REAL(wp), DIMENSION(1:jpjef,jpk) ::   &  !:
     uedta1 !: array used for interpolating monthly data on the east boundary

!DB jpjed <---- 1
  REAL(wp), DIMENSION(1:jpjef,jpk,5,jptobc) ::   &  !:
     tedta, sedta, tedta1, sedta1 !: array used for interpolating monthly data on the east boundary

  !!-------------------------------
  !! Arrays for radiative East OBC:
  !!-------------------------------
  REAL(wp), DIMENSION(jpj,jpk,3,3) ::   &  !:
     uebnd, vebnd                  !: baroclinic u & v component of the velocity over 3 rows
                                   ! and 3 time step (now, before, and before before)

  REAL(wp), DIMENSION(jpj,jpk,2,2) ::   &  !:
     tebnd, sebnd                  !: East boundary temperature and salinity over 2 rows
                                   ! and 2 time step (now and before)

  REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
     u_cxebnd, v_cxebnd            !: Zonal component of the phase speed ratio computed with
                                   ! radiation of u and v velocity (respectively) at the
                                   ! east open boundary (u_cxebnd = cx rdt )

  REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
     uemsk, vemsk, temsk           !: 2D mask for the East OB

!byoung
  REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
     temsk5

  ! Note that those arrays are optimized for mpp case
  ! (hence the dimension jpj is the size of one processor subdomain)

  !!--------------------
  !! West open boundary
  !!--------------------
  INTEGER ::   niw0  , niw1       !: do loop index in mpp case for jpiwob
  INTEGER ::   niw0p1, niw1p1     !: do loop index in mpp case for jpiwob+1
  INTEGER ::   njw0  , njw1       !: do loop index in mpp case for jpjwd, jpjwf
  INTEGER ::   njw0p1, njw1m1     !: do loop index in mpp case for jpjwdp1,jpjwfm1
  INTEGER ::   njw1m2, njw0m1     !: do loop index in mpp case for jpjwfm2,jpjwd

  REAL(wp), DIMENSION(jpj) ::   &  !:
     bsfwob              !: now barotropic stream fuction computed at the OBC. The corres-
     !                   !  ponding bsfn will be computed by the forward time step in dynspg.

  REAL(wp), DIMENSION(jpj,3,3) ::   &  !:
     bwbnd               !: West boundary barotropic streamfunction over
     !                   !  3 rows and 3 time step (now, before, and before before)

!DB
  REAL(wp), DIMENSION(1:jpjwf) ::   &  !:
     bfow,             & !: now climatology of the west boundary barotropic stream function
     sshfow,           & !: now climatology of the west boundary sea surface height
     eta_w,       & !: elevation computed from uwdta
     ubtfow,vbtfow       !: now climatology of the west boundary barotropic transport
     

  REAL(wp), DIMENSION(jpj,jpk) ::   &  !:
     ufow, vfow,       & !: now climatology of the west velocities
!     tfow, sfow,       & !: now climatology of the west temperature and salinity
     ucliw               !: baroclinic componant of the zonal velocity after the radiation
     !                   !  in the obcdyn.F90 routine

!byoung
  REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
     tfow, sfow        !: now climatology of the west temperature and salinity
     
  REAL(wp), DIMENSION(1:jpjwf,jpj) ::   &  !:
     sshfow_b            !: west boundary ssh correction averaged over the barotropic loop
                         !: (if Flather's algoritm applied at open boundary)

  REAL(wp), DIMENSION(1:jpjwf,0:jptobc+1) ::   &  !:
     sshwdta, ubtwdta    !: array used for interpolating monthly data on the west boundary

!ylu
!  REAL(wp), DIMENSION(1:jpjwf,4) ::    &    !:
  REAL(wp), DIMENSION(1:jpjwf,5) ::    &    !:
     tidesshwmag, tidesshwpha, tidevbtwmag, tidevbtwpha    

  REAL(wp), DIMENSION(1:jpjwf,jpk,jptobc) ::   &  !:
     uwdta !, uwdta1 !: array used for interpolating monthly data on the west boundary
!     uwdta, twdta, swdta !: array used for interpolating monthly data on the west boundary

!DBG
  REAL(wp), DIMENSION(1:jpjwf,jpk) ::   &  !:
     uwdta1 !: array used for interpolating monthly data on the west boundary


!byoung
  REAL(wp), DIMENSION(1:jpjwf,jpk,5,jptobc) ::   &  !:
     twdta, swdta, twdta1, swdta1 !: array used for interpolating monthly data on the west boundary

  !!-------------------------------
  !! Arrays for radiative West OBC
  !!-------------------------------
  REAL(wp), DIMENSION(jpj,jpk,3,3) ::   &  !:
     uwbnd, vwbnd                  !: baroclinic u & v components of the velocity over 3 rows
     !                             !  and 3 time step (now, before, and before before)

  REAL(wp), DIMENSION(jpj,jpk,2,2) ::   &  !:
     twbnd, swbnd                  !: west boundary temperature and salinity over 2 rows and
     !                             !  2 time step (now and before)

  REAL(wp), DIMENSION(jpj,jpk) ::    &  !:
     u_cxwbnd, v_cxwbnd            !: Zonal component of the phase speed ratio computed with
     !                             !  radiation of zonal and meridional velocity (respectively)
     !                             !  at the west open boundary (u_cxwbnd = cx rdt )

  REAL(wp), DIMENSION(jpj,jpk) ::    &  !:
     uwmsk, vwmsk, twmsk           !: 2D mask for the West OB

!byoung
  REAL(wp), DIMENSION(jpj,jpk,5) ::   &  !:
     twmsk5

  ! Note that those arrays are optimized for mpp case
  ! (hence the dimension jpj is the size of one processor subdomain)

  !!---------------------
  !! North open boundary
  !!---------------------
  INTEGER ::   nin0  , nin1       !: do loop index in mpp case for jpind, jpinf
  INTEGER ::   nin0p1, nin1m1     !: do loop index in mpp case for jpindp1, jpinfm1
  INTEGER ::   nin1m2, nin0m1     !: do loop index in mpp case for jpinfm1-1,jpind
  INTEGER ::   njn0  , njn1       !: do loop index in mpp case for jpnob
  INTEGER ::   njn0p1, njn1p1     !: do loop index in mpp case for jpnob+1
  INTEGER ::   njn0m1, njn1m1     !: do loop index in mpp case for jpnob-1

  REAL(wp), DIMENSION(jpi) ::   &  !:
     bsfnob              !: now barotropic stream fuction computed at the OBC. The corres-
     !                   !  ponding bsfn will be computed by the forward time step in dynspg.

  REAL(wp), DIMENSION(jpi,3,3) ::   &  !:
     bnbnd               !: north boundary barotropic streamfunction over
     !                   !  3 rows and 3 time step (now, before, and before before)

  REAL(wp), DIMENSION(1:jpinf) ::   &  !:
     bfon,             & !: now climatology of the north boundary barotropic stream function
     sshfon,           & !: now climatology of the north boundary sea surface height
     ubtfon,vbtfon       !: now climatology of the north boundary barotropic transport

  REAL(wp), DIMENSION(jpi,jpk) ::   &    !:
     ufon, vfon,       & !: now climatology of the north boundary velocities
!      tfon, sfon,       & !: now climatology of the north boundary temperature and salinity
     vclin               !: baroclinic componant of the meridian velocity after the radiation
     !                   !  in yhe obcdyn.F90 routine

!sujie
  REAL(wp), DIMENSION(jpi,jpk,5) ::   &    !:
     tfon, sfon        !: now climatology of the north boundary temperature and salinity

  REAL(wp), DIMENSION(1:jpinf,jpj) ::   &  !:
     sshfon_b            !: north boundary ssh correction averaged over the barotropic loop
                         !: (if Flather's algoritm applied at open boundary)

  REAL(wp), DIMENSION(1:jpinf,0:jptobc+1) ::   &  !:
     sshndta, vbtndta   !: array used for interpolating monthly data on the north boundary

!ylu
  REAL(wp), DIMENSION(1:jpinf,5) ::    &    !:
     tidesshnmag, tidesshnpha, tidevbtnmag, tidevbtnpha    

  REAL(wp), DIMENSION(1:jpinf,jpk,jptobc) ::   &  !:
     vndta !sujie, tndta, sndta !: array used for interpolating monthly data on the north boundary
!sujie
  REAL(wp), DIMENSION(1:jpinf,jpk,5,jptobc) ::   &  !:
     tndta, sndta !: array used for interpolating monthly data on the north boundary

  !!--------------------------------
  !! Arrays for radiative North OBC
  !!--------------------------------
  !!
  REAL(wp), DIMENSION(jpi,jpk,3,3) ::   &   !:
     unbnd, vnbnd                  !: baroclinic u & v components of the velocity over 3
     !                             !  rows and 3 time step (now, before, and before before)

  REAL(wp), DIMENSION(jpi,jpk,2,2) ::   &   !:
     tnbnd, snbnd                  !: north boundary temperature and salinity over
     !                             !  2 rows and 2 time step (now and before)

  REAL(wp), DIMENSION(jpi,jpk) ::   &     !:
     u_cynbnd, v_cynbnd            !: Meridional component of the phase speed ratio compu-
     !                             !  ted with radiation of zonal and meridional velocity
     !                             !  (respectively) at the north OB (u_cynbnd = cx rdt )

  REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
     unmsk, vnmsk, tnmsk           !: 2D mask for the North OB

!sujie
  REAL(wp), DIMENSION(jpi,jpk,5) ::   &  !:
     tnmsk5

  ! Note that those arrays are optimized for mpp case
  ! (hence the dimension jpj is the size of one processor subdomain)

  !!---------------------
  !! South open boundary
  !!---------------------
  INTEGER ::   nis0  , nis1       !: do loop index in mpp case for jpisd, jpisf
  INTEGER ::   nis0p1, nis1m1     !: do loop index in mpp case for jpisdp1, jpisfm1
  INTEGER ::   nis1m2, nis0m1     !: do loop index in mpp case for jpisfm1-1,jpisd
  INTEGER ::   njs0  , njs1       !: do loop index in mpp case for jpsob
  INTEGER ::   njs0p1, njs1p1     !: do loop index in mpp case for jpsob+1

  REAL(wp), DIMENSION(jpi) ::    &   !:
     bsfsob              !: now barotropic stream fuction computed at the OBC.The corres-
     !                   !  ponding bsfn will be computed by the forward time step in dynspg.
  REAL(wp), DIMENSION(jpi,3,3) ::   &   !:
     bsbnd               !: south boundary barotropic stream function over
     !                   !  3 rows and 3 time step (now, before, and before before)

!DB
  REAL(wp), DIMENSION(1:jpisf) ::    &   !:
     bfos,             & !: now climatology of the south boundary barotropic stream function
     sshfos,           & !: now climatology of the south boundary sea surface height
     eta_s,            & !: AD3Aug07:  elev. geostropic elev computed from vsdta
     ubtfos,vbtfos       !: now climatology of the south boundary barotropic transport

  REAL(wp), DIMENSION(jpi,jpk) ::    &   !:
     ufos, vfos ,       & !: now climatology of the south boundary velocities
!      tfos, sfos,       & !: now climatology of the south boundary temperature and salinity
     vclis               !: baroclinic componant of the meridian velocity after the radiation
     !                   !  in the obcdyn.F90 routine

!sujie
  REAL(wp), DIMENSION(jpi,jpk,5) ::    &   !:
     tfos, sfos        !: now climatology of the south boundary temperature and salinity

!DB
  REAL(wp), DIMENSION(1:jpisf,jpj) ::   &  !:
     sshfos_b            !: south boundary ssh correction averaged over the barotropic loop
                         !: (if Flather's algoritm applied at open boundary)

  REAL(wp), DIMENSION(1:jpisf,0:jptobc+1) ::    &    !:
     sshsdta, vbtsdta    !: array used for interpolating monthly data on the south boundary

!ylu
  REAL(wp), DIMENSION(1:jpisf,5) ::    &    !:
     tidesshsmag, tidesshspha, tidevbtsmag, tidevbtspha    

  REAL(wp), DIMENSION(1:jpisf,jpk,jptobc) ::    &    !:
     vsdta !sujie, tsdta, ssdta   !: array used for interpolating monthly data on the south boundary

!DB
  REAL(wp), DIMENSION(1:jpisf,jpk) ::    &    !:
     vsdta1  !AD: array used for interpolatated global obc array
     
!sujie
  REAL(wp), DIMENSION(1:jpisf,jpk,5,jptobc) ::    &    !:
     tsdta, ssdta, tsdta1, ssdta1   !: array used for interpolating monthly data on the south boundary

  !!--------------------------------
  !! Arrays for radiative South OBC
  !!--------------------------------
  !!                        computed by the forward time step in dynspg.
  REAL(wp), DIMENSION(jpi,jpk,3,3) ::   &   !:
     usbnd, vsbnd                  !: baroclinic u & v components of the velocity over 3
     !                             !  rows and 3 time step (now, before, and before before)

  REAL(wp), DIMENSION(jpi,jpk,2,2) ::   &  !:
     tsbnd, ssbnd                  !: south boundary temperature and salinity over
     !                             !  2 rows and 2 time step (now and before)

  REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
     u_cysbnd, v_cysbnd            !: Meridional component of the phase speed ratio compu-
     !                             !  ted with radiation of zonal and meridional velocity
     !                             !  (repsectively) at the south OB (u_cynbnd = cx rdt )

  REAL(wp), DIMENSION(jpi,jpk) ::   &  !:
     usmsk, vsmsk, tsmsk           !: 2D mask for the South OB

!sujie
  REAL(wp), DIMENSION(jpi,jpk,5) ::   &  !:
     tsmsk5

  ! Note that those arrays are optimized for mpp case
  ! (hence the dimension jpj is the size of one processor subdomain)

#else
  !!----------------------------------------------------------------------
  !!   Default option :                                       Empty module
  !!----------------------------------------------------------------------
#endif

  !!======================================================================
END MODULE obc_oce
