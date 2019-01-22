


PGF90 (Version     10.5)          09/23/2010  09:49:40      page 1

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

Filename: dynvor.F90

(    1) MODULE dynvor
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynvor  ***
(    4)    !! Ocean dynamics: Update the momentum trend with the relative and
(    5)    !!                 planetary vorticity trends
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   dyn_vor_enstrophy: enstrophy conserving scheme       (ln_dynvor_ens=T)
(   10)    !!   dyn_vor_energy   : energy conserving scheme          (ln_dynvor_ene=T)
(   11)    !!   dyn_vor_mixed    : mixed enstrophy/energy conserving (ln_dynvor_mix=T)
(   12)    !!   dyn_vor_ene_ens  : energy and enstrophy conserving   (ln_dynvor_een=T)
(   13)    !!   dyn_vor_ctl      : control of the different vorticity option
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE trdmod          ! ocean dynamics trends 
(   20)    USE trdmod_oce      ! ocean variables trends
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22)    USE prtctl          ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC dyn_vor_enstrophy      ! routine called by step.F90
(   29)    PUBLIC dyn_vor_energy         ! routine called by step.F90
(   30)    PUBLIC dyn_vor_mixed          ! routine called by step.F90
(   31)    PUBLIC dyn_vor_ene_ens        ! routine called by step.F90
(   32)    PUBLIC dyn_vor_ctl            ! routine called by step.F90
(   33) 
(   34)    !! * Shared module variables
(   35)    LOGICAL, PUBLIC ::   ln_dynvor_ene = .FALSE.   !: energy conserving scheme
(   36)    LOGICAL, PUBLIC ::   ln_dynvor_ens = .TRUE.    !: enstrophy conserving scheme
(   37)    LOGICAL, PUBLIC ::   ln_dynvor_mix = .FALSE.   !: mixed scheme
(   38)    LOGICAL, PUBLIC ::   ln_dynvor_een = .FALSE.    !: energy and enstrophy conserving scheme
(   39) 
(   40)    !! * Substitutions
(   41) #  include "domzgr_substitute.h90"
(   42) #  include "vectopt_loop_substitute.h90"
(   43)    !!----------------------------------------------------------------------
(   44)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   45)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynvor.F90,v 1.11 2006/03/21 08:26:27 opalod Exp $ 
(   46)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 2

(   47)    !!----------------------------------------------------------------------
(   48) 
(   49) CONTAINS
(   50) 
(   51)    SUBROUTINE dyn_vor_energy( kt )
(   52)       !!----------------------------------------------------------------------
(   53)       !!                  ***  ROUTINE dyn_vor_energy  ***
(   54)       !!
(   55)       !! ** Purpose :   Compute the now total vorticity trend and add it to 
(   56)       !!      the general trend of the momentum equation.
(   57)       !!
(   58)       !! ** Method  :   Trend evaluated using now fields (centered in time) 
(   59)       !!      and the Sadourny (1975) flux form formulation : conserves the
(   60)       !!      horizontal kinetic energy.
(   61)       !!      The trend of the vorticity term is given by:
(   62)       !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
(   63)       !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f  mi(e1v*e3v vn) ]
(   64)       !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f  mj(e2u*e3u un) ]
(   65)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(   66)       !!          voru = 1/e1u  mj-1[ (rotn+f)  mi(e1v vn) ]
(   67)       !!          vorv = 1/e2v  mi-1[ (rotn+f)  mj(e2u un) ]
(   68)       !!      Add this trend to the general momentum trend (ua,va):
(   69)       !!          (ua,va) = (ua,va) + ( voru , vorv )
(   70)       !!
(   71)       !! ** Action : - Update (ua,va) with the now vorticity term trend
(   72)       !!             - save the trends in (utrd,vtrd) in 2 parts (relative
(   73)       !!               and planetary vorticity trends) ('key_trddyn')
(   74)       !!
(   75)       !! References :
(   76)       !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
(   77)       !! History :
(   78)       !!   5.0  !  91-11  (G. Madec)  Original code
(   79)       !!   6.0  !  96-01  (G. Madec)  s-coord, suppress work arrays
(   80)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   81)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   82)       !!----------------------------------------------------------------------
(   83)       !! * Modules used     
(   84)       USE oce, ONLY :    ztdua => ta,     & ! use ta as 3D workspace   
(   85)                          ztdva => sa        ! use sa as 3D workspace   
(   86) 
(   87)       !! * Arguments
(   88)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   89) 
(   90)       !! * Local declarations
(   91)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(   92)       REAL(wp) ::   &
(   93)          zfact2, zua, zva,               &  ! temporary scalars
(   94)          zx1, zx2, zy1, zy2                 !    "         "
(   95)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   96)          zwx, zwy, zwz                      ! temporary workspace
(   97)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   98)          zcu, zcv                           !    "         "
(   99)       !!----------------------------------------------------------------------
(  100) 
(  101)       IF( kt == nit000 ) THEN
(  102)          IF(lwp) WRITE(numout,*)
(  103)          IF(lwp) WRITE(numout,*) 'dyn_vor_energy : vorticity term: energy conserving scheme'
(  104)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 3

(  105)       ENDIF
(  106) 
(  107)       ! Local constant initialization
(  108)       zfact2 = 0.5 * 0.5
(  109)       zcu (:,:,:) = 0.e0
(  110)       zcv (:,:,:) = 0.e0
(  111) 
(  112)       ! Save ua and va trends
(  113)       IF( l_trddyn )   THEN
(  114)          ztdua(:,:,:) = ua(:,:,:) 
(  115)          ztdva(:,:,:) = va(:,:,:) 
(  116)          zcu(:,:,:) = 0.e0
(  117)          zcv(:,:,:) = 0.e0
(  118)       ENDIF
(  119)       
(  120)       !                                                ! ===============
(  121)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  122)          !                                             ! ===============
(  123)          
(  124)          ! Potential vorticity and horizontal fluxes
(  125)          ! -----------------------------------------
(  126)          IF( lk_sco ) THEN
(  127)             zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) / fse3f(:,:,jk)
(  128)             zwx(:,:) = e2u(:,:) * fse3u(:,:,jk) * un(:,:,jk)
(  129)             zwy(:,:) = e1v(:,:) * fse3v(:,:,jk) * vn(:,:,jk)
(  130)          ELSE
(  131)             zwz(:,:) = rotn(:,:,jk) + ff(:,:)
(  132)             zwx(:,:) = e2u(:,:) * un(:,:,jk)
(  133)             zwy(:,:) = e1v(:,:) * vn(:,:,jk)
(  134)          ENDIF
(  135) 
(  136)          ! Compute and add the vorticity term trend
(  137)          ! ----------------------------------------
(  138)          DO jj = 2, jpjm1
(  139)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  140)                zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
(  141)                zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
(  142)                zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
(  143)                zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
(  144)                zua = zfact2 / e1u(ji,jj) * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
(  145)                zva =-zfact2 / e2v(ji,jj) * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
(  146)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  147)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  148)             END DO  
(  149)          END DO  
(  150)          !                                             ! ===============
(  151)       END DO                                           !   End of slab
(  152)       !                                                ! ===============
(  153) 
(  154)       ! save the relative & planetary vorticity trends for diagnostic
(  155)       ! momentum trends
(  156)       IF( l_trddyn )   THEN
(  157)          ! Compute the planetary vorticity term trend
(  158)          !                                                ! ===============
(  159)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  160)             !                                             ! ===============
(  161)             DO jj = 2, jpjm1
(  162)                DO ji = fs_2, fs_jpim1   ! vector opt.






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 4

(  163)                   zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
(  164)                   zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
(  165)                   zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
(  166)                   zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
(  167) # if defined key_s_coord
(  168)                  zcu(ji,jj,jk) = zfact2 / e1u(ji,jj) * ( ff(ji  ,jj-1) / fse3f(ji,jj-1,jk) * zy1  &
(  169)                    &                                   + ff(ji  ,jj  ) / fse3f(ji,jj  ,jk) * zy2 )
(  170)                  zcv(ji,jj,jk) =-zfact2 / e2v(ji,jj) * ( ff(ji-1,jj  ) / fse3f(ji-1,jj,jk) * zx1  &
(  171)                    &                                   + ff(ji  ,jj  ) / fse3f(ji  ,jj,jk) * zx2 )
(  172) # else
(  173)                  zcu(ji,jj,jk) = zfact2 / e1u(ji,jj) * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
(  174)                  zcv(ji,jj,jk) =-zfact2 / e2v(ji,jj) * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
(  175) # endif
(  176)                END DO  
(  177)             END DO  
(  178)             !                                             ! ===============
(  179)          END DO                                           !   End of slab
(  180)          !                                                ! ===============
(  181) 
(  182)          ! Compute the relative vorticity term trend
(  183)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
(  184)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 
(  185) 
(  186)          CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
(  187)          CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
(  188)          CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
(  189)       ENDIF
(  190) 
(  191)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  192)          CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
(  193)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  194)       ENDIF
(  195) 
(  196)    END SUBROUTINE dyn_vor_energy
(  197) 
(  198) 
(  199)    SUBROUTINE dyn_vor_mixed( kt )
(  200)       !!----------------------------------------------------------------------
(  201)       !!                 ***  ROUTINE dyn_vor_mixed  ***
(  202)       !!
(  203)       !! ** Purpose :   Compute the now total vorticity trend and add it to
(  204)       !!      the general trend of the momentum equation.
(  205)       !!
(  206)       !! ** Method  :   Trend evaluated using now fields (centered in time)
(  207)       !!      Mixte formulation : conserves the potential enstrophy of a hori-
(  208)       !!      zontally non-divergent flow for (rotzu x uh), the relative vor-
(  209)       !!      ticity term and the horizontal kinetic energy for (f x uh), the
(  210)       !!      coriolis term. the now trend of the vorticity term is given by:
(  211)       !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
(  212)       !!          voru = 1/e1u  mj-1(rotn/e3f) mj-1[ mi(e1v*e3v vn) ]
(  213)       !!              +1/e1u  mj-1[ f/e3f          mi(e1v*e3v vn) ]
(  214)       !!          vorv = 1/e2v  mi-1(rotn/e3f) mi-1[ mj(e2u*e3u un) ]
(  215)       !!              +1/e2v  mi-1[ f/e3f          mj(e2u*e3u un) ]
(  216)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(  217)       !!          voru = 1/e1u  mj-1(rotn) mj-1[ mi(e1v vn) ]
(  218)       !!              +1/e1u  mj-1[ f          mi(e1v vn) ]
(  219)       !!          vorv = 1/e2v  mi-1(rotn) mi-1[ mj(e2u un) ]
(  220)       !!              +1/e2v  mi-1[ f          mj(e2u un) ]






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 5

(  221)       !!      Add this now trend to the general momentum trend (ua,va):
(  222)       !!          (ua,va) = (ua,va) + ( voru , vorv )
(  223)       !!
(  224)       !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
(  225)       !!             - Save the trends in (utrd,vtrd) in 2 parts (relative
(  226)       !!               and planetary vorticity trends) ('key_trddyn')
(  227)       !!
(  228)       !! References :
(  229)       !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
(  230)       !! History :
(  231)       !!   5.0  !  91-11  (G. Madec) Original code, enstrophy-energy-combined schemes
(  232)       !!   6.0  !  96-01  (G. Madec) s-coord, suppress work arrays
(  233)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  234)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(  235)       !!----------------------------------------------------------------------
(  236)       !! * Modules used     
(  237)       USE oce, ONLY :    ztdua => ta,     & ! use ta as 3D workspace   
(  238)                          ztdva => sa        ! use sa as 3D workspace   
(  239)       !! * Arguments
(  240)       INTEGER, INTENT( in ) ::   kt         ! ocean timestep index
(  241) 
(  242)       !! * Local declarations
(  243)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(  244)       REAL(wp) ::   &
(  245)          zfact1, zfact2, zua, zva,       &  ! temporary scalars
(  246)          zcua, zcva, zx1, zx2, zy1, zy2
(  247)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  248)          zwx, zwy, zwz, zww                 ! temporary workspace
(  249)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  250)          zcu, zcv                           !    "         "
(  251)       !!----------------------------------------------------------------------
(  252) 
(  253)       IF( kt == nit000 ) THEN
(  254)          IF(lwp) WRITE(numout,*)
(  255)          IF(lwp) WRITE(numout,*) 'dyn_vor_mixed : vorticity term: mixed energy/enstrophy conserving scheme'
(  256)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
(  257)       ENDIF
(  258) 
(  259)       ! Local constant initialization
(  260)       zfact1 = 0.5 * 0.25
(  261)       zfact2 = 0.5 * 0.5
(  262) 
(  263)       ! Save ua and va trends
(  264)       IF( l_trddyn )   THEN
(  265)          ztdua(:,:,:) = ua(:,:,:) 
(  266)          ztdva(:,:,:) = va(:,:,:) 
(  267)          zcu(:,:,:) = 0.e0
(  268)          zcv(:,:,:) = 0.e0
(  269)       ENDIF
(  270) 
(  271)       !                                                ! ===============
(  272)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  273)          !                                             ! ===============
(  274) 
(  275)          ! Relative and planetary potential vorticity and horizontal fluxes
(  276)          ! ----------------------------------------------------------------
(  277)          IF( lk_sco ) THEN        
(  278)             zwz(:,:) = ff  (:,:)    / fse3f(:,:,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 6

(  279)             zww(:,:) = rotn(:,:,jk) / fse3f(:,:,jk)
(  280)             zwx(:,:) = e2u(:,:) * fse3u(:,:,jk) * un(:,:,jk)
(  281)             zwy(:,:) = e1v(:,:) * fse3v(:,:,jk) * vn(:,:,jk)
(  282)          ELSE
(  283)             zwz(:,:) = ff(:,:)
(  284)             zww(:,:) = rotn(:,:,jk)
(  285)             zwx(:,:) = e2u(:,:) * un(:,:,jk)
(  286)             zwy(:,:) = e1v(:,:) * vn(:,:,jk)
(  287)          ENDIF
(  288) 
(  289)          ! Compute and add the vorticity term trend
(  290)          ! ----------------------------------------
(  291)          DO jj = 2, jpjm1
(  292)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  293)                zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  294)                zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  295)                zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  296)                zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  297)                ! enstrophy conserving formulation for relative vorticity term
(  298)                zua = zfact1 * ( zww(ji  ,jj-1) + zww(ji,jj) ) * ( zy1 + zy2 )
(  299)                zva =-zfact1 * ( zww(ji-1,jj  ) + zww(ji,jj) ) * ( zx1 + zx2 )
(  300)                ! energy conserving formulation for planetary vorticity term
(  301)                zcua = zfact2 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
(  302)                zcva =-zfact2 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
(  303) 
(  304)                ua(ji,jj,jk) = ua(ji,jj,jk) + zcua + zua
(  305)                va(ji,jj,jk) = va(ji,jj,jk) + zcva + zva
(  306)             END DO  
(  307)          END DO  
(  308)          !                                             ! ===============
(  309)       END DO                                           !   End of slab
(  310)       !                                                ! ===============
(  311) 
(  312)       ! save the relative & planetary vorticity trends for diagnostic
(  313)       ! momentum trends
(  314)       IF( l_trddyn )   THEN
(  315)          ! Compute the planetary vorticity term trend
(  316)          !                                                ! ===============
(  317)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  318)             !                                             ! ===============
(  319)             DO jj = 2, jpjm1
(  320)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  321)                   zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  322)                   zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  323)                   zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  324)                   zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  325)                   ! energy conserving formulation for planetary vorticity term
(  326)                   zcu(ji,jj,jk) = zfact2 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
(  327)                   zcv(ji,jj,jk) =-zfact2 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
(  328)                END DO  
(  329)             END DO  
(  330)             !                                             ! ===============
(  331)          END DO                                           !   End of slab
(  332)          !                                                ! ===============
(  333) 
(  334)          ! Compute the relative vorticity term trend
(  335)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
(  336)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 7

(  337) 
(  338)          CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
(  339)          CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
(  340)          CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
(  341)       ENDIF
(  342) 
(  343)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  344)          CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
(  345)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  346)       ENDIF
(  347) 
(  348)    END SUBROUTINE dyn_vor_mixed
(  349) 
(  350) 
(  351)    SUBROUTINE dyn_vor_enstrophy( kt )
(  352)       !!----------------------------------------------------------------------
(  353)       !!                ***  ROUTINE dyn_vor_enstrophy  ***
(  354)       !!
(  355)       !! ** Purpose :   Compute the now total vorticity trend and add it to
(  356)       !!      the general trend of the momentum equation.
(  357)       !!
(  358)       !! ** Method  :   Trend evaluated using now fields (centered in time)
(  359)       !!      and the Sadourny (1975) flux FORM formulation : conserves the
(  360)       !!      potential enstrophy of a horizontally non-divergent flow. the
(  361)       !!      trend of the vorticity term is given by:
(  362)       !!       * s-coordinate (lk_sco=T), the e3. are inside the derivative:
(  363)       !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f ]  mj-1[ mi(e1v*e3v vn) ]
(  364)       !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f ]  mi-1[ mj(e2u*e3u un) ]
(  365)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(  366)       !!          voru = 1/e1u  mj-1[ rotn+f ]  mj-1[ mi(e1v vn) ]
(  367)       !!          vorv = 1/e2v  mi-1[ rotn+f ]  mi-1[ mj(e2u un) ]
(  368)       !!      Add this trend to the general momentum trend (ua,va):
(  369)       !!          (ua,va) = (ua,va) + ( voru , vorv )
(  370)       !!
(  371)       !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
(  372)       !!             - Save the trends in (utrd,vtrd) in 2 parts (relative 
(  373)       !!               and planetary vorticity trends) ('key_trddyn')
(  374)       !!
(  375)       !! References :
(  376)       !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
(  377)       !! History :
(  378)       !!   5.0  !  91-11  (G. Madec)  Original code
(  379)       !!   6.0  !  96-01  (G. Madec)  s-coord, suppress work arrays
(  380)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  381)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(  382)       !!----------------------------------------------------------------------
(  383)       !! * modules used
(  384)       USE oce, ONLY:   zwx  => ta,   & ! use ta as 3D workspace
(  385)                        zwy  => sa      ! use sa as 3D workspace
(  386)       !! * Arguments
(  387)       INTEGER, INTENT( in ) ::   kt    ! ocean timestep
(  388) 
(  389)       !! * Local declarations
(  390)       INTEGER ::   ji, jj, jk          ! dummy loop indices
(  391)       REAL(wp) ::   &
(  392)          zfact1, zua, zva, zuav, zvau  ! temporary scalars
(  393)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  394)          zcu, zcv, zwz,              & ! temporary workspace






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 8

(  395)          ztdua, ztdva                  ! temporary workspace
(  396)       !!----------------------------------------------------------------------
(  397)       
(  398)       IF( kt == nit000 ) THEN
(  399)          IF(lwp) WRITE(numout,*)
(  400)          IF(lwp) WRITE(numout,*) 'dyn_vor_enstrophy : vorticity term: enstrophy conserving scheme'
(  401)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~'
(  402)       ENDIF
(  403) 
(  404)       ! Local constant initialization
(  405)       zfact1 = 0.5 * 0.25
(  406) 
(  407)       ! Save ua and va trends
(  408)       IF( l_trddyn )   THEN
(  409)          ztdua(:,:,:) = ua(:,:,:) 
(  410)          ztdva(:,:,:) = va(:,:,:) 
(  411)          zcu(:,:,:) = 0.e0
(  412)          zcv(:,:,:) = 0.e0
(  413)       ENDIF
(  414) 
(  415)       !                                                ! ===============
(  416)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  417)          !                                             ! ===============
(  418) 
(  419)          ! Potential vorticity and horizontal fluxes
(  420)          ! -----------------------------------------
(  421)          IF( lk_sco ) THEN
(  422)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop 
(  423)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  424)                   zwz(ji,jj,jk) = ( rotn(ji,jj,jk) + ff(ji,jj) ) / fse3f(ji,jj,jk)
(  425)                   zwx(ji,jj,jk) = e2u(ji,jj) * fse3u(ji,jj,jk) * un(ji,jj,jk)
(  426)                   zwy(ji,jj,jk) = e1v(ji,jj) * fse3v(ji,jj,jk) * vn(ji,jj,jk)
(  427)                END DO
(  428)             END DO
(  429)          ELSE
(  430)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop 
(  431)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  432)                   zwz(ji,jj,jk) = rotn(ji,jj,jk) + ff(ji,jj)
(  433)                   zwx(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
(  434)                   zwy(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
(  435)                END DO
(  436)             END DO
(  437)          ENDIF
(  438) 
(  439) 
(  440)          ! Compute and add the vorticity term trend
(  441)          ! ----------------------------------------
(  442)          DO jj = 2, jpjm1
(  443)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  444)                zuav = zfact1 / e1u(ji,jj) * ( zwy(ji  ,jj-1,jk) + zwy(ji+1,jj-1,jk)   &
(  445)                                             + zwy(ji  ,jj  ,jk) + zwy(ji+1,jj  ,jk) )
(  446)                zvau =-zfact1 / e2v(ji,jj) * ( zwx(ji-1,jj  ,jk) + zwx(ji-1,jj+1,jk)   &
(  447)                                             + zwx(ji  ,jj  ,jk) + zwx(ji  ,jj+1,jk) )
(  448) 
(  449)                zua  = zuav * ( zwz(ji  ,jj-1,jk) + zwz(ji,jj,jk) )
(  450)                zva  = zvau * ( zwz(ji-1,jj  ,jk) + zwz(ji,jj,jk) )
(  451) 
(  452)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 9

(  453)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  454)             END DO  
(  455)          END DO  
(  456)          !                                             ! ===============
(  457)       END DO                                           !   End of slab
(  458)       !                                                ! ===============
(  459) 
(  460) 
(  461)       ! save the relative & planetary vorticity trends for diagnostic
(  462)       ! momentum trends
(  463)       IF( l_trddyn )   THEN
(  464)          ! Compute the planetary vorticity term trend
(  465)          !                                                ! ===============
(  466)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  467)             !                                             ! ===============
(  468)             DO jj = 2, jpjm1
(  469)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  470)                   zuav = zfact1 / e1u(ji,jj) * ( zwy(ji  ,jj-1,jk) + zwy(ji+1,jj-1,jk)   &
(  471)                      &                         + zwy(ji  ,jj  ,jk) + zwy(ji+1,jj  ,jk) )
(  472)                   zvau =-zfact1 / e2v(ji,jj) * ( zwx(ji-1,jj  ,jk) + zwx(ji-1,jj+1,jk)   &
(  473)                      &                         + zwx(ji  ,jj  ,jk) + zwx(ji  ,jj+1,jk) )
(  474) # if defined key_s_coord
(  475)                   zcu(ji,jj,jk)  = zuav * ( ff(ji  ,jj-1) / fse3f(ji  ,jj-1,jk)   &
(  476)                     &                     + ff(ji  ,jj  ) / fse3f(ji  ,jj  ,jk) )
(  477)                   zcv(ji,jj,jk)  = zvau * ( ff(ji-1,jj  ) / fse3f(ji-1,jj  ,jk)   &
(  478)                     &                     + ff(ji  ,jj  ) / fse3f(ji  ,jj  ,jk) )
(  479) # else
(  480)                   zcu(ji,jj,jk) = zuav * ( ff(ji  ,jj-1) + ff(ji,jj) )
(  481)                   zcv(ji,jj,jk) = zvau * ( ff(ji-1,jj  ) + ff(ji,jj) )
(  482) # endif
(  483)                END DO  
(  484)             END DO  
(  485)             !                                             ! ===============
(  486)          END DO                                           !   End of slab
(  487)          !                                                ! ===============
(  488) 
(  489)          ! Compute the relative vorticity term trend
(  490)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
(  491)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 
(  492) 
(  493)          CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
(  494)          CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
(  495)          CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
(  496)       ENDIF
(  497) 
(  498)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  499)          CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
(  500)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  501)       ENDIF
(  502) 
(  503)    END SUBROUTINE dyn_vor_enstrophy
(  504) 
(  505) 
(  506)    SUBROUTINE dyn_vor_ene_ens( kt )
(  507)       !!----------------------------------------------------------------------
(  508)       !!                ***  ROUTINE dyn_vor_ene_ens  ***
(  509)       !!
(  510)       !! ** Purpose :   Compute the now total vorticity trend and add it to 






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 10

(  511)       !!      the general trend of the momentum equation.
(  512)       !!
(  513)       !! ** Method  :   Trend evaluated using now fields (centered in time) 
(  514)       !!      and the Arakawa and Lamb (19XX) flux form formulation : conserves 
(  515)       !!      both the horizontal kinetic energy and the potential enstrophy
(  516)       !!      when horizontal divergence is zero.
(  517)       !!      The trend of the vorticity term is given by:
(  518)       !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
(  519)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(  520)       !!      Add this trend to the general momentum trend (ua,va):
(  521)       !!          (ua,va) = (ua,va) + ( voru , vorv )
(  522)       !!
(  523)       !! ** Action : - Update (ua,va) with the now vorticity term trend
(  524)       !!             - save the trends in (utrd,vtrd) in 2 parts (relative
(  525)       !!               and planetary vorticity trends) ('key_trddyn')
(  526)       !!
(  527)       !! References :
(  528)       !!      Arakawa and Lamb 1980, A potential enstrophy and energy conserving
(  529)       !!                             scheme for the Shallow water equations, 
(  530)       !!                             Monthly Weather Review, vol. 109, p 18-36
(  531)       !!
(  532)       !! History :
(  533)       !!   8.5  !  04-02  (G. Madec)  Original code
(  534)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(  535)       !!----------------------------------------------------------------------
(  536)       !! * Modules used     
(  537)       USE oce, ONLY :    ztdua => ta,    & ! use ta as 3D workspace   
(  538)                          ztdva => sa       ! use sa as 3D workspace   
(  539) 
(  540)       !! * Arguments
(  541)       INTEGER, INTENT( in ) ::   kt        ! ocean time-step index
(  542) 
(  543)       !! * Local declarations
(  544)       INTEGER ::   ji, jj, jk              ! dummy loop indices
(  545)       REAL(wp) ::   &
(  546)          zfac12, zua, zva                  ! temporary scalars
(  547)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  548)          zwx, zwy, zwz,                 &  ! temporary workspace
(  549)          ztnw, ztne, ztsw, ztse,        &  !    "           "
(  550)          zcor                              ! potential planetary vorticity (f/e3)
(  551)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  552)          zcu, zcv                          ! temporary workspace  
(  553)       REAL(wp), DIMENSION(jpi,jpj,jpk), SAVE ::   &
(  554)          ze3f
(  555)       !!----------------------------------------------------------------------
(  556) 
(  557)       IF( kt == nit000 ) THEN
(  558)          IF(lwp) WRITE(numout,*)
(  559)          IF(lwp) WRITE(numout,*) 'dyn_vor_ene_ens : vorticity term: energy and enstrophy conserving scheme'
(  560)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
(  561) 
(  562)          DO jk = 1, jpk
(  563)             DO jj = 1, jpjm1
(  564)                DO ji = 1, jpim1
(  565)                   ze3f(ji,jj,jk) = ( fse3t(ji,jj+1,jk)*tmask(ji,jj+1,jk) + fse3t(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)   &
(  566)                      &             + fse3t(ji,jj  ,jk)*tmask(ji,jj  ,jk) + fse3t(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) ) * 0.25
(  567) !!!               ze3f(ji,jj,jk) = MAX( ze3f(ji,jj,jk) , 1.e-20)
(  568)                   IF( ze3f(ji,jj,jk) /= 0.e0 )   ze3f(ji,jj,jk) = 1.e0 / ze3f(ji,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 11

(  569)                END DO
(  570)             END DO
(  571)          END DO
(  572)          CALL lbc_lnk( ze3f, 'F', 1. )
(  573)       ENDIF
(  574) 
(  575)       ! Local constant initialization
(  576)       zfac12 = 1.e0 / 12.e0
(  577) 
(  578)       ! Save ua and va trends
(  579)       IF( l_trddyn )   THEN
(  580)          ztdua(:,:,:) = ua(:,:,:) 
(  581)          ztdva(:,:,:) = va(:,:,:) 
(  582)          zcu(:,:,:) = 0.e0
(  583)          zcv(:,:,:) = 0.e0
(  584)       ENDIF
(  585)       
(  586)       !                                                ! ===============
(  587)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  588)          !                                             ! ===============
(  589)          
(  590)          ! Potential vorticity and horizontal fluxes
(  591)          ! -----------------------------------------
(  592) !!!bug   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) / fse3f(:,:,jk)
(  593)          zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) * ze3f(:,:,jk)
(  594)          zwx(:,:) = e2u(:,:) * fse3u(:,:,jk) * un(:,:,jk)
(  595)          zwy(:,:) = e1v(:,:) * fse3v(:,:,jk) * vn(:,:,jk)
(  596)          zcor(:,:) = ff(:,:) * ze3f(:,:,jk)
(  597) 
(  598)          ! Compute and add the vorticity term trend
(  599)          ! ----------------------------------------
(  600)          jj=2
(  601)          ztne(1,:) = 0 ; ztnw(1,:) = 0 ; ztse(1,:) = 0 ; ztsw(1,:) = 0
(  602)          DO ji = 2, jpi   
(  603)                ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
(  604)                ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
(  605)                ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
(  606)                ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
(  607)          END DO
(  608)          DO jj = 3, jpj
(  609)             DO ji = fs_2, jpi   ! vector opt.
(  610)                ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
(  611)                ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
(  612)                ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
(  613)                ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
(  614)             END DO
(  615)          END DO
(  616)          
(  617)          DO jj = 2, jpjm1
(  618)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  619)                zua = + zfac12 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
(  620)                   &                           + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
(  621)                zva = - zfac12 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
(  622)                   &                           + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
(  623)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  624)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  625)             END DO  
(  626)          END DO  






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 12

(  627)          !                                             ! ===============
(  628)       END DO                                           !   End of slab
(  629)       !                                                ! ===============
(  630) 
(  631)       ! save the relative & planetary vorticity trends for diagnostic
(  632)       ! momentum trends
(  633)       IF( l_trddyn )   THEN
(  634)          ! Compute the planetary vorticity term trend
(  635)          !                                                ! ===============
(  636)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  637)             !                                             ! ===============
(  638)             DO jj = 2, jpjm1
(  639)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  640)                   zcu(ji,jj,jk) = + zfac12 / e1u(ji,jj) * (  zcor(ji,jj  ) * zwy(ji  ,jj  ) + zcor(ji+1,jj) * zwy(ji+1,jj  )
(  641)                     &                                      + zcor(ji,jj  ) * zwy(ji  ,jj-1) + zcor(ji+1,jj) * zwy(ji+1,jj-1)
(  642)                   zcv(ji,jj,jk) = - zfac12 / e2v(ji,jj) * (  zcor(ji,jj+1) * zwx(ji-1,jj+1) + zcor(ji,jj+1) * zwx(ji  ,jj+1)
(  643)                     &                                      + zcor(ji,jj  ) * zwx(ji-1,jj  ) + zcor(ji,jj  ) * zwx(ji  ,jj  )
(  644)                END DO  
(  645)             END DO  
(  646)             !                                             ! ===============
(  647)          END DO                                           !   End of slab
(  648)          !                                                ! ===============
(  649) 
(  650)          ! Compute the relative vorticity term trend
(  651)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
(  652)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 
(  653) 
(  654)          CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
(  655)          CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
(  656)          CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
(  657)       ENDIF
(  658) 
(  659)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  660)          CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
(  661)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  662)       ENDIF
(  663) 
(  664)    END SUBROUTINE dyn_vor_ene_ens
(  665) 
(  666) 
(  667)    SUBROUTINE dyn_vor_ctl
(  668)       !!---------------------------------------------------------------------
(  669)       !!                  ***  ROUTINE dyn_vor_ctl  ***
(  670)       !!
(  671)       !! ** Purpose :   Control the consistency between cpp options for
(  672)       !!      tracer advection schemes
(  673)       !!
(  674)       !! History :
(  675)       !!   9.0  !  03-08  (G. Madec)  Original code
(  676)       !!----------------------------------------------------------------------
(  677)       !! * Local declarations
(  678)       INTEGER ::   ioptio          ! temporary integer
(  679) 
(  680)       NAMELIST/nam_dynvor/ ln_dynvor_ens, ln_dynvor_ene, ln_dynvor_mix, ln_dynvor_een
(  681)       !!----------------------------------------------------------------------
(  682) 
(  683)       ! Read Namelist nam_dynvor : Vorticity scheme options
(  684)       ! ------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 13

(  685)       REWIND ( numnam )
(  686)       READ   ( numnam, nam_dynvor )
(  687) 
(  688)       ! Control of vorticity scheme options
(  689)       ! -----------------------------------
(  690)       ! Control print
(  691)       IF(lwp) THEN
(  692)          WRITE(numout,*)
(  693)          WRITE(numout,*) 'dyn_vor_ctl : vorticity term : read namelist and control the consistency'
(  694)          WRITE(numout,*) '~~~~~~~~~~~'
(  695)          WRITE(numout,*) '              Namelist nam_dynvor : oice of the vorticity term scheme'
(  696)          WRITE(numout,*) '                 enstrophy conserving scheme                ln_dynvor_ens = ', ln_dynvor_ens
(  697)          WRITE(numout,*) '                 energy    conserving scheme                ln_dynvor_ene = ', ln_dynvor_ene
(  698)          WRITE(numout,*) '                 mixed enstrophy/energy conserving scheme   ln_dynvor_mix = ', ln_dynvor_mix
(  699)          WRITE(numout,*) '                 enstrophy and energy conserving scheme     ln_dynvor_een = ', ln_dynvor_een
(  700)       ENDIF
(  701) 
(  702)       ioptio = 0
(  703) 
(  704)       IF( ln_dynvor_ens ) THEN
(  705)          IF(lwp) WRITE(numout,*)
(  706)          IF(lwp) WRITE(numout,*) '              vorticity term : enstrophy conserving scheme'
(  707)          ioptio = ioptio + 1
(  708)       ENDIF
(  709)       IF( ln_dynvor_ene ) THEN
(  710)          IF(lwp) WRITE(numout,*)
(  711)          IF(lwp) WRITE(numout,*) '              vorticity term : energy conserving scheme'
(  712)          ioptio = ioptio + 1
(  713)       ENDIF
(  714)       IF( ln_dynvor_mix ) THEN
(  715)          IF(lwp) WRITE(numout,*)
(  716)          IF(lwp) WRITE(numout,*) '              vorticity term : mixed enstrophy/energy conserving scheme'
(  717)          ioptio = ioptio + 1
(  718)       ENDIF
(  719)       IF( ln_dynvor_een ) THEN
(  720)          IF(lwp) WRITE(numout,*)
(  721)          IF(lwp) WRITE(numout,*) '              vorticity term : energy and enstrophy conserving scheme'
(  722)          ioptio = ioptio + 1
(  723)       ENDIF
(  724)       IF ( ioptio /= 1 .AND. .NOT. lk_esopa ) THEN
(  725)           if(lwp) WRITE(numout,cform_err)
(  726)           IF(lwp) WRITE(numout,*) ' use ONE and ONLY one vorticity scheme'
(  727)           nstop = nstop + 1
(  728)       ENDIF
(  729) 
(  730)    END SUBROUTINE dyn_vor_ctl
(  731) 
(  732) !!==============================================================================
(  733) END MODULE dynvor















PGF90 (Version     10.5)          09/23/2010  09:49:40      page 14

(    1) # 1 "dynvor.F90"
(    1) MODULE dynvor
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynvor  ***
(    4)    !! Ocean dynamics: Update the momentum trend with the relative and
(    5)    !!                 planetary vorticity trends
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   dyn_vor_enstrophy: enstrophy conserving scheme       (ln_dynvor_ens=T)
(   10)    !!   dyn_vor_energy   : energy conserving scheme          (ln_dynvor_ene=T)
(   11)    !!   dyn_vor_mixed    : mixed enstrophy/energy conserving (ln_dynvor_mix=T)
(   12)    !!   dyn_vor_ene_ens  : energy and enstrophy conserving   (ln_dynvor_een=T)
(   13)    !!   dyn_vor_ctl      : control of the different vorticity option
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE trdmod          ! ocean dynamics trends 
(   20)    USE trdmod_oce      ! ocean variables trends
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22)    USE prtctl          ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC dyn_vor_enstrophy      ! routine called by step.F90
(   29)    PUBLIC dyn_vor_energy         ! routine called by step.F90
(   30)    PUBLIC dyn_vor_mixed          ! routine called by step.F90
(   31)    PUBLIC dyn_vor_ene_ens        ! routine called by step.F90
(   32)    PUBLIC dyn_vor_ctl            ! routine called by step.F90
(   33) 
(   34)    !! * Shared module variables
(   35)    LOGICAL, PUBLIC ::   ln_dynvor_ene = .FALSE.   !: energy conserving scheme
(   36)    LOGICAL, PUBLIC ::   ln_dynvor_ens = .TRUE.    !: enstrophy conserving scheme
(   37)    LOGICAL, PUBLIC ::   ln_dynvor_mix = .FALSE.   !: mixed scheme
(   38)    LOGICAL, PUBLIC ::   ln_dynvor_een = .FALSE.    !: energy and enstrophy conserving scheme
(   39) 
(   40)    !! * Substitutions
(   41) # 1 "./domzgr_substitute.h90"
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
(   48) # 42 "dynvor.F90"
(   42) # 1 "./vectopt_loop_substitute.h90"






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 15

(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 43 "dynvor.F90"
(   43) # 43
(   43)    !!----------------------------------------------------------------------
(   44)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   45)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynvor.F90,v 1.11 2006/03/21 08:26:27 opalod Exp $ 
(   46)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   47)    !!----------------------------------------------------------------------
(   48) 
(   49) CONTAINS
(   50) 
(   51)    SUBROUTINE dyn_vor_energy( kt )
(   52)       !!----------------------------------------------------------------------
(   53)       !!                  ***  ROUTINE dyn_vor_energy  ***
(   54)       !!
(   55)       !! ** Purpose :   Compute the now total vorticity trend and add it to 
(   56)       !!      the general trend of the momentum equation.
(   57)       !!
(   58)       !! ** Method  :   Trend evaluated using now fields (centered in time) 
(   59)       !!      and the Sadourny (1975) flux form formulation : conserves the
(   60)       !!      horizontal kinetic energy.
(   61)       !!      The trend of the vorticity term is given by:
(   62)       !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
(   63)       !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f  mi(e1v*e3v vn) ]
(   64)       !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f  mj(e2u*e3u un) ]
(   65)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(   66)       !!          voru = 1/e1u  mj-1[ (rotn+f)  mi(e1v vn) ]
(   67)       !!          vorv = 1/e2v  mi-1[ (rotn+f)  mj(e2u un) ]
(   68)       !!      Add this trend to the general momentum trend (ua,va):
(   69)       !!          (ua,va) = (ua,va) + ( voru , vorv )
(   70)       !!
(   71)       !! ** Action : - Update (ua,va) with the now vorticity term trend
(   72)       !!             - save the trends in (utrd,vtrd) in 2 parts (relative
(   73)       !!               and planetary vorticity trends) ('key_trddyn')
(   74)       !!
(   75)       !! References :
(   76)       !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
(   77)       !! History :
(   78)       !!   5.0  !  91-11  (G. Madec)  Original code
(   79)       !!   6.0  !  96-01  (G. Madec)  s-coord, suppress work arrays
(   80)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   81)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   82)       !!----------------------------------------------------------------------
(   83)       !! * Modules used     
(   84)       USE oce, ONLY :    ztdua => ta,     & ! use ta as 3D workspace   
(   85)                          ztdva => sa        ! use sa as 3D workspace   
(   86) 
(   87)       !! * Arguments






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 16

(   88)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   89) 
(   90)       !! * Local declarations
(   91)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(   92)       REAL(wp) ::   &
(   93)          zfact2, zua, zva,               &  ! temporary scalars
(   94)          zx1, zx2, zy1, zy2                 !    "         "
(   95)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   96)          zwx, zwy, zwz                      ! temporary workspace
(   97)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   98)          zcu, zcv                           !    "         "
(   99)       !!----------------------------------------------------------------------
(  100) 
(  101)       IF( kt == nit000 ) THEN
(  102)          IF(lwp) WRITE(numout,*)
(  103)          IF(lwp) WRITE(numout,*) 'dyn_vor_energy : vorticity term: energy conserving scheme'
(  104)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
(  105)       ENDIF
(  106) 
(  107)       ! Local constant initialization
(  108)       zfact2 = 0.5 * 0.5
(  109)       zcu (:,:,:) = 0.e0
(  110)       zcv (:,:,:) = 0.e0
(  111) 
(  112)       ! Save ua and va trends
(  113)       IF( l_trddyn )   THEN
(  114)          ztdua(:,:,:) = ua(:,:,:) 
(  115)          ztdva(:,:,:) = va(:,:,:) 
(  116)          zcu(:,:,:) = 0.e0
(  117)          zcv(:,:,:) = 0.e0
(  118)       ENDIF
(  119)       
(  120)       !                                                ! ===============
(  121)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  122)          !                                             ! ===============
(  123)          
(  124)          ! Potential vorticity and horizontal fluxes
(  125)          ! -----------------------------------------
(  126)          IF( lk_sco ) THEN
(  127)             zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) / e3t(jk)
(  128)             zwx(:,:) = e2u(:,:) * e3t(jk) * un(:,:,jk)
(  129)             zwy(:,:) = e1v(:,:) * e3t(jk) * vn(:,:,jk)
(  130)          ELSE
(  131)             zwz(:,:) = rotn(:,:,jk) + ff(:,:)
(  132)             zwx(:,:) = e2u(:,:) * un(:,:,jk)
(  133)             zwy(:,:) = e1v(:,:) * vn(:,:,jk)
(  134)          ENDIF
(  135) 
(  136)          ! Compute and add the vorticity term trend
(  137)          ! ----------------------------------------
(  138)          DO jj = 2, jpjm1
(  139)             DO ji = 2, jpim1   ! vector opt.
(  140)                zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
(  141)                zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
(  142)                zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
(  143)                zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
(  144)                zua = zfact2 / e1u(ji,jj) * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
(  145)                zva =-zfact2 / e2v(ji,jj) * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 17

(  146)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  147)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  148)             END DO  
(  149)          END DO  
(  150)          !                                             ! ===============
(  151)       END DO                                           !   End of slab
(  152)       !                                                ! ===============
(  153) 
(  154)       ! save the relative & planetary vorticity trends for diagnostic
(  155)       ! momentum trends
(  156)       IF( l_trddyn )   THEN
(  157)          ! Compute the planetary vorticity term trend
(  158)          !                                                ! ===============
(  159)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  160)             !                                             ! ===============
(  161)             DO jj = 2, jpjm1
(  162)                DO ji = 2, jpim1   ! vector opt.
(  163)                   zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
(  164)                   zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
(  165)                   zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
(  166)                   zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
(  167) # 173
(  173)                  zcu(ji,jj,jk) = zfact2 / e1u(ji,jj) * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
(  174)                  zcv(ji,jj,jk) =-zfact2 / e2v(ji,jj) * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
(  175) # 176
(  176)                END DO  
(  177)             END DO  
(  178)             !                                             ! ===============
(  179)          END DO                                           !   End of slab
(  180)          !                                                ! ===============
(  181) 
(  182)          ! Compute the relative vorticity term trend
(  183)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
(  184)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 
(  185) 
(  186)          CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
(  187)          CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
(  188)          CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
(  189)       ENDIF
(  190) 
(  191)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  192)          CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
(  193)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  194)       ENDIF
(  195) 
(  196)    END SUBROUTINE dyn_vor_energy
(  197) 
(  198) 
(  199)    SUBROUTINE dyn_vor_mixed( kt )
(  200)       !!----------------------------------------------------------------------
(  201)       !!                 ***  ROUTINE dyn_vor_mixed  ***
(  202)       !!
(  203)       !! ** Purpose :   Compute the now total vorticity trend and add it to
(  204)       !!      the general trend of the momentum equation.
(  205)       !!
(  206)       !! ** Method  :   Trend evaluated using now fields (centered in time)
(  207)       !!      Mixte formulation : conserves the potential enstrophy of a hori-
(  208)       !!      zontally non-divergent flow for (rotzu x uh), the relative vor-






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 18

(  209)       !!      ticity term and the horizontal kinetic energy for (f x uh), the
(  210)       !!      coriolis term. the now trend of the vorticity term is given by:
(  211)       !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
(  212)       !!          voru = 1/e1u  mj-1(rotn/e3f) mj-1[ mi(e1v*e3v vn) ]
(  213)       !!              +1/e1u  mj-1[ f/e3f          mi(e1v*e3v vn) ]
(  214)       !!          vorv = 1/e2v  mi-1(rotn/e3f) mi-1[ mj(e2u*e3u un) ]
(  215)       !!              +1/e2v  mi-1[ f/e3f          mj(e2u*e3u un) ]
(  216)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(  217)       !!          voru = 1/e1u  mj-1(rotn) mj-1[ mi(e1v vn) ]
(  218)       !!              +1/e1u  mj-1[ f          mi(e1v vn) ]
(  219)       !!          vorv = 1/e2v  mi-1(rotn) mi-1[ mj(e2u un) ]
(  220)       !!              +1/e2v  mi-1[ f          mj(e2u un) ]
(  221)       !!      Add this now trend to the general momentum trend (ua,va):
(  222)       !!          (ua,va) = (ua,va) + ( voru , vorv )
(  223)       !!
(  224)       !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
(  225)       !!             - Save the trends in (utrd,vtrd) in 2 parts (relative
(  226)       !!               and planetary vorticity trends) ('key_trddyn')
(  227)       !!
(  228)       !! References :
(  229)       !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
(  230)       !! History :
(  231)       !!   5.0  !  91-11  (G. Madec) Original code, enstrophy-energy-combined schemes
(  232)       !!   6.0  !  96-01  (G. Madec) s-coord, suppress work arrays
(  233)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  234)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(  235)       !!----------------------------------------------------------------------
(  236)       !! * Modules used     
(  237)       USE oce, ONLY :    ztdua => ta,     & ! use ta as 3D workspace   
(  238)                          ztdva => sa        ! use sa as 3D workspace   
(  239)       !! * Arguments
(  240)       INTEGER, INTENT( in ) ::   kt         ! ocean timestep index
(  241) 
(  242)       !! * Local declarations
(  243)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(  244)       REAL(wp) ::   &
(  245)          zfact1, zfact2, zua, zva,       &  ! temporary scalars
(  246)          zcua, zcva, zx1, zx2, zy1, zy2
(  247)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  248)          zwx, zwy, zwz, zww                 ! temporary workspace
(  249)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  250)          zcu, zcv                           !    "         "
(  251)       !!----------------------------------------------------------------------
(  252) 
(  253)       IF( kt == nit000 ) THEN
(  254)          IF(lwp) WRITE(numout,*)
(  255)          IF(lwp) WRITE(numout,*) 'dyn_vor_mixed : vorticity term: mixed energy/enstrophy conserving scheme'
(  256)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
(  257)       ENDIF
(  258) 
(  259)       ! Local constant initialization
(  260)       zfact1 = 0.5 * 0.25
(  261)       zfact2 = 0.5 * 0.5
(  262) 
(  263)       ! Save ua and va trends
(  264)       IF( l_trddyn )   THEN
(  265)          ztdua(:,:,:) = ua(:,:,:) 
(  266)          ztdva(:,:,:) = va(:,:,:) 






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 19

(  267)          zcu(:,:,:) = 0.e0
(  268)          zcv(:,:,:) = 0.e0
(  269)       ENDIF
(  270) 
(  271)       !                                                ! ===============
(  272)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  273)          !                                             ! ===============
(  274) 
(  275)          ! Relative and planetary potential vorticity and horizontal fluxes
(  276)          ! ----------------------------------------------------------------
(  277)          IF( lk_sco ) THEN        
(  278)             zwz(:,:) = ff  (:,:)    / e3t(jk)
(  279)             zww(:,:) = rotn(:,:,jk) / e3t(jk)
(  280)             zwx(:,:) = e2u(:,:) * e3t(jk) * un(:,:,jk)
(  281)             zwy(:,:) = e1v(:,:) * e3t(jk) * vn(:,:,jk)
(  282)          ELSE
(  283)             zwz(:,:) = ff(:,:)
(  284)             zww(:,:) = rotn(:,:,jk)
(  285)             zwx(:,:) = e2u(:,:) * un(:,:,jk)
(  286)             zwy(:,:) = e1v(:,:) * vn(:,:,jk)
(  287)          ENDIF
(  288) 
(  289)          ! Compute and add the vorticity term trend
(  290)          ! ----------------------------------------
(  291)          DO jj = 2, jpjm1
(  292)             DO ji = 2, jpim1   ! vector opt.
(  293)                zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  294)                zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  295)                zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  296)                zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  297)                ! enstrophy conserving formulation for relative vorticity term
(  298)                zua = zfact1 * ( zww(ji  ,jj-1) + zww(ji,jj) ) * ( zy1 + zy2 )
(  299)                zva =-zfact1 * ( zww(ji-1,jj  ) + zww(ji,jj) ) * ( zx1 + zx2 )
(  300)                ! energy conserving formulation for planetary vorticity term
(  301)                zcua = zfact2 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
(  302)                zcva =-zfact2 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
(  303) 
(  304)                ua(ji,jj,jk) = ua(ji,jj,jk) + zcua + zua
(  305)                va(ji,jj,jk) = va(ji,jj,jk) + zcva + zva
(  306)             END DO  
(  307)          END DO  
(  308)          !                                             ! ===============
(  309)       END DO                                           !   End of slab
(  310)       !                                                ! ===============
(  311) 
(  312)       ! save the relative & planetary vorticity trends for diagnostic
(  313)       ! momentum trends
(  314)       IF( l_trddyn )   THEN
(  315)          ! Compute the planetary vorticity term trend
(  316)          !                                                ! ===============
(  317)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  318)             !                                             ! ===============
(  319)             DO jj = 2, jpjm1
(  320)                DO ji = 2, jpim1   ! vector opt.
(  321)                   zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  322)                   zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  323)                   zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  324)                   zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 20

(  325)                   ! energy conserving formulation for planetary vorticity term
(  326)                   zcu(ji,jj,jk) = zfact2 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
(  327)                   zcv(ji,jj,jk) =-zfact2 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
(  328)                END DO  
(  329)             END DO  
(  330)             !                                             ! ===============
(  331)          END DO                                           !   End of slab
(  332)          !                                                ! ===============
(  333) 
(  334)          ! Compute the relative vorticity term trend
(  335)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
(  336)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 
(  337) 
(  338)          CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
(  339)          CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
(  340)          CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
(  341)       ENDIF
(  342) 
(  343)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  344)          CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
(  345)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  346)       ENDIF
(  347) 
(  348)    END SUBROUTINE dyn_vor_mixed
(  349) 
(  350) 
(  351)    SUBROUTINE dyn_vor_enstrophy( kt )
(  352)       !!----------------------------------------------------------------------
(  353)       !!                ***  ROUTINE dyn_vor_enstrophy  ***
(  354)       !!
(  355)       !! ** Purpose :   Compute the now total vorticity trend and add it to
(  356)       !!      the general trend of the momentum equation.
(  357)       !!
(  358)       !! ** Method  :   Trend evaluated using now fields (centered in time)
(  359)       !!      and the Sadourny (1975) flux FORM formulation : conserves the
(  360)       !!      potential enstrophy of a horizontally non-divergent flow. the
(  361)       !!      trend of the vorticity term is given by:
(  362)       !!       * s-coordinate (lk_sco=T), the e3. are inside the derivative:
(  363)       !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f ]  mj-1[ mi(e1v*e3v vn) ]
(  364)       !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f ]  mi-1[ mj(e2u*e3u un) ]
(  365)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(  366)       !!          voru = 1/e1u  mj-1[ rotn+f ]  mj-1[ mi(e1v vn) ]
(  367)       !!          vorv = 1/e2v  mi-1[ rotn+f ]  mi-1[ mj(e2u un) ]
(  368)       !!      Add this trend to the general momentum trend (ua,va):
(  369)       !!          (ua,va) = (ua,va) + ( voru , vorv )
(  370)       !!
(  371)       !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
(  372)       !!             - Save the trends in (utrd,vtrd) in 2 parts (relative 
(  373)       !!               and planetary vorticity trends) ('key_trddyn')
(  374)       !!
(  375)       !! References :
(  376)       !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
(  377)       !! History :
(  378)       !!   5.0  !  91-11  (G. Madec)  Original code
(  379)       !!   6.0  !  96-01  (G. Madec)  s-coord, suppress work arrays
(  380)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  381)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(  382)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 21

(  383)       !! * modules used
(  384)       USE oce, ONLY:   zwx  => ta,   & ! use ta as 3D workspace
(  385)                        zwy  => sa      ! use sa as 3D workspace
(  386)       !! * Arguments
(  387)       INTEGER, INTENT( in ) ::   kt    ! ocean timestep
(  388) 
(  389)       !! * Local declarations
(  390)       INTEGER ::   ji, jj, jk          ! dummy loop indices
(  391)       REAL(wp) ::   &
(  392)          zfact1, zua, zva, zuav, zvau  ! temporary scalars
(  393)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  394)          zcu, zcv, zwz,              & ! temporary workspace
(  395)          ztdua, ztdva                  ! temporary workspace
(  396)       !!----------------------------------------------------------------------
(  397)       
(  398)       IF( kt == nit000 ) THEN
(  399)          IF(lwp) WRITE(numout,*)
(  400)          IF(lwp) WRITE(numout,*) 'dyn_vor_enstrophy : vorticity term: enstrophy conserving scheme'
(  401)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~'
(  402)       ENDIF
(  403) 
(  404)       ! Local constant initialization
(  405)       zfact1 = 0.5 * 0.25
(  406) 
(  407)       ! Save ua and va trends
(  408)       IF( l_trddyn )   THEN
(  409)          ztdua(:,:,:) = ua(:,:,:) 
(  410)          ztdva(:,:,:) = va(:,:,:) 
(  411)          zcu(:,:,:) = 0.e0
(  412)          zcv(:,:,:) = 0.e0
(  413)       ENDIF
(  414) 
(  415)       !                                                ! ===============
(  416)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  417)          !                                             ! ===============
(  418) 
(  419)          ! Potential vorticity and horizontal fluxes
(  420)          ! -----------------------------------------
(  421)          IF( lk_sco ) THEN
(  422)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop 
(  423)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  424)                   zwz(ji,jj,jk) = ( rotn(ji,jj,jk) + ff(ji,jj) ) / e3t(jk)
(  425)                   zwx(ji,jj,jk) = e2u(ji,jj) * e3t(jk) * un(ji,jj,jk)
(  426)                   zwy(ji,jj,jk) = e1v(ji,jj) * e3t(jk) * vn(ji,jj,jk)
(  427)                END DO
(  428)             END DO
(  429)          ELSE
(  430)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop 
(  431)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  432)                   zwz(ji,jj,jk) = rotn(ji,jj,jk) + ff(ji,jj)
(  433)                   zwx(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
(  434)                   zwy(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
(  435)                END DO
(  436)             END DO
(  437)          ENDIF
(  438) 
(  439) 
(  440)          ! Compute and add the vorticity term trend






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 22

(  441)          ! ----------------------------------------
(  442)          DO jj = 2, jpjm1
(  443)             DO ji = 2, jpim1   ! vector opt.
(  444)                zuav = zfact1 / e1u(ji,jj) * ( zwy(ji  ,jj-1,jk) + zwy(ji+1,jj-1,jk)   &
(  445)                                             + zwy(ji  ,jj  ,jk) + zwy(ji+1,jj  ,jk) )
(  446)                zvau =-zfact1 / e2v(ji,jj) * ( zwx(ji-1,jj  ,jk) + zwx(ji-1,jj+1,jk)   &
(  447)                                             + zwx(ji  ,jj  ,jk) + zwx(ji  ,jj+1,jk) )
(  448) 
(  449)                zua  = zuav * ( zwz(ji  ,jj-1,jk) + zwz(ji,jj,jk) )
(  450)                zva  = zvau * ( zwz(ji-1,jj  ,jk) + zwz(ji,jj,jk) )
(  451) 
(  452)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  453)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  454)             END DO  
(  455)          END DO  
(  456)          !                                             ! ===============
(  457)       END DO                                           !   End of slab
(  458)       !                                                ! ===============
(  459) 
(  460) 
(  461)       ! save the relative & planetary vorticity trends for diagnostic
(  462)       ! momentum trends
(  463)       IF( l_trddyn )   THEN
(  464)          ! Compute the planetary vorticity term trend
(  465)          !                                                ! ===============
(  466)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  467)             !                                             ! ===============
(  468)             DO jj = 2, jpjm1
(  469)                DO ji = 2, jpim1   ! vector opt.
(  470)                   zuav = zfact1 / e1u(ji,jj) * ( zwy(ji  ,jj-1,jk) + zwy(ji+1,jj-1,jk)   &
(  471)                      &                         + zwy(ji  ,jj  ,jk) + zwy(ji+1,jj  ,jk) )
(  472)                   zvau =-zfact1 / e2v(ji,jj) * ( zwx(ji-1,jj  ,jk) + zwx(ji-1,jj+1,jk)   &
(  473)                      &                         + zwx(ji  ,jj  ,jk) + zwx(ji  ,jj+1,jk) )
(  474) # 480
(  480)                   zcu(ji,jj,jk) = zuav * ( ff(ji  ,jj-1) + ff(ji,jj) )
(  481)                   zcv(ji,jj,jk) = zvau * ( ff(ji-1,jj  ) + ff(ji,jj) )
(  482) # 483
(  483)                END DO  
(  484)             END DO  
(  485)             !                                             ! ===============
(  486)          END DO                                           !   End of slab
(  487)          !                                                ! ===============
(  488) 
(  489)          ! Compute the relative vorticity term trend
(  490)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
(  491)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 
(  492) 
(  493)          CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
(  494)          CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
(  495)          CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
(  496)       ENDIF
(  497) 
(  498)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  499)          CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
(  500)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  501)       ENDIF
(  502) 
(  503)    END SUBROUTINE dyn_vor_enstrophy






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 23

(  504) 
(  505) 
(  506)    SUBROUTINE dyn_vor_ene_ens( kt )
(  507)       !!----------------------------------------------------------------------
(  508)       !!                ***  ROUTINE dyn_vor_ene_ens  ***
(  509)       !!
(  510)       !! ** Purpose :   Compute the now total vorticity trend and add it to 
(  511)       !!      the general trend of the momentum equation.
(  512)       !!
(  513)       !! ** Method  :   Trend evaluated using now fields (centered in time) 
(  514)       !!      and the Arakawa and Lamb (19XX) flux form formulation : conserves 
(  515)       !!      both the horizontal kinetic energy and the potential enstrophy
(  516)       !!      when horizontal divergence is zero.
(  517)       !!      The trend of the vorticity term is given by:
(  518)       !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
(  519)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(  520)       !!      Add this trend to the general momentum trend (ua,va):
(  521)       !!          (ua,va) = (ua,va) + ( voru , vorv )
(  522)       !!
(  523)       !! ** Action : - Update (ua,va) with the now vorticity term trend
(  524)       !!             - save the trends in (utrd,vtrd) in 2 parts (relative
(  525)       !!               and planetary vorticity trends) ('key_trddyn')
(  526)       !!
(  527)       !! References :
(  528)       !!      Arakawa and Lamb 1980, A potential enstrophy and energy conserving
(  529)       !!                             scheme for the Shallow water equations, 
(  530)       !!                             Monthly Weather Review, vol. 109, p 18-36
(  531)       !!
(  532)       !! History :
(  533)       !!   8.5  !  04-02  (G. Madec)  Original code
(  534)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(  535)       !!----------------------------------------------------------------------
(  536)       !! * Modules used     
(  537)       USE oce, ONLY :    ztdua => ta,    & ! use ta as 3D workspace   
(  538)                          ztdva => sa       ! use sa as 3D workspace   
(  539) 
(  540)       !! * Arguments
(  541)       INTEGER, INTENT( in ) ::   kt        ! ocean time-step index
(  542) 
(  543)       !! * Local declarations
(  544)       INTEGER ::   ji, jj, jk              ! dummy loop indices
(  545)       REAL(wp) ::   &
(  546)          zfac12, zua, zva                  ! temporary scalars
(  547)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  548)          zwx, zwy, zwz,                 &  ! temporary workspace
(  549)          ztnw, ztne, ztsw, ztse,        &  !    "           "
(  550)          zcor                              ! potential planetary vorticity (f/e3)
(  551)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  552)          zcu, zcv                          ! temporary workspace  
(  553)       REAL(wp), DIMENSION(jpi,jpj,jpk), SAVE ::   &
(  554)          ze3f
(  555)       !!----------------------------------------------------------------------
(  556) 
(  557)       IF( kt == nit000 ) THEN
(  558)          IF(lwp) WRITE(numout,*)
(  559)          IF(lwp) WRITE(numout,*) 'dyn_vor_ene_ens : vorticity term: energy and enstrophy conserving scheme'
(  560)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
(  561) 






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 24

(  562)          DO jk = 1, jpk
(  563)             DO jj = 1, jpjm1
(  564)                DO ji = 1, jpim1
(  565)                   ze3f(ji,jj,jk) = ( e3t(jk)*tmask(ji,jj+1,jk) + e3t(jk)*tmask(ji+1,jj+1,jk)   &
(  566)                      &             + e3t(jk)*tmask(ji,jj  ,jk) + e3t(jk)*tmask(ji+1,jj  ,jk) ) * 0.25
(  567) !!!               ze3f(ji,jj,jk) = MAX( ze3f(ji,jj,jk) , 1.e-20)
(  568)                   IF( ze3f(ji,jj,jk) /= 0.e0 )   ze3f(ji,jj,jk) = 1.e0 / ze3f(ji,jj,jk)
(  569)                END DO
(  570)             END DO
(  571)          END DO
(  572)          CALL lbc_lnk( ze3f, 'F', 1. )
(  573)       ENDIF
(  574) 
(  575)       ! Local constant initialization
(  576)       zfac12 = 1.e0 / 12.e0
(  577) 
(  578)       ! Save ua and va trends
(  579)       IF( l_trddyn )   THEN
(  580)          ztdua(:,:,:) = ua(:,:,:) 
(  581)          ztdva(:,:,:) = va(:,:,:) 
(  582)          zcu(:,:,:) = 0.e0
(  583)          zcv(:,:,:) = 0.e0
(  584)       ENDIF
(  585)       
(  586)       !                                                ! ===============
(  587)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  588)          !                                             ! ===============
(  589)          
(  590)          ! Potential vorticity and horizontal fluxes
(  591)          ! -----------------------------------------
(  592) !!!bug   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) / e3t(jk)
(  593)          zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) * ze3f(:,:,jk)
(  594)          zwx(:,:) = e2u(:,:) * e3t(jk) * un(:,:,jk)
(  595)          zwy(:,:) = e1v(:,:) * e3t(jk) * vn(:,:,jk)
(  596)          zcor(:,:) = ff(:,:) * ze3f(:,:,jk)
(  597) 
(  598)          ! Compute and add the vorticity term trend
(  599)          ! ----------------------------------------
(  600)          jj=2
(  601)          ztne(1,:) = 0 ; ztnw(1,:) = 0 ; ztse(1,:) = 0 ; ztsw(1,:) = 0
(  602)          DO ji = 2, jpi   
(  603)                ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
(  604)                ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
(  605)                ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
(  606)                ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
(  607)          END DO
(  608)          DO jj = 3, jpj
(  609)             DO ji = 2, jpi   ! vector opt.
(  610)                ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
(  611)                ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
(  612)                ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
(  613)                ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
(  614)             END DO
(  615)          END DO
(  616)          
(  617)          DO jj = 2, jpjm1
(  618)             DO ji = 2, jpim1   ! vector opt.
(  619)                zua = + zfac12 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 25

(  620)                   &                           + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
(  621)                zva = - zfac12 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
(  622)                   &                           + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
(  623)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  624)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  625)             END DO  
(  626)          END DO  
(  627)          !                                             ! ===============
(  628)       END DO                                           !   End of slab
(  629)       !                                                ! ===============
(  630) 
(  631)       ! save the relative & planetary vorticity trends for diagnostic
(  632)       ! momentum trends
(  633)       IF( l_trddyn )   THEN
(  634)          ! Compute the planetary vorticity term trend
(  635)          !                                                ! ===============
(  636)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  637)             !                                             ! ===============
(  638)             DO jj = 2, jpjm1
(  639)                DO ji = 2, jpim1   ! vector opt.
(  640)                   zcu(ji,jj,jk) = + zfac12 / e1u(ji,jj) * (  zcor(ji,jj  ) * zwy(ji  ,jj  ) + zcor(ji+1,jj) * zwy(ji+1,jj  )   &
(  641)                     &                                      + zcor(ji,jj  ) * zwy(ji  ,jj-1) + zcor(ji+1,jj) * zwy(ji+1,jj-1) )
(  642)                   zcv(ji,jj,jk) = - zfac12 / e2v(ji,jj) * (  zcor(ji,jj+1) * zwx(ji-1,jj+1) + zcor(ji,jj+1) * zwx(ji  ,jj+1)   &
(  643)                     &                                      + zcor(ji,jj  ) * zwx(ji-1,jj  ) + zcor(ji,jj  ) * zwx(ji  ,jj  ) )
(  644)                END DO  
(  645)             END DO  
(  646)             !                                             ! ===============
(  647)          END DO                                           !   End of slab
(  648)          !                                                ! ===============
(  649) 
(  650)          ! Compute the relative vorticity term trend
(  651)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
(  652)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 
(  653) 
(  654)          CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
(  655)          CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
(  656)          CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
(  657)       ENDIF
(  658) 
(  659)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  660)          CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
(  661)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  662)       ENDIF
(  663) 
(  664)    END SUBROUTINE dyn_vor_ene_ens
(  665) 
(  666) 
(  667)    SUBROUTINE dyn_vor_ctl
(  668)       !!---------------------------------------------------------------------
(  669)       !!                  ***  ROUTINE dyn_vor_ctl  ***
(  670)       !!
(  671)       !! ** Purpose :   Control the consistency between cpp options for
(  672)       !!      tracer advection schemes
(  673)       !!
(  674)       !! History :
(  675)       !!   9.0  !  03-08  (G. Madec)  Original code
(  676)       !!----------------------------------------------------------------------
(  677)       !! * Local declarations






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 26

(  678)       INTEGER ::   ioptio          ! temporary integer
(  679) 
(  680)       NAMELIST/nam_dynvor/ ln_dynvor_ens, ln_dynvor_ene, ln_dynvor_mix, ln_dynvor_een
(  681)       !!----------------------------------------------------------------------
(  682) 
(  683)       ! Read Namelist nam_dynvor : Vorticity scheme options
(  684)       ! ------------------------
(  685)       REWIND ( numnam )
(  686)       READ   ( numnam, nam_dynvor )
(  687) 
(  688)       ! Control of vorticity scheme options
(  689)       ! -----------------------------------
(  690)       ! Control print
(  691)       IF(lwp) THEN
(  692)          WRITE(numout,*)
(  693)          WRITE(numout,*) 'dyn_vor_ctl : vorticity term : read namelist and control the consistency'
(  694)          WRITE(numout,*) '~~~~~~~~~~~'
(  695)          WRITE(numout,*) '              Namelist nam_dynvor : oice of the vorticity term scheme'
(  696)          WRITE(numout,*) '                 enstrophy conserving scheme                ln_dynvor_ens = ', ln_dynvor_ens
(  697)          WRITE(numout,*) '                 energy    conserving scheme                ln_dynvor_ene = ', ln_dynvor_ene
(  698)          WRITE(numout,*) '                 mixed enstrophy/energy conserving scheme   ln_dynvor_mix = ', ln_dynvor_mix
(  699)          WRITE(numout,*) '                 enstrophy and energy conserving scheme     ln_dynvor_een = ', ln_dynvor_een
(  700)       ENDIF
(  701) 
(  702)       ioptio = 0
(  703) 
(  704)       IF( ln_dynvor_ens ) THEN
(  705)          IF(lwp) WRITE(numout,*)
(  706)          IF(lwp) WRITE(numout,*) '              vorticity term : enstrophy conserving scheme'
(  707)          ioptio = ioptio + 1
(  708)       ENDIF
(  709)       IF( ln_dynvor_ene ) THEN
(  710)          IF(lwp) WRITE(numout,*)
(  711)          IF(lwp) WRITE(numout,*) '              vorticity term : energy conserving scheme'
(  712)          ioptio = ioptio + 1
(  713)       ENDIF
(  714)       IF( ln_dynvor_mix ) THEN
(  715)          IF(lwp) WRITE(numout,*)
(  716)          IF(lwp) WRITE(numout,*) '              vorticity term : mixed enstrophy/energy conserving scheme'
(  717)          ioptio = ioptio + 1
(  718)       ENDIF
(  719)       IF( ln_dynvor_een ) THEN
(  720)          IF(lwp) WRITE(numout,*)
(  721)          IF(lwp) WRITE(numout,*) '              vorticity term : energy and enstrophy conserving scheme'
(  722)          ioptio = ioptio + 1
(  723)       ENDIF
(  724)       IF ( ioptio /= 1 .AND. .NOT. lk_esopa ) THEN
(  725)           if(lwp) WRITE(numout,cform_err)
(  726)           IF(lwp) WRITE(numout,*) ' use ONE and ONLY one vorticity scheme'
(  727)           nstop = nstop + 1
(  728)       ENDIF
(  729) 
(  730)    END SUBROUTINE dyn_vor_ctl
(  731) 
(  732) !!==============================================================================
(  733) END MODULE dynvor
