


PGF90 (Version     10.5)          09/23/2010  09:49:44      page 1

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

Filename: dynspg_ts.F90

(    1) MODULE dynspg_ts
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_ts  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) #if ( defined key_dynspg_ts && ! defined key_autotasking ) ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_dynspg_ts'     free surface cst volume with time splitting
(    9)    !!   NOT 'key_autotasking'                      k-j-i loop (vector opt.)
(   10)    !!----------------------------------------------------------------------
(   11)    !!   dyn_spg_ts  : compute surface pressure gradient trend using a time-
(   12)    !!                 splitting scheme and add to the general trend 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE phycst          ! physical constants
(   18)    USE ocesbc          ! ocean surface boundary condition
(   19)    USE obcdta          ! open boundary condition data     
(   20)    USE obcfla          ! Flather open boundary condition  
(   21)    USE dynvor          ! vorticity term
(   22)    USE obc_oce         ! Lateral open boundary condition
(   23)    USE obc_par         ! open boundary condition parameters
(   24)    USE lib_mpp         ! distributed memory computing library
(   25)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   26)    USE prtctl          ! Print control
(   27)    USE dynspg_oce      ! surface pressure gradient variables
(   28)    USE in_out_manager  ! I/O manager
(   29) #ifdef key_RIVER_INPUT
(   30) !!DB
(   31)    USE rivers
(   32) #endif
(   33) 
(   34) 
(   35)    IMPLICIT NONE
(   36)    PRIVATE
(   37) 
(   38)    !! * Accessibility
(   39)    PUBLIC dyn_spg_ts  ! routine called by step.F90
(   40) 
(   41)    !! * Substitutions
(   42) #  include "domzgr_substitute.h90"
(   43) #  include "vectopt_loop_substitute.h90"
(   44)    !!----------------------------------------------------------------------
(   45)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   46)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_ts.F90,v 1.6 2006/01/03 15:04:14 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 2

(   47)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   48)    !!----------------------------------------------------------------------
(   49) 
(   50) CONTAINS
(   51) 
(   52)    SUBROUTINE dyn_spg_ts( kt )
(   53)       !!----------------------------------------------------------------------
(   54)       !!                  ***  routine dyn_spg_ts  ***
(   55)       !!
(   56)       !! ** Purpose :   Compute the now trend due to the surface pressure
(   57)       !!      gradient in case of free surface formulation with time-splitting.
(   58)       !!      Add it to the general trend of momentum equation.
(   59)       !!      Compute the free surface.
(   60)       !!
(   61)       !! ** Method  :   Free surface formulation with time-splitting
(   62)       !!      -1- Save the vertically integrated trend. This general trend is
(   63)       !!          held constant over the barotropic integration.
(   64)       !!          The Coriolis force is removed from the general trend as the
(   65)       !!          surface gradient and the Coriolis force are updated within
(   66)       !!          the barotropic integration.
(   67)       !!      -2- Barotropic loop : updates of sea surface height (ssha_e) and 
(   68)       !!          barotropic transports (ua_e and va_e) through barotropic 
(   69)       !!          momentum and continuity integration. Barotropic former 
(   70)       !!          variables are time averaging over the full barotropic cycle
(   71)       !!          (= 2 * baroclinic time step) and saved in zsshX_b, zuX_b 
(   72)       !!          and zvX_b (X specifying after, now or before).
(   73)       !!      -3- Update of sea surface height from time averaged barotropic 
(   74)       !!          variables.
(   75)       !!        - apply lateral boundary conditions on sshn.
(   76)       !!      -4- The new general trend becomes :
(   77)       !!          ua = ua - sum_k(ua)/H + ( zua_b - sum_k(ub) )/H
(   78)       !!
(   79)       !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
(   80)       !!
(   81)       !! References :
(   82)       !!   Griffies et al., (2003): A technical guide to MOM4. NOAA/GFDL
(   83)       !!
(   84)       !! History :
(   85)       !!   9.0  !  04-12  (L. Bessieres, G. Madec)  Original code
(   86)       !!        !  05-11  (V. Garnier, G. Madec)  optimization
(   87)       !!---------------------------------------------------------------------
(   88) 
(   89) 
(   90)       !! * Arguments
(   91)      INTEGER, INTENT( in )  ::   kt           ! ocean time-step index
(   92) 
(   93)       !! * Local declarations
(   94)       INTEGER  ::  ji, jj, jk, jit             ! dummy loop indices
(   95)       INTEGER  ::  icycle                      ! temporary scalar
(   96)       REAL(wp) ::                           &
(   97)          zraur, zcoef, z2dt_e, z2dt_b, zfac25,   &  ! temporary scalars
(   98)          zfact1, zspgu, zcubt, zx1, zy1,    &  !     "        "
(   99)          zfact2, zspgv, zcvbt, zx2, zy2        !     "        "
(  100)       REAL(wp), DIMENSION(jpi,jpj) ::       &
(  101)          zcu, zcv, zwx, zwy, zhdiv,         &  ! temporary arrays
(  102)          zua, zva, zub, zvb,                &  !     "        "
(  103)          zssha_b, zua_b, zva_b,             &  !     "        "
(  104)          zsshb_e, zub_e, zvb_e,             &  !     "        "






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 3

(  105)          zun_e, zvn_e                          !     "        "
(  106)       REAL(wp), DIMENSION(jpi,jpj),SAVE ::  &
(  107)          ztnw, ztne, ztsw, ztse
(  108)       !!----------------------------------------------------------------------
(  109) 
(  110)       ! Arrays initialization
(  111)       ! ---------------------
(  112)       zua_b(:,:) = 0.e0   ;   zub_e(:,:) = 0.e0   ;   zun_e(:,:) = 0.e0
(  113)       zva_b(:,:) = 0.e0   ;   zvb_e(:,:) = 0.e0   ;   zvn_e(:,:) = 0.e0
(  114)       zhdiv(:,:) = 0.e0
(  115) 
(  116) 
(  117)       IF( kt == nit000 ) THEN
(  118) 
(  119)          IF(lwp) WRITE(numout,*)
(  120)          IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
(  121)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
(  122)          IF(lwp) WRITE(numout,*) ' Number of sub cycle in 1 time-step (2 rdt) : icycle = ', FLOOR( 2*rdt/rdtbt )
(  123) 
(  124)          IF( .NOT. ln_rstart ) THEN
(  125)             ! initialize barotropic specific arrays
(  126)             sshb_b(:,:) = sshb(:,:)
(  127)             sshn_b(:,:) = sshn(:,:)
(  128)             un_b(:,:)   = 0.e0
(  129)             vn_b(:,:)   = 0.e0
(  130)             ! vertical sum
(  131)             IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
(  132)                DO jk = 1, jpkm1
(  133)                   DO ji = 1, jpij
(  134)                      un_b(ji,1) = un_b(ji,1) + fse3u(ji,1,jk) * un(ji,1,jk)
(  135)                      vn_b(ji,1) = vn_b(ji,1) + fse3v(ji,1,jk) * vn(ji,1,jk)
(  136)                   END DO
(  137)                END DO
(  138)             ELSE                             ! No  vector opt.
(  139)                DO jk = 1, jpkm1
(  140)                   un_b(:,:) = un_b(:,:) + fse3u(:,:,jk) * un(:,:,jk)
(  141)                   vn_b(:,:) = vn_b(:,:) + fse3v(:,:,jk) * vn(:,:,jk)
(  142)                END DO
(  143)             ENDIF
(  144)          ENDIF
(  145)          ssha_e(:,:) = sshn(:,:)
(  146)          ua_e(:,:)   = un_b(:,:)
(  147)          va_e(:,:)   = vn_b(:,:)
(  148) 
(  149)          IF( ln_dynvor_een ) THEN
(  150)             ztne(1,:) = 0.e0   ;   ztnw(1,:) = 0.e0   ;   ztse(1,:) = 0.e0   ;   ztsw(1,:) = 0.e0
(  151)             DO jj = 2, jpj
(  152)                DO ji = fs_2, jpi   ! vector opt.
(  153)                   ztne(ji,jj) = ( ff(ji-1,jj  ) + ff(ji  ,jj  ) + ff(ji  ,jj-1) ) / 3.
(  154)                   ztnw(ji,jj) = ( ff(ji-1,jj-1) + ff(ji-1,jj  ) + ff(ji  ,jj  ) ) / 3.
(  155)                   ztse(ji,jj) = ( ff(ji  ,jj  ) + ff(ji  ,jj-1) + ff(ji-1,jj-1) ) / 3.
(  156)                   ztsw(ji,jj) = ( ff(ji  ,jj-1) + ff(ji-1,jj-1) + ff(ji-1,jj  ) ) / 3.
(  157)                END DO
(  158)             END DO
(  159)          ENDIF
(  160) 
(  161)       ENDIF      !!kt==nit000
(  162)     






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 4

(  163)       ! Local constant initialization
(  164)       ! --------------------------------
(  165)       z2dt_b = 2.0 * rdt                                    ! baroclinic time step
(  166)       IF ( neuler == 0 .AND. kt == nit000 ) z2dt_b = rdt
(  167)       zfact1 = 0.5 * 0.25                                   ! coefficient for vorticity estimates
(  168)       zfact2 = 0.5 * 0.5
(  169)       zraur  = 1. / rauw                                    ! 1 / volumic mass of pure water
(  170)       
(  171)       ! -----------------------------------------------------------------------------
(  172)       !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
(  173)       ! -----------------------------------------------------------------------------
(  174) 
(  175)       ! Vertically integrated quantities
(  176)       ! --------------------------------
(  177)       zua(:,:) = 0.e0
(  178)       zva(:,:) = 0.e0
(  179)       zub(:,:) = 0.e0
(  180)       zvb(:,:) = 0.e0
(  181)       zwx(:,:) = 0.e0
(  182)       zwy(:,:) = 0.e0
(  183) 
(  184)       ! vertical sum
(  185)       IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
(  186)          DO jk = 1, jpkm1
(  187)             DO ji = 1, jpij
(  188)                !                                                           ! Vertically integrated momentum trends
(  189)                zua(ji,1) = zua(ji,1) + fse3u(ji,1,jk) * umask(ji,1,jk) * ua(ji,1,jk)
(  190)                zva(ji,1) = zva(ji,1) + fse3v(ji,1,jk) * vmask(ji,1,jk) * va(ji,1,jk)
(  191)                !                                                           ! Vertically integrated transports (before)
(  192)                zub(ji,1) = zub(ji,1) + fse3u(ji,1,jk) * ub(ji,1,jk)
(  193)                zvb(ji,1) = zvb(ji,1) + fse3v(ji,1,jk) * vb(ji,1,jk)
(  194)                !                                                           ! Planetary vorticity transport fluxes (now)
(  195)                zwx(ji,1) = zwx(ji,1) + e2u(ji,1) * fse3u(ji,1,jk) * un(ji,1,jk)
(  196)                zwy(ji,1) = zwy(ji,1) + e1v(ji,1) * fse3v(ji,1,jk) * vn(ji,1,jk)
(  197)             END DO
(  198)          END DO
(  199)       ELSE                             ! No  vector opt.
(  200)          DO jk = 1, jpkm1
(  201)             !                                                           ! Vertically integrated momentum trends
(  202)             zua(:,:) = zua(:,:) + fse3u(:,:,jk) * umask(:,:,jk) * ua(:,:,jk)
(  203)             zva(:,:) = zva(:,:) + fse3v(:,:,jk) * vmask(:,:,jk) * va(:,:,jk)
(  204)             !                                                           ! Vertically integrated transports (before)
(  205)             zub(:,:) = zub(:,:) + fse3u(:,:,jk) * ub(:,:,jk)
(  206)             zvb(:,:) = zvb(:,:) + fse3v(:,:,jk) * vb(:,:,jk)
(  207)             !                                                           ! Planetary vorticity (now)
(  208)             zwx(:,:) = zwx(:,:) + e2u(:,:) * fse3u(:,:,jk) * un(:,:,jk)
(  209)             zwy(:,:) = zwy(:,:) + e1v(:,:) * fse3v(:,:,jk) * vn(:,:,jk)
(  210)          END DO
(  211)       ENDIF
(  212) 
(  213) #ifdef key_RIVER_INPUT
(  214) !!DB: Adjust the above quantities to reflect runoff at head of SLE
(  215) !!     Works, but not clear whether this frag makes a difference -- see below
(  216)       call river_SLE_02(zua,zub,zwx)
(  217) #endif
(  218) 
(  219)       IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
(  220)          DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 5

(  221)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  222)                zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  223)                zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  224)                zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  225)                zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  226)                ! energy conserving formulation for planetary vorticity term
(  227)                zcu(ji,jj) = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
(  228)                zcv(ji,jj) =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
(  229)             END DO
(  230)          END DO
(  231) 
(  232)       ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
(  233)          DO jj = 2, jpjm1
(  234)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  235)                zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
(  236)                               + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  237)                zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
(  238)                               + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  239)                zcu(ji,jj)  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
(  240)                zcv(ji,jj)  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
(  241)             END DO
(  242)          END DO
(  243) 
(  244)       ELSEIF ( ln_dynvor_een ) THEN                    ! enstrophy and energy conserving scheme
(  245)          zfac25 = 0.25
(  246)          DO jj = 2, jpjm1
(  247)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  248)                zcu(ji,jj) = + zfac25 / e1u(ji,jj)   &
(  249)                   &       * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
(  250)                   &          + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
(  251)                zcv(ji,jj) = - zfac25 / e2v(ji,jj)   &
(  252)                   &       * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
(  253)                   &          + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
(  254)             END DO
(  255)          END DO
(  256) 
(  257)       ENDIF
(  258) 
(  259) 
(  260)       ! Remove barotropic trend from general momentum trend
(  261)       ! ---------------------------------------------------
(  262)       DO jk = 1 , jpkm1
(  263)          DO jj = 2, jpjm1
(  264)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  265)                ua(ji,jj,jk) = ua(ji,jj,jk) - zua(ji,jj) * hur(ji,jj)
(  266)                va(ji,jj,jk) = va(ji,jj,jk) - zva(ji,jj) * hvr(ji,jj)
(  267)             END DO
(  268)          END DO
(  269)       END DO
(  270) 
(  271)       ! Remove coriolis term from barotropic trend
(  272)       ! ------------------------------------------
(  273)       DO jj = 2, jpjm1
(  274)          DO ji = fs_2, fs_jpim1
(  275)             zua(ji,jj) = zua(ji,jj) - zcu(ji,jj)
(  276)             zva(ji,jj) = zva(ji,jj) - zcv(ji,jj)
(  277)          END DO
(  278)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 6

(  279) 
(  280)       ! -----------------------------------------------------------------------
(  281)       !  Phase 2 : Integration of the barotropic equations with time splitting
(  282)       ! -----------------------------------------------------------------------
(  283) 
(  284)       ! Initialisations
(  285)       !----------------
(  286)       ! Number of iteration of the barotropic loop
(  287)       icycle = FLOOR( z2dt_b / rdtbt )
(  288) 
(  289)       ! variables for the barotropic equations
(  290)       zsshb_e(:,:) = sshn_b(:,:)       ! (barotropic) sea surface height (before and now)
(  291)       sshn_e (:,:) = sshn_b(:,:)
(  292)       zub_e  (:,:) = un_b  (:,:)       ! barotropic transports issued from the barotropic equations (before and now)
(  293)       zvb_e  (:,:) = vn_b  (:,:)
(  294)       zun_e  (:,:) = un_b  (:,:)
(  295)       zvn_e  (:,:) = vn_b  (:,:)
(  296)       zssha_b(:,:) = sshn  (:,:)       ! time averaged variables over all sub-timesteps
(  297)       zua_b  (:,:) = un_b  (:,:)   
(  298)       zva_b  (:,:) = vn_b  (:,:)
(  299) 
(  300)       ! set ssh corrections to 0
(  301)       ! ssh corrections are applied to normal velocities (Flather's algorithm) and averaged over the barotropic loop
(  302) #if defined key_obc
(  303)       IF( lp_obc_east  )   sshfoe_b(:,:) = 0.e0
(  304)       IF( lp_obc_west  )   sshfow_b(:,:) = 0.e0
(  305)       IF( lp_obc_south )   sshfos_b(:,:) = 0.e0
(  306)       IF( lp_obc_north )   sshfon_b(:,:) = 0.e0
(  307) #endif
(  308) 
(  309)       ! Barotropic integration over 2 baroclinic time steps
(  310)       ! ---------------------------------------------------
(  311) 
(  312)       !                                                    ! ==================== !
(  313)       DO jit = 1, icycle                                   !  sub-time-step loop  !
(  314)          !                                                 ! ==================== !
(  315) 
(  316)          z2dt_e = 2. * rdtbt
(  317)          IF ( jit == 1 )   z2dt_e = rdtbt
(  318) 
(  319)          ! Time interpolation of open boundary condition data
(  320)          IF( lk_obc )   CALL obc_dta_bt( kt, jit )
(  321) 
(  322)          ! Horizontal divergence of barotropic transports
(  323)          !--------------------------------------------------
(  324)          DO jj = 2, jpjm1
(  325)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  326)                zhdiv(ji,jj) = ( e2u(ji  ,jj  ) * zun_e(ji  ,jj)              &
(  327)                   &            -e2u(ji-1,jj  ) * zun_e(ji-1,jj)              &
(  328)                   &            +e1v(ji  ,jj  ) * zvn_e(ji  ,jj)              &
(  329)                   &            -e1v(ji  ,jj-1) * zvn_e(ji  ,jj-1) )          &
(  330)                   &           / (e1t(ji,jj)*e2t(ji,jj))
(  331)             END DO
(  332)          END DO
(  333) 
(  334) #if defined key_obc
(  335)          ! open boundaries (div must be zero behind the open boundary)
(  336)          !  mpp remark: The zeroing of hdiv can probably be extended to 1->jpi/jpj for the correct row/column






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 7

(  337)          IF( lp_obc_east  )   zhdiv(nie0p1:nie1p1,nje0  :nje1)   = 0.e0      ! east
(  338)          IF( lp_obc_west  )   zhdiv(niw0  :niw1  ,njw0  :njw1)   = 0.e0      ! west
(  339)          IF( lp_obc_north )   zhdiv(nin0  :nin1  ,njn0p1:njn1p1) = 0.e0      ! north
(  340)          IF( lp_obc_south )   zhdiv(nis0  :nis1  ,njs0  :njs1)   = 0.e0      ! south
(  341) #endif
(  342) 
(  343)          ! Sea surface height from the barotropic system
(  344)          !----------------------------------------------
(  345)          DO jj = 2, jpjm1
(  346)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  347)                ssha_e(ji,jj) = ( zsshb_e(ji,jj) - z2dt_e *  ( zraur * emp(ji,jj)  &
(  348)             &            +  zhdiv(ji,jj) ) ) * tmask(ji,jj,1)
(  349)             END DO
(  350)          END DO
(  351) 
(  352)          ! evolution of the barotropic transport ( following the vorticity scheme used)
(  353)          ! ----------------------------------------------------------------------------
(  354)          zwx(:,:) = e2u(:,:) * zun_e(:,:)
(  355)          zwy(:,:) = e1v(:,:) * zvn_e(:,:)
(  356) 
(  357)          IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
(  358)             DO jj = 2, jpjm1
(  359)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  360)                   ! surface pressure gradient
(  361)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  362)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  363)                   ! energy conserving formulation for planetary vorticity term
(  364)                   zy1 = ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  365)                   zy2 = ( zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  366)                   zx1 = ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  367)                   zx2 = ( zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  368)                   zcubt = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
(  369)                   zcvbt =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
(  370)                   ! after transports
(  371)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  372)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  373)                END DO
(  374)             END DO
(  375) 
(  376)          ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
(  377)             DO jj = 2, jpjm1
(  378)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  379)                   ! surface pressure gradient
(  380)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  381)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  382)                   ! enstrophy conserving formulation for planetary vorticity term
(  383)                   zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
(  384)                                  + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  385)                   zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
(  386)                                  + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  387)                   zcubt  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
(  388)                   zcvbt  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
(  389)                   ! after transports
(  390)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  391)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  392)                END DO
(  393)             END DO
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 8

(  395)          ELSEIF ( ln_dynvor_een ) THEN                    ! energy and enstrophy conserving scheme
(  396)             zfac25 = 0.25
(  397)             DO jj = 2, jpjm1
(  398)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  399)                   ! surface pressure gradient
(  400)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  401)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  402)                   ! energy/enstrophy conserving formulation for planetary vorticity term
(  403)                   zcubt = + zfac25 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
(  404)                      &                             + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
(  405)                   zcvbt = - zfac25 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
(  406)                      &                             + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
(  407)                   ! after transports
(  408)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  409)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  410)                END DO
(  411)             END DO
(  412) 
(  413)          ENDIF
(  414) 
(  415) !!DB: add St. Lawrence River runoff to barotropic velocity
(  416) #ifdef key_RIVER_INPUT
(  417) !!DB: add trans to SLE via ua_e
(  418)          call river_SLE_03(ua_e)
(  419) #endif
(  420) 
(  421) 
(  422)          ! Flather's boundary condition for the barotropic loop :
(  423)          !         - Update sea surface height on each open boundary
(  424)          !         - Correct the barotropic transports
(  425)          IF( lk_obc )   CALL obc_fla_ts
(  426) 
(  427) 
(  428)          ! ... Boundary conditions on ua_e, va_e, ssha_e
(  429)          CALL lbc_lnk( ua_e  , 'U', -1. )
(  430)          CALL lbc_lnk( va_e  , 'V', -1. )
(  431)          CALL lbc_lnk( ssha_e, 'T',  1. )
(  432) 
(  433)          ! temporal sum
(  434)          !-------------
(  435)          zssha_b(:,:) = zssha_b(:,:) + ssha_e(:,:)
(  436)          zua_b  (:,:) = zua_b  (:,:) + ua_e  (:,:)
(  437)          zva_b  (:,:) = zva_b  (:,:) + va_e  (:,:) 
(  438) 
(  439)          ! Time filter and swap of dynamics arrays
(  440)          ! ---------------------------------------
(  441)          IF( neuler == 0 .AND. kt == nit000 ) THEN   ! Euler (forward) time stepping
(  442)             zsshb_e(:,:) = sshn_e(:,:)
(  443)             zub_e  (:,:) = zun_e (:,:)
(  444)             zvb_e  (:,:) = zvn_e (:,:)
(  445)             sshn_e (:,:) = ssha_e(:,:)
(  446)             zun_e  (:,:) = ua_e  (:,:)
(  447)             zvn_e  (:,:) = va_e  (:,:)
(  448)          ELSE                                        ! Asselin filtering
(  449)             zsshb_e(:,:) = atfp * ( zsshb_e(:,:) + ssha_e(:,:) ) + atfp1 * sshn_e(:,:)
(  450)             zub_e  (:,:) = atfp * ( zub_e  (:,:) + ua_e  (:,:) ) + atfp1 * zun_e  (:,:)
(  451)             zvb_e  (:,:) = atfp * ( zvb_e  (:,:) + va_e  (:,:) ) + atfp1 * zvn_e  (:,:)
(  452)             sshn_e (:,:) = ssha_e(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 9

(  453)             zun_e  (:,:) = ua_e  (:,:)
(  454)             zvn_e  (:,:) = va_e  (:,:)
(  455)          ENDIF
(  456) 
(  457)          !                                               ! ==================== !
(  458)       END DO       !   END OF  do jit = 1, icycle        !  sub-time-step loop  !
(  459)       !                                                  ! ==================== !
(  460) 
(  461) 
(  462)       ! Time average of after barotropic variables
(  463)       zcoef =  1.e0 / (  FLOAT( icycle +1 )  )
(  464)       zssha_b(:,:) = zcoef * zssha_b(:,:) 
(  465)       zua_b  (:,:) = zcoef *  zua_b (:,:) 
(  466)       zva_b  (:,:) = zcoef *  zva_b (:,:) 
(  467) #if defined key_obc
(  468)          IF( lp_obc_east  )   sshfoe_b(:,:) = zcoef * sshfoe_b(:,:)
(  469)          IF( lp_obc_west  )   sshfow_b(:,:) = zcoef * sshfow_b(:,:)
(  470)          IF( lp_obc_north )   sshfon_b(:,:) = zcoef * sshfon_b(:,:)
(  471)          IF( lp_obc_south )   sshfos_b(:,:) = zcoef * sshfos_b(:,:)
(  472) #endif
(  473)      
(  474) 
(  475)       ! ---------------------------------------------------------------------------
(  476)       ! Phase 3 : Update sea surface height from time averaged barotropic variables
(  477)       ! ---------------------------------------------------------------------------
(  478) 
(  479)  
(  480)       ! Horizontal divergence of time averaged barotropic transports
(  481)       !-------------------------------------------------------------
(  482)       DO jj = 2, jpjm1
(  483)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  484)             zhdiv(ji,jj) = ( e2u(ji,jj) * un_b(ji,jj) - e2u(ji-1,jj  ) * un_b(ji-1,jj  )     &
(  485)            &                +e1v(ji,jj) * vn_b(ji,jj) - e1v(ji  ,jj-1) * vn_b(ji  ,jj-1) )   &
(  486)            &             / ( e1t(ji,jj) * e2t(ji,jj) )
(  487)          END DO
(  488)       END DO
(  489) 
(  490) #if defined key_obc
(  491)       ! open boundaries (div must be zero behind the open boundary)
(  492)       !  mpp remark: The zeroing of hdiv can probably be extended to 1->jpi/jpj for the correct row/column
(  493)       IF( lp_obc_east  )   zhdiv(nie0p1:nie1p1,nje0  :nje1)   = 0.e0    ! east
(  494)       IF( lp_obc_west  )   zhdiv(niw0  :niw1  ,njw0  :njw1)   = 0.e0    ! west
(  495)       IF( lp_obc_north )   zhdiv(nin0  :nin1  ,njn0p1:njn1p1) = 0.e0    ! north
(  496)       IF( lp_obc_south )   zhdiv(nis0  :nis1  ,njs0  :njs1)   = 0.e0    ! south
(  497) #endif
(  498) 
(  499)       ! sea surface height
(  500)       !-------------------
(  501)       sshb(:,:) = sshn(:,:)
(  502)       sshn(:,:) = (  sshb_b(:,:) - z2dt_b * ( zraur * emp(:,:) + zhdiv(:,:) )  ) * tmask(:,:,1)
(  503) 
(  504)       ! ... Boundary conditions on sshn
(  505)       IF( .NOT. lk_obc ) CALL lbc_lnk( sshn, 'T', 1. )
(  506) 
(  507) 
(  508)       ! -----------------------------------------------------------------------------
(  509)       ! Phase 4. Coupling between general trend and barotropic estimates - (2nd step)
(  510)       ! -----------------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 10

(  511) 
(  512)       ! Swap on time averaged barotropic variables
(  513)       ! ------------------------------------------
(  514)       sshb_b(:,:) = sshn_b (:,:)
(  515)       sshn_b(:,:) = zssha_b(:,:)
(  516)       un_b  (:,:) = zua_b  (:,:) 
(  517)       vn_b  (:,:) = zva_b  (:,:) 
(  518)    
(  519)       ! add time averaged barotropic coriolis and surface pressure gradient
(  520)       ! terms to the general momentum trend
(  521)       ! --------------------------------------------------------------------
(  522)       DO jk=1,jpkm1
(  523)          ua(:,:,jk) = ua(:,:,jk) + hur(:,:) * ( zua_b(:,:) - zub(:,:) ) / z2dt_b
(  524)          va(:,:,jk) = va(:,:,jk) + hvr(:,:) * ( zva_b(:,:) - zvb(:,:) ) / z2dt_b
(  525)       END DO
(  526) 
(  527) 
(  528)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  529)          CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
(  530)       ENDIF
(  531)       
(  532)    END SUBROUTINE dyn_spg_ts
(  533) #else
(  534)    !!----------------------------------------------------------------------
(  535)    !!   Default case :   Empty module   No standart free surface cst volume
(  536)    !!----------------------------------------------------------------------
(  537)    USE in_out_manager
(  538) CONTAINS
(  539)    SUBROUTINE dyn_spg_ts( kt )       ! Empty routine
(  540)       if(lwp) WRITE(numout,*) 'dyn_spg_ts: You should not have seen this print! error?', kt
(  541)    END SUBROUTINE dyn_spg_ts
(  542) #endif
(  543)    
(  544)    !!======================================================================
(  545) END MODULE dynspg_ts





























PGF90 (Version     10.5)          09/23/2010  09:49:44      page 11

(    1) # 1 "dynspg_ts.F90"
(    1) MODULE dynspg_ts
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_ts  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_dynspg_ts'     free surface cst volume with time splitting
(    9)    !!   NOT 'key_autotasking'                      k-j-i loop (vector opt.)
(   10)    !!----------------------------------------------------------------------
(   11)    !!   dyn_spg_ts  : compute surface pressure gradient trend using a time-
(   12)    !!                 splitting scheme and add to the general trend 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE phycst          ! physical constants
(   18)    USE ocesbc          ! ocean surface boundary condition
(   19)    USE obcdta          ! open boundary condition data     
(   20)    USE obcfla          ! Flather open boundary condition  
(   21)    USE dynvor          ! vorticity term
(   22)    USE obc_oce         ! Lateral open boundary condition
(   23)    USE obc_par         ! open boundary condition parameters
(   24)    USE lib_mpp         ! distributed memory computing library
(   25)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   26)    USE prtctl          ! Print control
(   27)    USE dynspg_oce      ! surface pressure gradient variables
(   28)    USE in_out_manager  ! I/O manager
(   29) # 30
(   30) !!DB
(   31)    USE rivers
(   32) 
(   33) 
(   34) # 35
(   35)    IMPLICIT NONE
(   36)    PRIVATE
(   37) 
(   38)    !! * Accessibility
(   39)    PUBLIC dyn_spg_ts  ! routine called by step.F90
(   40) 
(   41)    !! * Substitutions
(   42) # 1 "./domzgr_substitute.h90"
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
(   48) # 43 "dynspg_ts.F90"






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 12

(   43) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 44 "dynspg_ts.F90"
(   44) # 44
(   44)    !!----------------------------------------------------------------------
(   45)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   46)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_ts.F90,v 1.6 2006/01/03 15:04:14 opalod Exp $ 
(   47)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   48)    !!----------------------------------------------------------------------
(   49) 
(   50) CONTAINS
(   51) 
(   52)    SUBROUTINE dyn_spg_ts( kt )
(   53)       !!----------------------------------------------------------------------
(   54)       !!                  ***  routine dyn_spg_ts  ***
(   55)       !!
(   56)       !! ** Purpose :   Compute the now trend due to the surface pressure
(   57)       !!      gradient in case of free surface formulation with time-splitting.
(   58)       !!      Add it to the general trend of momentum equation.
(   59)       !!      Compute the free surface.
(   60)       !!
(   61)       !! ** Method  :   Free surface formulation with time-splitting
(   62)       !!      -1- Save the vertically integrated trend. This general trend is
(   63)       !!          held constant over the barotropic integration.
(   64)       !!          The Coriolis force is removed from the general trend as the
(   65)       !!          surface gradient and the Coriolis force are updated within
(   66)       !!          the barotropic integration.
(   67)       !!      -2- Barotropic loop : updates of sea surface height (ssha_e) and 
(   68)       !!          barotropic transports (ua_e and va_e) through barotropic 
(   69)       !!          momentum and continuity integration. Barotropic former 
(   70)       !!          variables are time averaging over the full barotropic cycle
(   71)       !!          (= 2 * baroclinic time step) and saved in zsshX_b, zuX_b 
(   72)       !!          and zvX_b (X specifying after, now or before).
(   73)       !!      -3- Update of sea surface height from time averaged barotropic 
(   74)       !!          variables.
(   75)       !!        - apply lateral boundary conditions on sshn.
(   76)       !!      -4- The new general trend becomes :
(   77)       !!          ua = ua - sum_k(ua)/H + ( zua_b - sum_k(ub) )/H
(   78)       !!
(   79)       !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
(   80)       !!
(   81)       !! References :
(   82)       !!   Griffies et al., (2003): A technical guide to MOM4. NOAA/GFDL
(   83)       !!
(   84)       !! History :
(   85)       !!   9.0  !  04-12  (L. Bessieres, G. Madec)  Original code
(   86)       !!        !  05-11  (V. Garnier, G. Madec)  optimization
(   87)       !!---------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 13

(   88) 
(   89) 
(   90)       !! * Arguments
(   91)      INTEGER, INTENT( in )  ::   kt           ! ocean time-step index
(   92) 
(   93)       !! * Local declarations
(   94)       INTEGER  ::  ji, jj, jk, jit             ! dummy loop indices
(   95)       INTEGER  ::  icycle                      ! temporary scalar
(   96)       REAL(wp) ::                           &
(   97)          zraur, zcoef, z2dt_e, z2dt_b, zfac25,   &  ! temporary scalars
(   98)          zfact1, zspgu, zcubt, zx1, zy1,    &  !     "        "
(   99)          zfact2, zspgv, zcvbt, zx2, zy2        !     "        "
(  100)       REAL(wp), DIMENSION(jpi,jpj) ::       &
(  101)          zcu, zcv, zwx, zwy, zhdiv,         &  ! temporary arrays
(  102)          zua, zva, zub, zvb,                &  !     "        "
(  103)          zssha_b, zua_b, zva_b,             &  !     "        "
(  104)          zsshb_e, zub_e, zvb_e,             &  !     "        "
(  105)          zun_e, zvn_e                          !     "        "
(  106)       REAL(wp), DIMENSION(jpi,jpj),SAVE ::  &
(  107)          ztnw, ztne, ztsw, ztse
(  108)       !!----------------------------------------------------------------------
(  109) 
(  110)       ! Arrays initialization
(  111)       ! ---------------------
(  112)       zua_b(:,:) = 0.e0   ;   zub_e(:,:) = 0.e0   ;   zun_e(:,:) = 0.e0
(  113)       zva_b(:,:) = 0.e0   ;   zvb_e(:,:) = 0.e0   ;   zvn_e(:,:) = 0.e0
(  114)       zhdiv(:,:) = 0.e0
(  115) 
(  116) 
(  117)       IF( kt == nit000 ) THEN
(  118) 
(  119)          IF(lwp) WRITE(numout,*)
(  120)          IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
(  121)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
(  122)          IF(lwp) WRITE(numout,*) ' Number of sub cycle in 1 time-step (2 rdt) : icycle = ', FLOOR( 2*rdt/rdtbt )
(  123) 
(  124)          IF( .NOT. ln_rstart ) THEN
(  125)             ! initialize barotropic specific arrays
(  126)             sshb_b(:,:) = sshb(:,:)
(  127)             sshn_b(:,:) = sshn(:,:)
(  128)             un_b(:,:)   = 0.e0
(  129)             vn_b(:,:)   = 0.e0
(  130)             ! vertical sum
(  131)             IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
(  132)                DO jk = 1, jpkm1
(  133)                   DO ji = 1, jpij
(  134)                      un_b(ji,1) = un_b(ji,1) + e3t(jk) * un(ji,1,jk)
(  135)                      vn_b(ji,1) = vn_b(ji,1) + e3t(jk) * vn(ji,1,jk)
(  136)                   END DO
(  137)                END DO
(  138)             ELSE                             ! No  vector opt.
(  139)                DO jk = 1, jpkm1
(  140)                   un_b(:,:) = un_b(:,:) + e3t(jk) * un(:,:,jk)
(  141)                   vn_b(:,:) = vn_b(:,:) + e3t(jk) * vn(:,:,jk)
(  142)                END DO
(  143)             ENDIF
(  144)          ENDIF
(  145)          ssha_e(:,:) = sshn(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 14

(  146)          ua_e(:,:)   = un_b(:,:)
(  147)          va_e(:,:)   = vn_b(:,:)
(  148) 
(  149)          IF( ln_dynvor_een ) THEN
(  150)             ztne(1,:) = 0.e0   ;   ztnw(1,:) = 0.e0   ;   ztse(1,:) = 0.e0   ;   ztsw(1,:) = 0.e0
(  151)             DO jj = 2, jpj
(  152)                DO ji = 2, jpi   ! vector opt.
(  153)                   ztne(ji,jj) = ( ff(ji-1,jj  ) + ff(ji  ,jj  ) + ff(ji  ,jj-1) ) / 3.
(  154)                   ztnw(ji,jj) = ( ff(ji-1,jj-1) + ff(ji-1,jj  ) + ff(ji  ,jj  ) ) / 3.
(  155)                   ztse(ji,jj) = ( ff(ji  ,jj  ) + ff(ji  ,jj-1) + ff(ji-1,jj-1) ) / 3.
(  156)                   ztsw(ji,jj) = ( ff(ji  ,jj-1) + ff(ji-1,jj-1) + ff(ji-1,jj  ) ) / 3.
(  157)                END DO
(  158)             END DO
(  159)          ENDIF
(  160) 
(  161)       ENDIF      !!kt==nit000
(  162)     
(  163)       ! Local constant initialization
(  164)       ! --------------------------------
(  165)       z2dt_b = 2.0 * rdt                                    ! baroclinic time step
(  166)       IF ( neuler == 0 .AND. kt == nit000 ) z2dt_b = rdt
(  167)       zfact1 = 0.5 * 0.25                                   ! coefficient for vorticity estimates
(  168)       zfact2 = 0.5 * 0.5
(  169)       zraur  = 1. / rauw                                    ! 1 / volumic mass of pure water
(  170)       
(  171)       ! -----------------------------------------------------------------------------
(  172)       !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
(  173)       ! -----------------------------------------------------------------------------
(  174) 
(  175)       ! Vertically integrated quantities
(  176)       ! --------------------------------
(  177)       zua(:,:) = 0.e0
(  178)       zva(:,:) = 0.e0
(  179)       zub(:,:) = 0.e0
(  180)       zvb(:,:) = 0.e0
(  181)       zwx(:,:) = 0.e0
(  182)       zwy(:,:) = 0.e0
(  183) 
(  184)       ! vertical sum
(  185)       IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
(  186)          DO jk = 1, jpkm1
(  187)             DO ji = 1, jpij
(  188)                !                                                           ! Vertically integrated momentum trends
(  189)                zua(ji,1) = zua(ji,1) + e3t(jk) * umask(ji,1,jk) * ua(ji,1,jk)
(  190)                zva(ji,1) = zva(ji,1) + e3t(jk) * vmask(ji,1,jk) * va(ji,1,jk)
(  191)                !                                                           ! Vertically integrated transports (before)
(  192)                zub(ji,1) = zub(ji,1) + e3t(jk) * ub(ji,1,jk)
(  193)                zvb(ji,1) = zvb(ji,1) + e3t(jk) * vb(ji,1,jk)
(  194)                !                                                           ! Planetary vorticity transport fluxes (now)
(  195)                zwx(ji,1) = zwx(ji,1) + e2u(ji,1) * e3t(jk) * un(ji,1,jk)
(  196)                zwy(ji,1) = zwy(ji,1) + e1v(ji,1) * e3t(jk) * vn(ji,1,jk)
(  197)             END DO
(  198)          END DO
(  199)       ELSE                             ! No  vector opt.
(  200)          DO jk = 1, jpkm1
(  201)             !                                                           ! Vertically integrated momentum trends
(  202)             zua(:,:) = zua(:,:) + e3t(jk) * umask(:,:,jk) * ua(:,:,jk)
(  203)             zva(:,:) = zva(:,:) + e3t(jk) * vmask(:,:,jk) * va(:,:,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 15

(  204)             !                                                           ! Vertically integrated transports (before)
(  205)             zub(:,:) = zub(:,:) + e3t(jk) * ub(:,:,jk)
(  206)             zvb(:,:) = zvb(:,:) + e3t(jk) * vb(:,:,jk)
(  207)             !                                                           ! Planetary vorticity (now)
(  208)             zwx(:,:) = zwx(:,:) + e2u(:,:) * e3t(jk) * un(:,:,jk)
(  209)             zwy(:,:) = zwy(:,:) + e1v(:,:) * e3t(jk) * vn(:,:,jk)
(  210)          END DO
(  211)       ENDIF
(  212) 
(  213) # 214
(  214) !!DB: Adjust the above quantities to reflect runoff at head of SLE
(  215) !!     Works, but not clear whether this frag makes a difference -- see below
(  216)       call river_SLE_02(zua,zub,zwx)
(  217) 
(  218) # 219
(  219)       IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
(  220)          DO jj = 2, jpjm1
(  221)             DO ji = 2, jpim1   ! vector opt.
(  222)                zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  223)                zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  224)                zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  225)                zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  226)                ! energy conserving formulation for planetary vorticity term
(  227)                zcu(ji,jj) = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
(  228)                zcv(ji,jj) =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
(  229)             END DO
(  230)          END DO
(  231) 
(  232)       ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
(  233)          DO jj = 2, jpjm1
(  234)             DO ji = 2, jpim1   ! vector opt.
(  235)                zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
(  236)                               + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  237)                zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
(  238)                               + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  239)                zcu(ji,jj)  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
(  240)                zcv(ji,jj)  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
(  241)             END DO
(  242)          END DO
(  243) 
(  244)       ELSEIF ( ln_dynvor_een ) THEN                    ! enstrophy and energy conserving scheme
(  245)          zfac25 = 0.25
(  246)          DO jj = 2, jpjm1
(  247)             DO ji = 2, jpim1   ! vector opt.
(  248)                zcu(ji,jj) = + zfac25 / e1u(ji,jj)   &
(  249)                   &       * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
(  250)                   &          + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
(  251)                zcv(ji,jj) = - zfac25 / e2v(ji,jj)   &
(  252)                   &       * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
(  253)                   &          + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
(  254)             END DO
(  255)          END DO
(  256) 
(  257)       ENDIF
(  258) 
(  259) 
(  260)       ! Remove barotropic trend from general momentum trend
(  261)       ! ---------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 16

(  262)       DO jk = 1 , jpkm1
(  263)          DO jj = 2, jpjm1
(  264)             DO ji = 2, jpim1   ! vector opt.
(  265)                ua(ji,jj,jk) = ua(ji,jj,jk) - zua(ji,jj) * hur(ji,jj)
(  266)                va(ji,jj,jk) = va(ji,jj,jk) - zva(ji,jj) * hvr(ji,jj)
(  267)             END DO
(  268)          END DO
(  269)       END DO
(  270) 
(  271)       ! Remove coriolis term from barotropic trend
(  272)       ! ------------------------------------------
(  273)       DO jj = 2, jpjm1
(  274)          DO ji = 2, jpim1
(  275)             zua(ji,jj) = zua(ji,jj) - zcu(ji,jj)
(  276)             zva(ji,jj) = zva(ji,jj) - zcv(ji,jj)
(  277)          END DO
(  278)       END DO
(  279) 
(  280)       ! -----------------------------------------------------------------------
(  281)       !  Phase 2 : Integration of the barotropic equations with time splitting
(  282)       ! -----------------------------------------------------------------------
(  283) 
(  284)       ! Initialisations
(  285)       !----------------
(  286)       ! Number of iteration of the barotropic loop
(  287)       icycle = FLOOR( z2dt_b / rdtbt )
(  288) 
(  289)       ! variables for the barotropic equations
(  290)       zsshb_e(:,:) = sshn_b(:,:)       ! (barotropic) sea surface height (before and now)
(  291)       sshn_e (:,:) = sshn_b(:,:)
(  292)       zub_e  (:,:) = un_b  (:,:)       ! barotropic transports issued from the barotropic equations (before and now)
(  293)       zvb_e  (:,:) = vn_b  (:,:)
(  294)       zun_e  (:,:) = un_b  (:,:)
(  295)       zvn_e  (:,:) = vn_b  (:,:)
(  296)       zssha_b(:,:) = sshn  (:,:)       ! time averaged variables over all sub-timesteps
(  297)       zua_b  (:,:) = un_b  (:,:)   
(  298)       zva_b  (:,:) = vn_b  (:,:)
(  299) 
(  300)       ! set ssh corrections to 0
(  301)       ! ssh corrections are applied to normal velocities (Flather's algorithm) and averaged over the barotropic loop
(  302) # 303
(  303)       IF( lp_obc_east  )   sshfoe_b(:,:) = 0.e0
(  304)       IF( lp_obc_west  )   sshfow_b(:,:) = 0.e0
(  305)       IF( lp_obc_south )   sshfos_b(:,:) = 0.e0
(  306)       IF( lp_obc_north )   sshfon_b(:,:) = 0.e0
(  307) 
(  308) # 309
(  309)       ! Barotropic integration over 2 baroclinic time steps
(  310)       ! ---------------------------------------------------
(  311) 
(  312)       !                                                    ! ==================== !
(  313)       DO jit = 1, icycle                                   !  sub-time-step loop  !
(  314)          !                                                 ! ==================== !
(  315) 
(  316)          z2dt_e = 2. * rdtbt
(  317)          IF ( jit == 1 )   z2dt_e = rdtbt
(  318) 
(  319)          ! Time interpolation of open boundary condition data






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 17

(  320)          IF( lk_obc )   CALL obc_dta_bt( kt, jit )
(  321) 
(  322)          ! Horizontal divergence of barotropic transports
(  323)          !--------------------------------------------------
(  324)          DO jj = 2, jpjm1
(  325)             DO ji = 2, jpim1   ! vector opt.
(  326)                zhdiv(ji,jj) = ( e2u(ji  ,jj  ) * zun_e(ji  ,jj)              &
(  327)                   &            -e2u(ji-1,jj  ) * zun_e(ji-1,jj)              &
(  328)                   &            +e1v(ji  ,jj  ) * zvn_e(ji  ,jj)              &
(  329)                   &            -e1v(ji  ,jj-1) * zvn_e(ji  ,jj-1) )          &
(  330)                   &           / (e1t(ji,jj)*e2t(ji,jj))
(  331)             END DO
(  332)          END DO
(  333) 
(  334) # 335
(  335)          ! open boundaries (div must be zero behind the open boundary)
(  336)          !  mpp remark: The zeroing of hdiv can probably be extended to 1->jpi/jpj for the correct row/column
(  337)          IF( lp_obc_east  )   zhdiv(nie0p1:nie1p1,nje0  :nje1)   = 0.e0      ! east
(  338)          IF( lp_obc_west  )   zhdiv(niw0  :niw1  ,njw0  :njw1)   = 0.e0      ! west
(  339)          IF( lp_obc_north )   zhdiv(nin0  :nin1  ,njn0p1:njn1p1) = 0.e0      ! north
(  340)          IF( lp_obc_south )   zhdiv(nis0  :nis1  ,njs0  :njs1)   = 0.e0      ! south
(  341) 
(  342) # 343
(  343)          ! Sea surface height from the barotropic system
(  344)          !----------------------------------------------
(  345)          DO jj = 2, jpjm1
(  346)             DO ji = 2, jpim1   ! vector opt.
(  347)                ssha_e(ji,jj) = ( zsshb_e(ji,jj) - z2dt_e *  ( zraur * emp(ji,jj)  &
(  348)             &            +  zhdiv(ji,jj) ) ) * tmask(ji,jj,1)
(  349)             END DO
(  350)          END DO
(  351) 
(  352)          ! evolution of the barotropic transport ( following the vorticity scheme used)
(  353)          ! ----------------------------------------------------------------------------
(  354)          zwx(:,:) = e2u(:,:) * zun_e(:,:)
(  355)          zwy(:,:) = e1v(:,:) * zvn_e(:,:)
(  356) 
(  357)          IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
(  358)             DO jj = 2, jpjm1
(  359)                DO ji = 2, jpim1   ! vector opt.
(  360)                   ! surface pressure gradient
(  361)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  362)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  363)                   ! energy conserving formulation for planetary vorticity term
(  364)                   zy1 = ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  365)                   zy2 = ( zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  366)                   zx1 = ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  367)                   zx2 = ( zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  368)                   zcubt = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
(  369)                   zcvbt =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
(  370)                   ! after transports
(  371)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  372)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  373)                END DO
(  374)             END DO
(  375) 
(  376)          ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
(  377)             DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 18

(  378)                DO ji = 2, jpim1   ! vector opt.
(  379)                   ! surface pressure gradient
(  380)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  381)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  382)                   ! enstrophy conserving formulation for planetary vorticity term
(  383)                   zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
(  384)                                  + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  385)                   zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
(  386)                                  + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  387)                   zcubt  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
(  388)                   zcvbt  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
(  389)                   ! after transports
(  390)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  391)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  392)                END DO
(  393)             END DO
(  394) 
(  395)          ELSEIF ( ln_dynvor_een ) THEN                    ! energy and enstrophy conserving scheme
(  396)             zfac25 = 0.25
(  397)             DO jj = 2, jpjm1
(  398)                DO ji = 2, jpim1   ! vector opt.
(  399)                   ! surface pressure gradient
(  400)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  401)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  402)                   ! energy/enstrophy conserving formulation for planetary vorticity term
(  403)                   zcubt = + zfac25 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
(  404)                      &                             + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
(  405)                   zcvbt = - zfac25 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
(  406)                      &                             + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
(  407)                   ! after transports
(  408)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  409)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  410)                END DO
(  411)             END DO
(  412) 
(  413)          ENDIF
(  414) 
(  415) !!DB: add St. Lawrence River runoff to barotropic velocity
(  416) # 417
(  417) !!DB: add trans to SLE via ua_e
(  418)          call river_SLE_03(ua_e)
(  419) 
(  420) 
(  421) # 422
(  422)          ! Flather's boundary condition for the barotropic loop :
(  423)          !         - Update sea surface height on each open boundary
(  424)          !         - Correct the barotropic transports
(  425)          IF( lk_obc )   CALL obc_fla_ts
(  426) 
(  427) 
(  428)          ! ... Boundary conditions on ua_e, va_e, ssha_e
(  429)          CALL lbc_lnk( ua_e  , 'U', -1. )
(  430)          CALL lbc_lnk( va_e  , 'V', -1. )
(  431)          CALL lbc_lnk( ssha_e, 'T',  1. )
(  432) 
(  433)          ! temporal sum
(  434)          !-------------
(  435)          zssha_b(:,:) = zssha_b(:,:) + ssha_e(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 19

(  436)          zua_b  (:,:) = zua_b  (:,:) + ua_e  (:,:)
(  437)          zva_b  (:,:) = zva_b  (:,:) + va_e  (:,:) 
(  438) 
(  439)          ! Time filter and swap of dynamics arrays
(  440)          ! ---------------------------------------
(  441)          IF( neuler == 0 .AND. kt == nit000 ) THEN   ! Euler (forward) time stepping
(  442)             zsshb_e(:,:) = sshn_e(:,:)
(  443)             zub_e  (:,:) = zun_e (:,:)
(  444)             zvb_e  (:,:) = zvn_e (:,:)
(  445)             sshn_e (:,:) = ssha_e(:,:)
(  446)             zun_e  (:,:) = ua_e  (:,:)
(  447)             zvn_e  (:,:) = va_e  (:,:)
(  448)          ELSE                                        ! Asselin filtering
(  449)             zsshb_e(:,:) = atfp * ( zsshb_e(:,:) + ssha_e(:,:) ) + atfp1 * sshn_e(:,:)
(  450)             zub_e  (:,:) = atfp * ( zub_e  (:,:) + ua_e  (:,:) ) + atfp1 * zun_e  (:,:)
(  451)             zvb_e  (:,:) = atfp * ( zvb_e  (:,:) + va_e  (:,:) ) + atfp1 * zvn_e  (:,:)
(  452)             sshn_e (:,:) = ssha_e(:,:)
(  453)             zun_e  (:,:) = ua_e  (:,:)
(  454)             zvn_e  (:,:) = va_e  (:,:)
(  455)          ENDIF
(  456) 
(  457)          !                                               ! ==================== !
(  458)       END DO       !   END OF  do jit = 1, icycle        !  sub-time-step loop  !
(  459)       !                                                  ! ==================== !
(  460) 
(  461) 
(  462)       ! Time average of after barotropic variables
(  463)       zcoef =  1.e0 / (  FLOAT( icycle +1 )  )
(  464)       zssha_b(:,:) = zcoef * zssha_b(:,:) 
(  465)       zua_b  (:,:) = zcoef *  zua_b (:,:) 
(  466)       zva_b  (:,:) = zcoef *  zva_b (:,:) 
(  467) # 468
(  468)          IF( lp_obc_east  )   sshfoe_b(:,:) = zcoef * sshfoe_b(:,:)
(  469)          IF( lp_obc_west  )   sshfow_b(:,:) = zcoef * sshfow_b(:,:)
(  470)          IF( lp_obc_north )   sshfon_b(:,:) = zcoef * sshfon_b(:,:)
(  471)          IF( lp_obc_south )   sshfos_b(:,:) = zcoef * sshfos_b(:,:)
(  472) # 473
(  473)      
(  474) 
(  475)       ! ---------------------------------------------------------------------------
(  476)       ! Phase 3 : Update sea surface height from time averaged barotropic variables
(  477)       ! ---------------------------------------------------------------------------
(  478) 
(  479)  
(  480)       ! Horizontal divergence of time averaged barotropic transports
(  481)       !-------------------------------------------------------------
(  482)       DO jj = 2, jpjm1
(  483)          DO ji = 2, jpim1   ! vector opt.
(  484)             zhdiv(ji,jj) = ( e2u(ji,jj) * un_b(ji,jj) - e2u(ji-1,jj  ) * un_b(ji-1,jj  )     &
(  485)            &                +e1v(ji,jj) * vn_b(ji,jj) - e1v(ji  ,jj-1) * vn_b(ji  ,jj-1) )   &
(  486)            &             / ( e1t(ji,jj) * e2t(ji,jj) )
(  487)          END DO
(  488)       END DO
(  489) 
(  490) # 491
(  491)       ! open boundaries (div must be zero behind the open boundary)
(  492)       !  mpp remark: The zeroing of hdiv can probably be extended to 1->jpi/jpj for the correct row/column
(  493)       IF( lp_obc_east  )   zhdiv(nie0p1:nie1p1,nje0  :nje1)   = 0.e0    ! east






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 20

(  494)       IF( lp_obc_west  )   zhdiv(niw0  :niw1  ,njw0  :njw1)   = 0.e0    ! west
(  495)       IF( lp_obc_north )   zhdiv(nin0  :nin1  ,njn0p1:njn1p1) = 0.e0    ! north
(  496)       IF( lp_obc_south )   zhdiv(nis0  :nis1  ,njs0  :njs1)   = 0.e0    ! south
(  497) 
(  498) # 499
(  499)       ! sea surface height
(  500)       !-------------------
(  501)       sshb(:,:) = sshn(:,:)
(  502)       sshn(:,:) = (  sshb_b(:,:) - z2dt_b * ( zraur * emp(:,:) + zhdiv(:,:) )  ) * tmask(:,:,1)
(  503) 
(  504)       ! ... Boundary conditions on sshn
(  505)       IF( .NOT. lk_obc ) CALL lbc_lnk( sshn, 'T', 1. )
(  506) 
(  507) 
(  508)       ! -----------------------------------------------------------------------------
(  509)       ! Phase 4. Coupling between general trend and barotropic estimates - (2nd step)
(  510)       ! -----------------------------------------------------------------------------
(  511) 
(  512)       ! Swap on time averaged barotropic variables
(  513)       ! ------------------------------------------
(  514)       sshb_b(:,:) = sshn_b (:,:)
(  515)       sshn_b(:,:) = zssha_b(:,:)
(  516)       un_b  (:,:) = zua_b  (:,:) 
(  517)       vn_b  (:,:) = zva_b  (:,:) 
(  518)    
(  519)       ! add time averaged barotropic coriolis and surface pressure gradient
(  520)       ! terms to the general momentum trend
(  521)       ! --------------------------------------------------------------------
(  522)       DO jk=1,jpkm1
(  523)          ua(:,:,jk) = ua(:,:,jk) + hur(:,:) * ( zua_b(:,:) - zub(:,:) ) / z2dt_b
(  524)          va(:,:,jk) = va(:,:,jk) + hvr(:,:) * ( zva_b(:,:) - zvb(:,:) ) / z2dt_b
(  525)       END DO
(  526) 
(  527) 
(  528)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  529)          CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
(  530)       ENDIF
(  531)       
(  532)    END SUBROUTINE dyn_spg_ts
(  533) # 543
(  543)    
(  544)    !!======================================================================
(  545) END MODULE dynspg_ts
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (dynspg_ts.F90: 463)
