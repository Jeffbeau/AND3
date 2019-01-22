


PGF90 (Version     10.5)          09/23/2010  09:49:43      page 1

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

Filename: dynspg_ts_jki.F90

(    1) MODULE dynspg_ts_jki
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_ts_jki  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) #if ( defined key_dynspg_ts && defined key_autotasking )   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_dynspg_ts'                    free surface with time splitting
(    9)    !!   'key_autotasking'                          j-k-i loop (vector opt.)
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
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33)    !! * Accessibility
(   34)    PUBLIC dyn_spg_ts_jki  ! routine called by step.F90
(   35) 
(   36)    !! * Substitutions
(   37) #  include "domzgr_substitute.h90"
(   38) #  include "vectopt_loop_substitute.h90"
(   39)    !!----------------------------------------------------------------------
(   40)    !!   OPA 9.0 , LODYC-IPSL  (2005)
(   41)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_ts_jki.F90,v 1.2 2006/01/03 15:04:14 opalod Exp $ 
(   42)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   43)    !!----------------------------------------------------------------------
(   44) 
(   45) CONTAINS
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 2

(   47)    SUBROUTINE dyn_spg_ts_jki( kt )
(   48)       !!----------------------------------------------------------------------
(   49)       !!                  ***  routine dyn_spg_ts_jki  ***
(   50)       !!
(   51)       !! ** Purpose :   Compute the now trend due to the surface pressure
(   52)       !!      gradient in case of free surface formulation with time-splitting.
(   53)       !!      Add it to the general trend of momentum equation.
(   54)       !!      Compute the free surface.
(   55)       !!
(   56)       !! ** Method  :   Free surface formulation with time-splitting
(   57)       !!      -1- Save the vertically integrated trend. This general trend is
(   58)       !!          held constant over the barotropic integration.
(   59)       !!          The Coriolis force is removed from the general trend as the
(   60)       !!          surface gradient and the Coriolis force are updated within
(   61)       !!          the barotropic integration.
(   62)       !!      -2- Barotropic loop : updates of sea surface height (ssha_e) and 
(   63)       !!          barotropic transports (ua_e and va_e) through barotropic 
(   64)       !!          momentum and continuity integration. Barotropic former 
(   65)       !!          variables are time averaging over the full barotropic cycle
(   66)       !!          (= 2 * baroclinic time step) and saved in zsshX_b, zuX_b 
(   67)       !!          and zvX_b (X specifying after, now or before).
(   68)       !!      -3- Update of sea surface height from time averaged barotropic 
(   69)       !!          variables.
(   70)       !!        - apply lateral boundary conditions on sshn.
(   71)       !!      -4- The new general trend becomes :
(   72)       !!          ua = ua - sum_k(ua)/H + ( zua_b - sum_k(ub) )/H
(   73)       !!
(   74)       !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
(   75)       !!
(   76)       !! References :
(   77)       !!   Griffies et al., (2003): A technical guide to MOM4. NOAA/GFDL
(   78)       !!
(   79)       !! History :
(   80)       !!   9.0  !  04-12  (L. Bessieres, G. Madec)  Original code
(   81)       !!        !  05-11  (V. Garnier, G. Madec)  optimization
(   82)       !!---------------------------------------------------------------------
(   83)       !! * Arguments
(   84)       INTEGER, INTENT( in )  ::   kt           ! ocean time-step index
(   85) 
(   86)       !! * Local declarations
(   87)       INTEGER  ::  ji, jj, jk, jit             ! dummy loop indices
(   88)       INTEGER  ::  icycle                      ! temporary scalar
(   89)       REAL(wp) ::                           &
(   90)          zraur, zcoef, z2dt_e, z2dt_b, zfac25,   &  ! temporary scalars
(   91)          zfact1, zspgu, zcubt, zx1, zy1,    &  !     "        "
(   92)          zfact2, zspgv, zcvbt, zx2, zy2        !     "        "
(   93)       REAL(wp), DIMENSION(jpi,jpj) ::       &
(   94)          zcu, zcv, zwx, zwy, zhdiv,         &  ! temporary arrays
(   95)          zua, zva, zub, zvb,                &  !     "        "
(   96)          zssha_b, zua_b, zva_b,             &  !     "        "
(   97)          zsshb_e, zub_e, zvb_e,             &  !     "        "
(   98)          zun_e, zvn_e                          !     "        "
(   99)       REAL(wp), DIMENSION(jpi,jpj),SAVE ::  &
(  100)          ztnw, ztne, ztsw, ztse
(  101)       !!----------------------------------------------------------------------
(  102) 
(  103)       ! Arrays initialization
(  104)       ! ---------------------






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 3

(  105)       zua_b(:,:) = 0.e0   ;   zub_e(:,:) = 0.e0   ;   zun_e(:,:) = 0.e0
(  106)       zva_b(:,:) = 0.e0   ;   zvb_e(:,:) = 0.e0   ;   zvn_e(:,:) = 0.e0
(  107)       zhdiv(:,:) = 0.e0
(  108) 
(  109)       IF( kt == nit000 ) THEN
(  110) 
(  111)          IF(lwp) WRITE(numout,*)
(  112)          IF(lwp) WRITE(numout,*) 'dyn_spg_ts_jki : surface pressure gradient trend'
(  113)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   free surface with time splitting with j-k-i loop'
(  114)          IF(lwp) WRITE(numout,*) ' Number of sub cycle in 1 time-step (2 rdt) : icycle = ', FLOOR( 2*rdt/rdtbt )
(  115) 
(  116)          IF( .NOT. ln_rstart ) THEN
(  117)             ! initialize barotropic specific arrays
(  118)             sshb_b(:,:) = sshb(:,:)
(  119)             sshn_b(:,:) = sshn(:,:)
(  120)             un_b(:,:)   = 0.e0
(  121)             vn_b(:,:)   = 0.e0
(  122)             ! vertical sum
(  123)             DO jk = 1, jpkm1
(  124)                un_b(:,:) = un_b(:,:) + fse3u(:,:,jk) * un(:,:,jk)
(  125)                vn_b(:,:) = vn_b(:,:) + fse3v(:,:,jk) * vn(:,:,jk)
(  126)             END DO
(  127)          ENDIF
(  128)          ssha_e(:,:) = sshn(:,:)
(  129)          ua_e  (:,:) = un_b(:,:)
(  130)          va_e  (:,:) = vn_b(:,:)
(  131) 
(  132)          IF( ln_dynvor_een ) THEN
(  133)             ztne(1,:) = 0.e0   ;   ztnw(1,:) = 0.e0   ;   ztse(1,:) = 0.e0   ;   ztsw(1,:) = 0.e0
(  134)             DO jj = 2, jpj
(  135)                DO ji = 2, jpi
(  136)                   ztne(ji,jj) = ( ff(ji-1,jj  ) + ff(ji  ,jj  ) + ff(ji  ,jj-1) ) / 3.
(  137)                   ztnw(ji,jj) = ( ff(ji-1,jj-1) + ff(ji-1,jj  ) + ff(ji  ,jj  ) ) / 3.
(  138)                   ztse(ji,jj) = ( ff(ji  ,jj  ) + ff(ji  ,jj-1) + ff(ji-1,jj-1) ) / 3.
(  139)                   ztsw(ji,jj) = ( ff(ji  ,jj-1) + ff(ji-1,jj-1) + ff(ji-1,jj  ) ) / 3.
(  140)                END DO
(  141)             END DO
(  142)          ENDIF
(  143) 
(  144)       ENDIF
(  145)     
(  146)       ! Local constant initialization
(  147)       ! --------------------------------
(  148)       z2dt_b = 2.0 * rdt                                    ! baroclinic time step
(  149)       IF ( neuler == 0 .AND. kt == nit000 ) z2dt_b = rdt
(  150)       zfact1 = 0.5 * 0.25                                   ! coefficient for vorticity estimates
(  151)       zfact2 = 0.5 * 0.5
(  152)       zraur  = 1. / rauw                                    ! 1 / volumic mass of pure water
(  153)       
(  154)       ! -----------------------------------------------------------------------------
(  155)       !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
(  156)       ! -----------------------------------------------------------------------------
(  157) 
(  158)       DO jj = 1, jpj
(  159) 
(  160)          ! variables for the barotropic equations
(  161)          zsshb_e(:,jj) = sshn_b(:,jj)       ! (barotropic) sea surface height (before and now)
(  162)          sshn_e (:,jj) = sshn_b(:,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 4

(  163)          zub_e  (:,jj) = un_b  (:,jj)       ! barotropic transports issued from the barotropic equations (before and now)
(  164)          zvb_e  (:,jj) = vn_b  (:,jj)
(  165)          zun_e  (:,jj) = un_b  (:,jj)
(  166)          zvn_e  (:,jj) = vn_b  (:,jj)
(  167)          zssha_b(:,jj) = sshn  (:,jj)        ! time averaged variables over all sub-timesteps
(  168)          zua_b  (:,jj) = un_b  (:,jj)   
(  169)          zva_b  (:,jj) = vn_b  (:,jj)
(  170) 
(  171)          ! Vertically integrated quantities
(  172)          ! --------------------------------
(  173)          zua(:,jj) = 0.e0
(  174)          zva(:,jj) = 0.e0
(  175)          zub(:,jj) = 0.e0
(  176)          zvb(:,jj) = 0.e0
(  177)          zwx(:,jj) = 0.e0
(  178)          zwy(:,jj) = 0.e0
(  179) 
(  180)          ! vertical sum
(  181)          DO jk = 1, jpkm1
(  182)             !                                                           ! Vertically integrated momentum trends
(  183)             zua(:,jj) = zua(:,jj) + fse3u(:,jj,jk) * umask(:,jj,jk) * ua(:,jj,jk)
(  184)             zva(:,jj) = zva(:,jj) + fse3v(:,jj,jk) * vmask(:,jj,jk) * va(:,jj,jk)
(  185)             !                                                           ! Vertically integrated transports (before)
(  186)             zub(:,jj) = zub(:,jj) + fse3u(:,jj,jk) * ub(:,jj,jk)
(  187)             zvb(:,jj) = zvb(:,jj) + fse3v(:,jj,jk) * vb(:,jj,jk)
(  188)             !                                                           ! Planetary vorticity (now)
(  189)             zwx(:,jj) = zwx(:,jj) + e2u(:,jj) * fse3u(:,jj,jk) * un(:,jj,jk)
(  190)             zwy(:,jj) = zwy(:,jj) + e1v(:,jj) * fse3v(:,jj,jk) * vn(:,jj,jk)
(  191)          END DO
(  192) 
(  193)       END DO
(  194) 
(  195)       DO jj = 2, jpjm1
(  196) 
(  197)          IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
(  198)             DO ji = 2, jpim1
(  199)                zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  200)                zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  201)                zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  202)                zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  203)                ! energy conserving formulation for planetary vorticity term
(  204)                zcu(ji,jj) = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
(  205)                zcv(ji,jj) =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
(  206)             END DO
(  207) 
(  208)          ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
(  209)             DO ji = 2, jpim1
(  210)                zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
(  211)                               + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  212)                zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
(  213)                               + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  214)                zcu(ji,jj)  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
(  215)                zcv(ji,jj)  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
(  216)             END DO
(  217) 
(  218)          ELSEIF ( ln_dynvor_een ) THEN                    ! enstrophy and energy conserving scheme
(  219)          zfac25 = 0.25
(  220)             DO ji = 2, jpim1






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 5

(  221)                zcu(ji,jj) = + zfac25 / e1u(ji,jj)   &
(  222)                   &       * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
(  223)                   &          + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
(  224)                zcv(ji,jj) = - zfac25 / e2v(ji,jj)   &
(  225)                   &       * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
(  226)                   &          + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
(  227)             END DO
(  228) 
(  229)          ENDIF
(  230) 
(  231) 
(  232)          ! Remove barotropic trend from general momentum trend
(  233)          DO jk = 1 , jpkm1
(  234)             DO ji = 2, jpim1
(  235)                ua(ji,jj,jk) = ua(ji,jj,jk) - zua(ji,jj) * hur(ji,jj)
(  236)                va(ji,jj,jk) = va(ji,jj,jk) - zva(ji,jj) * hvr(ji,jj)
(  237)             END DO
(  238)          END DO
(  239) 
(  240)          ! Remove coriolis term from barotropic trend
(  241)          ! ------------------------------------------
(  242)          DO ji = 2, jpim1
(  243)             zua(ji,jj) = zua(ji,jj) - zcu(ji,jj)
(  244)             zva(ji,jj) = zva(ji,jj) - zcv(ji,jj)
(  245)          END DO
(  246) 
(  247)       END DO
(  248) 
(  249)       ! -----------------------------------------------------------------------
(  250)       !  Phase 2 : Integration of the barotropic equations with time splitting
(  251)       ! -----------------------------------------------------------------------
(  252) 
(  253)       ! Initialisations
(  254)       !----------------
(  255)       ! Number of iteration of the barotropic loop
(  256)       icycle = FLOOR( z2dt_b / rdtbt )
(  257) 
(  258)       ! set ssh corrections to 0
(  259)       ! ssh corrections are applied to normal velocities (Flather's algorithm) and averaged over the barotropic loop
(  260) #if defined key_obc
(  261)       IF( lp_obc_east  )   sshfoe_b(:,:) = 0.e0
(  262)       IF( lp_obc_west  )   sshfow_b(:,:) = 0.e0
(  263)       IF( lp_obc_south )   sshfos_b(:,:) = 0.e0
(  264)       IF( lp_obc_north )   sshfon_b(:,:) = 0.e0
(  265) #endif
(  266) 
(  267)       ! Barotropic integration over 2 baroclinic time steps
(  268)       ! ---------------------------------------------------
(  269) 
(  270)       !                                                    ! ==================== !
(  271)       DO jit = 1, icycle                                   !  sub-time-step loop  !
(  272)          !                                                 ! ==================== !
(  273) 
(  274)          z2dt_e = 2. * rdtbt
(  275)          IF ( jit == 1 )   z2dt_e = rdtbt
(  276) 
(  277)          ! Time interpolation of open boundary condition data
(  278)          IF( lk_obc )   CALL obc_dta_bt( kt, jit )






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 6

(  279) 
(  280)          DO jj = 2, jpjm1
(  281) 
(  282)          ! Horizontal divergence of barotropic transports
(  283)          !--------------------------------------------------
(  284)             DO ji = 2, jpim1
(  285)                zhdiv(ji,jj) = ( e2u(ji  ,jj  ) * zun_e(ji  ,jj)              &
(  286)                   &            -e2u(ji-1,jj  ) * zun_e(ji-1,jj)              &
(  287)                   &            +e1v(ji  ,jj  ) * zvn_e(ji  ,jj)              &
(  288)                   &            -e1v(ji  ,jj-1) * zvn_e(ji  ,jj-1) )          &
(  289)                   &           / (e1t(ji,jj)*e2t(ji,jj))
(  290)             END DO
(  291) 
(  292) #if defined key_obc
(  293)          ! open boundaries (div must be zero behind the open boundary)
(  294)          !  mpp remark: The zeroing of zhdiv can probably be extended to 1->jpi/jpj for the correct row/column
(  295)          IF( lp_obc_east  ) THEN
(  296)             IF( nje0   <= jj .AND. jj <= nje1   )   zhdiv(nie0p1:nie1p1,jj) = 0.e0      ! east
(  297)          ENDIF
(  298)          IF( lp_obc_west  ) THEN
(  299)             IF( njw0   <= jj .AND. jj <= njw1   )   zhdiv(niw0  :niw1  ,jj) = 0.e0      ! west
(  300)          ENDIF
(  301)          IF( lp_obc_north ) THEN
(  302)             IF( njn0p1 <= jj .AND. jj <= njn1p1 )   zhdiv(nin0  :nin1  ,jj) = 0.e0      ! north
(  303)          ENDIF
(  304)          IF( lp_obc_south ) THEN
(  305)             IF( njs0   <= jj .AND. jj <= njs1   )   zhdiv(nis0  :nis1  ,jj) = 0.e0      ! south
(  306)          ENDIF
(  307) #endif
(  308) 
(  309)          ! Sea surface height from the barotropic system
(  310)          !----------------------------------------------
(  311)             DO ji = 2, jpim1
(  312)                ssha_e(ji,jj) = ( zsshb_e(ji,jj) - z2dt_e *  ( zraur * emp(ji,jj)  &
(  313)                   &            +  zhdiv(ji,jj) ) ) * tmask(ji,jj,1)
(  314)             END DO
(  315) 
(  316)          END DO
(  317) 
(  318)          ! evolution of the barotropic transport ( following the vorticity scheme used)
(  319)          ! ----------------------------------------------------------------------------
(  320)          zwx(:,:) = e2u(:,:) * zun_e(:,:)
(  321)          zwy(:,:) = e1v(:,:) * zvn_e(:,:)
(  322) 
(  323)          DO jj = 2, jpjm1
(  324) 
(  325)             IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
(  326)                DO ji = 2, jpim1
(  327)                   ! surface pressure gradient
(  328)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  329)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  330)                   ! energy conserving formulation for planetary vorticity term
(  331)                   zy1 = ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
(  332)                   zy2 = ( zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  333)                   zx1 = ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
(  334)                   zx2 = ( zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  335)                   zcubt = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
(  336)                   zcvbt =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 7

(  337)                   ! after transports
(  338)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  339)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  340)                END DO
(  341) 
(  342)             ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
(  343)                DO ji = 2, jpim1
(  344)                   ! surface pressure gradient
(  345)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  346)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  347)                   ! enstrophy conserving formulation for planetary vorticity term
(  348)                   zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
(  349)                                  + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
(  350)                   zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
(  351)                                  + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
(  352)                   zcubt  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
(  353)                   zcvbt  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
(  354)                   ! after transports
(  355)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  356)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  357)                END DO
(  358) 
(  359)             ELSEIF ( ln_dynvor_een ) THEN                    ! energy and enstrophy conserving scheme
(  360)                zfac25 = 0.25
(  361)                DO ji = 2, jpim1
(  362)                   ! surface pressure gradient
(  363)                   zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
(  364)                   zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
(  365)                   ! energy/enstrophy conserving formulation for planetary vorticity term
(  366)                   zcubt = + zfac25 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
(  367)                      &                             + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
(  368)                   zcvbt = - zfac25 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
(  369)                      &                             + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
(  370)                   ! after transports
(  371)                   ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
(  372)                   va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
(  373)                END DO
(  374)             ENDIF
(  375) 
(  376)          END DO
(  377) 
(  378) 
(  379)          ! Flather's boundary condition for the barotropic loop :
(  380)          !         - Update sea surface height on each open boundary
(  381)          !         - Correct the barotropic transports
(  382)          IF( lk_obc )   CALL obc_fla_ts
(  383) 
(  384)          ! ... Boundary conditions on ua_e, va_e, ssha_e
(  385)          CALL lbc_lnk( ua_e  , 'U', -1. )
(  386)          CALL lbc_lnk( va_e  , 'V', -1. )
(  387)          CALL lbc_lnk( ssha_e, 'T',  1. )
(  388) 
(  389)          DO jj = 1, jpj
(  390) 
(  391)             ! temporal sum
(  392)             !-------------
(  393)             zssha_b(:,jj) = zssha_b(:,jj) + ssha_e(:,jj)
(  394)             zua_b  (:,jj) = zua_b  (:,jj) + ua_e  (:,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 8

(  395)             zva_b  (:,jj) = zva_b  (:,jj) + va_e  (:,jj) 
(  396) 
(  397)             ! Time filter and swap of dynamics arrays
(  398)             ! ---------------------------------------
(  399)             IF( neuler == 0 .AND. kt == nit000 ) THEN   ! Euler (forward) time stepping
(  400)                zsshb_e(:,jj) = sshn_e(:,jj)
(  401)                zub_e  (:,jj) = zun_e (:,jj)
(  402)                zvb_e  (:,jj) = zvn_e (:,jj)
(  403)                sshn_e (:,jj) = ssha_e(:,jj)
(  404)                zun_e  (:,jj) = ua_e  (:,jj)
(  405)                zvn_e  (:,jj) = va_e  (:,jj)
(  406)             ELSE                                        ! Asselin filtering
(  407)                zsshb_e(:,jj) = atfp * ( zsshb_e(:,jj) + ssha_e(:,jj) ) + atfp1 * sshn_e (:,jj)
(  408)                zub_e  (:,jj) = atfp * ( zub_e  (:,jj) + ua_e  (:,jj) ) + atfp1 * zun_e  (:,jj)
(  409)                zvb_e  (:,jj) = atfp * ( zvb_e  (:,jj) + va_e  (:,jj) ) + atfp1 * zvn_e  (:,jj)
(  410)                sshn_e (:,jj) = ssha_e(:,jj)
(  411)                zun_e  (:,jj) = ua_e  (:,jj)
(  412)                zvn_e  (:,jj) = va_e  (:,jj)
(  413)             ENDIF
(  414) 
(  415)          END DO
(  416) 
(  417)          !                                                 ! ==================== !
(  418)       END DO                                               !        end loop      !
(  419)       !                                                    ! ==================== !
(  420) 
(  421) 
(  422)       ! Time average of after barotropic variables
(  423)       zcoef =  1.e0 / (  FLOAT( icycle +1 )  )
(  424)       zssha_b(:,:) = zcoef * zssha_b(:,:) 
(  425)       zua_b  (:,:) = zcoef * zua_b  (:,:) 
(  426)       zva_b  (:,:) = zcoef * zva_b  (:,:) 
(  427) #if defined key_obc
(  428)          IF( lp_obc_east  )   sshfoe_b(:,:) = zcoef * sshfoe_b(:,:)
(  429)          IF( lp_obc_west  )   sshfow_b(:,:) = zcoef * sshfow_b(:,:)
(  430)          IF( lp_obc_north )   sshfon_b(:,:) = zcoef * sshfon_b(:,:)
(  431)          IF( lp_obc_south )   sshfos_b(:,:) = zcoef * sshfos_b(:,:)
(  432) #endif
(  433)      
(  434) 
(  435)       ! ---------------------------------------------------------------------------
(  436)       ! Phase 3 : Update sea surface height from time averaged barotropic variables
(  437)       ! ---------------------------------------------------------------------------
(  438) 
(  439)       sshb(:,:) = sshn(:,:)
(  440)  
(  441)       ! Horizontal divergence of time averaged barotropic transports
(  442)       !-------------------------------------------------------------
(  443)       DO jj = 2, jpjm1
(  444)          DO ji = 2, jpim1
(  445)             zhdiv(ji,jj) = ( e2u(ji,jj) * un_b(ji,jj) - e2u(ji-1,jj  ) * un_b(ji-1,jj  )     &
(  446)            &                +e1v(ji,jj) * vn_b(ji,jj) - e1v(ji  ,jj-1) * vn_b(ji  ,jj-1) )   &
(  447)            &             / ( e1t(ji,jj) * e2t(ji,jj) )
(  448)          END DO
(  449) 
(  450) #if defined key_obc
(  451)          ! open boundaries (div must be zero behind the open boundary)
(  452)          !  mpp remark: The zeroing of zhdiv can probably be extended to 1->jpi/jpj for the correct row/column






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 9

(  453)          IF( lp_obc_east  ) THEN
(  454)             IF( nje0   <= jj .AND. jj <= nje1   )   zhdiv(nie0p1:nie1p1,jj) = 0.e0      ! east
(  455)          ENDIF
(  456)          IF( lp_obc_west  ) THEN
(  457)             IF( njw0   <= jj .AND. jj <= njw1   )   zhdiv(niw0  :niw1  ,jj) = 0.e0      ! west
(  458)          ENDIF
(  459)          IF( lp_obc_north ) THEN
(  460)             IF( njn0p1 <= jj .AND. jj <= njn1p1 )   zhdiv(nin0  :nin1  ,jj) = 0.e0      ! north
(  461)          ENDIF
(  462)          IF( lp_obc_south ) THEN
(  463)             IF( njs0   <= jj .AND. jj <= njs1   )   zhdiv(nis0  :nis1  ,jj) = 0.e0      ! south
(  464)          ENDIF
(  465) #endif
(  466) 
(  467)          ! sea surface height
(  468)          !-------------------
(  469)          DO ji = 2, jpim1
(  470)             sshn(ji,jj) = (  sshb_b(ji,jj) - z2dt_b * ( zraur * emp(ji,jj) + zhdiv(ji,jj) )  ) * tmask(ji,jj,1)
(  471)          END DO
(  472) 
(  473)       END DO
(  474) 
(  475)       ! ... Boundary conditions on sshn
(  476)       IF( .NOT. lk_obc ) CALL lbc_lnk( sshn, 'T', 1. )
(  477) 
(  478) 
(  479)       ! -----------------------------------------------------------------------------
(  480)       ! Phase 4. Coupling between general trend and barotropic estimates - (2nd step)
(  481)       ! -----------------------------------------------------------------------------
(  482) 
(  483)       DO jj = 1, jpj
(  484) 
(  485)          ! Swap on time averaged barotropic variables
(  486)          ! ------------------------------------------
(  487)          sshb_b(:,jj) = sshn_b (:,jj)
(  488)          sshn_b(:,jj) = zssha_b(:,jj)
(  489)          un_b  (:,jj) = zua_b  (:,jj) 
(  490)          vn_b  (:,jj) = zva_b  (:,jj) 
(  491)    
(  492)          ! add time averaged barotropic coriolis and surface pressure gradient
(  493)          ! terms to the general momentum trend
(  494)          ! --------------------------------------------------------------------
(  495)          DO jk = 1, jpkm1
(  496)             ua(:,jj,jk) = ua(:,jj,jk) + hur(:,jj) * ( zua_b(:,jj) - zub(:,jj) ) / z2dt_b
(  497)             va(:,jj,jk) = va(:,jj,jk) + hvr(:,jj) * ( zva_b(:,jj) - zvb(:,jj) ) / z2dt_b
(  498)          END DO
(  499)       END DO
(  500) 
(  501)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  502)          CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
(  503)       ENDIF
(  504)       
(  505)    END SUBROUTINE dyn_spg_ts_jki
(  506) #else
(  507)    !!----------------------------------------------------------------------
(  508)    !!   Default case :   Empty module   No standart free surface cst volume
(  509)    !!----------------------------------------------------------------------
(  510)    USE in_out_manager






PGF90 (Version     10.5)          09/23/2010  09:49:43      page 10

(  511) CONTAINS
(  512)    SUBROUTINE dyn_spg_ts_jki( kt )       ! Empty routine
(  513)       if(lwp) WRITE(numout,*) 'dyn_spg_ts_jki: You should not have seen this print! error?', kt
(  514)    END SUBROUTINE dyn_spg_ts_jki
(  515) #endif
(  516)    
(  517)    !!======================================================================
(  518) END MODULE dynspg_ts_jki
























































PGF90 (Version     10.5)          09/23/2010  09:49:43      page 11

(    1) # 1 "dynspg_ts_jki.F90"
(    1) MODULE dynspg_ts_jki
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_ts_jki  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) # 507
(  507)    !!----------------------------------------------------------------------
(  508)    !!   Default case :   Empty module   No standart free surface cst volume
(  509)    !!----------------------------------------------------------------------
(  510)    USE in_out_manager
(  511) CONTAINS
(  512)    SUBROUTINE dyn_spg_ts_jki( kt )       ! Empty routine
(  513)       if(lwp) WRITE(numout,*) 'dyn_spg_ts_jki: You should not have seen this print! error?', kt
(  514)    END SUBROUTINE dyn_spg_ts_jki
(  515) # 516
(  516)    
(  517)    !!======================================================================
(  518) END MODULE dynspg_ts_jki
