


PGF90 (Version     10.5)          09/23/2010  09:50:09      page 1

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

Filename: sopa_mc.F90

(    1) !!DB 2009.08
(    2) !!Module that contains various (diagnostic) routines specific to the
(    3) !!Maritime Canada Shelf OPA model
(    4) !!Quantities computed once per day
(    5) 
(    6) !!Started as flxfwb; 
(    7) !! kept many of the variables although they may not be used
(    8) 
(    9) MODULE sopa_mc
(   10) 
(   11)   !! * Modules used
(   12)   USE oce             ! ocean dynamics and tracers
(   13)   USE dom_oce         ! ocean space and time domain
(   14)   USE cpl_oce         ! coupled atmosphere/ocean
(   15)   USE phycst          ! physical constants
(   16)   USE in_out_manager  ! I/O manager
(   17)   USE lib_mpp         ! distribued memory computing library
(   18)   USE flxrnf          ! ocean runoffs
(   19)   USE ocesbc          ! ocean surface boudaries conditions
(   20)   USE blk_oce
(   21)   USE flxblk          ! bulk formulea
(   22)   USE daymod          ! calendar
(   23) #if defined key_obc
(   24)   USE obc_oce
(   25) #endif
(   26) #if defined key_ice_lim
(   27)   USE ice
(   28) #endif
(   29) 
(   30)   IMPLICIT NONE
(   31)   PRIVATE
(   32) 
(   33)   !! * Routine accessibility
(   34)   PUBLIC sopa_mc_diagnostics      ! routine called by step
(   35) 
(   36)   !! * Share module variables
(   37)   REAL(wp), PUBLIC ::   &  !:
(   38)        vol_tot, sfce_area,  &
(   39)        glb_emp   ,  & ! domain averaged evaporation minus precipitation
(   40)        glb_precip,  & ! domain averaged precipitation
(   41)        glb_rnf   ,  & ! domain averaged runoff
(   42)        glb_ssh  ,  & ! domain averaged sea surface height
(   43)        glb_S   ,  & ! domain averaged ocean salinity
(   44)        glb_T        ! domain averaged ocean temperature
(   45) 
(   46)   REAL(wp), DIMENSION(jpi,jpj) ::  &






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 2

(   47)        e1e2_i    ! area of the interior domain (e1t*e2t*tmask_i)
(   48) 
(   49)   !! * Substitutions
(   50) #  include "domzgr_substitute.h90"
(   51) #  include "vectopt_loop_substitute.h90"
(   52) 
(   53) CONTAINS
(   54) 
(   55) 
(   56) !!DB: Once-per-day, compute globally-averaged (ie domain-averaged):
(   57) !!    T, S, emp, ssh, rnf, OBC TS flux, ice_area, ice_vol,  ...
(   58) !! Every dt: Output a check velocity from the Bay-of-Fundy
(   59) !! NB: Volume calc is done every call in case free sfce ever becomes part
(   60) !!     of fse3t(i,j,1). 
(   61) !!     Also cell_volume(i,j,k) could be computed at nit000 and saved but 
(   62) !!     this is not done as calc occurs only once-per-day
(   63) 
(   64) 
(   65)   SUBROUTINE sopa_mc_diagnostics( kt )
(   66)     !!---------------------------------------------------------------------
(   67)     !!
(   68)     !! ** Purpose : Compute some domain-averaged quantities 
(   69)     !!              Output to text file(s)
(   70)     !!
(   71)     !!	
(   72) 
(   73)     !!----------------------------------------------------------------------
(   74)     !! * Arguments
(   75)     INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   76) 
(   77)     !! * Local declarations
(   78)     INTEGER  ::   ji, jj, jk
(   79)     Logical ::  BoF_output
(   80)     CHARACTER*175 :: header
(   81)     REAL(wp) ::   zwei, zsrau, & 
(   82)          uT_east, uT_west, vT_south, vT_north, &
(   83)          uS_east, uS_west, vS_south, vS_north, &
(   84)          dT_sfce, dS_sfce, net_T_flux, net_S_flux, &
(   85)          ice_area, ice_vol, ave_T_sfce, ave_S_sfce, vol_sfce
(   86) 
(   87)     !!----------------------------------------------------------------------
(   88) 
(   89) 
(   90) !!DB: The dT & dS_sfce calcs mimic the ones found in trasbc.F90, viz.
(   91) !            zta = ro0cpr * ( qt(ji,jj) - qsr(ji,jj) ) * zse3t
(   92) !            zsa = emps(ji,jj) * zsrau * sn(ji,jj,1) * zse3t
(   93) 
(   94) 
(   95)     BoF_output = .true.    ! false ===> do not output the test vel in the BoF
(   96) 
(   97)     if( kt == nit000 ) then
(   98)        IF(lwp) THEN
(   99)           WRITE(numout2,*)
(  100)           WRITE(numout2,*) 'DB: SOPA_MC_DIAGNOSTICS'
(  101)           WRITE(numout2,*) '~~~~~~~'
(  102)        ENDIF
(  103) 
(  104)        ice_area = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 3

(  105)        ice_vol = 0.e0
(  106)        dT_sfce = 0.e0    !K/s
(  107)        dS_sfce = 0.e0    !S/s
(  108)        glb_emp    = 0.e0
(  109)        glb_precip = 0.e0
(  110)        glb_rnf    = 0.e0
(  111)        glb_ssh   = 0.e0   ! averaged sea surface height 
(  112)        glb_S = 0.e0   ! averaged ocean salinity 
(  113)        glb_T = 0.e0   ! averaged ocean temperature
(  114) !!DB OBC variables; assigned even if key_obc = off 
(  115)        uT_east = 0.e0
(  116)        uT_west = 0.e0
(  117)        vT_south = 0.e0
(  118)        vT_north = 0.e0
(  119)        uS_east = 0.e0
(  120)        uS_west = 0.e0
(  121)        vS_south = 0.e0
(  122)        vS_north = 0.e0
(  123) 
(  124) 
(  125) 
(  126) !! da(i,j) & surface area
(  127)        e1e2_i(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:)
(  128)        sfce_area    = SUM(e1e2_i(:,:))
(  129)        IF( lk_mpp )   CALL  mpp_sum( sfce_area    )   ! sum over the global domain
(  130) 
(  131)        if(BoF_output) then
(  132) !!       isolate the correct processor to open the file
(  133)           if(mi0(58)==mi1(58) .AND. mj0(106)== mj1(106)) then
(  134)              open(999,file='BoF_vel.dat')
(  135)           endif
(  136)        endif
(  137) 
(  138)        if(lwp) then
(  139)           open(998,file='sopa_mc_diagnostics.dat')
(  140)           header = 'year day kt  [ave_T ave_S ave_SSH] [OBC(VT(E,W,S) VS(E,W,S))] &
(  141)                dT_sfce dS_sfce net_T_flux net_S_flux (TS/sec) [ave(emp,precip,rnf)] ice_area ice_vol &
(  142)                [ave_T ave_S k=1,2]'
(  143)           write(998,'(A175)')header
(  144)        endif
(  145) 
(  146)     endif  !!kt = nit000 loop
(  147) 
(  148)     if(BoF_output) then
(  149) !!DBG: Output a velocity from the Bay of Fundy 
(  150) !!isolate processor that owns this global location
(  151)        ji = 58; jj = 106
(  152)        if(mi0(ji)==mi1(ji) .AND. mj0(jj)== mj1(jj)) then
(  153) !!OLD+NEW
(  154)           write(999,'(i7,2x,12(f8.3,1x))')kt, &
(  155)                (un(ji,jj,1)+un(ji,jj,2))/(fse3t(ji,jj,1)+fse3t(ji,jj,2)), &
(  156)                (vn(ji,jj,1)+vn(ji,jj,2))/(fse3t(ji,jj,1)+fse3t(ji,jj,2))
(  157) 
(  158) !          un(ji,jj,1), vn(ji,jj,1),un(ji,jj,3),vn(ji,jj,3)
(  159) !!NEW: decimal year, kt, Vel(i,j,1:2)/dz(1:2) ... in progress ...
(  160) !             write(999,'(i7,2x,12(f8.3,1x))')nyear, nday_year, 
(  161) !                  (un(ji,jj,1)+un(ji,jj,2))/(fse3t(ji,jj,1)+fse3t(ji,jj,2)), &
(  162) !                  (vn(ji,jj,1)+vn(ji,jj,2))/(fse3t(ji,jj,1)+fse3t(ji,jj,2))






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 4

(  163)           endif
(  164)        !!DB:END
(  165)     endif
(  166) 
(  167) 
(  168) !!DB: Do once-per-day -- Hardwired
(  169)       if(mod(kt-nit000,int(rday/rdt)) ==0) then
(  170) 
(  171)          zsrau = 1. / rauw
(  172)          dT_sfce = sum(e1e2_i(:,:)*qt(:,:)*ro0cpr)  ! >0 ===> dT > 0
(  173)          if(lk_mpp) call mpp_sum(dT_sfce)
(  174) !         dS_sfce = sum(e1e2_i(:,:)*emps(:,:)*sn(:,:,1)*zsrau)  ! >0 ===> dS > 0
(  175) !!Later must divide by glb_S        not           |||
(  176)          dS_sfce = sum(e1e2_i(:,:)*emps(:,:)*zsrau)  ! >0 ===> dS > 0
(  177)          if(lk_mpp) call mpp_sum(dS_sfce)
(  178) 
(  179)          glb_emp    = SUM( e1e2_i(:,:) * emp   (:,:) )
(  180)          IF( lk_mpp )   CALL  mpp_sum( glb_emp    )   ! sum over the global domain
(  181) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  182)          glb_precip = SUM( e1e2_i(:,:) * watm  (:,:) )
(  183)          IF( lk_mpp )   CALL  mpp_sum( glb_precip )   ! sum over the global domain
(  184) #endif
(  185)          glb_rnf    = SUM( e1e2_i(:,:) * runoff(:,:) )
(  186)          IF( lk_mpp )   CALL  mpp_sum( glb_rnf    )   ! sum over the global domain
(  187)          
(  188)          glb_ssh = SUM( e1e2_i(:,:) * sshn(:,:) )
(  189)          if( lk_mpp ) call mpp_sum(glb_ssh)
(  190) 
(  191) #if defined key_ice_lim
(  192)          ice_area = SUM(e1e2_i(:,:)*(1.0-frld(:,:)))  !m2
(  193)          ice_vol = SUM(e1e2_i(:,:)*(1.0-frld(:,:))*hicif(:,:)) !m3
(  194)          if( lk_mpp ) then 
(  195)             call mpp_sum(ice_area)
(  196)             call mpp_sum(ice_vol)
(  197)          endif
(  198) #endif
(  199)          vol_tot  = 0.e0
(  200)          glb_S = 0.0
(  201)          glb_T = 0.0
(  202)          do jk = 1, jpkm1   
(  203)             do jj = 2, jpjm1
(  204)                do ji = fs_2, fs_jpim1   ! vector opt.
(  205)                   zwei  = e1e2_i(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)
(  206)                   glb_S = glb_S + ( sn(ji,jj,jk) ) * zwei
(  207)                   glb_T = glb_T + ( tn(ji,jj,jk) ) * zwei
(  208)                   vol_tot  = vol_tot  + zwei
(  209)                enddo
(  210)             enddo
(  211)          enddo
(  212)          vol_sfce = 0.e0
(  213)          ave_T_sfce = 0.e0
(  214)          ave_S_sfce = 0.e0
(  215)          do jk = 1, 2
(  216)             do jj = 2, jpjm1
(  217)                do ji = fs_2, fs_jpim1   ! vector opt.
(  218)                   zwei  = e1e2_i(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)
(  219)                   ave_S_sfce = ave_S_sfce + ( sn(ji,jj,jk) ) * zwei
(  220)                   ave_T_sfce = ave_T_sfce + ( tn(ji,jj,jk) ) * zwei






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 5

(  221)                   vol_sfce  = vol_sfce + zwei
(  222)                enddo
(  223)             enddo
(  224)          enddo
(  225)          IF( lk_mpp ) THEN
(  226)             call  mpp_sum( glb_S )        
(  227)             call  mpp_sum( glb_T )        
(  228)             call  mpp_sum( vol_tot )        
(  229)             call  mpp_sum( ave_S_sfce )        
(  230)             call  mpp_sum( ave_T_sfce )        
(  231)             call  mpp_sum( vol_sfce )        
(  232)          ENDIF
(  233) 
(  234) #ifdef key_obc
(  235) !!DB -- OB T & S flux 
(  236) !!NB: because this routine is called after fields are updated, the
(  237) !!    now variables _should_ contain the correct values
(  238) 
(  239)          if( lp_obc_east  ) then
(  240)             uT_east = 0.e0;  uS_east = 0.e0
(  241)             do ji = nie0, nie1
(  242)                do jj = 1, jpj
(  243)                   do jk = 1, jpk
(  244)                      uT_east = uT_east +  & 
(  245)                           un(ji,jj,jk)*tn(ji+1,jj,jk)*e2u(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
(  246)                      uS_east = uS_east +  & 
(  247)                           un(ji,jj,jk)*sn(ji+1,jj,jk)*e2u(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
(  248)                   enddo
(  249)                enddo
(  250)             enddo
(  251)             if( lk_mpp ) then
(  252)                call mpp_sum(uT_east)
(  253)                call mpp_sum(uS_east)
(  254)             endif
(  255) 
(  256)          endif
(  257)          if( lp_obc_west  ) then
(  258)             uT_west = 0.e0;  uS_west = 0.e0
(  259)             do ji = niw0, niw1
(  260)                do jj = 1, jpj
(  261)                   do jk = 1, jpk
(  262)                      uT_west = uT_west + &
(  263)                           un(ji,jj,jk)*tn(ji,jj,jk)*e2u(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
(  264)                      uS_west = uS_west + &
(  265)                           un(ji,jj,jk)*sn(ji,jj,jk)*e2u(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
(  266)                   enddo
(  267)                enddo
(  268)             enddo
(  269)             if( lk_mpp ) then
(  270)                call mpp_sum(uT_west)
(  271)                call mpp_sum(uS_west)
(  272)             endif
(  273)          endif
(  274) 
(  275)          if( lp_obc_south ) then
(  276)             vT_south = 0.e0;  vS_south = 0.e0
(  277)             do jj = njs0, njs1
(  278)                do ji = 1, jpi






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 6

(  279)                   do jk = 1, jpk
(  280)                      vT_south = vT_south + &
(  281)                           vn(ji,jj,jk)*tn(ji,jj,jk)*e1v(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
(  282)                      vS_south = vS_south + &
(  283)                           vn(ji,jj,jk)*sn(ji,jj,jk)*e1v(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
(  284)                   enddo
(  285)                enddo
(  286)             enddo
(  287)             if( lk_mpp ) then
(  288)                call mpp_sum(vT_south)
(  289)                call mpp_sum(vS_south)
(  290)             endif
(  291)          endif
(  292) 
(  293) !!DB -- leave blank as SOPA-MC north is closed (lazy)
(  294) !         if( lp_obc_north ) then
(  295)          !write code later
(  296) !           vT_north = 0.e0
(  297) !         endif
(  298) 
(  299) #endif   !!key_obc
(  300) 
(  301) 
(  302) !! OUTPUT: average quantities over vol or area
(  303) !! NB: OBC fluxes are divided-by volume ---> units deg/s 
(  304) !! NB: leave emp, precip, runoff in whatever units they exist in 
(  305)          if(lwp) then
(  306)             glb_T = glb_T/vol_tot
(  307)             glb_S = glb_S/vol_tot
(  308)             glb_ssh = glb_ssh/sfce_area
(  309)             net_T_flux = -(uT_east/vol_tot-uT_west/vol_tot)  &
(  310)                  -(0.e0 - vT_south/vol_tot) + dT_sfce/vol_tot
(  311)             net_S_flux = -(uS_east/vol_tot-uS_west/vol_tot)  &
(  312)                  -(0.e0 - vS_south/vol_tot) + dS_sfce/(vol_tot*glb_S)
(  313)             write(998,'(2(i4,1x),i9,1x,3(f10.4,1x),19(e10.4,1x))')     &
(  314)                  nyear, nday_year, kt, & 
(  315)                  glb_T, glb_S, glb_ssh, &
(  316)                  uT_east/vol_tot,uT_west/vol_tot, vT_south/vol_tot, &
(  317)                  uS_east/vol_tot,uS_west/vol_tot, vS_south/vol_tot, &
(  318)                  dT_sfce/vol_tot, dS_sfce/(vol_tot*glb_S), &
(  319)                  net_T_flux, net_S_flux, &
(  320)                  glb_emp/sfce_area, glb_precip/sfce_area, glb_rnf/sfce_area, &
(  321)                  ice_area/1.e6, ice_vol/1.e9, &
(  322)                  ave_T_sfce/vol_sfce, ave_S_sfce/vol_sfce
(  323)          endif
(  324) 
(  325) !AD:
(  326)          call vertical_mixing_diagnostics(kt)
(  327) 
(  328)       endif     !! END (kt-nit000) % rday/rdt
(  329) 
(  330) 
(  331) 
(  332) !!Old code, but keep around for awhile
(  333)          ! Conversion in m3
(  334) !         glb_emp    = glb_emp    * rdttra(1) * 1.e-3 
(  335) !         glb_precip = glb_precip * rdttra(1) * 1.e-3 / rday
(  336) !         glb_rnf    = glb_rnf    * rdttra(1) * 1.e-3






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 7

(  337)     ! Ecriture des diagnostiques 
(  338)     ! --------------------------
(  339) 
(  340) !    IF( kt == nitend ) THEN
(  341) 
(  342) !       WRITE(inum,*)    'Net freshwater budget '
(  343) !       WRITE(inum,9010) '  emp    = ', a_emp   , ' m3 =', a_emp   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  344) !       WRITE(inum,9010) '  precip = ', a_precip, ' m3 =', a_precip/((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  345) !       WRITE(inum,9010) '  a_rnf  = ', a_rnf   , ' m3 =', a_rnf   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  346) !       WRITE(inum,*)    'Mean sea level : '
(  347) !       WRITE(inum,9010) '  at nit000 = ',a_sshb           ,' m3 '
(  348) !       WRITE(inum,9010) '  at nitend = ',a_sshend         ,' m3 '
(  349) !       WRITE(inum,9010) '  diff      = ',(a_sshend-a_sshb),' m3 =', (a_sshend-a_sshb)/((nitend-nit000+1)*rdt) * 1.e-6,' Sv'
(  350) !       WRITE(inum,9020) '  mean sea level elevation    =', a_sshend/zarea,' m'
(  351) !    ENDIF
(  352) 
(  353) 
(  354)   END SUBROUTINE sopa_mc_diagnostics
(  355) 
(  356) 
(  357)   SUBROUTINE vertical_mixing_diagnostics(kt)
(  358) 
(  359)     !!---------------------------------------------------------------------
(  360)     !!
(  361)     !! ** Purpose : Diagnostic output of various vertical mixing parameters
(  362)     !!              (Currently just the vertical viscosity coeff)
(  363)     !!  
(  364)     !! ** Method
(  365)     !! 
(  366)     !!   domain which contains the global OUTi, OUTj will compute the 
(  367)     !!   vertical maxval. Other domains remain at     maxV = 0.
(  368)     !!   call mpp_max(maxV) ensures all domains are updated. Only lwp 
(  369)     !!   writes output. For the global maximum, the generic mpp_maxloc is used.
(  370)     !!   It always returns a value from boundary. Might need to modify it to
(  371)     !!   look at internal domain only. 
(  372)     !!     
(  373)     !!	
(  374)     !!----------------------------------------------------------------------
(  375) 
(  376)     USE zdf_oce
(  377) 
(  378)     INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  379)     
(  380)     !! * Local declarations
(  381)     INTEGER  ::  jj,ji,fid
(  382)     REAL(wp) :: maxV,maxV_Tot
(  383)     INTEGER  :: MAXi,MAXj,MAXk
(  384)     INTEGER  :: OUTi_local, OUTj_local,maxV_local
(  385)     INTEGER  :: OUTi, OUTj
(  386)     parameter (OUTi=67,OUTj=63) !global location for desired output
(  387) 
(  388)     !!----------------------------------------------------------------------
(  389) 
(  390)     fid=numout2
(  391)     
(  392)     
(  393)     maxV = 0
(  394)     OUTj_local=OUTj + 1 -  njmpp






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 8

(  395)     OUTi_local=OUTi + 1 -  nimpp    
(  396)     IF  ( OUTi_local .le. jpi .AND. OUTi_local .ge. 1 .AND. &
(  397)             OUTj_local .le. jpj .AND. OUTj_local .ge. 1) then
(  398) 	  maxV = MAXVAL( avmu(OUTi_local,OUTj_local,:), & 
(  399) 	  mask = umask(OUTi_local,OUTj_local,:) == 1.e0)
(  400) !!DBG
(  401) !          write(3001,*) kt,avmu(OUTi_local,OUTj_local,1:mbathy(OUTi_local,OUTj_local))
(  402)     ENDIF
(  403) 
(  404)     IF(lk_mpp) call mpp_max(maxV)
(  405)     CALL mpp_maxloc( avmu,tmask,maxV_Tot,MAXi,MAXj,MAXk )
(  406) 
(  407)     IF(lwp) THEN 
(  408)        write(fid,'( "VERTICAL_MIXING_DIAGNOSTICS:"  )' )
(  409)        write(fid,*) '     @ kt,OUTi,OUTj,OUTi_local,OUTj_local', &
(  410)                        kt,OUTi,OUTj,OUTi_local,OUTj_local
(  411)        write(fid,*) '   glamu,gphiu = ', glamu(OUTi_local,OUTj_local) &
(  412)                                        , gphiu(OUTi_local,OUTj_local)  
(  413) ! AD. NOTE: glamu,gphiu not working!!!
(  414) !		       
(  415) 		       
(  416)        write(fid,'( "     vertical viscosity coeff. at uw-point"  )')
(  417)        write(fid,'( "MAXVAL avmu at: ",I4,I4," = ",E10.4,", &
(  418)                     Overall MAXVAL at: " ,I4,I4,I4," = ",E10.4)')  &  
(  419)               	    OUTi, OUTj,maxV,MAXi,MAXj,MAXk,maxV_Tot           
(  420)     ENDIF
(  421) 
(  422)   END SUBROUTINE vertical_mixing_diagnostics
(  423) 
(  424) 
(  425) END MODULE sopa_mc

































PGF90 (Version     10.5)          09/23/2010  09:50:09      page 9

(    1) # 1 "sopa_mc.F90"
(    1) !!DB 2009.08
(    2) !!Module that contains various (diagnostic) routines specific to the
(    3) !!Maritime Canada Shelf OPA model
(    4) !!Quantities computed once per day
(    5) 
(    6) !!Started as flxfwb; 
(    7) !! kept many of the variables although they may not be used
(    8) 
(    9) MODULE sopa_mc
(   10) 
(   11)   !! * Modules used
(   12)   USE oce             ! ocean dynamics and tracers
(   13)   USE dom_oce         ! ocean space and time domain
(   14)   USE cpl_oce         ! coupled atmosphere/ocean
(   15)   USE phycst          ! physical constants
(   16)   USE in_out_manager  ! I/O manager
(   17)   USE lib_mpp         ! distribued memory computing library
(   18)   USE flxrnf          ! ocean runoffs
(   19)   USE ocesbc          ! ocean surface boudaries conditions
(   20)   USE blk_oce
(   21)   USE flxblk          ! bulk formulea
(   22)   USE daymod          ! calendar
(   23) # 24
(   24)   USE obc_oce
(   25) # 27
(   27)   USE ice
(   28) 
(   29) # 30
(   30)   IMPLICIT NONE
(   31)   PRIVATE
(   32) 
(   33)   !! * Routine accessibility
(   34)   PUBLIC sopa_mc_diagnostics      ! routine called by step
(   35) 
(   36)   !! * Share module variables
(   37)   REAL(wp), PUBLIC ::   &  !:
(   38)        vol_tot, sfce_area,  &
(   39)        glb_emp   ,  & ! domain averaged evaporation minus precipitation
(   40)        glb_precip,  & ! domain averaged precipitation
(   41)        glb_rnf   ,  & ! domain averaged runoff
(   42)        glb_ssh  ,  & ! domain averaged sea surface height
(   43)        glb_S   ,  & ! domain averaged ocean salinity
(   44)        glb_T        ! domain averaged ocean temperature
(   45) 
(   46)   REAL(wp), DIMENSION(jpi,jpj) ::  &
(   47)        e1e2_i    ! area of the interior domain (e1t*e2t*tmask_i)
(   48) 
(   49)   !! * Substitutions
(   50) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 10

(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 51 "sopa_mc.F90"
(   51) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 52 "sopa_mc.F90"
(   52) 
(   53) # 53
(   53) CONTAINS
(   54) 
(   55) 
(   56) !!DB: Once-per-day, compute globally-averaged (ie domain-averaged):
(   57) !!    T, S, emp, ssh, rnf, OBC TS flux, ice_area, ice_vol,  ...
(   58) !! Every dt: Output a check velocity from the Bay-of-Fundy
(   59) !! NB: Volume calc is done every call in case free sfce ever becomes part
(   60) !!     of e3t(1). 
(   61) !!     Also cell_volume(i,j,k) could be computed at nit000 and saved but 
(   62) !!     this is not done as calc occurs only once-per-day
(   63) 
(   64) 
(   65)   SUBROUTINE sopa_mc_diagnostics( kt )
(   66)     !!---------------------------------------------------------------------
(   67)     !!
(   68)     !! ** Purpose : Compute some domain-averaged quantities 
(   69)     !!              Output to text file(s)
(   70)     !!
(   71)     !!	
(   72) 
(   73)     !!----------------------------------------------------------------------
(   74)     !! * Arguments
(   75)     INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   76) 
(   77)     !! * Local declarations
(   78)     INTEGER  ::   ji, jj, jk
(   79)     Logical ::  BoF_output
(   80)     CHARACTER*175 :: header
(   81)     REAL(wp) ::   zwei, zsrau, & 
(   82)          uT_east, uT_west, vT_south, vT_north, &
(   83)          uS_east, uS_west, vS_south, vS_north, &
(   84)          dT_sfce, dS_sfce, net_T_flux, net_S_flux, &
(   85)          ice_area, ice_vol, ave_T_sfce, ave_S_sfce, vol_sfce
(   86) 
(   87)     !!----------------------------------------------------------------------
(   88) 






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 11

(   89) 
(   90) !!DB: The dT & dS_sfce calcs mimic the ones found in trasbc.F90, viz.
(   91) !            zta = ro0cpr * ( qt(ji,jj) - qsr(ji,jj) ) * zse3t
(   92) !            zsa = emps(ji,jj) * zsrau * sn(ji,jj,1) * zse3t
(   93) 
(   94) 
(   95)     BoF_output = .true.    ! false ===> do not output the test vel in the BoF
(   96) 
(   97)     if( kt == nit000 ) then
(   98)        IF(lwp) THEN
(   99)           WRITE(numout2,*)
(  100)           WRITE(numout2,*) 'DB: SOPA_MC_DIAGNOSTICS'
(  101)           WRITE(numout2,*) '~~~~~~~'
(  102)        ENDIF
(  103) 
(  104)        ice_area = 0.e0
(  105)        ice_vol = 0.e0
(  106)        dT_sfce = 0.e0    !K/s
(  107)        dS_sfce = 0.e0    !S/s
(  108)        glb_emp    = 0.e0
(  109)        glb_precip = 0.e0
(  110)        glb_rnf    = 0.e0
(  111)        glb_ssh   = 0.e0   ! averaged sea surface height 
(  112)        glb_S = 0.e0   ! averaged ocean salinity 
(  113)        glb_T = 0.e0   ! averaged ocean temperature
(  114) !!DB OBC variables; assigned even if 1 = off 
(  115)        uT_east = 0.e0
(  116)        uT_west = 0.e0
(  117)        vT_south = 0.e0
(  118)        vT_north = 0.e0
(  119)        uS_east = 0.e0
(  120)        uS_west = 0.e0
(  121)        vS_south = 0.e0
(  122)        vS_north = 0.e0
(  123) 
(  124) 
(  125) 
(  126) !! da(i,j) & surface area
(  127)        e1e2_i(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:)
(  128)        sfce_area    = SUM(e1e2_i(:,:))
(  129)        IF( lk_mpp )   CALL  mpp_sum( sfce_area    )   ! sum over the global domain
(  130) 
(  131)        if(BoF_output) then
(  132) !!       isolate the correct processor to open the file
(  133)           if(mi0(58)==mi1(58) .AND. mj0(106)== mj1(106)) then
(  134)              open(999,file='BoF_vel.dat')
(  135)           endif
(  136)        endif
(  137) 
(  138)        if(lwp) then
(  139)           open(998,file='sopa_mc_diagnostics.dat')
(  140)           header = 'year day kt  [ave_T ave_S ave_SSH] [OBC(VT(E,W,S) VS(E,W,S))] &
(  141)                dT_sfce dS_sfce net_T_flux net_S_flux (TS/sec) [ave(emp,precip,rnf)] ice_area ice_vol &
(  142)                [ave_T ave_S k=1,2]'
(  143)           write(998,'(A175)')header
(  144)        endif
(  145) 
(  146)     endif  !!kt = nit000 loop






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 12

(  147) 
(  148)     if(BoF_output) then
(  149) !!DBG: Output a velocity from the Bay of Fundy 
(  150) !!isolate processor that owns this global location
(  151)        ji = 58; jj = 106
(  152)        if(mi0(ji)==mi1(ji) .AND. mj0(jj)== mj1(jj)) then
(  153) !!OLD+NEW
(  154)           write(999,'(i7,2x,12(f8.3,1x))')kt, &
(  155)                (un(ji,jj,1)+un(ji,jj,2))/(e3t(1)+e3t(2)), &
(  156)                (vn(ji,jj,1)+vn(ji,jj,2))/(e3t(1)+e3t(2))
(  157) 
(  158) !          un(ji,jj,1), vn(ji,jj,1),un(ji,jj,3),vn(ji,jj,3)
(  159) !!NEW: decimal year, kt, Vel(i,j,1:2)/dz(1:2) ... in progress ...
(  160) !             write(999,'(i7,2x,12(f8.3,1x))')nyear, nday_year, 
(  161) !                  (un(ji,jj,1)+un(ji,jj,2))/(e3t(1)+e3t(2)), &
(  162) !                  (vn(ji,jj,1)+vn(ji,jj,2))/(e3t(1)+e3t(2))
(  163)           endif
(  164)        !!DB:END
(  165)     endif
(  166) 
(  167) 
(  168) !!DB: Do once-per-day -- Hardwired
(  169)       if(mod(kt-nit000,int(rday/rdt)) ==0) then
(  170) 
(  171)          zsrau = 1. / rauw
(  172)          dT_sfce = sum(e1e2_i(:,:)*qt(:,:)*ro0cpr)  ! >0 ===> dT > 0
(  173)          if(lk_mpp) call mpp_sum(dT_sfce)
(  174) !         dS_sfce = sum(e1e2_i(:,:)*emps(:,:)*sn(:,:,1)*zsrau)  ! >0 ===> dS > 0
(  175) !!Later must divide by glb_S        not           |||
(  176)          dS_sfce = sum(e1e2_i(:,:)*emps(:,:)*zsrau)  ! >0 ===> dS > 0
(  177)          if(lk_mpp) call mpp_sum(dS_sfce)
(  178) 
(  179)          glb_emp    = SUM( e1e2_i(:,:) * emp   (:,:) )
(  180)          IF( lk_mpp )   CALL  mpp_sum( glb_emp    )   ! sum over the global domain
(  181) # 182
(  182)          glb_precip = SUM( e1e2_i(:,:) * watm  (:,:) )
(  183)          IF( lk_mpp )   CALL  mpp_sum( glb_precip )   ! sum over the global domain
(  184) # 185
(  185)          glb_rnf    = SUM( e1e2_i(:,:) * runoff(:,:) )
(  186)          IF( lk_mpp )   CALL  mpp_sum( glb_rnf    )   ! sum over the global domain
(  187)          
(  188)          glb_ssh = SUM( e1e2_i(:,:) * sshn(:,:) )
(  189)          if( lk_mpp ) call mpp_sum(glb_ssh)
(  190) 
(  191) # 192
(  192)          ice_area = SUM(e1e2_i(:,:)*(1.0-frld(:,:)))  !m2
(  193)          ice_vol = SUM(e1e2_i(:,:)*(1.0-frld(:,:))*hicif(:,:)) !m3
(  194)          if( lk_mpp ) then 
(  195)             call mpp_sum(ice_area)
(  196)             call mpp_sum(ice_vol)
(  197)          endif
(  198) # 199
(  199)          vol_tot  = 0.e0
(  200)          glb_S = 0.0
(  201)          glb_T = 0.0
(  202)          do jk = 1, jpkm1   
(  203)             do jj = 2, jpjm1
(  204)                do ji = 2, jpim1   ! vector opt.






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 13

(  205)                   zwei  = e1e2_i(ji,jj) * e3t(jk) * tmask(ji,jj,jk)
(  206)                   glb_S = glb_S + ( sn(ji,jj,jk) ) * zwei
(  207)                   glb_T = glb_T + ( tn(ji,jj,jk) ) * zwei
(  208)                   vol_tot  = vol_tot  + zwei
(  209)                enddo
(  210)             enddo
(  211)          enddo
(  212)          vol_sfce = 0.e0
(  213)          ave_T_sfce = 0.e0
(  214)          ave_S_sfce = 0.e0
(  215)          do jk = 1, 2
(  216)             do jj = 2, jpjm1
(  217)                do ji = 2, jpim1   ! vector opt.
(  218)                   zwei  = e1e2_i(ji,jj) * e3t(jk) * tmask(ji,jj,jk)
(  219)                   ave_S_sfce = ave_S_sfce + ( sn(ji,jj,jk) ) * zwei
(  220)                   ave_T_sfce = ave_T_sfce + ( tn(ji,jj,jk) ) * zwei
(  221)                   vol_sfce  = vol_sfce + zwei
(  222)                enddo
(  223)             enddo
(  224)          enddo
(  225)          IF( lk_mpp ) THEN
(  226)             call  mpp_sum( glb_S )        
(  227)             call  mpp_sum( glb_T )        
(  228)             call  mpp_sum( vol_tot )        
(  229)             call  mpp_sum( ave_S_sfce )        
(  230)             call  mpp_sum( ave_T_sfce )        
(  231)             call  mpp_sum( vol_sfce )        
(  232)          ENDIF
(  233) 
(  234) # 235
(  235) !!DB -- OB T & S flux 
(  236) !!NB: because this routine is called after fields are updated, the
(  237) !!    now variables _should_ contain the correct values
(  238) 
(  239)          if( lp_obc_east  ) then
(  240)             uT_east = 0.e0;  uS_east = 0.e0
(  241)             do ji = nie0, nie1
(  242)                do jj = 1, jpj
(  243)                   do jk = 1, jpk
(  244)                      uT_east = uT_east +  & 
(  245)                           un(ji,jj,jk)*tn(ji+1,jj,jk)*e2u(ji,jj)*e3t(jk)*tmask(ji,jj,jk)
(  246)                      uS_east = uS_east +  & 
(  247)                           un(ji,jj,jk)*sn(ji+1,jj,jk)*e2u(ji,jj)*e3t(jk)*tmask(ji,jj,jk)
(  248)                   enddo
(  249)                enddo
(  250)             enddo
(  251)             if( lk_mpp ) then
(  252)                call mpp_sum(uT_east)
(  253)                call mpp_sum(uS_east)
(  254)             endif
(  255) 
(  256)          endif
(  257)          if( lp_obc_west  ) then
(  258)             uT_west = 0.e0;  uS_west = 0.e0
(  259)             do ji = niw0, niw1
(  260)                do jj = 1, jpj
(  261)                   do jk = 1, jpk
(  262)                      uT_west = uT_west + &






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 14

(  263)                           un(ji,jj,jk)*tn(ji,jj,jk)*e2u(ji,jj)*e3t(jk)*tmask(ji,jj,jk)
(  264)                      uS_west = uS_west + &
(  265)                           un(ji,jj,jk)*sn(ji,jj,jk)*e2u(ji,jj)*e3t(jk)*tmask(ji,jj,jk)
(  266)                   enddo
(  267)                enddo
(  268)             enddo
(  269)             if( lk_mpp ) then
(  270)                call mpp_sum(uT_west)
(  271)                call mpp_sum(uS_west)
(  272)             endif
(  273)          endif
(  274) 
(  275)          if( lp_obc_south ) then
(  276)             vT_south = 0.e0;  vS_south = 0.e0
(  277)             do jj = njs0, njs1
(  278)                do ji = 1, jpi
(  279)                   do jk = 1, jpk
(  280)                      vT_south = vT_south + &
(  281)                           vn(ji,jj,jk)*tn(ji,jj,jk)*e1v(ji,jj)*e3t(jk)*tmask(ji,jj,jk)
(  282)                      vS_south = vS_south + &
(  283)                           vn(ji,jj,jk)*sn(ji,jj,jk)*e1v(ji,jj)*e3t(jk)*tmask(ji,jj,jk)
(  284)                   enddo
(  285)                enddo
(  286)             enddo
(  287)             if( lk_mpp ) then
(  288)                call mpp_sum(vT_south)
(  289)                call mpp_sum(vS_south)
(  290)             endif
(  291)          endif
(  292) 
(  293) !!DB -- leave blank as SOPA-MC north is closed (lazy)
(  294) !         if( lp_obc_north ) then
(  295)          !write code later
(  296) !           vT_north = 0.e0
(  297) !         endif
(  298) 
(  299) 
(  300) 
(  301) # 302
(  302) !! OUTPUT: average quantities over vol or area
(  303) !! NB: OBC fluxes are divided-by volume ---> units deg/s 
(  304) !! NB: leave emp, precip, runoff in whatever units they exist in 
(  305)          if(lwp) then
(  306)             glb_T = glb_T/vol_tot
(  307)             glb_S = glb_S/vol_tot
(  308)             glb_ssh = glb_ssh/sfce_area
(  309)             net_T_flux = -(uT_east/vol_tot-uT_west/vol_tot)  &
(  310)                  -(0.e0 - vT_south/vol_tot) + dT_sfce/vol_tot
(  311)             net_S_flux = -(uS_east/vol_tot-uS_west/vol_tot)  &
(  312)                  -(0.e0 - vS_south/vol_tot) + dS_sfce/(vol_tot*glb_S)
(  313)             write(998,'(2(i4,1x),i9,1x,3(f10.4,1x),19(e10.4,1x))')     &
(  314)                  nyear, nday_year, kt, & 
(  315)                  glb_T, glb_S, glb_ssh, &
(  316)                  uT_east/vol_tot,uT_west/vol_tot, vT_south/vol_tot, &
(  317)                  uS_east/vol_tot,uS_west/vol_tot, vS_south/vol_tot, &
(  318)                  dT_sfce/vol_tot, dS_sfce/(vol_tot*glb_S), &
(  319)                  net_T_flux, net_S_flux, &
(  320)                  glb_emp/sfce_area, glb_precip/sfce_area, glb_rnf/sfce_area, &






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 15

(  321)                  ice_area/1.e6, ice_vol/1.e9, &
(  322)                  ave_T_sfce/vol_sfce, ave_S_sfce/vol_sfce
(  323)          endif
(  324) 
(  325) !AD:
(  326)          call vertical_mixing_diagnostics(kt)
(  327) 
(  328)       endif     !! END (kt-nit000) % rday/rdt
(  329) 
(  330) 
(  331) 
(  332) !!Old code, but keep around for awhile
(  333)          ! Conversion in m3
(  334) !         glb_emp    = glb_emp    * rdttra(1) * 1.e-3 
(  335) !         glb_precip = glb_precip * rdttra(1) * 1.e-3 / rday
(  336) !         glb_rnf    = glb_rnf    * rdttra(1) * 1.e-3
(  337)     ! Ecriture des diagnostiques 
(  338)     ! --------------------------
(  339) 
(  340) !    IF( kt == nitend ) THEN
(  341) 
(  342) !       WRITE(inum,*)    'Net freshwater budget '
(  343) !       WRITE(inum,9010) '  emp    = ', a_emp   , ' m3 =', a_emp   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  344) !       WRITE(inum,9010) '  precip = ', a_precip, ' m3 =', a_precip/((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  345) !       WRITE(inum,9010) '  a_rnf  = ', a_rnf   , ' m3 =', a_rnf   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
(  346) !       WRITE(inum,*)    'Mean sea level : '
(  347) !       WRITE(inum,9010) '  at nit000 = ',a_sshb           ,' m3 '
(  348) !       WRITE(inum,9010) '  at nitend = ',a_sshend         ,' m3 '
(  349) !       WRITE(inum,9010) '  diff      = ',(a_sshend-a_sshb),' m3 =', (a_sshend-a_sshb)/((nitend-nit000+1)*rdt) * 1.e-6,' Sv'
(  350) !       WRITE(inum,9020) '  mean sea level elevation    =', a_sshend/zarea,' m'
(  351) !    ENDIF
(  352) 
(  353) 
(  354)   END SUBROUTINE sopa_mc_diagnostics
(  355) 
(  356) 
(  357)   SUBROUTINE vertical_mixing_diagnostics(kt)
(  358) 
(  359)     !!---------------------------------------------------------------------
(  360)     !!
(  361)     !! ** Purpose : Diagnostic output of various vertical mixing parameters
(  362)     !!              (Currently just the vertical viscosity coeff)
(  363)     !!  
(  364)     !! ** Method
(  365)     !! 
(  366)     !!   domain which contains the global OUTi, OUTj will compute the 
(  367)     !!   vertical maxval. Other domains remain at     maxV = 0.
(  368)     !!   call mpp_max(maxV) ensures all domains are updated. Only lwp 
(  369)     !!   writes output. For the global maximum, the generic mpp_maxloc is used.
(  370)     !!   It always returns a value from boundary. Might need to modify it to
(  371)     !!   look at internal domain only. 
(  372)     !!     
(  373)     !!	
(  374)     !!----------------------------------------------------------------------
(  375) 
(  376)     USE zdf_oce
(  377) 
(  378)     INTEGER, INTENT( in ) ::   kt      ! ocean time-step index






PGF90 (Version     10.5)          09/23/2010  09:50:09      page 16

(  379)     
(  380)     !! * Local declarations
(  381)     INTEGER  ::  jj,ji,fid
(  382)     REAL(wp) :: maxV,maxV_Tot
(  383)     INTEGER  :: MAXi,MAXj,MAXk
(  384)     INTEGER  :: OUTi_local, OUTj_local,maxV_local
(  385)     INTEGER  :: OUTi, OUTj
(  386)     parameter (OUTi=67,OUTj=63) !global location for desired output
(  387) 
(  388)     !!----------------------------------------------------------------------
(  389) 
(  390)     fid=numout2
(  391)     
(  392)     
(  393)     maxV = 0
(  394)     OUTj_local=OUTj + 1 -  njmpp
(  395)     OUTi_local=OUTi + 1 -  nimpp    
(  396)     IF  ( OUTi_local .le. jpi .AND. OUTi_local .ge. 1 .AND. &
(  397)             OUTj_local .le. jpj .AND. OUTj_local .ge. 1) then
(  398) 	  maxV = MAXVAL( avmu(OUTi_local,OUTj_local,:), & 
(  399) 	  mask = umask(OUTi_local,OUTj_local,:) == 1.e0)
(  400) !!DBG
(  401) !          write(3001,*) kt,avmu(OUTi_local,OUTj_local,1:mbathy(OUTi_local,OUTj_local))
(  402)     ENDIF
(  403) 
(  404)     IF(lk_mpp) call mpp_max(maxV)
(  405)     CALL mpp_maxloc( avmu,tmask,maxV_Tot,MAXi,MAXj,MAXk )
(  406) 
(  407)     IF(lwp) THEN 
(  408)        write(fid,'( "VERTICAL_MIXING_DIAGNOSTICS:"  )' )
(  409)        write(fid,*) '     @ kt,OUTi,OUTj,OUTi_local,OUTj_local', &
(  410)                        kt,OUTi,OUTj,OUTi_local,OUTj_local
(  411)        write(fid,*) '   glamu,gphiu = ', glamu(OUTi_local,OUTj_local) &
(  412)                                        , gphiu(OUTi_local,OUTj_local)  
(  413) ! AD. NOTE: glamu,gphiu not working!!!
(  414) !		       
(  415) 		       
(  416)        write(fid,'( "     vertical viscosity coeff. at uw-point"  )')
(  417)        write(fid,'( "MAXVAL avmu at: ",I4,I4," = ",E10.4,", &
(  418)                     Overall MAXVAL at: " ,I4,I4,I4," = ",E10.4)')  &  
(  419)               	    OUTi, OUTj,maxV,MAXi,MAXj,MAXk,maxV_Tot           
(  420)     ENDIF
(  421) 
(  422)   END SUBROUTINE vertical_mixing_diagnostics
(  423) 
(  424) 
(  425) END MODULE sopa_mc
