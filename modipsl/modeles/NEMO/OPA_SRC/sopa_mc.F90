!!DB 2009.08
!!Module that contains various (diagnostic) routines specific to the
!!Maritime Canada Shelf OPA model
!!Quantities computed once per day

!!Started as flxfwb; 
!! kept many of the variables although they may not be used

MODULE sopa_mc

  !! * Modules used
  USE oce             ! ocean dynamics and tracers
  USE dom_oce         ! ocean space and time domain
  USE cpl_oce         ! coupled atmosphere/ocean
  USE phycst          ! physical constants
  USE in_out_manager  ! I/O manager
  USE lib_mpp         ! distribued memory computing library
  USE flxrnf          ! ocean runoffs
  USE ocesbc          ! ocean surface boudaries conditions
  USE blk_oce
  USE flxblk          ! bulk formulea
  USE daymod          ! calendar
#if defined key_obc
  USE obc_oce
#endif
#if defined key_ice_lim
  USE ice
#endif

  IMPLICIT NONE
  PRIVATE

  !! * Routine accessibility
  PUBLIC sopa_mc_diagnostics      ! routine called by step

  !! * Share module variables
  REAL(wp), PUBLIC ::   &  !:
       vol_tot, sfce_area,  &
       glb_emp   ,  & ! domain averaged evaporation minus precipitation
       glb_precip,  & ! domain averaged precipitation
       glb_rnf   ,  & ! domain averaged runoff
       glb_ssh  ,  & ! domain averaged sea surface height
       glb_S   ,  & ! domain averaged ocean salinity
       glb_T        ! domain averaged ocean temperature

  REAL(wp), DIMENSION(jpi,jpj) ::  &
       e1e2_i    ! area of the interior domain (e1t*e2t*tmask_i)

  !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"

CONTAINS


!!DB: Once-per-day, compute globally-averaged (ie domain-averaged):
!!    T, S, emp, ssh, rnf, OBC TS flux, ice_area, ice_vol,  ...
!! Every dt: Output a check velocity from the Bay-of-Fundy
!! NB: Volume calc is done every call in case free sfce ever becomes part
!!     of fse3t(i,j,1). 
!!     Also cell_volume(i,j,k) could be computed at nit000 and saved but 
!!     this is not done as calc occurs only once-per-day


  SUBROUTINE sopa_mc_diagnostics( kt )
    !!---------------------------------------------------------------------
    !!
    !! ** Purpose : Compute some domain-averaged quantities 
    !!              Output to text file(s)
    !!
    !!	

    !!----------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

    !! * Local declarations
    INTEGER  ::   ji, jj, jk
    Logical ::  BoF_output
    CHARACTER*175 :: header
    REAL(wp) ::   zwei, zsrau, & 
         uT_east, uT_west, vT_south, vT_north, &
         uS_east, uS_west, vS_south, vS_north, &
         dT_sfce, dS_sfce, net_T_flux, net_S_flux, &
         ice_area, ice_vol, ave_T_sfce, ave_S_sfce, vol_sfce

    !!----------------------------------------------------------------------


!!DB: The dT & dS_sfce calcs mimic the ones found in trasbc.F90, viz.
!            zta = ro0cpr * ( qt(ji,jj) - qsr(ji,jj) ) * zse3t
!            zsa = emps(ji,jj) * zsrau * sn(ji,jj,1) * zse3t


    BoF_output = .false.    ! false ===> do not output the test vel in the BoF

    if( kt == nit000 ) then
       IF(lwp) THEN
          WRITE(numout2,*)
          WRITE(numout2,*) 'DB: SOPA_MC_DIAGNOSTICS'
          WRITE(numout2,*) '~~~~~~~'
       ENDIF

       ice_area = 0.e0
       ice_vol = 0.e0
       dT_sfce = 0.e0    !K/s
       dS_sfce = 0.e0    !S/s
       glb_emp    = 0.e0
       glb_precip = 0.e0
       glb_rnf    = 0.e0
       glb_ssh   = 0.e0   ! averaged sea surface height 
       glb_S = 0.e0   ! averaged ocean salinity 
       glb_T = 0.e0   ! averaged ocean temperature
!!DB OBC variables; assigned even if key_obc = off 
       uT_east = 0.e0
       uT_west = 0.e0
       vT_south = 0.e0
       vT_north = 0.e0
       uS_east = 0.e0
       uS_west = 0.e0
       vS_south = 0.e0
       vS_north = 0.e0



!! da(i,j) & surface area
       e1e2_i(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:)
       sfce_area    = SUM(e1e2_i(:,:))
       IF( lk_mpp )   CALL  mpp_sum( sfce_area    )   ! sum over the global domain

       if(BoF_output) then
!!       isolate the correct processor to open the file
          if(mi0(58)==mi1(58) .AND. mj0(106)== mj1(106)) then
             open(999,file='BoF_vel.dat')
          endif
       endif

       if(lwp) then
          open(998,file='sopa_mc_diagnostics.dat')
          header = 'year day kt  [ave_T ave_S ave_SSH] [OBC(VT(E,W,S) VS(E,W,S))] &
               dT_sfce dS_sfce net_T_flux net_S_flux (TS/sec) [ave(emp,precip,rnf)] ice_area ice_vol &
               [ave_T ave_S k=1,2]'
          write(998,'(A175)')header
       endif

    endif  !!kt = nit000 loop

    if(BoF_output) then
!!DBG: Output a velocity from the Bay of Fundy 
!!isolate processor that owns this global location
       ji = 58; jj = 106
       if(mi0(ji)==mi1(ji) .AND. mj0(jj)== mj1(jj)) then
!!OLD+NEW
          write(999,'(i7,2x,12(f8.3,1x))')kt, &
               (un(ji,jj,1)+un(ji,jj,2))/(fse3t(ji,jj,1)+fse3t(ji,jj,2)), &
               (vn(ji,jj,1)+vn(ji,jj,2))/(fse3t(ji,jj,1)+fse3t(ji,jj,2))

!          un(ji,jj,1), vn(ji,jj,1),un(ji,jj,3),vn(ji,jj,3)
!!NEW: decimal year, kt, Vel(i,j,1:2)/dz(1:2) ... in progress ...
!             write(999,'(i7,2x,12(f8.3,1x))')nyear, nday_year, 
!                  (un(ji,jj,1)+un(ji,jj,2))/(fse3t(ji,jj,1)+fse3t(ji,jj,2)), &
!                  (vn(ji,jj,1)+vn(ji,jj,2))/(fse3t(ji,jj,1)+fse3t(ji,jj,2))
          endif
       !!DB:END
    endif


!!DB: Do once-per-day -- Hardwired
      if(mod(kt-nit000,int(rday/rdt)) ==0) then

         zsrau = 1. / rauw
         dT_sfce = sum(e1e2_i(:,:)*qt(:,:)*ro0cpr)  ! >0 ===> dT > 0
         if(lk_mpp) call mpp_sum(dT_sfce)
!         dS_sfce = sum(e1e2_i(:,:)*emps(:,:)*sn(:,:,1)*zsrau)  ! >0 ===> dS > 0
!!Later must divide by glb_S        not           |||
         dS_sfce = sum(e1e2_i(:,:)*emps(:,:)*zsrau)  ! >0 ===> dS > 0
         if(lk_mpp) call mpp_sum(dS_sfce)

         glb_emp    = SUM( e1e2_i(:,:) * emp   (:,:) )
         IF( lk_mpp )   CALL  mpp_sum( glb_emp    )   ! sum over the global domain
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
         glb_precip = SUM( e1e2_i(:,:) * watm  (:,:) )
         IF( lk_mpp )   CALL  mpp_sum( glb_precip )   ! sum over the global domain
#endif
         glb_rnf    = SUM( e1e2_i(:,:) * runoff(:,:) )
         IF( lk_mpp )   CALL  mpp_sum( glb_rnf    )   ! sum over the global domain
         
         glb_ssh = SUM( e1e2_i(:,:) * sshn(:,:) )
         if( lk_mpp ) call mpp_sum(glb_ssh)

#if defined key_ice_lim
         ice_area = SUM(e1e2_i(:,:)*(1.0-frld(:,:)))  !m2
         ice_vol = SUM(e1e2_i(:,:)*(1.0-frld(:,:))*hicif(:,:)) !m3
         if( lk_mpp ) then 
            call mpp_sum(ice_area)
            call mpp_sum(ice_vol)
         endif
#endif
         vol_tot  = 0.e0
         glb_S = 0.0
         glb_T = 0.0
         do jk = 1, jpkm1   
            do jj = 2, jpjm1
               do ji = fs_2, fs_jpim1   ! vector opt.
                  zwei  = e1e2_i(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)
                  glb_S = glb_S + ( sn(ji,jj,jk) ) * zwei
                  glb_T = glb_T + ( tn(ji,jj,jk) ) * zwei
                  vol_tot  = vol_tot  + zwei
               enddo
            enddo
         enddo
         vol_sfce = 0.e0
         ave_T_sfce = 0.e0
         ave_S_sfce = 0.e0
         do jk = 1, 2
            do jj = 2, jpjm1
               do ji = fs_2, fs_jpim1   ! vector opt.
                  zwei  = e1e2_i(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)
                  ave_S_sfce = ave_S_sfce + ( sn(ji,jj,jk) ) * zwei
                  ave_T_sfce = ave_T_sfce + ( tn(ji,jj,jk) ) * zwei
                  vol_sfce  = vol_sfce + zwei
               enddo
            enddo
         enddo
         IF( lk_mpp ) THEN
            call  mpp_sum( glb_S )        
            call  mpp_sum( glb_T )        
            call  mpp_sum( vol_tot )        
            call  mpp_sum( ave_S_sfce )        
            call  mpp_sum( ave_T_sfce )        
            call  mpp_sum( vol_sfce )        
         ENDIF

#ifdef key_obc
!!DB -- OB T & S flux 
!!NB: because this routine is called after fields are updated, the
!!    now variables _should_ contain the correct values

         if( lp_obc_east  ) then
            uT_east = 0.e0;  uS_east = 0.e0
            do ji = nie0, nie1
               do jj = 1, jpj
                  do jk = 1, jpk
                     uT_east = uT_east +  & 
                          un(ji,jj,jk)*tn(ji+1,jj,jk)*e2u(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
                     uS_east = uS_east +  & 
                          un(ji,jj,jk)*sn(ji+1,jj,jk)*e2u(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
                  enddo
               enddo
            enddo
            if( lk_mpp ) then
               call mpp_sum(uT_east)
               call mpp_sum(uS_east)
            endif

         endif
         if( lp_obc_west  ) then
            uT_west = 0.e0;  uS_west = 0.e0
            do ji = niw0, niw1
               do jj = 1, jpj
                  do jk = 1, jpk
                     uT_west = uT_west + &
                          un(ji,jj,jk)*tn(ji,jj,jk)*e2u(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
                     uS_west = uS_west + &
                          un(ji,jj,jk)*sn(ji,jj,jk)*e2u(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
                  enddo
               enddo
            enddo
            if( lk_mpp ) then
               call mpp_sum(uT_west)
               call mpp_sum(uS_west)
            endif
         endif

         if( lp_obc_south ) then
            vT_south = 0.e0;  vS_south = 0.e0
            do jj = njs0, njs1
               do ji = 1, jpi
                  do jk = 1, jpk
                     vT_south = vT_south + &
                          vn(ji,jj,jk)*tn(ji,jj,jk)*e1v(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
                     vS_south = vS_south + &
                          vn(ji,jj,jk)*sn(ji,jj,jk)*e1v(ji,jj)*fse3t(ji,jj,jk)*tmask(ji,jj,jk)
                  enddo
               enddo
            enddo
            if( lk_mpp ) then
               call mpp_sum(vT_south)
               call mpp_sum(vS_south)
            endif
         endif

!!DB -- leave blank as SOPA-MC north is closed (lazy)
!         if( lp_obc_north ) then
         !write code later
!           vT_north = 0.e0
!         endif

#endif   !!key_obc


!! OUTPUT: average quantities over vol or area
!! NB: OBC fluxes are divided-by volume ---> units deg/s 
!! NB: leave emp, precip, runoff in whatever units they exist in 
         if(lwp) then
            glb_T = glb_T/vol_tot
            glb_S = glb_S/vol_tot
            glb_ssh = glb_ssh/sfce_area
            net_T_flux = -(uT_east/vol_tot-uT_west/vol_tot)  &
                 -(0.e0 - vT_south/vol_tot) + dT_sfce/vol_tot
            net_S_flux = -(uS_east/vol_tot-uS_west/vol_tot)  &
                 -(0.e0 - vS_south/vol_tot) + dS_sfce/(vol_tot*glb_S)
            write(998,'(2(i4,1x),i9,1x,3(f10.4,1x),19(e10.4,1x))')     &
                 nyear, nday_year, kt, & 
                 glb_T, glb_S, glb_ssh, &
                 uT_east/vol_tot,uT_west/vol_tot, vT_south/vol_tot, &
                 uS_east/vol_tot,uS_west/vol_tot, vS_south/vol_tot, &
                 dT_sfce/vol_tot, dS_sfce/(vol_tot*glb_S), &
                 net_T_flux, net_S_flux, &
                 glb_emp/sfce_area, glb_precip/sfce_area, glb_rnf/sfce_area, &
                 ice_area/1.e6, ice_vol/1.e9, &
                 ave_T_sfce/vol_sfce, ave_S_sfce/vol_sfce
         endif

!AD:
         call vertical_mixing_diagnostics(kt)

      endif     !! END (kt-nit000) % rday/rdt



!!Old code, but keep around for awhile
         ! Conversion in m3
!         glb_emp    = glb_emp    * rdttra(1) * 1.e-3 
!         glb_precip = glb_precip * rdttra(1) * 1.e-3 / rday
!         glb_rnf    = glb_rnf    * rdttra(1) * 1.e-3
    ! Ecriture des diagnostiques 
    ! --------------------------

!    IF( kt == nitend ) THEN

!       WRITE(inum,*)    'Net freshwater budget '
!       WRITE(inum,9010) '  emp    = ', a_emp   , ' m3 =', a_emp   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
!       WRITE(inum,9010) '  precip = ', a_precip, ' m3 =', a_precip/((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
!       WRITE(inum,9010) '  a_rnf  = ', a_rnf   , ' m3 =', a_rnf   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
!       WRITE(inum,*)    'Mean sea level : '
!       WRITE(inum,9010) '  at nit000 = ',a_sshb           ,' m3 '
!       WRITE(inum,9010) '  at nitend = ',a_sshend         ,' m3 '
!       WRITE(inum,9010) '  diff      = ',(a_sshend-a_sshb),' m3 =', (a_sshend-a_sshb)/((nitend-nit000+1)*rdt) * 1.e-6,' Sv'
!       WRITE(inum,9020) '  mean sea level elevation    =', a_sshend/zarea,' m'
!    ENDIF


  END SUBROUTINE sopa_mc_diagnostics


  SUBROUTINE vertical_mixing_diagnostics(kt)

    !!---------------------------------------------------------------------
    !!
    !! ** Purpose : Diagnostic output of various vertical mixing parameters
    !!              (Currently just the vertical viscosity coeff)
    !!  
    !! ** Method
    !! 
    !!   domain which contains the global OUTi, OUTj will compute the 
    !!   vertical maxval. Other domains remain at     maxV = 0.
    !!   call mpp_max(maxV) ensures all domains are updated. Only lwp 
    !!   writes output. For the global maximum, the generic mpp_maxloc is used.
    !!   It always returns a value from boundary. Might need to modify it to
    !!   look at internal domain only. 
    !!     
    !!	
    !!----------------------------------------------------------------------

    USE zdf_oce

    INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
    
    !! * Local declarations
    INTEGER  ::  jj,ji,fid
    REAL(wp) :: maxV,maxV_Tot
    INTEGER  :: MAXi,MAXj,MAXk
    INTEGER  :: OUTi_local, OUTj_local,maxV_local
    INTEGER  :: OUTi, OUTj
    parameter (OUTi=67,OUTj=63) !global location for desired output

    !!----------------------------------------------------------------------

    fid=numout2
    
    
    maxV = 0
    OUTj_local=OUTj + 1 -  njmpp
    OUTi_local=OUTi + 1 -  nimpp    
    IF  ( OUTi_local .le. jpi .AND. OUTi_local .ge. 1 .AND. &
            OUTj_local .le. jpj .AND. OUTj_local .ge. 1) then
	  maxV = MAXVAL( avmu(OUTi_local,OUTj_local,:), & 
	  mask = umask(OUTi_local,OUTj_local,:) == 1.e0)
!!DBG
!          write(3001,*) kt,avmu(OUTi_local,OUTj_local,1:mbathy(OUTi_local,OUTj_local))
    ENDIF

    IF(lk_mpp) call mpp_max(maxV)
    CALL mpp_maxloc( avmu,tmask,maxV_Tot,MAXi,MAXj,MAXk )

    IF(lwp) THEN 
       write(fid,'( "VERTICAL_MIXING_DIAGNOSTICS:"  )' )
       write(fid,*) '     @ kt,OUTi,OUTj,OUTi_local,OUTj_local', &
                       kt,OUTi,OUTj,OUTi_local,OUTj_local
       write(fid,*) '   glamu,gphiu = ', glamu(OUTi_local,OUTj_local) &
                                       , gphiu(OUTi_local,OUTj_local)  
! AD. NOTE: glamu,gphiu not working!!!
!		       
		       
       write(fid,'( "     vertical viscosity coeff. at uw-point"  )')
       write(fid,'( "MAXVAL avmu at: ",I4,I4," = ",E10.4,", &
                    Overall MAXVAL at: " ,I4,I4,I4," = ",E10.4)')  &  
              	    OUTi, OUTj,maxV,MAXi,MAXj,MAXk,maxV_Tot           
    ENDIF

  END SUBROUTINE vertical_mixing_diagnostics


END MODULE sopa_mc
