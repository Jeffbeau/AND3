!!DB
!!2008.05 -- Added tmask_i as in (latest) obcvol to account for inter CPU
!!overlap problems. This fixed a volume conservation problem which affected
!!my results. 
!! **(NB: should check shape factors e2t etc in obc transport calcs) 
!!I have added a lot of DBG code. Instead of deleting it in this version, I
!!just comment it out, as it is useful to retain for future debugging
Module obcctl
#ifdef	key_obc
  !!=================================================================================
  !!                       ***  MODULE  obcctl  ***
  !! Ocean dynamics:    Balancing  the water mass in the domain to avoid potential
  !!                    unbalance due to OBC 
  !!=================================================================================


  !!----------------------------------------------------------------------------------
  !! * Modules used
  USE oce             ! ocean dynamics and tracers
  USE dom_oce         ! ocean space and time domain
  USE dynspg_oce      ! surface pressure gradient variables
  USE obc_oce         ! ocean open boundary conditions
  USE phycst          ! physical constants
  USE lib_mpp
  USE in_out_manager  ! I/O logical units

  IMPLICIT NONE
  PRIVATE
# include "domzgr_substitute.h90"
  !! * Accessibility
  PUBLIC obc_ctl  ! 

  !!---------------------------------------------------------------------------------
  !!  OPA XX , BIO 
  !!---------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE obc_ctl(kt)
     !!------------------------------------------------------------------------------
     !!                      SUBROUTINE obc_ctl
     !!                     **********************
     !! ** Purpose :
     !!
     !! History : 
     !! xxx !  07-07  (Z. Wang) original
     !! xxx !  07-27  DB modified
     !! DB: 
     !! (1) Adjusted sum_obc to mirror divergence calc, i.e. div = (east-west)+(north-south)
     !! (2) Changed adjustment to vels to reflect (1) 
     !!------------------------------------------------------------------------------

!!DB
    INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
    REAL(wp):: transE,transW,transS,transN,areaE,areaW,areaS,areaN

    REAL(wp):: sum_obc, area, resid  ! temporary var.
    INTEGER :: ji,jj,jk,ij,ii
    COMPLEX(wp) :: cobcd,carea
    INTEGER :: itype
    cobcd = cmplx(0.e0, 0.e0, wp)
    carea = cmplx(0.e0, 0.e0, wp)

    sum_obc=0.0
    transE=0.0;    transW=0.0;    transS=0.0;    transN=0.0;
    areaE=0.0;     areaW=0.0;     areaS=0.0;     areaN=0.0; 

    IF( lp_obc_east  ) THEN 
       DO ji = nie0, nie1 
          DO  jk = 1, jpkm1
             DO jj = nje0p1, nje1m1
                ij = jj -1 + njmpp
                CALL DDPDD( cmplx(ufoe(jj,jk)*e2t(ji,jj)*fse3t(ji,jj,jk)*uemsk(jj,jk)*tmask_i(ji+1,jj) ,0.e0, wp), cobcd, 1, itype)
                CALL DDPDD( cmplx(e2t(ji,jj)*fse3t(ji,jj,jk)*uemsk(jj,jk)*tmask_i(ji+1,jj) ,0.e0, wp), carea, 1, itype) 
             END DO
          END DO
       END DO

    ENDIF

    IF( lp_obc_west  )  THEN 
       DO ji=niw0,niw1
          DO jk = 1, jpkm1
             DO jj = njw0p1, njw1m1
                ij = jj -1 + njmpp
                CALL DDPDD( cmplx(-ufow(jj,jk)*e2t(ji,jj)*fse3t(ji,jj,jk)*uwmsk(jj,jk)*tmask_i(ji,jj) ,0.e0, wp), cobcd, 1, itype)
                CALL DDPDD( cmplx(e2t(ji,jj)*fse3t(ji,jj,jk)*uwmsk(jj,jk)*tmask_i(ji,jj),0.e0, wp), carea, 1, itype) 
             END DO
          END DO
       END DO

    ENDIF
    
    IF( lp_obc_north )  THEN 
       DO jj = njn0,njn1
          DO jk = 1, jpkm1
             DO ji = nin0p1, nin1m1
                ii = ji -1 + nimpp
                CALL DDPDD( cmplx(vfon(ji,jk)*e2t(ji,jj)*fse3t(ji,jj,jk)*vnmsk(ji,jk) ,0.e0, wp), cobcd, 1, itype)
                CALL DDPDD( cmplx(e2t(ji,jj)*fse3t(ji,jj,jk)*vnmsk(ji,jk) ,0.e0, wp), carea, 1, itype)              
             END DO
          END DO
       END DO
       
    ENDIF
    
    IF( lp_obc_south )  THEN 
       DO  jj = njs0,njs1
          DO jk = 1, jpkm1
             DO ji = nis0p1, nis1m1
                ii = ji -1 + nimpp
                CALL DDPDD( cmplx(-vfos(ji,jk)*e2t(ji,jj)*fse3t(ji,jj,jk)*vsmsk(ji,jk)* tmask_i(ji,jj) ,0.e0, wp), cobcd, 1, itype)
                CALL DDPDD( cmplx(e2t(ji,jj)*fse3t(ji,jj,jk)*vsmsk(ji,jk)* tmask_i(ji,jj),0.e0, wp), carea, 1, itype)              

             END DO
          END DO
       END DO

    ENDIF


    IF(lk_mpp) call  mpp_sum(cobcd)
    IF(lk_mpp) call  mpp_sum(carea)
    
    sum_obc=real(cobcd)
    area=real(carea)
    
    resid=sum_obc/area

!!DBG: Output once-per-day
    if(lwp .AND. mod(kt-nit000,int(rday/rdt)) == 0) then
       write(numout2,'(A35,2x,i10,2x,3(e12.6,1x))')  &
            'DB -- obc_vol_ctl: kt, residual: ',kt,resid
    endif
    
!!DB: changed adjustments to reflect div calc
    IF( lp_obc_east  ) THEN 
       DO ji = nie0, nie1 
          DO  jk = 1, jpkm1
             DO jj = nje0p1, nje1m1
                ij = jj -1 + njmpp
                 ufoe(jj,jk)=(ufoe(jj,jk)-resid)*uemsk(jj,jk) 
             END DO
          END DO

       END DO

    ENDIF

    IF( lp_obc_west  )  THEN 
       DO ji=niw0,niw1
          DO jk = 1, jpkm1
             DO jj = njw0p1, njw1m1
                ij = jj -1 + njmpp
                 ufow(jj,jk)=(ufow(jj,jk)+resid)*uwmsk(jj,jk)
             END DO
          END DO

       END DO
    ENDIF

    IF( lp_obc_north )  THEN 
       DO jj = njn0,njn1
          DO jk = 1, jpkm1
             DO ji = nin0p1, nin1m1
                ii = ji -1 + nimpp
                 vfon(ji,jk)=(vfon(ji,jk)-resid)*vnmsk(ji,jk)
             END DO
          END DO
       END DO

    ENDIF

    IF( lp_obc_south )  THEN 
       DO  jj = njs0,njs1
          DO jk = 1, jpkm1
             DO ji = nis0p1, nis1m1
                ii = ji -1 + nimpp
                 vfos(ji,jk)=(vfos(ji,jk)+resid)*vsmsk(ji,jk)
             END DO
          END DO

       END DO
    ENDIF


     
   END SUBROUTINE obc_ctl
   
   
#else
   !!=================================================================================
   !!                       ***  MODULE  obcfla  ***
   !! Ocean dynamics:   Flather's algorithm at open boundaries for the time-splitting
   !!=================================================================================
 CONTAINS
   
   SUBROUTINE obc_ctl
   END SUBROUTINE obc_ctl
#endif
   
 END MODULE obcctl
 
