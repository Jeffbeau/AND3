!!DB: 2009.08.17
!!Eliminated some routines that are no longer necessary. If needed look in
!!bgcm_01_model.F91 (if it exists locally) or some achived version of the code
!!NB: RE Open Boundary Conditions:
!!NB: trc_trp() calls trc_nxt() which calls update_boundary_vals()
!!which is the correct location for updating the tracer arrays with OB values


!!DB: 2009.06.09
!!Latest version which uses "copy-to-boundary" OBC scheme for PZ and 
!!N_obc routine for N
!!Note that without further modifying the code structure (which is not
!!necessarily a bad idea) --
!!  All calls related to assignment of OBCs must be done from trcnxt.F90
!!Thus look there for OBC scheme currently used


!!OLD comments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!ca. 2009.04
!!It also contains OBC-related fixes relative to previous versions
!!(see also TRP/trcnxt.F90 for OBC-related change)
!!DB 2008.08.26 ...
!!Replace trcstp with 1 simpler module specific to key_BGCM_01
!!See TRP/trcstp.F90 to see how this module is used (i.e. how some
!!backward compatibility is maintained)

!!The BioGeoChemical Model works by (1) calling a Source-Minus-Sink (SMS)
!!routine that does the vertical biophysics only; and (2) calling the
!!regular tracer advection/diffusion schemes. The tracer array(s) default 
!!to trn(:,:,:,jptra) (see trc.F90 and par_trc ----> par_trc_trp). 
!!jptra is the number of 3D BGCM fields used. 


!!Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!-1- The main BGCM is in bgcm_sms() -- at bottom of this module. It is currently
!!    v.simple. This is where the major modifications need to be done.
!!-2- The SMS routine would get its various model and runtime parameters in the module
!!    bgcm_01_initrc.F90. At this time there is only a place-holder for this
!!    code in that module. 
!!END Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!END OLD comments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined key_passivetrc
MODULE bgcm_01_model
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================


   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE trc              ! ocean passive tracers variables 
   USE trctrp           ! passive tracers transport

   USE lib_bgcm_01
   USE lib_ncdf
!!DB
   USE oce, ONLY : sw_rad

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_stp           ! called by step

!!DB
   PUBLIC bgcm_sms

!!DB: bgcm_01.inc file variables used by step_npz_1D
  INTEGER :: last_lev
  REAL, DIMENSION(jpk) :: NN,PP,ZZ, NN0,PP0,ZZ0
  REAL(wp), DIMENSION(jpk+1) ::  wmask_1D


CONTAINS

   SUBROUTINE trc_stp( kt, kindic )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!
      !! History :
      !!   9.0  !  04-03  (C. Ethe)  Original
      !!-------------------------------------------------------------------

     USE lbclnk

      !! * Arguments
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      INTEGER, INTENT( in ) ::  kindic
      CHARACTER (len=25) :: charout

!!DB
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zfield  
      integer :: ji,jj,jk,jn, rec_num, status
      CHARACTER (len=80) :: fname
      REAL(wp), DIMENSION(100) :: bgcm_params
      REAL(wp) :: tmp1


!!DB -- 2009.10.07 -- write params to file (HARDWIRED):
      if(lwp .AND. kt==nit000) then
!!output to numout2
         write(numout2,*)'----------------------'
         write(numout2,*)'BGCM Parameters '
         write(numout2,'(a20,2x,e12.6)')'Vm = ',Vm
         write(numout2,'(a20,2x,e12.6)')'Rm = ',Rm
         write(numout2,'(a20,2x,e12.6)')'ks = ',ks
         write(numout2,'(a20,2x,e12.6)')'lambda = ',lambda
         write(numout2,'(a20,2x,e12.6)')'gamma = ',gamma
         write(numout2,'(a20,2x,e12.6)')'M_ZOO = ', M_ZOO
         write(numout2,'(a20,2x,e12.6)')'M_PHY = ', M_PHY
         write(numout2,'(a20,2x,e12.6)')'kext = ', kext
         write(numout2,'(a20,2x,e12.6)')'aaa = ', aaa
         write(numout2,'(a20,2x,e12.6)')'PAR = ', PAR
         write(numout2,'(a20,2x,e12.6)')'w_sink = ', w_sink
         write(numout2,*)'----------------------'
      endif


!!DB 2009.10.23 -- need sw_rad variable qsr() which is computed elsewhere
!!(method depends on the KEY). A problem arises in that tra_sbc() zeros it depending
!!on penetrative radiation flag = ln_trasqr. To guarantee that this routine has this
!!variable I search for a non-zero value and only assign sw_rad = qsr if qsr /= 0.
!!Note that I do this at every kt as qsr is not necessarily computed every kt so by
!!bad luck if I did this after the "return" statement below I may never get non-zero values.

!!NB: qsr() may not exist at all timesteps so only assign sw_rad when it does
!!NB: I use HARDWIRED position (10,10) which should always be OK 
    tmp1 = 0.0
    tmp1 = qsr(10,10)
    call mpp_sum(tmp1)
    if(tmp1 /= 0.0) then
       sw_rad(:,:) = qsr(:,:)
!Limit notification to 1/day if this works; if not it's no big deal
       if(lwp .AND. mod(kt-nit000,int(rday/rdt)) ==0) then
          write(numout2,*)'BGCM Updating sw_rad  at kt = ', kt
       endif
    endif

      ! this ROUTINE is called only every ndttrc time step
      IF( MOD( kt , ndttrc ) /= 0 ) RETURN

!!This is the routine that calls the BGCM Source-Minus-Sink code
      CALL bgcm_sms( kt )

      ! transport of passive tracers
!!DB: No need to change this; see trctrp.F90
!!NB: trc_trp() calls trc_nxt() which calls update_boundary_vals()
!!which is the correct location for updating the tracer arrays with OB values
      CALL trc_trp( kt )

!!DB: Open Boundary code 
!!This is where special OBC-related operations could also be done
#ifdef key_obc

#endif


!!DB: write restart -- not yet written
!!Restart file will be v.similar to regular output file 
!!the last frame of which can be used for the time-being. 
!      CALL trc_wri( kt )            ! outputs

!!DB: output model variables
      if(mod( kt, nwritetrc ) == 0 ) then
!!DB: write to ncdf file
         rec_num = (kt-nit000+1)/nwritetrc
         CALL ncdf_write(bgcm_fname, 'time_counter', REAL(kt * rdt), rec_num, status)

!!DB: 2008.10.17: Both of the below work 
!!DB: I use the NPZ version but keep the trn setup in the netcdf file in case
!!DB: it is desired to use it. (See lib_bgcm_01.F90 variable jf) 
!!DB: Note that to force a write of a 2D,3D,4D variable directly to a record number
!!    use a -ve value for rec_num

!         CALL ncdf_write(bgcm_fname, 'trn', trn, jptra, -rec_num, status)
!!DB: or
         zfield(:,:,:) = trn(:,:,:,1)
         CALL ncdf_write(bgcm_fname, 'N', zfield, -rec_num, status)
         zfield(:,:,:) = trn(:,:,:,2)
         CALL ncdf_write(bgcm_fname, 'P', zfield, -rec_num, status)
         zfield(:,:,:) = trn(:,:,:,3)
         CALL ncdf_write(bgcm_fname, 'Z', zfield, -rec_num, status)
!!AD/DB 2009.09.30
         CALL ncdf_write(bgcm_fname, 'ndastp',REAL(ndastp), rec_num, status)
         CALL ncdf_write(bgcm_fname, 'model_time_step',REAL(kt), rec_num, status)
         CALL ncdf_write(bgcm_fname, 'model_time',model_time, rec_num, status)

      endif

!!DB: tra is done in trc_trp ---> trc_nxt. 
!!To be safe do all of them (as done in p4zprg)
!!NB: Likely do not have to do this
      DO jn=1 , jptra
        CALL lbc_lnk(trn(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk(trb(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk(tra(:,:,:,jn), 'T', 1. )
      END DO



   END SUBROUTINE trc_stp

!!DB: BGCM_01 SMS model
   SUBROUTINE bgcm_sms( kt )
      !!===========================================================================================

      !! * Arguments
      !! -----------
     INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
     
!!DB
     real(wp) :: dk_tau
!!DB -- example of using pointers to make NPZ code more explicit
     real(wp), pointer :: N(:,:,:), P(:,:,:), Z(:,:,:)
     real(wp), pointer :: N0(:,:,:), P0(:,:,:), Z0(:,:,:)
     

!!DB -- assign pointers to default trn, trb arrays
     N => trn(:,:,:,1)
     P => trn(:,:,:,2)
     Z => trn(:,:,:,3)
     N0 => trb(:,:,:,1)
     P0 => trb(:,:,:,2)
     Z0 => trb(:,:,:,3)
     
     !! this ROUTINE is called only every ndttrc time step
     !! --------------------------------------------------
     IF ( MOD(kt,ndttrc) /= 0) RETURN
     
     call step_npz(kt)
     
   END SUBROUTINE bgcm_sms



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Simple P Franks routine, working on *2 versions of variables
! Call this after calls to STEP_TRACER which handles the adv/diff part
! of field updates. 
!
! REM: I0(d) = 125 + 100*cos(2*pi/360 * (day-180))
!            = fit by DB to radiation at ~ 45N
!
!REM: These model parameters are set at top of module
!    Vm = 2.0/86400.0  !! convert per day to per second
!    Rm = 0.5/86400.0  
!    ks = 1.0          !! \mu Mol N per liter
!    lambda = 0.2      !! grazing
!    gamma = 0.7
!    w_sink = -1.5/86400.0 !!sinking rate in m/day
!    M_ZOO = 0.04/86400.0 !! Zoo mort
!    M_PHY = 0.05/86400.0   !! M_PHY
!    lambda = 0.3      !! grazing
!    kext = 0.1        !!in m
!    aaa = 0.025       !! initial slope of P-I curve
!    PAR = 0.43        !!from Oschlies

  SUBROUTINE step_npz(kt)

#if ! defined key_flx_bulk_daily
    USE flxblk_solar
#endif

!!DB
    USE oce, ONLY : sw_rad
!!DBG
!    USE lib_mpp

    INTEGER, INTENT( in ) ::   kt

    REAL(wp) :: DTIN, tmp2, tmp1, fac, fac0
    REAL(wp), DIMENSION(jpi,jpj) :: C2
    REAL(wp), DIMENSION(jpi,jpj,2) ::  PL
    integer  :: k,l,j,lt,lb,L0,i

!!DB -- example of using pointers to make NPZ code more explicit
    real(wp), pointer :: N(:,:,:), P(:,:,:), Z(:,:,:)
    real(wp), pointer :: N0(:,:,:), P0(:,:,:), Z0(:,:,:)

!!DB -- assign pointers to default trn, trb arrays
    N => trn(:,:,:,1)
    P => trn(:,:,:,2)
    Z => trn(:,:,:,3)
    N0 => trb(:,:,:,1)
    P0 => trb(:,:,:,2)
    Z0 => trb(:,:,:,3)
    
    !===================================
    ! NPZ interaction equations
    !===================================

!!DB: Kept for reference purposes; should never be used anymore
!!If not calculating solar forcing elsewhere then must use the next lines
!    I_sfce = 125.0 + 100.0*cos(2.*rpi/360.0 * ((nday_year)-180.0)) !! W/m^2
!    I_sfce = PAR * aaa * I_sfce
!    C2(:,:) = (Vm*86400.0/I_sfce)**2    
!!Also must comment-out the I_sfce line below

#if ! defined key_flx_bulk_daily 
!!DB -- Use AD solar flux ---> ad_qsr_oce(i,j) 
   call flx_blk_solar(kt)
    sw_rad(:,:) = ad_qsr_oce(:,:)
#endif

!!     Do sfce layer first to simplify calc of sinking term: w_sink * dP/dz
!!NB C2(i,j) is computed here so it is not necessary to do so elsewhere
    lt = 1
    lb = 2
    do k = 1, 1
       L0 = k + 1
       do j = 1, jpj        
          do i = 1, jpi
!!DB 2009.10.23 -- below is the desired line-of-code
!!  However, sw_rad can contain ~zeros and depending on the machine this can result
!!  in a massive slowdown of execution. So I add 0.1 to sw_rad to avoid this problem.
!!  As of the above date, this seems to work
!             I_sfce = PAR * aaa * sw_rad(i,j)  !current desired line-of-code
             I_sfce = PAR * aaa * (sw_rad(i,j)+0.1)  !DBG
             C2(i,j) = (Vm*86400.0/I_sfce)**2    
             
             DTIN =  rdt * tmask(i,j,k)
             !                
             PL(i,j,lt) = 0.0
!!     Upwind diff
                fac0 = sign(1.0,w_sink)
                PL(i,j,lb) = ((1.+fac0)*P0(i,j,L0) +  &
                     (1.-fac0)*P0(i,j,k))/2.
!                
!!      this is z-integral of actual f(z,t); 
                Ibar = 0.0 +            &
                     (1.0/e3t(k)) * 1.0/kext *                           &
                                      (                                  &
                                      log( (exp(-kext*gdept(k))+         &
                                      sqrt(C2(i,j)+exp(-2.*kext*gdept(k)))) / &
                                      (exp(-kext*gdept(k+1))+            &
                                      sqrt(C2(i,j)+exp(-2.*kext*gdept(k+1))))  ) &
                                      ) 
!                
                P(i,j,k) = P(i,j,k) + DTIN * (          &
                     Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)    & 
                     - Z0(i,j,k)*Rm*(1.0-exp(-lambda*P0(i,j,k)))   & 
                     - M_PHY*P0(i,j,k)                  & 
!!DB: use  wmask here (computed in lib_bgcm_01)
                  -(wmask(i,j,k)*w_sink*PL(i,j,lt)-     &
                  wmask(i,j,L0)*w_sink*PL(i,j,lb))/e3t(k) )  !!sinking term
!                
                Z(i,j,k) = Z(i,j,k) + DTIN * (           &
                     gamma*Rm*Z0(i,j,k)*(1.0-exp(-lambda*P0(i,j,k)))        &
                     - M_ZOO*Z0(i,j,k) )
!                
                N(i,j,k) = N(i,j,k) + DTIN * (            &
                     -Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)        &
                     +(1.0-gamma)*Rm*Z0(i,j,k)*           &
                     (1.0-exp(-lambda*P0(i,j,k)))         &
                     + M_PHY*P0(i,j,k) + M_ZOO*Z0(i,j,k) )
!                
          enddo
       enddo
    enddo
    
!    
!!DB: as w_sink is always of the same sign (-ve) compute fac0 outside of loop
!!Use wmask to stop flux through bottom
    fac0 = sign(1.0,w_sink)
    lt = 2
    lb = 1
    do k = 2, jpkm1
       L0 = k + 1
       do j = 1, jpj
          do i = 1, jpi

             DTIN =  rdt * tmask(i,j,k)
!!     Upwind diff
             PL(i,j,lb) = ((1.+fac0)*P0(i,j,L0) +  &
                  (1.-fac0)*P0(i,j,k))/2.
!                
!!      this is z-integral of actual f(z,t);
             Ibar = 0.0 +            &
                  (1.0/e3t(k)) * 1.0/kext *                           &
                  (                                  &
                  log( (exp(-kext*gdept(k))+         &
                  sqrt(C2(i,j)+exp(-2.*kext*gdept(k)))) / &
                  (exp(-kext*gdept(k+1))+            &
                  sqrt(C2(i,j)+exp(-2.*kext*gdept(k+1))))  ) &
                  ) 
!                
             P(i,j,k) = P(i,j,k) + DTIN * (          &
                  Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)    & 
                  - Z0(i,j,k)*Rm*(1.0-exp(-lambda*P0(i,j,k)))   & 
                  - M_PHY*P0(i,j,k)                  & 
!!DB: use  wmask here (computed in lib_bgcm_01)
                  -(wmask(i,j,k)*w_sink*PL(i,j,lt)-     &
                  wmask(i,j,L0)*w_sink*PL(i,j,lb))/e3t(k) )  !!sinking term

             Z(i,j,k) = Z(i,j,k) + DTIN * (           &
                  gamma*Rm*Z0(i,j,k)*(1.0-exp(-lambda*P0(i,j,k)))        &
                  - M_ZOO*Z0(i,j,k) )
!                
             N(i,j,k) = N(i,j,k) + DTIN * (            &
                  -Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)        &
                  +(1.0-gamma)*Rm*Z0(i,j,k)*           &
                  (1.0-exp(-lambda*P0(i,j,k)))         &
                  + M_PHY*P0(i,j,k) + M_ZOO*Z0(i,j,k) )
             !                
          enddo
       enddo
    
       i = lt          !!switch indices
       lt = lb
       lb = i
       
    enddo       !! k-loop

!!DB: special loop for jpk layer (in case it exists)
!!NB: w_sink code frag
    do k = jpk, jpk
       L0 = k + 1
       do j = 1, jpj
          do i = 1, jpi

             DTIN =  rdt * tmask(i,j,k)
!                
!!      this is z-integral of actual f(z,t); calc here as it will become f(x,y,z,t)
             Ibar = 0.0 +            &
                  (1.0/e3t(k)) * 1.0/kext *                           &
                  (                                  &
                  log( (exp(-kext*gdept(k))+         &
                  sqrt(C2(i,j)+exp(-2.*kext*gdept(k)))) / &
                  (exp(-kext*gdept(k+1))+            &
                  sqrt(C2(i,j)+exp(-2.*kext*gdept(k+1))))  ) &
                  ) 
!                
             P(i,j,k) = P(i,j,k) + DTIN * (          &
                  Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)    & 
                  - Z0(i,j,k)*Rm*(1.0-exp(-lambda*P0(i,j,k)))   & 
                  - M_PHY*P0(i,j,k)                  & 
!!DB: NEW  use wmask here (computed in lib_bgcm_01)
                  -(wmask(i,j,k)*w_sink*P0(i,j,k-1)-     &
                  0.0)/e3t(k) )  !!sinking term for bottom layer

             Z(i,j,k) = Z(i,j,k) + DTIN * (           &
                  gamma*Rm*Z0(i,j,k)*(1.0-exp(-lambda*P0(i,j,k)))        &
                  - M_ZOO*Z0(i,j,k) )
!                

             N(i,j,k) = N(i,j,k) + DTIN * (            &
                  -Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)        &
                  +(1.0-gamma)*Rm*Z0(i,j,k)*           &
                  (1.0-exp(-lambda*P0(i,j,k)))         &
                  + M_PHY*P0(i,j,k) + M_ZOO*Z0(i,j,k) )
             !                
          enddo
       enddo
    
       i = lt          !!switch indices
       lt = lb
       lb = i
       
    enddo       !! k-loop


!!DB: 2009.04.27 add limiter to N,P,Z
    P(:,:,:) = max(P(:,:,:),1.e-4)
    Z(:,:,:) = max(Z(:,:,:),1.e-4)
    N(:,:,:) = max(N(:,:,:),1.e-4)

    
  END SUBROUTINE step_npz


END MODULE bgcm_01_model

#endif
