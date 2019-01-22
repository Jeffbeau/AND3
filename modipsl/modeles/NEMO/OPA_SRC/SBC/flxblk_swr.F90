!! Nicolas Lambert, feb. 2014


MODULE flxblk_swr
   !!======================================================================
   !!                       ***  MODULE  flxblk_swr  ***
   !!
   !!=================== SUMMARY =========================================
   !!Purpose:         Stand alone Short wave flux calcuation at kt 
   !!
   !!Inputs:         None / CMC humidity cloud cover and air temp
   !!
   !!Outputs:             qsr_ocep,     &        !: solar flux over ocean (at kt)
   !! 
   !! Usage:   -compilie with key_bulk_solar. No other keys required
   !!          -SUBROUTINE flx_blk_swr( kt) must be called each hour (for now in likemaker)
   !!          -Access qsr_ocep with "USE flxblk_solar"
   !!================= DESCRIPTION =======================================
   !!
   !!    SUBROUTINE flx_blk_swr( kt) computes shortwave radiation based
   !!     on the Zillman (1972) by computing the hours of the day. It using
   !!     the atmospheric data loaded by flx_blk (data in blk_oce) from CMC 
   !!     or NCEP (or others sources).
   !!   =====================================================================

#if defined key_BGCM_02
    USE in_out_manager
    USE daymod, only: nmonth,nday,nyear,nday_year !AD: to get opa calendar
    USE phycst, only:  rpi, rt0,rad ,raajj   ! physical constants (pi, zeroK,rad,ndaysyear !NL#8
    USE blk_oce, only: &     ! atmospheric parameters
         catm     ,      &    ! fraction of cloud
         tatm     ,      &    !: atmospheric temperature
         hatm                 !: relative humidity
    USE daymod, only: model_time   ! date : yyyymmdd.ddd
    USE oce_trc, only: gphit       !Latitude
    USE par_oce, only:jpi, jpj     ! domain dimension
    USE albedo, only:flx_blk_albedo ! to get the albedo


    IMPLICIT NONE
     PRIVATE
      
     !! * Accessibility
     PUBLIC flx_blk_swr 

     real, public :: qsr_ocep(jpi,jpj)          !short wave radiation NL#8


CONTAINS

      SUBROUTINE flx_blk_swr( kt )
      !NL#8 feb. 2014
      ! Calcul the short wave radiation (qsr_oce) for the Par_surf
      ! kt : time step (input) 
      ! qsr_oce_bio : short wave radiation (output)


      implicit none
      integer kt
      real net(jpi,jpj)
      integer jj,ji
      real a0,a1,a2,a3,a4,b1,b2,b3,b4
      real zp1,zdecl,zdaycor,jth,hour_GMT
      real zps,zpc,zlsrise,zlmunoon,zcldcor,zlha,zcmue,zsqsro,zev
      integer sign_zps,sign_ztamr
      REAL(wp), DIMENSION(jpi,jpj) ::   & 
         zalbocsd         ,  &   ! albedo of ocean
         zalboos          ,  &   ! albedo of ocean under overcast sky
         zalbics          ,  &   ! albedo of ice under clear sky
         zalbios          ,  &   ! albedo of ice under overcast sky
         zalbomu                 ! albedo of ocean when zcmue is 0.4
      REAL(wp)  ::            &
         zalbo             ,  &  ! albedo of sea-water
         zalbocs

      ! hours form the day fraction
      jth= mod(model_time,1.) * 24.

      !if(lwp) print*,'Calcul in swr_calcul (kt,hours),',kt,jth

      ! coefficients et constantes
      a0  =  0.39507671  
      a1  = 22.85684301  
      a2  = -0.38637317  
      a3  =  0.15096535  
      a4  = -0.00961411  
      b1  = -4.29692073  
      b2  =  0.05702074  
      b3  = -0.09028607  
      b4  =  0.00592797  
       
      zp1 = rpi * ( 2.0 * nday_year - 367.0 ) / raajj

      zdecl  = a0 + a1 * cos( zp1    ) + a2 * cos( 2.*zp1 ) &
                  + a3 * cos( 3.*zp1 ) + a4 * cos( 4.*zp1 ) &
                  + b1 * sin( zp1    ) + b2 * sin( 2.*zp1 ) &
                  + b3 * sin( 3.*zp1 ) + b4 * sin( 4.*zp1 )

      zdaycor = 1.0 + 0.0013 * sin(nday_year * 2.*rpi / raajj) & 
                    + 0.0342 * cos(nday_year * 2.*rpi / raajj)

      hour_GMT=-5.

      zlha=cos((12-jth-hour_GMT)*3.1416/12.)  ! (-5 GMT)

      ! call the get albedo
      CALL flx_blk_albedo( zalbios, zalboos, zalbics, zalbomu )
     
      do jj=1,jpj
          do ji=1,jpi
              
              zps  = sin(gphit(ji,jj)*rad) * sin(zdecl * rad)
              zpc  = cos(gphit(ji,jj)*rad) * cos(zdecl * rad)
              if (zps.lt.0.) then;sign_zps=-1;else;sign_zps=1;endif 
              zlsrise  = acos (-sign_zps * min(1.,sign_zps * zps/zpc ))
              
              zlmunoon = asin(( zps + zpc)) / rad
              zcldcor=min(1.,(1. - 0.62 *catm(ji,jj) + 0.0019 * zlmunoon))
              
              if (tatm(ji,jj).lt. 0) then ;sign_ztamr=-1;else;sign_ztamr=1;endif 
 
              zev = hatm(ji,jj) *611.0 * exp( abs(tatm(ji,jj)) * min(17.269 * sign_ztamr,21.875 * sign_ztamr) &
                    / (tatm(ji,jj) + rt0 - 35.86 + max (0.,28.200 * -sign_ztamr)))
              
              zcmue = max (0.,zps + zpc * zlha)
              zsqsro =   (1368.0 *zcmue**2) /((zcmue +2.7) * zev* 1.0e-05 + 1.085 * zcmue +0.10)

              ! calcul of albedo
              zalbocs  = 1.2 * 0.05 / ( 1.1 * zcmue**1.4 + 0.15 )
              zalbo = ( 1.0 - catm(ji,jj) ) * zalbocs + catm(ji,jj) * zalboos(ji,jj)
              
              net(ji,jj)=zsqsro * 0.9 * zcldcor* (1-zalbo)*zdaycor;        
          enddo
      enddo

      qsr_ocep = net

   END SUBROUTINE flx_blk_swr

#else
   !!----------------------------------------------------------------------
   !!   Default option :           Empty module                     NO bulk
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flx_blk_swr            ! Empty routine
   END SUBROUTINE flx_blk_swr
#endif


END MODULE flxblk_swr