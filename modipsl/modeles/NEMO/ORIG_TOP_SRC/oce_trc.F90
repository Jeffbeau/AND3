MODULE oce_trc
   !!======================================================================
   !!                      ***  MODULE  oce_trc  ***
   !! Ocean passive tracer  :  share ocean-passive tracers variables
   !!======================================================================
   !! History :
   !!   9.0  !  04-03  (C. Ethe)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!  TOP 1.0,  LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/oce_trc.F90,v 1.13 2006/04/10 15:40:28 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   !! Domain characteristics
   USE par_oce , ONLY :       &
      cp_cfg   =>   cp_cfg,  & !: name of the configuration
      jp_cfg   =>   jp_cfg,  & !: resolution of the configuration
      jpiglo   =>   jpiglo,  & !: first  dimension of global domain --> i
      jpjglo   =>   jpjglo,  & !: second dimension of global domain --> j
      jpi      =>   jpi   ,  & !: first  dimension of grid --> i 
      jpj      =>   jpj   ,  & !: second dimension of grid --> j  
      jpk      =>   jpk   ,  & !: number of levels  
      jpim1    =>   jpim1 ,  & !: jpi - 1
      jpjm1    =>   jpjm1 ,  & !: jpj - 1 
      jpkm1    =>   jpkm1 ,  & !: jpk - 1  
      jpij     =>   jpij  ,  & !: jpi x jpj
      jpidta   =>   jpidta,  & !: first horizontal dimension  > or = jpi
      jpjdta   =>   jpjdta,  & !: second horizontal dimension > or = jpj
      jpkdta   =>   jpkdta,  & !: number of levels            > or = jpk
      lk_esopa =>   lk_esopa   !: flag to activate the all option


   !! namelist parameters      
   USE in_out_manager , ONLY :  &    
      cexper   =>   cexper,  & !: experience name for vairmer format 
      no       =>   no    ,  & !: job number
      nrstdt   =>   nrstdt,  & !: control of the time step (0,  & 1 or 2)
      nit000   =>   nit000,  & !: number of the first time step
      nitend   =>   nitend,  & !: number of the last time step
      nleapy   =>   nleapy,  & !: Leap year calendar (0/1)
      nwrite   =>   nwrite,  & !: frequency of OUTPUT file
      nstock   =>   nstock,  & !: frequency of restart file
      nprint   =>   nprint,  & !: level of print (0 no print)
      lwp      =>   lwp   ,  & !: boolean term for mpp output
      ln_ctl   =>   ln_ctl,  & !: = ln_ctl.AND.lwp (print control on the 1st proc)
      nictls   =>   nictls,  & !: Start i indice for the SUM control
      nictle   =>   nictle,  & !: End   i indice for the SUM control
      njctls   =>   njctls,  & !: Start j indice for the SUM control
      njctle   =>   njctle,  & !: End   j indice for the SUM control
      isplt    =>   isplt ,  & !: number of processors following i
      jsplt    =>   jsplt ,  & !: number of processors following j
      ijsplt   =>   ijsplt,  & !: nb of local domain = nb of processors
      numout   =>   numout     !: logical unit for output print


   !! run controm   
   USE in_out_manager , ONLY :  &  
      nstop     =>   nstop    ,  &  !: e r r o r  flag (=number of reason for a
      !                             !                   prematurely stop the run)
      nwarn     =>   nwarn    ,  &  !: w a r n i n g  flag (=number of warning
      !                             !                   found during the run) 
      cform_err =>   cform_err,  &  !:
      cform_war =>   cform_war      !:  
      
   USE dom_oce , ONLY :           &            
      lzoom      => lzoom     ,  & !: zoom flag
      lzoom_e    => lzoom_e   ,  & !: East  zoom type flag
      lzoom_w    => lzoom_w   ,  & !: West  zoom type flag
      lzoom_s    => lzoom_s   ,  & !: South zoom type flag
      lzoom_n    => lzoom_n   ,  & !: North zoom type flag
      lzoom_arct => lzoom_arct,  & !: ORCA    arctic zoom flag
      lzoom_anta => lzoom_anta     !: ORCA antarctic zoom flag



   USE dom_oce , ONLY :       & 
      nperio   =>   nperio,  & !: type of lateral boundary condition       
!      nlci     =>   nlci  ,  & !: index i for the sub domain left bottom 
!      nlcj     =>   nlcj  ,  & !: index j for the sub domain left bottom 
      nimpp    =>   nimpp ,  & !: i index for mpp-subdomain left bottom
      njmpp    =>   njmpp ,  & !: j index for mpp-subdomain left bottom
      nproc    =>   nproc ,  & !: number for local processor
      narea    =>   narea ,  & !: number for local area
      mig      =>   mig   ,  & !: local  ==> global  domain i-indice
      mjg      =>   mjg   ,  & !: local  ==> global  domain i-indice
      mi0      =>   mi0   ,  & !: global ==> local domain i-indice 
      mi1      =>   mi1   ,  & !: (mi0=1 and mi1=0 if the global indice is not in the local domain)
      mj0      =>   mj0   ,  & !: global ==> local domain j-indice 
      mj1      =>   mj1   ,  & !: (mj0=1 and mj1=0 if the global indice is not in the local domain)
      nidom    =>   nidom
 
   USE dom_oce , ONLY :       & 
      nimppt   => nimppt  ,  & !:i-indexes for each processor
      njmppt   => njmppt  ,  & !:j-indexes for each processor
      ibonit   => ibonit  ,  & !:i-processor neighbour existence
      ibonjt   => ibonjt  ,  & !:j- processor neighbour existence 
      nlci     => nlci    ,  & !:i- & j-dimensions of the local subdomain
      nlcj     => nlcj    ,  & !:
      nldi     => nldi    ,  & !:first and last indoor i- and j-indexes
      nlei     => nlei    ,  & !:
      nldj     => nldj    ,  & !:
      nlej     => nlej    ,  & !:
      nlcit    => nlcit   ,  & !:dimensions of every i-subdomain
      nlcjt    => nlcjt   ,  & !:dimensions of every j-subdomain
      nldit    => nldit   ,  & !:first indoor index for each i-domain 
      nleit    => nleit   ,  & !:last indoor index for each i-domain 
      nldjt    => nldjt   ,  & !:first indoor index for each j-domain 
      nlejt    => nlejt        !:last indoor index for each j-domain 

    
      !! horizontal curvilinear coordinate and scale factors
   USE dom_oce , ONLY :            &    
      glamt    =>   glamt ,  & !: longitude of t-point (degre)  
      glamu    =>   glamu ,  & !: longitude of t-point (degre)  
      glamv    =>   glamv ,  & !: longitude of t-point (degre)  
      glamf    =>   glamf ,  & !: longitude of t-point (degre)  
      gphit    =>   gphit ,  & !: latitude  of t-point (degre)   
      gphiu    =>   gphiu ,  & !: latitude  of t-point (degre)   
      gphiv    =>   gphiv ,  & !: latitude  of t-point (degre)   
      gphif    =>   gphif ,  & !: latitude  of t-point (degre)   
      e1t      =>   e1t   ,  & !: horizontal scale factors at t-point (m)  
      e2t      =>   e2t   ,  & !: horizontal scale factors at t-point (m)   
      e1u      =>   e1u   ,  & !: horizontal scale factors at u-point (m)
      e2u      =>   e2u   ,  & !: horizontal scale factors at u-point (m)
      e1v      =>   e1v   ,  & !: horizontal scale factors at v-point (m)
      e2v      =>   e2v        !: horizontal scale factors at v-point (m)  

   !! vertical coordinate and scale factors
   USE dom_oce , ONLY :              &   
      gdept    =>   gdept ,  & !: reference depth of t-points (m)
      e3t      =>   e3t   ,  & !: reference depth of t-points (m)  
      e3w      =>   e3w   ,  & !: reference depth of w-points (m)
      gdepw    =>   gdepw      !: reference depth of w-points (m)

   USE dom_oce ,   ONLY :            &      
      lk_zps   =>  lk_zps ,  & !: partial steps flag
      lk_sco   =>  lk_sco ,  & !: s-coordinate flag
      lk_zco   =>  lk_zco      !: z-coordinate flag

   USE lib_mpp ,   ONLY :            &     
      lk_mpp   =>  lk_mpp      !: Mpp flag

   USE dynspg_oce ,   ONLY :            &     
      lk_dynspg_rl   =>  lk_dynspg_rl      !: rigid lid flag

#if defined key_partial_steps
   !! Partial steps ('key_partial_steps')
   !! -----------------------------------
   USE dom_oce , ONLY :                & 
      gdep3w   =>  gdep3w  ,  & !: ???
      gdept_ps =>  gdept_ps,  & !: depth of t-points (m)
      gdepw_ps =>  gdepw_ps,  & !: depth of t-points (m)
      e3t_ps   =>  e3t_ps  ,  & !: vertical scale factors at t-
      e3u_ps   =>  e3u_ps  ,  & !: vertical scale factors at u-
      e3v_ps   =>  e3v_ps  ,  & !: vertical scale factors v-
      e3w_ps   =>  e3w_ps  ,  & !: w-points (m)
      e3f_ps   =>  e3f_ps  ,  & !: f-points (m)
      e3uw_ps  =>  e3uw_ps ,  & !: uw-points (m)
      e3vw_ps  =>  e3vw_ps      !: vw-points (m)

   USE oce , ONLY :                &
      gtu   =>  gtu  ,  & !: t- horizontal gradient at u-
      gtv   =>  gtv       !: and v-points at bottom ocean level
#endif

#if defined key_s_coord
   USE dom_oce , ONLY :              &   
      hbatt   =>   hbatt  ,  & !: ocean depth at the vertical of  t-point (m)
      hbatu   =>   hbatu  ,  & !: ocean depth at the vertical of  u-point (m)
      hbatv   =>   hbatv  ,  & !: ocean depth at the vertical of w-point (m)
      gsigt   =>   gsigt  ,  & !: model level depth coefficient at t-,  & w-levelsvertical scale factors at u-
      gsigw   =>   gsigw  ,  & !: model level depth coefficient at t-,  & w-levelsvertical scale factors v-
      gsi3w   =>   gsi3w  ,  & !: model level depth coef at w-levels (defined as the sum of e3w)
      esigt   =>   esigt  ,  & !: vertical scale factor coef. at t-levels
      esigw   =>   esigw       !: vertical scale factor coef. at w-levels
#endif

   !! masks, bathymetry
   USE dom_oce , ONLY :             &    
      mbathy   =>   mbathy,  & !: number of ocean level (=0,  & 1, ... , jpk-1) 
      tmask_i  =>   tmask_i, & !: Interior mask at t-points
      tmask    =>   tmask ,  & !: land/ocean mask at t-points
      umask    =>   umask ,  & !: land/ocean mask at u-points   
      vmask    =>   vmask ,  & !: land/ocean mask at v-points 
      fmask    =>   fmask      !: land/ocean mask at f-points 

   USE dom_oce , ONLY :         &
      n_cla   =>   n_cla       !: flag (0/1) for cross land advection 

   !! time domain
   USE dom_oce , ONLY :                 &
      neuler   =>   neuler,  & !: restart euler forward option (0=Euler)
      rdt      =>   rdt   ,  & !: time step for the dynamics 
      atfp     =>   atfp  ,  & !: asselin time filter parameter
      atfp1    =>   atfp1 ,  & !: asselin time filter coeff. (atfp1= 1-2*atfp)
      rdttra   =>   rdttra     !: vertical profile of tracer time step

   USE daymod , ONLY :                 &
      ndastp    =>   ndastp,  &    !: time step date in year/month/day aammjj
      nday_year =>   nday_year, &  !: curent day counted from jan 1st of the current year
      nyear     =>   nyear,   &  !: Current year
      nmonth    =>   nmonth,  &  !: Current month
      nday      =>   nday        !: Current day

   !! physical constants
   USE phycst ,   ONLY :                &  
      ra       =>   ra    ,  & !: earth radius
      rpi      =>   rpi   ,  & !: pi
      rday     =>   rday  ,  & !: day
      rauw     =>   rauw  ,  & !: density of pure water kg/m3
      ro0cpr   =>   ro0cpr,  & !: = 1. / ( rau0 * rcp )
      rad      =>   rad   ,  & !: conversion coeff. from degre into radian
      raass    =>   raass ,  & !: number of seconds in one year
      rmoss    =>   rmoss ,  & !: number of seconds in one month
      rjjss    =>   rjjss      !: number of seconds in one day

   !! present fields (now)
   USE oce , ONLY :            &     
      ua      =>    ua    ,  & !: i-horizontal velocity (m s-1) 
      va      =>    va    ,  & !: j-horizontal velocity (m s-1)
      un      =>    un    ,  & !: i-horizontal velocity (m s-1) 
      vn      =>    vn    ,  & !: j-horizontal velocity (m s-1)
      wn      =>    wn    ,  & !: vertical velocity (m s-1)  
      tn      =>    tn    ,  & !: pot. temperature (celsius)
      sn      =>    sn    ,  & !: salinity (psu)
      rhop    =>    rhop  ,  & !: potential volumic mass (kg m-3) 
      rhd     =>    rhd        !: in situ density anomalie rhd=(rho-rau0)/rau0 (no units)

#if defined key_trc_diatrd
   USE oce , ONLY :          &
      hdivn   =>    hdivn      !: horizontal divergence (1/s)
#endif

#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
   !! wind speed
   USE blk_oce , ONLY :        &     
      vatm    =>    vatm       !: wind speed at sea surface (m s-1)
#endif

   !! wind speed
   USE taumod , ONLY :        &     
      taux    =>    taux ,  &  !: i-surface stress component
      tauy    =>    tauy       !: j-surface stress component

#if   defined key_trabbl_dif   ||   defined key_trabbl_adv
   USE trabbl , ONLY :           &      
      atrbbl   =>   atrbbl     !: lateral coeff. for bottom boundary layer scheme (m2/s)
#  if defined key_off_tra
   USE trabbl, ONLY :            &
      bblx   => bblx,       &
      bbly   => bbly
#  endif
#endif

   !! lateral diffusivity (tracers)
   USE ldftra_oce ,   ONLY :             &    
      aht0    =>   aht0  ,  &  !: horizontal eddy diffusivity for tracers (m2/s)
      ahtb0   =>   ahtb0 ,  &  !: background eddy diffusivity for isopycnal diff. (m2/s)
      ahtu    =>   ahtu  ,  &  !: lateral diffusivity coef. at u-points 
      ahtv    =>   ahtv  ,  &  !: lateral diffusivity coef. at v-points 
      ahtw    =>   ahtw  ,  &  !: lateral diffusivity coef. at w-points 
      ahtt    =>   ahtt  ,  &  !: lateral diffusivity coef. at t-points
      aeiv0   =>   aeiv0 ,  &  !: eddy induced velocity coefficient (m2/s) 
      aeiu    =>   aeiu  ,  &  !: eddy induced velocity coef. at u-points (m2/s)   
      aeiv    =>   aeiv  ,  &  !: eddy induced velocity coef. at v-points (m2/s) 
      aeiw    =>   aeiw        !: eddy induced velocity coef. at w-points (m2/s) 

   !! vertical diffusion
   USE zdf_oce , ONLY :      &    
      avt            =>   avt          ,  & !: vert. diffusivity coef. at w-point for temp  
      avt0           =>   avt0         ,  & !: vertical eddy diffusivity for tracers (m2/s)
#if ! defined key_off_tra
      l_trazdf_exp   =>   l_trazdf_exp ,  & !: explicit vertical diffusion scheme flag
#endif
      ln_zdfnpc      =>   ln_zdfnpc         !: convection: non-penetrative convection flag


#if defined key_zdfddm
   USE zdfddm , ONLY :             &     
      avs     =>    avs        !: salinity vertical diffusivity coeff. at w-point
#endif

   !! penetrative solar radiation
   USE traqsr , ONLY :            &      
      xsi1   =>   xsi1         !: first depth of extinction 

   !! surface fluxes
   USE ocesbc , ONLY :             &   
      qt      =>    qt    ,  & !: total surface heat flux (w m-2)   
      qsr     =>    qsr   ,  & !: penetrative solar radiation (w m-2)  
      emp     =>    emp   ,  & !: evaporation minus precipitation (kg m-2 s-2) 
      emps    =>    emps       !: evaporation minus precipitation (kg m-2 s-2)

   !! freezing area
   USE ocfzpt , ONLY :            &      
      freeze  =>    freeze,  & !: ice mask (0 or 1)  
      fzptn   =>    fzptn      !: now freezing temperature at ocean surface  


   !! mixing layer depth (turbocline)
   USE zdfmxl , ONLY :             &    
      hmld    =>   hmld   ,  & !: mixing layer depth (turbocline)
      hmlp    =>   hmlp   ,  & !: mixed layer depth  (rho=rho0+zdcrit) (m)
      hmlpt   =>   hmlpt       !: mixed layer depth at t-points (m)

   USE ldfslp , ONLY :              & 
      lk_ldfslp  =>  lk_ldfslp     !: slopes flag
#if   defined key_ldfslp
   !! direction of lateral diffusion (momentum  tracers) 
   USE ldfslp , ONLY :              & 
      uslp       =>   uslp    ,  & !: i-direction slope at u-, w-points
      vslp       =>   vslp    ,  & !: j-direction slope at v-, w-points
      wslpi      =>   wslpi   ,  & !: i-direction slope at u-, w-points
      wslpj      =>   wslpj        !: j-direction slope at v-, w-points
#endif

   !! ocean forcings runoff
   USE flxrnf , ONLY :              &   
      upsrnfh =>   upsrnfh ,  & !: mixed adv scheme in runoffs vicinity (hori.) 
      upsrnfz =>   upsrnfz ,  & !: mixed adv scheme in runoffs vicinity (vert.)
      upsadv  =>   upsadv       !: mixed adv scheme in straits vicinity (hori.)

END MODULE oce_trc
