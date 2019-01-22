


PGF90 (Version     10.5)          09/23/2010  09:49:16      page 1

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

Filename: oce_trc.F90

(    1) MODULE oce_trc
(    2)    !!======================================================================
(    3)    !!                      ***  MODULE  oce_trc  ***
(    4)    !! Ocean passive tracer  :  share ocean-passive tracers variables
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   9.0  !  04-03  (C. Ethe)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!  TOP 1.0,  LOCEAN-IPSL (2005)
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/oce_trc.F90,v 1.13 2006/04/10 15:40:28 opalod Exp $
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    !! Domain characteristics
(   15)    USE par_oce , ONLY :       &
(   16)       cp_cfg   =>   cp_cfg,  & !: name of the configuration
(   17)       jp_cfg   =>   jp_cfg,  & !: resolution of the configuration
(   18)       jpiglo   =>   jpiglo,  & !: first  dimension of global domain --> i
(   19)       jpjglo   =>   jpjglo,  & !: second dimension of global domain --> j
(   20)       jpi      =>   jpi   ,  & !: first  dimension of grid --> i 
(   21)       jpj      =>   jpj   ,  & !: second dimension of grid --> j  
(   22)       jpk      =>   jpk   ,  & !: number of levels  
(   23)       jpim1    =>   jpim1 ,  & !: jpi - 1
(   24)       jpjm1    =>   jpjm1 ,  & !: jpj - 1 
(   25)       jpkm1    =>   jpkm1 ,  & !: jpk - 1  
(   26)       jpij     =>   jpij  ,  & !: jpi x jpj
(   27)       jpidta   =>   jpidta,  & !: first horizontal dimension  > or = jpi
(   28)       jpjdta   =>   jpjdta,  & !: second horizontal dimension > or = jpj
(   29)       jpkdta   =>   jpkdta,  & !: number of levels            > or = jpk
(   30)       lk_esopa =>   lk_esopa   !: flag to activate the all option
(   31) 
(   32) 
(   33)    !! namelist parameters      
(   34)    USE in_out_manager , ONLY :  &    
(   35)       cexper   =>   cexper,  & !: experience name for vairmer format 
(   36)       no       =>   no    ,  & !: job number
(   37)       nrstdt   =>   nrstdt,  & !: control of the time step (0,  & 1 or 2)
(   38)       nit000   =>   nit000,  & !: number of the first time step
(   39)       nitend   =>   nitend,  & !: number of the last time step
(   40)       nleapy   =>   nleapy,  & !: Leap year calendar (0/1)
(   41)       nwrite   =>   nwrite,  & !: frequency of OUTPUT file
(   42)       nstock   =>   nstock,  & !: frequency of restart file
(   43)       nprint   =>   nprint,  & !: level of print (0 no print)
(   44)       lwp      =>   lwp   ,  & !: boolean term for mpp output
(   45)       ln_ctl   =>   ln_ctl,  & !: = ln_ctl.AND.lwp (print control on the 1st proc)
(   46)       nictls   =>   nictls,  & !: Start i indice for the SUM control






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 2

(   47)       nictle   =>   nictle,  & !: End   i indice for the SUM control
(   48)       njctls   =>   njctls,  & !: Start j indice for the SUM control
(   49)       njctle   =>   njctle,  & !: End   j indice for the SUM control
(   50)       isplt    =>   isplt ,  & !: number of processors following i
(   51)       jsplt    =>   jsplt ,  & !: number of processors following j
(   52)       ijsplt   =>   ijsplt,  & !: nb of local domain = nb of processors
(   53)       numout   =>   numout     !: logical unit for output print
(   54) 
(   55) 
(   56)    !! run controm   
(   57)    USE in_out_manager , ONLY :  &  
(   58)       nstop     =>   nstop    ,  &  !: e r r o r  flag (=number of reason for a
(   59)       !                             !                   prematurely stop the run)
(   60)       nwarn     =>   nwarn    ,  &  !: w a r n i n g  flag (=number of warning
(   61)       !                             !                   found during the run) 
(   62)       cform_err =>   cform_err,  &  !:
(   63)       cform_war =>   cform_war      !:  
(   64)       
(   65)    USE dom_oce , ONLY :           &            
(   66)       lzoom      => lzoom     ,  & !: zoom flag
(   67)       lzoom_e    => lzoom_e   ,  & !: East  zoom type flag
(   68)       lzoom_w    => lzoom_w   ,  & !: West  zoom type flag
(   69)       lzoom_s    => lzoom_s   ,  & !: South zoom type flag
(   70)       lzoom_n    => lzoom_n   ,  & !: North zoom type flag
(   71)       lzoom_arct => lzoom_arct,  & !: ORCA    arctic zoom flag
(   72)       lzoom_anta => lzoom_anta     !: ORCA antarctic zoom flag
(   73) 
(   74) 
(   75) 
(   76)    USE dom_oce , ONLY :       & 
(   77)       nperio   =>   nperio,  & !: type of lateral boundary condition       
(   78) !      nlci     =>   nlci  ,  & !: index i for the sub domain left bottom 
(   79) !      nlcj     =>   nlcj  ,  & !: index j for the sub domain left bottom 
(   80)       nimpp    =>   nimpp ,  & !: i index for mpp-subdomain left bottom
(   81)       njmpp    =>   njmpp ,  & !: j index for mpp-subdomain left bottom
(   82)       nproc    =>   nproc ,  & !: number for local processor
(   83)       narea    =>   narea ,  & !: number for local area
(   84)       mig      =>   mig   ,  & !: local  ==> global  domain i-indice
(   85)       mjg      =>   mjg   ,  & !: local  ==> global  domain i-indice
(   86)       mi0      =>   mi0   ,  & !: global ==> local domain i-indice 
(   87)       mi1      =>   mi1   ,  & !: (mi0=1 and mi1=0 if the global indice is not in the local domain)
(   88)       mj0      =>   mj0   ,  & !: global ==> local domain j-indice 
(   89)       mj1      =>   mj1   ,  & !: (mj0=1 and mj1=0 if the global indice is not in the local domain)
(   90)       nidom    =>   nidom
(   91)  
(   92)    USE dom_oce , ONLY :       & 
(   93)       nimppt   => nimppt  ,  & !:i-indexes for each processor
(   94)       njmppt   => njmppt  ,  & !:j-indexes for each processor
(   95)       ibonit   => ibonit  ,  & !:i-processor neighbour existence
(   96)       ibonjt   => ibonjt  ,  & !:j- processor neighbour existence 
(   97)       nlci     => nlci    ,  & !:i- & j-dimensions of the local subdomain
(   98)       nlcj     => nlcj    ,  & !:
(   99)       nldi     => nldi    ,  & !:first and last indoor i- and j-indexes
(  100)       nlei     => nlei    ,  & !:
(  101)       nldj     => nldj    ,  & !:
(  102)       nlej     => nlej    ,  & !:
(  103)       nlcit    => nlcit   ,  & !:dimensions of every i-subdomain
(  104)       nlcjt    => nlcjt   ,  & !:dimensions of every j-subdomain






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 3

(  105)       nldit    => nldit   ,  & !:first indoor index for each i-domain 
(  106)       nleit    => nleit   ,  & !:last indoor index for each i-domain 
(  107)       nldjt    => nldjt   ,  & !:first indoor index for each j-domain 
(  108)       nlejt    => nlejt        !:last indoor index for each j-domain 
(  109) 
(  110)     
(  111)       !! horizontal curvilinear coordinate and scale factors
(  112)    USE dom_oce , ONLY :            &    
(  113)       glamt    =>   glamt ,  & !: longitude of t-point (degre)  
(  114)       glamu    =>   glamu ,  & !: longitude of t-point (degre)  
(  115)       glamv    =>   glamv ,  & !: longitude of t-point (degre)  
(  116)       glamf    =>   glamf ,  & !: longitude of t-point (degre)  
(  117)       gphit    =>   gphit ,  & !: latitude  of t-point (degre)   
(  118)       gphiu    =>   gphiu ,  & !: latitude  of t-point (degre)   
(  119)       gphiv    =>   gphiv ,  & !: latitude  of t-point (degre)   
(  120)       gphif    =>   gphif ,  & !: latitude  of t-point (degre)   
(  121)       e1t      =>   e1t   ,  & !: horizontal scale factors at t-point (m)  
(  122)       e2t      =>   e2t   ,  & !: horizontal scale factors at t-point (m)   
(  123)       e1u      =>   e1u   ,  & !: horizontal scale factors at u-point (m)
(  124)       e2u      =>   e2u   ,  & !: horizontal scale factors at u-point (m)
(  125)       e1v      =>   e1v   ,  & !: horizontal scale factors at v-point (m)
(  126)       e2v      =>   e2v        !: horizontal scale factors at v-point (m)  
(  127) 
(  128)    !! vertical coordinate and scale factors
(  129)    USE dom_oce , ONLY :              &   
(  130)       gdept    =>   gdept ,  & !: reference depth of t-points (m)
(  131)       e3t      =>   e3t   ,  & !: reference depth of t-points (m)  
(  132)       e3w      =>   e3w   ,  & !: reference depth of w-points (m)
(  133)       gdepw    =>   gdepw      !: reference depth of w-points (m)
(  134) 
(  135)    USE dom_oce ,   ONLY :            &      
(  136)       lk_zps   =>  lk_zps ,  & !: partial steps flag
(  137)       lk_sco   =>  lk_sco ,  & !: s-coordinate flag
(  138)       lk_zco   =>  lk_zco      !: z-coordinate flag
(  139) 
(  140)    USE lib_mpp ,   ONLY :            &     
(  141)       lk_mpp   =>  lk_mpp      !: Mpp flag
(  142) 
(  143)    USE dynspg_oce ,   ONLY :            &     
(  144)       lk_dynspg_rl   =>  lk_dynspg_rl      !: rigid lid flag
(  145) 
(  146) #if defined key_partial_steps
(  147)    !! Partial steps ('key_partial_steps')
(  148)    !! -----------------------------------
(  149)    USE dom_oce , ONLY :                & 
(  150)       gdep3w   =>  gdep3w  ,  & !: ???
(  151)       gdept_ps =>  gdept_ps,  & !: depth of t-points (m)
(  152)       gdepw_ps =>  gdepw_ps,  & !: depth of t-points (m)
(  153)       e3t_ps   =>  e3t_ps  ,  & !: vertical scale factors at t-
(  154)       e3u_ps   =>  e3u_ps  ,  & !: vertical scale factors at u-
(  155)       e3v_ps   =>  e3v_ps  ,  & !: vertical scale factors v-
(  156)       e3w_ps   =>  e3w_ps  ,  & !: w-points (m)
(  157)       e3f_ps   =>  e3f_ps  ,  & !: f-points (m)
(  158)       e3uw_ps  =>  e3uw_ps ,  & !: uw-points (m)
(  159)       e3vw_ps  =>  e3vw_ps      !: vw-points (m)
(  160) 
(  161)    USE oce , ONLY :                &
(  162)       gtu   =>  gtu  ,  & !: t- horizontal gradient at u-






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 4

(  163)       gtv   =>  gtv       !: and v-points at bottom ocean level
(  164) #endif
(  165) 
(  166) #if defined key_s_coord
(  167)    USE dom_oce , ONLY :              &   
(  168)       hbatt   =>   hbatt  ,  & !: ocean depth at the vertical of  t-point (m)
(  169)       hbatu   =>   hbatu  ,  & !: ocean depth at the vertical of  u-point (m)
(  170)       hbatv   =>   hbatv  ,  & !: ocean depth at the vertical of w-point (m)
(  171)       gsigt   =>   gsigt  ,  & !: model level depth coefficient at t-,  & w-levelsvertical scale factors at u-
(  172)       gsigw   =>   gsigw  ,  & !: model level depth coefficient at t-,  & w-levelsvertical scale factors v-
(  173)       gsi3w   =>   gsi3w  ,  & !: model level depth coef at w-levels (defined as the sum of e3w)
(  174)       esigt   =>   esigt  ,  & !: vertical scale factor coef. at t-levels
(  175)       esigw   =>   esigw       !: vertical scale factor coef. at w-levels
(  176) #endif
(  177) 
(  178)    !! masks, bathymetry
(  179)    USE dom_oce , ONLY :             &    
(  180)       mbathy   =>   mbathy,  & !: number of ocean level (=0,  & 1, ... , jpk-1) 
(  181)       tmask_i  =>   tmask_i, & !: Interior mask at t-points
(  182)       tmask    =>   tmask ,  & !: land/ocean mask at t-points
(  183)       umask    =>   umask ,  & !: land/ocean mask at u-points   
(  184)       vmask    =>   vmask ,  & !: land/ocean mask at v-points 
(  185)       fmask    =>   fmask      !: land/ocean mask at f-points 
(  186) 
(  187)    USE dom_oce , ONLY :         &
(  188)       n_cla   =>   n_cla       !: flag (0/1) for cross land advection 
(  189) 
(  190)    !! time domain
(  191)    USE dom_oce , ONLY :                 &
(  192)       neuler   =>   neuler,  & !: restart euler forward option (0=Euler)
(  193)       rdt      =>   rdt   ,  & !: time step for the dynamics 
(  194)       atfp     =>   atfp  ,  & !: asselin time filter parameter
(  195)       atfp1    =>   atfp1 ,  & !: asselin time filter coeff. (atfp1= 1-2*atfp)
(  196)       rdttra   =>   rdttra     !: vertical profile of tracer time step
(  197) 
(  198)    USE daymod , ONLY :                 &
(  199)       ndastp    =>   ndastp,  &    !: time step date in year/month/day aammjj
(  200)       nday_year =>   nday_year, &  !: curent day counted from jan 1st of the current year
(  201)       nyear     =>   nyear,   &  !: Current year
(  202)       nmonth    =>   nmonth,  &  !: Current month
(  203)       nday      =>   nday        !: Current day
(  204) 
(  205)    !! physical constants
(  206)    USE phycst ,   ONLY :                &  
(  207)       ra       =>   ra    ,  & !: earth radius
(  208)       rpi      =>   rpi   ,  & !: pi
(  209)       rday     =>   rday  ,  & !: day
(  210)       rauw     =>   rauw  ,  & !: density of pure water kg/m3
(  211)       ro0cpr   =>   ro0cpr,  & !: = 1. / ( rau0 * rcp )
(  212)       rad      =>   rad   ,  & !: conversion coeff. from degre into radian
(  213)       raass    =>   raass ,  & !: number of seconds in one year
(  214)       rmoss    =>   rmoss ,  & !: number of seconds in one month
(  215)       rjjss    =>   rjjss      !: number of seconds in one day
(  216) 
(  217)    !! present fields (now)
(  218)    USE oce , ONLY :            &     
(  219)       ua      =>    ua    ,  & !: i-horizontal velocity (m s-1) 
(  220)       va      =>    va    ,  & !: j-horizontal velocity (m s-1)






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 5

(  221)       un      =>    un    ,  & !: i-horizontal velocity (m s-1) 
(  222)       vn      =>    vn    ,  & !: j-horizontal velocity (m s-1)
(  223)       wn      =>    wn    ,  & !: vertical velocity (m s-1)  
(  224)       tn      =>    tn    ,  & !: pot. temperature (celsius)
(  225)       sn      =>    sn    ,  & !: salinity (psu)
(  226)       rhop    =>    rhop  ,  & !: potential volumic mass (kg m-3) 
(  227)       rhd     =>    rhd        !: in situ density anomalie rhd=(rho-rau0)/rau0 (no units)
(  228) 
(  229) #if defined key_trc_diatrd
(  230)    USE oce , ONLY :          &
(  231)       hdivn   =>    hdivn      !: horizontal divergence (1/s)
(  232) #endif
(  233) 
(  234) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  235)    !! wind speed
(  236)    USE blk_oce , ONLY :        &     
(  237)       vatm    =>    vatm       !: wind speed at sea surface (m s-1)
(  238) #endif
(  239) 
(  240)    !! wind speed
(  241)    USE taumod , ONLY :        &     
(  242)       taux    =>    taux ,  &  !: i-surface stress component
(  243)       tauy    =>    tauy       !: j-surface stress component
(  244) 
(  245) #if   defined key_trabbl_dif   ||   defined key_trabbl_adv
(  246)    USE trabbl , ONLY :           &      
(  247)       atrbbl   =>   atrbbl     !: lateral coeff. for bottom boundary layer scheme (m2/s)
(  248) #  if defined key_off_tra
(  249)    USE trabbl, ONLY :            &
(  250)       bblx   => bblx,       &
(  251)       bbly   => bbly
(  252) #  endif
(  253) #endif
(  254) 
(  255)    !! lateral diffusivity (tracers)
(  256)    USE ldftra_oce ,   ONLY :             &    
(  257)       aht0    =>   aht0  ,  &  !: horizontal eddy diffusivity for tracers (m2/s)
(  258)       ahtb0   =>   ahtb0 ,  &  !: background eddy diffusivity for isopycnal diff. (m2/s)
(  259)       ahtu    =>   ahtu  ,  &  !: lateral diffusivity coef. at u-points 
(  260)       ahtv    =>   ahtv  ,  &  !: lateral diffusivity coef. at v-points 
(  261)       ahtw    =>   ahtw  ,  &  !: lateral diffusivity coef. at w-points 
(  262)       ahtt    =>   ahtt  ,  &  !: lateral diffusivity coef. at t-points
(  263)       aeiv0   =>   aeiv0 ,  &  !: eddy induced velocity coefficient (m2/s) 
(  264)       aeiu    =>   aeiu  ,  &  !: eddy induced velocity coef. at u-points (m2/s)   
(  265)       aeiv    =>   aeiv  ,  &  !: eddy induced velocity coef. at v-points (m2/s) 
(  266)       aeiw    =>   aeiw        !: eddy induced velocity coef. at w-points (m2/s) 
(  267) 
(  268)    !! vertical diffusion
(  269)    USE zdf_oce , ONLY :      &    
(  270)       avt            =>   avt          ,  & !: vert. diffusivity coef. at w-point for temp  
(  271)       avt0           =>   avt0         ,  & !: vertical eddy diffusivity for tracers (m2/s)
(  272) #if ! defined key_off_tra
(  273)       l_trazdf_exp   =>   l_trazdf_exp ,  & !: explicit vertical diffusion scheme flag
(  274) #endif
(  275)       ln_zdfnpc      =>   ln_zdfnpc         !: convection: non-penetrative convection flag
(  276) 
(  277) 
(  278) #if defined key_zdfddm






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 6

(  279)    USE zdfddm , ONLY :             &     
(  280)       avs     =>    avs        !: salinity vertical diffusivity coeff. at w-point
(  281) #endif
(  282) 
(  283)    !! penetrative solar radiation
(  284)    USE traqsr , ONLY :            &      
(  285)       xsi1   =>   xsi1         !: first depth of extinction 
(  286) 
(  287)    !! surface fluxes
(  288)    USE ocesbc , ONLY :             &   
(  289)       qt      =>    qt    ,  & !: total surface heat flux (w m-2)   
(  290)       qsr     =>    qsr   ,  & !: penetrative solar radiation (w m-2)  
(  291)       emp     =>    emp   ,  & !: evaporation minus precipitation (kg m-2 s-2) 
(  292)       emps    =>    emps       !: evaporation minus precipitation (kg m-2 s-2)
(  293) 
(  294)    !! freezing area
(  295)    USE ocfzpt , ONLY :            &      
(  296)       freeze  =>    freeze,  & !: ice mask (0 or 1)  
(  297)       fzptn   =>    fzptn      !: now freezing temperature at ocean surface  
(  298) 
(  299) 
(  300)    !! mixing layer depth (turbocline)
(  301)    USE zdfmxl , ONLY :             &    
(  302)       hmld    =>   hmld   ,  & !: mixing layer depth (turbocline)
(  303)       hmlp    =>   hmlp   ,  & !: mixed layer depth  (rho=rho0+zdcrit) (m)
(  304)       hmlpt   =>   hmlpt       !: mixed layer depth at t-points (m)
(  305) 
(  306)    USE ldfslp , ONLY :              & 
(  307)       lk_ldfslp  =>  lk_ldfslp     !: slopes flag
(  308) #if   defined key_ldfslp
(  309)    !! direction of lateral diffusion (momentum  tracers) 
(  310)    USE ldfslp , ONLY :              & 
(  311)       uslp       =>   uslp    ,  & !: i-direction slope at u-, w-points
(  312)       vslp       =>   vslp    ,  & !: j-direction slope at v-, w-points
(  313)       wslpi      =>   wslpi   ,  & !: i-direction slope at u-, w-points
(  314)       wslpj      =>   wslpj        !: j-direction slope at v-, w-points
(  315) #endif
(  316) 
(  317)    !! ocean forcings runoff
(  318)    USE flxrnf , ONLY :              &   
(  319)       upsrnfh =>   upsrnfh ,  & !: mixed adv scheme in runoffs vicinity (hori.) 
(  320)       upsrnfz =>   upsrnfz ,  & !: mixed adv scheme in runoffs vicinity (vert.)
(  321)       upsadv  =>   upsadv       !: mixed adv scheme in straits vicinity (hori.)
(  322) 
(  323) END MODULE oce_trc



















PGF90 (Version     10.5)          09/23/2010  09:49:16      page 7

(    1) # 1 "oce_trc.F90"
(    1) MODULE oce_trc
(    2)    !!======================================================================
(    3)    !!                      ***  MODULE  oce_trc  ***
(    4)    !! Ocean passive tracer  :  share ocean-passive tracers variables
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   9.0  !  04-03  (C. Ethe)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!  TOP 1.0,  LOCEAN-IPSL (2005)
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/oce_trc.F90,v 1.13 2006/04/10 15:40:28 opalod Exp $
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    !! Domain characteristics
(   15)    USE par_oce , ONLY :       &
(   16)       cp_cfg   =>   cp_cfg,  & !: name of the configuration
(   17)       jp_cfg   =>   jp_cfg,  & !: resolution of the configuration
(   18)       jpiglo   =>   jpiglo,  & !: first  dimension of global domain --> i
(   19)       jpjglo   =>   jpjglo,  & !: second dimension of global domain --> j
(   20)       jpi      =>   jpi   ,  & !: first  dimension of grid --> i 
(   21)       jpj      =>   jpj   ,  & !: second dimension of grid --> j  
(   22)       jpk      =>   jpk   ,  & !: number of levels  
(   23)       jpim1    =>   jpim1 ,  & !: jpi - 1
(   24)       jpjm1    =>   jpjm1 ,  & !: jpj - 1 
(   25)       jpkm1    =>   jpkm1 ,  & !: jpk - 1  
(   26)       jpij     =>   jpij  ,  & !: jpi x jpj
(   27)       jpidta   =>   jpidta,  & !: first horizontal dimension  > or = jpi
(   28)       jpjdta   =>   jpjdta,  & !: second horizontal dimension > or = jpj
(   29)       jpkdta   =>   jpkdta,  & !: number of levels            > or = jpk
(   30)       lk_esopa =>   lk_esopa   !: flag to activate the all option
(   31) 
(   32) 
(   33)    !! namelist parameters      
(   34)    USE in_out_manager , ONLY :  &    
(   35)       cexper   =>   cexper,  & !: experience name for vairmer format 
(   36)       no       =>   no    ,  & !: job number
(   37)       nrstdt   =>   nrstdt,  & !: control of the time step (0,  & 1 or 2)
(   38)       nit000   =>   nit000,  & !: number of the first time step
(   39)       nitend   =>   nitend,  & !: number of the last time step
(   40)       nleapy   =>   nleapy,  & !: Leap year calendar (0/1)
(   41)       nwrite   =>   nwrite,  & !: frequency of OUTPUT file
(   42)       nstock   =>   nstock,  & !: frequency of restart file
(   43)       nprint   =>   nprint,  & !: level of print (0 no print)
(   44)       lwp      =>   lwp   ,  & !: boolean term for mpp output
(   45)       ln_ctl   =>   ln_ctl,  & !: = ln_ctl.AND.lwp (print control on the 1st proc)
(   46)       nictls   =>   nictls,  & !: Start i indice for the SUM control
(   47)       nictle   =>   nictle,  & !: End   i indice for the SUM control
(   48)       njctls   =>   njctls,  & !: Start j indice for the SUM control
(   49)       njctle   =>   njctle,  & !: End   j indice for the SUM control
(   50)       isplt    =>   isplt ,  & !: number of processors following i
(   51)       jsplt    =>   jsplt ,  & !: number of processors following j
(   52)       ijsplt   =>   ijsplt,  & !: nb of local domain = nb of processors
(   53)       numout   =>   numout     !: logical unit for output print
(   54) 
(   55) 
(   56)    !! run controm   
(   57)    USE in_out_manager , ONLY :  &  






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 8

(   58)       nstop     =>   nstop    ,  &  !: e r r o r  flag (=number of reason for a
(   59)       !                             !                   prematurely stop the run)
(   60)       nwarn     =>   nwarn    ,  &  !: w a r n i n g  flag (=number of warning
(   61)       !                             !                   found during the run) 
(   62)       cform_err =>   cform_err,  &  !:
(   63)       cform_war =>   cform_war      !:  
(   64)       
(   65)    USE dom_oce , ONLY :           &            
(   66)       lzoom      => lzoom     ,  & !: zoom flag
(   67)       lzoom_e    => lzoom_e   ,  & !: East  zoom type flag
(   68)       lzoom_w    => lzoom_w   ,  & !: West  zoom type flag
(   69)       lzoom_s    => lzoom_s   ,  & !: South zoom type flag
(   70)       lzoom_n    => lzoom_n   ,  & !: North zoom type flag
(   71)       lzoom_arct => lzoom_arct,  & !: ORCA    arctic zoom flag
(   72)       lzoom_anta => lzoom_anta     !: ORCA antarctic zoom flag
(   73) 
(   74) 
(   75) 
(   76)    USE dom_oce , ONLY :       & 
(   77)       nperio   =>   nperio,  & !: type of lateral boundary condition       
(   78) !      nlci     =>   nlci  ,  & !: index i for the sub domain left bottom 
(   79) !      nlcj     =>   nlcj  ,  & !: index j for the sub domain left bottom 
(   80)       nimpp    =>   nimpp ,  & !: i index for mpp-subdomain left bottom
(   81)       njmpp    =>   njmpp ,  & !: j index for mpp-subdomain left bottom
(   82)       nproc    =>   nproc ,  & !: number for local processor
(   83)       narea    =>   narea ,  & !: number for local area
(   84)       mig      =>   mig   ,  & !: local  ==> global  domain i-indice
(   85)       mjg      =>   mjg   ,  & !: local  ==> global  domain i-indice
(   86)       mi0      =>   mi0   ,  & !: global ==> local domain i-indice 
(   87)       mi1      =>   mi1   ,  & !: (mi0=1 and mi1=0 if the global indice is not in the local domain)
(   88)       mj0      =>   mj0   ,  & !: global ==> local domain j-indice 
(   89)       mj1      =>   mj1   ,  & !: (mj0=1 and mj1=0 if the global indice is not in the local domain)
(   90)       nidom    =>   nidom
(   91)  
(   92)    USE dom_oce , ONLY :       & 
(   93)       nimppt   => nimppt  ,  & !:i-indexes for each processor
(   94)       njmppt   => njmppt  ,  & !:j-indexes for each processor
(   95)       ibonit   => ibonit  ,  & !:i-processor neighbour existence
(   96)       ibonjt   => ibonjt  ,  & !:j- processor neighbour existence 
(   97)       nlci     => nlci    ,  & !:i- & j-dimensions of the local subdomain
(   98)       nlcj     => nlcj    ,  & !:
(   99)       nldi     => nldi    ,  & !:first and last indoor i- and j-indexes
(  100)       nlei     => nlei    ,  & !:
(  101)       nldj     => nldj    ,  & !:
(  102)       nlej     => nlej    ,  & !:
(  103)       nlcit    => nlcit   ,  & !:dimensions of every i-subdomain
(  104)       nlcjt    => nlcjt   ,  & !:dimensions of every j-subdomain
(  105)       nldit    => nldit   ,  & !:first indoor index for each i-domain 
(  106)       nleit    => nleit   ,  & !:last indoor index for each i-domain 
(  107)       nldjt    => nldjt   ,  & !:first indoor index for each j-domain 
(  108)       nlejt    => nlejt        !:last indoor index for each j-domain 
(  109) 
(  110)     
(  111)       !! horizontal curvilinear coordinate and scale factors
(  112)    USE dom_oce , ONLY :            &    
(  113)       glamt    =>   glamt ,  & !: longitude of t-point (degre)  
(  114)       glamu    =>   glamu ,  & !: longitude of t-point (degre)  
(  115)       glamv    =>   glamv ,  & !: longitude of t-point (degre)  






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 9

(  116)       glamf    =>   glamf ,  & !: longitude of t-point (degre)  
(  117)       gphit    =>   gphit ,  & !: latitude  of t-point (degre)   
(  118)       gphiu    =>   gphiu ,  & !: latitude  of t-point (degre)   
(  119)       gphiv    =>   gphiv ,  & !: latitude  of t-point (degre)   
(  120)       gphif    =>   gphif ,  & !: latitude  of t-point (degre)   
(  121)       e1t      =>   e1t   ,  & !: horizontal scale factors at t-point (m)  
(  122)       e2t      =>   e2t   ,  & !: horizontal scale factors at t-point (m)   
(  123)       e1u      =>   e1u   ,  & !: horizontal scale factors at u-point (m)
(  124)       e2u      =>   e2u   ,  & !: horizontal scale factors at u-point (m)
(  125)       e1v      =>   e1v   ,  & !: horizontal scale factors at v-point (m)
(  126)       e2v      =>   e2v        !: horizontal scale factors at v-point (m)  
(  127) 
(  128)    !! vertical coordinate and scale factors
(  129)    USE dom_oce , ONLY :              &   
(  130)       gdept    =>   gdept ,  & !: reference depth of t-points (m)
(  131)       e3t      =>   e3t   ,  & !: reference depth of t-points (m)  
(  132)       e3w      =>   e3w   ,  & !: reference depth of w-points (m)
(  133)       gdepw    =>   gdepw      !: reference depth of w-points (m)
(  134) 
(  135)    USE dom_oce ,   ONLY :            &      
(  136)       lk_zps   =>  lk_zps ,  & !: partial steps flag
(  137)       lk_sco   =>  lk_sco ,  & !: s-coordinate flag
(  138)       lk_zco   =>  lk_zco      !: z-coordinate flag
(  139) 
(  140)    USE lib_mpp ,   ONLY :            &     
(  141)       lk_mpp   =>  lk_mpp      !: Mpp flag
(  142) 
(  143)    USE dynspg_oce ,   ONLY :            &     
(  144)       lk_dynspg_rl   =>  lk_dynspg_rl      !: rigid lid flag
(  145) 
(  146) 
(  147) 
(  148) # 178
(  178)    !! masks, bathymetry
(  179)    USE dom_oce , ONLY :             &    
(  180)       mbathy   =>   mbathy,  & !: number of ocean level (=0,  & 1, ... , jpk-1) 
(  181)       tmask_i  =>   tmask_i, & !: Interior mask at t-points
(  182)       tmask    =>   tmask ,  & !: land/ocean mask at t-points
(  183)       umask    =>   umask ,  & !: land/ocean mask at u-points   
(  184)       vmask    =>   vmask ,  & !: land/ocean mask at v-points 
(  185)       fmask    =>   fmask      !: land/ocean mask at f-points 
(  186) 
(  187)    USE dom_oce , ONLY :         &
(  188)       n_cla   =>   n_cla       !: flag (0/1) for cross land advection 
(  189) 
(  190)    !! time domain
(  191)    USE dom_oce , ONLY :                 &
(  192)       neuler   =>   neuler,  & !: restart euler forward option (0=Euler)
(  193)       rdt      =>   rdt   ,  & !: time step for the dynamics 
(  194)       atfp     =>   atfp  ,  & !: asselin time filter parameter
(  195)       atfp1    =>   atfp1 ,  & !: asselin time filter coeff. (atfp1= 1-2*atfp)
(  196)       rdttra   =>   rdttra     !: vertical profile of tracer time step
(  197) 
(  198)    USE daymod , ONLY :                 &
(  199)       ndastp    =>   ndastp,  &    !: time step date in year/month/day aammjj
(  200)       nday_year =>   nday_year, &  !: curent day counted from jan 1st of the current year
(  201)       nyear     =>   nyear,   &  !: Current year
(  202)       nmonth    =>   nmonth,  &  !: Current month






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 10

(  203)       nday      =>   nday        !: Current day
(  204) 
(  205)    !! physical constants
(  206)    USE phycst ,   ONLY :                &  
(  207)       ra       =>   ra    ,  & !: earth radius
(  208)       rpi      =>   rpi   ,  & !: pi
(  209)       rday     =>   rday  ,  & !: day
(  210)       rauw     =>   rauw  ,  & !: density of pure water kg/m3
(  211)       ro0cpr   =>   ro0cpr,  & !: = 1. / ( rau0 * rcp )
(  212)       rad      =>   rad   ,  & !: conversion coeff. from degre into radian
(  213)       raass    =>   raass ,  & !: number of seconds in one year
(  214)       rmoss    =>   rmoss ,  & !: number of seconds in one month
(  215)       rjjss    =>   rjjss      !: number of seconds in one day
(  216) 
(  217)    !! present fields (now)
(  218)    USE oce , ONLY :            &     
(  219)       ua      =>    ua    ,  & !: i-horizontal velocity (m s-1) 
(  220)       va      =>    va    ,  & !: j-horizontal velocity (m s-1)
(  221)       un      =>    un    ,  & !: i-horizontal velocity (m s-1) 
(  222)       vn      =>    vn    ,  & !: j-horizontal velocity (m s-1)
(  223)       wn      =>    wn    ,  & !: vertical velocity (m s-1)  
(  224)       tn      =>    tn    ,  & !: pot. temperature (celsius)
(  225)       sn      =>    sn    ,  & !: salinity (psu)
(  226)       rhop    =>    rhop  ,  & !: potential volumic mass (kg m-3) 
(  227)       rhd     =>    rhd        !: in situ density anomalie rhd=(rho-rau0)/rau0 (no units)
(  228) 
(  229) 
(  230) # 235
(  235)    !! wind speed
(  236)    USE blk_oce , ONLY :        &     
(  237)       vatm    =>    vatm       !: wind speed at sea surface (m s-1)
(  238) 
(  239) # 240
(  240)    !! wind speed
(  241)    USE taumod , ONLY :        &     
(  242)       taux    =>    taux ,  &  !: i-surface stress component
(  243)       tauy    =>    tauy       !: j-surface stress component
(  244) 
(  245) 
(  246) # 255
(  255)    !! lateral diffusivity (tracers)
(  256)    USE ldftra_oce ,   ONLY :             &    
(  257)       aht0    =>   aht0  ,  &  !: horizontal eddy diffusivity for tracers (m2/s)
(  258)       ahtb0   =>   ahtb0 ,  &  !: background eddy diffusivity for isopycnal diff. (m2/s)
(  259)       ahtu    =>   ahtu  ,  &  !: lateral diffusivity coef. at u-points 
(  260)       ahtv    =>   ahtv  ,  &  !: lateral diffusivity coef. at v-points 
(  261)       ahtw    =>   ahtw  ,  &  !: lateral diffusivity coef. at w-points 
(  262)       ahtt    =>   ahtt  ,  &  !: lateral diffusivity coef. at t-points
(  263)       aeiv0   =>   aeiv0 ,  &  !: eddy induced velocity coefficient (m2/s) 
(  264)       aeiu    =>   aeiu  ,  &  !: eddy induced velocity coef. at u-points (m2/s)   
(  265)       aeiv    =>   aeiv  ,  &  !: eddy induced velocity coef. at v-points (m2/s) 
(  266)       aeiw    =>   aeiw        !: eddy induced velocity coef. at w-points (m2/s) 
(  267) 
(  268)    !! vertical diffusion
(  269)    USE zdf_oce , ONLY :      &    
(  270)       avt            =>   avt          ,  & !: vert. diffusivity coef. at w-point for temp  
(  271)       avt0           =>   avt0         ,  & !: vertical eddy diffusivity for tracers (m2/s)
(  273)       l_trazdf_exp   =>   l_trazdf_exp ,  & !: explicit vertical diffusion scheme flag






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 11

(  275)       ln_zdfnpc      =>   ln_zdfnpc         !: convection: non-penetrative convection flag
(  276) 
(  277) 
(  278) 
(  279) # 283
(  283)    !! penetrative solar radiation
(  284)    USE traqsr , ONLY :            &      
(  285)       xsi1   =>   xsi1         !: first depth of extinction 
(  286) 
(  287)    !! surface fluxes
(  288)    USE ocesbc , ONLY :             &   
(  289)       qt      =>    qt    ,  & !: total surface heat flux (w m-2)   
(  290)       qsr     =>    qsr   ,  & !: penetrative solar radiation (w m-2)  
(  291)       emp     =>    emp   ,  & !: evaporation minus precipitation (kg m-2 s-2) 
(  292)       emps    =>    emps       !: evaporation minus precipitation (kg m-2 s-2)
(  293) 
(  294)    !! freezing area
(  295)    USE ocfzpt , ONLY :            &      
(  296)       freeze  =>    freeze,  & !: ice mask (0 or 1)  
(  297)       fzptn   =>    fzptn      !: now freezing temperature at ocean surface  
(  298) 
(  299) 
(  300)    !! mixing layer depth (turbocline)
(  301)    USE zdfmxl , ONLY :             &    
(  302)       hmld    =>   hmld   ,  & !: mixing layer depth (turbocline)
(  303)       hmlp    =>   hmlp   ,  & !: mixed layer depth  (rho=rho0+zdcrit) (m)
(  304)       hmlpt   =>   hmlpt       !: mixed layer depth at t-points (m)
(  305) 
(  306)    USE ldfslp , ONLY :              & 
(  307)       lk_ldfslp  =>  lk_ldfslp     !: slopes flag
(  308) # 309
(  309)    !! direction of lateral diffusion (momentum  tracers) 
(  310)    USE ldfslp , ONLY :              & 
(  311)       uslp       =>   uslp    ,  & !: i-direction slope at u-, w-points
(  312)       vslp       =>   vslp    ,  & !: j-direction slope at v-, w-points
(  313)       wslpi      =>   wslpi   ,  & !: i-direction slope at u-, w-points
(  314)       wslpj      =>   wslpj        !: j-direction slope at v-, w-points
(  315) 
(  316) # 317
(  317)    !! ocean forcings runoff
(  318)    USE flxrnf , ONLY :              &   
(  319)       upsrnfh =>   upsrnfh ,  & !: mixed adv scheme in runoffs vicinity (hori.) 
(  320)       upsrnfz =>   upsrnfz ,  & !: mixed adv scheme in runoffs vicinity (vert.)
(  321)       upsadv  =>   upsadv       !: mixed adv scheme in straits vicinity (hori.)
(  322) 
(  323) END MODULE oce_trc
