MODULE in_out_manager   

   USE lib_print         ! formated print library
   USE par_kind
   USE par_oce

   PUBLIC
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/in_out_manager.F90,v 1.1.1.1 2005/11/14 10:41:07 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! namelist parameters
   !! -------------------------------------
   ! namrun:  parameters of the run
   CHARACTER (len=16) ::    &   !:
      cexper = "exp0"           !: experiment name used for output filename
   
   LOGICAL ::   &              !!: * namelist namrun *
      ln_rstart = .FALSE. ,  &  !: start from (F) rest or (T) a restart file
      ln_ctl    = .FALSE.       !: run control for debugging
   
   INTEGER ::                & !!: * namelist namrun *
      no     = 0        ,    &  !: job number
      nrstdt = 0        ,    &  !: control of the time step (0, 1 or 2)
      nit000 = 1        ,    &  !: index of the first time step
      nitend = 10       ,    &  !: index of the last time step
      ndate0 = 961115   ,    &  !: initial calendar date aammjj
      nleapy = 0        ,    &  !: Leap year calendar flag (0/1 or 30)
      ninist = 0        ,    &  !: initial state output flag (0/1)
      nbench = 0                !: benchmark parameter (0/1)
   !!----------------------------------------------------------------------
   !!                          Run control  
   !!----------------------------------------------------------------------
   
   INTEGER ::                &  !:
      nstop = 0 ,            &  !: e r r o r  flag (=number of reason for a
      !                         !                   prematurely stop the run)
      nwarn = 0                 !: w a r n i n g  flag (=number of warning
      !                         !                       found during the run)

   
   CHARACTER (len=64) ::        &                                                    !:
      cform_err="(/,' ===>>> : E R R O R',     /,'         ===========',/)"    ,   & !:
      cform_war="(/,' ===>>> : W A R N I N G', /,'         ===============',/)"      !:
   !!----------------------------------------------------------------------
   !! output monitoring
   !! -----------------------------------

   LOGICAL ::   &               !:
      lwp                ,   &  !: boolean : true on the 1st processor only
      lsp_area = .TRUE.         !: to make a control print over a specific area

   INTEGER ::                &  !:
      nstock =   10 ,        &  !: restart file frequency
      nprint =    0 ,        &  !: level of print (0 no print)
      nwrite =   10 ,        &  !: restart file frequency
      nictls =    0 ,        &  !: Start i indice for the SUM control
      nictle =    0 ,        &  !: End   i indice for the SUM control
      njctls =    0 ,        &  !: Start j indice for the SUM control
      njctle =    0 ,        &  !: End   j indice for the SUM control
      isplt  =    1 ,        &  !: number of processors following i
      jsplt  =    1 ,        &  !: number of processors following j
      ijsplt =    1             !: nb of local domain = nb of processors

   !!----------------------------------------------------------------------
   !! logical units
   !! ------------------------------
   INTEGER ::                &  !:
      numstp     =  1 ,      &  !: logical unit for time step
      numout     =  2 ,      &  !: logical unit for output print
      numnam     =  3 ,      &  !: logical unit for namelist
      numnam_ice =  4 ,      &  !: logical unit for ice namelist
      numevo_ice = 17 ,      &  !: logical unit for ice variables (temp. evolution)
      numsol     = 25 ,      &  !: logical unit for solver statistics
      numwri     = 40 ,      &  !: logical unit for output write
      numisp     = 41 ,      &  !: logical unit for island statistics
      numgap     = 45 ,      &  !: logical unit for differences diagnostic
      numwrs     = 46 ,      &  !: logical unit for output restart
      numtdt     = 62 ,      &  !: logical unit for data temperature
      numsdt     = 63 ,      &  !: logical unit for data salinity
      numrnf     = 64 ,      &  !: logical unit for runoff data
      numwso     = 71 ,      &  !: logical unit for 2d output write
      numwvo     = 72 ,      &  !: logical unit for 3d output write
      numsst     = 65 ,      &  !: logical unit for surface temperature data
      numbol     = 67 ,      &  !: logical unit for "bol" diagnostics
      numptr     = 68 ,      &  !: logical unit for Poleward TRansports
      numflo     = 69 ,      &  !: logical unit for drifting floats
      !                         !: * coupled units
      numlhf     = 71 ,      &  !: unit to transfer fluxes
      numlws     = 72 ,      &  !: unit to transfer stress
      numlts     = 73 ,      &  !: unit to transfer sst
      numlic     = 74           !: unit to transfer ice cover


   !! Contral/debugging
   !! -----------------
   REAL(wp) ::               &  !:
      u_ctl, v_ctl,          &  !: sum of ua and va trend
      t_ctl, s_ctl              !: sum of ta and sa trend

END MODULE in_out_manager
