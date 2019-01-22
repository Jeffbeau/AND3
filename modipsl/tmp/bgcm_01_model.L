


PGF90 (Version     10.5)          09/23/2010  09:49:24      page 1

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

Filename: bgcm_01_model.F90

(    1) !!DB: 2009.08.17
(    2) !!Eliminated some routines that are no longer necessary. If needed look in
(    3) !!bgcm_01_model.F91 (if it exists locally) or some achived version of the code
(    4) !!NB: RE Open Boundary Conditions:
(    5) !!NB: trc_trp() calls trc_nxt() which calls update_boundary_vals()
(    6) !!which is the correct location for updating the tracer arrays with OB values
(    7) 
(    8) 
(    9) !!DB: 2009.06.09
(   10) !!Latest version which uses "copy-to-boundary" OBC scheme for PZ and 
(   11) !!N_obc routine for N
(   12) !!Note that without further modifying the code structure (which is not
(   13) !!necessarily a bad idea) --
(   14) !!  All calls related to assignment of OBCs must be done from trcnxt.F90
(   15) !!Thus look there for OBC scheme currently used
(   16) 
(   17) 
(   18) !!OLD comments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   19) !!ca. 2009.04
(   20) !!It also contains OBC-related fixes relative to previous versions
(   21) !!(see also TRP/trcnxt.F90 for OBC-related change)
(   22) !!DB 2008.08.26 ...
(   23) !!Replace trcstp with 1 simpler module specific to key_BGCM_01
(   24) !!See TRP/trcstp.F90 to see how this module is used (i.e. how some
(   25) !!backward compatibility is maintained)
(   26) 
(   27) !!The BioGeoChemical Model works by (1) calling a Source-Minus-Sink (SMS)
(   28) !!routine that does the vertical biophysics only; and (2) calling the
(   29) !!regular tracer advection/diffusion schemes. The tracer array(s) default 
(   30) !!to trn(:,:,:,jptra) (see trc.F90 and par_trc ----> par_trc_trp). 
(   31) !!jptra is the number of 3D BGCM fields used. 
(   32) 
(   33) 
(   34) !!Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   35) !!-1- The main BGCM is in bgcm_sms() -- at bottom of this module. It is currently
(   36) !!    v.simple. This is where the major modifications need to be done.
(   37) !!-2- The SMS routine would get its various model and runtime parameters in the module
(   38) !!    bgcm_01_initrc.F90. At this time there is only a place-holder for this
(   39) !!    code in that module. 
(   40) !!END Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   41) 
(   42) !!END OLD comments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   43) 
(   44) #if defined key_passivetrc
(   45) MODULE bgcm_01_model
(   46)    !!======================================================================






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 2

(   47)    !!                       ***  MODULE trcstp  ***
(   48)    !! Time-stepping    : time loop of opa for passive tracer
(   49)    !!======================================================================
(   50) 
(   51) 
(   52)    !!----------------------------------------------------------------------
(   53)    !!   trc_stp      : passive tracer system time-stepping
(   54)    !!----------------------------------------------------------------------
(   55)    !! * Modules used
(   56)    USE oce_trc          ! ocean dynamics and active tracers variables
(   57)    USE trc              ! ocean passive tracers variables 
(   58)    USE trctrp           ! passive tracers transport
(   59) 
(   60)    USE lib_bgcm_01
(   61)    USE lib_ncdf
(   62) !!DB
(   63)    USE oce, ONLY : sw_rad
(   64) 
(   65)    IMPLICIT NONE
(   66)    PRIVATE
(   67) 
(   68)    !! * Routine accessibility
(   69)    PUBLIC trc_stp           ! called by step
(   70) 
(   71) !!DB
(   72)    PUBLIC bgcm_sms
(   73) 
(   74) !!DB: bgcm_01.inc file variables used by step_npz_1D
(   75)   INTEGER :: last_lev
(   76)   REAL, DIMENSION(jpk) :: NN,PP,ZZ, NN0,PP0,ZZ0
(   77)   REAL(wp), DIMENSION(jpk+1) ::  wmask_1D
(   78) 
(   79) 
(   80) CONTAINS
(   81) 
(   82)    SUBROUTINE trc_stp( kt, kindic )
(   83)       !!-------------------------------------------------------------------
(   84)       !!                     ***  ROUTINE trc_stp  ***
(   85)       !!                      
(   86)       !! ** Purpose : Time loop of opa for passive tracer
(   87)       !! 
(   88)       !! ** Method  : 
(   89)       !!              Compute the passive tracers trends 
(   90)       !!              Update the passive tracers
(   91)       !!
(   92)       !! History :
(   93)       !!   9.0  !  04-03  (C. Ethe)  Original
(   94)       !!-------------------------------------------------------------------
(   95) 
(   96)      USE lbclnk
(   97) 
(   98)       !! * Arguments
(   99)       INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
(  100)       INTEGER, INTENT( in ) ::  kindic
(  101)       CHARACTER (len=25) :: charout
(  102) 
(  103) !!DB
(  104)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zfield  






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 3

(  105)       integer :: ji,jj,jk,jn, rec_num, status
(  106)       CHARACTER (len=80) :: fname
(  107)       REAL(wp), DIMENSION(100) :: bgcm_params
(  108)       REAL(wp) :: tmp1
(  109) 
(  110) 
(  111) !!DB -- 2009.10.07 -- write params to file (HARDWIRED):
(  112)       if(lwp .AND. kt==nit000) then
(  113) !!output to numout2
(  114)          write(numout2,*)'----------------------'
(  115)          write(numout2,*)'BGCM Parameters '
(  116)          write(numout2,'(a20,2x,e12.6)')'Vm = ',Vm
(  117)          write(numout2,'(a20,2x,e12.6)')'Rm = ',Rm
(  118)          write(numout2,'(a20,2x,e12.6)')'ks = ',ks
(  119)          write(numout2,'(a20,2x,e12.6)')'lambda = ',lambda
(  120)          write(numout2,'(a20,2x,e12.6)')'gamma = ',gamma
(  121)          write(numout2,'(a20,2x,e12.6)')'M_ZOO = ', M_ZOO
(  122)          write(numout2,'(a20,2x,e12.6)')'M_PHY = ', M_PHY
(  123)          write(numout2,'(a20,2x,e12.6)')'kext = ', kext
(  124)          write(numout2,'(a20,2x,e12.6)')'aaa = ', aaa
(  125)          write(numout2,'(a20,2x,e12.6)')'PAR = ', PAR
(  126)          write(numout2,'(a20,2x,e12.6)')'w_sink = ', w_sink
(  127)          write(numout2,*)'----------------------'
(  128)       endif
(  129) 
(  130) 
(  131) !!DB 2009.10.23 -- need sw_rad variable qsr() which is computed elsewhere
(  132) !!(method depends on the KEY). A problem arises in that tra_sbc() zeros it depending
(  133) !!on penetrative radiation flag = ln_trasqr. To guarantee that this routine has this
(  134) !!variable I search for a non-zero value and only assign sw_rad = qsr if qsr /= 0.
(  135) !!Note that I do this at every kt as qsr is not necessarily computed every kt so by
(  136) !!bad luck if I did this after the "return" statement below I may never get non-zero values.
(  137) 
(  138) !!NB: qsr() may not exist at all timesteps so only assign sw_rad when it does
(  139) !!NB: I use HARDWIRED position (10,10) which should always be OK 
(  140)     tmp1 = 0.0
(  141)     tmp1 = qsr(10,10)
(  142)     call mpp_sum(tmp1)
(  143)     if(tmp1 /= 0.0) then
(  144)        sw_rad(:,:) = qsr(:,:)
(  145) !Limit notification to 1/day if this works; if not it's no big deal
(  146)        if(lwp .AND. mod(kt-nit000,int(rday/rdt)) ==0) then
(  147)           write(numout2,*)'BGCM Updating sw_rad  at kt = ', kt
(  148)        endif
(  149)     endif
(  150) 
(  151)       ! this ROUTINE is called only every ndttrc time step
(  152)       IF( MOD( kt , ndttrc ) /= 0 ) RETURN
(  153) 
(  154) !!This is the routine that calls the BGCM Source-Minus-Sink code
(  155)       CALL bgcm_sms( kt )
(  156) 
(  157)       ! transport of passive tracers
(  158) !!DB: No need to change this; see trctrp.F90
(  159) !!NB: trc_trp() calls trc_nxt() which calls update_boundary_vals()
(  160) !!which is the correct location for updating the tracer arrays with OB values
(  161)       CALL trc_trp( kt )
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 4

(  163) !!DB: Open Boundary code 
(  164) !!This is where special OBC-related operations could also be done
(  165) #ifdef key_obc
(  166) 
(  167) #endif
(  168) 
(  169) 
(  170) !!DB: write restart -- not yet written
(  171) !!Restart file will be v.similar to regular output file 
(  172) !!the last frame of which can be used for the time-being. 
(  173) !      CALL trc_wri( kt )            ! outputs
(  174) 
(  175) !!DB: output model variables
(  176)       if(mod( kt, nwritetrc ) == 0 ) then
(  177) !!DB: write to ncdf file
(  178)          rec_num = (kt-nit000+1)/nwritetrc
(  179)          CALL ncdf_write(bgcm_fname, 'time_counter', REAL(kt * rdt), rec_num, status)
(  180) 
(  181) !!DB: 2008.10.17: Both of the below work 
(  182) !!DB: I use the NPZ version but keep the trn setup in the netcdf file in case
(  183) !!DB: it is desired to use it. (See lib_bgcm_01.F90 variable jf) 
(  184) !!DB: Note that to force a write of a 2D,3D,4D variable directly to a record number
(  185) !!    use a -ve value for rec_num
(  186) 
(  187) !         CALL ncdf_write(bgcm_fname, 'trn', trn, jptra, -rec_num, status)
(  188) !!DB: or
(  189)          zfield(:,:,:) = trn(:,:,:,1)
(  190)          CALL ncdf_write(bgcm_fname, 'N', zfield, -rec_num, status)
(  191)          zfield(:,:,:) = trn(:,:,:,2)
(  192)          CALL ncdf_write(bgcm_fname, 'P', zfield, -rec_num, status)
(  193)          zfield(:,:,:) = trn(:,:,:,3)
(  194)          CALL ncdf_write(bgcm_fname, 'Z', zfield, -rec_num, status)
(  195) !!AD/DB 2009.09.30
(  196)          CALL ncdf_write(bgcm_fname, 'ndastp',REAL(ndastp), rec_num, status)
(  197)          CALL ncdf_write(bgcm_fname, 'model_time_step',REAL(kt), rec_num, status)
(  198)          CALL ncdf_write(bgcm_fname, 'model_time',model_time, rec_num, status)
(  199) 
(  200)       endif
(  201) 
(  202) !!DB: tra is done in trc_trp ---> trc_nxt. 
(  203) !!To be safe do all of them (as done in p4zprg)
(  204) !!NB: Likely do not have to do this
(  205)       DO jn=1 , jptra
(  206)         CALL lbc_lnk(trn(:,:,:,jn), 'T', 1. )
(  207)         CALL lbc_lnk(trb(:,:,:,jn), 'T', 1. )
(  208)         CALL lbc_lnk(tra(:,:,:,jn), 'T', 1. )
(  209)       END DO
(  210) 
(  211) 
(  212) 
(  213)    END SUBROUTINE trc_stp
(  214) 
(  215) !!DB: BGCM_01 SMS model
(  216)    SUBROUTINE bgcm_sms( kt )
(  217)       !!===========================================================================================
(  218) 
(  219)       !! * Arguments
(  220)       !! -----------






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 5

(  221)      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
(  222)      
(  223) !!DB
(  224)      real(wp) :: dk_tau
(  225) !!DB -- example of using pointers to make NPZ code more explicit
(  226)      real(wp), pointer :: N(:,:,:), P(:,:,:), Z(:,:,:)
(  227)      real(wp), pointer :: N0(:,:,:), P0(:,:,:), Z0(:,:,:)
(  228)      
(  229) 
(  230) !!DB -- assign pointers to default trn, trb arrays
(  231)      N => trn(:,:,:,1)
(  232)      P => trn(:,:,:,2)
(  233)      Z => trn(:,:,:,3)
(  234)      N0 => trb(:,:,:,1)
(  235)      P0 => trb(:,:,:,2)
(  236)      Z0 => trb(:,:,:,3)
(  237)      
(  238)      !! this ROUTINE is called only every ndttrc time step
(  239)      !! --------------------------------------------------
(  240)      IF ( MOD(kt,ndttrc) /= 0) RETURN
(  241)      
(  242)      call step_npz(kt)
(  243)      
(  244)    END SUBROUTINE bgcm_sms
(  245) 
(  246) 
(  247) 
(  248) !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
(  249) ! Simple P Franks routine, working on *2 versions of variables
(  250) ! Call this after calls to STEP_TRACER which handles the adv/diff part
(  251) ! of field updates. 
(  252) !
(  253) ! REM: I0(d) = 125 + 100*cos(2*pi/360 * (day-180))
(  254) !            = fit by DB to radiation at ~ 45N
(  255) !
(  256) !REM: These model parameters are set at top of module
(  257) !    Vm = 2.0/86400.0  !! convert per day to per second
(  258) !    Rm = 0.5/86400.0  
(  259) !    ks = 1.0          !! \mu Mol N per liter
(  260) !    lambda = 0.2      !! grazing
(  261) !    gamma = 0.7
(  262) !    w_sink = -1.5/86400.0 !!sinking rate in m/day
(  263) !    M_ZOO = 0.04/86400.0 !! Zoo mort
(  264) !    M_PHY = 0.05/86400.0   !! M_PHY
(  265) !    lambda = 0.3      !! grazing
(  266) !    kext = 0.1        !!in m
(  267) !    aaa = 0.025       !! initial slope of P-I curve
(  268) !    PAR = 0.43        !!from Oschlies
(  269) 
(  270)   SUBROUTINE step_npz(kt)
(  271) 
(  272) #if ! defined key_flx_bulk_daily
(  273)     USE flxblk_solar
(  274) #endif
(  275) 
(  276) !!DB
(  277)     USE oce, ONLY : sw_rad
(  278) !!DBG






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 6

(  279) !    USE lib_mpp
(  280) 
(  281)     INTEGER, INTENT( in ) ::   kt
(  282) 
(  283)     REAL(wp) :: DTIN, tmp2, tmp1, fac, fac0
(  284)     REAL(wp), DIMENSION(jpi,jpj) :: C2
(  285)     REAL(wp), DIMENSION(jpi,jpj,2) ::  PL
(  286)     integer  :: k,l,j,lt,lb,L0,i
(  287) 
(  288) !!DB -- example of using pointers to make NPZ code more explicit
(  289)     real(wp), pointer :: N(:,:,:), P(:,:,:), Z(:,:,:)
(  290)     real(wp), pointer :: N0(:,:,:), P0(:,:,:), Z0(:,:,:)
(  291) 
(  292) !!DB -- assign pointers to default trn, trb arrays
(  293)     N => trn(:,:,:,1)
(  294)     P => trn(:,:,:,2)
(  295)     Z => trn(:,:,:,3)
(  296)     N0 => trb(:,:,:,1)
(  297)     P0 => trb(:,:,:,2)
(  298)     Z0 => trb(:,:,:,3)
(  299)     
(  300)     !===================================
(  301)     ! NPZ interaction equations
(  302)     !===================================
(  303) 
(  304) !!DB: Kept for reference purposes; should never be used anymore
(  305) !!If not calculating solar forcing elsewhere then must use the next lines
(  306) !    I_sfce = 125.0 + 100.0*cos(2.*rpi/360.0 * ((nday_year)-180.0)) !! W/m^2
(  307) !    I_sfce = PAR * aaa * I_sfce
(  308) !    C2(:,:) = (Vm*86400.0/I_sfce)**2    
(  309) !!Also must comment-out the I_sfce line below
(  310) 
(  311) #if ! defined key_flx_bulk_daily 
(  312) !!DB -- Use AD solar flux ---> ad_qsr_oce(i,j) 
(  313)    call flx_blk_solar(kt)
(  314)     sw_rad(:,:) = ad_qsr_oce(:,:)
(  315) #endif
(  316) 
(  317) !!     Do sfce layer first to simplify calc of sinking term: w_sink * dP/dz
(  318) !!NB C2(i,j) is computed here so it is not necessary to do so elsewhere
(  319)     lt = 1
(  320)     lb = 2
(  321)     do k = 1, 1
(  322)        L0 = k + 1
(  323)        do j = 1, jpj        
(  324)           do i = 1, jpi
(  325) !!DB 2009.10.23 -- below is the desired line-of-code
(  326) !!  However, sw_rad can contain ~zeros and depending on the machine this can result
(  327) !!  in a massive slowdown of execution. So I add 0.1 to sw_rad to avoid this problem.
(  328) !!  As of the above date, this seems to work
(  329) !             I_sfce = PAR * aaa * sw_rad(i,j)  !current desired line-of-code
(  330)              I_sfce = PAR * aaa * (sw_rad(i,j)+0.1)  !DBG
(  331)              C2(i,j) = (Vm*86400.0/I_sfce)**2    
(  332)              
(  333)              DTIN =  rdt * tmask(i,j,k)
(  334)              !                
(  335)              PL(i,j,lt) = 0.0
(  336) !!     Upwind diff






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 7

(  337)                 fac0 = sign(1.0,w_sink)
(  338)                 PL(i,j,lb) = ((1.+fac0)*P0(i,j,L0) +  &
(  339)                      (1.-fac0)*P0(i,j,k))/2.
(  340) !                
(  341) !!      this is z-integral of actual f(z,t); 
(  342)                 Ibar = 0.0 +            &
(  343)                      (1.0/e3t(k)) * 1.0/kext *                           &
(  344)                                       (                                  &
(  345)                                       log( (exp(-kext*gdept(k))+         &
(  346)                                       sqrt(C2(i,j)+exp(-2.*kext*gdept(k)))) / &
(  347)                                       (exp(-kext*gdept(k+1))+            &
(  348)                                       sqrt(C2(i,j)+exp(-2.*kext*gdept(k+1))))  ) &
(  349)                                       ) 
(  350) !                
(  351)                 P(i,j,k) = P(i,j,k) + DTIN * (          &
(  352)                      Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)    & 
(  353)                      - Z0(i,j,k)*Rm*(1.0-exp(-lambda*P0(i,j,k)))   & 
(  354)                      - M_PHY*P0(i,j,k)                  & 
(  355) !!DB: use  wmask here (computed in lib_bgcm_01)
(  356)                   -(wmask(i,j,k)*w_sink*PL(i,j,lt)-     &
(  357)                   wmask(i,j,L0)*w_sink*PL(i,j,lb))/e3t(k) )  !!sinking term
(  358) !                
(  359)                 Z(i,j,k) = Z(i,j,k) + DTIN * (           &
(  360)                      gamma*Rm*Z0(i,j,k)*(1.0-exp(-lambda*P0(i,j,k)))        &
(  361)                      - M_ZOO*Z0(i,j,k) )
(  362) !                
(  363)                 N(i,j,k) = N(i,j,k) + DTIN * (            &
(  364)                      -Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)        &
(  365)                      +(1.0-gamma)*Rm*Z0(i,j,k)*           &
(  366)                      (1.0-exp(-lambda*P0(i,j,k)))         &
(  367)                      + M_PHY*P0(i,j,k) + M_ZOO*Z0(i,j,k) )
(  368) !                
(  369)           enddo
(  370)        enddo
(  371)     enddo
(  372)     
(  373) !    
(  374) !!DB: as w_sink is always of the same sign (-ve) compute fac0 outside of loop
(  375) !!Use wmask to stop flux through bottom
(  376)     fac0 = sign(1.0,w_sink)
(  377)     lt = 2
(  378)     lb = 1
(  379)     do k = 2, jpkm1
(  380)        L0 = k + 1
(  381)        do j = 1, jpj
(  382)           do i = 1, jpi
(  383) 
(  384)              DTIN =  rdt * tmask(i,j,k)
(  385) !!     Upwind diff
(  386)              PL(i,j,lb) = ((1.+fac0)*P0(i,j,L0) +  &
(  387)                   (1.-fac0)*P0(i,j,k))/2.
(  388) !                
(  389) !!      this is z-integral of actual f(z,t);
(  390)              Ibar = 0.0 +            &
(  391)                   (1.0/e3t(k)) * 1.0/kext *                           &
(  392)                   (                                  &
(  393)                   log( (exp(-kext*gdept(k))+         &
(  394)                   sqrt(C2(i,j)+exp(-2.*kext*gdept(k)))) / &






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 8

(  395)                   (exp(-kext*gdept(k+1))+            &
(  396)                   sqrt(C2(i,j)+exp(-2.*kext*gdept(k+1))))  ) &
(  397)                   ) 
(  398) !                
(  399)              P(i,j,k) = P(i,j,k) + DTIN * (          &
(  400)                   Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)    & 
(  401)                   - Z0(i,j,k)*Rm*(1.0-exp(-lambda*P0(i,j,k)))   & 
(  402)                   - M_PHY*P0(i,j,k)                  & 
(  403) !!DB: use  wmask here (computed in lib_bgcm_01)
(  404)                   -(wmask(i,j,k)*w_sink*PL(i,j,lt)-     &
(  405)                   wmask(i,j,L0)*w_sink*PL(i,j,lb))/e3t(k) )  !!sinking term
(  406) 
(  407)              Z(i,j,k) = Z(i,j,k) + DTIN * (           &
(  408)                   gamma*Rm*Z0(i,j,k)*(1.0-exp(-lambda*P0(i,j,k)))        &
(  409)                   - M_ZOO*Z0(i,j,k) )
(  410) !                
(  411)              N(i,j,k) = N(i,j,k) + DTIN * (            &
(  412)                   -Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)        &
(  413)                   +(1.0-gamma)*Rm*Z0(i,j,k)*           &
(  414)                   (1.0-exp(-lambda*P0(i,j,k)))         &
(  415)                   + M_PHY*P0(i,j,k) + M_ZOO*Z0(i,j,k) )
(  416)              !                
(  417)           enddo
(  418)        enddo
(  419)     
(  420)        i = lt          !!switch indices
(  421)        lt = lb
(  422)        lb = i
(  423)        
(  424)     enddo       !! k-loop
(  425) 
(  426) !!DB: special loop for jpk layer (in case it exists)
(  427) !!NB: w_sink code frag
(  428)     do k = jpk, jpk
(  429)        L0 = k + 1
(  430)        do j = 1, jpj
(  431)           do i = 1, jpi
(  432) 
(  433)              DTIN =  rdt * tmask(i,j,k)
(  434) !                
(  435) !!      this is z-integral of actual f(z,t); calc here as it will become f(x,y,z,t)
(  436)              Ibar = 0.0 +            &
(  437)                   (1.0/e3t(k)) * 1.0/kext *                           &
(  438)                   (                                  &
(  439)                   log( (exp(-kext*gdept(k))+         &
(  440)                   sqrt(C2(i,j)+exp(-2.*kext*gdept(k)))) / &
(  441)                   (exp(-kext*gdept(k+1))+            &
(  442)                   sqrt(C2(i,j)+exp(-2.*kext*gdept(k+1))))  ) &
(  443)                   ) 
(  444) !                
(  445)              P(i,j,k) = P(i,j,k) + DTIN * (          &
(  446)                   Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)    & 
(  447)                   - Z0(i,j,k)*Rm*(1.0-exp(-lambda*P0(i,j,k)))   & 
(  448)                   - M_PHY*P0(i,j,k)                  & 
(  449) !!DB: NEW  use wmask here (computed in lib_bgcm_01)
(  450)                   -(wmask(i,j,k)*w_sink*P0(i,j,k-1)-     &
(  451)                   0.0)/e3t(k) )  !!sinking term for bottom layer
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 9

(  453)              Z(i,j,k) = Z(i,j,k) + DTIN * (           &
(  454)                   gamma*Rm*Z0(i,j,k)*(1.0-exp(-lambda*P0(i,j,k)))        &
(  455)                   - M_ZOO*Z0(i,j,k) )
(  456) !                
(  457) 
(  458)              N(i,j,k) = N(i,j,k) + DTIN * (            &
(  459)                   -Vm*N0(i,j,k)/(ks+N0(i,j,k))*Ibar*P0(i,j,k)        &
(  460)                   +(1.0-gamma)*Rm*Z0(i,j,k)*           &
(  461)                   (1.0-exp(-lambda*P0(i,j,k)))         &
(  462)                   + M_PHY*P0(i,j,k) + M_ZOO*Z0(i,j,k) )
(  463)              !                
(  464)           enddo
(  465)        enddo
(  466)     
(  467)        i = lt          !!switch indices
(  468)        lt = lb
(  469)        lb = i
(  470)        
(  471)     enddo       !! k-loop
(  472) 
(  473) 
(  474) !!DB: 2009.04.27 add limiter to N,P,Z
(  475)     P(:,:,:) = max(P(:,:,:),1.e-4)
(  476)     Z(:,:,:) = max(Z(:,:,:),1.e-4)
(  477)     N(:,:,:) = max(N(:,:,:),1.e-4)
(  478) 
(  479)     
(  480)   END SUBROUTINE step_npz
(  481) 
(  482) 
(  483) END MODULE bgcm_01_model
(  484) 
(  485) #endif































PGF90 (Version     10.5)          09/23/2010  09:49:24      page 10

(    1) # 1 "bgcm_01_model.F90"
(    1) !!DB: 2009.08.17
(    2) !!Eliminated some routines that are no longer necessary. If needed look in
(    3) !!bgcm_01_model.F91 (if it exists locally) or some achived version of the code
(    4) !!NB: RE Open Boundary Conditions:
(    5) !!NB: trc_trp() calls trc_nxt() which calls update_boundary_vals()
(    6) !!which is the correct location for updating the tracer arrays with OB values
(    7) 
(    8) 
(    9) !!DB: 2009.06.09
(   10) !!Latest version which uses "copy-to-boundary" OBC scheme for PZ and 
(   11) !!N_obc routine for N
(   12) !!Note that without further modifying the code structure (which is not
(   13) !!necessarily a bad idea) --
(   14) !!  All calls related to assignment of OBCs must be done from trcnxt.F90
(   15) !!Thus look there for OBC scheme currently used
(   16) 
(   17) 
(   18) !!OLD comments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   19) !!ca. 2009.04
(   20) !!It also contains OBC-related fixes relative to previous versions
(   21) !!(see also TRP/trcnxt.F90 for OBC-related change)
(   22) !!DB 2008.08.26 ...
(   23) !!Replace trcstp with 1 simpler module specific to key_BGCM_01
(   24) !!See TRP/trcstp.F90 to see how this module is used (i.e. how some
(   25) !!backward compatibility is maintained)
(   26) 
(   27) !!The BioGeoChemical Model works by (1) calling a Source-Minus-Sink (SMS)
(   28) !!routine that does the vertical biophysics only; and (2) calling the
(   29) !!regular tracer advection/diffusion schemes. The tracer array(s) default 
(   30) !!to trn(:,:,:,jptra) (see trc.F90 and par_trc ----> par_trc_trp). 
(   31) !!jptra is the number of 3D BGCM fields used. 
(   32) 
(   33) 
(   34) !!Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   35) !!-1- The main BGCM is in bgcm_sms() -- at bottom of this module. It is currently
(   36) !!    v.simple. This is where the major modifications need to be done.
(   37) !!-2- The SMS routine would get its various model and runtime parameters in the module
(   38) !!    bgcm_01_initrc.F90. At this time there is only a place-holder for this
(   39) !!    code in that module. 
(   40) !!END Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   41) 
(   42) !!END OLD comments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   43) 
PGF90-W-0006-Input file empty (bgcm_01_model.F90)
