


PGF90 (Version     10.5)          09/23/2010  09:49:17      page 1

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

Filename: trcbbl.F90

(    1) MODULE trcbbl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcbbl  ***
(    4)    !! Ocean passive tracers physics :  advective and/or diffusive bottom boundary 
(    5)    !!                                  layer scheme
(    6)    !!==============================================================================
(    7) #if  defined key_passivetrc && ( defined key_trcbbl_dif   ||   defined key_trcbbl_adv ) && ! defined key_cfg_1d
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_trcbbl_dif'   or            diffusive bottom boundary layer
(   10)    !!   'key_trcbbl_adv'                 advective bottom boundary layer
(   11)    !!----------------------------------------------------------------------
(   12)    !!   trc_bbl_dif  : update the passive tracer trends due to the bottom
(   13)    !!                  boundary layer (diffusive only)
(   14)    !!   trc_bbl_adv  : update the passive tracer trends due to the bottom
(   15)    !!                  boundary layer (advective and/or diffusive)
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE oce_trc             ! ocean dynamics and active tracers variables
(   19)    USE trc                 ! ocean passive tracers variables
(   20)    USE prtctl_trc          ! Print control for debbuging
(   21)    USE eosbn2
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC trc_bbl_dif    ! routine called by step.F90
(   27)    PUBLIC trc_bbl_adv    ! routine called by step.F90
(   28) 
(   29)    !! * Shared module variables
(   30) # if defined key_trcbbl_dif
(   31)    LOGICAL, PUBLIC, PARAMETER ::    &  !:
(   32)       lk_trcbbl_dif = .TRUE.   !: advective bottom boundary layer flag
(   33) 
(   34) # else
(   35)    LOGICAL, PUBLIC, PARAMETER ::    &  !:
(   36)       lk_trcbbl_dif = .FALSE.  !: advective bottom boundary layer flag
(   37) # endif
(   38) 
(   39) # if defined key_trcbbl_adv
(   40)    LOGICAL, PUBLIC, PARAMETER ::    &  !:
(   41)       lk_trcbbl_adv = .TRUE.   !: advective bottom boundary layer flag
(   42)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   43)        u_trc_bbl, v_trc_bbl, &  !: velocity involved in exhanges in the advective BBL
(   44)        w_trc_bbl                !: vertical increment of velocity due to advective BBL
(   45)        !                        !  only affect tracer vertical advection
(   46) # else






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 2

(   47)    LOGICAL, PUBLIC, PARAMETER ::    &  !:
(   48)       lk_trcbbl_adv = .FALSE.  !: advective bottom boundary layer flag
(   49) # endif
(   50) 
(   51)    !! * Module variables
(   52)    INTEGER, DIMENSION(jpi,jpj) ::   &  !:
(   53)       mbkt, mbku, mbkv                 ! ???
(   54) 
(   55)    REAL(wp) ::        &  !!! * trcbbl namelist *
(   56)       atrcbbl = 1.e+3      ! lateral coeff. for bottom boundary layer scheme (m2/s)
(   57) 
(   58)    !! * Substitutions
(   59) #  include "passivetrc_substitute.h90"
(   60)    !!----------------------------------------------------------------------
(   61)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   62)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcbbl.F90,v 1.11 2006/04/11 13:49:00 opalod Exp $ 
(   63)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   64)    !!----------------------------------------------------------------------
(   65) 
(   66) CONTAINS
(   67) 
(   68)    SUBROUTINE trc_bbl_dif( kt )
(   69)       !!----------------------------------------------------------------------
(   70)       !!                   ***  ROUTINE trc_bbl_dif  ***
(   71)       !!
(   72)       !! ** Purpose :   Compute the before tracer trend associated 
(   73)       !!      with the bottom boundary layer and add it to the general trend 
(   74)       !!      of tracer equations. The bottom boundary layer is supposed to be
(   75)       !!      a purely diffusive bottom boundary layer.
(   76)       !!
(   77)       !! ** Method  :   When the product grad( rho) * grad(h) < 0 (where grad 
(   78)       !!      is an along bottom slope gradient) an additional lateral diffu-
(   79)       !!      sive trend along the bottom slope is added to the general tracer
(   80)       !!      trend, otherwise nothing is done.
(   81)       !!      Second order operator (laplacian type) with variable coefficient
(   82)       !!      computed as follow for temperature (idem on s): 
(   83)       !!         difft = 1/(e1t*e2t*e3t) { di-1[ ahbt e2u*e3u/e1u di[ztb] ]
(   84)       !!                                 + dj-1[ ahbt e1v*e3v/e2v dj[ztb] ] }
(   85)       !!      where ztb is a 2D array: the bottom ocean temperature and ahtb
(   86)       !!      is a time and space varying diffusive coefficient defined by:
(   87)       !!         ahbt = zahbp    if grad(rho).grad(h) < 0
(   88)       !!              = 0.       otherwise.
(   89)       !!      Note that grad(.) is the along bottom slope gradient. grad(rho)
(   90)       !!      is evaluated using the local density (i.e. referenced at the
(   91)       !!      local depth). Typical value of ahbt is 2000 m2/s (equivalent to
(   92)       !!      a downslope velocity of 20 cm/s if the condition for slope
(   93)       !!      convection is satified)
(   94)       !!      Add this before trend to the general trend tra of the 
(   95)       !!      botton ocean tracer point:
(   96)       !!         tra = tra + difft
(   97)       !!
(   98)       !! ** Action  : - update tra at the bottom level with the bottom
(   99)       !!                boundary layer trend
(  100)       !!
(  101)       !! References :
(  102)       !!     Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
(  103)       !!
(  104)       !! History :






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 3

(  105)       !!   8.0  !  96-06  (L. Mortier)  Original code
(  106)       !!   8.0  !  97-11  (G. Madec)  Optimization
(  107)       !!   8.5  !  02-08  (G. Madec)  free form + modules
(  108)       !!   9.0  !  04-03  (C. Ethe)   Adaptation for passive tracers
(  109)       !!----------------------------------------------------------------------
(  110)       !! * Arguments 
(  111)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(  112) 
(  113)       !! * Local declarations
(  114)       INTEGER ::   ji, jj,jn                ! dummy loop indices
(  115)       INTEGER ::   ik
(  116)       INTEGER ::   ii0, ii1, ij0, ij1       ! temporary integers
(  117) #  if defined key_partial_steps
(  118)       INTEGER  ::   iku1, iku2, ikv1,ikv2   ! temporary intergers
(  119)       REAL(wp) ::   ze3u, ze3v              ! temporary scalars
(  120) #  else
(  121)       INTEGER ::   iku, ikv
(  122) #  endif
(  123)       REAL(wp) ::   &
(  124)          zsign, zt, zs, zh, zalbet,      &  ! temporary scalars
(  125)          zgdrho, zbtr, ztra
(  126)       REAL(wp), DIMENSION(jpi,jpj) ::    &
(  127)         zki, zkj, zkx, zky,    &  ! temporary workspace arrays
(  128)         ztnb, zsnb, zdep,                &
(  129)         ztrb, zahu, zahv
(  130)       CHARACTER (len=22) :: charout
(  131)       REAL(wp) ::   &
(  132)          fsalbt, pft, pfs, pfh              ! statement function
(  133)       !!----------------------------------------------------------------------
(  134)       ! ratio alpha/beta
(  135)       ! ================
(  136)       !  fsalbt: ratio of thermal over saline expension coefficients
(  137)       !       pft :  potential temperature in degrees celcius
(  138)       !       pfs :  salinity anomaly (s-35) in psu
(  139)       !       pfh :  depth in meters
(  140) 
(  141)       fsalbt( pft, pfs, pfh ) =                                              &
(  142)          ( ( ( -0.255019e-07 * pft + 0.298357e-05 ) * pft                    &
(  143)                                    - 0.203814e-03 ) * pft                    &
(  144)                                    + 0.170907e-01 ) * pft                    &
(  145)                                    + 0.665157e-01                            &
(  146)          +(-0.678662e-05 * pfs - 0.846960e-04 * pft + 0.378110e-02 ) * pfs   &
(  147)          +  ( ( - 0.302285e-13 * pfh                                         &
(  148)                 - 0.251520e-11 * pfs                                         &
(  149)                 + 0.512857e-12 * pft * pft          ) * pfh                  &
(  150)                                      - 0.164759e-06   * pfs                  &
(  151)              +(   0.791325e-08 * pft - 0.933746e-06 ) * pft                  &
(  152)                                      + 0.380374e-04 ) * pfh   
(  153)       !!----------------------------------------------------------------------
(  154) 
(  155) 
(  156)       IF( kt == nittrc000 )   CALL trc_bbl_init
(  157) 
(  158) 
(  159)       ! 0. 2D fields of bottom temperature and salinity, and bottom slope
(  160)       ! -----------------------------------------------------------------
(  161)       ! mbathy= number of w-level, minimum value=1 (cf dommsk.F)
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 4

(  163) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  164)       jj = 1
(  165)       DO ji = 1, jpij   ! vector opt. (forced unrolling)
(  166) #  else
(  167)       DO jj = 1, jpj
(  168)          DO ji = 1, jpi
(  169) #  endif
(  170)             ik = mbkt(ji,jj)                              ! index of the bottom ocean T-level
(  171)             ztnb(ji,jj) = tn(ji,jj,ik) * tmask(ji,jj,1)   ! masked now T and S at ocean bottom 
(  172)             zsnb(ji,jj) = sn(ji,jj,ik) * tmask(ji,jj,1)
(  173)             zdep(ji,jj) = fsdept(ji,jj,ik)                ! depth of the ocean bottom T-level
(  174) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  175)          END DO
(  176) #  endif
(  177)       END DO
(  178) 
(  179) #  if defined key_partial_steps
(  180)       ! partial steps correction 
(  181) #   if defined key_vectopt_loop   &&   ! defined key_autotasking
(  182)       jj = 1
(  183)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  184) #   else
(  185)       DO jj = 1, jpjm1
(  186)          DO ji = 1, jpim1
(  187) #   endif
(  188)             iku1 = MAX( mbathy(ji+1,jj  )-1, 1 )
(  189)             iku2 = MAX( mbathy(ji  ,jj  )-1, 1 )
(  190)             ikv1 = MAX( mbathy(ji  ,jj+1)-1, 1 )
(  191)             ikv2 = MAX( mbathy(ji  ,jj  )-1, 1 )
(  192)             ze3u = MIN( fse3u(ji,jj,iku1), fse3u(ji,jj,iku2) ) 
(  193)             ze3v = MIN( fse3v(ji,jj,ikv1), fse3v(ji,jj,ikv2) ) 
(  194)             zahu(ji,jj) = atrcbbl * e2u(ji,jj) * ze3u / e1u(ji,jj) * umask(ji,jj,1)
(  195)             zahv(ji,jj) = atrcbbl * e1v(ji,jj) * ze3v / e2v(ji,jj) * vmask(ji,jj,1)
(  196) #   if ! defined key_vectopt_loop   ||   defined key_autotasking
(  197)          END DO
(  198) #   endif
(  199)       END DO
(  200) #  else
(  201) #   if defined key_vectopt_loop   &&   ! defined key_autotasking
(  202)       jj = 1
(  203)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  204) #   else
(  205)       DO jj = 1, jpjm1
(  206)          DO ji = 1, jpim1
(  207) #   endif
(  208)             iku = mbku(ji,jj)
(  209)             ikv = mbkv(ji,jj)
(  210)             zahu(ji,jj) = atrcbbl * e2u(ji,jj) * fse3u(ji,jj,iku) / e1u(ji,jj) * umask(ji,jj,1)
(  211)             zahv(ji,jj) = atrcbbl * e1v(ji,jj) * fse3v(ji,jj,ikv) / e2v(ji,jj) * vmask(ji,jj,1)
(  212) #   if ! defined key_vectopt_loop   ||   defined key_autotasking
(  213)          END DO
(  214) #   endif
(  215)       END DO
(  216) #  endif
(  217) 
(  218) !!
(  219) !!     OFFLINE VERSION OF DIFFUSIVE BBL
(  220) !!






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 5

(  221) #if defined key_off_tra
(  222) 
(  223)       ! 2. Additional second order diffusive trends
(  224)       ! -------------------------------------------
(  225) 
(  226)       DO jn = 1, jptra
(  227)          ! first derivative (gradient)
(  228)          
(  229) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  230)          jj = 1
(  231)          DO ji = 1, jpij   ! vector opt. (forced unrolling)
(  232) #  else
(  233)          DO jj = 1, jpj
(  234)             DO ji = 1, jpi
(  235) #  endif
(  236)                ik = mbkt(ji,jj) 
(  237)                ztrb(ji,jj) = trb(ji,jj,ik,jn) * tmask(ji,jj,1)
(  238) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  239)             END DO
(  240) #  endif
(  241)          END DO
(  242) 
(  243) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  244)          jj = 1
(  245)          DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  246) #  else
(  247)          DO jj = 1, jpjm1
(  248)             DO ji = 1, jpim1
(  249) #  endif
(  250)                zkx(ji,jj) = bblx(ji,jj) * zahu(ji,jj) * ( ztrb(ji+1,jj) - ztrb(ji,jj) )
(  251)                zky(ji,jj) = bbly(ji,jj) * zahv(ji,jj) * ( ztrb(ji,jj+1) - ztrb(ji,jj) )
(  252) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  253)             END DO
(  254) #  endif
(  255)          END DO
(  256) !!
(  257) !!  ONLINE VERSION OF DIFFUSIVE BBL
(  258) !!
(  259) #else
(  260)       ! 1. Criteria of additional bottom diffusivity: grad(rho).grad(h)<0
(  261)       ! --------------------------------------------
(  262)       ! Sign of the local density gradient along the i- and j-slopes
(  263)       ! multiplied by the slope of the ocean bottom
(  264) 	SELECT CASE ( neos )
(  265) 
(  266)       	CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  267) 
(  268) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  269)       jj = 1
(  270)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  271) #  else
(  272)       DO jj = 1, jpjm1
(  273)          DO ji = 1, jpim1
(  274) #  endif
(  275)             ! temperature, salinity anomalie and depth
(  276)             zt = 0.5 * ( ztnb(ji,jj) + ztnb(ji+1,jj) )
(  277)             zs = 0.5 * ( zsnb(ji,jj) + zsnb(ji+1,jj) ) - 35.0
(  278)             zh = 0.5 * ( zdep(ji,jj) + zdep(ji+1,jj) )






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 6

(  279)             ! masked ratio alpha/beta
(  280)             zalbet = fsalbt( zt, zs, zh )*umask(ji,jj,1)
(  281)             ! local density gradient along i-bathymetric slope
(  282)             zgdrho = zalbet * ( ztnb(ji+1,jj) - ztnb(ji,jj) )   &
(  283)                    -          ( zsnb(ji+1,jj) - zsnb(ji,jj) )
(  284)             ! sign of local i-gradient of density multiplied by the i-slope
(  285)             zsign = SIGN( 0.5, - zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
(  286)             zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
(  287) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  288)          END DO
(  289) #  endif
(  290)       END DO
(  291) 
(  292) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  293)       jj = 1
(  294)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  295) #  else
(  296)       DO jj = 1, jpjm1
(  297)          DO ji = 1, jpim1
(  298) #  endif
(  299)             ! temperature, salinity anomalie and depth
(  300)             zt = 0.5 * ( ztnb(ji,jj+1) + ztnb(ji,jj) )
(  301)             zs = 0.5 * ( zsnb(ji,jj+1) + zsnb(ji,jj) ) - 35.0
(  302)             zh = 0.5 * ( zdep(ji,jj+1) + zdep(ji,jj) )
(  303)             ! masked ratio alpha/beta
(  304)             zalbet = fsalbt( zt, zs, zh )*vmask(ji,jj,1)
(  305)             ! local density gradient along j-bathymetric slope
(  306)             zgdrho = zalbet * ( ztnb(ji,jj+1) - ztnb(ji,jj) )   &
(  307)                    -          ( zsnb(ji,jj+1) - zsnb(ji,jj) )
(  308)             ! sign of local j-gradient of density multiplied by the j-slope
(  309)             zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
(  310)             zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
(  311) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  312)          END DO
(  313) #  endif
(  314)       END DO
(  315) 
(  316) 	CASE ( 1 )               ! Linear formulation function of temperature only
(  317) 
(  318) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  319)       jj = 1
(  320)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  321) #  else
(  322)       DO jj = 1, jpjm1
(  323)          DO ji = 1, jpim1
(  324) #  endif
(  325)             ! local density gradient along i-bathymetric slope
(  326)             zgdrho =  ( ztnb(ji+1,jj) - ztnb(ji,jj) )
(  327)             ! sign of local i-gradient of density multiplied by the i-slope
(  328)             zsign = SIGN( 0.5, - zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
(  329)             zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
(  330) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  331)          END DO
(  332) #  endif
(  333)       END DO
(  334) 
(  335) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  336)       jj = 1






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 7

(  337)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  338) #  else
(  339)       DO jj = 1, jpjm1
(  340)          DO ji = 1, jpim1
(  341) #  endif
(  342)             ! local density gradient along j-bathymetric slope
(  343)             zgdrho =  ( ztnb(ji,jj+1) - ztnb(ji,jj) )
(  344)             ! sign of local j-gradient of density multiplied by the j-slope
(  345)             zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
(  346)             zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
(  347) 
(  348) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  349)          END DO
(  350) #  endif
(  351)       END DO
(  352) 
(  353)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  354) 
(  355)       DO jj = 1, jpjm1
(  356)         DO ji = 1, fs_jpim1   ! vector opt.
(  357)             ! local density gradient along i-bathymetric slope
(  358)             zgdrho = - ( rbeta*( zsnb(ji+1,jj) - zsnb(ji,jj) )   &
(  359)                      -  ralpha*( ztnb(ji+1,jj) - ztnb(ji,jj) ) )
(  360)             ! sign of local i-gradient of density multiplied by the i-slope
(  361)             zsign = SIGN( 0.5, - zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
(  362)  	    zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
(  363)         END DO
(  364)       END DO
(  365) 
(  366)       DO jj = 1, jpjm1
(  367)         DO ji = 1, fs_jpim1   ! vector opt.
(  368)             ! local density gradient along j-bathymetric slope
(  369)             zgdrho = - ( rbeta*( zsnb(ji,jj+1) - zsnb(ji,jj) )   &
(  370)                      -  ralpha*( ztnb(ji,jj+1) - ztnb(ji,jj) ) )
(  371)             ! sign of local j-gradient of density multiplied by the j-slope
(  372)             zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
(  373)             zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
(  374)          END DO
(  375)       END DO
(  376) 
(  377) 
(  378)       CASE DEFAULT
(  379) 
(  380)          IF(lwp) WRITE(numout,cform_err)
(  381)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  382)          nstop = nstop + 1
(  383) 
(  384)       END SELECT
(  385)       
(  386)       ! 2. Additional second order diffusive trends
(  387)       ! -------------------------------------------
(  388) 
(  389)       DO jn = 1, jptra
(  390)          ! first derivative (gradient)
(  391) 
(  392) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  393)          jj = 1
(  394)          DO ji = 1, jpij   ! vector opt. (forced unrolling)






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 8

(  395) #  else
(  396)          DO jj = 1, jpj
(  397)             DO ji = 1, jpi
(  398) #  endif
(  399)                ik = mbkt(ji,jj)
(  400)                ztrb(ji,jj) = trb(ji,jj,ik,jn) * tmask(ji,jj,1)
(  401) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  402)             END DO
(  403) #  endif
(  404)          END DO
(  405) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  406)          jj = 1
(  407)          DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  408) #  else
(  409)          DO jj = 1, jpjm1
(  410)             DO ji = 1, jpim1
(  411) #  endif
(  412)                zkx(ji,jj) = zki(ji,jj) * ( ztrb(ji+1,jj) - ztrb(ji,jj) )
(  413)                zky(ji,jj) = zkj(ji,jj) * ( ztrb(ji,jj+1) - ztrb(ji,jj) )
(  414) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  415)             END DO
(  416) #  endif
(  417)          END DO
(  418) #endif
(  419) 
(  420) !!DB: deleted ORCA
(  421) !         IF( cp_cfg == "orca" ) THEN
(  422)          
(  423)          ! second derivative (divergence) and add to the general tracer trend
(  424) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  425)          jj = 1
(  426)          DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  427) #  else
(  428)          DO jj = 2, jpjm1
(  429)             DO ji = 2, jpim1
(  430) #  endif
(  431)                ik = MAX( mbathy(ji,jj)-1, 1 )
(  432)                zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,ik) )
(  433)                ztra = (  zkx(ji,jj) - zkx(ji-1,jj  )    &
(  434)                   &    + zky(ji,jj) - zky(ji  ,jj-1)  ) * zbtr
(  435)                tra(ji,jj,ik,jn) = tra(ji,jj,ik,jn) + ztra
(  436) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  437)             END DO
(  438) #  endif
(  439)          END DO
(  440) 
(  441)       END DO
(  442) 
(  443)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  444)          WRITE(charout, FMT="('bbl - dif')")
(  445)          CALL prt_ctl_trc_info(charout)
(  446)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  447)       ENDIF
(  448) 
(  449)    END SUBROUTINE trc_bbl_dif
(  450) 
(  451) # if defined key_trcbbl_adv
(  452)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 9

(  453)    !!   'key_trcbbl_adv'                    advective bottom boundary layer
(  454)    !!----------------------------------------------------------------------
(  455) #  include "trcbbl_adv.h90"
(  456) # else
(  457)    !!----------------------------------------------------------------------
(  458)    !!   Default option :                 NO advective bottom boundary layer
(  459)    !!----------------------------------------------------------------------
(  460)    SUBROUTINE trc_bbl_adv (kt )              ! Empty routine
(  461)       INTEGER, INTENT(in) :: kt
(  462) !      WRITE(*,*) 'trc_bbl_adv: You should not have seen this print! error?', kt
(  463)    END SUBROUTINE trc_bbl_adv
(  464) # endif
(  465) 
(  466)    SUBROUTINE trc_bbl_init
(  467)       !!----------------------------------------------------------------------
(  468)       !!                  ***  ROUTINE trc_bbl_init  ***
(  469)       !!
(  470)       !! ** Purpose :   Initialization for the bottom boundary layer scheme.
(  471)       !!
(  472)       !! ** Method  :   Read the namtrcbbl namelist and check the parameters
(  473)       !!      called by tra_bbl at the first timestep (nittrc000)
(  474)       !!
(  475)       !! History :
(  476)       !!    8.5  !  02-08  (G. Madec)  Original code
(  477)       !!----------------------------------------------------------------------
(  478)       !! * Local declarations
(  479)       INTEGER ::   ji, jj      ! dummy loop indices
(  480)       INTEGER :: numnat=80
(  481)       NAMELIST/namtrcbbl/ atrcbbl
(  482) 
(  483)       !!----------------------------------------------------------------------
(  484)       ! Read Namelist namtrcbbl : bottom boundary layer scheme
(  485)       ! --------------------
(  486) 
(  487)       OPEN(numnat,FILE='namelist.trp.cfc')
(  488)       REWIND ( numnat )
(  489)       READ   ( numnat, namtrcbbl )
(  490)       CLOSE(numnat)
(  491) 
(  492) 
(  493)       ! Parameter control and print
(  494)       ! ---------------------------
(  495)       IF(lwp) THEN
(  496)          WRITE(numout,*)
(  497)          WRITE(numout,*) 'trc_bbl_init : * Diffusive Bottom Boundary Layer'
(  498)          WRITE(numout,*) '~~~~~~~~~~~~'
(  499)          WRITE(numout,*) ' bottom boundary layer coef.    atrcbbl = ', atrcbbl
(  500) # if defined key_trcbbl_adv
(  501)             WRITE(numout,*) '               * Advective Bottom Boundary Layer'
(  502) # endif
(  503)          WRITE(numout,*)
(  504)       ENDIF
(  505)  
(  506)       DO jj = 1, jpj
(  507)          DO ji = 1, jpi
(  508)             mbkt(ji,jj) = MAX( mbathy(ji,jj) - 1, 1 )   ! vertical index of the bottom ocean T-level
(  509)          END DO
(  510)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 10

(  511)       DO jj = 1, jpjm1
(  512)          DO ji = 1, jpim1
(  513)             mbku(ji,jj) = MAX( MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) ) - 1, 1 )
(  514)             mbkv(ji,jj) = MAX( MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) ) - 1, 1 )
(  515)          END DO
(  516)       END DO
(  517) !!bug ???
(  518) !!bug Caution : define the vakue of mbku & mbkv everywhere!!! but lbc mpp lnk : pb when closed (0)
(  519) 
(  520) # if defined key_trcbbl_adv
(  521)       w_trc_bbl(:,:,:) = 0.e0    ! initialisation of w_trc_bbl to zero
(  522) # endif
(  523) 
(  524)    END SUBROUTINE trc_bbl_init
(  525) 
(  526) #else
(  527)    !!----------------------------------------------------------------------
(  528)    !!   Dummy module :                      No bottom boundary layer scheme
(  529)    !!----------------------------------------------------------------------
(  530)    LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbl_dif = .FALSE.   !: diff bbl flag
(  531)    LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbl_adv = .FALSE.   !: adv  bbl flag
(  532) CONTAINS
(  533)    SUBROUTINE trc_bbl_dif (kt )              ! Empty routine
(  534)       INTEGER, INTENT(in) :: kt
(  535) !      WRITE(*,*) 'trc_bbl_dif: You should not have seen this print! error?', kt
(  536)    END SUBROUTINE trc_bbl_dif
(  537)    SUBROUTINE trc_bbl_adv (kt )              ! Empty routine
(  538)       INTEGER, INTENT(in) :: kt
(  539) !      WRITE(*,*) 'trc_bbl_adv: You should not have seen this print! error?', kt
(  540)    END SUBROUTINE trc_bbl_adv
(  541) #endif
(  542) 
(  543)    !!======================================================================
(  544) END MODULE trcbbl






























PGF90 (Version     10.5)          09/23/2010  09:49:17      page 11

(    1) # 1 "trcbbl.F90"
(    1) MODULE trcbbl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcbbl  ***
(    4)    !! Ocean passive tracers physics :  advective and/or diffusive bottom boundary 
(    5)    !!                                  layer scheme
(    6)    !!==============================================================================
(    7) # 527
(  527)    !!----------------------------------------------------------------------
(  528)    !!   Dummy module :                      No bottom boundary layer scheme
(  529)    !!----------------------------------------------------------------------
(  530)    LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbl_dif = .FALSE.   !: diff bbl flag
(  531)    LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbl_adv = .FALSE.   !: adv  bbl flag
(  532) CONTAINS
(  533)    SUBROUTINE trc_bbl_dif (kt )              ! Empty routine
(  534)       INTEGER, INTENT(in) :: kt
(  535) !      WRITE(*,*) 'trc_bbl_dif: You should not have seen this print! error?', kt
(  536)    END SUBROUTINE trc_bbl_dif
(  537)    SUBROUTINE trc_bbl_adv (kt )              ! Empty routine
(  538)       INTEGER, INTENT(in) :: kt
(  539) !      WRITE(*,*) 'trc_bbl_adv: You should not have seen this print! error?', kt
(  540)    END SUBROUTINE trc_bbl_adv
(  541) 
(  542) # 543
(  543)    !!======================================================================
(  544) END MODULE trcbbl
