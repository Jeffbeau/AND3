


PGF90 (Version     10.5)          09/23/2010  09:49:23      page 1

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

Filename: trctrp.F90

(    1) MODULE trctrp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE trctrp  ***
(    4)    !! Ocean Physics    : manage the passive tracer transport
(    5)    !!======================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_trp        : passive tracer transport
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce_trc         ! ocean dynamics and active tracers variables
(   12)    USE trc             ! ocean passive tracers variables 
(   13) 
(   14) !   USE trctrp_lec      ! passive tracers transport
(   15) 
(   16)    USE trcbbl          ! bottom boundary layer               (trc_bbl routine)
(   17)    USE trcdmp          ! internal damping                    (trc_dmp routine)
(   18) 
(   19)    USE trcldf_bilapg   ! lateral mixing               (trc_ldf_bilapg routine)
(   20)    USE trcldf_bilap    ! lateral mixing                (trc_ldf_bilap routine)
(   21)    USE trcldf_iso      ! lateral mixing                  (trc_ldf_iso routine)
(   22)    USE trcldf_iso_zps  ! lateral mixing              (trc_ldf_iso_zps routine)
(   23)    USE trcldf_lap      ! lateral mixing                  (trc_ldf_lap routine)
(   24)  
(   25)    USE trcnxt          ! time-stepping                       (trc_nxt routine)
(   26)    USE trcrad          ! positivity                          (trc_rad routine)
(   27) 
(   28)    USE trcadv_cen2     ! 2nd order centered advection   (trc_adv_cen2 routine)
(   29)    USE trcadv_muscl    ! MUSCL advection               (trc_adv_muscl routine)
(   30)    USE trcadv_muscl2   ! MUSCL2 advection             (trc_adv_muscl2 routine)
(   31)    USE trcadv_tvd      ! TVD advection                   (trc_adv_tvd routine)
(   32)    USE trcadv_smolar   ! SMOLAR advection             (trc_adv_smolar routine)
(   33) 
(   34)    USE trczdf_exp      ! vertical diffusion              (trc_zdf_exp routine)
(   35)    USE trczdf_imp      ! vertical diffusion              (trc_zdf_exp routine)
(   36)    USE trczdf_iso      ! vertical diffusion              (trc_zdf_exp routine)
(   37)    USE trczdf_iso_vopt ! vertical diffusion              (trc_zdf_exp routine)
(   38)    USE trcsbc          ! surface boundary condition          (trc_sbc routine)
(   39) 
(   40)    USE zpshde_trc      ! partial step: hor. derivative   (zps_hde_trc routine)
(   41) 
(   42) #ifdef key_BGCM_01
(   43) !   USE bgcm_01_initrc
(   44)    USE lib_bgcm_01
(   45) #endif
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 2

(   47) 
(   48) 
(   49)    IMPLICIT NONE
(   50)    PRIVATE
(   51) 
(   52)    !! * Routine accessibility
(   53)    PUBLIC trc_trp            ! called by trc_stp
(   54) 
(   55)    !! * Substitutions
(   56) #  include "domzgr_substitute.h90"
(   57)    !!----------------------------------------------------------------------
(   58)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   59)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trctrp.F90,v 1.10 2006/04/11 13:49:00 opalod Exp $ 
(   60)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   61)    !!----------------------------------------------------------------------
(   62) 
(   63) CONTAINS
(   64) 
(   65)    SUBROUTINE trc_trp( kt )
(   66)       !!----------------------------------------------------------------------
(   67)       !!                     ***  ROUTINE trc_trp  ***
(   68)       !!                      
(   69)       !! ** Purpose : Management of passive tracers transport
(   70)       !! 
(   71)       !! ** Method  : 
(   72)       !!              Compute the passive tracers trends 
(   73)       !!              Update the passive tracers
(   74)       !!
(   75)       !! History :
(   76)       !!   9.0  !  04-03  (C. Ethe)  Original
(   77)       !!----------------------------------------------------------------------
(   78)       !! * Arguments
(   79)       INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
(   80)       !! ---------------------------------------------------------------------
(   81) 
(   82)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(   83)       ! Passive tracers
(   84)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(   85)       !-----------------------------------------------------------------------
(   86) 
(   87)       CALL trc_sbc( kt )            ! surface boundary condition
(   88) 
(   89) # if defined key_trcbbc
(   90)       IF(lwp) WRITE(numout,cform_err)
(   91)       IF(lwp) WRITE(numout,*) ' Bottom heat flux not yet implemented'
(   92)       IF(lwp) WRITE(numout,*) ' With passive tracers. '
(   93)       IF(lwp) WRITE(numout,*) ' Check trc_trp routine'
(   94)       nstop = nstop + 1
(   95) # endif 
(   96)       !                                                      ! bottom boundary condition
(   97)       IF( lk_trcbbl_dif    )   CALL trc_bbl_dif( kt )                ! diffusive bottom boundary layer scheme
(   98)       IF( lk_trcbbl_adv    )   CALL trc_bbl_adv( kt )                ! advective (and/or diffusive) bottom boundary layer sc
(   99) 
(  100)       IF( lk_trcdmp        )   CALL trc_dmp( kt )            ! internal damping trends
(  101) 
(  102) !!DB 2008.11.19 -- Enclose changes in ifdef just in case
(  103) #ifdef key_BGCM_01
(  104) !!DBG: 10.03 -- force this call as for some reason the tvd scheme --> nan






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 3

(  105) !!All of the below work, but cen2 seems to be slightly faster
(  106) !!DBG: no advection
(  107)       CALL trc_adv_cen2  ( kt )             ! 2nd order centered scheme -- works
(  108) 
(  109) !      CALL trc_adv_muscl ( kt )             ! MUSCL scheme -- works
(  110) !      CALL trc_adv_muscl2 ( kt )             ! MUSCL scheme -- works
(  111) #else
(  112) !!DB: NB same as above but other possibilities are more apparent
(  113)       !                                                      ! horizontal & vertical advection
(  114) !      IF( ln_trcadv_cen2   )   CALL trc_adv_cen2  ( kt )             ! 2nd order centered scheme
(  115) !      IF( ln_trcadv_muscl  )   CALL trc_adv_muscl ( kt )             ! MUSCL scheme
(  116) !      IF( ln_trcadv_muscl2 )   CALL trc_adv_muscl2( kt )             ! MUSCL2 scheme
(  117) !      IF( ln_trcadv_tvd    )   CALL trc_adv_tvd   ( kt )             ! TVD scheme
(  118) !      IF( ln_trcadv_smolar )   CALL trc_adv_smolar( kt )             ! SMOLARKIEWICZ scheme
(  119) 
(  120) !!DBG: 10.03 -- force this call as for some reason the tvd scheme --> nan
(  121) !!All of the below work, but cen2 seems to be slightly faster
(  122)       CALL trc_adv_cen2  ( kt )             ! 2nd order centered scheme -- works
(  123) !      CALL trc_adv_muscl ( kt )             ! MUSCL scheme -- works
(  124) !      CALL trc_adv_muscl2 ( kt )             ! MUSCL scheme -- works
(  125) #endif
(  126) 
(  127) 
(  128)       IF( n_cla == 1   )   THEN 
(  129)          IF(lwp) WRITE(numout,cform_err)
(  130)          IF(lwp) WRITE(numout,*) '          Cross Land Advection not yet implemented'
(  131)          IF(lwp) WRITE(numout,*) '          With Passive tracers. n_cla = ', n_cla
(  132)          IF(lwp) WRITE(numout,*) '          Check trc_trp routine'
(  133)          nstop = nstop + 1
(  134)       ENDIF
(  135) 
(  136)       !                                                      ! lateral mixing 
(  137)       IF( l_trcldf_bilapg  )   CALL trc_ldf_bilapg ( kt )            ! s-coord. horizontal bilaplacian
(  138)       IF( l_trcldf_bilap   )   CALL trc_ldf_bilap  ( kt )            ! iso-level bilaplacian 
(  139)       IF( l_trcldf_iso     )   CALL trc_ldf_iso    ( kt )            ! iso-neutral laplacian 
(  140)       IF( l_trcldf_iso_zps )   CALL trc_ldf_iso_zps( kt )            ! partial step iso-neutral laplacian
(  141)       IF( l_trcldf_lap     )   CALL trc_ldf_lap    ( kt )            ! iso-level laplacian
(  142) 
(  143) #ifdef key_BGCM_01
(  144) !!DB -- force implicit scheme
(  145)       call trc_zdf_imp( kt )                ! implicit time stepping (euler backward)
(  146) #else
(  147)       !                                                      ! vertical diffusion
(  148)       IF( l_trczdf_exp     )   CALL trc_zdf_exp( kt )                ! explicit time stepping (time splitting scheme)
(  149)       IF( l_trczdf_imp     )   CALL trc_zdf_imp( kt )                ! implicit time stepping (euler backward)
(  150)       IF( l_trczdf_iso     )   CALL trc_zdf_iso( kt )                ! isopycnal
(  151)       IF( l_trczdf_iso_vo  )   CALL trc_zdf_iso_vopt( kt )           ! vector opt. isopycnal
(  152) #endif
(  153) 
(  154)       CALL trc_nxt( kt )            ! tracer fields at next time step
(  155)  
(  156) !!DB: 2009.09.09 routine disabled
(  157) !      CALL trc_rad( kt )            ! Correct artificial negative concentrations for isopycnal scheme
(  158)       !                                                      
(  159) 
(  160)       IF( lk_zps .AND. .NOT. lk_trccfg_1d ) &
(  161)          &                     CALL zps_hde_trc( kt, trb, gtru, gtrv )  ! Partial steps: now horizontal gradient
(  162)       !                                                                 ! of passive tracers at the bottom ocean level






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 4

(  163) 
(  164) 
(  165)     END SUBROUTINE trc_trp
(  166) 
(  167) #else
(  168)    !!----------------------------------------------------------------------
(  169)    !!   Dummy module :                      NO passive tracers
(  170)    !!----------------------------------------------------------------------
(  171) CONTAINS
(  172)    SUBROUTINE trc_trp (kt )              ! Empty routine
(  173)       INTEGER, INTENT(in) :: kt
(  174) !      WRITE(*,*) 'trc_trp: You should not have seen this print! error?', kt
(  175)    END SUBROUTINE trc_trp
(  176) #endif
(  177)    
(  178)    !!======================================================================
(  179) END MODULE trctrp















































PGF90 (Version     10.5)          09/23/2010  09:49:23      page 5

(    1) # 1 "trctrp.F90"
(    1) MODULE trctrp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE trctrp  ***
(    4)    !! Ocean Physics    : manage the passive tracer transport
(    5)    !!======================================================================
(    6) # 168
(  168)    !!----------------------------------------------------------------------
(  169)    !!   Dummy module :                      NO passive tracers
(  170)    !!----------------------------------------------------------------------
(  171) CONTAINS
(  172)    SUBROUTINE trc_trp (kt )              ! Empty routine
(  173)       INTEGER, INTENT(in) :: kt
(  174) !      WRITE(*,*) 'trc_trp: You should not have seen this print! error?', kt
(  175)    END SUBROUTINE trc_trp
(  176) # 177
(  177)    
(  178)    !!======================================================================
(  179) END MODULE trctrp
