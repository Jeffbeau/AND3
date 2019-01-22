


PGF90 (Version     10.5)          09/23/2010  09:49:44      page 1

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

Filename: dynspg.F90

(    1) MODULE dynspg
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynspg  ***
(    4)    !! Ocean dynamics:  surface pressure gradient control
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_spg     : update the dynamics trend with the lateral diffusion
(    9)    !!   dyn_spg_ctl : initialization, namelist read, and parameters control
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce            ! ocean dynamics and tracers variables
(   13)    USE dom_oce        ! ocean space and time domain variables
(   14)    USE obc_oce        ! ocean open boundary conditions
(   15)    USE dynspg_oce     ! surface pressure gradient variables
(   16)    USE dynspg_exp     ! surface pressure gradient     (dyn_spg_exp routine)
(   17)    USE dynspg_ts      ! surface pressure gradient     (dyn_spg_ts  routine)
(   18)    USE dynspg_flt     ! surface pressure gradient     (dyn_spg_flt routine)
(   19)    USE dynspg_rl      ! surface pressure gradient     (dyn_spg_rl  routine)
(   20)    USE dynspg_exp_jki ! surface pressure gradient (dyn_spg_exp_jki routine)
(   21)    USE dynspg_ts_jki  ! surface pressure gradient (dyn_spg_ts_jki  routine)
(   22)    USE dynspg_flt_jki ! surface pressure gradient (dyn_spg_flt_jki routine)
(   23)    USE trdmod         ! ocean dynamics trends
(   24)    USE trdmod_oce     ! ocean variables trends
(   25)    USE prtctl         ! Print control                     (prt_ctl routine)
(   26)    USE in_out_manager ! I/O manager
(   27) 
(   28)    IMPLICIT NONE
(   29)    PRIVATE
(   30) 
(   31)    !! * Accessibility
(   32)    PUBLIC dyn_spg         ! routine called by step module
(   33) 
(   34)    !! * module variables
(   35)    INTEGER ::                        &
(   36)       nspg = 0                         ! type of surface pressure gradient scheme
(   37)       !                                ! defined from lk_dynspg_... 
(   38) 
(   39)    !! * Substitutions
(   40) #  include "domzgr_substitute.h90"
(   41) #  include "vectopt_loop_substitute.h90"
(   42)    !!----------------------------------------------------------------------
(   43)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   44)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg.F90,v 1.4 2005/12/29 10:51:26 opalod Exp $ 
(   45)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   46)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 2

(   47) 
(   48) CONTAINS
(   49) 
(   50)    SUBROUTINE dyn_spg( kt, kindic )
(   51)       !!----------------------------------------------------------------------
(   52)       !!                  ***  ROUTINE dyn_spg  ***
(   53)       !!
(   54)       !! ** Purpose :   compute the lateral ocean dynamics physics.
(   55)       !!
(   56)       !! History :
(   57)       !!   9.0  !  05-12  (C. Talandier, G. Madec)  Original code
(   58)       !!----------------------------------------------------------------------
(   59)       !! * Arguments
(   60)       INTEGER, INTENT( in  ) ::   kt     ! ocean time-step index
(   61)       INTEGER, INTENT( out ) ::   kindic ! solver flag
(   62) 
(   63)       !! * local declarations
(   64)       REAL(wp) ::   z2dt                      ! temporary scalar
(   65)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   66)          ztrdu, ztrdv                         ! 3D temporary workspace
(   67)       !!----------------------------------------------------------------------
(   68) 
(   69)       IF( kt == nit000 )   CALL dyn_spg_ctl      ! initialisation & control of options
(   70) 
(   71)       IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
(   72)          ztrdu(:,:,:) = ua(:,:,:)
(   73)          ztrdv(:,:,:) = va(:,:,:)
(   74)       ENDIF
(   75) 
(   76)       SELECT CASE ( nspg )                       ! compute surf. pressure gradient trend and add it to the general trend
(   77)       CASE ( -1 )                                       ! esopa: test all possibility with control print
(   78)          CALL dyn_spg_exp    ( kt )
(   79)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg0 - Ua: ', mask1=umask, &
(   80)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   81)          CALL dyn_spg_ts     ( kt )
(   82)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg1 - Ua: ', mask1=umask, &
(   83)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   84)          CALL dyn_spg_flt  ( kt, kindic )
(   85)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg2 - Ua: ', mask1=umask, &
(   86)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   87)          CALL dyn_spg_exp_jki( kt )
(   88)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg10- Ua: ', mask1=umask, &
(   89)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   90)          CALL dyn_spg_ts_jki ( kt )
(   91)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg12- Ua: ', mask1=umask, &
(   92)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   93)          CALL dyn_spg_flt_jki( kt, kindic )
(   94)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg13- Ua: ', mask1=umask, &
(   95)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   96)       CASE (  0 )                                       ! explicit
(   97)          CALL dyn_spg_exp    ( kt )
(   98)       CASE (  1 )                                       ! time-splitting
(   99)          CALL dyn_spg_ts     ( kt )
(  100)       CASE (  2 )                                       ! filtered
(  101)          CALL dyn_spg_flt    ( kt, kindic )
(  102)       CASE (  3 )                                       ! rigid lid
(  103)          CALL dyn_spg_rl     ( kt, kindic )
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 3

(  105)       CASE ( 10 )                                       ! explicit with j-k-i loop
(  106)          CALL dyn_spg_exp_jki( kt )
(  107)       CASE ( 11 )                                       ! time-splitting with j-k-i loop
(  108)          CALL dyn_spg_ts_jki ( kt )
(  109)       CASE ( 12 )                                       ! filtered with j-k-i loop
(  110)          CALL dyn_spg_flt_jki( kt, kindic )
(  111)       END SELECT
(  112) 
(  113)       !                                          ! save the horizontal diffusive trends for further diagnostics
(  114)       IF( l_trddyn )   THEN
(  115)          SELECT CASE ( nspg )
(  116)          CASE ( 0, 1, 3, 10, 11 )
(  117)             ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
(  118)             ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
(  119)          CASE( 2, 12 )
(  120)             z2dt = 2. * rdt
(  121)             IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt
(  122)             ztrdu(:,:,:) = ( ua(:,:,:) - ub(:,:,:) ) / z2dt - ztrdu(:,:,:)
(  123)             ztrdv(:,:,:) = ( va(:,:,:) - vb(:,:,:) ) / z2dt - ztrdv(:,:,:)
(  124)          END SELECT
(  125)          CALL trd_mod( ztrdu, ztrdv, jpdtdspg, 'DYN', kt )
(  126)       ENDIF
(  127) 
(  128)       !                                          ! print mean trends (used for debugging)
(  129)       IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg  - Ua: ', mask1=umask, &
(  130)          &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(  131) 
(  132)    END SUBROUTINE dyn_spg
(  133) 
(  134) 
(  135)    SUBROUTINE dyn_spg_ctl
(  136)       !!---------------------------------------------------------------------
(  137)       !!                  ***  ROUTINE dyn_spg_ctl  ***
(  138)       !!                
(  139)       !! ** Purpose :   Control the consistency between cpp options for 
(  140)       !!      surface pressure gradient schemes
(  141)       !!
(  142)       !! History :
(  143)       !!   9.0  !  05-10  (V. Garnier)  Original code : spg re-organization
(  144)       !!----------------------------------------------------------------------
(  145)       !! * Local declarations
(  146)       INTEGER ::   ioptio
(  147)       !!----------------------------------------------------------------------
(  148) 
(  149)       ! Parameter control and print
(  150)       ! ---------------------------
(  151)       ! Control print
(  152)       IF(lwp) THEN
(  153)          WRITE(numout,*)
(  154)          WRITE(numout,*) 'dyn_spg_ctl : choice of the surface pressure gradient scheme'
(  155)          WRITE(numout,*) '~~~~~~~~~~~'
(  156)          WRITE(numout,*) '     Explicit free surface                  lk_dynspg_exp = ', lk_dynspg_exp
(  157)          WRITE(numout,*) '     Free surface with time splitting       lk_dynspg_ts  = ', lk_dynspg_ts
(  158)          WRITE(numout,*) '     Filtered free surface cst volume       lk_dynspg_flt = ', lk_dynspg_flt
(  159)          WRITE(numout,*) '     Rigid-lid case                         lk_dynspg_rl  = ', lk_dynspg_rl
(  160)       ENDIF
(  161) 
(  162)       ! Control of surface pressure gradient scheme options






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 4

(  163)       ! ---------------------------------------------------
(  164)       ioptio = 0
(  165)       IF(lk_dynspg_exp)   ioptio = ioptio + 1
(  166)       IF(lk_dynspg_ts )   ioptio = ioptio + 1
(  167)       IF(lk_dynspg_flt)   ioptio = ioptio + 1
(  168)       IF(lk_dynspg_rl )   ioptio = ioptio + 1
(  169) 
(  170)       IF( ( ioptio > 1 .AND. .NOT. lk_esopa ) .OR. ioptio == 0 ) THEN
(  171)          IF(lwp) WRITE(numout,cform_err)
(  172)          IF(lwp) WRITE(numout,*) ' Choose only one surface pressure gradient scheme with a key cpp'
(  173)          nstop = nstop + 1
(  174)       ENDIF
(  175) 
(  176)       IF( lk_esopa     )   nspg = -1
(  177)       IF( lk_dynspg_exp)   nspg =  0
(  178)       IF( lk_dynspg_ts )   nspg =  1
(  179)       IF( lk_dynspg_flt)   nspg =  2
(  180)       IF( lk_dynspg_rl )   nspg =  3
(  181)       IF( lk_jki       )   nspg =  nspg + 10
(  182)       IF( nspg == 13   )   nspg =  3
(  183) 
(  184)       IF( lk_esopa     )   nspg = -1
(  185) 
(  186)      IF(lwp) THEN
(  187)          WRITE(numout,*)
(  188)          IF( nspg == -1 )   WRITE(numout,*) '     ESOPA test All scheme used except rigid-lid'
(  189)          IF( nspg ==  0 )   WRITE(numout,*) '     explicit free surface'
(  190)          IF( nspg ==  1 )   WRITE(numout,*) '     free surface with time splitting scheme'
(  191)          IF( nspg ==  2 )   WRITE(numout,*) '     filtered free surface'
(  192)          IF( nspg ==  3 )   WRITE(numout,*) '     rigid-lid'
(  193)          IF( nspg == 10 )   WRITE(numout,*) '     explicit free surface with j-k-i loop'
(  194)          IF( nspg == 11 )   WRITE(numout,*) '     time splitting free surface with j-k-i loop'
(  195)          IF( nspg == 12 )   WRITE(numout,*) '     filtered free surface with j-k-i loop'
(  196)       ENDIF
(  197) 
(  198)       ! Control of timestep choice
(  199)       ! --------------------------
(  200)       IF( lk_dynspg_ts ) THEN
(  201)          IF( MOD( rdt , rdtbt ) /= 0. ) THEN
(  202)             IF(lwp) WRITE(numout,cform_err)
(  203)             IF(lwp) WRITE(numout,*) ' The barotropic timestep must be an integer divisor of the baroclinic timestep'
(  204)             nstop = nstop + 1
(  205)          ENDIF
(  206)       ENDIF
(  207) 
(  208) #if key_obc
(  209)       ! Conservation of ocean volume (key_dynspg_flt)
(  210)       ! ---------------------------------------------
(  211)       IF( lk_dynspg_flt ) ln_vol_cst = .true.
(  212) 
(  213)       ! Application of Flather's algorithm at open boundaries
(  214)       ! -----------------------------------------------------
(  215)       IF( lk_dynspg_flt ) ln_obc_fla = .false.
(  216)       IF( lk_dynspg_exp ) ln_obc_fla = .true.
(  217)       IF( lk_dynspg_ts  ) ln_obc_fla = .true.
(  218) #endif
(  219) 
(  220)    END SUBROUTINE dyn_spg_ctl






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 5

(  221) 
(  222)   !!======================================================================
(  223) END MODULE dynspg





























































PGF90 (Version     10.5)          09/23/2010  09:49:44      page 6

(    1) # 1 "dynspg.F90"
(    1) MODULE dynspg
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynspg  ***
(    4)    !! Ocean dynamics:  surface pressure gradient control
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_spg     : update the dynamics trend with the lateral diffusion
(    9)    !!   dyn_spg_ctl : initialization, namelist read, and parameters control
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce            ! ocean dynamics and tracers variables
(   13)    USE dom_oce        ! ocean space and time domain variables
(   14)    USE obc_oce        ! ocean open boundary conditions
(   15)    USE dynspg_oce     ! surface pressure gradient variables
(   16)    USE dynspg_exp     ! surface pressure gradient     (dyn_spg_exp routine)
(   17)    USE dynspg_ts      ! surface pressure gradient     (dyn_spg_ts  routine)
(   18)    USE dynspg_flt     ! surface pressure gradient     (dyn_spg_flt routine)
(   19)    USE dynspg_rl      ! surface pressure gradient     (dyn_spg_rl  routine)
(   20)    USE dynspg_exp_jki ! surface pressure gradient (dyn_spg_exp_jki routine)
(   21)    USE dynspg_ts_jki  ! surface pressure gradient (dyn_spg_ts_jki  routine)
(   22)    USE dynspg_flt_jki ! surface pressure gradient (dyn_spg_flt_jki routine)
(   23)    USE trdmod         ! ocean dynamics trends
(   24)    USE trdmod_oce     ! ocean variables trends
(   25)    USE prtctl         ! Print control                     (prt_ctl routine)
(   26)    USE in_out_manager ! I/O manager
(   27) 
(   28)    IMPLICIT NONE
(   29)    PRIVATE
(   30) 
(   31)    !! * Accessibility
(   32)    PUBLIC dyn_spg         ! routine called by step module
(   33) 
(   34)    !! * module variables
(   35)    INTEGER ::                        &
(   36)       nspg = 0                         ! type of surface pressure gradient scheme
(   37)       !                                ! defined from lk_dynspg_... 
(   38) 
(   39)    !! * Substitutions
(   40) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 41 "dynspg.F90"
(   41) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 7

(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 42 "dynspg.F90"
(   42) # 42
(   42)    !!----------------------------------------------------------------------
(   43)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   44)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg.F90,v 1.4 2005/12/29 10:51:26 opalod Exp $ 
(   45)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   46)    !!----------------------------------------------------------------------
(   47) 
(   48) CONTAINS
(   49) 
(   50)    SUBROUTINE dyn_spg( kt, kindic )
(   51)       !!----------------------------------------------------------------------
(   52)       !!                  ***  ROUTINE dyn_spg  ***
(   53)       !!
(   54)       !! ** Purpose :   compute the lateral ocean dynamics physics.
(   55)       !!
(   56)       !! History :
(   57)       !!   9.0  !  05-12  (C. Talandier, G. Madec)  Original code
(   58)       !!----------------------------------------------------------------------
(   59)       !! * Arguments
(   60)       INTEGER, INTENT( in  ) ::   kt     ! ocean time-step index
(   61)       INTEGER, INTENT( out ) ::   kindic ! solver flag
(   62) 
(   63)       !! * local declarations
(   64)       REAL(wp) ::   z2dt                      ! temporary scalar
(   65)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   66)          ztrdu, ztrdv                         ! 3D temporary workspace
(   67)       !!----------------------------------------------------------------------
(   68) 
(   69)       IF( kt == nit000 )   CALL dyn_spg_ctl      ! initialisation & control of options
(   70) 
(   71)       IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
(   72)          ztrdu(:,:,:) = ua(:,:,:)
(   73)          ztrdv(:,:,:) = va(:,:,:)
(   74)       ENDIF
(   75) 
(   76)       SELECT CASE ( nspg )                       ! compute surf. pressure gradient trend and add it to the general trend
(   77)       CASE ( -1 )                                       ! esopa: test all possibility with control print
(   78)          CALL dyn_spg_exp    ( kt )
(   79)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg0 - Ua: ', mask1=umask, &
(   80)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   81)          CALL dyn_spg_ts     ( kt )
(   82)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg1 - Ua: ', mask1=umask, &
(   83)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   84)          CALL dyn_spg_flt  ( kt, kindic )
(   85)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg2 - Ua: ', mask1=umask, &
(   86)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   87)          CALL dyn_spg_exp_jki( kt )






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 8

(   88)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg10- Ua: ', mask1=umask, &
(   89)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   90)          CALL dyn_spg_ts_jki ( kt )
(   91)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg12- Ua: ', mask1=umask, &
(   92)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   93)          CALL dyn_spg_flt_jki( kt, kindic )
(   94)          IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg13- Ua: ', mask1=umask, &
(   95)                &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(   96)       CASE (  0 )                                       ! explicit
(   97)          CALL dyn_spg_exp    ( kt )
(   98)       CASE (  1 )                                       ! time-splitting
(   99)          CALL dyn_spg_ts     ( kt )
(  100)       CASE (  2 )                                       ! filtered
(  101)          CALL dyn_spg_flt    ( kt, kindic )
(  102)       CASE (  3 )                                       ! rigid lid
(  103)          CALL dyn_spg_rl     ( kt, kindic )
(  104) 
(  105)       CASE ( 10 )                                       ! explicit with j-k-i loop
(  106)          CALL dyn_spg_exp_jki( kt )
(  107)       CASE ( 11 )                                       ! time-splitting with j-k-i loop
(  108)          CALL dyn_spg_ts_jki ( kt )
(  109)       CASE ( 12 )                                       ! filtered with j-k-i loop
(  110)          CALL dyn_spg_flt_jki( kt, kindic )
(  111)       END SELECT
(  112) 
(  113)       !                                          ! save the horizontal diffusive trends for further diagnostics
(  114)       IF( l_trddyn )   THEN
(  115)          SELECT CASE ( nspg )
(  116)          CASE ( 0, 1, 3, 10, 11 )
(  117)             ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
(  118)             ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
(  119)          CASE( 2, 12 )
(  120)             z2dt = 2. * rdt
(  121)             IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt
(  122)             ztrdu(:,:,:) = ( ua(:,:,:) - ub(:,:,:) ) / z2dt - ztrdu(:,:,:)
(  123)             ztrdv(:,:,:) = ( va(:,:,:) - vb(:,:,:) ) / z2dt - ztrdv(:,:,:)
(  124)          END SELECT
(  125)          CALL trd_mod( ztrdu, ztrdv, jpdtdspg, 'DYN', kt )
(  126)       ENDIF
(  127) 
(  128)       !                                          ! print mean trends (used for debugging)
(  129)       IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg  - Ua: ', mask1=umask, &
(  130)          &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
(  131) 
(  132)    END SUBROUTINE dyn_spg
(  133) 
(  134) 
(  135)    SUBROUTINE dyn_spg_ctl
(  136)       !!---------------------------------------------------------------------
(  137)       !!                  ***  ROUTINE dyn_spg_ctl  ***
(  138)       !!                
(  139)       !! ** Purpose :   Control the consistency between cpp options for 
(  140)       !!      surface pressure gradient schemes
(  141)       !!
(  142)       !! History :
(  143)       !!   9.0  !  05-10  (V. Garnier)  Original code : spg re-organization
(  144)       !!----------------------------------------------------------------------
(  145)       !! * Local declarations






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 9

(  146)       INTEGER ::   ioptio
(  147)       !!----------------------------------------------------------------------
(  148) 
(  149)       ! Parameter control and print
(  150)       ! ---------------------------
(  151)       ! Control print
(  152)       IF(lwp) THEN
(  153)          WRITE(numout,*)
(  154)          WRITE(numout,*) 'dyn_spg_ctl : choice of the surface pressure gradient scheme'
(  155)          WRITE(numout,*) '~~~~~~~~~~~'
(  156)          WRITE(numout,*) '     Explicit free surface                  lk_dynspg_exp = ', lk_dynspg_exp
(  157)          WRITE(numout,*) '     Free surface with time splitting       lk_dynspg_ts  = ', lk_dynspg_ts
(  158)          WRITE(numout,*) '     Filtered free surface cst volume       lk_dynspg_flt = ', lk_dynspg_flt
(  159)          WRITE(numout,*) '     Rigid-lid case                         lk_dynspg_rl  = ', lk_dynspg_rl
(  160)       ENDIF
(  161) 
(  162)       ! Control of surface pressure gradient scheme options
(  163)       ! ---------------------------------------------------
(  164)       ioptio = 0
(  165)       IF(lk_dynspg_exp)   ioptio = ioptio + 1
(  166)       IF(lk_dynspg_ts )   ioptio = ioptio + 1
(  167)       IF(lk_dynspg_flt)   ioptio = ioptio + 1
(  168)       IF(lk_dynspg_rl )   ioptio = ioptio + 1
(  169) 
(  170)       IF( ( ioptio > 1 .AND. .NOT. lk_esopa ) .OR. ioptio == 0 ) THEN
(  171)          IF(lwp) WRITE(numout,cform_err)
(  172)          IF(lwp) WRITE(numout,*) ' Choose only one surface pressure gradient scheme with a key cpp'
(  173)          nstop = nstop + 1
(  174)       ENDIF
(  175) 
(  176)       IF( lk_esopa     )   nspg = -1
(  177)       IF( lk_dynspg_exp)   nspg =  0
(  178)       IF( lk_dynspg_ts )   nspg =  1
(  179)       IF( lk_dynspg_flt)   nspg =  2
(  180)       IF( lk_dynspg_rl )   nspg =  3
(  181)       IF( lk_jki       )   nspg =  nspg + 10
(  182)       IF( nspg == 13   )   nspg =  3
(  183) 
(  184)       IF( lk_esopa     )   nspg = -1
(  185) 
(  186)      IF(lwp) THEN
(  187)          WRITE(numout,*)
(  188)          IF( nspg == -1 )   WRITE(numout,*) '     ESOPA test All scheme used except rigid-lid'
(  189)          IF( nspg ==  0 )   WRITE(numout,*) '     explicit free surface'
(  190)          IF( nspg ==  1 )   WRITE(numout,*) '     free surface with time splitting scheme'
(  191)          IF( nspg ==  2 )   WRITE(numout,*) '     filtered free surface'
(  192)          IF( nspg ==  3 )   WRITE(numout,*) '     rigid-lid'
(  193)          IF( nspg == 10 )   WRITE(numout,*) '     explicit free surface with j-k-i loop'
(  194)          IF( nspg == 11 )   WRITE(numout,*) '     time splitting free surface with j-k-i loop'
(  195)          IF( nspg == 12 )   WRITE(numout,*) '     filtered free surface with j-k-i loop'
(  196)       ENDIF
(  197) 
(  198)       ! Control of timestep choice
(  199)       ! --------------------------
(  200)       IF( lk_dynspg_ts ) THEN
(  201)          IF( MOD( rdt , rdtbt ) /= 0. ) THEN
(  202)             IF(lwp) WRITE(numout,cform_err)
(  203)             IF(lwp) WRITE(numout,*) ' The barotropic timestep must be an integer divisor of the baroclinic timestep'






PGF90 (Version     10.5)          09/23/2010  09:49:44      page 10

(  204)             nstop = nstop + 1
(  205)          ENDIF
(  206)       ENDIF
(  207) 
(  208) # 209
(  209)       ! Conservation of ocean volume (key_dynspg_flt)
(  210)       ! ---------------------------------------------
(  211)       IF( lk_dynspg_flt ) ln_vol_cst = .true.
(  212) 
(  213)       ! Application of Flather's algorithm at open boundaries
(  214)       ! -----------------------------------------------------
(  215)       IF( lk_dynspg_flt ) ln_obc_fla = .false.
(  216)       IF( lk_dynspg_exp ) ln_obc_fla = .true.
(  217)       IF( lk_dynspg_ts  ) ln_obc_fla = .true.
(  218) 
(  219) # 220
(  220)    END SUBROUTINE dyn_spg_ctl
(  221) 
(  222)   !!======================================================================
(  223) END MODULE dynspg
