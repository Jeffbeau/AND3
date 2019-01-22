


PGF90 (Version     10.5)          09/23/2010  09:49:40      page 1

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

Filename: dynspg_exp_jki.F90

(    1) MODULE dynspg_exp_jki
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_exp_jki  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) #if ( defined key_dynspg_exp && defined key_autotasking )   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_dynspg_exp'                              explicit free surface
(    9)    !!   'key_autotasking'                                        j-k-i loop
(   10)    !!----------------------------------------------------------------------
(   11)    !!   dyn_spg_exp_jki  : update the momentum trend with the surface 
(   12)    !!                      pressure gradient in the free surface constant  
(   13)    !!                      volume case with vector optimization
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers 
(   17)    USE dom_oce         ! ocean space and time domain 
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE phycst          ! physical constants
(   20)    USE ocesbc          ! ocean surface boundary condition
(   21)    USE obc_oce         ! Lateral open boundary condition
(   22)    USE obc_par         ! open boundary condition parameters
(   23)    USE obcdta          ! open boundary condition data     (obc_dta_bt routine)
(   24)    USE lib_mpp         ! distributed memory computing library
(   25)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   26)    USE prtctl          ! Print control
(   27) 
(   28)    IMPLICIT NONE
(   29)    PRIVATE
(   30) 
(   31)    !! * Accessibility
(   32)    PUBLIC dyn_spg_exp_jki  ! routine called by step.F90
(   33) 
(   34)    !! * Substitutions
(   35) #  include "domzgr_substitute.h90"
(   36) #  include "vectopt_loop_substitute.h90"
(   37)    !!----------------------------------------------------------------------
(   38)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_exp_jki.F90,v 1.1 2005/12/29 10:44:29 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------
(   42) 
(   43) CONTAINS
(   44) 
(   45)    SUBROUTINE dyn_spg_exp_jki( kt )
(   46)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 2

(   47)       !!                  ***  routine dyn_spg_exp_jki  ***
(   48)       !!
(   49)       !! ** Purpose :   Compute the now trend due to the surface pressure
(   50)       !!      gradient in case of explicit free surface formulation and 
(   51)       !!      add it to the general trend of momentum equation. Compute
(   52)       !!      the free surface.
(   53)       !!
(   54)       !! ** Method  :   Explicit free surface formulation. The surface pressure
(   55)       !!      gradient is given by:
(   56)       !!         spgu = 1/rau0 d/dx(ps) =  g/e1u di( sshn )
(   57)       !!         spgv = 1/rau0 d/dy(ps) =  g/e2v dj( sshn )
(   58)       !!      -1- Compute the now surface pressure gradient
(   59)       !!      -2- Add it to the general trend
(   60)       !!      -3- Compute the horizontal divergence of velocities
(   61)       !!      - the now divergence is given by :
(   62)       !!         zhdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
(   63)       !!      - integrate the horizontal divergence from the bottom 
(   64)       !!         to the surface
(   65)       !!      - apply lateral boundary conditions on zhdivn
(   66)       !!      -4- Estimate the after sea surface elevation from the kinematic
(   67)       !!         surface boundary condition:
(   68)       !!         zssha = sshb - 2 rdt ( zhdiv + emp )
(   69)       !!      - Time filter applied on now sea surface elevation to avoid
(   70)       !!         the divergence of two consecutive time-steps and swap of free
(   71)       !!         surface arrays to start the next time step:
(   72)       !!         sshb = sshn + atfp * [ sshb + zssha - 2 sshn ]
(   73)       !!         sshn = zssha
(   74)       !!      - apply lateral boundary conditions on sshn
(   75)       !!
(   76)       !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
(   77)       !!
(   78)       !! References :
(   79)       !!
(   80)       !! History :
(   81)       !!   9.0  !  05-11  (V. Garnier, G. Madec, L. Bessieres) Original code
(   82)       !!
(   83)       !!---------------------------------------------------------------------
(   84)       !! * Arguments
(   85)       INTEGER, INTENT( in )  ::   kt         ! ocean time-step index
(   86) 
(   87)       !! * Local declarations
(   88)       INTEGER  ::   ji, jj, jk               ! dummy loop indices
(   89)       REAL(wp) ::   z2dt, zraur, zssha       ! temporary scalars 
(   90)       REAL(wp), DIMENSION(jpi,jpj)    ::  &  ! temporary arrays
(   91)          &         zhdiv
(   92)       !!----------------------------------------------------------------------
(   93) 
(   94)       IF( kt == nit000 ) THEN
(   95)          IF(lwp) WRITE(numout,*)
(   96)          IF(lwp) WRITE(numout,*) 'dyn_spg_exp_jki : surface pressure gradient trend'
(   97)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   (explicit free surface, j-k-i loop)'
(   98) 
(   99)          ! set to zero free surface specific arrays
(  100)          spgu(:,:) = 0.e0                     ! surface pressure gradient (i-direction)
(  101)          spgv(:,:) = 0.e0                     ! surface pressure gradient (j-direction)
(  102)       ENDIF
(  103) 
(  104)       ! 0. Local constant initialization






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 3

(  105)       ! --------------------------------
(  106)       ! read or estimate sea surface height and vertically integrated velocities
(  107)       IF( lk_obc )   CALL obc_dta_bt( kt, 0 )
(  108)       z2dt = 2. * rdt                                       ! time step: leap-frog
(  109)       IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt       ! time step: Euler if restart from rest
(  110)       zraur = 1.e0 / rauw
(  111) 
(  112) !CDIR PARALLEL DO
(  113) !$OMP PARALLEL DO
(  114)       !                                                ! =============== !
(  115)       DO jj = 2, jpjm1                                 !  Vertical slab  !
(  116)          !                                             ! =============== !
(  117) 
(  118)          ! Surface pressure gradient (now)
(  119)          DO ji = 2, jpim1
(  120)             spgu(ji,jj) = - grav * ( sshn(ji+1,jj) - sshn(ji,jj) ) / e1u(ji,jj)
(  121)             spgv(ji,jj) = - grav * ( sshn(ji,jj+1) - sshn(ji,jj) ) / e2v(ji,jj)
(  122)          END DO
(  123) 
(  124)          ! Add the surface pressure trend to the general trend
(  125)          DO jk = 1, jpkm1
(  126)             DO ji = 2, jpim1
(  127)                ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
(  128)                va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
(  129)             END DO
(  130)          END DO
(  131)      
(  132)          !  Vertical integration of horizontal divergence of velocities
(  133)          zhdiv(:,jj) = 0.e0
(  134)          DO jk = jpkm1, 1, -1
(  135)             DO ji = 2, jpim1
(  136)                zhdiv(ji,jj) = zhdiv(ji,jj) + (  e2u(ji  ,jj  ) * fse3u(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)      &
(  137)                   &                           - e2u(ji-1,jj  ) * fse3u(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)      &
(  138)                   &                           + e1v(ji  ,jj  ) * fse3v(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)      &
(  139)                   &                           - e1v(ji  ,jj-1) * fse3v(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)  )   &
(  140)                   &                        / ( e1t(ji,jj) * e2t(ji,jj) )
(  141)             END DO
(  142)          END DO
(  143) 
(  144) #if defined key_obc
(  145)          ! open boundaries (div must be zero behind the open boundary)
(  146)          !  mpp remark: The zeroing of zhdiv can probably be extended to 1->jpi/jpj for the correct row/column
(  147)          IF( lp_obc_east  ) THEN
(  148)             IF( nje0   <= jj .AND. jj <= nje1   )   zhdiv(nie0p1:nie1p1,jj) = 0.e0      ! east
(  149)          ENDIF
(  150)          IF( lp_obc_west  ) THEN
(  151)             IF( njw0   <= jj .AND. jj <= njw1   )   zhdiv(niw0  :niw1  ,jj) = 0.e0      ! west
(  152)          ENDIF
(  153)          IF( lp_obc_north ) THEN
(  154)             IF( njn0p1 <= jj .AND. jj <= njn1p1 )   zhdiv(nin0  :nin1  ,jj) = 0.e0      ! north
(  155)          ENDIF
(  156)          IF( lp_obc_south ) THEN
(  157)             IF( njs0   <= jj .AND. jj <= njs1   )   zhdiv(nis0  :nis1  ,jj) = 0.e0      ! south
(  158)          ENDIF
(  159) #endif
(  160) 
(  161)          ! Sea surface elevation time stepping
(  162)          IF( neuler == 0 .AND. kt == nit000 ) THEN       ! Euler (forward) time stepping, no time filter






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 4

(  163)             DO ji = 2, jpim1
(  164)                ! after free surface elevation
(  165)                zssha = sshb(ji,jj) - rdt * ( zraur * emp(ji,jj) + zhdiv(ji,jj) ) * tmask(ji,jj,1)
(  166)                ! swap of arrays
(  167)                sshb(ji,jj) = sshn(ji,jj)
(  168)                sshn(ji,jj) = zssha
(  169)             END DO
(  170)          ELSE                                            ! Leap-frog time stepping and time filter
(  171)             DO ji = 2, jpim1
(  172)                ! after free surface elevation
(  173)                zssha = sshb(ji,jj) - z2dt * ( zraur * emp(ji,jj) + zhdiv(ji,jj) ) * tmask(ji,jj,1)
(  174)                ! time filter and array swap
(  175)                sshb(ji,jj) = atfp * ( sshb(ji,jj) + zssha ) + atfp1 * sshn(ji,jj)
(  176)                sshn(ji,jj) = zssha
(  177)             END DO
(  178)          ENDIF
(  179)          !                                             ! =============== !
(  180)       END DO                                           !     end slab    !
(  181)       !                                                ! =============== !
(  182) 
(  183)       ! Boundary conditions on sshn
(  184)       IF( .NOT. lk_obc ) CALL lbc_lnk( sshn, 'T', 1. )
(  185)  
(  186)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  187)          CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
(  188)       ENDIF
(  189)       
(  190)    END SUBROUTINE dyn_spg_exp_jki
(  191) 
(  192) #else
(  193)    !!----------------------------------------------------------------------
(  194)    !!   Default case :   Empty module   No standart explicit free surface 
(  195)    !!----------------------------------------------------------------------
(  196)    USE in_out_manager
(  197) CONTAINS
(  198)    SUBROUTINE dyn_spg_exp_jki( kt )       ! Empty routine
(  199)       if(lwp) WRITE(numout,*) 'dyn_spg_exp_jki: You should not have seen this print! error?', kt
(  200)    END SUBROUTINE dyn_spg_exp_jki
(  201) #endif
(  202)    
(  203)    !!======================================================================
(  204) END MODULE dynspg_exp_jki






















PGF90 (Version     10.5)          09/23/2010  09:49:40      page 5

(    1) # 1 "dynspg_exp_jki.F90"
(    1) MODULE dynspg_exp_jki
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_exp_jki  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) # 193
(  193)    !!----------------------------------------------------------------------
(  194)    !!   Default case :   Empty module   No standart explicit free surface 
(  195)    !!----------------------------------------------------------------------
(  196)    USE in_out_manager
(  197) CONTAINS
(  198)    SUBROUTINE dyn_spg_exp_jki( kt )       ! Empty routine
(  199)       if(lwp) WRITE(numout,*) 'dyn_spg_exp_jki: You should not have seen this print! error?', kt
(  200)    END SUBROUTINE dyn_spg_exp_jki
(  201) # 202
(  202)    
(  203)    !!======================================================================
(  204) END MODULE dynspg_exp_jki
