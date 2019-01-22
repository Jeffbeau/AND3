


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

Filename: dynspg_exp.F90

(    1) MODULE dynspg_exp
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_exp  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) #if ( defined key_dynspg_exp && ! defined key_autotasking ) ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_dynspg_exp'       free sfce cst vol. without filter nor ts
(    9)    !!   NOT 'key_autotasking'                      k-j-i loop (vector opt.)
(   10)    !!----------------------------------------------------------------------
(   11)    !!   dyn_spg_exp  : update the momentum trend with the surface 
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
(   32)    PUBLIC dyn_spg_exp  ! routine called by step.F90
(   33) 
(   34)    !! * Substitutions
(   35) #  include "domzgr_substitute.h90"
(   36) #  include "vectopt_loop_substitute.h90"
(   37)    !!----------------------------------------------------------------------
(   38)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_exp.F90,v 1.4 2005/12/29 10:44:29 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------
(   42) 
(   43) 
(   44) CONTAINS
(   45) 
(   46)    SUBROUTINE dyn_spg_exp( kt )






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 2

(   47)       !!----------------------------------------------------------------------
(   48)       !!                  ***  routine dyn_spg_exp  ***
(   49)       !!
(   50)       !! ** Purpose :   Compute the now trend due to the surface pressure
(   51)       !!      gradient in case of explicit free surface formulation and 
(   52)       !!      add it to the general trend of momentum equation. Compute
(   53)       !!      the free surface.
(   54)       !!
(   55)       !! ** Method  :   Explicit free surface formulation. The surface pressure
(   56)       !!      gradient is given by:
(   57)       !!         spgu = 1/rau0 d/dx(ps) =  g/e1u di( sshn )
(   58)       !!         spgv = 1/rau0 d/dy(ps) =  g/e2v dj( sshn )
(   59)       !!      -1- Compute the now surface pressure gradient
(   60)       !!      -2- Add it to the general trend
(   61)       !!      -3- Compute the horizontal divergence of velocities
(   62)       !!      - the now divergence is given by :
(   63)       !!         zhdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
(   64)       !!      - integrate the horizontal divergence from the bottom 
(   65)       !!         to the surface
(   66)       !!      - apply lateral boundary conditions on zhdivn
(   67)       !!      -4- Estimate the after sea surface elevation from the kinematic
(   68)       !!         surface boundary condition:
(   69)       !!         zssha = sshb - 2 rdt ( zhdiv + emp )
(   70)       !!      - Time filter applied on now sea surface elevation to avoid
(   71)       !!         the divergence of two consecutive time-steps and swap of free
(   72)       !!         surface arrays to start the next time step:
(   73)       !!         sshb = sshn + atfp * [ sshb + zssha - 2 sshn ]
(   74)       !!         sshn = zssha
(   75)       !!      - apply lateral boundary conditions on sshn
(   76)       !!
(   77)       !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
(   78)       !!
(   79)       !! References :
(   80)       !!
(   81)       !! History :
(   82)       !!   9.0  !  05-11  (V. Garnier, G. Madec, L. Bessieres) Original code
(   83)       !!
(   84)       !!---------------------------------------------------------------------
(   85)       !! * Arguments
(   86)       INTEGER, INTENT( in )  ::   kt         ! ocean time-step index
(   87) 
(   88)       !! * Local declarations
(   89)       INTEGER  ::   ji, jj, jk               ! dummy loop indices
(   90)       REAL(wp) ::   z2dt, zraur, zssha       ! temporary scalars 
(   91)       REAL(wp), DIMENSION(jpi,jpj)    ::  &  ! temporary arrays
(   92)          &         zhdiv
(   93)       !!----------------------------------------------------------------------
(   94) 
(   95)       IF( kt == nit000 ) THEN
(   96)          IF(lwp) WRITE(numout,*)
(   97)          IF(lwp) WRITE(numout,*) 'dyn_spg_exp : surface pressure gradient trend'
(   98)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   (explicit free surface)'
(   99) 
(  100)          ! set to zero free surface specific arrays
(  101)          spgu(:,:) = 0.e0                     ! surface pressure gradient (i-direction)
(  102)          spgv(:,:) = 0.e0                     ! surface pressure gradient (j-direction)
(  103)       ENDIF
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 3

(  105)       ! 0. Initialization
(  106)       ! -----------------
(  107)       ! read or estimate sea surface height and vertically integrated velocities
(  108)       IF( lk_obc )   CALL obc_dta_bt( kt, 0 )
(  109)       z2dt = 2. * rdt                                       ! time step: leap-frog
(  110)       IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt       ! time step: Euler if restart from rest
(  111)       zraur = 1.e0 / rauw
(  112) 
(  113) 
(  114)       ! 1. Surface pressure gradient (now)
(  115)       ! ----------------------------
(  116)       DO jj = 2, jpjm1
(  117)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  118)             spgu(ji,jj) = - grav * ( sshn(ji+1,jj) - sshn(ji,jj) ) / e1u(ji,jj)
(  119)             spgv(ji,jj) = - grav * ( sshn(ji,jj+1) - sshn(ji,jj) ) / e2v(ji,jj)
(  120)          END DO
(  121)       END DO
(  122) 
(  123)       ! 2. Add the surface pressure trend to the general trend
(  124)       ! ------------------------------------------------------
(  125)       DO jk = 1, jpkm1
(  126)          DO jj = 2, jpjm1
(  127)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  128)                ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
(  129)                va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
(  130)             END DO
(  131)          END DO
(  132)       END DO
(  133)      
(  134)       ! 3. Vertical integration of horizontal divergence of velocities
(  135)       ! --------------------------------
(  136)       zhdiv(:,:) = 0.e0
(  137)       DO jk = jpkm1, 1, -1
(  138)          DO jj = 2, jpjm1
(  139)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  140)                zhdiv(ji,jj) = zhdiv(ji,jj) + (  e2u(ji  ,jj  ) * fse3u(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)      &
(  141)                   &                           - e2u(ji-1,jj  ) * fse3u(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)      &
(  142)                   &                           + e1v(ji  ,jj  ) * fse3v(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)      &
(  143)                   &                           - e1v(ji  ,jj-1) * fse3v(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)  )   &
(  144)                   &                        / ( e1t(ji,jj) * e2t(ji,jj) )
(  145)             END DO
(  146)          END DO
(  147)       END DO
(  148) 
(  149) #if defined key_obc
(  150)       ! open boundaries (div must be zero behind the open boundary)
(  151)       !  mpp remark: The zeroing of zhdiv can probably be extended to 1->jpi/jpj for the correct row/column
(  152)       IF( lp_obc_east  )   zhdiv(nie0p1:nie1p1,nje0  :nje1  ) = 0.e0      ! east
(  153)       IF( lp_obc_west  )   zhdiv(niw0  :niw1  ,njw0  :njw1  ) = 0.e0      ! west
(  154)       IF( lp_obc_north )   zhdiv(nin0  :nin1  ,njn0p1:njn1p1) = 0.e0      ! north
(  155)       IF( lp_obc_south )   zhdiv(nis0  :nis1  ,njs0  :njs1  ) = 0.e0      ! south
(  156) #endif
(  157) 
(  158) 
(  159)       ! 4. Sea surface elevation time stepping
(  160)       ! --------------------------------------
(  161)       ! Euler (forward) time stepping, no time filter
(  162)       IF( neuler == 0 .AND. kt == nit000 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:40      page 4

(  163)          DO jj = 1, jpj
(  164)             DO ji = 1, jpi
(  165)                ! after free surface elevation
(  166)                zssha = sshb(ji,jj) - rdt * ( zraur * emp(ji,jj) + zhdiv(ji,jj) ) * tmask(ji,jj,1)
(  167)                ! swap of arrays
(  168)                sshb(ji,jj) = sshn(ji,jj)
(  169)                sshn(ji,jj) = zssha
(  170)             END DO
(  171)          END DO
(  172)       ELSE
(  173)          ! Leap-frog time stepping and time filter
(  174)          DO jj = 1, jpj
(  175)             DO ji = 1, jpi
(  176)                ! after free surface elevation
(  177)                zssha = sshb(ji,jj) - z2dt * ( zraur * emp(ji,jj) + zhdiv(ji,jj) ) * tmask(ji,jj,1)
(  178)                ! time filter and array swap
(  179)                sshb(ji,jj) = atfp * ( sshb(ji,jj) + zssha ) + atfp1 * sshn(ji,jj)
(  180)                sshn(ji,jj) = zssha
(  181)             END DO
(  182)          END DO
(  183)       ENDIF
(  184) 
(  185)       ! Boundary conditions on sshn
(  186)       IF( .NOT. lk_obc ) CALL lbc_lnk( sshn, 'T', 1. )
(  187)  
(  188)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  189)          CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
(  190)       ENDIF
(  191)       
(  192)    END SUBROUTINE dyn_spg_exp
(  193) 
(  194) #else
(  195)    !!----------------------------------------------------------------------
(  196)    !!   Default case :   Empty module   No standart explicit free surface 
(  197)    !!----------------------------------------------------------------------
(  198)    USE in_out_manager
(  199) CONTAINS
(  200)    SUBROUTINE dyn_spg_exp( kt )       ! Empty routine
(  201)       if(lwp) WRITE(numout,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
(  202)    END SUBROUTINE dyn_spg_exp
(  203) #endif
(  204)    
(  205)    !!======================================================================
(  206) END MODULE dynspg_exp




















PGF90 (Version     10.5)          09/23/2010  09:49:40      page 5

(    1) # 1 "dynspg_exp.F90"
(    1) MODULE dynspg_exp
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dynspg_exp  ***
(    4)    !! Ocean dynamics:  surface pressure gradient trend
(    5)    !!======================================================================
(    6) # 195
(  195)    !!----------------------------------------------------------------------
(  196)    !!   Default case :   Empty module   No standart explicit free surface 
(  197)    !!----------------------------------------------------------------------
(  198)    USE in_out_manager
(  199) CONTAINS
(  200)    SUBROUTINE dyn_spg_exp( kt )       ! Empty routine
(  201)       if(lwp) WRITE(numout,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
(  202)    END SUBROUTINE dyn_spg_exp
(  203) # 204
(  204)    
(  205)    !!======================================================================
(  206) END MODULE dynspg_exp
