


PGF90 (Version     10.5)          09/23/2010  09:49:35      page 1

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

Filename: dyncor1d.F90

(    1) MODULE dyncor1d
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  ini1D  ***
(    4)    !! Ocean state   :  1D initialization
(    5)    !!=====================================================================
(    6) #if defined key_cfg_1d
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_cfg_1d'               1D Configuration
(    9)    !!----------------------------------------------------------------------
(   10)    !!----------------------------------------------------------------------
(   11)    !!   fcorio_1d   : Coriolis factor at T-point
(   12)    !!   dyn_cor_1d  : vorticity trend due to Coriolis
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce            ! ocean dynamics and tracers
(   16)    USE dom_oce        ! ocean space and time domain
(   17)    USE phycst         ! physical constants
(   18)    USE in_out_manager ! I/O manager
(   19)    USE prtctl         ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC fcorio_1d   ! routine called by OPA.F90
(   26)    PUBLIC dyn_cor_1d  ! routine called by step1d.F90
(   27) 
(   28)    !! * Substitutions
(   29) #  include "vectopt_loop_substitute.h90"
(   30)    !!----------------------------------------------------------------------
(   31)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/dyncor1d.F90,v 1.3 2005/10/03 09:20:35 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE fcorio_1d
(   39)       !!----------------------------------------------------------------------
(   40)       !!                   ***  ROUTINE fcorio_1d  ***
(   41)       !! 
(   42)       !! ** Purpose : Compute the Coriolis factor at T-point
(   43)       !!
(   44)       !! ** Method  :
(   45)       !!
(   46)       !! History :






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 2

(   47)       !!   9.0  !  04-09  (C. Ethe) 1D configuration
(   48)       !!----------------------------------------------------------------------
(   49)       !! * Local declarations
(   50)       !!----------------------------------------------------------------------
(   51)       REAL(wp) ::   &
(   52)          zphi0, zbeta, zf0         !  temporary scalars
(   53)  
(   54) 
(   55)       !!----------------------------------------------------------------------
(   56) 
(   57)       ! ================= !
(   58)       !  Coriolis factor  !
(   59)       ! ================= !
(   60)       IF(lwp) WRITE(numout,*)
(   61)       IF(lwp) WRITE(numout,*) 'fcorio_1d : Coriolis factor at T-point'
(   62)       IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
(   63) 
(   64)       SELECT CASE( jphgr_msh )   ! type of horizontal mesh
(   65) 
(   66)       CASE ( 0, 1, 4 )               ! mesh on the sphere
(   67) 
(   68)          ff(:,:) = 2. * omega * SIN( rad * gphit(:,:) ) 
(   69) 
(   70)       CASE ( 2 )                     ! f-plane at ppgphi0 
(   71) 
(   72)          ff(:,:) = 2. * omega * SIN( rad * ppgphi0 )
(   73) 
(   74)          IF(lwp) WRITE(numout,*) '          f-plane: Coriolis parameter = constant = ', ff(1,1)
(   75) 
(   76)       CASE ( 3 )                     ! beta-plane
(   77) 
(   78)          zbeta   = 2. * omega * COS( rad * ppgphi0 ) / ra                     ! beta at latitude ppgphi0
(   79)          zphi0   = ppgphi0 - FLOAT( jpjglo/2) * ppe2_m *1.e-3  / ( ra * rad ) ! latitude of the first row F-points
(   80)          zf0     = 2. * omega * SIN( rad * zphi0 )                            ! compute f0 1st point south
(   81) 
(   82)          ff(:,:) = ( zf0  + zbeta * gphit(:,:) * 1.e+3 )                      ! f = f0 +beta* y ( y=0 at south)
(   83) 
(   84)          IF(lwp) WRITE(numout,*) '          Beta-plane: Beta parameter = constant = ', ff(1,1)
(   85)          IF(lwp) WRITE(numout,*) '                      Coriolis parameter varies from ', ff(1,1),' to ', ff(1,jpj)
(   86) 
(   87)       CASE ( 5 )                     ! beta-plane and rotated domain
(   88) 
(   89)          zbeta = 2. * omega * COS( rad * ppgphi0 ) / ra                     ! beta at latitude ppgphi0
(   90)          zphi0 = 15.e0                                                      ! latitude of the first row F-points
(   91)          zf0   = 2. * omega * SIN( rad * zphi0 )                            ! compute f0 1st point south
(   92) 
(   93)          ff(:,:) = ( zf0 + zbeta * ABS( gphit(:,:) - zphi0 ) * rad * ra )   ! f = f0 +beta* y ( y=0 at south)
(   94) 
(   95)          IF(lwp) WRITE(numout,*) '          Beta-plane: Beta parameter = constant = ', ff(1,1)
(   96)          IF(lwp) WRITE(numout,*) '                      Coriolis parameter varies from ', ff(1,1),' to ', ff(1,jpj)
(   97) 
(   98)       END SELECT
(   99) 
(  100)    END SUBROUTINE fcorio_1d
(  101) 
(  102) 
(  103)    SUBROUTINE dyn_cor_1d( kt )
(  104)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 3

(  105)       !!                   ***  ROUTINE dyn_cor_1d  ***
(  106)       !! 
(  107)       !! ** Purpose :   Compute the now total vorticity trend and add it to 
(  108)       !!               the general trend of the momentum equation
(  109)       !!
(  110)       !! ** Method  :
(  111)       !!
(  112)       !! History :
(  113)       !!   9.0  !  04-09  (C. Ethe) 1D configuration
(  114)       !!----------------------------------------------------------------------
(  115)       !! * Arguments
(  116)       INTEGER, INTENT( in ) ::   kt        ! ocean time-step index
(  117) 
(  118)       !! * Local declarations
(  119)       INTEGER ::   ji, jj, jk              ! dummy loop indices
(  120)       REAL(wp) ::   &
(  121)          zua, zva                          ! temporary scalars
(  122) 
(  123)       !!----------------------------------------------------------------------
(  124) 
(  125)       IF( kt == nit000 ) THEN
(  126)          IF(lwp) WRITE(numout,*)
(  127)          IF(lwp) WRITE(numout,*) 'dyn_cor_1d : total vorticity trend in 1D'
(  128)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
(  129)       ENDIF
(  130) 
(  131)       DO jk = 1, jpkm1
(  132)          DO jj = 2, jpjm1
(  133)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  134)                zua =    ff(ji,jj) * vn(ji,jj,jk)
(  135)                zva =  - ff(ji,jj) * un(ji,jj,jk)
(  136)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  137)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  138)             END DO
(  139)          END DO
(  140)       END DO   
(  141) 
(  142)       IF(ln_ctl)   THEN
(  143)          CALL prt_ctl(tab3d_1=ua, clinfo1=' cor  - Ua: ', mask1=umask, &
(  144)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask)
(  145)       ENDIF
(  146) 
(  147) !     IF(l_ctl) THEN         ! print sum trends (used for debugging)
(  148) !        zua = SUM( ua(2:nictl,2:njctl,1:jpkm1) * umask(2:nictl,2:njctl,1:jpkm1) )
(  149) !        zva = SUM( va(2:nictl,2:njctl,1:jpkm1) * vmask(2:nictl,2:njctl,1:jpkm1) )
(  150) !        WRITE(numout,*) ' cor  - Ua: ', zua-u_ctl, ' Va: ', zva-v_ctl
(  151) !        u_ctl = zua   ;   v_ctl = zva
(  152) !     ENDIF
(  153) 
(  154)    END SUBROUTINE dyn_cor_1d
(  155) 
(  156) #else
(  157)    !!----------------------------------------------------------------------
(  158)    !!   Default key                                     NO 1D Config
(  159)    !!----------------------------------------------------------------------
(  160)    USE in_out_manager
(  161) CONTAINS
(  162)    SUBROUTINE fcorio_1d      ! Empty routine






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 4

(  163)    END SUBROUTINE fcorio_1d   
(  164)    SUBROUTINE dyn_cor_1d ( kt )
(  165)       if(lwp)WRITE(numout,*) 'dyn_cor_1d: You should not have seen this print! error?', kt
(  166)    END SUBROUTINE dyn_cor_1d
(  167) #endif
(  168) 
(  169)    !!=====================================================================
(  170) END MODULE dyncor1d
























































PGF90 (Version     10.5)          09/23/2010  09:49:35      page 5

(    1) # 1 "dyncor1d.F90"
(    1) MODULE dyncor1d
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  ini1D  ***
(    4)    !! Ocean state   :  1D initialization
(    5)    !!=====================================================================
(    6) # 157
(  157)    !!----------------------------------------------------------------------
(  158)    !!   Default key                                     NO 1D Config
(  159)    !!----------------------------------------------------------------------
(  160)    USE in_out_manager
(  161) CONTAINS
(  162)    SUBROUTINE fcorio_1d      ! Empty routine
(  163)    END SUBROUTINE fcorio_1d   
(  164)    SUBROUTINE dyn_cor_1d ( kt )
(  165)       if(lwp)WRITE(numout,*) 'dyn_cor_1d: You should not have seen this print! error?', kt
(  166)    END SUBROUTINE dyn_cor_1d
(  167) 
(  168) # 169
(  169)    !!=====================================================================
(  170) END MODULE dyncor1d
