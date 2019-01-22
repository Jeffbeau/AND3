


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

Filename: zpshde_trc.F90

(    1) MODULE zpshde_trc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE zpshde_trc   ***
(    4)    !! Ocean passive tracers: 
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc && ( defined key_partial_steps || defined key_esopa )
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_partial_steps' :               z-coordinate with partial steps
(    9)    !!----------------------------------------------------------------------
(   10)    !!   zps_hde_trc  :  Horizontal DErivative of passive tracers at the last
(   11)    !!                   ocean level (Z-coord. with Partial Steps)
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce_trc         ! ocean dynamics and tracers variables
(   15)    USE trc             ! ocean passive tracers variables
(   16)    USE lbclnk          ! lateral boundary conditions (or mpp link)
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC zps_hde_trc          ! routine called by step.F90
(   23) 
(   24)    !! * module variables
(   25)    INTEGER, DIMENSION(jpi,jpj) ::   &
(   26)       mbatu, mbatv      ! bottom ocean level index at U- and V-points
(   27) 
(   28)    !! * Substitutions
(   29) #  include "passivetrc_substitute.h90"
(   30)    !!----------------------------------------------------------------------
(   31)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/zpshde_trc.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE zps_hde_trc ( kt, ptra, pgtru, pgtrv )
(   39)       !!----------------------------------------------------------------------
(   40)       !!                     ***  ROUTINE zps_hde_trc  ***
(   41)       !!                    
(   42)       !! ** Purpose :   Compute the horizontal derivative of passive tracers
(   43)       !!      TRA at u- and v-points with a linear interpolation for z-coordinate
(   44)       !!      with partial steps.
(   45)       !!
(   46)       !! ** Method  :   In z-coord with partial steps, scale factors on last 






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 2

(   47)       !!      levels are different for each grid point, so that TRA points 
(   48)       !!      are not at the same depth as in z-coord. To have horizontal
(   49)       !!      gradients again, we interpolate TRA at the good depth : 
(   50)       !!      Linear interpolation of TRA  
(   51)       !!         Computation of di(trb) and dj(trb) by vertical interpolation:
(   52)       !!          di(tra) = tra~ - tra(i,j,k) or tra(i+1,j,k) - tra~
(   53)       !!          dj(tra) = tra~ - tra(i,j,k) or tra(i,j+1,k) - tra~
(   54)       !!         This formulation computes the two cases:
(   55)       !!                 CASE 1                   CASE 2  
(   56)       !!         k-1  ___ ___________   k-1   ___ ___________
(   57)       !!                  TRAi  TRA~             TRA~  TRAi+1
(   58)       !!                  _____                        _____
(   59)       !!         k        |   |TRAi+1   k         TRAi |   |
(   60)       !!                  |   |____                ____|   |
(   61)       !!              ___ |   |   |           ___  |   |   |
(   62)       !!                  
(   63)       !!      case 1->   e3w(i+1) >= e3w(i) ( and e3w(j+1) >= e3w(j) ) then
(   64)       !!      tra~ = tra(i+1,j  ,k) + (e3w(i+1) - e3w(i)) * dk(TRAi+1)/e3w(i+1)
(   65)       !!    ( tra~ = tra(i  ,j+1,k) + (e3w(j+1) - e3w(j)) * dk(TRAj+1)/e3w(j+1))
(   66)       !!          or
(   67)       !!      case 2->   e3w(i+1) <= e3w(i) ( and e3w(j+1) <= e3w(j) ) then
(   68)       !!       tra~ = tra(i,j,k) + (e3w(i) - e3w(i+1)) * dk(TRAi)/e3w(i )
(   69)       !!     ( tra~ = tra(i,j,k) + (e3w(j) - e3w(j+1)) * dk(TRAj)/e3w(j ) )
(   70)       !!      
(   71)       !!
(   72)       !! ** Action  : - pgtru : horizontal gradient of TRA at U-points 
(   73)       !!              - pgtrv : horizontal gradient of TRA at V-points 
(   74)       !!
(   75)       !! History :
(   76)       !!   8.5  !  02-04  (A. Bozec)  Original code
(   77)       !!   8.5  !  02-08  (G. Madec E. Durand)  Optimization and Free form
(   78)       !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
(   79)       !!----------------------------------------------------------------------
(   80)       !! * Arguments
(   81)       INTEGER, INTENT( in ) ::   kt ! ocean time-step index
(   82)       REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT( in ) ::   &
(   83)          ptra                       ! 4D tracers fields
(   84)       REAL(wp), DIMENSION(jpi,jpj,jptra), INTENT( out ) ::   &
(   85)          pgtru,                 &  ! horizontal grad. of TRA u- and v-points 
(   86)          pgtrv                     ! of the partial step level
(   87) 
(   88)       !! * Local declarations
(   89)       INTEGER ::   ji, jj,jn,     &  ! Dummy loop indices
(   90)                    iku,ikv          ! partial step level at u- and v-points
(   91)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   92)          zti, ztj                   ! tempory arrays
(   93) 
(   94)       REAL(wp), DIMENSION(jpi,jpj,jptra) ::   &
(   95)          ztrai, ztraj               ! interpolated value of TRA
(   96) 
(   97)       REAL(wp) ::   &
(   98)          ze3wu, ze3wv,           &  ! temporary scalars
(   99)          zmaxu1, zmaxu2,         &  !    "         "
(  100)          zmaxv1, zmaxv2             !    "         "
(  101)       !!----------------------------------------------------------------------
(  102) 
(  103)       ! Initialization (first time-step only): compute mbatu and mbatv
(  104)       IF( kt == nittrc000 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 3

(  105)          mbatu(:,:) = 0
(  106)          mbatv(:,:) = 0
(  107)          DO jj = 1, jpjm1
(  108)             DO ji = 1, fs_jpim1   ! vector opt.
(  109)                mbatu(ji,jj) = MAX( MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1, 2 )
(  110)                mbatv(ji,jj) = MAX( MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1, 2 )
(  111)             END DO
(  112)          END DO
(  113)          zti(:,:) = FLOAT( mbatu(:,:) )
(  114)          ztj(:,:) = FLOAT( mbatv(:,:) )
(  115)          ! lateral boundary conditions: T-point, sign unchanged
(  116)          CALL lbc_lnk( zti , 'U', 1. )
(  117)          CALL lbc_lnk( ztj , 'V', 1. )
(  118)          mbatu(:,:) = MAX( INT( zti(:,:) ), 2 )
(  119)          mbatv(:,:) = MAX( INT( ztj(:,:) ), 2 )
(  120)       ENDIF
(  121)       
(  122) 
(  123)       DO jn = 1, jptra
(  124)          ! Interpolation of passive tracers at the last ocean level
(  125) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  126)          jj = 1
(  127)          DO ji = 1, jpij-jpi   ! vector opt. (forced unrolled)
(  128) # else
(  129)          DO jj = 1, jpjm1
(  130)             DO ji = 1, jpim1
(  131) # endif
(  132)                ! last level
(  133)                iku = mbatu(ji,jj)
(  134)                ikv = mbatv(ji,jj)
(  135) 
(  136)                ze3wu  = fse3w(ji+1,jj  ,iku) - fse3w(ji,jj,iku)
(  137)                ze3wv  = fse3w(ji  ,jj+1,ikv) - fse3w(ji,jj,ikv)
(  138)                zmaxu1 =  ze3wu / fse3w(ji+1,jj  ,iku)
(  139)                zmaxu2 = -ze3wu / fse3w(ji  ,jj  ,iku)
(  140)                zmaxv1 =  ze3wv / fse3w(ji  ,jj+1,ikv)
(  141)                zmaxv2 = -ze3wv / fse3w(ji  ,jj  ,ikv)
(  142) 
(  143)                ! i- direction
(  144) 
(  145)                IF( ze3wu >= 0. ) THEN      ! case 1
(  146)                   ! interpolated values of passive tracers
(  147)                   ztrai(ji,jj,jn) = ptra(ji+1,jj,iku,jn) + zmaxu1 * ( ptra(ji+1,jj,iku-1,jn) - ptra(ji+1,jj,iku,jn) )
(  148)                   ! gradient of passive tracers
(  149)                   pgtru(ji,jj,jn) = umask(ji,jj,1) * ( ztrai(ji,jj,jn) - ptra(ji,jj,iku,jn) )
(  150)                ELSE                        ! case 2
(  151)                   ! interpolated values of passive tracers
(  152)                   ztrai(ji,jj,jn) = ptra(ji,jj,iku,jn) + zmaxu2 * ( ptra(ji,jj,iku-1,jn) - ptra(ji,jj,iku,jn) )
(  153)                   ! gradient of passive tracers
(  154)                   pgtru(ji,jj,jn) = umask(ji,jj,1) * ( ptra(ji+1,jj,iku,jn) - ztrai (ji,jj,jn) )
(  155)                ENDIF
(  156) 
(  157)                ! j- direction
(  158) 
(  159)                IF( ze3wv >= 0. ) THEN      ! case 1
(  160)                   ! interpolated values of passive tracers
(  161)                   ztraj(ji,jj,jn) = ptra(ji,jj+1,ikv,jn) + zmaxv1 * ( ptra(ji,jj+1,ikv-1,jn) - ptra(ji,jj+1,ikv,jn) )
(  162)                   ! gradient of passive tracers






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 4

(  163)                   pgtrv(ji,jj,jn) = vmask(ji,jj,1) * ( ztraj(ji,jj,jn) - ptra(ji,jj,ikv,jn) )
(  164)                ELSE                        ! case 2
(  165)                   ! interpolated values of passive tracers
(  166)                   ztraj(ji,jj,jn) = ptra(ji,jj,ikv,jn) + zmaxv2 * ( ptra(ji,jj,ikv-1,jn) - ptra(ji,jj,ikv,jn) )
(  167)                   ! gradient of passive tracers
(  168)                   pgtrv(ji,jj,jn) = vmask(ji,jj,1) * ( ptra(ji,jj+1,ikv,jn) - ztraj(ji,jj,jn) )
(  169)                ENDIF
(  170) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  171)             END DO
(  172) # endif
(  173)          END DO
(  174) 
(  175)          ! Lateral boundary conditions on each gradient
(  176)          CALL lbc_lnk( pgtru(:,:,jn) , 'U', -1. ) 
(  177)          CALL lbc_lnk( pgtrv(:,:,jn) , 'V', -1. )
(  178) 
(  179)       END DO
(  180) 
(  181)    END SUBROUTINE zps_hde_trc
(  182) 
(  183) #else
(  184)    !!----------------------------------------------------------------------
(  185)    !!   Default option                                         Empty module
(  186)    !!----------------------------------------------------------------------
(  187)    USE par_kind
(  188) CONTAINS
(  189)    SUBROUTINE zps_hde_trc ( kt, ptra, pgtru, pgtrv ) ! Empty routine
(  190)       INTEGER, INTENT( in) :: kt
(  191)       REAL(wp), DIMENSION(:,:,:,:) :: ptra
(  192)       REAL(wp) :: pgtru, pgtrv
(  193) !      WRITE(*,*) 'zps_hde_trc: You should not have seen this print! error?',   &
(  194) !         kt, ptra, pgtru, pgtrv
(  195)    END SUBROUTINE zps_hde_trc
(  196) #endif
(  197) 
(  198)    !!======================================================================
(  199) END MODULE zpshde_trc



























PGF90 (Version     10.5)          09/23/2010  09:49:23      page 5

(    1) # 1 "zpshde_trc.F90"
(    1) MODULE zpshde_trc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE zpshde_trc   ***
(    4)    !! Ocean passive tracers: 
(    5)    !!==============================================================================
(    6) # 184
(  184)    !!----------------------------------------------------------------------
(  185)    !!   Default option                                         Empty module
(  186)    !!----------------------------------------------------------------------
(  187)    USE par_kind
(  188) CONTAINS
(  189)    SUBROUTINE zps_hde_trc ( kt, ptra, pgtru, pgtrv ) ! Empty routine
(  190)       INTEGER, INTENT( in) :: kt
(  191)       REAL(wp), DIMENSION(:,:,:,:) :: ptra
(  192)       REAL(wp) :: pgtru, pgtrv
(  193) !      WRITE(*,*) 'zps_hde_trc: You should not have seen this print! error?',   &
(  194) !         kt, ptra, pgtru, pgtrv
(  195)    END SUBROUTINE zps_hde_trc
(  196) 
(  197) # 198
(  198)    !!======================================================================
(  199) END MODULE zpshde_trc
