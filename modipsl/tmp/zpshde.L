


PGF90 (Version     10.5)          09/23/2010  09:50:21      page 1

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

Filename: zpshde.F90

(    1) MODULE zpshde
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE zpshde   ***
(    4)    !! Ocean active tracers: 
(    5)    !!==============================================================================
(    6) #if defined key_partial_steps || defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_partial_steps' :               z-coordinate with partial steps
(    9)    !!----------------------------------------------------------------------
(   10)    !!   zps_hde      :  Horizontal DErivative of T, S and rd at the last
(   11)    !!                   ocean level (Z-coord. with Partial Steps)
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE dom_oce         ! ocean space domain variables
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE phycst          ! physical constants
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE eosbn2          ! ocean equation of state
(   19)    USE lbclnk          ! lateral boundary conditions (or mpp link)
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC zps_hde          ! routine called by step.F90
(   26) 
(   27)    !! * module variables
(   28)    INTEGER, DIMENSION(jpi,jpj) ::   &
(   29)       mbatu, mbatv      ! bottom ocean level index at U- and V-points
(   30) 
(   31)    !! * Substitutions
(   32) #  include "domzgr_substitute.h90"
(   33) #  include "vectopt_loop_substitute.h90"
(   34)    !!----------------------------------------------------------------------
(   35)    !!----------------------------------------------------------------------
(   36)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/zpshde.F90,v 1.4 2005/03/27 18:35:20 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!----------------------------------------------------------------------
(   40) CONTAINS
(   41) 
(   42)    SUBROUTINE zps_hde ( kt, ptem, psal, prd ,   &
(   43)                             pgtu, pgsu, pgru,   &
(   44)                             pgtv, pgsv, pgrv  )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                     ***  ROUTINE zps_hde  ***






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 2

(   47)       !!                    
(   48)       !! ** Purpose :   Compute the horizontal derivative of T, S and rd
(   49)       !!      at u- and v-points with a linear interpolation for z-coordinate
(   50)       !!      with partial steps.
(   51)       !!
(   52)       !! ** Method  :   In z-coord with partial steps, scale factors on last 
(   53)       !!      levels are different for each grid point, so that T, S and rd 
(   54)       !!      points are not at the same depth as in z-coord. To have horizontal
(   55)       !!      gradients again, we interpolate T and S at the good depth : 
(   56)       !!      Linear interpolation of T, S   
(   57)       !!         Computation of di(tb) and dj(tb) by vertical interpolation:
(   58)       !!          di(t) = t~ - t(i,j,k) or t(i+1,j,k) - t~
(   59)       !!          dj(t) = t~ - t(i,j,k) or t(i,j+1,k) - t~
(   60)       !!         This formulation computes the two cases:
(   61)       !!                 CASE 1                   CASE 2  
(   62)       !!         k-1  ___ ___________   k-1   ___ ___________
(   63)       !!                    Ti  T~                  T~  Ti+1
(   64)       !!                  _____                        _____
(   65)       !!         k        |   |Ti+1     k           Ti |   |
(   66)       !!                  |   |____                ____|   |
(   67)       !!              ___ |   |   |           ___  |   |   |
(   68)       !!                  
(   69)       !!      case 1->   e3w(i+1) >= e3w(i) ( and e3w(j+1) >= e3w(j) ) then
(   70)       !!          t~ = t(i+1,j  ,k) + (e3w(i+1) - e3w(i)) * dk(Ti+1)/e3w(i+1)
(   71)       !!        ( t~ = t(i  ,j+1,k) + (e3w(j+1) - e3w(j)) * dk(Tj+1)/e3w(j+1)  )
(   72)       !!          or
(   73)       !!      case 2->   e3w(i+1) <= e3w(i) ( and e3w(j+1) <= e3w(j) ) then
(   74)       !!          t~ = t(i,j,k) + (e3w(i) - e3w(i+1)) * dk(Ti)/e3w(i )
(   75)       !!        ( t~ = t(i,j,k) + (e3w(j) - e3w(j+1)) * dk(Tj)/e3w(j ) )
(   76)       !!          Idem for di(s) and dj(s)          
(   77)       !!
(   78)       !!      For rho, we call eos_insitu_2d which will compute rd~(t~,s~) at 
(   79)       !!      the good depth zh from interpolated T and S for the different
(   80)       !!      formulation of the equation of state (eos).
(   81)       !!      Gradient formulation for rho :
(   82)       !!          di(rho) = rd~ - rd(i,j,k) or rd (i+1,j,k) - rd~
(   83)       !!
(   84)       !! ** Action  : - pgtu, pgsu, pgru: horizontal gradient of T, S
(   85)       !!                and rd at U-points 
(   86)       !!              - pgtv, pgsv, pgrv: horizontal gradient of T, S
(   87)       !!                and rd at V-points 
(   88)       !!
(   89)       !! History :
(   90)       !!   8.5  !  02-04  (A. Bozec)  Original code
(   91)       !!   8.5  !  02-08  (G. Madec E. Durand)  Optimization and Free form
(   92)       !!----------------------------------------------------------------------
(   93)       !! * Arguments
(   94)       INTEGER, INTENT( in ) ::   kt ! ocean time-step index
(   95)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(   96)          ptem, psal, prd            ! 3D T, S and rd fields
(   97)       REAL(wp), DIMENSION(jpi,jpj), INTENT( out ) ::   &
(   98)          pgtu, pgsu, pgru,       &  ! horizontal grad. of T, S and rd at u- 
(   99)          pgtv, pgsv, pgrv           ! and v-points of the partial step level
(  100) 
(  101)       !! * Local declarations
(  102)       INTEGER ::   ji, jj,       &  ! Dummy loop indices
(  103)                    iku,ikv          ! partial step level at u- and v-points
(  104)       REAL(wp), DIMENSION(jpi,jpj) ::   &






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 3

(  105)          zti, ztj, zsi, zsj,     &  ! interpolated value of T, S 
(  106)          zri, zrj,               &  ! and rd
(  107)          zhgi, zhgj                 ! depth of interpolation for eos2d
(  108)       REAL(wp) ::   &
(  109)          ze3wu, ze3wv,           &  ! temporary scalars
(  110)          zmaxu1, zmaxu2,         &  !    "         "
(  111)          zmaxv1, zmaxv2             !    "         "
(  112) 
(  113)       ! Initialization (first time-step only): compute mbatu and mbatv
(  114)       IF( kt == nit000 ) THEN
(  115)          mbatu(:,:) = 0
(  116)          mbatv(:,:) = 0
(  117)          DO jj = 1, jpjm1
(  118)             DO ji = 1, fs_jpim1   ! vector opt.
(  119)                mbatu(ji,jj) = MAX( MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1, 2 )
(  120)                mbatv(ji,jj) = MAX( MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1, 2 )
(  121)             END DO
(  122)          END DO
(  123)          zti(:,:) = FLOAT( mbatu(:,:) )
(  124)          ztj(:,:) = FLOAT( mbatv(:,:) )
(  125)          ! lateral boundary conditions: T-point, sign unchanged
(  126)          CALL lbc_lnk( zti , 'U', 1. )
(  127)          CALL lbc_lnk( ztj , 'V', 1. )
(  128)          mbatu(:,:) = MAX( INT( zti(:,:) ), 2 )
(  129)          mbatv(:,:) = MAX( INT( ztj(:,:) ), 2 )
(  130)       ENDIF
(  131)       
(  132) 
(  133)       ! Interpolation of T and S at the last ocean level
(  134) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  135)          jj = 1
(  136)          DO ji = 1, jpij-jpi   ! vector opt. (forced unrolled)
(  137) # else
(  138)       DO jj = 1, jpjm1
(  139)          DO ji = 1, jpim1
(  140) # endif
(  141)             ! last level
(  142)             iku = mbatu(ji,jj)
(  143)             ikv = mbatv(ji,jj)
(  144) 
(  145)             ze3wu  = fse3w(ji+1,jj  ,iku) - fse3w(ji,jj,iku)
(  146)             ze3wv  = fse3w(ji  ,jj+1,ikv) - fse3w(ji,jj,ikv)
(  147)             zmaxu1 =  ze3wu / fse3w(ji+1,jj  ,iku)
(  148)             zmaxu2 = -ze3wu / fse3w(ji  ,jj  ,iku)
(  149)             zmaxv1 =  ze3wv / fse3w(ji  ,jj+1,ikv)
(  150)             zmaxv2 = -ze3wv / fse3w(ji  ,jj  ,ikv)
(  151) 
(  152)             ! i- direction
(  153) 
(  154)             IF( ze3wu >= 0. ) THEN      ! case 1
(  155)                ! interpolated values of T and S
(  156)                zti(ji,jj) = ptem(ji+1,jj,iku) + zmaxu1 * ( ptem(ji+1,jj,iku-1) - ptem(ji+1,jj,iku) )
(  157)                zsi(ji,jj) = psal(ji+1,jj,iku) + zmaxu1 * ( psal(ji+1,jj,iku-1) - psal(ji+1,jj,iku) )
(  158)                ! depth of the partial step level
(  159)                zhgi(ji,jj) = fsdept(ji,jj,iku)
(  160)                ! gradient of T and S
(  161)                pgtu(ji,jj) = umask(ji,jj,1) * ( zti(ji,jj) - ptem(ji,jj,iku) )
(  162)                pgsu(ji,jj) = umask(ji,jj,1) * ( zsi(ji,jj) - psal(ji,jj,iku) )






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 4

(  163) 
(  164)             ELSE                        ! case 2
(  165)                ! interpolated values of T and S
(  166)                zti(ji,jj) = ptem(ji,jj,iku) + zmaxu2 * ( ptem(ji,jj,iku-1) - ptem(ji,jj,iku) )
(  167)                zsi(ji,jj) = psal(ji,jj,iku) + zmaxu2 * ( psal(ji,jj,iku-1) - psal(ji,jj,iku) )
(  168)                ! depth of the partial step level
(  169)                zhgi(ji,jj) = fsdept(ji+1,jj,iku)
(  170)                ! gradient of T and S 
(  171)                pgtu(ji,jj) = umask(ji,jj,1) * ( ptem(ji+1,jj,iku) - zti (ji,jj) )
(  172)                pgsu(ji,jj) = umask(ji,jj,1) * ( psal(ji+1,jj,iku) - zsi (ji,jj) )
(  173)             ENDIF
(  174) 
(  175)             ! j- direction
(  176) 
(  177)             IF( ze3wv >= 0. ) THEN      ! case 1
(  178)                ! interpolated values of T and S
(  179)                ztj(ji,jj) = ptem(ji,jj+1,ikv) + zmaxv1 * ( ptem(ji,jj+1,ikv-1) - ptem(ji,jj+1,ikv) )
(  180)                zsj(ji,jj) = psal(ji,jj+1,ikv) + zmaxv1 * ( psal(ji,jj+1,ikv-1) - psal(ji,jj+1,ikv) )
(  181)                ! depth of the partial step level
(  182)                zhgj(ji,jj) = fsdept(ji,jj,ikv) 
(  183)                ! gradient of T and S
(  184)                pgtv(ji,jj) = vmask(ji,jj,1) * ( ztj(ji,jj) - ptem(ji,jj,ikv) )
(  185)                pgsv(ji,jj) = vmask(ji,jj,1) * ( zsj(ji,jj) - psal(ji,jj,ikv) )
(  186) 
(  187)             ELSE                        ! case 2
(  188)                ! interpolated values of T and S
(  189)                ztj(ji,jj) = ptem(ji,jj,ikv) + zmaxv2 * ( ptem(ji,jj,ikv-1) - ptem(ji,jj,ikv) )
(  190)                zsj(ji,jj) = psal(ji,jj,ikv) + zmaxv2 * ( psal(ji,jj,ikv-1) - psal(ji,jj,ikv) ) 
(  191)                ! depth of the partial step level
(  192)                zhgj(ji,jj) = fsdept(ji,jj+1,ikv) 
(  193)                ! gradient of T and S
(  194)                pgtv(ji,jj) = vmask(ji,jj,1) * ( ptem(ji,jj+1,ikv) - ztj(ji,jj) )
(  195)                pgsv(ji,jj) = vmask(ji,jj,1) * ( psal(ji,jj+1,ikv) - zsj(ji,jj) )
(  196)             ENDIF
(  197) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  198)          END DO
(  199) # endif
(  200)       END DO
(  201) 
(  202)       ! Compute interpolated rd from zti, zsi, ztj, zsj for the 2 cases at the depth of the partial
(  203)       ! step and store it in  zri, zrj for each  case
(  204)       CALL eos( zti, zsi, zhgi, zri )
(  205)       CALL eos( ztj, zsj, zhgj, zrj )
(  206) 
(  207) 
(  208)       ! Gradient of density at the last level 
(  209) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  210)          jj = 1
(  211)          DO ji = 1, jpij-jpi   ! vector opt. (forced unrolled)
(  212) # else
(  213)       DO jj = 1, jpjm1
(  214)          DO ji = 1, jpim1
(  215) # endif
(  216)             iku = mbatu(ji,jj)
(  217)             ikv = mbatv(ji,jj)
(  218)             ze3wu  = fse3w(ji+1,jj  ,iku) - fse3w(ji,jj,iku)
(  219)             ze3wv  = fse3w(ji  ,jj+1,ikv) - fse3w(ji,jj,ikv)
(  220)             IF( ze3wu >= 0. ) THEN    ! i-direction: case 1






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 5

(  221)                pgru(ji,jj) = umask(ji,jj,1) * ( zri(ji,jj) - prd(ji,jj,iku) )
(  222)             ELSE                      ! i-direction: case 2
(  223)                pgru(ji,jj) = umask(ji,jj,1) * ( prd(ji+1,jj,iku) - zri(ji,jj) )
(  224)             ENDIF
(  225)             IF( ze3wv >= 0. ) THEN    ! j-direction: case 1
(  226)                pgrv(ji,jj) = vmask(ji,jj,1) * ( zrj(ji,jj) - prd(ji,jj,ikv) )  
(  227)             ELSE                      ! j-direction: case 2
(  228)                pgrv(ji,jj) = vmask(ji,jj,1) * ( prd(ji,jj+1,ikv) - zrj(ji,jj) )
(  229)             ENDIF
(  230) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  231)          END DO
(  232) # endif
(  233)       END DO
(  234) 
(  235)       ! Lateral boundary conditions on each gradient
(  236)       CALL lbc_lnk( pgtu , 'U', -1. )   ;   CALL lbc_lnk( pgtv , 'V', -1. )
(  237)       CALL lbc_lnk( pgsu , 'U', -1. )   ;   CALL lbc_lnk( pgsv , 'V', -1. )
(  238)       CALL lbc_lnk( pgru , 'U', -1. )   ;   CALL lbc_lnk( pgrv , 'V', -1. )
(  239) 
(  240)    END SUBROUTINE zps_hde
(  241) 
(  242) #else
(  243)    !!----------------------------------------------------------------------
(  244)    !!   Default option                                         Empty module
(  245)    !!----------------------------------------------------------------------
(  246)    USE par_kind
(  247) CONTAINS
(  248)    SUBROUTINE zps_hde ( kt, ptem, psal, prd ,   &      ! Empty routine
(  249)                             pgtu, pgsu, pgru,   &
(  250)                             pgtv, pgsv, pgrv  )
(  251)       REAL(wp), DIMENSION(:,:,:) :: ptem, psal, prd
(  252)       REAL(wp) :: pgtu, pgsu, pgru, pgtv, pgsv, pgrv
(  253) !      WRITE(*,*) 'zps_hde: You should not have seen this print! error?',   &
(  254) !         kt, ptem, psal, prd, pgtu, pgsu, pgru, pgtv, pgsv, pgrv
(  255)    END SUBROUTINE zps_hde
(  256) #endif
(  257) 
(  258)    !!======================================================================
(  259) END MODULE zpshde

























PGF90 (Version     10.5)          09/23/2010  09:50:21      page 6

(    1) # 1 "zpshde.F90"
(    1) MODULE zpshde
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE zpshde   ***
(    4)    !! Ocean active tracers: 
(    5)    !!==============================================================================
(    6) # 243
(  243)    !!----------------------------------------------------------------------
(  244)    !!   Default option                                         Empty module
(  245)    !!----------------------------------------------------------------------
(  246)    USE par_kind
(  247) CONTAINS
(  248)    SUBROUTINE zps_hde ( kt, ptem, psal, prd ,   &      ! Empty routine
(  249)                             pgtu, pgsu, pgru,   &
(  250)                             pgtv, pgsv, pgrv  )
(  251)       REAL(wp), DIMENSION(:,:,:) :: ptem, psal, prd
(  252)       REAL(wp) :: pgtu, pgsu, pgru, pgtv, pgsv, pgrv
(  253) !      WRITE(*,*) 'zps_hde: You should not have seen this print! error?',   &
(  254) !         kt, ptem, psal, prd, pgtu, pgsu, pgru, pgtv, pgsv, pgrv
(  255)    END SUBROUTINE zps_hde
(  256) 
(  257) # 258
(  258)    !!======================================================================
(  259) END MODULE zpshde
