


PGF90 (Version     10.5)          09/23/2010  09:49:36      page 1

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

Filename: dynldf_bilap.F90

(    1) MODULE dynldf_bilap
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  dynldf_bilap  ***
(    4)    !! Ocean dynamics:  lateral viscosity trend
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_ldf_bilap : update the momentum trend with the lateral diffusion
(    9)    !!                   using an iso-level bilaplacian operator
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE ldfdyn_oce      ! ocean dynamics: lateral physics
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE trdmod          ! ocean dynamics trends 
(   17)    USE trdmod_oce      ! ocean variables trends
(   18)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC dyn_ldf_bilap  ! called by step.F90
(   26) 
(   27)    !! * Substitutions
(   28) #  include "domzgr_substitute.h90"
(   29) #  include "ldfdyn_substitute.h90"
(   30) #  include "vectopt_loop_substitute.h90"
(   31)    !!----------------------------------------------------------------------
(   32)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_bilap.F90,v 1.8 2005/09/02 15:45:23 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36) 
(   37) CONTAINS
(   38) 
(   39)    SUBROUTINE dyn_ldf_bilap( kt )
(   40)       !!----------------------------------------------------------------------
(   41)       !!                     ***  ROUTINE dyn_ldf_bilap  ***
(   42)       !!
(   43)       !! ** Purpose :   Compute the before trend of the lateral momentum
(   44)       !!      diffusion and add it to the general trend of momentum equation.
(   45)       !!
(   46)       !! ** Method  :   The before horizontal momentum diffusion trend is a 






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 2

(   47)       !!      bi-harmonic operator (bilaplacian type) which separates the
(   48)       !!      divergent and rotational parts of the flow.
(   49)       !!      Its horizontal components are computed as follow:
(   50)       !!      laplacian:
(   51)       !!          zlu = 1/e1u di[ hdivb ] - 1/(e2u*e3u) dj-1[ e3f rotb ]
(   52)       !!          zlv = 1/e2v dj[ hdivb ] + 1/(e1v*e3v) di-1[ e3f rotb ]
(   53)       !!      third derivative:
(   54)       !!       * multiply by the eddy viscosity coef. at u-, v-point, resp.
(   55)       !!          zlu = ahmu * zlu
(   56)       !!          zlv = ahmv * zlv
(   57)       !!       * curl and divergence of the laplacian
(   58)       !!          zuf = 1/(e1f*e2f) ( di[e2v zlv] - dj[e1u zlu] )
(   59)       !!          zut = 1/(e1t*e2t*e3t) ( di[e2u*e3u zlu] + dj[e1v*e3v zlv] )
(   60)       !!      bilaplacian:
(   61)       !!              diffu = 1/e1u di[ zut ] - 1/(e2u*e3u) dj-1[ e3f zuf ]
(   62)       !!              diffv = 1/e2v dj[ zut ] + 1/(e1v*e3v) di-1[ e3f zuf ]
(   63)       !!      If lk_sco=F and lk_zps=F, the vertical scale factors in the
(   64)       !!      rotational part of the diffusion are simplified
(   65)       !!      Add this before trend to the general trend (ua,va):
(   66)       !!            (ua,va) = (ua,va) + (diffu,diffv)
(   67)       !!      'key_trddyn' defined: the two components of the horizontal
(   68)       !!                               diffusion trend are saved.
(   69)       !!
(   70)       !! ** Action : - Update (ua,va) with the before iso-level biharmonic
(   71)       !!               mixing trend.
(   72)       !!             - Save in (ztdua,ztdva) the trends ('key_trddyn')
(   73)       !!
(   74)       !! History :
(   75)       !!        !  90-09  (G. Madec)  Original code
(   76)       !!        !  91-11  (G. Madec)
(   77)       !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
(   78)       !!        !  96-01  (G. Madec)  statement function for e3
(   79)       !!        !  97-07  (G. Madec)  lbc calls
(   80)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   81)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   82)       !!----------------------------------------------------------------------
(   83)       !! * Modules used     
(   84)       USE oce, ONLY :    ztdua => ta,      & ! use ta as 3D workspace   
(   85)                          ztdva => sa         ! use sa as 3D workspace   
(   86) 
(   87)       !! * Arguments
(   88)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   89) 
(   90)       !! * Local declarations
(   91)       INTEGER  ::   ji, jj, jk                ! dummy loop indices
(   92)       REAL(wp) ::   zua, zva, zbt, ze2u, ze2v ! temporary scalar
(   93)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   94)          zuf, zut, zlu, zlv, zcu, zcv         ! temporary workspace
(   95)       !!----------------------------------------------------------------------
(   96)       !!  OPA 8.5, LODYC-IPSL (2002)
(   97)       !!----------------------------------------------------------------------
(   98) 
(   99)       IF( kt == nit000 ) THEN
(  100)          IF(lwp) WRITE(numout,*)
(  101)          IF(lwp) WRITE(numout,*) 'dyn_ldf_bilap : iso-level bilaplacian operator'
(  102)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
(  103)       ENDIF
(  104)       zuf(:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 3

(  105)       zut(:,:) = 0.e0
(  106)       zlu(:,:) = 0.e0
(  107)       zlv(:,:) = 0.e0
(  108) 
(  109)       ! Save ua and va trends
(  110)       IF( l_trddyn )   THEN
(  111)          ztdua(:,:,:) = ua(:,:,:) 
(  112)          ztdva(:,:,:) = va(:,:,:) 
(  113)       ENDIF
(  114)       !                                                ! ===============
(  115)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  116)          !                                             ! ===============
(  117)          ! Laplacian
(  118)          ! ---------
(  119) 
(  120)          IF( lk_sco .OR. lk_zps ) THEN   ! s-coordinate or z-coordinate with partial steps
(  121)             zuf(:,:) = rotb(:,:,jk) * fse3f(:,:,jk)
(  122)             DO jj = 2, jpjm1
(  123)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  124)                   zlu(ji,jj) = - ( zuf(ji,jj) - zuf(ji,jj-1) ) / ( e2u(ji,jj) * fse3u(ji,jj,jk) )   &
(  125)                      &         + ( hdivb(ji+1,jj,jk) - hdivb(ji,jj,jk) ) / e1u(ji,jj)
(  126)    
(  127)                   zlv(ji,jj) = + ( zuf(ji,jj) - zuf(ji-1,jj) ) / ( e1v(ji,jj) * fse3v(ji,jj,jk) )   &
(  128)                      &         + ( hdivb(ji,jj+1,jk) - hdivb(ji,jj,jk) ) / e2v(ji,jj)
(  129)                END DO
(  130)             END DO
(  131)          ELSE                            ! z-coordinate
(  132)             DO jj = 2, jpjm1
(  133)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  134)                   zlu(ji,jj) = - ( rotb (ji  ,jj,jk) - rotb (ji,jj-1,jk) ) / e2u(ji,jj)   &
(  135)                      &         + ( hdivb(ji+1,jj,jk) - hdivb(ji,jj  ,jk) ) / e1u(ji,jj)
(  136)    
(  137)                   zlv(ji,jj) = + ( rotb (ji,jj  ,jk) - rotb (ji-1,jj,jk) ) / e1v(ji,jj)   &
(  138)                      &         + ( hdivb(ji,jj+1,jk) - hdivb(ji  ,jj,jk) ) / e2v(ji,jj)
(  139)                END DO  
(  140)             END DO  
(  141)          ENDIF
(  142) 
(  143)          ! Boundary conditions on the laplacian  (zlu,zlv)
(  144)          CALL lbc_lnk( zlu, 'U', -1. )
(  145)          CALL lbc_lnk( zlv, 'V', -1. )
(  146)          
(  147)          
(  148)          ! Third derivative
(  149)          ! ----------------
(  150)          
(  151)          ! Multiply by the eddy viscosity coef. (at u- and v-points)
(  152)          zlu(:,:) = zlu(:,:) * fsahmu(:,:,jk)
(  153)          zlv(:,:) = zlv(:,:) * fsahmv(:,:,jk)
(  154)          
(  155)          ! Contravariant "laplacian"
(  156)          zcu(:,:) = e1u(:,:) * zlu(:,:)
(  157)          zcv(:,:) = e2v(:,:) * zlv(:,:)
(  158)          
(  159)          ! Laplacian curl ( * e3f if s-coordinates or z-coordinate with partial steps)
(  160)          DO jj = 1, jpjm1
(  161)             DO ji = 1, fs_jpim1   ! vector opt.
(  162)                zuf(ji,jj) = fmask(ji,jj,jk) * (  zcv(ji+1,jj  ) - zcv(ji,jj)      &






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 4

(  163)                   &                            - zcu(ji  ,jj+1) + zcu(ji,jj)  )   &
(  164) #if defined key_s_coord || defined key_partial_steps
(  165)                   &       * fse3f(ji,jj,jk) / ( e1f(ji,jj)*e2f(ji,jj) )
(  166) #else
(  167)                   &                         / ( e1f(ji,jj)*e2f(ji,jj) )
(  168) #endif
(  169)             END DO  
(  170)          END DO  
(  171) 
(  172)          ! Laplacian Horizontal fluxes
(  173)          DO jj = 1, jpjm1
(  174)             DO ji = 1, fs_jpim1   ! vector opt.
(  175) #if defined key_s_coord || defined key_partial_steps
(  176)                zlu(ji,jj) = e2u(ji,jj) * fse3u(ji,jj,jk) * zlu(ji,jj)
(  177)                zlv(ji,jj) = e1v(ji,jj) * fse3v(ji,jj,jk) * zlv(ji,jj)
(  178) #else
(  179)                zlu(ji,jj) = e2u(ji,jj) * zlu(ji,jj)
(  180)                zlv(ji,jj) = e1v(ji,jj) * zlv(ji,jj)
(  181) #endif
(  182)             END DO
(  183)          END DO
(  184) 
(  185)          ! Laplacian divergence
(  186)          DO jj = 2, jpj
(  187)             DO ji = fs_2, jpi   ! vector opt.
(  188) #if defined key_s_coord || defined key_partial_steps
(  189)                zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
(  190) #else
(  191)                zbt = e1t(ji,jj) * e2t(ji,jj)
(  192) #endif
(  193)                zut(ji,jj) = (  zlu(ji,jj) - zlu(ji-1,jj  )   &
(  194)                   &          + zlv(ji,jj) - zlv(ji  ,jj-1)  ) / zbt
(  195)             END DO
(  196)          END DO
(  197) 
(  198) 
(  199)       ! boundary conditions on the laplacian curl and div (zuf,zut)
(  200)       CALL lbc_lnk( zuf, 'F', 1. )
(  201)       CALL lbc_lnk( zut, 'T', 1. )
(  202) 
(  203)          
(  204)          ! Bilaplacian
(  205)          ! -----------
(  206) 
(  207)          DO jj = 2, jpjm1
(  208)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  209) #if defined key_s_coord || defined key_partial_steps
(  210)                ze2u = e2u(ji,jj) * fse3u(ji,jj,jk)
(  211)                ze2v = e1v(ji,jj) * fse3v(ji,jj,jk)
(  212) #else
(  213)                ze2u = e2u(ji,jj)
(  214)                ze2v = e1v(ji,jj)
(  215) #endif
(  216)                ! horizontal biharmonic diffusive trends
(  217)                zua = - ( zuf(ji  ,jj) - zuf(ji,jj-1) ) / ze2u   &
(  218)                   &  + ( zut(ji+1,jj) - zut(ji,jj  ) ) / e1u(ji,jj)
(  219) 
(  220)                zva = + ( zuf(ji,jj  ) - zuf(ji-1,jj) ) / ze2v   &






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 5

(  221)                   &  + ( zut(ji,jj+1) - zut(ji  ,jj) ) / e2v(ji,jj)
(  222)                ! add it to the general momentum trends
(  223)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  224)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  225)             END DO
(  226)          END DO
(  227) 
(  228)          !                                             ! ===============
(  229)       END DO                                           !   End of slab
(  230)       !                                                ! ===============
(  231)       ! save the lateral diffusion trends for diagnostic
(  232)       ! momentum trends
(  233)       IF( l_trddyn )   THEN
(  234)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  235)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  236) 
(  237)          CALL trd_mod(ztdua, ztdva, jpdtdldf, 'DYN', kt)
(  238)       ENDIF
(  239) 
(  240)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  241)          CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
(  242)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  243)       ENDIF
(  244) 
(  245)    END SUBROUTINE dyn_ldf_bilap
(  246) 
(  247)    !!======================================================================
(  248) END MODULE dynldf_bilap




































PGF90 (Version     10.5)          09/23/2010  09:49:36      page 6

(    1) # 1 "dynldf_bilap.F90"
(    1) MODULE dynldf_bilap
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  dynldf_bilap  ***
(    4)    !! Ocean dynamics:  lateral viscosity trend
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_ldf_bilap : update the momentum trend with the lateral diffusion
(    9)    !!                   using an iso-level bilaplacian operator
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE ldfdyn_oce      ! ocean dynamics: lateral physics
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE trdmod          ! ocean dynamics trends 
(   17)    USE trdmod_oce      ! ocean variables trends
(   18)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC dyn_ldf_bilap  ! called by step.F90
(   26) 
(   27)    !! * Substitutions
(   28) # 1 "./domzgr_substitute.h90"
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
(   48) # 29 "dynldf_bilap.F90"
(   29) # 1 "./ldfdyn_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  ldfdyn_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsahm., the lateral eddy viscosity coeff. 
(    5)    !!      with a constant, or 1D, or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn_substitute.h90,v 1.3 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19
(   19) !   ' key_dynldf_c3d' :                  3D coefficient






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 7

(   20) # 30 "dynldf_bilap.F90"
(   30) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
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
(   12) # 31 "dynldf_bilap.F90"
(   31) # 31
(   31)    !!----------------------------------------------------------------------
(   32)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_bilap.F90,v 1.8 2005/09/02 15:45:23 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36) 
(   37) CONTAINS
(   38) 
(   39)    SUBROUTINE dyn_ldf_bilap( kt )
(   40)       !!----------------------------------------------------------------------
(   41)       !!                     ***  ROUTINE dyn_ldf_bilap  ***
(   42)       !!
(   43)       !! ** Purpose :   Compute the before trend of the lateral momentum
(   44)       !!      diffusion and add it to the general trend of momentum equation.
(   45)       !!
(   46)       !! ** Method  :   The before horizontal momentum diffusion trend is a 
(   47)       !!      bi-harmonic operator (bilaplacian type) which separates the
(   48)       !!      divergent and rotational parts of the flow.
(   49)       !!      Its horizontal components are computed as follow:
(   50)       !!      laplacian:
(   51)       !!          zlu = 1/e1u di[ hdivb ] - 1/(e2u*e3u) dj-1[ e3f rotb ]
(   52)       !!          zlv = 1/e2v dj[ hdivb ] + 1/(e1v*e3v) di-1[ e3f rotb ]
(   53)       !!      third derivative:
(   54)       !!       * multiply by the eddy viscosity coef. at u-, v-point, resp.
(   55)       !!          zlu = ahmu * zlu
(   56)       !!          zlv = ahmv * zlv
(   57)       !!       * curl and divergence of the laplacian
(   58)       !!          zuf = 1/(e1f*e2f) ( di[e2v zlv] - dj[e1u zlu] )
(   59)       !!          zut = 1/(e1t*e2t*e3t) ( di[e2u*e3u zlu] + dj[e1v*e3v zlv] )
(   60)       !!      bilaplacian:
(   61)       !!              diffu = 1/e1u di[ zut ] - 1/(e2u*e3u) dj-1[ e3f zuf ]
(   62)       !!              diffv = 1/e2v dj[ zut ] + 1/(e1v*e3v) di-1[ e3f zuf ]
(   63)       !!      If lk_sco=F and lk_zps=F, the vertical scale factors in the
(   64)       !!      rotational part of the diffusion are simplified
(   65)       !!      Add this before trend to the general trend (ua,va):
(   66)       !!            (ua,va) = (ua,va) + (diffu,diffv)
(   67)       !!      'key_trddyn' defined: the two components of the horizontal
(   68)       !!                               diffusion trend are saved.
(   69)       !!
(   70)       !! ** Action : - Update (ua,va) with the before iso-level biharmonic
(   71)       !!               mixing trend.
(   72)       !!             - Save in (ztdua,ztdva) the trends ('key_trddyn')
(   73)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 8

(   74)       !! History :
(   75)       !!        !  90-09  (G. Madec)  Original code
(   76)       !!        !  91-11  (G. Madec)
(   77)       !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
(   78)       !!        !  96-01  (G. Madec)  statement function for e3
(   79)       !!        !  97-07  (G. Madec)  lbc calls
(   80)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   81)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   82)       !!----------------------------------------------------------------------
(   83)       !! * Modules used     
(   84)       USE oce, ONLY :    ztdua => ta,      & ! use ta as 3D workspace   
(   85)                          ztdva => sa         ! use sa as 3D workspace   
(   86) 
(   87)       !! * Arguments
(   88)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   89) 
(   90)       !! * Local declarations
(   91)       INTEGER  ::   ji, jj, jk                ! dummy loop indices
(   92)       REAL(wp) ::   zua, zva, zbt, ze2u, ze2v ! temporary scalar
(   93)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   94)          zuf, zut, zlu, zlv, zcu, zcv         ! temporary workspace
(   95)       !!----------------------------------------------------------------------
(   96)       !!  OPA 8.5, LODYC-IPSL (2002)
(   97)       !!----------------------------------------------------------------------
(   98) 
(   99)       IF( kt == nit000 ) THEN
(  100)          IF(lwp) WRITE(numout,*)
(  101)          IF(lwp) WRITE(numout,*) 'dyn_ldf_bilap : iso-level bilaplacian operator'
(  102)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
(  103)       ENDIF
(  104)       zuf(:,:) = 0.e0
(  105)       zut(:,:) = 0.e0
(  106)       zlu(:,:) = 0.e0
(  107)       zlv(:,:) = 0.e0
(  108) 
(  109)       ! Save ua and va trends
(  110)       IF( l_trddyn )   THEN
(  111)          ztdua(:,:,:) = ua(:,:,:) 
(  112)          ztdva(:,:,:) = va(:,:,:) 
(  113)       ENDIF
(  114)       !                                                ! ===============
(  115)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  116)          !                                             ! ===============
(  117)          ! Laplacian
(  118)          ! ---------
(  119) 
(  120)          IF( lk_sco .OR. lk_zps ) THEN   ! s-coordinate or z-coordinate with partial steps
(  121)             zuf(:,:) = rotb(:,:,jk) * e3t(jk)
(  122)             DO jj = 2, jpjm1
(  123)                DO ji = 2, jpim1   ! vector opt.
(  124)                   zlu(ji,jj) = - ( zuf(ji,jj) - zuf(ji,jj-1) ) / ( e2u(ji,jj) * e3t(jk) )   &
(  125)                      &         + ( hdivb(ji+1,jj,jk) - hdivb(ji,jj,jk) ) / e1u(ji,jj)
(  126)    
(  127)                   zlv(ji,jj) = + ( zuf(ji,jj) - zuf(ji-1,jj) ) / ( e1v(ji,jj) * e3t(jk) )   &
(  128)                      &         + ( hdivb(ji,jj+1,jk) - hdivb(ji,jj,jk) ) / e2v(ji,jj)
(  129)                END DO
(  130)             END DO
(  131)          ELSE                            ! z-coordinate






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 9

(  132)             DO jj = 2, jpjm1
(  133)                DO ji = 2, jpim1   ! vector opt.
(  134)                   zlu(ji,jj) = - ( rotb (ji  ,jj,jk) - rotb (ji,jj-1,jk) ) / e2u(ji,jj)   &
(  135)                      &         + ( hdivb(ji+1,jj,jk) - hdivb(ji,jj  ,jk) ) / e1u(ji,jj)
(  136)    
(  137)                   zlv(ji,jj) = + ( rotb (ji,jj  ,jk) - rotb (ji-1,jj,jk) ) / e1v(ji,jj)   &
(  138)                      &         + ( hdivb(ji,jj+1,jk) - hdivb(ji  ,jj,jk) ) / e2v(ji,jj)
(  139)                END DO  
(  140)             END DO  
(  141)          ENDIF
(  142) 
(  143)          ! Boundary conditions on the laplacian  (zlu,zlv)
(  144)          CALL lbc_lnk( zlu, 'U', -1. )
(  145)          CALL lbc_lnk( zlv, 'V', -1. )
(  146)          
(  147)          
(  148)          ! Third derivative
(  149)          ! ----------------
(  150)          
(  151)          ! Multiply by the eddy viscosity coef. (at u- and v-points)
(  152)          zlu(:,:) = zlu(:,:) * ahm3(:,:,jk)
(  153)          zlv(:,:) = zlv(:,:) * ahm4(:,:,jk)
(  154)          
(  155)          ! Contravariant "laplacian"
(  156)          zcu(:,:) = e1u(:,:) * zlu(:,:)
(  157)          zcv(:,:) = e2v(:,:) * zlv(:,:)
(  158)          
(  159)          ! Laplacian curl ( * e3f if s-coordinates or z-coordinate with partial steps)
(  160)          DO jj = 1, jpjm1
(  161)             DO ji = 1, jpim1   ! vector opt.
(  162)                zuf(ji,jj) = fmask(ji,jj,jk) * (  zcv(ji+1,jj  ) - zcv(ji,jj)      &
(  163)                   &                            - zcu(ji  ,jj+1) + zcu(ji,jj)  )   &
(  167)                   &                         / ( e1f(ji,jj)*e2f(ji,jj) )
(  168) # 169
(  169)             END DO  
(  170)          END DO  
(  171) 
(  172)          ! Laplacian Horizontal fluxes
(  173)          DO jj = 1, jpjm1
(  174)             DO ji = 1, jpim1   ! vector opt.
(  175) # 179
(  179)                zlu(ji,jj) = e2u(ji,jj) * zlu(ji,jj)
(  180)                zlv(ji,jj) = e1v(ji,jj) * zlv(ji,jj)
(  181) # 182
(  182)             END DO
(  183)          END DO
(  184) 
(  185)          ! Laplacian divergence
(  186)          DO jj = 2, jpj
(  187)             DO ji = 2, jpi   ! vector opt.
(  188) # 191
(  191)                zbt = e1t(ji,jj) * e2t(ji,jj)
(  192) # 193
(  193)                zut(ji,jj) = (  zlu(ji,jj) - zlu(ji-1,jj  )   &
(  194)                   &          + zlv(ji,jj) - zlv(ji  ,jj-1)  ) / zbt
(  195)             END DO
(  196)          END DO
(  197) 






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 10

(  198) 
(  199)       ! boundary conditions on the laplacian curl and div (zuf,zut)
(  200)       CALL lbc_lnk( zuf, 'F', 1. )
(  201)       CALL lbc_lnk( zut, 'T', 1. )
(  202) 
(  203)          
(  204)          ! Bilaplacian
(  205)          ! -----------
(  206) 
(  207)          DO jj = 2, jpjm1
(  208)             DO ji = 2, jpim1   ! vector opt.
(  209) # 213
(  213)                ze2u = e2u(ji,jj)
(  214)                ze2v = e1v(ji,jj)
(  215) # 216
(  216)                ! horizontal biharmonic diffusive trends
(  217)                zua = - ( zuf(ji  ,jj) - zuf(ji,jj-1) ) / ze2u   &
(  218)                   &  + ( zut(ji+1,jj) - zut(ji,jj  ) ) / e1u(ji,jj)
(  219) 
(  220)                zva = + ( zuf(ji,jj  ) - zuf(ji-1,jj) ) / ze2v   &
(  221)                   &  + ( zut(ji,jj+1) - zut(ji  ,jj) ) / e2v(ji,jj)
(  222)                ! add it to the general momentum trends
(  223)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  224)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  225)             END DO
(  226)          END DO
(  227) 
(  228)          !                                             ! ===============
(  229)       END DO                                           !   End of slab
(  230)       !                                                ! ===============
(  231)       ! save the lateral diffusion trends for diagnostic
(  232)       ! momentum trends
(  233)       IF( l_trddyn )   THEN
(  234)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  235)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  236) 
(  237)          CALL trd_mod(ztdua, ztdva, jpdtdldf, 'DYN', kt)
(  238)       ENDIF
(  239) 
(  240)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  241)          CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
(  242)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  243)       ENDIF
(  244) 
(  245)    END SUBROUTINE dyn_ldf_bilap
(  246) 
(  247)    !!======================================================================
(  248) END MODULE dynldf_bilap
