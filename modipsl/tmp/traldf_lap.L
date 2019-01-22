


PGF90 (Version     10.5)          09/23/2010  09:50:15      page 1

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

Filename: traldf_lap.F90

(    1) MODULE traldf_lap
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traldf_lap  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_ldf_lap  : update the tracer trend with the horizontal diffusion
(    9)    !!                 using a iso-level harmonic (laplacien) operator.
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   15)    USE trdmod          ! ocean active tracers trends 
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE diaptr          ! poleward transport diagnostics
(   19)    USE prtctl          ! Print control
(   20) 
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC tra_ldf_lap  ! routine called by step.F90
(   27) 
(   28)    !! * Substitutions
(   29) #  include "domzgr_substitute.h90"
(   30) #  include "ldftra_substitute.h90"
(   31) #  include "vectopt_loop_substitute.h90"
(   32)    !!----------------------------------------------------------------------
(   33)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   34)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_lap.F90,v 1.7 2005/09/02 15:45:33 opalod Exp $ 
(   35)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   36)    !!----------------------------------------------------------------------
(   37)    
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE tra_ldf_lap( kt )
(   41)       !!----------------------------------------------------------------------
(   42)       !!                  ***  ROUTINE tra_ldf_lap  ***
(   43)       !!                   
(   44)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   45)       !!      trend and add it to the general trend of tracer equation.
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 2

(   47)       !! ** Method  :   Second order diffusive operator evaluated using before
(   48)       !!      fields (forward time scheme). The horizontal diffusive trends of 
(   49)       !!      temperature (idem for salinity) is given by:
(   50)       !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
(   51)       !!      factors e3. are inside the derivatives:
(   52)       !!          difft = 1/(e1t*e2t*e3t) {  di-1[ aht e2u*e3u/e1u di(tb) ]
(   53)       !!                                   + dj-1[ aht e1v*e3v/e2v dj(tb) ] }
(   54)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(   55)       !!          difft = 1/(e1t*e2t) {  di-1[ aht e2u/e1u di(tb) ]
(   56)       !!                               + dj-1[ aht e1v/e2v dj(tb) ] }
(   57)       !!      Add this trend to the general tracer trend (ta,sa):
(   58)       !!          (ta,sa) = (ta,sa) + ( difft , diffs )
(   59)       !!
(   60)       !! ** Action  : - Update (ta,sa) arrays with the before iso-level 
(   61)       !!                harmonic mixing trend.
(   62)       !!              - Save the trends in (ztdta,ztdsa) ('key_trdtra')
(   63)       !!
(   64)       !! History :
(   65)       !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
(   66)       !!        !  91-11  (G. Madec)
(   67)       !!        !  95-11  (G. Madec)  suppress volumetric scale factors
(   68)       !!        !  96-01  (G. Madec)  statement function for e3
(   69)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   70)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   71)       !!----------------------------------------------------------------------
(   72)       USE oce              , ztu => ua,  &  ! use ua as workspace
(   73)          &                   zsu => va      ! use va as workspace
(   74) 
(   75)       !! * Arguments
(   76)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   77)       
(   78)       !! * Local save
(   79)       REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
(   80)          ze1ur, ze2vr, zbtr2              ! scale factor coefficients
(   81)       
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(   84)       REAL(wp) ::   &
(   85)          zabe1, zabe2, zbtr, zta, zsa     ! temporary scalars
(   86)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   87)          ztv, zsv,                      &  ! temporary workspace arrays
(   88)          ztdta, ztdsa                      !    "         "
(   89)       !!----------------------------------------------------------------------
(   90)       
(   91)       IF( kt == nit000 ) THEN
(   92)          IF(lwp) WRITE(numout,*)
(   93)          IF(lwp) WRITE(numout,*) 'tra_ldf_lap : iso-level laplacian diffusion'
(   94)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(   95)          ze1ur(:,:) = e2u(:,:) / e1u(:,:)
(   96)          ze2vr(:,:) = e1v(:,:) / e2v(:,:)
(   97)          zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
(   98)       ENDIF
(   99)       
(  100)       ! Save ta and sa trends
(  101)       IF( l_trdtra )   THEN
(  102)          ztdta(:,:,:) = ta(:,:,:) 
(  103)          ztdsa(:,:,:) = sa(:,:,:) 
(  104)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 3

(  105) 
(  106)       !                                                  ! =============
(  107)       DO jk = 1, jpkm1                                   ! Vertical slab
(  108)          !                                               ! =============
(  109)          ! 1. First derivative (gradient)
(  110)          ! -------------------
(  111)          DO jj = 1, jpjm1
(  112)             DO ji = 1, fs_jpim1   ! vector opt.
(  113) #if defined key_s_coord
(  114)                zabe1 = fsahtu(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj) * fse3u(ji,jj,jk)
(  115)                zabe2 = fsahtv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj) * fse3v(ji,jj,jk)
(  116) #else
(  117)                zabe1 = fsahtu(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj)
(  118)                zabe2 = fsahtv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj)
(  119) #endif
(  120)                ztu(ji,jj,jk) = zabe1 * ( tb(ji+1,jj  ,jk) - tb(ji,jj,jk) )
(  121)                zsu(ji,jj,jk) = zabe1 * ( sb(ji+1,jj  ,jk) - sb(ji,jj,jk) )
(  122)                ztv(ji,jj,jk) = zabe2 * ( tb(ji  ,jj+1,jk) - tb(ji,jj,jk) )
(  123)                zsv(ji,jj,jk) = zabe2 * ( sb(ji  ,jj+1,jk) - sb(ji,jj,jk) )
(  124)             END DO  
(  125)          END DO  
(  126)          
(  127)          
(  128)          ! 2. Second derivative (divergence)
(  129)          ! --------------------
(  130)          DO jj = 2, jpjm1
(  131)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  132) #if defined key_s_coord
(  133)                zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  134) #else
(  135)                zbtr = zbtr2(ji,jj)
(  136) #endif
(  137)                ! horizontal diffusive trends
(  138)                zta = zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
(  139)                   &          + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
(  140)                zsa = zbtr * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk)   &
(  141)                   &          + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
(  142)                ! add it to the general tracer trends
(  143)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  144)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  145)             END DO  
(  146)          END DO  
(  147)          !                                               ! =============
(  148)       END DO                                             !  End of slab  
(  149)       !                                                  ! =============
(  150) 
(  151)       ! save the trends for diagnostic
(  152)       ! save the horizontal diffusive trends
(  153)       IF( l_trdtra )   THEN
(  154)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  155)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  156) 
(  157)          CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
(  158)       ENDIF
(  159) 
(  160)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  161)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  162)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 4

(  163)       ENDIF
(  164) 
(  165)       ! "zonal" mean lateral diffusive heat and salt transport 
(  166)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  167) # if defined key_s_coord || defined key_partial_steps
(  168)          pht_ldf(:) = ptr_vj( ztv(:,:,:) )
(  169)          pst_ldf(:) = ptr_vj( zsv(:,:,:) )
(  170) # else
(  171)          DO jk = 1, jpkm1
(  172)             DO jj = 2, jpjm1
(  173)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  174)                  ztv(ji,jj,jk) = ztv(ji,jj,jk) * fse3v(ji,jj,jk)
(  175)                  zsv(ji,jj,jk) = zsv(ji,jj,jk) * fse3v(ji,jj,jk)
(  176)                END DO
(  177)             END DO
(  178)          END DO
(  179)          pht_ldf(:) = ptr_vj( ztv(:,:,:) )
(  180)          pst_ldf(:) = ptr_vj( zsv(:,:,:) )
(  181) # endif
(  182)       ENDIF
(  183) 
(  184)    END SUBROUTINE tra_ldf_lap
(  185) 
(  186)    !!==============================================================================
(  187) END MODULE traldf_lap







































PGF90 (Version     10.5)          09/23/2010  09:50:15      page 5

(    1) # 1 "traldf_lap.F90"
(    1) MODULE traldf_lap
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traldf_lap  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_ldf_lap  : update the tracer trend with the horizontal diffusion
(    9)    !!                 using a iso-level harmonic (laplacien) operator.
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   15)    USE trdmod          ! ocean active tracers trends 
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE diaptr          ! poleward transport diagnostics
(   19)    USE prtctl          ! Print control
(   20) 
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC tra_ldf_lap  ! routine called by step.F90
(   27) 
(   28)    !! * Substitutions
(   29) # 1 "./domzgr_substitute.h90"
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
(   48) # 30 "traldf_lap.F90"
(   30) # 1 "./ldftra_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** ldftra_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_substitute.h90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 6

(   19) !   'key_traldf_c3d' :                 aht: 3D coefficient
(   20) # 31 "traldf_lap.F90"
(   31) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 32 "traldf_lap.F90"
(   32) # 32
(   32)    !!----------------------------------------------------------------------
(   33)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   34)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_lap.F90,v 1.7 2005/09/02 15:45:33 opalod Exp $ 
(   35)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   36)    !!----------------------------------------------------------------------
(   37)    
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE tra_ldf_lap( kt )
(   41)       !!----------------------------------------------------------------------
(   42)       !!                  ***  ROUTINE tra_ldf_lap  ***
(   43)       !!                   
(   44)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   45)       !!      trend and add it to the general trend of tracer equation.
(   46)       !!
(   47)       !! ** Method  :   Second order diffusive operator evaluated using before
(   48)       !!      fields (forward time scheme). The horizontal diffusive trends of 
(   49)       !!      temperature (idem for salinity) is given by:
(   50)       !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
(   51)       !!      factors e3. are inside the derivatives:
(   52)       !!          difft = 1/(e1t*e2t*e3t) {  di-1[ aht e2u*e3u/e1u di(tb) ]
(   53)       !!                                   + dj-1[ aht e1v*e3v/e2v dj(tb) ] }
(   54)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(   55)       !!          difft = 1/(e1t*e2t) {  di-1[ aht e2u/e1u di(tb) ]
(   56)       !!                               + dj-1[ aht e1v/e2v dj(tb) ] }
(   57)       !!      Add this trend to the general tracer trend (ta,sa):
(   58)       !!          (ta,sa) = (ta,sa) + ( difft , diffs )
(   59)       !!
(   60)       !! ** Action  : - Update (ta,sa) arrays with the before iso-level 
(   61)       !!                harmonic mixing trend.
(   62)       !!              - Save the trends in (ztdta,ztdsa) ('key_trdtra')
(   63)       !!
(   64)       !! History :
(   65)       !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
(   66)       !!        !  91-11  (G. Madec)
(   67)       !!        !  95-11  (G. Madec)  suppress volumetric scale factors
(   68)       !!        !  96-01  (G. Madec)  statement function for e3
(   69)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   70)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   71)       !!----------------------------------------------------------------------
(   72)       USE oce              , ztu => ua,  &  ! use ua as workspace
(   73)          &                   zsu => va      ! use va as workspace






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 7

(   74) 
(   75)       !! * Arguments
(   76)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   77)       
(   78)       !! * Local save
(   79)       REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
(   80)          ze1ur, ze2vr, zbtr2              ! scale factor coefficients
(   81)       
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(   84)       REAL(wp) ::   &
(   85)          zabe1, zabe2, zbtr, zta, zsa     ! temporary scalars
(   86)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   87)          ztv, zsv,                      &  ! temporary workspace arrays
(   88)          ztdta, ztdsa                      !    "         "
(   89)       !!----------------------------------------------------------------------
(   90)       
(   91)       IF( kt == nit000 ) THEN
(   92)          IF(lwp) WRITE(numout,*)
(   93)          IF(lwp) WRITE(numout,*) 'tra_ldf_lap : iso-level laplacian diffusion'
(   94)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(   95)          ze1ur(:,:) = e2u(:,:) / e1u(:,:)
(   96)          ze2vr(:,:) = e1v(:,:) / e2v(:,:)
(   97)          zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
(   98)       ENDIF
(   99)       
(  100)       ! Save ta and sa trends
(  101)       IF( l_trdtra )   THEN
(  102)          ztdta(:,:,:) = ta(:,:,:) 
(  103)          ztdsa(:,:,:) = sa(:,:,:) 
(  104)       ENDIF
(  105) 
(  106)       !                                                  ! =============
(  107)       DO jk = 1, jpkm1                                   ! Vertical slab
(  108)          !                                               ! =============
(  109)          ! 1. First derivative (gradient)
(  110)          ! -------------------
(  111)          DO jj = 1, jpjm1
(  112)             DO ji = 1, jpim1   ! vector opt.
(  113) # 117
(  117)                zabe1 = ahtu(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj)
(  118)                zabe2 = ahtv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj)
(  119) # 120
(  120)                ztu(ji,jj,jk) = zabe1 * ( tb(ji+1,jj  ,jk) - tb(ji,jj,jk) )
(  121)                zsu(ji,jj,jk) = zabe1 * ( sb(ji+1,jj  ,jk) - sb(ji,jj,jk) )
(  122)                ztv(ji,jj,jk) = zabe2 * ( tb(ji  ,jj+1,jk) - tb(ji,jj,jk) )
(  123)                zsv(ji,jj,jk) = zabe2 * ( sb(ji  ,jj+1,jk) - sb(ji,jj,jk) )
(  124)             END DO  
(  125)          END DO  
(  126)          
(  127)          
(  128)          ! 2. Second derivative (divergence)
(  129)          ! --------------------
(  130)          DO jj = 2, jpjm1
(  131)             DO ji = 2, jpim1   ! vector opt.
(  132) # 135
(  135)                zbtr = zbtr2(ji,jj)
(  136) # 137






PGF90 (Version     10.5)          09/23/2010  09:50:15      page 8

(  137)                ! horizontal diffusive trends
(  138)                zta = zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
(  139)                   &          + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
(  140)                zsa = zbtr * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk)   &
(  141)                   &          + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
(  142)                ! add it to the general tracer trends
(  143)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  144)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  145)             END DO  
(  146)          END DO  
(  147)          !                                               ! =============
(  148)       END DO                                             !  End of slab  
(  149)       !                                                  ! =============
(  150) 
(  151)       ! save the trends for diagnostic
(  152)       ! save the horizontal diffusive trends
(  153)       IF( l_trdtra )   THEN
(  154)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  155)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  156) 
(  157)          CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
(  158)       ENDIF
(  159) 
(  160)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  161)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  162)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  163)       ENDIF
(  164) 
(  165)       ! "zonal" mean lateral diffusive heat and salt transport 
(  166)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  167) # 171
(  171)          DO jk = 1, jpkm1
(  172)             DO jj = 2, jpjm1
(  173)                DO ji = 2, jpim1   ! vector opt.
(  174)                  ztv(ji,jj,jk) = ztv(ji,jj,jk) * e3t(jk)
(  175)                  zsv(ji,jj,jk) = zsv(ji,jj,jk) * e3t(jk)
(  176)                END DO
(  177)             END DO
(  178)          END DO
(  179)          pht_ldf(:) = ptr_vj( ztv(:,:,:) )
(  180)          pst_ldf(:) = ptr_vj( zsv(:,:,:) )
(  181) # 182
(  182)       ENDIF
(  183) 
(  184)    END SUBROUTINE tra_ldf_lap
(  185) 
(  186)    !!==============================================================================
(  187) END MODULE traldf_lap
