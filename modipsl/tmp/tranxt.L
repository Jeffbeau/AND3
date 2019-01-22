


PGF90 (Version     10.5)          09/23/2010  09:50:16      page 1

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

Filename: tranxt.F90

(    1) MODULE tranxt
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  tranxt  ***
(    4)    !! Ocean active tracers:  time stepping on temperature and salinity
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_nxt     : time stepping on temperature and salinity
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE zdf_oce         ! ???
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   16)    USE obctra          ! open boundary condition (obc_tra routine)
(   17)    USE prtctl          ! Print control
(   18)    USE agrif_opa_update
(   19)    USE agrif_opa_interp
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC tra_nxt          ! routine called by step.F90
(   26)    !!----------------------------------------------------------------------
(   27)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tranxt.F90,v 1.7 2006/03/10 10:55:45 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!----------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE tra_nxt( kt )
(   35)       !!----------------------------------------------------------------------
(   36)       !!                   ***  ROUTINE tranxt  ***
(   37)       !!
(   38)       !! ** Purpose :   Compute the temperature and salinity fields at the 
(   39)       !!      next time-step from their temporal trends and swap the fields.
(   40)       !! 
(   41)       !! ** Method  :   Apply lateral boundary conditions on (ta,ta) through 
(   42)       !!      call to lbc_lnk routine
(   43)       !!      After t and s are compute using a leap-frog scheme environment:
(   44)       !!         ta = tb + 2 rdttra(k) * ta
(   45)       !!         sa = sb + 2 rdttra(k) * sa
(   46)       !!      Compute and save in (ta,sa) an average over three time levels






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 2

(   47)       !!      (before,now and after) of temperature and salinity which is
(   48)       !!      used to compute rhd in eos routine and thus the hydrostatic 
(   49)       !!      pressure gradient (ln_dynhpg_imp = T)
(   50)       !!      Apply an Asselin time filter on now tracers (tn,sn) to avoid
(   51)       !!      the divergence of two consecutive time-steps and swap tracer
(   52)       !!      arrays to prepare the next time_step:
(   53)       !!         (zt,zs) = (ta+2tn+tb,sa+2sn+sb)/4       (ln_dynhpg_imp = T)
(   54)       !!         (zt,zs) = (0,0)                            (default option)
(   55)       !!         (tb,sb) = (tn,vn) + atfp [ (tb,sb) + (ta,sa) - 2 (tn,sn) ]
(   56)       !!         (tn,sn) = (ta,sa) 
(   57)       !!         (ta,sa) = (zt,zs)  (NB: reset to 0 after use in eos.F)
(   58)       !!
(   59)       !! ** Action  : - update (tb,sb) and (tn,sn) 
(   60)       !!              - (ta,sa) time averaged (t,s)      (ln_dynhpg_imp = T)
(   61)       !!
(   62)       !! History :
(   63)       !!   7.0  !  91-11  (G. Madec)  Original code
(   64)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   65)       !!        !  96-02  (G. Madec & M. Imbard)  opa release 8.0
(   66)       !!   8.0  !  96-04  (A. Weaver)  Euler forward step
(   67)       !!   8.2  !  99-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
(   68)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   69)       !!        !  02-11  (C. Talandier, A-M Treguier) Open boundaries
(   70)       !!----------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   73) 
(   74)       !! * Local declarations
(   75)       INTEGER ::   ji, jj, jk   ! dummy loop indices
(   76)       REAL(wp) ::   zt, zs      ! temporary scalars
(   77)       REAL(wp) ::   zfact       ! temporary scalar
(   78)       !!----------------------------------------------------------------------
(   79) 
(   80) 
(   81)       ! 0. Lateral boundary conditions on ( ta, sa )   (T-point, unchanged sign)
(   82)       ! ---------------------------------============
(   83)       CALL lbc_lnk( ta, 'T', 1. )   
(   84)       CALL lbc_lnk( sa, 'T', 1. )
(   85) 
(   86) 
(   87)       !                                                ! ===============
(   88)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   89)          !                                             ! ===============
(   90) 
(   91)          ! 1. Leap-frog scheme (only in explicit case, otherwise the 
(   92)          ! -------------------  time stepping is already done in trazdf)
(   93)          IF( l_trazdf_exp ) THEN
(   94)             zfact = 2. * rdttra(jk)
(   95)             IF( neuler == 0 .AND. kt == nit000 ) zfact = rdttra(jk)
(   96)             ta(:,:,jk) = ( tb(:,:,jk) + zfact * ta(:,:,jk) ) * tmask(:,:,jk)
(   97)             sa(:,:,jk) = ( sb(:,:,jk) + zfact * sa(:,:,jk) ) * tmask(:,:,jk)
(   98) !           sa(:,:,jk) = max(sa(:,:,jk),5d0)
(   99)          ENDIF
(  100) 
(  101)          DO ji=1,jpi
(  102)          DO jj=1,jpj
(  103) 
(  104)             IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0 !BIO






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 3

(  105)          END DO
(  106)          END DO
(  107) 
(  108) 
(  109) #if defined key_obc
(  110)          !                                             ! ===============
(  111)       END DO                                           !   End of slab
(  112)       !                                                ! ===============
(  113) 
(  114)       ! Update tracers on open boundaries.
(  115)       CALL obc_tra( kt )
(  116) 
(  117)       !                                                ! ===============
(  118)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  119)          !                                             ! ===============
(  120) #endif
(  121) #if defined key_agrif
(  122)          !                                             ! ===============
(  123)       END DO                                           !   End of slab
(  124)       !                                                ! ===============
(  125) 
(  126)       ! Update tracers on open boundaries.
(  127)       CALL Agrif_tra( kt )
(  128) 
(  129)       !                                                ! ===============
(  130)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  131)          !                                             ! ===============
(  132) #endif
(  133) 
(  134) 
(  135)          ! 2. Time filter and swap of arrays
(  136)          ! ---------------------------------
(  137)  
(  138)          IF( ln_dynhpg_imp ) THEN                       ! semi-implicite hpg
(  139)             IF( neuler == 0 .AND. kt == nit000 ) THEN
(  140)                DO jj = 1, jpj
(  141)                   DO ji = 1, jpi
(  142)                      zt = ( ta(ji,jj,jk) + 2. * tn(ji,jj,jk) + tb(ji,jj,jk) ) * 0.25
(  143)                      zs = ( sa(ji,jj,jk) + 2. * sn(ji,jj,jk) + sb(ji,jj,jk) ) * 0.25
(  144)                      tb(ji,jj,jk) = tn(ji,jj,jk)
(  145)                      sb(ji,jj,jk) = sn(ji,jj,jk)
(  146)                      tn(ji,jj,jk) = ta(ji,jj,jk)
(  147)                      IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
(  148)                      sn(ji,jj,jk) = sa(ji,jj,jk)
(  149) !                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
(  150)                      ta(ji,jj,jk) = zt
(  151)                      IF(zs.lt.0.0) zs=0.0 !BIO
(  152)                      sa(ji,jj,jk) = zs
(  153) !                    sa(ji,jj,jk) = max(zs,5d0)
(  154)                   END DO
(  155)                END DO
(  156)             ELSE
(  157)                DO jj = 1, jpj
(  158)                   DO ji = 1, jpi
(  159)                      zt = ( ta(ji,jj,jk) + 2. * tn(ji,jj,jk) + tb(ji,jj,jk) ) * 0.25
(  160)                      zs = ( sa(ji,jj,jk) + 2. * sn(ji,jj,jk) + sb(ji,jj,jk) ) * 0.25
(  161)                      tb(ji,jj,jk) = atfp  * ( tb(ji,jj,jk) + ta(ji,jj,jk) ) + atfp1 * tn(ji,jj,jk)
(  162)                      sb(ji,jj,jk) = atfp  * ( sb(ji,jj,jk) + sa(ji,jj,jk) ) + atfp1 * sn(ji,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 4

(  163)                      tn(ji,jj,jk) = ta(ji,jj,jk)
(  164)                     IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
(  165)                      sn(ji,jj,jk) = sa(ji,jj,jk)
(  166) !                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
(  167)                      ta(ji,jj,jk) = zt
(  168)                      IF(zs.lt.0.0) zs=0.0 !BIO
(  169)                      sa(ji,jj,jk) = zs
(  170) !                    sa(ji,jj,jk) = max(zs,5d0)
(  171)                   END DO
(  172)                END DO
(  173)             ENDIF
(  174)          ELSE                                          ! Default case
(  175)             IF( neuler == 0 .AND. kt == nit000 ) THEN
(  176)                DO jj = 1, jpj
(  177)                   DO ji = 1, jpi
(  178)                      tb(ji,jj,jk) = tn(ji,jj,jk)
(  179)                      sb(ji,jj,jk) = sn(ji,jj,jk)
(  180)                      tn(ji,jj,jk) = ta(ji,jj,jk)
(  181)                      IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
(  182)                      sn(ji,jj,jk) = sa(ji,jj,jk)
(  183) !                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
(  184)                   END DO
(  185)                END DO
(  186)             ELSE
(  187)                DO jj = 1, jpj
(  188)                   DO ji = 1, jpi
(  189)                      tb(ji,jj,jk) = atfp  * ( tb(ji,jj,jk) + ta(ji,jj,jk) ) + atfp1 * tn(ji,jj,jk)
(  190)                      sb(ji,jj,jk) = atfp  * ( sb(ji,jj,jk) + sa(ji,jj,jk) ) + atfp1 * sn(ji,jj,jk)
(  191)                      tn(ji,jj,jk) = ta(ji,jj,jk)
(  192)                      IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
(  193)                      sn(ji,jj,jk) = sa(ji,jj,jk)
(  194) !                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
(  195)                   END DO
(  196)                END DO
(  197)             ENDIF
(  198)          ENDIF
(  199)          !                                             ! ===============
(  200)       END DO                                           !   End of slab
(  201)       !                                                ! ===============
(  202) 
(  203)       IF(ln_ctl) THEN         ! print mean field (used for debugging)
(  204)          CALL prt_ctl(tab3d_1=tn, clinfo1=' nxt  - Tn: ', mask1=tmask, &
(  205)             &         tab3d_2=sn, clinfo2=' Sn: ', mask2=tmask)
(  206)       ENDIF
(  207)       
(  208) #if defined key_agrif
(  209)       IF (.NOT.Agrif_Root())    CALL Agrif_Update_Tra( kt )
(  210) #endif      
(  211) 
(  212)    END SUBROUTINE tra_nxt
(  213) 
(  214)    !!======================================================================
(  215) END MODULE tranxt











PGF90 (Version     10.5)          09/23/2010  09:50:16      page 5

(    1) # 1 "tranxt.F90"
(    1) MODULE tranxt
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  tranxt  ***
(    4)    !! Ocean active tracers:  time stepping on temperature and salinity
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_nxt     : time stepping on temperature and salinity
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE zdf_oce         ! ???
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   16)    USE obctra          ! open boundary condition (obc_tra routine)
(   17)    USE prtctl          ! Print control
(   18)    USE agrif_opa_update
(   19)    USE agrif_opa_interp
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC tra_nxt          ! routine called by step.F90
(   26)    !!----------------------------------------------------------------------
(   27)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tranxt.F90,v 1.7 2006/03/10 10:55:45 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!----------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE tra_nxt( kt )
(   35)       !!----------------------------------------------------------------------
(   36)       !!                   ***  ROUTINE tranxt  ***
(   37)       !!
(   38)       !! ** Purpose :   Compute the temperature and salinity fields at the 
(   39)       !!      next time-step from their temporal trends and swap the fields.
(   40)       !! 
(   41)       !! ** Method  :   Apply lateral boundary conditions on (ta,ta) through 
(   42)       !!      call to lbc_lnk routine
(   43)       !!      After t and s are compute using a leap-frog scheme environment:
(   44)       !!         ta = tb + 2 rdttra(k) * ta
(   45)       !!         sa = sb + 2 rdttra(k) * sa
(   46)       !!      Compute and save in (ta,sa) an average over three time levels
(   47)       !!      (before,now and after) of temperature and salinity which is
(   48)       !!      used to compute rhd in eos routine and thus the hydrostatic 
(   49)       !!      pressure gradient (ln_dynhpg_imp = T)
(   50)       !!      Apply an Asselin time filter on now tracers (tn,sn) to avoid
(   51)       !!      the divergence of two consecutive time-steps and swap tracer
(   52)       !!      arrays to prepare the next time_step:
(   53)       !!         (zt,zs) = (ta+2tn+tb,sa+2sn+sb)/4       (ln_dynhpg_imp = T)
(   54)       !!         (zt,zs) = (0,0)                            (default option)
(   55)       !!         (tb,sb) = (tn,vn) + atfp [ (tb,sb) + (ta,sa) - 2 (tn,sn) ]
(   56)       !!         (tn,sn) = (ta,sa) 
(   57)       !!         (ta,sa) = (zt,zs)  (NB: reset to 0 after use in eos.F)






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 6

(   58)       !!
(   59)       !! ** Action  : - update (tb,sb) and (tn,sn) 
(   60)       !!              - (ta,sa) time averaged (t,s)      (ln_dynhpg_imp = T)
(   61)       !!
(   62)       !! History :
(   63)       !!   7.0  !  91-11  (G. Madec)  Original code
(   64)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   65)       !!        !  96-02  (G. Madec & M. Imbard)  opa release 8.0
(   66)       !!   8.0  !  96-04  (A. Weaver)  Euler forward step
(   67)       !!   8.2  !  99-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
(   68)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   69)       !!        !  02-11  (C. Talandier, A-M Treguier) Open boundaries
(   70)       !!----------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   73) 
(   74)       !! * Local declarations
(   75)       INTEGER ::   ji, jj, jk   ! dummy loop indices
(   76)       REAL(wp) ::   zt, zs      ! temporary scalars
(   77)       REAL(wp) ::   zfact       ! temporary scalar
(   78)       !!----------------------------------------------------------------------
(   79) 
(   80) 
(   81)       ! 0. Lateral boundary conditions on ( ta, sa )   (T-point, unchanged sign)
(   82)       ! ---------------------------------============
(   83)       CALL lbc_lnk( ta, 'T', 1. )   
(   84)       CALL lbc_lnk( sa, 'T', 1. )
(   85) 
(   86) 
(   87)       !                                                ! ===============
(   88)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   89)          !                                             ! ===============
(   90) 
(   91)          ! 1. Leap-frog scheme (only in explicit case, otherwise the 
(   92)          ! -------------------  time stepping is already done in trazdf)
(   93)          IF( l_trazdf_exp ) THEN
(   94)             zfact = 2. * rdttra(jk)
(   95)             IF( neuler == 0 .AND. kt == nit000 ) zfact = rdttra(jk)
(   96)             ta(:,:,jk) = ( tb(:,:,jk) + zfact * ta(:,:,jk) ) * tmask(:,:,jk)
(   97)             sa(:,:,jk) = ( sb(:,:,jk) + zfact * sa(:,:,jk) ) * tmask(:,:,jk)
(   98) !           sa(:,:,jk) = max(sa(:,:,jk),5d0)
(   99)          ENDIF
(  100) 
(  101)          DO ji=1,jpi
(  102)          DO jj=1,jpj
(  103) 
(  104)             IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0 !BIO
(  105)          END DO
(  106)          END DO
(  107) 
(  108) 
(  109) # 110
(  110)          !                                             ! ===============
(  111)       END DO                                           !   End of slab
(  112)       !                                                ! ===============
(  113) 
(  114)       ! Update tracers on open boundaries.
(  115)       CALL obc_tra( kt )






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 7

(  116) 
(  117)       !                                                ! ===============
(  118)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  119)          !                                             ! ===============
(  120) 
(  121) 
(  122) # 135
(  135)          ! 2. Time filter and swap of arrays
(  136)          ! ---------------------------------
(  137)  
(  138)          IF( ln_dynhpg_imp ) THEN                       ! semi-implicite hpg
(  139)             IF( neuler == 0 .AND. kt == nit000 ) THEN
(  140)                DO jj = 1, jpj
(  141)                   DO ji = 1, jpi
(  142)                      zt = ( ta(ji,jj,jk) + 2. * tn(ji,jj,jk) + tb(ji,jj,jk) ) * 0.25
(  143)                      zs = ( sa(ji,jj,jk) + 2. * sn(ji,jj,jk) + sb(ji,jj,jk) ) * 0.25
(  144)                      tb(ji,jj,jk) = tn(ji,jj,jk)
(  145)                      sb(ji,jj,jk) = sn(ji,jj,jk)
(  146)                      tn(ji,jj,jk) = ta(ji,jj,jk)
(  147)                      IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
(  148)                      sn(ji,jj,jk) = sa(ji,jj,jk)
(  149) !                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
(  150)                      ta(ji,jj,jk) = zt
(  151)                      IF(zs.lt.0.0) zs=0.0 !BIO
(  152)                      sa(ji,jj,jk) = zs
(  153) !                    sa(ji,jj,jk) = max(zs,5d0)
(  154)                   END DO
(  155)                END DO
(  156)             ELSE
(  157)                DO jj = 1, jpj
(  158)                   DO ji = 1, jpi
(  159)                      zt = ( ta(ji,jj,jk) + 2. * tn(ji,jj,jk) + tb(ji,jj,jk) ) * 0.25
(  160)                      zs = ( sa(ji,jj,jk) + 2. * sn(ji,jj,jk) + sb(ji,jj,jk) ) * 0.25
(  161)                      tb(ji,jj,jk) = atfp  * ( tb(ji,jj,jk) + ta(ji,jj,jk) ) + atfp1 * tn(ji,jj,jk)
(  162)                      sb(ji,jj,jk) = atfp  * ( sb(ji,jj,jk) + sa(ji,jj,jk) ) + atfp1 * sn(ji,jj,jk)
(  163)                      tn(ji,jj,jk) = ta(ji,jj,jk)
(  164)                     IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
(  165)                      sn(ji,jj,jk) = sa(ji,jj,jk)
(  166) !                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
(  167)                      ta(ji,jj,jk) = zt
(  168)                      IF(zs.lt.0.0) zs=0.0 !BIO
(  169)                      sa(ji,jj,jk) = zs
(  170) !                    sa(ji,jj,jk) = max(zs,5d0)
(  171)                   END DO
(  172)                END DO
(  173)             ENDIF
(  174)          ELSE                                          ! Default case
(  175)             IF( neuler == 0 .AND. kt == nit000 ) THEN
(  176)                DO jj = 1, jpj
(  177)                   DO ji = 1, jpi
(  178)                      tb(ji,jj,jk) = tn(ji,jj,jk)
(  179)                      sb(ji,jj,jk) = sn(ji,jj,jk)
(  180)                      tn(ji,jj,jk) = ta(ji,jj,jk)
(  181)                      IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
(  182)                      sn(ji,jj,jk) = sa(ji,jj,jk)
(  183) !                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
(  184)                   END DO
(  185)                END DO






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 8

(  186)             ELSE
(  187)                DO jj = 1, jpj
(  188)                   DO ji = 1, jpi
(  189)                      tb(ji,jj,jk) = atfp  * ( tb(ji,jj,jk) + ta(ji,jj,jk) ) + atfp1 * tn(ji,jj,jk)
(  190)                      sb(ji,jj,jk) = atfp  * ( sb(ji,jj,jk) + sa(ji,jj,jk) ) + atfp1 * sn(ji,jj,jk)
(  191)                      tn(ji,jj,jk) = ta(ji,jj,jk)
(  192)                      IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
(  193)                      sn(ji,jj,jk) = sa(ji,jj,jk)
(  194) !                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
(  195)                   END DO
(  196)                END DO
(  197)             ENDIF
(  198)          ENDIF
(  199)          !                                             ! ===============
(  200)       END DO                                           !   End of slab
(  201)       !                                                ! ===============
(  202) 
(  203)       IF(ln_ctl) THEN         ! print mean field (used for debugging)
(  204)          CALL prt_ctl(tab3d_1=tn, clinfo1=' nxt  - Tn: ', mask1=tmask, &
(  205)             &         tab3d_2=sn, clinfo2=' Sn: ', mask2=tmask)
(  206)       ENDIF
(  207)       
(  208) 
(  209) # 212
(  212)    END SUBROUTINE tra_nxt
(  213) 
(  214)    !!======================================================================
(  215) END MODULE tranxt
