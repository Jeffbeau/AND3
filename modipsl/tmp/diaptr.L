


PGF90 (Version     10.5)          09/23/2010  09:49:28      page 1

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

Filename: diaptr.F90

(    1) MODULE diaptr
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  diaptr  ***
(    4)    !! Ocean physics:  brief description of the purpose of the module
(    5)    !!                 (please no more than 2 lines)
(    6)    !!=====================================================================
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dia_ptr      : Poleward Transport Diagnostics module
(    9)    !!   dia_ptr_init : Initialization, namelist read
(   10)    !!   dia_ptr_wri  : Output of poleward fluxes
(   11)    !!   ptr_vjk      : "zonal" sum computation of a "meridional" flux array
(   12)    !!   ptr_vtjk     : "zonal" mean computation of a tracer field
(   13)    !!   ptr_vj       : "zonal" and vertical sum computation of a "meridional"
(   14)    !!                : flux array; Generic interface: ptr_vj_3d, ptr_vj_2d
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce           ! ocean dynamics and active tracers
(   18)    USE dom_oce       ! ocean space and time domain
(   19)    USE ldftra_oce    ! ???
(   20)    USE lib_mpp
(   21)    USE in_out_manager
(   22)    USE dianam
(   23)    USE phycst
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    INTERFACE ptr_vj
(   29)       MODULE PROCEDURE ptr_vj_3d, ptr_vj_2d
(   30)    END INTERFACE
(   31) 
(   32)    !! *  Routine accessibility
(   33)    PUBLIC dia_ptr_init   ! call in opa module
(   34)    PUBLIC dia_ptr        ! call in step module
(   35)    PUBLIC ptr_vj         ! call by tra_ldf & tra_adv routines
(   36)    PUBLIC ptr_vjk        ! call by tra_ldf & tra_adv routines
(   37) 
(   38)    !! * Share Module variables
(   39)    LOGICAL, PUBLIC ::       & !!! ** init namelist (namptr) **
(   40)       ln_diaptr = .FALSE.,  &  !: Poleward transport flag (T) or not (F)
(   41)       ln_subbas = .FALSE.      !: Atlantic/Pacific/Indian basins calculation
(   42)    INTEGER, PUBLIC ::       & !!: ** ptr namelist (namptr) **
(   43)       nf_ptr = 15              !: frequency of ptr computation
(   44)    REAL(wp), PUBLIC, DIMENSION(jpj) ::   &   !!: poleward transport
(   45)       pht_adv, pst_adv,     &  !: heat and salt: advection
(   46)       pht_ove, pst_ove,     &  !: heat and salt: overturning






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 2

(   47)       pht_ldf, pst_ldf,     &  !: heat and salt: lateral diffusion
(   48) #if defined key_diaeiv
(   49)       pht_eiv, pst_eiv,     &  !: heat and salt: bolus advection
(   50) #endif
(   51)       ht_atl,ht_ind,ht_pac, &  !: heat
(   52)       st_atl,st_ind,st_pac     !: salt
(   53)    REAL(wp),DIMENSION(jpi,jpj) ::   &
(   54)       abasin,pbasin,ibasin     !: return function value
(   55)      
(   56) 
(   57)    !! Module variables
(   58)    REAL(wp), DIMENSION(jpj,jpk) ::   &  
(   59)       tn_jk  , sn_jk  ,  &  !: "zonal" mean temperature and salinity
(   60)       v_msf_atl       ,  &  !: "meridional" Stream-Function
(   61)       v_msf_glo       ,  &  !: "meridional" Stream-Function
(   62)       v_msf_ipc       ,  &  !: "meridional" Stream-Function
(   63) #if defined key_diaeiv
(   64)       v_msf_eiv       ,  &  !: bolus "meridional" Stream-Function
(   65) #endif
(   66)       surf_jk_r             !: inverse of the ocean "zonal" section surface
(   67) 
(   68)    !! * Substitutions
(   69) #  include "domzgr_substitute.h90"
(   70) #  include "vectopt_loop_substitute.h90"
(   71)    !!----------------------------------------------------------------------
(   72)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   73)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diaptr.F90,v 1.10 2006/03/20 16:05:53 opalod Exp $ 
(   74)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   75)    !!----------------------------------------------------------------------
(   76) 
(   77) CONTAINS
(   78) 
(   79)    FUNCTION ptr_vj_3d( pva )   RESULT ( p_fval )
(   80)       !!----------------------------------------------------------------------
(   81)       !!                    ***  ROUTINE ptr_vj_3d  ***
(   82)       !!
(   83)       !! ** Purpose :   "zonal" and vertical sum computation of a "meridional"
(   84)       !!      flux array
(   85)       !!
(   86)       !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
(   87)       !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
(   88)       !!
(   89)       !! ** Action  : - p_fval: i-k-mean poleward flux of pva
(   90)       !!
(   91)       !! History :
(   92)       !!   9.0  !  03-09  (G. Madec)  Original code
(   93)       !!----------------------------------------------------------------------
(   94)       !! * arguments
(   95)       REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
(   96)          pva                         ! mask flux array at V-point
(   97) 
(   98)       !! * local declarations
(   99)       INTEGER  ::   ji, jj, jk        ! dummy loop arguments
(  100) #if ! defined key_agrif
(  101)       INTEGER  ::   ijpj = jpj        ! ???
(  102) #else
(  103)       INTEGER  ::   ijpj             ! ???
(  104) #endif      






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 3

(  105)       REAL(wp),DIMENSION(jpj) ::   &
(  106)          p_fval                       ! function value
(  107)       !!--------------------------------------------------------------------
(  108) #if defined key_agrif
(  109)       ijpj = jpj
(  110) #endif      
(  111) 
(  112)       p_fval(:) = 0.e0
(  113)       DO jk = 1, jpkm1
(  114)          DO jj = 2, jpjm1
(  115)             DO ji = fs_2, fs_jpim1   ! Vector opt.
(  116)                p_fval(jj) = p_fval(jj) + pva(ji,jj,jk) * tmask_i(ji,jj+1) * tmask_i(ji,jj) 
(  117)             END DO
(  118)          END DO
(  119)       END DO
(  120) 
(  121)       IF( lk_mpp )   CALL mpp_sum( p_fval, ijpj )     !!bug  I presume
(  122) 
(  123)    END FUNCTION ptr_vj_3d
(  124) 
(  125) 
(  126) 
(  127)    FUNCTION ptr_vj_2d( pva )   RESULT ( p_fval )
(  128)       !!----------------------------------------------------------------------
(  129)       !!                    ***  ROUTINE ptr_vj_2d  ***
(  130)       !!
(  131)       !! ** Purpose :   "zonal" and vertical sum computation of a "meridional"
(  132)       !!      flux array
(  133)       !!
(  134)       !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
(  135)       !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
(  136)       !!
(  137)       !! ** Action  : - p_fval: i-k-mean poleward flux of pva
(  138)       !!
(  139)       !! History :
(  140)       !!   9.0  !  03-09  (G. Madec)  Original code
(  141)       !!----------------------------------------------------------------------
(  142)       !! * arguments
(  143)       REAL(wp) , INTENT(in), DIMENSION(jpi,jpj) ::   &
(  144)          pva                         ! mask flux array at V-point
(  145) 
(  146)       !! * local declarations
(  147)       INTEGER  ::   ji,jj             ! dummy loop arguments
(  148) #if ! defined key_agrif
(  149)       INTEGER  ::   ijpj = jpj        ! ???
(  150) #else
(  151)       INTEGER  ::   ijpj             ! ???
(  152) #endif      
(  153)       REAL(wp),DIMENSION(jpj) ::   &
(  154)          p_fval                       ! function value
(  155)       !!--------------------------------------------------------------------
(  156) #if defined key_agrif
(  157)       ijpj = jpj
(  158) #endif      
(  159)  
(  160)       p_fval(:) = 0.e0
(  161)       DO jj = 2, jpjm1
(  162)          DO ji = fs_2, fs_jpim1   ! Vector opt.






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 4

(  163)             p_fval(jj) = p_fval(jj) + pva(ji,jj) * tmask_i(ji,jj+1) * tmask_i(ji,jj)
(  164)          END DO
(  165)       END DO
(  166) 
(  167)       IF( lk_mpp )   CALL mpp_sum( p_fval, ijpj )     !!bug  I presume
(  168)  
(  169)     END FUNCTION ptr_vj_2d
(  170) 
(  171) 
(  172) 
(  173)    FUNCTION ptr_vjk( pva )   RESULT ( p_fval )
(  174)       !!----------------------------------------------------------------------
(  175)       !!                    ***  ROUTINE ptr_vjk  ***
(  176)       !!
(  177)       !! ** Purpose :   "zonal" sum computation of a "meridional" flux array
(  178)       !!
(  179)       !! ** Method  : - i-sum of pva using the interior 2D vmask (vmask_i).
(  180)       !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
(  181)       !!
(  182)       !! ** Action  : - p_fval: i-k-mean poleward flux of pva
(  183)       !!
(  184)       !! History :
(  185)       !!   9.0  !  03-09  (G. Madec)  Original code
(  186)       !!----------------------------------------------------------------------
(  187)       !! * arguments
(  188)       REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
(  189)          pva                         ! mask flux array at V-point
(  190) 
(  191)       !! * local declarations
(  192)       INTEGER  ::   ji, jj, jk        ! dummy loop arguments
(  193)       INTEGER, DIMENSION (1) :: ish
(  194)       INTEGER, DIMENSION (2) :: ish2
(  195)       REAL(wp),DIMENSION(jpj*jpk) ::   &
(  196)          zwork                        ! temporary vector for mpp_sum
(  197)       REAL(wp),DIMENSION(jpj,jpk) ::   &
(  198)          p_fval                       ! return function value
(  199)       !!--------------------------------------------------------------------
(  200)  
(  201)       p_fval(:,:) = 0.e0
(  202) 
(  203)       DO jk = 1, jpkm1
(  204)          DO jj = 2, jpjm1
(  205)            DO ji = fs_2, fs_jpim1
(  206)             p_fval(jj,jk) = p_fval(jj,jk) + pva(ji,jj,jk) * e1v(ji,jj) * fse3v(ji,jj,jk) &
(  207)                &            * tmask_i(ji,jj+1) * tmask_i(ji,jj)
(  208)            END DO
(  209)          END DO
(  210)       END DO
(  211) 
(  212)       IF(lk_mpp)   THEN
(  213)          ish(1) = jpj*jpk ; ish2(1)=jpj ; ish2(2)=jpk
(  214)          zwork(:)= RESHAPE(p_fval, ish )
(  215)          CALL mpp_sum(zwork, jpj*jpk )
(  216)          p_fval(:,:)= RESHAPE(zwork,ish2)
(  217)       END IF
(  218) 
(  219)    END FUNCTION ptr_vjk
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 5

(  221)    FUNCTION ptr_vtjk( pva )   RESULT ( p_fval )
(  222)       !!----------------------------------------------------------------------
(  223)       !!                    ***  ROUTINE ptr_vtjk  ***
(  224)       !!
(  225)       !! ** Purpose :   "zonal" mean computation of a tracer field
(  226)       !!
(  227)       !! ** Method  : - i-sum of mj(pva) using the interior 2D vmask (vmask_i)
(  228)       !!      multiplied by the inverse of the surface of the "zonal" ocean
(  229)       !!      section
(  230)       !!
(  231)       !! ** Action  : - p_fval: i-k-mean poleward flux of pva
(  232)       !!
(  233)       !! History :
(  234)       !!   9.0  !  03-09  (G. Madec)  Original code
(  235)       !!   9.0  !  06-01  (A. Biastoch)  Allow sub-basins computation
(  236)       !!----------------------------------------------------------------------
(  237)       !! * arguments
(  238)       REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
(  239)          pva                         ! mask flux array at V-point
(  240)  
(  241)       !! * local declarations
(  242)       INTEGER  ::   ji, jj, jk        ! dummy loop arguments
(  243)       INTEGER, DIMENSION (1) :: ish
(  244)       INTEGER, DIMENSION (2) :: ish2
(  245)       REAL(wp),DIMENSION(jpj*jpk) ::   &
(  246)          zwork                        ! temporary vector for mpp_sum
(  247)       REAL(wp),DIMENSION(jpj,jpk) ::   &
(  248)          p_fval                       ! return function value
(  249)       !!-------------------------------------------------------------------- 
(  250) 
(  251)       p_fval(:,:) = 0.e0
(  252)       DO jk = 1, jpkm1
(  253)          DO jj = 2, jpjm1
(  254)             DO ji = fs_2, fs_jpim1   ! Vector opt.
(  255)                p_fval(jj,jk) = p_fval(jj,jk) + ( pva(ji,jj,jk) + pva(ji,jj+1,jk) )              &
(  256)                   &                          * e1v(ji,jj) * fse3v(ji,jj,jk) * vmask(ji,jj,jk)   &
(  257)                   &                          * tmask_i(ji,jj+1) * tmask_i(ji,jj)
(  258)             END DO
(  259)          END DO
(  260)       END DO
(  261)       p_fval(:,:) = p_fval(:,:) * 0.5
(  262)       IF(lk_mpp)   THEN
(  263)          ish(1) = jpj*jpk ; ish2(1)=jpj ; ish2(2)=jpk
(  264)          zwork(:)= RESHAPE(p_fval, ish )
(  265)          CALL mpp_sum(zwork, jpj*jpk )
(  266)          p_fval(:,:)= RESHAPE(zwork,ish2)
(  267)       END IF
(  268) 
(  269)    END FUNCTION ptr_vtjk
(  270) 
(  271) 
(  272)    SUBROUTINE dia_ptr( kt )
(  273)       !!----------------------------------------------------------------------
(  274)       !!                  ***  ROUTINE dia_ptr  ***
(  275)       !!----------------------------------------------------------------------
(  276)       !! * Moudules used
(  277)       USE ioipsl
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 6

(  279)       !! * Argument
(  280)       INTEGER, INTENT(in) ::   kt   ! ocean time step index
(  281) 
(  282)       !! * Local variables
(  283)       INTEGER ::   jk,jj,ji               ! dummy loop
(  284)       REAL(wp) ::    &
(  285)          zsverdrup,  &              ! conversion from m3/s to Sverdrup
(  286)          zpwatt,     &              ! conversion from W    to PW
(  287)          zggram                     ! conversion from g    to Pg
(  288) 
(  289)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &
(  290)          v_atl , v_ipc,                    &
(  291)          vt_atl, vt_pac, vt_ind,           &
(  292)          vs_atl, vs_pac, vs_ind,           &
(  293)          zv_eiv
(  294)       CHARACTER (len=32) ::   &
(  295)          clnam = 'subbasins.nc'                
(  296)       INTEGER ::  itime,inum,ipi,ipj,ipk       ! temporary integer
(  297)       INTEGER, DIMENSION (1) ::   istep
(  298)       REAL(wp) ::    zdate0,zsecond,zdt        ! temporary scalars
(  299)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(  300)          zlamt, zphit, zdta             ! temporary workspace (NetCDF read)
(  301)       REAL(wp), DIMENSION(jpk) ::   &
(  302)          zdept                          ! temporary workspace (NetCDF read)
(  303)       !!----------------------------------------------------------------------
(  304) 
(  305)       IF( kt == nit000 .OR. MOD( kt, nf_ptr ) == 0 )   THEN
(  306) 
(  307)          zsverdrup = 1.e-6
(  308)          zpwatt    = 1.e-15
(  309)          zggram    = 1.e-6
(  310)          ipi       = jpidta
(  311)          ipj       = jpjdta
(  312)          ipk       = 1
(  313)          itime     = 1
(  314)          zsecond   = 0.e0
(  315)          zdate0    = 0.e0
(  316)    
(  317) # if defined key_diaeiv
(  318)          zv_eiv(:,:,:) = v_eiv(:,:,:)
(  319) # else
(  320)          zv_eiv(:,:,:) = 0.e0
(  321) # endif
(  322) 
(  323)          ! "zonal" mean temperature and salinity at V-points
(  324)          tn_jk(:,:) = ptr_vtjk( tn(:,:,:) ) * surf_jk_r(:,:)
(  325)          sn_jk(:,:) = ptr_vtjk( sn(:,:,:) ) * surf_jk_r(:,:)
(  326) 
(  327)          !--------------------------------------------------------
(  328)          ! overturning calculation:
(  329)  
(  330)          IF( ln_subbas ) THEN              ! Basins computation
(  331) 
(  332)             IF( kt == nit000 ) THEN                ! load basin mask
(  333)                itime = 1
(  334)                ipi   = jpidta
(  335)                ipj   = jpjdta
(  336)                ipk   = 1






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 7

(  337)                zdt   = 0.e0
(  338)                istep = 0
(  339)                clnam = 'subbasins.nc'
(  340) 
(  341)                CALL flinopen(clnam,1,jpidta,1,jpjdta,.FALSE.,ipi,ipj, &
(  342)                   &          ipk,zlamt,zphit,zdept,itime,istep,zdate0,zdt,inum)
(  343) 
(  344)                ! get basins:
(  345)                abasin (:,:) = 0.e0
(  346)                pbasin (:,:) = 0.e0
(  347)                ibasin (:,:) = 0.e0
(  348) 
(  349)                ! Atlantic basin
(  350)                CALL flinget(inum,'atlmsk',jpidta,jpjdta,1,itime,1,   &
(  351)                   &         0,1,jpidta,1,jpjdta,zdta(:,:))
(  352)                DO jj = 1, nlcj                                 ! interior values
(  353)                   DO ji = 1, nlci
(  354)                      abasin (ji,jj) = zdta( mig(ji), mjg(jj) )
(  355)                   END DO
(  356)                END DO
(  357) 
(  358)                ! Pacific basin
(  359)                CALL flinget(inum,'pacmsk',jpidta,jpjdta,1,itime,1,   &
(  360)                   &         0,1,jpidta,1,jpjdta,zdta(:,:))
(  361)                DO jj = 1, nlcj                                 ! interior values
(  362)                   DO ji = 1, nlci
(  363)                      pbasin (ji,jj) = zdta( mig(ji), mjg(jj) )
(  364)                   END DO
(  365)                END DO
(  366) 
(  367)                ! Indian basin
(  368)                CALL flinget(inum,'indmsk',jpidta,jpjdta,1,itime,1,   &
(  369)                   &         0,1,jpidta,1,jpjdta,zdta(:,:))
(  370)                DO jj = 1, nlcj                                 ! interior values
(  371)                   DO ji = 1, nlci
(  372)                      ibasin (ji,jj) = zdta( mig(ji), mjg(jj) )
(  373)                   END DO
(  374)                END DO
(  375) 
(  376)                CALL flinclo(inum)
(  377) 
(  378)             ENDIF
(  379) 
(  380)             ! basin separation:
(  381)             DO jj = 1, jpj
(  382)                DO ji = 1, jpi
(  383)                   ! basin separated velocity
(  384)                   v_atl(ji,jj,:) = (vn(ji,jj,:)+zv_eiv(ji,jj,:))*abasin(ji,jj)
(  385)                   v_ipc(ji,jj,:) = (vn(ji,jj,:)+zv_eiv(ji,jj,:))*(pbasin(ji,jj)+ibasin(ji,jj))
(  386) 
(  387)                   ! basin separated T times V on T points
(  388)                   vt_ind(ji,jj,:) = tn(ji,jj,:) *                                 &
(  389)                      &              ( (vn    (ji,jj,:) + vn    (ji,jj-1,:))*0.5   &
(  390)                      &              + (zv_eiv(ji,jj,:) + zv_eiv(ji,jj-1,:))*0.5 ) 
(  391)                   vt_atl(ji,jj,:) = vt_ind(ji,jj,:) * abasin(ji,jj)
(  392)                   vt_pac(ji,jj,:) = vt_ind(ji,jj,:) * pbasin(ji,jj)
(  393)                   vt_ind(ji,jj,:) = vt_ind(ji,jj,:) * ibasin(ji,jj)
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 8

(  395)                   ! basin separated S times V on T points
(  396)                   vs_ind(ji,jj,:) = sn(ji,jj,:) *                                 &
(  397)                      &              ( (vn    (ji,jj,:) + vn    (ji,jj-1,:))*0.5   &
(  398)                      &              + (zv_eiv(ji,jj,:) + zv_eiv(ji,jj-1,:))*0.5 ) 
(  399)                   vs_atl(ji,jj,:) = vs_ind(ji,jj,:) * abasin(ji,jj)
(  400)                   vs_pac(ji,jj,:) = vs_ind(ji,jj,:) * pbasin(ji,jj)
(  401)                   vs_ind(ji,jj,:) = vs_ind(ji,jj,:) * ibasin(ji,jj)
(  402)                END DO
(  403)             END DO
(  404) 
(  405)          ENDIF
(  406) 
(  407)          ! horizontal integral and vertical dz 
(  408)          v_msf_glo(:,:) = ptr_vjk( vn(:,:,:) ) 
(  409) #if defined key_diaeiv
(  410)          v_msf_eiv(:,:) = ptr_vjk( v_eiv(:,:,:) ) 
(  411) #endif
(  412)          IF( ln_subbas ) THEN
(  413)             v_msf_atl(:,:) = ptr_vjk( v_atl(:,:,:) ) 
(  414)             v_msf_ipc(:,:) = ptr_vjk( v_ipc(:,:,:) ) 
(  415)             ht_atl(:) = SUM(ptr_vjk( vt_atl(:,:,:)),2 )
(  416)             ht_pac(:) = SUM(ptr_vjk( vt_pac(:,:,:)),2 )
(  417)             ht_ind(:) = SUM(ptr_vjk( vt_ind(:,:,:)),2 )
(  418)             st_atl(:) = SUM(ptr_vjk( vs_atl(:,:,:)),2 )
(  419)             st_pac(:) = SUM(ptr_vjk( vs_pac(:,:,:)),2 )
(  420)             st_ind(:) = SUM(ptr_vjk( vs_ind(:,:,:)),2 )
(  421)          ENDIF
(  422) 
(  423)          ! poleward tracer transports: 
(  424)          ! overturning components:
(  425)          pht_ove(:) = SUM( v_msf_glo(:,:) * tn_jk(:,:), 2 )   ! SUM over jk
(  426)          pst_ove(:) = SUM( v_msf_glo(:,:) * sn_jk(:,:), 2 )   ! SUM over jk
(  427) #if defined key_diaeiv
(  428)          pht_eiv(:) = SUM( v_msf_eiv(:,:) * tn_jk(:,:), 2 )   ! SUM over jk
(  429)          pst_eiv(:) = SUM( v_msf_eiv(:,:) * sn_jk(:,:), 2 )   ! SUM over jk
(  430) #endif
(  431)       
(  432)          ! conversion in PW and G g
(  433)          zpwatt = zpwatt * rau0 * rcp
(  434)          pht_adv(:) = pht_adv(:) * zpwatt  
(  435)          pht_ove(:) = pht_ove(:) * zpwatt
(  436)          pht_ldf(:) = pht_ldf(:) * zpwatt
(  437)          pst_adv(:) = pst_adv(:) * zggram
(  438)          pst_ove(:) = pst_ove(:) * zggram
(  439)          pst_ldf(:) = pst_ldf(:) * zggram
(  440) #if defined key_diaeiv
(  441)          pht_eiv(:) = pht_eiv(:) * zpwatt
(  442)          pst_eiv(:) = pst_eiv(:) * zggram
(  443) #endif
(  444)          IF( ln_subbas ) THEN
(  445)             ht_atl(:) = ht_atl(:) * zpwatt
(  446)             ht_pac(:) = ht_pac(:) * zpwatt
(  447)             ht_ind(:) = ht_ind(:) * zpwatt
(  448)             st_atl(:) = st_atl(:) * zggram 
(  449)             st_pac(:) = st_pac(:) * zggram
(  450)             st_ind(:) = st_ind(:) * zggram
(  451)          ENDIF
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 9

(  453)          ! "Meridional" Stream-Function
(  454)          DO jk = 2,jpk 
(  455)             v_msf_glo(:,jk) = v_msf_glo(:,jk-1) + v_msf_glo(:,jk)
(  456)          END DO
(  457)          v_msf_glo(:,:) = v_msf_glo(:,:) * zsverdrup
(  458) 
(  459) #if defined key_diaeiv
(  460)          ! Bolus "Meridional" Stream-Function
(  461)          DO jk = 2,jpk 
(  462)             v_msf_eiv(:,jk) = v_msf_eiv(:,jk-1) + v_msf_eiv(:,jk)
(  463)          END DO
(  464)          v_msf_eiv(:,:) = v_msf_eiv(:,:) * zsverdrup
(  465) #endif
(  466) 
(  467)          IF( ln_subbas ) THEN
(  468)             DO jk = 2,jpk 
(  469)                v_msf_atl(:,jk) = v_msf_atl(:,jk-1) + v_msf_atl(:,jk)
(  470)                v_msf_ipc(:,jk) = v_msf_ipc(:,jk-1) + v_msf_ipc(:,jk)
(  471)             END DO
(  472)             v_msf_atl(:,:) = v_msf_atl(:,:) * zsverdrup
(  473)             v_msf_ipc(:,:) = v_msf_ipc(:,:) * zsverdrup
(  474)          ENDIF
(  475) 
(  476)          ! outputs
(  477)          CALL dia_ptr_wri( kt )
(  478) 
(  479)       ENDIF
(  480) 
(  481)       ! Close the file
(  482)       IF( kt == nitend ) CALL histclo( numptr )
(  483) 
(  484)    END SUBROUTINE dia_ptr
(  485) 
(  486) 
(  487)    SUBROUTINE dia_ptr_init
(  488)       !!----------------------------------------------------------------------
(  489)       !!                  ***  ROUTINE dia_ptr_init  ***
(  490)       !!                   
(  491)       !! ** Purpose :   Initialization, namelist read
(  492)       !!
(  493)       !! ** Method  :   
(  494)       !!
(  495)       !! ** input   :   Namlist namptr
(  496)       !!
(  497)       !! ** Action  :  
(  498)       !!
(  499)       !! history :
(  500)       !!   9.0  !  03-08  (Autor Names)  Original code
(  501)       !!----------------------------------------------------------------------
(  502)       !! * local declarations
(  503)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z_1         ! temporary workspace
(  504) 
(  505)       NAMELIST/namptr/ ln_diaptr, ln_subbas, nf_ptr
(  506)       !!----------------------------------------------------------------------
(  507) 
(  508)       ! Read Namelist namptr : poleward transport parameters
(  509)       REWIND ( numnam )
(  510)       READ   ( numnam, namptr )






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 10

(  511) 
(  512) 
(  513)       ! Control print
(  514)       IF(lwp) THEN
(  515)          WRITE(numout,*)
(  516)          WRITE(numout,*) 'dia_ptr_init : poleward transport and msf initialization'
(  517)          WRITE(numout,*) '~~~~~~~~~~~~'
(  518)          WRITE(numout,*) '          Namelist namptr : set ptr parameters'
(  519)          WRITE(numout,*) '             Switch for ptr diagnostic (T) or not (F) ln_diaptr = ', ln_diaptr
(  520)          WRITE(numout,*) '             Atla/Paci/Ind basins computation         ln_subbas = ', ln_subbas
(  521)          WRITE(numout,*) '             Frequency of computation                    nf_ptr = ', nf_ptr
(  522)       ENDIF
(  523) 
(  524)       ! inverse of the ocean "zonal" v-point section
(  525)       z_1(:,:,:) = 1.e0
(  526)       surf_jk_r(:,:) = ptr_vtjk( z_1(:,:,:) )
(  527)       WHERE( surf_jk_r(:,:) /= 0.e0 )   surf_jk_r(:,:) = 1.e0 / surf_jk_r(:,:)
(  528) 
(  529)    END SUBROUTINE dia_ptr_init
(  530) 
(  531)    !!---------------------------------------------------------------------
(  532)    !!   Default option :                                       NetCDF file
(  533)    !!---------------------------------------------------------------------
(  534) 
(  535)    SUBROUTINE dia_ptr_wri( kt )
(  536)       !!---------------------------------------------------------------------
(  537)       !!                ***  ROUTINE dia_ptr_wri  ***
(  538)       !!
(  539)       !! ** Purpose :   output of poleward fluxes
(  540)       !!
(  541)       !! ** Method  :   NetCDF file
(  542)       !!
(  543)       !! History :
(  544)       !!   9.0  !  03-09  (G. Madec)  Original code
(  545)       !!----------------------------------------------------------------------
(  546)       USE ioipsl          ! NetCDF IPSL library
(  547)       USE daymod
(  548) 
(  549)       !! * Arguments
(  550)       INTEGER, INTENT(in) ::   kt   ! ocean time-step index
(  551) 
(  552)       !! * Save variables   
(  553)       INTEGER, SAVE ::   nhoridz, ndepidzt, ndepidzw, ndex(1)
(  554) 
(  555)       !! * Local variables
(  556)       CHARACTER (len=40) ::   &
(  557)          clhstnam, clop             ! temporary names
(  558)       INTEGER ::   iline, it, ji    !
(  559)       REAL(wp) ::   &
(  560)          zsto, zout, zdt, zmax, &   ! temporary scalars
(  561)          zjulian
(  562)       REAL(wp), DIMENSION(jpj) ::   zphi, zfoo
(  563)       !!----------------------------------------------------------------------
(  564)       
(  565)       ! Define frequency of output and means
(  566)       zdt = rdt
(  567)       IF( nacc == 1 ) zdt = rdtmin
(  568) #if defined key_diainstant






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 11

(  569)          zsto = nf_ptr * zdt
(  570)          clop = "inst(x)"               ! no use of the mask value (require less cpu time)
(  571)          !!! clop="inst(only(x))"       ! put 1.e+20 on land (very expensive!!)
(  572) #else
(  573)          zsto = zdt
(  574)          clop = "ave(x)"                ! no use of the mask value (require less cpu time)
(  575)          !!! clop="ave(only(x))"        ! put 1.e+20 on land (very expensive!!)
(  576) #endif
(  577)       zout = nf_ptr * zdt
(  578)       zmax = ( nitend - nit000 + 1 ) * zdt
(  579)       
(  580)          
(  581)       ! define time axis
(  582)       it = kt - nit000 + 1
(  583) 
(  584)       ! Initialization
(  585)       ! --------------
(  586)       IF( kt == nit000 ) THEN
(  587)       
(  588)       zdt = rdt
(  589)       IF( nacc == 1 ) zdt = rdtmin
(  590) 
(  591)          ! Reference latitude
(  592)          ! ------------------
(  593)          !                                           ! =======================
(  594) !!DB: Orca-related
(  595) !         IF( cp_cfg == "orca" ) THEN                 !   ORCA configurations
(  596) !            !                                        ! =======================!
(  597) !
(  598) !            IF( jp_cfg == 05  )   iline = 192   ! i-line that passes near the North Pole
(  599) !            IF( jp_cfg == 025 )   iline = 384   ! i-line that passes near the North Pole
(  600) !            IF( jp_cfg == 2   )   iline =  48   ! i-line that passes near the North Pole
(  601) !            IF( jp_cfg == 4   )   iline =  24   ! i-line that passes near the North Pole
(  602) !            zphi(:) = 0.e0
(  603) !            DO ji = mi0(iline), mi1(iline) 
(  604) !               zphi(:) = gphiv(ji,:)         ! if iline is in the local domain
(  605) !               ! correct highest latitude for ORCA05
(  606) !               IF( jp_cfg == 05  ) zphi(jpj) = zphi(jpjm1) + (zphi(jpjm1)-zphi(jpj-2))/2.
(  607) !               IF( jp_cfg == 05  ) zphi(jpj) = MIN( zphi(jpj), 90.)!!!
(  608) !
(  609) !            END DO
(  610) !            ! provide the correct zphi to all local domains
(  611) !            IF( lk_mpp )   CALL mpp_sum( zphi, jpj )        !
(  612) !
(  613) !            !                                        ! =======================
(  614) !         ELSE                                        !   OTHER configurations!
(  615)             !                                        ! =======================!
(  616)             zphi(:) = gphiv(1,:)             ! assume lat/lon coordinate, select the first i-line
(  617)             !
(  618) !         ENDIF
(  619) 
(  620)          ! OPEN netcdf file 
(  621)          ! ----------------
(  622)          ! Define frequency of output and means
(  623)          zsto = nf_ptr * zdt
(  624)          clop = "ave(x)"
(  625)          zout = nf_ptr * zdt
(  626)          zfoo(:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 12

(  627) 
(  628)          ! Compute julian date from starting date of the run
(  629) 
(  630)          CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
(  631) 
(  632)          CALL dia_nam( clhstnam, nf_ptr, 'diaptr' )
(  633)          IF(lwp)WRITE( numout,*)" Name of diaptr NETCDF file ",clhstnam
(  634) 
(  635)          ! Horizontal grid : zphi()
(  636)          CALL histbeg(clhstnam, 1, zfoo, jpj, zphi,   &
(  637)             1, 1, 1, jpj, 0, zjulian, zdt, nhoridz, numptr, domain_id=nidom )
(  638)          ! Vertical grids : gdept, gdepw
(  639)          CALL histvert( numptr, "deptht", "Vertical T levels",   &
(  640)             "m", jpk, gdept, ndepidzt )
(  641)          CALL histvert( numptr, "depthw", "Vertical W levels",   &
(  642)             "m", jpk, gdepw, ndepidzw )
(  643)          
(  644)          !  Zonal mean T and S
(  645)          
(  646)          CALL histdef( numptr, "zotemglo", "Zonal Mean Temperature","C" ,   &
(  647)             1, jpj, nhoridz, jpk, 1, jpk, ndepidzt, 32, clop, zsto, zout )
(  648)          CALL histdef( numptr, "zosalglo", "Zonal Mean Salinity","PSU"  ,   &
(  649)             1, jpj, nhoridz, jpk, 1, jpk, ndepidzt, 32, clop, zsto, zout )
(  650) 
(  651)          !  Meridional Stream-Function (eulerian and bolus)
(  652)          
(  653)          CALL histdef( numptr, "zomsfglo", "Meridional Stream-Function: Global","Sv" ,   &
(  654)             1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
(  655)          IF( ln_subbas ) THEN
(  656)             CALL histdef( numptr, "zomsfatl", "Meridional Stream-Function: Atlantic","Sv" ,   &
(  657)                1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
(  658)             CALL histdef( numptr, "zomsfipc", "Meridional Stream-Function: Indo-Pacific","Sv" ,&
(  659)                1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
(  660)          ENDIF
(  661) 
(  662)          !  Heat transport 
(  663) 
(  664)          CALL histdef( numptr, "sophtadv", "Advective Heat Transport"      ,   &
(  665)             "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  666)          CALL histdef( numptr, "sophtldf", "Diffusive Heat Transport"      ,   &
(  667)             "PW",1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  668)          CALL histdef( numptr, "sophtove", "Overturning Heat Transport"    ,   &
(  669)             "PW",1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  670)          IF( ln_subbas ) THEN
(  671)             CALL histdef( numptr, "sohtatl", "Heat Transport Atlantic"      ,  &
(  672)                "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  673)             CALL histdef( numptr, "sohtpac", "Heat Transport Pacific"      ,   &
(  674)                "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  675)             CALL histdef( numptr, "sohtind", "Heat Transport Indic"      ,     &
(  676)                "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  677)          ENDIF
(  678) 
(  679) 
(  680)          !  Salt transport 
(  681) 
(  682)          CALL histdef( numptr, "sopstadv", "Advective Salt Transport"      ,   &
(  683)             "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  684)          CALL histdef( numptr, "sopstldf", "Diffusive Salt Transport"      ,   &






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 13

(  685)             "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  686)          CALL histdef( numptr, "sopstove", "Overturning Salt Transport"    ,   &
(  687)             "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  688) 
(  689) #if defined key_diaeiv
(  690)          ! Eddy induced velocity
(  691)          CALL histdef( numptr, "zomsfeiv", "Bolus Meridional Stream-Function: global",   &
(  692)             "Sv"      , 1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
(  693)          CALL histdef( numptr, "sophteiv", "Bolus Advective Heat Transport",   &
(  694)             "PW"      , 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  695)          CALL histdef( numptr, "sopsteiv", "Bolus Advective Salt Transport",   &
(  696)             "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  697) #endif
(  698)          IF( ln_subbas ) THEN
(  699)             CALL histdef( numptr, "sostatl", "Salt Transport Atlantic"      ,    &
(  700)                "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  701)             CALL histdef( numptr, "sostpac", "Salt Transport Pacific"      ,     &
(  702)                "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  703)             CALL histdef( numptr, "sostind", "Salt Transport Indic"      ,       &
(  704)                "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  705)          ENDIF
(  706)          
(  707) 
(  708)          CALL histend( numptr )
(  709) 
(  710)       ENDIF
(  711) 
(  712)       IF( MOD( kt, nf_ptr ) == 0 ) THEN
(  713) 
(  714)          IF(lwp) THEN
(  715)             WRITE(numout,*)
(  716)             WRITE(numout,*) 'dia_ptr : write Poleward Transports at time-step : ', kt
(  717)             WRITE(numout,*) '~~~~~~~~'
(  718)             WRITE(numout,*)
(  719)          ENDIF
(  720) 
(  721)          ! define time axis
(  722)          it= kt - nit000 + 1
(  723)          ndex(1) = 0
(  724)          CALL histwrite( numptr, "zotemglo", it, tn_jk    , jpj*jpk, ndex )
(  725)          CALL histwrite( numptr, "zosalglo", it, sn_jk    , jpj*jpk, ndex )
(  726)          ! overturning outputs:
(  727)          CALL histwrite( numptr, "zomsfglo", it, v_msf_glo , jpj*jpk, ndex )
(  728)          IF( ln_subbas ) THEN
(  729)             CALL histwrite( numptr, "zomsfatl", it, v_msf_atl , jpj*jpk, ndex )
(  730)             CALL histwrite( numptr, "zomsfipc", it, v_msf_ipc , jpj*jpk, ndex )
(  731)          ENDIF
(  732)          ! heat transport outputs:
(  733)          IF( ln_subbas ) THEN
(  734)             CALL histwrite( numptr, "sohtatl", it, ht_atl  , jpj, ndex )
(  735)             CALL histwrite( numptr, "sohtpac", it, ht_pac  , jpj, ndex )
(  736)             CALL histwrite( numptr, "sohtind", it, ht_ind  , jpj, ndex )
(  737)             CALL histwrite( numptr, "sostatl", it, st_atl  , jpj, ndex )
(  738)             CALL histwrite( numptr, "sostpac", it, st_pac  , jpj, ndex )
(  739)             CALL histwrite( numptr, "sostind", it, st_ind  , jpj, ndex )
(  740)          ENDIF
(  741) 
(  742)          CALL histwrite( numptr, "sophtadv", it, pht_adv  , jpj, ndex )






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 14

(  743)          CALL histwrite( numptr, "sophtldf", it, pht_ldf  , jpj, ndex )
(  744)          CALL histwrite( numptr, "sophtove", it, pht_ove  , jpj, ndex )
(  745)          CALL histwrite( numptr, "sopstadv", it, pst_adv  , jpj, ndex )
(  746)          CALL histwrite( numptr, "sopstldf", it, pst_ldf  , jpj, ndex )
(  747)          CALL histwrite( numptr, "sopstove", it, pst_ove  , jpj, ndex )
(  748) #if defined key_diaeiv
(  749)          CALL histwrite( numptr, "zomsfeiv", it, v_msf_eiv, jpj*jpk, ndex )
(  750)          CALL histwrite( numptr, "sophteiv", it, pht_eiv  , jpj    , ndex )
(  751)          CALL histwrite( numptr, "sopsteiv", it, pst_eiv  , jpj    , ndex )
(  752) #endif
(  753)  
(  754)       ENDIF
(  755) 
(  756)    END SUBROUTINE dia_ptr_wri
(  757) 
(  758)    !!======================================================================
(  759) END MODULE diaptr















































PGF90 (Version     10.5)          09/23/2010  09:49:28      page 15

(    1) # 1 "diaptr.F90"
(    1) MODULE diaptr
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  diaptr  ***
(    4)    !! Ocean physics:  brief description of the purpose of the module
(    5)    !!                 (please no more than 2 lines)
(    6)    !!=====================================================================
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dia_ptr      : Poleward Transport Diagnostics module
(    9)    !!   dia_ptr_init : Initialization, namelist read
(   10)    !!   dia_ptr_wri  : Output of poleward fluxes
(   11)    !!   ptr_vjk      : "zonal" sum computation of a "meridional" flux array
(   12)    !!   ptr_vtjk     : "zonal" mean computation of a tracer field
(   13)    !!   ptr_vj       : "zonal" and vertical sum computation of a "meridional"
(   14)    !!                : flux array; Generic interface: ptr_vj_3d, ptr_vj_2d
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce           ! ocean dynamics and active tracers
(   18)    USE dom_oce       ! ocean space and time domain
(   19)    USE ldftra_oce    ! ???
(   20)    USE lib_mpp
(   21)    USE in_out_manager
(   22)    USE dianam
(   23)    USE phycst
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    INTERFACE ptr_vj
(   29)       MODULE PROCEDURE ptr_vj_3d, ptr_vj_2d
(   30)    END INTERFACE
(   31) 
(   32)    !! *  Routine accessibility
(   33)    PUBLIC dia_ptr_init   ! call in opa module
(   34)    PUBLIC dia_ptr        ! call in step module
(   35)    PUBLIC ptr_vj         ! call by tra_ldf & tra_adv routines
(   36)    PUBLIC ptr_vjk        ! call by tra_ldf & tra_adv routines
(   37) 
(   38)    !! * Share Module variables
(   39)    LOGICAL, PUBLIC ::       & !!! ** init namelist (namptr) **
(   40)       ln_diaptr = .FALSE.,  &  !: Poleward transport flag (T) or not (F)
(   41)       ln_subbas = .FALSE.      !: Atlantic/Pacific/Indian basins calculation
(   42)    INTEGER, PUBLIC ::       & !!: ** ptr namelist (namptr) **
(   43)       nf_ptr = 15              !: frequency of ptr computation
(   44)    REAL(wp), PUBLIC, DIMENSION(jpj) ::   &   !!: poleward transport
(   45)       pht_adv, pst_adv,     &  !: heat and salt: advection
(   46)       pht_ove, pst_ove,     &  !: heat and salt: overturning
(   47)       pht_ldf, pst_ldf,     &  !: heat and salt: lateral diffusion
(   51)       ht_atl,ht_ind,ht_pac, &  !: heat
(   52)       st_atl,st_ind,st_pac     !: salt
(   53)    REAL(wp),DIMENSION(jpi,jpj) ::   &
(   54)       abasin,pbasin,ibasin     !: return function value
(   55)      
(   56) 
(   57)    !! Module variables
(   58)    REAL(wp), DIMENSION(jpj,jpk) ::   &  
(   59)       tn_jk  , sn_jk  ,  &  !: "zonal" mean temperature and salinity
(   60)       v_msf_atl       ,  &  !: "meridional" Stream-Function






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 16

(   61)       v_msf_glo       ,  &  !: "meridional" Stream-Function
(   62)       v_msf_ipc       ,  &  !: "meridional" Stream-Function
(   66)       surf_jk_r             !: inverse of the ocean "zonal" section surface
(   67) 
(   68)    !! * Substitutions
(   69) # 1 "./domzgr_substitute.h90"
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
(   48) # 70 "diaptr.F90"
(   70) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 71 "diaptr.F90"
(   71) # 71
(   71)    !!----------------------------------------------------------------------
(   72)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   73)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diaptr.F90,v 1.10 2006/03/20 16:05:53 opalod Exp $ 
(   74)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   75)    !!----------------------------------------------------------------------
(   76) 
(   77) CONTAINS
(   78) 
(   79)    FUNCTION ptr_vj_3d( pva )   RESULT ( p_fval )
(   80)       !!----------------------------------------------------------------------
(   81)       !!                    ***  ROUTINE ptr_vj_3d  ***
(   82)       !!
(   83)       !! ** Purpose :   "zonal" and vertical sum computation of a "meridional"
(   84)       !!      flux array
(   85)       !!
(   86)       !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
(   87)       !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
(   88)       !!
(   89)       !! ** Action  : - p_fval: i-k-mean poleward flux of pva
(   90)       !!
(   91)       !! History :
(   92)       !!   9.0  !  03-09  (G. Madec)  Original code
(   93)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 17

(   94)       !! * arguments
(   95)       REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
(   96)          pva                         ! mask flux array at V-point
(   97) 
(   98)       !! * local declarations
(   99)       INTEGER  ::   ji, jj, jk        ! dummy loop arguments
(  100) # 101
(  101)       INTEGER  ::   ijpj = jpj        ! ???
(  102) # 105
(  105)       REAL(wp),DIMENSION(jpj) ::   &
(  106)          p_fval                       ! function value
(  107)       !!--------------------------------------------------------------------
(  108) 
(  109) # 112
(  112)       p_fval(:) = 0.e0
(  113)       DO jk = 1, jpkm1
(  114)          DO jj = 2, jpjm1
(  115)             DO ji = 2, jpim1   ! Vector opt.
(  116)                p_fval(jj) = p_fval(jj) + pva(ji,jj,jk) * tmask_i(ji,jj+1) * tmask_i(ji,jj) 
(  117)             END DO
(  118)          END DO
(  119)       END DO
(  120) 
(  121)       IF( lk_mpp )   CALL mpp_sum( p_fval, ijpj )     !!bug  I presume
(  122) 
(  123)    END FUNCTION ptr_vj_3d
(  124) 
(  125) 
(  126) 
(  127)    FUNCTION ptr_vj_2d( pva )   RESULT ( p_fval )
(  128)       !!----------------------------------------------------------------------
(  129)       !!                    ***  ROUTINE ptr_vj_2d  ***
(  130)       !!
(  131)       !! ** Purpose :   "zonal" and vertical sum computation of a "meridional"
(  132)       !!      flux array
(  133)       !!
(  134)       !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
(  135)       !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
(  136)       !!
(  137)       !! ** Action  : - p_fval: i-k-mean poleward flux of pva
(  138)       !!
(  139)       !! History :
(  140)       !!   9.0  !  03-09  (G. Madec)  Original code
(  141)       !!----------------------------------------------------------------------
(  142)       !! * arguments
(  143)       REAL(wp) , INTENT(in), DIMENSION(jpi,jpj) ::   &
(  144)          pva                         ! mask flux array at V-point
(  145) 
(  146)       !! * local declarations
(  147)       INTEGER  ::   ji,jj             ! dummy loop arguments
(  148) # 149
(  149)       INTEGER  ::   ijpj = jpj        ! ???
(  150) # 153
(  153)       REAL(wp),DIMENSION(jpj) ::   &
(  154)          p_fval                       ! function value
(  155)       !!--------------------------------------------------------------------
(  156) # 159
(  159)  






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 18

(  160)       p_fval(:) = 0.e0
(  161)       DO jj = 2, jpjm1
(  162)          DO ji = 2, jpim1   ! Vector opt.
(  163)             p_fval(jj) = p_fval(jj) + pva(ji,jj) * tmask_i(ji,jj+1) * tmask_i(ji,jj)
(  164)          END DO
(  165)       END DO
(  166) 
(  167)       IF( lk_mpp )   CALL mpp_sum( p_fval, ijpj )     !!bug  I presume
(  168)  
(  169)     END FUNCTION ptr_vj_2d
(  170) 
(  171) 
(  172) 
(  173)    FUNCTION ptr_vjk( pva )   RESULT ( p_fval )
(  174)       !!----------------------------------------------------------------------
(  175)       !!                    ***  ROUTINE ptr_vjk  ***
(  176)       !!
(  177)       !! ** Purpose :   "zonal" sum computation of a "meridional" flux array
(  178)       !!
(  179)       !! ** Method  : - i-sum of pva using the interior 2D vmask (vmask_i).
(  180)       !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
(  181)       !!
(  182)       !! ** Action  : - p_fval: i-k-mean poleward flux of pva
(  183)       !!
(  184)       !! History :
(  185)       !!   9.0  !  03-09  (G. Madec)  Original code
(  186)       !!----------------------------------------------------------------------
(  187)       !! * arguments
(  188)       REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
(  189)          pva                         ! mask flux array at V-point
(  190) 
(  191)       !! * local declarations
(  192)       INTEGER  ::   ji, jj, jk        ! dummy loop arguments
(  193)       INTEGER, DIMENSION (1) :: ish
(  194)       INTEGER, DIMENSION (2) :: ish2
(  195)       REAL(wp),DIMENSION(jpj*jpk) ::   &
(  196)          zwork                        ! temporary vector for mpp_sum
(  197)       REAL(wp),DIMENSION(jpj,jpk) ::   &
(  198)          p_fval                       ! return function value
(  199)       !!--------------------------------------------------------------------
(  200)  
(  201)       p_fval(:,:) = 0.e0
(  202) 
(  203)       DO jk = 1, jpkm1
(  204)          DO jj = 2, jpjm1
(  205)            DO ji = 2, jpim1
(  206)             p_fval(jj,jk) = p_fval(jj,jk) + pva(ji,jj,jk) * e1v(ji,jj) * e3t(jk) &
(  207)                &            * tmask_i(ji,jj+1) * tmask_i(ji,jj)
(  208)            END DO
(  209)          END DO
(  210)       END DO
(  211) 
(  212)       IF(lk_mpp)   THEN
(  213)          ish(1) = jpj*jpk ; ish2(1)=jpj ; ish2(2)=jpk
(  214)          zwork(:)= RESHAPE(p_fval, ish )
(  215)          CALL mpp_sum(zwork, jpj*jpk )
(  216)          p_fval(:,:)= RESHAPE(zwork,ish2)
(  217)       END IF






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 19

(  218) 
(  219)    END FUNCTION ptr_vjk
(  220) 
(  221)    FUNCTION ptr_vtjk( pva )   RESULT ( p_fval )
(  222)       !!----------------------------------------------------------------------
(  223)       !!                    ***  ROUTINE ptr_vtjk  ***
(  224)       !!
(  225)       !! ** Purpose :   "zonal" mean computation of a tracer field
(  226)       !!
(  227)       !! ** Method  : - i-sum of mj(pva) using the interior 2D vmask (vmask_i)
(  228)       !!      multiplied by the inverse of the surface of the "zonal" ocean
(  229)       !!      section
(  230)       !!
(  231)       !! ** Action  : - p_fval: i-k-mean poleward flux of pva
(  232)       !!
(  233)       !! History :
(  234)       !!   9.0  !  03-09  (G. Madec)  Original code
(  235)       !!   9.0  !  06-01  (A. Biastoch)  Allow sub-basins computation
(  236)       !!----------------------------------------------------------------------
(  237)       !! * arguments
(  238)       REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
(  239)          pva                         ! mask flux array at V-point
(  240)  
(  241)       !! * local declarations
(  242)       INTEGER  ::   ji, jj, jk        ! dummy loop arguments
(  243)       INTEGER, DIMENSION (1) :: ish
(  244)       INTEGER, DIMENSION (2) :: ish2
(  245)       REAL(wp),DIMENSION(jpj*jpk) ::   &
(  246)          zwork                        ! temporary vector for mpp_sum
(  247)       REAL(wp),DIMENSION(jpj,jpk) ::   &
(  248)          p_fval                       ! return function value
(  249)       !!-------------------------------------------------------------------- 
(  250) 
(  251)       p_fval(:,:) = 0.e0
(  252)       DO jk = 1, jpkm1
(  253)          DO jj = 2, jpjm1
(  254)             DO ji = 2, jpim1   ! Vector opt.
(  255)                p_fval(jj,jk) = p_fval(jj,jk) + ( pva(ji,jj,jk) + pva(ji,jj+1,jk) )              &
(  256)                   &                          * e1v(ji,jj) * e3t(jk) * vmask(ji,jj,jk)   &
(  257)                   &                          * tmask_i(ji,jj+1) * tmask_i(ji,jj)
(  258)             END DO
(  259)          END DO
(  260)       END DO
(  261)       p_fval(:,:) = p_fval(:,:) * 0.5
(  262)       IF(lk_mpp)   THEN
(  263)          ish(1) = jpj*jpk ; ish2(1)=jpj ; ish2(2)=jpk
(  264)          zwork(:)= RESHAPE(p_fval, ish )
(  265)          CALL mpp_sum(zwork, jpj*jpk )
(  266)          p_fval(:,:)= RESHAPE(zwork,ish2)
(  267)       END IF
(  268) 
(  269)    END FUNCTION ptr_vtjk
(  270) 
(  271) 
(  272)    SUBROUTINE dia_ptr( kt )
(  273)       !!----------------------------------------------------------------------
(  274)       !!                  ***  ROUTINE dia_ptr  ***
(  275)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 20

(  276)       !! * Moudules used
(  277)       USE ioipsl
(  278) 
(  279)       !! * Argument
(  280)       INTEGER, INTENT(in) ::   kt   ! ocean time step index
(  281) 
(  282)       !! * Local variables
(  283)       INTEGER ::   jk,jj,ji               ! dummy loop
(  284)       REAL(wp) ::    &
(  285)          zsverdrup,  &              ! conversion from m3/s to Sverdrup
(  286)          zpwatt,     &              ! conversion from W    to PW
(  287)          zggram                     ! conversion from g    to Pg
(  288) 
(  289)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &
(  290)          v_atl , v_ipc,                    &
(  291)          vt_atl, vt_pac, vt_ind,           &
(  292)          vs_atl, vs_pac, vs_ind,           &
(  293)          zv_eiv
(  294)       CHARACTER (len=32) ::   &
(  295)          clnam = 'subbasins.nc'                
(  296)       INTEGER ::  itime,inum,ipi,ipj,ipk       ! temporary integer
(  297)       INTEGER, DIMENSION (1) ::   istep
(  298)       REAL(wp) ::    zdate0,zsecond,zdt        ! temporary scalars
(  299)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(  300)          zlamt, zphit, zdta             ! temporary workspace (NetCDF read)
(  301)       REAL(wp), DIMENSION(jpk) ::   &
(  302)          zdept                          ! temporary workspace (NetCDF read)
(  303)       !!----------------------------------------------------------------------
(  304) 
(  305)       IF( kt == nit000 .OR. MOD( kt, nf_ptr ) == 0 )   THEN
(  306) 
(  307)          zsverdrup = 1.e-6
(  308)          zpwatt    = 1.e-15
(  309)          zggram    = 1.e-6
(  310)          ipi       = jpidta
(  311)          ipj       = jpjdta
(  312)          ipk       = 1
(  313)          itime     = 1
(  314)          zsecond   = 0.e0
(  315)          zdate0    = 0.e0
(  316)    
(  317) # 320
(  320)          zv_eiv(:,:,:) = 0.e0
(  321) 
(  322) # 323
(  323)          ! "zonal" mean temperature and salinity at V-points
(  324)          tn_jk(:,:) = ptr_vtjk( tn(:,:,:) ) * surf_jk_r(:,:)
(  325)          sn_jk(:,:) = ptr_vtjk( sn(:,:,:) ) * surf_jk_r(:,:)
(  326) 
(  327)          !--------------------------------------------------------
(  328)          ! overturning calculation:
(  329)  
(  330)          IF( ln_subbas ) THEN              ! Basins computation
(  331) 
(  332)             IF( kt == nit000 ) THEN                ! load basin mask
(  333)                itime = 1
(  334)                ipi   = jpidta
(  335)                ipj   = jpjdta






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 21

(  336)                ipk   = 1
(  337)                zdt   = 0.e0
(  338)                istep = 0
(  339)                clnam = 'subbasins.nc'
(  340) 
(  341)                CALL flinopen(clnam,1,jpidta,1,jpjdta,.FALSE.,ipi,ipj, &
(  342)                   &          ipk,zlamt,zphit,zdept,itime,istep,zdate0,zdt,inum)
(  343) 
(  344)                ! get basins:
(  345)                abasin (:,:) = 0.e0
(  346)                pbasin (:,:) = 0.e0
(  347)                ibasin (:,:) = 0.e0
(  348) 
(  349)                ! Atlantic basin
(  350)                CALL flinget(inum,'atlmsk',jpidta,jpjdta,1,itime,1,   &
(  351)                   &         0,1,jpidta,1,jpjdta,zdta(:,:))
(  352)                DO jj = 1, nlcj                                 ! interior values
(  353)                   DO ji = 1, nlci
(  354)                      abasin (ji,jj) = zdta( mig(ji), mjg(jj) )
(  355)                   END DO
(  356)                END DO
(  357) 
(  358)                ! Pacific basin
(  359)                CALL flinget(inum,'pacmsk',jpidta,jpjdta,1,itime,1,   &
(  360)                   &         0,1,jpidta,1,jpjdta,zdta(:,:))
(  361)                DO jj = 1, nlcj                                 ! interior values
(  362)                   DO ji = 1, nlci
(  363)                      pbasin (ji,jj) = zdta( mig(ji), mjg(jj) )
(  364)                   END DO
(  365)                END DO
(  366) 
(  367)                ! Indian basin
(  368)                CALL flinget(inum,'indmsk',jpidta,jpjdta,1,itime,1,   &
(  369)                   &         0,1,jpidta,1,jpjdta,zdta(:,:))
(  370)                DO jj = 1, nlcj                                 ! interior values
(  371)                   DO ji = 1, nlci
(  372)                      ibasin (ji,jj) = zdta( mig(ji), mjg(jj) )
(  373)                   END DO
(  374)                END DO
(  375) 
(  376)                CALL flinclo(inum)
(  377) 
(  378)             ENDIF
(  379) 
(  380)             ! basin separation:
(  381)             DO jj = 1, jpj
(  382)                DO ji = 1, jpi
(  383)                   ! basin separated velocity
(  384)                   v_atl(ji,jj,:) = (vn(ji,jj,:)+zv_eiv(ji,jj,:))*abasin(ji,jj)
(  385)                   v_ipc(ji,jj,:) = (vn(ji,jj,:)+zv_eiv(ji,jj,:))*(pbasin(ji,jj)+ibasin(ji,jj))
(  386) 
(  387)                   ! basin separated T times V on T points
(  388)                   vt_ind(ji,jj,:) = tn(ji,jj,:) *                                 &
(  389)                      &              ( (vn    (ji,jj,:) + vn    (ji,jj-1,:))*0.5   &
(  390)                      &              + (zv_eiv(ji,jj,:) + zv_eiv(ji,jj-1,:))*0.5 ) 
(  391)                   vt_atl(ji,jj,:) = vt_ind(ji,jj,:) * abasin(ji,jj)
(  392)                   vt_pac(ji,jj,:) = vt_ind(ji,jj,:) * pbasin(ji,jj)
(  393)                   vt_ind(ji,jj,:) = vt_ind(ji,jj,:) * ibasin(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 22

(  394) 
(  395)                   ! basin separated S times V on T points
(  396)                   vs_ind(ji,jj,:) = sn(ji,jj,:) *                                 &
(  397)                      &              ( (vn    (ji,jj,:) + vn    (ji,jj-1,:))*0.5   &
(  398)                      &              + (zv_eiv(ji,jj,:) + zv_eiv(ji,jj-1,:))*0.5 ) 
(  399)                   vs_atl(ji,jj,:) = vs_ind(ji,jj,:) * abasin(ji,jj)
(  400)                   vs_pac(ji,jj,:) = vs_ind(ji,jj,:) * pbasin(ji,jj)
(  401)                   vs_ind(ji,jj,:) = vs_ind(ji,jj,:) * ibasin(ji,jj)
(  402)                END DO
(  403)             END DO
(  404) 
(  405)          ENDIF
(  406) 
(  407)          ! horizontal integral and vertical dz 
(  408)          v_msf_glo(:,:) = ptr_vjk( vn(:,:,:) ) 
(  409) # 412
(  412)          IF( ln_subbas ) THEN
(  413)             v_msf_atl(:,:) = ptr_vjk( v_atl(:,:,:) ) 
(  414)             v_msf_ipc(:,:) = ptr_vjk( v_ipc(:,:,:) ) 
(  415)             ht_atl(:) = SUM(ptr_vjk( vt_atl(:,:,:)),2 )
(  416)             ht_pac(:) = SUM(ptr_vjk( vt_pac(:,:,:)),2 )
(  417)             ht_ind(:) = SUM(ptr_vjk( vt_ind(:,:,:)),2 )
(  418)             st_atl(:) = SUM(ptr_vjk( vs_atl(:,:,:)),2 )
(  419)             st_pac(:) = SUM(ptr_vjk( vs_pac(:,:,:)),2 )
(  420)             st_ind(:) = SUM(ptr_vjk( vs_ind(:,:,:)),2 )
(  421)          ENDIF
(  422) 
(  423)          ! poleward tracer transports: 
(  424)          ! overturning components:
(  425)          pht_ove(:) = SUM( v_msf_glo(:,:) * tn_jk(:,:), 2 )   ! SUM over jk
(  426)          pst_ove(:) = SUM( v_msf_glo(:,:) * sn_jk(:,:), 2 )   ! SUM over jk
(  427) # 431
(  431)       
(  432)          ! conversion in PW and G g
(  433)          zpwatt = zpwatt * rau0 * rcp
(  434)          pht_adv(:) = pht_adv(:) * zpwatt  
(  435)          pht_ove(:) = pht_ove(:) * zpwatt
(  436)          pht_ldf(:) = pht_ldf(:) * zpwatt
(  437)          pst_adv(:) = pst_adv(:) * zggram
(  438)          pst_ove(:) = pst_ove(:) * zggram
(  439)          pst_ldf(:) = pst_ldf(:) * zggram
(  440) # 444
(  444)          IF( ln_subbas ) THEN
(  445)             ht_atl(:) = ht_atl(:) * zpwatt
(  446)             ht_pac(:) = ht_pac(:) * zpwatt
(  447)             ht_ind(:) = ht_ind(:) * zpwatt
(  448)             st_atl(:) = st_atl(:) * zggram 
(  449)             st_pac(:) = st_pac(:) * zggram
(  450)             st_ind(:) = st_ind(:) * zggram
(  451)          ENDIF
(  452) 
(  453)          ! "Meridional" Stream-Function
(  454)          DO jk = 2,jpk 
(  455)             v_msf_glo(:,jk) = v_msf_glo(:,jk-1) + v_msf_glo(:,jk)
(  456)          END DO
(  457)          v_msf_glo(:,:) = v_msf_glo(:,:) * zsverdrup
(  458) 
(  459) 






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 23

(  460) # 467
(  467)          IF( ln_subbas ) THEN
(  468)             DO jk = 2,jpk 
(  469)                v_msf_atl(:,jk) = v_msf_atl(:,jk-1) + v_msf_atl(:,jk)
(  470)                v_msf_ipc(:,jk) = v_msf_ipc(:,jk-1) + v_msf_ipc(:,jk)
(  471)             END DO
(  472)             v_msf_atl(:,:) = v_msf_atl(:,:) * zsverdrup
(  473)             v_msf_ipc(:,:) = v_msf_ipc(:,:) * zsverdrup
(  474)          ENDIF
(  475) 
(  476)          ! outputs
(  477)          CALL dia_ptr_wri( kt )
(  478) 
(  479)       ENDIF
(  480) 
(  481)       ! Close the file
(  482)       IF( kt == nitend ) CALL histclo( numptr )
(  483) 
(  484)    END SUBROUTINE dia_ptr
(  485) 
(  486) 
(  487)    SUBROUTINE dia_ptr_init
(  488)       !!----------------------------------------------------------------------
(  489)       !!                  ***  ROUTINE dia_ptr_init  ***
(  490)       !!                   
(  491)       !! ** Purpose :   Initialization, namelist read
(  492)       !!
(  493)       !! ** Method  :   
(  494)       !!
(  495)       !! ** input   :   Namlist namptr
(  496)       !!
(  497)       !! ** Action  :  
(  498)       !!
(  499)       !! history :
(  500)       !!   9.0  !  03-08  (Autor Names)  Original code
(  501)       !!----------------------------------------------------------------------
(  502)       !! * local declarations
(  503)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z_1         ! temporary workspace
(  504) 
(  505)       NAMELIST/namptr/ ln_diaptr, ln_subbas, nf_ptr
(  506)       !!----------------------------------------------------------------------
(  507) 
(  508)       ! Read Namelist namptr : poleward transport parameters
(  509)       REWIND ( numnam )
(  510)       READ   ( numnam, namptr )
(  511) 
(  512) 
(  513)       ! Control print
(  514)       IF(lwp) THEN
(  515)          WRITE(numout,*)
(  516)          WRITE(numout,*) 'dia_ptr_init : poleward transport and msf initialization'
(  517)          WRITE(numout,*) '~~~~~~~~~~~~'
(  518)          WRITE(numout,*) '          Namelist namptr : set ptr parameters'
(  519)          WRITE(numout,*) '             Switch for ptr diagnostic (T) or not (F) ln_diaptr = ', ln_diaptr
(  520)          WRITE(numout,*) '             Atla/Paci/Ind basins computation         ln_subbas = ', ln_subbas
(  521)          WRITE(numout,*) '             Frequency of computation                    nf_ptr = ', nf_ptr
(  522)       ENDIF
(  523) 






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 24

(  524)       ! inverse of the ocean "zonal" v-point section
(  525)       z_1(:,:,:) = 1.e0
(  526)       surf_jk_r(:,:) = ptr_vtjk( z_1(:,:,:) )
(  527)       WHERE( surf_jk_r(:,:) /= 0.e0 )   surf_jk_r(:,:) = 1.e0 / surf_jk_r(:,:)
(  528) 
(  529)    END SUBROUTINE dia_ptr_init
(  530) 
(  531)    !!---------------------------------------------------------------------
(  532)    !!   Default option :                                       NetCDF file
(  533)    !!---------------------------------------------------------------------
(  534) 
(  535)    SUBROUTINE dia_ptr_wri( kt )
(  536)       !!---------------------------------------------------------------------
(  537)       !!                ***  ROUTINE dia_ptr_wri  ***
(  538)       !!
(  539)       !! ** Purpose :   output of poleward fluxes
(  540)       !!
(  541)       !! ** Method  :   NetCDF file
(  542)       !!
(  543)       !! History :
(  544)       !!   9.0  !  03-09  (G. Madec)  Original code
(  545)       !!----------------------------------------------------------------------
(  546)       USE ioipsl          ! NetCDF IPSL library
(  547)       USE daymod
(  548) 
(  549)       !! * Arguments
(  550)       INTEGER, INTENT(in) ::   kt   ! ocean time-step index
(  551) 
(  552)       !! * Save variables   
(  553)       INTEGER, SAVE ::   nhoridz, ndepidzt, ndepidzw, ndex(1)
(  554) 
(  555)       !! * Local variables
(  556)       CHARACTER (len=40) ::   &
(  557)          clhstnam, clop             ! temporary names
(  558)       INTEGER ::   iline, it, ji    !
(  559)       REAL(wp) ::   &
(  560)          zsto, zout, zdt, zmax, &   ! temporary scalars
(  561)          zjulian
(  562)       REAL(wp), DIMENSION(jpj) ::   zphi, zfoo
(  563)       !!----------------------------------------------------------------------
(  564)       
(  565)       ! Define frequency of output and means
(  566)       zdt = rdt
(  567)       IF( nacc == 1 ) zdt = rdtmin
(  568) # 573
(  573)          zsto = zdt
(  574)          clop = "ave(x)"                ! no use of the mask value (require less cpu time)
(  575)          !!! clop="ave(only(x))"        ! put 1.e+20 on land (very expensive!!)
(  576) # 577
(  577)       zout = nf_ptr * zdt
(  578)       zmax = ( nitend - nit000 + 1 ) * zdt
(  579)       
(  580)          
(  581)       ! define time axis
(  582)       it = kt - nit000 + 1
(  583) 
(  584)       ! Initialization
(  585)       ! --------------






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 25

(  586)       IF( kt == nit000 ) THEN
(  587)       
(  588)       zdt = rdt
(  589)       IF( nacc == 1 ) zdt = rdtmin
(  590) 
(  591)          ! Reference latitude
(  592)          ! ------------------
(  593)          !                                           ! =======================
(  594) !!DB: Orca-related
(  595) !         IF( cp_cfg == "orca" ) THEN                 !   ORCA configurations
(  596) !            !                                        ! =======================!
(  597) !
(  598) !            IF( jp_cfg == 05  )   iline = 192   ! i-line that passes near the North Pole
(  599) !            IF( jp_cfg == 025 )   iline = 384   ! i-line that passes near the North Pole
(  600) !            IF( jp_cfg == 2   )   iline =  48   ! i-line that passes near the North Pole
(  601) !            IF( jp_cfg == 4   )   iline =  24   ! i-line that passes near the North Pole
(  602) !            zphi(:) = 0.e0
(  603) !            DO ji = mi0(iline), mi1(iline) 
(  604) !               zphi(:) = gphiv(ji,:)         ! if iline is in the local domain
(  605) !               ! correct highest latitude for ORCA05
(  606) !               IF( jp_cfg == 05  ) zphi(jpj) = zphi(jpjm1) + (zphi(jpjm1)-zphi(jpj-2))/2.
(  607) !               IF( jp_cfg == 05  ) zphi(jpj) = MIN( zphi(jpj), 90.)!!!
(  608) !
(  609) !            END DO
(  610) !            ! provide the correct zphi to all local domains
(  611) !            IF( lk_mpp )   CALL mpp_sum( zphi, jpj )        !
(  612) !
(  613) !            !                                        ! =======================
(  614) !         ELSE                                        !   OTHER configurations!
(  615)             !                                        ! =======================!
(  616)             zphi(:) = gphiv(1,:)             ! assume lat/lon coordinate, select the first i-line
(  617)             !
(  618) !         ENDIF
(  619) 
(  620)          ! OPEN netcdf file 
(  621)          ! ----------------
(  622)          ! Define frequency of output and means
(  623)          zsto = nf_ptr * zdt
(  624)          clop = "ave(x)"
(  625)          zout = nf_ptr * zdt
(  626)          zfoo(:) = 0.e0
(  627) 
(  628)          ! Compute julian date from starting date of the run
(  629) 
(  630)          CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
(  631) 
(  632)          CALL dia_nam( clhstnam, nf_ptr, 'diaptr' )
(  633)          IF(lwp)WRITE( numout,*)" Name of diaptr NETCDF file ",clhstnam
(  634) 
(  635)          ! Horizontal grid : zphi()
(  636)          CALL histbeg(clhstnam, 1, zfoo, jpj, zphi,   &
(  637)             1, 1, 1, jpj, 0, zjulian, zdt, nhoridz, numptr, domain_id=nidom )
(  638)          ! Vertical grids : gdept, gdepw
(  639)          CALL histvert( numptr, "deptht", "Vertical T levels",   &
(  640)             "m", jpk, gdept, ndepidzt )
(  641)          CALL histvert( numptr, "depthw", "Vertical W levels",   &
(  642)             "m", jpk, gdepw, ndepidzw )
(  643)          






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 26

(  644)          !  Zonal mean T and S
(  645)          
(  646)          CALL histdef( numptr, "zotemglo", "Zonal Mean Temperature","C" ,   &
(  647)             1, jpj, nhoridz, jpk, 1, jpk, ndepidzt, 32, clop, zsto, zout )
(  648)          CALL histdef( numptr, "zosalglo", "Zonal Mean Salinity","PSU"  ,   &
(  649)             1, jpj, nhoridz, jpk, 1, jpk, ndepidzt, 32, clop, zsto, zout )
(  650) 
(  651)          !  Meridional Stream-Function (eulerian and bolus)
(  652)          
(  653)          CALL histdef( numptr, "zomsfglo", "Meridional Stream-Function: Global","Sv" ,   &
(  654)             1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
(  655)          IF( ln_subbas ) THEN
(  656)             CALL histdef( numptr, "zomsfatl", "Meridional Stream-Function: Atlantic","Sv" ,   &
(  657)                1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
(  658)             CALL histdef( numptr, "zomsfipc", "Meridional Stream-Function: Indo-Pacific","Sv" ,&
(  659)                1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
(  660)          ENDIF
(  661) 
(  662)          !  Heat transport 
(  663) 
(  664)          CALL histdef( numptr, "sophtadv", "Advective Heat Transport"      ,   &
(  665)             "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  666)          CALL histdef( numptr, "sophtldf", "Diffusive Heat Transport"      ,   &
(  667)             "PW",1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  668)          CALL histdef( numptr, "sophtove", "Overturning Heat Transport"    ,   &
(  669)             "PW",1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  670)          IF( ln_subbas ) THEN
(  671)             CALL histdef( numptr, "sohtatl", "Heat Transport Atlantic"      ,  &
(  672)                "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  673)             CALL histdef( numptr, "sohtpac", "Heat Transport Pacific"      ,   &
(  674)                "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  675)             CALL histdef( numptr, "sohtind", "Heat Transport Indic"      ,     &
(  676)                "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  677)          ENDIF
(  678) 
(  679) 
(  680)          !  Salt transport 
(  681) 
(  682)          CALL histdef( numptr, "sopstadv", "Advective Salt Transport"      ,   &
(  683)             "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  684)          CALL histdef( numptr, "sopstldf", "Diffusive Salt Transport"      ,   &
(  685)             "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  686)          CALL histdef( numptr, "sopstove", "Overturning Salt Transport"    ,   &
(  687)             "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  688) 
(  689) # 698
(  698)          IF( ln_subbas ) THEN
(  699)             CALL histdef( numptr, "sostatl", "Salt Transport Atlantic"      ,    &
(  700)                "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  701)             CALL histdef( numptr, "sostpac", "Salt Transport Pacific"      ,     &
(  702)                "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  703)             CALL histdef( numptr, "sostind", "Salt Transport Indic"      ,       &
(  704)                "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
(  705)          ENDIF
(  706)          
(  707) 
(  708)          CALL histend( numptr )
(  709) 






PGF90 (Version     10.5)          09/23/2010  09:49:28      page 27

(  710)       ENDIF
(  711) 
(  712)       IF( MOD( kt, nf_ptr ) == 0 ) THEN
(  713) 
(  714)          IF(lwp) THEN
(  715)             WRITE(numout,*)
(  716)             WRITE(numout,*) 'dia_ptr : write Poleward Transports at time-step : ', kt
(  717)             WRITE(numout,*) '~~~~~~~~'
(  718)             WRITE(numout,*)
(  719)          ENDIF
(  720) 
(  721)          ! define time axis
(  722)          it= kt - nit000 + 1
(  723)          ndex(1) = 0
(  724)          CALL histwrite( numptr, "zotemglo", it, tn_jk    , jpj*jpk, ndex )
(  725)          CALL histwrite( numptr, "zosalglo", it, sn_jk    , jpj*jpk, ndex )
(  726)          ! overturning outputs:
(  727)          CALL histwrite( numptr, "zomsfglo", it, v_msf_glo , jpj*jpk, ndex )
(  728)          IF( ln_subbas ) THEN
(  729)             CALL histwrite( numptr, "zomsfatl", it, v_msf_atl , jpj*jpk, ndex )
(  730)             CALL histwrite( numptr, "zomsfipc", it, v_msf_ipc , jpj*jpk, ndex )
(  731)          ENDIF
(  732)          ! heat transport outputs:
(  733)          IF( ln_subbas ) THEN
(  734)             CALL histwrite( numptr, "sohtatl", it, ht_atl  , jpj, ndex )
(  735)             CALL histwrite( numptr, "sohtpac", it, ht_pac  , jpj, ndex )
(  736)             CALL histwrite( numptr, "sohtind", it, ht_ind  , jpj, ndex )
(  737)             CALL histwrite( numptr, "sostatl", it, st_atl  , jpj, ndex )
(  738)             CALL histwrite( numptr, "sostpac", it, st_pac  , jpj, ndex )
(  739)             CALL histwrite( numptr, "sostind", it, st_ind  , jpj, ndex )
(  740)          ENDIF
(  741) 
(  742)          CALL histwrite( numptr, "sophtadv", it, pht_adv  , jpj, ndex )
(  743)          CALL histwrite( numptr, "sophtldf", it, pht_ldf  , jpj, ndex )
(  744)          CALL histwrite( numptr, "sophtove", it, pht_ove  , jpj, ndex )
(  745)          CALL histwrite( numptr, "sopstadv", it, pst_adv  , jpj, ndex )
(  746)          CALL histwrite( numptr, "sopstldf", it, pst_ldf  , jpj, ndex )
(  747)          CALL histwrite( numptr, "sopstove", it, pst_ove  , jpj, ndex )
(  748) # 753
(  753)  
(  754)       ENDIF
(  755) 
(  756)    END SUBROUTINE dia_ptr_wri
(  757) 
(  758)    !!======================================================================
(  759) END MODULE diaptr
