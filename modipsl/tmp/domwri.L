


PGF90 (Version     10.5)          09/23/2010  09:49:33      page 1

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

Filename: domwri.F90

(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) MODULE domwri
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE domwri  ***
(    5)    !! Ocean initialization : write the ocean domain mesh ask file(s)
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   dom_wri        : create and write mesh and mask file(s)
(   10)    !!                    nmsh = 1  :   mesh_mask file
(   11)    !!                         = 2  :   mesh and mask file
(   12)    !!                         = 3  :   mesh_hgr, mesh_zgr and mask
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE in_out_manager
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC dom_wri        ! routine called by inidom.F90
(   23)    !!----------------------------------------------------------------------
(   24)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   25)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domwri.F90,v 1.9 2006/03/10 10:55:38 opalod Exp $ 
(   26)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   27)    !!----------------------------------------------------------------------
(   28) 
(   29) CONTAINS
(   30) 
(   31) 
(   32)    !!----------------------------------------------------------------------
(   33)    !!   Default option :                                        NetCDF file
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36)    SUBROUTINE dom_wri
(   37)       !!----------------------------------------------------------------------
(   38)       !!                  ***  ROUTINE dom_wri  ***
(   39)       !!                   
(   40)       !! ** Purpose :   Create the NetCDF file(s) which contain(s) all the
(   41)       !!      ocean domain informations (mesh and mask arrays). This (these)
(   42)       !!      file(s) is (are) used for visualisation (SAXO software) and
(   43)       !!      diagnostic computation.
(   44)       !!
(   45)       !! ** Method  :   Write in a file all the arrays generated in routines
(   46)       !!      domhgr, domzgr, and dommsk. Note: the file contain depends on






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 2

(   47)       !!      the vertical coord. used (z-coord, partial steps, s-coord)
(   48)       !!                    nmsh = 1  :   'mesh_mask.nc' file
(   49)       !!                         = 2  :   'mesh.nc' and mask.nc' files
(   50)       !!                         = 3  :   'mesh_hgr.nc', 'mesh_zgr.nc' and
(   51)       !!                                  'mask.nc' files
(   52)       !!      For huge size domain, use option 2 or 3 depending on your 
(   53)       !!      vertical coordinate.
(   54)       !!
(   55)       !! ** output file : 
(   56)       !!      meshmask.nc  : domain size, horizontal grid-point position,
(   57)       !!                     masks, depth and vertical scale factors
(   58)       !!
(   59)       !! History :
(   60)       !!        !  97-02  (G. Madec)  Original code
(   61)       !!        !  99-11  (M. Imbard)  NetCDF FORMAT with IOIPSL
(   62)       !!   9.0  !  02-08  (G. Madec)  F90 and several file
(   63)       !!----------------------------------------------------------------------
(   64)       !! * Modules used
(   65)       USE ioipsl
(   66) 
(   67)       !! * Local declarations
(   68)       INTEGER  ::                & !!! * temprary units for :
(   69)          inum0 ,                 &  ! 'mesh_mask.nc' file
(   70)          inum1 ,                 &  ! 'mesh.nc'      file
(   71)          inum2 ,                 &  ! 'mask.nc'      file
(   72)          inum3 ,                 &  ! 'mesh_hgr.nc'  file
(   73)          inum4                      ! 'mesh_zgr.nc'  file
(   74)       INTEGER  ::   itime           !  output from restini ???
(   75)       REAL(wp) ::   zdate0
(   76)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   77)          zprt                       ! temporary array for bathymetry 
(   78) 
(   79)       CHARACTER (len=21) ::      &
(   80)          clnam0  ,   &  ! filename (mesh and mask informations)
(   81)          clnam1  ,   &  ! filename (mesh informations)
(   82)          clnam2  ,   &  ! filename (mask informations)
(   83)          clnam3  ,   &  ! filename (horizontal mesh informations)
(   84)          clnam4         ! filename (vertical   mesh informations)
(   85)       !!----------------------------------------------------------------------
(   86) 
(   87)        IF(lwp) WRITE(numout,*)
(   88)        IF(lwp) WRITE(numout,*) 'dom_wri : create NetCDF mesh and mask information file(s)'
(   89)        IF(lwp) WRITE(numout,*) '~~~~~~~'
(   90) 
(   91)          clnam0 = 'mesh_mask'  ! filename (mesh and mask informations)
(   92)          clnam1 = 'mesh'       ! filename (mesh informations)
(   93)          clnam2 = 'mask'       ! filename (mask informations)
(   94)          clnam3 = 'mesh_hgr'   ! filename (horizontal mesh informations)
(   95)          clnam4 = 'mesh_zgr'   ! filename (vertical   mesh informations)
(   96) 
(   97) #if defined key_agrif
(   98)       if ( .NOT. Agrif_Root() ) then
(   99)         clnam0 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam0)
(  100)         clnam1 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam1)
(  101)         clnam2 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam2)
(  102)         clnam3 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam3)
(  103)         clnam4 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam4)
(  104)       endif






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 3

(  105) #endif
(  106) 
(  107)       CALL ymds2ju( 0, 1, 1, 0.e0, zdate0 )    ! calendar initialization
(  108) 
(  109) !       note that mbathy has been modified in dommsk or in solver.
(  110) !       it is the number of non-zero "w" levels in the water, and the minimum 
(  111) !       value (on land) is 2. We define zprt as the number of "T" points in the ocean 
(  112) !       at any location, and zero on land. 
(  113) !
(  114)       zprt = tmask(:,:,1)*(mbathy-1)
(  115) 
(  116)       SELECT CASE (nmsh)
(  117)          !                                     ! ============================
(  118)          CASE ( 1 )                            !  create 'mesh_mask.nc' file
(  119)             !                                  ! ============================
(  120) 
(  121)             IF(lwp) WRITE(numout,*) '          one file in "mesh_mask.nc" '
(  122)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_mask.nc' file
(  123)             &             jpk   , gdept , trim(clnam0)        ,  &   ! in unit inum0
(  124)             &             itime , zdate0, rdt   , inum0 , domain_id=nidom )
(  125)             inum2 = inum0                                            ! put all the informations
(  126)             inum3 = inum0                                            ! in unit inum0
(  127)             inum4 = inum0
(  128) 
(  129)             !                                  ! ============================
(  130)          CASE ( 2 )                            !  create 'mesh.nc' and 
(  131)             !                                  !         'mask.nc' files
(  132)             !                                  ! ============================
(  133) 
(  134)             IF(lwp) WRITE(numout,*) '          two files in "mesh.nc" and "mask.nc" '
(  135)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh.nc' file 
(  136)             &             jpk   , gdept , trim(clnam1)        ,  &   ! in unit inum1 
(  137)             &             itime , zdate0, rdt   , inum1, domain_id=nidom )
(  138)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mask.nc' file 
(  139)             &             jpk   , gdept , trim(clnam2)        ,  &   ! in unit inum2 
(  140)             &             itime , zdate0, rdt   , inum2, domain_id=nidom )
(  141)             inum3 = inum1                                            ! put mesh informations 
(  142)             inum4 = inum1                                            ! in unit inum1 
(  143) 
(  144)             !                                  ! ============================
(  145)          CASE ( 3 )                            !  create 'mesh_hgr.nc'
(  146)             !                                  !         'mesh_zgr.nc' and
(  147)             !                                  !         'mask.nc'     files
(  148)             !                                  ! ============================
(  149) 
(  150)             IF(lwp) WRITE(numout,*) '          three files in "mesh_hgr.nc" , mesh_zgr.nc" and "mask.nc" '
(  151)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_hgr.nc' file
(  152)             &             jpk   , gdept , trim(clnam3)        ,  &   ! in unit inum3
(  153)             &             itime , zdate0, rdt   , inum3, domain_id=nidom )
(  154)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_zgr.nc' file
(  155)             &             jpk   , gdept , trim(clnam4)        ,  &   ! in unit inum4
(  156)             &             itime , zdate0, rdt   , inum4, domain_id=nidom )
(  157)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mask.nc' file
(  158)             &             jpk   , gdept , trim(clnam2)        ,  &   ! in unit inum2
(  159)             &             itime , zdate0, rdt   , inum2, domain_id=nidom )
(  160) 
(  161)          END SELECT
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 4

(  163)          !                                                         ! masks (inum2) 
(  164)          CALL restput( inum2, 'tmask', jpi, jpj, jpk, 0, tmask ) 
(  165)          CALL restput( inum2, 'umask', jpi, jpj, jpk, 0, umask )
(  166)          CALL restput( inum2, 'vmask', jpi, jpj, jpk, 0, vmask )
(  167)          CALL restput( inum2, 'fmask', jpi, jpj, jpk, 0, fmask )
(  168) 
(  169)          !                                                         ! horizontal mesh (inum3)
(  170)          CALL restput( inum3, 'glamt', jpi, jpj, 1, 0, glamt )     !    ! latitude
(  171)          CALL restput( inum3, 'glamu', jpi, jpj, 1, 0, glamu )
(  172)          CALL restput( inum3, 'glamv', jpi, jpj, 1, 0, glamv )
(  173)          CALL restput( inum3, 'glamf', jpi, jpj, 1, 0, glamf )
(  174) 
(  175)          CALL restput( inum3, 'gphit', jpi, jpj, 1, 0, gphit )     !    ! longitude
(  176)          CALL restput( inum3, 'gphiu', jpi, jpj, 1, 0, gphiu )
(  177)          CALL restput( inum3, 'gphiv', jpi, jpj, 1, 0, gphiv )
(  178)          CALL restput( inum3, 'gphif', jpi, jpj, 1, 0, gphif )
(  179) 
(  180)          CALL restput( inum3, 'e1t', jpi, jpj, 1, 0, e1t )         !    ! e1 scale factors
(  181)          CALL restput( inum3, 'e1u', jpi, jpj, 1, 0, e1u )
(  182)          CALL restput( inum3, 'e1v', jpi, jpj, 1, 0, e1v )
(  183)          CALL restput( inum3, 'e1f', jpi, jpj, 1, 0, e1f )
(  184) 
(  185)          CALL restput( inum3, 'e2t', jpi, jpj, 1, 0, e2t )         !    ! e2 scale factors
(  186)          CALL restput( inum3, 'e2u', jpi, jpj, 1, 0, e2u )
(  187)          CALL restput( inum3, 'e2v', jpi, jpj, 1, 0, e2v )
(  188)          CALL restput( inum3, 'e2f', jpi, jpj, 1, 0, e2f )
(  189) 
(  190)          CALL restput( inum3, 'ff', jpi, jpj, 1, 0, ff )           !    ! coriolis factor
(  191) 
(  192)          CALL restput( inum4, 'mbathy', jpi, jpj, 1, 0, zprt )
(  193) 
(  194) # if defined key_s_coord
(  195)          !                                                         ! s-coordinate
(  196)          CALL restput( inum4, 'hbatt', jpi, jpj, 1, 0, hbatt )      !    ! depth
(  197)          CALL restput( inum4, 'hbatu', jpi, jpj, 1, 0, hbatu ) 
(  198)          CALL restput( inum4, 'hbatv', jpi, jpj, 1, 0, hbatv )
(  199)          CALL restput( inum4, 'hbatf', jpi, jpj, 1, 0, hbatf )
(  200) 
(  201)          CALL restput( inum4, 'gsigt', 1, 1, jpk, 0, gsigt )        !    ! scaling coef.
(  202)          CALL restput( inum4, 'gsigw', 1, 1, jpk, 0, gsigw )  
(  203)          CALL restput( inum4, 'gsi3w', 1, 1, jpk, 0, gsi3w )
(  204)          CALL restput( inum4, 'esigt', 1, 1, jpk, 0, esigt )
(  205)          CALL restput( inum4, 'esigw', 1, 1, jpk, 0, esigw )
(  206) 
(  207) # elif defined key_partial_steps
(  208)          !                                                          ! z-coordinate with partial steps
(  209)          CALL restput( inum4, 'hdept' , jpi, jpj, 1, 0, hdept  )    !    ! depth
(  210)          CALL restput( inum4, 'hdepw' , jpi, jpj, 1, 0, hdepw  ) 
(  211) 
(  212)          CALL restput( inum4, 'e3t_ps', jpi, jpj, jpk, 0, e3t_ps )  !    ! scale factors
(  213)          CALL restput( inum4, 'e3u_ps', jpi, jpj, jpk, 0, e3u_ps )
(  214)          CALL restput( inum4, 'e3v_ps', jpi, jpj, jpk, 0, e3v_ps )
(  215)          CALL restput( inum4, 'e3w_ps', jpi, jpj, jpk, 0, e3w_ps )
(  216) 
(  217)          CALL restput( inum4, 'gdept' , 1, 1, jpk, 0, gdept )       !    ! reference z-coord.
(  218)          CALL restput( inum4, 'gdepw' , 1, 1, jpk, 0, gdepw )
(  219)          CALL restput( inum4, 'e3t'   , 1, 1, jpk, 0, e3t   )
(  220)          CALL restput( inum4, 'e3w'   , 1, 1, jpk, 0, e3w   )






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 5

(  221) 
(  222) # else
(  223)          !                                                          ! z-coordinate 
(  224)          CALL restput( inum4, 'gdept', 1, 1, jpk, 0, gdept )        !    ! depth
(  225)          CALL restput( inum4, 'gdepw', 1, 1, jpk, 0, gdepw )
(  226)          CALL restput( inum4, 'e3t'  , 1, 1, jpk, 0, e3t   )        !    ! scale factors
(  227)          CALL restput( inum4, 'e3w'  , 1, 1, jpk, 0, e3w   )
(  228) # endif
(  229) 
(  230)          !                                     ! ============================
(  231)          !                                     !        close the files 
(  232)          !                                     ! ============================
(  233)          SELECT CASE ( nmsh )
(  234)             CASE ( 1 )                
(  235)                CALL restclo( inum0 )
(  236)             CASE ( 2 )
(  237)                CALL restclo( inum1 )
(  238)                CALL restclo( inum2 )
(  239)             CASE ( 3 )
(  240)                CALL restclo( inum2 )
(  241)                CALL restclo( inum3 )
(  242)                CALL restclo( inum4 )
(  243)          END SELECT
(  244) 
(  245)    END SUBROUTINE dom_wri
(  246) 
(  247) 
(  248)    !!======================================================================
(  249) END MODULE domwri



































PGF90 (Version     10.5)          09/23/2010  09:49:33      page 6

(    1) # 1 "domwri.F90"
(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) MODULE domwri
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE domwri  ***
(    5)    !! Ocean initialization : write the ocean domain mesh ask file(s)
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   dom_wri        : create and write mesh and mask file(s)
(   10)    !!                    nmsh = 1  :   mesh_mask file
(   11)    !!                         = 2  :   mesh and mask file
(   12)    !!                         = 3  :   mesh_hgr, mesh_zgr and mask
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE in_out_manager
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC dom_wri        ! routine called by inidom.F90
(   23)    !!----------------------------------------------------------------------
(   24)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   25)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domwri.F90,v 1.9 2006/03/10 10:55:38 opalod Exp $ 
(   26)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   27)    !!----------------------------------------------------------------------
(   28) 
(   29) CONTAINS
(   30) 
(   31) 
(   32)    !!----------------------------------------------------------------------
(   33)    !!   Default option :                                        NetCDF file
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36)    SUBROUTINE dom_wri
(   37)       !!----------------------------------------------------------------------
(   38)       !!                  ***  ROUTINE dom_wri  ***
(   39)       !!                   
(   40)       !! ** Purpose :   Create the NetCDF file(s) which contain(s) all the
(   41)       !!      ocean domain informations (mesh and mask arrays). This (these)
(   42)       !!      file(s) is (are) used for visualisation (SAXO software) and
(   43)       !!      diagnostic computation.
(   44)       !!
(   45)       !! ** Method  :   Write in a file all the arrays generated in routines
(   46)       !!      domhgr, domzgr, and dommsk. Note: the file contain depends on
(   47)       !!      the vertical coord. used (z-coord, partial steps, s-coord)
(   48)       !!                    nmsh = 1  :   'mesh_mask.nc' file
(   49)       !!                         = 2  :   'mesh.nc' and mask.nc' files
(   50)       !!                         = 3  :   'mesh_hgr.nc', 'mesh_zgr.nc' and
(   51)       !!                                  'mask.nc' files
(   52)       !!      For huge size domain, use option 2 or 3 depending on your 
(   53)       !!      vertical coordinate.
(   54)       !!
(   55)       !! ** output file : 
(   56)       !!      meshmask.nc  : domain size, horizontal grid-point position,
(   57)       !!                     masks, depth and vertical scale factors






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 7

(   58)       !!
(   59)       !! History :
(   60)       !!        !  97-02  (G. Madec)  Original code
(   61)       !!        !  99-11  (M. Imbard)  NetCDF FORMAT with IOIPSL
(   62)       !!   9.0  !  02-08  (G. Madec)  F90 and several file
(   63)       !!----------------------------------------------------------------------
(   64)       !! * Modules used
(   65)       USE ioipsl
(   66) 
(   67)       !! * Local declarations
(   68)       INTEGER  ::                & !!! * temprary units for :
(   69)          inum0 ,                 &  ! 'mesh_mask.nc' file
(   70)          inum1 ,                 &  ! 'mesh.nc'      file
(   71)          inum2 ,                 &  ! 'mask.nc'      file
(   72)          inum3 ,                 &  ! 'mesh_hgr.nc'  file
(   73)          inum4                      ! 'mesh_zgr.nc'  file
(   74)       INTEGER  ::   itime           !  output from restini ???
(   75)       REAL(wp) ::   zdate0
(   76)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   77)          zprt                       ! temporary array for bathymetry 
(   78) 
(   79)       CHARACTER (len=21) ::      &
(   80)          clnam0  ,   &  ! filename (mesh and mask informations)
(   81)          clnam1  ,   &  ! filename (mesh informations)
(   82)          clnam2  ,   &  ! filename (mask informations)
(   83)          clnam3  ,   &  ! filename (horizontal mesh informations)
(   84)          clnam4         ! filename (vertical   mesh informations)
(   85)       !!----------------------------------------------------------------------
(   86) 
(   87)        IF(lwp) WRITE(numout,*)
(   88)        IF(lwp) WRITE(numout,*) 'dom_wri : create NetCDF mesh and mask information file(s)'
(   89)        IF(lwp) WRITE(numout,*) '~~~~~~~'
(   90) 
(   91)          clnam0 = 'mesh_mask'  ! filename (mesh and mask informations)
(   92)          clnam1 = 'mesh'       ! filename (mesh informations)
(   93)          clnam2 = 'mask'       ! filename (mask informations)
(   94)          clnam3 = 'mesh_hgr'   ! filename (horizontal mesh informations)
(   95)          clnam4 = 'mesh_zgr'   ! filename (vertical   mesh informations)
(   96) 
(   97) 
(   98) # 107
(  107)       CALL ymds2ju( 0, 1, 1, 0.e0, zdate0 )    ! calendar initialization
(  108) 
(  109) !       note that mbathy has been modified in dommsk or in solver.
(  110) !       it is the number of non-zero "w" levels in the water, and the minimum 
(  111) !       value (on land) is 2. We define zprt as the number of "T" points in the ocean 
(  112) !       at any location, and zero on land. 
(  113) !
(  114)       zprt = tmask(:,:,1)*(mbathy-1)
(  115) 
(  116)       SELECT CASE (nmsh)
(  117)          !                                     ! ============================
(  118)          CASE ( 1 )                            !  create 'mesh_mask.nc' file
(  119)             !                                  ! ============================
(  120) 
(  121)             IF(lwp) WRITE(numout,*) '          one file in "mesh_mask.nc" '
(  122)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_mask.nc' file
(  123)             &             jpk   , gdept , trim(clnam0)        ,  &   ! in unit inum0






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 8

(  124)             &             itime , zdate0, rdt   , inum0 , domain_id=nidom )
(  125)             inum2 = inum0                                            ! put all the informations
(  126)             inum3 = inum0                                            ! in unit inum0
(  127)             inum4 = inum0
(  128) 
(  129)             !                                  ! ============================
(  130)          CASE ( 2 )                            !  create 'mesh.nc' and 
(  131)             !                                  !         'mask.nc' files
(  132)             !                                  ! ============================
(  133) 
(  134)             IF(lwp) WRITE(numout,*) '          two files in "mesh.nc" and "mask.nc" '
(  135)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh.nc' file 
(  136)             &             jpk   , gdept , trim(clnam1)        ,  &   ! in unit inum1 
(  137)             &             itime , zdate0, rdt   , inum1, domain_id=nidom )
(  138)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mask.nc' file 
(  139)             &             jpk   , gdept , trim(clnam2)        ,  &   ! in unit inum2 
(  140)             &             itime , zdate0, rdt   , inum2, domain_id=nidom )
(  141)             inum3 = inum1                                            ! put mesh informations 
(  142)             inum4 = inum1                                            ! in unit inum1 
(  143) 
(  144)             !                                  ! ============================
(  145)          CASE ( 3 )                            !  create 'mesh_hgr.nc'
(  146)             !                                  !         'mesh_zgr.nc' and
(  147)             !                                  !         'mask.nc'     files
(  148)             !                                  ! ============================
(  149) 
(  150)             IF(lwp) WRITE(numout,*) '          three files in "mesh_hgr.nc" , mesh_zgr.nc" and "mask.nc" '
(  151)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_hgr.nc' file
(  152)             &             jpk   , gdept , trim(clnam3)        ,  &   ! in unit inum3
(  153)             &             itime , zdate0, rdt   , inum3, domain_id=nidom )
(  154)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_zgr.nc' file
(  155)             &             jpk   , gdept , trim(clnam4)        ,  &   ! in unit inum4
(  156)             &             itime , zdate0, rdt   , inum4, domain_id=nidom )
(  157)             CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mask.nc' file
(  158)             &             jpk   , gdept , trim(clnam2)        ,  &   ! in unit inum2
(  159)             &             itime , zdate0, rdt   , inum2, domain_id=nidom )
(  160) 
(  161)          END SELECT
(  162) 
(  163)          !                                                         ! masks (inum2) 
(  164)          CALL restput( inum2, 'tmask', jpi, jpj, jpk, 0, tmask ) 
(  165)          CALL restput( inum2, 'umask', jpi, jpj, jpk, 0, umask )
(  166)          CALL restput( inum2, 'vmask', jpi, jpj, jpk, 0, vmask )
(  167)          CALL restput( inum2, 'fmask', jpi, jpj, jpk, 0, fmask )
(  168) 
(  169)          !                                                         ! horizontal mesh (inum3)
(  170)          CALL restput( inum3, 'glamt', jpi, jpj, 1, 0, glamt )     !    ! latitude
(  171)          CALL restput( inum3, 'glamu', jpi, jpj, 1, 0, glamu )
(  172)          CALL restput( inum3, 'glamv', jpi, jpj, 1, 0, glamv )
(  173)          CALL restput( inum3, 'glamf', jpi, jpj, 1, 0, glamf )
(  174) 
(  175)          CALL restput( inum3, 'gphit', jpi, jpj, 1, 0, gphit )     !    ! longitude
(  176)          CALL restput( inum3, 'gphiu', jpi, jpj, 1, 0, gphiu )
(  177)          CALL restput( inum3, 'gphiv', jpi, jpj, 1, 0, gphiv )
(  178)          CALL restput( inum3, 'gphif', jpi, jpj, 1, 0, gphif )
(  179) 
(  180)          CALL restput( inum3, 'e1t', jpi, jpj, 1, 0, e1t )         !    ! e1 scale factors
(  181)          CALL restput( inum3, 'e1u', jpi, jpj, 1, 0, e1u )






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 9

(  182)          CALL restput( inum3, 'e1v', jpi, jpj, 1, 0, e1v )
(  183)          CALL restput( inum3, 'e1f', jpi, jpj, 1, 0, e1f )
(  184) 
(  185)          CALL restput( inum3, 'e2t', jpi, jpj, 1, 0, e2t )         !    ! e2 scale factors
(  186)          CALL restput( inum3, 'e2u', jpi, jpj, 1, 0, e2u )
(  187)          CALL restput( inum3, 'e2v', jpi, jpj, 1, 0, e2v )
(  188)          CALL restput( inum3, 'e2f', jpi, jpj, 1, 0, e2f )
(  189) 
(  190)          CALL restput( inum3, 'ff', jpi, jpj, 1, 0, ff )           !    ! coriolis factor
(  191) 
(  192)          CALL restput( inum4, 'mbathy', jpi, jpj, 1, 0, zprt )
(  193) 
(  194) # 223
(  223)          !                                                          ! z-coordinate 
(  224)          CALL restput( inum4, 'gdept', 1, 1, jpk, 0, gdept )        !    ! depth
(  225)          CALL restput( inum4, 'gdepw', 1, 1, jpk, 0, gdepw )
(  226)          CALL restput( inum4, 'e3t'  , 1, 1, jpk, 0, e3t   )        !    ! scale factors
(  227)          CALL restput( inum4, 'e3w'  , 1, 1, jpk, 0, e3w   )
(  228) 
(  229) # 230
(  230)          !                                     ! ============================
(  231)          !                                     !        close the files 
(  232)          !                                     ! ============================
(  233)          SELECT CASE ( nmsh )
(  234)             CASE ( 1 )                
(  235)                CALL restclo( inum0 )
(  236)             CASE ( 2 )
(  237)                CALL restclo( inum1 )
(  238)                CALL restclo( inum2 )
(  239)             CASE ( 3 )
(  240)                CALL restclo( inum2 )
(  241)                CALL restclo( inum3 )
(  242)                CALL restclo( inum4 )
(  243)          END SELECT
(  244) 
(  245)    END SUBROUTINE dom_wri
(  246) 
(  247) 
(  248)    !!======================================================================
(  249) END MODULE domwri
