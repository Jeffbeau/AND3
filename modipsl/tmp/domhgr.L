


PGF90 (Version     10.5)          09/23/2010  09:49:32      page 1

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

Filename: domhgr.F90

(    1) !!DB 2008.09.30 -- Eliminated IOIPSL in hgr_read
(    2) !!   2009.08.31 -- Eliminated GYRE (jphgr_msh=5)
(    3) MODULE domhgr
(    4)    !!==============================================================================
(    5)    !!                       ***  MODULE domhgr   ***
(    6)    !! Ocean initialization : domain initialization
(    7)    !!==============================================================================
(    8) 
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dom_hgr        : initialize the horizontal mesh 
(   11)    !!   hgr_read       : read "coordinate" NetCDF file 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE phycst          ! physical constants
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE lib_mpp
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Module variables
(   23)    REAL(wp) ::   glam0, gphi0           ! variables corresponding to parameters
(   24)       !                                 ! ppglam0 ppgphi0 set in par_oce
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC dom_hgr        ! called by domain.F90
(   28)    !!----------------------------------------------------------------------
(   29)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domhgr.F90,v 1.15 2006/04/10 15:46:07 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36)    SUBROUTINE dom_hgr
(   37)       !!----------------------------------------------------------------------
(   38)       !!                  ***  ROUTINE dom_hgr  ***
(   39)       !!
(   40)       !! ** Purpose :   Compute the geographical position (in degre) of the 
(   41)       !!      model grid-points,  the horizontal scale factors (in meters) and 
(   42)       !!      the Coriolis factor (in s-1).
(   43)       !!
(   44)       !! ** Method  :   The geographical position of the model grid-points is
(   45)       !!      defined from analytical functions, fslam and fsphi, the deriva-
(   46)       !!      tives of which gives the horizontal scale factors e1,e2.






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 2

(   47)       !!      Defining two function fslam and fsphi and their derivatives in 
(   48)       !!      the two horizontal directions (fse1 and fse2), the model grid-
(   49)       !!      point position and scale factors are given by:
(   50)       !!         t-point:
(   51)       !!      glamt(i,j) = fslam(i    ,j    )   e1t(i,j) = fse1(i    ,j    )
(   52)       !!      gphit(i,j) = fsphi(i    ,j    )   e2t(i,j) = fse2(i    ,j    )
(   53)       !!         u-point:
(   54)       !!      glamu(i,j) = fslam(i+1/2,j    )   e1u(i,j) = fse1(i+1/2,j    )
(   55)       !!      gphiu(i,j) = fsphi(i+1/2,j    )   e2u(i,j) = fse2(i+1/2,j    )
(   56)       !!         v-point:
(   57)       !!      glamv(i,j) = fslam(i    ,j+1/2)   e1v(i,j) = fse1(i    ,j+1/2)
(   58)       !!      gphiv(i,j) = fsphi(i    ,j+1/2)   e2v(i,j) = fse2(i    ,j+1/2)
(   59)       !!            f-point:
(   60)       !!      glamf(i,j) = fslam(i+1/2,j+1/2)   e1f(i,j) = fse1(i+1/2,j+1/2)
(   61)       !!      gphif(i,j) = fsphi(i+1/2,j+1/2)   e2f(i,j) = fse2(i+1/2,j+1/2)
(   62)       !!      Where fse1 and fse2 are defined by:
(   63)       !!         fse1(i,j) = ra * rad * SQRT( (cos(phi) di(fslam))**2
(   64)       !!                                     +          di(fsphi) **2 )(i,j)
(   65)       !!         fse2(i,j) = ra * rad * SQRT( (cos(phi) dj(fslam))**2
(   66)       !!                                     +          dj(fsphi) **2 )(i,j)
(   67)       !!
(   68)       !!        The coriolis factor is given at z-point by:
(   69)       !!                     ff = 2.*omega*sin(gphif)      (in s-1)
(   70)       !!
(   71)       !!        This routine is given as an example, it must be modified
(   72)       !!      following the user s desiderata. nevertheless, the output as
(   73)       !!      well as the way to compute the model grid-point position and
(   74)       !!      horizontal scale factors must be respected in order to insure
(   75)       !!      second order accuracy schemes.
(   76)       !!
(   77)       !! N.B. If the domain is periodic, verify that scale factors are also
(   78)       !!      periodic, and the coriolis term again.
(   79)       !!
(   80)       !! ** Action  : - define  glamt, glamu, glamv, glamf: longitude of t-, 
(   81)       !!                u-, v- and f-points (in degre)
(   82)       !!              - define  gphit, gphiu, gphiv, gphit: latitude  of t-,
(   83)       !!               u-, v-  and f-points (in degre)
(   84)       !!        define e1t, e2t, e1u, e2u, e1v, e2v, e1f, e2f: horizontal
(   85)       !!      scale factors (in meters) at t-, u-, v-, and f-points.
(   86)       !!        define ff: coriolis factor at f-point
(   87)       !!
(   88)       !! References :
(   89)       !!      Marti, Madec and Delecluse, 1992, j. geophys. res., in press.
(   90)       !!
(   91)       !! History :
(   92)       !!        !  88-03  (G. Madec)
(   93)       !!        !  91-11  (G. Madec)
(   94)       !!        !  92-06  (M. Imbard)
(   95)       !!        !  96-01  (G. Madec)  terrain following coordinates
(   96)       !!        !  97-02  (G. Madec)  print mesh informations
(   97)       !!        !  01-09  (M. Levy)  eel config: grid in km, beta-plane
(   98)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module, namelist
(   99)       !!   9.0  !  04-01  (A.M. Treguier, J.M. Molines) Case 4 (Mercator mesh)
(  100)       !!                  use of parameters in par_CONFIG-Rxx.h90, not in namelist
(  101)       !!        !  04-05  (A. Koch-Larrouy) Add Gyre configuration 
(  102)       !!----------------------------------------------------------------------
(  103)       !! * local declarations
(  104)       INTEGER  ::   ji, jj              ! dummy loop indices






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 3

(  105)       INTEGER  ::   ii0, ii1, ij0, ij1  ! temporary integers
(  106)       INTEGER  ::   ijeq                ! index of equator T point (used in case 4)
(  107)       REAL(wp) ::   &
(  108)          zti, zui, zvi, zfi,         &  ! temporary scalars
(  109)          ztj, zuj, zvj, zfj,         &  !
(  110)          zphi0, zbeta, znorme,       &  !
(  111)          zarg, zf0, zminff, zmaxff
(  112)       REAL(wp) ::   &
(  113)          zlam1, zcos_alpha, zim1 , zjm1 , ze1, ze1deg,   &
(  114)          zphi1, zsin_alpha, zim05, zjm05
(  115)       !!----------------------------------------------------------------------
(  116) 
(  117)       IF(lwp) THEN
(  118)          WRITE(numout,*)
(  119)          WRITE(numout,*) 'dom_hgr : define the horizontal mesh from ithe following par_oce parameters '
(  120)          WRITE(numout,*) '~~~~~~~      type of horizontal mesh           jphgr_msh = ', jphgr_msh
(  121)          WRITE(numout,*) '             position of the first row and     ppglam0  = ', ppglam0
(  122)          WRITE(numout,*) '             column grid-point (degrees)       ppgphi0  = ', ppgphi0
(  123)          WRITE(numout,*) '             zonal      grid-spacing (degrees) ppe1_deg = ', ppe1_deg
(  124)          WRITE(numout,*) '             meridional grid-spacing (degrees) ppe2_deg = ', ppe2_deg
(  125)          WRITE(numout,*) '             zonal      grid-spacing (meters)  ppe1_m   = ', ppe1_m  
(  126)          WRITE(numout,*) '             meridional grid-spacing (meters)  ppe2_m   = ', ppe2_m  
(  127)       ENDIF
(  128) 
(  129) 
(  130)       SELECT CASE( jphgr_msh )   ! type of horizontal mesh
(  131) 
(  132)       CASE ( 0 )                     !  curvilinear coordinate on the sphere read in coordinate.nc file
(  133) 
(  134)          IF(lwp) WRITE(numout,*)
(  135)          IF(lwp) WRITE(numout,*) '          curvilinear coordinate on the sphere read in "coordinate" file'
(  136) 
(  137)          CALL hgr_read           ! Defaultl option  :   NetCDF file
(  138) 
(  139)          !                                                ! =====================
(  140) !!DB: delete ORCA
(  141) !         IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN    ! ORCA R2 configuration
(  142) 
(  143)          !                                                ! ======================
(  144) !         IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN   ! ORCA R05 configuration
(  145)          ! N.B. :  General case, lat and long function of both i and j indices:
(  146)          !     e1t(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphit(ji,jj) ) * fsdila( zti, ztj ) )**2   &
(  147)          !                                  + (                           fsdiph( zti, ztj ) )**2  )
(  148)          !     e1u(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiu(ji,jj) ) * fsdila( zui, zuj ) )**2   &
(  149)          !                                  + (                           fsdiph( zui, zuj ) )**2  )
(  150)          !     e1v(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiv(ji,jj) ) * fsdila( zvi, zvj ) )**2   &
(  151)          !                                  + (                           fsdiph( zvi, zvj ) )**2  )
(  152)          !     e1f(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphif(ji,jj) ) * fsdila( zfi, zfj ) )**2   &
(  153)          !                                  + (                           fsdiph( zfi, zfj ) )**2  )
(  154)          !
(  155)          !     e2t(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphit(ji,jj) ) * fsdjla( zti, ztj ) )**2   &
(  156)          !                                  + (                           fsdjph( zti, ztj ) )**2  )
(  157)          !     e2u(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiu(ji,jj) ) * fsdjla( zui, zuj ) )**2   &
(  158)          !                                  + (                           fsdjph( zui, zuj ) )**2  )
(  159)          !     e2v(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiv(ji,jj) ) * fsdjla( zvi, zvj ) )**2   &
(  160)          !                                  + (                           fsdjph( zvi, zvj ) )**2  )
(  161)          !     e2f(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphif(ji,jj) ) * fsdjla( zfi, zfj ) )**2   &
(  162)          !                                  + (                           fsdjph( zfi, zfj ) )**2  )






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 4

(  163) 
(  164) 
(  165)       CASE ( 1 )                     ! geographical mesh on the sphere with regular grid-spacing
(  166) 
(  167)          IF(lwp) WRITE(numout,*)
(  168)          IF(lwp) WRITE(numout,*) '          geographical mesh on the sphere with regular grid-spacing'
(  169)          IF(lwp) WRITE(numout,*) '          given by ppe1_deg and ppe2_deg' 
(  170) 
(  171)          DO jj = 1, jpj
(  172)             DO ji = 1, jpi
(  173)                zti = FLOAT( ji - 1 + nimpp - 1 )         ;   ztj = FLOAT( jj - 1 + njmpp - 1 )
(  174)                zui = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zuj = FLOAT( jj - 1 + njmpp - 1 )
(  175)                zvi = FLOAT( ji - 1 + nimpp - 1 )         ;   zvj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5
(  176)                zfi = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zfj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5
(  177)          ! Longitude
(  178)                glamt(ji,jj) = ppglam0 + ppe1_deg * zti
(  179)                glamu(ji,jj) = ppglam0 + ppe1_deg * zui
(  180)                glamv(ji,jj) = ppglam0 + ppe1_deg * zvi
(  181)                glamf(ji,jj) = ppglam0 + ppe1_deg * zfi
(  182)          ! Latitude
(  183)                gphit(ji,jj) = ppgphi0 + ppe2_deg * ztj
(  184)                gphiu(ji,jj) = ppgphi0 + ppe2_deg * zuj
(  185)                gphiv(ji,jj) = ppgphi0 + ppe2_deg * zvj
(  186)                gphif(ji,jj) = ppgphi0 + ppe2_deg * zfj
(  187)          ! e1
(  188)                e1t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
(  189)                e1u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
(  190)                e1v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
(  191)                e1f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
(  192)          ! e2
(  193)                e2t(ji,jj) = ra * rad * ppe2_deg
(  194)                e2u(ji,jj) = ra * rad * ppe2_deg
(  195)                e2v(ji,jj) = ra * rad * ppe2_deg
(  196)                e2f(ji,jj) = ra * rad * ppe2_deg
(  197)             END DO
(  198)          END DO
(  199) 
(  200) 
(  201)       CASE ( 2:3 )                   ! f- or beta-plane with regular grid-spacing
(  202) 
(  203)          IF(lwp) WRITE(numout,*)
(  204)          IF(lwp) WRITE(numout,*) '          f- or beta-plane with regular grid-spacing'
(  205)          IF(lwp) WRITE(numout,*) '          given by ppe1_m and ppe2_m' 
(  206) 
(  207)          ! Position coordinates (in kilometers)
(  208)          !                          ==========
(  209)          glam0 = 0.e0
(  210)          gphi0 = - ppe2_m * 1.e-3
(  211)          
(  212) #if defined key_agrif && defined key_eel_r6
(  213)          IF (.Not.Agrif_Root()) THEN
(  214)            glam0  = Agrif_Parent(glam0) + (Agrif_ix())*Agrif_Parent(ppe1_m) * 1.e-3
(  215)            gphi0  = Agrif_Parent(gphi0) + (Agrif_iy())*Agrif_Parent(ppe2_m) * 1.e-3
(  216)            ppe1_m = Agrif_Parent(ppe1_m)/Agrif_Rhox()
(  217)            ppe2_m = Agrif_Parent(ppe2_m)/Agrif_Rhoy()          
(  218)          ENDIF
(  219) #endif         
(  220)          DO jj = 1, jpj






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 5

(  221)             DO ji = 1, jpi
(  222)                glamt(ji,jj) = glam0 + ppe1_m * 1.e-3 * ( FLOAT( ji - 1 + nimpp - 1 )       )
(  223)                glamu(ji,jj) = glam0 + ppe1_m * 1.e-3 * ( FLOAT( ji - 1 + nimpp - 1 ) + 0.5 )
(  224)                glamv(ji,jj) = glamt(ji,jj)
(  225)                glamf(ji,jj) = glamu(ji,jj)
(  226)    
(  227)                gphit(ji,jj) = gphi0 + ppe2_m * 1.e-3 * ( FLOAT( jj - 1 + njmpp - 1 )       )
(  228)                gphiu(ji,jj) = gphit(ji,jj)
(  229)                gphiv(ji,jj) = gphi0 + ppe2_m * 1.e-3 * ( FLOAT( jj - 1 + njmpp - 1 ) + 0.5 )
(  230)                gphif(ji,jj) = gphiv(ji,jj)
(  231)             END DO
(  232)          END DO
(  233) 
(  234)          ! Horizontal scale factors (in meters)
(  235)          !                              ======
(  236)          e1t(:,:) = ppe1_m      ;      e2t(:,:) = ppe2_m
(  237)          e1u(:,:) = ppe1_m      ;      e2u(:,:) = ppe2_m
(  238)          e1v(:,:) = ppe1_m      ;      e2v(:,:) = ppe2_m
(  239)          e1f(:,:) = ppe1_m      ;      e2f(:,:) = ppe2_m
(  240) 
(  241)       CASE ( 4 )                     ! geographical mesh on the sphere, isotropic MERCATOR type
(  242) 
(  243)          IF(lwp) WRITE(numout,*)
(  244)          IF(lwp) WRITE(numout,*) '          geographical mesh on the sphere, MERCATOR type'
(  245)          IF(lwp) WRITE(numout,*) '          longitudinal/latitudinal spacing given by ppe1_deg'
(  246)          IF ( ppgphi0 == -90 ) THEN
(  247)                 IF(lwp) WRITE(numout,*) ' Mercator grid cannot start at south pole !!!! '
(  248)                 IF(lwp) WRITE(numout,*) ' We stop '
(  249)                 STOP
(  250)          ENDIF
(  251) 
(  252)          !  Find index corresponding to the equator, given the grid spacing e1_deg
(  253)          !  and the (approximate) southern latitude ppgphi0.
(  254)          !  This way we ensure that the equator is at a "T / U" point, when in the domain.
(  255)          !  The formula should work even if the equator is outside the domain.
(  256)          zarg = rpi / 4. - rpi / 180. * ppgphi0 / 2.
(  257)          ijeq = ABS( 180./rpi * LOG( COS( zarg ) / SIN( zarg ) ) / ppe1_deg )
(  258)          IF(  ppgphi0 > 0 )  ijeq = -ijeq
(  259) 
(  260)          IF(lwp) WRITE(numout,*) '          Index of the equator on the MERCATOR grid:', ijeq
(  261) 
(  262)          DO jj = 1, jpj
(  263)             DO ji = 1, jpi
(  264)                zti = FLOAT( ji - 1 + nimpp - 1 )         ;   ztj = FLOAT( jj - ijeq + njmpp - 1 )
(  265)                zui = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zuj = FLOAT( jj - ijeq + njmpp - 1 )
(  266)                zvi = FLOAT( ji - 1 + nimpp - 1 )         ;   zvj = FLOAT( jj - ijeq + njmpp - 1 ) + 0.5
(  267)                zfi = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zfj = FLOAT( jj - ijeq + njmpp - 1 ) + 0.5
(  268)          ! Longitude
(  269)                glamt(ji,jj) = ppglam0 + ppe1_deg * zti
(  270)                glamu(ji,jj) = ppglam0 + ppe1_deg * zui
(  271)                glamv(ji,jj) = ppglam0 + ppe1_deg * zvi
(  272)                glamf(ji,jj) = ppglam0 + ppe1_deg * zfi
(  273)          ! Latitude
(  274)                gphit(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* ztj ) )
(  275)                gphiu(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zuj ) )
(  276)                gphiv(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zvj ) )
(  277)                gphif(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zfj ) )
(  278)          ! e1






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 6

(  279)                e1t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
(  280)                e1u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
(  281)                e1v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
(  282)                e1f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
(  283)          ! e2
(  284)                e2t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
(  285)                e2u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
(  286)                e2v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
(  287)                e2f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
(  288)             END DO
(  289)          END DO
(  290) 
(  291)       CASE DEFAULT
(  292)          IF(lwp) WRITE(numout,cform_err)
(  293)          IF(lwp) WRITE(numout,*) '          bad flag value for jphgr_msh = ', jphgr_msh
(  294)          nstop = nstop + 1
(  295) 
(  296)       END SELECT
(  297) 
(  298) 
(  299)       ! Control printing : Grid informations (if not restart)
(  300)       ! ----------------
(  301) 
(  302)       IF(lwp .AND. .NOT.ln_rstart ) THEN
(  303)          WRITE(numout,*)
(  304)          WRITE(numout,*) '          longitude and e1 scale factors'
(  305)          WRITE(numout,*) '          ------------------------------'
(  306)          WRITE(numout,9300) ( ji, glamt(ji,1), glamu(ji,1),   &
(  307)             glamv(ji,1), glamf(ji,1),   &
(  308)             e1t(ji,1), e1u(ji,1),   &
(  309)             e1v(ji,1), e1f(ji,1), ji = 1, jpi,10)
(  310) 9300     FORMAT( 1x, i4, f8.2,1x, f8.2,1x, f8.2,1x, f8.2, 1x,    &
(  311)             f19.10, 1x, f19.10, 1x, f19.10, 1x, f19.10 )
(  312)          
(  313)          WRITE(numout,*)
(  314)          WRITE(numout,*) '          latitude and e2 scale factors'
(  315)          WRITE(numout,*) '          -----------------------------'
(  316)          WRITE(numout,9300) ( jj, gphit(1,jj), gphiu(1,jj),   &
(  317)             &                     gphiv(1,jj), gphif(1,jj),   &
(  318)             &                     e2t  (1,jj), e2u  (1,jj),   &
(  319)             &                     e2v  (1,jj), e2f  (1,jj), jj = 1, jpj, 10 )
(  320)       ENDIF
(  321) 
(  322)       
(  323)       IF( nprint == 1 .AND. lwp ) THEN
(  324)          WRITE(numout,*) '          e1u e2u '
(  325)          CALL prihre( e1u,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  326)          CALL prihre( e2u,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  327)          WRITE(numout,*) '          e1v e2v  '
(  328)          CALL prihre( e1v,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  329)          CALL prihre( e2v,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  330)          WRITE(numout,*) '          e1f e2f  '
(  331)          CALL prihre( e1f,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  332)          CALL prihre( e2f,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  333)       ENDIF
(  334) 
(  335) 
(  336)       ! ================= !






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 7

(  337)       !  Coriolis factor  !
(  338)       ! ================= !
(  339) 
(  340)       SELECT CASE( jphgr_msh )   ! type of horizontal mesh
(  341) 
(  342)       CASE ( 0, 1, 4 )               ! mesh on the sphere
(  343) 
(  344)          ff(:,:) = 2. * omega * SIN( rad * gphif(:,:) ) 
(  345) 
(  346)       CASE ( 2 )                     ! f-plane at ppgphi0 
(  347) 
(  348)          ff(:,:) = 2. * omega * SIN( rad * ppgphi0 )
(  349) 
(  350)          IF(lwp) WRITE(numout,*) '          f-plane: Coriolis parameter = constant = ', ff(1,1)
(  351) 
(  352)       CASE ( 3 )                     ! beta-plane
(  353) 
(  354)          zbeta   = 2. * omega * COS( rad * ppgphi0 ) / ra                       ! beta at latitude ppgphi0
(  355)          zphi0   = ppgphi0 - FLOAT( jpjglo/2) * ppe2_m / ( ra * rad )           ! latitude of the first row F-points
(  356)          
(  357) #if defined key_agrif && defined key_eel_r6
(  358)          IF (.Not.Agrif_Root()) THEN
(  359)            zphi0 = ppgphi0 - FLOAT( Agrif_Parent(jpjglo)/2)*Agrif_Parent(ppe2_m) / (ra * rad)
(  360)          ENDIF
(  361) #endif         
(  362)          zf0     = 2. * omega * SIN( rad * zphi0 )                              ! compute f0 1st point south
(  363) 
(  364)          ff(:,:) = ( zf0  + zbeta * gphif(:,:) * 1.e+3 )                        ! f = f0 +beta* y ( y=0 at south)
(  365)          
(  366)          IF(lwp) WRITE(numout,*) 
(  367)          IF(lwp) WRITE(numout,*) ' Beta-plane: Beta parameter = constant = ', ff(nldi,nldj)
(  368)          IF(lwp) WRITE(numout,*) ' Coriolis parameter varies from ', ff(nldi,nldj),' to ', ff(nldi,nlej)
(  369)          IF( lk_mpp ) THEN 
(  370)             zminff=ff(nldi,nldj)
(  371)             zmaxff=ff(nldi,nlej)
(  372)             CALL mpp_min( zminff )   ! min over the global domain
(  373)             CALL mpp_max( zmaxff )   ! max over the global domain
(  374)             IF(lwp) WRITE(numout,*) ' Coriolis parameter varies globally from ', zminff,' to ', zmaxff
(  375)          END IF
(  376) 
(  377) 
(  378)       END SELECT
(  379) 
(  380) 
(  381)       ! Control of domain for symetrical condition
(  382)       ! ------------------------------------------
(  383)       ! The equator line must be the latitude coordinate axe
(  384) 
(  385)       IF( nperio == 2 ) THEN
(  386)          znorme = SQRT( SUM( gphiu(:,2) * gphiu(:,2) ) ) / FLOAT( jpi )
(  387)          IF( znorme > 1.e-13 ) THEN
(  388)             IF(lwp) WRITE(numout,cform_err)
(  389)             IF(lwp) WRITE(numout,*) ' ===>>>> : symmetrical condition: rerun with good equator line'
(  390)             nstop = nstop + 1
(  391)          ENDIF
(  392)       ENDIF
(  393) 
(  394)    END SUBROUTINE dom_hgr






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 8

(  395) 
(  396) 
(  397)    SUBROUTINE hgr_read
(  398)       !!---------------------------------------------------------------------
(  399)       !!              ***  ROUTINE hgr_read  ***
(  400)       !!
(  401)       !! ** Purpose :   Read a coordinate file in NetCDF format 
(  402)       !!
(  403)       !! ** Method  :   The mesh file has been defined trough a analytical 
(  404)       !!      or semi-analytical method. It is read in a NetCDF file. 
(  405)       !!     
(  406)       !! References :
(  407)       !!      Marti, Madec and Delecluse, 1992, JGR, 97, 12,763-12,766.
(  408)       !!      Madec, Imbard, 1996, Clim. Dyn., 12, 381-388.
(  409)       !!
(  410)       !! History :
(  411)       !!        !         (O. Marti)  Original code
(  412)       !!        !  91-03  (G. Madec)
(  413)       !!        !  92-07  (M. Imbard)
(  414)       !!        !  99-11  (M. Imbard) NetCDF format with IOIPSL
(  415)       !!        !  00-08  (D. Ludicone) Reduced section at Bab el Mandeb
(  416)       !!   8.5  !  02-06  (G. Madec)  F90: Free form
(  417)       !!----------------------------------------------------------------------
(  418)       !! * Modules used
(  419) !!DB
(  420)       USE lib_ncdf
(  421) 
(  422) 
(  423)       !! * Local declarations
(  424)       LOGICAL ::   llog = .FALSE.
(  425)       CHARACTER(len=21) ::   clname
(  426)       INTEGER  ::   ji, jj              ! dummy loop indices
(  427)       INTEGER  ::   inum                ! temporary logical unit
(  428)       INTEGER  ::   ilev, itime         ! temporary integers
(  429)       REAL(wp) ::   zdt, zdate0         ! temporary scalars
(  430)       REAL(wp) ::   zdept(1)            ! temporary workspace
(  431)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(  432)          zlamt, zphit, zdta             ! temporary workspace (NetCDF read)
(  433) !!DB
(  434)       INTEGER :: status
(  435)       REAL(wp), DIMENSION(jpidta,jpjdta,1) ::   data_3d
(  436) 
(  437)       !!----------------------------------------------------------------------
(  438) 
(  439) !!DB
(  440)       clname = 'coordinates.nc'
(  441) 
(  442) 
(  443) #if defined key_agrif
(  444)       if ( .NOT. Agrif_Root() ) then
(  445)          clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
(  446)       endif
(  447) #endif         
(  448) 
(  449) 
(  450)       ! 1. Read of the grid coordinates and scale factors
(  451)       ! -------------------------------------------------
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 9

(  453)       IF(lwp) THEN
(  454)          WRITE(numout,*)
(  455)          WRITE(numout,*) 'hgr_read : read the horizontal coordinates'
(  456)          WRITE(numout,*) '~~~~~~~~~~~      jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpk = ', jpk
(  457)       ENDIF
(  458) 
(  459)       ! read the file
(  460)       itime = 0
(  461)       ilev = 1   
(  462)       zlamt(:,:) = 0.e0
(  463)       zphit(:,:) = 0.e0
(  464) 
(  465) !!DB: replace IOIPSL with lib_ncdf
(  466)       call ncdf_read_global(clname,'glamt',data_3d,-1,status)
(  467)       zdta(:,:) = data_3d(:,:,1)
(  468)       DO jj = 1, nlcj
(  469)          DO ji = 1, nlci
(  470)             glamt(ji,jj) =  zdta(mig(ji),mjg(jj))
(  471)          END DO
(  472)       END DO
(  473)       call ncdf_read_global(clname,'glamu',data_3d,-1,status)
(  474)       zdta(:,:) = data_3d(:,:,1)
(  475)       DO jj = 1, nlcj
(  476)          DO ji = 1, nlci
(  477)             glamu(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  478)          END DO
(  479)       END DO
(  480)       call ncdf_read_global(clname,'glamv',data_3d,-1,status)
(  481)       zdta(:,:) = data_3d(:,:,1)
(  482)       DO jj = 1, nlcj
(  483)          DO ji = 1, nlci
(  484)             glamv(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  485)          END DO
(  486)       END DO
(  487)       call ncdf_read_global(clname,'glamf',data_3d,-1,status)
(  488)       zdta(:,:) = data_3d(:,:,1)
(  489)       DO jj = 1, nlcj
(  490)          DO ji = 1, nlci
(  491)             glamf(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  492)          END DO
(  493)       END DO
(  494)       call ncdf_read_global(clname,'gphit',data_3d,-1,status)
(  495)       zdta(:,:) = data_3d(:,:,1)
(  496)       DO jj = 1, nlcj
(  497)          DO ji = 1, nlci
(  498)             gphit(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  499)          END DO
(  500)       END DO
(  501) !!DB: Initialize coriolis arrays along OBC
(  502) !!Uses f at T-pt
(  503) !!REM:   ff(:,:) = 2. * omega * SIN( rad * gphif(:,:) ) 
(  504)       do ji = 1, jpidta
(  505)          ff_s(ji) = 2. * omega * sin(rad*zdta(ji,1))
(  506)          ff_n(ji) = 2. * omega * sin(rad*zdta(ji,jpjdta))  !!for future use
(  507)       enddo
(  508)       do jj = 1, jpjdta
(  509)          ff_e(jj) = 2. * omega * sin(rad*zdta(jpidta,jj))
(  510)          ff_w(jj) =  2. * omega * sin(rad*zdta(1,jj))






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 10

(  511)       enddo
(  512) 
(  513)       call ncdf_read_global(clname,'gphiu',data_3d,-1,status)
(  514)       zdta(:,:) = data_3d(:,:,1)
(  515)       DO jj = 1, nlcj
(  516)          DO ji = 1, nlci
(  517)             gphiu(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  518)          END DO
(  519)       END DO
(  520)       call ncdf_read_global(clname,'gphiv',data_3d,-1,status)
(  521)       zdta(:,:) = data_3d(:,:,1)
(  522)       DO jj = 1, nlcj
(  523)          DO ji = 1, nlci
(  524)             gphiv(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  525)          END DO
(  526)       END DO
(  527)       call ncdf_read_global(clname,'gphif',data_3d,-1,status)
(  528)       zdta(:,:) = data_3d(:,:,1)
(  529)       DO jj = 1, nlcj
(  530)          DO ji = 1, nlci
(  531)             gphif(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  532)          END DO
(  533)       END DO
(  534) 
(  535)       call ncdf_read_global(clname,'e1t',data_3d,-1,status)
(  536)       zdta(:,:) = data_3d(:,:,1)
(  537)       DO jj = 1, nlcj
(  538)          DO ji = 1, nlci
(  539)             e1t  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  540)          END DO
(  541)       END DO
(  542)       call ncdf_read_global(clname,'e1u',data_3d,-1,status)
(  543)       zdta(:,:) = data_3d(:,:,1)
(  544)       DO jj = 1, nlcj
(  545)          DO ji = 1, nlci
(  546)             e1u  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  547)          END DO
(  548)       END DO
(  549)       call ncdf_read_global(clname,'e1v',data_3d,-1,status)
(  550)       zdta(:,:) = data_3d(:,:,1)
(  551)       DO jj = 1, nlcj
(  552)          DO ji = 1, nlci
(  553)             e1v  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  554)          END DO
(  555)       END DO
(  556)       call ncdf_read_global(clname,'e1f',data_3d,-1,status)
(  557)       zdta(:,:) = data_3d(:,:,1)
(  558)       DO jj = 1, nlcj
(  559)          DO ji = 1, nlci
(  560)             e1f  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  561)          END DO
(  562)       END DO
(  563)       call ncdf_read_global(clname,'e2t',data_3d,-1,status)
(  564)       zdta(:,:) = data_3d(:,:,1)
(  565)       DO jj = 1, nlcj
(  566)          DO ji = 1, nlci
(  567)             e2t  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  568)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 11

(  569)       END DO
(  570)       call ncdf_read_global(clname,'e2u',data_3d,-1,status)
(  571)       zdta(:,:) = data_3d(:,:,1)
(  572)       DO jj = 1, nlcj
(  573)          DO ji = 1, nlci
(  574)             e2u  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  575)          END DO
(  576)       END DO
(  577) !!DB: 2007.12.11 -- Initialize boundary array
(  578)       do ji = 1, jpidta
(  579)          e2u_s(ji) = zdta(ji,1)
(  580)          e2u_n(ji) = zdta(ji,jpjdta)
(  581)       enddo
(  582) 
(  583)       call ncdf_read_global(clname,'e2v',data_3d,-1,status)
(  584)       zdta(:,:) = data_3d(:,:,1)
(  585)       DO jj = 1, nlcj
(  586)          DO ji = 1, nlci
(  587)             e2v  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  588)          END DO
(  589)       END DO
(  590) !!DB: 2007.12.11 -- Initialize boundary arrays
(  591) !!   Not sure at this time exactly what i I should use, but 2 should be OK
(  592)       do jj = 1, jpjdta
(  593)          e2v_e(jj) = zdta(jpidta,jj)
(  594)          e2v_w(jj) = zdta(2,jj)
(  595)       enddo
(  596) 
(  597)       call ncdf_read_global(clname,'e2f',data_3d,-1,status)
(  598)       zdta(:,:) = data_3d(:,:,1)
(  599)       DO jj = 1, nlcj
(  600)          DO ji = 1, nlci
(  601)             e2f  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  602)          END DO
(  603)       END DO
(  604) 
(  605) 
(  606)       ! set extra rows add in mpp to none zero values
(  607)       DO jj = nlcj+1, jpj
(  608)          DO ji = 1, nlci
(  609)             glamt(ji,jj) = glamt(ji,1)   ;   gphit(ji,jj) = gphit(ji,1)
(  610)             glamu(ji,jj) = glamu(ji,1)   ;   gphiu(ji,jj) = gphiu(ji,1)
(  611)             glamv(ji,jj) = glamv(ji,1)   ;   gphiv(ji,jj) = gphiv(ji,1)
(  612)             glamf(ji,jj) = glamf(ji,1)   ;   gphif(ji,jj) = gphif(ji,1)
(  613)             e1t  (ji,jj) = e1t  (ji,1)   ;   e2t  (ji,jj) = e2t  (ji,1)
(  614)             e1u  (ji,jj) = e1u  (ji,1)   ;   e2u  (ji,jj) = e2u  (ji,1)
(  615)             e1v  (ji,jj) = e1v  (ji,1)   ;   e2v  (ji,jj) = e2v  (ji,1)
(  616)             e1f  (ji,jj) = e1f  (ji,1)   ;   e2f  (ji,jj) = e2f  (ji,1)
(  617)          END DO
(  618)       END DO
(  619) 
(  620)       ! set extra columns add in mpp to none zero values
(  621)       DO ji = nlci+1, jpi
(  622)          glamt(ji,:) = glamt(1,:)   ;   gphit(ji,:) = gphit(1,:)
(  623)          glamu(ji,:) = glamu(1,:)   ;   gphiu(ji,:) = gphiu(1,:)
(  624)          glamv(ji,:) = glamv(1,:)   ;   gphiv(ji,:) = gphiv(1,:)
(  625)          glamf(ji,:) = glamf(1,:)   ;   gphif(ji,:) = gphif(1,:)
(  626)          e1t  (ji,:) = e1t  (1,:)   ;   e2t  (ji,:) = e2t  (1,:)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 12

(  627)          e1u  (ji,:) = e1u  (1,:)   ;   e2u  (ji,:) = e2u  (1,:)
(  628)          e1v  (ji,:) = e1v  (1,:)   ;   e2v  (ji,:) = e2v  (1,:)
(  629)          e1f  (ji,:) = e1f  (1,:)   ;   e2f  (ji,:) = e2f  (1,:)
(  630)       END DO
(  631) 
(  632)    END SUBROUTINE hgr_read
(  633) 
(  634)    !!======================================================================
(  635) END MODULE domhgr























































PGF90 (Version     10.5)          09/23/2010  09:49:32      page 13

(    1) # 1 "domhgr.F90"
(    1) !!DB 2008.09.30 -- Eliminated IOIPSL in hgr_read
(    2) !!   2009.08.31 -- Eliminated GYRE (jphgr_msh=5)
(    3) MODULE domhgr
(    4)    !!==============================================================================
(    5)    !!                       ***  MODULE domhgr   ***
(    6)    !! Ocean initialization : domain initialization
(    7)    !!==============================================================================
(    8) 
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dom_hgr        : initialize the horizontal mesh 
(   11)    !!   hgr_read       : read "coordinate" NetCDF file 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE phycst          ! physical constants
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE lib_mpp
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Module variables
(   23)    REAL(wp) ::   glam0, gphi0           ! variables corresponding to parameters
(   24)       !                                 ! ppglam0 ppgphi0 set in par_oce
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC dom_hgr        ! called by domain.F90
(   28)    !!----------------------------------------------------------------------
(   29)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domhgr.F90,v 1.15 2006/04/10 15:46:07 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36)    SUBROUTINE dom_hgr
(   37)       !!----------------------------------------------------------------------
(   38)       !!                  ***  ROUTINE dom_hgr  ***
(   39)       !!
(   40)       !! ** Purpose :   Compute the geographical position (in degre) of the 
(   41)       !!      model grid-points,  the horizontal scale factors (in meters) and 
(   42)       !!      the Coriolis factor (in s-1).
(   43)       !!
(   44)       !! ** Method  :   The geographical position of the model grid-points is
(   45)       !!      defined from analytical functions, fslam and fsphi, the deriva-
(   46)       !!      tives of which gives the horizontal scale factors e1,e2.
(   47)       !!      Defining two function fslam and fsphi and their derivatives in 
(   48)       !!      the two horizontal directions (fse1 and fse2), the model grid-
(   49)       !!      point position and scale factors are given by:
(   50)       !!         t-point:
(   51)       !!      glamt(i,j) = fslam(i    ,j    )   e1t(i,j) = fse1(i    ,j    )
(   52)       !!      gphit(i,j) = fsphi(i    ,j    )   e2t(i,j) = fse2(i    ,j    )
(   53)       !!         u-point:
(   54)       !!      glamu(i,j) = fslam(i+1/2,j    )   e1u(i,j) = fse1(i+1/2,j    )
(   55)       !!      gphiu(i,j) = fsphi(i+1/2,j    )   e2u(i,j) = fse2(i+1/2,j    )
(   56)       !!         v-point:
(   57)       !!      glamv(i,j) = fslam(i    ,j+1/2)   e1v(i,j) = fse1(i    ,j+1/2)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 14

(   58)       !!      gphiv(i,j) = fsphi(i    ,j+1/2)   e2v(i,j) = fse2(i    ,j+1/2)
(   59)       !!            f-point:
(   60)       !!      glamf(i,j) = fslam(i+1/2,j+1/2)   e1f(i,j) = fse1(i+1/2,j+1/2)
(   61)       !!      gphif(i,j) = fsphi(i+1/2,j+1/2)   e2f(i,j) = fse2(i+1/2,j+1/2)
(   62)       !!      Where fse1 and fse2 are defined by:
(   63)       !!         fse1(i,j) = ra * rad * SQRT( (cos(phi) di(fslam))**2
(   64)       !!                                     +          di(fsphi) **2 )(i,j)
(   65)       !!         fse2(i,j) = ra * rad * SQRT( (cos(phi) dj(fslam))**2
(   66)       !!                                     +          dj(fsphi) **2 )(i,j)
(   67)       !!
(   68)       !!        The coriolis factor is given at z-point by:
(   69)       !!                     ff = 2.*omega*sin(gphif)      (in s-1)
(   70)       !!
(   71)       !!        This routine is given as an example, it must be modified
(   72)       !!      following the user s desiderata. nevertheless, the output as
(   73)       !!      well as the way to compute the model grid-point position and
(   74)       !!      horizontal scale factors must be respected in order to insure
(   75)       !!      second order accuracy schemes.
(   76)       !!
(   77)       !! N.B. If the domain is periodic, verify that scale factors are also
(   78)       !!      periodic, and the coriolis term again.
(   79)       !!
(   80)       !! ** Action  : - define  glamt, glamu, glamv, glamf: longitude of t-, 
(   81)       !!                u-, v- and f-points (in degre)
(   82)       !!              - define  gphit, gphiu, gphiv, gphit: latitude  of t-,
(   83)       !!               u-, v-  and f-points (in degre)
(   84)       !!        define e1t, e2t, e1u, e2u, e1v, e2v, e1f, e2f: horizontal
(   85)       !!      scale factors (in meters) at t-, u-, v-, and f-points.
(   86)       !!        define ff: coriolis factor at f-point
(   87)       !!
(   88)       !! References :
(   89)       !!      Marti, Madec and Delecluse, 1992, j. geophys. res., in press.
(   90)       !!
(   91)       !! History :
(   92)       !!        !  88-03  (G. Madec)
(   93)       !!        !  91-11  (G. Madec)
(   94)       !!        !  92-06  (M. Imbard)
(   95)       !!        !  96-01  (G. Madec)  terrain following coordinates
(   96)       !!        !  97-02  (G. Madec)  print mesh informations
(   97)       !!        !  01-09  (M. Levy)  eel config: grid in km, beta-plane
(   98)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module, namelist
(   99)       !!   9.0  !  04-01  (A.M. Treguier, J.M. Molines) Case 4 (Mercator mesh)
(  100)       !!                  use of parameters in par_CONFIG-Rxx.h90, not in namelist
(  101)       !!        !  04-05  (A. Koch-Larrouy) Add Gyre configuration 
(  102)       !!----------------------------------------------------------------------
(  103)       !! * local declarations
(  104)       INTEGER  ::   ji, jj              ! dummy loop indices
(  105)       INTEGER  ::   ii0, ii1, ij0, ij1  ! temporary integers
(  106)       INTEGER  ::   ijeq                ! index of equator T point (used in case 4)
(  107)       REAL(wp) ::   &
(  108)          zti, zui, zvi, zfi,         &  ! temporary scalars
(  109)          ztj, zuj, zvj, zfj,         &  !
(  110)          zphi0, zbeta, znorme,       &  !
(  111)          zarg, zf0, zminff, zmaxff
(  112)       REAL(wp) ::   &
(  113)          zlam1, zcos_alpha, zim1 , zjm1 , ze1, ze1deg,   &
(  114)          zphi1, zsin_alpha, zim05, zjm05
(  115)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 15

(  116) 
(  117)       IF(lwp) THEN
(  118)          WRITE(numout,*)
(  119)          WRITE(numout,*) 'dom_hgr : define the horizontal mesh from ithe following par_oce parameters '
(  120)          WRITE(numout,*) '~~~~~~~      type of horizontal mesh           jphgr_msh = ', jphgr_msh
(  121)          WRITE(numout,*) '             position of the first row and     ppglam0  = ', ppglam0
(  122)          WRITE(numout,*) '             column grid-point (degrees)       ppgphi0  = ', ppgphi0
(  123)          WRITE(numout,*) '             zonal      grid-spacing (degrees) ppe1_deg = ', ppe1_deg
(  124)          WRITE(numout,*) '             meridional grid-spacing (degrees) ppe2_deg = ', ppe2_deg
(  125)          WRITE(numout,*) '             zonal      grid-spacing (meters)  ppe1_m   = ', ppe1_m  
(  126)          WRITE(numout,*) '             meridional grid-spacing (meters)  ppe2_m   = ', ppe2_m  
(  127)       ENDIF
(  128) 
(  129) 
(  130)       SELECT CASE( jphgr_msh )   ! type of horizontal mesh
(  131) 
(  132)       CASE ( 0 )                     !  curvilinear coordinate on the sphere read in coordinate.nc file
(  133) 
(  134)          IF(lwp) WRITE(numout,*)
(  135)          IF(lwp) WRITE(numout,*) '          curvilinear coordinate on the sphere read in "coordinate" file'
(  136) 
(  137)          CALL hgr_read           ! Defaultl option  :   NetCDF file
(  138) 
(  139)          !                                                ! =====================
(  140) !!DB: delete ORCA
(  141) !         IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN    ! ORCA R2 configuration
(  142) 
(  143)          !                                                ! ======================
(  144) !         IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN   ! ORCA R05 configuration
(  145)          ! N.B. :  General case, lat and long function of both i and j indices:
(  146)          !     e1t(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphit(ji,jj) ) * fsdila( zti, ztj ) )**2   &
(  147)          !                                  + (                           fsdiph( zti, ztj ) )**2  )
(  148)          !     e1u(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiu(ji,jj) ) * fsdila( zui, zuj ) )**2   &
(  149)          !                                  + (                           fsdiph( zui, zuj ) )**2  )
(  150)          !     e1v(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiv(ji,jj) ) * fsdila( zvi, zvj ) )**2   &
(  151)          !                                  + (                           fsdiph( zvi, zvj ) )**2  )
(  152)          !     e1f(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphif(ji,jj) ) * fsdila( zfi, zfj ) )**2   &
(  153)          !                                  + (                           fsdiph( zfi, zfj ) )**2  )
(  154)          !
(  155)          !     e2t(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphit(ji,jj) ) * fsdjla( zti, ztj ) )**2   &
(  156)          !                                  + (                           fsdjph( zti, ztj ) )**2  )
(  157)          !     e2u(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiu(ji,jj) ) * fsdjla( zui, zuj ) )**2   &
(  158)          !                                  + (                           fsdjph( zui, zuj ) )**2  )
(  159)          !     e2v(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiv(ji,jj) ) * fsdjla( zvi, zvj ) )**2   &
(  160)          !                                  + (                           fsdjph( zvi, zvj ) )**2  )
(  161)          !     e2f(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphif(ji,jj) ) * fsdjla( zfi, zfj ) )**2   &
(  162)          !                                  + (                           fsdjph( zfi, zfj ) )**2  )
(  163) 
(  164) 
(  165)       CASE ( 1 )                     ! geographical mesh on the sphere with regular grid-spacing
(  166) 
(  167)          IF(lwp) WRITE(numout,*)
(  168)          IF(lwp) WRITE(numout,*) '          geographical mesh on the sphere with regular grid-spacing'
(  169)          IF(lwp) WRITE(numout,*) '          given by ppe1_deg and ppe2_deg' 
(  170) 
(  171)          DO jj = 1, jpj
(  172)             DO ji = 1, jpi
(  173)                zti = FLOAT( ji - 1 + nimpp - 1 )         ;   ztj = FLOAT( jj - 1 + njmpp - 1 )






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 16

(  174)                zui = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zuj = FLOAT( jj - 1 + njmpp - 1 )
(  175)                zvi = FLOAT( ji - 1 + nimpp - 1 )         ;   zvj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5
(  176)                zfi = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zfj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5
(  177)          ! Longitude
(  178)                glamt(ji,jj) = ppglam0 + ppe1_deg * zti
(  179)                glamu(ji,jj) = ppglam0 + ppe1_deg * zui
(  180)                glamv(ji,jj) = ppglam0 + ppe1_deg * zvi
(  181)                glamf(ji,jj) = ppglam0 + ppe1_deg * zfi
(  182)          ! Latitude
(  183)                gphit(ji,jj) = ppgphi0 + ppe2_deg * ztj
(  184)                gphiu(ji,jj) = ppgphi0 + ppe2_deg * zuj
(  185)                gphiv(ji,jj) = ppgphi0 + ppe2_deg * zvj
(  186)                gphif(ji,jj) = ppgphi0 + ppe2_deg * zfj
(  187)          ! e1
(  188)                e1t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
(  189)                e1u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
(  190)                e1v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
(  191)                e1f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
(  192)          ! e2
(  193)                e2t(ji,jj) = ra * rad * ppe2_deg
(  194)                e2u(ji,jj) = ra * rad * ppe2_deg
(  195)                e2v(ji,jj) = ra * rad * ppe2_deg
(  196)                e2f(ji,jj) = ra * rad * ppe2_deg
(  197)             END DO
(  198)          END DO
(  199) 
(  200) 
(  201)       CASE ( 2:3 )                   ! f- or beta-plane with regular grid-spacing
(  202) 
(  203)          IF(lwp) WRITE(numout,*)
(  204)          IF(lwp) WRITE(numout,*) '          f- or beta-plane with regular grid-spacing'
(  205)          IF(lwp) WRITE(numout,*) '          given by ppe1_m and ppe2_m' 
(  206) 
(  207)          ! Position coordinates (in kilometers)
(  208)          !                          ==========
(  209)          glam0 = 0.e0
(  210)          gphi0 = - ppe2_m * 1.e-3
(  211)          
(  212) # 220
(  220)          DO jj = 1, jpj
(  221)             DO ji = 1, jpi
(  222)                glamt(ji,jj) = glam0 + ppe1_m * 1.e-3 * ( FLOAT( ji - 1 + nimpp - 1 )       )
(  223)                glamu(ji,jj) = glam0 + ppe1_m * 1.e-3 * ( FLOAT( ji - 1 + nimpp - 1 ) + 0.5 )
(  224)                glamv(ji,jj) = glamt(ji,jj)
(  225)                glamf(ji,jj) = glamu(ji,jj)
(  226)    
(  227)                gphit(ji,jj) = gphi0 + ppe2_m * 1.e-3 * ( FLOAT( jj - 1 + njmpp - 1 )       )
(  228)                gphiu(ji,jj) = gphit(ji,jj)
(  229)                gphiv(ji,jj) = gphi0 + ppe2_m * 1.e-3 * ( FLOAT( jj - 1 + njmpp - 1 ) + 0.5 )
(  230)                gphif(ji,jj) = gphiv(ji,jj)
(  231)             END DO
(  232)          END DO
(  233) 
(  234)          ! Horizontal scale factors (in meters)
(  235)          !                              ======
(  236)          e1t(:,:) = ppe1_m      ;      e2t(:,:) = ppe2_m
(  237)          e1u(:,:) = ppe1_m      ;      e2u(:,:) = ppe2_m
(  238)          e1v(:,:) = ppe1_m      ;      e2v(:,:) = ppe2_m






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 17

(  239)          e1f(:,:) = ppe1_m      ;      e2f(:,:) = ppe2_m
(  240) 
(  241)       CASE ( 4 )                     ! geographical mesh on the sphere, isotropic MERCATOR type
(  242) 
(  243)          IF(lwp) WRITE(numout,*)
(  244)          IF(lwp) WRITE(numout,*) '          geographical mesh on the sphere, MERCATOR type'
(  245)          IF(lwp) WRITE(numout,*) '          longitudinal/latitudinal spacing given by ppe1_deg'
(  246)          IF ( ppgphi0 == -90 ) THEN
(  247)                 IF(lwp) WRITE(numout,*) ' Mercator grid cannot start at south pole !!!! '
(  248)                 IF(lwp) WRITE(numout,*) ' We stop '
(  249)                 STOP
(  250)          ENDIF
(  251) 
(  252)          !  Find index corresponding to the equator, given the grid spacing e1_deg
(  253)          !  and the (approximate) southern latitude ppgphi0.
(  254)          !  This way we ensure that the equator is at a "T / U" point, when in the domain.
(  255)          !  The formula should work even if the equator is outside the domain.
(  256)          zarg = rpi / 4. - rpi / 180. * ppgphi0 / 2.
(  257)          ijeq = ABS( 180./rpi * LOG( COS( zarg ) / SIN( zarg ) ) / ppe1_deg )
(  258)          IF(  ppgphi0 > 0 )  ijeq = -ijeq
(  259) 
(  260)          IF(lwp) WRITE(numout,*) '          Index of the equator on the MERCATOR grid:', ijeq
(  261) 
(  262)          DO jj = 1, jpj
(  263)             DO ji = 1, jpi
(  264)                zti = FLOAT( ji - 1 + nimpp - 1 )         ;   ztj = FLOAT( jj - ijeq + njmpp - 1 )
(  265)                zui = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zuj = FLOAT( jj - ijeq + njmpp - 1 )
(  266)                zvi = FLOAT( ji - 1 + nimpp - 1 )         ;   zvj = FLOAT( jj - ijeq + njmpp - 1 ) + 0.5
(  267)                zfi = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zfj = FLOAT( jj - ijeq + njmpp - 1 ) + 0.5
(  268)          ! Longitude
(  269)                glamt(ji,jj) = ppglam0 + ppe1_deg * zti
(  270)                glamu(ji,jj) = ppglam0 + ppe1_deg * zui
(  271)                glamv(ji,jj) = ppglam0 + ppe1_deg * zvi
(  272)                glamf(ji,jj) = ppglam0 + ppe1_deg * zfi
(  273)          ! Latitude
(  274)                gphit(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* ztj ) )
(  275)                gphiu(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zuj ) )
(  276)                gphiv(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zvj ) )
(  277)                gphif(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zfj ) )
(  278)          ! e1
(  279)                e1t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
(  280)                e1u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
(  281)                e1v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
(  282)                e1f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
(  283)          ! e2
(  284)                e2t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
(  285)                e2u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
(  286)                e2v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
(  287)                e2f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
(  288)             END DO
(  289)          END DO
(  290) 
(  291)       CASE DEFAULT
(  292)          IF(lwp) WRITE(numout,cform_err)
(  293)          IF(lwp) WRITE(numout,*) '          bad flag value for jphgr_msh = ', jphgr_msh
(  294)          nstop = nstop + 1
(  295) 
(  296)       END SELECT






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 18

(  297) 
(  298) 
(  299)       ! Control printing : Grid informations (if not restart)
(  300)       ! ----------------
(  301) 
(  302)       IF(lwp .AND. .NOT.ln_rstart ) THEN
(  303)          WRITE(numout,*)
(  304)          WRITE(numout,*) '          longitude and e1 scale factors'
(  305)          WRITE(numout,*) '          ------------------------------'
(  306)          WRITE(numout,9300) ( ji, glamt(ji,1), glamu(ji,1),   &
(  307)             glamv(ji,1), glamf(ji,1),   &
(  308)             e1t(ji,1), e1u(ji,1),   &
(  309)             e1v(ji,1), e1f(ji,1), ji = 1, jpi,10)
(  310) 9300     FORMAT( 1x, i4, f8.2,1x, f8.2,1x, f8.2,1x, f8.2, 1x,    &
(  311)             f19.10, 1x, f19.10, 1x, f19.10, 1x, f19.10 )
(  312)          
(  313)          WRITE(numout,*)
(  314)          WRITE(numout,*) '          latitude and e2 scale factors'
(  315)          WRITE(numout,*) '          -----------------------------'
(  316)          WRITE(numout,9300) ( jj, gphit(1,jj), gphiu(1,jj),   &
(  317)             &                     gphiv(1,jj), gphif(1,jj),   &
(  318)             &                     e2t  (1,jj), e2u  (1,jj),   &
(  319)             &                     e2v  (1,jj), e2f  (1,jj), jj = 1, jpj, 10 )
(  320)       ENDIF
(  321) 
(  322)       
(  323)       IF( nprint == 1 .AND. lwp ) THEN
(  324)          WRITE(numout,*) '          e1u e2u '
(  325)          CALL prihre( e1u,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  326)          CALL prihre( e2u,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  327)          WRITE(numout,*) '          e1v e2v  '
(  328)          CALL prihre( e1v,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  329)          CALL prihre( e2v,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  330)          WRITE(numout,*) '          e1f e2f  '
(  331)          CALL prihre( e1f,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  332)          CALL prihre( e2f,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
(  333)       ENDIF
(  334) 
(  335) 
(  336)       ! ================= !
(  337)       !  Coriolis factor  !
(  338)       ! ================= !
(  339) 
(  340)       SELECT CASE( jphgr_msh )   ! type of horizontal mesh
(  341) 
(  342)       CASE ( 0, 1, 4 )               ! mesh on the sphere
(  343) 
(  344)          ff(:,:) = 2. * omega * SIN( rad * gphif(:,:) ) 
(  345) 
(  346)       CASE ( 2 )                     ! f-plane at ppgphi0 
(  347) 
(  348)          ff(:,:) = 2. * omega * SIN( rad * ppgphi0 )
(  349) 
(  350)          IF(lwp) WRITE(numout,*) '          f-plane: Coriolis parameter = constant = ', ff(1,1)
(  351) 
(  352)       CASE ( 3 )                     ! beta-plane
(  353) 
(  354)          zbeta   = 2. * omega * COS( rad * ppgphi0 ) / ra                       ! beta at latitude ppgphi0






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 19

(  355)          zphi0   = ppgphi0 - FLOAT( jpjglo/2) * ppe2_m / ( ra * rad )           ! latitude of the first row F-points
(  356)          
(  357) # 362
(  362)          zf0     = 2. * omega * SIN( rad * zphi0 )                              ! compute f0 1st point south
(  363) 
(  364)          ff(:,:) = ( zf0  + zbeta * gphif(:,:) * 1.e+3 )                        ! f = f0 +beta* y ( y=0 at south)
(  365)          
(  366)          IF(lwp) WRITE(numout,*) 
(  367)          IF(lwp) WRITE(numout,*) ' Beta-plane: Beta parameter = constant = ', ff(nldi,nldj)
(  368)          IF(lwp) WRITE(numout,*) ' Coriolis parameter varies from ', ff(nldi,nldj),' to ', ff(nldi,nlej)
(  369)          IF( lk_mpp ) THEN 
(  370)             zminff=ff(nldi,nldj)
(  371)             zmaxff=ff(nldi,nlej)
(  372)             CALL mpp_min( zminff )   ! min over the global domain
(  373)             CALL mpp_max( zmaxff )   ! max over the global domain
(  374)             IF(lwp) WRITE(numout,*) ' Coriolis parameter varies globally from ', zminff,' to ', zmaxff
(  375)          END IF
(  376) 
(  377) 
(  378)       END SELECT
(  379) 
(  380) 
(  381)       ! Control of domain for symetrical condition
(  382)       ! ------------------------------------------
(  383)       ! The equator line must be the latitude coordinate axe
(  384) 
(  385)       IF( nperio == 2 ) THEN
(  386)          znorme = SQRT( SUM( gphiu(:,2) * gphiu(:,2) ) ) / FLOAT( jpi )
(  387)          IF( znorme > 1.e-13 ) THEN
(  388)             IF(lwp) WRITE(numout,cform_err)
(  389)             IF(lwp) WRITE(numout,*) ' ===>>>> : symmetrical condition: rerun with good equator line'
(  390)             nstop = nstop + 1
(  391)          ENDIF
(  392)       ENDIF
(  393) 
(  394)    END SUBROUTINE dom_hgr
(  395) 
(  396) 
(  397)    SUBROUTINE hgr_read
(  398)       !!---------------------------------------------------------------------
(  399)       !!              ***  ROUTINE hgr_read  ***
(  400)       !!
(  401)       !! ** Purpose :   Read a coordinate file in NetCDF format 
(  402)       !!
(  403)       !! ** Method  :   The mesh file has been defined trough a analytical 
(  404)       !!      or semi-analytical method. It is read in a NetCDF file. 
(  405)       !!     
(  406)       !! References :
(  407)       !!      Marti, Madec and Delecluse, 1992, JGR, 97, 12,763-12,766.
(  408)       !!      Madec, Imbard, 1996, Clim. Dyn., 12, 381-388.
(  409)       !!
(  410)       !! History :
(  411)       !!        !         (O. Marti)  Original code
(  412)       !!        !  91-03  (G. Madec)
(  413)       !!        !  92-07  (M. Imbard)
(  414)       !!        !  99-11  (M. Imbard) NetCDF format with IOIPSL
(  415)       !!        !  00-08  (D. Ludicone) Reduced section at Bab el Mandeb
(  416)       !!   8.5  !  02-06  (G. Madec)  F90: Free form






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 20

(  417)       !!----------------------------------------------------------------------
(  418)       !! * Modules used
(  419) !!DB
(  420)       USE lib_ncdf
(  421) 
(  422) 
(  423)       !! * Local declarations
(  424)       LOGICAL ::   llog = .FALSE.
(  425)       CHARACTER(len=21) ::   clname
(  426)       INTEGER  ::   ji, jj              ! dummy loop indices
(  427)       INTEGER  ::   inum                ! temporary logical unit
(  428)       INTEGER  ::   ilev, itime         ! temporary integers
(  429)       REAL(wp) ::   zdt, zdate0         ! temporary scalars
(  430)       REAL(wp) ::   zdept(1)            ! temporary workspace
(  431)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(  432)          zlamt, zphit, zdta             ! temporary workspace (NetCDF read)
(  433) !!DB
(  434)       INTEGER :: status
(  435)       REAL(wp), DIMENSION(jpidta,jpjdta,1) ::   data_3d
(  436) 
(  437)       !!----------------------------------------------------------------------
(  438) 
(  439) !!DB
(  440)       clname = 'coordinates.nc'
(  441) 
(  442) 
(  443) 
(  444) 
(  445) # 450
(  450)       ! 1. Read of the grid coordinates and scale factors
(  451)       ! -------------------------------------------------
(  452) 
(  453)       IF(lwp) THEN
(  454)          WRITE(numout,*)
(  455)          WRITE(numout,*) 'hgr_read : read the horizontal coordinates'
(  456)          WRITE(numout,*) '~~~~~~~~~~~      jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpk = ', jpk
(  457)       ENDIF
(  458) 
(  459)       ! read the file
(  460)       itime = 0
(  461)       ilev = 1   
(  462)       zlamt(:,:) = 0.e0
(  463)       zphit(:,:) = 0.e0
(  464) 
(  465) !!DB: replace IOIPSL with lib_ncdf
(  466)       call ncdf_read_global(clname,'glamt',data_3d,-1,status)
(  467)       zdta(:,:) = data_3d(:,:,1)
(  468)       DO jj = 1, nlcj
(  469)          DO ji = 1, nlci
(  470)             glamt(ji,jj) =  zdta(mig(ji),mjg(jj))
(  471)          END DO
(  472)       END DO
(  473)       call ncdf_read_global(clname,'glamu',data_3d,-1,status)
(  474)       zdta(:,:) = data_3d(:,:,1)
(  475)       DO jj = 1, nlcj
(  476)          DO ji = 1, nlci
(  477)             glamu(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  478)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 21

(  479)       END DO
(  480)       call ncdf_read_global(clname,'glamv',data_3d,-1,status)
(  481)       zdta(:,:) = data_3d(:,:,1)
(  482)       DO jj = 1, nlcj
(  483)          DO ji = 1, nlci
(  484)             glamv(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  485)          END DO
(  486)       END DO
(  487)       call ncdf_read_global(clname,'glamf',data_3d,-1,status)
(  488)       zdta(:,:) = data_3d(:,:,1)
(  489)       DO jj = 1, nlcj
(  490)          DO ji = 1, nlci
(  491)             glamf(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  492)          END DO
(  493)       END DO
(  494)       call ncdf_read_global(clname,'gphit',data_3d,-1,status)
(  495)       zdta(:,:) = data_3d(:,:,1)
(  496)       DO jj = 1, nlcj
(  497)          DO ji = 1, nlci
(  498)             gphit(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  499)          END DO
(  500)       END DO
(  501) !!DB: Initialize coriolis arrays along OBC
(  502) !!Uses f at T-pt
(  503) !!REM:   ff(:,:) = 2. * omega * SIN( rad * gphif(:,:) ) 
(  504)       do ji = 1, jpidta
(  505)          ff_s(ji) = 2. * omega * sin(rad*zdta(ji,1))
(  506)          ff_n(ji) = 2. * omega * sin(rad*zdta(ji,jpjdta))  !!for future use
(  507)       enddo
(  508)       do jj = 1, jpjdta
(  509)          ff_e(jj) = 2. * omega * sin(rad*zdta(jpidta,jj))
(  510)          ff_w(jj) =  2. * omega * sin(rad*zdta(1,jj))
(  511)       enddo
(  512) 
(  513)       call ncdf_read_global(clname,'gphiu',data_3d,-1,status)
(  514)       zdta(:,:) = data_3d(:,:,1)
(  515)       DO jj = 1, nlcj
(  516)          DO ji = 1, nlci
(  517)             gphiu(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  518)          END DO
(  519)       END DO
(  520)       call ncdf_read_global(clname,'gphiv',data_3d,-1,status)
(  521)       zdta(:,:) = data_3d(:,:,1)
(  522)       DO jj = 1, nlcj
(  523)          DO ji = 1, nlci
(  524)             gphiv(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  525)          END DO
(  526)       END DO
(  527)       call ncdf_read_global(clname,'gphif',data_3d,-1,status)
(  528)       zdta(:,:) = data_3d(:,:,1)
(  529)       DO jj = 1, nlcj
(  530)          DO ji = 1, nlci
(  531)             gphif(ji,jj) = zdta(mig(ji),mjg(jj))                    
(  532)          END DO
(  533)       END DO
(  534) 
(  535)       call ncdf_read_global(clname,'e1t',data_3d,-1,status)
(  536)       zdta(:,:) = data_3d(:,:,1)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 22

(  537)       DO jj = 1, nlcj
(  538)          DO ji = 1, nlci
(  539)             e1t  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  540)          END DO
(  541)       END DO
(  542)       call ncdf_read_global(clname,'e1u',data_3d,-1,status)
(  543)       zdta(:,:) = data_3d(:,:,1)
(  544)       DO jj = 1, nlcj
(  545)          DO ji = 1, nlci
(  546)             e1u  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  547)          END DO
(  548)       END DO
(  549)       call ncdf_read_global(clname,'e1v',data_3d,-1,status)
(  550)       zdta(:,:) = data_3d(:,:,1)
(  551)       DO jj = 1, nlcj
(  552)          DO ji = 1, nlci
(  553)             e1v  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  554)          END DO
(  555)       END DO
(  556)       call ncdf_read_global(clname,'e1f',data_3d,-1,status)
(  557)       zdta(:,:) = data_3d(:,:,1)
(  558)       DO jj = 1, nlcj
(  559)          DO ji = 1, nlci
(  560)             e1f  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  561)          END DO
(  562)       END DO
(  563)       call ncdf_read_global(clname,'e2t',data_3d,-1,status)
(  564)       zdta(:,:) = data_3d(:,:,1)
(  565)       DO jj = 1, nlcj
(  566)          DO ji = 1, nlci
(  567)             e2t  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  568)          END DO
(  569)       END DO
(  570)       call ncdf_read_global(clname,'e2u',data_3d,-1,status)
(  571)       zdta(:,:) = data_3d(:,:,1)
(  572)       DO jj = 1, nlcj
(  573)          DO ji = 1, nlci
(  574)             e2u  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  575)          END DO
(  576)       END DO
(  577) !!DB: 2007.12.11 -- Initialize boundary array
(  578)       do ji = 1, jpidta
(  579)          e2u_s(ji) = zdta(ji,1)
(  580)          e2u_n(ji) = zdta(ji,jpjdta)
(  581)       enddo
(  582) 
(  583)       call ncdf_read_global(clname,'e2v',data_3d,-1,status)
(  584)       zdta(:,:) = data_3d(:,:,1)
(  585)       DO jj = 1, nlcj
(  586)          DO ji = 1, nlci
(  587)             e2v  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  588)          END DO
(  589)       END DO
(  590) !!DB: 2007.12.11 -- Initialize boundary arrays
(  591) !!   Not sure at this time exactly what i I should use, but 2 should be OK
(  592)       do jj = 1, jpjdta
(  593)          e2v_e(jj) = zdta(jpidta,jj)
(  594)          e2v_w(jj) = zdta(2,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 23

(  595)       enddo
(  596) 
(  597)       call ncdf_read_global(clname,'e2f',data_3d,-1,status)
(  598)       zdta(:,:) = data_3d(:,:,1)
(  599)       DO jj = 1, nlcj
(  600)          DO ji = 1, nlci
(  601)             e2f  (ji,jj) = zdta(mig(ji),mjg(jj))                    
(  602)          END DO
(  603)       END DO
(  604) 
(  605) 
(  606)       ! set extra rows add in mpp to none zero values
(  607)       DO jj = nlcj+1, jpj
(  608)          DO ji = 1, nlci
(  609)             glamt(ji,jj) = glamt(ji,1)   ;   gphit(ji,jj) = gphit(ji,1)
(  610)             glamu(ji,jj) = glamu(ji,1)   ;   gphiu(ji,jj) = gphiu(ji,1)
(  611)             glamv(ji,jj) = glamv(ji,1)   ;   gphiv(ji,jj) = gphiv(ji,1)
(  612)             glamf(ji,jj) = glamf(ji,1)   ;   gphif(ji,jj) = gphif(ji,1)
(  613)             e1t  (ji,jj) = e1t  (ji,1)   ;   e2t  (ji,jj) = e2t  (ji,1)
(  614)             e1u  (ji,jj) = e1u  (ji,1)   ;   e2u  (ji,jj) = e2u  (ji,1)
(  615)             e1v  (ji,jj) = e1v  (ji,1)   ;   e2v  (ji,jj) = e2v  (ji,1)
(  616)             e1f  (ji,jj) = e1f  (ji,1)   ;   e2f  (ji,jj) = e2f  (ji,1)
(  617)          END DO
(  618)       END DO
(  619) 
(  620)       ! set extra columns add in mpp to none zero values
(  621)       DO ji = nlci+1, jpi
(  622)          glamt(ji,:) = glamt(1,:)   ;   gphit(ji,:) = gphit(1,:)
(  623)          glamu(ji,:) = glamu(1,:)   ;   gphiu(ji,:) = gphiu(1,:)
(  624)          glamv(ji,:) = glamv(1,:)   ;   gphiv(ji,:) = gphiv(1,:)
(  625)          glamf(ji,:) = glamf(1,:)   ;   gphif(ji,:) = gphif(1,:)
(  626)          e1t  (ji,:) = e1t  (1,:)   ;   e2t  (ji,:) = e2t  (1,:)
(  627)          e1u  (ji,:) = e1u  (1,:)   ;   e2u  (ji,:) = e2u  (1,:)
(  628)          e1v  (ji,:) = e1v  (1,:)   ;   e2v  (ji,:) = e2v  (1,:)
(  629)          e1f  (ji,:) = e1f  (1,:)   ;   e2f  (ji,:) = e2f  (1,:)
(  630)       END DO
(  631) 
(  632)    END SUBROUTINE hgr_read
(  633) 
(  634)    !!======================================================================
(  635) END MODULE domhgr
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (domhgr.F90: 173)
