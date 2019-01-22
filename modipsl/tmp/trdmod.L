


PGF90 (Version     10.5)          09/23/2010  09:49:10      page 1

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

Filename: trdmod.F90

(    1) MODULE trdmod
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trdmod  ***
(    4)    !! Ocean diagnostics:  ocean tracers and dynamic trends
(    5)    !!=====================================================================
(    6) #if  defined key_trdtra || defined key_trddyn || defined key_trdmld || defined key_trdvor || defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trd_mod          : Call the trend to be computed
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce                     ! ocean dynamics and tracers variables
(   12)    USE dom_oce                 ! ocean space and time domain variables
(   13)    USE trdmod_oce              ! ocean variables trends
(   14)    USE trdvor                  ! ocean vorticity trends 
(   15)    USE trdicp                  ! ocean bassin integral constraints properties
(   16)    USE trdmld                  ! ocean active mixed layer tracers trends 
(   17)    USE trabbl                  ! bottom boundary layer variables
(   18)    USE in_out_manager          ! I/O manager
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Routine accessibility
(   24)    PUBLIC trd_mod        ! called by all dynXX or traXX modules
(   25) 
(   26)    !! * Substitutions
(   27) #  include "domzgr_substitute.h90"
(   28) #  include "vectopt_loop_substitute.h90"
(   29)    !!----------------------------------------------------------------------
(   30)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmod.F90,v 1.2 2005/03/27 18:35:24 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35) CONTAINS
(   36) 
(   37)    SUBROUTINE trd_mod(ptrdx, ptrdy, ktrd, ctype, kt)
(   38)       !!---------------------------------------------------------------------
(   39)       !!                  ***  ROUTINE trd_mod  ***
(   40)       !! 
(   41)       !! ** Purpose : Dispatch all trends computation, e.g. vorticity, mld or 
(   42)       !!              integral constrains
(   43)       !!
(   44)       !! ** Method :
(   45)       !!
(   46)       !! History :






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 2

(   47)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   48)       !!----------------------------------------------------------------------
(   49)       !! * Modules used
(   50) #if defined key_trabbl_adv
(   51)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &  ! temporary arrays
(   52)          &         zun, zvn
(   53) #else
(   54)       USE oce                , zun => un,  &  ! When no bbl, zun == un
(   55)          &                     zvn => vn      ! When no bbl, zvn == vn
(   56) #endif
(   57) 
(   58)       !! * Arguments
(   59)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
(   60)          ptrdx,                      &   ! Temperature or U trend 
(   61)          ptrdy                           ! Salinity    or V trend
(   62) 
(   63)       INTEGER, INTENT( in ) ::   &
(   64)          kt  ,                   & ! time step
(   65)          ktrd                      ! tracer trend index
(   66) 
(   67)       CHARACTER(len=3), INTENT( in ) ::   &
(   68)          ctype                             ! momentum or tracers trends type
(   69)          !                                 ! 'DYN' or 'TRA'
(   70) 
(   71)       !! * Local save
(   72)       REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
(   73)          zbtr2
(   74) 
(   75)       !! * Local declarations
(   76)       INTEGER ::   ji, jj, jk    ! loop indices
(   77)       REAL(wp) ::   &
(   78)          zbtr,            &  ! temporary scalars
(   79)          zfui, zfvj,           &  !    "         "
(   80)          zfui1, zfvj1             !    "         "
(   81)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   82)          z2dx, z2dy                        ! workspace arrays
(   83)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   84)          z3dx, z3dy                            ! workspace arrays
(   85)       !!----------------------------------------------------------------------
(   86) 
(   87)       ! Initialization of workspace arrays
(   88)       z3dx(:,:,:) = 0.e0
(   89)       z3dy(:,:,:) = 0.e0
(   90)       z2dx(:,:) = 0.e0
(   91)       z2dy(:,:) = 0.e0
(   92) 
(   93)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(   94)       ! I. Bassin averaged properties for momentum and/or tracers trends
(   95)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(   96) 
(   97)       IF( ( mod(kt,ntrd) == 0 .OR. kt == nit000 .OR. kt == nitend) )   THEN
(   98) 
(   99)          ! Active tracers trends 
(  100)          IF( lk_trdtra .AND. ctype == 'TRA' )   THEN
(  101) 
(  102)             IF( ktrd == jpttdnsr )   THEN
(  103)                ! 2D array tracers surface forcing
(  104)                z2dx(:,:) = ptrdx(:,:,1)






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 3

(  105)                z2dy(:,:) = ptrdy(:,:,1)
(  106) 
(  107)                CALL trd(z2dx, z2dy, ktrd, ctype)
(  108)             ELSE
(  109)                ! 3D array
(  110)                CALL trd(ptrdx, ptrdy, ktrd, ctype)
(  111)             ENDIF
(  112) 
(  113)          ENDIF
(  114) 
(  115)          ! Momentum trends 
(  116)          IF( lk_trddyn .AND. ctype == 'DYN' )   THEN
(  117) 
(  118)             IF( ktrd == jpdtdswf .OR. ktrd == jpdtdbfr )   THEN
(  119)                ! momentum surface forcing/bottom friction  2D array
(  120)                z2dx(:,:) = ptrdx(:,:,1)
(  121)                z2dy(:,:) = ptrdy(:,:,1)
(  122) 
(  123)                CALL trd(z2dx, z2dy, ktrd, ctype)
(  124)             ELSE
(  125)                ! 3D array
(  126)                CALL trd(ptrdx, ptrdy, ktrd, ctype)
(  127)             ENDIF
(  128) 
(  129)          ENDIF
(  130) 
(  131)       ENDIF
(  132) 
(  133)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  134)       ! II. Vorticity trends
(  135)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  136) 
(  137)       IF( lk_trdvor .AND. ctype == 'DYN' )   THEN
(  138) 
(  139)          SELECT CASE ( ktrd )
(  140) 
(  141)          ! Pressure Gradient trend
(  142)          CASE ( jpdtdhpg )      
(  143)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorprg)
(  144) 
(  145)          ! KE Gradient trend
(  146)          CASE ( jpdtdkeg )      
(  147)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorkeg)
(  148) 
(  149)          ! Relative Vorticity trend
(  150)          CASE ( jpdtdrvo )      
(  151)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorrvo)
(  152) 
(  153)          ! Planetary Vorticity Term trend
(  154)          CASE ( jpdtdpvo )      
(  155)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorpvo)
(  156) 
(  157)          ! Horizontal Diffusion trend
(  158)          CASE ( jpdtdldf )      
(  159)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorldf)
(  160) 
(  161)          ! Vertical Advection trend
(  162)          CASE ( jpdtdzad )      






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 4

(  163)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorzad)
(  164) 
(  165)          ! Vertical Diffusion trend
(  166)          CASE ( jpdtdzdf )      
(  167)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorzdf)
(  168) 
(  169)          ! Surface Pressure Grad. trend
(  170)          CASE ( jpdtdspg )      
(  171)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorspg)
(  172) 
(  173)          ! Beta V trend 
(  174)          CASE ( jpdtddat )      
(  175)             CALL trd_vor_zint(ptrdx, ptrdy, jpvorbev)
(  176) 
(  177)          ! Wind stress forcing term
(  178)          CASE ( jpdtdswf )      
(  179)             z2dx(:,:) = ptrdx(:,:,1)
(  180)             z2dy(:,:) = ptrdy(:,:,1)
(  181) 
(  182)             CALL trd_vor_zint(z2dx, z2dy, jpvorswf)
(  183) 
(  184)          ! Bottom friction term
(  185)          CASE ( jpdtdbfr )      
(  186)             z2dx(:,:) = ptrdx(:,:,1)
(  187)             z2dy(:,:) = ptrdy(:,:,1)
(  188) 
(  189)             CALL trd_vor_zint(z2dx, z2dy, jpvorbfr)
(  190) 
(  191)          END SELECT
(  192) 
(  193)       ENDIF
(  194) 
(  195)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  196)       ! III. Mixed layer trends
(  197)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  198) 
(  199)       IF( lk_trdmld .AND. ctype == 'TRA' )   THEN
(  200)          
(  201)          SELECT CASE ( ktrd )
(  202) 
(  203)          ! horizontal advection trends
(  204)          CASE ( jpttdlad )      
(  205) 
(  206) #if defined key_trabbl_adv
(  207)             ! Advective bottom boundary layer 
(  208)             ! -------------------------------
(  209)             zun(:,:,:) = un(:,:,:) - u_bbl(:,:,:)
(  210)             zvn(:,:,:) = vn(:,:,:) - v_bbl(:,:,:)
(  211) #endif
(  212)             IF( kt == nit000 )   zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
(  213) 
(  214)             SELECT CASE ( l_adv )
(  215) 
(  216)             CASE ( 'ce2' )
(  217) 
(  218)                ! Split horizontal trends into i- and j- compnents for trdmld case 
(  219)                ! ----------------------------------------------------------------
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 5

(  221)                ! i- advective trend computed as Uh gradh(T)
(  222)                DO jk = 1, jpkm1
(  223)                   DO jj = 2, jpjm1
(  224)                      DO ji = fs_2, fs_jpim1   ! vector opt.
(  225) # if defined key_s_coord || defined key_partial_steps
(  226)                         zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  227) 
(  228)                         zfui = 0.5 * e2u(ji  ,jj) * fse3u(ji,  jj,jk) * zun(ji,  jj,jk)
(  229)                         zfui1= 0.5 * e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * zun(ji-1,jj,jk)
(  230) # else         
(  231)                         zbtr = zbtr2(ji,jj)
(  232) 
(  233)                         zfui = 0.5 * e2u(ji  ,jj) * zun(ji,  jj,jk)
(  234)                         zfui1= 0.5 * e2u(ji-1,jj) * zun(ji-1,jj,jk)
(  235) # endif
(  236)                         ! save i- advective trend 
(  237)                         z3dx(ji,jj,jk) = - zbtr * ( zfui  * ( tn(ji+1,jj,jk) - tn(ji  ,jj,jk) )    &
(  238)                             &                     + zfui1 * ( tn(ji  ,jj,jk) - tn(ji-1,jj,jk) ) )
(  239)                         z3dy(ji,jj,jk) = - zbtr * ( zfui  * ( sn(ji+1,jj,jk) - sn(ji  ,jj,jk) )    &
(  240)                             &                     + zfui1 * ( sn(ji  ,jj,jk) - sn(ji-1,jj,jk) ) )
(  241)                      END DO
(  242)                   END DO
(  243)                END DO
(  244) 
(  245)                ! save the i- horizontal trends for diagnostic
(  246)                CALL trd_mld_zint(z3dx, z3dy, jpmldxad, '3D')
(  247) 
(  248)                ! j- advective trend computed as Uh gradh(T)
(  249)                DO jk = 1, jpkm1
(  250)                   DO jj = 2, jpjm1
(  251)                      DO ji = fs_2, fs_jpim1   ! vector opt.
(  252) # if defined key_s_coord || defined key_partial_steps
(  253)                         zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  254) 
(  255)                         zfvj = 0.5 * e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * zvn(ji,jj  ,jk)
(  256)                         zfvj1= 0.5 * e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * zvn(ji,jj-1,jk)
(  257) # else         
(  258)                         zbtr = zbtr2(ji,jj)
(  259) 
(  260)                         zfvj = 0.5 * e1v(ji,jj  ) * zvn(ji,jj  ,jk)
(  261)                         zfvj1= 0.5 * e1v(ji,jj-1) * zvn(ji,jj-1,jk)
(  262) # endif
(  263)                         ! save j- advective trend 
(  264)                         z3dx(ji,jj,jk) = - zbtr * ( zfvj  * ( tn(ji,jj+1,jk) - tn(ji,jj  ,jk) )   &
(  265)                             &                     + zfvj1 * ( tn(ji,jj  ,jk) - tn(ji,jj-1,jk) ) )
(  266)                         z3dy(ji,jj,jk) = - zbtr * ( zfvj  * ( sn(ji,jj+1,jk) - sn(ji,jj  ,jk) )   &
(  267)                             &                     + zfvj1 * ( sn(ji,jj  ,jk) - sn(ji,jj-1,jk) ) )
(  268)                      END DO
(  269)                   END DO
(  270)                END DO
(  271) 
(  272)                ! save the j- horizontal trend for diagnostic
(  273)                CALL trd_mld_zint(z3dx, z3dy, jpmldyad, '3D')
(  274) 
(  275)             CASE ( 'tvd' )
(  276) 
(  277)                ! Recompute the horizontal advection term Div(Uh.T) term 
(  278)                z3dx(:,:,:) = ptrdx(:,:,:) - tn(:,:,:) * hdivn(:,:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 6

(  279)                z3dy(:,:,:) = ptrdy(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  280) 
(  281)                ! Deduce the i- horizontal advection in substracting the j- one.
(  282)                ! tladj()/sladj() are computed in traadv_tvd.F90 module
(  283)                z3dx(:,:,:) = z3dx(:,:,:) - tladj(:,:,:)
(  284)                z3dy(:,:,:) = z3dy(:,:,:) - sladj(:,:,:)
(  285) 
(  286)                DO jk = 1, jpkm1
(  287)                   DO jj = 2, jpjm1
(  288)                      DO ji = fs_2, fs_jpim1
(  289)                         zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  290) 
(  291)                         ! Compute the zonal et meridional divergence
(  292)                         zfui = e2u(ji  ,jj) * fse3u(ji  ,jj,jk) * zun(ji  ,jj,jk)  &
(  293)                              - e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * zun(ji-1,jj,jk)
(  294)                         zfvj = e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * zvn(ji,jj  ,jk)  &
(  295)                              - e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * zvn(ji,jj-1,jk)
(  296) 
(  297)                         ! i- advective trend computed as U gradx(T/S)
(  298)                         z3dx(ji,jj,jk) = z3dx(ji,jj,jk) + tn(ji,jj,jk) * zfui * zbtr
(  299)                         z3dy(ji,jj,jk) = z3dy(ji,jj,jk) + sn(ji,jj,jk) * zfui * zbtr
(  300) 
(  301)                         ! j- advective trend computed as V grady(T/S)
(  302)                         tladj(ji,jj,jk) = tladj(ji,jj,jk) + tn(ji,jj,jk) * zfvj * zbtr
(  303)                         sladj(ji,jj,jk) = sladj(ji,jj,jk) + sn(ji,jj,jk) * zfvj * zbtr
(  304) 
(  305)                      END DO
(  306)                   END DO
(  307)                END DO
(  308) 
(  309)                ! save the i- horizontal trend for diagnostic
(  310)                CALL trd_mld_zint(z3dx, z3dy, jpmldxad, '3D')
(  311) 
(  312)                ! save the j- horizontal trend for diagnostic
(  313)                CALL trd_mld_zint(tladj, sladi, jpmldyad, '3D')
(  314) 
(  315)             CASE ( 'mus', 'mu2' )
(  316) 
(  317)                !  Split horizontal trends in i- and j- direction for trdmld case 
(  318)                ! ----------------------------------------------------------------
(  319) 
(  320)                ! i- advective trend computed as U gradx(T/S)
(  321)                DO jk = 1, jpkm1
(  322)                   DO jj = 2, jpjm1      
(  323)                      DO ji = fs_2, fs_jpim1   ! vector opt.
(  324) # if defined key_s_coord || defined key_partial_steps
(  325)                         zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  326)                         zfui =  e2u(ji  ,jj) * fse3u(ji,  jj,jk) * zun(ji,  jj,jk)   &
(  327)                            & -  e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * zun(ji-1,jj,jk)
(  328) # else      
(  329)                         zbtr = zbtr2(ji,jj)
(  330)                         zfui = e2u(ji  ,jj) * zun(ji,  jj,jk)   &
(  331)                            & - e2u(ji-1,jj) * zun(ji-1,jj,jk)
(  332) # endif
(  333)                         ! save i- advective trend 
(  334)                         z3dx(ji,jj,jk) = - zbtr * ( tladi(ji,jj,jk) - tladi(ji-1,jj,jk) )   &
(  335)                             &                      + tn(ji,jj,jk) * zfui * zbtr
(  336)                         z3dy(ji,jj,jk) = - zbtr * ( sladi(ji,jj,jk) - sladi(ji-1,jj,jk) )  &






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 7

(  337)                             &                      + sn(ji,jj,jk) * zfui * zbtr
(  338)                      END DO
(  339)                   END DO
(  340)                END DO        
(  341) 
(  342)                ! save the i- horizontal trends for diagnostic
(  343)                CALL trd_mld_zint(z3dx, z3dy, jpmldxad, '3D')
(  344) 
(  345)                ! j- advective trend computed as V grady(T/S)
(  346)                DO jk = 1, jpkm1
(  347)                   DO jj = 2, jpjm1      
(  348)                      DO ji = fs_2, fs_jpim1   ! vector opt.
(  349) # if defined key_s_coord || defined key_partial_steps
(  350)                         zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  351)                         zfvj =  e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * zvn(ji,jj  ,jk)   &
(  352)                            & -  e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * zvn(ji,jj-1,jk)
(  353) # else      
(  354)                         zbtr = zbtr2(ji,jj)
(  355)                         zfvj = e1v(ji,jj  ) * zvn(ji,jj  ,jk)   &
(  356)                            & - e1v(ji,jj-1) * zvn(ji,jj-1,jk)
(  357) # endif
(  358)                         ! save j- advective trend 
(  359)                         z3dx(ji,jj,jk) =  - zbtr * ( tladj(ji,jj,jk) - tladj(ji,jj-1,jk) )   &
(  360)                             &                       + tn(ji,jj,jk) * zfvj * zbtr
(  361)                         z3dy(ji,jj,jk) =  - zbtr * ( sladj(ji,jj,jk) - sladj(ji,jj-1,jk) )   &
(  362)                             &                       + sn(ji,jj,jk) * zfvj * zbtr
(  363)                      END DO
(  364)                   END DO
(  365)                END DO        
(  366) 
(  367)                ! save the j- horizontal trends for diagnostic
(  368)                CALL trd_mld_zint(z3dx, z3dy, jpmldyad, '3D')
(  369) 
(  370)             END SELECT
(  371) 
(  372)          ! vertical advection trends
(  373)          CASE ( jpttdzad )      
(  374)             CALL trd_mld_zint(ptrdx, ptrdy, jpmldzad, '3D')
(  375) 
(  376)          ! lateral diffusion trends
(  377)          CASE ( jpttdldf )      
(  378)             CALL trd_mld_zint(ptrdx, ptrdy, jpmldldf, '3D')
(  379) # if defined key_traldf_eiv
(  380)             ! Save the i- and j- eddy induce velocity trends
(  381)             CALL trd_mld_zint(tladi, sladi, jpmldxei, '3D')
(  382)             CALL trd_mld_zint(tladj, sladj, jpmldyei, '3D')
(  383) # endif
(  384)             IF( lk_trabbl_dif )   THEN
(  385)                z3dx(:,:,:) = 0.e0
(  386)                z3dy(:,:,:) = 0.e0
(  387)                z3dx(:,:,1) = tldfbbl(:,:)
(  388)                z3dy(:,:,1) = sldfbbl(:,:)
(  389)                CALL trd_mld_zint(z3dx, z3dy, jpmldldf, '2D')
(  390)             ENDIF
(  391) 
(  392)          ! vertical diffusion trends
(  393)          CASE ( jpttdzdf )      
(  394)             CALL trd_mld_zint(ptrdx, ptrdy, jpmldzdf, '3D')






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 8

(  395) 
(  396)          ! vertical diffusion trends
(  397)          CASE ( jpttddoe )      
(  398)             CALL trd_mld_zint(ptrdx, ptrdy, jpmldzei, '3D')
(  399) 
(  400)          ! penetrative solar radiation trends
(  401)          CASE ( jpttdqsr )      
(  402)             CALL trd_mld_zint(ptrdx, ptrdy, jpmldfor, '3D')
(  403) 
(  404)          ! non penetrative solar radiation trends
(  405)          CASE ( jpttdnsr )
(  406)             ptrdx(:,:,2:jpk) = 0.e0
(  407)             ptrdy(:,:,2:jpk) = 0.e0
(  408)             CALL trd_mld_zint(ptrdx, ptrdy, jpmldfor, '2D')
(  409) 
(  410)          END SELECT   
(  411) 
(  412)       ENDIF
(  413) 
(  414) 
(  415)    END SUBROUTINE trd_mod
(  416) 
(  417) #   else
(  418)    !!----------------------------------------------------------------------
(  419)    !!   Default case :                                         Empty module
(  420)    !!----------------------------------------------------------------------
(  421)    USE trdmod_oce      ! ocean variables trends
(  422) 
(  423) CONTAINS
(  424)    SUBROUTINE trd_mod(ptrd3dx, ptrd3dy, ktrd , ctype, kt)       ! Empty routine
(  425)       REAL, DIMENSION(:,:,:), INTENT( in ) ::   &
(  426)           ptrd3dx,                     &   ! Temperature or U trend 
(  427)           ptrd3dy                          ! Salinity    or V trend
(  428)       INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
(  429)       INTEGER, INTENT( in ) ::   kt        ! Time step
(  430)       CHARACTER(len=3), INTENT( in ) ::   &
(  431)          ctype                             ! momentum or tracers trends type
(  432) !      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd3dx(1,1,1)
(  433) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd3dy(1,1,1)
(  434) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
(  435) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype
(  436) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', kt
(  437)    END SUBROUTINE trd_mod
(  438) #   endif
(  439) 
(  440)    !!======================================================================
(  441) END MODULE trdmod

















PGF90 (Version     10.5)          09/23/2010  09:49:10      page 9

(    1) # 1 "trdmod.F90"
(    1) MODULE trdmod
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trdmod  ***
(    4)    !! Ocean diagnostics:  ocean tracers and dynamic trends
(    5)    !!=====================================================================
(    6) # 418
(  418)    !!----------------------------------------------------------------------
(  419)    !!   Default case :                                         Empty module
(  420)    !!----------------------------------------------------------------------
(  421)    USE trdmod_oce      ! ocean variables trends
(  422) 
(  423) CONTAINS
(  424)    SUBROUTINE trd_mod(ptrd3dx, ptrd3dy, ktrd , ctype, kt)       ! Empty routine
(  425)       REAL, DIMENSION(:,:,:), INTENT( in ) ::   &
(  426)           ptrd3dx,                     &   ! Temperature or U trend 
(  427)           ptrd3dy                          ! Salinity    or V trend
(  428)       INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
(  429)       INTEGER, INTENT( in ) ::   kt        ! Time step
(  430)       CHARACTER(len=3), INTENT( in ) ::   &
(  431)          ctype                             ! momentum or tracers trends type
(  432) !      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd3dx(1,1,1)
(  433) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd3dy(1,1,1)
(  434) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
(  435) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype
(  436) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', kt
(  437)    END SUBROUTINE trd_mod
(  438) 
(  439) # 440
(  440)    !!======================================================================
(  441) END MODULE trdmod
