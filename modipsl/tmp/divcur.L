


PGF90 (Version     10.5)          09/23/2010  09:49:31      page 1

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

Filename: divcur.F90

(    1) MODULE divcur
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  divcur  ***
(    4)    !! Ocean diagnostic variable : horizontal divergence and relative vorticity
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   div_cur    : Compute the horizontal divergence and relative
(    9)    !!                vorticity fields
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE obc_oce         ! ocean lateral open boundary condition
(   16)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC div_cur    ! routine called by step.F90 and istate.F90
(   23) 
(   24)    !! * Substitutions
(   25) #  include "domzgr_substitute.h90"
(   26) #  include "vectopt_loop_substitute.h90"
(   27)    !!----------------------------------------------------------------------
(   28)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   29)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/divcur.F90,v 1.6 2006/03/10 10:55:41 opalod Exp $ 
(   30)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   31)    !!----------------------------------------------------------------------
(   32) 
(   33) CONTAINS
(   34) 
(   35) #if defined key_noslip_accurate
(   36)    !!----------------------------------------------------------------------
(   37)    !!   'key_noslip_accurate'                     2nd order centered scheme
(   38)    !!                                                4th order at the coast
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41)    SUBROUTINE div_cur( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                  ***  ROUTINE div_cur  ***
(   44)       !!
(   45)       !! ** Purpose :   compute the horizontal divergence and the relative
(   46)       !!      vorticity at before and now time-step






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 2

(   47)       !!
(   48)       !! ** Method  : 
(   49)       !!      I.  divergence :
(   50)       !!         - save the divergence computed at the previous time-step
(   51)       !!      (note that the Asselin filter has not been applied on hdivb)
(   52)       !!         - compute the now divergence given by :
(   53)       !!            * s-coordinate ('key_s_coord' defined)
(   54)       !!         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
(   55)       !!         * z-coordinate (default key)
(   56)       !!         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
(   57)       !!         - apply lateral boundary conditions on hdivn 
(   58)       !!      II. vorticity :
(   59)       !!         - save the curl computed at the previous time-step
(   60)       !!            rotb = rotn
(   61)       !!      (note that the Asselin time filter has not been applied to rotb)
(   62)       !!         - compute the now curl in tensorial formalism:
(   63)       !!            rotn = 1/(e1f*e2f) ( di[e2v vn] - dj[e1u un] )
(   64)       !!         - apply lateral boundary conditions on rotn through a call
(   65)       !!      of lbc_lnk routine.
(   66)       !!         - Coastal boundary condition: 'key_noslip_accurate' defined,
(   67)       !!      the no-slip boundary condition is computed using Schchepetkin
(   68)       !!      and O'Brien (1996) scheme (i.e. 4th order at the coast).
(   69)       !!      For example, along east coast, the one-sided finite difference
(   70)       !!      approximation used for di[v] is:
(   71)       !!         di[e2v vn] =  1/(e1f*e2f)
(   72)       !!                    * ( (e2v vn)(i) + (e2v vn)(i-1) + (e2v vn)(i-2) )
(   73)       !!
(   74)       !! ** Action  : - update hdivb, hdivn, the before & now hor. divergence
(   75)       !!              - update rotb , rotn , the before & now rel. vorticity
(   76)       !!
(   77)       !! History :
(   78)       !!   8.2  !  00-03  (G. Madec)  no slip accurate
(   79)       !!   9.0  !  03-08  (G. Madec)  merged of cur and div, free form, F90
(   80)       !!----------------------------------------------------------------------
(   81)       !! * Arguments
(   82)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   83)       
(   84)       !! * Local declarations
(   85)       INTEGER ::   ji, jj, jk     ! dummy loop indices
(   86)       INTEGER ::   ii, ij, jl     ! temporary integer
(   87)       INTEGER ::   ijt, iju       ! temporary integer
(   88)       REAL(wp) ::   zdiv, zdju
(   89)       REAL(wp), DIMENSION(   jpi  ,1:jpj+2) ::   zwu   ! workspace
(   90)       REAL(wp), DIMENSION(-1:jpi+2,  jpj  ) ::   zwv   ! workspace
(   91)       !!----------------------------------------------------------------------
(   92) 
(   93)       IF( kt == nit000 ) THEN
(   94)          IF(lwp) WRITE(numout,*)
(   95)          IF(lwp) WRITE(numout,*) 'div_cur : horizontal velocity divergence and relative vorticity'
(   96)          IF(lwp) WRITE(numout,*) '~~~~~~~   NOT optimal for auto-tasking case'
(   97)       ENDIF
(   98) 
(   99)       !                                                ! ===============
(  100)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  101)          !                                             ! ===============
(  102) 
(  103)          hdivb(:,:,jk) = hdivn(:,:,jk)    ! time swap of div arrays
(  104)          rotb (:,:,jk) = rotn (:,:,jk)    ! time swap of rot arrays






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 3

(  105) 
(  106)          !                                             ! --------
(  107)          ! Horizontal divergence                       !   div
(  108)          !                                             ! --------
(  109)          DO jj = 2, jpjm1
(  110)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  111) #if defined key_s_coord || defined key_partial_steps
(  112)                hdivn(ji,jj,jk) =   &
(  113)                   (  e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk) - e2u(ji-1,jj  )*fse3u(ji-1,jj  ,jk)  * un(ji-1,jj  ,jk)     
(  114)                    + e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk) - e1v(ji  ,jj-1)*fse3v(ji  ,jj-1,jk)  * vn(ji  ,jj-1,jk)  )  
(  115)                   / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  116) #else
(  117)                hdivn(ji,jj,jk) = (  e2u(ji,jj) * un(ji,jj,jk) - e2u(ji-1,jj  ) * un(ji-1,jj  ,jk)      &
(  118)                   &               + e1v(ji,jj) * vn(ji,jj,jk) - e1v(ji  ,jj-1) * vn(ji  ,jj-1,jk)  )   &
(  119)      &            / ( e1t(ji,jj) * e2t(ji,jj) )
(  120) #endif
(  121)             END DO
(  122)          END DO
(  123) 
(  124) #if defined key_obc
(  125)          ! open boundaries (div must be zero behind the open boundary)
(  126)          !  mpp remark: The zeroing of hdivn can probably be extended to 1->jpi/jpj for the correct row/column
(  127)          IF( lp_obc_east  )   hdivn(nie0p1:nie1p1,nje0  :nje1  ,jk) = 0.e0      ! east
(  128)          IF( lp_obc_west  )   hdivn(niw0  :niw1  ,njw0  :njw1  ,jk) = 0.e0      ! west
(  129)          IF( lp_obc_north )   hdivn(nin0  :nin1  ,njn0p1:njn1p1,jk) = 0.e0      ! north
(  130)          IF( lp_obc_south )   hdivn(nis0  :nis1  ,njs0  :njs1  ,jk) = 0.e0      ! south
(  131) #endif         
(  132) #if defined key_agrif
(  133)          if ( .NOT. AGRIF_Root() ) then
(  134)             IF ((nbondi ==  1).OR.(nbondi == 2)) hdivn(nlci-1 , :     ,jk) = 0.e0      ! east
(  135)             IF ((nbondi == -1).OR.(nbondi == 2)) hdivn(2      , :     ,jk) = 0.e0      ! west
(  136)             IF ((nbondj ==  1).OR.(nbondj == 2)) hdivn(:      ,nlcj-1 ,jk) = 0.e0      ! north
(  137)             IF ((nbondj == -1).OR.(nbondj == 2)) hdivn(:      ,2      ,jk) = 0.e0      ! south
(  138)          endif
(  139) #endif       
(  140) 
(  141)          !                                             ! --------
(  142)          ! relative vorticity                          !   rot 
(  143)          !                                             ! --------
(  144)          ! contravariant velocity (extended for lateral b.c.)
(  145)          ! inside the model domain
(  146)          DO jj = 1, jpj
(  147)             DO ji = 1, jpi
(  148)                zwu(ji,jj) = e1u(ji,jj) * un(ji,jj,jk)
(  149)                zwv(ji,jj) = e2v(ji,jj) * vn(ji,jj,jk)
(  150)             END DO  
(  151)          END DO  
(  152)  
(  153)          ! East-West boundary conditions
(  154)          IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6) THEN
(  155)             zwv(  0  ,:) = zwv(jpi-2,:)
(  156)             zwv( -1  ,:) = zwv(jpi-3,:)
(  157)             zwv(jpi+1,:) = zwv(  3  ,:)
(  158)             zwv(jpi+2,:) = zwv(  4  ,:)
(  159)          ELSE
(  160)             zwv(  0  ,:) = 0.e0
(  161)             zwv( -1  ,:) = 0.e0
(  162)             zwv(jpi+1,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 4

(  163)             zwv(jpi+2,:) = 0.e0
(  164)          ENDIF
(  165) 
(  166)          ! North-South boundary conditions
(  167)          IF( nperio == 3 .OR. nperio == 4 ) THEN
(  168)             ! north fold ( Grid defined with a T-point pivot) ORCA 2 degre
(  169)             zwu(jpi,jpj+1) = 0.e0
(  170)             zwu(jpi,jpj+2) = 0.e0
(  171)             DO ji = 1, jpi-1
(  172)                iju = jpi - ji + 1
(  173)                zwu(ji,jpj+1) = - zwu(iju,jpj-3)
(  174)                zwu(ji,jpj+2) = - zwu(iju,jpj-4)
(  175)             END DO
(  176)          ELSEIF( nperio == 5 .OR. nperio == 6 ) THEN
(  177)             ! north fold ( Grid defined with a F-point pivot) ORCA 0.5 degre\
(  178)             zwu(jpi,jpj+1) = 0.e0
(  179)             zwu(jpi,jpj+2) = 0.e0
(  180)             DO ji = 1, jpi-1
(  181)                iju = jpi - ji
(  182)                zwu(ji,jpj  ) = - zwu(iju,jpj-1)
(  183)                zwu(ji,jpj+1) = - zwu(iju,jpj-2)
(  184)                zwu(ji,jpj+2) = - zwu(iju,jpj-3)
(  185)             END DO
(  186)             DO ji = -1, jpi+2
(  187)                ijt = jpi - ji + 1
(  188)                zwv(ji,jpj) = - zwv(ijt,jpj-2)
(  189)             END DO
(  190)             DO ji = jpi/2+1, jpi+2
(  191)                ijt = jpi - ji + 1
(  192)                zwv(ji,jpjm1) = - zwv(ijt,jpjm1)
(  193)             END DO
(  194)          ELSE
(  195)             ! closed
(  196)             zwu(:,jpj+1) = 0.e0
(  197)             zwu(:,jpj+2) = 0.e0
(  198)          ENDIF
(  199) 
(  200)          ! relative vorticity (vertical component of the velocity curl) 
(  201)          DO jj = 1, jpjm1
(  202)             DO ji = 1, fs_jpim1   ! vector opt.
(  203)                rotn(ji,jj,jk) = (  zwv(ji+1,jj  ) - zwv(ji,jj)      &
(  204)                                  - zwu(ji  ,jj+1) + zwu(ji,jj)  )   &
(  205)                               * fmask(ji,jj,jk) / ( e1f(ji,jj)*e2f(ji,jj) )
(  206)             END DO
(  207)          END DO
(  208) 
(  209)          ! second order accurate scheme along straight coast
(  210)          DO jl = 1, npcoa(1,jk)
(  211)             ii = nicoa(jl,1,jk)
(  212)             ij = njcoa(jl,1,jk)
(  213)             rotn(ii,ij,jk) = 1. / ( e1f(ii,ij) * e2f(ii,ij) )   &
(  214)                            * ( + 4. * zwv(ii+1,ij) - zwv(ii+2,ij) + 0.2 * zwv(ii+3,ij) )
(  215)          END DO
(  216)          DO jl = 1, npcoa(2,jk)
(  217)             ii = nicoa(jl,2,jk)
(  218)             ij = njcoa(jl,2,jk)
(  219)             rotn(ii,ij,jk) = 1./(e1f(ii,ij)*e2f(ii,ij))   &
(  220)                *(-4.*zwv(ii,ij)+zwv(ii-1,ij)-0.2*zwv(ii-2,ij))






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 5

(  221)          END DO
(  222)          DO jl = 1, npcoa(3,jk)
(  223)             ii = nicoa(jl,3,jk)
(  224)             ij = njcoa(jl,3,jk)
(  225)             rotn(ii,ij,jk) = -1. / ( e1f(ii,ij)*e2f(ii,ij) )   &
(  226)                * ( +4. * zwu(ii,ij+1) - zwu(ii,ij+2) + 0.2 * zwu(ii,ij+3) )
(  227)          END DO
(  228)          DO jl = 1, npcoa(4,jk)
(  229)             ii = nicoa(jl,4,jk)
(  230)             ij = njcoa(jl,4,jk)
(  231)             rotn(ii,ij,jk) = -1. / ( e1f(ii,ij)*e2f(ii,ij) )   &
(  232)                * ( -4. * zwu(ii,ij) + zwu(ii,ij-1) - 0.2 * zwu(ii,ij-2) )
(  233)          END DO
(  234) 
(  235)          !                                             ! ===============
(  236)       END DO                                           !   End of slab
(  237)       !                                                ! ===============
(  238)       
(  239)       ! 4. Lateral boundary conditions on hdivn and rotn
(  240)       ! ---------------------------------=======---======
(  241)       CALL lbc_lnk( hdivn, 'T', 1. )     ! T-point, no sign change
(  242)       CALL lbc_lnk( rotn , 'F', 1. )     ! F-point, no sign change
(  243) 
(  244)    END SUBROUTINE div_cur
(  245)    
(  246) #else
(  247)    !!----------------------------------------------------------------------
(  248)    !!   Default option                           2nd order centered schemes
(  249)    !!----------------------------------------------------------------------
(  250) 
(  251)    SUBROUTINE div_cur( kt )
(  252)       !!----------------------------------------------------------------------
(  253)       !!                  ***  ROUTINE div_cur  ***
(  254)       !!                    
(  255)       !! ** Purpose :   compute the horizontal divergence and the relative
(  256)       !!      vorticity at before and now time-step
(  257)       !!
(  258)       !! ** Method  : - Divergence:
(  259)       !!      - save the divergence computed at the previous time-step
(  260)       !!      (note that the Asselin filter has not been applied on hdivb)
(  261)       !!      - compute the now divergence given by :
(  262)       !!         * s-coordinate ('key_s_coord' defined)
(  263)       !!         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
(  264)       !!         * z-coordinate (default key)
(  265)       !!         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
(  266)       !!      - apply lateral boundary conditions on hdivn 
(  267)       !!              - Relavtive Vorticity :
(  268)       !!      - save the curl computed at the previous time-step (rotb = rotn)
(  269)       !!      (note that the Asselin time filter has not been applied to rotb)
(  270)       !!      - compute the now curl in tensorial formalism:
(  271)       !!            rotn = 1/(e1f*e2f) ( di[e2v vn] - dj[e1u un] )
(  272)       !!      - apply lateral boundary conditions on rotn through a call to
(  273)       !!      routine lbc_lnk routine.
(  274)       !!      Note: Coastal boundary condition: lateral friction set through
(  275)       !!      the value of fmask along the coast (see dommsk.F90) and shlat
(  276)       !!      (namelist parameter)
(  277)       !!
(  278)       !! ** Action  : - update hdivb, hdivn, the before & now hor. divergence






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 6

(  279)       !!              - update rotb , rotn , the before & now rel. vorticity
(  280)       !!
(  281)       !! History :
(  282)       !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
(  283)       !!   4.0  !  91-11  (G. Madec)
(  284)       !!   6.0  !  93-03  (M. Guyon)  symetrical conditions
(  285)       !!   7.0  !  96-01  (G. Madec)  s-coordinates
(  286)       !!   8.0  !  97-06  (G. Madec)  lateral boundary cond., lbc
(  287)       !!   8.1  !  97-08  (J.M. Molines)  Open boundaries
(  288)       !!   9.0  !  02-09  (G. Madec, E. Durand)  Free form, F90
(  289)       !!----------------------------------------------------------------------
(  290)       !! * Arguments
(  291)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(  292)       
(  293)       !! * Local declarations
(  294)       INTEGER  ::   ji, jj, jk          ! dummy loop indices
(  295)       !!----------------------------------------------------------------------
(  296) 
(  297)       IF( kt == nit000 ) THEN
(  298)          IF(lwp) WRITE(numout,*)
(  299)          IF(lwp) WRITE(numout,*) 'div_cur : horizontal velocity divergence and'
(  300)          IF(lwp) WRITE(numout,*) '~~~~~~~   relative vorticity'
(  301)       ENDIF
(  302) 
(  303)       !                                                ! ===============
(  304)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  305)          !                                             ! ===============
(  306) 
(  307)          hdivb(:,:,jk) = hdivn(:,:,jk)    ! time swap of div arrays
(  308)          rotb (:,:,jk) = rotn (:,:,jk)    ! time swap of rot arrays
(  309) 
(  310)          !                                             ! --------
(  311)          ! Horizontal divergence                       !   div 
(  312)          !                                             ! --------
(  313)          DO jj = 2, jpjm1
(  314)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  315) #if defined key_s_coord || defined key_partial_steps
(  316)                hdivn(ji,jj,jk) =   &
(  317)                   (  e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk) - e2u(ji-1,jj  )*fse3u(ji-1,jj  ,jk)  * un(ji-1,jj  ,jk)     
(  318)                    + e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk) - e1v(ji  ,jj-1)*fse3v(ji  ,jj-1,jk)  * vn(ji  ,jj-1,jk)  )  
(  319)                   / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  320) #else
(  321)                hdivn(ji,jj,jk) = (  e2u(ji,jj) * un(ji,jj,jk) - e2u(ji-1,jj  ) * un(ji-1,jj  ,jk)      &
(  322)                   &               + e1v(ji,jj) * vn(ji,jj,jk) - e1v(ji  ,jj-1) * vn(ji  ,jj-1,jk)  )   & 
(  323)                   / ( e1t(ji,jj) * e2t(ji,jj) )
(  324) #endif
(  325)             END DO  
(  326)          END DO  
(  327) 
(  328) #if defined key_obc
(  329)          ! open boundaries (div must be zero behind the open boundary)
(  330)          !  mpp remark: The zeroing of hdivn can probably be extended to 1->jpi/jpj for the correct row/column
(  331)          IF( lp_obc_east  )   hdivn(nie0p1:nie1p1,nje0  :nje1  ,jk) = 0.e0      ! east
(  332)          IF( lp_obc_west  )   hdivn(niw0  :niw1  ,njw0  :njw1  ,jk) = 0.e0      ! west
(  333)          IF( lp_obc_north )   hdivn(nin0  :nin1  ,njn0p1:njn1p1,jk) = 0.e0      ! north
(  334)          IF( lp_obc_south )   hdivn(nis0  :nis1  ,njs0  :njs1  ,jk) = 0.e0      ! south
(  335) #endif         
(  336) #if defined key_agrif






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 7

(  337)          if ( .NOT. AGRIF_Root() ) then
(  338)             IF ((nbondi ==  1).OR.(nbondi == 2)) hdivn(nlci-1 , :     ,jk) = 0.e0      ! east
(  339)             IF ((nbondi == -1).OR.(nbondi == 2)) hdivn(2      , :     ,jk) = 0.e0      ! west
(  340)             IF ((nbondj ==  1).OR.(nbondj == 2)) hdivn(:      ,nlcj-1 ,jk) = 0.e0      ! north
(  341)             IF ((nbondj == -1).OR.(nbondj == 2)) hdivn(:      ,2      ,jk) = 0.e0      ! south
(  342)          endif
(  343) #endif       
(  344)          !                                             ! --------
(  345)          ! relative vorticity                          !   rot 
(  346)          !                                             ! --------
(  347)          DO jj = 1, jpjm1
(  348)             DO ji = 1, fs_jpim1   ! vector opt.
(  349)                rotn(ji,jj,jk) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ,jk) - e2v(ji,jj) * vn(ji,jj,jk)    &
(  350)                   &              - e1u(ji  ,jj+1) * un(ji  ,jj+1,jk) + e1u(ji,jj) * un(ji,jj,jk)  ) &
(  351)                   &           * fmask(ji,jj,jk) / ( e1f(ji,jj) * e2f(ji,jj) )
(  352)             END DO
(  353)          END DO
(  354)          !                                             ! ===============
(  355)       END DO                                           !   End of slab
(  356)       !                                                ! ===============
(  357)       
(  358)       ! 4. Lateral boundary conditions on hdivn and rotn
(  359)       ! ---------------------------------=======---======
(  360)       CALL lbc_lnk( hdivn, 'T', 1. )       ! T-point, no sign change
(  361)       CALL lbc_lnk( rotn , 'F', 1. )       ! F-point, no sign change
(  362) 
(  363)    END SUBROUTINE div_cur
(  364)    
(  365) #endif
(  366)    !!======================================================================
(  367) END MODULE divcur

































PGF90 (Version     10.5)          09/23/2010  09:49:31      page 8

(    1) # 1 "divcur.F90"
(    1) MODULE divcur
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  divcur  ***
(    4)    !! Ocean diagnostic variable : horizontal divergence and relative vorticity
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   div_cur    : Compute the horizontal divergence and relative
(    9)    !!                vorticity fields
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE obc_oce         ! ocean lateral open boundary condition
(   16)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC div_cur    ! routine called by step.F90 and istate.F90
(   23) 
(   24)    !! * Substitutions
(   25) # 1 "./domzgr_substitute.h90"
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
(   48) # 26 "divcur.F90"
(   26) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 27 "divcur.F90"
(   27) # 27
(   27)    !!----------------------------------------------------------------------
(   28)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   29)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/divcur.F90,v 1.6 2006/03/10 10:55:41 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 9

(   30)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   31)    !!----------------------------------------------------------------------
(   32) 
(   33) CONTAINS
(   34) 
(   35) # 247
(  247)    !!----------------------------------------------------------------------
(  248)    !!   Default option                           2nd order centered schemes
(  249)    !!----------------------------------------------------------------------
(  250) 
(  251)    SUBROUTINE div_cur( kt )
(  252)       !!----------------------------------------------------------------------
(  253)       !!                  ***  ROUTINE div_cur  ***
(  254)       !!                    
(  255)       !! ** Purpose :   compute the horizontal divergence and the relative
(  256)       !!      vorticity at before and now time-step
(  257)       !!
(  258)       !! ** Method  : - Divergence:
(  259)       !!      - save the divergence computed at the previous time-step
(  260)       !!      (note that the Asselin filter has not been applied on hdivb)
(  261)       !!      - compute the now divergence given by :
(  262)       !!         * s-coordinate ('key_s_coord' defined)
(  263)       !!         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
(  264)       !!         * z-coordinate (default key)
(  265)       !!         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
(  266)       !!      - apply lateral boundary conditions on hdivn 
(  267)       !!              - Relavtive Vorticity :
(  268)       !!      - save the curl computed at the previous time-step (rotb = rotn)
(  269)       !!      (note that the Asselin time filter has not been applied to rotb)
(  270)       !!      - compute the now curl in tensorial formalism:
(  271)       !!            rotn = 1/(e1f*e2f) ( di[e2v vn] - dj[e1u un] )
(  272)       !!      - apply lateral boundary conditions on rotn through a call to
(  273)       !!      routine lbc_lnk routine.
(  274)       !!      Note: Coastal boundary condition: lateral friction set through
(  275)       !!      the value of fmask along the coast (see dommsk.F90) and shlat
(  276)       !!      (namelist parameter)
(  277)       !!
(  278)       !! ** Action  : - update hdivb, hdivn, the before & now hor. divergence
(  279)       !!              - update rotb , rotn , the before & now rel. vorticity
(  280)       !!
(  281)       !! History :
(  282)       !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
(  283)       !!   4.0  !  91-11  (G. Madec)
(  284)       !!   6.0  !  93-03  (M. Guyon)  symetrical conditions
(  285)       !!   7.0  !  96-01  (G. Madec)  s-coordinates
(  286)       !!   8.0  !  97-06  (G. Madec)  lateral boundary cond., lbc
(  287)       !!   8.1  !  97-08  (J.M. Molines)  Open boundaries
(  288)       !!   9.0  !  02-09  (G. Madec, E. Durand)  Free form, F90
(  289)       !!----------------------------------------------------------------------
(  290)       !! * Arguments
(  291)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(  292)       
(  293)       !! * Local declarations
(  294)       INTEGER  ::   ji, jj, jk          ! dummy loop indices
(  295)       !!----------------------------------------------------------------------
(  296) 
(  297)       IF( kt == nit000 ) THEN
(  298)          IF(lwp) WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 10

(  299)          IF(lwp) WRITE(numout,*) 'div_cur : horizontal velocity divergence and'
(  300)          IF(lwp) WRITE(numout,*) '~~~~~~~   relative vorticity'
(  301)       ENDIF
(  302) 
(  303)       !                                                ! ===============
(  304)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  305)          !                                             ! ===============
(  306) 
(  307)          hdivb(:,:,jk) = hdivn(:,:,jk)    ! time swap of div arrays
(  308)          rotb (:,:,jk) = rotn (:,:,jk)    ! time swap of rot arrays
(  309) 
(  310)          !                                             ! --------
(  311)          ! Horizontal divergence                       !   div 
(  312)          !                                             ! --------
(  313)          DO jj = 2, jpjm1
(  314)             DO ji = 2, jpim1   ! vector opt.
(  315) # 321
(  321)                hdivn(ji,jj,jk) = (  e2u(ji,jj) * un(ji,jj,jk) - e2u(ji-1,jj  ) * un(ji-1,jj  ,jk)      &
(  322)                   &               + e1v(ji,jj) * vn(ji,jj,jk) - e1v(ji  ,jj-1) * vn(ji  ,jj-1,jk)  )   & 
(  323)                   / ( e1t(ji,jj) * e2t(ji,jj) )
(  324) # 325
(  325)             END DO  
(  326)          END DO  
(  327) 
(  328) # 329
(  329)          ! open boundaries (div must be zero behind the open boundary)
(  330)          !  mpp remark: The zeroing of hdivn can probably be extended to 1->jpi/jpj for the correct row/column
(  331)          IF( lp_obc_east  )   hdivn(nie0p1:nie1p1,nje0  :nje1  ,jk) = 0.e0      ! east
(  332)          IF( lp_obc_west  )   hdivn(niw0  :niw1  ,njw0  :njw1  ,jk) = 0.e0      ! west
(  333)          IF( lp_obc_north )   hdivn(nin0  :nin1  ,njn0p1:njn1p1,jk) = 0.e0      ! north
(  334)          IF( lp_obc_south )   hdivn(nis0  :nis1  ,njs0  :njs1  ,jk) = 0.e0      ! south
(  335) # 344
(  344)          !                                             ! --------
(  345)          ! relative vorticity                          !   rot 
(  346)          !                                             ! --------
(  347)          DO jj = 1, jpjm1
(  348)             DO ji = 1, jpim1   ! vector opt.
(  349)                rotn(ji,jj,jk) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ,jk) - e2v(ji,jj) * vn(ji,jj,jk)    &
(  350)                   &              - e1u(ji  ,jj+1) * un(ji  ,jj+1,jk) + e1u(ji,jj) * un(ji,jj,jk)  ) &
(  351)                   &           * fmask(ji,jj,jk) / ( e1f(ji,jj) * e2f(ji,jj) )
(  352)             END DO
(  353)          END DO
(  354)          !                                             ! ===============
(  355)       END DO                                           !   End of slab
(  356)       !                                                ! ===============
(  357)       
(  358)       ! 4. Lateral boundary conditions on hdivn and rotn
(  359)       ! ---------------------------------=======---======
(  360)       CALL lbc_lnk( hdivn, 'T', 1. )       ! T-point, no sign change
(  361)       CALL lbc_lnk( rotn , 'F', 1. )       ! F-point, no sign change
(  362) 
(  363)    END SUBROUTINE div_cur
(  364)    
(  365) # 366
(  366)    !!======================================================================
(  367) END MODULE divcur
