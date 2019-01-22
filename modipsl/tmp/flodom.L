


PGF90 (Version     10.5)          09/23/2010  09:49:48      page 1

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

Filename: flodom.F90

(    1) !!DB 2008.05.16
(    2) !!Modified routine for more specific DB particle tracking needs
(    3) !!For older code, see OLD_CODE/ or look elsewhere for an older version
(    4) 
(    5) 
(    6) MODULE flodom
(    7)    !!======================================================================
(    8)    !!                       ***  MODULE  flodom  ***
(    9)    !! Ocean floats :   domain
(   10)    !!======================================================================
(   11) #if   defined key_floats   ||   defined key_esopa
(   12)    !!----------------------------------------------------------------------
(   13)    !!   'key_floats'                                     float trajectories
(   14)    !!----------------------------------------------------------------------
(   15)    !!   flo_dom        : initialization of floats
(   16)    !!   findmesh       : compute index of position 
(   17)    !!   dstnce         : compute distance between face mesh and floats 
(   18)    !!----------------------------------------------------------------------
(   19)    !! * Modules used
(   20)    USE oce             ! ocean dynamics and tracers
(   21)    USE dom_oce         ! ocean space and time domain
(   22)    USE flo_oce         ! ocean drifting floats
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE lib_mpp         ! distribued memory computing library
(   25) 
(   26)    IMPLICIT NONE
(   27) 
(   28)    !! * Accessibility
(   29)    PRIVATE  dstnce
(   30)    PUBLIC flo_dom     ! routine called by floats.F90
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/flodom.F90,v 1.5 2005/09/22 10:24:53 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42)   SUBROUTINE flo_dom
(   43)       !! ---------------------------------------------------------------------
(   44)       !!                  ***  ROUTINE flo_dom  ***
(   45)       !!                 
(   46)       !!  ** Purpose :   Initialisation of floats






PGF90 (Version     10.5)          09/23/2010  09:49:48      page 2

(   47)       !!
(   48)       !!  ** Method  :   We put the floats  in the domain with the latitude,
(   49)       !!       the longitude (degree) and the depth (m).
(   50)       !!
(   51)       !!----------------------------------------------------------------------      
(   52)       !! * Local declarations
(   53)      LOGICAL  :: llinmesh
(   54)      CHARACTER (len=21) ::  clname
(   55)      INTEGER  :: ji, jj, jk               ! DO loop index on 3 directions
(   56)      INTEGER  :: jfl, jfl1                ! number of floats   
(   57)      INTEGER  :: inum = 11                ! logical unit for file read
(   58) !!DB
(   59)      INTEGER, PARAMETER :: max_num_floats = 100000
(   60) !     INTEGER, DIMENSION (jpnfl    )  ::   &
(   61)      INTEGER, DIMENSION (:),ALLOCATABLE  ::   &
(   62)           iimfl, ijmfl, ikmfl,    &          ! index mesh of floats
(   63)           idomfl,  ivtest, ihtest
(   64)      REAL(wp) :: zdxab, zdyad
(   65) !!DB
(   66)      INTEGER :: num_floats
(   67)      INTEGER, DIMENSION(:),ALLOCATABLE  ::    &
(   68)           tmp_nisobfl,    &  ! 0 for a isobar float
(   69)           !              ! 1 for a float following the w velocity
(   70)           tmp_ngrpfl         ! number to identify searcher group
(   71)      REAL(wp), DIMENSION(:),ALLOCATABLE ::    &
(   72)           tmp_flxx,       &  ! longitude of float (decimal degree)
(   73)           tmp_flyy,       &  ! latitude of float (decimal degree)
(   74)           tmp_flzz 
(   75) 
(   76)      !!---------------------------------------------------------------------
(   77)      
(   78)      ! Initialisation with the geographical position or restart
(   79)      
(   80)      IF(lwp) WRITE(numout,*) 'flo_dom : compute initial position of floats'
(   81)      IF(lwp) WRITE(numout,*) '~~~~~~~~'
(   82)      IF(lwp) WRITE(numout,*) '           jpnfl determined from init_floats'
(   83)      
(   84)      
(   85)      IF(lwp) WRITE(numout,*) '                     init_float read '
(   86) 
(   87) ! First initialisation of floats
(   88) ! the initials positions of floats are written in a file
(   89) ! with a variable to know if it is a isobar float a number 
(   90) ! to identified who want the trajectories of this float and 
(   91) ! an index for the number of the float         
(   92) ! open the init file 
(   93) !!DB: I keep the above, although it is not clear why
(   94) !!Change so that float file is: flxx  flyy  flzz  nisobfl  ngrpfl -- 1 for each jpnfl floats
(   95)      ALLOCATE(tmp_nisobfl(max_num_floats)); ALLOCATE(tmp_ngrpfl(max_num_floats))
(   96)      ALLOCATE(tmp_flxx(max_num_floats));ALLOCATE(tmp_flyy(max_num_floats))
(   97)      ALLOCATE(tmp_flzz(max_num_floats))
(   98)      
(   99)      clname='init_float'
(  100)      num_floats = 0
(  101)      open(inum,FILE=clname,STATUS='old',ERR=222)   !!Allow for no-file-found
(  102)      do while(num_floats >= 0 .AND. num_floats <= max_num_floats)
(  103)         num_floats = num_floats + 1
(  104)         jfl = num_floats      !lazy






PGF90 (Version     10.5)          09/23/2010  09:49:48      page 3

(  105)         read(inum,*,END=223) tmp_flxx(jfl),tmp_flyy(jfl),tmp_flzz(jfl),tmp_nisobfl(jfl),tmp_ngrpfl(jfl)
(  106)      enddo
(  107) 223  close(inum)
(  108)      num_floats = num_floats - 1
(  109)      if(lwp)write(numout,*)'Found ',num_floats,' floats'
(  110)      goto 224
(  111) !!DB file-not-found
(  112) 222  if(lwp)write(numout,*)'NB: Float file init_floats not found'
(  113) 
(  114) !!DB Allocate arrays here
(  115) 224  if(num_floats == 0) then
(  116) !!DB adjust allocation -- or not as it seems like allocating zero-size arrays is OK
(  117)         jpnfl = 0 ! Force correct value for jpnfl as it is needed elsewhere
(  118)      else
(  119) !!DB allocate arrays to num_floats elements and set jpnfl = num_floats
(  120)         jpnfl = num_floats
(  121)         ALLOCATE(iimfl(jpnfl));   ALLOCATE(ijmfl(jpnfl));   ALLOCATE(ikmfl(jpnfl))
(  122)         ALLOCATE(idomfl(jpnfl));  ALLOCATE(ivtest(jpnfl));  ALLOCATE(ihtest(jpnfl))
(  123)         ALLOCATE(nisobfl(jpnfl)); ALLOCATE(ngrpfl(jpnfl))
(  124)         ALLOCATE(flxx(jpnfl));    ALLOCATE(flyy(jpnfl));    ALLOCATE(flzz(jpnfl))
(  125)         ALLOCATE(tpifl(jpnfl));   ALLOCATE(tpjfl(jpnfl));   ALLOCATE(tpkfl(jpnfl))
(  126) !!Assign working arrays
(  127)         nisobfl(1:jpnfl) = tmp_nisobfl(1:jpnfl); ngrpfl(1:jpnfl) = tmp_ngrpfl(1:jpnfl); 
(  128)         flxx(1:jpnfl) = tmp_flxx(1:jpnfl);flyy(1:jpnfl) = tmp_flyy(1:jpnfl);flzz(1:jpnfl) = tmp_flzz(1:jpnfl);
(  129)      endif
(  130)      DEALLOCATE(tmp_nisobfl); DEALLOCATE(tmp_ngrpfl);
(  131)      DEALLOCATE(tmp_flxx);DEALLOCATE(tmp_flyy);DEALLOCATE(tmp_flzz)
(  132) 
(  133)      
(  134)      ! Test to find the grid point coordonate with the geographical position         
(  135)      DO jfl = 1, jpnfl
(  136) 
(  137)         ihtest(jfl) = 0
(  138)         ivtest(jfl) = 0
(  139)         ikmfl(jfl) = 0
(  140) # if   defined key_mpp_mpi   ||   defined key_mpp_shmem
(  141)         DO ji = MAX(nldi,2), nlei
(  142)            DO jj = MAX(nldj,2), nlej
(  143) # else
(  144)         DO ji = 2, jpi
(  145)            DO jj = 2, jpj
(  146) # endif                  
(  147)                     ! for each float we find the indexes of the mesh 
(  148) 
(  149)                     !!DB: Recalling:
(  150)                     !                               Vij
(  151)                     !                       ---------------Fij
(  152)                     !                       |              |
(  153)                     !                       |              |
(  154)                     !                U i-1,j|       Tij    |Uij
(  155)                     !                       |              |
(  156)                     !                       |              |
(  157)                     !                       ----------------
(  158)                     !                               Vi,j-1
(  159)                     !
(  160)                     ! ====> if .true. then float is in T-cell(ji,jj); ji,jj local indices
(  161) 
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:48      page 4

(  163)               CALL findmesh(glamf(ji-1,jj-1),gphif(ji-1,jj-1),   &
(  164)                    glamf(ji-1,jj  ),gphif(ji-1,jj  ),   &
(  165)                    glamf(ji  ,jj  ),gphif(ji  ,jj  ),   &
(  166)                    glamf(ji  ,jj-1),gphif(ji  ,jj-1),   &
(  167)                    flxx(jfl)       ,flyy(jfl)       ,   &
(  168)                    glamt(ji  ,jj  ),gphit(ji  ,jj  ), llinmesh)
(  169)               IF(llinmesh) THEN
(  170)                  iimfl(jfl)  = ji
(  171)                  ijmfl(jfl)  = jj
(  172)                  ihtest(jfl) = ihtest(jfl)+1
(  173)                  DO jk = 1, jpk-1
(  174)                     IF( (fsdepw(ji,jj,jk) <= flzz(jfl)) .AND. (fsdepw(ji,jj,jk+1) >  flzz(jfl)) ) THEN
(  175)                        ikmfl(jfl)  = jk
(  176)                        ivtest(jfl) = ivtest(jfl) + 1
(  177)                     ENDIF
(  178)                  END DO
(  179)               ENDIF
(  180)            END DO
(  181)         END DO
(  182)         
(  183)         ! If the float is in a mesh computed by an other processor we put iimfl=ijmfl=-1            
(  184)         IF( ihtest(jfl) == 0 ) THEN
(  185)            iimfl(jfl) = -1
(  186)            ijmfl(jfl) = -1
(  187)         ENDIF
(  188)      END DO
(  189)      
(  190)      ! A zero in the sum of the arrays "ihtest" and "ivtest"          
(  191)      IF( lk_mpp )   CALL mpp_sum(ihtest,jpnfl)   ! sums over the global domain
(  192)      IF( lk_mpp )   CALL mpp_sum(ivtest,jpnfl)
(  193)      
(  194)      DO jfl = 1, jpnfl
(  195)         IF( (ihtest(jfl) > 1 ) .OR. ( ivtest(jfl) > 1 )) THEN
(  196)            IF(lwp) WRITE(numout,*) 'THE FLOAT',jfl,' IS NOT IN ONLY ONE MESH'
(  197)         ENDIF
(  198)         IF( ihtest(jfl) == 0 ) THEN 
(  199)            IF(lwp) WRITE(numout,*)'THE FLOAT',jfl,' IS IN NO MESH'
(  200)         ENDIF
(  201)      END DO
(  202)      
(  203)      ! We compute the distance between the float and the face of  the mesh         
(  204)      DO jfl = 1, jpnfl
(  205)         ! Made only if the float is in the domain of the processor
(  206)         IF( (iimfl(jfl) >= 0 ) .AND. ( ijmfl(jfl) >= 0 ) ) THEN
(  207)            
(  208)            ! TEST TO KNOW IF THE FLOAT IS NOT INITIALISED IN THE COAST
(  209)            
(  210)            idomfl(jfl) = 0
(  211)            IF( tmask(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) == 0. ) idomfl(jfl)=1
(  212)            
(  213)            ! Computation of the distance between the float
(  214)            ! and the faces of the mesh
(  215)            !            zdxab
(  216)            !             .
(  217)            !        B----.---------C
(  218)            !        |    .         |
(  219)            !        |<------>flo   |
(  220)            !        |        ^     |






PGF90 (Version     10.5)          09/23/2010  09:49:48      page 5

(  221)            !        |        |.....|....zdyad
(  222)            !        |        |     |
(  223)            !        A--------|-----D
(  224)            
(  225)            zdxab = dstnce(flxx(jfl),flyy(jfl),glamf(iimfl(jfl)-1,ijmfl(jfl)-1),flyy(jfl))                
(  226)            zdyad = dstnce(flxx(jfl),flyy(jfl),flxx(jfl),gphif(iimfl(jfl)-1,ijmfl(jfl)-1))
(  227)            
(  228)            ! Translation of this distances (in meter) in indexes
(  229)            
(  230)            !!DB: If I take tp?fl to denote T-cell, then T-cell(i,j) goes from
(  231)            !!     (i-1/2 ... i+1/2 , j-1/2 ... j+1/2) 
(  232)            
(  233)            tpifl(jfl) = (iimfl(jfl)-0.5)+zdxab/ e1u(iimfl(jfl)-1,ijmfl(jfl))+(mig(1)-jpizoom)
(  234)            tpjfl(jfl) = (ijmfl(jfl)-0.5)+zdyad/ e2v(iimfl(jfl),ijmfl(jfl)-1)+(mjg(1)-jpjzoom)
(  235)            tpkfl(jfl) = (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - flzz(jfl))*(ikmfl(jfl))                     &
(  236)                 / (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))   &
(  237)                 + (flzz(jfl) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))*(ikmfl(jfl)+1)                     &
(  238)                 / (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))
(  239)         ELSE
(  240)            tpifl (jfl) = 0.e0
(  241)            tpjfl (jfl) = 0.e0
(  242)            tpkfl (jfl) = 0.e0
(  243)            idomfl(jfl) = 0
(  244)         ENDIF
(  245)      END DO
(  246)      
(  247)      ! The sum of all the arrays tpifl, tpjfl, tpkfl give 3 arrays with the positions of all the floats. 
(  248)      IF( lk_mpp )   CALL mpp_sum( tpifl , jpnfl )   ! sums over the global domain
(  249)      IF( lk_mpp )   CALL mpp_sum( tpjfl , jpnfl )
(  250)      IF( lk_mpp )   CALL mpp_sum( tpkfl , jpnfl )
(  251)      IF( lk_mpp )   CALL mpp_sum( idomfl, jpnfl )
(  252)      
(  253)      !!DB
(  254)      ! Print the initial positions of the floats
(  255)      if( .NOT. ln_rstflo ) THEN 
(  256)         !!DB         ! WARNING : initial position not in the sea .OR. it's ijk position
(  257)         DO jfl = 1, jpnfl
(  258)            IF( idomfl(jfl) == 1 ) THEN
(  259)               IF(lwp) WRITE(numout,*)'*****************************'
(  260)               IF(lwp) WRITE(numout,*)'!!!!!!!  WARNING   !!!!!!!!!!'
(  261)               IF(lwp) WRITE(numout,*)'*****************************'
(  262)               IF(lwp) WRITE(numout,*)'The float number',jfl,'is out of the sea.'
(  263)               IF(lwp) WRITE(numout,*)'geographical position',flxx(jfl),flyy(jfl),flzz(jfl)
(  264)               IF(lwp) WRITE(numout,*)'index position',tpifl(jfl),tpjfl(jfl),tpkfl(jfl)
(  265)            else
(  266)               IF(lwp) WRITE(numout,*)'DBG: Float number ',jfl
(  267)               IF(lwp) WRITE(numout,*)'Geographical position',flxx(jfl),flyy(jfl),flzz(jfl)
(  268)               IF(lwp) WRITE(numout,*)'Index position',tpifl(jfl),tpjfl(jfl),tpkfl(jfl)
(  269)            endif
(  270)            
(  271)         END DO
(  272)      endif
(  273)      
(  274)    END SUBROUTINE flo_dom
(  275)          
(  276) 
(  277)    SUBROUTINE findmesh( pax, pay, pbx, pby,   &
(  278)                         pcx, pcy, pdx, pdy,   &






PGF90 (Version     10.5)          09/23/2010  09:49:48      page 6

(  279)                         px  ,py  ,ptx, pty, ldinmesh )
(  280)       !! -------------------------------------------------------------
(  281)       !!                ***  ROUTINE findmesh  ***
(  282)       !!      
(  283)       !! ** Purpose :   Find the index of mesh for the point spx spy
(  284)       !!
(  285)       !! ** Method  : 
(  286)       !!
(  287)       !! History :
(  288)       !!   8.0  !  98-07 (Y.Drillet)  Original code
(  289)       !!----------------------------------------------------------------------
(  290)       !! * Arguments
(  291)       REAL(wp) ::   &
(  292)          pax, pay, pbx, pby,    &     ! ???
(  293)          pcx, pcy, pdx, pdy,    &     ! ???
(  294)          px, py,                &     ! longitude and latitude
(  295)          ptx, pty                     ! ???
(  296)       LOGICAL ::  ldinmesh            ! ???
(  297) 
(  298)       !! * local declarations
(  299)       REAL(wp) ::   &
(  300)          zabt, zbct, zcdt, zdat, zabpt, zbcpt, zcdpt, zdapt,  &
(  301)          psax,psay,psbx,psby,psx,psy
(  302)       REAL(wp) ::  fsline                ! Statement function
(  303) 
(  304)       !! * Substitutions
(  305)       fsline(psax, psay, psbx, psby, psx, psy) = psy  * ( psbx - psax )   &
(  306)                                                - psx  * ( psby - psay )   &
(  307)                                                + psax *   psby - psay * psbx
(  308)       !!---------------------------------------------------------------------
(  309)       
(  310)       ! 4 semi plane defined by the 4 points and including the T point
(  311)       zabt = fsline(pax,pay,pbx,pby,ptx,pty)
(  312)       zbct = fsline(pbx,pby,pcx,pcy,ptx,pty)
(  313)       zcdt = fsline(pcx,pcy,pdx,pdy,ptx,pty)
(  314)       zdat = fsline(pdx,pdy,pax,pay,ptx,pty)
(  315)       
(  316)       ! 4 semi plane defined by the 4 points and including the extrememity
(  317)       zabpt = fsline(pax,pay,pbx,pby,px,py)
(  318)       zbcpt = fsline(pbx,pby,pcx,pcy,px,py)
(  319)       zcdpt = fsline(pcx,pcy,pdx,pdy,px,py)
(  320)       zdapt = fsline(pdx,pdy,pax,pay,px,py)
(  321)        
(  322)       ! We compare the semi plane T with the semi plane including the point
(  323)       ! to know if it is in this  mesh.
(  324)       ! For numerical reasons it is possible that for a point which is on
(  325)       ! the line we don't have exactly zero with fsline function. We want 
(  326)       ! that a point can't be in 2 mesh in the same time, so we put the 
(  327)       ! coefficient to zero if it is smaller than 1.E-12
(  328)       
(  329)       IF( ABS(zabpt) <= 1.E-12 ) zabpt = 0.
(  330)       IF( ABS(zbcpt) <= 1.E-12 ) zbcpt = 0.
(  331)       IF( ABS(zcdpt) <= 1.E-12 ) zcdpt = 0.
(  332)       IF( ABS(zdapt) <= 1.E-12 ) zdapt = 0.
(  333)       IF( (zabt*zabpt >  0.) .AND. (zbct*zbcpt >= 0. ) .AND. ( zcdt*zcdpt >= 0. ) .AND. ( zdat*zdapt > 0. )   &
(  334)          .AND. ( px <= MAX(pcx,pdx) ) .AND. ( px >= MIN(pax,pbx) )    &
(  335)          .AND. ( py <= MAX(pby,pcy) ) .AND. ( py >= MIN(pay,pdy) ) ) THEN
(  336)          ldinmesh=.TRUE.






PGF90 (Version     10.5)          09/23/2010  09:49:48      page 7

(  337)       ELSE
(  338)          ldinmesh=.FALSE.
(  339)       ENDIF
(  340) 
(  341)    END SUBROUTINE findmesh
(  342) 
(  343) 
(  344)    FUNCTION dstnce( pla1, phi1, pla2, phi2 )
(  345)       !! -------------------------------------------------------------
(  346)       !!                 ***  Function dstnce  ***
(  347)       !!          
(  348)       !! ** Purpose :   returns distance (in m) between two geographical
(  349)       !!                points
(  350)       !! ** Method  : 
(  351)       !!         
(  352)       !!----------------------------------------------------------------------
(  353)       !! * Arguments
(  354)       REAL(wp), INTENT(in) ::   pla1, phi1, pla2, phi2   ! ???
(  355) 
(  356)       !! * Local variables
(  357)       REAL(wp) ::   dly1, dly2, dlx1, dlx2, dlx, dls, dld, dpi
(  358)       REAL(wp) ::   dstnce
(  359)       !!---------------------------------------------------------------------
(  360)       
(  361)       dpi  = 2.* ASIN(1.)
(  362)       dls  = dpi / 180.
(  363)       dly1 = phi1 * dls
(  364)       dly2 = phi2 * dls
(  365)       dlx1 = pla1 * dls
(  366)       dlx2 = pla2 * dls
(  367) 
(  368)       dlx = SIN(dly1) * SIN(dly2) + COS(dly1) * COS(dly2) * COS(dlx2-dlx1)
(  369)  
(  370)       IF( ABS(dlx) > 1.0 ) dlx = 1.0
(  371) 
(  372)       dld = ATAN(DSQRT( ( 1-dlx )/( 1+dlx ) )) * 222.24 / dls
(  373)       dstnce = dld * 1000.
(  374) 
(  375)    END FUNCTION dstnce
(  376) 
(  377) #  else
(  378)    !!----------------------------------------------------------------------
(  379)    !!   Default option                                         Empty module
(  380)    !!----------------------------------------------------------------------
(  381) CONTAINS
(  382)    SUBROUTINE flo_dom                 ! Empty routine
(  383)    END SUBROUTINE flo_dom
(  384) #endif
(  385) 
(  386)    !!======================================================================
(  387) END MODULE flodom













PGF90 (Version     10.5)          09/23/2010  09:49:48      page 8

(    1) # 1 "flodom.F90"
(    1) !!DB 2008.05.16
(    2) !!Modified routine for more specific DB particle tracking needs
(    3) !!For older code, see OLD_CODE/ or look elsewhere for an older version
(    4) 
(    5) 
(    6) MODULE flodom
(    7)    !!======================================================================
(    8)    !!                       ***  MODULE  flodom  ***
(    9)    !! Ocean floats :   domain
(   10)    !!======================================================================
(   11) # 378
(  378)    !!----------------------------------------------------------------------
(  379)    !!   Default option                                         Empty module
(  380)    !!----------------------------------------------------------------------
(  381) CONTAINS
(  382)    SUBROUTINE flo_dom                 ! Empty routine
(  383)    END SUBROUTINE flo_dom
(  384) 
(  385) # 386
(  386)    !!======================================================================
(  387) END MODULE flodom
