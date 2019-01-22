


PGF90 (Version     10.5)          09/23/2010  09:49:35      page 1

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

Filename: domcfg.F90

(    1) MODULE domcfg
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE domcfg   ***
(    4)    !! Ocean initialization : domain configuration initialization
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_cfg        : initialize the domain configuration
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE phycst          ! physical constants
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE lib_mpp         ! distributed memory computing library
(   15)    USE solisl          ! ???
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! * Routine accessibility
(   21)    PUBLIC dom_cfg        ! called by opa.F90
(   22)    !!----------------------------------------------------------------------
(   23)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   24)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domcfg.F90,v 1.4 2006/04/10 15:46:07 opalod Exp $ 
(   25)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   26)    !!----------------------------------------------------------------------
(   27) 
(   28) CONTAINS
(   29) 
(   30)    SUBROUTINE dom_cfg
(   31)       !!----------------------------------------------------------------------
(   32)       !!                  ***  ROUTINE dom_cfg  ***
(   33)       !!                    
(   34)       !! ** Purpose :   set the domain configuration
(   35)       !!
(   36)       !! ** Method  :
(   37)       !!
(   38)       !! History :
(   39)       !!   9.0  !  03-09  (G. Madec)  Original code
(   40)       !!----------------------------------------------------------------------
(   41)       !! * Local declarations
(   42)       INTEGER ::   iconf = 0         ! temporary integers
(   43)       !!----------------------------------------------------------------------
(   44) 
(   45)       IF(lwp) THEN
(   46)          WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 2

(   47)          WRITE(numout,*) 'dom_cfg : set the ocean configuration'
(   48)          WRITE(numout,*) '~~~~~~~      ocean model configuration used :',   &
(   49)             &                             ' cp_cfg = ', cp_cfg, ' jp_cfg = ', jp_cfg
(   50)       ENDIF
(   51) 
(   52)       ! Global domain boundary conditions
(   53)       ! ---------------------------------
(   54)       IF(lwp) THEN
(   55)          WRITE(numout,*) '          global domain lateral boundaries'
(   56) 
(   57)          IF( jperio == 0 ) WRITE(numout,*) '             jperio= 0, closed'
(   58)          IF( jperio == 1 ) WRITE(numout,*) '             jperio= 1, cyclic east-west'
(   59)          IF( jperio == 2 ) WRITE(numout,*) '             jperio= 2, equatorial symmetric'
(   60)          IF( jperio == 3 ) WRITE(numout,*) '             jperio= 3, north fold with T-point pivot'
(   61)          IF( jperio == 4 ) WRITE(numout,*) '             jperio= 4, cyclic east-west and',   &
(   62)                                                                   ' north fold with T-point pivot'
(   63)          IF( jperio == 5 ) WRITE(numout,*) '             jperio= 5, north fold with F-point pivot'
(   64)          IF( jperio == 6 ) WRITE(numout,*) '             jperio= 6, cyclic east-west and',   &
(   65)                                                                   ' north fold with F-point pivot'
(   66)       ENDIF
(   67)       IF( jperio <  0 .OR. jperio > 6 ) THEN
(   68)           IF(lwp) WRITE(numout,cform_err)
(   69)           IF(lwp) WRITE(numout,*) 'jperio is out of range'
(   70)           nstop = nstop + 1
(   71)       ENDIF
(   72) 
(   73) 
(   74)       ! global domain versus zoom and/or local domain
(   75)       ! ---------------------------------------------
(   76) 
(   77)       CALL dom_glo 
(   78) 
(   79)    END SUBROUTINE dom_cfg
(   80) 
(   81) 
(   82)    SUBROUTINE dom_glo
(   83)       !!----------------------------------------------------------------------
(   84)       !!                     ***  ROUTINE dom_glo  ***
(   85)       !!
(   86)       !! ** Purpose :   initialization for global domain, zoom and local domain
(   87)       !!
(   88)       !! ** Method  :   
(   89)       !!
(   90)       !! ** Action  : - mig  , mjg : 
(   91)       !!              - mi0  , mi1   :
(   92)       !!              - mj0, , mj1   :
(   93)       !!
(   94)       !! History :
(   95)       !!   8.5  !  02-08  (G. Madec)    Original code
(   96)       !!----------------------------------------------------------------------
(   97)       !! * Local declarations
(   98)       INTEGER ::   ji, jj            ! dummy loop argument
(   99)       !!----------------------------------------------------------------------
(  100) 
(  101)       ! Local domain 
(  102)       ! ============
(  103) 
(  104)       ! local domain indices ==> data domain indices






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 3

(  105)       DO ji = 1, jpi
(  106)         mig(ji) = ji + jpizoom - 1 + nimpp - 1
(  107)       END DO
(  108)       DO jj = 1, jpj
(  109)         mjg(jj) = jj + jpjzoom - 1 + njmpp - 1
(  110)       END DO
(  111) 
(  112)       ! data domain indices ==> local domain indices
(  113)       ! (return (m.0,m.1)=(1,0) if data domain gridpoint is to the west/south of the 
(  114)       ! local domain, or (m.0,m.1)=(jp.+1,jp.) to the east/north of local domain. 
(  115)       DO ji = 1, jpidta
(  116)         mi0(ji) = MAX( 1, MIN( ji - jpizoom + 1 - nimpp + 1, jpi+1 ) )
(  117)         mi1(ji) = MAX( 0, MIN( ji - jpizoom + 1 - nimpp + 1, jpi   ) )
(  118)       END DO
(  119)       DO jj = 1, jpjdta
(  120)         mj0(jj) = MAX( 1, MIN( jj - jpjzoom + 1 - njmpp + 1, jpj+1 ) )
(  121)         mj1(jj) = MAX( 0, MIN( jj - jpjzoom + 1 - njmpp + 1, jpj   ) )
(  122)       END DO
(  123) 
(  124)       IF(lwp) THEN
(  125)          WRITE(numout,*)
(  126)          WRITE(numout,*) 'dom_glo : domain: data / local '
(  127)          WRITE(numout,*) '~~~~~~~ '
(  128)          WRITE(numout,*) '          data input domain    : jpidta = ', jpidta,   &
(  129)             &                                            ' jpjdta = ', jpjdta, ' jpkdta = ', jpkdta
(  130)          WRITE(numout,*) '          global or zoom domain: jpiglo = ', jpiglo,   &
(  131)             &                                            ' jpjglo = ', jpjglo, ' jpk    = ', jpk
(  132)          WRITE(numout,*) '          local domain         : jpi    = ', jpi   ,   &
(  133)             &                                            ' jpj    = ', jpj   , ' jpk    = ', jpk
(  134)          WRITE(numout,*)
(  135)          WRITE(numout,*) '          south-west indices    jpizoom = ', jpizoom,   &
(  136)             &                                           ' jpjzoom = ', jpjzoom
(  137)          WRITE(numout,*)
(  138)          WRITE(numout,*) '          conversion local  ==> data i-index domain'
(  139)          WRITE(numout,25)              (mig(ji),ji = 1,jpi)
(  140)          WRITE(numout,*)
(  141)          WRITE(numout,*) '          conversion data   ==> local  i-index domain'
(  142)          WRITE(numout,*) '             starting index'
(  143)          WRITE(numout,25)              (mi0(ji),ji = 1,jpidta)
(  144)          WRITE(numout,*) '             ending index'
(  145)          WRITE(numout,25)              (mi1(ji),ji = 1,jpidta)
(  146)          WRITE(numout,*)
(  147)          WRITE(numout,*) '          conversion local  ==> data j-index domain'
(  148)          WRITE(numout,25)              (mjg(jj),jj = 1,jpj)
(  149)          WRITE(numout,*)
(  150)          WRITE(numout,*) '          conversion data  ==> local j-index domain'
(  151)          WRITE(numout,*) '             starting index'
(  152)          WRITE(numout,25)              (mj0(jj),jj = 1,jpjdta)
(  153)          WRITE(numout,*) '             ending index'
(  154)          WRITE(numout,25)              (mj1(jj),jj = 1,jpjdta)
(  155)       ENDIF
(  156)  25   FORMAT( 100(10x,19i4,/) )
(  157) 
(  158)       ! Zoom domain
(  159)       ! ===========
(  160) 
(  161)       ! zoom control
(  162)       IF( jpiglo + jpizoom - 1  >  jpidta .OR.   &






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 4

(  163)           jpjglo + jpjzoom - 1  >  jpjdta      ) THEN
(  164)          IF(lwp)WRITE(numout,cform_err)
(  165)          IF(lwp)WRITE(numout,*)' global or zoom domain exceed the data domain ! '
(  166)          nstop = nstop + 1
(  167)       ENDIF
(  168) 
(  169)       ! set zoom flag
(  170)       IF ( jpiglo < jpidta .OR. jpjglo < jpjdta )   lzoom = .TRUE.
(  171) 
(  172)       ! set zoom type flags
(  173)       IF( lzoom .AND. jpizoom /= 1 )   lzoom_w = .TRUE.                     ! 
(  174)       IF( lzoom .AND. jpjzoom /= 1 )   lzoom_s = .TRUE.
(  175)       IF( lzoom .AND. jpiglo + jpizoom -1 /= jpidta )   lzoom_e = .TRUE.
(  176)       IF( lzoom .AND. jpjglo + jpjzoom -1 /= jpjdta )   lzoom_n = .TRUE.
(  177) 
(  178)       IF(lwp) THEN
(  179)          WRITE(numout,*)
(  180)          WRITE(numout,*) '          zoom flags : '
(  181)          WRITE(numout,*) '             lzoom   = ', lzoom  , ' (T = zoom, F = global )'
(  182)          WRITE(numout,*) '             lzoom_e = ', lzoom_e, ' (T = forced closed east  boundary)'
(  183)          WRITE(numout,*) '             lzoom_w = ', lzoom_w, ' (T = forced closed west  boundary)'
(  184)          WRITE(numout,*) '             lzoom_s = ', lzoom_s, ' (T = forced closed South boundary)'
(  185)          WRITE(numout,*) '             lzoom_n = ', lzoom_n, ' (T = forced closed North boundary)'
(  186)       ENDIF
(  187)       IF(  ( lzoom_e .OR. lzoom_w )  .AND.  ( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 )  ) THEN
(  188)          IF(lwp)WRITE(numout,cform_err)
(  189)          IF(lwp)WRITE(numout,*)' Your zoom choice is inconsistent with east-west cyclic boundary condition'
(  190)          nstop = nstop + 1
(  191)       ENDIF
(  192)       IF(  lzoom_n  .AND.  ( 3 <= jperio .AND. jperio <= 6 )  ) THEN
(  193)          IF(lwp)WRITE(numout,cform_err)
(  194)          IF(lwp)WRITE(numout,*)' Your zoom choice is inconsistent with North fold boundary condition'
(  195)          nstop = nstop + 1
(  196)       ENDIF
(  197)       IF(  lzoom  .AND.  lk_isl  ) THEN
(  198)          IF(lwp)WRITE(numout,cform_err)
(  199)          IF(lwp)WRITE(numout,*)' key_islands and zoom are not allowed'
(  200)          nstop = nstop + 1
(  201)       ENDIF
(  202) 
(  203) !!DB: ORCA-related
(  204)       ! Pre-defined arctic/antarctic zoom of ORCA configuration flag
(  205) !      IF( cp_cfg == "orca" ) THEN
(  206) !         SELECT CASE ( jp_cfg )
(  207) !         !                                        ! =======================
(  208) !         CASE ( 2 )                               !  ORCA_R2 configuration
(  209) !            !                                     ! =======================
(  210) !            IF(  jpiglo  == 142    .AND. jpjglo  ==  53 .AND.   &
(  211) !               & jpizoom ==  21    .AND. jpjzoom ==  97         )   lzoom_arct = .TRUE.
(  212) !            IF(  jpiglo  == jpidta .AND. jpjglo  ==  50 .AND.   &
(  213) !               & jpizoom ==   1    .AND. jpjzoom ==   1         )   lzoom_anta = .TRUE.
(  214) !            !                                     ! =======================
(  215) !         CASE ( 05 )                              !  ORCA_R05 configuration
(  216) !            !                                     ! =======================
(  217) !            IF(  jpiglo  == 562    .AND. jpjglo  == 202 .AND.   &
(  218) !               & jpizoom ==  81    .AND. jpjzoom == 301         )   lzoom_arct = .TRUE.
(  219) !            IF(  jpiglo  == jpidta .AND. jpjglo  == 187 .AND.   &
(  220) !               & jpizoom ==   1    .AND. jpjzoom ==   1         )   lzoom_anta = .TRUE.






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 5

(  221) !         END SELECT
(  222) !         !
(  223) !         IF(lwp) WRITE(numout,*) '          ORCA configuration: antarctic/arctic zoom flags : '
(  224) !         IF(lwp) WRITE(numout,*) '             lzoom_arct = ', lzoom_arct, ' (T=   arctic zoom, F=global)'
(  225) !         IF(lwp) WRITE(numout,*) '             lzoom_anta = ', lzoom_anta, ' (T=antarctic zoom, F=global)'
(  226) !         !
(  227) !      ENDIF
(  228)          
(  229)    END SUBROUTINE dom_glo
(  230) 
(  231)    !!======================================================================
(  232) END MODULE domcfg




















































PGF90 (Version     10.5)          09/23/2010  09:49:35      page 6

(    1) # 1 "domcfg.F90"
(    1) MODULE domcfg
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE domcfg   ***
(    4)    !! Ocean initialization : domain configuration initialization
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_cfg        : initialize the domain configuration
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE phycst          ! physical constants
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE lib_mpp         ! distributed memory computing library
(   15)    USE solisl          ! ???
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! * Routine accessibility
(   21)    PUBLIC dom_cfg        ! called by opa.F90
(   22)    !!----------------------------------------------------------------------
(   23)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   24)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domcfg.F90,v 1.4 2006/04/10 15:46:07 opalod Exp $ 
(   25)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   26)    !!----------------------------------------------------------------------
(   27) 
(   28) CONTAINS
(   29) 
(   30)    SUBROUTINE dom_cfg
(   31)       !!----------------------------------------------------------------------
(   32)       !!                  ***  ROUTINE dom_cfg  ***
(   33)       !!                    
(   34)       !! ** Purpose :   set the domain configuration
(   35)       !!
(   36)       !! ** Method  :
(   37)       !!
(   38)       !! History :
(   39)       !!   9.0  !  03-09  (G. Madec)  Original code
(   40)       !!----------------------------------------------------------------------
(   41)       !! * Local declarations
(   42)       INTEGER ::   iconf = 0         ! temporary integers
(   43)       !!----------------------------------------------------------------------
(   44) 
(   45)       IF(lwp) THEN
(   46)          WRITE(numout,*)
(   47)          WRITE(numout,*) 'dom_cfg : set the ocean configuration'
(   48)          WRITE(numout,*) '~~~~~~~      ocean model configuration used :',   &
(   49)             &                             ' cp_cfg = ', cp_cfg, ' jp_cfg = ', jp_cfg
(   50)       ENDIF
(   51) 
(   52)       ! Global domain boundary conditions
(   53)       ! ---------------------------------
(   54)       IF(lwp) THEN
(   55)          WRITE(numout,*) '          global domain lateral boundaries'
(   56) 
(   57)          IF( jperio == 0 ) WRITE(numout,*) '             jperio= 0, closed'






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 7

(   58)          IF( jperio == 1 ) WRITE(numout,*) '             jperio= 1, cyclic east-west'
(   59)          IF( jperio == 2 ) WRITE(numout,*) '             jperio= 2, equatorial symmetric'
(   60)          IF( jperio == 3 ) WRITE(numout,*) '             jperio= 3, north fold with T-point pivot'
(   61)          IF( jperio == 4 ) WRITE(numout,*) '             jperio= 4, cyclic east-west and',   &
(   62)                                                                   ' north fold with T-point pivot'
(   63)          IF( jperio == 5 ) WRITE(numout,*) '             jperio= 5, north fold with F-point pivot'
(   64)          IF( jperio == 6 ) WRITE(numout,*) '             jperio= 6, cyclic east-west and',   &
(   65)                                                                   ' north fold with F-point pivot'
(   66)       ENDIF
(   67)       IF( jperio <  0 .OR. jperio > 6 ) THEN
(   68)           IF(lwp) WRITE(numout,cform_err)
(   69)           IF(lwp) WRITE(numout,*) 'jperio is out of range'
(   70)           nstop = nstop + 1
(   71)       ENDIF
(   72) 
(   73) 
(   74)       ! global domain versus zoom and/or local domain
(   75)       ! ---------------------------------------------
(   76) 
(   77)       CALL dom_glo 
(   78) 
(   79)    END SUBROUTINE dom_cfg
(   80) 
(   81) 
(   82)    SUBROUTINE dom_glo
(   83)       !!----------------------------------------------------------------------
(   84)       !!                     ***  ROUTINE dom_glo  ***
(   85)       !!
(   86)       !! ** Purpose :   initialization for global domain, zoom and local domain
(   87)       !!
(   88)       !! ** Method  :   
(   89)       !!
(   90)       !! ** Action  : - mig  , mjg : 
(   91)       !!              - mi0  , mi1   :
(   92)       !!              - mj0, , mj1   :
(   93)       !!
(   94)       !! History :
(   95)       !!   8.5  !  02-08  (G. Madec)    Original code
(   96)       !!----------------------------------------------------------------------
(   97)       !! * Local declarations
(   98)       INTEGER ::   ji, jj            ! dummy loop argument
(   99)       !!----------------------------------------------------------------------
(  100) 
(  101)       ! Local domain 
(  102)       ! ============
(  103) 
(  104)       ! local domain indices ==> data domain indices
(  105)       DO ji = 1, jpi
(  106)         mig(ji) = ji + jpizoom - 1 + nimpp - 1
(  107)       END DO
(  108)       DO jj = 1, jpj
(  109)         mjg(jj) = jj + jpjzoom - 1 + njmpp - 1
(  110)       END DO
(  111) 
(  112)       ! data domain indices ==> local domain indices
(  113)       ! (return (m.0,m.1)=(1,0) if data domain gridpoint is to the west/south of the 
(  114)       ! local domain, or (m.0,m.1)=(jp.+1,jp.) to the east/north of local domain. 
(  115)       DO ji = 1, jpidta






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 8

(  116)         mi0(ji) = MAX( 1, MIN( ji - jpizoom + 1 - nimpp + 1, jpi+1 ) )
(  117)         mi1(ji) = MAX( 0, MIN( ji - jpizoom + 1 - nimpp + 1, jpi   ) )
(  118)       END DO
(  119)       DO jj = 1, jpjdta
(  120)         mj0(jj) = MAX( 1, MIN( jj - jpjzoom + 1 - njmpp + 1, jpj+1 ) )
(  121)         mj1(jj) = MAX( 0, MIN( jj - jpjzoom + 1 - njmpp + 1, jpj   ) )
(  122)       END DO
(  123) 
(  124)       IF(lwp) THEN
(  125)          WRITE(numout,*)
(  126)          WRITE(numout,*) 'dom_glo : domain: data / local '
(  127)          WRITE(numout,*) '~~~~~~~ '
(  128)          WRITE(numout,*) '          data input domain    : jpidta = ', jpidta,   &
(  129)             &                                            ' jpjdta = ', jpjdta, ' jpkdta = ', jpkdta
(  130)          WRITE(numout,*) '          global or zoom domain: jpiglo = ', jpiglo,   &
(  131)             &                                            ' jpjglo = ', jpjglo, ' jpk    = ', jpk
(  132)          WRITE(numout,*) '          local domain         : jpi    = ', jpi   ,   &
(  133)             &                                            ' jpj    = ', jpj   , ' jpk    = ', jpk
(  134)          WRITE(numout,*)
(  135)          WRITE(numout,*) '          south-west indices    jpizoom = ', jpizoom,   &
(  136)             &                                           ' jpjzoom = ', jpjzoom
(  137)          WRITE(numout,*)
(  138)          WRITE(numout,*) '          conversion local  ==> data i-index domain'
(  139)          WRITE(numout,25)              (mig(ji),ji = 1,jpi)
(  140)          WRITE(numout,*)
(  141)          WRITE(numout,*) '          conversion data   ==> local  i-index domain'
(  142)          WRITE(numout,*) '             starting index'
(  143)          WRITE(numout,25)              (mi0(ji),ji = 1,jpidta)
(  144)          WRITE(numout,*) '             ending index'
(  145)          WRITE(numout,25)              (mi1(ji),ji = 1,jpidta)
(  146)          WRITE(numout,*)
(  147)          WRITE(numout,*) '          conversion local  ==> data j-index domain'
(  148)          WRITE(numout,25)              (mjg(jj),jj = 1,jpj)
(  149)          WRITE(numout,*)
(  150)          WRITE(numout,*) '          conversion data  ==> local j-index domain'
(  151)          WRITE(numout,*) '             starting index'
(  152)          WRITE(numout,25)              (mj0(jj),jj = 1,jpjdta)
(  153)          WRITE(numout,*) '             ending index'
(  154)          WRITE(numout,25)              (mj1(jj),jj = 1,jpjdta)
(  155)       ENDIF
(  156)  25   FORMAT( 100(10x,19i4,/) )
(  157) 
(  158)       ! Zoom domain
(  159)       ! ===========
(  160) 
(  161)       ! zoom control
(  162)       IF( jpiglo + jpizoom - 1  >  jpidta .OR.   &
(  163)           jpjglo + jpjzoom - 1  >  jpjdta      ) THEN
(  164)          IF(lwp)WRITE(numout,cform_err)
(  165)          IF(lwp)WRITE(numout,*)' global or zoom domain exceed the data domain ! '
(  166)          nstop = nstop + 1
(  167)       ENDIF
(  168) 
(  169)       ! set zoom flag
(  170)       IF ( jpiglo < jpidta .OR. jpjglo < jpjdta )   lzoom = .TRUE.
(  171) 
(  172)       ! set zoom type flags
(  173)       IF( lzoom .AND. jpizoom /= 1 )   lzoom_w = .TRUE.                     ! 






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 9

(  174)       IF( lzoom .AND. jpjzoom /= 1 )   lzoom_s = .TRUE.
(  175)       IF( lzoom .AND. jpiglo + jpizoom -1 /= jpidta )   lzoom_e = .TRUE.
(  176)       IF( lzoom .AND. jpjglo + jpjzoom -1 /= jpjdta )   lzoom_n = .TRUE.
(  177) 
(  178)       IF(lwp) THEN
(  179)          WRITE(numout,*)
(  180)          WRITE(numout,*) '          zoom flags : '
(  181)          WRITE(numout,*) '             lzoom   = ', lzoom  , ' (T = zoom, F = global )'
(  182)          WRITE(numout,*) '             lzoom_e = ', lzoom_e, ' (T = forced closed east  boundary)'
(  183)          WRITE(numout,*) '             lzoom_w = ', lzoom_w, ' (T = forced closed west  boundary)'
(  184)          WRITE(numout,*) '             lzoom_s = ', lzoom_s, ' (T = forced closed South boundary)'
(  185)          WRITE(numout,*) '             lzoom_n = ', lzoom_n, ' (T = forced closed North boundary)'
(  186)       ENDIF
(  187)       IF(  ( lzoom_e .OR. lzoom_w )  .AND.  ( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 )  ) THEN
(  188)          IF(lwp)WRITE(numout,cform_err)
(  189)          IF(lwp)WRITE(numout,*)' Your zoom choice is inconsistent with east-west cyclic boundary condition'
(  190)          nstop = nstop + 1
(  191)       ENDIF
(  192)       IF(  lzoom_n  .AND.  ( 3 <= jperio .AND. jperio <= 6 )  ) THEN
(  193)          IF(lwp)WRITE(numout,cform_err)
(  194)          IF(lwp)WRITE(numout,*)' Your zoom choice is inconsistent with North fold boundary condition'
(  195)          nstop = nstop + 1
(  196)       ENDIF
(  197)       IF(  lzoom  .AND.  lk_isl  ) THEN
(  198)          IF(lwp)WRITE(numout,cform_err)
(  199)          IF(lwp)WRITE(numout,*)' key_islands and zoom are not allowed'
(  200)          nstop = nstop + 1
(  201)       ENDIF
(  202) 
(  203) !!DB: ORCA-related
(  204)       ! Pre-defined arctic/antarctic zoom of ORCA configuration flag
(  205) !      IF( cp_cfg == "orca" ) THEN
(  206) !         SELECT CASE ( jp_cfg )
(  207) !         !                                        ! =======================
(  208) !         CASE ( 2 )                               !  ORCA_R2 configuration
(  209) !            !                                     ! =======================
(  210) !            IF(  jpiglo  == 142    .AND. jpjglo  ==  53 .AND.   &
(  211) !               & jpizoom ==  21    .AND. jpjzoom ==  97         )   lzoom_arct = .TRUE.
(  212) !            IF(  jpiglo  == jpidta .AND. jpjglo  ==  50 .AND.   &
(  213) !               & jpizoom ==   1    .AND. jpjzoom ==   1         )   lzoom_anta = .TRUE.
(  214) !            !                                     ! =======================
(  215) !         CASE ( 05 )                              !  ORCA_R05 configuration
(  216) !            !                                     ! =======================
(  217) !            IF(  jpiglo  == 562    .AND. jpjglo  == 202 .AND.   &
(  218) !               & jpizoom ==  81    .AND. jpjzoom == 301         )   lzoom_arct = .TRUE.
(  219) !            IF(  jpiglo  == jpidta .AND. jpjglo  == 187 .AND.   &
(  220) !               & jpizoom ==   1    .AND. jpjzoom ==   1         )   lzoom_anta = .TRUE.
(  221) !         END SELECT
(  222) !         !
(  223) !         IF(lwp) WRITE(numout,*) '          ORCA configuration: antarctic/arctic zoom flags : '
(  224) !         IF(lwp) WRITE(numout,*) '             lzoom_arct = ', lzoom_arct, ' (T=   arctic zoom, F=global)'
(  225) !         IF(lwp) WRITE(numout,*) '             lzoom_anta = ', lzoom_anta, ' (T=antarctic zoom, F=global)'
(  226) !         !
(  227) !      ENDIF
(  228)          
(  229)    END SUBROUTINE dom_glo
(  230) 
(  231)    !!======================================================================






PGF90 (Version     10.5)          09/23/2010  09:49:35      page 10

(  232) END MODULE domcfg
