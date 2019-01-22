


PGF90 (Version     10.5)          09/23/2010  09:49:21      page 1

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

Filename: trctrp_lec.F90

(    1) MODULE trctrp_lec
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trctrp_lec  ***
(    4)    !! Ocean passive tracers:  namelist read options for transport
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_trp_lec  : read the passive tracer namelist for transport
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce_trc             ! ocean dynamics and active tracers variables
(   12)    USE trc                 ! ocean passive tracers variables
(   13) 
(   14)    IMPLICIT NONE
(   15)    PRIVATE
(   16) 
(   17)    !! * Accessibility
(   18)    PUBLIC trc_trp_lec     ! routine called by step module
(   19)  
(   20)    !! * Share module variables
(   21) 
(   22)    !! Advection
(   23)    LOGICAL, PUBLIC ::   &
(   24)       ln_trcadv_cen2   = .FALSE. ,  & !!: 2nd order centered scheme flag
(   25)       ln_trcadv_tvd    = .FALSE. ,  &  !: TVD scheme flag
(   26)       ln_trcadv_muscl  = .FALSE. ,  &  !: MUSCL scheme flag
(   27)       ln_trcadv_muscl2 = .FALSE. ,  &  !: MUSCL2 scheme flag
(   28)       ln_trcadv_smolar = .TRUE.        !: Smolarkiewicz scheme flag
(   29) 
(   30)    !! Lateral diffusion
(   31)    LOGICAL , PUBLIC ::              & !!: ** lateral mixing namelist (nam_trcldf) **
(   32)       ln_trcldf_diff  = .FALSE. ,   &  !: flag of perform or not the lateral diff.
(   33)       ln_trcldf_lap   = .TRUE.  ,   &  !: laplacian operator
(   34)       ln_trcldf_bilap = .FALSE. ,   &  !: bilaplacian operator
(   35)       ln_trcldf_level = .FALSE. ,   &  !: iso-level direction
(   36)       ln_trcldf_hor   = .FALSE. ,   &  !: horizontal (geopotential) direction
(   37)       ln_trcldf_iso   = .TRUE.         !: iso-neutral direction
(   38) 
(   39)    LOGICAL , PUBLIC ::              & !!: flag of the lateral diff. scheme used
(   40)       l_trcldf_lap         ,        &  !: iso-level laplacian operator
(   41)       l_trcldf_bilap       ,        &  !: iso-level bilaplacian operator
(   42)       l_trcldf_bilapg      ,        &  !: geopotential bilap. (s-coord)
(   43)       l_trcldf_iso         ,        &  !: iso-neutral laplacian or horizontal lapacian (s-coord)
(   44)       l_trczdf_iso         ,        &  !: idem for the vertical component
(   45)       l_trczdf_iso_vo      ,        &  !: idem with vectopt_memory
(   46)       l_trcldf_iso_zps                 !: iso-neutral laplacian (partial steps)






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 2

(   47) 
(   48)    !! Vertical diffusion
(   49)    LOGICAL , PUBLIC ::              & !!: nam_trczdf: vertical diffusion
(   50)       ln_trczdf_exp = .FALSE.          !: explicit vertical diffusion scheme flag
(   51) 
(   52)    INTEGER, PUBLIC ::               & !!: namzdf:  vertical diffusion
(   53)       n_trczdf_exp = 3                 !: number of sub-time step (explicit time stepping)
(   54) 
(   55)    LOGICAL, PUBLIC ::               &  !:
(   56)       l_trczdf_exp     = .FALSE. ,  &  !: explicit vertical diffusion
(   57)       l_trczdf_imp     = .FALSE.       !: implicit vertical diffusion
(   58) 
(   59) #if defined key_trcdmp
(   60)    !! Newtonian damping
(   61)    INTEGER  , PUBLIC ::             & !!: * newtonian damping namelist (nam_trcdmp) *
(   62)       ndmptr   =   -1 ,             &  !: = 0/-1/'latitude' for damping over tracers
(   63)       ndmpftr  =    2 ,             &  !: = 1 create a damping.coeff NetCDF file 
(   64)       nmldmptr =    0                  !: = 0/1/2 flag for damping in the mixed layer
(   65) 
(   66)    REAL(wp) , PUBLIC ::             & !!:  * newtonian damping namelist *
(   67)       sdmptr   =   50.,             &  !: surface time scale for internal damping (days)
(   68)       bdmptr   =  360.,             &  !: bottom time scale for internal damping (days)
(   69)       hdmptr   =  800.                 !: depth of transition between sdmp and bdmp (meters)
(   70) #endif
(   71)    !!----------------------------------------------------------------------
(   72)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   73)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trctrp_lec.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
(   74)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   75)    !!----------------------------------------------------------------------
(   76) 
(   77) CONTAINS
(   78) 
(   79)    SUBROUTINE trc_trp_lec
(   80)       !!---------------------------------------------------------------------
(   81)       !!                  ***  ROUTINE trc_trp_lec  ***
(   82)       !!                
(   83)       !! ** Purpose :   Read Namelist for tracer transport option
(   84)       !!
(   85)       !! History :
(   86)       !!   9.0  !  04-03  (C. Ethe) 
(   87)       !!----------------------------------------------------------------------
(   88)       !! * Local declarations
(   89) 
(   90)       NAMELIST/namtrcadv/ ln_trcadv_cen2 , ln_trcadv_tvd,   &
(   91)          &                 ln_trcadv_muscl, ln_trcadv_muscl2, ln_trcadv_smolar
(   92) 
(   93)       NAMELIST/namtrcldf/  ln_trcldf_diff  , ln_trcldf_lap  , ln_trcldf_bilap, &
(   94)          &                 ln_trcldf_level, ln_trcldf_hor, ln_trcldf_iso,   &
(   95)          &                 ahtrc0, ahtrb0, aeivtr0, trcrat
(   96) 
(   97)       NAMELIST/namtrczdf/ ln_trczdf_exp, n_trczdf_exp
(   98) 
(   99) #if defined key_trcdmp
(  100)       NAMELIST/namtrcdmp/ ndmptr, ndmpftr, nmldmptr, sdmptr, bdmptr, hdmptr
(  101) #endif
(  102)       !!----------------------------------------------------------------------
(  103) 
(  104)       ! Read Namelist namtrcadv : tracer advection scheme






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 3

(  105)       ! -------------------------
(  106)       REWIND ( numnat )
(  107)       READ   ( numnat, namtrcadv )
(  108) 
(  109)       ! Parameter control and print
(  110)       ! ---------------------------
(  111)       ! Control print
(  112)       IF(lwp) THEN
(  113)          WRITE(numout,*)
(  114)          WRITE(numout,*) 'choice/control of the tracer advection scheme'
(  115)          WRITE(numout,*) '~~~~~~~~~~~'
(  116)          WRITE(numout,*) '          Namelist namtrcadv : chose a advection scheme for tracers'
(  117)          WRITE(numout,*)
(  118)          WRITE(numout,*) '             2nd order advection scheme     ln_trcadv_cen2   = ', ln_trcadv_cen2
(  119)          WRITE(numout,*) '             TVD advection scheme           ln_trcadv_tvd    = ', ln_trcadv_tvd
(  120)          WRITE(numout,*) '             MUSCL  advection scheme        ln_trcadv_muscl  = ', ln_trcadv_muscl
(  121)          WRITE(numout,*) '             MUSCL2 advection scheme        ln_trcadv_muscl2 = ', ln_trcadv_muscl2
(  122)          WRITE(numout,*) '             SMOLARKIEWICZ advection scheme ln_trcadv_smolar = ', ln_trcadv_smolar
(  123)       ENDIF
(  124) 
(  125)       !  Define the lateral tracer physics parameters
(  126)       ! =============================================
(  127)     
(  128)       ! Read Namelist namtrcldf : Lateral physics on tracers
(  129)       REWIND( numnat )
(  130)       READ  ( numnat, namtrcldf )
(  131) 
(  132)       IF(lwp) THEN
(  133)          WRITE(numout,*)
(  134)          WRITE(numout,*) 'lateral passive tracer physics'
(  135)          WRITE(numout,*) '~~~~~~~'
(  136)          WRITE(numout,*) '   Namelist namtrcldf : set lateral mixing parameters (type, direction, coefficients)'
(  137)          WRITE(numout,*) '     perform lateral diffusion or not               ln_trcldf_diff  = ', ln_trcldf_diff
(  138)          WRITE(numout,*) '     laplacian operator                             ln_trcldf_lap   = ', ln_trcldf_lap
(  139)          WRITE(numout,*) '     bilaplacian operator                           ln_trcldf_bilap = ', ln_trcldf_bilap
(  140)          WRITE(numout,*) '     iso-level                                      ln_trcldf_level = ', ln_trcldf_level
(  141)          WRITE(numout,*) '     horizontal (geopotential)                      ln_trcldf_hor   = ', ln_trcldf_hor
(  142)          WRITE(numout,*) '     iso-neutral                                    ln_trcldf_iso   = ', ln_trcldf_iso
(  143)          WRITE(numout,*) '     lateral eddy diffusivity                              ahtrc0   = ', ahtrc0
(  144)          WRITE(numout,*) '     background hor. diffusivity                            ahtrb0  = ', ahtrb0
(  145)          WRITE(numout,*) '     eddy induced velocity coef.                           aeivtr0  = ', aeivtr0
(  146)          WRITE(numout,*) '     ratio between passive and active tracer diffusion coef  trcrat = ', trcrat
(  147)       ENDIF
(  148) 
(  149)       ! Read namtrczdf namelist : vertical mixing parameters
(  150)       ! --------------------
(  151)       REWIND( numnat )
(  152)       READ  ( numnat, namtrczdf )
(  153) 
(  154)       ! Parameter print
(  155)       ! ---------------
(  156)       IF(lwp) THEN
(  157)          WRITE(numout,*)
(  158)          WRITE(numout,*) 'vertical physics'
(  159)          WRITE(numout,*) '~~~~~~~~'
(  160)          WRITE(numout,*) '          Namelist namtrczdf : set vertical diffusion parameters'
(  161)          WRITE(numout,*) '             time splitting / backward scheme ln_trczdf_exp = ', ln_trczdf_exp
(  162)          WRITE(numout,*) '             number of time step               n_trczdf_exp = ', n_trczdf_exp






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 4

(  163)       ENDIF
(  164) 
(  165) # if defined key_trcdmp
(  166)       ! Read Namelist namtdp : passive tracres damping term
(  167)       ! --------------------
(  168)       REWIND ( numnat )
(  169)       READ   ( numnat, namtrcdmp )
(  170)       IF( lzoom ) THEN
(  171)          nmldmptr = 0           ! restoring to climatology at closed north or south boundaries
(  172)       ENDIF
(  173) 
(  174)       ! Parameter control and print
(  175)       ! ---------------------------
(  176)       IF(lwp) THEN
(  177)          WRITE(numout,*)
(  178)          WRITE(numout,*) 'newtonian damping'
(  179)          WRITE(numout,*) '~~~~~~~'
(  180)          WRITE(numout,*) '          Namelist namtrcdmp : set damping parameter'
(  181)          WRITE(numout,*)
(  182)          WRITE(numout,*) '             tracers damping option         ndmptr   = ', ndmptr
(  183)          WRITE(numout,*) '             create a damping.coeff file    ndmpftr  = ', ndmpftr
(  184)          WRITE(numout,*) '             mixed layer damping option     nmldmptr = ', nmldmptr, '(zoom: forced to 0)'
(  185)          WRITE(numout,*) '             surface time scale (days)      sdmptr   = ', sdmptr
(  186)          WRITE(numout,*) '             bottom time scale (days)       bdmptr   = ', bdmptr
(  187)          WRITE(numout,*) '             depth of transition (meters)   hdmptr   = ', hdmptr
(  188)          WRITE(numout,*)
(  189)       ENDIF
(  190) 
(  191) #endif
(  192) 
(  193)    END SUBROUTINE trc_trp_lec
(  194) #else
(  195)    !!----------------------------------------------------------------------
(  196)    !!   Dummy module :                      NO passive tracer
(  197)    !!----------------------------------------------------------------------
(  198) CONTAINS
(  199)    SUBROUTINE trc_trp_lec              ! Empty routine
(  200)    END SUBROUTINE trc_trp_lec
(  201) #endif
(  202)   !!======================================================================
(  203) END MODULE trctrp_lec























PGF90 (Version     10.5)          09/23/2010  09:49:21      page 5

(    1) # 1 "trctrp_lec.F90"
(    1) MODULE trctrp_lec
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trctrp_lec  ***
(    4)    !! Ocean passive tracers:  namelist read options for transport
(    5)    !!==============================================================================
(    6) # 195
(  195)    !!----------------------------------------------------------------------
(  196)    !!   Dummy module :                      NO passive tracer
(  197)    !!----------------------------------------------------------------------
(  198) CONTAINS
(  199)    SUBROUTINE trc_trp_lec              ! Empty routine
(  200)    END SUBROUTINE trc_trp_lec
(  201) # 202
(  202)   !!======================================================================
(  203) END MODULE trctrp_lec
