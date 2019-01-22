


PGF90 (Version     10.5)          09/23/2010  09:49:16      page 1

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

Filename: traadv_ctl.F90

(    1) MODULE traadv_ctl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_ctl  ***
(    4)    !! Ocean active tracers:  advection scheme control
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_ctl  : control the different options of advection scheme
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE in_out_manager  ! I/O manager
(   13) 
(   14)    IMPLICIT NONE
(   15)    PRIVATE
(   16) 
(   17)    !! * Accessibility
(   18)    PUBLIC tra_adv_ctl     ! routine called by step module
(   19)  
(   20)    !! * Share module variables
(   21)    LOGICAL, PUBLIC ::   &
(   22)       ln_traadv_cen2   = .TRUE.  ,   &  ! 2nd order centered scheme flag
(   23)       ln_traadv_tvd    = .FALSE. ,   &  ! TVD scheme flag
(   24)       ln_traadv_muscl  = .FALSE. ,   &  ! MUSCL scheme flag
(   25)       ln_traadv_muscl2 = .FALSE.        ! MUSCL2 scheme flag
(   26) 
(   27)    !!----------------------------------------------------------------------
(   28)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   29)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_ctl.F90,v 1.3 2005/12/21 10:46:45 opalod Exp $ 
(   30)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   31)    !!----------------------------------------------------------------------
(   32) 
(   33) CONTAINS
(   34) 
(   35)    SUBROUTINE tra_adv_ctl
(   36)       !!---------------------------------------------------------------------
(   37)       !!                  ***  ROUTINE tra_adv_ctl  ***
(   38)       !!                
(   39)       !! ** Purpose :   Control the consistency between cpp options for 
(   40)       !!      tracer advection schemes
(   41)       !!
(   42)       !! History :
(   43)       !!   8.5  !  02-11  (G. Madec)  Original code
(   44)       !!----------------------------------------------------------------------
(   45)       !! * Local declarations
(   46)       INTEGER ::   ioptio






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 2

(   47) 
(   48)       NAMELIST/nam_traadv/ ln_traadv_cen2 , ln_traadv_tvd,   &
(   49)          &                 ln_traadv_muscl, ln_traadv_muscl2
(   50)       !!----------------------------------------------------------------------
(   51) 
(   52)       ! Read Namelist nam_traadv : tracer advection scheme
(   53)       ! -------------------------
(   54)       REWIND ( numnam )
(   55)       READ   ( numnam, nam_traadv )
(   56) 
(   57)       ! Parameter control and print
(   58)       ! ---------------------------
(   59)       ! Control print
(   60)       IF(lwp) THEN
(   61)          WRITE(numout,*)
(   62)          WRITE(numout,*) 'tra_adv_ctl : choice/control of the tracer advection scheme'
(   63)          WRITE(numout,*) '~~~~~~~~~~~'
(   64)          WRITE(numout,*) '          Namelist nam_tra_adv : chose a advection scheme for tracers'
(   65)          WRITE(numout,*)
(   66)          WRITE(numout,*) '             2nd order advection scheme     ln_traadv_cen2   = ', ln_traadv_cen2
(   67)          WRITE(numout,*) '             TVD advection scheme           ln_traadv_tvd    = ', ln_traadv_tvd
(   68)          WRITE(numout,*) '             MUSCL  advection scheme        ln_traadv_muscl  = ', ln_traadv_muscl
(   69)          WRITE(numout,*) '             MUSCL2 advection scheme        ln_traadv_muscl2 = ', ln_traadv_muscl2
(   70)       ENDIF
(   71) 
(   72)       ! Control of Advection scheme options
(   73)       ! -----------------------------------
(   74)       ioptio = 0
(   75)       IF( ln_traadv_cen2   )   ioptio = ioptio + 1
(   76)       IF( ln_traadv_tvd    )   ioptio = ioptio + 1
(   77)       IF( ln_traadv_muscl  )   ioptio = ioptio + 1
(   78)       IF( ln_traadv_muscl2 )   ioptio = ioptio + 1
(   79) 
(   80)       IF( lk_esopa ) THEN
(   81)          IF(lwp) WRITE(numout,*) ' esopa control : the use of all scheme is forced'
(   82)          ln_traadv_cen2   = .TRUE.
(   83)          ln_traadv_tvd    = .TRUE.
(   84)          ln_traadv_muscl  = .TRUE.
(   85)          ln_traadv_muscl2 = .TRUE.
(   86)       ELSEIF( ioptio > 1 .OR. ioptio == 0 ) THEN
(   87)          IF(lwp) WRITE(numout,cform_err)
(   88)          IF(lwp) WRITE(numout,*) ' Choose one advection scheme in namelist nam_traadv'
(   89)          IF(lwp) WRITE(numout,*) '        ***                              ***********'
(   90)          nstop = nstop + 1
(   91)       ENDIF
(   92) 
(   93)       IF( n_cla == 1 .AND. .NOT. ln_traadv_cen2 ) THEN
(   94)          IF(lwp) WRITE(numout,cform_err)
(   95)          IF(lwp) WRITE(numout,*) '     cross-land advection only with 2nd order advection scheme'
(   96)          nstop = nstop + 1
(   97)       ENDIF
(   98) 
(   99)    END SUBROUTINE tra_adv_ctl
(  100) 
(  101)   !!======================================================================
(  102) END MODULE traadv_ctl








PGF90 (Version     10.5)          09/23/2010  09:49:16      page 3

(    1) # 1 "traadv_ctl.F90"
(    1) MODULE traadv_ctl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_ctl  ***
(    4)    !! Ocean active tracers:  advection scheme control
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_ctl  : control the different options of advection scheme
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE in_out_manager  ! I/O manager
(   13) 
(   14)    IMPLICIT NONE
(   15)    PRIVATE
(   16) 
(   17)    !! * Accessibility
(   18)    PUBLIC tra_adv_ctl     ! routine called by step module
(   19)  
(   20)    !! * Share module variables
(   21)    LOGICAL, PUBLIC ::   &
(   22)       ln_traadv_cen2   = .TRUE.  ,   &  ! 2nd order centered scheme flag
(   23)       ln_traadv_tvd    = .FALSE. ,   &  ! TVD scheme flag
(   24)       ln_traadv_muscl  = .FALSE. ,   &  ! MUSCL scheme flag
(   25)       ln_traadv_muscl2 = .FALSE.        ! MUSCL2 scheme flag
(   26) 
(   27)    !!----------------------------------------------------------------------
(   28)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   29)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_ctl.F90,v 1.3 2005/12/21 10:46:45 opalod Exp $ 
(   30)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   31)    !!----------------------------------------------------------------------
(   32) 
(   33) CONTAINS
(   34) 
(   35)    SUBROUTINE tra_adv_ctl
(   36)       !!---------------------------------------------------------------------
(   37)       !!                  ***  ROUTINE tra_adv_ctl  ***
(   38)       !!                
(   39)       !! ** Purpose :   Control the consistency between cpp options for 
(   40)       !!      tracer advection schemes
(   41)       !!
(   42)       !! History :
(   43)       !!   8.5  !  02-11  (G. Madec)  Original code
(   44)       !!----------------------------------------------------------------------
(   45)       !! * Local declarations
(   46)       INTEGER ::   ioptio
(   47) 
(   48)       NAMELIST/nam_traadv/ ln_traadv_cen2 , ln_traadv_tvd,   &
(   49)          &                 ln_traadv_muscl, ln_traadv_muscl2
(   50)       !!----------------------------------------------------------------------
(   51) 
(   52)       ! Read Namelist nam_traadv : tracer advection scheme
(   53)       ! -------------------------
(   54)       REWIND ( numnam )
(   55)       READ   ( numnam, nam_traadv )
(   56) 
(   57)       ! Parameter control and print






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 4

(   58)       ! ---------------------------
(   59)       ! Control print
(   60)       IF(lwp) THEN
(   61)          WRITE(numout,*)
(   62)          WRITE(numout,*) 'tra_adv_ctl : choice/control of the tracer advection scheme'
(   63)          WRITE(numout,*) '~~~~~~~~~~~'
(   64)          WRITE(numout,*) '          Namelist nam_tra_adv : chose a advection scheme for tracers'
(   65)          WRITE(numout,*)
(   66)          WRITE(numout,*) '             2nd order advection scheme     ln_traadv_cen2   = ', ln_traadv_cen2
(   67)          WRITE(numout,*) '             TVD advection scheme           ln_traadv_tvd    = ', ln_traadv_tvd
(   68)          WRITE(numout,*) '             MUSCL  advection scheme        ln_traadv_muscl  = ', ln_traadv_muscl
(   69)          WRITE(numout,*) '             MUSCL2 advection scheme        ln_traadv_muscl2 = ', ln_traadv_muscl2
(   70)       ENDIF
(   71) 
(   72)       ! Control of Advection scheme options
(   73)       ! -----------------------------------
(   74)       ioptio = 0
(   75)       IF( ln_traadv_cen2   )   ioptio = ioptio + 1
(   76)       IF( ln_traadv_tvd    )   ioptio = ioptio + 1
(   77)       IF( ln_traadv_muscl  )   ioptio = ioptio + 1
(   78)       IF( ln_traadv_muscl2 )   ioptio = ioptio + 1
(   79) 
(   80)       IF( lk_esopa ) THEN
(   81)          IF(lwp) WRITE(numout,*) ' esopa control : the use of all scheme is forced'
(   82)          ln_traadv_cen2   = .TRUE.
(   83)          ln_traadv_tvd    = .TRUE.
(   84)          ln_traadv_muscl  = .TRUE.
(   85)          ln_traadv_muscl2 = .TRUE.
(   86)       ELSEIF( ioptio > 1 .OR. ioptio == 0 ) THEN
(   87)          IF(lwp) WRITE(numout,cform_err)
(   88)          IF(lwp) WRITE(numout,*) ' Choose one advection scheme in namelist nam_traadv'
(   89)          IF(lwp) WRITE(numout,*) '        ***                              ***********'
(   90)          nstop = nstop + 1
(   91)       ENDIF
(   92) 
(   93)       IF( n_cla == 1 .AND. .NOT. ln_traadv_cen2 ) THEN
(   94)          IF(lwp) WRITE(numout,cform_err)
(   95)          IF(lwp) WRITE(numout,*) '     cross-land advection only with 2nd order advection scheme'
(   96)          nstop = nstop + 1
(   97)       ENDIF
(   98) 
(   99)    END SUBROUTINE tra_adv_ctl
(  100) 
(  101)   !!======================================================================
(  102) END MODULE traadv_ctl
