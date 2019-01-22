


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

Filename: domstp.F90

(    1) MODULE domstp
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE domstp   ***
(    4)    !! Ocean initialization : time domain
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_stp        : ocean time domain initialization
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE in_out_manager  ! I/O manager
(   14) 
(   15)    IMPLICIT NONE
(   16)    PRIVATE
(   17) 
(   18)    !! * routine accessibility
(   19)    PUBLIC dom_stp        ! routine called by inidom.F90
(   20) 
(   21)    !! * Substitutions
(   22) #  include "domzgr_substitute.h90"
(   23)    !!----------------------------------------------------------------------
(   24)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   25)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domstp.F90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   26)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   27)    !!----------------------------------------------------------------------
(   28) 
(   29) CONTAINS
(   30) 
(   31)    SUBROUTINE dom_stp
(   32)       !!----------------------------------------------------------------------
(   33)       !!                    ***  ROUTINE dom_stp  ***
(   34)       !!          
(   35)       !! ** Purpose :   Intialize ocean time step for the run
(   36)       !!
(   37)       !! ** Method  : - Initialization of a coef. use in the Asselin time
(   38)       !!      filter:  atfp1 = 1 - 2 * atfp  where atfp is the Asselin time
(   39)       !!      filter parameter read in namelist
(   40)       !!              - Model time step:
(   41)       !!      nacc = 0 : synchronous time intergration. 
(   42)       !!      There is one time step only, defined by: rdt, rdttra(k)=rdt
(   43)       !!      nacc = 1 : accelerating the convergence. There is 2 different
(   44)       !!      time steps for dynamics and tracers:
(   45)       !!        rdt      : dynamical part
(   46)       !!        rdttra(k): temperature and salinity






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 2

(   47)       !!      The tracer time step is a function of vertical level. the model
(   48)       !!      reference time step ( i.e. for wind stress, surface heat and
(   49)       !!      salt fluxes) is the surface tracer time step is rdttra(1).
(   50)       !!         N.B. depth dependent acceleration of convergence is not im-
(   51)       !!      plemented for s-coordinate.
(   52)       !!
(   53)       !! ** Action  : - rdttra   : vertical profile of tracer time step
(   54)       !!              - atfp1    : = 1 - 2*atfp
(   55)       !!
(   56)       !! References :
(   57)       !!      Bryan, K., 1984, J. Phys. Oceanogr., 14, 666-673.
(   58)       !!
(   59)       !! History :
(   60)       !!        !  90-10  (O. Marti)  Original code
(   61)       !!        !  96-01  (G. Madec)  terrain following coordinates
(   62)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   63)       !!----------------------------------------------------------------------
(   64)       !! * Local declarations
(   65)       INTEGER ::   jk              ! dummy loop indice
(   66)       !!----------------------------------------------------------------------
(   67) 
(   68)       IF(lwp) THEN
(   69)          WRITE(numout,*)
(   70)          WRITE(numout,*) 'dom_stp : time stepping setting'
(   71)          WRITE(numout,*) '~~~~~~~'
(   72)       ENDIF
(   73) 
(   74)       ! 0. Asselin Time filter
(   75)       ! ----------------------
(   76)       
(   77)       atfp1 = 1. - 2. * atfp
(   78) 
(   79) 
(   80)       SELECT CASE ( nacc )
(   81) 
(   82)          CASE ( 0 )                ! Synchronous time stepping
(   83)             IF(lwp) WRITE(numout,*)'               synchronous time stepping'
(   84)             IF(lwp) WRITE(numout,*)'               dynamics and tracer time step = ', rdt/3600., ' hours'
(   85) 
(   86)             rdttra(:) = rdt
(   87) 
(   88)          CASE ( 1 )                ! Accelerating the convergence
(   89)             IF(lwp) WRITE(numout,*) '              no tracer damping in the turbocline'
(   90)             IF(lwp) WRITE(numout,*)'               accelerating the convergence'
(   91)             IF(lwp) WRITE(numout,*)'               dynamics time step = ', rdt/3600., ' hours'
(   92) #if defined key_s_coord
(   93)             IF( rdtmin /= rdtmax ) THEN
(   94)                IF(lwp) WRITE(numout,cform_err)
(   95)                IF(lwp) WRITE(numout,*)' depth dependent acceleration of &
(   96)                                       &convergence not implemented in s-coordinates'
(   97)                nstop = nstop + 1
(   98)             ENDIF
(   99) #endif
(  100) #if defined key_partial_steps
(  101)             IF( rdtmin /= rdtmax ) THEN
(  102)                IF(lwp) WRITE(numout,cform_err)
(  103)                IF(lwp) WRITE(numout,*)' depth dependent acceleration of &
(  104)                                       &convergence not implemented for partial steps case'






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 3

(  105)                nstop = nstop + 1
(  106)             ENDIF
(  107) #endif
(  108)             IF(lwp) WRITE(numout,*)'         tracers   time step :  dt (hours)  level'
(  109) 
(  110)             DO jk = 1, jpk
(  111)                IF( fsdept(1,1,jk) <= rdth ) rdttra(jk) = rdtmin
(  112)                IF( fsdept(1,1,jk) >  rdth ) THEN
(  113)                   rdttra(jk) = rdtmin + ( rdtmax - rdtmin )   &
(  114)                                       * ( EXP( ( fsdept(1,1,jk ) - rdth ) / rdth ) - 1. )   &
(  115)                                       / ( EXP( ( fsdept(1,1,jpk) - rdth ) / rdth ) - 1. )
(  116)                ENDIF
(  117)                IF(lwp) WRITE(numout,9200) rdttra(jk)/3600., jk
(  118)             END DO  
(  119)  9200       FORMAT(36x,f5.2,'     ',i3)
(  120) 
(  121)          CASE DEFAULT              ! E R R O R 
(  122) 
(  123)             IF(lwp) WRITE(numout,cform_err)
(  124)             IF(lwp) WRITE(numout,*) ' nacc value e r r o r, nacc= ',nacc
(  125)             IF(lwp) WRITE(numout,*) ' we stop'
(  126)             nstop = nstop + 1
(  127) 
(  128)       END SELECT
(  129) 
(  130)    END SUBROUTINE dom_stp
(  131) 
(  132)    !!======================================================================
(  133) END MODULE domstp



































PGF90 (Version     10.5)          09/23/2010  09:49:33      page 4

(    1) # 1 "domstp.F90"
(    1) MODULE domstp
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE domstp   ***
(    4)    !! Ocean initialization : time domain
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_stp        : ocean time domain initialization
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE in_out_manager  ! I/O manager
(   14) 
(   15)    IMPLICIT NONE
(   16)    PRIVATE
(   17) 
(   18)    !! * routine accessibility
(   19)    PUBLIC dom_stp        ! routine called by inidom.F90
(   20) 
(   21)    !! * Substitutions
(   22) # 1 "./domzgr_substitute.h90"
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
(   48) # 23 "domstp.F90"
(   23) # 23
(   23)    !!----------------------------------------------------------------------
(   24)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   25)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domstp.F90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   26)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   27)    !!----------------------------------------------------------------------
(   28) 
(   29) CONTAINS
(   30) 
(   31)    SUBROUTINE dom_stp
(   32)       !!----------------------------------------------------------------------
(   33)       !!                    ***  ROUTINE dom_stp  ***
(   34)       !!          
(   35)       !! ** Purpose :   Intialize ocean time step for the run
(   36)       !!
(   37)       !! ** Method  : - Initialization of a coef. use in the Asselin time
(   38)       !!      filter:  atfp1 = 1 - 2 * atfp  where atfp is the Asselin time
(   39)       !!      filter parameter read in namelist
(   40)       !!              - Model time step:
(   41)       !!      nacc = 0 : synchronous time intergration. 






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 5

(   42)       !!      There is one time step only, defined by: rdt, rdttra(k)=rdt
(   43)       !!      nacc = 1 : accelerating the convergence. There is 2 different
(   44)       !!      time steps for dynamics and tracers:
(   45)       !!        rdt      : dynamical part
(   46)       !!        rdttra(k): temperature and salinity
(   47)       !!      The tracer time step is a function of vertical level. the model
(   48)       !!      reference time step ( i.e. for wind stress, surface heat and
(   49)       !!      salt fluxes) is the surface tracer time step is rdttra(1).
(   50)       !!         N.B. depth dependent acceleration of convergence is not im-
(   51)       !!      plemented for s-coordinate.
(   52)       !!
(   53)       !! ** Action  : - rdttra   : vertical profile of tracer time step
(   54)       !!              - atfp1    : = 1 - 2*atfp
(   55)       !!
(   56)       !! References :
(   57)       !!      Bryan, K., 1984, J. Phys. Oceanogr., 14, 666-673.
(   58)       !!
(   59)       !! History :
(   60)       !!        !  90-10  (O. Marti)  Original code
(   61)       !!        !  96-01  (G. Madec)  terrain following coordinates
(   62)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   63)       !!----------------------------------------------------------------------
(   64)       !! * Local declarations
(   65)       INTEGER ::   jk              ! dummy loop indice
(   66)       !!----------------------------------------------------------------------
(   67) 
(   68)       IF(lwp) THEN
(   69)          WRITE(numout,*)
(   70)          WRITE(numout,*) 'dom_stp : time stepping setting'
(   71)          WRITE(numout,*) '~~~~~~~'
(   72)       ENDIF
(   73) 
(   74)       ! 0. Asselin Time filter
(   75)       ! ----------------------
(   76)       
(   77)       atfp1 = 1. - 2. * atfp
(   78) 
(   79) 
(   80)       SELECT CASE ( nacc )
(   81) 
(   82)          CASE ( 0 )                ! Synchronous time stepping
(   83)             IF(lwp) WRITE(numout,*)'               synchronous time stepping'
(   84)             IF(lwp) WRITE(numout,*)'               dynamics and tracer time step = ', rdt/3600., ' hours'
(   85) 
(   86)             rdttra(:) = rdt
(   87) 
(   88)          CASE ( 1 )                ! Accelerating the convergence
(   89)             IF(lwp) WRITE(numout,*) '              no tracer damping in the turbocline'
(   90)             IF(lwp) WRITE(numout,*)'               accelerating the convergence'
(   91)             IF(lwp) WRITE(numout,*)'               dynamics time step = ', rdt/3600., ' hours'
(   92) # 108
(  108)             IF(lwp) WRITE(numout,*)'         tracers   time step :  dt (hours)  level'
(  109) 
(  110)             DO jk = 1, jpk
(  111)                IF( gdept(jk) <= rdth ) rdttra(jk) = rdtmin
(  112)                IF( gdept(jk) >  rdth ) THEN
(  113)                   rdttra(jk) = rdtmin + ( rdtmax - rdtmin )   &
(  114)                                       * ( EXP( ( gdept(jk ) - rdth ) / rdth ) - 1. )   &






PGF90 (Version     10.5)          09/23/2010  09:49:33      page 6

(  115)                                       / ( EXP( ( gdept(jpk) - rdth ) / rdth ) - 1. )
(  116)                ENDIF
(  117)                IF(lwp) WRITE(numout,9200) rdttra(jk)/3600., jk
(  118)             END DO  
(  119)  9200       FORMAT(36x,f5.2,'     ',i3)
(  120) 
(  121)          CASE DEFAULT              ! E R R O R 
(  122) 
(  123)             IF(lwp) WRITE(numout,cform_err)
(  124)             IF(lwp) WRITE(numout,*) ' nacc value e r r o r, nacc= ',nacc
(  125)             IF(lwp) WRITE(numout,*) ' we stop'
(  126)             nstop = nstop + 1
(  127) 
(  128)       END SELECT
(  129) 
(  130)    END SUBROUTINE dom_stp
(  131) 
(  132)    !!======================================================================
(  133) END MODULE domstp
