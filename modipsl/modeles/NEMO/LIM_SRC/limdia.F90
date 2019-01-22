MODULE limdia
   !!======================================================================
   !!                       ***  MODULE limdia   ***
   !!                      diagnostics of ice model 
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_dia      : computation of the time evolution of keys var.
   !!   lim_dia_init : initialization and namelist read
   !!----------------------------------------------------------------------
   !! * Modules used
   USE phycst          ! 
   USE par_ice         ! ice parameters
   USE ice_oce         ! ice variables
   USE daymod          !
   USE dom_ice         !
   USE ice             !
   USE iceini          !
   USE limistate       !
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC lim_dia       ! called by ice_step

   !! * Shared module variables
   INTEGER, PUBLIC  ::  &  !:
      ntmoy   = 1 ,     &  !: instantaneous values of ice evolution or averaging ntmoy
      ninfo   = 1          !: frequency of ouputs on file ice_evolu in case of averaging

   !! * Module variables
   INTEGER, PARAMETER ::   &  ! Parameters for outputs to files "evolu"
      jpinfmx = 100         ,    &  ! maximum number of key variables
      jpchinf = 5           ,    &  ! ???
      jpchsep = jpchinf + 2         ! ???

   INTEGER ::   &
      nfrinf  = 4 ,     &  ! number of variables written in one line 
      nferme ,          &  ! last time step at which the var. are written on file
      nvinfo ,          &  ! number of total variables 
      nbvt   ,          &  ! number of time variables
      naveg                ! number of step for accumulation before averaging

   CHARACTER(len=8) ::   &
      fmtinf  = '1PE13.5 ' ! format of the output values  
   CHARACTER(len=30) ::   &
      fmtw  ,           &  ! formats
      fmtr  ,           &  ! ???
      fmtitr               ! ???
   CHARACTER(len=jpchsep), DIMENSION(jpinfmx) ::   &
      titvar               ! title of key variables
 
   REAL(wp) ::   &
      epsi06 = 1.e-06      ! ???
   REAL(wp), DIMENSION(jpinfmx) ::  &
      vinfom               ! temporary working space
   REAL(wp), DIMENSION(jpi,jpj) ::   &
      aire                 ! masked grid cell area

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdia.F90,v 1.5 2005/03/27 18:34:41 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_dia
      !!--------------------------------------------------------------------
      !!                  ***  ROUTINE lim_dia  ***
      !!   
      !! ** Purpose : Computation and outputs on file ice.evolu 
      !!      the temporal evolution of some key variables
      !!
      !! History :
      !!   8.0  !  97-06  (Louvain-La-Neuve)  Original code
      !!   8.5  !  02-09  (C. Ethe , G. Madec )  F90: Free form and module
      !!-------------------------------------------------------------------
      !! * Local variables
       INTEGER  ::   jv,ji, jj   ! dummy loop indices
       INTEGER  ::   nv          ! indice of variable 
       REAL(wp), DIMENSION(jpinfmx) ::  & 
          vinfor           ! temporary working space 
       REAL(wp) ::    &
          zarea    ,    &  ! sea ice area
          zldarea  ,    &  ! leads area
          zextent15,    &  ! sea ice extent (15%)
          zextent85,    &  ! sea ice extent (85%)
          zicevol  ,    &  ! sea ice volume
          zsnwvol  ,    &  ! snow volume over sea ice
          zicespd          ! sea ice velocity
       !!-------------------------------------------------------------------

       IF( numit == nstart )   CALL lim_dia_init   ! initialisation of ice_evolu file      

       ! computation of key variables at each time step   

       nv = 1 
       vinfor(nv) = REAL( numit )
       nv = nv + 1
       vinfor(nv) = nyear
 
       DO jv = nbvt + 1, nvinfo
          vinfor(jv) = 0.e0
       END DO

       zextent15 = 0.e0
       zextent85 = 0.e0
       ! variables in northern Hemis
       DO jj = njeq, jpjm1
          DO ji = fs_2, fs_jpim1   ! vector opt.
             IF( tms(ji,jj) == 1 ) THEN
                zarea = ( 1.0 - frld(ji,jj) ) * aire(ji,jj)
                IF (frld(ji,jj) <= 0.15 ) zextent15 = aire(ji,jj)    
                IF (frld(ji,jj) <= 0.85 ) zextent85 = aire(ji,jj)   
                zldarea = zarea   / MAX( ( 1 - frld(ji,jj) ) , epsi06 )
                zicevol = zarea   * hicif(ji,jj)
                zsnwvol = zarea   * hsnif(ji,jj)
                zicespd = zicevol * ( u_ice(ji,jj) * u_ice(ji,jj)   &
                   &                + v_ice(ji,jj) * v_ice(ji,jj) )
                vinfor(nv+ 1) = vinfor(nv+ 1) + zarea
                vinfor(nv+ 3) = vinfor(nv+ 3) + zextent15
                vinfor(nv+ 5) = vinfor(nv+ 5) + zextent85
                vinfor(nv+ 7) = vinfor(nv+ 7) + zldarea
                vinfor(nv+ 9) = vinfor(nv+ 9) + zicevol
                vinfor(nv+11) = vinfor(nv+11) + zsnwvol
                vinfor(nv+13) = vinfor(nv+13) + zicespd
             ENDIF
          END DO
       END DO
       vinfor(nv+13) = SQRT( vinfor(nv+13) / MAX( vinfor(nv+9) , epsi06 ) )


      ! variables in southern Hemis
       nv = nv + 1
       DO jj = 2, njeqm1
          DO ji = fs_2, fs_jpim1   ! vector opt.
             IF( tms(ji,jj) == 1 ) THEN
                zarea = ( 1.0 - frld(ji,jj) ) * aire(ji,jj)
                IF (frld(ji,jj) <= 0.15 ) zextent15 = aire(ji,jj)    
                IF (frld(ji,jj) <= 0.85 ) zextent85 = aire(ji,jj)   
                zldarea = zarea   / MAX( ( 1 - frld(ji,jj) ) , epsi06 )
                zicevol = zarea   * hicif(ji,jj)
                zsnwvol = zarea   * hsnif(ji,jj)
                zicespd = zicevol * ( u_ice(ji,jj) * u_ice(ji,jj)   &
                   &                + v_ice(ji,jj) * v_ice(ji,jj) )
                vinfor(nv+ 1) = vinfor(nv+ 1) + zarea
                vinfor(nv+ 3) = vinfor(nv+ 3) + zextent15
                vinfor(nv+ 5) = vinfor(nv+ 5) + zextent85
                vinfor(nv+ 7) = vinfor(nv+ 7) + zldarea
                vinfor(nv+ 9) = vinfor(nv+ 9) + zicevol
                vinfor(nv+11) = vinfor(nv+11) + zsnwvol
                vinfor(nv+13) = vinfor(nv+13) + zicespd
             ENDIF
          END DO
       END DO
       vinfor(nv+13) = SQRT( vinfor(nv+13) / MAX( vinfor(nv+9) , epsi06 ) )    

       !  Accumulation before averaging 
       DO jv = 1, nvinfo
          vinfom(jv) = vinfom(jv) + vinfor(jv)
       END DO
       naveg = naveg + 1  
    
       ! oututs on file ice_evolu    
       IF( MOD( numit , ninfo ) == 0 ) THEN
          WRITE(numevo_ice,fmtw) ( titvar(jv), vinfom(jv)/naveg, jv = 1, nvinfo )
          naveg = 0
          DO jv = 1, nvinfo
             vinfom(jv) = 0.e0
          END DO
       ENDIF
  
    END SUBROUTINE lim_dia
 

    SUBROUTINE lim_dia_init
       !!-------------------------------------------------------------------
       !!                  ***  ROUTINE lim_dia_init  ***
       !!             
       !! ** Purpose : Preparation of the file ice_evolu for the output of
       !!      the temporal evolution of key variables
       !!
       !! ** input   : Namelist namicedia
       !!
       !! history :
       !!  8.5  ! 03-08 (C. Ethe) original code
       !!-------------------------------------------------------------------
       NAMELIST/namicedia/fmtinf, nfrinf, ninfo, ntmoy

       INTEGER  ::   jv   ,     &  ! dummy loop indice
          &          ntot ,     &
          &          ndeb ,     &
          &          irecl

       INTEGER  ::   nv            ! indice of variable 

       REAL(wp) ::   zxx0, zxx1    ! temporary scalars

       CHARACTER(len=jpchinf) ::   titinf
       !!-------------------------------------------------------------------

       ! Read Namelist namicedia
       REWIND ( numnam_ice )
       READ   ( numnam_ice  , namicedia )
       IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'lim_dia_init : ice parameters for ice diagnostics '
          WRITE(numout,*) '~~~~~~~~~~~~'
          WRITE(numout,*) '   format of the output values                                 fmtinf = ', fmtinf
          WRITE(numout,*) '   number of variables written in one line                     nfrinf = ', nfrinf 
          WRITE(numout,*) '   Instantaneous values of ice evolution or averaging          ntmoy  = ', ntmoy
          WRITE(numout,*) '   frequency of ouputs on file ice_evolu in case of averaging  ninfo  = ', ninfo
       ENDIF

       ! masked grid cell area
       aire(:,:) = area(:,:) * tms(:,:)

       ! Titles of ice key variables :
       nv = 1
       titvar(nv) = 'NoIt'  ! iteration number
       nv = nv + 1
       titvar(nv) = 'T yr'  ! time step in years
       nv = nv + 1

       nbvt = nv - 1

       titvar(nv) = 'AEFN' ! sea ice area in the northern Hemisp.(10^12 km2)
       nv = nv + 1
       titvar(nv) = 'AEFS' ! sea ice area in the southern Hemisp.(10^12 km2)
       nv = nv + 1
       titvar(nv) = 'A15N'  ! sea ice extent (15%) in the northern Hemisp.(10^12 km2)
       nv = nv + 1
       titvar(nv) = 'A15S'  ! sea ice extent (15%) in the southern Hemisp.(10^12 km2)
       nv = nv + 1
       titvar(nv) = 'A85N'  ! sea ice extent (85%) in the northern Hemisp.(10^12 km2)
       nv = nv + 1
       titvar(nv) = 'A85S'  ! sea ice extent (85%) in the southern Hemisp.(10^12 km2)
       nv = nv + 1
       titvar(nv) = 'ALEN'  ! leads area in the northern Hemisp.(10^12 km2)
       nv = nv + 1
       titvar(nv) = 'ALES'  ! leads area in the southern Hemisp.(10^12 km2)
       nv = nv + 1
       titvar(nv) = 'VOLN'  ! sea ice volume in the northern Hemisp.(10^3 km3)
       nv = nv + 1
       titvar(nv) = 'VOLS'  ! sea ice volume in the southern Hemisp.(10^3 km3)
       nv = nv + 1
       titvar(nv) = 'VONN'  ! snow volume over sea ice in the northern Hemisp.(10^3 km3)
       nv = nv + 1
       titvar(nv) = 'VONS'  ! snow volume over sea ice in the southern Hemisp.(10^3 km3)
       nv = nv + 1
       titvar(nv) = 'ECGN'  ! mean sea ice velocity in the northern Hemisp.(m/s)
       nv = nv + 1
       titvar(nv) = 'ECGS'  ! mean sea ice velocity in the southern Hemisp.(m/s)

       nvinfo = nv

       ! Definition et Ecriture de l'entete : nombre d'enregistrements 
       ndeb   = ( nstart - 1 ) / ninfo
       IF( nstart == 1 ) ndeb = -1

       nferme = ( nstart - 1 + nitrun) / ninfo
       ntot   = nferme - ndeb
       ndeb   = ninfo * ( 1 + ndeb )
       nferme = ninfo * nferme

       ! definition of formats 
       WRITE( fmtw  , '(A,I3,A2,I1,A)' )  '(', nfrinf, '(A', jpchsep, ','//fmtinf//'))'
       WRITE( fmtr  , '(A,I3,A,I1,A)'  )  '(', nfrinf, '(', jpchsep, 'X,'//fmtinf//'))'
       WRITE( fmtitr, '(A,I3,A,I1,A)'  )  '(', nvinfo, 'A', jpchinf, ')'

       ! opening  "ice_evolu" file
       irecl = ( jpchinf + 1 ) * nvinfo 
       OPEN( numevo_ice, file='ice.evolu', status='unknown', RECL = irecl)
       OPEN( numevo_ice, file='ice.evolu', status='unknown')

       !- ecriture de 2 lignes d''entete :
       WRITE(numevo_ice,1000) fmtr, fmtw, fmtitr, nvinfo, ntot, 0, nfrinf
       zxx0 = 0.001 * REAL( ninfo )
       zxx1 = 0.001 * REAL( ndeb  )
       WRITE(numevo_ice,1111) REAL(jpchinf), 0., zxx1, zxx0, 0., 0., 0

       !- ecriture de 2 lignes de titre :
       WRITE(numevo_ice,'(A,I8,A,I8,A,I5)')                                      &
          'Evolution chronologique - Experience '//cexper   &
          //'   de', ndeb, ' a', nferme, ' pas', ninfo
       WRITE(numevo_ice,fmtitr) ( titvar(jv), jv = 1, nvinfo )


       !--preparation de "titvar" pour l''ecriture parmi les valeurs numeriques :
       DO  jv = 2 , nvinfo
          titinf     = titvar(jv)(:jpchinf)
          titvar(jv) = '  '//titinf
       END DO

       !--Initialisation of the arrays for the accumulation
       DO  jv = 1, nvinfo
          vinfom(jv) = 0.
       END DO
       naveg = 0

1000   FORMAT( 3(A20),4(1x,I6) )
1111   FORMAT( 3(F7.1,1X,F7.3,1X),I3,A )  

    END SUBROUTINE lim_dia_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :                               NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_dia         ! Empty routine
   END SUBROUTINE lim_dia
#endif

   !!======================================================================
END MODULE limdia
