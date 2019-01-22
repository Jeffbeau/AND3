MODULE trctrp_lec
   !!==============================================================================
   !!                       ***  MODULE  trctrp_lec  ***
   !! Ocean passive tracers:  namelist read options for transport
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_trp_lec  : read the passive tracer namelist for transport
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_trp_lec     ! routine called by step module
 
   !! * Share module variables

   !! Advection
   LOGICAL, PUBLIC ::   &
      ln_trcadv_cen2   = .FALSE. ,  & !!: 2nd order centered scheme flag
      ln_trcadv_tvd    = .FALSE. ,  &  !: TVD scheme flag
      ln_trcadv_muscl  = .FALSE. ,  &  !: MUSCL scheme flag
      ln_trcadv_muscl2 = .FALSE. ,  &  !: MUSCL2 scheme flag
      ln_trcadv_smolar = .TRUE.        !: Smolarkiewicz scheme flag

   !! Lateral diffusion
   LOGICAL , PUBLIC ::              & !!: ** lateral mixing namelist (nam_trcldf) **
      ln_trcldf_diff  = .FALSE. ,   &  !: flag of perform or not the lateral diff.
      ln_trcldf_lap   = .TRUE.  ,   &  !: laplacian operator
      ln_trcldf_bilap = .FALSE. ,   &  !: bilaplacian operator
      ln_trcldf_level = .FALSE. ,   &  !: iso-level direction
      ln_trcldf_hor   = .FALSE. ,   &  !: horizontal (geopotential) direction
      ln_trcldf_iso   = .TRUE.         !: iso-neutral direction

   LOGICAL , PUBLIC ::              & !!: flag of the lateral diff. scheme used
      l_trcldf_lap         ,        &  !: iso-level laplacian operator
      l_trcldf_bilap       ,        &  !: iso-level bilaplacian operator
      l_trcldf_bilapg      ,        &  !: geopotential bilap. (s-coord)
      l_trcldf_iso         ,        &  !: iso-neutral laplacian or horizontal lapacian (s-coord)
      l_trczdf_iso         ,        &  !: idem for the vertical component
      l_trczdf_iso_vo      ,        &  !: idem with vectopt_memory
      l_trcldf_iso_zps                 !: iso-neutral laplacian (partial steps)

   !! Vertical diffusion
   LOGICAL , PUBLIC ::              & !!: nam_trczdf: vertical diffusion
      ln_trczdf_exp = .FALSE.          !: explicit vertical diffusion scheme flag

   INTEGER, PUBLIC ::               & !!: namzdf:  vertical diffusion
      n_trczdf_exp = 3                 !: number of sub-time step (explicit time stepping)

   LOGICAL, PUBLIC ::               &  !:
      l_trczdf_exp     = .FALSE. ,  &  !: explicit vertical diffusion
      l_trczdf_imp     = .FALSE.       !: implicit vertical diffusion

#if defined key_trcdmp
   !! Newtonian damping
   INTEGER  , PUBLIC ::             & !!: * newtonian damping namelist (nam_trcdmp) *
      ndmptr   =   -1 ,             &  !: = 0/-1/'latitude' for damping over tracers
      ndmpftr  =    2 ,             &  !: = 1 create a damping.coeff NetCDF file 
      nmldmptr =    0                  !: = 0/1/2 flag for damping in the mixed layer

   REAL(wp) , PUBLIC ::             & !!:  * newtonian damping namelist *
      sdmptr   =   50.,             &  !: surface time scale for internal damping (days)
      bdmptr   =  360.,             &  !: bottom time scale for internal damping (days)
      hdmptr   =  800.                 !: depth of transition between sdmp and bdmp (meters)
#endif
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trctrp_lec.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_trp_lec
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trc_trp_lec  ***
      !!                
      !! ** Purpose :   Read Namelist for tracer transport option
      !!
      !! History :
      !!   9.0  !  04-03  (C. Ethe) 
      !!----------------------------------------------------------------------
      !! * Local declarations

      NAMELIST/namtrcadv/ ln_trcadv_cen2 , ln_trcadv_tvd,   &
         &                 ln_trcadv_muscl, ln_trcadv_muscl2, ln_trcadv_smolar

      NAMELIST/namtrcldf/  ln_trcldf_diff  , ln_trcldf_lap  , ln_trcldf_bilap, &
         &                 ln_trcldf_level, ln_trcldf_hor, ln_trcldf_iso,   &
         &                 ahtrc0, ahtrb0, aeivtr0, trcrat

      NAMELIST/namtrczdf/ ln_trczdf_exp, n_trczdf_exp

#if defined key_trcdmp
      NAMELIST/namtrcdmp/ ndmptr, ndmpftr, nmldmptr, sdmptr, bdmptr, hdmptr
#endif
      !!----------------------------------------------------------------------

      ! Read Namelist namtrcadv : tracer advection scheme
      ! -------------------------
      REWIND ( numnat )
      READ   ( numnat, namtrcadv )

      ! Parameter control and print
      ! ---------------------------
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist namtrcadv : chose a advection scheme for tracers'
         WRITE(numout,*)
         WRITE(numout,*) '             2nd order advection scheme     ln_trcadv_cen2   = ', ln_trcadv_cen2
         WRITE(numout,*) '             TVD advection scheme           ln_trcadv_tvd    = ', ln_trcadv_tvd
         WRITE(numout,*) '             MUSCL  advection scheme        ln_trcadv_muscl  = ', ln_trcadv_muscl
         WRITE(numout,*) '             MUSCL2 advection scheme        ln_trcadv_muscl2 = ', ln_trcadv_muscl2
         WRITE(numout,*) '             SMOLARKIEWICZ advection scheme ln_trcadv_smolar = ', ln_trcadv_smolar
      ENDIF

      !  Define the lateral tracer physics parameters
      ! =============================================
    
      ! Read Namelist namtrcldf : Lateral physics on tracers
      REWIND( numnat )
      READ  ( numnat, namtrcldf )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lateral passive tracer physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namtrcldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '     perform lateral diffusion or not               ln_trcldf_diff  = ', ln_trcldf_diff
         WRITE(numout,*) '     laplacian operator                             ln_trcldf_lap   = ', ln_trcldf_lap
         WRITE(numout,*) '     bilaplacian operator                           ln_trcldf_bilap = ', ln_trcldf_bilap
         WRITE(numout,*) '     iso-level                                      ln_trcldf_level = ', ln_trcldf_level
         WRITE(numout,*) '     horizontal (geopotential)                      ln_trcldf_hor   = ', ln_trcldf_hor
         WRITE(numout,*) '     iso-neutral                                    ln_trcldf_iso   = ', ln_trcldf_iso
         WRITE(numout,*) '     lateral eddy diffusivity                              ahtrc0   = ', ahtrc0
         WRITE(numout,*) '     background hor. diffusivity                            ahtrb0  = ', ahtrb0
         WRITE(numout,*) '     eddy induced velocity coef.                           aeivtr0  = ', aeivtr0
         WRITE(numout,*) '     ratio between passive and active tracer diffusion coef  trcrat = ', trcrat
      ENDIF

      ! Read namtrczdf namelist : vertical mixing parameters
      ! --------------------
      REWIND( numnat )
      READ  ( numnat, namtrczdf )

      ! Parameter print
      ! ---------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'vertical physics'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '          Namelist namtrczdf : set vertical diffusion parameters'
         WRITE(numout,*) '             time splitting / backward scheme ln_trczdf_exp = ', ln_trczdf_exp
         WRITE(numout,*) '             number of time step               n_trczdf_exp = ', n_trczdf_exp
      ENDIF

# if defined key_trcdmp
      ! Read Namelist namtdp : passive tracres damping term
      ! --------------------
      REWIND ( numnat )
      READ   ( numnat, namtrcdmp )
      IF( lzoom ) THEN
         nmldmptr = 0           ! restoring to climatology at closed north or south boundaries
      ENDIF

      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'newtonian damping'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '          Namelist namtrcdmp : set damping parameter'
         WRITE(numout,*)
         WRITE(numout,*) '             tracers damping option         ndmptr   = ', ndmptr
         WRITE(numout,*) '             create a damping.coeff file    ndmpftr  = ', ndmpftr
         WRITE(numout,*) '             mixed layer damping option     nmldmptr = ', nmldmptr, '(zoom: forced to 0)'
         WRITE(numout,*) '             surface time scale (days)      sdmptr   = ', sdmptr
         WRITE(numout,*) '             bottom time scale (days)       bdmptr   = ', bdmptr
         WRITE(numout,*) '             depth of transition (meters)   hdmptr   = ', hdmptr
         WRITE(numout,*)
      ENDIF

#endif

   END SUBROUTINE trc_trp_lec
#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_trp_lec              ! Empty routine
   END SUBROUTINE trc_trp_lec
#endif
  !!======================================================================
END MODULE trctrp_lec
