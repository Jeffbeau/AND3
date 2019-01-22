MODULE trabbl
   !!==============================================================================
   !!                       ***  MODULE  trabbl  ***
   !! Ocean physics :  advective and/or diffusive bottom boundary layer scheme
   !!==============================================================================
#if   defined key_trabbl_dif   ||   defined key_trabbl_adv   || defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_trabbl_dif'   or            diffusive bottom boundary layer
   !!----------------------------------------------------------------------
   !!   tra_bbl_dif  : update the active tracer trends due to the bottom
   !!                  boundary layer (diffusive only)
   !!   tra_bbl_init : initialization, namlist read, parameters control
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/TRA/trabbl.F90,v 1.2 2005/11/16 16:15:21 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_bbl_dif    ! routine called by step.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::            &  !!: * bbl namelist *
      atrbbl = 1.e+3                  !: lateral coeff. for bottom boundary 
      !                               !  layer scheme (m2/s)
   REAL(wp) , PUBLIC , DIMENSION(jpi,jpj) :: &
      bblx, bbly    ! Bottom boundary layer coefficients

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"

CONTAINS

   SUBROUTINE tra_bbl_dif( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_bbl_dif  ***
      !!
      !! ** Purpose :   Compute the before tracer (t & s) trend associated 
      !!      with the bottom boundary layer and add it to the general trend 
      !!      of tracer equations. The bottom boundary layer is supposed to be
      !!      a purely diffusive bottom boundary layer.
      !!
      !! ** Method  :   When the product grad( rho) * grad(h) < 0 (where grad 
      !!      is an along bottom slope gradient) an additional lateral diffu-
      !!      sive trend along the bottom slope is added to the general tracer
      !!      trend, otherwise nothing is done.
      !!      Second order operator (laplacian type) with variable coefficient
      !!      computed as follow for temperature (idem on s): 
      !!         difft = 1/(e1t*e2t*e3t) { di-1[ ahbt e2u*e3u/e1u di[ztb] ]
      !!                                 + dj-1[ ahbt e1v*e3v/e2v dj[ztb] ] }
      !!      where ztb is a 2D array: the bottom ocean temperature and ahtb
      !!      is a time and space varying diffusive coefficient defined by:
      !!         ahbt = zahbp    if grad(rho).grad(h) < 0
      !!              = 0.       otherwise.
      !!      Note that grad(.) is the along bottom slope gradient. grad(rho)
      !!      is evaluated using the local density (i.e. referenced at the
      !!      local depth). Typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!      Add this before trend to the general trend (ta,sa) of the 
      !!      botton ocean tracer point:
      !!         ta = ta + difft
      !!
      !! ** Action  : - update (ta,sa) at the bottom level with the bottom
      !!                boundary layer trend
      !!              - save the trends in bbltrd ('key_diatrends')
      !!
      !! References :
      !!     Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!
      !! History :
      !!   8.0  !  96-06  (L. Mortier)  Original code
      !!   8.0  !  97-11  (G. Madec)  Optimization
      !!   8.5  !  02-08  (G. Madec)  free form + modules
      !!----------------------------------------------------------------------
      !! * Arguments 
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !!----------------------------------------------------------------------

      IF( kt == nit000 )   CALL tra_bbl_init

   END SUBROUTINE tra_bbl_dif

   SUBROUTINE tra_bbl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_init  ***
      !!
      !! ** Purpose :   Initialization for the bottom boundary layer scheme.
      !!
      !! ** Method  :   Read the nambbl namelist and check the parameters
      !!      called by tra_bbl at the first timestep (nit000)
      !!
      !! History :
      !!    8.5  !  02-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      NAMELIST/nambbl/ atrbbl
      !!----------------------------------------------------------------------

      ! Read Namelist nambbl : bottom boundary layer scheme
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, nambbl )

      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*)
         WRITE(numout,*) 'tra_bbl_init : * Diffusive Bottom Boundary Layer'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist nambbl : set bbl parameters'
         WRITE(numout,*)
         WRITE(numout,*) '          bottom boundary layer coef.    atrbbl = ', atrbbl
         WRITE(numout,*)
      ENDIF
 
   END SUBROUTINE tra_bbl_init

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      No bottom boundary layer scheme
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE tra_bbl_dif (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'tra_bbl_dif: You should not have seen this print! error?', kt
   END SUBROUTINE tra_bbl_dif
#endif

   !!======================================================================
END MODULE trabbl
