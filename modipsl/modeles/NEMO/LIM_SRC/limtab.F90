MODULE limtab
   !!======================================================================
   !!                       ***  MODULE limtab   ***
   !!              transform 1D (2D) array to a 2D (1D) table
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   tab_2d_1d  : 2-D to 1-D
   !!   tab_1d_2d  : 1-D to 2-D
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tab_2d_1d  ! called by lim_ther
   PUBLIC tab_1d_2d  ! called by lim_ther

   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limtab.F90,v 1.2 2005/03/27 18:34:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tab_2d_1d ( ndim1d, tab1d, tab2d, ndim2d_x, ndim2d_y, tab_ind )

      INTEGER, INTENT(in) :: &
         ndim1d, ndim2d_x, ndim2d_y

      REAL(wp), DIMENSION (ndim2d_x, ndim2d_y), INTENT(in) ::  &
         tab2d

      INTEGER, DIMENSION ( ndim1d), INTENT ( in) :: &
         tab_ind

      REAL(wp), DIMENSION(ndim1d), INTENT ( out) ::  & 
         tab1d

      INTEGER ::  &
         jn , jid, jjd
        
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1, ndim2d_x ) + 1
         jjd        = ( tab_ind(jn) - 1 ) / ndim2d_x + 1
         tab1d( jn) = tab2d( jid, jjd)
      END DO 

   END SUBROUTINE tab_2d_1d


   SUBROUTINE tab_1d_2d ( ndim1d, tab2d, tab_ind, tab1d, ndim2d_x, ndim2d_y )

      INTEGER, INTENT ( in) :: &
         ndim1d, ndim2d_x, ndim2d_y

      INTEGER, DIMENSION (ndim1d) , INTENT (in) :: &
         tab_ind

      REAL(wp), DIMENSION(ndim1d), INTENT (in) ::  &
         tab1d  

      REAL(wp), DIMENSION (ndim2d_x, ndim2d_y), INTENT ( out) :: &
         tab2d

      INTEGER :: &
         jn, jid, jjd

      DO jn = 1, ndim1d
         jid             = MOD( tab_ind(jn) - 1, ndim2d_x) + 1
         jjd             =    ( tab_ind(jn) - 1 ) / ndim2d_x  + 1
         tab2d(jid, jjd) = tab1d( jn)
      END DO

   END SUBROUTINE tab_1d_2d

#endif
END MODULE limtab
