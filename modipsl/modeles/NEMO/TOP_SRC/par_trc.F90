MODULE par_trc
   !!======================================================================
   !!                        ***  par_trc  ***
   !! passive tracers :   set the passive tracers parameters
   !!======================================================================
   !! History :
   !!   8.2  !  96-01  (M. Levy)  Original code
   !!        !  99-07  (M. Levy)  for LOBSTER1 or NPZD model
   !!        !  00-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!   9.0  !  04-03  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------
   !!  TOP 1.0,  LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/par_trc.F90,v 1.4 2005/09/12 09:04:53 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
#if defined key_passivetrc

   USE par_trc_trp

   IMPLICIT NONE
   PUBLIC


#if defined key_trc_diatrd

!! number of dynamical trends
#  if defined key_trc_ldfeiv
!! we keep 3 more trends for eddy induced flux (gent velocity)
#    if defined key_trcdmp
   INTEGER , PARAMETER :: jpdiatrc = 10
#    else
   INTEGER , PARAMETER :: jpdiatrc = 9
#    endif
#  else
#    if defined key_trcdmp
   INTEGER , PARAMETER :: jpdiatrc = 7
#    else
   INTEGER , PARAMETER :: jpdiatrc = 6
#    endif
#  endif
# endif

#else
   !!======================================================================
   !!  Empty module : No passive tracer 
   !!======================================================================
#endif

END MODULE par_trc
