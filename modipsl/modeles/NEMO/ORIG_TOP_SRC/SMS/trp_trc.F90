MODULE trp_trc

#if defined key_passivetrc
   !!======================================================================
   !! Module trp_trc
   !!======================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/trp_trc.F90,v 1.4 2005/09/12 09:05:00 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! passive tracers number
   USE par_trc_trp , ONLY : &
      jptra    =>   jptra       !!: number of passive tracers

#if defined key_trc_diaadd
   USE par_trc_trp , ONLY : &
      jpdia2d  =>  jpdia2d , &  !!: number of passive tracers
      jpdia3d  =>  jpdia3d
#endif

   !! passive tracers fields 
   USE trc , ONLY :  &
      trai   =>   trai , &  !!: initial total tracer
      trb    =>   trb  , &  !!: tracer field (before)
      tra    =>   tra  , &  !!: tracer field (now)
      trn    =>   trn       !!: tracer field (after)

#if defined key_trc_diaadd
   USE trc , ONLY :  &
      trc2d   =>   trc2d , &  !!: additional 2D variable for ouputs
      trc3d   =>   trc3d      !!: additional 3D variable for ouputs
#endif
   !! time step
   USE trc , ONLY :  &
      ndttrc =>   ndttrc    !!: frequency of step on passive tracers (NAMELIST)

   !! non-centered advection scheme (smolarkiewicz)
   USE trc , ONLY : &
      rtrn   =>   rtrn      !!: value for truncation (NAMELIST)

   USE trc , ONLY : &
      ctrcnm   =>   ctrcnm      !!: value for truncation (NAMELIST)
#else
   !!======================================================================
   !!  Empty module : No passive tracer 
   !!======================================================================
#endif

END MODULE trp_trc
