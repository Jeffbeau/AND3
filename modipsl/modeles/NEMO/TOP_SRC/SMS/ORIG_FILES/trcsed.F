CC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/trcsed.F,v 1.5 2005/11/14 12:30:28 opalod Exp $
! CDIR$ LIST
      SUBROUTINE trcsed(kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcsed
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     compute the now trend due to the vertical sedimentation of
CCC     detritus and add it to the general trend of detritus equations.
CCC
CCC
CC   METHOD :
CC   -------
CC      this ROUTINE compute not exactly the advection but the
CC      transport term, i.e.  dz(wt) and dz(ws)., dz(wtr)
CC      using an upstream scheme
CC
CC	the now vertical advection of tracers is given by:
CC
CC      	dz(trn wn) = 1/bt dk+1( e1t e2t vsed (trn) )
CC
CC	add this trend now to the general trend of tracer (ta,sa,tra):
CC
CC                     tra = tra + dz(trn wn)
CC
CC      IF 'key_trc_diabio' key is activated, the now vertical advection
CC      trend of passive tracers is saved for futher diagnostics.
CC
CC	multitasked on vertical slab (jj-loop)
CC
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC            /cottrp/          : passive tracer fields
CC            /comtsk/          : multitasking
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cottrp/tra       : general tracer trend increased by the
CC				  now vertical tracer advection trend
CC            /cottbd/ trbio    : now vertical passive tracer advection
CC                                trend
CC                                (IF 'key_trc_diabio' key is activated)
CC
CC   WORKSPACE :
CC   ---------
CC	local
CC		ze1e2w, ze3tr, ztra
CC      COMMON
CC
CC   EXTERNAL :                   no
CC   --------
CC
CC   REFERENCES :                 no
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC       original : 95-06 (M. Levy)
CC       additions: 00-12 (E. Kestenare): clean up 
CC----------------------------------------------------------------------
! CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      USE lbclnk
      IMPLICIT NONE
! CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kt

#if defined key_passivetrc && defined key_trc_lobster1

      INTEGER ji,jj,jk
      REAL ze3tr,ztra
      REAL zwork(jpi,jpk),ze3t(jpk)
CC----------------------------------------------------------------------
CC statement functions
CC ===================
! CDIR$ NOLIST
#include "domzgr_substitute.h90"
! CDIR$ LIST
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
#     if defined key_trc_diaadd
      DO jk=1,jpkbm1
        ze3t(jk)=e3t(jk)*86400.
      END DO 
      DO jk=jpkb,jpk
        ze3t(jk)=0.
      END DO 
#     endif
C
C vertical slab
C =============
C
      DO 1000 jj = 1,jpj

C
C
C 1. sedimentation of detritus  : upstream scheme
C -----------------------------------------------
C
C
C for detritus sedimentation only - jpdet
C
C 1.1 initialisation needed for bottom and surface value
C
              DO jk=1,jpk
                DO  ji = 1,jpi
                  zwork(ji,jk) = 0.
                END DO 
              END DO
C
C 1.2 tracer flux at w-point: we use -vsed (downward flux)
C with simplification : no e1*e2
C
              DO  jk = 2,jpk
                DO  ji = 1,jpi
                  zwork(ji,jk) = -vsed * trn(ji,jj,jk - 1,jpdet)
                END DO
              END DO
C
C 1.3 tracer flux divergence at t-point added to the general trend
C
              DO  jk = 1,jpkm1
                DO  ji = 1,jpi
                  ze3tr = 1./fse3t(ji,jj,jk)
                  ztra = -ze3tr * (zwork(ji,jk) - zwork(ji,jk + 1))
                  tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + ztra
#            if defined key_trc_diabio
                  trbio(ji,jj,jk,8) = ztra
#            endif
#if defined key_trc_diaadd
            trc2d(ji,jj,8)=trc2d(ji,jj,8)+ztra*ze3t(jk)
#endif
                END DO
              END DO
C
C END of slab
C ===========

 1000 CONTINUE
C
#if defined key_trc_diabio
C Lateral boundary conditions on trcbio
      CALL lbc_lnk (trbio(:,:,1,8), 'T', 1. )
#endif
#if defined key_trc_diaadd
C Lateral boundary conditions on trc2d
      CALL lbc_lnk (trc2d(:,:,8), 'T', 1. )
#endif
C

#else
C
C     no passive tracer
C
#endif
C
      RETURN
      END
