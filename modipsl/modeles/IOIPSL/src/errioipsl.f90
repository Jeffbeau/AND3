!$Header: /home/ioipsl/CVSROOT/IOIPSL/src/errioipsl.f90,v 2.2 2005/02/22 10:14:14 adm Exp $
!-
MODULE errioipsl
!---------------------------------------------------------------------
IMPLICIT NONE
!-
PRIVATE
!-
PUBLIC :: ipslnlf, ipslerr, histerr, ipsldbg
!-
  INTEGER :: n_l=6
  LOGICAL :: ioipsl_debug=.FALSE.
!-
!===
CONTAINS
!===
SUBROUTINE ipslnlf (new_number,old_number)
!!--------------------------------------------------------------------
!! The "ipslnlf" routine allows to know and modify
!! the current logical number for the messages,
!!
!! SUBROUTINE ipslnlf (new_number,old_number)
!!
!! Optional INPUT argument
!!
!! (I) new_number : new logical number of the file
!!
!! Optional OUTPUT argument
!!
!! (I) old_number : current logical number of the file
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,OPTIONAL,INTENT(IN)  :: new_number
  INTEGER,OPTIONAL,INTENT(OUT) :: old_number
!---------------------------------------------------------------------
  IF (PRESENT(old_number)) THEN
    old_number = n_l
  ENDIF
  IF (PRESENT(new_number)) THEN
    n_l = new_number
  ENDIF
!---------------------
END SUBROUTINE ipslnlf
!===
SUBROUTINE ipslerr (plev,pcname,pstr1,pstr2,pstr3)
!---------------------------------------------------------------------
!! The "ipslerr" routine
!! allows to handle the messages to the user.
!!
!! INPUT
!!
!! plev   : Category of message to be reported to the user
!!          1 = Note to the user
!!          2 = Warning to the user
!!          3 = Fatal error
!! pcname : Name of subroutine which has called ipslerr
!! pstr1   
!! pstr2  : Strings containing the explanations to the user
!! pstr3
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER :: plev
   CHARACTER(LEN=*) :: pcname,pstr1,pstr2,pstr3
!-
   CHARACTER(LEN=30),DIMENSION(3) :: pemsg = &
  &  (/ "NOTE TO THE USER FROM ROUTINE ", &
  &     "WARNING FROM ROUTINE          ", &
  &     "FATAL ERROR FROM ROUTINE      " /)
!---------------------------------------------------------------------
   IF ( (plev >= 1).AND.(plev <= 3) ) THEN
     WRITE(n_l,'(/,A," ",A)') TRIM(pemsg(plev)),TRIM(pcname)
     WRITE(n_l,'(3(" --> ",A,/))') TRIM(pstr1),TRIM(pstr2),TRIM(pstr3)
   ENDIF
   IF (plev == 3) THEN
     STOP 'Fatal error from IOIPSL. See stdout for more details'
   ENDIF
!---------------------
END SUBROUTINE ipslerr
!===
SUBROUTINE histerr (plev,pcname,pstr1,pstr2,pstr3)
!---------------------------------------------------------------------
!- INPUT
!- plev   : Category of message to be reported to the user
!-          1 = Note to the user
!-          2 = Warning to the user
!-          3 = Fatal error
!- pcname : Name of subroutine which has called histerr
!- pstr1   
!- pstr2  : String containing the explanations to the user
!- pstr3
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER :: plev
   CHARACTER(LEN=*) :: pcname,pstr1,pstr2,pstr3
!-
   CHARACTER(LEN=30),DIMENSION(3) :: pemsg = &
  &  (/ "NOTE TO THE USER FROM ROUTINE ", &
  &     "WARNING FROM ROUTINE          ", &
  &     "FATAL ERROR FROM ROUTINE      " /)
!---------------------------------------------------------------------
   IF ( (plev >= 1).AND.(plev <= 3) ) THEN
     WRITE(*,'("     ")')
     WRITE(*,'(A," ",A)') TRIM(pemsg(plev)),TRIM(pcname)
     WRITE(*,'(" --> ",A)') pstr1
     WRITE(*,'(" --> ",A)') pstr2
     WRITE(*,'(" --> ",A)') pstr3
   ENDIF
   IF (plev == 3) THEN
     STOP 'Fatal error from IOIPSL. See stdout for more details'
   ENDIF
!---------------------
END SUBROUTINE histerr
!===
SUBROUTINE ipsldbg (new_status,old_status)
!!--------------------------------------------------------------------
!! The "ipsldbg" routine
!! allows to activate or deactivate the debug,
!! and to know the current status of the debug.
!!
!! SUBROUTINE ipsldbg (new_status,old_status)
!!
!! Optional INPUT argument
!!
!! (L) new_status : new status of the debug
!!
!! Optional OUTPUT argument
!!
!! (L) old_status : current status of the debug
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  LOGICAL,OPTIONAL,INTENT(IN)  :: new_status
  LOGICAL,OPTIONAL,INTENT(OUT) :: old_status
!---------------------------------------------------------------------
  IF (PRESENT(old_status)) THEN
    old_status = ioipsl_debug
  ENDIF
  IF (PRESENT(new_status)) THEN
    ioipsl_debug = new_status
  ENDIF
!---------------------
END SUBROUTINE ipsldbg
!===
!-------------------
END MODULE errioipsl
