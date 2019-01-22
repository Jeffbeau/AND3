!!DB: 2009.08.31 -- Removed ORCA-code
MODULE diafwb
   !!======================================================================
   !!                       ***  MODULE  diafwb  ***
   !! Ocean diagnostics: freshwater budget
   !!======================================================================
!#if ( defined key_orca_r2 || defined  key_orca_r4 ) && ! defined key_dynspg_rl && ! defined key_coupled

   LOGICAL, PUBLIC, PARAMETER ::   lk_diafwb = .FALSE.    !: fresh water budget flag
CONTAINS
   SUBROUTINE dia_fwb( kt )        ! Empty routine
     USE in_out_manager
      
     if(lwp .AND. kt==nit000) write(numout,*) 'DB: dia_fwb: disabled '

   END SUBROUTINE dia_fwb
END MODULE diafwb
