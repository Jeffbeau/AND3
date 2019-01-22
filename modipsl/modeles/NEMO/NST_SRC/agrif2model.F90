#if defined key_agrif
!     ************************************************************************** 
!!!   Subroutine   Agrif_Set_numberofcells
!     ************************************************************************** 
! 
      Subroutine Agrif_Set_numberofcells(Agrif_Gr)
      USE Agrif_Types
      Implicit none
      Type(Agrif_Grid), Pointer :: Agrif_Gr
      if ( associated(Agrif_Curgrid) )then
#include "SetNumberofcells.h"
      endif
      End Subroutine Agrif_Set_numberofcells
!
!     ************************************************************************** 
!!!   Subroutine   Agrif_Get_numberofcells
!     ************************************************************************** 
      Subroutine Agrif_Get_numberofcells(Agrif_Gr)
      USE Agrif_Types
      Implicit none
      Type(Agrif_Grid), Pointer :: Agrif_Gr
#include "GetNumberofcells.h"     
      End Subroutine Agrif_Get_numberofcells
!
!     ************************************************************************** 
!!!   Subroutine Agrif_Allocationcalls
!     ************************************************************************** 
      Subroutine Agrif_Allocationcalls(Agrif_Gr)
      USE Agrif_Types 
#include "include_agrif.h"
#include "include_use_instance_agrif.h"
      Implicit none
      Type(Agrif_Grid), Pointer :: Agrif_Gr
#include "allocations_calls_agrif.h"
      End Subroutine Agrif_Allocationcalls
!
!     **************************************************************************  
!!!   Subroutine Agrif_probdim_modtype_def
!     **************************************************************************
      Subroutine Agrif_probdim_modtype_def()
      Use Agrif_Types
      Implicit none
#include "modtype_agrif.h"
#include "probdim_agrif.h"
#include "keys_agrif.h"
      Return
      End Subroutine Agrif_probdim_modtype_def
!
!     **************************************************************************  
!!!   Subroutine Agrif_clustering_def
!     **************************************************************************  
      Subroutine Agrif_clustering_def()
      Use Agrif_Types
      Implicit none
#include "clustering_agrif.h"      
      Return
      End Subroutine Agrif_clustering_def
#else
      subroutine Agrif2Model
         write(*,*) 'Impossible to bet here'
      end subroutine Agrif2model
#endif
