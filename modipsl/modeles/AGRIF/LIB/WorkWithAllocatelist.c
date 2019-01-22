/******************************************************************************/
/*                                                                            */
/*     CONV (converter) for Agrif (Adaptive Grid Refinement In Fortran)       */
/*                                                                            */
/*     Copyright (C) 2005 Laurent Debreu (Laurent.Debreu@imag.fr)             */
/*                        Cyril Mazauric (Cyril.Mazauric@imag.fr)             */
/*                                                                            */
/*     This program is free software; you can redistribute it and/or modify   */
/*    it                                                                      */
/*                                                                            */
/*    This program is distributed in the hope that it will be useful,         */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of         */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          */
/*    GNU General Public License for more details.                            */
/*                                                                            */
/******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "decl.h"


/******************************************************************************/
/*                       OPTI_1_AddIdentToTheAllocateList                     */
/******************************************************************************/
/* Firstpass 1                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_1_AddIdentToTheAllocateList(char *nom)
{
   listallocate *newvar;
   listallocate *parcours;
   int out;

   if ( firstpass == 1 ) 
   {
      if ( !AllocateList )
      {
         newvar = (listallocate *)malloc(sizeof(listallocate));
         strcpy(newvar->nomvar,nom);
         strcpy(newvar->subroutine,subroutinename);
         strcpy(newvar->module,curmodulename);
         newvar->suiv = NULL;
         AllocateList = newvar;
      }
      else
      {
         parcours = AllocateList;
         out = 0 ; 
         while ( parcours->suiv && out == 0 )
         {
            if (  !strcasecmp(parcours->nomvar,nom) &&
                  !strcasecmp(parcours->subroutine,subroutinename) &&
                  !strcasecmp(parcours->module,curmodulename) ) out = 1;
            else
               parcours=parcours->suiv;               
         }
         if ( out == 0 ) 
         {
            if (  !strcasecmp(parcours->nomvar,nom) &&
                  !strcasecmp(parcours->subroutine,subroutinename) &&
                  !strcasecmp(parcours->module,curmodulename) ) out = 1;
            else
            {
               /* add the record                                              */
              newvar = (listallocate *)malloc(sizeof(listallocate));
              strcpy(newvar->nomvar,nom);
              strcpy(newvar->subroutine,subroutinename);
              strcpy(newvar->module,curmodulename);
              newvar->suiv = NULL;
              parcours->suiv = newvar;
            }
         }
      }
   }
}

/******************************************************************************/
/*                       OPTI_0_IsAllocateInThisSubroutine                    */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int OPTI_0_IsAllocateInThisSubroutine()
{
   listallocate *parcours;
   int out;

   out = 0 ; 
   if ( firstpass == 0 ) 
   {
      parcours = AllocateList;
      while ( parcours && out == 0 )
      {
         if ( !strcasecmp(parcours->subroutine,subroutinename)  )
         {
            out = 1 ;
         }
         else parcours=parcours->suiv;               
      }
   }
   return out;
}

/******************************************************************************/
/*                       OPTI_0_IsVarAllocatable                              */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int OPTI_0_IsVarAllocatable(char *ident)
{
   listallocate *parcours;
   int out;

   out = 0 ; 
   if ( firstpass == 0 ) 
   {
      parcours = AllocateList;
      while ( parcours && out == 0 )
      {
         if ( !strcasecmp(parcours->nomvar,ident)  ) out = 1 ;
         else parcours=parcours->suiv;               
      }
   }
   return out;
}


/******************************************************************************/
/*                     OPTI_0_varisallocatable                                */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int OPTI_0_varisallocatable(char *ident)
{
   listvar *newvar;
   listallocate *newvaralloc;
   int out;

   out =0;
   if (firstpass == 0 )
   {
      newvar = globalvarofusefile;
      while ( newvar && out == 0 )
      {
         if ( !strcmp(ident,newvar->var->nomvar) && 
              newvar->var->allocatable == 1 )  out = 1;
         else newvar = newvar->suiv;
      }
      if ( out == 0 )
      {
         newvar = globliste;
         while ( newvar && out == 0 )
         {
            if ( !strcmp(ident,newvar->var->nomvar) && 
                 newvar->var->allocatable == 1 )  out = 1;
            else newvar = newvar->suiv;
         }      
      }
      if ( out == 0 )
      {
         newvaralloc = AllocateList;
         while ( newvar && out == 0 )
         {
            if ( !strcasecmp(ident,newvaralloc->nomvar) && 
                 !strcasecmp(newvaralloc->subroutine,subroutinename)
               )  out = 1;
            else newvaralloc = newvaralloc->suiv;
         }      
      }
   }
   return out;
}
