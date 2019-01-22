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
/*                       COM_1_AddvartoParamlist                              */
/******************************************************************************/
/*  This subroutines is used to add the variable defined in common in the     */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void COM_1_AddvartoParamlist(listvar *listin)
{
   listvar *parcours;
   
   if ( firstpass == 1 )
   {
      if ( !parameterlist )
      {
         parameterlist = listin;
      }
      else
      {
         parcours = parameterlist;
         while (parcours->suiv) parcours=parcours->suiv;
      
         parcours->suiv = listin;
      }
   }
}

/******************************************************************************/
/*                   COM_1_UpdateparameterlistWithlistvarindoloop             */
/******************************************************************************/
/*  This subroutines is used to add the variable defined in common in the     */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void COM_1_UpdateparameterlistWithlistvarindoloop()
{
   listvar *parcours;
   listvar *parcours2;
   listvar *parcours3;
   listvar *parcoursprec;
   int out;
   
   parcours = parameterlist;
   while ( parcours )
   {
   if ( !strcasecmp(parcours->var->subroutinename,subroutinename) )
   { 
      /* We should look in the listvarindoloop if this variable is present    */
      parcours2=listvarindoloop;
      out=0;
      while( parcours2 && out == 0 )
      {
         if ( !strcasecmp(parcours->var->nomvar,parcours2->var->nomvar) &&
              !strcasecmp(parcours->var->subroutinename,
                                            parcours2->var->modulename) 
            )
         {
            parcours->var->VariableIsParameter = 1;
            /* we should find in the globliste the type of this variable      */
            parcours3 = globliste;
            while ( parcours3 && out == 0 )
            {
               if ( !strcasecmp(parcours3->var->nomvar,parcours->var->nomvar) )
               {
                  out = 1 ;
                  strcpy(parcours->var->typevar,parcours3->var->typevar);
               }
               else
               {
                  parcours3 = parcours3->suiv;
               }
            }            
            out = 1 ;
         }
         else
         {
            parcours2 = parcours2->suiv;
         }
      }
      if ( out == 0 )
      {
         /* we did not find it                                                */
         /* we should remove the variable from the globliste                  */
         if ( parcours ==  parameterlist)
         {
            parameterlist = parameterlist->suiv;
            parcours = parameterlist;
         }
         else
         {
            parcoursprec->suiv = parcours->suiv;
            parcours = parcoursprec->suiv;
         }
      }
      else
      {
         parcoursprec = parcours;
         parcours = parcours->suiv;
      }
   }
   else
   {
      parcoursprec = parcours;
      parcours = parcours->suiv;
   }
   }
}
