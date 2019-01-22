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
/*                    OPTI_1_ajoutvarofsubroutine                             */
/******************************************************************************/
/* Firstpass 1                                                                */
/* We should complete the listvarofsubroutine                                 */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void OPTI_1_ajoutvarofsubroutine(listvar *listtoadd)
{
   listvar *tmplist;
   
   tmplist = (listvar *)NULL;
   if ( firstpass == 1 && VariableIsParameter == 0 && SaveDeclare == 0)
   {
      tmplist = duplicatelistvar(listtoadd);
      varofsubroutineliste = AddListvarToListvar
                                               (tmplist,varofsubroutineliste,1);
   }
}

/******************************************************************************/
/*                 CleanThelistvarofsubroutineliste                           */
/******************************************************************************/
/* This subroutine is to remove from the varofsubroutineliste                 */
/* all variables which are not located in the subroutine argument             */
/******************************************************************************/
void CleanThelistvarofsubroutineliste()
{
  listvar *newvar;
  listvar *newvarprec;
  listvar *tmpglobvar;
  int out;

  newvarprec = (listvar *)NULL;
  newvar = varofsubroutineliste;
  while ( newvar )
  {

     out = 0;
     tmpglobvar = listargsubroutine;
     while ( tmpglobvar && out == 0 )
     {
        if ( !strcasecmp(newvar->var->nomvar,tmpglobvar->var->nomvar) &&
             !strcasecmp(newvar->var->modulename,subroutinename) )
        {
           out = 1;
	}
	else
	{
           tmpglobvar = tmpglobvar->suiv;	
	}	
     }
     /*  if the variable has not be found we should remove it                 */
     if ( out == 0 && !strcasecmp(newvar->var->modulename,subroutinename) )
     {
        /* remove the variable in the  varofsubroutineliste                   */
	if ( newvar == varofsubroutineliste )
	{
	   varofsubroutineliste = varofsubroutineliste->suiv;
	   newvar = varofsubroutineliste;
	}
	else
	{
	   newvarprec->suiv = newvar->suiv;
	   newvar = newvarprec->suiv;
	}
     }
     else
     {
         newvarprec= newvar;
	 newvar = newvar->suiv;
     }
  }
}


/******************************************************************************/
/*             UpdatevarofsubroutinelisteWithcommonlist                       */
/******************************************************************************/
/*  This subroutines is used to add the variable defined in common in the     */
/*    varofsubroutineliste                                                    */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void UpdatevarofsubroutinelisteWithcommonlist()
{
   listvarcommon *parcours;
   listvar *parcours2;
   listvar *parcoursvar;
   listvar *parcoursvarprec;
   int out;
   
   parcoursvar = varofsubroutineliste;
   parcoursvarprec = (listvar *)NULL;
   while ( parcoursvar )
   {
      /* We should look in the commonlist if this variable is present         */
      parcours=commonlist;
      out=0;
      while( parcours && out == 0 )
      {
         if ( !strcasecmp(parcoursvar->var->nomvar,parcours->nomvar) &&
              !strcasecmp(parcoursvar->var->subroutinename,
                                           parcours->subroutinename)
            )
         {
            out = 1 ;
         }
         else
         {
            parcours = parcours->suiv;
         }
      }
      if ( out == 1 )
      {
         /* we found it                                                       */
         /* we should remove the variable from the varofsubroutineliste       */
         if ( parcoursvar == varofsubroutineliste )
         {
            varofsubroutineliste = varofsubroutineliste->suiv;
            parcoursvar = varofsubroutineliste ;
         }
         else
         {
            parcoursvarprec->suiv = parcoursvar->suiv;
            parcoursvar = parcoursvarprec->suiv;
         }
      }
      else
      {
         parcoursvarprec = parcoursvar;
         parcoursvar = parcoursvar->suiv;
      }
   }
   
   /* now we should remove all parameters                                     */
   parcoursvar = varofsubroutineliste;
   while ( parcoursvar )
   {
      /* We should look in the commonlist if this variable is present         */
      parcours2=parameterlist;
      out=0;
      while( parcours2 && out == 0 )
      {
         if ( !strcasecmp(parcoursvar->var->nomvar,parcours2->var->nomvar) &&
              !strcasecmp(parcoursvar->var->subroutinename,
                                           parcours2->var->subroutinename) 
            )
         {
            out = 1 ;
            /*                                                                */
         }
         else
         {
            parcours2 = parcours2->suiv;
         }
      }
      if ( out == 1 )
      {
         /* we did find it                                                    */
         /* we should remove the variable from the varofsubroutineliste       */
         if ( parcoursvar == varofsubroutineliste )
         {
            varofsubroutineliste = varofsubroutineliste->suiv;
            parcoursvar = varofsubroutineliste;
         }
         else
         {
            parcoursvarprec->suiv = parcoursvar->suiv;
            parcoursvar = parcoursvarprec->suiv;
         }
      }
      else
      {
         parcoursvarprec = parcoursvar;
         parcoursvar = parcoursvar->suiv;
      }
   }
}


/******************************************************************************/
/*                COM_1_UpdatevarsubroutineWithvarofsubroutinelist            */
/******************************************************************************/
/*  This subroutines is used to add the variable defined in common in the     */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void COM_1_UpdatevarsubroutineWithvarofsubroutinelist()
{
   listvar *parcours;
   listvar *parcours2;
   listvar *parcoursprec;
   int out;
   
   parcours = varsubroutine;
   while ( parcours )
   {
      /* We should look in the varofsubroutineliste if this variable is       */
      /*    present                                                           */
      parcours2=varofsubroutineliste;
      out=0;
      while( parcours2 && out == 0 )
      {
         if ( !strcasecmp(parcours->var->nomvar,parcours2->var->nomvar) &&
              !strcasecmp(parcours->var->subroutinename,
                                            parcours2->var->modulename) 
            )
         {
            out = 1 ;
         }
         else
         {
            parcours2 = parcours2->suiv;
         }
      }
      if ( out == 1 )
      {
         /* we did not find it                                                */
         /* we should remove the variable from the varsubroutine              */
         if ( parcours ==  varsubroutine)
         {
            varsubroutine = varsubroutine->suiv;
            parcours = varsubroutine;
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
}
