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
/*                            initdimprob                                     */
/******************************************************************************/
/* This subroutine is used to initialized grid dimension variable             */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void initdimprob(int dimprobmod, char * nx, char * ny,char* nz)
{
  dimprob = dimprobmod;

  strcpy(nbmaillesX,nx);
  strcpy(nbmaillesY,ny);
  strcpy(nbmaillesZ,nz);
}

/******************************************************************************/
/*                      Variableshouldberemove                                */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/*               Agrif_<toto>(variable) ====>     Agrif_<toto>(variable)      */
/*                                                                            */
/******************************************************************************/
int Variableshouldberemove(char *nom)
{

   int remove;
   
   remove = 0 ; 
   
   if ( remove == 0 && !strcasecmp(nom,"RESHAPE") ) remove = 1 ;  
   if ( remove == 0 && AGRIF_n_Agrif_in_Tok_NAME(nom) == 1 ) remove = 1 ;  

   return remove;   
}

/******************************************************************************/
/*                          variableisglobal                                  */
/******************************************************************************/
/* This subroutine is to know if a variable is global                         */
/******************************************************************************/
int variableisglobal(listvar *curvar, listvar *listin)
{
  int Globalite;
  listvar *newvar;


  Globalite = 0;
  newvar = listin;
  while ( newvar && Globalite == 0 )
  {
     if ( !strcasecmp(newvar->var->nomvar,curvar->var->nomvar) )
     {
        Globalite = 1;
        /* Now we should give the definition of the variable in the           */
        /* table listvarindoloop                                              */
        strcpy(curvar->var->typevar,newvar->var->typevar);
        strcpy(curvar->var->dimchar,newvar->var->dimchar);
        curvar->var->nbdim = newvar->var->nbdim;
        curvar->var->dimensiongiven = newvar->var->dimensiongiven;
        curvar->var->typegiven = newvar->var->typegiven;
        curvar->var->allocatable = newvar->var->allocatable;
        curvar->var->pointerdeclare = newvar->var->pointerdeclare;
        curvar->var->indicetabvars = newvar->var->indicetabvars;
        strcpy(curvar->var->precision,newvar->var->precision);
        strcpy(curvar->var->readedlistdimension,
                                              newvar->var->readedlistdimension);
     }
     else
     {
         newvar = newvar->suiv;
     }
  }

  return Globalite ;
}

/******************************************************************************/
/*                     variableisparameterglobal                              */
/******************************************************************************/
/* This subroutine is to know if a variable is global                         */
/******************************************************************************/
int variableisparameterglobal(listvar *curvar, listparameter *listin)
{
  int Globalite;
  listparameter *newvar;

  Globalite = 0;
  newvar = listin;
  while ( newvar && Globalite == 0 )
  {
     if ( !strcasecmp(newvar->name,curvar->var->nomvar) ) Globalite = 1;
     else newvar = newvar->suiv;
  }

  return Globalite ;
}


/******************************************************************************/
/*                 addi_1_addsubroutine_inst_back_alloc                       */
/******************************************************************************/
/* Firstpass 0                                                                */
/* We should add subroutine of instanciation, back instaciation and           */
/* allocation                                                                 */
/* if moduleorcontains = 1 we are at the end module keyword                   */
/* if moduleorcontains = 0 we are at the contains   keyword                   */
/******************************************************************************/
void addi_0_addsubroutine_inst_back_alloc(int moduleorcontains)
{
   char ligne[LONGNOM];
   int Allocisempty;
   listvar *newvar;
   

   if ( firstpass == 0)
   {
     /* It is necessary to know if this subroutine is not empty               */
     Allocisempty = 0;
     newvar = globliste;
     while ( newvar && Allocisempty == 0 )
     {
        if ( !strcasecmp(newvar->var->modulename,curmodulename)) Allocisempty=1;
        else newvar = newvar->suiv;
     }
     if ( Allocisempty == 1 )
     {
         while ( newvar &&
                 !strcasecmp(newvar->var->modulename,curmodulename) &&
                 Allocisempty == 1 )
         {
            if ( (newvar->var->nbdim !=0          &&
                  newvar->var->allocatable != 1 ) ||
                 (newvar->var->nbdim == 0         &&
                  strcmp(newvar->var->initialvalue,"")) ) Allocisempty = 0;
            else newvar = newvar->suiv;
         }
     }
     if ( Allocisempty == 0 )
     {
      if ( MOD_n_InstanceInModule() == 1)
      {
         /* we should remove end module <name>                                */
         if ( moduleorcontains == 1 )
         {
            RemoveWordCUR(fortranout,(long)(-strlen(curmodulename)-12),
                                          strlen(curmodulename)+11);
         }
         /* we should remove contains                                         */
         if ( moduleorcontains == 0 )
         {
            RemoveWordCUR(fortranout,(long)(-8),8);
         }
         strcpy (ligne, "\n      PUBLIC Alloc_agrif_");
         strcat (ligne, curmodulename);
         strcat (ligne, "\n");
         fprintf(fortranout,ligne);
      }
      if (MOD_n_InstanceInModule() == 1)
      {      
         fprintf(fortranout,"\n      contains\n"); 
         strcpy (ligne, "\n#include \"alloc_agrif_");
         strcat (ligne, curmodulename);
         strcat (ligne, ".h\"\n");
         fprintf(fortranout,ligne);
         /* On reecrit la mot cle end module qui a ete efface du fichier      */
         /*    d'origine                                                      */
         if ( moduleorcontains == 1 ) fprintf(fortranout,"\n      end module %s"
                                                                ,curmodulename);
      }
     }
   }
}


/******************************************************************************/
/*                          OPTI_0_IsTabvarsUseInArgument                     */
/******************************************************************************/
/* Firstpass 1                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int OPTI_0_IsTabvarsUseInArgument()
{
   int out;
   int doloopout;
   listvar *parcours;   

   out=1;
  
   if ( listvarindoloop )
   {
      doloopout = 0;
      parcours = listvarindoloop;
      while ( parcours && doloopout == 0 )   
      {
         if ( !strcasecmp(parcours->var->modulename,subroutinename) ) 
                                                                  doloopout = 1;
         else parcours = parcours->suiv;
      }
      if (  doloopout == 0 ) out = 0;
      else out = 1 ;
   }
   else out = 0;

   return out;
}


/******************************************************************************/
/*                        ImplicitNoneInSubroutine                            */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int ImplicitNoneInSubroutine()
{
  listname *parcours;
  int out;

  parcours= listimplicitnone;
  out = 0 ;
  while ( parcours && out == 0 )
  {
     if ( !strcasecmp(parcours->name,subroutinename) ) out = 1;
     else parcours = parcours->suiv;
  
  }
  return out;
}

/******************************************************************************/
/*                          OPTI_0_varispointer                               */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int OPTI_0_varispointer(char *ident)
{
   listvar *newvar;
   int out;

   out =0;
   if (firstpass == 0 )
   {
      newvar = globalvarofusefile;
      while ( newvar && out == 0 )
      {
         if ( !strcmp(ident,newvar->var->nomvar) && 
              newvar->var->pointerdeclare == 1 )  out = 1;
         else newvar = newvar->suiv;
      }
   }
   return out;
}
