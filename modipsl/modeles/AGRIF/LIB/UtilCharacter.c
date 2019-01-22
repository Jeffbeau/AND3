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
/*                         FindAndChangeNameToTabvars                         */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
/* if  whichone = 0 ----> Agrif_tabvars(i) % var % array2                     */
/*                                                                            */
/* if  whichone = 1 ----> Agrif_tabvars(i) % parentvar % var % array2         */
/*                                                                            */
/******************************************************************************/
void FindAndChangeNameToTabvars(char name[LONGNOM],char toprint[LONGNOM],
                                listvar * listtosee, int whichone)
{
   listvar *newvar;
   int out;
   
   if ( strcasecmp(name,"") )
   {
      newvar=listtosee;
      out=0;
      while( newvar && out == 0 )
      {
         if ( !strcasecmp(newvar->var->nomvar,name) )
         {
            out = 1;
            strcat(toprint,vargridcurgridtabvars(newvar->var,whichone));
         }
         else newvar=newvar->suiv;
      }
      if ( out == 0 ) strcat(toprint,name);
   }
}


/******************************************************************************/
/*                     ChangeTheInitalvaluebyTabvarsName                      */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
char *ChangeTheInitalvaluebyTabvarsName(char *nom,listvar *listtoread, int whichone)
{
   char toprinttmp[LONGNOM];
   int i;
   char chartmp[2];
   
   i=0;
   strcpy(toprintglob,"");
   strcpy(toprinttmp,"");
   /*                                                                         */
   while ( i < strlen(nom) )
   {
      if ( nom[i] == '+' ) 
      {
         FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
         strcpy(toprinttmp,"");
         strcat(toprintglob,"+");
      }
      else if ( nom[i] == '-' ) 
      {
         FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
         strcpy(toprinttmp,"");
         strcat(toprintglob,"-");
      }
      else if ( nom[i] == '*' )
      {
         FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
         strcpy(toprinttmp,"");
         strcat(toprintglob,"*");
      }
      else if ( nom[i] == '/' )
      {
         FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
         strcpy(toprinttmp,"");
         strcat(toprintglob,"/");
      }
      else if ( nom[i] == '(' )
      {
         FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
         strcpy(toprinttmp,"");
         strcat(toprintglob,"(");
      }
      else if ( nom[i] == ')' )
      {
         FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
         strcpy(toprinttmp,"");
         strcat(toprintglob,")");
      }
      else if ( nom[i] == ':' )
      {
         FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
         strcpy(toprinttmp,"");
         strcat(toprintglob,":");
      }
      else if ( nom[i] == ',' )
      {
         FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
         strcpy(toprinttmp,"");
         strcat(toprintglob,",");
      }
      else
      {
         sprintf(chartmp,"%c",nom[i]);        
         strcat(toprinttmp,chartmp);
      }
      /*                                                                      */
      i=i+1;
   }
   FindAndChangeNameToTabvars(toprinttmp,toprintglob,listtoread,whichone);
   strcpy(toprinttmp,"");
   
   /*                                                                         */
   return toprintglob;
}


void IsVarInUseFile(char *nom)
{
   listvar *parcours;
   listparameter *parcoursparam;
   int out;

   out = 0;

   parcours = globliste;
   while( parcours && out == 0 )
   {
      if ( !strcasecmp(nom,parcours->var->nomvar) ) out =1 ;
     else parcours=parcours->suiv;
   }
   if ( out == 0 )
   {
      parcours = globparam;
      while( parcours && out == 0 )
      {
         if ( !strcasecmp(nom,parcours->var->nomvar) ) out =1 ;
        else parcours=parcours->suiv;
      }
   }
   if ( out == 0 )
   {
      parcours = parameterlist;
      while( parcours && out == 0 )
      {
         if ( !strcasecmp(nom,parcours->var->nomvar) ) out =1 ;
        else parcours=parcours->suiv;
      }
   }
   if ( out == 0 )
   {
      parcoursparam = tmpparameterlocallist;
      while( parcoursparam && out == 0 )
      {
         if ( !strcasecmp(nom,parcoursparam->name) ) out =2 ;
         else parcoursparam=parcoursparam->suiv;
      }
   }
   if ( out == 0 )
   {
      parcours = globalvarofusefile;
      while( parcours && out == 0 )
      {
         if ( !strcasecmp(nom,parcours->var->nomvar) ) out =2 ;
        else parcours=parcours->suiv;
      }
   }
   if ( out == 0 || out == 2 )
   {
      parcoursparam = tmpparameterlocallist2;
      while( parcoursparam && out != 1 )
      {
         if ( !strcasecmp(nom,parcoursparam->name) ) out =1 ;
         else parcoursparam=parcoursparam->suiv;
      }
      if ( out == 1 ) 
      {
         strcpy(charusemodule,parcoursparam->modulename);
         Addmoduletothelist(parcoursparam->modulename);
      }
   }
   if ( out == 0 ) printf("--- in UtilCharacter we do not found the \n");
   if ( out == 0 ) printf("---  variable %s, the module where this \n",nom);
   if ( out == 0 ) printf("---  variable has been defined has not been\n");
   if ( out == 0 ) printf("---  found.\n");
}

/******************************************************************************/
/*                      DecomposeTheNameinlistnom                             */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/******************************************************************************/
listnom *DecomposeTheNameinlistnom(char *nom, listnom * listout)
{
   char toprinttmp[LONGNOM];
   int i;
   char chartmp[2];
   
   i=0;
   strcpy(toprinttmp,"");
   /*                                                                         */
   while ( i < strlen(nom) )
   {
      if ( nom[i] == '+' ||
           nom[i] == '-' ||
           nom[i] == '*' ||
           nom[i] == '/' ||
           nom[i] == ')' ||
           nom[i] == '(' ||
           nom[i] == ',' ||
           nom[i] == ':' 
         ) 
      {
         if (strcasecmp(toprinttmp,"") && ( toprinttmp[0] >= 'A' ) )
         { 
             listout = Addtolistnom(toprinttmp,listout);
             
         }
         strcpy(toprinttmp,"");
      }
      else
      {
         sprintf(chartmp,"%c",nom[i]);        
         strcat(toprinttmp,chartmp);
      }
      /*                                                                      */
      i=i+1;
   }
   if (strcasecmp(toprinttmp,"") && ( toprinttmp[0] >= 'A' ) ) 
   { 
      listout = Addtolistnom(toprinttmp,listout);
   }
   strcpy(toprinttmp,"");
 
   return listout;   
}


/******************************************************************************/
/*                      DecomposeTheName                                      */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/*               Agrif_<toto>(variable) ====>     Agrif_<toto>(variable)      */
/*                                                                            */
/******************************************************************************/
void DecomposeTheName(char *nom)
{
   char toprinttmp[LONGNOM];
   int i;
   char chartmp[2];
   
   i=0;
   strcpy(toprinttmp,"");
   /*                                                                         */
   while ( i < strlen(nom) )
   {
      if ( nom[i] == '+' ||
           nom[i] == '-' ||
           nom[i] == '*' ||
           nom[i] == '/' ||
           nom[i] == ')' ||
           nom[i] == '(' ||
           nom[i] == ',' ||
           nom[i] == ':' 
         ) 
      {
         if (strcasecmp(toprinttmp,"") && ( toprinttmp[0] >= 'A' ) )
         { 
            ajoutevarindoloop_definedimension (toprinttmp);
            /* Is this variable present in globvarofusefile                   */
            IsVarInUseFile(toprinttmp);
         }
         strcpy(toprinttmp,"");
      }
      else
      {
         sprintf(chartmp,"%c",nom[i]);        
         strcat(toprinttmp,chartmp);
      }
      /*                                                                      */
      i=i+1;
   }
   if (strcasecmp(toprinttmp,"") && ( toprinttmp[0] >= 'A' ) ) 
   { 
      ajoutevarindoloop_definedimension (toprinttmp);
      /* Is this variable present in globvarofusefile                         */
      IsVarInUseFile(toprinttmp);
   }
   strcpy(toprinttmp,"");
   
}
