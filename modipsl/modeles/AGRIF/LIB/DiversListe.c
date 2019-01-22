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
/*                      COM_1_AddCommonvartolist                              */
/******************************************************************************/
/*  This subroutines is used to add the variable defined in common in the     */
/*     commonlist                                                             */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void COM_1_AddCommonvartolist()
{
   listvarcommon *newvar;
   
   newvar = (listvarcommon *)malloc(sizeof(listvarcommon));
   strcpy(newvar->nomvar,commonvar);
   strcpy(newvar->commonname,commonblockname);
   strcpy(newvar->subroutinename,subroutinename);
   newvar->positioninblock= positioninblock;

   newvar->suiv = NULL;

   if ( !commonlist )
   {
      commonlist = newvar;
   }
   else
   {
      newvar->suiv = commonlist;
      commonlist = newvar;
   }
}

/******************************************************************************/
/*                     Addtolistnom                                           */
/******************************************************************************/
/* This subroutine is used to add a variable to the list                      */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
listnom *Addtolistnom(char *nom, listnom *listin)
{
   listnom *newnom;
   listnom *parcours;
   int out;

   if ( !listin )
   {
      newnom=(listnom *) malloc (sizeof (listnom));
      strcpy(newnom->nom,nom);
      newnom->suiv = NULL;
      listin = newnom;
   }
   else
   {
      parcours = listin;
      out = 0 ;
      while ( parcours && out == 0 )
      {
         if ( !strcasecmp(parcours->nom,nom) ) out = 1 ;
         else parcours=parcours->suiv;
      }
      if ( out == 0 ) 
      {
          newnom=(listnom *) malloc (sizeof (listnom));
          strcpy(newnom->nom,nom);
          newnom->suiv = listin;
          listin = newnom;
      }
   }
   return listin;
}

/******************************************************************************/
/*                           Add_listname                                     */
/******************************************************************************/
/* This subroutine is used to add a        variable to the list               */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ glob +--->+ glob +--->+ glob +--->+ glob +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
listname *Add_listname(char *nom,listname *input)
{
   listname *newnom;
   listname *parcours;
   int out;

   if ( !input )
   {
      newnom=(listname *) malloc (sizeof (listname));
      strcpy(newnom->name,nom);
      newnom->suiv = NULL;
      input = newnom;
   }
   else
   {
      parcours = input;
      out = 0 ;
      while ( parcours && out == 0 )
      {
         if ( !strcasecmp(parcours->name,nom) ) out = 1;
         else parcours=parcours->suiv;         
      }
      if ( out == 0 )
      {
         newnom=(listname *) malloc (sizeof (listname));
         strcpy(newnom->name,nom);
         newnom->suiv = input;
         input = newnom;
      }
   }
   return input;
}

/******************************************************************************/
/*                     Add_ModuleTo_listofmodules                             */
/******************************************************************************/
/* This subroutine is used to add a        variable to the list               */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ glob +--->+ glob +--->+ glob +--->+ glob +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void Add_ModuleTo_listofmodules(char *nom)
{
   listnom *newnom;

   newnom=(listnom *) malloc (sizeof (listnom));
   strcpy(newnom->nom,nom);
   newnom->suiv = listofmodules;
   listofmodules = newnom;   
}

/******************************************************************************/
/*                    ModuleIsDefineInInputFile                               */
/******************************************************************************/
/* This subroutine is used to know if the module is defined in the input file */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/******************************************************************************/
int ModuleIsDefineInInputFile(char *name)
{
   listnom *newnom;
   int out;
   
   out = 0;
   if ( listofmodules ) 
   {
      newnom = listofmodules;
      while( newnom && out == 0 )
      {
         if ( !strcasecmp(newnom->nom,name) ) out = 1 ;
         else newnom=newnom->suiv;
      }
   }
   return out;
}

/******************************************************************************/
/*                         AddNameToListNamelist                              */
/******************************************************************************/
/* This subroutine is used to add a listvar l at the end of a listvar         */
/* glob.                                                                      */
/*                                                                            */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + glob +--->+ glob +--->+ glob +--->+ glob +--->+  l   +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/******************************************************************************/
void AddNameToListNamelist(char * name)
{
   listnamelist *newvar;

   if ( strcasecmp(name,"") )
   {
      newvar =(listnamelist*)malloc(sizeof(listnamelist));
      strcpy(newvar->name,name);
      newvar->suiv = listenamelist;
      listenamelist = newvar;
   }
}

/******************************************************************************/
/*                      Addmoduletothelisttmp                                 */
/******************************************************************************/
/* This subroutine is used to add a record to a list of struct                */
/* listusemodule                                                              */
/******************************************************************************/
/*                                                                            */
/*       subroutine sub ... USE mod1 ===> insert in list                      */
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ list +--->+ list +--->+ list +--->+ list +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*       list =  listofmoduletmp                                              */
/*                                                                            */
/******************************************************************************/
void Addmoduletothelisttmp(char *name)
{
  listusemodule *newmodule;
  listusemodule *parcours;
  int out;

  if ( !listofmoduletmp)
  {
    newmodule =(listusemodule *)malloc(sizeof(listusemodule));
    strcpy(newmodule->usemodule,name);
    strcpy(newmodule->cursubroutine,subroutinename);  
    newmodule->suiv = NULL;
    listofmoduletmp = newmodule ;
  }
  else
  {
    parcours = listofmoduletmp;
    out = 0;
    while( parcours && out == 0 )
    {
       if ( !strcasecmp(parcours->usemodule,name) ) out = 1;
       else parcours = parcours->suiv;
    }
    if ( out == 0 )
    {
       newmodule =(listusemodule *)malloc(sizeof(listusemodule));
       strcpy(newmodule->usemodule,name);
       strcpy(newmodule->cursubroutine,subroutinename);  
       newmodule->suiv = listofmoduletmp;
       listofmoduletmp = newmodule;
    }
  }
}

/******************************************************************************/
/*                     Add_ModuleTo_Modulelist                                */
/******************************************************************************/
/* This subroutine is used to add a        variable to the list               */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ glob +--->+ glob +--->+ glob +--->+ glob +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void Add_ModuleTo_Modulelist(char *nom)
{
   listnom *newnom;

   newnom=(listnom *) malloc (sizeof (listnom));
   strcpy(newnom->nom,nom);
   newnom->suiv = modulelist;
   modulelist = newnom;   
}

/******************************************************************************/
/*                 OPTI_1_completelistvarpointtovar                           */
/******************************************************************************/
/* Firstpass 1                                                                */
/* We should complete the listvarpointtovar                                   */ 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_1_completelistvarpointtovar(char *namemodule,listcouple *couple)
{
   listvarpointtovar *pointtmp;
   
   if ( firstpass == 1 ) 
   {
      /* we should complete the Listofvarpointtovar                           */
      pointtmp=(listvarpointtovar *)malloc(sizeof(listvarpointtovar));
      strcpy(pointtmp->usemodule,namemodule);
      strcpy(pointtmp->cursubroutine,subroutinename);
      pointtmp->couple = couple;
      if ( Listofvarpointtovar )
      {
         pointtmp->suiv = Listofvarpointtovar;
         Listofvarpointtovar = pointtmp;
      }
      else
      {
         pointtmp->suiv = NULL;
         Listofvarpointtovar = pointtmp;    
      }
   }
}

/******************************************************************************/
/*                        Addincludetothelist                                 */
/******************************************************************************/
/* This subroutine is used to add a record to a list of struct                */
/*  listofincludebysubroutine                                                 */
/******************************************************************************/
/*                                                                            */
/*       subroutine sub ... USE mod1 ===> insert in list                      */
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ list +--->+ list +--->+ list +--->+ list +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*       list =  listofmodulebysubroutine                                     */
/*                                                                            */
/******************************************************************************/
void Addincludetothelist(char *name)
{
  listusemodule *newinclude;

  newinclude =(listusemodule *)malloc(sizeof(listusemodule));
  strcpy(newinclude->usemodule,name);
  strcpy(newinclude->cursubroutine,subroutinename);  
  newinclude->suiv = NULL;

  if ( !listofincludebysubroutine)
  {
     listofincludebysubroutine  = newinclude ;
  }
  else
  {
    newinclude->suiv = listofincludebysubroutine;
    listofincludebysubroutine = newinclude;
  }
}


/******************************************************************************/
/*                        WriteIncludeDeclaration                             */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void WriteIncludeDeclaration()
{
  listusemodule *newinclude;

  newinclude = listofincludebysubroutine;
  fprintf(fortranout,"\n");
  while ( newinclude )
  {
     if ( !strcasecmp(newinclude->cursubroutine,subroutinename) )
     {
        fprintf(fortranout,"      INCLUDE %s \n",newinclude->usemodule);
     }
     newinclude = newinclude ->suiv;  
  }
}
