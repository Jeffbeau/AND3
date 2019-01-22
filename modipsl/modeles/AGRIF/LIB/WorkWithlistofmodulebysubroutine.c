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
/*                    RecordUseModulesVariables                               */
/******************************************************************************/
/******************************************************************************/
void RecordUseModulesVariables()
{
  listusemodule *tmplistmodule;

  /* we should record all variables defined in modules used in this file      */
  if ( listofmodulebysubroutine )
  {
     globalvarofusefile = (listvar *)NULL;
     tmpparameterlocallist = (listparameter *)NULL;
     listofvarofusemodulecreated = 1;
     tmplistmodule = listofmodulebysubroutine;
     while ( tmplistmodule )
     {
        if ( tmplistmodule->firstuse == 1 )
        {
           /* check if the file .depend<usemodule> exist                      */
           globalvarofusefile = Recordglobalvarofusefile
                                  (tmplistmodule->usemodule,globalvarofusefile);

           tmpparameterlocallist = ReaddependParameterList
                               (tmplistmodule->usemodule,tmpparameterlocallist);

        }

        tmplistmodule = tmplistmodule->suiv;
     }
  }
}

/******************************************************************************/
/*                RecordUseModulesUseModulesVariables                         */
/******************************************************************************/
/******************************************************************************/
void  RecordUseModulesUseModulesVariables()
{
  listusemodule *tmplistmodule;

  /* we should record all variables defined in modules used in this file      */
  if ( listofmodulebysubroutine )
  {
     /* and we should read the .depend of the module used by the module used  */
     tmplistmodule = listofmodulebysubroutine;
     while ( tmplistmodule )
     {
        Readthedependlistofmoduleused(tmplistmodule->usemodule);
        while( tmpuselocallist )
        {
           Addmoduletothelisttmp(tmpuselocallist->usemodule);
           tmpuselocallist = tmpuselocallist->suiv;
        }
        tmplistmodule = tmplistmodule->suiv;
     }
           
     globalvarofusefile2 = (listvar *)NULL;
     tmpparameterlocallist2 = (listparameter *)NULL;
     tmplistmodule = listofmoduletmp;
     while ( tmplistmodule )
     {
        /* check if the file .depend<usemodule> exist                         */
        globalvarofusefile2 = Recordglobalvarofusefile
                                 (tmplistmodule->usemodule,globalvarofusefile2);

        tmpparameterlocallist2 = ReaddependParameterList
                              (tmplistmodule->usemodule,tmpparameterlocallist2);
        
        tmplistmodule = tmplistmodule->suiv;
     }
  }
}


/******************************************************************************/
/*                        Addmoduletothelist                                  */
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
/*       list =  listofmodulebysubroutine                                     */
/*                                                                            */
/******************************************************************************/
void Addmoduletothelist(char *name)
{
  listusemodule *newmodule;
  listusemodule *parcours;
  int out;

  newmodule =(listusemodule *)malloc(sizeof(listusemodule));
  strcpy(newmodule->usemodule,name);
  strcpy(newmodule->charusemodule,charusemodule);
  strcpy(newmodule->cursubroutine,subroutinename);  
  newmodule->firstuse = 1 ;  
  newmodule->suiv = NULL;

  if ( !listofmodulebysubroutine)
  {
      listofmodulebysubroutine = newmodule ;
  }
  else
  {
    parcours = listofmodulebysubroutine;
    while ( parcours && newmodule->firstuse == 1 )
    {
       if ( !strcasecmp(name,parcours->usemodule) ) 
       {
          newmodule->firstuse = 0 ;
       }
       parcours=parcours->suiv;
    }
    /* we can not add the same module twice for the same subroutine           */
    parcours = listofmodulebysubroutine;
    out = 0 ;
    while ( parcours && out == 0 )
    {
       if ( !strcasecmp(name,parcours->usemodule) &&
            !strcasecmp(subroutinename,parcours->cursubroutine)
           )
       {
          out = 1 ;
          free(newmodule);
       }
       else parcours=parcours->suiv;
    }
    if ( out == 0 )
    {
       newmodule->suiv = listofmodulebysubroutine;
       listofmodulebysubroutine = newmodule;
    }
  }
}


/******************************************************************************/
/*                        WriteUsemoduleDeclaration                           */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void WriteUsemoduleDeclaration()
{
  listusemodule *newmodule;

  newmodule = listofmodulebysubroutine;
  fprintf(fortranout,"\n");
  while ( newmodule )
  {
     if ( !strcasecmp(newmodule->cursubroutine,subroutinename) )
     {
        if ( strcasecmp(newmodule->charusemodule,"Agrif_Util") ||
            adduseagrifutil != 1 ) fprintf(fortranout,"      USE %s \n"
                                                     ,newmodule->charusemodule);
     }
        newmodule = newmodule ->suiv;  
  }
}
