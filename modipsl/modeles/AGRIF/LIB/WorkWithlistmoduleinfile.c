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
/*                   MOD_1_FillInlistmodule                                   */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void MOD_1_FillInlistmodule()
{
   listmodule *tmplist;
   
   
   if (firstpass == 1) 
   {   
      tmplist = (listmodule *)malloc(sizeof(listmodule));
      strcpy(tmplist->module,curmodulename);
      tmplist->InstanceShouldMade = 0;
      tmplist->Instance = 0;
      /*         */
      if ( !listmoduleinfile)
      {
         listmoduleinfile = tmplist;  
         tmplist->suiv = NULL;
      }
      else
      {
         tmplist->suiv = listmoduleinfile;
         listmoduleinfile = tmplist;      
      }
   }
}


/******************************************************************************/
/*                MOD_1_InstanceShouldMadeTo0InModule                         */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void MOD_1_InstanceShouldMadeTo0InModule()
{
   listmodule *tmplist;
   
   
   if (firstpass == 1 && listmoduleinfile ) 
   {
      tmplist=listmoduleinfile;
      /* we should find the module in the listmoduleinfile                    */
      while ( strcasecmp(tmplist->module,curmodulename) ) tmplist=tmplist->suiv;
      /* and turn the flag to 0                                               */
      tmplist->InstanceShouldMade = 0 ;
   }
}


/******************************************************************************/
/*                MOD_1_InstanceShouldMadeTo1InModule                         */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void MOD_1_InstanceShouldMadeTo1InModule()
{
   listmodule *tmplist;
   
   
   if (firstpass == 1 && listmoduleinfile ) 
   {
      tmplist=listmoduleinfile;
      /* we should find the module in the listmoduleinfile                    */
      while ( strcasecmp(tmplist->module,curmodulename) ) tmplist=tmplist->suiv;
      /* and turn the flag to 0                                               */
      tmplist->InstanceShouldMade = 1 ;
   }
}

/******************************************************************************/
/*                     MOD_1_InstanceTo1InModule                              */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void MOD_1_InstanceTo1InModule()
{
   listmodule *tmplist;
   
   
   if (firstpass == 1 && listmoduleinfile ) 
   {
      tmplist=listmoduleinfile;
      /* we should find the module in the listmoduleinfile                    */
      while ( strcasecmp(tmplist->module,curmodulename) ) tmplist=tmplist->suiv;
      /* and turn the flag to 0                                               */
      tmplist->Instance = 1 ;
   }
}

/******************************************************************************/
/*                     MOD_n_InstanceShouldMadeInModule                       */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int MOD_n_InstanceShouldMadeInModule()
{
   listmodule *tmplist;
   
   
   if ( listmoduleinfile ) 
   {
      tmplist=listmoduleinfile;
      /* we should find the module in the listmoduleinfile                    */
      while ( strcasecmp(tmplist->module,curmodulename) ) tmplist=tmplist->suiv;
      /* and turn the flag to 0                                               */
      return tmplist->InstanceShouldMade;
   }
   else
   {
      return 0;
   }
}

/******************************************************************************/
/*                          MOD_n_InstanceInModule                            */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int MOD_n_InstanceInModule()
{
   listmodule *tmplist;
   
   
   if ( listmoduleinfile ) 
   {
      tmplist=listmoduleinfile;
      /* we should find the module in the listmoduleinfile                    */
      while ( strcasecmp(tmplist->module,curmodulename) ) tmplist=tmplist->suiv;
      /* and turn the flag to 0                                               */
      return tmplist->Instance;
   }
   else
   {
      return 0;
   }
}
