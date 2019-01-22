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
/*                            AddListvartolistvar                             */
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
listvar * AddListvarToListvar(listvar *l,listvar *glob,int ValueFirstpass)
{
   listvar *newvar;

   if ( firstpass == ValueFirstpass )
   {
      if ( !glob) glob = l ;
      else
      {
         newvar=glob;
         while (newvar->suiv) newvar = newvar->suiv;
         newvar->suiv = l;
      }
   }
   return glob;
}

/******************************************************************************/
/*                       CreateAndFillin_Curvar                               */
/******************************************************************************/
/* This subroutine is used to create the record corresponding to the          */
/* list of declaration                                                        */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void CreateAndFillin_Curvar(char *type,char *tokname,
                            listdim *dims,variable *curvar)
{
   if (!strcasecmp(type,"character") && CharacterSizeGiven == 1 )    
                            strcpy(curvar->dimchar,CharacterSize);

  /* On donne la precision de la variable si elle a ete donnee                */
  curvar->c_star = 0;
  if ( c_star == 1 ) curvar->c_star = 1;
  /*                                                                          */
  if ( lengspecgiven == 1 ) strcpy(curvar->vallengspec,vallengspec);
  curvar->lengspecgiven=0;
  if ( lengspecgiven == 1 ) curvar->lengspecgiven=1;

  if ( PrecisionGiven == 1 )  strcpy(curvar->precision,NamePrecision);
  /* Si cette variable a ete declaree dans un module on met curvar->module=1  */
  if ( inmoduledeclare == 1 || SaveDeclare == 1)
  {
      curvar->module = 1;
      /* Puis on donne le nom du module dans curvar->modulename               */
      strcpy(curvar->modulename,curmodulename);
   }
   else if (insubroutinedeclare == 1 )
   /* we give the name of the subroutine to the modulename                    */
   {
      strcpy(curvar->modulename,subroutinename);
   }
   /* Si cette variable a ete initialisee                                     */
   if (InitialValueGiven == 1 ) strcpy(curvar->initialvalue,InitValue); 
   /* Si cette variable est declaree en save                                  */
   if (SaveDeclare == 1 ) curvar->save = 1;
   /* Si cette variable est allocatable                                       */
   if (Allocatabledeclare == 1 ) curvar->allocatable=1;
   /* if INTENT spec has been given                                           */
   if ( IntentDeclare == 1 ) strcpy(curvar->IntentSpec,IntentSpec);
}


/******************************************************************************/
/*                        duplicatelistvar                                    */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
listvar * duplicatelistvar(listvar * orig)
{
   listvar *newlist;
   listvar *parcours;
   listvar *tmplistvar;
   listvar *tmplistvarprec;
   listdim *tmplistdim;
   variable *tmpvar;

   tmplistvarprec = (listvar *)NULL;
   newlist = (listvar *)NULL;
   parcours = orig;
   while ( parcours )
   {
      tmplistvar = (listvar *)malloc(sizeof(listvar));
      tmpvar = (variable *)malloc(sizeof(variable));
      /*                                                                      */
      strcpy(tmpvar->typevar,parcours->var->typevar);      
      strcpy(tmpvar->nomvar,parcours->var->nomvar);      
      strcpy(tmpvar->oldname,parcours->var->oldname);      
      strcpy(tmpvar->dimchar,parcours->var->dimchar);      
      if ( parcours->var->dimension )
      {
         tmplistdim = (listdim *)malloc(sizeof(listdim));
         tmplistdim = parcours->var->dimension;
         tmpvar->dimension = tmplistdim;
      }
      tmpvar->nbdim=parcours->var->nbdim;
      tmpvar->common=parcours->var->common;
      tmpvar->positioninblock=parcours->var->positioninblock;
      tmpvar->module=parcours->var->module;
      tmpvar->save=parcours->var->save;
      tmpvar->VariableIsParameter=parcours->var->VariableIsParameter;
      strcpy(tmpvar->modulename,parcours->var->modulename);      
      strcpy(tmpvar->commonname,parcours->var->commonname);      
      strcpy(tmpvar->vallengspec,parcours->var->vallengspec);
      strcpy(tmpvar->nameinttypename,parcours->var->nameinttypename);
      tmpvar->lengspecgiven=parcours->var->lengspecgiven;
      tmpvar->pointedvar=parcours->var->pointedvar;
      strcpy(tmpvar->commoninfile,parcours->var->commoninfile);      
      strcpy(tmpvar->subroutinename,parcours->var->subroutinename);      
      tmpvar->dimensiongiven=parcours->var->dimensiongiven;
      tmpvar->c_star=parcours->var->c_star;
      tmpvar->typegiven=parcours->var->typegiven;
      tmpvar->isparameter=parcours->var->isparameter;
      strcpy(tmpvar->precision,parcours->var->precision);
      strcpy(tmpvar->initialvalue,parcours->var->initialvalue);
      tmpvar->pointerdeclare=parcours->var->pointerdeclare;
      tmpvar->optionaldeclare=parcours->var->optionaldeclare;
      tmpvar->allocatable=parcours->var->allocatable;
      strcpy(tmpvar->IntentSpec,parcours->var->IntentSpec);
      tmpvar->dimsempty=parcours->var->dimsempty;
      strcpy(tmpvar->readedlistdimension,parcours->var->readedlistdimension);
      /*                                                                      */
      tmplistvar->var = tmpvar;
      tmplistvar->suiv = NULL;
      /*                                                                      */
      if ( !newlist )
      {
         newlist = tmplistvar;
         tmplistvarprec = newlist;
      }
      else
      {
         tmplistvarprec->suiv = tmplistvar;
         tmplistvarprec = tmplistvar;
      }
      /*                                                                      */
      parcours = parcours->suiv;
   }
   return newlist;
}

/******************************************************************************/
/*                           insertdim                                        */
/******************************************************************************/
/* This subroutine is used to insert a record in a list of                    */
/* struct : listdim                                                           */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ lin  +--->+ lin  +--->+ lin  +--->+  lin +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/******************************************************************************/
listdim * insertdim(listdim *lin,typedim nom)
{
   listdim *newdim ;

   newdim=(listdim *) malloc (sizeof (listdim));
   newdim->dim=nom;
   newdim->suiv=lin;
   
   return newdim;
}

/******************************************************************************/
/*                           reverse                                          */
/******************************************************************************/
/* This subroutine is used to reverse a list                                  */
/******************************************************************************/
/*        _______     _______                 _______     _______             */
/*       +      +    +      +                +      +    +      +             */
/*       +  A   +--->+   B  +--------------->+  B   +--->+   A  +             */
/*       +______+    +______+                +______+    +______+             */
/*                                                                            */
/******************************************************************************/
listdim *reverse(listdim *lin)
{
   listdim *newdim1;
   listdim *newdim2;
   listdim *lout;

   lout=(listdim *) NULL;

   newdim1=lin;
   while (newdim1)
   {
      newdim2=(listdim *) malloc(sizeof(listdim));
      newdim2->dim=newdim1->dim;
      newdim2->suiv=lout;
      lout=newdim2;
      newdim1=newdim1->suiv;
   }
   return lout;
}

/******************************************************************************/
/*                            change_dim_char                                 */
/******************************************************************************/
/* This subroutine is used to change the dimension in the list lin            */
/******************************************************************************/
/*        _______     _______                 _______     _______             */
/*       +  l   +    +  l   +                +  l   +    +   l  +             */
/*       + old  +--->+ old  +--------------->+ lin  +--->+  lin +             */
/*       +______+    +______+                +______+    +______+             */
/*                                                                            */
/******************************************************************************/
void change_dim_char(listdim *lin,listvar * l)
{
   listvar *parcours_var;
   variable *v;
  
   
   parcours_var=l;
   while(parcours_var)
   { 
      v=parcours_var->var;
      strcpy(v->dimchar,(lin->dim).last);
      parcours_var=parcours_var->suiv;
   }
}


/******************************************************************************/
/*                                num_dims                                    */
/******************************************************************************/
/* This subroutine is used to know the dimension of a table                   */
/******************************************************************************/
/*                                                                            */
/*             Dimension(jpi,jpj,jpk) ----------> num_dims = 3                */
/*                                                                            */
/******************************************************************************/
int num_dims(listdim *d)
{
   listdim *parcours;
   int compteur = 0;

   parcours = d;
   while(parcours)
   {
     compteur++;
     parcours=parcours->suiv;
   }
   return compteur;  
}


/******************************************************************************/
/*                          CREATEVAR                                         */
/******************************************************************************/
/* This subroutine is used to create and initialized a record of the          */
/*      struct : variable                                                     */
/******************************************************************************/
variable * createvar(char *nom,listdim *d)
{
  variable *var;
  listdim *dims;
  char ligne[LONGNOM];
  char listdimension[LONGNOM];

   var=(variable *) malloc(sizeof(variable));
   strcpy(var->nomvar,nom);
   /* Definition of the number of this variable in the table tabvars          */
   var->indicetabvars = 0;
   if ( firstpass == 1 && ( aftercontainsdeclare == 0 || 
                            SaveDeclare == 1          ||
                            fortran77 == 1 ) 
      )
   {
      indicemaxtabvars = indicemaxtabvars + 1;
      var->indicetabvars = indicemaxtabvars;
   }
   /*                                                                         */
   var->pointerdeclare=0;
   var->dimsempty=0;
   var->optionaldeclare=0;
   var->dimensiongiven=0;
   var->isparameter=0;
   var->positioninblock=0;
   var->VariableIsParameter = 0;
   var->PublicDeclare = 0;
   var->PrivateDeclare = 0;
   var->ExternalDeclare = 0;
   var->common=0;
   var->allocatable=0;
   var->module=0; 
   var->typegiven=0;
   var->save=0;
   /*                                                                         */
   strcpy(var->nameinttypename,"");
   strcpy(listdimension,"");
   strcpy(var->modulename,"");
   strcpy(var->commonname,"");
   strcpy(var->commoninfile,mainfile);
   strcpy(var->subroutinename,subroutinename);
   strcpy(var->dimchar,"");
   strcpy(var->oldname,"");
   strcpy(var->precision,""); 
   strcpy(var->initialvalue,""); 
   strcpy(var->IntentSpec,""); 
   /*                                                                         */
   if ( inttypename         == 1 ) strcpy(var->nameinttypename,nameinttypename);
   if ( optionaldeclare     == 1 ) var->optionaldeclare = 1;
   if ( pointerdeclare      == 1 ) var->pointerdeclare = 1;
   if ( VariableIsParameter == 1 ) var->isparameter = 1;
   if ( VariableIsParameter == 1 ) var->VariableIsParameter = 1 ;
   if ( PublicDeclare       == 1 ) var->PublicDeclare = 1 ;
   if ( PrivateDeclare      == 1 ) var->PrivateDeclare = 1;
   if ( ExternalDeclare     == 1 ) var->ExternalDeclare = 1; 
   /*                                                                         */
   var->dimension=d;
   /* Creation of the string for the dimension of this variable               */
   dimsempty = 1;
   if ( d )
   {
      var->dimensiongiven=1;
      dims = d;
      while (dims)
      {
         if ( strcasecmp(dims->dim.first,"") || strcasecmp(dims->dim.last,""))
                                                                  dimsempty = 0;
         sprintf(ligne,"%s:%s",dims->dim.first,dims->dim.last);
         strcat(listdimension,ligne);
         if ( dims->suiv )
         {
            strcat(listdimension,",");	     
         }
         dims = dims->suiv;
      }
      if ( dimsempty == 1 ) var->dimsempty=1;
   }
   strcpy(var->readedlistdimension,listdimension);
   /*                                                                         */
   var->nbdim=num_dims(d);
   /*                                                                         */
   return var;
}

/******************************************************************************/
/*                            INSERTVAR                                       */
/******************************************************************************/
/* This subroutine is used to insert a record in a list of the                */
/*      struct : listvar                                                      */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       +  lin +--->+  lin +--->+ lin  +--->+ lin  +--->+ NEW  +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
listvar * insertvar(listvar *lin,variable *v)
{
   listvar *newvar ;
   listvar *tmpvar ;

   newvar=(listvar *) malloc (sizeof (listvar));
   newvar->var=v;
   newvar->suiv = NULL;
   if (!lin)
   {
      newvar->suiv=NULL;
      lin = newvar;
   }
   else
   {
      tmpvar = lin ;
      while (tmpvar->suiv)
      {
         tmpvar = tmpvar ->suiv ;
      }
      tmpvar -> suiv = newvar;   
   }
   return lin;
}

/******************************************************************************/
/*                             SETTYPE                                        */
/******************************************************************************/
/* This subroutine is used to give the same variable type at each             */
/*      record of the list of the struct : listvar                            */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       + REAL +    + REAL +    + REAL +    + REAL +    + REAL +             */
/*       +  lin +--->+  lin +--->+ lin  +--->+ lin  +--->+ lin  +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
listvar *settype(char *nom,listvar *lin)
{
   listvar *newvar;
   variable *v;

   newvar=lin;
   while (newvar)
   {
      v=newvar->var;
      strcpy(v->typevar,nom);
      v->typegiven=1;
      newvar=newvar->suiv;
   }
   newvar=lin;
   return newvar ;
}

