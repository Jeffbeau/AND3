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
/*                     CompareNewparsingandoldone                             */
/******************************************************************************/
/* this subroutine is used to compare the old treatement with                 */
/* the new one                                                                */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + glob +--->+ glob +--->+ glob +--->+ glob +--->+ glob +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*           =         not=         =                                         */
/*        _______     _______     _______                                     */
/*       +      +    +      +    +      +                                     */
/*       + tmp  +--->+ tmp  +--->+ tmp  +                                     */
/*       +______+    +______+    +______+                                     */
/*                                                                            */
/******************************************************************************/
void CompareNewparsingandoldone()
{
   listvar *NewTreated;
   listvar *OldTreated;
   listvar *OldTreatedPrec;
   int Out;
   listindice *newindice;

   OldTreatedPrec = (listvar *)NULL;
   NewTreated = globliste;
   while ( NewTreated )
   {
         /* we are going to compare the two treatement                        */
         /* each time that we meet in the oldlist the same record we          */
         /* remove it from the tmplocallist                                   */
         OldTreated = tmplocallist;
         Out = 0;
         while ( OldTreated && Out == 0 )
         { 
            if ( !strcasecmp(NewTreated->var->nomvar, OldTreated->var->nomvar))
            {
               /* We should keep the same indice for tabvars table than       */
               /* the old one                                                 */
               NewTreated->var->indicetabvars = OldTreated->var->indicetabvars;
               /* we remove it from the tmplocallist                          */
               if ( OldTreated == tmplocallist ) 
               {
                  tmplocallist = tmplocallist -> suiv;
               }
               else
               {
                  OldTreatedPrec->suiv = OldTreated -> suiv;
               }  
               /* We go out of the loop because we find two variables         */
               /* with the same name in the same file                         */
               Out = 1;
            }
            else
            {
               OldTreatedPrec = OldTreated;
               OldTreated = OldTreated -> suiv;
            }
      }
      
      if ( !OldTreated  && Out == 0) 
      {
         /* if this tmplocallist has been readed without finding the          */
         /* variable of the globliste                                         */
         /* it means that this variable has been added                        */
         /* in this case we choose a new tabvars indice                       */
	 /* for this variable                                                 */
         if ( Listofavailableindices )
         {
            NewTreated->var->indicetabvars = Listofavailableindices -> indice;
            Listofavailableindices = Listofavailableindices ->suiv;
         }
         else
         {
            indicemaxtabvars = indicemaxtabvars + 1;
            NewTreated->var->indicetabvars = indicemaxtabvars;
         }
      }
      /*  On passe a l'enregistrement suivant */
      NewTreated = NewTreated -> suiv;
   }
   if ( tmplocallist )
   {
      /* if the tmplocallist is not empty it means that some variables        */
      /* has been removed in the new version of the filetoparse               */
      /* in this case we should record the indice of the tabvars              */
      /* to know that this field is empty                                     */
      while (tmplocallist)
      {    
         if ( tmplocallist -> var -> indicetabvars != 0 )
         {
            newindice=(listindice *) malloc (sizeof (listindice));
            newindice -> indice = tmplocallist -> var -> indicetabvars;
            newindice -> suiv = Listofavailableindices;
            Listofavailableindices = newindice;
         }
         tmplocallist = tmplocallist -> suiv;
      }
   }
}


/******************************************************************************/
/*                          decl_1_ajoutevar                                  */
/******************************************************************************/
/* Firstpass 1                                                                */
/* We should add this declaration to the globliste                            */
/******************************************************************************/
void decl_1_ajoutevar(listvar *listtoadd)
{
   if ( firstpass == 1 && VariableIsParameter == 0 )
   if ( aftercontainsdeclare == 0 || fortran77 == 1 )
   {
      globliste = AddListvarToListvar(listtoadd,globliste,1);
   }
}

/******************************************************************************/
/*                          decl_1_ajoutevarsave                              */
/******************************************************************************/
/* Firstpass 1                                                                */
/* We should add this declaration to the globliste. case SAVE                 */
/******************************************************************************/
void decl_1_ajoutevarsave(listvar *listtoadd)
{
   if ( VariableIsParameter == 0 && SaveDeclare == 1 && firstpass == 1 )
   {
      globliste = AddListvarToListvar(listtoadd,globliste,1);
   }      
}

/******************************************************************************/
/*                       UpdateIndiceTabvarsofGlobliste                      */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void UpdateIndiceTabvarsofGlobliste()
{

   int indicetmp;
   listvar *NewTreated;
   listvar *OldTreated;
   listvar *OldTreatedPrec;
   int Out;
   FILE *dependglobaloutput;
   
   
   if ( Did_filetoparse_treated == 1 )
   {

   Recordtmplocallist(NameTamponfile);

   if ((dependglobaloutput=fopen(".dependglobal","r"))!=NULL) 
   {
      fscanf(dependglobaloutput,"%d\n",&indicemaxtabvars);
      fclose(dependglobaloutput);
      oldindicemaxtabvars = indicemaxtabvars;
   }
   /* Read the list of available indice                                       */
   Readthedependavailablefile();
   /*                                                                         */
   indicetmp = indicemaxtabvars;
   OldTreatedPrec = (listvar *)NULL;
   NewTreated = globliste;
   while ( NewTreated )
   {
         /* we are going to compare the two treatement                        */
         /* each time that we meet in the oldlist the same record we          */
         /* remove it from the tmplocallist                                   */
         OldTreated = tmplocallist;
         Out = 0;
         while ( OldTreated && Out == 0 )
         { 
            if ( !strcasecmp(NewTreated->var->nomvar, OldTreated->var->nomvar))
            {
               /* We should keep the same indice for tabvars table than       */
               /* the old one                                                 */
               NewTreated->var->indicetabvars = OldTreated->var->indicetabvars;
               /* we remove it from the tmplocallist                          */
               if ( OldTreated == tmplocallist ) 
               {
                  tmplocallist = tmplocallist -> suiv;
               }
               else
               {
                  OldTreatedPrec->suiv = OldTreated -> suiv;
               }  
               /* We go out of the loop because we find two variables         */
               /* with the same name in the same file                         */
               Out = 1;
            }
            else
            {
               OldTreatedPrec = OldTreated;
               OldTreated = OldTreated -> suiv;
            }
      }
      
      if ( !OldTreated  && Out == 0) 
      {
         /* if this tmplocallist has been readed without finding the          */
         /* variable of the globliste                                         */
         /* it means that this variable has been added                        */
         /* in this case we choose a new tabvars indice                       */
	 /* for this variable                                                 */
         if ( Listofavailableindices )
         {
            NewTreated->var->indicetabvars = Listofavailableindices -> indice;
            Listofavailableindices = Listofavailableindices ->suiv;
         }
         else
         {
            indicetmp = indicetmp + 1;
            NewTreated->var->indicetabvars = indicetmp;
         }
      }
      /*  On passe a l'enregistrement suivant */
      NewTreated = NewTreated -> suiv;
   }
   tmplocallist = NULL;
   
   } /* end of Did_filetoparse_treated == 1                                   */
}

/******************************************************************************/
/*                  UpdateIndiceTabvarsofGloblisteFromCommon                  */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void UpdateIndiceTabvarsofGloblisteFromCommon()
{

   listnom *parcours;
   listvar *parcours2;
   listvar *listtmp;
   listvar *parcoursglob;
   listvar *parcoursglob2;
   int find;
   int NumberofNewVar;
   int NumberofVar;
   
   listtmp = ( listvar *)NULL ;
   NumberofNewVar=0;
   NumberofVar=0;
   parcours = modulelist;
   while( parcours )
   {
      listtmp = Recordglobalvarofusefile( parcours->nom  , listtmp);
      parcours=parcours->suiv;
   }
   /*                                                                         */
   parcoursglob = globliste;
   /* if this common has been ever read, we should update the tabvars         */
   /*    indices                                                              */
   while ( parcoursglob )
   {
      NumberofVar = NumberofVar +1 ;
      parcours2 = listtmp;
      find = 0 ;
      while ( parcours2 && find == 0 )
      {
         if ( !strcasecmp(parcoursglob->var->nomvar,parcours2->var->nomvar) &&
              !strcasecmp(parcoursglob->var->commonname,
                                                parcours2->var->commonname)
            )
         {
            parcoursglob->var->indicetabvars = parcours2->var->indicetabvars;
            strcpy(parcoursglob->var->commoninfile,
                                               parcours2->var->commoninfile);
               find = 1;
         }
         parcours2 = parcours2->suiv;
      }
      /* if we did not find it, it means that it is a new variable            */
      /*    we should see if this variable has not been defined twice later   */
      if ( find == 0 && (parcoursglob->var->indicetabvars > 
           (oldindicemaxtabvars + NumberofNewVar)) )
      {
         NumberofNewVar = NumberofNewVar +1 ;
         parcoursglob->var->indicetabvars = oldindicemaxtabvars 
                                                               + NumberofNewVar;
         parcoursglob2 = parcoursglob;
         while ( parcoursglob2 )
         {
            if ( !strcasecmp(parcoursglob->var->nomvar,
                                                parcoursglob2->var->nomvar) &&
                 !strcasecmp(parcoursglob->var->commonname,
                                            parcoursglob2->var->commonname)
               )
            {
               parcoursglob2->var->indicetabvars = oldindicemaxtabvars 
                                                               + NumberofNewVar;
            }
            parcoursglob2 = parcoursglob2->suiv;
         }
      }
      /*                                                                      */
      parcoursglob  = parcoursglob  ->suiv;
   }
   indicemaxtabvars = oldindicemaxtabvars + NumberofNewVar;
}


/******************************************************************************/
/*                   COM_1_UpdateGloblisteWithcommonlist                      */
/******************************************************************************/
/*  This subroutines is used to add the variable defined in common in the     */
/*     commonlist                                                             */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void COM_1_UpdateGloblisteWithcommonlist()
{
   listvarcommon *parcours;
   listvar *parcoursindic;
   listvar *parcoursvar;
   listvar *parcoursvarprec;
   int out;
   
   parcoursvar = globliste;
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
            /* we should update the globliste                                 */
            strcpy(parcoursvar->var->commonname,parcours->commonname);
            parcoursvar->var->positioninblock = parcours->positioninblock;
            parcoursvar->var->common = 1;
         }
         else
         {
            parcours = parcours->suiv;
         }
      }
      if ( out == 0 )
      {
         /* We should update the tabvarsindic of the following variable       */
         /*    present in the globliste                                       */
         parcoursindic = parcoursvar;
         indicemaxtabvars=indicemaxtabvars-1;
         while(parcoursindic)
         {
            parcoursindic->var->indicetabvars =
                                            parcoursindic->var->indicetabvars-1;
            parcoursindic = parcoursindic->suiv;
         }
         /* we did not find it                                                */
         /* we should remove the variable from the globliste                  */
         if ( parcoursvar == globliste )
         {
            globliste = globliste->suiv;
            parcoursvar = globliste;
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
