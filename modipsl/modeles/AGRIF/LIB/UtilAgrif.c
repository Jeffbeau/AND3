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
/*                      AGRIF_n_Vartonumber                                   */
/******************************************************************************/
/* This subroutine is used to know if Agrif_ is locate in the char            */
/* tokname                                                                    */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int AGRIF_n_Vartonumber(char *tokname)
{
   int agrifintheword;
   
   agrifintheword = 0;
        if ( !strcasecmp(tokname,"Agrif_parent")         ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_set_type")       ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_set_raf")        ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_set_bc")         ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_set_bcinterp")   ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_bc_variable")    ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_set_parent")     ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_interp_variable")) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_init_variable")  ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_update_variable")) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_Set_interp")     ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_Set_Update")     ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_Set_UpdateType") ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"Agrif_Set_restore")    ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"agrif_init_grids")     ) agrifintheword = 1;
   else if ( !strcasecmp(tokname,"agrif_step")           ) agrifintheword = 1;

   return agrifintheword;
}

/******************************************************************************/
/*                      AGRIF_n_Agrif_in_Tok_NAME                             */
/******************************************************************************/
/* This subroutine is used to know if Agrif_ is locate in the char            */
/* tokname                                                                    */
/******************************************************************************/
/*                                                                            */
/*                 Agrif_name --------------> Agrif_in_Tok_NAME = 1           */
/*                       name --------------> Agrif_in_Tok_NAME = 0           */
/*                                                                            */
/******************************************************************************/
int AGRIF_n_Agrif_in_Tok_NAME(char *tokname)
{
   int agrifintheword;
   
   if ( strncasecmp(tokname,"Agrif_",6) == 0 )  agrifintheword = 1;
   else agrifintheword = 0;

   return agrifintheword;
}


/******************************************************************************/
/*          AGRIF_1_completeListofvariableinagriffunction                     */
/******************************************************************************/
/* Firstpass 1                                                                */
/* We should complete the Listofvariableinagriffunction                       */
/******************************************************************************/
/*                                                                            */
/*               Agrif_Parent(variable) ====>         variable                */
/*                                                        \                   */
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + list +--->+ list +--->+ list +--->+ list +--->+ NEW  +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*       list =  Listofvariableinagriffunction                                */
/*                                                                            */
/******************************************************************************/
void AGRIF_1_completeListofvariableinagriffunction(char *ident)
{
   listnom *listnomtmp;
   listnom *parcours;

   if ( firstpass == 1 ) 
   {
      if ( Listofvariableinagriffunction )
      {
         parcours = Listofvariableinagriffunction;
	 while ( parcours && strcasecmp(parcours->nom,ident) )
	 {
	    parcours = parcours->suiv;
	 }
	 if ( !parcours )
	 {
            listnomtmp=(listnom *)malloc(sizeof(listnom));
            strcpy(listnomtmp->nom,ident);
	    listnomtmp->suiv = NULL;	
	    listnomtmp->suiv = Listofvariableinagriffunction;
            Listofvariableinagriffunction = listnomtmp;
	 }
      }
      else
      {
         listnomtmp=(listnom *)malloc(sizeof(listnom));
         strcpy(listnomtmp->nom,ident);
	 listnomtmp->suiv = NULL;
	 Listofvariableinagriffunction = listnomtmp;
      }
   }
}

/******************************************************************************/
/*                     AGRIF_0_ModifyTheVariableName                          */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/*               Agrif_<toto>(variable) ====>     Agrif_<toto>(variable)      */
/*                                                                            */
/******************************************************************************/
void AGRIF_0_ModifyTheVariableName(char *ident)
{
   listvar *newvar;
   int out;

   if ( firstpass == 0 ) 
   {
      /* looking for the ident in the listvarindoloop                         */
      if ( inagrifcallargument == 1 ) 
      {
         if ( fortran77 == 0 ) newvar = globalvarofusefile;
         else newvar = globliste;
      }
      else newvar = globliste;
      out=0;
      while ( newvar && out == 0 ) 
      {
         if ( !strcasecmp(newvar->var->nomvar,ident) ) out = 1;
         else newvar=newvar->suiv;      
      }
      if ( out == 1 ) 
      {
         /* remove the variable                                               */
         RemoveWordCUR(fortranout,(long)(-strlen(ident)),
                               strlen(ident));
         fseek(fortranout,(long)(-strlen(ident)),SEEK_CUR);
         /* then write the new name                                           */
         if ( inagrifcallargument == 1 && agrif_parentcall == 0 )
            fprintf(fortranout,"%d",newvar->var->indicetabvars);
         else
         {
            fprintf(fortranout,"%s",vargridcurgridtabvars(newvar->var,0));
            colnum = strlen(vargridcurgridtabvars(newvar->var,0));
            if ( colnum >= 25 ) 
            {
              if ( fortran77 == 0 )
              {
                 fprintf(fortranout," & \n      ");
              }
              else
              {
                 fprintf(fortranout," \n     & ");                 
              }
            }
         }
      }
      else
      {
         /* we should look in the globalvarofusefile                          */
         if ( inagrifcallargument != 1 )
         {
            newvar = globalvarofusefile;
            while ( newvar && out == 0 ) 
            {
               if ( !strcasecmp(newvar->var->nomvar,ident) ) out = 1;
               else newvar=newvar->suiv;      
            }
            if ( out == 1 ) 
            {
               /* remove the variable                                         */
               RemoveWordCUR(fortranout,(long)(-strlen(ident)),
                                     strlen(ident));
               fseek(fortranout,(long)(-strlen(ident)),SEEK_CUR);
               /* then write the new name                                     */
               if ( colnum >= 25 ) 
               {
                 if ( fortran77 == 0 )
                 {
                    fprintf(fortranout," & \n      ");
                 }
                 else
                 {
                    fprintf(fortranout," \n     & ");                 
                 }
               }
               fprintf(fortranout,"%s",vargridcurgridtabvars(newvar->var,0));
               colnum = strlen(vargridcurgridtabvars(newvar->var,0));
            }
         }
      }
   }
}


/******************************************************************************/
/*             AGRIF_n_AddsubroutineTolistsubwhereagrifused                   */
/******************************************************************************/
/* This subroutine is used to add a record to                                 */
/* listofsubroutinewhereagrifisused                                           */
/******************************************************************************/
/*                                                                            */
/*       subroutine sub ... Agrif_<something>                                 */
/*                                                                            */
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + list +--->+ list +--->+ list +--->+ list +--->+ sub  +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*       list = listofsubroutinewhereagrifisused                              */
/*                                                                            */
/******************************************************************************/
void  AGRIF_n_AddsubroutineTolistsubwhereagrifused()
{
  listnom *listnomtmp;
  listnom *parcours;

  if ( !listofsubroutinewhereagrifisused )
  {
     listnomtmp=(listnom *)malloc(sizeof(listnom));
     strcpy(listnomtmp->nom,subroutinename);
     listnomtmp->suiv = NULL;	
     listofsubroutinewhereagrifisused  =  listnomtmp;
  }
  else
  {
    parcours = listofsubroutinewhereagrifisused;
    while ( parcours && strcasecmp(parcours->nom,subroutinename) )
    {
       parcours = parcours->suiv;
    }
    if ( !parcours )
    {
       listnomtmp=(listnom *)malloc(sizeof(listnom));
       strcpy(listnomtmp->nom,subroutinename);
       listnomtmp->suiv = listofsubroutinewhereagrifisused;	
       listofsubroutinewhereagrifisused  =  listnomtmp;       
    }
  }
}

/******************************************************************************/
/*                          AGRIF_n_AddUseAgrifUtil                           */
/******************************************************************************/
/* Add use Agrif_Util at the beginning of the subroutine definition           */
/* if it is necessary                                                         */
/******************************************************************************/
/*                                                                            */
/*       subroutine sub            |  subroutine sub                          */
/*                                 |  USE Agrif_Util                          */
/*       implicit none             |  implicit none                           */
/*       ...                       |  ...                                     */
/*       ... Agrif_<something>     |  ... Agrif_<something>                   */
/*       ...                       |  ...                                     */
/*       end                       |  end                                     */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void  AGRIF_n_AddUseAgrifUtil()
{
  listnom *parcours;
  listusemodule *newmodule;
  int out;

  parcours = listofsubroutinewhereagrifisused;
  while ( parcours && strcasecmp(parcours->nom,subroutinename) )  
  {
     parcours = parcours -> suiv;
  }
  if ( parcours )
  {
     /* we should add the use agrif_util if it is necessary                   */
     newmodule = listofmodulebysubroutine;
     out=0;
     while( newmodule && out == 0)
     {
        if ( !strcasecmp(newmodule->cursubroutine,subroutinename) ||
             !strcasecmp(newmodule->cursubroutine," ")  )
        {
           if ( !strcasecmp(newmodule->charusemodule,"Agrif_Util") ) out = 1 ;
        }
        newmodule = newmodule ->suiv;
     }
     
     if ( out == 0 && inmodulemeet == 0 ) 
     {
        fprintf(fortranout,"\n      USE Agrif_Util \n");
        adduseagrifutil = 1 ;
     }
  }
}


/******************************************************************************/
/*                         AGRIF_0_AgrifParentNotify                          */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/*               Agrif_<toto>(variable) ====>     Agrif_<toto>(variable)      */
/*                                                                            */
/******************************************************************************/
void AGRIF_0_NotifyAgrifFunction(char *ident)
{
   if ( firstpass == 0 ) 
   {
      if ( !strcasecmp(ident,"Agrif_parent") ) 
      {
         InAgrifParentDef = 1;
         pos_curagrifparent = setposcur()-12;
      }
      else if ( !strcasecmp(ident,"Agrif_Get_Coarse_grid") ) 
      {
         InAgrifParentDef = 2;
         pos_curagrifparent = setposcur()-21;
      }
      else if ( !strcasecmp(ident,"Agrif_Rhox") ) 
      {
         InAgrifParentDef = 3;
         pos_curagrifparent = setposcur()-10;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_Rhox") ) 
      {
         InAgrifParentDef = 4;
         pos_curagrifparent = setposcur()-17;
      }
      else if ( !strcasecmp(ident,"Agrif_IRhox") ) 
      {
         InAgrifParentDef = 5;
         pos_curagrifparent = setposcur()-11;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_IRhox") ) 
      {
         InAgrifParentDef = 6;
         pos_curagrifparent = setposcur()-18;
      }
      else if ( !strcasecmp(ident,"Agrif_Rhoy") ) 
      {
         InAgrifParentDef = 7;
         pos_curagrifparent = setposcur()-10;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_Rhoy") ) 
      {
         InAgrifParentDef = 8;
         pos_curagrifparent = setposcur()-17;
      }
      else if ( !strcasecmp(ident,"Agrif_IRhoy") ) 
      {
         InAgrifParentDef = 9;
         pos_curagrifparent = setposcur()-11;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_IRhoy") ) 
      {
         InAgrifParentDef = 10;
         pos_curagrifparent = setposcur()-18;
      }
      else if ( !strcasecmp(ident,"Agrif_Rhoz") ) 
      {
         InAgrifParentDef = 11;
         pos_curagrifparent = setposcur()-10;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_Rhoz") ) 
      {
         InAgrifParentDef = 12;
         pos_curagrifparent = setposcur()-17;
      }
      else if ( !strcasecmp(ident,"Agrif_IRhoz") ) 
      {
         InAgrifParentDef = 13;
         pos_curagrifparent = setposcur()-11;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_IRhoz") )
      {
         InAgrifParentDef = 14;
         pos_curagrifparent = setposcur()-18;
      }
      else if ( !strcasecmp(ident,"Agrif_NearCommonBorderX") )
      {
         InAgrifParentDef = 15;
         pos_curagrifparent = setposcur()-23;
      }
      else if ( !strcasecmp(ident,"Agrif_NearCommonBorderY") )
      {
         InAgrifParentDef = 16;
         pos_curagrifparent = setposcur()-23;
      }
      else if ( !strcasecmp(ident,"Agrif_NearCommonBorderZ") )
      {
         InAgrifParentDef = 17;
         pos_curagrifparent = setposcur()-23;
      }
      else if ( !strcasecmp(ident,"Agrif_DistantCommonBorderX") )
      {
         InAgrifParentDef = 18;
         pos_curagrifparent = setposcur()-26;
      }
      else if ( !strcasecmp(ident,"Agrif_DistantCommonBorderY") )
      {
         InAgrifParentDef = 19;
         pos_curagrifparent = setposcur()-26;
      }
      else if ( !strcasecmp(ident,"Agrif_DistantCommonBorderZ") )
      {
         InAgrifParentDef = 20;
         pos_curagrifparent = setposcur()-26;
      }
      else if ( !strcasecmp(ident,"Agrif_Get_parent_id") )
      {
         InAgrifParentDef = 21;
         pos_curagrifparent = setposcur()-19;
      }
      else if ( !strcasecmp(ident,"Agrif_Get_grid_id") )
      {
         InAgrifParentDef = 22;
         pos_curagrifparent = setposcur()-17;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_Iz") )
      {
         InAgrifParentDef = 23;
         pos_curagrifparent = setposcur()-15;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_Iy") )
      {
         InAgrifParentDef = 24;
         pos_curagrifparent = setposcur()-15;
      }
      else if ( !strcasecmp(ident,"Agrif_Parent_Ix") )
      {
         InAgrifParentDef = 25;
         pos_curagrifparent = setposcur()-15;
      }
      else if ( !strcasecmp(ident,"Agrif_Iz") )
      {
         InAgrifParentDef = 26;
         pos_curagrifparent = setposcur()-8;
      }
      else if ( !strcasecmp(ident,"Agrif_Iy") )
      {
         InAgrifParentDef = 27;
         pos_curagrifparent = setposcur()-8;
      }
      else if ( !strcasecmp(ident,"Agrif_Ix") )
      {
         InAgrifParentDef = 28;
         pos_curagrifparent = setposcur()-8;
      }
      else if ( !strcasecmp(ident,"Agrif_Nb_Fine_Grids") )
      {
         InAgrifParentDef = 29;
         pos_curagrifparent = setposcur()-19;
      }
      else if ( !strcasecmp(ident,"AGRIF_Nb_Step") )
      {
         InAgrifParentDef = 30;
         pos_curagrifparent = setposcur()-13;
      }
   }
}

/******************************************************************************/
/*                       AGRIF_0_ModifyTheAgrifFunction                       */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/*               Agrif_<toto>(variable) ====>     Agrif_<toto>(variable)      */
/*                                                                            */
/******************************************************************************/
void AGRIF_0_ModifyTheAgrifFunction(char *ident)
{
   if ( InAgrifParentDef != 0 )
          AGRIF_0_AgriffunctionModify(ident,InAgrifParentDef);
   /*                                                                         */
   InAgrifParentDef = 0;
}


/******************************************************************************/
/*                         AGRIF_0_AgriffunctionModify                        */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/* if whichone = 1 Agrif_parent ===>                                          */
/*                                                                            */
/* if whichone = 2 Agrif_Get_coarse_grid ===>                                 */
/*                                                                            */
/* if whichone = 3 Agrif_Rhox ===>                                            */
/*                                                                            */
/* if whichone = 4 Agrif_Parent_Rhox ===>                                     */
/*                                                                            */
/* if whichone = 5 Agrif_IRhox ===>                                           */
/*                                                                            */
/* if whichone = 6 Agrif_Parent_IRhox ===>                                    */
/*                                                                            */
/* if whichone = 7 Agrif_Rhoy ===>                                            */
/*                                                                            */
/* if whichone = 8 Agrif_Parent_Rhoy ===>                                     */
/*                                                                            */
/* if whichone = 9 Agrif_IRhoy ===>                                           */
/*                                                                            */
/* if whichone = 10 Agrif_Parent_IRhoy ===>                                   */
/*                                                                            */
/* if whichone = 11 Agrif_Rhoz ===>                                           */
/*                                                                            */
/* if whichone = 12 Agrif_Parent_Rhoz ===>                                    */
/*                                                                            */
/* if whichone = 13 Agrif_IRhoz ===>                                          */
/*                                                                            */
/* if whichone = 14 Agrif_Parent_IRhoz ===>                                   */
/*                                                                            */
/* if whichone = 15 Agrif_NearCommonBorderX ===>                              */
/*                                                                            */
/* if whichone = 16 Agrif_NearCommonBorderX ===>                              */
/*                                                                            */
/* if whichone = 17 Agrif_NearCommonBorderX ===>                              */
/*                                                                            */
/* if whichone = 18 Agrif_DistantCommonBorderX ===>                           */
/*                                                                            */
/* if whichone = 19 Agrif_DistantCommonBorderY ===>                           */
/*                                                                            */
/* if whichone = 20 Agrif_DistantCommonBorderZ ===>                           */
/*                                                                            */
/* if whichone = 21 Agrif_Get_parent_id ===>                                  */
/*                                                                            */
/* if whichone = 22 Agrif_Get_grid_id ===>                                    */
/*                                                                            */
/* if whichone = 23 Agrif_Parent_Iz ===>                                      */
/*                                                                            */
/* if whichone = 24 Agrif_Parent_Iy ===>                                      */
/*                                                                            */
/* if whichone = 25 Agrif_Parent_Ix ===>                                      */
/*                                                                            */
/* if whichone = 26 Agrif_Iz ===>                                             */
/*                                                                            */
/* if whichone = 27 Agrif_Iy ===>                                             */
/*                                                                            */
/* if whichone = 28 Agrif_Ix ===>                                             */
/*                                                                            */
/* if whichone = 29 Agrif_Nb_Fine_Grids ===>                                  */
/*                                                                            */
/* if whichone = 30 AGRIF_Nb_Step ===>                                        */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void AGRIF_0_AgriffunctionModify(char *ident,int whichone)
{
   char toprint[LONGNOM];

   if ( firstpass == 0 ) 
   {
      strcpy(toprint,"");
      pos_end = setposcur();
      fseek(fortranout,pos_curagrifparent,SEEK_SET);
      if ( whichone == 1 || whichone == 2 ) 
      {
         /*                                                                   */
         FindAndChangeNameToTabvars(ident,toprint,globliste,1);
         if ( !strcasecmp(ident,toprint) )
         {
            if ( ! globalvarofusefile ) RecordUseModulesVariables();
            /* la liste des use de cette subroutine                           */
            strcpy(toprint,"");
            FindAndChangeNameToTabvars(ident,
                                          toprint,globalvarofusefile,whichone);
         }
      }
      else if ( whichone == 3 ) /* Agrif_Rhox                                 */
      {
         sprintf(toprint,"REAL(Agrif_Curgrid %% spaceref(1))");
      }
      else if ( whichone == 4 ) /* Agrif_Parent_Rhox                          */
      {
         sprintf(toprint,"REAL(Agrif_Curgrid %% parent %% spaceref(1))");
      }
      else if ( whichone == 5 ) /* Agrif_Rhox                                 */
      {
         sprintf(toprint,"Agrif_Curgrid %% spaceref(1)");
      }
      else if ( whichone == 6 ) /* Agrif_Parent_Rhox                          */
      {
         sprintf(toprint,"Agrif_Curgrid %% parent %% spaceref(1)");
      }
      else if ( whichone == 7 ) /* Agrif_Rhoy                                 */
      {
         sprintf(toprint,"REAL(Agrif_Curgrid %% spaceref(2))");
      }
      else if ( whichone == 8 ) /* Agrif_Parent_Rhoy                          */
      {
         sprintf(toprint,"REAL(Agrif_Curgrid %% parent %% spaceref(2))");
      }
      else if ( whichone == 9 ) /* Agrif_Rhoy                                 */
      {
         sprintf(toprint,"Agrif_Curgrid %% spaceref(2)");
      }
      else if ( whichone == 10 ) /* Agrif_Parent_Rhoy                         */
      {
         sprintf(toprint,"Agrif_Curgrid %% parent %% spaceref(2)");
      }
      else if ( whichone == 11 ) /* Agrif_Rhoz                                */
      {
         sprintf(toprint,"REAL(Agrif_Curgrid %% spaceref(3))");
      }
      else if ( whichone == 12 ) /* Agrif_Parent_Rhoz                         */
      {
         sprintf(toprint,"REAL(Agrif_Curgrid %% parent %% spaceref(3))");
      }
      else if ( whichone == 13 ) /* Agrif_Rhoz                                */
      {
         sprintf(toprint,"Agrif_Curgrid %% spaceref(3)");
      }
      else if ( whichone == 14 ) /* Agrif_Parent_Rhoz                         */
      {
         sprintf(toprint,"Agrif_Curgrid %% parent %% spaceref(3)");
      }
      else if ( whichone == 15 ) /* Agrif_NearCommonBorderX                   */
      {
         sprintf(toprint,"Agrif_Curgrid %% NearRootBorder(1)");
      }
      else if ( whichone == 16 ) /* Agrif_NearCommonBorderY                   */
      {
         sprintf(toprint,"Agrif_Curgrid %% NearRootBorder(2)");
      }
      else if ( whichone == 17 ) /* Agrif_NearCommonBorderZ                   */
      {
         sprintf(toprint,"Agrif_Curgrid %% NearRootBorder(3)");
      }
      else if ( whichone == 18 ) /* Agrif_NearCommonBorderX                   */
      {
         sprintf(toprint,"Agrif_Curgrid %% DistantRootBorder(1)");
      }
      else if ( whichone == 19 ) /* Agrif_NearCommonBorderY                   */
      {
         sprintf(toprint,"Agrif_Curgrid %% DistantRootBorder(2)");
      }
      else if ( whichone == 20 ) /* Agrif_NearCommonBorderZ                   */
      {
         sprintf(toprint,"Agrif_Curgrid %% DistantRootBorder(3)");
      }
      else if ( whichone == 21 ) /* Agrif_Get_parent_id                       */
      {
         sprintf(toprint,"Agrif_Curgrid %% parent %% grid_id");
      }
      else if ( whichone == 22 ) /*  Agrif_Get_grid_id                        */
      {
         sprintf(toprint,"Agrif_Curgrid %% grid_id");
      }
      else if ( whichone == 23 ) /*  Agrif_Parent_Iz                          */
      {
         sprintf(toprint,"Agrif_Curgrid %% parent %% ix(3)");
      }
      else if ( whichone == 24 ) /*  Agrif_Parent_Iy                          */
      {
         sprintf(toprint,"Agrif_Curgrid %% parent %% ix(2)");
      }
      else if ( whichone == 25 ) /*  Agrif_Parent_Ix                          */
      {
         sprintf(toprint,"Agrif_Curgrid %% parent %% ix(1)");
      }
      else if ( whichone == 26 ) /* Agrif_Iz                                  */
      {
         sprintf(toprint,"Agrif_Curgrid %% ix(3)");
      }
      else if ( whichone == 27 ) /* Agrif_Iy                                  */
      {
         sprintf(toprint,"Agrif_Curgrid %% ix(2)");
      }
      else if ( whichone == 28 ) /* Agrif_Ix                                  */
      {
         sprintf(toprint,"Agrif_Curgrid %% ix(1)");
      }
      else if ( whichone == 29 ) /* Agrif_Nb_Fine_Grids                       */
      {
         sprintf(toprint,"Agrif_nbfixedgrids");
      }
      else if ( whichone == 30 ) /* AGRIF_Nb_Step                             */
      {
         sprintf(toprint,"Agrif_Curgrid %% ngridstep");
      }
      if ( whichone == 1 || whichone == 2 ) 
      {
         tofich(fortranout,toprint,0);
      }
      else
      {
         if( fortran77 == 0 ) fprintf(fortranout," & \n");
         else fprintf(fortranout,"\n     & ");
         fprintf(fortranout,"%s",toprint);
      }
   }
}



/******************************************************************************/
/*                          AGRIF_0_AddUseAgrifInModuleDeclaration            */
/******************************************************************************/
/* Add use Agrif_Util at the beginning of the subroutine definition           */
/* if it is necessary                                                         */
/******************************************************************************/
/*                                                                            */
/*       subroutine sub            |  subroutine sub                          */
/*                                 |  USE Agrif_Util                          */
/*       implicit none             |  implicit none                           */
/*       ...                       |  ...                                     */
/*       ... Agrif_<something>     |  ... Agrif_<something>                   */
/*       ...                       |  ...                                     */
/*       end                       |  end                                     */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void  AGRIF_0_AddUseAgrifInModuleDeclaration()
{
  listusemodule *newmodule;
  int out;

   if ( firstpass == 0 ) 
   {
      out = 1 ;
      /* We should see if agrif_tabvars is the only                           */
      /*    necessary tools in the agrif librairy                             */
      newmodule = listofmodulebysubroutine;
      while( newmodule && out == 1 && !listofsubroutinewhereagrifisused )
      {
         if ( !strcasecmp(newmodule->cursubroutine,subroutinename) ||
              !strcasecmp(newmodule->cursubroutine," ")  )
         {
            if ( !strcasecmp(newmodule->charusemodule,"Agrif_Util") ) out = 0 ;
         }
         newmodule = newmodule ->suiv;
      }

      if ( out == 0 || listofsubroutinewhereagrifisused ) 
                               fprintf(fortranout,"\n       USE Agrif_Util \n");
      else fprintf(fortranout,
                           "\n       USE Agrif_types\n");

   }
}
