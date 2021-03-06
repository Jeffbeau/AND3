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
/*                    OPTI_1_cleanlistvarfordoloop                            */
/******************************************************************************/
/* Firstpass 1                                                                */
/* We should clean all the list used for the do loop OPTImization             */
/* if endsuborfunc = 1 we are at the end of the subroutine                    */
/* if endsuborfunc = 0 we are at the end of the function                      */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_1_cleanlistvarfordoloop(int endsuborfunc)
{
   listvar *tmplist;
   
   if ( firstpass == 1 ) 
   {
      if ( fortran77 == 1 ) UpdatevarofsubroutinelisteWithcommonlist();
      if ( fortran77 == 1 ) COM_1_UpdateparameterlistWithlistvarindoloop();
      if ( fortran77 == 1 ) COM_1_UpdateGloblisteWithcommonlist();
      if ( endsuborfunc == 1 ) CompleteThelistvarindoloop();
      if ( fortran77 == 0 ) UpdateIndiceTabvarsofGlobliste();
      else UpdateIndiceTabvarsofGloblisteFromCommon();
      CleanThelistvarindoloop();
      CleanFromThelistvarindoloopTheAgrifSubArguments();
      tmplist = (listvar *)NULL;
      if ( fortran77 == 1 ) tmplist = duplicatelistvar(varofsubroutineliste);
      if ( fortran77 == 1 ) varsubroutine = AddListvarToListvar
                                                      (tmplist,varsubroutine,1);
      CleanThelistvarofsubroutineliste();
      if ( fortran77 == 1 ) COM_1_UpdatevarsubroutineWithvarofsubroutinelist();
   }
}

/******************************************************************************/
/*                    OPTI_1_ajoutevarindoloop                                */
/******************************************************************************/
/* Firstpass 1                                                                */
/* We should complete the listvarindoloop                                     */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void OPTI_1_ajoutevarindoloop(char *ident)
{
   /* In the first pass we record all variables presents in the do loop       */
   if (firstpass == 1 && insubroutinedeclare == 1 ) ajoutevarindoloop(ident);
}

/******************************************************************************/
/*                        AJOUTEVARINDOLOOP                                   */
/******************************************************************************/
/* This subroutine is used to add a listvar to  listvarindoloop               */
/******************************************************************************/
void ajoutevarindoloop (char *name)
{
  listvar *newvar;
  listvar *tmpvar;
  int out;
  
  if ( !listvarindoloop )
  {
      newvar=(listvar *)malloc(sizeof(listvar));
      newvar->var=(variable *)malloc(sizeof(variable));
      newvar->suiv = NULL;
      strcpy(newvar->var->oldname,"");
      strcpy(newvar->var->nomvar,name);
      strcpy(newvar->var->modulename,subroutinename);
      newvar->var->pointedvar=pointedvar;
      newvar->var->indicetabvars=0;
      listvarindoloop = newvar ;
  }
  else
  {
      /* We should verify that this variable did not added                    */
      tmpvar = listvarindoloop;
      out = 0 ;
      while (tmpvar && out == 0 )
      {
         if ( !strcasecmp(tmpvar->var->nomvar,name) && 
              !strcasecmp(tmpvar->var->modulename,subroutinename)) out  = 1 ; 
         else tmpvar = tmpvar->suiv;
      }
      if ( out == 0 ) 
      {
         newvar=(listvar *)malloc(sizeof(listvar));
         newvar->var=(variable *)malloc(sizeof(variable));
         strcpy(newvar->var->oldname,"");
         strcpy(newvar->var->nomvar,name);
         strcpy(newvar->var->modulename,subroutinename);
         newvar->var->pointedvar=pointedvar;
         newvar->var->indicetabvars=0;
         newvar->suiv = listvarindoloop;
         listvarindoloop = newvar;
      }
  }
}

/******************************************************************************/
/*                        AJOUTEVARINDOLOOP_DEFINEDIMENSION                   */
/******************************************************************************/
/* This subroutine is used to add a listvar to  listvarindoloop               */
/******************************************************************************/
void ajoutevarindoloop_definedimension (char *name)
{
  listvar *newvar;
  listvar *tmpvar;
  listvar *tmpvarprec;
  int out;
  int tablemeet;
  
  if ( !listvarindoloop )
  {
      newvar=(listvar *)malloc(sizeof(listvar));
      newvar->var=(variable *)malloc(sizeof(variable));
      newvar->suiv = NULL;
      strcpy(newvar->var->oldname,"");
      strcpy(newvar->var->nomvar,name);
      strcpy(newvar->var->modulename,subroutinename);
      newvar->var->indicetabvars=0;
      newvar->var->pointedvar=pointedvar;
      listvarindoloop = newvar ;
  }
  else
  {
      /* We should verify that this variable did not added                    */
      tmpvarprec = (listvar *)NULL;
      tmpvar = listvarindoloop;
      out = 0 ;
      tablemeet = 0 ;
      while (tmpvar && out == 0 )
      {
         if ( tablemeet == 0 && tmpvar->var->nbdim != 0 ) tablemeet = 1 ;
         /*                                                                   */
         if ( !strcasecmp(tmpvar->var->nomvar,name) && 
              !strcasecmp(tmpvar->var->modulename,subroutinename)) 
         {
            out  = 1 ;
            /* if this variable has been define before a table we doi nothing */
            /*    else we should remove it                                    */
            if ( tablemeet == 1 )
            {
               tmpvarprec->suiv = tmpvar -> suiv;
               out = 2;
            }
         }
         else 
         {
            tmpvarprec = tmpvar;
            tmpvar = tmpvar->suiv;
         }
      }
      if ( out == 2 || out == 0 ) 
      {
         newvar=(listvar *)malloc(sizeof(listvar));
         newvar->var=(variable *)malloc(sizeof(variable));
         strcpy(newvar->var->nomvar,name);
         strcpy(newvar->var->oldname,"");
         newvar->var->indicetabvars=0;
         strcpy(newvar->var->modulename,subroutinename);
         newvar->var->pointedvar=pointedvar;
         /* we should find this new variable to know the tabvars indice       */
         if ( variableisglobal(newvar, globliste) == 1 )
         {
            newvar->suiv = listvarindoloop;
            listvarindoloop = newvar;
         }
         else if ( variableisglobal(newvar, globalvarofusefile) == 1 )
         {
            newvar->suiv = listvarindoloop;
            listvarindoloop = newvar;
         }
         else
         {
            free(newvar);
         }
     }
  }
}

/******************************************************************************/
/*        CleanFromThelistvarindoloopTheAgrifSubArguments                     */
/******************************************************************************/
/* This subroutine is to remove from the listvarindoloop all variables        */
/* which has been used in Agrif argument in order to avoid the                */
/* optimization code on Agrif function or subroutines                         */
/******************************************************************************/
void  CleanFromThelistvarindoloopTheAgrifSubArguments()
{
   listnom *parcours;
   listvar *parcoursvar;
   listvar *parcoursvarprec;
   
   parcoursvarprec = (listvar *)NULL;
   parcoursvar = listvarindoloop;
   while ( parcoursvar )
   {
      if ( !strcasecmp(parcoursvar->var->modulename,subroutinename) )
      {
         parcours = Listofvariableinagriffunction;
         while (parcours && strcasecmp(parcoursvar->var->nomvar,parcours->nom) )
         {
            parcours = parcours->suiv;
         }
         if ( parcours )
         {
            /* if we found the name in the listvarindoloop and                */
            /* Listofvariableinagriffunction we should remove it from         */
            /* listvarindoloop                                                */
            if ( parcoursvar == listvarindoloop )
            {
               listvarindoloop = listvarindoloop -> suiv;
               parcoursvar = listvarindoloop;
            }
            else
            {
               parcoursvarprec->suiv = parcoursvar->suiv;
               parcoursvar = parcoursvar->suiv;
            }
         }
         else
         {
            parcoursvarprec = parcoursvar;
            parcoursvar = parcoursvar ->suiv;
         }
      }
      else
      {
         parcoursvarprec = parcoursvar;
         parcoursvar = parcoursvar ->suiv;
      }
   }
   
}


/******************************************************************************/
/*                      CleanThelistvarindoloop                               */
/******************************************************************************/
/* This subroutine is to remove from the listvarindoloop all variables        */
/* which has not been declared as table in the globliste                      */
/******************************************************************************/
void CleanThelistvarindoloop ()
{
  listvar *newvar;
  listvar *newvarPrec;
  listvar *tmpglobvar;
  listallocate *parcoursallocate;
  listnamelist *newnamelist;
  int not_remove;

  RecordUseModulesVariables();
  RecordUseModulesUseModulesVariables();
  /*                                                                          */
  not_remove = 0 ;
  newvarPrec = (listvar *)NULL;
  newvar = listvarindoloop;
  while ( newvar )
  {
  if ( !strcasecmp(newvar->var->modulename,subroutinename))
  {
     not_remove = 0;
     if ( Variableshouldberemove(newvar->var->nomvar) == 0 )
     {
/******************************************************************************/
/*                      look in the globliste                                 */
/******************************************************************************/
/******************************************************************************/
/*                      look in the varofsubroutineliste                      */
/******************************************************************************/
        tmpglobvar = varofsubroutineliste;
        while ( tmpglobvar && not_remove == 0 )
        {
           if ( !strcasecmp(tmpglobvar->var->nomvar,newvar->var->nomvar) &&
                !strcasecmp
                   (tmpglobvar->var->modulename,newvar->var->modulename)
              )
               not_remove = 2;
          else tmpglobvar = tmpglobvar->suiv;
       }

     if (not_remove == 0 ) tmpglobvar = globliste;
     else tmpglobvar = (listvar *)NULL;

     while ( tmpglobvar && not_remove == 0 )
     {
        if ( !strcasecmp(tmpglobvar->var->nomvar,newvar->var->nomvar) )
        {
           not_remove = 1;
           /* Now we should give the definition of the variable in the        */
           /*    table listvarindoloop                                        */
           strcpy(newvar->var->typevar,tmpglobvar->var->typevar);
           strcpy(newvar->var->dimchar,tmpglobvar->var->dimchar);
           newvar->var->nbdim = tmpglobvar->var->nbdim;
           newvar->var->dimensiongiven = tmpglobvar->var->dimensiongiven;
           newvar->var->typegiven = tmpglobvar->var->typegiven;
           newvar->var->allocatable = tmpglobvar->var->allocatable;
           newvar->var->pointerdeclare = tmpglobvar->var->pointerdeclare;
           newvar->var->indicetabvars = tmpglobvar->var->indicetabvars;
           strcpy(newvar->var->precision,tmpglobvar->var->precision);
           strcpy(newvar->var->readedlistdimension,
                                          tmpglobvar->var->readedlistdimension);
           DecomposeTheName(newvar->var->readedlistdimension);
        }
        else tmpglobvar = tmpglobvar->suiv;
     }
     
/******************************************************************************/
/*                      look in the globparam                                 */
/******************************************************************************/
     if ( not_remove == 0 )
     {
        tmpglobvar = globparam;
        while ( tmpglobvar && not_remove == 0 )
        {
           if ( !strcasecmp(tmpglobvar->var->nomvar,newvar->var->nomvar) &&
                !strcasecmp(tmpglobvar->var->subroutinename,
                                                newvar->var->modulename) 
               ) not_remove = 2;
           else tmpglobvar = tmpglobvar->suiv;
        }
     }
     
     if ( not_remove == 0 )
     {
/******************************************************************************/
/*                      look in the listenamelist                             */
/******************************************************************************/
        newnamelist = listenamelist;
        while ( newnamelist && not_remove == 0 )
        {
           if ( !strcasecmp(newnamelist->name,newvar->var->nomvar)) not_remove = 2;
           else newnamelist = newnamelist->suiv;
        }
     }

     if ( not_remove == 0 )
     {
/******************************************************************************/
/*                      look in the varofsubroutineliste                      */
/******************************************************************************/
        tmpglobvar = varofsubroutineliste;
        while ( tmpglobvar && not_remove == 0 )
        {
           if ( !strcasecmp(tmpglobvar->var->nomvar,newvar->var->nomvar) &&
                !strcasecmp(tmpglobvar->var->modulename,
                                                newvar->var->modulename)
              ) not_remove = 2;
          else tmpglobvar = tmpglobvar->suiv;
       }
     }
/******************************************************************************/
/*            look in the .dependfile and .dependparameterfile                */
/******************************************************************************/
     if ( not_remove == 0 && not_remove == 0 )
     {
        /* la liste des use de cette subroutine                               */
        not_remove = 0 ;

        if ( variableisparameterglobal(newvar,tmpparameterlocallist) == 1 )
        {
           not_remove = 2 ;
        }
        else if ( variableisglobal(newvar,globalvarofusefile) == 1 )
        {
           not_remove = 1 ;
           DecomposeTheName(newvar->var->readedlistdimension);
        }

/******************************************************************************/
/*    look in the .dependfile and .dependparameterfile of USE modules         */
/******************************************************************************/
        if ( not_remove == 0 )
        {
           if ( variableisparameterglobal(newvar,tmpparameterlocallist2) == 1 )
           {
              not_remove = 2 ;
           }
           else if ( variableisglobal(newvar, globalvarofusefile2) == 1 )
           {
              not_remove = 1 ;
              DecomposeTheName(newvar->var->readedlistdimension);
           }
        }
     }
/******************************************************************************/
/*                          look if pointer variable                          */
/******************************************************************************/
     /* if this variable is a pointer we should remove it                     */
     if ( not_remove == 1 && newvar->var->pointerdeclare == 1 )
     {
        not_remove = 2;
     }
     /* if this variable is an allocatable var we should remove it            */
     if ( not_remove == 1 && newvar->var->allocatable == 1 )
     {
        not_remove = 2;
     }
/******************************************************************************/
/*                          look in the AllocateList                          */
/******************************************************************************/
     /* if this variable has been used in a allocate we should remove it      */
     if ( not_remove == 1 && newvar->var->nbdim != 0 )
     {
        parcoursallocate = AllocateList;
        while ( parcoursallocate && not_remove == 1 )
        {
           if ( !strcasecmp(parcoursallocate->nomvar,newvar->var->nomvar) &&
                !strcasecmp(parcoursallocate->subroutine,subroutinename)
           ) not_remove = 2;
           else parcoursallocate = parcoursallocate->suiv;
        }
     }
     /*                                                                       */
     } /* end of strcasecmp(newvar->var->nomvar,"") */
     else
     {
        not_remove = 2;
     }
/******************************************************************************/
/*                          REMOVE                                            */
/******************************************************************************/
     if ( (   not_remove == 0 || not_remove == 2 ) && 
              newvar->var->pointedvar == 0 
        )
     {
        if ( newvar == listvarindoloop )
        {
           listvarindoloop = listvarindoloop->suiv;
           newvar = listvarindoloop;
        }
        else
        {
           newvarPrec->suiv = newvar->suiv;
           newvar = newvarPrec->suiv;
        }
     }
     else
     {        
        /*                                                                    */
        newvarPrec = newvar;
        newvar = newvar->suiv ;
     }
  }
  else
  {
     newvarPrec = newvar;
     newvar = newvar->suiv;
  }
  }
}


/******************************************************************************/
/*                        ModifyThelistvarindoloop                            */
/******************************************************************************/
/* This subroutine is to give the old name to the which has been              */
/* declared as USE MOD, U => V in this case we should replace in the          */
/* name V by the old name U in the listvarindoloop                            */
/******************************************************************************/
void  ModifyThelistvarindoloop()
{
  listvar *newvar;
      
  newvar = listvarindoloop;
  while ( newvar )
  {
     if ( strcasecmp(newvar->var->oldname,"") )
     {
        strcpy(newvar->var->nomvar,newvar->var->oldname);
     }
     newvar = newvar->suiv;
  }
}

/******************************************************************************/
/*                          CompleteThelistvarindoloop                        */
/******************************************************************************/
/* This subroutine is to add to the listvarindoloop all variables which       */
/* has been declared as USE MOD, U => V in this case we should replace        */
/* in the listvarindoloop the word U by the word V                            */
/******************************************************************************/
void  CompleteThelistvarindoloop()
{
  listvar *newvar;
  listvarpointtovar *pointtmplist;
  listcouple *coupletmp;
  int outvar;
      
  pointtmplist = Listofvarpointtovar;
  
  while ( pointtmplist )
  {
      coupletmp = pointtmplist->couple;
      while ( coupletmp )
      {
         newvar = listvarindoloop;
         outvar = 0 ;
         while ( newvar && outvar == 0)
         {
           /* we should find the same variable name in the same subroutine    */
           if ( !strcasecmp(newvar->var->nomvar,coupletmp->namevar) &&
                !strcasecmp(newvar->var->modulename,
                                       pointtmplist->cursubroutine) &&
                 strcasecmp(coupletmp->namepointedvar,"") 
              )
           {
              outvar = 1;
              strcpy(newvar->var->oldname,newvar->var->nomvar);
              strcpy(newvar->var->nomvar,coupletmp->namepointedvar);
           }
           else
           {      
              newvar = newvar->suiv;
           }
         }
         coupletmp = coupletmp->suiv;     
     }
     pointtmplist = pointtmplist->suiv;
  }
}
