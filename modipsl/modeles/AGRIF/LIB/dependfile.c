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
/*            Creation and modification of .dependfile                        */
/******************************************************************************/
/*  .dependnbxnby            : this file contains tabvars indices of variables*/
/*                                given in agrif.in as number of cells        */
/*  .dependuse<module>       : this file contains all modules used in the     */
/*                                current file                                */
/*  .dependparameter<module> : this file contains all parmeters defined in    */
/*                                the current file                            */
/*  .depend<module name>     : this file contains all globals variables       */
/*                                informations (name, dim, etc ...)           */
/*  .dependavailable         : this file contains all tabvars indices which   */
/*                                are not used.                               */
/******************************************************************************/


/******************************************************************************/
/*                 Writethedependnbxnbyfile                                   */
/******************************************************************************/
/* This subroutine is used to create the .dependnbxnby                        */
/******************************************************************************/
/*                                                                            */
/*                     .dependnbxnby                                          */
/*                                                                            */
/*                     nbmaillesX                                             */
/*                     nbmaillesY                                             */
/*                     nbmaillesZ                                             */
/*                                                                            */
/******************************************************************************/
void Writethedependnbxnbyfile()
{
  FILE *dependfileoutput;
  listvar *parcours;
  int out;
   
  /* We are looking for all the variable of the current filetoparse file      */
  /*    in the globliste                                                      */
  parcours =globliste;
  out = 0;
  while (parcours && out == 0 )
  {
     if ( !strcasecmp(parcours->var->nomvar,nbmaillesX) ) out = 1;
     else parcours = parcours->suiv;
  }
  if ( out == 1 ) 
  {
     dependfileoutput = fopen(".dependnbxnby","w");
     fprintf(dependfileoutput,"%d\n",parcours->var->indicetabvars);
     IndicenbmaillesX = parcours->var->indicetabvars;

     if ( dimprob > 1 )
     {
        parcours =globliste;
        out = 0;
        while (parcours && out == 0 )
        {
           if ( !strcasecmp(parcours->var->nomvar,nbmaillesY) ) out = 1;
           else parcours = parcours->suiv;
        }
        if ( out == 1 ) 
        {
           fprintf(dependfileoutput,"%d\n",parcours->var->indicetabvars);
           IndicenbmaillesY = parcours->var->indicetabvars;
        }
     }

     if ( dimprob > 2 )
     {
        parcours =globliste;
        out = 0;
        while (parcours && out == 0 )
        {
           if ( !strcasecmp(parcours->var->nomvar,nbmaillesZ) ) out = 1;
           else parcours = parcours->suiv;
        }
        if ( out == 1 ) 
        {
           fprintf(dependfileoutput,"%d\n",parcours->var->indicetabvars);
           IndicenbmaillesZ = parcours->var->indicetabvars;
        }
     }

     if ( out == 1 ) fclose(dependfileoutput);
   }
}

/******************************************************************************/
/*                 Readthedependnbxnbyfile                                    */
/******************************************************************************/
/* This subroutine is used to create the .dependnbxnby                        */
/******************************************************************************/
/*                                                                            */
/*                     .dependnbxnby                                          */
/*                                                                            */
/*                     nbmaillesX                                             */
/*                     nbmaillesY                                             */
/*                     nbmaillesZ                                             */
/*                                                                            */
/******************************************************************************/
void Readthedependnbxnbyfile()
{
  FILE *dependfileoutput;
   
  if ((dependfileoutput = fopen(".dependnbxnby","r"))!=NULL) 
  {
     fscanf(dependfileoutput,"%d\n",&IndicenbmaillesX);
     if ( dimprob > 1 ) fscanf(dependfileoutput,"%d\n",&IndicenbmaillesY);
     if ( dimprob > 2 ) fscanf(dependfileoutput,"%d\n",&IndicenbmaillesZ);
     fclose(dependfileoutput);
  }
}

/******************************************************************************/
/*                     Writethedependlistofmoduleused                         */
/******************************************************************************/
/* This subroutine is used to create the .dependuse<module>                   */
/******************************************************************************/
/*                                                                            */
/*               .dependuse<name>                                             */
/*                                                                            */
/*               mod1                                                         */
/*               mod2                                                         */
/*                                                                            */
/******************************************************************************/
void Writethedependlistofmoduleused(char *NameTampon )
{
  FILE *dependfileoutput;
  listusemodule *parcours;
  char ligne[LONGNOM];
   
  if ( listofmodulebysubroutine )
  {
     sprintf(ligne,".dependuse%s",NameTampon);
     dependfileoutput = fopen(ligne,"w");
     /* We are looking for all the variable of the current filetoparse file   */
     /*    in the globliste                                                   */
     parcours = listofmodulebysubroutine;
     while (parcours)
     {
        fprintf(dependfileoutput,"%s\n",parcours->usemodule);
        parcours = parcours->suiv;
     }
     fclose(dependfileoutput);
  }
}

/******************************************************************************/
/*                    Readthedependlistofmoduleused                           */
/******************************************************************************/
/* This subroutine is used to create the .dependuse<module>                   */
/******************************************************************************/
/*                                                                            */
/*               .dependuse<name>                                             */
/*                                                                            */
/*               mod1                                                         */
/*               mod2                                                         */
/*                                                                            */
/******************************************************************************/
void Readthedependlistofmoduleused(char *NameTampon)
{
  FILE *dependfileoutput;
  listusemodule *parcours;
  char ligne[LONGNOM];

  sprintf(ligne,".dependuse%s",NameTampon);
  
  tmpuselocallist = (listusemodule *)NULL;
  if ((dependfileoutput = fopen(ligne,"r"))==NULL) 
  {
  }
  else
  {
    /* if the file exist we should verify that this file has changed          */
      while (!feof(dependfileoutput))
      {
         parcours=(listusemodule *)malloc(sizeof(listusemodule));
         fscanf(dependfileoutput,"%s\n",parcours->usemodule);

         parcours->suiv = tmpuselocallist;
         tmpuselocallist = parcours;

         parcours = NULL;
      }
      fclose(dependfileoutput);
  }
}


/******************************************************************************/
/*                        WritedependParameterList                            */
/******************************************************************************/
/* This subroutine is used to create the .dependparameter<name>               */
/******************************************************************************/
/*                                                                            */
/*               .dependparameter<name>                                       */
/*                                                                            */
/*               mod1                                                         */
/*               mod2                                                         */
/*                                                                            */
/******************************************************************************/
void WritedependParameterList(char *NameTampon )
{
  FILE *dependfileoutput;
  listvar *parcours;
  char ligne[LONGNOM];
   
  if ( globparam )
  {
  sprintf(ligne,".dependparameter%s",NameTampon);
  dependfileoutput = fopen(ligne,"w");

  parcours =globparam;
  while (parcours)
  {
     fprintf(dependfileoutput,"%s\n",parcours->var->nomvar);
     fprintf(dependfileoutput,"%s\n",parcours->var->modulename);
     parcours = parcours->suiv;
  }
  fclose(dependfileoutput);
  }
}


/******************************************************************************/
/*                         ReaddependParameterList                            */
/******************************************************************************/
/* This subroutine is used to create the .dependparameter<name>               */
/******************************************************************************/
/*                                                                            */
/*               .dependparameter<name>                                       */
/*                                                                            */
/*               mod1                                                         */
/*               mod2                                                         */
/*                                                                            */
/******************************************************************************/
listparameter *ReaddependParameterList(char *NameTampon,listparameter *listout)
{
  FILE *dependfileoutput;
  listparameter *parcours;
  char ligne[LONGNOM];

  sprintf(ligne,".dependparameter%s",NameTampon);

  if ((dependfileoutput = fopen(ligne,"r"))==NULL) 
  {
  }
  else
  {
    /* if the file exist we should verify that this file has changed          */
      while (!feof(dependfileoutput))
      {
         parcours=(listparameter *)malloc(sizeof(listparameter));
         fscanf(dependfileoutput,"%s\n",parcours->name);
         fscanf(dependfileoutput,"%s\n",parcours->modulename);

         parcours->suiv = listout;
         listout = parcours;

         parcours = NULL;
      }
      fclose(dependfileoutput);
  }
  return listout;
}

/******************************************************************************/
/*                   Writethedependfile                                       */
/******************************************************************************/
/* This subroutine is used to create the .depend<name>                        */
/******************************************************************************/
/*                                                                            */
/*                     .depend<name>                                          */
/*                                                                            */
/*                      REAL                                                  */
/*                      Variable                                              */
/*                      char dimension or T                                   */
/*                      table dimension                                       */
/*                      is type given                                         */
/*                      precision or T                                        */
/*                      initial value or T                                    */
/*                      indice in the tabvars                                 */
/*                      listdimension or T                                    */
/*                      -------------------------                             */
/*                                                                            */
/******************************************************************************/
void Writethedependfile(char *NameTampon, listvar *input )
{
  FILE *dependfileoutput;
  listvar *parcours;
  listdim *dims;
  char ligne[LONGNOM];
  char listdimension[LONGNOM];
  char curname[LONGNOM];
  int out;
   
  if ( input )
  {
  sprintf(ligne,".depend%s",NameTampon);
  dependfileoutput = fopen(ligne,"w");
  /* We are looking for all the variable of the current filetoparse file      */
  /*    in the input                                                          */
  parcours =input;
  out = 0;
  strcpy(curname,"");
  while (parcours && out == 0 )
  {
     if ( !strcasecmp(parcours->var->modulename,NameTampon) ||
          !strcasecmp(parcours->var->commonname,NameTampon) )
     {
        /*                                                                    */
        if ( strcmp(curname,"") && 
             !strcmp(curname,parcours->var->nomvar) ) out = 1 ;
        if ( !strcmp(curname,"") ) strcpy(curname,parcours->var->nomvar);
        /*                                                                    */
        if ( out == 0 )
        {
           fprintf(dependfileoutput,"%s\n",parcours->var->typevar);
           fprintf(dependfileoutput,"%s\n",parcours->var->nomvar);
           if ( strcasecmp(parcours->var->dimchar,"") )
           {
              fprintf(dependfileoutput,"%s\n",parcours->var->dimchar);
           }
           else
           {
              fprintf(dependfileoutput,"T\n");
           }
           if ( strcasecmp(parcours->var->commoninfile,"") )
           {
              fprintf(dependfileoutput,"%s\n",parcours->var->commoninfile);
           }
           else
           {
              fprintf(dependfileoutput,"T\n");
           }
           if ( strcasecmp(parcours->var->commonname,"") )
           {
              fprintf(dependfileoutput,"%s\n",parcours->var->commonname);
           }
           else
           {
              fprintf(dependfileoutput,"T\n");
           }
           fprintf(dependfileoutput,"%d\n",parcours->var->nbdim);
           fprintf(dependfileoutput,"%d\n",parcours->var->dimensiongiven);
           fprintf(dependfileoutput,"%d\n",parcours->var->typegiven);
           fprintf(dependfileoutput,"%d\n",parcours->var->allocatable);
           fprintf(dependfileoutput,"%d\n",parcours->var->pointerdeclare);
           if ( strcasecmp(parcours->var->precision,"") )
           {
              fprintf(dependfileoutput,"%s\n",parcours->var->precision); 
           }
           else
           {
              fprintf(dependfileoutput,"T\n");
           }
           if ( strcasecmp(parcours->var->initialvalue,"") )
           {
              fprintf(dependfileoutput,"%s\n",parcours->var->initialvalue); 
           }
           else
           {
              fprintf(dependfileoutput,"T\n");
           }
           fprintf(dependfileoutput,"%d\n",parcours->var->indicetabvars);
           if ( parcours->var->dimensiongiven == 1 )
           {
              dims = parcours->var->dimension;
              strcpy(listdimension,"");
              while (dims)
              {
                 sprintf(ligne,"%s:%s",dims->dim.first,dims->dim.last);
                 strcat(listdimension,ligne);
                 if ( dims->suiv )
                 {
                    strcat(listdimension,",");     
                 }
                 dims = dims->suiv;
              }
              fprintf(dependfileoutput,"%s\n",listdimension);    
           }
           else
           {
              fprintf(dependfileoutput,"T\n");
           }
           fprintf(dependfileoutput,"------------------------\n");
        }
     }
     parcours = parcours->suiv;
  }
  fclose(dependfileoutput);
  }
}

/******************************************************************************/
/*                        Recordtmplocallist                                  */
/******************************************************************************/
/* This subroutine is used to read the .depend<name> if the                   */
/* filetoparse has been read before                                           */
/* if it is the case we record the .depend<name> in the tmplocallist          */
/******************************************************************************/
/*                                                                            */
/*           .depend<name> -------->  tmplocallist = list of var              */
/*                                                                            */
/*        not.depend<name> -------->  tmplocallist = NULL                     */
/*                                                                            */
/******************************************************************************/
void Recordtmplocallist( char *NameTampon)
{
  FILE *dependfileoutput;
  listvar *parcours;
  listvar *parcoursprec;
  char ligne[LONGNOM];
  char nothing[LONGNOM];

  parcoursprec = (listvar *)NULL;
  sprintf(ligne,".depend%s",NameTampon);
  if ((dependfileoutput = fopen(ligne,"r"))==NULL) 
  {
    /* if the file doesn't exist it means that it is the first time           */
    /*    we tried to parse this file                                         */
  }
  else
  {
    /* if the file exist we should verify that this file has changed          */
      while (!feof(dependfileoutput))
      {
         parcours=(listvar *)malloc(sizeof(listvar));
         parcours->var=(variable *)malloc(sizeof(variable));
         fscanf(dependfileoutput,"%s\n",parcours->var->typevar);
         fscanf(dependfileoutput,"%s\n",parcours->var->nomvar);
         fscanf(dependfileoutput,"%s\n",parcours->var->dimchar);
         if ( !strcasecmp(parcours->var->dimchar,"T") )
         { 
            strcpy(parcours->var->dimchar,"");
         }
         fscanf(dependfileoutput,"%s\n",parcours->var->commoninfile);
         if ( !strcasecmp(parcours->var->commoninfile,"T") )
         { 
            strcpy(parcours->var->commoninfile,"");
         }
         fscanf(dependfileoutput,"%s\n",parcours->var->commonname);
         if ( !strcasecmp(parcours->var->commonname,"T") )
         { 
            strcpy(parcours->var->commonname,"");
         }
         fscanf(dependfileoutput,"%d\n",&parcours->var->nbdim);
         fscanf(dependfileoutput,"%d\n",&parcours->var->dimensiongiven);
         fscanf(dependfileoutput,"%d\n",&parcours->var->typegiven);
         fscanf(dependfileoutput,"%d\n",&parcours->var->allocatable);
         fscanf(dependfileoutput,"%d\n",&parcours->var->pointerdeclare);
         fscanf(dependfileoutput,"%[^\n] \n",parcours->var->precision);
         if ( !strcasecmp(parcours->var->precision,"T") )
         {
            strcpy(parcours->var->precision,"");
         }
         fscanf(dependfileoutput,"%[^\n] \n",parcours->var->initialvalue);
         if ( !strcasecmp(parcours->var->initialvalue,"T") )
         {
            strcpy(parcours->var->initialvalue,"");
         }
         fscanf(dependfileoutput,"%d\n",&parcours->var->indicetabvars);
         fscanf(dependfileoutput,"%s\n",parcours->var->readedlistdimension);
         if ( !strcasecmp(parcours->var->readedlistdimension,"T") )
         {
            strcpy(parcours->var->readedlistdimension,"");
         }
         fscanf(dependfileoutput,"%s\n",nothing);
         parcours->suiv = NULL; 
         if ( !tmplocallist )
         { 
            tmplocallist = parcours;
            parcoursprec = parcours;
         }
         else
         {
            parcoursprec->suiv = parcours;
            parcoursprec = parcours;
         }
         parcours = NULL;
      }
      fclose(dependfileoutput);
  }
}


/******************************************************************************/
/*                     Recordglobalvarofusefile                              */
/******************************************************************************/
/* This subroutine is used to read the .dependfile<name> and to insert new    */
/*    information in the listout list.                                        */
/******************************************************************************/
/*                                                                            */
/*           .dependmodule -------->  globalvarofusefile = list of var        */
/*                                                                            */
/*        not.dependmodule -------->                                          */
/*                                                                            */
/******************************************************************************/
listvar *Recordglobalvarofusefile( char *NameTampon , listvar *listout)
{
  char ligne[LONGNOM];
  FILE *dependfileoutput;
  listvar *parcours0;
  listvar *parcours;
  listvar *parcoursprec;
  char nothing[LONGNOM];

  parcoursprec = (listvar *)NULL;
  /* we should free the listvar globalvarofusefile                            */
  sprintf(ligne,".depend%s",NameTampon);
  if ((dependfileoutput = fopen(ligne,"r"))==NULL) 
  {
    /* if the file doesn't exist it means that it is the first time           */
    /*    we tried to parse this file                                         */
  }
  else
  {
    /* if the file exist we should verify that this file has changed          */
      while (!feof(dependfileoutput))
      {
         parcours=(listvar *)malloc(sizeof(listvar));
         parcours->var=(variable *)malloc(sizeof(variable));
         fscanf(dependfileoutput,"%s\n",parcours->var->typevar);
         fscanf(dependfileoutput,"%s\n",parcours->var->nomvar);
         fscanf(dependfileoutput,"%s\n",parcours->var->dimchar);
         if ( !strcasecmp(parcours->var->dimchar,"T") )
         { 
            strcpy(parcours->var->dimchar,"");
         }
         fscanf(dependfileoutput,"%s\n",parcours->var->commoninfile);
         if ( !strcasecmp(parcours->var->commoninfile,"T") )
         { 
            strcpy(parcours->var->commoninfile,"");
         }
         fscanf(dependfileoutput,"%s\n",parcours->var->commonname);
         if ( !strcasecmp(parcours->var->commonname,"T") )
         { 
            strcpy(parcours->var->commonname,"");
         }
         fscanf(dependfileoutput,"%d\n",&parcours->var->nbdim);
         fscanf(dependfileoutput,"%d\n",&parcours->var->dimensiongiven);
         fscanf(dependfileoutput,"%d\n",&parcours->var->typegiven);
         fscanf(dependfileoutput,"%d\n",&parcours->var->allocatable);
         fscanf(dependfileoutput,"%d\n",&parcours->var->pointerdeclare);
         fscanf(dependfileoutput,"%[^\n] \n",parcours->var->precision);
         if ( !strcasecmp(parcours->var->precision,"T") )
         {
            strcpy(parcours->var->precision,"");
         }
         fscanf(dependfileoutput,"%[^\n] \n",parcours->var->initialvalue);
         if ( !strcasecmp(parcours->var->initialvalue,"T") )
         {
            strcpy(parcours->var->initialvalue,"");
         }
         fscanf(dependfileoutput,"%d\n",&parcours->var->indicetabvars);
         fscanf(dependfileoutput,"%s\n",parcours->var->readedlistdimension);
         if ( !strcasecmp(parcours->var->readedlistdimension,"T") )
         {
            strcpy(parcours->var->readedlistdimension,"");
         }
         fscanf(dependfileoutput,"%s\n",nothing);
         parcours->suiv = NULL; 
         if ( !listout )
         { 
            listout = parcours;
            parcoursprec = parcours;
         }
         else
         {
            if ( parcoursprec )
            {
               parcoursprec->suiv = parcours;
               parcoursprec = parcours;
            }
            else
            {
               parcours0 = listout;
               while ( parcours0->suiv ) parcours0=parcours0->suiv;
               parcours0->suiv = parcours;
               parcoursprec = parcours0->suiv;
            }
         }
         parcours = NULL;
      }
      fclose(dependfileoutput);
  }
  return listout;
}

/******************************************************************************/
/*                        Writethedependavailablefile                         */
/******************************************************************************/
/* This subroutine is used to write the .dependfileavailable file             */
/******************************************************************************/
/*                                                                            */
/*                                  .dependavailable                          */
/*     tabvars(1) = var1                                                      */
/*     tabvars(3) = var1                  2                                   */
/*     tabvars(4) = var1         =====>   5                                   */
/*     tabvars(6) = var1                                                      */
/*     tabvars(7) = var1                                                      */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void Writethedependavailablefile()
{
  FILE *dependfileoutput;
  listindice *parcours;

  if ( Listofavailableindices )
  {
  if ((dependfileoutput=fopen(".dependavailable","w"))!=NULL) 
  {
     /* We are looking for all the indices of the Listofavailableindices      */
     parcours = Listofavailableindices;
     while (parcours)
     {
        if ( parcours->indice != 0 )
        {
           fprintf(dependfileoutput,"%d\n",parcours->indice);
        }
        parcours = parcours->suiv;
     }
     fclose(dependfileoutput);
  }
  }
}

/******************************************************************************/
/*                        Readthedependavailablefile                          */
/******************************************************************************/
/* This subroutine is used to read the .dependfileavailable file              */
/******************************************************************************/
/*                                                                            */
/*                                  .dependavailable                          */
/*     tabvars(1) = var1                                                      */
/*     tabvars(3) = var1                  2                                   */
/*     tabvars(4) = var1         =====>   5  ==> Listofavailableindices       */
/*     tabvars(6) = var1                                                      */
/*     tabvars(7) = var1                                                      */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void Readthedependavailablefile()
{
  FILE *dependfileoutput;
  listindice *parcours;

  if ((dependfileoutput=fopen(".dependavailable","r"))!=NULL) 
  {
     /* We are looking for all the indices of the Listofavailableindices      */
     Listofavailableindices = (listindice *)NULL;
     while (!feof(dependfileoutput))
     {
        parcours=(listindice *)malloc(sizeof(listindice));
        fscanf(dependfileoutput,"%d\n",&parcours->indice);
        if ( parcours->indice != 0 && parcours->indice < 10000000 )
        {
           parcours -> suiv = Listofavailableindices;
           Listofavailableindices = parcours;
        }
        else
        {
           free(parcours);
        }
     }
     fclose(dependfileoutput);
  }
}


/******************************************************************************/
/*                      Did_filetoparse_readed                                */
/******************************************************************************/
/* This subroutine is used to know if the .depend<NameTampon> exist           */
/*    it means if the file has been ever parsed                               */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int Did_filetoparse_readed(char *NameTampon)
{
  FILE *dependfileoutput;
  char ligne[LONGNOM];
  int out;

  sprintf(ligne,".depend%s",NameTampon);
  if ((dependfileoutput = fopen(ligne,"r"))==NULL) 
  {
      out = 0;
  }
  else
  {
      out = 1;
      fclose(dependfileoutput);
  }
  return out;
}
