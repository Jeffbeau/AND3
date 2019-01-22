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
%{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "decl.h"
%}

%union {
       int ival;
       char na[LONGNOM];
       listnom * ln;
       }

%token TOK_SEP
%token TOK_USE
%token TOK_MINWIDTH        /* minimum width for the adaptive refinement       */
%token TOK_REGRIDDING
%token TOK_COEFFRAFX       /* space refinement factor in the x direction      */
%token TOK_COEFFRAFY       /* space refinement factor in the y direction      */
%token TOK_COEFFRAFZ       /* space refinement factor in the z direction      */
%token TOK_COEFFRAFTX      /* time refinement factor in the x direction       */
%token TOK_COEFFRAFTY      /* time refinement factor in the y direction       */
%token TOK_COEFFRAFTZ      /* time refinement factor in the z direction       */
%token TOK_MODULEMAIN      /* name of the module                              */
%token TOK_NOTGRIDDEP      /* Variable which are not grid dependent           */
%token <ival> TOK_NUM
%token <na> TOK_FILENAME
%token <na> TOK_USEITEM
%token <na> TOK_NAME
%token <na> TOK_PROBTYPE   /* dimension of the problem                        */
%token <na> TOK_EFFICIENCY /* efficiency for the adaptive refinement          */
%token <na> TOK_RAFMAX
%token <na> TOK_RAFMAXX
%token <na> TOK_RAFMAXY
%token <na> TOK_RAFMAXZ
%token ','
%token ';'
%token ':'
%token '('
%token ')'   
%token '['
%token ']'
%%
input :
      | input line
;
line :'\n'
      | TOK_PROBTYPE TOK_NAME ';'                  {initdimprob(1,$2,"0","0");}
      | TOK_PROBTYPE TOK_NAME ',' TOK_NAME ';'     {initdimprob(2,$2, $4,"0");} 
      | TOK_PROBTYPE TOK_NAME ',' TOK_NAME ',' TOK_NAME ';'
                                                   {initdimprob(3,$2, $4, $6);}
      | TOK_COEFFRAFX TOK_NUM ';'                              {coeffrafx =$2;}
      | TOK_REGRIDDING TOK_NUM ';'                             {regridding=$2;}
      | TOK_COEFFRAFY TOK_NUM ';'                              {coeffrafy =$2;}
      | TOK_COEFFRAFZ TOK_NUM ';'                              {coeffrafz =$2;}
      | TOK_COEFFRAFTX TOK_NUM ';'                             {coeffraftx=$2;}
      | TOK_COEFFRAFTY TOK_NUM ';'                             {coeffrafty=$2;}
      | TOK_COEFFRAFTZ TOK_NUM ';'                             {coeffraftz=$2;}
      | TOK_MODULEMAIN TOK_NAME ';'            {Add_ModuleTo_listofmodules($2);
                                                       Addmoduletothelist($2);}
      | TOK_NOTGRIDDEP TOK_SEP TOK_NAME ';'             {ajoutenotgriddep($3);}
      | TOK_EFFICIENCY TOK_NUM ';'            {userefficiency=1;efficiency=$2;}
      | TOK_RAFMAX TOK_NUM ';'              {rafmaxx=$2,rafmaxy=$2,rafmaxz=$2;}
      | TOK_RAFMAXX TOK_NUM ';'                                {rafmaxx  =$2;}
      | TOK_RAFMAXY TOK_NUM ';'                                {rafmaxy  =$2;}
      | TOK_RAFMAXZ TOK_NUM ';'                                {rafmaxz  =$2;}
      | TOK_MINWIDTH TOK_NUM ';'                               {minwidth =$2;}
      | TOK_USE TOK_USEITEM ';'  {
                                    if (!strcasecmp($2,"FIXED_GRIDS"))
                                                                 fixedgrids=1;
                                    if (!strcasecmp($2,"ONLY_FIXED_GRIDS"))
                                                             onlyfixedgrids=1;
                                    if (!strcasecmp($2,"DEBUG"))    todebug=1;
                                 }
      ;
%%

int main(int argc,char *argv[])
{
   extern FILE * yyin ;
   FILE *dependglobaloutput;
   char *tmp;
   int i;
   listnom *parcours;

/******************************************************************************/
/*  1-  Variables initialization                                              */
/******************************************************************************/
   globliste=(listvar *)NULL;
   listenamelist=(listnamelist *)NULL;
   globparam=(listvar *)NULL;
   globvarforcommon=(listvar *)NULL;
   AllocateList=(listallocate *)NULL;
   commonlist=(listvarcommon *)NULL;
   listofsubroutinewhereagrifisused=(listnom *)NULL;
   listofincludebysubroutine=(listusemodule *)NULL;
   listofmodulebysubroutine=(listusemodule *)NULL;
   listofmoduletmp=(listusemodule *)NULL;
   listmoduleinfile=(listmodule *)NULL;
   varofsubroutineliste=(listvar *)NULL;
   varsubroutine=(listvar *)NULL;
   listvarindoloop=(listvar *)NULL;
   listenotgriddepend=(listvar *)NULL;
   Listofavailableindices=(listindice *)NULL;
   Listofvarpointtovar=(listvarpointtovar *)NULL;
   globalvarofusefile = (listvar *)NULL;
   globalvarofusefile2 = (listvar *)NULL;
   
   strcpy(mainfile,argv[1]);    
   strcpy(nomdir,"AGRIF_INC");
   strcpy(commondirin,".");
   strcpy(commondirout,".");
   strcpy(filetoparse," "); 
   strcpy(subofagrifinitgrids,""); 
   strcpy(meetagrifinitgrids,"");
   strcpy(meetmpiinit,"");
   strcpy(mpiinitvar,"");

   listofvarofusemodulecreated=0;
   checkexistcommon=1;
   fortran77 = 0 ;
   Did_filetoparse_treated = 0 ;
   userefficiency=0;
   todebug=0;
   onlyfixedgrids=0;
   fixedgrids=0;
   InAgrifParentDef = 0;
   rafmaxx=0;
   rafmaxy=0;
   rafmaxz=0;
   IntegerIShouldBeAdd=0;
   IndicenbmaillesX=0;
   IndicenbmaillesY=0;
   IndicenbmaillesZ=0;
   minwidth=1;
   coeffrafx=0;
   coeffrafy=0;
   coeffrafz=0;
   indicemaxtabvars = 0;   /* current indice in the table tabvars             */
   oldindicemaxtabvars = 0;/* current indice in the table tabvars             */
   SubloopScalar = 0;
/******************************************************************************/
/*  2-  Program arguments                                                     */
/******************************************************************************/
   if (argc < 2) 
   {
       printf("usage : conv <file> [-rm] [-incdir <directory>] \n");
       printf(" [-comdirin   <directory>] [-comdirout <directory>]\n");
       printf(" [-convfile  <FILENAME >] -SubloopScalar -SubloopScalar1 \n"); 
       printf(" -f77\n"); 
       exit(0);
   }


   if ((yyin=fopen(argv[1],"r"))==NULL) 
   {
      printf("the file %s doesn't exist \n",argv[1]);
      exit(0);    
   }

   i=2;
   while (i<argc)
   {
      if (!strcasecmp(argv[i],"-incdir")) 
      {
         strcpy(nomdir,argv[i+1]);
         i++;
      }
      else if (!strcasecmp(argv[i],"-comdirin")) /* input directory           */
      {     
         strcpy(commondirin,argv[i+1]);
         i++;
      }
      else if (!strcasecmp(argv[i],"-comdirout")) /* output directory         */
      {
         strcpy(commondirout,argv[i+1]);
         i++;
      }      
      else if (!strcasecmp(argv[i],"-convfile")) /* file to parse             */
      {     
         strcpy(filetoparse,argv[i+1]);
         i++;
      }   
      else if (!strcasecmp(argv[i],"-f77")) /* fortran 77 file to parse       */
      {     
         fortran77 = 1 ;
      }   
      else if (!strcasecmp(argv[i],"-SubloopScalar")) /* file to parse        */
      {     
         SubloopScalar = 1 ;
      }   
      else if (!strcasecmp(argv[i],"-SubloopScalar1")) /* file to parse       */
      {     
         SubloopScalar = 2 ;
      }   
      else if (!strcasecmp(argv[i],"-rm")) 
      {     
         checkexistcommon=0;
      }
      else 
      {
         printf("Unkwon option : %s\n",argv[i]);
         exit(0);
      }
      i++;       
   }  

/******************************************************************************/
/*  3-  Parsing of the  conv file <name>.in                                   */
/******************************************************************************/

   if ((yyin=fopen(argv[1],"r"))==NULL) 
   {
       printf("the file %s doesn't exist \n",argv[1]);
       exit(0);    
   }
   strcpy(mainfile,argv[1]);    

   yyparse();

/******************************************************************************/
/*  4-  Preparation of the file parsing                                       */
/******************************************************************************/
   if ((yyin=fopen(filetoparse,"r"))==NULL) /* Is the file to parse exist ?   */
   {
      printf("the file %s doesn't exist \n",filetoparse);
      exit(0);    
   }
   /* NameTamponfile : the name of the model file extract from the name       */
   /*    of agrif_module_<NameTamponfile>                                     */
   tmp = strchr(filetoparse, '.');
   NameTamponfile=(char *)malloc(
                              (strlen(filetoparse)-strlen(tmp)+1)*sizeof(char));
   strncpy(NameTamponfile,filetoparse,strlen(filetoparse)-strlen(tmp)+1);
   strcpy (&NameTamponfile[strlen(filetoparse)-strlen(tmp)], "\0");
   /* mainfile : the name of the file to parse                                */
   strcpy(mainfile,filetoparse);    
   /* We should verify that this file has not been read before                */
   /* if it is the case we record the old globliste in the tmplocallist       */
   tmplocallist = NULL;
   tmpuselocallist = NULL;
   Did_filetoparse_treated = Did_filetoparse_readed(NameTamponfile);
   /* if  Did_filetoparse_treated = 1 then the file to parse has been treated */
   if ( Did_filetoparse_treated == 0 ) 
   {
     /* if the filetoparse has not been treated, we should know the last      */
     /*    tabvars indices which has been used                                */
     if ((dependglobaloutput=fopen(".dependglobal","r"))!=NULL) 
     {
        fscanf(dependglobaloutput,"%d\n",&indicemaxtabvars);
        fclose(dependglobaloutput);
        oldindicemaxtabvars = indicemaxtabvars;
     }
   }   
   /* Write the .dependnbxnby file which contains indices of nbmaillsX,       */
   /*    nbmailleY and nbmailleZ                                              */
   Readthedependnbxnbyfile();

/******************************************************************************/
/*  4-  Parsing of the input file (2 times)                                   */
/******************************************************************************/

   firstpass = 1; 
   processfortran(filetoparse); 
   firstpass = 0; 
   processfortran(filetoparse);

/******************************************************************************/
/*  5-  Write informations in output files                                    */
/******************************************************************************/

   if ( Did_filetoparse_treated == 0 ) /* if the file has never been treated  */
   {
      /* Write the .dependglobal file which contain the max indice            */
      /*    of the tabvars table                                              */
      dependglobaloutput = fopen(".dependglobal","w");
      fprintf(dependglobaloutput,"%d\n",indicemaxtabvars);
      fclose(dependglobaloutput);
      /* Write the .depend<namefile> file which contain general informations  */
      /*    about variable of this file                                       */
      parcours = modulelist;
      while( parcours )
      {
         Writethedependfile(parcours->nom,globliste);
         parcours=parcours->suiv;
      }
   }

/******************************************************************************/
/*  7-  Remove the non grid dependent variables                               */
/******************************************************************************/

   /* we should remove from the globliste the non grid dependent variables    */
   RemoveNotgriddependFromGlobliste();

/******************************************************************************/
/*  8-  Write informations in output files                                    */
/******************************************************************************/

   /* if this file has been treated in past called,                           */
   /*    we should compare the old parsing (record in the tmplocallist)       */
   /*    and the new one contained in the globliste                           */
   if ( Did_filetoparse_treated == 1 ) 
   {
      parcours = modulelist;
      while( parcours )
      {
         Recordtmplocallist(parcours->nom);
         parcours=parcours->suiv;
      }
      /* if the filetoparse has not been treated, we should know              */
      /*    the last tabvars indices which has been used                      */
     if ((dependglobaloutput=fopen(".dependglobal","r"))!=NULL) 
     {
        fscanf(dependglobaloutput,"%d\n",&indicemaxtabvars);
        fclose(dependglobaloutput);
        oldindicemaxtabvars = indicemaxtabvars;
     }
     /* Read the list of available indice                                     */
     Readthedependavailablefile();
     /* the old treatement has been recorded in the tmplocallist              */
     /* Now we should compare the old treatement with the new one             */
/*mazauric for each module */
     CompareNewparsingandoldone();
     /* Write the .dependglobal file which contain general informations       */
     /*    about globlist                                                     */
     dependglobaloutput = fopen(".dependglobal","w");
     fprintf(dependglobaloutput,"%d\n",indicemaxtabvars);
     fclose(dependglobaloutput);
     /* Write the list of available indice                                    */
     Writethedependavailablefile();  
     /* Write the .depend<namefile> file which contain general                */
     /*    informations about variable of this file                           */
     parcours = modulelist;
     while( parcours )
     {
        Writethedependfile(parcours->nom,globliste);
        parcours=parcours->suiv;
     }
     /* Write the .dependnbxnby file which contains indices of nbmaillsX,     */
     /*    nbmailleY and nbmailleZ                                            */
     Writethedependnbxnbyfile();
   }
   /* Write the .dependnbxnby file which contains indices of nbmaillsX,       */
   /*    nbmailleY and nbmailleZ                                              */
   Writethedependnbxnbyfile();
   /* Write the list of module used in this file                              */
   Writethedependlistofmoduleused(NameTamponfile);
   WritedependParameterList(NameTamponfile);
/******************************************************************************/
/*  8-  Create files in AGRIF_INC directory                                   */
/******************************************************************************/
   creefichieramr(NameTamponfile);
   return 0;
}
