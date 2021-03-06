#ifndef lint
static char yysccsid[] = "@(#)yaccpar	1.9 (Berkeley) 02/21/93";
#endif
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define yyclearin (yychar=(-1))
#define yyerrok (yyerrflag=0)
#define YYRECOVERING (yyerrflag!=0)
#define YYPREFIX "yy"
#line 18 "convert.y"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "decl.h"
#line 24 "convert.y"
typedef union {
       int ival;
       char na[LONGNOM];
       listnom * ln;
       } YYSTYPE;
#line 23 "y.tab.c"
#define TOK_SEP 257
#define TOK_USE 258
#define TOK_MINWIDTH 259
#define TOK_REGRIDDING 260
#define TOK_COEFFRAFX 261
#define TOK_COEFFRAFY 262
#define TOK_COEFFRAFZ 263
#define TOK_COEFFRAFTX 264
#define TOK_COEFFRAFTY 265
#define TOK_COEFFRAFTZ 266
#define TOK_MODULEMAIN 267
#define TOK_NOTGRIDDEP 268
#define TOK_NUM 269
#define TOK_FILENAME 270
#define TOK_USEITEM 271
#define TOK_NAME 272
#define TOK_PROBTYPE 273
#define TOK_EFFICIENCY 274
#define TOK_RAFMAX 275
#define TOK_RAFMAXX 276
#define TOK_RAFMAXY 277
#define TOK_RAFMAXZ 278
#define YYERRCODE 256
short yylhs[] = {                                        -1,
    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,
};
short yylen[] = {                                         2,
    0,    2,    1,    3,    5,    7,    3,    3,    3,    3,
    3,    3,    3,    3,    4,    3,    3,    3,    3,    3,
    3,    3,
};
short yydefred[] = {                                      1,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    3,    2,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,   22,   21,    8,
    7,    9,   10,   11,   12,   13,   14,    0,    0,    4,
   16,   17,   18,   19,   20,   15,    0,    0,    5,    0,
    6,
};
short yydgoto[] = {                                       1,
   20,
};
short yysindex[] = {                                      0,
  -10, -268, -265, -264, -263, -262, -261, -260, -259, -258,
 -257, -245, -254, -256, -255, -250, -249, -248,    0,    0,
  -37,  -36,  -35,  -34,  -33,  -32,  -31,  -30,  -29,  -28,
 -240,  -43,  -26,  -25,  -24,  -23,  -22,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,  -21, -233,    0,
    0,    0,    0,    0,    0,    0,  -42, -232,    0,  -18,
    0,
};
short yyrindex[] = {                                      0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,
};
short yygindex[] = {                                      0,
    0,
};
#define YYTABLESIZE 268
short yytable[] = {                                      19,
   49,   58,   21,   22,   23,   24,   25,   26,   27,   28,
   29,   31,   33,   34,   30,   50,   59,   32,   35,   36,
   37,   38,   39,   40,   41,   42,   43,   44,   45,   46,
   47,   48,   51,   52,   53,   54,   55,   56,   57,   60,
   61,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    2,    3,    4,
    5,    6,    7,    8,    9,   10,   11,   12,    0,    0,
    0,    0,   13,   14,   15,   16,   17,   18,
};
short yycheck[] = {                                      10,
   44,   44,  271,  269,  269,  269,  269,  269,  269,  269,
  269,  257,  269,  269,  272,   59,   59,  272,  269,  269,
  269,   59,   59,   59,   59,   59,   59,   59,   59,   59,
   59,  272,   59,   59,   59,   59,   59,   59,  272,  272,
   59,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,  258,  259,  260,
  261,  262,  263,  264,  265,  266,  267,  268,   -1,   -1,
   -1,   -1,  273,  274,  275,  276,  277,  278,
};
#define YYFINAL 1
#ifndef YYDEBUG
#define YYDEBUG 1
#endif
#define YYMAXTOKEN 278
#if YYDEBUG
char *yyname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,"'\\n'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,"'('","')'",0,0,"','",0,0,0,0,0,0,0,0,0,0,0,0,0,"':'","';'",0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"'['",0,"']'",0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
"TOK_SEP","TOK_USE","TOK_MINWIDTH","TOK_REGRIDDING","TOK_COEFFRAFX",
"TOK_COEFFRAFY","TOK_COEFFRAFZ","TOK_COEFFRAFTX","TOK_COEFFRAFTY",
"TOK_COEFFRAFTZ","TOK_MODULEMAIN","TOK_NOTGRIDDEP","TOK_NUM","TOK_FILENAME",
"TOK_USEITEM","TOK_NAME","TOK_PROBTYPE","TOK_EFFICIENCY","TOK_RAFMAX",
"TOK_RAFMAXX","TOK_RAFMAXY","TOK_RAFMAXZ",
};
char *yyrule[] = {
"$accept : input",
"input :",
"input : input line",
"line : '\\n'",
"line : TOK_PROBTYPE TOK_NAME ';'",
"line : TOK_PROBTYPE TOK_NAME ',' TOK_NAME ';'",
"line : TOK_PROBTYPE TOK_NAME ',' TOK_NAME ',' TOK_NAME ';'",
"line : TOK_COEFFRAFX TOK_NUM ';'",
"line : TOK_REGRIDDING TOK_NUM ';'",
"line : TOK_COEFFRAFY TOK_NUM ';'",
"line : TOK_COEFFRAFZ TOK_NUM ';'",
"line : TOK_COEFFRAFTX TOK_NUM ';'",
"line : TOK_COEFFRAFTY TOK_NUM ';'",
"line : TOK_COEFFRAFTZ TOK_NUM ';'",
"line : TOK_MODULEMAIN TOK_NAME ';'",
"line : TOK_NOTGRIDDEP TOK_SEP TOK_NAME ';'",
"line : TOK_EFFICIENCY TOK_NUM ';'",
"line : TOK_RAFMAX TOK_NUM ';'",
"line : TOK_RAFMAXX TOK_NUM ';'",
"line : TOK_RAFMAXY TOK_NUM ';'",
"line : TOK_RAFMAXZ TOK_NUM ';'",
"line : TOK_MINWIDTH TOK_NUM ';'",
"line : TOK_USE TOK_USEITEM ';'",
};
#endif
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 500
#define YYMAXDEPTH 500
#endif
#endif
int yydebug;
int yynerrs;
int yyerrflag;
int yychar;
short *yyssp;
YYSTYPE *yyvsp;
YYSTYPE yyval;
YYSTYPE yylval;
short yyss[YYSTACKSIZE];
YYSTYPE yyvs[YYSTACKSIZE];
#define yystacksize YYSTACKSIZE
#line 93 "convert.y"

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
#line 497 "y.tab.c"
#define YYABORT goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR goto yyerrlab
int
yyparse()
{
    register int yym, yyn, yystate;
#if YYDEBUG
    register char *yys;
    extern char *getenv();

    if (yys = getenv("YYDEBUG"))
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = (-1);

    yyssp = yyss;
    yyvsp = yyvs;
    *yyssp = yystate = 0;

yyloop:
    if (yyn = yydefred[yystate]) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yyssp >= yyss + yystacksize - 1)
        {
            goto yyoverflow;
        }
        *++yyssp = yystate = yytable[yyn];
        *++yyvsp = yylval;
        yychar = (-1);
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;
#ifdef lint
    goto yynewerror;
#endif
yynewerror:
    yyerror("syntax error");
#ifdef lint
    goto yyerrlab;
#endif
yyerrlab:
    ++yynerrs;
yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yyssp]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yyssp, yytable[yyn]);
#endif
                if (yyssp >= yyss + yystacksize - 1)
                {
                    goto yyoverflow;
                }
                *++yyssp = yystate = yytable[yyn];
                *++yyvsp = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yyssp);
#endif
                if (yyssp <= yyss) goto yyabort;
                --yyssp;
                --yyvsp;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = (-1);
        goto yyloop;
    }
yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    yyval = yyvsp[1-yym];
    switch (yyn)
    {
case 4:
#line 64 "convert.y"
{initdimprob(1,yyvsp[-1].na,"0","0");}
break;
case 5:
#line 65 "convert.y"
{initdimprob(2,yyvsp[-3].na, yyvsp[-1].na,"0");}
break;
case 6:
#line 67 "convert.y"
{initdimprob(3,yyvsp[-5].na, yyvsp[-3].na, yyvsp[-1].na);}
break;
case 7:
#line 68 "convert.y"
{coeffrafx =yyvsp[-1].ival;}
break;
case 8:
#line 69 "convert.y"
{regridding=yyvsp[-1].ival;}
break;
case 9:
#line 70 "convert.y"
{coeffrafy =yyvsp[-1].ival;}
break;
case 10:
#line 71 "convert.y"
{coeffrafz =yyvsp[-1].ival;}
break;
case 11:
#line 72 "convert.y"
{coeffraftx=yyvsp[-1].ival;}
break;
case 12:
#line 73 "convert.y"
{coeffrafty=yyvsp[-1].ival;}
break;
case 13:
#line 74 "convert.y"
{coeffraftz=yyvsp[-1].ival;}
break;
case 14:
#line 75 "convert.y"
{Add_ModuleTo_listofmodules(yyvsp[-1].na);
                                                       Addmoduletothelist(yyvsp[-1].na);}
break;
case 15:
#line 77 "convert.y"
{ajoutenotgriddep(yyvsp[-1].na);}
break;
case 16:
#line 78 "convert.y"
{userefficiency=1;efficiency=yyvsp[-1].ival;}
break;
case 17:
#line 79 "convert.y"
{rafmaxx=yyvsp[-1].ival,rafmaxy=yyvsp[-1].ival,rafmaxz=yyvsp[-1].ival;}
break;
case 18:
#line 80 "convert.y"
{rafmaxx  =yyvsp[-1].ival;}
break;
case 19:
#line 81 "convert.y"
{rafmaxy  =yyvsp[-1].ival;}
break;
case 20:
#line 82 "convert.y"
{rafmaxz  =yyvsp[-1].ival;}
break;
case 21:
#line 83 "convert.y"
{minwidth =yyvsp[-1].ival;}
break;
case 22:
#line 84 "convert.y"
{
                                    if (!strcasecmp(yyvsp[-1].na,"FIXED_GRIDS"))
                                                                 fixedgrids=1;
                                    if (!strcasecmp(yyvsp[-1].na,"ONLY_FIXED_GRIDS"))
                                                             onlyfixedgrids=1;
                                    if (!strcasecmp(yyvsp[-1].na,"DEBUG"))    todebug=1;
                                 }
break;
#line 721 "y.tab.c"
    }
    yyssp -= yym;
    yystate = *yyssp;
    yyvsp -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yyssp = YYFINAL;
        *++yyvsp = yyval;
        if (yychar < 0)
        {
            if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yyssp, yystate);
#endif
    if (yyssp >= yyss + yystacksize - 1)
    {
        goto yyoverflow;
    }
    *++yyssp = yystate;
    *++yyvsp = yyval;
    goto yyloop;
yyoverflow:
    yyerror("yacc stack overflow");
yyabort:
    return (1);
yyaccept:
    return (0);
}
#line 2 "convert.yy.c"
/* A lexical scanner generated by flex */

/* Scanner skeleton version:
 * $Header: /home/opalod/NEMOCVSROOT/AGRIF/LIB/main.c,v 1.1.1.1 2006/03/10 17:58:35 opalod Exp $
 */

#define FLEX_SCANNER
#define YY_FLEX_MAJOR_VERSION 2
#define YY_FLEX_MINOR_VERSION 5

#include <stdio.h>
#include <unistd.h>


/* cfront 1.2 defines "c_plusplus" instead of "__cplusplus" */
#ifdef c_plusplus
#ifndef __cplusplus
#define __cplusplus
#endif
#endif


#ifdef __cplusplus

#include <stdlib.h>

/* Use prototypes in function declarations. */
#define YY_USE_PROTOS

/* The "const" storage-class-modifier is valid. */
#define YY_USE_CONST

#else	/* ! __cplusplus */

#if __STDC__

#define YY_USE_PROTOS
#define YY_USE_CONST

#endif	/* __STDC__ */
#endif	/* ! __cplusplus */

#ifdef __TURBOC__
 #pragma warn -rch
 #pragma warn -use
#include <io.h>
#include <stdlib.h>
#define YY_USE_CONST
#define YY_USE_PROTOS
#endif

#ifdef YY_USE_CONST
#define yyconst const
#else
#define yyconst
#endif


#ifdef YY_USE_PROTOS
#define YY_PROTO(proto) proto
#else
#define YY_PROTO(proto) ()
#endif

/* Returned upon end-of-file. */
#define YY_NULL 0

/* Promotes a possibly negative, possibly signed char to an unsigned
 * integer for use as an array index.  If the signed char is negative,
 * we want to instead treat it as an 8-bit unsigned char, hence the
 * double cast.
 */
#define YY_SC_TO_UI(c) ((unsigned int) (unsigned char) c)

/* Enter a start condition.  This macro really ought to take a parameter,
 * but we do it the disgusting crufty way forced on us by the ()-less
 * definition of BEGIN.
 */
#define BEGIN yy_start = 1 + 2 *

/* Translate the current start state into a value that can be later handed
 * to BEGIN to return to the state.  The YYSTATE alias is for lex
 * compatibility.
 */
#define YY_START ((yy_start - 1) / 2)
#define YYSTATE YY_START

/* Action number for EOF rule of a given start state. */
#define YY_STATE_EOF(state) (YY_END_OF_BUFFER + state + 1)

/* Special action meaning "start processing a new file". */
#define YY_NEW_FILE yyrestart( yyin )

#define YY_END_OF_BUFFER_CHAR 0

/* Size of default input buffer. */
#define YY_BUF_SIZE 16384

typedef struct yy_buffer_state *YY_BUFFER_STATE;

extern int yyleng;
extern FILE *yyin, *yyout;

#define EOB_ACT_CONTINUE_SCAN 0
#define EOB_ACT_END_OF_FILE 1
#define EOB_ACT_LAST_MATCH 2

/* The funky do-while in the following #define is used to turn the definition
 * int a single C statement (which needs a semi-colon terminator).  This
 * avoids problems with code like:
 *
 * 	if ( condition_holds )
 *		yyless( 5 );
 *	else
 *		do_something_else();
 *
 * Prior to using the do-while the compiler would get upset at the
 * "else" because it interpreted the "if" statement as being all
 * done when it reached the ';' after the yyless() call.
 */

/* Return all but the first 'n' matched characters back to the input stream. */

#define yyless(n) \
	do \
		{ \
		/* Undo effects of setting up yytext. */ \
		*yy_cp = yy_hold_char; \
		YY_RESTORE_YY_MORE_OFFSET \
		yy_c_buf_p = yy_cp = yy_bp + n - YY_MORE_ADJ; \
		YY_DO_BEFORE_ACTION; /* set up yytext again */ \
		} \
	while ( 0 )

#define unput(c) yyunput( c, yytext_ptr )

/* The following is because we cannot portably get our hands on size_t
 * (without autoconf's help, which isn't available because we want
 * flex-generated scanners to compile on their own).
 */
typedef unsigned int yy_size_t;


struct yy_buffer_state
	{
	FILE *yy_input_file;

	char *yy_ch_buf;		/* input buffer */
	char *yy_buf_pos;		/* current position in input buffer */

	/* Size of input buffer in bytes, not including room for EOB
	 * characters.
	 */
	yy_size_t yy_buf_size;

	/* Number of characters read into yy_ch_buf, not including EOB
	 * characters.
	 */
	int yy_n_chars;

	/* Whether we "own" the buffer - i.e., we know we created it,
	 * and can realloc() it to grow it, and should free() it to
	 * delete it.
	 */
	int yy_is_our_buffer;

	/* Whether this is an "interactive" input source; if so, and
	 * if we're using stdio for input, then we want to use getc()
	 * instead of fread(), to make sure we stop fetching input after
	 * each newline.
	 */
	int yy_is_interactive;

	/* Whether we're considered to be at the beginning of a line.
	 * If so, '^' rules will be active on the next match, otherwise
	 * not.
	 */
	int yy_at_bol;

	/* Whether to try to fill the input buffer when we reach the
	 * end of it.
	 */
	int yy_fill_buffer;

	int yy_buffer_status;
#define YY_BUFFER_NEW 0
#define YY_BUFFER_NORMAL 1
	/* When an EOF's been seen but there's still some text to process
	 * then we mark the buffer as YY_EOF_PENDING, to indicate that we
	 * shouldn't try reading from the input source any more.  We might
	 * still have a bunch of tokens to match, though, because of
	 * possible backing-up.
	 *
	 * When we actually see the EOF, we change the status to "new"
	 * (via yyrestart()), so that the user can continue scanning by
	 * just pointing yyin at a new input file.
	 */
#define YY_BUFFER_EOF_PENDING 2
	};

static YY_BUFFER_STATE yy_current_buffer = 0;

/* We provide macros for accessing buffer states in case in the
 * future we want to put the buffer states in a more general
 * "scanner state".
 */
#define YY_CURRENT_BUFFER yy_current_buffer


/* yy_hold_char holds the character lost when yytext is formed. */
static char yy_hold_char;

static int yy_n_chars;		/* number of characters read into yy_ch_buf */


int yyleng;

/* Points to current character in buffer. */
static char *yy_c_buf_p = (char *) 0;
static int yy_init = 1;		/* whether we need to initialize */
static int yy_start = 0;	/* start state number */

/* Flag which is used to allow yywrap()'s to do buffer switches
 * instead of setting up a fresh yyin.  A bit of a hack ...
 */
static int yy_did_buffer_switch_on_eof;

void yyrestart YY_PROTO(( FILE *input_file ));

void yy_switch_to_buffer YY_PROTO(( YY_BUFFER_STATE new_buffer ));
void yy_load_buffer_state YY_PROTO(( void ));
YY_BUFFER_STATE yy_create_buffer YY_PROTO(( FILE *file, int size ));
void yy_delete_buffer YY_PROTO(( YY_BUFFER_STATE b ));
void yy_init_buffer YY_PROTO(( YY_BUFFER_STATE b, FILE *file ));
void yy_flush_buffer YY_PROTO(( YY_BUFFER_STATE b ));
#define YY_FLUSH_BUFFER yy_flush_buffer( yy_current_buffer )

YY_BUFFER_STATE yy_scan_buffer YY_PROTO(( char *base, yy_size_t size ));
YY_BUFFER_STATE yy_scan_string YY_PROTO(( yyconst char *yy_str ));
YY_BUFFER_STATE yy_scan_bytes YY_PROTO(( yyconst char *bytes, int len ));

static void *yy_flex_alloc YY_PROTO(( yy_size_t ));
static void *yy_flex_realloc YY_PROTO(( void *, yy_size_t ));
static void yy_flex_free YY_PROTO(( void * ));

#define yy_new_buffer yy_create_buffer

#define yy_set_interactive(is_interactive) \
	{ \
	if ( ! yy_current_buffer ) \
		yy_current_buffer = yy_create_buffer( yyin, YY_BUF_SIZE ); \
	yy_current_buffer->yy_is_interactive = is_interactive; \
	}

#define yy_set_bol(at_bol) \
	{ \
	if ( ! yy_current_buffer ) \
		yy_current_buffer = yy_create_buffer( yyin, YY_BUF_SIZE ); \
	yy_current_buffer->yy_at_bol = at_bol; \
	}

#define YY_AT_BOL() (yy_current_buffer->yy_at_bol)

typedef unsigned char YY_CHAR;
FILE *yyin = (FILE *) 0, *yyout = (FILE *) 0;
typedef int yy_state_type;
extern char *yytext;
#define yytext_ptr yytext

static yy_state_type yy_get_previous_state YY_PROTO(( void ));
static yy_state_type yy_try_NUL_trans YY_PROTO(( yy_state_type current_state ));
static int yy_get_next_buffer YY_PROTO(( void ));
static void yy_fatal_error YY_PROTO(( yyconst char msg[] ));

/* Done after the current pattern has been matched and before the
 * corresponding action - sets up yytext.
 */
#define YY_DO_BEFORE_ACTION \
	yytext_ptr = yy_bp; \
	yyleng = (int) (yy_cp - yy_bp); \
	yy_hold_char = *yy_cp; \
	*yy_cp = '\0'; \
	yy_c_buf_p = yy_cp;

#define YY_NUM_RULES 27
#define YY_END_OF_BUFFER 28
static yyconst short int yy_accept[133] =
    {   0,
        0,    0,    0,    0,   28,   27,   26,   25,   27,   24,
       23,   23,   23,   23,   24,   22,   22,   22,   22,   22,
       22,   22,   22,   22,   22,   22,   26,    0,   17,   23,
       21,   18,    0,   22,   22,   22,   22,   22,   22,   22,
       22,   22,   22,   22,   22,   19,   22,   22,   22,   22,
       22,   22,   22,   22,   22,   22,   15,   19,   22,   22,
       22,   22,   22,   22,   22,   22,   22,   22,   22,   20,
       22,   22,   22,   22,   22,   22,   22,   22,   22,   22,
       22,   22,   22,   22,   22,   10,   22,   22,   22,   22,
       22,   22,   22,   22,   11,   12,   13,   22,   22,   22,

       22,   16,   22,   22,   22,   22,   22,    2,    3,    4,
       22,   22,   22,   22,   22,   22,    5,    6,    7,    9,
       22,   14,   22,   22,    1,   22,    8,   22,   22,   22,
       22,    0
    } ;

static yyconst int yy_ec[256] =
    {   0,
        1,    1,    1,    1,    1,    1,    1,    1,    2,    3,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    2,    1,    1,    1,    1,    4,    1,    1,    5,
        6,    1,    1,    7,    1,    8,    1,    9,   10,   11,
       12,    9,    9,    9,    9,    9,    9,   13,   14,    1,
        1,    1,    1,    1,   18,   19,   20,   21,   22,   23,
       24,   25,   26,   27,   27,   28,   29,   30,   31,   32,
       27,   33,   34,   35,   36,   27,   37,   38,   39,   40,
       15,    1,   16,    1,   17,    1,   18,   19,   20,   21,

       22,   23,   24,   25,   26,   27,   27,   28,   29,   30,
       31,   32,   27,   33,   34,   35,   36,   27,   37,   38,
       39,   40,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,

        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1
    } ;

static yyconst int yy_meta[41] =
    {   0,
        1,    1,    2,    1,    1,    1,    1,    3,    4,    4,
        4,    4,    1,    1,    1,    1,    5,    5,    5,    5,
        5,    5,    5,    5,    5,    5,    5,    5,    5,    5,
        5,    5,    5,    5,    5,    5,    5,    5,    5,    5
    } ;

static yyconst short int yy_base[137] =
    {   0,
        0,    0,    0,    0,  284,  285,  281,  285,  278,  285,
       32,   36,   40,   44,  268,  272,   50,   51,   52,   54,
       56,   55,   58,   59,   61,   60,  277,  274,  273,   80,
      285,  285,    0,  268,   62,   66,   64,   63,   68,   70,
       85,   87,   88,   91,   92,    0,   89,   94,   95,   96,
       98,   99,  100,  101,  102,  108,  267,    0,  109,  114,
      116,  119,  117,  118,  120,  121,  126,  134,  125,  266,
      137,  138,  140,  139,  141,  144,  145,  146,  148,  149,
      151,  154,  158,  160,  161,  164,  166,  168,  169,  172,
      170,  173,  174,  176,  265,  264,  263,  180,  182,  188,

      185,  262,  192,  193,  177,  199,  201,  253,  250,  249,
      208,  202,  210,  211,  216,  219,  248,  247,  246,  242,
      217,  240,  218,  223,  230,  222,  229,  220,  226,  228,
      225,  285,  259,  262,  193,  264
    } ;

static yyconst short int yy_def[137] =
    {   0,
      132,    1,    1,    1,  132,  132,  132,  132,  133,  132,
      132,  132,  132,  132,  132,  134,  134,  134,  134,  134,
      134,  134,  134,  134,  134,  134,  132,  133,  133,  132,
      132,  132,  135,  134,  134,  134,  134,  134,  134,  134,
      134,  134,  134,  134,  134,  136,  134,  134,  134,  134,
      134,  134,  134,  134,  134,  134,  134,  136,  134,  134,
      134,  134,  134,  134,  134,  134,  134,  134,  134,  134,
      134,  134,  134,  134,  134,  134,  134,  134,  134,  134,
      134,  134,  134,  134,  134,  134,  134,  134,  134,  134,
      134,  134,  134,  134,  134,  134,  134,  134,  134,  134,

      134,  134,  134,  134,  134,  134,  134,  134,  134,  134,
      134,  134,  134,  134,  134,  134,  134,  134,  134,  134,
      134,  134,  134,  134,  134,  134,  134,  134,  134,  134,
      134,    0,  132,  132,  132,  132
    } ;

static yyconst short int yy_nxt[326] =
    {   0,
        6,    7,    8,    9,   10,   10,   10,    6,   11,   12,
       13,   14,   15,   10,   10,   10,   16,   16,   16,   17,
       18,   19,   20,   16,   16,   16,   16,   16,   21,   22,
       23,   24,   25,   16,   16,   26,   16,   16,   16,   16,
       30,   30,   30,   30,   30,   30,   30,   30,   30,   30,
       30,   30,   30,   30,   30,   30,   31,   33,   33,   33,
       31,   33,   33,   33,   31,   33,   33,   33,   33,   33,
       33,   33,   36,   33,   37,   33,   42,   33,   43,   38,
       35,   39,   44,   47,   48,   40,   49,   41,   30,   30,
       30,   30,   33,   45,   33,   33,   33,   51,   33,   33,

       50,   33,   33,   33,   52,   33,   33,   33,   33,   33,
       55,   59,   53,   57,   56,   33,   33,   62,   66,   54,
       61,   33,   64,   33,   33,   33,   33,   33,   33,   60,
       67,   69,   33,   33,   63,   71,   75,   70,   65,   72,
       68,   33,   73,   77,   33,   33,   33,   33,   33,   76,
       74,   33,   33,   33,   81,   33,   33,   79,   33,   78,
       82,   33,   80,   84,   83,   33,   87,   33,   33,   88,
       89,   33,   85,   33,   90,   33,   33,   33,   92,   33,
       33,   33,   86,   33,   33,   93,   98,   33,   91,   33,
       99,   94,   33,  103,  102,   33,  105,   46,  100,   33,

       33,   95,   96,   97,  101,  106,   33,  111,   33,   33,
      112,  104,  115,  113,  114,   33,  107,   33,   33,  108,
      109,  110,  121,   33,   33,   33,   33,   33,  116,   33,
       33,  123,   33,   33,  126,   33,   33,   33,  117,  118,
      119,  122,  125,  124,  127,  128,  120,   33,  131,   33,
       70,  130,  129,   33,   33,   33,   33,   33,   70,   28,
       33,   28,   28,   28,   34,   34,   34,   58,   58,   33,
       33,   33,   33,   33,   33,   33,   29,   29,   27,   33,
       32,   29,   27,  132,    5,  132,  132,  132,  132,  132,
      132,  132,  132,  132,  132,  132,  132,  132,  132,  132,

      132,  132,  132,  132,  132,  132,  132,  132,  132,  132,
      132,  132,  132,  132,  132,  132,  132,  132,  132,  132,
      132,  132,  132,  132,  132
    } ;

static yyconst short int yy_chk[326] =
    {   0,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
        1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
       11,   11,   11,   11,   12,   12,   12,   12,   13,   13,
       13,   13,   14,   14,   14,   14,   12,   17,   18,   19,
       13,   20,   22,   21,   14,   23,   24,   26,   25,   35,
       38,   37,   18,   36,   19,   39,   24,   40,   25,   20,
       17,   21,   25,   35,   36,   22,   37,   23,   30,   30,
       30,   30,   41,   26,   42,   43,   47,   39,   44,   45,

       38,   48,   49,   50,   40,   51,   52,   53,   54,   55,
       43,   47,   41,   45,   44,   56,   59,   50,   54,   42,
       49,   60,   52,   61,   63,   64,   62,   65,   66,   48,
       55,   59,   69,   67,   51,   61,   65,   60,   53,   62,
       56,   68,   63,   67,   71,   72,   74,   73,   75,   66,
       64,   76,   77,   78,   72,   79,   80,   69,   81,   68,
       73,   82,   71,   75,   74,   83,   78,   84,   85,   79,
       80,   86,   76,   87,   81,   88,   89,   91,   83,   90,
       92,   93,   77,   94,  105,   84,   87,   98,   82,   99,
       88,   85,  101,   92,   91,  100,   94,  135,   89,  103,

      104,   86,   86,   86,   90,   98,  106,  100,  107,  112,
      101,   93,  105,  103,  104,  111,   99,  113,  114,   99,
       99,   99,  112,  115,  121,  123,  116,  128,  106,  126,
      124,  114,  131,  129,  123,  130,  127,  125,  107,  107,
      107,  113,  116,  115,  124,  126,  111,  122,  130,  120,
      121,  129,  128,  119,  118,  117,  110,  109,  131,  133,
      108,  133,  133,  133,  134,  134,  134,  136,  136,  102,
       97,   96,   95,   70,   57,   34,   29,   28,   27,   16,
       15,    9,    7,    5,  132,  132,  132,  132,  132,  132,
      132,  132,  132,  132,  132,  132,  132,  132,  132,  132,

      132,  132,  132,  132,  132,  132,  132,  132,  132,  132,
      132,  132,  132,  132,  132,  132,  132,  132,  132,  132,
      132,  132,  132,  132,  132
    } ;

static yy_state_type yy_last_accepting_state;
static char *yy_last_accepting_cpos;

/* The intent behind this definition is that it'll catch
 * any uses of REJECT which flex missed.
 */
#define REJECT reject_used_but_not_detected
#define yymore() yymore_used_but_not_detected
#define YY_MORE_ADJ 0
#define YY_RESTORE_YY_MORE_OFFSET
char *yytext;
#line 1 "convert.lex"
#define INITIAL 0
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
#define character 1

#line 19 "convert.lex"
#include <math.h>
#include <stdlib.h>
#include <string.h>
int line_num=1;
extern FILE * yyin;
#define MAX_INCLUDE_DEPTH 30
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
#line 505 "convert.yy.c"

/* Macros after this point can all be overridden by user definitions in
 * section 1.
 */

#ifndef YY_SKIP_YYWRAP
#ifdef __cplusplus
extern "C" int yywrap YY_PROTO(( void ));
#else
extern int yywrap YY_PROTO(( void ));
#endif
#endif

#ifndef YY_NO_UNPUT
static void yyunput YY_PROTO(( int c, char *buf_ptr ));
#endif

#ifndef yytext_ptr
static void yy_flex_strncpy YY_PROTO(( char *, yyconst char *, int ));
#endif

#ifdef YY_NEED_STRLEN
static int yy_flex_strlen YY_PROTO(( yyconst char * ));
#endif

#ifndef YY_NO_INPUT
#ifdef __cplusplus
static int yyinput YY_PROTO(( void ));
#else
static int input YY_PROTO(( void ));
#endif
#endif

#if YY_STACK_USED
static int yy_start_stack_ptr = 0;
static int yy_start_stack_depth = 0;
static int *yy_start_stack = 0;
#ifndef YY_NO_PUSH_STATE
static void yy_push_state YY_PROTO(( int new_state ));
#endif
#ifndef YY_NO_POP_STATE
static void yy_pop_state YY_PROTO(( void ));
#endif
#ifndef YY_NO_TOP_STATE
static int yy_top_state YY_PROTO(( void ));
#endif

#else
#define YY_NO_PUSH_STATE 1
#define YY_NO_POP_STATE 1
#define YY_NO_TOP_STATE 1
#endif

#ifdef YY_MALLOC_DECL
YY_MALLOC_DECL
#else
#if __STDC__
#ifndef __cplusplus
#include <stdlib.h>
#endif
#else
/* Just try to get by without declaring the routines.  This will fail
 * miserably on non-ANSI systems for which sizeof(size_t) != sizeof(int)
 * or sizeof(void*) != sizeof(int).
 */
#endif
#endif

/* Amount of stuff to slurp up with each read. */
#ifndef YY_READ_BUF_SIZE
#define YY_READ_BUF_SIZE 8192
#endif

/* Copy whatever the last rule matched to the standard output. */

#ifndef ECHO
/* This used to be an fputs(), but since the string might contain NUL's,
 * we now use fwrite().
 */
#define ECHO (void) fwrite( yytext, yyleng, 1, yyout )
#endif

/* Gets input and stuffs it into "buf".  number of characters read, or YY_NULL,
 * is returned in "result".
 */
#ifndef YY_INPUT
#define YY_INPUT(buf,result,max_size) \
	if ( yy_current_buffer->yy_is_interactive ) \
		{ \
		int c = '*', n; \
		for ( n = 0; n < max_size && \
			     (c = getc( yyin )) != EOF && c != '\n'; ++n ) \
			buf[n] = (char) c; \
		if ( c == '\n' ) \
			buf[n++] = (char) c; \
		if ( c == EOF && ferror( yyin ) ) \
			YY_FATAL_ERROR( "input in flex scanner failed" ); \
		result = n; \
		} \
	else if ( ((result = fread( buf, 1, max_size, yyin )) == 0) \
		  && ferror( yyin ) ) \
		YY_FATAL_ERROR( "input in flex scanner failed" );
#endif

/* No semi-colon after return; correct usage is to write "yyterminate();" -
 * we don't want an extra ';' after the "return" because that will cause
 * some compilers to complain about unreachable statements.
 */
#ifndef yyterminate
#define yyterminate() return YY_NULL
#endif

/* Number of entries by which start-condition stack grows. */
#ifndef YY_START_STACK_INCR
#define YY_START_STACK_INCR 25
#endif

/* Report a fatal error. */
#ifndef YY_FATAL_ERROR
#define YY_FATAL_ERROR(msg) yy_fatal_error( msg )
#endif

/* Default declaration of generated scanner - a define so the user can
 * easily add parameters.
 */
#ifndef YY_DECL
#define YY_DECL int yylex YY_PROTO(( void ))
#endif

/* Code executed at the beginning of each rule, after yytext and yyleng
 * have been set up.
 */
#ifndef YY_USER_ACTION
#define YY_USER_ACTION
#endif

/* Code executed at the end of each rule. */
#ifndef YY_BREAK
#define YY_BREAK break;
#endif

#define YY_RULE_SETUP \
	YY_USER_ACTION

YY_DECL
	{
	register yy_state_type yy_current_state;
	register char *yy_cp, *yy_bp;
	register int yy_act;

#line 39 "convert.lex"

#line 658 "convert.yy.c"

	if ( yy_init )
		{
		yy_init = 0;

#ifdef YY_USER_INIT
		YY_USER_INIT;
#endif

		if ( ! yy_start )
			yy_start = 1;	/* first start state */

		if ( ! yyin )
			yyin = stdin;

		if ( ! yyout )
			yyout = stdout;

		if ( ! yy_current_buffer )
			yy_current_buffer =
				yy_create_buffer( yyin, YY_BUF_SIZE );

		yy_load_buffer_state();
		}

	while ( 1 )		/* loops until end-of-file is reached */
		{
		yy_cp = yy_c_buf_p;

		/* Support of yytext. */
		*yy_cp = yy_hold_char;

		/* yy_bp points to the position in yy_ch_buf of the start of
		 * the current run.
		 */
		yy_bp = yy_cp;

		yy_current_state = yy_start;
yy_match:
		do
			{
			register YY_CHAR yy_c = yy_ec[YY_SC_TO_UI(*yy_cp)];
			if ( yy_accept[yy_current_state] )
				{
				yy_last_accepting_state = yy_current_state;
				yy_last_accepting_cpos = yy_cp;
				}
			while ( yy_chk[yy_base[yy_current_state] + yy_c] != yy_current_state )
				{
				yy_current_state = (int) yy_def[yy_current_state];
				if ( yy_current_state >= 133 )
					yy_c = yy_meta[(unsigned int) yy_c];
				}
			yy_current_state = yy_nxt[yy_base[yy_current_state] + (unsigned int) yy_c];
			++yy_cp;
			}
		while ( yy_base[yy_current_state] != 285 );

yy_find_action:
		yy_act = yy_accept[yy_current_state];
		if ( yy_act == 0 )
			{ /* have to back up */
			yy_cp = yy_last_accepting_cpos;
			yy_current_state = yy_last_accepting_state;
			yy_act = yy_accept[yy_current_state];
			}

		YY_DO_BEFORE_ACTION;


do_action:	/* This label is used only to access EOF actions. */


		switch ( yy_act )
	{ /* beginning of action switch */
			case 0: /* must back up */
			/* undo the effects of YY_DO_BEFORE_ACTION */
			*yy_cp = yy_hold_char;
			yy_cp = yy_last_accepting_cpos;
			yy_current_state = yy_last_accepting_state;
			goto yy_find_action;

case 1:
YY_RULE_SETUP
#line 40 "convert.lex"
return TOK_REGRIDDING; /* period of regridding                    */
	YY_BREAK
case 2:
YY_RULE_SETUP
#line 41 "convert.lex"
return TOK_COEFFRAFX;  /* space refinement in the x direction     */
	YY_BREAK
case 3:
YY_RULE_SETUP
#line 42 "convert.lex"
return TOK_COEFFRAFY;  /* space refinement in the y direction     */
	YY_BREAK
case 4:
YY_RULE_SETUP
#line 43 "convert.lex"
return TOK_COEFFRAFZ;  /* space refinement in the z direction     */
	YY_BREAK
case 5:
YY_RULE_SETUP
#line 44 "convert.lex"
return TOK_COEFFRAFTX; /* time refinement in the x direction      */
	YY_BREAK
case 6:
YY_RULE_SETUP
#line 45 "convert.lex"
return TOK_COEFFRAFTY; /* time refinement in the y direction      */
	YY_BREAK
case 7:
YY_RULE_SETUP
#line 46 "convert.lex"
return TOK_COEFFRAFTZ; /* time refinement in the z direction      */
	YY_BREAK
case 8:
YY_RULE_SETUP
#line 47 "convert.lex"
return TOK_MODULEMAIN; /* name of the module                      */
	YY_BREAK
case 9:
YY_RULE_SETUP
#line 48 "convert.lex"
return TOK_EFFICIENCY; /* efficiency for the adaptive refinement  */
	YY_BREAK
case 10:
YY_RULE_SETUP
#line 49 "convert.lex"
return TOK_RAFMAX;     /* minimum size in all directions          */
	YY_BREAK
case 11:
YY_RULE_SETUP
#line 50 "convert.lex"
return TOK_RAFMAXX;    /* minimum size in x direction             */
	YY_BREAK
case 12:
YY_RULE_SETUP
#line 51 "convert.lex"
return TOK_RAFMAXY;    /* minimum size in y direction             */
	YY_BREAK
case 13:
YY_RULE_SETUP
#line 52 "convert.lex"
return TOK_RAFMAXZ;    /* minimum size in z direction             */
	YY_BREAK
case 14:
YY_RULE_SETUP
#line 53 "convert.lex"
return TOK_NOTGRIDDEP; /* variable which are not grid dependent   */
	YY_BREAK
case 15:
YY_RULE_SETUP
#line 54 "convert.lex"
return TOK_USE;
	YY_BREAK
case 16:
YY_RULE_SETUP
#line 55 "convert.lex"
return TOK_MINWIDTH;   /* minimum width of rectangles for the     */
	YY_BREAK
/*    adaptive refinement                  */
case 17:
YY_RULE_SETUP
#line 57 "convert.lex"
{}
	YY_BREAK
case 18:
YY_RULE_SETUP
#line 58 "convert.lex"
return TOK_SEP;
	YY_BREAK
case 19:
YY_RULE_SETUP
#line 59 "convert.lex"
{strcpy(yylval.na,yytext); return TOK_FILENAME;}
	YY_BREAK
case 20:
YY_RULE_SETUP
#line 60 "convert.lex"
{strcpy(yylval.na,yytext); return TOK_USEITEM;}
	YY_BREAK
case 21:
YY_RULE_SETUP
#line 61 "convert.lex"
{strcpy(yylval.na,yytext); return TOK_PROBTYPE;}
	YY_BREAK
/* dimension of the problem                */
case 22:
YY_RULE_SETUP
#line 63 "convert.lex"
{strcpy(yylval.na,yytext); return TOK_NAME;}
	YY_BREAK
case 23:
YY_RULE_SETUP
#line 64 "convert.lex"
{yylval.ival=atoi(yytext); return TOK_NUM;}
	YY_BREAK
case 24:
YY_RULE_SETUP
#line 65 "convert.lex"
{return (int) *yytext;}
	YY_BREAK
case 25:
YY_RULE_SETUP
#line 66 "convert.lex"
{line_num++;return (int) *yytext;}
	YY_BREAK
case 26:
YY_RULE_SETUP
#line 67 "convert.lex"
;
	YY_BREAK
case 27:
YY_RULE_SETUP
#line 68 "convert.lex"
ECHO;
	YY_BREAK
#line 878 "convert.yy.c"
case YY_STATE_EOF(INITIAL):
case YY_STATE_EOF(character):
	yyterminate();

	case YY_END_OF_BUFFER:
		{
		/* Amount of text matched not including the EOB char. */
		int yy_amount_of_matched_text = (int) (yy_cp - yytext_ptr) - 1;

		/* Undo the effects of YY_DO_BEFORE_ACTION. */
		*yy_cp = yy_hold_char;
		YY_RESTORE_YY_MORE_OFFSET

		if ( yy_current_buffer->yy_buffer_status == YY_BUFFER_NEW )
			{
			/* We're scanning a new file or input source.  It's
			 * possible that this happened because the user
			 * just pointed yyin at a new source and called
			 * yylex().  If so, then we have to assure
			 * consistency between yy_current_buffer and our
			 * globals.  Here is the right place to do so, because
			 * this is the first action (other than possibly a
			 * back-up) that will match for the new input source.
			 */
			yy_n_chars = yy_current_buffer->yy_n_chars;
			yy_current_buffer->yy_input_file = yyin;
			yy_current_buffer->yy_buffer_status = YY_BUFFER_NORMAL;
			}

		/* Note that here we test for yy_c_buf_p "<=" to the position
		 * of the first EOB in the buffer, since yy_c_buf_p will
		 * already have been incremented past the NUL character
		 * (since all states make transitions on EOB to the
		 * end-of-buffer state).  Contrast this with the test
		 * in input().
		 */
		if ( yy_c_buf_p <= &yy_current_buffer->yy_ch_buf[yy_n_chars] )
			{ /* This was really a NUL. */
			yy_state_type yy_next_state;

			yy_c_buf_p = yytext_ptr + yy_amount_of_matched_text;

			yy_current_state = yy_get_previous_state();

			/* Okay, we're now positioned to make the NUL
			 * transition.  We couldn't have
			 * yy_get_previous_state() go ahead and do it
			 * for us because it doesn't know how to deal
			 * with the possibility of jamming (and we don't
			 * want to build jamming into it because then it
			 * will run more slowly).
			 */

			yy_next_state = yy_try_NUL_trans( yy_current_state );

			yy_bp = yytext_ptr + YY_MORE_ADJ;

			if ( yy_next_state )
				{
				/* Consume the NUL. */
				yy_cp = ++yy_c_buf_p;
				yy_current_state = yy_next_state;
				goto yy_match;
				}

			else
				{
				yy_cp = yy_c_buf_p;
				goto yy_find_action;
				}
			}

		else switch ( yy_get_next_buffer() )
			{
			case EOB_ACT_END_OF_FILE:
				{
				yy_did_buffer_switch_on_eof = 0;

				if ( yywrap() )
					{
					/* Note: because we've taken care in
					 * yy_get_next_buffer() to have set up
					 * yytext, we can now set up
					 * yy_c_buf_p so that if some total
					 * hoser (like flex itself) wants to
					 * call the scanner after we return the
					 * YY_NULL, it'll still work - another
					 * YY_NULL will get returned.
					 */
					yy_c_buf_p = yytext_ptr + YY_MORE_ADJ;

					yy_act = YY_STATE_EOF(YY_START);
					goto do_action;
					}

				else
					{
					if ( ! yy_did_buffer_switch_on_eof )
						YY_NEW_FILE;
					}
				break;
				}

			case EOB_ACT_CONTINUE_SCAN:
				yy_c_buf_p =
					yytext_ptr + yy_amount_of_matched_text;

				yy_current_state = yy_get_previous_state();

				yy_cp = yy_c_buf_p;
				yy_bp = yytext_ptr + YY_MORE_ADJ;
				goto yy_match;

			case EOB_ACT_LAST_MATCH:
				yy_c_buf_p =
				&yy_current_buffer->yy_ch_buf[yy_n_chars];

				yy_current_state = yy_get_previous_state();

				yy_cp = yy_c_buf_p;
				yy_bp = yytext_ptr + YY_MORE_ADJ;
				goto yy_find_action;
			}
		break;
		}

	default:
		YY_FATAL_ERROR(
			"fatal flex scanner internal error--no action found" );
	} /* end of action switch */
		} /* end of scanning one token */
	} /* end of yylex */


/* yy_get_next_buffer - try to read in a new buffer
 *
 * Returns a code representing an action:
 *	EOB_ACT_LAST_MATCH -
 *	EOB_ACT_CONTINUE_SCAN - continue scanning from current position
 *	EOB_ACT_END_OF_FILE - end of file
 */

static int yy_get_next_buffer()
	{
	register char *dest = yy_current_buffer->yy_ch_buf;
	register char *source = yytext_ptr;
	register int number_to_move, i;
	int ret_val;

	if ( yy_c_buf_p > &yy_current_buffer->yy_ch_buf[yy_n_chars + 1] )
		YY_FATAL_ERROR(
		"fatal flex scanner internal error--end of buffer missed" );

	if ( yy_current_buffer->yy_fill_buffer == 0 )
		{ /* Don't try to fill the buffer, so this is an EOF. */
		if ( yy_c_buf_p - yytext_ptr - YY_MORE_ADJ == 1 )
			{
			/* We matched a single character, the EOB, so
			 * treat this as a final EOF.
			 */
			return EOB_ACT_END_OF_FILE;
			}

		else
			{
			/* We matched some text prior to the EOB, first
			 * process it.
			 */
			return EOB_ACT_LAST_MATCH;
			}
		}

	/* Try to read more data. */

	/* First move last chars to start of buffer. */
	number_to_move = (int) (yy_c_buf_p - yytext_ptr) - 1;

	for ( i = 0; i < number_to_move; ++i )
		*(dest++) = *(source++);

	if ( yy_current_buffer->yy_buffer_status == YY_BUFFER_EOF_PENDING )
		/* don't do the read, it's not guaranteed to return an EOF,
		 * just force an EOF
		 */
		yy_current_buffer->yy_n_chars = yy_n_chars = 0;

	else
		{
		int num_to_read =
			yy_current_buffer->yy_buf_size - number_to_move - 1;

		while ( num_to_read <= 0 )
			{ /* Not enough room in the buffer - grow it. */
#ifdef YY_USES_REJECT
			YY_FATAL_ERROR(
"input buffer overflow, can't enlarge buffer because scanner uses REJECT" );
#else

			/* just a shorter name for the current buffer */
			YY_BUFFER_STATE b = yy_current_buffer;

			int yy_c_buf_p_offset =
				(int) (yy_c_buf_p - b->yy_ch_buf);

			if ( b->yy_is_our_buffer )
				{
				int new_size = b->yy_buf_size * 2;

				if ( new_size <= 0 )
					b->yy_buf_size += b->yy_buf_size / 8;
				else
					b->yy_buf_size *= 2;

				b->yy_ch_buf = (char *)
					/* Include room in for 2 EOB chars. */
					yy_flex_realloc( (void *) b->yy_ch_buf,
							 b->yy_buf_size + 2 );
				}
			else
				/* Can't grow it, we don't own it. */
				b->yy_ch_buf = 0;

			if ( ! b->yy_ch_buf )
				YY_FATAL_ERROR(
				"fatal error - scanner input buffer overflow" );

			yy_c_buf_p = &b->yy_ch_buf[yy_c_buf_p_offset];

			num_to_read = yy_current_buffer->yy_buf_size -
						number_to_move - 1;
#endif
			}

		if ( num_to_read > YY_READ_BUF_SIZE )
			num_to_read = YY_READ_BUF_SIZE;

		/* Read in more data. */
		YY_INPUT( (&yy_current_buffer->yy_ch_buf[number_to_move]),
			yy_n_chars, num_to_read );

		yy_current_buffer->yy_n_chars = yy_n_chars;
		}

	if ( yy_n_chars == 0 )
		{
		if ( number_to_move == YY_MORE_ADJ )
			{
			ret_val = EOB_ACT_END_OF_FILE;
			yyrestart( yyin );
			}

		else
			{
			ret_val = EOB_ACT_LAST_MATCH;
			yy_current_buffer->yy_buffer_status =
				YY_BUFFER_EOF_PENDING;
			}
		}

	else
		ret_val = EOB_ACT_CONTINUE_SCAN;

	yy_n_chars += number_to_move;
	yy_current_buffer->yy_ch_buf[yy_n_chars] = YY_END_OF_BUFFER_CHAR;
	yy_current_buffer->yy_ch_buf[yy_n_chars + 1] = YY_END_OF_BUFFER_CHAR;

	yytext_ptr = &yy_current_buffer->yy_ch_buf[0];

	return ret_val;
	}


/* yy_get_previous_state - get the state just before the EOB char was reached */

static yy_state_type yy_get_previous_state()
	{
	register yy_state_type yy_current_state;
	register char *yy_cp;

	yy_current_state = yy_start;

	for ( yy_cp = yytext_ptr + YY_MORE_ADJ; yy_cp < yy_c_buf_p; ++yy_cp )
		{
		register YY_CHAR yy_c = (*yy_cp ? yy_ec[YY_SC_TO_UI(*yy_cp)] : 1);
		if ( yy_accept[yy_current_state] )
			{
			yy_last_accepting_state = yy_current_state;
			yy_last_accepting_cpos = yy_cp;
			}
		while ( yy_chk[yy_base[yy_current_state] + yy_c] != yy_current_state )
			{
			yy_current_state = (int) yy_def[yy_current_state];
			if ( yy_current_state >= 133 )
				yy_c = yy_meta[(unsigned int) yy_c];
			}
		yy_current_state = yy_nxt[yy_base[yy_current_state] + (unsigned int) yy_c];
		}

	return yy_current_state;
	}


/* yy_try_NUL_trans - try to make a transition on the NUL character
 *
 * synopsis
 *	next_state = yy_try_NUL_trans( current_state );
 */

#ifdef YY_USE_PROTOS
static yy_state_type yy_try_NUL_trans( yy_state_type yy_current_state )
#else
static yy_state_type yy_try_NUL_trans( yy_current_state )
yy_state_type yy_current_state;
#endif
	{
	register int yy_is_jam;
	register char *yy_cp = yy_c_buf_p;

	register YY_CHAR yy_c = 1;
	if ( yy_accept[yy_current_state] )
		{
		yy_last_accepting_state = yy_current_state;
		yy_last_accepting_cpos = yy_cp;
		}
	while ( yy_chk[yy_base[yy_current_state] + yy_c] != yy_current_state )
		{
		yy_current_state = (int) yy_def[yy_current_state];
		if ( yy_current_state >= 133 )
			yy_c = yy_meta[(unsigned int) yy_c];
		}
	yy_current_state = yy_nxt[yy_base[yy_current_state] + (unsigned int) yy_c];
	yy_is_jam = (yy_current_state == 132);

	return yy_is_jam ? 0 : yy_current_state;
	}


#ifndef YY_NO_UNPUT
#ifdef YY_USE_PROTOS
static void yyunput( int c, register char *yy_bp )
#else
static void yyunput( c, yy_bp )
int c;
register char *yy_bp;
#endif
	{
	register char *yy_cp = yy_c_buf_p;

	/* undo effects of setting up yytext */
	*yy_cp = yy_hold_char;

	if ( yy_cp < yy_current_buffer->yy_ch_buf + 2 )
		{ /* need to shift things up to make room */
		/* +2 for EOB chars. */
		register int number_to_move = yy_n_chars + 2;
		register char *dest = &yy_current_buffer->yy_ch_buf[
					yy_current_buffer->yy_buf_size + 2];
		register char *source =
				&yy_current_buffer->yy_ch_buf[number_to_move];

		while ( source > yy_current_buffer->yy_ch_buf )
			*--dest = *--source;

		yy_cp += (int) (dest - source);
		yy_bp += (int) (dest - source);
		yy_current_buffer->yy_n_chars =
			yy_n_chars = yy_current_buffer->yy_buf_size;

		if ( yy_cp < yy_current_buffer->yy_ch_buf + 2 )
			YY_FATAL_ERROR( "flex scanner push-back overflow" );
		}

	*--yy_cp = (char) c;


	yytext_ptr = yy_bp;
	yy_hold_char = *yy_cp;
	yy_c_buf_p = yy_cp;
	}
#endif	/* ifndef YY_NO_UNPUT */


#ifdef __cplusplus
static int yyinput()
#else
static int input()
#endif
	{
	int c;

	*yy_c_buf_p = yy_hold_char;

	if ( *yy_c_buf_p == YY_END_OF_BUFFER_CHAR )
		{
		/* yy_c_buf_p now points to the character we want to return.
		 * If this occurs *before* the EOB characters, then it's a
		 * valid NUL; if not, then we've hit the end of the buffer.
		 */
		if ( yy_c_buf_p < &yy_current_buffer->yy_ch_buf[yy_n_chars] )
			/* This was really a NUL. */
			*yy_c_buf_p = '\0';

		else
			{ /* need more input */
			int offset = yy_c_buf_p - yytext_ptr;
			++yy_c_buf_p;

			switch ( yy_get_next_buffer() )
				{
				case EOB_ACT_LAST_MATCH:
					/* This happens because yy_g_n_b()
					 * sees that we've accumulated a
					 * token and flags that we need to
					 * try matching the token before
					 * proceeding.  But for input(),
					 * there's no matching to consider.
					 * So convert the EOB_ACT_LAST_MATCH
					 * to EOB_ACT_END_OF_FILE.
					 */

					/* Reset buffer status. */
					yyrestart( yyin );

					/* fall through */

				case EOB_ACT_END_OF_FILE:
					{
					if ( yywrap() )
						return EOF;

					if ( ! yy_did_buffer_switch_on_eof )
						YY_NEW_FILE;
#ifdef __cplusplus
					return yyinput();
#else
					return input();
#endif
					}

				case EOB_ACT_CONTINUE_SCAN:
					yy_c_buf_p = yytext_ptr + offset;
					break;
				}
			}
		}

	c = *(unsigned char *) yy_c_buf_p;	/* cast for 8-bit char's */
	*yy_c_buf_p = '\0';	/* preserve yytext */
	yy_hold_char = *++yy_c_buf_p;


	return c;
	}


#ifdef YY_USE_PROTOS
void yyrestart( FILE *input_file )
#else
void yyrestart( input_file )
FILE *input_file;
#endif
	{
	if ( ! yy_current_buffer )
		yy_current_buffer = yy_create_buffer( yyin, YY_BUF_SIZE );

	yy_init_buffer( yy_current_buffer, input_file );
	yy_load_buffer_state();
	}


#ifdef YY_USE_PROTOS
void yy_switch_to_buffer( YY_BUFFER_STATE new_buffer )
#else
void yy_switch_to_buffer( new_buffer )
YY_BUFFER_STATE new_buffer;
#endif
	{
	if ( yy_current_buffer == new_buffer )
		return;

	if ( yy_current_buffer )
		{
		/* Flush out information for old buffer. */
		*yy_c_buf_p = yy_hold_char;
		yy_current_buffer->yy_buf_pos = yy_c_buf_p;
		yy_current_buffer->yy_n_chars = yy_n_chars;
		}

	yy_current_buffer = new_buffer;
	yy_load_buffer_state();

	/* We don't actually know whether we did this switch during
	 * EOF (yywrap()) processing, but the only time this flag
	 * is looked at is after yywrap() is called, so it's safe
	 * to go ahead and always set it.
	 */
	yy_did_buffer_switch_on_eof = 1;
	}


#ifdef YY_USE_PROTOS
void yy_load_buffer_state( void )
#else
void yy_load_buffer_state()
#endif
	{
	yy_n_chars = yy_current_buffer->yy_n_chars;
	yytext_ptr = yy_c_buf_p = yy_current_buffer->yy_buf_pos;
	yyin = yy_current_buffer->yy_input_file;
	yy_hold_char = *yy_c_buf_p;
	}


#ifdef YY_USE_PROTOS
YY_BUFFER_STATE yy_create_buffer( FILE *file, int size )
#else
YY_BUFFER_STATE yy_create_buffer( file, size )
FILE *file;
int size;
#endif
	{
	YY_BUFFER_STATE b;

	b = (YY_BUFFER_STATE) yy_flex_alloc( sizeof( struct yy_buffer_state ) );
	if ( ! b )
		YY_FATAL_ERROR( "out of dynamic memory in yy_create_buffer()" );

	b->yy_buf_size = size;

	/* yy_ch_buf has to be 2 characters longer than the size given because
	 * we need to put in 2 end-of-buffer characters.
	 */
	b->yy_ch_buf = (char *) yy_flex_alloc( b->yy_buf_size + 2 );
	if ( ! b->yy_ch_buf )
		YY_FATAL_ERROR( "out of dynamic memory in yy_create_buffer()" );

	b->yy_is_our_buffer = 1;

	yy_init_buffer( b, file );

	return b;
	}


#ifdef YY_USE_PROTOS
void yy_delete_buffer( YY_BUFFER_STATE b )
#else
void yy_delete_buffer( b )
YY_BUFFER_STATE b;
#endif
	{
	if ( ! b )
		return;

	if ( b == yy_current_buffer )
		yy_current_buffer = (YY_BUFFER_STATE) 0;

	if ( b->yy_is_our_buffer )
		yy_flex_free( (void *) b->yy_ch_buf );

	yy_flex_free( (void *) b );
	}



#ifdef YY_USE_PROTOS
void yy_init_buffer( YY_BUFFER_STATE b, FILE *file )
#else
void yy_init_buffer( b, file )
YY_BUFFER_STATE b;
FILE *file;
#endif


	{
	yy_flush_buffer( b );

	b->yy_input_file = file;
	b->yy_fill_buffer = 1;

#if YY_ALWAYS_INTERACTIVE
	b->yy_is_interactive = 1;
#else
#if YY_NEVER_INTERACTIVE
	b->yy_is_interactive = 0;
#else
	b->yy_is_interactive = file ? (isatty( fileno(file) ) > 0) : 0;
#endif
#endif
	}


#ifdef YY_USE_PROTOS
void yy_flush_buffer( YY_BUFFER_STATE b )
#else
void yy_flush_buffer( b )
YY_BUFFER_STATE b;
#endif

	{
	if ( ! b )
		return;

	b->yy_n_chars = 0;

	/* We always need two end-of-buffer characters.  The first causes
	 * a transition to the end-of-buffer state.  The second causes
	 * a jam in that state.
	 */
	b->yy_ch_buf[0] = YY_END_OF_BUFFER_CHAR;
	b->yy_ch_buf[1] = YY_END_OF_BUFFER_CHAR;

	b->yy_buf_pos = &b->yy_ch_buf[0];

	b->yy_at_bol = 1;
	b->yy_buffer_status = YY_BUFFER_NEW;

	if ( b == yy_current_buffer )
		yy_load_buffer_state();
	}


#ifndef YY_NO_SCAN_BUFFER
#ifdef YY_USE_PROTOS
YY_BUFFER_STATE yy_scan_buffer( char *base, yy_size_t size )
#else
YY_BUFFER_STATE yy_scan_buffer( base, size )
char *base;
yy_size_t size;
#endif
	{
	YY_BUFFER_STATE b;

	if ( size < 2 ||
	     base[size-2] != YY_END_OF_BUFFER_CHAR ||
	     base[size-1] != YY_END_OF_BUFFER_CHAR )
		/* They forgot to leave room for the EOB's. */
		return 0;

	b = (YY_BUFFER_STATE) yy_flex_alloc( sizeof( struct yy_buffer_state ) );
	if ( ! b )
		YY_FATAL_ERROR( "out of dynamic memory in yy_scan_buffer()" );

	b->yy_buf_size = size - 2;	/* "- 2" to take care of EOB's */
	b->yy_buf_pos = b->yy_ch_buf = base;
	b->yy_is_our_buffer = 0;
	b->yy_input_file = 0;
	b->yy_n_chars = b->yy_buf_size;
	b->yy_is_interactive = 0;
	b->yy_at_bol = 1;
	b->yy_fill_buffer = 0;
	b->yy_buffer_status = YY_BUFFER_NEW;

	yy_switch_to_buffer( b );

	return b;
	}
#endif


#ifndef YY_NO_SCAN_STRING
#ifdef YY_USE_PROTOS
YY_BUFFER_STATE yy_scan_string( yyconst char *yy_str )
#else
YY_BUFFER_STATE yy_scan_string( yy_str )
yyconst char *yy_str;
#endif
	{
	int len;
	for ( len = 0; yy_str[len]; ++len )
		;

	return yy_scan_bytes( yy_str, len );
	}
#endif


#ifndef YY_NO_SCAN_BYTES
#ifdef YY_USE_PROTOS
YY_BUFFER_STATE yy_scan_bytes( yyconst char *bytes, int len )
#else
YY_BUFFER_STATE yy_scan_bytes( bytes, len )
yyconst char *bytes;
int len;
#endif
	{
	YY_BUFFER_STATE b;
	char *buf;
	yy_size_t n;
	int i;

	/* Get memory for full buffer, including space for trailing EOB's. */
	n = len + 2;
	buf = (char *) yy_flex_alloc( n );
	if ( ! buf )
		YY_FATAL_ERROR( "out of dynamic memory in yy_scan_bytes()" );

	for ( i = 0; i < len; ++i )
		buf[i] = bytes[i];

	buf[len] = buf[len+1] = YY_END_OF_BUFFER_CHAR;

	b = yy_scan_buffer( buf, n );
	if ( ! b )
		YY_FATAL_ERROR( "bad buffer in yy_scan_bytes()" );

	/* It's okay to grow etc. this buffer, and we should throw it
	 * away when we're done.
	 */
	b->yy_is_our_buffer = 1;

	return b;
	}
#endif


#ifndef YY_NO_PUSH_STATE
#ifdef YY_USE_PROTOS
static void yy_push_state( int new_state )
#else
static void yy_push_state( new_state )
int new_state;
#endif
	{
	if ( yy_start_stack_ptr >= yy_start_stack_depth )
		{
		yy_size_t new_size;

		yy_start_stack_depth += YY_START_STACK_INCR;
		new_size = yy_start_stack_depth * sizeof( int );

		if ( ! yy_start_stack )
			yy_start_stack = (int *) yy_flex_alloc( new_size );

		else
			yy_start_stack = (int *) yy_flex_realloc(
					(void *) yy_start_stack, new_size );

		if ( ! yy_start_stack )
			YY_FATAL_ERROR(
			"out of memory expanding start-condition stack" );
		}

	yy_start_stack[yy_start_stack_ptr++] = YY_START;

	BEGIN(new_state);
	}
#endif


#ifndef YY_NO_POP_STATE
static void yy_pop_state()
	{
	if ( --yy_start_stack_ptr < 0 )
		YY_FATAL_ERROR( "start-condition stack underflow" );

	BEGIN(yy_start_stack[yy_start_stack_ptr]);
	}
#endif


#ifndef YY_NO_TOP_STATE
static int yy_top_state()
	{
	return yy_start_stack[yy_start_stack_ptr - 1];
	}
#endif

#ifndef YY_EXIT_FAILURE
#define YY_EXIT_FAILURE 2
#endif

#ifdef YY_USE_PROTOS
static void yy_fatal_error( yyconst char msg[] )
#else
static void yy_fatal_error( msg )
char msg[];
#endif
	{
	(void) fprintf( stderr, "%s\n", msg );
	exit( YY_EXIT_FAILURE );
	}



/* Redefine yyless() so it works in section 3 code. */

#undef yyless
#define yyless(n) \
	do \
		{ \
		/* Undo effects of setting up yytext. */ \
		yytext[yyleng] = yy_hold_char; \
		yy_c_buf_p = yytext + n; \
		yy_hold_char = *yy_c_buf_p; \
		*yy_c_buf_p = '\0'; \
		yyleng = n; \
		} \
	while ( 0 )


/* Internal utility routines. */

#ifndef yytext_ptr
#ifdef YY_USE_PROTOS
static void yy_flex_strncpy( char *s1, yyconst char *s2, int n )
#else
static void yy_flex_strncpy( s1, s2, n )
char *s1;
yyconst char *s2;
int n;
#endif
	{
	register int i;
	for ( i = 0; i < n; ++i )
		s1[i] = s2[i];
	}
#endif

#ifdef YY_NEED_STRLEN
#ifdef YY_USE_PROTOS
static int yy_flex_strlen( yyconst char *s )
#else
static int yy_flex_strlen( s )
yyconst char *s;
#endif
	{
	register int n;
	for ( n = 0; s[n]; ++n )
		;

	return n;
	}
#endif


#ifdef YY_USE_PROTOS
static void *yy_flex_alloc( yy_size_t size )
#else
static void *yy_flex_alloc( size )
yy_size_t size;
#endif
	{
	return (void *) malloc( size );
	}

#ifdef YY_USE_PROTOS
static void *yy_flex_realloc( void *ptr, yy_size_t size )
#else
static void *yy_flex_realloc( ptr, size )
void *ptr;
yy_size_t size;
#endif
	{
	/* The cast to (char *) in the following accommodates both
	 * implementations that use char* generic pointers, and those
	 * that use void* generic pointers.  It works with the latter
	 * because both ANSI C and C++ allow castless assignment from
	 * any pointer type to void*, and deal with argument conversions
	 * as though doing an assignment.
	 */
	return (void *) realloc( (char *) ptr, size );
	}

#ifdef YY_USE_PROTOS
static void yy_flex_free( void *ptr )
#else
static void yy_flex_free( ptr )
void *ptr;
#endif
	{
	free( ptr );
	}

#if YY_MAIN
int main()
	{
	yylex();
	return 0;
	}
#endif
#line 68 "convert.lex"



int yywrap()
{
}


yyerror(char *s)
{
if (!strcasecmp(curfile,mainfile))
{
   printf("Dans convert %s line %d, fichier %s\n",s,line_num,curfile);
}
else
{
   printf("Dans convert %s line %d, fichier %s\n",s,line_num,curfile);
}
exit(0);
}
