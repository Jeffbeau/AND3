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
#define LONGNOM 800
#define LONGLIGNE 800

/******************************************************************************/
/*********** Declaration of structures used in conv ***************************/
/******************************************************************************/

typedef struct 
{
   char first[LONGNOM];
   char last[LONGNOM];
} typedim ;                /* fortran dimension as 'ndeb:nfin'                */

typedef struct listdim 
{   
   typedim dim;
   struct listdim *suiv;
} listdim;                 /* list of the dimensions of a variable            */
      
typedef struct variable 
{
   char typevar[LONGNOM];
   char nomvar[LONGNOM] ;
   char oldname[LONGNOM] ;
   char dimchar[LONGNOM];
   listdim *dimension;
   int nbdim;
   struct variable *vinit;
   int common;
   int positioninblock;
   int module; 
   int save;
   int VariableIsParameter;
   int PublicDeclare;
   int PrivateDeclare;
   int ExternalDeclare;
   char modulename[LONGNOM]; 
   char commonname[LONGNOM];
   char vallengspec[LONGNOM];
   char nameinttypename[LONGNOM];
   int lengspecgiven;
   int pointedvar;
   char commoninfile[LONGNOM];
   char subroutinename[LONGNOM];
   int dimensiongiven;
   int c_star;
   int typegiven;
   int isparameter;
   char precision[LONGNOM]; 
   char initialvalue[LONGNOM]; 
   int indicetabvars; 
   int pointerdeclare; 
   int optionaldeclare;
   int allocatable; 
   char IntentSpec[LONGNOM];    
   int dimsempty;
   char readedlistdimension[LONGNOM];    
} variable ;               /* type of a variable                              */
                           /* typevar : type (integer, real, ...)             */
                           /* nomvar : name of the variable                   */
                           /* dimension : list of dimensions of the variable  */ 
                           /* nbdim: 1 if the variable is 1d, etc ...         */
                           /* precision : Name of the variable which          */
                           /* determine the precision. example : wp in the    */
                           /* case where REAL(wp)                             */
         
typedef struct listvar
{
   variable *var ;
   struct listvar * suiv;
} listvar ;                /* list of variables                               */
  

typedef struct listvarcommon
{
   char nomvar[LONGNOM] ;
   char commonname[LONGNOM];
   char commoninfile[LONGNOM];
   char subroutinename[LONGNOM];
   int dimensiongiven;
   int nbdim;
   int indicetabvars;
   int positioninblock;
   listdim *dimension;
   char readedlistdimension[LONGNOM];    
   struct listvarcommon * suiv;
} listvarcommon ;          /* list of variables in common block               */
  

typedef struct listusemodule 
{
   char usemodule[LONGNOM];
   char charusemodule[LONGNOM];
   char cursubroutine[LONGNOM];
   int firstuse;
   struct listusemodule * suiv;
} listusemodule;           /* list of names                                   */

typedef struct listparameter
{
   char name[LONGNOM];
   char modulename[LONGNOM];
   struct listparameter * suiv;
} listparameter ;           /* list of names                                  */

typedef struct listnamelist
{
   char name[LONGNOM];
   struct listnamelist * suiv;
} listnamelist ;            /* list of names                                  */

typedef struct listname
{
   char name[LONGNOM];
   struct  listname* suiv;
} listname ;            /* list of names                                  */

listname *listimplicitnone;

typedef struct listmodule 
{
   char module[LONGNOM];
   int InstanceShouldMade;
   int Instance;
   struct listmodule * suiv;
} listmodule;              /* list of names                                   */



typedef struct listcouple 
{
   char namevar[LONGNOM];
   char namepointedvar[LONGNOM];
   struct listcouple * suiv;
} listcouple;              /* list of names                                   */


typedef struct listnom 
{
   char nom[LONGNOM];
   listcouple *couple;
   struct listnom * suiv;
} listnom;                 /* list of names                                   */


typedef struct listallocate 
{
   char nomvar[LONGNOM];
   char subroutine[LONGNOM];
   char module[LONGNOM];
   struct listallocate * suiv;
} listallocate ;

 
typedef struct listvarpointtovar 
{
   char usemodule[LONGNOM];
   char cursubroutine[LONGNOM];
   listcouple *couple;
   struct  listvarpointtovar* suiv;
}listvarpointtovar ;       /* list of names                                   */


typedef struct listindice 
{
   int indice;
   struct  listindice * suiv;
} listindice;              /* list of indiced                                 */
 

typedef struct listparameters 
{
   char initial[LONGNOM];
   char nom[LONGNOM];
   char subroutinename[LONGNOM];
   struct listparameters * suiv;
 } listparameters;         /* list of parameters                              */

typedef struct 
{
   char name[LONGNOM];
   char filename[LONGNOM];
   int inclcommon;
   listvar *lvar;
} common;

typedef struct listcommons
{
 common *curcommon ;
 struct listcommons *suiv;
} listcommons ;


typedef struct listcommonsname
{
 char nom[LONGNOM] ;
 listnom *subroutine;
 struct listcommonsname *suiv;
} listcommonsname ;



 int fortran77;            /* = 1; the code has been writen in                */
                           /*    fortran77 else in fortran 90                 */
/******************************************************************************/
/****************   *** COMMON Variables ***  *********************************/
/******************************************************************************/

 int positioninblock;
 char commonvar[LONGNOM];
 char commonblockname[LONGNOM];
 listdim *commondim;

/******************************************************************************/
/****************   *** AGRIF Variables ***   *********************************/
/******************************************************************************/
 int inagrifcallargument;
 int adduseagrifutil;

 int InAgrifParentDef;

/******************************************************************************/
/****************   *** VAR DEF Variables ***   *******************************/
/******************************************************************************/
 int oldindicemaxtabvars;  /* Number of variables in the model i.e. last      */
 int indicemaxtabvars;     /* Number of variables in the model i.e. last      */
                           /*    indice used in  the tabvars table            */
 int PublicDeclare;        /* Variable has been declared as PUBLIC */ 
 int PrivateDeclare;       /* Variable has been declared as PRIVATE */ 
 int ExternalDeclare;      /* Variable has been declared as EXTERNAL */ 
 int PrecisionGiven;       /* A precision has been given for the variable */ 
 int CharacterSizeGiven;   /* A size for the character has been given */ 
 int InitialValueGiven;    /* An initial value has been given */ 
 int formatdeclare;
 int inttypename;
 int Allocatabledeclare;
 int lengspecgiven;
 int SaveDeclare;
 int pointerdeclare;
 int optionaldeclare;
 int VariableIsParameter; 
 int dimsgiven;
 int IntentDeclare;
 int c_star;
 char DeclType[LONGNOM]; 
 char nameinttypename[LONGNOM]; 
 char InitValue[LONGNOM*2]; 
 char IntentSpec[LONGNOM];
 char NamePrecision[LONGNOM]; 
 char CharacterSize[LONGNOM]; 
 char curmodulename[LONGNOM];
 char vallengspec[LONGNOM];
 char subroutinename[LONGNOM];

/******************************************************************************/
/****************   *** TOAMR Variables ***   *********************************/
/******************************************************************************/
 char Alloctreatedname[LONGNOM];

/******************************************************************************/
/****************   *** CONV Variables ***   **********************************/
/******************************************************************************/
 int coeffrafx;            /* space refinement factor in the x                */
 int coeffrafy;            /* space refinement factor in the y                */ 
 int coeffrafz;            /* space refinement factor in the z                */
 int coeffraftx;           /* time refinement factor in the x                 */
 int coeffrafty;           /* time refinement factor in the y                 */
 int coeffraftz;           /* time refinement factor in the z                 */
 int regridding;           /* number of time steps between two regridding     */
 int dimprob ;             /* dimension of the problem : 1 for 1D,2 for 2D,   */
                           /*    3 for 3D                                     */
 int rafmaxx;
 int rafmaxy;
 int rafmaxz;
 int userefficiency;       /* = 1 efficiency is given in the input            */
 int minwidth;             /* minimum width of the rectangles in the          */
                           /*    clustering algorithm                         */
 int onlyfixedgrids;       /* = 1 if onlyfixedgrids is true                   */
 int todebug;
 int fixedgrids;           /* = 1 if fixedgrids is true                       */
 int efficiency;           /* efficacity of the clustering algorithm          */
                           /*    (by percent)                                 */
 char nbmaillesX[LONGNOM]; /* number of cells in the x direction              */
 char nbmaillesY[LONGNOM]; /* number of cells in the y direction              */
 char nbmaillesZ[LONGNOM]; /* number of cells in the z direction              */
 int IndicenbmaillesX;
 int IndicenbmaillesY;
 int IndicenbmaillesZ;

 listvar *listvartempo;
 variable *variabletempo;

 listdim *curdim;
 variable *curvar;
 listvar *globliste;
 listvar *globvarforcommon;
 listvar *globparam;
 listvar *listdatavariable;
 listvar *listargsubroutine;
 listvar *varofsubroutineliste;
 listvar *varsubroutine;
 listvar *listvarindoloop;
 listvar *finglobliste;
 listvar *tmplocallist;
 listvar *parameterlist;
 listvar *globalvarofusefile2;
 listvar *globalvarofusefile;
 listvar *functionlistvar;
 listvar *listenotgriddepend; /* List of the variables which are not grid dependent */
 listvarcommon *commonlist;
 listusemodule *listofmodulebysubroutine;
 listusemodule *listofincludebysubroutine;
 listusemodule *listofmoduletmp;
 listusemodule *tmpuselocallist;
 listparameter *tmpparameterlocallist2;
 listparameter *tmpparameterlocallist;
 listmodule *listmoduleinfile;
 listnamelist *listenamelist;
 listnom *NewModuleList;
 listnom *listofmodules;
 listnom *modulelist;
 listnom *Listofvariableinagriffunction;
 listnom *listofsubroutinewhereagrifisused;
 listallocate *AllocateList;
 listvarpointtovar *Listofvarpointtovar; 
                           /*  variables which are pointed to an other one    */
 listindice *Listofavailableindices; 
                           /* List of available indices in the tabvars table  */
 int indeclarationvar;
 int inmodulemeet;
 int incalldeclare;
 int Did_filetoparse_treated;
 int aftercontainsdeclare; /* Signale si l'on vient d'un contains ou non */
 int colnum;
 int callagrifinitgrids;
 int callmpiinit;
 int firstpass;
 int listofvarofusemodulecreated;
 int couldaddvariable;
 int Savemeet;
 int pointedvar;
 int agrif_parentcall;
 int didvariableadded;
 int infunctiondeclare;
 int SubloopScalar;        /* =1 we should put in argument of sub_loop        */
                           /*    only                                         */
                           /*    scalar and not table u(1,1,1) in place of u  */
 int checkexistcommon;
 int insubroutinedeclare;
 int tmpdeclaration_everdone;
 char meetagrifinitgrids[LONGNOM];
 char meetmpiinit[LONGNOM];
 char mpiinitvar[LONGNOM];
 int paramdeclaration_everdone;
 int inmoduledeclare;
 int IntegerIShouldBeAdd;
 int AllocEmpty;
 int dimsempty;
 char *NameTamponfile;
 char toprintglob[LONGNOM];
 char recorddimension[LONGNOM];
 char tmpvargridname[LONGLIGNE];
 char curdimchar[10];
 char OriginalFileName[LONGNOM]; /* Name of the parsing file*/ 
 char EmptyChar[LONGNOM];        /* An empty char */ 
 char commonfile[LONGNOM];
 char curfilename[LONGNOM];
 char commonfile_main[LONGNOM];  /* name of the common file */
 char nomfileoutput[LONGNOM];
 char curbuf[100*LONGNOM];
 char motparse[LONGNOM];
 char charusemodule[LONGNOM];
 char subofagrifinitgrids[LONGNOM];
 char motparse1[LONGNOM];
 char curfile[LONGNOM];         /* name of the current file */
 char mainfile[LONGNOM];        /* name of the configuration file */
 char nomdir[LONGNOM];          /* name of the directory where include files are put */
 char commondirout[LONGNOM];    /* name of the directory where comon files are put */
 char commondirin[LONGNOM];     /* name of the directory containing the common files */
 char filetoparse[LONGNOM];     /* name of the file where all the module file are listed */ 

 FILE *fortranout; /* Output File */
 FILE *fortranin; /* Input File */
 FILE *oldfortranout;
 FILE *subloop;
 FILE *commontomoduleout;
 FILE *paramtomoduleout;
 FILE *inputmodule;
 FILE *allocationagrif;

 long int pos_cur;         /* current position in the output file             */
 long int pos_curagrifparent;
                           /* current position in the output file             */
 long int pos_curcall;     /* current position in the output file             */
 long int pos_curuse;      /* current position in the output file             */
 long int pos_cur_decl;    /* current position in the output file             */
 long int pos_curdata;     /* current position in the output file             */
 long int pos_curparameter;/* current position in the output file             */
 long int pos_curcommon;   /* current position in the output file             */
 long int pos_curinit;     /* current position in the output file             */
 long int pos_curinclude;  /* final position of a line in file                */
 long int pos_end;         /* final position of a line in file                */


/******************************************************************************/
/*********** Declaration of externals subroutines *****************************/
/***************************************************** ************************/

/******************************************************************************/
/*********** UtilNotGridDep.c *************************************************/
/******************************************************************************/
extern void ajoutenotgriddep (char *name);
extern void RemoveNotgriddependFromGlobliste();
extern int VarIsNonGridDepend(char *name);
extern void DECL_0_NonGridDepDeclaration(listvar *listtomodify);
/******************************************************************************/
/*********** WriteInFile.c ****************************************************/
/******************************************************************************/
extern void tofich_reste (FILE * filout, char *s,int returnlineornot);
extern void tofich (FILE * filout, char *s,int returnlineornot);
extern void tofich_blanc (FILE * filout, int size);
extern void RemoveWordCUR(FILE * filout, long int position, 
                                         long int sizetoremove);
extern void RemoveWordSET(FILE * filout, long int position, 
                                         long int sizetoremove);
/******************************************************************************/
/*********** Writedeclarations.c **********************************************/
/******************************************************************************/
extern void WriteBeginDeclaration(variable *v,char ligne[LONGLIGNE]);
extern void WriteScalarDeclaration(variable *v,char ligne[LONGLIGNE]);
extern void WriteTableDeclaration(variable * v,char ligne[LONGLIGNE],int tmpok);
extern void ModifTableDeclaration(variable * v,char ligne[LONGLIGNE]);
extern void writevardeclaration (listvar * var_record, FILE *fileout);
extern void NonGridDepDeclaration(listvar * deb_common);
extern void writedeclaration (listvar * deb_common, FILE *fileout, 
                                                    listvar *presentinthislist);
extern void writesub_loopdeclaration (listvar * deb_common, FILE *fileout);
extern void writedeclarationintoamr (listvar * deb_common, FILE *fileout,
                                    listvar *listin , char commonname[LONGNOM]);
extern void  writedeclarationsubroutinedeclaration(listvar * deb_common,
                                                 FILE *fileout,listvar *listin);
/******************************************************************************/
/*********** WorkWithvarofsubroutineliste.c ***********************************/
/******************************************************************************/
extern void CleanThelistvarofsubroutineliste();
extern void UpdatevarofsubroutinelisteWithcommonlist();
extern void OPTI_1_ajoutvarofsubroutine(listvar *listtoadd);
extern void COM_1_UpdatevarsubroutineWithvarofsubroutinelist();
/******************************************************************************/
/*********** toamr.c **********************************************************/
/******************************************************************************/
extern char *variablenameroottabvars (variable * var);
extern char *variablenametabvars (variable * var, int iorindice);
extern char *variablecurgridtabvars (variable * var,int ParentOrCurgrid);
extern char *vargridnametabvars (variable * var,int iorindice);
extern char *vargridcurgridtabvars (variable * var,int ParentOrCurgrid);
extern char *vargridparam (variable * v, int whichone);
extern void write_probdimagrif_file();
extern void write_includeagrif_file();
extern void write_keysagrif_file();
extern void write_clusteringagrif_file();
extern void write_modtypeagrif_file();
extern void write_createvarnameagrif_file(variable *v,FILE *createvarname,
                                                                int *InitEmpty);
extern void write_Setnumberofcells_file();
extern void write_Getnumberofcells_file();
extern void write_initialisationsagrif_file(variable *v,FILE *initproc,
                                     int *VarnameEmpty);
extern listnom *write_allocation(listvar *newvar,variable *v,
                          listnom *listedesnoms,
                          FILE *alloccalls,
                          FILE *instanceUSE,
                          FILE *modulealloc,
                          int *IndiceMax);
extern void creefichieramr (char *NameTampon);
/******************************************************************************/
/*********** dependfile.c *****************************************************/
/******************************************************************************/
extern void Writethedependnbxnbyfile_fromgloliste();
extern void Readthedependnbxnbyfile_fromgloliste();
extern void Writethedependlistofmoduleused(char *NameTampon );
extern void Readthedependlistofmoduleused(char *NameTampon);
extern void WritedependParameterList(char *NameTampon );
extern listparameter *ReaddependParameterList(char *NameTampon, 
                                                        listparameter *listout);
extern void Writethedependfile(char *NameTampon, listvar *input );
extern void Recordtmplocallist( char *NameTampon);
extern listvar *Recordglobalvarofusefile( char *NameTampon , listvar *listout);
extern void Writethedependavailablefile();
extern void Readthedependavailablefile();
extern int Did_filetoparse_readed(char *NameTampon);
/******************************************************************************/
/*********** SubLoopCreation.c ************************************************/
/******************************************************************************/
extern void OPTI_0_writeheadnewsubforsub();
extern void OPTI_0_writeheadnewsubforfunc();
extern void OPTI_0_writesubroutinedeclaration(listvar *listtomodify);
extern void  WriteVariablelist_subloop(FILE *outputfile);
extern void  WriteVariablelist_subloop_Call(FILE *outputfile);
extern void  WriteVariablelist_subloop_Def(FILE *outputfile);
extern void  WriteHeadofSubroutineLoop();
extern void OPTI_0_closeandcallsubloopandincludeit(int suborfun, 
                                   char endsub[LONGNOM], char optname[LONGNOM]);
/******************************************************************************/
/*********** WorkWithglobliste.c **********************************************/
/******************************************************************************/
extern void CompareNewparsingandoldone();
extern void decl_1_ajoutevar(listvar *listtoadd);
extern void decl_1_ajoutevarsave(listvar *listtoadd);
extern void UpdateIndiceTabvarsofGlobliste();
extern void UpdateIndiceTabvarsofGloblisteFromCommon();
extern void COM_1_UpdateGloblisteWithcommonlist();
/******************************************************************************/
/*********** WorkWithlistvarindoloop.c ****************************************/
/******************************************************************************/
extern void OPTI_1_cleanlistvarfordoloop(int endsuborfunc);
extern void OPTI_1_ajoutevarindoloop(char *ident);
extern void ajoutevarindoloop (char *name);
extern void ajoutevarindoloop_definedimension (char *name);
extern void CleanFromThelistvarindoloopTheAgrifSubArguments();
extern void CleanThelistvarindoloop ();
extern void ModifyThelistvarindoloop();
extern void CompleteThelistvarindoloop();
/******************************************************************************/
/*********** test.c ***********************************************************/
/******************************************************************************/
extern int tests_entrees();
/******************************************************************************/
/*********** WorkWithlistdatavariable.c ***************************************/
/******************************************************************************/
extern void DATA_n_CompleteDataList (char *name,char *values);
extern void DATA_1_CompleteGlobListeWithDatalist ();
/******************************************************************************/
/*********** UtilAgrif.c ******************************************************/
/******************************************************************************/
extern int AGRIF_n_Vartonumber(char *tokname);
extern int AGRIF_n_Agrif_in_Tok_NAME(char *tokname);
extern void AGRIF_1_completeListofvariableinagriffunction(char *ident);
extern void AGRIF_0_ModifyTheVariableName(char *ident);
extern void AGRIF_n_AddsubroutineTolistsubwhereagrifused();
extern void AGRIF_n_AddUseAgrifUtil();
extern void AGRIF_0_NotifyAgrifFunction(char *ident);
extern void AGRIF_0_ModifyTheAgrifFunction(char *ident);
extern void AGRIF_0_AgriffunctionModify(char *ident,int whichone);
extern void AGRIF_0_AddUseAgrifInModuleDeclaration();
/******************************************************************************/
/*********** WorkWithParameterlist.c ******************************************/
/******************************************************************************/
extern void COM_1_AddvartoParamlist(listvar *listin);
extern void COM_1_UpdateparameterlistWithlistvarindoloop();
/******************************************************************************/
/*********** WorkWithAllocatelist.c *******************************************/
/******************************************************************************/
extern void OPTI_1_AddIdentToTheAllocateList(char *nom);
extern int OPTI_0_IsAllocateInThisSubroutine();
extern int OPTI_0_IsVarAllocatable(char *ident);
extern int OPTI_0_varisallocatable(char *ident);
/******************************************************************************/
/*********** UtilCharacter.c **************************************************/
/******************************************************************************/
extern void FindAndChangeNameToTabvars(char name[LONGNOM],
                char toprint[LONGNOM],listvar * listtosee, int ParentOrCurgrid);
extern char *ChangeTheInitalvaluebyTabvarsName(char *nom,listvar *listtoread,
                                                                  int whichone);
extern void IsVarInUseFile(char *nom);
extern listnom *DecomposeTheNameinlistnom(char *nom, listnom * listout);
extern void DecomposeTheName(char *nom);
/******************************************************************************/
/*********** UtilListe.c ******************************************************/
/******************************************************************************/
extern listvar *AddListvarToListvar(listvar *l,listvar *glob,
                                                            int ValueFirstpass);
extern void CreateAndFillin_Curvar(char *type,char *tokname,
                                                listdim *dims,variable *curvar);
extern listvar *duplicatelistvar(listvar * orig);
extern listdim *insertdim(listdim *lin,typedim nom);
extern listdim *reverse(listdim *lin);
extern void change_dim_char(listdim *lin,listvar * l);
extern int num_dims(listdim *d);
extern variable *createvar(char *nom,listdim *d);
extern listvar *insertvar(listvar *lin,variable *v);
extern listvar *settype(char *nom,listvar *lin);
/******************************************************************************/
/*********** UtilFile.c *******************************************************/
/******************************************************************************/
extern FILE * associate (char *filename);
extern FILE * associateaplus (char *filename);
extern long int setposcur();
extern long int setposcurinoldfortranout();
extern void decl_0_modifdeclarationssave(listvar *listtomodify);
extern void OPTI_0_copyuse(char *namemodule);
extern void OPTI_0_copyuseonly(char *namemodule);
/******************************************************************************/
/*********** WorkWithlistofmodulebysubroutine.c *******************************/
/******************************************************************************/
extern void RecordUseModulesVariables();
extern void  RecordUseModulesUseModulesVariables();
extern void Addmoduletothelist(char *name);
extern void WriteUsemoduleDeclaration();
/******************************************************************************/
/*********** WorkWithlistmoduleinfile.c ***************************************/
/******************************************************************************/
extern void MOD_1_FillInlistmodule();
extern void MOD_1_InstanceShouldMadeTo0InModule();
extern void MOD_1_InstanceShouldMadeTo1InModule();
extern void MOD_1_InstanceTo1InModule();
extern int MOD_n_InstanceShouldMadeInModule();
extern int MOD_n_InstanceInModule();
/******************************************************************************/
/*********** UtilFortran.c ****************************************************/
/******************************************************************************/
extern void initdimprob(int dimprobmod, char * nx, char * ny,char* nz);
extern int Variableshouldberemove(char *nom);
extern int variableisglobal(listvar *curvar, listvar *listin);
extern int variableisparameterglobal(listvar *curvar, listparameter *listin);
extern void addi_0_addsubroutine_inst_back_alloc(int moduleorcontains);
extern int OPTI_0_IsTabvarsUseInArgument();
extern int ImplicitNoneInSubroutine();
extern int OPTI_0_varispointer(char *ident);
/******************************************************************************/
/*********** DiversListe.c ****************************************************/
/******************************************************************************/
extern void COM_1_AddCommonvartolist();
extern listnom *Addtolistnom(char *nom, listnom *listin);
extern listname *Add_listname(char *nom,listname *input);
extern void Add_ModuleTo_listofmodules(char *nom);
extern int ModuleIsDefineInInputFile(char *name);
extern void AddNameToListNamelist(char * name);
extern void Addmoduletothelisttmp(char *name);
extern void Add_ModuleTo_Modulelist(char *nom);
extern void OPTI_1_completelistvarpointtovar(char *namemodule,
                                                            listcouple *couple);
extern void Addincludetothelist(char *name);
extern void WriteIncludeDeclaration();
/******************************************************************************/
extern void processfortran(char *fichier_entree);
