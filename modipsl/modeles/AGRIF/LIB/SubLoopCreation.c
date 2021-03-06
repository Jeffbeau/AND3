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
/*     preparation and write of the argument list of a subroutine             */
/******************************************************************************/
 
 
/******************************************************************************/
/*                        OPTI_0_writeheadnewsubforsub                        */
/******************************************************************************/
/* Firstpass 0                                                                */
/* We should write the head of the subroutine sub_loop_<subroutinename>       */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_0_writeheadnewsubforsub()
{
   int out;
   listusemodule *newmodule;
   
   if ( firstpass == 0 && OPTI_0_IsTabvarsUseInArgument() == 1 ) 
   {

      /* we should add the use agrif_util if it is necessary                  */
      newmodule = listofmodulebysubroutine;
      out = 0 ;
      if ( adduseagrifutil != 1 && inmodulemeet == 0 )
      {
         while ( newmodule && out == 0 )
         {
             if ( !strcasecmp(newmodule->cursubroutine,subroutinename) ||
                  !strcasecmp(newmodule->cursubroutine," ")
                 )
             {
                if ( !strcasecmp(newmodule->charusemodule,"Agrif_Util") ) 
                                                                        out = 1;
             }
             newmodule = newmodule ->suiv;  
         }
         if ( out == 0 ) tofich(fortranout,
                             "\n      USE Agrif_Types, ONLY : Agrif_tabvars",1);
      }
      /* we should modify the name of the variable in the                     */
      /* listvarindoloop which has been declared by the way of the            */
      /* USE fortran  function : USE MOD, U => V                              */
      ModifyThelistvarindoloop();
      /* And write the head of the new subroutine                             */
      WriteHeadofSubroutineLoop(); 
      if ( OPTI_0_IsAllocateInThisSubroutine() == 1 && inmodulemeet == 0 )
      {
         tofich(fortranout,"\n      USE Agrif_Types, ONLY : Agrif_tabvars",1); 
      }
      else if ( out == 1 && fortran77 == 1 ) tofich(fortranout,
                                                    "\n      USE Agrif_Util",1);
      WriteUsemoduleDeclaration();
      WriteIncludeDeclaration();
      if ( ImplicitNoneInSubroutine() == 1 ) fprintf(fortranout,
                                                       "      IMPLICIT NONE\n");
      /* We should write once the declaration of tables (extract              */
      /*    from pointer) in the new subroutine                               */
      tmpdeclaration_everdone = 1;
      if ( todebug == 1 ) fprintf(fortranout,"!!! 111111111111111 \n");
      if ( fortran77 == 1 ) writesub_loopdeclaration
                                                     (parameterlist,fortranout);
      if ( todebug == 1 ) fprintf(fortranout,"!!! 222222222222222 \n");
      writesub_loopdeclaration(listvarindoloop,fortranout);
      if ( todebug == 1 ) fprintf(fortranout,"!!! 333333333333333 \n");
      if ( fortran77 == 1 ) writesub_loopdeclaration
                                              (varofsubroutineliste,fortranout);
      if ( todebug == 1 ) fprintf(fortranout,"!!! 444444444444444 \n");
      if ( fortran77 == 1 ) writesub_loopdeclaration
                                                     (varsubroutine,fortranout);
      if ( todebug == 1 ) fprintf(fortranout,"!!! 555555555555555 \n");
      /* now we should write the function declaration                         */
      /*    case if it is the                                                 */
      writedeclaration (functionlistvar, fortranout,varofsubroutineliste);
      if ( todebug == 1 ) fprintf(fortranout,"!!! 666666666666666 \n");
   }
   else if ( firstpass == 0 )
   {
      WriteUsemoduleDeclaration();
      WriteIncludeDeclaration();
   }
   
}


/******************************************************************************/
/*                        OPTI_0_writeheadnewsubforfunc                       */
/******************************************************************************/
/* Firstpass 0                                                                */
/* We should write the head of the subroutine sub_loop_<subroutinename>       */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_0_writeheadnewsubforfunc()
{
   int out;
   listusemodule *newmodule;

   if ( firstpass == 0 && OPTI_0_IsTabvarsUseInArgument() == 1 )
   {
      /* we should add the use agrif_util if it is necessary                  */
      newmodule = listofmodulebysubroutine;
      out = 0 ;
      if ( adduseagrifutil != 1 && inmodulemeet == 0 )
      {
         while ( newmodule && out == 0 )
         {
             if ( !strcasecmp(newmodule->cursubroutine,subroutinename) ||
                  !strcasecmp(newmodule->cursubroutine," ")
                 )
             {
                if ( !strcasecmp(newmodule->charusemodule,"Agrif_Util") ) 
                                                                        out = 1;
             }
             newmodule = newmodule ->suiv;  
         }
         if ( out == 0 ) tofich(fortranout,
                             "\n      USE Agrif_Types, ONLY : Agrif_tabvars",1);
      }
      WriteHeadofSubroutineLoop();   
      if ( OPTI_0_IsAllocateInThisSubroutine() == 1 && inmodulemeet == 0 )
           tofich(fortranout,"\n      USE Agrif_Types, ONLY : Agrif_tabvars",1);
      WriteUsemoduleDeclaration();
      WriteIncludeDeclaration();
      if ( ImplicitNoneInSubroutine() == 1 ) fprintf(fortranout,
                                                       "      IMPLICIT NONE\n");
      /* We should write once the declaration of tables (extract              */
      /*    from pointer) in the new subroutine                               */
      tmpdeclaration_everdone = 1;
      if ( fortran77 == 1 ) writesub_loopdeclaration
                                                     (parameterlist,fortranout);
      writesub_loopdeclaration(listvarindoloop,fortranout);
      if ( fortran77 == 1 ) writesub_loopdeclaration
                                                     (varsubroutine,fortranout);
      if ( fortran77 == 1 ) writesub_loopdeclaration
                                              (varofsubroutineliste,fortranout);
      /* now we should write the function declaration                         */
      /*    case if it is the                                                 */
      writedeclaration (functionlistvar, fortranout,varofsubroutineliste);
   }
   else if ( firstpass == 0 ) 
   {
      WriteUsemoduleDeclaration();
      WriteIncludeDeclaration();
   }
}


/******************************************************************************/
/*                    OPTI_0_writesubroutinedeclaration                       */
/******************************************************************************/
/* Firstpass 0                                                                */
/* We should write the declaration of the subroutine in order to              */
/* create the new sub_loop subroutine                                         */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_0_writesubroutinedeclaration(listvar *listtomodify)
{
   if ( VariableIsParameter == 0 && SaveDeclare == 0)
   {
      if (firstpass == 0 && OPTI_0_IsTabvarsUseInArgument() == 1 )
      {
         /* We should write this declaration into the original                */
         /*    subroutine too                                                 */
        if ( fortran77 == 1                 && 
             paramdeclaration_everdone == 0 && 
             varofsubroutineliste 
           )
        {
           paramdeclaration_everdone = 1;
           writedeclarationsubroutinedeclaration
                             (parameterlist,oldfortranout,varofsubroutineliste);
        }
        writedeclaration (listtomodify, oldfortranout,varofsubroutineliste);
      }
   }
}

/******************************************************************************/
/*                    WriteVariablelist_subloop                               */
/******************************************************************************/
/* This subroutine is used to write the list of the variable which            */
/* should be called by the sub_loop_<name> subroutine                         */
/* The first part is composed by the list of the local variables              */
/******************************************************************************/
/*                                                                            */
/*    varofsubroutineliste              a,b,c,  &                             */
/*                                      d,e,f,  &                             */
/*     a,b,c,d,e,f,g,h     ========>    g,h                                   */
/*                                                                            */
/******************************************************************************/
void WriteVariablelist_subloop(FILE *outputfile)
{
   listvar *parcours;   
   char ligne[LONGNOM];
   int compteur;
   
   parcours = varofsubroutineliste;
   didvariableadded = 0;
   compteur = 0 ;
   
   while ( parcours )   
   {
      /* if the readed variable is a variable of the subroutine               */
      /*    subroutinename we should write the name of this variable          */
      /*    in the output file                                                */
      if ( !strcasecmp(parcours->var->modulename,subroutinename) )
      {      
         if ( didvariableadded == 0 )
         {
            strcpy(ligne,"");
         }
         else
         {
            if ( compteur == 0 ) strcpy(ligne,"");
            strcat(ligne,",");
         }
         strcat(ligne,parcours->var->nomvar);
         didvariableadded = 1;
         compteur = compteur + 1;
         if ( compteur == 3 ) 
         {
            if ( fortran77 == 0 ) strcat(ligne," &");
            if ( fortran77 == 0 ) fprintf(outputfile,"\n      %s",ligne);
            else fprintf(outputfile,"\n     & %s",ligne);
            compteur = 0;
         }
      }
      parcours = parcours -> suiv;   
   }
   if ( compteur != 3 && compteur != 0 ) 
   {
      if ( fortran77 == 0 ) fprintf(outputfile,"\n      %s &",ligne);
      else fprintf(outputfile,"\n     & %s ",ligne);
   }
}


/******************************************************************************/
/*                     WriteVariablelist_subloop_Call                         */
/******************************************************************************/
/* This subroutine is used to write the list of the variable which            */
/* should be called by the sub_loop_<name> subroutine into the called         */
/* The second part is composed by the list of the global table                */
/******************************************************************************/
/*                                                                            */
/*      listvarindoloop        SubloopScalar = 0 | SubloopScalar = 1          */
/*                                a,b,c,  &      |  a,b(1,1),c,      &        */
/*     a,b,c,d,e,f,g,h  =====>    d,e,f,  &      |  d(1),e(1,1,1),f, &        */
/*                                g,h            |  g,h(1,1)                  */
/*                                                                            */
/******************************************************************************/
void  WriteVariablelist_subloop_Call(FILE *outputfile)
{
   listvar *parcours;   
   char ligne[LONGNOM*100];
   char ligne2[10];
   int i;
   
   parcours = listvarindoloop;
   while ( parcours )   
   {
      /* if the readed variable is a variable of the subroutine               */
      /*    subroutinename we should write the name of this variable          */
      /*    in the output file                                                */
      if ( !strcasecmp(parcours->var->modulename,subroutinename) )
      {
         if ( didvariableadded == 0 )
         {
            strcpy(ligne,"");
         }
         else
         {
            strcpy(ligne,"");
            strcat(ligne,",");
         }
         strcat(ligne,vargridcurgridtabvars(parcours->var,0));
         /* if it is asked in the call of the conv we should give             */
         /* scalar in argument, so we should put (1,1,1) after the            */
         /* the name of the variable                                          */
         if (  SubloopScalar != 0 && 
               (OPTI_0_IsVarAllocatable(parcours->var->nomvar) == 0 && 
	       parcours->var->pointerdeclare == 0 ) &&
               parcours->var->nbdim != 0 )
         {
             i = 1;
             while ( i <=  parcours->var->nbdim )
             {
                if ( i == 1 ) strcat(ligne,"( ");
                if ( SubloopScalar == 2 )
                {
                   strcat(ligne,":");
                   if ( i != parcours->var->nbdim ) strcat(ligne,",");
                }
                else
                {
                   strcat(ligne," lbound( ");
                   strcat(ligne,vargridcurgridtabvars(parcours->var,0));
                   strcat(ligne,",");
                   strcpy(ligne2,"");
                   sprintf(ligne2,"%d",i);
                   strcat(ligne,ligne2);
                   if ( i != parcours->var->nbdim ) strcat(ligne,"),");
                }
                if ( i == parcours->var->nbdim ) strcat(ligne,"))");
                i++;
             }
         }
         didvariableadded = 1;
         if ( fortran77 == 0 ) strcat(ligne," &");
         if ( fortran77 == 0 ) fprintf(outputfile,"\n");
         else fprintf(outputfile,"\n     & ");
         tofich(outputfile,ligne,0);
      }
      parcours = parcours -> suiv;   
   }
   /* Now we should replace the last ", &" by " &"                            */
   if ( didvariableadded != 0 && fortran77 == 0 ) fseek(outputfile,-1,SEEK_CUR);
   if ( didvariableadded == 0 ) fseek(outputfile,-2,SEEK_CUR);
}


/******************************************************************************/
/*                       WriteVariablelist_subloop_Def                        */
/******************************************************************************/
/* This subroutine is used to write the list of the variable which            */
/* should be called by the sub_loop_<name> subroutine into the def            */
/* The second part is composed by the list of the global table                */
/* <name>_tmp                                                                 */
/******************************************************************************/
/*                                                                            */
/*       listvarindoloop                                                      */
/*                                a-tmp,b-tmp,c_tmp, &                        */
/*     a,b,c,d,e,f,g,h  =====>    d_tmp,e_tmp,f_tmp, &                        */
/*                                g_tmp,h_tmp                                 */
/*                                                                            */
/******************************************************************************/
void  WriteVariablelist_subloop_Def(FILE *outputfile)
{
   listvar *parcours;   
   char ligne[LONGNOM];
   int compteur;

   parcours = listvarindoloop;
   compteur = 0 ;
   while ( parcours )   
   {
      /* if the readed variable is a variable of the subroutine               */
      /*    subrotinename we should write the name of this variable           */
      /*    in the output file                                                */
      if ( !strcasecmp(parcours->var->modulename,subroutinename) )
      {
         if ( didvariableadded == 0 )
         {
            strcpy(ligne,"");
         }
         else
         {
            if ( compteur == 0 ) strcpy(ligne,"");
            strcat(ligne,",");
         }
         strcat(ligne,parcours->var->nomvar);
         compteur = compteur + 1;
         didvariableadded = 1;
         if ( compteur == 3 ) 
         {
            if ( fortran77 == 0 ) strcat(ligne," &");
            if ( fortran77 == 0 ) fprintf(outputfile,"\n      %s",ligne);
            else fprintf(outputfile,"\n     & %s",ligne);
            compteur = 0;
         }
      }
      parcours = parcours -> suiv;   
   }   
   if ( compteur != 3 && compteur != 0 )
   {
      if ( fortran77 == 0 ) fprintf(outputfile,"\n      %s &",ligne);
      else fprintf(outputfile,"\n     & %s",ligne);
   }

   /* Now we should replace the last ", &" by " &"                            */
   if ( didvariableadded != 0 && fortran77 == 0 ) fseek(outputfile,-1,SEEK_CUR);
   if ( didvariableadded == 0 ) fseek(outputfile,-1,SEEK_CUR);
}


/******************************************************************************/
/*                      WriteHeadofSubroutineLoop                             */
/******************************************************************************/
/* This subroutine is used to write the head of the subroutine                */
/* Sub_Loop_<name>                                                            */
/******************************************************************************/
/*                 Sub_loop_subroutine.h                                      */
/*                                                                            */
/*                 subroutine Sub_Loop_subroutine ( &                         */
/*                 a,b,c, &                                                   */
/* SubLoopScalar   d,e(1,1),f(1,1,1), &                                       */
/*                 g,h  &                                                     */
/*                 )                                                          */
/* adduseagrifutil USE Agrif_Util                                             */
/******************************************************************************/
void  WriteHeadofSubroutineLoop()
{
   char ligne[LONGNOM];
   FILE * subloop;


   tofich(fortranout,"\n",1);
   /* Open this newfile                                                       */
   sprintf(ligne,"Sub_Loop_%s.h",subroutinename);
   subloop = associate(ligne);
   /*                                                                         */
   if ( fortran77 == 0 ) sprintf(ligne,"      subroutine Sub_Loop_%s( &"
                                                               ,subroutinename);
   else sprintf(ligne,"      subroutine Sub_Loop_%s( ",subroutinename);
   fprintf(subloop,ligne);
   /*                                                                         */
   WriteVariablelist_subloop(subloop);
   WriteVariablelist_subloop_Def(subloop);
   /*                                                                         */
   sprintf(ligne,")");
   tofich(subloop,ligne,1);   
   /* if USE agrif_Util should be add                                         */
   if ( adduseagrifutil == 1 ) fprintf(subloop,"\n      USE Agrif_Util\n");
   /*                                                                         */
   oldfortranout = fortranout;
   fortranout = subloop;
}

/******************************************************************************/
/*                OPTI_0_closeandcallsubloopandincludeit                      */
/******************************************************************************/
/* Firstpass 0                                                                */
/* We should close the sub_loop subroutine, call it and close the             */
/* function (suborfun = 0)                                                    */
/* subroutine (suborfun = 1)                                                  */
/* end (suborfun = 2)                                                         */
/* end program (suborfun = 3)                                                 */
/* and include the sub_loop subroutine after                                  */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_0_closeandcallsubloopandincludeit(int suborfun, char endsub[LONGNOM],
                                                          char optname[LONGNOM])
{
   char ligne[LONGNOM];

   if ( firstpass == 0 ) 
   {
   if ( OPTI_0_IsTabvarsUseInArgument() == 1 )
   {
      /* We should remove the key word end subroutine                         */
      RemoveWordCUR(fortranout,(long)(-strlen(optname)-strlen(endsub)-1),
                                       strlen(optname)+strlen(endsub)+1);
      /* We should close the loop subroutine                                  */
      sprintf(ligne,"\n      end subroutine Sub_Loop_%s",subroutinename);
      tofich(fortranout,ligne,1);
      fclose(fortranout);  
      fortranout = oldfortranout;
      /* Now we add the call af the new subroutine                            */
      if ( fortran77 == 0 ) sprintf(ligne,"\n      Call Sub_Loop_%s( &"
                                                               ,subroutinename);
      else sprintf(ligne,"\n      Call Sub_Loop_%s( ",subroutinename);
      fprintf(fortranout,ligne);
      /* Write the list of the local variables used in this new subroutine    */
      WriteVariablelist_subloop(fortranout);
      /* Write the list of the global tables used in this new subroutine      */
      /*    in doloop                                                         */
      WriteVariablelist_subloop_Call(fortranout);  
      /* Close the parenthesis of the new subroutine called                   */
      sprintf(ligne,")");
      tofich(fortranout,ligne,1);  
      /* We should close the original subroutine                              */
      if ( !strcasecmp(subofagrifinitgrids,subroutinename) )
      {
/*         fprintf(fortranout,"      CALL Agrif_Deallocation\n");*/
      }
      if ( suborfun == 3 ) sprintf(ligne,"\n      end program %s"
                                                               ,subroutinename);
      if ( suborfun == 2 ) sprintf(ligne,"\n      end");
      if ( suborfun == 1 ) sprintf(ligne,"\n      end subroutine %s"
                                                               ,subroutinename);
      if ( suborfun == 0 ) sprintf(ligne,"\n      end function %s"
                                                               ,subroutinename);
      tofich(fortranout,ligne,1);
      /* we should include the above file in the original code                */
      sprintf(ligne,"\n#include \"Sub_Loop_%s.h\" \n",subroutinename);
      tofich(fortranout,ligne,1);
      }
   }
}
