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
/*                            associate                                       */
/******************************************************************************/
/* This subroutine is used to open a file                                     */
/******************************************************************************/
FILE * associate (char *filename)
{
  char filefich[LONGNOM];
  sprintf(filefich,"%s/%s",nomdir,filename);
  return fopen (filefich, "w");
}


/******************************************************************************/
/*                          associateaplus                                    */
/******************************************************************************/
/* This subroutine is used to open a file with option a+                      */
/******************************************************************************/
FILE * associateaplus (char *filename)
{
  char filefich[LONGNOM];
  sprintf(filefich,"%s/%s",nomdir,filename);
  return fopen (filefich, "a+");
}


/******************************************************************************/
/*                           setposcur                                        */
/******************************************************************************/
/* This subroutine is used to know the current position in the file           */
/******************************************************************************/
/*                                                                            */
/*                      setposcur ---------> position in file                 */
/*                                                                            */
/******************************************************************************/
long int setposcur()
{
   fflush(fortranout);
   return ftell(fortranout);
}

/******************************************************************************/
/*                      setposcurinoldfortranout                              */
/******************************************************************************/
/* This subroutine is used to know the position in the oldfortranout         */
/******************************************************************************/
/*                                                                            */
/*             setposcurinoldfortranout ---------> position in file           */
/*                                                                            */
/******************************************************************************/
long int setposcurinoldfortranout()
{
   fflush(oldfortranout);
   return ftell(oldfortranout);
}

/******************************************************************************/
/*                       decl_0_modifdeclarationssave                         */
/******************************************************************************/
/* Firstpass 0                                                                */
/* We should modify this declaration in the file fortranout. case SAVE        */
/******************************************************************************/
void decl_0_modifdeclarationssave(listvar *listtomodify)
{
   if ( VariableIsParameter == 0 && SaveDeclare == 1) 
   {
      if (firstpass == 0)
      {
         pos_end = setposcur();
         RemoveWordSET(fortranout,pos_cur,
                               pos_end-pos_cur);
      }
   }      
}

/******************************************************************************/
/*                    OPTI_0_copyuse                                          */
/******************************************************************************/
/* Firstpass 0                                                                */
/* We should write in the fortranout the USE tok_name                         */
/* read in the original file                                                  */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_0_copyuse(char *namemodule)
{
   if (firstpass == 0 && OPTI_0_IsTabvarsUseInArgument() == 1 )
   {
      /* We should write this declaration into the original subroutine too    */
      fprintf(oldfortranout,"      USE %s \n",namemodule);
   }
}

/******************************************************************************/
/*                    OPTI_0_copyuseonly                                      */
/******************************************************************************/
/* Firstpass 0                                                                */
/* We should write in the fortranout the USE tok_name, only                   */
/* read in the original file                                                  */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void OPTI_0_copyuseonly(char *namemodule)
{
   if (firstpass == 0 && OPTI_0_IsTabvarsUseInArgument() == 1 )
   {
      /* We should write this declaration into the original subroutine too    */
      fprintf(oldfortranout,"      USE %s , ONLY : \n",namemodule);
   }
}
