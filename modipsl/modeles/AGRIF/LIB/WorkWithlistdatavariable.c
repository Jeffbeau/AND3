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
/*                          DATA_n_COMPLETEDATALIST                           */
/******************************************************************************/
/* This subroutine is used to add a record to listdatavariable                */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ data +--->+ data +--->+ data +--->+  data+             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/******************************************************************************/
void DATA_n_CompleteDataList (char *name,char *values)
{
  listvar *newvar;
  
  newvar=(listvar *)malloc(sizeof(listvar));
  newvar->var=(variable *)malloc(sizeof(variable));
  strcpy(newvar->var->nomvar,name);  
  strcpy(newvar->var->initialvalue,values); 
  newvar->suiv = NULL; 

  if ( !listdatavariable )
  {
     listdatavariable  = newvar ;
  }
  else
  {
     newvar->suiv = listdatavariable;
     listdatavariable = newvar;
  }
}


/******************************************************************************/
/*                    DATA_1_COMPLETEGLOBLISTEWITHDATALIST                    */
/******************************************************************************/
/* This subroutine is used to complete the variable initialisation            */
/* in the globliste with the listdatavariable                                 */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + data +--->+ data +--->+ data +--->+ data +--->+ data +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                  ||                                        */
/*                                  ||                                        */
/*                                  ||                                        */
/*                                  ||                                        */
/*                            initialvalue                                    */
/*                                  ||                                        */
/*                                  ||                                        */
/*        _______     _______     __\/___     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + glob +--->+ glob +--->+ glob +--->+ glob +--->+ glob +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/******************************************************************************/
void DATA_1_CompleteGlobListeWithDatalist ()
{
  listvar *newvar;
  listvar *globlistetmp;
  int out;
  
  /* We are looking for each variable of the listdatavariable where           */
  /* are they located in the globliste                                        */
  newvar = listdatavariable;
  while ( newvar )
  {
     globlistetmp = globliste;
     out = 0 ;
     while( globlistetmp && out == 0 )
     {
        if ( !strcasecmp(newvar->var->nomvar,globlistetmp->var->nomvar) )
	{
           out = 1;
	   if ( strcasecmp(globlistetmp->var->initialvalue,"") )
	   {
              printf("The variable %s has ever a initial value \n"
                                                         , newvar->var->nomvar);
              printf("Error in the CompleteGlobListeWithDatalist routine \n");
              exit(0);	      
	   }
	   else
	   {
              strcpy(globlistetmp->var->initialvalue,newvar->var->initialvalue);
	   }
	}
	else
	{
        globlistetmp = globlistetmp->suiv;
	}
     }
     if ( !globlistetmp )
     {
        printf("The variable %s has not bee found in the globliste \n"
                                                         , newvar->var->nomvar);
        printf("Error in the CompleteGlobListeWithDatalist routine \n");
	exit(0);
     }
     newvar = newvar->suiv;
  }
}

