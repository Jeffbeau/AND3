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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "decl.h"

/*Tests de correction*/

int tests_entrees() 
{
   int erreur=0;

   if (onlyfixedgrids == 0)
   {
      if(regridding==0)
      {
         printf("ERROR: Regridding interval is missing or equals 0\n");
         erreur=1;
      }
  
      if (userefficiency==1)  
      {
         if ((efficiency>100) || (efficiency < 0)) 
         {	 
            printf("ERROR : invalid efficiency (must be between 0 and 100) \n");
            erreur=1;
         }
      }
      else 
      {
         printf("WARNING : clustering efficiency not given: set to 70%% \n");
         efficiency=70;
      }  
   }
  

   if(erreur) return 1;
 
   return erreur;
}
