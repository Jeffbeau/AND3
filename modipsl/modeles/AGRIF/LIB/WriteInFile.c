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
/*                            tofich_reste                                    */
/******************************************************************************/
/* This subroutine is used to write the string s into the fileout             */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void tofich_reste (FILE * filout, char *s,int returnlineornot)	
{
  char temp[61];
  char *tmp;
  int size;


  if (strlen (s) <= 60)
    {
      if ( returnlineornot == 0 ) fprintf (filout, "     & %s", s);
      else fprintf (filout, "     & %s\n", s);
      if ( returnlineornot == 0 ) colnum=colnum+strlen(s)+6;
      else colnum=0;
    }
  else
    {
      strncpy (temp, s, 60);
      strcpy (&temp[60], "\0");

      tmp = strrchr(temp, '+');

      if ( !tmp || strlen(tmp) == 60 ) tmp = strrchr(temp, '-');
      if ( !tmp || strlen(tmp) == 60  ) tmp = strrchr(temp, '/');
      if ( !tmp || strlen(tmp) == 60  ) tmp = strrchr(temp, '*');
      if ( !tmp || strlen(tmp) == 60  ) tmp = strrchr(temp, '%');
      if ( !tmp || strlen(tmp) == 60  ) tmp = strrchr(temp, ':');
      if ( !tmp || strlen(tmp) == 60  ) tmp = strrchr(temp, ')');
      if ( tmp )
      {
         size = strlen(tmp);
      }
      else
      {
         size = 0 ;
      }

      strcpy (&temp[60-size], "\0");

      if ( fortran77 == 0 ) fprintf (filout, "     & %s  &\n", temp);
      else fprintf (filout, "     & %s  \n", temp);
      colnum=0;
      tofich_reste (filout, (char *) &s[60-size],returnlineornot);
    }
}

/******************************************************************************/
/*                            tofich                                          */
/******************************************************************************/
/* This subroutine is used to write the string s into the fileout             */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void tofich (FILE * filout, char *s, int returnlineornot)	
{
  char temp[61];
  char *tmp;
  int size;
  
  if (strlen (s) <= 60)
    {
      if ( returnlineornot == 0 ) fprintf (filout, "      %s", s);
      else fprintf (filout, "      %s\n", s);
      if ( returnlineornot == 0 ) colnum=colnum+strlen(s)+6;
      else colnum=0;
    }
  else
    {
      strncpy (temp, s, 60);
      strcpy (&temp[60], "\0");

      tmp = strrchr(temp, ',');
      if ( !tmp  || strlen(tmp) == strlen(temp) ) tmp = strrchr(temp, '+');
      if ( !tmp  || strlen(tmp) == strlen(temp)  ) tmp = strrchr(temp, '-');
      if ( !tmp  || strlen(tmp) == strlen(temp)  ) tmp = strrchr(temp, '/');
      if ( !tmp  || strlen(tmp) == strlen(temp)  ) tmp = strrchr(temp, '*');
      if ( !tmp  || strlen(tmp) == strlen(temp)  ) tmp = strrchr(temp, '%');
      if ( !tmp  || strlen(tmp) == strlen(temp)  ) tmp = strrchr(temp, ':');
      if ( !tmp  || strlen(tmp) == strlen(temp)  ) tmp = strrchr(temp, ')');
      if ( tmp )
      {
         size = strlen(tmp);
      }
      else
      {
         size = 0 ;
      }

      strcpy (&temp[60-size], "\0");

      if ( fortran77 == 0 ) fprintf (filout, "      %s  &\n", temp);
      else fprintf (filout, "      %s  \n", temp);
      colnum=0;
      tofich_reste (filout, (char *) &s[60-size], returnlineornot);
    }
}

/******************************************************************************/
/*                       tofich_blanc                                         */
/******************************************************************************/
/* This subroutine is used to write size blank into the fileout               */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void tofich_blanc (FILE * filout, int size)
{
  int i;
  
  if (size <= 65)
    {
      fprintf (filout, "%*s\n",size,EmptyChar);
    }
  else
    {
      i=0;
      do
      {
         fprintf (filout, "%*s\n",65,EmptyChar);
         i = i+1;
      } while ( i <= size / 65 );
         fprintf (filout, "%*s\n",size%65,EmptyChar);
    }
}


/******************************************************************************/
/*                         RemoveWord                                         */
/******************************************************************************/
/* This subroutine is used to remove a sentence in the file filout            */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void RemoveWordSET(FILE * filout, long int position, long int sizetoremove)
{
   fseek(filout,position,SEEK_SET);
   tofich_blanc(filout,sizetoremove);
}


/******************************************************************************/
/*                         RemoveWord                                         */
/******************************************************************************/
/* This subroutine is used to remove a sentence in the file filout            */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void RemoveWordCUR(FILE * filout, long int position, long int sizetoremove)
{
   fseek(filout,position,SEEK_CUR);
   tofich_blanc(filout,sizetoremove);
}
