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
%s character
%{
#include <math.h>
#include <stdlib.h>
#include <string.h>
int line_num=1;
extern FILE * yyin;
#define MAX_INCLUDE_DEPTH 30
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
%}

COMMENT "%"  
SEPARATEUR "::"
NIMPORTEQUOI .
COMMENTAIRES1 {COMMENT}{NIMPORTEQUOI}*{COMMENT}
PROBTYPE "1D"|"2D"|"3D"
USEITEM "FIXED_GRIDS"|"ONLY_FIXED_GRIDS"|"DEBUG"
NAME [a-zA-Z\_][a-zA-Z0-9\_]*
DIGIT [0-9]+
NUM {DIGIT}
NEXTLINE \n+[ \t]+"$"|\n+[ \t]+"&"
FILENAME {NAME}"."{NAME}
%%
regridding  return TOK_REGRIDDING; /* period of regridding                    */
coeffrefx   return TOK_COEFFRAFX;  /* space refinement in the x direction     */
coeffrefy   return TOK_COEFFRAFY;  /* space refinement in the y direction     */
coeffrefz   return TOK_COEFFRAFZ;  /* space refinement in the z direction     */
coeffreftx  return TOK_COEFFRAFTX; /* time refinement in the x direction      */
coeffrefty  return TOK_COEFFRAFTY; /* time refinement in the y direction      */
coeffreftz  return TOK_COEFFRAFTZ; /* time refinement in the z direction      */
parammodule return TOK_MODULEMAIN; /* name of the module                      */
efficiency  return TOK_EFFICIENCY; /* efficiency for the adaptive refinement  */
rafmax      return TOK_RAFMAX;     /* minimum size in all directions          */
rafmaxx     return TOK_RAFMAXX;    /* minimum size in x direction             */
rafmaxy     return TOK_RAFMAXY;    /* minimum size in y direction             */
rafmaxz     return TOK_RAFMAXZ;    /* minimum size in z direction             */
notgriddep  return TOK_NOTGRIDDEP; /* variable which are not grid dependent   */
use         return TOK_USE;
minwidth    return TOK_MINWIDTH;   /* minimum width of rectangles for the     */
                                   /*    adaptive refinement                  */
{COMMENTAIRES1}    {}
{SEPARATEUR}        return TOK_SEP;
{FILENAME}         {strcpy(yylval.na,yytext); return TOK_FILENAME;}
{USEITEM}          {strcpy(yylval.na,yytext); return TOK_USEITEM;}
{PROBTYPE}         {strcpy(yylval.na,yytext); return TOK_PROBTYPE;}
                                   /* dimension of the problem                */
{NAME}             {strcpy(yylval.na,yytext); return TOK_NAME;}
{NUM}              {yylval.ival=atoi(yytext); return TOK_NUM;}
;|\,|\(|\)|:|\[|\] {return (int) *yytext;}
\n                 {line_num++;return (int) *yytext;}
[ \t]+ ;
%%


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
