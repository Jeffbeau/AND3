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
%x parameter
%s character
%{
#include <math.h>
#include <stdlib.h>
#include <string.h>
extern FILE * yyin;
#define MAX_INCLUDE_DEPTH 30
#define tabsize 6
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
int line_num_fortran=1;
int line_num_fortran_common=1;
int newlinef90 = 0;
char *tmp;
/******************************************************************************/
/**************PETITS PB NON PREVUS *******************************************/
/******************************************************************************/
/* NEXTLINF77 un ligne fortran 77 peut commencer par -      &a=b or on        */
/*            a prevu seulement       & a=b avec l'espace entre le symbole    */
/*            de la 7eme et le debut de la ligne de commande                  */
/*            le ! est aussi interdit comme symbole de la 7 eme colonne       */
/*            Normalement NEXTLINEF77 \n+[ ]{5}[^ ]                           */
/******************************************************************************/
#define YY_USER_ACTION \
	{\
	   if (firstpass == 0) \
           {\
              strcat(curbuf,yytext); \
	      strcpy(motparse,yytext);\
              colnum = colnum + strlen(motparse);\
	      /*printf("motparse = %s %d\n",motparse,strlen(motparse));*/\
              ECHO; \
           }\
           strcpy(motparse1,yytext);\
	/*if ( firstpass == 1 ) 
                      printf("yytext = %s %d\n",yytext,strlen(yytext));*/\
	}
%}
AGRIFDEB "Agrif_debut"
AGRIFFIN "Agrif_fin"
NOTTREAT Agrif_do_not_treat
ENDNOTTREAT Agrif_end_do_not_treat

NIMPORTEQUOI .
SAME_LINE ";"
SLASH "/"
DSLASH "/"[ \t]*"/"
NAME [a-zA-Z\_][a-zA-Z0-9\_]*
DIGIT [0-9]+
INT {DIGIT}
EXPONENT e[-+]?{DIGIT}
DEXPONENT d[-+]?{DIGIT}
QEXPONENT q[-+]?{DIGIT}
REAL (({DIGIT}\.[0-9]+|[0-9]*\.{DIGIT}){EXPONENT}?)|{DIGIT}\.{EXPONENT}
REALDP (({DIGIT}\.[0-9]+|[0-9]*\.{DIGIT}){DEXPONENT}?)|{DIGIT}\.{DEXPONENT}
REALQP (({DIGIT}\.[0-9]+|[0-9]*\.{DIGIT}){QEXPONENT}?)|{DIGIT}\.{QEXPONENT}
ENDFUNCTION end[ \t]*function
DOUBLEPRECISION double[ \t]*precision
DOUBLECOMPLEX double[ \t]*complex

COMMENTAIRESFORTRAN77 (^[Cc]{NIMPORTEQUOI}*)
COMMENTAIRESFORTRAN90 ^([ \t]*!{NIMPORTEQUOI}*\n)
COMMENTAIRESFORTRAN90_2 (!{NIMPORTEQUOI}*)
NEXTLINEF90 "&"{NIMPORTEQUOI}*[\n]*
NEXTLINEF77 \n[ \t]*"&"
%%
{NAME}\= {
		if (firstpass == 0)
		{
		fseek(fortranout,-1,1);
		strcpy(&curbuf[strlen(curbuf)-1],"\0");
		}
		yyless(yyleng-1);
		strcpy(yylval.na,yytext);
		return TOK_NAME;
		}
^C${NOTTREAT}           {return TOK_DONOTTREAT;}
^C${ENDNOTTREAT}        {return TOK_ENDDONOTTREAT;}
^C${AGRIFDEB}            return TOK_DEBUT;
^C${AGRIFFIN}            return TOK_FIN;
^C$OMP[ \t]*{NIMPORTEQUOI}* return TOK_OMP;
^C$[ \t]*{NIMPORTEQUOI}* return TOK_DOLLAR;

subroutine              {return TOK_SUBROUTINE;}
program                 {return TOK_PROGRAM;}
allocate                {return TOK_ALLOCATE;}
deallocate              {return TOK_DEALLOCATE;}
result                  {return TOK_RESULT;}
function                {return TOK_FUNCTION;}
end[ \t]*subroutine     {strcpy(yylval.na,yytext);return TOK_ENDSUBROUTINE;}
end[ \t]*program        {strcpy(yylval.na,yytext);return TOK_ENDPROGRAM;}
end[ \t]*function       {strcpy(yylval.na,yytext);return TOK_ENDFUNCTION;}
end                     {strcpy(yylval.na,yytext);return TOK_ENDUNIT;}
include                  return TOK_INCLUDE;
use                     {return TOK_USE;}
rewind                  {return TOK_REWIND;}
implicit                 return TOK_IMPLICIT;
none                     return TOK_NONE;
call                     return TOK_CALL;
.true.                   return TOK_TRUE;
.false.                  return TOK_FALSE;
\=\>                    {return TOK_POINT_TO;}
\*\*                    {strcpy(yylval.na,yytext);return TOK_DASTER;}
\.[ \t]*eq\.            {strcpy(yylval.na,yytext);return TOK_EQ;}
\.[ \t]*gt\.            {strcpy(yylval.na,yytext);return TOK_GT;}
\.[ \t]*ge\.            {strcpy(yylval.na,yytext);return TOK_GE;}
\.[ \t]*lt\.            {strcpy(yylval.na,yytext);return TOK_LT;}
\.[ \t]*le\.            {strcpy(yylval.na,yytext);return TOK_LE;}
\.[ \t]*ne\.            {strcpy(yylval.na,yytext);return TOK_NE;}
\.[ \t]*not\.           {strcpy(yylval.na,yytext);return TOK_NOT;}
\.[ \t]*or\.            {strcpy(yylval.na,yytext);return TOK_OR;}
\.[ \t]*xor\.           {strcpy(yylval.na,yytext);return TOK_XOR;}
\.[ \t]*and\.           {strcpy(yylval.na,yytext);return TOK_AND;}
module                  {return TOK_MODULE;}
do[ \t]*while           {return TOK_DOWHILE;}
end[ \t]*module          return TOK_ENDMODULE;
end[ \t]*do              return TOK_ENDDO;
do                      {return TOK_PLAINDO;}
real                    {strcpy(yylval.na,yytext);return TOK_REAL;}
integer                 {strcpy(yylval.na,yytext);return TOK_INTEGER;}
logical                 {strcpy(yylval.na,yytext);return TOK_LOGICAL;}
character               {strcpy(yylval.na,yytext);return TOK_CHARACTER;}
allocatable             {return TOK_ALLOCATABLE;}
close                    return TOK_CLOSE;
inquire                  return TOK_INQUIRE;
dimension               {return TOK_DIMENSION;}
pause                    return TOK_PAUSE;
equivalence              return TOK_EQUIVALENCE;
stop                     return TOK_STOP;
where                    return TOK_WHERE;
end[ \t]*where           return TOK_ENDWHERE;
else[ \t]*where          return TOK_ELSEWHERE;
complex                 {return TOK_COMPLEX;}
^[ \t]*contains         {return TOK_CONTAINS;}
only                    {return TOK_ONLY;}
parameter               {return TOK_PARAMETER;}
common                  {return TOK_COMMON;}
external                {return TOK_EXTERNAL;}
intent                  {return TOK_INTENT;}
kind                    {return TOK_KIND;}
pointer                 {return TOK_POINTER;}
optional                {return TOK_OPTIONAL;}
save                    {return TOK_SAVE;}
^[ \t]*type              {return TOK_TYPE;}
end[ \t]*type           {return TOK_ENDTYPE;}
open                     return TOK_OPEN;
return                   return TOK_RETURN;
exit                     return TOK_EXIT;
print                    return TOK_PRINT;
module[ \t]*procedure   {return TOK_PROCEDURE;}
read                    {return TOK_READ;}
namelist                {return TOK_NAMELIST;}
write                   {return TOK_WRITE;}
target                  {return TOK_TARGET;}
public                  {return TOK_PUBLIC;}
private                 {return TOK_PRIVATE;}
in                      {return TOK_IN;}
data                    {return TOK_DATA;} 
continue                 return TOK_CONTINUE;
go[ \t]*to              {return TOK_PLAINGOTO;}
out                     {return TOK_OUT;}
inout                   {return TOK_INOUT;}
intrinsic               {return TOK_INTRINSIC;}
then                    {return TOK_THEN;}
else[ \t]*if            {return TOK_ELSEIF;}
else                    {return TOK_ELSE;}
end[ \t]*if             {return TOK_ENDIF;}
if                      {return TOK_LOGICALIF;}
sum[ \t]*\(              {return TOK_SUM;}
max                     {return TOK_MAX;}
tanh                    {return TOK_TANH;}
maxval                  {return TOK_MAXVAL;}
trim                    {return TOK_TRIM;}
sqrt                    {return TOK_SQRT;}
select[ \t]*case        {return TOK_SELECTCASE;}
case                    {return TOK_CASE;}
case[ \t]*default       {return TOK_CASEDEFAULT;}
end[ \t]*select         {return TOK_ENDSELECT;}
file[ \t]*\=            {return TOK_FILE;}
exist[ \t]*\=           {return TOK_EXIST;}
min[ \t]*\(             {return TOK_MIN;}
int                     {return TOK_INT;}
nint                    {return TOK_NINT;}
float                   {return TOK_FLOAT;}
exp                     {return TOK_EXP;}
cos                     {return TOK_COS;}
cosh                    {return TOK_COSH;}
acos                    {return TOK_ACOS;}
sin                     {return TOK_SIN;}
sinh                    {return TOK_SINH;}
asin                    {return TOK_ASIN;}
log                     {return TOK_LOG;}
tan                     {return TOK_TAN;}
atan                    {return TOK_ATAN;}
abs                     {return TOK_ABS;}
mod                     {return TOK_MOD;}
sign                    {return TOK_SIGN;}
minloc                  {return TOK_MINLOC;}
maxloc                  {return TOK_MAXLOC;}
minval                  {return TOK_MINVAL;}
interface               {return TOK_INTERFACE;}
end[ \t]*interface      {return TOK_ENDINTERFACE;}
\({SLASH}               {return TOK_LEFTAB;}
{SLASH}\)               {return TOK_RIGHTAB;}
format                  {return TOK_FORMAT;}
{DOUBLEPRECISION}       {strcpy(yylval.na,yytext);return TOK_DOUBLEPRECISION;}
{DOUBLECOMPLEX}         {strcpy(yylval.na,yytext);return TOK_DOUBLECOMPLEX;}
{SAME_LINE}             {return '\n';}
{SLASH}                 {strcpy(yylval.na,yytext);return TOK_SLASH;}
DSLASH                  {strcpy(yylval.na,yytext);return TOK_DSLASH;}
(\')[^']*&{0,1}\n[ \t]*&{0,1}[^']*(\')         {strcpy(yylval.na,yytext);return TOK_CHAR_CUT;}
(\')[^\n']*(\')         {strcpy(yylval.na,yytext);return TOK_CHAR_CONSTANT;}
(\")[^\n]*(\")          {strcpy(yylval.na,yytext);return TOK_CHAR_MESSAGE;}
({NAME}{REAL})          {strcpy(yylval.na,yytext);return TOK_CHAR_INT;}
{NAME}                  {strcpy(yylval.na,yytext);return TOK_NAME;}
{REAL}                  {strcpy(yylval.na,yytext);return TOK_CSTREAL;}
{REALDP}                {strcpy(yylval.na,yytext);return TOK_CSTREALDP;}
{REALQP}                {strcpy(yylval.na,yytext);return TOK_CSTREALQP;}
({DIGIT}\.)/[^{NAME}|"and."|"false."|"true."|"eq."|"or."|"gt."|"ge."|"lt."|"le."|"not."|"ne."] {strcpy(yylval.na,yytext);return TOK_CSTREAL;}
\.                      {return TOK_POINT;}
{INT}                   {strcpy(yylval.na,yytext);return TOK_CSTINT;}
\$ {}
\'|\"                   {return TOK_QUOTE;}
;|\(|\)|:|\[|\]|\+|\-|\*|\% {strcpy(yylval.na,yytext);return (int) *yytext;} 
\,                      {return (int) *yytext;}
\=                      {return (int) *yytext;}
\<                      {return (int) *yytext;}
\>                      {return (int) *yytext;}
\n                      {colnum=0;line_num_fortran++;line_num_fortran_common++; return (int) *yytext;}
^[ ]*$
^(((" "|[0-9]){1,5})|([ \t]{1,5}))[ &]+ {if (newlinef90 == 0) return TOK_LABEL; else newlinef90 = 0;}
[ ]+
[\t]+                   {colnum=colnum-1+tabsize;}
[ \t]+ ;
{NEXTLINEF90}           {line_num_fortran++;line_num_fortran_common++;newlinef90=1;colnum=0;}
{NEXTLINEF77}           {line_num_fortran++;line_num_fortran_common++;colnum=0;}
{COMMENTAIRESFORTRAN77} {
                                       tmp =  strstr(motparse1,"contains");
                           if ( !tmp ) tmp =  strstr(motparse1,"CONTAINS");
                           if ( !tmp ) tmp =  strstr(motparse1,"Contains");
                           if (  tmp ) 
                           {
                              if ( strlen(motparse1) == strlen(tmp)+1 ) 
                              {
                                 return TOK_CONTAINS;
                              }
                              else
                              {
                                 colnum=0;line_num_fortran++;line_num_fortran_common++;
                              }
                           }
                           else
                           {
                              colnum=0;line_num_fortran++;line_num_fortran_common++;                           
                           }
                         }
{COMMENTAIRESFORTRAN90}   {
                             if ( !strcasecmp(motparse1,"!$AGRIF_DO_NOT_TREAT\n")) return TOK_DONOTTREAT; 
                             if ( !strcasecmp(motparse1,"!$AGRIF_END_DO_NOT_TREAT\n")) return TOK_ENDDONOTTREAT; 
                          }
{COMMENTAIRESFORTRAN90_2} {
                             if ( !strcasecmp(motparse1,"!$AGRIF_DO_NOT_TREAT\n")) return TOK_DONOTTREAT; 
                             if ( !strcasecmp(motparse1,"!$AGRIF_END_DO_NOT_TREAT\n")) return TOK_ENDDONOTTREAT; 
                          }
%%

fortranerror(char *s)
{
   if (!strcasecmp(curfile,mainfile))
   {
      printf("%s line %d, fichier %s\n",s,line_num_fortran,curfile);
   }
   else
   {
      printf("%s line %d, fichier %s\n",s,line_num_fortran_common,curfile);
   }
}

int fortranwrap()
{
}
