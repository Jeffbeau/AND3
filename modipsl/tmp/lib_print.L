


PGF90 (Version     10.5)          09/23/2010  18:36:34      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /distrib/local/netcdf/pgf/include/
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: lib_print.f90

(    1) MODULE lib_print
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  lib_print  ***
(    4)    !! print librairy :  formated real and integer array print
(    5)    !!=====================================================================
(    6)       
(    7)    !!----------------------------------------------------------------------
(    8)    !!   prihin       : print an integer 2D horizontal field
(    9)    !!   prihre       : print an real 2D horizontal field
(   10)    !!   prizre       : print an real 2D vertical field
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Mules used
(   13)    USE par_kind      ! kind parameters
(   14) 
(   15)    IMPLICIT NONE
(   16)    PRIVATE
(   17) 
(   18)    !! * Routine accessibility
(   19)    PUBLIC prihin, prihre, prizre
(   20)    !!----------------------------------------------------------------------
(   21)    !!----------------------------------------------------------------------
(   22)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   23)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_print.f90,v 1.2 2005/03/27 18:34:47 opalod Exp $ 
(   24)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   25)    !!----------------------------------------------------------------------
(   26) CONTAINS
(   27) 
(   28)    SUBROUTINE prihin( ktab, ki   , kj   , kideb, kifin ,   &
(   29)                       kind, kjdeb, kjfin, kjnd , kscale, kumout )
(   30)       !!----------------------------------------------------------------------
(   31)       !!                   ***  SUBROUTINE  prihre  ***
(   32)       !!  
(   33)       !! ** purpose :   Print an integer field
(   34)       !!
(   35)       !! ** method :   format of print is selected with the dummy argument kscale
(   36)       !!
(   37)       !! History :
(   38)       !!        !  90-04 (0. Marti)  Original code
(   39)       !!        !  92-02 (M. Imbard)
(   40)       !!        !  03-07 (G. Madec)  F90, free form
(   41)       !!----------------------------------------------------------------------
(   42)       !! * Arguments
(   43)       INTEGER, INTENT( in ) ::   &
(   44)          ki, kj,                 &  ! array dimensions
(   45)          kideb, kifin, kind,     &  ! first and last index, increment for i 
(   46)          kjdeb, kjfin, kjnd,     &  ! first and last index, increment for j
(   47)          kscale,                 &  ! kscale=0 or > 5  print ktab with format i8
(   48)       !                             !         kscale=1 print ktab with format i1
(   49)       !                             !         kscale=2 print ktab with format i2






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 2

(   50)       !                             !         kscale=3 print ktab with format i3
(   51)       !                             !         kscale=4 print ktab with format i4
(   52)       !                             !         kscale=5 print ktab with format i5
(   53)          kumout                     ! unit in which print
(   54)       INTEGER, DIMENSION(ki,kj), INTENT( in ) ::   &
(   55)          ktab                       ! integer 2D array to be print
(   56) 
(   57)       !! * local declarations
(   58)       INTEGER ::   ji, jj, jn       ! dummy loop indices
(   59)       INTEGER ::   isca, il1, il2   ! temporary integers
(   60)       INTEGER ::   iind, ijnd       ! temporary integers
(   61) 
(   62)       isca = 10
(   63)       IF( kscale == 0 )   isca = 10
(   64)       IF( kscale == 1 )   isca = 100
(   65)       IF( kscale == 2 )   isca = 60
(   66)       IF( kscale == 3 )   isca = 40
(   67)       IF( kscale == 4 )   isca = 30
(   68)       IF( kscale == 5 )   isca = 20
(   69) 
(   70)       iind = MAX( 1, kind )
(   71)       ijnd = MAX( 1, kjnd )
(   72) 
(   73)       il1 = kideb
(   74) 
(   75)       DO jn = 1, (kifin-kideb+1)/(iind*isca) + 1
(   76) 
(   77)         IF( il1 > kifin ) RETURN
(   78)         WRITE(kumout,'(/)')
(   79)         il2 = il1+iind*(isca-1)
(   80)         IF( il2 > kifin )   il2 = kifin
(   81) 
(   82)         IF( kscale == 1 ) THEN
(   83)             WRITE(kumout,'(4x,i14," to ",1i4," each ",1i4,/)') il1, il2, iind
(   84)             DO jj = kjfin, kjdeb, -ijnd
(   85)               WRITE (kumout,'(1x,i3,100i1)') jj, ( ktab(ji,jj), ji = il1, il2, iind )
(   86)             END DO  
(   87)         ELSEIF( kscale == 2 ) THEN
(   88)             WRITE(kumout,'(4x,i14," to ",1i4," each ",1i4,/)')il1, il2, iind
(   89)             DO jj = kjfin, kjdeb, -ijnd
(   90)               WRITE (kumout,'(1x,i3,60i2)') jj, ( ktab(ji,jj), ji = il1, il2, iind )
(   91)             END DO  
(   92)         ELSEIF( kscale == 3 ) THEN
(   93)             WRITE(kumout,'(4x,i14," to ",1i4," each ",1i4,/)')il1, il2, iind
(   94)             DO jj = kjfin, kjdeb, -ijnd
(   95)               WRITE (kumout,'(1x,i3,40i3)') jj, ( ktab(ji,jj), ji = il1, il2, iind )
(   96)             END DO  
(   97)         ELSEIF( kscale == 4 ) THEN
(   98)             WRITE(kumout,'(4x,30i4,/)') ( ji, ji = il1, il2, iind )
(   99)             DO jj = kjfin, kjdeb, -ijnd
(  100)               WRITE (kumout,'(1x,i3,30i4)') jj, ( ktab(ji,jj), ji = il1, il2, iind )
(  101)             END DO  
(  102)         ELSEIF( kscale == 5 ) THEN
(  103)             WRITE(kumout,'(4x,20i5,/)') ( ji, ji = il1, il2, iind )
(  104)             DO jj = kjfin, kjdeb, -ijnd
(  105)               WRITE (kumout,'(1x,i3,20i5)') jj, ( ktab(ji,jj), ji = il1, il2, iind )
(  106)             END DO  
(  107)         ELSE






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 3

(  108)             WRITE(kumout,'(4x,10i8,/)') ( ji, ji = il1, il2, iind )
(  109)             DO jj = kjfin, kjdeb, -ijnd
(  110)               WRITE (kumout,'(1x,i3,10i8)') jj, ( ktab(ji,jj), ji = il1, il2, iind )
(  111)             END DO  
(  112)         ENDIF
(  113) 
(  114)         il1 = il1 + iind * isca
(  115)       END DO  
(  116) 
(  117)    END SUBROUTINE prihin
(  118) 
(  119) 
(  120)    SUBROUTINE prihre( ptab, ki   , kj   , kideb, kifin ,   &
(  121)                       kind, kjdeb, kjfin, kjnd , pscale, kumout )
(  122)       !!----------------------------------------------------------------------
(  123)       !!                  ***  ROUTINE prihre  ***
(  124)       !!      
(  125)       !! ** purpose :   Print a real field with the format 10e12.4 or 20f6.2
(  126)       !!
(  127)       !! ** method  :   the print format is selected with the pscale argument
(  128)       !!
(  129)       !! History :
(  130)       !!   1.0  !  86-01  (P. Andrich)  Original code
(  131)       !!        !  89-11  (C. Levy)
(  132)       !!        !  92-02  (M. Imbard)
(  133)       !!        !  92-06  (M. Imbard)
(  134)       !!----------------------------------------------------------------------
(  135)       !! * Arguments
(  136)       INTEGER, INTENT( in ) ::   &
(  137)          ki, kj,                 &  ! array dimensions
(  138)          kideb, kifin, kind,     &  ! first and last index, increment for i
(  139)          kjdeb, kjfin, kjnd,     &  ! first and last index, increment for j
(  140)          kumout                     ! unit in which print
(  141)       REAL(wp), INTENT( in ) ::   &
(  142)          pscale                     ! = 0  print        ptab with e13.5 format
(  143)          !                          ! else print pscale*ptab with f6.2 format
(  144)       REAL(wp), DIMENSION(ki,kj), INTENT( in ) ::   &
(  145)          ptab                       ! integer 2D array to be print
(  146) 
(  147)       !! * Local variables
(  148)       INTEGER ::   ji, jj, jn       ! dummy loop indices
(  149)       INTEGER ::   isca, il1, il2   ! temporary integers
(  150) 
(  151)       isca = 10
(  152)       IF( pscale /= 0. )   isca=20
(  153) 
(  154)       il1 = kideb
(  155) 
(  156)       DO jn = 1, (kifin-kideb+1)/(kind*isca) + 1
(  157) 
(  158)         IF( il1 > kifin )   RETURN
(  159) 
(  160)         WRITE(kumout,9100)
(  161) 
(  162)         il2 = il1+kind*(isca-1)
(  163)         IF(il2 > kifin) il2 = kifin
(  164)         IF( pscale == 0.) THEN
(  165)             WRITE(kumout,9101) ( ji, ji = il1, il2, kind )






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 4

(  166)             DO jj = kjfin, kjdeb, -kjnd
(  167)               WRITE(kumout,9102) jj, ( ptab(ji,jj), ji = il1, il2, kind )
(  168)             END DO  
(  169)         ELSE
(  170)             WRITE(kumout,9103) ( ji, ji = il1, il2, kind )
(  171)             DO jj = kjfin, kjdeb, -kjnd
(  172)               WRITE(kumout,9104) jj, ( pscale*ptab(ji,jj), ji = il1, il2, kind )
(  173)             END DO  
(  174)         ENDIF
(  175)         il1 = il1+kind*isca
(  176) 
(  177)       END DO  
(  178) 
(  179)       ! formats
(  180)  9100 FORMAT(/)
(  181)  9101 FORMAT(10i12, /)
(  182)  9102 FORMAT(1x, i3, 10(1pe12.4))
(  183)  9103 FORMAT(3x, 20i6, /)
(  184)  9104 FORMAT(1x, i3, 1x, 20f6.2)
(  185) 
(  186)    END SUBROUTINE prihre
(  187) 
(  188) 
(  189)    SUBROUTINE prizre( ptab , ki   , kj   , kk   , kjcut ,   &
(  190)                       kideb, kifin, kid  , kkdeb, kkfin ,   &
(  191)                       kkd  , pscale, kumout )
(  192)       !!----------------------------------------------------------------------
(  193)       !!                      ***  ROUTINE prizre  ***
(  194)       !!
(  195)       !! ** purpose :   Print a vertical slab from a tridimentional real field
(  196)       !!
(  197)       !!   METHOD :
(  198)       !! ** method  :   the print format is selected with the argument pscale
(  199)       !!
(  200)       !! History :
(  201)       !!      original : 86-01 (o. Marti)
(  202)       !!      addition : 92-02 (M. Imbard)
(  203)       !!      addition : 92-06 (M. Imbard)
(  204)       !!----------------------------------------------------------------------
(  205)       !! * Arguments
(  206)       INTEGER, INTENT( in ) ::   &
(  207)          ki, kj, kk,             &  ! array dimensions
(  208)          kjcut,                  &  ! index j for the vertical slab
(  209)          kideb, kifin, kid,      &  ! first and last index, increment for i
(  210)          kkdeb, kkfin, kkd,      &  ! first and last index, increment for k
(  211)          kumout                     ! unit in which print
(  212)       REAL(wp), INTENT( in ) ::   &
(  213)          pscale                     ! = 0  print        ptab with e12.4 format
(  214)          !                          ! else print pscale*ptab with f6.2 format
(  215)       REAL(wp), DIMENSION(ki,kj,kk), INTENT( in ) ::   &
(  216)          ptab                       ! integer 3D array to be print
(  217) 
(  218)       !! * Local variables
(  219)       INTEGER ::   ji, jn, jk       ! dummy loop indices
(  220)       INTEGER ::   isca, il1, il2   ! temporary integers
(  221)       INTEGER ::   iind, iknd       !    "          "
(  222) 
(  223) 






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 5

(  224)       iind = kid
(  225)       iknd = kkd
(  226)       isca = 10
(  227)       IF( pscale /= 0.) isca = 20
(  228) 
(  229)       IF (iind == 0) iind = 1
(  230)       IF (iknd == 0) iknd = 1
(  231) 
(  232)       il1 = kideb
(  233) 
(  234)       DO jn = 1, (kifin-kideb+1)/(iind*isca) + 1
(  235) 
(  236)         IF(il1 > kifin) RETURN
(  237)         WRITE(kumout,9100)
(  238)         il2 = il1+iind*(isca-1)
(  239)         IF(il2 > kifin) il2 = kifin
(  240) 
(  241)         IF( pscale == 0.) THEN
(  242)             WRITE(kumout,9101) ( ji, ji = il1, il2, iind )
(  243)             DO jk = kkdeb, kkfin, iknd
(  244)               WRITE (kumout,9102) jk, ( ptab(ji,kjcut,jk), ji = il1, il2, iind )
(  245)             END DO  
(  246)         ELSE
(  247)             WRITE (kumout,9103) ( ji, ji = il1, il2, iind )
(  248)             DO jk = kkdeb, kkfin, iknd
(  249)               WRITE(kumout,9104)jk, ( pscale*ptab(ji,kjcut,jk), ji = il1, il2, iind )
(  250)             END DO  
(  251)         ENDIF
(  252) 
(  253)         il1 = il1+iind*isca
(  254)       END DO      
(  255) 
(  256)  9100 FORMAT(/)
(  257)  9101 FORMAT(10i12, /)
(  258)  9102 FORMAT(1x, i3, 10(1pe12.4))
(  259)  9103 FORMAT(3x, 20i6, /)
(  260)  9104 FORMAT(1x, i3, 1x, 20f6.1)
(  261) 
(  262)       END SUBROUTINE prizre
(  263) 
(  264)    !!======================================================================
(  265) END MODULE lib_print
