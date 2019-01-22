


PGF90 (Version     10.5)          09/23/2010  09:49:41      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir /usr/local/include
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /usr/local/include
          -idir /usr/include/mpich2-x86_64
          -idir /usr/include/mpich2-x86_64
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: lib_feti.F90

(    1) MODULE lib_feti
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE lib_feti   ***
(    4)    !! Ocean solver : FETI library
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_feti.F90,v 1.2 2005/03/27 18:34:47 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) #if defined key_feti
(   13)    !!  feti_creadr
(   14)    !!  feti_prext
(   15)    !!  feti_init
(   16)    !!  feti_iclr
(   17)    !!  feti_vclr
(   18)    !!  feti_vmov
(   19)    !!  feti_vneg    (used in lib_feti)
(   20)    !!  feti_vsub
(   21)    !!  feti_vadd    (used in lib_feti)
(   22)    !!  feti_ldlt    (used in lib_feti)
(   23)    !!  feti_desremo (used in lib_feti)
(   24)    !!  feti_mxv     (used in lib_feti)
(   25)    !!  feti_mxvadd  (used in lib_feti)
(   26)    !!  feti_saut    (used in lib_feti)
(   27)    !!  feti_moyen   (used in lib_feti)
(   28)    !!  feti_assemb  (used in lib_feti)
(   29)    !!  feti_poids
(   30)    !!  feti_probit  (used in lib_feti)
(   31)    !!  feti_inisub
(   32)    !!  feti_subound
(   33)    !!  feti_subdir
(   34)    !!  feti_listdir
(   35)    !!  feti_blodir
(   36)    !!  feti_numblo
(   37)    !!  feti_blomat
(   38)    !!  feti_blomat1
(   39)    !!  feti_nullsp
(   40)    !!  feti_project (used in lib_feti)
(   41)    !!  feti_etesolv (used in lib_feti)
(   42)    !!  feti_osaxpy  (used in lib_feti) 
(   43)    !!  feti_proe    (used in lib_feti)
(   44)    !!  feti_proet   (used in lib_feti)
(   45)    !!  feti_consrhs (used in lib_feti)
(   46)    !!  feti_dualschur






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 2

(   47)    !!  feti_extend  (used in lib_feti)
(   48)    !!  feti_restri  (used in lib_feti)
(   49)    !!  feti_halfint
(   50)    !!  feti_front
(   51)    !!  feti_calinv  (used in lib_feti)
(   52)    !!  feti_nbdia   (used in lib_feti)
(   53)    !!  feti_liblod  (used in lib_feti)
(   54)    !!  feti_liblos  (used in lib_feti)
(   55)    !!  feti_resloc
(   56)    !!  feti_pblos   (used in lib_feti)
(   57)    !!  feti_pbloi   (used in lib_feti) 
(   58)    !!  feti_proax
(   59)    !!
(   60)    !!  EXTERNAL used (BLAS): sgemv, strsv, saxpy, irecv
(   61)    !!
(   62)    !! References :
(   63)    !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
(   64)    !!      A domain decomposition solver to compute the barotropic 
(   65)    !!      component of an OGCM in the parallel processing field.
(   66)    !!      Ocean Modelling, issue 105, december 94.
(   67)    !!
(   68)    !! Modifications :
(   69)    !!      original :  98-02 (M. Guyon)
(   70)    !!------------------------------------------------------------------------
(   71)    !! * Modules used
(   72)    USE lib_mpp                 ! distribued memory computing
(   73) 
(   74) CONTAINS
(   75) 
(   76)       subroutine feti_creadr(ial,ialmax,lmax,ltab,iatab,chtab)
(   77) !  ial = derniere adresse libre dans le super-tableau de longueur lmax .
(   78) !  ltab = taille de la structure de donnee tab , d'adresse iatab .
(   79)       implicit none
(   80)       external mynode
(   81) 
(   82)       integer ial,ialmax,lmax,ltab,iatab
(   83)       integer iialmax,iibidon,mynode
(   84)       character(len=12) chtab
(   85) 
(   86)       iatab=ial
(   87)       ial=iatab+ltab
(   88)       ialmax=max0(ial,ialmax)
(   89)       if(ial > lmax) then 
(   90)           write(0,*)'  creation de la structure de donnee ',chtab
(   91)           write(0,*)'  depassement taille super-tableau :'
(   92)           write(0,*)'  taille maximale , taille demandee :',lmax,ial
(   93)           stop
(   94)       endif
(   95) 
(   96) !... test de la memoire
(   97) 
(   98)       iialmax=-ialmax
(   99)       call mpp_min(iialmax,1,iibidon)
(  100)       iialmax=-iialmax
(  101) !     if(mynode() == 0)then
(  102) !         write(0,*)'  creation de la structure de donnee ',chtab
(  103) !         write(0,*)'taille max vaut : ',iialmax
(  104) !     endif






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 3

(  105) 
(  106)       end subroutine feti_creadr
(  107) !***********************************************************************
(  108)       subroutine feti_prext(n,x)
(  109)       implicit none
(  110)       integer n
(  111)       REAL(kind=8)  x(n)
(  112)       REAL(kind=8)   xmin,xmax
(  113)       integer i
(  114) 
(  115)       xmin=x(1)
(  116)       xmax=x(1)
(  117)       do i=2,n
(  118)         xmin=min(xmin,x(i))
(  119)         xmax=max(xmax,x(i))
(  120)       end do
(  121) !     write(0,100)xmin,xmax
(  122)  100  format(8d10.3)
(  123) 
(  124)       end subroutine feti_prext
(  125) !*********************************************************************
(  126)       subroutine feti_init(noeuds,x)
(  127)       implicit none
(  128)       integer noeuds
(  129)       REAL(kind=8)  x(noeuds)
(  130)       integer ji
(  131) 
(  132)       do ji=1,noeuds
(  133)         x(ji)=1.
(  134)       end do
(  135) 
(  136)       end subroutine feti_init
(  137) !***********************************************************************
(  138)       subroutine feti_iclr(n,x)
(  139)       implicit none
(  140)       integer n
(  141)       integer x(n)
(  142)       integer i
(  143) 
(  144)       do i=1,n
(  145)         x(i)=0
(  146)       end do
(  147) 
(  148)       end subroutine feti_iclr
(  149) !***********************************************************************
(  150)       subroutine feti_vclr(n,x)
(  151)       implicit none
(  152)       integer n
(  153)       REAL(kind=8)  x(n)
(  154)       integer i
(  155) 
(  156)       do i=1,n
(  157)         x(i)=0.e0
(  158)       end do
(  159) 
(  160)       return
(  161)       end
(  162) !***********************************************************************






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 4

(  163)       subroutine feti_vmov(n,x,y)
(  164)       implicit none
(  165)       integer n
(  166)       REAL(kind=8)  x(n),y(n)
(  167)       integer i
(  168)       do i=1,n
(  169)         y(i)=x(i)
(  170)       end do
(  171)       return
(  172)       end
(  173) !***********************************************************************
(  174)       subroutine feti_vneg(n,x,y)
(  175)       implicit none
(  176)       integer n
(  177)       REAL(kind=8)  x(n),y(n)
(  178)       integer i
(  179) 
(  180)       do i=1,n
(  181)         y(i)=-x(i)
(  182)       end do
(  183) 
(  184)       return
(  185)       end
(  186) !***********************************************************************
(  187)       subroutine feti_vsub(n,x,y,z)
(  188)       implicit none
(  189)       integer n
(  190)       REAL(kind=8)  x(n),y(n),z(n)
(  191)       integer i
(  192)       do i=1,n
(  193)         z(i)=x(i)-y(i)
(  194)       end do
(  195)       return
(  196)       end
(  197) !***********************************************************************
(  198)       subroutine feti_vadd(n,x,y,z)
(  199)       implicit none
(  200)       integer n
(  201)       REAL(kind=8)  x(n),y(n),z(n)
(  202)       integer i
(  203) 
(  204)       do i=1,n
(  205)         z(i)=y(i)+x(i)
(  206)       end do
(  207) 
(  208)       return
(  209)       end
(  210) !***********************************************************************
(  211)       subroutine feti_ldlt(n,a,d,v,ndlblo,lisblo,nmdlblo,j0)
(  212) 
(  213) !...ndlblo nombre de liberte bloques = nombre de pivots nuls = dim(Ker(Ep))
(  214) 
(  215) !...nmdlblo = Max (ndlblo(Ep)) <= superstructure management in f77 context
(  216) !            p=1,Np
(  217) 
(  218)       implicit none
(  219)       integer n,ndlblo,nmdlblo
(  220)       integer j0






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 5

(  221)       REAL(kind=8)  a(n,n),v(n),d(n)
(  222)       integer lisblo(nmdlblo)
(  223)       integer i,j,k
(  224)       REAL(kind=8)  piv,epsilo
(  225) 
(  226) !  pivot de reference .
(  227)       piv=a(1,1)
(  228) !  epsilon .
(  229)       epsilo=1.e-10
(  230) !  factorisation de Crout d'une matrice pleine .
(  231)       do j = 1 , n
(  232) !  calcul du produit de la ligne numero j par la diagonale .
(  233)         do k = 1 , j - 1
(  234)           v(k) = a(j,k) * a(k,k)
(  235)         end do
(  236) !  calcul des termes de la partie inferieure de la colonne .
(  237)         call sgemv('N',n-j+1,j-1,-1.e0,a(j,1),n,v,1,1.e0,a(j,j),1)
(  238) !  detection des pivots nuls .
(  239)         if(a(j,j) < (piv*epsilo)) then
(  240)             ndlblo=ndlblo+1
(  241)             lisblo(ndlblo)=j0+j
(  242)             do i=1,j-1
(  243)               a(i,j)=0.e0
(  244)             end do
(  245)             do i=j+1,n
(  246)               a(j,i)=0.e0
(  247)             end do
(  248)             a(j,j)=0.e0
(  249)             d(j)=0.e0
(  250)         else
(  251)             piv=a(j,j)
(  252)             d(j)=1.e0/a(j,j)
(  253)             do i = j+1 , n
(  254)               a(i,j) = a(i,j) * d(j)
(  255)             end do
(  256)         endif
(  257)       end do
(  258) !  recopie de la partie triangulaire superieure .
(  259)       do i = 1, n
(  260)         do j = i+1, n
(  261)           a(i,j)=a(j,i)
(  262)         end do
(  263)       end do
(  264)       do i=1,n
(  265)         a(i,i)=d(i)
(  266)       end do
(  267) 
(  268)       return
(  269)       end
(  270) !**********************************************************************c
(  271)       subroutine feti_desremo(n,j0,a,y,x) 
(  272)       implicit none
(  273)       integer n,i,j,j0
(  274)       REAL(kind=8)  a(n,n),x(n),y(n)
(  275) 
(  276) !  descente remontee du systeme .
(  277) !  initialisation du vecteur resultat .
(  278)       do i=1,n






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 6

(  279)         x(i)=y(i)
(  280)       end do
(  281) !  descente du systeme L x = x , L est a diagonale unitaire .
(  282)       call strsv('L','N','U',n-j0,a(j0+1,j0+1),n,x(j0+1),1)
(  283) !  division par la diagonale .
(  284)       do j=j0+1,n
(  285)         x(j)=x(j)*a(j,j)
(  286)       end do          
(  287) !  remontee du systeme U x = x , U est a diagonale unitaire . 
(  288)       call strsv('U','N','U',n-j0,a(j0+1,j0+1),n,x(j0+1),1)
(  289) 
(  290)       return
(  291)       end
(  292) !***********************************************************************
(  293)       subroutine feti_mxv(n,m,a,x,y) 
(  294) !  calcul du produit matrice-vecteur :  y = A x
(  295)       implicit none
(  296)       integer n,m
(  297)       REAL(kind=8)  a(n,m),x(m),y(n)
(  298)       REAL(kind=8)  alpha,beta
(  299) !  initialisation du vecteur y .
(  300)       alpha=1.e0
(  301)       beta=0.e0
(  302)       call sgemv('N',n,m,alpha,a,n,x,1,beta,y,1)
(  303)       return
(  304)       end
(  305) !**********************************************************************c
(  306)       subroutine feti_mxvadd(n,m,a,x,y) 
(  307) !  calcul du produit matrice-vecteur :  y = y + A x
(  308)       implicit none
(  309)       integer n,m
(  310)       REAL(kind=8)  a(n,m),x(m),y(n)
(  311)       REAL(kind=8)  alpha,beta
(  312) !  initialisation du vecteur y .
(  313) 
(  314)       alpha=1.e0
(  315)       beta=1.e0
(  316)       if((n > 0).and.(m > 0))then
(  317)           call sgemv('N',n,m,alpha,a,n,x,1,beta,y,1)
(  318)       endif
(  319) 
(  320)       return
(  321)       end
(  322) !**********************************************************************c
(  323)       subroutine feti_saut(numnes,listnes,    &
(  324)           plistin,numins,listin,              &
(  325)                 numnos,v,w,bufin,bufout)
(  326) !!        narea,numnos,v,w,bufin,bufout)
(  327) 
(  328) !  computation of the jump of a field on the interfaces, in the case 
(  329) !  where subdomain number i is allocated to processor number i-1 .
(  330) 
(  331) !     INPUTS : numnes  =  number of neighbouring subdomains
(  332) !               listnes =  list of neighbouring subdomains
(  333) !               lint    =  list of inner interfaces
(  334) !               lext    =  list of extern interfaces
(  335) !               plistin =  pointer of sublists of interface nodes
(  336) !               numins  =  number of interface nodes






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 7

(  337) !               listin  =  list of interface nodes
(  338) !               narea     =  subdomain number
(  339) !               numnos  =  number of nodes in the subdomain
(  340) !               v       =  local field
(  341) 
(  342) !     OUTPUTS : w  =  jump of v
(  343) 
(  344) !     WORKSPACE : bufin   =  contributions of local field on interface
(  345) !                 bufout  =  contributions of outer fields on interface
(  346) !                 work
(  347)       implicit none 
(  348) !      external irecv
(  349) !!    integer irecv, work(100)
(  350) !!    integer numnes,numins,numnos,narea
(  351)       integer numnes,numins,numnos
(  352)       integer listnes(numnes,3), plistin(numnes+1), listin(numins)
(  353)       REAL(kind=8)  v(numnos),w(numins)
(  354)       REAL(kind=8)  bufin(numins),bufout(numins)
(  355)       integer i,j,l
(  356)       integer mesg
(  357) 
(  358)       mesg=128
(  359) !  receiving the values of the outer fields on the interface
(  360) !     do i=1,numnes
(  361) !       l=(plistin(i+1)-plistin(i))*8
(  362) !       work(i)=irecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
(  363) !     end do
(  364) !  gathering the values of the local field on the interface
(  365)       do i=1,numnes
(  366)         do j=plistin(i)+1,plistin(i+1)
(  367)           bufin(j)=v(listin(j))
(  368)         end do
(  369)       end do
(  370) !  sending the values of the local field on the interface
(  371)       do i=1,numnes
(  372)         l=(plistin(i+1)-plistin(i))*8
(  373) !       call mppsend(mesg+narea-1,bufin(plistin(i)+1),l,listnes(i)-1,0)
(  374)         call mppsend(mesg+listnes(i,3)-1,bufin(plistin(i)+1),l,   &
(  375)             listnes(i,1)-1,0)
(  376)       end do
(  377) !  computing the jump on each interface 
(  378)       do i=1,numnes
(  379) !       call msgwait(work(i))
(  380)         l=(plistin(i+1)-plistin(i))*8
(  381) !       call mpprecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
(  382)         call mpprecv(mesg+listnes(i,2)-1,bufout(plistin(i)+1),l)
(  383)       end do
(  384)       do j=1,numins
(  385)         w(j)=bufin(j)-bufout(j)
(  386)       end do
(  387) 
(  388)       return
(  389)       end
(  390) !**********************************************************************c
(  391)       subroutine feti_moyen(numnes,listnes,   &
(  392)           plistin,numins,listin,              &
(  393)           narea,numnos,weight,v,w,bufin,bufout)
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 8

(  395) !  averaging a field on the interfaces, in the case where subdomain
(  396) !  number i is allocated to processor number i-1 .
(  397) 
(  398) !     INPUTS : numnes  =  number of neighbouring subdomains
(  399) !               listnes =  list of neighbouring subdomains
(  400) !               lint    =  list of inner interfaces
(  401) !               lext    =  list of extern interfaces
(  402) !               plistin =  pointer of sublists of interface nodes
(  403) !               numins  =  number of interface nodes
(  404) !               listin  =  list of interface nodes
(  405) !               narea    =  subdomain number
(  406) !               numnos  =  number of nodes in the subdomain
(  407) !               weight  =  weighting vector
(  408) !               v       =  local field
(  409) 
(  410) !     OUTPUTS :      w  =  averaged field
(  411) 
(  412) !     WORKSPACE : bufin   =  contributions of local field on interface
(  413) !                 bufout  =  contributions of outer fields on interface
(  414) !                 work
(  415)       implicit none 
(  416) !     external irecv
(  417) !!    integer irecv, work(100)
(  418)       integer numnes,numins,numnos,narea
(  419)       integer listnes(numnes,3), plistin(numnes+1), listin(numins)
(  420)       REAL(kind=8)  weight(numnos),v(numnos),w(numnos)
(  421)       REAL(kind=8)  bufin(numins),bufout(numins)
(  422)       integer i,j,l,n
(  423)       integer mesg
(  424) 
(  425)       mesg=128*2
(  426) !  initialization
(  427)       call feti_vmov(numnos,v,w)
(  428) !  receiving the values of the outer fields on the interface
(  429) !     do i=1,numnes
(  430) !       l=(plistin(i+1)-plistin(i))*8
(  431) !       work(i)=irecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
(  432) !     end do
(  433) !  gathering the values of the local field on the interface
(  434)       do i=1,numnes
(  435)         do j=plistin(i)+1,plistin(i+1)
(  436)           bufin(j)=v(listin(j))
(  437)         end do
(  438)       end do
(  439) !  sending the values of the local field on the interface
(  440)       do i=1,numnes
(  441)         l=(plistin(i+1)-plistin(i))*8
(  442) !       call mppsend(mesg+narea-1,bufin(plistin(i)+1),l,listnes(i)-1,0)
(  443)         call mppsend(mesg+listnes(i,3)-1,bufin(plistin(i)+1),l,   &
(  444)             listnes(i,1)-1,0)
(  445)       end do
(  446) !  computing the jump on each interface 
(  447)       do i=1,numnes
(  448) !       call msgwait(work(i))
(  449)         l=(plistin(i+1)-plistin(i))*8
(  450) !       call mpprecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
(  451)         call mpprecv(mesg+listnes(i,2)-1,bufout(plistin(i)+1),l)
(  452)       end do






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 9

(  453)       do j=1,numins
(  454)         n=listin(j)
(  455)         w(n)=w(n)+bufout(j)
(  456)       end do
(  457)       do n=1,numnos
(  458)         w(n)=w(n)*weight(n)
(  459)       end do
(  460) 
(  461)       return
(  462)       end
(  463) !**********************************************************************c
(  464)       subroutine feti_assemb(numnes,listnes,   &
(  465)           plistin,numins,listin,               &
(  466)           narea,numnos,v,w,bufin,bufout)
(  467) 
(  468) !  assembling a field on the interfaces, in the case where subdomain
(  469) !  number i is allocated to processor number i-1 .
(  470) 
(  471) !     INPUTS :      numnes  =  number of neighbouring subdomains
(  472) !               listnes =  list of neighbouring subdomains
(  473) !               lint    =  list of inner interfaces
(  474) !               lext    =  list of extern interfaces
(  475) !               plistin =  pointer of sublists of interface nodes
(  476) !               numins  =  number of interface nodes
(  477) !               listin  =  list of interface nodes
(  478) !               narea     =  subdomain number
(  479) !               numnos  =  number of nodes in the subdomain
(  480) !               v       =  unassembled field
(  481) 
(  482) !     OUTPUTS :      w  =  assembled field
(  483) 
(  484) !     WORKSPACE : bufin   =  contributions of local field on interface
(  485) !                 bufout  =  contributions of outer fields on interface
(  486) !                 work
(  487)       implicit none 
(  488) !      external irecv
(  489) !!    integer irecv, work(100)
(  490)       integer numnes,numins,numnos,narea
(  491)       integer listnes(numnes,3), plistin(numnes+1), listin(numins)
(  492)       REAL(kind=8)  v(numnos),w(numnos)
(  493)       REAL(kind=8)  bufin(numins),bufout(numins)
(  494)       integer i,j,l,n
(  495)       integer mesg
(  496) 
(  497)       mesg=128*3
(  498) !  initialization
(  499)       call feti_vmov(numnos,v,w)
(  500) !  receiving the values of the outer fields on the interface
(  501) !     do i=1,numnes
(  502) !       l=(plistin(i+1)-plistin(i))*8
(  503) !       work(i)=irecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
(  504) !     end do
(  505) !  gathering the values of the local field on the interface
(  506)       do i=1,numnes
(  507)         do j=plistin(i)+1,plistin(i+1)
(  508)           bufin(j)=v(listin(j))
(  509)         end do
(  510)       end do






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 10

(  511) !  sending the values of the local field on the interface
(  512)       do i=1,numnes
(  513)         l=(plistin(i+1)-plistin(i))*8
(  514) !       write(90,*)'envoie a ',listnes(i),'de ',l,'octets'
(  515) !       write(90,*)'tag est : ',listnes(i,3)
(  516) !       call mppsend(mesg+narea-1,bufin(plistin(i)+1),l,listnes(i)-1,0)
(  517)         call mppsend(mesg+listnes(i,3)-1,bufin(plistin(i)+1),l,   &
(  518)             listnes(i,1)-1,0)
(  519)       end do
(  520) !  computing the jump on each interface 
(  521)       do i=1,numnes
(  522) !       call msgwait(work(i))
(  523)         l=(plistin(i+1)-plistin(i))*8
(  524) !       write(90,*)'recoie de ',listnes(i),'de ',l,'octets'
(  525) !       write(90,*)'recois avec tag :',listnes(i,2)
(  526) !       call mpprecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
(  527)         call mpprecv(mesg+listnes(i,2)-1,bufout(plistin(i)+1),l)
(  528)       end do
(  529)       do j=1,numins
(  530)         n=listin(j)
(  531)         w(n)=w(n)+bufout(j)
(  532)       end do
(  533) 
(  534)       return
(  535)       end
(  536) !**********************************************************************c
(  537) !!    subroutine feti_poids(numnes,listnes,plistin,numins,listin,   &
(  538) !!        narea,numnos,w,bufin,bufout)
(  539)       subroutine feti_poids(numnes,                numins,listin,   &
(  540)                 numnos,w             )
(  541) 
(  542) !  computation of the subdomain weighting vector, in the case where 
(  543) !  subdomain number i is allocated to processor number i-1 .
(  544) 
(  545) !     INPUTS :      numnes  =  number of neighbouring subdomains
(  546) !               listnes =  list of neighbouring subdomains
(  547) !               plistin =  pointer of sublists of interface nodes
(  548) !               numins  =  number of interface nodes
(  549) !               listin  =  list of interface nodes
(  550) !               narea     =  subdomain number
(  551) !               numnos  =  number of nodes in the subdomain
(  552) 
(  553) !     OUTPUTS :      w  =  weighting vector
(  554) 
(  555) !     WORKSPACE : bufin   =  contributions of local field on interface
(  556) !                 bufout  =  contributions of outer fields on interface
(  557) !                 work
(  558)       implicit none 
(  559) !!    integer numnes,numins,numnos,narea
(  560)       integer numnes,numins,numnos
(  561) !!    integer listnes(numnes), plistin(numnes+1), listin(numins)
(  562)       integer                                     listin(numins)
(  563)       REAL(kind=8)  w(numnos)
(  564) !!    REAL(kind=8)  bufin(numins),bufout(numins)
(  565)       integer j,n
(  566) 
(  567) !  initialization
(  568)       do n=1,numnos






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 11

(  569)         w(n)=1.e0
(  570)       end do
(  571) !  computing the weight on each interface 
(  572)       do j=1,numins
(  573)         n=listin(j)
(  574)         w(n)=w(n)+1.e0
(  575)       end do
(  576)       do n=1,numnos
(  577)         w(n)=1.e0/w(n)
(  578)       end do
(  579) 
(  580)       return
(  581)       end
(  582) !**********************************************************************c
(  583)       subroutine feti_probit(numnes,listnes,plistin,numins,listin,   &
(  584)           numnos,w,bitw)
(  585) 
(  586) !  computing the local right-hand-side associated with a field
(  587) !  on the interface .
(  588) 
(  589) !     INPUTS :      numnes  =  number of neighbouring subdomains
(  590) !               listnes =  list of neighbouring subdomains
(  591) !               plistin =  pointer of sublists of interface nodes
(  592) !               numins  =  number of interface nodes
(  593) !               listin  =  list of interface nodes
(  594) !               numnos  =  number of nodes in the subdomain
(  595) !               w       =  interface field
(  596) 
(  597) !     OUTPUTS :      bitw  =  local right-hand side
(  598) 
(  599)       implicit none 
(  600)       integer numnes,numins,numnos
(  601)       integer listnes(numnes,3), plistin(numnes+1), listin(numins)
(  602)       REAL(kind=8)  w(numins)
(  603)       REAL(kind=8)  bitw(numnos)
(  604)       integer i,j,n
(  605) 
(  606) !  initialization
(  607)       do n=1,numnos
(  608)         bitw(n)=0.e0
(  609)       end do
(  610) !  assembling the values on the interfaces 
(  611)       do i=1,numnes
(  612)         do j=plistin(i)+1,plistin(i+1)
(  613)           n=listin(j)
(  614)           bitw(n)=bitw(n)+w(j)
(  615)         end do
(  616)       end do
(  617) 
(  618)       return
(  619)       end
(  620) !**********************************************************************c
(  621)       subroutine feti_inisub(ni,nj,ibondi,ibondj,iperio,   &
(  622)           ibsw,ibnw,ibse,ibne,                             &
(  623)           ninterf,ninterfc,nni,nnic)
(  624)       implicit none
(  625)       integer ni,nj,ibondi,ibondj,iperio
(  626)       integer ibsw,ibnw,ibse,ibne






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 12

(  627)       integer ibsw2,ibnw2,ibse2,ibne2
(  628)       integer ninterf,ninterfc,nni,nnic
(  629) !     integer i,j
(  630) 
(  631) !  determination de la position du sous-domaine dans la grille
(  632) !  cela depnd de ibondi et ibondj
(  633) 
(  634) !  nombre de sous-domaines voisins et de la longueur du tableau des
(  635) !  noeuds voisins
(  636) 
(  637) !  initialisation
(  638) 
(  639)       ninterf=0
(  640)       nni=0
(  641) 
(  642) !...first, the "segment-points"
(  643) 
(  644)       if((ibondi == 2).AND.(iperio /= 1))then
(  645) 
(  646) !  la periodicite : pas d'interface east / west si nbondi = 2
(  647) !  ET nperio != 1
(  648) 
(  649) !      elseif(iperio == 1) then
(  650) 
(  651) !...ca doit marcher!?
(  652) 
(  653)       else
(  654) 
(  655) !  west
(  656) 
(  657)           if(ibondi /= -1)then
(  658)               ninterf=ninterf+1
(  659)               nni=nni+nj+1
(  660)           endif
(  661) 
(  662) !  east
(  663) 
(  664)           if(ibondi /= 1)then
(  665)               ninterf=ninterf+1
(  666)               nni=nni+nj+1
(  667)           endif
(  668) 
(  669)       endif
(  670) 
(  671) !  south
(  672) 
(  673)       if(ibondj /= -1.AND.ibondj /= 2)then
(  674)           ninterf=ninterf+1
(  675)           nni=nni+ni+1
(  676)       endif
(  677) 
(  678) !  north
(  679) 
(  680)       if(ibondj /= 1.AND.ibondj /= 2)then
(  681)           ninterf=ninterf+1
(  682)           nni=nni+ni+1
(  683)       endif
(  684) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 13

(  685)       nnic=nni
(  686)       ninterfc=ninterf
(  687) 
(  688) 
(  689) !...second, the "corner-points"
(  690) 
(  691) !...determination of the neighboorings, they depends on iperio
(  692) !   boundary condition => local description : ipXX becomes ipXX2...
(  693) 
(  694)       ibnw2 = ibnw
(  695)       ibne2 = ibne
(  696)       ibsw2 = ibsw
(  697)       ibse2 = ibse
(  698) 
(  699) !...iperio boundary condition effect
(  700) 
(  701)       if(iperio == 1) then
(  702) 
(  703) !...general case : Earth == infinite tube
(  704) 
(  705)           ibnw2 = 1
(  706)           ibne2 = 1
(  707)           ibsw2 = 1
(  708)           ibse2 = 1
(  709) 
(  710) !...real boundary condition
(  711) 
(  712)           if(ibondj == -1.OR.ibondj == 2) then
(  713)               ibsw2 = 0
(  714)               ibse2 = 0
(  715)           endif
(  716) 
(  717)           if(ibondj == 1.OR.ibondj == 2) then
(  718)               ibnw2 = 0
(  719)               ibne2 = 0
(  720)           endif
(  721) 
(  722)       endif
(  723) 
(  724) !  la periodicite : pas de coin si nbondi = 2
(  725) !  ET nperio != 1
(  726) 
(  727) !  north-west
(  728) 
(  729)       if(ibnw == 1)then
(  730)           ninterfc=ninterfc+1
(  731)           nnic=nnic+1
(  732)       endif
(  733) 
(  734) !  north-east
(  735) 
(  736)       if(ibne == 1)then
(  737)           ninterfc=ninterfc+1
(  738)           nnic=nnic+1
(  739)       endif
(  740) 
(  741) !  south-west
(  742) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 14

(  743)       if(ibsw == 1)then
(  744)           ninterfc=ninterfc+1
(  745)           nnic=nnic+1
(  746)       endif
(  747) 
(  748) !  south-east
(  749) 
(  750)       if(ibse == 1)then
(  751)           ninterfc=ninterfc+1
(  752)           nnic=nnic+1
(  753)       endif
(  754) 
(  755) !       write(0,*)'  nbre de voisins =',ninterf
(  756) !       write(0,*)'  longueur des interfaces =',nni
(  757) !       write(0,*)'  longueur des interfaces =',nnic
(  758) 
(  759)       return
(  760)       end
(  761) !**********************************************************************c
(  762)       subroutine feti_subound(ni,nj,ildi,ilei,ildj,ilej,   &
(  763)           imoi,ibondi,ibondj,iperio,   &
(  764)           ninterf,ninterfc,   &
(  765)           iowe,ioea,ioso,iono,   &
(  766)           ibsw,ibnw,ibse,ibne,   &
(  767)           ipsw,ipnw,ipse,ipne,   &
(  768)           nsdvois,nsdvoisc,   &
(  769)           plistin,nni,listin)
(  770) 
(  771)       implicit none
(  772)       integer ni,nj,ildi,ilei,ildj,ilej
(  773)       integer imoi,ibondi,ibondj,iperio
(  774)       integer iowe,ioea,ioso,iono
(  775)       integer ibsw,ibnw,ibse,ibne
(  776)       integer ipsw,ipnw,ipse,ipne
(  777)       integer ninterf,ninterfc,nni
(  778)       integer nsdvois(ninterf,3),nsdvoisc(ninterfc,3)
(  779)       integer plistin(ninterfc+1),listin(nni)
(  780)       integer     nint,nni0,ii,jj
(  781) !     integer i,j,nint,nni0,ii,jj
(  782) 
(  783)       external mynode
(  784)       integer  mynode
(  785) 
(  786) !  determination de la position du sous-domaine dans la grille
(  787) !  cela depnd de ibondi et ibondj
(  788) 
(  789) ! liste des sous-domaines voisins et des noeuds interface
(  790) 
(  791) ! be carefull!!! in the FETI algorithm, area number is considered
(  792) ! instead of process number
(  793) 
(  794) !  initialisation
(  795) 
(  796)       nint=0
(  797)       nni0=0
(  798)       plistin(1)=0
(  799) 
(  800) !...first, the "segment-points"






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 15

(  801) 
(  802)       if((ibondi == 2).AND.(iperio /= 1))then
(  803) 
(  804) !  la periodicite : pas d'interface east / west si nbondi = 2
(  805) !  ET nperio != 1
(  806) 
(  807)       elseif(iperio == 1) then
(  808) 
(  809) !  west
(  810) 
(  811)           nint=nint+1
(  812)           nsdvoisc(nint,1)=imoi
(  813)           do jj=1,nj+1
(  814)             listin(nni0+jj)=(jj-1)*(ni+1)+1
(  815)           end do
(  816)           nni0=nni0+nj+1
(  817)           plistin(nint+1)=nni0
(  818)           nsdvoisc(nint,2)=1
(  819)           nsdvoisc(nint,3)=2
(  820) 
(  821) !  east
(  822) 
(  823)           nint=nint+1
(  824)           nsdvoisc(nint,1)=imoi
(  825)           do jj=1,nj+1
(  826)             listin(nni0+jj)=(jj-1)*(ni+1)+ilei
(  827)           end do
(  828)           nni0=nni0+nj+1
(  829)           plistin(nint+1)=nni0
(  830)           nsdvoisc(nint,2)=2
(  831)           nsdvoisc(nint,3)=1
(  832) 
(  833)       else
(  834) 
(  835) !  west
(  836) 
(  837)           if(ibondi /= -1)then
(  838)               nint=nint+1
(  839)               nsdvoisc(nint,1)=(iowe+1)
(  840)               do jj=1,nj+1
(  841)                 listin(nni0+jj)=(jj-1)*(ni+1)+1
(  842)               end do
(  843)               nni0=nni0+nj+1
(  844)               plistin(nint+1)=nni0
(  845)               nsdvoisc(nint,2)=1
(  846)               nsdvoisc(nint,3)=2
(  847)           endif
(  848) 
(  849) !  east
(  850) 
(  851)           if(ibondi /= 1)then
(  852)               nint=nint+1
(  853)               nsdvoisc(nint,1)=(ioea+1)
(  854)               do jj=1,nj+1
(  855)                 listin(nni0+jj)=(jj-1)*(ni+1)+ilei
(  856)               end do
(  857)               nni0=nni0+nj+1
(  858)               plistin(nint+1)=nni0






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 16

(  859)               nsdvoisc(nint,2)=2
(  860)               nsdvoisc(nint,3)=1
(  861)           endif
(  862) 
(  863)       endif
(  864) 
(  865) !  south
(  866) 
(  867)       if(ibondj /= -1.AND.ibondj /= 2)then
(  868)           nint=nint+1
(  869)           nsdvoisc(nint,1)=(ioso+1)
(  870)           do ii=1,ni+1
(  871)             listin(nni0+ii)=ii
(  872)           end do
(  873)           nni0=nni0+ni+1
(  874)           plistin(nint+1)=nni0
(  875)           nsdvoisc(nint,2)=3
(  876)           nsdvoisc(nint,3)=4
(  877)       endif
(  878) 
(  879) !  north
(  880) 
(  881)       if(ibondj /= 1.AND.ibondj /= 2)then
(  882)           nint=nint+1
(  883)           nsdvoisc(nint,1)=(iono+1)
(  884)           do ii=1,ni+1
(  885)             listin(nni0+ii)=(ilej-1)*(ni+1)+ii
(  886)           end do
(  887)           nni0=nni0+ni+1
(  888)           plistin(nint+1)=nni0
(  889)           nsdvoisc(nint,2)=4
(  890)           nsdvoisc(nint,3)=3
(  891)       endif
(  892) 
(  893) !...second, the "corner-points"
(  894) 
(  895) !...determination of the neighboorings, they depends on iperio 
(  896) !   boundary condition => local description : inimpp or inimpp2...
(  897) 
(  898) !  la periodicite : pas de coin si nbondi = 2
(  899) !  ET nperio != 1
(  900) 
(  901) !!!!!!!!!verifier la definition des voisins-coin
(  902) 
(  903) !  north-west
(  904) 
(  905)       if(ibnw == 1)then
(  906)           nint=nint+1
(  907)           nsdvoisc(nint,1)=(ipnw+1)
(  908)           listin(nni0+1)=(ilej-1)*(ni+1)+1
(  909)           nni0=nni0+1
(  910)           plistin(nint+1)=nni0
(  911)           nsdvoisc(nint,2)=5
(  912)           nsdvoisc(nint,3)=8
(  913)       endif
(  914) 
(  915) !  north-east
(  916) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 17

(  917)       if(ibne == 1)then
(  918)           nint=nint+1
(  919)           nsdvoisc(nint,1)=(ipne+1)
(  920)           listin(nni0+1)=(ilej-1)*(ni+1)+ilei
(  921)           nni0=nni0+1
(  922)           plistin(nint+1)=nni0
(  923)           nsdvoisc(nint,2)=6
(  924)           nsdvoisc(nint,3)=7
(  925)       endif
(  926) 
(  927) !  south-west
(  928) 
(  929)       if(ibsw == 1)then
(  930)           nint=nint+1
(  931)           nsdvoisc(nint,1)=(ipsw+1)
(  932)           listin(nni0+1)=1
(  933)           nni0=nni0+1
(  934)           plistin(nint+1)=nni0
(  935)           nsdvoisc(nint,2)=7
(  936)           nsdvoisc(nint,3)=6
(  937)       endif
(  938) 
(  939) !  south-east
(  940) 
(  941)       if(ibse == 1)then
(  942)           nint=nint+1
(  943)           nsdvoisc(nint,1)=(ipse+1)
(  944)           listin(nni0+1)=ilei
(  945)           nni0=nni0+1
(  946)           plistin(nint+1)=nni0
(  947)           nsdvoisc(nint,2)=8
(  948)           nsdvoisc(nint,3)=5
(  949)       endif
(  950) 
(  951) !...je compte les coins
(  952) 
(  953)       do ii=1,ninterf
(  954)         nsdvois(ii,1) = nsdvoisc(ii,1)
(  955)         nsdvois(ii,2) = nsdvoisc(ii,2)
(  956)         nsdvois(ii,3) = nsdvoisc(ii,3)
(  957)       enddo
(  958) 
(  959) 
(  960) !  impressions
(  961) 
(  962) !     write(*,*)'  sous-domaine :',imoi
(  963) !     write(*,*)'  nombre de sous-domaines voisins :',ninterfc
(  964) !      do ii=1,ninterfc
(  965) !      write(*,*)'  sous-domaine voisin :',nsdvoisc(ii,1)
(  966) !     write(*,*)'  nombre de noeuds interface :',   &
(  967) !         plistin(ii+1)-plistin(ii)
(  968) !     write(*,*)'  noeuds interface :',(listin(jj),jj=plistin(ii)+1,   &
(  969) !         plistin(ii+1))
(  970) !     end do
(  971) 
(  972)       return
(  973)       end
(  974) !**********************************************************************c






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 18

(  975)       subroutine feti_subdir(ni,nj,noeuds,ndir,mgcnum)
(  976)       implicit none 
(  977)       integer ni,nj,noeuds,ndir
(  978)       integer mgcnum(ni+1,nj+1)
(  979)       integer i,j
(  980) !     integer i,j,k
(  981) 
(  982) !  determination de la position du sous-domaine dans la grille
(  983) !  narea = (j-1)*iglo + i , avec 1<=i<=iglo et 1<=j<=jglo 
(  984) 
(  985) 
(  986) 
(  987) ! 1. NUMBER THE GRID-POINTS USING BMASK
(  988) ! -------------------------------------
(  989) 
(  990) 
(  991)       ndir=0
(  992)       do j=1,nj+1
(  993)         do i=1,ni+1
(  994)           ndir=ndir+1-mgcnum(i,j)
(  995) !      do k=1,noeuds
(  996) !        ndir=ndir+1-mgcnum(k,1)
(  997)         enddo
(  998)       enddo
(  999) 
( 1000)       return
( 1001)       end
( 1002) !**********************************************************************c
( 1003)       subroutine feti_listdir(ni,nj,logdir,ndir,lisdir)
( 1004)       implicit none 
( 1005)       integer ni,nj,ndir
( 1006)       integer logdir(ni,nj),lisdir(ndir)
( 1007)       integer i0,ji,jj
( 1008) 
( 1009) !  creation de la liste des degres de liberte bloques
( 1010)       i0=0
( 1011)       do ji=1,ni
( 1012)         do jj=1,nj
( 1013)           if(logdir(ji,jj) == 0) then
( 1014)               i0=i0+1
( 1015)               lisdir(i0)=ji+(jj-1)*ni
( 1016)           endif
( 1017)         end do
( 1018)       end do
( 1019)       if(i0 /= ndir) then
( 1020) !      write(*,*)'  nombre de ddl bloques/prevus :',i0,ndir
( 1021) !      write(*,*)'liste des ddl Dirichlet :',(lisdir(ji),ji=1,ndir)
( 1022)           stop
( 1023)       endif
( 1024) 
( 1025)       return
( 1026)       end
( 1027) !***********************************************************************
( 1028)       subroutine feti_blodir(n,x,ndlblo,list)
( 1029)       implicit none
( 1030)       integer n,ndlblo
( 1031)       integer list(ndlblo)
( 1032)       REAL(kind=8)  x(n)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 19

( 1033)       integer i
( 1034) 
( 1035) !  remise a zero des ddl bloques
( 1036)       do i=1,ndlblo
( 1037)         x(list(i))=0.e0
( 1038)       end do
( 1039) 
( 1040)       return
( 1041)       end
( 1042) !***********************************************************************
( 1043)       subroutine feti_numblo(ndlblo,lisblo)
( 1044)       implicit none
( 1045)       external mynode
( 1046) !!    integer mynode,indlblog,iibidon
( 1047)       integer mynode
( 1048)       integer ndlblo
( 1049)       integer lisblo(ndlblo)
( 1050) !!    integer n
( 1051) 
( 1052) !      indlblog=-ndlblo
( 1053) !      call mpp_min(indlblog,1,iibidon)
( 1054) !      indlblog=-indlblog
( 1055) !      if(mynode() == 0)THEN
( 1056) !      write(0,*)'  nombre de degres de liberte flottants :',ndlblo
( 1057) !      write(0,*)'  liste :',(lisblo(n),n=1,ndlblo)
( 1058) !      endif
( 1059) 
( 1060)       return
( 1061)       end
( 1062) !***********************************************************************
( 1063)       subroutine feti_blomat(ni,nj,a,ndlblo,lisblo)
( 1064)       implicit none
( 1065)       integer ni,nj
( 1066)       REAL(kind=8)  a(ni,nj,5)
( 1067)       integer ndlblo
( 1068)       integer lisblo(ndlblo)
( 1069) !!    integer i,j,k,l,n
( 1070)       integer i,j,k,n
( 1071) 
( 1072)       do n=1,ndlblo
( 1073) !  degre de liberte bloque : lisblo(n) = (j-1)*ni + i
( 1074)         i=mod(lisblo(n)-1,ni)+1
( 1075)         j=((lisblo(n)-1)/ni)+1
( 1076) !  annulation des coefficients sur la ligne
( 1077)         do k=1,5
( 1078)           a(i,j,k)=0.e0
( 1079)         end do
( 1080)         a(i,j,3)=1.e0
( 1081) !  annulation des coefficients dans les colonnes correspondantes
( 1082)         if(i > 1) then
( 1083)             a(i-1,j,4)=0.e0
( 1084)         endif
( 1085)         if(i < ni) then
( 1086)             a(i+1,j,2)=0.e0
( 1087)         endif
( 1088)         if(j > 1) then
( 1089)             a(i,j-1,5)=0.e0
( 1090)         endif






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 20

( 1091)         if(j < nj) then
( 1092)             a(i,j+1,1)=0.e0
( 1093)         endif
( 1094)       end do
( 1095) 
( 1096)       return
( 1097)       end
( 1098) !***********************************************************************
( 1099)       subroutine feti_blomat1(ni,nj,a,ndlblo,lisblo,nsp)
( 1100)       implicit none
( 1101)       integer ni,nj,ndlblo
( 1102)       REAL(kind=8)  a(ni,nj,5),nsp(ni,nj,ndlblo)
( 1103)       integer lisblo(ndlblo)
( 1104) !!    integer i,j,k,l,n
( 1105)       integer i,j,n
( 1106) 
( 1107)       do n=1,ndlblo
( 1108) !  degre de liberte bloque : lisblo(n) = (j-1)*ni + i
( 1109)         i=mod(lisblo(n)-1,ni)+1
( 1110)         j=((lisblo(n)-1)/ni)+1
( 1111) !  annulation des coefficients dans les colonnes correspondantes
( 1112)         if(i > 1) then
( 1113)             nsp(i-1,j,n)=-a(i-1,j,4)
( 1114)         endif
( 1115)         if(i < ni) then
( 1116)             nsp(i+1,j,n)=-a(i+1,j,2)
( 1117)         endif
( 1118)         if(j > 1) then
( 1119)             nsp(i,j-1,n)=-a(i,j-1,5)
( 1120)         endif
( 1121)         if(j < nj) then
( 1122)             nsp(i,j+1,n)=-a(i,j+1,1)
( 1123)         endif
( 1124)       end do
( 1125) 
( 1126)       return
( 1127)       end
( 1128) !***********************************************************************
( 1129)       subroutine feti_nullsp(noeuds,ni,nj,lpblo,blo,a,ndlblo,lisblo,nsp,z)
( 1130) !  calcul du noyau de la matrice
( 1131)       implicit none
( 1132)       integer noeuds,ni,nj,lpblo,ndlblo
( 1133)       REAL(kind=8)  a(ni,nj,5)
( 1134)       integer lisblo(ndlblo)
( 1135)       REAL(kind=8)  blo(lpblo),nsp(noeuds,ndlblo),z(noeuds)
( 1136)       external sdot
( 1137)       REAL(kind=8)  sdot,res
( 1138)       integer j,i
( 1139) 
( 1140) !  calcul de [Aii]-1 * Aib
( 1141)       do j=1,ndlblo
( 1142)         call feti_resloc(noeuds,ni,nj,a,lpblo,blo,nsp(1,j),nsp(1,j),z)
( 1143)       end do
( 1144) !  remise a identite du bloc diagonal correspondant aux modes bloques
( 1145)       do i=1,ndlblo
( 1146)         do j=1,ndlblo
( 1147)           nsp(lisblo(i),j)=0.e0
( 1148)         end do






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 21

( 1149)         nsp(lisblo(i),i)=1.e0
( 1150)       end do
( 1151) !  verifications
( 1152)       do j=1,ndlblo
( 1153)         call feti_proax(noeuds,ni,nj,a,nsp(1,j),z)
( 1154)         res=sqrt(sdot(noeuds,z,1,z,1)/   &
( 1155)             sdot(noeuds,nsp(1,j),1,nsp(1,j),1))
( 1156) !       write(0,*)'  residu de la colonne ',j,' du noyau :',res
( 1157) !       write(0,*)'  vecteur colonne :'
( 1158) !       call prvec(noeuds,nsp(1,j))
( 1159)       end do
( 1160) 
( 1161)       return
( 1162)       end
( 1163) !**********************************************************************c
( 1164)       subroutine feti_project(gint,pgint,nspdim,x,b,nmspdim,numit0   &
( 1165)           ,nitmax,gvec,agvec,d,ad,add,gamm,numnes,listnes,plistin,numins   &
( 1166)           ,listin,narea,numnos,nsp,v,w,bufin,bufout,work)
( 1167) 
( 1168)       implicit none 
( 1169) 
( 1170)       integer nitmax, nspdim
( 1171)       REAL(kind=8)  x(nspdim), b(nspdim)
( 1172) 
( 1173) !...the Krylov optimisation used in PCG associated to the Coarse Solver
( 1174) !   involve to keep di & adi vectors in etesolve procedure
( 1175) !   the nmspdim parameter is introduced to avoid message "out of bound"
( 1176) !   arrising qhen the memory is checked!
( 1177) 
( 1178) !   nmspdim  =   Max (nspdim(Ep)) <= memory management (f77 context)
( 1179) !               p=1,Np
( 1180) 
( 1181) !      REAL(kind=8)  gvec(nspdim), agvec(nspdim),
( 1182) !     &       d(nspdim,0:nitmax-1), ad(nspdim,0:nitmax-1), 
( 1183) !     &       add(0:nitmax-1), gamm(0:nitmax-1)
( 1184) 
( 1185)       integer nmspdim
( 1186)       REAL(kind=8)  gvec(nspdim), agvec(nspdim),   &
( 1187)           d(nmspdim,0:nitmax-1), ad(nmspdim,0:nitmax-1),    &
( 1188)           add(0:nitmax-1), gamm(0:nitmax-1)
( 1189) 
( 1190)       integer numnes,numins,narea,  numnos, numit0
( 1191)       integer listnes(numnes,3), plistin(numnes+1), listin(numins)
( 1192)       REAL(kind=8)  nsp(numnos,nspdim), v(numnos),    &
( 1193)           w(numins), gint(numins), pgint(numins) 
( 1194)       REAL(kind=8)  bufin(numins), bufout(numins)
( 1195)       REAL(kind=8)  work(0:nitmax-1)
( 1196) 
( 1197) 
( 1198)       call feti_proet(numnes,listnes,plistin,numins,listin,numnos,   &
( 1199)           gint,nspdim,nsp,b,v)
( 1200)       call feti_etesolv(nspdim,x,b,nmspdim,numit0,nitmax,gvec,agvec,d,ad   &
( 1201)           ,add,gamm,numnes,listnes,plistin,numins,listin,narea,numnos   &
( 1202)           ,nsp,v,w,bufin,bufout,work)
( 1203)       call feti_proe(numnes,listnes,plistin,numins,listin,narea,numnos,   &
( 1204)           nspdim,nsp,x,v,w,bufin,bufout)
( 1205)       call feti_vsub(numins,gint,w,pgint)
( 1206) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 22

( 1207)       return
( 1208)       end
( 1209) !**********************************************************************c
( 1210)       subroutine feti_etesolv(nspdim,x,b,nmspdim,numit0,nitmax,gvec    &
( 1211)           ,agvec,d,ad,add,gamm,numnes,listnes,plistin,numins,listin    &
( 1212)           ,narea,numnos,nsp,v,w,bufin,bufout,work)
( 1213) 
( 1214) !     INPUTS :   nspdim  =   dimension of null space : LOCAL
( 1215) !                x       =   initial solution
( 1216) !                b       =   right hand side
( 1217) !               nmspdim  =   Max (nspdim(Ep)) <= memory management
( 1218) !                            p=1,Np
( 1219) !               numit0   =   number of previous direction vectors, at most
( 1220) !                            nitmax
( 1221) !               nitmax   =   number of maximum direction vectors for 
( 1222) !                            reconjugation  
( 1223) !               gvec     =   gradient vector
( 1224) !               agvec    =   A * g
( 1225) !               d        =   direction vector
( 1226) !               ad       =   A * d
( 1227) !               add      =   dot products (Ad,d)
( 1228) !               gamm     =   GAMMA coefficients
( 1229) !               numnes   =   number of neighbouring subdomains
( 1230) !               listnes  =   list of neighbouring subdomains
( 1231) !               plistin  =   pointer of sublists of interface nodes
( 1232) !               numins   =   number of interface nodes
( 1233) !               listin   =   list of interface nodes
( 1234) !               narea    =   subdomain number
( 1235) !               numnos   =   number of nodes in the subdomain
( 1236) 
( 1237) !     WORKSPACE : bufin  =  contributions of local fields on interface
( 1238) !                 bufout =  contributions of outer fields on interface
( 1239) !                 work   =  workspace for mpp_sum calls
( 1240) 
( 1241)       implicit none 
( 1242) 
( 1243)       integer nitmax, nspdim
( 1244)       REAL(kind=8)  x(nspdim), b(nspdim)
( 1245) 
( 1246) !...the Krylov optimisation involve to keep di & adi vectors
( 1247) !   the nmspdim parameter is introduce to avoid message "out of bound"
( 1248) !   arrising qhen the memory is checked!
( 1249) 
( 1250) !      REAL(kind=8)  gvec(nspdim), agvec(nspdim),
( 1251) !     &       d(nspdim,0:nitmax-1), ad(nspdim,0:nitmax-1), 
( 1252) !     &       rho, temp1(2), temp2(2), add(0:nitmax-1),
( 1253) !     &       gamm(0:nitmax-1), facg0, facgn, facst, eps
( 1254) 
( 1255)       integer nmspdim
( 1256)       REAL(kind=8)  gvec(nmspdim), agvec(nmspdim),    &
( 1257)           d(nspdim,0:nitmax-1), ad(nspdim,0:nitmax-1),     &
( 1258)           rho, temp1(2), temp2(2), add(0:nitmax-1),    &
( 1259)           gamm(0:nitmax-1), facg0, facgn, facst, eps
( 1260) 
( 1261)       integer numnes,numins,narea,  numnos
( 1262)       integer listnes(numnes,3), plistin(numnes+1), listin(numins)
( 1263)       REAL(kind=8)  nsp(numnos,nspdim), v(numnos), w(numins)
( 1264)       REAL(kind=8)  bufin(numins), bufout(numins)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 23

( 1265)       REAL(kind=8)  work(0:nitmax-1)
( 1266) 
( 1267)       integer i, j, n, numit0, nn, n1, numit1
( 1268) 
( 1269)       external sdot
( 1270)       REAL(kind=8)  sdot
( 1271) 
( 1272) !    initialisation 
( 1273) 
( 1274)       do 1 i=1,nspdim
( 1275)         x(i) = 0.e0
( 1276)         gvec(i) = -b(i)
( 1277)  1    continue
( 1278)       facg0 = sdot( nspdim, gvec, 1, gvec , 1)
( 1279)       call mpp_sum(facg0,1,work)
( 1280) 
( 1281)       if( facg0  ==  0.e0 ) return
( 1282) 
( 1283)       eps=1.e-24
( 1284)       facst=facg0*eps
( 1285) 
( 1286) !    initial projection 
( 1287) 
( 1288)       if(numit0 > 0) then
( 1289)           do 91 j = 0 , numit0
( 1290)             gamm(j) = - sdot(nspdim,gvec,1,d(1,j),1) / add(j)
( 1291)  91       continue
( 1292)           call mpp_sum(gamm,numit0+1,work)
( 1293)           call feti_mxvadd( nspdim, numit0+1, d, gamm, x )
( 1294)           call feti_mxvadd( nspdim, numit0+1, ad, gamm, gvec )
( 1295) 
( 1296) !  test residual after initial projection
( 1297) 
( 1298)           facgn = sdot( nspdim, gvec, 1, gvec, 1 )
( 1299)           call mpp_sum(facgn,1,work)
( 1300)           if( facgn <= facst ) then
( 1301) !             write(0,*) 'residual after initial projection for feti_etesolv='   &
( 1302) !                 ,facgn/facg0
( 1303)               return
( 1304)           endif
( 1305) 
( 1306)       endif
( 1307) 
( 1308) !     compute (Et E) * gvec
( 1309) 
( 1310)       call feti_proe(numnes,listnes,plistin,numins,listin,narea,numnos,   &
( 1311)           nspdim,nsp,gvec,v,w,bufin,bufout)
( 1312)       call feti_proet(numnes,listnes,plistin,numins,listin,numnos,w,   &
( 1313)           nspdim,nsp,agvec,v)
( 1314) 
( 1315) !   conjugate d .
( 1316) 
( 1317)       if(numit0 > 0) then
( 1318)           numit1=min0(numit0+1,nitmax-1)
( 1319)           do 94 i=0,numit0
( 1320)             gamm(i) = -sdot(nspdim,gvec,1,ad(1,i),1)/ add(i)
( 1321)  94       continue
( 1322)           call mpp_sum(gamm,numit0+1,work)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 24

( 1323) 
( 1324)           call feti_osaxpy( nspdim,gamm(numit0),d(1,numit0),   &
( 1325)               gvec,d(1,numit1))
( 1326)           call feti_osaxpy( nspdim,gamm(numit0),ad(1,numit0),   &
( 1327)               agvec,ad(1,numit1))
( 1328) 
( 1329)           call feti_mxvadd( nspdim, numit1-1, d, gamm, d(1,numit1) )
( 1330)           call feti_mxvadd( nspdim, numit1-1, ad, gamm, ad(1,numit1) )
( 1331)       else
( 1332)           numit1=0
( 1333) 
( 1334)           do 20 i = 1, nspdim
( 1335)             d(i,numit1) = gvec(i)
( 1336)  20       continue
( 1337) 
( 1338)           do 30 i = 1, nspdim
( 1339)             ad(i,numit1) = agvec(i)
( 1340)  30       continue
( 1341)       endif
( 1342) 
( 1343) !     computing rho
( 1344) 
( 1345)       do 50 nn = numit1, nitmax
( 1346)         n=min0(nn,nitmax-1)
( 1347)         temp1(1) = sdot( nspdim, gvec, 1, d(1,n), 1 )
( 1348)         temp1(2) = sdot( nspdim, ad(1,n), 1, d(1,n), 1 )
( 1349)         call mpp_sum(temp1,2,temp2)
( 1350)         add(n) = temp1(2)
( 1351)         rho = - temp1(1) / add(n)
( 1352) 
( 1353)         call feti_osaxpy( nspdim, rho, d(1,n), x, x )
( 1354)         call feti_osaxpy( nspdim, rho, ad(1,n), gvec, gvec )
( 1355) 
( 1356) !        test residual . 
( 1357) 
( 1358)         facgn = sdot( nspdim, gvec, 1, gvec, 1 )
( 1359)         call mpp_sum(facgn,1,work)
( 1360) !        write(0,*) 'iteration =', nn+1, 'residual = ', facgn/facg0
( 1361)         if( facgn <= facst ) then
( 1362)             goto 999
( 1363)         endif
( 1364) 
( 1365) !   reconjugation of gvec
( 1366) 
( 1367)         do 60 i=0,n
( 1368)           gamm(i) = -sdot(nspdim,gvec,1,ad(1,i),1)/ add(i)
( 1369)  60     continue
( 1370)         call mpp_sum(gamm,n+1,work)
( 1371)         n1=min0(nn+1,nitmax-1)
( 1372) 
( 1373)         call feti_osaxpy( nspdim,gamm(n),d(1,n),gvec,d(1,n1))
( 1374) 
( 1375)         call feti_mxvadd( nspdim, n1-1, d, gamm, d(1,n1) )
( 1376) 
( 1377) !     compute (Et E) * d
( 1378) 
( 1379)         call feti_proe(numnes,listnes,plistin,numins,listin,narea,   &
( 1380)             numnos,nspdim,nsp,d(1,n1),v,w,bufin,bufout)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 25

( 1381)         call feti_proet(numnes,listnes,plistin,numins,listin,numnos,w,   &
( 1382)             nspdim,nsp,ad(1,n1),v)
( 1383) 
( 1384)  50   continue
( 1385) 
( 1386)       if(narea == 1) then
( 1387)           write(0,*) 'No convergence for ete in ',nitmax,' iterations'
( 1388)           stop
( 1389)       endif
( 1390) 
( 1391) 
( 1392)  999  continue
( 1393)       numit0=min0(n,nitmax-1)
( 1394) !     write(0,*) 'number of iterations for ete :', nn+1, ' ,  residual =', facgn
( 1395) 
( 1396)       return
( 1397)       end
( 1398) !**********************************************************************c
( 1399)       subroutine feti_osaxpy( n, a, x, y, z )
( 1400) 
( 1401)       implicit none
( 1402) 
( 1403)       integer n
( 1404)       REAL(kind=8)  a, x(n), y(n), z(n)
( 1405) 
( 1406)       integer i
( 1407) 
( 1408)       do 10 i = 1, n
( 1409)         z(i) = y(i) + a * x(i)
( 1410)  10   continue
( 1411) 
( 1412)       return
( 1413)       end
( 1414) !**********************************************************************c
( 1415)       subroutine feti_proe(numnes,listnes,plistin,numins,listin,narea,   &
( 1416)           numnos,nspdim,nsp,alpha,v,w,bufin,   &
( 1417)           bufout)
( 1418) 
( 1419) !  jump of a velocity field on the interfaces, in the case where 
( 1420) !  subdomain number i is allocated to processor number i-1 .
( 1421) !  the velocity field is Nsp * alpha, where Nsp is the nullspace .
( 1422) 
( 1423) !     INPUTS :      numnes  =  number of neighbouring subdomains
( 1424) !               listnes =  list of neighbouring subdomains
( 1425) !               plistin =  pointer of sublists of interface nodes
( 1426) !               numins  =  number of interface nodes
( 1427) !               listin  =  list of interface nodes
( 1428) !               narea     =  subdomain number
( 1429) !               numnos  =  number of nodes in the subdomain
( 1430) !               nspdim  =  dimension of the null space
( 1431) !               nsp     =  null space
( 1432) !               alpha   =  components of the velocity field in the null
( 1433) !                          space
( 1434) !               v       =  corresponding velocity field
( 1435) 
( 1436) !     OUTPUTS :      w  =  jump on the interface
( 1437) 
( 1438) !     WORKSPACE : bufin   =  contributions of local field on interface






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 26

( 1439) !                 bufout  =  contributions of outer fields on interface
( 1440) !                 work
( 1441)       implicit none 
( 1442)       integer numnes,numins,numnos,narea,nspdim
( 1443)       integer listnes(numnes,3), plistin(numnes+1), listin(numins)
( 1444)       REAL(kind=8)  v(numnos),w(numins)
( 1445)       REAL(kind=8)  bufin(numins),bufout(numins)
( 1446)       REAL(kind=8)  alpha(nspdim), nsp(numnos,nspdim)
( 1447)       integer i
( 1448) 
( 1449) !  computing the values of the local velocity field 
( 1450)       do i=1,numnos
( 1451)         v(i)=0.e0
( 1452)       end do
( 1453)       call feti_mxvadd( numnos, nspdim, nsp, alpha, v )
( 1454) !  computing the jump on the interface
( 1455)       call feti_saut(numnes,listnes,plistin,numins,listin,   &
( 1456)                 numnos,v,w,bufin,bufout)
( 1457) !!        narea,numnos,v,w,bufin,bufout)
( 1458) 
( 1459)       return
( 1460)       end
( 1461) !**********************************************************************c
( 1462)       subroutine feti_proet(numnes,listnes,plistin,numins,listin,   &
( 1463)           numnos,w,nspdim,nsp,alpha,v)
( 1464) 
( 1465) !  projection in the null space of the right-hand-side associated with
( 1466) !  a Lagrange multiplier w .
( 1467) 
( 1468) !     INPUTS :      numins  =  number of interface nodes
( 1469) !               listin  =  list of interface nodes
( 1470) !               numnos  =  number of nodes in the subdomain
( 1471) !               w       =  Lagrange multiplier
( 1472) !               nspdim  =  dimension of the null space
( 1473) !               nsp     =  null space
( 1474) 
( 1475) !     OUTPUTS :      v     =  local right-hand-side
( 1476) !               alpha =  projection of v in the null space
( 1477) 
( 1478)       implicit none 
( 1479)       external sdot
( 1480)       REAL(kind=8)  sdot
( 1481)       integer numins,numnos,nspdim,numnes
( 1482)       integer listnes(numnes,3),plistin(numnes+1),listin(numins)
( 1483)       REAL(kind=8)  v(numnos),w(numins)
( 1484)       REAL(kind=8)  alpha(nspdim), nsp(numnos,nspdim)
( 1485)       integer i
( 1486) 
( 1487) !  assembling the values of the Lagrange multiplier on the interface
( 1488)       call feti_probit(numnes,listnes,plistin,numins,listin,numnos,w,v)
( 1489) !  computing the projection of v in the null space .
( 1490)       do 1 i=1,nspdim
( 1491)         alpha(i)=sdot(numnos,v,1,nsp(1,i),1)
( 1492)  1    continue
( 1493) 
( 1494)       return
( 1495)       end
( 1496) !**********************************************************************c






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 27

( 1497)       subroutine feti_consrhs( vdim,nspdim,rhs,nsp,cb )
( 1498) 
( 1499) 
( 1500)       implicit none
( 1501) 
( 1502)       external sdot
( 1503)       REAL(kind=8)  sdot
( 1504)       integer vdim, nspdim
( 1505) 
( 1506)       REAL(kind=8)   rhs(vdim), nsp(vdim,nspdim), cb(nspdim)
( 1507) 
( 1508) 
( 1509)       integer i
( 1510) 
( 1511) !  compute the right hand side for constraint 
( 1512) 
( 1513)       do i=1,nspdim
( 1514)         cb(i)=-sdot(vdim,rhs,1,nsp(1,i),1)
( 1515)       end do
( 1516) 
( 1517)       return
( 1518)       end
( 1519) !**********************************************************************c
( 1520)       subroutine feti_dualschur(noeuds,ni,nj,a,lpblo,blo,   &
( 1521)           ninterf,ninterfc,nni,nnic,   &
( 1522)           ndvois,ndvoisc,plistin,listin,   &
( 1523)           poids,u,v,f,bitw,utilu,   &
( 1524)           lambda,g,pg,mg,nitmax,nmaxd,j0,wj,dwj,dwwj,   &
( 1525)           gamm,work,bufin,bufout,narea,epsilo,ndlblo,   &
( 1526)           lisblo,ndkerep,   &
( 1527)           xnul,ynul,numit0ete,nitmaxete,eteg,   &
( 1528)           eteag,eted,etead,eteadd,etegamm,nsp,etev,   &
( 1529)           etew,nnih,plistih,gh,w,dw,   &
( 1530)           residu,indic,jn)
( 1531)       implicit none
( 1532) !  nombre de ddl par noeuds, nombre de noeuds interne et interface
( 1533)       integer noeuds,ni,nj,nni,nnic
( 1534) !      integer jpj,jpi
( 1535) !  numero du sous-domaine = numero de processeur + 1 
( 1536)       integer narea
( 1537) !  tableaux descripteurs de la matrice locale
( 1538)       REAL(kind=8)  a(ni,nj,5)
( 1539) !  tableaux descripteurs de l'inverse de la matrice du probleme local
( 1540)       integer lpblo
( 1541)       REAL(kind=8)  blo(lpblo)
( 1542) !  tableaux descripteurs de l'interface
( 1543)       integer ninterf,ninterfc
( 1544)       integer ndvois(ninterf,3),ndvoisc(ninterfc,3)
( 1545)       integer plistin(ninterfc+1),listin(nnic)
( 1546) !  tableaux descripteurs de l'interface ajoutes pour nperio == 1
( 1547)       REAL(kind=8)  poids(noeuds)
( 1548) !  utilitaires pour les assemblages aux interfaces
( 1549)       REAL(kind=8)  bufin(nnic),bufout(nnic)
( 1550) !  vecteurs locaux
( 1551)       REAL(kind=8)  u(noeuds),v(noeuds),f(noeuds)
( 1552)       REAL(kind=8)  bitw(noeuds),utilu(noeuds)
( 1553) !  vecteurs inconnues aux interfaces
( 1554) !  mg est le gradient preconditionne : M = sum(Bi Ai Bit)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 28

( 1555)       REAL(kind=8)  lambda(nni),g(nni),pg(nni),mg(nni)
( 1556) !  tableaux servant a gerer le stockage de la moitie des coefficients
( 1557) !  des tableaux des directions de descente sur chaque interface
( 1558)       integer nnih
( 1559)       integer plistih(ninterf+1)
( 1560)       REAL(kind=8)  gh(nnih)
( 1561) !  tableaux des directions de descente
( 1562) !  on appelle D la matrice condensee a l'interface : sum(Bi [Ai]-1 Bit)
( 1563) !  j0 est le nombre de vecteurs D orthogonaux deja stockes
( 1564) !  j0 est egal au plus a nmaxd-1
( 1565) !  lorsque le nomtre total d'iterations atteint nmaxd, seule la derniere
( 1566) !  direction calculee est stockee dans la colonne numero nmaxd
( 1567)       integer nitmax,nmaxd,j0
( 1568)       REAL(kind=8)  wj(nnih,nmaxd),dwj(nnih,nmaxd),dwwj(nmaxd)
( 1569)       REAL(kind=8)  w(nni),dw(nni)
( 1570) !  utilitaires pour les calculs de produits scalaires globaux
( 1571)       REAL(kind=8)  work(nmaxd),gamm(nmaxd)
( 1572) !  tableaux descripteurs du noyau de la matrice locale et tableaux 
( 1573) !  servant a la projection sur le noyau
( 1574) 
( 1575) !...ndlblo is the dimension of the local nullspace .=<. the size of the
( 1576) !   memory of the superstructure associated to the nullspace : ndkerep
( 1577) !   indeed ndkerep = Max ndlblo = Max dim(Ker(Ep))
( 1578) !                                p=1,Np
( 1579) !   ndkerep is introduced to avoid messages "out of bounds" when memory
( 1580) !   is checked
( 1581) 
( 1582)       integer ndkerep
( 1583)       integer ndlblo,numit0ete,nitmaxete
( 1584)       integer lisblo(ndlblo)
( 1585)       REAL(kind=8)  xnul(ndlblo),ynul(ndlblo),eteg(ndlblo),eteag(ndlblo)
( 1586) 
( 1587) !...the Krylov optimisation used in PCG associated to the Coarse Solver
( 1588) !   involve to keep di & adi vectors in etesolve procedure
( 1589) !   the ndkerep parameter is introduced to avoid message "out of bound"
( 1590) !   arrising when the memory is checked!
( 1591) 
( 1592) !   ndkerep  =   Max (ndlblo(Ep)) <= memory management (f77 context)
( 1593) !               p=1,Np
( 1594) 
( 1595) !      REAL(kind=8)  eted(ndlblo,nitmaxete),etead(ndlblo,nitmaxete)
( 1596) 
( 1597)       REAL(kind=8)  eted(ndkerep,nitmaxete),etead(ndkerep,nitmaxete)
( 1598)       REAL(kind=8)  eteadd(nitmaxete),etegamm(nitmaxete)
( 1599)       REAL(kind=8)  etev(noeuds),etew(nni),nsp(noeuds,ndlblo)
( 1600) !  indicateur de convergence
( 1601) 
( 1602)       integer indic
( 1603)       REAL(kind=8)  residu
( 1604)       integer jn
( 1605) !  utilitaires
( 1606) !!    integer j,k,j1,ji,jj,jk
( 1607)       integer j,k,j1
( 1608)       REAL(kind=8)  scal(2),ff,roj,epsilo,epsilo2
( 1609) !  fonctions externes
( 1610)       external sdot
( 1611)       REAL(kind=8)  sdot
( 1612) !  booleen pour caracterisation






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 29

( 1613) !  numerique si llconv=1
( 1614)       integer llconv
( 1615) 
( 1616)       llconv=0
( 1617) 
( 1618) !...critere d'arret epsilo2 = epsilo **2 => evite sqrt
( 1619) 
( 1620)       epsilo2 = epsilo * epsilo
( 1621) 
( 1622) !  calcul du carre scalaire du second membre global
( 1623) 
( 1624)       call feti_assemb(ninterfc,ndvoisc,plistin,nnic,listin,narea,noeuds   &
( 1625)           ,f,utilu,bufin,bufout)
( 1626)       ff=sdot(noeuds,utilu,1,utilu,1)
( 1627)       call mpp_sum(ff,1,work)
( 1628) !     write(90,*)'valeur du carre du second membre ds le Schur Dual :',ff
( 1629) 
( 1630) !  calcul du second membre contraint 
( 1631) 
( 1632)       call feti_consrhs(noeuds,ndlblo,f,nsp,ynul)
( 1633) 
( 1634) !  calcul de lambda0 satisfaisant la contrainte
( 1635)  
( 1636) !     call feti_vclr(nni,lambda)
( 1637)       call feti_etesolv(ndlblo,xnul,ynul,ndkerep,numit0ete,nitmaxete   &
( 1638)           ,eteg,eteag,eted,etead,eteadd,etegamm,ninterf,ndvois,plistin   &
( 1639)           ,nni,listin,narea,noeuds,nsp,etev,etew,bufin,bufout,work)
( 1640) 
( 1641)       call feti_proe(ninterf,ndvois,plistin,nni,listin,narea,noeuds,   &
( 1642)           ndlblo,nsp,xnul,etev,lambda,bufin,bufout)
( 1643) 
( 1644) !  calcul du nouveau second membre
( 1645) 
( 1646)       call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds,lambda,bitw)
( 1647) 
( 1648)       call feti_vadd(noeuds,bitw,f,bitw)
( 1649) 
( 1650) 
( 1651) !  calcul du champ local initial u0
( 1652) !     call feti_blodir(noeuds,bitw,ndlblo,lisblo)
( 1653)       call feti_resloc(noeuds,ni,nj,a,lpblo,blo,bitw,u,utilu)
( 1654) 
( 1655) !  calcul du gradient initial g0 = feti_saut(u0)
( 1656) !!    call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds,u,g   &
( 1657)       call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds,u,g   &
( 1658)           ,bufin,bufout)
( 1659) 
( 1660) !  reinitialisation de lambdaj0 par reconjugaison
( 1661) 
( 1662)       call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea,g, gh)
( 1663)       do k=1,j0
( 1664)         gamm(k)=sdot(nnih,gh,1,wj(1,k),1)
( 1665)       end do
( 1666)       if(j0 >= 1) call mpp_sum(gamm,j0,work)
( 1667)       do k=1,j0
( 1668)         gamm(k)=-gamm(k)/dwwj(k)
( 1669)       end do
( 1670)       call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, lambda,gh)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 30

( 1671)       call feti_mxvadd(nnih,j0,wj,gamm,gh)
( 1672)       call feti_extend(ninterf,ndvois,plistin,nni,plistih,nnih,narea,gh, lambda)
( 1673) !  calcul du second membre associe a lambdaj0
( 1674)       if(j0 > 0) then
( 1675)           call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds, lambda,bitw)
( 1676)           call feti_vadd(noeuds,bitw,f,bitw)
( 1677) !  calcul du champ local associe uj0
( 1678) !         call feti_blodir(noeuds,bitw,ndlblo,lisblo)
( 1679)           call feti_resloc(noeuds,ni,nj,a,lpblo,blo,bitw,u,utilu)
( 1680) !  calcul du gradient initial gj0= feti_saut(uj0)
( 1681) !!        call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds   &
( 1682)           call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds   &
( 1683)               ,u,g,bufin,bufout)
( 1684)       endif
( 1685) !  calcul du gradient projete 
( 1686) !     call feti_vmov(nni,g,pg)
( 1687)       call feti_project(g,pg,ndlblo,xnul,ynul,ndkerep,numit0ete   &
( 1688)           ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf,ndvois   &
( 1689)           ,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin,bufout,work)
( 1690) !  calcul du gradient preconditionne 
( 1691) !     call feti_vmov(nni,pg,mg)
( 1692) 
( 1693)       call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds, pg,bitw)
( 1694) 
( 1695)       call feti_proax(noeuds,ni,nj,a,bitw,v) 
( 1696) !!    call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds,v, mg,bufin,bufout)
( 1697)       call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds,v, mg,bufin,bufout)
( 1698) !  verification du residu initial 
( 1699)       residu=sdot(nni,mg,1,mg,1)
( 1700)       call mpp_sum(residu,1,work)
( 1701)       residu = (4. * residu) / ff
( 1702)       if(residu <= epsilo2)then
( 1703) !         residu=10.*sqrt(residu/ff)
( 1704) !         if(residu < epsilo)then
( 1705)           if(llconv == 1)then
( 1706)               if(narea == 1)then
( 1707)                   write(0,*)'  residu carre approche initial :',residu
( 1708)               endif
( 1709)           endif
( 1710) !  calcul du champ global moyenne aux interfaces et du residu global 
( 1711)           call feti_vsub(noeuds,u,etev,utilu)
( 1712)           call feti_moyen(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
( 1713)               noeuds,poids,utilu,utilu,bufin,bufout)
( 1714)           call feti_vmov(noeuds,utilu,u)
( 1715)           call feti_proax(noeuds,ni,nj,a,utilu,v) 
( 1716)           call feti_vsub(noeuds,v,f,v)
( 1717)           call feti_assemb(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
( 1718)               noeuds,v,utilu,bufin,bufout)
( 1719)           residu=sdot(noeuds,utilu,1,utilu,1)
( 1720)           call mpp_sum(residu,1,work)
( 1721)           residu=sqrt(residu/ff)
( 1722) !         residu=residu/ff
( 1723)           if(llconv == 1)then
( 1724)               if(narea == 1)then
( 1725)                   write(0,*)'  residu exact initial :',residu
( 1726)               endif
( 1727)           endif
( 1728) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 31

( 1729)  101      format(10d9.2)
( 1730) 
( 1731)           return
( 1732) 
( 1733)       endif
( 1734) 
( 1735) 
( 1736) !  calcul de la premiere direction de descente
( 1737) !  calcul du gradient preconditionne projete 
( 1738)       call feti_project(mg,mg,ndlblo,xnul,ynul,ndkerep,numit0ete   &
( 1739)           ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf,ndvois   &
( 1740)           ,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin,bufout ,work)
( 1741)       call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea,   &
( 1742)           mg,wj(1,j0+1))
( 1743)       do k=1,j0
( 1744)         gamm(k)=sdot(nnih,wj(1,j0+1),1,dwj(1,k),1)
( 1745)       end do
( 1746)       if(j0 >= 1) call mpp_sum(gamm,j0,work)
( 1747)       do k=1,j0
( 1748)         gamm(k)=-gamm(k)/dwwj(k)
( 1749)       end do
( 1750)       call feti_mxvadd(nnih,j0,wj,gamm,wj(1,j0+1))
( 1751) 
( 1752) !  iterations de gradient conjugue a l'interface
( 1753) 
( 1754) 
( 1755) !debug
( 1756)       do j=j0+1,j0+nitmax
( 1757) 
( 1758)         if(llconv == 1.and.narea == 1) write(0,*)'etape numero',j-j0
( 1759) !  j1 est le numero de colonne de la derniere direction de descente
( 1760) !  j1 = j, numero d'iteration total, tant que j < nmaxd
( 1761) !  au-dela, j1 = nmaxd
( 1762)         j1=min0(j,nmaxd)
( 1763) !  calcul de Dwi , D = sum ( Bi [Ai]-1 tBi )
( 1764)         call feti_extend(ninterf,ndvois,plistin,nni,plistih,nnih,narea, wj(1,j1),w)
( 1765)         call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds, w,bitw)
( 1766) !       call feti_blodir(noeuds,bitw,ndlblo,lisblo)
( 1767)         call feti_resloc(noeuds,ni,nj,a,lpblo,blo,bitw,v,utilu)
( 1768) !!      call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds,v   &
( 1769)         call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds,v   &
( 1770)             ,dw,bufin,bufout)
( 1771)         call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, dw,dwj(1,j1))
( 1772) !  calcul du coefficient de descente
( 1773)         call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, g,gh)
( 1774)         scal(1)=sdot(nnih,gh,1,wj(1,j1),1)
( 1775)         scal(2)=sdot(nnih,dwj(1,j1),1,wj(1,j1),1)
( 1776)         call mpp_sum(scal,2,work)
( 1777)         roj=-scal(1)/scal(2)
( 1778)         dwwj(j1)=scal(2)
( 1779) !  remise a jour de lambda, g et u
( 1780)         call saxpy(nni,roj,w,1,lambda,1)
( 1781)         call saxpy(nni,roj,dw,1,g,1)
( 1782)         call saxpy(noeuds,roj,v,1,u,1)
( 1783) !  calcul du gradient projete 
( 1784) !       call feti_vmov(nni,g,pg)
( 1785)         call feti_project(g,pg,ndlblo,xnul,ynul,ndkerep,numit0ete   &
( 1786)             ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf   &






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 32

( 1787)             ,ndvois,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin   &
( 1788)             ,bufout,work)
( 1789) !  calcul du gradient preconditionne Mg , M = sum ( Bi Ai tBi )
( 1790) !       call feti_vmov(nni,pg,mg)
( 1791)         call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds, pg,bitw)
( 1792)         call feti_proax(noeuds,ni,nj,a,bitw,v) 
( 1793) !!      call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds,v   &
( 1794)         call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds,v   &
( 1795)             ,mg,bufin,bufout)
( 1796) !  verification du residu global approche
( 1797)         residu=sdot(nni,mg,1,mg,1)
( 1798)         call mpp_sum(residu,1,work)
( 1799) !       residu=10.*sqrt(residu/ff)
( 1800)         residu = (4. * residu) / ff
( 1801)         if(llconv == 1.and.narea == 1)then
( 1802)             write(0,*)'  residu carre global approche apres ',   &
( 1803)                 j-j0,' iterations :',residu
( 1804)         endif
( 1805) 
( 1806) !       if(residu <  epsilo) then
( 1807)         if(residu <= epsilo2)then
( 1808)             if(llconv == 1.and.narea == 1)then
( 1809)                 write(0,*)'residu carre global approche apres ',   &
( 1810)                     j-j0,' iterations :',residu
( 1811)             endif
( 1812) !  calcul du champ global moyenne aux interfaces et du residu global 
( 1813)             call feti_vsub(noeuds,u,etev,utilu)
( 1814)             call feti_moyen(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
( 1815)                 noeuds,poids,utilu,utilu,bufin,bufout)
( 1816)             call feti_vmov(noeuds,utilu,u)
( 1817)             call feti_proax(noeuds,ni,nj,a,utilu,v) 
( 1818)             call feti_vsub(noeuds,v,f,v)
( 1819)             call feti_assemb(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
( 1820)                 noeuds,v,utilu,bufin,bufout)
( 1821)             residu=sdot(noeuds,utilu,1,utilu,1)
( 1822)             call mpp_sum(residu,1,work)
( 1823) !           residu=residu/ff
( 1824)             residu=sqrt(residu/ff)
( 1825)             if(llconv == 1.and.narea == 1)then
( 1826)                 write(0,*)'  residu global exact apres ',   &
( 1827)                     j-j0,' iterations :',residu
( 1828)             endif
( 1829) 
( 1830)             jn = j - j0
( 1831)             j0=min0(j,nmaxd-1)
( 1832)             if(llconv == 1.and.narea == 1)then
( 1833)                 write(0,*)'  nombre de directions orthogonales conservees :',j0
( 1834)             endif
( 1835) 
( 1836) 
( 1837)             return
( 1838)         endif
( 1839) !##################################################################
( 1840) !  calcul du champ global moyenne aux interfaces et du residu global 
( 1841)         call feti_vsub(noeuds,u,etev,utilu)
( 1842)         call feti_moyen(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
( 1843)             noeuds,poids,utilu,utilu,bufin,bufout)
( 1844)         call feti_proax(noeuds,ni,nj,a,utilu,v) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 33

( 1845)         call feti_vsub(noeuds,v,f,v)
( 1846)         call feti_assemb(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
( 1847)             noeuds,v,utilu,bufin,bufout)
( 1848)         residu=sdot(noeuds,utilu,1,utilu,1)
( 1849)         call mpp_sum(residu,1,work)
( 1850) !       residu=sqrt(residu/ff)
( 1851)         residu=residu/ff
( 1852)         if(llconv == 1.and.narea == 1)then
( 1853)             write(0,*)'  residu carre global exact apres ',   &
( 1854)                 j-j0,' iterations :',residu
( 1855)         endif
( 1856) 
( 1857) !#####################################################################
( 1858) !  calcul de la nouvelle direction de descente par reconjugaison
( 1859) !  calcul du gradient preconditionne projete 
( 1860)         call feti_project(mg,mg,ndlblo,xnul,ynul,ndkerep,numit0ete   &
( 1861)             ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf   &
( 1862)             ,ndvois,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin   &
( 1863)             ,bufout,work)
( 1864)         call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, mg,gh)
( 1865)         do k=1,j1
( 1866)           gamm(k)=sdot(nnih,gh,1,dwj(1,k),1)
( 1867)         end do
( 1868)         call mpp_sum(gamm,j1,work)
( 1869)         do k=1,j1
( 1870)           gamm(k)=-gamm(k)/dwwj(k)
( 1871)         end do
( 1872)         call feti_mxvadd(nnih,j1,wj,gamm,gh)
( 1873)         j1=min0(j1+1,nmaxd)
( 1874) 
( 1875) !%%%%%%%%%%stabilisation numerique %%%%%%%%%%%%%%%%%%%%%%%%%
( 1876)         call feti_extend(ninterf,ndvois,plistin,nni,plistih,nnih,narea, gh,w)
( 1877)         call feti_project(w,w,ndlblo,xnul,ynul,ndkerep,numit0ete   &
( 1878)             ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf   &
( 1879)             ,ndvois,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin   &
( 1880)             ,bufout,work)
( 1881)         call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, w,gh)
( 1882) !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
( 1883) 
( 1884)         call feti_vmov(nnih,gh,wj(1,j1))
( 1885) 
( 1886)       end do
( 1887) 
( 1888)       residu = sqrt(residu)
( 1889)       if(narea == 1) then
( 1890)           write(0,*)'  residu global approche apres ',nitmax,   &
( 1891)               '  iterations :',residu
( 1892)       endif
( 1893)       indic=-2
( 1894) 
( 1895)       return
( 1896)       end
( 1897) !**********************************************************************c
( 1898)       subroutine feti_extend(numnes,listnes,   &
( 1899)           plistin,numins,plistih,   &
( 1900)           numinh,narea,wh,w)
( 1901) 
( 1902) !  reassembling an half interface field 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 34

( 1903) 
( 1904) !     INPUTS :      numnes  =  number of neighbouring subdomains
( 1905) !               listnes =  list of neighbouring subdomains
( 1906) !               lint    =  list of inner interfaces
( 1907) !               lext    =  list of extern ibnterfaces
( 1908) !               plistin =  pointer of sublists of interface nodes
( 1909) !               plistih =  pointer of half-sublists of interface nodes
( 1910) !               numins  =  number of interface nodes
( 1911) !               numinh  =  half number of interface nodes
( 1912) !               narea    =  subdomain number
( 1913) !               wh      =  half interface vector 
( 1914) 
( 1915) !     OUTPUTS :      w  =  complete interface vector
( 1916) 
( 1917) !     WORKSPACE : work
( 1918)       implicit none 
( 1919) !      external irecv
( 1920) !!    integer irecv, work(100)
( 1921)       integer numnes,numins,numinh,narea
( 1922)       integer listnes(numnes,3), plistin(numnes+1), plistih(numnes+1)
( 1923)       REAL(kind=8)  wh(numinh),w(numins)
( 1924)       integer i,l,lin,lout
( 1925)       integer mesg
( 1926) 
( 1927)       mesg=128*10
( 1928) !  if subdomains i < j are neighbours, subdomain i manage the first
( 1929) !  half of the interface, and subdomain j the remaining part
( 1930) 
( 1931) !  receiving the values of the outer fields on the interface
( 1932) !  sending the values of the inner field on the interface
( 1933) 
( 1934)       do i=1,numnes
( 1935)         l=(plistin(i+1)-plistin(i))
( 1936) !        if(narea < listnes(i)) then
( 1937)         if(listnes(i,3) < listnes(i,2))then
( 1938)             lout=(l-(l/2))
( 1939)             lin=(l/2)
( 1940) !           work(i)=irecv(mesg+listnes(i)-1,w(plistin(i)+(l/2)+1), lout*8)
( 1941) !           call mppsend(mesg+narea-1,wh(plistih(i)+1),lin*8,listnes(i) -1,0)
( 1942)             call mppsend(mesg+listnes(i,3)-1,wh(plistih(i)+1),lin*8, listnes(i,1)-1,0)
( 1943)             call feti_vmov(lin,wh(plistih(i)+1),w(plistin(i)+1))
( 1944)         else
( 1945)             lout=(l/2)
( 1946)             lin=(l-(l/2))
( 1947) !           work(i)=irecv(mesg+listnes(i)-1,w(plistin(i)+1),lout*8)
( 1948) !           call mppsend(mesg+narea-1,wh(plistih(i)+1),lin*8,listnes(i)-1,0)
( 1949)             call mppsend(mesg+listnes(i,3)-1,wh(plistih(i)+1),lin*8, listnes(i,1)-1,0)
( 1950)             call feti_vmov(lin,wh(plistih(i)+1),w(plistin(i)+(l/2)+1))
( 1951)         endif
( 1952)       enddo
( 1953) !  waiting for the completion on each interface 
( 1954)       do i=1,numnes
( 1955)         l=(plistin(i+1)-plistin(i))
( 1956) !       if(narea < listnes(i)) then
( 1957)         if(listnes(i,3) < listnes(i,2))then
( 1958)             lout=(l-(l/2))
( 1959) !           call mpprecv(mesg+listnes(i)-1,w(plistin(i)+(l/2)+1), lout*8)
( 1960)             call mpprecv(mesg+listnes(i,2)-1,w(plistin(i)+(l/2)+1), lout*8)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 35

( 1961)             call feti_vneg(lout,w(plistin(i)+(l/2)+1),   &
( 1962)                 w(plistin(i)+(l/2)+1))
( 1963)         else
( 1964)             lout=(l/2)
( 1965) !           call mpprecv(mesg+listnes(i)-1,w(plistin(i)+1),lout*8)
( 1966)             call mpprecv(mesg+listnes(i,2)-1,w(plistin(i)+1),lout*8)
( 1967)             call feti_vneg(lout,w(plistin(i)+1),w(plistin(i)+1))
( 1968)         endif
( 1969)       end do
( 1970) 
( 1971)       return
( 1972)       end
( 1973) !**********************************************************************c
( 1974)       subroutine feti_restri(numnes,listnes,plistin,numins,plistih,   &
( 1975)           numinh,narea,w,wh)
( 1976) 
( 1977) !  restriction of an interface field to one half on each interface
( 1978) 
( 1979) !     INPUTS :      numnes  =  number of neighbouring subdomains
( 1980) !               listnes =  list of neighbouring subdomains
( 1981) !               plistin =  pointer of sublists of interface nodes
( 1982) !               plistih =  pointer of half-sublists of interface nodes
( 1983) !               numins  =  number of interface nodes
( 1984) !               numinh  =  half number of interface nodes
( 1985) !               narea    =  subdomain number
( 1986) !               w       =  complete interface vector 
( 1987) 
( 1988) !     OUTPUTS :      wh  =  half interface vector
( 1989) 
( 1990) !     WORKSPACE : work
( 1991)       implicit none 
( 1992)       integer numnes,numins,numinh,narea
( 1993)       integer listnes(numnes,3), plistin(numnes+1), plistih(numnes+1)
( 1994)       REAL(kind=8)  wh(numinh),w(numins)
( 1995)       integer i,l,lin
( 1996) 
( 1997) !  if subdomains i < j are neighbours, subdomain i manage the first
( 1998) !  half of the interface, and subdomain j the remaining part
( 1999) 
( 2000)       do i=1,numnes
( 2001)         l=(plistin(i+1)-plistin(i))
( 2002) !       if(narea < listnes(i)) then
( 2003)         if(listnes(i,3) < listnes(i,2))then
( 2004)             lin=(l/2)
( 2005)             call feti_vmov(lin,w(plistin(i)+1),wh(plistih(i)+1))
( 2006)         else
( 2007)             lin=(l-(l/2))
( 2008)             call feti_vmov(lin,w(plistin(i)+(l/2)+1),wh(plistih(i)+1))
( 2009)         endif
( 2010)       end do
( 2011) 
( 2012)       return
( 2013)       end
( 2014) !**********************************************************************c
( 2015) !!    subroutine feti_halfint(numnes,listnes,plistin,numins,plistih   &
( 2016) !!        ,numinh,narea)
( 2017)       subroutine feti_halfint(numnes,listnes,plistin,       plistih   &
( 2018)           ,numinh      )






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 36

( 2019) 
( 2020) !  construction of the pointer of the restriction of an interface 
( 2021) !  field to one half on each interface
( 2022) 
( 2023) !     INPUTS :      numnes  =  number of neighbouring subdomains
( 2024) !               listnes =  list of neighbouring subdomains
( 2025) !               plistin =  pointer of sublists of interface nodes
( 2026) !               numins  =  number of interface nodes
( 2027) !               narea    =  subdomain number
( 2028) 
( 2029) !     OUTPUTS :      numinh  =  half number of interface nodes
( 2030) !               plistih =  pointer of half-sublists of interface nodes
( 2031) 
( 2032) !     WORKSPACE : work
( 2033)       implicit none 
( 2034) !!    integer numnes,numins,numinh,narea
( 2035)       integer numnes,       numinh
( 2036)       integer listnes(numnes,3), plistin(numnes+1), plistih(numnes+1)
( 2037)       integer i,l
( 2038) 
( 2039) !  if subdomains i < j are neighbours, subdomain i manage the first
( 2040) !  half of the interface, and subdomain j the remaining part
( 2041) 
( 2042)       plistih(1)=0
( 2043)       numinh=0
( 2044)       do i=1,numnes
( 2045)         l=(plistin(i+1)-plistin(i))
( 2046) !        if(narea < listnes(i,1)) then
( 2047)         if(listnes(i,3) < listnes(i,2)) then
( 2048) !          write(*,*)'je gere interface avec : ',listnes(i,1)
( 2049)             plistih(i+1)=plistih(i)+(l/2)
( 2050)         else
( 2051) !          write(*,*)'il gere interface avec : ',listnes(i,1)
( 2052)             plistih(i+1)=plistih(i)+(l-(l/2))
( 2053)         endif
( 2054)       end do
( 2055)       numinh=plistih(numnes+1)
( 2056) 
( 2057)       return
( 2058)       end
( 2059) !***********************************************************************
( 2060)       subroutine feti_front(n,ni,nj,a,lcmat,cmat,d,v,w,ndlblo,lisblo   &
( 2061)           ,nmdlblo)
( 2062)       implicit REAL(kind=8) (a-h,o-z)
( 2063)       dimension a(n,5)
( 2064) !  routine de calcul de la decomposition frontale par blocs de
( 2065) !  la matrice pentadiagonale a .
( 2066) !  cmat contient les termes des blocs diagonaux factorises .
( 2067)       dimension cmat(lcmat),d(n)
( 2068) 
( 2069) !...ndlblo nombre de liberte bloques = nombre de pivots nuls = dim(Ker(Ep))
( 2070) 
( 2071) !...nmdlblo = Max (ndlblo(Ep)) <= superstructure management in f77 context
( 2072) !            p=1,Np
( 2073) 
( 2074)       integer ndlblo,nmdlblo
( 2075)       integer lisblo(nmdlblo)
( 2076) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 37

( 2077) !  d ne contient pas, en sortie, l'inverse de la diagonale de Crout .
( 2078) !  v est un tableau servant a stocker la partie triangulaire inferieure
( 2079) !  stricte d'un bloc diagonal , et w un tableau servant a stocker un
( 2080) !  bloc sur-diagonal plein .
( 2081) !  de dimension egale a la largeur maximale d'un front .
( 2082)       dimension v(ni,ni),w(ni)
( 2083) !!!
( 2084)       ndimd=ni
( 2085)       lbd=ndimd*ndimd
( 2086) 
( 2087) !      lecture du tout premier bloc diagonal
( 2088) 
( 2089)       nb=1
( 2090) !  adresse du bloc diagonal courant
( 2091)       iac=1
( 2092)       call feti_liblod(n,ni,a,ndimd,nb,cmat(iac))
( 2093) !  inversion du premier bloc diagonal .
( 2094) !     write(*,*)'  inversion du premier bloc diagonal '
( 2095)       j0=0
( 2096)       call feti_ldlt(ndimd,cmat(iac),d(1),w,ndlblo,lisblo,nmdlblo,j0)
( 2097)       call feti_calinv(ndimd,cmat(iac),v)
( 2098) 
( 2099) !  contraction puis inversion des blocs diagonaux de la la matrice
( 2100) !  front par front .
( 2101) 
( 2102)       do 2 nb=2,nj
( 2103) !       write(*,*)'  front numero ',nb
( 2104) !  adresse du bloc diagonal precedent
( 2105)         iap=iac
( 2106)         iac=(nb-1)*lbd+1
( 2107) 
( 2108) !      lecture du bloc sur-diagonal precedent: passage bande-vecteur
( 2109) 
( 2110)         call feti_liblos(n,a,ndimd,nb-1,w)
( 2111) 
( 2112) !  elimination du bloc extra-diagonal .
( 2113) 
( 2114) !   lecture du bloc diagonal courant: passage bande-plein
( 2115) 
( 2116)         call feti_liblod(n,ni,a,ndimd,nb,cmat(iac))
( 2117) 
( 2118) !  calcul du nouveau bloc diagonal courant .
( 2119) 
( 2120) !       write(*,*)'  calcul du nouveau bloc diagonal courant '
( 2121)         call feti_nbdia(ndimd,w,cmat(iap),cmat(iac))
( 2122) 
( 2123) !  inversion du nouveau bloc diagonal courant .
( 2124) 
( 2125) !       write(*,*)'  factorisation du nouveau bloc diagonal courant '
( 2126)         j0=(nb-1)*ndimd
( 2127)         call feti_ldlt(ndimd,cmat(iac),d((nb-1)*ndimd+1),w,ndlblo,   &
( 2128)             lisblo,nmdlblo,j0)
( 2129)         call feti_calinv(ndimd,cmat(iac),v)
( 2130) 
( 2131)  2    continue
( 2132) 
( 2133)       return
( 2134)       end






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 38

( 2135) !**********************************************************************c
( 2136)       subroutine feti_calinv(dim,a,b)
( 2137) !  calcul de la partie triangulaire inferieure du bloc :
( 2138) !  b = [a]-1 puis recopie dans a complet .
( 2139)       implicit none
( 2140)       integer dim
( 2141)       REAL(kind=8)  a(dim,dim),b(dim,dim)
( 2142)       integer i,j
( 2143) 
( 2144)       call feti_vclr(dim*dim,b)
( 2145)       do j=1,dim
( 2146)         b(j,j)=1.e0
( 2147)         call feti_desremo(dim,j-1,a,b(1,j),b(1,j))
( 2148)       end do
( 2149)       do j=1,dim
( 2150)         do i=j,dim
( 2151)           a(i,j)=b(i,j)
( 2152)           a(j,i)=b(i,j)
( 2153)         end do
( 2154)       end do
( 2155) 
( 2156)       return
( 2157)       end
( 2158) !***********************************************************************
( 2159)       subroutine feti_nbdia(n,d,w,v)
( 2160)       implicit REAL(kind=8) (a-h,o-z)
( 2161)       dimension d(n),v(n,n),w(n,n)
( 2162) 
( 2163)       do 1 j=1,n
( 2164)         do 2 i=1,n
( 2165)           v(i,j)=v(i,j)-d(i)*w(i,j)*d(j)
( 2166)  2      continue
( 2167)  1    continue
( 2168) 
( 2169)       return
( 2170)       end
( 2171) !***********************************************************************
( 2172)       subroutine feti_liblod(n,ni,a,ndimd,nb,db)
( 2173) 
( 2174)       implicit REAL(kind=8) (a-h,o-z)
( 2175)       dimension db(ndimd,ndimd)
( 2176) !  tableaux descripteurs de la structure morse de la matrice .
( 2177)       dimension a(n,5)
( 2178) 
( 2179)       call feti_vclr(ndimd*ndimd,db)
( 2180) !       do i=1,ndimd*ndimd
( 2181) !         db(i,1)=0.e0
( 2182) !8     continue
( 2183)       n0=(nb-1)*ndimd
( 2184) !  lecture de la diagonale db(i,i) .
( 2185)       do 3 i=1,ndimd
( 2186)         db(i,i)=a(n0+i,3)
( 2187)  3    continue
( 2188) !  lecture de la premiere sur-diagonale db(i,i+1) .
( 2189)       do 4 i=1,ndimd-1
( 2190)         db(i,i+1)=a(n0+i,4)
( 2191)  4    continue
( 2192) !  lecture de la premiere sous-diagonale db(i,i-1) .






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 39

( 2193)       do 2 i=2,ndimd
( 2194)         db(i,i-1)=a(n0+i,2)
( 2195)  2    continue
( 2196) 
( 2197)       return
( 2198)       end
( 2199) !***********************************************************************
( 2200)       subroutine feti_liblos(n,a,ndimd,nb,v)
( 2201) 
( 2202)       implicit REAL(kind=8) (a-h,o-z)
( 2203)       dimension v(ndimd)
( 2204) !  tableaux descripteurs de la structure morse de la matrice .
( 2205)       dimension a(n,5)
( 2206) 
( 2207)       n0=(nb-1)*ndimd
( 2208) !  lecture de la deuxieme sur-diagonale v(i) .
( 2209)       do 5 i=1,ndimd
( 2210)         v(i)=a(n0+i,5)
( 2211)  5    continue
( 2212) 
( 2213)       return
( 2214)       end
( 2215) !***********************************************************************
( 2216)       subroutine feti_resloc(n,ni,nj,a,lcmat,cmat,y,x,z)
( 2217)       implicit REAL(kind=8) (a-h,o-z)
( 2218)       dimension a(n,5)
( 2219) !  routine de resolution utilisant la decomposition frontale par blocs
( 2220) !  de la matrice heptadiagonale a .
( 2221) !  cmat contient les termes des blocs diagonaux factorises .
( 2222)       dimension cmat(lcmat)
( 2223) !  d contient l'inverse de la diagonale de Crout .
( 2224)       dimension y(n),x(n),z(n)
( 2225) !
( 2226)       call feti_vmov(n,y,z)
( 2227)       ndimd=ni
( 2228)       lbd=ndimd*ndimd
( 2229) !  descente du systeme .
( 2230)       do nb=1,nj-1
( 2231) !  calcul du second membre condense : z(i+1)=y(i+1)- FD-1 z(i) .
( 2232)         iacmat=(nb-1)*lbd+1
( 2233)         n0=(nb-1)*ndimd
( 2234)         call feti_mxv(ndimd,ndimd,cmat(iacmat),z(n0+1),x(n0+1))
( 2235)         call feti_pbloi(n,ni,nj,a,x,z,nb)
( 2236)       end do
( 2237) !  inversion du dernier bloc diagonal .
( 2238)       iacmat=(nj-1)*lbd+1
( 2239)       n0=(nj-1)*ndimd
( 2240)       call feti_mxv(ndimd,ndimd,cmat(iacmat),z(n0+1),x(n0+1))
( 2241) !  remontee du systeme par bloc .
( 2242)       do nb=nj-1,1,-1
( 2243) !  calcul du second membre condense : z(i)=z(i)- E x(i+1) .
( 2244)         call feti_pblos(n,ni,nj,a,x,z,nb)
( 2245) !  calcul de la solution pour le bloc diagonal courant .
( 2246)         iacmat=(nb-1)*lbd+1
( 2247)         n0=(nb-1)*ndimd
( 2248)         call feti_mxv(ndimd,ndimd,cmat(iacmat),z(n0+1),x(n0+1))
( 2249)       end do
( 2250) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 40

( 2251)       end
( 2252) !***********************************************************************
( 2253)       subroutine feti_pblos(n,ni,nj,a,x,y,nb)
( 2254) !  calcul du produit matrice-vecteur morse : y = y - A x .
( 2255) !  produit par le bloc superieur numero nb .
( 2256)       implicit REAL(kind=8) (a-h,o-z)
( 2257)       dimension a(n,5)
( 2258)       dimension x(n),y(n)
( 2259) 
( 2260)       n0=(nb-1)*ni
( 2261) !  produit par la deuxieme sur-diagonale a(i,i+ni) .
( 2262)       do i=n0+1,n0+ni
( 2263)         y(i)=y(i)-x(i+ni)*a(i,5)
( 2264)       end do
( 2265) 
( 2266)       end subroutine feti_pblos
( 2267) !***********************************************************************
( 2268)       subroutine feti_pbloi(n,ni,nj,a,x,y,nb)
( 2269) !  calcul du produit matrice-vecteur morse : y = y - A x .
( 2270) !  produit par le bloc inferieur numero nb .
( 2271)       implicit REAL(kind=8) (a-h,o-z)
( 2272)       dimension a(n,5)
( 2273)       dimension x(n),y(n)
( 2274) 
( 2275)       n0=nb*ni
( 2276) !  produit par la deuxieme sous-diagonale a(i,i-ni) .
( 2277)       do i=n0+1,n0+ni
( 2278)         y(i)=y(i)-x(i-ni)*a(i,1)
( 2279)       end do
( 2280) 
( 2281)       end subroutine feti_pbloi
( 2282) !***********************************************************************
( 2283)       subroutine feti_proax(n,ni,nj,a,x,y)
( 2284) !  calcul du produit matrice-vecteur morse par la matrice heptadiagonale
( 2285) !  a stockee par bandes .
( 2286)       implicit REAL(kind=8) (a-h,o-z)
( 2287)       dimension a(n,5)
( 2288)       dimension x(n),y(n)
( 2289) 
( 2290) !  produit par la diagonale .
( 2291)       do i=1,n
( 2292)         y(i)=x(i)*a(i,3)
( 2293)       end do
( 2294) !  produit par la premiere sur-diagonale a(i,i+1) .
( 2295)       do i=1,n-1
( 2296)         y(i)=y(i)+x(i+1)*a(i,4)
( 2297)       end do
( 2298) !  produit par la deuxieme sur-diagonale a(i,i+ni) .
( 2299)       do i=1,n-ni
( 2300)         y(i)=y(i)+x(i+ni)*a(i,5)
( 2301)       end do
( 2302) ! produit par la premiere sous-diagonale a(i,i-1)
( 2303)       do i=2,n
( 2304)           y(i)=y(i)+x(i-1)*a(i,2)
( 2305)       end do
( 2306) ! produit par la deuxieme sous-diagonale a(i,i-ni)
( 2307)       do i=1+ni,n
( 2308)           y(i)=y(i)+x(i-ni)*a(i,1)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 41

( 2309)       end do
( 2310) 
( 2311)       end subroutine feti_proax
( 2312) 
( 2313) #else
( 2314)    !! no use of FETI librairy
( 2315) #endif
( 2316) 
( 2317)    !!======================================================================
( 2318) END MODULE lib_feti






















































PGF90 (Version     10.5)          09/23/2010  09:49:41      page 42

(    1) # 1 "lib_feti.F90"
(    1) MODULE lib_feti
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE lib_feti   ***
(    4)    !! Ocean solver : FETI library
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_feti.F90,v 1.2 2005/03/27 18:34:47 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 2314
( 2314)    !! no use of FETI librairy
( 2315) 
( 2316) # 2317
( 2317)    !!======================================================================
( 2318) END MODULE lib_feti
