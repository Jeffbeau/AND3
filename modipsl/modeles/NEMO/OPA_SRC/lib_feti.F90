MODULE lib_feti
   !!==============================================================================
   !!                       ***  MODULE lib_feti   ***
   !! Ocean solver : FETI library
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_feti.F90,v 1.2 2005/03/27 18:34:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_feti
   !!  feti_creadr
   !!  feti_prext
   !!  feti_init
   !!  feti_iclr
   !!  feti_vclr
   !!  feti_vmov
   !!  feti_vneg    (used in lib_feti)
   !!  feti_vsub
   !!  feti_vadd    (used in lib_feti)
   !!  feti_ldlt    (used in lib_feti)
   !!  feti_desremo (used in lib_feti)
   !!  feti_mxv     (used in lib_feti)
   !!  feti_mxvadd  (used in lib_feti)
   !!  feti_saut    (used in lib_feti)
   !!  feti_moyen   (used in lib_feti)
   !!  feti_assemb  (used in lib_feti)
   !!  feti_poids
   !!  feti_probit  (used in lib_feti)
   !!  feti_inisub
   !!  feti_subound
   !!  feti_subdir
   !!  feti_listdir
   !!  feti_blodir
   !!  feti_numblo
   !!  feti_blomat
   !!  feti_blomat1
   !!  feti_nullsp
   !!  feti_project (used in lib_feti)
   !!  feti_etesolv (used in lib_feti)
   !!  feti_osaxpy  (used in lib_feti) 
   !!  feti_proe    (used in lib_feti)
   !!  feti_proet   (used in lib_feti)
   !!  feti_consrhs (used in lib_feti)
   !!  feti_dualschur
   !!  feti_extend  (used in lib_feti)
   !!  feti_restri  (used in lib_feti)
   !!  feti_halfint
   !!  feti_front
   !!  feti_calinv  (used in lib_feti)
   !!  feti_nbdia   (used in lib_feti)
   !!  feti_liblod  (used in lib_feti)
   !!  feti_liblos  (used in lib_feti)
   !!  feti_resloc
   !!  feti_pblos   (used in lib_feti)
   !!  feti_pbloi   (used in lib_feti) 
   !!  feti_proax
   !!
   !!  EXTERNAL used (BLAS): sgemv, strsv, saxpy, irecv
   !!
   !! References :
   !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
   !!      A domain decomposition solver to compute the barotropic 
   !!      component of an OGCM in the parallel processing field.
   !!      Ocean Modelling, issue 105, december 94.
   !!
   !! Modifications :
   !!      original :  98-02 (M. Guyon)
   !!------------------------------------------------------------------------
   !! * Modules used
   USE lib_mpp                 ! distribued memory computing

CONTAINS

      subroutine feti_creadr(ial,ialmax,lmax,ltab,iatab,chtab)
!  ial = derniere adresse libre dans le super-tableau de longueur lmax .
!  ltab = taille de la structure de donnee tab , d'adresse iatab .
      implicit none
      external mynode

      integer ial,ialmax,lmax,ltab,iatab
      integer iialmax,iibidon,mynode
      character(len=12) chtab

      iatab=ial
      ial=iatab+ltab
      ialmax=max0(ial,ialmax)
      if(ial > lmax) then 
          write(0,*)'  creation de la structure de donnee ',chtab
          write(0,*)'  depassement taille super-tableau :'
          write(0,*)'  taille maximale , taille demandee :',lmax,ial
          stop
      endif

!... test de la memoire

      iialmax=-ialmax
      call mpp_min(iialmax,1,iibidon)
      iialmax=-iialmax
!     if(mynode() == 0)then
!         write(0,*)'  creation de la structure de donnee ',chtab
!         write(0,*)'taille max vaut : ',iialmax
!     endif

      end subroutine feti_creadr
!***********************************************************************
      subroutine feti_prext(n,x)
      implicit none
      integer n
      REAL(kind=8)  x(n)
      REAL(kind=8)   xmin,xmax
      integer i

      xmin=x(1)
      xmax=x(1)
      do i=2,n
        xmin=min(xmin,x(i))
        xmax=max(xmax,x(i))
      end do
!     write(0,100)xmin,xmax
 100  format(8d10.3)

      end subroutine feti_prext
!*********************************************************************
      subroutine feti_init(noeuds,x)
      implicit none
      integer noeuds
      REAL(kind=8)  x(noeuds)
      integer ji

      do ji=1,noeuds
        x(ji)=1.
      end do

      end subroutine feti_init
!***********************************************************************
      subroutine feti_iclr(n,x)
      implicit none
      integer n
      integer x(n)
      integer i

      do i=1,n
        x(i)=0
      end do

      end subroutine feti_iclr
!***********************************************************************
      subroutine feti_vclr(n,x)
      implicit none
      integer n
      REAL(kind=8)  x(n)
      integer i

      do i=1,n
        x(i)=0.e0
      end do

      return
      end
!***********************************************************************
      subroutine feti_vmov(n,x,y)
      implicit none
      integer n
      REAL(kind=8)  x(n),y(n)
      integer i
      do i=1,n
        y(i)=x(i)
      end do
      return
      end
!***********************************************************************
      subroutine feti_vneg(n,x,y)
      implicit none
      integer n
      REAL(kind=8)  x(n),y(n)
      integer i

      do i=1,n
        y(i)=-x(i)
      end do

      return
      end
!***********************************************************************
      subroutine feti_vsub(n,x,y,z)
      implicit none
      integer n
      REAL(kind=8)  x(n),y(n),z(n)
      integer i
      do i=1,n
        z(i)=x(i)-y(i)
      end do
      return
      end
!***********************************************************************
      subroutine feti_vadd(n,x,y,z)
      implicit none
      integer n
      REAL(kind=8)  x(n),y(n),z(n)
      integer i

      do i=1,n
        z(i)=y(i)+x(i)
      end do

      return
      end
!***********************************************************************
      subroutine feti_ldlt(n,a,d,v,ndlblo,lisblo,nmdlblo,j0)

!...ndlblo nombre de liberte bloques = nombre de pivots nuls = dim(Ker(Ep))

!...nmdlblo = Max (ndlblo(Ep)) <= superstructure management in f77 context
!            p=1,Np

      implicit none
      integer n,ndlblo,nmdlblo
      integer j0
      REAL(kind=8)  a(n,n),v(n),d(n)
      integer lisblo(nmdlblo)
      integer i,j,k
      REAL(kind=8)  piv,epsilo

!  pivot de reference .
      piv=a(1,1)
!  epsilon .
      epsilo=1.e-10
!  factorisation de Crout d'une matrice pleine .
      do j = 1 , n
!  calcul du produit de la ligne numero j par la diagonale .
        do k = 1 , j - 1
          v(k) = a(j,k) * a(k,k)
        end do
!  calcul des termes de la partie inferieure de la colonne .
        call sgemv('N',n-j+1,j-1,-1.e0,a(j,1),n,v,1,1.e0,a(j,j),1)
!  detection des pivots nuls .
        if(a(j,j) < (piv*epsilo)) then
            ndlblo=ndlblo+1
            lisblo(ndlblo)=j0+j
            do i=1,j-1
              a(i,j)=0.e0
            end do
            do i=j+1,n
              a(j,i)=0.e0
            end do
            a(j,j)=0.e0
            d(j)=0.e0
        else
            piv=a(j,j)
            d(j)=1.e0/a(j,j)
            do i = j+1 , n
              a(i,j) = a(i,j) * d(j)
            end do
        endif
      end do
!  recopie de la partie triangulaire superieure .
      do i = 1, n
        do j = i+1, n
          a(i,j)=a(j,i)
        end do
      end do
      do i=1,n
        a(i,i)=d(i)
      end do

      return
      end
!**********************************************************************c
      subroutine feti_desremo(n,j0,a,y,x) 
      implicit none
      integer n,i,j,j0
      REAL(kind=8)  a(n,n),x(n),y(n)

!  descente remontee du systeme .
!  initialisation du vecteur resultat .
      do i=1,n
        x(i)=y(i)
      end do
!  descente du systeme L x = x , L est a diagonale unitaire .
      call strsv('L','N','U',n-j0,a(j0+1,j0+1),n,x(j0+1),1)
!  division par la diagonale .
      do j=j0+1,n
        x(j)=x(j)*a(j,j)
      end do          
!  remontee du systeme U x = x , U est a diagonale unitaire . 
      call strsv('U','N','U',n-j0,a(j0+1,j0+1),n,x(j0+1),1)

      return
      end
!***********************************************************************
      subroutine feti_mxv(n,m,a,x,y) 
!  calcul du produit matrice-vecteur :  y = A x
      implicit none
      integer n,m
      REAL(kind=8)  a(n,m),x(m),y(n)
      REAL(kind=8)  alpha,beta
!  initialisation du vecteur y .
      alpha=1.e0
      beta=0.e0
      call sgemv('N',n,m,alpha,a,n,x,1,beta,y,1)
      return
      end
!**********************************************************************c
      subroutine feti_mxvadd(n,m,a,x,y) 
!  calcul du produit matrice-vecteur :  y = y + A x
      implicit none
      integer n,m
      REAL(kind=8)  a(n,m),x(m),y(n)
      REAL(kind=8)  alpha,beta
!  initialisation du vecteur y .

      alpha=1.e0
      beta=1.e0
      if((n > 0).and.(m > 0))then
          call sgemv('N',n,m,alpha,a,n,x,1,beta,y,1)
      endif

      return
      end
!**********************************************************************c
      subroutine feti_saut(numnes,listnes,    &
          plistin,numins,listin,              &
                numnos,v,w,bufin,bufout)
!!        narea,numnos,v,w,bufin,bufout)

!  computation of the jump of a field on the interfaces, in the case 
!  where subdomain number i is allocated to processor number i-1 .

!     INPUTS : numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               lint    =  list of inner interfaces
!               lext    =  list of extern interfaces
!               plistin =  pointer of sublists of interface nodes
!               numins  =  number of interface nodes
!               listin  =  list of interface nodes
!               narea     =  subdomain number
!               numnos  =  number of nodes in the subdomain
!               v       =  local field

!     OUTPUTS : w  =  jump of v

!     WORKSPACE : bufin   =  contributions of local field on interface
!                 bufout  =  contributions of outer fields on interface
!                 work
      implicit none 
!      external irecv
!!    integer irecv, work(100)
!!    integer numnes,numins,numnos,narea
      integer numnes,numins,numnos
      integer listnes(numnes,3), plistin(numnes+1), listin(numins)
      REAL(kind=8)  v(numnos),w(numins)
      REAL(kind=8)  bufin(numins),bufout(numins)
      integer i,j,l
      integer mesg

      mesg=128
!  receiving the values of the outer fields on the interface
!     do i=1,numnes
!       l=(plistin(i+1)-plistin(i))*8
!       work(i)=irecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
!     end do
!  gathering the values of the local field on the interface
      do i=1,numnes
        do j=plistin(i)+1,plistin(i+1)
          bufin(j)=v(listin(j))
        end do
      end do
!  sending the values of the local field on the interface
      do i=1,numnes
        l=(plistin(i+1)-plistin(i))*8
!       call mppsend(mesg+narea-1,bufin(plistin(i)+1),l,listnes(i)-1,0)
        call mppsend(mesg+listnes(i,3)-1,bufin(plistin(i)+1),l,   &
            listnes(i,1)-1,0)
      end do
!  computing the jump on each interface 
      do i=1,numnes
!       call msgwait(work(i))
        l=(plistin(i+1)-plistin(i))*8
!       call mpprecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
        call mpprecv(mesg+listnes(i,2)-1,bufout(plistin(i)+1),l)
      end do
      do j=1,numins
        w(j)=bufin(j)-bufout(j)
      end do

      return
      end
!**********************************************************************c
      subroutine feti_moyen(numnes,listnes,   &
          plistin,numins,listin,              &
          narea,numnos,weight,v,w,bufin,bufout)

!  averaging a field on the interfaces, in the case where subdomain
!  number i is allocated to processor number i-1 .

!     INPUTS : numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               lint    =  list of inner interfaces
!               lext    =  list of extern interfaces
!               plistin =  pointer of sublists of interface nodes
!               numins  =  number of interface nodes
!               listin  =  list of interface nodes
!               narea    =  subdomain number
!               numnos  =  number of nodes in the subdomain
!               weight  =  weighting vector
!               v       =  local field

!     OUTPUTS :      w  =  averaged field

!     WORKSPACE : bufin   =  contributions of local field on interface
!                 bufout  =  contributions of outer fields on interface
!                 work
      implicit none 
!     external irecv
!!    integer irecv, work(100)
      integer numnes,numins,numnos,narea
      integer listnes(numnes,3), plistin(numnes+1), listin(numins)
      REAL(kind=8)  weight(numnos),v(numnos),w(numnos)
      REAL(kind=8)  bufin(numins),bufout(numins)
      integer i,j,l,n
      integer mesg

      mesg=128*2
!  initialization
      call feti_vmov(numnos,v,w)
!  receiving the values of the outer fields on the interface
!     do i=1,numnes
!       l=(plistin(i+1)-plistin(i))*8
!       work(i)=irecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
!     end do
!  gathering the values of the local field on the interface
      do i=1,numnes
        do j=plistin(i)+1,plistin(i+1)
          bufin(j)=v(listin(j))
        end do
      end do
!  sending the values of the local field on the interface
      do i=1,numnes
        l=(plistin(i+1)-plistin(i))*8
!       call mppsend(mesg+narea-1,bufin(plistin(i)+1),l,listnes(i)-1,0)
        call mppsend(mesg+listnes(i,3)-1,bufin(plistin(i)+1),l,   &
            listnes(i,1)-1,0)
      end do
!  computing the jump on each interface 
      do i=1,numnes
!       call msgwait(work(i))
        l=(plistin(i+1)-plistin(i))*8
!       call mpprecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
        call mpprecv(mesg+listnes(i,2)-1,bufout(plistin(i)+1),l)
      end do
      do j=1,numins
        n=listin(j)
        w(n)=w(n)+bufout(j)
      end do
      do n=1,numnos
        w(n)=w(n)*weight(n)
      end do

      return
      end
!**********************************************************************c
      subroutine feti_assemb(numnes,listnes,   &
          plistin,numins,listin,               &
          narea,numnos,v,w,bufin,bufout)

!  assembling a field on the interfaces, in the case where subdomain
!  number i is allocated to processor number i-1 .

!     INPUTS :      numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               lint    =  list of inner interfaces
!               lext    =  list of extern interfaces
!               plistin =  pointer of sublists of interface nodes
!               numins  =  number of interface nodes
!               listin  =  list of interface nodes
!               narea     =  subdomain number
!               numnos  =  number of nodes in the subdomain
!               v       =  unassembled field

!     OUTPUTS :      w  =  assembled field

!     WORKSPACE : bufin   =  contributions of local field on interface
!                 bufout  =  contributions of outer fields on interface
!                 work
      implicit none 
!      external irecv
!!    integer irecv, work(100)
      integer numnes,numins,numnos,narea
      integer listnes(numnes,3), plistin(numnes+1), listin(numins)
      REAL(kind=8)  v(numnos),w(numnos)
      REAL(kind=8)  bufin(numins),bufout(numins)
      integer i,j,l,n
      integer mesg

      mesg=128*3
!  initialization
      call feti_vmov(numnos,v,w)
!  receiving the values of the outer fields on the interface
!     do i=1,numnes
!       l=(plistin(i+1)-plistin(i))*8
!       work(i)=irecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
!     end do
!  gathering the values of the local field on the interface
      do i=1,numnes
        do j=plistin(i)+1,plistin(i+1)
          bufin(j)=v(listin(j))
        end do
      end do
!  sending the values of the local field on the interface
      do i=1,numnes
        l=(plistin(i+1)-plistin(i))*8
!       write(90,*)'envoie a ',listnes(i),'de ',l,'octets'
!       write(90,*)'tag est : ',listnes(i,3)
!       call mppsend(mesg+narea-1,bufin(plistin(i)+1),l,listnes(i)-1,0)
        call mppsend(mesg+listnes(i,3)-1,bufin(plistin(i)+1),l,   &
            listnes(i,1)-1,0)
      end do
!  computing the jump on each interface 
      do i=1,numnes
!       call msgwait(work(i))
        l=(plistin(i+1)-plistin(i))*8
!       write(90,*)'recoie de ',listnes(i),'de ',l,'octets'
!       write(90,*)'recois avec tag :',listnes(i,2)
!       call mpprecv(mesg+listnes(i)-1,bufout(plistin(i)+1),l)
        call mpprecv(mesg+listnes(i,2)-1,bufout(plistin(i)+1),l)
      end do
      do j=1,numins
        n=listin(j)
        w(n)=w(n)+bufout(j)
      end do

      return
      end
!**********************************************************************c
!!    subroutine feti_poids(numnes,listnes,plistin,numins,listin,   &
!!        narea,numnos,w,bufin,bufout)
      subroutine feti_poids(numnes,                numins,listin,   &
                numnos,w             )

!  computation of the subdomain weighting vector, in the case where 
!  subdomain number i is allocated to processor number i-1 .

!     INPUTS :      numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               plistin =  pointer of sublists of interface nodes
!               numins  =  number of interface nodes
!               listin  =  list of interface nodes
!               narea     =  subdomain number
!               numnos  =  number of nodes in the subdomain

!     OUTPUTS :      w  =  weighting vector

!     WORKSPACE : bufin   =  contributions of local field on interface
!                 bufout  =  contributions of outer fields on interface
!                 work
      implicit none 
!!    integer numnes,numins,numnos,narea
      integer numnes,numins,numnos
!!    integer listnes(numnes), plistin(numnes+1), listin(numins)
      integer                                     listin(numins)
      REAL(kind=8)  w(numnos)
!!    REAL(kind=8)  bufin(numins),bufout(numins)
      integer j,n

!  initialization
      do n=1,numnos
        w(n)=1.e0
      end do
!  computing the weight on each interface 
      do j=1,numins
        n=listin(j)
        w(n)=w(n)+1.e0
      end do
      do n=1,numnos
        w(n)=1.e0/w(n)
      end do

      return
      end
!**********************************************************************c
      subroutine feti_probit(numnes,listnes,plistin,numins,listin,   &
          numnos,w,bitw)

!  computing the local right-hand-side associated with a field
!  on the interface .

!     INPUTS :      numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               plistin =  pointer of sublists of interface nodes
!               numins  =  number of interface nodes
!               listin  =  list of interface nodes
!               numnos  =  number of nodes in the subdomain
!               w       =  interface field

!     OUTPUTS :      bitw  =  local right-hand side

      implicit none 
      integer numnes,numins,numnos
      integer listnes(numnes,3), plistin(numnes+1), listin(numins)
      REAL(kind=8)  w(numins)
      REAL(kind=8)  bitw(numnos)
      integer i,j,n

!  initialization
      do n=1,numnos
        bitw(n)=0.e0
      end do
!  assembling the values on the interfaces 
      do i=1,numnes
        do j=plistin(i)+1,plistin(i+1)
          n=listin(j)
          bitw(n)=bitw(n)+w(j)
        end do
      end do

      return
      end
!**********************************************************************c
      subroutine feti_inisub(ni,nj,ibondi,ibondj,iperio,   &
          ibsw,ibnw,ibse,ibne,                             &
          ninterf,ninterfc,nni,nnic)
      implicit none
      integer ni,nj,ibondi,ibondj,iperio
      integer ibsw,ibnw,ibse,ibne
      integer ibsw2,ibnw2,ibse2,ibne2
      integer ninterf,ninterfc,nni,nnic
!     integer i,j

!  determination de la position du sous-domaine dans la grille
!  cela depnd de ibondi et ibondj

!  nombre de sous-domaines voisins et de la longueur du tableau des
!  noeuds voisins

!  initialisation

      ninterf=0
      nni=0

!...first, the "segment-points"

      if((ibondi == 2).AND.(iperio /= 1))then

!  la periodicite : pas d'interface east / west si nbondi = 2
!  ET nperio != 1

!      elseif(iperio == 1) then

!...ca doit marcher!?

      else

!  west

          if(ibondi /= -1)then
              ninterf=ninterf+1
              nni=nni+nj+1
          endif

!  east

          if(ibondi /= 1)then
              ninterf=ninterf+1
              nni=nni+nj+1
          endif

      endif

!  south

      if(ibondj /= -1.AND.ibondj /= 2)then
          ninterf=ninterf+1
          nni=nni+ni+1
      endif

!  north

      if(ibondj /= 1.AND.ibondj /= 2)then
          ninterf=ninterf+1
          nni=nni+ni+1
      endif

      nnic=nni
      ninterfc=ninterf


!...second, the "corner-points"

!...determination of the neighboorings, they depends on iperio
!   boundary condition => local description : ipXX becomes ipXX2...

      ibnw2 = ibnw
      ibne2 = ibne
      ibsw2 = ibsw
      ibse2 = ibse

!...iperio boundary condition effect

      if(iperio == 1) then

!...general case : Earth == infinite tube

          ibnw2 = 1
          ibne2 = 1
          ibsw2 = 1
          ibse2 = 1

!...real boundary condition

          if(ibondj == -1.OR.ibondj == 2) then
              ibsw2 = 0
              ibse2 = 0
          endif

          if(ibondj == 1.OR.ibondj == 2) then
              ibnw2 = 0
              ibne2 = 0
          endif

      endif

!  la periodicite : pas de coin si nbondi = 2
!  ET nperio != 1

!  north-west

      if(ibnw == 1)then
          ninterfc=ninterfc+1
          nnic=nnic+1
      endif

!  north-east

      if(ibne == 1)then
          ninterfc=ninterfc+1
          nnic=nnic+1
      endif

!  south-west

      if(ibsw == 1)then
          ninterfc=ninterfc+1
          nnic=nnic+1
      endif

!  south-east

      if(ibse == 1)then
          ninterfc=ninterfc+1
          nnic=nnic+1
      endif

!       write(0,*)'  nbre de voisins =',ninterf
!       write(0,*)'  longueur des interfaces =',nni
!       write(0,*)'  longueur des interfaces =',nnic

      return
      end
!**********************************************************************c
      subroutine feti_subound(ni,nj,ildi,ilei,ildj,ilej,   &
          imoi,ibondi,ibondj,iperio,   &
          ninterf,ninterfc,   &
          iowe,ioea,ioso,iono,   &
          ibsw,ibnw,ibse,ibne,   &
          ipsw,ipnw,ipse,ipne,   &
          nsdvois,nsdvoisc,   &
          plistin,nni,listin)

      implicit none
      integer ni,nj,ildi,ilei,ildj,ilej
      integer imoi,ibondi,ibondj,iperio
      integer iowe,ioea,ioso,iono
      integer ibsw,ibnw,ibse,ibne
      integer ipsw,ipnw,ipse,ipne
      integer ninterf,ninterfc,nni
      integer nsdvois(ninterf,3),nsdvoisc(ninterfc,3)
      integer plistin(ninterfc+1),listin(nni)
      integer     nint,nni0,ii,jj
!     integer i,j,nint,nni0,ii,jj

      external mynode
      integer  mynode

!  determination de la position du sous-domaine dans la grille
!  cela depnd de ibondi et ibondj

! liste des sous-domaines voisins et des noeuds interface

! be carefull!!! in the FETI algorithm, area number is considered
! instead of process number

!  initialisation

      nint=0
      nni0=0
      plistin(1)=0

!...first, the "segment-points"

      if((ibondi == 2).AND.(iperio /= 1))then

!  la periodicite : pas d'interface east / west si nbondi = 2
!  ET nperio != 1

      elseif(iperio == 1) then

!  west

          nint=nint+1
          nsdvoisc(nint,1)=imoi
          do jj=1,nj+1
            listin(nni0+jj)=(jj-1)*(ni+1)+1
          end do
          nni0=nni0+nj+1
          plistin(nint+1)=nni0
          nsdvoisc(nint,2)=1
          nsdvoisc(nint,3)=2

!  east

          nint=nint+1
          nsdvoisc(nint,1)=imoi
          do jj=1,nj+1
            listin(nni0+jj)=(jj-1)*(ni+1)+ilei
          end do
          nni0=nni0+nj+1
          plistin(nint+1)=nni0
          nsdvoisc(nint,2)=2
          nsdvoisc(nint,3)=1

      else

!  west

          if(ibondi /= -1)then
              nint=nint+1
              nsdvoisc(nint,1)=(iowe+1)
              do jj=1,nj+1
                listin(nni0+jj)=(jj-1)*(ni+1)+1
              end do
              nni0=nni0+nj+1
              plistin(nint+1)=nni0
              nsdvoisc(nint,2)=1
              nsdvoisc(nint,3)=2
          endif

!  east

          if(ibondi /= 1)then
              nint=nint+1
              nsdvoisc(nint,1)=(ioea+1)
              do jj=1,nj+1
                listin(nni0+jj)=(jj-1)*(ni+1)+ilei
              end do
              nni0=nni0+nj+1
              plistin(nint+1)=nni0
              nsdvoisc(nint,2)=2
              nsdvoisc(nint,3)=1
          endif

      endif

!  south

      if(ibondj /= -1.AND.ibondj /= 2)then
          nint=nint+1
          nsdvoisc(nint,1)=(ioso+1)
          do ii=1,ni+1
            listin(nni0+ii)=ii
          end do
          nni0=nni0+ni+1
          plistin(nint+1)=nni0
          nsdvoisc(nint,2)=3
          nsdvoisc(nint,3)=4
      endif

!  north

      if(ibondj /= 1.AND.ibondj /= 2)then
          nint=nint+1
          nsdvoisc(nint,1)=(iono+1)
          do ii=1,ni+1
            listin(nni0+ii)=(ilej-1)*(ni+1)+ii
          end do
          nni0=nni0+ni+1
          plistin(nint+1)=nni0
          nsdvoisc(nint,2)=4
          nsdvoisc(nint,3)=3
      endif

!...second, the "corner-points"

!...determination of the neighboorings, they depends on iperio 
!   boundary condition => local description : inimpp or inimpp2...

!  la periodicite : pas de coin si nbondi = 2
!  ET nperio != 1

!!!!!!!!!verifier la definition des voisins-coin

!  north-west

      if(ibnw == 1)then
          nint=nint+1
          nsdvoisc(nint,1)=(ipnw+1)
          listin(nni0+1)=(ilej-1)*(ni+1)+1
          nni0=nni0+1
          plistin(nint+1)=nni0
          nsdvoisc(nint,2)=5
          nsdvoisc(nint,3)=8
      endif

!  north-east

      if(ibne == 1)then
          nint=nint+1
          nsdvoisc(nint,1)=(ipne+1)
          listin(nni0+1)=(ilej-1)*(ni+1)+ilei
          nni0=nni0+1
          plistin(nint+1)=nni0
          nsdvoisc(nint,2)=6
          nsdvoisc(nint,3)=7
      endif

!  south-west

      if(ibsw == 1)then
          nint=nint+1
          nsdvoisc(nint,1)=(ipsw+1)
          listin(nni0+1)=1
          nni0=nni0+1
          plistin(nint+1)=nni0
          nsdvoisc(nint,2)=7
          nsdvoisc(nint,3)=6
      endif

!  south-east

      if(ibse == 1)then
          nint=nint+1
          nsdvoisc(nint,1)=(ipse+1)
          listin(nni0+1)=ilei
          nni0=nni0+1
          plistin(nint+1)=nni0
          nsdvoisc(nint,2)=8
          nsdvoisc(nint,3)=5
      endif

!...je compte les coins

      do ii=1,ninterf
        nsdvois(ii,1) = nsdvoisc(ii,1)
        nsdvois(ii,2) = nsdvoisc(ii,2)
        nsdvois(ii,3) = nsdvoisc(ii,3)
      enddo


!  impressions

!     write(*,*)'  sous-domaine :',imoi
!     write(*,*)'  nombre de sous-domaines voisins :',ninterfc
!      do ii=1,ninterfc
!      write(*,*)'  sous-domaine voisin :',nsdvoisc(ii,1)
!     write(*,*)'  nombre de noeuds interface :',   &
!         plistin(ii+1)-plistin(ii)
!     write(*,*)'  noeuds interface :',(listin(jj),jj=plistin(ii)+1,   &
!         plistin(ii+1))
!     end do

      return
      end
!**********************************************************************c
      subroutine feti_subdir(ni,nj,noeuds,ndir,mgcnum)
      implicit none 
      integer ni,nj,noeuds,ndir
      integer mgcnum(ni+1,nj+1)
      integer i,j
!     integer i,j,k

!  determination de la position du sous-domaine dans la grille
!  narea = (j-1)*iglo + i , avec 1<=i<=iglo et 1<=j<=jglo 



! 1. NUMBER THE GRID-POINTS USING BMASK
! -------------------------------------


      ndir=0
      do j=1,nj+1
        do i=1,ni+1
          ndir=ndir+1-mgcnum(i,j)
!      do k=1,noeuds
!        ndir=ndir+1-mgcnum(k,1)
        enddo
      enddo

      return
      end
!**********************************************************************c
      subroutine feti_listdir(ni,nj,logdir,ndir,lisdir)
      implicit none 
      integer ni,nj,ndir
      integer logdir(ni,nj),lisdir(ndir)
      integer i0,ji,jj

!  creation de la liste des degres de liberte bloques
      i0=0
      do ji=1,ni
        do jj=1,nj
          if(logdir(ji,jj) == 0) then
              i0=i0+1
              lisdir(i0)=ji+(jj-1)*ni
          endif
        end do
      end do
      if(i0 /= ndir) then
!      write(*,*)'  nombre de ddl bloques/prevus :',i0,ndir
!      write(*,*)'liste des ddl Dirichlet :',(lisdir(ji),ji=1,ndir)
          stop
      endif

      return
      end
!***********************************************************************
      subroutine feti_blodir(n,x,ndlblo,list)
      implicit none
      integer n,ndlblo
      integer list(ndlblo)
      REAL(kind=8)  x(n)
      integer i

!  remise a zero des ddl bloques
      do i=1,ndlblo
        x(list(i))=0.e0
      end do

      return
      end
!***********************************************************************
      subroutine feti_numblo(ndlblo,lisblo)
      implicit none
      external mynode
!!    integer mynode,indlblog,iibidon
      integer mynode
      integer ndlblo
      integer lisblo(ndlblo)
!!    integer n

!      indlblog=-ndlblo
!      call mpp_min(indlblog,1,iibidon)
!      indlblog=-indlblog
!      if(mynode() == 0)THEN
!      write(0,*)'  nombre de degres de liberte flottants :',ndlblo
!      write(0,*)'  liste :',(lisblo(n),n=1,ndlblo)
!      endif

      return
      end
!***********************************************************************
      subroutine feti_blomat(ni,nj,a,ndlblo,lisblo)
      implicit none
      integer ni,nj
      REAL(kind=8)  a(ni,nj,5)
      integer ndlblo
      integer lisblo(ndlblo)
!!    integer i,j,k,l,n
      integer i,j,k,n

      do n=1,ndlblo
!  degre de liberte bloque : lisblo(n) = (j-1)*ni + i
        i=mod(lisblo(n)-1,ni)+1
        j=((lisblo(n)-1)/ni)+1
!  annulation des coefficients sur la ligne
        do k=1,5
          a(i,j,k)=0.e0
        end do
        a(i,j,3)=1.e0
!  annulation des coefficients dans les colonnes correspondantes
        if(i > 1) then
            a(i-1,j,4)=0.e0
        endif
        if(i < ni) then
            a(i+1,j,2)=0.e0
        endif
        if(j > 1) then
            a(i,j-1,5)=0.e0
        endif
        if(j < nj) then
            a(i,j+1,1)=0.e0
        endif
      end do

      return
      end
!***********************************************************************
      subroutine feti_blomat1(ni,nj,a,ndlblo,lisblo,nsp)
      implicit none
      integer ni,nj,ndlblo
      REAL(kind=8)  a(ni,nj,5),nsp(ni,nj,ndlblo)
      integer lisblo(ndlblo)
!!    integer i,j,k,l,n
      integer i,j,n

      do n=1,ndlblo
!  degre de liberte bloque : lisblo(n) = (j-1)*ni + i
        i=mod(lisblo(n)-1,ni)+1
        j=((lisblo(n)-1)/ni)+1
!  annulation des coefficients dans les colonnes correspondantes
        if(i > 1) then
            nsp(i-1,j,n)=-a(i-1,j,4)
        endif
        if(i < ni) then
            nsp(i+1,j,n)=-a(i+1,j,2)
        endif
        if(j > 1) then
            nsp(i,j-1,n)=-a(i,j-1,5)
        endif
        if(j < nj) then
            nsp(i,j+1,n)=-a(i,j+1,1)
        endif
      end do

      return
      end
!***********************************************************************
      subroutine feti_nullsp(noeuds,ni,nj,lpblo,blo,a,ndlblo,lisblo,nsp,z)
!  calcul du noyau de la matrice
      implicit none
      integer noeuds,ni,nj,lpblo,ndlblo
      REAL(kind=8)  a(ni,nj,5)
      integer lisblo(ndlblo)
      REAL(kind=8)  blo(lpblo),nsp(noeuds,ndlblo),z(noeuds)
      external sdot
      REAL(kind=8)  sdot,res
      integer j,i

!  calcul de [Aii]-1 * Aib
      do j=1,ndlblo
        call feti_resloc(noeuds,ni,nj,a,lpblo,blo,nsp(1,j),nsp(1,j),z)
      end do
!  remise a identite du bloc diagonal correspondant aux modes bloques
      do i=1,ndlblo
        do j=1,ndlblo
          nsp(lisblo(i),j)=0.e0
        end do
        nsp(lisblo(i),i)=1.e0
      end do
!  verifications
      do j=1,ndlblo
        call feti_proax(noeuds,ni,nj,a,nsp(1,j),z)
        res=sqrt(sdot(noeuds,z,1,z,1)/   &
            sdot(noeuds,nsp(1,j),1,nsp(1,j),1))
!       write(0,*)'  residu de la colonne ',j,' du noyau :',res
!       write(0,*)'  vecteur colonne :'
!       call prvec(noeuds,nsp(1,j))
      end do

      return
      end
!**********************************************************************c
      subroutine feti_project(gint,pgint,nspdim,x,b,nmspdim,numit0   &
          ,nitmax,gvec,agvec,d,ad,add,gamm,numnes,listnes,plistin,numins   &
          ,listin,narea,numnos,nsp,v,w,bufin,bufout,work)

      implicit none 

      integer nitmax, nspdim
      REAL(kind=8)  x(nspdim), b(nspdim)

!...the Krylov optimisation used in PCG associated to the Coarse Solver
!   involve to keep di & adi vectors in etesolve procedure
!   the nmspdim parameter is introduced to avoid message "out of bound"
!   arrising qhen the memory is checked!

!   nmspdim  =   Max (nspdim(Ep)) <= memory management (f77 context)
!               p=1,Np

!      REAL(kind=8)  gvec(nspdim), agvec(nspdim),
!     &       d(nspdim,0:nitmax-1), ad(nspdim,0:nitmax-1), 
!     &       add(0:nitmax-1), gamm(0:nitmax-1)

      integer nmspdim
      REAL(kind=8)  gvec(nspdim), agvec(nspdim),   &
          d(nmspdim,0:nitmax-1), ad(nmspdim,0:nitmax-1),    &
          add(0:nitmax-1), gamm(0:nitmax-1)

      integer numnes,numins,narea,  numnos, numit0
      integer listnes(numnes,3), plistin(numnes+1), listin(numins)
      REAL(kind=8)  nsp(numnos,nspdim), v(numnos),    &
          w(numins), gint(numins), pgint(numins) 
      REAL(kind=8)  bufin(numins), bufout(numins)
      REAL(kind=8)  work(0:nitmax-1)


      call feti_proet(numnes,listnes,plistin,numins,listin,numnos,   &
          gint,nspdim,nsp,b,v)
      call feti_etesolv(nspdim,x,b,nmspdim,numit0,nitmax,gvec,agvec,d,ad   &
          ,add,gamm,numnes,listnes,plistin,numins,listin,narea,numnos   &
          ,nsp,v,w,bufin,bufout,work)
      call feti_proe(numnes,listnes,plistin,numins,listin,narea,numnos,   &
          nspdim,nsp,x,v,w,bufin,bufout)
      call feti_vsub(numins,gint,w,pgint)

      return
      end
!**********************************************************************c
      subroutine feti_etesolv(nspdim,x,b,nmspdim,numit0,nitmax,gvec    &
          ,agvec,d,ad,add,gamm,numnes,listnes,plistin,numins,listin    &
          ,narea,numnos,nsp,v,w,bufin,bufout,work)

!     INPUTS :   nspdim  =   dimension of null space : LOCAL
!                x       =   initial solution
!                b       =   right hand side
!               nmspdim  =   Max (nspdim(Ep)) <= memory management
!                            p=1,Np
!               numit0   =   number of previous direction vectors, at most
!                            nitmax
!               nitmax   =   number of maximum direction vectors for 
!                            reconjugation  
!               gvec     =   gradient vector
!               agvec    =   A * g
!               d        =   direction vector
!               ad       =   A * d
!               add      =   dot products (Ad,d)
!               gamm     =   GAMMA coefficients
!               numnes   =   number of neighbouring subdomains
!               listnes  =   list of neighbouring subdomains
!               plistin  =   pointer of sublists of interface nodes
!               numins   =   number of interface nodes
!               listin   =   list of interface nodes
!               narea    =   subdomain number
!               numnos   =   number of nodes in the subdomain

!     WORKSPACE : bufin  =  contributions of local fields on interface
!                 bufout =  contributions of outer fields on interface
!                 work   =  workspace for mpp_sum calls

      implicit none 

      integer nitmax, nspdim
      REAL(kind=8)  x(nspdim), b(nspdim)

!...the Krylov optimisation involve to keep di & adi vectors
!   the nmspdim parameter is introduce to avoid message "out of bound"
!   arrising qhen the memory is checked!

!      REAL(kind=8)  gvec(nspdim), agvec(nspdim),
!     &       d(nspdim,0:nitmax-1), ad(nspdim,0:nitmax-1), 
!     &       rho, temp1(2), temp2(2), add(0:nitmax-1),
!     &       gamm(0:nitmax-1), facg0, facgn, facst, eps

      integer nmspdim
      REAL(kind=8)  gvec(nmspdim), agvec(nmspdim),    &
          d(nspdim,0:nitmax-1), ad(nspdim,0:nitmax-1),     &
          rho, temp1(2), temp2(2), add(0:nitmax-1),    &
          gamm(0:nitmax-1), facg0, facgn, facst, eps

      integer numnes,numins,narea,  numnos
      integer listnes(numnes,3), plistin(numnes+1), listin(numins)
      REAL(kind=8)  nsp(numnos,nspdim), v(numnos), w(numins)
      REAL(kind=8)  bufin(numins), bufout(numins)
      REAL(kind=8)  work(0:nitmax-1)

      integer i, j, n, numit0, nn, n1, numit1

      external sdot
      REAL(kind=8)  sdot

!    initialisation 

      do 1 i=1,nspdim
        x(i) = 0.e0
        gvec(i) = -b(i)
 1    continue
      facg0 = sdot( nspdim, gvec, 1, gvec , 1)
      call mpp_sum(facg0,1,work)

      if( facg0  ==  0.e0 ) return

      eps=1.e-24
      facst=facg0*eps

!    initial projection 

      if(numit0 > 0) then
          do 91 j = 0 , numit0
            gamm(j) = - sdot(nspdim,gvec,1,d(1,j),1) / add(j)
 91       continue
          call mpp_sum(gamm,numit0+1,work)
          call feti_mxvadd( nspdim, numit0+1, d, gamm, x )
          call feti_mxvadd( nspdim, numit0+1, ad, gamm, gvec )

!  test residual after initial projection

          facgn = sdot( nspdim, gvec, 1, gvec, 1 )
          call mpp_sum(facgn,1,work)
          if( facgn <= facst ) then
!             write(0,*) 'residual after initial projection for feti_etesolv='   &
!                 ,facgn/facg0
              return
          endif

      endif

!     compute (Et E) * gvec

      call feti_proe(numnes,listnes,plistin,numins,listin,narea,numnos,   &
          nspdim,nsp,gvec,v,w,bufin,bufout)
      call feti_proet(numnes,listnes,plistin,numins,listin,numnos,w,   &
          nspdim,nsp,agvec,v)

!   conjugate d .

      if(numit0 > 0) then
          numit1=min0(numit0+1,nitmax-1)
          do 94 i=0,numit0
            gamm(i) = -sdot(nspdim,gvec,1,ad(1,i),1)/ add(i)
 94       continue
          call mpp_sum(gamm,numit0+1,work)

          call feti_osaxpy( nspdim,gamm(numit0),d(1,numit0),   &
              gvec,d(1,numit1))
          call feti_osaxpy( nspdim,gamm(numit0),ad(1,numit0),   &
              agvec,ad(1,numit1))

          call feti_mxvadd( nspdim, numit1-1, d, gamm, d(1,numit1) )
          call feti_mxvadd( nspdim, numit1-1, ad, gamm, ad(1,numit1) )
      else
          numit1=0

          do 20 i = 1, nspdim
            d(i,numit1) = gvec(i)
 20       continue

          do 30 i = 1, nspdim
            ad(i,numit1) = agvec(i)
 30       continue
      endif

!     computing rho

      do 50 nn = numit1, nitmax
        n=min0(nn,nitmax-1)
        temp1(1) = sdot( nspdim, gvec, 1, d(1,n), 1 )
        temp1(2) = sdot( nspdim, ad(1,n), 1, d(1,n), 1 )
        call mpp_sum(temp1,2,temp2)
        add(n) = temp1(2)
        rho = - temp1(1) / add(n)

        call feti_osaxpy( nspdim, rho, d(1,n), x, x )
        call feti_osaxpy( nspdim, rho, ad(1,n), gvec, gvec )

!        test residual . 

        facgn = sdot( nspdim, gvec, 1, gvec, 1 )
        call mpp_sum(facgn,1,work)
!        write(0,*) 'iteration =', nn+1, 'residual = ', facgn/facg0
        if( facgn <= facst ) then
            goto 999
        endif

!   reconjugation of gvec

        do 60 i=0,n
          gamm(i) = -sdot(nspdim,gvec,1,ad(1,i),1)/ add(i)
 60     continue
        call mpp_sum(gamm,n+1,work)
        n1=min0(nn+1,nitmax-1)

        call feti_osaxpy( nspdim,gamm(n),d(1,n),gvec,d(1,n1))

        call feti_mxvadd( nspdim, n1-1, d, gamm, d(1,n1) )

!     compute (Et E) * d

        call feti_proe(numnes,listnes,plistin,numins,listin,narea,   &
            numnos,nspdim,nsp,d(1,n1),v,w,bufin,bufout)
        call feti_proet(numnes,listnes,plistin,numins,listin,numnos,w,   &
            nspdim,nsp,ad(1,n1),v)

 50   continue

      if(narea == 1) then
          write(0,*) 'No convergence for ete in ',nitmax,' iterations'
          stop
      endif


 999  continue
      numit0=min0(n,nitmax-1)
!     write(0,*) 'number of iterations for ete :', nn+1, ' ,  residual =', facgn

      return
      end
!**********************************************************************c
      subroutine feti_osaxpy( n, a, x, y, z )

      implicit none

      integer n
      REAL(kind=8)  a, x(n), y(n), z(n)

      integer i

      do 10 i = 1, n
        z(i) = y(i) + a * x(i)
 10   continue

      return
      end
!**********************************************************************c
      subroutine feti_proe(numnes,listnes,plistin,numins,listin,narea,   &
          numnos,nspdim,nsp,alpha,v,w,bufin,   &
          bufout)

!  jump of a velocity field on the interfaces, in the case where 
!  subdomain number i is allocated to processor number i-1 .
!  the velocity field is Nsp * alpha, where Nsp is the nullspace .

!     INPUTS :      numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               plistin =  pointer of sublists of interface nodes
!               numins  =  number of interface nodes
!               listin  =  list of interface nodes
!               narea     =  subdomain number
!               numnos  =  number of nodes in the subdomain
!               nspdim  =  dimension of the null space
!               nsp     =  null space
!               alpha   =  components of the velocity field in the null
!                          space
!               v       =  corresponding velocity field

!     OUTPUTS :      w  =  jump on the interface

!     WORKSPACE : bufin   =  contributions of local field on interface
!                 bufout  =  contributions of outer fields on interface
!                 work
      implicit none 
      integer numnes,numins,numnos,narea,nspdim
      integer listnes(numnes,3), plistin(numnes+1), listin(numins)
      REAL(kind=8)  v(numnos),w(numins)
      REAL(kind=8)  bufin(numins),bufout(numins)
      REAL(kind=8)  alpha(nspdim), nsp(numnos,nspdim)
      integer i

!  computing the values of the local velocity field 
      do i=1,numnos
        v(i)=0.e0
      end do
      call feti_mxvadd( numnos, nspdim, nsp, alpha, v )
!  computing the jump on the interface
      call feti_saut(numnes,listnes,plistin,numins,listin,   &
                numnos,v,w,bufin,bufout)
!!        narea,numnos,v,w,bufin,bufout)

      return
      end
!**********************************************************************c
      subroutine feti_proet(numnes,listnes,plistin,numins,listin,   &
          numnos,w,nspdim,nsp,alpha,v)

!  projection in the null space of the right-hand-side associated with
!  a Lagrange multiplier w .

!     INPUTS :      numins  =  number of interface nodes
!               listin  =  list of interface nodes
!               numnos  =  number of nodes in the subdomain
!               w       =  Lagrange multiplier
!               nspdim  =  dimension of the null space
!               nsp     =  null space

!     OUTPUTS :      v     =  local right-hand-side
!               alpha =  projection of v in the null space

      implicit none 
      external sdot
      REAL(kind=8)  sdot
      integer numins,numnos,nspdim,numnes
      integer listnes(numnes,3),plistin(numnes+1),listin(numins)
      REAL(kind=8)  v(numnos),w(numins)
      REAL(kind=8)  alpha(nspdim), nsp(numnos,nspdim)
      integer i

!  assembling the values of the Lagrange multiplier on the interface
      call feti_probit(numnes,listnes,plistin,numins,listin,numnos,w,v)
!  computing the projection of v in the null space .
      do 1 i=1,nspdim
        alpha(i)=sdot(numnos,v,1,nsp(1,i),1)
 1    continue

      return
      end
!**********************************************************************c
      subroutine feti_consrhs( vdim,nspdim,rhs,nsp,cb )


      implicit none

      external sdot
      REAL(kind=8)  sdot
      integer vdim, nspdim

      REAL(kind=8)   rhs(vdim), nsp(vdim,nspdim), cb(nspdim)


      integer i

!  compute the right hand side for constraint 

      do i=1,nspdim
        cb(i)=-sdot(vdim,rhs,1,nsp(1,i),1)
      end do

      return
      end
!**********************************************************************c
      subroutine feti_dualschur(noeuds,ni,nj,a,lpblo,blo,   &
          ninterf,ninterfc,nni,nnic,   &
          ndvois,ndvoisc,plistin,listin,   &
          poids,u,v,f,bitw,utilu,   &
          lambda,g,pg,mg,nitmax,nmaxd,j0,wj,dwj,dwwj,   &
          gamm,work,bufin,bufout,narea,epsilo,ndlblo,   &
          lisblo,ndkerep,   &
          xnul,ynul,numit0ete,nitmaxete,eteg,   &
          eteag,eted,etead,eteadd,etegamm,nsp,etev,   &
          etew,nnih,plistih,gh,w,dw,   &
          residu,indic,jn)
      implicit none
!  nombre de ddl par noeuds, nombre de noeuds interne et interface
      integer noeuds,ni,nj,nni,nnic
!      integer jpj,jpi
!  numero du sous-domaine = numero de processeur + 1 
      integer narea
!  tableaux descripteurs de la matrice locale
      REAL(kind=8)  a(ni,nj,5)
!  tableaux descripteurs de l'inverse de la matrice du probleme local
      integer lpblo
      REAL(kind=8)  blo(lpblo)
!  tableaux descripteurs de l'interface
      integer ninterf,ninterfc
      integer ndvois(ninterf,3),ndvoisc(ninterfc,3)
      integer plistin(ninterfc+1),listin(nnic)
!  tableaux descripteurs de l'interface ajoutes pour nperio == 1
      REAL(kind=8)  poids(noeuds)
!  utilitaires pour les assemblages aux interfaces
      REAL(kind=8)  bufin(nnic),bufout(nnic)
!  vecteurs locaux
      REAL(kind=8)  u(noeuds),v(noeuds),f(noeuds)
      REAL(kind=8)  bitw(noeuds),utilu(noeuds)
!  vecteurs inconnues aux interfaces
!  mg est le gradient preconditionne : M = sum(Bi Ai Bit)
      REAL(kind=8)  lambda(nni),g(nni),pg(nni),mg(nni)
!  tableaux servant a gerer le stockage de la moitie des coefficients
!  des tableaux des directions de descente sur chaque interface
      integer nnih
      integer plistih(ninterf+1)
      REAL(kind=8)  gh(nnih)
!  tableaux des directions de descente
!  on appelle D la matrice condensee a l'interface : sum(Bi [Ai]-1 Bit)
!  j0 est le nombre de vecteurs D orthogonaux deja stockes
!  j0 est egal au plus a nmaxd-1
!  lorsque le nomtre total d'iterations atteint nmaxd, seule la derniere
!  direction calculee est stockee dans la colonne numero nmaxd
      integer nitmax,nmaxd,j0
      REAL(kind=8)  wj(nnih,nmaxd),dwj(nnih,nmaxd),dwwj(nmaxd)
      REAL(kind=8)  w(nni),dw(nni)
!  utilitaires pour les calculs de produits scalaires globaux
      REAL(kind=8)  work(nmaxd),gamm(nmaxd)
!  tableaux descripteurs du noyau de la matrice locale et tableaux 
!  servant a la projection sur le noyau

!...ndlblo is the dimension of the local nullspace .=<. the size of the
!   memory of the superstructure associated to the nullspace : ndkerep
!   indeed ndkerep = Max ndlblo = Max dim(Ker(Ep))
!                                p=1,Np
!   ndkerep is introduced to avoid messages "out of bounds" when memory
!   is checked

      integer ndkerep
      integer ndlblo,numit0ete,nitmaxete
      integer lisblo(ndlblo)
      REAL(kind=8)  xnul(ndlblo),ynul(ndlblo),eteg(ndlblo),eteag(ndlblo)

!...the Krylov optimisation used in PCG associated to the Coarse Solver
!   involve to keep di & adi vectors in etesolve procedure
!   the ndkerep parameter is introduced to avoid message "out of bound"
!   arrising when the memory is checked!

!   ndkerep  =   Max (ndlblo(Ep)) <= memory management (f77 context)
!               p=1,Np

!      REAL(kind=8)  eted(ndlblo,nitmaxete),etead(ndlblo,nitmaxete)

      REAL(kind=8)  eted(ndkerep,nitmaxete),etead(ndkerep,nitmaxete)
      REAL(kind=8)  eteadd(nitmaxete),etegamm(nitmaxete)
      REAL(kind=8)  etev(noeuds),etew(nni),nsp(noeuds,ndlblo)
!  indicateur de convergence

      integer indic
      REAL(kind=8)  residu
      integer jn
!  utilitaires
!!    integer j,k,j1,ji,jj,jk
      integer j,k,j1
      REAL(kind=8)  scal(2),ff,roj,epsilo,epsilo2
!  fonctions externes
      external sdot
      REAL(kind=8)  sdot
!  booleen pour caracterisation
!  numerique si llconv=1
      integer llconv

      llconv=0

!...critere d'arret epsilo2 = epsilo **2 => evite sqrt

      epsilo2 = epsilo * epsilo

!  calcul du carre scalaire du second membre global

      call feti_assemb(ninterfc,ndvoisc,plistin,nnic,listin,narea,noeuds   &
          ,f,utilu,bufin,bufout)
      ff=sdot(noeuds,utilu,1,utilu,1)
      call mpp_sum(ff,1,work)
!     write(90,*)'valeur du carre du second membre ds le Schur Dual :',ff

!  calcul du second membre contraint 

      call feti_consrhs(noeuds,ndlblo,f,nsp,ynul)

!  calcul de lambda0 satisfaisant la contrainte
 
!     call feti_vclr(nni,lambda)
      call feti_etesolv(ndlblo,xnul,ynul,ndkerep,numit0ete,nitmaxete   &
          ,eteg,eteag,eted,etead,eteadd,etegamm,ninterf,ndvois,plistin   &
          ,nni,listin,narea,noeuds,nsp,etev,etew,bufin,bufout,work)

      call feti_proe(ninterf,ndvois,plistin,nni,listin,narea,noeuds,   &
          ndlblo,nsp,xnul,etev,lambda,bufin,bufout)

!  calcul du nouveau second membre

      call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds,lambda,bitw)

      call feti_vadd(noeuds,bitw,f,bitw)


!  calcul du champ local initial u0
!     call feti_blodir(noeuds,bitw,ndlblo,lisblo)
      call feti_resloc(noeuds,ni,nj,a,lpblo,blo,bitw,u,utilu)

!  calcul du gradient initial g0 = feti_saut(u0)
!!    call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds,u,g   &
      call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds,u,g   &
          ,bufin,bufout)

!  reinitialisation de lambdaj0 par reconjugaison

      call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea,g, gh)
      do k=1,j0
        gamm(k)=sdot(nnih,gh,1,wj(1,k),1)
      end do
      if(j0 >= 1) call mpp_sum(gamm,j0,work)
      do k=1,j0
        gamm(k)=-gamm(k)/dwwj(k)
      end do
      call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, lambda,gh)
      call feti_mxvadd(nnih,j0,wj,gamm,gh)
      call feti_extend(ninterf,ndvois,plistin,nni,plistih,nnih,narea,gh, lambda)
!  calcul du second membre associe a lambdaj0
      if(j0 > 0) then
          call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds, lambda,bitw)
          call feti_vadd(noeuds,bitw,f,bitw)
!  calcul du champ local associe uj0
!         call feti_blodir(noeuds,bitw,ndlblo,lisblo)
          call feti_resloc(noeuds,ni,nj,a,lpblo,blo,bitw,u,utilu)
!  calcul du gradient initial gj0= feti_saut(uj0)
!!        call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds   &
          call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds   &
              ,u,g,bufin,bufout)
      endif
!  calcul du gradient projete 
!     call feti_vmov(nni,g,pg)
      call feti_project(g,pg,ndlblo,xnul,ynul,ndkerep,numit0ete   &
          ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf,ndvois   &
          ,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin,bufout,work)
!  calcul du gradient preconditionne 
!     call feti_vmov(nni,pg,mg)

      call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds, pg,bitw)

      call feti_proax(noeuds,ni,nj,a,bitw,v) 
!!    call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds,v, mg,bufin,bufout)
      call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds,v, mg,bufin,bufout)
!  verification du residu initial 
      residu=sdot(nni,mg,1,mg,1)
      call mpp_sum(residu,1,work)
      residu = (4. * residu) / ff
      if(residu <= epsilo2)then
!         residu=10.*sqrt(residu/ff)
!         if(residu < epsilo)then
          if(llconv == 1)then
              if(narea == 1)then
                  write(0,*)'  residu carre approche initial :',residu
              endif
          endif
!  calcul du champ global moyenne aux interfaces et du residu global 
          call feti_vsub(noeuds,u,etev,utilu)
          call feti_moyen(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
              noeuds,poids,utilu,utilu,bufin,bufout)
          call feti_vmov(noeuds,utilu,u)
          call feti_proax(noeuds,ni,nj,a,utilu,v) 
          call feti_vsub(noeuds,v,f,v)
          call feti_assemb(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
              noeuds,v,utilu,bufin,bufout)
          residu=sdot(noeuds,utilu,1,utilu,1)
          call mpp_sum(residu,1,work)
          residu=sqrt(residu/ff)
!         residu=residu/ff
          if(llconv == 1)then
              if(narea == 1)then
                  write(0,*)'  residu exact initial :',residu
              endif
          endif

 101      format(10d9.2)

          return

      endif


!  calcul de la premiere direction de descente
!  calcul du gradient preconditionne projete 
      call feti_project(mg,mg,ndlblo,xnul,ynul,ndkerep,numit0ete   &
          ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf,ndvois   &
          ,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin,bufout ,work)
      call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea,   &
          mg,wj(1,j0+1))
      do k=1,j0
        gamm(k)=sdot(nnih,wj(1,j0+1),1,dwj(1,k),1)
      end do
      if(j0 >= 1) call mpp_sum(gamm,j0,work)
      do k=1,j0
        gamm(k)=-gamm(k)/dwwj(k)
      end do
      call feti_mxvadd(nnih,j0,wj,gamm,wj(1,j0+1))

!  iterations de gradient conjugue a l'interface


!debug
      do j=j0+1,j0+nitmax

        if(llconv == 1.and.narea == 1) write(0,*)'etape numero',j-j0
!  j1 est le numero de colonne de la derniere direction de descente
!  j1 = j, numero d'iteration total, tant que j < nmaxd
!  au-dela, j1 = nmaxd
        j1=min0(j,nmaxd)
!  calcul de Dwi , D = sum ( Bi [Ai]-1 tBi )
        call feti_extend(ninterf,ndvois,plistin,nni,plistih,nnih,narea, wj(1,j1),w)
        call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds, w,bitw)
!       call feti_blodir(noeuds,bitw,ndlblo,lisblo)
        call feti_resloc(noeuds,ni,nj,a,lpblo,blo,bitw,v,utilu)
!!      call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds,v   &
        call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds,v   &
            ,dw,bufin,bufout)
        call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, dw,dwj(1,j1))
!  calcul du coefficient de descente
        call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, g,gh)
        scal(1)=sdot(nnih,gh,1,wj(1,j1),1)
        scal(2)=sdot(nnih,dwj(1,j1),1,wj(1,j1),1)
        call mpp_sum(scal,2,work)
        roj=-scal(1)/scal(2)
        dwwj(j1)=scal(2)
!  remise a jour de lambda, g et u
        call saxpy(nni,roj,w,1,lambda,1)
        call saxpy(nni,roj,dw,1,g,1)
        call saxpy(noeuds,roj,v,1,u,1)
!  calcul du gradient projete 
!       call feti_vmov(nni,g,pg)
        call feti_project(g,pg,ndlblo,xnul,ynul,ndkerep,numit0ete   &
            ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf   &
            ,ndvois,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin   &
            ,bufout,work)
!  calcul du gradient preconditionne Mg , M = sum ( Bi Ai tBi )
!       call feti_vmov(nni,pg,mg)
        call feti_probit(ninterf,ndvois,plistin,nni,listin,noeuds, pg,bitw)
        call feti_proax(noeuds,ni,nj,a,bitw,v) 
!!      call feti_saut(ninterf,ndvois,plistin,nni,listin,narea,noeuds,v   &
        call feti_saut(ninterf,ndvois,plistin,nni,listin,      noeuds,v   &
            ,mg,bufin,bufout)
!  verification du residu global approche
        residu=sdot(nni,mg,1,mg,1)
        call mpp_sum(residu,1,work)
!       residu=10.*sqrt(residu/ff)
        residu = (4. * residu) / ff
        if(llconv == 1.and.narea == 1)then
            write(0,*)'  residu carre global approche apres ',   &
                j-j0,' iterations :',residu
        endif

!       if(residu <  epsilo) then
        if(residu <= epsilo2)then
            if(llconv == 1.and.narea == 1)then
                write(0,*)'residu carre global approche apres ',   &
                    j-j0,' iterations :',residu
            endif
!  calcul du champ global moyenne aux interfaces et du residu global 
            call feti_vsub(noeuds,u,etev,utilu)
            call feti_moyen(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
                noeuds,poids,utilu,utilu,bufin,bufout)
            call feti_vmov(noeuds,utilu,u)
            call feti_proax(noeuds,ni,nj,a,utilu,v) 
            call feti_vsub(noeuds,v,f,v)
            call feti_assemb(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
                noeuds,v,utilu,bufin,bufout)
            residu=sdot(noeuds,utilu,1,utilu,1)
            call mpp_sum(residu,1,work)
!           residu=residu/ff
            residu=sqrt(residu/ff)
            if(llconv == 1.and.narea == 1)then
                write(0,*)'  residu global exact apres ',   &
                    j-j0,' iterations :',residu
            endif

            jn = j - j0
            j0=min0(j,nmaxd-1)
            if(llconv == 1.and.narea == 1)then
                write(0,*)'  nombre de directions orthogonales conservees :',j0
            endif


            return
        endif
!##################################################################
!  calcul du champ global moyenne aux interfaces et du residu global 
        call feti_vsub(noeuds,u,etev,utilu)
        call feti_moyen(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
            noeuds,poids,utilu,utilu,bufin,bufout)
        call feti_proax(noeuds,ni,nj,a,utilu,v) 
        call feti_vsub(noeuds,v,f,v)
        call feti_assemb(ninterfc,ndvoisc,plistin,nnic,listin,narea,   &
            noeuds,v,utilu,bufin,bufout)
        residu=sdot(noeuds,utilu,1,utilu,1)
        call mpp_sum(residu,1,work)
!       residu=sqrt(residu/ff)
        residu=residu/ff
        if(llconv == 1.and.narea == 1)then
            write(0,*)'  residu carre global exact apres ',   &
                j-j0,' iterations :',residu
        endif

!#####################################################################
!  calcul de la nouvelle direction de descente par reconjugaison
!  calcul du gradient preconditionne projete 
        call feti_project(mg,mg,ndlblo,xnul,ynul,ndkerep,numit0ete   &
            ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf   &
            ,ndvois,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin   &
            ,bufout,work)
        call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, mg,gh)
        do k=1,j1
          gamm(k)=sdot(nnih,gh,1,dwj(1,k),1)
        end do
        call mpp_sum(gamm,j1,work)
        do k=1,j1
          gamm(k)=-gamm(k)/dwwj(k)
        end do
        call feti_mxvadd(nnih,j1,wj,gamm,gh)
        j1=min0(j1+1,nmaxd)

!%%%%%%%%%%stabilisation numerique %%%%%%%%%%%%%%%%%%%%%%%%%
        call feti_extend(ninterf,ndvois,plistin,nni,plistih,nnih,narea, gh,w)
        call feti_project(w,w,ndlblo,xnul,ynul,ndkerep,numit0ete   &
            ,nitmaxete,eteg,eteag,eted,etead,eteadd,etegamm,ninterf   &
            ,ndvois,plistin,nni,listin,narea,noeuds,nsp,etev,etew,bufin   &
            ,bufout,work)
        call feti_restri(ninterf,ndvois,plistin,nni,plistih,nnih,narea, w,gh)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call feti_vmov(nnih,gh,wj(1,j1))

      end do

      residu = sqrt(residu)
      if(narea == 1) then
          write(0,*)'  residu global approche apres ',nitmax,   &
              '  iterations :',residu
      endif
      indic=-2

      return
      end
!**********************************************************************c
      subroutine feti_extend(numnes,listnes,   &
          plistin,numins,plistih,   &
          numinh,narea,wh,w)

!  reassembling an half interface field 

!     INPUTS :      numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               lint    =  list of inner interfaces
!               lext    =  list of extern ibnterfaces
!               plistin =  pointer of sublists of interface nodes
!               plistih =  pointer of half-sublists of interface nodes
!               numins  =  number of interface nodes
!               numinh  =  half number of interface nodes
!               narea    =  subdomain number
!               wh      =  half interface vector 

!     OUTPUTS :      w  =  complete interface vector

!     WORKSPACE : work
      implicit none 
!      external irecv
!!    integer irecv, work(100)
      integer numnes,numins,numinh,narea
      integer listnes(numnes,3), plistin(numnes+1), plistih(numnes+1)
      REAL(kind=8)  wh(numinh),w(numins)
      integer i,l,lin,lout
      integer mesg

      mesg=128*10
!  if subdomains i < j are neighbours, subdomain i manage the first
!  half of the interface, and subdomain j the remaining part

!  receiving the values of the outer fields on the interface
!  sending the values of the inner field on the interface

      do i=1,numnes
        l=(plistin(i+1)-plistin(i))
!        if(narea < listnes(i)) then
        if(listnes(i,3) < listnes(i,2))then
            lout=(l-(l/2))
            lin=(l/2)
!           work(i)=irecv(mesg+listnes(i)-1,w(plistin(i)+(l/2)+1), lout*8)
!           call mppsend(mesg+narea-1,wh(plistih(i)+1),lin*8,listnes(i) -1,0)
            call mppsend(mesg+listnes(i,3)-1,wh(plistih(i)+1),lin*8, listnes(i,1)-1,0)
            call feti_vmov(lin,wh(plistih(i)+1),w(plistin(i)+1))
        else
            lout=(l/2)
            lin=(l-(l/2))
!           work(i)=irecv(mesg+listnes(i)-1,w(plistin(i)+1),lout*8)
!           call mppsend(mesg+narea-1,wh(plistih(i)+1),lin*8,listnes(i)-1,0)
            call mppsend(mesg+listnes(i,3)-1,wh(plistih(i)+1),lin*8, listnes(i,1)-1,0)
            call feti_vmov(lin,wh(plistih(i)+1),w(plistin(i)+(l/2)+1))
        endif
      enddo
!  waiting for the completion on each interface 
      do i=1,numnes
        l=(plistin(i+1)-plistin(i))
!       if(narea < listnes(i)) then
        if(listnes(i,3) < listnes(i,2))then
            lout=(l-(l/2))
!           call mpprecv(mesg+listnes(i)-1,w(plistin(i)+(l/2)+1), lout*8)
            call mpprecv(mesg+listnes(i,2)-1,w(plistin(i)+(l/2)+1), lout*8)
            call feti_vneg(lout,w(plistin(i)+(l/2)+1),   &
                w(plistin(i)+(l/2)+1))
        else
            lout=(l/2)
!           call mpprecv(mesg+listnes(i)-1,w(plistin(i)+1),lout*8)
            call mpprecv(mesg+listnes(i,2)-1,w(plistin(i)+1),lout*8)
            call feti_vneg(lout,w(plistin(i)+1),w(plistin(i)+1))
        endif
      end do

      return
      end
!**********************************************************************c
      subroutine feti_restri(numnes,listnes,plistin,numins,plistih,   &
          numinh,narea,w,wh)

!  restriction of an interface field to one half on each interface

!     INPUTS :      numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               plistin =  pointer of sublists of interface nodes
!               plistih =  pointer of half-sublists of interface nodes
!               numins  =  number of interface nodes
!               numinh  =  half number of interface nodes
!               narea    =  subdomain number
!               w       =  complete interface vector 

!     OUTPUTS :      wh  =  half interface vector

!     WORKSPACE : work
      implicit none 
      integer numnes,numins,numinh,narea
      integer listnes(numnes,3), plistin(numnes+1), plistih(numnes+1)
      REAL(kind=8)  wh(numinh),w(numins)
      integer i,l,lin

!  if subdomains i < j are neighbours, subdomain i manage the first
!  half of the interface, and subdomain j the remaining part

      do i=1,numnes
        l=(plistin(i+1)-plistin(i))
!       if(narea < listnes(i)) then
        if(listnes(i,3) < listnes(i,2))then
            lin=(l/2)
            call feti_vmov(lin,w(plistin(i)+1),wh(plistih(i)+1))
        else
            lin=(l-(l/2))
            call feti_vmov(lin,w(plistin(i)+(l/2)+1),wh(plistih(i)+1))
        endif
      end do

      return
      end
!**********************************************************************c
!!    subroutine feti_halfint(numnes,listnes,plistin,numins,plistih   &
!!        ,numinh,narea)
      subroutine feti_halfint(numnes,listnes,plistin,       plistih   &
          ,numinh      )

!  construction of the pointer of the restriction of an interface 
!  field to one half on each interface

!     INPUTS :      numnes  =  number of neighbouring subdomains
!               listnes =  list of neighbouring subdomains
!               plistin =  pointer of sublists of interface nodes
!               numins  =  number of interface nodes
!               narea    =  subdomain number

!     OUTPUTS :      numinh  =  half number of interface nodes
!               plistih =  pointer of half-sublists of interface nodes

!     WORKSPACE : work
      implicit none 
!!    integer numnes,numins,numinh,narea
      integer numnes,       numinh
      integer listnes(numnes,3), plistin(numnes+1), plistih(numnes+1)
      integer i,l

!  if subdomains i < j are neighbours, subdomain i manage the first
!  half of the interface, and subdomain j the remaining part

      plistih(1)=0
      numinh=0
      do i=1,numnes
        l=(plistin(i+1)-plistin(i))
!        if(narea < listnes(i,1)) then
        if(listnes(i,3) < listnes(i,2)) then
!          write(*,*)'je gere interface avec : ',listnes(i,1)
            plistih(i+1)=plistih(i)+(l/2)
        else
!          write(*,*)'il gere interface avec : ',listnes(i,1)
            plistih(i+1)=plistih(i)+(l-(l/2))
        endif
      end do
      numinh=plistih(numnes+1)

      return
      end
!***********************************************************************
      subroutine feti_front(n,ni,nj,a,lcmat,cmat,d,v,w,ndlblo,lisblo   &
          ,nmdlblo)
      implicit REAL(kind=8) (a-h,o-z)
      dimension a(n,5)
!  routine de calcul de la decomposition frontale par blocs de
!  la matrice pentadiagonale a .
!  cmat contient les termes des blocs diagonaux factorises .
      dimension cmat(lcmat),d(n)

!...ndlblo nombre de liberte bloques = nombre de pivots nuls = dim(Ker(Ep))

!...nmdlblo = Max (ndlblo(Ep)) <= superstructure management in f77 context
!            p=1,Np

      integer ndlblo,nmdlblo
      integer lisblo(nmdlblo)

!  d ne contient pas, en sortie, l'inverse de la diagonale de Crout .
!  v est un tableau servant a stocker la partie triangulaire inferieure
!  stricte d'un bloc diagonal , et w un tableau servant a stocker un
!  bloc sur-diagonal plein .
!  de dimension egale a la largeur maximale d'un front .
      dimension v(ni,ni),w(ni)
!!!
      ndimd=ni
      lbd=ndimd*ndimd

!      lecture du tout premier bloc diagonal

      nb=1
!  adresse du bloc diagonal courant
      iac=1
      call feti_liblod(n,ni,a,ndimd,nb,cmat(iac))
!  inversion du premier bloc diagonal .
!     write(*,*)'  inversion du premier bloc diagonal '
      j0=0
      call feti_ldlt(ndimd,cmat(iac),d(1),w,ndlblo,lisblo,nmdlblo,j0)
      call feti_calinv(ndimd,cmat(iac),v)

!  contraction puis inversion des blocs diagonaux de la la matrice
!  front par front .

      do 2 nb=2,nj
!       write(*,*)'  front numero ',nb
!  adresse du bloc diagonal precedent
        iap=iac
        iac=(nb-1)*lbd+1

!      lecture du bloc sur-diagonal precedent: passage bande-vecteur

        call feti_liblos(n,a,ndimd,nb-1,w)

!  elimination du bloc extra-diagonal .

!   lecture du bloc diagonal courant: passage bande-plein

        call feti_liblod(n,ni,a,ndimd,nb,cmat(iac))

!  calcul du nouveau bloc diagonal courant .

!       write(*,*)'  calcul du nouveau bloc diagonal courant '
        call feti_nbdia(ndimd,w,cmat(iap),cmat(iac))

!  inversion du nouveau bloc diagonal courant .

!       write(*,*)'  factorisation du nouveau bloc diagonal courant '
        j0=(nb-1)*ndimd
        call feti_ldlt(ndimd,cmat(iac),d((nb-1)*ndimd+1),w,ndlblo,   &
            lisblo,nmdlblo,j0)
        call feti_calinv(ndimd,cmat(iac),v)

 2    continue

      return
      end
!**********************************************************************c
      subroutine feti_calinv(dim,a,b)
!  calcul de la partie triangulaire inferieure du bloc :
!  b = [a]-1 puis recopie dans a complet .
      implicit none
      integer dim
      REAL(kind=8)  a(dim,dim),b(dim,dim)
      integer i,j

      call feti_vclr(dim*dim,b)
      do j=1,dim
        b(j,j)=1.e0
        call feti_desremo(dim,j-1,a,b(1,j),b(1,j))
      end do
      do j=1,dim
        do i=j,dim
          a(i,j)=b(i,j)
          a(j,i)=b(i,j)
        end do
      end do

      return
      end
!***********************************************************************
      subroutine feti_nbdia(n,d,w,v)
      implicit REAL(kind=8) (a-h,o-z)
      dimension d(n),v(n,n),w(n,n)

      do 1 j=1,n
        do 2 i=1,n
          v(i,j)=v(i,j)-d(i)*w(i,j)*d(j)
 2      continue
 1    continue

      return
      end
!***********************************************************************
      subroutine feti_liblod(n,ni,a,ndimd,nb,db)

      implicit REAL(kind=8) (a-h,o-z)
      dimension db(ndimd,ndimd)
!  tableaux descripteurs de la structure morse de la matrice .
      dimension a(n,5)

      call feti_vclr(ndimd*ndimd,db)
!       do i=1,ndimd*ndimd
!         db(i,1)=0.e0
!8     continue
      n0=(nb-1)*ndimd
!  lecture de la diagonale db(i,i) .
      do 3 i=1,ndimd
        db(i,i)=a(n0+i,3)
 3    continue
!  lecture de la premiere sur-diagonale db(i,i+1) .
      do 4 i=1,ndimd-1
        db(i,i+1)=a(n0+i,4)
 4    continue
!  lecture de la premiere sous-diagonale db(i,i-1) .
      do 2 i=2,ndimd
        db(i,i-1)=a(n0+i,2)
 2    continue

      return
      end
!***********************************************************************
      subroutine feti_liblos(n,a,ndimd,nb,v)

      implicit REAL(kind=8) (a-h,o-z)
      dimension v(ndimd)
!  tableaux descripteurs de la structure morse de la matrice .
      dimension a(n,5)

      n0=(nb-1)*ndimd
!  lecture de la deuxieme sur-diagonale v(i) .
      do 5 i=1,ndimd
        v(i)=a(n0+i,5)
 5    continue

      return
      end
!***********************************************************************
      subroutine feti_resloc(n,ni,nj,a,lcmat,cmat,y,x,z)
      implicit REAL(kind=8) (a-h,o-z)
      dimension a(n,5)
!  routine de resolution utilisant la decomposition frontale par blocs
!  de la matrice heptadiagonale a .
!  cmat contient les termes des blocs diagonaux factorises .
      dimension cmat(lcmat)
!  d contient l'inverse de la diagonale de Crout .
      dimension y(n),x(n),z(n)
!
      call feti_vmov(n,y,z)
      ndimd=ni
      lbd=ndimd*ndimd
!  descente du systeme .
      do nb=1,nj-1
!  calcul du second membre condense : z(i+1)=y(i+1)- FD-1 z(i) .
        iacmat=(nb-1)*lbd+1
        n0=(nb-1)*ndimd
        call feti_mxv(ndimd,ndimd,cmat(iacmat),z(n0+1),x(n0+1))
        call feti_pbloi(n,ni,nj,a,x,z,nb)
      end do
!  inversion du dernier bloc diagonal .
      iacmat=(nj-1)*lbd+1
      n0=(nj-1)*ndimd
      call feti_mxv(ndimd,ndimd,cmat(iacmat),z(n0+1),x(n0+1))
!  remontee du systeme par bloc .
      do nb=nj-1,1,-1
!  calcul du second membre condense : z(i)=z(i)- E x(i+1) .
        call feti_pblos(n,ni,nj,a,x,z,nb)
!  calcul de la solution pour le bloc diagonal courant .
        iacmat=(nb-1)*lbd+1
        n0=(nb-1)*ndimd
        call feti_mxv(ndimd,ndimd,cmat(iacmat),z(n0+1),x(n0+1))
      end do

      end
!***********************************************************************
      subroutine feti_pblos(n,ni,nj,a,x,y,nb)
!  calcul du produit matrice-vecteur morse : y = y - A x .
!  produit par le bloc superieur numero nb .
      implicit REAL(kind=8) (a-h,o-z)
      dimension a(n,5)
      dimension x(n),y(n)

      n0=(nb-1)*ni
!  produit par la deuxieme sur-diagonale a(i,i+ni) .
      do i=n0+1,n0+ni
        y(i)=y(i)-x(i+ni)*a(i,5)
      end do

      end subroutine feti_pblos
!***********************************************************************
      subroutine feti_pbloi(n,ni,nj,a,x,y,nb)
!  calcul du produit matrice-vecteur morse : y = y - A x .
!  produit par le bloc inferieur numero nb .
      implicit REAL(kind=8) (a-h,o-z)
      dimension a(n,5)
      dimension x(n),y(n)

      n0=nb*ni
!  produit par la deuxieme sous-diagonale a(i,i-ni) .
      do i=n0+1,n0+ni
        y(i)=y(i)-x(i-ni)*a(i,1)
      end do

      end subroutine feti_pbloi
!***********************************************************************
      subroutine feti_proax(n,ni,nj,a,x,y)
!  calcul du produit matrice-vecteur morse par la matrice heptadiagonale
!  a stockee par bandes .
      implicit REAL(kind=8) (a-h,o-z)
      dimension a(n,5)
      dimension x(n),y(n)

!  produit par la diagonale .
      do i=1,n
        y(i)=x(i)*a(i,3)
      end do
!  produit par la premiere sur-diagonale a(i,i+1) .
      do i=1,n-1
        y(i)=y(i)+x(i+1)*a(i,4)
      end do
!  produit par la deuxieme sur-diagonale a(i,i+ni) .
      do i=1,n-ni
        y(i)=y(i)+x(i+ni)*a(i,5)
      end do
! produit par la premiere sous-diagonale a(i,i-1)
      do i=2,n
          y(i)=y(i)+x(i-1)*a(i,2)
      end do
! produit par la deuxieme sous-diagonale a(i,i-ni)
      do i=1+ni,n
          y(i)=y(i)+x(i-ni)*a(i,1)
      end do

      end subroutine feti_proax

#else
   !! no use of FETI librairy
#endif

   !!======================================================================
END MODULE lib_feti
