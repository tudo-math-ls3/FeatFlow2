      subroutine  vanca_full 
     *    (u1,u2,p,f1,f2,fp,
     *    a,kcola,klda,b1,b2,kcolb,kldb,nu,np,
     *    kmbd,kvert,kmid,knpr,nmbd)
************************************************************************
c
c one iteration of vanka smother (block gauss-seidel) on the system
c
c    A  0  B1 | U1    F1
c    0  A  B2 | U2  = F2
c    D1 D2 0  | PR    FP
c
c Steps:
c loop over all elements
c   compute local residuum ff
c   extract local matrix aa
c   solve for xx= aa^(-1) ff
c   update global solution vector (u1,u2,p)=(u1,u2,p)+ rlxsm*xx
c end loop
c
c This version will be able to solve the real full system of the form
c
c    A1 A2 B1 | U1    F1
c    A3 A4 B2 | U2  = F2
c    D1 D2 0  | PR    FP
c
c once we have full coupling and even if coupled with additional 
c equations as
c
      implicit none
      
c common blocks      
      
      include 'cout.inc'
      
      include 'cbasictria.inc'
      include 'ctria.inc'
      
      include 'cinidat.inc'
      
c parameters

      double precision u1(*),u2(*),p(*)
      double precision f1(*),f2(*),fp(*)
      real a(*),b1(*),b2(*)
      integer kcola(*),klda(*),kcolb(*),kldb(*)
      integer kmbd(*),kvert(nnve,*),kmid(nnve,*),knpr(*)
      integer nu,np,nmbd
      
c *** local arrays for informations about one element
      
      integer nnld
      parameter(nnld=2*4+1)
      double precision aa(nnld,nnld),ff(nnld)
      integer ipiv(nnld)
      integer iu(4)
      logical bdbc(4),bnbc(4)
      
c local variables

      integer iel,imid,ia1,ia2,ia,ib,ib1,ib2,info
      integer ii,jj,i,j,k,jp,jp1,jp2,iu1,iu2,ipr,nn
      double precision ffp,pjp,aoff,a1,a2,aux1,aux2,aux

c-----------------------------------------------------------------------
c pointers to the solution components in the local vector
      iu1=0
      iu2=iu1+4
      ipr=iu2+4
c total number of unknowns on each element
      nn=nnld

c=======================================================================
c     block gauss-seidel on schur complement
c=======================================================================
c *** loop over all elements
c
      do iel=1,nel

c get the local numbering and the BC info for dirichlet nodes
        do ii=1,4
          imid=kmid(ii,iel)
          i=imid-nvt
          iu(ii)=i
          if (knpr(imid).eq.0)  then
            bdbc(ii)=.false.
          else
            bdbc(ii)=.true.
          endif
        enddo

c compute the local residuum ff=f-Ax
c * get the local rhs
        ff(ipr+1)=fp(iel)
        do ii=1,4
          i=iu(ii)
          ff(iu1+ii)=f1(i)
          ff(iu2+ii)=f2(i)
        enddo

C * compute ff=ff - A x for the local unknowns
        do ii=1,4          
          i=iu(ii)
          ia1=klda(i)
          ia2=klda(i+1)-1
          do ia=ia1,ia2
            j=kcola(ia)
            aux=dble(a(ia))
            ff(1,iu1+ii)=ff(1,iu1+ii)-aux*u1(j)
            ff(1,iu2+ii)=ff(1,iu2+ii)-aux*u2(j)
          enddo
          ib1=kldb(i)
          ib2=kldb(i+1)-1
          do ib=ib1,ib2
            k=kcolb(ib)
            aux1=dble(b1(ib))
            aux2=dble(b2(ib))
            if(iel.eq.k) ff(1,ipr+1)=ff(1,ipr+1)-aux1*u1(i)-aux2*u2(i)
            if(bdbc(ii)) aux1=0d0 
            if(bdbc(ii)) aux2=0d0
            ff(1,iu1+ii)=ff(1,iu1+ii)-aux1*p(k)
            ff(1,iu2+ii)=ff(1,iu2+ii)-aux2*p(k)
          enddo
        enddo
          
c extract the local matrix into one array aa(nn,nn)
        aa(ipr+1,ipr+1)=0d0
        do ii=1,4
          i=iu(ii)
          do jj=1,4
            j=iu(jj)
            aa(iu1+ii,iu1+jj)=0d0
            aa(iu2+ii,iu2+jj)=0d0
            aa(iu1+ii,iu2+jj)=0d0
            aa(iu2+ii,iu1+jj)=0d0
            ia1=klda(i)
            ia2=klda(i+1)-1
            do ia=ia1,ia2
              if(j.eq.kcola(ia)) then
                aux=dble(a(ia))
                aa(iu1+ii,iu1+jj)=aux
                aa(iu2+ii,iu2+jj)=aux
              endif
            enddo
          enddo
          ib1=kldb(i)
          ib2=kldb(i+1)-1
          do ib=ib1,ib2
            k=kcolb(ib)
            aux1=dble(b1(ib))
            aux2=dble(b2(ib))
            if(iel.eq.k) then
              aa(ipr+1,iu1+ii)=aux1
              aa(ipr+1,iu2+ii)=aux2
              if(bdbc(ii)) aux1=0d0 
              if(bdbc(ii)) aux2=0d0
              aa(iu1+ii,ipr+1)=aux1
              aa(iu2+ii,ipr+1)=aux2           
            endif
          enddo
        enddo

c solve the local system aa xx = ff using lapack
        call DGETRF( nn, nn, aa, nn, ipiv, info )
        if(info.ne.0) print *,'ERROR: LAPACK(DGETRF) LU decomposition'
        call DGETRS('N', nn, 1, aa, nn, ipiv, ff, nn, info )
        if(info.ne.0) print *,'ERROR: LAPACK(DGETRS) back substitution'

c global update/correction (u1,u2,p)=(u1,u2,p) + rlxsm*xx
        do ii=1,4
          i=iu(ii)
          u1(i)=u1(i)+rlxsm*ff(1,iu1+ii)
          u2(i)=u2(i)+rlxsm*ff(1,iu2+ii)
        enddo
        p(iel)=p(iel)+rlxsm*ff(1,ipr+1)

      enddo ! over all elements

      end

*=======================================================================
*=======================================================================

      subroutine  vanca_diag 
     *    (u1,u2,p,f1,f2,fp,
     *    a,kcola,klda,b1,b2,kcolb,kldb,nu,np,
     *    kmbd,kvert,kmid,knpr,nmbd)
************************************************************************
c
c one iteration of vanka smother (block gauss-seidel) on the system
c
c    A  0  B1 | U1    F1
c    0  A  B2 | U2  = F2
c    D1 D2 0  | PR    FP
c
c Steps:
c loop over all elements
c   compute local residuum ff
c   extract local matrix aa with only the diagonal part of the block A
c   solve for xx= aa^(-1) ff
c   update global solution vector (u1,u2,p)=(u1,u2,p)+ rlxsm*xx
c end loop
c
c
c
      implicit none
      
c common blocks      
      
      include 'cout.inc'
      
      include 'cbasictria.inc'
      include 'ctria.inc'
      
      include 'cinidat.inc'
      
c parameters

      double precision u1(*),u2(*),p(*)
      double precision f1(*),f2(*),fp(*)
      real a(*),b1(*),b2(*)
      integer kcola(*),klda(*),kcolb(*),kldb(*)
      integer kmbd(*),kvert(nnve,*),kmid(nnve,*),knpr(*)
      integer nu,np,nmbd
      
c *** local arrays for informations about one element
      
      integer nnld
      parameter(nnld=2*4)
      double precision aa(nnld),bb(nnld),dd(nnld)
      double precision xx(nnld+1), ff(nnld+1)
      integer ipiv(nnld)
      integer iu(4)
      logical bdbc(4),bnbc(4)
      
c local variables

      integer iel,imid,ia1,ia2,ia,ib,ib1,ib2,info
      integer ii,jj,i,j,k,jp,jp1,jp2,iu1,iu2,ipr,nn
      double precision ffp,pjp,aoff,a1,a2,aux1,aux2,aux

c-----------------------------------------------------------------------
      iu1=0
      iu2=iu1+4
      ipr=iu2+4
      nn=nnld
c=======================================================================
c     block gauss-seidel on schur complement
c=======================================================================
c *** loop over all elements
c
      do iel=1,nel

c get the local numbering and the BC info for dirichlet nodes
        do ii=1,4
          imid=kmid(ii,iel)
          i=imid-nvt
          iu(ii)=i
          if (knpr(imid).eq.0)  then
            bdbc(ii)=.false.
          else
            bdbc(ii)=.true.
          endif
        enddo

c compute the local residuum
        ff(ipr+1)=fp(iel)
        do ii=1,4          
          i=iu(ii)
          ff(iu1+ii)=f1(i)
          ff(iu2+ii)=f2(i)

          ia=klda(i)
          aux=1d0/dble(a(ia))
          aa(iu1+ii)=aux
          aa(iu2+ii)=aux

          ia1=klda(i)
          ia2=klda(i+1)-1
          do ia=ia1,ia2
            j=kcola(ia)
            aux=dble(a(ia))
            ff(iu1+ii)=ff(iu1+ii)-aux*u1(j)
            ff(iu2+ii)=ff(iu2+ii)-aux*u2(j)
          enddo

          ib1=kldb(i)
          ib2=kldb(i+1)-1
          do ib=ib1,ib2
            k=kcolb(ib)
            aux1=dble(b1(ib))
            aux2=dble(b2(ib))
            if(iel.eq.k) then
              ff(ipr+1)=ff(ipr+1)-aux1*u1(i)-aux2*u2(i)
              dd(iu1+ii)=aux1
              dd(iu2+ii)=aux2
              if(.not.bdbc(ii)) then
                bb(iu1+ii)=aux1
                bb(iu2+ii)=aux2           
                else 
                bb(iu1+ii)=0d0
                bb(iu2+ii)=0d0
              endif
            endif
            if(.not.bdbc(ii)) then
              ff(iu1+ii)=ff(iu1+ii)-aux1*p(k)
              ff(iu2+ii)=ff(iu2+ii)-aux2*p(k)
            endif
          enddo
        enddo
          

c solve the local system aa(nn,nn) xx = ff of the form
c
c     |aa 0  bb|
c     |0  aa bb| xx = ff 
c     |dd dd 0 |
c
c by schur complement method
c
        aux1=0d0
        aux2=0d0
        do ii=1,8
          aux1=aux1+dd(ii)*ff(ii)*aa(ii)
          aux2=aux2+dd(ii)*bb(ii)*aa(ii)
        enddo
        xx(ipr+1)=(aux1-ff(ipr+1))/aux2
        do ii=1,8
          xx(ii)=(ff(ii)-bb(ii)*xx(ipr+1))*aa(ii)
        enddo

c global update (u1,u2,p)=(u1,u2,p) + rlxsm*xx
        do ii=1,4
          i=iu(ii)
          u1(i)=u1(i)+rlxsm*xx(iu1+ii)
          u2(i)=u2(i)+rlxsm*xx(iu2+ii)
        enddo
        p(iel)=p(iel)+rlxsm*xx(ipr+1)

      enddo ! over all elements

      end



C
