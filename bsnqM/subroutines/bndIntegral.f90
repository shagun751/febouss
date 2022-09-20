subroutine bndIntegral1(npl,npt,nele,nbnd,conn,mabnd,Sz,ivl,ivq,&
  linkl,linkq,invJ,bndS,dep,gFPP,gFPQ,gFQP,gFQQ,gFW)
use bsnqGlobVars
implicit none

  integer(kind=C_K1),intent(in)::npl,npt,nele,nbnd
  integer(kind=C_K1),intent(in)::Sz(4),conn(nele,6)
  integer(kind=C_K1),intent(in)::ivl(0:npt),ivq(0:npt)
  integer(kind=C_K1),intent(in)::linkl(Sz(3))
  integer(kind=C_K1),intent(in)::linkq(Sz(4))
  integer(kind=C_K1),intent(in)::mabnd(nbnd,6)
  integer(kind=C_K1)::i,j,k,ele,en(6),gRow,gCol,lRow,lCol,sTyp
  integer(kind=C_K1)::nlinkl(ivl(0)),nlinkq(ivq(0))
  
  real(kind=C_K2),intent(in)::invJ(nele,5),bndS(nbnd,3),dep(npt)
  real(kind=C_K2),intent(out)::gFPP(Sz(4)),gFPQ(Sz(4))
  real(kind=C_K2),intent(out)::gFQP(Sz(4)),gFQQ(Sz(4))
  real(kind=C_K2),intent(out)::gFW(Sz(1))
  real(kind=C_K2)::cnst(10)
  real(kind=C_K2)::nx,ny,sL
  real(kind=C_K2)::lFPP(6,6),lFPQ(6,6),tmp(6,6)
  real(kind=C_K2)::lFQP(6,6),lFQQ(6,6),lFW(3,3),tmp3x3(3,3)
  real(kind=C_K2)::b11,b12,b21,b22


  cnst(1)=grav                    !gravity
  cnst(2)=BsqC                    !Bsnq constant
  cnst(3)=BsqC+(1d0/3d0)          !dervied 1
  cnst(4)=2d0*BsqC+(1d0/3d0)      !derived 2
  cnst(5)=2d0*BsqC+(1d0/2d0)      !derived 3
  cnst(6)=grav*BsqC               !derived 4

  gFPP=0d0
  gFPQ=0d0
  gFQP=0d0
  gFQQ=0d0
  gFW=0d0

  do i=1,nbnd
    ele=mabnd(i,3)
    en=conn(ele,:)

    b11=invJ(ele,1)
    b12=invJ(ele,2)
    b21=invJ(ele,3)
    b22=invJ(ele,4)
    nx=bndS(i,1)
    ny=bndS(i,2)
    sL=bndS(i,3)
    sTyp=mabnd(i,4)    

    if(mabnd(i,6).eq.1) then
      call bsnqBnd12(lFPP,lFPQ,lFQP,lFQQ,cnst,b11,b12,&
        b21,b22,nx,ny,dep(en(1)),dep(en(2)),dep(en(3)),tmp)
    else if(mabnd(i,6).eq.2) then
      call bsnqBnd23(lFPP,lFPQ,lFQP,lFQQ,cnst,b11,b12,&
        b21,b22,nx,ny,dep(en(1)),dep(en(2)),dep(en(3)),tmp)
    else if(mabnd(i,6).eq.3) then
      call bsnqBnd31(lFPP,lFPQ,lFQP,lFQQ,cnst,b11,b12,&
        b21,b22,nx,ny,dep(en(1)),dep(en(2)),dep(en(3)),tmp)
    endif

    lFW=0d0
    if((sTyp.eq.11).or.(sTyp.eq.14))then
      call femBnd_N3i_Sc3_dN3jdx(tmp3x3,dep(en(1)),dep(en(2)),&
        dep(en(3)),b11,b12,mabnd(i,6))
      lFW = tmp3x3*nx*sL
      call femBnd_N3i_Sc3_dN3jdx(tmp3x3,dep(en(1)),dep(en(2)),&
        dep(en(3)),b21,b22,mabnd(i,6))
      lFW = lFW + tmp3x3*ny*sL
    endif

    lFPP=-lFPP*sL
    lFPQ=-lFPQ*sL
    lFQP=-lFQP*sL
    lFQQ=-lFQQ*sL

    !6x6
    do lRow=1,6
      gRow=en(lRow)
      k=(gRow-1)*ivq(0)
      nlinkq=linkq(k+1:k+ivq(0))
      do lCol=1,6
        gCol=en(lCol)
        do j=1,ivq(gRow)
          if(nlinkq(j).eq.gCol) goto 11
        enddo
        write(9,*)"[Err] node conn missing in Bsnq at",gRow
        stop
        11 gFPP(k+j)=gFPP(k+j)+lFPP(lRow,lCol)
        gFPQ(k+j)=gFPQ(k+j)+lFPQ(lRow,lCol)
        gFQP(k+j)=gFQP(k+j)+lFQP(lRow,lCol)
        gFQQ(k+j)=gFQQ(k+j)+lFQQ(lRow,lCol)
      enddo
    enddo    

    if((sTyp.eq.11).or.(sTyp.eq.14))then
      continue
    else
      cycle
    endif

    !3x3
    do lRow=1,3
      gRow=en(lRow)
      k=(gRow-1)*ivl(0)
      nlinkl=linkl(k+1:k+ivl(0))
      do lCol=1,3
        gCol=en(lCol)
        do j=1,ivl(gRow)
          if(nlinkl(j).eq.gCol) goto 14
        enddo
        write(9,*)"[Err] node conn missing in gFW at",gRow
        stop
        14 gFW(k+j)=gFW(k+j) + lFW(lRow,lCol)
      enddo
    enddo

  enddo

end subroutine bndIntegral1



subroutine bsnqBnd12(lFPP,lFPQ,lFQP,lFQQ,cnst,&
  b11,b12,b21,b22,nx,ny,h1,h2,h3,tmp)
use bsnqGlobVars
implicit none
  
  real(kind=C_K2),intent(out)::lFPP(6,6),lFPQ(6,6),tmp(6,6)
  real(kind=C_K2),intent(out)::lFQP(6,6),lFQQ(6,6)
  real(kind=C_K2),intent(in)::cnst(10),b11,b12,b21,b22
  real(kind=C_K2),intent(in)::nx,ny,h1,h2,h3  


  call bsnqBnd_N6i_hsq_dN6jdx_S12(tmp,b11,b12,h1,h2,h3)
  lFPP=(cnst(3)*nx*tmp)
  lFQP=(cnst(3)*ny*tmp)

  call bsnqBnd_N6i_hsq_dN6jdx_S12(tmp,b21,b22,h1,h2,h3)
  lFPQ=(cnst(3)*nx*tmp)
  lFQQ=(cnst(3)*ny*tmp)

end subroutine bsnqBnd12



subroutine bsnqBnd23(lFPP,lFPQ,lFQP,lFQQ,cnst,&
  b11,b12,b21,b22,nx,ny,h1,h2,h3,tmp)
use bsnqGlobVars
implicit none
  
  real(kind=C_K2),intent(out)::lFPP(6,6),lFPQ(6,6),tmp(6,6)
  real(kind=C_K2),intent(out)::lFQP(6,6),lFQQ(6,6)
  real(kind=C_K2),intent(in)::cnst(10),b11,b12,b21,b22
  real(kind=C_K2),intent(in)::nx,ny,h1,h2,h3  


  call bsnqBnd_N6i_hsq_dN6jdx_S23(tmp,b11,b12,h1,h2,h3)
  lFPP=(cnst(3)*nx*tmp)
  lFQP=(cnst(3)*ny*tmp)

  call bsnqBnd_N6i_hsq_dN6jdx_S23(tmp,b21,b22,h1,h2,h3)
  lFPQ=(cnst(3)*nx*tmp)
  lFQQ=(cnst(3)*ny*tmp)

end subroutine bsnqBnd23



subroutine bsnqBnd31(lFPP,lFPQ,lFQP,lFQQ,cnst,&
  b11,b12,b21,b22,nx,ny,h1,h2,h3,tmp)
use bsnqGlobVars
implicit none
  
  real(kind=C_K2),intent(out)::lFPP(6,6),lFPQ(6,6),tmp(6,6)
  real(kind=C_K2),intent(out)::lFQP(6,6),lFQQ(6,6)
  real(kind=C_K2),intent(in)::cnst(10),b11,b12,b21,b22
  real(kind=C_K2),intent(in)::nx,ny,h1,h2,h3  


  call bsnqBnd_N6i_hsq_dN6jdx_S31(tmp,b11,b12,h1,h2,h3)
  lFPP=(cnst(3)*nx*tmp)
  lFQP=(cnst(3)*ny*tmp)

  call bsnqBnd_N6i_hsq_dN6jdx_S31(tmp,b21,b22,h1,h2,h3)
  lFPQ=(cnst(3)*nx*tmp)
  lFQQ=(cnst(3)*ny*tmp)

end subroutine bsnqBnd31



subroutine bsnqBnd_N6i_hsq_dN6jdx_S12(lF,b11,b12,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::lF(6,6)
  real(kind=C_K2),intent(in)::b11,b12,h1,h2,h3

  lF=0d0

  lF(1,1)=(-(1d0/60d0))*(b11 + b12)*(23d0*h1**2 + 6d0*h1*h2 + h2**2) 
  lF(1,2)=-((b11*h1**2)/12d0) - (b11*h1*h2)/30d0 - (b11*h2**2)/20d0 
  lF(1,3)=-((3d0*b12*h1**2)/20d0) - (b12*h1*h2)/30d0 &
    + (b12*h2**2)/60d0
  lF(1,4)=(7d0*b11*h1**2)/15d0 - (b12*h1**2)/15d0 &
    + (2d0*b11*h1*h2)/15d0 + (b11*h2**2)/15d0 + (b12*h2**2)/15d0 
  lF(1,5)=(b12*h1**2)/15d0 - (b12*h2**2)/15d0 
  lF(1,6)=(8d0*b12*h1**2)/15d0 + (2d0*b12*h1*h2)/15d0 

  lF(2,1)=(1d0/60d0)*(b11 + b12)*(3d0*h1**2 + 2d0*h1*h2 + 5d0*h2**2)
  lF(2,2)=(b11*h1**2)/60d0 + (b11*h1*h2)/10d0 + (23d0*b11*h2**2)/60d0
  lF(2,3)=(b12*h1**2)/60d0 - (b12*h1*h2)/30d0 - (3d0*b12*h2**2)/20d0 
  lF(2,4)=-((b11*h1**2)/15d0) - (2d0*b11*h1*h2)/15d0 &
    - (2d0*b12*h1*h2)/15d0 - (7d0*b11*h2**2)/15d0 &
    - (8d0*b12*h2**2)/15d0 
  lF(2,5)=(2d0*b12*h1*h2)/15d0 + (8d0*b12*h2**2)/15d0 
  lF(2,6)=-((b12*h1**2)/15d0) + (b12*h2**2)/15d0    

  lF(4,1)=(-(1d0/15d0))*(b11 + b12)*(5d0*h1**2 + 4d0*h1*h2 + h2**2)
  lF(4,2)=(b11*h1**2)/15d0 + (4d0*b11*h1*h2)/15d0 + (b11*h2**2)/3d0 
  lF(4,3)=-((b12*h1**2)/5d0) - (4d0*b12*h1*h2)/15d0 - (b12*h2**2)/5d0 
  lF(4,4)=(4d0*b11*h1**2)/15d0 - (4d0*b12*h1**2)/15d0 &
    - (8d0*b12*h1*h2)/15d0 - (4d0*b11*h2**2)/15d0 &
    - (8d0*b12*h2**2)/15d0
  lF(4,5)=(4d0*b12*h1**2)/15d0 + (8d0*b12*h1*h2)/15d0 &
    + (8d0*b12*h2**2)/15d0 
  lF(4,6)=(8d0*b12*h1**2)/15d0 + (8d0*b12*h1*h2)/15d0 &
    + (4d0*b12*h2**2)/15d0 

end subroutine bsnqBnd_N6i_hsq_dN6jdx_S12



subroutine bsnqBnd_N6i_hsq_dN6jdx_S23(lF,b11,b12,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::lF(6,6)
  real(kind=C_K2),intent(in)::b11,b12,h1,h2,h3

  lF=0d0

  lF(2,1)=(1d0/60d0)*(b11 + b12)*(9d0*h2**2 + 2d0*h2*h3 - h3**2) 
  lF(2,2)=(23d0*b11*h2**2)/60d0 + (b11*h2*h3)/10d0 + (b11*h3**2)/60d0 
  lF(2,3)=-((b12*h2**2)/12d0) - (b12*h2*h3)/30d0 - (b12*h3**2)/20d0 
  lF(2,4)=(-(2d0/15d0))*(b11 + b12)*h2*(4d0*h2 + h3)
  lF(2,5)=(b11*h2**2)/15d0 + (8d0*b12*h2**2)/15d0 &
      + (2d0*b12*h2*h3)/15d0 - (b11*h3**2)/15d0 
  lF(2,6)=(1d0/15d0)*(b11 + b12)*(-h2**2 + h3**2) 

  lF(3,1)=(-(1d0/60d0))*(b11 + b12)*(h2**2 - 2d0*h2*h3 - 9d0*h3**2)
  lF(3,2)=-((b11*h2**2)/20d0) - (b11*h2*h3)/30d0 - (b11*h3**2)/12d0 
  lF(3,3)=(b12*h2**2)/60d0 + (b12*h2*h3)/10d0 + (23d0*b12*h3**2)/60d0 
  lF(3,4)=(1d0/15d0)*(b11 + b12)*(h2**2 - h3**2) 
  lF(3,5)=-((b12*h2**2)/15d0) + (2d0*b11*h2*h3)/15d0 &
    + (8d0*b11*h3**2)/15d0 + (b12*h3**2)/15d0 
  lF(3,6)=(-(2d0/15d0))*(b11 + b12)*h3*(h2 + 4d0*h3)

  lF(5,1)=(1d0/15d0)*(b11 + b12)*(3d0*h2**2 + 4d0*h2*h3 + 3d0*h3**2)
  lF(5,2)=(b11*h2**2)/3d0 + (4d0*b11*h2*h3)/15d0 + (b11*h3**2)/15d0 
  lF(5,3)=(b12*h2**2)/15d0 + (4d0*b12*h2*h3)/15d0 + (b12*h3**2)/3d0 
  lF(5,4)=(-(4d0/15d0))*(b11 + b12)*(2d0*h2**2 + 2d0*h2*h3 + h3**2) 
  lF(5,5)=(4d0*b11*h2**2)/15d0 + (8d0*b12*h2**2)/15d0 &
    + (8d0*b11*h2*h3)/15d0 + (8d0*b12*h2*h3)/15d0 &
    + (8d0*b11*h3**2)/15d0 + (4d0*b12*h3**2)/15d0 
  lF(5,6)=(-(4d0/15d0))*(b11 + b12)*(h2**2 + 2d0*h2*h3 + 2d0*h3**2)

end subroutine bsnqBnd_N6i_hsq_dN6jdx_S23



subroutine bsnqBnd_N6i_hsq_dN6jdx_S31(lF,b11,b12,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::lF(6,6)
  real(kind=C_K2),intent(in)::b11,b12,h1,h2,h3

  lF=0d0

  lF(1,1)=(-(1d0/60d0))*(b11 + b12)*(23d0*h1**2 + 6d0*h1*h3 + h3**2)
  lF(1,2)=-((3d0*b11*h1**2)/20d0) - (b11*h1*h3)/30d0 &
    + (b11*h3**2)/60d0 
  lF(1,3)=-((b12*h1**2)/12d0) - (b12*h1*h3)/30d0 - (b12*h3**2)/20d0 
  lF(1,4)=(8d0*b11*h1**2)/15d0 + (2d0*b11*h1*h3)/15d0 
  lF(1,5)=(b11*h1**2)/15d0 - (b11*h3**2)/15d0 
  lF(1,6)=-((b11*h1**2)/15d0) + (7d0*b12*h1**2)/15d0 &
    + (2d0*b12*h1*h3)/15d0 + (b11*h3**2)/15d0 + (b12*h3**2)/15d0 

  lF(3,1)=(1d0/60d0)*(b11 + b12)*(3d0*h1**2 + 2d0*h1*h3 + 5d0*h3**2)
  lF(3,2)=(b11*h1**2)/60d0 - (b11*h1*h3)/30d0 - (3d0*b11*h3**2)/20d0 
  lF(3,3)=(b12*h1**2)/60d0 + (b12*h1*h3)/10d0 + (23d0*b12*h3**2)/60d0 
  lF(3,4)=-((b11*h1**2)/15d0) + (b11*h3**2)/15d0 
  lF(3,5)=(2d0*b11*h1*h3)/15d0 + (8d0*b11*h3**2)/15d0 
  lF(3,6)=-((b12*h1**2)/15d0) - (2d0*b11*h1*h3)/15d0 &
    - (2d0*b12*h1*h3)/15d0 - (8d0*b11*h3**2)/15d0 &
    - (7d0*b12*h3**2)/15d0 

  lF(6,1)=(-(1d0/15d0))*(b11 + b12)*(5d0*h1**2 + 4d0*h1*h3 + h3**2)
  lF(6,2)=-((b11*h1**2)/5d0) - (4d0*b11*h1*h3)/15d0 - (b11*h3**2)/5d0 
  lF(6,3)=(b12*h1**2)/15d0 + (4d0*b12*h1*h3)/15d0 + (b12*h3**2)/3d0 
  lF(6,4)=(8d0*b11*h1**2)/15d0 + (8d0*b11*h1*h3)/15d0 &
    + (4d0*b11*h3**2)/15d0 
  lF(6,5)=(4d0*b11*h1**2)/15d0 + (8d0*b11*h1*h3)/15d0 &
    + (8d0*b11*h3**2)/15d0 
  lF(6,6)=-((4d0*b11*h1**2)/15d0) + (4d0*b12*h1**2)/15d0 &
    - (8d0*b11*h1*h3)/15d0 - (8d0*b11*h3**2)/15d0 &
    - (4d0*b12*h3**2)/15d0

end subroutine bsnqBnd_N6i_hsq_dN6jdx_S31