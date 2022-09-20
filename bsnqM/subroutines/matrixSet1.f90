subroutine matrixSet1(npoinl,npoint,nelem,conn,Sz,ivl,ivq,&
  linkl,linkq,invJ,depth,por,mass1,mass2,gBs1,gBs2,gBs3,gBs4,&
  gCxFlux,gCyFlux,gDMat,gBs5,gBs6,ele6x6,ele6x3)
use bsnqGlobVars
implicit none

  !Mass matrices M1(6x6) and M2(3x3)
  !Bsnq matrices 1,2,3,4
  !c Flux gradient mtrices x,y
  !D matrix
  
  integer(kind=C_K1),intent(in)::npoinl,npoint,nelem,Sz(4)
  integer(kind=C_K1),intent(in)::conn(nelem,6)
  integer(kind=C_K1),intent(in)::ivl(0:npoint),linkl(Sz(3))
  integer(kind=C_K1),intent(in)::ivq(0:npoint),linkq(Sz(4))
  integer(kind=C_K1),intent(out)::ele6x6(nelem,36)
  integer(kind=C_K1),intent(out)::ele6x3(nelem,18)
  integer(kind=C_K1)::i,j,k,i2,j2,k2,n(6),gRow,gCol,lRow,lCol
  integer(kind=C_K1)::nlinkl(ivl(0)),nlinkq(ivq(0))

  real(kind=C_K2),intent(in)::invJ(nelem,5),depth(npoint)  
  real(kind=C_K2),intent(in)::por(npoint)
  real(kind=C_K2),intent(out)::mass1(Sz(4))
  real(kind=C_K2),intent(out)::mass2(Sz(1))  
  real(kind=C_K2),intent(out)::gBs1(Sz(4)),gBs2(Sz(4))
  real(kind=C_K2),intent(out)::gBs3(Sz(4)),gBs4(Sz(4))
  real(kind=C_K2),intent(out)::gCxFlux(Sz(2))
  real(kind=C_K2),intent(out)::gCyFlux(Sz(2))
  real(kind=C_K2),intent(out)::gDMat(Sz(1))    
  real(kind=C_K2),intent(out)::gBs5(Sz(3)),gBs6(Sz(3))
  real(kind=C_K2)::bs1t1(6,6),bs1t2(6,6)
  real(kind=C_K2)::bs2t1(6,6),bs2t2(6,6),bs2t3(6,6)
  real(kind=C_K2)::bs3t1(6,6),bs3t2(6,6),bs3t3(6,6)
  real(kind=C_K2)::bs4t1(6,6),bs4t2(6,6)  
  real(kind=C_K2)::cxFlux(3,6),cyFlux(3,6)
  real(kind=C_K2)::dMat(3,3)
  real(kind=C_K2)::lBs5(6,3),lBs6(6,3)
  real(kind=C_K2)::cnst(10)
  real(kind=C_K2)::locM1(6,6),locM2(3,3),lScN3(3)

  cnst(1)=grav                    !gravity
  cnst(2)=BsqC                    !Bsnq constant
  cnst(3)=BsqC+(1d0/3d0)          !dervied 1
  cnst(4)=2d0*BsqC+(1d0/3d0)      !derived 2
  cnst(5)=2d0*BsqC+(1d0/2d0)      !derived 3
  cnst(6)=grav*BsqC               !derived 4

  mass1=0d0
  mass2=0d0
  gBs1=0d0
  gBs2=0d0
  gBs3=0d0
  gBs4=0d0
  gCxFlux=0d0
  gCyFlux=0d0
  gDMat=0d0
  gBs5=0d0
  gBs6=0d0

  locM1=0d0
  locM1(1,:)=(/ 1d0/60d0,-1d0/360d0,-1d0/360d0,0d0,-1d0/90d0,0d0 /)
  locM1(2,:)=(/ -1d0/360d0,1d0/60d0,-1d0/360d0,0d0,0d0,-1d0/90d0 /)
  locM1(3,:)=(/ -1d0/360d0,-1d0/360d0,1d0/60d0,-1d0/90d0,0d0,0d0 /)
  locM1(4,:)=(/ 0d0,0d0,-1d0/90d0,4d0/45d0,2d0/45d0,2d0/45d0 /)
  locM1(5,:)=(/ -1d0/90d0,0d0,0d0,2d0/45d0,4d0/45d0,2d0/45d0 /)
  locM1(6,:)=(/ 0d0,-1d0/90d0,0d0,2d0/45d0,2d0/45d0,4d0/45d0 /)

  locM2=0d0
  locM2(1,:)=(/ 1d0/12d0,1d0/24d0,1d0/24d0 /)
  locM2(2,:)=(/ 1d0/24d0,1d0/12d0,1d0/24d0 /)
  locM2(3,:)=(/ 1d0/24d0,1d0/24d0,1d0/12d0 /) 

  do i=1,nelem
    n=conn(i,:)    
    call bsnq1Term1(bs1t1,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq1Term2(bs1t2,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq2Term1(bs2t1,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq2Term2(bs2t2,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq2Term3(bs2t3,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq3Term1(bs3t1,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    bs3t2=bs2t2
    bs3t3=bs2t3

    call bsnq4Term1(bs4t1,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq4Term2(bs4t2,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call cFluxMat(cxFlux,cyFlux,invJ(i,1),invJ(i,2),&
      invJ(i,3),invJ(i,4))

    call remainMat(dMat,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))  

    lScN3=grav*BsqC*por(n(1:3))*(depth(n(1:3))**2)
    call fem_N6i_Sc3_dN3jdx(lBs5,lScN3(1),lScN3(2),lScN3(3),&
      invJ(i,1),invJ(i,2))
    call fem_N6i_Sc3_dN3jdx(lBs6,lScN3(1),lScN3(2),lScN3(3),&
      invJ(i,3),invJ(i,4))
    

    bs1t1=invJ(i,5)*((cnst(3)*bs1t1)+(cnst(4)*bs1t2))
    bs2t1=invJ(i,5)*((cnst(3)*bs2t1)+(cnst(5)*bs2t2)-(bs2t3/6d0))
    bs3t1=invJ(i,5)*((cnst(3)*bs3t1)-(bs3t2/6d0)+(cnst(5)*bs3t3))
    bs4t1=invJ(i,5)*((cnst(3)*bs4t1)+(cnst(4)*bs4t2))
    cxFlux=-invJ(i,5)*cxFlux
    cyFlux=-invJ(i,5)*cyFlux
    dMat=-invJ(i,5)*dMat   
    lBs5=invJ(i,5)*lBs5 
    lBs6=invJ(i,5)*lBs6 

    !6x6
    do lRow=1,6
      gRow=n(lRow)
      k=(gRow-1)*ivq(0)
      nlinkq=linkq(k+1:k+ivq(0))
      do lCol=1,6
        gCol=n(lCol)
        do j=1,ivq(gRow)
          if(nlinkq(j).eq.gCol) goto 11
        enddo
        write(9,*)"[Err] node conn missing in Bsnq at",gRow
        stop
        11 mass1(k+j)=mass1(k+j)+(locM1(lRow,lCol)*invJ(i,5))
        gBs1(k+j)=gBs1(k+j)+bs1t1(lRow,lCol)
        gBs2(k+j)=gBs2(k+j)+bs2t1(lRow,lCol)
        gBs3(k+j)=gBs3(k+j)+bs3t1(lRow,lCol)
        gBs4(k+j)=gBs4(k+j)+bs4t1(lRow,lCol)        
        ele6x6(i,(lRow-1)*6+lCol)=k+j
      enddo
    enddo

    !6x3
    do lRow=1,6
      gRow=n(lRow)
      k=(gRow-1)*ivl(0)
      nlinkl=linkl(k+1:k+ivl(0))
      do lCol=1,3
        gCol=n(lCol)
        do j=1,ivl(gRow)
          if(nlinkl(j).eq.gCol) goto 12
        enddo
        write(9,*)"[Err] node conn missing in Bsnq at",gRow
        stop
        12 gBs5(k+j)=gBs5(k+j)+lBs5(lRow,lCol)
        gBs6(k+j)=gBs6(k+j)+lBs6(lRow,lCol)
        ele6x3(i,(lRow-1)*3+lCol)=k+j
      enddo
    enddo

    !3x6    
    do lRow=1,3
      gRow=n(lRow)
      k=(gRow-1)*ivq(0)
      nlinkq=linkq(k+1:k+ivq(0))
      do lCol=1,6
        gCol=n(lCol)
        do j=1,ivq(gRow)
          if(nlinkq(j).eq.gCol) goto 13
        enddo
        write(9,*)"[Err] node conn missing in C Flux at",gRow
        stop
        13 gCxFlux(k+j)=gCxFlux(k+j)+cxFlux(lRow,lCol)
        gCyFlux(k+j)=gCyFlux(k+j)+cyFlux(lRow,lCol)
      enddo
    enddo    
    
    !3x3
    do lRow=1,3
      gRow=n(lRow)
      k=(gRow-1)*ivl(0)
      nlinkl=linkl(k+1:k+ivl(0))
      do lCol=1,3
        gCol=n(lCol)
        do j=1,ivl(gRow)
          if(nlinkl(j).eq.gCol) goto 14
        enddo
        write(9,*)"[Err] node conn missing in dMat at",gRow
        stop
        14 mass2(k+j)=mass2(k+j)+(locM2(lRow,lCol)*invJ(i,5))
        gDMat(k+j)=gDMat(k+j)+dMat(lRow,lCol)        
      enddo
    enddo
  enddo

end subroutine matrixSet1


subroutine bsnq1Term1(bs1t1,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::bs1t1(6,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs1t1(1,1)=(1d0/180d0)*(b11 + b12)**2d0 *(39d0*h1**2d0 &
    + 15d0*h1*(h2 + h3) + 7d0*(h2**2d0 + h2*h3 + h3**2d0))

  bs1t1(1,2)=(1d0/180d0)*b11*(b11 + b12)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs1t1(1,3)=(1d0/180d0)*b12*(b11 + b12)*(9d0*h1**2d0 &
    + h2**2d0 + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs1t1(1,4)=(-(1d0/45d0))*(b11 + b12)*(b11*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b12*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs1t1(1,5)=(1d0/45d0)*(b11 + b12)*(b12*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs1t1(1,6)=(-(1d0/45d0))*(b11 + b12)*(b12*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs1t1(2,1)=(1d0/180d0)*b11*(b11 + b12)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs1t1(2,2)=(1d0/180d0)*b11**2d0*(7d0*h1**2d0 + 15d0*h1*h2 &
    + 39d0*h2**2d0 + 7d0*h1*h3 + 15d0*h2*h3 + 7d0*h3**2d0)

  bs1t1(2,3)=(-(1d0/180d0))*b11*b12*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs1t1(2,4)=(-(1d0/45d0))*b11*(b11*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs1t1(2,5)=(1d0/45d0)*b11*(b11*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs1t1(2,6)=(-(1d0/45d0))*b11*(b12*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b11*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs1t1(3,1)=(1d0/180d0)*b12*(b11 + b12)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs1t1(3,2)=(-(1d0/180d0))*b11*b12*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs1t1(3,3)=(1d0/180d0)*b12**2d0*(7d0*(h1**2d0 + h1*h2 &
    + h2**2d0) + 15d0*(h1 + h2)*h3 + 39d0*h3**2d0)

  bs1t1(3,4)=(-(1d0/45d0))*b12*(b11*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b12*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs1t1(3,5)=(-(1d0/45d0))*b12*(b12*(h1**2d0 + 2d0*h1*h2 &
    + 3d0*h2**2d0 - h1*h3 - 2d0*h2*h3 - 3d0*h3**2d0) &
    - b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs1t1(3,6)=(-(1d0/45d0))*b12*(b12*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs1t1(4,1)=(-(1d0/45d0))*(b11 + b12)*(b11*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b12*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs1t1(4,2)=(-(1d0/45d0))*b11*(b11*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs1t1(4,3)=(-(1d0/45d0))*b12*(b11*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b12*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs1t1(4,4)=(4d0/45d0)*(b11*b12*(-h1**2d0 + 2d0*h1*h2 &
    + 9d0*h2**2d0 + 4d0*h2*h3 + h3**2d0) + b11**2d0*(4d0*h1**2d0 &
    + 4d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(h2 + h3)) &
    + b12**2d0*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3)))

  bs1t1(4,5)=(-(4d0/45d0))*((-b11**2d0)*(h1 - h2)*(h1 + h2 + h3) &
    + b11*b12*(2d0*h1*h2 + 6d0*h2**2d0 + h1*h3 + 4d0*h2*h3 &
    + 2d0*h3**2d0) + b12**2d0*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 &
    + h3**2d0 + h1*(3d0*h2 + h3)))

  bs1t1(4,6)=(4d0/45d0)*((-b11**2d0)*(h1 - h2)*(h1 + h2 + h3) &
    - b12**2d0*(h1 - h3)*(h1 + h2 + h3) + b11*b12*(4d0*h1**2d0 &
    + 2d0*h2**2d0 + 3d0*h2*h3 + 2d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  bs1t1(5,1)=(1d0/45d0)*(b11 + b12)*(b12*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs1t1(5,2)=(1d0/45d0)*b11*(b11*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs1t1(5,3)=(-(1d0/45d0))*b12*(b12*(h1**2d0 + 2d0*h1*h2 &
    + 3d0*h2**2d0 - h1*h3 - 2d0*h2*h3 - 3d0*h3**2d0) &
    - b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs1t1(5,4)=(-(4d0/45d0))*((-b11**2d0)*(h1 - h2)*(h1 + h2 + h3) &
    + b11*b12*(2d0*h1*h2 + 6d0*h2**2d0 + h1*h3 + 4d0*h2*h3 &
    + 2d0*h3**2d0) +  b12**2d0*(h1**2d0 + 6d0*h2**2d0 &
    + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3)))

  bs1t1(5,5)=(4d0/45d0)*(b11*b12*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + b12**2d0*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3)) + b11**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs1t1(5,6)=(-(4d0/45d0))*((-b12**2d0)*(h1 - h3)*(h1 + h2 &
    + h3) + b11*b12*(h1*h2 + 2d0*h2**2d0 + 2d0*h1*h3 &
    + 4d0*h2*h3 + 6d0*h3**2) +  b11**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs1t1(6,1)=(-(1d0/45d0))*(b11 + b12)*(b12*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs1t1(6,2)=(-(1d0/45d0))*b11*(b12*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b11*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs1t1(6,3)=(-(1d0/45d0))*b12*(b12*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs1t1(6,4)=(4d0/45d0)*((-b11**2d0)*(h1 - h2)*(h1 + h2 + h3) &
    - b12**2d0*(h1 - h3)*(h1 + h2 + h3) + b11*b12*(4d0*h1**2d0 &
    + 2d0*h2**2d0 + 3d0*h2*h3 + 2d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  bs1t1(6,5)=(-(4d0/45d0))*((-b12**2d0)*(h1 - h3)*(h1 + h2 + h3) &
    + b11*b12*(h1*h2 + 2d0*h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 6d0*h3**2d0) + b11**2d0*(h1**2d0 + h2**2d0 + 3d0*h2*h3 &
    + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs1t1(6,6)=(4d0/45d0)*(b11*b12*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 &
    + 4d0*h2*h3 + 9d0*h3**2) + b12**2d0*(4d0*h1**2d0 + h2**2d0 &
    + 2d0*h2*h3 + 4d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + b11**2d0*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))) 

end subroutine bsnq1Term1


subroutine bsnq1Term2(bs1t2,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::bs1t2(6,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs1t2(1,1)=(1d0/120d0)*(b11 + b12)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + h3)

  bs1t2(1,2)=(1d0/360d0)*b11*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*h1 + 5d0*h2 + h3)

  bs1t2(1,3)=(1d0/360d0)*b12*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*h1 + h2 + 5d0*h3)

  bs1t2(1,4)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(2d0*h2 + h3) + b11*(6d0*h1 + 2d0*h2 + h3))

  bs1t2(1,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(2d0*h2 + h3) + b11*(h2 + 2d0*h3))

  bs1t2(1,6)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b11*(h2 + 2d0*h3) + b12*(6d0*h1 + h2 + 2d0*h3))

  bs1t2(2,1)=(1d0/360d0)*(b11 + b12)*(5d0*h1 + 6d0*h2 &
    + h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs1t2(2,2)=(-(1d0/120d0))*b11*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + h3)

  bs1t2(2,3)=(1d0/360d0)*b12*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + 5d0*h3)

  bs1t2(2,4)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*b12*h2 + b11*(2d0*h1 + 6d0*h2 + h3))

  bs1t2(2,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(-6d0*b12*h2 + b11*(h1 + 2d0*h3))

  bs1t2(2,6)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(-h1 + h3) + b11*(h1 + 2d0*h3))

  bs1t2(3,1)=(1d0/360d0)*(b11 + b12)*(5d0*h1 + h2 &
    + 6d0*h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs1t2(3,2)=(1d0/360d0)*b11*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 5d0*h2 + 6d0*h3)

  bs1t2(3,3)=(1d0/120d0)*b12*(h1 + h2 + 6d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs1t2(3,4)=(1d0/90d0)*(b11*(h1 - h2) - b12*(h1 &
    + 2d0*h2))*(b11*(h1 - h2) + b12*(h1 - h3))

  bs1t2(3,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(h1 + 2d0*h2) - 6d0*b11*h3)

  bs1t2(3,6)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*b11*h3 + b12*(2d0*h1 + h2 + 6d0*h3))

  bs1t2(4,1)=(1d0/90d0)*(b11 + b12)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + 2d0*h2 + h3)

  bs1t2(4,2)=(-(1d0/90d0))*b11*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(2d0*h1 + 6d0*h2 + h3)

  bs1t2(4,3)=(1d0/90d0)*b12*(-2d0*h1 - 2d0*h2 &
    + h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs1t2(4,4)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b11*(h1 - h2) - b12*(2d0*h1 + 3d0*h2 + h3))

  bs1t2(4,5)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b11*(h1 + h2 + h3) + b12*(2d0*h1 + 3d0*h2 + h3))

  bs1t2(4,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*((-b12)*(2d0*h1 + h2) + b11*(h1 + h2 + h3))

  bs1t2(5,1)=(1d0/90d0)*(b11 + b12)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(h1 - 2d0*(h2 + h3))

  bs1t2(5,2)=(-(1d0/90d0))*b11*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + 2d0*h3)

  bs1t2(5,3)=(1d0/90d0)*b12*(h1 + 2d0*h2 + 6d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs1t2(5,4)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b11*(2d0*h2 + h3) + b12*(h1 + 3d0*h2 + 2d0*h3))

  bs1t2(5,5)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b12*(h1 + 3d0*h2 + 2d0*h3)+b11*(h1 + 2d0*h2 + 3d0*h3))

  bs1t2(5,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(h2 + 2d0*h3) + b11*(h1 + 2d0*h2 + 3d0*h3))

  bs1t2(6,1)=(1d0/90d0)*(b11 + b12)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + 2d0*h3)

  bs1t2(6,2)=(1d0/90d0)*b11*(-2d0*h1 + h2 - 2d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs1t2(6,3)=(-(1d0/90d0))*b12*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(2d0*h1 + h2 + 6d0*h3)

  bs1t2(6,4)=(-(2d0/45d0))*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*((-b11)*(2d0*h1 + h3) + b12*(h1 + h2 + h3))

  bs1t2(6,5)=(-(2d0/45d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(h1 + h2 + h3) + b11*(2d0*h1 + h2 + 3d0*h3))

  bs1t2(6,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(-h1 + h3) + b11*(2d0*h1 + h2 + 3d0*h3))

end subroutine bsnq1Term2

subroutine bsnq2Term1(bs2t1,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::bs2t1(6,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs2t1(1,1)=(1d0/180d0)*(b11 + b12)*(b21 + b22)*(39d0*h1**2d0 &
    + 15d0*h1*(h2 + h3) + 7d0*(h2**2d0 + h2*h3 + h3**2d0))

  bs2t1(1,2)=(1d0/180d0)*(b11 + b12)*b21*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs2t1(1,3)=(1d0/180d0)*(b11 + b12)*b22*(9d0*h1**2d0 &
    + h2**2d0 + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs2t1(1,4)=(-(1d0/45d0))*(b11 + b12)*(b21*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b22*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs2t1(1,5)=(1d0/45d0)*(b11 + b12)*(b22*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs2t1(1,6)=(-(1d0/45d0))*(b11 + b12)*(b22*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs2t1(2,1)=(1d0/180d0)*b11*(b21 + b22)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs2t1(2,2)=(1d0/180d0)*b11*b21*(7d0*h1**2d0 + 15d0*h1*h2 &
    + 39d0*h2**2d0 + 7d0*h1*h3 + 15d0*h2*h3 + 7d0*h3**2d0)

  bs2t1(2,3)=(-(1d0/180d0))*b11*b22*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs2t1(2,4)=(-(1d0/45d0))*b11*(b21*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs2t1(2,5)=(1d0/45d0)*b11*(b21*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs2t1(2,6)=(-(1d0/45d0))*b11*(b22*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b21*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs2t1(3,1)=(1d0/180d0)*b12*(b21 + b22)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs2t1(3,2)=(-(1d0/180d0))*b12*b21*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs2t1(3,3)=(1d0/180d0)*b12*b22*(7d0*h1**2d0 + 7d0*h1*h2 &
    + 7d0*h2**2d0 + 15d0*h1*h3 + 15d0*h2*h3 + 39d0*h3**2d0)

  bs2t1(3,4)=(-(1d0/45d0))*b12*(b21*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b22*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs2t1(3,5)=(1d0/45d0)*b12*(b22*(-h1**2d0 - 3d0*h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs2t1(3,6)=(-(1d0/45d0))*b12*(b22*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs2t1(4,1)=(-(1d0/45d0))*(b21 + b22)*(b11*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b12*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs2t1(4,2)=(-(1d0/45d0))*b21*(b11*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs2t1(4,3)=(-(1d0/45d0))*b22*(b11*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b12*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs2t1(4,4)=(2d0/45d0)*(b12*(2d0*b22*(h1**2d0 + 3d0*h1*h2 &
    + 6d0*h2**2d0 + h1*h3 + 3d0*h2*h3 + h3**2d0) &
    + b21*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b11*(b22*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 &
    + 4d0*h2*h3 + h3**2d0) + 2d0*b21*(4d0*h1**2d0 + 4d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(h2 + h3))))

  bs2t1(4,5)=(-(2d0/45d0))*(b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b12*(b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) + 2d0*b22*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3))))

  bs2t1(4,6)=(2d0/45d0)*(b12*(-2d0*b22*(h1 - h3)*(h1 + h2 + h3) &
    + b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))) + b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(7d0*h1**2d0 + 2d0*h1*(h2 + h3) + (h2 + h3)**2d0)))

  bs2t1(5,1)=(1d0/45d0)*(b21 + b22)*(b12*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs2t1(5,2)=(1d0/45d0)*b21*(b11*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) + b12*(h1**2d0 &
    + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 + h1*(6d0*h2 + h3)))

  bs2t1(5,3)=(1d0/45d0)*b22*(b12*(-h1**2d0 - 3d0*h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs2t1(5,4)=(-(2d0/45d0))*(b12*(2d0*b22*(h1**2d0 + 3d0*h1*h2 &
    + 6d0*h2**2d0 + h1*h3 + 3d0*h2*h3 + h3**2d0) + b21*(-h1**2d0 &
    + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 + h3**2d0)) &
    + b11*(-2d0*b21*(h1 - h2)*(h1 + h2 + h3) + b22*(h1**2d0 &
    + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3))))

  bs2t1(5,5)=(2d0/45d0)*(b12*(b21*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b22*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3))) + b11*(b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b21*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))))

  bs2t1(5,6)=(-(2d0/45d0))*(b11*(2d0*b21*(h1**2d0 + h1*h2 &
    + h2**2d0 + 3d0*h1*h3 + 3d0*h2*h3 + 6d0*h3**2d0) &
    + b22*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0)) + b12*(-2d0*b22*(h1 - h3)*(h1 + h2 + h3) &
    + b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))))

  bs2t1(6,1)=(-(1d0/45d0))*(b21 + b22)*(b12*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs2t1(6,2)=(-(1d0/45d0))*b21*(b12*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b11*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs2t1(6,3)=(-(1d0/45d0))*b22*(b12*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs2t1(6,4)=(2d0/45d0)*(b11*(-2d0*b21*(h1 - h2)*(h1 + h2 + h3) &
    + b22*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))) + b12*(-2d0*b22*(h1 - h3)*(h1 + h2 &
    + h3) + b21*(7d0*h1**2d0 + 2d0*h1*(h2 + h3) + (h2 + h3)**2d0)))

  bs2t1(6,5)=(-(2d0/45d0))*(b12*(-2d0*b22*(h1 - h3)*(h1 + h2 &
    + h3) + b21*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0)) + b11*(b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b21*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))))

  bs2t1(6,6)=(2d0/45d0)*(b11*(2d0*b21*(h1**2d0 + h1*h2 + h2**2d0 &
    + 3d0*h1*h3 + 3d0*h2*h3 + 6d0*h3**2d0) + b22*(-h1**2d0 &
    + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 + 9d0*h3**2d0)) &
    + b12*(b21*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0) + 2d0*b22*(4d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 4d0*h3**2d0 + 2d0*h1*(h2 + h3))))

end subroutine bsnq2Term1

subroutine bsnq2Term2(bs2t2,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::bs2t2(6,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs2t2(1,1)=(1d0/120d0)*(b21 + b22)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + h3)

  bs2t2(1,2)=(1d0/360d0)*b21*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + 5d0*h2 + h3)

  bs2t2(1,3)=(1d0/360d0)*b22*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + 5d0*h3)

  bs2t2(1,4)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(2d0*h2 + h3) + b21*(6d0*h1 + 2d0*h2 + h3))

  bs2t2(1,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(2d0*h2 + h3) + b21*(h2 + 2d0*h3))

  bs2t2(1,6)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b21*(h2 + 2d0*h3) + b22*(6d0*h1 + h2 + 2d0*h3))

  bs2t2(2,1)=(1d0/360d0)*(b21 + b22)*(5d0*h1 + 6d0*h2 &
    + h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs2t2(2,2)=(-(1d0/120d0))*b21*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + h3)

  bs2t2(2,3)=(1d0/360d0)*b22*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + 5d0*h3)

  bs2t2(2,4)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*b22*h2 + b21*(2d0*h1 + 6d0*h2 + h3))

  bs2t2(2,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(-6d0*b22*h2 + b21*(h1 + 2d0*h3))

  bs2t2(2,6)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(-h1 + h3) + b21*(h1 + 2d0*h3))

  bs2t2(3,1)=(1d0/360d0)*(b21 + b22)*(5d0*h1 + h2 &
    + 6d0*h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs2t2(3,2)=(1d0/360d0)*b21*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(h1 + 5d0*h2 + 6d0*h3)

  bs2t2(3,3)=(1d0/120d0)*b22*(h1 + h2 + 6d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs2t2(3,4)=(1d0/90d0)*(b21*(h1 - h2) - b22*(h1 &
    + 2d0*h2))*(b11*(h1 - h2) + b12*(h1 - h3))

  bs2t2(3,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(h1 + 2d0*h2) - 6d0*b21*h3)

  bs2t2(3,6)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*b21*h3 + b22*(2d0*h1 + h2 + 6d0*h3))

  bs2t2(4,1)=(1d0/90d0)*(b21 + b22)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + 2d0*h2 + h3)

  bs2t2(4,2)=(-(1d0/90d0))*b21*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(2d0*h1 + 6d0*h2 + h3)

  bs2t2(4,3)=(1d0/90d0)*b22*(-2d0*h1 - 2d0*h2 + h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs2t2(4,4)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b21*(h1 - h2) - b22*(2d0*h1 + 3d0*h2 + h3))

  bs2t2(4,5)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b21*(h1 + h2 + h3) + b22*(2d0*h1 + 3d0*h2 + h3))

  bs2t2(4,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*((-b22)*(2d0*h1 + h2) + b21*(h1 + h2 + h3))

  bs2t2(5,1)=(1d0/90d0)*(b21 + b22)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(h1 - 2d0*(h2 + h3))

  bs2t2(5,2)=(-(1d0/90d0))*b21*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(h1 + 6d0*h2 + 2d0*h3)

  bs2t2(5,3)=(1d0/90d0)*b22*(h1 + 2d0*h2 + 6d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs2t2(5,4)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b21*(2d0*h2 + h3) + b22*(h1 + 3d0*h2 + 2d0*h3))

  bs2t2(5,5)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b22*(h1 + 3d0*h2 + 2d0*h3) + b21*(h1 + 2d0*h2 &
    + 3d0*h3))

  bs2t2(5,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 - h3))*(b22*(h2 &
    + 2d0*h3) + b21*(h1 + 2d0*h2 + 3d0*h3))

  bs2t2(6,1)=(1d0/90d0)*(b21 + b22)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + 2d0*h3)

  bs2t2(6,2)=(1d0/90d0)*b21*(-2d0*h1 + h2 - 2d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs2t2(6,3)=(-(1d0/90d0))*b22*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(2d0*h1 + h2 + 6d0*h3)

  bs2t2(6,4)=(-(2d0/45d0))*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*((-b21)*(2d0*h1 + h3) + b22*(h1 + h2 + h3))

  bs2t2(6,5)=(-(2d0/45d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(h1 + h2 + h3) + b21*(2d0*h1 + h2 + 3d0*h3))

  bs2t2(6,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(-h1 + h3) + b21*(2d0*h1 + h2 + 3d0*h3))

end subroutine bsnq2Term2

subroutine bsnq2Term3(bs2t3,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::bs2t3(6,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs2t3(1,1)=(1d0/120d0)*(b11 + b12)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + h2 + h3)

  bs2t3(1,2)=(1d0/360d0)*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*h1 + 5d0*h2 + h3)

  bs2t3(1,3)=(1d0/360d0)*b12*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*h1 + h2 + 5d0*h3)

  bs2t3(1,4)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(2d0*h2 + h3) + b11*(6d0*h1 + 2d0*h2 + h3))

  bs2t3(1,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(2d0*h2 + h3) + b11*(h2 + 2d0*h3))

  bs2t3(1,6)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b11*(h2 + 2d0*h3) + b12*(6d0*h1 + h2 + 2d0*h3))

  bs2t3(2,1)=(1d0/360d0)*(b11 + b12)*(5d0*h1 + 6d0*h2 &
    + h3)*(b21*(-h1 + h2) + b22*(-h1 + h3))

  bs2t3(2,2)=(-(1d0/120d0))*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + h3)

  bs2t3(2,3)=(1d0/360d0)*b12*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + 5d0*h3)

  bs2t3(2,4)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*b12*h2 + b11*(2d0*h1 + 6d0*h2 + h3))

  bs2t3(2,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(-6d0*b12*h2 + b11*(h1 + 2d0*h3))

  bs2t3(2,6)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(-h1 + h3) + b11*(h1 + 2d0*h3))

  bs2t3(3,1)=(1d0/360d0)*(b11 + b12)*(5d0*h1 + h2 &
    + 6d0*h3)*(b21*(-h1 + h2) + b22*(-h1 + h3))

  bs2t3(3,2)=(1d0/360d0)*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 5d0*h2 + 6d0*h3)

  bs2t3(3,3)=(1d0/120d0)*b12*(h1 + h2 + 6d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs2t3(3,4)=(1d0/90d0)*(b11*(h1 - h2) - b12*(h1 &
    + 2d0*h2))*(b21*(h1 - h2) + b22*(h1 - h3))

  bs2t3(3,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(h1 + 2d0*h2) - 6d0*b11*h3)

  bs2t3(3,6)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*b11*h3 + b12*(2d0*h1 + h2 + 6d0*h3))

  bs2t3(4,1)=(1d0/90d0)*(b11 + b12)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + 2d0*h2 + h3)

  bs2t3(4,2)=(-(1d0/90d0))*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(2d0*h1 + 6d0*h2 + h3)

  bs2t3(4,3)=(1d0/90d0)*b12*(-2d0*h1 - 2d0*h2 + h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs2t3(4,4)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b11*(h1 - h2) - b12*(2d0*h1 + 3d0*h2 + h3))

  bs2t3(4,5)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b11*(h1 + h2 + h3) + b12*(2d0*h1 + 3d0*h2 + h3))

  bs2t3(4,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*((-b12)*(2d0*h1 + h2) + b11*(h1 + h2 + h3))

  bs2t3(5,1)=(1d0/90d0)*(b11 + b12)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(h1 - 2d0*(h2 + h3))

  bs2t3(5,2)=(-(1d0/90d0))*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + 2d0*h3)

  bs2t3(5,3)=(1d0/90d0)*b12*(h1 + 2d0*h2 + 6d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs2t3(5,4)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b11*(2d0*h2 + h3) + b12*(h1 + 3d0*h2 + 2d0*h3))

  bs2t3(5,5)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b12*(h1 + 3d0*h2 + 2d0*h3) + b11*(h1 + 2d0*h2 &
    + 3d0*h3))

  bs2t3(5,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(h2 + 2d0*h3) + b11*(h1 + 2d0*h2 + 3d0*h3))

  bs2t3(6,1)=(1d0/90d0)*(b11 + b12)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + h2 + 2d0*h3)

  bs2t3(6,2)=(1d0/90d0)*b11*(-2d0*h1 + h2 - 2d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs2t3(6,3)=(-(1d0/90d0))*b12*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(2d0*h1 + h2 + 6d0*h3)

  bs2t3(6,4)=(-(2d0/45d0))*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*((-b11)*(2d0*h1 + h3) + b12*(h1 + h2 + h3))

  bs2t3(6,5)=(-(2d0/45d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(h1 + h2 + h3) + b11*(2d0*h1 + h2 + 3d0*h3))

  bs2t3(6,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(-h1 + h3) + b11*(2d0*h1 + h2 + 3d0*h3))

end subroutine bsnq2Term3

subroutine bsnq3Term1(bs3t1,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::bs3t1(6,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3

  bs3t1(1,1)=(1d0/180d0)*(b11 + b12)*(b21 + b22)*(39d0*h1**2d0 &
    + 15d0*h1*(h2 + h3) + 7d0*(h2**2d0 + h2*h3 + h3**2d0))

  bs3t1(1,2)=(1d0/180d0)*b11*(b21 + b22)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs3t1(1,3)=(1d0/180d0)*b12*(b21 + b22)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs3t1(1,4)=(-(1d0/45d0))*(b21 + b22)*(b11*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b12*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs3t1(1,5)=(1d0/45d0)*(b21 + b22)*(b12*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 &
    + h3)) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs3t1(1,6)=(-(1d0/45d0))*(b21 + b22)*(b12*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs3t1(2,1)=(1d0/180d0)*(b11 + b12)*b21*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs3t1(2,2)=(1d0/180d0)*b11*b21*(7d0*h1**2d0 + 15d0*h1*h2 &
    + 39d0*h2**2d0 + 7d0*h1*h3 + 15d0*h2*h3 + 7d0*h3**2d0)

  bs3t1(2,3)=(-(1d0/180d0))*b12*b21*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs3t1(2,4)=(-(1d0/45d0))*b21*(b11*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs3t1(2,5)=(1d0/45d0)*b21*(b11*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs3t1(2,6)=(-(1d0/45d0))*b21*(b12*(h1 - h3)*(2d0*h1 &
    - h2 + 2d0*h3) + b11*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 &
    - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0))

  bs3t1(3,1)=(1d0/180d0)*(b11 + b12)*b22*(9d0*h1**2d0 &
    + h2**2d0 + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs3t1(3,2)=(-(1d0/180d0))*b11*b22*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs3t1(3,3)=(1d0/180d0)*b12*b22*(7d0*h1**2d0 + 7d0*h1*h2 &
    + 7d0*h2**2d0 + 15d0*h1*h3 + 15d0*h2*h3 + 39d0*h3**2d0)

  bs3t1(3,4)=(-(1d0/45d0))*b22*(b11*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b12*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs3t1(3,5)=(1d0/45d0)*b22*(b12*(-h1**2d0 - 3d0*h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs3t1(3,6)=(-(1d0/45d0))*b22*(b12*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs3t1(4,1)=(-(1d0/45d0))*(b11 + b12)*(b21*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b22*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs3t1(4,2)=(-(1d0/45d0))*b11*(b21*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs3t1(4,3)=(-(1d0/45d0))*b12*(b21*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b22*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs3t1(4,4)=(2d0/45d0)*(b12*(2d0*b22*(h1**2d0 + 3d0*h1*h2 &
    + 6d0*h2**2d0 + h1*h3 + 3d0*h2*h3 + h3**2d0) &
    + b21*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b11*(b22*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 &
    + 4d0*h2*h3 + h3**2d0) + 2d0*b21*(4d0*h1**2d0 + 4d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(h2 + h3))))

  bs3t1(4,5)=(-(2d0/45d0))*(b12*(2d0*b22*(h1**2d0 + 3d0*h1*h2 &
    + 6d0*h2**2d0 + h1*h3 + 3d0*h2*h3 + h3**2d0) &
    + b21*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b11*(-2d0*b21*(h1 - h2)*(h1 + h2 + h3) &
    + b22*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))))

  bs3t1(4,6)=(2d0/45d0)*(b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + h3))) + b12*(-2d0*b22*(h1 &
    - h3)*(h1 + h2 + h3) + b21*(7d0*h1**2d0 + 2d0*h1*(h2 + h3) &
    + (h2 + h3)**2d0)))

  bs3t1(5,1)=(1d0/45d0)*(b11 + b12)*(b22*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs3t1(5,2)=(1d0/45d0)*b11*(b21*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs3t1(5,3)=(1d0/45d0)*b12*(b22*(-h1**2d0 - 3d0*h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs3t1(5,4)=(-(2d0/45d0))*(b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b12*(b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) + 2d0*b22*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3))))

  bs3t1(5,5)=(2d0/45d0)*(b12*(b21*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b22*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3))) + b11*(b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b21*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))))

  bs3t1(5,6)=(-(2d0/45d0))*(b12*(-2d0*b22*(h1 - h3)*(h1 + h2 &
    + h3) + b21*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0)) + b11*(b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b21*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))))

  bs3t1(6,1)=(-(1d0/45d0))*(b11 + b12)*(b22*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs3t1(6,2)=(-(1d0/45d0))*b11*(b22*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b21*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs3t1(6,3)=(-(1d0/45d0))*b12*(b22*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs3t1(6,4)=(2d0/45d0)*(b12*(-2d0*b22*(h1 - h3)*(h1 + h2 &
    + h3) + b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))) + b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(7d0*h1**2d0 + 2d0*h1*(h2 + h3) &
    + (h2 + h3)**2d0)))

  bs3t1(6,5)=(-(2d0/45d0))*(b11*(2d0*b21*(h1**2d0 + h1*h2 &
    + h2**2d0 + 3d0*h1*h3 + 3d0*h2*h3 + 6d0*h3**2d0) &
    + b22*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0)) + b12*(-2d0*b22*(h1 - h3)*(h1 + h2 + h3) &
    + b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))))

  bs3t1(6,6)=(2d0/45d0)*(b11*(2d0*b21*(h1**2d0 + h1*h2 + h2**2d0 &
    + 3d0*h1*h3 + 3d0*h2*h3 + 6d0*h3**2d0) + b22*(-h1**2d0 &
    + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 + 9d0*h3**2d0)) &
    + b12*(b21*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0) + 2d0*b22*(4d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 4d0*h3**2d0 + 2d0*h1*(h2 + h3))))

end subroutine bsnq3Term1

subroutine bsnq4Term1(bs4t1,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::bs4t1(6,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs4t1(1,1)=(1d0/180d0)*(b21 + b22)**2d0*(39d0*h1**2d0 &
    + 15d0*h1*(h2 + h3) + 7d0*(h2**2d0 + h2*h3 + h3**2d0))

  bs4t1(1,2)=(1d0/180d0)*b21*(b21 + b22)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs4t1(1,3)=(1d0/180d0)*b22*(b21 + b22)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs4t1(1,4)=(-(1d0/45d0))*(b21 + b22)*(b21*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b22*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs4t1(1,5)=(1d0/45d0)*(b21 + b22)*(b22*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs4t1(1,6)=(-(1d0/45d0))*(b21 + b22)*(b22*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs4t1(2,1)=(1d0/180d0)*b21*(b21 + b22)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs4t1(2,2)=(1d0/180d0)*b21**2d0*(7d0*h1**2d0 + 15d0*h1*h2 &
    + 39d0*h2**2d0 + 7d0*h1*h3 + 15d0*h2*h3 + 7d0*h3**2d0)

  bs4t1(2,3)=(-(1d0/180d0))*b21*b22*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs4t1(2,4)=(-(1d0/45d0))*b21*(b21*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs4t1(2,5)=(1d0/45d0)*b21*(b21*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs4t1(2,6)=(-(1d0/45d0))*b21*(b22*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b21*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs4t1(3,1)=(1d0/180d0)*b22*(b21 + b22)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs4t1(3,2)=(-(1d0/180d0))*b21*b22*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs4t1(3,3)=(1d0/180d0)*b22**2d0*(7d0*(h1**2d0 + h1*h2 &
    + h2**2d0) + 15d0*(h1 + h2)*h3 + 39d0*h3**2d0)

  bs4t1(3,4)=(-(1d0/45d0))*b22*(b21*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b22*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs4t1(3,5)=(-(1d0/45d0))*b22*(b22*(h1**2d0 + 2d0*h1*h2 &
    + 3d0*h2**2d0 - h1*h3 - 2d0*h2*h3 - 3d0*h3**2d0) &
    - b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs4t1(3,6)=(-(1d0/45d0))*b22*(b22*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs4t1(4,1)=(-(1d0/45d0))*(b21 + b22)*(b21*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b22*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs4t1(4,2)=(-(1d0/45d0))*b21*(b21*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs4t1(4,3)=(-(1d0/45d0))*b22*(b21*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b22*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs4t1(4,4)=(4d0/45d0)*(b21*b22*(-h1**2d0 + 2d0*h1*h2 &
    + 9d0*h2**2d0 + 4d0*h2*h3 + h3**2d0) &
    + b21**2d0*(4d0*h1**2d0 + 4d0*h2**2d0 + 2d0*h2*h3 &
    + h3**2d0 + 2d0*h1*(h2 + h3)) + b22**2d0*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3)))

  bs4t1(4,5)=(-(4d0/45d0))*((-b21**2d0)*(h1 - h2)*(h1 + h2 &
    + h3) + b21*b22*(2d0*h1*h2 + 6d0*h2**2d0 + h1*h3 &
    + 4d0*h2*h3 + 2d0*h3**2d0) + b22**2d0*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3)))

  bs4t1(4,6)=(4d0/45d0)*((-b21**2d0)*(h1 - h2)*(h1 + h2 &
    + h3) - b22**2d0*(h1 - h3)*(h1 + h2 + h3) &
    + b21*b22*(4d0*h1**2d0 + 2d0*h2**2d0 + 3d0*h2*h3 &
    + 2d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  bs4t1(5,1)=(1d0/45d0)*(b21 + b22)*(b22*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs4t1(5,2)=(1d0/45d0)*b21*(b21*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs4t1(5,3)=(-(1d0/45d0))*b22*(b22*(h1**2d0 + 2d0*h1*h2 &
    + 3d0*h2**2d0 - h1*h3 - 2d0*h2*h3 - 3d0*h3**2d0) &
    - b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs4t1(5,4)=(-(4d0/45d0))*((-b21**2d0)*(h1 - h2)*(h1 + h2 &
    + h3) + b21*b22*(2d0*h1*h2 + 6d0*h2**2d0 + h1*h3 &
    + 4d0*h2*h3 + 2d0*h3**2d0) + b22**2d0*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3)))

  bs4t1(5,5)=(4d0/45d0)*(b21*b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + b22**2d0*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3)) + b21**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs4t1(5,6)=(-(4d0/45d0))*((-b22**2d0)*(h1 - h3)*(h1 + h2 &
    + h3) + b21*b22*(h1*h2 + 2d0*h2**2d0 + 2d0*h1*h3 &
    + 4d0*h2*h3 + 6d0*h3**2d0) + b21**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs4t1(6,1)=(-(1d0/45d0))*(b21 + b22)*(b22*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs4t1(6,2)=(-(1d0/45d0))*b21*(b22*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b21*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs4t1(6,3)=(-(1d0/45d0))*b22*(b22*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs4t1(6,4)=(4d0/45d0)*((-b21**2d0)*(h1 - h2)*(h1 + h2 &
    + h3) - b22**2d0*(h1 - h3)*(h1 + h2 + h3) &
    + b21*b22*(4d0*h1**2d0 + 2d0*h2**2d0 + 3d0*h2*h3 &
    + 2d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  bs4t1(6,5)=(-(4d0/45d0))*((-b22**2d0)*(h1 - h3)*(h1 + h2 &
    + h3) + b21*b22*(h1*h2 + 2d0*h2**2d0 + 2d0*h1*h3 &
    + 4d0*h2*h3 + 6d0*h3**2d0) + b21**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs4t1(6,6)=(4d0/45d0)*(b21*b22*(-h1**2d0 + h2**2d0 &
    + 2d0*h1*h3 + 4d0*h2*h3 + 9d0*h3**2d0) &
    + b22**2d0*(4d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 4d0*h3**2d0 &
    + 2d0*h1*(h2 + h3)) + b21**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

end subroutine bsnq4Term1

subroutine bsnq4Term2(bs4t2,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::bs4t2(6,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs4t2(1,1)=(1d0/120d0)*(b21 + b22)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + h2 + h3)

  bs4t2(1,2)=(1d0/360d0)*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*h1 + 5d0*h2 + h3)

  bs4t2(1,3)=(1d0/360d0)*b22*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*h1 + h2 + 5d0*h3)

  bs4t2(1,4)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(2d0*h2 + h3) + b21*(6d0*h1 + 2d0*h2 + h3))

  bs4t2(1,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(2d0*h2 + h3) + b21*(h2 + 2d0*h3))

  bs4t2(1,6)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b21*(h2 + 2d0*h3) + b22*(6d0*h1 + h2 + 2d0*h3))

  bs4t2(2,1)=(1d0/360d0)*(b21 + b22)*(5d0*h1 + 6d0*h2 &
    + h3)*(b21*(-h1 + h2) + b22*(-h1 + h3))

  bs4t2(2,2)=(-(1d0/120d0))*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + h3)

  bs4t2(2,3)=(1d0/360d0)*b22*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + 5d0*h3)

  bs4t2(2,4)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*b22*h2 + b21*(2d0*h1 + 6d0*h2 + h3))

  bs4t2(2,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(-6d0*b22*h2 + b21*(h1 + 2d0*h3))

  bs4t2(2,6)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(-h1 + h3) + b21*(h1 + 2d0*h3))

  bs4t2(3,1)=(1d0/360d0)*(b21 + b22)*(5d0*h1 + h2 &
    + 6d0*h3)*(b21*(-h1 + h2) + b22*(-h1 + h3))

  bs4t2(3,2)=(1d0/360d0)*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 5d0*h2 + 6d0*h3)

  bs4t2(3,3)=(1d0/120d0)*b22*(h1 + h2 + 6d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs4t2(3,4)=(1d0/90d0)*(b21*(h1 - h2) - b22*(h1 &
    + 2d0*h2))*(b21*(h1 - h2) + b22*(h1 - h3))

  bs4t2(3,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(h1 + 2d0*h2) - 6d0*b21*h3)

  bs4t2(3,6)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*b21*h3 + b22*(2d0*h1 + h2 + 6d0*h3))

  bs4t2(4,1)=(1d0/90d0)*(b21 + b22)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + 2d0*h2 + h3)

  bs4t2(4,2)=(-(1d0/90d0))*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(2d0*h1 + 6d0*h2 + h3)

  bs4t2(4,3)=(1d0/90d0)*b22*(-2d0*h1 - 2d0*h2 + h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs4t2(4,4)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b21*(h1 - h2) - b22*(2d0*h1 + 3d0*h2 + h3))

  bs4t2(4,5)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b21*(h1 + h2 + h3) + b22*(2d0*h1 + 3d0*h2 + h3))

  bs4t2(4,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*((-b22)*(2d0*h1 + h2) + b21*(h1 + h2 + h3))

  bs4t2(5,1)=(1d0/90d0)*(b21 + b22)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(h1 - 2d0*(h2 + h3))

  bs4t2(5,2)=(-(1d0/90d0))*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + 2d0*h3)

  bs4t2(5,3)=(1d0/90d0)*b22*(h1 + 2d0*h2 + 6d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs4t2(5,4)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b21*(2d0*h2 + h3) + b22*(h1 + 3d0*h2 + 2d0*h3))

  bs4t2(5,5)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b22*(h1 + 3d0*h2 + 2d0*h3) + b21*(h1 &
    + 2d0*h2 + 3d0*h3))

  bs4t2(5,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(h2 + 2d0*h3) + b21*(h1 + 2d0*h2 + 3d0*h3))

  bs4t2(6,1)=(1d0/90d0)*(b21 + b22)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + h2 + 2d0*h3)

  bs4t2(6,2)=(1d0/90d0)*b21*(-2d0*h1 + h2 - 2d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs4t2(6,3)=(-(1d0/90d0))*b22*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(2d0*h1 + h2 + 6d0*h3)

  bs4t2(6,4)=(-(2d0/45d0))*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*((-b21)*(2d0*h1 + h3) + b22*(h1 + h2 + h3))

  bs4t2(6,5)=(-(2d0/45d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(h1 + h2 + h3) + b21*(2d0*h1 + h2 + 3d0*h3))

  bs4t2(6,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(-h1 + h3) + b21*(2d0*h1 + h2 + 3d0*h3))

end subroutine bsnq4Term2

subroutine cFluxMat(cxFlux,cyFlux,b11,b12,b21,b22)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::cxFlux(3,6),cyFlux(3,6)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22

  ! !As per paper dN3idx N6j
  ! cxFlux=0d0
  ! cxFlux(1,4)=(1d0/6d0)*(-b11 - b12)
  ! cxFlux(1,5)=(1d0/6d0)*(-b11 - b12)
  ! cxFlux(1,6)=(1d0/6d0)*(-b11 - b12)
  ! cxFlux(2,4)=b11/6d0
  ! cxFlux(2,5)=b11/6d0
  ! cxFlux(2,6)=b11/6d0
  ! cxFlux(3,4)=b12/6d0
  ! cxFlux(3,5)=b12/6d0
  ! cxFlux(3,6)=b12/6d0

  ! !As per paper dN3idy N6j
  ! cyFlux=0d0
  ! cyFlux(1,4)=(1d0/6d0)*(-b21 - b22)
  ! cyFlux(1,5)=(1d0/6d0)*(-b21 - b22)
  ! cyFlux(1,6)=(1d0/6d0)*(-b21 - b22)
  ! cyFlux(2,4)=b21/6d0
  ! cyFlux(2,5)=b21/6d0
  ! cyFlux(2,6)=b21/6d0
  ! cyFlux(3,4)=b22/6d0
  ! cyFlux(3,5)=b22/6d0
  ! cyFlux(3,6)=b22/6d0

  ! As per Shagun N3i dN6jdx
  cxFlux(1,1)=(1d0/6d0)*(-b11 - b12)
  cxFlux(1,2)=0d0
  cxFlux(1,3)=0d0
  cxFlux(1,4)=(b11 - b12)/6d0 
  cxFlux(1,5)=(b11 + b12)/6d0
  cxFlux(1,6)=(1d0/6d0)*(-b11 + b12)

  cxFlux(2,1)=0d0 
  cxFlux(2,2)=b11/6d0 
  cxFlux(2,3)=0d0
  cxFlux(2,4)=(1d0/6d0)*(-b11 - 2d0*b12)
  cxFlux(2,5)=(1d0/6d0)*(b11 + 2d0*b12) 
  cxFlux(2,6)=-(b11/6d0) 

  cxFlux(3,1)=0d0 
  cxFlux(3,2)=0d0 
  cxFlux(3,3)=b12/6d0 
  cxFlux(3,4)=-(b12/6d0)
  cxFlux(3,5)=(1d0/6d0)*(2d0*b11 + b12)
  cxFlux(3,6)=(1d0/6d0)*(-2d0*b11 - b12)

  ! As per Shagun N3i dN6jdy
  cyFlux(1,1)=(1d0/6d0)*(-b21 - b22)
  cyFlux(1,2)=0d0
  cyFlux(1,3)=0d0
  cyFlux(1,4)=(b21 - b22)/6d0 
  cyFlux(1,5)=(b21 + b22)/6d0
  cyFlux(1,6)=(1d0/6d0)*(-b21 + b22)

  cyFlux(2,1)=0d0
  cyFlux(2,2)=b21/6d0
  cyFlux(2,3)=0d0
  cyFlux(2,4)=(1d0/6d0)*(-b21 - 2d0*b22) 
  cyFlux(2,5)=(1d0/6d0)*(b21 + 2d0*b22)
  cyFlux(2,6)=-(b21/6d0)

  cyFlux(3,1)=0d0
  cyFlux(3,2)=0d0
  cyFlux(3,3)=b22/6d0 
  cyFlux(3,4)=-(b22/6d0)
  cyFlux(3,5)=(1d0/6d0)*(2d0*b21 + b22)
  cyFlux(3,6)=(1d0/6d0)*(-2d0*b21 - b22)

end subroutine cFluxMat

subroutine remainMat(dMat,b11,b12,b21,b22,h1,h2,h3)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::dMat(3,3)
  real(kind=C_K2),intent(in)::b11,b12,b21,b22,h1,h2,h3   


  !Taking aux var equation as
  !w = d/dy(h d(eta)/dx) + d/dx(h d(eta)/dy)
  ! dMat(1,1)=(1d0/3d0)*(b11 + b12)*(b21 + b22)*(h1 + h2 + h3)

  ! dMat(1,2)=(-(1d0/6d0))*(b12*b21 + b11*(2d0*b21 &
  !   + b22))*(h1 + h2 + h3)

  ! dMat(1,3)=(-(1d0/6d0))*(b11*b22 + b12*(b21 &
  !   + 2d0*b22))*(h1 + h2 + h3)

  ! dMat(2,1)=(-(1d0/6d0))*(b12*b21 + b11*(2d0*b21 &
  !   + b22))*(h1 + h2 + h3)

  ! dMat(2,2)=(1d0/3d0)*b11*b21*(h1 + h2 + h3)

  ! dMat(2,3)=(1d0/6d0)*(b12*b21 + b11*b22)*(h1 + h2 + h3)

  ! dMat(3,1)=(-(1d0/6d0))*(b11*b22 + b12*(b21 &
  !   + 2d0*b22))*(h1 + h2 + h3)

  ! dMat(3,2)=(1d0/6d0)*(b12*b21 + b11*b22)*(h1 + h2 + h3)

  ! dMat(3,3)=(1d0/3d0)*b12*b22*(h1 + h2 + h3)

  !Taking aux var equation as
  !w = d/dx(h d(eta)/dx) + d/dy(h d(eta)/dy)
  dMat(1,1)=(1d0/6d0)*((b11 + b12)**2d0 + (b21 &
    + b22)**2d0)*(h1 + h2 + h3)

  dMat(1,2)=(-(1d0/6d0))*(b11**2d0 + b11*b12 + b21*(b21 &
    + b22))*(h1 + h2 + h3)

  dMat(1,3)=(-(1d0/6d0))*(b11*b12 + b12**2d0 + b22*(b21 &
    + b22))*(h1 + h2 + h3)

  dMat(2,1)=(-(1d0/6d0))*(b11**2d0 + b11*b12 + b21*(b21 &
    + b22))*(h1 + h2 + h3)

  dMat(2,2)=(1d0/6d0)*(b11**2d0 + b21**2d0)*(h1 + h2 + h3)

  dMat(2,3)=(1d0/6d0)*(b11*b12 + b21*b22)*(h1 + h2 + h3)

  dMat(3,1)=(-(1d0/6d0))*(b11*b12 + b12**2d0 + b22*(b21 &
    + b22))*(h1 + h2 + h3)

  dMat(3,2)=(1d0/6d0)*(b11*b12 + b21*b22)*(h1 + h2 + h3)

  dMat(3,3)=(1d0/6d0)*(b12**2d0 + b22**2d0)*(h1 + h2 + h3)

end subroutine remainMat
