!!---------------------------meshFreeMod---------------------------!!
module meshFreeMod
use bsnqGlobVars
implicit none
  
  type, public :: mfPoiTyp
    !! Random location with FEM mesh neghs
    integer(kind=C_K1)::nn,nnMax,bsnqId
    integer(kind=C_K1),allocatable::neid(:)
    real(kind=C_K2)::cx,cy
    real(kind=C_K2)::rad !Influence radius of the point
    real(kind=C_K2),allocatable::phi(:),phiDx(:),phiDy(:)
    ! bsnqId is set = 0 if the point is not a bsnq point
  contains
    procedure :: initPoi
    procedure :: setPoi
    procedure :: setPoiNZOnly
  end type mfPoiTyp

  type, extends(mfPoiTyp), public :: mfPoiVertVelTyp
    !! Inheriting the base class mfPoiTyp from meshFreeMod
    real(kind=C_K2)::ux,uxx
    real(kind=C_K2)::px,pxx
  !contains
    !procedure :: calcDxx
  end type mfPoiVertVelTyp

contains


!!---------------------------initPoi---------------------------!!
  subroutine initPoi(m,nnMax,cx,cy,rad)
  implicit none

    class(mfPoiTyp),intent(inout)::m
    integer(kind=C_K1),intent(in)::nnMax
    real(kind=C_K2),intent(in)::cx,cy,rad

    if(allocated(m%neid)) deallocate(m%neid)
    if(allocated(m%phi)) deallocate(m%phi)
    if(allocated(m%phiDx)) deallocate(m%phiDx)
    if(allocated(m%phiDy)) deallocate(m%phiDy)
    allocate(m%neid(nnMax),m%phi(nnMax))
    allocate(m%phiDx(nnMax),m%phiDy(nnMax))
    m%nnMax=nnMax
    m%nn=0
    m%cx=cx
    m%cy=cy
    m%rad=rad
    m%bsnqId=0

  end subroutine initPoi
!!-------------------------End initPoi-------------------------!!



!!---------------------------setPoi----------------------------!!
  subroutine setPoi(m,nn,nnMax,bsnqId,cx,cy,rad,&
    nei,phi,phiDx,phiDy)
  implicit none

    class(mfPoiTyp),intent(inout)::m
    integer(kind=C_K1),intent(in)::nn,nnMax,nei(nn),bsnqId
    real(kind=C_K2),intent(in)::cx,cy,rad
    real(kind=C_K2),intent(in)::phi(nn),phiDx(nn),phiDy(nn)

    if(allocated(m%neid)) deallocate(m%neid)
    if(allocated(m%phi)) deallocate(m%phi)
    if(allocated(m%phiDx)) deallocate(m%phiDx)
    if(allocated(m%phiDy)) deallocate(m%phiDy)
    allocate(m%neid(nnMax),m%phi(nnMax))
    allocate(m%phiDx(nnMax),m%phiDy(nnMax))

    if(nn.gt.nnMax)then
      write(9,'(" [ERR] numNei > numNeiMax",2I10)')nn,nnMax
      stop
    endif

    m%nnMax=nnMax
    m%nn=nn
    m%bsnqId=bsnqId
    m%cx=cx
    m%cy=cy
    m%rad=rad
    m%neid(1:nn)=nei
    m%phi(1:nn)=phi
    m%phiDx(1:nn)=phiDx
    m%phiDy(1:nn)=phiDy

  end subroutine setPoi
!!-------------------------End setPoi--------------------------!!



!!---------------------setPoisetPoiNZOnly----------------------!!
  subroutine setPoiNZOnly(m,nn,nnMax,bsnqId,cx,cy,rad,&
    nei,phi,phiDx,phiDy, wrki)
  implicit none

    ! Only stores the neighs with non-zero value for
    ! atleast one of phi, phiDx and phiDy
    ! This minimises the size of the object

    class(mfPoiTyp),intent(inout)::m
    integer(kind=C_K1),intent(in)::nn,nnMax,nei(nn),bsnqId
    real(kind=C_K2),intent(in)::cx,cy,rad
    real(kind=C_K2),intent(in)::phi(nn),phiDx(nn),phiDy(nn)

    integer(kind=C_K1),intent(inout)::wrki

    m%bsnqId=bsnqId
    m%cx=cx
    m%cy=cy
    m%rad=rad
    

    m%nn = 0
    do wrki = 1, nn
      if( abs(phi(wrki)).lt.1d-10 .and. &
        abs(phiDx(wrki)).lt.1d-10 .and. &
        abs(phiDy(wrki)).lt.1d-10 ) cycle

      m%nn = m%nn + 1      
    enddo

    m%nnMax = m%nn

    if(allocated(m%neid)) deallocate(m%neid)
    if(allocated(m%phi)) deallocate(m%phi)
    if(allocated(m%phiDx)) deallocate(m%phiDx)
    if(allocated(m%phiDy)) deallocate(m%phiDy)
    allocate( m%neid(m%nnMax), m%phi(m%nnMax))
    allocate( m%phiDx(m%nnMax), m%phiDy(m%nnMax))


    m%nn = 0
    do wrki = 1, nn
      if( abs(phi(wrki)).lt.1d-10 .and. &
        abs(phiDx(wrki)).lt.1d-10 .and. &
        abs(phiDy(wrki)).lt.1d-10 ) cycle
      
      m%nn = m%nn + 1      
      m%neid(m%nn) = nei(wrki)
      m%phi(m%nn) = phi(wrki)
      m%phiDx(m%nn)=phiDx(wrki)
      m%phiDy(m%nn)=phiDy(wrki)
    enddo       

    if(m%nn.ne.m%nnMax)then
      write(9,'(" [ERR] numNei .ne. numNeiMax in setPoiNZOnly",3I10)') &
        m%bsnqId, m%nn, m%nnMax
      stop
    endif 

  end subroutine setPoiNZOnly
!!-------------------End setPoisetPoiNZOnly--------------------!!



!!---------------------------calcGrad--------------------------!!
  subroutine calcGrad(nn,phiDx,f,fx,i)
  implicit none

    integer(kind=C_K1),intent(in)::nn
    integer(kind=C_K1),intent(out)::i
    real(kind=C_K2),intent(in)::phiDx(nn),f(nn)
    real(kind=C_K2),intent(out)::fx

    fx=0d0
    do i=1,nn
      fx = fx + phiDx(i) * f(i)
    enddo

  end subroutine calcGrad
!!---------------------------calcGrad--------------------------!!



!!---------------------------mls2DDx---------------------------!!
  subroutine mls2DDx(x,y,nn,corx,cory,neiR,phi,phiDx,phiDy,err)
  implicit none

    integer(kind=C_K1),intent(in)::nn        
    real(kind=C_K2),intent(in)::x,y,corx(nn),cory(nn),neiR(nn)
    integer(kind=C_K1),intent(out)::err
    real(kind=C_K2),intent(out)::phi(nn),phiDx(nn),phiDy(nn)

    integer(kind=C_K1)::j
    real(kind=C_K2)::rMax,dr,dx,dy,wj,wjx,wjy,drr,xj,yj,rj
    real(kind=C_K2)::A(3,3),AInv(3,3),ADet
    real(kind=C_K2)::Ax(3,3),Ay(3,3),M(3,3)
    real(kind=C_K2)::pT(3),pTx(3),pTy(3)
    real(kind=C_K2)::Bj(3),Bjx(3),Bjy(3),tmpr1,tmpr2
    real(kind=C_K2)::gT(3),gTx(3),gTy(3),gTemp(3)

    phi=0d0
    phiDx=0d0
    phiDy=0d0
    err=-1
   
    A=0d0    
    Ax=0d0
    Ay=0d0
    pT(1)=1d0
    pT(2)=x
    pT(3)=y      
    pTx(1)=0d0
    pTx(2)=1d0
    pTx(3)=0d0    
    pTy(1)=0d0
    pTy(2)=0d0
    pTy(3)=1d0    

    do j=1,nn
      xj=corx(j)
      yj=cory(j)
      rj=neiR(j)
      dx=xj-x
      dy=yj-y
      dr=dsqrt(dx**2 + dy**2)
      drr=dr/rj
      
      ! drr.gt.1 taken care of in weightFnc
      call weightFnc(dx,dy,dr,drr,rj,wj,wjx,wjy,1)              

      A(1,1)=A(1,1) + wj
      A(1,2)=A(1,2) + wj*xj
      A(1,3)=A(1,3) + wj*yj
      A(2,2)=A(2,2) + wj*xj*xj
      A(2,3)=A(2,3) + wj*xj*yj
      A(3,3)=A(3,3) + wj*yj*yj        

      Ax(1,1)=Ax(1,1) + wjx
      Ax(1,2)=Ax(1,2) + wjx*xj
      Ax(1,3)=Ax(1,3) + wjx*yj
      Ax(2,2)=Ax(2,2) + wjx*xj*xj
      Ax(2,3)=Ax(2,3) + wjx*xj*yj
      Ax(3,3)=Ax(3,3) + wjx*yj*yj        

      Ay(1,1)=Ay(1,1) + wjy
      Ay(1,2)=Ay(1,2) + wjy*xj
      Ay(1,3)=Ay(1,3) + wjy*yj
      Ay(2,2)=Ay(2,2) + wjy*xj*xj
      Ay(2,3)=Ay(2,3) + wjy*xj*yj
      Ay(3,3)=Ay(3,3) + wjy*yj*yj        
    enddo    

    A(2,1)=A(1,2)
    A(3,1)=A(1,3)
    A(3,2)=A(2,3)

    Ax(2,1)=Ax(1,2)
    Ax(3,1)=Ax(1,3)
    Ax(3,2)=Ax(2,3)

    Ay(2,1)=Ay(1,2)
    Ay(3,1)=Ay(1,3)
    Ay(3,2)=Ay(2,3)

    call findInvSymm3x3(A,AInv,ADet)    

    if(abs(ADet).lt.1e-10)then
      write(9,'(" [ERR] ADet is too small",F15.6)')ADet
      write(9,'(" [---] Location X Y",2F15.6)')x,y
      write(9,'(" [---] Num Neigh",I10)')nn
      return
    endif    

    !! gammaTranspose
    call matMul_V13_A33(pT,AInv,gT)
    
    !! gammaTranspose dx
    call matMul_ASym33_BSym33(Ax,AInv,M)
    ! Note : M is not symmetric
    call matMul_V13_A33(gT,M,gTemp)
    call matMul_V13_A33(pTx,AInv,gTx)
    gTx=gTx-gTemp

    !! gammaTranspose dx
    call matMul_ASym33_BSym33(Ay,AInv,M)
    ! Note : M is not symmetric
    call matMul_V13_A33(gT,M,gTemp)
    call matMul_V13_A33(pTy,AInv,gTy)
    gTy=gTy-gTemp

    do j=1,nn
      xj=corx(j)
      yj=cory(j)
      rj=neiR(j)
      dx=xj-x
      dy=yj-y
      dr=dsqrt(dx**2 + dy**2)
      drr=dr/rj

      ! drr.gt.1 taken care of in weightFnc
      call weightFnc(dx,dy,dr,drr,rj,wj,wjx,wjy,1)      

      Bj(1)=wj
      Bj(2)=wj*xj
      Bj(3)=wj*yj
      Bjx(1)=wjx
      Bjx(2)=wjx*xj
      Bjx(3)=wjx*yj
      Bjy(1)=wjy
      Bjy(2)=wjy*xj
      Bjy(3)=wjy*yj

      !! phi
      call matMul_V13_V31(gT,Bj,phi(j))

      !! phiDx
      call matMul_V13_V31(gTx,Bj,tmpr1)
      call matMul_V13_V31(gT,Bjx,tmpr2)
      phiDx(j)=tmpr1+tmpr2

      !! phiDy
      call matMul_V13_V31(gTy,Bj,tmpr1)
      call matMul_V13_V31(gT,Bjy,tmpr2)
      phiDy(j)=tmpr1+tmpr2

      err=0
    enddo    

  end subroutine mls2DDx
!!------------------------End mls2DDx-------------------------!!



!!-------------------------testMls2DDx-------------------------!!
  subroutine testMls2DDx
  implicit none

    integer(kind=C_K1)::nn,i,err
    real(kind=C_K2)::xi,yi,tmpr1,tmpr2,tmpr3
    real(kind=C_K2),allocatable::corx(:),cory(:),nerad(:)
    real(kind=C_K2),allocatable::phi(:),phiDx(:),phiDy(:)
    real(kind=C_K2)::A(3,3),AInv(3,3),ADet


    write(*,*)' \|/ '
    write(*,*)'-- --'
    write(*,*)' /|\ '
    xi=10d0
    yi=89.55d0
    nn=9    

    allocate(corx(nn),cory(nn),nerad(nn))
    allocate(phi(nn),phiDx(nn),phiDy(nn))

    corx = (/ 1d0, 1d0, 0d0, -1d0, -1d0, -1d0, 0d0, 1d0, 0d0 /)
    cory = (/ 0d0, 1d0, 1d0, 1d0, 0d0, -1d0, -1d0, -1d0, 0d0 /)    
    corx=corx+xi
    cory=cory+yi
    nerad=dsqrt(2d0)/0.8d0    

    call mls2DDx(xi, yi, nn, corx, cory, nerad, &
        phi, phiDx, phiDy, err)          

    write(*,'(I10)')nn
    do i=1,nn
      write(*,'(5F15.6)')corx(i),cory(i),phi(i),phiDx(i),phiDx(i)
    enddo

    write(*,'(3F15.6)')sum(phi),sum(phiDx),sum(phiDy)
    tmpr1=0d0
    tmpr2=0d0
    tmpr3=0d0
    do i=1,nn
      tmpr1 = tmpr1 + phi(i)*(1+corx(i))*(10+cory(i))
      tmpr2 = tmpr2 + phiDx(i)*(1+corx(i))
      tmpr3 = tmpr3 + phiDy(i)*(10+cory(i))
    enddo
    write(*,'(3F15.6)')tmpr1,tmpr2,tmpr3


    deallocate(corx,cory,nerad,phi,phiDx,phiDy)
    write(*,*)
    write(*,*)

    write(*,*)'  |  '
    write(*,*)'-- --'
    write(*,*)'  |  '
    xi=10d0
    yi=89.550d0
    nn=5  

    allocate(corx(nn),cory(nn),nerad(nn))
    allocate(phi(nn),phiDx(nn),phiDy(nn))

    corx = (/ 1d0, 0d0, -1d0, 0d0, 0d0 /)
    cory = (/ 0d0, 1d0, 0d0, -1d0, 0d0 /)
    corx=corx+xi
    cory=cory+yi
    nerad=1d0/0.8d0

    call mls2DDx(xi, yi, nn, corx, cory, nerad, &
        phi, phiDx, phiDy, err)          

    write(*,'(I10)')nn
    do i=1,nn
      write(*,'(5F15.6)')corx(i),cory(i),phi(i),phiDx(i),phiDx(i)
    enddo

    write(*,'(3F15.6)')sum(phi),sum(phiDx),sum(phiDy)
    tmpr1=0d0
    tmpr2=0d0
    tmpr3=0d0
    do i=1,nn
      tmpr1 = tmpr1 + phi(i)*(1+corx(i))*(10+cory(i))
      tmpr2 = tmpr2 + phiDx(i)*(1+corx(i))
      tmpr3 = tmpr3 + phiDy(i)*(10+cory(i))
    enddo
    write(*,'(3F15.6)')tmpr1,tmpr2,tmpr3


    deallocate(corx,cory,nerad,phi,phiDx,phiDy)
    write(*,*)
    write(*,*)
    
    write(*,*)'  | /  '
    write(*,*)'  |/   '
    write(*,*)'   ----'
    xi=10d0
    yi=89.55d0
    nn=6    

    allocate(corx(nn),cory(nn),nerad(nn))
    allocate(phi(nn),phiDx(nn),phiDy(nn))

    corx = (/ 1d0, 0d0, 2d0, 0d0, 1d0, 0d0 /)
    cory = (/ 0d0, 1d0, 0d0, 2d0, 1d0, 0d0 /)
    corx=corx+xi
    cory=cory+yi
    nerad=1d0/0.48d0

    call mls2DDx(xi, yi, nn, corx, cory, nerad, &
        phi, phiDx, phiDy, err)          

    write(*,'(I10)')nn
    do i=1,nn
      write(*,'(5F15.6)')corx(i),cory(i),phi(i),phiDx(i),phiDx(i)
    enddo

    write(*,'(3F15.6)')sum(phi),sum(phiDx),sum(phiDy)
    tmpr1=0d0
    tmpr2=0d0
    tmpr3=0d0
    do i=1,nn
      tmpr1 = tmpr1 + phi(i)*(1+corx(i))*(10+cory(i))
      tmpr2 = tmpr2 + phiDx(i)*(1+corx(i))
      tmpr3 = tmpr3 + phiDy(i)*(10+cory(i))
    enddo
    write(*,'(3F15.6)')tmpr1,tmpr2,tmpr3


    deallocate(corx,cory,nerad,phi,phiDx,phiDy)
    write(*,*)
    write(*,*)
    
    write(*,*)'  |  '
    write(*,*)'   --'
    write(*,*)'  |  '
    xi=10d0
    yi=89.55d0
    nn=4    

    allocate(corx(nn),cory(nn),nerad(nn))
    allocate(phi(nn),phiDx(nn),phiDy(nn))

    corx = (/ 1d0, 0d0, 0d0, 0d0 /)
    cory = (/ 0d0, 1d0, -1d0, 0d0 /)
    corx=corx+xi
    cory=cory+yi
    nerad=1d0/0.8d0

    call mls2DDx(xi, yi, nn, corx, cory, nerad, &
        phi, phiDx, phiDy, err)          

    write(*,'(I10)')nn
    do i=1,nn
      write(*,'(5F15.6)')corx(i),cory(i),phi(i),phiDx(i),phiDx(i)
    enddo

    write(*,'(3F15.6)')sum(phi),sum(phiDx),sum(phiDy)
    tmpr1=0d0
    tmpr2=0d0
    tmpr3=0d0
    do i=1,nn
      tmpr1 = tmpr1 + phi(i)*(1+corx(i))*(10+cory(i))
      tmpr2 = tmpr2 + phiDx(i)*(1+corx(i))
      tmpr3 = tmpr3 + phiDy(i)*(10+cory(i))
    enddo
    write(*,'(3F15.6)')tmpr1,tmpr2,tmpr3

    deallocate(corx,cory,nerad,phi,phiDx,phiDy)
    write(*,*)
    write(*,*)


    write(*,*)' \ / '
    write(*,*)'     '
    write(*,*)' / \ '
    xi=10d0
    yi=89.55d0
    nn=5    

    allocate(corx(nn),cory(nn),nerad(nn))
    allocate(phi(nn),phiDx(nn),phiDy(nn))

    corx = (/ 1d0, -0.8d0, -1d0, 1d0, 0d0 /)
    cory = (/ 1d0, 0.8d0, -1d0, -1d0, 0d0 /)    
    corx=corx+xi
    cory=cory+yi
    nerad=dsqrt(2d0)/0.8d0    

    call mls2DDx(xi, yi, nn, corx, cory, nerad, &
        phi, phiDx, phiDy, err)          

    write(*,'(I10)')nn
    do i=1,nn
      write(*,'(5F15.6)')corx(i),cory(i),phi(i),phiDx(i),phiDx(i)
    enddo

    write(*,'(3F15.6)')sum(phi),sum(phiDx),sum(phiDy)
    tmpr1=0d0
    tmpr2=0d0
    tmpr3=0d0
    do i=1,nn
      tmpr1 = tmpr1 + phi(i)*(1+corx(i))*(10+cory(i))
      tmpr2 = tmpr2 + phiDx(i)*(1+corx(i))
      tmpr3 = tmpr3 + phiDy(i)*(10+cory(i))
    enddo
    write(*,'(3F15.6)')tmpr1,tmpr2,tmpr3

    deallocate(corx,cory,nerad,phi,phiDx,phiDy)
    write(*,*)
    write(*,*)
    

    ! A(1,:) = (/ 1d0,0.2d0,3d0 /)
    ! A(2,:) = (/ 0.2d0,4d0,5d0 /)
    ! A(3,:) = (/ 3d0,5d0,6d0 /)

    ! call findInvSymm3x3(A,AInv,ADet)
    ! write(*,*)
    ! write(*,'(3F15.6)')AInv
    ! write(*,'(F15.6)')ADet

  end subroutine testMls2DDx
!!-----------------------End testMls2DDx-----------------------!!



! !!-----------------------mls2DDxSAThesis-----------------------!!
!   subroutine mls2DDxSAThesis(xi,yi,nn,R,corx,cory,phi,phiDx,phiDy)
!   implicit none

!     integer(kind=C_K1),intent(in)::nn    
!     real(kind=C_K2),intent(in)::xi,yi,corx(nn),cory(nn),R
!     real(kind=C_K2),intent(out)::phi(nn),phiDx(nn),phiDy(nn)

!     integer(kind=C_K1)::j,k1,k2,i2
!     real(kind=C_K2)::rMax,dr,dx,dy,wj,wjx,wjy,drr
!     real(kind=C_K2)::m(0:5),mx(0:5),my(0:5)
!     real(kind=C_K2)::c0,c1,c21,c22
!     real(kind=C_K2)::c0x,c1x,c21x,c22x
!     real(kind=C_K2)::c0y,c1y,c21y,c22y

!     !! [Note] : The thesis derivation is probably wrong
!     !! It was taken as the summation instead of integration 
!     !! of the RKPM formulation. But on closer investigation
!     !! it seems MLS is not the same as doing integration 
!     !! instead of summation for RKPM. Anyway MLS is not too difficult.
    
!     !m(0:5) = (/ m0 m1 m2 m11 m12 m22 /)

!     phi=0d0
!     phiDx=0d0
!     phiDy=0d0
   
!     m=0d0
!     mx=0d0
!     my=0d0
!     do j=1,nn
!       dx=corx(j)-xi
!       dy=cory(j)-yi
!       dr=dsqrt(dx**2 + dy**2)
!       drr=dr/R

!       wj=0d0      
!       if(drr.le.1d0)then
!         !! Exp (ending at r/R = 1)
!         ! wj=0.5d0/pi*dexp(-4.5d0*drr*drr)
!         ! wjx=9d0*dx/R/R*wj
!         ! wjy=9d0*dy/R/R*wj

!         !! Biquadratic
!         wj = 1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4
!         wjx= 12d0*dx/(R**2) - 24d0*dx*dr/(R**3) &
!           + 12d0*dx*dr*dr/(R**4)
!         wjy= 12d0*dy/(R**2) - 24d0*dy*dr/(R**3) &
!           + 12d0*dy*dr*dr/(R**4)
!       endif

!       !if(wj.gt.1e-15)then
!         m(0)=m(0)+wj
!         m(1)=m(1)+wj*dx/R
!         m(2)=m(2)+wj*dy/R
!         m(3)=m(3)+wj*dx*dx/R/R
!         m(4)=m(4)+wj*dx*dy/R/R
!         m(5)=m(5)+wj*dy*dy/R/R

!         mx(0)=mx(0)+wjx
!         mx(1)=mx(1)+wjx*dx/R - wj/R
!         mx(2)=mx(2)+wjx*dy/R
!         mx(3)=mx(3)+wjx*dx*dx/R/R - 2d0*wj*dx/R/R
!         mx(4)=mx(4)+wjx*dx*dy/R/R - wj*dy/R/R
!         mx(5)=mx(5)+wjx*dy*dy/R/R

!         my(0)=my(0)+wjy
!         my(1)=my(1)+wjy*dx/R
!         my(2)=my(2)+wjy*dy/R - wj/R
!         my(3)=my(3)+wjy*dx*dx/R/R
!         my(4)=my(4)+wjy*dx*dy/R/R - wj*dx/R/R
!         my(5)=my(5)+wjy*dy*dy/R/R - 2d0*wj*dy/R/R
!       !endif
!     enddo    

!     c0 = m(0)*( m(4)*m(6) - m(5)**2 ) - ( m(1)*m(1)*m(6) &
!       - 2d0*m(1)*m(2)*m(5) + m(2)*m(2)*m(4) )

!     if(abs(c0).lt.1e-10)then
!       write(9,'(" [ERR] C0 is too small",F15.6)')c0
!       write(9,'(" [---] Location X Y",2F15.6)')xi,yi
!       write(9,'(" [---] Num Neigh",I10)')nn
!       stop
!     endif

!     c1 = (m(4)*m(6) - m(5)**2)/c0
!     c21 = (m(2)*m(5) - m(1)*m(6))/c0
!     c22 = (m(1)*m(5) - m(2)*m(4))/c0

!     c0x = m(0)*( mx(4)*m(6) + m(4)*mx(6) - 2d0*m(5)*mx(5) ) &
!       + mx(0)*( m(4)*m(6) - m(5)**2 ) &
!       - ( 2d0*m(1)*mx(1)*m(6) + m(1)*m(1)*mx(6) ) &
!       - ( 2d0*m(2)*mx(2)*m(4) + m(2)*m(2)*mx(4) ) &
!       + ( 2d0*mx(1)*m(2)*m(5) + 2d0*m(1)*mx(2)*m(5) &
!         + 2d0*m(1)*m(2)*mx(5) )
!     c1x = ( mx(4)*m(6) + m(4)*mx(6) - 2d0*mx(5)*m(5) )/c0 &
!       - ( m(4)*m(6) - m(5)**2 )*c0x/c0/c0
!     c21x = ( mx(2)*m(5) + m(2)*mx(5) - mx(1)*m(6) - m(1)*mx(6) )/c0 &
!       - ( m(2)*m(5) - m(1)*m(6) )*c0x/c0/c0
!     c22x = ( mx(1)*m(5) + m(1)*mx(5) - mx(2)*m(4) - m(2)*mx(4) )/c0 &
!       - ( m(1)*m(5) - m(2)*m(4) )*c0x/c0/c0

!     c0y = m(0)*( my(4)*m(6) + m(4)*my(6) - 2d0*m(5)*my(5) ) &
!       + my(0)*( m(4)*m(6) - m(5)**2 ) &
!       - ( 2d0*m(1)*my(1)*m(6) + m(1)*m(1)*my(6) ) &
!       - ( 2d0*m(2)*my(2)*m(4) + m(2)*m(2)*my(4) ) &
!       + ( 2d0*my(1)*m(2)*m(5) + 2d0*m(1)*my(2)*m(5) &
!         + 2d0*m(1)*m(2)*my(5) )
!     c1y = ( my(4)*m(6) + m(4)*my(6) - 2d0*my(5)*m(5) )/c0 &
!       - ( m(4)*m(6) - m(5)**2 )*c0y/c0/c0
!     c21y = ( my(2)*m(5) + m(2)*my(5) - my(1)*m(6) - m(1)*my(6) )/c0 &
!       - ( m(2)*m(5) - m(1)*m(6) )*c0y/c0/c0
!     c22y = ( my(1)*m(5) + m(1)*my(5) - my(2)*m(4) - m(2)*my(4) )/c0 &
!       - ( m(1)*m(5) - m(2)*m(4) )*c0y/c0/c0


!     do j=1,nn
!       dx=corx(j)-xi
!       dy=cory(j)-yi
!       dr=dsqrt(dx**2 + dy**2)
!       drr=dr/R

!       wj=0d0      
!       if(drr.le.1d0)then
!         !! Exp (ending at r/R = 1)
!         ! wj=0.5d0/pi*dexp(-4.5d0*drr*drr)
!         ! wjx=9d0*dx/R/R*wj
!         ! wjy=9d0*dy/R/R*wj

!         !! Biquadratic
!         wj = 1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4
!         wjx= 12d0*dx/(R**2) - 24d0*dx*dr/(R**3) &
!           + 12d0*dx*dr*dr/(R**4)
!         wjy= 12d0*dy/(R**2) - 24d0*dy*dr/(R**3) &
!           + 12d0*dy*dr*dr/(R**4)
!       endif

!       !if(wj.gt.1e-15)then
!         phi(j)=wj*( c1 + c21*dx/R + c22*dy/R )
!         phiDx(j)=wj*( c1x + c21x*dx/R + c22x*dy/R -c21/R ) &
!           + wjx*( c1 + c21*dx/R + c22*dy/R )
!         phiDy(j)=wj*( c1y + c21y*dx/R + c22y*dy/R -c22/R ) &
!           + wjy*( c1 + c21*dx/R + c22*dy/R )
!       !endif
!     enddo    

!   end subroutine mls2DDxSAThesis
! !!---------------------End mls2DDxSAThesis---------------------!!



!!--------------------------weightFnc--------------------------!!
  subroutine weightFnc(dx,dy,dr,drr,R,wj,wjx,wjy,typ)
  implicit none

    integer(kind=C_K1),intent(in)::typ
    real(kind=C_K2),intent(in)::dx,dy,dr,drr,R
    real(kind=C_K2),intent(out)::wj,wjx,wjy

    select case(typ)
      case(1)
        !! Exp (ending at r/R = 1)
        if(drr.gt.1d0)then
          wj=0d0
          wjx=0d0
          wjy=0d0          
        else
          wj=0.5d0/pi*dexp(-4.5d0*drr*drr)
          wjx=9d0*dx/R/R*wj
          wjy=9d0*dy/R/R*wj
        endif

      case(2)
        !! Biquadratic
        if(drr.gt.1d0)then
          wj=0d0
          wjx=0d0
          wjy=0d0          
        else
          wj = 1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4        
          wjx= -12d0*dx/(R**2) + 24d0*dx*dr/(R**3) &
            - 12d0*dx*dr*dr/(R**4)
          wjy= 12d0*dy/(R**2) - 24d0*dy*dr/(R**3) &
            + 12d0*dy*dr*dr/(R**4)    
        endif
    end select

  end subroutine weightFnc
!!------------------------End weightFnc------------------------!!




!!-----------------------findInvSymm3x3------------------------!!
  subroutine findInvSymm3x3(A,AInv,ADet)
  implicit none

    real(kind=C_K2),intent(in)::A(3,3)
    real(kind=C_K2),intent(out)::AInv(3,3),ADet

    AInv(1,1)=-A(2,3)**2 + A(2,2)*A(3,3) 
    AInv(1,2)=A(1,3)*A(2,3) - A(1,2)*A(3,3) 
    AInv(1,3)=-A(1,3)*A(2,2) + A(1,2)*A(2,3) 

    AInv(2,1)=AInv(1,2) 
    AInv(2,2)=-A(1,3)**2 + A(1,1)*A(3,3) 
    AInv(2,3)=A(1,2)*A(1,3) - A(1,1)*A(2,3) 
   
    AInv(3,1)=AInv(1,3)
    AInv(3,2)=AInv(2,3)
    AInv(3,3)=-A(1,2)**2 + A(1,1)*A(2,2)

    ADet = -A(1,3)*A(1,3)*A(2,2) + 2d0*A(1,2)*A(1,3)*A(2,3) &
      - A(1,1)*A(2,3)*A(2,3) - A(1,2)*A(1,2)*A(3,3) &
      + A(1,1)*A(2,2)*A(3,3)

    if(ADet.ne.0d0)then
      AInv = AInv/ADet
    else
      AInv=0d0
    endif

  end subroutine findInvSymm3x3
!!---------------------End findInvSymm3x3----------------------!!



!!---------------------------matMul----------------------------!!
  subroutine matMul_V13_M33_V31(V13,A,V31,res)
  implicit none

    real(kind=C_K2),intent(in)::V13(3),V31(3),A(3,3)
    real(kind=C_K2),intent(out)::res

    res=(A(1,1)*V13(1) + A(2,1)*V13(2) + A(3,1)*V13(3))*V31(1) &
      + (A(1,2)*V13(1) + A(2,2)*V13(2) + A(3,2)*V13(3))*V31(2) &
      + (A(1,3)*V13(1) + A(2,3)*V13(2) + A(3,3)*V13(3))*V31(3)

  end subroutine matMul_V13_M33_V31

  subroutine matMul_ASym33_BSym33(A,B,R)
  implicit none

    real(kind=C_K2),intent(in)::A(3,3),B(3,3)
    real(kind=C_K2),intent(out)::R(3,3)


    R(1,1)=A(1,1)*B(1,1) + A(1,2)*B(1,2) + A(1,3)*B(1,3)
    R(1,2)=A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(2,3)
    R(1,3)=A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)

    R(2,1)=A(1,2)*B(1,1) + A(2,2)*B(1,2) + A(2,3)*B(1,3)
    R(2,2)=A(1,2)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(2,3)
    R(2,3)=A(1,2)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)

    R(3,1)=A(1,3)*B(1,1) + A(2,3)*B(1,2) + A(3,3)*B(1,3)
    R(3,2)=A(1,3)*B(1,2) + A(2,3)*B(2,2) + A(3,3)*B(2,3)
    R(3,3)=A(1,3)*B(1,3) + A(2,3)*B(2,3) + A(3,3)*B(3,3)

  end subroutine matMul_ASym33_BSym33

  subroutine matMul_V13_A33(V,A,res)
  implicit none

    real(kind=C_K2),intent(in)::V(3),A(3,3)
    real(kind=C_K2),intent(out)::res(3)
    
    res(1)=A(1,1)*V(1) + A(2,1)*V(2) + A(3,1)*V(3)
    res(2)=A(1,2)*V(1) + A(2,2)*V(2) + A(3,2)*V(3) 
    res(3)=A(1,3)*V(1) + A(2,3)*V(2) + A(3,3)*V(3)
  end subroutine matMul_V13_A33

  subroutine matMul_V13_V31(V13,V31,res)
  implicit none

    real(kind=C_K2),intent(in)::V13(3),V31(3)
    real(kind=C_K2),intent(out)::res
    
    res = V13(1)*V31(1) + V13(2)*V31(2) + V13(3)*V31(3)
  end subroutine matMul_V13_V31
!!-------------------------End matMul--------------------------!!

end module meshFreeMod
!!-------------------------End meshFreeMod-------------------------!!



