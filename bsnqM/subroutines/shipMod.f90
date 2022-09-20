!!-----------------------------shipMod-----------------------------!!
module shipMod
use bsnqGlobVars
use meshFreeMod
implicit none

  type, public :: shipType    
    integer(kind=C_K1)::posN,posI,totNShip,presFnc
    integer(kind=C_K1)::gridNL,gridNB,gridNN
    real(kind=C_K2)::cl,cb,al,x0,y0,thDeg    
    real(kind=C_K2)::L,B,T
    real(kind=C_K2),allocatable::posData(:,:),posDadtt(:,:)
    real(kind=C_K2),allocatable::gP(:,:) !X major form
    real(kind=C_K2),allocatable::gPPres(:), gPNatCor(:,:)
    integer(kind=C_K1),allocatable::gPFEMele(:)
    logical(kind=C_LG)::dragFlag, initEtaFlag
    type(mfPoiTyp)::gPObj    

    !! gPFEMele is the FEM ele whose nodes neighbours will 
    !! be used by the point-cloud nodes 
  contains
    procedure ::  getPress    
    procedure ::  getLoc
    procedure ::  getLocCS
    procedure ::  calcDrag
    procedure ::  generatePointCloud
    procedure ::  setCubicSpline
  end type shipType  

  interface shipType
    procedure :: initShip
  end interface shipType
  
contains



!!------------------------initShipType-------------------------!!
  type(shipType) function initShip(mf,numShip,simEnd)
  implicit none

    integer(kind=C_K1),intent(in)::mf,numShip
    real(kind=C_K2),intent(in)::simEnd
    integer(kind=C_K1)::i,j,nnMax,k,k2
    real(kind=C_K2)::tmpr1,tmpr2
    character(len=C_KSTR)::bqtxt  

    initShip%totNShip=numShip
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%presFnc
    select case(initShip%presFnc)
      case (1)        
        read(mf,*,end=83,err=83)bqtxt
        read(mf,*,end=83,err=83)initShip%cl,initShip%cb

      case (2) 
        read(mf,*,end=83,err=83)bqtxt
        read(mf,*,end=83,err=83)initShip%cl,initShip%cb,initShip%al
      
      case (3)
        initShip%cl=0
        initShip%cb=0
        initShip%al=0
      
      case DEFAULT
        write(9,*)'[ERR] Invalid ship type'
        stop
    end select

    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%L,initShip%B,initShip%T

    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%initEtaFlag

    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%dragFlag
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)nnMax
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%gridNL,initShip%gridNB

    if(mod(initShip%gridNL,2).eq.0) initShip%gridNL = initShip%gridNL+1
    if(mod(initShip%gridNB,2).eq.0) initShip%gridNB = initShip%gridNB+1

    initShip%gridNN = initShip%gridNL*initShip%gridNB
    allocate( initShip%gP( initShip%gridNN, 2 ) )
    allocate( initShip%gPPres( initShip%gridNN ) )
    allocate( initShip%gPFEMele( initShip%gridNN ) )
    allocate( initShip%gPNatCor( initShip%gridNN, 2 ) )
    call initShip%gPObj%initPoi(nnMax,0d0,0d0,0d0)    

    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%posN        
    read(mf,*,end=83,err=83)bqtxt
    allocate(initShip%posData( initShip%posN, 4))
    allocate(initShip%posDadtt( initShip%posN, 4))
    do j=1,initShip%posN
      read(mf,*,end=83,err=84)initShip%posData(j,1:4)      
    enddo        
    initShip%posDadtt(:,1)=initShip%posData(:,1)
    initShip%posI=1
    call initShip%setCubicSpline

    if(simEnd.gt.initShip%posData( initShip%posN,1 ))then
      write(9,*)"[ERR] Insufficient ship position time information"
      stop
    endif    

    write(9,*)
    write(9,'(" [INF] Ship Read")')
    write(9,'(" [---] NN NL NB (made odd if were given as even)")')
    write(9,'(" [---] ",3I10)') initShip%gridNN, & 
      initShip%gridNL, initShip%gridNB
    write(9,*)

    goto 84
    83 write(9,*) "[ERR] Check ship file format"
    stop
    84 continue

  end function initShip
!!----------------------End initShipType-----------------------!!



!!-----------------------setCubicSpline------------------------!!
  subroutine setCubicSpline(sh)
  implicit none

    class(shipType),intent(inout)::sh    
    integer(kind=C_K1)::i,nn,k    
    real(kind=C_K2),allocatable::mA(:,:),rhs(:,:)

    nn=sh%posN

    ! mA (n x n)
    ! [b1 c1 ----------------------]
    ! [a2 b2 c2 -------------------]
    ! [-- a3 b3 c3 ----------------]
    ! [----- a4 b4 c4 -------------]
    ! [----------------------------]
    ! [---------------------- an bn]
    !
    ! Convert to
    ! 
    ! mA (n x 3)
    ! [ 0 b1 c1]
    ! [a2 b2 c2]
    ! [a3 b3 c3]
    ! [a4 b4 c4]
    ! [--------]
    ! [an bn 0 ]

    ! rhs
    ! corresponding to 
    ! [x, y, th]

    write(9,*)
    write(9,'(" [INF] Ship Cublic Spine Coeffs")')    
    
    allocate(mA(nn,3), rhs(nn,3))
    mA=0d0
    rhs=0d0

    ! natural spline setup
    mA(1,2)=1d0 !b1
    mA(nn,2)=1d0 !bn
    do i=2,nn-1
      mA(i,1) = (sh%posData(i,1)-sh%posData(i-1,1))/6d0
      mA(i,2) = (sh%posData(i+1,1)-sh%posData(i-1,1))/3d0
      mA(i,3) = (sh%posData(i+1,1)-sh%posData(i,1))/6d0

      k=2
      rhs(i,k-1) = (sh%posData(i+1,k) - sh%posData(i,k)) &
        / (sh%posData(i+1,1) - sh%posData(i,1)) - &
        (sh%posData(i,k) - sh%posData(i-1,k)) &
        / (sh%posData(i,1) - sh%posData(i-1,1)) 

      k=3
      rhs(i,k-1) = (sh%posData(i+1,k) - sh%posData(i,k)) &
        / (sh%posData(i+1,1) - sh%posData(i,1)) - &
        (sh%posData(i,k) - sh%posData(i-1,k)) &
        / (sh%posData(i,1) - sh%posData(i-1,1)) 

      k=4
      rhs(i,k-1) = (sh%posData(i+1,k) - sh%posData(i,k)) &
        / (sh%posData(i+1,1) - sh%posData(i,1)) - &
        (sh%posData(i,k) - sh%posData(i-1,k)) &
        / (sh%posData(i,1) - sh%posData(i-1,1)) 
    enddo    

    call triDiagSolver(nn, mA, rhs(:,1), sh%posDadtt(:,2))
    call triDiagSolver(nn, mA, rhs(:,2), sh%posDadtt(:,3))
    call triDiagSolver(nn, mA, rhs(:,3), sh%posDadtt(:,4))
    
    ! write(*,*)sh%posN
    ! do i = 1,nn
    !   !write(*,'(3F15.6," | ",F15.6)')mA(i,1:3), rhs(i,3)
    !   write(*,'(3F15.6," | ",F15.6)')mA(i,1:3), sh%posDadtt(i,4)
    !   !write(*,'(3F15.6," | ",F15.6)')sh%posData(i,1:4)
    ! enddo            

    deallocate(mA, rhs)
    write(9,'(" [---] Done")')    
    write(9,*)
    
  end subroutine setCubicSpline
!!---------------------End setCubicSpline----------------------!!



!!---------------------------getLoc----------------------------!!
  subroutine getLoc(sh,rTime,x0,y0,thDeg)
  implicit none

    class(shipType),intent(inout)::sh    
    integer(kind=C_K1)::i,k,lPosI
    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::x0,y0,thDeg

    if((rTime.ge.sh%posData(sh%posI,1)) .and. &
      (rTime.lt.sh%posData(sh%posI+1,1)))then
      lPosI=sh%posI      
    
    else
      do i=1,sh%posN
        if(sh%posData(i,1).gt.rTime) exit      
      enddo
      lPosI=i-1
      if((lPosI.eq.0) .or. (rTime.gt.sh%posData(sh%posN,1)) )then
        write(9,*)"[ERR] Ship position unavailable at time ",rTime
        stop
      endif
      sh%posI=lPosI
    endif

    k=2
    x0=sh%posData(lPosI,k) &
      +(sh%posData(lPosI+1,k)-sh%posData(lPosI,k)) &
      /(sh%posData(lPosI+1,1)-sh%posData(lPosI,1)) &
      *(rTime-sh%posData(lPosI,1))
    k=3
    y0=sh%posData(lPosI,k) &
      +(sh%posData(lPosI+1,k)-sh%posData(lPosI,k)) &
      /(sh%posData(lPosI+1,1)-sh%posData(lPosI,1)) &
      *(rTime-sh%posData(lPosI,1))
    k=4
    thDeg=sh%posData(lPosI,k) &
      +(sh%posData(lPosI+1,k)-sh%posData(lPosI,k)) &
      /(sh%posData(lPosI+1,1)-sh%posData(lPosI,1)) &
      *(rTime-sh%posData(lPosI,1))

  end subroutine getLoc
!!-------------------------End getLoc--------------------------!!



!!--------------------------getLocCS---------------------------!!
  subroutine getLocCS(sh,rTime,x0,y0,thDeg)
  implicit none

    class(shipType),intent(inout)::sh    
    integer(kind=C_K1)::i,k,lPosI
    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::x0,y0,thDeg

    real(kind=C_K2)::cA, cB, cC, cD, dx

    if((rTime.ge.sh%posData(sh%posI,1)) .and. &
      (rTime.lt.sh%posData(sh%posI+1,1)))then
      lPosI=sh%posI      
    
    else
      do i=1,sh%posN
        if(sh%posData(i,1).gt.rTime) exit      
      enddo
      lPosI=i-1
      if((lPosI.eq.0) .or. (rTime.gt.sh%posData(sh%posN,1)) )then
        write(9,*)"[ERR] Ship position unavailable at time ",rTime
        stop
      endif
      sh%posI=lPosI
    endif

    dx = (sh%posData(lPosI+1,1) - sh%posData(lPosI,1))
    cA = (sh%posData(lPosI+1,1) - rTime) / dx
    cB = 1d0 - cA
    cC = (cA**3 - cA)*dx*dx/6d0
    cD = (cB**3 - cB)*dx*dx/6d0

    k=2
    x0 = (cA * sh%posData(lposI,k)) + (cB * sh%posData(lposI+1,k)) &
      + (cC * sh%posDadtt(lposI,k)) + (cD * sh%posDadtt(lposI+1,k))
  
    k=3
    y0 = (cA * sh%posData(lposI,k)) + (cB * sh%posData(lposI+1,k)) &
      + (cC * sh%posDadtt(lposI,k)) + (cD * sh%posDadtt(lposI+1,k))      

    k=4
    thDeg = (cA * sh%posData(lposI,k)) + (cB * sh%posData(lposI+1,k)) &
      + (cC * sh%posDadtt(lposI,k)) + (cD * sh%posDadtt(lposI+1,k))
  end subroutine getLocCS
!!------------------------End getLocCS-------------------------!!




!!--------------------------getPress---------------------------!!
  subroutine getPress(sh,rTime,npt,cor,pres)
  implicit none

    class(shipType),intent(inout)::sh
    integer(kind=C_K1),intent(in)::npt
    real(kind=C_K2),intent(in)::rTime,cor(npt,2)
    real(kind=C_K2),intent(out)::pres(npt)
    
    integer(kind=C_K1)::i
    real(kind=C_K2)::x0,y0,thDeg,thRad,rRef,cost,sint
    real(kind=C_K2)::x,y,dr2,lx,ly,px,py
    
    call sh%getLocCS(rTime,x0,y0,thDeg)
    thRad=thDeg*deg2rad

    pres=0d0

    select case(sh%presFnc)
    case (1)
        rRef=(max(sh%L,sh%B)/2d0*1.5d0)**2
        cost=dcos(thRad)
        sint=dsin(thRad)

        do i=1,npt
          x=cor(i,1)-x0
          y=cor(i,2)-y0
          dr2=(x**2 + y**2)
          if(dr2.lt.rRef)then
            lx=abs(+x*cost + y*sint)/sh%L
            ly=abs(-x*sint + y*cost)/sh%B
            ! lx=(-x*cost - y*sint)/L
            ! ly=(+x*sint - y*cost)/B
            if(lx.gt.0.5d0)cycle
            if(ly.gt.0.5d0)cycle
            
            if(lx.lt.0.5d0*sh%cl)then
              px=1        
            else
              px=dcos(pi*(lx-0.5d0*sh%cl)/(1d0-sh%cl))      
            endif
            
            if(ly.lt.0.5d0*sh%cb)then
              py=1
            else
              py=dcos(pi*(ly-0.5d0*sh%cb)/(1d0-sh%cb))      
            endif
            
            pres(i)=sh%T*px*px*py*py
          endif
        enddo

      case (2)
        rRef=(max(sh%L,sh%B)/2d0*1.5d0)**2
        cost=dcos(thRad)
        sint=dsin(thRad)

        do i=1,npt
          x=cor(i,1)-x0
          y=cor(i,2)-y0
          dr2=(x**2 + y**2)
          if(dr2.lt.rRef)then
            lx=(+x*cost + y*sint)/sh%L
            ly=(-x*sint + y*cost)/sh%B
            ! lx=(-x*cost - y*sint)/L
            ! ly=(+x*sint - y*cost)/B
            if(abs(lx).gt.0.5d0)cycle
            if(abs(ly).gt.0.5d0)cycle        
            pres(i)=sh%T*(1-sh%cl*(lx**4))*(1-sh%cb*(ly**2)) &
              *exp(-sh%al*(ly**2))
          endif
        enddo

      case (3)        
        do i=1,npt
          x=abs(cor(i,1)-x0)/sh%L
          if(x.gt.0.5d0)cycle
          lx=dcos(pi*x)
          pres(i)=sh%T*lx**2
        enddo

    end select

  end subroutine getPress
!!------------------------End getPress-------------------------!!



!!---------------------generatePointCloud----------------------!!
  subroutine generatePointCloud(sh, rTime)
  implicit none

    class(shipType),intent(inout)::sh
    real(kind=C_K2),intent(in)::rTime

    integer(kind=C_K1)::i, j, k, k2
    real(kind=C_K2)::x0, y0, thDeg, csth, snth, dx, dy
    real(kind=C_K2)::x, y

    call sh%getLocCS(rTime,x0,y0,thDeg)
    csth=dcos(thDeg*deg2rad)
    snth=dsin(thDeg*deg2rad)
    sh%x0 = x0
    sh%y0 = y0
    sh%thDeg = thDeg

    dx = sh%L / (sh%gridNL-1)
    dy = sh%B / (sh%gridNB-1)
    k2=sh%gridNB
    do i=1,sh%gridNL
      do j=1,sh%gridNB
        k=(i-1)*k2+j
        x= -sh%L/2d0 + (i-1)*dx
        y= -sh%B/2d0 + (j-1)*dy
        sh%gP(k,1) = x0 + x*csth - y*snth  
        sh%gP(k,2) = y0 + x*snth + y*csth  
      enddo
    enddo

  end subroutine generatePointCloud
!!-------------------End generatePointCloud--------------------!!



!!--------------------------calcDrag---------------------------!!
  subroutine calcDrag(sh,rTime,nele,np,conn,corx,cory,&
    femPObf,eta,fx,fy,fm,wrkNeiRad)
  implicit none

    class(shipType),intent(inout)::sh
    integer(kind=C_K1),intent(in)::np,nele, conn(nele,6)
    real(kind=C_K2),intent(in)::rTime,corx(np),cory(np)
    real(kind=C_K2),intent(in)::eta(np)
    type(mfPoiTyp),intent(in),target::femPObf(np)
    real(kind=C_K2),intent(out)::fx,fy,fm,wrkNeiRad(np)

    integer(kind=C_K1)::i,j,j2,k,k2,shI,shJ,nn,err
    integer(kind=C_K1)::femNode, femEle
    real(kind=C_K2)::tmpr1,rj
    real(kind=C_K2)::dx,dy,x,y,etaDx,etaDy, lfx, lfy
    type(mfPoiTyp),pointer::femPThis

    
    ! do j=sh%gridNB,1,-1
    !   do i=1,sh%gridNL
    !     k=(i-1)*k2+j
    !     write(*,'(F6.2,",",F6.2," ")',advance='no')sh%gP(k,1), &
    !       sh%gP(k,2)
    !   enddo
    !   write(*,*)
    ! enddo
    ! stop

    dx = sh%L / (sh%gridNL-1)
    dy = sh%B / (sh%gridNB-1)
    
    call sh%getPress(rTime,sh%gridNN,sh%gP,sh%gPPres)    

    fx=0d0
    fy=0d0
    fm=0d0
    k2=sh%gridNB
    do shI=1,sh%gridNL
      do shJ=1,sh%gridNB

        k=(shI-1)*k2+shJ
        sh%gPObj%cx = sh%gP(k,1)
        sh%gPObj%cy = sh%gP(k,2)
        femEle = sh%gPFEMele(k)            
        
        sh%gPObj%rad=0        

        nn=0
        do i = 1,3
          femPThis => femPObf(conn(femEle, i))
          
          do j = 1,femPThis%nn            
            j2 = femPThis%neid(j)
            if( count(sh%gPObj%neid(1:nn).eq.j2).eq.0 )then
              tmpr1 = (corx(j2)-sh%gPObj%cx)**2 &
                + (cory(j2)-sh%gPObj%cy)**2
              rj = femPObf(j2)%rad
              if(tmpr1.gt.(rj*rj))cycle

              nn = nn + 1
              if(nn .gt. sh%gPObj%nnMax)then
                write(9,'("      |",a6,a)')"[ERR]",&
                  "ship%CalcDrag| Increase ship max numOfNegh"
                fx=0d0
                fy=0d0
                fm=0d0
                return
              endif                            

              sh%gPObj%neid(nn)=j2
              wrkNeiRad(nn) = rj
            endif
          enddo
        enddo
        sh%gPObj%nn = nn
                        
        if(nn .lt. 4)then
          write(9,'("      |",a6,a)')"[ERR]",&
            "ship%CalcDrag| Insufficient neghs"
          ! write(*,'(2I5,2F15.6)')shI,shJ,sh%gPObj%cx,sh%gPObj%cy
          ! write(*,'(I5)')nn
          fx=0d0
          fy=0d0
          fm=0d0
          return
        endif
                
        call mls2DDx(sh%gPObj%cx, sh%gPObj%cy, nn, &
          corx( sh%gPObj%neid(1:nn) ), &
          cory( sh%gPObj%neid(1:nn) ), &
          wrkNeiRad( 1:nn ), &
          sh%gPObj%phi(1:nn), sh%gPObj%phiDx(1:nn), &
          sh%gPObj%phiDy(1:nn), err)        

        if(err.ne.0) then !If any err in mls2DDx then return fx=0
          fx=0d0
          fy=0d0
          fm=0d0
          write(9,'("      |",a6,a)')"[ERR]",&
            "ship%CalcDrag| Error in mls2DDx"
          return
        endif

        !! etaDx
        etaDx=0d0
        etaDy=0d0
        do i=1,nn
          j = sh%gPObj%neid(i)
          etaDx=etaDx+ sh%gPObj%phiDx(i) * eta(j)
          etaDy=etaDy+ sh%gPObj%phiDy(i) * eta(j)
        enddo

        !! Simpson's integration
        tmpr1=1d0
        if(mod(shI,2).ne.0)then
          if((shI.eq.1).or.(shI.eq.sh%gridNL))then
            tmpr1=tmpr1*1d0
          else
            tmpr1=tmpr1*2d0
          endif
        else
          tmpr1=tmpr1*4d0
        endif
        if(mod(shJ,2).ne.0)then
          if((shJ.eq.1).or.(shJ.eq.sh%gridNB))then
            tmpr1=tmpr1*1d0
          else
            tmpr1=tmpr1*2d0
          endif
        else
          tmpr1=tmpr1*4d0
        endif
        tmpr1=-tmpr1*dx*dy/9d0 !! correct sign for inward normal to area
        lfx = tmpr1*sh%gPPres(k)*etaDx
        lfy = tmpr1*sh%gPPres(k)*etaDy        
        fx = fx + lfx
        fy = fy + lfy
        fm = fm + (sh%gP(k,1)-sh%x0)*lfy - (sh%gP(k,2)-sh%y0)*lfx
      enddo
    enddo

  end subroutine calcDrag
!!------------------------End calcDrag-------------------------!!



!!------------------------triDiagSolver------------------------!!
  subroutine triDiagSolver(n,A,B,x)
  implicit none

    integer(kind=C_K1),intent(in)::n
    real(kind=C_K2),intent(in)::A(n,3),B(n)
    real(kind=C_K2),intent(out)::x(n)

    integer(kind=C_K1)::i    
    real(kind=C_K2),allocatable::prC(:),prD(:)    

    ! A (n x n)
    ! [b1 c1 ----------------------]
    ! [a2 b2 c2 -------------------]
    ! [-- a3 b3 c3 ----------------]
    ! [----- a4 b4 c4 -------------]
    ! [----------------------------]
    ! [---------------------- an bn]
    !
    ! Converted to following in input itself
    ! 
    ! A (n x 3)
    ! [ 0 b1 c1]
    ! [a2 b2 c2]
    ! [a3 b3 c3]
    ! [a4 b4 c4]
    ! [--------]
    ! [an bn 0 ]

    ! [ai bi ci][xi] = [di]      

    allocate(prC(n), prD(n))

    prC(1) = A(1,3)/A(1,2)
    prD(1) = B(1)/A(1,2)

    do i = 2,n-1
      prC(i) = A(i,3) / ( A(i,2) - A(i,1)*prC(i-1) )      
    enddo

    do i = 2,n      
      prD(i) = ( B(i) - A(i,1)*prD(i-1) ) &
        / ( A(i,2) - A(i,1)*prC(i-1) )
    enddo

    x(n) = prD(n)

    do i = n-1,1,-1
      x(i) = prD(i) - prC(i)*x(i+1)
    enddo

    deallocate(prC, prD)

  end subroutine triDiagSolver
!!----------------------End triDiagSolver----------------------!!


end module shipMod
!!---------------------------End shipMod---------------------------!!