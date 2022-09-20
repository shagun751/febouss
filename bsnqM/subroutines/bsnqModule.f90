include 'mods.f90'
include 'modsInletBC.f90'
include 'modsMFree.f90'
include 'modsVertVel.f90'
include 'shipMod.f90'
include 'bndCondition_v1.0.f90'
include 'bndIntegral.f90'
include 'femAnalyticalTri_v1.1.f90'
include 'findNeiFromLinkList.f90'
include 'geometry_v1.0.f90'
include 'matrixSet1.f90'
include 'matrixSet2.f90'
include 'mergeSort.f90'
include 'nodeConnAll.f90'
include 'rk4Interpolation.f90'
include 'solver_v1.0.f90'
include 'modCellStore.f90'

module bsnqModule
use bsnqGlobVars
use waveFileModule
use waveSourceFunction
use outAbsModule
use shipMod
use rk4InterpMod
use meshFreeMod
use vertVelMod
use cellStoreMod
implicit none

  interface
    subroutine paralution_init(nthreads) BIND(C)
      use, intrinsic :: ISO_C_BINDING, only : C_INT
      integer(kind=C_INT), value, intent(in)  :: nthreads
    end subroutine paralution_init

    

    subroutine paralution_stop() BIND(C)
    end subroutine paralution_stop        

    

    function createLSObj() result(obj) BIND(C, name="createLSObj")
      use, intrinsic :: ISO_C_BINDING, only : C_PTR

      type(C_PTR) :: obj      
    end function createLSObj



    subroutine checkLSObj(obj) BIND(C, name="checkLSObj")
      use, intrinsic :: ISO_C_BINDING, only : C_PTR 

      type(C_PTR), value, intent(in) :: obj
    end subroutine checkLSObj


    subroutine initLSSys(lsObj, np,nnz,iv,jv,gA, &
      atol, rtol, div, maxIter)

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE
      use bsnqGlobVars

      type(C_PTR)::lsObj
      integer(kind=C_K1),intent(in)::np,nnz, maxiter
      integer(kind=C_INT),intent(in),target::iv(np+1),jv(nnz)
      real(kind=C_K2),intent(in)::atol, rtol, div
      real(kind=C_DOUBLE),intent(in),target::gA(nnz)
    end subroutine initLSSys


    subroutine solveSys2(obj,np,nnz,gB,gX,&
      iter,resnorm,ier)

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE
      use bsnqGlobVars

      type(C_PTR),intent(in)::obj
      integer(kind=C_K1),intent(in)::np,nnz  
      integer(kind=C_K1),intent(out)::iter,ier  
      real(kind=C_DOUBLE),intent(in),target::gB(np)  
      real(kind=C_DOUBLE),intent(out),target::gX(np)
      real(kind=C_K2),intent(out)::resnorm  

    end subroutine solveSys2


  end interface
  
  private
  integer(kind=C_K1)::maxNePoi=30
  integer(kind=C_K1)::maxNeEle=10
  !! Number of elements in sysC and sysT
  integer(kind=C_K1),parameter::nSysC=10 

  type, public :: bsnqVars
    integer(kind=C_K1)::npl,npt
    real(kind=C_K2)::rtm
    real(kind=C_K2),allocatable::p(:),q(:),e(:),tD(:)
  contains
    procedure ::  initBsnqVars
  end type bsnqVars

  
  type, public :: bsnqCase
    
    character(len=C_KSTR)::probname,resumeFile
    integer(kind=C_K1)::npl,npq,npt,nele,nbnd,nbndtyp,nedg
    integer(kind=C_K1)::maxNePoi,maxNeEle,nbndp,nthrd,Sz(4)
    integer(kind=C_K1)::maxIter,fileOut,resumeOut
    integer(kind=C_K1)::tStepMethod !(0=RK4, 1=AdBa3, 2=SSPRK3)
    integer(kind=C_K1)::nnzl,nnzf,tStep,nTOb,nSOb
    integer(kind=C_K1),allocatable::conn(:,:),mabnd(:,:)    
    integer(kind=C_K1),allocatable::npoisur(:,:),bndP(:)
    integer(kind=C_K1),allocatable::bndPT(:)
    integer(kind=C_K1),allocatable::wpEle(:)
    integer(kind=C_K1),allocatable::ivl(:),linkl(:),ivq(:),linkq(:)            
    integer(kind=C_INT),allocatable::ivsl(:),jvsl(:)
    integer(kind=C_INT),allocatable::ivsf(:),jvsf(:)
    integer(kind=C_INT),allocatable::ele6x6(:,:),ele6x3(:,:)

    real(kind=C_K2)::dt,errLim,endTime,wvHReset
    real(kind=C_K2)::sysRate,sysT(nSysC)
    real(kind=C_K2)::botFricCd=0d0
    real(kind=C_K2),allocatable::cor(:,:),dep(:)    
    real(kind=C_K2),allocatable::invJ(:,:),bndS(:,:),bndPN(:,:)    
    real(kind=C_K2),allocatable::por(:),presr(:)
    real(kind=C_K2),allocatable::vec6Tmp(:), vec6Tmp2(:)
    real(kind=C_K2),allocatable::ur(:),vr(:),pbpr(:),qbpr(:)    
    real(kind=C_K2),allocatable::uhr(:),vhr(:)
    real(kind=C_K2),allocatable::mass1(:),mass2(:)    
    real(kind=C_K2),allocatable::rowMaxW(:),rowMaxE(:),rowMaxPQ(:)
    real(kind=C_K2),allocatable::gBs5(:),gBs6(:),absC(:)
    real(kind=C_K2),allocatable::gGx(:),gGy(:),gNAdv(:)    
    real(kind=C_K2),allocatable::gPGx(:),gPGy(:)
    real(kind=C_K2),allocatable::gCxF(:),gCyF(:),gDMat(:)
    real(kind=C_K2),allocatable::gBs1(:),gBs2(:),gBs3(:),gBs4(:)
    real(kind=C_K2),allocatable::etaMax(:),etaMin(:),wpLoc(:,:)
    real(kind=C_K2),allocatable::wvAng(:), botFricN6(:)
    real(kind=C_DOUBLE),allocatable::gXW(:),gXE(:),gXPQ(:)
    real(kind=C_DOUBLE),allocatable::gRE(:),gRPQ(:)
    real(kind=C_DOUBLE),allocatable::gMW(:),gME(:),gMPQ(:)
    type(C_PTR)::paralsW, paralsE, paralsPQ

    !! Deallocated in destructR1
    integer(kind=C_K1),allocatable::p2p(:,:),p2e(:,:)
    integer(kind=C_K1),allocatable::ivf(:),linkf(:)  
    real(kind=C_K2),allocatable::massW(:),massE(:)
    real(kind=C_K2),allocatable::gFBs1(:),gFBs2(:),gFBs3(:),gFBs4(:)
    real(kind=C_K2),allocatable::aFull(:),gFW(:)

    integer(kind=C_KCLK)::sysC(nSysC)
    logical(kind=C_LG)::resume, presOn, absOn, botFricOn
    type(wvFileType)::wvF
    type(wvSourceType)::wvS
    type(shipType),allocatable::sh(:)
    type(absTyp),allocatable::absOb(:)
    type(bsnqVars),allocatable::tOb(:),sOb(:)
    !type(bsnqVars)::siOb(:),tiOb !tiOb and siOb For time-interpolation
    !type(rk4InterpTyp)::rk4intp
    type(cellofEleTyp)::cofe
    type(vertVelProbes)::vvPrb, vvMsh    

    !! Optional initialisation
    type(mfPoiTyp),allocatable::pObf(:)
    type(vertVelDerv3D),allocatable::bDf(:), iDf(:)
    

  contains    

    procedure ::  initMat
    procedure ::  meshRead
    procedure ::  femInit
    procedure ::  setRun
    procedure ::  statMatrices
    procedure ::  dynaMatrices
    procedure ::  CSRMatrices
    procedure ::  destructR1    
    procedure ::  outputXML    
    procedure ::  preInstructsRK4
    procedure ::  postInstructs
    procedure ::  solveAll
    procedure ::  diriBCEtaDt
    procedure ::  diriBCPQDt
    procedure ::  diriBCEta
    procedure ::  diriBCPQ
    procedure ::  updateSoln
    procedure ::  writeWaveProbe !and write vertVelProbes
    procedure ::  getEtaPQForXY
    procedure ::  writeResume
    procedure ::  readResume
    procedure ::  caseOutputs
    procedure ::  timeStepRK4
    procedure ::  timeStepSSPRK3
    procedure ::  timeStepEuEx1
    procedure ::  setParalutionLS
    !procedure ::  destructor
    
    procedure ::  setMFree
    procedure ::  calcDepResDerivAll
    procedure ::  calcDepResDeriv
    procedure ::  findEleForLocXY1    !For one location. No OpenMP
    procedure ::  findEleForLocXY2    !For a matrix of locs. OpenMP
    procedure ::  findEleForLocXY3    !For a matrix of locs. OpenMP + Cells
    procedure ::  getVertVel          !using vertVelExp to calculate
    procedure ::  testGetVertVel          
    procedure ::  locWvAng
    !procedure ::  setDepResAtTime     !Set depth resolved quants at time    

  end type bsnqCase

contains

  include 'bsnqModuleFncs.f90'
  include 'bsnqModuleFncs2.f90'
  include 'initialCondition.f90'
  include 'outputXML.f90'

!!-------------------------preInstructsRK4-------------------------!!
  subroutine preInstructsRK4(b)
  implicit none

    class(bsnqCase),intent(inout)::b    
    integer(kind=C_K1)::i

    call system_clock(b%sysC(3)) 
    
    b%tStep=b%tStep+1            

    b%sysT(1)=0d0 ! To time PQ soln in Predictor + Corrector    
    b%sysT(2)=0d0 ! To time solveAll

    do i=b%nTOb-1,1,-1
      b%tOb(i)%rtm = b%tOb(i-1)%rtm
      b%tOb(i)%e = b%tOb(i-1)%e
      b%tOb(i)%p = b%tOb(i-1)%p
      b%tOb(i)%q = b%tOb(i-1)%q
      b%tOb(i)%tD = b%tOb(i-1)%tD
    enddo
    b%tOb(0)%rtm = b%tOb(1)%rtm + b%dt

    write(9,*)
    write(9,'(" ------Time : ",F20.6,"------")') b%tOb(0)%rtm

  end subroutine preInstructsRK4
!!-----------------------End preInstructsRK4-----------------------!!



!!--------------------------postInstructs--------------------------!!
  subroutine postInstructs(b, ltStepMethod)
  implicit none

    class(bsnqCase),intent(inout),target::b    
    integer(kind=C_K1),intent(in)::ltStepMethod
    !(0=RK4), (1=AdBa3), (11=EuEx1)

    type(shipType),pointer::shThis

    integer(kind=C_K1)::i, ishp, j, k
    real(kind=C_K2)::tmpr1,tmpr2,tmpr3, wei(3)

    select case (ltStepMethod)
      case (0) !RK4
        b%tOb(0)%e = b%tOb(1)%e + b%dt/6d0*(b%sOb(1)%e &
          + 2d0*b%sOb(2)%e + 2d0*b%sOb(3)%e + b%sOb(4)%e)
        b%tOb(0)%p = b%tOb(1)%p + b%dt/6d0*(b%sOb(1)%p &
          + 2d0*b%sOb(2)%p + 2d0*b%sOb(3)%p + b%sOb(4)%p)
        b%tOb(0)%q = b%tOb(1)%q + b%dt/6d0*(b%sOb(1)%q &
          + 2d0*b%sOb(2)%q + 2d0*b%sOb(3)%q + b%sOb(4)%q)

      
      case (1) !AdBa3
        b%tOb(0)%e = b%tOb(1)%e + b%dt/12d0*(23d0*b%sOb(1)%e &
          - 16d0*b%sOb(2)%e + 5d0*b%sOb(3)%e)
        b%tOb(0)%p = b%tOb(1)%p + b%dt/12d0*(23d0*b%sOb(1)%p &
          - 16d0*b%sOb(2)%p + 5d0*b%sOb(3)%p)
        b%tOb(0)%q = b%tOb(1)%q + b%dt/12d0*(23d0*b%sOb(1)%q &
          - 16d0*b%sOb(2)%q + 5d0*b%sOb(3)%q)


      case (2) !SSPRK3
        b%tOb(0)%e = b%tOb(1)%e + b%dt/6d0 * ( b%sOb(1)%e &
          + b%sOb(2)%e + 4d0*b%sOb(3)%e )
        b%tOb(0)%p = b%tOb(1)%p + b%dt/6d0 * ( b%sOb(1)%p &
          + b%sOb(2)%p + 4d0*b%sOb(3)%p )
        b%tOb(0)%q = b%tOb(1)%q + b%dt/6d0 * ( b%sOb(1)%q &
          + b%sOb(2)%q + 4d0*b%sOb(3)%q )        


      case (11) !EuEx1
        b%tOb(0)%e = b%tOb(1)%e + b%dt*b%sOb(1)%e
        b%tOb(0)%p = b%tOb(1)%p + b%dt*b%sOb(1)%p
        b%tOb(0)%q = b%tOb(1)%q + b%dt*b%sOb(1)%q

      case DEFAULT
        write(9,'(" [ERR] Wrong time-step method!")')
        stop
    end select

    !! Forcing Dirichlet BC
    call b%diriBCEta(b%tOb(0)%e,b%tOb(0)%rtm)
    call b%diriBCPQ(b%tOb(0)%p,b%tOb(0)%q,b%tOb(0)%rtm)
    
    b%tOb(0)%tD(1:b%npl) = b%dep(1:b%npl) + b%tOb(0)%e
    call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD)

    do i=1,b%npl
      tmpr1=b%tOb(0)%e(i)
      if(b%etaMin(i).gt.tmpr1) b%etaMin(i)=tmpr1
      if(b%etaMax(i).lt.tmpr1) b%etaMax(i)=tmpr1
    enddo

    tmpr1 = mod( b%tOb(0)%rtm, b%wvHReset )
    if( ( tmpr1 .lt. 0.1*b%dt ) .or. &
        ( tmpr1 .gt. b%wvHReset-0.1*b%dt ) )then 
      !the above coz possible mod(11.9999,3) or mod(12.01,3)
      b%etaMin=b%tOb(0)%e
      b%etaMax=b%tOb(0)%e
    endif

    ! do i=1,b%nbndp
    !   j2=b%bndPT(i)
    !   if(j2.eq.11)then
    !     i2=b%bndP(i)
    !     call b%wvIn%getEta(b%tOb(0)%rtm,b%cor(i2,1),b%cor(i2,2),tmpr1)
    !     call b%wvIn%getEta(b%tOb(1)%rtm,b%cor(i2,1),b%cor(i2,2),tmpr2)
    !     write(9,302)"InWv",b%tOb(0)%rtm,b%tOb(0)%e(i2),&
    !       b%tOb(0)%p(i2),b%tOb(0)%q(i2)
    !     write(9,302)"InEt",b%tOb(1)%e(i2),tmpr1,tmpr2,tmpr1-tmpr2
    !     exit
    !   endif
    ! enddo

    ! b%vec6Tmp = b%tOb(0)%tD - b%dep
    ! call b%locWvAng( b%npt, b%vec6Tmp ) !-pi to pi

    if(allocated(b%bDf))then 
      call b%calcDepResDerivAll( b%npt, b%tOb(0)%p, b%tOb(0)%q, &
        b%tOb(0)%tD, b%bDf )
      !call b%testGetVertVel
    endif
        
    !! Ship drag calculation and position reporting
    if(b%presOn)then
      b%vec6Tmp = b%tOb(0)%tD - b%dep
      
      do ishp = 1, b%sh(1)%totNShip
        if(b%sh(ishp)%dragFlag)then !Optional to calc drag          
          
          shThis => b%sh(ishp)
          ! Ship centre location updated in this function
          ! along with the generation of the point cloud
          call shThis%generatePointCloud( b%tOb(0)%rtm )

          write(9,303) 'shP', ishp, b%tOb(0)%rtm, shThis%x0, &
            shThis%y0, shThis%thDeg
          
          ! get the element that each point in pointCloud is
          call b%findEleForLocXY3( shThis%gridNN, &
            shThis%gP(:,1), shThis%gP(:,2), &
            shThis%gPFEMele, shThis%gPNatCor(:,1), & 
            shThis%gPNatCor(:,2) )

          if(minval(shThis%gPFEMele).eq.-1) then
            write(9,'("      |",a6,a)')"[ERR]",&
              "ship%CalcDrag| Ship not fully in domain "
            write(9,303)'shF',ishp,b%tOb(0)%rtm,0d0,0d0,0d0
            cycle
          endif                    

          i=b%npt
          call shThis%calcDrag(b%tOb(0)%rtm, b%nele, i, &
            b%conn, b%cor(1:i,1), b%cor(1:i,2), &
            b%pObf(1:i), b%vec6Tmp(1:i), &
            tmpr1, tmpr2, tmpr3, b%vec6Tmp2(1:i))
          write(9,303)'shF',ishp,b%tOb(0)%rtm,tmpr1,tmpr2,tmpr3
        endif
      enddo
    endif

    !write(201,*)    

    call system_clock(b%sysC(4))
    ! bq%sysT(1) = To time PQ soln in Predictor + Corrector    
    ! bq%sysT(2) = To time solveAll    
    ! bq%sysT(3) = To time the current time-step
    b%sysT(3)=1d0*(b%sysC(4)-b%sysC(3))/b%sysRate    
    write(9,301)"[TIL]", b%sysT(3), b%sysT(1)/b%sysT(3)*100d0, &
      b%sysT(2)/b%sysT(3)*100d0
    !write(9,*)
    301 format('      |',a6,3F15.4)
    302 format('      |',a6,4F15.4)
    303 format('      |',a6,I10,4F15.4)

  end subroutine postInstructs
!!------------------------End postInstructs------------------------!!



!!---------------------------caseOutputs---------------------------!!
  subroutine caseOutputs(b)
  implicit none

    class(bsnqCase),intent(inout)::b    
    real(kind=C_K2)::tmpr1

    call system_clock(b%sysC(7))

    if(mod(b%tStep,b%fileOut).eq.0) then
      call b%getVertVel( b%vvMsh%np, b%vvMsh%x, b%vvMsh%y, &
        b%vvMsh%z, b%vvMsh%u, b%vvMsh%v, b%vvMsh%w, &
        b%vvMsh%p, b%vvMsh%eta, &
        b%vvMsh%wrki, b%vvMsh%wrkr, b%vvMsh%err, tmpr1 )
      call b%outputXML
    endif

    if(mod(b%tStep,b%resumeOut).eq.0) then
      if(b%tStep.ne.0)then
        call b%writeResume
      endif
    endif

    call b%writeWaveProbe    !and write vertVelProbes

    call system_clock(b%sysC(8))

    write(9,301) "[TIO]", &
      1d0*(b%sysC(8) - b%sysC(7))/b%sysRate    

    301 format('      |',a6,F15.4)

  end subroutine caseOutputs
!!-------------------------End caseOutputs-------------------------!!



!!---------------------------timeStepRK4---------------------------!!
  subroutine timeStepRK4(bq)
  implicit none

    class(bsnqCase),intent(inout)::bq
    real(kind=C_K2)::rTime

    !!-------------RK4 S1--------------!!
    bq%ur = bq%tOb(1)%p / bq%tOb(1)%tD
    bq%vr = bq%tOb(1)%q / bq%tOb(1)%tD
    bq%pbpr = bq%tOb(1)%p / bq%por
    bq%qbpr = bq%tOb(1)%q / bq%por   
    rTime=bq%tOb(1)%rtm
    call bq%dynaMatrices(rTime,bq%tOb(1)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(1)%p, bq%tOb(1)%q, &
      bq%pbpr, bq%qbpr, bq%presr, bq%tOb(1)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ, bq%sysC)
    call bq%updateSoln(1)
    !!-----------End RK4 S1------------!!

    !!-------------RK4 S2--------------!!
    bq%ur = bq%tOb(0)%p / bq%tOb(0)%tD
    bq%vr = bq%tOb(0)%q / bq%tOb(0)%tD
    bq%pbpr = bq%tOb(0)%p / bq%por
    bq%qbpr = bq%tOb(0)%q / bq%por   
    rTime=(bq%tOb(0)%rtm + bq%tOb(1)%rtm)/2d0
    call bq%dynaMatrices(rTime,bq%tOb(0)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(0)%p, bq%tOb(0)%q, &
      bq%pbpr, bq%qbpr, bq%presr, bq%tOb(0)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ, bq%sysC)
    call bq%updateSoln(2)
    !!-----------End RK4 S2------------!!

    !!-------------RK4 S3--------------!!
    bq%ur = bq%tOb(0)%p / bq%tOb(0)%tD
    bq%vr = bq%tOb(0)%q / bq%tOb(0)%tD
    bq%pbpr = bq%tOb(0)%p / bq%por
    bq%qbpr = bq%tOb(0)%q / bq%por   
    rTime=(bq%tOb(0)%rtm + bq%tOb(1)%rtm)/2d0
    call bq%dynaMatrices(rTime,bq%tOb(0)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(0)%p, bq%tOb(0)%q, &
      bq%pbpr, bq%qbpr, bq%presr, bq%tOb(0)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ, bq%sysC)
    call bq%updateSoln(3)
    !!-----------End RK4 S3------------!!

    !!-------------RK4 S4--------------!!
    bq%ur = bq%tOb(0)%p / bq%tOb(0)%tD
    bq%vr = bq%tOb(0)%q / bq%tOb(0)%tD
    bq%pbpr = bq%tOb(0)%p / bq%por
    bq%qbpr = bq%tOb(0)%q / bq%por   
    rTime=bq%tOb(0)%rtm
    call bq%dynaMatrices(rTime,bq%tOb(0)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(0)%p, bq%tOb(0)%q, &
      bq%pbpr, bq%qbpr, bq%presr, bq%tOb(0)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ, bq%sysC)
    call bq%updateSoln(4)
    !!-----------End RK4 S4------------!!

  end subroutine timeStepRK4
!!-------------------------End timeStepRK4-------------------------!!



!!----------------------------updateSoln---------------------------!!
  subroutine updateSoln(b,step)
  implicit none

    class(bsnqCase),intent(inout)::b    
    integer(kind=4),intent(in)::step
    integer(kind=4)::i,j

    b%sysT(1)=b%sysT(1)+1d0*(b%sysC(8)-b%sysC(7))/b%sysRate
    b%sysT(2)=b%sysT(2)+1d0*(b%sysC(6)-b%sysC(5))/b%sysRate
    
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j)
    select case(step)

      case (1) !RK4 - Step1
        !$OMP DO SCHEDULE(dynamic,1000)
        do i = 1, b%npl          
          b%sOb(1)%e(i) = b%gXE(i)          
          
          b%tOb(0)%e(i) = b%tOb(1)%e(i) + b%dt * b%gXE(i)/2d0          
          b%tOb(0)%tD(i) = b%dep(i) + b%tOb(0)%e(i)
        enddo
        !$OMP END DO NOWAIT
        
        !$OMP DO SCHEDULE(dynamic,1000)
        do i = 1, b%npt
          j = b%npt + i          
          b%sOb(1)%p(i) = b%gXPQ(i)
          b%sOb(1)%q(i) = b%gXPQ(j)        
                    
          b%tOb(0)%p(i) = b%tOb(1)%p(i) + b%dt * b%gXPQ(i)/2d0
          b%tOb(0)%q(i) = b%tOb(1)%q(i) + b%dt * b%gXPQ(j)/2d0          
        enddo        
        !$OMP END DO NOWAIT
        call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD) 

      case (2) !RK4 - Step2
        !$OMP DO SCHEDULE(dynamic,1000)
        do i = 1, b%npl          
          b%sOb(2)%e(i) = b%gXE(i)          
          
          b%tOb(0)%e(i) = b%tOb(1)%e(i) + b%dt * b%gXE(i)/2d0          
          b%tOb(0)%tD(i) = b%dep(i) + b%tOb(0)%e(i)
        enddo
        !$OMP END DO NOWAIT
        
        !$OMP DO SCHEDULE(dynamic,1000)
        do i = 1, b%npt
          j = b%npt + i          
          b%sOb(2)%p(i) = b%gXPQ(i)
          b%sOb(2)%q(i) = b%gXPQ(j)        
                    
          b%tOb(0)%p(i) = b%tOb(1)%p(i) + b%dt * b%gXPQ(i)/2d0
          b%tOb(0)%q(i) = b%tOb(1)%q(i) + b%dt * b%gXPQ(j)/2d0          
        enddo        
        !$OMP END DO NOWAIT
        call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD) 

      case (3) !RK4 - Step3
        !$OMP DO SCHEDULE(dynamic,1000)
        do i = 1, b%npl          
          b%sOb(3)%e(i) = b%gXE(i)          
          
          b%tOb(0)%e(i) = b%tOb(1)%e(i) + b%dt * b%gXE(i)
          b%tOb(0)%tD(i) = b%dep(i) + b%tOb(0)%e(i)
        enddo
        !$OMP END DO NOWAIT
        
        !$OMP DO SCHEDULE(dynamic,1000)
        do i = 1, b%npt
          j = b%npt + i          
          b%sOb(3)%p(i) = b%gXPQ(i)
          b%sOb(3)%q(i) = b%gXPQ(j)        
                    
          b%tOb(0)%p(i) = b%tOb(1)%p(i) + b%dt * b%gXPQ(i)
          b%tOb(0)%q(i) = b%tOb(1)%q(i) + b%dt * b%gXPQ(j)
        enddo        
        !$OMP END DO NOWAIT
        call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD) 

      case (4) !RK4 - Step4
        !$OMP DO SCHEDULE(dynamic,1000)
        do i = 1, b%npl          
          b%sOb(4)%e(i) = b%gXE(i)                  
        enddo
        !$OMP END DO NOWAIT
        
        !$OMP DO SCHEDULE(dynamic,1000)
        do i = 1, b%npt
          j = b%npt + i          
          b%sOb(4)%p(i) = b%gXPQ(i)
          b%sOb(4)%q(i) = b%gXPQ(j)        
        enddo        
        !$OMP END DO NOWAIT

    end select        
    !$OMP END PARALLEL          

  end subroutine updateSoln
!!--------------------------End updateSoln-------------------------!!



!!--------------------------timeStepEuEx1--------------------------!!
  subroutine timeStepEuEx1(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=4)::i,j
    real(kind=C_K2)::rTime

    !---------preInstruct----------!
    call system_clock(b%sysC(3)) 
    
    b%tStep=b%tStep+1            

    b%sysT(1)=0d0 ! To time PQ soln in Predictor + Corrector    
    b%sysT(2)=0d0 ! To time solveAll

    do i=b%nSOb,2,-1
      b%sOb(i)%rtm = b%sOb(i-1)%rtm
      b%sOb(i)%e = b%sOb(i-1)%e
      b%sOb(i)%p = b%sOb(i-1)%p
      b%sOb(i)%q = b%sOb(i-1)%q
      b%sOb(i)%tD = b%sOb(i-1)%tD
    enddo

    do i=b%nTOb-1,1,-1
      b%tOb(i)%rtm = b%tOb(i-1)%rtm
      b%tOb(i)%e = b%tOb(i-1)%e
      b%tOb(i)%p = b%tOb(i-1)%p
      b%tOb(i)%q = b%tOb(i-1)%q
      b%tOb(i)%tD = b%tOb(i-1)%tD
    enddo
    b%tOb(0)%rtm = b%tOb(1)%rtm + b%dt

    write(9,*)
    write(9,'(" ------Time : ",F20.6,"------")') b%tOb(0)%rtm
    !-------End preInstruct--------!
    
    !-----------solution-----------!
    b%ur = b%tOb(1)%p / b%tOb(1)%tD
    b%vr = b%tOb(1)%q / b%tOb(1)%tD
    b%pbpr = b%tOb(1)%p / b%por
    b%qbpr = b%tOb(1)%q / b%por   
    rTime=b%tOb(1)%rtm
    call b%dynaMatrices(rTime,b%tOb(1)%tD, b%ur, b%vr)
    call b%solveAll(rTime, b%tOb(1)%p, b%tOb(1)%q, &
      b%pbpr, b%qbpr, b%presr, b%tOb(1)%e, &
      b%gXW, b%gXE, b%gXPQ, b%gRE, b%gRPQ, b%sysC)   
    !---------End solution---------!

    !----------updateSoln----------!
    b%sysT(1)=b%sysT(1)+1d0*(b%sysC(8)-b%sysC(7))/b%sysRate
    b%sysT(2)=b%sysT(2)+1d0*(b%sysC(6)-b%sysC(5))/b%sysRate

    b%sOb(1)%e = b%gXE
    b%sOb(1)%p = b%gXPQ(1:b%npt)
    b%sOb(1)%q = b%gXPQ(b%npt+1:2*b%npt)        
    !--------End updateSoln--------!

  end subroutine timeStepEuEx1
!!------------------------End timeStepEuEx1------------------------!!



!!-------------------------timeStepSSPRK3--------------------------!!
  subroutine timeStepSSPRK3(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=4)::i,j
    real(kind=C_K2)::rTime


    !!------------preInstructs-----------!!
    call system_clock(b%sysC(3)) 
    
    b%tStep=b%tStep+1            

    b%sysT(1)=0d0 ! To time PQ soln in Predictor + Corrector    
    b%sysT(2)=0d0 ! To time solveAll

    do i=b%nTOb-1,1,-1
      b%tOb(i)%rtm = b%tOb(i-1)%rtm
      b%tOb(i)%e = b%tOb(i-1)%e
      b%tOb(i)%p = b%tOb(i-1)%p
      b%tOb(i)%q = b%tOb(i-1)%q
      b%tOb(i)%tD = b%tOb(i-1)%tD
    enddo
    b%tOb(0)%rtm = b%tOb(1)%rtm + b%dt

    write(9,*)
    write(9,'(" ------Time : ",F20.6,"------")') b%tOb(0)%rtm
    !!----------End preInstructs---------!!


    !!-------------SSPRK3 S1-------------!!
    b%ur = b%tOb(1)%p / b%tOb(1)%tD
    b%vr = b%tOb(1)%q / b%tOb(1)%tD
    b%pbpr = b%tOb(1)%p / b%por
    b%qbpr = b%tOb(1)%q / b%por   
    rTime = b%tOb(1)%rtm
    call b%dynaMatrices(rTime,b%tOb(1)%tD, b%ur, b%vr)
    call b%solveAll(rTime, b%tOb(1)%p, b%tOb(1)%q, &
      b%pbpr, b%qbpr, b%presr, b%tOb(1)%e, &
      b%gXW, b%gXE, b%gXPQ, b%gRE, b%gRPQ, b%sysC)    
    !!-----------End SSPRK3 S1-----------!!


    !!-----------updateSoln S1-----------!!
    b%sysT(1) = b%sysT(1) + 1d0*(b%sysC(8)-b%sysC(7))/b%sysRate
    b%sysT(2) = b%sysT(2) + 1d0*(b%sysC(6)-b%sysC(5))/b%sysRate

    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j)
    !$OMP DO SCHEDULE(dynamic,1000)
    do i = 1, b%npl          
      b%sOb(1)%e(i) = b%gXE(i)          
      
      b%tOb(0)%e(i) = b%tOb(1)%e(i) + b%dt * b%gXE(i)
      b%tOb(0)%tD(i) = b%dep(i) + b%tOb(0)%e(i)
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO SCHEDULE(dynamic,1000)
    do i = 1, b%npt
      j = b%npt + i          
      b%sOb(1)%p(i) = b%gXPQ(i)
      b%sOb(1)%q(i) = b%gXPQ(j)        
                
      b%tOb(0)%p(i) = b%tOb(1)%p(i) + b%dt * b%gXPQ(i)
      b%tOb(0)%q(i) = b%tOb(1)%q(i) + b%dt * b%gXPQ(j)          
    enddo        
    !$OMP END DO NOWAIT
    call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD) 
    !$OMP END PARALLEL          
    !!---------End updateSoln S1---------!!


    !!-------------SSPRK3 S2-------------!!
    b%ur = b%tOb(0)%p / b%tOb(0)%tD
    b%vr = b%tOb(0)%q / b%tOb(0)%tD
    b%pbpr = b%tOb(0)%p / b%por
    b%qbpr = b%tOb(0)%q / b%por   
    rTime = b%tOb(1)%rtm + b%dt
    call b%dynaMatrices(rTime,b%tOb(0)%tD, b%ur, b%vr)
    call b%solveAll(rTime, b%tOb(0)%p, b%tOb(0)%q, &
      b%pbpr, b%qbpr, b%presr, b%tOb(0)%e, &
      b%gXW, b%gXE, b%gXPQ, b%gRE, b%gRPQ, b%sysC)    
    !!-----------End SSPRK3 S2-----------!!


    !!-----------updateSoln S2-----------!!
    b%sysT(1) = b%sysT(1) + 1d0*(b%sysC(8)-b%sysC(7))/b%sysRate
    b%sysT(2) = b%sysT(2) + 1d0*(b%sysC(6)-b%sysC(5))/b%sysRate

    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j)
    !$OMP DO SCHEDULE(dynamic,1000)
    do i = 1, b%npl          
      b%sOb(2)%e(i) = b%gXE(i)          
      
      b%tOb(0)%e(i) = b%tOb(1)%e(i) &
        + b%dt/4d0 * ( b%sOb(1)%e(i) + b%sOb(2)%e(i) )
      b%tOb(0)%tD(i) = b%dep(i) + b%tOb(0)%e(i)
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO SCHEDULE(dynamic,1000)
    do i = 1, b%npt
      j = b%npt + i          
      b%sOb(2)%p(i) = b%gXPQ(i)
      b%sOb(2)%q(i) = b%gXPQ(j)        
                
      b%tOb(0)%p(i) = b%tOb(1)%p(i) &
        + b%dt/4d0 * ( b%sOb(1)%p(i) + b%sOb(2)%p(i) )
      b%tOb(0)%q(i) = b%tOb(1)%q(i) &
        + b%dt/4d0 * ( b%sOb(1)%q(i) + b%sOb(2)%q(i) )
    enddo        
    !$OMP END DO NOWAIT
    call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD) 
    !$OMP END PARALLEL          
    !!---------End updateSoln S2---------!!


    !!-------------SSPRK3 S3-------------!!
    b%ur = b%tOb(0)%p / b%tOb(0)%tD
    b%vr = b%tOb(0)%q / b%tOb(0)%tD
    b%pbpr = b%tOb(0)%p / b%por
    b%qbpr = b%tOb(0)%q / b%por   
    rTime = b%tOb(1)%rtm + b%dt/2d0
    call b%dynaMatrices(rTime,b%tOb(0)%tD, b%ur, b%vr)
    call b%solveAll(rTime, b%tOb(0)%p, b%tOb(0)%q, &
      b%pbpr, b%qbpr, b%presr, b%tOb(0)%e, &
      b%gXW, b%gXE, b%gXPQ, b%gRE, b%gRPQ, b%sysC)    
    !!-----------End SSPRK3 S3-----------!!


    !!-----------updateSoln S3-----------!!
    b%sysT(1) = b%sysT(1) + 1d0*(b%sysC(8)-b%sysC(7))/b%sysRate
    b%sysT(2) = b%sysT(2) + 1d0*(b%sysC(6)-b%sysC(5))/b%sysRate

    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j)
    !$OMP DO SCHEDULE(dynamic,1000)
    do i = 1, b%npl          
      b%sOb(3)%e(i) = b%gXE(i)                    
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO SCHEDULE(dynamic,1000)
    do i = 1, b%npt
      j = b%npt + i          
      b%sOb(3)%p(i) = b%gXPQ(i)
      b%sOb(3)%q(i) = b%gXPQ(j)                              
    enddo        
    !$OMP END DO NOWAIT    
    !$OMP END PARALLEL          
    !!---------End updateSoln S3---------!!
    

  end subroutine timeStepSSPRK3
!!-----------------------End timeStepSSPRK3------------------------!!



!!----------------------------diriBCEta----------------------------!!
  subroutine diriBCEta(b,mat,rTt0)
  implicit none

    class(bsnqCase),intent(in)::b    
    real(kind=C_K2),intent(in)::rTt0
    real(kind=C_K2),intent(inout)::mat(b%npl)    
    integer(kind=C_K1)::i,i2,j2
    real(kind=C_K2)::leta, letadt

    !! Note : Consistent with SemiDirect only
    call b%wvF%getEta(rTt0, leta, letadt)    
    do i=1,b%nbndp
      i2=b%bndP(i)
      j2=b%bndPT(i)
      if(i2.gt.b%npl)cycle !Linear only
      if(j2.eq.11)then        
        mat(i2)=leta

      elseif(j2.eq.14)then
        mat(i2)=0d0

      endif
    enddo
  end subroutine diriBCEta
!!--------------------------End diriBCEta--------------------------!!



!!----------------------------diriBCPQ-----------------------------!!
  subroutine diriBCPQ(b,p,q,rTt0)
  implicit none

    class(bsnqCase),intent(in)::b    
    real(kind=C_K2),intent(in)::rTt0
    real(kind=C_K2),intent(inout)::p(b%npt),q(b%npt)
    integer(kind=C_K1)::i,i2,j2
    real(kind=C_K2)::lp,lq,lpdt,lqdt

    !! Note : Consistent with SemiDirect only
    call b%wvF%getPQ(rTt0, lp, lq, lpdt, lqdt)    
    do i=1,b%nbndp
      i2=b%bndP(i)
      j2=b%bndPT(i)
      if(j2.eq.11)then        
        p(i2)=lp
        q(i2)=lq

      elseif((j2.eq.12).or.(j2.eq.14))then
        p(i2)=0d0
        q(i2)=0d0

      elseif(j2.eq.13)then
        if( abs(b%bndPN(i2,1)) .gt. 0.86d0 )then !degLT30
          p(i2) = -q(i2)*b%bndPN(i2,2)/b%bndPN(i2,1)
          !write(*,*) '!!!!222',i2,abs(b%bndPN(i2,1))
        else
          q(i2) = -p(i2)*b%bndPN(i2,1)/b%bndPN(i2,2)
        endif        

      endif
    enddo
  end subroutine diriBCPQ
!!--------------------------End diriBCPQ---------------------------!!



!!---------------------------diriBCEtaDt---------------------------!!
  subroutine diriBCEtaDt(b,mat,rTt0)
  implicit none

    class(bsnqCase),intent(in)::b
    real(kind=C_K2),intent(in)::rTt0
    real(kind=C_K2),intent(inout)::mat(b%npl)    
    integer(kind=C_K1)::i,i2,j2
    real(kind=C_K2)::leta, letadt

    !! Note : Consistent with SemiDirect only
    call b%wvF%getEta(rTt0, leta, letadt)    
    !write(201,'(3F15.6)',advance='no')rTt0,tmpr1,tmpr2
    do i=1,b%nbndp
      i2=b%bndP(i)
      j2=b%bndPT(i)
      if(i2.gt.b%npl)cycle !Linear only
      if(j2.eq.11)then        
        mat(i2)=letadt

      elseif(j2.eq.14)then
        mat(i2)=0d0

      endif
    enddo
  end subroutine diriBCEtaDt
!!-------------------------End diriBCEtaDt-------------------------!!



!!---------------------------diriBCPQDt----------------------------!!
  subroutine diriBCPQDt(b,mat,rTt0)
  implicit none

    class(bsnqCase),intent(in)::b
    real(kind=C_K2),intent(in)::rTt0
    real(kind=C_K2),intent(inout)::mat(2*b%npt)
    integer(kind=C_K1)::i,i2,j2
    real(kind=C_K2)::lp, lq, lpdt, lqdt

    !! Note : Consistent with SemiDirect only
    call b%wvF%getPQ(rTt0, lp, lq, lpdt, lqdt)             
    !write(201,'(2F15.6)',advance='no')tmpr1,tmpr3
    do i=1,b%nbndp
      i2=b%bndP(i)
      j2=b%bndPT(i)
      if(j2.eq.11)then        
        mat(i2)=lpdt
        mat(b%npt+i2)=lqdt

      elseif((j2.eq.12).or.(j2.eq.14))then
        mat(i2)=0d0
        mat(b%npt+i2)=0d0

      elseif(j2.eq.13)then
        if( abs(b%bndPN(i2,1)) .gt. 0.86d0 )then !degLT30
          mat(i2) = 0d0
          !write(*,*) '!!!!111',i2,abs(b%bndPN(i2,1))
        else
          mat(b%npt+i2)=0d0
        endif        

      endif
    enddo
  end subroutine diriBCPQDt
!!-------------------------End diriBCPQDt--------------------------!!



!!-----------------------------initMat-----------------------------!!
  subroutine initMat(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::i,i1,j,j1    

    i=b%npl
    j=b%npt
    i1=b%ivl(0)
    j1=b%ivq(0)

    !Allocations
    allocate(b%gXE(i),b%gXPQ(2*j),b%gXW(i))
    allocate(b%por(j), b%vec6Tmp(j), b%vec6Tmp2(j))
    allocate(b%ur(j),b%vr(j),b%pbpr(j),b%qbpr(j))
    allocate(b%uhr(j),b%vhr(j))

    allocate(b%massW(i1*i),b%massE(i1*i))
    allocate(b%mass1(j1*j),b%mass2(i1*i))    
    allocate(b%gBs1(j1*j),b%gBs2(j1*j))
    allocate(b%gBs3(j1*j),b%gBs4(j1*j))
    allocate(b%gCxF(j1*i),b%gCyF(j1*i),b%gDMat(i1*i))
    allocate(b%gBs5(i1*j),b%gBs6(i1*j),b%absC(j))
    allocate(b%gGx(i1*j),b%gGy(i1*j),b%gNAdv(j1*j))
    allocate(b%gPGx(j1*j),b%gPGy(j1*j))
    allocate(b%gFBs1(j1*j),b%gFBs2(j1*j))
    allocate(b%gFBs3(j1*j),b%gFBs4(j1*j),b%gFW(i1*i))
    allocate(b%aFull(b%ivf(0)*2*j))
    allocate(b%rowMaxW(i),b%rowMaxE(i),b%rowMaxPQ(2*j))
    allocate(b%gRE(i),b%gRPQ(2*j))
    allocate(b%etaMax(i),b%etaMin(i),b%presr(j))
    allocate(b%ele6x6(b%nele,36),b%ele6x3(b%nele,18))
    allocate( b%wvAng(j), b%botFricN6(j) )

    b%Sz(1)=i1*i ![3x3] ![ivl(0) * npl]
    b%Sz(2)=j1*i ![3x6] ![ivq(0) * npl]
    b%Sz(3)=i1*j ![6x3] ![ivl(0) * npt]
    b%Sz(4)=j1*j ![6x6] ![ivq(0) * npt]
    b%por=1d0
    
    b%absC=0d0
    if(b%absOn)then
      do i=1,b%absOb(1)%N
        call b%absOb(i)%calcAbsC(b%npt,b%cor,b%absC)
      enddo  
    endif

    b%nTOb=2
    allocate(b%tOb(0:b%nTOb-1))
    do i=0,b%nTOb-1
      call b%tOb(i)%initBsnqVars(b%npl,b%npt)
    enddo

    b%nSOb=4
    allocate( b%sOb(b%nSOb))
    do i=1,b%nSOb
      call b%sOb(i)%initBsnqVars(b%npl,b%npt)
    enddo    

    b%tStep=0

    if(b%resume)then
      ! [Note]:
      ! Check bsnqM/Dev - Logs/log_bsnqM_v0003.md
      ! to understand why the variables 
      ! gXW, gXE, gXPQ are being written in 
      ! resume file.
      call b%readResume
    
    else
      b%tOb(0)%rtm=0d0
      b%tOb(0)%e=0d0
      b%tOb(0)%p=0d0
      b%tOb(0)%q=0d0      

      b%gXW = 0d0
      b%gXE = 0d0
      b%gXPQ = 0d0

      ! Initialis eta to ship shape if instructed by user
      if(b%presOn)then
        b%presr=0d0
        do i=1,b%sh(1)%totNShip
          if(b%sh(i)%initEtaFlag)then
            call b%sh(i)%getPress(b%tOb(0)%rtm,b%npt,b%cor,b%vec6Tmp)
            b%presr=b%presr+b%vec6Tmp
          endif
        enddo
        b%tOb(0)%e = -b%presr
      endif
    endif

    ! call solitICFromFile(b%npl, b%npt, b%cor, &
    !   b%tOb(0)%e, b%tOb(0)%p, b%tOb(0)%q)
    
    b%tOb(0)%tD(1:b%npl)=b%dep(1:b%npl)+b%tOb(0)%e
    call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD)        
    
    b%etaMin=b%tOb(0)%e
    b%etaMax=b%tOb(0)%e    
    b%wvAng=0d0
    b%botFricN6 = 0d0

    b%presr=0d0    

    ! call testMls2DDx
    ! stop    

    call paralution_init(b%nthrd)    

    !call b%outputXML

    write(9,*)"[MSG] Done initMat"
    write(9,*)

  end subroutine initMat
!!---------------------------End initMat---------------------------!!



!!-------------------------writeWaveProbe--------------------------!!
  subroutine writeWaveProbe(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::nq(6),i,i2,k
    real(kind=C_K2)::N3(3),N6(6),tmpr1,tmpr2,tmpr3

    !Writing wave probes value
    k=b%wpEle(-1)
    write(k,'(F15.6)',advance='no')b%tOb(0)%rtm
    do i=1,b%wpEle(0)

      tmpr1=0d0
      tmpr2=0d0
      tmpr3=0d0

      if(b%wpEle(i).ne.-1)then

        nq=b%conn(b%wpEle(i),:)
        call fem_N6i(b%wpLoc(i,3), b%wpLoc(i,4), N6)
        call fem_N3i(b%wpLoc(i,3), b%wpLoc(i,4), N3)

        do i2=1,3
          tmpr1 = tmpr1 + b%tOb(0)%e(nq(i2)) * N3(i2) 
        enddo
        do i2=1,6
          tmpr2 = tmpr2 + b%tOb(0)%p(nq(i2)) * N6(i2) 
          tmpr3 = tmpr3 + b%tOb(0)%q(nq(i2)) * N6(i2)                     
        enddo

      endif

      write(k,'(5F15.6)',advance='no')b%wpLoc(i,1), b%wpLoc(i,2),&
        tmpr1, tmpr2, tmpr3
    enddo
    write(k,*)



    !Write vertVel probes value
    if(b%vvPrb%np .gt. 0)then
      k = b%vvPrb%fileid
      write(k,'(F15.6)',advance='no')b%tOb(0)%rtm
      
      call b%getVertVel( b%vvPrb%np, b%vvPrb%x, b%vvPrb%y, &
        b%vvPrb%z, b%vvPrb%u, b%vvPrb%v, b%vvPrb%w, &
        b%vvPrb%p, b%vvPrb%eta, &
        b%vvPrb%wrki, b%vvPrb%wrkr, b%vvPrb%err, tmpr1 )

      do i = 1, b%vvPrb%np
        if(b%vvPrb%err(i).eq.0)then   !noError
          write(k,'(8F15.6)',advance='no')b%vvPrb%x(i), &
            b%vvPrb%y(i), b%vvPrb%z(i), b%vvPrb%p(i), &
            b%vvPrb%u(i), b%vvPrb%v(i), b%vvPrb%w(i), &
            b%vvPrb%eta(i)
        else
          write(k,'(3F15.6, 4A15, F15.6)',advance='no')&
            b%vvPrb%x(i), b%vvPrb%y(i), b%vvPrb%z(i), &
            '---', '---', '---', '---', b%vvPrb%eta(i)
        endif 
      enddo
      write(k,*)
    endif


  end subroutine writeWaveProbe
!!-----------------------End writeWaveProbe------------------------!!



!!--------------------------getEtaPQForXY--------------------------!!
  subroutine getEtaPQForXY(b,np,xin,yin,eta,p,q,wrki,wrkr,err,rTime)
  implicit none

    class(bsnqCase),intent(in)::b    
    integer(kind=C_K1),intent(in)::np
    real(kind=C_K2),intent(in)::xin(np),yin(np)
    integer(kind=C_K1),intent(out)::wrki(np),err(np)
    real(kind=C_K2),intent(out)::eta(np),p(np),q(np),wrkr(np,2)
    real(kind=C_K2),intent(out)::rTime

    integer(kind=C_KCLK)::lsysC(2)
    integer(kind=C_K1)::nq(6),i,k
    real(kind=C_K2)::wei(6),etaLoc,pLoc,qLoc,hLoc

    call system_clock(lsysC(1)) 

    call b%findEleForLocXY3(np,xin,yin,wrki,wrkr(:,1),wrkr(:,2))

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i, nq, wei, hLoc, etaLoc, pLoc, qLoc, k)
    !$OMP DO SCHEDULE(dynamic,10)
    do i=1,np

      if(wrki(i).eq.-1)then
        err(i)=1
        eta(i)=0d0
        p(i)=0d0
        q(i)=0d0   
        cycle     
      endif

      nq = b%conn(wrki(i),:)
      call fem_N6i(wrkr(i,1),wrkr(i,2),wei)

      hLoc=0d0
      etaLoc=0d0
      pLoc=0d0
      qLoc=0d0

      do k=1,6
        hLoc = hLoc + wei(k) * b%dep(nq(k))
        etaLoc = etaLoc + wei(k) * b%tOb(0)%tD(nq(k)) !Using totalDep 
        pLoc = pLoc + wei(k) * b%tOb(0)%p(nq(k))
        qLoc = qLoc + wei(k) * b%tOb(0)%q(nq(k))
      enddo
      etaLoc = etaLoc - hLoc

      err(i)=0
      eta(i)=etaLoc
      p(i)=pLoc
      q(i)=qLoc

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

    call system_clock(lsysC(2)) 
    rTime = 1d0 * (lsysC(2) - lsysC(1)) / b%sysRate   

  end subroutine getEtaPQForXY

!!------------------------End getEtaPQForXY------------------------!!



!!-------------------------setParalutionLS-------------------------!!
subroutine setParalutionLS(b)
  implicit none

    class(bsnqCase),intent(inout)::b

    b%paralsW  = createLSObj();
    b%paralsE  = createLSObj();
    b%paralsPQ = createLSObj();

    ! errLim is setting the absolute tolerence
    call initLSSys(b%paralsW, b%npl, b%nnzl, b%ivsl, b%jvsl, &
      b%gMW, b%errLim, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, b%maxIter)

    call initLSSys(b%paralsE, b%npl, b%nnzl, b%ivsl, b%jvsl, &
      b%gME, b%errLim, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, b%maxIter)

    call initLSSys(b%paralsPQ, 2*b%npt, b%nnzf, b%ivsf, b%jvsf, &
      b%gMPQ, b%errLim, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, b%maxIter)

    ! call checkLSObj(b%paralsW)
    ! call checkLSObj(b%paralsE)
    ! call checkLSObj(b%paralsPQ)            

    write(9,*)
    write(9,'(" [MSG] Created Linear Solver objects in C++")')
    write(9,*)
    !stop
end subroutine setParalutionLS
!!-----------------------End setParalutionLS-----------------------!!



!!--------------------------initBsnqVars---------------------------!!
  subroutine initBsnqVars(b,npl,npt)
  use bsnqGlobVars
  implicit none

    class(bsnqVars),intent(inout)::b
    integer(kind=C_K1),intent(in)::npl,npt
    logical(kind=C_LG)::tmp

    tmp=allocated(b%e).or.allocated(b%p).or.&
      allocated(b%q).or.allocated(b%tD)    
    if(tmp)then
      write(9,'(" [ERR] Already allocated bsnqVars")')
      stop
    endif
    allocate(b%e(npl),b%p(npt),b%q(npt),b%tD(npt))
    b%npl=npl
    b%npt=npt
    b%e=0d0
    b%p=0d0
    b%q=0d0
    b%tD=0d0

  end subroutine initBsnqVars
!!------------------------End initBsnqVars-------------------------!!


end module bsnqModule