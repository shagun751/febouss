!!------------------------ Version M 1.x.x ------------------------!!
!!  -> Quadratic + Linear
!!  -> Boundary - SemiDirect + Penalty + Gauss Seidel
!!    -> 11 - Inlet - No Absorbing
!!    -> 12 - NoSlip Wall 
!!    -> 13 - Slip Wall - Rectangular wall only
!!    -> 14 - Sponge - BOUSS2D approach generalised input
!!    -> 15 - Outlet  - Not coded
!!  -> Porosity - Emergent structure only (not done yet)
!!    -> Generalised input
!!  -> Solver - Normalised
!!  -> Generalised input - v3
!!  -> GMRES
!!  -> Dirichlet BndCond for eta, p, q
!!  -> Paralution CSR
!!	-> XML output
!!-----------------------------------------------------------------!!
!! Time-Stepping : RK4 

program boussinesqQuad
use bsnqGlobVars
use bsnqModule
implicit none

!!--------------------------Declarations---------------------------!!
  type(bsnqCase)::bq  
  character(len=C_KSTR)::bqtxt
!!------------------------End Declarations-------------------------!!
  
  call getarg(1,bq%probname)  
  do while(len_trim(bq%probname).lt.1)
    write(*,'(A)',advance='no')"Enter Problem Name:"
    read(*,*)bq%probname
  enddo
  write(*,*)"Problem Name: "//trim(bq%probname)

  bqtxt=trim(bq%probname)//'.rout'
  open(9,file=trim(bqtxt))  

  call system_clock(bq%sysC(1))
  call bq%meshRead
  call bq%femInit
  call bq%setRun 
  call bq%setMFree 
  call bq%initMat  
  call bq%statMatrices    
  call bq%setParalutionLS  

  call bq%caseOutputs

  select case (bq%tStepMethod)
  case (0) !RK4
  
    do while(abs(bq%tOb(0)%rtm-bq%endTime).gt.bq%dt/2d0)
    
      call bq%preInstructsRK4

      call bq%timeStepRK4
    
      call bq%postInstructs(0) !RK4

      call bq%caseOutputs

    enddo


  case (1) !AdBa3

    do while(bq%tStep.le.2) 
      call bq%timeStepEuEx1
      call bq%postInstructs(11) !EuEx1
      call bq%caseOutputs      
    enddo

    do while(abs(bq%tOb(0)%rtm-bq%endTime).gt.bq%dt/2d0) 
      call bq%timeStepEuEx1
      call bq%postInstructs(1) !AdBa3
      call bq%caseOutputs
    enddo


  case (2) !SSPRK3
  
    do while(abs(bq%tOb(0)%rtm-bq%endTime).gt.bq%dt/2d0)        

      call bq%timeStepSSPRK3
    
      call bq%postInstructs(2) !SSPRK3

      call bq%caseOutputs

    enddo

  case default
    write(9,'(" [ERR] Wrong time-step method!")')
    stop
  end select

  call system_clock(bq%sysC(2))
  close(bq%wpEle(-1))
  close(bq%vvPrb%fileid)
  write(9,*)"[MSG] boussinesqQuad End"
  write(9,'(" [TIM] ",F15.4)')1d0*(bq%sysC(2)-bq%sysC(1))/bq%sysRate
  close(9)  
end program boussinesqQuad