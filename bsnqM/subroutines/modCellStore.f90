!!--------------------------cellStoreMod---------------------------!!
module cellStoreMod
use bsnqGlobVars
implicit none
  
  type, private ::cellTyp
    integer(kind=C_K1)::nm=0, nmax=0
    integer(kind=C_K1),allocatable::m(:)
  end type cellTyp



  type, public :: cellofEleTyp

    integer(kind=C_K1)::ntot, nx, ny
    real(kind=C_K2)::xmin, xmax, ymin, ymax, cellR
    type(cellTyp),allocatable::cell(:)
  contains
    procedure :: setEleCells      
    procedure :: getCellNo
  end type cellofEleTyp

contains



!!-------------------------setEleCells-------------------------!!
  subroutine setEleCells(f, xmin, ymin, xmax, ymax, cellR, &
    np, nele, corx, cory, conn)
  implicit none

    class(cellofEleTyp),intent(inout)::f    
    integer(kind=C_K1),intent(in)::np, nele, conn(nele,6)
    real(kind=C_K2),intent(in)::xmin, xmax, ymin, ymax, cellR    
    real(kind=C_K2),intent(in)::corx(np), cory(np)

    integer(kind=C_KCLK)::lsysC(3)
    integer(kind=C_K1)::iele, nl(3), ip, i, i11, i12, nm
    integer(kind=C_K1)::icell
    integer(kind=C_K1),allocatable::tmpi(:)
    real(kind=C_K2)::exy(3,2)

    call system_clock(lsysC(1), lsysC(3))

    f%xmin = xmin
    f%ymin = ymin
    f%xmax = xmax
    f%ymax = ymax
    f%cellR = cellR

    f%nx = floor((f%xmax - f%xmin)/f%cellR) + 11
    f%ny = floor((f%ymax - f%ymin)/f%cellR) + 11
    f%ntot = f%nx * f%ny 

    allocate(f%cell(f%ntot))
    allocate(tmpi(np))

    do i = 1, f%ntot
      f%cell(i)%nm = 0
      f%cell(i)%nmax = nele
      allocate( f%cell(i)%m(nele) )
    enddo


    do iele = 1, nele
      nl = conn(iele, 1:3)
      do ip = 1, 3
        call f%getCellNo( corx(nl(ip)), cory(nl(ip)), &
          i11, i12, icell )
        nm = f%cell(icell)%nm
        if( count(f%cell(icell)%m(1:nm) .eq. iele) .eq. 0 )then
          f%cell(icell)%nm = f%cell(icell)%nm + 1
          f%cell(icell)%m(nm+1) = iele
        endif
      enddo
    enddo


    do i = 1, f%ntot
      nm = f%cell(i)%nm
      
      if(nm.le.0)then
        deallocate(f%cell(i)%m)
        f%cell(i)%nm = 0
        f%cell(i)%nmax = 0
        cycle
      endif

      tmpi( 1:nm ) = f%cell(i)%m( 1:nm )
      deallocate(f%cell(i)%m)
      allocate( f%cell(i)%m(1:nm) )
      f%cell(i)%m( 1:nm ) = tmpi( 1:nm )
    enddo

    deallocate(tmpi)

    call system_clock(lsysC(2))
    write(9,*)
    write(9,'(" [INF] setEleCells")')
    write(9,'(" [---] xmin, xmax : ",2F15.6)') f%xmin, f%xmax
    write(9,'(" [---] ymin, ymax : ",2F15.6)') f%ymin, f%ymax
    write(9,'(" [---] cellR, NX, NY, NTOT : ",F15.6, 3I10)') &
      f%cellR, f%nx-10, f%ny-10, f%ntot
    write(9,'(" [TIM] setEleCells Time : ", F15.6)') &
      1d0*(lsysC(2)-lsysC(1))/lsysC(3)
    write(9,*)

  end subroutine setEleCells
!!-----------------------End setEleCells-----------------------!!



!!--------------------------getCellNo--------------------------!!
  subroutine getCellNo(f, xin, yin, cellx, celly, celli)
  implicit none

    class(cellofEleTyp),intent(in)::f    
    real(kind=C_K2),intent(in)::xin, yin
    integer(kind=C_K1),intent(out)::cellx, celly, celli

    cellx = floor( (xin - f%xmin) / f%cellR )
    celly = floor( (yin - f%ymin) / f%cellR )

    if(cellx.lt.0) cellx=0
    if(celly.lt.0) celly=0
    if(cellx.gt.f%nx-1) cellx=f%nx-1
    if(celly.gt.f%ny-1) celly=f%ny-1

    celli = cellx + (f%nx * celly) + 1

  end subroutine getCellNo

!!------------------------End getCellNo------------------------!!



end module cellStoreMod
!!------------------------End cellStoreMod-------------------------!!