! **************************************************************************
!
!    PARALUTION   www.paralution.com
!
!    Copyright (C) 2015  PARALUTION Labs UG (haftungsbeschränkt) & Co. KG
!                        Am Hasensprung 6, 76571 Gaggenau
!                        Handelsregister: Amtsgericht Mannheim, HRA 706051
!                        Vertreten durch:
!                        PARALUTION Labs Verwaltungs UG (haftungsbeschränkt)
!                        Am Hasensprung 6, 76571 Gaggenau
!                        Handelsregister: Amtsgericht Mannheim, HRB 721277
!                        Geschäftsführer: Dimitar Lukarski, Nico Trost
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! **************************************************************************

!
!
! PARALUTION version 1.1.0 
!
!
program paralution_solver

  use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

  implicit none

  interface
    ! Shagun modify 2017_08_14
    subroutine paralution_init(nthreads) BIND(C)
      use, intrinsic :: ISO_C_BINDING, only : C_INT
      integer(kind=C_INT), value, intent(in)  :: nthreads
    end subroutine paralution_init

    subroutine paralution_stop() BIND(C)
    end subroutine paralution_stop    

    subroutine paralution_fortran_solve_coo( n, m, nnz, solver, mformat, preconditioner, pformat,    &
    &                                        rows, cols, rval, rhs, atol, rtol, div, maxiter, basis, &
    &                                        p, q, x, iter, resnorm, ierr ) BIND(C)

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

      integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
      real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
      integer(kind=C_INT),        intent(out) :: iter, ierr
      real(kind=C_DOUBLE),        intent(out) :: resnorm
      type(C_PTR),         value, intent(in)  :: rows, cols, rval, rhs
      type(C_PTR),         value              :: x
      character(kind=C_CHAR)                  :: solver, mformat, preconditioner, pformat

    end subroutine paralution_fortran_solve_coo
  end interface

  integer, parameter    :: infile = 10
  integer(kind=C_INT)   :: n, m, nnz, fnz, i, j, iter, ierr
  real(kind=C_DOUBLE)   :: resnorm
  integer, dimension(8) :: tbegin, tend
  real(kind=8)          :: tsolver

  logical               :: sym = .false.

  integer(kind=C_INT), allocatable, target :: rows(:), cols(:)
  real(kind=C_DOUBLE), allocatable, target :: ival(:), rval(:), cval(:)
  real(kind=C_DOUBLE), allocatable, target :: rhs(:), x(:)

  character(len=10)  :: rep
  character(len=7)   :: field
  character(len=19)  :: symm
  character(len=128) :: arg,text


  ! Get command line option to read file
  do i = 1, iargc()
    call getarg(i, arg)
  end do
  open( unit = infile, file = arg, action = 'read', status = 'old' )

  call mminfo( infile, rep, field, symm, n, m, fnz )

  nnz = fnz
  if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
    nnz = 2 * ( fnz - n ) + n
    sym = .true.
  end if

  ! Allocate memory for COO format specific arrays
  allocate( rows(nnz), cols(nnz) )
  allocate( ival(nnz), rval(nnz), cval(nnz) )

  ! Read COO matrix in matrix market format
  call mmread( infile, rep, field, symm, n, m, fnz, nnz, &
  &            rows, cols, ival, rval, cval )

  ! Close file
  close( unit = infile )

  ! Deallocate arrays we do not use
  deallocate( ival, cval )

  ! Fill 2nd half of matrix if symmetric
  if ( sym ) then
    j = fnz + 1
    do i = 1, fnz
      if ( rows(i) .ne. cols(i) ) then
        rows(j) = cols(i)
        cols(j) = rows(i)
        rval(j) = rval(i)
        j = j + 1
      end if
    end do
  end if

  ! Allocate and initialize rhs and solution vector
  allocate( rhs(n), x(n) )
  do i = 1, n
    rhs(i) = 1._C_DOUBLE
    x(i)   = 0._C_DOUBLE
  end do

  text=arg(1:len_trim(arg))//'2'
  open(111,file=text(1:len_trim(text)))

  do i=1,n
    read(111,*)rhs(i)
  enddo
  close(111)
  write(*,*)'[Msg] RHS imported!',sum(rhs)

  ! Print L2 norm of solution vector
  write(*,fmt='(A,F0.2)') '(Fortran) Initial L2 Norm(x) = ', sqrt( sum( x**2 ) )

  ! Shagun modify 2017_08_14
  ! Initialize PARALUTION backend
  call paralution_init(2)

  call date_and_time(values = tbegin)

  ! Run paralution C function for COO matrices
  ! Doing a GMRES with MultiColored ILU(1,2) preconditioner
  ! Check paralution documentation for a detailed argument explanation
  call paralution_fortran_solve_coo( n, m, nnz,                                          &
  &                                  'GMRES' // C_NULL_CHAR,                                &
  &                                  'CSR' // C_NULL_CHAR,                               &
  &                                  'None' // C_NULL_CHAR,                   &
  &                                  'CSR' // C_NULL_CHAR,                               &
  &                                  C_LOC(rows), C_LOC(cols), C_LOC(rval), C_LOC(rhs),  &
  &                                  1e-15_C_DOUBLE, 1e-8_C_DOUBLE, 1e+8_C_DOUBLE, 5000, &
  &                                  30, 0, 1, C_LOC(x), iter, resnorm, ierr )

  call date_and_time(values = tend)

  tbegin = tend - tbegin
  tsolver = 0.001 * tbegin(8) + tbegin(7) + 60 * tbegin(6) + 3600 * tbegin(5)
  write(*,fmt='(A,F0.2,A)') '(Fortran) Solver ended after ', tsolver,'sec.'

  ! Print solver details
  if ( ierr .eq. 0 ) then
    write(*,fmt='(A,I0,A,E11.5,A)') '(Fortran) Solver took ', iter, ' iterations with residual norm ', resnorm, '.'
    write(*,fmt='(A,F0.2)') '(Fortran) Final L2 Norm(x)   = ', sqrt( sum( x**2 ) )
  else
    write(*,fmt='(A,I0)') '(Fortran) Solver returned status code ', ierr
  end if

  do i=1,n
    write(10,*) x(i)
  end do

  deallocate( rows, cols, rval, rhs, x )

  ! Stop PARALUTION backend
  call paralution_stop

end program paralution_solver

