!!-----------------------------Custom------------------------------!!
subroutine solveSys2(obj,np,nnz,gB,gX,&
  iter,resnorm,ier)
use bsnqGlobVars
implicit none

  interface    
    subroutine solveLSObj( obj, rhs, x, iter, resnorm, &
      ierr ) BIND(C, name="solveLSObj")

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE

      type(C_PTR), value, intent(in) :: obj      
      integer(kind=C_INT),        intent(out) :: iter, ierr
      real(kind=C_DOUBLE),        intent(out) :: resnorm
      type(C_PTR),         value, intent(in)  :: rhs
      type(C_PTR),         value              :: x

    end subroutine solveLSObj
  end interface
  
  type(C_PTR),intent(in)::obj
  integer(kind=C_K1),intent(in)::np,nnz  
  integer(kind=C_K1),intent(out)::iter,ier  
  real(kind=C_DOUBLE),intent(in),target::gB(np)  
  real(kind=C_DOUBLE),intent(out),target::gX(np)
  real(kind=C_K2),intent(out)::resnorm  
  
  
  call solveLSObj(obj,C_LOC(gB), &    
    C_LOC(gX), iter, resnorm, ier)  

end subroutine solveSys2



subroutine initLSSys(lsObj, np,nnz,iv,jv,gA, &
  atol, rtol, div, maxIter)
use bsnqGlobVars
implicit none

  interface
    subroutine initLSObj(obj, np, nnz, ivCSR, jvCSR, matCSR, & 
      atol, rtol, div, maxiter) BIND(C, name="initLSObj")
      
      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE
      
      type(C_PTR), value, intent(in) :: obj
      integer(C_INT), value, intent(in)  :: np, nnz
      type(C_PTR), value, intent(in)  :: ivCSR, jvCSR, matCSR
      integer(kind=C_INT), value, intent(in)  :: maxiter
      real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
    end subroutine initLSObj
  end interface

  type(C_PTR)::lsObj
  integer(kind=C_K1),intent(in)::np,nnz, maxiter
  integer(kind=C_INT),intent(in),target::iv(np+1),jv(nnz)
  real(kind=C_K2),intent(in)::atol, rtol, div
  real(kind=C_DOUBLE),intent(in),target::gA(nnz)
  
  call initLSObj(lsObj, np, nnz, C_LOC(iv), C_LOC(jv),&
    C_LOC(gA), atol, rtol, div, maxiter)
end subroutine initLSSys
!!---------------------------End Custom----------------------------!!



!!-----------------------------Default-----------------------------!!
subroutine solveSys(np,nnz,iv,jv,gA,gB,gX,errLim,maxiter,&
  iter,resnorm,ier)
use bsnqGlobVars
implicit none

  interface    
    subroutine paralution_fortran_solve_csr( n, m, nnz, solver, &
      mformat, preconditioner, pformat, &
      ivCSR, jvCSR, rval, rhs, atol, rtol, div, maxiter, basis, &
      p, q, x, iter, resnorm, ierr ) BIND(C, name ="paralution_fortran_solve_csr")

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

      integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
      real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
      integer(kind=C_INT),        intent(out) :: iter, ierr
      real(kind=C_DOUBLE),        intent(out) :: resnorm
      type(C_PTR),         value, intent(in)  :: ivCSR, jvCSR, rval, rhs
      type(C_PTR),         value              :: x
      character(kind=C_CHAR)                  :: solver, mformat, preconditioner, pformat

    end subroutine paralution_fortran_solve_csr
  end interface
  
  integer(kind=C_K1),intent(in)::np,nnz,maxIter
  integer(kind=C_INT),intent(in),target::iv(np+1),jv(nnz)
  integer(kind=C_K1),intent(out)::iter,ier
  real(kind=C_K2),intent(in)::errLim  
  real(kind=C_DOUBLE),intent(in),target::gA(nnz),gB(np)  
  real(kind=C_DOUBLE),intent(out),target::gX(np)
  real(kind=C_K2),intent(out)::resnorm  

  call paralution_fortran_solve_csr(np,np,nnz,&
    'BiCGStab' // C_NULL_CHAR, &
    'CSR' // C_NULL_CHAR, &
    'None' // C_NULL_CHAR, &
    'CSR' // C_NULL_CHAR, &
    C_LOC(iv), C_LOC(jv), &
    C_LOC(gA), C_LOC(gB), &
    errLim, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, maxIter, &
    30, 0, 1, C_LOC(gX), iter, resnorm, ier)

end subroutine solveSys
!!---------------------------End Default---------------------------!!