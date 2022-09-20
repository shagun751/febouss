!!------------------------rk4Interpolation-------------------------!!
module rk4InterpMod
use bsnqGlobVars
implicit none

  type, public :: rk4InterpTyp        
    real(kind=C_K2)::M(4,4)
  
  contains
    procedure ::  setrk4InterpMatrix
    procedure ::  get_ktilde
  end type rk4InterpTyp
  
contains

  subroutine setrk4InterpMatrix(f, dtC, dtF)
  implicit none

    class(rk4InterpTyp),intent(inout)::f
    real(kind=C_K2),intent(in)::dtC, dtF
    real(kind=C_K2)::tr

    tr = dtF/dtC
    
    f%M(1,:) = (/ 1d0, 0d0, 0d0, 0d0 /)

    f%M(2,:) = (/ 1d0 - tr, tr, 0d0, 0d0 /) 

    f%M(3,:) = (/ 1d0 - tr, tr - tr**2, tr**2, 0d0 /)
   
    f%M(4,:) = (/ 1d0 - 2d0*tr + tr**3, 2d0*tr*(1d0 - tr),&
      2d0*tr*tr*(1d0 - tr), tr**3 /)

  end subroutine setrk4InterpMatrix



  subroutine get_ktilde(f, k, ktil)
  implicit none

    class(rk4InterpTyp),intent(in)::f
    real(kind=C_K2),intent(in)::k(4)
    real(kind=C_K2),intent(out)::ktil(4)

    ktil(1) = k(1)
    ktil(2) = k(1)*f%M(2,1) + k(2)*f%M(2,2)
    ktil(3) = k(1)*f%M(3,1) + k(2)*f%M(3,2) + k(3)*f%M(3,3)
    ktil(4) = k(1)*f%M(4,1) + k(2)*f%M(4,2) + k(3)*f%M(4,3) + k(4)*f%M(4,4)
    
  end subroutine get_ktilde

end module rk4InterpMod
!!----------------------End rk4Interpolation-----------------------!!