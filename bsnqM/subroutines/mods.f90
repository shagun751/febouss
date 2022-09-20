!!----------------------------basicVars----------------------------!!
module bsnqGlobVars 
use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE
use, intrinsic :: ISO_C_BINDING, only : C_CHAR, C_NULL_CHAR, C_LOC 
implicit none
  
  integer,parameter::C_K1=C_INT,C_K2=C_DOUBLE,C_KSTR=256,C_LG=1  
  integer,parameter::C_KCLK=8

  !!-----------------------Constants-----------------------!!
  !! Gravity
  real(kind=C_K2),parameter::grav=9.81d0          
  !! PI
  real(kind=C_K2),parameter::pi=atan(1d0)*4d0       
  real(kind=C_K2),parameter::deg2rad=pi/180d0
  real(kind=C_K2),parameter::rad2deg=180d0/pi
  !! Water density
  real(kind=C_K2),parameter::rhoW=1000d0  
  !! Bsnq constant
  real(kind=C_K2),parameter::BsqC=1d0/15d0
  !!---------------------End Constants---------------------!!


  !!  mafi defintions
  !!  mafi(1)     Mesh File
  !!  mafi(2)     Paraview output
  !!  mafi(3)     Volume output
  !!  mafi(4)     <Unknown>
  !!  mafi(5)     Input file
  !!  mafi(6)     Porosity file
  !!  mafi(7)     Wave probes files

contains
  
  subroutine vecCross2D(v1x,v1y,v2x,v2y,res)
  implicit none

    real(kind=C_K2),intent(in)::v1x,v1y,v2x,v2y
    real(kind=C_K2),intent(out)::res

    res = v1x*v2y - v1y*v2x

  end subroutine vecCross2D


end module bsnqGlobVars
!!--------------------------End basicVars--------------------------!!






