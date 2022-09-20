subroutine fem_N6i(ep,et,shF)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(in)::ep,et
  real(kind=C_K2),intent(out)::shF(6)

  shF(1) = 1d0 - 3d0*ep - 3d0*et + 2d0*ep*ep &
    + 2d0*et*et + 4d0*ep*et
  shF(2) = 2d0*ep*ep - ep
  shF(3) = 2d0*et*et - et
  shF(4) = 4d0*ep - 4d0*ep*ep - 4d0*ep*et
  shF(5) = 4d0*ep*et
  shF(6) = 4d0*et - 4d0*et*et - 4d0*ep*et

end subroutine fem_N6i

subroutine fem_N3i(ep,et,shF)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(in)::ep,et
  real(kind=C_K2),intent(out)::shF(3)

  shF(1) = 1d0 - ep - et
  shF(2) = ep
  shF(3) = et  

end subroutine fem_N3i


subroutine fem_N6i_Sc3_dN3jdx(mat,h1,h2,h3,b11,b12)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::mat(6,3)
  real(kind=C_K2),intent(in)::h1,h2,h3,b11,b12

  mat(1,1)=(-(1d0/120d0))*(b11 + b12)*(2d0*h1 - h2 - h3)
  mat(1,2)=(1d0/120d0)*b11*(2d0*h1 - h2 - h3)
  mat(1,3)=(1d0/120d0)*b12*(2d0*h1 - h2 - h3)

  mat(2,1)=(1d0/120d0)*(b11 + b12)*(h1 - 2d0*h2 + h3)
  mat(2,2)=(-(1d0/120d0))*b11*(h1 - 2d0*h2 + h3)
  mat(2,3)=(-(1d0/120d0))*b12*(h1 - 2d0*h2 + h3)

  mat(3,1)=(1d0/120d0)*(b11 + b12)*(h1 + h2 - 2d0*h3)
  mat(3,2)=(-(1d0/120d0))*b11*(h1 + h2 - 2d0*h3)
  mat(3,3)=(-(1d0/120d0))*b12*(h1 + h2 - 2d0*h3)

  mat(4,1)=(-(1d0/30d0))*(b11 + b12)*(2d0*h1 + 2d0*h2 + h3)
  mat(4,2)=(1d0/30d0)*b11*(2d0*h1 + 2d0*h2 + h3)
  mat(4,3)=(1d0/30d0)*b12*(2d0*h1 + 2d0*h2 + h3)

  mat(5,1)=(-(1d0/30d0))*(b11 + b12)*(h1 + 2d0*(h2 + h3))
  mat(5,2)=(1d0/30d0)*b11*(h1 + 2d0*(h2 + h3))
  mat(5,3)=(1d0/30d0)*b12*(h1 + 2d0*(h2 + h3))

  mat(6,1)=(-(1d0/30d0))*(b11 + b12)*(2d0*h1 + h2 + 2d0*h3)
  mat(6,2)=(1d0/30d0)*b11*(2d0*h1 + h2 + 2d0*h3)
  mat(6,3)=(1d0/30d0)*b12*(2d0*h1 + h2 + 2d0*h3)

end subroutine fem_N6i_Sc3_dN3jdx

subroutine fem_N6i_Sc6_dN6jdx(mat,h1,h2,h3,h4,h5,h6,&
  b11,b12)
use bsnqGlobVars
implicit none

  ! Integration term
  ! phi_i h d(phi_j)/dx

  real(kind=C_K2),intent(out)::mat(6,6)
  real(kind=C_K2),intent(in)::h1,h2,h3,h4,h5,h6
  real(kind=C_K2),intent(in)::b11,b12

  mat(1,1)=(-(1d0/840d0))*(b11 + b12)*(26d0*h1 - 3d0*h2 &
    - 3d0*h3 + 16d0*h4 + 4d0*h5 + 16d0*h6)
  mat(1,2)=-((b11*(18d0*h1 + 9d0*h2 - 11d0*h3 + 32d0*h4 &
    + 20d0*h5 + 16d0*h6))/2520d0)
  mat(1,3)=-((b12*(18d0*h1 - 11d0*h2 + 9d0*h3 + 16d0*h4 &
    + 20d0*h5 + 32d0*h6))/2520d0)
  mat(1,4)=(1d0/630d0)*(b12*(-6d0*h1 + 4d0*h2 - h3 &
    + 8d0*h4 + 12d0*h5 + 4d0*h6) + b11*(24d0*h1 - 5d0*h3 &
    + 20d0*h4 + 8d0*h5 + 16d0*h6))
  mat(1,5)=(1d0/630d0)*(b12*(6d0*h1 - 4d0*h2 + h3 - 8d0*h4 &
    - 12d0*h5 - 4d0*h6) + b11*(6d0*h1 + h2 &
    - 4d0*(h3 + h4 + 3d0*h5 + 2d0*h6)))
  mat(1,6)=(1d0/630d0)*(b12*(24d0*h1 - 5d0*h2 + 16d0*h4 &
    + 8d0*h5 + 20d0*h6) - b11*(6d0*h1 + h2 &
    - 4d0*(h3 + h4 + 3d0*h5 + 2d0*h6)))

  mat(2,1)=((b11 + b12)*(9d0*h1 + 18d0*h2 - 11d0*h3 &
    + 32d0*h4 + 16d0*h5 + 20d0*h6))/2520d0
  mat(2,2)=(1d0/840d0)*b11*(-3d0*h1 + 26d0*h2 - 3d0*h3 &
    + 16d0*h4 + 16d0*h5 + 4d0*h6)
  mat(2,3)=(b12*(11d0*h1 - 18d0*h2 - 9d0*h3 - 16d0*h4 &
    - 32d0*h5 - 20d0*h6))/2520d0
  mat(2,4)=(1d0/630d0)*(4d0*b12*h1 - 30d0*b12*h2 &
    + 4d0*b12*(h3 - 3d0*h4 - 3d0*h5 + h6) &
    - b11*(24d0*h2 - 5d0*h3 + 20d0*h4 + 16d0*h5 + 8d0*h6))
  mat(2,5)=(1d0/630d0)*(2d0*b12*(-2d0*h1 + 15d0*h2 &
    - 2d0*(h3 - 3d0*h4 - 3d0*h5 + h6)) &
    + b11*(h1 + 6d0*h2 - 4d0*(h3 + h4 + 2*h5 + 3*h6)))
  mat(2,6)=(1d0/630d0)*(b12*(-5d0*h1 + 5d0*h3 &
    - 4d0*h4 + 4d0*h5) + b11*(-h1 - 6d0*h2 &
    + 4d0*(h3 + h4 + 2d0*h5 + 3d0*h6)))

  mat(3,1)=((b11 + b12)*(9d0*h1 - 11d0*h2 + 18d0*h3 &
    + 20d0*h4 + 16d0*h5 + 32d0*h6))/2520d0
  mat(3,2)=(b11*(11d0*h1 - 9d0*h2 - 2d0*(9d0*h3 &
    + 10d0*h4 + 16d0*h5 + 8d0*h6)))/2520d0
  mat(3,3)=(1d0/840d0)*b12*(-3d0*h1 - 3d0*h2 + 26d0*h3 &
    + 4d0*h4 + 16d0*h5 + 16d0*h6)
  mat(3,4)=(1d0/630d0)*(b11*(-5d0*h1 + 5d0*h2 + 4d0*h5 &
    - 4d0*h6) + b12*(-h1 + 4d0*h2 - 6d0*h3 + 12d0*h4 &
    + 8d0*h5 + 4d0*h6))
  mat(3,5)=(1d0/630d0)*(-2d0*b11*(2d0*h1 + 2d0*h2 - 15d0*h3 &
    + 2d0*h4 - 6d0*h5 - 6d0*h6) + b12*(h1 &
    - 2d0*(2d0*h2 - 3d0*h3 + 6d0*h4 + 4d0*h5 + 2d0*h6)))
  mat(3,6)=(1d0/630d0)*(5d0*b12*h2 + 2d0*b11*(2d0*h1 &
    + 2d0*h2 - 15d0*h3 + 2d0*h4 - 6d0*h5 - 6d0*h6) &
    - 4d0*b12*(6d0*h3 + 2d0*h4 + 4d0*h5 + 5d0*h6))

  mat(4,1)=(-(1d0/630d0))*(b11 + b12)*(12d0*h1 - 8d0*h2 &
    - 5d0*h3 + 40d0*h4 + 4d0*h5 + 20d0*h6)
  mat(4,2)=(1d0/630d0)*b11*(-8d0*h1 + 12d0*h2 - 5d0*h3 &
    + 40d0*h4 + 20d0*h5 + 4d0*h6)
  mat(4,3)=(1d0/630d0)*b12*(-4d0*h1 - 4d0*h2 + 3d0*h3 &
    - 24d0*h4 + 4d0*h5 + 4d0*h6)
  mat(4,4)=(2d0/315d0)*(b12*(2d0*h1 - 3d0*h2 + 3d0*h3 &
    - 24d0*h4 - 12d0*h5 - 8d0*h6) &
    + b11*(5d0*h1 - 5d0*h2 - 4d0*h5 + 4d0*h6))
  mat(4,5)=(-(2d0/315d0))*(b12*(2d0*h1 - 3d0*h2 + 3d0*h3 &
    - 24d0*h4 - 12d0*h5 - 8d0*h6) + b11*(h1 + h2 + h3 &
    - 8d0*h4 - 8d0*h5 - 8d0*h6))
  mat(4,6)=(2d0/315d0)*(b12*(4d0*h1 - h2 - 2d0*h3 + 16d0*h4 &
    + 4d0*h6) + b11*(h1 + h2 + h3 - 8d0*(h4 + h5 + h6)))

  mat(5,1)=(1d0/630d0)*(b11 + b12)*(-3d0*h1 &
    + 4d0*(h2 + h3 - h4 + 6d0*h5 - h6))
  mat(5,2)=(1d0/630d0)*b11*(-5d0*h1 &
    + 4d0*(3d0*h2 - 2d0*h3 + 5d0*h4 + 10d0*h5 + h6))
  mat(5,3)=(1d0/630d0)*b12*(-5d0*h1 + 4d0*(-2d0*h2 &
    + 3d0*h3 + h4 + 10d0*h5 + 5d0*h6))
  mat(5,4)=(2d0/315d0)*(b11*(2d0*h1 - 4d0*h2 + h3 - 4d0*h4 &
    - 16d0*h5) + b12*(3d0*h1 - 3d0*h2 + 2d0*h3 - 12d0*h4 &
    - 24d0*h5 - 8d0*h6))
  mat(5,5)=(1d0/315d0)*(2d0*b12*(-3d0*h1 + 3d0*h2 &
    - 2d0*h3 + 12d0*h4 + 24d0*h5 + 8d0*h6) &
    + b11*(-6d0*h1 - 4d0*h2 + 6d0*h3 + 16d0*h4 &
    + 48d0*h5 + 24d0*h6))
  mat(5,6)=(2d0/315d0)*(b11*(3d0*h1 + 2d0*h2 - 3d0*h3 &
    - 8d0*h4 - 24d0*h5 - 12d0*h6) + b12*(2d0*h1 + h2 &
    - 4d0*(h3 + 4d0*h5 + h6)))

  mat(6,1)=(-(1d0/630d0))*(b11 + b12)*(12d0*h1 - 5d0*h2 &
    + 4d0*(-2d0*h3 + 5d0*h4 + h5 + 10d0*h6))
  mat(6,2)=(1d0/630d0)*b11*(-4d0*h1 + 3d0*h2 &
    + 4d0*(-h3 + h4 + h5 - 6d0*h6))
  mat(6,3)=(1d0/630d0)*b12*(-8d0*h1 - 5d0*h2 &
    + 4d0*(3d0*h3 + h4 + 5d0*h5 + 10d0*h6))
  mat(6,4)=(2d0/315d0)*(b11*(4d0*h1 - 2d0*h2 - h3 &
    + 4d0*h4 + 16d0*h6) + b12*(h1 + h2 + h3 - 8d0*(h4 + h5 + h6)))
  mat(6,5)=(-(2d0/315d0))*(b11*(2d0*h1 + 3d0*h2 - 3d0*h3 &
    - 8d0*h4 - 12d0*h5 - 24d0*h6) + b12*(h1 + h2 &
    + h3 - 8d0*h4 - 8d0*h5 - 8d0*h6))
  mat(6,6)=(2d0/315d0)*(b12*(5d0*h1 - 5d0*h3 + 4d0*h4 &
    - 4d0*h5) + b11*(2d0*h1 + 3d0*h2 - 3d0*h3 - 8d0*h4 &
    - 12d0*h5 - 24d0*h6))

end subroutine fem_N6i_Sc6_dN6jdx



subroutine fem_N6i_du6N6jdx(mat,u1,u2,u3,u4,u5,u6,b11,b12)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::mat(6,6)
  real(kind=C_K2),intent(in)::u1,u2,u3,u4,u5,u6,b11,b12

  mat(1,1)=(1d0/840d0)*b11*(-52d0*u1 - 3d0*u2 + 3d0*u3 + 16d0*u4 &
    + 4d0*u5 - 24d0*u6) + (1d0/840d0)*b12*(-52d0*u1 + 3d0*u2 &
    - 3d0*u3 - 24d0*u4 + 4d0*u5 + 16d0*u6)

  mat(1,2)=(b11*(-9d0*u1 - 18d0*u2 + 11d0*u3 - 32d0*u4 - 16d0*u5 &
    - 20d0*u6))/2520d0 + (b12*(9d0*u1 + 11d0*u3 + 16d0*u4 &
    - 16d0*u5 - 20d0*u6))/2520d0

  mat(1,3)=(b12*(-9d0*u1 + 11d0*u2 - 18d0*u3 - 20d0*u4 - 16d0*u5 &
    - 32d0*u6))/2520d0 + (b11*(9d0*u1 + 11d0*u2 - 20d0*u4 &
    - 16d0*u5 + 16d0*u6))/2520d0

  mat(1,4)=(1d0/630d0)*b12*(-18d0*u1 + 4d0*u2 - 5d0*u3 + 16d0*u4 &
    + 4d0*u5 + 20d0*u6) + (1d0/630d0)*b11*(12d0*u1 - 8d0*u2 &
    - 5d0*u3 + 40d0*u4 + 4d0*u5 + 20d0*u6)

  mat(1,5)=(1d0/630d0)*b11*(3d0*u1 - 4d0*(u2 + u3 - u4 + 6d0*u5 &
    - u6)) + (1d0/630d0)*b12*(3d0*u1 &
    - 4d0*(u2 + u3 - u4 + 6d0*u5 - u6))

  mat(1,6)=(1d0/630d0)*b11*(-18d0*u1 - 5d0*u2 &
    + 4d0*(u3 + 5d0*u4 + u5 + 4d0*u6)) &
    + (1d0/630d0)*b12*(12d0*u1 - 5d0*u2 &
    + 4d0*(-2d0*u3 + 5d0*u4 + u5 + 10d0*u6))


  mat(2,1)=(1d0/420d0)*b12*(3d0*u1 + 3d0*u2 + 8d0*u4) &
    + (b11*(18d0*u1 + 9d0*u2 - 11d0*u3 + 32d0*u4 &
    + 20d0*u5 + 16d0*u6))/2520d0

  mat(2,2)=(1d0/420d0)*b12*(3d0*u1 - 3d0*u3 - 20d0*u4 &
    + 20d0*u5) + (1d0/840d0)*b11*(3d0*u1 + 52d0*u2 - 3d0*u3 &
    - 16d0*u4 + 24d0*u5 - 4d0*u6)

  mat(2,3)=(-(1d0/420d0))*b12*(3d0*u2 + 3d0*u3 + 8d0*u5) &
    + (b11*(-11d0*u1 - 9d0*u2 + 20d0*u4 &
    - 16d0*u5 + 16d0*u6))/2520d0

  mat(2,4)=(1d0/105d0)*b12*(2d0*u1 - 5d0*u2 - 4d0*u4) &
    + (1d0/630d0)*b11*(8d0*u1 - 12d0*u2 + 5d0*u3 &
    - 40d0*u4 - 20d0*u5 - 4d0*u6)

  mat(2,5)=(1d0/105d0)*b12*(5d0*u2 - 2d0*u3 + 4d0*u5) &
    + (1d0/630d0)*b11*(5d0*u1 + 18d0*u2 &
    - 4d0*(u3 + 5d0*u4 + 4d0*u5 + u6))

  mat(2,6)=(1d0/630d0)*b11*(4d0*u1 - 3d0*u2 &
    + 4d0*(u3 - u4 - u5 + 6d0*u6))


  mat(3,1)=(1d0/420d0)*b11*(3d0*u1 + 3d0*u3 + 8d0*u6) &
    + (b12*(18d0*u1 - 11d0*u2 + 9d0*u3 + 16d0*u4 &
    + 20d0*u5 + 32d0*u6))/2520d0

  mat(3,2)=(-(1d0/420d0))*b11*(3d0*u2 + 3d0*u3 + 8d0*u5) &
    + (b12*(-11d0*u1 - 9d0*u3 + 16d0*u4 &
    - 16d0*u5 + 20d0*u6))/2520d0

  mat(3,3)=(1d0/840d0)*b11*(6d0*u1 - 6d0*u2 + 40d0*u5 &
    - 40d0*u6) + (1d0/840d0)*b12*(3d0*u1 - 3d0*u2 &
    + 52d0*u3 - 4d0*u4 + 24d0*u5 - 16d0*u6)

  mat(3,4)=(1d0/630d0)*b12*(4d0*u1 + 4d0*u2 - 3d0*u3 &
    + 24d0*u4 - 4d0*u5 - 4d0*u6)
 
  mat(3,5)=(1d0/105d0)*b11*(-2d0*u2 + 5d0*u3 + 4d0*u5) &
    + (1d0/630d0)*b12*(5d0*u1 &
    - 2d0*(2d0*u2 - 9d0*u3 + 2d0*(u4 + 4d0*u5 + 5d0*u6)))

  mat(3,6)=(1d0/105d0)*b11*(2d0*u1 - 5d0*u3 - 4d0*u6) &
    + (1d0/630d0)*b12*(8d0*u1 + 5d0*u2 &
    - 4d0*(3d0*u3 + u4 + 5d0*u5 + 10d0*u6))


  mat(4,1)=(1d0/630d0)*b11*(-24d0*u1 + 5d0*u3 &
    - 20d0*u4 - 8d0*u5 - 16d0*u6) &
    + (1d0/630d0)*b12*(-24d0*u1 + 8d0*u2 + u3 &
    - 32d0*u4 - 12d0*u5 - 4d0*u6)

  mat(4,2)=(2d0/315d0)*b12*(2d0*u1 - u3 - 3d0*u4 &
    + 3d0*u5 - u6) + (1d0/630d0)*b11*(24d0*u2 - 5d0*u3 &
    + 20d0*u4 + 16d0*u5 + 8d0*u6)

  mat(4,3)=(1d0/630d0)*b11*(5d0*u1 - 5d0*u2 - 4d0*u5 &
    + 4d0*u6) + (1d0/630d0)*b12*(u1 &
    - 2d0*(2d0*u2 - 3d0*u3 + 6d0*u4 + 4d0*u5 + 2d0*u6))

  mat(4,4)=(-(2d0/315d0))*b12*(8d0*u1 + 3d0*u2 + 3d0*u3 &
    + 48d0*u4 - 12d0*u5 - 8d0*u6) &
    - (2d0/315d0)*b11*(5d0*u1 - 5d0*u2 - 4d0*u5 + 4d0*u6)

  mat(4,5)=(1d0/315d0)*b11*(-4d0*u1 + 8d0*u2 - 2d0*u3 &
    + 8d0*u4 + 32d0*u5) + (2d0/315d0)*b12*(-3d0*u1 &
    + 3d0*u2 - 2d0*u3 + 12d0*u4 + 24d0*u5 + 8d0*u6)

  mat(4,6)=(-(2d0/315d0))*b12*(u1 + u2 + u3 - 8d0*u4 &
    - 8d0*u5 - 8d0*u6) - (2d0/315d0)*b11*(4d0*u1 &
    - 2d0*u2 - u3 + 4d0*u4 + 16d0*u6)


  mat(5,1)=(1d0/630d0)*b12*(-6d0*u1 + 4d0*u2 - u3 &
    + 8d0*u4 + 12d0*u5 + 4d0*u6) + (1d0/630d0)*b11*(-6d0*u1 &
    - u2 + 4d0*(u3 + u4 + 3d0*u5 + 2d0*u6))

  mat(5,2)=(2d0/315d0)*b12*(u1 - 2d0*u3 - 3d0*u4 + 3d0*u5 &
    + u6) + (1d0/630d0)*b11*(-u1 &
    + 4d0*(6d0*u2 - 2d0*u3 + u4 + 8d0*u5 + 3d0*u6))

  mat(5,3)=(2d0/315d0)*b11*(u1 - 2d0*u2 + u4 + 3d0*u5 &
    - 3d0*u6) + (1d0/630d0)*b12*(-u1 + 4d0*(-2d0*u2 &
    + 6d0*u3 + 3d0*u4 + 8d0*u5 + u6))

  mat(5,4)=(2d0/315d0)*b12*(2d0*u1 - 3d0*u2 + 3d0*u3 &
    - 24d0*u4 - 12d0*u5 - 8d0*u6) &
    + (2d0/315d0)*b11*(u1 + u2 + u3 &
    - 8d0*u4 - 8d0*u5 - 8d0*u6)

  mat(5,5)=(2d0/315d0)*b11*(3d0*u1 + 8d0*u2 + 3d0*u3 &
    - 8d0*u4 + 48d0*u5 - 12d0*u6) + (2d0/315d0)*b12*(3d0*u1 &
    + 3d0*u2 + 8d0*u3 - 12d0*u4 + 48d0*u5 - 8d0*u6)

  mat(5,6)=(2d0/315d0)*b11*(2d0*u1 + 3d0*u2 - 3d0*u3 &
    - 8d0*u4 - 12d0*u5 - 24d0*u6) + (2d0/315d0)*b12*(u1 &
    + u2 + u3 - 8d0*u4 - 8d0*u5 - 8d0*u6)


  mat(6,1)=(1d0/630d0)*b11*(-24d0*u1 + u2 + 8d0*u3 - 4d0*u4 &
    - 12d0*u5 - 32d0*u6) + (1d0/630d0)*b12*(-24d0*u1 &
    + 5d0*u2 - 16d0*u4 - 8d0*u5 - 20d0*u6)

  mat(6,2)=(1d0/630d0)*b12*(5d0*u1 - 5d0*u3 + 4d0*u4 - 4d0*u5) &
    + (1d0/630d0)*b11*(u1 + 6d0*u2 &
    - 4d0*(u3 + u4 + 2d0*u5 + 3d0*u6))

  mat(6,3)=(2d0/315d0)*b11*(2d0*u1 - u2 - u4 &
    + 3d0*u5 - 3d0*u6) + (1d0/630d0)*b12*(-5d0*u2 &
    + 4d0*(6d0*u3 + 2d0*u4 + 4d0*u5 + 5d0*u6))

  mat(6,4)=(-(2d0/315d0))*b11*(u1 + u2 + u3 - 8d0*u4 &
    - 8d0*u5 - 8d0*u6) - (2d0/315d0)*b12*(4d0*u1 - u2 &
    - 2d0*u3 + 16d0*u4 + 4d0*u6)

  mat(6,5)=(-(2d0/315d0))*b11*(3d0*u1 + 2d0*u2 - 3d0*u3 &
    - 8d0*u4 - 24d0*u5 - 12d0*u6) &
    - (2d0/315d0)*b12*(2d0*u1 + u2 - 4d0*(u3 + 4d0*u5 + u6))

  mat(6,6)=(-(2d0/315d0))*b12*(5d0*u1 - 5d0*u3 + 4d0*u4 &
    - 4d0*u5) - (2d0/315d0)*b11*(8d0*u1 + 3d0*u2 + 3d0*u3 &
    - 8d0*u4 - 12d0*u5 + 48d0*u6)

end subroutine fem_N6i_du6N6jdx



subroutine fem_dN6iSc6dx_N6j(mat,h1,h2,h3,h4,h5,h6,&
  b11,b12)
use bsnqGlobVars
implicit none

  ! Integration term
  ! d(phi_i h)/dx phi_j

  real(kind=C_K2),intent(out)::mat(6,6)
  real(kind=C_K2),intent(in)::h1,h2,h3,h4,h5,h6
  real(kind=C_K2),intent(in)::b11,b12


  mat(1,1) = (1d0/840d0)*(b11*(-52d0*h1 - 3d0*h2 + 3d0*h3 + 16d0*h4 &
    + 4d0*h5 - 24d0*h6) + b12*(-52d0*h1 + 3d0*h2 - 3d0*h3 - 24d0*h4 &
    + 4d0*h5 + 16d0*h6))

  mat(1,2) = (6d0*b12*(3d0*h1 + 3d0*h2 + 8d0*h4) + b11*(18d0*h1 &
    + 9d0*h2 - 11d0*h3 + 32d0*h4 + 20d0*h5 + 16d0*h6))/2520d0

  mat(1,3) = (6d0*b11*(3d0*h1 + 3d0*h3 + 8d0*h6) + b12*(18d0*h1 &
    - 11d0*h2 + 9d0*h3 + 16d0*h4 + 20d0*h5 + 32d0*h6))/2520d0

  mat(1,4) = (1d0/630d0)*(b12*(-24d0*h1 + 8d0*h2 + h3 - 32d0*h4 &
    - 12d0*h5 - 4d0*h6) - b11*(24d0*h1 - 5d0*h3 + 20d0*h4 + 8d0*h5 &
    + 16d0*h6))

  mat(1,5) = (1d0/630d0)*(b12*(-6d0*h1 + 4d0*h2 - h3 + 8d0*h4 &
    + 12d0*h5 + 4d0*h6) + b11*(-6d0*h1 - h2 &
    + 4d0*(h3 + h4 + 3d0*h5 + 2d0*h6)))

  mat(1,6) = (1d0/630d0)*(b11*(-24d0*h1 + h2 + 8d0*h3 - 4d0*h4 &
    - 12d0*h5 - 32d0*h6) - b12*(24d0*h1 - 5d0*h2 + 16d0*h4 &
    + 8d0*h5 + 20d0*h6))


  mat(2,1) = (b12*(9d0*h1 + 11d0*h3 + 16d0*h4 - 16d0*h5 &
    - 20d0*h6) - b11*(9d0*h1 + 18d0*h2 - 11d0*h3 + 32d0*h4 &
    + 16d0*h5 + 20d0*h6))/2520d0

  mat(2,2) = (1d0/840d0)*(2d0*b12*(3d0*h1 - 3d0*h3 - 20d0*h4 &
    + 20d0*h5) + b11*(3d0*h1 + 52d0*h2 - 3d0*h3 - 16d0*h4 &
    + 24d0*h5 - 4d0*h6))

  mat(2,3) = (-6d0*b11*(3d0*h2 + 3d0*h3 + 8d0*h5) &
    + b12*(-11d0*h1 - 9d0*h3 + 16d0*h4 - 16d0*h5 + 20d0*h6))/2520d0

  mat(2,4) = (1d0/630d0)*(-4d0*b12*(-2d0*h1 + h3 + 3d0*h4 &
    - 3d0*h5 + h6) + b11*(24d0*h2 - 5d0*h3 + 20d0*h4 + 16d0*h5 &
    + 8d0*h6))

  mat(2,5) = (1d0/630d0)*((-b11)*h1 + 4d0*b12*(h1 - 2d0*h3 &
    - 3d0*h4 + 3d0*h5 + h6) + 4d0*b11*(6d0*h2 - 2d0*h3 + h4 &
    + 8d0*h5 + 3d0*h6))

  mat(2,6) = (1d0/630d0)*(b12*(5d0*h1 - 5d0*h3 + 4d0*h4 &
    - 4d0*h5) + b11*(h1 + 6d0*h2 &
    - 4d0*(h3 + h4 + 2d0*h5 + 3d0*h6)))


  mat(3,1) = (b11*(9d0*h1 + 11d0*h2 - 20d0*h4 - 16d0*h5 &
    + 16d0*h6) - b12*(9d0*h1 - 11d0*h2 + 18d0*h3 + 20d0*h4 &
    + 16d0*h5 + 32d0*h6))/2520d0

  mat(3,2) = (-6d0*b12*(3d0*h2 + 3d0*h3 + 8d0*h5) &
    + b11*(-11d0*h1 - 9d0*h2 + 20d0*h4 - 16d0*h5 + 16d0*h6))/2520d0

  mat(3,3) = (1d0/840d0)*(b11*(6d0*h1 - 6d0*h2 + 40d0*h5 &
    - 40d0*h6) + b12*(3d0*h1 - 3d0*h2 + 52d0*h3 - 4d0*h4 &
    + 24d0*h5 - 16d0*h6))

  mat(3,4) = (1d0/630d0)*(b11*(5d0*h1 - 5d0*h2 - 4d0*h5 &
    + 4d0*h6) + b12*(h1 &
    - 2d0*(2d0*h2 - 3d0*h3 + 6d0*h4 + 4d0*h5 + 2d0*h6)))

  mat(3,5) = (1d0/630d0)*((-b12)*h1 &
    + 4d0*b11*(h1 - 2d0*h2 + h4 + 3d0*h5 - 3d0*h6) &
    + 4d0*b12*(-2d0*h2 + 6d0*h3 + 3d0*h4 + 8d0*h5 + h6))

  mat(3,6) = (1d0/630d0)*(-5d0*b12*h2 &
    + 4d0*b11*(2d0*h1 - h2 - h4 + 3d0*h5 - 3d0*h6) &
    + 4d0*b12*(6d0*h3 + 2d0*h4 + 4d0*h5 + 5d0*h6))


  mat(4,1) = (1d0/630d0)*(b12*(-18d0*h1 + 4d0*h2 - 5d0*h3 &
    + 16d0*h4 + 4d0*h5 + 20d0*h6) + b11*(12d0*h1 - 8d0*h2 &
    - 5d0*h3 + 40d0*h4 + 4d0*h5 + 20d0*h6))

  mat(4,2) = (1d0/630d0)*(6d0*b12*(2d0*h1 - 5d0*h2 - 4d0*h4) &
    + b11*(8d0*h1 - 12d0*h2 + 5d0*h3 - 40d0*h4 - 20d0*h5 - 4d0*h6))

  mat(4,3) = (1d0/630d0)*b12*(4d0*h1 + 4d0*h2 - 3d0*h3 &
    + 24d0*h4 - 4d0*h5 - 4d0*h6)

  mat(4,4) = (-(2d0/315d0))*(b12*(8d0*h1 + 3d0*h2 + 3d0*h3 &
    + 48d0*h4 - 12d0*h5 - 8d0*h6) + b11*(5d0*h1 &
    - 5d0*h2 - 4d0*h5 + 4d0*h6))

  mat(4,5) = (2d0/315d0)*(b12*(2d0*h1 - 3d0*h2 + 3d0*h3 &
    - 24d0*h4 - 12d0*h5 - 8d0*h6) + b11*(h1 + h2 + h3 &
    - 8d0*h4 - 8d0*h5 - 8d0*h6))

  mat(4,6) = (2d0/315d0)*((-b11)*(h1 + h2 + h3 - 8d0*h4 &
    - 8d0*h5 - 8d0*h6) + b12*(-4d0*h1 + h2 + 2d0*h3 &
    - 16d0*h4 - 4d0*h6))


  mat(5,1) = (1d0/630d0)*(b11 + b12)*(3d0*h1 &
    - 4d0*(h2 + h3 - h4 + 6d0*h5 - h6))

  mat(5,2) = (1d0/630d0)*(6d0*b12*(5d0*h2 - 2d0*h3 + 4d0*h5) &
    + b11*(5d0*h1 + 18d0*h2 - 4d0*(h3 + 5d0*h4 + 4d0*h5 + h6)))

  mat(5,3) = (1d0/630d0)*(5d0*b12*h1 &
    + 6d0*b11*(-2d0*h2 + 5d0*h3 + 4d0*h5) &
    - 2d0*b12*(2d0*h2 - 9d0*h3 + 2d0*(h4 + 4d0*h5 + 5d0*h6)))

  mat(5,4) = (1d0/315d0)*(b11*(-4d0*h1 + 8d0*h2 - 2d0*h3 &
    + 8d0*h4 + 32d0*h5) + 2d0*b12*(-3d0*h1 + 3d0*h2 - 2d0*h3 &
    + 12d0*h4 + 24d0*h5 + 8d0*h6))

  mat(5,5) = (2d0/315d0)*(b11*(3d0*h1 + 8d0*h2 + 3d0*h3 &
    - 8d0*h4 + 48d0*h5 - 12d0*h6) + b12*(3d0*h1 + 3d0*h2 &
    + 8d0*h3 - 12d0*h4 + 48d0*h5 - 8d0*h6))

  mat(5,6) = (-(2d0/315d0))*(b11*(3d0*h1 + 2d0*h2 - 3d0*h3 &
    - 8d0*h4 - 24d0*h5 - 12d0*h6) + b12*(2d0*h1 + h2 &
    - 4d0*(h3 + 4d0*h5 + h6)))


  mat(6,1) = (1d0/630d0)*(b11*(-18d0*h1 - 5d0*h2 &
    + 4d0*(h3 + 5d0*h4 + h5 + 4d0*h6)) + b12*(12d0*h1 &
    - 5d0*h2 + 4d0*(-2d0*h3 + 5d0*h4 + h5 + 10d0*h6)))

  mat(6,2) = (1d0/630d0)*b11*(4d0*h1 - 3d0*h2 &
    + 4d0*(h3 - h4 - h5 + 6d0*h6))

  mat(6,3) = (1d0/630d0)*(6d0*b11*(2d0*h1 - 5d0*h3 - 4d0*h6) &
    + b12*(8d0*h1 + 5d0*h2 - 4d0*(3d0*h3 + h4 + 5d0*h5 + 10d0*h6)))

  mat(6,4) = (-(2d0/315d0))*(b11*(4d0*h1 - 2d0*h2 - h3 &
    + 4d0*h4 + 16d0*h6) + b12*(h1 + h2 + h3 - 8d0*(h4 + h5 + h6)))

  mat(6,5) = (2d0/315d0)*(b11*(2d0*h1 + 3d0*h2 - 3d0*h3 &
    - 8d0*h4 - 12d0*h5 - 24d0*h6) + b12*(h1 + h2 + h3 &
    - 8d0*h4 - 8d0*h5 - 8d0*h6))

  mat(6,6) = (-(2d0/315d0))*(b12*(5d0*h1 - 5d0*h3 &
    + 4d0*h4 - 4d0*h5) + b11*(8d0*h1 + 3d0*h2 + 3d0*h3 &
    - 8d0*h4 - 12d0*h5 + 48d0*h6))

end subroutine fem_dN6iSc6dx_N6j



subroutine femBnd_N3i_Sc3_N3j(mat,h1,h2,h3,side)
use bsnqGlobVars
implicit none
  
  integer(kind=C_K1),intent(in)::side
  real(kind=C_K2),intent(out)::mat(3,3)
  real(kind=C_K2),intent(in)::h1,h2,h3

  mat=0d0

  if(side.eq.1)then !S12
    mat(1,1)=h1/4d0 + h2/12d0 
    mat(1,2)=h1/12d0 + h2/12d0

    mat(2,1)=h1/12d0 + h2/12d0
    mat(2,2)=h1/12d0 + h2/4d0 

  
  elseif(side.eq.2)then !S23
    mat(2,2)=h2/4d0 + h3/12d0
    mat(2,3)=h2/12d0 + h3/12d0

    mat(3,2)=h2/12d0 + h3/12d0 
    mat(3,3)=h2/12d0 + h3/4d0

  
  elseif(side.eq.3)then !S31
    mat(1,1)=h1/4d0 + h3/12d0
    mat(1,3)=h1/12d0 + h3/12d0

    mat(3,1)=h1/12d0 + h3/12d0
    mat(3,3)=h1/12d0 + h3/4d0


  else
    write(*,'(" [ERR] Wrong side number in femBnd_N3i_Sc3_N3j")')
    stop
  endif

end subroutine femBnd_N3i_Sc3_N3j



subroutine femBnd_N3i_Sc3_dN3jdx(mat,h1,h2,h3,b11,b12,side)
use bsnqGlobVars
implicit none
  
  integer(kind=C_K1),intent(in)::side
  real(kind=C_K2),intent(out)::mat(3,3)
  real(kind=C_K2),intent(in)::h1,h2,h3,b11,b12

  mat=0d0

  if(side.eq.1)then !S12
    mat(1,1)=-((b11*h1)/3d0) - (b12*h1)/3d0 &
      - (b11*h2)/6d0 - (b12*h2)/6d0 
    mat(1,2)=(b11*h1)/3d0 + (b11*h2)/6d0
    mat(1,3)=(b12*h1)/3d0 + (b12*h2)/6d0

    mat(2,1)=-((b11*h1)/6d0) - (b12*h1)/6d0 &
      - (b11*h2)/3d0 - (b12*h2)/3d0
    mat(2,2)=(b11*h1)/6d0 + (b11*h2)/3d0
    mat(2,3)=(b12*h1)/6d0 + (b12*h2)/3d0

  
  elseif(side.eq.2)then !S23
    mat(2,1)=-((b11*h2)/3d0) - (b12*h2)/3d0 &
      - (b11*h3)/6d0 - (b12*h3)/6d0
    mat(2,2)=(b11*h2)/3d0 + (b11*h3)/6d0
    mat(2,3)=(b12*h2)/3d0 + (b12*h3)/6d0

    mat(3,1)=-((b11*h2)/6d0) - (b12*h2)/6d0 &
      - (b11*h3)/3d0 - (b12*h3)/3d0
    mat(3,2)=(b11*h2)/6d0 + (b11*h3)/3d0
    mat(3,3)=(b12*h2)/6d0 + (b12*h3)/3d0

  
  elseif(side.eq.3)then !S31
    mat(1,1)=-((b11*h1)/3d0) - (b12*h1)/3d0 &
      - (b11*h3)/6d0 - (b12*h3)/6d0
    mat(1,2)=(b11*h1)/3d0 + (b11*h3)/6d0
    mat(1,3)=(b12*h1)/3d0 + (b12*h3)/6d0

    mat(3,1)=-((b11*h1)/6d0) - (b12*h1)/6d0 &
      - (b11*h3)/3d0 - (b12*h3)/3d0 
    mat(3,2)=(b11*h1)/6d0 + (b11*h3)/3d0
    mat(3,3)=(b12*h1)/6d0 + (b12*h3)/3d0


  else
    write(*,'(" [ERR] Wrong side number in femBnd_N3i_Sc3_N3j")')
    stop
  endif

end subroutine femBnd_N3i_Sc3_dN3jdx