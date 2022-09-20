!!---------------------------vertVelMod----------------------------!!
module vertVelMod
use bsnqGlobVars
implicit none

  type, public :: vertVelDerv      
    real(kind=C_K2)::u,ux,uxx,uxxx
    real(kind=C_K2)::uh,uhx,uhxx,uhxxx
    real(kind=C_K2)::hx
    ! d2(U)/dxdt and d2(Uh)/dxdt
    real(kind=C_K2)::uxt,uhxt      
    ! d(U)/dx and d(Uh)/dx at t(n-1), t(n-2)   
    real(kind=C_K2)::ux_tn(2),uhx_tn(2)    
  contains
  !   procedure ::  init => initVertVelDerv
    procedure ::  initByInterp
  end type vertVelDerv



  type, public :: vertVelDerv3D

    ! P = uh , note its not u(h+eta)
    real(kind=C_K2)::u, v, P, Q
    ! ug = del.(u i + v j) . Pg = del.(uh i + vh j)
    real(kind=C_K2)::ug, Pg
    ! ugx = d( del.(u) ) / dx
    real(kind=C_K2)::ugx, ugy, Pgx, Pgy
    ! uggg = del.( del( del.(u) ) )
    real(kind=C_K2)::uggg, Pggg
    ! hx = d(h)/dx 
    real(kind=C_K2):: hx, hy
    ! ugt = d( del.(u) )/dt
    real(kind=C_K2)::ugt, Pgt
    ! del.u and del.P at t(n-1), t(n-2)   
    real(kind=C_K2)::ug_tn(2),Pg_tn(2)    
  
  contains  

    procedure ::  initByInterp3D

  end type vertVelDerv3D



  type, public :: vertVelProbes
    integer(kind=C_K1)::np, fileid
    real(kind=C_K2),allocatable::x(:), y(:), z(:), wrkr(:,:)
    real(kind=C_K2),allocatable::u(:), v(:), w(:), p(:), eta(:)
    integer(kind=C_K1),allocatable::wrki(:), err(:)
  contains
    procedure :: initvvProbes    
  end type  vertVelProbes

contains


!!------------------------initByInterp-------------------------!!
  subroutine initByInterp(b,nNei,nei,wei,i)
  implicit none

    class(vertVelDerv),intent(out)::b
    integer(kind=C_K1),intent(in)::nNei
    real(kind=C_K2),intent(in)::wei(nNei)
    type(vertVelDerv),intent(in)::nei(nNei)

    integer(kind=C_K1),intent(out)::i

    b%u=0d0
    b%ux=0d0
    b%uxx=0d0
    b%uxxx=0d0
    b%uh=0d0
    b%uhx=0d0
    b%uhxx=0d0
    b%uhxxx=0d0
    b%hx=0d0
    b%uxt=0d0
    b%uhxt=0d0    

    do i=1,nNei
      b%u = b%u + wei(i) * nei(i)%u
      b%ux =  b%ux + wei(i) * nei(i)%ux
      b%uxx = b%uxx +  wei(i) * nei(i)%uxx
      b%uxxx =  b%uxxx + wei(i) * nei(i)%uxxx
      b%uh =  b%uh + wei(i) * nei(i)%uh
      b%uhx = b%uhx +  wei(i) * nei(i)%uhx
      b%uhxx =  b%uhxx + wei(i) * nei(i)%uhxx
      b%uhxxx = b%uhxxx +  wei(i) * nei(i)%uhxxx
      b%hx =  b%hx + wei(i) * nei(i)%hx
      b%uxt = b%uxt +  wei(i) * nei(i)%uxt
      b%uhxt =  b%uhxt + wei(i) * nei(i)%uhxt      
    enddo

  end subroutine initByInterp
!!----------------------End initByInterp-----------------------!!



!!-----------------------initByInterp3D------------------------!!
  subroutine initByInterp3D(b,nNei,nei,wei,i)
  implicit none

    class(vertVelDerv3D),intent(out)::b
    integer(kind=C_K1),intent(in)::nNei
    real(kind=C_K2),intent(in)::wei(nNei)
    type(vertVelDerv3D),intent(in)::nei(nNei)

    integer(kind=C_K1),intent(out)::i


    b%u = 0d0
    b%v = 0d0
    b%P = 0d0
    b%Q = 0d0

    b%ug = 0d0
    b%Pg = 0d0
    b%hx = 0d0
    b%hy = 0d0

    b%ugx = 0d0
    b%ugy = 0d0
    b%Pgx = 0d0
    b%Pgy = 0d0

    b%uggg = 0d0
    b%Pggg = 0d0
    b%ugt = 0d0
    b%Pgt = 0d0
    
    do i = 1, nNei
      b%u = b%u + wei(i) * nei(i)%u
      b%v = b%v + wei(i) * nei(i)%v
      b%P = b%P + wei(i) * nei(i)%P
      b%Q = b%Q + wei(i) * nei(i)%Q

      b%ug = b%ug + wei(i) * nei(i)%ug
      b%Pg = b%Pg + wei(i) * nei(i)%Pg
      b%hx = b%hx + wei(i) * nei(i)%hx
      b%hy = b%hy + wei(i) * nei(i)%hy

      b%ugx = b%ugx + wei(i) * nei(i)%ugx
      b%ugy = b%ugy + wei(i) * nei(i)%ugy
      b%Pgx = b%Pgx + wei(i) * nei(i)%Pgx
      b%Pgy = b%Pgy + wei(i) * nei(i)%Pgy

      b%uggg = b%uggg + wei(i) * nei(i)%uggg
      b%Pggg = b%Pggg + wei(i) * nei(i)%Pggg
      b%ugt = b%ugt + wei(i) * nei(i)%ugt
      b%Pgt = b%Pgt + wei(i) * nei(i)%Pgt
    enddo

  end subroutine initByInterp3D
!!---------------------End initByInterp3D----------------------!!



!!-------------------------vertVelExp--------------------------!!
  subroutine vertVelExp(z,h,eta,hx,u,ux,uxx,uxxx,uhx,uhxx,uhxxx,&
    uxt,uhxt,uc,wc,pc)
  implicit none

    real(kind=C_K2),intent(in)::z,h,eta,hx,u,ux,uxx,uxxx
    real(kind=C_K2),intent(in)::uhx,uhxx,uhxxx,uxt,uhxt
    real(kind=C_K2),intent(out)::uc,wc,pc

    uc = u - (0.5d0*h*uhxx - h*h*uxx/6d0) - (z*uhxx + 0.5d0*z*z*uxx)
    wc = -uhx - z*ux + z/2d0*( hx*uhxx + h*uhxxx ) &
      - z/6d0*( 2d0*h*hx*uxx + h*h*uxxx ) &
      + z*z/2d0*uhxxx + z*z*z/6d0*uxxx
    pc = rhoW * ( grav*( -z + eta ) + z*(uhxt + 0.5d0*z*uxt) )    

  end subroutine vertVelExp
!!-----------------------End vertVelExp------------------------!!



!!------------------------vertVelExp3D-------------------------!!
  subroutine vertVelExp3D(z, h, eta, hx, hy, u, v, P, Q,&
    ug, Pg, ugx, ugy, Pgx, Pgy, uggg, Pggg, ugt, Pgt, &
    uc, vc, wc, pc)
  implicit none

    real(kind=C_K2),intent(in)::z, h, eta, hx, hy, u, v, P, Q
    real(kind=C_K2),intent(in)::ug, Pg, ugx, ugy, Pgx, Pgy
    real(kind=C_K2),intent(in)::uggg, Pggg, ugt, Pgt    
    real(kind=C_K2),intent(out)::uc, vc, wc, pc


    uc = u + (h*h/6d0 - z*z/2d0)*ugx - (h/2d0 + z)*Pgx

    vc = v + (h*h/6d0 - z*z/2d0)*ugy - (h/2d0 + z)*Pgy

    wc = -Pg - z*ug - z*h/3d0*(hx*ugx + hy*ugy) &
      + z/2d0*(hx*Pgx + hy*Pgy) &
      - z/6d0*(h*h - z*z)*uggg &
      + z/2d0*(h + z)*Pggg
    
    pc = rhoW * ( grav*(eta - z) + z*(Pgt + 0.5d0*z*ugt) )    

  end subroutine vertVelExp3D
!!----------------------End vertVelExp3D-----------------------!!



!!------------------------initvvProbes-------------------------!!
  subroutine initvvProbes(b,np)
  implicit none

    class(vertVelProbes),intent(inout)::b
    integer(kind=C_K1),intent(in)::np

    b%np = np
    allocate(b%x(np), b%y(np), b%z(np), b%p(np), b%wrkr(np,2))
    allocate(b%u(np), b%v(np), b%w(np), b%wrki(np), b%err(np))
    allocate(b%eta(np))

    b%x = 0d0
    b%y = 0d0
    b%z = 0d0
    b%err = 0

  end subroutine initvvProbes
!!----------------------End initvvProbes-----------------------!!


end module vertVelMod
!!---------------------------End shipMod---------------------------!!