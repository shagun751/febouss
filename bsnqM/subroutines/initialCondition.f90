!! Apply initial conditions inside the subroutine initMat() inside 
!! bsnqModule.f90

!!-----------------------------solitIC-----------------------------!!
  subroutine solitIC(npl,npt,cor,er,pr,qr)
  use bsnqGlobVars
  implicit none
    
    integer(kind=C_K1),intent(in)::npl,npt
    real(kind=C_K2),intent(in)::cor(npt,2)
    real(kind=C_K2),intent(out)::er(npl),pr(npt),qr(npt)

    integer(kind=C_K1)::i
    real(kind=C_K2)::tmpr1,tmpr2,tmpr3,tmpr4,tmpr5

    ! pr=0d0
    ! qr=0d0
    ! !er=0.045d0*dexp(-2d0*( (cor(1:npl,1)-18.288d0)**2 ))
    ! do i=1,npl
    !   er(i)=(cor(i,1)-18.288d0)**2
    !   if(er(i).gt.25)then
    !     er(i)=0d0
    !   else
    !     er(i)=0.045*exp(-2d0*er(i))
    !   endif
    ! enddo

    er=0d0
    pr=0d0
    qr=0d0
    tmpr1=0.45d0
    tmpr2=0.045d0
    tmpr3=dsqrt(grav*(tmpr1+tmpr2))
    tmpr4=dsqrt(3*tmpr2/(4*(tmpr1**3)))
    do i=1,npt
      if((cor(i,1).ge.3d0).and.(cor(i,1).le.19d0)) then
        tmpr5=tmpr2/(dcosh(tmpr4*(cor(i,1)-(11d0)))**2)
        pr(i)=tmpr3*tmpr5
        if(i.le.npl) then
          er(i)=tmpr5
        endif
      endif
    enddo  

  end subroutine solitIC
!!---------------------------End solitIC---------------------------!!



!!-------------------------solitICFromFile-------------------------!!
  subroutine solitICFromFile(npl,npt,cor,er,pr,qr)
  use bsnqGlobVars
  implicit none
    
    integer(kind=C_K1),intent(in)::npl,npt
    real(kind=C_K2),intent(in)::cor(npt,2)
    real(kind=C_K2),intent(out)::er(npl),pr(npt),qr(npt)

    integer(kind=C_K1)::i, ff, nn, j
    real(kind=C_K2)::tmpr1, x, dx
    real(kind=C_K2),allocatable::data(:,:)
    character(len=C_KSTR)::fileName

    er = 0d0
    pr = 0d0
    qr = 0d0
  
    write(9,'(" [INF] Solitary from file")')    

    !Opening mesh file  
    fileName='solitaryWave_dx0010.dat'    
    open(newunit=ff,file=trim(fileName))       

    read(ff,*,end=11,err=11)nn
    allocate(data(nn,3))

    do i = 1, nn
      read(ff,*,end=11,err=11) data(i,1:3)
    enddo

    dx = data(2,1) - data(1,1)
    dx = dx/10d0

    do i = 1, npt
      x = cor(i,1)

      j=1
      do while(.true.)        
        if(abs(data(j,1)-x).lt.dx) exit
        j = j+1
        if(j.gt.nn) goto 21
      enddo

      if(i.le.npl) er(i) = data(j,2)
      pr(i) = data(j,3)

      21 continue
    enddo

    write(9,'(" [---] End Solitary from file")')    
    write(9,*)

    return
    11 write(9,'(" [ERR] Error in solitary input file")')    
    stop

  end subroutine solitICFromFile
!!-----------------------End solitICFromFile-----------------------!!