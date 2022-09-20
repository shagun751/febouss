!!-----------------------------setMFree----------------------------!!
  subroutine setMFree(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::nn,err,i,i2,j,k1,j2
    integer(kind=C_K1),allocatable::neid(:),newrk(:)
    real(kind=C_K2)::cx,cy,rad,tmpr3,tmpr5
    real(kind=C_K2),allocatable::nedr(:),nerad(:)
    real(kind=C_K2),allocatable::phi(:),phiDx(:),phiDy(:)
    real(kind=C_K2)::n3wei(3,3)

    n3wei(1,:) = (/ 0.5d0, 0.5d0, 0.0d0 /)
    n3wei(2,:) = (/ 0.0d0, 0.5d0, 0.5d0 /)
    n3wei(3,:) = (/ 0.5d0, 0.0d0, 0.5d0 /)

    call system_clock(b%sysC(5))
    allocate(b%pObf(b%npt))
    allocate(neid(b%npt), newrk(b%npt), nedr(b%npt))
    allocate(nerad(b%npt))
    allocate(phi(b%npt),phiDx(b%npt),phiDy(b%npt))

    phi=0d0
    phiDx=0d0
    phiDy=0d0

    !tmpr5=1.2d0 !Coeff of multipliciation to max rad in linkList    
    tmpr5=1.8d0 !Coeff of multipliciation to min rad in linkList    

    ! First find influence radius
    do i=1,b%npl !vertex nodes first
      call findRadLinkList(i, b%npl, b%Sz(1), b%ivl, b%linkl, &
        b%cor(1:b%npl,:), tmpr5, rad, i2, j, k1, tmpr3, cx, cy)    
      b%pObf(i)%cx = cx
      b%pObf(i)%cy = cy
      b%pObf(i)%rad = rad
    enddo  
    do k1 = 1, b%nele !For side-midpoint nodes
      do i2 = 4,6
        i = b%conn(k1,i2)
        b%pObf(i)%cx = b%cor(i,1)
        b%pObf(i)%cy = b%cor(i,2)
        rad=0d0
        do j2 = 1,3
          rad = rad + n3wei(i2-3,j2) * b%pObf(b%conn(k1,j2))%rad
        enddo
        b%pObf(i)%rad = rad
      enddo
    enddo


    ! Find neigh and set the mf object
    ! Find MLS shape fnc for interp and 1st derivative
    do i=1,b%npt
      cx = b%pObf(i)%cx
      cy = b%pObf(i)%cy
      rad = b%pObf(i)%rad

      call findNeiLinkList(i, rad, b%npt, b%Sz(4), b%ivq, &
        b%linkq, b%cor, b%npt, nn, neid, newrk, nedr)          

      ! call findNeiBruteForce(i, rad, b%npt, b%Sz(4), b%ivq, &
      !   b%linkq, b%cor, b%npt, nn, neid)
          
      nerad(1:nn) = b%pObf(neid(1:nn))%rad

      call mls2DDx(cx, cy, nn, &
        b%cor(neid(1:nn),1), b%cor(neid(1:nn),2), nerad(1:nn), &
        phi(1:nn), phiDx(1:nn), phiDy(1:nn), err)

      if(err.ne.0)then
        write(9,'(" [ERR] No MFree at node ", I10)')i
        write(9,'(" [---] Cx, Cy ",2F15.6)')cx,cy
      endif      

      call b%pObf(i)%setPoiNZOnly(nn, nn, i, cx, cy, rad, &
        neid(1:nn), phi(1:nn), phiDx(1:nn), phiDy(1:nn), j)      
      

      ! !if(i.ne.2341) cycle !rect2d vertex
      ! if(i.ne.15151) cycle !rect2d side-mid
      ! write(*,'(I15,2F15.6)')b%pObf(i)%bsnqId,b%pObf(i)%cx,b%pObf(i)%cy
      ! write(*,'(2I15,F15.6)')b%pObf(i)%nn,b%pObf(i)%nnMax,b%pObf(i)%rad
      ! nn = b%pObf(i)%nn
      ! do j=1,nn
      !   j2 = b%pObf(i)%neid(j)        
      !   tmpr3 = dsqrt( (b%cor(j2,1) - b%pObf(i)%cx)**2 + &
      !    (b%cor(j2,2) - b%pObf(i)%cy)**2 )
      !   write(*,'(I15,8F15.6)')j2, b%cor(j2,1:2), tmpr3/0.2032d0, &
      !     b%pObf(j2)%rad, tmpr3/b%pObf(j2)%rad, & 
      !     b%pObf(i)%phi(j), b%pObf(i)%phiDx(j), b%pObf(i)%phiDy(j)
      ! enddo
      ! write(*,'(F15.6)')sum(b%pObf(i)%phi(1:nn))
    enddo      

    deallocate(neid,newrk,nedr,phi,phiDx,phiDy,nerad)

    !call b%bDf%init( b%npt )

    ! iDf will be used for time-interpolated values.
    ! will only do velocities.
    ! wont do pressure in iDf
    allocate( b%bDf(b%npt), b%iDf(b%npt))
    do i=1,b%npt
      b%bDf(i)%ug=0d0
      b%bDf(i)%Pg=0d0    
      b%bDf(i)%ug_tn=0d0
      b%bDf(i)%Pg_tn=0d0

      b%iDf(i)%ug=0d0
      b%iDf(i)%Pg=0d0        
    enddo

    ! call testMls2DDx
    ! stop    

    ! Set vvMsh to output vel at -h/3 for all lin nodes
    call b%vvMsh%initvvProbes( b%npl )
    do i =1, b%npl
      b%vvMsh%x(i) = b%cor(i,1)
      b%vvMsh%y(i) = b%cor(i,2)
      b%vvMsh%z(i) = -b%dep(i)/3d0
    enddo

    call system_clock(b%sysC(6))
    write(9,'(" [MSG] Size of mFree : npt, sum(nn), MB ",2I15,F15.6)') &
      b%npt, sum(b%pObf(:)%nn), &
      sum(b%pObf(:)%nn)*(8d0*3+4d0*1)/(1024d0*1024d0)
    write(9,*)"[MSG] Done setMFree"
    write(9,'(" [TIM] ",F15.4)')1d0*(b%sysC(6)-b%sysC(5))/b%sysRate
    write(9,*)    

  end subroutine setMFree
!!---------------------------End setMFree--------------------------!!



!!-----------------------calcDepResDerivAll------------------------!!
  subroutine calcDepResDerivAll(b, npt, p, q, tD, rDf)
  implicit none

    class(bsnqCase),intent(in)::b
    integer(kind=C_K1),intent(in)::npt
    real(kind=C_K2),intent(in)::p(npt), q(npt), tD(npt)
    type(vertVelDerv3D),intent(out)::rDf(npt)

    integer(kind=C_K1)::i,i2,j,k2
    real(kind=C_K2)::tmpx, tmpy    


    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i,i2,j,k2,tmpx,tmpy)
    !$OMP DO SCHEDULE(dynamic,1000)    
    do i=1,npt
      rDf(i)%u = p(i) / tD(i)
      rDf(i)%v = q(i) / tD(i)
      rDf(i)%P = rDf(i)%u * b%dep(i) ! P = uh , note its not u(h+eta)
      rDf(i)%Q = rDf(i)%v * b%dep(i)

      ! Storing ! d(U)/dx and d(Uh)/dx at t(n-1), t(n-2)    
      rDf(i)%ug_tn(2) = rDf(i)%ug_tn(1)  
      rDf(i)%ug_tn(1) = rDf(i)%ug
      rDf(i)%Pg_tn(2) = rDf(i)%Pg_tn(1)    
      rDf(i)%Pg_tn(1) = rDf(i)%Pg
    enddo
    !$OMP END DO NOWAIT    
        
    !b%ur = dsin(b%cor(:,1))

    !! First derivative    
    !$OMP DO SCHEDULE(dynamic,1000)    
    do i=1, npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! dudx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf( b%pObf(i)%neid )%u, tmpx, j )

      ! dvdy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf( b%pObf(i)%neid )%v, tmpy, j )

      ! du
      rDf(i2)%ug = tmpx + tmpy

      ! dPdx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf( b%pObf(i)%neid )%P, tmpx, j )

      ! dQdy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf( b%pObf(i)%neid )%Q, tmpy, j )

      ! dP
      rDf(i2)%Pg = tmpx + tmpy

      ! dhdx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%dep( b%pObf(i)%neid ), rDf(i2)%hx, j )

      ! dhdy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        b%dep( b%pObf(i)%neid ), rDf(i2)%hy, j )      

    enddo
    !$OMP END DO NOWAIT


    !! Second derivative    
    !$OMP DO SCHEDULE(dynamic,1000)    
    do i=1, npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! d( del.(u) ) / dx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf(b%pObf(i)%neid)%ug, rDf(i2)%ugx, j )

      ! d( del.(u) ) / dy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf(b%pObf(i)%neid)%ug, rDf(i2)%ugy, j )

      ! d( del.(P) ) / dx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf(b%pObf(i)%neid)%Pg, rDf(i2)%Pgx, j )

      ! d( del.(P) ) / dy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf(b%pObf(i)%neid)%Pg, rDf(i2)%Pgy, j )

    enddo
    !$OMP END DO NOWAIT


    !! Third derivative and time-derivative
    !$OMP DO SCHEDULE(dynamic,1000)    
    do i=1, npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! ugxx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf(b%pObf(i)%neid)%ugx, tmpx, j )

      !ugyy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf(b%pObf(i)%neid)%ugy, tmpy, j )

      ! del.( del( del.(u) ) )
      rDf(i2)%uggg = tmpx + tmpy      

      ! Pgxx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf(b%pObf(i)%neid)%Pgx, tmpx, j )

      !Pgyy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf(b%pObf(i)%neid)%Pgy, tmpy, j )

      ! del.( del( del.(P) ) )
      rDf(i2)%Pggg = tmpx + tmpy            


      ! Time-derivative for pressure
      rDf(i)%ugt = ( 3d0*rDf(i)%ug - 4d0*rDf(i)%ug_tn(1) &
        + rDf(i)%ug_tn(2) )/ (2d0*b%dt)
      rDf(i)%Pgt = ( 3d0*rDf(i)%Pg - 4d0*rDf(i)%Pg_tn(1) &
        + rDf(i)%Pg_tn(2) )/ (2d0*b%dt)

    enddo
    !$OMP END DO NOWAIT          
    !$OMP END PARALLEL


  end subroutine calcDepResDerivAll
!!---------------------End calcDepResDerivAll----------------------!!



!!-------------------------calcDepResDeriv-------------------------!!
  subroutine calcDepResDeriv(b, npt, p, q, tD, rDf)
  implicit none

    class(bsnqCase),intent(in)::b    
    integer(kind=C_K1),intent(in)::npt
    real(kind=C_K2),intent(in)::p(npt), q(npt), tD(npt)
    type(vertVelDerv3D),intent(out)::rDf(npt)

    integer(kind=C_K1)::i,i2,j,k2
    real(kind=C_K2)::tmpx, tmpy    


    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i,i2,j,k2,tmpx,tmpy)
    do i = 1, npt
      rDf(i)%u = p(i) / tD(i)
      rDf(i)%v = q(i) / tD(i)
      rDf(i)%P = rDf(i)%u * b%dep(i) ! P = uh , note its not u(h+eta)
      rDf(i)%Q = rDf(i)%v * b%dep(i)
    enddo
      
    !! First derivative    
    !$OMP DO SCHEDULE(dynamic,1000)
    do i = 1, npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! dudx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf( b%pObf(i)%neid )%u, tmpx, j )

      ! dvdy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf( b%pObf(i)%neid )%v, tmpy, j )

      ! du
      rDf(i2)%ug = tmpx + tmpy

      ! dPdx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf( b%pObf(i)%neid )%P, tmpx, j )

      ! dQdy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf( b%pObf(i)%neid )%Q, tmpy, j )

      ! dP
      rDf(i2)%Pg = tmpx + tmpy

      ! dhdx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%dep( b%pObf(i)%neid ), rDf(i2)%hx, j )

      ! dhdy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        b%dep( b%pObf(i)%neid ), rDf(i2)%hy, j )      

    enddo
    !$OMP END DO NOWAIT    


    !! Second derivative    
    !$OMP DO SCHEDULE(dynamic,1000)
    do i=1,npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! d( del.(u) ) / dx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf(b%pObf(i)%neid)%ug, rDf(i2)%ugx, j )

      ! d( del.(u) ) / dy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf(b%pObf(i)%neid)%ug, rDf(i2)%ugy, j )

      ! d( del.(P) ) / dx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf(b%pObf(i)%neid)%Pg, rDf(i2)%Pgx, j )

      ! d( del.(P) ) / dy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf(b%pObf(i)%neid)%Pg, rDf(i2)%Pgy, j )

    enddo
    !$OMP END DO NOWAIT    


    !! Third derivative    
    !$OMP DO SCHEDULE(dynamic,1000)    
    do i=1,npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! ugxx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf(b%pObf(i)%neid)%ugx, tmpx, j )

      !ugyy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf(b%pObf(i)%neid)%ugy, tmpy, j )

      ! del.( del( del.(u) ) )
      rDf(i2)%uggg = tmpx + tmpy      

      ! Pgxx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        rDf(b%pObf(i)%neid)%Pgx, tmpx, j )

      !Pgyy
      call calcGrad( k2, b%pObf(i)%phiDy, &
        rDf(b%pObf(i)%neid)%Pgy, tmpy, j )

      ! del.( del( del.(P) ) )
      rDf(i2)%Pggg = tmpx + tmpy            

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL    
    
    ! Pressure related terms not calculated.  


  end subroutine calcDepResDeriv
!!-----------------------End calcDepResDeriv-----------------------!!



!!-------------------------findEleForLocXY-------------------------!!
  subroutine findEleForLocXY1(b,xin,yin,eleOut,ep,et)
  implicit none

    class(bsnqCase),intent(in)::b
    real(kind=C_K2),intent(in)::xin,yin
    integer(kind=C_K1),intent(out)::eleOut
    real(kind=C_K2),intent(out)::ep,et
    integer(kind=C_K1)::ind(3,2),p1,p2,dir,liel,li,nl(3)
    real(kind=C_K2)::xy(3,2),vec1(2),vec2(2),res

    ind(1,:)=(/ 1,2 /)
    ind(2,:)=(/ 2,3 /)
    ind(3,:)=(/ 3,1 /)

    do liel = 1, b%nele
      nl=b%conn(liel,1:3)
      do li=1,3
        xy(li,:)=b%cor(nl(li),:)
      enddo

      dir=1
      do li=1,3
        p1=ind(li,1)
        p2=ind(li,2)
        vec1(1)=xy(p2,1)-xy(p1,1)
        vec1(2)=xy(p2,2)-xy(p1,2)
        vec2(1)=xin-xy(p1,1)
        vec2(2)=yin-xy(p1,2)
        call vecCross2D(vec1(1),vec1(2),vec2(1),vec2(2),res)
        if(res.lt.0d0) dir=-1        
      enddo

      if(dir.eq.1)then
        eleOut=liel
        ep = ( b%invJ(liel,1)*(xin-xy(1,1)) ) &
          + ( b%invJ(liel,3)*(yin-xy(1,2)) )
        et = ( b%invJ(liel,2)*(xin-xy(1,1)) ) &
          + ( b%invJ(liel,4)*(yin-xy(1,2)) )
        return
      endif

    enddo

    eleOut=-1
    ep=-1
    et=-1

  end subroutine findEleForLocXY1


  subroutine findEleForLocXY2(b,np,xin,yin,eleOut,ep,et)
  implicit none
    
    class(bsnqCase),intent(in)::b
    integer(kind=C_K1),intent(in)::np
    real(kind=C_K2),intent(in)::xin(np),yin(np)
    integer(kind=C_K1),intent(out)::eleOut(np)
    real(kind=C_K2),intent(out)::ep(np),et(np)
    integer(kind=C_K1)::ind(3,2),p1,p2,dir,pp,liel,li,nl(3)
    real(kind=C_K2)::xy(3,2),vec1(2),vec2(2),res,xp,yp

    ind(1,:)=(/ 1,2 /)
    ind(2,:)=(/ 2,3 /)
    ind(3,:)=(/ 3,1 /)


    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(pp,xp,yp,liel,nl,xy,dir,li,p1,p2,vec1,vec2,res)
    !$OMP DO SCHEDULE(dynamic,100)
    do pp=1,np
      xp=xin(pp)
      yp=yin(pp)

      do liel = 1, b%nele
        nl=b%conn(liel,1:3)
        do li=1,3
          xy(li,:)=b%cor(nl(li),:)
        enddo

        dir=1
        do li=1,3
          p1=ind(li,1)
          p2=ind(li,2)
          vec1(1)=xy(p2,1)-xy(p1,1)
          vec1(2)=xy(p2,2)-xy(p1,2)
          vec2(1)=xp-xy(p1,1)
          vec2(2)=yp-xy(p1,2)
          call vecCross2D(vec1(1),vec1(2),vec2(1),vec2(2),res)
          if(res.lt.0d0) dir=-1        
        enddo

        if(dir.eq.1)then
          eleOut(pp)=liel
          ep(pp) = ( b%invJ(liel,1)*(xp-xy(1,1)) ) &
            + ( b%invJ(liel,3)*(yp-xy(1,2)) )
          et(pp) = ( b%invJ(liel,2)*(xp-xy(1,1)) ) &
            + ( b%invJ(liel,4)*(yp-xy(1,2)) )
          exit
        endif
      enddo
      
      if(dir.eq.1) cycle

      eleOut(pp)=-1
      ep(pp)=-1
      et(pp)=-1
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end subroutine findEleForLocXY2


  subroutine findEleForLocXY3(b,np,xin,yin,eleOut,ep,et)
  implicit none
    
    class(bsnqCase),intent(in)::b
    integer(kind=C_K1),intent(in)::np
    real(kind=C_K2),intent(in)::xin(np),yin(np)
    integer(kind=C_K1),intent(out)::eleOut(np)
    real(kind=C_K2),intent(out)::ep(np),et(np)
    integer(kind=C_K1)::ind(3,2),p1,p2,dir,pp,liel,li,nl(3)
    integer(kind=C_K1)::celli, cellx, celly, i1
    real(kind=C_K2)::xy(3,2),vec1(2),vec2(2),res,xp,yp

    ind(1,:)=(/ 1,2 /)
    ind(2,:)=(/ 2,3 /)
    ind(3,:)=(/ 3,1 /)


    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(pp,xp,yp,liel,nl,xy,dir,li,&
    !$OMP     p1,p2,vec1,vec2,res, celli, cellx, celly, i1)
    !$OMP DO SCHEDULE(dynamic,100)
    do pp=1,np
      xp=xin(pp)
      yp=yin(pp)

      call b%cofe%getCellNo( xp, yp, cellx, celly, celli)

      dir=-1
      do i1 = 1, b%cofe%cell(celli)%nm
        liel = b%cofe%cell(celli)%m(i1)
        nl=b%conn(liel,1:3)
        do li=1,3
          xy(li,:)=b%cor(nl(li),:)
        enddo

        dir=1
        do li=1,3
          p1=ind(li,1)
          p2=ind(li,2)
          vec1(1)=xy(p2,1)-xy(p1,1)
          vec1(2)=xy(p2,2)-xy(p1,2)
          vec2(1)=xp-xy(p1,1)
          vec2(2)=yp-xy(p1,2)
          call vecCross2D(vec1(1),vec1(2),vec2(1),vec2(2),res)
          if(res.lt.0d0) dir=-1        
        enddo

        if(dir.eq.1)then
          eleOut(pp)=liel
          ep(pp) = ( b%invJ(liel,1)*(xp-xy(1,1)) ) &
            + ( b%invJ(liel,3)*(yp-xy(1,2)) )
          et(pp) = ( b%invJ(liel,2)*(xp-xy(1,1)) ) &
            + ( b%invJ(liel,4)*(yp-xy(1,2)) )
          exit
        endif
      enddo
      
      if(dir.eq.1) cycle

      eleOut(pp)=-1
      ep(pp)=-1
      et(pp)=-1
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end subroutine findEleForLocXY3
!!-----------------------End findEleForLocXY-----------------------!!



!!----------------------------getVertVel---------------------------!!
  subroutine getVertVel(b,np,xin,yin,zin,uOut,vOut,wOut,pOut,&
    etaOut,wrki,wrkr,err,rTime)
  implicit none

    class(bsnqCase),intent(in)::b
    integer(kind=C_K1),intent(in)::np
    real(kind=C_K2),intent(in)::xin(np),yin(np),zin(np)

    integer(kind=C_K1),intent(out)::wrki(np),err(np)
    real(kind=C_K2),intent(out)::uOut(np),vOut(np),wOut(np)
    real(kind=C_K2),intent(out)::pOut(np),etaOut(np),wrkr(np,2)    
    real(kind=C_K2),intent(out)::rTime

    integer(kind=C_KCLK)::lsysC(2)
    integer(kind=C_K1)::nq(6),i,k
    real(kind=C_K2)::wei(6),hLoc,zLoc,etaLoc
    real(kind=C_K2)::uLoc,vLoc,wLoc,pLoc    
    type(vertVelDerv3D)::tmp

    !Note: pLoc here is pressure. It is not depth-integrated vel-x

    call system_clock(lsysC(1)) 

    call b%findEleForLocXY3(np,xin,yin,wrki,wrkr(:,1),wrkr(:,2))

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i, nq, wei, tmp, k, hLoc, etaLoc, zLoc, &
    !$OMP     uLoc, vLoc, wLoc, pLoc)
    !$OMP DO SCHEDULE(dynamic,100)
    do i = 1,np

      if(wrki(i).eq.-1)then
        err(i)=1
        uOut(i)=0d0
        vOut(i)=0d0
        wOut(i)=0d0
        pOut(i)=0d0     
        etaOut(i)=0d0
        cycle   
      endif

      nq = b%conn(wrki(i),:)
      call fem_N6i(wrkr(i,1),wrkr(i,2),wei)      

      call tmp%initByInterp3D( 6, b%bDf(nq), wei, k )

      hLoc=0d0
      etaLoc=0d0      
      do k=1,6
        hLoc = hLoc + wei(k)*b%dep(nq(k))
        etaLoc = etaLoc + wei(k)*b%tOb(0)%tD(nq(k))
      enddo
      etaLoc = etaLoc - hLoc
      zLoc = zin(i) !zRef is mean sea level

      ! call vertVelExp(zLoc, hLoc, etaLoc, tmp%hx, tmp%hy, &
      !   tmp%u, tmp%ux, tmp%uxx, tmp%uxxx, &
      !   tmp%uhx, tmp%uhxx, tmp%uhxxx, &
      !   tmp%uxt, tmp%uhxt, uLoc, wLoc, pLoc)  

      call vertVelExp3D(zLoc, hLoc, etaLoc, tmp%hx, tmp%hy, &
        tmp%u, tmp%v, tmp%P, tmp%Q, tmp%ug, tmp%Pg, &
        tmp%ugx, tmp%ugy, tmp%Pgx, tmp%Pgy, &
        tmp%uggg, tmp%Pggg, tmp%ugt, tmp%Pgt, &
        uLoc, vLoc, wLoc, pLoc)      

      err(i)=0
      uOut(i)=uLoc
      vOut(i)=vLoc
      wOut(i)=wLoc
      pOut(i)=pLoc
      etaOut(i)=etaLoc

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL    

    call system_clock(lsysC(2)) 
    rTime = 1d0 * (lsysC(2) - lsysC(1)) / b%sysRate
    
  end subroutine getVertVel
!!--------------------------End getVertVel-------------------------!!



!!--------------------------testGetVertVel-------------------------!!
  subroutine testGetVertVel(b)
  implicit none

    class(bsnqCase),intent(in)::b
    
    integer,parameter::np=4
    integer(kind=C_K1)::wrki(np),err(np),i
    real(kind=C_K2)::xin(np),yin(np),zin(np),wrkr(np,2)
    real(kind=C_K2)::uOut(np),vOut(np),wOut(np),pOut(np)    
    real(kind=C_K2)::etaOut(np), lRtime
    
    xin(1)=15d0
    yin(1)=2.00d0
    zin(1)=-0.20d0

    xin(2)=15d0
    yin(2)=2.00d0
    zin(2)=-0.35d0

    xin(3)=16.02d0
    yin(3)=2.02d0
    zin(3)=-0.20d0

    xin(4)=16.02d0
    yin(4)=2.02d0
    zin(4)=-0.35d0

    call b%getVertVel(np,xin,yin,zin,uOut,vOut,wOut,pOut,&
      etaOut,wrki,wrkr,err,lRtime)

    write(120,'(F15.6)',advance='no')b%tOb(0)%rtm
    do i=1,np
      write(120,'(5F15.6)',advance='no')pOut(i),uOut(i),&
        vOut(i),wOut(i),etaOut(i)
    enddo
    write(120,*)

  end subroutine testGetVertVel
!!------------------------End testGetVertVel-----------------------!!



!!----------------------------locWvAng-----------------------------!!
  ! subroutine locWvAng(b)
  ! implicit none

  !   class(bsnqCase),intent(inout)::b

  !   integer(kind=C_K1)::i
  !   real(kind=C_K2)::p,q,pMag2
  !   real(kind=C_K2),parameter::velLowLimit2=1d-20

  !   !$OMP PARALLEL DEFAULT(shared) PRIVATE(i, p, q, pMag2)
  !   !$OMP DO SCHEDULE(dynamic,100)
  !   do i = 1, b%npt
  !     p = b%tOb(0)%p(i)
  !     q = b%tOb(0)%q(i)
  !     pMag2 = p**2 + q**2
  !     if(pMag2.lt.velLowLimit2) then 
  !       b%wvAng(i) = 0d0
  !     else  
  !       ! RESULT = ATAN2(Y, X) -pi to pi
  !       b%wvAng(i) = atan2(q, p)
  !     endif
  !   enddo
  !   !$OMP END DO NOWAIT
  !   !$OMP END PARALLEL    

  ! end subroutine locWvAng

  
  subroutine locWvAng(b, npt, eta6)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1),intent(in)::npt
    real(kind=C_K2),intent(in)::eta6(npt)

    integer(kind=C_K1)::i, i2, nn, j
    real(kind=C_K2)::etx,ety

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i, i2, nn, j, etx, ety)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i = 1, b%npt
      i2 = b%pObf(i)%bsnqId
      nn = b%pObf(i)%nn

      call calcGrad( nn, b%pObf(i)%phiDx, &
        eta6( b%pObf(i)%neid ), etx, j )      

      call calcGrad( nn, b%pObf(i)%phiDy, &
        eta6( b%pObf(i)%neid ), ety, j )      
      
      ! RESULT = ATAN2(Y, X) -pi to pi
      b%wvAng(i2) = atan2(ety, etx)
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

  end subroutine locWvAng
!!--------------------------End locWvAng---------------------------!!



!!-------------------------setDepResAtTime-------------------------!!
  ! subroutine setDepResAtTime(b, tin, err)
  !   implicit none

  !   class(bsnqCase),intent(inout)::b

  !   real(kind=C_K2),intent(in)::tin
  !   integer(kind=C_K1),intent(out)::err

  !   integer(kind=C_K1)::i
  !   real(kind=C_K2)::tr, dtin

  !   ! b%tOb(1) = t(n-1)
  !   ! b%tOb(0) = t(n)

  !   ! err = 0 no error
  !   ! err = 1 some error

  !   err = 0

  !   dtin = (tin-b%tOb(1)%rtm)
  !   tr =  dtin / b%dt

  !   if((tr.lt.0d0) .or. (tr.gt.1d0))then
  !     write(9,'(" [ERR] setDepResAtTime")')
  !     write(9,'(" [---] Tin not within t(n-1) and t(n)")')
  !     write(9,'(" [---] ",3a15)')'Tin', 't(n-1)', 't(n)'
  !     write(9,'(" [---] ",3f15.6)')tin, b%tOb(1)%rtm, b%tOb(0)%rtm
  !     err = 1
  !     return
  !   endif

  !   call b%rk4intp%setrk4InterpMatrix( b%dt, dtin)

  !   ! Cant use for now
  !   ! Made probably a humungous blunder in RK4 time-stepping
  !   !$OMP PARALLEL DEFAULT(shared) PRIVATE(i)
  !   !$OMP DO SCHEDULE(dynamic,1000)
  !   do i = 1, b%npl
  !     call b%rk4intp%get_ktilde( b%sOb(1:4)%e(i), &
  !       b%siOb(1:4)%e(i) )
  !   enddo
  !   !$OMP END DO NOWAIT

  !   !$OMP DO SCHEDULE(dynamic,1000)
  !   do i = 1, b%npt
  !     call b%rk4intp%get_ktilde( b%sOb(1:4)%p(i), &
  !       b%siOb(1:4)%p(i) )

  !     call b%rk4intp%get_ktilde( b%sOb(1:4)%q(i), &
  !       b%siOb(1:4)%q(i) )
  !   enddo
  !   !$OMP END DO NOWAIT
  !   !$OMP END PARALLEL          


  !   b%tiOb%rtm = tin
  !   b%tiOb%e = b%tOb(1)%e + 1d0/6d0*(b%siOb(1)%e &
  !     + 2d0*b%siOb(2)%e + 2d0*b%siOb(3)%e + b%siOb(4)%e)
  !   b%tiOb%p = b%tOb(1)%p + 1d0/6d0*(b%siOb(1)%p &
  !     + 2d0*b%siOb(2)%p + 2d0*b%siOb(3)%p + b%siOb(4)%p)
  !   b%tiOb%q = b%tOb(1)%q + 1d0/6d0*(b%siOb(1)%q &
  !     + 2d0*b%siOb(2)%q + 2d0*b%siOb(3)%q + b%siOb(4)%q)

  !   !! Forcing Dirichlet BC
  !   call b%diriBCEta(b%tiOb%e,b%tiOb%rtm)
  !   call b%diriBCPQ(b%tiOb%p,b%tiOb%q,b%tiOb%rtm)
    
  !   b%tiOb%tD(1:b%npl) = b%dep(1:b%npl) + b%tiOb%e
  !   call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tiOb%tD)

  !   call b%calcDepResDeriv(b%npt, b%tiOb%p, b%tiOb%q, &
  !     b%tiOb%tD, b%iDf)

  ! end subroutine setDepResAtTime
!!-----------------------End setDepResAtTime-----------------------!!



