!!----------------------------outputXML----------------------------!!
  subroutine outputXML(b)
  implicit none

    class(bsnqCase),intent(inout)::b 

    integer(kind=C_K1)::mf,i,j,k, i2, k2
    real(kind=C_K2)::tmpr1,tmpr2,tmpr3
    character(len=C_KSTR)::bqtxt

    write(bqtxt,'(I15)')int(b%tOb(0)%rtm*1000)
    bqtxt=adjustl(bqtxt)
    bqtxt="Output/"//trim(b%probname)//"_"//trim(bqtxt)//".vtu"
    open(newunit=mf,file=trim(bqtxt))

    write(mf,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
    write(mf,'(T3,a)')'<UnstructuredGrid>'
    write(mf,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',b%npl,'" NumberOfCells="',b%nele,'">'

    ! PointData
    write(mf,'(T5,a)')'<PointData Scalars="eta" Vectors="vel">'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="eta" format="ascii">'
    write(mf,'(F20.6)')b%tOb(0)%e
    write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="waveH" format="ascii">'
    write(mf,'(F20.6)')b%etaMax-b%etaMin
    write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="press" format="ascii">'
    write(mf,'(F15.6)')b%presr(1:b%npl)
    write(mf,'(T7,a)')'</DataArray>'    
    
    write(mf,'(T7,a)')'<DataArray type="Float64" Name="depth" format="ascii">'
    write(mf,'(F15.6)')-b%dep(1:b%npl)
    write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="mlsRad" format="ascii">'
    write(mf,'(F15.6)')b%pObf(1:b%npl)%rad
    write(mf,'(T7,a)')'</DataArray>'

    ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="wvAngDeg" format="ascii">'
    ! write(mf,'(F15.6)')b%wvAng(1:b%npl)*rad2deg
    ! write(mf,'(T7,a)')'</DataArray>'

    ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="porH" format="ascii">'
    ! write(mf,*)porH-depth(1:npoinl)
    ! write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="vel" NumberOfComponents="3" format="ascii">'  
    do i=1,b%npl
      write(mf,'(2F20.6,F5.1)')b%tOb(0)%p(i), b%tOb(0)%q(i), 0d0
      !write(mf,'(2F20.6,F5.1)')b%bndPN(i,1), b%bndPN(i,2), 0d0
    enddo
    write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="velDepRes" NumberOfComponents="3" format="ascii">'  
    do i=1,b%npl
      write(mf,'(3F20.6)')b%vvMsh%u(i), b%vvMsh%v(i), b%vvMsh%w(i)
    enddo
    write(mf,'(T7,a)')'</DataArray>'

    !sin(x)
    write(mf,'(T7,a)')'<DataArray type="Float64" Name="sinx" NumberOfComponents="3" format="ascii">'  
    b%bDf%u = dsin(b%cor(:,1)*2d0*pi/(2.032*3d0))
    do i=1,b%npt      
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! dudx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf( b%pObf(i)%neid )%u, b%bDf(i2)%ug, j )
      
    enddo

    do i=1,b%npt      
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! dudx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf( b%pObf(i)%neid )%ug, b%bDf(i2)%ugx, j )
      
    enddo

    do i=1,b%npt      
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      ! dudx
      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf( b%pObf(i)%neid )%ugx, b%bDf(i2)%uggg, j )
      
    enddo

    do i=1,b%npl      
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      write(mf,'(3F20.6)')b%bDf(i2)%ug, b%bDf(i2)%ugx, b%bDf(i2)%uggg
      
    enddo    
    write(mf,'(T7,a)')'</DataArray>'

    ! gradEta
    ! if(allocated(b%pObf))then
    !   write(mf,'(T7,a)')'<DataArray type="Float64" Name="gradEta" NumberOfComponents="3" format="ascii">'  
    !   do i=1,b%npl
    !     tmpr1=0d0
    !     tmpr2=0d0
    !     do j=1,b%pObf(i)%nn
    !       k = b%pObf(i)%neid(j)
    !       tmpr3 = b%tOb(0)%tD( k ) - b%dep( k )
    !       tmpr1=tmpr1 + b%pObf(i)%phiDx(j) * tmpr3
    !       tmpr2=tmpr2 + b%pObf(i)%phiDy(j) * tmpr3
    !     enddo
    !     if(abs(tmpr1).gt.100d0)tmpr1=-10
    !     if(abs(tmpr2).gt.100d0)tmpr2=-10
    !     write(mf,'(2F15.6,F5.1)')tmpr1,tmpr2,0d0
    !   enddo
    !   write(mf,'(T7,a)')'</DataArray>'

    !   ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="uDerv" NumberOfComponents="3" format="ascii">'  
    !   ! do i=1,b%npl
    !   !   write(mf,'(3F15.6)')b%bDf(i)%ux,b%bDf(i)%uxx,b%bDf(i)%uxxx
    !   ! enddo
    !   ! write(mf,'(T7,a)')'</DataArray>'

    !   ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="pDerv" NumberOfComponents="3" format="ascii">'  
    !   ! do i=1,b%npl
    !   !   write(mf,'(3F15.6)')b%bDf(i)%px,b%bDf(i)%pxx,b%bDf(i)%pxxx
    !   ! enddo
    !   ! write(mf,'(T7,a)')'</DataArray>'
    ! endif

    write(mf,'(T5,a)')'</PointData>'


    ! CellData
    write(mf,'(T5,a)')'<CellData>'
    write(mf,'(T5,a)')'</CellData>'


    ! Mesh
    write(mf,'(T5,a)')'<Points>'
    write(mf,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
    do i=1,b%npl
      write(mf,'(2F15.4,F5.1)')b%cor(i,:),0d0
    enddo
    write(mf,'(T7,a)')'</DataArray>'
    write(mf,'(T5,a)')'</Points>'

    write(mf,'(T5,a)')'<Cells>'
    write(mf,'(T7,a)')'<DataArray type="Int32" Name="connectivity" format="ascii">'  
    do i=1,b%nele
      write(mf,'(3I12)')b%conn(i,1:3)-1
    enddo
    write(mf,'(T7,a)')'</DataArray>'
    write(mf,'(T7,a)')'<DataArray type="Int32" Name="offsets" format="ascii">'  
    do i=1,b%nele
      write(mf,'(I12)')3*i
    enddo
    write(mf,'(T7,a)')'</DataArray>'
    write(mf,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
    do i=1,b%nele
      write(mf,'(I12)')5
    enddo
    write(mf,'(T7,a)')'</DataArray>'
    write(mf,'(T5,a)')'</Cells>'

    write(mf,'(T5,a)')'</Piece>'
    write(mf,'(T3,a)')'</UnstructuredGrid>'
    write(mf,'(a)')'</VTKFile>'
    close(mf)

    write(9,'(" [MSG] OutXML at ",F15.4)')b%tOb(0)%rtm

  end subroutine outputXML
!!--------------------------End outputXML--------------------------!!



!!---------------------------writeResume---------------------------!!
  subroutine writeResume(b)
  implicit none

    ! [Note]:
    ! Check bsnqM/Dev - Logs/log_bsnqM_v0003.md
    ! to understand why the variables 
    ! gXW, gXE, gXPQ are being written in 
    ! resume file.

    class(bsnqCase),intent(in)::b 
    integer(kind=C_K1)::mf,i
    character(len=C_KSTR)::bqtxt

    write(bqtxt,'(I15)')int(b%tOb(0)%rtm*1000)
    bqtxt=adjustl(bqtxt)
    bqtxt="Output/Resume_"//trim(b%probname)//"_"//trim(bqtxt)//".dat"
    open(newunit=mf, file=trim(bqtxt), form='unformatted')

    write(mf) b%tOb(0)%rtm

    do i = 1, b%npl
      write(mf) b%tOb(0)%e(i)
    enddo

    do i = 1, b%npt
      write(mf) b%tOb(0)%p(i)
    enddo

    do i = 1, b%npt
      write(mf) b%tOb(0)%q(i)
    enddo    

    do i = 1, b%npl
      write(mf) b%gXW(i)
    enddo

    do i = 1, b%npl
      write(mf) b%gXE(i)
    enddo

    do i = 1, 2*b%npt
      write(mf) b%gXPQ(i)
    enddo

    close(mf)

  end subroutine writeResume
!!-------------------------End writeResume-------------------------!!



!!---------------------------readResume----------------------------!!
  subroutine readResume(b)
  implicit none

    class(bsnqCase),intent(inout)::b 
    integer(kind=C_K1)::mf,i
    logical(kind=C_LG)::ex

    inquire(file=trim(b%resumeFile),exist=ex)
    if(ex) then
      open(newunit=mf, file=trim(b%resumeFile), &
        form='unformatted')
    else
      write(9,*)"[ERR] Missing resume file"
      write(9,*)"[---]",trim(b%resumeFile)
      stop
    endif      

    read(mf) b%tOb(0)%rtm
    
    do i = 1, b%npl
      read(mf) b%tOb(0)%e(i)
    enddo

    do i = 1, b%npt
      read(mf) b%tOb(0)%p(i)
    enddo

    do i = 1, b%npt
      read(mf) b%tOb(0)%q(i)
    enddo    

    do i = 1, b%npl
      read(mf) b%gXW(i)
    enddo

    do i = 1, b%npl
      read(mf) b%gXE(i)
    enddo

    do i = 1, 2*b%npt
      read(mf) b%gXPQ(i)
    enddo

    close(mf)    

  end subroutine readResume
!!-------------------------End readResume--------------------------!!  
