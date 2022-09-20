program genTriMesh
implicit none

  integer(kind=4)::i,j,k,l,i2,j2,k2,nx,ny,en(3)
  integer(kind=4)::npoi,nele,nbnd,nedg,bndt1,bndt2
  integer(kind=4),allocatable::conn(:,:),mabnd(:,:)
  real(kind=8)::domx(2),domy(2)
  real(kind=8)::tmpr1,tmpr2,tmpr3,tmpr4,tmpr5,dx,dy
  real(kind=8),allocatable::corx(:),cory(:)
  real(kind=8),allocatable::dep(:)
  character(len=100)::probname

  probname='test'

  !Bottom left X Y
  domx(1)=0d0
  domy(1)=0d0

  !Top right X Y
  domx(2)=150.4d0
  domy(2)=4.8d0

  !Mesh dx and dy
  dx=0.16d0
  dy=0.16d0

  !!-------------- Do not change --------------!!
  open(11,file = "test.msh")

  nx=int((domx(2)-domx(1))/dx,4)
  ny=int((domy(2)-domy(1))/dy,4)
  write(*,*)nx,ny

  npoi=(nx+1)*(ny+1)
  nele=nx*ny*2
  nedg=3*nx*ny+nx+ny
  nbnd=2*(nx+ny)

  allocate(corx(npoi),cory(npoi),conn(nele,3))
  allocate(mabnd(nbnd,3),dep(npoi))

  !! Nodes
  k=0
  do i=1,nx+1
    do j=1,ny+1
      k=k+1
      corx(k)=(i-1)*dx+domx(1)
      cory(k)=(j-1)*dy+domy(1)
    enddo
  enddo
  write(*,*)k,npoi
  ! do i=1,npoi
  !   write(*,'(2F15.6)')corx(i),cory(i)
  ! enddo

  !! Elements
  k=0 
  do i=1,nx
    do j=1,ny

      if((j*1d0 - ny/2d0).le.0)then
        !   /|
        !  / |
        ! /__|
        k=k+1
        k2=(i-1)*(ny+1)+j
        en(1)=k2
        en(2)=k2+ny+1
        en(3)=k2+ny+2      
        conn(k,:)=en

        ! |  /
        ! | /
        ! |/
        k=k+1
        k2=(i-1)*(ny+1)+j
        en(1)=k2
        en(2)=k2+ny+2
        en(3)=k2+1      
        conn(k,:)=en

      else
        ! |\
        ! | \
        ! |__\
        k=k+1
        k2=(i-1)*(ny+1)+j
        en(1)=k2
        en(2)=k2+ny+1
        en(3)=k2+1
        conn(k,:)=en

        ! \  |
        !  \ |
        !   \|
        k=k+1
        k2=(i-1)*(ny+1)+j
        en(1)=k2+ny+1
        en(2)=k2+ny+2
        en(3)=k2+1      
        conn(k,:)=en
      endif

    enddo
  enddo
  write(*,*)k,nele
  ! do i=1,nele
  !   write(*,'(4I15)')conn(i,:)
  ! enddo

  !! Boundaries
  !Left
  k=0
  do i=1,ny
    k=k+1
    mabnd(k,1)=i+1
    mabnd(k,2)=i
    if((i*1d0 - ny/2d0).le.0)then
      mabnd(k,3)=i*2
    else
      mabnd(k,3)=i*2-1
    endif
  enddo
  bndt1=k

  !Bottom
  do i=1,nx
    k=k+1
    j = (i-1)*(ny+1)+1
    mabnd(k,1)=j
    mabnd(k,2)=j+ny+1
    mabnd(k,3)=(i-1)*ny*2+1    
  enddo
  !Top
  do i=nx,1,-1
    k=k+1
    j = i*(ny+1)
    mabnd(k,1)=j+ny+1
    mabnd(k,2)=j
    mabnd(k,3)=i*ny*2
  enddo
  bndt2 = k

  !Right
  do i=1,ny
    k=k+1
    j = nx*(ny+1) + i
    mabnd(k,1)=j
    mabnd(k,2)=j+1
    if((i*1d0 - ny/2d0).le.0)then
      mabnd(k,3)= (nx-1)*ny*2 + i*2 - 1
    else
      mabnd(k,3)= (nx-1)*ny*2 + i*2 
    endif
  enddo  
  write(*,*)k,nbnd
  !!------------ End Do not change ------------!!

  

  !!---------- Can change depth code-----------!!

  ! Constant depth
  dep=1

  ! ! Depth code  with slope
  ! tmpr5 = 1d0/20d0  !! Slope
  ! tmpr3 = 300d0
  ! tmpr4 = 376d0
  ! do i=1,npoi
  !   if( corx(i) .le. tmpr3)then
  !     dep(i) = 4
    
  !   elseif( ( corx(i) .gt. tmpr3) .and. ( corx(i) .le. tmpr4) )then
  !     dep(i) = 4 - (corx(i) - tmpr3)*tmpr5
    
  !   else
  !     dep(i) = 0.2
  !   endif    
  ! enddo

  ! Whalin
  
  !!---------End Can change depth code---------!!
  
  

  !!-------------- Do not change --------------!!
  call outputXML(npoi,nele,conn,corx,cory,dep)
  ! call out4NXML(probname,npoi,npoi,nele,21,0,&
  ! conn,corx,cory,dep,dep,dep,dep)

  !! Mesh File
  write(11,'(3A15)')"Elements","Nodes", "Edges"
  write(11,'(3I15)')nele,npoi,nedg
  write(11,'(2A15)')"BndSides","BndTypes"
  write(11,'(2I15)')nbnd,3
  write(11,'("Nodes")')
  do i=1,npoi
    write(11,'(2F20.8)')corx(i),cory(i)
  enddo
  write(11,'("Elements")')
  do i=1,nele
    write(11,'(3I15)')conn(i,:)
  enddo
  write(11,'("Boundary")')
  write(11,'(2I15)')11,bndt1
  do i=1,bndt1
    write(11,'(3I15)')mabnd(i,:)
  enddo
  write(11,'("Boundary")')
  write(11,'(2I15)')13,bndt2-bndt1
  do i=bndt1+1,bndt2
    write(11,'(3I15)')mabnd(i,:)
  enddo
  write(11,'("Boundary")')
  write(11,'(2I15)')12,nbnd - bndt2
  do i=bndt2+1,nbnd
    write(11,'(3I15)')mabnd(i,:)
  enddo  
  write(11,'("Depth")')
  do i=1,npoi
    write(11,'(F15.6)')dep(i)
  enddo
  write(11,'("####")')
  close(11)
  !!------------ End Do not change ------------!!


end program genTriMesh



subroutine outputXML(npoinl,nelem,conn,corx,cory,depth)
implicit none
  
  integer(kind=4),intent(in)::npoinl,nelem
  integer(kind=4),intent(in)::conn(nelem,3)  
  integer(kind=4)::i,j,k,code

  real(kind=8),intent(in)::corx(npoinl), cory(npoinl)
  real(kind=8),intent(in)::depth(npoinl)  
  
  character(len=100)::text

  code=12

  text="test.vtu"
  open(code,file=text(1:len_trim(text)))
  write(code,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
  write(code,'(T3,a)')'<UnstructuredGrid>'
  write(code,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',npoinl,'" NumberOfCells="',nelem,'">'

  write(code,'(T5,a)')'<PointData Scalars="dep">'
  
  write(code,'(T7,a)')'<DataArray type="Float64" Name="dep" format="ascii">'
  write(code,*)-depth
  write(code,'(T7,a)')'</DataArray>'  
  
  write(code,'(T5,a)')'</PointData>'  

  write(code,'(T5,a)')'<CellData>'
  write(code,'(T5,a)')'</CellData>'

  write(code,'(T5,a)')'<Points>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,*)corx(i),cory(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Points>'

  write(code,'(T5,a)')'<Cells>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="connectivity" format="ascii">'  
  do i=1,nelem
    write(code,*)conn(i,1:3)-1
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="offsets" format="ascii">'  
  do i=1,nelem
    write(code,*)3*i
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
  do i=1,nelem
    write(code,*)5
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Cells>'

  write(code,'(T5,a)')'</Piece>'
  write(code,'(T3,a)')'</UnstructuredGrid>'
  write(code,'(a)')'</VTKFile>'
  close(code)

end subroutine outputXML
