subroutine matrixSet2(npoinl,npoint,nelem,conn,Sz,ivl,ivq,linkl,&
  linkq,invJ,ele6x6,ele6x3,depth,por,tDr,ur,vr,gGx,gGy,gNAdv,gPGx,gPGy)
use bsnqGlobVars
implicit none

  integer(kind=C_K1),intent(in)::npoinl,npoint,nelem
  integer(kind=C_K1),intent(in)::Sz(4),conn(nelem,6)
  integer(kind=C_K1),intent(in)::ivl(0:npoint),linkl(Sz(3))
  integer(kind=C_K1),intent(in)::ivq(0:npoint),linkq(Sz(4))
  integer(kind=C_K1),intent(in)::ele6x6(nelem,36)
  integer(kind=C_K1),intent(in)::ele6x3(nelem,18)
  integer(kind=C_K1)::i,j,k,i2,j2,k2,n(6),gRow,gCol,lRow,lCol
  integer(kind=C_K1)::nlinkl(ivl(0)),nlinkq(ivq(0))

  real(kind=C_K2),intent(in)::invJ(nelem,5),depth(npoint)  
  real(kind=C_K2),intent(in)::por(npoint),tDr(npoint)
  real(kind=C_K2),intent(in)::ur(npoint),vr(npoint)
  real(kind=C_K2),intent(out)::gGx(Sz(3)),gPGx(Sz(4))
  real(kind=C_K2),intent(out)::gGy(Sz(3)),gPGy(Sz(4))
  real(kind=C_K2),intent(out)::gNAdv(Sz(4))
  real(kind=C_K2)::lGx(6,3),lGy(6,3),lNAd(6,6),tmp(6,6)
  real(kind=C_K2)::lPGx(6,6),lPGy(6,6)
  real(kind=C_K2)::lScN3(3),lScN6(6)
  

  gGx=0d0
  gGy=0d0  
  gNAdv=0d0
  gPGx=0d0
  gPGy=0d0

  !$OMP PARALLEL DEFAULT(shared) &
  !$OMP   PRIVATE(i, n, lScN3, lScN6, lGx, lGy, tmp, lNAd, &
  !$OMP     lPGx, lPGy, lRow, lCol, j2)
  !$OMP DO SCHEDULE(dynamic,100)  
  do i=1,nelem
    n=conn(i,:)    
    
    lScN3=grav*por(n(1:3))*tDr(n(1:3))
    call fem_N6i_Sc3_dN3jdx(lGx,lScN3(1),lScN3(2),lScN3(3),&
      invJ(i,1),invJ(i,2))
    call fem_N6i_Sc3_dN3jdx(lGy,lScN3(1),lScN3(2),lScN3(3),&
      invJ(i,3),invJ(i,4))

    call fem_N6i_du6N6jdx(tmp,ur(n(1)),ur(n(2)),ur(n(3)),&
      ur(n(4)),ur(n(5)),ur(n(6)),invJ(i,1),invJ(i,2))
    lNAd=tmp
    call fem_N6i_du6N6jdx(tmp,vr(n(1)),vr(n(2)),vr(n(3)),&
      vr(n(4)),vr(n(5)),vr(n(6)),invJ(i,3),invJ(i,4))
    lNAd=lNAd+tmp

    lScN6=grav*por(n(1:6))*tDr(n(1:6))
    call fem_N6i_Sc6_dN6jdx(lPGx,lScN6(1),lScN6(2),lScN6(3),&
      lScN6(4),lScN6(5),lScN6(6),invJ(i,1),invJ(i,2))
    call fem_N6i_Sc6_dN6jdx(lPGy,lScN6(1),lScN6(2),lScN6(3),&
      lScN6(4),lScN6(5),lScN6(6),invJ(i,3),invJ(i,4))


    lGx=-invJ(i,5)*lGx
    lGy=-invJ(i,5)*lGy
    lNAd=-invJ(i,5)*lNAd
    lPGx=-invJ(i,5)*lPGx
    lPGy=-invJ(i,5)*lPGy
    
    
    !$OMP CRITICAL
    do lRow=1,6
      !6x6
      do lCol=1,6        
        j2=ele6x6(i,(lRow-1)*6+lCol)
        gNAdv(j2)=gNAdv(j2)+lNAd(lRow,lCol)
        gPGx(j2)=gPGx(j2)+lPGx(lRow,lCol)
        gPGy(j2)=gPGy(j2)+lPGy(lRow,lCol)
      enddo

      !6x3
      do lCol=1,3        
        j2=ele6x3(i,(lRow-1)*3+lCol)
        gGx(j2)=gGx(j2)+lGx(lRow,lCol)
        gGy(j2)=gGy(j2)+lGy(lRow,lCol)
      enddo
    enddo
    !$OMP END CRITICAL

    ! !6x3
    ! do lRow=1,6
    !   do lCol=1,3        
    !     j2=ele6x3(i,(lRow-1)*3+lCol)
    !     gGx(j2)=gGx(j2)+lGx(lRow,lCol)
    !     gGy(j2)=gGy(j2)+lGy(lRow,lCol)
    !   enddo
    ! enddo

    ! !6x6
    ! do lRow=1,6
    !   gRow=n(lRow)
    !   k=(gRow-1)*ivq(0)
    !   nlinkq=linkq(k+1:k+ivq(0))
    !   do lCol=1,6
    !     gCol=n(lCol)
    !     do j=1,ivq(gRow)
    !       if(nlinkq(j).eq.gCol) goto 11
    !     enddo
    !     write(9,*)"[Err] node conn missing in Bsnq at",gRow
    !     stop
    !     11 gNAdv(k+j)=gNAdv(k+j)+lNAd(lRow,lCol)
    !     gPGx(k+j)=gPGx(k+j)+lPGx(lRow,lCol)
    !     gPGy(k+j)=gPGy(k+j)+lPGy(lRow,lCol)
    !   enddo
    ! enddo

    ! !6x3
    ! do lRow=1,6
    !   gRow=n(lRow)
    !   k=(gRow-1)*ivl(0)
    !   nlinkl=linkl(k+1:k+ivl(0))
    !   do lCol=1,3
    !     gCol=n(lCol)
    !     do j=1,ivl(gRow)
    !       if(nlinkl(j).eq.gCol) goto 12
    !     enddo
    !     write(9,*)"[Err] node conn missing in Bsnq at",gRow
    !     stop
    !     12 gGx(k+j)=gGx(k+j)+lGx(lRow,lCol)
    !     gGy(k+j)=gGy(k+j)+lGy(lRow,lCol)
    !   enddo
    ! enddo

    ! !3x6    
    ! do lRow=1,3
    !   gRow=n(lRow)
    !   k=(gRow-1)*ivq(0)
    !   nlinkq=linkq(k+1:k+ivq(0))
    !   do lCol=1,6
    !     gCol=n(lCol)
    !     do j=1,ivq(gRow)
    !       if(nlinkq(j).eq.gCol) goto 12
    !     enddo
    !     write(9,*)"[Err] node conn missing in C Flux at",gRow
    !     stop
    !     12 gCxFlux(k+j)=gCxFlux(k+j)+cxFlux(lRow,lCol)
    !     gCyFlux(k+j)=gCyFlux(k+j)+cyFlux(lRow,lCol)
    !   enddo
    ! enddo    
    
    ! !3x3
    ! do lRow=1,3
    !   gRow=n(lRow)
    !   k=(gRow-1)*ivl(0)
    !   nlinkl=linkl(k+1:k+ivl(0))
    !   do lCol=1,3
    !     gCol=n(lCol)
    !     do j=1,ivl(gRow)
    !       if(nlinkl(j).eq.gCol) goto 14
    !     enddo
    !     write(9,*)"[Err] node conn missing in dMat at",gRow
    !     stop
    !     14 mass2(k+j)=mass2(k+j)+(locM2(lRow,lCol)*invJ(i,5))
    !     gDMat(k+j)=gDMat(k+j)+dMat(lRow,lCol)        
    !   enddo
    ! enddo

  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

end subroutine matrixSet2
