subroutine mergeSort(maxNePoi,nNePoi,poi2poi)
implicit none
  integer(kind=4),intent(in)::maxNePoi,nNePoi
  integer(kind=4),intent(inout)::poi2poi(maxNePoi)
  integer(kind=4)::temp(maxNePoi),i,j,k,l,m
  integer(kind=4)::istart,imiddle,iend,width

  temp=0
  width=1
  do while(width.lt.nNePoi)
    do i=1,nNePoi,(2*width)
      istart=i
      imiddle=min(i+width-1,nNePoi)
      iend=min(i+2*width-1,nNePoi)

      j=istart
      k=imiddle+1
      do l=istart,iend
        if((j.le.imiddle).and.(k.gt.iend .or. poi2poi(j).le.poi2poi(k))) then
          temp(l)=poi2poi(j)
          j=j+1
        else
          temp(l)=poi2poi(k)
          k=k+1
        endif
      enddo  
      !write(9,*)(temp(m),m=1,maxNePoi)
    enddo    
    poi2poi=temp
    width=2*width
  enddo

end subroutine mergeSort