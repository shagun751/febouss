subroutine nodeConnVSR(npoin,nelem,maxNePoi,maxNeEle,conn,&
  poi2poi,poi2ele,npoisur)

  ![Note]:
  ! This one has connectivity matrix with each row for a node
  ! containing nodes excluding itself. This is because in the
  ! VSR storage the value corresponding to the node itself in
  ! the node row of the matrix has to be placed at the end.
  ! Hence excluding it from this table makes it easier so that
  ! this table can be used directly for both linear and quadratic
  ! linkTables and then node itself can be added while looping
  ! for the linkTable

implicit none
  integer(kind=4),intent(in)::npoin,nelem,maxNePoi,maxNeEle
  integer(kind=4)::i,j,k,l,n(6),i2,j2,k2,l2
  integer(kind=4),intent(in)::conn(nelem,6)
  integer(kind=4),intent(out)::poi2poi(npoin,maxNePoi),poi2ele(npoin,maxNeEle),npoisur(npoin,3)

  !Initialisations
  poi2poi=0
  poi2ele=0
  npoisur=0

  do i=1,nelem
    n=(/ (conn(i,k),k=1,6) /)
    do j=1,6      
      do l=1,6
        if(n(j).ne.n(l)) then
          do i2=1,maxNePoi
            if(poi2poi(n(j),i2).eq.0) then
              npoisur(n(j),1)=npoisur(n(j),1)+1
              poi2poi(n(j),i2)=n(l)                        
              exit    
            elseif(poi2poi(n(j),i2).eq.n(l)) then
              exit
            elseif(i2.eq.maxNePoi) then
              write(9,*) "[Error] Increase maxNePoi at node ",n(j)
              stop
            endif      
          enddo
        endif
      enddo
      npoisur(n(j),2)=npoisur(n(j),2)+1
      if(npoisur(n(j),2).gt.maxNeEle) then
        write(9,*) "[Error] Increase maxNeEle at node ",n(j)
        stop
      else  
        poi2ele(n(j),npoisur(n(j),2))=i
      endif
    enddo
  enddo

end subroutine nodeConnVSR

subroutine nodeConnInit(npoin,nelem,maxNeEle,conn,poi2ele)

  ![Note]:
  ! This one has connectivity matrix with each row for a node
  ! containing nodes excluding itself. This is because in the
  ! VSR storage the value corresponding to the node itself in
  ! the node row of the matrix has to be placed at the end.
  ! Hence excluding it from this table makes it easier so that
  ! this table can be used directly for both linear and quadratic
  ! linkTables and then node itself can be added while looping
  ! for the linkTable

implicit none
  integer(kind=4),intent(in)::npoin,nelem,maxNeEle
  integer(kind=4)::i,j,k,l,n(3),i2,j2,k2,l2
  integer(kind=4),intent(in)::conn(nelem,6)
  integer(kind=4),intent(out)::poi2ele(npoin,maxNeEle)
  integer(kind=4)::npoisur(npoin)

  !Initialisations  
  poi2ele=0
  npoisur=0

  do i=1,nelem
    n=conn(i,1:3)
    do j=1,3      
      npoisur(n(j))=npoisur(n(j))+1
      if(npoisur(n(j)).gt.maxNeEle) then
        write(9,*) "[Error] Increase maxNeEle at node ",n(j)
        stop
      else  
        poi2ele(n(j),npoisur(n(j)))=i
      endif
    enddo
  enddo

end subroutine nodeConnInit

subroutine nodeConnYale(npoin,nelem,maxNePoi,maxNeEle,conn,poi2poi,poi2ele,npoisur)

  ![Note]:
  ! This one has connectiity matrix with each row for a node
  ! containing nodes including itslef. This is because Yale
  ! sotrage doesnt require any special storage for the value
  ! corresponding to the node itself in the row in the matrix
  ! unlike the VSR storage.

implicit none
  integer(kind=4),intent(in)::npoin,nelem,maxNePoi,maxNeEle
  integer(kind=4)::i,j,k,l,n(6),i2,j2,k2,l2
  integer(kind=4),intent(in)::conn(nelem,6)
  integer(kind=4),intent(out)::poi2poi(npoin,maxNePoi),poi2ele(npoin,maxNeEle),npoisur(npoin,3)

  !Initialisations
  poi2poi=0
  poi2ele=0
  npoisur=0

  do i=1,nelem
    n=(/ (conn(i,k),k=1,6) /)
    do j=1,6      
      do l=1,6        
        do i2=1,maxNePoi
          if(poi2poi(n(j),i2).eq.0) then
            npoisur(n(j),1)=npoisur(n(j),1)+1
            poi2poi(n(j),i2)=n(l)                        
            exit    
          elseif(poi2poi(n(j),i2).eq.n(l)) then
            exit
          elseif(i2.eq.maxNePoi) then
            write(9,*) "[Error] Increase maxNePoi at node ",n(j)
            stop
          endif      
        enddo
      enddo
      npoisur(n(j),2)=npoisur(n(j),2)+1
      if(npoisur(n(j),2).gt.maxNeEle) then
        write(9,*) "[Error] Increase maxNeEle at node ",n(j)
        stop
      else  
        poi2ele(n(j),npoisur(n(j),2))=i
      endif
    enddo
  enddo

end subroutine nodeConnYale