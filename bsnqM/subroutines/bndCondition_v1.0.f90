!!---------------------------diriBCMass----------------------------!!
subroutine diriBCMass(npl,npt,nbndp,bndP,bndPT,Sz,ivl,ivq,&
 linkl,linkq,bndPN,gBs1,gBs2,gBs3,gBs4,massE)
use bsnqGlobVars
implicit none

  integer(kind=C_K1),intent(in)::npl,npt,nbndp,Sz(4)
  integer(kind=C_K1),intent(in)::ivl(0:npt),ivq(0:npt)
  integer(kind=C_K1),intent(in)::linkl(Sz(3))
  integer(kind=C_K1),intent(in)::linkq(Sz(4))
  integer(kind=C_K1),intent(in)::bndP(nbndp),bndPT(nbndp)
  integer(kind=C_K1)::i,j,k,i1,j1,k1
  
  real(kind=C_K2),intent(in)::bndPN(npt,2)
  real(kind=C_K2),intent(inout)::gBs1(Sz(4)),gBs2(Sz(4))
  real(kind=C_K2),intent(inout)::gBs3(Sz(4)),gBs4(Sz(4))
  real(kind=C_K2),intent(inout)::massE(Sz(1))

  !! BndCond vel eta SemiDirect
  do i=1,nbndp
    i1=bndP(i)
    j1=bndPT(i)
    k=(i1-1)*ivq(0)
    k1=ivq(i1)
    if((j1.eq.11).or.(j1.eq.12).or.(j1.eq.14))then
      gBs1(k+k1)=1d0
      gBs1(k+1:k+k1-1)=0d0
      gBs2(k+1:k+k1)=0d0
      
      gBs4(k+k1)=1d0
      gBs4(k+1:k+k1-1)=0d0      
      gBs3(k+1:k+k1)=0d0
    
    elseif(j1.eq.13)then
      if( abs(bndPN(i1,1)) .gt. 0.86d0 )then !degLT30
        gBs1(k+k1)=1d0
        gBs1(k+1:k+k1-1)=0d0
        gBs2(k+1:k+k1-1)=0d0
        gBs2(k+k1)=bndPN(i1,2)/bndPN(i1,1)
        !write(*,*) '!!!!000',i1
      
      else        
        gBs4(k+k1)=1d0
        gBs4(k+1:k+k1-1)=0d0      
        gBs3(k+1:k+k1-1)=0d0
        gBs3(k+k1)=bndPN(i1,1)/bndPN(i1,2)
        !write(*,*) '!!!!555',i1
      endif

    ! else
    !   write(9,'(" [ERR] Not coded bndType",I10)')bndPT(i)
    !   stop
    endif

    ! Linear nodes only
    if(i1.gt.npl)cycle
    k=(i1-1)*ivl(0)
    k1=ivl(i1)
    if((j1.eq.11).or.(j1.eq.14))then
      massE(k+k1)=1d0
      massE(k+1:k+k1-1)=0d0
    endif
  enddo

end subroutine diriBCMass
!!-------------------------End diriBCMass--------------------------!!