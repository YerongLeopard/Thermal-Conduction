
program Main
INTEGER::N1=6,N2=6,N3,i,j,k,s=0
INTEGER::first1, first2, first3,end1,end2,end3
N3=6
first1=-N1+1
first2=-N2+1
first3=-N3+1
end1=N1
end2=N2
end3=N3
do N3=6, 100
do i=first1,end1
  do j=first2,end2
    do k=first3,end3
      if (((-N1<2*i).AND.(2*i<N1+1)).AND.((-N2<2*j).AND.(2*j<N2+1)).AND.((-N3<2*k).AND.(2*k<N3+1))) then
         s=s+1
!result=result+((v**2)/sq)*((Debyek/k)**p)/(1+4*((np.sin(q0/2))**2)*(MinPath*(v/k)*((Debyek/k)**p)+(MinPath*(v/k)*((Debyek/k)**p))**2))
!result=result/(N3*N1*N2*4)*(3-p)*((np.cos(q0/2))**2)
      endif
    enddo
  enddo
enddo
if (s /= N1*N2*N3) then
  write(*,*) s,N1*N2*N3
endif
enddo
end program Main
