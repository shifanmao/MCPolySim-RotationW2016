!---------------------------------------------------------------*
      
PROGRAM test

DOUBLE PRECISION U(10,3)  ! Tangent vectors

do i=1:m
  do j=1:n
    write(*,"(A)",advance="no") U(i,j)
  end do
end do

END PROGRAM
      
!---------------------------------------------------------------*
