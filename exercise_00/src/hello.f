program hello

use constants
use kind_parameters
      
      
implicit none

real :: numbers(5)
integer :: i
real(RNP) :: MoreNumbers(5)
      
write(*,*) 'Howdy World :)'
write(*,*) PI

do i=1,5

      numbers(i) = i
      write(*,*) numbers(i)

      MoreNumbers(i) = sin(numbers(i))
      write(*,*) MoreNumbers(i)

end do

      
      
      
         
      

end program hello
