program Calculate_e
  use Kind_Parameters ! Higher precision kinds
  use Constants       ! constants like PI
  implicit none

  integer :: i
  integer :: n_max = 6
  real(RNP) :: e_estimation = 1.0
  real(RNP) :: nn = 1.0    

  do i=1,n_max
    nn = i * nn      
    e_estimation = e_estimation + (1.0 / nn)
  end do

  write(*,*) 'e estimation: ', e_estimation
    
  write(*,*) 'e error: ', (exp(1.0) - e_estimation)
      
      
end program Calculate_e
