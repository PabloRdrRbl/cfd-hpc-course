program Calculate_Pi
  use Kind_Parameters ! Higher precision kinds
  use Constants       ! constants like PI
  implicit none

  integer :: i
  integer :: max_n = 1E4
  real(RNP) :: pi_estimation = 0.0
      
  do i=0,max_n
      pi_estimation = pi_estimation + ((-1.0)**i / (2.0*i + 1.0))
  end do
      
  pi_estimation = 4 * pi_estimation
  
  write(*,*) 'PI estimation: ', pi_estimation
    
  write(*,*) 'PI error: ', (PI - pi_estimation)
      
end program Calculate_Pi
