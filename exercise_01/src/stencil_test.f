!> \file      stencil_test.f
!> \brief     A small program running a stencil operation
!> \author    Immo Huismann
!> \date      2016/11/21
!> \copyright Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================

program Stencil_Test
  use Kind_Parameters,  only: RDP
  use Constants,        only: ZERO, ONE, PI

  use Stencil_Baseline, only: StencilBaseline
  use Stencil_Direct,   only: StencilDirect
  use Stencil_Blocked,  only: StencilBlocked
  

  use MPI_F08
  implicit none

  integer,   parameter   :: N = 200             ! number of intervals
  real(RDP), parameter   :: L = ONE             ! domain length
  real(RDP), parameter   :: H = L / (N - 1)     ! interval width
  integer,   parameter   :: N_REPITITION = 1000 ! number of repititions for test  

  real(RDP), allocatable :: u(:,:), v(:,:), x(:,:), y(:,:)
  integer :: i, j

  integer :: i_repitition

  real(RDP) :: start_time, stop_time, time_taken, time_per_application

  !.............................................................................
  ! Allocation of data
  
  allocate(u(N,N), source = ZERO)
  allocate(v(N,N), source = ZERO)
  allocate(x(N,N), source = ZERO)
  allocate(y(N,N), source = ZERO)

  !.............................................................................
  ! Initialization

  call MPI_Init()
  
  do j = 1, N
  do i = 1, N
    x(i,j) = (i - 1) * h
    y(i,j) = (j - 1) * h
  end do
  end do

  u = sin(2 * PI * x) * cos(2 * PI * x)

  !.............................................................................
  ! Baseline implementation

  start_time = MPI_Wtime()
  
  do i_repitition = 1, N_REPITITION
    call StencilBaseline(N,h,u,v)
  end do

  stop_time            = MPI_Wtime()

  time_taken           = stop_time - start_time
  time_per_application = time_taken / N_REPITITION
  
  print '(A,1X,E15.7)', 'Time taken:',            time_per_application
  print '(A,1X,E15.7)', 'MLUPS:     ', (N-1)**2 / time_per_application
  print *

  !.............................................................................
  ! Direct implementation

  start_time = MPI_Wtime()
  
  do i_repitition = 1, N_REPITITION
    call StencilDirect(N,h,u,v)
  end do

  stop_time            = MPI_Wtime()
  
  time_taken           = stop_time - start_time
  time_per_application = time_taken / N_REPITITION
  
  print '(A,1X,E15.7)', 'Time taken:',            time_per_application
  print '(A,1X,E15.7)', 'MLUPS:     ', (N-1)**2 / time_per_application
  print *

  !.............................................................................
  ! Blocked variant

  start_time = MPI_Wtime()
  
  do i_repitition = 1, N_REPITITION
    call StencilBlocked(N,h,u,v)
  end do

  stop_time            = MPI_Wtime()

  time_taken           = stop_time - start_time
  time_per_application = time_taken / N_REPITITION
  
  print '(A,1X,E15.7)', 'Time taken:',            time_per_application
  print '(A,1X,E15.7)', 'MLUPS:     ', (N-1)**2 / time_per_application
  print *

  !.............................................................................
  ! Cleanup
  
  call MPI_Finalize()

!===============================================================================

end program Stencil_Test
