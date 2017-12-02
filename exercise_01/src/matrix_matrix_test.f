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

  use Matrix_Matrix_Simple
  use Matrix_Matrix_Blocked

  use MPI_F08
  implicit none

  integer,   parameter   :: N = 200            ! number of intervals
  integer,   parameter   :: N_REPITITION = 100 ! number of repititions for test

  real(RDP), parameter   :: FLOP          = 2 * real(N, RDP)**3
  real(RDP), parameter   :: FLOP_TO_GFLOP = 1e-9_RDP
  real(RDP), parameter   :: ALLOWED_ERROR = 1e-9_RDP
  real(RDP), allocatable :: A(:,:), B(:,:), C(:,:), A_ref(:,:)
  integer :: i, j

  integer :: i_repitition

  real(RDP) :: start_time, stop_time, time_taken, time_per_application

  !.............................................................................
  ! Allocation of data

  allocate(A(N,N),     source = ZERO)
  allocate(B(N,N),     source = ZERO)
  allocate(C(N,N),     source = ZERO)
  allocate(A_ref(N,N), source = ZERO)
  
  !.............................................................................
  ! Initialization

  call MPI_Init()

  do j = 1, N
  do i = 1, N
    B(i,j) = 5 * i - 3 * j
    C(i,j) = i + j
  end do
  end do

  A_ref = matmul(B, C)

  !.............................................................................
  ! Baseline implementation

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    A = matmul(B, C)
  end do

  stop_time            = MPI_Wtime()

  time_taken           = stop_time - start_time
  time_per_application = time_taken / N_REPITITION

  print *, 'Baseline implementation'    
  if( maxval(abs(A_ref - A)) > ALLOWED_ERROR) then
    print '(A,1X,E15.7)', 'Error:  ', maxval(abs(A_ref - A))
  end if
  print '(A,1X,E15.7)', 'Time:   ',                        time_per_application
  print '(A,1X,E15.7)', 'GFLOPS: ', FLOP * FLOP_TO_GFLOP / time_per_application
  print *

  !.............................................................................
  ! Direct implementation

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call MatrixMatrixSimple(N, A, B, C)
  end do

  stop_time            = MPI_Wtime()

  time_taken           = stop_time - start_time
  time_per_application = time_taken / N_REPITITION

  print *, 'Direct implementation'         
  if( maxval(abs(A_ref - A)) > ALLOWED_ERROR) then
    print '(A,1X,E15.7)', 'Error:  ', maxval(abs(A_ref - A))
  end if
  print '(A,1X,E15.7)', 'Time:   ',                        time_per_application
  print '(A,1X,E15.7)', 'GFLOPS: ', FLOP * FLOP_TO_GFLOP / time_per_application
  print *

  !.............................................................................
  ! Blocked variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call MatrixMatrixBlocked(N, A, B, C)
  end do

  stop_time            = MPI_Wtime()

  time_taken           = stop_time - start_time
  time_per_application = time_taken / N_REPITITION

  print *, 'Blocked implementation'        
  if( maxval(abs(A_ref - A)) > ALLOWED_ERROR) then
    print '(A,1X,E15.7)', 'Error:  ', maxval(abs(A_ref - A))
  end if
  print '(A,1X,E15.7)', 'Time:   ',                        time_per_application
  print '(A,1X,E15.7)', 'GFLOPS: ', FLOP * FLOP_TO_GFLOP / time_per_application
  print *

  !.............................................................................
  ! Cleanup

  call MPI_Finalize()

!===============================================================================

end program Stencil_Test
