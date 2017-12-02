!> \file      poutpouri_a1.f
!> \brief     A small test program
!> \author    Immo Huismann
!> \date      2016/11/22
!> \copyright Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================

program Poutpouri_A1
  use Kind_Parameters,       only: RDP

  use Potpourri_Print_MLUPS, only: PrintMLUPS
  use MPI_F08
  implicit none

  integer, parameter :: N            = 2000
  integer, parameter :: N_DOF        = N**2
  integer, parameter :: N_REPITITION = 100

  real(RDP), allocatable :: u(:,:), v(:,:)
  real(RDP) :: start_time, stop_time
  integer   :: i_repitition

  !.............................................................................
  ! Allocation of data

  allocate(u(N,N), v(N,N))

  !.............................................................................
  ! Initialization

  u = 7
  v = 0
  call MPI_Init()

  !.............................................................................
  ! Test with baseline variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call Baseline(N,u,v)
  end do

  stop_time            = MPI_Wtime()

  call PrintMLUPS(n_dof,n_repitition,start_time,stop_time)

  !.............................................................................
  ! Test with improved variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call Improved(N,u,v)
  end do

  stop_time            = MPI_Wtime()

  call PrintMLUPS(n_dof,n_repitition,start_time,stop_time)

  !.............................................................................
  ! Cleanup

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> \brief   Copy operator -- Baseline
!> \author  Immo Huismann

subroutine Baseline(n_point,u,v)
  integer,   intent(in)  :: n_point
  real(RDP), intent(in)  :: u(n_point,n_point)
  real(RDP), intent(out) :: v(n_point,n_point)

  integer :: i,j

  do i = 1, n_point
  do j = 1, n_point
    v(i,j) = u(i,j)
  end do
  end do

end subroutine Baseline

!-------------------------------------------------------------------------------
!> \brief   Copy operator -- Improved
!> \author  Immo Huismann

subroutine Improved(n_point,u,v)
  integer,   intent(in)  :: n_point
  real(RDP), intent(in)  :: u(n_point,n_point)
  real(RDP), intent(out) :: v(n_point,n_point)

end subroutine Improved

!===============================================================================

end program Poutpouri_A1
