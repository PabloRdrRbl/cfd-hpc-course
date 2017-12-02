!> \file      potpouri_a3.f
!> \brief     3D interpolation operator for Spectral-Element / DG methods
!> \author    Immo Huismann
!> \date      2016/11/22
!> \copyright Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================

program Potpouri_A3
  use Kind_Parameters, only: RDP
  use Constants,       only: ONE, ZERO

  use Potpourri_Print_MLUPS, only: PrintMLUPS
  use MPI_F08
  implicit none

  integer, parameter :: N            = 100000
  integer, parameter :: N_DOF        = N
  integer, parameter :: N_REPITITION = 100

  real(RDP), allocatable :: u(:), v(:), w(:)
  real(RDP) :: start_time, stop_time
  integer   :: i_repitition
  real(RDP) :: h

  !.............................................................................
  ! Allocation of data

  allocate(u(N), v(N), w(N))

  !.............................................................................
  ! Initialization

  u = 7
  v = [((-1)**i_repitition, i_repitition = 1, N)]
  h = ONE / (N - 1)

  call MPI_Init()

  !.............................................................................
  ! Test with baseline variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call Baseline(N,h,u,v,w)
  end do

  stop_time            = MPI_Wtime()

  call PrintMLUPS(n_dof,n_repitition,start_time,stop_time)

  !.............................................................................
  ! Test with improved variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call Improved(N,h,u,v,w)
  end do

  stop_time            = MPI_Wtime()

  call PrintMLUPS(n_dof,n_repitition,start_time,stop_time)

  !.............................................................................
  ! Cleanup

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> \brief   1D upwind scheme for convection -- Baseline
!> \author  Immo Huismann

subroutine Baseline(n,h,u,v,w)
  integer,   intent(in)  :: n    !< points per direction
  real(RDP), intent(in)  :: h    !< mesh width
  real(RDP), intent(in)  :: u(n) !< values to convect
  real(RDP), intent(in)  :: v(n) !< convection velocity
  real(RDP), intent(out) :: w(n) !< convection term

  integer :: i

  do i = 2, n-1

    if(v(i) > 0) then
      w(i) = (u(i  ) - u(i-1)) / h
    else
      w(i) = (u(i+1) - u(i  )) / h
    end if

  end do

end subroutine Baseline

!-------------------------------------------------------------------------------
!> \brief   1D upwind scheme for convection -- Improved
!> \author  Immo Huismann

subroutine Improved(n,h,u,v,w)
  integer,   intent(in)  :: n    !< points per direction
  real(RDP), intent(in)  :: h    !< mesh width
  real(RDP), intent(in)  :: u(n) !< values to convect
  real(RDP), intent(in)  :: v(n) !< convection velocity
  real(RDP), intent(out) :: w(n) !< convection term

end subroutine Improved

!===============================================================================

end program Potpouri_A3
