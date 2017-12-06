!> \file       test_interpolation.f
!> \brief      A runtime test for the interpolation on spectral elements
!> \author     Immo Huismann
!> \date       2016/11/29
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Test_Interpolation
  use Kind_Parameters,                only: RDP
  use Spectral_Element_Interpolation

  use MPI_F08
  implicit none

  !.............................................................................
  ! Field sizes

  integer :: po                             !< polynomial order
  integer :: ne_x1                          !< number of elements, x_1 direction
  integer :: ne_x2                          !< number of elements, x_2 direction
  integer :: ni                             !< number of nodes to interpolate to

  !.............................................................................
  ! Arrays

  real(RDP), allocatable :: A(:,:,:,:)      !< 2D to 2D interpolation matrix
  real(RDP), allocatable :: u  (:,:,:,:)    !< variable to interpolate
  real(RDP), allocatable :: u_0(:,:,:,:)    !< initial value of u
  real(RDP), allocatable :: v  (:,:,:,:)    !< result of interpolation

  !.............................................................................
  ! measurements
  
  real(RDP) :: start_time, stop_time                  ! time stamps
  real(RDP) :: time_per_application, mlups            ! resulting measurements
  real(RDP) :: max_error                              ! error in implementation

  !.............................................................................
  ! input, loop variables, etc.

  character(len=*), parameter :: INPUT_FILE = 'input' ! input file to read from
  integer,          parameter :: N_REPITITION = 100   ! number of repititions
  
  integer :: uni, i_repitition, n_dof
  integer :: i, j

  !.............................................................................
  ! read input

  open(newunit = uni, file = INPUT_FILE)
  read(uni,*) po
  read(uni,*) ne_x1
  read(uni,*) ne_x2  
  close(uni)
  
  ni    = po + 1
  n_dof = ni**3 * ne_x1 * ne_x2

  !.............................................................................
  ! allocate fields

  allocate(A  (1:ni,1:ni,0:po,0:po))
  allocate(u  (0:po,0:po,1:ne_x1,1:ne_x2))
  allocate(u_0(0:po,0:po,1:ne_x1,1:ne_x2))
  allocate(v  (1:ni,1:ni,1:ne_x1,1:ne_x2))

  call Random_Number(u)
  u_0 = u

  ! Initialize A as the identity matrix
  A = 0
  do j = 1, ni
  do i = 1, ni
    A(i,j,i-1,j-1) = 1
  end do
  end do

  !.............................................................................
  ! Run the interpolation N_REPITITION times
  !
  ! The test is divided into two parts to ensure that the subroutine call is not
  ! optimized out of the program.
  
  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION

    if(mod(i_repitition,2) == 1) then
      call SpectralElementInterpolation(A, u, v)
    else
      call SpectralElementInterpolation(A, v, u)
    end if

  end do

  stop_time = MPI_Wtime()
  
  !.............................................................................
  ! Print the time per application and the Mega Lattice Updates Per Second

  max_error = maxval(abs(u-u_0))

  time_per_application = (stop_time - start_time) / N_REPITITION
  mlups                = n_dof / time_per_application / (10**6)

  print '(A,1X,I15  )', 'Number of degrees of freedom: ', n_dof
  print '(A,1X,E15.7)', 'Time per application [s]:     ', time_per_application
  print '(A,1X,E15.7)', 'MLUPS:                        ', mlups
  print '(A,1X,E15.7)', 'Maximum error:                ', max_error

!===============================================================================

end program Test_Interpolation
