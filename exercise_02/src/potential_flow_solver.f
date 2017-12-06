!> \file       potential_flow_solver.f
!> \brief      2d solver for potential flows
!> \author     Immo Huismann
!> \date       2014/01/09
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Potential_Flow_Solver
  use Kind_Parameters, only: RDP
  use Constants,       only: ZERO, ONE, TWO, PI
  use Helmholtz_2D
  implicit none

  !.............................................................................
  ! GLL variables in standard element

  real(RDP), allocatable :: xi(:)          ! collocation points in [-1,1]
  real(RDP), allocatable :: w(:)           ! integration weights
  real(RDP), allocatable :: diff(:,:)      ! differentiation matrix
  real(RDP), allocatable :: stiff(:,:)     ! stiffness matrix
  real(RDP), allocatable :: id(:)          ! identity matrix

  !.............................................................................
  ! field variables

  real(RDP), allocatable :: x  (:,:,:,:,:) ! domain coordinates
  real(RDP), allocatable :: u  (:,:,:,:,:) ! velocity
  real(RDP), allocatable :: phi(:,:,:,:  ) ! potential to solve for
  real(RDP), allocatable :: err(:,:,:,:  ) ! solution error
  real(RDP), allocatable :: f  (:,:,:,:)   ! discrete RHS
  real(RDP), allocatable :: dbc(:,:,:,:)   ! Dirichlet BC

  !.............................................................................
  ! parameters

  integer   :: po = 2                      ! polynomial degrees
  integer   :: n_element(2) = [3, 3]       ! number of elements
  logical   :: cyclic_bc(2)                ! switch for cyclic BC

  !.............................................................................
  ! element widths & loop bounds

  real(RDP) :: h(2)                        ! element width
  integer   :: i,j,m,o                     ! loop indices
  integer   :: uni                         ! output unit

  !.............................................................................
  ! IO parameters

  character(len=*), parameter ::  INPUT_FILE = 'input'
  character(len=*), parameter :: OUTPUT_FILE = 'output'

  !............................................................................
  ! read polynomial degree and number of elements

  open(newunit = uni, file = INPUT_FILE)
  read(uni,*) po
  read(uni,*) n_element(1)
  read(uni,*) n_element(2)
  close(uni)

  !.............................................................................
  ! Allocation and Initialization

  call GetGLLVariables(po, xi, w, diff, stiff)

  allocate(id(0:po), source = ONE)

  allocate(x  (0:po,0:po,1:n_element(1),1:n_element(2),2), source = ZERO)
  allocate(u  (0:po,0:po,1:n_element(1),1:n_element(2),2), source = ZERO)
  allocate(phi(0:po,0:po,1:n_element(1),1:n_element(2)  ), source = ZERO)
  allocate(err(0:po,0:po,1:n_element(1),1:n_element(2)  ), source = ZERO)
  allocate(f  (0:po,0:po,1:n_element(1),1:n_element(2)  ), source = ZERO)
  allocate(dbc(0:po,0:po,1:n_element(1),1:n_element(2)  ), source = ZERO)

  ! Set the domain to map from [-1;1]^2
  h = TWO / n_element

  x = CartesianGrid(xi,[-ONE,-ONE],n_element,h)

  !.............................................................................
  ! Setup the test case

  ! setup a mask array. 1: no Dirichlet BC, 0: Dirichlet BC
  dbc                      = 1
  dbc( 0,:,           1,:) = 0
  dbc(po,:,n_element(1),:) = 0
  dbc(:, 0,:,           1) = 0
  dbc(:,po,:,n_element(2)) = 0

  ! Set potential to the correct values on the boundary for the Dirichlet BC
  phi = StagnationFlowPotential([-ONE,-ONE],x(:,:,:,:,1),x(:,:,:,:,2))
  phi = phi * (1 - dbc)

  ! the discrete RHS is zero as a homogeneous Laplace equation is solved
  ! no Neumann boundary conditions added
  f = 0

  cyclic_bc = [.false.,.false.]

  !.............................................................................
  ! Solution process

  ! A scalar laplace equation is solved using the CG solver
  call ConjugatedGradients(ZERO,h,w,stiff,cyclic_bc,dbc,f,phi)

  !.............................................................................
  ! Calculate the velocities from the potential

  ! derivative in first direction for the first velocity
  u(:,:,:,:,1) = TensorProduct2D(id  ,diff*h(2)/h(1),phi)

  ! derivative in second direction for the second velocity
  u(:,:,:,:,2) = TensorProduct2D(diff/h(2)*h(1),id  ,phi)

  ! ...........................................................................
  ! Output

  err = phi - StagnationFlowPotential([-ONE,-ONE],x(:,:,:,:,1),x(:,:,:,:,2))

  print '(A,E22.7)', 'Maximum error    : ',maxval(abs(err))

  open(newunit = uni, file = OUTPUT_FILE)
  write(uni, *) po, po, n_element,0,0
  do j = 1, n_element(2)
  do o = 0, po
  do i = 1, n_element(1)
  do m = 0, po
    write(uni,'(6(1X,E22.7))') &
         x(m,o,i,j,:), u(m,o,i,j,:),phi(m,o,i,j),err(m,o,i,j)
  end do
  end do
  end do
  end do
  close(uni)

contains

!-------------------------------------------------------------------------------
!> \brief   The potential of the stagnation flow
!> \author  Immo Huismann

pure function StagnationFlowPotential(center,x_1,x_2) result(phi)
  real(RDP), intent(in) :: center(2) !< stagnation point
  real(RDP), intent(in) :: x_1(:,:,:,:) !< first  coordinate
  real(RDP), intent(in) :: x_2(:,:,:,:) !< second coordinate
  real(RDP) :: phi(size(x_1,1),size(x_1,2),size(x_1,3),size(x_1,4))

  phi = &
       + (x_1 - center(1)) ** 3 * (x_2 - center(2)) &
       - (x_1 - center(1))      * (x_2 - center(2)) ** 3

end function StagnationFlowPotential

!===============================================================================

end program Potential_Flow_Solver
