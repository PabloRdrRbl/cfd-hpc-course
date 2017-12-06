!> \file       helmholtz_2d.f
!> \brief      Routines for the two-dimensional Helmholtz equation
!> \author     Immo Huismann
!> \date       2014/01/09
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Helmholtz_2D
  use Kind_Parameters, only: RDP
  use Constants,       only: ZERO
  use Gauss_Jacobi,    only: GLL_Points, GLL_Weights, GLL_DiffMatrix
  implicit none

  private

  public :: GetGLLVariables
  public :: AssembleResidual
  public :: TensorProduct2D
  public :: HelmholtzResidual
  public :: HelmholtzOperator
  public :: ConjugatedGradients
  public :: CartesianGrid

  interface TensorProduct2D
    module procedure TensorProduct2D_DxM
    module procedure TensorProduct2D_MxD
    module procedure TensorProduct2D_DxD
  end interface

contains

!-------------------------------------------------------------------------------
!> \brief   This routine sets the GLL Variables
!> \author  Immo Huismann
!>
!> \details
!> This routine gets the node locations, the integration weights, and the
!> differentiation matrix depending on the polynomial degree p.

subroutine GetGLLVariables(po, xi, wi, diff, stiff)
  integer,                intent(in)  :: po         !< desired polynomial degree
  real(RDP), allocatable, intent(out) :: xi   (:)   !< node locations
  real(RDP), allocatable, intent(out) :: wi   (:)   !< weights
  real(RDP), allocatable, intent(out) :: diff (:,:) !< differentiation matrix
  real(RDP), allocatable, intent(out) :: stiff(:,:) !< differentiation matrix

  if(allocated(xi   )) deallocate(xi   )
  if(allocated(wi   )) deallocate(wi   )
  if(allocated(diff )) deallocate(diff )
  if(allocated(stiff)) deallocate(stiff)

  ! allocate and set the GLL variables for the polynomial degree po
  allocate(xi   (0:po)     , source = GLL_Points     (po))
  allocate(wi   (0:po)     , source = GLL_Weights    (xi))
  allocate(diff (0:po,0:po), source = GLL_DiffMatrix (xi))
  allocate(stiff(0:po,0:po), source = StiffnessMatrix(wi,diff))

end subroutine GetGLLVariables

!-------------------------------------------------------------------------------
!> \brief   This routine sets the GLL stiffness matrix
!> \author  Immo Huismann

pure function StiffnessMatrix(wi, diff) result(stiff)
  real(RDP), intent(in ) :: wi  (0:)    !< weights
  real(RDP), intent(in ) :: diff(0:,0:) !< differentiation matrix
  real(RDP) :: stiff(0:ubound(diff,1),0:ubound(diff,2)) !< stiffness matrix

  integer :: i, po
  real(RDP) :: mass(0:ubound(diff,1),0:ubound(diff,2)) ! mass matrix

  po = ubound(diff,1)

  mass = 0
  do i = 0, po
    mass(i,i) = wi(i)
  end do

  stiff = matmul(transpose(diff),matmul(mass,diff))

end function StiffnessMatrix

!-------------------------------------------------------------------------------
!> \brief   Returns a cartesian grid
!> \author  Immo Huismann

pure function CartesianGrid(xi,x_start,n_ele,h) result(x)
  real(RDP), intent(in) :: xi(0:)     !< 1D collocation point distribution
  real(RDP), intent(in) :: x_start(2) !< starting point
  integer  , intent(in) :: n_ele(2)   !< number of elements
  real(RDP), intent(in) :: h(2)
  real(RDP) :: x(0:ubound(xi,1),0:ubound(xi,1),n_ele(1),n_ele(2),2)

  integer :: i,j,m,o
  integer :: po

  po = ubound(xi,1)

  ! first dimension
  do j = 1, n_ele(2)
  do i = 1, n_ele(1)
  do o = 0, po
    x(:,o,i,j,1) = (xi + 1) / 2 * h(1) + (i-1) * h(1) + x_start(1)
  end do
  end do
  end do

  ! second dimension
  do j = 1, n_ele(2)
  do i = 1, n_ele(1)
  do m = 0, po
    x(m,:,i,j,2) = (xi + 1) / 2 * h(2) + (j-1) * h(2) + x_start(2)
  end do
  end do
  end do

end function CartesianGrid

!-------------------------------------------------------------------------------
!> \brief   Elementwise Helmholtz residual r_e = f_e - H u_e
!> \author  Immo Huismann

subroutine HelmholtzResidual(lambda, h, w, stiff, u, f, r)
  real(RDP), intent(in)  :: lambda          !< Helmholtz parameter
  real(RDP), intent(in)  :: h(2)            !< element width
  real(RDP), intent(in)  :: w(0:)           !< integration weights
  real(RDP), intent(in)  :: stiff(0:,0:)    !< stiffness matrix
  real(RDP), intent(in)  :: u(0:,0:,1:,1:)  !< values to operate on
  real(RDP), intent(in)  :: f(0:,0:,1:,1:)  !< discrete RHS
  real(RDP), intent(out) :: r(0:,0:,1:,1:)  !< discrete RHS

  integer :: n_element(2)
  integer :: po
  integer :: e_1, e_2, i, j

  po           = ubound(u,1)
  n_element(1) = ubound(r,3)
  n_element(2) = ubound(r,4)

  call HelmholtzOperator(lambda, h, w, stiff, u, r)

  do e_2 = 1, n_element(2)
  do e_1 = 1, n_element(1)
    do j = 0, po
    do i = 0, po
      r(i,j,e_1,e_2) = f(i,j,e_1,e_2) - r(i,j,e_1,e_2)
    end do
    end do
  end do
  end do

end subroutine HelmholtzResidual

!-------------------------------------------------------------------------------
!> Evaluates the Helmholtz Operator elementwise - no assembly yet

subroutine HelmholtzOperator(lambda, h, w, stiff, u, r)
  real(RDP), intent(in)  :: lambda          !< Helmholtz parameter
  real(RDP), intent(in)  :: h(2)            !< element width (cartesian meshes)
  real(RDP), intent(in)  :: w(0:)           !< GLL integration weights
  real(RDP), intent(in)  :: stiff(0:,0:)    !< stiffness matrix
  real(RDP), intent(in)  :: u(0:,0:,1:,1:)  !< values to operate on
  real(RDP), intent(out) :: r(0:,0:,1:,1:)  !< values to operate on

  integer   :: n_element(2)
  integer   :: po
  integer   :: i, j, k, l, m
  real(RDP) :: tmp_sum, fact(0:2)


  po           = ubound(u,1)
  n_element(1) = ubound(r,3)
  n_element(2) = ubound(r,4)

  ! compute scaling factors
  fact = [lambda * h(1) * h(2) / 4, h(2) / h(1), h(1) / h(2)]

  do l = 1, n_element(2)
  do k = 1, n_element(1)
    do j = 0, po
    do i = 0, po
      r(i,j,k,l) = w(i) * w(j) * u(i,j,k,l) * fact(0)
    end do
    end do

    do j = 0, po
    do i = 0, po
      tmp_sum = 0
      do m = 0, po
        tmp_sum = tmp_sum + stiff(i,m) * u(m,j,k,l)
      end do
      r(i,j,k,l) = r(i,j,k,l) + tmp_sum * w(j) * fact(1)
    end do
    end do

    do j = 0, po
    do i = 0, po
      tmp_sum = 0
      do m = 0, po
        tmp_sum = tmp_sum + stiff(j,m) * u(i,m,k,l)
      end do
      r(i,j,k,l) = r(i,j,k,l) + tmp_sum * w(i) * fact(2)
    end do
    end do

  end do
  end do

end subroutine HelmholtzOperator

!-------------------------------------------------------------------------------
!> Assembles a residual on a 2D structured mesh

pure subroutine AssembleResidual(cyclic_bc,res)
  logical  , intent(in   ) :: cyclic_bc(2)     !< are the dimensions cyclic?
  real(RDP), intent(inout) :: res(0:,0:,1:,1:) !< the residual to assemble

  integer :: po, n_element(2) ! polynomial degree and number of elements
  integer :: i, j             ! loop indices

  po        =  ubound(res,1)
  n_element = [ubound(res,3), ubound(res,4)]

  ! work on the first one
  do i = 2, n_element(1)
    res( 0,:,i  ,:) = res(0,:,i,:) + res(po,:,i-1,:)
    res(po,:,i-1,:) = res(0,:,i,:)
  end do

  ! are there cyclic BC?
  if(cyclic_bc(1)) then
    res( 0,:,          1, :) = res(0,:,1,:) + res(po,:,n_element(1),:)
    res(po,:,n_element(1),:) = res(0,:,1,:)
  end if

  ! work on the second one
  do j = 2, n_element(2)
    res(:, 0,:,j  ) = res(:,0,:,j) + res(:,po,:,j-1)
    res(:,po,:,j-1) = res(:,0,:,j)
  end do

  ! are there cyclic BC?
  if(cyclic_bc(2)) then
    res(:, 0,:,          1 ) = res(:,0,:,1) + res(:,po,:,n_element(2))
    res(:,po,:,n_element(2)) = res(:,0,:,1)
  end if

end subroutine AssembleResidual

!===============================================================================
! Solver

!-------------------------------------------------------------------------------
!> \brief   CG solver for the Helmholtz problem
!> \author  Immo Huismann

subroutine ConjugatedGradients(lambda,h,w,stiff,cyclic_bc,dbc,f,u)
  real(RDP), intent(in)    :: lambda           !< Helmholtz parameter
  real(RDP), intent(in)    :: h(2)             !< element width
  real(RDP), intent(in)    :: w(0:)            !< weights for first direction
  real(RDP), intent(in)    :: stiff(0:,0:)     !< stiffness matrix first dir
  logical  , intent(in)    :: cyclic_bc(2)     !< are there cyclic bc?
  real(RDP), intent(in)    :: dbc(0:,0:,1:,1:) !< Dirichlet BC
  real(RDP), intent(in)    :: f  (0:,0:,1:,1:) !< discrete RHS
  real(RDP), intent(inout) :: u  (0:,0:,1:,1:) !< variable to solve for

  real(RDP), parameter :: r_max    = 10.0_RDP**(-12)
  integer,   parameter :: iter_max = 2000

  real(RDP), allocatable, dimension(:,:,:,:) :: r,p,q,multi_inv
  real(RDP) :: alpha, beta, rho, rho_0
  integer   :: iter

  integer   :: po, n_element(2)

  po        =  ubound(f,1)
  n_element = [ubound(f,3),ubound(f,4)]

  allocate(r        (0:po,0:po,1:n_element(1),1:n_element(2)), source = ZERO)
  allocate(p        (0:po,0:po,1:n_element(1),1:n_element(2)), source = ZERO)
  allocate(q        (0:po,0:po,1:n_element(1),1:n_element(2)), source = ZERO)
  allocate(multi_inv(0:po,0:po,1:n_element(1),1:n_element(2)), source = ZERO)

  ! Calculate inverted multiplicity of local collocation points
  multi_inv = 1
  call AssembleResidual(cyclic_bc,multi_inv)
  multi_inv = 1 / multi_inv

  ! Get the residual of the Helmholtz equation, mind that the Neumann BC are needed
  call HelmholtzResidual(lambda,h,w,stiff,u,f,r)

  call AssembleResidual(cyclic_bc,r)

  r = r * dbc

  ! set the first basis vector (first search direction) as the residual
  p   = r

  ! exit if the residual is too low to achieve a better solution
  rho = WeightedScalarProduct(multi_inv,r,r)

  if(sqrt(rho) == 0) then

    return
  end if

  do iter = 1, iter_max

    !...........................................................................
    ! Calculate helmholtz operator

    ! get the impact of direction p on the residual, store in q
    call HelmholtzOperator(lambda,h,w,stiff,p,q)

    ! assemble q to ensure the continuity (global vs. local formulation)
    call AssembleResidual(cyclic_bc,q)

    ! Apply Dirichlet boundary conditions
    q = q * dbc

    !...........................................................................
    ! calculate update

    ! Get the distance to move in direction p
    alpha = rho / WeightedScalarProduct(multi_inv,q,p)

    ! advance distance alpha in direction p (q = H p for r)
    u     = u + alpha * p
    r     = r - alpha * q

    ! exit if the residual is too small

    rho_0 = rho
    rho   = WeightedScalarProduct(multi_inv,r,r)

    if(sqrt(rho) < r_max) exit

    !...........................................................................
    ! Calculate next search vector p

    ! get the residual reduction by p
    beta = rho / rho_0
    p    = r + beta * p

  end do

  print '(A,E22.7)', 'Residual:          ', sqrt(rho)
  print '(A,I22,A)', 'Nr. of iterations: ', iter - 1

end subroutine ConjugatedGradients

!-------------------------------------------------------------------------------
!> \brief   Computes a scalar product with a weight
!> \author  Immo Huismann

function WeightedScalarProduct(w,x,y) result(prod)
  real(RDP), intent(in) :: w(0:,0:,1:,1:) !< weights
  real(RDP), intent(in) :: x(0:,0:,1:,1:) !< first  input vector
  real(RDP), intent(in) :: y(0:,0:,1:,1:) !< second input vector
  real(RDP)             :: prod           !< resulting scalar product

  integer :: n_element(2)
  integer :: po
  integer :: i, j, k, l

  po           = ubound(x,1)
  n_element(1) = ubound(x,3)
  n_element(2) = ubound(x,4)

  prod = 0

  do l = 1, n_element(2)
  do k = 1, n_element(1)
    do j = 0, po
    do i = 0, po
      prod = prod + x(i,j,k,l) * w(i,j,k,l) * y(i,j,k,l)
    end do
    end do
  end do
  end do

end function WeightedScalarProduct

!===============================================================================
! Tensor products

!-------------------------------------------------------------------------------
!> \brief   Tensor product M_2 x M_1 working with a diagonal matrix M_2
!> \author  Immo Huismann

pure function TensorProduct2D_DxM(matrix_2, matrix_1, val) result(res)
  real(RDP), intent(in) :: matrix_2  (:) !< matrix for the second dimension
  real(RDP), intent(in) :: matrix_1(:,:) !< matrix for the first  dimension
  real(RDP), intent(in) :: val (:,:,:,:) !< values to work on
  real(RDP) :: res(size(val,1),size(val,2),size(val,3),size(val,4))

  integer :: n(2), n_ele(2) ! number of points and number of elements
  integer :: j,k,o,m        ! loop indices

  n     = [size(matrix_1,1), size(matrix_2,1)]
  n_ele = [size(val     ,3), size(val     ,4)]

  do k = 1, n_ele(2)
  do j = 1, n_ele(1)
    ! work on the first dimension utilizing the matrix-matrix product
    res(:,:,j,k) = matmul(matrix_1,val(:,:,j,k))

    ! work on the secong dimension using the diagonal matrix
    do o = 1, n(2)
    do m = 1, n(1)
      res(m,o,j,k) = matrix_2(o) * res(m,o,j,k)
    end do
    end do

  end do
  end do

end function TensorProduct2D_DxM

!-------------------------------------------------------------------------------
!> \brief   Tensor product M_1 x M_2 working with a diagonal matrix M_1
!> \author  Immo Huismann

pure function TensorProduct2D_MxD(matrix_2, matrix_1, val) result(res)
  real(RDP), intent(in) :: matrix_2(:,:)
  real(RDP), intent(in) :: matrix_1  (:)
  real(RDP), intent(in) :: val   (:,:,:,:)

  real(RDP) :: res(size(val,1),size(val,2),size(val,3),size(val,4))

  integer :: n(2), n_ele(2) ! number of points and number of elements
  integer :: j,k,o,m        ! loop indices

  n     = [size(matrix_1,1), size(matrix_2,1)]
  n_ele = [size(val     ,3), size(val     ,4)]

  do k = 1, n_ele(2)
  do j = 1, n_ele(1)

    ! work on the first dimension using the diagonal matrix
    do o = 1, n(2)
    do m = 1, n(1)
      res(m,o,j,k) = matrix_1(m) * val(m,o,j,k)
    end do
    end do

    ! work on the second one with the matrix-matrix product
    res(:,:,j,k) = matmul(res(:,:,j,k), transpose(matrix_2))

  end do
  end do

end function TensorProduct2D_MxD

!-------------------------------------------------------------------------------
!> \brief   Tensor product M_1 x M_2 working with a diagonal matrices M_1, M_2
!> \author  Immo Huismann

pure function TensorProduct2D_DxD(matrix_2, matrix_1, val) result(res)
  real(RDP), intent(in) :: matrix_2(:)
  real(RDP), intent(in) :: matrix_1(:)
  real(RDP), intent(in) :: val   (:,:,:,:)

  real(RDP) :: res(size(val,1),size(val,2),size(val,3),size(val,4))

  integer :: n(2), n_ele(2) ! number of points and number of elements
  integer :: j,k,o,m        ! loop indices

  n     = [size(matrix_1,1), size(matrix_2,1)]
  n_ele = [size(val     ,3), size(val     ,4)]

  do k = 1, n_ele(2)
  do j = 1, n_ele(1)

    ! work on the both dimensions using the diagonal matrices
    do o = 1, n(2)
    do m = 1, n(1)
      res(m,o,j,k) = matrix_1(m) * matrix_2(o) * val(m,o,j,k)
    end do
    end do

  end do
  end do

end function TensorProduct2D_DxD

!===============================================================================

end module Helmholtz_2D
