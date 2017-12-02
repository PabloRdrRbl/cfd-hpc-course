
!     0;95;0c> \file      jacobi_iteration.f
!> \brief     Program solving the Poisson equation
!> \author    Immo Huismann
!> \date      2014/11/27
!> \copyright Institute of Fluid Mechanics
!>            TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================

program Jacobi_Iteration
  use Kind_Parameters, only: RNP
  implicit none

  integer, parameter :: n = 20
  integer, parameter :: iter = 30
  integer :: i, j
  integer :: h    
  real(RNP) :: x_h(n+1, 1), f_h(n+1, 1), u_h(n+1, 1)    
  real(RNP) :: A(n-1, n-1)

  h = 1.0 * 2.0 / n

  x_h(:, 1) = (/((i*h),i=0,n)/)    

  f_h(:, 1) = 1.0 * sin(x_h)    

  u_h(:, 1) = (/(0.0, i=0,n)/)    

  A = reshape((/(0.0, i=1,((n-1)*(n-1)))/), shape(A))

  A(1, 1) = 2.0
  A(n-1, n-1) = 2.0
  A(1+1, 1) = -1.0
  A(n-2, n-1) = -1.0
      
  do i=2,n-2    
    A(i, i) = 2.0
    A(i+1, i) = -1.0
    A(i-1, i) = -1.0
  end do    

  do i=1,iter  
      u_h(:, 1) = (1.0/2.0) * (-f_h(:, 1) - matmul(A, u_h))
  end do
      
end program Jacobi_Iteration
