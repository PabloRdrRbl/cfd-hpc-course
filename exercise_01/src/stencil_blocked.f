!> \file      stencil_blocked.f
!> \brief     A baseline version of the 5 point stencil
!> \author    Immo Huismann
!> \date      2016/11/21
!> \copyright Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Stencil_Blocked
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: StencilBlocked

contains

!-------------------------------------------------------------------------------
!> \brief   A baseline version of the five point stencil
!> \author  Immo Huismann

subroutine StencilBlocked(n,h,a,b)
  integer,   intent(in)  :: n
  real(RDP), intent(in)  :: h
  real(RDP), intent(in)  :: a(n,n)
  real(RDP), intent(out) :: b(n,n)

  integer   :: i,j
  real(RDP) :: s

  s = 1 / h**2

  do j = 2, n - 1
  do i = 2, n - 1
    b(i,j) = (a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1)) * s
  end do
  end do

end subroutine StencilBlocked

!===============================================================================

end module Stencil_Blocked
