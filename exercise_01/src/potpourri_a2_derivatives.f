!> \file      potpourri_a2_derivatives.f
!> \brief     Derivatives
!> \author    Immo Huismann
!> \date      2016/11/23
!> \copyright Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================

module Potpourri_A2_Derivatives
  use Kind_Parameters, only: RDP
  implicit none
  private

  public ::  FirstDerivative
  public :: SecondDerivative

contains

pure function FirstDerivative(h,u_left,u_right) result(du)
  real(RDP), intent(in) :: h
  real(RDP), intent(in) :: u_left
  real(RDP), intent(in) :: u_right
  real(RDP)             :: du

  du = (u_right - u_left) / (2 * h)

end function FirstDerivative

pure function SecondDerivative(h,u_left,u_middle,u_right) result(du)
  real(RDP), intent(in) :: h
  real(RDP), intent(in) :: u_left
  real(RDP), intent(in) :: u_middle
  real(RDP), intent(in) :: u_right
  real(RDP)             :: du

  du = (u_right - 2 * u_middle + u_left) / (h**2)

end function SecondDerivative

!===============================================================================

end module Potpourri_A2_Derivatives
