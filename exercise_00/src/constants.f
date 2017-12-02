!> \file       constants.f
!> \brief      Definition of common constants
!> \author     Joerg Stiller
!> \date       2013/10/17
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Constants
  use Kind_Parameters
  implicit none
  private

  !-----------------------------------------------------------------------------
  ! numeric constants

  real(RNP), parameter, public ::  ZERO  = 0
  real(RNP), parameter, public ::  ONE   = 1
  real(RNP), parameter, public ::  TWO   = 2
  real(RNP), parameter, public ::  THREE = 3
  real(RNP), parameter, public ::  FOUR  = 4
  real(RNP), parameter, public ::  HALF  = 0.5_RNP
  real(RNP), parameter, public ::  THIRD = ONE/3

  real(RHP), parameter, public ::  PI = 3.1415926535897932384626433832795029_RHP

  !==============================================================================

end module Constants