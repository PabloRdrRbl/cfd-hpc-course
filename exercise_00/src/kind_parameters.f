!> \file       kind_parameters.f
!> \brief      Definition of intrinsic type kind parameters
!> \author     Joerg Stiller
!> \date       2014/03/12
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \todo
!> Once available with the Intel compiler, define RHP via
!>      integer, parameter :: RHP = merge(REAL128, RNP, REAL128 > 0)
!===============================================================================

module Kind_Parameters
  implicit none
  public

  !-----------------------------------------------------------------------------
  ! integer kinds
  integer, parameter :: IXS = selected_int_kind( 2) !<  8 Bit, range > 10^2
  integer, parameter :: IXL = selected_int_kind(18) !< 64 Bit, range > 10^18

  !-----------------------------------------------------------------------------
  ! real kinds
  integer, parameter :: RDP = kind(1D0) !< real double precision
  integer, parameter :: RNP = RDP       !< real normal precision
  integer, parameter :: RHP = RDP       !< real high precision

  !=============================================================================

end module Kind_Parameters