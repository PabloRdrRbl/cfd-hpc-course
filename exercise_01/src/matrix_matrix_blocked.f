!> \file      matrix_matrix_blocked.f
!> \brief     Blocked implementation of a matrix-matrix product
!> \author    Immo Huismann
!> \date      2017/11/ 7
!> \copyright Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================

module Matrix_Matrix_Blocked
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: MatrixMatrixBlocked

contains

subroutine MatrixMatrixBlocked(n, A, B, C)
  integer,   intent(in)  :: n
  real(RDP), intent(out) :: A(n,n)
  real(RDP), intent(in)  :: B(n,n)
  real(RDP), intent(in)  :: C(n,n)

  ! this is definitly incorrect, change this
  A = B + C

end subroutine MatrixMatrixBlocked

!===============================================================================

end module Matrix_Matrix_Blocked
