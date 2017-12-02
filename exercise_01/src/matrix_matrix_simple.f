!> \file      matrix_matrix_simple.f
!> \brief     Simple implementation of a matrix-matrix product
!> \author    Immo Huismann
!> \date      2017/11/07
!> \copyright Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================

module Matrix_Matrix_Simple
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: MatrixMatrixSimple

contains

subroutine MatrixMatrixSimple(n, A, B, C)
  integer,   intent(in)  :: n
  real(RDP), intent(out) :: A(n,n)
  real(RDP), intent(in)  :: B(n,n)
  real(RDP), intent(in)  :: C(n,n)

  integer   :: i,j,k
  
  A = 0
  
  do j=1,n  ! Fortran stores colums
      do i=1,n
         !A(i,j) = 0  ! Just in case A = 0 does not work
         do k=1,n      
            A(i,j) = A(i,j) +  B(i,k) * C(k,j)
         end do
      end do
  end do


end subroutine MatrixMatrixSimple

!===============================================================================

end module Matrix_Matrix_Simple
