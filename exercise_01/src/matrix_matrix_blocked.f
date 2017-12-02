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

    integer :: nn  ! Block size
    integer :: i, j, k, ii, jj, kk
    
    nn = 50
    
    A = 0

    ! Iterate blocks
    do k=1,n,nn
       do j=1,n,nn
          do i=1,n,nn
             do jj=j,min(j+(nn-1), n)
                do ii=i,min(i+(nn-1), n)
                   do kk=k,min(k+(nn-1), n)
                      A(ii, jj) = A(ii,jj) + B(ii, kk) * C(kk, jj)
                   end do
                end do
             end do
          end do
       end do
    end do
    
  end subroutine MatrixMatrixBlocked
        

!===============================================================================

end module Matrix_Matrix_Blocked
