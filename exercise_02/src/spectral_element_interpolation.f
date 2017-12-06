!> \file       spectral_element_interpolation.f
!> \brief      2d spectral element interpolation
!> \author     Joerg Stiller
!> \date       2016/10/22
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Element_Interpolation
  use Kind_Parameters, only: RDP ! Imported real kind
  implicit none                  ! No implicit declaration of variables
  private                        ! No module objects visible by default

  public :: SpectralElementInterpolation

contains

!-------------------------------------------------------------------------------
!> Loop-based implementation of the interpolation operator

subroutine SpectralElementInterpolation(A, u, v)
  real(RDP), intent(in)  :: A(1:,1:,0:,0:) !< Interpolation matrix
  real(RDP), intent(in)  :: u(0:,0:,1:,1:) !< values to interpolate
  real(RDP), intent(out) :: v(1:,1:,1:,1:) !< result from interpolation

  ! Internal Variables
  integer   :: po, ne1, ne2, ni
  integer   :: i, j, m, n, p, q
  real(RDP) :: s

  !.............................................................................
  ! Inquire data sizes

  po  = ubound(u,1)
  ne1 = ubound(u,3)
  ne2 = ubound(u,4)
  ni  = ubound(v,1)

  !.............................................................................
  ! Compute interpolation

  do n = 1, ne2
  do m = 1, ne1
    do q = 1, ni
    do p = 1, ni
      s = 0
      do j = 0, po
      do i = 0, po
         s = s + A(p,q,i,j) * u(i,j,m,n)
      end do
      end do
      v(p,q,m,n) = s
    end do
    end do
  end do
  end do

end subroutine SpectralElementInterpolation

!===============================================================================

end module Spectral_Element_Interpolation
