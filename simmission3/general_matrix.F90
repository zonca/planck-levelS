!-----------------------------------------------------------------------------
!
!  Copyright (C) 1999-2013 Daniel Mortlock
!
!  This file is part of the "simmission" component of the Planck simulation
!  package.
!
!  This code is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This code is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this code; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!-----------------------------------------------------------------------------

module general_matrix
  use planck_config
  use general_error
  implicit none
  private

  public :: gn_jacobian_sm

contains

  ! Returns the number of rows of the symmetric matrix stored column by
  ! column in a vector of length nelem. If nelem cannot be expressed
  ! as nelem = (nrow (nrow + 1)) / 2 then a fatal error occurs.
  function gn_nrow_sm(nelem) result(nrow)
    integer, intent(in) :: nelem
    integer :: nrow

    call gn_assert (nelem>0, 'gn_nrow_sm: nelem <= 0', nelem)
    nrow = nint(0.5 * (sqrt(real(8*nelem+1,dp)) - 1.0))
    call gn_assert ((nrow * (nrow + 1)) / 2 == nelem, &
      'gn_nrow_sm: cannot have a symmetric matrix with nelem', nelem)
  end function gn_nrow_sm

  ! Given a symmetric matrix SM (stored as a vector, column by column),
  ! and a square matrix (in 2-dimensional array format), R, this returns
  ! R^T SM R -- the Jacobian or similarity transform.

  ! At present it only works for 3 x 3 matrices, but this will be changed
  ! in the future.
  function gn_jacobian_sm(r, sm) result(smtrans)
    real(dp), intent(in) :: r(:, :), sm(:)
    real(dp) :: smtrans(size(sm))

    real(dp) :: tm(3, 3)
    integer :: i, nrow

    ! First of all check to see that only 3 x 3 matrices are input
    ! and that the 2-dimensional matrix is compatible with this.

    nrow = gn_nrow_sm(size(sm))
    call gn_assert (nrow==3,'gn_jacobian_sm: nrow /= 3 not yet programmed')
    call gn_assert (size(r)==nrow**2, &
      'gn_jacobian_sm: r is incompatible with for sm')

    ! Now compute the transformation in two parts.

    do i = 1, 3
      tm(1, i) = r(1, i) * sm(1) + r(2, i) * sm(2) + r(3, i) * sm(4)
      tm(2, i) = r(1, i) * sm(2) + r(2, i) * sm(3) + r(3, i) * sm(5)
      tm(3, i) = r(1, i) * sm(4) + r(2, i) * sm(5) + r(3, i) * sm(6)
    end do

    smtrans(1) = tm(1, 1) * r(1, 1) + tm(2, 1) * r(2, 1) + tm(3, 1) * r(3, 1)
    smtrans(2) = tm(1, 1) * r(1, 2) + tm(2, 1) * r(2, 2) + tm(3, 1) * r(3, 2)
    smtrans(3) = tm(1, 2) * r(1, 2) + tm(2, 2) * r(2, 2) + tm(3, 2) * r(3, 2)
    smtrans(4) = tm(1, 1) * r(1, 3) + tm(2, 1) * r(2, 3) + tm(3, 1) * r(3, 3)
    smtrans(5) = tm(1, 2) * r(1, 3) + tm(2, 2) * r(2, 3) + tm(3, 2) * r(3, 3)
    smtrans(6) = tm(1, 3) * r(1, 3) + tm(2, 3) * r(2, 3) + tm(3, 3) * r(3, 3)
  end function gn_jacobian_sm

end module general_matrix
