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

module general_maths
  use planck_config
  use general_const
  implicit none
  private

  public:: gn_atan

contains

  ! Returns the angle (in the range 0.0 to 2 pi) in the direction of the
  ! point (x, y), or 0.0 if x = y = 0.0 (in which case an error message
  ! is also generated).
  function gn_atan(x, y) result(phi)
    real(dp), intent(in) :: x, y
    real(dp) :: phi

    if ((x == 0.0) .and. (y == 0.0)) then
      phi = 0.0
    else
      phi = atan2(y, x)
      if (phi < 0.0) then
        phi = phi + twopi
      end if
    end if
  end function gn_atan

end module general_maths
