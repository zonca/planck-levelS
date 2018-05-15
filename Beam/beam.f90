!-----------------------------------------------------------------------------
!
!  Copyright (C) 2002-2013 Mark Ashdown
!
!  This file is part of the "Beam" component of the Planck simulation
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

module beam

  ! This is the master module for the beam package. Using this module
  ! uses all of the other modules in the package.
  !
  ! Mark Ashdown, CPAC.

  use beam_grid
  use beam_cut
  use beam_polar
  use beam_square
  use beam_bessel
  use beam_alm
  use beam_convert
  use beam_interpolate
  use beam_crosspol
  use beam_transform

end module beam
