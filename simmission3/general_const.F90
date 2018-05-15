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

module general_const
  use planck_config
  implicit none
  private

  ! Size limits.
  ! FIXME simmission3 only
  real(dp), public, parameter :: GNDP_MAX = 0.99 * huge(0.0_dp)

  real(dp), public, parameter :: &
    GNRAD_DEG = 57.29577951308232087679815481410517033241_dp, &
    GNDEG_RAD = 0.01745329251994329576923690768488612713443_dp, &
    GNARCMIN_RAD = 0.0002908882086657215961539484614147687855738_dp, &
    GNRAD_ARCSEC = 206264.8062470963551564733573307786131967_dp

end module general_const
