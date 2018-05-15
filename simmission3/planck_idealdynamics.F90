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

module planck_idealdynamics
  use planck_config
  use general_vector
  implicit none
  private

  public :: pl_rotm_ideal

contains

  ! Given the the rotation matrix, rotm_0, that takes the satellite from its
  ! reference orientation to that at reference time t_0, return the
  ! rotation matrix, rotm, appropriate at time t_0 + delta_t, given that the
  ! satellite is spinning with about its own z-axis with rate rate_rot_z.
  function pl_rotm_ideal(delta_t, rotm_0, rate_rot_z) result(rotm)
    real(dp), intent(in) :: delta_t
    real(dp), intent(in) :: rotm_0(3, 3)
    real(dp), intent(in) :: rate_rot_z
    real(dp) :: rotm(3, 3)

    real(dp) :: rotm_rot(3, 3), phase_rot

    ! The rotation of the satellite about its own z-axis.
    phase_rot = rate_rot_z * delta_t
    rotm_rot = gn_rotm(phase_rot, 3)

    ! Then combine this with the reference rotation matrix, remembering
    ! that it is the satellite's revolution that must be applied first.
    rotm = matmul(rotm_0, rotm_rot)
  end function pl_rotm_ideal

end module planck_idealdynamics
