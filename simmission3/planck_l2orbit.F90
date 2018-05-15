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

module planck_l2orbit
  use planck_config
  use general_const
  use general_time
  use solarsystem_planet
  use solarsystem_l2orbit
  use ls_paramfile_io
  implicit none
  private

  public :: pl_l2orbit_init

contains

  ! Initialise the satellite orbit about L2 of the Earth-Sun system,
  ! just calling the more general routine in the Solar system package.
  ! It is assumed that the input time is the mission start time, although
  ! it need not be.
  function pl_l2orbit_init(params, tsd, tst, earth) result(l2orbit)
    type(paramfile_handle), intent(inout) :: params
    character(len=*), intent(in) :: tsd,tst
    type(ssplanet), intent(in) :: earth
    type(ssl2orbit) :: l2orbit

    real(dp) :: pos_rel_0(3), phase_0
    type(gnsec) :: t_0
    character(len=filenamelen) :: string_date, string_time

    write(*, '(/,a,/)') 'Planck L2 orbit parameters.'

    ! Set the reference time for the orbit (using the mission start time
    ! as the default).
    string_date = parse_string(params,'date_l2_0_orbit',tsd)
    string_time = parse_string(params,'time_l2_0_orbit',tst)
    t_0 = gn_strings2sec(string_date, 'yyyymmdd', string_time, 'hhmmssdsss')

    ! Get the reference phase and location in the Lissajous orbit.
    pos_rel_0(1) = parse_double(params,'pos_l2_0_x_orbit')
    pos_rel_0(2) = parse_double(params,'pos_l2_0_y_orbit')
    pos_rel_0(3) = parse_double(params,'pos_l2_0_z_orbit')

    phase_0 = GNDEG_RAD * parse_double(params,'phase_l2_0_orbit', &
      vmin=0.0_dp, vmax=360.0_dp)

    ! Initialise the Lissajous orbit parameters using these inputs (and
    ! copying the inputs into the orbit structure).
    l2orbit = ss_l2orbit_init(t_0, pos_rel_0, phase_0, earth)
  end function pl_l2orbit_init

end module planck_l2orbit
