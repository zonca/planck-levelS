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

module solarsystem_planet
! Basic simulation of a Solar system planet.
  use planck_config
  use general_const
  use general_error
  use general_time
  use general_vector
  use solarsystem_orbit
  use solarsystem_l2
  implicit none
  private

  public :: ssplanet, ss_planet_init_0, ss_planet_obs

  ! Structure containing information on a single planet, combining the
  ! common English name, mass, m (in kg), physical radius, r (in m) and
  ! the planet's orbit (with respect to the central star). Finally, the
  ! generic properties of the star-planet L2 point are stored in the l2
  ! structure; aside from the location of the L2 point, the local potential
  ! can be precalculated as well.
  type ssplanet
    character(len=filenamelen) :: name
    real(dp) :: m
    real(dp) :: r
    type(ssorbit) :: orbit
    type(ssl2) :: l2
  end type ssplanet

  ! Newton's gravitational constant [in m^3 sec^(-2) kg^(-1)].
  real(dp), parameter :: GNG_KGMSEC = 6.67259e-11_dp

contains

  ! Initialise a planet from a full list of basic parameters, most of which
  ! pertain to the planet's orbit.
  function ss_planet_init_0(name, m, r, t_0, a, e, inc, ascnode, lonperi, &
    meanlon, m_star) result(planet)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: m, r
    type(gnsec), intent(in) :: t_0
    real(dp), intent(in) :: a, e, inc, ascnode, lonperi, meanlon, m_star
    type(ssplanet) :: planet

    real(dp) :: mu

    planet%name = name

    ! Physical properties
    call gn_assert (m>0.0, 'ss_planet_init_0: m <= 0.0', m)
    planet%m = m
    call gn_assert (r>0.0, 'ss_planet_init_0: r <= 0.0', r)
    planet%r = r

    ! Orbital properties.
    mu = GNG_KGMSEC * (m + m_star)
    planet%orbit = ss_orbit_init(t_0, a, e, inc, ascnode, lonperi, meanlon, mu)

    ! L2 properties.
    planet%l2 = ss_l2_init(m_star, m)
  end function ss_planet_init_0

  ! Given the time (in s since 1970.0), t, and the observer's position,
  ! (relative to the central star, in ecliptic coordinates and with the
  ! three Cartesian components measured in m), pos_obs, a number of
  ! observable quantities of the input planet are calculated. These are:
  ! the distance from the central star to the planet (in m), dist_starpl;
  ! the distance from the observer to the planet (in m), dist_pl; the
  ! angular position of the planet as viewed from the observer (in ecliptic
  ! coordinates), ang_pl; the angular radius of the planet as seen by the
  ! observer, angradius_pl; and the angle between the central star and the
  ! observer subtended at the planet, theta_starobs. The rationale for the
  ! choice of observables is that these are sufficient to calculate the
  ! thermal emission from the planet as seen by the observer, including
  ! the possibility that the planet's temperature (and hence thermal
  ! emission) is directly determined by its distance from the central star
  ! and also whether it is the night or day side that is being seen.
  subroutine ss_planet_obs(t, pos_obs, planet, dist_starpl, dist_pl, ang_pl, &
    angradius_pl, theta_starobs)
    type(gnsec), intent(in) :: t
    real(dp), intent(in) :: pos_obs(3)
    type(ssplanet), intent(in) :: planet
    real(dp), intent(out) :: dist_starpl, dist_pl
    type(gnsphangle_double), intent(out) :: ang_pl
    real(dp), intent(out) :: angradius_pl, theta_starobs

    real(dp) :: pos_pl(3), relpos_pl(3), r_pl, costheta_starobs

    ! Calculate the position of the planet, relative to the central star, in
    ! ecliptic coordinates, and the position relative to the observer.
    pos_pl = ss_plpos_orbit(t, planet%orbit)
    relpos_pl = pos_pl - pos_obs

    ! Calculate the distance from the planet to the central star.
    dist_starpl = gn_absv(pos_pl)

    ! Convert the relative position into angular coordinates and distance;
    ! the former are output directly (giving the position of the planet on
    ! the ``sky'' of the observer).
    ang_pl = gn_v2ang(relpos_pl)
    dist_pl = gn_absv(relpos_pl)

    r_pl = planet%r
    if (dist_pl < r_pl) then
      call gn_warning('ss_planet_obs: dist_pl < r_pl: observer inside planet')
      angradius_pl = GNDP_MAX
    else
      angradius_pl = asin(r_pl / dist_pl)
    end if

    ! Finally calculate the angle between the central star and the observer
    ! as seen from the planet. (Note that the inner product should be
    ! computed between the negative of the two vectors, but the double
    ! negatives cancel out.)
    if (dist_starpl == 0.0) then
      call gn_warning( &
        'ss_planet_obs: dist_starpl = 0.0: setting theta_starobs = 0.0')
      theta_starobs = 0.0
    else if (dist_pl == 0.0) then
      call gn_warning( &
        'ss_planet_obs: dist_pl = 0.0: setting theta_starobs = 0.0')
      theta_starobs = 0.0
    else
      costheta_starobs = dot_product(pos_pl, relpos_pl) / (dist_starpl * dist_pl)
      theta_starobs = acos(costheta_starobs)
    end if
  end subroutine ss_planet_obs

end module solarsystem_planet
