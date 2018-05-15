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

module solarsystem_l2orbit
  use planck_config
  use general_time
  use general_vector
  use solarsystem_orbit
  use solarsystem_planet
  use solarsystem_l2
  implicit none
  private

  public :: ssl2orbit, ss_l2orbit_init, ss_l2orbitpos, ss_l2vel_planet, &
    ss_l2pos_planet

  ! Orbital parameters for motion in a Lissajous pattern around the L2 point
  ! of a star-planet system). Most of the parameters (aa, c1, c2, oxy, oz, lxy)
  ! are derived from the properties of the space around L2, themselves
  ! related to the mass of the two bodies. The initial conditions of the
  ! orbit are given by the starting time, t_0, the position and velocity
  ! (relative to L2, in m) and the orbital phase. These are not all
  ! indepdendent, and the vector a is calculated from them.
  type ssl2orbit
    type(gnsec) :: t_0
    real(dp) :: pos_rel_0(3)
    real(dp) :: vel_rel_0(3)
    real(dp) :: phase_0
    real(dp) :: a(4)
  end type ssl2orbit

contains

  ! Initialise an orbit around the L2 Lagrange point of the star-planet
  ! system given a reference time (in years since zero), position (in metres
  ! from L2) and orbital phase (in radians).
  function ss_l2orbit_init(t_0, pos_rel_0, phase_0, planet) &
    result(l2orbit)
    type(gnsec), intent(in) :: t_0
    real(dp), intent(in) :: pos_rel_0(3)
    real(dp), intent(in) :: phase_0
    type(ssplanet), intent(in) :: planet
    type(ssl2orbit) :: l2orbit

    real(dp) :: vel_rel_0(3), t1, t2, q1, q2, r1, r2, aa(4, 4)

    ! Extract pre-calculated quantities from the planet structure.
    aa = planet%l2%aa

    ! Calculate the initial velocity from the phase information.
    t1 = aa(1, 1) * pos_rel_0(1) + aa(2, 1) * pos_rel_0(2)
    t2 = aa(1, 2) * pos_rel_0(1) + aa(2, 2) * pos_rel_0(2)
    q1 = aa(3, 1)
    q2 = aa(3, 2)
    r1 = aa(4, 1)
    r2 = aa(4, 2)

    vel_rel_0(3) = pos_rel_0(3) * cos(phase_0)
    vel_rel_0(2) = - (t2 - q2 * t1 / q1) / (r2 - q2 * r1 / q1)
    vel_rel_0(1) = (- t1 - r1 * vel_rel_0(2)) / q1

    ! Copy orbital quantities into the structure.
    l2orbit%t_0 = t_0;
    l2orbit%pos_rel_0 = pos_rel_0
    l2orbit%vel_rel_0 = vel_rel_0
    l2orbit%phase_0 = phase_0

    ! Finally compute the a vector, depending on the initial position and
    ! velocity.
    call ss_l2orbit_init_a(l2orbit, planet%l2)
  end function ss_l2orbit_init

  ! Initialise the initial conditions array a(1: 4) in the L2 orbit
  ! structure.
  subroutine ss_l2orbit_init_a(l2orbit, l2)
    type(ssl2orbit), intent(in out) :: l2orbit
    type(ssl2), intent(in) :: l2

    real(dp) :: pos_rel_0(3), vel_rel_0(3), aa(4, 4), a(4)
    integer :: i

    ! Copy over structure variables.
    pos_rel_0 = l2orbit%pos_rel_0
    vel_rel_0 = l2orbit%vel_rel_0
    aa = l2%aa

    ! Calculate the a values, checking to make sure that tiny values are
    ! set to zero.
    do i = 1, 4

      a(i) = pos_rel_0(1) * aa(1, i) + pos_rel_0(2) * aa(2, i) &
        + vel_rel_0(1) * aa(3, i) + vel_rel_0(2) * aa(4, i)

      if (abs(a(i)) < 1.0e-12) then
        a(i) = 0.0
      end if

    end do

    ! Copy the values across to the structure.
    l2orbit%a = a
  end subroutine ss_l2orbit_init_a

  ! Returns the relative position of the body orbiting around L2 as a
  ! function of time. The position returned is relative to L2, with
  ! coordinates measured in metres and oriented so that the positive
  ! x-axis is from the planet to L2 (and beyond), the positive z-axis
  ! is parallel to the z-axis of the ecliptic and the positive y-axis
  ! is defined such that the three constitute a standard right-handed
  ! triad.
  function ss_l2orbitpos_rel(t, l2orbit, l2) result(l2orbitpos_rel)
    type(gnsec), intent(in) :: t
    type(ssl2orbit), intent(in) :: l2orbit
    type(ssl2), intent(in) :: l2
    real(dp) :: l2orbitpos_rel(3)

    real(dp) :: sinxy, cosxy, cosz, e_plus, e_minus, phase, &
      pos_rel_0(3), phase_0, oxy, oz, lxy, a(4), c1, c2
    type(gnsec) :: t_0

    ! Copy out orbital parameters to local variables.
    t_0 = l2orbit%t_0
    phase_0 = l2orbit%phase_0
    pos_rel_0 = l2orbit%pos_rel_0
    a = l2orbit%a

    oxy = l2%oxy
    oz = l2%oz
    lxy = l2%lxy
    c1 = l2%c1
    c2 = l2%c2

    ! Orbital phase (not sure about this).
    phase = twopi * (t - t_0) / GNYR_SEC

    ! Precalculated quantities.
    e_minus = exp(- phase * lxy)
    e_plus = exp(phase * lxy)

    sinxy = sin(phase * oxy)
    cosxy = cos(phase * oxy)

    cosz = cos(phase * oz + phase_0)

    ! Position relative to L2.
    l2orbitpos_rel(1) = e_plus * a(1) + e_minus * a(2) &
      + cosxy * a(3) + sinxy * a(4)
    l2orbitpos_rel(2) = c1 * (e_plus * a(1) - e_minus * a(2)) &
      - c2 * (sinxy * a(3) - cosxy * a(4))
    l2orbitpos_rel(3) = cosz * pos_rel_0(3)
  end function ss_l2orbitpos_rel

  ! Returns the absolute position of the body orbiting around L2 as a
  ! function of time (relative to the central star).
  function ss_l2orbitpos(t, l2orbit, planet) result(l2orbitpos)
    type(gnsec), intent(in) :: t
    type(ssl2orbit), intent(in) :: l2orbit
    type(ssplanet), intent(in) :: planet
    real(dp) :: l2orbitpos(3)

    real(dp) :: l2orbitpos_rel(3), l2pos_hat(3), l2pos(3)

    ! Firstly find the position of the L2 orbit relative to L2 itself
    ! and the position of L2 separately.
    l2orbitpos_rel = ss_l2orbitpos_rel(t, l2orbit, planet%l2)
    l2pos = ss_l2pos_planet(t, planet)
    l2pos_hat = gn_vhat(l2pos)

    ! Add these two vectors, remembering that l2orbitpos_rel is defined
    ! using a different choice of axes, oriented with respect to the
    ! star-planet pair, rather than the absolute system coordinates.
    ! To effect this transformation simply resolve the components of
    ! the relative displacement vector onto the global reference frame
    ! using the l2pos unit vector. It should be possible to do this more
    ! elegantly using vector/matrix methods, but for the moment an
    ! explicit transformation is used.
    l2orbitpos(1) = l2pos(1) &
      + l2pos_hat(1) * l2orbitpos_rel(1) - l2pos_hat(2) * l2orbitpos_rel(2)
    l2orbitpos(2) = l2pos(2) &
      + l2pos_hat(2) * l2orbitpos_rel(1) + l2pos_hat(1) * l2orbitpos_rel(2)
    l2orbitpos(3) = l2pos(3) &
      + l2orbitpos_rel(3)
  end function ss_l2orbitpos

  ! Return the absolute position of the star-planet L2 point.
  function ss_l2pos_planet(t, planet) result(l2pos)
    type(gnsec), intent(in) :: t
    type(ssplanet), intent(in) :: planet
    real(dp) :: l2pos(3)

    real(dp) :: plpos(3)

    ! Position of the planet.
    plpos = ss_plpos_orbit(t, planet%orbit)

    ! L2 position is proportional the planet's position, but with a
    ! precalculated factor just above unity determined by the planet's
    ! mass.
    l2pos = (1.0 + planet%l2%dratio) * plpos
  end function ss_l2pos_planet

  ! Return the absolute velocity of the star-planet L2 point.
  function ss_l2vel_planet(t, planet) result(l2vel)
    type(gnsec), intent(in) :: t
    type(ssplanet), intent(in) :: planet
    real(dp) :: l2vel(3)

    real(dp) :: plvel(3)

    ! Velocity of the planet.
    plvel = ss_plvel_orbit(t, planet%orbit)

    ! L2 velocity is proportional the planet's velocity, but with a
    ! precalculated factor just above unity determined by the planet's
    ! mass.
    l2vel = (1.0 + planet%l2%dratio) * plvel
  end function ss_l2vel_planet

end module solarsystem_l2orbit
