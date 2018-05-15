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

module solarsystem_orbit
  use planck_config
  use general_error
  use general_time
  use general_vector
  implicit none
  private

  public :: ssorbit, ss_orbit_init, ss_plpos_orbit, ss_plvel_orbit

  ! Structure containing all the information required to specify the
  ! classical Keplerian elliptical orbit for a single planet at a single
  ! time. The rate of change of the planetary orbit is assumed to be
  ! sufficiently slow that the information here is a good approximation
  ! over a range of times, but t_0, included in the structure, gives the
  ! time at which it's a formal best fit to the planet's orbit. The
  ! six independent Keplerian orbital parameters are taken to be: the
  ! orbit's semi-major axis, a (in m); the eccentricity, e; the orbital
  ! inclination, inc (in radians); the longitude of the ascending node,
  ! ascnode (in radians); the longitude of perihelion, lonperi (in radians),
  ! and the mean longitude, meanlon (in radians). The timescale of the
  ! orbit is characterised by mu = G (M_star + M_planet) [in m^3 s^(-2)].
  !
  ! Also stored are several derivable quantities (which should thus never
  ! be calculated externally, rather being calculated ``automatically''
  ! when initialising a planetary orbit). These are: the mean motion,
  ! n = (mu / a^3)^(1/2) [in s^(-1)] and the Laplace vectors, phat and qhat,
  ! which give the coordinate directions of the pericentre point, P, and
  ! the orthogonal vector in the orbit plane, Q. (These make an orthogonal
  ! triad with the angular momentum vector of the orbit, c, such that
  ! Q = c x P.
  type ssorbit
    type(gnsec) :: t_0
    real(dp) :: a
    real(dp) :: e
    real(dp) :: inc
    real(dp) :: ascnode
    real(dp) :: lonperi
    real(dp) :: meanlon
    real(dp) :: mu
    real(dp) :: n
    type(gnsec) :: t_peri
    real(dp) :: phat(3)
    real(dp) :: qhat(3)
  end type ssorbit

contains

  ! Initialise a planetary orbit from a full list of orbital parameters.
  function ss_orbit_init(t_0, a, e, inc, ascnode, lonperi, meanlon, mu) &
    result(orbit)
    type(gnsec), intent(in) :: t_0
    real(dp), intent(in) :: a
    real(dp), intent(in) :: e
    real(dp), intent(in) :: inc
    real(dp), intent(in) :: ascnode
    real(dp), intent(in) :: lonperi
    real(dp), intent(in) :: meanlon
    real(dp), intent(in) :: mu
    type(ssorbit) :: orbit

    real(dp) :: meananom, sininc, cosinc, sinascnode, cosascnode, &
      sinargperi, cosargperi, argperi

    ! Reference time at which the orbit is the best fit to the planet;
    ! in the ideal case this is at all times, but note that this value
    ! is only utilised for the calculation of the planetary phase, via
    ! the formula for t_peri; it makes no difference to the overall
    ! orbit followed.
    orbit%t_0 = t_0

    ! Classic Keplerian orbital parameters, each of which is checked to
    ! see that it is in a sensible range.
    call gn_assert(a>0.0, 'ss_orbit_init: a <= 0.0', a)
    orbit%a = a

    call gn_assert(e>=0.0, 'ss_orbit_init: e < 0.0', e)
    call gn_assert(e<=1.0, 'ss_orbit_init: e > 1.0', e)
    orbit%e = e

    call gn_assert(inc>=-halfpi, 'ss_orbit_init: inc < - pi / 2.0', inc)
    call gn_assert(inc<= halfpi, 'ss_orbit_init: inc > pi / 2.0', inc)
    orbit%inc = inc

    if (ascnode < 0.0) then
      call gn_warning('ss_orbit_init: ascnode < 0.0', ascnode)
    else if (ascnode >= twopi) then
      call gn_warning('ss_orbit_init: ascnode >= 2.0 pi', ascnode)
    end if
    orbit%ascnode = mod(ascnode, twopi)

    if (lonperi < 0.0) then
      call gn_warning('ss_orbit_init: lonperi < 0.0', lonperi)
    else if (lonperi > twopi) then
      call gn_warning('ss_orbit_init: lonperi >= 2.0 pi', lonperi)
    end if
    orbit%lonperi = mod(lonperi, twopi)

    if (meanlon < 0.0) then
      call gn_warning('ss_orbit_init: meanlon < 0.0', meanlon)
    else if (meanlon > twopi) then
      call gn_warning('ss_orbit_init: meanlon >= 2.0 pi', meanlon)
    end if
    orbit%meanlon = mod(meanlon, twopi)

    call gn_assert(mu>0.0,'ss_orbit_init: mu <= 0.0', mu)
    orbit%mu = mu

    ! Calculate other basic/alternative orbital elements: the mean motion,
    ! n; the mean anomaly, meananom; and hence the time of pericentre,
    ! t_peri.
    orbit%n = sqrt(orbit%mu / (orbit%a)**3)
    meananom = mod(orbit%meanlon - orbit%lonperi, twopi)
    orbit%t_peri = t_0 - meananom / orbit%n

    ! Calculate the Laplace vectors, P and Q, precomputing trigonometric
    ! functions, from Appendix A of `Modern Astrodynamics', with the
    ! argument of pericentre, argperi, calculated first.
    argperi = orbit%lonperi - orbit%ascnode

    sininc = sin(orbit%inc)
    cosinc = cos(orbit%inc)
    sinascnode = sin(orbit%ascnode)
    cosascnode = cos(orbit%ascnode)
    sinargperi = sin(argperi)
    cosargperi = cos(argperi)

    orbit%phat(1) &
      = cosargperi * cosascnode - sinargperi * sinascnode * cosinc
    orbit%phat(2) &
      = cosargperi * sinascnode + sinargperi * cosascnode * cosinc
    orbit%phat(3) = sinargperi * sininc

    orbit%qhat(1) &
      = - sinargperi * cosascnode - cosargperi * sinascnode * cosinc
    orbit%qhat(2) &
      = - sinargperi * sinascnode + cosargperi * cosascnode * cosinc
    orbit%qhat(3) = cosargperi * sininc
  end function ss_orbit_init

  ! Return the position of a planet (relative to the central star).
  function ss_plpos_orbit(t, orbit) result(plpos)
    type(gnsec), intent(in) :: t
    type(ssorbit), intent(in) :: orbit
    real(dp) :: plpos(3)

    real(dp) :: meananom, accuracy, semilatus, phat(3), qhat(3), &
      a, e, n, mu, eccenanom
    type(gnsec) :: t_peri

    ! Local copies of orbit variables.
    t_peri = orbit%t_peri
    a = orbit%a
    e = orbit%e
    mu = orbit%mu
    n = orbit%n
    phat = orbit%phat
    qhat = orbit%qhat

    ! First calculate the mean anomaly.
    meananom = n * (t - t_peri)

    ! Then solve Kepler's equation for the eccentric anomaly.
    accuracy = 1.0e-6
    eccenanom = ss_eccenanom_kepler(e, meananom, accuracy)

    ! Next compute the semi-latus rectum, p, of the orbit. Note that there are
    ! two equally applicable formulae for this: p = a (1 - e^2) and p = P / mu
    ! where P is the magnitude of the Laplace vector P.
    semilatus = a * (1.0 - e**2)

    ! Finally calculate the position from Eq. (4.15) of `Modern
    ! Astrodynamics', doing all three components at once.
    plpos = a * (cos(eccenanom) - e) * phat &
      + sqrt(a * semilatus) * sin(eccenanom) * qhat
  end function ss_plpos_orbit

  ! Return the velocity of a planet (relative to that of the central star).
  function ss_plvel_orbit(t, orbit) result(plvel)
    type(gnsec), intent(in) :: t
    type(ssorbit), intent(in) :: orbit
    real(dp) :: plvel(3)

    real(dp) :: plpos(3), eccenanom, accuracy, semilatus, phat(3), &
      qhat(3), a, e, n, mu, meananom, pldist
    type(gnsec) :: t_peri

    ! Local copies of orbit variables.
    t_peri = orbit%t_peri
    a = orbit%a
    e = orbit%e
    mu = orbit%mu
    n = orbit%n
    phat = orbit%phat
    qhat = orbit%qhat

    ! First calculate the mean anomaly.
    meananom = n * (t - t_peri)

    ! Then solve Kepler's equation for the eccentric anomaly.
    accuracy = 1.0e-6
    eccenanom = ss_eccenanom_kepler(e, meananom, accuracy)

    ! Next compute the semi-latus rectum, p, of the orbit. Note that there are
    ! two equally applicable formulae for this: p = a (1 - e^2) and p = P / mu
    ! where P is the magnitude of the Laplace vector P.
    semilatus = a * (1.0 - e**2)

    ! Calculate the position from Eq. (4.15) of `Modern Astrodynamics',
    ! doing all three components at once.
    plpos = a * (cos(eccenanom) - e) * phat &
      + sqrt(a * semilatus) * sin(eccenanom) * qhat

    ! Calculate the distance from the star to the planet (just the magnitude
    ! of the position vector).
    pldist = gn_absv(plpos)

    ! Then calculate the velocity of the planet from Eq. (4.16) of `Modern
    ! Astrodynamics', doing all three components at once.
    plvel = - sqrt(mu * a) * sin(eccenanom) / pldist * phat &
      + sqrt(mu * semilatus) * cos(eccenanom) / pldist * qhat
  end function ss_plvel_orbit

  ! Solve Kepler's equation for the eccentric anomaly, eccenanom, given
  ! the eccentricity of the orbit, e, the mean anomaly (implicitly at
  ! the current time). The solution accuracy is specified by the accuracy
  ! parameter [something like 10^(-6) being appropriate]. The algorithm to
  ! do this is based on that in Section 4.5.2 of `Modern Astrodynamics',
  ! the major change being that the mean anomaly is given as an input
  ! parameter, rather than being calculated in the function. Further, the
  ! starting point for the eccentric anomaly is chosen to be the mean
  ! anomaly for simplicity.
  function ss_eccenanom_kepler(e, meananom, accuracy) result(eccenanom)
    real(dp), intent(in) :: e
    real(dp), intent(in) :: meananom
    real(dp), intent(in) :: accuracy
    real(dp) :: eccenanom

    integer :: i, imax
    real(dp) :: f_eccenanom, dfdeccenanom, deccenanom

    call gn_assert (e>=0.0,'ss_eccenanom_kepler: e < 0.0',e)
    call gn_assert (e<=1.0,'ss_eccenanom_kepler: e > 1.0',e)
    call gn_assert (accuracy>0.0,'ss_eccenanom_kepler: accuracy <= 0.0',&
      accuracy)

    ! As a starting point for the iteration, choose the mean anomaly.
    eccenanom = meananom

    ! Then continue adding terms until the desired accuracy has been
    ! achieved.
    imax = 10
    do i = 1, imax

      ! Value and derivative of the equation that is being solved.
      f_eccenanom = eccenanom - e * sin(eccenanom) - meananom
      dfdeccenanom = 1.0 - e * cos(eccenanom)

      ! Difference between the new value and the old, and the new value.
      deccenanom = - f_eccenanom / dfdeccenanom
      eccenanom = eccenanom + deccenanom

      ! If the difference is sufficiently small, assume that the desired
      ! accuracy has been reached; otherwise iterate again. But if we've
      ! reached the imax'th iteration then exit anyway (with a warning
      ! that the iteration hasn't converged).
      if (abs(deccenanom) < accuracy) exit

      if (i == imax) call gn_warning( &
          'ss_eccenanom_kepler: exiting with deccenanom / accuracy', &
          abs(deccenanom)/accuracy)

    end do
  end function ss_eccenanom_kepler

end module solarsystem_orbit
