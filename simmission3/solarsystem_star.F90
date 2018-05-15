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

module solarsystem_star
  use planck_config
  use general_const
  use general_error
  use general_vector
  implicit none
  private

  public :: ssstar, ss_star_init, ss_radpressure_star, ss_star_obs

  ! Structure containing information about the central star of a Solar
  ! system (i.e., the Sun, although this is kept general to allow for
  ! other applications). The parameters are mass, m (in kg), radius,
  ! r (in m), the total radiation force (integrated over a sphere centred
  ! on the star), radforce (in N).
  type ssstar
    real(dp) :: m
    real(dp) :: r
    real(dp) :: radforce
  end type ssstar

contains

  ! Initialise a star from the full list of basic parameters: radius,
  ! r (in kg) and radius, r (in m).
  function ss_star_init(m, r, radforce) result(star)
    real(dp), intent(in) :: m, r, radforce
    type(ssstar) :: star

    ! Physical properties.
    call gn_assert(m>0.0,'ss_star_init: m <= 0.0', m)
    star%m = m
    call gn_assert(r>0.0,'ss_star_init: r <= 0.0', r)
    star%r = r
    call gn_assert(radforce>=0.0,'ss_star_init: radforce < 0.0', radforce)
    star%radforce = radforce
  end function ss_star_init

  ! Given the observer's position in the Solar system, calculate the
  ! observables of the central star: its angular position in the sky of
  ! the observer (in ecliptic coordinates), its distance from the observer,
  ! and its angular size as seen by the observer.
  subroutine ss_star_obs(pos_obs, star, ang_star, dist_star, angradius_star)
    real(dp), intent(in) :: pos_obs(3)
    type(ssstar), intent(in) :: star
    type(gnsphangle_double), intent(out) :: ang_star
    real(dp), intent(out) :: dist_star, angradius_star

    real(dp) :: relpos_star(3), r_star

    ! The position of the star relative to the observer is just the negative
    ! of the star-centric position of the observer.
    relpos_star = - pos_obs

    ! Then convert the relative position into angular coordinates and
    ! distance; the former are output directly (giving the position of
    ! the planet on the ``sky'' of the observer); the latter are also output
    ! (to allow for correct flux calculations) as well as being converted
    ! into angular sizes.
    ang_star = gn_v2ang(relpos_star)
    dist_star = gn_absv(relpos_star)

    r_star = star%r
    if (dist_star < r_star) then
      call gn_warning('ss_star_obs: dist_star < r_star: observer inside star')
      angradius_star = GNDP_MAX
    else
      angradius_star = asin(r_star / dist_star)
    end if
  end subroutine ss_star_obs

  ! Stellar radiation pressure [in N m^(-2)] felt at a distance dist (in m)
  ! from the central star, calculated using the normalisation constant of
  ! the radiation pressure an AU from the Sun.
  function ss_radpressure_star(dist, star) result(radpressure)
    real(dp), intent(in) :: dist
    type(ssstar), intent(in) :: star
    real(dp) :: radpressure

    call gn_assert(dist>=0.0,'ss_radpressure_star: dist < 0.0', dist)
    if (dist < star%r) then
      call gn_warning( &
        'ss_radpressure_star: dist < r_star: assuming infinite', dist)
      radpressure = GNDP_MAX
    else
      radpressure = star%radforce / (fourpi * dist**2)
    end if
  end function ss_radpressure_star

end module solarsystem_star
