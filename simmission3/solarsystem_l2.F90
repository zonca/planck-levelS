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

module solarsystem_l2
  use planck_config
  use general_error
  implicit none
  private

  public :: ssl2, ss_l2_init

  ! Information on an L2 second Lagrange point, including the distance
  ! ratio between L2 to primary and secondary to primary and a variety
  ! of quantities relating to the potential in the vicinity of L2.
  type ssl2
    real(dp) :: dratio
    real(dp) :: aa(4, 4)
    real(dp) :: c1
    real(dp) :: c2
    real(dp) :: oxy
    real(dp) :: oz
    real(dp) :: lxy
  end type ssl2

  ! Maximum mass fraction for accurate calculation of the star-planet L2
  ! position.
  real(dp), parameter :: SSM_RATIO_MAX = 0.01_dp

contains

  ! Perform generic initialisations of the L2 orbit structure that depend
  ! only on the properties of the gravitational potential in the vicinity
  ! of L2. Local copies of all the variables are used for the calculation
  ! before the values are copied into the structure.
  function ss_l2_init(m_star, m_planet) result(l2)
    real(dp), intent(in) :: m_star
    real(dp), intent(in) :: m_planet
    type(ssl2) :: l2

    real(dp) :: mratio, dratio, k, d1, d2, oxy, oz, lxy, c1, c2, &
      aa(4, 4)

    call gn_assert (m_star>0.0,'ss_l2_init: m_star <= 0.0', m_star)
    call gn_assert (m_planet>0.0,'ss_l2_init: m_planet <= 0.0', m_planet)

    ! Calculate the mass ratio of the star and planet.
    mratio = m_planet / m_star

    ! Calculate the ratio of the distance from planet to L2 to that of
    ! the star to L2.
    dratio = ss_l2_dratio(m_star, m_planet)

    k = mratio / (dratio**3) + (1.0_dp - mratio) / (1.0 + dratio)**3

    ! Calculations to do with the orbital period.
    oxy = sqrt((- k + 2.0_dp + sqrt((9.0_dp * k - 8.0_dp) * k)) &
      / 2.0_dp)
    oz = sqrt(k)
    lxy = sqrt((k - 2.0_dp + sqrt((9.0_dp * k - 8.0_dp) * k)) &
      / 2.0_dp)

    c1 = (lxy * lxy - 1.0_dp - 2.0_dp * k) / (2.0_dp * lxy)
    c2 = (oxy * oxy + 1.0_dp + 2.0_dp * k) / (2.0_dp * oxy)

    d1 = c1 * lxy + c2 * oxy
    d2 = c1 * oxy - c2 * lxy

    ! Initialise the aa matrix.
    aa(1, 1) = 0.5_dp * c2 * oxy / d1
    aa(2, 1) = 0.5_dp * oxy / d2
    aa(3, 1) = - 0.5_dp * c2 / d2
    aa(4, 1) = 0.5_dp / d1

    aa(1, 2) = aa(1, 1)
    aa(2, 2) = - aa(2, 1)
    aa(3, 2) = - aa(3, 1)
    aa(4, 2) = aa(4, 1)

    aa(1, 3) = c1 * lxy / d1
    aa(2, 3) = 0.0_dp
    aa(3, 3) = 0.0_dp
    aa(4, 3) = - 1.0_dp / d1

    aa(1, 4) = 0.0_dp
    aa(2, 4) = - lxy / d2
    aa(3, 4) = c1 / d2
    aa(4, 4) = 0.0_dp

    ! Copy the calculated values into the structure.
    l2%dratio = dratio
    l2%aa = aa
    l2%oxy = oxy
    l2%oz = oz
    l2%lxy = lxy
    l2%c1 = c1
    l2%c2 = c2
  end function ss_l2_init

  ! Calculates the L2 point of the star-planet system. The return value,
  ! dratio, is the ratio of the distance between the L2 point and the
  ! planet to the distance between the L2 point and the star.
  function ss_l2_dratio(m_star, m_planet) result(dratio)
    real(dp), intent(in) :: m_star
    real(dp), intent(in) :: m_planet
    real(dp) :: dratio

    real(dp) :: mratio, mratio3

    call gn_assert (m_star>0.0,'ss_l2_dratio: m_star <= 0.0', m_star)
    call gn_assert (m_planet>0.0,'ss_l2_dratio: m_planet <= 0.0', m_planet)

    ! Mass ratio quantities, with warnings if the planet is too heavy
    ! for the assumptions here to be used.
    mratio = m_planet / (m_star + m_planet)
    mratio3 = (mratio/3.0_dp)** (1.0_dp/3.0_dp)

    if (mratio > 1.0) then
      call gn_warning( &
        'ss_l2_dratio: planet heavier than star: mratio > 1.0', mratio)
    else if (mratio > 0.1) then
      call gn_warning('ss_l2_dratio: mratio << 1.0 assumption broken', &
        mratio)
    end if

    ! Calculate dratio, using a power law expansion of the fifth order
    ! equation for the L2 points that is accurate for small planets.
    dratio = mratio3 + mratio3**2 / 3.0 &
      - mratio3**3 / 9.0 + 50.0 * mratio3**4 / 81.0
  end function ss_l2_dratio

end module solarsystem_l2
