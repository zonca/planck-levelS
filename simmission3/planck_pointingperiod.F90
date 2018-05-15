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

module planck_pointingperiod
  use planck_config
  use general_const
  use general_vector
  use solarsystem_star
  use planck_satellite
  use planck_analyticaldynamics
  implicit none
  private

  public :: pl_pointingperiod_init

contains

  ! Initialise a pointing period, calculating all the needed quantities
  ! at the start of the period. This routine is completely deterministic --
  ! any errors in the initial pointing of the satellite are implicit in
  ! the input values for the satellite axes, Tait-Bryan angles and the
  ! angular rotation. The only effect of this routine is to calculate
  ! a number of quantities that are used later on; no global or permanent
  ! variables are affected.
  subroutine pl_pointingperiod_init(pos_sat, rotm_nom, tb, rate_rot_z, &
    star, satellite, rate_rot, nutation)
    real(dp), intent(in) :: pos_sat(3)
    real(dp), intent(in) :: rotm_nom(3, 3)
    real(dp), intent(in) :: tb(3)
    real(dp), intent(in) :: rate_rot_z
    type(ssstar), intent(in) :: star
    type(plsatellite), intent(in) :: satellite
    real(dp), intent(out) :: rate_rot(3)
    type(plnutation), optional, intent(out) :: nutation

    real(dp) :: pos_cog(3), rotm_srs2irs(3, 3), i_irs(6), f_star, &
      specref_panel, diffref_panel, area_panel, f_star_max, cosxi_star, &
      xi_star, dist_star, g, h, drate_rot_zdt, f_1, f_2, r_1, r_2, gamma, &
      a_x, a_y, b_x, b_y, b, torque_irs_z

    ! Make local copies of the *input* satellite variables: the centre of
    ! gravity, the inertia tensor and the solar panel properties.
    pos_cog = satellite%pos_cog
    rotm_srs2irs = satellite%rotm_srs2irs
    i_irs = satellite%i_irs

    specref_panel = satellite%solarpanel%specref
    diffref_panel = satellite%solarpanel%diffref
    area_panel = satellite%solarpanel%area

    ! Precalculate the angle between the satellite's z-axis and the central
    ! star, xi_star, as well as the distance from satellite to star. Note
    ! that the satellite's z-axis is just given by the third column of the
    ! rotation matrix giving its overall orientation; hence in the nominal
    ! position it is the third column of rotm_nom. For most scan strategies
    ! xi_star < 0.05 radians.
    cosxi_star = dot_product(gn_vhat(pos_sat), gn_vhat(rotm_nom(:, 3)))
    xi_star = acos(cosxi_star)
    dist_star = gn_absv(pos_sat)

    ! Calculate the force of the Solar radiation pressure on the satellite's
    ! Solar panel (assuming the satellite was face on) and the actual
    ! force of the radiation pressure (given the orientation of the
    ! satellite). The total force on the satellite should be about
    ! f_star = 7 * 10^(-5) N.
    f_star_max = ss_radpressure_star(dist_star, star) * area_panel
    f_star = f_star_max * cosxi_star

    ! Solar panel scattering coefficients [as introduced just before
    ! equation (A8) in van Leeuwen et al. (2002)].
    g = 1.0 - specref_panel
    h = (1.0 + specref_panel) * cosxi_star + 2.0 / 3.0 * diffref_panel

    ! Calculate the rate of change of the satellite rotation at the
    ! start of the pointing period (to be used in the analytic simulations
    ! in which it is assumed to stay constant over the pointing period).
    ! The expression for the torque is the third component of equation (A14)
    ! of van Leeuwen et al. (2002), which then gives the rate of change
    ! of the satellite's spin as per equation (A4). For the nominal
    ! satellite parameters (and Solar parameters) it should be around
    ! 2 * 10^(-6) in the x- and y-directions and 10^(-7) in the z-direction,
    ! although this does depend sensitively on how the satellite is
    ! oriented towards the Sun. This in turn implies that the rate of
    ! change of the rotation rate should be around 10^(-6) arcsec s^(-2).
    torque_irs_z = f_star * (g - h) &
      * (pos_cog(1) * rotm_srs2irs(2, 3) + pos_cog(2) * rotm_srs2irs(3, 1))
    drate_rot_zdt = torque_irs_z / i_irs(6)

    ! Factors f_1 and f_2, defined in equation (A5) of van Leeuwen et al.
    ! (2002); obviously ratios of inertia tensor elements, but not sure
    ! what intuitive significance they have.
    f_1 = (i_irs(6) - i_irs(3)) / i_irs(1)
    f_2 = (i_irs(6) - i_irs(1)) / i_irs(3)

    ! Gamma, as defined at the end of Section A3 of van Leeuwen et al. (2002).
    ! This should have a value of approximately gamma = 0.3 for the Planck
    ! satellite.
    gamma = sqrt(f_1 * f_2)

    ! Factors r_1 and r_2, introduced in equation (A28) of van Leeuwen
    ! et al. (2002).
    r_1 = (1.0 - f_1) / (1.0 - gamma**2) / rate_rot_z
    r_2 = - (1.0 - f_2) / (1.0 - gamma**2) / rate_rot_z

    ! Another initialisation quantity used in the analytical approximation
    ! to the dynamics, these are implicitly defined in equation (A23) of
    ! van Leeuwen et al. (2002), with the second line of each formula
    ! given explicitly in equation (A14).
    a_x = f_star / (i_irs(3) * rate_rot_z * f_2) &
      * (pos_cog(2) * h + pos_cog(3) * g * rotm_srs2irs(3, 2) &
      - pos_cog(1) * h * rotm_srs2irs(3, 1))
    a_y = f_star / (i_irs(1) * rate_rot_z * f_1) &
      * (- pos_cog(1) * h - pos_cog(3) * g * rotm_srs2irs(3, 1) &
      - pos_cog(2) * h * rotm_srs2irs(2, 1))

    ! Calculate the coefficients b_x and b_y as defined in equation (A16)
    ! of van Leeuwen et al. (2002). These are used in the calculation of
    ! the intial rotation rates about the x- and y-axes.
    b = 0.5 * pos_cog(3) * g * sin(2.0 * xi_star) * f_star_max &
      / (i_irs(1) + i_irs(3) - i_irs(6)) / rate_rot_z
    b_y = b * (2.0 * i_irs(1) / i_irs(6) - 1.0)
    b_x = b * (2.0 * i_irs(3) / i_irs(6) - 1.0)

    ! Finally, initialise the x- and y- rotation rates from the initial
    ! values of the Tait-Bryan angles for those axes, tb(1) and tb(2),
    ! a_x and a_y, the nominal rotation rate about the z-axis, rate_rot_z,
    ! and r_1 and r_2.
    rate_rot(1) = (tb(2) - a_x / rate_rot_z) / r_2 - a_x
    rate_rot(2) = (tb(1) - a_y / rate_rot_z) / r_1 + a_y
    rate_rot(3) = rate_rot_z

    ! Finally write the necessary calculated values into the structure
    ! used in doing the analytical dynamical simulation, if these are
    ! needed.
    if (present(nutation)) &
      nutation &
        = pl_nutation_init(tb, rate_rot, drate_rot_zdt, a_x, a_y, f_1, f_2)
  end subroutine pl_pointingperiod_init

end module planck_pointingperiod
