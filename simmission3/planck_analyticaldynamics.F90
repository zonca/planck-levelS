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

module planck_analyticaldynamics
  use planck_config
  use general_const
  use general_error
  use general_vector
  use general_maths
  use planck_idealdynamics
  implicit none
  private

  public :: plnutation, pl_nutation_init, pl_rotm_analytical

  ! Dynamical quantities for a pointing period as used in the analytical
  ! reconstruction of the satellite's dynamics: the rate of change of the
  ! satellite's rotation period, drate_rot_zdt [in s^(-2)]; the nutation
  ! rate of the satellite, rate_nut [in s^(-1)]; and the six coefficients
  ! that determine the initial values and phase of the nutation, A_x, A_y,
  ! D_x, D_y, E_x, E_y, as defined in equation (A25) of van Leeuwen et al.
  ! (2002). Note that the values of the Tait-Bryan angles at the reference
  ! time of the scan are encoded in these coefficients (specifically D_x
  ! and D_y, the coefficients of the cosine terms, with compensating
  ! values of A_x and A_y, the coefficients of the constant terms).
  type plnutation
    real(dp) :: drate_rot_zdt
    real(dp) :: rate_nut
    real(dp) :: a_x, a_y, d_x, d_y, e_x, e_y
  end type plnutation

contains

  ! Initialise the analytical dynamics structure from input parameters.
  function pl_nutation_init(tb_0, rate_rot_0, drate_rot_zdt, a_x_in, &
    a_y_in, f_1, f_2) result(nutation)
    real(dp), intent(in) :: tb_0(3), rate_rot_0(3)
    real(dp), intent(in) :: drate_rot_zdt, a_x_in, a_y_in, f_1, f_2
    type(plnutation) :: nutation

    real(dp) :: a_x, a_y, d_x, d_y, e_x, e_y, w_1, w_2, &
      gamma, rate_nut, amp_nut_1, amp_nut_2, phase_max_1, phase_max_2

! In trying to debug the simulations run with the anomalous fifteenth
! pointing period the difference between pointing period 15 and, say,
! the previous period, has been narrowed down to the first two components
! of either tb_0 or rate_rot_0. Right; we've now established that the
! problem lies in the translation from the Tait-Bryan angles, tb_0,
! possibly when one is much larger than the other, and this results
! in a nutation amplitude that is much larger than the initial offset.
! In the example here, with the anomalous fifteenth pointing period
! the maximum Tait-Bryan angles should be 2 arcmin, but somehow the
! nutation amplitude gets up to 10 arcmin. Doesn't make any sense to me,
! I have to say ... so what's the answer?

    ! Nutation rate.

    gamma = sqrt(f_1 * f_2)
    rate_nut = gamma * rate_rot_0(3)

    ! Integration constants defined implicitly around equation (A15) of
    ! van Leeuwen et al. (2002). However note that the correspondence
    ! is not obvious.

    w_1 = rate_rot_0(1) - a_x_in
    w_2 = gamma / f_1 * (rate_rot_0(2) - a_y_in)

    ! Further intialisations, A_x, A_y, D_x, D_y, E_x and E_y being given
    ! in various forms in equation (A25) of van Leeuwen et al. (2002). Note
    ! the presence of the initial Tait-Bryan angles in D_x and D_y; these
    ! have a cosine dependence in the evolution of the Tait-Bryan angles
    ! during the pointing period and their presence in D_x and D_y ensure
    ! the correct starting values are maintained. Note that the original
    ! code and van Leeuwen et al. (2002) differ (as is so often the case)
    ! on these lines -- it's possible that a_x_in and a_y_in should be
    ! swapped in the next two lines, which would also mean a_y and a_x
    ! should be swapped in the next two lines giving d_x and d_y.

    a_x = a_x_in / rate_rot_0(3)
    a_y = a_y_in / rate_rot_0(3)

    d_x = tb_0(1) - a_y
    d_y = tb_0(2) - a_x

    e_x = w_1 / w_2 * d_x
    e_y = - w_2 / w_1 * d_y

    ! Work out the two nutation amplitudes (mainly for user output, but
    ! also to check that they aren't ridiculously large).

    phase_max_1 = gn_atan(d_x, e_x)
    amp_nut_1 = d_x * cos(phase_max_1) + e_x * sin(phase_max_1)
    phase_max_2 = gn_atan(d_y, e_y)
    amp_nut_2 = d_y * cos(phase_max_2) + e_y * sin(phase_max_2)

    ! There have been problems with unusually large nutation amplitudes,
    ! much larger than would be expected from the initial Tait-Bryan
    ! angles.  So we've installed a warning if the nutation amplitudes
    ! are more than three times the initial Tait-Bryan angles.

    if (amp_nut_1 > 3.0 * sqrt(tb_0(1)**2 + tb_0(2)**2)) then
      write(*,'(a)') 'pl_nutation_init: amp_nut_1 >> initial TB angle'
      write(*, '(a, f10.5, a, f10.5, a)') &
        '  TB_0 = (', GNRAD_ARCSEC * tb_0(1), ' arcsec, ', &
        GNRAD_ARCSEC * tb_0(2), ' arcsec).'
      write(*, '(a, f10.5, a)') &
        '  nutation amplitude = ', GNRAD_ARCSEC * amp_nut_1, ' arcsec.'
    end if

    if (amp_nut_2 > 3.0 * sqrt(tb_0(1)**2 + tb_0(2)**2)) then
      write(*,'(a)') 'pl_nutation_init: amp_nut_2 >> initial TB angle'
      write(*, '(a, f10.5, a, f10.5, a)') &
        '  TB_0 = (', GNRAD_ARCSEC * tb_0(1), ' arcsec, ', &
        GNRAD_ARCSEC * tb_0(2), ' arcsec).'
      write(*, '(a, f10.5, a)') &
        '  nutation amplitude = ', GNRAD_ARCSEC * amp_nut_2, ' arcsec.'
    end if

    ! Copy these quantities across to the structure.
    nutation%drate_rot_zdt = drate_rot_zdt

    nutation%rate_nut = rate_nut

    nutation%a_x = a_x
    nutation%a_y = a_y
    nutation%d_x = d_x
    nutation%d_y = d_y
    nutation%e_x = e_x
    nutation%e_y = e_y
  end function pl_nutation_init

  ! Given the the rotation matrix, rotm_0, that takes the satellite from its
  ! reference orientation to that at reference time t_0, return the
  ! rotation matrix, rotm, appropriate at time t_0 + deltat, given that the
  ! satellite is spinning about its own z-axis with rate rate_rot_z, and has
  ! the nutation properties summarised by the nutation structure. (Thev
  ! nutation is assumed to be sufficiently non-extreme that the rotation
  ! rate and the axis of rotation can be assumed constant.)
  function pl_rotm_analytical(deltat, rotm_0, rate_rot_z, nutation) &
    result(rotm)
    real(dp), intent(in) :: deltat
    real(dp), intent(in) :: rotm_0(3, 3)
    real(dp), intent(in) :: rate_rot_z
    type(plnutation), intent(in) :: nutation
    real(dp) :: rotm(3, 3)

    real(dp) :: tb(3), rotm_tb(3, 3), rotm_ideal(3, 3)

    ! Calculate the nominal rotation matrix for the scan, using the
    ! ideal dynamics function; if there was no nutation then this
    ! would be the full answer.
    rotm_ideal = pl_rotm_ideal(deltat, rotm_0, rate_rot_z)

    ! Calculate the Tait-Bryan angles at the time in question and from these
    ! the perturbed rotation matrix.
    tb = pl_tb_analytical(deltat, nutation)
    rotm_tb = gn_rotm(tb)

    ! Then combine this with the reference rotation matrix, remembering
    ! that it is the perturbation that must be applied first.
    rotm = matmul(rotm_ideal, rotm_tb)
  end function pl_rotm_analytical

  ! Calculate the Tait-Bryan angles taking the satellite from its reference
  ! orientation to that at time deltat in the given pointing period, using the
  ! analytical dynamical model. The overall model for the dynamics is
  ! described in Section A6 of van Leeuwen et al. (2002).
  function pl_tb_analytical(deltat, nutation) result(tb)
    real(dp), intent(in) :: deltat
    type(plnutation), intent(in) :: nutation
    real(dp) :: tb(3)

    real(dp) :: phase_nut, rate_nut, drate_rot_zdt, &
      a_x, a_y, d_x, d_y, e_x, e_y

    ! Nutation phase, rate_nut = gamma * rate_rot_z = gamma * rate_rot(3).
    rate_nut = nutation%rate_nut

    phase_nut = rate_nut * deltat

    ! Calculate the evolving Tait-Bryan angles at the time in question,
    ! using a simplified version of equation (A24) -- specifically
    ! equation (B1) -- from van Leeuwen et al. (2002).
    a_x = nutation%a_x
    a_y = nutation%a_y
    d_x = nutation%d_x
    d_y = nutation%d_x
    e_x = nutation%e_x
    e_y = nutation%e_y

    tb(1) = a_x + d_x * cos(phase_nut) + e_x * sin(phase_nut)
    tb(2) = a_y + d_y * cos(phase_nut) + e_y * sin(phase_nut)

    ! The perturbation to the azimuthal orientation of the satellite
    ! is simply quadratic with time, a constant acceleration provided
    ! by the torque of the Solar radiation pressure. [Remember that
    ! drate_rot_zdt has units of (radians) s^(-2) so that this is
    ! dimensionally correct.]
    drate_rot_zdt = nutation%drate_rot_zdt

    tb(3) = 0.5 * drate_rot_zdt * deltat**2
  end function pl_tb_analytical

end module planck_analyticaldynamics
