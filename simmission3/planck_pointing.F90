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

module planck_pointing
  use planck_config
  use general_const
  use general_error
  use general_rand
  use ls_paramfile_io
  implicit none
  private

  public :: plpointing, pl_pointing_init, pl_monte_tb, &
    pl_monte_delta_rate_rot, pl_monte_delta_phase_rot

  ! Planck pointing details, mainly whether or not the pointing can be
  ! accomplished ideally, and, in the case that pointing errors are
  ! included, the magnitude of these errors. Three separate errors are
  ! modelled in the pointing process: the error in the z-axis direction,
  ! with the Monte Carlo mode, mean error and maximum error as specified;
  ! the mode of the phase error at the start of the pointing period; and
  ! the error in the rotation rate, specified by the mode of the error,
  ! the mean error and the maximum error. The latter are given in units
  ! of rotations per second, i.e., s^(-1).
  type plpointing
    character(len=filenamelen) :: mode_zaxis
    real(dp) :: sigma_zaxis_x
    real(dp) :: sigma_zaxis_y
    real(dp) :: delta_zaxis_max_x
    real(dp) :: delta_zaxis_max_y
    character(len=filenamelen) :: mode_phase
    character(len=filenamelen) :: mode_rate_rot
    real(dp) :: sigma_rate_rot
    real(dp) :: delta_rate_rot_max
  end type plpointing

contains

  ! Initialise the pointing properties of the Planck satellite.
  function pl_pointing_init(params) result(pointing)
    type(paramfile_handle), intent(inout) :: params
    type(plpointing) :: pointing

    write(*, '(/,a,/)') 'Planck pointing parameters.'

    ! Z-AXIS POINTING ERRORS
    pointing%mode_zaxis = parse_string(params,'mode_zaxis_pointing')

    ! The mean and maximum pointing error
    ! (dependent on which pointing mode is used).
    select case (pointing%mode_zaxis)
    case ('ideal')
      pointing%sigma_zaxis_x = 0.0_dp
      pointing%sigma_zaxis_y = 0.0_dp
      pointing%delta_zaxis_max_x = 0.0_dp
      pointing%delta_zaxis_max_y = 0.0_dp
    case ('gaussian')
      pointing%sigma_zaxis_x &
        = GNARCMIN_RAD * parse_double(params,'sigma_zaxis_x_pointing', &
        vmin = 0.0_dp)
      pointing%sigma_zaxis_y &
        = GNARCMIN_RAD * parse_double(params,'sigma_zaxis_y_pointing', &
        vmin = 0.0_dp)
      ! Check to see that the old version of this parameter has not been used.
      call gn_assert (.not. key_present (params,'sigma_zaxis_pointing'), &
        'pl_pointing_init: use of sigma_zaxis_pointing is deprecated')
      pointing%delta_zaxis_max_x &
        = GNARCMIN_RAD * parse_double(params,'delta_zaxis_max_x_pointing', &
        vmin = 0.0_dp)
      pointing%delta_zaxis_max_y &
        = GNARCMIN_RAD * parse_double(params,'delta_zaxis_max_y_pointing', &
        vmin = 0.0_dp)
    case default
      call gn_fatal('pl_pointing_init: mode_zaxis unknown', &
        pointing%mode_zaxis)
    end select

    ! Mode of setting the initial scan phase of a pointing period.
    pointing%mode_phase = parse_string(params,'mode_phase_pointing')
    call gn_assert ((pointing%mode_phase=='ideal') .or. &
      (pointing%mode_phase=='random'), &
      'pl_pointing_init: mode_phase unknown', pointing%mode_phase)

    ! Rotation rate errors.
    pointing%mode_rate_rot = parse_string(params,'mode_rate_rot_pointing')

    ! Mean and maximum rotation rate error.
    select case (pointing%mode_rate_rot)
    case ('ideal')
      pointing%sigma_rate_rot = 0.0_dp
      pointing%delta_rate_rot_max = 0.0_dp
    case ('gaussian')
      pointing%sigma_rate_rot = GNDEG_RAD * &
        parse_double(params,'sigma_rate_rot_pointing', vmin = 0.0_dp)
      pointing%delta_rate_rot_max &
        = GNDEG_RAD * parse_double(params,'delta_rate_rot_max_pointing', &
        vmin = 0.0_dp)
    case default
      call gn_fatal('pl_pointing_init: mode_rate_rot unknown', &
        pointing%mode_rate_rot)
    end select
  end function pl_pointing_init

  ! Return random values for the Tait-Bryan angles [with the z-angle,
  ! tb(3), set to zero in all cases], using the error values given in
  ! the pointing structure to determine the variance.
  function pl_monte_tb(pointing) result(tb)
    type(plpointing), intent(in) :: pointing
    real(dp) :: tb(3)

    real(dp) :: delta_zaxis_x, delta_zaxis_y

    select case (pointing%mode_zaxis)
    case ('ideal')
      ! No errors, just set the Tait-Bryan angles to zero.
      tb = 0.0_dp
    case ('gaussian')
      ! Gaussian errors on the z-axis position (but not on the rotation phase).
      delta_zaxis_x = gn_monte_gauss_clip( &
        pointing%sigma_zaxis_x, pointing%delta_zaxis_max_x)
      delta_zaxis_y = gn_monte_gauss_clip( &
        pointing%sigma_zaxis_y, pointing%delta_zaxis_max_y)

      tb(1) = delta_zaxis_x
      tb(2) = delta_zaxis_y
      tb(3) = 0.0_dp
    end select
  end function pl_monte_tb

  ! Return a random value for the error in the scanning rotation rate,
  ! with variance determined by the pointing structure.
  function pl_monte_delta_rate_rot(pointing) result(delta_rate_rot)
    type(plpointing), intent(in) :: pointing
    real(dp) :: delta_rate_rot

    delta_rate_rot = 0.0_dp
    select case (pointing%mode_rate_rot)
    case ('ideal')
      delta_rate_rot = 0.0_dp
    case ('gaussian')
      delta_rate_rot = gn_monte_gauss_clip(pointing%sigma_rate_rot, &
        pointing%delta_rate_rot_max)
    end select
  end function pl_monte_delta_rate_rot

  ! Return a random starting scan phase_rot, phase, generally in the range
  ! 0.0 to 2.0 * pi, although possibly predetermined to coincide with
  ! the nominal scan simulations.
  function pl_monte_delta_phase_rot(pointing) result(delta_phase_rot)
    type(plpointing), intent(in) :: pointing
    real(dp) :: delta_phase_rot

    delta_phase_rot = 0.0_dp
    select case (pointing%mode_phase)
    case ('ideal')
      delta_phase_rot = 0.0_dp
    case ('random')
      delta_phase_rot = gn_monte_rand(0.0_dp, twopi)
    end select
  end function pl_monte_delta_phase_rot

end module planck_pointing
