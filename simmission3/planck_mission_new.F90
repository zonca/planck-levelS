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

module planck_mission_new
  use planck_config
  use general_const
  use general_error
  use general_vector
  use general_time
  use solarsystem_star
  use planck_satellite
  use planck_pointing
  use planck_scanning_new
  use planck_idealdynamics
  use planck_analyticaldynamics
  use planck_pointingperiod
  use planck_io
  use levels_output
  use ls_paramfile_io
  use ephemerides_f
  implicit none
  private

  public :: plmission, pl_mission_init, pl_mission_finish, pl_simmission

  ! Model of the Planck satellite/mission including a physical model of
  ! the satellite, parameters for its orbit, the scan strategy, information
  ! on pointing errors and mission times.
  type plmission
    type(plsatellite) :: satellite
    type(plpointing) :: pointing
    type(plscanstrategy) :: scanstrategy
    type(eph_handle) :: eph_sun
  end type plmission

  ! Length conversions between astronomical units (AU) and metres (M).
  real(dp), parameter :: GNM_AU = 6.6845862e-12_dp

contains

  ! Initialise the overall Planck simulation.
  function pl_mission_init(params) result(mission)
    type(paramfile_handle), intent(inout) :: params
    type(plmission) :: mission

    character (len=filenamelen) :: eph_object

    write(*, '(/,a,/)') 'Planck mission parameters.'

    ! First the (nominal) scan strategy.
    mission%scanstrategy = pl_scanstrategy_init(params)

    ! Then initialise the physical properties of the satellite itself
    mission%satellite = pl_satellite_init(params)

    ! Then the satellite's pointing characteristics.
    mission%pointing = pl_pointing_init(params)

    ! Then the ephemerides
    eph_object = parse_string(params,"ephemerides")
    call eph_open (mission%eph_sun, eph_object)
    call eph_load_body (mission%eph_sun, "Sun")

    write(*, '(/,a,/)') 'Planck mission initialised.'
  end function pl_mission_init

  subroutine pl_mission_finish(mission)
    type(plmission), intent(inout) :: mission

    call eph_close (mission%eph_sun)
    call pl_scanstrategy_finish(mission%scanstrategy)
  end subroutine pl_mission_finish

  ! Simulate the Planck mission. This function simply loops over the
  ! (possibly arbitrary) pointing periods of the mission, simulating
  ! each of these independently. Moreover, due to the expected size
  ! of the simulation, the only output of this routine is i) written
  ! to disk rather than returned to the calling program and ii) a
  ! meta-file containing only header information detailing the files
  ! produced with information on each pointing period.
  subroutine pl_simmission(mission, missionoutput)
    type(plmission), intent(in out) :: mission
    type(ploutput), intent(in) :: missionoutput

    integer :: pt, n_pt

    n_pt = mission%scanstrategy%n_pt
    do pt = 1, n_pt
      call pl_simpointingperiod_scan(pt, mission, missionoutput)
    end do
  end subroutine pl_simmission

  ! Simulate a pointing period starting from the scan strategy and
  ! overall mission parameters.
  subroutine pl_simpointingperiod_scan(pt, mission, missionoutput)
    integer, intent(in) :: pt
    type(plmission), intent(in out) :: mission
    type(ploutput), intent(in) :: missionoutput

    real(dp) :: rotm_nom_startpt(3, 3), period_pt, &
      alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt
    type(gnsec) :: t_startpt
    integer :: n_pt

    n_pt = mission%scanstrategy%n_pt

    call gn_assert(pt>0,'pl_simpointingperiod_scan: pt <= 0',pt)
    call gn_assert(pt<=n_pt,'pl_simpointingperiod_scan: pt > n_pt',pt)

    ! Work out times for this pointing.
    t_startpt = pl_t_startpt(pt, mission%scanstrategy)
    period_pt = pl_period_pt(pt, mission%scanstrategy)

    ! Get the nominal orientation from the scan strategy.
    rotm_nom_startpt = pl_rotm_nom(pt, mission%scanstrategy)
    call gn_rotm2ang(rotm_nom_startpt, &
      alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt)

    call pl_simpointingperiod_nom(pt, t_startpt, period_pt, &
      alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt, &
      mission, missionoutput)
  end subroutine pl_simpointingperiod_scan

  ! Almost all of the functionality of the mission simulation is contained
  ! in this massive, inelegant subroutine.  To some degree it is divided into
  ! initialisations, simulations, output and clean-up, but unfortunately
  ! I can't think of a useful way to compartmentalise the routine further
  ! due to the number of locally initialised variables in use.  An
  ! understanding of how it works can only be obtained by reading through
  ! the next 500 or so lines of comments and code -- sorry.
  !
  ! There is also some ambiguity in the use of both the explicit variables
  ! giving the nominal pointing time and initial Euler angles, as these
  ! are also contained in the mission%scanstrategy.
  subroutine pl_simpointingperiod_nom(pt, t_startpt, period_pt, &
    alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt, &
    mission, missionoutput)
    integer, intent(in) :: pt
    type(gnsec), intent(in) :: t_startpt
    real(dp), intent(in) :: period_pt
    real(dp), intent(in) :: alpha_nom_startpt
    real(dp), intent(in) :: beta_nom_startpt
    real(dp), intent(in) :: gamma_nom_startpt
    type(plmission), intent(in out) :: mission
    type(ploutput), intent(in) :: missionoutput

    integer :: n_pt, n_sm_pt, sm
    type(gnsec) :: t_endpt
    real(dp) :: t_midsm
    type(gnsphangle_double) :: axes(3), axes_startpt(3), zaxis_nom_startpt
    real(dp) :: satpos_startpt(3), &
      satpos_endpt(3), rotm_nom_startpt(3, 3), &
      tb_startpt(3), delta_rate_z_rot_startpt, delta_phase_rot_startpt, &
      rotm(3, 3), rotm_startpt(3, 3), &
      rate_rot_startpt(3), rate_z_rot_startpt, rate_z_rot_nom_startpt, &
      rotm_delta_phase_rot(3, 3), alpha_startpt, beta_startpt, gamma_startpt
    real(dp), pointer :: sataxes(:, :) => null()
    type(plnutation) :: nutation

real(dp) :: t_tmp
type(ssstar):: mystar

    n_pt = mission%scanstrategy%n_pt

    n_sm_pt = -1 ! to shut up compiler warnings

    call gn_assert(pt>0,'pl_simpointingperiod_nom: pt <= 0',pt)
    call gn_assert(pt<=n_pt,'pl_simpointingperiod_nom: pt > n_pt',pt)

    ! INITIALISATIONS.

    ! REFERENCE TIMES.

    ! Reference times for the pointing period: the pointing period,
    ! period_pt; the starting time of the pointing period, t_startpt;
    ! the middle point of the pointing period, t_midpt; and the end
    ! time of the pointing period, t_endpt.

    call gn_assert(period_pt>=0.0, &
      'pl_simpointingperiod_nom: period_pt < 0.0', period_pt)
    t_endpt = t_startpt + period_pt

    if (missionoutput%simulatepointings) then

      ! Reference sampling period. Note that there is a potential difficulty
      ! calculating the number of samples per pointing period as there is
      ! no guarantee that an integer number of samples fit into a pointing
      ! period -- this doesn't matter for the simulation of a single
      ! pointing period, but does matter for the synthesis of an entire
      ! mission. One solution is to doctor the pointing period to be divible
      ! exactly by the sample period; a more realistic solution is to allow
      ! for the repositioning time, the gaps in the data collection this
      ! causes being more than sufficient to get around this problem.

      n_sm_pt = nint(period_pt / missionoutput%period_sm)

      if (abs(real(n_sm_pt,dp) * missionoutput%period_sm - period_pt) > 1e-6_dp) &
        write (*,'(a)') &
          'Warning: pl_simpointingperiod_nom: non-integer number of samples in period'

    end if

    ! SATELLITE POSITION.

    ! Calculate the position of the satellite (in the orbit around the
    ! L2 point of the Sun-Earth system) at the start and the end of the
    ! pointing period.

t_tmp = t_startpt%int + t_startpt%frac + 378691200._dp
satpos_startpt = -eph_posRelSat_m (mission%eph_sun, t_tmp)
t_tmp = t_endpt%int + t_endpt%frac + 378691200._dp
satpos_endpt = -eph_posRelSat_m (mission%eph_sun, t_tmp)

    ! INITIAL CONDITIONS OF THE POINTING PERIOD.

    ! Nominal starting values for the satellite's orientation (expressed
    ! as a rotation matrix to take the satellite from its reference
    ! orientation) and its rotation rate (just the nominal rotation rate
    ! from the scanning strategy).

    rotm_nom_startpt &
      = gn_rotm(alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt)
    zaxis_nom_startpt = gn_sphang(beta_nom_startpt, alpha_nom_startpt)

    rate_z_rot_nom_startpt = mission%scanstrategy%rate_rot

    ! Errors at the start of the pointing period: the Tait-Bryan angles,
    ! tb_startpt, give the initial mispointing, the perturbation to the
    ! rotation rate is given by delta_rate_z_rot_startpt, and the
    ! perturbation to the initial rotation phase, delta_phase_startpt,
    ! is also calculated here. These are the only aspects of the pointing
    ! that are non-deterministic; they are combined here to give the
    ! rotation matrix of the satellite at the start of the pointing
    ! period.

    tb_startpt = pl_monte_tb(mission%pointing)

    delta_phase_rot_startpt = pl_monte_delta_phase_rot(mission%pointing)
    rotm_delta_phase_rot = gn_rotm(delta_phase_rot_startpt, 3)

    rotm_startpt = matmul(matmul(rotm_nom_startpt, rotm_delta_phase_rot), &
      gn_rotm(tb_startpt))
    axes_startpt = gn_rotm2axes(rotm_startpt)
    call gn_rotm2ang(rotm_startpt, &
      alpha_startpt, beta_startpt, gamma_startpt)

    delta_rate_z_rot_startpt = pl_monte_delta_rate_rot(mission%pointing)
    rate_z_rot_startpt = rate_z_rot_nom_startpt + delta_rate_z_rot_startpt

    ! Initialise the pointing period, primarily by calculating various
    ! quantitites related to the inertia tensor of the satellite and
    ! the initial inertial rates (having taken the initial z-axis rotation
    ! rate as an input). Unlike the previous routine, this routine is
    ! deterministic. The output values are the updated satellite structure
    ! and the three inertial rates, along with the nutation structure,
    ! although if we are doing an ideal dynamics simulation then
    ! this is not supplied and that part of the code is avoided.

!FIXME: quickly create the Sun out of nothing
mystar=ss_star_init(1._dp,1._dp,4.5d-6 * (fourpi /  (GNM_AU**2)))

    if (mission%satellite%dynamics == 'ideal') then
      call pl_pointingperiod_init(satpos_startpt, rotm_nom_startpt, &
        tb_startpt, rate_z_rot_startpt, mystar, mission%satellite, &
        rate_rot_startpt)
    else
      call pl_pointingperiod_init(satpos_startpt, rotm_nom_startpt, &
        tb_startpt, rate_z_rot_startpt, mystar, mission%satellite, &
        rate_rot_startpt, nutation)
    end if

    ! POINTING PERIOD SIMULATIONS.

    ! These entire loops (over the samples of a given pointing period)
    ! are only entered if simulatepointings = .true., i.e., that
    ! we actually want to generate pointing period files with outputs
    ! for each sample.

    if (missionoutput%simulatepointings) then

      ! Allocate memory for the angle arrays, either six angles for the
      ! (theta, phi) for each satellite axis or three angles for the three
      ! Euler angles (alpha, beta, gamma).

      allocate(sataxes(n_sm_pt, 6))

      ! SAMPLE LOOPS.

      ! The rest of the pointing period simulation varies greatly depending
      ! on the assumed behaviour of the satellite over the period in
      ! question. The three main options are for an ideal evolution, a
      ! analytical approximation to the dynamics, or a more accurate but
      ! slower numerical integration. At this point all the initial
      ! quantities are ``fundamental'' and physically-meainingful, and
      ! not specific to either method of integration.

      if (mission%satellite%dynamics == 'ideal') then

        ! In the ideal case simply use the scan strategy information to
        ! repeatedly evaluate the orientation of the satellite at each
        ! sample time over the pointing period.

        do sm = 1, n_sm_pt
          ! Reference time for the sampling period.
          t_midsm = (sm-1+0.5_dp) * missionoutput%period_sm

          ! Get the rotation matrix for the sample time and convert to
          ! Euler angles.
          rotm = pl_rotm_ideal(t_midsm, rotm_startpt, rate_rot_startpt(3))

          ! Convert from the rotation matrix to the more compact notion
          ! of angles (assuming they are required for output).
          axes = gn_rotm2axes(rotm)
          sataxes(sm, 1) = axes(1)%theta
          sataxes(sm, 2) = axes(1)%phi
          sataxes(sm, 3) = axes(2)%theta
          sataxes(sm, 4) = axes(2)%phi
          sataxes(sm, 5) = axes(3)%theta
          sataxes(sm, 6) = axes(3)%phi
        end do

      else if (mission%satellite%dynamics == 'analytical') then

        ! Only ring-based scanning strategies can be evolved dynamically

        ! In the analytical case get the orientation of the satellite
        ! at the start of the pointing period and then evolve it using
        ! a simple analytical/approximate integration of the Euler
        ! equation.

        do sm = 1, n_sm_pt
          ! Reference time for the sampling period.
          t_midsm = (sm-1+0.5_dp) * missionoutput%period_sm

          ! Get the rotation matrix for the sample time (given its
          ! value at the reference time at the start of the pointing
          ! period) and convert to Euler angles.
          rotm = pl_rotm_analytical(t_midsm, rotm_startpt, &
            rate_rot_startpt(3), nutation)

          axes = gn_rotm2axes(rotm)
          sataxes(sm, 1) = axes(1)%theta
          sataxes(sm, 2) = axes(1)%phi
          sataxes(sm, 3) = axes(2)%theta
          sataxes(sm, 4) = axes(2)%phi
          sataxes(sm, 5) = axes(3)%theta
          sataxes(sm, 6) = axes(3)%phi
        end do

      else

        call gn_fatal('pl_simpointingperiod_nom: method_dynamics unknown', &
          mission%satellite%dynamics)

      end if

      ! This is the end of the ``if (simulatepointings)'' section of code.

    end if

    ! OUTPUTS

    ! POINTING PERIOD SUMMARY

    ! Write out a copy of the general geometrical information to the
    ! main mission file (only output that is always done).
    call add_summary (pt, t_startpt, axes_startpt, rate_z_rot_startpt, t_endpt)
    if (missionoutput%simulatepointings) call add_pointing(pt, sataxes)

    ! CLEAR-UP

    ! Deallocate temporary memory.
    if (missionoutput%simulatepointings) deallocate(sataxes)

  end subroutine pl_simpointingperiod_nom

end module planck_mission_new
