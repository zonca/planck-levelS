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

module planck_mission
  use planck_config
  use general_const
  use general_error
  use general_vector
  use general_time
  use general_matrix
  use solarsystem_l2orbit
  use solarsystem_solsys
  use solarsystem_planet
  use solarsystem_orbit
  use solarsystem_star
  use planck_satellite
  use planck_l2orbit
  use planck_pointing
  use planck_scanning
  use planck_idealdynamics
  use planck_analyticaldynamics
  use planck_pointingperiod
  use planck_io
  use levels_output
  use ls_paramfile_io
  implicit none
  private

  public :: plmission, pl_mission_init, pl_simmission

  ! Model of the Planck satellite/mission including a physical model of
  ! the satellite, parameters for its orbit, the scan strategy, information
  ! on pointing errors and mission times.
  type plmission
    type(gnsec) :: t_start
    type(gnsec) :: t_end
    type(plsatellite) :: satellite
    type(ssl2orbit) :: l2orbit
    type(plpointing) :: pointing
    type(plscanstrategy) :: scanstrategy
  end type plmission

  ! Length conversions between astronomical units (AU), kilometres (KM)
  ! and metres (M).
  real(dp), parameter :: GNM_AU = 6.6845862e-12_dp

contains

  ! Initialise the overall Planck simulation.
  function pl_mission_init(params,solsys) result(mission)
    type(paramfile_handle), intent(inout) :: params
    type(sssolsys), intent(in) :: solsys
    type(plmission) :: mission

    type(ssplanet) :: earth
    character(len=20) :: tsd, tst, ted, tet

    write(*, '(/,a,/)') 'Planck mission parameters.'

    ! Firstly get the mission start and end times.
    tsd = trim(parse_string(params,'date_start'))
    tst = trim(parse_string(params,'time_start'))
    ted = trim(parse_string(params,'date_end'))
    tet = trim(parse_string(params,'time_end'))
    mission%t_start = gn_strings2sec(tsd, 'YYYYMMDD', tst, 'hhmmssdsss')
    mission%t_end = gn_strings2sec(ted, 'YYYYMMDD', tet, 'hhmmssdsss')

    ! Make sure that t_end is after t_start; if not exit with humourous
    ! message to the user.
    call gn_assert(mission%t_end>mission%t_start, &
      'pl_missiontimes_init: t_start >= t_end (detonation on launchpad!)')

    ! Then initialise the physical properties of the satellite itself
    mission%satellite = pl_satellite_init(params)

    ! Next comes the satellite's orbit around the Sun-Earth L2 point.
    earth = ss_planet_name('Earth', solsys)
    mission%l2orbit = pl_l2orbit_init(params,tsd,tst,earth)

    ! Then the satellite's pointing characteristics.
    mission%pointing = pl_pointing_init(params)

    ! Finally the (nominal) scan strategy.
    mission%scanstrategy = pl_scanstrategy_init(params, tsd, tst, &
      mission%t_start, mission%t_end)

    write(*, '(/,a,/)') 'Planck mission initialised.'
  end function pl_mission_init

  ! Trim the mission end time to ensure that an integer number of pointing
  ! periods occur over the entire mission.
  subroutine pl_trimmission(mission)
    type(plmission), intent(in out) :: mission

    integer :: n_pt
    real(dp) :: period_pt
    type(gnsec) :: t_start, t_end

    t_start = mission%t_start
    n_pt = mission%scanstrategy%n_pt
    period_pt = mission%scanstrategy%period_pt

    ! Calculate the end time implied by the pointing period and the number
    ! of pointings.
    t_end = t_start + real(n_pt,dp) * period_pt

    ! Replace the mission end time with this value, but write out a
    ! message to the user if there is a significant difference.
    if (abs(t_end - mission%t_end) > 1.0) then
      write(*, '(a)') 'Trimming mission end time to fit pointing period.'
    end if

    mission%t_end = t_end
  end subroutine pl_trimmission

  ! Simulate the Planck mission. This function simply loops over the
  ! (possibly arbitrary) pointing periods of the mission, simulating
  ! each of these independently. Moreover, due to the expected size
  ! of the simulation, the only output of this routine is i) written
  ! to disk rather than returned to the calling program and ii) a
  ! meta-file containing only header information detailing the files
  ! produced with information on each pointing period.
  subroutine pl_simmission(mission, solsys, missionoutput)
    type(plmission), intent(in out) :: mission
    type(sssolsys), intent(in) :: solsys
    type(ploutput), intent(in) :: missionoutput

    integer :: pt, n_pt

    ! MISSION FILE

    ! Write out overall mission file -- essentially just headers with
    ! information about how the actual mission simulation has been output.

    ! POINTING PERIODS

    ! Then loop over each pointing period in turn.

    n_pt = mission%scanstrategy%n_pt

    ! This loop can usefully be parallelised -- in fact doing so is
    ! likely to make pl_simmission become close to maximally efficient
    ! in any multiple processor situation. The large mission structure
    ! with the satellite data etc. must be copied due to the use of
    ! temporary variables within them However one disadvantage of this
    ! approach is that, in cases where the random number generator
    ! is used for pointing errors, etc., the run becomes unrepeatable,
    ! even in the case that the random number seeds are known as the
    ! order of the calls to the random number generator are mixed up.

! $omp parallel do &
! $omp   schedule(dynamic) &
! $omp   shared(n_pt, solsys, missionoutput) &
! $omp   private(pt, mission)

    do pt = 1, n_pt
      call pl_simpointingperiod_scan(pt, mission, solsys, missionoutput)
    end do

! $omp end parallel do

  end subroutine pl_simmission

  ! Simulate a pointing period starting from the scan strategy and
  ! overall mission parameters.
  subroutine pl_simpointingperiod_scan(pt, mission, solsys, missionoutput)
    integer, intent(in) :: pt
    type(plmission), intent(in out) :: mission
    type(sssolsys), intent(in) :: solsys
    type(ploutput), intent(in) :: missionoutput

    real(dp) :: rotm_nom_startpt(3, 3), period_pt, &
      alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt
    type(gnsec) :: t_startpt
    integer :: n_pt
    type(ssplanet) :: earth

    n_pt = mission%scanstrategy%n_pt

    call gn_assert(pt>0,'pl_simpointingperiod_scan: pt <= 0',pt)
    call gn_assert(pt<=n_pt,'pl_simpointingperiod_scan: pt > n_pt',pt)

    ! Get the Earth's orbital properties -- needed for some scan
    ! strategies.
    earth = ss_planet_name('Earth', solsys)

    ! Work out times for this pointing.
    t_startpt = pl_t_startpt(pt, mission%t_start, mission%scanstrategy)
    period_pt = pl_period_pt(pt, mission%scanstrategy)

    ! Get the nominal orientation from the scan strategy.
    rotm_nom_startpt = pl_rotm_nom( &
      pt, mission%t_start, mission%scanstrategy, earth)
    call gn_rotm2ang(rotm_nom_startpt, &
      alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt)

    call pl_simpointingperiod_nom(pt, t_startpt, period_pt, &
      alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt, &
      mission, solsys, missionoutput)
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
    mission, solsys, missionoutput)
    integer, intent(in) :: pt
    type(gnsec), intent(in) :: t_startpt
    real(dp), intent(in) :: period_pt
    real(dp), intent(in) :: alpha_nom_startpt
    real(dp), intent(in) :: beta_nom_startpt
    real(dp), intent(in) :: gamma_nom_startpt
    type(plmission), intent(in out) :: mission
    type(sssolsys), intent(in) :: solsys
    type(ploutput), intent(in) :: missionoutput

    type(ssplanet) :: earth
    integer :: n_pt, n_sm_pt, sm
    type(gnsec) :: t_endpt, t_midpt, t_startsm, t_endsm, t_midsm
    type(gnsphangle_double) :: axes(3), axes_startpt(3), zaxis_nom_startpt
    real(dp) :: period_sm, satpos_startpt(3), satpos_endpt(3), &
      rotm_nom_startpt(3, 3), tb_startpt(3), delta_rate_z_rot_startpt, &
      delta_phase_rot_startpt, rotm(3, 3), rotm_startpt(3, 3), &
      rate_rot_startpt(3), rate_z_rot_startpt, rate_z_rot_nom_startpt, &
      rotm_delta_phase_rot(3, 3), angle_sun, angle_earth, alpha_startpt, &
      beta_startpt, gamma_startpt, earthpos_startpt(3), earthpos_endpt(3)
    real(dp), pointer :: sataxes(:, :) => null()
    type(plnutation) :: nutation

    n_sm_pt = -1 ! to shut up the compiler
    n_pt = mission%scanstrategy%n_pt

    call gn_assert(pt>0,'pl_simpointingperiod_nom: pt <= 0',pt)
    call gn_assert(pt<=n_pt,'pl_simpointingperiod_nom: pt > n_pt',pt)

    ! INITIALISATIONS.

    ! REFERENCE TIMES.

    ! Reference times for the pointing period: the pointing period,
    ! period_pt; the starting time of the pointing period, t_startpt;
    ! the middle point of the pointing period, t_midpt; and the end
    ! time of the pointing period, t_endpt.

!   period_pt = mission%scanstrategy%period_pt
!   t_startpt = mission%t_start + real(pt-1,dp) * period_pt

    call gn_assert(period_pt>=0.0, &
      'pl_simpointingperiod_nom: period_pt < 0.0', period_pt)
    t_midpt = t_startpt + period_pt / 2.0
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

      period_sm = missionoutput%period_sm
      n_sm_pt = nint(period_pt / period_sm)

      if (abs(real(n_sm_pt,dp) * period_sm - period_pt) > 1.0e-6_dp) &
        call gn_warning( &
          'pl_simpointingperiod_nom: non-integer number of samples in period')

    end if

    ! SATELLITE POSITION.

    ! Calculate the position of the satellite (in the orbit around the
    ! L2 point of the Sun-Earth system) at the start and the end of the
    ! pointing period.

    earth = ss_planet_name('Earth', solsys)
    satpos_startpt = ss_l2orbitpos(t_startpt, mission%l2orbit, earth)
    satpos_endpt = ss_l2orbitpos(t_endpt, mission%l2orbit, earth)
    earthpos_startpt = ss_plpos_orbit(t_startpt, earth%orbit)
    earthpos_endpt = ss_plpos_orbit(t_endpt, earth%orbit)

    ! SOLAR SYSTEM.

    ! INITIAL CONDITIONS OF THE POINTING PERIOD.

    ! Nominal starting values for the satellite's orientation (expressed
    ! as a rotation matrix to take the satellite from its reference
    ! orientation) and its rotation rate (just the nominal rotation rate
    ! from the scanning strategy).

    rotm_nom_startpt &
      = gn_rotm(alpha_nom_startpt, beta_nom_startpt, gamma_nom_startpt)
    zaxis_nom_startpt = gn_sphang(beta_nom_startpt, alpha_nom_startpt)

    rate_z_rot_nom_startpt = mission%scanstrategy%rate_rot

    ! Test that the Solar and Earth aspect angles are not outside
    ! the acceptable limits.  This could be done at both the start and
    ! the end of the pointing period, and for both the nominal and
    ! actual satellite orientations, although testing the actual satellite
    ! orientation at the end of the pointing period would be rather
    ! difficult as a full pointing period simulation would be required
    ! to determine the orientation of the satellite at the end of the
    ! pointing period.  Hence the tests are just done at the start of
    ! pointing period, here for the nominal values and, at the later
    ! appropriate point, for the actual values.

    angle_sun = pl_sataspectangle(rotm_nom_startpt, satpos_startpt, &
      (/0.0_dp, 0.0_dp, 0.0_dp/));
    call pl_testangle_sun(angle_sun, mission%scanstrategy)

!   angle_sun = pl_sataspectangle(rotm_nom_startpt, satpos_endpt, &
!     (/0.0_dp, 0.0_dp, 0.0_dp/));
!   call pl_testangle_sun(angle_sun, mission%scanstrategy)

    angle_earth = pl_sataspectangle(rotm_nom_startpt, satpos_startpt, &
      ss_plpos_orbit(t_startpt, earth%orbit))
    call pl_testangle_earth(angle_earth, mission%scanstrategy)

!   angle_earth = pl_sataspectangle(rotm_nom_startpt, satpos_endpt, &
!     ss_plpos_orbit(t_endpt, earth%orbit))
!   call pl_testangle_earth(angle_earth, mission%scanstrategy)

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

    ! Test that the Solar and Earth aspect angles are not outside
    ! the acceptable limits.  This could be done at both the start and
    ! the end of the pointing period, and for both the nominal and
    ! actual satellite orientations, although testing the actual satellite
    ! orientation at the end of the pointing period would be rather
    ! difficult as a full pointing period simulation would be required
    ! to determine the orientation of the satellite at the end of the
    ! pointing period.  Hence the tests are just done at the start of
    ! pointing period, here for the nominal values and, at the later
    ! appropriate point, for the actual values.

    angle_sun = pl_sataspectangle(rotm_startpt, satpos_startpt, &
      (/0.0_dp, 0.0_dp, 0.0_dp/));
    call pl_testangle_sun(angle_sun, mission%scanstrategy)

    angle_earth = pl_sataspectangle(rotm_startpt, satpos_startpt, &
      ss_plpos_orbit(t_startpt, earth%orbit))
    call pl_testangle_earth(angle_earth, mission%scanstrategy)

    ! Initialise the pointing period, primarily by calculating various
    ! quantitites related to the inertia tensor of the satellite and
    ! the initial inertial rates (having taken the initial z-axis rotation
    ! rate as an input). Unlike the previous routine, this routine is
    ! deterministic. The output values are the updated satellite structure
    ! and the three inertial rates, along with the nutation structure,
    ! although if we are doing an ideal dynamics simulation then
    ! this is not supplied and that part of the code is avoided.

    if (mission%satellite%dynamics == 'ideal') then
      call pl_pointingperiod_init(satpos_startpt, rotm_nom_startpt, &
        tb_startpt, rate_z_rot_startpt, solsys%star, mission%satellite, &
        rate_rot_startpt)
    else
      call pl_pointingperiod_init(satpos_startpt, rotm_nom_startpt, &
        tb_startpt, rate_z_rot_startpt, solsys%star, mission%satellite, &
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
          period_sm = missionoutput%period_sm
          t_startsm = t_startpt + real(sm-1,dp) * period_sm
          t_midsm = t_startsm + 0.5 * period_sm
          t_endsm = t_startsm + period_sm

          ! Get the rotation matrix for the sample time and convert to
          ! Euler angles.
          if ((mission%scanstrategy%method == 'function') &
            .and. (mission%scanstrategy%mode == 'fastprecession')) then
            rotm = pl_rotm_nom(t_midsm, mission%scanstrategy, earth)
          else
            rotm = pl_rotm_ideal(t_midsm-t_startpt, rotm_startpt, &
              rate_rot_startpt(3))
          end if

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

        ! Only ring-based scanning strategies can be evolved dynamically,

        if ((mission%scanstrategy%method == 'function') &
          .and. (mission%scanstrategy%mode == 'fastprecession')) then

          call gn_fatal( &
          'pl_simpointingperiod_nom: no dynamical model for fast precession')

        else

          ! In the analytical case get the orientation of the satellite
          ! at the start of the pointing period and then evolve it using
          ! a simple analytical/approximate integration of the Euler
          ! equation.

          do sm = 1, n_sm_pt
            ! Reference time for the sampling period.
            period_sm = missionoutput%period_sm
            t_startsm = t_startpt + real(sm-1,dp) * period_sm
            t_midsm = t_startsm + 0.5 * period_sm
            t_endsm = t_startsm + period_sm

            ! Get the rotation matrix for the sample time (given its
            ! value at the reference time at the start of the pointing
            ! period) and convert to Euler angles.
            rotm = pl_rotm_analytical(t_midsm-t_startpt, rotm_startpt, &
              rate_rot_startpt(3), nutation)

            axes = gn_rotm2axes(rotm)
            sataxes(sm, 1) = axes(1)%theta
            sataxes(sm, 2) = axes(1)%phi
            sataxes(sm, 3) = axes(2)%theta
            sataxes(sm, 4) = axes(2)%phi
            sataxes(sm, 5) = axes(3)%theta
            sataxes(sm, 6) = axes(3)%phi
          end do

        end if

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
    if (missionoutput%simulatepointings) then
      deallocate(sataxes)
    end if

  end subroutine pl_simpointingperiod_nom

end module planck_mission
