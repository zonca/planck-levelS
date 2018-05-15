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

module planck_scanning
  use planck_config
  use general_const
  use general_error
  use general_vector
  use general_time
  use solarsystem_planet
  use solarsystem_orbit
  use solarsystem_solsys
  use planck_missionio
  use ls_paramfile_io
  use ls_misc_utils
  implicit none
  private

  public :: plscanstrategy, pl_scanstrategy_init, pl_rotm_nom, pl_t_startpt, &
    pl_period_pt, pl_sataspectangle, pl_testangle_sun, pl_testangle_earth

  ! Planck scan strategy, defined primarily by mode (which determines the
  ! type of scan strategy), and in all cases parameterised in ecliptic
  ! coordinates, but using standard mathematical angles (theta, phi) to
  ! specify a position on the sphere. Also note that, except at input and
  ! output, all angles are measured in radians; similarly all times
  ! periods are measured in seconds.
  !
  ! The overall scan strategy is given by mode; the rate at which the
  ! azimuthal angle of the scan progresses is given by phimode (either
  ! linear or keeping in line with the sun-earth line). The reference
  ! time is t_0; the base latitude is theta_z_0; the magnitude of the
  ! latitude variation is delta_theta_z. The period with which the
  ! satellite rotates is period_rot; the rate at which the satellite
  ! rotates is rate_rot. The pointing period (defined to be the period
  ! over which the satellite spins passively) is period_pt; the repointing
  ! period (defined to be the period between active repointings of the
  ! satellite, as distinct from the ``jitter'' that may be added to the
  ! various pointing periods within one repointing period); the period of
  ! the overall scanning motion (if there is such variation) is
  ! period_motion, and the reference phase at time t_0 of the motion
  ! is phase_motion_0. The number of pointing periods over the entire
  ! mission is n_pt.
  !
  ! Finally, the restrictions placed on the scan strategy through the
  ! maximum Solar aspect angle (i.e., the angle between the *negative*
  ! z-axis and the Sun, which should be close to zero), angle_sun_max,
  ! and the Earth aspect angle (i.e., the angle between the *negative*
  ! z-axis and the Earth, which should also be close to zero),
  ! angle_earth_max. If these limits are violated (i.e., the back of
  ! satellite makes too much of an angle with either the Sun, in which
  ! case the Solar panel might not be sufficiently shielding, or with
  ! the Earth, in which case telemetry is problematic in addition to
  ! any heating effects) then either a warning message is generated
  ! or execution is halted, depending on whether angle_sunearth_status
  ! is ``fatal'' (in which case execution is halted), ``warning'' (in which
  ! case only a warning message is generated) or ``none'' (in which case
  ! no such tests are made).
  !
  ! A 2007 modification was to allow a choice of overall method: if the
  ! scan pattern is to be calculated then method = 'function'; if it is to
  ! be read in from a file then method = 'table'.  In the latter case all
  ! the parameters

  type plscanstrategy

    ! Generic parametrs
    character(len=filenamelen) :: method
    integer :: n_pt
    real(dp) :: period_rot
    real(dp) :: rate_rot

    ! Functional scan strategy parameters
    character(len=filenamelen) :: mode
    character(len=filenamelen) :: phimode
    type(gnsec) :: t_0
    real(dp) :: theta_z_0
    real(dp) :: delta_theta_z
    real(dp) :: period_pt
    real(dp) :: period_rept
    real(dp) :: period_motion
    real(dp) :: phase_motion_0

    ! Look-up table for non-simple scan strategies
    type(gnsphangle_double), pointer :: zaxis_pts(:) => null()
    type(gnsec), pointer :: t_startpts(:) => null()
    real(dp), pointer :: period_pts(:) => null()

    ! Solar and terrestrial angle warnings
    character(len=filenamelen) :: angle_sunearth_status
    real(dp) :: angle_sun_max
    real(dp) :: angle_earth_max

  end type plscanstrategy

  interface pl_rotm_nom
    module procedure pl_rotm_nom_t, pl_rotm_nom_pt
  end interface

contains

  ! Initialise the scan strategy, also ensuring that the start and end
  ! times for the mission are consistent with the number and duration
  ! of pointing periods.
  function pl_scanstrategy_init(params,tsd,tst,t_start,t_end) &
    result(scanstrategy)
    type(paramfile_handle), intent(inout) :: params
    character(len=*), intent(in) :: tsd,tst
    type(gnsec), intent(inout) ::  t_start, t_end
    type(plscanstrategy) :: scanstrategy

    real(dp) :: n_rot_pt, n_motion_yr
    integer :: n_pt_rept, pt
    character(len=filenamelen) :: scanfilename, scandataformat, &
      string_date, string_time

    write(*, '(/,a,/)') 'Planck scan strategy'

    ! Overall scan method
    scanstrategy%method = parse_string(params,'method_scan')

    ! The (primary) rotational period of the satellite (and the *angular*
    ! rotation rate implied by this).  These have to be read in at this
    ! point because the PPL format which specifies a scan strategy does
    ! not specify a rotation rate.
    scanstrategy%period_rot = GNMIN_SEC * parse_double(params, &
      'period_rot_scan', vmin=0.0_dp)
    call gn_assert(scanstrategy%period_rot>0.0, &
      'pl_scanstrategy_init: period_rot = 0.0')
    scanstrategy%rate_rot = twopi / scanstrategy%period_rot

    ! Now split up according to which method is being used.
    if (scanstrategy%method == 'function') then

      ! FUNCTIONAL SCAN STRATEGY

      ! Set the reference time for the scan strategy (using the mission
      ! start time as the default).
      string_date = parse_string(params,'date_0_scan',default=tsd)
      string_time = parse_string(params,'time_0_scan',default=tst)
      scanstrategy%t_0 = gn_strings2sec(string_date, 'yyyymmdd', &
        string_time, 'hhmmssdsss')

      ! Read in the overall scanning mode.
      scanstrategy%mode = parse_string(params,'mode_scan')

      ! Read in the azimuthal scanning mode.
      scanstrategy%phimode = parse_string(params,'phimode_scan')

      ! Read in the fiducial colatitude and colatitude variation of the
      ! satellite's z-axis.
      scanstrategy%theta_z_0 = GNDEG_RAD * parse_double(params, &
        'theta_z_0_scan', default=90.0_dp, vmin=0.0_dp, vmax=180.0_dp)

      if (scanstrategy%mode == 'constantlatitude') then
        scanstrategy%delta_theta_z = 0.0_dp
      else
        scanstrategy%delta_theta_z &
          = GNDEG_RAD * parse_double(params,'delta_theta_z_scan', &
          default=0.0_dp, vmin=0.0_dp, vmax=90.0_dp)
      end if

      ! Read in the pointing period of the satellite (in terms of rotational
      ! periods). For most scanning strategies this should be an integer
      ! value (i.e., an integer number of rotations for each output pointing
      ! period), although this is not a formal requirement. The exception
      ! is the fast precession scan strategy, in which case there is no
      ! reason to strongly prefer an integer value.
      n_rot_pt = parse_double(params,'n_rot_pt_scan',vmin = 0.0_dp)
      call gn_assert(n_rot_pt>0.0,'pl_scanstrategy_init: n_rot_pt = 0.0')
      scanstrategy%period_pt = n_rot_pt * scanstrategy%period_rot

      ! Read in the repointing period of the satellite (in terms of pointing
      ! periods). In most cases this is likely to be unity (i.e., the
      ! repointing period between different satellite orientations is
      ! the same as the pointing period used for outputs (and between
      ! ``jitter''-style perturbations). The minimum value here is unity
      ! (as we can't have pointing periods that stretch over more than
      ! one repointing period) and the input value is an integer (for
      ! similar reasons). However this has not been thought through for
      ! long enough to preclude a change to more flexible real values
      ! at some stage.
      n_pt_rept = parse_int(params,'n_pt_rept_scan', vmin=1)
      scanstrategy%period_rept = real(n_pt_rept,dp) * scanstrategy%period_pt

      ! Read in the relevant angular and temporal scales for each scanning
      ! stratgey.
      if (scanstrategy%mode == 'constantlatitude') then
        scanstrategy%period_motion = 0.0_dp
        scanstrategy%phase_motion_0 = 0.0_dp
      else
        n_motion_yr = parse_double(params,'n_motion_yr_scan', vmin=0.0_dp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code left in to process old parameter files.  Note that, because this is
! called after the new parameter, n_motion_yr_scan, has already been
! processed, it will have no effect unless the parameter n_prec_yr_scan
! is also present.
n_motion_yr = parse_double(params,'n_prec_yr_scan', &
  default=n_motion_yr, vmin=0.0_dp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (n_motion_yr == 0.0) then
          call gn_warning( &
          'pl_scanstrategy_init: n_motion_yr = 0.0: constant latitude scanning')
          scanstrategy%mode = 'constantlatitude'
          scanstrategy%period_motion = 0.0_dp
          scanstrategy%phase_motion_0 = 0.0_dp
        else
          scanstrategy%period_motion = GNYR_SEC / n_motion_yr
          scanstrategy%phase_motion_0 = parse_double(params, &
            'phase_motion_0_scan',default=0.0_dp,vmin=0.0_dp,vmax=1.0_dp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code left in to process old parameter files.  Note that, because this is
! called after the new parameter, n_motion_yr_scan, has already been
! processed, it will have no effect unless the parameter phase_prec_0_scan
! is also present.

scanstrategy%phase_motion_0 = parse_double(params,'phase_prec_0_scan', &
  default=scanstrategy%phase_motion_0, vmin=0.0_dp, vmax=1.0_dp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end if

      end if

      ! Calculate the total number of pointing periods during the mission,
      ! remembering that the nominal mission time might correspond to a
      ! non-integer number of pointing periods.  Note also that, due to
      ! various silly restrictions (lengths of FITSIO header keywords, etc.)
      ! there is a hard upper limit that the number of pointing periods
      ! must be expressable as a six digit number (i.e., n_pt <= 999,999).
      ! This may seem like an enormous number, but given roughly 10,000
      ! hour-long repointing periods and 60 potential pointing periods
      ! per repointing period (if one wants each satellite revolution to
      ! be separate) then that's 600,000 already, and almost up to this
      ! limit.  It's enough to make one hope the bloody thing blows up
      ! on the launchpad ...
      call gn_assert(t_end>=t_start,'pl_scanstrategy_init: t_start after t_end')

      scanstrategy%n_pt = floor((t_end - t_start) / scanstrategy%period_pt)
      t_end = t_start + real(scanstrategy%n_pt,dp) * scanstrategy%period_pt

      write(*, '(/,a, f5.1, a, i6, a, f7.3, a)') &
        'Mission duration: ', (t_end - t_start)/GNDAY_SEC, ' days; ', &
        scanstrategy%n_pt, ' pointing periods of ', &
        scanstrategy%period_pt/GNHR_SEC, ' hr.'

      call gn_assert(scanstrategy%n_pt<1000000, &
        'pl_scanstrategy_init: n_pt >= 1000000')

    else if (scanstrategy%method == 'table') then

      ! LOOK-UP TABLE OF POINTINGS

      ! File details.
      scanfilename = parse_string(params,'scanfile')
      call assert_present(scanfilename)

      scandataformat = parse_string(params,'scandataformat')

      ! Read in the pointing times and axis angles from the scan
      ! file.
      call pl_readpointings_nom( &
        scanfilename, scandataformat, &
        t_start, t_end, scanstrategy%n_pt, scanstrategy%zaxis_pts, &
        scanstrategy%t_startpts, scanstrategy%period_pts)

    end if

    ! SOLAR AND TERRESTRIAL ANGLE LIMITS

    ! Read in the maximum Solar and Earth aspect angles and whether
    ! violations of these are critical (and if these errors are not
    ! checked don't bother reading in values).

    scanstrategy%angle_sunearth_status = parse_string(params, &
      'angle_sunearth_status', default='warning')
    if (scanstrategy%angle_sunearth_status /= 'none') then
      scanstrategy%angle_sun_max = GNDEG_RAD &
        * parse_double(params,'angle_sun_max', default=10.0_dp, &
        vmin=0.0_dp, vmax=180.0_dp)
      scanstrategy%angle_earth_max = GNDEG_RAD &
        * parse_double(params,'angle_earth_max', default=15.0_dp, &
        vmin=0.0_dp, vmax=180.0_dp)
    else
      scanstrategy%angle_sun_max = pi
      scanstrategy%angle_earth_max = pi
    end if

    ! CHECK THE SCAN STRATEGY

! Ideally there wouldn't be the 1e-4 at the end of this if statement
! to test whether pointing periods overlap, however it's proved
! necessary to avoid erroneous errors when there is overlap due
! to numerical imprecision.  Martin Reinecke has this flagged for
! correction, although it's not particularly clear how to fix this
! in any other way, short of ``snapping'' the pointing periods together
! earlier on in the code.

    if (scanstrategy%method == 'table') then
      do pt = 1, scanstrategy%n_pt - 1
        call gn_assert ((scanstrategy%t_startpts(pt + 1) + 1e-4_dp) &
          > (scanstrategy%t_startpts(pt)+scanstrategy%period_pts(pt)), &
          'pl_scanstrategy_init: pointing periods overlap')
      end do
    end if
  end function pl_scanstrategy_init

  ! Returns the rotation matrix that takes the satellite from its reference
  ! orientation to its nominal orientation at time t of the given scan
  ! strategy. This includes the rotation of the satellite about its own
  ! z-axis and, for ideal error-free scanning, is sufficient to define the
  ! entire mission. If there is pointing noise this must be added; if there
  ! is to be a full dynamical simulation of the passive rotation periods
  ! of the satellite then this function is only useful to work out the
  ! mean spin-axis position of a pointing period.
  function pl_rotm_nom_t(t, scanstrategy, earth) result(rotm)
    type(gnsec), intent(in) :: t
    type(plscanstrategy), intent(in) :: scanstrategy
    type(ssplanet), intent(in) :: earth
    real(dp) :: rotm(3, 3)

    integer :: n_rept_t
    type(gnsec) :: t_0, t_midrept, tpluseps
    real(dp) :: theta_z_0, delta_theta_z, period_rot, period_rept, &
      period_motion, phase_motion_0, phase_rot, phase_motion, theta_z, phi_z
    character(len=filenamelen) :: mode, method

    method = scanstrategy%method

    if (method == 'function') then

      ! INITIALISATION
      !
      ! First up make some local copies of the various parameters that
      ! define the scan strategy and calculate some derived quantities.

      ! Overall scanning mode.
      mode = scanstrategy%mode

      ! Angular scales that define the various possibilities within a
      ! given scan strategy.
      theta_z_0 = scanstrategy%theta_z_0
      delta_theta_z = scanstrategy%delta_theta_z

      ! Scan strategy reference time (not necessarily the start time of the
      ! mission, but that is the default value), from which the start, end
      ! and midpoint of any repointing period can be calculated.

      ! These times also allow various phases to be calculated which are
      ! the rotation angles that must be converted into rotation matrices.
      ! Note that mod(phase, 2 pi) part of these definitions is optional as
      ! the phases are always arguments of rotations, in which case the
      ! additional angle of the form 2 pi n is irrelevant. Note also that
      ! the phases are defined with respect to the overall mission t_0;
      ! it might be more reasonable to define zero phases for the various
      ! subcomponents of the scan strategy (e.g., the main satellite rotation)
      ! in terms of the time since, e.g., the start of the pointing period
      ! and this could easily be implemented if it was deemed desireable.
      t_0 = scanstrategy%t_0

      ! The time at the middle of the current repointing period. This is
      ! needed for most scan strategies as the azimuthal position (and
      ! latitude variation) of the satellite's nominal z-axis orientation
      ! should be calculated at this time (and not the start of the
      ! pointing period). This time, t_midrept, is calculated by first
      ! finding out which pointing period we are in, rept_t, and then
      ! subtracting the the time at the start of this repointing period
      ! from the current time. Thus
      !
      !   t_midrept = t_0 + (n_rept_t + 1/2) period_rept.
      !
      ! Note, however, that there is the possibility for error here if the
      ! time t is at the exact start of a pointing period -- rounding errors
      ! could result in the previous pointing period being assumed instead.
      ! To circumvent this, an ugliness has been introduced into the code,
      ! with a small time added to t in this function. It is chosen to be
      ! of order 10^(-6) s as this is far above the numerical resolution of
      ! the gnsec variable type, but much smaller than any timescales likely
      ! to be relevant to the Planck mission.
      tpluseps = t + 1.0e-6_dp

      period_rept = scanstrategy%period_rept
      n_rept_t = floor((tpluseps - t_0) / period_rept)
      t_midrept = t_0 + (real(n_rept_t,dp) + 0.5) * period_rept
      ! The phase of the main scanning motion. It is defined such that it has
      ! angular phase 2 pi phase_motion_0 at time t_0 and that it evolves
      ! through angle 2 pi in each period of period_motion. The definition
      ! of the zero phase is thus unambiguous for the main scanning motion.
      ! This also has to be adjusted depending on whether discrete pointing
      ! periods are used (in which case the phase should be evaluated at
      ! the middle of the pointing period) or whether a fast precession
      ! scan strategy is used (in which case the evaluation is at the
      ! present time, t).

      period_motion = scanstrategy%period_motion
      phase_motion_0 = scanstrategy%phase_motion_0
      if (period_motion == 0.0) then
        phase_motion = phase_motion_0
      else if (mode == 'fastprecession') then
        phase_motion = twopi * ((t - t_0) / period_motion + phase_motion_0)
      else
        phase_motion &
          = twopi * ((t_midrept - t_0) / period_motion + phase_motion_0)
      end if
      phase_motion = mod(phase_motion, twopi)

      ! Rotation period and phase. This is defined so that the satellite's
      ! rotation is ``continuous'' in the sense that it continues independent
      ! of any repointing manouevers; hence phase_rot = 0 at t = t_0. It is
      ! this choice of zero phase that seems most arbitrary; choosing instead
      ! phase_rot = 0 (or some other canonical value) at the start of each
      ! repointing period would seem to make just as much sense.
      period_rot = scanstrategy%period_rot
      phase_rot = mod(twopi * (t - t_0) / period_rot, twopi)

      ! SATELLITE ROTATION
      !
      ! First part of the rotation is the simple (and very nominal) rotation
      ! about the z-axis of the satellite. This only changes the x- and y-axis
      ! orientations, and gets ``superceded'' in the more detailed/realistic
      ! simulations which only let the nominal rotm act on the z-axis. In
      ! more physical terms this line sets the satellite spinning, remembering
      ! that the coordinate axes are fixed in space (so that this rotation
      ! does not affect the sense of the others).  If only this rotation
      ! matrix were applied the satellite would rotate at a steady rate
      ! about the (and its) z-axis.
      rotm = gn_rotm(phase_rot, 3)

      ! SCAN STRATEGY
      !
      ! The rest of the rotations pre-multiply the above rotation and
      ! affect all the axes, thus encoding the ``chosen'' part of the
      ! scan strategy as opposed to the more fundamental fact that the
      ! satellite rotates.
      if (mode == 'constantlatitude') then

        ! Satellite z-axis stays at constant latitude.
        theta_z = theta_z_0
        rotm = gn_rotm(theta_z, 2, rotm)

        ! Satellite z-axis moves in a step-wise fashion, remaining
        ! constant for each pointing period.
        phi_z = pl_phi_z_smooth_t(t_midrept, scanstrategy, earth)
        rotm = gn_rotm(phi_z, 3, rotm)

      else if (mode == 'sinusoidal') then

        ! The colatitude varies sinusoidally with longitude. Note that
        ! the nominal z-axis position of the satellite is in the mean
        ! postition (theta_z_0) if phase_motion is any integrer multiple
        ! of pi.
        theta_z = theta_z_0 + delta_theta_z * sin(phase_motion)
        rotm = gn_rotm(theta_z, 2, rotm)

        ! Satellite z-axis moves in a step-wise fashion, remaining
        ! constant for each pointing period.
        phi_z = pl_phi_z_smooth_t(t_midrept, scanstrategy, earth)
        rotm = gn_rotm(phi_z, 3, rotm)

      else if (mode == 'cycloidal') then

        ! The spin axis of the satellite traces out a cycloid centred on
        ! the mean colatitude, theta_z_0.  This is most simply thought
        ! of as the addition or superposition of two independent motions:
        ! the fairly steady longitudinal movement, encoded in phi_z; and
        ! a circular motion around this ``mean'' position.  The only
        ! possible ambiguity comes about due to the azimuthal motion. If
        ! this is linear then the motion traced out is a true cycloid, but
        ! if the Solar longitude is followed then the cycloid locus would
        ! only be traced if the angular rotation of the satellite about its
        ! mean orientation changes rate.  At the moment the opposite choice
        ! is taken for computational convenience: this scanning strategy
        ! is synthesized by first ``making'' the satellite's z-axis trace out
        ! a cone and then adding in the longitudinal motion. This means that
        ! the cycloidal motion must be added first, followed by the
        ! azimuthal motion.

        ! First rotation of the z-axis takes describes the magnitude and
        ! rate of the cycloidal motion.
        rotm = gn_rotm(delta_theta_z, 2, rotm)
        rotm = gn_rotm(phase_motion, 3, rotm)

        ! Then apply the same rotations as for the constant latitude case,
        ! taking the z-axis to the ecliptic plane.
        theta_z = theta_z_0
        rotm = gn_rotm(theta_z, 2, rotm)

        ! Then move around the ecliptic plane at the standard rate.
        phi_z = pl_phi_z_smooth_t(t_midrept, scanstrategy, earth)
        rotm = gn_rotm(phi_z, 3, rotm)

!       theta_z = theta_z_0 + asin(sin(delta_theta_z) * sin(phase_motion))
!       rotm = gn_rotm(theta_z, 2, rotm)
!
!       phi_z = pl_phi_z_smooth_t(t_midrept, scanstrategy, earth)
!       if ((phase_motion < GNPION2) .or. (phase_motion > GN3PION2)) then
!         phi_z = phi_z + acos(cos(delta_theta_z) / sin(theta_z))
!       else
!         phi_z = phi_z - acos(cos(delta_theta_z) / sin(theta_z))
!       end if
!       rotm = gn_rotm(phi_z, 3, rotm)

      else if (mode == 'step') then

        ! The latitude is either ``up'' or ``down'' for each half of
        ! the scanning motion period.
        if (phase_motion < pi) then
          theta_z = theta_z_0 + delta_theta_z
        else
          theta_z = theta_z_0 - delta_theta_z
        end if
        rotm = gn_rotm(theta_z, 2, rotm)

        ! Satellite z-axis moves in a step-wise fashion, remaining
        ! constant for each pointing period.
        phi_z = pl_phi_z_smooth_t(t_midrept, scanstrategy, earth)
        rotm = gn_rotm(phi_z, 3, rotm)

      else if (mode == 'fastprecession') then

        ! First rotation of the z-axis takes describes the magnitude and
        ! rate of the ``fast'' precession.
        rotm = gn_rotm(delta_theta_z, 2, rotm)
        rotm = gn_rotm(phase_motion, 3, rotm)

        ! Then apply the same rotations as for the constant latitude case,
        ! taking the z-axis to the ecliptic plane.
        theta_z = theta_z_0
        rotm = gn_rotm(theta_z, 2, rotm)

        ! Then move around the ecliptic plane at the standard rate.
        phi_z = pl_phi_z_smooth_t(t, scanstrategy, earth)
        rotm = gn_rotm(phi_z, 3, rotm)

      else

        call gn_fatal('pl_rotm_nom_t: scan strategy unknown', mode)

      end if

    else if (method == 'table') then

      ! For a table scan strategy we need a two-step method where we
      ! first convert from the input time to the pointing period which
      ! this time lies within (if it does at all); the second step is
      ! then to make a call to pl_rotm_nom_pt.
      call gn_fatal('pl_rotm_nom_t: table scan strategy not yet programmed')

    else

      call gn_fatal('pl_rotm_nom_t: method unknown', method)

    end if
  end function pl_rotm_nom_t

  ! Returns the ``average'' ecliptic longitude of the satellite's z-axis
  ! as a function of time, depending on the scanning mode. The term
  ! ``average'' is used here as i) pointing noise is not accounted for;
  ! ii) the change in z-axis during a pointing period is not accounted for;
  ! and iii) the more rapid variation in e.g., the fast precession mode is
  ! not accounted for. In all cases the satellite's z-axis points roughly
  ! away from the Sun, so it is given by the ecliptic longitude of the
  ! Sun + pi (to make the satellite point away from, rather than towards)
  ! the Sun.
  function pl_phi_z_smooth_t(t, scanstrategy, earth) result(phi_z)
    type(gnsec), intent(in) :: t
    type(plscanstrategy), intent(in) :: scanstrategy
    type(ssplanet), intent(in) :: earth
    real(dp) :: phi_z

    type(gnsec) :: t_0
    real(dp) :: pos_earth(3)
    type(gnsphangle_double) :: ang_earth

    phi_z = GNDP_MAX ! Just to shut up compiler warnings
    if (scanstrategy%method == 'function') then
      if (scanstrategy%phimode == 'linear') then
        ! For linear phimode just evolve the longitude of the satellite's
        ! z-axis linearly with time.
        t_0 = scanstrategy%t_0
        phi_z = ss_lambda_star(t_0) + pi + twopi * (t - t_0) / GNYR_SEC
      else if (scanstrategy%phimode == 'followsun') then
        ! Orient the satellite so that its z-axis is in the plane through
        ! the Sun and the satellite that is perpendicular to the ecliptic.
        pos_earth = ss_plpos_orbit(t, earth%orbit)
        ang_earth = gn_v2ang(pos_earth)
        phi_z = ang_earth%phi
      else
        call gn_fatal('pl_phi_z_smooth_t: phimode unknown', &
          scanstrategy%phimode)
      end if
    else
      call gn_fatal('pl_phi_z_smooth_t: scan method unknown', &
        scanstrategy%method)
    end if
  end function pl_phi_z_smooth_t

  ! Calculate the aspect angle between the satellite's negative z-axis
  ! and position pos_0. This is defined such that the angle is 0 if
  ! the body is directly behind the satellite's Solar panel.
  function pl_sataspectangle(rotm_sat, pos_sat, pos_0) result(aspectangle)
    real(dp), intent(in) :: rotm_sat(3, 3)
    real(dp), intent(in) :: pos_sat(3)
    real(dp), intent(in) :: pos_0(3)
    real(dp) :: aspectangle

    real(dp) :: pos_sat_rel(3)

    pos_sat_rel = pos_sat - pos_0
    aspectangle = acos(dot_product(gn_vhat(pos_sat_rel), &
      gn_vhat(rotm_sat(:, 3))))
  end function pl_sataspectangle

  ! Test whether the current position and orientation of the satellite
  ! violate the Solar aspect angle criterion. Ths Sun position is not
  ! included in this call as the Sun is the origin of the coordinate system
  ! used in the Planck package.
  subroutine pl_testangle_sun(angle_sun, scanstrategy)
    real(dp), intent(in) :: angle_sun
    type(plscanstrategy), intent(in) :: scanstrategy

    character(len=filenamelen) :: angle_sunearth_status
    real(dp) :: angle_sun_max

    angle_sunearth_status = scanstrategy%angle_sunearth_status
    angle_sun_max = scanstrategy%angle_sun_max

    ! Only do anything if the error status is not ``none''.
    if (angle_sunearth_status /= 'none') then
      ! Test to see if the Sun angle is legitimate.
      if (angle_sun > angle_sun_max) then
        if (angle_sunearth_status == 'warning') then
          call gn_warning('pl_testangle_sun: angle_sun > angle_sun_max', &
            GNRAD_DEG * angle_sun)
        else if (angle_sunearth_status == 'fatal') then
          call gn_fatal('pl_testangle_sun: angle_sun > angle_sun_max', &
            GNRAD_DEG * angle_sun)
        else
          call gn_fatal('pl_testangle_sun: angle_sunearth_status unknown', &
            angle_sunearth_status)
        end if
      end if
    end if
  end subroutine pl_testangle_sun

  ! Test whether the current position and orientation of the satellite
  ! violate the Earth aspect angle criterion.
  subroutine pl_testangle_earth(angle_earth, scanstrategy)
    real(dp), intent(in) :: angle_earth
    type(plscanstrategy), intent(in) :: scanstrategy

    character(len=filenamelen) :: angle_sunearth_status
    real(dp) :: angle_earth_max

    angle_sunearth_status = scanstrategy%angle_sunearth_status
    angle_earth_max = scanstrategy%angle_earth_max

    ! Only do anything if the error status is not ``none''.
    if (angle_sunearth_status /= 'none') then
      ! Then test to see if the Earth angle is legitimate.
      if (angle_earth > angle_earth_max) then
        if (angle_sunearth_status == 'warning') then
          call gn_warning( &
            'pl_testangle_earth: angle_earth > angle_earth_max', &
            GNRAD_DEG * angle_earth)
        else if (angle_sunearth_status == 'fatal') then
          call gn_fatal( &
            'pl_testangle_earth: angle_earth > angle_earth_max', &
            GNRAD_DEG * angle_earth)
        else
          call gn_fatal('pl_testangle_earth: angle_sunearth_status unknown', &
            angle_sunearth_status)
        end if
      end if
    end if
  end subroutine pl_testangle_earth

  ! Returns the time at which the pt'th pointing period starts.
  function pl_t_startpt(pt, t_start, scanstrategy) result(t_startpt)
    integer, intent(in) :: pt
    type(gnsec), intent(in) :: t_start
    type(plscanstrategy), intent(in) :: scanstrategy
    type(gnsec) :: t_startpt

    real(dp) :: period_pt

    call gn_assert(pt>=1,'pl_t_startpt: pt < 1', pt)
    call gn_assert(pt<=scanstrategy%n_pt,'pl_t_startpt: pt > n_pt', pt)
    if (scanstrategy%method == 'function') then
      period_pt = pl_period_pt(pt, scanstrategy)
      t_startpt = t_start + real(pt-1,dp) * period_pt
    else if (scanstrategy%method == 'table') then
      t_startpt = scanstrategy%t_startpts(pt)
    else
      call gn_fatal('pl_t_startpt: method unknown', scanstrategy%method)
      t_startpt = t_start ! Just to shut up compiler warnings
    end if
  end function pl_t_startpt

  ! Returns the duration of the pt'th pointing period.
  function pl_period_pt(pt, scanstrategy) result(period_pt)
    integer, intent(in) :: pt
    type(plscanstrategy), intent(in) :: scanstrategy
    real(dp) :: period_pt

    call gn_assert(pt>=1,'pl_period_pt: pt < 1', pt)
    call gn_assert(pt<=scanstrategy%n_pt,'pl_period_pt: pt > n_pt', pt)
    if (scanstrategy%method == 'function') then
      period_pt = scanstrategy%period_pt
    else if (scanstrategy%method == 'table') then
      period_pt = scanstrategy%period_pts(pt)
    else
      call gn_fatal('pl_period_pt: method unknown', scanstrategy%method)
      period_pt = scanstrategy%period_pt ! Just to shut up compiler warnings
    end if
  end function pl_period_pt

  ! Returns the nominal rotation matrix at the start of the pt'th
  ! pointing period.
  function pl_rotm_nom_pt(pt, t_start, scanstrategy, earth) &
    result(rotm_nom_startpt)
    integer, intent(in) :: pt
    type(gnsec), intent(in) :: t_start
    type(plscanstrategy), intent(in) :: scanstrategy
    type(ssplanet), intent(in) :: earth
    real(dp) :: rotm_nom_startpt(3, 3)

    type(gnsec) :: t_startpt

    call gn_assert(pt>=1,'pl_rotm_nom_pt: pt < 1', pt)
    call gn_assert(pt<=scanstrategy%n_pt,'pl_rotm_nom_pt: pt > n_pt', pt)
    if (scanstrategy%method == 'function') then
      t_startpt = pl_t_startpt(pt, t_start, scanstrategy)
      rotm_nom_startpt = pl_rotm_nom(t_startpt, scanstrategy, earth)
    else if (scanstrategy%method == 'table') then
      rotm_nom_startpt = gn_rotm(scanstrategy%zaxis_pts(pt)%phi, &
        scanstrategy%zaxis_pts(pt)%theta, 0.0_dp)
    else
      call gn_fatal('pl_rotm_nom_pt: method unknown', scanstrategy%method)
    end if
  end function pl_rotm_nom_pt

end module planck_scanning
