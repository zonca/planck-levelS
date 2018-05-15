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

module planck_scanning_new
  use planck_config
  use general_const
  use general_error
  use general_vector
  use general_time
  use planck_missionio_new
  use ls_paramfile_io
  use ls_misc_utils
  implicit none
  private

  public :: plscanstrategy, pl_scanstrategy_init, pl_scanstrategy_finish, &
    pl_rotm_nom, pl_t_startpt, pl_period_pt

  ! Planck scan strategy, parameterised in ecliptic
  ! coordinates, but using standard mathematical angles (theta, phi) to
  ! specify a position on the sphere. Also note that, except at input and
  ! output, all angles are measured in radians; similarly all times
  ! are measured in seconds.
  !
  ! The rate at which the satellite rotates is rate_rot.
  ! The number of pointing periods over the entire mission is n_pt.

  type plscanstrategy

    ! Generic parameters
    integer :: n_pt
    real(dp) :: rate_rot

    ! Look-up table for non-simple scan strategies
    type(gnsphangle_double), pointer :: zaxis_pts(:) => null()
    type(gnsec), pointer :: t_startpts(:) => null()
    real(dp), pointer :: period_pts(:) => null()

  end type plscanstrategy

contains

  ! Initialise the scan strategy, also ensuring that the start and end
  ! times for the mission are consistent with the number and duration
  ! of pointing periods.
  function pl_scanstrategy_init(params) &
    result(scanstrategy)
    type(paramfile_handle), intent(inout) :: params
    type(plscanstrategy) :: scanstrategy

    integer :: pt
    character(len=filenamelen) :: scanfilename, scandataformat
    real(dp) :: period_rot

    write(*, '(/,a,/)') 'Planck scan strategy'

    ! The (primary) rotational period of the satellite (and the *angular*
    ! rotation rate implied by this).  These have to be read in at this
    ! point because the PPL format which specifies a scan strategy does
    ! not specify a rotation rate.
    period_rot = GNMIN_SEC * parse_double(params,'period_rot_scan', vmin=0.0_dp)
    call gn_assert(period_rot>0.0, 'pl_scanstrategy_init: period_rot = 0.0')
    scanstrategy%rate_rot = twopi / period_rot

    ! LOOK-UP TABLE OF POINTINGS

    ! File details.
    scanfilename = parse_string(params,'scanfile')
    call assert_present(scanfilename)

    scandataformat = parse_string(params,'scandataformat')

    ! Read in the pointing times and axis angles from the scan file.
    call pl_readpointings_nom(scanfilename, scandataformat, params, &
      scanstrategy%n_pt, scanstrategy%zaxis_pts, &
      scanstrategy%t_startpts, scanstrategy%period_pts)

    ! CHECK THE SCAN STRATEGY

! Ideally there wouldn't be the 1e-4 at the end of this if statement
! to test whether pointing periods overlap, however it's proved
! necessary to avoid erroneous errors when there is overlap due
! to numerical imprecision.  Martin Reinecke has this flagged for
! correction, although it's not particularly clear how to fix this
! in any other way, short of ``snapping'' the pointing periods together
! earlier on in the code.

    do pt = 1, scanstrategy%n_pt - 1
      call gn_assert ((scanstrategy%t_startpts(pt + 1) + 1e-4_dp) &
        > (scanstrategy%t_startpts(pt)+scanstrategy%period_pts(pt)), &
        'pl_scanstrategy_init: pointing periods overlap')
    end do
  end function pl_scanstrategy_init

  ! Returns the time at which the pt'th pointing period starts.
  function pl_t_startpt(pt, scanstrategy) result(t_startpt)
    integer, intent(in) :: pt
    type(plscanstrategy), intent(in) :: scanstrategy
    type(gnsec) :: t_startpt

    call gn_assert(pt>=1,'pl_t_startpt: pt < 1', pt)
    call gn_assert(pt<=scanstrategy%n_pt,'pl_t_startpt: pt > n_pt', pt)
    t_startpt = scanstrategy%t_startpts(pt)
  end function pl_t_startpt

  ! Returns the duration of the pt'th pointing period.
  function pl_period_pt(pt, scanstrategy) result(period_pt)
    integer, intent(in) :: pt
    type(plscanstrategy), intent(in) :: scanstrategy
    real(dp) :: period_pt

    call gn_assert(pt>=1,'pl_period_pt: pt < 1', pt)
    call gn_assert(pt<=scanstrategy%n_pt,'pl_period_pt: pt > n_pt', pt)
    period_pt = scanstrategy%period_pts(pt)
  end function pl_period_pt

  ! Returns the nominal rotation matrix at the start of the pt'th
  ! pointing period.
  function pl_rotm_nom(pt, scanstrategy) &
    result(rotm_nom_startpt)
    integer, intent(in) :: pt
    type(plscanstrategy), intent(in) :: scanstrategy
    real(dp) :: rotm_nom_startpt(3, 3)

    call gn_assert(pt>=1,'pl_rotm_nom: pt < 1', pt)
    call gn_assert(pt<=scanstrategy%n_pt,'pl_rotm_nom: pt > n_pt', pt)
    rotm_nom_startpt = gn_rotm(scanstrategy%zaxis_pts(pt)%phi, &
      scanstrategy%zaxis_pts(pt)%theta, 0.0_dp)
  end function pl_rotm_nom

  subroutine pl_scanstrategy_finish (scanstrategy)
    type(plscanstrategy), intent(inout) :: scanstrategy

    if (associated(scanstrategy%zaxis_pts)) &
      deallocate (scanstrategy%zaxis_pts, scanstrategy%t_startpts, &
        scanstrategy%period_pts)
  end subroutine pl_scanstrategy_finish


end module planck_scanning_new
