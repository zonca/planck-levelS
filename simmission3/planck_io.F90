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

module planck_io
  use planck_config
  use ls_paramfile_io
  implicit none
  private

  public :: ploutput, pl_missionoutput_init

  ! OUTPUT PARAMETERS
  !
  ! Most important is simulatepointings (which determines whether any
  ! detailed simulation of the pointing periods, as opposed to just their
  ! basic properties is performed)
  !
  ! The number of pointings and the sample period are included
  ! in the ploutput structure as well, atlthough these are
  ! generally just derived from the scan strategy and the detector
  ! sampling frequencies.  A related piece of information is the number
  ! of pointing periods' data per pointing file (in the case of full
  ! missions), the natural default for which is one, but which can be
  ! higher if needed.

  type ploutput
    logical :: simulatepointings
    real(dp) :: period_sm
    integer :: n_pt
  end type ploutput

contains

  ! Initialises the mission output parameters, with the only inputs being
  ! the number of pointing periods and the default sampling period. The
  ! latter can be changed (and in fact is irrelevant if simulatepointings
  ! == .false.); but n_pt is rigidly set.
  function pl_missionoutput_init(params,n_pt) result(missionoutput)
    type(paramfile_handle), intent(inout) :: params
    integer, intent(in) :: n_pt
    type(ploutput) :: missionoutput

    write(*, '(/,a,/)') 'Simulation output details.'

    missionoutput%simulatepointings = parse_lgt(params, 'simulatepointings', &
      default=.false.)
    missionoutput%n_pt = n_pt
    missionoutput%period_sm = parse_double(params,'period_sm_mission', &
      default=1.0_dp, vmin=0.0_dp)
  end function pl_missionoutput_init

end module planck_io
