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

module simmission4_main

contains

function simmission4_module2(args)
  use planck_config
  use general_rand
  use planck_mission_new
  use planck_io
  use levels_output
  use announce_mod
  use ls_paramfile_io
  use dmc_io
  implicit none

  character(len=*), intent(in) :: args(:)
  integer(i4b) :: simmission4_module2

  type(plmission) :: mission
  type(ploutput) :: missionoutput
  type(paramfile_handle) :: params

  call module_startup('simmission4',args,.true.)
  call dmc_init(args)
  params = dmc_getrunparams()

  ! Initialise the program.
  write(*, '(a,/)') 'GENERAL INITIALISATIONS'
  call gn_rand_init(params)
  write(*, '(a)')

  ! Initialise the mission model.
  mission = pl_mission_init(params)

  ! Calculate the default sampling period and number of pointing periods
  ! and use it to initialse the output options.
  missionoutput = pl_missionoutput_init(params,mission%scanstrategy%n_pt)

  call lsoutput_init(params,mission%scanstrategy%n_pt)

  ! Now simulate the whole mission, writing out the desired data whilst
  ! going along.
  call pl_simmission(mission, missionoutput)

  call pl_mission_finish(mission)

  call lsoutput_finish
  call parse_finish(params)

  simmission4_module2=0
end function simmission4_module2

function simmission4_module (init_object)
  use planck_config

  character(len=*), intent(in) :: init_object
  integer(i4b) :: simmission4_module
  simmission4_module = simmission4_module2 ( (/init_object/) )
end function simmission4_module

end module simmission4_main
