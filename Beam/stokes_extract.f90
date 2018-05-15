!-----------------------------------------------------------------------------
!
!  Copyright (C) 2002-2013 Mark Ashdown
!
!  This file is part of the "Beam" component of the Planck simulation
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

program stokes_extract

  use planck_config
  use announce_mod
  use ls_misc_utils
  use ls_paramfile_io
  use dmc_io

  use beam

  implicit none

  character(len=filenamelen) :: beam_format, beam_file, text_file
  character(len=filenamelen),allocatable :: args(:)
  type(paramfile_handle) :: params
  type(bmpolar) :: beamp
  type(bmsquare) :: beams

  args=getcmdline()

  ! Announce the program.

  call module_startup('stokes_extract', args, .true.)

  ! Initialise DMC and get parameters.

  call dmc_init(args)
  params = dmc_getrunparams()

  beam_format = 'polar'
  beam_file = parse_string(params, 'beam_polar', default=' ', &
      descr='Input polar beam file')

  if (beam_file == ' ') then

    beam_format = 'square'
    beam_file = parse_string(params, 'beam_square', &
        descr='Input square beam file')

  end if

  text_file = parse_string(params, 'text_file', descr='Output beam file')

  select case (beam_format)
  case ('polar')

    write(*, '("Reading Stokes parameters from ", a)') trim(beam_file)

    call bm_polar_read(beamp, beam_file)

    write(*, '("Writing Stokes parameters to ", a)') trim(text_file)

    call bm_polar_write_txt(beamp, text_file)

  case ('square')

    write(*, '("Reading Stokes parameters from ", a)') trim(beam_file)

    call bm_square_read(beams, beam_file)

    write(*, '("Writing Stokes parameters to ", a)') trim(text_file)

    call bm_square_write_txt(beams, text_file)

  end select

  ! Shutdown DMC.

  call dmc_shutdown

end program stokes_extract
