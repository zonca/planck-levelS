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

module grasp2stokes_main

  implicit none

contains

  !======================================================================

  ! Reads polarised beam amplitudes in Grasp 8 "cut" or "grd" format
  ! and converts them to Stokes parameters.
  !
  ! Mark Ashdown, CPAC.

  function grasp2stokes_module2(args)

    use planck_config
    use ls_misc_utils
    use announce_mod
    use ls_paramfile_io
    use dmc_io

    use beam

    character(len=*), intent(in) :: args(:)
    integer(i4b) :: grasp2stokes_module2
    character(len=20) :: grasp_format, grasp_copol, grasp_norm
    character(len=filenamelen) :: grasp_file, stokes_file
    type(paramfile_handle) :: params
    type(bmgrid) :: grasp_g
    type(bmcut) :: grasp_c
    type(bmpolar) :: stokes_p
    type(bmsquare) :: stokes_s

    ! Announce the program.

    call module_startup('grasp2stokes',args,.true.)

    ! Initialise DMC and get parameters.

    call dmc_init(args)

    params = dmc_getrunparams()

    grasp_file = parse_string(params, 'grasp_file', descr='Input Grasp file')
    grasp_format = trim(parse_string(params, 'grasp_format', &
        descr='Grasp file format'))
    grasp_copol = trim(parse_string(params, 'grasp_copol', &
        descr='Grasp file co-polar direction'))
    grasp_norm = trim(parse_string(params, 'grasp_norm', &
        descr='Grasp beam normalisation'))

    write(*, '("")')

    select case (grasp_format)

    case ('grd_polar')

      ! Grasp file is in polar grid format. Read .grd file and convert
      ! to Stokes parameters on polar grid.

      stokes_file = parse_string(params, 'stokes_file_polar', &
          descr='Output Stokes parameters file')

      write(*, '("Reading Grasp beam from ", a)') trim(grasp_file)

      call bm_grid_read(grasp_g, grasp_file)

      write(*, '("Converting Grasp beam to Stokes parameters")')

      call bm_grid2polar(grasp_g, stokes_p, grasp_copol)

      write(*, '("Normalising beam")')

      call bm_polar_normalise(stokes_p, grasp_norm)

      write(*, '("Writing Stokes beam to ", a)') trim(stokes_file)

      call bm_polar_write(stokes_p, stokes_file)

      call bm_polar_free(stokes_p)

      call bm_grid_free(grasp_g)

    case ('grd_square')

      ! Grasp file is in square grid format. Read .grd file and
      ! convert to Stokes parameters on square grid.

      stokes_file = parse_string(params, 'stokes_file_square', &
          descr='Output Stokes parameters file')

      write(*, '("Reading Grasp beam from ", a)') trim(grasp_file)

      call bm_grid_read(grasp_g, grasp_file)

      write(*, '("Converting Grasp beam to Stokes parameters")')

      call bm_grid2square(grasp_g, stokes_s, grasp_copol)

      write(*, '("Normalising beam")')

      call bm_square_normalise(stokes_s, grasp_norm)

      write(*, '("Writing Stokes beam to ", a)') trim(stokes_file)

      call bm_square_write(stokes_s, stokes_file)

      call bm_square_free(stokes_s)

      call bm_grid_free(grasp_g)

    case ('cut')

      ! Grasp file is in cut format. Read .cut file and convert to
      ! Stokes parameters on polar grid.

      stokes_file = parse_string(params, 'stokes_file_polar', &
          descr='Output Stokes parameters file')

      write(*, '("Reading Grasp beam from ", a)') trim(grasp_file)

      call bm_cut_read(grasp_c, grasp_file)

      write(*, '("Converting Grasp beam to Stokes parameters")')

      call bm_cut2polar(grasp_c, stokes_p, grasp_copol)

      write(*, '("Normalising beam")')

      call bm_polar_normalise(stokes_p, grasp_norm)

      write(*, '("Writing Stokes beam to ", a)') trim(stokes_file)

      call bm_polar_write(stokes_p, stokes_file)

      call bm_polar_free(stokes_p)

      call bm_cut_free(grasp_c)

    case default

      ! Unknown Grasp sub-format.

      call exit_with_status(1, 'Unknown format of Grasp beam')

    end select

    call parse_finish(params)

    ! Shutdown DMC.

    call dmc_shutdown

    ! Set return value.

    grasp2stokes_module2 = 0

  end function grasp2stokes_module2

  function grasp2stokes_module (init_object)
    use planck_config

    character(len=*), intent(in) :: init_object
    integer(i4b) :: grasp2stokes_module
    grasp2stokes_module = grasp2stokes_module2 ( (/init_object/) )
  end function grasp2stokes_module

  !======================================================================

end module grasp2stokes_main
