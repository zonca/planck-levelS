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

module crosspol_main

  implicit none

contains

  !======================================================================

  ! Reads beam Stokes parameters and creates an effective beam
  ! incorporating the cross-polar leakage due to the detector.  It is
  ! done by one of two possible methods:
  !
  ! (1) If only the co-polar beam is supplied, it scales the Stokes
  !     parameters of the beam to introduce the cross-polar leakage.
  !
  ! (2) If both the co-polar and cross-polar beams are supplied, it
  !     combines them to introduce the cross-polar leakage.
  !
  ! If the shapes of the co-polar and cross-polar beam are identical
  ! such that they only differ in their polarisation orientations, the
  ! two methods will produce the same result.
  !
  ! Mark Ashdown, CPAC

  function crosspol_module2(args)

    use planck_config
    use ls_misc_utils
    use announce_mod
    use ls_paramfile_io
    use focalplanemod
    use dmc_io

    use beam

    character(len=*), intent(in) :: args(:)
    integer(i4b) :: crosspol_module2
    logical :: use_cross
    real(dp) :: epsilon, angle
    character(len=20) :: beam_format
    character(len=filenamelen) :: detdbfile, co_file, cross_file, eff_file
    character(len=9) :: detid
    type(paramfile_handle) :: params
    type(bmpolar) :: co_p, cross_p, eff_p
    type(bmsquare) :: co_s, cross_s, eff_s

    ! Announce the program.

    call module_startup('crosspol',args,.true.)

    ! Initialise DMC and get parameters.

    call dmc_init(args)

    params = dmc_getrunparams()

    beam_format = 'polar'
    if (parse_string(params, 'co_file_square', default='', &
        descr='Input co-polar beam file')/='') then
      beam_format = 'square'
    endif

    co_file = parse_string(params, 'co_file_'//beam_format, &
        descr='Input co-polar beam file')
    cross_file = parse_string(params, 'cross_file_'//beam_format, &
        descr='Input cross-polar beam file', default='')
    eff_file = parse_string(params, 'eff_file_'//beam_format, &
        descr='Output effective beam file')

    detdbfile = parse_string(params, 'focalplane_db', &
        descr='focal plane database file')
    detid = trim(parse_string(params, 'detector_id', descr='detector ID'))
    call fpdb_open(detdbfile)
    epsilon = fpdb_value_real8(detid, 'epsilon')
    call fpdb_close

    use_cross = (cross_file /= ' ')

    if (use_cross) then
      angle = parse_double(params, 'angle', default=0.0_dp, &
          descr='Angle between co- and cross-polar beam coordinate systems')
    end if

    call parse_finish(params)

    write(*, '("")')

    select case (beam_format)
    case('polar')

      write(*, '("Reading co-polar beam from ", a)') trim(co_file)
      call bm_polar_read(co_p, co_file)

      if (use_cross) then

        write(*, '("Reading cross-polar beam from ", a)') trim(cross_file)
        call bm_polar_read(cross_p, cross_file)

        write(*, '("Combining beams to create effective beam")')
        call bm_polar_crosspol2(co_p, cross_p, epsilon, angle, eff_p)

      else

        write(*, '("Adjusting beam to create effective beam")')
        call bm_polar_crosspol1(co_p, epsilon, eff_p)

      end if

      write(*, '("Writing effective beam to ", a)') trim(eff_file)
      call bm_polar_write(eff_p, eff_file)

      call bm_polar_free(eff_p)
      if (use_cross) call bm_polar_free(cross_p)
      call bm_polar_free(co_p)

    case('square')

      write(*, '("Reading co-polar beam from ", a)') trim(co_file)
      call bm_square_read(co_s, co_file)

      if (use_cross) then

        write(*, '("Reading cross-polar beam from ", a)') trim(cross_file)
        call bm_square_read(cross_s, cross_file)

        write(*, '("Combining beams to create effective beam")')
        call bm_square_crosspol2(co_s, cross_s, epsilon, angle, eff_s)

      else

        write(*, '("Adjusting beam to create effective beam")')
        call bm_square_crosspol1(co_s, epsilon, eff_s)

      end if

      write(*, '("Writing effective beam to ", a)') trim(eff_file)
      call bm_square_write(eff_s, eff_file)

      call bm_square_free(eff_s)
      if (use_cross) call bm_square_free(cross_s)
      call bm_square_free(co_s)

    case default

      call exit_with_status(1, 'Unknown beam format')

    end select

    ! Shutdown DMC.

    call dmc_shutdown

    ! Set return value.

    crosspol_module2 = 0

  end function crosspol_module2

  function crosspol_module (init_object)
    use planck_config

    character(len=*), intent(in) :: init_object
    integer(i4b) :: crosspol_module
    crosspol_module = crosspol_module2 ( (/init_object/) )
  end function crosspol_module

  !======================================================================

end module crosspol_main
