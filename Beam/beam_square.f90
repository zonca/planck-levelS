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

module beam_square

  use planck_config
  use dmc_io

  implicit none

  ! Type to store Stokes parameters of a beam on a square (x-y) grid.
  !
  ! Internally, this package uses the 'Ludwig 3' definition for the
  ! polarisation basis with the co-polar (positive Q) direction
  ! aligned with the y-axis.

  type bmsquare
    integer :: nx, ny, xcentre, ycentre
    real(dp) :: xdelta, ydelta
    real(dp), dimension(:,:,:), pointer :: stokes
  end type bmsquare

contains

  !======================================================================

  subroutine bm_square_init(beam, nx, ny, xdelta, ydelta, xcentre, ycentre)
    type(bmsquare), intent(inout) :: beam
    integer, intent(in) :: nx, ny, xcentre, ycentre
    real(dp), intent(in) :: xdelta, ydelta

    beam%nx = nx
    beam%ny = ny
    beam%xdelta = xdelta
    beam%ydelta = ydelta
    beam%xcentre = xcentre
    beam%ycentre = ycentre
    allocate(beam%stokes(4, nx, ny))

  end subroutine bm_square_init

  !======================================================================

  subroutine bm_square_free(beam)
    type(bmsquare), intent(inout) :: beam

    beam%nx = 0
    beam%ny = 0
    beam%xdelta = 0.0_dp
    beam%ydelta = 0.0_dp
    beam%xcentre = 0
    beam%ycentre = 0
    deallocate(beam%stokes)

  end subroutine bm_square_free

  !======================================================================

  subroutine bm_square_read(beam, filename)
    type(bmsquare), intent(inout) :: beam
    character(len=*), intent(in) :: filename

    integer :: nx, ny, xcentre, ycentre
    real(dp) :: xdelta, ydelta
    real(dp), dimension(:), allocatable :: tmp
    type(dmc_handle) :: inp

    call dmc_open(inp, filename, 'beam.LS_cart_beammap_pol')

    call dmc_get_key(inp, 'Nx', nx)
    call dmc_get_key(inp, 'Ny', ny)
    call dmc_get_key(inp, 'Xcentre', xcentre)
    call dmc_get_key(inp, 'Ycentre', ycentre)
    call dmc_get_key(inp, 'Xdelta', xdelta)
    call dmc_get_key(inp, 'Ydelta', ydelta)

    call bm_square_init(beam, nx, ny, xdelta, ydelta, xcentre, ycentre)

    allocate(tmp(nx*ny))
    call dmc_read_column(inp, dmc_colnum(inp, 'Beamdata'), tmp, 0_i8b)
    beam%stokes(1,:,:) = reshape(tmp, (/nx,ny/))
    call dmc_read_column(inp, dmc_colnum(inp, 'BeamdataQ'), tmp, 0_i8b)
    beam%stokes(2,:,:) = reshape(tmp, (/nx,ny/))
    call dmc_read_column(inp, dmc_colnum(inp, 'BeamdataU'), tmp, 0_i8b)
    beam%stokes(3,:,:) = reshape(tmp, (/nx,ny/))
    call dmc_read_column(inp, dmc_colnum(inp, 'BeamdataV'), tmp, 0_i8b)
    beam%stokes(4,:,:) = reshape(tmp, (/nx,ny/))
    deallocate(tmp)

    call dmc_close(inp)

  end subroutine bm_square_read

  !======================================================================

  subroutine bm_square_write(beam, filename)
    type(bmsquare), intent(in) :: beam
    character(len=*), intent(in) :: filename

    integer :: npix
    real(dp), dimension(:), allocatable :: tmp
    type(dmc_handle) :: out

    call dmc_create(out, filename, 'beam.LS_cart_beammap_pol')

    call dmc_set_key(out, 'Nx', beam%nx)
    call dmc_set_key(out, 'Ny', beam%ny)
    call dmc_set_key(out, 'Xcentre', beam%xcentre)
    call dmc_set_key(out, 'Ycentre', beam%ycentre)
    call dmc_set_key(out, 'Xdelta', beam%xdelta)
    call dmc_set_key(out, 'Ydelta', beam%ydelta)

    npix = beam%nx*beam%ny
    allocate(tmp(npix))
    tmp = reshape(beam%stokes(1,:,:), (/npix/))
    call dmc_append_column(out, dmc_colnum(out, 'Beamdata'), tmp)
    tmp = reshape(beam%stokes(2,:,:), (/npix/))
    call dmc_append_column(out, dmc_colnum(out, 'BeamdataQ'), tmp)
    tmp = reshape(beam%stokes(3,:,:), (/npix/))
    call dmc_append_column(out, dmc_colnum(out, 'BeamdataU'), tmp)
    tmp = reshape(beam%stokes(4,:,:), (/npix/))
    call dmc_append_column(out, dmc_colnum(out, 'BeamdataV'), tmp)
    deallocate(tmp)

    call dmc_close(out)

  end subroutine bm_square_write

  !======================================================================

  ! Write bmsquare type to test file. Useful for visualisation!

  subroutine bm_square_write_txt(beam, filename)
    type(bmsquare), intent(in) :: beam
    character(len=*), intent(in) :: filename

    integer :: i, j
    real(dp) :: x, y

    open(1, file=filename, status='replace')

    do j = 1, beam%ny
      y = real(j-beam%ycentre, dp) * beam%ydelta
      write(1, *)
      do i = 1, beam%nx
        x = real(i-beam%xcentre, dp) * beam%xdelta
        write(1, *) x, y, beam%stokes(:, i, j)
      end do
    end do

    close(1)

  end subroutine bm_square_write_txt

  !======================================================================

  ! Normalise beam.  The normalisation convention of the beam on input
  ! is specified by norm.  The options are "unity", "four_pi" or
  ! "eight_pi" which should be self-explanatory.  On output, the beam
  ! will be normalised to 1/2.

  subroutine bm_square_normalise(beam, norm)
    type(bmsquare), intent(inout) :: beam
    character(len=20) :: norm

    select case (norm)
    case('unity')

      beam%stokes = beam%stokes * 0.5_dp

    case ('four_pi')

      beam%stokes = beam%stokes * 0.5_dp / fourpi

    case ('eight_pi')

      beam%stokes = beam%stokes * 0.25_dp / fourpi

    case default

      call exit_with_status(1, 'Unknown normalisation convention')

    end select

  end subroutine bm_square_normalise

  !======================================================================

  ! Normalise beam by integration.  Note this should only be used for
  ! a main beam (and only then if you are ignoring the slidelobe beam).

  subroutine bm_square_normalise_int(beam)
    type(bmsquare), intent(inout) :: beam

    real(dp) :: omega, dxdy

    dxdy = beam%xdelta*beam%ydelta
    omega = sum(beam%stokes(1, :, :)) * dxdy
    beam%stokes = beam%stokes * 0.5_dp / omega

  end subroutine bm_square_normalise_int

  !======================================================================

end module beam_square
