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

module beam_polar

  use planck_config
  use dmc_io

  implicit none

  ! Type to store Stokes parameters of a beam on a spherical polar
  ! (theta-phi) grid.  This type is an input to the spherical harmonic
  ! transform.
  !
  ! Internally, this package uses the 'Ludwig 3' definition for the
  ! polarisation basis with the co-polar (positive Q) direction
  ! aligned with the y-axis.

  type bmpolar
    integer :: nphi, ntheta
    real(dp) :: theta_min, theta_max ! N.B. these are in radians!
    real(dp), dimension(:,:,:), pointer :: stokes
  end type bmpolar

contains

  !======================================================================

  ! Initialise bmpolar type.

  subroutine bm_polar_init(beam, nphi, ntheta, theta_min, theta_max)
    type(bmpolar), intent(out) :: beam
    integer, intent(in) :: nphi, ntheta
    real(dp), intent(in) :: theta_min, theta_max

    beam%nphi = nphi
    beam%ntheta = ntheta
    beam%theta_min = theta_min
    beam%theta_max = theta_max
    allocate(beam%stokes(4, nphi, ntheta))

  end subroutine bm_polar_init

  !======================================================================

  ! Free bmpolar type.

  subroutine bm_polar_free(beam)
    type(bmpolar), intent(inout) ::  beam

    beam%nphi = 0
    beam%ntheta = 0
    beam%theta_min = 0.0_dp
    beam%theta_max = 0.0_dp
    deallocate(beam%stokes)

  end subroutine bm_polar_free

  !======================================================================

  ! Read bmpolar type.

  subroutine bm_polar_read(beam, filename)
    type(bmpolar), intent(inout) :: beam
    character(len=*), intent(in) :: filename

    integer :: ntheta, nphi
    real(dp) :: theta_min, theta_max
    real(dp), dimension(:), allocatable :: tmp
    type(dmc_handle) inp

    call dmc_open(inp, filename, 'beam.LS_beammap_pol')

    call dmc_get_key(inp, 'Ntheta', ntheta)
    call dmc_get_key(inp, 'Nphi', nphi)
    !call dmc_get_key(inp, 'Mintheta', theta_min)
    call dmc_get_key(inp, 'Maxtheta', theta_max)
    theta_min = 0.0_dp

    ! Initialise beam type (allocate array).

    call bm_polar_init(beam, nphi, ntheta, theta_min, theta_max)

    allocate(tmp(ntheta*nphi))
    call dmc_read_column(inp, dmc_colnum(inp, 'Beamdata'), tmp, 0_i8b)
    beam%stokes(1,:,:) = reshape(tmp, (/nphi, ntheta/))
    call dmc_read_column(inp, dmc_colnum(inp, 'BeamdataQ'), tmp, 0_i8b)
    beam%stokes(2,:,:) = reshape(tmp, (/nphi, ntheta/))
    call dmc_read_column(inp, dmc_colnum(inp, 'BeamdataU'), tmp, 0_i8b)
    beam%stokes(3,:,:) = reshape(tmp, (/nphi, ntheta/))
    call dmc_read_column(inp, dmc_colnum(inp, 'BeamdataV'), tmp, 0_i8b)
    beam%stokes(4,:,:) = reshape(tmp, (/nphi, ntheta/))
    deallocate(tmp)

    call dmc_close(inp)

  end subroutine bm_polar_read

  !======================================================================

  ! Write bmpolar type.

  subroutine bm_polar_write(beam, filename)
    type(bmpolar), intent(in) :: beam
    character(len=*), intent(in) :: filename

    integer :: npix
    real(dp), dimension(:), allocatable :: tmp(:)
    type(dmc_handle) :: out

    call dmc_create(out, filename, 'beam.LS_beammap_pol')

    call dmc_set_key(out, 'Ntheta', beam%ntheta)
    call dmc_set_key(out, 'Nphi', beam%nphi)
    call dmc_set_key(out, 'Mintheta', beam%theta_min)
    call dmc_set_key(out, 'Maxtheta', beam%theta_max)

    npix = beam%nphi*beam%ntheta
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

  end subroutine bm_polar_write

  !======================================================================

  ! Write bmpolar type to text file. Useful for visualisation!

  subroutine bm_polar_write_txt(beam, filename)
    type(bmpolar), intent(in) :: beam
    character(len=*), intent(in) :: filename

    integer :: itheta, iphi
    real(dp) :: theta, phi, dtheta, dphi

    dphi = twopi / real(beam%nphi, dp)
    dtheta = (beam%theta_max - beam%theta_min) / real(beam%ntheta-1, dp)

    open(1, file=filename, status='replace')

    ! Note angles are converted to degrees.

    do itheta = 1, beam%ntheta
      theta = beam%theta_min + real(itheta-1, dp)*dtheta
      theta = theta * 180.0_dp / pi
      write(1, *)
      do iphi = 1, beam%nphi
        phi = real(iphi-1, dp)*dphi
        phi = phi * 180.0 / pi
        write(1, *) theta, phi, beam%stokes(:, iphi, itheta)
      end do
    end do

    close(1)

  end subroutine bm_polar_write_txt

  !======================================================================

  ! Normalise beam.  The normalisation convention of the beam on input
  ! is specified by norm.  The options are "unity", "four_pi" or
  ! "eight_pi" which should be self-explanatory.  On output, the beam
  ! will be normalised to 1/2.

  subroutine bm_polar_normalise(beam, norm)
    type(bmpolar), intent(inout) :: beam
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

  end subroutine bm_polar_normalise

  !======================================================================

  ! Normalise beam by integration.  Note that this should only be used
  ! for a full-sky beam or a main beam (and only in the latter case if
  ! you are ignoring the sidelobe beam).

  subroutine bm_polar_normalise_int(beam)
    type(bmpolar), intent(inout) :: beam

    integer :: itheta, iphi
    real(dp) :: theta, sinth, dtheta, dphi, omega

    dphi = twopi / real(beam%nphi, dp)
    dtheta = (beam%theta_max - beam%theta_min) / real(beam%ntheta-1, dp)

    omega = 0.0_dp

    do itheta = 1, beam%ntheta
      theta = beam%theta_min + real(itheta-1, dp)*dtheta
      sinth = sin(theta)
      do iphi = 1, beam%nphi
        omega = omega + beam%stokes(1, iphi, itheta)*sinth*dtheta*dphi
      end do
    end do

    beam%stokes = beam%stokes * 0.5_dp / omega

  end subroutine bm_polar_normalise_int

  !======================================================================

  ! Rotates Q and U Stokes parameters from the co-cross basis to the
  ! polar basis.
  !
  ! The Q and U Stokes parameters are usually represented in the
  ! co-cross basis, where the co-polar direction is aligned with the
  ! y-axis (consistent with Ludwig 3 convention).  For the purposes of
  ! extracting the spherical harmonic coefficients, it is more useful
  ! to represent them in the polar basis.  This routine should only be
  ! called just before the spherical transform routines.

  subroutine bm_polar_stokesrotate(beam)
    type(bmpolar), intent(inout) :: beam

    integer :: iphi, itheta
    real(dp) :: phi, dphi, theta, dtheta
    real(dp) :: cos2phi, sin2phi
    real(dp) :: q, u, qp, up

    dphi = twopi / real(beam%nphi, dp)
    dtheta = (beam%theta_max - beam%theta_min) / &
        real(beam%ntheta-1, dp)

    do itheta = 1, beam%ntheta

      theta = beam%theta_min + (itheta-1) * dtheta
      if (theta == 0.0_dp) cycle

      do iphi = 1, beam%nphi

        phi = (iphi-1) * dphi
        cos2phi = cos(2.0_dp * phi)
        sin2phi = sin(2.0_dp * phi)

        ! Change co-polar basis from y to x.

        q = -beam%stokes(2, iphi, itheta)
        u = -beam%stokes(3, iphi, itheta)

        ! Rotate to polar basis.

        beam%stokes(2, iphi, itheta) =  q * cos2phi + u * sin2phi
        beam%stokes(3, iphi, itheta) = -q * sin2phi + u * cos2phi

      end do
    end do

  end subroutine bm_polar_stokesrotate

  !======================================================================

end module beam_polar
