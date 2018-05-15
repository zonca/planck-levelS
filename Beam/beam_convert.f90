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

module beam_convert

  ! Beam conversion module. This contains routines for converting the
  ! beams from the raw-ish Grasp amplitudes in the bmgrid and bmcut
  ! types to Stokes parameters in the bmpolar and bmsquare types.
  ! bmpolar contains the Stokes parameters on a spherical polar
  ! (theta, phi) grid, and bmsquare on a square (u, v) grid. These
  ! routines do no interpolation; the conversions are exact.
  !
  ! Internally, this package uses the 'Ludwig 3' definition for the
  ! polarisation basis with the co-polar (positive Q) direction
  ! aligned with the y-axis.
  !
  ! Mark Ashdown, CPAC

  use planck_config

  use ls_misc_utils
  use beam_grid
  use beam_cut
  use beam_polar
  use beam_square

  implicit none

contains

  !======================================================================

  ! Converts beam in square grid format into Stokes parameters on a
  ! square grid.  The value of copol specifies the alignment of the
  ! co-polar basis ('x' or 'y') of the input Grasp file.

  subroutine bm_grid2square(beam_g, beam_s, copol)
    type(bmgrid), intent(in) :: beam_g
    type(bmsquare), intent(inout) :: beam_s
    character, intent(in) :: copol

    integer :: nx, ny, xcentre, ycentre, i, j, sign
    real(dp) :: xdelta, ydelta, modc2, modx2
    complex(dp) :: c, x, acaxs

    ! Check metadata of input beam.

    call assert(beam_g%ncomp==2, &
      'Error in bm_grid2square: beam has wrong number of components')

    call assert(beam_g%kgrid==1, &
      'Error in bm_grid2square: beam is not on u-v grid')

    nx = beam_g%nx
    ny = beam_g%ny

    xdelta = (beam_g%xe - beam_g%xs) / real(nx-1, dp)
    ydelta = (beam_g%ye - beam_g%ys) / real(ny-1, dp)

    xcentre = beam_g%ix - nint(beam_g%xs / xdelta) + 1
    ycentre = beam_g%iy - nint(beam_g%ys / ydelta) + 1

    call bm_square_init(beam_s, nx, ny, xdelta, ydelta, xcentre, &
        ycentre)

    select case (beam_g%kcomp)
    case (3)

      ! Beam is expressed in linear co and cx components.

      select case (copol)
      case ('x')
        sign = -1
      case('y')
        sign = 1
      case default
        call exit_with_status(1, &
            'Error in bm_grid2square: unknown value for copol')
      end select

      do j = 1, ny
        do i = 1, nx

          c = beam_g%amp(1, i, j)
          x = beam_g%amp(2, i, j)

          modc2 = abs(c)**2
          modx2 = abs(x)**2
          acaxs = c * conjg(x)

          beam_s%stokes(1, i, j) = modc2 + modx2
          beam_s%stokes(2, i, j) = sign * (modc2 - modx2)
          beam_s%stokes(3, i, j) = sign * 2.0_dp * real(acaxs)
          beam_s%stokes(4, i, j) = 2.0_dp * aimag(acaxs)

          end do
        end do

    case (9)

      ! Beam is expressed in |E| and sqrt(rhc/lhc) format.  Ignoring
      ! polarisation!

      do j = 1, ny
        do i = 1, nx

          c = beam_g%amp(1, i, j)
          modc2 = abs(c)**2

          beam_s%stokes(1, i, j) = modc2
          beam_s%stokes(2, i, j) = 0.0_dp
          beam_s%stokes(3, i, j) = 0.0_dp
          beam_s%stokes(4, i, j) = 0.0_dp

        end do
      end do

    case default

      call exit_with_status(1,'Error in bm_grid2square:&
          & beam is not in supported grid sub-format')

    end select

  end subroutine bm_grid2square

  !======================================================================

  ! Converts beam in polar grid format into Stokes parameters on a
  ! polar grid.  The value of copol specifies the alignment of the
  ! co-polar basis ('x' or 'y') of the input Grasp file.

  subroutine bm_grid2polar(beam_g, beam_p, copol)
    type(bmgrid), intent(in) :: beam_g
    type(bmpolar), intent(inout) :: beam_p
    character, intent(in) :: copol

    logical :: swaptheta
    integer :: nphi, ntheta, iphi, itheta, itheta2, sign
    real(dp) :: theta_min, theta_max, dphi
    real(dp) :: modc2, modx2
    complex(dp) :: c, x, acaxs

    ! Check metadata of input beam.

    call assert(beam_g%ncomp==2, &
      'Error in bm_grid2polar: beam has wrong number of components')

    call assert(beam_g%kgrid==7, &
      'Error in bm_grid2polar: beam is not on theta-phi grid')

    call assert(abs(beam_g%xs)<=1e-5_dp, &
      'Error in bm_grid2polar: phi coordinate does not start at zero')

    ! Amplitudes at phi = 0 degrees are repeated at phi = 360 degrees,
    ! so nphi is set to nx - 1.

    call assert(abs(beam_g%xe-beam_g%xs-360._dp)<=1e-5_dp, &
      'Error in bm_grid2polar: phi range is not 360 degrees')

    nphi = beam_g%nx-1
    ntheta = beam_g%ny

    ! Note theta min and max angles are in radians.

    theta_min = beam_g%ys * pi / 180.0_dp
    theta_max = beam_g%ye * pi / 180.0_dp

    swaptheta = theta_min>theta_max

    if (swaptheta) then
      print *,"Warning: swapping theta direction"
      theta_min = beam_g%ye * pi / 180.0_dp
      theta_max = beam_g%ys * pi / 180.0_dp
    endif

    call bm_polar_init(beam_p, nphi, ntheta, theta_min, theta_max)

    dphi = twopi / real(nphi, dp)

    select case (beam_g%kcomp)
    case (3)

      ! Beam is expressed in linear co and cx components.

      select case (copol)
      case ('x')
        sign = -1
      case ('y')
        sign = 1
      case default
        call exit_with_status(1, &
            'Error in bm_grid2polar: unknown value for copol')
      end select

      do itheta = 1, ntheta
        itheta2=itheta
        if (swaptheta) itheta2=ntheta+1-itheta
        do iphi = 1, nphi

          c = beam_g%amp(1, iphi, itheta)
          x = beam_g%amp(2, iphi, itheta)

          modc2 = abs(c)**2
          modx2 = abs(x)**2
          acaxs = c * conjg(x)

          beam_p%stokes(1, iphi, itheta2) = modc2 + modx2
          beam_p%stokes(2, iphi, itheta2) = sign * (modc2 - modx2)
          beam_p%stokes(3, iphi, itheta2) = sign * 2.0_dp * real(acaxs)
          beam_p%stokes(4, iphi, itheta2) = 2.0_dp * aimag(acaxs)

        end do
      end do

    case (9)

      ! Beam is expressed in |E| and sqrt(rhc/lhc) format.  Ignoring
      ! polarisation!

      do itheta = 1, ntheta
        itheta2=itheta
        if (swaptheta) itheta2=ntheta+1-itheta
        do iphi = 1, nphi

          c = beam_g%amp(1, iphi, itheta)
          modc2 = abs(c)**2

          beam_p%stokes(1, iphi, itheta2) = modc2
          beam_p%stokes(2, iphi, itheta2) = 0.0_dp
          beam_p%stokes(3, iphi, itheta2) = 0.0_dp
          beam_p%stokes(4, iphi, itheta2) = 0.0_dp

        end do
      end do

    case default

      call exit_with_status(1,'Error in bm_grid2square:&
          & beam is not in supported grid sub-format')

    end select

  end subroutine bm_grid2polar

  !======================================================================

  ! Converts beam in "cut" format to Stokes parameters on a polar
  ! grid.  Assumes that cuts are evenly spaced in theta.  The value of
  ! copol specifies the alignment of the co-polar basis ('x' or 'y')
  ! of the input Grasp file.

  subroutine bm_cut2polar(beam_c, beam_p, copol)
    type(bmcut), intent(in) :: beam_c
    type(bmpolar), intent(inout) :: beam_p
    character, intent(in) :: copol

    integer :: ncut, np, nphi, ntheta, icut, iphi, itheta, sign
    real(dp) :: theta_min, theta_max, dphi
    real(dp) :: modc2, modx2
    complex(dp) :: c, x, acaxs
    complex(dp), dimension(:,:,:), allocatable :: amp_tmp

    ! Check metadata of input beam.

    call assert(beam_c%icomp==3, &
      'Error in bm_cut2polar: beam is not in linear co and cx components')

    call assert(beam_c%icon==1, &
      'Error in bm_cut2polar: beam is not in phi cuts')

    call assert(beam_c%ncomp==2, &
      'Error in bm_cut2polar: beam has the wrong number of components')

    ! Work out the number of theta and phi values in the cuts.

    ncut = beam_c%ncut
    np = beam_c%np

    nphi = 2 * ncut
    ntheta = (np + 1) / 2

    ! Note theta min and max angles are in radians.

    theta_min = 0.0_dp
    theta_max = abs(beam_c%sa) * pi / 180.0_dp

    call bm_polar_init(beam_p, nphi, ntheta, theta_min, theta_max)

    ! Reshape amplitude array into a temporary array.

    allocate(amp_tmp(2, nphi, ntheta))

    do icut = 1, ncut
      amp_tmp(:, icut, :) = beam_c%amp(:, ntheta:np, icut)
      amp_tmp(:, ncut + icut, :) = beam_c%amp(:, ntheta:1:-1, icut)
    end do

    ! Convert amplitudes into Stokes parameters.

    dphi = twopi / real(nphi, dp)

    select case (copol)
    case ('x')
      sign = -1
    case ('y')
      sign = 1
    case default
      call exit_with_status(1, &
          'Error in bm_cut2polar: unknown value for copol')
    end select

    do itheta = 1, ntheta
      do iphi = 1, nphi

        c = amp_tmp(1, iphi, itheta)
        x = amp_tmp(2, iphi, itheta)

        modc2 = abs(c)**2
        modx2 = abs(x)**2
        acaxs = c * conjg(x)

        beam_p%stokes(1, iphi, itheta) = modc2 + modx2
        beam_p%stokes(2, iphi, itheta) = sign * (modc2 - modx2)
        beam_p%stokes(3, iphi, itheta) = sign * 2.0_dp * real(acaxs)
        beam_p%stokes(4, iphi, itheta) = 2.0_dp * aimag(acaxs)

      end do
    end do

    deallocate(amp_tmp)

  end subroutine bm_cut2polar

  !======================================================================

end module beam_convert
