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

module beam_crosspol

  use ls_misc_utils
  use beam_polar
  use beam_square

  implicit none

contains

  !======================================================================

  ! Adjusts co-polar beam to create an effective beam with the
  ! required degree of cross-polar leakage.

  subroutine bm_polar_crosspol1(beam_co, epsilon, beam_eff)
    type(bmpolar), intent(in) :: beam_co
    real(dp), intent(in) :: epsilon
    type(bmpolar), intent(out) :: beam_eff

    call bm_polar_init(beam_eff, beam_co%nphi, beam_co%ntheta, &
        beam_co%theta_min, beam_co%theta_max)

    beam_eff%stokes(1, :, :) = (1.0_dp+epsilon)*beam_co%stokes(1, :, :)
    beam_eff%stokes(2, :, :) = (1.0_dp-epsilon)*beam_co%stokes(2, :, :)
    beam_eff%stokes(3, :, :) = (1.0_dp-epsilon)*beam_co%stokes(3, :, :)
    beam_eff%stokes(4, :, :) = (1.0_dp+epsilon)*beam_co%stokes(4, :, :)

  end subroutine bm_polar_crosspol1

  !======================================================================

  ! Combines co-polar beam and cross-polar beam to create an effective
  ! beam with the required degree of cross-polar leakage.  The
  ! cross-polar beam is rotated through angle "angle" (in degrees) in
  ! a right-handed sense about its z-axis before being combined with
  ! the co-polar beam.

  subroutine bm_polar_crosspol2(beam_co, beam_cx, epsilon, angle, beam_eff)
    type(bmpolar), intent(in) :: beam_co, beam_cx
    real(dp), intent(in) :: epsilon, angle
    type(bmpolar), intent(out) :: beam_eff

    integer :: itheta, iphi, iphi_dest, iphi_offset
    real(dp) :: i, q, u, v, dphi
    real(dp) :: angle_tmp, angle_error, angle_rad, cos2ang, sin2ang

    ! Check co- and croos-polar beams have same metadata.

    call assert(beam_co%nphi==beam_cx%nphi, &
      'Co- and cross-polar beams do not have same nphi')

    call assert(beam_co%ntheta==beam_cx%ntheta, &
      'Co- and cross-polar beams do not have same ntheta')

    call assert(beam_co%theta_min==beam_cx%theta_min, &
      'Co- and cross-polar beams do not have same theta_min')

    call assert(beam_co%theta_max==beam_cx%theta_max, &
      'Co- and cross-polar beams do not have same theta_max')

    ! Make angle non-negative.

    angle_tmp = angle

    do
      if (angle_tmp < 0.0_dp) then
        angle_tmp = angle_tmp + 360.0_dp
      else
        exit
      end if
    end do

    ! Check that angle corresponds to a integer multiple of dphi,
    ! where dphi is the grid spacing in the phi direction.  N.B. this
    ! calculation is done in degrees.

    dphi = 360.0_dp / real(beam_co%nphi, dp)
    iphi_offset = nint(angle_tmp/dphi)
    angle_error = angle_tmp - real(iphi_offset, dp) * dphi

    !print *, 'angle = ', angle_tmp
    !print *, 'dphi = ', dphi
    !print *, 'iphi_offset = ', iphi_offset
    !print *, 'angle_error = ', angle_error

    call assert(angle_error==0.0_dp, &
      'Angle must correspond to a whole number of steps in phi')

    ! Initialise output beam type.

    call bm_polar_init(beam_eff, beam_co%nphi, beam_co%ntheta, &
        beam_co%theta_min, beam_co%theta_max)

    ! N.B. convert angle to radians!

    angle_rad = angle_tmp * pi/180.0_dp
    cos2ang = cos(2.0_dp * angle_rad)
    sin2ang = sin(2.0_dp * angle_rad)

    ! Rotate cross-polar beam through angle.  The basis vectors rotate
    ! too, so the Q and U components are modified appropriately.  The
    ! output beam is used to store the rotated cross-polar beam.

    do itheta = 1, beam_eff%ntheta
      do iphi = 1, beam_eff%nphi

        iphi_dest = iphi + iphi_offset
        if (iphi_dest > beam_eff%nphi) then
          iphi_dest = iphi_dest - beam_eff%nphi
        end if

        i = beam_cx%stokes(1, iphi, itheta)
        q = beam_cx%stokes(2, iphi, itheta)
        u = beam_cx%stokes(3, iphi, itheta)
        v = beam_cx%stokes(4, iphi, itheta)

        beam_eff%stokes(1, iphi_dest, itheta) = i
        beam_eff%stokes(2, iphi_dest, itheta) = q*cos2ang + u*sin2ang
        beam_eff%stokes(3, iphi_dest, itheta) = u*cos2ang - q*sin2ang
        beam_eff%stokes(4, iphi_dest, itheta) = v

      end do
    end do

    ! Add co-polar beam and rotated cross-polar beam with epsilon
    ! factors.

    beam_eff%stokes = beam_co%stokes + epsilon * beam_eff%stokes

  end subroutine bm_polar_crosspol2

  !======================================================================

  subroutine bm_square_crosspol1(beam_co, epsilon, beam_eff)
    type(bmsquare), intent(in) :: beam_co
    real(dp), intent(in) :: epsilon
    type(bmsquare), intent(out) :: beam_eff

    call bm_square_init(beam_eff, beam_co%nx, beam_co%ny, &
        beam_co%xdelta, beam_co%ydelta, beam_co%xcentre, &
        beam_co%ycentre)

    beam_eff%stokes(1, :, :) = (1.0_dp+epsilon)*beam_co%stokes(1, :, :)
    beam_eff%stokes(2, :, :) = (1.0_dp-epsilon)*beam_co%stokes(2, :, :)
    beam_eff%stokes(3, :, :) = (1.0_dp-epsilon)*beam_co%stokes(3, :, :)
    beam_eff%stokes(4, :, :) = (1.0_dp+epsilon)*beam_co%stokes(4, :, :)

  end subroutine bm_square_crosspol1

  !======================================================================

  subroutine bm_square_crosspol2(beam_co, beam_cx, epsilon, angle, beam_eff)
    type(bmsquare), intent(in) :: beam_co, beam_cx
    real(dp), intent(in) :: epsilon, angle
    type(bmsquare), intent(out) :: beam_eff

    write(*, '("Warning: cross-polar leakage calculation using cross-polar ")')
    write(*, '("beam has not yet been implemented for square beams. Effective")')
    write(*, '("beam will be calculated by adjusting co-polar beam.")')

    call bm_square_crosspol1(beam_co, epsilon, beam_eff)

  end subroutine bm_square_crosspol2

  !======================================================================

end module beam_crosspol
