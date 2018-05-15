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

module general_vector
  use planck_config
  use general_const
  use general_error
  use general_maths
  implicit none
  private

  public :: gnsphangle_double, gn_sphang, gn_rotm, gn_absv, gn_v2ang, &
    gn_vhat, gn_rotm2axes, gn_rotm2ang

  ! Angular position on the sphere with colatitude theta (radians) and
  ! azimuthal angle phi (radians).
  type gnsphangle_double
    real(dp) :: theta, phi
  end type gnsphangle_double

  ! Calculates a rotation matrix.
  interface gn_rotm
    module procedure gn_rotm_euler, gn_rotm_axis, gn_rotm_multaxis, gn_rotm_tb
  end interface

contains

  ! Returns a double precision spherical angle for which 0 <= theta <= pi
  ! and 0 <= phi < 2 pi.
  function gn_ang2ang(ang_in) result(ang_out)
    type(gnsphangle_double), intent(in) :: ang_in
    type(gnsphangle_double) :: ang_out

    real(dp) :: theta, phi

    ! Get values from the structure.
    theta = ang_in%theta
    phi = ang_in%phi

    ! First of all get theta into the range 0 <= theta < 2 pi.
    theta = theta - twopi * real(floor(theta/twopi),dp)

    ! Then get it into the range 0 <= theta < pi (which implies changing
    ! phi as well).
    if (theta > pi) then
      theta = twopi - theta
      phi = phi + pi
    end if

    ! Finally get phi into the range 0 <= phi < 2 pi.
    phi = phi - twopi * real(floor(phi/twopi),dp)

    ! Put the values into the output structure.
    ang_out%theta = theta
    ang_out%phi = phi
  end function gn_ang2ang

  ! Returns the double precision unit vector vhat associated with the
  ! double precision vector of arbitrary length, v. If |v| = 0.0 then
  ! a warning message results and vhat is set to (1.0, 0.0, 0.0, ...).
  function gn_vhat(v) result(vhat)
    real(dp), intent(in) :: v(:)
    real(dp) :: vhat(size(v))

    real(dp) :: absv

    absv = gn_absv(v)
    call gn_assert(absv>0, 'gn_vhat: |v| = 0.0')
    vhat = v/absv
  end function gn_vhat

  ! Calculates the angle on the sphere, ang, of a 3-dimensional double
  ! precision vector, v. If |v| = 0.0 then a warning message is generated
  ! and (ang%theta, ang%phi) is set to (0.0, 0.0).
  function gn_v2ang(v) result(ang)
    real(dp), intent(in) :: v(:)
    type(gnsphangle_double) :: ang

    integer :: nrow
    real(dp) :: absv

    nrow = size(v)
    call gn_assert (nrow==3, 'gn_v2ang: nrow /= 3', nrow)
    absv = gn_absv(v)
    call gn_assert(absv>0, 'gn_v2ang: |v| = 0.0')
    ang%theta = acos(v(3) / absv)
    ang%phi = gn_atan(v(1), v(2))
  end function gn_v2ang

  ! Returns the absolute magnitude of a double precision vector, v, of
  ! arbitrary length.
  function gn_absv(v) result(absv)
    real(dp), intent(in) :: v(:)
    real(dp) :: absv

    absv = sqrt(dot_product(v,v))
  end function gn_absv

  ! Creates the double precision 3 x 3 rotation matrix that rotates a
  ! 3-dimensional vector through Euler angles alpha, beta and gamma.
  function gn_rotm_euler(alpha, beta, gamma) result(rotm)
    real(dp), intent(in) :: alpha, beta, gamma
    real(dp) :: rotm(3, 3)

    real(dp) :: cosalpha, cosbeta, cosgamma, sinalpha, sinbeta, singamma

    sinalpha = sin(alpha)
    sinbeta = sin(beta)
    singamma = sin(gamma)
    cosalpha = cos(alpha)
    cosbeta = cos(beta)
    cosgamma = cos(gamma)

    rotm(1, 1) = cosalpha * cosbeta * cosgamma - sinalpha * singamma
    rotm(1, 2) = - cosalpha * cosbeta * singamma - sinalpha * cosgamma
    rotm(1, 3) = cosalpha * sinbeta

    rotm(2, 1) = sinalpha * cosbeta * cosgamma + cosalpha * singamma
    rotm(2, 2) = - sinalpha * cosbeta * singamma + cosalpha * cosgamma
    rotm(2, 3) = sinalpha * sinbeta

    rotm(3, 1) = - sinbeta * cosgamma
    rotm(3, 2) = sinbeta * singamma
    rotm(3, 3) = cosbeta
  end function gn_rotm_euler

  ! Creates the double precision 3 x 3 rotation matrix that rotates by
  ! angle phi (in a right-handed sense) around the specified axis (1 = x,
  ! 2 = y, 3 = z).
  function gn_rotm_axis(phi, axis) result(rotm)
    real(dp), intent(in) :: phi
    integer, intent(in) :: axis
    real(dp) :: rotm(3, 3)

    real(dp) :: sinphi, cosphi

    sinphi = sin(phi)
    cosphi = cos(phi)

    rotm = 0.0_dp

    if (axis == 1) then
      rotm(1, 1) = 1.0
      rotm(2, 2) = cosphi
      rotm(2, 3) = - sinphi
      rotm(3, 2) = sinphi
      rotm(3, 3) = cosphi
    else if (axis == 2) then
      rotm(1, 1) = cosphi
      rotm(1, 3) = sinphi
      rotm(2, 2) = 1.0
      rotm(3, 1) = - sinphi
      rotm(3, 3) = cosphi
    else if (axis == 3) then
      rotm(1, 1) = cosphi
      rotm(1, 2) = - sinphi
      rotm(2, 1) = sinphi
      rotm(2, 2) = cosphi
      rotm(3, 3) = 1.0
    else
      call gn_fatal('gn_rotm_axis: axis incorrect', axis)
    end if
  end function gn_rotm_axis

  ! Returns the modified input 3 x 3 rotation matrix, rotm_in, by
  ! premultiplying it with the 3 x 3 rotation matirx corresponding to
  ! a rotation of angle phi (in a right-handed sense) around the specified
  ! axis (1 = x, 2 = y, 3 = z). (In other words this rotation is physically
  ! applied to a vector after those encoded in rotm_in.)
  function gn_rotm_multaxis(phi, axis, rotm_in) result(rotm)
    real(dp), intent(in) :: phi
    integer, intent(in) :: axis
    real(dp), intent(in) :: rotm_in(3, 3)
    real(dp) :: rotm(3, 3)

    ! Calculate the rotation matrix corresponding to the new rotation.
    rotm = gn_rotm(phi, axis)

    ! Postmultiply it by the old rotation.
    rotm = matmul(rotm, rotm_in)
  end function gn_rotm_multaxis

  ! Calculates the Euler angles encoded in a double precision 3 x 3 rotation
  ! matrix, rotm. In the case of ambiguity (e.g., when beta = 0.0) it is
  ! assumed that gamma = 0.0 and all the azimuthal rotation is encoded in
  ! alpha. Note that it is assumed that the input matrix is a rotation
  ! matrix and hence that the elements have values in the range - 1.0 to
  ! 1.0, facillitating the use of inverse trignonometric functions; if
  ! an arbitrary matrix is input anything could happen.
  subroutine gn_rotm2ang(rotm, alpha, beta, gamma)
    real(dp), intent(in) :: rotm(3, 3)
    real(dp), intent(out) :: alpha, beta, gamma

    real(dp) :: cosbeta, sinbeta, cosalpha, sinalpha, cosgamma, singamma

    ! First get beta unambiguously from rotm(3, 3) = cos(beta) and
    ! the fact that 0.0 <= beta <= pi (unlike the other two Euler
    ! angles which run from 0.0 <= {alpha, gamma} < 2 pi.
    cosbeta = rotm(3, 3)
    beta = acos(cosbeta)
    sinbeta = sin(beta)

    if (abs(sinbeta) <= 1.0e-6) then

      ! If sin(beta) = 0.0 then the aziumthal angles are ambiguous.
      if (cosbeta > 0.0) then
        ! With beta = 0.0 there is no polar rotation at all, just an
        ! azimuthal rotation about the z-axis which, if alpha = 0.0 is
        ! assumed, gives gamma unambiguously.
        alpha = 0.0
        cosgamma = rotm(1, 1)
        singamma = rotm(2, 1)
        gamma = gn_atan(cosgamma, singamma)
      else
        ! If beta = pi then the sky is flipped, top to bottom, between
        ! the two azimuthal rotations and so it is the difference
        ! alpha - gamma that is constrained, which again gives gamma if
        ! alpha = 0.0 is assumed.
        alpha = 0.0
        cosgamma = - rotm(1, 1)
        singamma = rotm(1, 2)
        gamma = gn_atan(cosgamma, singamma)
      end if

    else
      ! In the more general case the elements of the rotation matrix
      ! can be ``inverted'' to give the Euler angles. But the inversion
      ! requires two elements (sin and cos) brought together to obtain
      ! the values for the azimuthal angles. These give
      !
      !   cos(alpha) = rotm(1, 3) / sin(beta)
      !   sin(alpha) = rotm(2, 3) / sin(beta)
      !
      ! and
      !
      !   cos(gamma) = - rotm(3, 1) / sin(beta)
      !   sin(gamma) = rotm(3, 2) / sin(beta)
      !
      ! but hence can be combined to give alpha and gamma without including
      ! the sin(beta) terms. It is also possible that these expressions
      ! might become innacurate for certain combinations angles; this has
      ! not been investigated.
      cosalpha = rotm(1, 3)
      sinalpha = rotm(2, 3)
      alpha = gn_atan(cosalpha, sinalpha)

      cosgamma = - rotm(3, 1)
      singamma = rotm(3, 2)
      gamma = gn_atan(cosgamma, singamma)
    end if
  end subroutine gn_rotm2ang

  ! Calculate a single precision rotation matrix from the Tait-Bryan angles.
  function gn_rotm_tb(tb) result(rotm)
    real(dp), intent(in) :: tb(3)
    real(dp) :: rotm(3, 3)

    ! The Tait-Bryan angles correspond to consecutive rotations about the
    ! three priniciple axes.
    rotm = gn_rotm(tb(1), 1)
    rotm = gn_rotm(tb(2), 2, rotm)
    rotm = gn_rotm(tb(3), 3, rotm)
  end function gn_rotm_tb

  ! Return the spherical angles of the rotated coordinate axes from a
  ! rotation matrix.
  function gn_rotm2axes(rotm) result(axes)
    real(dp), intent(in) :: rotm(3, 3)
    type(gnsphangle_double) :: axes(3)

    axes(1) = gn_v2ang(rotm(:, 1))
    axes(2) = gn_v2ang(rotm(:, 2))
    axes(3) = gn_v2ang(rotm(:, 3))
  end function gn_rotm2axes

  ! Returns a spherical angle in double precision.
  function gn_sphang(theta, phi) result(ang)
    real(dp), intent(in) :: theta, phi
    type(gnsphangle_double) :: ang

    ang%theta = theta
    ang%phi = phi
    ang = gn_ang2ang(ang)
  end function gn_sphang

end module general_vector
