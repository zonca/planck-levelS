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

module beam_transform

  ! This module contains routines to perform the spherical harmonic
  ! analysis of a beam on a polar grid.
  !
  ! Mark Ashdown, CPAC.

  use planck_config
  use sharpf_mod
  use beam_polar
  use beam_alm

  implicit none

contains

  !======================================================================

  ! Subroutine to extract the multipoles from the Stokes parameters of
  ! beam on a polar grid by quadrature.  The input bmpolar type can be
  ! a full-sky beam or just a main beam. In the latter case it is
  ! assumed that the beam is zero elsewhere on the sphere.  Note that
  ! this subroutine does not set the multipoles in the bmalm type to
  ! zero before doing the transform; instead it adds the extracted
  ! multipoles to the input values.

  subroutine bm_polar2alm(beam, alm)
    type(bmpolar), intent(in) :: beam
    type(bmalm), intent(inout) :: alm

    integer :: i, s_phi
    real(dp) :: dtheta
    integer, dimension(beam%ntheta) :: nph, thetaofs, phistride
    integer, dimension(alm%mmax+1) :: alm_m_start
    real(dp), dimension(beam%ntheta) :: phi0, theta, weight

    dtheta = (beam%theta_max-beam%theta_min)/(beam%ntheta-1)
    s_phi = size(beam%stokes,1)
    call assert(s_phi==4, 'unexpected number of Stokes parameters in beam')

    if ((2*alm%lmax+1)>(pi/dtheta + 10)) &
      print *,"WARNING: lmax too high for this grid spacing"
    if ((2*alm%mmax+1)>beam%nphi) &
      print *,"WARNING: mmax too high for this grid spacing"

    do i = 1, beam%ntheta
      nph(i) = beam%nphi
      thetaofs(i) = (i-1)*s_phi*beam%nphi
      phistride(i) = s_phi
      phi0(i) = 0.0_dp
      theta(i) = beam%theta_min + (i-1.0_dp)*dtheta
      if (theta(i)==0.0_dp) theta(i) = 0.001_dp*dtheta
      weight(i) = twopi*sin(theta(i))*dtheta/beam%nphi
    end do

    do i = 1, alm%mmax+1
      alm_m_start(i) = (i-1)*(alm%lmax+1)*alm%nstokes
    end do

    select case (alm%nstokes)
    case (1)
      call sharpf_do_job_d(beam%ntheta, nph, thetaofs, phistride, phi0, &
          theta, weight, beam%stokes, 1, alm%lmax, alm%mmax, alm%a, 1, &
          alm_m_start, alm%nstokes, MAP2ALM, 1)
    case (3)
      call sharpf_do_job_d(beam%ntheta, nph, thetaofs, phistride, phi0, &
          theta, weight, beam%stokes, 1, alm%lmax, alm%mmax, alm%a, 1, &
          alm_m_start, alm%nstokes, MAP2ALM_POL, 1)
    case default
      call exit_with_status(1, 'incorrect nstokes in bm_polar2alm')
    end select

  end subroutine bm_polar2alm

  !======================================================================

end module beam_transform
