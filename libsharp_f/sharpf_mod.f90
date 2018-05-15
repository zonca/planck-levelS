!-----------------------------------------------------------------------------
!
!  This file is part of the "libsharp_f" component of the Planck simulation
!  package.
!
!  libsharp_f is free software; you can redistribute it and/or modify
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
!  along with libsharp_f; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!-----------------------------------------------------------------------------

!  Copyright (C) 2012-2013 Max-Planck-Society
!  \author Martin Reinecke

module sharpf_mod
use planck_config
implicit none
private

integer, parameter, public:: MAP2ALM=1, MAP2ALM_POL=2

public sharpf_do_job_s, sharpf_do_job_d

interface

subroutine sharps_do_job_f (ntheta, nphi, theta_start, phi_stride, phi_0, &
  theta, weight, map, map_stride, lmax, mmax, alm, alm_stride, alm_m_start, &
  alm_l_stride, jobtype, add) bind(c,name="sharps_do_job_f")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), value  :: ntheta, map_stride, lmax, mmax, alm_stride, &
    alm_l_stride, jobtype, add
  integer(c_int), dimension(*), intent(in) :: nphi, theta_start, phi_stride, &
    alm_m_start
  real(c_double), dimension(*), intent(in) :: phi_0, theta, weight
  real(c_float), dimension(*) :: map
  complex(c_float_complex), dimension(*) :: alm
end subroutine

subroutine sharpd_do_job_f (ntheta, nphi, theta_start, phi_stride, phi_0, &
  theta, weight, map, map_stride, lmax, mmax, alm, alm_stride, alm_m_start, &
  alm_l_stride, jobtype, add) bind(c,name="sharpd_do_job_f")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), value  :: ntheta, map_stride, lmax, mmax, alm_stride, &
    alm_l_stride, jobtype, add
  integer(c_int), dimension(*), intent(in) :: nphi, theta_start, phi_stride, &
    alm_m_start
  real(c_double), dimension(*), intent(in) :: phi_0, theta, weight
  real(c_double), dimension(*) :: map
  complex(c_double_complex), dimension(*) :: alm
end subroutine

end interface

contains

! Parameter definition
! ntheta: number of thetas (or rings) in this pixelisation
! nphi(ntheta): number of pixels for each individual theta
! theta_start(ntheta): index of the first pixel of each ring minus
!                      index of the first pixel of the first ring
! phi_stride: index difference in the map array between two ring pixels
! phi_0(ntheta): longitude(rad) of the first pixel in each ring
! theta(ntheta): colatitude(rad) of each ring
! weight(ntheta): weight of an individual _pixel_ in each ring
!                 i.e. area of a single ring pixel / 4*pi
! map(*): array containing all pixels of all maps
! map_stride: index difference in the map array between the first pixel of
!             the second map and the first pixel of the first map
! alm(*): array containing all coefficients of all a_lm components
! alm_stride: index difference in the alm array between the first coeff of
!             the second a_lm set and the first coeff of the second a_lm set

! Example: The index the third a_lm element in the a_lm array is (zero-based):
! i = 2*alm_stride + alm_m_start(m) + l*alm_l_stride

subroutine sharpf_do_job_s (ntheta, nphi, theta_start, phi_stride, phi_0, &
  theta, weight, map, map_stride, lmax, mmax, alm, alm_stride, alm_m_start, &
  alm_l_stride, jobtype, add)
  integer, intent(in) :: ntheta, map_stride, lmax, mmax, alm_stride, &
    alm_l_stride, jobtype, add
  integer, dimension(:), intent(in) :: nphi, theta_start, phi_stride, &
    alm_m_start
  real(dp), dimension(:), intent(in) :: phi_0, theta, weight
  real(sp), dimension(*) :: map
  complex(sp), dimension(*) :: alm

  call sharps_do_job_f (ntheta, nphi, theta_start, phi_stride, phi_0, &
    theta, weight, map, map_stride, lmax, mmax, alm, alm_stride, alm_m_start, &
    alm_l_stride, jobtype, add)
end subroutine sharpf_do_job_s

subroutine sharpf_do_job_d (ntheta, nphi, theta_start, phi_stride, phi_0, &
  theta, weight, map, map_stride, lmax, mmax, alm, alm_stride, alm_m_start, &
  alm_l_stride, jobtype, add)
  integer, intent(in) :: ntheta, map_stride, lmax, mmax, alm_stride, &
    alm_l_stride, jobtype, add
  integer, dimension(:), intent(in) :: nphi, theta_start, phi_stride, &
    alm_m_start
  real(dp), dimension(:), intent(in) :: phi_0, theta, weight
  real(dp), dimension(*) :: map
  complex(dp), dimension(*) :: alm

  call sharpd_do_job_f (ntheta, nphi, theta_start, phi_stride, phi_0, &
    theta, weight, map, map_stride, lmax, mmax, alm, alm_stride, alm_m_start, &
    alm_l_stride, jobtype, add)
end subroutine sharpf_do_job_d

end module sharpf_mod
