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

module general_rand
  use planck_config
  use general_error
  use ls_paramfile_io
  implicit none
  private

  public :: gn_rand_init, gn_monte_rand, gn_monte_gauss, gn_monte_gauss_clip

  ! Initialises the random number generator in the specified way.
  interface gn_rand_init
    module procedure gn_rand_init_manual, gn_rand_init_par
  end interface

  ! Returns a random number drawn from a Gaussian distribution.
  interface gn_monte_gauss
    module procedure gn_monte_gauss_0, gn_monte_gauss_mean
  end interface

  ! Returns a random number drawn from a clipped Gaussian distribution.
  interface gn_monte_gauss_clip
    module procedure gn_monte_gauss_clip_0, gn_monte_gauss_clip_mean
  end interface

  ! The random number generator used here makes extensive use of both
  ! global and saved variables. These are declared here, but are in
  ! all cases private and cannot be accessed from outside this module.
  real(dp), save :: &
    gnc = 362436.0 / 16777216.0, &
    gncd = 7654321.0 / 16777216.0, &
    gncm = 16777213.0 / 16777216.0
  integer, save :: &
    gni97 = 97, &
    gnj97 = 33

  integer, parameter, public :: &
    GNIJ_MAX = 31328, &
    GNKL_MAX = 30081

  integer, save :: &
    ij_save = 0, &
    kl_save = 0

  logical, save :: &
    gnrand_init = .false., &
    gnrand_set = .false.
  real(dp), save :: gnu(97)
  real(dp), save :: gnmonte_save

contains

  ! Seed the random number generator by reading in the two seeds via
  ! parameter parsing (but leave the default values as the clock seeds).
  subroutine gn_rand_init_par(params)
    type(paramfile_handle), intent(inout) :: params
    integer :: ij, kl

    ! Read in the two seeds, with default values set from the system
    ! clock.  Note that the clock calls may have to be modified to
    ! ensure that they are independent in the case that the return
    ! value of the clock function is below both GNIJ_MAX and GNKL_MAX
    ! (i.e., the situation in which the two would have the same value).
    ij = parse_int(params,'randseed_1', vmin=0, vmax=GNIJ_MAX)
    kl = parse_int(params,'randseed_2', vmin=0, vmax=GNKL_MAX)

    ! Call the manual seeding subroutine.
    call gn_rand_init(ij, kl)
  end subroutine gn_rand_init_par

  ! Seed the random number generator manually and set the gn_rand_init
  ! flag to true so the program knows that this call has been made.
  subroutine gn_rand_init_manual(ij, kl)
    integer, intent(in) :: ij, kl

    integer :: i, j, k, l, m, ii, jj
    real(dp) :: s, t

    ! Check they are within acceptable limits.
    call gn_assert (ij>=0,'gn_rand_init_manual: ij < 0: seed out of range', ij)
    call gn_assert (ij<=GNIJ_MAX, &
      'gn_rand_init_manual: ij > GNIJ_MAX: seed out of range', ij)
    call gn_assert (kl>=0,'gn_rand_init_manual: kl < 0: seed out of range', kl)
    call gn_assert (kl<=GNKL_MAX, &
      'gn_rand_init_manual: kl > GNKL_MAX: seed out of range', kl)

    ! Perform initialisation calculations.
    i = mod(ij / 177, 177) + 2
    j = mod(ij , 177) + 2
    k = mod(kl / 169, 178) + 1
    l = mod(kl, 169)

    ii = 1
    do
      if (ii > 97) then
        exit
      end if

        s = 0.0
        t = 0.5

        jj = 1
        do
          if (jj > 24) exit

          m = mod(mod(i * j, 179) * k, 179)
          i = j
          j = k
          k = m
          l = mod(53 * l + 1, 169)
          if (mod(l * m, 64) >= 32) then
             s = s + t
          end if
          t = 0.5 * t

          jj = jj + 1
        end do

        gnu(ii) = s

      ii = ii + 1
    end do

    ! Copy the values of ij and kl to prviate module values that can be
    ! accessed using the gn_rand_seeds routine for immortalisation in a
    ! data file.
    ij_save = ij
    kl_save = kl

    ! Set gnrand_init to true so that calls to gn_monte_rand know that the
    ! initialisation has been called.
    gnrand_init = .true.
  end subroutine gn_rand_init_manual

  ! Returns a uniform deviate random number between min and max.
  function gn_monte_rand(min, max) result(monte)
    real(dp), intent(in) :: min, max
    real(dp) :: monte

    real(dp) :: monte_temp

    call gn_assert (gnrand_init, &
      'gn_monte_rand: gn_rand_init must be called first')
    call gn_assert(min<max,'gn_monte_rand: min >= max')

    monte_temp = gnu(gni97) - gnu(gnj97)
    if (monte_temp < 0.0) monte_temp = monte_temp + 1.0

    gnu(gni97) = monte_temp
    gni97 = gni97 - 1
    if (gni97 == 0) gni97 = 97

    gnj97 = gnj97 - 1
    if (gnj97 == 0) gnj97 = 97

    gnc = gnc - gncd
    if (gnc < 0.0) gnc = gnc + gncm

    monte_temp = monte_temp - gnc
    if (monte_temp < 0.0) monte_temp = monte_temp + 1.0

    ! Scale monte to the required range, and convert to GNDP, in order
    ! to interface with the rest of the program.
    monte = min + monte_temp * (max - min)
  end function gn_monte_rand

  ! Returns a random number generated from a Gaussian distribution with
  ! expectation value mean and standard deviation sigma. This is generated
  ! using the Box-Muller method.
  function gn_monte_gauss_mean(mean, sigma) result(monte)
    real(dp), intent(in) :: mean, sigma
    real(dp) :: monte

    real(dp) :: fac, rsq, v1, v2

    call gn_assert(sigma>=0, 'gn_monte_gauss_mean: sigma < 0.0', sigma)
    if (sigma == 0.0) then

      monte = mean

    else

      if (.not. gnrand_set) then

        do
          v1 = gn_monte_rand(- 1.0_dp, 1.0_dp)
          v2 = gn_monte_rand(- 1.0_dp, 1.0_dp)
          rsq = v1**2 + v2**2

          if ((rsq < 1.0) .and. (rsq /= 0.0)) exit
        end do

        fac = sqrt(- 2.0 * log(rsq) / rsq)
        monte = v1 * fac
        gnmonte_save = v2 * fac
        gnrand_set = .true.

      else

        monte = gnmonte_save
        gnrand_set = .false.

      end if

    end if

    monte = mean + sigma * monte
  end function gn_monte_gauss_mean

  ! Returns a random number generated from a Gaussian distribution with
  ! zero mean and standard deviation sigma; this is just a wrapper to the
  ! function gn_monte_gauss_x.
  function gn_monte_gauss_0(sigma) result(monte)
    real(dp), intent(in) :: sigma
    real(dp) :: monte

    monte = gn_monte_gauss(0.0_dp, sigma)
  end function gn_monte_gauss_0

  ! Returns a random number generated from a Gaussian distribution with
  ! expectation value mean and standard deviation sigma, but with the
  ! distribution ``clipped'' so that deviations of more than clip from
  ! mean are rejected.
  function gn_monte_gauss_clip_mean(mean, sigma, clip) result(monte)
    real(dp), intent(in) :: mean, sigma, clip
    real(dp) :: monte

    call gn_assert (clip>=0.0,'gn_monte_gauss_clip_mean: clip < 0.0', clip)
    if (clip == 0.0) then
      monte = mean
    else
      ! Call the Gaussian random number generator repeatedly, rejecting
      ! each value that is too far from the mean.
      do
        monte = gn_monte_gauss(mean, sigma)
        if (abs(monte - mean) <= clip) exit
      end do
    end if
  end function gn_monte_gauss_clip_mean

  ! Returns a random number generated from a Gaussian distribution with
  ! zero mean and standard deviation sigma, but with the distribution
  ! ``clipped'' to have absolute value less than clip.
  function gn_monte_gauss_clip_0(sigma, clip) result(monte)
    real(dp), intent(in) :: sigma, clip
    real(dp) :: monte

    monte = gn_monte_gauss_clip(0.0_dp, sigma, clip)
  end function gn_monte_gauss_clip_0

end module general_rand
