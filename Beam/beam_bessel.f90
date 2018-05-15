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

module beam_bessel

  ! Bessel functions for calculating elliptical beam multipoles.

  use planck_config

  implicit none

  private :: bessi0, bessi1, bessi, poly

contains

  !======================================================================

  ! Modified Bessel function.

  function bessel_i(n, x)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: bessel_i

    select case (n)
    case(0)
      bessel_i = bessi0(x)
    case(1)
      bessel_i = bessi1(x)
    case(2:)
      bessel_i = bessi(n, x)
    case default
      bessel_i = 0
      call exit_with_status (1, 'bad order for Bessel function')
    end select

  end function bessel_i

  !======================================================================

  function bessi0(x)
    real(dp), intent(in) :: x
    real(dp) :: bessi0

    real(dp) :: ax
    real(dp), dimension(7) :: p = (/ 1.0_dp, 3.5156229_dp, &
        3.0899424_dp, 1.2067492_dp, 0.2659732_dp, 0.360768e-1_dp, &
        0.45813e-2_dp /)
    real(dp), dimension(9) :: q = (/ 0.39894228_dp, 0.1328592e-1_dp, &
        0.225319e-2_dp, -0.157565e-2_dp, 0.916281e-2_dp, &
        -0.2057706e-1_dp, 0.2635537e-1_dp, -0.1647633e-1_dp, &
        0.392377e-2_dp /)

    ax = abs(x)
    if (ax < 3.75_dp) then
      bessi0 = poly(real((x/3.75_dp)**2, dp), p)
    else
      bessi0 = (exp(ax)/sqrt(ax))*poly(real(3.75_dp/ax, dp), q)
    end if

  end function bessi0

  !======================================================================

  function bessi1(x)
    real(dp), intent(in) :: x
    real(dp) :: bessi1

    real(dp) :: ax
    real(dp), dimension(7) :: p = (/ 0.5_dp, 0.87890594_dp, &
        0.51498869_dp, 0.15084934_dp, 0.2658733e-1_dp, &
        0.301532e-2_dp, 0.32411e-3_dp /)
    real(dp), dimension(9) :: q = (/ 0.39894228_dp, -0.3988024e-1_dp, &
        -0.362018e-2_dp, 0.163801e-2_dp, -0.1031555e-1_dp, &
        0.2282967e-1_dp, -0.2895312e-1_dp, 0.1787654e-1_dp, &
        -0.420059e-2_dp /)

    ax = abs(x)
    if (ax < 3.75_dp) then
      bessi1 = ax * poly(real((x/3.75_dp)**2, dp), p)
    else
      bessi1 = (exp(ax)/sqrt(ax))*poly(real(3.75_dp/ax, dp), q)
    end if
    if (x < 0.0) bessi1 = -bessi1

  end function bessi1

  !======================================================================

  function bessi(n, x)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: bessi

    integer, parameter :: iacc = 40, iexp = maxexponent(x)/2
    integer :: j, m
    real(dp) :: bi, bim, bip, tox
    bessi = 0.0_dp
    if (x*x <= 8.0_dp*tiny(x)) return
    tox = 2.0_dp/abs(x)
    bip = 0.0_dp
    bi = 1.0_dp
    m = 2*(n+int(sqrt(real(iacc*n, dp))))
    do j = m, 1, -1
      bim = bip + j*tox*bi
      bip = bi
      bi = bim
      if (exponent(bi) > iexp) then
        bessi = scale(bessi, -iexp)
        bi = scale(bi, -iexp)
        bip = scale(bip, -iexp)
      end if
      if (j == n) bessi = bip
    end do
    bessi = bessi*bessi0(x)/bi
    if (x < 0.0_dp .and. mod(n,2) == 1) bessi = -bessi

  end function bessi

  !======================================================================

  ! Evaluates polynomial given argument and coefficients

  function poly(x, coeffs)
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: coeffs
    real(dp) :: poly

    integer :: i, n

    n = size(coeffs)
    if (n <= 0) then
      poly = 0.0_dp
    else
      poly = coeffs(n)
      do i = n-1, 1, -1
        poly = x*poly + coeffs(i)
      end do
    end if

  end function poly

  !======================================================================

end module beam_bessel
