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

module beam_interpolate

  ! Beam interpolation module. This contains routines for
  ! interpolating beams from square grids to polar grids and
  ! interpolating the values of the beams on a polar grid in the theta
  ! direction.
  !
  ! Mark Ashdown, CPAC

  use planck_config
  use ls_misc_utils

  use beam_polar
  use beam_square

  implicit none

  private :: bm_spline, bm_splint, bm_locate, bm_tridiag

contains

  !======================================================================

  ! Interpolates beam from square grid to polar grid.

  subroutine bm_square2polar(beam_s, beam_p, nphi, ntheta)
    type(bmsquare), intent(in) :: beam_s
    type(bmpolar), intent(inout) :: beam_p
    integer, intent(in) :: nphi, ntheta

    integer :: i, j, itheta, iphi
    real(dp) :: theta_min, theta_max
    real(dp) :: theta, phi, dtheta, dphi, sinth
    real(dp) :: s, t, x, y, fx, fy

    ! Note theta min and max angles are in radians.

    theta_min = 0.0_dp
    theta_max = asin(beam_s%xdelta * real((beam_s%nx-1)/2, dp))

    call bm_polar_init(beam_p, nphi, ntheta, theta_min, theta_max)

    ! Interpolate Stokes parameters onto polar grid using bilinear
    ! interpolation.

    dtheta = theta_max / real(ntheta-1, dp)
    dphi = twopi / real(nphi, dp)

    do itheta = 1, ntheta

      theta = real(itheta-1, dp) * dtheta
      sinth = sin(theta)

      do iphi = 1, nphi

        phi = real(iphi-1, dp) * dphi
        x = sinth * cos(phi)
        y = sinth * sin(phi)

        fx = x / beam_s%xdelta
        fy = y / beam_s%ydelta
        i = floor(fx)
        j = floor(fy)
        s = fx - real(i, dp)
        t = fy - real(j, dp)

        ! Add correct offset to centre of grid.

        i = i + beam_s%xcentre
        j = j + beam_s%ycentre
!FIXME
if ((i<0).or.(i>beam_s%nx).or.(j<0).or.(j>beam_s%ny)) &
  print *,"Problem: ",i,j
        if (i >= beam_s%nx) then
          beam_p%stokes(:, iphi, itheta) = &
              (1.0_dp-t)*beam_s%stokes(:, beam_s%nx, j) + &
              t*beam_s%stokes(:, beam_s%nx, j+1)
        else if (i <= 0) then
          beam_p%stokes(:, iphi, itheta) = &
              (1.0_dp-t)*beam_s%stokes(:, 1, j) + &
              t*beam_s%stokes(:, 1, j+1)
        else if (j >= beam_s%ny) then
          beam_p%stokes(:, iphi, itheta) = &
              (1.0_dp-s)*beam_s%stokes(:, i, beam_s%ny) + &
              s*beam_s%stokes(:, i+1, beam_s%ny)
        else if (j <= 0) then
          beam_p%stokes(:, iphi, itheta) = &
              (1.0_dp-s)*beam_s%stokes(:, i, 1) + &
              s*beam_s%stokes(:, i+1, 1)
        else
          beam_p%stokes(:, iphi, itheta) = &
              (1.0_dp-s)*(1.0_dp-t)*beam_s%stokes(:, i, j) + &
              s*(1.0_dp-t)*beam_s%stokes(:, i+1, j) + &
              (1.0_dp-s)*t*beam_s%stokes(:, i, j+1) + &
              s*t*beam_s%stokes(:, i+1, j+1)
        end if

      end do
    end do

  end subroutine bm_square2polar

  !======================================================================

  ! Interpolates beam on polar grid in the theta direction using
  ! linear interpolation.

  subroutine bm_interp_linear(beam, tminnew, tmaxnew, ntnew)
    type(bmpolar), intent(inout) :: beam
    real(dp), intent(in) :: tminnew, tmaxnew
    integer, intent(in) :: ntnew

    real(dp), dimension(:,:,:), pointer :: valnew, valold
    real(dp), dimension(:), allocatable :: tnew, told
    real(dp) :: t, dtheta, v1, v2, t1, t2
    real(dp) :: tminold, tmaxold
    integer :: itheta, iphi, icmpt, i
    integer :: ntold

    ntold = beam%ntheta
    tminold = beam%theta_min
    tmaxold = beam%theta_max

    allocate(told(ntold)) ! Old theta values
    allocate(tnew(ntnew)) ! New theta values
    allocate(valnew(4, beam%nphi, ntnew)) ! New beam values
    valold => beam%stokes

    ! Calculate old theta values

    dtheta = (tmaxold-tminold)/real(ntold-1, kind=dp)
    do itheta = 1, ntold
      told(itheta) = tminold + real(itheta-1, kind=dp)*dtheta
    end do

    ! Calculate new theta values

    dtheta =  (tmaxnew-tminnew)/real(ntnew-1, kind=dp)
    do itheta = 1, ntnew
      tnew(itheta) = tminnew + real(itheta-1, kind=dp)*dtheta
    end do

    ! Linear interpolation

    do iphi = 1, beam%nphi
      do icmpt = 1, 4
        do itheta = 1, ntnew
          t = tnew(itheta)
          i = bm_locate(told, t)
          if (i == 0) then
            valnew(icmpt, iphi, itheta) = valold(icmpt, iphi, 1)
          else if (i == ntold) then
            valnew(icmpt, iphi, itheta) = valold(icmpt, iphi, ntold)
          else
            v1 = valold(icmpt, iphi, i)
            v2 = valold(icmpt, iphi, i+1)
            t1 = told(i)
            t2 = told(i+1)
            valnew(icmpt, iphi, itheta) = (v1*(t2-t)+v2*(t-t1))/(t2-t1)
          end if
        end do
      end do
    end do

    ! Put the new values into beam

    beam%theta_min = tminnew
    beam%theta_max = tmaxnew
    beam%ntheta = ntnew
    beam%stokes => valnew

    ! Deallocate the old pointers and the temporary arrays

    deallocate(tnew, told, valold)

  end subroutine bm_interp_linear

  !======================================================================

  ! Interpolates beam on polar grid in the theta direction using
  ! cubic spline interpolation.

  subroutine bm_interp_spline(beam, tminnew, tmaxnew, ntnew)
    type(bmpolar), intent(inout) :: beam
    real(dp), intent(in) :: tminnew, tmaxnew
    integer, intent(in) :: ntnew

    real(dp), dimension(:,:,:), pointer :: valnew, valold
    real(dp), dimension(:), allocatable :: tnew, told
    real(dp), dimension(:), allocatable :: v2, valtemp
    real(dp) :: theta, dtheta
    real(dp) :: tminold, tmaxold
    integer :: itheta, iphi, icmpt
    integer :: ntold

    ntold = beam%ntheta
    tminold = beam%theta_min
    tmaxold = beam%theta_max

    allocate(v2(ntold))      ! Second derivatives of beam values
    allocate(valtemp(ntold)) ! Temporary storage of beam values
    allocate(told(ntold))    ! Old theta values
    allocate(tnew(ntnew))    ! New theta values
    allocate(valnew(4, beam%nphi, ntnew)) ! New beam values
    valold => beam%stokes

    ! Calculate old theta values

    dtheta = (tmaxold-tminold)/real(ntold-1, kind=dp)
    do itheta = 1, ntold
      told(itheta) = tminold + real(itheta-1, kind=dp)*dtheta
    end do

    ! Calculate new theta values

    dtheta = (tmaxnew-tminnew)/real(ntnew-1, kind=dp)
    do itheta = 1, ntnew
      tnew(itheta) = tminnew + real(itheta-1, kind=dp)*dtheta
    end do

    ! Do the interpolation with spline

    do iphi = 1, beam%nphi
      do icmpt = 1, 4
        valtemp =  valold(icmpt, iphi, :)
        call bm_spline(told, valtemp, 0.0_dp, 0.0_dp, v2)
        do itheta = 1, ntnew
          theta = tnew(itheta)
          valnew(icmpt, iphi, itheta) = bm_splint(told, valtemp, v2, theta)
        end do
      end do
    end do

    ! Put the new values into beam

    beam%theta_min = tminnew
    beam%theta_max = tmaxnew
    beam%ntheta = ntnew
    beam%stokes => valnew

    ! Deallocate the old pointers and the temporary arrays

    deallocate(tnew, told, valold, valtemp, v2)

  end subroutine bm_interp_spline

  !======================================================================

  ! Given arrays x and y of length n containing a tabulated function,
  ! i.e. y(i)=f(x(i)), with x(1) < x(2) < ... <x(n), and given values
  ! yp1 and ypn for the first derivative of the interpolating function
  ! at points 1 and n, respectively, this routine returns an array y2 of
  ! length n that contains the second derivatives of the interpolating
  ! function at the tabulated points x(i). If yp1 or ypn is equal to
  ! 10e30 or larger, the routine is signaled to set the corresponding
  ! boundary condition for a natural spline, with zero second derivative
  ! on that boundary.

  subroutine bm_spline(x, y, yp1, ypn, y2)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), intent(in) :: yp1, ypn
    real(dp), dimension(:), intent(out) :: y2

    integer :: n
    real(dp), dimension(size(x)) :: a, b, c, r

    n = size(x)

    ! Set up the tridiagonal equations.
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0

    if (yp1 > 0.99e30_dp) then
      ! The lower boundary condition is set either to be natural
      r(1)=0.0
      c(1)=0.0
    else
      ! or else to have the specified first derivative.
      r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      c(1)=0.5
    end if
    if (ypn > 0.99e30_dp) then
      ! The upper boundary condition is set either to be natural
      r(n)=0.0
      a(n)=0.0
    else
      !or else to have the specified first derivative.
      r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
      a(n)=0.5
    end if

    call bm_tridiag(a(2:n), b(1:n), c(1:n-1), r(1:n), y2(1:n))

  end subroutine bm_spline

  !======================================================================

  ! Given the arrays xa and ya, which tabulate a function (with the
  ! elements of xa in increasing or decreasing order),and given the
  ! array y2a which is the output from bm_spline above, and given
  ! a value of x, this routine returns a cubic-spline interpolated
  ! value. The arrays xa, ya and y2a are of the same size.

  function bm_splint(xa, ya, y2a, x)
    real(dp), dimension(:), intent(in) :: xa, ya, y2a
    real(dp), intent(in) :: x
    real(dp):: bm_splint

    integer :: khi,klo,n
    real(dp) :: a,b,h

    n = size(xa)

    ! Find the values klo and khi which bracket the input value of x
    ! using bm_locate.

    klo = max(min(bm_locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    ! Check whether the xa's are distinct
    call assert(h/=0.0,'bad xa input in bm_splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    ! Cubic spline polynomial is now evaluated.
    bm_splint=a*ya(klo)+b*ya(khi)+ &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
  end function bm_splint

  !======================================================================

  function bm_locate(xx, x)
    real(dp), dimension(:), intent(in) :: xx
    real(dp), intent(in) :: x

    integer :: bm_locate
    integer :: n, jl, jm, ju
    logical :: ascnd

    n = size(xx)
    ascnd = (xx(n)>=xx(1))
    jl=0
    ju=n+1
    do
      if (ju-jl <= 1) exit
      jm=(ju+jl)/2
      if (ascnd .eqv. (x >=xx(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do
    if (x == xx(1)) then
      bm_locate=1
    else if (x == xx(n)) then
      bm_locate=n-1
    else
      bm_locate=jl
    end if
  end function bm_locate

  !======================================================================

  subroutine bm_tridiag(a,b,c,r,u)
    real(dp), dimension(:), intent(in) :: a, b, c, r
    real(dp), dimension(:), intent(out) :: u

    real(dp), dimension(size(b)) :: gam
    integer :: n, j
    real(dp) :: bet

    n = size(b)
    bet=b(1)

    ! If the following happens then you should rewrite your equations as a
    ! set of order n-1 ,with u(2) trivially eliminated.
    call assert (bet/=0.0_dp, 'bm_tridiag: Error in stage 1')

    u(1)=r(1)/bet
    do j = 2, n
      ! Decomposition and forward substitution.
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j-1)*gam(j)
      call assert (bet/=0.0_dp, 'bm_tridiag: Error in stage 2')
      u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j = n-1, 1, -1
      ! Back substitution.
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do

  end subroutine bm_tridiag

  !======================================================================

end module beam_interpolate
