!-----------------------------------------------------------------------------
!
!  This file is part of the Planck simulation package
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version.

!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.

!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!-----------------------------------------------------------------------------

!  The Planck simulation package is being developed at the Max-Planck-Institut
!  fuer Astrophysik and financially supported by the Deutsches Zentrum fuer
!  Luft- und Raumfahrt (DLR).

!  Copyright (C) 2002-2013 Max-Planck-Society
!  \author Martin Reinecke

module minihealpix
use planck_config
use ls_misc_utils
implicit none
private

public pix2ang_ring, ang2pix_ring, ang2vec, vec2ang

real(dp), parameter :: ns_max=8192

contains

subroutine pix2ang_ring(nside, ipix, theta, phi)
!=======================================================================
!     renders theta and phi coordinates of the nominal pixel center
!     for the pixel number ipix (RING scheme)
!     given the map resolution parameter nside
!=======================================================================
  integer(i4b), intent(in) :: ipix, nside
  real(dp), intent(out) :: theta, phi

  integer(i4b) ::  nl2, nl4, ncap, iring, iphi, ip, ipix1
  real(dp) :: fodd, hip, fihip
!-----------------------------------------------------------------------

  ipix1 = ipix + 1 ! in {1, npix}
  nl2 = 2*nside
  ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1

  if (ipix1 <= ncap) then ! North Polar cap -------------
    hip   = ipix1*0.5_dp
    fihip = aint (hip, kind=dp)
    iring = int(sqrt(hip-sqrt(fihip))) + 1 ! counted from North pole
    iphi  = ipix1 - 2*iring*(iring-1)

    theta = acos(1.0_dp - iring**2/(3.0_dp*nside**2))
    phi   = (real(iphi,kind=dp) - 0.5_dp) * pi/(2.0_dp*iring)
  else if (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------
    ip    = ipix1 - ncap - 1
    nl4   = 4*nside
    iring = int(ip/nl4) + nside ! counted from North pole
    iphi  = modulo(ip,nl4) + 1

! 1 if iring+nside is odd, 1/2 otherwise
    fodd  = 0.5_dp * (1 + modulo(iring+nside,2))
    theta = acos((nl2-iring)/(1.5_dp*nside))
    phi   = (real(iphi,kind=dp) - fodd) * pi/(real(nl2,kind=dp))
  else ! South Polar cap -----------------------------------
    ip    = 12*nside*nside - ipix1 + 1
    hip   = ip*0.5_dp
    fihip = aint (hip, kind=dp)
    iring = int(sqrt(hip-sqrt(fihip))) + 1 ! counted from South pole
    iphi  = 4*iring + 1 - (ip-2*iring*(iring-1))

    theta = acos(-1.0_dp + iring**2/(3.0_dp*nside*nside))
    phi   = (real(iphi,kind=dp) - 0.5_dp) * pi/(2.0_dp*iring)
  endif
end subroutine

subroutine ang2pix_ring(nside, theta, phi, ipix)
!=======================================================================
!     renders the pixel number ipix (RING scheme) for a pixel which contains
!     a point on a sphere at coordinates theta and phi, given the map
!     resolution parameter nside
!=======================================================================
  integer(i4b), intent(in) :: nside
  integer(i4b), intent(out) :: ipix
  real(dp), intent(in) ::  theta, phi

  integer(i4b) ::  nl4, jp, jm
  real(dp) ::  z, za, tt, tp, tmp, temp1, temp2
  integer(i4b) ::  ir, ip, kshift

!-----------------------------------------------------------------------
  z = cos(theta)
  za = abs(z)
  tt = modulo(phi,twopi)/halfpi  ! in [0,4)

  if (za<=twothird) then ! Equatorial region
    temp1 = nside*(.5_dp+tt)
    temp2 = nside*.75_dp*z
    jp = int(temp1-temp2) ! index of  ascending edge line
    jm = int(temp1+temp2) ! index of descending edge line

    ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
    kshift = 1 - modulo(ir,2) ! kshift=1 if ir even, 0 otherwise

    nl4 = 4*nside
    ip = int((jp+jm-nside+kshift+1)/2) ! in {0,4n-1}
    if (ip >= nl4) ip = ip - nl4

    ipix = 2*nside*(nside-1) + nl4*(ir-1) + ip
  else ! North & South polar caps -----------------------------
    tp = tt - int(tt)      !MODULO(tt,1.0_dp)
    tmp = nside*sqrt(3.0_dp*(1.0_dp-za))

    jp = int(tp          * tmp ) ! increasing edge line index
    jm = int((1.0_dp-tp) * tmp ) ! decreasing edge line index

    ir = jp + jm + 1        ! ring number counted from the closest pole
    ip = int(tt*ir) ! in {0,4*ir-1}
    if (ip >= 4*ir) ip = ip - 4*ir

    if (z>0._dp) then
      ipix = 2*ir*(ir-1)+ip
    else
      ipix = 12*nside*nside-2*ir*(ir+1)+ip
    endif
  endif
end subroutine

subroutine ang2vec(theta, phi, vector)
!=======================================================================
!     renders the vector (x,y,z) corresponding to angles
!     theta (co-latitude measured from North pole, in [0,Pi] radians)
!     and phi (longitude measured eastward, in radians)
!     North pole is (x,y,z)=(0,0,1)
!=======================================================================
  real(kind=dp), intent(in) :: theta, phi
  real(kind=dp), intent(out), dimension(1:) :: vector

  real(kind=dp) :: sintheta
!-----------------------------------------------------------------------

  call assert ((theta>=0.0_dp) .and. (theta<=pi), &
    "ANG2VEC: theta is out of range [0, Pi]")
  sintheta = SIN(theta)

  vector(1) = sintheta * COS(phi)
  vector(2) = sintheta * SIN(phi)
  vector(3) = COS(theta)

  return
end subroutine ang2vec

subroutine vec2ang(vector, theta, phi)
!=======================================================================
!     renders the angles theta, phi corresponding to vector (x,y,z)
!     theta (co-latitude measured from North pole, in [0,Pi] radians)
!     and phi (longitude measured eastward, in [0,2Pi[ radians)
!     North pole is (x,y,z)=(0,0,1)
!=======================================================================
  real(kind=dp), intent(in), dimension(1:) :: vector
  real(kind=dp), intent(out) :: theta, phi

  real(kind=dp) :: dnorm, z
!-----------------------------------------------------------------------

  dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)

  z = vector(3) / dnorm
  theta = ACOS(z)

  phi = 0.0_dp
  if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
       &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]
  if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[

  return
end subroutine vec2ang

end module minihealpix
