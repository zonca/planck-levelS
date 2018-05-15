!-----------------------------------------------------------------------------
!
!  Copyright (C) 2000-2013 Michele Maris
!
!  This file is part of the zodiacal light component of the Planck simulation
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

MODULE Interpolation
implicit none

Character(len=*), Parameter :: Interpolation_Version = 'Interpolation - 2.0 - By M.Maris - 23 May 2002 -'
integer, parameter::dp=kind(1d0)
integer, parameter::i4b=kind(1)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! In this code KEEP Blocks have been transformed in module files using
! keep2mod.pl - 0.0 - By M. Maris - 23 Nov 2000 -
!
! Original Name : interpolation.f90
! Final    Name : interpolation\interpolation.f90
! Date          : Fri Nov 24 18:06:55 2000
! Operated by   : M.Maris
!
! List of replacements:
!
! Number of Warnings: 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This code has been formatted in F90 + Fortran 90 Guidelines style using
! f772f90.pl - 0.0 - By M. Maris - 23 May 2000 -
!
! Original Name : f77\interpolation.for
! Final    Name : f772f90\interpolation.f90
! Date          : Thu Nov 23 18:01:58 2000
! Operated by   :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CMZ :          07/12/96  15.16.24  by  M.Maris
!-- Author :    M.Maris   19/06/96

!
! INTERPOLATION unit collects many general purposes interpolation procedures
!

!cccccccccccccccccccccccccccc
! BEGIN UNIT: INTERPOLATION c
!cccccccccccccccccccccccccccc

!2..!+...1....!....2....!....3....!....4....!....5....!....6....!....7.@..!....8

subroutine poly_interpol (xa,ya,n,x,y,dy)
  integer, intent(in) :: n
  real(dp), intent(in) :: xa(n), ya(n), x
  real(dp), intent(out) :: y,dy

  real(dp) dif,den,c(n),d(n)
  integer i, ns, m

  c=ya
  d=ya
  ns=1
  dif=abs(x-xa(1))
  do i=2,n
    if (abs(x-xa(i))<dif) then
      ns=i
      dif=abs(x-xa(i))
    endif
  end do

  y=ya(ns)
  ns=ns-1
  do m=1,n-1
    do i=1,n-m
      den=(c(i+1)-d(i))/(xa(i)-xa(i+m))
      d(i)=(xa(i+m)-x)*den
      c(i)=(xa(i)-x)*den
    end do
    if ((2*ns)<(n-m)) then
      dy=c(ns+1)
    else
      dy=d(ns)
      ns=ns-1
    endif
    y=y+dy
  end do
end subroutine poly_interpol

!cccccccccccccccccccccccccc
! END UNIT: INTERPOLATION c
!cccccccccccccccccccccccccc

END MODULE Interpolation
