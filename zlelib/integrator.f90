!-----------------------------------------------------------------------------
!
!  Copyright (C) 1996-2013 Michele Maris
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

MODULE Integrator
!
! Integrator.f90 - 2.2 - By M. Maris - 19 Jun 1996 - 18 Nov 1997 - 23 Nov 2000 - 24 Nov 2000 - 31 Jul 2002 - 1 Aug 2002 - 30 May 2003 - 23 Aug 2004 -
!
!
! :Update: 31 Jul 2002 : M. Maris :
!   The subroutine TRAPZD has a potential bug since IT is not declared
! with the attribute SAVE it may work properly as it may not, depending
! on the compiler. To be safe SAVE has been added.
!
! :Update: 31 Jul 2002 : M. Maris :
!   A new version of TRAPZD have been added: TRAPZD_ERR it estimated also
! the integration error.
! To support this new version a new, private, data structure T_TRAP_ERR
! instantiated as: P_TRAP_ERR accessible through Integrator_Trap_Err(),
! has been created.
!
! :Warning: 31 Jul 2002 : M. Maris
!   Consecutive calls of TRAPZD_ERR will affect the same P_TRAP_ERR
! variable so, interleaved executions are not allowed.
!
! :Warning: 31 Jul 2002 : M. Maris
!     This is a beta release, after validation clean all the unused parts
!
! :Update: 1 Aug 2002 : M. Maris
!     Added QSIMP_ERR to handle QTRAP_ERR.
!
! :Update: 30 May 2003 : M. Maris
!     REAL -> DOUBLE PRECISION
!
! :Update: 23 Aug 2004 : M. Maris
!     Some correction to allow compatibility with IFOR compiler
!
USE Interpolation
implicit none
private
public romberg,simpson

Character(len=*), Parameter :: INTEGRATOR_VERSION = 'INTEGRATOR.F90 2.3 By M. Maris - 8 Mar 2001 - 23 Aug 2004 - '

CONTAINS

subroutine trapezoidal(func,a,b,res,n)
  integer(i4b), intent(in) :: n
  real(dp), intent(in) :: a,b
  real(dp), intent(inout) :: res
  real(dp), external:: func

  integer(i4b) :: it,j
  real(dp) del,sum
  if (n==1) then
    res=0.5_dp*(b-a)*(func(a)+func(b))
  else
    it=2**(n-2)
    del=(b-a)/real(it,kind=dp)
    sum=0._dp
    do j=1,it
      sum=sum+func(a + del*(j-0.5_dp))
    end do
    res=0.5_dp*(res+(b-a)*sum/real(it,kind=dp))
  endif
end subroutine

function simpson (func,a,b,eps,nmax)
! Estimation of an integral with a control over error.
! The procedure stops when the integration error is < eps
! or when the refinement becomes too high.
!
! The error considered in integration is: eps > 0, absolute error
!                                         eps < 0, relative error
  real(dp), external :: func
  real(dp), intent(in) :: a,b,eps
  integer, intent(in) :: nmax
  real(dp) :: simpson

  integer, parameter :: k=5
  real(dp) res_old, ds
  integer(i4b) n
  logical relative

  res_old = 0._dp
  relative = (eps<=0._dp)

  do n=1,nmax
    call trapezoidal(func,a,b,simpson,n)
    if (n>=k) then
      ds = abs(simpson-res_old)
      if (relative) then
        if (ds<(-eps*0.5_dp*abs(simpson+res_old))) return
      else
        if (ds<eps) return
      endif
    endif
    res_old = simpson
  enddo
  print *,'simpson did not converge'
end function

function romberg(func,a,b,eps,nmax)
! Estimation of an integral with a control over error.
! The procedure stops when the integration error is < eps
! or when the refinement becomes too high.
!
! The error considered in integration is: eps > 0, absolute error
!                                         eps < 0, relative error

  real(dp), external :: func
  real(dp), intent(in) :: a,b, eps
  integer, intent(in) :: nmax
  real(dp) :: romberg

  real(dp) dres

  integer(i4b), parameter :: k=5
  real(dp) s(nmax+1), h(nmax+1)
  logical relative
  integer(i4b) j

  relative = (eps<=0._dp)

  h(1)=1.
  do j=1,nmax
    call trapezoidal(func,a,b,s(j),j)
    if (j>=k) then
      call poly_interpol(h(j-k+1),s(j-k+1),k,0._dp,romberg,dres)
      if (relative) then
        if (abs(dres)<(-eps*abs(romberg))) return
      else
        if (abs(dres)<eps) return
      endif
    endif
    s(j+1)=s(j)
    h(j+1)=0.25_dp*h(j)
  end do
  print*,'romberg did not converge'
end function

END MODULE Integrator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! In this code KEEP Blocks have been transformed in module files using
! keep2mod.pl - 0.0 - By M. Maris - 23 Nov 2000 -
!
! Original Name : integrator.f90
! Final    Name : integrator\integrator.f90
! Date          : Fri Nov 24 18:06:57 2000
! Operated by   : M.Maris
!
! List of replacements:
! !KEEP,integralcom replaced by MODULE integralcom_obj at 159
! KEND. Replaced at 178
!
! Number of Warnings: 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This code has been formatted in F90 + Fortran 90 Guidelines style using
! f772f90.pl - 0.0 - By M. Maris - 23 May 2000 -
!
! Original Name : f77\integrator.for
! Final    Name : f772f90\integrator.f90
! Date          : Thu Nov 23 18:01:58 2000
! Operated by   :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!CMZ :          18/11/97  11.33.58  by  M.Maris
!-- Author :    M.Maris   19/06/96
