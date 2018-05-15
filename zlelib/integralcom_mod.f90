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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is the _MOD.F90 file replacing the KEEP Blocks, transformed in module
! file using
! keep2mod.pl - 0.0 - By M. Maris - 23 Nov 2000 - 30 May 2003 -
!
! Original Name : re.f90
! Final    Name : re/integralcom_obj.f90
! Date          : Fri Nov 24 18:07:07 2000
! Operated by   : M.Maris
!
! :Update: 30 MAY 2003 :
!   REAL converted to DOUBLE PRECISION.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE integralcom_mod

!KEEP,INTEGRALCOM.
!
! Contains definitions and commons to pass values to integral arguments
!

!
! IntegralEps is the integration precision (real)
! IntegralN is the Number of loops (integer)
! IntegralMAxN is the Maximum number of loops (integer)
!
      Double Precision IntegralEps
      Integer IntegralN, IntegralMaxN

END MODULE
