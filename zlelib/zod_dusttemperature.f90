!-----------------------------------------------------------------------------
!
!  Copyright (C) 2001-2013 Michele Maris
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
MODULE ZOD_DustTemperature
!
! ZOD_DustTemperature 1.0 By M. Maris - 8 Mar 2001 - 10 Sep 2002 - 8 Jul 2003 -
!
! Class to handle the Dust Temperature for ZODIACAL LIGHT
!
! :Update: 10 Sep 2002 : M. Maris :
!  The ZOD_DustTemperature_Show hase been updated
!
! :BUG: 26 June 2003 : M. Maris :
!  The Dust Temperature function is wrongth!!!! T increases with Rh instead of
! to decrease!!!!
!
! :BUGFIXED: 26 June 2003 : M. Maris :
!    The bug have been fixed
!
! :Update: 8 Jul 2003 : M. Maris :
!   All REAL data have been converted to Double Precision
!

Use ZOD_Debug
IMPLICIT NONE

Character(len=*), Parameter :: ZOD_DustTemperature_Version &
 & = 'ZOD_DustTemperature 1.0 By M. Maris - 8 Mar 2001 - 10 Sep 2002 - 8 July 2003 - '

Type T_DustTemperature
!
! Base Dust Temperature
!
    Double Precision :: T0           ! Temperature at 1 AU
    Double Precision :: SigmaT0      ! Temperature at 1 AU Uncertainty (1sigma)
    Double Precision :: Delta        ! Temperature Radial Dependence Exponent
    Double Precision :: SigmaDelta   ! Temperature Radial Dependence Exponent Uncertainty (1sigma)
End Type T_DustTemperature

CONTAINS

Double Precision function TDust(R,DTPar)
!
! Computes the IPD Dust Temperature as a function of R
!
!   R in AU
!
    Type(T_DustTemperature), Intent(IN) :: DTPar ! The temeprature distribution parameters
    Double Precision, Intent(IN) :: R                        ! The eliocentric distance in AU

! Local Variables
    Double Precision TD

    if (R.le.0) Call Die('ZOD_DustTemperature::TDUST','Negative or Null R')

!    TD = DTPar%T0 * power(R,DTPar%delta)

! :BUG: 26 June 2003 : M. Maris :
!  The Dust Temperature function is wrongth!!!! T increases with Rh instead of
! to decrease!!!!
!
!!! This is the old - incriminated line
!       TD = DTPar%T0*R**DTPar%delta
!
!   This is the bug fixed

    TD = DTPar%T0/R**DTPar%delta
    TDust = TD
end function TDust

Double Precision function SigmaTDust(R,DTPar)
!
! Computes the IPD Dust Temperature Uncertainty as a function of R
!
!   R in AU
!
    Type(T_DustTemperature), Intent(IN) :: DTPar ! The temeprature distribution parameters
    Double Precision, Intent(IN) :: R                  ! The eliocentric distance in AU

! Local Variables
    Double Precision TD
    Double Precision STD
    Double Precision VCR1,VCR2

    if (R.le.0) Call Die('ZOD_DustTemperature::SigmaTDUST','Negative or Null R')

    TD = TDust(R,DTPar)

    VCR1 = (DTPar%SigmaT0/DTPar%T0)
    VCR1 = VCR1*VCR1

! :BUG: 26 June 2003 : M. Maris :
!  The Dust Temperature function is wrongth!!!! T increases with Rh instead of
! to decrease!!!!
!
!!! This is the old - incriminated line
!    VCR2 = log(R) * DTPar%SigmaDelta
!   This is the bug fixed

    VCR2 = log(R) * (-DTPar%SigmaDelta)

    VCR2 = VCR2*VCR2

    STD = sqrt(VCR1+VCR2)

    SigmaTDust = STD
end function SigmaTDust

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! DEFAULT OBJECTS HANDLERS !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZOD_DustTemperature_Set(This,T0,SigmaT0,Delta,SigmaDelta)
!
! Sets parameters to dust temperature
!
IMPLICIT NONE

  Type (T_DustTemperature), Intent(OUT) :: This
  Double Precision, Optional,           Intent(IN)  :: T0
  Double Precision, Optional,           Intent(IN)  :: SigmaT0
  Double Precision, Optional,           Intent(IN)  :: Delta
  Double Precision, Optional,           Intent(IN)  :: SigmaDelta

  if (present(T0)) This%T0 = T0
  if (present(SigmaT0)) This%SigmaT0 = SigmaT0
  if (present(Delta)) This%Delta = Delta
  if (present(SigmaDelta)) This%SigmaDelta = SigmaDelta

end subroutine ZOD_DustTemperature_Set

subroutine ZOD_DustTemperature_New(This)
!
! Initializes a new DustTemperature type
!
    Type (T_DustTemperature), Intent(OUT) :: This

! Type (T_DustTemperature), Parameter  :: T_DustTemperature_ZERO = /278.,0.0,-0.5,0.0/
!    TDustPar = T_DustTemperature_ZERO

    This%T0 = 278.
    This%SigmaT0 = 0.
    This%Delta = -0.5
    This%SigmaDelta = 0.

end subroutine ZOD_DustTemperature_New

Subroutine ZOD_DustTemperature_Destroy(This)
!
! Destroy a dust temperature object
!
    Type (T_DustTemperature), Intent(OUT) :: This

    This%T0    = -1.

End Subroutine ZOD_DustTemperature_Destroy

Subroutine ZOD_DustTemperature_Show(This,Unit,Title)
!
! Shows a dust temperature object
!
    Type (T_DustTemperature), Intent(IN) :: This
    Integer, Optional, Intent(IN) :: Unit ! The unit number
    Character(len=*) , Optional, Intent(IN) :: Title ! The Title

    If (present(UNIT)) then
       if (present(Title)) write(UNIT,*) '! ',Title
       write(UNIT,*) '&ZOD_DustTemperature'
       write(UNIT,*) 'T0         = ',This%T0,','
       write(UNIT,*) 'SigmaT0    = ',This%SigmaT0,','
       write(UNIT,*) 'Delta      = ',This%Delta,','
       write(UNIT,*) 'SigmaDelta = ',This%SigmaDelta
       write(UNIT,*) '/'
    else
       if (present(Title)) print*,Title
       print*,'T0         = ',This%T0
       print*,'SigmaT0    = ',This%SigmaT0
       print*,'Delta      = ',This%Delta
       print*,'SigmaDelta = ',This%SigmaDelta
    EndIF
End Subroutine ZOD_DustTemperature_Show

END MODULE ZOD_DustTemperature
