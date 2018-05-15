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
MODULE ZOD_ColorCorrection
!
! ZOD_ColorCorrection 0.0 By M. Maris - 8 Mar 2001 -
!
! Class to handle the Color Correction for ZODIACAL LIGHT
!
IMPLICIT NONE

Character(len=*), Parameter :: ZOD_ColorCorrection_Version = 'ZOD_ColorCorrection 0.0 By M. Maris - 8 Mar 2001 -'

Type T_ColorCorrection
!
! Base Color Correction
!
  Real :: P
End Type T_ColorCorrection

Type T_ColorCorrection_COBE
!
! COBE Color Correction
!
    Type (T_ColorCorrection) :: CC
End Type T_ColorCorrection_COBE

INTERFACE ColorCorrection_New
    MODULE PROCEDURE ColorCorrection_Base_New
    MODULE PROCEDURE ColorCorrection_COBE_New
END INTERFACE

CONTAINS

Subroutine ColorCorrection_Base_New(This)
!
! Inits a New Color Correction object
!

    Type (T_ColorCorrection), Intent(OUT) :: This
    This%P = 0.0
End Subroutine ColorCorrection_Base_New

Subroutine ColorCorrection_COBE_New(This)
!
! Inits a New Color Correction object of COBE Type
!

    Type (T_ColorCorrection_COBE), Intent(OUT) :: This

    Call ColorCorrection_Base_New(This%CC)

End Subroutine ColorCorrection_COBE_New

END MODULE ZOD_ColorCorrection
