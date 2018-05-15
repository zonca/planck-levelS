!-----------------------------------------------------------------------------
!
!  Copyright (C) 2004-2013 Michele Maris
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
MODULE ZOD_BI_BASE
!
! ZOD_BI_BASE 0.0 By M. Maris - 12 Jul 2004 - 23 Aug 2004 -
!
! Class base to handle the Brightness Integrals
!
! This is a sublibrary of ZOD_BRINT.F90
!
! :Beware:
! ========
!    All the private attributes and methods has P_ as a prefix
!
! :Update: 23 Aug 2004 : M. Maris
!    File Renamed to ZOD_BI_BASE from ZOD_BRINT_BASE
!

USE Physical_Parameters_CGS
USE ZOD_Debug
USE ZOD_BlackBody
USE ZOD_DustTemperature
USE ZOD_DensDstr
IMPLICIT NONE

Character(len=*), Parameter :: ZOD_BI_BASE_VERSION = 'ZOD_BI_BASE 0.0 By M. Maris - 12 Jul 2004 - 23 Aug 2004 -'


!
! Private Members used to pass parameters to a Function as a side effect
!
  Type (T_DensDstr_COBE_SMOOTH) :: P_MY_DensDstr_COBE_SMOOTH
  Type (T_DensDstr_COBE_BAND  ) :: P_MY_DensDstr_COBE_BAND
  Type (T_DensDstr_COBE_CRING ) :: P_MY_DensDstr_COBE_CRING
  Type (T_DensDstr_COBE_TBLOB ) :: P_MY_DensDstr_COBE_TBLOB

!
! Used by Func X to pass parameters from and to sub functions
!
  TYPE :: P_T_MY
    Double Precision :: PLong     ! Pointing Longitude (deg)
    Double Precision :: PLat      ! Pointing Latitude  (deg)

    Double Precision :: SLong     ! Spacecraft Longitude (deg)
    Double Precision :: SLat      ! Spacecraft Latitude  (deg)
    Double Precision :: SDistance ! Spacecraft Heliocentric Distance

    Double Precision, Dimension(1:3) :: Sun ! Solar Baricentric Position Vector
    Double Precision, Dimension(1:3) :: SpaceCraft ! Spacecraft Baricentric Position Vector
    Double Precision, Dimension(1:3) :: HSpaceCraft ! Spacecraft Heliocentric Position Vector
    Double Precision, Dimension(1:3) :: SC_Center ! Spacecraft IDP CLoud Center Position

    Double Precision, Dimension(1:3) :: Pointing ! Pointing

    Integer :: Power_of_T ! Power of T^N by which the density distribution has to be scaled

    Double Precision :: MeanEarthLong     ! [deg] The Earth Mean Longitude
    Double Precision :: MeanEarthLongRad  ! [rad] The Earth Mean Longitude in Radiants

!
! These quantities are usually filled at the first call
!
    Type (T_DustTemperature) :: DustTemperaturePar ! Parameters for the dust temperature
    Integer :: MaxPower ! Maximum Power of T
    Integer :: IMethod ! Integration Method
    Integer :: MaxNLoops ! Maximum Number of loops
    Double Precision :: Eps ! Integration accuracy
    Double Precision :: MinDistance ! Minimum distance over which to integrate
    Double Precision :: MaxDistance ! Maximum distance over which to integrate
    Double Precision  :: FrequencyGHz ! Frequency in GHz
    Double Precision  :: FrequencyHz  ! Frequency in Hz
  END TYPE P_T_MY

! Used by Func X to pass parameters from and to sub functions
  Type (P_T_MY) :: P_MY

!CONTAINS

END MODULE ZOD_BI_BASE
