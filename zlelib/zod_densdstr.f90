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
MODULE ZOD_DensDstr
!
! ZOD_DensDstr 0.3 By M. Maris - 8 Mar 2001 - 20 May 2004 - 24 May 2004 - 21 June 2004 -
!
! Class to handle the Density Distribution for ZODIACAL LIGHT
!
! A density distribution is created with default values        with ZOD_DensDstr_NEW,
! parameters different from default are readed from files      with ZOD_DensDstr_GET,
! density distributions precalculated parameters are completed with ZOD_DensDstr_COMPLETE.
!
! :Beware:
! ========
!    All the private attributes and methods has P_ as a prefix
!
! :Note: M .Maris : 19 Jun 2001
! ======
!   The Flight Simulator generates vectorial positions for Planck,
!   the Sun and other Solar System objects in the Baricentric
!   Reference System
!
! :Note: M. Maris : 30 Jul 2002
!   Added the parameter CutOff_Radius_Squared to fix from where, in the Solar System,
!   the density distribution has to be considered null.
!   The parameter is squared to avoid the need of calculate a square root.
!
!
! :Update: M. Maris, S. Fogliani : 9 Jul 2003 :
!   Added densities for COBE Bands, Circumsolar Ring and Trailing Blob.
!
! :Note: M. Maris : 9 Jul 2003
! To compute ZOD_density_COBE_TBLOB one has to calculate the mean earth longitude
! for the given epoch. Shall arrive from the Flight Simulator or shall be
! computed using the pointing list? See what done with MatLab
!
! :Update: M. Maris : 20 May 2004
! Added SMOOTH0 density distribuzione: SMOOTH with a not tilted, not shifted plane
!
! :Update: M. Maris : 24 May 2004
! In the calculation of the SMOOTH component, the calculation of sines and
! cosines have been moved from calculation at integration time to calculation
! at initialization time.
!
! :Update: M. Maris : 21 June 2004
! Added comments and identification strings to density objects
!
! :Update: M. Maris : 21 June 2004
! Added the possiblity to handle cloud parameters loaded from outside.
!
! :Update: M. Maris : 6 July 2004
!    Sublibraries scorporated and created
!

USE ZOD_Debug
USE ZOD_DensDstr_BASE
USE ZOD_DD_COBE_SMOOTH
USE ZOD_DD_COBE_BAND
USE ZOD_DD_COBE_CRING
USE ZOD_DD_COBE_TBLOB
IMPLICIT NONE

Character(len=*), Parameter :: ZOD_DensDstr_Version = 'ZOD_DensDstr 0.4 By M. Maris - 8 Mar 2001 - 21 June 2004 - 6 July 2004 -'

!
! Generic interface to DensDstr_New
!
INTERFACE ZOD_DensDstr_New
    MODULE PROCEDURE ZOD_DensDstr_base_New
    MODULE PROCEDURE ZOD_DensDstr_COBE_New
    MODULE PROCEDURE ZOD_DensDstr_COBE_SMOOTH_New
    MODULE PROCEDURE ZOD_DensDstr_COBE_BAND_New
    MODULE PROCEDURE ZOD_DensDstr_COBE_CRING_New
    MODULE PROCEDURE ZOD_DensDstr_COBE_TBLOB_New
END INTERFACE

!
! Generic inteface to DensDstr_GET
!
INTERFACE ZOD_DensDstr_GET
    MODULE PROCEDURE ZOD_DensDstr_COBE_SMOOTH_GETP
    MODULE PROCEDURE ZOD_DensDstr_COBE_BAND_GETP
    MODULE PROCEDURE ZOD_DensDstr_COBE_CRING_GETP
    MODULE PROCEDURE ZOD_DensDstr_COBE_TBLOB_GETP
END INTERFACE

!
! Generic inteface to DensDstr_SHOW
!
INTERFACE ZOD_DensDstr_SHOW
    MODULE PROCEDURE ZOD_DensDstr_COBE_SMOOTH_SHP
    MODULE PROCEDURE ZOD_DensDstr_COBE_BAND_SHP
    MODULE PROCEDURE ZOD_DensDstr_COBE_CRING_SHP
    MODULE PROCEDURE ZOD_DensDstr_COBE_TBLOB_SHP
END INTERFACE

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PRIVATE MEMBERS !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!


!
! Private Members used to pass parameters to a Function as a side effect
!

 ! Used by Func X to pass parameters from and to sub functions
  TYPE, PRIVATE :: P_T_MY
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
  END TYPE P_T_MY

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PRIVATE METHODS !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

PRIVATE :: ZOD_DensDstr_base_New
PRIVATE :: ZOD_DensDstr_COBE_New
PRIVATE :: ZOD_DensDstr_COBE_SMOOTH_New
PRIVATE :: ZOD_DensDstr_COBE_BAND_New
PRIVATE :: ZOD_DensDstr_COBE_CRING_New
PRIVATE :: ZOD_DensDstr_COBE_TBLOB_New

PRIVATE :: ZOD_DensDstr_COBE_SMOOTH_GETP
PRIVATE :: ZOD_DensDstr_COBE_BAND_GETP
PRIVATE :: ZOD_DensDstr_COBE_CRING_GETP
PRIVATE :: ZOD_DensDstr_COBE_TBLOB_GETP

PRIVATE :: ZOD_DensDstr_COBE_SMOOTH_SHP
PRIVATE :: ZOD_DensDstr_COBE_BAND_SHP
PRIVATE :: ZOD_DensDstr_COBE_CRING_SHP
PRIVATE :: ZOD_DensDstr_COBE_TBLOB_SHP

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MISCELLANEOUS SERVICES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZOD_DENSDSTR_STARTUP()
end subroutine ZOD_DENSDSTR_STARTUP


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PRIVATE METHODS !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ZOD_DensDstr
