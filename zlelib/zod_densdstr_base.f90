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
MODULE ZOD_DensDstr_BASE
!
! ZOD_DensDstr 0.0 By M. Maris - 6 July 2004 -
!
! Class base to handle the Density Distribution for ZODIACAL LIGHT
!
! Subfile of ZOD_DensDstr.f90
!
USE Physical_Parameters_CGS
USE ZOD_Debug
IMPLICIT NONE

Character(len=*), Parameter :: ZOD_DensDstr_BASE_Version = 'ZOD_DensDstr_BASE 0.0 By M. Maris - 6 July 2004 -'

!
! Private Members used as constants
!

Type T_DensDstr
!
! Base Density Distribution
!
    Double Precision :: MinR ! Minimum Heliocentric Distance (AU)
    Double Precision :: MaxR ! Maximum Heliocentric Distance (AU)
End Type T_DensDstr

Type T_DensDstr_COBE
!
! COBE Density Distribution
!
    Type (T_DensDstr) :: DD
End Type T_DensDstr_COBE


CONTAINS

!
! Subroutines
!

!
! Creators of ZOD structures
!

subroutine ZOD_DensDstr_base_New(This)
!
! Initializes a new density distribution
!

    Type (T_DensDstr), Intent(OUT) :: This

    This%MinR = -1. ! Minimum Heliocentric Distance (AU)
    This%MaxR = -1. ! Maximum Heliocentric Distance (AU)
end subroutine ZOD_DensDstr_base_New

subroutine ZOD_DensDstr_COBE_New(This)
!
! Initializes a new density distribution of COBE type
!
    Type (T_DensDstr_COBE), Intent(OUT) :: This
!    Character(len=70), Optional, Intent(OUT) :: CNAME

    Call ZOD_DensDstr_base_New(This%DD)
end subroutine ZOD_DensDstr_COBE_New

END MODULE ZOD_DensDstr_BASE

MODULE ZOD_P_DensDstr_COMMON
implicit none
!
! Private Local Common Block for ZOD_DensStr
!

    Character(64)  :: Name         ! Name of the distribution
    Character(64)  :: Comment1     ! First  line of comment
    Character(64)  :: Comment2     ! Second line of comment
    Character(64)  :: Comment3     ! Third  line of comment
    Integer :: ICOMPONENT

    Double Precision :: MinR, MaxR, n0, alpha, beta, gamma, mu
    Double Precision :: i, Omega, X0, Y0, Z0, CutOff_Radius_Squared
    Double Precision :: T0, delta

    Double Precision :: delta_RB, nB, deg_delta_gzetaB, vB, pB
    Integer :: IBAND

    Double Precision :: nSR              ! [UA^-1] density at 3 AU
    Double Precision :: RSR              ! [UA] Radius of peak density
    Double Precision :: sigma_RSR        ! [UA] Radial Dispersion
    Double Precision :: sigma_zSR        ! [UA] Vertical Dispersion

    Double Precision :: nTB              ! [UA^-1] density at 3 AU
    Double Precision :: RTB              ! [UA] Radius of peak density
    Double Precision :: sigma_RTB        ! [UA] Radial Dispersion
    Double Precision :: sigma_zTB        ! [UA] Vertical Dispersion
    Double Precision :: theta_TB         ! [deg] longitude w.r.t. Earth
    Double Precision :: sigma_theta      ! [deg] longitude dispersion

NAMELIST/DensDstr_COBE_SMOOTH/ &
        ICOMPONENT, Name, Comment1, Comment2, Comment3, &
        MinR, MaxR, &
        i, Omega, X0, Y0, Z0, CutOff_Radius_Squared, &
        T0, delta, &
        n0, alpha, beta, gamma, mu

NAMELIST/DensDstr_COBE_BAND/ &
        ICOMPONENT, Name, Comment1, Comment2, Comment3, &
        MinR, MaxR, &
        i, Omega, X0, Y0, Z0, &
        T0, delta, &
        IBAND, nB, deg_delta_gzetaB, vB, pB, delta_RB

NAMELIST/DensDstr_COBE_CRING/ &
        ICOMPONENT, Name, Comment1, Comment2, Comment3, &
        MinR, MaxR, &
        i, Omega, X0, Y0, Z0, &
        T0, delta, &
        nSR, RSR, sigma_RSR, sigma_zSR

NAMELIST/DensDstr_COBE_TBLOB/ &
        ICOMPONENT, Name, Comment1, Comment2, Comment3, &
        MinR, MaxR, &
        i, Omega, X0, Y0, Z0, &
        T0, delta, &
        nTB, RTB, sigma_RTB, sigma_zTB, theta_TB, sigma_theta

CONTAINS

SUBROUTINE ZOD_P_DensDstr_COMMON_RESET()
!
! resets the content of SUBROUTINE ZOD_P_DensDstr_COMMON
!

Double Precision :: DINFTY = -1d150
Integer :: INFTY = -1000000

    Name = ''
    Comment1 = ''
    Comment2 = ''
    Comment3 = ''

!!!!!!!!!!!!!!!!!!!!
    ICOMPONENT = INFTY

    MinR = DINFTY
    MaxR = DINFTY

    i = DINFTY
    Omega = DINFTY
    X0 = DINFTY
    Y0 = DINFTY
    Z0 = DINFTY
    CutOff_Radius_Squared = DINFTY

    T0 = DINFTY
    delta  = DINFTY

!!!!!!!!!!!!!!!!!!!!
    n0 = DINFTY
    alpha = DINFTY
    beta = DINFTY
    gamma = DINFTY
    mu = DINFTY

    delta_RB  = DINFTY
    nB  = DINFTY
    deg_delta_gzetaB  = DINFTY
    vB  = DINFTY
    pB  = DINFTY
    IBAND = INFTY

    nSR = DINFTY
    RSR = DINFTY
    sigma_RSR = DINFTY
    sigma_zSR = DINFTY

    nTB = DINFTY
    RTB = DINFTY
    sigma_RTB = DINFTY
    sigma_zTB  = DINFTY
    theta_TB  = DINFTY
    sigma_theta  = DINFTY

RETURN
END SUBROUTINE ZOD_P_DensDstr_COMMON_RESET

END MODULE ZOD_P_DensDstr_COMMON
