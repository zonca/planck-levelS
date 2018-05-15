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
MODULE ZOD_BI_COBE_SMOOTH
!
! ZOD_BI_COBE_SMOOTH 0.0 By M. Maris - 12 Jul 2004 - 23 Aug 2004 -
!
! Class to handle the Brightness Integral for COBE SMOOTH
!
! This is a sublibrary of ZOD_BRINT.F90
!
! :Beware:
! ========
!    All the private attributes and methods has P_ as a prefix
!
!    Despite the declarations:
!        i) Raileight - Jeans approx is not used
!       ii) Frequency is expected to be passed instead of lambda
!      iii) Results are Scaled in MJ/sterad
!
!   pi is provided through Physical_Parameters_CGS
!
! :Update: 23 Aug 2004 : M. Maris
!    File renamed to ZOD_BI_COBE_SMOOT fromZOD_BRINT_COBE_SMOOT

USE Physical_Parameters_CGS
USE ZOD_Debug
USE ZOD_BlackBody
USE ZOD_DustTemperature
USE ZOD_DensDstr
IMPLICIT NONE

Character(len=*), Parameter :: ZOD_BI_COBE_SMOOTH_VERSION = 'ZOD_BI_COBE_SMOOTH 0.1 By M. Maris - 12 Jul 2004 - 23 Aug 2004 -'

CONTAINS

!
! Smooth component
!

Double Precision function P_FX_COBE_SMOOTH_DENSITY(Distance)
!
! for an heliocentric DISTANCE in AU computes the local 3D optical density
!
! Global Variables
!
! Parameters are passed through Private Members
USE ZOD_BI_BASE

    Double Precision, Intent(IN) :: Distance

! local variables

    Double Precision, Dimension(1:3) :: HPos, P
    Double Precision :: d

    P = Distance * P_MY%Pointing
    HPos = P + P_MY%HSpaceCraft

    d=ZOD_density_COBE_SMOOTH(HPos,P_MY_DensDstr_COBE_SMOOTH)

    P_FX_COBE_SMOOTH_DENSITY = d
End Function P_FX_COBE_SMOOTH_DENSITY

Double Precision function P_FX_COBE_SMOOTH_DENS(Distance)
!
! for an heliocentric DISTANCE in AU computes the local 3D optical density
! scaled by the Black Body brightness .
!
! Parameters are passed through Private Members
!
! The BB is in cgs units and it is assumes Lambda=1cm
!
USE ZOD_BI_BASE

    Double Precision, Intent(IN) :: Distance

! local variables

    Double Precision, Dimension(1:3) :: HPos, P
    Double Precision :: R, d, TR, bb_rj

    ! :Step: Computes the heliocentric position for the given pointing
    P = Distance * P_MY%Pointing
    HPos = P + P_MY%HSpaceCraft

    ! :Step: Computes the dust temperature
    R = Dot_Product(HPos,HPos)
    R = sqrt(R)
    TR = TDust(R=R,DTPar=P_MY%DustTemperaturePar)

    ! :Step: Computes the BB brightness (RJ approx) at lambda=1.0cm
    !    bb_rj = bbl_rj_cgs(1.0,TR)
    !print*,P_MY%FrequencyHz
    bb_rj = bbn_cgs(P_MY%FrequencyHz,TR)
    bb_rj = bb_rj*ScaleToMJ
    !print*,bb_rj

    ! :Step: Computes the density
    d=ZOD_density_COBE_SMOOTH(HPos,P_MY_DensDstr_COBE_SMOOTH)

    P_FX_COBE_SMOOTH_DENS = d*bb_rj
End Function P_FX_COBE_SMOOTH_DENS

Double Precision function P_FX_COBE_SMOOTH_DENS_TN(Distance)
!
! for an heliocentric DISTANCE in AU computes the local 3D optical density
! scaled by the Black Body brightness and T^N.
!
! Parameters are passed through Private Members
!
! The BB is in cgs units and it is assumes Lambda=1cm
!
USE ZOD_BI_BASE

    Double Precision, Intent(IN) :: Distance

! local variables

    Double Precision, Dimension(1:3) :: HPos, P
    Double Precision :: R, d, TR, bb_rj, TN

    ! :Step: Computes the heliocentric position for the given pointing
    P = Distance * P_MY%Pointing
    HPos = P + P_MY%HSpaceCraft

    ! :Step: Computes the dust temperature
    R = Dot_Product(HPos,HPos)
    R = sqrt(R)
    TR = TDust(R=R,DTPar=P_MY%DustTemperaturePar)
    TN = TR**P_MY%Power_of_T

    ! :Step: Computes the BB brightness (RJ approx) at lambda=1.0cm
!    bb_rj = bbl_rj_cgs(1.0,TR)
!bb_rj = bbn_rj_cgs(1.0,TR)
bb_rj = bbn_cgs(P_MY%FrequencyHz,TR)
bb_rj = bb_rj*ScaleToMJ

    ! :Step: Computes the density
    d=ZOD_density_COBE_SMOOTH(HPos,P_MY_DensDstr_COBE_SMOOTH)

    P_FX_COBE_SMOOTH_DENS_TN = d*bb_rj*TN
End Function P_FX_COBE_SMOOTH_DENS_TN

Double Precision function P_FX_COBE_T_SMOOTH_DENS(Distance)
!
! for an heliocentric DISTANCE in AU computes the local 3D optical density
! scaled by T
!
! Parameters are passed through Private Members
!
USE ZOD_BI_BASE

    Double Precision, Intent(IN) :: Distance

! local variables

    Double Precision, Dimension(1:3) :: HPos, P
    Double Precision :: R, d, TR

    ! :Step: Computes the heliocentric position for the given pointing
    P = Distance * P_MY%Pointing
    HPos = P + P_MY%HSpaceCraft

    ! :Step: Computes the dust temperature
    R = Dot_Product(HPos,HPos)
    R = sqrt(R)
    TR = TDust(R=R,DTPar=P_MY%DustTemperaturePar)

    ! :Step: Computes the density
    d=ZOD_density_COBE_SMOOTH(HPos,P_MY_DensDstr_COBE_SMOOTH)

    P_FX_COBE_T_SMOOTH_DENS = d*TR
End Function P_FX_COBE_T_SMOOTH_DENS

END MODULE ZOD_BI_COBE_SMOOTH
