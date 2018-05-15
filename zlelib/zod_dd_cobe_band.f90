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
MODULE ZOD_DD_COBE_BAND
!
! ZOD_DD_COBE_BAND 0.0 By M. Maris - 6 July 2004 -
!
! Class to handle the Density Distribution, COBE, BAND
!
! Subfile of ZOD_DensDstr.f90
!
! :Note: M .Maris : 6 July 2004
! ======
!   Library created.
!

USE ZOD_Debug
USE ZOD_DensDstr_BASE
IMPLICIT NONE

Character(len=*), Parameter :: ZOD_DD_COBE_BAND_Version = 'ZOD_DD_COBE_BAND 0.0 By M. Maris - 6 July 2004 -'

Type T_DensDstr_COBE_BAND
!
! COBE Density Distribution
! BAND Component
!

    ! Description Parameters
    Integer        :: ICOMPONENT   ! Code of the distribution
    Character(64)  :: Name         ! Name of the distribution
    Character(64)  :: Comment1     ! First  line of comment
    Character(64)  :: Comment2     ! Second line of comment
    Character(64)  :: Comment3     ! Third  line of comment

    ! Parameters of the distribution
  Integer          :: IBAND            ! The index of the band for which parameters are stored
  Double Precision :: nB               ! [UA^-1] density at 3 AU
  Double Precision :: deg_delta_gzetaB ! [deg] Shape parameter
  Double Precision :: vB               ! Shape parameter
  Double Precision :: pB               ! Shape parameter
  Double Precision :: i                ! [deg] Inclination
  Double Precision :: Omega            ! [deg] Ascending node
  Double Precision :: delta_RB         ! [AU] Inner radial cut off
  Double Precision :: X0               ! [AU]
  Double Precision :: Y0               ! [AU]
  Double Precision :: Z0               ! [AU]
  Double Precision :: T0               ! [K]
  Double Precision :: delta            !

  ! Processing Flags
  Logical          :: Completed        ! True if Derived Parameters are generated

  ! Derived Parameters
  Double Precision :: sin_delta_gzetaB ! the sin of deg_delta_gzetaB
  Double Precision :: i_rad            ! i in radiants
  Double Precision :: omega_rad        ! omega in radiants
  Double Precision :: cos_omega        ! cos(omega)
  Double Precision :: sin_omega        ! sin(omega)
  Double Precision :: cos_i            ! cos(i)
  Double Precision :: sin_i            ! sin(i
  Double Precision :: COSI             ! cos(omega)*sin(omega)
  Double Precision :: SOSI             ! sin(omega)*sin(omega)

  ! Generic parameters for a COBE model
  Type (T_DensDstr_COBE) :: DD
  Double Precision, Dimension(1:3) :: Center ! The cloud center Center = (/X0,Y0,Z0/)
  Double Precision :: CutOff_Radius_Squared  ! The squared eliocentric distance over which to force
                                             ! the BAND density to be zero
End Type T_DensDstr_COBE_BAND

CONTAINS

subroutine ZOD_DensDstr_COBE_BAND_GC(This)
!
! Completes the initializzation for the COBE density distribution BAND
!
    Type (T_DensDstr_COBE_BAND), Intent(INOUT) :: This

! Initializes secondary members used to accelerate calculations
    This%sin_delta_gzetaB = sin(THIS%deg_delta_gzetaB*pi/180.d0)
    This%i_rad = This%i*pi/180.d0
    This%Omega_rad = This%Omega*pi/180.d0

    This%Cos_Omega = cos(This%Omega_rad)
    This%Sin_Omega = sin(This%Omega_rad)

    This%Cos_i = cos(This%i_rad)
    This%Sin_i = sin(This%i_rad)

    This%SOSI = This%Sin_Omega*This%Sin_i
    This%COSI = This%Cos_Omega*This%Sin_i

! This is just a repeated definition usefull for vectorial astrometry
    This%Center(1) = This%X0
    This%Center(2) = This%Y0
    This%Center(3) = This%Z0

    ! Signals that the resource has been completed
    This%Completed = .true.
end subroutine ZOD_DensDstr_COBE_BAND_GC


subroutine ZOD_DensDstr_COBE_BAND_New(This,iband,CNAME)
!
! Initializes a new density distribution of COBE type for BANDS
!
! IBAND = 1, 2, 3 for the three bands
!
    Integer, Intent(IN) :: IBAND
    Type (T_DensDstr_COBE_BAND), Intent(OUT) :: This
    Character(len=70), Optional, Intent(OUT) :: CNAME

    ! Declares uncompleted parameters
    This%Completed = .false.

    Call ZOD_DensDstr_COBE_New(This%DD)

    This%Comment1 = 'Default Model defined by: ZOD_DensDstr_COBE_BAND_New'
    This%Comment2 = ''
    This%Comment3 = ''

    This%DD%DD%MinR = 0.1   ! AU
    This%DD%DD%MaxR = 5.2   ! AU

    THIS%IBAND            = IBAND

    BandSelection : SELECT CASE (IBAND)
        CASE (1)
           ! The First Band
           This%Name             = 'BAND1'
           This%Comment1 = 'Default Model defined by: ZOD_DensDstr_COBE_CBAND_New'
           This%Comment2 = 'BAND 1'
           This%Comment3 = ''
           This%IComponent       = IBAND
           THIS%nB               = 5.59d-10 ! [UA^-1] density at 3 AU
           THIS%deg_delta_gzetaB = 8.78     ! [deg] Shape parameter
           THIS%vB               = 0.10     ! Shape parameter
           THIS%pB               = 4.       ! Shape parameter
           THIS%i                = 0.56     ! [deg] Inclination
           THIS%Omega            = 80.      ! [deg] Ascending node
           THIS%delta_RB         = 1.5      ! [AU] Inner radial cut off
           THIS%X0               = 0.0d0    ! [AU]
           THIS%Y0               = 0.0d0    ! [AU]
           THIS%Z0               = 0.0d0    ! [AU]
           THIS%T0               = 286.     ! [K]
           THIS%delta            = 0.467    !
        CASE (2)
           ! The Second Band
           This%Name             = 'BAND2'
           This%Comment1 = 'Default Model defined by: ZOD_DensDstr_COBE_CBAND_New'
           This%Comment2 = 'BAND 2'
           This%Comment3 = ''
           This%IComponent       = IBAND
           THIS%nB               = 1.99d-9  ! [UA^-1] density at 3 AU
           THIS%deg_delta_gzetaB = 1.99     ! [deg] Shape parameter
           THIS%vB               = 0.90     ! Shape parameter
           THIS%pB               = 4.       ! Shape parameter
           THIS%i                = 1.2      ! [deg] Inclination
           THIS%Omega            = 30.3     ! [deg] Ascending node
           THIS%delta_RB         = 0.94     ! [AU] Inner radial cut off
           THIS%X0               = 0.0d0    ! [AU]
           THIS%Y0               = 0.0d0    ! [AU]
           THIS%Z0               = 0.0d0    ! [AU]
           THIS%T0               = 286.     ! [K]
           THIS%delta            = 0.467    !
        CASE (3)
           ! The Third Band
           This%Name             = 'BAND3'
           This%Comment1 = 'Default Model defined by: ZOD_DensDstr_COBE_CBAND_New'
           This%Comment2 = 'BAND 3'
           This%Comment3 = ''
           This%IComponent       = IBAND
           THIS%nB               = 1.44d-10 ! [UA^-1] density at 3 AU
           THIS%deg_delta_gzetaB = 15.0     ! [deg] Shape parameter
           THIS%vB               = 0.05     ! Shape parameter
           THIS%pB               = 4.       ! Shape parameter
           THIS%i                = 0.8      ! [deg] Inclination
           THIS%Omega            = 80.0     ! [deg] Ascending node
           THIS%delta_RB         = 1.5      ! [AU] Inner radial cut off
           THIS%X0               = 0.0d0    ! [AU]
           THIS%Y0               = 0.0d0    ! [AU]
           THIS%Z0               = 0.0d0    ! [AU]
           THIS%T0               = 286.     ! [K]
           THIS%delta            = 0.467    !
        CASE DEFAULT  ! Unrecongnized Band
             print*
             print*,'Condition Critique ... le IBAND est ne pas reconnuit IBAND = ',IBAND
             print*,'Les ordinateur ferme'
             print*
             stop
    END SELECT BandSelection

! if required generates the name
    if (present(CNAME)) write(CNAME,1000) IBAND
1000 FORMAT('COBE BAND ',I1)

    This%CutOff_Radius_Squared = 5.2*5.2 ! AU^2

   ! Completes the definitions
   call ZOD_DensDstr_COBE_BAND_GC(This)

end subroutine ZOD_DensDstr_COBE_BAND_New

Double Precision Function ZOD_density_COBE_BAND(HPos,DistPar)
!
! for an heliocentr position HPos = (X, Y, Z) (UA) computes the local 3D optical density
! for a given set of model parameters P for the density distribution
!
! If the squared heliocentric distance is greater than
! DistPar%CutOff_Radius_Squared, density = 0.
!
IMPLICIT NONE

    Double Precision, Dimension(3), Intent(IN)  :: HPos
    Type (T_DensDstr_COBE_BAND),    Intent(IN)  :: DistPar

! local variables
    Double Precision, Dimension(3) :: HPosP

    Double Precision :: R, Rc, RcS
    Double Precision :: zB, gzetaB, gzetaB_deltazetaB
    Double Precision :: F1, F2,F3, F4, density

    HPosP = HPos - DistPar%Center
    RcS = sum(HPosP*HPosP)
    if (RcS.gt.DistPar%CutOff_Radius_Squared) then
       ZOD_density_COBE_BAND = 0.d0
    else
       R = sqrt(sum(HPos*HPos))
       Rc  = sqrt(RcS)

       ! This is the eq.(5) of the paper
       zB =     HPosP(1)*DistPar%SOSI  &
            & - HPosP(2)*DistPar%COSI  &
            & + HPosP(3)*DistPar%Cos_i

       gzetaB = abs(zB/Rc)

       gzetaB_deltazetaB = gzetaB/DistPar%sin_delta_gzetaB

       F1 = 3.d0*DistPar%nB/R
       F2 = exp(-gzetaB_deltazetaB**6.d0)
       F3 = DistPar%vB + gzetaB_deltazetaB**DistPar%pB
       F4 = 1.d0 - exp(-(R/DistPar%delta_RB)**20.d0)

       density = F1*F2*F3*F4

       ZOD_density_COBE_BAND = density
    endif
end function ZOD_density_COBE_BAND

subroutine ZOD_DensDstr_COBE_BAND_GETP(This,IUNIT,VERBOSE)
!
! Gets a new density distribution of COBE BAND Type
!
USE ZOD_P_DensDstr_COMMON

    Type (T_DensDstr_COBE_BAND), Intent(OUT) :: This
    Integer, Optional, Intent(IN) :: IUNIT
    Integer, Optional, Intent(IN) :: VERBOSE

    Call ZOD_P_DensDstr_COMMON_RESET

    ! Declares uncompleted parameters
    This%Completed = .false.

    Call ZOD_DensDstr_COBE_New(This%DD)

    if (present(IUNIT)) then
       read(IUNIT,NML = DensDstr_COBE_BAND)
    else
       read(*,NML = DensDstr_COBE_BAND)
    endif

    This%Name     = Name
    This%Comment1 = Comment1
    This%Comment2 = Comment2
    This%Comment3 = Comment3

    This%ICOMPONENT = ICOMPONENT

    This%DD%DD%MinR = MinR  ! AU
    This%DD%DD%MaxR = MaxR  ! AU

    This%i     = i       ! deg
    This%Omega = Omega       ! deg

    This%X0    = X0     ! AU
    This%Y0    = Y0    ! AU
    This%Z0    = Z0   ! AU

    This%T0    = T0       ! K
    This%delta = delta      ! ad

    This%CutOff_Radius_Squared = CutOff_Radius_Squared ! AU^2

    ! component related values
    THIS%nB               = nB ! [UA^-1] density at 3 AU
    THIS%deg_delta_gzetaB = deg_delta_gzetaB ! [deg] Shape parameter
    THIS%vB               = vB ! Shape parameter
    THIS%pB               = pB ! Shape parameter
    THIS%IBAND            = IBAND ! Band Index

    print*,'BAND   Cloud Geometry'
    print*,'====================='
    print*,'Inclination = ',This%i
    print*,'Omega       = ',This%Omega
    print*
    print*,'X0     [AU] = ',This%X0
    print*,'Y0     [AU] = ',This%Y0
    print*,'Z0     [AU] = ',This%Z0
    print*

    ! Completes the definitions
    call ZOD_DensDstr_COBE_BAND_GC(This)

    ! if verbose shows what is readed
    if (present(VERBOSE)) then
       if (VERBOSE.gt.0) then
          call ZOD_DensDstr_COBE_BAND_SHP(This)
       endif
    endif

end subroutine ZOD_DensDstr_COBE_BAND_GETP

subroutine ZOD_DensDstr_COBE_BAND_SHP(This,IOUT)
!
! Shows a new density distribution of COBE type
!
! If IOUT is specified
!
USE ZOD_P_DensDstr_COMMON

    Type (T_DensDstr_COBE_BAND), Intent(IN) :: This
    Integer, Intent(IN), Optional             :: IOUT

    Call ZOD_P_DensDstr_COMMON_RESET

    MinR = This%DD%DD%MinR  ! AU
    MaxR = This%DD%DD%MaxR  ! AU

    i = This%i              ! deg
    Omega = This%Omega      ! deg
    X0    = This%X0         ! AU
    Y0    = This%Y0         ! AU
    Z0    = This%Z0         ! AU
    T0    = This%T0         ! K
    delta = This%delta      ! ad

    CutOff_Radius_Squared = This%CutOff_Radius_Squared   ! AU^2

    NAME = This%Name
    Comment1 = This%Comment1
    Comment2 = This%Comment2
    Comment3 = This%Comment3
    ICOMPONENT = This%ICOMPONENT

    ! Component related section
    nB = THIS%nB ! [UA^-1] density at 3 AU
    deg_delta_gzetaB = THIS%deg_delta_gzetaB ! [deg] Shape parameter
    vB = THIS%vB ! Shape parameter
    pB = THIS%pB ! Shape parameter

    IBAND = This%IBAND

    if (.not.present(IOUT)) then
        print*
        print*,'Density Distribution COBE BAND  '
        print*,'================================'
        print*
        if (.not.This%Completed) then
            print*,'Beware: '
            print*,'Beware: NOT COMPLETED OBJECT'
            print*,'Beware: '
        endif
        write(*,NML = DensDstr_COBE_BAND)
    else
        if (.not.This%Completed) then
            write(IOUT,*) 'Beware: '
            write(IOUT,*) 'Beware: NOT COMPLETED OBJECT'
            write(IOUT,*) 'Beware: '
        endif
        write(IOUT,NML = DensDstr_COBE_BAND)
    endif

end subroutine ZOD_DensDstr_COBE_BAND_SHP

END MODULE ZOD_DD_COBE_BAND
