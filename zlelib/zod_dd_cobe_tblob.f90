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
MODULE ZOD_DD_COBE_TBLOB
!
! ZOD_DD_COBE_TBLOB 0.0 By M. Maris - 6 July 2004 -
!
! Class to handle the Density Distribution, COBE, TBLOB
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

Character(len=*), Parameter :: ZOD_DD_COBE_TBLOB_Version = 'ZOD_DD_COBE_TBLOB 0.0 By M. Maris - 6 July 2004 -'

Type T_DensDstr_COBE_TBLOB
!
! COBE Density Distribution
! TBLOB Component
!
    ! Description Parameters
    Integer        :: ICOMPONENT   ! Code of the distribution
    Character(64)  :: Name         ! Name of the distribution
    Character(64)  :: Comment1     ! First  line of comment
    Character(64)  :: Comment2     ! Second line of comment
    Character(64)  :: Comment3     ! Third  line of comment

    ! Parameters of the distribution
  Double Precision :: nTB              ! [UA^-1] density at 3 AU
  Double Precision :: RTB              ! [UA] Radius of peak density
  Double Precision :: sigma_RTB        ! [UA] Radial Dispersion
  Double Precision :: sigma_zTB        ! [UA] Vertical Dispersion
  Double Precision :: theta_TB         ! [deg] longitude w.r.t. Earth
  Double Precision :: sigma_theta      ! [deg] longitude dispersion
  Double Precision :: i                ! [deg] Inclination
  Double Precision :: Omega            ! [deg] Ascending node
  Double Precision :: X0               ! [AU]
  Double Precision :: Y0               ! [AU]
  Double Precision :: Z0               ! [AU]
  Double Precision :: T0               ! [K]
  Double Precision :: delta            !

  ! Processing Flags
  Logical          :: Completed        ! True if Derived Parameters are generated

  ! Derived Parameters
  Double Precision :: theta_rad        ! [rad] theta in rad
  Double Precision :: sigma_theta_rad  ! [rad] longitude dispersion in rad
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
                                             ! the CRING density to be zero
End Type T_DensDstr_COBE_TBLOB

CONTAINS

subroutine ZOD_DensDstr_COBE_TBLOB_GC(This)
!
! Completes the initializzation for the COBE density distribution TBLOB
!
    Type (T_DensDstr_COBE_TBLOB), Intent(OUT) :: This

! Initializes secondary members used to accelerate calculations
    This%i_rad = This%i*pi/180.d0
    This%Omega_rad = This%Omega*pi/180.d0

    This%theta_rad = This%theta_TB*pi/180.d0
    This%sigma_theta_rad = This%sigma_theta*pi/180.d0

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
end subroutine ZOD_DensDstr_COBE_TBLOB_GC


subroutine ZOD_DensDstr_COBE_TBLOB_New(This,CNAME)
!
! Initializes a new density distribution of COBE type for TBLOBS
!
    Type (T_DensDstr_COBE_TBLOB), Intent(OUT) :: This
    Character(len=70), Optional, Intent(OUT) :: CNAME

    Call ZOD_DensDstr_COBE_New(This%DD)

    CNAME = 'TBLOB'

    This%DD%DD%MinR = 0.1   ! AU
    This%DD%DD%MaxR = 5.2   ! AU

    ! The Cobe Model
      THIS%nTB = 1.9d-8                ! [AU^-1] density at 1 AU
      THIS%RTB = 1.06d0                ! [UA] Radius of peak density
      THIS%sigma_RTB = 0.010           ! [UA] Radial Dispersion
      THIS%sigma_zTB = 0.091           ! [UA] Vertical Dispersion
      THIS%theta_TB  = -10.            ! [deg] lognitude w.r.t. Earth
      THIS%sigma_theta  = 12.1         ! [deg] lognitude dispersion
      THIS%i = 0.49                    ! [deg] Inclination
      THIS%Omega  = 22.3               ! [deg] Ascending node
      THIS%X0     = 0.0                ! [AU]
      THIS%Y0     = 0.0                ! [AU]
      THIS%Z0     = 0.0                ! [AU]
      THIS%T0     = 286.               ! [K]
      THIS%delta  = 0.467              !

    This%CutOff_Radius_Squared = 5.2*5.2 ! AU^2

   ! Completes the definitions
   call ZOD_DensDstr_COBE_TBLOB_GC(This)

end subroutine ZOD_DensDstr_COBE_TBLOB_New

Double Precision Function ZOD_density_COBE_TBLOB(HPos,MeanEarthLongRad,DistPar)
!
! for an heliocentr position HPos = (X, Y, Z) (UA) computes the local 3D optical density
! for a given set of model parameters P for the density distribution
!
! If the squared heliocentric distance is greater than
! DistPar%CutOff_Radius_Squared, density = 0.
!
IMPLICIT NONE

    Double Precision, Dimension(3), Intent(IN)  :: HPos
    Double Precision, Intent(IN) :: MeanEarthLongRad
    Type (T_DensDstr_COBE_TBLOB),    Intent(IN)  :: DistPar

! local variables
    Double Precision, Dimension(3) :: HPosP

    Double Precision :: R, RcS,Rc
    Double Precision :: zR, csiB, Long,theta,absdt
    Double Precision :: F1, F2,F3, density

    HPosP = HPos - DistPar%Center
    RcS = sum(HPosP*HPosP)
    if (RcS.gt.DistPar%CutOff_Radius_Squared) then
       ZOD_density_COBE_TBLOB = 0.d0
    else
       R = sqrt(sum(HPos*HPos))
       Rc  = sqrt(RcS)

       ! This is the eq.(5) of the paper
       zR =     HPosP(1)*DistPar%SOSI  &
            & - HPosP(2)*DistPar%COSI  &
            & + HPosP(3)*DistPar%Cos_i

       csiB = abs(zR/DistPar%sigma_zTB)

       F1 = ((Rc - DistPar%RTB)/DistPar%sigma_RTB)**2/2d0
       F1 = exp(-F1)
       F2 = exp(-csiB)

       Long = atan2(HPosP(2),HPosP(1));
       if (Long .lt. 0) Long = 2.d0*pi + Long

       theta = Long - MeanEarthLongRad

       absdt = (theta - DistPar%theta_Rad)
       absdt = absdt - (2.d0*pi)*floor(absdt/(2.d0*pi));
       if (absdt .gt. pi) absdt = 2.d0*pi - absdt

       F3 = exp(-(absdt/DistPar%sigma_theta_Rad)**2/2.d0)

       density = DistPar%nTB*F1*F2*F3

       ZOD_density_COBE_TBLOB = density
    endif
end function ZOD_density_COBE_TBLOB

subroutine ZOD_DensDstr_COBE_TBLOB_GETP(This,IUNIT,VERBOSE)
!
! Gets a new density distribution of COBE TBLOB Type
!
USE ZOD_P_DensDstr_COMMON

    Type (T_DensDstr_COBE_TBLOB), Intent(OUT) :: This
    Integer, Optional, Intent(IN) :: IUNIT
    Integer, Optional, Intent(IN) :: VERBOSE

    Call ZOD_P_DensDstr_COMMON_RESET

    ! Declares uncompleted parameters
    This%Completed = .false.

    Call ZOD_DensDstr_COBE_New(This%DD)

    if (present(IUNIT)) then
       read(IUNIT,NML = DensDstr_COBE_TBLOB)
    else
       read(*,NML = DensDstr_COBE_TBLOB)
    endif

    This%DD%DD%MinR = MinR  ! AU
    This%DD%DD%MaxR = MaxR  ! AU

    This%Name     = Name
    This%Comment1 = Comment1
    This%Comment2 = Comment2
    This%Comment3 = Comment3

    This%ICOMPONENT = ICOMPONENT

    This%i     = i       ! deg
    This%Omega = Omega       ! deg

    This%X0    = X0     ! AU
    This%Y0    = Y0    ! AU
    This%Z0    = Z0   ! AU

    This%T0    = T0       ! K
    This%delta = delta      ! ad

    This%CutOff_Radius_Squared = CutOff_Radius_Squared ! AU^2

    ! component related values
    THIS%nTB = nTB                ! [AU^-1] density at 1 AU
    THIS%RTB = RTB                ! [UA] Radius of peak density
    THIS%sigma_RTB = sigma_RTB    ! [UA] Radial Dispersion
    THIS%sigma_zTB = sigma_zTB    ! [UA] Vertical Dispersion
    THIS%theta_TB  = theta_TB     ! [deg] lognitude w.r.t. Earth
    THIS%sigma_theta  = sigma_theta   ! [deg] lognitude dispersion

    print*,'TBLOB Cloud Geometry'
    print*,'====================='
    print*,'Inclination = ',This%i
    print*,'Omega       = ',This%Omega
    print*
    print*,'X0     [AU] = ',This%X0
    print*,'Y0     [AU] = ',This%Y0
    print*,'Z0     [AU] = ',This%Z0
    print*


    ! Completes the definitions
    call ZOD_DensDstr_COBE_TBLOB_GC(This)

    ! if verbose shows what is readed
    if (present(VERBOSE)) then
       if (VERBOSE.gt.0) then
          call ZOD_DensDstr_COBE_TBLOB_SHP(This)
       endif
    endif

end subroutine ZOD_DensDstr_COBE_TBLOB_GETP

subroutine ZOD_DensDstr_COBE_TBLOB_SHP(This,IOUT)
!
! Shows a new density distribution of COBE type
!
! If IOUT is specified
!
USE ZOD_P_DensDstr_COMMON

    Type (T_DensDstr_COBE_TBLOB), Intent(IN) :: This
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

    Name = This%Name
    Comment1 = This%Comment1
    Comment2 = This%Comment2
    Comment3 = This%Comment3
    ICOMPONENT = This%ICOMPONENT

    ! Component related section
    nTB = THIS%nTB                ! [AU^-1] density at 1 AU
    RTB = THIS%RTB                ! [UA] Radius of peak density
    sigma_RTB = THIS%sigma_RTB    ! [UA] Radial Dispersion
    sigma_zTB = THIS%sigma_zTB     ! [UA] Vertical Dispersion
    theta_TB = THIS%theta_TB       ! [deg] lognitude w.r.t. Earth
    sigma_theta = THIS%sigma_theta   ! [deg] lognitude dispersion

    if (.not.present(IOUT)) then
        print*
        print*,'Density Distribution COBE TBLOB'
        print*,'================================'
        print*
        if (.not.This%Completed) then
            print*,'Beware: '
            print*,'Beware: NOT COMPLETED OBJECT'
            print*,'Beware: '
        endif
        write(*,NML = DensDstr_COBE_TBLOB)
    else
        if (.not.This%Completed) then
            write(IOUT,*) 'Beware: '
            write(IOUT,*) 'Beware: NOT COMPLETED OBJECT'
            write(IOUT,*) 'Beware: '
        endif
        write(IOUT,NML = DensDstr_COBE_TBLOB)
    endif

end subroutine ZOD_DensDstr_COBE_TBLOB_SHP

END MODULE ZOD_DD_COBE_TBLOB
