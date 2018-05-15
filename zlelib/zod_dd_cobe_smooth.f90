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
MODULE ZOD_DD_COBE_SMOOTH
!
! ZOD_DD_COBE_SMOOTH 0.0 By M. Maris - 6 July 2004 -
!
! Class to handle the Density Distribution, COBE, SMOOTH
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

Character(len=*), Parameter :: ZOD_DD_COBE_SMOOTH_Version = 'ZOD_DD_COBE_SMOOTH 0.0 By M. Maris - 6 July 2004 -'

Type T_DensDstr_COBE_SMOOTH
!
! COBE Density Distribution
! SMOOTH Component
!
    ! Description Parameters
    Integer        :: ICOMPONENT   ! Code of the distribution
    Character(64)  :: Name         ! Name of the distribution
    Character(64)  :: Comment1     ! First  line of comment
    Character(64)  :: Comment2     ! Second line of comment
    Character(64)  :: Comment3     ! Third  line of comment

    ! Parameters of the distribution
    Double Precision :: n0      ! Density 0
    Double Precision :: alpha
    Double Precision :: mu
    Double Precision :: beta
    Double Precision :: gamma
    Double Precision :: X0
    Double Precision :: Y0
    Double Precision :: Z0
    Double Precision :: Omega  ! deg
    Double Precision :: i      ! deg
    Double Precision :: T0
    Double Precision :: Delta

  ! Processing Flags
  Logical          :: Completed        ! True if Derived Parameters are generated

  ! Derived Parameters
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
                                               ! the SMOOTH density to be zero
End Type T_DensDstr_COBE_SMOOTH

CONTAINS

subroutine ZOD_DensDstr_COBE_SMOOTH_GC(This)
!
! Completes the initializzation for the COBE density distribution SMOOTH
!
    Type (T_DensDstr_COBE_SMOOTH), Intent(OUT) :: This

    ! This is just a repeated definition usefull for vectorial astrometry
    This%Center(1) = This%X0
    This%Center(2) = This%Y0
    This%Center(3) = This%Z0

    ! Initializes secondary members used to accelerate calculations
    This%i_rad = This%i*pi/180.d0
    This%Omega_rad = This%Omega*pi/180.d0

    This%Cos_Omega = cos(This%Omega_rad)
    This%Sin_Omega = sin(This%Omega_rad)

    This%Cos_i = cos(This%i_rad)
    This%Sin_i = sin(This%i_rad)

    This%SOSI = This%Sin_Omega*This%Sin_i
    This%COSI = This%Cos_Omega*This%Sin_i

    ! Signals that the resource has been completed
    This%Completed = .true.
end subroutine ZOD_DensDstr_COBE_SMOOTH_GC

subroutine ZOD_DensDstr_COBE_Smooth_New(This,ICOMPONENT,CNAME)
!
! Initializes a new density distribution of COBE type
!
! ICOMPONENT = 001 : the model with shifted and tilted dust plane
! ICOMPONENT = 002 : the model with dust plane not shifted and not tilted
!

    Type (T_DensDstr_COBE_SMOOTH), Intent(OUT) :: This
    Integer, Intent(IN) :: ICOMPONENT
    Character(len=70), Optional, Intent(OUT) :: CNAME

    ! Declares uncompleted parameters
    This%Completed = .false.

    if ((ICOMPONENT.ne.001).and.(ICOMPONENT.ne.002)) then
       print*,'Error: in ZOD_DensDstr_COBE_Smooth_New'
       print*,'ICOMPONENT nor 001 neither 002'
       print*,'Program Stops'
       stop
    endif

    Call ZOD_DensDstr_COBE_New(This%DD)

    This%DD%DD%MinR = 0.1   ! AU
    This%DD%DD%MaxR = 6.  ! AU

    This%n0    = 1.13e-7    ! AU^-1
    This%alpha = 1.34       ! ad
    This%beta  = 4.14       ! ad
    This%gamma = 0.942      ! ad
    This%mu    = 0.189      ! ad
    This%i     = 2.03       ! deg
    This%Omega = 77.7       ! deg
    This%X0    = 0.0119     ! AU
    This%Y0    = 0.00548    ! AU
! BEWARE: in previous versions This%T0 whas set at -0.00125 AU!!!!
    This%Z0    = -0.00215   ! AU
    This%T0    = 286.       ! K
    This%delta = 0.467      ! ad

    This%CutOff_Radius_Squared = 5.2*5.2 ! AU^2

! if required generates the name
    This%Name     = 'COBE SMOOTH'
    if (present(CNAME)) CNAME = This%Name

    This%Comment1 = 'Default Model defined by: ZOD_DensDstr_COBE_SMOOTH_New'
    This%Comment2 = ''
    This%Comment3 = ''
    This%ICOMPONENT = ICOMPONENT

! put to zero the geometrical tild and shift of the cloud
   if (ICOMPONENT.eq.002) then
       This%i     = 0.         ! deg
       This%Omega = 0.         ! deg

       This%X0    = 0.         ! AU
       This%Y0    = 0.         ! AU
       This%Z0    = 0.         ! AU

       This%Name     = 'SMOOTH0'
       This%Comment2 = 'SMOOTH Cloud with null tilt and shift'

       if (present(CNAME)) CNAME = 'COBE SMOOTH0'

       print*
       print*,'SMOOTH Cloud with null tilt and shift'
       print*
    endif

    print*,'Smooth Cloud Geometry'
    print*,'====================='
    print*,'Inclination = ',This%i
    print*,'Omega       = ',This%Omega
    print*
    print*,'X0     [AU] = ',This%X0
    print*,'Y0     [AU] = ',This%Y0
    print*,'Z0     [AU] = ',This%Z0
    print*


    ! Completes the definitions
    call ZOD_DensDstr_COBE_SMOOTH_GC(This)

end subroutine ZOD_DensDstr_COBE_Smooth_New

subroutine ZOD_DensDstr_COBE_Smooth_GETP(This,IUNIT,VERBOSE)
!
! Gets a new density distribution of COBE Smooth Type
!
USE ZOD_P_DensDstr_COMMON

    Type (T_DensDstr_COBE_SMOOTH), Intent(OUT) :: This
    Integer, Optional, Intent(IN) :: IUNIT
    Integer, Optional, Intent(IN) :: VERBOSE

    Call ZOD_P_DensDstr_COMMON_RESET

    ! Declares uncompleted parameters
    This%Completed = .false.

    Call ZOD_DensDstr_COBE_New(This%DD)

    if (present(IUNIT)) then
       read(IUNIT,NML = DensDstr_COBE_SMOOTH)
    else
       read(*,NML = DensDstr_COBE_SMOOTH)
    endif

    This%DD%DD%MinR = MinR  ! AU
    This%DD%DD%MaxR = MaxR  ! AU

    This%n0    = n0    ! AU^-1
    This%alpha = alpha       ! ad
    This%beta  = beta       ! ad
    This%gamma = gamma      ! ad
    This%mu    = mu      ! ad
    This%i     = i       ! deg
    This%Omega = Omega       ! deg
    This%X0    = X0     ! AU
    This%Y0    = Y0    ! AU
    This%Z0    = Z0   ! AU
    This%T0    = T0       ! K
    This%delta = delta      ! ad

    This%CutOff_Radius_Squared = CutOff_Radius_Squared ! AU^2

    This%Name     = Name
    This%Comment1 = Comment1
    This%Comment2 = Comment2
    This%Comment3 = Comment3
    This%ICOMPONENT = ICOMPONENT

    print*,'Smooth Cloud Geometry'
    print*,'====================='
    print*,'Inclination = ',This%i
    print*,'Omega       = ',This%Omega
    print*
    print*,'X0     [AU] = ',This%X0
    print*,'Y0     [AU] = ',This%Y0
    print*,'Z0     [AU] = ',This%Z0
    print*

    ! Completes the definitions
    call ZOD_DensDstr_COBE_SMOOTH_GC(This)

    ! if verbose shows what is readed
    if (present(VERBOSE)) then
       if (VERBOSE.gt.0) then
          call ZOD_DensDstr_COBE_SMOOTH_SHP(This)
       endif
    endif

end subroutine ZOD_DensDstr_COBE_Smooth_GETP

subroutine ZOD_DensDstr_COBE_Smooth_SHP(This,IOUT)
!
! Shows a new density distribution of COBE type
!
! If IOUT is specified
!
USE ZOD_P_DensDstr_COMMON

    Type (T_DensDstr_COBE_SMOOTH), Intent(IN) :: This
    Integer, Intent(IN), Optional             :: IOUT

    Call ZOD_P_DensDstr_COMMON_RESET

    Name = This%Name
    Comment1 = This%Comment1
    Comment2 = This%Comment2
    Comment3 = This%Comment3
    ICOMPONENT = This%ICOMPONENT

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

    n0 = This%n0     ! AU^-1
    alpha = This%alpha       ! ad
    beta = This%beta        ! ad
    gamma = This%gamma      ! ad
    mu = This%mu            ! ad

    if (.not.present(IOUT)) then
        print*
        print*,'Density Distribution COBE Smooth'
        print*,'================================'
        print*
        if (.not.This%Completed) then
            print*,'Beware: '
            print*,'Beware: NOT COMPLETED OBJECT'
            print*,'Beware: '
        endif
        write(*,NML = DensDstr_COBE_SMOOTH)
    else
        if (.not.This%Completed) then
            write(IOUT,*) 'Beware: '
            write(IOUT,*) 'Beware: NOT COMPLETED OBJECT'
            write(IOUT,*) 'Beware: '
        endif
        write(IOUT,NML = DensDstr_COBE_SMOOTH)
    endif

end subroutine ZOD_DensDstr_COBE_Smooth_SHP

Double Precision Function ZOD_density_COBE_smooth(HPos,DistPar)
!
! for an heliocentr position HPos = (X, Y, Z) (UA) computes the local 3D optical density
! for a given set of model parameters P
!
! If the squared heliocentric distance is greater than
! DistPar%CutOff_Radius_Squared, density = 0.
!
IMPLICIT NONE

    Double Precision, Dimension(3)                      , Intent(IN)  :: HPos
    Type (T_DensDstr_COBE_SMOOTH), Intent(IN)  :: DistPar

!    Type Double Precision                               , Intent(OUT) :: Dens

! local variables
    Double Precision :: Xp, Yp, Zp, Rc, Csi, Dens, Zc, gcsi, fcsi, RadRc

! execution
    Rc = sum(HPos*HPos)
    if (Rc.gt.DistPar%CutOff_Radius_Squared) then
      ZOD_density_COBE_SMOOTH = 0.
    else
      Xp = HPos(1) - DistPar%X0
      Yp = HPos(2) - DistPar%Y0
      Zp = HPos(3) - DistPar%Z0

      Rc = sqrt(Xp*Xp+Yp*Yp+Zp*Zp)

!
!      OmegaRad = DistPar%Omega * pi/180.
!      InclRad  = DistPar%i     * pi/180.
!
!      Zc = Xp*sin(OmegaRad)*sin(InclRad) - &
!           Yp*cos(OmegaRad)*sin(InclRad) + &
!           Zp*cos(InclRad)

      Zc = Xp*DistPar%SOSI - &
           Yp*DistPar%COSI + &
           Zp*DistPar%Cos_i

      csi = abs(Zc/Rc)

      if (csi < DistPar%mu) then
         gcsi = csi*csi/(2.*DistPar%mu)
      else
         gcsi = csi - DistPar%mu/2.
      endif

      fcsi = gcsi**DistPar%gamma
      fcsi = exp(-DistPar%beta * fcsi)

      RadRc = Rc**(-DistPar%alpha)

      dens = DistPar%n0 * RadRc * fcsi

      ZOD_density_COBE_SMOOTH = dens
    endif
end function ZOD_density_COBE_smooth

END MODULE ZOD_DD_COBE_SMOOTH
