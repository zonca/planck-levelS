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
MODULE ZOD_BRINT_FAST
!
! ZOD_BRINT_FAST 1.0 By M. Maris - 8 Mar 2001 - 1 Aug 2002 - 23 Aug 2004 - 07 Feb 2011 -
!
! Class to handle the Various components of the Brightness Integral, fast integration
! Fast integration just generates ZLE Flux
!    Density * Black Body emission and not ancillary quantities such as Temperature
!    moments.
!
! :Beware:
! ========
!    All the private attributes and methods has P_ as a prefix
!
!  :Beware:
!  ========
!    Despite the declarations:
!        i) Raileight - Jeans approx is not used
!       ii) Frequency is expected to be passed instead of lambda
!      iii) Results are Scaled in MJ/sterad
!
! :Update: 1 Aug 2002 : M. Maris
!   Added to ZOD_BRINT ... the possibility to handle integration error.
!
! :Update: 1 Aug 2002 : M. Maris
!   Added to ZOD_BRINT_SHOW_BASE the possibility to have one single-row records.
!
! :Update: 5 July 2003 : M. Maris, S. Fogliani
!   Added the bands components
!
! :Update: 8 July 2003 : M. Maris
!   Added the Circumsolar Ring (CRING) and Trailing Blob (TBlob) components
!
! :Update: 8 July 2004 : M. Maris
!   pi is provided through Physical_Parameters_CGS
!
! :Update: 12 July 2004 : M. Maris
!   The library has been splitted in
!       ZOD_BRINT.F90
!       ZOD_BRINT_BASE.F90
!       ZOD_BI_COBE_SMOOTH.F90
!       ZOD_BI_COBE_BAND.F90
!       ZOD_BI_COBE_CRING.F90
!       ZOD_BI_COBE_TBLOB.F90
!
! :Update: 23 Aug 2004 : M. Maris
!    Updated to take in account of changes in zod_brint.f90 sublibraries.
!    Depending on the compiler the line between !BEGIN_DUMMY and !END_DUMMY
!    has to be commented or uncommented
!

USE Physical_Parameters_CGS
USE ZOD_Debug
USE ZOD_BlackBody
USE ZOD_DustTemperature
USE ZOD_DensDstr

USE ZOD_BI_BASE
USE ZOD_BI_COBE_SMOOTH
USE ZOD_BI_COBE_BAND
USE ZOD_BI_COBE_CRING
USE ZOD_BI_COBE_TBLOB

IMPLICIT NONE

Character(len=*), Parameter :: ZOD_BRINT_FAST_VERSION = 'ZOD_BRINT_FAST 1.0 By M. Maris - 8 Mar 2001 - 2011 Feb 07 -'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! OBJECTS DEFINITIONS !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This is the maximum power of the T expansion allowed.
! For cleaness, despite it is privare, it is declared here
!
Integer, Parameter, Private :: P_UpperMaxTPower = 5

!
! The Brightness integral component class
!
TYPE T_BRINT_FAST_COMPONENT
   Double Precision :: Dens_BB ! Density x BlackBody = Flux
END TYPE T_BRINT_FAST_COMPONENT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! GENERIC INTERFACES !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTERFACE ZOD_BRINT_FAST_NEW
   MODULE PROCEDURE P_ZOD_BRINT_FAST_NEW_BASE
END INTERFACE ! ZOD_BRINT_FAST_NEW

INTERFACE ZOD_BRINT_FAST_DESTROY
   MODULE PROCEDURE P_ZOD_BRINT_FAST_DESTROY_BASE
END INTERFACE ! ZOD_BRINT_FAST_DESTROY

INTERFACE ZOD_BRINT_FAST_SHOW
   MODULE PROCEDURE P_ZOD_BRINT_FAST_SHOW_BASE
END INTERFACE ! ZOD_BRINT_FAST_SHOW

!INTERFACE ZOD_BRINT_FAST_COBE
!   MODULE PROCEDURE ZOD_BI_FAST_CC_SMOOTH
!   MODULE PROCEDURE ZOD_BI_FAST_CC_BAND
!   MODULE PROCEDURE ZOD_BI_FAST_CC_TBLOB
!   MODULE PROCEDURE ZOD_BI_FAST_CC_CRING
!END INTERFACE ! ZOD_BRINT_FAST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! PRIVATE DEFINITIONS !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PRIVATE P_ZOD_BRINT_FAST_NEW_BASE
PRIVATE P_ZOD_BRINT_FAST_DESTROY_BASE
PRIVATE P_ZOD_BRINT_FAST_SHOW_BASE

!
! Private Members used as constants
!

CONTAINS

subroutine ZOD_BI_FAST_CC_SMOOTH( &
          & THIS,                     &
          & Pointing,                 &
          & SpaceCraft,               &
          & Sun,                      &
          & MinDistance,              &
          & MaxDistance,              &
          & ParMod,                   &
          & DustTemperaturePar,       &
          & IMethod,                  &
          & Eps,                      &
          & MaxNLoops,                &
          & FrequencyGHz,             &
          & NoQuitAfterSetUp,         &
          & SetUp,                    &
          & FirstCall,                &
          & THIS_ERR                  &
          &)
!
! Computes the components for the brightness integral,
! along a given line of sight defined by the spacecraft Pointing
! and the spacecraft position SPACECRAFT.
!
! Results are saved in the THIS structure of type T_BRINT_FAST_COMPONENT.
!
! Positions are in baricentric coordinates.
!
! Black Body is approximated by the Raileight-Jeans formula in cgs units
!
! Optional parameters shall be passed at the first call, then they are saved
! and reused at each next call (see the subsequent comments).
! An internal counter checks that all the required parameters have been passed.
!
! This_Err is used to put in output the integration error
! when IMethod = -1 (trapezoidal integration TRAPZD_ERR)
! It is a T_BRINT_FAST_COMPONENT structure.
!
! :Update: 1 Aug 2002 : M. Maris
!   Added the possibility to handle the integration error.
!
USE Integrator
IMPLICIT NONE

    Type (T_BRINT_FAST_COMPONENT), Intent(OUT) :: THIS ! The data structure with the
                                                  ! information about components
    Double Precision, Dimension (1:3), Intent(IN) :: Pointing    ! Pointing vector
    Double Precision, Dimension (1:3), Intent(IN) :: SpaceCraft  ! SpaceCraft Position
    Double Precision, Dimension (1:3), Intent(IN) :: Sun         ! Sun Position

    Double Precision, Optional, Intent(IN)        :: MinDistance ! Lower limit for the integration
    Double Precision, Optional, Intent(IN)        :: MaxDistance ! Upper limit for the integration

    Type (T_BRINT_FAST_COMPONENT), Optional, Intent(OUT) :: THIS_ERR ! The data structure with the
                                                                ! information about the integration error

    Type (T_DensDstr_COBE_SMOOTH), Optional, Intent(IN) :: ParMod ! Density Distribution Parameters
    Type (T_DustTemperature), Optional, Intent(IN) :: DustTemperaturePar ! Dust Temperature Parameters

    Integer, Optional, Intent(IN) :: IMethod ! Integration Method
                                             !   0 for QSimp (default)
                                             !   1 for Roemberg

    Double Precision,    Optional, Intent(IN) :: FrequencyGHz ! Frequency in GHz

    Double Precision,    Optional, Intent(IN) :: Eps     ! Integration accuracy
                                             !    Eps < 0, relative accuracy (def Eps = -1e-3)
                                             !    Eps > 0, absolute accuracy

    Integer, Optional, Intent(IN) :: MaxNLoops ! Maximum Number of Loops
                                               !    MaxNLoop = 8 (default)

    Logical, Optional, Intent(IN) :: Setup ! Set this to .true. if you like to enter
                                           ! SetUp mode and quit just after the optional
                                           ! parameters setup

    Logical, Optional, Intent(IN) :: NoQuitAfterSetUp ! Set this to true to avoid return
                                                   ! after setup and allow calculation

    Logical, Optional, Intent(IN) :: FirstCall ! Set this to .true. if the current call
                                               ! is the first one (default .false.).
                                               !
                                               ! If .true. the internal counter of optional
                                               ! parameters is resetted and all
                                               ! the not passed optional
                                               ! values are set to their default values (if any).
                                               !
                                               ! It shall be used in conjunction with
                                               ! SetUp = .true.

! Internal Counter to count parameters set
! Uses the Internal Subroutine P_DecrmementCounter
!
! The initial value for T_MY%NFILLS
!  :Warning:
!  =========
!   Please update this constat each time you add a new optional parameter
!   to be setup at startup
!
  Integer, Parameter :: P_INITIAL_NFILLS = 7

! This variable is set to P_INITIAL_NFILLS by ZOD_BRINT_FAST_INIT and
! it is used to check that all the parameters to be passed at the
! first call have been really passed.
!
    Integer, SAVE :: NFills = P_INITIAL_NFILLS

!
! Local Variables
!
    Double Precision :: Result
    Logical :: LFirstCall

! The Execution Starts here

    ! If required sets the P_T_MY%NFILLS counter
    if (present(FirstCall)) then
       if (FirstCall.and.(.not.(present(SetUp)))) Call Die(&
             &  'Error in ZOD_BRINT_FAST::ZOD_BI_FAST_CC_SMOOTH' &
             & ,'A FirstCall is performed outside SetUp')
       if (FirstCall.and.(.not.(SetUp))) Call Die(&
             &  'Error in ZOD_BRINT_FAST::ZOD_BI_FAST_CC_SMOOTH' &
             & ,'A FirstCall is performed outside SetUp')
       if (FirstCall) NFILLS = P_INITIAL_NFILLS
    endif

!
! :Step: if required enters the setup mode
!
    if (present(FirstCall)) then
       LFirstCall = FirstCall
    else
       LFirstCall = .false.
    endif

    if (present(SetUp)) then
!
! Pass the other variables to the P_FX methods
! through the private members
!
! These quantities shall be passed at the first call only
!
    if (.not.PRESENT(IMethod)) then
        if (LFirstCall) then
           P_MY%IMETHOD = 0
           Call P_DecrementCounter()
        endif
    else
       if ((IMETHOD.ne.1).and.(IMETHOD.ne.0).and.(IMETHOD.ne.1)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_SMOOTH', &
          & 'Invalid IMETHOD', &
          & IMethod)
       P_MY%IMETHOD = IMETHOD
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(EPS)) then
        if (LFirstCall) then
           P_MY%EPS = -1e-3
           Call P_DecrementCounter()
        endif
    else
       if ((EPS.eq.0)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_SMOOTH', &
          & 'Null EPS')
       P_MY%EPS = EPS
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MaxNLoops)) then
        if (LFirstCall) then
           P_MY%MaxNLoops = 8
           Call P_DecrementCounter()
        EndIf
    else
       if ((MaxNLoops.le.0).or.(MaxNLoops.gt.15)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_SMOOTH', &
          & 'Null MaxNLoops or larger than 15')
        P_MY%MaxNLoops = MaxNLoops
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(DustTemperaturePar)) then
    else
       P_MY%DustTemperaturePar = DustTemperaturePar
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MinDistance)) then
    else
       P_MY%MinDistance = dble(MinDistance)
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MaxDistance)) then
    else
       P_MY%MaxDistance = dble(MaxDistance)
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(FrequencyGHz)) then
    else
       P_MY%FrequencyGHz = dble(FrequencyGHz)
       P_MY%FrequencyHz  = dble(FrequencyGHz*1.d9)
       Call P_DecrementCounter()
    endif


    if (NFILLS .ne. 0) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_SMOOTH', &
          & 'Not all the required parameters have been passed at 1st call', &
          & NFILLS)

    if (Present(ParMod)) P_MY_DensDstr_COBE_SMOOTH = ParMod

    ! Return unless forbidden
    if (.not.Present(NoQuitAfterSetUp)) Return

    if ( .not.(Present(THIS_ERR) .neqv. (P_MY%IMethod .eq. -1)) ) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_SMOOTH', &
          & 'The presence or absence of THIS_ERR parameter is not compatible with the selected IMethod', &
          & P_MY%IMethod)

   endif ! End of setup mode

!
! These quantities shall be passed at any call
!
! Geometric transforms
!
    P_MY%Sun         = Sun
    P_MY%SpaceCraft  = Spacecraft
    P_MY%HSpaceCraft = Spacecraft - Sun
    P_MY%SC_Center   = SpaceCraft - Sun - P_MY_DensDstr_COBE_SMOOTH%Center

!
! Pointing
!
    P_MY%Pointing    = Pointing

!
! The Interface Ends, The serious game begins
!
! Perform the integrations
!

    ! :Step: Integrates over the dust density scaled by BB in RJ approx
    if (P_MY%IMETHOD .eq. 0) then ! Uses QSIMP Integration
       result=simpson(P_FX_COBE_SMOOTH_DENS,P_MY%MinDistance,P_MY%MaxDistance,P_MY%Eps,P_MY%MaxNLoops)
    else                     ! Uses QROMB Integration
       result=romberg(P_FX_COBE_SMOOTH_DENS,P_MY%MinDistance,P_MY%MaxDistance,P_MY%EPS,P_MY%MaxNLoops)
    endif
     This%Dens_BB       = Result * 1.0

   ! :Beware: 23 Aug 2004 : M. Maris
   !     If present THIS_ERR passes a dummy value
   !     This is required by some compilers
!BEGIN_DUMMY
!   This_Err%Dens_T = -1
!END_DUMMY

CONTAINS
Subroutine P_DecrementCounter()
!
! Used to decrement the NFILLS counter
!
  if (NFills.gt.0) NFills = NFills-1
End Subroutine P_DecrementCounter

END SUBROUTINE ZOD_BI_FAST_CC_SMOOTH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZOD_BI_FAST_CC_BAND( &
          & THIS,                     &
          & Pointing,                 &
          & SpaceCraft,               &
          & Sun,                      &
          & MinDistance,              &
          & MaxDistance,              &
          & ParMod,                   &
          & DustTemperaturePar,       &
          & IMethod,                  &
          & Eps,                      &
          & MaxNLoops,                &
          & FrequencyGHz,             &
          & NoQuitAfterSetUp,         &
          & SetUp,                    &
          & FirstCall,                &
          & THIS_ERR                  &
          &)
!
! Computes the components for the brightness integral,
! along a given line of sight defined by the spacecraft Pointing
! and the spacecraft position SPACECRAFT.
!
! Results are saved in the THIS structure of type T_BRINT_FAST_COMPONENT.
!
! Positions are in baricentric coordinates.
!
! Black Body is approximated by the Raileight-Jeans formula in cgs units
!
! Optional parameters shall be passed at the first call, then they are saved
! and reused at each next call (see the subsequent comments).
! An internal counter checks that all the required parameters have been passed.
!
! This_Err is used to put in output the integration error
! when IMethod = -1 (trapezoidal integration TRAPZD_ERR)
! It is a T_BRINT_FAST_COMPONENT structure.
!
! :Update: 1 Aug 2002 : M. Maris
!   Added the possibility to handle the integration error.
!
! :Update: 3 June 2003 : M. Maris
!   This version of BRINT handles COBE BAND MODEL
!
USE Integrator
IMPLICIT NONE

    Type (T_BRINT_FAST_COMPONENT), Intent(OUT) :: THIS ! The data structure with the
                                                  ! information about components
    Double Precision, Dimension (1:3), Intent(IN) :: Pointing    ! Pointing vector
    Double Precision, Dimension (1:3), Intent(IN) :: SpaceCraft  ! SpaceCraft Position
    Double Precision, Dimension (1:3), Intent(IN) :: Sun         ! Sun Position

    Double Precision, Optional, Intent(IN)        :: MinDistance ! Lower limit for the integration
    Double Precision, Optional, Intent(IN)        :: MaxDistance ! Upper limit for the integration

    Type (T_BRINT_FAST_COMPONENT), Optional, Intent(OUT) :: THIS_ERR ! The data structure with the
                                                                ! information about the integration error

    Type (T_DensDstr_COBE_BAND), Optional, Intent(IN) :: ParMod ! Density Distribution Parameters
    Type (T_DustTemperature), Optional, Intent(IN) :: DustTemperaturePar ! Dust Temperature Parameters

    Integer, Optional, Intent(IN) :: IMethod ! Integration Method
                                             !   0 for QSimp (default)
                                             !   1 for Roemberg

    Double Precision,    Optional, Intent(IN) :: FrequencyGHz ! Frequency in GHz

    Double Precision,    Optional, Intent(IN) :: Eps     ! Integration accuracy
                                             !    Eps < 0, relative accuracy (def Eps = -1e-3)
                                             !    Eps > 0, absolute accuracy

    Integer, Optional, Intent(IN) :: MaxNLoops ! Maximum Number of Loops
                                               !    MaxNLoop = 8 (default)

    Logical, Optional, Intent(IN) :: Setup ! Set this to .true. if you like to enter
                                           ! SetUp mode and quit just after the optional
                                           ! parameters setup

    Logical, Optional, Intent(IN) :: NoQuitAfterSetUp ! Set this to true to avoid return
                                                   ! after setup and allow calculation

    Logical, Optional, Intent(IN) :: FirstCall ! Set this to .true. if the current call
                                               ! is the first one (default .false.).
                                               !
                                               ! If .true. the internal counter of optional
                                               ! parameters is resetted and all
                                               ! the not passed optional
                                               ! values are set to their default values (if any).
                                               !
                                               ! It shall be used in conjunction with
                                               ! SetUp = .true.

! Internal Counter to count parameters set
! Uses the Internal Subroutine P_DecrmementCounter
!
! The initial value for T_MY%NFILLS
!  :Warning:
!  =========
!   Please update this constat each time you add a new optional parameter
!   to be setup at startup
!
  Integer, Parameter :: P_INITIAL_NFILLS = 7

! This variable is set to P_INITIAL_NFILLS by ZOD_BRINT_FAST_INIT and
! it is used to check that all the parameters to be passed at the
! first call have been really passed.
!
    Integer, SAVE :: NFills = P_INITIAL_NFILLS

!
! Local Variables
!
    Double Precision :: Result, Result_err
    Logical :: LFirstCall

! The Execution Starts here

    ! If required sets the P_T_MY%NFILLS counter
    if (present(FirstCall)) then
       if (FirstCall.and.(.not.(present(SetUp)))) Call Die(&
             &  'Error in ZOD_BRINT_FAST::ZOD_BI_FAST_CC_BAND' &
             & ,'A FirstCall is performed outside SetUp')
       if (FirstCall.and.(.not.(SetUp))) Call Die(&
             &  'Error in ZOD_BRINT_FAST::ZOD_BI_FAST_CC_BAND' &
             & ,'A FirstCall is performed outside SetUp')
       if (FirstCall) NFILLS = P_INITIAL_NFILLS
    endif

!
! :Step: if required enters the setup mode
!
    if (present(FirstCall)) then
       LFirstCall = FirstCall
    else
       LFirstCall = .false.
    endif

    if (present(SetUp)) then
!
! Pass the other variables to the P_FX methods
! through the private members
!
! These quantities shall be passed at the first call only
!
    if (.not.PRESENT(IMethod)) then
        if (LFirstCall) then
           P_MY%IMETHOD = 0
           Call P_DecrementCounter()
        endif
    else
       if ((IMETHOD.ne.1).and.(IMETHOD.ne.0).and.(IMETHOD.ne.1)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_BAND', &
          & 'Invalid IMETHOD', &
          & IMethod)
       P_MY%IMETHOD = IMETHOD
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(EPS)) then
        if (LFirstCall) then
           P_MY%EPS = -1e-3
           Call P_DecrementCounter()
        endif
    else
       if ((EPS.eq.0)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_BAND', &
          & 'Null EPS')
       P_MY%EPS = EPS
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MaxNLoops)) then
        if (LFirstCall) then
           P_MY%MaxNLoops = 8
           Call P_DecrementCounter()
        EndIf
    else
       if ((MaxNLoops.le.0).or.(MaxNLoops.gt.15)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_BAND', &
          & 'Null MaxNLoops or larger than 15')
        P_MY%MaxNLoops = MaxNLoops
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(DustTemperaturePar)) then
    else
       P_MY%DustTemperaturePar = DustTemperaturePar
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MinDistance)) then
    else
       P_MY%MinDistance = dble(MinDistance)
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MaxDistance)) then
    else
       P_MY%MaxDistance = dble(MaxDistance)
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(FrequencyGHz)) then
    else
       P_MY%FrequencyGHz = dble(FrequencyGHz)
       P_MY%FrequencyHz  = dble(FrequencyGHz*1.d9)
       Call P_DecrementCounter()
    endif


    if (NFILLS .ne. 0) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_BAND', &
          & 'Not all the required parameters have been passed at 1st call', &
          & NFILLS)

    if (Present(ParMod)) P_MY_DensDstr_COBE_BAND = ParMod

    ! Return unless forbidden
    if (.not.Present(NoQuitAfterSetUp)) Return

    if ( .not.(Present(THIS_ERR) .neqv. (P_MY%IMethod .eq. -1)) ) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_BAND', &
          & 'The presence or absence of THIS_ERR parameter is not compatible with the selected IMethod', &
          & P_MY%IMethod)

   endif ! End of setup mode

!
! These quantities shall be passed at any call
!
! Geometric transforms
!
    P_MY%Sun         = Sun
    P_MY%SpaceCraft  = Spacecraft
    P_MY%HSpaceCraft = Spacecraft - Sun
    P_MY%SC_Center   = SpaceCraft - Sun - P_MY_DensDstr_COBE_BAND%Center

!
! Pointing
!
    P_MY%Pointing    = Pointing

!
! The Interface Ends, The serious game begins
!
! Perform the integrations
!


    ! :Step: Integrates over the dust density scaled by BB in RJ approx
    Result = -1.
    if (P_MY%IMETHOD .eq. 0) then ! Uses QSIMP Integration
       result=simpson(P_FX_COBE_BAND_DENS,P_MY%MinDistance,P_MY%MaxDistance,P_MY%Eps,P_MY%MaxNLoops)
    else                     ! Uses QROMB Integration
       result=romberg(P_FX_COBE_BAND_DENS,P_MY%MinDistance,P_MY%MaxDistance,P_MY%EPS,P_MY%MaxNLoops)
    endif
    This%Dens_BB       = Result * 1.0

!print*,'END OF COBE_T_BAND_DENS'

   ! :Beware: 23 Aug 2004 : M. Maris
   !     If present THIS_ERR passes a dummy value
   !     This is required by some compilers
   if (present(This_Err)) This_Err%Dens_BB = result_err

!print*,'QUITTING '

CONTAINS
Subroutine P_DecrementCounter()
!
! Used to decrement the NFILLS counter
!
  if (NFills.gt.0) NFills = NFills-1
End Subroutine P_DecrementCounter

END SUBROUTINE ZOD_BI_FAST_CC_BAND


!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZOD_BI_FAST_CC_CRING( &
          & THIS,                     &
          & Pointing,                 &
          & SpaceCraft,               &
          & Sun,                      &
          & MinDistance,              &
          & MaxDistance,              &
          & ParMod,                   &
          & DustTemperaturePar,       &
          & IMethod,                  &
          & Eps,                      &
          & MaxNLoops,                &
          & FrequencyGHz,             &
          & NoQuitAfterSetUp,         &
          & SetUp,                    &
          & FirstCall,                &
          & THIS_ERR                  &
          &)
!
! Computes the components for the brightness integral,
! along a given line of sight defined by the spacecraft Pointing
! and the spacecraft position SPACECRAFT.
!
! Results are saved in the THIS structure of type T_BRINT_FAST_COMPONENT.
!
! Positions are in baricentric coordinates.
!
! Black Body is approximated by the Raileight-Jeans formula in cgs units
!
! Optional parameters shall be passed at the first call, then they are saved
! and reused at each next call (see the subsequent comments).
! An internal counter checks that all the required parameters have been passed.
!
! This_Err is used to put in output the integration error
! when IMethod = -1 (trapezoidal integration TRAPZD_ERR)
! It is a T_BRINT_FAST_COMPONENT structure.
!
! :Update: 1 Aug 2002 : M. Maris
!   Added the possibility to handle the integration error.
!
! :Update: 8 July 2003 : M. Maris
!   This version of BRINT handles COBE CRING MODEL
!
USE Integrator
IMPLICIT NONE

    Type (T_BRINT_FAST_COMPONENT), Intent(OUT) :: THIS ! The data structure with the
                                                  ! information about components
    Double Precision, Dimension (1:3), Intent(IN) :: Pointing    ! Pointing vector
    Double Precision, Dimension (1:3), Intent(IN) :: SpaceCraft  ! SpaceCraft Position
    Double Precision, Dimension (1:3), Intent(IN) :: Sun         ! Sun Position

    Double Precision, Optional, Intent(IN)        :: MinDistance ! Lower limit for the integration
    Double Precision, Optional, Intent(IN)        :: MaxDistance ! Upper limit for the integration

    Type (T_BRINT_FAST_COMPONENT), Optional, Intent(OUT) :: THIS_ERR ! The data structure with the
                                                                ! information about the integration error

    Type (T_DensDstr_COBE_CRING), Optional, Intent(IN) :: ParMod ! Density Distribution Parameters
    Type (T_DustTemperature), Optional, Intent(IN) :: DustTemperaturePar ! Dust Temperature Parameters

    Integer, Optional, Intent(IN) :: IMethod ! Integration Method
                                             !   0 for QSimp (default)
                                             !   1 for Roemberg

    Double Precision,    Optional, Intent(IN) :: FrequencyGHz ! Frequency in GHz

    Double Precision,    Optional, Intent(IN) :: Eps     ! Integration accuracy
                                             !    Eps < 0, relative accuracy (def Eps = -1e-3)
                                             !    Eps > 0, absolute accuracy

    Integer, Optional, Intent(IN) :: MaxNLoops ! Maximum Number of Loops
                                               !    MaxNLoop = 8 (default)

    Logical, Optional, Intent(IN) :: Setup ! Set this to .true. if you like to enter
                                           ! SetUp mode and quit just after the optional
                                           ! parameters setup

    Logical, Optional, Intent(IN) :: NoQuitAfterSetUp ! Set this to true to avoid return
                                                   ! after setup and allow calculation

    Logical, Optional, Intent(IN) :: FirstCall ! Set this to .true. if the current call
                                               ! is the first one (default .false.).
                                               !
                                               ! If .true. the internal counter of optional
                                               ! parameters is resetted and all
                                               ! the not passed optional
                                               ! values are set to their default values (if any).
                                               !
                                               ! It shall be used in conjunction with
                                               ! SetUp = .true.

! Internal Counter to count parameters set
! Uses the Internal Subroutine P_DecrmementCounter
!
! The initial value for T_MY%NFILLS
!  :Warning:
!  =========
!   Please update this constat each time you add a new optional parameter
!   to be setup at startup
!
  Integer, Parameter :: P_INITIAL_NFILLS = 7

! This variable is set to P_INITIAL_NFILLS by ZOD_BRINT_FAST_INIT and
! it is used to check that all the parameters to be passed at the
! first call have been really passed.
!
    Integer, SAVE :: NFills = P_INITIAL_NFILLS

!
! Local Variables
!
    Double Precision :: Result
    Logical :: LFirstCall

! The Execution Starts here

    ! If required sets the P_T_MY%NFILLS counter
    if (present(FirstCall)) then
       if (FirstCall.and.(.not.(present(SetUp)))) Call Die(&
             &  'Error in ZOD_BRINT_FAST::ZOD_BI_FAST_CC_CRING' &
             & ,'A FirstCall is performed outside SetUp')
       if (FirstCall.and.(.not.(SetUp))) Call Die(&
             &  'Error in ZOD_BRINT_FAST::ZOD_BI_FAST_CC_CRING' &
             & ,'A FirstCall is performed outside SetUp')
       if (FirstCall) NFILLS = P_INITIAL_NFILLS
    endif

!
! :Step: if required enters the setup mode
!
    if (present(FirstCall)) then
       LFirstCall = FirstCall
    else
       LFirstCall = .false.
    endif

    if (present(SetUp)) then
!
! Pass the other variables to the P_FX methods
! through the private members
!
! These quantities shall be passed at the first call only
!
    if (.not.PRESENT(IMethod)) then
        if (LFirstCall) then
           P_MY%IMETHOD = 0
           Call P_DecrementCounter()
        endif
    else
       if ((IMETHOD.ne.1).and.(IMETHOD.ne.0).and.(IMETHOD.ne.1)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_CRING', &
          & 'Invalid IMETHOD', &
          & IMethod)
       P_MY%IMETHOD = IMETHOD
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(EPS)) then
        if (LFirstCall) then
           P_MY%EPS = -1e-3
           Call P_DecrementCounter()
        endif
    else
       if ((EPS.eq.0)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_CRING', &
          & 'Null EPS')
       P_MY%EPS = EPS
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MaxNLoops)) then
        if (LFirstCall) then
           P_MY%MaxNLoops = 8
           Call P_DecrementCounter()
        EndIf
    else
       if ((MaxNLoops.le.0).or.(MaxNLoops.gt.15)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_CRING', &
          & 'Null MaxNLoops or larger than 15')
        P_MY%MaxNLoops = MaxNLoops
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(DustTemperaturePar)) then
    else
       P_MY%DustTemperaturePar = DustTemperaturePar
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MinDistance)) then
    else
       P_MY%MinDistance = dble(MinDistance)
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MaxDistance)) then
    else
       P_MY%MaxDistance = dble(MaxDistance)
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(FrequencyGHz)) then
    else
       P_MY%FrequencyGHz = dble(FrequencyGHz)
       P_MY%FrequencyHz  = dble(FrequencyGHz*1.d9)
       Call P_DecrementCounter()
    endif


    if (NFILLS .ne. 0) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_CRING', &
          & 'Not all the required parameters have been passed at 1st call', &
          & NFILLS)

    if (Present(ParMod)) P_MY_DensDstr_COBE_CRING = ParMod

    ! Return unless forbidden
    if (.not.Present(NoQuitAfterSetUp)) Return

    if ( .not.(Present(THIS_ERR) .neqv. (P_MY%IMethod .eq. -1)) ) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_CRING', &
          & 'The presence or absence of THIS_ERR parameter is not compatible with the selected IMethod', &
          & P_MY%IMethod)

   endif ! End of setup mode

!
! These quantities shall be passed at any call
!
! Geometric transforms
!
    P_MY%Sun         = Sun
    P_MY%SpaceCraft  = Spacecraft
    P_MY%HSpaceCraft = Spacecraft - Sun
    P_MY%SC_Center   = SpaceCraft - Sun - P_MY_DensDstr_COBE_CRING%Center

!
! Pointing
!
    P_MY%Pointing    = Pointing

!
! The Interface Ends, The serious game begins
!
! Perform the integrations
!

    ! :Step: Integrates over the dust density scaled by BB in RJ approx
    if (P_MY%IMETHOD .eq. 0) then ! Uses QSIMP Integration
       result=simpson(P_FX_COBE_CRING_DENS,P_MY%MinDistance,P_MY%MaxDistance,P_MY%Eps,P_MY%MaxNLoops)
    else                     ! Uses QROMB Integration
       result=romberg(P_FX_COBE_CRING_DENS,P_MY%MinDistance,P_MY%MaxDistance,P_MY%EPS,P_MY%MaxNLoops)
    endif
     This%Dens_BB       = Result * 1.0

!print*,'END OF FLUX'

!print*,'END OF COBE_T_CRING_DENS'

   ! :Beware: 23 Aug 2004 : M. Maris
   !     If present THIS_ERR passes a dummy value
   !     This is required by some compilers
   if (present(This_Err)) This_Err%Dens_BB = -1

CONTAINS
Subroutine P_DecrementCounter()
!
! Used to decrement the NFILLS counter
!
  if (NFills.gt.0) NFills = NFills-1
End Subroutine P_DecrementCounter

END SUBROUTINE ZOD_BI_FAST_CC_CRING


!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZOD_BI_FAST_CC_TBLOB( &
          & THIS,                     &
          & Pointing,                 &
          & SpaceCraft,               &
          & Sun,                      &
          & MinDistance,              &
          & MaxDistance,              &
          & ParMod,                   &
          & DustTemperaturePar,       &
          & IMethod,                  &
          & Eps,                      &
          & MaxNLoops,                &
          & FrequencyGHz,             &
          & NoQuitAfterSetUp,         &
          & SetUp,                    &
          & FirstCall,                &
          & THIS_ERR                  &
          &)
!
! Computes the components for the brightness integral,
! along a given line of sight defined by the spacecraft Pointing
! and the spacecraft position SPACECRAFT.
!
! Results are saved in the THIS structure of type T_BRINT_FAST_COMPONENT.
!
! Positions are in baricentric coordinates.
!
! Black Body is approximated by the Raileight-Jeans formula in cgs units
!
! Optional parameters shall be passed at the first call, then they are saved
! and reused at each next call (see the subsequent comments).
! An internal counter checks that all the required parameters have been passed.
!
! This_Err is used to put in output the integration error
! when IMethod = -1 (trapezoidal integration TRAPZD_ERR)
! It is a T_BRINT_FAST_COMPONENT structure.
!
! :Update: 1 Aug 2002 : M. Maris
!   Added the possibility to handle the integration error.
!
! :Update: 3 June 2003 : M. Maris
!   This version of BRINT handles COBE TBLOB MODEL
!
USE Integrator
IMPLICIT NONE

    Type (T_BRINT_FAST_COMPONENT), Intent(OUT) :: THIS ! The data structure with the
                                                  ! information about components
    Double Precision, Dimension (1:3), Intent(IN) :: Pointing    ! Pointing vector
    Double Precision, Dimension (1:3), Intent(IN) :: SpaceCraft  ! SpaceCraft Position
    Double Precision, Dimension (1:3), Intent(IN) :: Sun         ! Sun Position

    Double Precision, Optional, Intent(IN)        :: MinDistance ! Lower limit for the integration
    Double Precision, Optional, Intent(IN)        :: MaxDistance ! Upper limit for the integration

    Type (T_BRINT_FAST_COMPONENT), Optional, Intent(OUT) :: THIS_ERR ! The data structure with the
                                                                ! information about the integration error

    Type (T_DensDstr_COBE_TBLOB), Optional, Intent(IN) :: ParMod ! Density Distribution Parameters
    Type (T_DustTemperature), Optional, Intent(IN) :: DustTemperaturePar ! Dust Temperature Parameters

    Integer, Optional, Intent(IN) :: IMethod ! Integration Method
                                             !   0 for QSimp (default)
                                             !   1 for Roemberg

    Double Precision,    Optional, Intent(IN) :: FrequencyGHz ! Frequency in GHz

    Double Precision,    Optional, Intent(IN) :: Eps     ! Integration accuracy
                                             !    Eps < 0, relative accuracy (def Eps = -1e-3)
                                             !    Eps > 0, absolute accuracy

    Integer, Optional, Intent(IN) :: MaxNLoops ! Maximum Number of Loops
                                               !    MaxNLoop = 8 (default)

    Logical, Optional, Intent(IN) :: Setup ! Set this to .true. if you like to enter
                                           ! SetUp mode and quit just after the optional
                                           ! parameters setup

    Logical, Optional, Intent(IN) :: NoQuitAfterSetUp ! Set this to true to avoid return
                                                   ! after setup and allow calculation

    Logical, Optional, Intent(IN) :: FirstCall ! Set this to .true. if the current call
                                               ! is the first one (default .false.).
                                               !
                                               ! If .true. the internal counter of optional
                                               ! parameters is resetted and all
                                               ! the not passed optional
                                               ! values are set to their default values (if any).
                                               !
                                               ! It shall be used in conjunction with
                                               ! SetUp = .true.

! Internal Counter to count parameters set
! Uses the Internal Subroutine P_DecrmementCounter
!
! The initial value for T_MY%NFILLS
!  :Warning:
!  =========
!   Please update this constat each time you add a new optional parameter
!   to be setup at startup
!
  Integer, Parameter :: P_INITIAL_NFILLS = 7

! This variable is set to P_INITIAL_NFILLS by ZOD_BRINT_FAST_INIT and
! it is used to check that all the parameters to be passed at the
! first call have been really passed.
!
    Integer, SAVE :: NFills = P_INITIAL_NFILLS

!
! Local Variables
!
    Double Precision :: Result
    Logical :: LFirstCall

! The Execution Starts here

    ! If required sets the P_T_MY%NFILLS counter
    if (present(FirstCall)) then
       if (FirstCall.and.(.not.(present(SetUp)))) Call Die(&
             &  'Error in ZOD_BRINT_FAST::ZOD_BI_FAST_CC_TBLOB' &
             & ,'A FirstCall is performed outside SetUp')
       if (FirstCall.and.(.not.(SetUp))) Call Die(&
             &  'Error in ZOD_BRINT_FAST::ZOD_BI_FAST_CC_TBLOB' &
             & ,'A FirstCall is performed outside SetUp')
       if (FirstCall) NFILLS = P_INITIAL_NFILLS
    endif

!
! :Step: if required enters the setup mode
!
    if (present(FirstCall)) then
       LFirstCall = FirstCall
    else
       LFirstCall = .false.
    endif

    if (present(SetUp)) then
!
! Pass the other variables to the P_FX methods
! through the private members
!
! These quantities shall be passed at the first call only
!
    if (.not.PRESENT(IMethod)) then
        if (LFirstCall) then
           P_MY%IMETHOD = 0
           Call P_DecrementCounter()
        endif
    else
       if ((IMETHOD.ne.1).and.(IMETHOD.ne.0).and.(IMETHOD.ne.1)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_TBLOB', &
          & 'Invalid IMETHOD', &
          & IMethod)
       P_MY%IMETHOD = IMETHOD
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(EPS)) then
        if (LFirstCall) then
           P_MY%EPS = -1e-3
           Call P_DecrementCounter()
        endif
    else
       if ((EPS.eq.0)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_TBLOB', &
          & 'Null EPS')
       P_MY%EPS = EPS
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MaxNLoops)) then
        if (LFirstCall) then
           P_MY%MaxNLoops = 8
           Call P_DecrementCounter()
        EndIf
    else
       if ((MaxNLoops.le.0).or.(MaxNLoops.gt.15)) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_TBLOB', &
          & 'Null MaxNLoops or larger than 15')
        P_MY%MaxNLoops = MaxNLoops
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(DustTemperaturePar)) then
    else
       P_MY%DustTemperaturePar = DustTemperaturePar
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MinDistance)) then
    else
       P_MY%MinDistance = dble(MinDistance)
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(MaxDistance)) then
    else
       P_MY%MaxDistance = dble(MaxDistance)
       Call P_DecrementCounter()
    endif

    if (.not.PRESENT(FrequencyGHz)) then
    else
       P_MY%FrequencyGHz = dble(FrequencyGHz)
       P_MY%FrequencyHz  = dble(FrequencyGHz*1.d9)
       Call P_DecrementCounter()
    endif


    if (NFILLS .ne. 0) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_TBLOB', &
          & 'Not all the required parameters have been passed at 1st call', &
          & NFILLS)

    if (Present(ParMod)) P_MY_DensDstr_COBE_TBLOB = ParMod

    ! Return unless forbidden
    if (.not.Present(NoQuitAfterSetUp)) Return

    if ( .not.(Present(THIS_ERR) .neqv. (P_MY%IMethod .eq. -1)) ) &
          & call die('ZOD_BRINT_FAST::ZOD_BI_FAST_CC_TBLOB', &
          & 'The presence or absence of THIS_ERR parameter is not compatible with the selected IMethod', &
          & P_MY%IMethod)

   endif ! End of setup mode

!
! These quantities shall be passed at any call
!
! Geometric transforms
!
    P_MY%Sun         = Sun
    P_MY%SpaceCraft  = Spacecraft
    P_MY%HSpaceCraft = Spacecraft - Sun
    P_MY%SC_Center   = SpaceCraft - Sun - P_MY_DensDstr_COBE_TBLOB%Center

!
! Pointing
!
    P_MY%Pointing    = Pointing

!
! The Interface Ends, The serious game begins
!
! Perform the integrations
!

    ! :Step: Integrates over the dust density scaled by BB in RJ approx
    if (P_MY%IMETHOD .eq. 0) then ! Uses QSIMP Integration
       result=simpson(P_FX_COBE_TBLOB_DENS,P_MY%MinDistance,P_MY%MaxDistance,P_MY%Eps,P_MY%MaxNLoops)
    else                     ! Uses QROMB Integration
       result=romberg(P_FX_COBE_TBLOB_DENS,P_MY%MinDistance,P_MY%MaxDistance,P_MY%EPS,P_MY%MaxNLoops)
    endif
     This%Dens_BB       = Result * 1.0

!print*,'END OF COBE_T_TBLOB_DENS'

   ! :Beware: 23 Aug 2004 : M. Maris
   !     If present THIS_ERR passes a dummy value
   !     This is required by some compilers
   if (present(This_Err)) This_Err%Dens_BB = -1

CONTAINS
Subroutine P_DecrementCounter()
!
! Used to decrement the NFILLS counter
!
  if (NFills.gt.0) NFills = NFills-1
End Subroutine P_DecrementCounter

END SUBROUTINE ZOD_BI_FAST_CC_TBLOB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MISCELLANEOUS SERVICES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZOD_BRINT_FAST_STARTUP()
!
! Intialization of ZOD_BRINT_FAST
!
! To be runned only once at start time
!

  print*
  print*,'ZOD_BRINT_FAST Startup Successfull'
  print*
end subroutine ZOD_BRINT_FAST_STARTUP

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PRIVATE METHODS !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Those are methods used to support the integration of
! the Zodiacal light density along the line of view
! with various scalings
!

Double Precision function P_FX_COBE_T(Distance)
!
! for an heliocentric DISTANCE in AU computes the local T of Dust
!
! Parameters are passed through Private Members
!

    Double Precision, Intent(IN) :: Distance

! local variables

    Double Precision, Dimension(1:3) :: HPos, P
    Double Precision :: R, TR

    ! :Step: Computes the heliocentric position for the given pointing
    P = Distance * P_MY%Pointing
    HPos = P + P_MY%HSpaceCraft

    ! :Step: Computes the dust temperature
    R = Dot_Product(HPos,HPos)
    R = sqrt(R)
    TR = TDust(R=R,DTPar=P_MY%DustTemperaturePar)

    P_FX_COBE_T = TR
End Function P_FX_COBE_T


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MISCELLANEOUS SERVICES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! DEFAULT OBJECTS HANDLERS !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE P_ZOD_BRINT_FAST_NEW_BASE(This)
!
! The creator for the basic ZOD_BRINT_FAST object
!
! At creation time all the integrated quantities are set to -1
!
IMPLICIT NONE
  Type (T_BRINT_FAST_COMPONENT), Intent(OUT) :: This

  This%Dens_BB  = -1.d0
END SUBROUTINE P_ZOD_BRINT_FAST_NEW_BASE

SUBROUTINE P_ZOD_BRINT_FAST_DESTROY_BASE(This)
!
! The DESTRUCTOR for the basic ZOD_BRINT_FAST object
!
IMPLICIT NONE
  Type (T_BRINT_FAST_COMPONENT), Intent(OUT) :: This

  call P_ZOD_BRINT_FAST_NEW_BASE(This)

END SUBROUTINE P_ZOD_BRINT_FAST_DESTROY_BASE

SUBROUTINE P_ZOD_BRINT_FAST_SHOW_BASE(This,Unit)
!
! The displayer for the basic ZOD_BRINT_FAST object
!
IMPLICIT NONE
  Type (T_BRINT_FAST_COMPONENT), Intent(IN) :: This
  Integer, Optional,        Intent(IN) :: Unit

!  Double Precision, Allocatable, Dimension(:) :: A
  Integer NCOLS
  Character(len=20) :: VCs1,vcs2

  IF (Present(UNIT)) THEN
     ! Creates the formatting string
     NCOLS = 1
     write(vcs1,*) NCOLS
     vcs1 = adjustl(vcs1)
     write(vcs2,*) '(',trim(vcs1),'g20.8)'
     write(UNIT,vcs2) This%Dens_BB
  ELSE
     print*,This%Dens_BB
  ENDIF

END SUBROUTINE P_ZOD_BRINT_FAST_SHOW_BASE


END MODULE ZOD_BRINT_FAST
