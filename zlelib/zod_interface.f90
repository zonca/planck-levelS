!-----------------------------------------------------------------------------
!
!  Copyright (C) 2010-2013 Michele Maris
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
module zle_interface
use planck_config
use ZOD_dd_cobe_smooth
use ZOD_dd_cobe_band
use ZOD_dd_cobe_cring
use ZOD_dd_cobe_tblob
use ZOD_brint_fast
use ZOD_ColorCorrection
implicit none
public

integer,allocatable, save :: comp(:) ! set of components to compute
real(dp),allocatable, save :: fact(:) ! factors for individual components
real(dp), save :: frequencyGHz ! The frequency in GHz
real(dp), save :: integration_mindistance, integration_maxdistance ! in AU
integer, save :: Integration_MaxNLoops ! Maximum number of loops in integration
real(dp), save :: Integration_Accuracy  ! Integration accuracy (>0: absolute, <0: relative)
integer, save :: Integration_Method    ! Integration Method 0, QSIMP; 1, QROMB

Type (T_DustTemperature), save :: DustTemperaturePar

Type (T_DensDstr_COBE_SMOOTH), save :: DensDstr_COBE_SMOOTH(2)
Type (T_DensDstr_COBE_BAND), save :: DensDstr_COBE_BAND(3)
Type (T_DensDstr_COBE_CRING), save :: DensDstr_COBE_CRING
Type (T_DensDstr_COBE_TBLOB), save :: DensDstr_COBE_TBLOB

Type (T_ColorCorrection_COBE), save :: ColorCorrection_COBE

end module zle_interface

subroutine zle_init(freq,component,fct,ncomp,accuracy,nloop,mindist,maxdist) &
  bind(c)
use iso_c_binding
use zle_interface
use ZOD_brint_fast
use ZOD_densdstr
implicit none
real(c_double), intent(in), value :: freq
integer(c_int), intent(in) :: component(*)
real(c_double), intent(in) :: fct(*)
integer(c_int), intent(in), value :: ncomp
real(c_double), intent(in), value :: accuracy
integer(c_int), intent(in), value :: nloop
real(c_double), intent(in), value :: mindist,maxdist

integer i, icomponent
character(len=70) :: compname

Call ZOD_DENSDSTR_STARTUP
!Call ZOD_BRINT_FAST_STARTUP ! No need to call, only prints a banner

call ZOD_DustTemperature_New(DustTemperaturePar)
call ColorCorrection_New(ColorCorrection_COBE)

frequencyGHz=freq
allocate(comp(ncomp),fact(ncomp))
comp=component(1:ncomp)
fact=fct(1:ncomp)
Integration_Accuracy=accuracy
Integration_MaxNLoops=nloop
integration_mindistance=mindist
integration_maxdistance=maxdist

do i = 1,size(comp) ! loop over all active components
  icomponent=comp(i)

  ComponentSelection: SELECT CASE (ABS(ICOMPONENT))
  case (001, 002) ! SMOOTH COMPONENT
    call ZOD_DensDstr_New(DensDstr_COBE_SMOOTH(icomponent), &
      ICOMPONENT=ICOMPONENT,CNAME=compname)

    call ZOD_DustTemperature_Set(DustTemperaturePar, &
      T0=DensDstr_COBE_SMOOTH(icomponent)%T0,SigmaT0 = 0.d0, &
      Delta=DensDstr_COBE_SMOOTH(icomponent)%Delta,SigmaDelta = 0.d0)

  case (011, 012, 013)  ! BANDS
    call ZOD_DensDstr_New(DensDstr_COBE_BAND(icomponent-10), &
      IBAND=ICOMPONENT-10,CNAME=compname)

    call ZOD_DustTemperature_Set(DustTemperaturePar, &
      T0=DensDstr_COBE_BAND(icomponent-10)%T0,SigmaT0 = 0.d0, &
      Delta=DensDstr_COBE_BAND(icomponent-10)%Delta,SigmaDelta = 0.d0)

  case (020)  ! CRING COMPONENT
    call ZOD_DensDstr_New(DensDstr_COBE_CRING,CNAME=compname)

    call ZOD_DustTemperature_Set(DustTemperaturePar, &
      T0=DensDstr_COBE_CRING%T0,SigmaT0 = 0.d0, &
      Delta=DensDstr_COBE_CRING%Delta,SigmaDelta = 0.d0)

  case (030)  ! TBLOB COMPONENT
    call ZOD_DensDstr_New(DensDstr_COBE_TBLOB,CNAME=compname)

    call ZOD_DustTemperature_Set(DustTemperaturePar, &
      T0=DensDstr_COBE_TBLOB%T0,SigmaT0 = 0.d0, &
      Delta=DensDstr_COBE_TBLOB%Delta,SigmaDelta = 0.d0)

  CASE DEFAULT  ! Unrecongnized Component
    call exit_with_status(1,"unknown component")
  END SELECT ComponentSelection

    ! :Step: Shows the defined/readed component
!  ComponentShowing: SELECT CASE (abs(ICOMPONENT))
!    case (001, 002      )  ! SMOOTH COMPONENT
!        CALL ZOD_DensDstr_COBE_SMOOTH_SHP(DensDstr_COBE_SMOOTH(icomponent))
!    case (011, 012, 013)  ! BANDS
!        CALL ZOD_DensDstr_COBE_BAND_SHP(DensDstr_COBE_BAND(icomponent-10))
!    case (020            )  ! CRING COMPONENT
!        CALL ZOD_DensDstr_COBE_CRING_SHP(DensDstr_COBE_CRING)
!    case (030            )  ! TBLOB COMPONENT
!        CALL ZOD_DensDstr_COBE_TBLOB_SHP(DensDstr_COBE_TBLOB)
!  END SELECT ComponentShowing
end do

call ZOD_DustTemperature_Show(DustTemperaturePar)
end subroutine zle_init

subroutine zle_compute (sun,bary,nsamples,directions,signal) bind(c)
use iso_c_binding
use zle_interface
use zod_brint_fast
implicit none
real(c_double), intent(in) :: sun(3), bary(3)
integer(c_int), intent(in), value :: nsamples
real(c_double), intent(in) :: directions(*)
real(c_double), intent(out) :: signal(*)

real(dp), parameter :: AU=1.49597870691e11_dp ! astronomical unit in m

integer i,icomponent,ipsi
real(dp) :: m_plk(3), m_sun(3), m_beam(3), antfct

real(dp), parameter :: speedOfLight=2.99792458e8_dp, &
                       kBoltzmann=1.3806504e-23_dp

Type (T_BRINT_FAST_COMPONENT) :: zbf_comp

CALL ZOD_BRINT_FAST_NEW(zbf_comp)

m_plk = -bary/AU
m_sun = (sun-bary)/AU

antfct= speedOfLight*speedOfLight &
        /(2*frequencyGHz*frequencyGHz*1e18_dp*kBoltzmann)*1e-20;

do i=1,size(comp)
  icomponent=comp(i)

  InitCmpSel: SELECT CASE (abs(ICOMPONENT))
  CASE (001,002)  ! SMOOTH COMPONENT
    Call ZOD_BI_FAST_CC_SMOOTH(SetUp=.true.,FirstCall= .true., &
      ParMod = DensDstr_COBE_SMOOTH(icomponent), &
      DustTemperaturePar = DustTemperaturePar, &
      FrequencyGHz = FrequencyGHz,MinDistance = Integration_MinDistance, &
      MaxDistance = Integration_MaxDistance, IMETHOD = Integration_Method, &
      Eps = Integration_Accuracy,MaxNLoops = Integration_MaxNLoops, &
      THIS=zbf_comp, Pointing=M_BEAM,SpaceCraft=M_PLK, &
      Sun=M_Sun)

    loop_ipsi_smooth: do ipsi=0,NSamples-1
      M_BEAM    = directions(3*ipsi+1:3*ipsi+3)
      Call ZOD_BI_FAST_CC_SMOOTH(THIS=zbf_comp, &
        Pointing=M_BEAM,SpaceCraft=M_PLK,Sun=M_Sun)
      signal(ipsi+1)=signal(ipsi+1)+fact(i)*antfct*zbf_comp%Dens_BB
    enddo loop_ipsi_smooth

  CASE (011, 012, 013)  ! BANDS
    Call ZOD_BI_FAST_CC_BAND(SetUp=.true.,FirstCall= .true., &
      ParMod = DensDstr_COBE_BAND(icomponent-10), &
      DustTemperaturePar = DustTemperaturePar, &
      FrequencyGHz = FrequencyGHz,MinDistance = Integration_MinDistance, &
      MaxDistance = Integration_MaxDistance, IMETHOD = Integration_Method, &
      Eps = Integration_Accuracy,MaxNLoops = Integration_MaxNLoops, &
      THIS=zbf_comp, Pointing=M_BEAM,SpaceCraft=M_PLK, &
      Sun=M_Sun)

    loop_ipsi_band: do ipsi=0,NSamples-1
      M_BEAM    = directions(3*ipsi+1:3*ipsi+3)
      Call ZOD_BI_FAST_CC_BAND(THIS=zbf_comp, &
        Pointing=M_BEAM,SpaceCraft=M_PLK,Sun=M_Sun)
      signal(ipsi+1)=signal(ipsi+1)+fact(i)*antfct*zbf_comp%Dens_BB
    enddo loop_ipsi_band

  CASE (020) ! CRING
    Call ZOD_BI_FAST_CC_CRING(SetUp=.true.,FirstCall= .true., &
      ParMod = DensDstr_COBE_CRING, &
      DustTemperaturePar = DustTemperaturePar, &
      FrequencyGHz = FrequencyGHz,MinDistance = Integration_MinDistance, &
      MaxDistance = Integration_MaxDistance, IMETHOD = Integration_Method, &
      Eps = Integration_Accuracy,MaxNLoops = Integration_MaxNLoops, &
      THIS=zbf_comp, Pointing=M_BEAM,SpaceCraft=M_PLK, &
      Sun=M_Sun)
    loop_ipsi_cring: do ipsi=0,NSamples-1
      M_BEAM    = directions(3*ipsi+1:3*ipsi+3)
      Call ZOD_BI_FAST_CC_CRING(THIS=zbf_comp, &
        Pointing=M_BEAM,SpaceCraft=M_PLK,Sun=M_Sun)
      signal(ipsi+1)=signal(ipsi+1)+fact(i)*antfct*zbf_comp%Dens_BB
    enddo loop_ipsi_cring

  CASE (030) ! TBLOB
    Call ZOD_BI_FAST_CC_TBLOB(SetUp=.true.,FirstCall= .true., &
      ParMod = DensDstr_COBE_TBLOB, &
      DustTemperaturePar = DustTemperaturePar, &
      FrequencyGHz = FrequencyGHz,MinDistance = Integration_MinDistance, &
      MaxDistance = Integration_MaxDistance, IMETHOD = Integration_Method, &
      Eps = Integration_Accuracy,MaxNLoops = Integration_MaxNLoops, &
      THIS=zbf_comp, Pointing=M_BEAM,SpaceCraft=M_PLK, &
      Sun=M_Sun)

    loop_ipsi_tblob: do ipsi=0,NSamples-1
      M_BEAM    = directions(3*ipsi+1:3*ipsi+3)
      Call ZOD_BI_FAST_CC_TBLOB(THIS=Zbf_comp, &
        Pointing=M_BEAM,SpaceCraft=M_PLK,Sun=M_Sun)
      signal(ipsi+1)=signal(ipsi+1)+fact(i)*antfct*zbf_comp%Dens_BB
    enddo loop_ipsi_tblob

  END SELECT InitCmpSel

end do

CALL ZOD_BRINT_FAST_DESTROY(zbf_comp)

end subroutine zle_compute

subroutine zle_destroy() bind(c)
use iso_c_binding
use zle_interface
implicit none
deallocate(comp)
CALL ZOD_DustTemperature_DESTROY(DustTemperaturePar)
end subroutine zle_destroy
