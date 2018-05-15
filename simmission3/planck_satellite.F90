!-----------------------------------------------------------------------------
!
!  Copyright (C) 1999-2013 Daniel Mortlock
!
!  This file is part of the "simmission" component of the Planck simulation
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

module planck_satellite
  use planck_config
  use general_vector
  use general_matrix
  use ls_paramfile_io
  implicit none
  private

  public :: plsolarpanel, plsatellite, pl_satellite_init

  ! The Solar panel of the Planck satellite, characterised by its radius,
  ! radius (in m), its area, area (in m^2), and its specular and diffuse
  ! reflection coefficients, specref and diffref.
  type plsolarpanel
    real(dp) :: radius
    real(dp) :: area
    real(dp) :: specref
    real(dp) :: diffref
  end type plsolarpanel

  ! Physical model of the Planck satellite. The dynamics keyword indicates
  ! the assumed satellite dynamics: ideal or analytical. The
  ! centre-of-gravity, pos_cog, is given in m relative to the origin of the
  ! satellite reference system coordinate axes (in m). The initial/fiducial
  ! value of the inertia tensor in the satellite reference system, i_srs_0,
  ! is given in kg m^2 and kg m^2 s^(-1). Note that the inertia tensor
  ! (and its derivatives, etc.), being symmetric, are stored as vectors,
  ! rather than matrices, with the ``packing'' defined by:
  !
  !  I(1) = I(1, 1)
  !  I(2) = I(2, 1) = I(1, 2)
  !  I(3) = I(2, 2)
  !  I(4) = I(3, 1) = I(1, 3)
  !  I(5) = I(3, 2) = I(2, 3)
  !  I(6) = I(3, 3)
  !
  ! Also stored are a variety of derived quantities that change over the
  ! duration of the mission: the actual inertia tensor in the satellite
  ! reference system, i_srs (in kg m^2); the inertia tensor in the rotated
  ! frame in which only the x-y off-diagonal couplings are retained,
  ! i_irs (in kg m^2); and a number of other quantities that are reused in the
  ! calculations.
  ! Also, the rotation matrices stored here act on vectors in the IRS
  ! and SRS, respectively, but are defined to rotate frames rather than
  ! rotate the vectors, so should be used in transpose.
  type plsatellite
    character(len=filenamelen) :: dynamics
    real(dp) :: pos_cog(3)
    real(dp) :: i_srs_0(6)
    real(dp) :: i_irs(6)
    real(dp) :: rotm_srs2irs(3, 3)
    type(plsolarpanel) :: solarpanel
  end type plsatellite

contains

  ! Initialise a model of the Planck satellite.
  function pl_satellite_init(params) result(satellite)
    type(paramfile_handle), intent(inout) :: params
    type(plsatellite) :: satellite

    write(*, '(/,a,/)') 'Planck satellite parameters.'

    ! Read in the dynamics.
    satellite%dynamics = parse_string(params,'dynamics')

    ! Read in the coordinates of the centre of gravity.
    write(*, '(a)')
    satellite%pos_cog(1) = parse_double(params,'pos_cog_sat_x')
    satellite%pos_cog(2) = parse_double(params,'pos_cog_sat_y')
    satellite%pos_cog(3) = parse_double(params,'pos_cog_sat_z')

    ! Read in the elements of the inertia tensor, copying the symmetric
    ! parts across the diagonal.
    write(*, '(a)')
    satellite%i_srs_0(1) = parse_double(params,'i_sat_xx', vmin=0.0_dp)
    satellite%i_srs_0(3) = parse_double(params,'i_sat_yy', vmin=0.0_dp)
    satellite%i_srs_0(6) = parse_double(params,'i_sat_zz', vmin=0.0_dp)
    satellite%i_srs_0(2) = parse_double(params,'i_sat_xy')
    satellite%i_srs_0(4) = parse_double(params,'i_sat_xz')
    satellite%i_srs_0(5) = parse_double(params,'i_sat_yz')

    call pl_satellite_t(satellite)

    ! Read in the properties of the satellite's solar panel.
    write(*, '(a)')
    satellite%solarpanel%radius = parse_double(params,'r_panel_sat',vmin=0.0_dp)
    satellite%solarpanel%area = pi * satellite%solarpanel%radius**2

    satellite%solarpanel%specref = parse_double(params,'specref_panel_sat', &
      vmin=0.0_dp, vmax=1.0_dp)
    satellite%solarpanel%diffref = parse_double(params,'diffrel_panel_sat', &
      vmin=0.0_dp, vmax=1.0_dp - satellite%solarpanel%specref)
  end function pl_satellite_init

  ! Calculate the inertia tensor, both in the satellite reference system (SRS)
  ! and the inertial reference system (IRS).
  subroutine pl_satellite_t(satellite)
    type(plsatellite), intent(in out) :: satellite

    real(dp) :: b, db, f, df, i_srs(6), i_irs(6), rotm_srs2irs(3, 3), tolerance
    integer :: i

    ! Calculate the current value of the inertia tensor in the SRS.
    i_srs = satellite%i_srs_0

    ! Iterative evaluation of the inertia tensor in the IRS.
    b = 0.0
    f = 0.0
    i_irs = i_srs
    tolerance = 1e-7_dp

    do i = 1, 5
      db = asin(i_irs(4) / (i_irs(1) - i_irs(6)))
      df = asin(i_irs(5) / (i_irs(6) - i_irs(3)))

      ! If these perturbations are below the tolerance (and we have
      ! already gone through one loop) save time by ceasing the calculation.
      if (i > 1) then
        if ((abs(db)<tolerance) .and. (abs(df)<tolerance)) exit
      end if

      b = b + db
      f = f + df

      ! Build up the rotation matrix, noting that the old calls used the
      ! opposite convention for the direction of the rotation. The old
      ! call used to be different in terms of the angle defintions:
      !
      ! am = gn_identitym(3)
      ! call dn_mrot3(2, b, am, an)
      ! call dn_mrot3(1, f, an, am)
      rotm_srs2irs = gn_rotm(- b, 2)
      rotm_srs2irs = gn_rotm(- f, 1, rotm_srs2irs)

      ! Apply the current best guess rotation and iterate.
      i_irs = gn_jacobian_sm(rotm_srs2irs, i_srs)
    end do

    ! Copy the newly calculated matrices/tensors into the satellite
    ! structure.
    satellite%rotm_srs2irs = rotm_srs2irs
    satellite%i_irs = i_irs
  end subroutine pl_satellite_t

end module planck_satellite
