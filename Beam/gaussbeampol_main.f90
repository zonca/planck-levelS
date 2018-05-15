!-----------------------------------------------------------------------------
!
!  Copyright (C) 2002-2013 Mark Ashdown
!
!  This file is part of the "Beam" component of the Planck simulation
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

module gaussbeampol_main

  implicit none

contains

  !======================================================================

  ! Calculates the multipoles of an linearly-polarised elliptical
  ! Gaussian beam with following properties -
  !
  ! - beam points in z-direction (towards "north pole").
  !
  ! - beam ellipse is described by three parameters:
  !   1) mean FWHM, defined as fwhm = sqrt(fwhm_max*fwhm_min);
  !   2) ellipticity, defined as fwhm_max/fwhm_min;
  !   3) orientation angle. Major axis is orinted at angle psi_ell to
  !      the x-axis;
  !
  ! - beam is polarised along direction at angle psi_pol to the
  !   x-axis.
  !
  ! Mark Ashdown, CPAC.

  function gaussbeampol_module2 (args)

    use planck_config
    use ls_misc_utils
    use announce_mod
    use ls_paramfile_io
    use dmc_io
    use focalplanemod
    use beam
    implicit none

    character(len=*), intent(in) :: args(:)
    integer(i4b) :: gaussbeampol_module2
    type(paramfile_handle) :: params
    character(len=filenamelen) :: detdbfile, almfile
    character(len=9) :: detid
    real(dp) :: fwhm, ell, psi_ell, psi_pol, epsilon
    integer :: lmax, mmax, nstokes
    logical :: elliptic
    type(bmalm) :: blm

    ! Announce the program

    call module_startup('gaussbeampol',args,.true.)

    ! Get name of parameter file

    call dmc_init(args)
    params = dmc_getrunparams()

    detdbfile = parse_string(params, 'focalplane_db', &
        descr='focal plane database file')
    detid = trim(parse_string(params, 'detector_id', descr='detector ID'))
    lmax = parse_int(params, 'beam_lmax', default=1024, vmin=2, &
        descr='maximum l for analysis')
    mmax = parse_int(params, 'beam_mmax', default=2, vmin=2, &
        descr='maximum m for analysis')
    elliptic = parse_lgt(params, 'beam_elliptic', default=.true., &
        descr='should an elliptic beam be modelled?')
    nstokes = parse_int(params, 'beam_nstokes', default=3, vmin=1, vmax=4, &
        descr='number of Stokes parameters (1 = Intensity only)')
    almfile = parse_string(params, 'beam_alm_file', &
        descr='output beam multipole file name')

    ! Read detector database

    call fpdb_open(detdbfile)

    ! Get beam parameters from detector database

    fwhm = fpdb_value_real8(detid, 'beamfwhm')
    if (elliptic) then
      ell = fpdb_value_real8(detid, 'ellipticity')
      psi_ell = fpdb_value_real8(detid, 'psi_ell')
    else
      ell = 1.0_dp
      psi_ell = 0.0_dp
    endif
    psi_pol = fpdb_value_real8(detid, 'psi_pol')
    epsilon = fpdb_value_real8(detid, 'epsilon')

    ! Close the detector database

    call fpdb_close

    call parse_finish(params)

    ! Convert fwhm and angles from degrees to radians

    fwhm = fwhm*pi/180.0_dp
    psi_ell = psi_ell*pi/180.0_dp
    psi_pol = psi_pol*pi/180.0_dp

    write(*,'("")')

    ! Allocate multipoles and create gaussian beam.

    write(*,'("Calculating beam multipoles")')

    call bm_alm_gaussbeam(blm, lmax, mmax, nstokes, fwhm, ell, psi_ell, &
        psi_pol, epsilon)

    ! Write multipoles to FITS file.

    write(*,'("Writing multipoles to ",a)') trim(almfile)

    call bm_alm_write(blm, almfile)

    ! Deallocate multipole type.

    call bm_alm_free(blm)

    ! Shutdown the DMC.

    call dmc_shutdown

    ! Set return value.

    gaussbeampol_module2 = 0

  end function gaussbeampol_module2

  function gaussbeampol_module (init_object)
    use planck_config

    character(len=*), intent(in) :: init_object
    integer(i4b) :: gaussbeampol_module
    gaussbeampol_module = gaussbeampol_module2 ( (/init_object/) )
  end function gaussbeampol_module

  !======================================================================

end module gaussbeampol_main
