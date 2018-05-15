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

module beam2alm_main

  implicit none

contains

  !======================================================================

  ! Reads beam Stokes parameters on a square or polar grid and
  ! calculates their spherical harmonic coefficients by quadrature.
  ! If the input beam is on a square grid, it is interpolated onto a
  ! polar grid before the analysis.
  !
  ! This program can take its input from two files, one containing the
  ! main beam at high resolution, the other the full-sky beam at lower
  ! resolution. In the main beam is present, the full beam is
  ! interpolated in the theta direction so that its resolution is the
  ! same as the main beam.
  !
  ! Mark Ashdown, CPAC.

  function beam2alm_module2(args)

    use planck_config
    use ls_misc_utils
    use announce_mod
    use ls_paramfile_io
    use dmc_io
    use minihealpix
    use beam

    character(len=*), intent(in) :: args(:)
    integer(i4b) :: beam2alm_module2
    logical :: main_present, full_present
    integer :: lmax, mmax, ntheta, nphi
    real(dp) :: tmin, tmax, dtheta
    character(len=20) :: main_format, full_format
    character(len=filenamelen) :: main_file, full_file, alm_file
    type(paramfile_handle) :: params
    type(bmsquare) :: main_s
    type(bmpolar) :: main, full
    type(bmalm) :: alm
    real(dp), allocatable :: map(:,:)
    character(len=filenamelen) :: healpix_file
    type(dmc_handle) :: out
    real(dp) :: theta,phi
    integer(i4b) :: nside,pix,x

    ! Announce the program.

    call module_startup('beam2alm',args,.true.)

    ! Initialise DMC and get parameters.

    call dmc_init(args)

    params = dmc_getrunparams()

    main_format = ' '
    if (parse_string(params, 'beam_main_file_polar', &
        default='', descr='polar main beam file') /= '') then
      main_format = 'polar'
    elseif (parse_string(params, 'beam_main_file_square', &
        default='', descr='square main beam file') /= '') then
      main_format = 'square'
    endif

    full_format = ' '
    if (parse_string(params, 'beam_full_file', &
        default='', descr='polar full beam file') /= '') then
      full_format = 'polar'
    endif

    ! Determine the type of input data

    main_present = .false.
    select case (main_format)
    case('square', 'polar')
      ! There is a main beam.
      main_present = .true.
    case(' ')
      ! No main beam, do nothing
    case default
      ! Print error message
      call exit_with_status(1, 'Unknown format for main beam')
    end select

    full_present = .false.
    select case (full_format)
    case('polar')
      ! There is a full-sky beam.
      full_present = .true.
    case(' ')
      ! No full beam, do nothing
    case default
      ! Print error message
      call exit_with_status(1, 'Unknown format for full-sky beam')
    end select

    if ((.not. main_present) .and. (.not. full_present)) then
      call exit_with_status(1, &
          'There must be at least one input beam data file')
    end if

    if (main_present) then
      main_file = parse_string(params, 'beam_main_file_'//main_format, &
          descr='main beam file name')
    end if

    if (full_present) then
      full_file = parse_string(params, 'beam_full_file', &
          descr='full-sky beam file name')
    end if

    alm_file = parse_string(params, 'beam_alm_file', &
        descr='output beam multipole file name')
    lmax = parse_int(params, 'beam_lmax', vmin=0, &
        descr='maximum l for analysis')
    mmax = parse_int(params, 'beam_mmax', vmin=0, vmax=lmax, &
        descr='maximum m for analysis')

    if (main_format == 'square') then

      ntheta = parse_int(params, 'beam_ntheta', vmin=2, &
          descr='number of theta values for main beam interpolation')
      nphi = parse_int(params, 'beam_nphi', vmin=1, &
          descr='number of phi values for main beam interpolation')

    end if

    healpix_file = parse_string(params, 'beam_healpix_file', &
        descr='optional Healpix map file name', default='')
    if (healpix_file/='') &
      nside = parse_int(params, 'beam_nside', descr='Nside parameter')

    call parse_finish(params)

    write(*, '("")')

    ! Read beam data from input files.

    if (main_present) then

      ! There is a main beam to be analysed.

      select case (main_format)
      case ('polar')

        write(*, '("Reading main beam from ",a)') trim(main_file)
        call bm_polar_read(main, main_file)

      case ('square')

        write(*, '("Reading main beam from ",a)') trim(main_file)

        call bm_square_read(main_s, main_file)

        write(*, '("Interpolating main beam onto polar grid")')

        call bm_square2polar(main_s, main, nphi, ntheta)
        call bm_square_free(main_s)

      end select

    end if

    if (full_present) then

      write(*, '("Reading full-sky beam from ",a)') trim(full_file)

      call bm_polar_read(full, full_file)

      if (main_present) then

        ! If main beam is also present then interpolate full sky
        ! beam values to same theta resolution.

        write(*, '("Interpolating full-sky beam values")')

        ! Work out theta_max, theta_min and ntheta values for the
        ! interpolation.

        dtheta = ((main%theta_max-main%theta_min) / &
            real(main%ntheta-1,dp))
        tmin = main%theta_max + dtheta
        tmax = pi
        ntheta = nint((tmax-tmin)/dtheta)+1

        ! Do the interpolation.

        call bm_interp_linear(full, tmin, tmax, ntheta)

      end if

    end if

    ! Initialise multipoles.

    call bm_alm_init(alm, lmax, mmax, 3)

    ! Perform the transform.

    write(*, '("Calculating multipoles up to lmax = ",i5,", mmax = ",i5)') &
        lmax, mmax

    ! Rotate beam Stokes parameters to polar basis.

    if (main_present) call bm_polar_stokesrotate(main)
    if (full_present) call bm_polar_stokesrotate(full)

    if (healpix_file/='') then
      allocate(map(12*nside*nside,3))
      map = -1.6375e30
      do pix=1,12*nside*nside
        call pix2ang_ring (nside,pix-1,theta,phi)
        if (main_present .and. (theta<main%theta_max)) then
          do x=1,3
            map(pix,x) = bm_polar_get_value(main,theta,phi,x)
          end do
        else if (full_present) then
          do x=1,3
            map(pix,x) = bm_polar_get_value(full,theta,phi,x)
          end do
        else
          exit
        endif
      end do

      call dmc_create(out, healpix_file, 'map.LS_map_pol')
      call dmc_set_key(out, 'Nside', nside)
      call dmc_set_key(out, 'Ordering', 'RING')
      call dmc_append_column(out, dmc_colnum(out, 'I_Stokes'), map(:,1))
      call dmc_append_column(out, dmc_colnum(out, 'Q_Stokes'), map(:,2))
      call dmc_append_column(out, dmc_colnum(out, 'U_Stokes'), map(:,3))
      call dmc_close(out)
      deallocate(map)
    endif

    ! Call the correct transform with the correct data types.

    if (main_present) call bm_polar2alm(main, alm)
    if (full_present) call bm_polar2alm(full, alm)

    ! Deallocate the polar beam types.

    if (main_present) call bm_polar_free(main)
    if (full_present) call bm_polar_free(full)

    ! Write multipoles to FITS file.

    write(*, '("Writing beam multipoles to ",a)') trim(alm_file)
    call bm_alm_write(alm, alm_file)

    ! Deallocate the multipoles.

    call bm_alm_free(alm)

    ! Shutdown the DMC.

    call dmc_shutdown

    ! Set return value.

    beam2alm_module2 = 0

contains

  function bm_polar_get_value(beam,theta,phi,x)
    type(bmpolar), intent(in) :: beam
    real(dp), intent(in) :: theta,phi
    integer, intent(in) :: x
    real(dp) bm_polar_get_value

    real(dp) dth,dph,wth,wph,th1,ph1
    integer(i4b) ith1,ith2,iph1,iph2

    dth=(beam%theta_max-beam%theta_min)/(beam%ntheta-1)
    dph=twopi/beam%nphi

    ith1=1+int((theta-beam%theta_min)/dth)
    ith1=max(1,min(beam%ntheta-1,ith1))
    ith2=ith1+1
    iph1=1+int(phi/dph)
    iph1=max(1,min(beam%nphi,iph1))
    iph2=iph1+1
    if (iph2>beam%nphi) iph2=1

    th1=beam%theta_min+(ith1-1)*dth
    wth=1._dp - (theta-th1)/dth

    ph1=(iph1-1)*dph
    wph=1._dp - (phi-ph1)/dph
    bm_polar_get_value = &
             wth * (wph*beam%stokes(x,iph1,ith1) + (1._dp-wph)*beam%stokes(x,iph2,ith1)) + &
      (1._dp-wth)* (wph*beam%stokes(x,iph1,ith2) + (1._dp-wph)*beam%stokes(x,iph2,ith2))
  end function

  end function beam2alm_module2

  function beam2alm_module (init_object)
    use planck_config

    character(len=*), intent(in) :: init_object
    integer(i4b) :: beam2alm_module
    beam2alm_module = beam2alm_module2 ( (/init_object/) )
  end function beam2alm_module

  !======================================================================

end module beam2alm_main
