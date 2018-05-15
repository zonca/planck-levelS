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

module beam_alm

  ! Beam multipole module.

  use planck_config
  use ls_misc_utils
  use dmc_io

  use beam_bessel

  implicit none

  real(dp), parameter, private :: &
      bmln2 = 0.693147180559945309417232121458176568_dp

  ! Beam multipoles type.
  !
  ! Internally, this module uses the conventions for G and C
  ! multipoles from Challinor et al. (astro-ph/0008228).  The
  ! normalisation and signs are changed in output routine bm_alm_write
  ! to be compatible with the new conventions used in HEALPix 1.2.

  type bmalm
    integer :: lmax, mmax, nstokes
    complex(dp), dimension(:,:,:), pointer :: a
  end type bmalm

contains

  !======================================================================

  ! Initialise a bmalm type with lmax, mmax, and nstokes Stokes
  ! parameters.

  subroutine bm_alm_init(alm, lmax, mmax, nstokes)
    type(bmalm), intent(inout) :: alm
    integer, intent(in) :: lmax, mmax, nstokes

    integer :: mmax_temp

    mmax_temp = mmax
    if (lmax < 0) then
      call exit_with_status(1, 'bm_alm_init: lmax < 0')
    else if (mmax < 0) then
      call exit_with_status(1, 'bm_alm_init: mmax < 0')
    else if (mmax > lmax) then
      write(*,'("Warning from bm_alm_init:&
          & mmax > lmax: setting mmax = lmax")')
      mmax_temp = lmax
    end if

    if (.not.(nstokes == 1 .or. nstokes == 3 .or. nstokes == 4 )) then
      call exit_with_status(1, 'bm_alm_init: nstokes must be  1, 3 or 4')
    end if

    alm%lmax = lmax
    alm%mmax = mmax_temp
    alm%nstokes = nstokes
    allocate(alm%a(nstokes, 0:lmax, 0:mmax_temp))
    alm%a = (0.0_dp, 0.0_dp)

  end subroutine bm_alm_init

  !======================================================================

  ! Free bmalm type.

  subroutine bm_alm_free(alm)
    type(bmalm), intent(inout) :: alm

    alm%lmax = 0
    alm%mmax = 0
    alm%nstokes = 0
    deallocate(alm%a)

  end subroutine bm_alm_free

  !======================================================================

  ! Calculate the multipoles of a normalised elliptical gaussian beam
  ! pointing in the z-direction. If there are 3 (IGC) or 4 (IGCV)
  ! Stokes parameters, it will create a linearly polarised beam.
  !
  ! The parameters are -
  !
  ! fwhm: full width at half-maximum / radians.
  ! ell: ellipticity (max fwhm/min fwhm)
  ! psi_ell: orientation of beam ellipse major axis wrt x-axis /
  !          radians
  ! psi_pol: orientation of electric polariation measurement wrt
  !          x-axis / radians
  !
  ! This routine assumes that fwhm is small.

  subroutine bm_alm_gaussbeam(alm, lmax, mmax, nstokes, fwhm, ell, &
      psi_ell, psi_pol, epsilon)
    type(bmalm), intent(inout) :: alm
    integer, intent(in) :: lmax, mmax, nstokes
    real(dp), intent(in) :: fwhm, ell, psi_ell, psi_pol, epsilon

    integer :: l, m
    real(dp) :: sigmasq, sigmaxsq, esq, tmp, tmp2, rho
    complex(dp) :: val, f1, f2, b1, b2
    logical :: polarised

    call bm_alm_init(alm, lmax, mmax, nstokes)

    if (psi_pol > 1.0e10) then
      ! Beam is unpolarised, set G and C multipoles to zero.
      polarised = .false.
    else
      polarised = .true.
    end if

    ! Work out whether to calculate elliptical expansion.

    if (ell == 1.0_dp) then

      ! Beam is circular

      sigmasq = fwhm*fwhm/(8.0_dp * bmln2)

      ! Calculate I component

      do l = 0, lmax
        alm%a(1, l, 0) = cmplx(sqrt(real(2*l+1, dp) / fourpi) * &
            exp(-0.5_dp*sigmasq*real(l*l, dp)), 0.0_dp)
      end do

      if (nstokes > 1 .and. polarised) then

        ! Do the G and C components

        f1 = cmplx(cos(2.0_dp*psi_pol), -sin(2.0_dp*psi_pol))

        do l = 2, lmax
          val = cmplx(sqrt(real(2*l+1, dp) / (8.0_dp*fourpi)) * &
              exp(-0.5_dp*sigmasq*real(l*l, dp)), 0.0_dp) * f1
          alm%a(2,l,2) = val
          alm%a(3,l,2) = val * (0.0_dp, 1.0_dp)
        end do
      end if

    else

      ! Beam is elliptical
      !
      ! Simulate beam with major axis of ellipse aligned with x-axis
      ! and polarisation direction aligned at angle (psi_pol - psi_ell)
      ! to x-axis.

      esq = 1.0_dp-1.0_dp/(ell**2)
      sigmaxsq = fwhm*fwhm*ell/(bmln2 * 8.0_dp)

      ! I component

      do l = 0, lmax
        tmp = real(l*l, dp)*sigmaxsq
        do m = 0, min(l, mmax), 2
          alm%a(1, l, m) = cmplx(sqrt(real(2*l+1, dp) / fourpi) * &
              exp(-0.5_dp*tmp*(1.0_dp-0.5_dp*esq)) * &
              bessel_i(m/2, 0.25_dp*tmp*esq), 0.0_dp)
        end do
      end do

      if (nstokes > 1 .and. polarised) then

        ! Do G and C components

        rho = psi_pol - psi_ell
        f1 = cmplx(cos(2.0_dp*rho), -sin(2.0_dp*rho))
        f2 = cmplx(cos(2.0_dp*rho), sin(2.0_dp*rho))

        do l = 2, lmax
          tmp = real(l*l, dp)*sigmaxsq
          tmp2 = 0.25_dp*tmp*esq

          ! m = 0

          val = sqrt(real(2*l+1, dp) / (2.0_dp*fourpi)) * &
              exp(-0.5_dp*tmp*(1.0_dp-0.5_dp*esq)) * &
              bessel_i(1, tmp2)

          alm%a(2, l, 0) = val*cos(2.0_dp*rho)
          alm%a(3, l, 0) = val*sin(2.0_dp*rho)

          ! m = 2, 4, and so on.

          do m = 2, min(l, mmax), 2
            val = sqrt(real(2*l+1, dp) / (8.0_dp*fourpi)) * &
                exp(-0.5_dp*tmp*(1.0_dp-0.5_dp*esq))
            b1 = f1*bessel_i((m-2)/2, tmp2)
            b2 = f2*bessel_i((m+2)/2, tmp2)
            alm%a(2,l,m) = val * (b1+b2)
            alm%a(3,l,m) = val * (b1-b2) * (0.0_dp, 1.0_dp)
          end do
        end do

      end if

      ! Rotate multipoles through angle psi_ell about z-axis, so the
      ! beam is in the right orientation (only need do this for even m).

      do m = 0, mmax, 2
        f1 = cmplx(cos(real(m, dp)*psi_ell), -sin(real(m, dp)*psi_ell))
        alm%a(:, :, m) = alm%a(:, :, m) * f1
      end do

    end if

    ! Adjust multipoles for cross-polar leakage.

    alm%a(1, :, :) = alm%a(1, :, :) * 0.5_dp * (1.0_dp + epsilon)

    if (nstokes > 1)  then
      alm%a(2:3, :, :) = alm%a(2:3, :, :) * 0.5_dp * (1.0_dp - epsilon)
    end if

    ! Adjust normalisation.

    if (nstokes > 1)  then
      alm%a(2:3, :, :) = -alm%a(2:3, :, :) * sqrt2
    end if

  end subroutine bm_alm_gaussbeam

  !======================================================================

  ! Normalise a set of multipoles.
  !
  ! This scales the harmonics so that a^I_00 = 1/sqrt(4pi) if a^I_00
  ! is non-zero.

  !subroutine bm_alm_normalise(alm)
  !  type(bmalm), intent(inout) :: alm
  !
  !  real(sp) :: norm
  !
  !  norm = real(alm%a(1,0,0), sp)*sqrt(fourpi)
  !  if (norm > tiny(0.0_sp)) alm%a = alm%a/norm
  !
  !end subroutine bm_alm_normalise

  !======================================================================

  ! Write out a set of multipoles (HEALPix-style) to a DMC object.

  subroutine bm_alm_write(alm, filename)
    type(bmalm), intent(in) :: alm
    character(len=*), intent(in) :: filename

    integer :: l, m, nalm, pol
    real(dp), pointer :: realtmp(:), imagtmp(:)
    integer(i4b), pointer :: lmtmp(:)
    character (len=1) :: suffix
    type(dmc_handle) :: out

    ! First write the a_lms from the data structure into the temporary
    ! array.

    if (alm%nstokes == 1) then
      call dmc_create(out, filename, 'alm.LS_alm')
    else if (alm%nstokes == 3) then
      call dmc_create(out, filename, 'alm.LS_alm_pol')
    else
      call exit_with_status(1, 'bad a_lm format for DMC')
    endif

    do pol = 1, alm%nstokes
      select case (pol)
        case (1)
          suffix='T'
        case (2)
          suffix='G'
        case (3)
          suffix='C'
      end select

      if (pol==1) then
        call dmc_set_key(out, 'lmaxT', alm%lmax)
        call dmc_set_key(out, 'mmaxT', alm%mmax)
        if (alm%nstokes>=3) then
          call dmc_set_key(out, 'lmaxG', alm%lmax)
          call dmc_set_key(out, 'mmaxG', alm%mmax)
          call dmc_set_key(out, 'lmaxC', alm%lmax)
          call dmc_set_key(out, 'mmaxC', alm%mmax)
        endif
      endif

      nalm = ((alm%mmax+1)*(alm%mmax+2))/2 + (alm%mmax+1)*(alm%lmax-alm%mmax)

      allocate(lmtmp(nalm), realtmp(nalm), imagtmp(nalm))

      nalm = 0
      do l = 0,alm%lmax
        do m = 0,min(l,alm%mmax)
          nalm = nalm+1
          lmtmp(nalm) = l*l+l+m+1
          realtmp(nalm) = real (alm%a(pol, l, m))
          imagtmp(nalm) = aimag(alm%a(pol, l, m))
        end do
      end do

      call dmc_append_column(out, dmc_colnum(out, 'Index'//suffix), lmtmp)
      call dmc_append_column(out, dmc_colnum(out, 'Real'//suffix), realtmp)
      call dmc_append_column(out, dmc_colnum(out, 'Imag'//suffix), imagtmp)

      deallocate(lmtmp, realtmp, imagtmp)

    end do

    call dmc_close(out)

  end subroutine bm_alm_write

  !======================================================================

end module beam_alm
