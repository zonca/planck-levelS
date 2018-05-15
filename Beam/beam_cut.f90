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

module beam_cut

  use planck_config
  use ls_misc_utils

  implicit none

  ! Beam type corresponding to Grasp "cut" file format.

  type bmcut
    integer :: ncomp, np, ncut
    integer :: icomp, icon
    real(dp) :: sa, da
    real(dp), dimension(:), pointer :: ca
    complex(dp), dimension(:,:,:), pointer :: amp
  end type bmcut

contains

  !======================================================================

  subroutine bm_cut_read(beam, filename)
    type(bmcut), intent(inout) :: beam
    character(len=*), intent(in) :: filename

    integer :: ncut, np, icomp, icon, ncomp
    integer :: icut, i
    real(dp) :: sa, da, ca, tmp1, tmp2, tmp3, tmp4
    character(len=80) :: line

    open(1, file=filename, status='old', action='read')

    ncut = 0

    ! Count the number of cuts.

    do

      ! Skip over header lines until line of correct format is encountered.

      do
        read(1, '(a)', end=4712) line
        read(line, *, err=4711) sa, da, np, ca, icomp, icon, ncomp
        exit
4711    continue
      end do

      ncut = ncut + 1

      call assert(ncomp<=2, 'Three field components present.&
            & Beam package can only handle two field components.')
      call assert(2*(np/2)/=np, 'The number of pixels in a cut must be odd.')

      ! Skip over data.

      do i = 1, np
        read(1, *)
      end do

    end do

4712 rewind(1)

    !print *, 'DEBUG: ', ncut, sa, da, np, icomp, icon, ncomp

    beam%ncut = ncut
    beam%sa = sa
    beam%da = da
    beam%np = np
    beam%icomp = icomp
    beam%icon = icon
    beam%ncomp = ncomp
    allocate(beam%ca(ncut))
    allocate(beam%amp(ncomp, np, ncut))

    do icut = 1, ncut

      ! Skip over header.

      do
        read(1, '(a)') line
        read(line, *, err=4713) sa, da, np, ca, icomp, icon, ncomp
        exit
4713    continue
      end do

      ! Set metadata.

      beam%ca(icut) = ca

      ! Read amplitudes.

      do i = 1, np
        read(1, *) tmp1, tmp2, tmp3, tmp4
        beam%amp(1, i, icut) = cmplx(tmp1, tmp2)
        beam%amp(2, i, icut) = cmplx(tmp3, tmp4)
      end do

    end do

    close(1)

  end subroutine bm_cut_read

  !======================================================================

  subroutine bm_cut_free(beam)
    type(bmcut), intent(inout) :: beam

    beam%ncomp = 0
    beam%np = 0
    beam%ncut = 0
    deallocate(beam%ca)
    deallocate(beam%amp)

  end subroutine bm_cut_free

  !======================================================================

end module beam_cut
