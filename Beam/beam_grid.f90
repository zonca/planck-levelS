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

module beam_grid

  use planck_config
  use ls_misc_utils

  implicit none

  ! Beam type corresponding to Grasp "grid" file format.

  type bmgrid
    integer :: ncomp, nx, ny
    integer :: kcomp, kgrid
    integer :: ix, iy
    real(dp) :: xs, ys, xe, ye
    complex(dp), dimension(:,:,:), pointer :: amp
  end type bmgrid

contains

  !======================================================================

  subroutine bm_grid_read(beam, filename)
    type(bmgrid), intent(inout) :: beam
    character(len=*), intent(in) :: filename

    integer :: i, j
    integer :: ktype, nset, klimit
    integer :: is, in, junk1, junk2
    real(dp) :: tmp1, tmp2, tmp3, tmp4
    character(len=80) :: line

    ! Set beam type from input.

    open(1, file=filename, status='old')

    ! Read lines until one starting '++++' is found.

    do
      read(1, '(a)') line
      if (line(1:4) == '++++') exit
    end do

    read(1, *) ktype

    call assert (ktype==1, 'Unknown Grasp grid format, ktype /= 1')

    read(1, *) nset, beam%kcomp, beam%ncomp, beam%kgrid

    if (nset > 1) &
      print *, 'Warning: nset > 1, only reading first beam in file.'

    call assert (beam%ncomp<=2, 'Three field components present.&
          & Beam package can only handle two field components.')

    read(1, *) beam%ix, beam%iy

    do i = 2, nset
      read(1, *) junk1, junk2
    end do

    ! Read grid limits.

    read(1, *) beam%xs, beam%ys, beam%xe, beam%ye
    read(1, *) beam%nx, beam%ny, klimit

    allocate(beam%amp(beam%ncomp, beam%nx, beam%ny))

    beam%amp = (0.0_dp, 0.0_dp)

    do j = 1, beam%ny

      if (klimit == 1) then
        read(1, *) is, in
      else
        is = 1
        in = beam%nx
      end if

      do i = is, is+in-1
        read(1, *) tmp1, tmp2, tmp3, tmp4
        beam%amp(1, i, j) = cmplx(tmp1, tmp2)
        beam%amp(2, i, j) = cmplx(tmp3, tmp4)
      end do

    end do

    close(1)

  end subroutine bm_grid_read

  !======================================================================

  subroutine bm_grid_free(beam)
    type(bmgrid), intent(inout) :: beam

    beam%ncomp = 0
    beam%nx = 0
    beam%ny = 0
    deallocate(beam%amp)

  end subroutine bm_grid_free

  !======================================================================

end module beam_grid
