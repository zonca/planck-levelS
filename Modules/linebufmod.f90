!-----------------------------------------------------------------------------
!
!  This file is part of the Planck simulation package
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version.

!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.

!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!-----------------------------------------------------------------------------

!  The Planck simulation package is being developed at the Max-Planck-Institut
!  fuer Astrophysik and financially supported by the Deutsches Zentrum fuer
!  Luft- und Raumfahrt (DLR).

!  Copyright (C) 2002-2013 Max-Planck-Society
!  \author Martin Reinecke

module linebufmod
implicit none
private

public linebuf, linebuf_add, linebuf_clear, linebuf_set

type linebuf
  character(len=80), dimension(:), pointer :: data=>NULL()
  integer :: size=0
end type linebuf

contains

subroutine linebuf_add (buf, line)
  type(linebuf), intent(inout) :: buf
  character(len=*), intent(in) :: line
  character(len=80), dimension(:), pointer :: data2
  data2=>NULL()

  if (.not. associated(buf%data)) allocate(buf%data(10))

  if (size(buf%data)==buf%size) then
    allocate (data2(2*size(buf%data)))
    data2(1:size(buf%data)) = buf%data
    deallocate (buf%data)
    buf%data => data2
  endif

  buf%size=buf%size+1
  buf%data(buf%size) = trim (line)
end subroutine linebuf_add

subroutine linebuf_clear (buf)
  type(linebuf), intent(inout) :: buf

  if (associated(buf%data)) deallocate(buf%data)
  buf%size=0
end subroutine linebuf_clear

subroutine linebuf_set (buf, orig)
  type(linebuf), intent(inout) :: buf
  character(len=80), dimension(:), intent(in) :: orig

  if (associated(buf%data)) deallocate(buf%data)

  allocate(buf%data(size(orig)))
  buf%data=orig
  buf%size=size(orig)
end subroutine linebuf_set

end module linebufmod
