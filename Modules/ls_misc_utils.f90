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

module ls_misc_utils
use planck_config
implicit none
private

public fatal_error, assert, assert_present, assert_not_present, &
       assert_alloc, file_present, toString, getcmdline, scan_args, &
       scan_parfile

interface toString
  module procedure toString_i4b, toString_i8b, toString_sp, toString_dp, &
                   toString_l
end interface

contains

subroutine fatal_error (msg)
  character(len=*), intent(in), optional :: msg

  if (present(msg)) then
    print *,'Fatal error: ', trim(msg)
  else
    print *,'Fatal error'
  endif
  call exit_with_status(1)
end subroutine fatal_error

function file_present (filename)
  character(len=*), intent(in) :: filename
  logical :: file_present

  inquire(file=trim(filename),exist=file_present)
end function

function getcmdline()
  character(len=filenamelen),dimension(:),allocatable :: getcmdline

  integer m

  allocate(getcmdline(command_argument_count()))
  do m=1,command_argument_count()
    call get_command_argument(m, getcmdline(m))
  end do
end function

subroutine scan_args(args,keys,values)
  character(len=*),intent(in) :: args(:)
  character(len=*),pointer, intent(out) :: keys(:), values(:)

  integer m,eqpos

  allocate(keys(size(args)),values(size(args)))
  do m=1,size(args)
    eqpos=scan(args(m),'=')
    call assert(eqpos/=0,"bad command line")
    keys(m)=trim(adjustl(args(m)(:eqpos-1)))
    values(m)=trim(adjustl(args(m)(eqpos+1:)))
  end do
end subroutine

subroutine scan_parfile(name,keys,values)
  character(len=*),intent(in) :: name
  character(len=filenamelen),pointer, intent(out) :: keys(:), values(:)

  integer i,cnt
  character(len=filenamelen) :: line

  call assert_present(name)
  if (len(name)>filenamelen) then
    print *,"Parser: error: file name too long"
    call exit_with_status(1)
  endif
  open (1, file=trim(name))
  cnt=0
  do
    read (1,'(a)',end=2) line
    i=scan(line,'#')
    if (i/=0) line=line(1:i-1)
    i=scan(line,'=')
    if (i/=0) cnt=cnt+1
  end do
2 close (1)
  allocate(keys(cnt), values(cnt))
  open (1, file=trim(name))
  cnt=0
  do
    read (1,'(a)',end=3) line
    i=scan(line,'#')
    if (i/=0) line=line(1:i-1)
    i=scan(line,'=')
    if (i/=0) then
      cnt=cnt+1
      keys(cnt)=trim(adjustl(line(:i-1)))
      values(cnt)=trim(adjustl(line(i+1:)))
    endif
  end do
3 close (1)
end subroutine

subroutine assert_present (filename)
  character(len=*), intent(in) :: filename

  if (.not. file_present(trim(filename))) then
    print *, "Error:  file '" // trim(filename) // "' does not exist!"
    call exit_with_status(1)
  end if
end subroutine assert_present

subroutine assert_not_present (filename)
  character(len=*), intent(in) :: filename

  if (file_present(trim(filename))) then
    print *, "Error:  file '" // trim(filename) // "' already exists!"
    call exit_with_status(1)
  end if
end subroutine assert_not_present

subroutine assert_alloc (stat,code,arr)
  integer, intent(in) :: stat
  character(len=*), intent(in) :: code, arr

  if (stat==0) return

  print *, trim(code)//'> cannot allocate memory for array: '//trim(arr)
  call exit_with_status(1)
end subroutine assert_alloc

subroutine assert (testval,msg,errcode)
  logical, intent(in) :: testval
  character(len=*), intent(in), optional :: msg
  integer, intent(in), optional :: errcode

  if (testval) return

  print *,"Assertion failed: "
  if (present(msg)) print *, trim(msg)
  if (present(errcode)) call exit_with_status (errcode)
  call exit_with_status(1)
end subroutine

function toString_i4b (num)
  integer(i4b), intent(in) :: num
  character (len=30) :: toString_i4b

  write(toString_i4b,*) num
  toString_i4b = trim(adjustl(toString_i4b))
end function

function toString_i8b (num)
  integer(i8b), intent(in) :: num
  character (len=30) :: toString_i8b

  write(toString_i8b,*) num
  toString_i8b = trim(adjustl(toString_i8b))
end function

function toString_sp (num)
  real(sp), intent(in) :: num
  character (len=30) :: toString_sp

  write(toString_sp,*) num
  toString_sp = trim(adjustl(toString_sp))
end function

function toString_dp (num)
  real(dp), intent(in) :: num
  character (len=30) :: toString_dp

  write(toString_dp,*) num
  toString_dp = trim(adjustl(toString_dp))
end function

function toString_l (val)
  logical, intent(in) :: val
  character (len=10) :: toString_l

  write(toString_l,*) val
  toString_l = trim(adjustl(toString_l))
end function

end module ls_misc_utils
