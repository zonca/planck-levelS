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

module ls_paramfile_io
  use planck_config
  use ls_misc_utils
  implicit none
  private

  public paramfile_handle, parse_init, parse_real, parse_double, parse_int, &
         parse_long, parse_lgt, parse_string, parse_finish, key_present

  type paramfile_handle
    character(len=filenamelen) filename
    character(len=filenamelen), pointer, dimension(:) :: keylist=>NULL()
    character(len=filenamelen), pointer, dimension(:) :: valuelist=>NULL()
    logical, pointer, dimension(:) :: usedlist=>NULL()
    logical interactive, verbose
  end type paramfile_handle

  interface parse_init
    module procedure parse_init_file, parse_init_arrays
  end interface

  character(len=*), parameter, public :: ret = achar(10)//' '

contains

subroutine notify_user (keyname, rdef, rmin, rmax, ddef, dmin, dmax, &
  idef, imin, imax, ldef, lmin, lmax, logdef, chdef, descr)
  character(len=*), intent(in) :: keyname
  real(sp), intent(in), optional :: rdef, rmin, rmax
  real(dp), intent(in), optional :: ddef, dmin, dmax
  integer(i4b), intent(in), optional :: idef, imin, imax
  integer(i8b), intent(in), optional :: ldef, lmin, lmax
  logical, intent(in), optional :: logdef
  character(len=*), intent(in), optional :: chdef, descr

  if (present(descr)) then
    print *, trim(descr)
  else
    print *, 'Please enter a value for the key ', keyname
  endif
  if (present(rmin)) print *, "min value: ", rmin
  if (present(dmin)) print *, "min value: ", dmin
  if (present(imin)) print *, "min value: ", imin
  if (present(lmin)) print *, "min value: ", lmin
  if (present(rmax)) print *, "max value: ", rmax
  if (present(dmax)) print *, "max value: ", dmax
  if (present(imax)) print *, "max value: ", imax
  if (present(lmax)) print *, "max value: ", lmax
  if (present(rdef)) print *, "default value: ", rdef
  if (present(ddef)) print *, "default value: ", ddef
  if (present(idef)) print *, "default value: ", idef
  if (present(ldef)) print *, "default value: ", ldef
  if (present(logdef)) print *, "default value: ", logdef
  if (present(chdef)) print *, "default value: ", trim(chdef)
end subroutine notify_user

function parse_init_file (fname, verbose)
  character(len=*), intent(in) :: fname
  logical, intent(in), optional :: verbose
  type(paramfile_handle) parse_init_file

  parse_init_file%verbose=.true.
  if (present(verbose)) parse_init_file%verbose=verbose
  if (len(trim(fname))==0) then
    parse_init_file%filename = ''
    parse_init_file%interactive=.true.
  else
    parse_init_file%filename = fname
    parse_init_file%interactive=.false.
    call scan_parfile(fname,parse_init_file%keylist,parse_init_file%valuelist)
    allocate(parse_init_file%usedlist(size(parse_init_file%keylist)))
    parse_init_file%usedlist=.false.
  endif
end function parse_init_file

function parse_init_arrays (key, val, verbose)
  character(len=*), intent(in) :: key(:), val(:)
  logical, intent(in), optional :: verbose
  type(paramfile_handle) parse_init_arrays

  integer cnt

  call assert(size(key)==size(val),"key and value arrays have different size")
  cnt = size(key)
  parse_init_arrays%filename = ''
  parse_init_arrays%verbose=.true.
  if (present(verbose)) parse_init_arrays%verbose=verbose
  parse_init_arrays%interactive=.false.
  allocate(parse_init_arrays%keylist(cnt),parse_init_arrays%valuelist(cnt), &
           parse_init_arrays%usedlist(cnt))
  parse_init_arrays%keylist = key
  parse_init_arrays%valuelist = val
  parse_init_arrays%usedlist = .false.
end function parse_init_arrays

subroutine parse_finish (handle, verbose)
  type(paramfile_handle), intent(inout) :: handle
  logical, intent(in), optional :: verbose
  integer i
  logical verb2

  verb2 = .true.
  if (present(verbose)) verb2 = verbose
  if (verb2) then
    do i=1,size(handle%keylist)
      if (.not. handle%usedlist(i)) &
        print *, "Parser warning: unused parameter '", &
                 trim(handle%keylist(i)),"'"
    enddo
  endif

  if (associated(handle%keylist)) &
    deallocate(handle%keylist, handle%valuelist, handle%usedlist)
end subroutine parse_finish

function key_present (handle, keyname)
  type(paramfile_handle), intent(in) :: handle
  character(len=*), intent(in) :: keyname
  logical key_present
  integer i

  if (handle%interactive) then
    key_present=.true.
  else
    key_present=.false.
    do i=1,size(handle%keylist)
      if (trim(handle%keylist(i))==keyname) then
        key_present=.true.
        return
      end if
    end do
  endif
end function key_present

subroutine find_param (handle,keyname,result,found,rdef,rmin,rmax, &
    ddef,dmin,dmax,idef,imin,imax,ldef,lmin,lmax,logdef,chdef,descr)
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  character(len=*), intent(out) :: result
  logical, intent(out) :: found
  real(sp), intent(in), optional :: rdef, rmin, rmax
  real(dp), intent(in), optional :: ddef, dmin, dmax
  integer(i4b), intent(in), optional :: idef, imin, imax
  integer(i8b), intent(in), optional :: ldef, lmin, lmax
  logical, intent(in), optional :: logdef
  character(len=*), intent(in), optional :: chdef, descr

  integer i

  found=.false.

  if (.not. handle%interactive) then
    do i=1,size(handle%keylist)
      if (trim(handle%keylist(i))==keyname) then
        result=trim(handle%valuelist(i))
        found=.true.
        handle%usedlist(i) = .true.
!        return
      end if
    end do
  else
    call notify_user (keyname,rdef,rmin,rmax,ddef,dmin,dmax, &
                      idef,imin,imax,ldef,lmin,lmax,logdef,chdef,descr)
    read (*,'(a)',err=5) result
    found = (trim(result)/='')
  end if

  return

5 print*,"Parser: could not parse '",trim(result), &
    "' as string value for key '",trim(keyname),"'."
  call exit_with_status(1)
end subroutine find_param

function parse_real (handle, keyname, default, vmin, vmax, descr)
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  real(sp), intent(in), optional :: default, vmin, vmax
  character(len=*), intent(in), optional :: descr
  real(sp) :: parse_real

  character(len=filenamelen) :: result
  logical found

10 continue

  call find_param (handle, trim(keyname), result, found, rdef=default, &
                   rmin=vmin, rmax=vmax, descr=descr)
  if (found) then
    read (result,*,err=5) parse_real
  else
    if (present(default)) then
      if (handle%verbose) &
        print *,"Parser: warning: using default value for '",trim(keyname),"'"
      parse_real = default
    else
      print *,"Parser: error: key '",trim(keyname),"' not found."
      goto 2
    endif
  endif
  if (handle%verbose) print *,"Parser: ",trim(keyname)," = ",parse_real
  if (present(vmin)) then
    if (parse_real<vmin) then
      print *,"Parser: error: value for '", trim(keyname),"' too small."
      goto 2
    endif
  endif
  if (present(vmax)) then
    if (parse_real>vmax) then
      print *,"Parser: error: value for '", trim(keyname),"' too large."
      goto 2
    endif
  endif

  return

5 print*,"Parser: could not parse '",trim(result), &
    "' as single precision value for key '",trim(keyname),"'."

2 if (handle%interactive) goto 10 ! try again

  call exit_with_status(1)
end function parse_real

function parse_double (handle, keyname, default, vmin, vmax, descr)
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  real(dp), intent(in), optional :: default, vmin, vmax
  character(len=*), intent(in), optional :: descr
  real(dp) :: parse_double

  character(len=filenamelen) :: result
  logical found

10 continue

  call find_param (handle, trim(keyname), result, found, ddef=default, &
                   dmin=vmin, dmax=vmax, descr=descr)
  if (found) then
    read (result,*,err=5) parse_double
  else
    if (present(default)) then
      if (handle%verbose) &
        print *,"Parser: warning: using default value for '",trim(keyname),"'"
      parse_double = default
    else
      print *,"Parser: error: key '",trim(keyname),"' not found."
      goto 2
    endif
  endif
  if (handle%verbose) print *,"Parser: ",trim(keyname)," = ",parse_double
  if (present(vmin)) then
    if (parse_double<vmin) then
      print *,"Parser: error: value for '", trim(keyname),"' too small."
      goto 2
    endif
  endif
  if (present(vmax)) then
    if (parse_double>vmax) then
      print *,"Parser: error: value for '", trim(keyname),"' too large."
      goto 2
    endif
  endif

  return

5 print*,"Parser: could not parse '",trim(result), &
    "' as double precision value for key '",trim(keyname),"'."

2 if (handle%interactive) goto 10 ! try again

  call exit_with_status(1)
end function parse_double

function parse_int (handle, keyname, default, vmin, vmax, descr)
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  integer(i4b), intent(in), optional :: default, vmin, vmax
  character(len=*), intent(in), optional :: descr
  integer(i4b) :: parse_int

  character(len=filenamelen) :: result
  logical found

10 continue

  call find_param (handle, trim(keyname), result, found, idef=default, &
                   imin=vmin, imax=vmax, descr=descr)
  if (found) then
    read (result,*,err=5) parse_int
  else
    if (present(default)) then
      if (handle%verbose) &
        print *,"Parser: warning: using default value for '",trim(keyname),"'"
      parse_int = default
    else
      print *,"Parser: error: key '",trim(keyname),"' not found."
      goto 2
    endif
  endif
  if (handle%verbose) print *,"Parser: ",trim(keyname)," = ",parse_int
  if (present(vmin)) then
    if (parse_int<vmin) then
      print *,"Parser: error: value for '", trim(keyname),"' too small."
      goto 2
    endif
  endif
  if (present(vmax)) then
    if (parse_int>vmax) then
      print *,"Parser: error: value for '", trim(keyname),"' too large."
      goto 2
    endif
  endif

  return

5 print*,"Parser: could not parse '",trim(result), &
    "' as 32bit integer value for key '",trim(keyname),"'."

2 if (handle%interactive) goto 10 ! try again

  call exit_with_status(1)
end function parse_int

function parse_long (handle, keyname, default, vmin, vmax, descr)
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  integer(i8b), intent(in), optional :: default, vmin, vmax
  character(len=*), intent(in), optional :: descr
  integer(i8b) :: parse_long

  character(len=filenamelen) :: result
  logical found

10 continue

  call find_param (handle, trim(keyname), result, found, ldef=default, &
                   lmin=vmin, lmax=vmax, descr=descr)
  if (found) then
    read (result,*,err=5) parse_long
  else
    if (present(default)) then
      if (handle%verbose) &
        print *,"Parser: warning: using default value for '",trim(keyname),"'"
      parse_long = default
    else
      print *,"Parser: error: key '",trim(keyname),"' not found."
      goto 2
    endif
  endif
  if (handle%verbose) print *,"Parser: ",trim(keyname)," = ",parse_long
  if (present(vmin)) then
    if (parse_long<vmin) then
      print *,"Parser: error: value for '", trim(keyname),"' too small."
      goto 2
    endif
  endif
  if (present(vmax)) then
    if (parse_long>vmax) then
      print *,"Parser: error: value for '", trim(keyname),"' too large."
      goto 2
    endif
  endif

  return

5 print*,"Parser: could not parse '",trim(result), &
    "' as 64bit integer value for key '",trim(keyname),"'."

2 if (handle%interactive) goto 10 ! try again

  call exit_with_status(1)
end function parse_long

function parse_lgt (handle, keyname, default, descr)
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  logical, intent(in), optional :: default
  character(len=*), intent(in), optional :: descr
  logical :: parse_lgt

  character(len=filenamelen) :: result
  logical found

10 continue

  call find_param (handle, trim(keyname), result, found, logdef=default, &
                   descr=descr)
  if (found) then
    select case (result)
      case ('y','Y','t','T','true','TRUE','.true.','.TRUE.')
        parse_lgt = .true.
      case ('n','N','f','F','false','FALSE','.false.','.FALSE.')
        parse_lgt= .false.
      case default
        goto 5
    end select
  else
    if (present(default)) then
      if (handle%verbose) &
        print *,"Parser: warning: using default value for '",trim(keyname),"'"
      parse_lgt = default
    else
      print *,"Parser: error: key '",trim(keyname),"' not found."
      goto 2
    endif
  endif
  if (handle%verbose) print *,"Parser: ",trim(keyname)," = ",parse_lgt

  return

5 print*,"Parser: could not parse '",trim(result), &
    "' as boolean value for key '",trim(keyname),"'."

2 if (handle%interactive) goto 10 ! try again

  call exit_with_status(1)
end function parse_lgt

function parse_string (handle, keyname, default, descr, filestatus)
  type(paramfile_handle), intent(inout) :: handle
  character(len=*), intent(in) :: keyname
  character(len=*), intent(in), optional :: default
  character(len=*), intent(in), optional :: descr
  character(len=*), intent(in), optional :: filestatus
  character(len=filenamelen) :: parse_string

  character(len=filenamelen) :: result
  logical found, there

  call find_param (handle, trim(keyname), result, found, chdef=default, &
                   descr=descr)
  if (found) then
    parse_string = trim(result)
  else
    if (present(default)) then
      if (handle%verbose) &
        print *,"Parser: warning: using default value for '",trim(keyname),"'"
      parse_string = trim(default)
    else
      print *,"Parser: error: key '",trim(keyname),"' not found."
      call exit_with_status(1)
    endif
  endif
  if (handle%verbose) &
    print *,"Parser: ",trim(keyname)," = '",trim(parse_string),"'"
  if (trim(adjustl(parse_string)) == "''")  parse_string = ''
  if (present(filestatus) .and. trim(parse_string) /= '') then
    if (trim(filestatus)=='new') then
      inquire(file=trim(parse_string),exist=there)
      if (there) then
        print *, "Parser: error: output file '" // trim(parse_string) // &
                 "' already exists!"
        call exit_with_status(1)
      end if
    else if (trim(filestatus)=='old') then
      inquire(file=trim(parse_string),exist=there)
      if (.not. there) then
        print *, "Parser: error: input file '" // trim(parse_string) // &
                 "' does not exist!"
        call exit_with_status(1)
      end if
    else
      print *, "Parser: error: wrong value for filestatus"
      call exit_with_status(1)
    endif
  endif
end function parse_string

end module ls_paramfile_io
