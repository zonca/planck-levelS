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

!  Copyright (C) 2005-2013 Max-Planck-Society
!  \author Martin Reinecke

! make sure at most one DMC is specified
#if ((defined(USE_TOODI)) && (defined(USE_HFIDMC)))
#error support for both DMCs requested
#elif ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
! DMC mode

!#define CHECK(X) print *,"calling ", 'X'; status = X; if (status/=0) print *,"status = ",status; call assert(status==0,"DMC error"); print *,"... done"
#define CHECK(X) status = X; call assert(status==0,"DMC error",status)

!! DMC I/O routines
module cll_io

use planck_config
use planck_types
use ls_misc_utils
use ls_paramfile_io
#ifdef USE_TOODI
use cllio
#else
use llio
#endif
implicit none
private

#ifdef USE_TOODI
integer, parameter :: CLL_HANDLE_KIND=i4b
integer, parameter :: CLL_STRLEN=filenamelen
#else
integer, parameter :: CLL_HANDLE_KIND=i8b
integer, parameter :: CLL_STRLEN=DMCSTRINGMAXLEN
#endif

integer(CLL_HANDLE_KIND), save :: inithandle, session
integer(i4b), save :: status

character(len=filenamelen),dimension(:),allocatable,save :: args

!! data type containing information about a single column in a DMC object
type cll_column
  !! column name
  character(len=CLL_STRLEN) :: name=''
  !! physical unit of the column
  character(len=CLL_STRLEN) :: unit=''
  !! data type of the column
  integer(i4b) :: type=PLANCK_INVALID_TYPE
end type

public cll_handle, cll_column, cll_open, cll_close, cll_flush, &
  cll_create, cll_set_key, cll_get_key, cll_delete_key, &
  cll_key_present, cll_append_column, cll_read_column, &
  cll_colnum, cll_collength, cll_init, cll_shutdown, cll_getrunparams, &
  cll_numcols, cll_get_colname, cll_column_present, cll_get_colunit, &
  cll_coltype, cll_eff_chunksize

interface cll_open
  module procedure cll_open_generic, cll_open_with_type
end interface

interface cll_set_key
  module procedure cll_set_key_int1, cll_set_key_int2, cll_set_key_int4, &
    cll_set_key_int8, cll_set_key_real4, cll_set_key_real8, cll_set_key_bool, &
    cll_set_key_string
end interface

interface cll_get_key
  module procedure cll_get_key_int1, cll_get_key_int2, cll_get_key_int4, &
    cll_get_key_int8, cll_get_key_real4, cll_get_key_real8, cll_get_key_bool, &
    cll_get_key_string
end interface

interface cll_append_column
  module procedure cll_append_column_real4, cll_append_column_real8, &
    cll_append_column_int1, cll_append_column_int2, cll_append_column_int4, &
    cll_append_column_int8, cll_append_column_string
end interface

interface cll_read_column
  module procedure cll_read_column_real4, cll_read_column_real8, &
    cll_read_column_int1, cll_read_column_int2, cll_read_column_int4, &
    cll_read_column_int8, cll_read_column_string
end interface

!! data type containing information about a DMC object
type cll_handle
  integer(CLL_HANDLE_KIND) :: hnd = -1
  logical :: connected = .false.
  logical :: readonly = .true.
  character(len=CLL_STRLEN) :: type = ''

  type(cll_column), dimension(:), pointer :: columns=>NULL()
end type

contains

subroutine assert_connected(handle)
  type(cll_handle), intent(inout) :: handle

  call assert(handle%connected,"trying to access an unconnected DMC handle")
end subroutine

subroutine assert_read(handle)
  type(cll_handle), intent(inout) :: handle

  call assert(handle%readonly,"trying to read from a write-only DMC handle")
end subroutine

subroutine assert_write(handle)
  type(cll_handle), intent(inout) :: handle

  call assert(.not. handle%readonly,"trying to write to a read-only DMC handle")
end subroutine

function cll_getrunparams (verbose)
  logical, intent(in), optional :: verbose
  type(paramfile_handle) :: cll_getrunparams

  character(len=CLL_STRLEN),dimension(:),pointer :: keys, values
  integer nk, m

  if ((size(args)==1) .and. (scan(args(1),'=')==0)) then
    nullify(keys)
    CHECK(dmcgetmetakeys(inithandle,nk,keys))
    allocate(values(nk))
    do m=1,nk
      CHECK(dmcgetmetavalueasstring(inithandle,keys(m),values(m)))
    end do
  else
    call scan_args(args,keys,values)
  endif
  cll_getrunparams = parse_init(keys,values,verbose)
  deallocate(keys,values)
end function

subroutine cll_init (xargs)
  character(len=*), intent(in) :: xargs(:)

  logical, save :: initialised = .false.

  if (.not. initialised) then
    call assert(size(xargs)>0,"bad argument format")
    allocate(args(size(xargs)))
    args=xargs
    if ((size(args)==1) .and. (scan(args(1),'=')==0)) then
#ifdef USE_TOODI
! FIXME: this is not conforming to the CLLIO specification
      if (args(1)=='') then
        CHECK(dmcinitializellio('TOODI%file', inithandle))
      else
        CHECK(dmcinitializellio(trim(args(1)), inithandle))
      endif
#else
      CHECK(dmcinitializellio(trim(args(1)), inithandle))
#endif
      CHECK(dmcopenstores(inithandle,session))
      CHECK(dmcbegintransaction(session))
    else
#ifdef USE_TOODI
! FIXME: this is not conforming to the CLLIO specification
      CHECK(dmcinitializellio('TOODI%file', inithandle))
#else
      CHECK(dmcinitializellio('', inithandle))
#endif
      CHECK(dmcopenstores(inithandle,session))
      CHECK(dmcbegintransaction(session))
    endif
    initialised = .true.
  else
    print *,"Warning: tried to initialise the DMC more than once"
  endif
end subroutine

subroutine cll_shutdown ()
  logical, save :: done = .false.

  if (.not. done) then
    CHECK(dmccommittransaction(session))
    CHECK(dmcclosestores(session))
    CHECK(dmcclosellio(inithandle))
    deallocate(args)
    done = .true.
  else
    print *,"Warning: tried to shut down the DMC more than once"
  endif
end subroutine

subroutine clean_all (handle)
  type(cll_handle), intent(inout) :: handle

  if (.not. handle%connected) return
  CHECK(dmccloseobject(handle%hnd))
  if (.not. handle%readonly) CHECK(dmccommitandresume(session))
  handle%hnd=-1
  handle%connected=.false.
  handle%readonly=.true.
  handle%type=''
  deallocate(handle%columns)
  nullify(handle%columns)
end subroutine

subroutine cll_fill_cache (handle)
  type(cll_handle), intent(inout) :: handle

  integer ncols, m
  character(len=CLL_STRLEN) :: ptype

  call assert_connected(handle)
#ifdef USE_TOODI
  CHECK(dmcgetmetavalueasstring(handle%hnd,"objType",handle%type))
#else
  CHECK(dmcgetobjecttype(handle%hnd,handle%type))
#endif
  CHECK(dmcgetnumberofcolumns(handle%hnd,ncols))
  allocate(handle%columns(0:ncols-1))
  do m=0,ncols-1
    CHECK(dmcgetnameofcolumn(handle%hnd,m,handle%columns(m)%name))
    CHECK(dmcgetunitofcolumn(handle%hnd,m,handle%columns(m)%unit))
    CHECK(dmcgettypeofcolumn(handle%hnd,m,ptype))
    if (ptype=="float64") then
      handle%columns(m)%type=PLANCK_FLOAT64
    else if (ptype=="float32") then
      handle%columns(m)%type=PLANCK_FLOAT32
    else if (ptype=="byte") then
      handle%columns(m)%type=PLANCK_INT8
    else if (ptype=="char") then
      handle%columns(m)%type=PLANCK_INT8
    else if (ptype=="int16") then
      handle%columns(m)%type=PLANCK_INT16
    else if (ptype=="int32") then
      handle%columns(m)%type=PLANCK_INT32
    else if (ptype=="int64") then
      handle%columns(m)%type=PLANCK_INT64
    else if (ptype=="string") then
      handle%columns(m)%type=PLANCK_STRING
    else if (ptype=="boolean") then
      handle%columns(m)%type=PLANCK_BOOL
    else
      print *,"unknown data type '", ptype, "'"
      call exit_with_status(1)
    endif
  end do

end subroutine

subroutine cll_close (handle)
  type(cll_handle), intent(inout) :: handle

  call assert_connected(handle)
  CHECK(dmccommitandresume(session))
  call clean_all(handle)
end subroutine

subroutine cll_flush (handle)
  type(cll_handle), intent(inout) :: handle

  call assert_connected(handle)
  CHECK(dmccommitandresume(session))
end subroutine

subroutine cll_open_generic (handle, name)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name

  call clean_all(handle)
  CHECK(dmcopenpersistentobject (session,"common.GenericCore",name,handle%hnd))
  handle%connected = .true.
  handle%readonly = .true.
  call cll_fill_cache (handle)
end subroutine

subroutine cll_open_with_type (handle, name, type)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name, type

  call clean_all(handle)
  CHECK(dmcopenpersistentobject (session,type,name,handle%hnd))
  handle%connected = .true.
  handle%readonly = .true.
  call cll_fill_cache (handle)
end subroutine

subroutine cll_create (handle, name, type)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name, type

  call clean_all(handle)
  CHECK(dmcopennewobject(session,type,name,handle%hnd))
  handle%connected = .true.
  handle%readonly = .false.
  call cll_fill_cache(handle)
end subroutine

subroutine cll_set_key_int1 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i1b), intent(in) :: value

  call assert_write(handle)
  CHECK(dmcsetmetavalueasbyte(handle%hnd,key,value))
end subroutine

subroutine cll_set_key_int2 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i2b), intent(in) :: value

  call assert_write(handle)
  CHECK(dmcsetmetavalueasint16(handle%hnd,key,value))
end subroutine

subroutine cll_set_key_int4 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i4b), intent(in) :: value

  call assert_write(handle)
  CHECK(dmcsetmetavalueasint32(handle%hnd,key,value))
end subroutine

subroutine cll_set_key_int8 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i8b), intent(in) :: value

  call assert_write(handle)
  CHECK(dmcsetmetavalueasint64(handle%hnd,key,value))
end subroutine

subroutine cll_set_key_real4 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(sp), intent(in) :: value

  call assert_write(handle)
  CHECK(dmcsetmetavalueasfloat(handle%hnd,key,value))
end subroutine

subroutine cll_set_key_real8 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(dp), intent(in) :: value

  call assert_write(handle)
  CHECK(dmcsetmetavalueasdouble(handle%hnd,key,value))
end subroutine

subroutine cll_set_key_bool (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical, intent(in) :: value

  logical(DMCBOOL) tmp

  call assert_write(handle)
  tmp = value
  CHECK(dmcsetmetavalueasboolean(handle%hnd,key,tmp))
end subroutine

subroutine cll_set_key_string (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in) :: value

  call assert_write(handle)
  CHECK(dmcsetmetavalueasstring(handle%hnd,key,value))
end subroutine

subroutine cll_delete_key (handle,key)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key

  call assert_write(handle)
  CHECK(dmcremovemetakey(handle%hnd,key))
end subroutine

function cll_key_present (handle, key)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical cll_key_present

  logical(DMCBOOL) tmp

  call assert_read(handle)
  CHECK(dmcmetakeyispresent(handle%hnd,key,tmp))
  cll_key_present = tmp
end function

subroutine cll_get_key_int1 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i1b), intent(out) :: value

  call assert_read(handle)
  CHECK(dmcgetmetavalueasbyte(handle%hnd,key,value))
end subroutine

subroutine cll_get_key_int2 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i2b), intent(out) :: value

  call assert_read(handle)
  CHECK(dmcgetmetavalueasint16(handle%hnd,key,value))
end subroutine

subroutine cll_get_key_int4 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i4b), intent(out) :: value

  call assert_read(handle)
  CHECK(dmcgetmetavalueasint32(handle%hnd,key,value))
end subroutine

subroutine cll_get_key_int8 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i8b), intent(out) :: value

  call assert_read(handle)
  CHECK(dmcgetmetavalueasint64(handle%hnd,key,value))
end subroutine

subroutine cll_get_key_real4 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(sp), intent(out) :: value

  call assert_read(handle)
  CHECK(dmcgetmetavalueasfloat(handle%hnd,key,value))
end subroutine

subroutine cll_get_key_real8 (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(dp), intent(out) :: value

  call assert_read(handle)
  CHECK(dmcgetmetavalueasdouble(handle%hnd,key,value))
end subroutine

subroutine cll_get_key_bool (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical, intent(out) :: value

  logical(DMCBOOL) tmp

  call assert_read(handle)
  CHECK(dmcgetmetavalueasboolean(handle%hnd,key,tmp))
  value = tmp
end subroutine

subroutine cll_get_key_string (handle,key,value)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(out) :: value

  call assert_read(handle)
  CHECK(dmcgetmetavalueasstring(handle%hnd,key,value))
end subroutine

subroutine cll_append_column_real4 (handle, colnum, data)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(sp), intent(in) :: data(:)

  real(dp), allocatable :: datad(:)
  integer(i8b) :: nelems

  call assert_write(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_FLOAT32)
      CHECK(dmcappendfloatvalues(handle%hnd,colnum,nelems,data,-1))
    case (PLANCK_FLOAT64)
      allocate(datad(size(data)))
      datad=data
      CHECK(dmcappenddoublevalues(handle%hnd,colnum,nelems,datad,-1))
      deallocate(datad)
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_append_column_real8 (handle, colnum, data)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(dp), intent(in) :: data(:)

  real(sp), allocatable :: datas(:)
  integer(i8b) :: nelems

  call assert_write(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_FLOAT64)
      CHECK(dmcappenddoublevalues(handle%hnd,colnum,nelems,data,-1))
    case (PLANCK_FLOAT32)
      allocate(datas(size(data)))
      datas=data
      CHECK(dmcappendfloatvalues(handle%hnd,colnum,nelems,datas,-1))
      deallocate(datas)
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_append_column_int1 (handle, colnum, data)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i1b), intent(in) :: data(:)

  integer(i8b) :: nelems

  call assert_write(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_INT8)
      CHECK(dmcappendbytevalues(handle%hnd,colnum,nelems,data,-1))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_append_column_int2 (handle, colnum, data)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i2b), intent(in) :: data(:)

  integer(i8b) :: nelems

  call assert_write(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_INT16)
      CHECK(dmcappendint16values(handle%hnd,colnum,nelems,data,-1))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_append_column_int4 (handle, colnum, data)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b), intent(in) :: data(:)

  real(sp), allocatable :: datas(:)
  real(dp), allocatable :: datad(:)
  integer(i8b) :: nelems

  call assert_write(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_INT32)
      CHECK(dmcappendint32values(handle%hnd,colnum,nelems,data,-1))
    case (PLANCK_FLOAT32)
      allocate(datas(size(data)))
      datas=data
      CHECK(dmcappendfloatvalues(handle%hnd,colnum,nelems,datas,-1))
      deallocate(datas)
    case (PLANCK_FLOAT64)
      allocate(datad(size(data)))
      datad=data
      CHECK(dmcappenddoublevalues(handle%hnd,colnum,nelems,datad,-1))
      deallocate(datad)
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_append_column_int8 (handle, colnum, data)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(in) :: data(:)

  integer(i8b) :: nelems

  call assert_write(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_INT64)
      CHECK(dmcappendint64values(handle%hnd,colnum,nelems,data,-1))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_append_column_string (handle, colnum, data)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(in) :: data(:)

  integer(i8b) :: nelems

  call assert_write(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_STRING)
      CHECK(dmcappendstringvalues(handle%hnd,colnum,nelems,data,-1))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_read_column_real4 (handle, colnum, data, offset)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(sp), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  integer(i4b), allocatable :: datai(:)
  real(dp), allocatable :: datad(:)
  integer(i8b) :: nelems

  call assert_read(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_FLOAT32)
      CHECK (dmcgetfloatvalues(handle%hnd,colnum,offset,nelems,data))
    case (PLANCK_INT32)
      allocate(datai(size(data)))
      CHECK (dmcgetint32values(handle%hnd,colnum,offset,nelems,datai))
      data = datai
      deallocate(datai)
    case (PLANCK_FLOAT64)
      allocate(datad(size(data)))
      CHECK (dmcgetdoublevalues(handle%hnd,colnum,offset,nelems,datad))
      data = datad
      deallocate(datad)
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_read_column_real8 (handle, colnum, data, offset)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(dp), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  integer(i4b), allocatable :: datai(:)
  real(sp), allocatable :: datas(:)
  integer(i8b) :: nelems

  call assert_read(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_FLOAT64)
      CHECK (dmcgetdoublevalues(handle%hnd,colnum,offset,nelems,data))
    case (PLANCK_INT32)
      allocate(datai(size(data)))
      CHECK (dmcgetint32values(handle%hnd,colnum,offset,nelems,datai))
      data = datai
      deallocate(datai)
    case (PLANCK_FLOAT32)
      allocate(datas(size(data)))
      CHECK (dmcgetfloatvalues(handle%hnd,colnum,offset,nelems,datas))
      data = datas
      deallocate(datas)
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_read_column_int1 (handle, colnum, data, offset)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i1b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  integer(i8b) :: nelems

  call assert_read(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_INT8)
      CHECK (dmcgetbytevalues(handle%hnd,colnum,offset,nelems,data))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_read_column_int2 (handle, colnum, data, offset)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i2b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  integer(i8b) :: nelems

  call assert_read(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_INT16)
      CHECK (dmcgetint16values(handle%hnd,colnum,offset,nelems,data))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_read_column_int4 (handle, colnum, data, offset)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  integer(i8b) :: nelems

  call assert_read(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_INT32)
      CHECK (dmcgetint32values(handle%hnd,colnum,offset,nelems,data))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_read_column_int8 (handle, colnum, data, offset)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  integer(i8b) :: nelems

  call assert_read(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_INT64)
      CHECK (dmcgetint64values(handle%hnd,colnum,offset,nelems,data))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

subroutine cll_read_column_string (handle, colnum, data, offset)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  integer(i8b) :: nelems

  call assert_read(handle)
  nelems = size(data)
  if (nelems==0) return
  select case (handle%columns(colnum)%type)
    case (PLANCK_STRING)
      CHECK (dmcgetstringvalues(handle%hnd,colnum,offset,nelems,data))
    case default
      call exit_with_status(1,'incompatible data types')
  end select
end subroutine

function cll_column_present (handle, colname)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  logical cll_column_present

  integer m

  call assert_connected(handle)
  do m=0,size(handle%columns)-1
    if (handle%columns(m)%name==colname) then
      cll_column_present=.true.
      return
    endif
  end do

  cll_column_present=.false.
end function

function cll_colnum (handle, colname)
  type(cll_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer cll_colnum

  integer m

  call assert_connected(handle)
  do m=0,size(handle%columns)-1
    if (handle%columns(m)%name==colname) then
      cll_colnum=m
      return
    endif
  end do

  cll_colnum=-1
  call fatal_error("column name not found")
end function

function cll_collength (handle, colnum)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b) cll_collength

  call assert_read(handle)
  CHECK(dmcgetsizeofcolumn(handle%hnd,colnum,cll_collength))
end function

function cll_numcols (handle)
  type(cll_handle), intent(inout) :: handle
  integer(i4b) cll_numcols

  call assert_connected(handle)
  cll_numcols = size(handle%columns)
end function

subroutine cll_get_colname (handle, colnum, colname)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: colname

  call assert_connected(handle)
  colname = handle%columns(colnum)%name
end subroutine

subroutine cll_get_colunit (handle, colnum, colunit)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: colunit

  call assert_connected(handle)
  colunit = handle%columns(colnum)%unit
end subroutine

function cll_coltype (handle, colnum)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b) :: cll_coltype

  call assert_connected(handle)
  cll_coltype = handle%columns(colnum)%type
end function

function cll_eff_chunksize (handle, colnum)
  type(cll_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b) :: cll_eff_chunksize

  call assert_connected(handle)
#ifdef USE_HFIDMC
  cll_eff_chunksize = 1000000000000_i8b ! avoid chunks at all costs for HFI
#else
  cll_eff_chunksize = 1024*1024 ! no better heuristics available at the moment
#endif
end function

end module cll_io

#endif

!! CLLIO FITS backend
module fts_io

use planck_config
use planck_types
use ls_misc_utils
use ls_paramfile_io
use fitsmod2
implicit none
private

character(len=filenamelen), allocatable, dimension(:), save :: args

!! data type containing information about a single column in a FITS object
type fts_column
  !! column name
  character(len=80) :: name=''
  !! physical unit of the column
  character(len=80) :: unit=''
  !! data type of the column
  integer(i4b) :: type=PLANCK_INVALID_TYPE
  integer(i4b) :: hdu, idx, width
  integer(i8b) :: size
end type

public fts_handle, fts_column, fts_open, fts_close, fts_flush, &
  fts_create, fts_set_key, fts_get_key, fts_delete_key, &
  fts_key_present, fts_append_column, fts_read_column, fts_reserve_column, &
  fts_colnum, fts_collength, fts_init, fts_shutdown, fts_getrunparams, &
  fts_numcols, fts_get_colname, fts_column_present, fts_get_colunit, &
  fts_coltype, fts_eff_chunksize

interface fts_open
  module procedure fts_open_generic, fts_open_with_type
end interface

interface fts_set_key
  module procedure fts_set_key_int1, fts_set_key_int2, fts_set_key_int4, &
    fts_set_key_int8, fts_set_key_real4, fts_set_key_real8, fts_set_key_bool, &
    fts_set_key_string
end interface

interface fts_get_key
  module procedure fts_get_key_int1, fts_get_key_int2, fts_get_key_int4, &
    fts_get_key_int8, fts_get_key_real4, fts_get_key_real8, fts_get_key_bool, &
    fts_get_key_string
end interface

interface fts_append_column
  module procedure fts_append_column_real4, fts_append_column_real8, &
    fts_append_column_int1, fts_append_column_int2, fts_append_column_int4, &
    fts_append_column_int8, fts_append_column_string
end interface

interface fts_read_column
  module procedure fts_read_column_real4, fts_read_column_real8, &
    fts_read_column_int1, fts_read_column_int2, fts_read_column_int4, &
    fts_read_column_int8, fts_read_column_string
end interface

!! data type containing information about a FITS object
type fts_handle
  type(fitshandle) :: hnd
  logical :: connected = .false.
  logical :: readonly = .true.
  character(len=80) :: type = 'INVALID'

  type(fts_column), dimension(:), pointer :: columns=>NULL()
end type

contains

subroutine assert_connected(handle)
  type(fts_handle), intent(inout) :: handle

  call assert(handle%connected,"trying to access an unconnected FITS handle")
end subroutine

subroutine assert_read(handle)
  type(fts_handle), intent(inout) :: handle

  call assert(handle%readonly,"trying to read from a write-only FITS handle")
end subroutine

subroutine assert_write(handle)
  type(fts_handle), intent(inout) :: handle

  call assert(.not.handle%readonly,"trying to write to a read-only FITS handle")
end subroutine

function fts_getrunparams (verbose)
  logical, intent(in), optional :: verbose
  type(paramfile_handle) :: fts_getrunparams

  character(len=filenamelen), pointer :: keys(:)=>NULL(),values(:)=>NULL()

  if ((size(args)==1) .and. (scan(args(1),'=')==0)) then
    fts_getrunparams = parse_init(args(1),verbose)
  else
    call scan_args(args,keys,values)
    fts_getrunparams = parse_init(keys,values,verbose)
    deallocate(keys,values)
  endif
end function

subroutine fts_init (xargs)
  character(len=*), intent(in) :: xargs(:)

  call assert(size(xargs)>0,"incorrect command line format");
  allocate (args(size(xargs)))
  args=xargs
end subroutine

subroutine fts_shutdown ()
  deallocate(args)
end subroutine

subroutine clean_all (handle)
  type(fts_handle), intent(inout) :: handle

  if (.not. handle%connected) return
  call fits_close(handle%hnd)
  handle%connected=.false.
  handle%readonly=.true.
  handle%type='INVALID'
  deallocate(handle%columns)
  nullify(handle%columns)
end subroutine

subroutine fts_fill_cache (handle)
  type(fts_handle), intent(inout) :: handle

  integer ncols, m, n

  call fits_goto_hdu(handle%hnd,2)
  if (fits_key_present(handle%hnd,"objType")) then
    call fits_get_key(handle%hnd,"objType",handle%type)
  else
    handle%type = 'UNKNOWN'
  endif

  ncols = 0
  do m=2, fits_num_hdus(handle%hnd)
    call fits_goto_hdu(handle%hnd,m)
    ncols=ncols+size(handle%hnd%columns)
  end do
  allocate (handle%columns(0:ncols-1))

  ncols=0
  do m=2, fits_num_hdus(handle%hnd)
    call fits_goto_hdu(handle%hnd,m)
    do n=1, size(handle%hnd%columns)
      handle%columns(ncols)%name=handle%hnd%columns(n)%name
      handle%columns(ncols)%unit=handle%hnd%columns(n)%unit
      handle%columns(ncols)%width=handle%hnd%columns(n)%repcount
      handle%columns(ncols)%hdu=m
      handle%columns(ncols)%idx=n
      handle%columns(ncols)%type=fitstype2type(handle%hnd%columns(n)%type)
      if (handle%columns(ncols)%type==PLANCK_STRING) then
        handle%columns(ncols)%size=handle%hnd%nrows
      else
        handle%columns(ncols)%size=handle%hnd%nrows &
                                  *handle%columns(ncols)%width
      endif
      ncols=ncols+1
    end do
  end do
end subroutine

subroutine fts_check_consistency (handle)
  type(fts_handle), intent(inout) :: handle

  integer m, hdu, idx
  character(len=200) tp

  call assert_read(handle)

  call fits_goto_hdu(handle%hnd,2)
  if (fits_key_present(handle%hnd,"objtype")) then
    call fits_get_key(handle%hnd,"objtype",tp)
    if (handle%type/=tp) &
      print *,"Warning: expected object type '",trim(handle%type), &
              "', but found '",trim(tp),"'."
  else
    print *, "Warning: FITS file does not contain type information"
  endif
  do m=0,size(handle%columns)-1
    hdu = handle%columns(m)%hdu
    idx = handle%columns(m)%idx
    call assert (fits_num_hdus(handle%hnd)>=hdu, &
      "FITS file has too few HDUs")
    call fits_goto_hdu(handle%hnd,hdu)
    call assert (size(handle%hnd%columns)>=idx, &
      "FITS HDU has too few columns")
    if (handle%columns(m)%name/=handle%hnd%columns(idx)%name) &
      print *,"Warning: HDU ",hdu,",column ",idx,": expected name '", &
              trim(handle%columns(m)%name),"', but found '", &
              trim(handle%hnd%columns(idx)%name), "'."
    if (handle%columns(m)%type/=fitstype2type(handle%hnd%columns(idx)%type)) &
      print *,"Warning: HDU ",hdu,",column ",idx,": expected type '", &
              trim(type2string(handle%columns(m)%type)),"', but found '", &
             trim(type2string(fitstype2type(handle%hnd%columns(idx)%type))),"'."
    if (handle%columns(m)%unit/='[arbitrary]') then
      if (handle%columns(m)%unit/=handle%hnd%columns(idx)%unit) &
        print *,"Warning: HDU ",hdu,",column ",idx,": expected unit '", &
                trim(handle%columns(m)%unit),"', but found '", &
                trim(handle%hnd%columns(idx)%unit), "'."
    endif
    if (handle%columns(m)%width/=handle%hnd%columns(idx)%repcount) &
      print *,"Warning: HDU ",hdu,",column ",idx,": expected repcount ", &
              handle%columns(m)%width," but found ", &
              handle%hnd%columns(idx)%repcount
  end do
end subroutine

subroutine fts_populate_fits (handle)
  type(fts_handle), intent(inout) :: handle

  integer maxhdu, maxidx, idx, m, n
  type(fitscolumn),allocatable :: tab(:)
  type(paramfile_handle) :: params

  maxhdu=0
  maxidx=0
  do m=0,size(handle%columns)-1
    if (handle%columns(m)%hdu > maxhdu) maxhdu = handle%columns(m)%hdu
    if (handle%columns(m)%idx > maxidx) maxidx = handle%columns(m)%idx
  end do

  allocate(tab(maxidx))
  do n=2,maxhdu
    maxidx=0
    do m=0,size(handle%columns)-1
      if (handle%columns(m)%hdu ==n) then
        idx = handle%columns(m)%idx
        tab(idx)%name = handle%columns(m)%name
        tab(idx)%unit = handle%columns(m)%unit
        tab(idx)%repcount = handle%columns(m)%width
        tab(idx)%type = type2fitstype(handle%columns(m)%type)
        if (idx>maxidx) maxidx=idx
      endif
    end do
    call fits_insert_bintab (handle%hnd,tab(1:maxidx))
  end do
  deallocate(tab)

!FIXME Doesn't work in a mixed scenario if the DMC holds the parameters
  if (args(1)/='') then
    params = fts_getrunparams (verbose=.false.)
    call fits_goto_hdu(handle%hnd,1)
    do n=1,size(params%keylist)
      call fits_update_key(handle%hnd,params%keylist(n),params%valuelist(n))
    end do
    call parse_finish(params, verbose=.false.)
  endif
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,"objType",handle%type)
end subroutine

subroutine fts_update_colsizes (handle)
  type(fts_handle), intent(inout) :: handle

  integer m

  call assert_read(handle)
  do m=0,size(handle%columns)-1
    call fits_goto_hdu(handle%hnd,handle%columns(m)%hdu)
    if (handle%columns(m)%type==PLANCK_STRING) then
      handle%columns(m)%size = handle%hnd%nrows
    else
      handle%columns(m)%size = handle%hnd%nrows*handle%columns(m)%width
    endif
  end do
end subroutine

#include "ddl.f90"

subroutine fts_close (handle)
  type(fts_handle), intent(inout) :: handle

  call assert_connected(handle)
  call clean_all(handle)
end subroutine

subroutine fts_flush (handle)
  type(fts_handle), intent(inout) :: handle
  ! do nothing for now
end subroutine

subroutine fts_open_generic (handle, name)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name

  call clean_all(handle)
  call fits_open(handle%hnd,name)
  handle%connected = .true.
  handle%readonly = .true.
  call fts_fill_cache (handle)
end subroutine

subroutine fts_open_with_type (handle, name, type)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name, type

  call clean_all(handle)
  handle%type = type
  call fts_fill_cache_from_ddl (handle)
  call fits_open(handle%hnd,name)
  handle%connected = .true.
  handle%readonly = .true.
  call fts_check_consistency (handle)
  call fts_update_colsizes (handle)
end subroutine

subroutine fts_create (handle, name, type)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name, type

  call clean_all(handle)
  handle%type = type
  call fts_fill_cache_from_ddl (handle)
  call fits_create(handle%hnd,name)
  handle%connected = .true.
  handle%readonly = .false.
  call fts_populate_fits(handle)
end subroutine

subroutine fts_set_key_int1 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i1b), intent(in) :: value

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,key,value)
end subroutine

subroutine fts_set_key_int2 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i2b), intent(in) :: value

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,key,value)
end subroutine

subroutine fts_set_key_int4 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i4b), intent(in) :: value

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,key,value)
end subroutine

subroutine fts_set_key_int8 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i8b), intent(in) :: value

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,key,value)
end subroutine

subroutine fts_set_key_real4 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(sp), intent(in) :: value

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,key,value)
end subroutine

subroutine fts_set_key_real8 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(dp), intent(in) :: value

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,key,value)
end subroutine

subroutine fts_set_key_bool (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical, intent(in) :: value

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,key,value)
end subroutine

subroutine fts_set_key_string (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in) :: value

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_update_key(handle%hnd,key,value)
end subroutine

subroutine fts_delete_key (handle,key)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_delete_key(handle%hnd,key)
end subroutine

function fts_key_present (handle, key)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical fts_key_present

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  fts_key_present=fits_key_present(handle%hnd,key)
end function

subroutine fts_get_key_int1 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i1b), intent(out) :: value

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_get_key(handle%hnd,key,value)
end subroutine

subroutine fts_get_key_int2 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i2b), intent(out) :: value

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_get_key(handle%hnd,key,value)
end subroutine

subroutine fts_get_key_int4 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i4b), intent(out) :: value

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_get_key(handle%hnd,key,value)
end subroutine

subroutine fts_get_key_int8 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i8b), intent(out) :: value

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_get_key(handle%hnd,key,value)
end subroutine

subroutine fts_get_key_real4 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(sp), intent(out) :: value

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_get_key(handle%hnd,key,value)
end subroutine

subroutine fts_get_key_real8 (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(dp), intent(out) :: value

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_get_key(handle%hnd,key,value)
end subroutine

subroutine fts_get_key_bool (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical, intent(out) :: value

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_get_key(handle%hnd,key,value)
end subroutine

subroutine fts_get_key_string (handle,key,value)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(out) :: value

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,2)
  call fits_get_key(handle%hnd,key,value)
end subroutine

subroutine fts_append_column_real4 (handle, colnum, data)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(sp), intent(in) :: data(:)

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_write_column(handle%hnd,handle%columns(colnum)%idx,data, &
                         handle%columns(colnum)%size)
  handle%columns(colnum)%size=handle%columns(colnum)%size+size(data)
end subroutine

subroutine fts_append_column_real8 (handle, colnum, data)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(dp), intent(in) :: data(:)

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_write_column(handle%hnd,handle%columns(colnum)%idx,data, &
                         handle%columns(colnum)%size)
  handle%columns(colnum)%size=handle%columns(colnum)%size+size(data)
end subroutine

subroutine fts_append_column_int1 (handle, colnum, data)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i1b), intent(in) :: data(:)

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_write_column(handle%hnd,handle%columns(colnum)%idx,data, &
                         handle%columns(colnum)%size)
  handle%columns(colnum)%size=handle%columns(colnum)%size+size(data)
end subroutine

subroutine fts_append_column_int2 (handle, colnum, data)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i2b), intent(in) :: data(:)

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_write_column(handle%hnd,handle%columns(colnum)%idx,data, &
                         handle%columns(colnum)%size)
  handle%columns(colnum)%size=handle%columns(colnum)%size+size(data)
end subroutine

subroutine fts_append_column_int4 (handle, colnum, data)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b), intent(in) :: data(:)

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_write_column(handle%hnd,handle%columns(colnum)%idx,data, &
                         handle%columns(colnum)%size)
  handle%columns(colnum)%size=handle%columns(colnum)%size+size(data)
end subroutine

subroutine fts_append_column_int8 (handle, colnum, data)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(in) :: data(:)

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_write_column(handle%hnd,handle%columns(colnum)%idx,data, &
                         handle%columns(colnum)%size)
  handle%columns(colnum)%size=handle%columns(colnum)%size+size(data)
end subroutine

subroutine fts_append_column_string (handle, colnum, data)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(in) :: data(:)

  call assert_write(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_write_column(handle%hnd,handle%columns(colnum)%idx,data, &
                         handle%columns(colnum)%size)
  handle%columns(colnum)%size=handle%columns(colnum)%size+size(data)
end subroutine

subroutine fts_read_column_real4 (handle, colnum, data, offset)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(sp), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_read_column(handle%hnd,handle%columns(colnum)%idx,data,offset)
end subroutine

subroutine fts_read_column_real8 (handle, colnum, data, offset)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(dp), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_read_column(handle%hnd,handle%columns(colnum)%idx,data,offset)
end subroutine

subroutine fts_read_column_int1 (handle, colnum, data, offset)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i1b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_read_column(handle%hnd,handle%columns(colnum)%idx,data,offset)
end subroutine

subroutine fts_read_column_int2 (handle, colnum, data, offset)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i2b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_read_column(handle%hnd,handle%columns(colnum)%idx,data,offset)
end subroutine

subroutine fts_read_column_int4 (handle, colnum, data, offset)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_read_column(handle%hnd,handle%columns(colnum)%idx,data,offset)
end subroutine

subroutine fts_read_column_int8 (handle, colnum, data, offset)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_read_column(handle%hnd,handle%columns(colnum)%idx,data,offset)
end subroutine

subroutine fts_read_column_string (handle, colnum, data, offset)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  call assert_read(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  call fits_read_column(handle%hnd,handle%columns(colnum)%idx,data,offset)
end subroutine

subroutine fts_reserve_column (handle, colnum, sz)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(in) :: sz

  integer dummy(1)

  call assert_write(handle)
  call assert (handle%columns(colnum)%size==0,"cannot reserve after writing")
  call assert (sz>=1,"size must be at least 1")
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  dummy=0
  call fits_write_column(handle%hnd,handle%columns(colnum)%idx,dummy,sz-1)
end subroutine

function fts_column_present (handle, colname)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  logical fts_column_present

  integer m

  call assert_connected(handle)
  do m=0,size(handle%columns)-1
    if (handle%columns(m)%name==colname) then
      fts_column_present=.true.
      return
    endif
  end do

  fts_column_present=.false.
end function

function fts_colnum (handle, colname)
  type(fts_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer fts_colnum

  integer m

  call assert_connected(handle)
  do m=0,size(handle%columns)-1
    if (handle%columns(m)%name==colname) then
      fts_colnum=m
      return
    endif
  end do

  fts_colnum=-1
  call fatal_error("column name not found")
end function

function fts_collength (handle, colnum)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b) fts_collength

  call assert_read(handle)
  fts_collength = handle%columns(colnum)%size
end function

function fts_numcols (handle)
  type(fts_handle), intent(inout) :: handle
  integer(i4b) fts_numcols

  call assert_connected(handle)
  fts_numcols = size(handle%columns)
end function

subroutine fts_get_colname (handle, colnum, colname)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: colname

  call assert_connected(handle)
  colname = handle%columns(colnum)%name
end subroutine

subroutine fts_get_colunit (handle, colnum, colunit)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: colunit

  call assert_connected(handle)
  colunit = handle%columns(colnum)%unit
end subroutine

function fts_coltype (handle, colnum)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b) :: fts_coltype

  call assert_connected(handle)
  fts_coltype = handle%columns(colnum)%type
end function

function fts_eff_chunksize (handle, colnum)
  type(fts_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b) :: fts_eff_chunksize

  call assert_connected(handle)
  call fits_goto_hdu(handle%hnd,handle%columns(colnum)%hdu)
  fts_eff_chunksize = fits_fast_nelms(handle%hnd,handle%columns(colnum)%idx)
end function

end module fts_io

module dmc_io
use planck_config
use ls_paramfile_io
use fts_io
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
use cll_io
#endif
implicit none
private

#ifdef USE_TOODI
character(len=*), parameter, public :: backend_type = "TOODI"
#else
#ifdef USE_HFIDMC
character(len=*), parameter, public :: backend_type = "HFIDMC"
#else
character(len=*), parameter, public :: backend_type = "FITS"
#endif
#endif

public dmc_handle, dmc_open, dmc_close, dmc_flush, &
  dmc_create, dmc_set_key, dmc_get_key, dmc_delete_key, &
  dmc_key_present, dmc_append_column, dmc_read_column, dmc_reserve_column, &
  dmc_colnum, dmc_collength, dmc_init, dmc_shutdown, dmc_getrunparams, &
  dmc_numcols, dmc_get_colname, dmc_column_present, dmc_get_colunit, &
  dmc_coltype, dmc_eff_chunksize

interface dmc_init
  module procedure dmc_init_old, dmc_init_new
end interface

interface dmc_open
  module procedure dmc_open_generic, dmc_open_with_type
end interface

interface dmc_set_key
  module procedure dmc_set_key_int1, dmc_set_key_int2, dmc_set_key_int4, &
    dmc_set_key_int8, dmc_set_key_real4, dmc_set_key_real8, dmc_set_key_bool, &
    dmc_set_key_string
end interface

interface dmc_get_key
  module procedure dmc_get_key_int1, dmc_get_key_int2, dmc_get_key_int4, &
    dmc_get_key_int8, dmc_get_key_real4, dmc_get_key_real8, dmc_get_key_bool, &
    dmc_get_key_string
end interface

interface dmc_append_column
  module procedure dmc_append_column_real4, dmc_append_column_real8, &
    dmc_append_column_int1, dmc_append_column_int2, dmc_append_column_int4, &
    dmc_append_column_int8, dmc_append_column_string, &
    dmc_append_column_real4_s, dmc_append_column_real8_s, &
    dmc_append_column_int1_s, dmc_append_column_int2_s, &
    dmc_append_column_int4_s, dmc_append_column_int8_s, &
    dmc_append_column_string_s
end interface

interface dmc_read_column
  module procedure dmc_read_column_real4, dmc_read_column_real8, &
    dmc_read_column_int1, dmc_read_column_int2, dmc_read_column_int4, &
    dmc_read_column_int8, dmc_read_column_string, &
    dmc_read_column_real4_s, dmc_read_column_real8_s, &
    dmc_read_column_int1_s, dmc_read_column_int2_s, dmc_read_column_int4_s, &
    dmc_read_column_int8_s, dmc_read_column_string_s
end interface

interface dmc_reserve_column
  module procedure dmc_reserve_column_n, dmc_reserve_column_s
end interface

interface dmc_collength
  module procedure dmc_collength_n, dmc_collength_s
end interface

interface dmc_get_colunit
  module procedure dmc_get_colunit_n, dmc_get_colunit_s
end interface

interface dmc_coltype
  module procedure dmc_coltype_n, dmc_coltype_s
end interface

interface dmc_eff_chunksize
  module procedure dmc_eff_chunksize_n, dmc_eff_chunksize_s
end interface

type dmc_handle
  type(fts_handle) :: fh
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  type(cll_handle) :: ch
#endif
  logical is_fits
end type

contains

subroutine prepare_name (name, name2, isfits)
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: name2
  logical, intent(out) :: isfits

#if defined(USE_TOODI)
  if (index (name,"TOODI:")==1) then
    name2=name(7:len(name))
    isfits=.false.
    return
  endif
  if (index (name,"FITS:")/=1) then
    name2=name
    isfits=.false.
    return
  endif
#endif
#if defined(USE_HFIDMC)
  if (index (name,"HFIDMC:")==1) then
    name2=name(8:len(name))
    isfits=.false.
    return
  endif
  if (index (name,"FITS:")/=1) then
    name2=name
    isfits=.false.
    return
  endif
#endif
  if (index (name,"FITS:")==1) then
    name2=name(6:len(name))
    isfits=.true.
    return
  endif
  name2=name
  isfits=.true.
end subroutine

subroutine dmc_init_new (args)
  character(len=*), intent(in) :: args(:)

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  call cll_init(args)
  call fts_init( (/''/) )
#else
  call fts_init(args)
#endif
end subroutine

subroutine dmc_init_old (args)
  character(len=*), intent(in) :: args

  call dmc_init_new( (/args/) )
end subroutine

subroutine dmc_shutdown ()
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  call cll_shutdown
#else
  call fts_shutdown
#endif
end subroutine

function dmc_getrunparams (verbose)
  logical, intent(in), optional :: verbose
  type(paramfile_handle) :: dmc_getrunparams

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  dmc_getrunparams = cll_getrunparams(verbose)
#else
  dmc_getrunparams = fts_getrunparams(verbose)
#endif
end function

subroutine dmc_close (handle)
  type(dmc_handle), intent(inout) :: handle

  if (handle%is_fits) then
    call fts_close(handle%fh)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_close(handle%ch)
#endif
  endif
end subroutine

subroutine dmc_flush (handle)
  type(dmc_handle), intent(inout) :: handle

  if (handle%is_fits) then
    call fts_flush(handle%fh)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_flush(handle%ch)
#endif
  endif
end subroutine

subroutine dmc_open_generic (handle, name)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name

  character(len(name)) :: name2

  call prepare_name (name, name2, handle%is_fits)
  if (handle%is_fits) then
    call fts_open(handle%fh,name2)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_open(handle%ch,name2)
#endif
  endif
end subroutine

subroutine dmc_open_with_type (handle, name, type)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name, type

  character(len(name)) :: name2

  call prepare_name (name, name2, handle%is_fits)
  if (handle%is_fits) then
    call fts_open(handle%fh,name2,type)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_open(handle%ch,name2,type)
#endif
  endif
end subroutine

subroutine dmc_create (handle, name, type)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name, type

  character(len(name)) :: name2

  call prepare_name (name, name2, handle%is_fits)
  if (handle%is_fits) then
    call fts_create(handle%fh,name2,type)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_create(handle%ch,name2,type)
#endif
  endif
end subroutine

subroutine dmc_set_key_int1 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i1b), intent(in) :: value

  if (handle%is_fits) then
    call fts_set_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_set_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_set_key_int2 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i2b), intent(in) :: value

  if (handle%is_fits) then
    call fts_set_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_set_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_set_key_int4 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i4b), intent(in) :: value

  if (handle%is_fits) then
    call fts_set_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_set_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_set_key_int8 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i8b), intent(in) :: value

  if (handle%is_fits) then
    call fts_set_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_set_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_set_key_real4 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(sp), intent(in) :: value

  if (handle%is_fits) then
    call fts_set_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_set_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_set_key_real8 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(dp), intent(in) :: value

  if (handle%is_fits) then
    call fts_set_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_set_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_set_key_bool (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical, intent(in) :: value

  if (handle%is_fits) then
    call fts_set_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_set_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_set_key_string (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in) :: value

  if (handle%is_fits) then
    call fts_set_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_set_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_get_key_int1 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i1b), intent(out) :: value

  if (handle%is_fits) then
    call fts_get_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_get_key_int2 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i2b), intent(out) :: value

  if (handle%is_fits) then
    call fts_get_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_get_key_int4 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i4b), intent(out) :: value

  if (handle%is_fits) then
    call fts_get_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_get_key_int8 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i8b), intent(out) :: value

  if (handle%is_fits) then
    call fts_get_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_get_key_real4 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(sp), intent(out) :: value

  if (handle%is_fits) then
    call fts_get_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_get_key_real8 (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(dp), intent(out) :: value

  if (handle%is_fits) then
    call fts_get_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_get_key_bool (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical, intent(out) :: value

  if (handle%is_fits) then
    call fts_get_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_get_key_string (handle,key,value)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(out) :: value

  if (handle%is_fits) then
    call fts_get_key(handle%fh,key,value)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_key(handle%ch,key,value)
#endif
  endif
end subroutine

subroutine dmc_delete_key (handle,key)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key

  if (handle%is_fits) then
    call fts_delete_key(handle%fh,key)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_delete_key(handle%ch,key)
#endif
  endif
end subroutine

function dmc_key_present (handle, key)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical dmc_key_present

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  if (handle%is_fits) then
    dmc_key_present=fts_key_present(handle%fh,key)
  else
    dmc_key_present=cll_key_present(handle%ch,key)
  endif
#else
  dmc_key_present=fts_key_present(handle%fh,key)
#endif
end function

subroutine dmc_append_column_real4 (handle, colnum, data)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(sp), intent(in) :: data(:)

  if (handle%is_fits) then
    call fts_append_column(handle%fh,colnum,data)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_append_column(handle%ch,colnum,data)
#endif
  endif
end subroutine

subroutine dmc_append_column_real8 (handle, colnum, data)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(dp), intent(in) :: data(:)

  if (handle%is_fits) then
    call fts_append_column(handle%fh,colnum,data)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_append_column(handle%ch,colnum,data)
#endif
  endif
end subroutine

subroutine dmc_append_column_int1 (handle, colnum, data)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i1b), intent(in) :: data(:)

  if (handle%is_fits) then
    call fts_append_column(handle%fh,colnum,data)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_append_column(handle%ch,colnum,data)
#endif
  endif
end subroutine

subroutine dmc_append_column_int2 (handle, colnum, data)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i2b), intent(in) :: data(:)

  if (handle%is_fits) then
    call fts_append_column(handle%fh,colnum,data)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_append_column(handle%ch,colnum,data)
#endif
  endif
end subroutine

subroutine dmc_append_column_int4 (handle, colnum, data)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b), intent(in) :: data(:)

  if (handle%is_fits) then
    call fts_append_column(handle%fh,colnum,data)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_append_column(handle%ch,colnum,data)
#endif
  endif
end subroutine

subroutine dmc_append_column_int8 (handle, colnum, data)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(in) :: data(:)

  if (handle%is_fits) then
    call fts_append_column(handle%fh,colnum,data)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_append_column(handle%ch,colnum,data)
#endif
  endif
end subroutine

subroutine dmc_append_column_string (handle, colnum, data)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(in) :: data(:)

  if (handle%is_fits) then
    call fts_append_column(handle%fh,colnum,data)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_append_column(handle%ch,colnum,data)
#endif
  endif
end subroutine

subroutine dmc_read_column_real4 (handle, colnum, data, offset)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(sp), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  if (handle%is_fits) then
    call fts_read_column(handle%fh,colnum,data,offset)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_read_column(handle%ch,colnum,data,offset)
#endif
  endif
end subroutine

subroutine dmc_read_column_real8 (handle, colnum, data, offset)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(dp), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  if (handle%is_fits) then
    call fts_read_column(handle%fh,colnum,data,offset)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_read_column(handle%ch,colnum,data,offset)
#endif
  endif
end subroutine

subroutine dmc_read_column_int1 (handle, colnum, data, offset)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i1b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  if (handle%is_fits) then
    call fts_read_column(handle%fh,colnum,data,offset)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_read_column(handle%ch,colnum,data,offset)
#endif
  endif
end subroutine

subroutine dmc_read_column_int2 (handle, colnum, data, offset)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i2b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  if (handle%is_fits) then
    call fts_read_column(handle%fh,colnum,data,offset)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_read_column(handle%ch,colnum,data,offset)
#endif
  endif
end subroutine

subroutine dmc_read_column_int4 (handle, colnum, data, offset)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  if (handle%is_fits) then
    call fts_read_column(handle%fh,colnum,data,offset)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_read_column(handle%ch,colnum,data,offset)
#endif
  endif
end subroutine

subroutine dmc_read_column_int8 (handle, colnum, data, offset)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  if (handle%is_fits) then
    call fts_read_column(handle%fh,colnum,data,offset)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_read_column(handle%ch,colnum,data,offset)
#endif
  endif
end subroutine

subroutine dmc_read_column_string (handle, colnum, data, offset)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset

  if (handle%is_fits) then
    call fts_read_column(handle%fh,colnum,data,offset)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_read_column(handle%ch,colnum,data,offset)
#endif
  endif
end subroutine

subroutine dmc_reserve_column_n (handle, colnum, sz)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(in) :: sz

  if (handle%is_fits) call fts_reserve_column(handle%fh,colnum,sz)
end subroutine

function dmc_column_present (handle, colname)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  logical dmc_column_present

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  if (handle%is_fits) then
    dmc_column_present=fts_column_present(handle%fh,colname)
  else
    dmc_column_present=cll_column_present(handle%ch,colname)
  endif
#else
  dmc_column_present=fts_column_present(handle%fh,colname)
#endif
end function

function dmc_colnum (handle, colname)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer dmc_colnum

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  if (handle%is_fits) then
    dmc_colnum=fts_colnum(handle%fh,colname)
  else
    dmc_colnum=cll_colnum(handle%ch,colname)
  endif
#else
  dmc_colnum=fts_colnum(handle%fh,colname)
#endif
end function

function dmc_collength_n (handle, colnum)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b) dmc_collength_n

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  if (handle%is_fits) then
    dmc_collength_n=fts_collength(handle%fh,colnum)
  else
    dmc_collength_n=cll_collength(handle%ch,colnum)
  endif
#else
  dmc_collength_n=fts_collength(handle%fh,colnum)
#endif
end function

function dmc_numcols (handle)
  type(dmc_handle), intent(inout) :: handle
  integer(i4b) dmc_numcols

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  if (handle%is_fits) then
    dmc_numcols=fts_numcols(handle%fh)
  else
    dmc_numcols=cll_numcols(handle%ch)
  endif
#else
  dmc_numcols=fts_numcols(handle%fh)
#endif
end function

subroutine dmc_get_colname (handle, colnum, colname)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: colname

  if (handle%is_fits) then
    call fts_get_colname(handle%fh,colnum,colname)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_colname(handle%ch,colnum,colname)
#endif
  endif
end subroutine

subroutine dmc_get_colunit_n(handle, colnum, colunit)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: colunit

  if (handle%is_fits) then
    call fts_get_colunit(handle%fh,colnum,colunit)
#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  else
    call cll_get_colunit(handle%ch,colnum,colunit)
#endif
  endif
end subroutine

function dmc_coltype_n (handle, colnum)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b) :: dmc_coltype_n

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  if (handle%is_fits) then
    dmc_coltype_n=fts_coltype(handle%fh,colnum)
  else
    dmc_coltype_n=cll_coltype(handle%ch,colnum)
  endif
#else
  dmc_coltype_n=fts_coltype(handle%fh,colnum)
#endif
end function

function dmc_eff_chunksize_n (handle, colnum)
  type(dmc_handle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b) :: dmc_eff_chunksize_n

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))
  if (handle%is_fits) then
    dmc_eff_chunksize_n=fts_eff_chunksize(handle%fh,colnum)
  else
    dmc_eff_chunksize_n=cll_eff_chunksize(handle%ch,colnum)
  endif
#else
  dmc_eff_chunksize_n=fts_eff_chunksize(handle%fh,colnum)
#endif
end function

subroutine dmc_append_column_real4_s (handle, colname, data)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  real(sp), intent(in) :: data(:)
  call dmc_append_column(handle,dmc_colnum(handle,colname),data)
end subroutine
subroutine dmc_append_column_real8_s (handle, colname, data)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  real(dp), intent(in) :: data(:)
  call dmc_append_column(handle,dmc_colnum(handle,colname),data)
end subroutine
subroutine dmc_append_column_int1_s (handle, colname, data)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i1b), intent(in) :: data(:)
  call dmc_append_column(handle,dmc_colnum(handle,colname),data)
end subroutine
subroutine dmc_append_column_int2_s (handle, colname, data)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i2b), intent(in) :: data(:)
  call dmc_append_column(handle,dmc_colnum(handle,colname),data)
end subroutine
subroutine dmc_append_column_int4_s (handle, colname, data)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i4b), intent(in) :: data(:)
  call dmc_append_column(handle,dmc_colnum(handle,colname),data)
end subroutine
subroutine dmc_append_column_int8_s (handle, colname, data)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i8b), intent(in) :: data(:)
  call dmc_append_column(handle,dmc_colnum(handle,colname),data)
end subroutine
subroutine dmc_append_column_string_s (handle, colname, data)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  character(len=*), intent(in) :: data(:)
  call dmc_append_column(handle,dmc_colnum(handle,colname),data)
end subroutine

subroutine dmc_read_column_real4_s (handle, colname, data, offset)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  real(sp), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset
  call dmc_read_column(handle,dmc_colnum(handle,colname),data,offset)
end subroutine
subroutine dmc_read_column_real8_s (handle, colname, data, offset)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  real(dp), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset
  call dmc_read_column(handle,dmc_colnum(handle,colname),data,offset)
end subroutine
subroutine dmc_read_column_int1_s (handle, colname, data, offset)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i1b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset
  call dmc_read_column(handle,dmc_colnum(handle,colname),data,offset)
end subroutine
subroutine dmc_read_column_int2_s (handle, colname, data, offset)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i2b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset
  call dmc_read_column(handle,dmc_colnum(handle,colname),data,offset)
end subroutine
subroutine dmc_read_column_int4_s (handle, colname, data, offset)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i4b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset
  call dmc_read_column(handle,dmc_colnum(handle,colname),data,offset)
end subroutine
subroutine dmc_read_column_int8_s (handle, colname, data, offset)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i8b), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset
  call dmc_read_column(handle,dmc_colnum(handle,colname),data,offset)
end subroutine
subroutine dmc_read_column_string_s (handle, colname, data, offset)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  character(len=*), intent(out) :: data(:)
  integer(i8b), intent(in) :: offset
  call dmc_read_column(handle,dmc_colnum(handle,colname),data,offset)
end subroutine

subroutine dmc_reserve_column_s (handle, colname, sz)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i8b), intent(in) :: sz

  call dmc_reserve_column_n(handle,dmc_colnum(handle,colname),sz)
end subroutine

function dmc_collength_s (handle, colname)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i8b) dmc_collength_s

  dmc_collength_s=dmc_collength(handle,dmc_colnum(handle,colname))
end function

subroutine dmc_get_colunit_s (handle, colname, colunit)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  character(len=*), intent(out) :: colunit
  call dmc_get_colunit(handle,dmc_colnum(handle,colname),colunit)
end subroutine

function dmc_coltype_s (handle, colname)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i4b) :: dmc_coltype_s
  dmc_coltype_s=dmc_coltype(handle,dmc_colnum(handle,colname))
end function

function dmc_eff_chunksize_s (handle, colname)
  type(dmc_handle), intent(inout) :: handle
  character(len=*), intent(in) :: colname
  integer(i8b) dmc_eff_chunksize_s

  dmc_eff_chunksize_s=dmc_eff_chunksize(handle,dmc_colnum(handle,colname))
end function

end module dmc_io
