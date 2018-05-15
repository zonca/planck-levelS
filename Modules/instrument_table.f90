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

!  Copyright (C) 2011-2013 Max-Planck-Society
!  \author Martin Reinecke

module instrument_table
use planck_config
use ls_misc_utils
use dmc_io
implicit none
private

public itable_open, itable_close, itable_quantity_present, &
       itable_quantity_unit, itable_value_real8, itable_value_int4, &
       itable_metavalue_present, itable_metavalue_real8, &
       itable_metavalue_string, itable_handle

type itable_handle
  character(len=100), allocatable, dimension(:) :: keys
  type(dmc_handle) :: inp
end type

contains

function itable_open (file,type,keycolumn) result(handle)
  character(len=*), intent(in) :: file,type,keycolumn
  type(itable_handle) :: handle

  if (type=="") then
    call dmc_open (handle%inp, file)
  else
    call dmc_open (handle%inp, file, type)
  endif
  allocate(handle%keys(0:dmc_collength(handle%inp,keycolumn)-1))
  call dmc_read_column (handle%inp, keycolumn, handle%keys, 0_i8b)
end function

subroutine itable_close(handle)
  type(itable_handle), intent(inout) :: handle
  call dmc_close(handle%inp)
  deallocate(handle%keys)
end subroutine

function itable_quantity_present(handle,quantity)
  type(itable_handle), intent(inout) :: handle
  character(len=*), intent(in) :: quantity
  logical :: itable_quantity_present

  itable_quantity_present = dmc_column_present(handle%inp,quantity)
end function

function itable_quantity_unit(handle,quantity)
  type(itable_handle), intent(inout) :: handle
  character(len=*), intent(in) :: quantity
  character(len=filenamelen) :: itable_quantity_unit

  call dmc_get_colunit(handle%inp,quantity,itable_quantity_unit)
end function

function itable_value_real8 (handle, name, quantity)
  type(itable_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name, quantity
  real(dp) :: itable_value_real8

  integer(i8b) m
  real(dp) :: tmp(1)

  do m=0,size(handle%keys)-1
    if (name==handle%keys(m)) then
      call dmc_read_column (handle%inp, dmc_colnum(handle%inp,quantity), tmp, m)
      itable_value_real8=tmp(1)
      return
    endif
  end do
  itable_value_real8=0
  call fatal_error ("requested entity not found in instrument table")
end function

function itable_value_int4 (handle, name, quantity)
  type(itable_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name, quantity
  integer(i4b) :: itable_value_int4

  integer(i8b) m
  integer(i4b) :: tmp(1)

  do m=0,size(handle%keys)-1
    if (name==handle%keys(m)) then
      call dmc_read_column (handle%inp, dmc_colnum(handle%inp,quantity), tmp, m)
      itable_value_int4=tmp(1)
      return
    endif
  end do
  itable_value_int4=0
  call fatal_error ("requested entity not found in instrument table")
end function

function itable_metavalue_present (handle, quantity)
  type(itable_handle), intent(inout) :: handle
  character(len=*), intent(in) :: quantity
  logical :: itable_metavalue_present

  itable_metavalue_present = dmc_key_present(handle%inp, quantity)
end function

function itable_metavalue_real8 (handle, quantity)
  type(itable_handle), intent(inout) :: handle
  character(len=*), intent(in) :: quantity
  real(dp) :: itable_metavalue_real8

  call dmc_get_key(handle%inp, quantity, itable_metavalue_real8)
end function

function itable_metavalue_string (handle, quantity)
  type(itable_handle), intent(inout) :: handle
  character(len=*), intent(in) :: quantity
  character(len=filenamelen) :: itable_metavalue_string

  call dmc_get_key(handle%inp, quantity, itable_metavalue_string)
end function

end module instrument_table
