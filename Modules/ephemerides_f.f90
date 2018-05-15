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

!  Copyright (C) 2010-2013 Max-Planck-Society
!  \author Martin Reinecke

module ephemerides_f
use planck_config
use ls_misc_utils
use dmc_io
implicit none
private

! This type is supposed to be a black box for the user
! User code must not read any data from inside this type
type eph_handle
  type(dmc_handle) :: hnd
  character(len=32), dimension(:), pointer :: &
    body=>NULL(), scalar=>NULL(), array=>NULL()
  real(dp), dimension(:), pointer :: data=>NULL()
  integer(i8b) nsamp
  real(dp) t0, dt
  integer datacol
  integer(i8b) o_x, o_y, o_z, o_angdiam, o_dsun, o_sto
end type

public eph_handle, eph_open, eph_close, eph_load_body, eph_posRelSat_m

contains

function arr_index (arr,val)
  character(len=*), dimension(:), intent(in) :: arr
  character(len=*), intent(in) :: val
  integer arr_index

  do arr_index = 1, size(arr)
    if (val==arr(arr_index)) return
  end do
  call exit_with_status(1,"value not found in array")
end function

subroutine eph_open (handle, name)
  type(eph_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name

  integer :: colnum
  integer(i8b) :: num

  call dmc_open (handle%hnd, name, 'ephemeris.LS_ephemeris')
  call dmc_get_key (handle%hnd, 'nsamples', handle%nsamp)
  call dmc_get_key (handle%hnd, 'T0', handle%t0)
  call dmc_get_key (handle%hnd, 'deltaT', handle%dt)
  colnum = dmc_colnum(handle%hnd,'bodies')
  num = dmc_collength(handle%hnd,colnum)
  call assert (num>0,"no bodies found in ephemeris object")
  allocate (handle%body(num))
  call dmc_read_column(handle%hnd,colnum,handle%body,0_i8b)
  colnum = dmc_colnum(handle%hnd,'scalars')
  num = dmc_collength(handle%hnd,colnum)
  allocate (handle%scalar(num))
  call dmc_read_column(handle%hnd,colnum,handle%scalar,0_i8b)
  colnum = dmc_colnum(handle%hnd,'arrays')
  num = dmc_collength(handle%hnd,colnum)
  allocate (handle%array(num))
  call dmc_read_column(handle%hnd,colnum,handle%array,0_i8b)

  num = size(handle%scalar)+handle%nsamp*size(handle%array)
  call assert (num>0, "no data found in ephemeris object")

  allocate (handle%data(num))
  handle%datacol = dmc_colnum(handle%hnd,'data')

  handle%o_x = eph_get_array_offset(handle,"x_rel_sat_ecl/m")
  handle%o_y = eph_get_array_offset(handle,"y_rel_sat_ecl/m")
  handle%o_z = eph_get_array_offset(handle,"z_rel_sat_ecl/m")
  handle%o_angdiam = eph_get_array_offset(handle,"angdiam_rel_sat/arcsec")
  handle%o_dsun = eph_get_array_offset(handle,"d_sun/m")
  handle%o_sto = eph_get_array_offset(handle,"angle_S_T_O/deg")
end subroutine

subroutine eph_close (handle)
  type(eph_handle), intent(inout) :: handle

  call dmc_close (handle%hnd)
  deallocate (handle%body, handle%scalar, handle%array, handle%data)
end subroutine

subroutine eph_load_body (handle, name)
  type(eph_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name

  integer(i8b) offset

  offset = size(handle%data)*(arr_index(handle%body,name)-1)
  call dmc_read_column (handle%hnd,handle%datacol,handle%data,offset)
end subroutine

function eph_get_scalar_offset (handle, name)
  type(eph_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name
  integer eph_get_scalar_offset

  eph_get_scalar_offset = arr_index(handle%scalar,name)-1
end function

function eph_get_array_offset (handle, name)
  type(eph_handle), intent(inout) :: handle
  character(len=*), intent(in) :: name
  integer(i8b) eph_get_array_offset

  eph_get_array_offset = &
    size(handle%scalar) + handle%nsamp*(arr_index(handle%array,name)-1)
end function

function eph_get_scalar (handle, offset)
  type(eph_handle), intent(inout) :: handle
  integer, intent(in) :: offset
  real(dp) eph_get_scalar

  eph_get_scalar = handle%data(1+offset)
end function

function eph_get_array_value (handle, offset, time)
  type(eph_handle), intent(inout) :: handle
  integer, intent(in) :: offset
  real(dp), intent(in) :: time
  real(dp) eph_get_array_value

  real(dp) t2, frac
  integer idx

  t2 = time - handle%t0
  call assert(t2>=0,'requested time < t0')
  frac = t2/handle%dt
  idx = int(frac)
  call assert(idx<handle%nsamp-1,'requested time too large')
  frac = frac-idx
  idx = idx + offset + 1
  eph_get_array_value = (1._dp-frac)*handle%data(idx) + frac*handle%data(idx+1)
end function

function eph_posRelSat_m (handle, time)
  type(eph_handle), intent(inout) :: handle
  real(dp), intent(in) :: time
  real(dp), dimension(3) :: eph_posRelSat_m

  real(dp) t2, frac
  integer(i8b) idx, i2

  t2 = time - handle%t0
  call assert(t2>=0,'requested time < t0')
  frac = t2/handle%dt
  idx = int(frac)
  call assert(idx<handle%nsamp-1,'requested time too large')
  frac = frac-idx
  i2 = idx + handle%o_x + 1
  eph_posRelSat_m(1) = (1._dp-frac)*handle%data(i2) + frac*handle%data(i2+1)
  i2 = idx + handle%o_y + 1
  eph_posRelSat_m(2) = (1._dp-frac)*handle%data(i2) + frac*handle%data(i2+1)
  i2 = idx + handle%o_z + 1
  eph_posRelSat_m(3) = (1._dp-frac)*handle%data(i2) + frac*handle%data(i2+1)
end function

end module ephemerides_f
