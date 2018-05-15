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

module focalplanemod
use planck_config
use instrument_table
implicit none
private

public fpdb_open, fpdb_close, fpdb_quantity_present, &
       fpdb_value_real8, fpdb_value_int4

character(len=100), public, save :: fpdb_version
real(dp), public, save :: fpdb_theta_b
type(itable_handle), save :: inp

contains

subroutine fpdb_open (file)
  character(len=*), intent(in) :: file

! FIXME: in principle this should better be opened without explicitly
!        specifying a DDL type. HFIDMC currently does not allow this,
!        however ...
  inp = itable_open (file,'focalplane.LS_focalplanedb','detector')
  fpdb_version = trim(itable_metavalue_string(inp,'fpdb_version'))
  fpdb_theta_b = itable_metavalue_real8(inp,'theta_b')
end subroutine

subroutine fpdb_close
  call itable_close(inp)
end subroutine

function fpdb_quantity_present(quantity)
  character(len=*), intent(in) :: quantity
  logical :: fpdb_quantity_present

  fpdb_quantity_present = itable_quantity_present(inp,quantity)
end function

function fpdb_value_real8 (name, quantity)
  character(len=*), intent(in) :: name, quantity
  real(dp) :: fpdb_value_real8

  fpdb_value_real8 = itable_value_real8(inp,name,quantity)
end function

function fpdb_value_int4 (name, quantity)
  character(len=*), intent(in) :: name, quantity
  integer(i4b) :: fpdb_value_int4

  fpdb_value_int4 = itable_value_int4(inp,name,quantity)
end function

end module focalplanemod
