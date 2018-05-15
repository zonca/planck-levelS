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

!  Copyright (C) 2006-2013 Max-Planck-Society
!  \author Martin Reinecke

module planck_types
use planck_config
implicit none
private

integer(i4b), parameter, public :: PLANCK_BOOL    = 0, &
                                   PLANCK_INT8    = 1, &
                                   PLANCK_INT16   = 2, &
                                   PLANCK_INT32   = 3, &
                                   PLANCK_INT64   = 4, &
                                   PLANCK_FLOAT32 = 5, &
                                   PLANCK_FLOAT64 = 6, &
                                   PLANCK_STRING  = 7, &
                                   PLANCK_INVALID_TYPE = -1

public string2type, type2string

contains

function string2type (stype)
  character(len=*), intent(in) :: stype
  integer(i4b) :: string2type

  if (stype=="INT8") then
    string2type = PLANCK_INT8
  elseif (stype=="INT16") then
    string2type = PLANCK_INT16
  elseif (stype=="INT32") then
    string2type = PLANCK_INT32
  elseif (stype=="INT64") then
    string2type = PLANCK_INT64
  elseif (stype=="FLOAT32") then
    string2type = PLANCK_FLOAT32
  elseif (stype=="FLOAT64") then
    string2type = PLANCK_FLOAT64
  elseif (stype=="BOOL") then
    string2type = PLANCK_BOOL
  elseif (stype=="STRING") then
    string2type = PLANCK_STRING
  else
    string2type = PLANCK_INVALID_TYPE
    call exit_with_status(1,"string2type: unknown type")
  endif
end function

function type2string (type)
  integer(i4b), intent(in) :: type
  character(len=10) :: type2string

  select case(type)
    case (PLANCK_INT8)
      type2string = "INT8"
    case (PLANCK_INT16)
      type2string = "INT16"
    case (PLANCK_INT32)
      type2string = "INT32"
    case (PLANCK_INT64)
      type2string = "INT64"
    case (PLANCK_FLOAT32)
      type2string = "FLOAT32"
    case (PLANCK_FLOAT64)
      type2string = "FLOAT64"
    case (PLANCK_BOOL)
      type2string = "BOOL"
    case (PLANCK_STRING)
      type2string = "STRING"
    case default
      call exit_with_status(1,"type2string: unknown type")
  end select
end function

end module
