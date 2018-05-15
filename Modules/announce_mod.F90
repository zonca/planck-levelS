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

!  Copyright (C) 2008-2013 Max-Planck-Society
!  \author Martin Reinecke

module announce_mod
use planck_config
use ls_misc_utils
use dmc_io
implicit none
private

public module_startup, announce

contains

subroutine openmp_status
#ifdef _OPENMP
  integer threads
  interface
    integer function omp_get_max_threads()
    end function
  end interface
  threads = omp_get_max_threads()
  if (threads == 1) then
    print *, "OpenMP active, but running with 1 thread only."
  else
    print *, "OpenMP active: max. ", threads, " threads."
  endif
#else
  print *, "OpenMP: not supported by this binary"
#endif
end subroutine

subroutine IO_backend_status
  print *, "Default I/O backend: ", backend_type
end subroutine

subroutine target_status
#ifndef LEVELS_TARGET
  print *, "LevelS compilation target: unknown"
#else
  print *, "LevelS compilation target: ", LEVELS_TARGET
#endif
end subroutine

subroutine announce (name, verbose)
  character(len=*), intent(in) :: name
  logical, intent(in), optional :: verbose
  integer nlen, i

  if (present(verbose).and. (.not. verbose)) return
  nlen=len(trim(name))
  print *
  print *,"+-",("-",i=1,nlen),"-+"
  print *,"| ",trim(name)," |"
  print *,"+-",("-",i=1,nlen),"-+"
  print *
  call openmp_status
  call IO_backend_status
  call target_status()
  print *
end subroutine

subroutine module_startup (name,args,verbose)
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: args(:)
  logical, intent(in) :: verbose

  if (verbose) call announce(name)
  if (size(args)<1) then
    if (verbose) then
      print *,"Usage:"
      print *,trim(name)," <parameter file / init object>"
      print *,"or:"
      print *,trim(name)," par1=val1 par2=val2 ..."
    endif
    call exit_with_status(1)
  endif
end subroutine

end module announce_mod
