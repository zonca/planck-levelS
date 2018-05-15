!-----------------------------------------------------------------------------
!
!  Copyright (C) 1999-2013 Daniel Mortlock
!
!  This file is part of the "simmission" component of the Planck simulation
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

module general_error
  use planck_config
  use ls_misc_utils
  implicit none
  private

  !FIXME everything warning-related is simmission3-only

  public :: gn_setwritewarnings, gn_warning, gn_assert, gn_fatal

  ! Print out a warning message with error value of specified type.
  interface gn_warning
    module procedure gn_warning_null, gn_warning_int, &
      gn_warning_string, gn_warning_double
  end interface

  ! Print out a fatal error message with error value of specified type, if
  ! a given condition does not hold
  interface gn_assert
    module procedure gn_assert_null, gn_assert_int, &
      gn_assert_string, gn_assert_double
  end interface

  ! Print out a fatal error message with error value of specified type.
  interface gn_fatal
    module procedure gn_fatal_null, gn_fatal_int, &
      gn_fatal_string, gn_fatal_double
  end interface

  ! These flags determine if low-grade warning messages and normal warning
  ! messages are written out to the screen. The normal mode of operation is
  ! to write out normal warnings (which shouldn't happen in the course of
  ! operations) and ignore messages (which can easily happen without anything
  ! going wrong).
  logical, save :: &
    gnwritemessages = .false., &
    gnwritewarnings = .true.

contains

  ! Set gnwritewarnings to the given value.
  subroutine gn_setwritewarnings(writewarnings)
    logical, intent(in) :: writewarnings
    gnwritewarnings = writewarnings
  end subroutine gn_setwritewarnings

  subroutine gn_diag (flag,label,text,vtext)
    logical, intent(in) :: flag
    character(len=*), intent(in) :: label, text
    character(len=*), intent(in), optional :: vtext

    if (flag) then
      if (present(vtext)) then
        write(*, '(/, a, a, a, a, a, a)') trim(label), ': ', trim(text), ': ', &
          trim(vtext), '.'
      else
        write(*, '(/, a, a, a, a)') trim(label), ': ', trim(text), '.'
      endif
    end if
  end subroutine gn_diag

  ! Print out warning message with no error value.
  subroutine gn_warning_null(message)
    character(len=*), intent(in) :: message
    call gn_diag (gnwritewarnings, 'Warning', message)
  end subroutine gn_warning_null

  ! Print out warning message with string error value.
  subroutine gn_warning_string(message, value)
    character(len=*), intent(in) :: message, value
    call gn_diag (gnwritewarnings, 'Warning', message, value)
  end subroutine gn_warning_string

  ! Print out warning message with integer error value.
  subroutine gn_warning_int(message, value)
    character(len=*), intent(in) :: message
    integer, intent(in) :: value
    call gn_diag (gnwritewarnings, 'Warning', message, toString(value))
  end subroutine gn_warning_int

  ! Print out warning message with real (double precision) error value.
  subroutine gn_warning_double(message, value)
    character(len=*), intent(in) :: message
    real(dp), intent(in) :: value
    call gn_diag (gnwritewarnings, 'Warning', message, toString(value))
  end subroutine gn_warning_double

  ! Print out fatal error message with no error value.
  subroutine gn_fatal_null(message)
    character(len=*), intent(in) :: message
    call gn_diag (.true., 'Fatal error', message)
    call exit_with_status(1)
  end subroutine gn_fatal_null

  ! Print out fatal error message with string error value.
  subroutine gn_fatal_string(message, value)
    character(len=*), intent(in) :: message, value
    call gn_diag (.true., 'Fatal error', message, value)
    call exit_with_status(1)
  end subroutine gn_fatal_string

  ! Print out fatal error message with integer error value.
  subroutine gn_fatal_int(message, value)
    character(len=*), intent(in) :: message
    integer, intent(in) :: value
    call gn_fatal (message, toString(value))
  end subroutine gn_fatal_int

  ! Print out fatal error message with real (double precision) error value.
  subroutine gn_fatal_double(message, value)
    character(len=*), intent(in) :: message
    real(dp), intent(in) :: value
    call gn_fatal (message, toString(value))
  end subroutine gn_fatal_double

  subroutine gn_assert_null(cond, message)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: message
    if (cond) return
    call gn_diag (.true., 'Assertion failed', message)
    call exit_with_status(1)
  end subroutine gn_assert_null

  subroutine gn_assert_string(cond, message, value)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: message, value
    if (cond) return
    call gn_diag (.true., 'Assertion failed', message, value)
    call exit_with_status(1)
  end subroutine gn_assert_string

  subroutine gn_assert_int(cond, message, value)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: message
    integer, intent(in) :: value
    if (cond) return
    call gn_assert (cond, message, toString(value))
  end subroutine gn_assert_int

  subroutine gn_assert_double(cond, message, value)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: message
    real(dp), intent(in) :: value
    if (cond) return
    call gn_assert (cond, message, toString(value))
  end subroutine gn_assert_double

end module general_error
