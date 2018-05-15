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

module general_time
  use planck_config
  use general_const
  use general_error
  implicit none
  private

  public :: gnsec, gn_string2sec, operator(+), operator(>)

  !FIXME simmission3 only
  public :: gn_strings2sec, gn_sec2cent, gn_cent2sec, gn_sec2year, &
    operator(-), operator(>=)

  ! The date, given by the three integers yr, mth and day.
  type gndate
    integer :: yr
    integer :: mth
    integer :: day
  end type gndate

  ! The time (of day), given by the two integers hr, min, and the
  ! real sec.
  type gntime
    integer :: hr
    integer :: min
    real(dp) :: sec
  end type gntime

  ! The time as specified by date and time (of day).
  type gndateandtime
    type(gndate) :: date
    type(gntime) :: time
  end type gndateandtime

  type gnsec
    integer(i8b) :: int
    real(dp) :: frac
  end type gnsec

  ! Addition interface for gnsec.
  interface operator(+)
    module procedure gn_addsec_double_post
  end interface

  ! Subtraction interface for gnsec.
  interface operator(-)
    module procedure gn_subtractsec, gn_subtractsec_double
  end interface

  ! Greater than or equal to relational operator for gnsec.
  interface operator(>=)
    module procedure gn_gesec
  end interface

  ! Greater than relational operator for gnsec.
  interface operator(>)
    module procedure gn_gtsec
  end interface

  ! Time conversions between centuries (CENT), years (YR), days (DAY),
  ! hours (HR), minutes (MIN), seconds (SEC) and milliseconds (MSEC).
  real(dp), public, parameter :: &
    GNYR_SEC = 31557600.0_dp, &
    GNDAY_SEC = 86400.0_dp, &
    GNHR_SEC = 3600.0_dp, &
    GNMIN_SEC = 60.0_dp

  ! Details of the conversion from calendar/Julian dates to ``Java'' time,
  ! which is defined as the integer number of milliseconds since 00:00:00
  ! on January 1, 1970 (GMT). Due to the ``resolution'' this must be
  ! measured in long integers (as defined in gn_types). GNJULIAN_1970 is
  ! the Julian date corresponding to this zero for time in milliseconds
  real(dp), parameter :: &
    GNJULIAN_1970 = 2440587.5_dp

  real(dp), parameter :: &
    GNCENT_DAY = 36525.0_dp, &
    GNYR_DAY = 365.25_dp

contains

  ! The number of days in the month of the calendar year specified. A fatal
  ! error occurs if mth < 1 or mth > 12.
  function gn_maxday_yrmth(yr, mth) result(nday)
    integer, intent(in) :: yr, mth
    integer :: nday

    call gn_assert(yr>=1950,'gn_maxday_yrmth: yr < 1950')
    call gn_assert(mth>=1,'gn_maxday_yrmth: mth < 1')
    call gn_assert(mth<=12,'gn_maxday_yrmth: mth > 12')

    if ((mth == 1) .or. (mth == 3) .or. (mth == 5) .or. (mth == 7) &
      .or. (mth == 8) .or. (mth == 10) .or. (mth == 12)) then
      nday = 31
    else if ((mth == 4) .or. (mth == 6) .or. (mth == 9) .or. (mth == 11)) then
      nday = 30
    else
      ! Here the various leap years must be dealt with. Prior to the
      ! Gregorian reform every year divisible by 4 was a leap year. In
      ! the time since the following exception has been introduced: every
      ! year that is divisible by 100 but not by 400 is *not* a leap year.
      if (mod(yr, 400)==0) then
        nday = 29
      else if (mod(yr, 100)==0) then
        nday = 28
      else if (mod(yr, 4)==0) then
        nday = 29
      else
        nday = 28
      end if
    end if
  end function gn_maxday_yrmth

  ! Test date and time to make sure it is valid.
  subroutine gn_testdateandtime(dateandtime)
    type(gndateandtime), intent(in) :: dateandtime

    ! Test date and time separately using dedicated routines.
    call gn_testdate(dateandtime%date)
    call gn_testtime(dateandtime%time)
  end subroutine gn_testdateandtime

  ! Test date to make sure it is valid.
  subroutine gn_testdate(date)
    type(gndate), intent(in) :: date

    call gn_assert(date%yr>=1950,'gn_testdate: yr < 1950')
    call gn_assert(date%mth>=1,'gn_testdate: mth < 1',date%mth)
    call gn_assert(date%mth<=12,'gn_testdate: mth > 12',date%mth)
    call gn_assert(date%day>=1,'gn_testdate: day < 1',date%day)
    call gn_assert(date%day<=gn_maxday_yrmth(date%yr,date%mth), &
      'gn_testdate: day > maxday_yrmth',date%day)
  end subroutine gn_testdate

  ! Test time (of day) to make sure it is valid.
  subroutine gn_testtime(time)
    type(gntime), intent(in) :: time

    call gn_assert(time%hr>=0,'gn_testtime: hr < 0',time%hr)
    call gn_assert(time%hr<=23,'gn_testtime: hr > 23',time%hr)
    call gn_assert(time%min>=0,'gn_testtime: min < 0',time%min)
    call gn_assert(time%min<=59,'gn_testtime: min > 59',time%min)
    call gn_assert(time%sec>=0,'gn_testtime: sec < 0.0',time%sec)
    call gn_assert(time%sec<=60.0,'gn_testtime: sec >= 60.0',time%sec)
  end subroutine gn_testtime

  ! Convert a string in the given format to time in seconds since 1970.
  function gn_string2sec(string, format) result(sec)
    character(len=*), intent(in) :: string, format
    type(gnsec) :: sec

    sec = gn_dateandtime2sec(gn_string2dateandtime(string,format))
  end function gn_string2sec

  ! Convert a string in the given format to time in seconds since 1970.
  function gn_strings2sec(string_date, format_date, string_time, &
    format_time) result(sec)
    character(len=*), intent(in) :: string_date, format_date
    character(len=*), intent(in) :: string_time, format_time
    type(gnsec) :: sec

    sec = gn_dateandtime2sec(gn_strings2dateandtime &
      (string_date, format_date, string_time, format_time))
  end function gn_strings2sec

  ! Convert a string in the given format to date and time.
  function gn_string2dateandtime(string, format) result(dateandtime)
    character(len=*), intent(in) :: string, format
    type(gndateandtime) :: dateandtime

    call gn_assert (format=='ppl', &
      'gn_string2dateandtime: format unknown', format)

    dateandtime = gn_strings2dateandtime(string(1:10),'ppl',string(12:19),'ppl')
  end function gn_string2dateandtime

  ! Convert date and time strings in the given formats to date and time.
  function gn_strings2dateandtime(string_date, format_date, string_time, &
    format_time) result(dateandtime)
    character(len=*), intent(in) :: string_date, format_date
    character(len=*), intent(in) :: string_time, format_time
    type(gndateandtime) :: dateandtime

    dateandtime%date = gn_string2date(string_date, format_date)
    dateandtime%time = gn_string2time(string_time, format_time)
  end function gn_strings2dateandtime

  ! Convert a date string in the given format to a date.
  function gn_string2date(string, format) result(date)
    character(len=*), intent(in) :: string, format
    type(gndate) :: date

    integer stat

    ! Split up according to the supported formats (or generate an error).
    if ((format == 'yyyymmdd') .or. (format == 'YYYYMMDD')) then
      read (string,'(I4,I2,I2)', iostat=stat) date%yr, date%mth, date%day
      call gn_assert(stat==0, 'gn_string2date: parse error')
    else if (format == 'ppl') then
      read (string,'(I4,1X,I2,1X,I2)', iostat=stat) date%yr, date%mth, date%day
      call gn_assert(stat==0, 'gn_string2date: parse error')
    else
      call gn_fatal('gn_string2date: format unknown', format)
    end if

    call gn_testdate(date)
  end function gn_string2date

  ! Convert a time string in the given format to a time.
  function gn_string2time(string, format) result (time)
    character(len=*), intent(in) :: string, format
    type(gntime) :: time

    integer stat

    ! Split up according to time format.
    if ((format == 'hhmmssdsss') .or. (format == 'HHMMSSDSSS')) then
      read (string,'(I2,I2,F6.0)', iostat=stat) time%hr, time%min, time%sec
      call gn_assert(stat==0, 'gn_string2time: parse error')
    else if (format == 'ppl') then
      read (string,'(I2,1X,I2,1X,F2.0)',iostat=stat) time%hr, time%min, time%sec
      call gn_assert(stat==0, 'gn_string2time: parse error')
    else
      call gn_fatal('gn_string2time: format unknown', format)
    end if

    call gn_testtime(time)
  end function gn_string2time

  ! Convert centuries to time in seconds since 1970.
  function gn_cent2sec(cent) result(sec)
    real(dp), intent(in) :: cent
    type(gnsec) :: sec

    real(dp) :: julian, t_sec

    julian = 1721045.0_dp + GNCENT_DAY * cent
    t_sec = GNDAY_SEC * (julian - GNJULIAN_1970)

    ! Break this into two parts.
    sec%int = nint(t_sec, i8b)
    sec%frac = t_sec - real(sec%int, dp)

    ! Correct the format to ensure 0.0 <= frac < 1.0.
    call gn_correctsec(sec)
  end function gn_cent2sec

  ! Convert time in seconds since 1970 to centuries.
  function gn_sec2cent(sec) result(cent)
    type(gnsec), intent(in) :: sec
    real(dp) :: cent

    real(dp) :: julian

    call gn_testsec(sec)
    julian = (real(sec%int, dp) + sec%frac) / GNDAY_SEC + GNJULIAN_1970
    cent = (julian - 1721045.0_dp) / GNCENT_DAY
  end function gn_sec2cent

  ! Convert time in seconds since 1970 to years.
  function gn_sec2year(sec) result(year)
    type(gnsec), intent(in) :: sec
    real(dp) :: year

    real(dp) :: julian

    call gn_testsec(sec)
    julian = (real(sec%int, dp) + sec%frac) / GNDAY_SEC + GNJULIAN_1970
    year = (julian - 1721045.0_dp) / GNYR_DAY
  end function gn_sec2year

  ! Convert date and time to time in seconds since 1970.
  function gn_dateandtime2sec(dateandtime) result(sec)
    type(gndateandtime), intent(in) :: dateandtime
    type(gnsec) :: sec

    call gn_testdateandtime(dateandtime)
    sec = gn_date2sec(dateandtime%date) + gn_time2sec(dateandtime%time)
  end function gn_dateandtime2sec

  ! Convert date to seconds since Jan 1, 1970 (at the start of that date).
  function gn_date2sec(date) result(sec)
    type(gndate), intent(in) :: date
    type(gnsec) :: sec

    integer :: a, mth, yr, julian_day
    real(dp) :: julian
    real(dp), parameter :: GNYR_MTH = 12.0_dp, &
                           GNMTH_DAY = 30.4375_dp

    call gn_testdate(date)

    ! First calculate the Julian date
    ! This is taken from the Press et al. (1992) algorithm.
    yr = date%yr
    if (yr < 0) yr = yr + 1

    if (date%mth > 2) then
      mth = date%mth + 1
    else
      yr = yr - 1
      mth = date%mth + 13
    end if

    julian_day = int(GNYR_DAY * yr) + int(30.6001 * mth) + date%day + 1720995

    if (date%day + floor(GNMTH_DAY + 1) * (date%mth + nint(GNYR_MTH) &
      * date%yr) >= 588829) then
      a = int(0.01 * yr)
      julian_day = julian_day + 2 - a + int(0.25 * a)
    end if

    julian = real(julian_day,dp)

    ! Finally, because of the definition here that the date is equivalent
    ! to the time at the start of that date, we need to subtract 0.5 days
    ! off this figure.
    julian = julian - 0.5

    sec%int = nint(julian-GNJULIAN_1970,kind=i8b)*nint(GNDAY_SEC,kind=i8b)
    sec%frac = 0
  end function gn_date2sec

  ! Convert time (of day) to seconds since midnight.
  function gn_time2sec(time) result(secs)
    type(gntime), intent(in) :: time
    real(dp) :: secs

    call gn_testtime(time)
    secs = real(time%hr,dp) * GNHR_SEC &
      + real(time%min,dp) * GNMIN_SEC &
      + time%sec
  end function gn_time2sec

  ! Test time of type gnsec to see if it is in correct format. The integer
  ! part can take on any value; the fractional part must be between 0.0 and
  ! 1.0. Note that rounding error can result in the fractional part being
  ! just below 0.0 or just above 1.0 very easily.
  subroutine gn_testsec(sec)
    type(gnsec), intent(in) :: sec

    call gn_assert(sec%frac>=0,'gn_testsec: frac < 0.0',sec%frac)
    call gn_assert(sec%frac<1.0,'gn_testsec: frac >= 1.0',sec%frac)
  end subroutine gn_testsec

  ! Correct time of type gnsec to convert to proper format by fixing
  ! rounding errors to ensure that frac is between 0.0 and 1.0.
  subroutine gn_correctsec(sec)
    type(gnsec), intent(in out) :: sec

    integer :: int_frac

    ! Calculate the integer part of the fractional component of sec.
    ! Note that floor only accepts a second ``i8b''-type of
    ! argument in Fortran95.
    int_frac = floor(sec%frac)
    sec%int = sec%int + int(int_frac, i8b)
    sec%frac = sec%frac - real(int_frac, dp)
  end subroutine gn_correctsec

  ! Add a double precision real time interval to a time of type gnsec.
  function gn_addsec_double_post(sec, interval) result(seci)
    type(gnsec), intent(in) :: sec
    real(dp), intent(in) :: interval
    type(gnsec) :: seci

    integer(i8b) :: interval_int

    ! Add the integer and fractional parts of the time separately.
    interval_int = nint(interval, i8b)
    seci%int = sec%int + interval_int
    seci%frac = sec%frac + (interval - real(interval_int, dp))

    ! Correct the format to ensure 0.0 <= frac < 1.0.
    call gn_correctsec(seci)
  end function gn_addsec_double_post

  ! Subtract two times of type gnsec (assuming that they are both
  ! correctly formatted to begin with).
  function gn_subtractsec(sec_1, sec_2) result(interval)
    type(gnsec), intent(in) :: sec_1, sec_2
    real(dp) :: interval

    interval = real(sec_1%int - sec_2%int, dp) + sec_1%frac - sec_2%frac
  end function gn_subtractsec

  ! Subtract a double precision real time interval from a time of type gnsec.
  function gn_subtractsec_double(sec, n) result(secn)
    type(gnsec), intent(in) :: sec
    real(dp), intent(in) :: n
    type(gnsec) :: secn

    ! Use the code for addition.
    secn = gn_addsec_double_post(sec, -n)
  end function gn_subtractsec_double

  ! Greater than or equal to relational operator for time of type gnsec.
  function gn_gesec(sec_1, sec_2) result(ge)
    type(gnsec), intent(in) :: sec_1, sec_2
    logical :: ge

    if (sec_1%int > sec_2%int) then
      ge = .true.
    else if (sec_1%int == sec_2%int) then
      ge = (sec_1%frac >= sec_2%frac)
    else
      ge = .false.
    end if
  end function gn_gesec

  ! Greater than relational operator for time of type gnsec.
  function gn_gtsec(sec_1, sec_2) result(gt)
    type(gnsec), intent(in) :: sec_1, sec_2
    logical :: gt

    if (sec_1%int > sec_2%int) then
      gt = .true.
    else if (sec_1%int == sec_2%int) then
      gt = (sec_1%frac>sec_2%frac)
    else
      gt = .false.
    end if
  end function gn_gtsec

end module general_time
