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

module planck_missionio_new
  use planck_config
  use general_const
  use general_error
  use general_vector
  use general_time
  use ls_paramfile_io
  implicit none
  private

  public :: pl_readpointings_nom

  ! Angular position on the sphere in ecliptic coordinates with longitude
  ! lambda (degrees) and latitude beta (degrees).
  type gneclangle
    real(dp) :: lambda, beta
  end type gneclangle

contains

  ! Convert ecliptic position to angular position on the sphere.
  function gn_eclang2ang(eclang) result(ang)
    type(gneclangle), intent(in) :: eclang
    type(gnsphangle_double) :: ang

    ! The longitude, lambda, must be between 0.0 and 360.0.
    call gn_assert(eclang%lambda>=0.0, &
      'gn_testeclang: lambda < 0.0', eclang%lambda)
    call gn_assert(eclang%lambda<360.0, &
      'gn_testeclang: lambda >= 360.0', eclang%lambda)

    ! The latitude, beta, must be between - 90.0 and 90.0.
    call gn_assert(eclang%beta>=-90.0, &
      'gn_testeclang: beta < -90.0', eclang%beta)
    call gn_assert(eclang%beta<=90.0, &
      'gn_testeclang: beta > 90.0', eclang%beta)
    ang%theta = halfpi * (1.0 - eclang%beta / 90.0)
    ang%phi = eclang%lambda / GNRAD_DEG
  end function gn_eclang2ang

  ! Read a file of nominal pointings, allocating the arrays for the
  ! zaxis angles, pointing start times and pointing durations.
  subroutine pl_readpointings_nom(filename, dataformat, params, &
    n_pt, zaxis_pts, t_startpts, period_pts)
    character(len=*), intent(in) :: filename, dataformat
    type(paramfile_handle), intent(inout) :: params
    integer, intent(out) :: n_pt
    type(gnsphangle_double), pointer :: zaxis_pts(:)
    type(gnsec), pointer :: t_startpts(:)
    real(dp), pointer :: period_pts(:)

    integer :: n_ln, pt, p1, p2
    type(gneclangle) :: zaxis_ecl
    character(len=20) :: st1
    character(len=filenamelen) :: line

    ! PRE-PROCESSING

    ! Have a first look at the file to assess how many pointings are in it.
    call gn_assert((dataformat=='ppl') .or. (dataformat=='appls') &
      .or. (dataformat=='plan'), &
      'pl_readpointings_nom: dataformat unknown', dataformat)

    ! obtain number of lines in the file
    open(1,file=filename,status='old',action='read',err=4711)
    n_ln=0
    do
      read(1,'(A)',err=4711,end=1)
      n_ln = n_ln+1
    end do
 1  close(1)

    if (dataformat == 'ppl') then
      call gn_assert(n_ln>2,'pl_readpointings_nom: ppl format: n_ln <= 2')
      n_pt = n_ln-2
      print *,"found ",n_pt," pointings in the PPL file"
    else if (dataformat == 'plan') then
      call gn_assert(n_ln>2,'pl_readpointings_nom: plan format: n_ln <= 2')
      n_pt = n_ln-2
      print *,"found ",n_pt," pointings in the plan file"
    else if (dataformat == 'appls') then
      call gn_assert(n_ln>52, &
        'pl_readpointings_nom: appls format: n_ln <= 52')
      open(1,file=filename,status='old',action='read',err=4711)
      read(1,'(51/,a)',err=4711)

      n_pt = 0
      do
        read (1,'(A)',err=4711) line
        if (line=='') exit
        n_pt = n_pt + 1
      end do
      close(1)
      print *,"found ",n_pt," pointings in the APPLS file"
    end if

    ! select subset of pointings
    p1 = parse_int(params,'first_period',vmin=-1,default=-1)
    p2 = parse_int(params,'last_period',vmin=-1,default=-1)
    if (p1<=0) p1=1
    if (p2<=0) p2=n_pt
    call gn_assert (p2<=n_pt, "last_period too large", p2)
    call gn_assert (p1<=p2, "first_period larger than last_period", p1)
    n_pt = p2-p1+1

    ! Allocate pointings.
    allocate(zaxis_pts(n_pt),t_startpts(n_pt),period_pts(n_pt))

    open(1,file=filename,status='old',action='read',err=4711)

    ! Skip header
    if (dataformat == 'ppl') then
      read (1,'(/)',err=4711)
    else if (dataformat == 'plan') then
      read (1,'(/)',err=4711)
    else if (dataformat == 'appls') then
      read(1,'(51/)',err=4711)
    endif
    ! Skip unwanted leading pointing periods
    do pt = 1,p1-1
      read (1,*)
    end do

    ! BODY
    do pt=1,n_pt
      if (dataformat == 'ppl') then
        read (1,'(9X,F8.0,1X,F8.0,1X,A20,43X,F6.0)',err=4711) &
          zaxis_ecl%lambda, zaxis_ecl%beta, st1, period_pts(pt)
      else if (dataformat == 'plan') then
        read (1,'(F7.0,1X,F7.0,1X,A20,1X,F8.0)',err=4711) &
          zaxis_ecl%lambda, zaxis_ecl%beta, st1, period_pts(pt)
      else if (dataformat == 'appls') then
        read (1,'(A20,10X,F8.0,1X,F8.0,18X,F5.0)',err=4711) &
          st1, zaxis_ecl%lambda, zaxis_ecl%beta, period_pts(pt)
      endif

      zaxis_pts(pt) = gn_eclang2ang(zaxis_ecl)
      t_startpts(pt) = gn_string2sec(st1, 'ppl')
    end do

    close(1)

! Fix period durations in case of "plan" format
    if (dataformat == 'plan') then
      print *,"Warning: 'plan' format, stitching pointing periods together"
      do pt=1,n_pt-1
        period_pts(pt) = t_startpts(pt+1)-t_startpts(pt)
      end do
    endif

    return

4711 call exit_with_status(1,"error reading input file")

  end subroutine pl_readpointings_nom

end module planck_missionio_new
