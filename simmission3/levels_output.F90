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

module levels_output
use planck_config
use dmc_io
use ls_misc_utils
use ls_paramfile_io
use general_const
use general_time
use general_vector
implicit none
private

type(dmc_handle), save :: mhandle
integer, save :: satcol(23)

real(dp), allocatable, save :: summary(:,:)
integer(i4b), allocatable, save :: nsamples(:)

integer, save :: samples_per_period
integer(i8b) :: eff_size

public lsoutput_init, lsoutput_finish, add_summary, add_pointing

contains

subroutine lsoutput_init (params,nptg)
  type(paramfile_handle), intent(inout) :: params
  integer, intent(in) :: nptg

  character (len=filenamelen) :: mfile

  mfile = parse_string(params,'file_mission')

  allocate(summary(17,nptg))
  allocate(nsamples(nptg))
  nsamples = 0

  call dmc_create(mhandle,mfile,"satelliteinfo.LS_satinfo")
  call dmc_set_key (mhandle, "sataxis_convention", "HFI")
  call dmc_set_key (mhandle, "TREF1958", .true.)
  call dmc_reserve_column (mhandle, "t_startpt_int", int(nptg,i8b))
  satcol = (/ dmc_colnum(mhandle,"t_startpt_int"), &
              dmc_colnum(mhandle,"t_startpt_frac"), &
              dmc_colnum(mhandle,"satpos_x_startpt"), &
              dmc_colnum(mhandle,"satpos_y_startpt"), &
              dmc_colnum(mhandle,"satpos_z_startpt"), &
              dmc_colnum(mhandle,"theta_x_startpt"), &
              dmc_colnum(mhandle,"phi_x_startpt"), &
              dmc_colnum(mhandle,"theta_y_startpt"), &
              dmc_colnum(mhandle,"phi_y_startpt"), &
              dmc_colnum(mhandle,"theta_z_startpt"), &
              dmc_colnum(mhandle,"phi_z_startpt"), &
              dmc_colnum(mhandle,"rate_rot_z_startpt"), &
              dmc_colnum(mhandle,"t_endpt_int"), &
              dmc_colnum(mhandle,"t_endpt_frac"), &
              dmc_colnum(mhandle,"satpos_x_endpt"), &
              dmc_colnum(mhandle,"satpos_y_endpt"), &
              dmc_colnum(mhandle,"satpos_z_endpt"), &
              dmc_colnum(mhandle,"theta_x"), &
              dmc_colnum(mhandle,"phi_x"), &
              dmc_colnum(mhandle,"theta_y"), &
              dmc_colnum(mhandle,"phi_y"), &
              dmc_colnum(mhandle,"theta_z"), &
              dmc_colnum(mhandle,"phi_z") /)

  samples_per_period = -1
  eff_size = dmc_eff_chunksize(mhandle,"theta_x")
end subroutine lsoutput_init

subroutine add_summary (pt, t_startpt, axes_startpt, rate_z_rot_startpt,t_endpt)
  integer, intent(in) :: pt
  type(gnsec), intent(in) :: t_startpt
  type(gnsphangle_double), intent(in) :: axes_startpt(3)
  real(dp), intent(in) :: rate_z_rot_startpt
  type(gnsec), intent(in) :: t_endpt

  real(dp),parameter :: silly_number=1e300_dp

  summary(:,pt) = (/ real(t_startpt%int,dp)+sec_58_70, &
                     t_startpt%frac, &
                     -silly_number, &
                     -silly_number, &
                     -silly_number, &
                     axes_startpt(3)%theta, &
                     axes_startpt(3)%phi, &
                     pi-axes_startpt(2)%theta, &
                     mod(axes_startpt(2)%phi+pi,twopi), &
                     axes_startpt(1)%theta, &
                     axes_startpt(1)%phi, &
                     rate_z_rot_startpt, &
                     real(t_endpt%int,dp)+sec_58_70, &
                     t_endpt%frac, &
                     silly_number, &
                     silly_number, &
                     silly_number /)

end subroutine add_summary

subroutine add_pointing(pt, sataxes)
  integer, intent(in) :: pt
  real(dp), pointer :: sataxes(:, :)

  integer(i8b) nchunks, i, i1, i2

  nsamples(pt) = size(sataxes,dim=1)
  nchunks = nsamples(pt)/eff_size
  if (eff_size*nchunks < nsamples(pt)) nchunks=nchunks+1
  do i=0,nchunks-1
    i1 = i*eff_size+1
    i2 = min((i+1)*eff_size,nsamples(pt))
    call dmc_append_column (mhandle,satcol(18),sataxes(i1:i2,5))
    call dmc_append_column (mhandle,satcol(19),sataxes(i1:i2,6))
    call dmc_append_column (mhandle,satcol(20),pi-sataxes(i1:i2,3))
    call dmc_append_column (mhandle,satcol(21),mod(sataxes(i1:i2,4)+pi,twopi))
    call dmc_append_column (mhandle,satcol(22),sataxes(i1:i2,1))
    call dmc_append_column (mhandle,satcol(23),sataxes(i1:i2,2))
  end do
end subroutine add_pointing

subroutine lsoutput_finish
  integer col

  do col=1,17
    call dmc_append_column (mhandle,satcol(col),summary(col,:))
  end do
  call dmc_append_column (mhandle,dmc_colnum(mhandle,"nsamples"),nsamples)
  deallocate(summary)
  deallocate(nsamples)
  call dmc_close(mhandle)

  call dmc_shutdown
end subroutine lsoutput_finish

end module levels_output
