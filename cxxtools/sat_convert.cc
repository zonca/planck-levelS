/*
 *  This file is part of the Planck simulation package.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  This code is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2009-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "iohandle_current.h"
#include "io_utils.h"
#include "pointing.h"
#include "arr.h"
#include "paramfile.h"

using namespace std;

namespace {

void read_entire_timeinfo (iohandle &inp, arr<double> &tstart,
  arr<double> &tend)
  {
  double toffset = sec_58_70;
  if (inp.keyPresent("TREF1958") && inp.getKey<bool>("TREF1958")) toffset=0.;
  inp.readEntireColumn("t_endpt_int",tend);
  arr<double> tmp;
  inp.readEntireColumn("t_endpt_frac",tmp);
  for (tsize i=0; i<tmp.size(); ++i) tend[i]+=tmp[i]+toffset;
  inp.readEntireColumn("t_startpt_int",tstart);
  inp.readEntireColumn("t_startpt_frac",tmp);
  for (tsize i=0; i<tmp.size(); ++i) tstart[i]+=tmp[i]+toffset;

  for (tsize i=0; i<tstart.size()-1; ++i)
    planck_assert (approx(tend[i],tstart[i+1]),
      "gap between pointing periods detected!");
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("sat_convert",argc,argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  safe_ptr<iohandle> inp (HandleManager.openObject
    (params.find<string>("infile"),"satelliteinfo.LS_satinfo"));
  safe_ptr<iohandle> out (HandleManager.createObject
    (params.find<string>("outfile"),"sat.LS_satpoint_real"));

  arr<double> tstart, tend;
  read_entire_timeinfo (*inp, tstart, tend);
  arr<int> nsamples;
  inp->readEntireColumn("nsamples",nsamples);

  tsize npp=tstart.size();

  arr<string> ppname(npp);
  for (tsize m=0; m<npp; ++m)
    ppname[m]=dataToString(m+1);
  arr<int64> ifirst(npp);
  ifirst[0]=0;
  for (tsize m=1; m<npp; ++m)
    ifirst[m]=ifirst[m-1]+nsamples[m-1];

  int thetaxcol = inp->columnNumber("theta_x"),
      phixcol = inp->columnNumber("phi_x"),
      thetaycol = inp->columnNumber("theta_y"),
      phiycol = inp->columnNumber("phi_y"),
      thetazcol = inp->columnNumber("theta_z"),
      phizcol = inp->columnNumber("phi_z");
  int qwcol = out->columnNumber("quat_w"),
      qxcol = out->columnNumber("quat_x"),
      qycol = out->columnNumber("quat_y"),
      qzcol = out->columnNumber("quat_z"),
      qtcol = out->columnNumber("quat_time"),
      qfcol = out->columnNumber("quat_flag");

  for (tsize period=0; period<npp; ++period)
    {
    arr<pointing> real_px(nsamples[period]),
                  real_py(nsamples[period]),
                  real_pz(nsamples[period]);
    readPointing (*inp,thetaxcol,phixcol,thetaycol,phiycol,thetazcol,phizcol,
                  real_px,real_py,real_pz,ifirst[period]);
    double dt = (tend[period]-tstart[period])/nsamples[period];

    arr<quaternion> rotquat (nsamples[period]);
    arr<int32> qflag (nsamples[period]);
    arr<float64> qtime (nsamples[period]);
    for (tsize m=0; m<tsize(nsamples[period]); ++m)
      {
      rotquat[m]=quaternion(rotmatrix(real_px[m].to_vec3(),real_py[m].to_vec3(),
                                      real_pz[m].to_vec3()));
      qflag[m] = 0;
      qtime[m] = tstart[period]+(m+0.5)*dt;
      }
    appendQuaternions (*out, qxcol, qycol, qzcol, qwcol, rotquat);
    out->appendColumn(qtcol,qtime);
    out->appendColumn(qfcol,qflag);
    }

  out->appendColumn("pp_name",ppname);
  out->appendColumn("pp_tstart",tstart);
  out->appendColumn("pp_tend",tend);
  out->appendColumn("pp_ifirst",ifirst);
  out->appendColumn("pp_nquat",nsamples);
PLANCK_DIAGNOSIS_END
  }
