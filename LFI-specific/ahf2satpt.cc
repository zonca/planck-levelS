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

#include <fstream>
#include "paramfile.h"
#include "iohandle_current.h"
#include "io_utils.h"
#include "announce.h"
#include "rotmatrix.h"
#include "lsconstants.h"

using namespace std;

namespace {

void parse_words (const string &filename, vector<string> &words)
  {
  words.clear();
  ifstream inp (filename.c_str());
  string dummy;
  while (!inp.eof())
    {
    dummy="";
    inp >> dummy;
    if (!dummy.empty()) words.push_back(dummy);
    }
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
  PLANCK_DIAGNOSIS_BEGIN
  module_startup ("ahf2satpt", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  safe_ptr<iohandle> out (HandleManager.createObject
    (params.find<string>("outfile"),"sat.LS_satpoint_real"));
  string wobble=params.find<string>("outfile_wobble");
  int cxo = out->columnNumber("quat_x");
  int cyo = out->columnNumber("quat_y");
  int czo = out->columnNumber("quat_z");
  int cwo = out->columnNumber("quat_w");

  vector<string> infiles;
  parse_words (params.find<string>("input_list"),infiles);

  vector<int64> pp_ifirst;
  vector<double> pp_tstart, pp_tstable, pp_tilt1, pp_tilt2;
  vector<string> pp_name;
  int64 glob_offset=0;
  double tlast=-1;
  tdiff nval=-1;
  double t1acc=0, t2acc=0;
  int period=params.find<int>("first_period_number",1)-1;

  for (tsize m=0; m<infiles.size(); ++m)
    {
    cout << "processing AHF " << infiles[m]
         << ", first period: " << period+1 << endl;
    safe_ptr<iohandle> inp (HandleManager.openObject
      (infiles[m],"toi.attitude.HighFrequency"));

    arr<double> times, tilt1, tilt2;
    arr<string> name;
    arr<quaternion> quat;
    int cx = inp->columnNumber("quaternionX");
    int cy = inp->columnNumber("quaternionY");
    int cz = inp->columnNumber("quaternionZ");
    int cw = inp->columnNumber("quaternionS");

    inp->readEntireColumn("timeDataValid",times);
// CAUTION! The tilt angles are stored in arcmin, even if the DDL says degrees!
    inp->readEntireColumn("tiltAngle1",tilt1);
    inp->readEntireColumn("tiltAngle2",tilt2);
    inp->readEntireColumn("pointingID",name);
    tsize nsamp = times.size();
    planck_assert(name.size()==nsamp,"column size mismatch");

    const double tfact = ldexp(1.,-16);
    for (tsize i=0; i<nsamp; ++i)
      times[i]*=tfact;
    tlast = times[times.size()-1];
    quat.alloc(nsamp);
    readQuaternions(*inp, cx, cy, cz, cw, quat, 0);
    planck_assert(quat.size()==nsamp,"column size mismatch");

    out->appendColumn("quat_time",times);
    string cur_name="foobar";

    for (tsize i=0; i<nsamp; ++i)
      {
      if (name[i]!=cur_name)
        {
        cur_name=name[i];
        ++period;
//        cout << "pointing period " << period << " starting at " << i << endl;
        pp_ifirst.push_back(i+glob_offset);
        pp_tstart.push_back(times[i]);
        pp_name.push_back(dataToString(period));
        if (nval>0)
          {
          pp_tilt1.push_back(t1acc/nval);
          pp_tilt2.push_back(t2acc/nval);
          }
        t1acc=t2acc=0;
        nval=0;
        }
      t1acc+=tilt1[i];
      t2acc+=tilt2[i];
      ++nval;
      }
    appendQuaternions(*out, cxo, cyo, czo, cwo, quat);
    glob_offset+=nsamp;
    }
  pp_tilt1.push_back(t1acc/nval);
  pp_tilt2.push_back(t2acc/nval);

  tsize npp=pp_name.size();
  planck_assert(pp_tilt1.size()==npp,"inconsistency in tilt angles");
  arr<int32> pp_nquat(npp);
  arr<double> pp_tend(npp);
  for (tsize i=0; i<npp-1; ++i)
    {
    pp_tend[i] = pp_tstart[i+1];
    pp_nquat[i] = pp_ifirst[i+1]-pp_ifirst[i];
    }
  pp_tend[npp-1]=tlast;
  pp_nquat[npp-1]=glob_offset-pp_ifirst[npp-1];
  out->appendColumn("pp_name",pp_name);
  out->appendColumn("pp_tstart",pp_tstart);
  out->appendColumn("pp_tend",pp_tend);
  out->appendColumn("pp_ifirst",pp_ifirst);
  out->appendColumn("pp_nquat",pp_nquat);

  if (wobble!="")
    {
    safe_ptr<iohandle> outw (HandleManager.createObject
      (wobble,"sat.LS_tiltAngles"));
    outw->appendColumn("tilt1",pp_tilt1);
    outw->appendColumn("tilt2",pp_tilt2);
    }
  PLANCK_DIAGNOSIS_END
  }
