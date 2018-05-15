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
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "writers.h"
#include "arr.h"
#include "io_utils.h"

using namespace std;

detpt_writer::detpt_writer (bool full_info, bool singleprec,
  const string &filename, int first, int last)
  {
  string type = "pointing.LS_detpoint";
  if (full_info) type += "_with_orientation";
  if (!singleprec) type += "_dp";
  out=HandleManager.createObject(filename,type);
  out->setKey("firstPointingID",first);
  out->setKey("lastPointingID",last);
  ctheta = out->columnNumber("theta");
  cphi = out->columnNumber("phi");
  cpsi = full_info ? out->columnNumber("psi") : -1;
  }

void detpt_writer::write (const arr<pointing> &ptg)
  {
  planck_assert (cpsi<0, "detpt_writer: wrong writing function for this type");
  appendPointing (*out,ctheta,cphi,ptg);
  }

void detpt_writer::write (const arr<pointing> &ptg, const arr<double> &heading)
  {
  planck_assert (cpsi>=0,"detpt_writer: wrong writing function for this type");
  appendPointing (*out,ctheta,cphi,cpsi,ptg,heading);
  }

tod_writer::tod_writer (const string &filename, int first, int last)
  {
  out=HandleManager.createObject(filename,"toi.LS_toi");
  out->setKey("firstPointingID",first);
  out->setKey("lastPointingID",last);
  todcol = out->columnNumber("signal");
  }

void tod_writer::write (const arr<double> &toddata)
  {
  out->appendColumn(todcol,toddata);
  }

time_writer::time_writer (const string &filename)
  {
  out=HandleManager.createObject(filename,"toi.LS_timestamps");
  timecol = out->columnNumber("timestamp");
  }

void time_writer::write (const arr<double> &times)
  {
  arr<double> timedata(times.size());
  for (tsize m=0; m<times.size(); ++m) timedata[m]=1e6*times[m];
  out->appendColumn(timecol,timedata);
  }


quat_writer::quat_writer (const string &filename, const string &fpdb)
  {
  out=HandleManager.createObject(filename,"quat.LS_satpt_quat");
  qxcol = out->columnNumber("quat_x");
  qycol = out->columnNumber("quat_y");
  qzcol = out->columnNumber("quat_z");
  qwcol = out->columnNumber("quat_w");
  out->setKey("TUPLENbVec",int16(4));
  out->setKey("RefFrame",string("ECLIPTIC"));
  out->setKey("FormatPtg",string("QUATERNION"));
  out->setKey("ID",string("SM"));
  out->setKey("IMOID",fpdb);
  out->setKey("IDBolo",string(""));
  out->setKey("IDBeam",string(""));
  out->setKey("FormatOri",string(""));
  out->setKey("ConvOri",string(""));
  }

void quat_writer::write (const arr<quaternion> &quatdata)
  {
  appendQuaternions (*out,qxcol,qycol,qzcol,qwcol,quatdata);
  }


flag_writer::flag_writer (const string &filename)
  {
  out=HandleManager.createObject(filename,"flag.LS_flag");
  flagcol = out->columnNumber("flag");
  }

void flag_writer::write (const arr<bool> &flagdata)
  {
  out->appendColumn(flagcol,flagdata);
  }
