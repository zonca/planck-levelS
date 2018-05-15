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
 *  Copyright (C) 2010-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

/* Utility for converting Horizon text files to ephemeris DDL objects.
   Send the batch file to horizons@ssd.jpl.nasa.gov with subject "JOB".
   The required Horizon settings are the following:

!$$SOF
COMMAND= '0' '10' '199' '299' '399' '499' '599' '699' '799' '899' '999'
CENTER= '500@-489'
MAKE_EPHEM= 'YES'
TABLE_TYPE= 'OBSERVER'
START_TIME= '2009-05-15'
STOP_TIME= '2013-10-01'
STEP_SIZE= '1 h'
CAL_FORMAT= 'CAL'
TIME_DIGITS= 'FRACSEC'
ANG_FORMAT= 'DEG'
OUT_UNITS= 'KM-S'
RANGE_UNITS= 'KM'
APPARENT= 'AIRLESS'
SOLAR_ELONG= '0,180'
SUPPRESS_RANGE_RATE= 'NO'
SKIP_DAYLT= 'NO'
EXTRA_PREC= 'NO'
R_T_S_ONLY= 'NO'
REF_SYSTEM= 'J2000'
CSV_FORMAT= 'NO'
OBJ_DATA= 'YES'
QUANTITIES= '10,13,19,20,23,24,31'
!$$EOF
*/

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "iohandle_current.h"
#include "pointing.h"
#include "time_utils.h"
#include "lsconstants.h"
#include "announce.h"

using namespace std;

namespace {

void assert_string_present (const vector<string> &lines,
  const string &key)
  {
  for (tsize i=0; i<lines.size(); ++i)
    {
    string::size_type pos=lines[i].find(key);
    if (pos!=string::npos) return;
    }
  planck_fail ("could not find string '"+key+"' in HORIZON header");
  }

template <typename T> T get_value_after (const vector<string> &lines,
  const string &key)
  {
  for (tsize i=0; i<lines.size(); ++i)
    {
    string::size_type pos=lines[i].find(key);
    if (pos!=string::npos)
      {
      istringstream tmp (lines[i].substr(pos+key.length()));
      T result;
      tmp >> result;
      return result;
      }
    }
  planck_fail ("could not find value after '"+key+"'");
  }

void skip (istream &inp, int nwords)
  {
  string dum;
  for (int i=0; i<nwords; ++i) inp >> dum;
  }

void add_planet (iohandle &out, const string &horifile, const string &name,
  int64 &nsamples)
  {
  ifstream inp (horifile.c_str());
  vector<string> hlines;
  string line;
  planck_assert(getline(inp,line),"read error");
  while (line != "$$SOE")
    {
    hlines.push_back(line);
    planck_assert(getline(inp,line),"read error");
    }
  assert_string_present (hlines, "Target body name: "+name);
  assert_string_present (hlines, "Center body name: Planck Space Observatory");
  assert_string_present (hlines, "Center-site name: BODYCENTRIC");
  assert_string_present (hlines,
    "Start time      : A.D. 2009-May-15 00:00:00.0000 UT");
  assert_string_present (hlines,
    "Stop  time      : A.D. 2013-Oct-01 00:00:00.0000 UT");
  assert_string_present (hlines, "Step-size       : 60 minutes");
  assert_string_present (hlines,
    "Date__(UT)__HR:MN:SC.fff       Illu% Ang-diam               r        "
    "rdot            delta      deldot    S-O-T    S-T-O      "
    "phi  PAB-LON  PAB-LAT    ObsEcLon    ObsEcLat");
  assert_string_present (hlines, "RA format       : DEG");
  assert_string_present (hlines, "Time format     : CAL");
  out.appendColumn("bodies",name);
  int dcol = out.columnNumber("data");
  vector<double> xsat,ysat,zsat,vdsun,vangdiam,vsto;
  int64 nsamp=0;
  planck_assert(getline(inp,line),"read error");
  while (line != "$$EOE")
    {
    double dsun, dsat, angdiam, sto;
    ++nsamp;
    istringstream tmp (line);
    skip(tmp,3);
    tmp >> angdiam;
    vangdiam.push_back(angdiam);
    tmp >> dsun;
    dsun *= 1000; // km -> m
    vdsun.push_back(dsun);
    skip(tmp,1);
    tmp >> dsat;
    dsat *= 1000; // km -> m
    skip(tmp,2);
    tmp >> sto;
    vsto.push_back(sto);
    pointing ptg;
    skip(tmp,3);
    tmp >> ptg.phi >> ptg.theta;
    ptg.phi*=degr2rad;
    ptg.theta = halfpi - ptg.theta*degr2rad;
    vec3 relpos = ptg;
    relpos *= dsat;
    xsat.push_back(relpos.x);
    ysat.push_back(relpos.y);
    zsat.push_back(relpos.z);
    planck_assert(getline(inp,line),"read error");
    }
  if (nsamples==-1) nsamples=nsamp;
  planck_assert(nsamples==nsamp,"mismatch of nsamples");
  out.appendColumnRaw(dcol,&xsat[0],nsamp);
  out.appendColumnRaw(dcol,&ysat[0],nsamp);
  out.appendColumnRaw(dcol,&zsat[0],nsamp);
  out.appendColumnRaw(dcol,&vangdiam[0],nsamp);
  out.appendColumnRaw(dcol,&vdsun[0],nsamp);
  out.appendColumnRaw(dcol,&vsto[0],nsamp);
  }

} // unnamed namespace

int main(int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("read_horizons", argc, argv);

  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  safe_ptr<iohandle> out (HandleManager.createObject
    (params.find<string>("outfile"),"ephemeris.LS_ephemeris"));

//  out->appendColumn("scalars",string("mass/kg"));
  out->appendColumn("arrays",string("x_rel_sat_ecl/m"));
  out->appendColumn("arrays",string("y_rel_sat_ecl/m"));
  out->appendColumn("arrays",string("z_rel_sat_ecl/m"));
  out->appendColumn("arrays",string("angdiam_rel_sat/arcsec"));
  out->appendColumn("arrays",string("d_sun/m"));
  out->appendColumn("arrays",string("angle_S_T_O/deg"));

  int64 nsamp=-1;
  add_planet(*out, params.find<string>("sun"), "Sun", nsamp);
  add_planet(*out, params.find<string>("mercury"), "Mercury", nsamp);
  add_planet(*out, params.find<string>("venus"), "Venus", nsamp);
  add_planet(*out, params.find<string>("earth"), "Earth", nsamp);
  add_planet(*out, params.find<string>("mars"), "Mars", nsamp);
  add_planet(*out, params.find<string>("jupiter"), "Jupiter", nsamp);
  add_planet(*out, params.find<string>("saturn"), "Saturn", nsamp);
  add_planet(*out, params.find<string>("uranus"), "Uranus", nsamp);
  add_planet(*out, params.find<string>("neptune"), "Neptune", nsamp);
  add_planet(*out, params.find<string>("pluto"), "Pluto", nsamp);
  add_planet(*out, params.find<string>("ssb"), "Solar System Barycenter", nsamp);
  out->setKey("nsamples",nsamp);
  out->setKey("T0",double(dateandtime2planck(2009,5,15,0,0,0)));
  out->setKey("deltaT",3600.);
PLANCK_DIAGNOSIS_END
  }
