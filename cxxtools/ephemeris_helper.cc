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

#include <fstream>
#include "iohandle_current.h"
#include "ephemerides.h"
#include "pointing.h"
#include "announce.h"
#include "string_utils.h"

using namespace std;

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("ephemeris_helper", argc, argv, 7,
    "<ephemeris file> <outfile> <body> <t1> <t2> <deltat>");

  iohandle_current::Manager mng ("");

  double t1=stringToData<double>(argv[4]);
  double t2=stringToData<double>(argv[5]);
  double dt=stringToData<double>(argv[6]);

  safe_ptr<ephemerides> eph (getEphemerides(argv[1]));
  eph->loadBody(argv[3]);
  ofstream out (argv[2]);
  for (double t=t1; t<t2; t+=dt)
    {
    vec3 pos = eph->posRelSat_m (t);
    pointing ptg(pos);
    out << dataToString(t) << " "
        << dataToString(ptg.theta) << " " << dataToString(ptg.phi) << endl;
    }
PLANCK_DIAGNOSIS_END
  }
