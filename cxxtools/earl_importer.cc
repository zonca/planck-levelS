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
 *  Copyright (C) 2011-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <fstream>
#include "iohandle_current.h"
#include "announce.h"
#include "pointing.h"
#include "string_utils.h"

using namespace std;

int main(int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("earl_importer", argc, argv);

  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  ifstream inp(params.find<string>("infile").c_str());

  vector<double> signal, noise, theta, phi;

  while (!inp.eof())
    {
    string line;
    getline(inp,line);
    vector<double> list;
    split(line,list);
    if (list.size()==0) continue; // empty line
    planck_assert(list.size()==5,"incorrect input file format");
    signal.push_back(list[0]);
    noise.push_back(list[1]);
    pointing ptg (vec3(list[2],list[3],list[4]));
    theta.push_back(ptg.theta);
    phi.push_back(ptg.phi);
    }

  safe_ptr<iohandle> out_tod(HandleManager.createObject
    (params.find<string>("out_tod"),"toi.LS_toi"));
  out_tod->appendColumn("signal",signal);
  safe_ptr<iohandle> out_noise(HandleManager.createObject
    (params.find<string>("out_noise"),"toi.LS_toi"));
  out_noise->appendColumn("signal",noise);
  safe_ptr<iohandle> out_ptg(HandleManager.createObject
    (params.find<string>("out_pointing"),"pointing.LS_detpoint"));
  out_ptg->appendColumn("theta",theta);
  out_ptg->appendColumn("phi",phi);
PLANCK_DIAGNOSIS_END
  }
