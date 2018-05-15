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
#include "powspec.h"
#include "powspec_dmcio.h"
#include "string_utils.h"

using namespace std;

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("text2powspec", argc, argv);
  iohandle_current::Manager mng (argc,argv);
  paramfile params (mng.getParams());

  ifstream inp(params.find<string>("in").c_str());
  planck_assert(inp,"could not open input file");

  arr<vector<double> > data;

  tsize nc=0;
  while(inp)
    {
    string line;
    getline(inp,line);
    vector<double> list;
    split(line,list);
    if (list.size()==0) continue; // empty line
    nc = list.size();
    if (data.size()==0)
      {
      planck_assert((nc==1)||(nc==4)||(nc==6),
        "bad number of spectrum components");
      data.alloc(nc);
      }
    else
      planck_assert(nc==data.size(),"number of components is changing");
    for (tsize i=0; i<nc; ++i)
      data[i].push_back(list[i]);
    }

  tsize lmax = data[0].size()-1;
  PowSpec ps(nc,lmax);
  for (tsize l=0; l<=lmax; ++l)
    {
    ps.tt(l)=data[0][l];
    if (nc==1) continue;
    ps.gg(l)=data[1][l];
    ps.cc(l)=data[2][l];
    ps.tg(l)=data[3][l];
    if (nc==4) continue;
    ps.tc(l)=data[4][l];
    ps.gc(l)=data[5][l];
    }
  write_powspec_to_dmc (params.find<string>("out"),ps,nc);
PLANCK_DIAGNOSIS_END
  }
