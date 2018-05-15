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
 *  Copyright (C) 2004-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "announce.h"
#include "iohandle_current.h"
#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "paramfile.h"
#include "string_utils.h"

using namespace std;

namespace {

int find_nmaps(paramfile &params)
  {
  planck_assert (!params.param_present("map0"),
    "Key 'map0' found, but first component must have index 1");
  planck_assert (!params.param_present("factor0"),
    "Key 'factor0' found, but first component must have index 1");
  int res=0;
  while (true)
    {
    if (params.param_present(string("map")+dataToString(res+1)))
      ++res;
    else
      return res;
    }
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup("addMaps",(argc==2)||(((argc&1)==1)&&(argc>=5)),
    "Usage: addMaps <init object>\n"
    "or:    addMaps <file1> <factor1> [<file2> <factor2>] [...] [output file] [pol]");

  string outfile;
  bool pol;
  arr<string> mapname;
  arr<double> factor;

  iohandle_current::Manager mng ((argc==2) ? argv[1] : "");

  if (argc==2)
    {
    paramfile params (mng.getParams());

    pol = params.find<bool>("polarisation",false);
    outfile = params.find<string>("outfile");

    int nmaps = find_nmaps(params);
    cout << "Co-adding " << nmaps << " maps." << endl;
    mapname.alloc(nmaps);
    factor.alloc(nmaps);
    for (int i=0; i<nmaps; ++i)
      {
      string suffix = dataToString(i+1);
      mapname[i] = params.find<string>(string("map")+suffix);
      factor[i] = params.find<double>(string("factor")+suffix);
      }
    }
  else
    {
    int nmaps = (argc-3)/2;
    mapname.alloc(nmaps);
    factor.alloc(nmaps);
    for (int i=0; i<nmaps; ++i)
      {
      mapname[i] = argv[2*i+1];
      factor[i] = stringToData<double>(argv[2*i+2]);
      }
    outfile = argv[argc-2];
    pol = stringToData<bool>(argv[argc-1]);
    }

  arr<string> cname(3);
  cname[0] = "I_Stokes"; cname[1] = "Q_Stokes"; cname[2] = "U_Stokes";

  Healpix_Map<double> outmap[3], tmap;
  int ncomp = pol ? 3:1;
  for (tsize n=0; n<mapname.size(); ++n)
    {
    for (int i=0; i<ncomp; ++i)
      {
      read_Healpix_map_from_dmc (mapname[n], tmap, cname[i], false);
      if (n==0)
        {
        outmap[i].SetNside(tmap.Nside(),tmap.Scheme());
        outmap[i].fill(0);
        }
      else
        planck_assert (outmap[i].conformable(tmap),
          "maps are not conformable");

      for (int m=0; m<tmap.Npix(); ++m)
        if (approx(tmap[m],Healpix_undef)
            || approx(outmap[i][m],Healpix_undef))
          outmap[i][m] = Healpix_undef;
        else
          outmap[i][m] += factor[n]*tmap[m];
      }
    }

  if (pol)
    write_Healpix_map_to_dmc (outfile, outmap[0], outmap[1], outmap[2]);
  else
    write_Healpix_map_to_dmc (outfile, outmap[0]);

PLANCK_DIAGNOSIS_END
  }
