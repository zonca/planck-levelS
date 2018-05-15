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
 *  Copyright (C) 2005-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "iohandle_current.h"
#include "io_utils.h"
#include "pointing.h"
#include "arr.h"
#include "paramfile.h"
#include "lsconstants.h"
#include "announce.h"

using namespace std;

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("pntsrc2map",argc,argv);
  iohandle_current::Manager mng (argc,argv);
  paramfile params (mng.getParams());

  int nside = params.find<int>("nside");
  string infile = params.find<string>("catalog");
  string outfile = params.find<string>("outfile");
  bool pol = params.find<bool>("polarisation",false);

  Healpix_Map<double> mapT (nside,RING,SET_NSIDE),
                      mapQ (nside,RING,SET_NSIDE),
                      mapU (nside,RING,SET_NSIDE);
  mapT.fill(0); mapQ.fill(0); mapU.fill(0);

  arr<pointing> pos;
  arr<double> flux, frac, angle;
  arr<string> name;
  readPointSourcesNew (infile, pol, pos, flux, angle, frac, name, 1, 0);

  double factor = mapT.Npix()/fourpi;

  for (tsize i=0; i<flux.size(); ++i)
    {
    tsize pix = mapT.ang2pix(pos[i]);
    mapT[pix] += flux[i]*factor;
    if (pol)
      {
      mapQ[pix] += flux[i]*frac[i]*factor*cos(2*angle[i]);
      mapU[pix] += flux[i]*frac[i]*factor*sin(2*angle[i]);
      }
    }

  if (pol)
    write_Healpix_map_to_dmc (outfile, mapT, mapQ, mapU);
  else
    write_Healpix_map_to_dmc (outfile, mapT);
PLANCK_DIAGNOSIS_END
  }
