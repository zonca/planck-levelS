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

#include "arr.h"
#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "paramfile.h"
#include "focalplane_db.h"
#include "iohandle_current.h"
#include "safe_ptr.h"
#include "levels_facilities.h"
#include "lsconstants.h"
#include "announce.h"
#include "string_utils.h"

#undef CALC_FLUX

using namespace std;

namespace {

#include "spectra.cc"
#include "detector_responses.cc"

//! converts from flux in MJy/sr to antenna K at \a freq.
double antennaFact (double freq)
  {
  return speedOfLight*speedOfLight/(2*freq*freq*kBoltzmann)*1e-20;
  }

void map_import (const string &name, Healpix_Map<float> &map,
  const string &colname)
  {
  Healpix_Map<float> tmap;
  read_Healpix_map_from_dmc (name,tmap,colname,false);
  map.Import(tmap);
  }

#include "emitters.cc"

} // unnamed namespace

int skymixer3_module (int argc, const char **argv)
  {
  module_startup ("skymixer3", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  int nside = params.find<int> ("nside");
  Healpix_Map<float> mapT(nside,RING,SET_NSIDE),
                     mapQ(nside,RING,SET_NSIDE),
                     mapU(nside,RING,SET_NSIDE);
  mapT.fill(0); mapQ.fill(0); mapU.fill(0);
  safe_ptr<Detector_Response> resp (make_response(params));
  vector<Emitter *> emitters;
  collect_emitters(params,emitters);
  for (unsigned long m=0; m<emitters.size(); ++m)
    {
    emitters[m]->add_emission(*resp,mapT,mapQ,mapU);
    delete emitters[m];
    }

  string outfname = params.find<string>("output_map");
  write_Healpix_map_to_dmc (outfname,mapT,mapQ,mapU);
  return 0;
  }
