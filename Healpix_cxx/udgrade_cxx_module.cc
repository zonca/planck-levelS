/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "paramfile.h"
#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "iohandle_current.h"
#include "levels_facilities.h"
#include "announce.h"

using namespace std;

namespace {

template<typename T> void udgrade_cxx (paramfile &params)
  {
  string infile = params.template find<string>("infile");
  string outfile = params.template find<string>("outfile");
  int order = Healpix_Base::nside2order (params.template find<int>("nside"));
  bool polarisation = params.template find<bool>("polarisation",false);
  bool pessimistic = params.template find<bool>("pessimistic",false);

  if (!polarisation)
    {
    Healpix_Map<T> inmap;
    read_Healpix_map_from_dmc(infile,inmap);
    Healpix_Map<T> outmap (order, inmap.Scheme());

    outmap.Import(inmap,pessimistic);
    write_Healpix_map_to_dmc (outfile,outmap);
    }
  else
    {
    Healpix_Map<T> inmap;
    read_Healpix_map_from_dmc(infile,inmap);
    Healpix_Map<T> outmapT (order, inmap.Scheme()),
                   outmapQ (order, inmap.Scheme()),
                   outmapU (order, inmap.Scheme());
    planck_assert(inmap.Nside()<=outmapT.Nside(),
      "degrading not supported for polarised maps");

    outmapT.Import(inmap,pessimistic);
    read_Healpix_map_from_dmc(infile,inmap,"Q_Stokes");
    outmapQ.Import(inmap,pessimistic);
    read_Healpix_map_from_dmc(infile,inmap,"U_Stokes");
    outmapU.Import(inmap,pessimistic);
    write_Healpix_map_to_dmc (outfile,outmapT,outmapQ,outmapU);
    }
  }

} // unnamed namespace

int udgrade_cxx_module (int argc, const char **argv)
  {
  module_startup ("udgrade_cxx", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  bool dp = params.find<bool> ("double_precision",false);
  dp ? udgrade_cxx<double>(params) : udgrade_cxx<float>(params);

  return 0;
  }
