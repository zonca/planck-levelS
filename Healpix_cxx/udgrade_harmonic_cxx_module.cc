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
 *  Copyright (C) 2012-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "paramfile.h"
#include "alm.h"
#include "xcomplex.h"
#include "healpix_data_io.h"
#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "iohandle_current.h"
#include "alm_healpix_tools.h"
#include "levels_facilities.h"
#include "announce.h"
#include "lsconstants.h"

using namespace std;

namespace {

template<typename T> void udgrade_harmonic_cxx (paramfile &params)
  {
  string infile = params.template find<string>("infile");
  string outfile = params.template find<string>("outfile");
  int nlmax = params.template find<int>("nlmax");
  int nside = params.template find<int>("nside");
  string pixwin_in = params.template find<string>("pixwin_in","");
  string pixwin_out = params.template find<string>("pixwin_out","");
  int num_iter = params.template find<int>("iter_order",0);
  bool polarisation = params.template find<bool>("polarisation",false);

  if (!polarisation)
    {
    Healpix_Map<T> map;
    read_Healpix_map_from_dmc(infile,map);

    double avg=map.average();
    map.Add(T(-avg));

    arr<double> weight;
    get_ring_weights (params,map.Nside(),weight);

    Alm<xcomplex<T> > alm(nlmax,nlmax);
    if (map.Scheme()==NEST) map.swap_scheme();
    map2alm_iter(map,alm,num_iter,weight);

    arr<double> temp(nlmax+1);
    if (pixwin_in!="")
      {
      read_pixwin(pixwin_in,temp);
      for (int l=0; l<=nlmax; ++l)
        temp[l] = 1/temp[l];
      alm.ScaleL (temp);
      }
    if (pixwin_out!="")
      {
      read_pixwin(pixwin_out,temp);
      alm.ScaleL (temp);
      }

    alm(0,0) += T(avg*sqrt(fourpi));
    map.SetNside(nside,RING);
    alm2map (alm,map);

    write_Healpix_map_to_dmc (outfile,map);
    }
  else
    {
    Healpix_Map<T> mapT, mapQ, mapU;
    read_Healpix_map_from_dmc(infile,mapT,mapQ,mapU);

    double avg=mapT.average();
    mapT.Add(T(-avg));

    arr<double> weight;
    get_ring_weights (params,mapT.Nside(),weight);

    Alm<xcomplex<T> > almT(nlmax,nlmax), almG(nlmax,nlmax), almC(nlmax,nlmax);
    if (mapT.Scheme()==NEST) mapT.swap_scheme();
    if (mapQ.Scheme()==NEST) mapQ.swap_scheme();
    if (mapU.Scheme()==NEST) mapU.swap_scheme();
    map2alm_pol_iter
      (mapT,mapQ,mapU,almT,almG,almC,num_iter,weight);

    arr<double> temp(nlmax+1), pol(nlmax+1);
    if (pixwin_in!="")
      {
      read_pixwin(pixwin_in,temp,pol);
      for (int l=0; l<=nlmax; ++l)
        { temp[l] = 1/temp[l]; pol[l] = 1/pol[l]; }
      almT.ScaleL(temp); almG.ScaleL(pol); almC.ScaleL(pol);
      }
    if (pixwin_out!="")
      {
      read_pixwin(pixwin_out,temp,pol);
      almT.ScaleL(temp); almG.ScaleL(pol); almC.ScaleL(pol);
      }

    almT(0,0) += T(avg*sqrt(fourpi));
    mapT.SetNside(nside,RING);
    mapQ.SetNside(nside,RING);
    mapU.SetNside(nside,RING);
    alm2map_pol (almT,almG,almC,mapT,mapQ,mapU);

    write_Healpix_map_to_dmc (outfile,mapT,mapQ,mapU);
    }
  }

} // unnamed namespace

int udgrade_harmonic_cxx_module (int argc, const char **argv)
  {
  module_startup ("udgrade_harmonic_cxx", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  bool dp = params.find<bool> ("double_precision",false);
  dp ? udgrade_harmonic_cxx<double>(params)
     : udgrade_harmonic_cxx<float> (params);

  return 0;
  }
