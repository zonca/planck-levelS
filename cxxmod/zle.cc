/*
 *  This file is part of libcxxmod.
 *
 *  libcxxmod is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxmod is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxmod; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxmod is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2007-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <map>
#include <vector>
#include "zle.h"
#include "arr.h"
#include "focalplane_db.h"
#include "paramfile.h"
#include "zod_interface.h"
#include "healpix_base.h"
#include "string_utils.h"

using namespace std;

#ifdef MIXBIN_SUPPORTED

ZLE::ZLE (focalplane_db &focal, const string &det_id, paramfile &params)
  {
  vector<int> comp_i;
  split(params.find<string>("zle_components"),comp_i);
  vector<double> fact_i(comp_i.size(),1.);
  if (params.param_present("zle_factors"))
    split(params.find<string>("zle_factors"),fact_i);
  planck_assert(fact_i.size()==comp_i.size(),
    "wrong number of factors in 'zle_factors'");
  zle_init (focal.getValue<double>(det_id,"nu_cen")*1e-9,&comp_i[0],&fact_i[0],
    comp_i.size(),-1e-3,13,0.,5.2);
  ephSun = getEphemerides(params.find<string>("zle_ephemeris"));
  ephBary = getEphemerides(params.find<string>("zle_ephemeris"));
  ephSun->loadBody("Sun");
  ephBary->loadBody("Solar System Barycenter");
  nside = params.find<int>("zle_nside",128);
  }
ZLE::~ZLE()
  { zle_destroy(); }

void ZLE::Add_Intensities (const arr<pointing> &detpt, const arr<vec3> &vdetpt,
  const arr<double> &times, arr<double> &intensity) const
  {
  tsize nsamp=vdetpt.size();
  planck_assert(multiequal(nsamp,times.size(),intensity.size()),
    "array size mismatch");
  double time=0.5*(times[0]+times[nsamp-1]);
  vec3 psun = ephSun->posRelSat_m(time),
       pbary = ephBary->posRelSat_m(time);
  Healpix_Base base(nside,NEST,SET_NSIDE);
  arr<int> ppix(nsamp);
  map<int,int> pixmap;
  vector<vec3> vpix;
  vpix.reserve(10*nside);
  tsize npixhit=0;
  int last_pix=-1, last_idx=-1;
  for (tsize i=0; i<nsamp; ++i)
    {
    int pix = base.zphi2pix(vdetpt[i].z,detpt[i].phi);
    if (pix==last_pix)
      ppix[i]=last_idx;
    else
      {
      map<int,int>::iterator entry=pixmap.find(pix);
      if (entry==pixmap.end())
        {
        pixmap[pix]=ppix[i]=npixhit++;
        vpix.push_back(base.pix2vec(pix));
        }
      else
        ppix[i]=entry->second;

      last_pix=pix;
      last_idx=ppix[i];
      }
    }
  arr<double> int2(npixhit,0.);
  zle_compute(&psun.x, &pbary.x, npixhit, &vpix[0].x, &int2[0]);
  for (tsize i=0; i<nsamp; ++i)
    intensity[i]+=int2[ppix[i]];
  }

#else

ZLE::ZLE (focalplane_db &, const string &, paramfile &)
  { planck_fail("unsupported call into Fortran routine"); }
ZLE::~ZLE()
  { planck_fail("unsupported call into Fortran routine"); }

void ZLE::Add_Intensities (const arr<pointing> &, const arr<vec3> &,
  const arr<double> &, arr<double> &) const
  { planck_fail("unsupported call into Fortran routine"); }

#endif
