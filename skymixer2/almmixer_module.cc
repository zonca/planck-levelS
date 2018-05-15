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
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "xcomplex.h"
#include "mix_utils.h"
#include "error_handling.h"
#include "iohandle_current.h"
#include "focalplane_db.h"
#include "paramfile.h"
#include "lsconstants.h"
#include "alm.h"
#include "alm_dmcio.h"
#include "safe_ptr.h"
#include "announce.h"
#include "levels_facilities.h"
#include "string_utils.h"

using namespace std;

int almmixer_module (int argc, const char **argv)
  {
  module_startup ("almmixer", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  // sanity check
  planck_assert (!params.param_present("amx_mapname0"),
    "Key 'amx_mapname0' found, but first component must have index 1");

  focalplane_db database(params);
  string det_id = params.find<string> ("detector_id");
  double detfrq=database.getValue<double>(det_id,"nu_cen");
  double minfrq=database.getValue<double>(det_id,"nu_min");
  double maxfrq=database.getValue<double>(det_id,"nu_max");
  double detbw=maxfrq-minfrq;
  string rtype = params.find<string> ("amx_response_type");
  safe_ptr<Response_function> det_response;
  if (rtype=="DELTA")
    det_response=new Delta(detfrq,detbw);
  else if (rtype=="TOPHAT")
    det_response=new Tophat(minfrq,maxfrq);
  else if (rtype=="GAUSS")
    det_response=new Gauss(detfrq,detbw);
  else
    planck_fail("Bad amx_response_type '"+rtype+"'.");

  double antfct = speedOfLight*speedOfLight/
         (2*detfrq*detfrq*kBoltzmann*det_response->area())*1e-20;

  int lmax = params.find<int> ("amx_lmax");
  int mmax = params.find<int> ("amx_mmax",lmax);
  bool polarisation = params.find<bool> ("amx_polar");
  int nmaps = params.find<int> ("amx_nmaps");
  bool monopole = params.find<bool> ("amx_monopole",false);
  arr<safe_ptr<iohandle> > inp(nmaps);
  arr<string> mapname(nmaps);
  arr<string> maptype(nmaps);
  arr<double> preoffset(nmaps);
  arr<double> prefactor(nmaps);
  arr<double> factor(nmaps);
  arr<double> dusttemp(nmaps);
  arr<bool> has_pol(nmaps);
  for (long m=0; m<nmaps; ++m)
    {
    string suffix = dataToString(m+1);
    mapname[m] = params.find<string> ("amx_mapname"+suffix);
    maptype[m] = params.find<string> ("amx_maptype"+suffix);

    preoffset[m]=0;
    prefactor[m]=1;
    factor[m]=1;
    if (maptype[m]=="CMB")
      {
      double x=(detfrq*hPlanck)/(kBoltzmann*tcmb);
      if (monopole) preoffset[m] = 1e20*sqrt(fourpi);
      prefactor[m] = 1e20*x/((1.0-exp(-x))*tcmb);
      factor[m] = det_response->convolve(Planck_Spectrum(tcmb));
      has_pol[m] = true;
      }
    else if (maptype[m]=="CMB_UNPOL")
      {
      double x=(detfrq*hPlanck)/(kBoltzmann*tcmb);
      if (monopole) preoffset[m] = 1e20*sqrt(fourpi);
      prefactor[m] = 1e20*x/((1.0-exp(-x))*tcmb);
      factor[m] = det_response->convolve(Planck_Spectrum(tcmb));
      has_pol[m] = false;
      }
    else if (maptype[m]=="DUST")
      {
      dusttemp[m] = params.find<double>("amx_dust_temp"+suffix);
      factor[m] = det_response->convolve
        (Dust_Spectrum(dusttemp[m], speedOfLight/1e-4, 2));
      has_pol[m] = false;
      }
    else if (maptype[m]=="DUST2")
      {
      factor[m] = det_response->convolve(Dust2_Spectrum());
      has_pol[m] = false;
      }
    else if (maptype[m]=="SZ")
      {
      prefactor[m] = 1e20;
      factor[m] = det_response->convolve(SZ_Spectrum(tcmb));
      has_pol[m] = false;
      }
    else if (maptype[m]=="SZKIN")
      {
      prefactor[m] = 1e20;
      factor[m] = det_response->convolve(SZkin_Spectrum(tcmb));
      has_pol[m] = false;
      }
    else if (maptype[m]=="SYNCHRO")
      {
      prefactor[m] = pow(22/0.408,1.25-0.75);
      factor[m] = det_response->convolve(Power_Spectrum(-1.25,0.408e9));
      has_pol[m] = true;
      }
    else if (maptype[m]=="FREEFREE")
      {
      prefactor[m] = 1e20;
      factor[m] = det_response->convolve(FreeFree_Spectrum(1e4));
      has_pol[m] = false;
      }
    else if (maptype[m]=="CO")
      {
      prefactor[m] = 1000/speedOfLight * 1e20 * det_response->area();
      double t_co = params.find<double>("amx_co_temp"+suffix);

      factor[m] = det_response->convolve(CO_Spectrum(t_co));
      has_pol[m] = false;
      }
    else if (maptype[m]=="MJY")
      {
      factor[m] = det_response->area();
      has_pol[m] = false;
      }
    else if (maptype[m]=="MJY_POL")
      {
      factor[m] = det_response->area();
      has_pol[m] = true;
      }
    else if (maptype[m]=="KRJ")
      {
      factor[m] = 1./antfct;
      has_pol[m] = false;
      }
    else if (maptype[m]=="KRJ_POL")
      {
      factor[m] = 1./antfct;
      has_pol[m] = true;
      }
    else
      planck_fail ("unknown map type '" + maptype[m] + "'");

    inp[m] = HandleManager.openObject(mapname[m],
      (has_pol[m]) ? "alm.LS_alm_pol" : "alm.LS_alm");
    }

  string outfname = params.find<string> ("amx_output_map");
  safe_ptr<iohandle> out (HandleManager.createObject(outfname,
    polarisation ? "alm.LS_alm_pol" : "alm.LS_alm"));

  Alm<xcomplex<float> > alm_in(lmax,mmax), alm_out(lmax,mmax);
  alm_out.SetToZero();
  for (long j=0; j<nmaps; ++j)
    {
    read_Alm_from_dmc (*inp[j], alm_in, lmax, mmax, 'T');
    alm_in.Scale (prefactor[j]);
    alm_in(0,0) += float(preoffset[j]);
    alm_in.Scale (factor[j]);
    for (int m=0; m<=mmax; ++m)
      for (int l=m; l<=lmax; ++l)
        alm_out(l,m)  += alm_in(l,m);
    }
  alm_out.Scale (antfct);
  write_Alm_to_dmc (*out, alm_out, lmax, mmax, 'T');

  if (polarisation)
    {
    for (int hdu=3; hdu<=4; ++hdu)
      {
      char suffix = (hdu==3) ? 'G' : 'C';
      alm_out.SetToZero();
      for (long j=0; j<nmaps; ++j)
        {
        if (has_pol[j])
          {
          read_Alm_from_dmc (*inp[j], alm_in, lmax, mmax, suffix);
          alm_in.Scale (prefactor[j]*factor[j]);
          for (int m=0; m<=mmax; ++m)
            for (int l=m; l<=lmax; ++l)
              alm_out(l,m)  += alm_in(l,m);
          }
        }
      alm_out.Scale (antfct);
      write_Alm_to_dmc (*out, alm_out, lmax, mmax, suffix);
      }
    }
  return 0;
  }
