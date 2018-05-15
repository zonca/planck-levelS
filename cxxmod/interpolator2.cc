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
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "healpix_map_dmcio.h"
#include "pointing.h"
#include "interpolator2.h"
#include "paramfile.h"
#include "focalplane_db.h"

using namespace std;

Interpolator2::Interpolator2 (paramfile &params, focalplane_db &fpdb,
  const string &det_id, double psi_pol)
  : psi(psi_pol), ecl2gal(2000,2000,Ecliptic,Galactic)
  {
  string outstr = params.find<string>("interpol2_output_type","SIGNAL");
  output=SIGNAL;
  if (outstr=="SIGNAL") output=SIGNAL;
  else if (outstr=="I") output=I;
  else if (outstr=="Q") output=Q;
  else if (outstr=="U") output=U;
  else planck_fail ("Unknown value '" + outstr + "' for interpol2_output_type");

  string mapname = params.find<string>("interpol2_map");
  polarisation = params.find<bool>("interpol2_polarisation",false);
  planck_assert ((psi_pol<1e10)||(!polarisation),
    "Detector is non-polarised, but polarisation was requested.");
  galactic = params.find<bool>("interpol2_galactic",false);

  if (polarisation)
    read_Healpix_map_from_dmc(mapname,mapT,mapQ,mapU);
  else
    read_Healpix_map_from_dmc(mapname,mapT);

  double epsilon = fpdb.getValue<double>(det_id,"epsilon");
  mapT.Scale(float(0.5*(1+epsilon)));
  if (polarisation)
    {
    mapQ.Scale(float(0.5*(1-epsilon)));
    mapU.Scale(float(0.5*(1-epsilon)));
    }
  }

void Interpolator2::Add_Intensities (const arr<pointing> &detpt,
  const arr<double> &heading, arr<double> &intensity) const
  {
  planck_assert(multiequal(intensity.size(),heading.size(),detpt.size()),
    "wrong array sizes");

  if ((!polarisation) && ((output==Q)||(output==U))) return;

#pragma omp parallel
{
  int m, sz=detpt.size();
#pragma omp for schedule (static)
  for (m=0; m<sz; ++m)
    {
    pointing ptg;
    double hdg;
    if (galactic)
      {
      ecl2gal.rotatefull (detpt[m],ptg,hdg);
      hdg += heading[m];
      }
    else
      {
      ptg = detpt[m];
      hdg = heading[m];
      }

    int pix = mapT.ang2pix(ptg);

    switch (output)
      {
      case SIGNAL:
        {
        intensity[m] += mapT[pix];
        if (polarisation)
          {
          double omega = 2*(hdg+psi);
          intensity[m] += mapQ[pix]*cos(omega) + mapU[pix]*sin(omega);
          }
        }
        break;
      case I:
        intensity[m] += mapT[pix];
        break;
      case Q:
        intensity[m] += mapQ[pix];
        break;
      case U:
        intensity[m] += mapU[pix];
        break;
      }
    }
}
  }
