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
 *  This file contains the sampler module for Planck originally written
 *  by Ian Grivell and Bob Mann. For details see the note "Sampler module
 *  for scan strategy simulations" in LiveLink.
 *
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "sampler.h"
#include "arr.h"
#include "focalplane_db.h"
#include "paramfile.h"

using namespace std;

Sampler::Sampler (paramfile &params, focalplane_db &fpdb, const string &det_id,
  int n_integ, double fsamp1)
  : f_samp1(fsamp1)
  {
  planck_assert(f_samp1>0.,"invalid f_samp1 in sampler");
  double tau_int = fpdb.getValue<double>(det_id,"tau_int");
  double tau_bol = fpdb.getValue<double>(det_id,"tau_bol");
  f_samp_det = fpdb.getValue<double>(det_id,"f_samp");
  int nread = fpdb.getValue<int>(det_id,"nread");
  timeshift = params.find<double>("sampler_timeshift",0.)
              * f_samp1/f_samp_det;
  // if time response is instantaneous, don't need to integrate against it
  if (tau_bol<1e-4) n_integ=1;
  int n_used = 1 + int (tau_int*f_samp_det*nread);

  arr<double> toffset(n_integ*n_used);
  deltat.alloc(n_integ*n_used);

  for (int i=0; i<n_used; ++i)
    for (int j=0; j<n_integ; ++j)
      toffset[j+i*n_integ] = f_samp1 *
        (i/(f_samp_det*nread) - tau_bol*log((double(j+1))/n_integ));
  toffset.sort();

  double avg=0;
  for (tsize m=0; m<toffset.size(); ++m)
    avg+=toffset[m];
  avg*=f_samp_det/(f_samp1*toffset.size());
  cout << "Sampler: avg offset = " << avg-timeshift*f_samp_det/f_samp1
       << " samples" << endl;

  deltat[0] = 0;
  for (tsize m=1; m<deltat.size(); ++m)
    deltat[m]=toffset[m]-toffset[m-1];

  scale = 1./deltat.size();
  }

void Sampler::sample (const arr<double> &in, const arr<double> &t1,
  const arr<double> &t2, arr<double> &out) const
  {
  planck_assert(in.size()==t1.size(),"incompatible array sizes");
  int insize = in.size();
  int nsamp = t2.size();
  out.alloc(nsamp);

#pragma omp parallel
{
  int m;
#pragma omp for schedule (static)
  for (m=0; m<nsamp; ++m)
    {
    double t = (t2[m]-t1[0])*f_samp1 + timeshift;
    double value = 0;
    // loop over all time offsets:
    int ind = max(0,min(int(t),insize-2));
    double frac = t-ind;

    for (tsize i=0; i<deltat.size(); ++i)
      {
      frac -= deltat[i];
      while ((frac<0.) && (ind>0))
        { --ind; frac+=1; }
      value += in[ind]*(1-frac) + in[ind+1]*frac;
      }
    out[m]=value*scale;
    }
}
  }
