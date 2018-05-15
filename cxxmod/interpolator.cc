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
 *  This file implements the angular interpolation of the totalconvolver
 *  output.
 *
 *  Copyright (C) 2003-2010 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "xcomplex.h"
#include "pointing.h"
#include "interpolator.h"
#include "iohandle.h"
#include "paramfile.h"
#include "trafos.h"
#include "lsconstants.h"

using namespace std;

Interpolator::Interpolator (paramfile &params, double psi_pol)
  : psi(psi_pol), ecl2gal(2000,2000,Ecliptic,Galactic)
  {
  string outstr = params.find<string>("output_type","SIGNAL");
  output=SIGNAL;
  if (outstr=="SIGNAL") output=SIGNAL;
  else if (outstr=="I") output=I;
  else if (outstr=="Q") output=Q;
  else if (outstr=="U") output=U;
  else planck_fail ("Unknown value '" + outstr + "' for output_type");

  string ringset = params.find<string>("ringset");
  int order = params.find<int>("interpol_order",1);
  galactic = params.find<bool>("interpol_galactic",false);
  init (ringset, order);
  }

Interpolator::Interpolator (const string &ringset, double psi_pol,
  int order)
  : psi(psi_pol),galactic(false),ecl2gal(2000,2000,Ecliptic,Galactic)
  {
  output=SIGNAL;
  init (ringset, order);
  }

void Interpolator::init (const string &ringset, int order)
  {
  safe_ptr<iohandle> inp
    (HandleManager.openObject(ringset,"ringset.LS_ringset"));
  inp->getKey ("beam_mmax",beammmax);
  npsi=0;
  psiarr.alloc(beammmax+1);
  psiarr.fill(-1);
  int col_m = inp->columnNumber("ringsets_present");
  int num_m = inp->columnLength(col_m);
  for (int i=0; i<num_m; ++i)
    {
    int m;
    inp->readColumn(col_m, m, i);
    psiarr[m] = npsi;
    (m==0) ? npsi+=1 : npsi+=2;
    }

  inp->getKey("nphi",nphi);
  double dphi=degr2rad*inp->getKey<double>("dphi");
  inv_delta_phi = 1./dphi;
  double phi0=degr2rad*inp->getKey<double>("phi0");
  phioffset = phi0/dphi;
  planck_assert(approx(abs(nphi/inv_delta_phi),twopi), "phi extent is not 2 pi");
  inp->getKey("ntheta",ntheta);
  double dtheta=degr2rad*inp->getKey<double>("dtheta");
  inv_delta_theta = 1./dtheta;
  double theta0=degr2rad*inp->getKey<double>("theta0");
  thetaoffset = theta0/dtheta;

  sky.alloc(nphi);
  for (tsize i=0; i<sky.size(); ++i)
    sky[i].alloc(ntheta,npsi);

  for (int k=0; k<npsi; ++k)
    {
    arr2<float> tmpsky(ntheta,nphi);
    inp->readColumnRaw("ringsetdata",&(tmpsky[0][0]),
      ntheta*nphi,int64(k)*ntheta*nphi);
    for (int i=0; i<ntheta; ++i)
      for (int j=0; j<nphi; ++j)
        sky[j][i][k] = (k==0) ? tmpsky[i][j] : 2*tmpsky[i][j];
    }

  planck_assert ((order>0) && (order&1) && (order<=max_order),
    "bad interpolation order");
  npoints = order+1;
  ioffset = order/2;
  base_wgt.alloc(npoints);
  base_wgt.fill(1);
  for (int m=0; m<npoints; ++m)
    for (int n=0; n<npoints; ++n)
      if (m!=n) base_wgt[m] *= m-n;
  for (int m=0; m<npoints; ++m) base_wgt[m] = 1./base_wgt[m];
  }

inline void Interpolator::weight_n (double x, double *wgt) const
  {
  for (int m=0; m<npoints; ++m)
    wgt[m] = base_wgt[m];
  double mul1=x, mul2=x-npoints+1;
  for (int m=1; m<npoints; ++m)
    {
    wgt[m]*=mul1;
    wgt[npoints-m-1]*=mul2;
    mul1*=x-m;
    mul2*=x-npoints+m+1;
    }
  }

void Interpolator::interpol_n (double theta, double phi,
  arr<double> &result) const
  {
  double wgt1[max_order+1];
  double frac = theta*inv_delta_theta - thetaoffset;
  int itheta0 = int (frac) - ioffset;
  if (itheta0>(ntheta-npoints)) itheta0 = ntheta-npoints;
  if (itheta0<0) itheta0 = 0;
  frac -= itheta0;
  weight_n (frac,wgt1);

  double wgt2[max_order+1];
  frac = phi*inv_delta_phi - phioffset;
  frac = fmodulo (frac,double(nphi));
  int iphi0 = int (frac) - ioffset;
  frac -= iphi0;
  if (iphi0 >= nphi) iphi0-=nphi;
  if (iphi0 < 0) iphi0+=nphi;
  weight_n (frac,wgt2);

  result.fill(0);
  int iphi = iphi0;
  for (int i=0; i<npoints; ++i)
    {
    for (int j=0; j<npoints; ++j)
      {
      double weight = wgt2[i]*wgt1[j];
      const float *ref = sky[iphi][itheta0+j];
      for (int k=0; k<npsi; ++k)
        result[k] += weight*ref[k];
      }
    if (++iphi>=nphi) iphi-=nphi;
    }
  }

double Interpolator::interpol_psi(double omega, const arr<double> &karr) const
  {
  double result = (psiarr[0]>=0) ? karr[0] : 0;
  if (npsi<=1) return result;

  double cosang=1, sinang=0;
  double sinomg=sin(omega), cosomg=cos(omega);

  for (int k=1; k<=beammmax; ++k)
    {
    const double tmp = sinang*cosomg + cosang*sinomg;
    cosang=cosang*cosomg - sinang*sinomg;
    sinang=tmp;
    if (psiarr[k]>=0)
      result += cosang*karr[psiarr[k]] - sinang*karr[psiarr[k]+1];
    }
  return result;
  }

void Interpolator::Get_Intensities (const arr<pointing> &detpt,
  arr<double> &ii, arr<double> &iq, arr<double> &iu) const
  {
  ii.alloc(detpt.size());
  iq.alloc(detpt.size());
  iu.alloc(detpt.size());

#pragma omp parallel
{
  arr<double> karr (npsi);
  int m, sz=detpt.size();
#pragma omp for schedule (static)
  for (m=0; m<sz; ++m)
    {
    interpol_n (detpt[m].theta, detpt[m].phi, karr);

    double val0  = interpol_psi(-psi, karr),
           val90 = interpol_psi(halfpi-psi,karr);
    ii[m] = 0.5 * (val0+val90);
    iq[m] = 0.5 * (val90-val0);
    iu[m] = 0.5 * (interpol_psi(1.5*halfpi-psi, karr)
                  -interpol_psi(0.5*halfpi-psi,karr));
    }
}
  }

void Interpolator::Add_Intensities (const arr<pointing> &detpt,
  const arr<double> &heading, arr<double> &intensity) const
  {
  planck_assert(multiequal(intensity.size(),heading.size(),detpt.size()),
    "wrong array sizes");

#pragma omp parallel
{
  arr<double> karr (npsi);
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
    interpol_n (ptg.theta, ptg.phi, karr);

    switch (output)
      {
      case SIGNAL:
        {
        double omega_wg = hdg-halfpi;
        intensity[m] += interpol_psi (omega_wg, karr);
        }
        break;
      case I:
        intensity[m] += 0.5 * (interpol_psi(0,karr)+interpol_psi(halfpi,karr));
        break;
      case Q:
        intensity[m] += 0.5 * (interpol_psi(halfpi-psi, karr)
                              -interpol_psi(-psi,karr));
        break;
      case U:
        intensity[m] += 0.5 * (interpol_psi(1.5*halfpi-psi, karr)
                              -interpol_psi(0.5*halfpi-psi, karr));
        break;
      }
    }
}
  }
