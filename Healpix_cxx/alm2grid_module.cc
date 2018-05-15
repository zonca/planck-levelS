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
 *  Copyright (C) 2004-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "xcomplex.h"
#include "paramfile.h"
#include "alm.h"
#include "alm_dmcio.h"
#include "arr.h"
#include "iohandle_current.h"
#include "sharp_cxx.h"
#include "levels_facilities.h"
#include "lsconstants.h"
#include "announce.h"

using namespace std;

namespace {

template<typename T> inline void stokesrotate (T &q, T &u, double phi)
  {
  double c2=cos(2*phi), s2=sin(2*phi);
  double q2=q, u2=u;
  q = T(-c2*q2 + s2*u2);
  u = T(-s2*q2 - c2*u2);
  }

template<typename T> void alm2grid (const Alm<xcomplex<T> > &almT,
  const Alm<xcomplex<T> > &almG, const Alm<xcomplex<T> > &almC,
  arr2<T> &gridI, arr2<T> &gridQ, arr2<T> &gridU, bool polarisation,
  double thetamax)
  {
  int ntheta = gridI.size1();
  int nphi_ = gridI.size2();
  arr<int> nphi(ntheta), stride(ntheta);
  arr<ptrdiff_t> ofs(ntheta);
  arr<double> phi0(ntheta), theta(ntheta), weight(ntheta);

  for (int ith=0; ith<ntheta; ++ith)
    {
    theta[ith] = double(ith)*thetamax/(ntheta-1);
    // move first ring a tiny little bit away from the pole
    if (ith==0) theta[ith] = 0.001*thetamax/(ntheta-1);
    phi0[ith]=0;
    ofs[ith]=&gridI[ith][0]-&gridI[0][0];
    stride[ith]=1;
    weight[ith]=0;
    nphi[ith]=nphi_;
    }

  sharp_cxxjob<T> job;
  job.set_general_geometry (ntheta, &nphi[0], &ofs[0], &stride[0],
    &phi0[0], &theta[0], &weight[0]);
  job.set_triangular_alm_info (almT.Lmax(), almT.Mmax());
  job.alm2map(&almT(0,0).re, &gridI[0][0], false);
  if (polarisation)
    job.alm2map_spin(&almG(0,0).re, &almC(0,0).re,
      &gridQ[0][0], &gridU[0][0], 2, false);

  // going from polar basis to co-cross basis
  if (polarisation)
    {
    double dphi = twopi / nphi_;
    for (int ith=0; ith<ntheta; ++ith)
      for (int iph=0; iph<nphi_; ++iph)
        stokesrotate(gridQ[ith][iph],gridU[ith][iph],iph*dphi);
    }
  }

} // unnamed namespace

int alm2grid_module (int argc, const char **argv)
  {
  module_startup ("alm2grid", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  int nlmax = params.find<int>("nlmax");
  int nmmax = params.find<int>("nmmax",nlmax);
  planck_assert(nmmax<=nlmax,"nmmax must not be larger than nlmax");
  string infile = params.find<string>("infile");
  string outfile = params.find<string>("outfile");
  int ntheta = params.find<int>("ntheta");
  int nphi = params.find<int>("nphi");
  double thetamax = params.find<double>("thetamax")*degr2rad;
  bool polarisation = params.find<bool>("polarisation",false);

  Alm<xcomplex<float> > almT, almG, almC;
  if (polarisation)
    read_Alm_from_dmc(infile,almT,almG,almC,nlmax,nmmax);
  else
    read_Alm_from_dmc(infile,almT,nlmax,nmmax);

  arr2<float> gridI(ntheta,nphi), gridQ, gridU;
  if (polarisation)
    { gridQ.alloc(ntheta,nphi); gridU.alloc(ntheta,nphi); }

  float offset = float(almT(0,0).real()/sqrt(fourpi));
  almT(0,0) = 0;
  alm2grid(almT,almG,almC,gridI,gridQ,gridU,polarisation,thetamax);
  for (tsize m=0; m<gridI.size1(); ++m)
    for (tsize n=0; n<gridI.size2(); ++n)
      gridI[m][n]+=offset;
  safe_ptr<iohandle> out (HandleManager.createObject(outfile,
    polarisation ? "beam.LS_beammap_pol" : "beam.LS_beammap"));

  out->setKey("Ntheta",ntheta);
  out->setKey("Nphi",nphi);
  out->setKey("Maxtheta",thetamax);
  out->appendColumnRaw("Beamdata",&(gridI[0][0]), gridI.size());
  if (polarisation)
    {
    out->appendColumnRaw("BeamdataQ",&(gridQ[0][0]), gridQ.size());
    out->appendColumnRaw("BeamdataU",&(gridU[0][0]), gridU.size());
    }
  return 0;
  }
