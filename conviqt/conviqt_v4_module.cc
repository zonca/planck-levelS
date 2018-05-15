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
 *  This file implements the total convolution algorithm suggested by
 *  Gary Prezeau (Prezeau & Reinecke 2010, http://arxiv.org/pdf/1002.1050)
 *
 *  Copyright (C) 2009-2013 Max-Planck-Society
 *  Authors: Martin Reinecke, Gary Prezeau
 */

#include "levels_facilities.h"
#include "alm.h"
#include "paramfile.h"
#include "mpi_support.h"
#include "iohandle_current.h"
#include "lsconstants.h"
#include "announce.h"
#include "sharp_cxx.h"
#include "convolver_helper.h"

using namespace std;

namespace {

struct ringpair
  {
  int r1,r2;
  ringpair (int i1, int i2) : r1(i1), r2(i2) {}
  };

template <typename T> void sort_theta(const arr<double> &th, arr<double> &th2,
  arr2<T> &img)
  {
  th2=th;
  for (tsize i=0; i<th2.size()-1; ++i)
    {
    tsize imin = i;
    for (tsize j=i+1; j<th2.size(); ++j)
      if (th2[j]>th2[imin]) imin = j; // decreasing thetas for TC compatibility

    if (imin!=i)
      {
      swap (th2[imin],th2[i]);
      for (tsize k=0; k<img.size2(); ++k) swap (img[imin][k],img[i][k]);
      }
    }
  }

template <typename T> void write_data (iohandle &out, arr2<T> &imgarr,
  const arr<double> &allthetas, bool &keys_written)
  {
  arr<double> th2;
  sort_theta (allthetas, th2, imgarr);
  out.appendColumnRaw("ringsetdata",&(imgarr[0][0]),imgarr.size());
  if (keys_written) return;
  out.setKey("nphi",int32(imgarr.size2()));
  out.setKey("phi0",90.);
  out.setKey("dphi",360./imgarr.size2());
  out.setKey("ntheta",int32(imgarr.size1()));
  out.setKey("theta0",th2[0]*rad2degr);
  out.setKey("dtheta",(th2[1]-th2[0])*rad2degr);
  keys_written = true;
  }

template <typename T> void write_file (arr2<T> &imgarr,
  iohandle &out, const arr<double> &allthetas, bool &keys_written)
  {
  if (mpiMgr.num_ranks()>1)
    {
    if (mpiMgr.master())
      {
      arr2<T> imgarr2;
      mpiMgr.gatherv_m(imgarr,imgarr2);
      write_data(out,imgarr2,allthetas,keys_written);
      }
    else
      mpiMgr.gatherv_s(imgarr);
    }
  else
    write_data(out,imgarr,allthetas,keys_written);
  }

void make_thetas (int ntheta, arr<double> &thetas, arr<double> &allthetas,
  vector<ringpair> &pair)
  {
  pair.clear();
  const int addrows=10;
  double dtheta=pi/(ntheta-1);

  double mydtheta=mpiMgr.num_ranks()*dtheta;
  double mytheta0=dtheta*(mpiMgr.rank()-addrows);
  double dntheta=1 + (halfpi-mytheta0)/mydtheta;
  int myntheta=int(dntheta-1e-7); // don't include a potential ring at pi/2
  if (approx(dntheta,floor(dntheta+1e-7),1e-10)) // ring at pi/2 exists
    {
    thetas.alloc(2*myntheta+1);
    thetas[2*myntheta]=halfpi;
    pair.push_back(ringpair(2*myntheta,-1));
    }
  else
    thetas.alloc(2*myntheta);

  for (int i=0; i<myntheta; ++i)
    {
    thetas[2*i] = mytheta0+i*mydtheta;
    thetas[2*i+1] = pi-thetas[2*i];
    pair.push_back(ringpair(2*i,2*i+1));
    }

  mpiMgr.gatherv(thetas,allthetas);
  }

template <typename T> void convolve (paramfile &params)
  {
  convolver_helper<T> chelper(params);
  int lmax = chelper.Lmax();
  int olmax = params.find<int>("lmax_out",lmax);
  if (olmax<0) olmax=lmax;
  int nphi = params.find<int>("nphi",2*olmax+1);
  if (nphi<0) nphi=2*olmax+1;
  olmax=(nphi-1)/2;
  planck_assert (olmax>=lmax, "lmax_out is smaller than conv_lmax");
  int ntheta = params.find<int>("ntheta",olmax+2);
  if (ntheta<0) ntheta=olmax+2;

  arr<double> thetas, allthetas;
  vector<ringpair> pairs;
  make_thetas (ntheta, thetas, allthetas, pairs);

  safe_ptr<iohandle> out;
  vector<int32> present_sets;
  bool dp_output = params.find<bool>("double_precision_output",false);
  if (mpiMgr.master())
    {
    out = HandleManager.createObject (params.find<string>("ringset"),
      dp_output ? "ringset.LS_ringset_dp" : "ringset.LS_ringset");
    out->setKey ("beam_mmax",int32(chelper.beamMmax()));
    }
  bool keys_written=false;

  arr2<T> map1(thetas.size(),nphi),map2;
  Alm<xcomplex<T> > a1,a2;

  tsize nth=thetas.size();
  sharp_cxxjob<T> job;

  arr<int> nph(nth,nphi), stride(nth,1);
  arr<ptrdiff_t> ofs(nth);
  arr<double> phi0(nth,halfpi);
  for (tsize i=0; i<nth; ++i)
    ofs[i]=i*nphi;
  job.set_general_geometry (nth, &nph[0], &ofs[0], &stride[0], &phi0[0],
    &thetas[0], NULL);

  job.set_triangular_alm_info (lmax, lmax);

  for(int k=0; k<=chelper.beamMmax(); ++k)
    {
    if (chelper.blZero(k)) continue; // all terms for this k are zero

    if (k>0) map2.alloc(thetas.size(),nphi);
    chelper.getAlm(k,a1,a2);
    (k==0) ? job.alm2map(&a1(0,0).re,&map1[0][0],false) :
             job.alm2map_spin(&a1(0,0).re,&a2(0,0).re,
               &map1[0][0],&map2[0][0],k,false);

    int sign1=1,sign2=-1;
/* rotate by 90 degrees in psi */
    int quadrant=imodulo(k,4);
    if (quadrant&1)
      { map1.swap(map2); swap(sign1,sign2); }

    if ((quadrant==1)||(quadrant==2))
      sign1=-sign1;
    if ((quadrant==2)||(quadrant==3))
      sign2=-sign2;

    if (sign1!=1) map1.scale(sign1);
    if (sign2!=1) map2.scale(sign2);

    if (mpiMgr.master())
      cout << "output: mbeam=" << k << endl;
    if (mpiMgr.master())
      present_sets.push_back(int(k));
    write_file (map1,*out,allthetas,keys_written);
    if (k!=0)
      write_file (map2,*out,allthetas,keys_written);
    }

  if (mpiMgr.master())
    out->appendColumnRaw("ringsets_present",&present_sets[0],
      present_sets.size());
  }

} // unnamed namespace

int conviqt_v4_module (int argc, const char **argv)
  {
  module_startup ("conviqt_v4", argc, argv, mpiMgr.master());
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams(mpiMgr.master()));

  params.find<bool>("double_precision",false) ?
    convolve<double>(params) : convolve<float>(params);

  return 0;
  }
