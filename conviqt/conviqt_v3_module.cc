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
 *  Copyright (C) 2009-2011 Max-Planck-Society
 *  Authors: Martin Reinecke, Gary Prezeau
 */

#include "levels_facilities.h"
#include "alm.h"
#include "alm_dmcio.h"
#include "alm_powspec_tools.h"
#include "paramfile.h"
#include "fftpack_support.h"
#include "mpi_support.h"
#include "iohandle_current.h"
#include "wigner.h"
#include "sse_utils_cxx.h"
#include "lsconstants.h"
#include "announce.h"

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

template <typename T> void write_file (const arr2<xcomplex<double> >&tmm,
  iohandle &out, const arr<double> &allthetas, bool realpart,
  bool &keys_written)
  {
  arr2<T> imgarr(tmm.size1(), tmm.size2());
  for (tsize mm=0; mm<imgarr.size1(); ++mm)
    for (tsize m=0; m<imgarr.size2(); ++m)
      imgarr[mm][m] = T (realpart ? tmm[mm][m].real() : tmm[mm][m].imag());

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

template <typename T> void checkZero (bool pol, const Alm<xcomplex<T> > &almT,
  const Alm<xcomplex<T> > &almG , const Alm<xcomplex<T> > &almC,
  arr<bool> &lzero)
  {
  lzero.allocAndFill (almT.Mmax()+1,true);
  for (int m=0; m<=almT.Mmax(); ++m)
    for (int l=m; l<=almT.Lmax(); ++l)
      {
      double sum = almT(l,m).norm();
      if (pol) sum += almG(l,m).norm()+almC(l,m).norm();
      if (sum>0) { lzero[m]=false; break; }
      }
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
  int lmax = params.find<int>("conv_lmax");
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

  int beammmax = params.find<int>("beammmax");
  bool pol = params.find<bool>("polarisation",true);
  double epsPow = params.find<double>("exponent",100.);
  safe_ptr<iohandle> out;
  vector<int32> present_sets;
  bool dp_output = params.find<bool>("double_precision_output",false);
  if (mpiMgr.master())
    {
    out = HandleManager.createObject (params.find<string>("ringset"),
      dp_output ? "ringset.LS_ringset_dp" : "ringset.LS_ringset");
    out->setKey ("beam_mmax",int32(beammmax));
    }
  bool keys_written=false;

  string infile_sky  = params.find<string>("sky_alm"),
         infile_beam = params.find<string>("beam_alm");
  Alm<xcomplex<T> > blmT,slmT,blmG,slmG,blmC,slmC;
  pol ? read_Alm_from_dmc (infile_sky, slmT, slmG, slmC, lmax, lmax)
      : read_Alm_from_dmc (infile_sky, slmT, lmax, lmax);
  pol ? read_Alm_from_dmc (infile_beam,blmT, blmG, blmC, lmax, beammmax)
      : read_Alm_from_dmc (infile_beam,blmT, lmax, beammmax);
  double fwhm = arcmin2rad*params.find<double>("fwhm_deconv",0.);
  pol ? smoothWithGauss (slmT, slmG, slmC, -fwhm)
      : smoothWithGauss (slmT, -fwhm);

  //Pre-calculation of sines and cosines used in rotating the output
  arr<xcomplex<double> > expmsky(nphi);
  for (int msky=-olmax; msky<nphi-olmax; ++msky)
    {
    double ang = -twopi*(((msky+olmax)*olmax)%nphi)/nphi;
    expmsky[msky+olmax].Set(cos(ang),sin(ang));
    }

  fix_arr<xcomplex<double>,4> rothalfpi;
  rothalfpi[0].Set( 1,0); rothalfpi[1].Set(0, 1);
  rothalfpi[2].Set(-1,0); rothalfpi[3].Set(0,-1);

  arr2<xcomplex<double> > imgarr(thetas.size(),nphi);

  arr<bool> blzero, slzero;
  checkZero (pol,blmT,blmG,blmC,blzero);
  checkZero (pol,slmT,slmG,slmC,slzero);

  arr<double> rthetas(pairs.size());
  for (tsize i=0; i<pairs.size(); ++i)
    rthetas[i]=thetas[pairs[i].r1];

  for(int mbeam=0; mbeam<=beammmax; ++mbeam)
    {
    if (blzero[mbeam]) continue; // all terms for this mbeam are zero

    imgarr.fill(0.);

#pragma omp parallel
{
    wignergen wgen(lmax,rthetas,1e-30);
    wigner_estimator estimator(lmax,epsPow);
#ifdef __SSE2__
    arr_align<V2df,16> val1(lmax+1), val2(lmax+1), val3(lmax+1), val4(lmax+1);
#else
    arr<double> val1(lmax+1), val2(lmax+1), val3(lmax+1), val4(lmax+1);
#endif

    int msky;
#pragma omp for schedule(dynamic,1)
    for(msky=-lmax; msky<=lmax; ++msky)
      {
      int absmsky = abs(msky);
      if (slzero[absmsky]) continue; // all terms for this msky are zero

      // s_l,-m = (-1)^m conj(s_lm)
      bool flipR=false, flipI=false;
      if (msky<0)
        {
        flipR = ((-msky)&1);
        flipI = !flipR;
        }
      // Pre-calculation of the sum over the product of b_lm and s_lm
      for (int l=max(absmsky,mbeam); l<=lmax; ++l)
        {
        xcomplex<T> sT=slmT(l,absmsky), bT=blmT(l,mbeam);
        double RR = sT.re*bT.re, II = sT.im*bT.im,
               IR = sT.im*bT.re, RI = sT.re*bT.im;
        if (pol)
          {
          xcomplex<T> sG=slmG(l,absmsky), bG=blmG(l,mbeam),
                      sC=slmC(l,absmsky), bC=blmC(l,mbeam);
          RR += sG.re*bG.re + sC.re*bC.re; II += sG.im*bG.im + sC.im*bC.im;
          IR += sG.im*bG.re + sC.im*bC.re; RI += sG.re*bG.im + sC.re*bC.im;
          }
        if (flipR) { RR=-RR; RI=-RI; }
        if (flipI) { IR=-IR; II=-II; }

#ifdef __SSE2__
        // s_lm*conj(b_lm')
        val1[l]=RR+II; val2[l]=IR-RI;
        // (-1)^(l+mbeam+msky) conj(s_lm*b_lm') (needed for ring pairs only)
        val3[l]=xpow(l+mbeam+msky,RR-II);
        val4[l]=xpow(l+mbeam+msky,-IR-RI);
#else
        // s_lm*conj(b_lm')
        val1[l]=RR+II; val2[l]=IR-RI;
        // (-1)^(l+mbeam+msky) conj(s_lm*b_lm') (needed for ring pairs only)
        val3[l]=xpow(l+mbeam+msky,RR-II); val4[l]=xpow(l+mbeam+msky,-IR-RI);
#endif
        }

      estimator.prepare_m (mbeam,msky);
      wgen.prepare(mbeam,msky);

#ifdef __SSE2__
      for (tsize ipair=0; ipair<pairs.size(); ipair+=2)
        {
        bool dual = ipair<(pairs.size()-1);
        double theta1=thetas[pairs[ipair].r1],
               theta2= dual ? thetas[pairs[ipair+1].r1] : theta1;
        if (estimator.canSkip(theta1) && estimator.canSkip(theta2))
          continue; // negligible dmm

        int firstl;
        const arr_align<V2df,16> &dmm = dual ?
          wgen.calc(ipair,ipair+1,firstl) : wgen.calc(ipair,ipair,firstl);

        V2df t1(0.), t2(0.), t3(0.), t4(0.);
        for (int l=firstl; l<=lmax; ++l)
          {
          t1+=val1[l]*dmm[l]; t2+=val2[l]*dmm[l];
          t3+=val3[l]*dmm[l]; t4+=val4[l]*dmm[l];
          }

        // sum_l d^l_mm'(theta) * s_lm * conj(b_lm')
        imgarr[pairs[ipair].r1][msky+olmax].Set(t1[0],t2[0]);
        if (pairs[ipair].r2>0) // we have a ring at (pi-theta1)
          imgarr[pairs[ipair].r2][-msky+olmax].Set(t3[0],t4[0]);
        if (dual)
          {
          imgarr[pairs[ipair+1].r1][msky+olmax].Set(t1[1],t2[1]);
          if (pairs[ipair+1].r2>0) // we have a ring at (pi-theta2)
            imgarr[pairs[ipair+1].r2][-msky+olmax].Set(t3[1],t4[1]);
          }
        }
#else
      for (tsize ipair=0; ipair<pairs.size(); ++ipair)
        {
        double theta1=thetas[pairs[ipair].r1];
        if (estimator.canSkip(theta1)) continue; // negligible dmm

        int firstl;
        const arr<double> &dmm=wgen.calc(ipair,firstl);
        double t1=0., t2=0., t3=0., t4=0.;

        if (pairs[ipair].r2>0) // pair of rings at theta1 and (pi-theta1)
          for (int l=firstl; l<=lmax; ++l)
            {
            t1 += val1[l]*dmm[l]; t2 += val2[l]*dmm[l];
            t3 += val3[l]*dmm[l]; t4 += val4[l]*dmm[l];
            }
        else // isolated ring at theta1
          for (int l=firstl; l<=lmax; ++l)
            { t1 += val1[l]*dmm[l]; t2 += val2[l]*dmm[l]; }

        // sum_l d^l_mm'(theta) * s_lm * conj(b_lm')
        imgarr[pairs[ipair].r1][msky+olmax].Set(t1,t2);
        if (pairs[ipair].r2>0) // we have a ring at (pi-theta1)
          imgarr[pairs[ipair].r2][-msky+olmax].Set(t3,t4);
        }
#endif
      }
} // end of parallel region

#pragma omp parallel
{
    cfft plan(nphi);
    int ith;
#pragma omp for schedule(static)
    for (ith=0; ith<int(imgarr.size1()); ++ith)
      {
      // rotate by -pi/2 in phi and psi direction for TC compatibility
      for (int msky=-lmax; msky<=lmax; ++msky)
        imgarr[ith][msky+olmax]*=rothalfpi[imodulo(-msky-mbeam,4)];
      plan.backward(&imgarr[ith][0]);
      for (int msky=-olmax; msky<nphi-olmax; ++msky)
        imgarr[ith][msky+olmax]*=expmsky[msky+olmax];
      }
} // end of parallel region

    if (mpiMgr.master())
      cout << "output: mbeam=" << mbeam << endl;
    if (mpiMgr.master())
      present_sets.push_back(int(mbeam));
    dp_output ? write_file<double> (imgarr,*out,allthetas,true,keys_written)
              : write_file<float>  (imgarr,*out,allthetas,true,keys_written);
    if (mbeam!=0)
      dp_output ? write_file<double> (imgarr,*out,allthetas,false,keys_written)
                : write_file<float>  (imgarr,*out,allthetas,false,keys_written);
    }

  if (mpiMgr.master())
    out->appendColumnRaw("ringsets_present",&present_sets[0],
      present_sets.size());
  }

} // unnamed namespace

int conviqt_v3_module (int argc, const char **argv)
  {
  module_startup ("conviqt_v3", argc, argv, mpiMgr.master());
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams(mpiMgr.master()));

  params.find<bool>("double_precision",false) ?
    convolve<double>(params) : convolve<float>(params);

  return 0;
  }
