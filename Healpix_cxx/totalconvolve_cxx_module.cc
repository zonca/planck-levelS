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
 *  This file implements the total convolution algorithm suggested by
 *  Wandelt & Gorski, Phys. Rev. D 63, 123002 (2001),
 *  with additions to include polarisation by
 *  Challinor et al., Phys. Rev. D 62, 123002 (2000)
 *
 *  It is based on the Fortran90 implementation by B. Wandelt
 *  with partial OpenMP parallelisation by S. Colombi.
 *
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm.h"
#include "iohandle_current.h"
#include "alm_dmcio.h"
#include "alm_powspec_tools.h"
#include "paramfile.h"
#include "fftpack_support.h"
#include "wigner.h"
#include "levels_facilities.h"
#include "lsconstants.h"
#include "announce.h"

using namespace std;

namespace {

void shrink_lmmax (const string &infile_sky, const string &infile_beam,
  int &lmax, int &beammmax, bool polarisation)
  {
  int olmax=lmax, obeammmax=beammmax;
  int tlmax, tmmax;
  polarisation ? get_almsize_dmc_pol (infile_sky, tlmax, tmmax)
               : get_almsize_dmc (infile_sky, tlmax, tmmax);
  lmax = min(lmax,tlmax);
  beammmax = min(beammmax,tmmax);
  polarisation ? get_almsize_dmc_pol (infile_beam, tlmax, tmmax)
               : get_almsize_dmc (infile_beam, tlmax, tmmax);
  lmax = min(lmax,tlmax);
  beammmax = min(beammmax,tmmax);
  if (lmax!=olmax)
    cout << "lmax reduced to " << lmax << endl;
  if (beammmax!=obeammmax)
    cout << "beammmax reduced to " << beammmax << endl;
  }

template<typename T> void write_file (const arr2<xcomplex<double> >&tmm,
  iohandle &out, bool realpart, bool &keys_set)
  {
// additional rows in theta-direction for higher-order interpolation
  tsize addrows = min(tmm.size1(),tsize(10));
  arr2<T> imgarr(tmm.size1()/2+1+2*addrows, tmm.size2());
  for (tsize mm=addrows; mm<imgarr.size1(); ++mm)
    for (tsize m=0; m<imgarr.size2(); ++m)
      imgarr[mm][m] = T(realpart ?
        tmm[mm-addrows][m].real() : tmm[mm-addrows][m].imag());
  for (tsize mm=0; mm<addrows; ++mm)
    for (tsize m=0; m<imgarr.size2(); ++m)
      imgarr[mm][m] = T (realpart ? tmm[tmm.size1()-addrows+mm][m].real()
                                  : tmm[tmm.size1()-addrows+mm][m].imag());

  if (!keys_set)
    {
    out.setKey("nphi",int32(imgarr.size2()));
    out.setKey("phi0",90.);
    out.setKey("dphi",360./imgarr.size2());
    out.setKey("ntheta",int32(imgarr.size1()));
    double dtheta = -360./tmm.size1();
    out.setKey("theta0",180.-addrows*dtheta);
    out.setKey("dtheta",dtheta);
    keys_set=true;
    }
  out.appendColumnRaw("ringsetdata",&(imgarr[0][0]),
    imgarr.size());
  }

template <typename T> void convolve (paramfile &params)
  {
  double fwhm = arcmin2rad*params.find<double>("fwhm_deconv",0.);
  bool polarisation = params.find<bool>("polarisation",true);

  string infile_beam = params.find<string>("beam_alm");
  string infile_sky = params.find<string>("sky_alm");
  string outfile = params.find<string>("ringset");

  int lmax = params.find<int>("conv_lmax");
  int lmax_out = params.find<int>("lmax_out",lmax);
  int beammmax = params.find<int>("beammmax");
// reduce lmax and beammmax if the input files allow it
  shrink_lmmax (infile_sky, infile_beam, lmax, beammmax, polarisation);

  int m_blocksize = params.find<int>("m_blocksize",beammmax+1);
  planck_assert(m_blocksize>0,"bad 'm_blocksize'");
  bool dp_output = params.find<bool>("double_precision_output",false);
  bool keys_set=false;

  Alm<xcomplex<T> > blmT,slmT,blmG,blmC,slmG,slmC;
  if (polarisation)
    {
    read_Alm_from_dmc (infile_sky, slmT, slmG, slmC, lmax, lmax);
    smoothWithGauss (slmT, slmG, slmC, -fwhm);
    read_Alm_from_dmc (infile_beam, blmT, blmG, blmC, lmax, lmax);
    }
  else
    {
    read_Alm_from_dmc (infile_sky, slmT, lmax, lmax);
    smoothWithGauss (slmT, -fwhm);
    read_Alm_from_dmc (infile_beam, blmT, lmax, lmax);
    }

  vector<int32> present_sets;
  safe_ptr<iohandle> out (HandleManager.createObject(outfile,
    dp_output ? "ringset.LS_ringset_dp" : "ringset.LS_ringset"));
  out->setKey ("beam_mmax",beammmax);

  int mmm_start=0;
  while (mmm_start<=beammmax)
    {
    int n_mmm=min(m_blocksize,beammmax+1-mmm_start);
    int mmm_stop=mmm_start+n_mmm-1;

    arr<arr2<xcomplex<double> > >tmmm(n_mmm);

{ // block for doing all pre-FFT calculations

    wigner_d_halfpi_risbo_openmp wigner(lmax);
    arr<double[4]> tmp2(lmax+1);

    cout << "starting convolution" << endl;
    for (int n=0; n<=lmax; ++n)
      {
      int flipn=(n%2==0) ? 1 : -1;
      const arr2<double> &d(wigner.recurse());

      for (int mmm=mmm_start; mmm<=min(n,mmm_stop); ++mmm)
        {
        int mmm0 = mmm-mmm_start;
        int flipmmm = (mmm%2==0) ? 1 : -1;
        xcomplex<double> tmp1T = xcomplex<double>(conj(blmT(n,mmm))),
                         tmp1G(0.,0.), tmp1C(0.,0.);
        double tmpnorm = norm(tmp1T);
        if (polarisation)
          {
          tmp1G = xcomplex<double>(conj(blmG(n,mmm)));
          tmp1C = xcomplex<double>(conj(blmC(n,mmm)));
          tmpnorm += norm(tmp1G)+norm(tmp1C);
          }

        if (tmpnorm>0)
          {
          if (tmmm[mmm0].size1()==0)
            {
            tmmm[mmm0].alloc(lmax+1,2*lmax+1);
#pragma omp parallel
{
            int mm;
#pragma omp for schedule(static)
            for (mm=0; mm<=lmax; ++mm)
              for (int m=0; m<=2*lmax; ++m)
                tmmm[mmm0][mm][m] = 0;
}
            }
          for (int m=0; m<=n; ++m)
            {
            double t1 = slmT(n,m).re*tmp1T.re;
            double t2 = slmT(n,m).re*tmp1T.im;
            double t3 =-slmT(n,m).im*tmp1T.im;
            double t4 = slmT(n,m).im*tmp1T.re;
            if (polarisation)
              {
              t1 += slmG(n,m).re*tmp1G.re + slmC(n,m).re*tmp1C.re;
              t2 += slmG(n,m).re*tmp1G.im + slmC(n,m).re*tmp1C.im;
              t3 -= slmG(n,m).im*tmp1G.im + slmC(n,m).im*tmp1C.im;
              t4 += slmG(n,m).im*tmp1G.re + slmC(n,m).im*tmp1C.re;
              }
            tmp2[m][0] = t1+t3;
            tmp2[m][1] = t2+t4;
            tmp2[m][2] = t1-t3;
            tmp2[m][3] = t2-t4;
            }

#pragma omp parallel
{
          int mm;
#pragma omp for schedule(static)
          for (mm=0; mm<=n; ++mm)
            {
            if (abs(d[n-mmm][n-mm])>1e-30)
              {
              int flipmm = (mm%2==0) ? 1 : -1;
              double tmpfact = (flipmmm*flipmm)*d[n-mmm][n-mm];

              tmmm[mmm0][mm][0].re += (flipn*d[n-mm][n])*tmpfact*tmp2[0][0];
              tmmm[mmm0][mm][0].im += (flipn*d[n-mm][n])*tmpfact*tmp2[0][1];
              int flip = -flipmm;

              xcomplex<double> *ptr = &tmmm[mmm0][mm][0];
              for (int m=1; m<=n; ++m)
                {
                if (abs(d[n-mm][n-m])>1e-30)
                  {
                  double fct=d[n-mm][n-m]*tmpfact;

                  if (flip<0) fct=-fct;
                  ptr[2*m-1].re += fct*tmp2[m][0];
                  ptr[2*m-1].im += fct*tmp2[m][1];

                  if (flipn!=flip) fct=-fct;
                  ptr[2*m].re += fct*tmp2[m][2];
                  ptr[2*m].im += fct*tmp2[m][3];
                  }
                flip=-flip;
                }
              }
            }
}
          }
        }
      }

} // end of block for doing all pre-FFT calculations

    arr2<xcomplex<double> > outarr(2*lmax_out+2,2*lmax_out+1);

    for (int mmm=mmm_start; mmm<=mmm_stop; ++mmm)
      {
      int mmm0=mmm-mmm_start;

      if (tmmm[mmm0].size1()!=0)
        {
// pad the array to the desired dimensions, apply symmetry condition
#pragma omp parallel
{
        int mm;
#pragma omp for schedule(static)
        for (mm=0; mm<=lmax_out+1; ++mm)
          for (int m=0; m<=lmax_out; ++m)
            {
            int xmm = 2*lmax_out+2-mm;
            int xm = 2*lmax_out+1-m;
            if ((m>lmax) || (mm>lmax))
              {
              outarr[mm][m] = 0;
              if (m>0) outarr[mm][xm] = 0;
              if (mm>0) outarr[xmm][m] = 0;
              if ((m>0)&&(mm>0)) outarr[xmm][xm] = 0;
              }
            else
              {
              int xxm1 = max(2*m-1,0);
              int xxm2 = 2*m;
              int flip = (((m+mmm)%2) != 0) ? -1 : 1;
              outarr[mm][m] = tmmm[mmm0][mm][xxm1];
              if (m>0) outarr[mm][xm] = tmmm[mmm0][mm][xxm2];
              if (mm>0) outarr[xmm][m] = tmmm[mmm0][mm][xxm1] * flip;
              if ((m>0)&&(mm>0)) outarr[xmm][xm] = tmmm[mmm0][mm][xxm2] * flip;
              }
            }
}

cout << "FFT phi" << endl;
#pragma omp parallel
{
        cfft p(2*lmax_out+1);

        int mm;
#pragma omp for schedule (static)
        for (mm=0; mm<=2*lmax_out+1; ++mm)
          p.backward (&(outarr[mm][0].re));
}

cout << "FFT theta" << endl;
#pragma omp parallel
{
        cfft p(2*lmax_out+2);
        arr<xcomplex<double> > tmparr(2*lmax_out+2);

        int m;
#pragma omp for schedule (static)
        for (m=0; m<=2*lmax_out; ++m)
          {
          for (int mm=0; mm<=2*lmax_out+1; ++mm) tmparr[mm]=outarr[mm][m];
          p.backward(tmparr);
          for (int mm=0; mm<=2*lmax_out+1; ++mm) outarr[mm][m]=tmparr[mm];
          }
}

cout << "output" << endl;
        present_sets.push_back(mmm);
        dp_output ? write_file<double> (outarr,*out,true,keys_set)
                  : write_file<float>  (outarr,*out,true,keys_set);
        if (mmm!=0)
          dp_output ? write_file<double> (outarr,*out,false,keys_set)
                    : write_file<float>  (outarr,*out,false,keys_set);
        }
      }

    mmm_start+=n_mmm;
    }

  out->appendColumnRaw ("ringsets_present",&present_sets[0],
    present_sets.size());
  }

} // unnamed namespace

int totalconvolve_cxx_module (int argc, const char **argv)
  {
  module_startup ("totalconvolve_cxx", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  params.find<bool>("double_precision",false) ?
    convolve<double>(params) : convolve<float>(params);

  return 0;
  }
