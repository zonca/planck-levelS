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

/*! \file convolver_helper.h
 *  Utilities for total convolution applications
 *
 *  Copyright (C) 2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_CONVOLVER_HELPER_H
#define PLANCK_CONVOLVER_HELPER_H

#include <string>
#include "paramfile.h"
#include "lsconstants.h"
#include "xcomplex.h"
#include "alm.h"
#include "alm_powspec_tools.h"
#include "alm_dmcio.h"
#include "arr.h"
#include "sharp_cxx.h"

template<typename T> class convolver_helper
  {
  private:
    bool pol;
    int lmax, beammmax;
    Alm<xcomplex<T> > slmT, slmG, slmC, blmT, blmG, blmC;
    arr<double> lnorm;
    arr<bool> blzero;

    void checkZero (bool pol_, const Alm<xcomplex<T> > &almT,
      const Alm<xcomplex<T> > &almG , const Alm<xcomplex<T> > &almC,
      arr<bool> &lzero)
      {
      lzero.allocAndFill (almT.Mmax()+1,true);
      for (int m=0; m<=almT.Mmax(); ++m)
        for (int l=m; l<=almT.Lmax(); ++l)
          {
          double sum = almT(l,m).norm();
          if (pol_) sum += almG(l,m).norm()+almC(l,m).norm();
          if (sum>0) { lzero[m]=false; break; }
          }
      }

    void init(paramfile &params, const std::string &infile_beam)
      {
      using namespace std;
      lmax = params.find<int>("conv_lmax");
      beammmax = params.find<int>("beammmax");
      pol = params.find<bool>("polarisation",true);

      pol ? read_Alm_from_dmc (infile_beam,blmT, blmG, blmC, lmax, beammmax)
          : read_Alm_from_dmc (infile_beam,blmT, lmax, beammmax);

      if (params.param_present("sky_alm"))
        {
        string infile_sky  = params.find<string>("sky_alm");
        pol ? read_Alm_from_dmc (infile_sky, slmT, slmG, slmC, lmax, lmax)
            : read_Alm_from_dmc (infile_sky, slmT, lmax, lmax);
        double fwhm = arcmin2rad*params.find<double>("fwhm_deconv",0.);
        pol ? smoothWithGauss (slmT, slmG, slmC, -fwhm)
            : smoothWithGauss (slmT, -fwhm);
        }
      else
        {
        // special mode: compute sky a_lm of a single point source with
        // user-defined parameters
        slmT.Set(lmax,lmax);
        if (pol)
          { slmG.Set(lmax,lmax); slmC.Set(lmax,lmax); }

        sharp_cxxjob<T> job;
        double theta, phi0;
        T si=1, sq=0, su=0;
        int nph=1;
        ptrdiff_t ofs=0;
        // FIXME: still missing potential conversions deg->rad,
        //        latitude->colatitude, galactic->ecliptic
        // FIXME: still missing flux input, conversion to K_RJ or K_CMB
        theta=params.find<double>("pntsrc_theta");
        phi0=params.find<double>("pntsrc_phi");
        if (pol)
          {
          double vpol = si*params.find<double>("pntsrc_polfrac");
          double polang = params.find<double>("pntsrc_polangle");
          sq = vpol*cos(2*polang);
          su = vpol*sin(2*polang);
          }
        job.set_general_geometry (1, &nph, &ofs, &nph, &phi0, &theta, NULL);
        job.set_triangular_alm_info (lmax, lmax);
        job.map2alm(&si,&slmT(0,0).re,false);
        if (pol)
          job.map2alm_spin(&sq,&su,&slmG(0,0).re,&slmC(0,0).re,2,false);
        }
      checkZero(pol,blmT,blmG,blmC,blzero);

      lnorm.alloc(lmax+1);
      for (int i=0; i<=lmax; ++i)
        lnorm[i]=sqrt(4*pi/(2*i+1.));
      }
    void init(paramfile &params)
      { init (params, params.find<std::string>("beam_alm")); }


  public:
    convolver_helper(paramfile &params)
      { init(params); }
    convolver_helper(paramfile &params, const std::string &infile_beam)
      { init(params,infile_beam); }

    int Lmax () const
      { return lmax; }
    int beamMmax () const
      { return beammmax; }
    bool blZero (int k) const
      {
      planck_assert((k>=0) && (k<=beammmax), "bad value for k");
      return blzero[k];
      }

    void getAlm (int k, Alm<xcomplex<T> > &a1, Alm<xcomplex<T> > &a2) const
      {
      planck_assert((k>=0) && (k<=beammmax), "bad value for k");
      a1.Set(lmax,lmax);
      if (k!=0) a2.Set(lmax,lmax);
      double spinsign = (k==0) ? 1. : -1.;
      for (int m=0; m<=lmax; ++m)
        {
        double mfac=(m&1) ? -1.:1.;
        for (int l=m; l<=lmax; ++l)
          {
          if (l<k)
            a1(l,m)=a2(l,m)=0.;
          else
            {
            xcomplex<T> v1=slmT(l,m)*blmT(l,k),
                        v2=conj(slmT(l,m))*blmT(l,k)*mfac;
            if (pol)
              {
              v1+=slmG(l,m)*blmG(l,k)+slmC(l,m)*blmC(l,k);
              v2+=(conj(slmG(l,m))*blmG(l,k)+conj(slmC(l,m))*blmC(l,k))*mfac;
              }
            a1(l,m) = (v1+conj(v2)*mfac)*(0.5*spinsign*lnorm[l]);
            if (k>0)
              a2(l,m) = ((v1-conj(v2)*mfac).times_i())*(-spinsign*0.5*lnorm[l]);
            }
          }
        }
      }
  };

#endif
