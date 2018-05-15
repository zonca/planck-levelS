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

#include "oofnoise.h"
#include "planck_rng.h"
#include "arr.h"
#include "focalplane_db.h"
#include "lsconstants.h"
#include "paramfile.h"

using namespace std;

class sde_process
  {
  private:
    double fknee, fsamp, sigma, alpha, beta, gamma, state;

    void recalc()
      {
      alpha = exp(-fknee/fsamp);
      beta = sigma*(1-alpha);
      gamma = 1/fknee;
      }

  public:
    sde_process (double sigma_, double fknee_, double fsamp_)
      : fknee(fknee_), fsamp(fsamp_), sigma(sigma_), state(0)
      { recalc(); }

    void reset()
      { state=0; }

    double get_next (double rand)
      {
      state = alpha*state + beta*rand;
      return state;
      }

    double power (double freq) const
      {
      double t1 = freq*gamma;
      return sigma*sigma/(1+t1*t1);
      }

    double Sigma() const
      { return sigma; }

    double Fknee() const
      { return fknee; }

    double Beta() const
      { return beta; }

    double Gamma() const
      { return gamma; }

    void update_sigma (double sigma_new)
      {
      state *= sigma_new/sigma;
      sigma=sigma_new;
      recalc();
      }
    void update_fsamp (double fsamp_new)
      {
      fsamp=fsamp_new;
      recalc();
      }
  };

class sde_set: public OofGenerator
  {
  private:
    planck_rng rng;
    double sigmawhite, f_knee, f_min, f_samp, slope;
    vector<sde_process> sde;
    bool do_white, do_oof;

    void warm_up (double f, int it)
      {
      for (unsigned int m=0; m<sde.size(); ++m)
        sde[m].update_fsamp(f);
      for (int i=0; i<it; ++i)
        nextSample();
      }

  public:
    sde_set (double sigmawhite_, double f_knee_, double f_min_, double f_samp_,
             double slope_, int seed, bool do_white_, bool do_oof_)
      : rng(seed), sigmawhite(sigmawhite_), f_knee(f_knee_), f_min(f_min_),
        f_samp(f_samp_), slope(slope_), do_white(do_white_), do_oof(do_oof_)
      {
      if (do_white)
        sde.push_back (sde_process (sigmawhite,1000*f_samp,f_samp));
      if (do_oof)
        {
        if (approx(slope,-2.0))
          {
          sde.push_back (sde_process (sigmawhite*pow(f_min/f_knee,
                         0.5*slope),f_min,f_samp));
          }
        else
          {
          int lbound = sde.size();
          double fhi=10*f_knee;
          double flo=0.1*f_min;
          int nproc = int(2*log10(fhi/flo));

          cout << "adding " << nproc << " processes" << endl;
          for (int m=0; m<nproc; ++m)
            {
            double freq = flo*pow(fhi/flo,double(m)/(nproc-1));
            double power = designpower(freq);
            double sigma = sqrt(0.5*power);
            sde.push_back (sde_process(sigma, freq, f_samp));
            }

          // iterate the sigmas until the correct total power is reached
          for (int q=0; q<50000; ++q)
            {
            for (unsigned int m=lbound; m<sde.size(); ++m)
              {
              double testpower = totalpower(sde[m].Fknee());
              double sigmagoal = sde[m].Sigma()
                               * sqrt(designpower(sde[m].Fknee())/testpower);
              sde[m].update_sigma(sde[m].Sigma()
                                  +0.01*(sigmagoal-sde[m].Sigma()));
              }
            }

          // "warm up" the coupled SDE processes by advancing them with
          // increasingly higher frequencies
          for (double freq=flo; freq<fhi; freq*=1.1)
            warm_up (freq,1000);
          warm_up (f_samp,1000);
          }
        }
      }

    virtual double nextSample()
      {
      double rand = rng.rand_gauss(), result=0;
      for (unsigned int m=0; m<sde.size(); ++m)
        result += sde[m].get_next(rand);
      return result;
      }

    virtual void addVec (arr<double> &result, int num_real)
      {
      for (int cnt=0; cnt<num_real; ++cnt)
        for (tsize n=0; n<result.size(); ++n)
          result[n] += nextSample();
      }

    virtual void reset()
      {
      for (unsigned int m=0; m<sde.size(); ++m)
        sde[m].reset();
      double fhi=10*f_knee;
      double flo=0.1*f_min;
      for (double freq=flo; freq<fhi; freq*=1.1)
        warm_up (freq,1000);
      warm_up (f_samp,1000);
      }

    double designpower (double freq) const
      {
      return sigmawhite*sigmawhite
        * (do_white
           + do_oof*pow(f_knee,-slope)/(pow(f_min,-slope)+pow(freq,-slope)));
      }

    double uncoupled_power(double freq) const
      {
      double result=0;
      for (unsigned int m=0; m<sde.size(); ++m)
        result += sde[m].power(freq);
      return result;
      }

    double crosspower(double freq) const
      {
      double result=0;
      double fsq=freq*freq;

      for (unsigned int i=0; i<sde.size(); ++i)
        for (unsigned int j=0; j<i; ++j)
          {
          double t1 = 1+fsq*sde[i].Gamma()*sde[j].Gamma();
          double t2 = sde[i].Gamma()-sde[j].Gamma();
          result += 2*sde[i].Sigma()*sde[j].Sigma()*t1
                  / (t1*t1 + fsq*t2*t2);
          }
      return result;
      }

    double totalpower(double freq) const
      { return uncoupled_power(freq) + crosspower(freq); }
  };


/*! A numeric filter which produces noise with the power spectrum

    P(f)=(1/fsample)^2*(f^2+fknee^2)/(f^2+fmin^2)

    when fed with Gaussian random numbers of sigma=1.
    \author Stephane Plaszczynski (plaszczy@lal.in2p3.fr) */
class oof2filter
  {
  private:
    double x1, y1, c0, c1, d0;

  public:
    oof2filter (double fmin, double fknee, double fsample)
      : x1(0), y1(0)
      {
      double w0 = pi*fmin/fsample, w1=pi*fknee/fsample;
      c0 = (1+w1)/(1+w0);
      c1 =-(1-w1)/(1+w0);
      d0 = (1-w0)/(1+w0);
      }

    void reset()
      { x1=y1=0; }

    double operator()(double x2)
      {
      double y2 = c0*x2 + c1*x1 + d0*y1;
      x1 = x2;
      y1 = y2;
      return y2;
      }
  };

/*! A numeric filter, based on superposition of 1/f2 filters.
    see : {Keshner,PROC-IEE,vol-70 (1982)}
    that approximates the power spectrum

    P(f)=(1/fsamp)^2[(f^2+fknee^2)/(f^2+fmin^2)]^(alpha/2)

    for 0<alpha<2, when fed with Gaussian random numbers of sigma=1.

    Errors should be below 1% for any alpha.

    \author Stephane Plaszczynski (plaszczy@lal.in2p3.fr) */
class oofafilter
  {
  private:
    vector<oof2filter> filter;

  public:
    oofafilter (double alpha, double fmin, double fknee, double fsample)
      {
      double lw0 = log10(twopi*fmin), lw1 = log10(twopi*fknee);

      int Nproc = max(1,int(2*(lw1-lw0)));
      double dp = (lw1-lw0)/Nproc;
      double p0 = lw0 + dp*0.5*(1+0.5*alpha);
      for (int i=0; i<Nproc; ++i)
        {
        double p_i = p0+i*dp;
        double z_i = p_i - 0.5*dp*alpha;

        filter.push_back
          (oof2filter(pow(10.,p_i)/twopi,pow(10.,z_i)/twopi,fsample));
        }
      }

    double operator()(double x2)
      {
      for (unsigned int i=0; i<filter.size(); ++i)
        x2 = filter[i](x2);
      return x2;
      }

    void reset()
      {
      for (unsigned int i=0; i<filter.size(); ++i)
        filter[i].reset();
      }
  };

class Oof2Noise: public OofGenerator
  {
  private:
    planck_rng rng;
    oof2filter filter;
    double sigma_;

  public:
    Oof2Noise(double sigma, double fmin, double fknee, double fsample,
      int seed)
      : rng(seed), filter(fmin, fknee, fsample), sigma_(sigma)
      {}

    virtual double nextSample()
      {
      return sigma_*filter(rng.rand_gauss());
      }

    virtual void addVec (arr<double> &result, int num_real)
      {
      for (int cnt=0; cnt<num_real; ++cnt)
        for (tsize n=0; n<result.size(); ++n)
          result[n] += nextSample();
      }

    virtual void reset()
      {
      filter.reset();
      }
  };

class OofaNoise : public OofGenerator
  {
  private:
    planck_rng rng;
    oofafilter filter;
    double sigma;

  public:
    OofaNoise(double sigmawhite_, double f_knee_, double f_min_,
      double f_samp_, double slope_, int seed)
      : rng(seed), filter(slope_, f_min_, f_knee_, f_samp_), sigma(sigmawhite_)
      {}

    virtual double nextSample()
      {
      return sigma*filter(rng.rand_gauss());
      }

    virtual void addVec (arr<double> &result, int num_real)
      {
      for (int cnt=0; cnt<num_real; ++cnt)
        for (tsize n=0; n<result.size(); ++n)
          result[n] += nextSample();
      }

    virtual void reset()
      {
      filter.reset();
      }
  };

class sde: public OofGenerator
  {
  private:
    arr<double> y0, ttau, weight;
    int nprocesses;
    double sigma;
    bool whiteQ;
    arr<planck_rng> rng;
    planck_rng rngwhite;

  public:
    sde (int seed, double sigma_in, double fknee, double slope,
         double f0, bool whiteQ_in, bool oofQ)
      : sigma(sigma_in), whiteQ(whiteQ_in)
      {
      const double tau1=1.7, tau2=1./f0;
      if (oofQ)
//        if (approx(slope,-2))
//          nprocesses = 1;
//        else
          nprocesses = int(2*log10(tau2/tau1));
      else
        nprocesses = 0;
      cout << "nprocesses = " << nprocesses << endl;
      arr<double> tau(nprocesses);
      ttau.alloc(nprocesses);
      y0.alloc(nprocesses);
      weight.alloc(nprocesses);
      rng.alloc(nprocesses);
      rngwhite.seed(seed+nprocesses);
      for (int i=0; i<nprocesses; ++i)
        {
        rng[i].seed(seed+i);
        double tmp = (nprocesses==1) ? 1 : (double(i)/(nprocesses-1));
        tau[i] = tau1* pow(tau2/tau1,tmp);
//        ttau[i] = exp(-1/tau[i]);
        ttau[i] = 1-1/tau[i];
        y0[i] = rng[i].rand_gauss()*sqrt(tau[i]/(2-1/tau[i]));
        weight[i] = pow(tau[i],(-0.5*slope)-1.);
        }
      double pps = 0;
      for (int i=0; i<nprocesses; ++i)
        pps += weight[i]*weight[i]*tau[i]*tau[i]/(1+fknee*fknee*tau[i]*tau[i]);
      double tmp1 = sigma/sqrt(pps);
      for (int i=0; i<nprocesses; ++i) weight[i] *= tmp1;
      }

    virtual double nextSample()
      {
      double result= whiteQ ? sigma*rngwhite.rand_gauss() : 0;
      for (int i=0; i<nprocesses; ++i)
        {
        y0[i] = y0[i]*ttau[i] + rng[i].rand_gauss();
        result += y0[i]*weight[i];
        }
      return result;
      }

#ifdef _OPENMP

    virtual void addVec (arr<double> &result, int num_real)
      {
      arr2<double> tmp(nprocesses,result.size());
#pragma omp parallel
{
      int i;
#pragma omp for schedule(static,1)
      for (i=0; i<=nprocesses; ++i)
        {
        if (i<nprocesses)
          {
          double y0l = y0[i];
          planck_rng rngl=rng[i];
          for (tsize n=0; n<result.size(); ++n) tmp[i][n]=0;
          for (int cnt=0; cnt<num_real; ++cnt)
            {
            for (tsize n=0; n<result.size(); ++n)
              {
              y0l = y0l*ttau[i] + rngl.rand_gauss();
              tmp[i][n] += y0l*weight[i];
              }
            }
          y0[i]=y0l;
          rng[i]=rngl;
          }
        else
          {
          if (whiteQ)
            for (int cnt=0; cnt<num_real; ++cnt)
              {
              for (tsize n=0; n<result.size(); ++n)
                result[n] += sigma*rngwhite.rand_gauss();
              }
          }
        }
} // end of parallel region

#pragma omp parallel
{
      int n, sz=result.size();
#pragma omp for schedule(static)
      for (n=0; n<sz; ++n)
        {
        double t1 = 0;
        for (int i=0; i<nprocesses; ++i) t1+=tmp[i][n];
        result[n] += t1;
        }
} //end of parallel region
      }

#else

    virtual void addVec (arr<double> &result, int num_real)
      {
      for (int cnt=0; cnt<num_real; ++cnt)
        {
        for (tsize n=0; n<result.size(); ++n)
          {
          for (int i=0; i<nprocesses; ++i)
            {
            y0[i] = y0[i]*ttau[i] + rng[i].rand_gauss();
            result[n] += y0[i]*weight[i];
            }
          if (whiteQ)
            result[n] += sigma*rngwhite.rand_gauss();
          }
        }
      }

#endif

    virtual void reset()
      {
      for (int i=0; i<nprocesses; ++i)
        {
        double tau = 1./(1-ttau[i]);
        y0[i] = rng[i].rand_gauss()*sqrt(tau/(2-1/tau));
        }
     }
  };

Oofnoise::Oofnoise (focalplane_db &fpdb, const string &det_id,
  int num_realisations, paramfile &params)
  : num_real(num_realisations),
    inv_num_real(1./num_realisations), curperiod(0)
  {
  int noise_seed = params.find<int>("oof_rand_seed");

  string mode = params.find<string>("oof_mode","BOTH");
  bool white=true, oof=true;
  if (mode=="WHITE") oof=false;
  else if (mode=="OOF") white=false;
  else planck_assert (mode=="BOTH", "unknown OOF mode");

  double f_samp = fpdb.getValue<double>(det_id,"f_samp");
  double f_knee = fpdb.getValue<double>(det_id,"f_knee");
  double net_rj = fpdb.getValue<double>(det_id,"net_rj");
  double alpha = fpdb.getValue<double>(det_id,"alpha");
  double f_min = fpdb.getValue<double>(det_id,"f_min");
  double epsilon = fpdb.getValue<double>(det_id,"epsilon");

  double sigma = 0.5*(1+epsilon)*net_rj*sqrt(f_samp);
  cout << "calculated noise rms: " << sigma << endl;

  string oofmethod = params.find<string>("oof_method","CLASSIC");
  stationary_periods = params.find<int>("oof_stationary_periods",-1);

  if (oofmethod=="CLASSIC")
    {
    double fknp = twopi*f_knee/f_samp;
    double f0 = twopi*f_min/f_samp;
    generator = new sde(noise_seed,sigma,fknp,-alpha,f0,white,oof);
    }
  else if (oofmethod=="NEWNOISE")
    generator = new sde_set (sigma, twopi*f_knee, twopi*f_min,
                             f_samp, -alpha, noise_seed, white, oof);
  else if (oofmethod=="OOF2NOISE")
    {
    planck_assert(approx(alpha,2.),
      "Oof2Noise can only be used with a slope of -2");
    planck_assert(mode=="BOTH",
      "Oof2Noise can only be used when oof_mode==BOTH");
    generator = new Oof2Noise (sigma, f_min, f_knee, f_samp, noise_seed);
    }
  else if (oofmethod=="OOFANOISE")
    {
    planck_assert(mode=="BOTH",
      "OofaNoise can only be used when oof_mode==BOTH");
    generator = new OofaNoise (sigma, f_knee, f_min, f_samp, -alpha,
                               noise_seed);
    }
  else
    planck_fail ("unknown oof_method '" + oofmethod + "' specified");
  }

void Oofnoise::Get_Noise (int /*period*/, int num_samp, arr<double> &noise)
  {
  noise.alloc(num_samp);
  noise.fill(0);

  if (++curperiod==stationary_periods)
    {
    generator->reset();
    curperiod=0;
    }

  generator->addVec (noise, num_real);
  for (int i=0; i<num_samp; ++i) noise[i]*=inv_num_real;
  }

void Oofnoise::Skip (int nskip)
  {
  for (int i=0; i<nskip; ++i) generator->nextSample ();
  }
