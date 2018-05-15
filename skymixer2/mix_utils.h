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

#ifndef PLANCK_MIX_UTILS_H
#define PLANCK_MIX_UTILS_H

#include "lsconstants.h"
#include "arr.h"
#include <algorithm>
#include <cmath>

class Spectrum
  {
  public:
    virtual double value(double f) const = 0;
    /// Returns the average spectrum value between \a flo and \a fhi.
    virtual double value_avg(double flo, double fhi) const
      { return value (0.5*(flo+fhi)); }
    virtual ~Spectrum() {}
  };

class Planck_Spectrum: public Spectrum
  {
  protected:
    double T;

  public:
    Planck_Spectrum (double temp)
      : T(temp) {}

    virtual double value (double f) const
      {
      using namespace std;
      double f2=2.0*(hPlanck*f)*(f/speedOfLight)*(f/speedOfLight);
      double e=(hPlanck*f)/(kBoltzmann*T);
      return f2/(exp(e)-1.0);
      }
  };

class Power_Spectrum: public Spectrum
  {
  private:
    double exponent, ref_freq;
//    double nu_0, ampl;

  public:
    Power_Spectrum (double expo, double fref)
      : exponent(expo), ref_freq(fref) {}

    virtual double value (double f) const
      {
      using namespace std;
      return pow(f/ref_freq,exponent);
//      return ampl*pow((f-nu_0)/ref_freq,exponent);
      }
  };

class Dust_Spectrum: public Planck_Spectrum
  {
  private:
    double a, ref;

  public:
    Dust_Spectrum (double temp, double f0, double alpha=2)
      : Planck_Spectrum(temp),
        a(alpha),
        ref(std::pow(f0,alpha)*Planck_Spectrum::value(f0))
      {}

    virtual double value (double f) const
      {
      using namespace std;
      return pow(f,a)*Planck_Spectrum::value(f)/ref;
      }
  };

class Dust2_Spectrum: public Spectrum
  {
  private:
    double a1, a2, frac1, frac2, q1q2, t1, t2, f_ref;
    Planck_Spectrum spec1, spec2;

  public:
    Dust2_Spectrum ()
      : a1(1.67), a2(2.70), frac1(0.0363), frac2(1-frac1), q1q2(13),
        t1(9.4), t2(16.2), f_ref(speedOfLight/1e-4), spec1(t1), spec2(t2)
      {}

    virtual double value (double f) const
      {
      using namespace std;
      return   (  frac1*q1q2*pow(f/f_ref,a1)*spec1.value(f)
                + frac2*     pow(f/f_ref,a2)*spec2.value(f))
             / (  frac1*q1q2*                spec1.value(f_ref)
                + frac2*                     spec2.value(f_ref));
      }
  };

class SZ_Spectrum: public Planck_Spectrum
  {
  public:
    SZ_Spectrum (double temp)
      : Planck_Spectrum(temp) {}

    virtual double value (double f) const
      {
      using namespace std;
      double x=hPlanck*f/(kBoltzmann*T);
      double e=exp(x);
      double d=x*e/(e-1)*(x*(e+1)/(e-1)-4.0);
      return d*Planck_Spectrum::value(f);
      }
  };

class SZkin_Spectrum: public Planck_Spectrum
  {
  public:
    SZkin_Spectrum (double temp)
      : Planck_Spectrum(temp) {}

    virtual double value (double f) const
      {
      using namespace std;
      double x=hPlanck*f/(kBoltzmann*T);
      double e=exp(x);
      double d=x*e/(e-1.0);
      return -d*Planck_Spectrum::value(f);
      }
  };

class FreeFree_Spectrum: public Spectrum
  {
  private:
    double f1, f2;

  public:
    FreeFree_Spectrum (double T=10000)
      {
      using namespace std;
      double e=sqrt(2.30708E-19); // electron charge in cgs
      double m_e=9.1094E-28;  // electron mass in g
      double k = 1.38066E-16; // Boltzmann constant in erg/K
      f1 = pow(2*k*T,1.5)/(pi*e*e*sqrt(m_e));
      double Tnorm = T/1e4;
      f2 = 14*pow(Tnorm,0.517)*pow(10,0.029/Tnorm);
      }

    virtual double value (double f) const
      {
      using namespace std;
      double gamma = 0.577215664901532; // Euler's gamma constant
      double g = sqrt(3.)/pi * (log(f1/f)-2.5*gamma);
      double Tbright = 1e-6*f2*g/(f*f*1e-20);
      return 2*f*f*kBoltzmann*Tbright/(speedOfLight*speedOfLight);
      }
  };

class CO_Spectrum: public Spectrum
  {
  private:
    arr<double> strength;

  public:
    CO_Spectrum (double T=20)
      : strength(7)
      {
      using namespace std;
      double tmp1 = exp(-(2*1.15e11*hPlanck/(kBoltzmann*T)));
      strength[0] = 1;
      for (int j=1; j<7; ++j)
        {
        double q = (2*j+3)/(2*j+1)*tmp1;
        strength[j] = q*strength[j-1];
        }
      for (int j=0; j<7; ++j)
        {
        double freq = (j+1)*1.15e11;
        strength[j] *= 2*kBoltzmann*freq*freq/(speedOfLight*speedOfLight);
        }
      }

    virtual double value (double) const
      {
      planck_fail ("CO_Spectrum::value() should not be called");
      }

    virtual double value_avg (double flo, double fhi) const
      {
      double res = 0;
      for (int j=0; j<7; ++j)
        {
        double freq = (j+1)*1.15e11;
        if ((freq>=flo) && (freq<fhi))
          res += strength[j]/(fhi-flo);
        }
      return res;
      }
  };


class Response_function
  {
  public:
    virtual double convolve (const Spectrum &spec) const = 0;
    virtual double area () const = 0;
    virtual ~Response_function() {}
  };

class Tophat: public Response_function
  {
  private:
    double flo, fhi, height;

  public:
    Tophat (double f1, double f2, double h=1)
      : flo(f1), fhi(f2), height(h) {}

    virtual double convolve (const Spectrum &spec) const
      {
      const int nsteps=10000;
      double conv=0;
      double deltaf=(fhi-flo)/nsteps;
      for (int m=0; m<nsteps; ++m)
        {
        double f1 = flo + m*deltaf;
        conv += spec.value_avg(f1,f1+deltaf);
        }
      return conv*deltaf*height;
      }
    virtual double area () const
      { return (fhi-flo)*height; }
  };

class Delta: public Response_function
  {
  private:
    double freq,integ;

  public:
    Delta (double f, double i=1)
      : freq(f), integ(i) {}

    virtual double convolve (const Spectrum &spec) const
      { return integ*spec.value_avg(freq*0.99999,freq*1.00001); }
    virtual double area () const
      { return integ; }
  };

class Gauss: public Response_function
  {
  private:
    double mean, sigma;

  public:
    Gauss (double Mean, double Sigma)
      : mean(Mean), sigma(Sigma) {}

    virtual double convolve (const Spectrum &spec) const
      {
      using namespace std;
      const int nsteps=10000;
      double conv=0;
      double flo=max(0.,mean-30*sigma),fhi=mean+30*sigma;
      double deltaf=(fhi-flo)/nsteps;
      for (int m=0; m<nsteps; ++m)
        {
        double fmid = flo + (m+.5)*deltaf;
        double f1 = flo + m*deltaf;
        conv += spec.value_avg(f1,f1+deltaf)*exp(-pow((fmid-mean)/sigma,2));
        }
      return conv*deltaf;
      }
    virtual double area () const
      { return std::sqrt(pi)*sigma; }
  };

#endif
