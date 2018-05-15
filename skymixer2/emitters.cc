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
 *  Copyright (C) 2007-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

class Emitter
  {
  public:
//  Emitter (paramfile &params);
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &mapQ,
                               Healpix_Map<float> &mapU) const = 0;
    virtual ~Emitter() {}
  };

class CmbEmitter_PRSM: public Emitter
  {
  private:
    string mapname;

  public:
    CmbEmitter_PRSM (paramfile &params)
      {
      mapname = params.find<string> ("CMB_PRSM_map");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &mapQ,
                               Healpix_Map<float> &mapU) const
      {
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor = 0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double fmid = 0.5*(freq[m]+freq[m+1]);
        double x=(fmid*hPlanck)/(kBoltzmann*tcmb);
        double exm1 = exp(x)-1;
        double nu_c = fmid/speedOfLight;
        double conv = 1e-6*x*x*exp(x)/(exm1*exm1)
                      *2e20*kBoltzmann*nu_c*nu_c;
        double dflux = wgt[m]*conv;
#ifdef CALC_FLUX
        factor += dflux;
#else
        factor += dflux*antennaFact(fmid);
#endif
        }
      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (mapname,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor*map[m]);
      map_import (mapname,map,"Q_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapQ[m] += float(factor*map[m]);
      map_import (mapname,map,"U_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapU[m] += float(factor*map[m]);
      }
  };

class DustEmitter_PRSM: public Emitter
  {
  private:
    string sfdmap, cosmap, sinmap;

  public:
    DustEmitter_PRSM (paramfile &params)
      {
      sfdmap = params.find<string> ("DUST_PRSM_sfdmap");
      cosmap = params.find<string> ("DUST_PRSM_cosmap");
      sinmap = params.find<string> ("DUST_PRSM_sinmap");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &mapQ,
                               Healpix_Map<float> &mapU) const
      {
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor = 0;
      Dust2_Spectrum spec;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double conv = spec.value_avg(freq[m],freq[m+1]);
        double dflux = wgt[m]*conv;
#ifdef CALC_FLUX
        factor += dflux;
#else
        factor += dflux*antennaFact(0.5*(freq[m]+freq[m+1]));
#endif
        }
      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (sfdmap,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor*map[m]);
      const double p=0.05;
      Healpix_Map<float> map2(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (cosmap,map2,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapQ[m] += float(p*factor*map[m]*map2[m]);
      map_import (sinmap,map2,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapU[m] += float(p*factor*map[m]*map2[m]);
      }
  };

class SynchroEmitter_PRSM: public Emitter
  {
  private:
    string n_map, n_idxmap, sinmap, cosmap;

  public:
    SynchroEmitter_PRSM (paramfile &params)
      {
      n_map = params.find<string> ("SYNCHRO_PRSM_map");
      n_idxmap = params.find<string> ("SYNCHRO_PRSM_idxmap");
      cosmap = params.find<string> ("SYNCHRO_PRSM_cosmap");
      sinmap = params.find<string> ("SYNCHRO_PRSM_sinmap");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &mapQ,
                               Healpix_Map<float> &mapU) const
      {
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double nu0=0.408e9; // Hz
      Healpix_Map<float> map(mapT.Nside(), mapT.Scheme(), SET_NSIDE),
                         idxmap(mapT.Nside(), mapT.Scheme(), SET_NSIDE);
      map_import(n_map,map,"I_Stokes");
      map_import (n_idxmap,idxmap,"I_Stokes");
      arr<double> midfreq(wgt.size());
      for (tsize i=0; i<wgt.size(); ++i)
        {
        midfreq[i]=0.5*(freq[i]+freq[i+1]);
#ifdef CALC_FLUX
        wgt[i]/=antennaFact(midfreq[i]);
#endif
        midfreq[i]/=nu0;
        }
      Healpix_Map<float> tmap(mapT.Nside(), mapT.Scheme(), SET_NSIDE);
      tmap.fill(0);
      for (int m=0; m<tmap.Npix(); ++m)
        {
        if ((m%100000)==0) cout << m << " pixels" << endl;
        for (tsize i=0; i<wgt.size(); ++i)
          {
          tmap[m] += float(map[m]*pow(midfreq[i],double(-idxmap[m]))*wgt[i]);
          }
        }
      for (int m=0; m<mapT.Npix(); ++m)
        mapT[m] += tmap[m];

      map_import (cosmap,map,"I_Stokes");
      for (int m=0; m<mapT.Npix(); ++m)
        {
        double frac=(3.*abs(idxmap[m])-3.)/(3.*abs(idxmap[m])-1.);
        mapQ[m] += float(frac*tmap[m]*map[m]);
        }
      map_import (sinmap,map,"I_Stokes");
      for (int m=0; m<mapT.Npix(); ++m)
        {
        double frac=(3.*abs(idxmap[m])-3.)/(3.*abs(idxmap[m])-1.);
        mapU[m] += float(frac*tmap[m]*map[m]);
        }
      }
  };

class SZEmitter_PRSM: public Emitter
  {
  private:
    string szmap;

    static double specfunc(double f)
      {
      double x = hPlanck*f/(kBoltzmann*tcmb);
      return x*x*x*x*exp(x)/pow(exp(x)-1,2)*(x/tanh(0.5*x)-4);
      }

  public:
    SZEmitter_PRSM (paramfile &params)
      {
      szmap = params.find<string> ("SZ_PRSM_map");
      }
   virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &,
                               Healpix_Map<float> &) const
      {
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor = 0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double fmid = 0.5*(freq[m]+freq[m+1]);
        double dflux = wgt[m]*specfunc(fmid);
#ifdef CALC_FLUX
        factor += dflux;
#else
        factor += dflux*antennaFact(fmid);
#endif
        }
      factor *= 1e-6/(4*pi)*12.*2048.*2048.; // Jy -> MJy/sr
      factor /= specfunc(4e11);
      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (szmap,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor*map[m]);
      }
  };

class RadioEmitter_PRSM: public Emitter
  {
  private:
    string flat, f_cos, f_sin, steep, s_cos, s_sin;

  public:
    RadioEmitter_PRSM (paramfile &params)
      {
      flat = params.find<string> ("RADIO_PRSM_flatmap");
      f_cos = params.find<string> ("RADIO_PRSM_flatcos");
      f_sin = params.find<string> ("RADIO_PRSM_flatsin");
      steep = params.find<string> ("RADIO_PRSM_steepmap");
      s_cos = params.find<string> ("RADIO_PRSM_steepcos");
      s_sin = params.find<string> ("RADIO_PRSM_steepsin");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &mapQ,
                               Healpix_Map<float> &mapU) const
      {
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor_flat=0, factor_steep=0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double fmid = 0.5*(freq[m]+freq[m+1]);
#ifdef CALC_FLUX
        factor_flat += wgt[m];
        factor_steep += pow(fmid/3e10,0.8)*wgt[m];
#else
        factor_flat += wgt[m]*antennaFact(fmid);
        factor_steep += pow(fmid/3e10,0.8)*wgt[m]*antennaFact(fmid);
#endif
        }
      factor_flat *= 1e-6/(4*pi)*12.*2048.*2048.; // Jy -> MJy/sr
      factor_steep *= 1e-6/(4*pi)*12.*2048.*2048.; // Jy -> MJy/sr

      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      Healpix_Map<float> map2(mapT.Nside(),mapT.Scheme(),SET_NSIDE);

      map_import (flat,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor_flat*map[m]);
      const double p_flat=0.027;
      map_import (f_cos,map2,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapQ[m] += float(factor_flat*p_flat*map[m]*map2[m]);
      map_import (f_sin,map2,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapU[m] += float(factor_flat*p_flat*map[m]*map2[m]);

      map_import (steep,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor_steep*map[m]);
      const double p_steep=0.048;
      map_import (s_cos,map2,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapQ[m] += float(factor_steep*p_steep*map[m]*map2[m]);
      map_import (s_sin,map2,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapU[m] += float(factor_steep*p_steep*map[m]*map2[m]);
      }
  };

class InfraredEmitter_PRSM: public Emitter
  {
  private:
    string path;
    enum { nmaps=6 };

  public:
    InfraredEmitter_PRSM (paramfile &params, const string &suffix)
      {
      path = params.find<string> ("data_path"+suffix);
      }

    void getmaps (Healpix_Map<float> &maplo, Healpix_Map<float> &maphi,
      double &flo, double &fhi, double f) const
      {
      const double mapfreq[] = { 100e9, 143e9, 217e9, 353e9, 545e9, 857e9 };
      const char *mapname[] = { "all_sky_PLANCK_3000mic.hpx2048.fits",
                                "all_sky_PLANCK_2100mic.hpx2048.fits",
                                "all_sky_PLANCK_1380mic.hpx2048.fits",
                                "all_sky_PLANCK_850mic.hpx2048.fits",
                                "all_sky_PLANCK_550mic.hpx2048.fits",
                                "all_sky_PLANCK_350mic.hpx2048.fits" };

      int mres=0;
      for (int m=1; m<nmaps-1; ++m)
        if (f>mapfreq[m]) mres=m;

      if (approx(flo,mapfreq[mres]))
        return;

      map_import(path+mapname[mres],maplo,"I_Stokes");
      map_import(path+mapname[mres+1],maphi,"I_Stokes");
      flo = mapfreq[mres];
      fhi = mapfreq[mres+1];
      }

    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &mapQ,
                               Healpix_Map<float> &mapU) const
      {
      const double p = 0.01;
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      Healpix_Map<float> maplo(mapT.Nside(), mapT.Scheme(), SET_NSIDE),
                         maphi(mapT.Nside(), mapT.Scheme(), SET_NSIDE),
                         cmap (mapT.Nside(), mapT.Scheme(), SET_NSIDE),
                         smap (mapT.Nside(), mapT.Scheme(), SET_NSIDE);

      map_import(path+"/Cos_2theta_2048_random.fits",cmap,"I_Stokes");
      map_import(path+"/Sin_2theta_2048_random.fits",smap,"I_Stokes");

      arr<double> midfreq(wgt.size());
      for (tsize i=0; i<wgt.size(); ++i)
        {
        midfreq[i]=0.5*(freq[i]+freq[i+1]);
#ifdef CALC_FLUX
        wgt[i]*=1e-6;
#else
        wgt[i]*=antennaFact(midfreq[i])*1e-6;
#endif
        }
      double flo=-1, fhi=-1;
      for (tsize i=0; i<wgt.size(); ++i)
        {
        getmaps (maplo,maphi,flo,fhi,midfreq[i]);
        double fact = 1./log(fhi/flo);
        for (int m=0; m<mapT.Npix(); ++m)
          {
          if ((maplo[m]>0) && (maphi[m]>0))
            {
            double beta = log(double(maphi[m]/maplo[m]))*fact;
            double val = wgt[i]*maplo[m]*pow(midfreq[i]/flo,beta);
            mapT[m] += float(val);
            mapQ[m] += float(p*val*cmap[m]);
            mapU[m] += float(p*val*smap[m]);
            }
          }
        }
      }
 };

class FreeFreeEmitter_PRSM: public Emitter
  {
  private:
    string ha_name, ebv_name;

    static double r2ff (double freq)
      {
      const double T4 = 0.7, Z=1.0;
      double gaunt_correct = 3.96*pow(T4,0.21)*pow(freq*1e-9/40,-.14);
      return 14 * pow(T4,.517)*pow(10.,.029/T4)*1.08*Z*Z*gaunt_correct
             /(freq*freq*1e-20);
      }

  public:
    FreeFreeEmitter_PRSM (paramfile &params)
      {
      ha_name = params.find<string> ("FREEFREE_PRSM_hamap");
      ebv_name = params.find<string> ("FREEFREE_PRSM_ebvmap");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &,
                               Healpix_Map<float> &) const
      {
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      Healpix_Map<float> ha(mapT.Nside(), mapT.Scheme(), SET_NSIDE),
                         dust(mapT.Nside(), mapT.Scheme(), SET_NSIDE);
      map_import(ha_name,ha,"I_Stokes");
      map_import (ebv_name,dust,"I_Stokes");

      double factor=0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double fmid = 0.5*(freq[m]+freq[m+1]);
        double dTan = wgt[m]*r2ff(fmid)*1e-6;
#ifdef CALC_FLUX
        factor += dTan/antennaFact(fmid);
#else
        factor += dTan;
#endif
        }
      const double scale=2.51, limit=1, f_d=0.33;
      for (int m=0; m<mapT.Npix(); ++m)
        {
        if (!(approx<double>(ha[m],Healpix_undef)
            ||approx<double>(dust[m],Healpix_undef)))
          {
          double abs_correct = min(limit,f_d*scale*dust[m]);
          double abst = pow(10.,abs_correct*0.4);
          double halpha_correct = ha[m]*abst;
          mapT[m] += float(halpha_correct*factor);
          }
        }
      }
  };

class CmbEmitter: public Emitter
  {
  private:
    string mapname;
    bool monopole;

  public:
    CmbEmitter (paramfile &params)
      {
      mapname = params.find<string> ("CMB_LS_map");
      monopole = params.find<bool> ("CMB_LS_monopole", false);
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &mapQ,
                               Healpix_Map<float> &mapU) const
      {
      Planck_Spectrum spec(tcmb);
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor = 0;
      double offset = 0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double fmid = 0.5*(freq[m]+freq[m+1]);
        double specval = spec.value_avg(freq[m],freq[m+1]);
        double x=(fmid*hPlanck)/(kBoltzmann*tcmb);
        double prefactor = 1e20*x/((1.0-exp(-x))*tcmb);
        double dflux = wgt[m]*prefactor*specval;
#ifdef CALC_FLUX
        factor += dflux;
        offset += wgt[m]*1e20*specval;
#else
        factor += dflux*antennaFact(fmid);
        offset += wgt[m]*1e20*specval*antennaFact(fmid);
#endif
        }
      if (!monopole) offset=0;
      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (mapname,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(offset+factor*map[m]);
      map_import (mapname,map,"Q_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapQ[m] += float(factor*map[m]);
      map_import (mapname,map,"U_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapU[m] += float(factor*map[m]);
      }
  };

class SZEmitter: public Emitter
  {
  private:
    string mapname;

  public:
    SZEmitter (paramfile &params)
      {
      mapname = params.find<string> ("SZ_LS_map");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &,
                               Healpix_Map<float> &) const
      {
      SZ_Spectrum spec(tcmb);
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor = 0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double specval = spec.value_avg(freq[m],freq[m+1]);
        double dflux = wgt[m]*1e20*specval;
#ifdef CALC_FLUX
        factor += dflux;
#else
        factor += dflux*antennaFact(0.5*(freq[m]+freq[m+1]));
#endif
        }
      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (mapname,map,"I_Stokes");
      planck_assert(map.Nside()==mapT.Nside(),"Nside mismatch");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor*map[m]);
      }
  };

class SZkinEmitter: public Emitter
  {
  private:
    string mapname;

  public:
    SZkinEmitter (paramfile &params)
      {
      mapname = params.find<string> ("SZKIN_LS_map");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &,
                               Healpix_Map<float> &) const
      {
      SZkin_Spectrum spec(tcmb);
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor = 0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double specval = spec.value_avg(freq[m],freq[m+1]);
        double dflux = wgt[m]*1e20*specval;
#ifdef CALC_FLUX
        factor += dflux;
#else
        factor += dflux*antennaFact(0.5*(freq[m]+freq[m+1]));
#endif
        }
      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (mapname,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor*map[m]);
      }
  };

class FreeFreeEmitter: public Emitter
  {
  private:
    string mapname;

  public:
    FreeFreeEmitter (paramfile &params)
      {
      mapname = params.find<string> ("FREEFREE_LS_map");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &,
                               Healpix_Map<float> &) const
      {
      FreeFree_Spectrum spec(1e4);
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor = 0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double specval = spec.value_avg(freq[m],freq[m+1]);
        double dflux = wgt[m]*1e20*specval;
#ifdef CALC_FLUX
        factor += dflux;
#else
        factor += dflux*antennaFact(0.5*(freq[m]+freq[m+1]));
#endif
        }
      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (mapname,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor*map[m]);
      }
  };

class COEmitter: public Emitter
  {
  private:
    string mapname;
    double t_co;

  public:
    COEmitter (paramfile &params)
      {
      mapname = params.find<string> ("CO_LS_map");
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &,
                               Healpix_Map<float> &) const
      {
      CO_Spectrum spec(t_co);
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      double factor = 0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double specval = spec.value_avg(freq[m],freq[m+1]);
        double dflux = wgt[m]*1000/speedOfLight*1e20*specval;
#ifdef CALC_FLUX
        factor += dflux;
#else
        factor += dflux*antennaFact(0.5*(freq[m]+freq[m+1]));
#endif
        }
      Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
      map_import (mapname,map,"I_Stokes");
      for (int m=0; m<map.Npix(); ++m)
        mapT[m] += float(factor*map[m]);
      }
  };

class TabularEmitter: public Emitter
  {
  private:
    arr<string> mapname;
    arr<double> mapfreq;
    bool polarised;

  public:
    TabularEmitter (paramfile &params, int idx)
      {
      string prefix = string("tabular_")+dataToString(idx)+"_";
      polarised = params.find<bool> (prefix+"polarised",true);
      int nmaps = params.find<int> (prefix+"nmaps");
      planck_assert(nmaps>=1,"TabularEmitter: need at least one map");
      mapname.alloc(nmaps);
      if (nmaps>1) mapfreq.alloc(nmaps);
      for (int m=0; m<nmaps; ++m)
        {
        mapname[m] = params.find<string>
                     (prefix+"map_"+dataToString(m+1));
        if (nmaps>1)
          mapfreq[m] = params.find<double>
                       (prefix+"freq_"+dataToString(m+1));
        }
      }
    virtual void add_emission (const Detector_Response &resp,
                               Healpix_Map<float> &mapT,
                               Healpix_Map<float> &mapQ,
                               Healpix_Map<float> &mapU) const
      {
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);
      tsize nmaps = mapname.size();
      arr<double> factor(nmaps,0);
      tsize il=0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double fmid = 0.5*(freq[m]+freq[m+1]);
        if (nmaps>1)
          {
          while ((il<(nmaps-2)) && (fmid>mapfreq[il+1]))
            ++il;
          double frac = (fmid-mapfreq[il])/(mapfreq[il+1]-mapfreq[il]);
          double dfluxl = (1-frac)*wgt[m];
          double dfluxr = frac*wgt[m];
#ifdef CALC_FLUX
          factor[il]   += dfluxl;
          factor[il+1] += dfluxr;
#else
          factor[il]   += dfluxl*antennaFact(fmid);
          factor[il+1] += dfluxr*antennaFact(fmid);
#endif
          }
        else
          {
#ifdef CALC_FLUX
          factor[0] += wgt[m];
#else
          factor[0] += wgt[m]*antennaFact(fmid);
#endif
          }
        }
      for (tsize i=0; i<nmaps; ++i)
        {
        Healpix_Map<float> map(mapT.Nside(),mapT.Scheme(),SET_NSIDE);
        if (factor[i]!=0)
          {
          map_import (mapname[i],map,"I_Stokes");
          for (int m=0; m<map.Npix(); ++m)
            mapT[m] += float(factor[i]*map[m]);
          if (polarised)
            {
            map_import (mapname[i],map,"Q_Stokes");
            for (int m=0; m<map.Npix(); ++m)
              mapQ[m] += float(factor[i]*map[m]);
            map_import (mapname[i],map,"U_Stokes");
            for (int m=0; m<map.Npix(); ++m)
              mapU[m] += float(factor[i]*map[m]);
            }
          }
        }
      }
  };

void collect_emitters (paramfile &params, vector<Emitter *> &emitters)
  {
  if (params.find<bool>("calc_CMB",false))
    emitters.push_back (new CmbEmitter(params));
  if (params.find<bool>("calc_CMB_PRSM",false))
    emitters.push_back (new CmbEmitter_PRSM(params));
  if (params.find<bool>("calc_DUST_PRSM",false))
    emitters.push_back (new DustEmitter_PRSM(params));
  if (params.find<bool>("calc_SYNCHRO_PRSM",false))
    emitters.push_back (new SynchroEmitter_PRSM(params));
  if (params.find<bool>("calc_SZ_PRSM",false))
    emitters.push_back (new SZEmitter_PRSM(params));
  if (params.find<bool>("calc_RADIO_PRSM",false))
    emitters.push_back (new RadioEmitter_PRSM(params));
//  if (params.find<bool>("calc_INFRARED_PRSM",false))
//    emitters.push_back (new InfraredEmitter_PRSM(params));
  if (params.find<bool>("calc_FREEFREE_PRSM",false))
    emitters.push_back (new FreeFreeEmitter_PRSM(params));
  if (params.find<bool>("calc_SZ",false))
    emitters.push_back (new SZEmitter(params));
  if (params.find<bool>("calc_SZKIN",false))
    emitters.push_back (new SZkinEmitter(params));
  if (params.find<bool>("calc_FREEFREE",false))
    emitters.push_back (new FreeFreeEmitter(params));
  if (params.find<bool>("calc_CO",false))
    emitters.push_back (new COEmitter(params));

  int ntab = params.find<int> ("n_tabular",0);
  planck_assert (!params.param_present("tabular_0_map_0"),
    "Key 'tabular_0_map_0' found, but first component must have index 1");
  for (int m=1; m<=ntab; ++m)
    emitters.push_back (new TabularEmitter(params,m));
  }
