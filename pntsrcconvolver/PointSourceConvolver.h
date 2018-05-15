/*
 *  This file is part of the Planck simulation package
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
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  This code is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  \author Martin Reinecke
 *  \author Reinhard Hell
 */

#ifndef PLANCK_POINTSOURCECONVOLVER_H
#define PLANCK_POINTSOURCECONVOLVER_H

#include <string>
#include <vector>
#include <set>
#include "healpix_map.h"
#include "paramfile.h"
#include "pointing.h"
#include "vec3.h"
#include "focalplane_db.h"
#include "sat_info.h"
#include "safe_ptr.h"
#include "geom_utils.h"
#include "BeamCharacteristic.h"
#include "lsconstants.h"

namespace psc {

class interval
  {
  public:
    int lo, hi;

    interval (int Lo, int Hi)
      : lo(Lo), hi(Hi) {}
    interval (int val)
      : lo(val), hi(val) {}

    bool operator< (const interval &other) const
      { return (lo<other.lo); }

    bool tryToAdd (int val)
      {
      if ((val>=lo) && (val<=hi)) return true;
      if (val==lo-1) { --lo; return true; }
      if (val==hi+1) { ++hi; return true; }
      return false;
      }
  };

class SmallPtSrc
  {
  protected:
    bool fixed;
    pointing locP1, locP2;

  public:
    double flux, polangle, polfrac;

    SmallPtSrc()
      : fixed(true), locP1(0,0), locP2(0,0), flux(0), polangle(0), polfrac(0) {}
    SmallPtSrc (const pointing &ptg, double fl, double polang, double polfr)
      : fixed(true), locP1(ptg), locP2(ptg), flux(fl), polangle(polang),
        polfrac(polfr) {}
    SmallPtSrc (const pointing &ptg1, const pointing &ptg2, double fl,
        double polang, double polfr)
      : fixed(false), locP1(ptg1), locP2(ptg2), flux(fl), polangle(polang),
        polfrac(polfr)
        {
        const double min_angle = degr2rad/(60.*60.); // 1 arcsec
        if (v_angle(locP1,locP2)<min_angle)
          {
          locP1 = locP2 = locP(0.5);
          fixed = true;
          }
        }

    pointing locP (double frac) const
      {
      if (fixed) return locP1;
      return pointing (vec3(locP1)*(1.-frac) + vec3(locP2)*frac);
      }
  };

class QuickPtSrc: public SmallPtSrc
  {
  private:
    vec3 locV1, locV2;

  public:
    unsigned int id;

    QuickPtSrc (const SmallPtSrc &sps, int id_)
      : SmallPtSrc(sps), locV1(locP1), locV2(locP2), id(id_) {}
    QuickPtSrc (const SmallPtSrc &sps, double factor, int id_)
      : SmallPtSrc(sps), locV1(locP1), locV2(locP2), id(id_)
      { flux*=factor; }

    vec3 locV (double frac) const
      {
      if (fixed) return locV1;
      vec3 tmp = locV1*(1.-frac) + locV2*frac;
      tmp.Normalize();
      return tmp;
      }
  };

class PtSrc
  {
  protected:
    std::string name_;

  public:
    PtSrc (const std::string &nm)
      : name_(nm) {}
    virtual ~PtSrc() {}
    const std::string &name() const { return name_; }
    virtual SmallPtSrc getData (int pp) = 0;
  };

class FixedPtSrc: public PtSrc
  {
  private:
    SmallPtSrc sps;

  public:
    FixedPtSrc (const SmallPtSrc &sps_, const std::string &nm)
      : PtSrc(nm), sps(sps_) {}
    virtual SmallPtSrc getData (int /*pp*/)
      { return sps; }
  };

class VariablePtSrc: public PtSrc
  {
  private:
    arr<SmallPtSrc> sps;

  public:
    VariablePtSrc(int npp, const std::string &nm)
      : PtSrc(nm), sps(npp) {}
    void setData (int pp, const SmallPtSrc &ps)
      { sps[pp] = ps; }

    virtual SmallPtSrc getData (int pp)
      { return sps[pp]; }
  };

typedef Healpix_Map<std::vector<QuickPtSrc> > PtSrcMapType;

class PtSrcMapServer
  {
  private:
    PtSrcMapType srcMap;
    std::vector<FixedPtSrc> fixedSources;
    std::vector<VariablePtSrc> variableSources;
    arr<int> variablePixels;
    double varPntsrcFactor;

    void readFixedSources (paramfile &params,bool polarisation);
    void readPlanets (paramfile &params,const Sat_Info &satinfo,double freqGHz);

  public:
    PtSrcMapServer (paramfile &params, const focalplane_db &fpdb,
      const std::string &det_id, const Sat_Info &satinfo);

    const PtSrcMapType &provideMap (int pp)
      {
      // clean up old variable sources
      for (tsize m=0; m<variablePixels.size();++m)
        srcMap[variablePixels[m]].pop_back();
      // add variable sources for this pointing period
      variablePixels.alloc(variableSources.size());
      for (tsize m=0; m<variablePixels.size();++m)
        {
        SmallPtSrc sps = variableSources[m].getData(pp);
        int id = m+fixedSources.size();
        variablePixels[m] = srcMap.ang2pix(sps.locP(0.5));
        srcMap[variablePixels[m]].push_back(QuickPtSrc(sps,varPntsrcFactor,id));
        }
      return srcMap;
      }
    int nSources() const
      { return fixedSources.size()+variableSources.size(); }
    const PtSrc &getPtSrc(unsigned int id) const
      {
      if (id<fixedSources.size()) return fixedSources[id];
      return variableSources[id-fixedSources.size()];
      }
  };

class PointSourceConvolver
  {
  private:
    std::string statFile;
    double beamRadius;
    arr<SmallPtSrc> maxSrc;
    arr<bool> seen;
    arr<int> nrOfObs, ppOfObsMax;
    arr<std::vector<interval> > ppsOfHits;
    PtSrcMapServer mapServer;
    safe_ptr<BeamCharacteristic> beamCharac;

  public:
    PointSourceConvolver(paramfile &params, const focalplane_db &fpdb,
      const std::string &det_id, const Sat_Info &satinfo);

    void Get_Intensities(int pp, const arr<pointing> &detpt,
      const arr<vec3> &vdetpt, const arr<double> &heading,
      arr<double> &intensity);
    void Write_Hitsfile() const;
  };

} //end of namespace psc

#endif
