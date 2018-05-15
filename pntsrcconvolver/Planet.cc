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

#include "Planet.h"
#include "iohandle.h"
#include "ephemerides.h"

namespace psc {

using namespace std;

double Planet::effectiveTemp(int pp) const
  {
  // effective temperature in K:
  return meaneffectivetemp*sqrt(meansundistance/sunDistance[pp]);
  }

double Planet::brightnessTemp(double frequency, int pp) const
  {
  if (brTemp[0]<0) return effectiveTemp(pp);

  const double ofreq[] // Frequencies of Planck satellite
    = { 30.0, 44.0, 70.0, 100.0, 143.0, 217.0, 353.0, 545.0, 857.0};

  tsize idx;
  double frac;
  interpol_helper(&ofreq[0],&ofreq[0]+9,frequency,idx,frac);

  return sqrt(meansundistance/sunDistance[pp])
    *(brTemp[idx] + frac*(brTemp[idx+1]-brTemp[idx]));
  }

void Planet::readEphem(ephemerides &eph,
      const arr<double> &tstart, const arr<double> &tend)
  {
  int nperiods = tstart.size();
  eph.loadBody(planet_name);

  sources.alloc(nperiods);
  sunDistance.allocAndFill(nperiods,-1.);
  satDistance.allocAndFill(nperiods,-1.);
  sunPlanetSatAngle.allocAndFill(nperiods,-1.);
  for (int m=0; m<nperiods; ++m)
    {
    if ((tstart[m]>=0.) && (tend[m]>=0.)) // pointing period fully valid
      {
      vec3 pos1(eph.posRelSat_m(tstart[m])),
           pos2(eph.posRelSat_m(tend[m]));
      satDistance[m] = 0.5*(pos1.Length()+pos2.Length());
      sunDistance[m] = 0.5*(eph.dSun_m(tstart[m]) + eph.dSun_m(tend[m]));
      sunPlanetSatAngle[m] = 0.5*
        (eph.angle_S_T_O_rad(tstart[m])+eph.angle_S_T_O_rad(tend[m]));
      sources[m] = SmallPtSrc (pointing(pos1), pointing(pos2), 0, 0, 0);
      }
    }
  }

void MarsPlanet::readTempdata (const string &objname, const arr<double> &tstart,
  const arr<double> &tend)
  {
  tempdata.allocAndFill(tstart.size(),-1e30);
  arr<double> temps;
  safe_ptr<iohandle> inp (HandleManager.openObject(objname,"toi.LS_toi"));
  inp->readEntireColumn("signal",temps);
  double t0 = inp->getKey<double>("tstart"),
         dt = inp->getKey<double>("timestep");
  for (tsize i=0; i<tstart.size(); ++i)
    {
    if ((tstart[i]>=0.) && (tend[i]>=0.)) // pointing period fully valid
      {
      double t = (0.5*(tstart[i]+tend[i])-t0)/dt;
      planck_assert ((t>=0.0)&&(t<temps.size()-1),"bad time");
      tsize it = tsize(t);
      double frac = t-it;
      tempdata[i] = (1.-frac)*temps[it] + frac*temps[it+1];
      }
    }
  }

} // end of namespace psc
