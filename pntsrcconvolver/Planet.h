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

#ifndef PLANCK_PLANET_H
#define PLANCK_PLANET_H

#include <string>
#include "lsconstants.h"
#include "arr.h"
#include "PointSourceConvolver.h"

class ephemerides;

/*
 *  written by  Dr. Reinhard M. Hell,  Max-Planck-Institute for Astrophysics
 */

namespace psc {

struct PlanetData
  {
  const char *planet_name;
  int rotationdirection; // Moore (2000)
  double rotationperiod; // [s] Moore (2000)
  double equatorialradius, polarradius; // [m] Moore (2000)
  double meansundistance; // [m] Moore (2000)
  double bondalbedo; // Moore (2000)
  double surfaceemissivity; // Neugebauer et al. (1971)
  double heatratio; // Unsoeld and Baschek (1999)
  double coolingfactor; // Milton, ed. (1989); Sinton and Strong (1959)
  double meaneffectivetemp; // [K] Lang (1991)
  double brTemp[9];
  };

const PlanetData MarsData=
  {
  "Mars",
  1,
  88642.6,
  3397.0e3, 3379.5e3,
  1.524*astronomicalUnit,
  0.16,
  0.9,
  1.0,
  0.70,
  212.0,
  { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  };
const PlanetData JupiterData=
  {
  "Jupiter",
  1,
  35729.0,
  71442.0e3, 66854.0e3,
  5.203*astronomicalUnit,
  0.43,
  1.0,
  1.7,
  1.0,
  129.0,
  // (Naselsky et al., 2002)
  {152.0, 158.0, 167.0, 173.0, 176.0, 179.0, 178.0, 133.0, 145.0},
  };
const PlanetData SaturnData=
  {
  "Saturn",
  1,
  38365.0,
  60268.0e3, 54364.0e3,
  9.359*astronomicalUnit,
  0.47,
  1.0,
  1.8,
  1.0,
  97.0,
  // (Naselsky et al., 2002)
  {133.0, 135.0, 138.0, 141.0, 144.0, 146.0, 141.0, 118.0, 114.0},
  };
const PlanetData UranusData=
  {
  "Uranus",
  1,
  62064.0,
  25559.0e3, 24973.0e3,
  19.181*astronomicalUnit,
  0.51,
  1.0,
  1.1,
  0.85,
  58.0,
  // (353-857GHz: Hildebrand et al., 1985; 30-217GHz: estimate using Griffin
  // et al., 1986, Naselsky et al., 2002, and Sinton and Strong, 1959)
  { 79.0,  82.0,  87.0,  90.0,  92.0,  93.0,  91.0,  77.0,  62.5},
  };
const PlanetData NeptuneData=
  {
  "Neptune",
  1,
  58020.0,
  25269.0e3, 24800.0e3,
  30.058*astronomicalUnit,
  0.35,
  1.0,
  2.6,
  1.0,
  56.0,
  // (353-857GHz: Hildebrand et al., 1985; 30-217GHz: estimate using Griffin
  // et al., 1986, Naselsky et al., 2002, and Sinton and Strong, 1959)
  { 63.0,  66.0,  69.0,  72.0,  73.0,  74.0,  87.0,  76.0,  63.0},
  };
const PlanetData PlutoData=
  {
  "Pluto",
  1,
  551858.4,
  1162.0e3, 1162.0e3,
  39.500*astronomicalUnit,
  0.55,
  0.9,
  1.0,
  0.5,
  60.0,
  // (Estimate using Griffin et al., 1986, Naselsky et al., 2002,
  // and Sinton and Strong, 1959)
  { 55.0,  57.0,  61.0,  63.0,  64.0,  65.0,  65.0,  48.0,  53.0},
  };
const PlanetData AsteroidData=
  {
  "Asteroid",
  1,
  21600.0,
  75.0e3, 75.0e3,
  3*astronomicalUnit,
  0.062816,
  0.9,
  1.0,
  0.8,
  151.0,
  { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  };

class Planet: public PlanetData
  {
  protected:
    arr<double> sunDistance, satDistance, sunPlanetSatAngle;
    arr<SmallPtSrc> sources;

    double effectiveTemp(int pp) const; // in K
    double brightnessTemp(double frequency, int pp) const; // in K

    double fact1(int pp) const
      { return pi*equatorialradius*polarradius/pow(satDistance[pp],2); }
    double fact2(int pp) const
      {
      return coolingfactor
            +(1.0-coolingfactor)*0.5*(1.0+cos(sunPlanetSatAngle[pp]));
      }

    double fluxBright (double frequency, int pp) const
      {
      return fact1(pp)*fact2(pp)*surfaceemissivity
            *brightnessTemp(frequency,pp);
      }

    void readEphem(ephemerides &eph,
      const arr<double> &tstart, const arr<double> &tend);

  public:
    Planet(const PlanetData &pdata)
      : PlanetData(pdata) {}

    const arr<SmallPtSrc> &getSources() const { return sources; }
  };

class GasPlanet: public Planet
  {
  public:
    GasPlanet(const PlanetData &pdata, ephemerides &eph,
      const arr<double> &tstart, const arr<double> &tend, double freq)
      : Planet(pdata)
      {
      readEphem(eph,tstart,tend);
      for (tsize i=0; i<sources.size(); ++i)
        sources[i].flux = (sunDistance[i]>=0.) ? fluxBright(freq,i) : -1.;
      }
  };

class RockyPlanet: public Planet
  {
  public:
    RockyPlanet(const PlanetData &pdata, ephemerides &eph,
      const arr<double> &tstart, const arr<double> &tend, double freq)
      : Planet(pdata)
      {
      readEphem(eph,tstart,tend);
      for (tsize i=0; i<sources.size(); ++i)
        sources[i].flux = (sunDistance[i]>=0.) ? fluxBright(freq,i) : -1.;
      }
  };

class MarsPlanet: public Planet
  {
  private:
    arr<double> tempdata;

    void readTempdata(const std::string &objname, const arr<double> &tstart,
      const arr<double> &tend);

    double fluxSpecial (int pp) const
      { return fact1(pp)*tempdata[pp]; }

  public:
    MarsPlanet(const PlanetData &pdata, ephemerides &eph,
      const arr<double> &tstart, const arr<double> &tend,
      const std::string &objname)
      : Planet(pdata)
      {
      readEphem(eph,tstart,tend);
      readTempdata(objname,tstart,tend);
      for (tsize i=0; i<sources.size(); ++i)
        sources[i].flux = (sunDistance[i]>=0.) ? fluxSpecial(i) : -1.;
      }
  };

} // end of namespace psc

#endif
