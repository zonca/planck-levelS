/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Utilities for accessing ephemeris data
 *
 *  Copyright (C) 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "ephemerides.h"
#include "lsconstants.h"

using namespace std;

class ephemerides_cllio: public ephemerides
  {
  private:
    arr<string> body, scalar, array;
    arr<double> data;
    tsize nsamp;
    double t0, inv_dt;
    safe_ptr<iohandle> inp;
    tsize o_x, o_y, o_z, o_angdiam, o_dsun, o_sto;

    tsize getScalarOffset (const string &name) const
      { return scalar.find(name); }

    double getScalar (tsize offset) const
      { return data[offset]; }

    double getScalar (const string &name) const
      { return getScalar(getScalarOffset(name)); }

    tsize getArrayOffset (const string &name) const
      { return scalar.size() + nsamp*array.find(name); }

    struct inter_struct
      {
      double frac;
      tsize idx;
      };

    inter_struct getArrayInterpol (double time) const
      {
      time -= t0;
      planck_assert(time>=0.,"requested time < t0");
      inter_struct inter;
      inter.frac = time*inv_dt;
      inter.idx = tsize(inter.frac);
      planck_assert(inter.idx<nsamp-1,"requested time too large");
      inter.frac -= double(inter.idx);
      return inter;
      }

    double getArrayValue (tsize offset, const inter_struct &inter) const
      {
      tsize ofs = offset + inter.idx;
      return (1.-inter.frac)*data[ofs] + inter.frac*data[ofs+1];
      }
    double getArrayValue (tsize offset, double time) const
      { return getArrayValue(offset,getArrayInterpol(time)); }
    double getArrayValue (const string &name, double time) const
      { return getArrayValue (getArrayOffset(name), time); }

  public:
    ephemerides_cllio (const string &obj)
      {
      inp = HandleManager.openObject(obj,"ephemeris.LS_ephemeris");
      nsamp = inp->getKey<int64>("nsamples");
      t0 = inp->getKey<double>("T0");
      inv_dt = 1./inp->getKey<double>("deltaT");
      inp->readEntireColumn("bodies",body);
      inp->readEntireColumn("scalars",scalar);
      inp->readEntireColumn("arrays",array);

      data.alloc (scalar.size() + nsamp*array.size());

      planck_assert(data.size()>0,"no data fields in ephemerides file");
      planck_assert(body.size()>0,"no bodies in ephemerides file");

      o_x = getArrayOffset("x_rel_sat_ecl/m");
      o_y = getArrayOffset("y_rel_sat_ecl/m");
      o_z = getArrayOffset("z_rel_sat_ecl/m");
      o_angdiam = getArrayOffset("angdiam_rel_sat/arcsec");
      o_dsun = getArrayOffset("d_sun/m");
      o_sto = getArrayOffset("angle_S_T_O/deg");
      }

    virtual const arr<string> &getBodies () const
      { return body; }

    virtual void loadBody (const string &name)
      { inp->readColumn("data",data,body.find(name)*data.size()); }

    virtual vec3 posRelSat_m (double time) const
      {
      inter_struct inter = getArrayInterpol (time);
      return vec3 (getArrayValue (o_x, inter),
                   getArrayValue (o_y, inter),
                   getArrayValue (o_z, inter));
      }
    virtual double dSun_m (double time) const
      { return getArrayValue (o_dsun, time); }
    virtual double angdiam_arcmin (double time) const
      { return 1./60.*getArrayValue (o_angdiam, time); }
    virtual double angle_S_T_O_rad (double time) const
      { return degr2rad*getArrayValue (o_sto, time); }
  };

ephemerides *getEphemerides (const string &obj)
  { return new ephemerides_cllio (obj); }
