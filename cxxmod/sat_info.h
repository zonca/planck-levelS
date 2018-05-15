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

#ifndef PLANCK_SAT_INFO_H
#define PLANCK_SAT_INFO_H

#include <string>
#include "arr.h"
#include "pointing.h"
#include "vec3.h"
#include "rotmatrix.h"
#include "iohandle.h"
#include "quaternion.h"
#include "ephemerides.h"

class paramfile;

class Sat_Info
  {
  protected:
    arr<double> tstart, tend, duration;
    safe_ptr<ephemerides> eph;
    bool nominal;

// general cached data:
    int current_period;

// cached data for realistic pointing:
    arr<int32> nval_ptg;
    arr<int64> offset_ptg;

    arr<double> qtime;
    arr<quaternion> rotquat;
    arr<double> rangle, rxsin;
    arr<bool> rotflip;

    safe_ptr<iohandle> inp;

    void computeInterpolation();

  public:
    Sat_Info (bool nompoint, const paramfile &parfile);
    virtual ~Sat_Info() {}

    virtual void gotoPeriod (int period) = 0;
    virtual vec3 velocity () const;

    virtual int numSamples (int period, double f_samp) const;
    virtual int64 startIndex (int period, double f_samp) const;
    virtual void getTimeStamps (int period, double f_samp, arr<double> &times)
      const;

    void getTransform (double time, rotmatrix &rmat) const;
    int numPeriods() const { return tstart.size(); }
    int curPeriod() const { return current_period; }

    double startTime (int period) const { return tstart[period]; }
    const arr<double> &startTimes () const { return tstart; }
    const arr<double> &endTimes () const { return tend; }
    const arr<double> &durations () const { return duration; }
  };

class Sat_Info_Simmission: public Sat_Info
  {
  private:
    arr<pointing> nom_px, nom_py, nom_pz;

// general cached data:
    int thetaxcol, phixcol, thetaycol, phiycol, thetazcol, phizcol;

  public:
    Sat_Info_Simmission (bool nompoint, const paramfile &parfile);
    virtual void gotoPeriod (int period);
  };

class Sat_Info_LFI: public Sat_Info
  {
  private:
    int xcol, ycol, zcol, wcol, tcol;

// data for realistic sampling
    arr<int32> nval_samp;
    arr<int64> offset_samp;

// data for time stamp generation
    arr<int64> s_idx;
    arr<double> s_time;

    double sampletime (int64 samp);

  public:
    Sat_Info_LFI (bool nompoint, const paramfile &parfile);
    virtual void gotoPeriod (int period);
    virtual int numSamples (int period, double f_samp) const;
    virtual int64 startIndex (int period, double f_samp) const;
    virtual void getTimeStamps (int period, double f_samp, arr<double> &times)
      const;
  };

class Sat_Info_HFI: public Sat_Info
  {
  private:
    int xcol, ycol, zcol, wcol, tcol;

    safe_ptr<iohandle> inp_time;

  public:
    Sat_Info_HFI (bool nompoint, paramfile &parfile);
    virtual void gotoPeriod (int period);
    virtual int numSamples (int period, double f_samp) const;
    virtual int64 startIndex (int period, double f_samp) const;
    virtual void getTimeStamps (int period, double f_samp, arr<double> &times)
      const;
  };

#endif
