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
 *  Copyright (C) 2003-2014 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "sat_info.h"
#include "io_utils.h"
#include "iohandle.h"
#include "quaternion.h"
#include "paramfile.h"
#include "ephemerides.h"
#include "lsconstants.h"

using namespace std;

namespace {

void read_entire_timeinfo (iohandle &inp, arr<double> &tstart,
  arr<double> &tend, arr<double> &duration)
  {
  double toffset = sec_58_70;
  if (inp.keyPresent("TREF1958") && inp.getKey<bool>("TREF1958")) toffset=0.;
  inp.readEntireColumn("t_endpt_int",tend);
  arr<double> tmp;
  inp.readEntireColumn("t_endpt_frac",tmp);
  for (tsize i=0; i<tmp.size(); ++i) tend[i]+=tmp[i]+toffset;
  inp.readEntireColumn("t_startpt_int",tstart);
  inp.readEntireColumn("t_startpt_frac",tmp);
  for (tsize i=0; i<tmp.size(); ++i) tstart[i]+=tmp[i]+toffset;

  for (tsize i=0; i<tstart.size()-1; ++i)
    planck_assert (approx(tend[i],tstart[i+1]),
      "gap between pointing periods detected!");

  duration.alloc(tstart.size());
  for (tsize i=0; i<duration.size(); ++i) duration[i]=tend[i]-tstart[i];
  }

} // unnamed namespace

Sat_Info::Sat_Info (bool nompoint, const paramfile &parfile)
  : nominal(nompoint), current_period(-1)
  {
  eph = getEphemerides(parfile.find<string>("satinfo_ephemeris"));
  eph->loadBody("Solar System Barycenter");
  }

void Sat_Info::computeInterpolation()
  {
  rangle.alloc(nval_ptg[current_period]-1);
  rxsin.alloc(nval_ptg[current_period]-1);
  rotflip.alloc(nval_ptg[current_period]-1);

#pragma omp parallel
{
  int m;
#pragma omp for schedule (static)
  for (m=0; m<nval_ptg[current_period]-1; ++m)
    {
    quaternion delta(rotquat[m+1]*rotquat[m].conj());
    rotflip[m]=false;
    if (delta.w < 0.)
      { rotflip[m]=true; delta.Flip(); }
    vec3 v;
    double omega;
    delta.toAxisAngle(v,omega);
    rangle[m]=omega*.5;
    rxsin[m]=1./sin(rangle[m]);
    }
}
  }

//virtual
vec3 Sat_Info::velocity() const
  {
  double tlight=462; // light time correction proposed by Reijo Keskitalo
  return (eph->posRelSat_m(tend[current_period]+tlight)
          -eph->posRelSat_m(tstart[current_period]+tlight))
         /(-duration[current_period]);
  }

//virtual
int Sat_Info::numSamples (int period, double f_samp) const
  {
  if (nominal) return nearest<int>(60*f_samp);

  double te = tend[period];
  if (period<(int(tstart.size())-1))
    te = tstart[period+1];

  int64 first_sample = nearest<int64>((tstart[period]-tstart[0])*f_samp);
  int64 last_sample = nearest<int64>((te-tstart[0])*f_samp);
  return last_sample-first_sample;
  }

//virtual
int64 Sat_Info::startIndex (int period, double f_samp) const
  {
  if (nominal) return period*nearest<int>(60*f_samp);

  return nearest<int64>((tstart[period]-tstart[0])*f_samp);
  }

//virtual
void Sat_Info::getTimeStamps (int period, double f_samp, arr<double> &times)
  const
  {
  if (nominal)
    {
    times.alloc(nearest<int>(60*f_samp));
    for (tsize m=0; m<times.size(); ++m)
      times[m] = tstart[period]+m/f_samp;
    }
  else
    {
    double te = tend[period];
    if (period<(int(tstart.size())-1))
      te = tstart[period+1];

    int64 first_sample = nearest<int64>((tstart[period]-tstart[0])*f_samp);
    int64 last_sample = nearest<int64>((te-tstart[0])*f_samp);
    times.alloc(last_sample-first_sample);
    for (tsize m=0; m<times.size(); ++m)
      times[m] = tstart[0]+(m+first_sample)/f_samp;
    }
  }

void Sat_Info::getTransform (double time, rotmatrix &rmat) const
  {
  tsize idx;
  double frac;
  qtime.interpol_helper(time,idx,frac);
  double omega = rangle[idx];
  double xsin = rxsin[idx];
  double w1 = sin((1.-frac)*omega)*xsin,
         w2 = sin(frac*omega)*xsin;
  if (rotflip[idx]) w1=-w1;
  const quaternion &q1(rotquat[idx]), &q2(rotquat[idx+1]);
  quaternion q (w1*q1.w + w2*q2.w,
                w1*q1.x + w2*q2.x,
                w1*q1.y + w2*q2.y,
                w1*q1.z + w2*q2.z);
  q.make_rotmatrix(rmat);
  }

Sat_Info_Simmission::Sat_Info_Simmission (bool nompoint,
  const paramfile &parfile)
  : Sat_Info (nompoint,parfile)
  {
  inp = HandleManager.openObject(parfile.find<string>("sat_info"),
    "satelliteinfo.LS_satinfo");

  string axcon = inp->getKey<string>("sataxis_convention");
  planck_assert((axcon=="hfi")||(axcon=="HFI"),
    "wrong axis convention in the satellite pointing file");

  read_entire_timeinfo (*inp, tstart, tend, duration);

  if (nominal)
    {
    readEntirePointing (*inp, "theta_x_startpt","phi_x_startpt",nom_px);
    readEntirePointing (*inp, "theta_y_startpt","phi_y_startpt",nom_py);
    readEntirePointing (*inp, "theta_z_startpt","phi_z_startpt",nom_pz);
    nval_ptg.alloc(nom_px.size());
    nval_ptg.fill(2);
    }
  else
    {
    thetaxcol = inp->columnNumber("theta_x");
    phixcol = inp->columnNumber("phi_x");
    thetaycol = inp->columnNumber("theta_y");
    phiycol = inp->columnNumber("phi_y");
    thetazcol = inp->columnNumber("theta_z");
    phizcol = inp->columnNumber("phi_z");
    planck_assert(inp->columnLength(thetaxcol)>0,
      "full pointing information is not available");
    inp->readEntireColumn("nsamples",nval_ptg);
    planck_assert (nval_ptg.size()>0,"nsamples column is empty");
    offset_ptg.alloc(nval_ptg.size());
    offset_ptg[0] = 0;
    for (tsize m=1; m<nval_ptg.size(); ++m)
      offset_ptg[m] = offset_ptg[m-1] + nval_ptg[m-1];
    }
  }

//virtual
void Sat_Info_Simmission::gotoPeriod (int period)
  {
  if (period==current_period) return;

  rotquat.alloc(nval_ptg[period]);
  qtime.alloc(nval_ptg[period]);

  if (nominal)
    {
    qtime[0]=tstart[period];
    qtime[1]=qtime[0]+1.;
    rotquat[0] = rotmatrix(nom_px[period].to_vec3(),
      nom_py[period].to_vec3(), nom_pz[period].to_vec3());
    quaternion t1 (vec3(1,0,0),twopi/60.);
    rotquat[1]=rotquat[0]*t1;
    }
  else
    {
    arr<pointing> real_px(nval_ptg[period]),
                  real_py(nval_ptg[period]),
                  real_pz(nval_ptg[period]);
    readPointing (*inp,thetaxcol,phixcol,thetaycol,phiycol,thetazcol,phizcol,
                  real_px,real_py,real_pz,offset_ptg[period]);
    double x_dt_in = nval_ptg[period]/duration[period];

#pragma omp parallel
{
    int m;
#pragma omp for schedule (static)
    for (m=0; m<nval_ptg[period]; ++m)
      {
      rotquat[m]=quaternion(rotmatrix(real_px[m].to_vec3(),real_py[m].to_vec3(),
                                      real_pz[m].to_vec3()));
      qtime[m] = tstart[period]+(m+0.5)/x_dt_in;
      }
}
    }

  current_period = period;
  computeInterpolation();
  }

double Sat_Info_LFI::sampletime (int64 samp)
  {
  tsize iv=std::lower_bound (s_idx.begin(),s_idx.end(),samp)-s_idx.begin();
  if ((samp==s_idx[s_idx.size()-1]) && (iv==s_idx.size()-1))
    --iv;
  planck_assert(iv<(s_idx.size()-1), "leaving the last interval");
// FIXME: sanity checks (only temporary)
  planck_assert(s_idx[iv]<=samp,"must not happen");
  planck_assert(s_idx[iv+1]>=samp,"must not happen");
  return s_time[iv]
    + ((samp-s_idx[iv])*(s_time[iv+1]-s_time[iv]))/(s_idx[iv+1]-s_idx[iv]);
  }

Sat_Info_LFI::Sat_Info_LFI (bool nompoint, const paramfile &parfile)
  : Sat_Info(false,parfile)
  {
  planck_assert(!nompoint, "Sat_Info_LFI does not support nominal pointing");

  inp = HandleManager.openObject(parfile.find<string>("satinfo_input"),
    "sat.LS_satpoint_real");

  xcol = inp->columnNumber("quat_x");
  ycol = inp->columnNumber("quat_y");
  zcol = inp->columnNumber("quat_z");
  wcol = inp->columnNumber("quat_w");
  tcol = inp->columnNumber("quat_time");
  uint64 nquat=inp->columnLength(xcol);
  planck_assert(multiequal(nquat,inp->columnLength(ycol),
    inp->columnLength(zcol),inp->columnLength(wcol),inp->columnLength(tcol)),
    "array size mismatch");

  inp->readEntireColumn ("pp_ifirst", offset_ptg);
  inp->readEntireColumn ("pp_nquat", nval_ptg);
  tsize npp = offset_ptg.size();
  planck_assert(npp==nval_ptg.size(), "array size mismatch");

  safe_ptr<iohandle> inp2 (HandleManager.openObject
    (parfile.find<string>("satinfo_sampleinfo"),"mission.LS_samples_real"));
  inp2->readEntireColumn ("pp_nsamp",nval_samp);
  planck_assert(nval_samp.size()==npp,"array size mismatch");
  offset_samp.alloc(npp);
  offset_samp[0]=0;
  for (tsize m=1; m<npp; ++m)
    offset_samp[m] = offset_samp[m-1]+nval_samp[m-1];
  inp2->readEntireColumn ("sample_index",s_idx);
  inp2->readEntireColumn ("sample_time",s_time);
  planck_assert(s_idx.size()==s_time.size(),"array size mismatch");

  tstart.alloc(npp); tend.alloc(npp); duration.alloc(npp);
  tstart[0]=sampletime(offset_samp[0]);
  for (tsize m=1; m<npp; ++m)
    tend[m-1]=tstart[m]=sampletime(offset_samp[m]);
  tend[npp-1]=sampletime(offset_samp[npp-1]+nval_samp[npp-1]-1);
  for (tsize m=0; m<npp; ++m)
    duration[m] = tend[m]-tstart[m];
  }

//virtual
void Sat_Info_LFI::gotoPeriod (int period)
  {
  if (period==current_period) return;

  rotquat.alloc(nval_ptg[period]);
  qtime.alloc(nval_ptg[period]);
  readQuaternions (*inp,xcol,ycol,zcol,wcol,rotquat,offset_ptg[period]);
  inp->readColumn(tcol,qtime,offset_ptg[period]);

  current_period = period;
  computeInterpolation();
  }

//virtual
int Sat_Info_LFI::numSamples (int period, double f_samp) const
  {
  if (f_samp<0) return nval_samp[period];
  return Sat_Info::numSamples(period,f_samp);
  }

//virtual
int64 Sat_Info_LFI::startIndex (int period, double f_samp) const
  {
  if (f_samp<0) return offset_samp[period];
  return Sat_Info::startIndex(period,f_samp);
  }

//virtual
void Sat_Info_LFI::getTimeStamps (int period, double f_samp,
  arr<double> &times) const
  {
  if (f_samp<0)
    {
    planck_assert(period==current_period,"period must be the current period");
    times.alloc(nval_samp[period]);
    tsize iv=std::lower_bound
      (s_idx.begin(),s_idx.end(),offset_samp[period])-s_idx.begin();
    planck_assert(iv<(s_idx.size()-1), "leaving the last interval");
// FIXME: sanity checks (only temporary)
    planck_assert(s_idx[iv]<=offset_samp[period],"must not happen");
    planck_assert(s_idx[iv+1]>offset_samp[period],"must not happen");
    for (tsize m=0; m<times.size(); ++m)
      {
      int64 idx = m+offset_samp[period];
      if (idx>s_idx[iv+1]) ++iv;
      planck_assert(iv<(s_idx.size()-1), "leaving the last interval");
      times[m] = s_time[iv]
        + ((idx-s_idx[iv])*(s_time[iv+1]-s_time[iv]))/(s_idx[iv+1]-s_idx[iv]);
      }
    return;
    }
  return Sat_Info::getTimeStamps(period,f_samp,times);
  }

Sat_Info_HFI::Sat_Info_HFI (bool nompoint, paramfile &parfile)
  : Sat_Info(false,parfile)
  {
  planck_assert(!nompoint, "Sat_Info_HFI does not support nominal pointing");

  inp = HandleManager.openObject(parfile.find<string>("satinfo_quaternions"),
    "quat.LS_satpt_quat");
  inp_time = HandleManager.openObject(parfile.find<string>("satinfo_ctr"),
    "ctr.LS_ctr");

  xcol = inp->columnNumber("quat_x");
  ycol = inp->columnNumber("quat_y");
  zcol = inp->columnNumber("quat_z");
  wcol = inp->columnNumber("quat_w");
  tcol = inp_time->columnNumber("CTR");

  {
  safe_ptr<iohandle> inp_idx (HandleManager.openObject
    (parfile.find<string>("satinfo_indexobject"),"table.LS_toi_index_table"));
  inp_idx->readEntireColumn ("first_sample",offset_ptg);
  arr<int64> tmp;
  inp_idx->readEntireColumn ("last_sample",tmp);
  planck_assert(offset_ptg.size()==tmp.size(),"array size mismatch");
  nval_ptg.alloc(offset_ptg.size());
  for (tsize m=0; m<tmp.size(); ++m)
    nval_ptg[m] = tmp[m]-offset_ptg[m]+1;
  }

  tsize nperiods = nval_ptg.size();

  tstart.allocAndFill(nperiods,-1e30);
  tend.allocAndFill(nperiods,-1e30);
  duration.allocAndFill(nperiods,-1e30);

//FIXME: this is still quite inefficient

  int first = parfile.find<int>("first_pointing",1);
  if (first==-1) first=1;
  --first;
  int last = parfile.find<int>("last_pointing",nperiods);
  if (last==-1) last=nperiods;
  --last;
  planck_assert ((first>=0)&&(last>=0)&&(first<=last)&&(last<int(nperiods)),
    "inconsistent first and last pointing periods!");

  int64 rawtime;
  rawtime = 1621174818000457763LL; // start of the real mission in nanosecs
  tstart[0] = rawtime*1e-9;
  for (int m=first; m<=last; ++m)
    {
    inp_time->readColumn (tcol,rawtime,offset_ptg[m]);
    tstart[m] = rawtime*1e-9;
    }
  for (int m=first; m<last; ++m)
    tend[m] = tstart[m+1];
  (last==int(nperiods)-1) ?
    inp_time->readColumn (tcol,rawtime,inp_time->columnLength(tcol)-1) :
    inp_time->readColumn (tcol,rawtime,offset_ptg[last+1]);
  tend[last] = rawtime*1e-9;
  if (last<int(nperiods)-1) tstart[last+1]=tend[last];
  for (int m=first; m<=last; ++m)
    duration[m] = tend[m]-tstart[m];
  }

//virtual
void Sat_Info_HFI::gotoPeriod (int period)
  {
  if (period==current_period) return;

  rotquat.alloc(nval_ptg[period]);
  qtime.alloc(nval_ptg[period]);
  readQuaternions (*inp,xcol,ycol,zcol,wcol,rotquat,offset_ptg[period]);
  arr<int64> rawtime(nval_ptg[period]);
  inp_time->readColumn(tcol,rawtime,offset_ptg[period]);
  for (int m=0; m<nval_ptg[period]; ++m)
    qtime[m] = 1e-9*rawtime[m];
  current_period = period;
  computeInterpolation();
  }

//virtual
int Sat_Info_HFI::numSamples (int period, double f_samp) const
  {
  if (f_samp<0) return nval_ptg[period];
  return Sat_Info::numSamples (period, f_samp);
  }

//virtual
int64 Sat_Info_HFI::startIndex (int period, double f_samp) const
  {
  if (f_samp<0) return offset_ptg[period];
  return Sat_Info::startIndex (period, f_samp);
  }

//virtual
void Sat_Info_HFI::getTimeStamps (int period, double f_samp,
  arr<double> &times) const
  {
  if (f_samp<0)
    {
    planck_assert(period==current_period,"period must be the current period");
    times=qtime;
    return;
    }
  return Sat_Info::getTimeStamps (period, f_samp, times);
  }
