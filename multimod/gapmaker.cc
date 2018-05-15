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
 *  Copyright (C) 2008-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "lsconstants.h"
#include "arr.h"
#include "announce.h"
#include "paramfile.h"
#include "sat_info.h"
#include "focalplane_db.h"
#include "iohandle_current.h"

using namespace std;

namespace {

template <typename T> class interval
  {
  public:
    T lo, hi;

    interval (T Lo, T Hi)
      : lo(Lo), hi(Hi) {}

    bool operator< (const interval &other) const
      { return (lo<other.lo); }
  };

template<typename T> void copyInterval (iohandle &inp, iohandle &out,
  int col, int64 lo, int64 hi)
  {
  if (lo<0) lo=0;
  if (hi>int64(inp.columnLength(col))) hi=inp.columnLength(col);
  if (hi<=lo) return;

  const uint64 chunksize = 1024*256;
  int64 offset=0;
  arr<T> data;

  while(offset<(hi-lo))
    {
    uint64 nsamples = min(chunksize,uint64(hi-lo-offset));
    data.alloc(nsamples);
    inp.readColumn(col,data,lo+offset);
    out.appendColumn(col, data);
    offset+=chunksize;
    }
  }

template<typename T> void cutColumn (iohandle &inp, iohandle &out,
  int col, int64 minidx, const vector<interval<int64> > &gapsample)
  {
  if (gapsample.size()==0)
    copyInterval<T>(inp, out, col, 0, inp.columnLength(col));
  else
    {
    copyInterval<T>(inp, out, col, 0, gapsample[0].lo-minidx);
    for (unsigned int i=0; i<gapsample.size()-1; ++i)
      {
      copyInterval<T>(inp, out, col, gapsample[i].hi-minidx,
        gapsample[i+1].lo-minidx);
      }
    copyInterval<T>(inp, out, col, gapsample[gapsample.size()-1].hi-minidx,
      inp.columnLength(col));
    }
  }

void appendData (iohandle &out, int col, int64 lo, int64 hi,
  int64 minidx, int64 maxidx, bool value)
  {
  if (lo<minidx) lo=minidx;
  if (hi>maxidx) hi=maxidx;
  if (hi<=lo) return;

  const uint64 chunksize = 1024*256;
  int64 offset=0;
  arr<bool> data;

  while(offset<(hi-lo))
    {
    uint64 nsamples = min(chunksize,uint64(hi-lo-offset));
    data.alloc(nsamples);
    data.fill(value);
    out.appendColumn(col, data);
    offset+=chunksize;
    }
  }

void write_flagfile (paramfile &par_file,
  const vector<interval<int64> > &gapsample, int64 firstsample,
  int64 lastsample)
  {
  string flagfile = par_file.find<string> ("flagfile","");
  if (flagfile=="") return;

  safe_ptr<iohandle> out (HandleManager.createObject(flagfile,"flag.LS_flag"));
  int colnum = out->columnNumber("flag");

  if (gapsample.size()==0)
    appendData (*out,colnum,firstsample,lastsample,
                firstsample,lastsample,true);
  else
    {
    appendData (*out,colnum,firstsample,gapsample[0].lo,
      firstsample,lastsample,true);
    for (unsigned int i=0; i<gapsample.size()-1; ++i)
      {
      appendData (*out,colnum,gapsample[i].lo,gapsample[i].hi,
        firstsample,lastsample,false);
      appendData (*out,colnum,gapsample[i].hi,gapsample[i+1].lo,
        firstsample,lastsample,true);
      }
    appendData (*out,colnum,gapsample[gapsample.size()-1].lo,
      gapsample[gapsample.size()-1].hi,firstsample,lastsample,false);
    appendData (*out,colnum,gapsample[gapsample.size()-1].hi,
      lastsample,firstsample,lastsample,true);
    }
  }

void cut_timestamps (paramfile &par_file,
  const vector<interval<int64> > &gapsample, int64 firstsample)
  {
  string input = par_file.find<string>("timestamp_in","");
  if (input=="") return;
  string output = par_file.find<string>("timestamp_out");
  safe_ptr<iohandle> inp (HandleManager.openObject(input,"toi.LS_timestamps"));
  safe_ptr<iohandle> out (HandleManager.createObject
    (output,"toi.LS_timestamps"));

  cutColumn<double>(*inp,*out,inp->columnNumber("timestamp"),
    firstsample,gapsample);
  }

void cut_toi (paramfile &par_file,
  const vector<interval<int64> > &gapsample, int64 firstsample)
  {
  string input = par_file.find<string>("toi_in","");
  if (input=="") return;
  string output = par_file.find<string>("toi_out");
  safe_ptr<iohandle> inp (HandleManager.openObject(input,"toi.LS_toi"));
  safe_ptr<iohandle> out (HandleManager.createObject(output,"toi.LS_toi"));

  cutColumn<float>(*inp,*out,inp->columnNumber("signal"),firstsample,gapsample);
  }

void cut_quaternions (paramfile &par_file,
  const vector<interval<int64> > &gapsample, int64 firstsample)
  {
  string input = par_file.find<string>("quaternions_in","");
  if (input=="") return;
  string output = par_file.find<string>("quaternions_out");
  safe_ptr<iohandle> inp (HandleManager.openObject(input,"quat.LS_satpt_quat"));
  safe_ptr<iohandle> out (HandleManager.createObject
    (output,"quat.LS_satpt_quat"));

  cutColumn<double>(*inp,*out,inp->columnNumber("quat_w"),firstsample,gapsample);
  cutColumn<double>(*inp,*out,inp->columnNumber("quat_x"),firstsample,gapsample);
  cutColumn<double>(*inp,*out,inp->columnNumber("quat_y"),firstsample,gapsample);
  cutColumn<double>(*inp,*out,inp->columnNumber("quat_z"),firstsample,gapsample);
  }

void cut_detpt (paramfile &par_file,
  const vector<interval<int64> > &gapsample, int64 firstsample)
  {
  string input = par_file.find<string>("detpt_in","");
  if (input=="") return;
  string output = par_file.find<string>("detpt_out");
  safe_ptr<iohandle> inp (HandleManager.openObject(input));
  string type = inp->objectType();
  bool fullptg;
  if ((type=="pointing.LS_detpoint") || (type=="pointing.LS_detpoint_dp"))
    fullptg=false;
  else if ((type=="pointing.LS_detpoint_with_orientation") ||
     (type=="pointing.LS_detpoint_with_orientation_dp"))
    fullptg=true;
  else
    planck_fail("bad detector pointing file type '"+type+"'.");

  safe_ptr<iohandle> out (HandleManager.createObject(output,type));

  cutColumn<double>(*inp,*out,inp->columnNumber("theta"),firstsample,gapsample);
  cutColumn<double>(*inp,*out,inp->columnNumber("phi"),firstsample,gapsample);
  if (fullptg)
    cutColumn<double>(*inp,*out,inp->columnNumber("psi"),firstsample,gapsample);
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("gapmaker", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile par_file(mng.getParams());

  Sat_Info_Simmission satinfo (false, par_file);
  focalplane_db fpdb(par_file);
  string det_id = par_file.find<string>("detector_id");
  double f_samp = fpdb.getValue<double>(det_id,"f_samp");

  int first = par_file.find<int>("first_pointing",1);
  if (first==-1) first=1;
  --first;
  int last = par_file.find<int>("last_pointing",satinfo.numPeriods());
  if (last==-1) last=satinfo.numPeriods();
  --last;
  planck_assert ((first>=0)&&(last>=0)&&(first<=last)
    &&(last<satinfo.numPeriods()),
    "inconsistent first and last pointing periods!");

  int64 first_sample = satinfo.startIndex(first,f_samp);
  int64 last_sample = satinfo.startIndex(last,f_samp)
                    + satinfo.numSamples(last,f_samp);

  safe_ptr<iohandle> inp_gap (HandleManager.openObject
    (par_file.find<string> ("gapfile")));
  arr<double> tmp, tmp2;
  inp_gap->readEntireColumn("gap_start",tmp);
  inp_gap->readEntireColumn("gap_end",tmp2);
  planck_assert(tmp.size()==tmp2.size(), "ill-formed gap object");

  double startTime = satinfo.startTime(0);
  vector<interval<int64> > gapsample;
  for (tsize i=0; i<tmp.size(); ++i)
    {
    int64 start = nearest<int64> ((tmp [i]-startTime)*f_samp);
    int64 end   = nearest<int64> ((tmp2[i]-startTime)*f_samp);
    if (start<=end)
      gapsample.push_back (interval<int64>(start,end));
    }

  sort (gapsample.begin(), gapsample.end());

  unsigned int cur=0;
  while (cur<gapsample.size()-1)
    {
    if (gapsample[cur].hi>=gapsample[cur+1].lo)
      {
      if (gapsample[cur].hi<gapsample[cur+1].hi)
        gapsample[cur].hi = gapsample[cur+1].hi;
      gapsample.erase (gapsample.begin()+cur+1);
      }
    else
      ++cur;
    }

  if (gapsample.size()==0)
    cout << "Warning: gap list does not contain valid gaps!" << endl;

  write_flagfile (par_file, gapsample, first_sample, last_sample);

  cut_timestamps (par_file, gapsample, first_sample);
  cut_toi (par_file, gapsample, first_sample);
  cut_detpt (par_file, gapsample, first_sample);
  cut_quaternions (par_file, gapsample, first_sample);
PLANCK_DIAGNOSIS_END
  }
