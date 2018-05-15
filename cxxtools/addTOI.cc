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
 *  Copyright (C) 2006-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "announce.h"
#include "iohandle_current.h"
#include "paramfile.h"
#include "share_utils.h"
#include "string_utils.h"

using namespace std;

namespace {

int find_ntoi(paramfile &params)
  {
  planck_assert (!params.param_present("toi0"),
    "Key 'toi0' found, but first component must have index 1");
  planck_assert (!params.param_present("factor0"),
    "Key 'factor0' found, but first component must have index 1");
  int res=0;
  while (true)
    {
    if (params.param_present(string("toi")+dataToString(res+1)))
      ++res;
    else
      return res;
    }
  }

void checkDoubleData (const iohandle &inp, double &value, const string &name)
  {
  if (inp.keyPresent(name))
    {
    double dtmp=inp.getKey<double>(name);
    if (value<0)
      value=dtmp;
    else
      planck_assert (approx(dtmp,value,1e-14),
        "inconsistent values for '"+name+"'");
    }
  else
    cout << "Warning: key '" << name << "' missing in input file" << endl;
  }

void assembleData (const iohandle &inp, int64 &first_sample,
  double &start_time_mission, double &start_time, double &f_samp, uint64 &nsamp)
  {
  int incol = inp.columnNumber("signal");
  uint64 utmp = inp.columnLength(incol);
  if (nsamp==0)
    nsamp=utmp;
  else
    planck_assert(utmp==nsamp,"trying to coadd TOD of different lengths");

  if (inp.keyPresent("first_sample"))
    {
    int64 itmp=inp.getKey<int64>("first_sample");
    if (first_sample<0)
      first_sample=itmp;
    else
      planck_assert (itmp==first_sample,
        "inconsistent values for 'first_sample'");
    }
  else
    cout << "Warning: key 'first_sample' missing in input file" << endl;
  checkDoubleData (inp, start_time_mission, "start_time_mission");
  checkDoubleData (inp, start_time, "start_time");
  checkDoubleData (inp, f_samp, "f_samp");
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup("addTOI",(argc==2)||(((argc&1)==0)&&(argc>=4)),
    "Usage: addTOI <init object>\n"
    "or:    addTOI <file1> <factor1> [<file2> <factor2>] [...] [output file]");

  string outfile;
  arr<string> toiname;
  arr<double> factor;
  int ntoi;

  iohandle_current::Manager mng ((argc==2) ? argv[1] : "");

  if (argc==2)
    {
    paramfile params (mng.getParams());

    outfile = params.find<string>("outfile");

    ntoi = find_ntoi(params);
    cout << "Co-adding " << ntoi << " TOIs." << endl;
    toiname.alloc(ntoi);
    factor.alloc(ntoi);
    for (int i=0; i<ntoi; ++i)
      {
      string suffix = dataToString(i+1);
      toiname[i] = params.find<string>(string("toi")+suffix);
      factor[i] = params.find<double>(string("factor")+suffix);
      }
    }
  else
    {
    ntoi = (argc-2)/2;
    toiname.alloc(ntoi);
    factor.alloc(ntoi);
    for (int i=0; i<ntoi; ++i)
      {
      toiname[i] = argv[2*i+1];
      factor[i] = stringToData<double>(argv[2*i+2]);
      }
    outfile = argv[argc-1];
    }

  arr<safe_ptr<iohandle> > inp(ntoi);
  arr<int> incol(ntoi);
  uint64 nsamp=0;
  int64 first_sample=-1;
  double start_time_mission=-1, start_time=-1, f_samp=-1;
  for (int i=0; i<ntoi; ++i)
    {
    inp[i]=HandleManager.openObject(toiname[i],"toi.LS_toi");
    incol[i] = inp[i]->columnNumber("signal");
    assembleData (*inp[i], first_sample, start_time_mission, start_time,
      f_samp, nsamp);
    }
  safe_ptr<iohandle> out (HandleManager.createObject(outfile,"toi.LS_toi"));
  int outcol=out->columnNumber("signal");
  if (first_sample>=0)
    out->setKey("first_sample",first_sample);
  if (start_time_mission>=0)
    out->setKey("start_time_mission",start_time_mission);
  if (start_time>=0)
    out->setKey("start_time",start_time);
  if (f_samp>=0)
    out->setKey("f_samp",f_samp);

  arr<float> inpdata,outdata;
  chunkMaker cm(nsamp,1024*256);
  uint64 offset,psamp;
  while (cm.getNext(offset,psamp))
    {
    inpdata.alloc(psamp);
    outdata.alloc(psamp);
    outdata.fill(0);

    for (long i=0; i<ntoi; ++i)
      {
      inp[i]->readColumn(incol[i],inpdata,offset);
      for (tsize m=0; m<inpdata.size(); ++m)
        outdata[m]+=float(inpdata[m]*factor[i]);
      }
    out->appendColumn(outcol,outdata);
    }

PLANCK_DIAGNOSIS_END
  }
