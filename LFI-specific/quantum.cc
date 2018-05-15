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
 * This file implements the LFI Onboard Processing and Ground Segment
 * simulator. Simulates the PType2 (mixing and quantization) and the
 * following reconstruction by the ground segment. Neither packaging nor
 * compression/decompression are included.
 *
 * This is based on the original work by M.Maris and M.Frailis from
 * OAT - 06/29/06.
 *
 * This is the version integrated into the LFI DPC pipeline.
 *
 * Copyright (C) 2006, LFI-DPC
 * Author: Davide Maino
 */

#include <cmath>
#include "arr.h"
#include "announce.h"
#include "paramfile.h"
#include "iohandle_current.h"

using namespace std;

namespace {

void lfi_obc (arr<int> &q1, arr<int> &q2, arr<float> &sky, arr<float> &load,
  float second_quant, float offset_adjust, float gmf1,
  float gmf2, int naver)
  {
  double xnaver=1./naver;
  for (tsize i=0;i<q1.size();i++)
    {
    double p1 = (sky[i]-gmf1*load[i])*xnaver;
    double p2 = (sky[i]-gmf2*load[i])*xnaver;
    q1[i] = nearest<int>(second_quant*(p1+offset_adjust));
    q2[i] = nearest<int>(second_quant*(p2+offset_adjust));
    }
  }

void lfi_gs (arr<float> &sky, arr<float> &load, arr<int> &q1,
  arr<int> &q2, float second_quant, float offset_adjust,
  float gmf1, float gmf2)
  {
  for (tsize i=0;i<sky.size();i++)
    {
    double delta = gmf2 - gmf1;
    double rp1 = q1[i]/second_quant-offset_adjust;
    double rp2 = q2[i]/second_quant-offset_adjust;
    sky[i] = float((gmf2*rp1 - gmf1*rp2)/delta);
    load[i] = float(-(rp2 - rp1)/delta);
    }
  }

void lfi_proc_err (arr<float> &esky, arr<float> &eload, arr<float> &rsky,
  arr<float> &rload, arr<float> &sky, arr<float> &load, int naver)
  {
  double xnaver=1./naver;
  for (tsize i=0;i<esky.size();i++)
    {
    esky[i]  = rsky[i]  - float(sky[i]*xnaver);
    eload[i] = rload[i] - float(load[i]*xnaver);
    }
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
  module_startup ("quantum", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  string fnameout = params.find<string>("fnameout");
  string fnamein = params.find<string>("fnamein");
  int n_aver = params.find<int>("n_aver",1);
  float offset = params.find<float>("offset");
  float quant = params.find<float>("quant");
  float gmf1 = params.find<float>("gmf1");
  float gmf2 = params.find<float>("gmf2");

  safe_ptr<iohandle>
    inp (HandleManager.openObject(fnamein,"toi.science.LFI_Data")),
    out (HandleManager.createObject(fnameout,"toi.science.LFI_Data"));
  // int obtcol = inp->columnNumber("sampleOBT");
  int skycol = inp->columnNumber("sky_adu");
  int refcol = inp->columnNumber("ref_adu");
  uint64 rows = inp->columnLength(skycol);

  int out_skycol = out->columnNumber("sky_adu");
  int out_refcol = out->columnNumber("ref_adu");

  // simulates obc processing
  arr<int> q1,q2;
  arr<float> fsky, fload, rfsky, rfload, efsky, efload;

  const uint64 chunksize=1024*256;
  uint64 offset_read = 0;
  while (offset_read<rows)
    {
    uint64 psamp=min(chunksize,rows-offset_read);
    q1.alloc(psamp); q2.alloc(psamp);
    fsky.alloc(psamp); fload.alloc(psamp);
    rfsky.alloc(psamp); rfload.alloc(psamp);
    efsky.alloc(psamp); efload.alloc(psamp);

    inp->readColumn(skycol,fsky,offset_read);
    inp->readColumn(refcol,fload,offset_read);

    lfi_obc(q1,q2,fsky,fload,quant,offset,gmf1,gmf2,n_aver);
    lfi_gs(rfsky, rfload, q1, q2, quant, offset, gmf1, gmf2);
    lfi_proc_err(efsky, efload, rfsky, rfload, fsky, fload, n_aver);

    out->appendColumn(out_skycol,rfsky);
    out->appendColumn(out_refcol,rfload);
    offset_read+=chunksize;
    }
  }
