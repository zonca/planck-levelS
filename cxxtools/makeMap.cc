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
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "announce.h"
#include "paramfile.h"
#include "iohandle_current.h"
#include "io_utils.h"
#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "string_utils.h"

using namespace std;

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup("makeMap",(argc==2)||(argc==5),
    "Usage: makeMap <init object>\n"
    "or:    makeMap <TOD file> <pointing file> <map file> <nside>");
  iohandle_current::Manager mng ((argc==2) ? argv[1] : "");

  string outfile;
  int nside;
  safe_ptr<iohandle> inp1, inp2;

  if (argc==2)
    {
    paramfile params (mng.getParams());

    outfile = params.find<string>("outfile");
    nside = params.find<int>("nside");
    inp1=HandleManager.openObject(params.find<string>("tod_file"));
    inp2=HandleManager.openObject(params.find<string>("pointing_file"));
    }
  else
    {
    outfile = argv[3];
    nside = stringToData<int>(argv[4]);
    inp1=HandleManager.openObject(argv[1]);
    inp2=HandleManager.openObject(argv[2]);
    }

  int toicol = inp1->columnNumber("signal");
  int thetacol = inp2->columnNumber("theta");
  int phicol = inp2->columnNumber("phi");

  planck_assert (inp1->columnLength(toicol)==inp2->columnLength(thetacol),
    "lengths of TOI and pointing are different");
  uint64 nsamples = inp1->columnLength(toicol);

  Healpix_Map<double> map(nside,RING,SET_NSIDE);
  map.fill(0);
  Healpix_Map<int> imap(nside,RING,SET_NSIDE);
  imap.fill(0);

  const uint64 chunksize = 1024*256;
  uint64 offset=0;
  arr<float> val;
  arr<pointing> ptg;
  while (offset<nsamples)
    {
    int sz=min(chunksize,nsamples-offset);
    val.alloc(sz);
    ptg.alloc(sz);
    readPointing(*inp2,thetacol,phicol,ptg,offset);
    inp1->readColumn(toicol,val,offset);
    for (int m=0; m<sz; ++m)
      {
      int pix = map.ang2pix(ptg[m]);
      map[pix] += val[m];
      ++imap[pix];
      }
    offset+=chunksize;
    }

  for (int m=0; m<map.Npix(); ++m)
    map[m] = (imap[m]>0) ? map[m]/imap[m] : Healpix_undef;

  write_Healpix_map_to_dmc(outfile,map);
PLANCK_DIAGNOSIS_END
  }
