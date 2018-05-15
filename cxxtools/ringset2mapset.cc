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
 *  Copyright (C) 2012-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "interpolator.h"
#include "healpix_base.h"
#include "iohandle_current.h"
#include "pointing.h"
#include "arr.h"
#include "paramfile.h"
#include "lsconstants.h"
#include "share_utils.h"
#include "announce.h"
#include "string_utils.h"

using namespace std;

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("ringset2mapset",argc,argv);
  iohandle_current::Manager mng (argc,argv);

  paramfile params(mng.getParams());
  tsize npsi = params.find<int>("npsi");
  string ringset = params.find<string>("ringset");
  int nside = params.find<int>("nside");
  string mapname = params.find<string>("mapname");
  int iorder = params.find<int>("interpol_order",1);

  for (tsize ipsi=0; ipsi<npsi; ++ipsi)
    {
    double psi=ipsi*(twopi/npsi);
    Interpolator interpol(ringset,psi,iorder);
    Healpix_Base hpbase (nside,RING,SET_NSIDE);

    safe_ptr<iohandle> out (HandleManager.createObject
      (mapname+intToString(ipsi,3),"map.LS_map"));
    int c1 = out->columnNumber("I_Stokes");
    out->setKey("Nside",nside);
    out->setKey("Ordering",string("RING"));

    arr<pointing> ptg;
    arr<double> hdg, res;

    chunkMaker cm(hpbase.Npix(),out->efficientChunkSize(c1));
    uint64 offset, ppix;
    while (cm.getNext(offset,ppix))
      {
      ptg.alloc(ppix);
      hdg.allocAndFill(ppix,psi);
      res.allocAndFill(ppix,0.);
      for (tsize m=0; m<ppix; ++m)
        ptg[m] = hpbase.pix2ang(m+offset);
      interpol.Add_Intensities (ptg, hdg, res);
      out->appendColumn(c1,res);
      }
    }
PLANCK_DIAGNOSIS_END
  }
