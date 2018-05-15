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
 *  Copyright (C) 2007-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "iohandle_current.h"
#include "pointing.h"
#include "arr.h"
#include "paramfile.h"
#include "quaternion.h"
#include "io_utils.h"
#include "announce.h"
#include "share_utils.h"

using namespace std;

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("sat2quat", argc, argv);

  iohandle_current::Manager mng (argc, argv);
  paramfile params(mng.getParams());

  safe_ptr<iohandle> inp (HandleManager.openObject
    (params.find<string>("infile"),"satelliteinfo.LS_satinfo"));
  safe_ptr<iohandle> out(HandleManager.createObject
    (params.find<string>("outfile"),"quat.LS_satpt_quat"));

  int thetaxcol = inp->columnNumber("theta_x");
  int phixcol = inp->columnNumber("phi_x");
  int thetaycol = inp->columnNumber("theta_y");
  int phiycol = inp->columnNumber("phi_y");
  int thetazcol = inp->columnNumber("theta_z");
  int phizcol = inp->columnNumber("phi_z");

  int quatwcol = out->columnNumber("quat_w");
  int quatxcol = out->columnNumber("quat_x");
  int quatycol = out->columnNumber("quat_y");
  int quatzcol = out->columnNumber("quat_z");

  arr<pointing> ptx,pty,ptz;
  arr<quaternion> quat;

  chunkMaker cm(inp->columnLength(thetaxcol),1024*256);
  uint64 offset, psamp;
  arr<double> tmp;
  while (cm.getNext(offset,psamp))
    {
    ptx.alloc(psamp);
    pty.alloc(psamp);
    ptz.alloc(psamp);
    quat.alloc(psamp);

    readPointing (*inp,thetaxcol,phixcol,thetaycol,phiycol,thetazcol,phizcol,
      ptx,pty,ptz,offset);

    for (tsize m=0; m<ptx.size(); ++m)
      quat[m] = quaternion (rotmatrix(ptx[m],pty[m],ptz[m]));

    appendQuaternions (*out,quatwcol,quatxcol,quatycol,quatzcol,quat);
    }

PLANCK_DIAGNOSIS_END
  }
