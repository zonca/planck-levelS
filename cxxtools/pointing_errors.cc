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
#include "planck_rng.h"
#include "iohandle_current.h"
#include "quaternion.h"
#include "io_utils.h"

using namespace std;

namespace {

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("pointing_errors", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params(mng.getParams());

  safe_ptr<iohandle> inp (HandleManager.openObject
    (params.find<string>("infile"),"quat.LS_satpt_quat"));
  safe_ptr<iohandle> out(HandleManager.createObject
    (params.find<string>("outfile"),"quat.LS_satpt_quat"));

  planck_rng rng(params.find<int>("rand_seed"));

  double sigma = params.find<double>("sigma");

  int wcoli = inp->columnNumber("quat_w");
  int xcoli = inp->columnNumber("quat_x");
  int ycoli = inp->columnNumber("quat_y");
  int zcoli = inp->columnNumber("quat_z");

  int wcolo = out->columnNumber("quat_w");
  int xcolo = out->columnNumber("quat_x");
  int ycolo = out->columnNumber("quat_y");
  int zcolo = out->columnNumber("quat_z");

  const uint64 chunksize=1024*256;
  uint64 nsamp = inp->columnLength(wcoli);
  uint64 offset = 0;
  arr<quaternion> quat;

  while (offset<nsamp)
    {
    uint64 psamp=min(chunksize,nsamp-offset);
    quat.alloc(psamp);

    readQuaternions (*inp,wcoli,xcoli,ycoli,zcoli,quat,offset);

    for (tsize m=0; m<quat.size(); ++m)
      {
      double alpha = rng.rand_gauss()*sigma;
      double beta = rng.rand_uni()*twopi;
      alpha*=0.5;
      double salpha=sin(alpha);
      quaternion perturb (cos(alpha),salpha*cos(beta),salpha*sin(beta),0);
      quat[m] *= perturb;
      }

    appendQuaternions (*out,wcolo,xcolo,ycolo,zcolo,quat);
    offset+=chunksize;
    }

PLANCK_DIAGNOSIS_END
  }
