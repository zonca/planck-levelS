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
 *  Copyright (C) 2011-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "wobble_correction.h"
#include "lsconstants.h"
#include "iohandle_current.h"

using namespace std;

wobble_correction::wobble_correction (paramfile &params)
  {
  safe_ptr<iohandle> inp(HandleManager.openObject
    (params.find<string>("wcorr_tilt_angles"),"sat.LS_tiltAngles"));
  inp->readEntireColumn("tilt1",psi1);
  inp->readEntireColumn("tilt2",psi2);
  planck_assert(psi1.size()==psi2.size(),"array size mismatch");
  for (tsize i=0; i<psi1.size(); ++i)
    { psi1[i]*=arcmin2rad; psi2[i]*=arcmin2rad; }
  }

rotmatrix wobble_correction::get_matrix (int period) const
  {
  planck_assert((period>=0)&&(period<int(psi1.size())), "bad pointing period");

  const double psi1Ref =   4.491 *arcmin2rad,
               psi2Ref = -28.2002*arcmin2rad;

  rotmatrix psi1R( cos(psi1Ref), sin(psi1Ref), 0,
                  -sin(psi1Ref), cos(psi1Ref), 0,
                              0,            0, 1);
  rotmatrix psi2R(cos(psi2Ref), 0, -sin(psi2Ref),
                             0, 1,             0,
                  sin(psi2Ref), 0,  cos(psi2Ref));
  rotmatrix psi1RT(cos(psi1[period]), -sin(psi1[period]), 0,
                   sin(psi1[period]),  cos(psi1[period]), 0,
                                   0,                  0, 1);
  rotmatrix R_prod1 = psi2R * psi1R;

  rotmatrix psi2RT( cos(psi2[period]), 0, sin(psi2[period]),
                                    0, 1,                 0,
                   -sin(psi2[period]), 0, cos(psi2[period]));
  rotmatrix R_prod2 = psi2RT * R_prod1;
  rotmatrix R_prod3 = psi1RT * R_prod2;
  return R_prod3;
  }
