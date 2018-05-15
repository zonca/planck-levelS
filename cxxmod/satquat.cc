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
 *  Copyright (C) 2007-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "satquat.h"
#include "rotmatrix.h"
#include "sat_info.h"
#include "arr.h"

using namespace std;

void Satquat::Get_Quaternions (const arr<double> &times, arr<quaternion> &quat)
  const
  {
  int sz=times.size();
  quat.alloc(sz);

#pragma omp parallel
{
  rotmatrix trans;
  int i;
#pragma omp for schedule (static)
  for (i=0; i<sz; ++i)
    {
    satptg.getTransform(times[i], trans);
    quat[i] = trans;
    }
}
  }
