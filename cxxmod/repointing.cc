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
 *  Copyright (C) 2009-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "repointing.h"
#include "arr.h"
#include "paramfile.h"

using namespace std;

repointing::repointing (paramfile &params, bool nominal)
  : rng(params.find<int>("repointing_rand_seed",4711))
  {
  planck_assert(!nominal,"repointing flag does not work with nominal pointing");
  }

void repointing::get_data (const arr<double> &times, arr<bool> &flag)
  {
  flag.alloc (times.size());

  const double gap_avg   = 180.,
               gap_sigma = 20.,
               gap_min   = 120.,
               gap_max   = 240.;

  double t_rept;
  do
    { t_rept = gap_avg + gap_sigma*rng.rand_gauss(); }
  while ((t_rept<gap_min) || (t_rept>gap_max));

  for (tsize i=0; i<times.size(); ++i)
    flag[i] = (times[i]-times[0])<t_rept;
  }
