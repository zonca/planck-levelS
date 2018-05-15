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
 *  This file contains the sampler module for Planck originally written
 *  by Ian Grivell and Bob Mann. For details see the note "Sampler module
 *  for scan strategy simulations" in LiveLink.
 *
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SAMPLER_H
#define PLANCK_SAMPLER_H

#include <string>
#include "arr.h"

class focalplane_db;
class paramfile;

class Sampler
  {
  private:
    double tstep, scale, f_samp_det, f_samp1, timeshift;
    arr<double> deltat;

  public:
    Sampler (paramfile &params, focalplane_db &fpdb, const std::string &det_id,
      int n_integ, double fsamp1);
    void sample (const arr<double> &in, const arr<double> &t1,
      const arr<double> &t2, arr<double> &out) const;
  };

#endif
