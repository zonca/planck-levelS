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
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_OOFNOISE_H
#define PLANCK_OOFNOISE_H

#include <string>
#include "safe_ptr.h"

template<typename T> class arr;
class focalplane_db;
class paramfile;

class OofGenerator
  {
  public:
    virtual ~OofGenerator() {}
    virtual double nextSample() = 0;
    virtual void addVec(arr<double> &result, int num_real) = 0;
    virtual void reset() = 0;
  };

class Oofnoise
  {
  private:
    safe_ptr<OofGenerator> generator;
    int num_real;
    double inv_num_real;
    int stationary_periods, curperiod;

  public:
    Oofnoise (focalplane_db &fpdb, const std::string &det_id,
      int num_realisations, paramfile &params);

    void Get_Noise (int period, int num_samp, arr<double> &noise);
    void Skip (int nskip);
  };

#endif
