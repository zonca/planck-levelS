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

#ifndef PLANCK_ZLE_H
#define PLANCK_ZLE_H

#include <string>
#include "ephemerides.h"
#include "pointing.h"

template<typename T> class arr;
class paramfile;
class focalplane_db;

class ZLE
  {
  private:
    safe_ptr<ephemerides> ephSun, ephBary;
    int nside;

  public:
    ZLE (focalplane_db &focal, const std::string &det_id, paramfile &params);
    ~ZLE();
    void Add_Intensities (const arr<pointing> &detpt, const arr<vec3> &vdetpt,
      const arr<double> &times, arr<double> &intensity) const;
  };

#endif
