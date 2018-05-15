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

#ifndef PLANCK_DETECTOR_POINTING_H
#define PLANCK_DETECTOR_POINTING_H

#include <string>
#include "vec3.h"
#include "safe_ptr.h"
#include "wobble_correction.h"
#include "ptcor.h"

template <typename T> class arr;
class pointing;
class Sat_Info;
class focalplane_db;
class paramfile;

class Detector_Pointing
  {
  private:
    Sat_Info &satptg;
    vec3 direction, xdir, vsldp;
    bool aberration;
    safe_ptr<wobble_correction> wcorr;
    safe_ptr<ptcor> ptcorr;

    void Get_Pointings (const arr<double> &times,
      arr<pointing> *detpt, arr<vec3> *vdetpt, arr<double> *heading,
      arr<vec3> *vsldppt) const;

  public:
    Detector_Pointing (paramfile &params, Sat_Info &sptg,
      const focalplane_db &fpdb, const std::string &det_id, bool aberration_);
    void Get_Pointings (const arr<double> &times, arr<pointing> &detpt,
      arr<vec3> &vdetpt, arr<double> &heading, arr<vec3> &vsldppt) const
      { Get_Pointings (times, &detpt, &vdetpt, &heading, &vsldppt); }
    void Get_Pointings (const arr<double> &times, arr<pointing> &detpt,
      arr<vec3> &vdetpt, arr<double> &heading) const
      { Get_Pointings (times, &detpt, &vdetpt, &heading, 0); }
    void Get_Pointings (const arr<double> &times, arr<pointing> &detpt,
      arr<vec3> &vdetpt) const
      { Get_Pointings (times, &detpt, &vdetpt, 0, 0); }
    void Get_Pointings (const arr<double> &times, arr<pointing> &detpt,
      arr<double> &heading) const
      { Get_Pointings (times, &detpt, 0, &heading, 0); }
    void Get_Pointings (const arr<double> &times, arr<vec3> &vdetpt) const
      { Get_Pointings (times, 0, &vdetpt, 0, 0); }
  };

#endif
