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

#ifndef PLANCK_DIPOLE_H
#define PLANCK_DIPOLE_H

#include <string>
#include "vec3.h"
template<typename T> class arr;
class paramfile;
class Sat_Info;
class focalplane_db;

class Dipole
  {
  private:
    enum speedtype { SOLSYS, SATELLITE, TOTAL };

    Sat_Info &satinfo;
    focalplane_db &fpdb;
    bool thermotemp, do_dipole, do_fsldp;
    int outputtype;
    vec3 solsysdir_v;
    speedtype speed;
    double dip_norm;

  public:
    Dipole (Sat_Info &info, focalplane_db &focal, paramfile &params,
            bool source_dipole, bool source_fsldp);
    void Add_Intensities (const std::string &det_id, const arr<vec3> &vdetpt,
      const arr<vec3> &vsldppt, arr<double> &intensity) const;
    void Add_Intensities (const std::string &det_id, const arr<vec3> &vdetpt,
      arr<double> &intensity) const;
  };

#endif
