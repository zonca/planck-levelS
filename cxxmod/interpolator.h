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
 *  Copyright (C) 2003, 2004 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_INTERPOLATOR_H
#define PLANCK_INTERPOLATOR_H

#include <string>
#include "arr.h"
#include "trafos.h"

class pointing;
class paramfile;

class Interpolator
  {
  private:
    enum outtype { SIGNAL, I, Q, U };
    enum { max_order=19 };

    double inv_delta_phi, inv_delta_theta, phioffset, thetaoffset, psi;
    arr<arr2<float> > sky;
    arr<double> base_wgt;
    arr<int> psiarr;
    int nphi, ntheta, npsi, beammmax;
    int npoints, ioffset;
    outtype output;
    bool galactic;
    Trafo ecl2gal;

    void init (const std::string &ringset, int order);

    inline void weight_n (double x, double *wgt) const;

    inline void interpol_n
      (double theta, double phi, arr<double> &result) const;
    inline double interpol_psi(double omega, const arr<double> &karr) const;

  public:
    Interpolator (paramfile &params, double psi_pol=0);
    Interpolator (const std::string &ringset, double psi_pol=0,
      int order=1);
    void Get_Intensities (const arr<pointing> &detpt,
      arr<double> &ii, arr<double> &iq, arr<double> &iu) const;
    void Add_Intensities (const arr<pointing> &detpt,
      const arr<double> &heading, arr<double> &intensity) const;
  };

#endif
