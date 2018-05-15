/*
 *  This file is part of the Planck simulation package
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
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  This code is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  \author Martin Reinecke
 *  \author Reinhard Hell
 */

#ifndef PLANCK_BEAMCHARACTERISTIC_H
#define PLANCK_BEAMCHARACTERISTIC_H

#include <string>
#include <algorithm>
#include "arr.h"

/*
 *  written by  Dr. Reinhard M. Hell,  Max-Planck-Institute for Astrophysics
 */

class BeamCharacteristic
  {
  public:
    virtual ~BeamCharacteristic() {}

    //! Scalar and polarised values of the beam at position (theta,phi) with beam orientation psi (angles in radian)
    virtual void beamValue(double theta, double phi, double psi, double &scalar,
      double &pol) const = 0;

    virtual bool asymmetric() const { return true; }

    virtual double max_theta() const = 0;
  };

class GridBeam : public BeamCharacteristic
  {
  private:
    arr2<float> grid, gridQ, gridU;
    int ntheta, nphi;
    double maxtheta, idphi, idtheta;
    bool do_pol;

  public:
    GridBeam(const std::string &beamfile, bool polarisation);

    virtual void beamValue(double theta, double phi, double, double &scalar,
      double &pol) const;

    virtual double max_theta() const { return maxtheta; }
  };

class CartBeam: public BeamCharacteristic
  {
  private:
    arr2<float> grid, gridQ, gridU;
    int nx, ny;
    double xc, yc, dx, dy, maxtheta;
    bool do_pol;

  public:
    CartBeam(const std::string &beamfile, bool polarisation);

    virtual void beamValue(double theta, double phi, double psi, double &scalar,
      double &pol) const;

    virtual double max_theta() const { return maxtheta; }
  };

class GaussBeam: public BeamCharacteristic
  {
  private:
    double sigma;
    bool do_pol;
    double pol_angle, f1, f2;

  public:
    GaussBeam(double fwhm, bool polarisation, double polangle, double epsilon);

    virtual void beamValue(double theta, double, double psi, double &scalar,
      double &pol) const;

    virtual bool asymmetric() const { return false; }

    virtual double max_theta() const { return 6*sigma; }
  };

class EllipticGaussBeam: public BeamCharacteristic
  {
  private:
    double sigma1, sigma2, rotAngle;
    bool do_pol;
    double pol_angle, f1, f2;

  public:
    EllipticGaussBeam (double fwhm1, double fwhm2, double rotangle,
      bool polarisation, double polangle, double epsilon);

    virtual void beamValue(double theta, double phi, double psi, double &scalar,
      double &pol) const;

    virtual double max_theta() const { return 6*std::max(sigma1,sigma2); }
  };

#endif
