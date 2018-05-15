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

#include "BeamCharacteristic.h"
#include "pointing.h"
#include "iohandle.h"
#include "lsconstants.h"
#include <algorithm>
#include <cmath>

using namespace std;

namespace {

template<typename T> inline void stokesrotate (T &q, T &u, double phi)
  {
  double c2=cos(2*phi), s2=sin(2*phi);
  double q2=q, u2=u;
  q = T(-c2*q2 - s2*u2);
  u = T(s2*q2 - c2*u2);
  }

}

GridBeam::GridBeam(const string &beamfile, bool polarisation)
  : do_pol(polarisation)
  {
  // Read in values of the beam:
  safe_ptr<iohandle> inp (HandleManager.openObject(beamfile,
    do_pol ? "beam.LS_beammap_pol" : "beam.LS_beammap"));

  inp->getKey("Nphi",nphi);
  inp->getKey("Ntheta",ntheta);
  inp->getKey("Maxtheta",maxtheta);
  double dtheta = maxtheta/(ntheta-1);
  double dphi = twopi/nphi;
  idtheta = 1./dtheta;
  idphi = 1./dphi;

  grid.alloc(ntheta,nphi);
  inp->readColumnRaw("Beamdata", &(grid[0][0]), grid.size());
  if (do_pol)
    {
    gridQ.alloc(ntheta,nphi);
    gridU.alloc(ntheta,nphi);
    inp->readColumnRaw("BeamdataQ", &(gridQ[0][0]), gridQ.size());
    inp->readColumnRaw("BeamdataU", &(gridU[0][0]), gridU.size());
    for (int ith=0; ith<ntheta; ++ith)
      for (int iph=0; iph<nphi; ++iph)
        stokesrotate(gridQ[ith][iph],gridU[ith][iph],iph*dphi);
    }
  }

void GridBeam::beamValue(double theta, double phi, double psi,
  double &scalar, double &pol) const
  {
  if (theta>=maxtheta) { scalar = pol = 0; return; }
  // theta,phi: position of the point source wrt the center of the beam
  // psi: rotation of the beam reference frame
  psi=-psi;
  psi-=phi;
  phi+=pi;
  double fphi = phi*idphi;
  int iphi = int(fphi);
  fphi -= iphi;
  iphi = imodulo(iphi,nphi);
  double ftheta = theta*idtheta;
  int itheta = min(ntheta-1,max(0,int(ftheta)));
  ftheta -= itheta;
  int iphi1=iphi+1;
  if (iphi1==nphi) iphi1=0;
  double w00 = (1-ftheta)*(1-fphi),
         w01 = (1-ftheta)*(  fphi),
         w10 = (  ftheta)*(1-fphi),
         w11 = (  ftheta)*(  fphi);
  scalar =  w00 * grid[itheta  ][iphi ]
          + w01 * grid[itheta  ][iphi1]
          + w10 * grid[itheta+1][iphi ]
          + w11 * grid[itheta+1][iphi1];
  if (do_pol)
    {
    double qbeam, ubeam;
    qbeam = w00 * gridQ[itheta  ][iphi ]
          + w01 * gridQ[itheta  ][iphi1]
          + w10 * gridQ[itheta+1][iphi ]
          + w11 * gridQ[itheta+1][iphi1];
    ubeam = w00 * gridU[itheta  ][iphi ]
          + w01 * gridU[itheta  ][iphi1]
          + w10 * gridU[itheta+1][iphi ]
          + w11 * gridU[itheta+1][iphi1];
    pol = qbeam * cos(2.*psi) + ubeam * sin(2.*psi);
    }
  else
    pol = 0;
  }

CartBeam::CartBeam(const string &beamfile, bool polarisation)
  : do_pol(polarisation)
  {
  // Read in values of the beam:
  safe_ptr<iohandle> inp (HandleManager.openObject(beamfile,
    do_pol ? "beam.LS_cart_beammap_pol" : "beam.LS_cart_beammap"));

  inp->getKey("Nx",nx);
  inp->getKey("Ny",ny);
  inp->getKey("Xcentre",xc);
  inp->getKey("Ycentre",yc);
  inp->getKey("Xdelta",dx);
  inp->getKey("Ydelta",dy);
  planck_assert(nx==ny,"CartBeam: nx!=ny");
  planck_assert(approx(dx,dy,1e-12),"CartBeam: grid spacings differ");
  maxtheta=sqrt(2.)*0.5*(nx-1)*dx;
cout << "maxtheta: " << maxtheta*rad2degr << endl;
  grid.alloc(ny,nx);
  inp->readColumnRaw("Beamdata", &(grid[0][0]), grid.size());
  if (do_pol)
    {
    gridQ.alloc(ny,nx);
    gridU.alloc(ny,nx);
    inp->readColumnRaw("BeamdataQ", &(gridQ[0][0]), gridQ.size());
    inp->readColumnRaw("BeamdataU", &(gridU[0][0]), gridU.size());
    for (int y=0; y<ny; ++y)
      for (int x=0; x<nx; ++x)
        stokesrotate(gridQ[y][x],gridU[y][x],
                     safe_atan2((y-.5*(ny-1))*dy,(x-.5*(nx-1))*dx));
    }
  }

void CartBeam::beamValue(double theta, double phi, double psi,
  double &scalar, double &pol) const
  {
  psi=-psi;
  psi-=phi;
  phi+=pi;
  // Determine position numbers of point source in beam array
  double xpos=theta*cos(phi)/dx + (nx-1)*0.5;
  double ypos=theta*sin(phi)/dy + (ny-1)*0.5;
  if ((xpos<0) || (xpos>nx-1) || (ypos<0) || (ypos>ny-1))
    { scalar = pol = 0; return; }
  int xposfloor=int(xpos);
  int xposceil=xposfloor+1;
  int yposfloor=int(ypos);
  int yposceil=yposfloor+1;
  double xdelta=xpos-xposfloor;
  double ydelta=ypos-yposfloor;
  double w00 = (1-ydelta)*(1-xdelta),
         w01 = (1-ydelta)*   xdelta ,
         w10 =    ydelta *(1-xdelta),
         w11 =    ydelta *   xdelta;
  // Interpolate in beam array
  scalar =
    grid[yposfloor][xposfloor] * w00 +
    grid[yposfloor][xposceil ] * w01 +
    grid[yposceil ][xposfloor] * w10 +
    grid[yposceil ][xposceil ] * w11;
  if (do_pol)
    {
    double qbeam, ubeam;
    qbeam = w00 * gridQ[yposfloor][xposfloor]
          + w01 * gridQ[yposfloor][xposceil ]
          + w10 * gridQ[yposceil ][xposfloor]
          + w11 * gridQ[yposceil ][xposceil ];
    ubeam = w00 * gridU[yposfloor][xposfloor]
          + w01 * gridU[yposfloor][xposceil ]
          + w10 * gridU[yposceil ][xposfloor]
          + w11 * gridU[yposceil ][xposceil ];
    pol = qbeam * cos(2.*psi) + ubeam * sin(2.*psi);
    }
  else
    pol = 0;
  }

GaussBeam::GaussBeam(double fwhm, bool polarisation, double polangle,
  double epsilon)
  : do_pol(polarisation), pol_angle(polangle)
  {
  sigma=fwhm*fwhm2sigma;
  double beamOmega=2.0*pi*sigma*sigma;
  f1=(0.5*(1+epsilon))/beamOmega;
  f2=(0.5*(1-epsilon))/beamOmega;
  }

void GaussBeam::beamValue(double theta, double, double psi,
  double &scalar, double &pol) const
  {
  double thesi=theta/sigma;
  double tmp=exp(-0.5*thesi*thesi);
  scalar = f1*tmp;
  if (do_pol)
    pol = f2*tmp*cos(2*(psi-pol_angle));
  else
    pol = 0;
  }

EllipticGaussBeam::EllipticGaussBeam(double fwhm1, double fwhm2,
  double rotangle, bool polarisation, double polangle, double epsilon)
  : rotAngle(rotangle), do_pol(polarisation), pol_angle(polangle)
  {
  sigma1=fwhm1*fwhm2sigma;
  sigma2=fwhm2*fwhm2sigma;
  double beamOmega=2.0*pi*sigma1*sigma2;
  f1=(0.5*(1+epsilon))/beamOmega;
  f2=(0.5*(1-epsilon))/beamOmega;
  }

void EllipticGaussBeam::beamValue(double theta, double phi, double psi,
  double &scalar, double &pol) const
  {
  phi-=rotAngle;
  double tmp1=theta*cos(phi)/sigma1;
  double tmp2=theta*sin(phi)/sigma2;
  double tmp=exp(-0.5*(tmp1*tmp1+tmp2*tmp2));
  scalar = f1*tmp;
  if (do_pol)
    pol = f2*tmp*cos(2*(psi-pol_angle));
  else
    pol = 0;
  }
