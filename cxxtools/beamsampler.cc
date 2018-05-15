/*
 *  This file is part of the Planck simulation package.
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
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  This code is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2007-2013 Max-Planck-Society
 *  Authors: Martin Reinecke, Mark Ashdown
 */

#include "alm.h"
#include "alm_dmcio.h"
#include "iohandle_current.h"
#include "alm_powspec_tools.h"
#include "focalplane_db.h"
#include "lsconstants.h"
#include "announce.h"

using namespace std;

int main(int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("beamsampler", argc, argv);

  iohandle_current::Manager mng (argc, argv);
  paramfile par_file(mng.getParams());

  string infile  = par_file.find<string>("beam_in");
  string outfile = par_file.find<string>("beam_out");
  string det_id = par_file.find<string>("detector_id");
  focalplane_db fpdb(par_file);
  double theta_b = fpdb.theta_b()*degr2rad;
  double tau_int = fpdb.getValue<double>(det_id,"tau_int");
  double tau_bol = fpdb.getValue<double>(det_id,"tau_bol");
  double f_samp = fpdb.getValue<double>(det_id,"f_samp");
  double phi_uv = fpdb.getValue<double>(det_id,"phi_uv")*degr2rad;
  double theta_uv = fpdb.getValue<double>(det_id,"theta_uv")*degr2rad;
  double psi_uv = fpdb.getValue<double>(det_id,"psi_uv")*degr2rad;
  int nread = fpdb.getValue<int>(det_id,"nread");
  // if time response is instantaneous, don't need to integrate against it
  int n_integ = 1;
  // otherwise get the number of integration points
  if (tau_bol>=1e-4)
    n_integ = par_file.find<int> ("n_integ", 5);
  int n_used = 1 + int (tau_int*f_samp*nread);
  double ang_vel = twopi/60*par_file.find<double>("sat_rpm",1.);
  int mmax_out = par_file.find<int>("mmax_out");

  arr<double> toffset(n_integ*n_used);
  for (int i=0; i<n_used; ++i)
    for (int j=0; j<n_integ; ++j)
      toffset[j+i*n_integ] = ang_vel *
        (i/(f_samp*nread) - tau_bol*log((double(j+1))/n_integ));
  toffset.sort();

  int lmax, dummy;
  get_almsize_dmc_pol (infile,lmax,dummy);
  Alm<xcomplex<double> > almT,almG,almC;
  read_Alm_from_dmc (infile, almT, almG, almC, lmax, lmax); // sic

  double v_xdet = cos(theta_b)*sin(theta_uv)*sin(psi_uv-phi_uv) -
    sin(theta_b) * (cos(theta_uv)*cos(phi_uv)*sin(psi_uv-phi_uv) +
    sin(phi_uv)*cos(psi_uv-phi_uv));
  double v_ydet = cos(theta_b)*sin(theta_uv)*cos(psi_uv-phi_uv) -
    sin(theta_b) * (cos(theta_uv)*cos(phi_uv)*cos(psi_uv-phi_uv) -
    sin(phi_uv)*sin(psi_uv-phi_uv));

  double azimuth = atan2(v_ydet, v_xdet);

  Alm<xcomplex<double> > almT2=almT, almG2=almG, almC2=almC;
  for (tsize i=1; i<toffset.size(); ++i)
    {
    Alm<xcomplex<double> > almT3=almT, almG3=almG, almC3=almC;
    rotate_alm(almT3,almG3,almC3,-azimuth,toffset[i],azimuth);
    almT2.Add(almT3); almG2.Add(almG3); almC2.Add(almC3);
    }
  almT2.Scale(1./toffset.size());
  almG2.Scale(1./toffset.size());
  almC2.Scale(1./toffset.size());

  write_Alm_to_dmc (outfile,almT2,almG2,almC2,lmax,mmax_out);

PLANCK_DIAGNOSIS_END
  }
