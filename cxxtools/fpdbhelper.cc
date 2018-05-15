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
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "iohandle_current.h"
#include "focalplane_db.h"
#include "lsconstants.h"

using namespace std;

namespace {

void usage()
  {
  cerr <<
    "Usage: fpdbhelper <database> <detector> <quantity>\n\n"
    "<database>: name of the focalplane data base file\n\n"
    "<detector>: name of the required detector\n\n"
    "<quantity>: determines the output of the program\n"
    "  fwhm_arcmin   : FWHM in arc minutes\n"
    "  freq_GHz      : central frequency in GHz\n"
    "                  (rounded to the nearest integer)\n"
    "  lmax          : a hint at a suitable lmax parameter for this detector\n"
    "  ringres       : number of samples taken per minute\n"
    "                  (rounded to the nearest integer)\n"
    "  f_samp        : sampling frequency (in Hz)\n"
    "  alpha         : slope of the noise spectrum (as a positive number)\n"
    "  f_knee        : knee frequency (in Hz)\n"
    "  f_min         : minimum noise frequency (in Hz)\n"
    "  net_rj        : noise equivalent antenna temperature (in K/sqrt(Hz))\n"
    "  tau_bol       : bolometer time constant (in s)\n"
    "  tau_int       : integration time for each sample (in s)\n"
    "  ellipticity   : (max FWHM)/(min FWHM)\n"
    "  theta_uv_deg  : the theta_uv angle (in degrees)\n"
    "  phi_uv_deg    : the phi_uv angle (in degrees)\n"
    "  psi_uv_deg    : the psi_uv angle (in degrees)\n"
    "  psi_pol_deg   : the psi_pol angle (in degrees)\n"
    "  psi_ell_deg   : the psi_ell angle (in degrees)\n\n";
  planck_fail_quietly("Incorrect usage");
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  if (argc!=4) usage();
  iohandle_current::Manager mng ("");

  focalplane_db fpdb (argv[1]);
  string det = argv[2];
  string quantity = argv[3];
  if (quantity == "fwhm_arcmin")
    { cout << (fpdb.getValue<double>(det,"beamfwhm")*60); return 0; }
  if (quantity == "freq_GHz")
    { cout << (nearest<int>(fpdb.getValue<double>(det,"nu_cen")*1e-9));return 0; }
  if (quantity == "lmax")
    {
    double sigma = fpdb.getValue<double>(det,"beamfwhm")*degr2rad*fwhm2sigma;
    cout << min (4096, 2*nearest<int>(3/sigma));
    return 0;
    }
  if (quantity == "ringres")
    { cout << nearest<int>(60*fpdb.getValue<double>(det,"f_samp")); return 0; }
  if (quantity == "f_samp")
    { cout << fpdb.getValue<double>(det,"f_samp"); return 0; }
  if (quantity == "alpha")
    { cout << fpdb.getValue<double>(det,"alpha"); return 0; }
  if (quantity == "f_knee")
    { cout << fpdb.getValue<double>(det,"f_knee"); return 0; }
  if (quantity == "f_min")
    { cout << fpdb.getValue<double>(det,"f_min"); return 0; }
  if (quantity == "net_rj")
    { cout << fpdb.getValue<double>(det,"net_rj"); return 0; }
  if (quantity == "tau_bol")
    { cout << fpdb.getValue<double>(det,"tau_bol"); return 0; }
  if (quantity == "tau_int")
    { cout << fpdb.getValue<double>(det,"tau_int"); return 0; }
  if (quantity == "ellipticity")
    { cout << fpdb.getValue<double>(det,"ellipticity"); return 0; }
  if (quantity == "theta_uv_deg")
    { cout << fpdb.getValue<double>(det,"theta_uv"); return 0; }
  if (quantity == "phi_uv_deg")
    { cout << fpdb.getValue<double>(det,"phi_uv"); return 0; }
  if (quantity == "psi_uv_deg")
    { cout << fpdb.getValue<double>(det,"psi_uv"); return 0; }
  if (quantity == "psi_pol_deg")
    { cout << fpdb.getValue<double>(det,"psi_pol"); return 0; }
  if (quantity == "psi_ell_deg")
    { cout << fpdb.getValue<double>(det,"psi_ell"); return 0; }
  usage();
PLANCK_DIAGNOSIS_END
  }
