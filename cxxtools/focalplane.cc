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

#include <sstream>
#include <fstream>
#include "paramfile.h"
#include "iohandle_current.h"
#include "announce.h"
#include "lsconstants.h"
#include "string_utils.h"

using namespace std;

int main(int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("focalplane", argc, argv);

  //read the parameter file:
  iohandle_current::Manager mng (argc, argv);
  paramfile par_file (mng.getParams());

  ifstream inp (par_file.find<string>("focalplane").c_str());
  planck_assert(inp, "input file could not be opened");

  //initialise the output file:
  safe_ptr<iohandle> out (HandleManager.createObject
    (par_file.find<string>("focalplane_db"),"focalplane.LS_focalplanedb"));
  out->setKey("ORIGIN",string("LPAC"));
  out->setKey("TELESCOP",string("PLANCK"));
  out->setKey("theta_b",par_file.find<double>("theta_b"));

  // find out the version of the focalplane database being used
  bool version_found = false;
  int ndet = 0;
  const int ndmax=100;
  vector<double> phi_uv(ndmax),theta_uv(ndmax),psi_uv(ndmax),psi_pol(ndmax),
    centre_freq(ndmax),min_freq(ndmax),max_freq(ndmax),f_knee(ndmax),
    alpha(ndmax),f_min(ndmax),sampling_freq(ndmax),tau_bol(ndmax),
    tau_integ(ndmax),beam_fwhm(ndmax),ellipt(ndmax),psi_ell(ndmax),
    net_rj(ndmax),epsilon(ndmax),sldp_x(ndmax),sldp_y(ndmax),sldp_z(ndmax);
  vector<int> nread(ndmax);
  vector<string> det_id(ndmax);
  while (inp)
    {
    string line;
    getline(inp,line);
    string::size_type pos = line.find("#");
    if (pos != string::npos) line = line.substr(0,pos);
    line=trim(line);
    if (line!="")
      {
      pos = line.find("=");
      if (pos != string::npos) // (equals sign found)
        {
        string par_name = trim(line.substr(0,pos));
        if (par_name=="Version")
          {
          string dbversion=trim(line.substr(pos+1));
          planck_assert
            ((dbversion.length()>=2)&&(dbversion.substr(0,2)=="6."),
            "incompatible DB version: need version 6.x");
          version_found=true;
          out->setKey("fpdb_version",dbversion);
          }
        }
      else // detector line found
        {
        planck_assert(ndet<ndmax,"too many detectors found");
        // reading focal plane positions for each detector:
        istringstream buf(line);
        buf >> det_id[ndet] >> phi_uv[ndet] >> theta_uv[ndet] >> psi_uv[ndet]
            >> psi_pol[ndet] >> epsilon[ndet] >> centre_freq[ndet]
            >> min_freq[ndet] >> max_freq[ndet] >> f_knee[ndet] >> alpha[ndet]
            >> f_min[ndet] >> sampling_freq[ndet] >> tau_bol[ndet]
            >> tau_integ[ndet] >> nread[ndet] >> beam_fwhm[ndet]
            >> ellipt[ndet] >> psi_ell[ndet] >> net_rj[ndet] >> sldp_x[ndet]
            >> sldp_y[ndet] >> sldp_z[ndet];

        beam_fwhm[ndet] /= 60.0; // convert to degrees

        cout << det_id[ndet] << ":" << phi_uv[ndet] << ":"
             << theta_uv[ndet] << ":" << psi_uv[ndet] << endl;
        ++ndet;
        }
      }
    }

  planck_assert (version_found, "did not find FPDB version in input file");
  planck_assert (ndet>0, "no detector lines found in input file");
  cout << ndet << " detectors found in input file" << endl;

  out->appendColumnRaw("detector",&det_id[0],ndet);
  out->appendColumnRaw("phi_uv",&phi_uv[0],ndet);
  out->appendColumnRaw("theta_uv",&theta_uv[0],ndet);
  out->appendColumnRaw("psi_uv",&psi_uv[0],ndet);
  out->appendColumnRaw("psi_pol",&psi_pol[0],ndet);
  out->appendColumnRaw("epsilon",&epsilon[0],ndet);
  out->appendColumnRaw("nu_cen",&centre_freq[0],ndet);
  out->appendColumnRaw("nu_min",&min_freq[0],ndet);
  out->appendColumnRaw("nu_max",&max_freq[0],ndet);
  out->appendColumnRaw("f_knee",&f_knee[0],ndet);
  out->appendColumnRaw("alpha",&alpha[0],ndet);
  out->appendColumnRaw("f_min",&f_min[0],ndet);
  out->appendColumnRaw("f_samp",&sampling_freq[0],ndet);
  out->appendColumnRaw("tau_bol",&tau_bol[0],ndet);
  out->appendColumnRaw("tau_int",&tau_integ[0],ndet);
  out->appendColumnRaw("nread",&nread[0],ndet);
  out->appendColumnRaw("beamfwhm",&beam_fwhm[0],ndet);
  out->appendColumnRaw("ellipticity",&ellipt[0],ndet);
  out->appendColumnRaw("psi_ell",&psi_ell[0],ndet);
  out->appendColumnRaw("net_rj",&net_rj[0],ndet);
  out->appendColumnRaw("sldp_x",&sldp_x[0],ndet);
  out->appendColumnRaw("sldp_y",&sldp_y[0],ndet);
  out->appendColumnRaw("sldp_z",&sldp_z[0],ndet);

PLANCK_DIAGNOSIS_END
    }
