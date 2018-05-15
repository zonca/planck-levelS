/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm.h"
#include "alm_dmcio.h"
#include "iohandle_current.h"
#include "announce.h"
#include "string_utils.h"
#include "xcomplex.h"
#include "lsconstants.h"

using namespace std;

namespace {

template<typename T> void rotalm_phi (paramfile &params, bool dp)
  {
  bool polarisation = params.template find<bool>("polarisation");
  string infile = params.template find<string>("infile");
  string outfile = params.template find<string>("outfile");
  double dphi = params.template find<double>("dphi_degrees")*degr2rad;
  int lmax, mmax;

  if (polarisation)
    {
    Alm<xcomplex<T> > almT,almG,almC;
    get_almsize_dmc_pol(infile, lmax, mmax, dp);
    read_Alm_from_dmc (infile, almT, almG, almC, lmax, mmax, dp);
    arr<xcomplex<T> > fact(mmax+1);
    for (tsize m=0; m<fact.size(); ++m)
      fact[m] = xcomplex<T>(cos(m*dphi),-sin(m*dphi));
    almT.ScaleM(fact); almG.ScaleM(fact); almC.ScaleM(fact);
    write_Alm_to_dmc (outfile, almT, almG, almC, lmax, mmax, dp);
    }
  else
    {
    Alm<xcomplex<T> > alm;
    get_almsize_dmc(infile, lmax, mmax, dp);
    read_Alm_from_dmc (infile, alm, lmax, mmax, dp);
    arr<xcomplex<T> > fact(mmax+1);
    for (tsize m=0; m<fact.size(); ++m)
      fact[m] = xcomplex<T>(cos(m*dphi),-sin(m*dphi));
    alm.ScaleM(fact);
    write_Alm_to_dmc (outfile, alm, lmax, mmax, dp);
    }
  }

} // unnamed namespace

int main(int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup("rotalm_phi",argc,argv);

  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  bool dp = params.find<bool> ("double_precision",false);
  dp ? rotalm_phi<double>(params,dp) : rotalm_phi<float>(params,dp);

PLANCK_DIAGNOSIS_END
  }
