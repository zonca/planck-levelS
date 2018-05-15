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
 *  Copyright (C) 2005-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm.h"
#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "iohandle_current.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"
#include "trafos.h"
#include "announce.h"
#include "string_utils.h"

using namespace std;

namespace {

Trafo maketrafo (int num)
  {
  switch (num)
    {
    case  1: return Trafo(2000,2000,Equatorial,Galactic);
    case  2: return Trafo(2000,2000,Galactic,Equatorial);
    case  3: return Trafo(2000,2000,Equatorial,Ecliptic);
    case  4: return Trafo(2000,2000,Ecliptic,Equatorial);
    case  5: return Trafo(2000,2000,Ecliptic,Galactic);
    case  6: return Trafo(2000,2000,Galactic,Ecliptic);
    case  7: return Trafo(1950,1950,Equatorial,Galactic);
    case  8: return Trafo(1950,1950,Galactic,Equatorial);
    case  9: return Trafo(1950,1950,Equatorial,Ecliptic);
    case 10: return Trafo(1950,1950,Ecliptic,Equatorial);
    case 11: return Trafo(1950,1950,Ecliptic,Galactic);
    case 12: return Trafo(1950,1950,Galactic,Ecliptic);
    default: planck_fail("Unsupported transformation "+dataToString(num));
    }
  }

} // unnamed namespace

int main(int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup("rotmap_cxx", argc==5,
    "Usage: rotmap_cxx <infile> <outfile> <itransform> <pol>\n\n"
    "Transform 1: Equatorial (2000) -> Galactic   (2000)\n"
    "          2: Galactic   (2000) -> Equatorial (2000)\n"
    "          3: Equatorial (2000) -> Ecliptic   (2000)\n"
    "          4: Ecliptic   (2000) -> Equatorial (2000)\n"
    "          5: Ecliptic   (2000) -> Galactic   (2000)\n"
    "          6: Galactic   (2000) -> Ecliptic   (2000)\n"
    "          7: Equatorial (1950) -> Galactic   (1950)\n"
    "          8: Galactic   (1950) -> Equatorial (1950)\n"
    "          9: Equatorial (1950) -> Ecliptic   (1950)\n"
    "         10: Ecliptic   (1950) -> Equatorial (1950)\n"
    "         11: Ecliptic   (1950) -> Galactic   (1950)\n"
    "         12: Galactic   (1950) -> Ecliptic   (1950)\n\n"
    "pol: T or F\n");

  iohandle_current::Manager mng ("");

  string infile  = argv[1];
  string outfile = argv[2];
  int trafo = stringToData<int>(argv[3]);
  bool polarisation = stringToData<bool>(argv[4]);

  Trafo tr(maketrafo(trafo));

  Healpix_Map<double> mapT,mapQ,mapU;
  Alm<xcomplex<double> > almT,almG,almC;

  if (!polarisation)
    {
    read_Healpix_map_from_dmc (infile, mapT);
    double avg=mapT.average();
    mapT.Add(-avg);
    int lmax = mapT.Nside()*3;
    almT.Set(lmax,lmax);
    map2alm_iter(mapT,almT,3);
    rotate_alm(almT,tr.Matrix());
    alm2map(almT,mapT);
    mapT.Add(avg);
    write_Healpix_map_to_dmc (outfile,mapT);
    }
  else
    {
    read_Healpix_map_from_dmc (infile, mapT, mapQ, mapU);
    double avg=mapT.average();
    mapT.Add(-avg);
    int lmax = mapT.Nside()*3;
    almT.Set(lmax,lmax); almG.Set(lmax,lmax); almC.Set(lmax,lmax);
    map2alm_pol_iter(mapT,mapQ,mapU,almT,almG,almC,3);
    rotate_alm(almT,almG,almC,tr.Matrix());
    alm2map_pol(almT,almG,almC,mapT,mapQ,mapU);
    mapT.Add(avg);
    write_Healpix_map_to_dmc (outfile,mapT,mapQ,mapU);
    }

PLANCK_DIAGNOSIS_END
  }
