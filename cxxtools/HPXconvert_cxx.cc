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

#include "trafos.h"
#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "iohandle_current.h"
#include "announce.h"
#include "string_utils.h"

using namespace std;

namespace {

Trafo maketrafo (int num)
  {
  switch (num)
    {
    case  2: return Trafo(2000,2000,Equatorial,Galactic);
    case  1: return Trafo(2000,2000,Galactic,Equatorial);
    case  4: return Trafo(2000,2000,Equatorial,Ecliptic);
    case  3: return Trafo(2000,2000,Ecliptic,Equatorial);
    case  6: return Trafo(2000,2000,Ecliptic,Galactic);
    case  5: return Trafo(2000,2000,Galactic,Ecliptic);
    case  8: return Trafo(1950,1950,Equatorial,Galactic);
    case  7: return Trafo(1950,1950,Galactic,Equatorial);
    case 10: return Trafo(1950,1950,Equatorial,Ecliptic);
    case  9: return Trafo(1950,1950,Ecliptic,Equatorial);
    case 12: return Trafo(1950,1950,Ecliptic,Galactic);
    case 11: return Trafo(1950,1950,Galactic,Ecliptic);
    default: planck_fail("Unsupported transformation "+dataToString(num));
    }
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup("HPXconvert_cxx", argc==5,
    "Usage: HPXconvert_cxx <infile> <outfile> <itransform> <pol>\n\n"
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

  Trafo tr (maketrafo(trafo));

  if (!polarisation)
    {
    Healpix_Map<float> inmap;
    read_Healpix_map_from_dmc(infile,inmap);
    Healpix_Map<float> outmap (inmap.Nside(), inmap.Scheme(), SET_NSIDE);

    for (int m=0; m<inmap.Npix(); ++m)
      {
      pointing pt = outmap.pix2ang(m);
      outmap[m] = inmap.interpolated_value(tr(pt));
      }
    write_Healpix_map_to_dmc(outfile, outmap);
    }
  else
    {
    Healpix_Map<float> inmapT, inmapQ, inmapU;
    read_Healpix_map_from_dmc(infile,inmapT,inmapQ,inmapU);
    Healpix_Map<float> outmapT (inmapT.Nside(), inmapT.Scheme(), SET_NSIDE),
                       outmapQ (inmapT.Nside(), inmapT.Scheme(), SET_NSIDE),
                       outmapU (inmapT.Nside(), inmapT.Scheme(), SET_NSIDE);

    for (int m=0; m<inmapT.Npix(); ++m)
      {
      pointing pt2;
      double dpsi;
      tr.rotatefull(outmapT.pix2ang(m),pt2,dpsi);
      fix_arr<int,4> pix;
      fix_arr<double,4> wgt;
      inmapT.get_interpol(pt2,pix,wgt);
      outmapT[m] = inmapT.interpolation(pix,wgt);
      float valq = inmapQ.interpolation(pix,wgt),
            valu = inmapU.interpolation(pix,wgt);
      outmapQ[m] = float(cos(2*dpsi)*valq + sin(2*dpsi)*valu);
      outmapU[m] = float(cos(2*dpsi)*valu - sin(2*dpsi)*valq);
      }
    write_Healpix_map_to_dmc(outfile, outmapT, outmapQ, outmapU);
    }
PLANCK_DIAGNOSIS_END
  }
