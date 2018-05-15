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
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "healpix_map_dmcio.h"
#include "healpix_map.h"
#include "iohandle.h"
#include "arr.h"
#include "share_utils.h"
#include "string_utils.h"

using namespace std;

namespace {

string mapType (bool pol, bool dp)
  {
  string res = "map.LS_map";
  if (pol) res += "_pol";
  if (dp) res += "_dp";
  return res;
  }

void checkNsideNpix (int64 nside, int64 npix)
  {
  planck_assert (npix==12*nside*nside,
    string("mismatch between number of map pixels ("
    +dataToString(npix)+") and Nside ("+dataToString(nside)+")"));
  }

} // unnamed namespace

template<typename T> void read_Healpix_map_from_dmc
  (iohandle &inp, Healpix_Map<T> &map, const string &colname)
  {
  arr<T> myarr;
  inp.readEntireColumn (colname, myarr);
  checkNsideNpix (inp.getKey<int>("Nside"), myarr.size());
  map.Set (myarr, string2HealpixScheme(inp.getKey<string>("Ordering")));
  }

template void read_Healpix_map_from_dmc
  (iohandle &inp, Healpix_Map<float> &map, const string &colname);
template void read_Healpix_map_from_dmc
  (iohandle &inp, Healpix_Map<double> &map, const string &colname);
template void read_Healpix_map_from_dmc
  (iohandle &inp, Healpix_Map<int> &map, const string &colname);

template<typename T> void read_Healpix_map_from_dmc
  (const string &name, Healpix_Map<T> &map, const string &colname, bool dp)
  {
  string type = (colname=="Hits") ? "map.LS_hitmap" :
                mapType(colname!="I_Stokes",dp);

  safe_ptr<iohandle> inp (HandleManager.openObject(name,type));
  read_Healpix_map_from_dmc (*inp,map,colname);
  }

template void read_Healpix_map_from_dmc (const string &name,
  Healpix_Map<float> &map, const string &colname, bool dp);
template void read_Healpix_map_from_dmc (const string &name,
  Healpix_Map<double> &map, const string &colname, bool dp);

template<typename T> void read_Healpix_map_from_dmc
  (const string &name, Healpix_Map<T> &map, bool dp)
  { read_Healpix_map_from_dmc (name,map,"I_Stokes",dp); }

template void read_Healpix_map_from_dmc (const string &name,
  Healpix_Map<float> &map, bool dp);
template void read_Healpix_map_from_dmc (const string &name,
  Healpix_Map<double> &map, bool dp);

template<typename T> void read_Healpix_map_from_dmc
  (const string &name, Healpix_Map<T> &mapT, Healpix_Map<T> &mapQ,
   Healpix_Map<T> &mapU, bool dp)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject(name,mapType(true,dp)));
  int icol=inp->columnNumber("I_Stokes"),
      qcol=inp->columnNumber("Q_Stokes"),
      ucol=inp->columnNumber("U_Stokes");

  tsize npix=inp->columnLength(icol);
  planck_assert ((npix==inp->columnLength(qcol)) &&
                 (npix==inp->columnLength(ucol)), "column length mismatch");
  tsize nside=inp->getKey<int>("Nside");
  checkNsideNpix (nside, npix);
  Healpix_Ordering_Scheme scheme
    = string2HealpixScheme(inp->getKey<string>("Ordering"));
  mapT.SetNside(nside,scheme);
  mapQ.SetNside(nside,scheme);
  mapU.SetNside(nside,scheme);
  chunkMaker cm (npix,inp->efficientChunkSize(icol));
  uint64 offset, ppix;
  while (cm.getNext (offset,ppix))
    {
    inp->readColumnRaw(icol,&mapT[offset],ppix,offset);
    inp->readColumnRaw(qcol,&mapQ[offset],ppix,offset);
    inp->readColumnRaw(ucol,&mapU[offset],ppix,offset);
    }
  }

template void read_Healpix_map_from_dmc
  (const string &name, Healpix_Map<float> &mapT, Healpix_Map<float> &mapQ,
   Healpix_Map<float> &mapU, bool dp);
template void read_Healpix_map_from_dmc
  (const string &name, Healpix_Map<double> &mapT, Healpix_Map<double> &mapQ,
   Healpix_Map<double> &mapU, bool dp);

namespace {

void write_Healpix_keys (iohandle &out, const Healpix_Base &base)
  {
  string ordering = (base.Scheme()==RING) ? "RING" : "NESTED";
  out.setKey ("PIXTYPE",string("HEALPIX"));
  out.setKey ("Ordering",ordering);
  out.setKey ("Nside",base.Nside());
  out.setKey ("FIRSTPIX",0);
  out.setKey ("LASTPIX",base.Npix()-1);
  out.setKey ("INDXSCHM",string("IMPLICIT"));
  }

}

template<typename T> void write_Healpix_map_to_dmc
  (const string &name, const Healpix_Map<T> &map, bool dp)
  {
  safe_ptr<iohandle> out (HandleManager.createObject(name,mapType(false,dp)));
  write_Healpix_keys (*out, map);
  out->appendColumn("I_Stokes",map.Map());
  }

template void write_Healpix_map_to_dmc
  (const string &name, const Healpix_Map<float> &map, bool dp);
template void write_Healpix_map_to_dmc
  (const string &name, const Healpix_Map<double> &map, bool dp);

template<typename T> void write_Healpix_hitmap_to_dmc
  (const string &name, const Healpix_Map<T> &map)
  {
  safe_ptr<iohandle> out (HandleManager.createObject(name,"map.LS_hitmap"));
  write_Healpix_keys (*out, map);
  out->appendColumn("Hits",map.Map());
  }

template void write_Healpix_hitmap_to_dmc
  (const string &name, const Healpix_Map<int> &map);

template<typename T> void write_Healpix_map_to_dmc
  (const string &name, const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ, const Healpix_Map<T> &mapU, bool dp)
  {
  planck_assert (mapT.conformable(mapQ) && mapT.conformable(mapU),
    "write_Healpix_map_to_dmc: maps are not conformable");

  safe_ptr<iohandle> out (HandleManager.createObject(name,mapType(true,dp)));
  write_Healpix_keys (*out, mapT);
  int icol=out->columnNumber("I_Stokes"),
      qcol=out->columnNumber("Q_Stokes"),
      ucol=out->columnNumber("U_Stokes");
  chunkMaker cm (mapT.Npix(),out->efficientChunkSize(icol));
  uint64 offset, ppix;
  while (cm.getNext (offset,ppix))
    {
    out->appendColumnRaw(icol,&mapT[offset],ppix);
    out->appendColumnRaw(qcol,&mapQ[offset],ppix);
    out->appendColumnRaw(ucol,&mapU[offset],ppix);
    }
  }

template void write_Healpix_map_to_dmc
  (const string &name, const Healpix_Map<float> &mapT,
   const Healpix_Map<float> &mapQ, const Healpix_Map<float> &mapU, bool dp);
template void write_Healpix_map_to_dmc
  (const string &name, const Healpix_Map<double> &mapT,
   const Healpix_Map<double> &mapQ, const Healpix_Map<double> &mapU, bool dp);
