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

/*! \file healpix_map_dmcio.h
 *  Copyright (C) 2003, 2004, 2005, 2006 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_MAP_DMCIO_H
#define HEALPIX_MAP_DMCIO_H

#include <string>
class iohandle;

template<typename T> class Healpix_Map;

template<typename T> void read_Healpix_map_from_dmc
  (iohandle &inp, Healpix_Map<T> &map, const std::string &colname);
template<typename T> void read_Healpix_map_from_dmc (const std::string &name,
  Healpix_Map<T> &map, const std::string &colname, bool dp);
template<typename T> void read_Healpix_map_from_dmc
  (const std::string &name, Healpix_Map<T> &map, bool dp=false);
template<typename T> void read_Healpix_map_from_dmc
  (const std::string &name, Healpix_Map<T> &mapT, Healpix_Map<T> &mapQ,
   Healpix_Map<T> &mapU, bool dp=false);

template<typename T> void write_Healpix_map_to_dmc
  (const std::string &name, const Healpix_Map<T> &map, bool dp=false);
template<typename T> void write_Healpix_hitmap_to_dmc
  (const std::string &name, const Healpix_Map<T> &map);
template<typename T> void write_Healpix_map_to_dmc
  (const std::string &name, const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ, const Healpix_Map<T> &mapU, bool dp=false);

#endif
