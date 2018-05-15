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
 *  Copyright (C) 2003, 2005, 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "healpix_data_io.h"
#include "arr.h"
#include "iohandle.h"
#include "paramfile.h"

using namespace std;

void read_weight_ring (const string &weightfile, arr<double> &weight)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject(weightfile,
    "constant.LS_ringweight"));
  inp->readColumn("Weight",weight);
  }

void get_ring_weights (paramfile &params, int nside, arr<double> &weight)
  {
  string weightfile = params.find<string>("ringweights","");
  weight.alloc (2*nside);
  if (weightfile!="")
    {
    read_weight_ring (weightfile, weight);
    for (tsize m=0; m<weight.size(); ++m) weight[m]+=1;
    }
  else
    weight.fill(1);
  }

void read_pixwin (const string &file, arr<double> &temp)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject(file,"constant.LS_pixwin"));
  if (temp.size()==0)
    inp->readEntireColumn("Temperature",temp);
  else
    inp->readColumn("Temperature",temp);
  }
void read_pixwin (const string &file, arr<double> &temp, arr<double> &pol)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject(file,"constant.LS_pixwin"));
  if (temp.size()==0)
    inp->readEntireColumn("Temperature",temp);
  else
    inp->readColumn("Temperature",temp);
  if (pol.size()==0)
    inp->readEntireColumn("Polarisation",pol);
  else
    inp->readColumn("Polarisation",pol);
  }

void get_pixwin (paramfile &params, int lmax, arr<double> &pixwin)
  {
  string windowfile = params.find<string>("windowfile","");
  pixwin.alloc(lmax+1);
  pixwin.fill(1);
  if (windowfile!="")
    read_pixwin (windowfile,pixwin);
  }
void get_pixwin (paramfile &params, int lmax, arr<double> &pixwin,
  arr<double> &pixwin_pol)
  {
  string windowfile = params.find<string>("windowfile","");
  pixwin.alloc(lmax+1);
  pixwin.fill(1);
  pixwin_pol.alloc(lmax+1);
  pixwin_pol.fill(1);
  if (windowfile!="")
    read_pixwin (windowfile,pixwin,pixwin_pol);
  }
