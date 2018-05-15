/*
 *  This file is part of libcxxmod.
 *
 *  libcxxmod is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxmod is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxmod; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxmod is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2012-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <fstream>
#include "ptcor.h"
#include "math_utils.h"
#include "lsconstants.h"

using namespace std;

ptcor::ptcor (paramfile &params, double theta_b)
  {
  matLos.Make_Axis_Rotation_Transform(vec3(0.,1.,0.),halfpi-theta_b);
  matLosInv.Make_Axis_Rotation_Transform(vec3(0.,1.,0.),theta_b-halfpi);
  ifstream inp(params.find<string>("ptcor_file").c_str());
  planck_assert(inp,"could not open ptcor CSV file");
  string dummy;
  getline(inp,dummy); // skip first line
  while (inp)
    {
    string token;
    getline(inp, token,',');
    if (trim(token)=="") break; // can happen on the last(empty) line
    time.push_back(stringToData<double>(token));
    getline(inp, token,',');
    incor.push_back(stringToData<double>(token));
    getline(inp, token);
    xcor.push_back(stringToData<double>(token));
    }
  }

rotmatrix ptcor::get_matrix (double obt) const
  {
  tsize idx;
  double frac;
  interpol_helper (time.begin(),time.end(),obt,idx,frac);
  // sanity check: warn if the time is  far outside the values in the CSV file
  // (this most likely means that the time conventions differ)
  if ((frac<-10.)||(frac>11.))
    cout << "Warning: weird interpolation values in ptcor: " << frac << endl;
  double xscan=xcor[idx]+frac*(xcor[idx+1]-xcor[idx]);
  double inscan=incor[idx]+frac*(incor[idx+1]-incor[idx]);
  rotmatrix matX,matIn;
  matX.Make_Axis_Rotation_Transform(vec3(0,1,0),xscan);
  matIn.Make_Axis_Rotation_Transform(vec3(1,0,0),inscan);
  return matLos*((matX*matIn)*matLosInv);
  }
