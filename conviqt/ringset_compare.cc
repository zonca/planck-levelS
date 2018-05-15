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
 *  Copyright (C) 2009-2013 Max-Planck-Society
 *  Authors: Martin Reinecke
 */

#include "paramfile.h"
#include "iohandle_current.h"
#include "announce.h"
#include "ls_image.h"

using namespace std;

int main (int argc, const char **argv)
  {
  module_startup ("ringset_compare", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  safe_ptr<iohandle> inp1 (HandleManager.openObject
    (params.find<string>("ringset1"),"ringset.LS_ringset_dp"));
  safe_ptr<iohandle> inp2 (HandleManager.openObject
    (params.find<string>("ringset2"),"ringset.LS_ringset_dp"));
  arr<double> d1,d2;
  inp1->readEntireColumn("ringsetdata",d1);
  inp2->readEntireColumn("ringsetdata",d2);
  planck_assert(d1.size()==d2.size(),"array length mismatch");
  {
  double maxdiff=0,maxval=0;
  for (tsize i=0; i<d1.size(); ++i)
    {
    if (abs(d1[i])>maxval) maxval=abs(d1[i]);
    if (abs(d1[i]-d2[i])>maxdiff) maxdiff=abs(d1[i]-d2[i]);
    }
  cout << "maxval: " << maxval << endl;
  cout << "maxdiff: " << dataToString(maxdiff) << endl;
  }
  int nphi=inp1->getKey<int>("nphi");
  int ntheta=inp1->getKey<int>("ntheta");
  int ncomp=inp1->columnLength("ringsets_present")*2-1;
  for (int j=0; j<ncomp; ++j)
    {
    cout << "component " << j << endl;
    double maxdiff=0,maxval=0;
    tsize ofs=j*nphi*ntheta;
    for (tsize i=ofs; i<ofs+nphi*ntheta; ++i)
      {
      if (abs(d1[i])>maxval) maxval=abs(d1[i]);
      if (abs(d1[i]-d2[i])>maxdiff) maxdiff=abs(d1[i]-d2[i]);
      }
    cout << "maxval: " << maxval << endl;
    cout << "maxdiff: " << dataToString(maxdiff) << endl;
    LS_Image img(nphi,ntheta);
    for (tsize i=ofs; i<ofs+nphi*ntheta; ++i)
      {
      int x=(i-ofs)%nphi, y=(i-ofs)/nphi;
      double val=abs(d1[i]-d2[i])/maxdiff;
      img.put_pixel(x,y,Colour(val,val,val));
      }
    img.write_TGA (string("blub_")+intToString(j,2));
    }
  }
