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
 *  Copyright (C) 2010-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "iohandle_current.h"
#include "announce.h"
#include "string_utils.h"

using namespace std;

namespace {

template<typename T> void showColumn (iohandle &inp, int num)
  {
  uint64 len = inp.columnLength(num);
  arr<T> data;
  if (len<20)
    {
    inp.readEntireColumn(num, data);
    for (tsize i=0; i<data.size(); ++i)
      cout << intToString(i,2) << ": " << dataToString(data[i]) << endl;
    }
  else
    {
    data.alloc(10);
    inp.readColumn(num,data,0);
    for (tsize i=0; i<data.size(); ++i)
      cout << intToString(i,2) << ": " << dataToString(data[i]) << endl;
    cout << "[...]" << endl;
    inp.readColumn(num,data,len-10);
    for (tsize i=0; i<data.size(); ++i)
      cout << len-10+i << ": " << dataToString(data[i]) << endl;
    }
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("object_inspector", argc, argv, 2, "<object name>");

  iohandle_current::Manager mng ("");

  safe_ptr<iohandle> inp(HandleManager.openObject(argv[1]));
  cout << "Inspecting the object " << argv[1] << " ("
       << inp->objectType() << ")" << endl;
  for (int i=0; i<inp->numColumns(); ++i)
    {
    cout << "\nColumn #" << i << " ('" << inp->columnName(i) << "', '"
         << inp->columnUnit(i) << "', " << type2string(inp->columnType(i))
         << "):" << endl;
    switch (inp->columnType(i))
      {
      case PLANCK_FLOAT32:
        showColumn<float32>(*inp,i);
        break;
      case PLANCK_FLOAT64:
        showColumn<float64>(*inp,i);
        break;
      case PLANCK_INT8:
        showColumn<int8>(*inp,i);
        break;
      case PLANCK_INT16:
        showColumn<int16>(*inp,i);
        break;
      case PLANCK_INT32:
        showColumn<int32>(*inp,i);
        break;
      case PLANCK_INT64:
        showColumn<int64>(*inp,i);
        break;
      case PLANCK_UINT8:
        showColumn<uint8>(*inp,i);
        break;
      case PLANCK_UINT16:
        showColumn<uint16>(*inp,i);
        break;
      case PLANCK_UINT32:
        showColumn<uint32>(*inp,i);
        break;
      case PLANCK_UINT64:
        showColumn<uint64>(*inp,i);
        break;
      case PLANCK_BOOL:
        showColumn<bool>(*inp,i);
        break;
      case PLANCK_STRING:
        showColumn<string>(*inp,i);
        break;
      default:
        cout << "unsupported data type" << endl;
        break;
      }
    }
PLANCK_DIAGNOSIS_END
  }
