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
 *  Copyright (C) 2006-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "iohandle_current.h"
#include "iohandle_fits.h"
#include "alm_dmcio.h"
#include "iohandle_helpers.h"
#include "announce.h"

using namespace std;

namespace {

void checkColumns (const iohandle &inp, const iohandle &out)
  {
  for (int i=0; i<inp.numColumns(); ++i)
    if (!out.columnPresent(inp.columnName(i)))
      cout << "Warning: input column '" + inp.columnName(i)
        + "' not present in output object!" << endl;
  for (int i=0; i<out.numColumns(); ++i)
    if (!inp.columnPresent(out.columnName(i)))
      cout << "Warning: output column '" + out.columnName(i)
        + "' not present in input object!" << endl;
  }

void copyKeyTypeless (const iohandle &inp, iohandle &out, const string &key)
  {
  cout << "copying key '" << key <<"'" << endl;
  out.setKeyTypeless(key,inp.getKeyTypeless(key));
  }

void copyColumnTypeless (const iohandle &inp, iohandle &out, const string &name)
  {
  if (!out.columnPresent(name))
    {
    cout << "skipping column '" << name << "'" << endl;
    return;
    }

  PDT type = inp.columnType(name);
  cout << "copying column '" << name << "' ("
       << type2string(type) << ")" << endl;
  switch (type)
    {
    case PLANCK_BOOL:
      copyEntireColumn<bool> (inp, out, name);
      break;
    case PLANCK_FLOAT32:
      copyEntireColumn<float32> (inp, out, name);
      break;
    case PLANCK_FLOAT64:
      copyEntireColumn<float64> (inp, out, name);
      break;
    case PLANCK_INT8:
      copyEntireColumn<int8> (inp, out, name);
      break;
    case PLANCK_UINT8:
      copyEntireColumn<uint8> (inp, out, name);
      break;
    case PLANCK_INT16:
      copyEntireColumn<int16> (inp, out, name);
      break;
    case PLANCK_UINT16:
      copyEntireColumn<uint16> (inp, out, name);
      break;
    case PLANCK_INT32:
      copyEntireColumn<int32> (inp, out, name);
      break;
    case PLANCK_UINT32:
      copyEntireColumn<uint32> (inp, out, name);
      break;
    case PLANCK_INT64:
      copyEntireColumn<int64> (inp, out, name);
      break;
    case PLANCK_UINT64:
      copyEntireColumn<uint64> (inp, out, name);
      break;
    case PLANCK_STRING:
      copyEntireColumn<string> (inp, out, name);
      break;
    default:
      planck_fail (string("unsupported column type '")+type2string(type)+"'");
      break;
    }
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("dmc_transfer", argc, argv, 6,
    "<name1> <DDL type/'auto'> <name2> <DDL type/'auto'> <DMC init string>");

  iohandle_current::Manager mng (argv[5]);

  safe_ptr<iohandle> inp, out;
  string type_in;
  if (string(argv[2])=="auto")
    {
    inp = HandleManager.openObject(argv[1]);
    type_in = inp->objectType();
    }
  else
    {
    type_in = argv[2];
    inp = HandleManager.openObject(argv[1],type_in);
    }

  string type_out = (string(argv[4])=="auto") ? type_in : argv[4];
  out = HandleManager.createObject(argv[3],type_out);

  cout << "transferring an object of type '" << type_in
       << "' to an object of type '" << type_out << "'." << endl;

  vector<string> keys;
  inp->getAllKeys (keys);
  for (tsize i=0; i<keys.size(); ++i)
    copyKeyTypeless (*inp, *out, keys[i]);
  checkColumns(*inp,*out);
  for (int i=0; i<inp->numColumns(); ++i)
    copyColumnTypeless (*inp, *out, inp->columnName(i));
PLANCK_DIAGNOSIS_END
  }
