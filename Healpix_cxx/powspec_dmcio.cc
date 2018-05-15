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
 *  Copyright (C) 2003-2010 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "powspec_dmcio.h"
#include "powspec.h"
#include "iohandle.h"

using namespace std;

namespace {

string pstype (int nspecs)
  {
  if (nspecs==1) return "powerspectrum.LS_powspec";
  if (nspecs==4) return "powerspectrum.LS_powspec_pol";
  if (nspecs==6) return "powerspectrum.LS_powspec_pol_full";
  planck_fail ("incorrect number of spectra");
  }

} // unnamed namespace

void read_powspec_from_dmc (iohandle &inp, PowSpec &powspec,
  int nspecs, int lmax)
  {
  planck_assert ((nspecs==1)||(nspecs==4)||(nspecs==6),
    "wrong number of spectra");
  planck_assert (inp.columnLength("TT")>=tsize(lmax),
    "trying to read more C_l than are available in the input object");
  arr<double> tt(lmax+1), gg(lmax+1), cc(lmax+1), tg(lmax+1),
    tc(lmax+1), gc(lmax+1);

  inp.readColumn ("TT",tt);
  if (nspecs>=4)
    {
    inp.readColumn ("EE",gg);
    inp.readColumn ("BB",cc);
    inp.readColumn ("TE",tg);
    }
  if (nspecs==6)
    {
    inp.readColumn ("TB",tc);
    inp.readColumn ("EB",gc);
    }

  if (nspecs==1) powspec.Set(tt);
  if (nspecs==4) powspec.Set(tt,gg,cc,tg);
  if (nspecs==6) powspec.Set(tt,gg,cc,tg,tc,gc);
  }

void read_powspec_from_dmc (const string &name, PowSpec &powspec,
  int nspecs, int lmax)
  {
  planck_assert ((nspecs==1)||(nspecs==4), "wrong number of spectra");
  safe_ptr<iohandle> inp (HandleManager.openObject(name, pstype(nspecs)));
  read_powspec_from_dmc (*inp, powspec, nspecs, lmax);
  }

void write_powspec_to_dmc (iohandle &out, const PowSpec &powspec, int nspecs)
  {
  planck_assert ((nspecs==1)||(nspecs==4)||(nspecs==6),
    "wrong number of spectra");
  out.appendColumn("TT",powspec.tt());
  if (nspecs>1)
    {
    out.appendColumn("EE",powspec.gg());
    out.appendColumn("BB",powspec.cc());
    out.appendColumn("TE",powspec.tg());
    }
  if (nspecs>4)
    {
    out.appendColumn("TB",powspec.tc());
    out.appendColumn("EB",powspec.gc());
    }
  }

void write_powspec_to_dmc (const string &name, const PowSpec &powspec,
  int nspecs)
  {
  planck_assert ((nspecs==1)||(nspecs==4)||(nspecs==6),
    "wrong number of spectra");
  safe_ptr<iohandle> out (HandleManager.createObject(name, pstype(nspecs)));
  write_powspec_to_dmc (*out, powspec, nspecs);
  }
