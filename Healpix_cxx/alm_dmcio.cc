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

#include <string>
#include "alm_dmcio.h"
#include "alm.h"
#include "iohandle.h"
#include "xcomplex.h"
#include "share_utils.h"

using namespace std;

namespace {

string almType (bool pol, bool dp)
  {
  string res = "alm.LS_alm";
  if (pol) res += "_pol";
  if (dp) res += "_dp";
  return res;
  }

} // unnamed namespace

void get_almsize_dmc (iohandle &inp, char suffix, int &lmax, int &mmax)
  {
  string lms = string("lmax")+suffix, mms = string("mmax")+suffix;

  if (inp.keyPresent(lms) && inp.keyPresent(mms))
    {
    lmax = inp.getKey<int> (lms);
    mmax = inp.getKey<int> (mms);
    return;
    }

  cerr << "Warning: keywords " << lms << " and " << mms << " not found.\n"
       << "Trying to determine the values from the index column." << endl;

  int icol = inp.columnNumber(string("Index")+suffix);
  int n_alms = inp.columnLength(icol);
  arr<int> index;
  lmax=mmax=-1;

  chunkMaker cm (n_alms,inp.efficientChunkSize(icol));
  uint64 offset, ppix;
  while (cm.getNext (offset,ppix))
    {
    index.alloc(ppix);
    inp.readColumn(icol,index,offset);

    for (tsize i=0; i<ppix; ++i)
      {
      int l = isqrt(index[i]-1);
      int m = index[i] - l*l - l - 1;
      if (l>lmax) lmax=l;
      if (m>mmax) mmax=m;
      }
    }
  cerr << lms << " = " << lmax << "; " << mms << " = " << mmax << endl;
  }

void get_almsize_dmc(const string &name, int &lmax, int &mmax, bool dp)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject(name,almType(false,dp)));
  get_almsize_dmc (*inp,'T',lmax,mmax);
  }

void get_almsize_dmc_pol(const string &name, int &lmax, int &mmax, bool dp)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject(name,almType(true,dp)));
  int lmax2, mmax2;
  get_almsize_dmc (*inp,'T',lmax,mmax);
  get_almsize_dmc (*inp,'G',lmax2,mmax2);
  lmax=max(lmax,lmax2); mmax=max(mmax,mmax2);
  get_almsize_dmc (*inp,'C',lmax2,mmax2);
  lmax=max(lmax,lmax2); mmax=max(mmax,mmax2);
  }

template<typename T> void read_Alm_from_dmc
  (iohandle &inp, Alm<xcomplex<T> >&alms, int lmax, int mmax,
   char suffix)
  {
  int col_idx = inp.columnNumber(string("Index")+suffix);
  int col_re = inp.columnNumber(string("Real")+suffix);
  int col_im = inp.columnNumber(string("Imag")+suffix);
  int n_alms = inp.columnLength(col_idx);

  arr<int> index;
  arr<T> re, im;

  alms.Set(lmax, mmax);
  alms.SetToZero();
  int max_index = lmax*lmax + lmax + mmax + 1;
  chunkMaker cm (n_alms,inp.efficientChunkSize(col_idx));
  uint64 offset, ppix;
  while (cm.getNext (offset,ppix))
    {
    index.alloc(ppix);
    re.alloc(ppix); im.alloc(ppix);
    inp.readColumn(col_idx,index,offset);
    inp.readColumn(col_re,re,offset);
    inp.readColumn(col_im,im,offset);

    for (tsize i=0; i<ppix; ++i)
      {
      if (index[i]>max_index) continue;

      int l = int(sqrt(double(index[i]-1)));
      int m = index[i] - l*l - l - 1;
      planck_assert(m>=0,"negative m encountered");
      planck_assert(l>=m, "wrong l,m combination");
      if ((l<=lmax) && (m<=mmax))
        alms(l,m).Set (re[i], im[i]);
      }
    }
  }

template void read_Alm_from_dmc (iohandle &inp, Alm<xcomplex<float> >&alms,
  int lmax, int mmax, char suffix);
template void read_Alm_from_dmc (iohandle &inp, Alm<xcomplex<double> >&alms,
  int lmax, int mmax, char suffix);

template<typename T> void read_Alm_from_dmc
  (const string &name, Alm<xcomplex<T> >&alms, int lmax, int mmax, bool dp)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject(name,almType(false,dp)));
  read_Alm_from_dmc (*inp,alms,lmax,mmax,'T');
  }

template void read_Alm_from_dmc (const string &name,
   Alm<xcomplex<float> >&alms, int lmax, int mmax, bool dp);
template void read_Alm_from_dmc (const string &name,
   Alm<xcomplex<double> >&alms, int lmax, int mmax, bool dp);

template<typename T> void read_Alm_from_dmc
  (const string &name, Alm<xcomplex<T> >&almT, Alm<xcomplex<T> >&almG,
   Alm<xcomplex<T> >&almC, int lmax, int mmax, bool dp)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject(name,almType(true,dp)));
  read_Alm_from_dmc (*inp,almT,lmax,mmax,'T');
  read_Alm_from_dmc (*inp,almG,lmax,mmax,'G');
  read_Alm_from_dmc (*inp,almC,lmax,mmax,'C');
  }

template void read_Alm_from_dmc (const string &name,
   Alm<xcomplex<float> >&almT, Alm<xcomplex<float> >&almG,
   Alm<xcomplex<float> >&almC, int lmax, int mmax, bool dp);
template void read_Alm_from_dmc (const string &name,
   Alm<xcomplex<double> >&almT, Alm<xcomplex<double> >&almG,
   Alm<xcomplex<double> >&almC, int lmax, int mmax, bool dp);

#if 0

template<typename T> void write_compressed_Alm_to_dmc
  (iohandle &out, const Alm<xcomplex<T> > &alms, int lmax, int mmax,
   char suffix)
  {
  int col_idx = out.columnNumber(string("Index")+suffix);
  int col_re = out.columnNumber(string("Real")+suffix);
  int col_im = out.columnNumber(string("Imag")+suffix);
  arr<int> index;
  arr<double> re, im;

  int n_alms = 0;
  for (int m=0; m<=mmax; ++m)
    for (int l=m; l<=lmax; ++l)
      if (alms(l,m).norm()>0) ++n_alms;

  int l=0, m=0;
  int real_lmax = 0, real_mmax=0;
  chunkMaker cm (n_alms,out.efficientChunkSize(col_idx));
  uint64 offset, ppix;
  while (cm.getNext (offset,ppix))
    {
    index.alloc(ppix);
    re.alloc(ppix); im.alloc(ppix);
    for (int i=0; i<ppix; ++i)
      {
      while (alms(l,m).norm()==0)
        {
        ++m;
        if ((m>l) || (m>mmax)) { ++l; m=0; }
        }
      index[i] = l*l + l + m + 1;
      re[i] = alms(l,m).re;
      im[i] = alms(l,m).im;
      if (l>real_lmax) real_lmax=l;
      if (m>real_mmax) real_mmax=m;
      ++m;
      if ((m>l) || (m>mmax)) { ++l; m=0; }
      }
    out.appendColumn(col_idx,index);
    out.appendColumn(col_re,re);
    out.appendColumn(col_im,im);
    }
  out.setKey(string("lmax")+suffix,real_lmax);
  out.setKey(string("mmax")+suffix,real_mmax);
  }

#endif

template<typename T> void write_Alm_to_dmc
  (iohandle &out, const Alm<xcomplex<T> > &alms, int lmax, int mmax,
   char suffix)
  {
  int col_idx = out.columnNumber(string("Index")+suffix);
  int col_re = out.columnNumber(string("Real")+suffix);
  int col_im = out.columnNumber(string("Imag")+suffix);
  arr<int> index;
  arr<double> re, im;

  int lm=alms.Lmax(), mm=alms.Mmax();
  int n_alms = ((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax);

  int l=0, m=0;
  chunkMaker cm (n_alms,out.efficientChunkSize(col_idx));
  uint64 offset, ppix;
  while (cm.getNext (offset,ppix))
    {
    index.alloc(ppix);
    re.alloc(ppix); im.alloc(ppix);
    for (tsize i=0; i<ppix; ++i)
      {
      index[i] = l*l + l + m + 1;
      if ((l<=lm) && (m<=mm))
        { re[i] = alms(l,m).re; im[i] = alms(l,m).im; }
      else
        { re[i] = 0; im[i] = 0; }
      ++m;
      if ((m>l) || (m>mmax)) { ++l; m=0; }
      }
    out.appendColumn(col_idx,index);
    out.appendColumn(col_re,re);
    out.appendColumn(col_im,im);
    }
  out.setKey(string("lmax")+suffix,lmax);
  out.setKey(string("mmax")+suffix,mmax);
  }

template void write_Alm_to_dmc (iohandle &out,
  const Alm<xcomplex<float> > &alms, int lmax, int mmax, char suffix);
template void write_Alm_to_dmc (iohandle &out,
  const Alm<xcomplex<double> > &alms, int lmax, int mmax, char suffix);

template<typename T> void write_Alm_to_dmc (const string &name,
  const Alm<xcomplex<T> > &alms, int lmax, int mmax, bool dp)
  {
  safe_ptr<iohandle> out (HandleManager.createObject(name,almType(false,dp)));
  write_Alm_to_dmc (*out,alms,lmax,mmax,'T');
  }

template void write_Alm_to_dmc (const string &name,
  const Alm<xcomplex<double> > &alms, int lmax, int mmax, bool dp);
template void write_Alm_to_dmc (const string &name,
  const Alm<xcomplex<float> > &alms, int lmax, int mmax, bool dp);

template<typename T> void write_Alm_to_dmc
  (const string &name, const Alm<xcomplex<T> > &almT,
   const Alm<xcomplex<T> > &almG, const Alm<xcomplex<T> > &almC,
   int lmax, int mmax, bool dp)
  {
  safe_ptr<iohandle> out (HandleManager.createObject(name,almType(true,dp)));
  write_Alm_to_dmc (*out,almT,lmax,mmax,'T');
  write_Alm_to_dmc (*out,almG,lmax,mmax,'G');
  write_Alm_to_dmc (*out,almC,lmax,mmax,'C');
  }

template void write_Alm_to_dmc
  (const string &name, const Alm<xcomplex<float> > &almT,
   const Alm<xcomplex<float> > &almG, const Alm<xcomplex<float> > &almC,
   int lmax, int mmax, bool dp);
template void write_Alm_to_dmc
  (const string &name, const Alm<xcomplex<double> > &almT,
   const Alm<xcomplex<double> > &almG, const Alm<xcomplex<double> > &almC,
   int lmax, int mmax, bool dp);
