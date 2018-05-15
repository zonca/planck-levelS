/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file io_utils.h
 *  Convenience functions for I/O with iohandles.
 *
 *  Copyright (C) 2009-2014 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_IO_UTILS_H
#define PLANCK_IO_UTILS_H

#include "iohandle.h"
#include "pointing.h"
#include "quaternion.h"
#include "arr.h"
#include "share_utils.h"

template<typename A, typename B> inline void writeMulticol(iohandle &out,
  int cn1, const A &a1, int cn2, const B &a2)
  {
  planck_assert(multiequal(a1.size(),a2.size()),"array size mismatch");
  chunkMaker cm (a1.size(),out.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    out.appendColumnRaw(cn1,&a1[offset2],nval);
    out.appendColumnRaw(cn2,&a2[offset2],nval);
    }
  }
template<typename A, typename B> inline void writeMulticol(iohandle &out,
  const std::string &c1, const A &a1,
  const std::string &c2, const B &a2)
  {
  writeMulticol (out, out.columnNumber(c1), a1, out.columnNumber(c2), a2);
  }

template<typename A, typename B, typename C> inline void writeMulticol
  (iohandle &out, int cn1, const A &a1, int cn2, const B &a2,
    int cn3, const C &a3)
  {
  planck_assert(multiequal(a1.size(),a2.size(),a3.size()),
    "array size mismatch");
  chunkMaker cm (a1.size(),out.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    out.appendColumnRaw(cn1,&a1[offset2],nval);
    out.appendColumnRaw(cn2,&a2[offset2],nval);
    out.appendColumnRaw(cn3,&a3[offset2],nval);
    }
  }
template<typename A, typename B, typename C> inline void writeMulticol
  (iohandle &out, const std::string &c1, const A &a1,
   const std::string &c2, const B &a2, const std::string &c3, const C &a3)
  {
  writeMulticol (out, out.columnNumber(c1), a1, out.columnNumber(c2), a2,
    out.columnNumber(c3), a3);
  }

template<typename A, typename B, typename C, typename D>
  inline void writeMulticol (iohandle &out,
  int cn1, const A &a1, int cn2, const B &a2,
  int cn3, const C &a3, int cn4, const D &a4)
  {
  planck_assert(multiequal(a1.size(),a2.size(),a3.size(),a4.size()),
    "array size mismatch");
  chunkMaker cm (a1.size(),out.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    out.appendColumnRaw(cn1,&a1[offset2],nval);
    out.appendColumnRaw(cn2,&a2[offset2],nval);
    out.appendColumnRaw(cn3,&a3[offset2],nval);
    out.appendColumnRaw(cn4,&a4[offset2],nval);
    }
  }
template<typename A, typename B, typename C, typename D>
  inline void writeMulticol (iohandle &out,
  const std::string &c1, const A &a1,
  const std::string &c2, const B &a2,
  const std::string &c3, const C &a3,
  const std::string &c4, const D &a4)
  {
  writeMulticol (out, out.columnNumber(c1), a1, out.columnNumber(c2), a2,
    out.columnNumber(c3), a3, out.columnNumber(c4), a4);
  }

template<typename A, typename B, typename C, typename D, typename E>
  inline void writeMulticol (iohandle &out,
  int cn1, const A &a1, int cn2, const B &a2,
  int cn3, const C &a3, int cn4, const D &a4,
  int cn5, const E &a5)
  {
  planck_assert(multiequal(a1.size(),a2.size(),a3.size(),a4.size(),a5.size()),
    "array size mismatch");
  chunkMaker cm (a1.size(),out.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    out.appendColumnRaw(cn1,&a1[offset2],nval);
    out.appendColumnRaw(cn2,&a2[offset2],nval);
    out.appendColumnRaw(cn3,&a3[offset2],nval);
    out.appendColumnRaw(cn4,&a4[offset2],nval);
    out.appendColumnRaw(cn5,&a5[offset2],nval);
    }
  }
template<typename A, typename B, typename C, typename D, typename E>
  inline void writeMulticol (iohandle &out,
  const std::string &c1, const A &a1,
  const std::string &c2, const B &a2,
  const std::string &c3, const C &a3,
  const std::string &c4, const D &a4,
  const std::string &c5, const E &a5)
  {
  writeMulticol (out, out.columnNumber(c1), a1, out.columnNumber(c2), a2,
    out.columnNumber(c3), a3, out.columnNumber(c4), a4,
    out.columnNumber(c5), a5);
  }
template<typename A, typename B, typename C, typename D, typename E, typename F>
  inline void writeMulticol (iohandle &out,
  int cn1, const A &a1, int cn2, const B &a2,
  int cn3, const C &a3, int cn4, const D &a4,
  int cn5, const E &a5, int cn6, const F &a6)
  {
  planck_assert(multiequal(a1.size(),a2.size(),a3.size(),a4.size(),a5.size(),
    a6.size()),"array size mismatch");
  chunkMaker cm (a1.size(),out.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    out.appendColumnRaw(cn1,&a1[offset2],nval);
    out.appendColumnRaw(cn2,&a2[offset2],nval);
    out.appendColumnRaw(cn3,&a3[offset2],nval);
    out.appendColumnRaw(cn4,&a4[offset2],nval);
    out.appendColumnRaw(cn5,&a5[offset2],nval);
    out.appendColumnRaw(cn6,&a6[offset2],nval);
    }
  }
template<typename A, typename B, typename C, typename D, typename E, typename F>
  inline void writeMulticol (iohandle &out,
  const std::string &c1, const A &a1,
  const std::string &c2, const B &a2,
  const std::string &c3, const C &a3,
  const std::string &c4, const D &a4,
  const std::string &c5, const E &a5,
  const std::string &c6, const F &a6)
  {
  writeMulticol (out, out.columnNumber(c1), a1, out.columnNumber(c2), a2,
    out.columnNumber(c3), a3, out.columnNumber(c4), a4,
    out.columnNumber(c5), a5, out.columnNumber(c6), a6);
  }

template<typename A, typename B> inline void readMulticol(iohandle &inp,
  int cn1, A &a1, int cn2, B &a2, tsize len, uint64 offset)
  {
  a1.resize(len); a2.resize(len);
  chunkMaker cm (len,inp.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    inp.readColumnRaw(cn1,&a1[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn2,&a2[offset2],nval,offset+offset2);
    }
  }
template<typename A, typename B> inline void readMulticol(iohandle &inp,
  const std::string &c1, A &a1,
  const std::string &c2, B &a2, tsize len, uint64 offset)
  {
  readMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    len, offset);
  }
template<typename A, typename B, typename C> inline void readMulticol
  (iohandle &inp, int cn1, A &a1, int cn2, B &a2, int cn3, C &a3,
  tsize len, uint64 offset)
  {
  a1.resize(len); a2.resize(len); a3.resize(len);
  chunkMaker cm (len,inp.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    inp.readColumnRaw(cn1,&a1[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn2,&a2[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn3,&a3[offset2],nval,offset+offset2);
    }
  }
template<typename A, typename B, typename C> inline void readMulticol
  (iohandle &inp, const std::string &c1, A &a1,
  const std::string &c2, B &a2, const std::string &c3, C &a3,
  tsize len, uint64 offset)
  {
  readMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    inp.columnNumber(c3), a3, len, offset);
  }
template<typename A, typename B, typename C, typename D>
  inline void readMulticol(iohandle &inp, int cn1, A &a1,
  int cn2, B &a2, int cn3, C &a3,
  int cn4, D &a4, tsize len, uint64 offset)
  {
  a1.resize(len); a2.resize(len); a3.resize(len); a4.resize(len);
  chunkMaker cm (len,inp.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    inp.readColumnRaw(cn1,&a1[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn2,&a2[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn3,&a3[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn4,&a4[offset2],nval,offset+offset2);
    }
  }
template<typename A, typename B, typename C, typename D>
  inline void readMulticol(iohandle &inp, const std::string &c1, A &a1,
  const std::string &c2, B &a2, const std::string &c3, C &a3,
  const std::string &c4, D &a4, tsize len, uint64 offset)
  {
  readMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    inp.columnNumber(c3), a3, inp.columnNumber(c4), a4, len, offset);
  }
template<typename A, typename B, typename C, typename D, typename E>
  inline void readMulticol(iohandle &inp, int cn1, A &a1,
  int cn2, B &a2, int cn3, C &a3, int cn4, D &a4,
  int cn5, E &a5, tsize len, uint64 offset)
  {
  a1.resize(len);a2.resize(len);a3.resize(len);a4.resize(len);a5.resize(len);
  chunkMaker cm (len,inp.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    inp.readColumnRaw(cn1,&a1[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn2,&a2[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn3,&a3[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn4,&a4[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn5,&a5[offset2],nval,offset+offset2);
    }
  }
template<typename A, typename B, typename C, typename D, typename E>
  inline void readMulticol(iohandle &inp, const std::string &c1, A &a1,
  const std::string &c2, B &a2, const std::string &c3, C &a3,
  const std::string &c4, D &a4, const std::string &c5, E &a5,
  tsize len, uint64 offset)
  {
  readMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    inp.columnNumber(c3), a3, inp.columnNumber(c4), a4,
    inp.columnNumber(c5), a5, len, offset);
  }
template<typename A, typename B, typename C, typename D, typename E, typename F>
  inline void readMulticol(iohandle &inp, int cn1, A &a1,
  int cn2, B &a2, int cn3, C &a3, int cn4, D &a4,
  int cn5, E &a5, int cn6, F &a6, tsize len, uint64 offset)
  {
  a1.resize(len);a2.resize(len);a3.resize(len);a4.resize(len);a5.resize(len);
  a6.resize(len);
  chunkMaker cm (len,inp.efficientChunkSize(cn1));
  uint64 offset2, nval;
  while (cm.getNext (offset2,nval))
    {
    inp.readColumnRaw(cn1,&a1[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn2,&a2[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn3,&a3[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn4,&a4[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn5,&a5[offset2],nval,offset+offset2);
    inp.readColumnRaw(cn6,&a6[offset2],nval,offset+offset2);
    }
  }
template<typename A, typename B, typename C, typename D, typename E, typename F>
  inline void readMulticol(iohandle &inp, const std::string &c1, A &a1,
  const std::string &c2, B &a2, const std::string &c3, C &a3,
  const std::string &c4, D &a4, const std::string &c5, E &a5,
  const std::string &c6, F &a6, tsize len, uint64 offset)
  {
  readMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    inp.columnNumber(c3), a3, inp.columnNumber(c4), a4,
    inp.columnNumber(c5), a5, inp.columnNumber(c6), a6, len, offset);
  }

template<typename A, typename B> inline void readEntireMulticol(iohandle &inp,
  int cn1, A &a1, int cn2, B &a2)
  {
  planck_assert(multiequal(inp.columnLength(cn1),inp.columnLength(cn2)),
    "column length mismatch");
  readMulticol(inp,cn1,a1,cn2,a2,inp.columnLength(cn1),0);
  }
template<typename A, typename B> inline void readEntireMulticol(iohandle &inp,
  const std::string &c1, A &a1, const std::string &c2, B &a2)
  {
  readEntireMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2);
  }
template<typename A, typename B, typename C> inline void readEntireMulticol
  (iohandle &inp, int cn1, A &a1, int cn2, B &a2, int cn3, C &a3)
  {
  planck_assert(multiequal(inp.columnLength(cn1),inp.columnLength(cn2),
    inp.columnLength(cn3)), "column length mismatch");
  readMulticol(inp,cn1,a1,cn2,a2,cn3,a3,inp.columnLength(cn1),0);
  }
template<typename A, typename B, typename C> inline void readEntireMulticol
  (iohandle &inp, const std::string &c1, A &a1,
   const std::string &c2, B &a2, const std::string &c3, C &a3)
  {
  readEntireMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    inp.columnNumber(c3), a3);
  }
template<typename A, typename B, typename C, typename D>
  inline void readEntireMulticol (iohandle &inp,
    int cn1, A &a1, int cn2, B &a2, int cn3, C &a3,
    int cn4, D &a4)
  {
  planck_assert(multiequal(inp.columnLength(cn1),inp.columnLength(cn2),
    inp.columnLength(cn3),inp.columnLength(cn4)), "column length mismatch");
  readMulticol(inp,cn1,a1,cn2,a2,cn3,a3,cn4,a4,inp.columnLength(cn1),0);
  }
template<typename A, typename B, typename C, typename D>
  inline void readEntireMulticol (iohandle &inp,
    const std::string &c1, A &a1, const std::string &c2, B &a2,
    const std::string &c3, C &a3, const std::string &c4, D &a4)
  {
  readEntireMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    inp.columnNumber(c3), a3, inp.columnNumber(c4), a4);
  }
template<typename A, typename B, typename C, typename D, typename E>
  inline void readEntireMulticol (iohandle &inp,
    int cn1, A &a1, int cn2, B &a2, int cn3, C &a3,
    int cn4, D &a4, int cn5, E &a5)
  {
  planck_assert(multiequal(inp.columnLength(cn1),inp.columnLength(cn2),
    inp.columnLength(cn3),inp.columnLength(cn4),inp.columnLength(cn5)),
    "column length mismatch");
  readMulticol(inp,cn1,a1,cn2,a2,cn3,a3,cn4,a4,cn5,a5,inp.columnLength(cn1),0);
  }
template<typename A, typename B, typename C, typename D, typename E>
  inline void readEntireMulticol (iohandle &inp,
    const std::string &c1, A &a1, const std::string &c2, B &a2,
    const std::string &c3, C &a3, const std::string &c4, D &a4,
    const std::string &c5, E &a5)
  {
  readEntireMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    inp.columnNumber(c3), a3, inp.columnNumber(c4), a4,
    inp.columnNumber(c5), a5);
  }
template<typename A, typename B, typename C, typename D, typename E, typename F>
  inline void readEntireMulticol (iohandle &inp,
    int cn1, A &a1, int cn2, B &a2, int cn3, C &a3,
    int cn4, D &a4, int cn5, E &a5, int cn6, F &a6)
  {
  planck_assert(multiequal(inp.columnLength(cn1),inp.columnLength(cn2),
    inp.columnLength(cn3),inp.columnLength(cn4),inp.columnLength(cn5),
    inp.columnLength(cn6)), "column length mismatch");
  readMulticol(inp,cn1,a1,cn2,a2,cn3,a3,cn4,a4,cn5,a5,cn6,a6,
    inp.columnLength(cn1),0);
  }
template<typename A, typename B, typename C, typename D, typename E, typename F>
  inline void readEntireMulticol (iohandle &inp,
    const std::string &c1, A &a1, const std::string &c2, B &a2,
    const std::string &c3, C &a3, const std::string &c4, D &a4,
    const std::string &c5, E &a5, const std::string &c6, F &a6)
  {
  readEntireMulticol (inp, inp.columnNumber(c1), a1, inp.columnNumber(c2), a2,
    inp.columnNumber(c3), a3, inp.columnNumber(c4), a4,
    inp.columnNumber(c5), a5, inp.columnNumber(c6), a6);
  }

/*! Reads colatitudes from the column with the number \a ctheta and longitudes
    from the column with the number \a cphi from \a inp, starting at the column
    entry with the number \a offset, and returns them in \a pt. The size of
    \a pt determines the number of values read. */
void readPointing (iohandle &inp, int ctheta, int cphi,
  arr<pointing> &pt, uint64 offset);

/*! Reads colatitudes from the column with the number \a ctheta, longitudes
    from the column with the number \a cphi and orientation angles from the
    column with the number \a cpsi from \a inp, starting at the column
    entry with the number \a offset, and returns them in \a pt and \a psi.
    The size of \a pt determines the number of values read. \a pt and \a psi
    must have the same number of entries. */
void readPointing (iohandle &inp, int ctheta, int cphi, int cpsi,
  arr<pointing> &pt, arr<double> &psi, uint64 offset);

/*! Reads colatitudes from the column with the name \a ntheta and longitudes
    from the column with the name \a nphi from \a inp, starting at the column
    entry with the number \a offset, and returns them in \a pt. The size of
    \a pt determines the number of values read. */
inline void readPointing (iohandle &inp, const std::string &ntheta,
  const std::string &nphi, arr<pointing> &pt, uint64 offset)
  {
  readPointing (inp, inp.columnNumber(ntheta), inp.columnNumber(nphi),
                pt, offset);
  }

/*! Reads colatitudes from the column with the name \a ntheta, longitudes
    from the column with the name \a nphi and orientation angles from the
    column with the name \a npsi from \a inp, starting at the column
    entry with the number \a offset, and returns them in \a pt and \a psi.
    The size of \a pt determines the number of values read. \a pt and \a psi
    must have the same number of entries. */
inline void readPointing (iohandle &inp, const std::string &ntheta,
  const std::string &nphi, const std::string &npsi, arr<pointing> &pt,
  arr<double> &psi, uint64 offset)
  {
  readPointing (inp, inp.columnNumber(ntheta), inp.columnNumber(nphi),
                inp.columnNumber(npsi), pt, psi, offset);
  }

/*! Reads three sets of pointings from 6 columns of \a inp, starting at the
    column entry with the number \a offset. */
void readPointing (iohandle &inp, int ctx, int cpx, int cty, int cpy,
  int ctz, int cpz, arr<pointing> &px, arr<pointing> &py, arr<pointing> &pz,
  uint64 offset);


/*! Reads all colatitudes from the column with the number \a ctheta and
    all longitudes from the column with the number \a cphi from \a inp, and
    returns them in \a pt (which is resized accordingly). */
void readEntirePointing (iohandle &inp, int ctheta, int cphi,
  arr<pointing> &pt);

/*! Reads all colatitudes from the column with the name \a ntheta and
    all longitudes from the column with the name \a nphi from \a inp, and
    returns them in \a pt (which is resized accordingly). */
inline void readEntirePointing (iohandle &inp, const std::string &ntheta,
  const std::string &nphi, arr<pointing> &pt)
  {
  readEntirePointing (inp,inp.columnNumber(ntheta),inp.columnNumber(nphi),pt);
  }

/*! Appends the colatitudes and longitudes from \a pt to the columns with the
    numbers \a ctheta and \a cphi, respectively, in \a out. */
void appendPointing (iohandle &out, int ctheta, int cphi,
  const arr<pointing> &pt);

/*! Appends the colatitudes and longitudes from \a pt to the columns with the
    numbers \a ctheta and \a cphi, respectively, in \a out. Also appends the
    orientations in \a psi to the column with the number \a cpsi. */
void appendPointing (iohandle &out, int ctheta, int cphi, int cpsi,
  const arr<pointing> &pt, const arr<double> &psi);

/*! Reads the x,y,z, and w quaternion components from the columns with the
    numbers \a cx, \a cy, \a cz and \a cw in \a inp, starting at the column
    entry with the number \a offset, and returns them in \a quat. The size of
    \a quat determines the number of values read.   */
void readQuaternions (iohandle &inp, int cx, int cy, int cz, int cw,
  arr<quaternion> &quat, uint64 offset);

/*! Appends the x,y,z, and w components of \a quat to the columns with the
    numbers \a cx, \a cy, \a cz and \a cw in \a out. */
void appendQuaternions (iohandle &out, int cx, int cy, int cz, int cw,
  const arr<quaternion> &quat);

void readPointSources (const std::string &infile, const std::string &suffix,
  bool polarisation, arr<pointing> &pos, arr<double> &flux,
  arr<double> &polang, arr<double> &polfrac, arr<std::string> &name,
  int nshares, int myshare);

void readPointSourcesNew (const std::string &infile,
  bool polarisation, arr<pointing> &pos, arr<double> &flux,
  arr<double> &polang, arr<double> &polfrac, arr<std::string> &name,
  int nshares, int myshare);

#endif
