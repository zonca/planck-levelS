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

/*! \file alm_dmcio.h
 *  DMC I/O for spherical harmonic coefficients
 *
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef ALM_DMCIO_H
#define ALM_DMCIO_H

#include <string>
#include "xcomplex.h"
class iohandle;

template<typename T> class Alm;

/*! \defgroup alm_dmcio_group DMC-based I/O of a_lm */
/*! \{ */

/*! Returns the maximum \a l and \a m multipole moments found in the object
    \a inp for the component \a suffix in \a lmax and \a mmax. */
void get_almsize_dmc (iohandle &inp, char suffix, int &lmax, int &mmax);
/*! Returns the maximum \a l and \a m multipole moments found in the object
    \a filename in \a lmax and \a mmax. */
void get_almsize_dmc(const std::string &name, int &lmax, int &mmax,
  bool dp=false);
/*! Returns the maximum \a l and \a m multipole moments (in all three
    polarisation components) found in the object \a filename in \a lmax
    and \a mmax. */
void get_almsize_dmc_pol(const std::string &name, int &lmax, int &mmax,
  bool dp=false);

template<typename T> void read_Alm_from_dmc
  (iohandle &inp, Alm<xcomplex<T> >&alms, int lmax, int mmax,
   char suffix);

template<typename T> void read_Alm_from_dmc (const std::string &name,
  Alm<xcomplex<T> > &alms, int lmax, int mmax, bool dp=false);

template<typename T> void read_Alm_from_dmc
  (const std::string &name, Alm<xcomplex<T> > &almT, Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC, int lmax, int mmax, bool dp=false);

template<typename T> void write_Alm_to_dmc
  (iohandle &out, const Alm<xcomplex<T> > &alms, int lmax, int mmax,
   char suffix);

template<typename T> void write_Alm_to_dmc (const std::string &name,
  const Alm<xcomplex<T> > &alms, int lmax, int mmax, bool dp=false);

template<typename T> void write_Alm_to_dmc
  (const std::string &name, const Alm<xcomplex<T> > &almT,
   const Alm<xcomplex<T> > &almG, const Alm<xcomplex<T> > &almC,
   int lmax, int mmax, bool dp=false);

/*! \} */

#endif
