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
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_WRITERS_H
#define PLANCK_WRITERS_H

#include <string>
#include "iohandle.h"
#include "pointing.h"
#include "quaternion.h"

template <typename T> class arr;

class generic_writer
  {
  protected:
    safe_ptr<iohandle> out;

  public:
    template<typename T>void setKey(const std::string &key,
      const T &value)
      { out->setKey(key, value); }
  };

class detpt_writer: public generic_writer
  {
  private:
    int ctheta, cphi, cpsi;

  public:
    detpt_writer (bool full_info, bool singleprec,
      const std::string &filename, int first, int last);

    void write (const arr<pointing> &ptg);
    void write (const arr<pointing> &ptg, const arr<double> &heading);
  };

class tod_writer: public generic_writer
  {
  private:
    int todcol;

  public:
    tod_writer (const std::string &filename, int first, int last);

    void write (const arr<double> &toddata);
  };

class time_writer: public generic_writer
  {
  private:
    int timecol;

  public:
    time_writer (const std::string &filename);

    void write (const arr<double> &times);
  };

class quat_writer: public generic_writer
  {
  private:
    int qxcol, qycol, qzcol, qwcol;

  public:
    quat_writer (const std::string &filename, const std::string &fpdb);

    void write (const arr<quaternion> &quatdata);
  };

class flag_writer: public generic_writer
  {
  private:
    int flagcol;

  public:
    flag_writer (const std::string &filename);

    void write (const arr<bool> &flagdata);
  };

#endif
