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
 *  Copyright (C) 2003-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "dipole.h"
#include "pointing.h"
#include "vec3.h"
#include "sat_info.h"
#include "arr.h"
#include "focalplane_db.h"
#include "lsconstants.h"
#include "paramfile.h"

using namespace std;

namespace {

//! Colatitude of the solar system motion relative to CMB
//! (ecliptical coordinates) (Hinshaw et al. 2009, ApJS, 180, 225).
const double solsysdir_ecl_theta = 1.765248346;

//! Longitude of the solar system motion relative to CMB
//! (ecliptical coordinates) (Hinshaw et al. 2009, ApJS, 180, 225).
const double solsysdir_ecl_phi = 2.995840906;

//! Speed of the solar system motion relative to CMB in m/s
//! (Hinshaw et al. 2009, ApJS, 180, 225).
const double solsysspeed = 369000.0;

inline double T_ant (double t_thermo, double hnydk)
  { return hnydk/(exp(hnydk/t_thermo)-1); }

} // unnamed namespace

Dipole::Dipole (Sat_Info &info, focalplane_db &focal, paramfile &params,
                bool source_dipole, bool source_fsldp)
  : satinfo(info), fpdb(focal), do_dipole(source_dipole),
    do_fsldp(source_fsldp)
  {
  string speedstr = params.find<string>("dipole_speed","TOTAL");
  speed = TOTAL;
  if (speedstr=="TOTAL") speed = TOTAL;
  else if (speedstr=="SOLSYS") speed = SOLSYS;
  else if (speedstr=="SATELLITE") speed = SATELLITE;
  else planck_fail ("Incorrect value '" + speedstr + "' for dipole_speed");
  thermotemp = params.find<bool>("dipole_thermotemp",false);
  outputtype = params.find<int>("dipole_type",2);
  planck_assert ((outputtype>=1) && (outputtype<=7), "dipole: wrong type");
  solsysdir_v = pointing(solsysdir_ecl_theta, solsysdir_ecl_phi).to_vec3();
  dip_norm = params.find<double>("dip_norm",1);
  }

void Dipole::Add_Intensities (const string &det_id, const arr<vec3> &vdetpt,
  const arr<vec3> &vsldppt, arr<double> &intensity) const
  {
  planck_assert (multiequal(vdetpt.size(),intensity.size(),vsldppt.size()),
    "array size mismatch");

  if (do_fsldp)
    {
    double detfreq = fpdb.getValue<double>(det_id,"nu_cen");
    double fact = 0.5 * (1 + fpdb.getValue<double>(det_id,"epsilon"))*dip_norm;
    double hnydk = hPlanck*detfreq/kBoltzmann;
    double x = hnydk/tcmb;
    double expx = exp(x);
    double planck = x/(expx-1);
    double thermo2ant = expx *planck*planck;
    vec3 vsat = satinfo.velocity();

    vec3 detspeed;
    if (speed==TOTAL)     detspeed = solsysdir_v*solsysspeed + vsat;
    if (speed==SOLSYS)    detspeed = solsysdir_v*solsysspeed;
    if (speed==SATELLITE) detspeed = vsat;
    vec3 speed_dir = detspeed; speed_dir.Normalize();
    double detspeed_value = detspeed.Length();

    double beta = detspeed_value/speedOfLight;

#pragma omp parallel
{
    int m, sz=vdetpt.size();
#pragma omp for schedule(static)
    for (m=0; m<sz; ++m)
      {
      /* Start by adding far sidelobe dipole since this is always a plain
         dipole irrespective of the outputtype switch */
      double fsldp_thermo_dp = tcmb * beta * dotprod(speed_dir,vsldppt[m]);
      intensity[m] += fact *(thermotemp ?
        fsldp_thermo_dp :
        thermo2ant * fsldp_thermo_dp);
      }
}
    }

  /* now add contribution through main beam */
  Add_Intensities (det_id, vdetpt, intensity);
  }

void Dipole::Add_Intensities (const string &det_id, const arr<vec3> &vdetpt,
  arr<double> &intensity) const
  {
  planck_assert (vdetpt.size()==intensity.size(), "array size mismatch");

  if (!do_dipole) return;

  double detfreq = fpdb.getValue<double>(det_id,"nu_cen");
  double fact = 0.5 * (1 + fpdb.getValue<double>(det_id,"epsilon"))*dip_norm;
  double hnydk = hPlanck*detfreq/kBoltzmann;
  double tcmbant = T_ant(tcmb, hnydk);

  vec3 vsat = satinfo.velocity();

  vec3 detspeed;
  if (speed==TOTAL)     detspeed = solsysdir_v*solsysspeed + vsat;
  if (speed==SOLSYS)    detspeed = solsysdir_v*solsysspeed;
  if (speed==SATELLITE) detspeed = vsat;
  vec3 speed_dir = detspeed; speed_dir.Normalize();
  double detspeed_value = detspeed.Length();

  double beta = detspeed_value/speedOfLight;
  double gamma = 1/sqrt(1-beta*beta);

#pragma omp parallel
{
  int m, sz=vdetpt.size();
#pragma omp for schedule(static)
  for (m=0; m<sz; ++m)
    {
    double cosdir = dotprod(speed_dir,vdetpt[m]);

    switch (outputtype)
      {
      case 1:
        {
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        intensity[m] += fact * (thermotemp ?
          t_thermo_all :
          T_ant(t_thermo_all,hnydk));
        }
        break;
      case 2:
        {
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        intensity[m] += fact * (thermotemp ?
          (t_thermo_all-tcmb) :
          (T_ant(t_thermo_all,hnydk)-tcmbant));
        }
        break;
      case 3:
        {
        double t_thermo_mp_dp = tcmb * (1+beta*cosdir);
        intensity[m] += fact * (thermotemp ?
          t_thermo_mp_dp :
          T_ant(t_thermo_mp_dp, hnydk));
        }
        break;
      case 4:
        {
        double t_thermo_mp_dp = tcmb * (1+beta*cosdir);
        intensity[m] += fact * (thermotemp ?
          (t_thermo_mp_dp-tcmb) :
          (T_ant(t_thermo_mp_dp, hnydk)-tcmbant));
        }
        break;
      case 5:
        {
        double t_thermo_mp_dp = tcmb * (1+beta*cosdir);
        double t_thermo_mp_dp_qp = tcmb *
          (1+beta*cosdir + beta*beta*(cosdir*cosdir-0.5));
        intensity[m] += fact * (thermotemp ?
          (t_thermo_mp_dp_qp-t_thermo_mp_dp) :
          (T_ant(t_thermo_mp_dp_qp,hnydk)-T_ant(t_thermo_mp_dp,hnydk)));
        }
        break;
      case 6:
        {
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        double t_thermo_mp_dp_qp = tcmb *
          (1+beta*cosdir + beta*beta*(cosdir*cosdir-0.5));
        intensity[m] += fact * (thermotemp ?
          (t_thermo_all-t_thermo_mp_dp_qp) :
          (T_ant(t_thermo_all,hnydk)-T_ant(t_thermo_mp_dp_qp,hnydk)));
        }
        break;
      case 7:
        {
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        double t_thermo_mp_dp = tcmb * (1+beta*cosdir);
        intensity[m] += fact * (thermotemp ?
          (t_thermo_all-t_thermo_mp_dp) :
          (T_ant(t_thermo_all,hnydk)-T_ant(t_thermo_mp_dp,hnydk)));
        }
        break;
      }
    }
}
  }
