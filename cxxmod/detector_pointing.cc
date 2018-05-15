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
 *  Copyright (C) 2003-2014 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "detector_pointing.h"
#include "rotmatrix.h"
#include "sat_info.h"
#include "focalplane_db.h"
#include "arr.h"
#include "pointing.h"
#include "geom_utils.h"
#include "lsconstants.h"

using namespace std;

Detector_Pointing::Detector_Pointing (paramfile &params, Sat_Info &sptg,
  const focalplane_db &fpdb, const string &det_id, bool aberration_)
  : satptg(sptg), aberration(aberration_)
  {
  if (params.find<bool>("detpt_wcorr",false))
    wcorr = new wobble_correction(params);
  double theta_b = fpdb.theta_b()*degr2rad,
         theta_uv = fpdb.getValue<double>(det_id,"theta_uv")*degr2rad,
         phi_uv = fpdb.getValue<double>(det_id,"phi_uv")*degr2rad,
         psi_uv = fpdb.getValue<double>(det_id,"psi_uv")*degr2rad,
         sldp_x = fpdb.getValue<double>(det_id,"sldp_x"),
         sldp_y = fpdb.getValue<double>(det_id,"sldp_y"),
         sldp_z = fpdb.getValue<double>(det_id,"sldp_z");
  if (params.find<bool>("detpt_ptcor",false))
    ptcorr = new ptcor(params,theta_b);

  /* Construct rotation matrix that makes an active transform
     from spin axis frame to detector frame: */
  rotmatrix mat1,mat2,mat3;
  mat1.Make_Axis_Rotation_Transform(vec3(0,1,0),halfpi-theta_b);
  mat2.Make_Axis_Rotation_Transform(vec3(-sin(phi_uv),cos(phi_uv),0),theta_uv);
  mat3.Make_Axis_Rotation_Transform(vec3(0,0,1),psi_uv);
  rotmatrix ax2det=mat1*(mat2*mat3);

  /* Find direction of detector in spin axis frame
     by transforming vector along spin axis: */
  direction = ax2det.Transform(vec3(0,0,1));
  /* Find direction of sidelobe dipole vector in spin axis frame by
     applying the same transform: */
  vsldp = ax2det.Transform(vec3(sldp_x,sldp_y,sldp_z));
  /* Transform spin axis "x" vector to detector frame as a way of
     recording detector orientation on the sky */
  xdir = ax2det.Transform(vec3(1,0,0));
  }

void Detector_Pointing::Get_Pointings (const arr<double> &times,
      arr<pointing> *detpt, arr<vec3> *vdetpt, arr<double> *heading,
      arr<vec3> *vsldppt) const
  {
  int num_samp=times.size();
  if (detpt) detpt->alloc(num_samp);
  if (vdetpt) vdetpt->alloc(num_samp);
  if (heading) heading->alloc(num_samp);
  if (vsldppt) vsldppt->alloc(num_samp);
  if (num_samp==0) return;
  vec3 vcorr = satptg.velocity()*(1./speedOfLight);
  rotmatrix wobble;
  if (wcorr) wobble = wcorr->get_matrix(satptg.curPeriod());

#pragma omp parallel
{
  rotmatrix trans;
  vec3 vdpt, vxpt;
  int i;
#pragma omp for schedule (static)
  for (i=0; i<num_samp; ++i)
    {
    satptg.getTransform(times[i], trans);
    if (wcorr) trans = trans*wobble;
    if (ptcorr)
      {
      rotmatrix mptcor;
      mptcor = ptcorr->get_matrix(times[i]);
      trans = trans*mptcor;
      }
    trans.Transform (direction,vdpt);
    if (aberration) vdpt=(vdpt+crossprod(vdpt,crossprod(vdpt,vcorr))).Norm();
    if (vdetpt) (*vdetpt)[i]=vdpt;
    if (detpt) (*detpt)[i].from_vec3(vdpt);
    if (heading)
      {
      trans.Transform (xdir,vxpt);
      double hdg = orientation (vdpt,vxpt)+pi;
      if (hdg>twopi) hdg-=twopi;
      (*heading)[i] = hdg;
      }
    if (vsldppt) trans.Transform(vsldp,(*vsldppt)[i]);
    }
}
  }
