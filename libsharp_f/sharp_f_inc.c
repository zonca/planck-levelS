/*
   This file is part of the "libsharp_f" component of the Planck simulation
   package.

   libsharp_f is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with libsharp_f; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
   Copyright (C) 2012-2013 Max-Planck-Society
   \author Martin Reinecke
*/

void X(do_job_f) (int ntheta, const int *nphi, const int *theta_start,
  const int *phi_stride, const double *phi_0, const double *theta,
  const double *weight, FLT *map, int map_stride, int lmax, int mmax, FLT *alm,
  int alm_stride, const int *alm_m_start, int alm_l_stride, int jobtype, int add);

void X(do_job_f) (int ntheta, const int *nphi, const int *theta_start,
  const int *phi_stride, const double *phi_0, const double *theta,
  const double *weight, FLT *map, int map_stride, int lmax, int mmax, FLT *alm,
  int alm_stride, const int *alm_m_start, int alm_l_stride, int jobtype, int add)
  {
  sharp_geom_info *tinfo;
  sharp_alm_info *alms;
  ptrdiff_t *theta_start2=RALLOC(ptrdiff_t,ntheta),
            *alm_m_start2=RALLOC(ptrdiff_t,mmax+1);
  COPY_ARRAY(theta_start,theta_start2,0,ntheta);
  COPY_ARRAY(alm_m_start,alm_m_start2,0,mmax+1);

  CHECK_STACK_ALIGN(8);

  sharp_make_geom_info (ntheta, nphi, theta_start2, phi_stride, phi_0, theta,
    weight, &tinfo);
  sharp_make_alm_info(lmax,mmax,alm_l_stride,alm_m_start2,&alms);

  int flags=0;
#ifdef DP_TRANS
  flags|=SHARP_DP;
#endif
  if (add) flags |= SHARP_ADD;

  switch(jobtype)
    {
    case 1:
      sharp_execute (SHARP_MAP2ALM, 0, &alm, &map, tinfo, alms, 1, flags, NULL, NULL);
      break;
    case 2:
      {
      sharp_execute (SHARP_MAP2ALM, 0, &alm, &map, tinfo, alms, 1, flags, NULL, NULL);
      FLT *mp[2];
      mp[0]=map+map_stride; mp[1]=map+2*map_stride;
      FLT *ap[2];
      ap[0]=alm+2*alm_stride; ap[1]=alm+4*alm_stride;
      sharp_execute (SHARP_MAP2ALM, 2, &ap[0], &mp[0], tinfo, alms, 1, flags, NULL, NULL);
      }
      break;
    default:
      UTIL_FAIL("bad job type in do_job()");
      break;
    }

  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);
  DEALLOC(theta_start2);
  DEALLOC(alm_m_start2);
  }
