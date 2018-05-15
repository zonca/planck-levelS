/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file svd_c.h
 *  Functions for singular value decomposition
 *
 *  Copyright (C) 2008, 2009 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SVD_C_H
#define PLANCK_SVD_C_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
  {
  double **u, **v, *w;
  int m, n;
  } svd_obj;

void svd_init (double **matrix, double epsilon, int m, int n,
  svd_obj *obj);
void svd_solve (const svd_obj *obj, double *b);
void svd_destroy (svd_obj *obj);

#ifdef __cplusplus
}
#endif

#endif
