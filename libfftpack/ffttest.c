/*
 *  This file is part of libfftpack.
 *
 *  libfftpack is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libfftpack is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libfftpack; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libfftpack is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Test codes for libfftpack.
 *
 *  Copyright (C) 2004-2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ls_fft.h"
#include "fftpack.h"
#include "ls_fft.c"
#include "bluestein.h"

#define maxlen 4096

static void fill_random (double *data, double *odata, size_t length)
  {
  size_t m;
  for (m=0; m<length; ++m)
    data[m] = odata[m] = rand()/(RAND_MAX+1.0)-0.5;
  }
static void fill_random_c (double *data, double *odata, size_t length)
  {
  size_t m;
  for (m=0; m<length; ++m)
    {
    data[2*m] = odata[2*m] = rand()/(RAND_MAX+1.0)-0.5;
    data[2*m+1] = odata[2*m+1] = 0.;
    }
  }

static void normalize (double *data, size_t length, double norm)
  {
  size_t m;
  for (m=0; m<length; ++m)
    data[m] /= norm;
  }

static double errcalc (double *data, double *odata, size_t length)
  {
  size_t m;
  double sum, errsum;
  sum = errsum = 0;
  for (m=0; m<length; ++m)
    {
    errsum += (data[m]-odata[m])*(data[m]-odata[m]);
    sum += odata[m]*odata[m];
    }
  return sqrt(errsum/sum);
  }

static void test_real(void)
  {
  double data[2*maxlen], odata[2*maxlen];
  int length;
  double err;
  const double epsilon=3e-15;
  real_plan plan;
  for (length=1; length<=maxlen; ++length)
    {
    fill_random_c (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_c(plan, data);
    real_plan_backward_c(plan, data);
    normalize (data, 2*length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, 2*length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_fftpack (plan, data);
    real_plan_backward_fftpack (plan, data);
    normalize (data, length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_fftw (plan, data);
    real_plan_backward_fftw (plan, data);
    normalize (data, length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_fftw (plan, data);
    halfcomplex2fftpack (data,plan->length);
    real_plan_backward_fftpack (plan, data);
    normalize (data, length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_fftpack (plan, data);
    fftpack2halfcomplex (data,plan->length);
    real_plan_backward_fftw (plan, data);
    normalize (data, length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  }

static void test_complex(void)
  {
  double data[2*maxlen], odata[2*maxlen];
  int length;
  double err;
  const double epsilon=3e-15;
  complex_plan plan;
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, 2*length);
    plan = make_complex_plan (length);
    complex_plan_forward(plan, data);
    complex_plan_backward(plan, data);
    normalize (data, 2*length, length);
    kill_complex_plan (plan);
    err = errcalc (data, odata, 2*length);
    if (err>epsilon) printf("problem at complex length %i: %e\n",length,err);
    }
  }

int main(void)
  {
  test_real();
  test_complex();
  return 0;
  }
