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

#include "c_utils.h"
#include "sharp_lowlevel.h"

#define CONCAT(a,b) a ## b

#define FLT double
#define DP_TRANS
#define X(arg) CONCAT(sharpd_,arg)
#include "sharp_f_inc.c"
#undef DP_TRANS
#undef FLT
#undef X

#define FLT float
#define X(arg) CONCAT(sharps_,arg)
#include "sharp_f_inc.c"
#undef FLT
#undef X

#undef CONCAT
