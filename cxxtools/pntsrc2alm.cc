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
 *  Copyright (C) 2009-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm.h"
#include "alm_dmcio.h"
#include "xcomplex.h"
#include "pointing.h"
#include "sharp_cxx.h"
#include "planck_rng.h"
#include "paramfile.h"
#include "iohandle_current.h"
#include "mpi_support.h"
#include "io_utils.h"
#include "lsconstants.h"
#include "announce.h"
#include "string_utils.h"

using namespace std;

namespace {

struct ptsrc
  {
  double theta, phi, i, q, u;
  };

void readSources (paramfile &params, bool polarisation,
  arr<ptsrc> &src)
  {
  string infile = params.find<string>("pntsrc_file","");
  arr<pointing> pos;
  arr<double> flux, polang, polfrac;
  arr<string> name;
  readPointSourcesNew (infile, polarisation, pos, flux, polang,
    polfrac, name, mpiMgr.num_ranks(), mpiMgr.rank());

  src.alloc(flux.size());
  for (tsize j=0; j<src.size(); ++j)
    {
    src[j].phi = fmodulo(pos[j].phi, twopi);
    src[j].theta = pos[j].theta;
    src[j].i = flux[j];
    if (polarisation)
      {
      double vpol = src[j].i*polfrac[j];
      src[j].q = vpol*cos(2*polang[j]);
      src[j].u = vpol*sin(2*polang[j]);
      }
    else
      src[j].q = src[j].u = 0;
    }
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("pntsrc2alm", argc, argv, mpiMgr.master());

  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams(mpiMgr.master()));

  bool polarisation = params.find<bool>("polarisation");
  arr<ptsrc> src;
  readSources (params, polarisation, src);

  int lmax = params.find<int>("lmax");
  int mmax = lmax;
  string outfile = params.find<string>("outFile");

  Alm<xcomplex<double> > almT(lmax,mmax), almG, almC;
  if (polarisation)
    { almG.Set(lmax,mmax); almC.Set(lmax,mmax); }

  sharp_cxxjob<double> job;
  arr<double> theta(src.size()), phi0(src.size()),
    si(src.size()), sq(src.size()), su(src.size());
  arr<int> nph(src.size(),1);
  arr<ptrdiff_t> ofs(src.size());
  for (tsize i=0; i<src.size(); ++i)
    {
    theta[i]=src[i].theta;
    phi0[i]=src[i].phi;
    ofs[i]=i;
    si[i]=src[i].i;
    sq[i]=src[i].q;
    su[i]=src[i].u;
    }
  job.set_general_geometry (src.size(), &nph[0], &ofs[0], &nph[0], &phi0[0],
    &theta[0], NULL);
  job.set_triangular_alm_info (lmax, lmax);
  job.map2alm(&si[0],&almT(0,0).re,false);
  if (polarisation)
    job.map2alm_spin(&sq[0],&su[0],&almG(0,0).re,&almC(0,0).re,2,false);

  if (mpiMgr.num_ranks()>1)
    {
    tsize nalm = almT.Alms().size();
    Alm<xcomplex<double> > tmp (almT);
    mpiMgr.allreduceRaw(&tmp(0,0).re,&almT(0,0).re,2*nalm,MPI_Manager::Sum);
    if (polarisation)
      {
      tmp=almG;
      mpiMgr.allreduceRaw(&tmp(0,0).re,&almG(0,0).re,2*nalm,MPI_Manager::Sum);
      tmp=almC;
      mpiMgr.allreduceRaw(&tmp(0,0).re,&almC(0,0).re,2*nalm,MPI_Manager::Sum);
      }
    }

  if (mpiMgr.master())
    {
    if (polarisation)
      write_Alm_to_dmc (outfile,almT,almG,almC,lmax,mmax,true);
    else
      write_Alm_to_dmc (outfile,almT,lmax,mmax,true);
    }

PLANCK_DIAGNOSIS_END
  }
