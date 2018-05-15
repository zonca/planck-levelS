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

#include "lsconstants.h"
#include "arr.h"
#include "announce.h"
#include "focalplane_db.h"
#include "paramfile.h"
#include "pointing.h"
#include "vec3.h"
#include "rotmatrix.h"
#include "detector_pointing.h"
#include "interpolator.h"
#include "interpolator2.h"
#include "sat_info.h"
#include "sampler.h"
#include "dipole.h"
#include "oofnoise.h"
#include "satquat.h"
#include "PointSourceConvolver.h"
#include "writers.h"
#include "healpix_map.h"
#include "healpix_map_dmcio.h"
#include "iohandle_current.h"
#include "zle.h"
#include "repointing.h"
#include "levels_facilities.h"
#include "walltimer.h"
#include "mpi_support.h"

#include "wobble_correction.h"

using namespace std;

namespace {

#define TIME(op,timer) \
  do { tstack_push(timer); op; tstack_pop(timer); } while(0)

void parallel_add (arr<double> &a1, const arr<double> &a2)
  {
#pragma omp parallel
{
  int m, sz=a1.size();
#pragma omp for schedule (static)
  for (m=0; m<sz; ++m)
    a1[m] += a2[m];
}
  }

void parallel_zero (arr<double> &a)
  {
#pragma omp parallel
{
  int m, sz=a.size();
#pragma omp for schedule (static)
  for (m=0; m<sz; ++m)
    a[m] = 0.;
}
  }

void get_params (const string &name, paramfile &par_file,
  string &outfile, bool &write)
  {
  outfile = par_file.find<string> (name+string("_file"),"");
  write = (outfile!="");
  }

void  ingestTimeData (int64 first_sample, double start_time_mission,
  double  start_time, double f_samp, generic_writer *writer)
  {
  if (!writer) return;
  writer->setKey("first_sample", first_sample);
  writer->setKey("start_time_mission", start_time_mission);
  writer->setKey("start_time", start_time);
  writer->setKey("f_samp", f_samp);
  }

void replaceKey (const string &key, paramfile &params, int idx)
  {
  if (params.param_present(key))
    params.setParam(key,params.find<string>(key)+intToString(idx,5));
  }

void special_mpi_treatment_for_helsinki (paramfile &params)
  {
  if (!params.param_present("helsinki_mpi_mode")) return;
  if (!params.find<bool>("helsinki_mpi_mode",false)) return;
  if (mpiMgr.num_ranks()<2) return;
  planck_assert(!params.find<bool>("source_oof",false),
    "cannot compute 1/f noise in Helsinki MPI mode");

// in Helsinki MPI mode, first_pointing and last_pointing MUST be specified!
  int first = params.find<int>("first_pointing");
  int last = params.find<int>("last_pointing");
  planck_assert ((first>=1)&&(last>=1)&&(first<=last),
    "inconsistent first and last pointing periods!");

  int64 myfirst, mylast;
  mpiMgr.calcShare (first, last+1, myfirst, mylast);
  --mylast;
  params.setParam("first_pointing",myfirst);
  params.setParam("last_pointing",mylast);
  replaceKey ("sampinfo_file",params,myfirst);
  replaceKey ("detpt_file",params,myfirst);
  replaceKey ("tod_file",params,myfirst);
// add calls for all other output file names that need to be replaced
  }

} // unnamed namespace

int multimod_module(int argc, const char **argv)
  {
  module_startup ("multimod", argc, argv);
  tstack_push("multimod");
  tstack_push("startup");
  iohandle_current::Manager mng (argc, argv);
  paramfile par_file(mng.getParams());
  special_mpi_treatment_for_helsinki(par_file);

  string detpt_file, cnt_file, map_file, tod_file, quat_file, index_file,
    time_file, repointing_flag_file, sampinfo_file;
  bool write_detpt, write_cnt, write_map, write_tod, write_quat, write_index,
    write_time, write_repointing_flag, write_sampinfo;
  get_params("detpt",par_file,detpt_file,write_detpt);
  get_params("cnt",par_file,cnt_file,write_cnt);
  get_params("map",par_file,map_file,write_map);
  get_params("tod",par_file,tod_file,write_tod);
  get_params("quaternions",par_file,quat_file,write_quat);
  get_params("timestamp",par_file,time_file,write_time);
  get_params("index",par_file,index_file,write_index);
  get_params("repointing_flag",par_file,repointing_flag_file,
    write_repointing_flag);
  get_params("sampinfo",par_file,sampinfo_file,write_sampinfo);
  bool detpt_sp = par_file.find<bool> ("single_precision_detpt", false);
  bool source_mixed = par_file.find<bool> ("source_mixed", false);
  bool source_mixed2 = par_file.find<bool> ("source_mixed2", false);
  bool source_dipole = par_file.find<bool> ("source_dipole", false);
  bool source_fsldp  = par_file.find<bool> ("source_fsldp", false);
  bool source_oof = par_file.find<bool> ("source_oof", false);
  bool source_pntsrc = par_file.find<bool> ("source_pntsrc", false);
  bool source_zle = par_file.find<bool> ("source_zle", false);

  bool do_detpt1 = false,
       do_detpt2 = write_detpt,
       do_cnt    = write_cnt,
       do_map    = write_map,
       do_rings  = false,
       do_tod    = write_tod,
       do_quat   = write_quat;

// dependencies
  if (do_map) do_detpt2 = do_cnt = do_tod = true;
  if (do_tod &&
      (source_mixed || source_mixed2 || source_dipole || source_fsldp ||
       source_pntsrc || source_zle))
    do_rings = true;
  if (do_cnt) do_detpt2 = true;
  if (do_rings &&
      (source_mixed || source_mixed2 || source_dipole || source_fsldp ||
       source_pntsrc || source_zle))
    do_detpt1=true;

  if (do_detpt1) cout << "doing detpt1" << endl;
  if (do_rings) cout << "doing rings" << endl;
  if (do_detpt2) cout << "doing detpt2" << endl;
  if (do_cnt) cout << "doing cnt" << endl;
  if (do_tod) cout << "doing tod" << endl;
  if (do_map) cout << "doing map" << endl;
  if (do_quat) cout << "doing quaternions" << endl;

  string det_id = par_file.find<string>("detector_id");
  bool nominal_pointing = par_file.find<bool>("nominal_pointing");

  focalplane_db fpdb(par_file);

  safe_ptr<Sat_Info> satinfo;
  string sattype = par_file.find<string>("satinfo_type","SIMMISSION");
  if (sattype=="SIMMISSION")
    satinfo = new Sat_Info_Simmission (nominal_pointing, par_file);
  else if (sattype=="LFI")
    satinfo = new Sat_Info_LFI (nominal_pointing, par_file);
  else if (sattype=="HFI")
    satinfo = new Sat_Info_HFI (nominal_pointing, par_file);
  else
    planck_fail ("Incorrect value '" +sattype+ "' for satellite pointing type");

  int first = par_file.find<int>("first_pointing",1);
  if (first==-1) first=1;
  --first;
  int last = par_file.find<int>("last_pointing",satinfo->numPeriods());
  if (last==-1) last=satinfo->numPeriods();
  --last;
  planck_assert ((first>=0)&&(last>=0)&&(first<=last)
    &&(last<satinfo->numPeriods()),
    "inconsistent first and last pointing periods!");

  Detector_Pointing pntgs(par_file,*satinfo,fpdb,det_id,
    par_file.find("detpt_aberration",false));

  double o_factor = 1.;
  if (do_rings || do_detpt1)
    {
    o_factor = par_file.find<double>("oversampling_factor");
    planck_assert(o_factor>0,"incorrect value for oversampling_factor");
    }

  double f_samp_det=fpdb.getValue<double>(det_id,"f_samp");
  double f_samp2 = (sattype=="HFI"||sattype=="LFI") ?
    -f_samp_det : f_samp_det;
  double f_samp1=f_samp_det*o_factor;

  bool bypass_sampler = par_file.find<bool>("bypass_sampler",false);
  if (bypass_sampler)
    {
    f_samp1=f_samp2;
    planck_assert(o_factor==1.,
      "if bypass_sampler==true, oversampling_factor must be 1.");
    }

  double calibration_factor = 1;
  if (do_tod)
    if (par_file.find<bool> ("calibrate_signal"))
      calibration_factor = 2./(1.+fpdb.getValue<double>(det_id,"epsilon"));

  safe_ptr<Interpolator> inter;
  double psi_pol = degr2rad*fpdb.getValue<double>(det_id,"psi_pol");
  if (do_rings && source_mixed)
    inter = new Interpolator(par_file, psi_pol);
  safe_ptr<Interpolator2> inter2;
  if (do_rings && source_mixed2)
    inter2 = new Interpolator2(par_file, fpdb, det_id, psi_pol);
  safe_ptr<Dipole> dipole;
  if (do_rings && (source_dipole || source_fsldp))
    dipole = new Dipole(*satinfo, fpdb, par_file, source_dipole, source_fsldp);
  safe_ptr<ZLE> zle;
  if (do_rings && source_zle)
    zle = new ZLE(fpdb,det_id,par_file);
  safe_ptr<Satquat> satquat;
  if (do_quat)
    satquat = new Satquat(*satinfo);

  safe_ptr<psc::PointSourceConvolver> pntsrcconv;
  if (do_rings && source_pntsrc)
    pntsrcconv = new psc::PointSourceConvolver(par_file,fpdb,det_id,*satinfo);

  safe_ptr<Oofnoise> oof;
  if (do_tod && source_oof)
    {
    int num_real=1;
    if (nominal_pointing)
      num_real = par_file.find<int>("oof_num_real",60);
    oof = new Oofnoise(fpdb,det_id,num_real,par_file);
    }

  safe_ptr<repointing> rept;
  if (write_repointing_flag)
    rept = new repointing (par_file, nominal_pointing);

  safe_ptr<Sampler> sampler;
  if (do_tod && do_rings && (!bypass_sampler))
    sampler = new Sampler(par_file,fpdb,det_id,5,f_samp1);
  Healpix_Map<double> mapper;
  Healpix_Map<int> mapcnt;

  int mapnside=-1;
  if (do_map || do_cnt)
    mapnside = par_file.find<int>("nside");

  if (do_cnt) { mapcnt.SetNside(mapnside,RING); mapcnt.fill(0); }
  if (do_map) { mapper.SetNside(mapnside,RING); mapper.fill(0); }

  safe_ptr<tod_writer> tod_out;
  if (write_tod)
    tod_out = new tod_writer (tod_file,first+1,last+1);

  safe_ptr<detpt_writer> pt_writer;
  if (write_detpt)
    pt_writer = new detpt_writer (true, detpt_sp, detpt_file, first+1, last+1);
  safe_ptr<quat_writer> quat_out;
  if (write_quat)
    quat_out = new quat_writer (quat_file,
      par_file.find<string>("focalplane_db"));
  safe_ptr<time_writer> time_out;
  if (write_time)
    time_out = new time_writer (time_file);
  safe_ptr<flag_writer> rept_out;
  if (write_repointing_flag)
    rept_out = new flag_writer (repointing_flag_file);

  arr<pointing> ptarr1,ptarr2_x;
  arr<vec3> vptarr1,vptarr2_x,vsldpptarr;
  arr<double> heading1,heading2_x;
  bool separate_pointings =
    (!do_detpt1) || (!approx(o_factor,1.)) || (f_samp1*f_samp2<0);
  arr<pointing> &ptarr2 (separate_pointings ? ptarr2_x : ptarr1);
  arr<vec3> &vptarr2 (separate_pointings ? vptarr2_x : vptarr1);
  arr<double> &heading2 (separate_pointings ? heading2_x : heading1);
  arr<double> intensity, tod, tmprings, tmptod, timearr;
  arr<quaternion> quatarr;
  arr<bool> rept_flag_arr;
  arr<int64> nsamp(last-first+1);
  arr<double> tstart(last-first+1),tend(last-first+1);
  arr<double> times1, times2;

  {
  int64 samp0=satinfo->startIndex(first,f_samp2);
  double start_time_mission=satinfo->startTime(0);
  double start_time=satinfo->startTime(first);
  ingestTimeData (samp0, start_time_mission, start_time, f_samp2, tod_out);
  ingestTimeData (samp0, start_time_mission, start_time, f_samp2, pt_writer);
  ingestTimeData (samp0, start_time_mission, start_time, f_samp2, quat_out);
  ingestTimeData (samp0, start_time_mission, start_time, f_samp2, time_out);
  }

  tstack_replace("startup","period_loop");

  for (int period=first; period<=last; ++period)
    {
    tstack_push("satinfo");
    TIME(satinfo->gotoPeriod(period),"goto_period");

    tstack_push("gettimes");
    satinfo->getTimeStamps (period, f_samp1, times1);
    satinfo->getTimeStamps (period, f_samp2, times2);
    tstack_pop("gettimes");
    tstack_pop("satinfo");
    int numsamples1 = times1.size(),
        numsamples2 = times2.size();
    nsamp[period-first]=numsamples2;
    if ((numsamples1<=0)||(numsamples2<=0))
      {
      cout << "WARNING: empty pointing period detected: " << period+1 << endl;
      tstart[period-first] = tend[period-first] = satinfo->startTime(period);
      }
    else
      {
      tstart[period-first]=times2[0];
      tend[period-first]=times2[numsamples2-1];
      }

    if (do_detpt1)
      TIME(pntgs.Get_Pointings (times1, ptarr1, vptarr1, heading1, vsldpptarr),
           "detpt");

    if (do_rings)
      {
      intensity.alloc(numsamples1);
      parallel_zero(intensity);
      if (source_mixed)
        TIME(inter->Add_Intensities(ptarr1, heading1, intensity),"interpol");
      if (source_mixed2)
        TIME(inter2->Add_Intensities(ptarr1, heading1, intensity),"interpol2");
      if (source_dipole || source_fsldp)
        TIME(dipole->Add_Intensities(det_id, vptarr1, vsldpptarr, intensity),
             "dipole");
      if (source_zle)
        TIME(zle->Add_Intensities(ptarr1, vptarr1, times1, intensity),"zle");

      if (source_pntsrc)
        {
        tstack_push("pntsrc");
        pntsrcconv->Get_Intensities(period,ptarr1,vptarr1,heading1,tmprings);
        parallel_add (intensity,tmprings);
        tstack_pop("pntsrc");
        }
      }

    if (do_detpt2 && separate_pointings)
      TIME(pntgs.Get_Pointings (times2, ptarr2, vptarr2, heading2),"detpt");

    if (do_tod)
      {
      tod.alloc(numsamples2);
      parallel_zero(tod);

      if (source_oof)
        {
        tstack_push("oof");
        oof->Get_Noise (period,numsamples2,tmptod);
        parallel_add(tod,tmptod);
        tstack_pop("oof");
        }

      if (do_rings)
        {
        if (bypass_sampler)
          parallel_add(tod,intensity);
        else
          {
          tstack_push("sampler");
          sampler->sample(intensity,times1,times2,tmptod);
          parallel_add(tod,tmptod);
          tstack_pop("sampler");
          }
        }
      }

    for (tsize m=0; m<tod.size(); ++m)
      tod[m]*=calibration_factor;

    if (write_tod)
      tod_out->write(tod);

    if (write_detpt)
      pt_writer->write(ptarr2, heading2);

    if (do_quat)
      TIME(satquat->Get_Quaternions (times2, quatarr),"quaternions");
    if (write_quat)
      quat_out->write(quatarr);

    if (write_time)
      time_out->write(times2);

    if (write_repointing_flag)
      {
      rept->get_data(times2, rept_flag_arr);
      rept_out->write(rept_flag_arr);
      }

    if (do_cnt)
      {
      tstack_push("cnt/map");
      arr<int> pix(ptarr2.size());
#pragma omp parallel
{
      int m, sz=pix.size();
#pragma omp for schedule(static)
      for (m=0; m<sz; ++m)
        pix[m] = mapcnt.zphi2pix(vptarr2[m].z,ptarr2[m].phi);
}
      for (tsize m=0; m<ptarr2.size(); ++m)
        {
        ++mapcnt[pix[m]];
        if (do_map) mapper[pix[m]] += tod[m];
        }
      tstack_pop("cnt/map");
      }
    }

  if (pntsrcconv) pntsrcconv->Write_Hitsfile();

  if (write_cnt)
    write_Healpix_hitmap_to_dmc (cnt_file,mapcnt);
  if (write_map)
    {
    for (int m=0; m<mapper.Npix(); ++m)
      mapper[m] = (mapcnt[m]>0) ? (mapper[m]/mapcnt[m]) : Healpix_undef;
    write_Healpix_map_to_dmc (map_file,mapper);
    }
  if (write_index)
    {
    arr<int64> first_sample(nsamp.size()), last_sample(nsamp.size());
    first_sample[0]=0;
    last_sample[0]=nsamp[0]-1;
    for (tsize i=1; i<nsamp.size(); ++i)
      {
      first_sample[i]=last_sample[i-1]+1;
      last_sample[i]=first_sample[i]+nsamp[i]-1;
      }
    safe_ptr<iohandle> out (HandleManager.createObject
      (index_file,"table.LS_toi_index_table"));
    out->appendColumn("first_sample",first_sample);
    out->appendColumn("last_sample",last_sample);
    }
  if (write_sampinfo)
    {
    arr<int> pid(nsamp.size()), nsamp2(nsamp.size());
    for (tsize i=0; i<pid.size(); ++i)
      {
      pid[i]=i+first+1;
      nsamp2[i]=nsamp[i];
      tstart[i]*=1e6;
      tend[i]*=1e6;
      }
    safe_ptr<iohandle> out (HandleManager.createObject
      (sampinfo_file,"table.elina_sampleinfo"));
    out->appendColumn("pointingID",pid);
    out->appendColumn("tot_samples",nsamp2);
    out->appendColumn("start_time_SCET",tstart);
    out->appendColumn("end_time_SCET",tend);
    }

  tstack_pop("period_loop");
  tstack_pop("multimod");
  tstack_report("multimod");
  return 0;
  }
