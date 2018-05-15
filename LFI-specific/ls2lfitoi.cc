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
 * This file implements a conversion for Level S signal TOI
 * to LFI data type toi.science.LFI_Data of LFI DMC.
 * This takes already generated T_sky signal (with or without
 * instrumental noise), creates T_load time stream (either with
 * provided variations (analytical formula) or stable 4 K load.
 *
 * For the time being it does not include sampleOBT for
 * simplicity. In the future will be inserted properly.
 *
 * 12-11-06 : Original version delivered also to Level S
 *            CVS repository in ESTEC for End-To-End tests
 *
 * 01-04-07 : Update after discussion with Paddy Leahy
 *            for proper generation of R accounting for
 *            total power 1/f noise (more realism in the
 *            estimation of R).
 *
 * 12-13-07 : Update version for inclusion of:
 *            - thermal effects (4 and 20 Kon both Sky and Load)
 *            - frequency spikes (due to H/K interference)
 *            - time stamp as produced by LevelS
 *
 * 03-26-08 : Update to include
 *            - thermometry for the telescope. Assuming to get
 *              data as base temp ref + flucts. Convert physical
 *              temperature into antenna temperature and adds them
 *              to singal from Sky horn
 *
 * Copyright (C) 2006-2007-2008, LFI-DPC
 * Author: Davide Maino
 */

#include "arr.h"
#include "lsconstants.h"
#include "math.h"
#include "focalplane_db.h"
#include "announce.h"
#include "paramfile.h"
#include "iohandle_current.h"

using namespace std;

namespace {

inline double T_ant (double t_thermo, double hnydk)
  { return hnydk/(exp(hnydk/t_thermo)-1); }

void get_params (const string &name, paramfile &par_file,
  string &outfile, bool &write)
  {
  outfile = par_file.find<string> (name+string("_file"),"");
  write = (outfile!="");
  }

void checkDoubleData (const iohandle &inp, double &value, const string &name)
  {
  if (inp.keyPresent(name))
    {
    double dtmp=inp.getKey<double>(name);
    if (value<0)
      value=dtmp;
    else
      planck_assert (approx(dtmp,value), "inconsistent values for '"+name+"'");
    }
  else
    cout << "Warning: key '" << name << "' missing in input file" << endl;
  }

void assembleData (const iohandle &inp, int64 &first_sample,
  double &start_time_mission, double &start_time, double &f_samp, uint64 &nsamp)
  {
  int incol = inp.columnNumber("signal");
  uint64 utmp = inp.columnLength(incol);
  if (nsamp==0)
    nsamp=utmp;
  else
    planck_assert(utmp==nsamp,"trying to coadd TOD of different lengths");

  if (inp.keyPresent("first_sample"))
    {
    int64 itmp=inp.getKey<int64>("first_sample");
    if (first_sample<0)
      first_sample=itmp;
    else
      planck_assert (itmp==first_sample,
        "inconsistent values for 'first_sample'");
    }
  else
    cout << "Warning: key 'first_sample' missing in input file" << endl;
  checkDoubleData (inp, start_time_mission, "start_time_mission");
  checkDoubleData (inp, start_time, "start_time");
  checkDoubleData (inp, f_samp, "f_samp");
  }

} // unnamed namespace

int main (int argc, const char **argv)
  {
  module_startup ("ls2lfitoi", argc, argv);
  iohandle_current::Manager mng(argc, argv);
  paramfile params (mng.getParams());

  // This section checks if thermal, spikes and time files
  // has to be read in and considered in the creation of LFI TOI
  string therm_sky_file, therm_load_file, fspike_file, time_file, tele_file;
  bool read_thermsky, read_thermload, read_fspike, read_time, read_tele;
  get_params("therm_sky",params,therm_sky_file,read_thermsky);
  get_params("therm_load",params,therm_load_file,read_thermload);
  get_params("fspike",params,fspike_file,read_fspike);
  get_params("time",params,time_file,read_time);
  get_params("telescope",params,tele_file,read_tele);

  bool do_time = read_time,
    do_therm_sky = read_thermsky,
    do_therm_load = read_thermload,
    do_fspike = read_fspike,
    do_tele = read_tele;

  cout << "do_time " << do_time << " do_tsky "<< do_therm_sky << endl;
  cout << "do_tload " << do_therm_load << " do_spike "<<do_fspike << endl;
  cout << "do_tele " << do_tele << endl;

  double sqrt2 = sqrt(2.);
  double tref_temp = params.find<double>("tref_temp",4.5);
  double t_telescope = params.find<double>("t_telescope",42.);
  double emissivity = params.find<double>("emissivity",1.);
  bool has_monopole = params.find<bool>("has_monopole");
  string det_id = params.find<string>("detector_id");
  string tsky = params.find<string>("tsky");
  string tsky_white_noise = params.find<string>("tsky_noise");
  string tref_white_noise = params.find<string>("tref_noise");
  string oof_noise = params.find<string>("oof_noise");
  string oof_noise_tp = params.find<string>("oof_noise_tp");
  string file_toi = params.find<string>("file_toi");
  bool verbosity = params.find<bool>("verbosity",false);

  /* Read data from FPDB to get the actual band, noise and frequency */
  focalplane_db fpdb(params);
  double freq = fpdb.getValue<double>(det_id,"nu_cen");
  double net_rj = fpdb.getValue<double>(det_id,"net_rj");
  double band = fpdb.getValue<double>(det_id,"nu_max")
              - fpdb.getValue<double>(det_id,"nu_min");

  double hnydk = hPlanck*freq/kBoltzmann;
  double ta=0;
  if (!has_monopole)
    ta = T_ant(tcmb,hnydk);
  double Tnoise = net_rj * sqrt(band)/sqrt2 - T_ant(tcmb,hnydk);
  //WARNING: R0 is computed assuming a physical fluctuation around
  // 4.5 K as in code Luca Terenzi's code. This is the reason why
  // 4.5 K is the default value for tref_temp. If you have
  // fluctuations with other mean value you have to insert it
  // in tref_temp so that at the final stage it is properly canceled
  // out. Otherwise you are left with an offset.
  double ta_ref = T_ant(tref_temp,hnydk);
  double ta_tele_mean = T_ant(t_telescope, hnydk);
  ta_tele_mean *= emissivity;
  double R0 = (ta_ref + Tnoise)/(T_ant(tcmb,hnydk) + Tnoise + ta_tele_mean);
  if (verbosity)
    cout << "R0 = " << 1./R0 << endl;

  safe_ptr<iohandle> inp_sky (HandleManager.openObject(tsky,"toi.LS_toi"));

  /* from this file get the time information. Then check
     that this is consistent with the same keywords in all
     the other files. For thermal and spikes which are supposed
     to be a single file (for the time being) we simply move
     along the file till finding the required sample */
  double start_time_mission=-1, start_time=-1, f_samp=-1;
  int64 first_sample = -1;
  uint64 nsamp=0;
  assembleData (*inp_sky, first_sample, start_time_mission, start_time,
    f_samp, nsamp);

  safe_ptr<iohandle> wn_sky (HandleManager.openObject(tsky_white_noise,
    "toi.LS_toi"));
  safe_ptr<iohandle> wn_load (HandleManager.openObject(tref_white_noise,
    "toi.LS_toi"));
  safe_ptr<iohandle> oof (HandleManager.openObject(oof_noise,"toi.LS_toi"));
  safe_ptr<iohandle> ooftp (HandleManager.openObject(oof_noise_tp,
    "toi.LS_toi"));

  // This part is for the new inclusion of systematic effects
  // Defines variables and then, if required, link them to the
  // proper handle
  safe_ptr<iohandle> inp_therm_sky, inp_therm_load, inp_fspike, inp_time,
    inp_tele;
  int thskycol=-1,thloadcol=-1,fspikecol=-1,ticol=-1,telecol=-1;

  if (do_fspike) {
    inp_fspike=HandleManager.openObject(fspike_file,"toi.LS_toi");
    fspikecol = inp_fspike->columnNumber("signal");
  }
  if (do_therm_sky) {
    inp_therm_sky=HandleManager.openObject(therm_sky_file,"toi.LS_toi");
    thskycol = inp_therm_sky->columnNumber("signal");
  }
  if (do_therm_load) {
    inp_therm_load=HandleManager.openObject(therm_load_file,"toi.LS_toi");
    thloadcol = inp_therm_load->columnNumber("signal");
  }
  if (do_time) {
    inp_time=HandleManager.openObject(time_file,"toi.LS_timestamps");
    ticol = inp_time->columnNumber("timestamp");
  }
  if (do_tele) {
    inp_tele=HandleManager.openObject(tele_file,"toi.LS_toi");
    telecol = inp_tele->columnNumber("signal");
  }

  int skycol = inp_sky->columnNumber("signal");
  int wnskycol = wn_sky->columnNumber("signal");
  int wnloadcol = wn_load->columnNumber("signal");
  int oofcol = oof->columnNumber("signal");
  int ooftpcol = ooftp->columnNumber("signal");

  uint64 nsamples_sky = inp_sky->columnLength(skycol);
  uint64 nsamples_ref = wn_load->columnLength(wnloadcol);
  uint64 nsamples = min(nsamples_sky,nsamples_ref);

  safe_ptr<iohandle> out (HandleManager.createObject(file_toi,
    "toi.science.LFI_Data"));

  int scetcol = out->columnNumber("sampleSCET");
  int oskycol = out->columnNumber("sky_volt");
  int orefcol = out->columnNumber("load_volt");

  arr<double> dval, dval1;
  arr<float> val, val1, val2, val3;
  uint64 offset =0;
  const uint64 chunksize=1024*256;
  while(offset<nsamples)
    {
    uint64 psamp=min(chunksize,nsamples-offset);
    val.alloc(psamp);
    val1.alloc(psamp);
    val2.alloc(psamp);
    val3.alloc(psamp);
    dval.alloc(psamp);

    // Start with Sky signal creation
    inp_sky->readColumn(skycol,val,offset);
    if (verbosity) cout << " signal[0] = " << val[0] << endl;
    for (uint64 m=0;m<psamp;++m) {
      if (do_tele) {
         dval[m] = double(val[m]) + ta + Tnoise;
      } else {
         dval[m] = double(val[m]) + ta + Tnoise + ta_tele_mean;
      }
    }
    wn_sky->readColumn(wnskycol,val,offset);
    if (verbosity) cout << " wn_sky[0] = " << val[0] << endl;
    for (uint64 m=0;m<psamp;++m)
      dval[m] += double(val[m])/sqrt2;
    oof->readColumn(oofcol,val,offset);
    if (verbosity) cout << " oof diff[0] = " << val[0] << endl;
    for (uint64 m=0;m<psamp;++m)
      dval[m] += double(val[m])/2.;
    ooftp->readColumn(ooftpcol,val1,offset);
    if (verbosity) cout << " oof tp[0] = " << val1[0] << endl;
    for (uint64 m=0;m<psamp;++m)
      dval[m] += double(val1[m]);

    //This part is for the inclusion of systematic effects
    if (do_therm_sky) {
      inp_therm_sky->readColumn(thskycol,val3,offset+first_sample);
      for (uint64 m=0;m<psamp;++m)
        dval[m] += double(val3[m]);
    }
    if (do_fspike) {
      inp_fspike->readColumn(fspikecol,val3,offset+first_sample);
      for (uint64 m=0;m<psamp;++m)
        dval[m] += double(val3[m]);
    }
    if (do_tele) {
      inp_tele->readColumn(telecol,val3,offset+first_sample);
      for (uint64 m=0;m<psamp;++m) {
        double ta_tele = T_ant(t_telescope + val3[m],hnydk);
      // dval[m] += emissivity*double(val3[m]);
        dval[m] += emissivity * ta_tele;
      }
    }
    if (verbosity) cout << " sky[0] = " << dval[0] << endl;
    out->appendColumn(oskycol,dval);

    // This includes SCET into the object
    // Definition compliant with LFI Level 1
    if (do_time) {
      dval1.alloc(psamp);
      inp_time->readColumn(ticol,dval1,offset);
      cout << "time[1] = " << dval1[1] << endl;
      out->appendColumn(scetcol,dval1);
    }

    // Now create the Load signal
    wn_load->readColumn(wnloadcol,val2,offset);
    if (verbosity) cout << "wn load[0] = " << val2[0] << endl;
    for (uint64 m=0;m<psamp;++m)
      dval[m] = ta_ref + Tnoise + R0*double(val2[m])/sqrt2 +
                R0*(double(val1[m]) - double(val[m]/2.));
    if (verbosity) cout << " load[0] = " << dval[0] << endl;

    // This part adds the thermal effect on the Load signal
    // and remove the mean value of the fluctuation which
    // is inserter in dval variable for the constant tref
    // case.
    if (do_therm_load) {
      inp_therm_load->readColumn(thloadcol,val3,offset+first_sample);
      for (uint64 m=0;m<psamp;++m)
        dval[m] = dval[m]+double(val3[m])-ta_ref;
    }
    out->appendColumn(orefcol,dval);
    offset+=chunksize;
    }
  }
