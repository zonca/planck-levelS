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
 *  Copyright (C) 2007-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "arr.h"
#include "paramfile.h"
#include "focalplane_db.h"
#include "iohandle_current.h"
#include "iohandle_helpers.h"
#include "fitshandle.h"
#include "safe_ptr.h"
#include "lsconstants.h"
#include "announce.h"
#include "string_utils.h"

using namespace std;

namespace {
//! converts from flux in Jy to antenna K * sr at \a freq.
double antennaFact (double freq)
  {
  return speedOfLight*speedOfLight/(2*freq*freq*kBoltzmann)*Jansky2SI;
  }

#include "detector_responses.cc"

class PntsrcEmitter
  {
  private:
    arr<string> catname, fluxcol, polangcol, polpercentcol;
    arr<double> catfreq;
    string cat_out, thetacol, phicol, namecol;
    bool polarised;

    void get_input(tsize m, arr<double> &fi, arr<double> &fq, arr<double> &fu,
      arr<double> &theta, arr<double> &phi, arr<string> &name)
      {
      safe_ptr<iohandle> inp (HandleManager.openObject(catname[m]));
      inp->readEntireColumn(thetacol,theta);
      inp->readEntireColumn(phicol,phi);
      inp->readEntireColumn(namecol,name);
      inp->readEntireColumn(fluxcol[m],fi);
      planck_assert(multiequal(theta.size(),phi.size(),name.size(),fi.size()),
        "inconsistent array sizes");
      if (polarised)
        {
        inp->readEntireColumn(polangcol[m],fq);
        inp->readEntireColumn(polpercentcol[m],fu);
        planck_assert(multiequal(fi.size(),fq.size(),fu.size()),
          "inconsistent array sizes");
        for (tsize i=0; i<fi.size(); ++i)
          {
          double polangle = fq[i]*degr2rad;
          double polfrac = fu[i]*0.01;
          fq[i] = fi[i]*polfrac*cos(2*polangle);
          fu[i] = fi[i]*polfrac*sin(2*polangle);
          }
        }
      }

  public:
    PntsrcEmitter (paramfile &params)
      {
      polarised = params.find<bool> ("polarised",false);
      int ncat = params.find<int> ("ncatalogs");
      planck_assert(ncat>=1,"PntsrcEmitter: need at least one catalog");
      cat_out = params.find<string> ("output_catalog");
      thetacol = params.find<string> ("colname_theta");
      phicol = params.find<string> ("colname_phi");
      namecol = params.find<string> ("colname_name");
      catname.alloc(ncat); fluxcol.alloc(ncat); catfreq.alloc(ncat);
      if (polarised) { polangcol.alloc(ncat); polpercentcol.alloc(ncat); }
      for (int m=0; m<ncat; ++m)
        {
        string suffix=dataToString(m+1);
        catname[m] = params.find<string>("catalog_"+suffix);
        fluxcol[m] = params.find<string>("colname_flux_"+suffix);
        catfreq[m] = params.find<double>("freq_"+suffix);
        if (polarised)
          {
          polangcol[m] = params.find<string>("colname_polang_"+suffix);
          polpercentcol[m] = params.find<string>("colname_polpercent_"+suffix);
          }
        }
      }

    virtual ~PntsrcEmitter() {}

    virtual void calc_factors (const Detector_Response &resp,
      arr<double> &factor) const
      {
      arr<double> freq, wgt;
      resp.get_weights (freq,wgt);

      tsize ncat = catfreq.size();
      factor.allocAndFill(ncat,0.);
      tsize il=0;
      for (tsize m=0; m<wgt.size(); ++m)
        {
        double fmid = 0.5*(freq[m]+freq[m+1]);
        if (ncat>1)
          {
          while ((il<(ncat-2)) && (fmid>catfreq[il+1]))
            ++il;
          double frac = (fmid-catfreq[il])/(catfreq[il+1]-catfreq[il]);
          double dfluxl = (1-frac)*wgt[m];
          double dfluxr = frac*wgt[m];
          factor[il]   += dfluxl*antennaFact(fmid);
          factor[il+1] += dfluxr*antennaFact(fmid);
          }
        else
          {
          factor[0] += wgt[m]*antennaFact(fmid);
          }
        }
      }
    virtual void calc_flux (const Detector_Response &resp)
      {
      arr<double> factor, fI, fQ, fU, theta, phi;
      arr<string> name;
      calc_factors (resp, factor);
      safe_ptr<iohandle> out (HandleManager.createObject(cat_out,
        polarised ? "catalog.LS_pointsource_catalog_pol_new"
                  : "catalog.LS_pointsource_catalog_new"));

      tsize ncat = catfreq.size();
      for (tsize m=0; m<ncat; ++m)
        {
        arr<double> fi_tmp, fq_tmp, fu_tmp, theta_tmp, phi_tmp;
        arr<string> name_tmp;
        get_input(m,fi_tmp,fq_tmp,fu_tmp,theta_tmp,phi_tmp,name_tmp);
        if (m==0)
          {
          fI.allocAndFill (fi_tmp.size(),0);
          theta=theta_tmp;
          phi=phi_tmp;
          name=name_tmp;
          if (polarised)
            {
            fQ.allocAndFill (fq_tmp.size(),0);
            fU.allocAndFill (fu_tmp.size(),0);
            }
          }
        else
          {
          planck_assert (fI.size()==fi_tmp.size(),
            "size mismatch in flux columns");
          planck_assert (theta_tmp.contentsEqual(theta),"theta mismatch");
          planck_assert (phi_tmp.contentsEqual(phi),"phi mismatch");
          planck_assert (name_tmp.contentsEqual(name),"source name mismatch");
          }
        for (tsize i=0; i<fI.size(); ++i)
          {
          fI[i]+=factor[m]*fi_tmp[i];
          if (polarised)
            {
            fQ[i]+=factor[m]*fq_tmp[i];
            fU[i]+=factor[m]*fu_tmp[i];
            }
          }
        }

      out->appendColumn ("theta_ecl",theta);
      out->appendColumn ("phi_ecl",phi);
      out->appendColumn ("source_name",name);
      out->appendColumn ("flux",fI);
      if (polarised)
        {
        arr<double> angle(fI.size()), percent(fI.size());
        for (tsize i=0; i<fI.size(); ++i)
          {
          angle[i] = rad2degr*0.5*atan2(fU[i],fQ[i]);
          percent[i] = 100.*sqrt(fQ[i]*fQ[i] + fU[i]*fU[i])/fI[i];
          }
        out->appendColumn ("polangle",angle);
        out->appendColumn ("polpercent",percent);
        }
      }
  };

} // unnamed namespace

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("pntsrcmixer", argc, argv);
  iohandle_current::Manager mng (argc, argv);
  paramfile params (mng.getParams());

  safe_ptr<Detector_Response> resp (make_response(params));
  PntsrcEmitter emitter (params);
  emitter.calc_flux (*resp);
PLANCK_DIAGNOSIS_END
  }
