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

class Detector_Response
  {
  public:
//  Detector_Response (paramfile &params);
    virtual ~Detector_Response() {}
//  sum of all weights must be 1
//  on exit: freq.size() == wgt.size()+1
    virtual void get_weights (arr<double> &freq, arr<double> &wgt) const = 0;
  };

class Tophat: public Detector_Response
  {
  private:
    double flo, fhi;

  public:
    Tophat (paramfile &params)
      {
      focalplane_db database(params);
      string det_id = params.find<string> ("detector_id");
      flo=database.getValue<double>(det_id,"nu_min");
      fhi=database.getValue<double>(det_id,"nu_max");
      }

    virtual void get_weights (arr<double> &freq, arr<double> &wgt) const
      {
      const int nsteps = 100;
      freq.alloc(nsteps+1);
      wgt.alloc(nsteps);
      double df = (fhi-flo)/nsteps;
      for (int m=0; m<nsteps; ++m)
        {
        freq[m] = flo+m*df;
        wgt[m] = 1./nsteps;
        }
      freq[nsteps] = fhi;
      }
  };

class Delta: public Detector_Response
  {
  private:
    double f;

  public:
    Delta (paramfile &params)
      {
      focalplane_db database(params);
      string det_id = params.find<string> ("detector_id");
      f=database.getValue<double>(det_id,"nu_cen");
      }

    virtual void get_weights (arr<double> &freq, arr<double> &wgt) const
      {
      freq.alloc(2);
      wgt.alloc(1);
      freq[0] = 0.99999*f;
      freq[1] = 1.00001*f;
      wgt[0] = 1;
      }
  };

class Gauss: public Detector_Response
  {
  private:
    double mean, sigma;

  public:
    Gauss (paramfile &params)
      {
      focalplane_db database(params);
      string det_id = params.find<string> ("detector_id");
      mean=database.getValue<double>(det_id,"nu_cen");
      sigma=database.getValue<double>(det_id,"nu_max")
           -database.getValue<double>(det_id,"nu_min");
      }

    virtual void get_weights (arr<double> &freq, arr<double> &wgt) const
      {
      const int nsteps = 100;
      freq.alloc(nsteps+1);
      wgt.alloc(nsteps);
      double flo = max(mean-5*sigma,0.001*mean); // do not go across f=0
      double fhi = mean+5*sigma;
      double df = (fhi-flo)/nsteps;
      double normfact = (fhi-flo)/(nsteps*sqrt(pi)*sigma);
      for (int m=0; m<nsteps; ++m)
        {
        freq[m] = flo+m*df;
        double fmid = freq[m]+0.5*df;
        wgt[m] = normfact*exp(-pow((fmid-mean)/sigma,2));
        }
      freq[nsteps] = fhi;
      }
  };

class Tabular: public Detector_Response
  {
  private:
    arr<double> frequency, response;

  public:
    Tabular (paramfile &params)
      {
      safe_ptr<iohandle> inp(HandleManager.openObject
        (params.find<string>("response_table"), "table.LS_detector_response"));
      inp->readEntireColumn("frequency",frequency);
      inp->readEntireColumn("response",response);
      planck_assert(frequency.size()==response.size()+1,
        "inconsistent column sizes in detector response object");
      planck_assert(frequency.size()>=2,
        "tabular frequency response: need at least 2 frequencies");
      bool do_norm = params.find<bool>("normalise_bandpass",false);
      if (do_norm)
        {
        cout << "normalising the bandpass" << endl;
        double norm = 0;
        for (tsize m=0; m<response.size(); ++m)
          norm+=response[m];
        norm = 1./norm;
        for (tsize m=0; m<response.size(); ++m)
          response[m]*=norm;
        }
      }

    virtual void get_weights (arr<double> &freq, arr<double> &wgt) const
      {
      freq = frequency;
      wgt = response;
      }
  };

Detector_Response *make_response (paramfile &params)
  {
  string rtype = params.find<string> ("response_type");
  if (rtype=="DELTA")
    return new Delta(params);
  else if (rtype=="TOPHAT")
    return new Tophat(params);
  else if (rtype=="TABULAR")
    return new Tabular(params);
  else if (rtype=="GAUSS")
// FIXME
    planck_fail ("response type GAUSS has been disabled for now.");
//    return new Gauss(params);
  else
    planck_fail("Bad response_type '"+rtype+"'.");
  }
