/*
 *  This file is part of the Planck simulation package
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
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  This code is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003-2014 Max-Planck-Society
 *  \author Martin Reinecke
 *  \author Reinhard Hell
 */

#include "PointSourceConvolver.h"
#include "iohandle_current.h"
#include "geom_utils.h"
#include "Planet.h"
#include "sat_info.h"
#include "ephemerides.h"
#include "io_utils.h"
#include "string_utils.h"

using namespace std;

namespace psc {

PtSrcMapServer::PtSrcMapServer (paramfile &params, const focalplane_db &fpdb,
  const string &det_id, const Sat_Info &satinfo)
  {
  double freqGHz = fpdb.getValue<double>(det_id,"nu_cen")/1e9;
  bool polarisation = params.find<bool>("pntsrc_polarisation",false);
  double brm = degr2rad*params.find<double>("beam_radius_max");
  varPntsrcFactor = params.find<double>("variable_pntsrc_factor",1.);
  int nside=nearest<int>(2.5*sqrt(fourpi/(brm*brm*pi*12)));
  nside=max(1,min(nside,128));
  cout << "determined srcnside is " << nside << endl;
  srcMap.SetNside(nside,RING);
  readFixedSources(params, polarisation);
  readPlanets(params, satinfo, freqGHz);
  }

void PtSrcMapServer::readFixedSources (paramfile &params, bool polarisation)
  {
  string infile = params.find<string>("pntsrc_file","");
  if (infile=="") return;

  arr<pointing> pos;
  arr<double> flux, polang, polfrac;
  arr<string> name;
  readPointSourcesNew (infile, polarisation, pos, flux, polang, polfrac, name,
    1, 0);

  for (tsize j=0; j<flux.size(); ++j)
    {
    pos[j].normalize();
    SmallPtSrc sps (pos[j],flux[j],polang[j],polfrac[j]);
    fixedSources.push_back(FixedPtSrc(sps,name[j]));
    int ipix=srcMap.ang2pix(pos[j]);
    srcMap[ipix].push_back(QuickPtSrc(fixedSources[j].getData(0),j));
    }
  }

namespace {

void readPlanet (Planet &planet, vector<VariablePtSrc> &varsrc)
  {
  const arr<SmallPtSrc> &origsrc (planet.getSources());
  int npp=origsrc.size();
  varsrc.push_back(VariablePtSrc(npp,planet.planet_name));
  VariablePtSrc &src(varsrc[varsrc.size()-1]);
  for (int m=0; m<npp; ++m)
    src.setData(m,origsrc[m]);
  }

} // unnamed namespace

void PtSrcMapServer::readPlanets (paramfile &params, const Sat_Info &satinfo,
  double freqGHz)
  {
  string infile = params.find<string>("planet_file","");
  if (infile=="") return;

  const arr<double> tstart(satinfo.startTimes()),
                    tend(satinfo.endTimes());

  safe_ptr<ephemerides> eph(getEphemerides(infile));

  string marstemp = params.find<string>("mars_temperatures","");
  if (marstemp=="")
    {
    RockyPlanet mars(MarsData, *eph, tstart, tend, freqGHz);
    readPlanet (mars, variableSources);
    }
  else
    {
    MarsPlanet mars(MarsData, *eph, tstart, tend, marstemp);
    readPlanet (mars, variableSources);
    }

  RockyPlanet pluto(PlutoData, *eph, tstart, tend, freqGHz);
  GasPlanet jupiter(JupiterData, *eph, tstart, tend, freqGHz),
            saturn(SaturnData, *eph, tstart, tend, freqGHz),
            uranus(UranusData, *eph, tstart, tend, freqGHz),
            neptune(NeptuneData, *eph, tstart, tend, freqGHz);

  readPlanet (jupiter, variableSources);
  readPlanet (saturn, variableSources);
  readPlanet (uranus, variableSources);
  readPlanet (neptune, variableSources);
  readPlanet (pluto, variableSources);
  }

PointSourceConvolver::PointSourceConvolver (paramfile &params,
  const focalplane_db &fpdb, const string &det_id, const Sat_Info &satinfo)
  : mapServer (params, fpdb, det_id, satinfo)
  {
  // Read in parameters:
  bool polarisation = params.find<bool>("pntsrc_polarisation",false);
  statFile = params.find<string>("psrchits_file","");

  double maxradius = params.find<double>("beam_radius_max");
  planck_assert (maxradius<=180.0,
    "Error: The maximum beam radius must not be larger than 180 degrees.");

  int beamType = params.find<int>("beam_file_type");
  switch (beamType)
    {
    case 2: // analytic Gauss beam
      {
      double epsilon = fpdb.getValue<double>(det_id,"epsilon");
      double fwhm = degr2rad*fpdb.getValue<double>(det_id,"beamfwhm");
      double polangle = degr2rad*fpdb.getValue<double>(det_id,"psi_pol");
      beamCharac = new GaussBeam(fwhm,polarisation,polangle,epsilon);
      break;
      }
    case 3: // elliptic Gauss beam
      {
      double epsilon = fpdb.getValue<double>(det_id,"epsilon");
      double fwhm = degr2rad*fpdb.getValue<double>(det_id,"beamfwhm");
      double ellipticity = fpdb.getValue<double>(det_id,"ellipticity");
      double fwhm1 = fwhm*sqrt(ellipticity);
      double fwhm2 = fwhm/sqrt(ellipticity);
      double rotangle = degr2rad*fpdb.getValue<double>(det_id,"psi_ell");
      double polangle = degr2rad*fpdb.getValue<double>(det_id,"psi_pol");
      beamCharac = new EllipticGaussBeam(fwhm1, fwhm2, rotangle, polarisation,
        polangle,epsilon);
      break;
      }
    case 4: // beam on ECP grid
      {
      string beamFile = params.find<string>("beam_file");
      beamCharac = new GridBeam (beamFile, polarisation);
      break;
      }
    case 5: // beam on Cartesian grid
      {
      string beamFile = params.find<string>("beam_file");
      beamCharac = new CartBeam (beamFile, polarisation);
      break;
      }
    case 0: case 1:
      planck_fail("obsolete beam type " + dataToString(beamType));
    default:
      planck_fail("wrong beam type " + dataToString(beamType));
    }

  double freqGHz = fpdb.getValue<double>(det_id,"nu_cen")/1e9;
  cout << "FocalPlaneDB: frequency = " << freqGHz << endl;

  beamRadius = maxradius*degr2rad;
  if (beamCharac->max_theta()<beamRadius)
    {
    cout << "reducing maximum beam radius to "
         << beamCharac->max_theta()*rad2degr << " degrees" << endl;
    beamRadius = beamCharac->max_theta();
    }

  if (statFile!="")
    {
    int nsrc = mapServer.nSources();
    ppsOfHits.alloc(nsrc);
    nrOfObs.alloc(nsrc);
    nrOfObs.fill(0);
    ppOfObsMax.alloc(nsrc);
    ppOfObsMax.fill(-1);
    maxSrc.alloc(nsrc);
    seen.alloc(nsrc);
    }
  }

void PointSourceConvolver::Get_Intensities(int pp, const arr<pointing> &detpt,
  const arr<vec3> &vdetpt, const arr<double> &heading, arr<double> &intensity)
  {
  intensity.alloc(detpt.size());

  const PtSrcMapType &srcMap (mapServer.provideMap(pp));
  bool asymmetric = beamCharac->asymmetric();
  bool do_stat = (statFile!="");
  if (do_stat) seen.fill(false);

  double beamRadius_incl = beamRadius+srcMap.max_pixrad();
  double cosBeamRadius = cos(beamRadius);
#ifdef _OPENMP
  int parflag = do_stat ? 0 : 1;
#pragma omp parallel if (parflag)
#endif
{
  vector<int> pixelList;

  int ppp, sz=detpt.size();
#pragma omp for schedule(static)
  for (ppp=0; ppp<sz; ++ppp)
    {
    double frac = ppp/(sz-1.);
    intensity[ppp]=0.0;
    vec3 detpntVec=vdetpt[ppp];

    // Find point sources near the scanning directions
    srcMap.query_disc(detpt[ppp], beamRadius_incl, pixelList);

    // Add intensities of neighbouring point sources
    double alphaSat=heading[ppp];
    for (unsigned int k=0; k<pixelList.size(); ++k)
      {
      int pixel=pixelList[k];
      for (unsigned int n=0; n<srcMap[pixel].size(); ++n)
        {
        // Calculate position angles of point source in beam array (focal plane)
        const QuickPtSrc &qps(srcMap[pixel][n]);
        vec3 qpos = qps.locV(frac);
        double cosAngleDistPnt2Src=min(dotprod(qpos,detpntVec),1.);
        if (cosAngleDistPnt2Src>=cosBeamRadius)
          {
          double angleDistPnt2Src=acos(cosAngleDistPnt2Src);
          double alpha=0;
          if (asymmetric)
            alpha=orientation(detpntVec,qpos-detpntVec)-alphaSat;

          // Count the number of observations
          if (do_stat)
            {
            ++nrOfObs[qps.id];
            if (!seen[qps.id])
              {
              seen[qps.id]=true;
              vector<interval> &hits(ppsOfHits[qps.id]);
              if ((hits.size()==0) || (!hits[hits.size()-1].tryToAdd(pp)))
                hits.push_back(interval(pp));
              if (qps.flux>maxSrc[qps.id].flux)
                {
                maxSrc[qps.id]=qps;
                ppOfObsMax[qps.id]=pp;
                }
              }
            }

          // Determine the intensity at the detector (in antenna K)
          vec3 eastdir_beam (-detpntVec.y, detpntVec.x, 0);
          double delta_psi=halfpi+orientation(qpos,eastdir_beam);
          double scalar, pol;
          beamCharac->beamValue(angleDistPnt2Src,alpha,
            delta_psi+alphaSat-qps.polangle,scalar,pol);
          intensity[ppp] += qps.flux * (scalar + pol*qps.polfrac);
          }
        }
      }
    }
} // end of parallel region

  }

void PointSourceConvolver::Write_Hitsfile() const
  {
  if (statFile=="") return;

  int nsrc=mapServer.nSources();
  arr<string> pntsrcNames(nsrc), ppstrOfHits(nsrc);
  arr<double> thetaAtMax(nsrc), phiAtMax(nsrc), fluxAtMax(nsrc);
  arr<int> pptNrOfHits(nsrc), cppNrOfHits(nsrc);

  for (int m=0; m<nsrc; ++m)
    {
    const PtSrc &psc=mapServer.getPtSrc(m);
    pntsrcNames[m]=psc.name();
    pointing ptg=maxSrc[m].locP(0.5);
    thetaAtMax[m]=ptg.theta*rad2degr;
    phiAtMax[m]=ptg.phi*rad2degr;
    fluxAtMax[m]=maxSrc[m].flux;
    cppNrOfHits[m]=0;

    string &ppstr (ppstrOfHits[m]);
    pptNrOfHits[m]=0;
    const vector<interval> &hits(ppsOfHits[m]);
    for (unsigned int i=0; i<hits.size(); ++i)
      {
      if (ppstr!="") ppstr +=" ";
      ppstr+=dataToString(hits[i].lo);
      if (hits[i].lo!=hits[i].hi)
        ppstr+="-"+dataToString(hits[i].hi);
      pptNrOfHits[m]+= hits[i].hi-hits[i].lo+1;
      ++cppNrOfHits[m];
      }
    }

  safe_ptr<iohandle> out(HandleManager.createObject (statFile,
    "module.psc.LS_pointsource_hits"));
  out->appendColumn ("source_name", pntsrcNames);
  out->appendColumn ("max_period", ppOfObsMax);
//  out->appendColumn ("max_intensity", fluxAtMax);
  out->appendColumn ("max_antenna_temp", fluxAtMax);
  out->appendColumn ("max_theta", thetaAtMax);
  out->appendColumn ("max_phi", phiAtMax);
  out->appendColumn ("number_of_hits", nrOfObs);
  out->appendColumn ("periods_of_hits", pptNrOfHits);
  out->appendColumn ("contiguous_periods", cppNrOfHits);
  out->appendColumn ("Times_of_hits", ppstrOfHits);
  }

} // end of namespace psc
