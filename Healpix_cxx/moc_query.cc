/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file moc_query.cc
 *  Copyright (C) 2014 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "moc_query.h"
#include "geom_utils.h"
#include "lsconstants.h"

using namespace std;

namespace {

template<typename I> class querulator
  {
  private:
    int order, omax;
    bool inclusive;
    vector<MocQueryComponent> comp;
    vector<T_Healpix_Base<I> > base;
    arr<double> cr;
    arr2<double> crmin;
    arr2<double> crmax;

    vector<pair<I,int> > stk; // stack for pixel numbers and their orders
    I pix;
    int o;
    int stacktop; // a place to save a stack position
    vec3 pv;

    void check_pixel (int zone, Moc<I> &pixset)
      {
      if (zone==0) return;
      if (o<order)
        {
        if (zone>=3)
          pixset.appendPixel(o,pix); // output all subpixels
        else // (zone>=1)
          for (int i=0; i<4; ++i)
            stk.push_back(make_pair(4*pix+3-i,o+1)); // add children
        }
      else if (o>order) // this implies that inclusive==true
        {
        if (zone>=2) // pixel center in shape
          {
          pixset.appendPixel(order,pix>>(2*(o-order))); // output parent pixel
          stk.resize(stacktop); // unwind the stack
          }
        else // (zone>=1): pixel center in safety range
          {
          if (o<omax) // check sublevels
            for (int i=0; i<4; ++i) // add children in reverse order
              stk.push_back(make_pair(4*pix+3-i,o+1));
          else // at resolution limit
            {
            pixset.appendPixel(order,pix>>(2*(o-order))); // output parent pixel
            stk.resize(stacktop); // unwind the stack
            }
          }
        }
      else // o==order
        {
        if (zone>=2)
          pixset.appendPixel(order,pix);
        else if (inclusive) // and (zone>=1)
          {
          if (order<omax) // check sublevels
            {
            stacktop=stk.size(); // remember current stack position
            for (int i=0; i<4; ++i) // add children in reverse order
              stk.push_back(make_pair(4*pix+3-i,o+1));
            }
          else // at resolution limit
            pixset.appendPixel(order,pix); // output the pixel
          }
        }
      }

    void correctLoc(int &loc)
      {
      int myloc=loc--;
      planck_assert((myloc>=0)&&(myloc<int(comp.size())),"inconsistency");
      switch (comp[myloc].op)
        {
        case AND: case OR: case XOR:
          correctLoc(loc);
        case NOT:
          correctLoc(loc);
        case NONE:;
        }
      }
    int getZone (int &loc, int zmin, int zmax)
      {
      if (zmin==zmax) { correctLoc(loc); return zmin; } // short-circuit
      int myloc=loc--;
      planck_assert((myloc>=0)&&(myloc<int(comp.size())),"inconsistency");
      switch (comp[myloc].op)
        {
        case AND:
          {
          int z1=getZone(loc,zmin,zmax);
          return getZone(loc,zmin,z1);
          }
        case OR:
          {
          int z1=getZone(loc,zmin,zmax);
          return getZone(loc,z1,zmax);
          }
        case XOR:
          {
          int z1=getZone(loc,0,3);
          int z2=getZone(loc,0,3);
          return max(zmin,min(zmax,max(min(z1,3-z2),min(3-z1,z2))));
          }
        case NOT:
          return 3-getZone(loc,3-zmax,3-zmin);
        case NONE:
          {
          int res=zmax;
          double crad=dotprod(pv,comp[myloc].center);
          if (crad<=crmax(o,myloc)) res=0;
          else if (crad<=cr[myloc]) res=1;
          else if (crad<=crmin(o,myloc)) res=2;
          return max(zmin,min(zmax,res));
          }
        }
      planck_fail("must not get here");
      }

  public:
    querulator (int order_, int omax_, bool inclusive_,
      const vector<MocQueryComponent> &comp_)
      : order(order_), omax(omax_), inclusive(inclusive_), comp(comp_),
        base(omax+1), cr(comp.size()), crmin(omax+1,comp.size()),
        crmax(omax+1,comp.size())
      {
      planck_assert(comp.size()>=1,"bad query component vector");
      planck_assert(order<=omax,"order>omax");
      if (!inclusive) planck_assert(order==omax,"inconsistency");
      planck_assert(omax<=T_Healpix_Base<I>::order_max,"omax too high");

      for (tsize i=0; i<comp.size(); ++i)
        if (comp[i].op==NONE) // it's a cap
          cr[i]=cos(comp[i].radius);
      for (o=0; o<=omax; ++o) // prepare data at the required orders
        {
        base[o].Set(o,NEST);
        double dr=base[o].max_pixrad(); // safety distance
        for (tsize i=0; i<comp.size(); ++i)
          if (comp[i].op==NONE) // it's a cap
            {
            crmax(o,i) = (comp[i].radius+dr>=pi) ? -1.01 : cos(comp[i].radius+dr);
            crmin(o,i) = (comp[i].radius-dr<=0.) ?  1.01 : cos(comp[i].radius-dr);
            }
        }
      }
    Moc<I> result()
      {
      Moc<I> pixset;
      stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
      for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
        stk.push_back(make_pair(I(11-i),0));

      stacktop=0; // a place to save a stack position

      while (!stk.empty()) // as long as there are pixels on the stack
        {
        // pop current pixel number and order from the stack
        pix=stk.back().first;
        o=stk.back().second;
        stk.pop_back();
        pv = base[o].pix2vec(pix);

        int loc=comp.size()-1;
        tsize zone=getZone(loc,0,3);
        check_pixel (zone, pixset);
        planck_assert(loc==-1,"stack not used up");
        }
      return pixset;
      }
  };

} // unnamed namespace

template<typename I> Moc<I> mocQuery (int order,
  const vector<MocQueryComponent> &comp)
  {
  querulator<I> quer(order,order,false,comp);
  return quer.result();
  }
template Moc<int> mocQuery (int order, const vector<MocQueryComponent> &comp);
template Moc<int64> mocQuery (int order, const vector<MocQueryComponent> &comp);

template<typename I> Moc<I> mocQueryInclusive (int order, int omax,
  const vector<MocQueryComponent> &comp)
  {
  querulator<I> quer(order,omax,true,comp);
  return quer.result();
  }
template Moc<int> mocQueryInclusive (int order, int omax,
  const vector<MocQueryComponent> &comp);
template Moc<int64> mocQueryInclusive (int order, int omax,
  const vector<MocQueryComponent> &comp);

bool isEar(tsize i, const vector<vec3> &vv, const vector<vec3> &normal, const vector<bool> &convex)
  {
  if (!convex[i]) return false;
  tsize nv=normal.size();
  tsize pred=(i+nv-1)%nv, succ=(i+1)%nv;
  vec3 ab=normal[i],bc=crossprod(vv[succ],vv[pred]).Norm(),ca=normal[pred];
  for (tsize j=(succ+1)%nv; j!=pred; j=(j+1)%nv)
    if ((dotprod(vv[j],ab)>0)
      &&(dotprod(vv[j],bc)>0)
      &&(dotprod(vv[j],ca)>0))
      return false;
  return true;
  }

// vertices are assumed to be normalized
void processVertices(vector<vec3> &vv, vector<vec3> &normal,
  vector<bool> &convex, vector<bool> &ear)
  {
  tsize nv=vv.size();
  do
    {
    nv=vv.size();
    for (tsize i=0; i<vv.size(); ++i)
      {
      tsize succ=(i+1)%vv.size();
      double v_ang = v_angle(vv[i],vv[succ]);
      if (v_ang<1e-8) // vertices are very close
        {
        vv.erase(vv.begin()+i);
        break;
        }
      if (v_ang>pi-1e-8) // vertices almost antipodal
        {
        planck_fail("degenerate polygon edge");
        break;
        }
      tsize pred=(i+vv.size()-1)%vv.size();
      vec3 npred=crossprod(vv[pred],vv[i]),
           nsucc=crossprod(vv[i],vv[succ]);
      double n_ang = v_angle(npred,nsucc);
      if (n_ang<1e-8) // vertices almost on a straight line
        {
        vv.erase(vv.begin()+i);
        break;
        }
      if (n_ang>pi-1e-8) // very thin polygon spike
        {
        vv.erase(vv.begin()+i);
        break;
        }
      }
    } while (nv!=vv.size());

  normal.resize(nv);
  convex.resize(nv);
  ear.resize(nv);
  for (tsize i=0; i<nv;++i)
    normal[i]=crossprod(vv[i],vv[(i+1)%nv]).Norm();
  for (tsize i=0; i<nv;++i)
    convex[i]=dotprod(normal[i],vv[(i+nv-1)%nv])>0;
  for (tsize i=0; i<nv;++i)
    ear[i]=isEar(i,vv,normal,convex);
  }

vector<MocQueryComponent> prepPolygon (const vector<vec3> &vertex)
  {
  planck_assert(vertex.size()>=3,"not enough vertices in polygon");
  vector<vec3> vv(vertex.size());
  for (tsize i=0; i<vertex.size(); ++i)
    vv[i]=vertex[i].Norm();
  vector<vec3> normal;
  vector<bool> convex, ear;

  processVertices(vv,normal,convex,ear);
  vector<MocQueryComponent> comp;

  tsize nconvex=0;
  for (tsize i=0; i<vv.size(); ++i)
    if (convex[i]) ++nconvex;
  if (nconvex==vv.size()) // fully convex polygon
    {
    for (tsize i=0; i<vv.size(); ++i)
      comp.push_back(MocQueryComponent(normal[i],halfpi));
    for (tsize i=1; i<vv.size(); ++i)
      comp.push_back(MocQueryComponent(AND));
    return comp;
    }
  if (nconvex==0) // complement of a fully convex polygon??
    {
    for (tsize i=0; i<vv.size(); ++i)
      comp.push_back(MocQueryComponent(-normal[i],halfpi));
    for (tsize i=1; i<vv.size(); ++i)
      comp.push_back(MocQueryComponent(AND));
    comp.push_back(MocQueryComponent(NOT));
    return comp;
    }

  int earcount=0;
  while (vv.size()>2) // try to clip ears
    {
    tsize nvlast=vv.size();
    for (tsize i=0; i<vv.size();++i)
      {
      tsize pred=(i+vv.size()-1)%vv.size(), succ=(i+1)%vv.size();

      if (ear[i])
        {
        comp.push_back(MocQueryComponent(normal[i],halfpi));
        comp.push_back(MocQueryComponent(crossprod(vv[succ],vv[pred]).Norm(),halfpi));
        comp.push_back(MocQueryComponent(AND));
        comp.push_back(MocQueryComponent(normal[pred],halfpi));
        comp.push_back(MocQueryComponent(AND));
        ++earcount;
        vv.erase(vv.begin()+i);
        processVertices(vv,normal,convex,ear);
        break;
        }
      }
    planck_assert(vv.size()<nvlast,"failed to clip an ear");
    }
  for (int i=0; i<earcount-1; ++i)
    comp.push_back(MocQueryComponent(OR));
  return comp;
  }
