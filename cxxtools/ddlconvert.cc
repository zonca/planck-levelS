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
 *  Copyright (C) 2006-2013 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include "announce.h"
#include "error_handling.h"

using namespace std;

namespace {

map<string,string> ancestor;

string f_typeconvert(const string &tp)
  {
  if (tp=="BOOL") return "PLANCK_BOOL";
  if (tp=="BYTE") return "PLANCK_INT8";
  if (tp=="INT8") return "PLANCK_INT8";
  if (tp=="INT16") return "PLANCK_INT16";
  if (tp=="INT32") return "PLANCK_INT32";
  if (tp=="INT64") return "PLANCK_INT64";
  if (tp=="FLOAT32") return "PLANCK_FLOAT32";
  if (tp=="FLOAT64") return "PLANCK_FLOAT64";
  if (tp=="STRING") return "PLANCK_STRING";
  planck_fail ("unexpected column datatype '" + tp + "'");
  }

void f_prologue(ostream &out)
  {
  out << "! AUTOMATICALLY GENERATED FILE. DO NOT EDIT!\n\n"
         "subroutine fts_fill_cache_from_ddl(handle)\n"
         "  type(fts_handle), intent(inout) :: handle\n\n";
  }

void cxx_prologue(ostream &out)
  {
  out << "// AUTOMATICALLY GENERATED FILE. DO NOT EDIT!\n\n"
         "const std::string ddlstring =\n";
  }

void f_epilogue(ostream &out)
  {
  out << "  call fatal_error('type '''//trim(handle%type)//''' "
           "not found in DDL')\n"
         "end subroutine\n";
  }

void cxx_epilogue(ostream &out)
  {
  out << ";\n";
  out << "const std::string ddl_inheritance =\n";
  for (map<string,string>::const_iterator it=ancestor.begin();
    it!=ancestor.end(); ++it)
    out << "\"" << it->first << " " << it->second << " \"\n";
  out << ";\n";
  }

void read_type (ifstream &inp, ostream &outf, ostream &outc)
  {
  string word;
  int ncols;

  inp >> word;
  if (!inp) return;
  if ((word.size()>0) && (word[word.size()-1]==':'))
    {
    word.erase(word.size()-1);
    string parent;
    inp >> parent;
    ancestor[word]=parent;
    }
  cout << word << endl;
  inp >> ncols;
  outf << "  if (handle%type=='"<<word<<"') then\n"
          "    allocate(handle%columns(0:"<<ncols-1<<"))\n";
  outc << "\"" << word << " " << ncols << " \"\n";
  for (int m=0; m<ncols; ++m)
    {
    string name, unit, type;
    int hdu, col, repc;
    inp >> name >> unit >> type >> hdu >> col >> repc;
    outf << "    handle%columns("<<m<<")=fts_column('"<<name<<"','"<<unit<<"',"
         << f_typeconvert(type)<<","<<hdu<<","<<col<<","<<repc<<",0)\n";
    outc << "\"" << name << " " << unit << " " << type << " " << hdu << " "
         << col << " " << repc << " \"\n";
    }
  outf << "    return\n"
          "  endif\n";
  outc << "\n";
  }

} // unnamed namespace

int main(int argc, const char **argv)
  {
  module_startup ("ddlconvert", argc, argv, 4,
    "<DDL description> <C++ file> <Fortran file>");
  ifstream inp(argv[1]);
  planck_assert (inp,"could not open DDL description");
  ostringstream outc, outf;
  f_prologue(outf);
  cxx_prologue(outc);
  while (!inp.eof())
    {
    read_type (inp,outf,outc);
    }
  f_epilogue(outf);
  cxx_epilogue(outc);

  ofstream outc2(argv[2]);
  planck_assert (outc2,"could not create C++ file");
  ofstream outf2(argv[3]);
  planck_assert (outf2,"could not create Fortran file");
  outc2 << outc.str();
  outf2 << outf.str();
  }
