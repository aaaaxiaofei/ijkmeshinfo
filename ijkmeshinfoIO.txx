/// \file ijkmeshinfoIO.txx
/// IO templates for ijkmeshinfo
/// Version 0.2.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2017 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#ifndef _IJKMESHINFO_IO_TXX_
#define _IJKMESHINFO_IO_TXX_

#include <iostream>
#include <string>

#include "ijkmeshinfo.h"

namespace IJKMESHINFO {

  /// Output minimum and maximum values.
  template <typename MDATA_TYPE, typename MTYPE, typename CTYPE,
            typename MINVAL_TYPE, typename MAXVAL_TYPE,
            typename PMIN_TYPE, typename PMAX_TYPE>
  void output_min_max_values
  (std::ostream & out,
   const MDATA_TYPE & mesh_data,
   const MTYPE & polymesh, const CTYPE * vertex_coord, 
   const bool flag_output_min,
   const bool flag_output_max,
   const bool flag_internal,
   const char * value_descriptor,
   void compute_min_max_values
   (const MDATA_TYPE &, const MTYPE &, const CTYPE *, const bool,
    MINVAL_TYPE &, MAXVAL_TYPE &, PMIN_TYPE &, PMAX_TYPE &),
   MINVAL_TYPE & min_val, MAXVAL_TYPE & max_val)
  {
    PMIN_TYPE poly_with_min_val;
    PMAX_TYPE poly_with_max_val;
    const char * internal_descriptor;

    using std::endl;

    compute_min_max_values
      (mesh_data, polymesh, vertex_coord, flag_internal,
       min_val, max_val, poly_with_min_val, poly_with_max_val);

    if (flag_internal) { internal_descriptor = "internal "; }
    else { internal_descriptor = ""; }

    if (flag_output_min) {
      out << "Min " << internal_descriptor << value_descriptor 
          << ": " << min_val << endl; 
    }

    if (flag_output_max) {
      out << "Max " << internal_descriptor << value_descriptor
          << ": " << max_val << endl; 
    }
  }


  /// Output minimum and maximum values.
  /// - Version which does not return min/max values.
  template <typename MDATA_TYPE, typename MTYPE, typename CTYPE,
            typename MINVAL_TYPE, typename MAXVAL_TYPE,
            typename PMIN_TYPE, typename PMAX_TYPE>
  void output_min_max_values
  (std::ostream & out,
   const MDATA_TYPE & mesh_data,
   const MTYPE & polymesh, const CTYPE * vertex_coord, 
   const bool flag_output_min,
   const bool flag_output_max,
   const bool flag_internal,
   const char * value_descriptor,
   void compute_min_max_values
   (const MDATA_TYPE &, const MTYPE &, const CTYPE *, const bool,
    MINVAL_TYPE &, MAXVAL_TYPE &, PMIN_TYPE &, PMAX_TYPE &))
  {
    MINVAL_TYPE min_value;
    MAXVAL_TYPE max_value;
    
    output_min_max_values
      (out, mesh_data, polymesh, vertex_coord,
       flag_output_min, flag_output_max, flag_internal, value_descriptor,
       compute_min_max_values, min_value, max_value);
  }

  /// Print message "(Additional polytopes not listed.)"
  template <typename NTYPE0, typename NTYPE1>
  void output_additional_not_listed_message
  (std::ostream & out,
   const NTYPE0 num_poly_out, const NTYPE1 max_num_poly_out,
   const char * poly_descriptor)
  {
    using std::endl;

    if (num_poly_out > max_num_poly_out) {
      out << "  ..." << endl;
      out << "  (Additional " << num_poly_out-max_num_poly_out
          << " " << poly_descriptor << " not listed.)" << endl;
    }
  }


  /// Output polytopes with min value.
  template <typename MDATA_TYPE, typename MTYPE, typename CTYPE,
            typename VTYPE, 
            typename NTYPE0, typename NTYPE1, typename NTYPE2,
            typename MINVAL_TYPE, typename MAXVAL_TYPE,
            typename PMIN_TYPE, typename PMAX_TYPE>
  void output_polytopes_with_min_value
  (std::ostream & out,
   const MDATA_TYPE & mesh_data,
   const MTYPE & polymesh, const CTYPE * vertex_coord, 
   const bool flag_internal,
   const char * poly_descriptor,
   const char * poly_short_descriptor,
   const char * value_descriptor,
   const NTYPE0 max_num_poly_out,
   void compute_min_max_values
   (const MDATA_TYPE &, const MTYPE &, const CTYPE *, const bool,
    MINVAL_TYPE &, MAXVAL_TYPE &, PMIN_TYPE &, PMAX_TYPE &),
   void compute_min_max_poly_values
   (const MDATA_TYPE &, const VTYPE *, const NTYPE1, const CTYPE *,
    MINVAL_TYPE &, MAXVAL_TYPE &, NTYPE2 &))
  {
    MINVAL_TYPE min_val;
    MAXVAL_TYPE max_val;
    PMIN_TYPE poly_with_min_val;
    PMAX_TYPE poly_with_max_val;
    std::string header_string;

    typedef typename MTYPE::NUMBER_TYPE NUM_POLY_TYPE;

    using std::endl;

    compute_min_max_values
      (mesh_data, polymesh, vertex_coord, flag_internal, 
       min_val, max_val, poly_with_min_val, poly_with_max_val);

    if (flag_internal) 
      { header_string = std::string("Internal ") + poly_descriptor; }
    else {
      header_string = poly_descriptor;
      if (header_string.length() > 0) 
        { header_string[0] = toupper(header_string[0]); }
    }

    out << header_string << " with minimum " << value_descriptor
        << " " << min_val << ":" << endl;

    NTYPE0 num_poly_out = 0;
    for (NUM_POLY_TYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

      if (flag_internal) {
        if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      MINVAL_TYPE min_val_i;
      MAXVAL_TYPE max_val_i;
      NTYPE1 num_val;
      compute_min_max_poly_values
        (mesh_data, polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly),
         vertex_coord, min_val_i, max_val_i, num_val);

      if (num_val > 0) {

        if (ipoly == poly_with_min_val || min_val_i <= min_val) {

          if (num_poly_out < max_num_poly_out) {
            out << "  " << poly_short_descriptor << " " << ipoly << ": ";
            IJK::print_list
              (out, polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly));
            out << "  Min " << value_descriptor << ": " << min_val_i
                << endl;
          }

          num_poly_out++;
        }
      }
    }

    output_additional_not_listed_message
      (out, num_poly_out, max_num_poly_out, poly_descriptor);
  }


  /// Output polytopes with max value.
  template <typename MDATA_TYPE, typename MTYPE, typename CTYPE,
            typename VTYPE, 
            typename NTYPE0, typename NTYPE1, typename NTYPE2,
            typename MINVAL_TYPE, typename MAXVAL_TYPE,
            typename PMIN_TYPE, typename PMAX_TYPE>
  void output_polytopes_with_max_value
  (std::ostream & out,
   const MDATA_TYPE & mesh_data,
   const MTYPE & polymesh, const CTYPE * vertex_coord, 
   const bool flag_internal,
   const char * poly_descriptor,
   const char * poly_short_descriptor,
   const char * value_descriptor,
   const NTYPE0 max_num_poly_out,
   void compute_min_max_values
   (const MDATA_TYPE &, const MTYPE &, const CTYPE *, const bool,
    MINVAL_TYPE &, MAXVAL_TYPE &, PMIN_TYPE &, PMAX_TYPE &),
   void compute_min_max_poly_values
   (const MDATA_TYPE &, const VTYPE *, const NTYPE1, const CTYPE *,
    MINVAL_TYPE &, MAXVAL_TYPE &, NTYPE2 &))
  {
    MINVAL_TYPE min_val;
    MAXVAL_TYPE max_val;
    PMIN_TYPE poly_with_min_val;
    PMAX_TYPE poly_with_max_val;
    std::string header_string;

    typedef typename MTYPE::NUMBER_TYPE NUM_POLY_TYPE;

    using std::endl;

    compute_min_max_values
      (mesh_data, polymesh, vertex_coord, flag_internal, 
       min_val, max_val, poly_with_min_val, poly_with_max_val);

    if (flag_internal) 
      { header_string = std::string("Internal ") + poly_descriptor; }
    else {
      header_string = poly_descriptor;
      if (header_string.length() > 0) 
        { header_string[0] = toupper(header_string[0]); }
    }

    out << header_string << " with maximum " << value_descriptor
        << " " << max_val << ":" << endl;

    NTYPE0 num_poly_out = 0;
    for (NUM_POLY_TYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

      if (flag_internal) {
        if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      MINVAL_TYPE min_val_i;
      MAXVAL_TYPE max_val_i;
      NTYPE1 num_val;
      compute_min_max_poly_values
        (mesh_data, polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly),
         vertex_coord, min_val_i, max_val_i, num_val);

      if (num_val > 0) {

        if (ipoly == poly_with_max_val || max_val_i >= max_val) {

          if (num_poly_out < max_num_poly_out) {
            out << "  " << poly_short_descriptor << " " << ipoly << ": ";
            IJK::print_list
              (out, polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly));
            out << "  Max " << value_descriptor << ": " << max_val_i
                << endl;
          }

          num_poly_out++;
        }
      }
    }

    output_additional_not_listed_message
      (out, num_poly_out, max_num_poly_out, poly_descriptor);
  }

}


#endif
