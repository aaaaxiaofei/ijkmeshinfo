/// \file ijkmeshinfo_compute.txx
/// Templates for computing angles, edge lengths, etc.
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


#ifndef _IJKMESHINFO_COMPUTE_TXX_
#define _IJKMESHINFO_COMPUTE_TXX_

#include "ijkmeshinfo.h"

namespace IJKMESHINFO {

  template <typename MDATA_TYPE,
            typename MTYPE, typename CTYPE,
            typename MINVAL_TYPE0, typename MAXVAL_TYPE0, 
            typename MINVAL_TYPE1, typename MAXVAL_TYPE1, 
            typename PTYPE0, typename PTYPE1,
            typename VTYPE, typename NTYPE0, typename NTYPE1>
  void compute_min_max_values
  (const MDATA_TYPE & mesh_data, const MTYPE & polymesh,
   const COORD_TYPE * vertex_coord, const bool flag_internal,
   const MINVAL_TYPE0 min_init, const MAXVAL_TYPE0 max_init,
   MINVAL_TYPE1 & min_value, MAXVAL_TYPE1 & max_value,
   PTYPE0 & poly_with_min_value, PTYPE1 & poly_with_max_value,
   void compute_min_max_poly_values
   (const MDATA_TYPE &, const VTYPE *, const NTYPE0, const CTYPE *,
    MINVAL_TYPE1 &, MAXVAL_TYPE1 &, NTYPE1 &))
  {
    min_value = min_init;
    max_value = max_init;
    poly_with_min_value = 0;
    poly_with_max_value = 0;

    if (flag_internal) {
      IJK::PROCEDURE_ERROR error("compute_min_max_values");
      if (!check_boundary_facets(mesh_data, error)) { throw error; } 
    }

    bool is_set(false);
    for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

      if (polymesh.poly_data[ipoly].IsDegenerate()) { continue; }

      if (flag_internal) {
        if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      MINVAL_TYPE1 min_value_i;
      MAXVAL_TYPE1 max_value_i;
      NTYPE1 num_values;
      compute_min_max_poly_values
        (mesh_data, polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly),
         vertex_coord, min_value_i, max_value_i, num_values);

      if (num_values > 0) {
        if (!is_set || min_value_i < min_value) { 
          min_value = min_value_i; 
          poly_with_min_value = ipoly;
        }

        if (!is_set || max_value_i > max_value) { 
          max_value = max_value_i; 
          poly_with_max_value = ipoly;
        }

        is_set = true;
      }
    }

  }

}

#endif

