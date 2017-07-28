/// \file ijkmeshinfo_compute.cxx
/// Compute angles, edge lengths, etc.
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

#include <cmath>

#include "ijkcoord.txx"

#include "ijkmeshinfo_compute.h"

// **************************************************
// COMPUTE ANGLE ROUTINES
// **************************************************

// Compute min/max polygon angles.
// @param flag_internal If true, compute angles for interior polygons.
// @param num_poly_edges If num_poly_edges > 0, compute angles only
//          for polygons with num_poly_edges.
// @param num_poly_edges If num_poly_edges = 0, compute angles 
//          for all polygons.
// @pre mesh_dimension == 2.
void IJKMESHINFO::compute_min_max_polygon_angles
(const int dimension, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  COORD_TYPE cos_min = -1;
  COORD_TYPE cos_max = 1;
  poly_with_min_angle = 0;
  poly_with_max_angle = 0;

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  int num_angle;
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate()) { continue; }

    if (num_poly_edges > 0) {
      if (polymesh.NumPolyVert(ipoly) != num_poly_edges) 
        { continue; }
    }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    IJK::compute_cos_min_max_polygon_angles
      (dimension, polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly),
       vertex_coord, cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      if (cos_min_i > cos_min) {
        cos_min = cos_min_i;
        poly_with_min_angle = ipoly;
      }

      if (cos_max_i < cos_max) {
        cos_max = cos_max_i;
        poly_with_max_angle = ipoly;
      }
    }
  }


  min_angle = std::acos(cos_min) * 180.0/M_PI;
  max_angle = std::acos(cos_max) * 180.0/M_PI;
}


// Compute min/max polygon angles.
// - Version with num_poly_edges set to 0.
void IJKMESHINFO::compute_min_max_polygon_angles
(const int dimension, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  compute_min_max_polygon_angles
    (dimension, polymesh, vertex_coord, flag_internal, 0, 
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute min/max polygon angles.
// - Version without return arguments poly_with_min_angle and
//   poly_with_max_angle.
void IJKMESHINFO::compute_min_max_polygon_angles
(const int dimension, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_polygon_angles
    (dimension, polymesh, vertex_coord, flag_internal, num_poly_edges, 
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute min/max polygon angles.
// - Version without return arguments poly_with_min_angle and
//   poly_with_max_angle.
// - Version with num_poly_edges set to 0.
void IJKMESHINFO::compute_min_max_polygon_angles
(const int dimension, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  compute_min_max_polygon_angles
    (dimension, polymesh, vertex_coord, flag_internal, 0, 
     min_angle, max_angle);
}


// Compute number of angles less than or equal to min_angle and
//   greater than or equal to max_angle.
// @pre mesh_dimension == 2.
void IJKMESHINFO::compute_num_polygon_angles
(const int dimension,
 const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, 
 const ANGLE_TYPE min_angle, const ANGLE_TYPE max_angle,
 int & num_le, int & num_ge)
{
  // Initialize to zero.
  num_le = 0;
  num_ge = 0;

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  COORD_TYPE cos_min = cos(min_angle*M_PI/180.0);
  COORD_TYPE cos_max = cos(max_angle*M_PI/180.0);
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    int num_angle;
    IJK::compute_cos_min_max_polygon_angles
      (dimension, polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly),
       vertex_coord, cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      if (cos_min_i >= cos_min) { num_le++; }
      if (cos_max_i <= cos_max) { num_ge++; }
    }

  }

}

// Compute number of dihedral angles less than or equal to min_angle and
//   greater than or equal to max_angle.
// @pre mesh_dimension == 3.
void IJKMESHINFO::compute_num_dihedral_angles
(const int dimension,
 const int mesh_dimension, 
 const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, 
 const ANGLE_TYPE min_angle, const ANGLE_TYPE max_angle,
 int & num_le, int & num_ge)
{
  // Initialize to zero.
  num_le = 0;
  num_ge = 0;

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  COORD_TYPE cos_min = cos(min_angle*M_PI/180.0);
  COORD_TYPE cos_max = cos(max_angle*M_PI/180.0);
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    int num_angle;
    IJK::compute_cos_min_max_tetrahedron_dihedral_angles
      (dimension, polymesh.VertexList(ipoly), vertex_coord,
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      if (cos_min_i >= cos_min) { num_le++; }
      if (cos_max_i <= cos_max) { num_ge++; }
    }

  }

}


// Compute cosine of min/max  angle over all tetrahedron facets.
void IJKMESHINFO::compute_cos_min_max_tetrahedron_facet_angles
(const int dimension,
 const VERTEX_INDEX tetrahedron_vert[], const COORD_TYPE * vertex_coord,
 COORD_TYPE & cos_min, COORD_TYPE & cos_max,
 int & num_angle)
{
  const int NUM_VERT_PER_TETRAHEDRON(4);

  cos_min = -1;
  cos_max = 1;
  num_angle = 0;

  for (int i0 = 0; i0 < NUM_VERT_PER_TETRAHEDRON; i0++) {
    const int iv0 = tetrahedron_vert[i0];

    const int i1 = (i0+1)%NUM_VERT_PER_TETRAHEDRON;
    const int i2 = (i0+2)%NUM_VERT_PER_TETRAHEDRON;
    const int i3 = (i0+3)%NUM_VERT_PER_TETRAHEDRON;

    const int iv1 = tetrahedron_vert[i1];
    const int iv2 = tetrahedron_vert[i2];
    const int iv3 = tetrahedron_vert[i3];

    double cos_angle;
    bool flag_duplicate_point;

    // Compute min/max angles of triangle (iv0, iv1, iv3).
    IJK::compute_cos_triangle_angle_coord_list
      (dimension, vertex_coord, iv0, iv1, iv2, cos_angle, 
       flag_duplicate_point);

    if (!flag_duplicate_point) {
      num_angle++;
      if (cos_angle > cos_min) { cos_min = cos_angle; }
      if (cos_angle < cos_max) { cos_max = cos_angle; }
    }
  }
}


// Compute number of tetrahedra facet angles less than or equal 
//   to min_angle and  greater than or equal to max_angle.
// @pre mesh_dimension == 3.
void IJKMESHINFO::compute_num_tetrahedra_facet_angles
(const int dimension, const int mesh_dimension, 
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal, 
 const ANGLE_TYPE min_angle, const ANGLE_TYPE max_angle,
 int & num_le, int & num_ge)
{
  // Initialize to zero.
  num_le = 0;
  num_ge = 0;

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  COORD_TYPE cos_min = cos(min_angle*M_PI/180.0);
  COORD_TYPE cos_max = cos(max_angle*M_PI/180.0);
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    int num_angle;
    compute_cos_min_max_tetrahedron_facet_angles
      (dimension, polymesh.VertexList(ipoly), vertex_coord,
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      if (cos_min_i >= cos_min) { num_le++; }
      if (cos_max_i <= cos_max) { num_ge++; }
    }

  }

}


// Compute min/max angle over all tetrahedra facets.
// @pre Tetrahedra are in 3D.
void IJKMESHINFO::compute_min_max_tetrahedra_facet_angles
(const int dimension, const int mesh_dimension,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  COORD_TYPE cos_min = -1;
  COORD_TYPE cos_max = 1;

  poly_with_min_angle = 0;
  poly_with_max_angle = 0;

  int num_angle;
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].is_degenerate) { continue; }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    compute_cos_min_max_tetrahedron_facet_angles
      (dimension, polymesh.VertexList(ipoly), vertex_coord, 
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      if (cos_min_i > cos_min) { 
        cos_min = cos_min_i;
        poly_with_min_angle = ipoly;
      }

      if (cos_max_i < cos_max) {
        cos_max = cos_max_i;
        poly_with_max_angle = ipoly;
      }
    }
  }

  min_angle = std::acos(cos_min) * 180.0/M_PI;
  max_angle = std::acos(cos_max) * 180.0/M_PI;
}


// Compute min/max angle over all tetrahedra facets.
void IJKMESHINFO::compute_min_max_tetrahedra_facet_angles
(const int dimension, const int mesh_dimension,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_tetrahedra_facet_angles
    (dimension, mesh_dimension, polymesh, vertex_coord, flag_internal,
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// **************************************************
// COMPUTE DIHEDRAL ANGLE ROUTINES
// **************************************************

// Compute min/max dihedral angles of simplices
void IJKMESHINFO::compute_min_max_dihedral_angles
(const int dimension,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle,
 const bool flag_internal)
{
  COORD_TYPE cos_min = -1;
  COORD_TYPE cos_max = 1;
  poly_with_min_angle = 0;
  poly_with_max_angle = 0;

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  int num_angle;
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate()) { continue; }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    IJK::compute_cos_min_max_tetrahedron_dihedral_angles
      (dimension, polymesh.VertexList(ipoly), vertex_coord,
       cos_min_i, cos_max_i, num_angle);


    if (num_angle > 0) {
      if (cos_min_i > cos_min) { 
        cos_min = cos_min_i; 
        poly_with_min_angle = ipoly;
      }

      if (cos_max_i < cos_max) {
        cos_max = cos_max_i;
        poly_with_max_angle = ipoly;
      }
    }
  }

  min_angle = std::acos(cos_min) * 180.0/M_PI;
  max_angle = std::acos(cos_max) * 180.0/M_PI;
}


void IJKMESHINFO::compute_min_max_dihedral_angles
(const int dimension,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 const bool flag_internal)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_dihedral_angles
    (dimension, polymesh, vertex_coord, min_angle, max_angle, 
     poly_with_min_angle, poly_with_max_angle, flag_internal);
}
