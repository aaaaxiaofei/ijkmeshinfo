/// \file ijkmeshinfo_compute.cxx
/// Compute angles, edge lengths, etc.
/// Version 0.2.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2017-2018 Rephael Wenger

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
#include "ijkmeshinfo_compute.txx"


// **************************************************
// COMPUTE ANGLE ROUTINES
// **************************************************

// Compute min/max polygon angles.
// @param flag_internal If true, compute angles for interior polygons.
// @param num_poly_vert If num_poly_vert > 0, compute angles only
//          for polygons with num_poly_vert.
// @param num_poly_edges If num_poly_vert = 0, compute angles 
//          for all polygons.
// @pre mesh_dimension == 2.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_vert,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  const int dimension = mesh_data.dimension;
  IJK::PROCEDURE_ERROR error("compute_min_max_polygon_angles");

  poly_with_min_angle = 0;
  poly_with_max_angle = 0;

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  COORD_TYPE min_cos, max_cos;
  if (num_poly_vert == 0) {
    // Note: Poly with max angle has min_cos.
    //       Poly with min angle has max_cos.
    compute_min_max_plist_values
      (mesh_data, polymesh, vertex_coord, flag_internal, 
       1, -1, min_cos, max_cos, poly_with_max_angle, poly_with_min_angle,
       compute_min_max_cos_polygon_angles);
  }
  else {
    // Note: Poly with max angle has min_cos.
    //       Poly with min angle has max_cos.
    compute_min_max_plist_values_select_poly_by_numv
      (mesh_data, polymesh, vertex_coord, flag_internal, 1, -1, num_poly_vert,
       min_cos, max_cos, poly_with_max_angle, poly_with_min_angle,
       compute_min_max_cos_polygon_angles);
  }

  min_angle = std::acos(max_cos) * 180.0/M_PI;
  max_angle = std::acos(min_cos) * 180.0/M_PI;
}


// Compute min/max polygon angles.
// - Version with num_poly_edges set to 0.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  compute_min_max_polygon_angles
    (mesh_data, polymesh, vertex_coord, flag_internal, 0, 
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute min/max polygon angles.
// - Version without return arguments poly_with_min_angle and
//   poly_with_max_angle.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_polygon_angles
    (mesh_data, polymesh, vertex_coord, flag_internal, num_poly_edges, 
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute min/max polygon angles.
// - Version without return arguments poly_with_min_angle and
//   poly_with_max_angle.
// - Version with num_poly_edges set to 0.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal, ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  compute_min_max_polygon_angles
    (mesh_data, polymesh, vertex_coord, flag_internal, 0, 
     min_angle, max_angle);
}


// Compute min/max polygon angles.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle, int & num_angle)
{
  const int dimension = mesh_data.dimension;

  // Initialize
  min_angle = 0;
  max_angle = 180;

  COORD_TYPE cos_min, cos_max;
  IJK::compute_cos_min_max_polygon_angles
    (dimension, poly_vert, num_vert, vertex_coord, 
     cos_min, cos_max, num_angle);

  if (num_angle > 0) {
    min_angle = std::acos(cos_min) * 180.0/M_PI;
    max_angle = std::acos(cos_max) * 180.0/M_PI;
  }
}


// Compute min/max cosine of the angles of a single polygon.
// - Note: acos(min_cos) is the max angle.
//         acos(max_cos) is the min angle.
void IJKMESHINFO::compute_min_max_cos_polygon_angles
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_cos, COORD_TYPE & max_cos, int & num_angle)
{
  const int dimension = mesh_data.dimension;

  IJK::compute_cos_min_max_polygon_angles
    (dimension, poly_vert, num_vert, vertex_coord, 
     max_cos, min_cos, num_angle);
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
void IJKMESHINFO::compute_num_tetmesh_dihedral_angles
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


// Compute cosine of min/max angle over all tetrahedron facets.
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


// Compute min/max angle over all tetrahedron facets.
void IJKMESHINFO::compute_min_max_tetrahedron_facet_angles
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX tetrahedron_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & num_angle)
{
  const int dimension = mesh_data.dimension;

  // Initialize
  min_angle = 0;
  max_angle = 180;

  COORD_TYPE cos_min, cos_max;
  compute_cos_min_max_tetrahedron_facet_angles
    (dimension, tetrahedron_vert, vertex_coord, cos_min, cos_max, num_angle);

  if (num_angle > 0) {
    min_angle = std::acos(cos_min) * 180.0/M_PI;
    max_angle = std::acos(cos_max) * 180.0/M_PI;
  }
}


// Compute min/max angle over all tetrahedra facets.
// @pre Tetrahedra are in 3D.
void IJKMESHINFO::compute_min_max_tetrahedra_facet_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  const int dimension = mesh_data.dimension;

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
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_tetrahedra_facet_angles
    (mesh_data, polymesh, vertex_coord, flag_internal,
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
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


// **************************************************
// COMPUTE DIHEDRAL ANGLE ROUTINES
// **************************************************

// Compute min/max dihedral angles of tetrahedra.
void IJKMESHINFO::compute_min_max_tetmesh_dihedral_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  const int dimension = mesh_data.dimension;
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


void IJKMESHINFO::compute_min_max_tetmesh_dihedral_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_tetmesh_dihedral_angles
    (mesh_data, polymesh, vertex_coord, flag_internal,
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute min/max dihedral angles of a tetrahedron.
void IJKMESHINFO::compute_min_max_tetrahedron_dihedral_angles
(const MESH_DATA & mesh_data, const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle, int & num_angle)
{
  const int dimension = mesh_data.dimension;

  // Initialize
  min_angle = 0;
  max_angle = 180;

  COORD_TYPE cos_min, cos_max;
  IJK::compute_cos_min_max_tetrahedron_dihedral_angles
    (dimension, poly_vert, vertex_coord, cos_min, cos_max, num_angle);

  if (num_angle > 0) {
    min_angle = std::acos(cos_min) * 180.0/M_PI;
    max_angle = std::acos(cos_max) * 180.0/M_PI;
  }
}


// **************************************************
// COMPUTE EDGE LENGTH ROUTINES
// **************************************************


// Compute min/max edge lengths of polygons.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  compute_min_max_plist_values
    (mesh_data, polymesh, vertex_coord, flag_internal,
     0, 0, min_length, max_length,
     poly_with_min_edge_length, poly_with_max_edge_length,
     compute_min_max_polygon_edge_lengths);
}

// Compute min/max edge lengths of polygons.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths_select_poly_by_numv
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_vert,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  if (num_poly_vert == 0) {
    compute_min_max_polygon_edge_lengths
      (mesh_data, polymesh, vertex_coord,
       flag_internal, min_length, max_length,
       poly_with_min_edge_length, poly_with_max_edge_length);
  }
  else {
    compute_min_max_plist_values_select_poly_by_numv
      (mesh_data, polymesh, vertex_coord, flag_internal,
       0, 0, num_poly_vert, min_length, max_length,
       poly_with_min_edge_length, poly_with_max_edge_length,
       compute_min_max_polygon_edge_lengths);
  }
}

// Compute min/max edge lengths of polygons.
// - Version which does not return polygons with min/max edge lengths.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths_select_poly_by_numv
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_vert,
 COORD_TYPE & min_length, COORD_TYPE & max_length)
{
  int poly_with_min_edge_length, poly_with_max_edge_length;

  compute_min_max_polygon_edge_lengths_select_poly_by_numv
    (mesh_data, polymesh, vertex_coord, flag_internal, num_poly_vert,
     min_length, max_length, 
     poly_with_min_edge_length, poly_with_max_edge_length);
}

// Compute min/max edge lengths of one polygon.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;
  COORD_TYPE edge_length_squared;

  // Initialize.
  min_length = 0;
  max_length = 0;

  COORD_TYPE min_edge_length_squared = 0;
  COORD_TYPE max_edge_length_squared = 0;
  num_lengths = 0;

  bool flag_set = false;
  for (int i0 = 0; i0 < num_vert; i0++) {
    const int i1 = (i0+1)%num_vert;

    const VERTEX_INDEX iv0 = poly_vert[i0];
    const VERTEX_INDEX iv1 = poly_vert[i1];

    const COORD_TYPE * v0coord = vertex_coord+iv0*dimension;
    const COORD_TYPE * v1coord = vertex_coord+iv1*dimension;

    IJK::compute_distance_squared
      (dimension, v0coord, v1coord, edge_length_squared);

    if (!flag_set || edge_length_squared < min_edge_length_squared) 
      { min_edge_length_squared = edge_length_squared; }
    if (!flag_set || edge_length_squared > max_edge_length_squared) 
      { max_edge_length_squared = edge_length_squared; }

    flag_set = true;
    num_lengths++;
  }

  if (num_lengths > 0) {
    min_length = std::sqrt(min_edge_length_squared);
    max_length = std::sqrt(max_edge_length_squared);
  }
}


// Compute min/max edge lengths of tetrahedra.
void IJKMESHINFO::compute_min_max_tetrahedra_edge_lengths
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  const int dimension = mesh_data.dimension;

  compute_min_max_plist_values
    (mesh_data, polymesh, vertex_coord, flag_internal,
     0, 0, min_length, max_length,
     poly_with_min_edge_length, poly_with_max_edge_length,
     compute_min_max_tetrahedron_edge_lengths);
}


// Compute min/max edge lengths of one tetrahedron.
void IJKMESHINFO::compute_min_max_tetrahedron_edge_lengths
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX tet_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord, 
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;
  const int NUM_VERT_PER_TETRAHEDRON(4);
  COORD_TYPE edge_length_squared;

  // Initialize.
  min_length = 0;
  max_length = 0;

  COORD_TYPE min_edge_length_squared = 0;
  COORD_TYPE max_edge_length_squared = 0;
  num_lengths = 0;

  bool flag_set = false;
  for (int i0 = 0; i0 < NUM_VERT_PER_TETRAHEDRON; i0++) {
    for (int i1 = i0+1; i1 < NUM_VERT_PER_TETRAHEDRON; i1++) {

      const VERTEX_INDEX iv0 = tet_vert[i0];
      const VERTEX_INDEX iv1 = tet_vert[i1];

      const COORD_TYPE * v0coord = vertex_coord+iv0*dimension;
      const COORD_TYPE * v1coord = vertex_coord+iv1*dimension;

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, edge_length_squared);

      if (!flag_set || edge_length_squared < min_edge_length_squared) 
        { min_edge_length_squared = edge_length_squared; }
      if (!flag_set || edge_length_squared > max_edge_length_squared) 
        { max_edge_length_squared = edge_length_squared; }

      flag_set = true;
      num_lengths++;
    }
  }

  if (num_lengths > 0) {
    min_length = std::sqrt(min_edge_length_squared);
    max_length = std::sqrt(max_edge_length_squared);
  }
}


// Compute min/max edge lengths of a polymesh of simplices.
void IJKMESHINFO::compute_min_max_simplices_edge_lengths
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  compute_min_max_plist_values
    (mesh_data, polymesh, vertex_coord, flag_internal,
     0, 0, min_length, max_length,
     poly_with_min_edge_length, poly_with_max_edge_length,
     compute_min_max_simplex_edge_lengths);
}


// Compute min/max edge lengths of one simplex.
void IJKMESHINFO::compute_min_max_simplex_edge_lengths
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX simplex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;
  COORD_TYPE edge_length_squared;

  // Initialize.
  min_length = 0;
  max_length = 0;

  COORD_TYPE min_edge_length_squared = 0;
  COORD_TYPE max_edge_length_squared = 0;
  num_lengths = 0;

  bool flag_set = false;
  for (int i0 = 0; i0 < num_vert; i0++) {
    for (int i1 = i0+1; i1 < num_vert; i1++) {

      const VERTEX_INDEX iv0 = simplex_vert[i0];
      const VERTEX_INDEX iv1 = simplex_vert[i1];

      const COORD_TYPE * v0coord = vertex_coord+iv0*dimension;
      const COORD_TYPE * v1coord = vertex_coord+iv1*dimension;

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, edge_length_squared);

      if (!flag_set || edge_length_squared < min_edge_length_squared) 
        { min_edge_length_squared = edge_length_squared; }
      if (!flag_set || edge_length_squared > max_edge_length_squared) 
        { max_edge_length_squared = edge_length_squared; }

      flag_set = true;
      num_lengths++;
    }
  }

  if (num_lengths > 0) {
    min_length = std::sqrt(min_edge_length_squared);
    max_length = std::sqrt(max_edge_length_squared);
  }
}


// Compute min/max edge lengths of hexahedra.
void IJKMESHINFO::compute_min_max_hexahedra_edge_lengths
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  compute_min_max_plist_values
    (mesh_data, polymesh, vertex_coord, flag_internal,
     0, 0, min_length, max_length,
     poly_with_min_edge_length, poly_with_max_edge_length,
     compute_min_max_hexahedron_edge_lengths);
}


// Compute min/max edge lengths of one hexahedron.
void IJKMESHINFO::compute_min_max_hexahedron_edge_lengths
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;
  const static CUBE_TYPE cube(DIM3);
  COORD_TYPE edge_length_squared;

  // Initialize.
  min_length = 0;
  max_length = 0;

  COORD_TYPE min_edge_length_squared = 0;
  COORD_TYPE max_edge_length_squared = 0;
  num_lengths = 0;

  bool flag_set = false;
  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    const int iend0 = cube.EdgeEndpoint(ie, 0);
    const int iend1 = cube.EdgeEndpoint(ie, 1);

    const VERTEX_INDEX iv0 = hex_vert[iend0];
    const VERTEX_INDEX iv1 = hex_vert[iend1];

    const COORD_TYPE * v0coord = vertex_coord+iv0*dimension;
    const COORD_TYPE * v1coord = vertex_coord+iv1*dimension;

    IJK::compute_distance_squared
      (dimension, v0coord, v1coord, edge_length_squared);

    if (!flag_set || edge_length_squared < min_edge_length_squared) 
      { min_edge_length_squared = edge_length_squared; }
    if (!flag_set || edge_length_squared > max_edge_length_squared) 
      { max_edge_length_squared = edge_length_squared; }

    flag_set = true;
    num_lengths++;
  }

  if (num_lengths > 0) {
    min_length = std::sqrt(min_edge_length_squared);
    max_length = std::sqrt(max_edge_length_squared);
  }

}


// **************************************************
// COMPUTE JACOBIAN DETERMINANTS
// **************************************************

// Compute min/max Jacobian matrix determinants of hexahedra.
void IJKMESHINFO::compute_min_max_hexahedra_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant)
{
 compute_min_max_plist_values
    (mesh_data, polymesh, vertex_coord, flag_internal,
     0, 0, min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant,
     compute_min_max_hexahedron_Jacobian_determinants);
}


// Compute min/max Jacobian matrix determinants of hexahedra.
// - Version which does not return hexahedra containing min/max 
//   Jacobian determinants.
void IJKMESHINFO::compute_min_max_hexahedra_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  int poly_with_min_Jacobian_determinant;
  int poly_with_max_Jacobian_determinant;

  compute_min_max_hexahedra_Jacobian_determinants
    (mesh_data, polymesh, vertex_coord, flag_internal,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant);
}

// Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hexahedron_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_determinant,
 COORD_TYPE & max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);

  IJK::compute_min_max_hexahedron_Jacobian_determinant_3D
    (hex_vert, mesh_data.orientation, vertex_coord, cube,
     min_Jacobian_determinant, max_Jacobian_determinant);

  num_Jacobian_determinants = 9;
}


// Compute min/max Jacobian matrix determinants of hexahedra vertices.
// - Version which returns vertices with min/max Jacobian determinants.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant,
 int & vert_with_min_Jacobian_determinant, 
 int & vert_with_max_Jacobian_determinant)
{
  compute_min_max_vlist_values
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     0, 0, min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant, vert_with_max_Jacobian_determinant,
     compute_min_max_hex_vert_Jacobian_determinants);
}


// Compute min/max of the Jacobian matrix determinants of hexahedra vertices.
// - Version which does not return vertices or polytopes with min/max values.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly, const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  int poly_with_min_Jacobian_determinant;
  int poly_with_max_Jacobian_determinant;
  int vert_with_min_Jacobian_determinant;
  int vert_with_max_Jacobian_determinant;

  compute_min_max_hex_vert_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant,
     poly_with_max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant,
     vert_with_max_Jacobian_determinant);
}


// Compute min/max Jacobian matrix determinants of a vertex 
//   in a hexahedral mesh.
// - Version with input vertex-poly incidence data structure.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const VERTEX_INDEX iv0,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);
  const int NUM_CUBE_VERT = cube.NumVertices();
  const int num_incident_poly = 
    vertex_poly_incidence.NumIncidentPoly(iv0);
  COORD_TYPE Jacobian_determinant;

  // Initialize
  num_Jacobian_determinants = 0;
  min_Jacobian_determinant = 0.0;
  max_Jacobian_determinant = 0.0;
  poly_with_min_Jacobian_determinant = 0;
  poly_with_max_Jacobian_determinant = 0;

  if (flag_internal_vert) {
    if (polymesh.vertex_data[iv0].OnBoundary()) { 
      // Vertex iv0 is not internal.
      return;
    }
  }

  for (int k = 0; k < num_incident_poly; k++) {

    const int kpoly = vertex_poly_incidence.IncidentPoly(iv0, k);
    const int kloc = vertex_poly_incidence.VertexLocInPolyVertexList(iv0, k);
    IJK::compute_Jacobian_determinant_at_hex_vertex_3D
      (polymesh.VertexList(kpoly), mesh_data.orientation,
       vertex_coord, cube, kloc, Jacobian_determinant);

    if (flag_internal_poly) {
      if (polymesh.poly_data[kpoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    if (num_Jacobian_determinants == 0) {
      min_Jacobian_determinant = Jacobian_determinant;
      max_Jacobian_determinant = Jacobian_determinant;
      poly_with_min_Jacobian_determinant = kpoly;
      poly_with_max_Jacobian_determinant = kpoly;
    }
    else {

      if (Jacobian_determinant < min_Jacobian_determinant) {
        min_Jacobian_determinant = Jacobian_determinant;
        poly_with_min_Jacobian_determinant = kpoly;
      }

      if (Jacobian_determinant > max_Jacobian_determinant) {
        max_Jacobian_determinant = Jacobian_determinant;
        poly_with_max_Jacobian_determinant = kpoly;
      }
    }

    num_Jacobian_determinants++;
  }

}


// **************************************************
// COMPUTE NORMALIZED JACOBIAN DETERMINANTS
// **************************************************

// Compute min/max normalized Jacobian matrix determinants of hexahedra.
void IJKMESHINFO::compute_min_max_hexahedra_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant)
{
 compute_min_max_plist_values
    (mesh_data, polymesh, vertex_coord, flag_internal,
     0, 0, min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant,
     compute_min_max_hexahedron_normalized_Jacobian_determinants);
}


// Compute min/max normalized Jacobian matrix determinants of hexahedra.
// - Version which does not return hexahedra containing min/max 
//   normalized Jacobian determinants.
void IJKMESHINFO::compute_min_max_hexahedra_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  int poly_with_min_value;
  int poly_with_max_value;

  compute_min_max_hexahedra_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_coord, flag_internal,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_value, poly_with_max_value);
}


// Compute min/max of the nine normalized Jacobian matrix determinants 
//   of a hexahedron.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hexahedron_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_determinant,
 COORD_TYPE & max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);
  const COORD_TYPE max_small_magnitude(0.0);

  IJK::compute_min_max_hexahedron_normalized_Jacobian_determinant_3D
    (hex_vert, mesh_data.orientation, vertex_coord, cube, max_small_magnitude,
     min_Jacobian_determinant, max_Jacobian_determinant, 
     num_Jacobian_determinants);
}


// Compute min/max of the eight normalized Jacobian matrix determinants 
//   at the eight vertices of a hexahedron.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_determinant,
 COORD_TYPE & max_Jacobian_determinant,
 VERTEX_INDEX & vert_with_min_Jacobian_determinant,
 VERTEX_INDEX & vert_with_max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);

  const COORD_TYPE max_small_magnitude(0.0);

  IJK::compute_min_max_hex_vert_normalized_Jacobian_determinant_3D
    (hex_vert, mesh_data.orientation, vertex_coord, cube, max_small_magnitude,
     min_Jacobian_determinant, max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant, 
     vert_with_max_Jacobian_determinant,
     num_Jacobian_determinants);
}


// Compute min/max normalized Jacobian matrix determinants 
//   of hexahedra vertices.
// - Version which returns vertices with min/max Jacobian determinants.
void IJKMESHINFO::compute_min_max_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, 
 COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant,
 int & vert_with_min_Jacobian_determinant, 
 int & vert_with_max_Jacobian_determinant)
{
  compute_min_max_vlist_values
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     0, 0, min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant, vert_with_max_Jacobian_determinant,
     compute_min_max_hex_vert_normalized_Jacobian_determinants);
}


// Compute min/max of the normalized Jacobian matrix determinants 
//   of hexahedra vertices.
// - Version which does not return vertices or polytopes with min/max values.
void IJKMESHINFO::compute_min_max_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly, const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  int poly_with_min_Jacobian_determinant;
  int poly_with_max_Jacobian_determinant;
  int vert_with_min_Jacobian_determinant;
  int vert_with_max_Jacobian_determinant;

  compute_min_max_hex_vert_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant,
     poly_with_max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant,
     vert_with_max_Jacobian_determinant);
}


// Compute min/max normalized Jacobian matrix determinants of a vertex 
//   in a hexahedral mesh.
// - Version with input vertex-poly incidence data structure.
void IJKMESHINFO::compute_min_max_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const VERTEX_INDEX iv0,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);
  const int NUM_CUBE_VERT = cube.NumVertices();
  const int num_incident_poly = 
    vertex_poly_incidence.NumIncidentPoly(iv0);
  const COORD_TYPE max_small_magnitude(0.0);
  bool flag_zero;
  COORD_TYPE Jacobian_determinant;

  // Initialize
  num_Jacobian_determinants = 0;
  min_Jacobian_determinant = 0.0;
  max_Jacobian_determinant = 0.0;
  poly_with_min_Jacobian_determinant = 0;
  poly_with_max_Jacobian_determinant = 0;

  if (flag_internal_vert) {
    if (polymesh.vertex_data[iv0].OnBoundary()) { 
      // Vertex iv0 is not internal.
      return;
    }
  }

  for (int k = 0; k < num_incident_poly; k++) {

    const int kpoly = vertex_poly_incidence.IncidentPoly(iv0, k);
    const int kloc = vertex_poly_incidence.VertexLocInPolyVertexList(iv0, k);
    IJK::compute_normalized_Jacobian_determinant_at_hex_vertex_3D
      (polymesh.VertexList(kpoly), mesh_data.orientation,
       vertex_coord, cube, kloc, max_small_magnitude,
       Jacobian_determinant, flag_zero);

    if (flag_internal_poly) {
      if (polymesh.poly_data[kpoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    if (flag_zero) { continue; }

    if (num_Jacobian_determinants == 0) {
      min_Jacobian_determinant = Jacobian_determinant;
      max_Jacobian_determinant = Jacobian_determinant;
      poly_with_min_Jacobian_determinant = kpoly;
      poly_with_max_Jacobian_determinant = kpoly;
    }
    else {

      if (Jacobian_determinant < min_Jacobian_determinant) {
        min_Jacobian_determinant = Jacobian_determinant;
        poly_with_min_Jacobian_determinant = kpoly;
      }

      if (Jacobian_determinant > max_Jacobian_determinant) {
        max_Jacobian_determinant = Jacobian_determinant;
        poly_with_max_Jacobian_determinant = kpoly;
      }
    }

    num_Jacobian_determinants++;
  }

}


// **************************************************
// COMPUTE SHAPE METRIC BASED ON JACOBIAN MATRICES
// **************************************************

// Compute min/max hexahedra shape metric based on Jacobian matrices.
void IJKMESHINFO::compute_min_max_hexahedra_Jacobian_shape
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape,
 int & poly_with_min_Jacobian_shape, 
 int & poly_with_max_Jacobian_shape)
{
 compute_min_max_plist_values
    (mesh_data, polymesh, vertex_coord, flag_internal,
     0, 0, min_Jacobian_shape, max_Jacobian_shape,
     poly_with_min_Jacobian_shape, poly_with_max_Jacobian_shape,
     compute_min_max_hexahedron_Jacobian_shape);
}


// Compute min/max hexahedra shape metric based on Jacobian matrices.
// - Version which does not return hexahedra containing min/max shape.
void IJKMESHINFO::compute_min_max_hexahedra_Jacobian_shape
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape)
{
  int poly_with_min_value;
  int poly_with_max_value;

  compute_min_max_hexahedra_Jacobian_shape
    (mesh_data, polymesh, vertex_coord, flag_internal,
     min_Jacobian_shape, max_Jacobian_shape,
     poly_with_min_value, poly_with_max_value);
}


// Compute min/max hexahedron shape metric based on Jacobian matrices.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hexahedron_Jacobian_shape
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_shape,
 COORD_TYPE & max_Jacobian_shape,
 int & num_Jacobian_shapes)
{
  const static CUBE_TYPE cube(DIM3);
  const COORD_TYPE max_small_magnitude(0.0);

  IJK::compute_min_max_hexahedron_Jacobian_shape_3D
    (hex_vert, vertex_coord, cube, max_small_magnitude,
     min_Jacobian_shape, max_Jacobian_shape,
     num_Jacobian_shapes);
}


// Compute min/max of the eight Jacobian shape metrics
//   at the eight vertices of a hexahedron.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_shape,
 COORD_TYPE & max_Jacobian_shape,
 VERTEX_INDEX & vert_with_min_Jacobian_shape,
 VERTEX_INDEX & vert_with_max_Jacobian_shape,
 int & num_Jacobian_shapes)
{
  const static CUBE_TYPE cube(DIM3);

  const COORD_TYPE max_small_magnitude(0.0);

  IJK::compute_min_max_hex_vert_Jacobian_shape_3D
    (hex_vert, vertex_coord, cube, max_small_magnitude,
     min_Jacobian_shape, max_Jacobian_shape,
     vert_with_min_Jacobian_shape, 
     vert_with_max_Jacobian_shape,
     num_Jacobian_shapes);
}


// Compute min/max shape metrics based on Jacobian matrices at hex vertices.
// - Version which returns vertices with min/max shape metric.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_shape, 
 COORD_TYPE & max_Jacobian_shape,
 int & poly_with_min_Jacobian_shape, 
 int & poly_with_max_Jacobian_shape,
 int & vert_with_min_Jacobian_shape, 
 int & vert_with_max_Jacobian_shape)
{
  compute_min_max_vlist_values
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     0, 0, min_Jacobian_shape, max_Jacobian_shape,
     poly_with_min_Jacobian_shape, poly_with_max_Jacobian_shape,
     vert_with_min_Jacobian_shape, vert_with_max_Jacobian_shape,
     compute_min_max_hex_vert_Jacobian_shape);
}


// Compute min/max shape metrics based on Jacobian matrices at hex vertices.
// - Version which does not return vertices or polytopes with min/max values.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly, const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape)
{
  int poly_with_min_Jacobian_shape;
  int poly_with_max_Jacobian_shape;
  int vert_with_min_Jacobian_shape;
  int vert_with_max_Jacobian_shape;

  compute_min_max_hex_vert_Jacobian_shape
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_shape, max_Jacobian_shape,
     poly_with_min_Jacobian_shape,
     poly_with_max_Jacobian_shape,
     vert_with_min_Jacobian_shape,
     vert_with_max_Jacobian_shape);
}


// Compute min/max shape metrics based on Jacobian matrices at hex vertices.
// - Version with input vertex-poly incidence data structure.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const VERTEX_INDEX iv0,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape,
 int & poly_with_min_Jacobian_shape, 
 int & poly_with_max_Jacobian_shape,
 int & num_Jacobian_shapes)
{
  const static CUBE_TYPE cube(DIM3);
  const int NUM_CUBE_VERT = cube.NumVertices();
  const int num_incident_poly = 
    vertex_poly_incidence.NumIncidentPoly(iv0);
  const COORD_TYPE max_small_magnitude(0.0);
  bool flag_zero;
  COORD_TYPE Jacobian_shape;

  // Initialize
  num_Jacobian_shapes = 0;
  min_Jacobian_shape = 0.0;
  max_Jacobian_shape = 0.0;
  poly_with_min_Jacobian_shape = 0;
  poly_with_max_Jacobian_shape = 0;

  if (flag_internal_vert) {
    if (polymesh.vertex_data[iv0].OnBoundary()) { 
      // Vertex iv0 is not internal.
      return;
    }
  }

  for (int k = 0; k < num_incident_poly; k++) {

    const int kpoly = vertex_poly_incidence.IncidentPoly(iv0, k);
    const int kloc = vertex_poly_incidence.VertexLocInPolyVertexList(iv0, k);
    IJK::compute_Jacobian_shape_at_hex_vertex_3D
      (polymesh.VertexList(kpoly), vertex_coord, cube, kloc, 
       max_small_magnitude, Jacobian_shape, flag_zero);

    if (flag_internal_poly) {
      if (polymesh.poly_data[kpoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    if (flag_zero) { continue; }

    if (num_Jacobian_shapes == 0) {
      min_Jacobian_shape = Jacobian_shape;
      max_Jacobian_shape = Jacobian_shape;
      poly_with_min_Jacobian_shape = kpoly;
      poly_with_max_Jacobian_shape = kpoly;
    }
    else {

      if (Jacobian_shape < min_Jacobian_shape) {
        min_Jacobian_shape = Jacobian_shape;
        poly_with_min_Jacobian_shape = kpoly;
      }

      if (Jacobian_shape > max_Jacobian_shape) {
        max_Jacobian_shape = Jacobian_shape;
        poly_with_max_Jacobian_shape = kpoly;
      }
    }

    num_Jacobian_shapes++;
  }

}
