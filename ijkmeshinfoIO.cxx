/// \file ijkmeshinfoIO.cxx
/// IO routines for ijkmeshinfo
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


#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ijkcommand_line.txx"
#include "ijkcoord.txx"
#include "ijkIO.txx"
#include "ijkprint.txx"
#include "ijkstring.txx"

#include "ijkmeshinfoIO.h"
#include "ijkmeshinfoIO.txx"
#include "ijkmeshinfo_compute.h"

using namespace std;
using namespace IJKMESHINFO;
using namespace IJK;


// **************************************************
// OUTPUT ROUTINES
// **************************************************

void IJKMESHINFO::output_poly_info
(const int dimension, const POLYMESH_TYPE & polymesh, 
 const COORD_TYPE * vertex_coord, const int poly_index)
{
  cout << "Poly: " << poly_index << endl;

  if (poly_index < 0 || poly_index >= polymesh.NumPoly()) {
    cout << "  Illegal poly index.  Poly index should be in range["
         << 0 << "," << polymesh.NumPoly()-1 << "]." << endl;
    return;
  };

  cout << "  Vertices:" << endl;
  for (int k = 0; k < polymesh.NumPolyVert(poly_index); k++) {
    int iv = polymesh.Vertex(poly_index, k);
    cout << "    " << setw(6) << iv << "  ";
    IJK::print_list(cout, vertex_coord+iv*dimension, dimension);
    cout << endl;
  }

  cout << endl;
}


void IJKMESHINFO::output_degenerate_poly
(const POLYMESH_TYPE & polymesh, const MESH_INFO & mesh_info)
{
  if (mesh_info.num_poly_with_duplicate_vertices > 0) {
    cout << "Degenerate poly:" << endl;
    for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      if (polymesh.poly_data[ipoly].IsDegenerate()) {
        cout << "  Poly " << ipoly << ": ";
        print_list(cout, polymesh.VertexList(ipoly), 
                   polymesh.NumPolyVert(ipoly));
        cout << endl;
      }
    }
  }
}


void IJKMESHINFO::output_duplicate_poly
(const POLYMESH_TYPE & polymesh, const MESH_INFO & mesh_info)
{
  if (mesh_info.num_duplicate_poly > 0) {
    cout << "Duplicate poly:" << endl;
    for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      if (polymesh.poly_data[ipoly].IsDuplicate()) {
        cout << "  Poly " << ipoly << ": ";
        print_list(cout, polymesh.VertexList(ipoly), 
                   polymesh.NumPolyVert(ipoly));
        cout << endl;
      }
    }
  }
}


void IJKMESHINFO::output_duplicate_vertices
(const int dimension, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE vertex_coord[])
{
  for (int jpoly = 0; jpoly < polymesh.NumPoly(); jpoly++) {

    const int * pvert = polymesh.VertexList(jpoly);

    bool flag_dup = false;

    for (int k = 0; k+1 < polymesh.NumPolyVert(jpoly); k++) {
      for (int k2 = k+1; k2 < polymesh.NumPolyVert(jpoly); k2++) {
        if (pvert[k] == pvert[k2]) { flag_dup = true; }
      }
    }

    if (flag_dup) {
      cout << "Duplicate vertices in poly " << jpoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(jpoly));
      cout << endl;
    }

    if (!flag_dup) {

      int i0, i1;
      for (int k = 0; k+1 < polymesh.NumPolyVert(jpoly); k++) {
        for (int k2 = k+1; k2 < polymesh.NumPolyVert(jpoly); k2++) {
          if (IJK::is_coord_equal
              (dimension, vertex_coord+dimension*pvert[k],
               vertex_coord+dimension*pvert[k2])) {
            flag_dup = true;
            i0 = k;
            i1 = k2;
          }
        }
      }

      if (flag_dup) {
        cout << "Vertices with identical coordinates in poly " << jpoly << ": ";
        print_list(cout, pvert, polymesh.NumPolyVert(jpoly));
        cout << endl;
        cout << "  Vertex " << pvert[i0] << ": ";
        print_list(cout, vertex_coord+pvert[i0]*dimension, dimension);
        cout << endl;
        cout << "  Vertex " << pvert[i1] << ": ";
        print_list(cout, vertex_coord+pvert[i1]*dimension, dimension);
        cout << endl;
      }
    }
  }

}


// **************************************************
// OUTPUT POLYGON ANGLE ROUTINES
// **************************************************

void IJKMESHINFO::output_min_max_polygon_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info,
 const bool flag_internal, const int num_poly_edges)
{
  const int dimension = mesh_data.dimension;
  const char * polygon_descriptor;
  IJK::PROCEDURE_ERROR error("output_min_max_polygon_angle");

  string polyname = "polygon";

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  ANGLE_TYPE min_angle, max_angle;
  compute_min_max_polygon_angles
    (mesh_data, polymesh, vertex_coord, flag_internal, num_poly_edges, 
     min_angle, max_angle);

  get_polygon_name(num_poly_edges, polyname);

  if (io_info.flag_output_min_angle) {
    cout << "Min ";
    if (flag_internal) { cout << "internal "; }
    cout << polyname << " angle: ";
    cout << min_angle << endl; 
  }

  if (io_info.flag_output_max_angle) {
    cout << "Max ";
    if (flag_internal) { cout << "internal "; }
    cout << polyname << " angle: ";
    cout << max_angle << endl; 
  }

  if (io_info.angle_le.IsSet() || io_info.angle_ge.IsSet()) {
    int num_le, num_ge;

    if (flag_internal) { polygon_descriptor = "internal polygons"; }
    else { polygon_descriptor = "polygons"; }

    compute_num_polygon_angles
      (dimension, polymesh, vertex_coord, flag_internal, 
       io_info.angle_le.Value(), io_info.angle_ge.Value(), num_le, num_ge);

    if (io_info.angle_le.IsSet() && io_info.flag_output_min_angle) {
      cout << "Number of " << polygon_descriptor << " with angles <= ";
      cout << io_info.angle_le.Value() << ": " << num_le << endl;
    }

    if (io_info.angle_ge.IsSet() && io_info.flag_output_max_angle) {
      cout << "Number of " << polygon_descriptor << " with angles >= ";
      cout << io_info.angle_ge.Value() << ": " << num_ge << endl;
    }
  }

}


// Output polygons with minimum angle
void IJKMESHINFO::output_polygons_with_min_angle
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_polygons_with_min_angle");

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_min_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "polygons", "Poly", "angle", io_info.max_num_poly_out,
     compute_min_max_polygon_angles,
     compute_min_max_polygon_angles);
}


// Output polygons with maximum angle
void IJKMESHINFO::output_polygons_with_max_angle
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_polygons_with_max_angle");

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_max_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "polygons", "Poly", "angle", io_info.max_num_poly_out,
     compute_min_max_polygon_angles,
     compute_min_max_polygon_angles);
}


// Output polygons with minimum and maximum angles
void IJKMESHINFO::output_polygons_with_min_max_angles
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  output_polygons_with_min_angle
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal);
  output_polygons_with_max_angle
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal);
}


void IJKMESHINFO::output_polygons_with_small_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const bool flag_internal, const ANGLE_TYPE angle_bound)
{
  const int dimension = mesh_data.dimension;
  IJK::PROCEDURE_ERROR error("output_polygons_with_small_angles");

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  int num_out = 0;
  COORD_TYPE cos_angle_bound = cos(angle_bound*M_PI/180.0);
  COORD_TYPE cos_min_i, cos_max_i;
  int num_angle;

  if (flag_internal) {
    cout << "Internal polygons with small angles: " << endl;
  }
  else {
    cout << "Polygons with small angles: " << endl;
  }

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].is_degenerate) { continue; }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    IJK::compute_cos_min_max_polygon_angles
      (dimension, polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly),
       vertex_coord, cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0 && cos_min_i >= cos_angle_bound) {
      num_out++;

      const int * pvert = polymesh.VertexList(ipoly);
      ANGLE_TYPE min_angle_i = std::acos(cos_min_i) * 180.0/M_PI;

      cout << "  Poly " << ipoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
      cout << "  Min angle: " << min_angle_i;
      cout << endl;
    }
  }

  if (num_out == 0) {
    if (flag_internal) {
      cout << "  No internal polygons with angles <= " 
           << angle_bound << "." << endl;
    }
    else {
      cout << "  No polygons with angles <= " << angle_bound << "." << endl;
    }
  }
  
}


// Output polygons with large angles
void IJKMESHINFO::output_polygons_with_large_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const bool flag_internal, const ANGLE_TYPE angle_bound)
{
  const int dimension = mesh_data.dimension;
  IJK::PROCEDURE_ERROR error("output_polygons_with_large_angles");

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  int num_out = 0;
  COORD_TYPE cos_angle_bound = cos(angle_bound*M_PI/180.0);
  COORD_TYPE cos_min_i, cos_max_i;
  int num_angle;

  if (flag_internal) {
    cout << "Internal polygons with large angles: " << endl;
  }
  else {
    cout << "Polygons with large angles: " << endl;
  }

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate()) { continue; }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    IJK::compute_cos_min_max_polygon_angles
      (dimension, polymesh.VertexList(ipoly), 
       polymesh.NumPolyVert(ipoly),
       vertex_coord, cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0 && cos_max_i <= cos_angle_bound) {
      num_out++;

      const int * pvert = polymesh.VertexList(ipoly);
      ANGLE_TYPE max_angle_i = std::acos(cos_max_i) * 180.0/M_PI;

      cout << "  Poly " << ipoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
      cout << "  Max angle: " << max_angle_i;
      cout << endl;
    }
  }

  if (num_out == 0) {

    if (flag_internal) {
      cout << "  No internal polygons with angles >= " 
           << angle_bound << "." << endl;
    }
    else {
      cout << "  No polygons with angles >= " << angle_bound << "." << endl;
    }
  }  
}


// **************************************************
// OUTPUT TETRAHEDRA FACET ANGLE ROUTINES
// **************************************************

void IJKMESHINFO::output_min_max_tetrahedra_facet_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_min_max_tetrahedra_facet_angle");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_min_max_poly_values
    (cout, mesh_data, polymesh, vertex_coord,
     io_info.flag_output_min_angle, io_info.flag_output_max_angle,
     flag_internal, "triangle angle", 
     compute_min_max_tetrahedra_facet_angles);
}


void IJKMESHINFO::output_tetrahedra_facet_angle_count
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  const char * tetrahedra_descriptor;
  IJK::PROCEDURE_ERROR error("output_tetrahedra_facet_angle_count");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  if (flag_internal) { tetrahedra_descriptor = "internal tetrahedra"; }
  else { tetrahedra_descriptor = "tetrahedra"; }

  int num_le, num_ge;
  if (io_info.facet_angle_le.IsSet() || io_info.facet_angle_ge.IsSet()) {
    compute_num_tetrahedra_facet_angles
      (dimension, mesh_dimension, polymesh, vertex_coord, flag_internal, 
       io_info.facet_angle_le.Value(), io_info.facet_angle_ge.Value(), 
       num_le, num_ge);

    if (io_info.facet_angle_le.IsSet() && io_info.flag_output_min_angle) {
      cout << "Number of " << tetrahedra_descriptor
           << " with facet angles <= ";
      cout << io_info.facet_angle_le.Value() << ": " << num_le << endl;
    }

    if (io_info.facet_angle_ge.IsSet() && io_info.flag_output_max_angle) {
      cout << "Number of " << tetrahedra_descriptor
           << " with facet angles >= ";
      cout << io_info.facet_angle_ge.Value() << ": " << num_ge << endl;
    }
  }

}


// Output tetrahedra with minimum facet angle
void IJKMESHINFO::output_tetrahedra_with_min_facet_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_min_facet_angle");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_min_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "tetrahedra", "Tet", "triangle angle", io_info.max_num_poly_out,
     compute_min_max_tetrahedra_facet_angles,
     compute_min_max_tetrahedron_facet_angles);
}


// Output tetrahedra with maximum facet angle
void IJKMESHINFO::output_tetrahedra_with_max_facet_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_max_facet_angle");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_max_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "tetrahedra", "Tet", "triangle angle", io_info.max_num_poly_out,
     compute_min_max_tetrahedra_facet_angles,
     compute_min_max_tetrahedron_facet_angles);
}


// Output tetrahedra with minimum and maximum facet angles
void IJKMESHINFO::output_tetrahedra_with_min_max_facet_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  output_tetrahedra_with_min_facet_angle
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal);
  output_tetrahedra_with_max_facet_angle
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal);
}


void IJKMESHINFO::output_tetrahedra_with_small_facet_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const bool flag_internal, const ANGLE_TYPE angle_bound)
{
  const int dimension = mesh_data.dimension;
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_small_facet_angles");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  int num_out = 0;
  COORD_TYPE cos_angle_bound = cos(angle_bound*M_PI/180.0);
  COORD_TYPE cos_min_i, cos_max_i;
  int num_angle;

  if (flag_internal) {
    cout << "Internal tetrahedra with small facet angles: " << endl;
  }
  else {
    cout << "Tetrahedra with small facet angles: " << endl;
  }

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].is_degenerate) { continue; }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    compute_cos_min_max_tetrahedron_facet_angles
      (dimension, polymesh.VertexList(ipoly),
       vertex_coord, cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0 && cos_min_i >= cos_angle_bound) {
      num_out++;

      const int * pvert = polymesh.VertexList(ipoly);
      ANGLE_TYPE min_angle_i = std::acos(cos_min_i) * 180.0/M_PI;

      cout << "  Tet " << ipoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
      cout << "  Min facet angle: " << min_angle_i;
      cout << endl;
    }
  }

  if (num_out == 0) {
    if (flag_internal) {
      cout << "  No internal tetrahedra with facet angles <= " 
           << angle_bound << "." << endl;
    }
    else {
      cout << "  No tetrahedra with facet angles <= " 
           << angle_bound << "." << endl;
    }
  }
  
}


// Output tetrahedra with large facet angles.
void IJKMESHINFO::output_tetrahedra_with_large_facet_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const bool flag_internal, const ANGLE_TYPE angle_bound)
{
  const int dimension = mesh_data.dimension;
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_large_facet_angles");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  int num_out = 0;
  COORD_TYPE cos_angle_bound = cos(angle_bound*M_PI/180.0);
  COORD_TYPE cos_min_i, cos_max_i;
  int num_angle;

  if (flag_internal) {
    cout << "Internal tetrahedra with large facet angles: " << endl;
  }
  else {
    cout << "Tetrahedra with large facet angles: " << endl;
  }

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate()) { continue; }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    compute_cos_min_max_tetrahedron_facet_angles
      (dimension, polymesh.VertexList(ipoly), 
       vertex_coord, cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0 && cos_max_i <= cos_angle_bound) {
      num_out++;

      const int * pvert = polymesh.VertexList(ipoly);
      ANGLE_TYPE max_angle_i = std::acos(cos_max_i) * 180.0/M_PI;

      cout << "  Tet " << ipoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
      cout << "  Max facet angle: " << max_angle_i;
      cout << endl;
    }
  }

  if (num_out == 0) {

    if (flag_internal) {
      cout << "  No internal tetrahedra with facet angles >= " 
           << angle_bound << "." << endl;
    }
    else {
      cout << "  No tetrahedra with facet angles >= " 
           << angle_bound << "." << endl;
    }
  }  
}


// **************************************************
// OUTPUT DIHEDRAL ANGLE ROUTINES
// **************************************************

void IJKMESHINFO::output_min_max_dihedral_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_min_max_dihedral_angle");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_min_max_poly_values
    (cout, mesh_data, polymesh, vertex_coord,
     io_info.flag_output_min_angle, io_info.flag_output_max_angle,
     flag_internal, "dihedral angle", compute_min_max_tetmesh_dihedral_angles);
}


void IJKMESHINFO::output_dihedral_angle_count
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  const char * tetrahedra_descriptor;
  IJK::PROCEDURE_ERROR error("output_dihedral_angle_count");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  if (flag_internal) { tetrahedra_descriptor = "internal tetrahedra"; }
  else { tetrahedra_descriptor = "tetrahedra"; }

  int num_le, num_ge;
  if (io_info.angle_le.IsSet() || io_info.angle_ge.IsSet()) {
    compute_num_tetmesh_dihedral_angles
      (dimension, mesh_dimension, polymesh, vertex_coord, flag_internal, 
       io_info.angle_le.Value(), io_info.angle_ge.Value(), num_le, num_ge);

    if (io_info.angle_le.IsSet() && io_info.flag_output_min_angle) {
      cout << "Number of " << tetrahedra_descriptor
           << " with dihedral angles <= ";
      cout << io_info.angle_le.Value() << ": " << num_le << endl;
    }

    if (io_info.angle_ge.IsSet() && io_info.flag_output_max_angle) {
      cout << "Number of " << tetrahedra_descriptor
           << " with dihedral angles >= ";
      cout << io_info.angle_ge.Value() << ": " << num_ge << endl;
    }
  }

}


// Output tetrahedra with minimum dihedral angle.
void IJKMESHINFO::output_tetrahedra_with_min_dihedral_angle
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_min_dihedral_angle");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_min_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "tetrahedra", "Tet", "dihedral angle", io_info.max_num_poly_out,
     compute_min_max_tetmesh_dihedral_angles,
     compute_min_max_tetrahedron_dihedral_angles);
}


// Output tetrahedra with maximum dihedral angle.
void IJKMESHINFO::output_tetrahedra_with_max_dihedral_angle
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_max_dihedral_angle");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_max_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "tetrahedra", "Tet", "dihedral angle", io_info.max_num_poly_out,
     compute_min_max_tetmesh_dihedral_angles,
     compute_min_max_tetrahedron_dihedral_angles);

}


// Output tetrahedra with minimum and maximum dihedral angles
void IJKMESHINFO::output_tetrahedra_with_min_max_dihedral_angles
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  output_tetrahedra_with_min_dihedral_angle
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal);
  output_tetrahedra_with_max_dihedral_angle
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal);
}


void IJKMESHINFO::output_tetrahedra_with_small_dihedral_angles
(const MESH_DATA & mesh_data, 
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const bool flag_internal, const ANGLE_TYPE angle_bound)
{
  const int dimension = mesh_data.dimension;
  IJK::PROCEDURE_ERROR error("output_polygons_with_small_angles");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  int num_out = 0;
  COORD_TYPE cos_angle_bound = cos(angle_bound*M_PI/180.0);
  COORD_TYPE cos_min_i, cos_max_i;
  int num_angle;

  if (flag_internal) 
    { cout << "Internal tetrahedra with small dihedral angles: " << endl; }
  else 
    { cout << "Tetrahedra  with small dihedral angles: " << endl; }

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].is_degenerate) { continue; }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    compute_cos_min_max_tetrahedron_dihedral_angles
      (dimension, polymesh.VertexList(ipoly), vertex_coord,
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0 && cos_min_i >= cos_angle_bound) {
      num_out++;

      const int * pvert = polymesh.VertexList(ipoly);
      ANGLE_TYPE min_angle_i = std::acos(cos_min_i) * 180.0/M_PI;

      cout << "  Tet " << ipoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
      cout << "  Min dihedral angle: " << min_angle_i;
      cout << endl;
    }
  }

  if (num_out == 0) {
    if (flag_internal) {
      cout << "  No internal tetrahedra with angles <= " 
           << angle_bound << "." << endl;
    }
    else {
      cout << "  No tetrahedra with angles <= " << angle_bound << "." << endl;
    }
  }
  
}


// Output tetrahedra with large dihedral angles
void IJKMESHINFO::output_tetrahedra_with_large_dihedral_angles
(const MESH_DATA & mesh_data, 
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const bool flag_internal, const ANGLE_TYPE angle_bound)
{
  const int dimension = mesh_data.dimension;
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_large_dihedral_angles");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  int num_out = 0;
  COORD_TYPE cos_angle_bound = cos(angle_bound*M_PI/180.0);
  COORD_TYPE cos_min_i, cos_max_i;
  int num_angle;

  if (flag_internal) 
    { cout << "Internal tetrahedra with large dihedral angles: " << endl; }
  else 
    { cout << "Tetrahedra with large angles: " << endl; }

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate()) { continue; }

    if (flag_internal) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    compute_cos_min_max_tetrahedron_dihedral_angles
      (dimension, polymesh.VertexList(ipoly), vertex_coord,
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0 && cos_max_i <= cos_angle_bound) {
      num_out++;

      const int * pvert = polymesh.VertexList(ipoly);
      ANGLE_TYPE max_angle_i = std::acos(cos_max_i) * 180.0/M_PI;

      cout << "  Poly " << ipoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
      cout << "  Max angle: " << max_angle_i;
      cout << endl;
    }
  }

  if (num_out == 0) {
    if (flag_internal) {
      cout << "  No internal tetrahedra with angles >= " 
           << angle_bound << "." << endl;
    }
    else {
      cout << "  No tetrahedra with angles >= " << angle_bound << "." << endl;
    }
  }  
}


// **************************************************
// OUTPUT EDGE LENGTH ROUTINES
// **************************************************

void IJKMESHINFO::output_min_max_polygon_edge_lengths
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal,
 const int num_poly_vert)
{
  const int dimension = mesh_data.dimension;
  const char * polygon_descriptor;
  IJK::PROCEDURE_ERROR error("output_min_max_polygon_edge_lengths");

  string polyname = "polygon";

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  COORD_TYPE min_edge_length, max_edge_length;

  compute_min_max_polygon_edge_lengths_select_poly_by_numv
    (mesh_data, polymesh, vertex_coord, flag_internal, num_poly_vert, 
     min_edge_length, max_edge_length);


  get_polygon_name(num_poly_vert, polyname);

  if (io_info.flag_output_min_edge_length) {
    cout << "Min ";
    if (flag_internal) { cout << "internal "; }
    cout << polyname << " edge length: ";
    cout << min_edge_length << endl; 
  }

  if (io_info.flag_output_max_edge_length) {
    cout << "Max ";
    if (flag_internal) { cout << "internal "; }
    cout << polyname << " edge length: ";
    cout << max_edge_length << endl; 
  }
}


void IJKMESHINFO::output_min_max_tetrahedra_edge_lengths
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  IJK::PROCEDURE_ERROR error("output_min_max_tetrahedra_edge_lengths");

  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_min_max_poly_values
    (cout, mesh_data, polymesh, vertex_coord,
     io_info.flag_output_min_edge_length, io_info.flag_output_max_edge_length,
     flag_internal, "edge length", compute_min_max_tetrahedra_edge_lengths);
}


void IJKMESHINFO::output_min_max_simplices_edge_lengths
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  IJK::PROCEDURE_ERROR error("output_min_max_simplices_edge_lengths");

  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_min_max_poly_values
    (cout, mesh_data, polymesh, vertex_coord,
     io_info.flag_output_min_edge_length, io_info.flag_output_max_edge_length,
     flag_internal, "edge length", compute_min_max_simplices_edge_lengths);
}


void IJKMESHINFO::output_min_max_hexahedra_edge_lengths
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  IJK::PROCEDURE_ERROR error("output_min_max_hexahedra_edge_lengths");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_min_max_poly_values
    (cout, mesh_data, polymesh, vertex_coord,
     io_info.flag_output_min_edge_length, io_info.flag_output_max_edge_length,
     flag_internal, "edge length", compute_min_max_hexahedra_edge_lengths);
}


// Output polygons with min edge lengths.
void IJKMESHINFO::output_polygons_with_min_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_polygons_with_min_edge_lengths");

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_min_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "polgyons", "Poly", "edge length", io_info.max_num_poly_out,
     compute_min_max_polygon_edge_lengths,
     compute_min_max_polygon_edge_lengths);
}


// Output polygons with max edge lengths.
void IJKMESHINFO::output_polygons_with_max_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_polygons_with_max_edge_lengths");

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_max_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "polgyons", "Poly", "edge length", io_info.max_num_poly_out,
     compute_min_max_polygon_edge_lengths,
     compute_min_max_polygon_edge_lengths);
}


// Output polygons with min and maximum edge lengths.
void IJKMESHINFO::output_polygons_with_min_max_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  output_polygons_with_min_edge_lengths
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal);
  output_polygons_with_max_edge_lengths
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal);
}


// Output tetrahedra with min edge lengths.
void IJKMESHINFO::output_tetrahedra_with_min_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_min_edge_lengths");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_min_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "tetrahedra", "Tet", "edge length", io_info.max_num_poly_out,
     compute_min_max_tetrahedra_edge_lengths,
     compute_min_max_tetrahedron_edge_lengths);
}


// Output tetrahedra with max edge lengths.
void IJKMESHINFO::output_tetrahedra_with_max_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_tetrahedra_with_max_edge_lengths");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_max_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "tetrahedra", "Tet", "edge length", io_info.max_num_poly_out,
     compute_min_max_tetrahedra_edge_lengths,
     compute_min_max_tetrahedron_edge_lengths);
}


// Output hexahedra with min edge lengths.
void IJKMESHINFO::output_hexahedra_with_min_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_hexahedra_with_min_edge_lengths");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_min_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "hexahedra", "Hex", "edge length", io_info.max_num_poly_out,
     compute_min_max_hexahedra_edge_lengths,
     compute_min_max_hexahedron_edge_lengths);
}


// Output hexahedra with max edge lengths.
void IJKMESHINFO::output_hexahedra_with_max_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error("output_hexahedra_with_max_edge_lengths");

  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_max_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "hexahedra", "Hex", "edge length", io_info.max_num_poly_out,
     compute_min_max_hexahedra_edge_lengths,
     compute_min_max_hexahedron_edge_lengths);
}


// **************************************************
// OUTPUT JACOBIAN ROUTINES
// **************************************************

void IJKMESHINFO::output_min_max_hexahedra_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, 
 COORD_TYPE & max_Jacobian_determinant)
{
  IJK::PROCEDURE_ERROR error("output_min_max_hexahedra_Jacobian_determinants");

  if (!check_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_min_max_poly_values
    (cout, mesh_data, polymesh, vertex_coord, 
     io_info.flag_output_min_Jacobian_determinant,
     io_info.flag_output_max_Jacobian_determinant,
     flag_internal, "Jacobian determinant", 
     compute_min_max_hexahedra_Jacobian_determinants,
     min_Jacobian_determinant, max_Jacobian_determinant);
}


void IJKMESHINFO::output_min_max_hexahedra_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  COORD_TYPE min_Jacobian_determinant, max_Jacobian_determinant;

  output_min_max_hexahedra_Jacobian_determinants
    (mesh_data, polymesh, vertex_coord, io_info,
     flag_internal, min_Jacobian_determinant, max_Jacobian_determinant);
}


void IJKMESHINFO::output_min_max_hexahedra_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, 
 COORD_TYPE & max_Jacobian_determinant)
{
  IJK::PROCEDURE_ERROR error
    ("output_min_max_hexahedra_normalized_Jacobian_determinants");

  if (!check_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_min_max_poly_values
    (cout, mesh_data, polymesh, vertex_coord, 
     io_info.flag_output_min_normalized_Jacobian_determinant,
     io_info.flag_output_max_normalized_Jacobian_determinant,
     flag_internal, "normalized Jacobian determinant", 
     compute_min_max_hexahedra_normalized_Jacobian_determinants,
     min_Jacobian_determinant, max_Jacobian_determinant);
}


void IJKMESHINFO::output_min_max_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, 
 const bool flag_internal_poly,
 const bool flag_internal_vertex,
 COORD_TYPE & min_Jacobian_determinant, 
 COORD_TYPE & max_Jacobian_determinant)
{
  IJK::PROCEDURE_ERROR error("output_min_max_hex_vert_Jacobian_determinants");

  if (!check_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal_poly || flag_internal_vertex) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  if (flag_internal_vertex) {
    output_min_max_vertex_values
      (cout, mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info.flag_output_min_Jacobian_determinant,
       io_info.flag_output_max_Jacobian_determinant,
       flag_internal_poly, flag_internal_vertex,
       "Jacobian determinant (at internal vert)", 
       compute_min_max_hex_vert_Jacobian_determinants,
       min_Jacobian_determinant, max_Jacobian_determinant);
  }
  else {
    output_min_max_vertex_values
      (cout, mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info.flag_output_min_Jacobian_determinant,
       io_info.flag_output_max_Jacobian_determinant,
       flag_internal_poly, flag_internal_vertex,
       "Jacobian determinant (at vert)", 
       compute_min_max_hex_vert_Jacobian_determinants,
       min_Jacobian_determinant, max_Jacobian_determinant);
  }
}


void IJKMESHINFO::output_min_max_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, 
 const bool flag_internal_poly,
 const bool flag_internal_vertex,
 COORD_TYPE & min_Jacobian_determinant, 
 COORD_TYPE & max_Jacobian_determinant)
{
  IJK::PROCEDURE_ERROR error
    ("output_min_max_hex_vert_normalized_Jacobian_determinants");

  if (!check_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal_poly || flag_internal_vertex) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  if (flag_internal_vertex) {
    output_min_max_vertex_values
      (cout, mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info.flag_output_min_Jacobian_determinant,
       io_info.flag_output_max_Jacobian_determinant,
       flag_internal_poly, flag_internal_vertex,
       "normalized Jacobian determinant (at internal vert)", 
       compute_min_max_hex_vert_normalized_Jacobian_determinants,
       min_Jacobian_determinant, max_Jacobian_determinant);
  }
  else {
    output_min_max_vertex_values
      (cout, mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info.flag_output_min_Jacobian_determinant,
       io_info.flag_output_max_Jacobian_determinant,
       flag_internal_poly, flag_internal_vertex,
       "normalized Jacobian determinant (at vert)", 
       compute_min_max_hex_vert_normalized_Jacobian_determinants,
       min_Jacobian_determinant, max_Jacobian_determinant);
  }
}


// Output hexahedra with min Jacobian determinants.
void IJKMESHINFO::output_hexahedra_with_min_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error
    ("output_hexahedra_with_min_Jacobian_determinants");

  if (!check_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_min_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "hexahedra", "Hex", "Jacobian determinant", io_info.max_num_poly_out,
     compute_min_max_hexahedra_Jacobian_determinants,
     compute_min_max_hexahedron_Jacobian_determinants);
}


// Output hexahedra with max Jacobian determinants.
void IJKMESHINFO::output_hexahedra_with_max_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, const bool flag_internal)
{
  IJK::PROCEDURE_ERROR error
    ("output_hexahedra_with_max_Jacobian_determinants");

  if (!check_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  output_polytopes_with_max_value
    (cout, mesh_data, polymesh, vertex_coord, flag_internal,
     "hexahedra", "Hex", "Jacobian determinant", io_info.max_num_poly_out,
     compute_min_max_hexahedra_Jacobian_determinants,
     compute_min_max_hexahedron_Jacobian_determinants);
}


// Output hex mesh vertices with min Jacobian determinants.
void IJKMESHINFO::output_hex_mesh_vertices_with_min_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, 
 const bool flag_internal_poly,
 const bool flag_internal_vertex)
{
  COORD_TYPE min_Jacobian_determinant, max_Jacobian_determinant;
  int poly_with_min_val, poly_with_max_val;
  VERTEX_INDEX vert_with_min_val, vert_with_max_val;
  const CUBE_TYPE cube(DIM3);
  IJK::PROCEDURE_ERROR error
    ("output_hex_mesh_vertices_with_min_Jacobian_determinants");

  if (!check_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal_poly || flag_internal_vertex) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  compute_min_max_hex_vert_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vertex,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_val, poly_with_max_val, 
     vert_with_min_val, vert_with_max_val);

  if (flag_internal_vertex) {
    cout << "Internal vertices with minimum Jacobian determinant " 
         << min_Jacobian_determinant <<  ":" << endl;
  }
  else {
    cout << "Vertices with minimum Jacobian determinant " 
         << min_Jacobian_determinant <<  ":" << endl;
  }

  for (VERTEX_INDEX iv = 0; iv < vertex_poly_incidence.NumVertices(); iv++) {
    bool flag_minv = false;

    if (flag_internal_vertex) {

      if (polymesh.vertex_data[iv].OnBoundary()) {
        // Vertex iv is not internal.
        continue;
      }
    }

    for (int j = 0; j < vertex_poly_incidence.NumIncidentPoly(iv); j++) {
      const int jpoly = vertex_poly_incidence.IncidentPoly(iv, j);
      
      if (flag_internal_poly) {
        if (polymesh.poly_data[jpoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      COORD_TYPE det;
      int kloc;
      if (polymesh.DoesPolyContainVertex(jpoly,iv,kloc)) {
        IJK::compute_hexahedron_Jacobian_determinant_3D
          (polymesh.VertexList(jpoly), mesh_data.orientation,
           vertex_coord, cube, kloc, det);
      }
      else {
        error.AddMessage
          ("Programming error.  Inconsistency between polymesh and vertex adjacency list.");
        error.AddMessage
          ("  Vertex adjacency list for vertex ", iv,
           " contains poly ", jpoly, ".");
        error.AddMessage
          ("  Polymesh poly ", jpoly, " does not have vertex ", iv, ".");
        throw error;
      }

      if ((iv == vert_with_min_val && jpoly == poly_with_min_val) ||
          (det <= min_Jacobian_determinant)) {
        if (flag_minv) { cout << "  " << jpoly; }
        else {
          cout << "  Vertex " << iv << ".  Min Jacobian det in polytopes: "
               << jpoly;
        }
        flag_minv = true;
      }
    }

    if (flag_minv) { cout << endl; }
  }

}


// Output hex mesh vertices with min normalized Jacobian determinants.
void IJKMESHINFO::
output_hex_mesh_vertices_with_min_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,
 const IO_INFO & io_info, 
 const bool flag_internal_poly,
 const bool flag_internal_vertex)
{
  COORD_TYPE min_Jacobian_determinant, max_Jacobian_determinant;
  int poly_with_min_val, poly_with_max_val;
  VERTEX_INDEX vert_with_min_val, vert_with_max_val;
  const CUBE_TYPE cube(DIM3);
  const COORD_TYPE max_small_magnitude(0.0);
  bool flag_zero;
  IJK::PROCEDURE_ERROR error
    ("output_hex_mesh_vertices_with_min_Jacobian_determinants");

  if (!check_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (!check_mesh_dimension<DIM3>(mesh_data, error)) { throw error; }
  if (flag_internal_poly || flag_internal_vertex) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  compute_min_max_hex_vert_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vertex,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_val, poly_with_max_val, 
     vert_with_min_val, vert_with_max_val);

  if (flag_internal_vertex) {
    cout << "Internal vertices with minimum normalized Jacobian determinant " 
         << min_Jacobian_determinant <<  ":" << endl;
  }
  else {
    cout << "Vertices with minimum normalized Jacobian determinant " 
         << min_Jacobian_determinant <<  ":" << endl;
  }

  for (VERTEX_INDEX iv = 0; iv < vertex_poly_incidence.NumVertices(); iv++) {
    bool flag_minv = false;

    if (flag_internal_vertex) {

      if (polymesh.vertex_data[iv].OnBoundary()) {
        // Vertex iv is not internal.
        continue;
      }
    }

    for (int j = 0; j < vertex_poly_incidence.NumIncidentPoly(iv); j++) {
      const int jpoly = vertex_poly_incidence.IncidentPoly(iv, j);
      
      if (flag_internal_poly) {
        if (polymesh.poly_data[jpoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      COORD_TYPE det;
      int kloc;
      if (polymesh.DoesPolyContainVertex(jpoly,iv,kloc)) {
        IJK::compute_hexahedron_normalized_Jacobian_determinant_3D
          (polymesh.VertexList(jpoly), mesh_data.orientation,
           vertex_coord, cube, kloc, max_small_magnitude, 
           det, flag_zero);
      }
      else {
        error.AddMessage
          ("Programming error.  Inconsistency between polymesh and vertex adjacency list.");
        error.AddMessage
          ("  Vertex adjacency list for vertex ", iv,
           " contains poly ", jpoly, ".");
        error.AddMessage
          ("  Polymesh poly ", jpoly, " does not have vertex ", iv, ".");
        throw error;
      }

      if ((iv == vert_with_min_val && jpoly == poly_with_min_val) ||
          (det <= min_Jacobian_determinant && !flag_zero)) {
        if (flag_minv) { cout << "  " << jpoly; }
        else {
          cout << "  Vertex " << iv 
               << ".  Min normalized Jacobian det in polytopes: " << jpoly;
        }
        flag_minv = true;
      }
    }

    if (flag_minv) { cout << endl; }
  }

}


// **************************************************
// POLYGON/POLYTOPE NAME ROUTINES
// **************************************************

void IJKMESHINFO::get_polygon_name
(const int num_polygon_edges, std::string & polygon_name)
{
  switch(num_polygon_edges) {

  case 0:
    polygon_name = "polygon";
    break;

  case 3:
    polygon_name = "triangle";
    break;

  case 4:
    polygon_name = "quadrilateral";
    break;

  case 5:
    polygon_name = "pentagon";
    break;

  case 6:
    polygon_name = "hexagon";
    break;

  default:
    {
      string s;
      val2string(num_polygon_edges, s);
      polygon_name = string("polygon (") + s + " edges)"; 
      break;
    }
  }
}


// **************************************************
// PRINT ROUTINES
// **************************************************

void IJKMESHINFO::print_poly_index_and_vert
(std::ostream & out, 
 const char * prefix, const char * separator, const char * suffix,
 const POLYMESH_TYPE & polymesh, const int ipoly)
{
  const VERTEX_INDEX * pvert = polymesh.VertexList(ipoly);

  out << prefix << ipoly << separator;
  print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
  out << suffix;
}



// **************************************************
// USAGE/HELP MESSAGES
// **************************************************

void IJKMESHINFO::usage_msg()
{
  cerr << "Usage: ijkmeshinfo [OPTIONS] {input file}" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  [-mesh_dim {mdim}] [-reverse_orient]" << endl;
  cerr << "  [-vertex {vnum}] [-simplex {snum}] [-poly {pnum}]"
       << endl;
  cerr << "  [-vlist | -vlist_min] [-plist]" << endl;
  cerr << "  [-manifold | -oriented_manifold | -check_facetI]" << endl;
  cerr << "  [-containsv {vnum}] [-containse {end0} {end1}]" << endl;
  cerr << "  [-minc \"min coord\"] [-maxc \"max coord\"]" << endl;
  cerr << "  [-min_numv <N>] [-max_numv <N>]" << endl;
  cerr << "  [-angle_le <A>] [-angle_ge <A>]" << endl;
  cerr << "  [-facet_angle_le <A>] [-facet_angle_ge <A>]" << endl;
  cerr << "  [-list_dup] [-internal] [-internal_vert]" << endl;
  cerr << "  [-selfI]" << endl;
  cerr << "  [-out_values] [-out_min_angle] [-out_max_angle]"
       << endl;
  cerr << "  [-out_values] [-out_min_jdet]"
       << endl;
  cerr << "  [-plot_angles] [-plot_edge_lengths] [-plot_jacobian]"
       << endl;
  cerr << "  [-report_deep] [-for_each_type]" << endl;
  cerr << "  [-max_out <N>] [-terse] [-help]" << endl;
}

void IJKMESHINFO::usage_error()
{
  usage_msg();
  exit(10);
}

void IJKMESHINFO::help_msg()
{
  cerr << "Usage: ijkmeshinfo [OPTIONS] {input file}" << endl;
  cerr << "  -mesh_dim {mdim}:   Mesh dimension." << endl;
  cerr << "  -reverse_orient:    Reverse the orientation of each mesh element." << endl;
  cerr << "  -vertex {vnum}:     Print coordinates of vertex {vnum}."
       << endl;
  cerr << "  -simplex {snum}:    Print vertices of simplex {snum}."
       << endl;
  cerr << "  -poly {pnum}:       Print vertices of poly {pnum}."
       << endl;
  cerr << "  -vlist:             Print list of vertex coordinates." << endl;
  cerr << "     Prints only vertices with coordinates between min coord and max coord."
       << endl;
  cerr << "  -vlist_min:         Print vertices with min values." << endl;
  cerr << "     (Only implemented for hexahedral meshes.)" << endl;
  cerr << "  -plist:             Print list of polygons/polytopes." << endl;
  cerr << "     Print polygons with min/max angles, edge lengths, and"
       << endl
       << "       Jacobian determinants." << endl;
  cerr << "     If -min_numv or -max_numv is set, prints only "
       << endl
       << "       polygons/polytopes with number of vertices between"
       << endl;
  cerr << "       <min_numv> and <max_numv>." << endl;
  cerr << "  -manifold:          Print manifold information." << endl;
  cerr << "  -oriented_manifold: Print manifold and orientation information."
       << endl;
  cerr << "  -check_facetI:      Check facet intersections." << endl
       << "     Check for hexahedra facets which intersect in exactly two edges." << endl
       << "     Options -manifold and -oriented_manifold always perform this check." << endl;
  cerr << "  -selfI:             Print self intersections." << endl;
  cerr << "  -containsv {vnum}:  Print list of simplices containing vertex {vnum}."
       << endl;
  cerr << "  -containse {end0} {end1}}:  Print list of simplices containing"
       << endl;
  cerr << "                          edge ({end0},{end1})." << endl;
  cerr << "  -minc \"min coord\":  Minimum coordinates." << endl;
  cerr << "  -maxc \"min coord\":  Maximum coordinates." << endl;
  cerr << "  -min_numv <N>:      Minimum number of poly vertices." << endl;
  cerr << "  -max_numv <N>:      Maximum number of poly vertices." << endl;
  cerr << "  -angle_le <A>:" << endl;
  cerr << "     Report number of polygons with angles less than or equal to <A>"
       << endl
       << "       or number of tetrahedra with dihedral angles less than" 
       << endl
       << "       or equal to <A>." << endl;
  cerr << "     Use with -plist to list polygons or tetrahedra with angles"
       << endl
       << "       or dihedral angles less than or equal to <A>." << endl;
  cerr << "  -angle_ge <A>:" << endl;
  cerr << "     Report number of polygons with angles greater than or equal to <A>" << endl
       << "       or number of tetrahedra with dihedral angles greater than" 
       << endl
       << "       or equal to <A>." << endl;
  cerr << "     Use with -plist to list polygons or tetrahedra with angles or"
       << endl
       << "       dihedral angles greater than or equal to <A>." << endl;
  cerr << "  -facet_angle_le <A>:" << endl;
  cerr << "     Report number of tetrahedra with facet angles less than or equal to <A>." << endl;
  cerr << "     Use with -plist to list tetrahedra with facet angles less than"
       << endl
       << "       or equal to <A>." << endl;
  cerr << "  -facet_angle_ge <A>:" << endl;
  cerr << "     Report number of tetrahedra with facet angles greater than or equal to <A>." << endl;
  cerr << "     Use with -plist to list tetrahedra with facet angles greater than" << endl
       << "       or equal to <A>." << endl;
  cerr << "  -list_dup:          List duplicate poly vertices or poly vertices" << endl;
  cerr << "                           with identical coordinates." << endl;
  cerr << "  -internal:          Report only min/max values for internal polytopes."
       << endl;
  cerr << "  -internal_vert:     Report only min/max values at internal mesh vertices."
       << endl;
  cerr << "  -out_values:        Output only values for -manifold option."
       << endl;
  cerr << "     Output number of non-manifold vertices, number of non-manifold vertices"
       << endl
       << "     at least 1 from bounding box, number of non-manifold edges," 
       << endl
       << "     number of internal boundary facets and number of facets"
       << endl
       << "     at least 1 from bounding box." << endl;
  cerr << "  -out_min_angle:     Output minimum triangle angle." << endl;
  cerr << "  -out_max_angle:     Output maximum triangle angle." << endl;
  cerr << "  -out_min_jdet:      Output minimum Jacobian determinant." << endl
       << "    (Use with -vlist to list vertices with minimum Jacobian" 
       << endl
       << "     determinant of vertices.  Note that the minimum Jacobian "
       << endl
       << "     determinant of vertices may be greater than the minimum"
       << endl
       << "     Jacobian determinant.)" << endl;
  cerr << "  -out_min_jdet:      Output maximum Jacobian determinant." 
       << endl;
  cerr << "  -plot_angles:       Create gnuplot (.gplt) files of min and max"
       << endl
       << "                          triangle angles." << endl;
  cerr << "  -plot_edge_lengths: Create gnuplot (.gplt) files of min and max"
       << endl
       << "                          polytope edge lengths." << endl;
  cerr << "  -plot_jacobian:     Create gnuplot (.gplt) files of min and max"
       << endl
       << "                          hexahedra Jacobian determinants." 
       << endl;
  cerr << "  -report_deep:       Report only boundary facets at least distance 1"
       << endl
       << "                          from bounding box." << endl;
  cerr << "  -for_each_type:     Reports angles separately for triangles, "
       << endl
       << "                          quadrilaterals, pentagons, ..." << endl;
  cerr << "  -max_out <N>:       Limit number of polytopes in output lists to <N>." << endl;
  cerr << "  -terse:             Terse output.  Affects -manifold and -containsv options." << endl;
  cerr << "  -help:              Print this help message." << endl;
  exit(20);
}


// **************************************************
// Class IO_INFO member functions
// **************************************************

void IO_INFO::Init()
{
  max_num_poly_out = 5;
  flag_output_min_angle = false;
  flag_output_max_angle = false;
  flag_output_min_edge_length = false;
  flag_output_max_edge_length = false;
  flag_output_all_min_max = false;
  flag_output_min_Jacobian_determinant = false;
  flag_output_max_Jacobian_determinant = false;
  flag_output_min_normalized_Jacobian_determinant = false;
  flag_output_max_normalized_Jacobian_determinant = false;
  flag_general_info = true;
};

