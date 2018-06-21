/// \file ijkmeshinfo.cxx
/// compute mesh information
/// Version 0.2.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2008-2018 Rephael Wenger

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


#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <algorithm>
#include <map>

#include "ijkcommand_line.txx"
#include "ijkcoord.txx"
#include "ijkIO.txx"
#include "ijkmerge.txx"
#include "ijkmesh.txx"
#include "ijkmesh_faces.txx"
#include "ijkprint.txx"
#include "ijkstring.txx"

#include "ijkmeshinfo.h"
#include "ijkmeshinfo_compute.h"
#include "ijkmeshinfoIO.h"

using namespace std;
using namespace IJKMESHINFO;
using namespace IJK;

// global variables
bool is_mesh_dimension_set(false);
int num_simplices = 0;
int num_edges = 0;
int num_poly = 0;
COORD_TYPE * vertex_coord = NULL;
int * simplex_vert = NULL;
int num_vert_per_simplex;
POLYMESH_TYPE polymesh, polymesh_sorted;
int num_vert_per_poly(0);  // Set only if all poly have same number of vertices.
bool standard_input = false;
char * input_filename = NULL;
char * output_filename = NULL;
BOUNDING_BOX bounding_box(3);
BOUNDING_BOX contracted_bounding_box(3);
COORD_TYPE contract_margin = 1.0;
bool flag_small_bounding_box = false;
bool is_min_coord_set = false;
bool is_max_coord_set = false;
int min_num_polyv_output = 0;
int max_num_polyv_output = 0;
bool is_min_num_polyv_output_set = false;
bool is_max_num_polyv_output_set = false;
bool flag_simplex_file(false);
bool flag_cube_file(false);
bool flag_polyfile(false);
std::vector<COORD_TYPE> min_coord;
std::vector<COORD_TYPE> max_coord;
bool flag_list_duplicate_vertices = false;
bool flag_list_duplicate_poly = false;
bool flag_report_deep = false;         // Report only facets "deep" in bounding box.
bool flag_output_only_values = false;  // Output only values without any text.
bool flag_internal_poly = false;       // Output properties of internal polytopes.
bool flag_internal_vert = false;       // Output properties of internal vertices.
bool flag_for_each_type;               // Report min/max for tri, quad, etc.
bool flag_report_self_intersections = false;  // Report self intersections.
double selfI_epsilon = 1.0e-10;
bool flag_use_grid_of_bins = true;
int num_bins_per_axis = 10;
bool is_num_bins_per_axis_set = false;
bool flag_plot_angles = false;               // Plot some angles.
bool flag_plot_min_polygon_angles = false;   // Plot min polygon angles.
bool flag_plot_max_polygon_angles = false;   // Plot max polygon angles.
bool flag_plot_min_edge_lengths = false;     // Plot min edges lengths.
bool flag_plot_max_edge_lengths = false;     // Plot min edges lengths.
bool flag_plot_min_jacobian = false;         // Plot min Jacobian determinant.
bool flag_plot_max_jacobian = false;         // Plot max Jacobian determinant.
bool flag_plot_min_jacobian_shape = false;   // Plot min Jacobian shape.
bool flag_plot_max_jacobian_shape = false;   // Plot max Jacobian shape.
bool flag_normalize = false;
bool flag_silent_write = false; // if true, suppress message "Writing table..."
int DEFAULT_TABLE_COLUMN_WIDTH = 8;
int DEFAULT_TABLE_PRECISION = 4;
COORD_TYPE edge_interval = 0.05;
COORD_TYPE Jacobian_table_interval = 0.05;
MESH_INFO mesh_info;
IO_INFO io_info;
MESH_DATA mesh_data;

// List of poly with orientation conflicts.
vector<int> orientation_conflict_list;  

vector<int> non_manifold_facet_vert; // list of non_manifold facet vertices
vector<int> non_manifold_edge_vert;  // list of non_manifold edge vertices
vector<int> boundary_facet_vert;     // list of boundary facet vertices

// List of pairs of hex facets which share exactly two edges.
FACET_INFO_PAIRS_ARRAY hex_facet_pairs_sharing_exactly_two_edges;

// list of vertices in facets with mismatched orientations.
vector<int> orientation_mismatch_facet_vert;

// true if boundary facet is inside bounding box
vector<bool> internal_boundary_facet;  

// true if boundary facet is far from bounding_box
vector<bool> far_from_bounding_box;

int num_internal_boundary_facets = 0;
int num_deep_boundary_facets = 0;
vector<bool> in_non_manifold_facet;  // true if vertex is in non_manifold facet
vector<bool> in_non_manifold_edge;   // true if vertex is in non_manifold edge
vector<bool> non_manifold_vert;      // true if vertex is non_manifold
vector<int> non_manifold_vert_list;  // list of non_manifold vertices

vector<int> sorted_poly;            // list of polytopes in sorted order

int vertex_index = 0;
int contains_vertex_index = 0;
int simplex_index = 0;
int poly_index = 0;
int edge_end0_index = 0;
int edge_end1_index = 0;

bool vertex_info_flag = false;
bool simplex_info_flag = false;
bool poly_info_flag = false;
bool manifold_flag = false;
bool oriented_manifold_flag = false;
bool check_facet_intersections_flag = false;
bool terse_flag = false;
bool vlist_flag = false;
bool vlist_min_flag = false;
bool plist_flag = false;
bool contains_vertex_flag = false;
bool contains_edge_flag = false;

// compute info routines
void compute_bounding_box();
int compute_num_edges();
int count_deep_vertices(const int dimension, const std::vector<int> & vlist);
void compute_facet_info();
int count_num_poly(const int num_poly_vert);
int count_num_internal_vertices(const POLYMESH_TYPE & polymesh);


// mesh processing routines
void sort_poly();
void identify_duplicates();
int identify_poly_with_duplicate_vertices();
int identify_duplicate_polytopes();
void identify_polygon_edges();
void set_in_non_manifold_facet();
void set_in_non_manifold_edge();
void identify_hex_sharing_exactly_two_facet_edges
(const POLYMESH_TYPE & polymesh, 
 FACET_INFO_PAIRS_ARRAY & facet_pairs_sharing_exactly_two_edges);
void set_boundary_vertices
(const vector<int> & boundary_facet_vert, POLYMESH_TYPE & polymesh);


// manifold routines
void identify_non_manifold();
void identify_non_manifold_and_boundary_facets();
int identify_non_manifold_vertices();
int identify_non_manifold_edges();
bool is_internal
(const BOUNDING_BOX & bounding_box, const vector<int> & facet_vlist, 
 const int numv_per_facet, const int jf);

// write table routines
void compute_polygon_angles(ANGLE_TABLE & angle_table);
void compute_hex_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const COORD_TYPE edge_interval,
 EDGE_LENGTH_TABLE & edge_length_table);
void compute_hex_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const COORD_TYPE min_table_Jacobian,
 const COORD_TYPE jacobian_interval,
 JACOBIAN_TABLE & jacobian_table);
void compute_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 const COORD_TYPE min_table_Jacobian,
 const COORD_TYPE jacobian_interval,
 JACOBIAN_TABLE & jacobian_table);
void compute_hex_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const COORD_TYPE min_table_Jacobian,
 const COORD_TYPE jacobian_interval,
 JACOBIAN_TABLE & jacobian_table);
void compute_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 const COORD_TYPE min_table_Jacobian,
 const COORD_TYPE jacobian_interval,
 JACOBIAN_TABLE & jacobian_table);
void compute_hex_Jacobian_shape
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const COORD_TYPE min_table_value,
 const COORD_TYPE jacobian_interval,
 JACOBIAN_SHAPE_TABLE & jacobian_shape_table);
void compute_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 const COORD_TYPE min_table_value,
 const COORD_TYPE jacobian_interval,
 JACOBIAN_SHAPE_TABLE & jacobian_shape_table);


void write_angle_table_gplt
(const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table);
void write_edge_length_table_gplt
(const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table);
void write_jacobian_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table);
void write_normalized_jacobian_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table);
void write_jacobian_shape_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table);


// output info routines
void output_general_info
(const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence);
void output_vertex_info(const int dimension, const int vertex_index);
void output_simplex_info(const int dimension, const int simplex_index);
void output_manifold_info();
void output_non_manifold_facets();
void output_non_manifold_edges();
void output_non_manifold_vertices();
void output_poly_with_orientation_conflicts();
void output_internal_boundary_facets();
void output_vertex_list();
void output_simplices();
void output_polytopes();
bool output_self_intersections();
bool output_self_intersections_using_grid_of_bins();
void output_manifold_and_boundary_counts();
void output_hex_facet_pairs_sharing_exactly_two_edge
(const FACET_INFO_PAIRS_ARRAY & facet_pairs_sharing_exactly_two_edges);
void output_poly_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal_poly);

void output_min_max_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal_poly);
void output_min_max_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal_poly, 
 const int num_poly_edges);
void output_min_max_polygon_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, const IO_INFO & io_info,
 const bool flag_internal_poly);
void output_min_max_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, const IO_INFO & io_info,
 const bool flag_internal_poly);
void output_poly_with_min_max_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly);
void output_poly_with_min_max_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly);
void output_min_max_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant);
void output_min_max_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert);
void output_min_max_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant);
void output_min_max_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert);
void output_min_max_Jacobian_shape
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert, 
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape);
void output_min_max_Jacobian_shape
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert);
void write_non_manifold_edges();



// intersection routines
int is_num_poly_vert_constant(const std::vector<int> & num_poly_vert);
bool intersect_edge_triangle
(const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
 const COORD_TYPE w0[DIM3], const COORD_TYPE w1[DIM3],
 const COORD_TYPE w2[DIM3], const COORD_TYPE epsilon,
 COORD_TYPE intersection_point[DIM3]);
bool intersect_triangle_triangle
(const int js1, const int js2, const double epsilon,
 COORD_TYPE intersection_point[DIM3]);

// misc routines
void read_input_file
(const char * input_filename, POLYMESH_TYPE & polymesh);
void memory_exhaustion();
void parse_command_line(int argc, char **argv);
void check_input();


// PARAMETER TYPE
typedef enum
  {POLYFILE_PARAM, MESH_DIM_PARAM, REVERSE_ORIENT_PARAM,
   VERTEX_PARAM, SIMPLEX_PARAM, POLY_PARAM,
   VLIST_PARAM, VLIST_MIN_PARAM,
   PLIST_PARAM, 
   CONTAINSV_PARAM, CONTAINSE_PARAM,
   MANIFOLD_PARAM, ORIENTED_MANIFOLD_PARAM,
   CHECK_FACET_INTERSECTIONS_PARAM,
   SELFI_PARAM, SELFI_NO_GRID_PARAM, GRID_LENGTH_PARAM,
   MULTIVERT_PARAM,
   MINC_PARAM, MAXC_PARAM, MIN_NUMV_PARAM, MAX_NUMV_PARAM,
   ANGLE_LE_PARAM, ANGLE_GE_PARAM,
   FACET_ANGLE_LE_PARAM, FACET_ANGLE_GE_PARAM,
   PJACOBIAN_PARAM, VJACOBIAN_PARAM,
   LIST_DUP_PARAM, INTERNAL_POLY_PARAM, INTERNAL_VERT_PARAM, 
   REPORT_DEEP_PARAM, OUT_VALUES_PARAM,
   OUT_MIN_ANGLE_PARAM, OUT_MAX_ANGLE_PARAM, 
   OUT_MIN_JACOBIAN_DET_PARAM, OUT_MAX_JACOBIAN_DET_PARAM,
   OUT_MIN_NORMALIZED_JACOBIAN_DET_PARAM, 
   OUT_MAX_NORMALIZED_JACOBIAN_DET_PARAM,
   PLOT_ANGLES_PARAM, PLOT_EDGE_LENGTHS_PARAM,
   PLOT_JACOBIAN_PARAM, PLOT_JACOBIAN_SHAPE_PARAM,
   FOR_EACH_TYPE_PARAM, 
   MAX_OUT_PARAM,
   TERSE_PARAM, HELP_PARAM, UNKNOWN_PARAM} PARAMETER;
const char * parameter_string[] = 
  {"-polyfile", "-mesh_dim", "-reverse_orient",
   "-vertex", "-simplex", "-poly",
   "-vlist", "-vlist_min", 
   "-plist",
   "-containsv", "-containse",
   "-manifold", "-oriented_manifold",
   "-check_facetI",
   "-selfI", "-selfI_no_grid", "-grid_length",
   "-multivert",
   "-minc", "-maxc", "-min_numv", "-max_numv",
   "-angle_le", "-angle_ge",
   "-facet_angle_le", "-facet_angle_ge",
   "-pJacobian", "-vJacobian",
   "-list_dup", "-internal", "-internal_vert",
   "-report_deep", "-out_values", "-out_min_angle", "-out_max_angle",
   "-out_min_jdet", "-out_max_jdet",
   "-out_min_normalized_jdet", "-out_max_normalized_jdet",
   "-plot_angles", "-plot_edge_lengths", 
   "-plot_jacobian", "-plot_jshape",
   "-for_each_type",
   "-max_out", "-terse", "-help", "-unknown"};


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  bool passed_all_manifold_tests = false;
  bool passed_boundary_test = false;
  bool flag_self_intersect = false;
  std::string input_filename_prefix;

  std::set_new_handler(memory_exhaustion);

  parse_command_line(argc, argv);
  get_filename_remove_suffix(input_filename, "off", input_filename_prefix);

  try {

    read_input_file(input_filename, polymesh);

    check_input();

    VERTEX_POLY_INCIDENCE_TYPE vertex_info(polymesh);

    sort_poly();
    identify_duplicates();

    if (flag_internal_poly || flag_internal_vert) { compute_facet_info(); }

    if (io_info.flag_general_info) 
      { output_general_info(polymesh, vertex_info); }
    else if (io_info.flag_output_min_angle || io_info.flag_output_max_angle) {
      output_min_max_angle
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly); 
    }
    else if (io_info.flag_output_min_Jacobian_determinant ||
             io_info.flag_output_max_Jacobian_determinant ||
             io_info.flag_output_min_normalized_Jacobian_determinant ||
             io_info.flag_output_max_normalized_Jacobian_determinant ||
             vlist_min_flag) {
      output_min_max_Jacobian_determinants
        (mesh_data, polymesh, vertex_info, vertex_coord, io_info, 
         flag_internal_poly, flag_internal_vert); 
      output_min_max_normalized_Jacobian_determinants
        (mesh_data, polymesh, vertex_info, vertex_coord, io_info, 
         flag_internal_poly, flag_internal_vert); 
      output_min_max_Jacobian_shape
        (mesh_data, polymesh, vertex_info, vertex_coord, io_info, 
         flag_internal_poly, flag_internal_vert); 
    }
    else if (io_info.flag_output_min_Jacobian_determinant ||
             io_info.flag_output_max_Jacobian_determinant) {
      output_min_max_Jacobian_determinants
        (mesh_data, polymesh, vertex_info, vertex_coord, io_info, 
         flag_internal_poly, flag_internal_vert); 
    }

    if (vertex_info_flag) 
      { output_vertex_info(mesh_data.dimension, vertex_index); };
    if (poly_info_flag) { 
      output_poly_info
        (mesh_data.dimension, polymesh, vertex_coord, poly_index); 
    };
    if (simplex_info_flag && !flag_simplex_file) {
      if (simplex_info_flag) { 
        output_poly_info
          (mesh_data.dimension, polymesh, vertex_coord, simplex_index); 
      };
    }
    else {
      if (simplex_info_flag) 
        { output_simplex_info(mesh_data.dimension, simplex_index); };
    }

    if (manifold_flag) { 
      identify_non_manifold();

      if (mesh_info.AreAllNonManifoldZero()) 
        { passed_all_manifold_tests = true; }

      if (flag_report_deep) {
        if (num_deep_boundary_facets == 0) 
          { passed_boundary_test = true; }
      }
      else {
        if (num_internal_boundary_facets == 0) 
          { passed_boundary_test = true; }
      }

      if (flag_output_only_values) {
        output_manifold_and_boundary_counts();
      }
      else {

        bool flag_passed_tests = mesh_info.AreAllZero();
        if (mesh_data.mesh_dimension < mesh_data.dimension) {
          flag_passed_tests = (flag_passed_tests && passed_boundary_test);
        }

        if (terse_flag && flag_passed_tests) {
          if (mesh_data.mesh_dimension < mesh_data.dimension) {
            if (oriented_manifold_flag) 
              { cout << "Passed all manifold, boundary and orientation tests." 
                     << endl; }
            else 
              { cout << "Passed all manifold and boundary tests." << endl; }
          }
          else {
            if (oriented_manifold_flag) 
              { cout << "Passed all manifold and orientation tests." << endl; }
            else
              { cout << "Passed all manifold tests." << endl; }
          }
        }
        else if (mesh_data.mesh_dimension == DIM2 && output_filename != NULL) {
          write_non_manifold_edges();
        }
        else {

          if (terse_flag) {
            if (!passed_all_manifold_tests) 
              { cout << "Failed manifold tests.  "; }
            if (mesh_data.mesh_dimension < mesh_data.dimension) {
              if (flag_report_deep) {
                if (num_deep_boundary_facets != 0) 
                  { cout << "Failed boundary test.";  }
              }
              else {
                if (num_internal_boundary_facets != 0) 
                  { cout << "Failed boundary test.";  }
              }
            }

            if (oriented_manifold_flag) {
              if (mesh_info.num_poly_with_orientation_conflicts != 0)
                { cout << "Failed orientation test.";  }
            }
            cout << endl;
          }

          output_manifold_info();
        };
      };
    }
    else if (check_facet_intersections_flag) {
      if (mesh_data.dimension == DIM3 && flag_cube_file) {

        if (!mesh_data.are_hex_sharing_exactly_two_facet_edges_identified) {
          identify_hex_sharing_exactly_two_facet_edges
            (polymesh, hex_facet_pairs_sharing_exactly_two_edges);
        }

        output_hex_facet_pairs_sharing_exactly_two_edge
          (hex_facet_pairs_sharing_exactly_two_edges);
      }
    }

    if (flag_report_self_intersections) {
      if (flag_simplex_file && num_vert_per_simplex == 3 &&
          mesh_data.dimension == DIM3) {
        flag_self_intersect = output_self_intersections();
      }
    }

    if (vlist_flag) {
      if (io_info.flag_output_min_Jacobian_determinant ||
          io_info.flag_output_min_normalized_Jacobian_determinant ||
          io_info.flag_output_min_Jacobian_shape) {

        if (io_info.flag_output_min_Jacobian_determinant) {
          output_hex_mesh_vertices_with_min_Jacobian_determinants
            (mesh_data, polymesh, vertex_info, vertex_coord,
             io_info, flag_internal_poly, flag_internal_vert);
        }

        if (io_info.flag_output_min_normalized_Jacobian_determinant) {
          output_hex_mesh_vertices_with_min_normalized_Jacobian_determinants
            (mesh_data, polymesh, vertex_info, vertex_coord,
             io_info, flag_internal_poly, flag_internal_vert);
        }

        if (io_info.flag_output_min_Jacobian_shape) {
          output_hex_mesh_vertices_with_min_Jacobian_shape
            (mesh_data, polymesh, vertex_info, vertex_coord,
             io_info, flag_internal_poly, flag_internal_vert);
        }
      }
      else {
        // Output all vertices.
        output_vertex_list();
      }
    }

    if (vlist_min_flag) {

      if (flag_cube_file && mesh_data.mesh_dimension == DIM3) {
        output_hex_mesh_vertices_with_min_Jacobian_determinants
          (mesh_data, polymesh, vertex_info, vertex_coord,
           io_info, flag_internal_poly, flag_internal_vert);
        output_hex_mesh_vertices_with_min_normalized_Jacobian_determinants
          (mesh_data, polymesh, vertex_info, vertex_coord,
           io_info, flag_internal_poly, flag_internal_vert);
          output_hex_mesh_vertices_with_min_Jacobian_shape
            (mesh_data, polymesh, vertex_info, vertex_coord,
             io_info, flag_internal_poly, flag_internal_vert);
      }
      else {
        cerr << "Error. Option -vlist_min only implemented for hexahedral mesh."
             << endl;
      }
    }

    if (flag_polyfile) {
      if (contains_vertex_flag)
        { output_polytopes(); }

      if (flag_list_duplicate_vertices) {
        output_duplicate_vertices(mesh_data.dimension, polymesh, vertex_coord); 
        cout << endl;
      }

      if (flag_list_duplicate_poly) {
        output_duplicate_poly(polymesh, mesh_info); 
        cout << endl;
      }
    }
    else {
      if (contains_vertex_flag || contains_edge_flag)
        { output_simplices(); }
    }

    if (plist_flag) {
      output_poly_angles
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      output_poly_with_min_max_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      output_poly_with_min_max_Jacobian_determinants
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }

    if (is_mesh_dimension_set || !flag_polyfile) {
      if (mesh_data.mesh_dimension == DIM2 && flag_plot_angles) {

        write_angle_table_gplt
          (input_filename_prefix, flag_plot_min_polygon_angles,
           flag_plot_max_polygon_angles);
      }
    }

    if (flag_plot_min_edge_lengths || flag_plot_max_edge_lengths) {
      if (mesh_data.mesh_dimension == DIM3 && flag_cube_file) {

        write_edge_length_table_gplt
          (input_filename_prefix, flag_plot_min_edge_lengths,
           flag_plot_max_edge_lengths);
      }
    }


    if (flag_plot_min_jacobian || flag_plot_max_jacobian) {
      if (mesh_data.mesh_dimension == DIM3 && flag_cube_file) {

        write_jacobian_table_gplt
          (mesh_data, polymesh, vertex_info, input_filename_prefix, 
           flag_plot_min_jacobian, flag_plot_max_jacobian);
        write_normalized_jacobian_table_gplt
          (mesh_data, polymesh, vertex_info, input_filename_prefix, 
           flag_plot_min_jacobian, flag_plot_max_jacobian);
      }
    }

    if (flag_plot_min_jacobian_shape || flag_plot_max_jacobian_shape) {
      if (mesh_data.mesh_dimension == DIM3 && flag_cube_file) {

        write_jacobian_shape_table_gplt
          (mesh_data, polymesh, vertex_info, input_filename_prefix, 
           flag_plot_min_jacobian_shape, 
           flag_plot_max_jacobian_shape);
      }
    }

  }
  catch (ERROR & error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  }

  delete [] simplex_vert;;
  delete [] vertex_coord;

  if (manifold_flag && !passed_all_manifold_tests) 
    { return(1); }

  if (oriented_manifold_flag &&
      mesh_info.num_poly_with_orientation_conflicts > 0)
    { return(1); }

  if (manifold_flag && (mesh_data.mesh_dimension < mesh_data.dimension) && 
      !passed_boundary_test) {
    return(1);
  }
  else if (flag_report_self_intersections && flag_self_intersect) {
    return(1);
  }
  else {
    return(0);
  }

}

// **************************************************
// READ INPUT FILE ROUTINES
// **************************************************

void read_input_file
(const char * input_filename, POLYMESH_TYPE & polymesh)
{
  const int NUMV_PER_CUBE = 8;
  IJK::PROCEDURE_ERROR error("read_input_file");

  if (input_filename == NULL) {
    error.AddMessage("Programming error. Input filename not set.");
    throw error;
  }

  ifstream in(input_filename, ios::in);
  if (!in.good()) {
    cerr << "Unable to open file " << input_filename << "." << endl;
    exit(30);
  };

  ijkinPolytopeOFF
    (in, mesh_data.dimension, vertex_coord, mesh_data.num_vertices, 
     polymesh.list_length, polymesh.element, polymesh.first_element);
  polymesh.poly_data.resize(polymesh.NumPoly());
  polymesh.vertex_data.resize(mesh_data.num_vertices);

  in.close();

  if (!polymesh.Check(error)) { throw error; }

  num_poly = polymesh.NumPoly();

  if (num_poly == 0) { return; }

  int min_num_vert = 
    *(min_element(polymesh.list_length.begin(), polymesh.list_length.end()));
  int max_num_vert = 
    *(max_element(polymesh.list_length.begin(), polymesh.list_length.end()));

  if (min_num_vert == max_num_vert) {

    num_vert_per_poly = min_num_vert;

    if (!is_mesh_dimension_set) {

      if (mesh_data.dimension == DIM3 && num_vert_per_poly == 4) {
        cerr << "Unable to determine mesh dimension." << endl;
        cerr << "  Input polytopes could be tetrahedra or quadrilaterals."
             << endl;
        cerr << "Use option -mesh_dim <mdim>." << endl;
        exit(20);
      }

      if (num_vert_per_poly <= mesh_data.dimension+1) {
        // Assume poly are all simplices.
        mesh_data.mesh_dimension = num_vert_per_poly-1;
      }
      else if (num_vert_per_poly == NUMV_PER_CUBE) {
        mesh_data.mesh_dimension = DIM3;
        flag_cube_file = true;
      }
      else {
        mesh_data.mesh_dimension = mesh_data.dimension;
      }
    }

    if (num_vert_per_poly <= mesh_data.mesh_dimension+1) {

      // Copy poly_vert into simplex_vert.
      simplex_vert = new int[polymesh.element.size()];
      num_vert_per_simplex = num_vert_per_poly;
      num_simplices = num_poly;

      std::copy(polymesh.element.begin(), polymesh.element.end(),
                simplex_vert);
      flag_simplex_file = true;
    }
    else if (num_vert_per_poly == NUMV_PER_CUBE) {
      flag_cube_file = true;
      flag_polyfile = true;
    }
    else {
      flag_polyfile = true;
    }

    if (flag_cube_file && is_mesh_dimension_set) {
      const int numv_per_cube = (1 << mesh_data.mesh_dimension);
      if (numv_per_cube != num_vert_per_poly) {
        const int cube_dimension = 
          IJK::compute_cube_dimension_from_num_vertices(num_vert_per_poly);
        cerr << "Usage error.  Mismatch between input file of cubes and -mesh_dim argument." << endl;
        cerr << "  Input file cubes have "
             << numv_per_cube << " vertices and dimension "
             << cube_dimension << "." << endl;
        cerr << "  Argument of -mesh_dim is " << mesh_data.mesh_dimension 
             << "." << endl;
        cerr << "  Change argument of -mesh_dim to " << cube_dimension
             << "." << endl;
        exit(20);
      }
    }

    if (is_mesh_dimension_set) {
      if (num_vert_per_poly <= mesh_data.mesh_dimension) {
        cerr << "Warning:  All polygons in input file are degenerate."
             << endl;
        cerr << "  Argument of -mesh_dim is " << mesh_data.mesh_dimension 
             << "." << endl;
        cerr << "  All polygons have "  << num_vert_per_poly
             << " vertices." << endl;
        cerr << endl;
      }
    }

  }
  else {
    flag_polyfile = true;

    if (!is_mesh_dimension_set) {
      mesh_data.mesh_dimension = mesh_data.dimension-1;

      cerr << "Warning: Unable to determine mesh dimension from input file."
           << endl;
      cerr << "  Assuming mesh dimension is " 
           << mesh_data.mesh_dimension << "." << endl;
      cerr << "  Use option \"-mesh_dim {mdim}\" to specify mesh dimension."
           << endl;
      cerr << endl;
    }
  }

}


int is_num_poly_vert_constant(const std::vector<int> & num_poly_vert)
{
  if (num_poly_vert.size() == 0) { 
    // Trivially true.
    return(true); 
  }

  const int num_poly_vert0 = num_poly_vert[0];

  for (int i = 1; i < num_poly_vert.size(); i++) {
    if (num_poly_vert0 != num_poly_vert[i])
      { return(false); }
  }

  return(true);
}


// **************************************************
// OUTPUT INFO ROUTINES
// **************************************************

void output_general_info
(const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  const int num_vertices = mesh_data.num_vertices;

  compute_bounding_box();
  if (!flag_polyfile) {
    num_edges = compute_num_edges();
  }

  cout << "Volume dimension: " << dimension << endl;
  cout << "Mesh dimension: " << mesh_dimension << endl;
  cout << "Number of mesh vertices: " << num_vertices << endl;
  if (flag_internal_vert) {
    const int num_internal_vertices = count_num_internal_vertices(polymesh);
    cout << "Number of internal mesh vertices: " 
         << num_internal_vertices << endl;
  }
  if (flag_polyfile) {
    if (mesh_dimension <= DIM2) {
      cout << "Number of mesh polygons: " << num_poly << endl;
    }
    else {
      cout << "Number of mesh polytopes: " << num_poly << endl;
    }
  }
  else {
    if (dimension > DIM2) {
      cout << "Number of mesh edges: " << num_edges << endl;
      if (mesh_dimension > DIM2) {
        cout << "Number of mesh simplices: " << num_simplices << endl;
      }
      else {
        cout << "Number of mesh triangles: " << num_simplices << endl;
      }
    }
    else if (dimension == DIM2) {
      cout << "Number of mesh edges: " << num_simplices << endl;
    }

    if (dimension == DIM3) {
      cout << "#V - #E + #F = " 
           << num_vertices - num_edges + num_simplices 
           << endl;
    }
  }

  if (mesh_info.num_poly_with_duplicate_vertices > 0) {
    cout << "Number of poly with duplicate vertices: " 
         << mesh_info.num_poly_with_duplicate_vertices << endl;
  }

  if (mesh_info.num_duplicate_poly > 0) {
    cout << "Number of duplicate poly: " 
         << mesh_info.num_duplicate_poly << endl;
  }

  io_info.flag_output_min_angle = true;
  io_info.flag_output_max_angle = true;
  output_min_max_angle(mesh_data, polymesh, vertex_coord, io_info, false);

  if (flag_internal_poly) {
    output_min_max_angle(mesh_data, polymesh, vertex_coord, io_info, true);
  }

  io_info.flag_output_min_edge_length = true;
  io_info.flag_output_max_edge_length = true;
  output_min_max_edge_lengths
    (mesh_data, polymesh, vertex_coord, io_info, false);
  if (flag_internal_poly) {
    output_min_max_edge_lengths
      (mesh_data, polymesh, vertex_coord, io_info, true);
  }

  io_info.flag_output_min_Jacobian_determinant = true;
  io_info.flag_output_max_Jacobian_determinant = true;
  output_min_max_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     io_info, false, false);
  if (flag_internal_poly) {
    output_min_max_Jacobian_determinants
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info, true, false);
  }
  if (flag_internal_vert) {
    output_min_max_Jacobian_determinants
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info, false, true);
  }

  io_info.flag_output_min_normalized_Jacobian_determinant = true;
  io_info.flag_output_max_normalized_Jacobian_determinant = true;
  output_min_max_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     io_info, false, false);
  if (flag_internal_poly) {
    output_min_max_normalized_Jacobian_determinants
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info, true, false);
  }
  if (flag_internal_vert) {
    output_min_max_normalized_Jacobian_determinants
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info, false, true);
  }

  io_info.flag_output_min_Jacobian_shape = true;
  io_info.flag_output_max_Jacobian_shape = true;
  output_min_max_Jacobian_shape
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     io_info, false, false);
  if (flag_internal_poly) {
    output_min_max_Jacobian_shape
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info, true, false);
  }
  if (flag_internal_vert) {
    output_min_max_Jacobian_shape
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       io_info, false, true);
  }

  cout << "Bounding box: (";
  IJK::print_list(cout, bounding_box.MinCoord(), bounding_box.Dimension());
  cout << " ";
  IJK::print_list(cout, bounding_box.MaxCoord(), bounding_box.Dimension());
  cout << ")" << endl;

  cout << endl;
}

void output_min_max_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal_poly)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (flag_simplex_file && mesh_dimension == DIM3) {
    output_min_max_dihedral_angle
      (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
  }

  output_min_max_angle
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly, 0);

  if (flag_for_each_type) {
    for (int num_poly_edges = 3; num_poly_edges < 10; num_poly_edges++) {
      int npoly = count_num_poly(num_poly_edges);
      if (npoly > 0) {
        output_min_max_angle
          (mesh_data, polymesh, vertex_coord, io_info, 
           flag_internal_poly, num_poly_edges);
      }
    }
  }

  if (flag_simplex_file && mesh_dimension == DIM3) {
    output_dihedral_angle_count
      (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    output_tetrahedra_facet_angle_count
      (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
  }

}


void output_min_max_angle
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info,
 const bool flag_internal_poly, const int num_poly_edges)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (mesh_dimension == DIM2) {
    output_min_max_polygon_angle
      (mesh_data, polymesh, vertex_coord, io_info, 
       flag_internal_poly, num_poly_edges);
  }
  else if (flag_simplex_file && mesh_dimension == DIM3) {
    output_min_max_tetrahedra_facet_angle
      (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
  }
}


// Output polytopes with min/max angles or angles less then and equal to
//   or greater than and equal to given bounds.
void output_poly_angles
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh, const COORD_TYPE * vertex_coord, 
 const IO_INFO & io_info, const bool flag_internal_poly)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (flag_internal_poly &&
      !mesh_data.are_boundary_facets_identified)
    { identify_non_manifold(); }

  if (io_info.angle_le.IsSet() || io_info.angle_ge.IsSet() ||
      io_info.facet_angle_le.IsSet() || io_info.facet_angle_ge.IsSet() ||
      io_info.flag_output_min_angle || io_info.flag_output_max_angle) {

    if (io_info.angle_le.IsSet()) {
      if (mesh_dimension == DIM2) {
        output_polygons_with_small_angles
          (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
           io_info.angle_le.Value());
      }
      else if (mesh_dimension == DIM3 && flag_simplex_file) {
        output_tetrahedra_with_small_dihedral_angles
          (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
           io_info.angle_le.Value());
      }
    }
    else if (io_info.flag_output_min_angle) {
      if (mesh_dimension == DIM2) {
        output_polygons_with_min_angle
          (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      }
      else if (mesh_dimension == DIM3 && flag_simplex_file) {
        output_tetrahedra_with_min_dihedral_angle
          (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      }
    }

    if (io_info.angle_ge.IsSet()) {
      if (mesh_dimension == DIM2) {
        output_polygons_with_large_angles
          (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
           io_info.angle_ge.Value());
      }
      else if (mesh_dimension == DIM3 && flag_simplex_file) {
        output_tetrahedra_with_large_dihedral_angles
          (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
           io_info.angle_ge.Value());
      }
    }
    else if (io_info.flag_output_max_angle) {
      if (mesh_dimension == DIM2) {
        output_polygons_with_max_angle
          (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      }
      else if (mesh_dimension == DIM3 && flag_simplex_file) {
        output_tetrahedra_with_max_dihedral_angle
          (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      }
    }

    if (io_info.facet_angle_le.IsSet()) {
      if (mesh_dimension == DIM3 && flag_simplex_file) {
        output_tetrahedra_with_small_facet_angles
          (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
           io_info.facet_angle_le.Value());
      }
    }
    else if (io_info.flag_output_min_angle) {
      if (mesh_dimension == DIM3 && flag_simplex_file) {
        output_tetrahedra_with_min_facet_angle
          (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      }
    }

    if (io_info.facet_angle_ge.IsSet()) {
      if (mesh_dimension == DIM3 && flag_simplex_file) {
        output_tetrahedra_with_large_facet_angles
          (mesh_data, polymesh, vertex_coord,
           flag_internal_poly, io_info.facet_angle_ge.Value());
      }
    }
    else if (io_info.flag_output_min_angle) {
      if (mesh_dimension == DIM3 && flag_simplex_file) {
        output_tetrahedra_with_max_facet_angle
          (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      }
    }

  }
  else if (io_info.flag_output_all_min_max) {

    if (mesh_dimension == DIM2) {
      output_polygons_with_min_max_angles
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }
    else if (mesh_dimension == DIM3 && flag_simplex_file) {
      output_tetrahedra_with_min_max_dihedral_angles
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
      output_tetrahedra_with_min_max_facet_angles
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }
  }
}


void output_min_max_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly)
{
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (mesh_dimension == DIM2) {
    output_min_max_polygon_edge_lengths
      (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
  }
  else if (flag_simplex_file) {
    if (mesh_dimension == DIM3) {
      output_min_max_tetrahedra_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }
    else {
      output_min_max_simplices_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }
  }
  else if (flag_cube_file && mesh_dimension == DIM3) {
    output_min_max_hexahedra_edge_lengths
      (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
  }
}

void output_poly_with_min_max_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (mesh_dimension == 2) {

    if (io_info.flag_output_min_edge_length || 
        io_info.flag_output_all_min_max) {
      output_polygons_with_min_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }

    if (io_info.flag_output_max_edge_length || 
        io_info.flag_output_all_min_max) {
      output_polygons_with_max_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }
  }
  else if (flag_simplex_file) {

    if (io_info.flag_output_min_edge_length || 
        io_info.flag_output_all_min_max) {
      output_tetrahedra_with_min_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }

    if (io_info.flag_output_max_edge_length || 
        io_info.flag_output_all_min_max) {
      output_tetrahedra_with_max_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }
  }
  else if (flag_cube_file && mesh_dimension == DIM3) {
    if (io_info.flag_output_min_edge_length || 
        io_info.flag_output_all_min_max) {
      output_hexahedra_with_min_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }

    if (io_info.flag_output_max_edge_length || 
        io_info.flag_output_all_min_max) {
      output_hexahedra_with_max_edge_lengths
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }
  }
}


void output_poly_with_min_max_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (flag_cube_file && mesh_dimension == DIM3 && dimension == DIM3) {
    if (io_info.flag_output_min_edge_length || 
        io_info.flag_output_all_min_max) {
      output_hexahedra_with_min_Jacobian_determinants
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }

    if (io_info.flag_output_max_edge_length || 
        io_info.flag_output_all_min_max) {
      output_hexahedra_with_max_Jacobian_determinants
        (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly);
    }
  }
}


void output_min_max_polygon_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord, const IO_INFO & io_info,
 const bool flag_internal_poly)
{
  output_min_max_polygon_edge_lengths
    (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly, 0);

  if (flag_for_each_type) {
    for (int num_poly_vert = 3; num_poly_vert < 10; num_poly_vert++) {
      int npoly = count_num_poly(num_poly_vert);
      if (npoly > 0) {
        output_min_max_polygon_edge_lengths
          (mesh_data, polymesh, vertex_coord, io_info, 
           flag_internal_poly, num_poly_vert);
      }
    }
  }
}


void output_min_max_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert, 
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (flag_cube_file && mesh_dimension == DIM3 && dimension == DIM3) {

    if (flag_internal_poly || !flag_internal_vert) {
      if (io_info.flag_pJacobian.IsSetAndTrue()) {
        output_min_max_hexahedra_Jacobian_determinants
          (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly, 
           min_Jacobian_determinant, max_Jacobian_determinant);
      }

      if (io_info.flag_vJacobian.IsSetAndTrue()) {
        output_min_max_hex_vert_Jacobian_determinants
          (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
           io_info, flag_internal_poly, flag_internal_vert,
           min_Jacobian_determinant, max_Jacobian_determinant);
      }
    }
    else if (flag_internal_vert) {
      output_min_max_hex_vert_Jacobian_determinants
        (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
         io_info, flag_internal_poly, flag_internal_vert,
         min_Jacobian_determinant, max_Jacobian_determinant);
    }
  }

}


void output_min_max_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert)
{
  COORD_TYPE min_Jacobian_determinant, max_Jacobian_determinant;

  output_min_max_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, io_info,
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_determinant, max_Jacobian_determinant);
}


void output_min_max_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert, 
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (flag_cube_file && mesh_dimension == DIM3 && dimension == DIM3) {

    if (flag_internal_poly || !flag_internal_vert) {
      if (io_info.flag_pJacobian.IsSetAndTrue()) {
        output_min_max_hexahedra_normalized_Jacobian_determinants
          (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly, 
           min_Jacobian_determinant, max_Jacobian_determinant);
      }

      if (io_info.flag_vJacobian.IsSetAndTrue()) {
        output_min_max_hex_vert_normalized_Jacobian_determinants
          (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
           io_info, flag_internal_poly, flag_internal_vert,
           min_Jacobian_determinant, max_Jacobian_determinant);
      }
    }
    else if (flag_internal_vert) {
      output_min_max_hex_vert_normalized_Jacobian_determinants
        (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
         io_info, flag_internal_poly, flag_internal_vert,
         min_Jacobian_determinant, max_Jacobian_determinant);
    }
  }

}


void output_min_max_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert)
{
  COORD_TYPE min_Jacobian_determinant, max_Jacobian_determinant;

  output_min_max_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, io_info,
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_determinant, max_Jacobian_determinant);
}


void output_min_max_Jacobian_shape
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert, 
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape)
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (flag_cube_file && mesh_dimension == DIM3 && dimension == DIM3) {

    if (flag_internal_poly || !flag_internal_vert) {
      if (io_info.flag_pJacobian.IsSet()) {
        if (io_info.flag_pJacobian.Value()) {
          output_min_max_hexahedra_Jacobian_shape
            (mesh_data, polymesh, vertex_coord, io_info, flag_internal_poly, 
             min_Jacobian_shape, max_Jacobian_shape);
        }
      }

      if (io_info.flag_vJacobian.IsSet()) {
        if (io_info.flag_vJacobian.Value()) {
          output_min_max_hex_vert_Jacobian_shape
            (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
             io_info, flag_internal_poly, flag_internal_vert,
             min_Jacobian_shape, max_Jacobian_shape);
        }
      }
    }
    else if (flag_internal_vert) {
      output_min_max_hex_vert_Jacobian_shape
        (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
         io_info, flag_internal_poly, flag_internal_vert,
         min_Jacobian_shape, max_Jacobian_shape);
    }
  }

}


void output_min_max_Jacobian_shape
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
 const COORD_TYPE * vertex_coord,  const IO_INFO & io_info,
 const bool flag_internal_poly, const bool flag_internal_vert)
{
  COORD_TYPE min_Jacobian_shape, max_Jacobian_shape;

  output_min_max_Jacobian_shape
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, io_info,
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_shape, max_Jacobian_shape);
}


void output_vertex_info(const int dimension, const int vertex_index)
{
  const int num_vertices = mesh_data.num_vertices;

  cout << "Vertex: " << vertex_index << endl;

  if (vertex_index < 0 || vertex_index >= num_vertices)
    {
      cout << "  Illegal vertex index " << vertex_index
           << ".  Vertex index should be in range["
           << 0 << "," << num_vertices-1 << "]." << endl;
      return;
    };

  cout << "  Coordinates: ";
  IJK::print_list
    (cout, vertex_coord+vertex_index*dimension, dimension);
  cout << endl;

  cout << endl;
}


void output_simplex_info(const int dimension, const int simplex_index)
{
  const int numv_per_simplex = num_vert_per_poly;

  cout << "Simplex: " << simplex_index << endl;

  if (simplex_index < 0 || simplex_index >= num_simplices) {
    cout << "  Illegal simplex index.  Simplex index should be in range["
         << 0 << "," << num_simplices-1 << "]." << endl;
    return;
  };

  cout << "  Vertices:" << endl;
  for (int k = 0; k < numv_per_simplex; k++) {
    int iv = simplex_vert[simplex_index*numv_per_simplex+k];
    cout << "    " << setw(6) << iv << "  ";
    IJK::print_list
      (cout, vertex_coord+iv*dimension, dimension);
    cout << endl;
  }

  cout << endl;
}


void output_manifold_info()
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  const int NUMV_PER_CUBE = 8;
  const int NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;
  const char * poly_str = "polytopes";
  const char * facet_str = "facets";
  bool flag_newline(false);

  if (mesh_dimension == DIM2) { 
    poly_str = "polygons"; 
    facet_str = "edges";
  }

  if (mesh_info.num_poly_with_duplicate_vertices == 0) {
    if (!terse_flag) {
      cout << "No degenerate " << poly_str << "." << endl;
    }
  }
  else {
    if (!terse_flag) { cout << endl; }
    cout << "Num degenerate " << poly_str << ":  " 
         << mesh_info.num_poly_with_duplicate_vertices << endl;

    if (!terse_flag) {
      output_degenerate_poly(polymesh, mesh_info);
      cout << endl;
      flag_newline = true;
    }
  }

  if (mesh_info.num_duplicate_poly == 0) {
    if (!terse_flag) {
      cout << "No duplicate " << poly_str << "." << endl;
    }
    flag_newline = false;
  }
  else {
    if (!flag_newline && !terse_flag) { cout << endl; }
    cout << "Num duplicate " << poly_str << ":  "
         << mesh_info.num_duplicate_poly << endl;

    if (!terse_flag) {
      output_duplicate_poly(polymesh, mesh_info);
      cout << endl;
      flag_newline = true;
    }
  }

  if (mesh_data.dimension == DIM3 && flag_cube_file) {

    if (hex_facet_pairs_sharing_exactly_two_edges.size() > 0) {
      if (!flag_newline && !terse_flag) { cout << endl; }
    }

    output_hex_facet_pairs_sharing_exactly_two_edge
      (hex_facet_pairs_sharing_exactly_two_edges);

    if (hex_facet_pairs_sharing_exactly_two_edges.size() > 0 &&
        !terse_flag) {
      cout << endl;
      flag_newline = true;
    }
  }

  if (mesh_dimension == DIM2 || flag_simplex_file || flag_cube_file) {

    int numv_per_facet;
    if (mesh_dimension == DIM2 || flag_simplex_file) 
      { numv_per_facet = mesh_dimension;  }
    else 
      { numv_per_facet = NUMV_PER_CUBE_FACET; }

    const int num_non_manifold_facets = 
      non_manifold_facet_vert.size()/numv_per_facet;

    if (num_non_manifold_facets == 0) {
      if (!terse_flag) {
        cout << "No non-manifold " << facet_str << "." << endl;
        flag_newline = false;
      }
    }
    else {
      if (!flag_newline && !terse_flag) { cout << endl; }
      cout << "Num non-manifold " << facet_str << ":  " 
           << num_non_manifold_facets << endl;
      if (!terse_flag) {
        output_non_manifold_facets();
        cout << endl;
        flag_newline = true;
      }
    }

    if (flag_cube_file && mesh_dimension > 2) {
      // If mesh_dimension == DIM2, then edges are reported as facets.
      // Non-manifold edge detection only implemented for cubes.

      if (mesh_info.num_non_manifold_edges == 0) {
        if (!terse_flag) {
          cout << "No non-manifold edges." << endl;
          flag_newline = false;
        }
      }
      else {
        if (!flag_newline && !terse_flag) { cout << endl; }
        cout << "Num non-manifold edges: " 
             << mesh_info.num_non_manifold_edges << endl;
        if (!terse_flag) {
          output_non_manifold_edges();
          cout << endl;
          flag_newline = true;
        }

      }
    }

    if (mesh_info.num_non_manifold_vertices == 0) {
      if (!terse_flag) {
        cout << "No non-manifold vertices." << endl;
        flag_newline = false;
      }
    }
    else {
      if (!flag_newline && !terse_flag) { cout << endl; }

      int num_deep = count_deep_vertices(dimension, non_manifold_vert_list);
      cout << "Num non-manifold vertices:  " 
           << mesh_info.num_non_manifold_vertices << endl;

      if (flag_report_deep) {
        cout << "Num non-manifold vertices at least " << contract_margin
             << " from bounding box boundary: " << num_deep << endl;
      }

      if (!terse_flag) {
        output_non_manifold_vertices(); 
        cout << endl;
        flag_newline = true;
      }
    }

    if (oriented_manifold_flag) {
      if (mesh_info.num_poly_with_orientation_conflicts == 0) {
        if (!terse_flag) {
          cout << "All polytope orientations match." << endl;
          flag_newline = false;
        }
      }
      else {
        if (!flag_newline && !terse_flag) { cout << endl; }

        cout << "Num " << poly_str
             << " with conflicting orientations: "
             << mesh_info.num_poly_with_orientation_conflicts << endl;

        if (!terse_flag) {
          output_poly_with_orientation_conflicts(); 
          cout << endl;
          flag_newline = true;
        }
      }
    }

    if (mesh_dimension < dimension) { output_internal_boundary_facets(); }
  }
  else {
    if (!terse_flag) {
      cout << "Unable to determine " << poly_str 
           << " facets to check manifold/boundary conditions." << endl;

    }
  }

}


void output_non_manifold_facets()
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  const int NUMV_PER_CUBE = 8;
  const int NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;
  const char * poly_str = "Polytopes";
  const char * facet_str = "facets";
  int numv_per_facet = mesh_dimension;

  if (mesh_dimension == DIM2 || flag_simplex_file) 
    { numv_per_facet = mesh_dimension;  }
  else 
    { numv_per_facet = NUMV_PER_CUBE_FACET; }

  if (mesh_dimension == DIM2) {
    poly_str = "Polygons";
    facet_str = "edges"; 
  }

  int num_non_manifold_facets = 
    non_manifold_facet_vert.size()/numv_per_facet;

  cout << "Non-manifold " << facet_str << ":" << endl;

  for (int jf = 0; jf < num_non_manifold_facets; jf++) {
    cout << "  ";
    print_list(cout, &non_manifold_facet_vert.front()+jf*numv_per_facet,
               numv_per_facet);

    if (mesh_dimension == 2) {
      cout << "  (";
      for (int i = 0; i < numv_per_facet; i++) {
        int iv = non_manifold_facet_vert[jf*numv_per_facet+i];
        print_list(cout, vertex_coord+iv*dimension, dimension);
        if (i+1 < numv_per_facet) { cout << ","; }
      }
      cout << ")";
    }

    cout << endl;
  }
  cout << endl;

  cout << poly_str << " containing non-manifold " << facet_str 
       << ": " << endl;

  int num_output = 0;
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {
    if (polymesh.poly_data[ipoly].ContainsNonManifoldFacet()) {
      cout << "  " << ipoly;
      num_output++;

      if (num_output%10 == 0)
        cout << endl;
    }
  }
  if (num_output%10 != 0) { cout << endl; };
}


void output_non_manifold_edges()
{
  const int DIM3(3);
  const int mesh_dimension = mesh_data.mesh_dimension;
  const char * poly_str = "Cubes";
  const char * edge_str = "edges";
  const int num_non_manifold_edges = non_manifold_edge_vert.size()/2;

  if (mesh_dimension == DIM3) 
    { poly_str = "Hexahedra"; }

  cout << "Non-manifold " << edge_str << ":" << endl;

  for (int je = 0; je < num_non_manifold_edges; je++) {
    cout << "  ";
    print_list(cout, &non_manifold_edge_vert.front()+je*2, 2);
    cout << endl;
  }
  cout << endl;

  cout << poly_str << " containing non-manifold edges:" << endl;

  int num_output = 0;
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {
    if (polymesh.poly_data[ipoly].ContainsNonManifoldEdge()) {
      cout << "  " << ipoly;
      num_output++;

      if (num_output%10 == 0)
        cout << endl;
    }
  }
  if (num_output%10 != 0) { cout << endl; };
}


void output_poly_with_orientation_conflicts()
{
  cout << "Poly whose orientation does not match some neighbor:" << endl;

  for (int i = 0; i < orientation_conflict_list.size(); i++) {
    cout << "  " << orientation_conflict_list[i];
  }
  cout << endl;
}


void output_internal_boundary_facets()
{
  const int mesh_dimension = mesh_data.mesh_dimension;
  int numv_per_facet;
  IJK::PROCEDURE_ERROR error("output_internal_boundary_facet");

  if (flag_cube_file) {
    CUBE_TYPE cube(mesh_dimension);
    numv_per_facet = cube.NumFacetVertices();
  }
  else if (flag_simplex_file) {
    numv_per_facet = mesh_dimension;
  }
  else if (mesh_dimension == DIM2) {
    numv_per_facet = 2;
  }
  else {
    error.AddMessage
      ("Programming error.  Unable to determine number of vertices per facet.");
    throw error;
  }

  if (!terse_flag) {
    cout << "Bounding box: (";
    cout << bounding_box.MinCoord(0) << ", "
         << bounding_box.MinCoord(1) << ", "
         << bounding_box.MinCoord(2) << ")  ("
         << bounding_box.MaxCoord(0) << ", "
         << bounding_box.MaxCoord(1) << ", "
         << bounding_box.MaxCoord(2) << ")" << endl;
  }

  if (boundary_facet_vert.size() == 0) {
    if (!terse_flag) {
      cout << "Surface has no boundary." << endl;
    }
  }
  else if (num_internal_boundary_facets == 0) {
    if (!terse_flag) {
      cout << "Surface boundary lies on bounding box boundary." << endl;
    }
  }
  else if (!flag_report_deep || num_deep_boundary_facets > 0) {
    if (!flag_report_deep) {
      cout << "Number of surface boundary facets inside bounding box: "
           << num_internal_boundary_facets << endl;
    }
    cout << "Number of surface boundary facets at least " 
         << contract_margin << " from bounding box boundary: "
         << num_deep_boundary_facets << endl;

    if (!terse_flag) {

      if (flag_report_deep) {
        cout << "Facets on surface boundary and FAR from bounding box boundary:"
             << endl;

        for (int jf = 0; jf < internal_boundary_facet.size(); jf++) {
          if (far_from_bounding_box[jf]) {
            for (int i = jf*numv_per_facet; 
                 i < (jf+1)*numv_per_facet; i++) {
              cout << "  " << boundary_facet_vert[i];
            }
            cout << endl;
          }
        }
      }
      else {
        cout << "Facets on surface boundary but not on bounding box boundary:"
             << endl;

        for (int jf = 0; jf < internal_boundary_facet.size(); jf++) {
          if (internal_boundary_facet[jf]) {
            for (int i = jf*numv_per_facet; 
                 i < (jf+1)*numv_per_facet; i++) {
              cout << "  " << boundary_facet_vert[i];
            }
            cout << endl;
          }
        }
      }


      cout << endl;
    }
  }
}

void output_hex_facet_pairs_sharing_exactly_two_edge
(const FACET_INFO_PAIRS_ARRAY & facet_pairs_sharing_exactly_two_edges)
{
  if (facet_pairs_sharing_exactly_two_edges.size() == 0) {
    if (!terse_flag) {
      cout << "No irregular hex facet intersections." << endl;
    }
  }
  else {
    cout << "Number of irregular hex facet intersections: "
         << facet_pairs_sharing_exactly_two_edges.size() << endl;

    if (!terse_flag) {
      cout << "Hex pairs whose facets intersect in exactly two edges:" << endl;
      for (int i = 0; i < facet_pairs_sharing_exactly_two_edges.size(); i++) {
        cout << "  (" << facet_pairs_sharing_exactly_two_edges[i].first.poly_containing_face
             << "," << facet_pairs_sharing_exactly_two_edges[i].second.poly_containing_face
             << ")";
      }
      cout << endl;
    }
  }
}


// Output intersection point between triangles
void output_triangle_triangle_intersection
(const int js1, const int js2, const COORD_TYPE intersection_point[DIM3])
{
  cout << "Triangles " << js1 << " and " << js2
       << " intersect at: ";
  print_coord3D(cout, intersection_point, "\n");
}

// Compute intersection between triangles, if one exists
// Return true if intersection found.
bool intersect_triangle_triangle
(const int js1, const int js2, const double epsilon,
 COORD_TYPE intersection_point[DIM3])
{
  COORD_TYPE_PTR w0, w1, w2, v0, v1;

  w0 = vertex_coord + DIM3*simplex_vert[js1*DIM3];
  w1 = vertex_coord + DIM3*simplex_vert[js1*DIM3+1];
  w2 = vertex_coord + DIM3*simplex_vert[js1*DIM3+2];

  for (int k = 0; k < 3; k++) {
    v0 = vertex_coord + DIM3*simplex_vert[js2*DIM3+k];
    int k2 = (k+1)%DIM3;
    v1 = vertex_coord + DIM3*simplex_vert[js2*DIM3+k2];
        
    if (intersect_edge_triangle
        (v0, v1, w0, w1, w2, epsilon, intersection_point))
      { return(true); }
  }

  w0 = vertex_coord + DIM3*simplex_vert[js2*DIM3];
  w1 = vertex_coord + DIM3*simplex_vert[js2*DIM3+1];
  w2 = vertex_coord + DIM3*simplex_vert[js2*DIM3+2];

  for (int k = 0; k < 3; k++) {
    v0 = vertex_coord + DIM3*simplex_vert[js1*DIM3+k];
    int k2 = (k+1)%DIM3;
    v1 = vertex_coord + DIM3*simplex_vert[js1*DIM3+k2];

    if (intersect_edge_triangle
        (v0, v1, w0, w1, w2, epsilon, intersection_point))
      { return(true); }
  }

  return(false);
}

bool output_self_intersections()
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE intersection_point[DIM3];

  if (dimension != DIM3) { return(false); }
  if (num_vert_per_simplex != 3) { return(false); }

  bool flag_intersect = false;
  if (flag_use_grid_of_bins) {
    flag_intersect = output_self_intersections_using_grid_of_bins();
  }
  else {
    for (int js1 = 0; js1 < num_simplices; js1++) {
      for (int js2 = js1+1; js2 < num_simplices; js2++) {

        if (intersect_triangle_triangle
            (js1, js2, selfI_epsilon, intersection_point)) {
          output_triangle_triangle_intersection(js1, js2, intersection_point);
          flag_intersect = true;
        }
      }
    }
  }

  if (!flag_intersect) {
    cout << "Surface has no self intersections." << endl;
  }

  return(flag_intersect);
}

bool output_self_intersections_using_grid_of_bins()
{
  INTEGER_LIST<int, int> triList(num_simplices);
  vector<int> sorted_list;
  COORD_TYPE_PTR w0, w1, w2;
  COORD_TYPE intersection_point[DIM3];
  int min_grid_coord[DIM3];
  int max_grid_coord[DIM3];
  int binCoord[DIM3];
  const int NUM_SIMPLICES_PER_TIC = 10000;
  bool flag_output_tics = false;
  int num_tics = 0;

  if (num_simplices >= 5*NUM_SIMPLICES_PER_TIC)
    { flag_output_tics = true; }

  if (!is_num_bins_per_axis_set) {
    if (num_simplices > 1000) {
      if (num_simplices < 5000) 
        { num_bins_per_axis = 20; }
      else if (num_simplices < 10000)
        { num_bins_per_axis = 50; }
      else if (num_simplices < 50000)
        { num_bins_per_axis = 100; }
      else if (num_simplices < 100000)
        { num_bins_per_axis = 200; }
      else
        { num_bins_per_axis = 500; }

      cout << "Setting number of bins per axis to: "
           << num_bins_per_axis << "." << endl;
    }
  }

  GRID_OF_BINS_3D binGrid(num_bins_per_axis);

  compute_bounding_box();
  binGrid.SetMinCoord(bounding_box.MinCoord());
  binGrid.SetMaxCoord(bounding_box.MaxCoord());

  for (int js = 0; js < num_simplices; js++) {

    w0 = vertex_coord + DIM3*simplex_vert[js*DIM3];
    w1 = vertex_coord + DIM3*simplex_vert[js*DIM3+1];
    w2 = vertex_coord + DIM3*simplex_vert[js*DIM3+2];

    binGrid.InsertTri(w0, w1, w2, js);
  }

  bool flag_intersect = false;
  for (int js1 = 0; js1 < num_simplices; js1++) {

    if ((js1%NUM_SIMPLICES_PER_TIC) == 0 && flag_output_tics && 
        !flag_intersect) {
      cout << "."; 
      cout.flush();
      num_tics++;
    }


    triList.ClearList();

    w0 = vertex_coord + DIM3*simplex_vert[js1*DIM3];
    w1 = vertex_coord + DIM3*simplex_vert[js1*DIM3+1];
    w2 = vertex_coord + DIM3*simplex_vert[js1*DIM3+2];

    binGrid.ComputeMinMaxBinCoord(w0, w1, w2, min_grid_coord, max_grid_coord);

    for (binCoord[2] = min_grid_coord[2]; 
         binCoord[2] <= max_grid_coord[2]; binCoord[2]++) {
      for (binCoord[1] = min_grid_coord[1]; 
           binCoord[1] <= max_grid_coord[1]; binCoord[1]++) {
        binCoord[0] = min_grid_coord[0];
        int ibin = binGrid.ComputeBinIndex(binCoord);

        for (int x = min_grid_coord[0]; x <= max_grid_coord[0]; x++) {

          for (int k = 0; k < binGrid.Bin(ibin)->size(); k++) 
            { triList.Insert((*binGrid.Bin(ibin))[k]); }
          ibin++;
        }
      }
    }

    sorted_list.resize(triList.ListLength());
    for (int i = 0; i < triList.ListLength(); i++) 
      { sorted_list[i] = triList.List(i); }
    std::sort(sorted_list.begin(), sorted_list.end());

    for (int k = 0; k < sorted_list.size(); k++) {
      int js2 = sorted_list[k];
      if (js1 < js2) {
        if (intersect_triangle_triangle
            (js1, js2, selfI_epsilon, intersection_point)) {

          if (num_tics > 0) {
            cout << endl;
            num_tics = 0;
          }

          output_triangle_triangle_intersection(js1, js2, intersection_point);
          flag_intersect = true;
        }
      }
    }
  }
  
  if (flag_output_tics)
    { cout << endl; }

  return(flag_intersect);
}


void output_manifold_and_boundary_counts()
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  const int numv_per_facet = mesh_dimension;
  const int num_non_manifold_facets = 
    non_manifold_facet_vert.size()/numv_per_facet;
  const int num_deep_vertices = 
    count_deep_vertices(dimension, non_manifold_vert_list);

  cout << non_manifold_vert_list.size() << " "
       << num_deep_vertices << " "
       << num_non_manifold_facets << " "
       << num_internal_boundary_facets << " "
       << num_deep_boundary_facets << endl;
}


void write_non_manifold_edges()
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;
  const int num_vertices = mesh_data.num_vertices;
  const int numv_per_facet = mesh_dimension;
  const int num_non_manifold_edges =
    non_manifold_facet_vert.size()/numv_per_facet;

  if (mesh_dimension != DIM2) {
    cerr << "Only able to write non-manifold edges for mesh dimension two."
         << endl;
    return;
  }

  if (num_non_manifold_edges == 0) {

    cout << "No non-manifold edges." << endl;
    return;
  }

  ofstream output_file;
  output_file.open(output_filename, ios::out);

  float rgba[4] = { 1, 0, 1, 1};

  IJK::ijkoutColorLINE
    (output_file, dimension, vertex_coord, num_vertices,
     vector2pointer(non_manifold_facet_vert), num_non_manifold_edges, rgba);

  output_file.close();
}

void output_non_manifold_vertices()
{
  if (!terse_flag) {
    int num_output = 0;
    for (int i = 0; i < non_manifold_vert_list.size(); i++) {
      cout << "  " << non_manifold_vert_list[i];
      num_output++;
      if (num_output%10 == 0) { cout << endl; };
    }
    if (num_output%10 != 0) { cout << endl; };
  }
}

/// Return true if region from (0,...) to region_max[] contains point p[].
bool region_contains
(const int dimension, const COORD_TYPE * region_max, const COORD_TYPE * p)
{
  for (int d = 0; d < dimension; d++) {
    if (p[d] > region_max[d]) { return(false); }
  }

  return(true);
}

void output_vertex_list()
{
  const int dimension = mesh_data.dimension;
  const int num_vertices = mesh_data.num_vertices;

  for (int iv = 0; iv < num_vertices; iv++) {

    const COORD_TYPE * vcoord = vertex_coord + iv*dimension;

    if (is_min_coord_set) {
      if (!region_contains(dimension, vcoord, &(min_coord[0]))) {
        continue;
      }
    }

    if (is_max_coord_set) {
      if (!region_contains(dimension, &(max_coord[0]), vcoord)) {
        continue;
      }
    }

    cout << "Vertex " << iv << ": ";
    cout << "(";
    for (int ic = 0; ic < dimension; ic++) {
      cout << vertex_coord[iv*dimension+ic];
      if (ic+1 < dimension) { cout << ","; }
    }
    cout << ")";
    cout << endl;
    
  }
}

bool simplex_contains_vertex
(const int numv_per_simplex, const int * svert, const int iv)
{
  for (int k = 0; k < numv_per_simplex; k++) {
    if (svert[k] == iv) { return(true); }
  }

  return(false);
}

/// @param iend0 = Edge endpoint 0.
/// @param iend1 = Edge endpoint 1.
bool simplex_contains_edge
(const int numv_per_simplex, const int * svert, 
 const int iend0, const int iend1)
{
  bool contains_end0 = false;
  bool contains_end1 = false;
  for (int k = 0; k < numv_per_simplex; k++) {

    if (svert[k] == iend0) 
      { contains_end0 = true; }

    if (svert[k] == iend1) 
      { contains_end1 = true; }
  }

  return (contains_end0 && contains_end1);
}


bool poly_contains_vertex
(const int num_poly_vert, const int * poly_vert, const int iv)
{
  for (int k = 0; k < num_poly_vert; k++) {
    if (poly_vert[k] == iv) { return(true); }
  }

  return(false);
}

void output_simplices()
{
  const int num_vertices = mesh_data.num_vertices;

  if (contains_vertex_flag) {
    if (contains_vertex_index < 0 || contains_vertex_index >= num_vertices)
      {
        cout << "Illegal vertex index " << contains_vertex_index
             << ".  Vertex index should be in range["
             << 0 << "," << num_vertices-1 << "]." << endl;
        return;
      }
  }

  int num_out = 0;
  for (int js = 0; js < num_simplices; js++) {

    const int * svert = simplex_vert + js*num_vert_per_poly;

    if (contains_vertex_flag) {
      if (!simplex_contains_vertex
          (num_vert_per_poly, svert, contains_vertex_index))
        { continue; }
    }

    if (contains_edge_flag) {
      if (!simplex_contains_edge
          (num_vert_per_poly, svert, 
           edge_end0_index, edge_end1_index))
        { continue; }
    }

    num_out++;
    cout << "Simplex " << js << ": ";
    cout << "(";
    for (int k = 0; k < num_vert_per_poly; k++) {
      cout << svert[k];
      if (k+1 < num_vert_per_poly) { cout << ","; }
    }
    cout << ")";
    cout << endl;
  }

  if (num_out == 0 && !terse_flag) {

    if (contains_edge_flag && contains_vertex_flag) {
      cout << "No simplices contain vertex " 
           << contains_vertex_index << " and edge (" 
           << edge_end0_index << ", " << edge_end1_index
           << ")." << endl;
    }
    else if (contains_vertex_flag) {
      cout << "No simplices contain vertex " 
           << contains_vertex_index << "." << endl;
    }
    else if (contains_edge_flag) {
      cout << "No simplices contain edge (" 
           << edge_end0_index << ", " << edge_end1_index
           << ")." << endl;
    }

  }

}


void output_polytopes()
{
  const int num_vertices = mesh_data.num_vertices;

  if (contains_vertex_flag) {
    if (contains_vertex_index < 0 || contains_vertex_index >= num_vertices)
      {
        cout << "Illegal vertex index " << contains_vertex_index
             << ".  Vertex index should be in range["
             << 0 << "," << num_vertices-1 << "]." << endl;
        return;
      }
  }

  int num_out = 0;
  for (int jpoly = 0; jpoly < num_poly; jpoly++) {

    const int * pvert = polymesh.VertexList(jpoly);

    if (contains_vertex_flag) {
      if (!poly_contains_vertex
          (polymesh.NumPolyVert(jpoly), pvert, contains_vertex_index))
        { continue; }
    }

    if (is_min_num_polyv_output_set) {
      if (polymesh.NumPolyVert(jpoly) < min_num_polyv_output) { continue; }
    }

    if (is_max_num_polyv_output_set) {
      if (polymesh.NumPolyVert(jpoly) > max_num_polyv_output) { continue; }
    }

    num_out++;
    cout << "Poly " << jpoly << ": ";
    print_list(cout, pvert, polymesh.NumPolyVert(jpoly));
    cout << endl;
  }

  if (contains_vertex_flag) {
    if (num_out == 0 && !terse_flag) {

      if (contains_vertex_flag) {
        cout << "No polytopes contain vertex " 
             << contains_vertex_index << "." << endl;
      }
    }
  }

}


// **************************************************
// Write Tables
// **************************************************

template <class TABLE_TYPE>
void write_table_gplt(ofstream & ofile, const TABLE_TYPE & table)
{

  ofile.precision(DEFAULT_TABLE_PRECISION);
  ofile.setf(ios::left);

  if (flag_normalize) {
    ofile << "# normalized values" << endl;
  }
  ofile << "#";
  table.WriteColumnLabels(ofile, "  ");
  ofile << endl;

  int width = DEFAULT_TABLE_COLUMN_WIDTH;
  if (flag_normalize) {
    table.WriteNormalizedColumnData(ofile, "  ", width, 1);
  }
  else {
    table.WriteColumnData(ofile, "  ", width);
  }
}


template <class TABLE_TYPE>
void write_table_gplt(const string & filename, const TABLE_TYPE & table)
{
  ofstream ofile(filename.c_str(), ios::out);
  if (!flag_silent_write) {
    cout << "Writing table: " << filename << endl;
  }
  write_table_gplt(ofile, table);
  ofile.close();
}


template <class TABLE_TYPE>
void write_table_gplt(const string & filename_prefix, 
                      const string & filename_suffix,
                      const TABLE_TYPE & table)
{
  string filename = filename_prefix + "." + filename_suffix;
  write_table_gplt(filename, table);
}


void write_angle_table_gplt
(const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  ANGLE_TABLE angle_table;
  angle_table.angle.Include();
  angle_table.angle.SetAtIntervals(0, 1);

  compute_polygon_angles(angle_table);

  if (flag_min_table) {
    angle_table.HideAllExceptAngleColumn();
    angle_table.min_polygon_angle_freq.Show();

    write_table_gplt
      (output_filename_prefix, "min_poly_angle_freq.gplt", angle_table);
  }

  if (flag_max_table) {
    angle_table.HideAllExceptAngleColumn();
    angle_table.max_polygon_angle_freq.Show();

    write_table_gplt
      (output_filename_prefix, "max_poly_angle_freq.gplt", angle_table);
  }
}


void write_edge_length_table_gplt
(const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  EDGE_LENGTH_TABLE edge_length_table;
  edge_length_table.edge_length.Include();
  edge_length_table.edge_length.SetAtIntervals(0, edge_interval);

  compute_hex_edge_lengths
    (mesh_data, polymesh, vertex_coord, edge_interval, edge_length_table);

  if (flag_min_table) {
    edge_length_table.HideAllExceptEdgeLengthColumn();
    edge_length_table.min_edge_length_freq.Show();

    write_table_gplt
      (output_filename_prefix, "min_edge_length_freq.gplt", 
       edge_length_table);
  }

  if (flag_max_table) {
    edge_length_table.HideAllExceptEdgeLengthColumn();
    edge_length_table.max_edge_length_freq.Show();

    write_table_gplt
      (output_filename_prefix, "max_edge_length_freq.gplt", 
       edge_length_table);
  }
}



 // Determine min (starting) value of Jacobian determinants
//    in the Jacobian determinant table.
void determine_min_table_Jacobian
(const COORD_TYPE min_Jacobian, COORD_TYPE & min_table_Jacobian)
{
  min_table_Jacobian = 0.0;

  if (min_Jacobian < -0.5) 
    { min_table_Jacobian = -1.0; }
  else if (min_Jacobian < -0.2) 
    { min_table_Jacobian = -0.5; }
  else if (min_Jacobian < -0.1) 
    { min_table_Jacobian = -0.2; }
  else if (min_Jacobian < 0.0) 
    { min_table_Jacobian = -0.1; }
}


// Determine min value of Jacobian determinants in table
//   and number of table rows.
void determine_jacobian_table_parameters
(const COORD_TYPE min_jdet, const COORD_TYPE max_jdet,
 const COORD_TYPE table_interval,
 COORD_TYPE & min_table_Jacobian,
 int & num_rows)
{
  using namespace IJKDATATABLE;

  determine_min_table_Jacobian(min_jdet, min_table_Jacobian);

  if (max_jdet > 2.0) {
    num_rows = compute_num_buckets<int>
      (min_table_Jacobian, 3.0, table_interval); 
  }
  else if (max_jdet > 1.0) {
    num_rows = compute_num_buckets<int>
      (min_table_Jacobian, 3.0, table_interval); 
  }
  else {
    num_rows = compute_num_buckets<int>
      (min_table_Jacobian, 1.0, table_interval); 
  }

}


// Determine min value of normalized Jacobian determinant in table
//   and number of table rows.
void determine_normalized_jacobian_table_parameters
(const COORD_TYPE min_normalized_jdet,
 const COORD_TYPE max_normalized_jdet,
 const COORD_TYPE table_interval,
 COORD_TYPE & min_table_normalized_jdet,
 int & num_rows)
{
  using namespace IJKDATATABLE;

  determine_min_table_Jacobian
    (min_normalized_jdet, min_table_normalized_jdet);

  if (min_table_normalized_jdet < -1.0) { 
    // normalized Jacobian determinant should be in range [-1,1];
    min_table_normalized_jdet = -1.0;
  }

  num_rows = compute_num_buckets<int>
    (min_table_normalized_jdet, 1.0, table_interval); 
}


void write_min_jacobian_table_gplt
(const string & filename_prefix, 
 const string & filename_suffix,
 JACOBIAN_TABLE & jacobian_table)
{
  jacobian_table.HideAllExceptJacobianColumn();
  jacobian_table.min_jacobian_freq.Show();

  write_table_gplt
    (filename_prefix, "min_"+filename_suffix, jacobian_table);
}


void write_max_jacobian_table_gplt
(const string & filename_prefix, 
 const string & filename_suffix,
 JACOBIAN_TABLE & jacobian_table)
{
  jacobian_table.HideAllExceptJacobianColumn();
  jacobian_table.max_jacobian_freq.Show();

  write_table_gplt
    (filename_prefix, "max_"+filename_suffix, jacobian_table);
}


void write_jacobian_table_gplt
(const string & filename_prefix, 
 const string & filename_suffix,
 const bool flag_min_table, const bool flag_max_table,
 JACOBIAN_TABLE & jacobian_table)
{
  if (flag_min_table) {
    write_min_jacobian_table_gplt
      (filename_prefix, filename_suffix, jacobian_table);
  }

  if (flag_max_table) {
    write_max_jacobian_table_gplt
      (filename_prefix, filename_suffix, jacobian_table);
  }
}


void write_jacobian_poly_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  std::string filename_suffix;
  COORD_TYPE min_Jacobian, max_Jacobian;
  COORD_TYPE min_table_Jacobian;
  int num_rows;

  compute_min_max_hexahedra_Jacobian_determinants
    (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
     min_Jacobian, max_Jacobian);

  determine_jacobian_table_parameters
    (min_Jacobian, max_Jacobian, Jacobian_table_interval, 
     min_table_Jacobian, num_rows);

  JACOBIAN_TABLE jacobian_table(num_rows);
  jacobian_table.jacobian.Include();
  jacobian_table.jacobian.SetAtIntervals
    (min_table_Jacobian, Jacobian_table_interval);

  compute_hex_Jacobian_determinants
    (mesh_data, polymesh, vertex_coord, flag_internal_poly,
     min_table_Jacobian, Jacobian_table_interval, jacobian_table);

  filename_suffix = ".gplt";
  if (flag_internal_poly)
    { filename_suffix = ".internalP" + filename_suffix; }
  filename_suffix = "jdet_freq" + filename_suffix;

  write_jacobian_table_gplt
    (output_filename_prefix, filename_suffix, 
     flag_min_table, flag_max_table, jacobian_table);
}


void write_jacobian_vert_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  std::string filename_suffix;
  COORD_TYPE min_Jacobian, max_Jacobian;
  COORD_TYPE min_table_Jacobian;
  int num_rows;

  compute_min_max_hex_vert_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     min_Jacobian, max_Jacobian);

  determine_jacobian_table_parameters
    (min_Jacobian, max_Jacobian, Jacobian_table_interval, 
     min_table_Jacobian, num_rows);

  JACOBIAN_TABLE jacobian_table(num_rows);
  jacobian_table.jacobian.Include();
  jacobian_table.jacobian.SetAtIntervals
    (min_table_Jacobian, Jacobian_table_interval);

  compute_hex_vert_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert, 
     min_table_Jacobian, Jacobian_table_interval, jacobian_table);

  filename_suffix = ".gplt";
  if (flag_internal_vert) {
    if (flag_internal_poly) 
      { filename_suffix = ".internalVP" + filename_suffix; }
    else
      { filename_suffix = ".internalV" + filename_suffix; }
  }
  else {
    filename_suffix = ".V" + filename_suffix;
  }
  filename_suffix = "jdet_freq" + filename_suffix;

  write_jacobian_table_gplt
    (output_filename_prefix, filename_suffix, 
     flag_min_table, flag_max_table, jacobian_table);
}


void write_jacobian_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  if (flag_internal_vert) {
    write_jacobian_vert_table_gplt
      (mesh_data, polymesh, vertex_poly_incidence,
       output_filename_prefix, flag_min_table, flag_max_table);
  }
  else {
    if (io_info.flag_pJacobian.IsSetAndTrue()) {
      write_jacobian_poly_table_gplt
        (mesh_data, polymesh, vertex_poly_incidence,
         output_filename_prefix, flag_min_table, flag_max_table);
    }

    if (io_info.flag_vJacobian.IsSetAndTrue()) {
      write_jacobian_vert_table_gplt
        (mesh_data, polymesh, vertex_poly_incidence,
         output_filename_prefix, flag_min_table, flag_max_table);
    }
  }
}


void write_normalized_jacobian_poly_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  COORD_TYPE min_Jacobian, max_Jacobian;
  COORD_TYPE min_table_Jacobian;
  int num_rows;

  compute_min_max_hexahedra_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
     min_Jacobian, max_Jacobian);

  determine_normalized_jacobian_table_parameters
    (min_Jacobian, max_Jacobian, Jacobian_table_interval, 
     min_table_Jacobian, num_rows);

  JACOBIAN_TABLE jacobian_table(num_rows);
  jacobian_table.jacobian.Include();
  jacobian_table.jacobian.SetAtIntervals
    (min_table_Jacobian, Jacobian_table_interval);
  std::string filename_suffix;

  compute_hex_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_coord, flag_internal_poly,
     min_table_Jacobian, Jacobian_table_interval, jacobian_table);

  filename_suffix = ".gplt";
  if (flag_internal_poly)
    { filename_suffix = ".internalP" + filename_suffix; }
  filename_suffix = "normalized_jdet_freq" + filename_suffix;

  write_jacobian_table_gplt
    (output_filename_prefix, filename_suffix, 
     flag_min_table, flag_max_table, jacobian_table);
}


void write_normalized_jacobian_vert_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  COORD_TYPE min_Jacobian, max_Jacobian;
  COORD_TYPE min_table_Jacobian;
  int num_rows;

  compute_min_max_hexahedra_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
     min_Jacobian, max_Jacobian);

  determine_normalized_jacobian_table_parameters
    (min_Jacobian, max_Jacobian, Jacobian_table_interval, 
     min_table_Jacobian, num_rows);

  JACOBIAN_TABLE jacobian_table(num_rows);
  jacobian_table.jacobian.Include();
  jacobian_table.jacobian.SetAtIntervals
    (min_table_Jacobian, Jacobian_table_interval);
  std::string filename_suffix;

  compute_hex_vert_normalized_Jacobian_determinants
    (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert, 
     min_table_Jacobian, Jacobian_table_interval, jacobian_table);

  filename_suffix = ".gplt";
  if (flag_internal_vert) {
    if (flag_internal_poly) 
      { filename_suffix = ".internalVP" + filename_suffix; }
    else
      { filename_suffix = ".internalV" + filename_suffix; }
  }
  else {
    filename_suffix = ".V" + filename_suffix;
  }
  filename_suffix = "normalized_jdet_freq" + filename_suffix;

  write_jacobian_table_gplt
    (output_filename_prefix, filename_suffix, 
     flag_min_table, flag_max_table, jacobian_table);
}


void write_normalized_jacobian_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  if (flag_internal_vert) {
    write_normalized_jacobian_vert_table_gplt
      (mesh_data, polymesh, vertex_poly_incidence,
       output_filename_prefix, flag_min_table, flag_max_table);
  }
  else {
    if (io_info.flag_pJacobian.IsSetAndTrue()) {
      write_normalized_jacobian_poly_table_gplt
        (mesh_data, polymesh, vertex_poly_incidence,
         output_filename_prefix, flag_min_table, flag_max_table);
    }

    if (io_info.flag_vJacobian.IsSetAndTrue()) {
      write_normalized_jacobian_vert_table_gplt
        (mesh_data, polymesh, vertex_poly_incidence,
         output_filename_prefix, flag_min_table, flag_max_table);
    }
  }
}


void write_jacobian_shape_table_gplt
(const MESH_DATA & mesh_data,
 const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const std::string & output_filename_prefix,
 const bool flag_min_table, const bool flag_max_table)
{
  COORD_TYPE min_Jacobian_shape, max_Jacobian_shape;
  const COORD_TYPE min_table_value = 0;
  int num_rows;

  if (flag_internal_vert) {
    compute_min_max_hex_vert_Jacobian_shape
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       flag_internal_poly, flag_internal_vert,
       min_Jacobian_shape, max_Jacobian_shape);
  }
  else {
    compute_min_max_hexahedra_Jacobian_shape
      (mesh_data, polymesh, vertex_coord, flag_internal_poly, 
       min_Jacobian_shape, max_Jacobian_shape);
  }

  JACOBIAN_SHAPE_TABLE jacobian_shape_table;
  jacobian_shape_table.jshape.Include();
  jacobian_shape_table.jshape.SetAtIntervals(0, Jacobian_table_interval);
  std::string filename_suffix;

  if (flag_internal_vert) {
    compute_hex_vert_Jacobian_shape
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, 
       flag_internal_poly, flag_internal_vert, 
       min_table_value, Jacobian_table_interval, jacobian_shape_table);
  }
  else {
    compute_hex_Jacobian_shape
      (mesh_data, polymesh, vertex_coord, flag_internal_poly,
       min_table_value, Jacobian_table_interval, jacobian_shape_table);
  }

  filename_suffix = ".gplt";
  if (flag_internal_vert) {
    if (flag_internal_poly) 
      { filename_suffix = ".internalVP" + filename_suffix; }
    else
      { filename_suffix = ".internalV" + filename_suffix; }
  }
  else if (flag_internal_poly)
    { filename_suffix = ".internalP" + filename_suffix; }
  filename_suffix = "jshape_freq" + filename_suffix;

  if (flag_min_table) {
    jacobian_shape_table.HideAllExceptShapeColumn();
    jacobian_shape_table.min_jshape_freq.Show();

    write_table_gplt
      (output_filename_prefix, "min_"+filename_suffix, 
       jacobian_shape_table);
  }

  if (flag_max_table) {
    jacobian_shape_table.HideAllExceptShapeColumn();
    jacobian_shape_table.max_jshape_freq.Show();

    write_table_gplt
      (output_filename_prefix, "max_"+filename_suffix, jacobian_shape_table);
  }

}


// **************************************************
// COMPUTE INFO ROUTINES
// **************************************************

void compute_bounding_box()
{
  const int dimension = mesh_data.dimension;
  const int num_vertices = mesh_data.num_vertices;

  bounding_box.SetDimension(dimension);

  if (num_vertices < 1) {
    bounding_box.SetAllMinCoord(0);
    bounding_box.SetAllMaxCoord(0);
    return;
  }

  for (int ic = 0; ic < dimension; ic++) {
    COORD_TYPE minc = vertex_coord[ic];
    COORD_TYPE maxc = vertex_coord[ic];

    for (int iv = 1; iv < num_vertices; iv++) {
      COORD_TYPE c = vertex_coord[iv*dimension+ic];
      if (c < minc) { minc = c; };
      if (c > maxc) { maxc = c; };
    }
    
    bounding_box.SetMinCoord(ic, minc);
    bounding_box.SetMaxCoord(ic, maxc);
  }

  contracted_bounding_box = bounding_box;
  for (int d = 0; d < dimension; d++) {
    COORD_TYPE minc = bounding_box.MinCoord(d);
    COORD_TYPE maxc = bounding_box.MaxCoord(d);

    minc += contract_margin;
    maxc -= contract_margin;

    if (minc < maxc) {
      contracted_bounding_box.SetMinCoord(d, minc);
      contracted_bounding_box.SetMaxCoord(d, maxc);
    }
    else {
      flag_small_bounding_box = true;
    }
  }
}

/// Compute number of edges.
int compute_num_edges()
{
  const int numv_per_simplex = num_vert_per_poly;
  const int nume_per_simplex = numv_per_simplex*(numv_per_simplex-1)/2;
  const int max_elist_length = 2*nume_per_simplex*num_simplices;

  // Store edges in elist.
  IJK::ARRAY<int> elist(max_elist_length);
  int elist_length = 0;
  for (int js = 0; js < num_simplices; js++) {
    int k = js * numv_per_simplex;
    for (int i0 = 0; i0+1 < numv_per_simplex; i0++) {
      for (int i1 = i0+1; i1 < numv_per_simplex; i1++) {
        int iv0 = simplex_vert[k+i0];
        int iv1 = simplex_vert[k+i1];
        if (iv0 > iv1) { std::swap(iv0,iv1); };
        if (iv0 != iv1) {
          elist[elist_length] = iv0;
          elist[elist_length+1] = iv1;
          elist_length += 2;
        }
      }
    }
  }

  IJK::ARRAY<int> elist_loc(max_elist_length);
  std::vector<int> elist_nodup;

  merge_pairs
    (elist.PtrConst(), elist_length/2, elist_nodup, elist_loc.Ptr());

  return(elist_nodup.size()/2);
}

// Compute facet information
void compute_facet_info()
{
  const int mesh_dimension = mesh_data.mesh_dimension;
  const int NUMV_PER_CUBE = 8;
  const int NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;
  const int DIM3(3);

  compute_bounding_box();
  if (flag_simplex_file || mesh_dimension == 2) {
    if (!mesh_data.are_non_manifold_facets_identified || 
        !mesh_data.are_boundary_facets_identified) 
      { identify_non_manifold_and_boundary_facets(); }

    const int num_vert_per_facet = mesh_dimension;
    mesh_info.num_non_manifold_facets =
      (non_manifold_facet_vert.size())/num_vert_per_facet;

    set_boundary_vertices(boundary_facet_vert, polymesh);
  }
  else if (flag_cube_file) {

    if (!mesh_data.are_non_manifold_facets_identified || 
        !mesh_data.are_boundary_facets_identified) 
      { identify_non_manifold_and_boundary_facets(); }

    mesh_info.num_non_manifold_facets =
      (non_manifold_facet_vert.size())/NUMV_PER_CUBE_FACET;

    set_boundary_vertices(boundary_facet_vert, polymesh);
  }

  if (mesh_data.dimension == DIM3 && flag_cube_file) {
    identify_hex_sharing_exactly_two_facet_edges
      (polymesh, hex_facet_pairs_sharing_exactly_two_edges);
    mesh_info.num_hex_facet_pairs_sharing_exactly_two_edges =
      hex_facet_pairs_sharing_exactly_two_edges.size();
  }

  if (oriented_manifold_flag) {
    mesh_info.num_poly_with_orientation_conflicts =
      orientation_conflict_list.size();
  }
}


// **************************************************
// MESH PROCESSING ROUTINES
// **************************************************

void sort_poly()
{
  polymesh_sorted.Copy(polymesh);
  polymesh_sorted.SortVert();
  polymesh_sorted.GetSortedPolytopeIndices(sorted_poly);
}

// Count number of poly with num_poly_vert vertices
int count_num_poly(const int num_poly_vert)
{
  int n = 0;
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
    if (polymesh.NumPolyVert(ipoly) == num_poly_vert) 
      { n++; }
  }

  return(n);
}

// Count number of internal mesh vertices.
// @pre set_boundary_vertices must be called before calling this routine.
int count_num_internal_vertices(const POLYMESH_TYPE & polymesh)
{
  int n = 0;
  for (int iv = 0; iv < polymesh.NumVertices(); iv++) {
    if (polymesh.vertex_data[iv].IsInternal()) 
      { n++; }
  }

  return(n);
}

// **************************************************
// IDENTIFY DUPLICATE POLY OR POLY VERTICES
// **************************************************

void identify_duplicates()
{
  mesh_info.num_poly_with_duplicate_vertices =
    identify_poly_with_duplicate_vertices();
  mesh_info.num_duplicate_poly =
    identify_duplicate_polytopes();
}

// Identify polytopes with duplicate vertices.
// Return number of poly with duplicate vertices.
// @pre polymesh_sorted is created and its vertices are sorted.
int identify_poly_with_duplicate_vertices()
{
  int num_poly_with_duplicate_vertices = 0;

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {
    if (polymesh_sorted.NumPolyVert(ipoly) < 2) { continue; }

    int iv0 = polymesh_sorted.Vertex(ipoly, 0);
    for (int j = 1; j < polymesh_sorted.NumPolyVert(ipoly); j++) {
      int iv1 = polymesh_sorted.Vertex(ipoly, j);
      if (iv0 == iv1) {
        polymesh.poly_data[ipoly].is_degenerate = true;
        num_poly_with_duplicate_vertices++;
        break;
      }
      iv0 = iv1;
    }
  }

  return(num_poly_with_duplicate_vertices);
}

// Identify polytopes with duplicate vertices.
// Return number of duplicate polytopes.
// @pre polymesh_sorted is created and its vertices are sorted.
int identify_duplicate_polytopes()
{
  POLYMESH_LESS_THAN<POLYMESH_TYPE> polymesh_lt(&polymesh_sorted);

  int num_duplicate_poly = 0;

  if (sorted_poly.size() != polymesh_sorted.NumPoly()) 
    { polymesh_sorted.GetSortedPolytopeIndices(sorted_poly); }

  for (int i0 = 0; i0+1 < sorted_poly.size(); i0++) {
    int ipoly0 = sorted_poly[i0];
    int ipoly1 = sorted_poly[i0+1];

    if (!polymesh_lt(ipoly0,ipoly1) && !polymesh_lt(ipoly1,ipoly0)) {
      polymesh.poly_data[ipoly0].is_duplicate = true;
      polymesh.poly_data[ipoly1].is_duplicate = true;
    }
  }

  for (int i = 0; i < polymesh.NumPoly(); i++) {
    if (polymesh.poly_data[i].IsDuplicate()) { num_duplicate_poly++; }
  }


  return(num_duplicate_poly);
}


// **************************************************
// MANIFOLD ROUTINES
// **************************************************

void identify_non_manifold()
{
  const int dimension = mesh_data.dimension;
  const int mesh_dimension = mesh_data.mesh_dimension;

  if (flag_simplex_file || mesh_dimension == 2 || flag_cube_file) {
    compute_facet_info();

    if (flag_cube_file && mesh_dimension > 2) {
      mesh_info.num_non_manifold_edges = 
        identify_non_manifold_edges(); 
    }

    mesh_info.num_non_manifold_vertices = 
      identify_non_manifold_vertices();
    mesh_info.num_deep_non_manifold_vertices = 
      count_deep_vertices(dimension, non_manifold_vert_list);
  }
}


// Return true if jf0 equals jf1.
// @pre facet vertices are listed in sorted order.
bool facet_equals
(const int * facet_vert, const int numv_per_facet, const int jf0, const int jf1)
{
  return(are_lists_equal
         (facet_vert+jf0*numv_per_facet, facet_vert+jf1*numv_per_facet, 2));
}


void set_in_non_manifold_facet()
{
  const int num_vertices = mesh_data.num_vertices;

  // initialize in_non_manifold_facet[iv] to false for all vertices iv
  in_non_manifold_facet.assign(num_vertices, false);

  for (int i = 0; i < non_manifold_facet_vert.size(); i++) {
    int iv = non_manifold_facet_vert[i];
    in_non_manifold_facet[iv] = true;
  }
}


void set_in_non_manifold_edge()
{
  const int num_vertices = mesh_data.num_vertices;

  // initialize in_non_manifold_edge[iv] to false for all vertices iv
  in_non_manifold_edge.assign(num_vertices, false);

  for (int i = 0; i < non_manifold_edge_vert.size(); i++) {
    int iv = non_manifold_edge_vert[i];
    in_non_manifold_edge[iv] = true;
  }
}


void set_non_manifold_vert
(const vector<int> & non_manifold_vert_list, 
 vector<bool> & non_manifold_vert)
{
  const int num_vertices = mesh_data.num_vertices;

  // Initialize non_manifold_vert[iv] to false for all vertices iv.
  non_manifold_vert.assign(num_vertices, false);

  for (int i = 0; i < non_manifold_vert_list.size(); i++) {
    int iv = non_manifold_vert_list[i];
    non_manifold_vert[iv] = true;
  }
}


void set_poly_containing_non_manifold_facets
(const FACET_INFO_ARRAY & non_manifold_facets, POLYMESH_TYPE & polymesh)
{
  // Flag polytopes containing non-manifold facets.
  for (int j = 0; j < non_manifold_facets.size(); j++) {
    const int jpoly = non_manifold_facets[j].poly_containing_face;
    polymesh.poly_data[jpoly].contains_non_manifold_facet = true;
  }
}


void set_poly_containing_boundary_facets
(const FACET_INFO_ARRAY & boundary_facets, POLYMESH_TYPE & polymesh)
{
  // Flag polytopes containing boundary facets.
  for (int j = 0; j < boundary_facets.size(); j++) {
    const int jpoly = boundary_facets[j].poly_containing_face;
    polymesh.poly_data[jpoly].contains_boundary_facet = true;
  }
}


void set_poly_containing_mismatched_orientations
(const FACET_INFO_ARRAY & orientation_mismatch_facets, 
 POLYMESH_TYPE & polymesh, vector<int> & orientation_conflict_list)
{
  // Flag polytopes containing orientation mismatch facets.
  for (int j = 0; j < orientation_mismatch_facets.size(); j++) {
    const int jpoly = orientation_mismatch_facets[j].poly_containing_face;
    polymesh.poly_data[jpoly].orientation_conflict = true;
  }

  // Get list of polytopes with mismatched orientations.
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
    if (polymesh.poly_data[ipoly].orientation_conflict) 
      { orientation_conflict_list.push_back(ipoly); }
  }
}

void determine_internal_boundary_facets
(const vector<int> & boundary_facet_vert,
 const int num_vert_per_facet,
 vector<bool> & internal_boundary_facet,
 vector<bool> & far_from_bounding_box)
{
  internal_boundary_facet.clear();
  far_from_bounding_box.clear();

  // Determine internal boundary facets;
  const int num_boundary_facets = 
    boundary_facet_vert.size()/num_vert_per_facet;

  for (int j = 0; j < num_boundary_facets; j++) {
    bool flag_internal_poly =
      is_internal
      (bounding_box, boundary_facet_vert, num_vert_per_facet, j);
    internal_boundary_facet.push_back(flag_internal_poly);
    if (flag_internal_poly) { num_internal_boundary_facets++; }

    if (flag_small_bounding_box) 
      { flag_internal_poly = false; }
    else {
      flag_internal_poly =
        is_internal(contracted_bounding_box, boundary_facet_vert,
                    num_vert_per_facet, j);
    }
    far_from_bounding_box.push_back(flag_internal_poly);
    if (flag_internal_poly) { num_deep_boundary_facets++; }
  }
}

void set_skip_degen_dup_poly_flag
(const POLYMESH_TYPE & polymesh, std::vector<bool> & flag_skip_poly)
{
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
    if (polymesh.poly_data[ipoly].IsDegenerate() ||
        polymesh.poly_data[ipoly].IsDuplicate()) 
      { flag_skip_poly[ipoly] = true; }
  }
}

// Set flagB to equal flagA or flagB.
void set_flag_or_flag
(const std::vector<bool> & flagA,
 std::vector<bool> & flagB)
{
  if (flagB.size() != flagA.size()) 
    { flagB.resize(flagA.size(), false); }

  for (int i = 0; i < flagB.size(); i++) {
    if (flagA[i]) { flagB[i] = true; }
  }
}


// Identify non-manifold edges and boundary edges of 2D mesh.
// An edge is non-manifold if it is in more than two polygons.
// An edge is boundary if it is in only one polygon.
void identify_non_manifold_and_boundary_edges_of_mesh2D()
{
  const int NUM_VERT_PER_EDGE(2);
  FACET_INFO_ARRAY non_manifold_edges;
  FACET_INFO_ARRAY boundary_edges;
  FACET_INFO_ARRAY orientation_mismatch_edges;
  std::vector<bool> flag_skip_poly(polymesh.NumPoly(), false);

  set_skip_degen_dup_poly_flag(polymesh, flag_skip_poly);

  get_non_manifold_and_boundary_edges_of_mesh2D
    (polymesh, flag_skip_poly, non_manifold_facet_vert, non_manifold_edges, 
     boundary_facet_vert, boundary_edges,
     orientation_mismatch_facet_vert, orientation_mismatch_edges);

  set_poly_containing_non_manifold_facets(non_manifold_edges, polymesh);
  set_poly_containing_boundary_facets(boundary_edges, polymesh);
  set_poly_containing_mismatched_orientations
    (orientation_mismatch_edges, polymesh, orientation_conflict_list);

  determine_internal_boundary_facets
    (boundary_facet_vert, NUM_VERT_PER_EDGE, internal_boundary_facet,
     far_from_bounding_box);
}


// Identify non-manifold facets and boundary facets of tet mesh.
// A facet is non-manifold if it is in more than two polytope.
// A facet is boundary if it is in only one polytope.
void identify_non_manifold_and_boundary_facets_of_tet_mesh()
{
  const int NUM_VERT_PER_TET(4);
  const int NUM_VERT_PER_FACET(NUM_VERT_PER_TET-1);
  FACET_INFO_ARRAY non_manifold_facets;
  FACET_INFO_ARRAY boundary_facets;
  FACET_INFO_ARRAY orientation_mismatch_facets;

  get_non_manifold_and_boundary_facets_of_tet_mesh
    (polymesh, non_manifold_facet_vert, non_manifold_facets, 
     boundary_facet_vert, boundary_facets,
     orientation_mismatch_facet_vert, orientation_mismatch_facets);

  set_poly_containing_non_manifold_facets(non_manifold_facets, polymesh);
  set_poly_containing_boundary_facets(boundary_facets, polymesh);
  set_poly_containing_mismatched_orientations
    (orientation_mismatch_facets, polymesh, orientation_conflict_list);

  determine_internal_boundary_facets
    (boundary_facet_vert, NUM_VERT_PER_FACET, internal_boundary_facet,
     far_from_bounding_box);
}


// Identify non-manifold facets and boundary facets of simplicial mesh.
// A facet is non-manifold if it is in more than two polytope.
// A facet is boundary if it is in only one polytope.
void identify_non_manifold_and_boundary_facets_of_simplicial_mesh()
{
  const int NUM_VERT_PER_TET(4);
  const int NUM_VERT_PER_FACET(NUM_VERT_PER_TET-1);
  FACET_INFO_ARRAY non_manifold_facets;
  FACET_INFO_ARRAY boundary_facets;
  FACET_INFO_ARRAY orientation_mismatch_facets;

  get_non_manifold_and_boundary_facets_of_simplicial_mesh
    (polymesh, non_manifold_facet_vert, non_manifold_facets, 
     boundary_facet_vert, boundary_facets,
     orientation_mismatch_facet_vert, orientation_mismatch_facets);

  set_poly_containing_non_manifold_facets(non_manifold_facets, polymesh);
  set_poly_containing_boundary_facets(boundary_facets, polymesh);
  set_poly_containing_mismatched_orientations
    (orientation_mismatch_facets, polymesh, orientation_conflict_list);

  determine_internal_boundary_facets
    (boundary_facet_vert, NUM_VERT_PER_FACET, internal_boundary_facet,
     far_from_bounding_box);
}


// Identify non-manifold facets and boundary facets of hex mesh.
// A facet is non-manifold if it is in more than two polytope.
// A facet is boundary if it is in only one polytope.
void identify_non_manifold_and_boundary_facets_of_hex_mesh()
{
  const CUBE_TYPE cube(DIM3);
  const int NUM_VERT_PER_FACET(cube.NumFacetVertices());
  FACET_INFO_ARRAY non_manifold_facets;
  FACET_INFO_ARRAY boundary_facets;
  FACET_INFO_ARRAY orientation_mismatch_facets;

  get_non_manifold_and_boundary_facets_of_hex_mesh
    (polymesh, cube, non_manifold_facet_vert, non_manifold_facets, 
     boundary_facet_vert, boundary_facets,
     orientation_mismatch_facet_vert, orientation_mismatch_facets);

  set_poly_containing_non_manifold_facets(non_manifold_facets, polymesh);
  set_poly_containing_boundary_facets(boundary_facets, polymesh);
  set_poly_containing_mismatched_orientations
    (orientation_mismatch_facets, polymesh, orientation_conflict_list);

  determine_internal_boundary_facets
    (boundary_facet_vert, NUM_VERT_PER_FACET, internal_boundary_facet,
     far_from_bounding_box);
}


// Identify non-manifold facets and boundary facets.
// A facet is non-manifold if it is in more than two polytope.
// A facet is boundary if it is in only one polytope.
void identify_non_manifold_and_boundary_facets()
{
  const int DIM3(3);
  const int mesh_dimension = mesh_data.mesh_dimension;
  IJK::PROCEDURE_ERROR error
    ("identify_non_manifold_and_boundary_facets");

  if (mesh_dimension == 2) 
    { identify_non_manifold_and_boundary_edges_of_mesh2D(); }
  else if (flag_simplex_file) {
    if (mesh_dimension == DIM3) 
      { identify_non_manifold_and_boundary_facets_of_tet_mesh(); }
    else 
      { identify_non_manifold_and_boundary_facets_of_simplicial_mesh(); }
  }
  else if (flag_cube_file) 
    { identify_non_manifold_and_boundary_facets_of_hex_mesh(); }
  else {
    error.AddMessage
      ("Programming error.  Unable to determine polytope type.");
    throw error;
  }

  mesh_data.are_boundary_facets_identified = true;
  mesh_data.are_non_manifold_facets_identified = true;
  set_in_non_manifold_facet();
}


// Return number of deep vertices in list
int count_deep_vertices
(const int dimension, const std::vector<int> & vlist)
{
  int count = 0;
  for (int i = 0; i < vlist.size(); i++) {
    int iv = vlist[i];
    COORD_TYPE * coord_ptr = vertex_coord+dimension*iv;
    if (contracted_bounding_box.Contains(coord_ptr)) 
      { count++; }
  }
  return(count);
}



bool are_poly_adjacent
(const POLYMESH_TYPE & polymesh, const int kpoly, const int jpoly)
// return true if polytopes kpoly and jpoly share facet
{
  const int mesh_dimension = mesh_data.mesh_dimension;
  const int NUMV_PER_CUBE = 8;
  static CUBE_TYPE cube(DIM3);
  IJK::PROCEDURE_ERROR error("are_poly_adjacent");

  if (mesh_dimension == 2) {
    return(are_poly2D_adjacent(polymesh, kpoly, jpoly));
  }
  else if (flag_simplex_file) {
    return(are_simplices_adjacent(polymesh, kpoly, jpoly));
  }
  else if (num_vert_per_poly == NUMV_PER_CUBE) {
    return(are_hexahedra_adjacent(polymesh, kpoly, jpoly, cube));
  }
  else {
    error.AddMessage
      ("Programming error. Polytopes are not simplices and not 2D.");
    error.AddMessage("Cannot determine if polytopes are adjacent.");
    throw error;
  }
}


/// Set on_boundary flag to tree for boundary polymesh vertices.
void set_boundary_vertices
(const vector<int> & boundary_facet_vert,
 POLYMESH_TYPE & polymesh)
{
  // Initialize on_boundary to false.
  for (int iv = 0; iv < polymesh.NumVertices(); iv++) 
    { polymesh.vertex_data[iv].on_boundary = false; }

  for (int i = 0; i < boundary_facet_vert.size(); i++) {
    const VERTEX_INDEX iv = boundary_facet_vert[i]; 
    polymesh.vertex_data[iv].on_boundary = true;
  }
}


/// Identify non-manifold vertices of 2D mesh.
int identify_non_manifold_vertices_of_mesh2D
  (const std::vector<bool> & flag_skip_vert)
{
  VERTEX_POLY_INCIDENCE_TYPE vertex_info(polymesh);

  get_non_manifold_vertices_of_mesh2D
    (polymesh, vertex_info, flag_skip_vert, non_manifold_vert_list);

  set_non_manifold_vert(non_manifold_vert_list, non_manifold_vert);

  return(non_manifold_vert_list.size());
}


/// Identify non-manifold vertices of simplicial mesh.
int identify_non_manifold_vertices_of_simplicial_mesh
  (const std::vector<bool> & flag_skip_vert)
{
  std::vector<bool> flag_skip_poly(polymesh.NumPoly(), false);
  VERTEX_POLY_INCIDENCE_TYPE vertex_info(polymesh);

  get_non_manifold_vertices_of_simplicial_mesh
    (polymesh, vertex_info, flag_skip_vert, non_manifold_vert_list);

  set_non_manifold_vert(non_manifold_vert_list, non_manifold_vert);

  return(non_manifold_vert_list.size());
}


/// Identify non-manifold vertices of tetrahedral mesh.
int identify_non_manifold_vertices_of_tet_mesh
  (const std::vector<bool> & flag_skip_vert)
{
  std::vector<bool> flag_skip_poly(polymesh.NumPoly(), false);
  VERTEX_POLY_INCIDENCE_TYPE vertex_info(polymesh);

  get_non_manifold_vertices_of_tet_mesh
    (polymesh, vertex_info, flag_skip_vert, 
     non_manifold_vert, non_manifold_vert_list);

  set_non_manifold_vert(non_manifold_vert_list, non_manifold_vert);

  return(non_manifold_vert_list.size());
}


/// Identify non-manifold vertices of hexahedral mesh.
int identify_non_manifold_vertices_of_hex_mesh
  (const std::vector<bool> & flag_skip_vert)
{
  const int DIM3(3);
  VERTEX_POLY_INCIDENCE_TYPE vertex_info(polymesh);
  CUBE_TYPE cube(DIM3);

  get_non_manifold_vertices_of_hex_mesh
    (polymesh, vertex_info, flag_skip_vert, cube, 
     non_manifold_vert, non_manifold_vert_list);

  set_non_manifold_vert(non_manifold_vert_list, non_manifold_vert);

  return(non_manifold_vert_list.size());
}


/// Identify non-manifold vertices.
int identify_non_manifold_vertices()
{
  const int DIM2(2);
  const int DIM3(3);
  const int mesh_dimension = mesh_data.mesh_dimension;
  const int num_vertices = mesh_data.num_vertices;
  std::vector<bool> flag_skip_vert(num_vertices, false);
  int num_non_manifold_vertices;
  IJK::PROCEDURE_ERROR error("identify_non_manifold_vertices");

  non_manifold_vert.assign(num_vertices, false);
  set_flag_or_flag(in_non_manifold_edge, flag_skip_vert);
  set_flag_or_flag(in_non_manifold_facet, flag_skip_vert);

  if (mesh_dimension == DIM2) {
    num_non_manifold_vertices =
      identify_non_manifold_vertices_of_mesh2D(flag_skip_vert);
  }
  else if (flag_simplex_file) {
    if (mesh_dimension == DIM3) {
      num_non_manifold_vertices =
        identify_non_manifold_vertices_of_tet_mesh(flag_skip_vert);
    }
    else {
      num_non_manifold_vertices =
        identify_non_manifold_vertices_of_simplicial_mesh(flag_skip_vert);
    }
  }
  else if (flag_cube_file) {
    num_non_manifold_vertices =
      identify_non_manifold_vertices_of_hex_mesh(flag_skip_vert);
  }
  else {
    error.AddMessage
      ("Programming error.  Unable to identify non_manifold vertices.");
    error.AddMessage("  Unknown mesh type.");
    throw error;
  }

  return(num_non_manifold_vertices);
}


/// Identify non-manifold edges
int identify_non_manifold_edges()
{
  const int num_vertices = mesh_data.num_vertices;
  const int NUMV_PER_CUBE = 8;
  const CUBE_TYPE cube(DIM3);
  std::vector<int> poly_list;
  int num_non_manifold_edges = 0;
  IJK::PROCEDURE_ERROR error("identify_non_manifold_edges");

  if (num_vert_per_poly != NUMV_PER_CUBE) {
      ("Programming error.  Mesh is not a mesh of hexahedra.");
    error.AddMessage
      ("  Non-manifold edges only implemented for mesh of hexahedra.");
    throw error;
  }

  non_manifold_edge_vert.clear();

  VERTEX_POLY_INCIDENCE<int,int> vertex_info(polymesh);
  VERTEX_ADJACENCY_LIST<int,int> adjacency_list;

  if (in_non_manifold_facet.size() != num_vertices) {
    error.AddMessage
      ("Programming error.  Array in_non_manifold_facet[] not set.");
    throw error;
  }

  adjacency_list.SetFromMeshOfCubes(polymesh, cube);

  for (int iv0 = 0; iv0 < adjacency_list.NumVertices(); iv0++) {

    if (in_non_manifold_facet[iv0]) { continue; }

    for (int i1 = 0; i1 < adjacency_list.NumAdjacent(iv0); i1++) {
      const int iv1 = adjacency_list.AdjacentVertex(iv0, i1);

      if (iv0 >= iv1) {
        // Skip if iv0 == iv1.
        // Skip if iv0 > iv1, since edge handled as neighbor of iv1.
        continue;
      }

      if (in_non_manifold_facet[iv1]) { continue; }

      vertex_info.GetPolyContaining(iv0, iv1, poly_list);
      const int num_incident_poly = poly_list.size();

      if (num_incident_poly == 0) { continue; }

      vector<bool> is_reachable(num_incident_poly, false);
      vector<int> reachable;
      reachable.push_back(0);

      while (reachable.size() > 0) {
        int n = reachable.size()-1;
        int j = reachable[n];
        reachable.pop_back();
        is_reachable[j] = true;
        int jpoly = poly_list[j];

        for (int k = 0; k < num_incident_poly; k++) {
          if (!is_reachable[k]) {
            int kpoly = poly_list[k];
            if (are_poly_adjacent(polymesh, kpoly, jpoly)) {
              reachable.push_back(k);
            }
          }
        }
      }

      for (int k = 0; k < num_incident_poly; k++) {
        if (!is_reachable[k]) {
          non_manifold_edge_vert.push_back(iv0);
          non_manifold_edge_vert.push_back(iv1);
          num_non_manifold_edges++;

          for (int k2 = 0; k2 < poly_list.size(); k2++) {
            const int kpoly = poly_list[k2];
            polymesh.poly_data[kpoly].contains_non_manifold_edge =
              true;
          }

          break;
        }
      }
    }
  }

  set_in_non_manifold_edge();

  return(num_non_manifold_edges);
}


// Return true if facet jf is internal
bool is_internal(const BOUNDING_BOX & bounding_box, 
                 const vector<int> & facet_vlist, 
                 const int numv_per_facet,
                 const int jf)
{
  const int dimension = mesh_data.dimension;

  for (int d = 0; d < dimension; d++) {
    COORD_TYPE minc = bounding_box.MinCoord(d);
    COORD_TYPE maxc = bounding_box.MaxCoord(d);

    bool greater_than_minc = false;
    bool less_than_maxc = false;
    for (int i = jf*numv_per_facet; i < (jf+1)*numv_per_facet; i++) {
      int iv = facet_vlist[i];
      COORD_TYPE c = vertex_coord[iv*dimension+d];
      if (c > minc) { greater_than_minc = true; };
      if (c < maxc) { less_than_maxc = true; };
    }
    if (!greater_than_minc || !less_than_maxc) { return(false); }
  }

  return(true);
}


// **************************************************
// POLYTOPE INTERSECTION ROUTNES
// **************************************************

namespace {

  /// Return true if first three vertices of vlistA and vlistB are equal.
  bool do_first_three_vertices_match
  (const VERTEX_INDEX vlistA[], const VERTEX_INDEX vlistB[]) {
    if (vlistA[0] == vlistB[0] &&
        vlistA[1] == vlistB[1] &&
        vlistA[2] == vlistB[2])
      { return(true); }

    return(false);
  }

}

// Identify hexahedra which share exactly two facet edges, but not a facet.
void identify_hex_sharing_exactly_two_facet_edges
(const POLYMESH_TYPE & polymesh, 
 FACET_INFO_PAIRS_ARRAY & facet_pairs_sharing_exactly_two_edges)
{
  const CUBE_TYPE cube(DIM3);
  const int NUM_VERT_PER_HEX_FACET(4);
  VERTEX_INDEX facet_vert[NUM_VERT_PER_HEX_FACET];
  VERTEX_INDEX ordered_facet_vert[NUM_VERT_PER_HEX_FACET];
  vector<int> index_sorted;

  // List of hex facets with each facet listed four times
  //   with four different last vertices.
  IJK::FACET_LIST_BASE<VERTEX_INDEX,int,FACET_INFO> ordered_facet;

  for (int ihex = 0; ihex < polymesh.NumPoly(); ihex++) {

    for (int jf = 0; jf < cube.NumFacets(); jf++) {
      facet_vert[0] = cube.FacetVertex(jf, 0);
      facet_vert[1] = cube.FacetVertex(jf, 1);

      // Reverse vertices 2 and 3 to list facets in clockwise or
      //   counter-clockwise order around facet.
      facet_vert[2] = cube.FacetVertex(jf, 3);
      facet_vert[3] = cube.FacetVertex(jf, 2);

      for (int k = 0; k < NUM_VERT_PER_HEX_FACET; k++) 
        { facet_vert[k] = polymesh.Vertex(ihex, facet_vert[k]); }

      for (int k0 = 0; k0 < NUM_VERT_PER_HEX_FACET; k0++) {
        for (int k1 = 0; k1 < NUM_VERT_PER_HEX_FACET; k1++) {
          int k2 = (k0+k1)%NUM_VERT_PER_HEX_FACET;
          ordered_facet_vert[k1] = facet_vert[k2];
        }

        if (ordered_facet_vert[0] > ordered_facet_vert[2]) {
          // Reverse facet orientation so that 
          //   ordered_facet_vert[0] <= ordered_facet_vert[2].
          std::swap(ordered_facet_vert[0], ordered_facet_vert[2]);
        }

        const int ifacet = ordered_facet.NumFacets();
        ordered_facet.AddPolytope(ordered_facet_vert, NUM_VERT_PER_HEX_FACET);
        ordered_facet.poly_data[ifacet].poly_containing_face = ihex;
        ordered_facet.poly_data[ifacet].face_index = jf;
      }
    }
  }

  index_sorted.resize(ordered_facet.NumFacets());
  ordered_facet.GetSortedPolytopeIndices(index_sorted);

  // Check for hexahedra sharing exactly two edges but no facets.
  for (int k1 = 1; k1 < ordered_facet.NumFacets(); k1++) {
    const int kf1 = index_sorted[k1];
    const int k0 = k1-1;
    const int kf0 = index_sorted[k0];

    if (do_first_three_vertices_match
        (ordered_facet.VertexList(kf0), ordered_facet.VertexList(kf1))) {

      if (ordered_facet.Vertex(kf0,3) != ordered_facet.Vertex(kf1,3)) {
        FACET_INFO facet0_info = ordered_facet.poly_data[kf0];
        FACET_INFO facet1_info = ordered_facet.poly_data[kf1];
        facet_pairs_sharing_exactly_two_edges.push_back
          (FACET_INFO_PAIR(facet0_info, facet1_info));
      }
    }

  }

  mesh_data.are_hex_sharing_exactly_two_facet_edges_identified = true;
}


// **************************************************
// COMPUTE TABLE ROUTINES
// **************************************************

void compute_polygon_angles(ANGLE_TABLE & angle_table)
{
  const int dimension = mesh_data.dimension;

  ANGLE_TYPE min_angle, max_angle;
  int num_angle;

  angle_table.min_polygon_angle_freq.Include();
  angle_table.min_polygon_angle_freq.SetAll(0);
  angle_table.max_polygon_angle_freq.Include();
  angle_table.max_polygon_angle_freq.SetAll(0);

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate())
      { continue; }

    IJK::compute_min_max_polygon_angles
      (dimension, polymesh.VertexList(ipoly), 
       polymesh.NumPolyVert(ipoly),
       vertex_coord, min_angle, max_angle, num_angle);

    if (num_angle < 3) { continue; }

    int imin = floor(min_angle+0.5);
    int imax = floor(max_angle+0.5);

    angle_table.min_polygon_angle_freq.Increment(imin);
    angle_table.max_polygon_angle_freq.Increment(imax);
  }
}


void compute_hex_edge_lengths
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const COORD_TYPE edge_interval,
 EDGE_LENGTH_TABLE & edge_length_table)
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE min_length, max_length;
  int num_lengths;

  edge_length_table.min_edge_length_freq.Include();
  edge_length_table.min_edge_length_freq.SetAll(0);
  edge_length_table.max_edge_length_freq.Include();
  edge_length_table.max_edge_length_freq.SetAll(0);

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate())
      { continue; }

    compute_min_max_hexahedron_edge_lengths
      (mesh_data, polymesh.VertexList(ipoly),
       polymesh.NumPolyVert(ipoly), vertex_coord,
       min_length, max_length, num_lengths);

    if (num_lengths < 3) { continue; }

    int imin = 
      IJKDATATABLE::compute_bucket
      (min_length, 0, edge_interval, edge_length_table.NumRows());
    int imax = 
      IJKDATATABLE::compute_bucket
      (max_length, 0, edge_interval, edge_length_table.NumRows());

    edge_length_table.min_edge_length_freq.Increment(imin);
    edge_length_table.max_edge_length_freq.Increment(imax);
  }
}


void compute_hex_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const COORD_TYPE min_table_Jacobian,
 const COORD_TYPE table_interval,
 JACOBIAN_TABLE & jacobian_table)
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE min_Jacobian, max_Jacobian;
  int num_rows;
  int num_Jacobian_determinants;

  using namespace IJKDATATABLE;

  jacobian_table.min_jacobian_freq.Include();
  jacobian_table.min_jacobian_freq.SetAll(0);
  jacobian_table.max_jacobian_freq.Include();
  jacobian_table.max_jacobian_freq.SetAll(0);

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate())
      { continue; }

    if (flag_internal_poly) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    compute_min_max_hexahedron_Jacobian_determinants
      (mesh_data, polymesh.VertexList(ipoly),
       polymesh.NumPolyVert(ipoly), vertex_coord,
       min_Jacobian, max_Jacobian, num_Jacobian_determinants);

    if (num_Jacobian_determinants < 1) { continue; }

    int imin = 
      compute_bucket
      (min_Jacobian, min_table_Jacobian, table_interval, 
       jacobian_table.NumRows());
    int imax = 
      compute_bucket
      (max_Jacobian, min_table_Jacobian, table_interval, 
       jacobian_table.NumRows());

    jacobian_table.min_jacobian_freq.Increment(imin);
    jacobian_table.max_jacobian_freq.Increment(imax);
  }
}


void compute_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 const COORD_TYPE min_table_Jacobian,
 const COORD_TYPE table_interval,
 JACOBIAN_TABLE & jacobian_table)
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE min_Jacobian, max_Jacobian;
  int ipoly_with_min, ipoly_with_max;
  int num_Jacobian_determinants;

  jacobian_table.min_jacobian_freq.Include();
  jacobian_table.min_jacobian_freq.SetAll(0);
  jacobian_table.max_jacobian_freq.Include();
  jacobian_table.max_jacobian_freq.SetAll(0);

  for (int iv = 0; iv < vertex_poly_incidence.NumVertices(); iv++) {

    if (flag_internal_vert) {
      if (polymesh.vertex_data[iv].OnBoundary()) { 
        // Vertex iv is not internal.
        continue;
      }
    }

    compute_min_max_hex_vert_Jacobian_determinants
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, iv,
       flag_internal_poly, flag_internal_vert, min_Jacobian, max_Jacobian,
       ipoly_with_min, ipoly_with_max, num_Jacobian_determinants);

    if (num_Jacobian_determinants < 1) { continue; }

    int imin = 
      IJKDATATABLE::compute_bucket
      (min_Jacobian, min_table_Jacobian, table_interval, 
       jacobian_table.NumRows());
    int imax = 
      IJKDATATABLE::compute_bucket
      (max_Jacobian, min_table_Jacobian, table_interval, 
       jacobian_table.NumRows());

    jacobian_table.min_jacobian_freq.Increment(imin);
    jacobian_table.max_jacobian_freq.Increment(imax);
  }
}


void compute_hex_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const COORD_TYPE min_table_Jacobian,
 const COORD_TYPE table_interval,
 JACOBIAN_TABLE & jacobian_table)
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE min_Jacobian, max_Jacobian;
  int num_Jacobian_determinants;

  jacobian_table.min_jacobian_freq.Include();
  jacobian_table.min_jacobian_freq.SetAll(0);
  jacobian_table.max_jacobian_freq.Include();
  jacobian_table.max_jacobian_freq.SetAll(0);

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate())
      { continue; }

    if (flag_internal_poly) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    compute_min_max_hexahedron_normalized_Jacobian_determinants
      (mesh_data, polymesh.VertexList(ipoly),
       polymesh.NumPolyVert(ipoly), vertex_coord,
       min_Jacobian, max_Jacobian, num_Jacobian_determinants);

    if (num_Jacobian_determinants < 1) { continue; }

    int imin = 
      IJKDATATABLE::compute_bucket
      (min_Jacobian, min_table_Jacobian, table_interval, 
       jacobian_table.NumRows());
    int imax = 
      IJKDATATABLE::compute_bucket
      (max_Jacobian, min_table_Jacobian, table_interval, 
       jacobian_table.NumRows());

    jacobian_table.min_jacobian_freq.Increment(imin);
    jacobian_table.max_jacobian_freq.Increment(imax);
  }
}


void compute_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 const COORD_TYPE min_table_Jacobian,
 const COORD_TYPE table_interval,
 JACOBIAN_TABLE & jacobian_table)
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE min_Jacobian, max_Jacobian;
  int ipoly_with_min, ipoly_with_max;
  int num_Jacobian_determinants;

  jacobian_table.min_jacobian_freq.Include();
  jacobian_table.min_jacobian_freq.SetAll(0);
  jacobian_table.max_jacobian_freq.Include();
  jacobian_table.max_jacobian_freq.SetAll(0);

  for (int iv = 0; iv < vertex_poly_incidence.NumVertices(); iv++) {

    if (flag_internal_vert) {
      if (polymesh.vertex_data[iv].OnBoundary()) { 
        // Vertex iv is not internal.
        continue;
      }
    }

    compute_min_max_hex_vert_normalized_Jacobian_determinants
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, iv,
       flag_internal_poly, flag_internal_vert, min_Jacobian, max_Jacobian,
       ipoly_with_min, ipoly_with_max, num_Jacobian_determinants);

    if (num_Jacobian_determinants < 1) { continue; }

    int imin = 
      IJKDATATABLE::compute_bucket
      (min_Jacobian, min_table_Jacobian, table_interval, 
       jacobian_table.NumRows());
    int imax = 
      IJKDATATABLE::compute_bucket
      (max_Jacobian, min_table_Jacobian, table_interval, 
       jacobian_table.NumRows());

    jacobian_table.min_jacobian_freq.Increment(imin);
    jacobian_table.max_jacobian_freq.Increment(imax);
  }
}


void compute_hex_Jacobian_shape
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const COORD_TYPE min_table_value,
 const COORD_TYPE table_interval,
 JACOBIAN_SHAPE_TABLE & jacobian_shape_table)
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE min_jshape, max_jshape;
  int num_jshape_determinants;

  jacobian_shape_table.min_jshape_freq.Include();
  jacobian_shape_table.min_jshape_freq.SetAll(0);
  jacobian_shape_table.max_jshape_freq.Include();
  jacobian_shape_table.max_jshape_freq.SetAll(0);

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (polymesh.poly_data[ipoly].IsDegenerate())
      { continue; }

    if (flag_internal_poly) {
      if (polymesh.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    compute_min_max_hexahedron_Jacobian_shape
      (mesh_data, polymesh.VertexList(ipoly),
       polymesh.NumPolyVert(ipoly), vertex_coord,
       min_jshape, max_jshape, num_jshape_determinants);

    if (num_jshape_determinants < 1) { continue; }

    int imin = 
      IJKDATATABLE::compute_bucket
      (min_jshape, min_table_value, table_interval, 
       jacobian_shape_table.NumRows());
    int imax = 
      IJKDATATABLE::compute_bucket
      (max_jshape, min_table_value, table_interval, 
       jacobian_shape_table.NumRows());

    jacobian_shape_table.min_jshape_freq.Increment(imin);
    jacobian_shape_table.max_jshape_freq.Increment(imax);
  }
}


void compute_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data, const POLYMESH_TYPE & polymesh,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 const COORD_TYPE min_table_value,
 const COORD_TYPE table_interval,
 JACOBIAN_SHAPE_TABLE & jacobian_shape_table)
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE min_jshape, max_jshape;
  int ipoly_with_min, ipoly_with_max;
  int num_jshape_determinants;

  jacobian_shape_table.min_jshape_freq.Include();
  jacobian_shape_table.min_jshape_freq.SetAll(0);
  jacobian_shape_table.max_jshape_freq.Include();
  jacobian_shape_table.max_jshape_freq.SetAll(0);

  for (int iv = 0; iv < vertex_poly_incidence.NumVertices(); iv++) {

    if (flag_internal_vert) {
      if (polymesh.vertex_data[iv].OnBoundary()) { 
        // Vertex iv is not internal.
        continue;
      }
    }

    compute_min_max_hex_vert_Jacobian_shape
      (mesh_data, polymesh, vertex_poly_incidence, vertex_coord, iv,
       flag_internal_poly, flag_internal_vert, min_jshape, max_jshape,
       ipoly_with_min, ipoly_with_max, num_jshape_determinants);

    if (num_jshape_determinants < 1) { continue; }

    int imin = 
      IJKDATATABLE::compute_bucket
      (min_jshape, min_table_value, table_interval, 
       jacobian_shape_table.NumRows());
    int imax = 
      IJKDATATABLE::compute_bucket
      (max_jshape, min_table_value, table_interval, 
       jacobian_shape_table.NumRows());

    jacobian_shape_table.min_jshape_freq.Increment(imin);
    jacobian_shape_table.max_jshape_freq.Increment(imax);
  }
}


// **************************************************
// INTERSECTION ROUTINES
// **************************************************

// Return true if x0 equals x1, or x0 or x1 equal zero.
template <typename TYPEA>
inline bool are_equal_or_zero(const TYPEA x0, const TYPEA x1)
{
  if ((x0 == x1) || (x0 == 0) || (x1 == 0))
    { return(true); }
  else
    { return(false); }
}

// Compute cross product of two 3D vectors.
void compute_cross_product
(const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
 COORD_TYPE v2[DIM3])
{
  determinant_2x2(v0[1], v1[1], v0[2], v1[2], v2[0]);
  determinant_2x2(v0[2], v1[2], v0[0], v1[0], v2[1]);
  determinant_2x2(v0[0], v1[0], v0[1], v1[1], v2[2]);
}

// Compute affine subspace spanning points p0, p1, p2.
// @param span_dimension Dimension of spanning subspace. 0, 1 or 2.
// @param w If span_dimension = 2, w is unit normal to the plane.
// @param w If span_dimension = 1, w is unit line direction.
// @param w If span_dimension = 0, w is p0.
void compute_span
(const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
 const COORD_TYPE p2[DIM3], 
 int & span_dimension, COORD_TYPE w[DIM3])
{
  COORD_TYPE v01[DIM3], v02[DIM3], v12[DIM3];
  COORD_TYPE normalized_v01[DIM3], normalized_v02[DIM3], 
    normalized_v12[DIM3];
  COORD_TYPE v01_mag, v02_mag, v12_mag;
  COORD_TYPE_PTR u0, u1, u2;
  COORD_TYPE u1_orth[DIM3], u2_orth[DIM3];
  COORD_TYPE normalized_u1_orth[DIM3], normalized_u2_orth[DIM3];
  COORD_TYPE u1_orth_mag, u2_orth_mag, w_mag;

  subtract_coord_3D(p1, p0, v01);
  subtract_coord_3D(p2, p0, v02);
  subtract_coord_3D(p2, p1, v12);

  normalize_vector_robust(DIM3, v01, 0, normalized_v01, v01_mag);
  normalize_vector_robust(DIM3, v02, 0, normalized_v02, v02_mag);
  normalize_vector_robust(DIM3, v12, 0, normalized_v12, v12_mag);

  if (v01_mag <= 0 && v02_mag <= 0 && v12_mag <= 0) {
    span_dimension = 0;
    copy_coord_3D(p0, w);
    return;
  }

  if (v01_mag >= v02_mag && v01_mag >= v12_mag) {
    u0 = normalized_v01;
    u1 = v02;
    u2 = v12;
  }
  else if (v02_mag >= v12_mag) {
    u0 = normalized_v02;
    u1 = v01;
    u2 = v12;
  }
  else {
    u0 = normalized_v12;
    u1 = v01;
    u2 = v02;
  }

  compute_orthogonal_vector(DIM3, u1, u0, u1_orth);
  compute_orthogonal_vector(DIM3, u2, u0, u2_orth);

  normalize_vector_robust(DIM3, u1_orth, 0, normalized_u1_orth, u1_orth_mag);
  normalize_vector_robust(DIM3, u2_orth, 0, normalized_u2_orth, u2_orth_mag);

  if (u1_orth_mag == 0 && u2_orth_mag == 0) {
    span_dimension = 1;
    copy_coord_3D(u0, w);
    return;
  }

  span_dimension = 2;
  if (u1_orth_mag >= u2_orth_mag) {
    compute_cross_product(u0, normalized_u1_orth, w);
  }
  else {
    compute_cross_product(u0, normalized_u2_orth, w);
  }

  // Renormalized, just in case.
  normalize_vector_robust(DIM3, w, 0, w, w_mag);
}

template <typename DTYPE, typename CTYPE>
CTYPE compute_max_coord(const DTYPE dimension, const CTYPE coord[])
{
  CTYPE max_coord;

  if (dimension <= 0) { return(0); }

  max_coord = coord[0];
  for (int d = 1; d < dimension; d++) {
    if (max_coord < coord[d]) 
      { max_coord = coord[d]; }
  }
}


inline void compute_projection_coefficient
(const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3], 
 const COORD_TYPE u[DIM3], COORD_TYPE & s)
{
  COORD_TYPE v01[DIM3];

  subtract_coord_3D(p1, p0, v01);
  compute_inner_product_3D(v01, u, s);
}

// Compute min/max coordinates of triangle vertices.
void compute_min_max_coord
(const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
 const COORD_TYPE p2[DIM3],
 COORD_TYPE min_coord[DIM3], COORD_TYPE max_coord[DIM3])
{
  COORD_TYPE minC, maxC;

  for (int d = 0; d < DIM3; d++) {
    minC = p0[d];
    maxC = p0[d];
    if (p1[d] < minC) { minC = p1[d]; }
    if (p1[d] > maxC) { maxC = p1[d]; }
    if (p2[d] < minC) { minC = p2[d]; }
    if (p2[d] > maxC) { maxC = p2[d]; }
    min_coord[d] = minC;
    max_coord[d] = maxC;
  }
}

// Compute orientation of (p1-p0), (p2-p0) and w
int compute_orientation
(const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
 const COORD_TYPE p2[DIM3], const COORD_TYPE w[DIM3],
 const COORD_TYPE epsilon)
{
  COORD_TYPE v01[DIM3], v02[DIM3];
  COORD_TYPE det;

  subtract_coord_3D(p1, p0, v01);
  subtract_coord_3D(p2, p0, v02);

  determinant_3x3(v01, v02, w, det);
  if (abs(det) < epsilon) { det = 0; }

  if (det < 0) { return(-1); }
  else if (det > 0) { return(1); }
  else { return(0); }
}


/// Return true if triangle contains point.
/// Assumes plane containes point and triangle.
/// @param normal[] Normal to plane containing point and triangle.
bool does_triangle_contain_point
(const COORD_TYPE w0[DIM3], const COORD_TYPE w1[DIM3],
 const COORD_TYPE w2[DIM3], const COORD_TYPE p[DIM3],
 const COORD_TYPE normal[DIM3], const COORD_TYPE epsilon)
{
  COORD_TYPE min_coord[DIM3], max_coord[DIM3];
  int sign0, sign1, sign2;

  // Handles cases where w0, w1, w2 and p are collinear.
  compute_min_max_coord(w0, w1, w2, min_coord, max_coord);
  for (int d = 0; d < DIM3; d++) {
    if (p[d] + epsilon < min_coord[d]) { return(false); }
    if (p[d] - epsilon > max_coord[d]) { return(false); }
  }

  // Check that intersection_point is in triangle (w0,w1,w2).
  sign0 = compute_orientation(w0, w1, p, normal, epsilon);
  sign1 = compute_orientation(w1, w2, p, normal, epsilon);
  sign2 = compute_orientation(w2, w0, p, normal, epsilon);

  if (!are_equal_or_zero(sign0, sign1)) { return(false); }
  if (!are_equal_or_zero(sign1, sign2)) { return(false); }
  if (!are_equal_or_zero(sign0, sign2)) { return(false); }

  return(true);
}

// Compute intersection of edge (v0,v1) and triangle (w0,w1,w2).
// Use vectors (w1-w0) and (w2-w0)
// Return true if edge intersects triangle.
bool intersect_edge_triangle
(const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
 const COORD_TYPE w0[DIM3], const COORD_TYPE w1[DIM3],
 const COORD_TYPE w2[DIM3], const COORD_TYPE epsilon,
 COORD_TYPE intersection_point[DIM3])
{
  COORD_TYPE w[DIM3];
  COORD_TYPE v01[DIM3], normalized_v01[DIM3];
  COORD_TYPE v0w0[DIM3];
  COORD_TYPE u0[DIM3], u1[DIM3];
  COORD_TYPE w_mag, v01_mag;
  COORD_TYPE s0, s1, t;
  COORD_TYPE product0, product1;
  int span_dimension;

  // No intersection if line segment and triangle share a vertex.
  if (is_coord_equal_3D(v0, w0)) { return(false); };
  if (is_coord_equal_3D(v0, w1)) { return(false); };
  if (is_coord_equal_3D(v0, w2)) { return(false); };
  if (is_coord_equal_3D(v1, w0)) { return(false); };
  if (is_coord_equal_3D(v1, w1)) { return(false); };
  if (is_coord_equal_3D(v1, w2)) { return(false); };

  compute_span(w0, w1, w2, span_dimension, w);

  if (span_dimension < 2) { return(false); };

  // Check w is non-zero.
  compute_magnitude_3D(w, w_mag);
  if (w_mag <= epsilon) { return(false); };

  subtract_coord_3D(v1, v0, v01);
  subtract_coord_3D(w0, v0, v0w0);

  compute_inner_product_3D(v01, w, s0);
  compute_inner_product_3D(v0w0, w, s1);

  if (s0 == 0) { return(false); }

  if (fabs(s1) <= (1+epsilon)*fabs(s0)) 
    { t = s1/s0; }
  else
    { return(false); }

  if (t < -epsilon || t > 1+epsilon) 
    { return(false); }

  if (t*w_mag < -epsilon || t*w_mag > 1+epsilon)
    { return(false); }

  add_scaled_coord_3D(t, v01, v0, intersection_point);

  // Check that intersection_point is on line segment (v0,v1).
  subtract_coord_3D(intersection_point, v0, u0);
  subtract_coord_3D(intersection_point, v1, u1);
  compute_inner_product_3D(u0, v01, product0);
  compute_inner_product_3D(u1, v01, product1);

  if (product0 < epsilon || product1 > -epsilon) { return(false); }

  bool flag = does_triangle_contain_point
    (w0, w1, w2, intersection_point, w, epsilon);

  return(flag);
}

// **************************************************
// Class MESH_INFO
// **************************************************

void MESH_INFO::Init()
{
  num_poly_with_duplicate_vertices = 0;
  num_duplicate_poly = 0;
  num_non_manifold_facets = 0;
  num_non_manifold_edges = 0;
  num_non_manifold_vertices = 0;
  num_deep_non_manifold_vertices = 0;
  num_poly_with_orientation_conflicts = 0;
  num_hex_facet_pairs_sharing_exactly_two_edges = 0;
}

// Return true if all non-manifold numbers are zero.
bool MESH_INFO::AreAllNonManifoldZero() const
{
  if (num_poly_with_duplicate_vertices > 0) { return(false); }
  if (num_duplicate_poly > 0) { return(false); }
  if (num_hex_facet_pairs_sharing_exactly_two_edges > 0) { return(false); }
  if (num_non_manifold_facets > 0) { return(false); }
  if (num_non_manifold_edges > 0) { return(false); }
  if (num_non_manifold_vertices > 0) { return(false); }
  if (num_deep_non_manifold_vertices > 0) { return(false); }

  return(true);
}

// Return true if all numbers are zero.
bool MESH_INFO::AreAllZero() const
{
  if (!AreAllNonManifoldZero()) { return(false); }
  if (num_poly_with_orientation_conflicts > 0) { return(false); }

  return(true);
}


// **************************************************
// Class MESH_DATA member functions
// **************************************************

void MESH_DATA::Init()
{
  dimension = DIM3;
  mesh_dimension = DIM2;
  num_vertices = 0;

  are_boundary_facets_identified = false;
  are_non_manifold_facets_identified = false;
  are_hex_sharing_exactly_two_facet_edges_identified = false;

  // Default to positive orientation.
  orientation = 1;
}


// **************************************************
// Class GRID_OF_BINS_3D
// **************************************************

void GRID_OF_BINS_3D::Init(const int k)
{
  int k2 = k;
  if (k2 < 1) { k2 = 1; }

  for (int d = 0; d < DIM3; d++)
    { num_bins_along_axis[d] = k2; }

  num_bins = 1;
  for (int d = 0; d < DIM3; d++)
    { num_bins = num_bins*num_bins_along_axis[d]; }

  typedef vector<int> * BIN_PTR;

  bin = new BIN_PTR[num_bins];
  for (int ibin = 0; ibin < num_bins; ibin++) 
    { bin[ibin] = NULL; }

  // initialize minC[] and maxC[]
  for (int d = 0; d < DIM3; d++) {
    minC[d] = 0;
    maxC[d] = 1;
  }

  bins_are_empty = true;
}

// Destructor
GRID_OF_BINS_3D::~GRID_OF_BINS_3D()
{
  for (int ibin = 0; ibin < num_bins; ibin++) {
    if (bin[ibin] != NULL) { delete bin[ibin]; }
    bin[ibin] = NULL;
  }

  delete [] bin;
  bin = NULL;

  num_bins = 0;
}

void GRID_OF_BINS_3D::SetMinCoord(const COORD_TYPE minC[DIM3])
{
  IJK::PROCEDURE_ERROR error("GRID_OF_BINS_3D::SetMinCoord");

  if (!bins_are_empty) {
    error.AddMessage("Programming error. Illegal to call SetMinCoord after inserting any elements.");
    throw error;
  }

  for (int d = 0; d < DIM3; d++)
    { this->minC[d] = minC[d]; }
}

void GRID_OF_BINS_3D::SetMaxCoord(const COORD_TYPE maxC[DIM3])
{
  IJK::PROCEDURE_ERROR error("GRID_OF_BINS_3D::SetMaxCoord");

  if (!bins_are_empty) {
    error.AddMessage("Programming error. Illegal to call SetMaxCoord after inserting any elements.");
    throw error;
  }

  for (int d = 0; d < DIM3; d++)
    { this->maxC[d] = maxC[d]; }
}

int GRID_OF_BINS_3D::LocateBinCoord(const int d, const COORD_TYPE c) const
{
  const int nbin = num_bins_along_axis[d];
  const COORD_TYPE bin_width = (maxC[d] - minC[d])/nbin;

  if (c <= minC[d] || bin_width <= 0) { return(0); }

  if (nbin*bin_width <= (c-minC[d])) { return(nbin-1); }

  int binCoord = int((c-minC[d])/bin_width);

  if (binCoord >= nbin) { binCoord = nbin-1; }
  
  return(binCoord);
}

int GRID_OF_BINS_3D::ComputeBinIndex(const int binCoord[DIM3]) const
{
  int index = binCoord[0] + 
    num_bins_along_axis[0]*(binCoord[1]+ binCoord[2]*num_bins_along_axis[1]);
  return(index);
}

void GRID_OF_BINS_3D::ComputeMinMaxBinCoord
(const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
 const COORD_TYPE v2[DIM3],
 int min_grid_coord[DIM3], int max_grid_coord[DIM3]) const
{
  BOUNDING_BOX bounding_box(DIM3);

  bounding_box.SetCoord(v0, v0);
  bounding_box.Extend(v1);
  bounding_box.Extend(v2);

  for (int d = 0; d < DIM3; d++) {
    min_grid_coord[d] = LocateBinCoord(d, bounding_box.MinCoord(d));
    max_grid_coord[d] = LocateBinCoord(d, bounding_box.MaxCoord(d));
  }
}

void GRID_OF_BINS_3D::InsertTri
(const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
 const COORD_TYPE v2[DIM3], const int triangle_index)
{
  int min_grid_coord[DIM3];
  int max_grid_coord[DIM3];

  ComputeMinMaxBinCoord(v0, v1, v2, min_grid_coord, max_grid_coord);

  int binCoord[DIM3];
  for (binCoord[2] = min_grid_coord[2]; 
       binCoord[2] <= max_grid_coord[2]; binCoord[2]++) {
    for (binCoord[1] = min_grid_coord[1]; 
         binCoord[1] <= max_grid_coord[1]; binCoord[1]++) {
      binCoord[0] = min_grid_coord[0];
      int ibin = ComputeBinIndex(binCoord);

      for (int x = min_grid_coord[0]; x <= max_grid_coord[0]; x++) {
        if (bin[ibin] == NULL)
          { bin[ibin] = new std::vector<int>; }
        bin[ibin]->push_back(triangle_index);
        ibin++;
      }
    }
  }

  bins_are_empty = false;
}


// **************************************************
// Class POLY_DATA member functions
// **************************************************

void IJKMESHINFO::POLY_DATA::Init()
{
  is_degenerate = false;
  is_duplicate = false;
  contains_non_manifold_facet = false;
  contains_non_manifold_edge = false;
  contains_boundary_facet = false;
  orientation_conflict = false;
}


// **************************************************
// Class VERTEX_DATA member functions
// **************************************************

void IJKMESHINFO::VERTEX_DATA::Init()
{
  on_boundary = false;
}


// **************************************************
// Class ANGLE_TABLE member functions
// **************************************************

void ANGLE_TABLE::HideAllExceptAngleColumn()
{
  angle.Show();
  min_polygon_angle_freq.Hide();
  max_polygon_angle_freq.Hide();
}

void ANGLE_TABLE::WriteColumnLabels
(std::ostream & out, const std::string & separator) const
{
  angle.WriteLabel(out, separator);
  min_polygon_angle_freq.WriteLabel(out, separator);
  max_polygon_angle_freq.WriteLabel(out, separator);
}

void ANGLE_TABLE::WriteColumnData
(std::ostream & out, const std::string & separator, 
 const NUM_TYPE width) const
{
  for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
    angle.WriteData(out, "", width, irow);
    min_polygon_angle_freq.WriteData(out, "  ", width, irow);
    max_polygon_angle_freq.WriteData(out, "  ", width, irow);
    out << endl;
  }
}

void ANGLE_TABLE::WriteNormalizedColumnData
(std::ostream & out, const std::string & separator, const NUM_TYPE width,
 const double normalization_factor) const
{
  for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
    angle.WriteData(out, "", width, irow);
    min_polygon_angle_freq.WriteNormalizedData
      (out, "  ", width, normalization_factor, irow);
    max_polygon_angle_freq.WriteNormalizedData
      (out, "  ", width, normalization_factor, irow);
    out << endl;
  }
}


// **************************************************
// Class EDGE_LENGTH_TABLE member functions
// **************************************************

void EDGE_LENGTH_TABLE::HideAllExceptEdgeLengthColumn()
{
  edge_length.Show();
  min_edge_length_freq.Hide();
  max_edge_length_freq.Hide();
}

void EDGE_LENGTH_TABLE::WriteColumnLabels
(std::ostream & out, const std::string & separator) const
{
  edge_length.WriteLabel(out, separator);
  min_edge_length_freq.WriteLabel(out, separator);
  max_edge_length_freq.WriteLabel(out, separator);
}

void EDGE_LENGTH_TABLE::WriteColumnData
(std::ostream & out, const std::string & separator, 
 const NUM_TYPE width) const
{
  for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
    edge_length.WriteData(out, "", width, irow);
    min_edge_length_freq.WriteData(out, "  ", width, irow);
    max_edge_length_freq.WriteData(out, "  ", width, irow);
    out << endl;
  }
}

void EDGE_LENGTH_TABLE::WriteNormalizedColumnData
(std::ostream & out, const std::string & separator, const NUM_TYPE width,
 const double normalization_factor) const
{
  for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
    edge_length.WriteData(out, "", width, irow);
    min_edge_length_freq.WriteNormalizedData
      (out, "  ", width, normalization_factor, irow);
    max_edge_length_freq.WriteNormalizedData
      (out, "  ", width, normalization_factor, irow);
    out << endl;
  }
}


// **************************************************
// Class JACOBIAN_TABLE member functions
// **************************************************

void JACOBIAN_TABLE::HideAllExceptJacobianColumn()
{
  jacobian.Show();
  min_jacobian_freq.Hide();
  max_jacobian_freq.Hide();
}

void JACOBIAN_TABLE::WriteColumnLabels
(std::ostream & out, const std::string & separator) const
{
  jacobian.WriteLabel(out, separator);
  min_jacobian_freq.WriteLabel(out, separator);
  max_jacobian_freq.WriteLabel(out, separator);
}

void JACOBIAN_TABLE::WriteColumnData
(std::ostream & out, const std::string & separator, 
 const NUM_TYPE width) const
{
  for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
    jacobian.WriteData(out, "", width, irow);
    min_jacobian_freq.WriteData(out, "  ", width, irow);
    max_jacobian_freq.WriteData(out, "  ", width, irow);
    out << endl;
  }
}

void JACOBIAN_TABLE::WriteNormalizedColumnData
(std::ostream & out, const std::string & separator, const NUM_TYPE width,
 const double normalization_factor) const
{
  for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
    jacobian.WriteData(out, "", width, irow);
    min_jacobian_freq.WriteNormalizedData
      (out, "  ", width, normalization_factor, irow);
    max_jacobian_freq.WriteNormalizedData
      (out, "  ", width, normalization_factor, irow);
    out << endl;
  }
}


// **************************************************
// Class JACOBIAN_SHAPE_TABLE member functions
// **************************************************

void JACOBIAN_SHAPE_TABLE::HideAllExceptShapeColumn()
{
  jshape.Show();
  min_jshape_freq.Hide();
  max_jshape_freq.Hide();
}

void JACOBIAN_SHAPE_TABLE::WriteColumnLabels
(std::ostream & out, const std::string & separator) const
{
  jshape.WriteLabel(out, separator);
  min_jshape_freq.WriteLabel(out, separator);
  max_jshape_freq.WriteLabel(out, separator);
}

void JACOBIAN_SHAPE_TABLE::WriteColumnData
(std::ostream & out, const std::string & separator, 
 const NUM_TYPE width) const
{
  for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
    jshape.WriteData(out, "", width, irow);
    min_jshape_freq.WriteData(out, "  ", width, irow);
    max_jshape_freq.WriteData(out, "  ", width, irow);
    out << endl;
  }
}

void JACOBIAN_SHAPE_TABLE::WriteNormalizedColumnData
(std::ostream & out, const std::string & separator, const NUM_TYPE width,
 const double normalization_factor) const
{
  for (NUM_TYPE irow = 0; irow < this->NumRows(); irow++) {
    jshape.WriteData(out, "", width, irow);
    min_jshape_freq.WriteNormalizedData
      (out, "  ", width, normalization_factor, irow);
    max_jshape_freq.WriteNormalizedData
      (out, "  ", width, normalization_factor, irow);
    out << endl;
  }
}


// **************************************************
// MISCELLANEOUS ROUTINES
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}


PARAMETER get_parameter_token(char * s)
// convert string s into parameter token
{
  for (int i = 0; i < int(UNKNOWN_PARAM); i++)
    if (string(s) == string(parameter_string[i]))
      return(PARAMETER(i));
  return(UNKNOWN_PARAM);
}

void get_coord(const char * s, vector<COORD_TYPE> & coord)
{
  istringstream coord_string;

  coord.clear();

  string s2 = s;
  // remove trailing blanks from s2
  size_t pos = 0;
  for (size_t i = 0; i < s2.length(); i++) {
    if (!isspace(s2[i])) { pos = i+1; }
  }
  if (pos < s2.length()) { s2.erase(pos); };

  coord_string.str(s2);
  while (coord_string.good()) {
    COORD_TYPE c;
    coord_string >> c;
    coord.push_back(c);
  }

  if (coord_string.fail() && !coord_string.eof()) {
    cerr << "Error reading coordinates: "
         << "\"" << s << "\"" << endl;
    cerr << "  Non-numeric character in coordinate string." << endl;
    exit(600);
  }

}

void parse_command_line(int argc, char **argv)
{
  float x;
  IJK::ERROR error;

  if (argc == 1) { usage_error(); };

  int iarg = 1;

  try {
    while (iarg < argc && argv[iarg][0] == '-') {
      PARAMETER param = get_parameter_token(argv[iarg]);
      if (param == UNKNOWN_PARAM) break;

      switch(param) {

      case POLYFILE_PARAM:
        flag_polyfile = true;
        cerr << "WARNING: Option -polyfile is deprecated." << endl;
        break;

      case MESH_DIM_PARAM:
        mesh_data.mesh_dimension = get_arg_int(iarg, argc, argv, error);
        is_mesh_dimension_set = true;
        iarg++;
        break;

      case REVERSE_ORIENT_PARAM:
        mesh_data.orientation = -1;
        break;

      case VERTEX_PARAM:
        vertex_index = get_arg_int(iarg, argc, argv, error);
        iarg++;
        vertex_info_flag = true;
        io_info.flag_general_info = false;
        break;

      case SIMPLEX_PARAM:
        simplex_index = get_arg_int(iarg, argc, argv, error);
        iarg++;
        simplex_info_flag = true;
        io_info.flag_general_info = false;
        break;

      case POLY_PARAM:
        poly_index = get_arg_int(iarg, argc, argv, error);
        iarg++;
        poly_info_flag = true;
        io_info.flag_general_info = false;
        break;

      case VLIST_PARAM:
        vlist_flag = true;
        io_info.flag_general_info = false;
        break;

      case VLIST_MIN_PARAM:
        vlist_min_flag = true;
        io_info.flag_output_min_Jacobian_determinant = true;
        io_info.flag_output_min_normalized_Jacobian_determinant = true;
        io_info.flag_output_min_Jacobian_shape = true;
        io_info.flag_general_info = false;
        break;

      case PLIST_PARAM:
        plist_flag = true;
        io_info.flag_general_info = false;
        break;

      case CONTAINSV_PARAM:
        contains_vertex_index = get_arg_int(iarg, argc, argv, error);
        iarg++;
        contains_vertex_flag = true;
        io_info.flag_general_info = false;
        break;

      case CONTAINSE_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%d", &edge_end0_index);
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%d", &edge_end1_index);
        contains_edge_flag = true;
        io_info.flag_general_info = false;
        break;

      case MANIFOLD_PARAM:
        manifold_flag = true;
        io_info.flag_general_info = false;
        break;

      case ORIENTED_MANIFOLD_PARAM:
        manifold_flag = true;
        oriented_manifold_flag = true;
        io_info.flag_general_info = false;
        break;

      case CHECK_FACET_INTERSECTIONS_PARAM:
        check_facet_intersections_flag = true;
        io_info.flag_general_info = false;
        break;

      case SELFI_PARAM:
        flag_report_self_intersections = true;
        flag_use_grid_of_bins = true;
        io_info.flag_general_info = false;
        break;

      case SELFI_NO_GRID_PARAM:
        flag_report_self_intersections = true;
        flag_use_grid_of_bins = false;
        io_info.flag_general_info = false;
        break;

      case GRID_LENGTH_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%d", &num_bins_per_axis);
        is_num_bins_per_axis_set = true;
        break;

      case MINC_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        get_coord(argv[iarg], min_coord);
        is_min_coord_set = true;
        break;

      case MAXC_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        get_coord(argv[iarg], max_coord);
        is_max_coord_set = true;
        break;

      case MIN_NUMV_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%d", &min_num_polyv_output);
        is_min_num_polyv_output_set = true;
        break;

      case MAX_NUMV_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%d", &max_num_polyv_output);
        is_max_num_polyv_output_set = true;
        break;

      case ANGLE_LE_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%f", &x);
        io_info.angle_le.Set(x);
        break;

      case ANGLE_GE_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%f", &x);
        io_info.angle_ge.Set(x);
        break;

      case FACET_ANGLE_LE_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%f", &x);
        io_info.facet_angle_le.Set(x);
        break;

      case FACET_ANGLE_GE_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%f", &x);
        io_info.facet_angle_ge.Set(x);
        break;

      case PJACOBIAN_PARAM:
        io_info.flag_pJacobian.Set(true);
        break;

      case VJACOBIAN_PARAM:
        io_info.flag_vJacobian.Set(true);
        break;

      case LIST_DUP_PARAM:
        flag_list_duplicate_vertices = true;
        flag_list_duplicate_poly = true;
        break;

      case INTERNAL_POLY_PARAM:
        flag_internal_poly = true;
        break;

      case INTERNAL_VERT_PARAM:
        flag_internal_vert = true;
        break;

      case REPORT_DEEP_PARAM:
        flag_report_deep = true;
        break;

      case OUT_VALUES_PARAM:
        flag_output_only_values = true;
        break;

      case OUT_MIN_ANGLE_PARAM:
        io_info.flag_output_min_angle = true;
        io_info.flag_general_info = false;
        break;

      case OUT_MAX_ANGLE_PARAM:
        io_info.flag_output_max_angle = true;
        io_info.flag_general_info = false;
        break;

      case OUT_MIN_JACOBIAN_DET_PARAM:
        io_info.flag_output_min_Jacobian_determinant = true;
        io_info.flag_output_min_normalized_Jacobian_determinant = true;
        io_info.flag_output_min_Jacobian_shape = true;
        io_info.flag_general_info = false;
        break;

      case OUT_MAX_JACOBIAN_DET_PARAM:
        io_info.flag_output_max_Jacobian_determinant = true;
        io_info.flag_general_info = false;
        break;

      case OUT_MIN_NORMALIZED_JACOBIAN_DET_PARAM:
        io_info.flag_output_min_normalized_Jacobian_determinant = true;
        io_info.flag_general_info = false;
        break;

      case OUT_MAX_NORMALIZED_JACOBIAN_DET_PARAM:
        io_info.flag_output_max_normalized_Jacobian_determinant = true;
        io_info.flag_general_info = false;
        break;

      case PLOT_ANGLES_PARAM:
        flag_plot_angles = true;
        flag_plot_min_polygon_angles = true;
        flag_plot_max_polygon_angles = true;
        break;

      case PLOT_EDGE_LENGTHS_PARAM:
        flag_plot_min_edge_lengths = true;
        flag_plot_max_edge_lengths = true;
        break;

      case PLOT_JACOBIAN_PARAM:
        flag_plot_min_jacobian = true;
        flag_plot_max_jacobian = true;
        break;

      case PLOT_JACOBIAN_SHAPE_PARAM:
        flag_plot_min_jacobian_shape = true;
        flag_plot_max_jacobian_shape = true;
        break;

      case FOR_EACH_TYPE_PARAM:
        flag_for_each_type = true;
        break;

      case MAX_OUT_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%d", &io_info.max_num_poly_out);
        break;

      case TERSE_PARAM:
        terse_flag = true;
        break;

      case HELP_PARAM:
        help_msg();
        break;

      default:
        usage_error();
      };

      iarg++;
    };
  }
  catch (ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    exit(10);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  }
  
  if (iarg >= argc) {
    cerr << "Error.  Missing input file name." << endl;
    usage_error();
  };

  input_filename = argv[iarg]; 
  iarg++;

  if (iarg < argc) { 
    output_filename = argv[iarg]; 
    iarg++;
  }

  if (iarg != argc) { usage_error(); }

  if (!io_info.flag_pJacobian.IsSet() && !io_info.flag_vJacobian.IsSet()) {

    if (vlist_flag || vlist_min_flag) { io_info.flag_vJacobian.Set(true); }
    else { io_info.flag_pJacobian.Set(true); }
  }

  if (flag_polyfile) {

    if (contains_edge_flag) {
      cerr << "Option -containse not implemented with -polyfile." << endl;
      exit(20);
    }

    if (manifold_flag) {
      cerr << "Option -manifold not implemented with -polyfile." << endl;
      exit(20);
    }
  }

  if (io_info.angle_le.IsSet()) {
    if (io_info.angle_le.Value() < 0 || io_info.angle_le.Value() > 180) {
      cerr << "Illegal value " << io_info.angle_le.Value() 
           << " for option -angle_le <A>." << endl;
      cerr << "  Angle <A> must be in range [0,180].";
      exit(20);
    }
  }

  if (io_info.angle_ge.IsSet()) {
    if (io_info.angle_ge.Value() < 0 || io_info.angle_ge.Value() > 180) {
      cerr << "Illegal value " << io_info.angle_ge.Value() 
           << " for option -angle_ge <A>." << endl;
      cerr << "  Angle <A> must be in range [0,180].";
      exit(20);
    }
  }

  if (!io_info.flag_output_min_angle && !io_info.flag_output_max_angle &&
      !io_info.angle_le.IsSet() && !io_info.angle_ge.IsSet() &&
      !io_info.facet_angle_le.IsSet() && !io_info.facet_angle_ge.IsSet() &&
      !io_info.flag_output_min_edge_length && 
      !io_info.flag_output_max_edge_length) 
    { io_info.flag_output_all_min_max = true; }
}

void check_input()
{
  if (!flag_simplex_file) {

    if (contains_edge_flag) {
      cerr << "Usage error.  Option \"-containse\" only implemented for mesh of simplices." << endl;
      exit(20);
    }
  }

}

