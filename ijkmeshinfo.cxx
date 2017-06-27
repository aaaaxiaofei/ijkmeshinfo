/// \file ijkmeshinfo.cxx
/// compute mesh information
/// Version 0.2.1

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2008-2017 Rephael Wenger

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
#include "ijkcube.txx"
#include "ijkIO.txx"
#include "ijkmerge.txx"
#include "ijkmesh.txx"
#include "ijkmesh_datastruct.txx"
#include "ijkprint.txx"
#include "ijkstring.txx"

#include "ijkmeshinfo.h"

using namespace std;
using namespace IJKMESHINFO;
using namespace IJK;

// types
typedef int VERTEX_INDEX;
typedef typename IJK::POLYMESH<VERTEX_INDEX,int> POLYMESH_TYPE;
typedef typename IJK::VERTEX_POLY_INCIDENCE<int,int> 
VERTEX_POLY_INCIDENCE_TYPE;
typedef typename IJK::CUBE_FACE_INFO<int,int,int> CUBE_TYPE;

// global variables
int dimension(DIM3);
int mesh_dimension;
bool is_mesh_dimension_set(false);
int num_vertices = 0;
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
ANGLE_TYPE angle_le = 0;
ANGLE_TYPE angle_ge = 180;
bool is_min_angle_set = false;
bool is_max_angle_set = false;
bool flag_simplex_file(false);
bool flag_cube_file(false);
bool flag_polyfile(false);
std::vector<COORD_TYPE> min_coord;
std::vector<COORD_TYPE> max_coord;
bool flag_list_duplicate_vertices = false;
bool flag_list_duplicate_poly = false;
bool flag_report_deep = false;         // Report only facets "deep" in bounding box.
bool flag_output_only_values = false;  // Output only values without any text.
bool flag_output_min_angle = false;    // Output minimum angle.
bool flag_output_max_angle = false;    // Output maximum angle.
bool flag_internal = false;            // Output angles of internal polygons.
bool flag_for_each_type;               // Report min/max for tri, quad, etc.
bool flag_report_self_intersections = false;  // Report self intersections.
double selfI_epsilon = 1.0e-10;
bool flag_use_grid_of_bins = true;
int num_bins_per_axis = 10;
bool is_num_bins_per_axis_set = false;
bool flag_plot_angles = false;               // Plot some angles.
bool flag_plot_min_polygon_angles = false;   // Plot min polygon angles.
bool flag_plot_max_polygon_angles = false;   // Plot max polygon angles.
bool flag_normalize = false;
bool flag_silent_write = false; // if true, suppress message "Writing table..."
int DEFAULT_TABLE_COLUMN_WIDTH = 8;
int DEFAULT_TABLE_PRECISION = 4;
MESH_INFO mesh_info;

vector<bool> is_degenerate;        // true if simplex/poly is degenerate
vector<bool> is_duplicate;         // true if simplex/poly is duplicate

// true if simplex contains a nonmanifold facet
vector<bool> contains_nonmanifold_facet; 

// true if simplex contains a boundary facet
vector<bool> contains_boundary_facet; 

vector<int> nonmanifold_facet_vert; // list of nonmanifold facet vertices
vector<int> nonmanifold_edge_vert;  // list of nonmanifold edge vertices
vector<int> boundary_facet_vert;    // list of boundary facet vertices

// true if boundary facet is inside bounding box
vector<bool> internal_boundary_facet;  

// true if boundary facet is far from bounding_box
vector<bool> far_from_bounding_box;

int num_internal_boundary_facets = 0;
int num_deep_boundary_facets = 0;
vector<bool> in_nonmanifold_facet;  // true if vertex is in nonmanifold facet
vector<bool> in_nonmanifold_edge;   // true if vertex is in nonmanifold edge
vector<bool> nonmanifold_vert;      // true if vertex is nonmanifold
vector<int> nonmanifold_vert_list;  // list of nonmanifold vertices

vector<int> sorted_poly;            // list of polytopes in sorted order

int vertex_index = 0;
int contains_vertex_index = 0;
int simplex_index = 0;
int poly_index = 0;
int edge_end0_index = 0;
int edge_end1_index = 0;

bool general_info_flag = true;
bool vertex_info_flag = false;
bool simplex_info_flag = false;
bool poly_info_flag = false;
bool manifold_flag = false;
bool terse_flag = false;
bool vlist_flag = false;
bool plist_flag = false;
bool contains_vertex_flag = false;
bool contains_edge_flag = false;

// compute info routines
void compute_bounding_box();
int compute_num_edges();
int count_deep_vertices(const std::vector<int> & vlist);
void compute_facet_info();
int count_num_poly(const int num_poly_vert);

// compute angle routines
void compute_min_max_polygon_angles
(const bool flag_internal, ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle);
void compute_min_max_polygon_angles
(const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle);
void compute_min_max_polygon_angles
(const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle);
void compute_min_max_polygon_angles
(const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle);
void compute_min_max_polygon_angles
(const int ipoly, ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & num_angle);
void compute_num_polygon_angles
(const bool flag_internal, 
 const ANGLE_TYPE & min_angle, const ANGLE_TYPE & max_angle,
 int & num_le, int & num_ge);
void compute_min_max_polygon_cos
(const int ipoly, COORD_TYPE & min_cos, COORD_TYPE & max_cos,
 int & num_angle);

// compute edge length routines
void compute_min_max_edge_lengths
(const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_edge_length, ANGLE_TYPE & max_edge_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length);
void compute_min_max_edge_lengths
(const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_edge_length, ANGLE_TYPE & max_edge_length);

// mesh processing routines
void sort_poly();
void identify_duplicates();
int identify_poly_with_duplicate_vertices();
int identify_duplicate_polytopes();
void identify_polygon_edges();
void set_in_nonmanifold_facet();
void set_in_nonmanifold_edge();

// manifold routines
void identify_nonmanifold();
void identify_nonmanifold_and_boundary_facets();
int identify_nonmanifold_vertices();
int identify_nonmanifold_edges();
bool is_internal
(const BOUNDING_BOX & bounding_box, const vector<int> & facet_vlist, 
 const int numv_per_facet, const int jf);

// plot routines
void compute_polygon_angles(ANGLE_TABLE & angle_table);
template <class TABLE_TYPE>
void write_table_gplt
(const string & filename_prefix, const string & filename_suffix,
 const TABLE_TYPE & table);
template <class TABLE_TYPE>
void write_table_gplt(ofstream & ofile, const TABLE_TYPE & table);


// output info routines
void output_general_info();
void output_vertex_info(const int vertex_index);
void output_simplex_info(const int simplex_index);
void output_poly_info(const int poly_index);
void output_manifold_info();
void output_degenerate_poly();
void output_duplicate_poly();
void output_duplicate_vertices();
void output_nonmanifold_facets();
void output_nonmanifold_edges();
void output_nonmanifold_vertices();
void output_internal_boundary_facets();
void output_poly(const int ipoly);
void output_simplex(const int * simplex);
void output_vertex_list();
void output_simplices();
void output_polytopes();
bool output_self_intersections();
bool output_self_intersections_using_grid_of_bins();
void output_manifold_and_boundary_counts();
void output_small_angles(const bool flag_internal, const ANGLE_TYPE min_angle);
void output_large_angles(const bool flag_internal, const ANGLE_TYPE max_angle);
void output_min_max_angle(const bool flag_internal);
void output_min_max_angle
(const bool flag_internal, const int num_poly_edges);
void output_poly_with_min_max_angles(const bool flag_internal);
void output_poly_with_min_angle(const bool flag_internal);
void output_poly_with_max_angle(const bool flag_internal);
void output_min_max_edge_lengths
(const bool flag_internal, const int num_poly_edges);
void output_min_max_edge_lengths(const bool flag_internal);
void write_nonmanifold_edges();

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
void read_input_file(const char * input_filename);
void memory_exhaustion();
void parse_command_line(int argc, char **argv);
void check_input();
void usage_error(), help_msg();

// PARAMETER TYPE
typedef enum
  {POLYFILE_PARAM, MESH_DIM_PARAM,
   VERTEX_PARAM, SIMPLEX_PARAM, POLY_PARAM,
   VLIST_PARAM, PLIST_PARAM, CONTAINSV_PARAM, CONTAINSE_PARAM,
   MANIFOLD_PARAM, 
   SELFI_PARAM, SELFI_NO_GRID_PARAM, GRID_LENGTH_PARAM,
   MULTIVERT_PARAM,
   MINC_PARAM, MAXC_PARAM, MIN_NUMV_PARAM, MAX_NUMV_PARAM,
   ANGLE_LE_PARAM, ANGLE_GE_PARAM,
   LIST_DUP_PARAM, INTERNAL_PARAM,
   REPORT_DEEP_PARAM, OUT_VALUES_PARAM,
   OUT_MIN_ANGLE_PARAM, OUT_MAX_ANGLE_PARAM, PLOT_ANGLES_PARAM,
   FOR_EACH_TYPE_PARAM, 
   TERSE_PARAM, HELP_PARAM, UNKNOWN_PARAM} PARAMETER;
const char * parameter_string[] = 
  {"-polyfile", "-mesh_dim", "-vertex", "-simplex", "-poly",
   "-vlist", "-plist", "-containsv", "-containse",
   "-manifold", 
   "-selfI", "-selfI_no_grid", "-grid_length",
   "-multivert",
   "-minc", "-maxc", "-min_numv", "-max_numv",
   "-angle_le", "-angle_ge",
   "-list_dup", "-internal",
   "-report_deep", "-out_values", "-out_min_angle", "-out_max_angle",
   "-plot_angles", "-for_each_type",
   "-terse", "-help", "-unknown"};


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  bool passed_all_manifold_tests = false;
  bool passed_boundary_test = false;
  bool flag_self_intersect = false;

  std::set_new_handler(memory_exhaustion);

  parse_command_line(argc, argv);

  try {

    read_input_file(input_filename);

    check_input();

    sort_poly();
    identify_duplicates();

    if (flag_internal && simplex_info_flag) { compute_facet_info(); }

    if (general_info_flag) { output_general_info(); }
    else if (flag_output_min_angle || flag_output_max_angle) {
      output_min_max_angle(flag_internal); 
    }
    if (vertex_info_flag) { output_vertex_info(vertex_index); };
    if (poly_info_flag) { output_poly_info(poly_index); };
    if (simplex_info_flag && !flag_simplex_file) {
      if (simplex_info_flag) { output_poly_info(simplex_index); };
    }
    else {
      if (simplex_info_flag) { output_simplex_info(simplex_index); };
    }

    if (manifold_flag) { 
      identify_nonmanifold();

      if (mesh_info.AreAllZero()) { passed_all_manifold_tests = true; }

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

        bool flag_passed_tests = passed_all_manifold_tests;
        if (mesh_dimension < dimension) {
          flag_passed_tests = 
            (passed_all_manifold_tests && passed_boundary_test);
        }

        if (terse_flag && flag_passed_tests) {
          if (mesh_dimension < dimension) {
            cout << "Passed all manifold and boundary tests." << endl;
          }
          else {
            cout << "Passed all manifold tests." << endl;
          }
        }
        else if (mesh_dimension == 2 && output_filename != NULL) {
          write_nonmanifold_edges();
        }
        else {

          if (terse_flag) {
            if (!passed_all_manifold_tests) 
              { cout << "Failed manifold tests.  "; }
            if (mesh_dimension < dimension) {
              if (flag_report_deep) {
                if (num_deep_boundary_facets != 0) 
                  { cout << "Failed boundary test.";  }
              }
              else {
                if (num_internal_boundary_facets != 0) 
                  { cout << "Failed boundary test.";  }
              }
            }
            cout << endl;
          }

          output_manifold_info();
        };
      };
    };

    if (flag_report_self_intersections) {
      if (flag_simplex_file && num_vert_per_simplex == 3 &&
          dimension == DIM3) {
        flag_self_intersect = output_self_intersections();
      }
    }

    if (vlist_flag) 
      { output_vertex_list(); }

    if (flag_polyfile) {
      if (contains_vertex_flag)
        { output_polytopes(); }

      if (flag_list_duplicate_vertices) {
        output_duplicate_vertices(); 
        cout << endl;
      }

      if (flag_list_duplicate_poly) {
        output_duplicate_poly(); 
        cout << endl;
      }
    }
    else {
      if (contains_vertex_flag || contains_edge_flag)
        { output_simplices(); }
    }

    if (is_mesh_dimension_set || !flag_polyfile) {
      if (mesh_dimension == 2) {

        if (plist_flag) {

          if (flag_internal &&
              contains_boundary_facet.size() != num_poly)
            { identify_nonmanifold(); }

          if (is_min_angle_set || is_max_angle_set ||
              flag_output_min_angle || flag_output_max_angle) {

            if (is_min_angle_set) {
              output_small_angles(flag_internal, angle_le);
            }
            else if (flag_output_min_angle) {
              output_poly_with_min_angle(flag_internal);
            }

            if (is_max_angle_set) {
              output_large_angles(flag_internal, angle_ge);
            }
            else if (flag_output_max_angle) {
              output_poly_with_max_angle(flag_internal);
            }
          }
          else {
            output_poly_with_min_max_angles(flag_internal);
          }
        }
      }
    }

    if (is_mesh_dimension_set || !flag_polyfile) {
      if (mesh_dimension == 2 && flag_plot_angles) {
        ANGLE_TABLE angle_table;
        angle_table.angle.Include();
        angle_table.angle.SetAtIntervals(0, 1);

        compute_polygon_angles(angle_table);

        string output_prefix;
        string output_suffix = ".gplt";

        // create output filename
        string fname = input_filename;

        // remove path from file name
        string prefix, suffix;
        split_string(fname, '/', prefix, suffix);
        if (suffix != "") { fname = suffix; }
        split_string(fname, '.', prefix, suffix);
        if (suffix == "off") {
          output_prefix = prefix;
        }
        else {
          output_prefix = input_filename;
        }

        if (flag_plot_min_polygon_angles) {
          angle_table.HideAllExceptAngleColumn();
          angle_table.min_polygon_angle_freq.Show();
          write_table_gplt
            (output_prefix, "min_poly_angle_freq"+output_suffix, angle_table);
        }

        if (flag_plot_max_polygon_angles) {
          angle_table.HideAllExceptAngleColumn();
          angle_table.max_polygon_angle_freq.Show();
          write_table_gplt
            (output_prefix, "max_poly_angle_freq"+output_suffix, angle_table);
        }

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

  if (manifold_flag && !passed_all_manifold_tests) {
    return(1);
  }
  if (manifold_flag && (mesh_dimension < dimension) && 
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

void read_input_file(const char * input_filename)
{
  const int DIM3(3);
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
    (in, dimension, vertex_coord, num_vertices, 
     polymesh.list_length, polymesh.element, polymesh.first_element);

  in.close();

  num_poly = polymesh.NumPoly();

  if (num_poly == 0) { return; }

  int min_num_vert = 
    *(min_element(polymesh.list_length.begin(), polymesh.list_length.end()));
  int max_num_vert = 
    *(max_element(polymesh.list_length.begin(), polymesh.list_length.end()));

  if (min_num_vert == max_num_vert) {

    num_vert_per_poly = min_num_vert;

    if (!is_mesh_dimension_set) {

      if (dimension == 3 && num_vert_per_poly == 4) {
        cerr << "Unable to determine mesh dimension." << endl;
        cerr << "  Input polytopes could be tetrahedra or quadrilaterals."
             << endl;
        cerr << "Use option -mesh_dim <mdim>." << endl;
        exit(20);
      }

      if (num_vert_per_poly <= dimension+1) {
        // Assume poly are all simplices.
        mesh_dimension = num_vert_per_poly-1;
      }
      else if (num_vert_per_poly == NUMV_PER_CUBE) {
        mesh_dimension = DIM3;
        flag_cube_file = true;
      }
      else {
        mesh_dimension = dimension;
      }
    }

    if (num_vert_per_poly <= mesh_dimension+1) {

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
      const int numv_per_cube = (1 << mesh_dimension);
      if (numv_per_cube != num_vert_per_poly) {
        const int cube_dimension = 
          IJK::compute_cube_dimension_from_num_vertices(num_vert_per_poly);
        cerr << "Usage error.  Mismatch between input file of cubes and -mesh_dim argument." << endl;
        cerr << "  Input file cubes have "
             << numv_per_cube << " vertices and dimension "
             << cube_dimension << "." << endl;
        cerr << "  Argument of -mesh_dim is " << mesh_dimension << "." << endl;
        cerr << "  Change argument of -mesh_dim to " << cube_dimension
             << "." << endl;
        exit(20);
      }
    }

    if (is_mesh_dimension_set) {
      if (num_vert_per_poly <= mesh_dimension) {
        cerr << "Warning:  All polygons in input file are degenerate."
             << endl;
        cerr << "  Argument of -mesh_dim is " << mesh_dimension << "." << endl;
        cerr << "  All polygons have "  << num_vert_per_poly
             << " vertices." << endl;
        cerr << endl;
      }
    }

  }
  else {
    flag_polyfile = true;

    if (!is_mesh_dimension_set) {
      mesh_dimension = dimension-1;

      cerr << "Warning: Unable to determine mesh dimension from input file."
           << endl;
      cerr << "  Assuming mesh dimension is " << mesh_dimension << "." << endl;
      cerr << "  Use option \"-mesh_dim {mdim}\" to specify mesh dimension."
           << endl;
      cerr << endl;
    }
  }

  if (flag_internal) {
    if (!flag_simplex_file && !mesh_dimension == 2) {
      cerr << "Option -internal only implemented for mesh of simplices or for mesh dimension 2." << endl;
      exit(35);
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

void output_general_info()
{
  compute_bounding_box();
  if (!flag_polyfile) {
    num_edges = compute_num_edges();
  }

  cout << "Volume dimension: " << dimension << endl;
  cout << "Mesh dimension: " << mesh_dimension << endl;
  cout << "Number of mesh vertices: " << num_vertices << endl;
  if (flag_polyfile) {
    if (mesh_dimension <= 2) {
      cout << "Number of mesh polygons: " << num_poly << endl;
    }
    else {
      cout << "Number of mesh polytopes: " << num_poly << endl;
    }
  }
  else {
    if (dimension > 2) {
      cout << "Number of mesh edges: " << num_edges << endl;
      cout << "Number of mesh simplices: " << num_simplices << endl;
    }
    else if (dimension == 2) {
      cout << "Number of mesh edges: " << num_simplices << endl;
    }

    if (dimension == 3) {
      cout << "#V - #E + #F = " << num_vertices - num_edges + num_simplices 
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

  flag_output_min_angle = true;
  flag_output_max_angle = true;
  output_min_max_angle(false);

  if (flag_internal) 
    { output_min_max_angle(true); }

  output_min_max_edge_lengths(false);
  if (flag_internal) 
    { output_min_max_edge_lengths(true); }

  cout << "Bounding box: (";
  IJK::print_list(cout, bounding_box.MinCoord(), bounding_box.Dimension());
  cout << " ";
  IJK::print_list(cout, bounding_box.MaxCoord(), bounding_box.Dimension());
  cout << ")" << endl;

  cout << endl;
}

void output_min_max_angle(const bool flag_internal)
{
  output_min_max_angle(flag_internal, 0);

  if (flag_for_each_type) {
    for (int num_poly_edges = 3; num_poly_edges < 10; num_poly_edges++) {
      int npoly = count_num_poly(num_poly_edges);
      if (npoly > 0) {
        output_min_max_angle(flag_internal, num_poly_edges);
      }
    }
  }
}

void output_min_max_angle
(const bool flag_internal, const int num_poly_edges)
{
  string polyname = "polygon";

  if (mesh_dimension == 2) {
    ANGLE_TYPE min_angle, max_angle;
    compute_min_max_polygon_angles
      (flag_internal, num_poly_edges, min_angle, max_angle);

    if (num_poly_edges == 3) { polyname = "triangle"; }
    else if (num_poly_edges == 4) { polyname = "quadrilateral"; }
    else if (num_poly_edges == 5) { polyname = "pentagon"; }
    else if (num_poly_edges == 6) { polyname = "hexagon"; }
    else if (num_poly_edges > 0) { 
      string s;
      val2string(num_poly_edges, s);
      polyname = string("polygon (") + s + " edges)"; };

    if (flag_output_min_angle) {
      cout << "Min ";
      if (flag_internal) { cout << "internal "; }
      cout << polyname << " angle: ";
      cout << min_angle << endl; 
    }
    if (flag_output_max_angle) {
      cout << "Max ";
      if (flag_internal) { cout << "internal "; }
      cout << polyname << " angle: ";
      cout << max_angle << endl; 
    }

    int num_le, num_ge;
    if (is_min_angle_set || is_max_angle_set) {
      compute_num_polygon_angles
        (flag_internal, angle_le, angle_ge, num_le, num_ge);
      if (is_min_angle_set && flag_output_min_angle) {
        if (flag_internal) {
          cout << "Number of internal polygons with angles <= ";
        }
        else {
          cout << "Number of polygons with angles <= ";
        }
        cout << angle_le << ": " << num_le << endl;
      }
      if (is_max_angle_set && flag_output_max_angle) {
        if (flag_internal) {
          cout << "Number of internal polygons with angles >= ";
        }
        else {
          cout << "Number of polygons with angles >= ";
        }
        cout << angle_ge << ": " << num_ge << endl;
      }
    }
  }

}


void output_min_max_edge_lengths
(const bool flag_internal, const int num_poly_edges)
{
  string polyname = "polygon";

  if (is_mesh_dimension_set || !flag_polyfile) {
    if (mesh_dimension == 2) {
      ANGLE_TYPE min_edge_length, max_edge_length;
      compute_min_max_edge_lengths
        (flag_internal, num_poly_edges, min_edge_length, max_edge_length);

      if (num_poly_edges == 3) { polyname = "triangle"; }
      else if (num_poly_edges == 4) { polyname = "quadrilateral"; }
      else if (num_poly_edges == 5) { polyname = "pentagon"; }
      else if (num_poly_edges == 6) { polyname = "hexagon"; }
      else if (num_poly_edges > 0) { 
        string s;
        val2string(num_poly_edges, s);
        polyname = string("polygon (") + s + " edges)"; };

      if (flag_output_min_angle) {
        cout << "Min ";
        if (flag_internal) { cout << "internal "; }
        cout << polyname << " edge length: ";
        cout << min_edge_length << endl; 
      }
      if (flag_output_max_angle) {
        cout << "Max ";
        if (flag_internal) { cout << "internal "; }
        cout << polyname << " edge length: ";
        cout << max_edge_length << endl; 
      }
    }
  }

}

void output_min_max_edge_lengths(const bool flag_internal)
{
  output_min_max_edge_lengths(flag_internal, 0);

  if (flag_for_each_type) {
    for (int num_poly_edges = 3; num_poly_edges < 10; num_poly_edges++) {
      int npoly = count_num_poly(num_poly_edges);
      if (npoly > 0) {
        output_min_max_edge_lengths(flag_internal, num_poly_edges);
      }
    }
  }
}

void output_vertex_info(const int vertex_index)
{
  cout << "Vertex: " << vertex_index << endl;

  if (vertex_index < 0 || vertex_index >= num_vertices)
    {
      cout << "  Illegal vertex index " << vertex_index
           << ".  Vertex index should be in range["
           << 0 << "," << num_vertices-1 << "]." << endl;
      return;
    };

  cout << "  Coordinates: ";
  IJK::print_list(cout, vertex_coord+vertex_index*dimension, dimension);
  cout << endl;

  cout << endl;
}


void output_simplex_info(const int simplex_index)
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
    IJK::print_list(cout, vertex_coord+iv*dimension, dimension);
    cout << endl;
  }

  cout << endl;
}

void output_poly_info(const int poly_index)
{
  cout << "Poly: " << poly_index << endl;

  if (poly_index < 0 || poly_index >= num_poly) {
    cout << "  Illegal poly index.  Poly index should be in range["
         << 0 << "," << num_poly-1 << "]." << endl;
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


bool coord_match(const COORD_TYPE * vertex_coord,
                 const int iv0, const int iv1)
{
  const COORD_TYPE * vcoord0 = vertex_coord + iv0*dimension;
  const COORD_TYPE * vcoord1 = vertex_coord + iv1*dimension;
  for (int d = 0; d < dimension; d++) {
    if (vcoord0[d] != vcoord1[d]) 
      { return(false); }
  }

  return(true);
}


void output_manifold_info()
{
  const int NUMV_PER_CUBE = 8;
  const int NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;
  const char * poly_str = "polytopes";
  const char * facet_str = "facets";
  bool flag_newline(false);

  if (mesh_dimension == 2) { 
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
      output_degenerate_poly();
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
      output_duplicate_poly();
      cout << endl;
      flag_newline = true;
    }
  }

  if (mesh_dimension == 2 || flag_simplex_file || flag_cube_file) {

    int numv_per_facet;
    if (mesh_dimension == 2 || flag_simplex_file) 
      { numv_per_facet = mesh_dimension;  }
    else 
      { numv_per_facet = NUMV_PER_CUBE_FACET; }

    const int num_nonmanifold_facets = 
      nonmanifold_facet_vert.size()/numv_per_facet;

    if (num_nonmanifold_facets == 0) {
      if (!terse_flag) {
        cout << "No non-manifold " << facet_str << "." << endl;
        flag_newline = false;
      }
    }
    else {
      if (!flag_newline && !terse_flag) { cout << endl; }
      cout << "Num non-manifold " << facet_str << ":  " 
           << num_nonmanifold_facets << endl;
      if (!terse_flag) {
        output_nonmanifold_facets();
        cout << endl;
        flag_newline = true;
      }
    }

    if (flag_cube_file && mesh_dimension > 2) {
      // If mesh_dimension == 2, then edges are reported as facets.
      // Non-manifold edge detection only implemented for cubes.

      if (mesh_info.num_nonmanifold_edges == 0) {
        if (!terse_flag) {
          cout << "No non-manifold edges." << endl;
          flag_newline = false;
        }
      }
      else {
        if (!flag_newline && !terse_flag) { cout << endl; }
        cout << "Num non-manifold edges " 
             << mesh_info.num_nonmanifold_edges << endl;
        if (!terse_flag) {
          output_nonmanifold_edges();
          cout << endl;
          flag_newline = true;
        }

      }
    }

    if (mesh_info.num_nonmanifold_vertices == 0) {
      if (!terse_flag) {
        cout << "No non-manifold vertices." << endl;
        flag_newline = false;
      }
    }
    else {
      if (!flag_newline && !terse_flag) { cout << endl; }

      int num_deep = count_deep_vertices(nonmanifold_vert_list);
      cout << "Num non-manifold vertices:  " 
           << mesh_info.num_nonmanifold_vertices << endl;

      if (flag_report_deep) {
        cout << "Num non-manifold vertices at least " << contract_margin
             << " from bounding box boundary: " << num_deep << endl;
      }

      if (!terse_flag) {
        output_nonmanifold_vertices(); 
        cout << endl;
        flag_newline = true;
      }
    }

    if (mesh_dimension < dimension)
      { output_internal_boundary_facets(); }
  }
  else {
    if (!terse_flag) {
      cout << "Unable to determine " << poly_str 
           << " facets to check manifold/boundary conditions." << endl;

    }
  }

}


void output_degenerate_poly()
{
  if (mesh_info.num_poly_with_duplicate_vertices > 0) {
    cout << "Degenerate poly:" << endl;
    for (int ipoly = 0; ipoly < is_degenerate.size(); ipoly++) {
      if (is_degenerate[ipoly]) {
        cout << "  Poly " << ipoly << ": ";
        print_list(cout, polymesh.VertexList(ipoly), 
                   polymesh.NumPolyVert(ipoly));
        cout << endl;
      }
    }
  }
}


void output_duplicate_poly()
{
  if (mesh_info.num_duplicate_poly > 0) {
    cout << "Duplicate poly:" << endl;
    for (int ipoly = 0; ipoly < is_duplicate.size(); ipoly++) {
      if (is_duplicate[ipoly]) {
        cout << "  Poly " << ipoly << ": ";
        print_list(cout, polymesh.VertexList(ipoly), 
                   polymesh.NumPolyVert(ipoly));
        cout << endl;
      }
    }
  }
}


void output_nonmanifold_facets()
{
  const int NUMV_PER_CUBE = 8;
  const int NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;
  const char * poly_str = "Polytopes";
  const char * facet_str = "facets";
  int numv_per_facet = mesh_dimension;
  const int num_nonmanifold_facets = 
    nonmanifold_facet_vert.size()/numv_per_facet;

  if (mesh_dimension == 2 || flag_simplex_file) 
    { numv_per_facet = mesh_dimension;  }
  else 
    { numv_per_facet = NUMV_PER_CUBE_FACET; }

  if (mesh_dimension == 2) {
    poly_str = "Polygons";
    facet_str = "edges"; 
  }

  cout << "Non-manifold " << facet_str << ":" << endl;

  for (int jf = 0; jf < num_nonmanifold_facets; jf++) {
    cout << "  ";
    print_list(cout, &nonmanifold_facet_vert.front()+jf*numv_per_facet,
               numv_per_facet);

    if (mesh_dimension == 2) {
      cout << "  (";
      for (int i = 0; i < numv_per_facet; i++) {
        int iv = nonmanifold_facet_vert[jf*numv_per_facet+i];
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
    if (contains_nonmanifold_facet[ipoly]) {
      cout << "  " << ipoly;
      num_output++;

      if (num_output%10 == 0)
        cout << endl;
    }
  }
  if (num_output%10 != 0) { cout << endl; };
}


void output_nonmanifold_edges()
{
  const char * poly_str = "Cubes";
  const char * edge_str = "edges";
  const int num_nonmanifold_edges = nonmanifold_edge_vert.size()/2;

  cout << "Non-manifold " << edge_str << ":" << endl;

  for (int je = 0; je < num_nonmanifold_edges; je++) {
    cout << "  ";
    print_list(cout, &nonmanifold_edge_vert.front()+je*2, 2);
    cout << endl;
  }
  cout << endl;

  /* DEBUG
  cout << poly_str << " containing non-manifold " << facet_str 
       << ": " << endl;

  int num_output = 0;
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {
    if (contains_nonmanifold_edge[ipoly]) {
      cout << "  " << ipoly;
      num_output++;

      if (num_output%10 == 0)
        cout << endl;
    }
  }
  if (num_output%10 != 0) { cout << endl; };
  */
}


void output_internal_boundary_facets()
{
  int numv_per_facet;
  IJK::PROCEDURE_ERROR error("output_internal_boundary_facet");

  if (flag_cube_file) {
    CUBE_TYPE cube(mesh_dimension);
    numv_per_facet = cube.NumFacetVertices();
  }
  else if (flag_simplex_file) {
    numv_per_facet = mesh_dimension;
  }
  else if (mesh_dimension == 2) {
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


void output_duplicate_vertices()
{
  for (int jpoly = 0; jpoly < num_poly; jpoly++) {

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
          if (is_coord_equal(dimension, vertex_coord+dimension*pvert[k],
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
  const int numv_per_facet = mesh_dimension;
  const int num_nonmanifold_facets = 
    nonmanifold_facet_vert.size()/numv_per_facet;
  const int num_deep_vertices = count_deep_vertices(nonmanifold_vert_list);

  cout << nonmanifold_vert_list.size() << " "
       << num_deep_vertices << " "
       << num_nonmanifold_facets << " "
       << num_internal_boundary_facets << " "
       << num_deep_boundary_facets << endl;
}


void write_nonmanifold_edges()
{
  const int numv_per_facet = mesh_dimension;
  const int num_nonmanifold_edges =
    nonmanifold_facet_vert.size()/numv_per_facet;

  if (mesh_dimension != 2) {
    cerr << "Only able to write non-manifold edges for mesh dimension two."
         << endl;
    return;
  }

  if (num_nonmanifold_edges == 0) {

    cout << "No non-manifold edges." << endl;
    return;
  }

  ofstream output_file;
  output_file.open(output_filename, ios::out);

  float rgba[4] = { 1, 0, 1, 1};

  IJK::ijkoutColorLINE
    (output_file, dimension, vertex_coord, num_vertices,
     vector2pointer(nonmanifold_facet_vert), num_nonmanifold_edges, rgba);

  output_file.close();
}

void output_nonmanifold_vertices()
{
  if (!terse_flag) {
    int num_output = 0;
    for (int i = 0; i < nonmanifold_vert_list.size(); i++) {
      cout << "  " << nonmanifold_vert_list[i];
      num_output++;
      if (num_output%10 == 0) { cout << endl; };
    }
    if (num_output%10 != 0) { cout << endl; };
  }
}

void output_poly(const int ipoly)
// output polytope
{
  cout << "(";
  for (int j = 0; j < polymesh.NumPolyVert(ipoly); j++) {
    cout << polymesh.Vertex(ipoly,j);
    if (j+1 < polymesh.NumPolyVert(ipoly))
      { cout << "  "; }
  }
  cout << ")";
}

void output_simplex(const int * simplex)
// output simplex
{
  cout << "(";
  for (int j = 0; j < num_vert_per_poly; j++) {
    cout << *(simplex+j);
    if (j+1 < num_vert_per_poly)
      { cout << "  "; }
  }
  cout << ")";
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
  for (int iv = 0; iv < num_vertices; iv++) {

    const COORD_TYPE * vcoord = vertex_coord + iv*dimension;

    if (is_min_coord_set) {
      if (!region_contains
          (dimension, vcoord, &(min_coord[0]))) {
        continue;
      }
    }

    if (is_max_coord_set) {
      if (!region_contains
          (dimension, &(max_coord[0]), vcoord)) {
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

/// Output polygons with small angles
void output_small_angles
(const bool flag_internal, const ANGLE_TYPE min_angle)
{
  IJK::PROCEDURE_ERROR error("output_small_angles");

  if (mesh_dimension != 2) {
    error.AddMessage("Programming error. Mesh dimension should be 2.");
    error.AddMessage("  Mesh dimension = ", mesh_dimension, ".");
    throw error;
  }

  if (flag_internal && contains_boundary_facet.size() != num_poly) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }

  int num_out = 0;

  COORD_TYPE max_cos = cos(min_angle*M_PI/180.0);
  COORD_TYPE min_cos_i, max_cos_i;
  int num_angle;

  if (flag_internal) {
    cout << "Internal polygons with small angles: " << endl;
  }
  else {
    cout << "Polygons with small angles: " << endl;
  }

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (is_degenerate[ipoly]) { continue; }

    if (flag_internal) 
      { if (contains_boundary_facet[ipoly]) { continue; } }

    compute_min_max_polygon_cos(ipoly, min_cos_i, max_cos_i, num_angle);
    if (num_angle > 0 && max_cos_i >= max_cos) {
      num_out++;

      const int * pvert = polymesh.VertexList(ipoly);

      // Note: min_angle_i is acos of max_cos_i.
      ANGLE_TYPE min_angle_i = std::acos(max_cos_i) * 180.0/M_PI;

      cout << "  Poly " << ipoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
      cout << "  Min angle: " << min_angle_i;
      cout << endl;
    }
  }

  if (num_out == 0) {
    cout << "  No polytopes with angles <= " << min_angle << endl;
  }
  
}

/// Output polygons with large angles
void output_large_angles
(const bool flag_internal, const ANGLE_TYPE max_angle)
{
  IJK::PROCEDURE_ERROR error("output_large_angles");

  if (mesh_dimension != 2) {
    error.AddMessage("Programming error. Mesh dimension should be 2.");
    error.AddMessage("  Mesh dimension = ", mesh_dimension, ".");
    throw error;
  }

  if (flag_internal && contains_boundary_facet.size() != num_poly) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }

  int num_out = 0;

  COORD_TYPE min_cos = cos(max_angle*M_PI/180.0);
  COORD_TYPE min_cos_i, max_cos_i;
  int num_angle;

  if (flag_internal) {
    cout << "Internal polygons with large angles: " << endl;
  }
  else {
    cout << "Polygons with large angles: " << endl;
  }

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (is_degenerate[ipoly]) { continue; }

    if (flag_internal) 
      { if (contains_boundary_facet[ipoly]) { continue; } }

    compute_min_max_polygon_cos(ipoly, min_cos_i, max_cos_i, num_angle);
    if (num_angle > 0 && min_cos_i < min_cos) {
      num_out++;

      const int * pvert = polymesh.VertexList(ipoly);

      // Note: max_angle_i is acos of min_cos_i.
      ANGLE_TYPE max_angle_i = std::acos(min_cos_i) * 180.0/M_PI;

      cout << "  Poly " << ipoly << ": ";
      print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
      cout << "  Max angle: " << max_angle_i;
      cout << endl;
    }
  }

  if (num_out == 0) {
    cout << "  No polytopes with angles >= " << max_angle << endl;
  }  
}

/// Output polygons with minimum and maximum angles
void output_poly_with_min_max_angles(const bool flag_internal)
{
  output_poly_with_min_angle(flag_internal);
  output_poly_with_max_angle(flag_internal);
}

/// Output polygons with minimum angle
void output_poly_with_min_angle(const bool flag_internal)
{
  ANGLE_TYPE min_angle, max_angle;
  int poly_with_min_angle, poly_with_max_angle;
  int num_angle;
  IJK::PROCEDURE_ERROR error("output_poly_with_min_min_angle");

  if (flag_internal && contains_boundary_facet.size() != num_poly) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }

  compute_min_max_polygon_angles
    (flag_internal, min_angle, max_angle, 
     poly_with_min_angle, poly_with_max_angle);

  if (flag_internal) {
    cout << "Internal polygons with minimum angle ";
  }
  else {
    cout << "Polygons with minimum angle ";
  }
  cout << min_angle << ":" << endl;

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (flag_internal) 
      { if (contains_boundary_facet[ipoly]) { continue; } }

    COORD_TYPE min_cos_i, max_cos_i;
    compute_min_max_polygon_cos(ipoly, min_cos_i, max_cos_i, num_angle);

    if (num_angle > 0) {

      // Note: min_angle_i is acos of max_cos_i.
      ANGLE_TYPE min_angle_i = std::acos(max_cos_i) * 180.0/M_PI;

      if (ipoly == poly_with_min_angle || min_angle_i == min_angle) {

        const int * pvert = polymesh.VertexList(ipoly);

        cout << "  Poly " << ipoly << ": ";
        print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
        cout << "  Min angle: " << min_angle_i;
        cout << endl;
      }
    }
  }
}

/// Output polygons with maximum angle
void output_poly_with_max_angle(const bool flag_internal)
{
  ANGLE_TYPE min_angle, max_angle;
  int poly_with_min_angle, poly_with_max_angle;
  int num_angle;
  IJK::PROCEDURE_ERROR error("output_poly_with_max_angle");

  if (flag_internal && contains_boundary_facet.size() != num_poly) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }

  compute_min_max_polygon_angles
    (flag_internal, min_angle, max_angle, 
     poly_with_min_angle, poly_with_max_angle);

  if (flag_internal) {
    cout << "Internal polygons with maximum angle ";
  }
  else {
    cout << "Polygons with maximum angle ";
  }
  cout << max_angle << ":" << endl;

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (flag_internal) 
      { if (contains_boundary_facet[ipoly]) { continue; } }

    COORD_TYPE min_cos_i, max_cos_i;
    compute_min_max_polygon_cos(ipoly, min_cos_i, max_cos_i, num_angle);

    if (num_angle > 0) {

      // Note: max_angle_i is acos of min_cos_i.
      ANGLE_TYPE max_angle_i = std::acos(min_cos_i) * 180.0/M_PI;

      if (ipoly == poly_with_max_angle || max_angle_i == max_angle) {

        const int * pvert = polymesh.VertexList(ipoly);

        cout << "  Poly " << ipoly << ": ";
        print_list(cout, pvert, polymesh.NumPolyVert(ipoly));
        cout << "  Max angle: " << max_angle_i;
        cout << endl;
      }
    }
  }

}


// **************************************************
// Write Tables
// **************************************************

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


// **************************************************
// COMPUTE INFO ROUTINES
// **************************************************

void compute_bounding_box()
{
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
  const int NUMV_PER_CUBE = 8;
  const int NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;

  compute_bounding_box();
  if (flag_simplex_file || mesh_dimension == 2) {
    identify_nonmanifold_and_boundary_facets();

    const int num_vert_per_facet = mesh_dimension;
    mesh_info.num_nonmanifold_facets =
      (nonmanifold_facet_vert.size())/num_vert_per_facet;
  }
  else if (flag_cube_file) {
    identify_nonmanifold_and_boundary_facets();
    mesh_info.num_nonmanifold_facets =
      (nonmanifold_facet_vert.size())/NUMV_PER_CUBE_FACET;
  }
}



// **************************************************
// COMPUTE ANGLE ROUTINES
// **************************************************

void compute_cos_angle
(const int iv0, const int iv1, const int iv2,
 double & cos_angle, bool & flag_duplicate_point);

/// Compute min/max polygon angles.
/// @pre mesh_dimension == 2.
void compute_min_max_polygon_angles
(const bool flag_internal, ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  compute_min_max_polygon_angles(flag_internal, 0, min_angle, max_angle);
}

/// Compute min/max polygon angles.
/// @pre mesh_dimension == 2.
void compute_min_max_polygon_angles
(const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_polygon_angles
    (flag_internal, num_poly_edges, min_angle, max_angle, 
     poly_with_min_angle, poly_with_max_angle);
}

/// Compute min/max polygon angles.
/// @pre mesh_dimension == 2.
void compute_min_max_polygon_angles
(const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  compute_min_max_polygon_angles
    (flag_internal, 0, min_angle, max_angle, 
     poly_with_min_angle, poly_with_max_angle);
}

/// Compute min/max polygon angles.
/// @param flag_internal If true, compute angles for interior polygons.
/// @param num_poly_edges If num_poly_edges > 0, compute angles only
///          for polygons with num_poly_edges.
/// @param num_poly_edges If num_poly_edges = 0, compute angles 
///          for all polygons
/// @pre mesh_dimension == 2.
void compute_min_max_polygon_angles
(const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  IJK::PROCEDURE_ERROR error("compute_min_max_polygon_angles");

  if (mesh_dimension != 2) {
    error.AddMessage("Programming error. Mesh dimension should be 2.");
    error.AddMessage("  Mesh dimension = ", mesh_dimension, ".");
    throw error;
  }

  if (flag_internal && contains_boundary_facet.size() != num_poly) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }

  COORD_TYPE min_cos = 1;
  COORD_TYPE max_cos = -1;
  poly_with_min_angle = 0;
  poly_with_max_angle = 0;

  int num_angle;
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (is_degenerate[ipoly]) { continue; }

    if (num_poly_edges > 0) {
      if (polymesh.NumPolyVert(ipoly) != num_poly_edges) 
        { continue; }
    }

    if (flag_internal) 
      { if (contains_boundary_facet[ipoly]) { continue; } }

    COORD_TYPE min_cos_i, max_cos_i;
    compute_min_max_polygon_cos(ipoly, min_cos_i, max_cos_i, num_angle);

    if (num_angle > 0) {
      if (min_cos_i < min_cos) { 
        min_cos = min_cos_i; 
        // Note: max_angle is acos of min_cos.
        poly_with_max_angle = ipoly;
      }
      if (max_cos_i > max_cos) { 
        max_cos = max_cos_i; 
        // Note: min_angle is acos of max_cos.
        poly_with_min_angle = ipoly;
      }
    }
  }


  // Note: min_angle is acos of max_cos.
  // Note: max_angle is acos of min_cos.
  min_angle = std::acos(max_cos) * 180.0/M_PI;
  max_angle = std::acos(min_cos) * 180.0/M_PI;
}

/// Compute number of angles less than or equal to min_angle and
///   greater than or equal to max_angle.
/// @pre mesh_dimension == 2.
void compute_num_polygon_angles
(const bool flag_internal, 
 const ANGLE_TYPE & min_angle, const ANGLE_TYPE & max_angle,
 int & num_le, int & num_ge)
{
  IJK::PROCEDURE_ERROR error("compute_num_polygon_angles");

  // Initialize to zero.
  num_le = 0;
  num_ge = 0;

  if (mesh_dimension != 2) {
    error.AddMessage("Programming error.  Mesh dimension should be 2.");
    error.AddMessage("  Mesh dimension = ", mesh_dimension, ".");
    throw error;
  }

  if (flag_internal && contains_boundary_facet.size() != num_poly) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }

  // Note: max_cos is cosine of min_angle.
  //       min_cos is cosine of max_angle. 
  COORD_TYPE max_cos = cos(min_angle*M_PI/180.0);
  COORD_TYPE min_cos = cos(max_angle*M_PI/180.0);
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (flag_internal) 
      { if (contains_boundary_facet[ipoly]) { continue; } }

    COORD_TYPE min_cos_i, max_cos_i;
    int num_angle;
    compute_min_max_polygon_cos(ipoly, min_cos_i, max_cos_i, num_angle);

    if (num_angle > 0) {
      if (min_cos_i <= min_cos) { num_ge++; }
      if (max_cos_i >= max_cos) { num_le++; }
    }

  }

}

/// Compute min/max cosine of polygon angles for polygon ipoly.
/// @param[out] num_angle Number of angles evaluated.
/// @pre mesh_dimension == 2.
void compute_min_max_polygon_cos
(const int ipoly, COORD_TYPE & min_cos, COORD_TYPE & max_cos,
 int & num_angle)
{
  min_cos = 1;
  max_cos = -1;
  num_angle = 0;
  // Ignore degenerate polygons with 1 or 2 vertices.
  if (polymesh.NumPolyVert(ipoly) > 2) {
    int nv = polymesh.NumPolyVert(ipoly);
    for (int i0 = 0; i0 < nv; i0++) {
      int i1 = (i0+1)%nv;
      int i2 = (i0+2)%nv;
      int iv0 = polymesh.Vertex(ipoly,i0);
      int iv1 = polymesh.Vertex(ipoly,i1);
      int iv2 = polymesh.Vertex(ipoly,i2);
      double cos_angle;
      bool flag_duplicate_point;
      compute_cos_angle(iv0, iv1, iv2, cos_angle, flag_duplicate_point);

      if (!flag_duplicate_point) {
        num_angle++;
        if (cos_angle < min_cos) { min_cos = cos_angle; }
        if (cos_angle > max_cos) { max_cos = cos_angle; }
      }
    }
  }

}

/// Compute min/max polygon angles for polygon ipoly.
/// @param[out] num_angle Number of angles evaluated.
/// @pre mesh_dimension == 2.
void compute_min_max_polygon_angles
(const int ipoly, ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & num_angle)
{
  COORD_TYPE min_cos, max_cos;

  compute_min_max_polygon_cos(ipoly, min_cos, max_cos, num_angle);

  min_angle = std::acos(max_cos) * 180.0/M_PI;
  max_angle = std::acos(min_cos) * 180.0/M_PI;
}

/// Compute cosine of the angle between (iv0-iv1) and (iv2-iv1).
void compute_cos_angle
(const int iv0, const int iv1, const int iv2,
 double & cos_angle, bool & flag_duplicate_point)
{
  ARRAY<COORD_TYPE> w0(dimension), w1(dimension);

  subtract_coord
    (dimension, vertex_coord+dimension*iv0, vertex_coord+dimension*iv1,
     w0.Ptr());
  subtract_coord
    (dimension, vertex_coord+dimension*iv2, vertex_coord+dimension*iv1,
     w1.Ptr());
  IJK::compute_cos_angle(dimension, w0.PtrConst(), w1.PtrConst(),
                         cos_angle, flag_duplicate_point);
}


// **************************************************
// COMPUTE EDGE LENGTH ROUTINES
// **************************************************

void compute_min_max_polygon_edge_length_squared
(const int ipoly, COORD_TYPE & min_edge_length_squared, 
 COORD_TYPE & max_edge_length_squared, int & num_edge);


/// Compute min/max edge lengths.
/// Skip polygons with duplicate vertices.
/// @param flag_internal If true, compute edge lengths for interior polygons.
/// @param num_poly_edges If num_poly_edges > 0, compute angles only
///          for polygons with num_poly_edges.
/// @param num_poly_edges If num_poly_edges = 0, compute angles 
///          for all polygons
/// @pre mesh_dimension == 2.
void compute_min_max_edge_lengths
(const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_edge_length, ANGLE_TYPE & max_edge_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  COORD_TYPE min_length_squared, max_length_squared;
  bool is_length_set;
  IJK::PROCEDURE_ERROR error("compute_min_max_edge_lengths");

  if (mesh_dimension != 2) {
    error.AddMessage("Programming error. Mesh dimension should be 2.");
    error.AddMessage("  Mesh dimension = ", mesh_dimension, ".");
    throw error;
  }

  if (flag_internal && contains_boundary_facet.size() != num_poly) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }

  is_length_set = false;
  min_length_squared = 0;
  max_length_squared = 0;
  poly_with_min_edge_length = 0;
  poly_with_max_edge_length = 0;

  int num_length;
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (is_degenerate[ipoly]) { continue; }

    if (num_poly_edges > 0) {
      if (polymesh.NumPolyVert(ipoly) != num_poly_edges) 
        { continue; }
    }

    if (flag_internal) 
      { if (contains_boundary_facet[ipoly]) { continue; } }

    COORD_TYPE min_length_squared_i, max_length_squared_i;
    int num_edges;
    compute_min_max_polygon_edge_length_squared
      (ipoly, min_length_squared_i, max_length_squared_i, num_edges);

    if (num_edges > 0) {
      if (!is_length_set || min_length_squared_i < min_length_squared) { 
        min_length_squared = min_length_squared_i; 
        poly_with_min_edge_length = ipoly;
      }
      if (!is_length_set || max_length_squared_i > max_length_squared) { 
        max_length_squared = max_length_squared_i;
        poly_with_max_edge_length = ipoly;
      }
      is_length_set = true;
    }
  }

  if (min_length_squared > 0) 
    { min_edge_length = std::sqrt(min_length_squared); }
  else
    { min_edge_length = 0; }

  if (max_length_squared > 0) 
    { max_edge_length = std::sqrt(max_length_squared); }
  else
    { max_edge_length = 0; }
}

void compute_min_max_edge_lengths
(const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_edge_length, ANGLE_TYPE & max_edge_length)
{
  int poly_with_min_edge_length, poly_with_max_edge_length;

  compute_min_max_edge_lengths
    (flag_internal, num_poly_edges, min_edge_length, max_edge_length,
     poly_with_min_edge_length, poly_with_max_edge_length);
}

/// Compute min/max edge length of polygon ipoly.
/// @param[out] num_edges Number of edges evaluated.
/// @pre mesh_dimension == 2.
void compute_min_max_polygon_edge_length_squared
(const int ipoly, COORD_TYPE & min_edge_length_squared, 
 COORD_TYPE & max_edge_length_squared, int & num_edge)
{
  COORD_TYPE edge_length_squared;

  min_edge_length_squared = 0;
  max_edge_length_squared = 0;
  num_edges = 0;
  // Ignore degenerate polygons with 1 or 2 vertices.
  if (polymesh.NumPolyVert(ipoly) > 2) {
    int nv = polymesh.NumPolyVert(ipoly);

    for (int i0 = 0; i0 < nv; i0++) {
      int i1 = (i0+1)%nv;
      const COORD_TYPE * v0 = vertex_coord+polymesh.Vertex(ipoly,i0)*dimension;
      const COORD_TYPE * v1 = vertex_coord+polymesh.Vertex(ipoly,i1)*dimension;

      IJK::compute_distance_squared(dimension, v0, v1, edge_length_squared);

      if (i0 == 0 || edge_length_squared < min_edge_length_squared)
        { min_edge_length_squared = edge_length_squared; }
      if (i0 == 0 || edge_length_squared > max_edge_length_squared)
        { max_edge_length_squared = edge_length_squared; }
    }

    num_edge = nv;
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

  // initialize is_degenerate[ipoly] to false for all polytopes
  // a polytope with duplicate vertices is degenerate
  is_degenerate.assign(num_poly, false);

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {
    if (polymesh_sorted.NumPolyVert(ipoly) < 2) { continue; }

    int iv0 = polymesh_sorted.Vertex(ipoly, 0);
    for (int j = 1; j < polymesh_sorted.NumPolyVert(ipoly); j++) {
      int iv1 = polymesh_sorted.Vertex(ipoly, j);
      if (iv0 == iv1) {
        is_degenerate[ipoly] = true;
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

  // initialize is_duplicate[ipoly] to false for all polytopes
  is_duplicate.assign(num_poly, false);

  for (int i0 = 0; i0+1 < sorted_poly.size(); i0++) {
    int ipoly0 = sorted_poly[i0];
    int ipoly1 = sorted_poly[i0+1];

    if (!polymesh_lt(ipoly0,ipoly1) && !polymesh_lt(ipoly1,ipoly0)) {
      is_duplicate[ipoly0] = true;
      is_duplicate[ipoly1] = true;
    }
  }

  for (int i = 0; i < is_duplicate.size(); i++) {
    if (is_duplicate[i]) { num_duplicate_poly++; }
  }

  return(num_duplicate_poly);
}


// **************************************************
// MANIFOLD ROUTINES
// **************************************************

void identify_nonmanifold()
{
  if (flag_simplex_file || mesh_dimension == 2 || flag_cube_file) {
    compute_facet_info();

    if (flag_cube_file && mesh_dimension > 2)
      { mesh_info.num_nonmanifold_edges = identify_nonmanifold_edges(); }
    mesh_info.num_nonmanifold_vertices = identify_nonmanifold_vertices();
    mesh_info.num_deep_nonmanifold_vertices = 
      count_deep_vertices(nonmanifold_vert_list);
  }
}


bool simplex_equals(const int * simplex_vert, const int num_simplices,
                    const int is0, const int is1)
// return true if is0 equals is1
// assumes simplex vertices are listed in sorted order
{
  const int numv_per_simplex = num_vert_per_poly;

  for (int k = 0; k < numv_per_simplex; k++) {
    if (simplex_vert[is0*numv_per_simplex+k] !=
        simplex_vert[is1*numv_per_simplex+k])
      return(false);
  }

  return(true);
}

// Return true if jf0 equals jf1.
// @pre facet vertices are listed in sorted order.
bool facet_equals
(const int * facet_vert, const int numv_per_facet, const int jf0, const int jf1)
{
  for (int k = 0; k < numv_per_facet; k++) {
    if (facet_vert[jf0*numv_per_facet+k] !=
        facet_vert[jf1*numv_per_facet+k])
      return(false);
  }

  return(true);
}

void store_polygon_edges
(const POLYMESH_TYPE & polymesh, const int ipoly,
 int & k, int & num_proper_facets, vector<int> & facet_vert_list)
{
  for (int j0 = 0; j0 < polymesh.NumPolyVert(ipoly); j0++) {
    const int j1 = (j0+1)%polymesh.NumPolyVert(ipoly);
    facet_vert_list[k] = polymesh.Vertex(ipoly, j0);
    facet_vert_list[k+1] = polymesh.Vertex(ipoly, j1);
    if (facet_vert_list[k] > facet_vert_list[k+1])
      { std::swap(facet_vert_list[k], facet_vert_list[k+1]); }
    num_proper_facets++;
    k += 2;
  }
}


void store_simplex_facets
(const POLYMESH_TYPE & polymesh, const int ipoly,
 const int num_facets_per_simplex,
 int & k, int & num_proper_facets, vector<int> & facet_vert_list)
{
  for (int jf = 0; jf < num_facets_per_simplex; jf++) {
    for (int jv = 0; jv < num_vert_per_simplex; jv++) {
      if (jv != jf) {
        facet_vert_list[k] = polymesh.Vertex(ipoly, jv);
        k++;
      }
    }
    num_proper_facets++;
  }
}


void store_cube_facets
(const POLYMESH_TYPE & polymesh, const int ipoly, const CUBE_TYPE & cube,
 int & k, int & num_proper_facets, vector<int> & facet_vert_list)
{  
  for (int jf = 0; jf < cube.NumFacets(); jf++) {
    for (int jv = 0; jv < cube.NumFacetVertices(); jv++) {
      int jv2 = cube.FacetVertex(jf, jv);
      facet_vert_list[k] = polymesh.Vertex(ipoly, jv2);
      k++;
    }
    num_proper_facets++;
  }
}


void set_in_nonmanifold_facet()
{
  // initialize in_nonmanifold_facet[iv] to false for all vertices iv
  in_nonmanifold_facet.assign(num_vertices, false);

  for (int i = 0; i < nonmanifold_facet_vert.size(); i++) {
    int iv = nonmanifold_facet_vert[i];
    in_nonmanifold_facet[iv] = true;
  }
}


void set_in_nonmanifold_edge()
{
  // initialize in_nonmanifold_edge[iv] to false for all vertices iv
  in_nonmanifold_edge.assign(num_vertices, false);

  for (int i = 0; i < nonmanifold_edge_vert.size(); i++) {
    int iv = nonmanifold_edge_vert[i];
    in_nonmanifold_edge[iv] = true;
  }
}


// Identify non-manifold facets and boundary facets.
// A facet is non-manifold if it is in more than two polytope.
// A facet is boundary if it is in only one polytope.
void identify_nonmanifold_and_boundary_facets()
{
  const int NUMV_PER_CUBE = 8;
  const int NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;
  const int NUM_FACETS_PER_CUBE = 6;
  const int num_facets_per_simplex = mesh_dimension+1;
  const int num_vert_per_simplex = mesh_dimension+1;
  const CUBE_TYPE cube(DIM3);
  int numv_per_facet = mesh_dimension;
  int num_facets;
  vector<int> facet_vert_list;
  vector<int> first_poly_facet(num_poly);
  vector<bool> is_facet_nonmanifold;
  vector<bool> is_facet_boundary;
  int num_proper_facets = 0;
  IJK::PROCEDURE_ERROR error("identify_nonmanifold_and_boundary_facets");

  if (mesh_dimension == 2) {
    num_facets = sum_num_poly_vert(polymesh);
  }
  else if (flag_simplex_file) {
    num_facets = polymesh.NumPoly()*num_facets_per_simplex;
  }
  else if (flag_cube_file) {
    // All mesh elements are cubes.
    num_facets = polymesh.NumPoly()*NUM_FACETS_PER_CUBE;
    numv_per_facet = NUMV_PER_CUBE_FACET;
  }
  else {
    error.AddMessage
      ("Programming error.  Mesh is not dimension 2 and not a mesh of simplices.");
    error.AddMessage("  Unable to determine polytope facets.");
    throw error;
  }

  facet_vert_list.resize(num_facets*numv_per_facet);

  // store facets in facet_vert_list
  int k = 0;
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {
    first_poly_facet[ipoly] = 0;
    // ignore degenerate polytopes and duplicate polytopes
    if (!is_degenerate[ipoly] && !is_duplicate[ipoly]) {
      first_poly_facet[ipoly] = num_proper_facets;

      if (mesh_dimension == 2) {
        store_polygon_edges
          (polymesh, ipoly, k, num_proper_facets, facet_vert_list);
      }
      else if (flag_simplex_file) {
        store_simplex_facets
          (polymesh, ipoly, num_facets_per_simplex, 
           k, num_proper_facets, facet_vert_list);
      }
      else if (flag_cube_file) {
        store_cube_facets
          (polymesh, ipoly, cube, k, num_proper_facets, facet_vert_list);
      }
      else {
        error.AddMessage
          ("Programming error.  Unable to determine polytope facets.");
        throw error;
      }
    }
  }

  if (k != numv_per_facet * num_proper_facets) {
    error.AddMessage
      ("Programming error.  Incorrect number of vertices in facet_vert_list.");
    throw error;
  }

  if (mesh_dimension > 2) {
    // Sort vertices in each facet
    // If mesh_dimension == 2, edge endpoints are already sorted.
    for (int jf = 0; jf < num_proper_facets; jf++) {
      sort(facet_vert_list.begin()+jf*numv_per_facet,
           facet_vert_list.begin()+(jf+1)*numv_per_facet);
    }
  }

  vector<int> index_sorted(num_proper_facets);

  for (int i = 0; i < num_proper_facets; i++) 
    { index_sorted[i] = i; }

  TUPLE_LESS_THAN<int,int> facet_less_than
    (numv_per_facet, &(facet_vert_list.front()));
  sort(index_sorted.begin(), index_sorted.end(), facet_less_than);

  is_facet_nonmanifold.assign(num_proper_facets, false);
  contains_nonmanifold_facet.assign(num_poly, false);
  is_facet_boundary.assign(num_proper_facets, false);
  contains_boundary_facet.assign(num_poly, false);

  int j = 0;
  while (j+1 < num_proper_facets) {
    int k = j+1;
    int jf = index_sorted[j];
    int kf = index_sorted[k];
    while (k < num_proper_facets && 
           facet_equals(&(facet_vert_list.front()), numv_per_facet, jf, kf)) {
      k++;
      kf = index_sorted[k];
    };

    int num_duplicate = k - j;

    if (num_duplicate > 2) {

      // store non-manifold facet
      for (int iv = jf*numv_per_facet; iv < (jf+1)*numv_per_facet; iv++)
        { nonmanifold_facet_vert.push_back(facet_vert_list[iv]); };

      // set facets to non-manifold
      for (int m = j; m <= k; m++) {
        int mf = index_sorted[m];
        is_facet_nonmanifold[mf] = true;
      };
    }
    else if (num_duplicate == 1) {

      // store boundary facet
      for (int iv = jf*numv_per_facet; iv < (jf+1)*numv_per_facet; iv++) {
        int iv2 = facet_vert_list[iv];
        boundary_facet_vert.push_back(iv2);
      }

      int num_boundary_facets = boundary_facet_vert.size()/numv_per_facet;
      bool flag_internal =
        is_internal(bounding_box, boundary_facet_vert, 
                    numv_per_facet, num_boundary_facets-1);
      internal_boundary_facet.push_back(flag_internal);
      if (flag_internal) { num_internal_boundary_facets++; }

      if (flag_small_bounding_box) {
        flag_internal = false;
      }
      else {
        flag_internal =
          is_internal(contracted_bounding_box, boundary_facet_vert,
                      numv_per_facet, num_boundary_facets-1);
      }
      far_from_bounding_box.push_back(flag_internal);
      if (flag_internal) { num_deep_boundary_facets++; }

      // set facets to boundary
      for (int m = j; m < k; m++) {
        int mf = index_sorted[m];
        is_facet_boundary[mf] = true;
      };
    }

    j = k;
  }

  // set contains_nonmanifold_facet[] and contains_boundary_facet[]
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {
    if (!is_degenerate[ipoly] && !is_duplicate[ipoly]) {
      int jf = first_poly_facet[ipoly];

      if (mesh_dimension == 2) {
        for (int k = 0; k < polymesh.NumPolyVert(ipoly); k++) {
          if (is_facet_nonmanifold[jf+k])
            contains_nonmanifold_facet[ipoly] = true;
          if (is_facet_boundary[jf+k])
            contains_boundary_facet[ipoly] = true;
        }
      }
      else {
        int num_facets_per_poly; 

        if (flag_cube_file) {
          num_facets_per_poly = NUM_FACETS_PER_CUBE;
        }
        else {
          num_facets_per_poly = num_facets_per_simplex;
        }

        // flag_simplex_input is true.
        for (int k = 0; k < num_facets_per_poly; k++) {
          if (is_facet_nonmanifold[jf+k])
            contains_nonmanifold_facet[ipoly] = true;
          if (is_facet_boundary[jf+k])
            contains_boundary_facet[ipoly] = true;
        }
      }
    }
  }

  set_in_nonmanifold_facet();
}


// Return number of deep vertices in list
int count_deep_vertices(const std::vector<int> & vlist)
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

bool are_simplices_adjacent(int ks, int js)
// return true if simplices ks and js share facet
{
  const int numv_per_simplex = num_vert_per_poly;
  const int numv_per_facet = mesh_dimension;
  int jfacet_vert[numv_per_simplex];
  int kfacet_vert[numv_per_simplex];

  for (int i = 0; i < numv_per_simplex; i++) {
    jfacet_vert[i] = simplex_vert[js*numv_per_simplex+i];
    kfacet_vert[i] = simplex_vert[ks*numv_per_simplex+i];
  }
			
  sort(jfacet_vert, jfacet_vert+numv_per_simplex);
  sort(kfacet_vert, kfacet_vert+numv_per_simplex);
  
  int j = 0;
  int k = 0;
  int num_match = 0;
  while (j < numv_per_simplex && k < numv_per_simplex) {
    if (jfacet_vert[j] == kfacet_vert[k]) {
      j++;
      k++;
      num_match++;
    }
    else if (jfacet_vert[j] > kfacet_vert[k]) {
      k++;
    }
    else {
      j++;
    }
  }

  if (num_match >= numv_per_facet)
    return(true);
  else
    return(false);
}

// return true if 2D polygons kpoly and jpoly are adjacent
bool are_poly2D_adjacent(int kpoly, int jpoly)
{
  for (int k0 = 0; k0 < polymesh.NumPolyVert(kpoly); k0++) {
    int kv0 = polymesh.Vertex(kpoly, k0);
    int k1 = (k0+1)%polymesh.NumPolyVert(kpoly);
    int kv1 = polymesh.Vertex(kpoly, k1);

    for (int j0 = 0; j0 < polymesh.NumPolyVert(jpoly); j0++) {
      int jv0 = polymesh.Vertex(jpoly, j0);
      int j1 = (j0+1)%polymesh.NumPolyVert(jpoly);
      int jv1 = polymesh.Vertex(jpoly, j1);
      if (kv0 == jv0 && kv1 == jv1) { return(true); }
      if (kv0 == jv1 && kv1 == jv0) { return(true); }
    }
  }

  return(false);
}


/// Return true if hexahedra ihexA and ihexB share facet.
bool are_hexahedra_adjacent(int ihexA, int ihexB)
{
  const int DIM3(3);
  const int NUMV_PER_CUBE = 8;
  const int NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;
  const int NUM_CUBE_FACETS = 6;
  int facetA_vert[NUMV_PER_CUBE_FACET];
  int facetB_vert[NUMV_PER_CUBE_FACET];
  static CUBE_TYPE cube(DIM3);

  for (int jfA = 0; jfA < NUM_CUBE_FACETS; jfA++) {
    
    for (int i = 0; i < NUMV_PER_CUBE_FACET; i++) {
      int iA = cube.FacetVertex(jfA, i);
      facetA_vert[i] = polymesh.Vertex(ihexA, iA);
    }
    sort(facetA_vert, facetA_vert+NUMV_PER_CUBE_FACET);

    for (int jfB = 0; jfB < NUM_CUBE_FACETS; jfB++) {

      for (int i = 0; i < NUMV_PER_CUBE_FACET; i++) {
        int iB = cube.FacetVertex(jfB, i);
        facetB_vert[i] = polymesh.Vertex(ihexB, iB);
      }
      sort(facetB_vert, facetB_vert+NUMV_PER_CUBE_FACET);

      if (equal(facetA_vert, facetA_vert+NUMV_PER_CUBE_FACET,
                facetB_vert)) 
        { return(true);  }
    }
  }

  return(false);
}


bool are_poly_adjacent(int kpoly, int jpoly)
// return true if polytopes kpoly and jpoly share facet
{
  const int NUMV_PER_CUBE = 8;
  IJK::PROCEDURE_ERROR error("are_poly_adjacent");

  if (mesh_dimension == 2) {
    return(are_poly2D_adjacent(kpoly, jpoly));
  }
  else if (flag_simplex_file) {
    return(are_simplices_adjacent(kpoly, jpoly));
  }
  else if (num_vert_per_poly == NUMV_PER_CUBE) {
    return(are_hexahedra_adjacent(kpoly, jpoly));
  }
  else {
    error.AddMessage
      ("Programming error. Polytopes are not simplices and not 2D.");
    error.AddMessage("Cannot determine if polytopes are adjacent.");
    throw error;
  }
}

bool are_edges_adjacent
(const int je, const int ke, 
 const std::vector<VERTEX_INDEX> & boundary_edge_endpoint)
{
  const int j = 2*je;
  const int k = 2*ke;

  if (boundary_edge_endpoint[j] == boundary_edge_endpoint[k])
    { return(true); }
  if (boundary_edge_endpoint[j+1] == boundary_edge_endpoint[k])
    { return(true); }
  if (boundary_edge_endpoint[j] == boundary_edge_endpoint[k+1])
    { return(true); }
  if (boundary_edge_endpoint[j+1] == boundary_edge_endpoint[k+1])
    { return(true); }

  return(false);
}


bool are_edges_connected
(const std::vector<VERTEX_INDEX> & boundary_edge_endpoint)
{
  const int num_edges = boundary_edge_endpoint.size()/2;

  if (num_edges == 0) { return(true); }

  vector<bool> is_reachable(num_edges,false);
  vector<int> reachable;
  reachable.push_back(0);

  while (reachable.size() > 0) {
    int je = reachable.back();
    reachable.pop_back();
    is_reachable[je] = true;

    for (int ke = 0; ke < num_edges; ke++) {
      if (!is_reachable[ke]) {
        if (are_edges_adjacent(je, ke, boundary_edge_endpoint)) 
          { reachable.push_back(ke); }
      }
    }
  }

  for (int ke = 0; ke < num_edges; ke++) {
    if (!is_reachable[ke]) 
      { return(false); }
  }

  return(true);
}

void compute_boundary_edges
(const POLYMESH_TYPE & polymesh, 
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
 const VERTEX_INDEX iv,
 const CUBE_TYPE & cube,
 std::vector<VERTEX_INDEX> & boundary_edge_endpoint)
{
  POLYMESH_TYPE link_mesh;

  boundary_edge_endpoint.clear();

  compute_vertex_link_in_cube_mesh
    (polymesh, vertex_info, iv, cube, link_mesh);

  // Convert quad vertices to cylic order around quadrilaterals.
  reorder_quad_vertices(link_mesh.element);

  const int num_edges = sum_num_poly_vert(link_mesh);
  std::vector<int> edge_endpoint(2*num_edges);

  int k = 0;
  int num_edges2;
  for (int ipoly = 0; ipoly < link_mesh.NumPoly(); ipoly++) {
    store_polygon_edges(link_mesh, ipoly, k, num_edges2, edge_endpoint);
  }

  std::vector<int> index_sorted(num_edges);

  for (int i = 0; i < num_edges; i++) 
    { index_sorted[i] = i; }

  TUPLE_LESS_THAN<int,int> edge_less_than(2, &(edge_endpoint.front()));
  sort(index_sorted.begin(), index_sorted.end(), edge_less_than);

  int j = 0;
  while (j+1 < num_edges) {
    int k = j+1;
    int jf = index_sorted[j];
    int kf = index_sorted[k];
    while (k < num_edges && 
           facet_equals(&(edge_endpoint.front()), 2, jf, kf)) {
      k++;
      kf = index_sorted[k];
    };

    int num_duplicate = k - j;

    if (num_duplicate == 1) {

      // store boundary edge
      for (int iv = jf*2; iv < (jf+1)*2; iv++) {
        int iv2 = edge_endpoint[iv];
        boundary_edge_endpoint.push_back(iv2);
      }
    }

    j = k;
  }

}


/// Identify non-manifold vertices.
int identify_nonmanifold_vertices()
{
  const int NUMV_PER_CUBE = 8;
  IJK::PROCEDURE_ERROR error("identify_nonmanifold_vertices");

  if (mesh_dimension != 2 && (!flag_simplex_file) &&
      (num_vert_per_poly != NUMV_PER_CUBE)) {
      ("Programming error.  Mesh is not dimension 2 and not a mesh of simplices.");
    error.AddMessage("  Unable to determine adjacent polytopes.");
    throw error;
  }

  if (in_nonmanifold_facet.size() != num_vertices) {
    error.AddMessage
      ("Programming error.  Array in_nonmanifold_facet[] not set.");
    throw error;
  }

  if (flag_cube_file) {
    if (in_nonmanifold_edge.size() != num_vertices) {
      error.AddMessage
        ("Programming error.  Array in_nonmanifold_edge[] not set.");
      throw error;
    }
  }

  nonmanifold_vert_list.clear();
  nonmanifold_vert.assign(num_vertices, false);

  VERTEX_POLY_INCIDENCE_TYPE vertex_info(polymesh);

  for (VERTEX_INDEX iv = 0; iv < vertex_info.NumVertices(); iv++) {
    const int num_incident = vertex_info.NumIncidentPoly(iv);
    if (!in_nonmanifold_facet[iv] && num_incident > 0) {

      if (flag_cube_file) {
        if (in_nonmanifold_edge[iv]) { continue; }
      }

      vector<bool> is_reachable(num_incident, false);
      vector<int> reachable;
      reachable.push_back(0);

      while (reachable.size() > 0) {
        int n = reachable.size()-1;
        int j = reachable[n];
        reachable.pop_back();
        is_reachable[j] = true;
        int jpoly = vertex_info.IncidentPoly(iv, j);

        for (int k = 0; k < num_incident; k++) {
          if (!is_reachable[k]) {
            int kpoly = vertex_info.IncidentPoly(iv, k);
            if (are_poly_adjacent(kpoly, jpoly)) {
              reachable.push_back(k);
            }
          }
        }
      }

      for (int k = 0; k < num_incident; k++) {
        if (!is_reachable[k]) {
          nonmanifold_vert_list.push_back(iv);
          nonmanifold_vert[iv] = true;
          break;
        }

      }
    }
  }


  if (flag_cube_file && mesh_dimension == DIM3) {

    CUBE_TYPE cube(mesh_dimension);

    // Check that the boundary of each link is connected.
    for (VERTEX_INDEX iv = 0; iv < vertex_info.NumVertices(); iv++) {

      if (nonmanifold_vert[iv] || in_nonmanifold_facet[iv] ||
          in_nonmanifold_edge[iv])
        { continue; }

      vector<int> boundary_edge_endpoint;

      compute_boundary_edges
        (polymesh, vertex_info, iv, cube, boundary_edge_endpoint);

      if (!are_edges_connected(boundary_edge_endpoint)) {
        nonmanifold_vert_list.push_back(iv);
        nonmanifold_vert[iv] = true;
      }
    }

  }

  return(nonmanifold_vert_list.size());
}


/// Identify non-manifold vertices in mesh of 2D polygons.
void identify_nonmanifold_vertices_in_dim2_mesh()
{
  nonmanifold_vert_list.clear();

  // initialize in_nonmanifold_facet[iv] to false for all vertices iv
  in_nonmanifold_facet.assign(num_vertices, false);

  for (int i = 0; i < nonmanifold_facet_vert.size(); i++) {
    int iv = nonmanifold_facet_vert[i];
    in_nonmanifold_facet[iv] = true;
  }

  VERTEX_POLY_INCIDENCE<int,int> vertex_info(polymesh);

  for (int iv = 0; iv < vertex_info.NumVertices(); iv++) {
    const int num_incident = vertex_info.NumIncidentPoly(iv);
    if (!in_nonmanifold_facet[iv] && num_incident > 0) {
      
      vector<bool> is_reachable(num_incident, false);
      vector<int> reachable;
      reachable.push_back(0);

      while (reachable.size() > 0) {
        int n = reachable.size()-1;
        int j = reachable[n];
        reachable.pop_back();
        is_reachable[j] = true;
        int jpoly = vertex_info.IncidentPoly(iv, j);

        for (int k = 0; k < num_incident; k++) {
          if (!is_reachable[k]) {
            int kpoly = vertex_info.IncidentPoly(iv, k);
            if (are_poly_adjacent(kpoly, jpoly)) {
              reachable.push_back(k);
            }
          }
        }
      }

      for (int k = 0; k < num_incident; k++) {
        if (!is_reachable[k]) {
          nonmanifold_vert_list.push_back(iv);
          break;
        }

      }
    }
  }

}


/// Identify non-manifold edges
int identify_nonmanifold_edges()
{
  const int DIM3(3);
  const int NUMV_PER_CUBE = 8;
  const CUBE_TYPE cube(DIM3);
  std::vector<int> poly_list;
  int num_nonmanifold_edges = 0;
  IJK::PROCEDURE_ERROR error("identify_nonmanifold_edges");

  if (num_vert_per_poly != NUMV_PER_CUBE) {
      ("Programming error.  Mesh is not a mesh of hexahedra.");
    error.AddMessage("  Non-manifold edges only implemented for mesh of hexahedra.");
    throw error;
  }

  nonmanifold_edge_vert.clear();

  VERTEX_POLY_INCIDENCE<int,int> vertex_info(polymesh);
  VERTEX_ADJACENCY_LIST<int,int> adjacency_list;

  if (in_nonmanifold_facet.size() != num_vertices) {
    error.AddMessage
      ("Programming error.  Array in_nonmanifold_facet[] not set.");
    throw error;
  }

  adjacency_list.SetFromMeshOfCubes(polymesh, cube);

  for (int iv0 = 0; iv0 < adjacency_list.NumVertices(); iv0++) {

    if (in_nonmanifold_facet[iv0]) { continue; }

    for (int i1 = 0; i1 < adjacency_list.NumAdjacent(iv0); i1++) {
      const int iv1 = adjacency_list.AdjacentVertex(iv0, i1);

      if (iv0 >= iv1) {
        // Skip if iv0 == iv1.
        // Skip if iv0 > iv1, since edge handled as neighbor of iv1.
        continue;
      }

      if (in_nonmanifold_facet[iv1]) { continue; }

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
            if (are_poly_adjacent(kpoly, jpoly)) {
              reachable.push_back(k);
            }
          }
        }
      }

      for (int k = 0; k < num_incident_poly; k++) {
        if (!is_reachable[k]) {
          nonmanifold_edge_vert.push_back(iv0);
          nonmanifold_edge_vert.push_back(iv1);
          num_nonmanifold_edges++;
          break;
        }
      }
    }
  }

  set_in_nonmanifold_edge();

  return(num_nonmanifold_edges);
}

// Return true if facet jf is internal
bool is_internal(const BOUNDING_BOX & bounding_box, 
                 const vector<int> & facet_vlist, 
                 const int numv_per_facet,
                 const int jf)
{
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
// PLOT ANGLE ROUTINES
// **************************************************

void compute_polygon_angles(ANGLE_TABLE & angle_table)
{
  ANGLE_TYPE min_angle, max_angle;
  int num_angle;

  angle_table.min_polygon_angle_freq.Include();
  angle_table.min_polygon_angle_freq.SetAll(0);
  angle_table.max_polygon_angle_freq.Include();
  angle_table.max_polygon_angle_freq.SetAll(0);

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    if (is_degenerate[ipoly]) { continue; }

    compute_min_max_polygon_angles(ipoly, min_angle, max_angle, num_angle);

    if (num_angle < 3) { continue; }

    int imin = floor(min_angle+0.5);
    int imax = floor(max_angle+0.5);

    angle_table.min_polygon_angle_freq.Increment(imin);
    angle_table.max_polygon_angle_freq.Increment(imax);
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
  num_nonmanifold_facets = 0;
  num_nonmanifold_edges = 0;
  num_nonmanifold_vertices = 0;
  num_deep_nonmanifold_vertices = 0;
}

// Return true if all numbers are zero.
bool MESH_INFO::AreAllZero() const
{
  if (num_poly_with_duplicate_vertices > 0) { return(false); }
  if (num_duplicate_poly > 0) { return(false); }
  if (num_nonmanifold_facets > 0) { return(false); }
  if (num_nonmanifold_edges > 0) { return(false); }
  if (num_nonmanifold_vertices > 0) { return(false); }
  if (num_deep_nonmanifold_vertices > 0) { return(false); }

  return(true);
}


// **************************************************
// Class SIMPLEX_INFO
// **************************************************

SIMPLEX_INFO::SIMPLEX_INFO
(const int numv_per_simplex, const int num_simplices,
 const int * simplex_vert)
{
  this->numv_per_simplex = numv_per_simplex;
  this->num_simplices = num_simplices;

  num_adjacent = new int[num_simplices];
  adjacent = new int[num_simplices*numv_per_simplex];

}

SIMPLEX_INFO::~SIMPLEX_INFO()
{
  delete [] num_adjacent;
  num_adjacent = NULL;
  delete [] adjacent;
  adjacent = NULL;
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
        mesh_dimension = get_arg_int(iarg, argc, argv, error);
        is_mesh_dimension_set = true;
        iarg++;
        break;

      case VERTEX_PARAM:
        vertex_index = get_arg_int(iarg, argc, argv, error);
        iarg++;
        vertex_info_flag = true;
        general_info_flag = false;
        break;

      case SIMPLEX_PARAM:
        simplex_index = get_arg_int(iarg, argc, argv, error);
        iarg++;
        simplex_info_flag = true;
        general_info_flag = false;
        break;

      case POLY_PARAM:
        poly_index = get_arg_int(iarg, argc, argv, error);
        iarg++;
        poly_info_flag = true;
        general_info_flag = false;
        break;

      case VLIST_PARAM:
        vlist_flag = true;
        general_info_flag = false;
        break;

      case PLIST_PARAM:
        plist_flag = true;
        general_info_flag = false;
        break;

      case CONTAINSV_PARAM:
        contains_vertex_index = get_arg_int(iarg, argc, argv, error);
        iarg++;
        contains_vertex_flag = true;
        general_info_flag = false;
        break;

      case CONTAINSE_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%d", &edge_end0_index);
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%d", &edge_end1_index);
        contains_edge_flag = true;
        general_info_flag = false;
        break;

      case MANIFOLD_PARAM:
        manifold_flag = true;
        general_info_flag = false;
        break;

      case SELFI_PARAM:
        flag_report_self_intersections = true;
        flag_use_grid_of_bins = true;
        general_info_flag = false;
        break;

      case SELFI_NO_GRID_PARAM:
        flag_report_self_intersections = true;
        flag_use_grid_of_bins = false;
        general_info_flag = false;
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
        sscanf(argv[iarg], "%f", &angle_le);
        is_min_angle_set = true;
        break;

      case ANGLE_GE_PARAM:
        iarg++;
        if (iarg >= argc) usage_error();
        sscanf(argv[iarg], "%f", &angle_ge);
        is_max_angle_set = true;
        break;

      case LIST_DUP_PARAM:
        flag_list_duplicate_vertices = true;
        flag_list_duplicate_poly = true;
        break;

      case INTERNAL_PARAM:
        flag_internal = true;
        break;

      case REPORT_DEEP_PARAM:
        flag_report_deep = true;
        break;

      case OUT_VALUES_PARAM:
        flag_output_only_values = true;
        break;

      case OUT_MIN_ANGLE_PARAM:
        flag_output_min_angle = true;
        general_info_flag = false;
        break;

      case OUT_MAX_ANGLE_PARAM:
        flag_output_max_angle = true;
        general_info_flag = false;
        break;

      case PLOT_ANGLES_PARAM:
        flag_plot_angles = true;
        flag_plot_min_polygon_angles = true;
        flag_plot_max_polygon_angles = true;
        break;

      case FOR_EACH_TYPE_PARAM:
        flag_for_each_type = true;
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

  if (is_min_angle_set) {
    if (angle_le < 0 || angle_le > 180) {
      cerr << "Illegal value " << angle_le << " for option -angle_le <A>."
           << endl;
      cerr << "  Angle <A> must be in range [0,180].";
      exit(20);
    }
  }

  if (is_max_angle_set) {
    if (angle_ge < 0 || angle_ge > 180) {
      cerr << "Illegal value " << angle_ge << " for option -angle_ge <A>."
           << endl;
      cerr << "  Angle <A> must be in range [0,180].";
      exit(20);
    }
  }
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

void usage_msg()
{
  cerr << "Usage: ijkmeshinfo [OPTIONS] {input file}" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  [-mesh_dim {mdim}] [-vertex {vnum}] [-simplex {snum}] [-poly {pnum}]"
       << endl;
  cerr << "  [-vlist] [-plist] [-manifold]" << endl;
  cerr << "  [-containsv {vnum}] [-containse {end0} {end1}]" << endl;
  cerr << "  [-minc \"min coord\"] [-maxc \"max coord\"]" << endl;
  cerr << "  [-min_numv <N>] [-max_numv <N>]" << endl;
  cerr << "  [-angle_le <A>] [-angle_ge <A>]" << endl;
  cerr << "  [-list_dup] [-internal]" << endl;
  cerr << "  [-selfI]" << endl;
  cerr << "  [-out_values] [-out_min_angle] [-out_max_angle] [-plot_angles]"
       << endl;
  cerr << "  [-report_deep] [-for_each_type]" << endl;
  cerr << "  [-terse] [-help]" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help_msg()
{
  cerr << "Usage: ijkmeshinfo [OPTIONS] {input file}" << endl;
  cerr << "  [-mesh_dim {mdim}]:   Mesh dimension." << endl;
  cerr << "  [-vertex {vnum}]:     Print coordinates of vertex {vnum}."
       << endl;
  cerr << "  [-simplex {snum}]:    Print vertices of simplex {snum}."
       << endl;
  cerr << "  [-poly {pnum}]:       Print vertices of poly {pnum}."
       << endl;
  cerr << "  [-vlist]:             Print list of vertex coordinates." << endl;
  cerr << "     Prints only vertices with coordinates between min coord and max coord."
       << endl;
  cerr << "  [-plist]:             Print list of polygons/polytopes." << endl;
  cerr << "     Prints only polygons/polytopes with number of vertices between"
       << endl;
  cerr << "       <min_numv> and <max_numv>." << endl;
  cerr << "  [-manifold]:          Print manifold information." << endl;
  cerr << "  [-selfI]:             Print self intersections." << endl;
  cerr << "  [-containsv {vnum}]:  Print list of simplices containing vertex {vnum}."
       << endl;
  cerr << "  [-containse {end0} {end1}}]:  Print list of simplices containing"
       << endl;
  cerr << "                          edge ({end0},{end1})." << endl;
  cerr << "  [-minc \"min coord\"]:  Minimum coordinates." << endl;
  cerr << "  [-maxc \"min coord\"]:  Maximum coordinates." << endl;
  cerr << "  [-min_numv <N>]:      Minimum number of poly vertices." << endl;
  cerr << "  [-max_numv <N>]:      Maximum number of poly vertices." << endl;
  cerr << "  [-angle_le <A>]:" << endl;
  cerr << "     Report number of polygons with angle less than or equal to <A>." << endl;
  cerr << "     Use with -plist to list polygons with angle less than or equal to <A>." << endl;
  cerr << "  [-angle_ge <A>]:" << endl;
  cerr << "     Report number of polygons with angle greater than or equal to <A>." << endl;
  cerr << "     Use with -plist to list polygons with angle greater than or equal to <A>." << endl;
  cerr << "  [-list_dup]:          List duplicate poly vertices or poly vertices" << endl;
  cerr << "                           with identical coordinates." << endl;
  cerr << "  [-internal]:          Report only angles for internal polygons."
       << endl;
  cerr << "  [-out_values]:        Output only values for -manifold option."
       << endl;
  cerr << "     Output number of non-manifold vertices, number of non-manifold vertices"
       << endl
       << "     at least 1 from bounding box, number of non-manifold edges," 
       << endl
       << "     number of internal boundary facets and number of facets"
       << endl
       << "     at least 1 from bounding box." << endl;
  cerr << "  [-out_min_angle]:     Output minimum triangle angle." << endl;
  cerr << "  [-out_max_angle]:     Output maximum triangle angle." << endl;
  cerr << "  [-plot_angles]:       Create gnuplot (.gplt) files of min and max"
       << endl
       << "                          triangle angles." << endl;
  cerr << "  [-report_deep]:       Report only boundary facets at least distance 1"
       << endl
       << "                          from bounding box." << endl;
  cerr << "  [-for_each_type]:     Reports angles separately for triangles, "
       << endl
       << "                          quadrilaterals, pentagons, ..." << endl;
  cerr << "  [-terse]:             Terse output.  Affects -manifold and -containsv options." << endl;
  cerr << "  [-help]:              Print this help message." << endl;
  exit(20);
}

