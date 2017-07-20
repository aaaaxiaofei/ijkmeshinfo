/// \file ijkmeshinfo.h
/// compute mesh information

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2010-2016 Rephael Wenger

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

#ifndef _IJKMESHINFO_
#define _IJKMESHINFO_

#include <iostream>
#include <string>

#include "ijk.txx"
#include "ijkdatatable.txx"


namespace IJKMESHINFO {

// **************************************************
// Global constants
// **************************************************

  const int DIM2(2);
  const int DIM3(3);

// **************************************************
// Types
// **************************************************

  typedef double COORD_TYPE;    // Use double to detect self-intersections.
  typedef COORD_TYPE * COORD_TYPE_PTR;
  typedef float ANGLE_TYPE;
  typedef IJK::BOX<COORD_TYPE> BOUNDING_BOX;
  typedef int NUM_TYPE;


// **************************************************
// Class MESH_INFO
// **************************************************

// mesh information
  class MESH_INFO {

  protected:
    void Init();

  public:
    MESH_INFO(){ Init(); };

    int num_poly_with_duplicate_vertices;
    int num_duplicate_poly;
    int num_nonmanifold_facets;
    int num_nonmanifold_edges;
    int num_nonmanifold_vertices;
    int num_deep_nonmanifold_vertices;
    int num_poly_with_orientation_conflicts;

    /// Return true if all non-manifold numbers are zero.
    bool AreAllNonManifoldZero() const;

    /// Return true if all numbers are zero.
    bool AreAllZero() const;
  };


// **************************************************
// Class SIMPLEX_INFO
// **************************************************

// list of simplices adjacent to a given simplex
  class SIMPLEX_INFO {
  protected:
    int numv_per_simplex;
    int num_simplices;
    int * num_adjacent;
    int * adjacent;

  public:
    SIMPLEX_INFO(const int numv_per_simplex, const int nums,
		 const int * simplex_vert);
    ~SIMPLEX_INFO();

    // get functions
    int NumAdjacent(const int is) const
    { return(num_adjacent[is]); };

    int AdjacentSimplex(const int iv, const int k)
    { return(adjacent[numv_per_simplex*iv + k]); };
  };


// **************************************************
// Class FACET_INFO
// **************************************************

  class FACET_INFO {
  public:
    int poly_containing_facet;
    int facet_index;

  public:
    FACET_INFO() {};
    FACET_INFO(const int ipoly, const int jf)
    { poly_containing_facet = ipoly, facet_index = jf; }
  };

  typedef std::vector<FACET_INFO> FACET_INFO_ARRAY;


// **************************************************
// Class GRID_OF_BINS_3D
// *************************************************

  class GRID_OF_BINS_3D {

  protected:
    std::vector<int> ** bin;
    int num_bins_along_axis[DIM3];
    int num_bins;
    COORD_TYPE minC[DIM3];
    COORD_TYPE maxC[DIM3];
    bool bins_are_empty;

    void Init(const int k);

  public:
    GRID_OF_BINS_3D(const int k) { Init(k); };
    ~GRID_OF_BINS_3D();

    int Dimension() const { return(DIM3); };
    int NumBinsAlongAxis(const int d) const 
    { return(num_bins_along_axis[d]); };
    int NumBins() const { return(num_bins); };
    const COORD_TYPE * MinCoord() const
    { return(minC); }
    const COORD_TYPE * MaxCoord() const
    { return(maxC); }

    int LocateBinCoord(const int d, const COORD_TYPE c) const;
    int ComputeBinIndex(const int binCoord[DIM3]) const;

    /// Compute min/max bin coordinates containing triangle (v0,v2,v2).
    void ComputeMinMaxBinCoord
      (const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
       const COORD_TYPE v2[DIM3],
       int min_grid_coord[DIM3], int max_grid_coord[DIM3]) const;

    const std::vector<int> * Bin(const int i)
    { return(bin[i]); };

    void SetMinCoord(const COORD_TYPE minC[DIM3]);
    void SetMaxCoord(const COORD_TYPE maxC[DIM3]);

    void InsertTri(const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
                   const COORD_TYPE v2[DIM3], const int triangle_index);

  };


  // **************************************************
  // Class ANGLE_TABLE
  // **************************************************

  template <class NTYPE, class DTYPE>
    class DATA_COLUMN: public IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE> {

    protected:
      bool is_included;
      bool is_hidden;
      DTYPE sum;

      void Init() { is_included = false; is_hidden = false; sum = 0; };

    public:
    DATA_COLUMN<NTYPE,DTYPE>() :
      IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>() { Init(); };
      DATA_COLUMN<NTYPE,DTYPE>
      (const std::string & label, const NTYPE num_rows) :
        IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>(label, num_rows) { Init(); };

    // set functions
    void Include();
    void Hide() { is_hidden = true; };
    void Show() { is_hidden = false; };
    void ComputeSum();

    // get functions
    bool IsIncluded() const { return(is_included); };
    bool IsHidden() const { return(is_hidden); };
    DTYPE Sum() const { return(sum); };

    // write functions
    void WriteLabel(std::ostream & out, const std::string & separator) const;
    void WriteData(std::ostream & out, const std::string & separator, 
		   const NTYPE width, const NTYPE irow) const;
    void WriteNormalizedData
      (std::ostream & out, const std::string & separator, const NTYPE width,
       const double normalization_factor, const NTYPE irow) const;
  };

  class ANGLE_TABLE:public IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE> {

  public:
    static const NUM_TYPE NUM_ROWS = 180;

    DATA_COLUMN<NUM_TYPE, NUM_TYPE> angle;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> min_polygon_angle_freq;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> max_polygon_angle_freq;

  public:
    ANGLE_TABLE():
      IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE>(NUM_ROWS),
      angle("angle", NUM_ROWS), 
      min_polygon_angle_freq("min-poly-angle", NUM_ROWS),
      max_polygon_angle_freq("max-poly-angle", NUM_ROWS)
        {};

    // set routines
    void HideAllExceptAngleColumn();
    void ComputeSum();

    // write routines
    void WriteColumnLabels
      (std::ostream & out, const std::string & separator) const;
    void WriteColumnData
      (std::ostream & out, const std::string & separator, 
       const NUM_TYPE width) const;
    void WriteNormalizedColumnData
      (std::ostream & out, const std::string & separator, const NUM_TYPE width,
       const double normalization_factor) const;

  };

// **************************************************
// DATA_COLUMN member functions
// **************************************************

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::Include()
    {
      is_included = true;

      if (this->NumRows() > 0 && this->data == NULL) 
        { this->data = new DTYPE[this->NumRows()]; }
    }

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::ComputeSum()
    {
      if (IsIncluded()) 
        { sum = IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>::Sum(); }
      else
        { sum = 0; }
    }

  template <class NTYPE, class DTYPE>
  void DATA_COLUMN<NTYPE,DTYPE>::WriteLabel
  (std::ostream & out, const std::string & separator) const
  {
    if (IsIncluded() && !IsHidden()) 
      { out << separator << this->Label(); }
  }

  template <class NTYPE, class DTYPE>
  void DATA_COLUMN<NTYPE,DTYPE>::WriteData
  (std::ostream & out, const std::string & separator, 
   const NTYPE width, const NTYPE irow) const
  {
    if (IsIncluded() && !IsHidden()) {
      out << separator;
      out.width(width);
      out << this->data[irow]; 
    }
  }

  template <class NTYPE, class DTYPE>
  void DATA_COLUMN<NTYPE,DTYPE>::WriteNormalizedData
  (std::ostream & out, const std::string & separator, const NTYPE width,
   const double normalization_factor, const NTYPE irow) const
  {
    if (IsIncluded() && !IsHidden()) {
      double v = normalization_factor*double(this->data[irow])/Sum();
      out << separator;
      out.width(width);
      out << v; 
    }
  }

};

#endif
