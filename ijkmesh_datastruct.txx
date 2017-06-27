/// \file ijkmesh_datastruct.txx
/// ijk template classes for polyhedral mesh data structures.
/// Version 0.2.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2010-2017 Rephael Wenger

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

#ifndef _IJKMESH_DATASTRUCT_
#define _IJKMESH_DATASTRUCT_

#include "ijk.txx"
#include "ijklist.txx"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>


namespace IJK {

  // **************************************************
  // Class POLYMESH 
  // **************************************************

  /// Mesh of polytopes
  template <typename VTYPE, typename NTYPE>
  class POLYMESH:public LIST_OF_LISTS<VTYPE,NTYPE> {

  public:

    /// Vertex index type.
    typedef VTYPE VERTEX_INDEX_TYPE;

    /// NUMBER type.
    typedef NTYPE NUMBER_TYPE;

  public:
    /// constructor
    POLYMESH(){}

    NTYPE NumPoly() const               ///< Number of polytopes.
    { return(this->NumLists()); }

    /// Number of vertices of polytope i.
    NTYPE NumPolyVert(const NTYPE ipoly) const
    { return(this->ListLength(ipoly)); }

    /// Return j'th vertex of polytope ipoly.
    VTYPE Vertex(const NTYPE ipoly, const NTYPE j) const
    { return(this->Element(ipoly,j)); }

    /// Return pointer to first vertex of polytope ipoly.
    const VTYPE * VertexList(const NTYPE ipoly) const
    { return(this->List(ipoly)); }

    /// Return true if polytope ipoly contains vertex iv.
    template <typename VTYPE2>
    bool DoesPolyContainVertex(const NTYPE ipoly, const VTYPE2 iv) const
    { return(this->DoesListContain(ipoly, iv)); }

    /// Return true if 2D polygon ipoly contains edge (iv0,iv1).
    /// @pre Polygon vertices are listed in cylic order around the polygon.
    template <typename VTYPE0, typename VTYPE1>
    bool DoesPoly2DContainEdge
    (const NTYPE ipoly, const VTYPE0 iv0, const VTYPE1 iv1) const;

    /// Add a polytope with list_length vertices.
    template <typename VTYPE2, typename NTYPE2>
    void AddPolytope(const VTYPE2 poly_vert_list[],
                     const NTYPE2 list_length)
    { this->AddList(poly_vert_list, list_length); }

    /// Add a polytope.
    /// - C++ STL vector type for array poly_vert_list[].
    /// @param poly_vert_list List of polytope vertices in C++ STL vector.
    template <typename VTYPE2>
    void AddPolytope(const std::vector<VTYPE2> & poly_vert_list)
    { AddPolytope(&poly_vert_list.front(), poly_vert_list.size()); }

    /// Add polytopes from list where each polytope has 
    ///   num_vert_per_poly vertices.
    /// @param num_vert Total number of vertices in poly_vert_list[].
    /// @param num_vert_per_poly Number of vertices in each polytope.
    template <typename VTYPE2, typename NTYPE1, typename NTYPE2>
    void AddPolytopes(const VTYPE2 poly_vert_list[],
                      const NTYPE1 num_vert,
                      const NTYPE2 num_vert_per_poly)
    { this->AddLists(poly_vert_list, num_vert, num_vert_per_poly); }

    /// Add polytopes from list where each polytope has 
    ///   num_vert_per_poly vertices.
    /// @param num_vert_per_poly Number of vertices in each polytope.
    template <typename VTYPE2>
    void AddPolytopes(const std::vector<VTYPE2> & poly_vert_list, 
                      const NTYPE num_vert_per_poly)
    { this->AddLists(poly_vert_list, num_vert_per_poly); }

    /// Sort the vertices of each polytope.
    /// Note: Sorting may change simplex orientation.
    void SortVert()
    { this->SortEachList(); }

    /// Return list of polytope indices in sorted order.
    /// @param sorted_poly[i] = Index of i'th poly in sorted order.
    /// @pre SortVert() should be called before GetSortedPolytopeIndices.
    template <typename NTYPE2>
    void GetSortedPolytopeIndices(std::vector<NTYPE2> & sorted_poly) const
    { this->GetSortedListIndices(sorted_poly); }

  };


  // **************************************************
  // Class VERTEX_POLY_INCIDENCE_ELEMENT
  // **************************************************

  /// Base element of VERTEX_POLY_INCIDENCE_BASE
  template <typename PTYPE>
  class VERTEX_POLY_INCIDENCE_ELEMENT {

  protected:
    PTYPE poly_index;

  public:
    typedef PTYPE POLY_INDEX_TYPE;

  public:
    PTYPE PolyIndex() const { return(poly_index); };
    void SetPolyIndex(const PTYPE ipoly) { poly_index = ipoly; }
  };


  // **************************************************
  // Class VERTEX_POLY_INCIDENCE_BASE and
  //   VERTEX_POLY_INCIDENCE
  // **************************************************

  /// List of polytopes incident on each vertex. Base class.
  /// @tparam ETYPE List element type.
  ///   Usually derived from VERTEX_POLY_INCIDENCE_ELEMENT
  /// @tparam NTYPE Number type.
  template <typename ETYPE, typename NTYPE>
  class VERTEX_POLY_INCIDENCE_BASE:protected LIST_OF_LISTS<ETYPE,NTYPE> {

  protected:

    /// Number of vertices.
    NTYPE num_vertices;

    /// Number of polytopes.
    NTYPE num_poly;

  protected:

    typedef typename ETYPE::POLY_INDEX_TYPE PTYPE;

    /// Compute number of polytopes incident on each vertex.
    /// Allocate arrays list_length[], first_element[] and element[].
    /// Set number of polytopes and number of vertices.
    /// Set arrays list_length[] and first_element[].
    template <typename VTYPE2, typename NTYPE2>
    void AllocateLists(const POLYMESH<VTYPE2,NTYPE2> & polymesh);


  public:

    // Constructor
    VERTEX_POLY_INCIDENCE_BASE()
    { Clear(); };

    template <typename VTYPE, typename NTYPE2>
    VERTEX_POLY_INCIDENCE_BASE
    (const POLYMESH<VTYPE,NTYPE2> & polymesh)
    { Set(polymesh); }

    NTYPE NumPoly() const
    { return(num_poly); }

    NTYPE NumVertices() const
    { return(num_vertices); }

    template <typename VTYPE>
    NTYPE NumIncidentPoly(const VTYPE iv) const
    { return(this->ListLength(iv)); }

    /// Return list element containing j'th poly incident on vertex iv.
    template <typename VTYPE, typename NTYPE2>
    const ETYPE & Element(const VTYPE iv, const NTYPE2 j) const
    { return(this->element[this->ElementIndex(iv,j)]); }

    /// Return j'th poly containing vertex iv.
    template <typename VTYPE, typename NTYPE2>
    PTYPE IncidentPoly(const VTYPE iv, const NTYPE2 j) const
    { return(Element(iv,j).PolyIndex()); }

    /// Return true if vertex iv is incident on polytope ipoly.
    /// @param[out] iloc Location (index) of ipoly in list
    ///   list of poly incident on iv.
    /// iloc is undefined if vertex iv is not incident on ipoly.
    template <typename VTYPE2, typename PTYPE2, typename NTYPE2>
    bool IsVertexIncidentOnPoly
    (const VTYPE2 iv, const PTYPE2 ipoly, NTYPE2 & iloc) const;

    /// Return true if vertex iv is incident on polytope ipoly.
    /// Version which does not return iloc.
    template <typename VTYPE2, typename PTYPE2>
    bool IsVertexIncidentOnPoly(const VTYPE2 iv, const PTYPE2 ipoly) const
    { NTYPE iloc; 
      return(IsVertexIncidentOnPoly(iv, ipoly, iloc)); }

    /// Get polytopes containing vertices iv0 and iv1.
    /// @param poly_list[] List of polytopes containing vertices iv0 and iv1.
    template <typename VTYPE0, typename VTYPE1, typename PTYPE2>
    void GetPolyContaining
    (const VTYPE0 iv0, const VTYPE1 iv1, std::vector<PTYPE2> & poly_list);

    /// Set vertex poly incidence from polymesh.
    template <typename VTYPE2, typename NTYPE2>
    void Set(const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    /// Clear all lists.
    void Clear();
  };

  /// List of polytopes incident on each vertex.
  /// @tparam PTYPE Polytope index type.
  /// @tparam NTYPE Number type.
  template <typename PTYPE, typename NTYPE>
  class VERTEX_POLY_INCIDENCE:
    public VERTEX_POLY_INCIDENCE_BASE
    <VERTEX_POLY_INCIDENCE_ELEMENT<PTYPE>,NTYPE> {

  public:
    VERTEX_POLY_INCIDENCE() {};

    template <typename VTYPE, typename NTYPE2>
    VERTEX_POLY_INCIDENCE
    (const POLYMESH<VTYPE,NTYPE2> & polymesh)
    { this->Set(polymesh); }
  };


  // **************************************************
  // Class VERTEX_ADJACENCY_LIST_ELEMENT
  // **************************************************


  /// Base element of VERTEX_ADJACENCY_LIST_BASE
  template <typename VTYPE>
  class VERTEX_ADJACENCY_LIST_ELEMENT {

  protected:
    VTYPE vertex;

  public:
    typedef VTYPE VERTEX_INDEX_TYPE;

  public:
    VTYPE Vertex() const { return(vertex); };
    void SetVertex(const VTYPE iv) { vertex = iv; }
  };


  // **************************************************
  // Classes VERTEX_ADJACENCY_LIST_BASE and
  //   VERTEX_ADJACENCY_LIST.
  // **************************************************

  /// List of vertices adjacent to each vertex. Base class.
  /// @tparam ETYPE List element type.
  ///   Usually derived from VERTEX_ADJACENCY_LIST_ELEMENT.
  /// @tparam NTYPE Number type.
  template <typename ETYPE, typename NTYPE>
  class VERTEX_ADJACENCY_LIST_BASE:protected LIST_OF_LISTS<ETYPE,NTYPE> {

  protected:

    /// Number of vertices.
    NTYPE num_vertices;

  protected:
    
    typedef typename ETYPE::VERTEX_INDEX_TYPE VTYPE;

    /// Allocate arrays list_length[], first_element[] and element[].
    /// Set number of vertices.
    /// Set arrays list_length[] and first_element[].
    template <typename NTYPE2>
    void AllocateLists(const std::vector<NTYPE2> & num_adjacent);


  public:

    // Constructor
    VERTEX_ADJACENCY_LIST_BASE() 
    { Clear(); };

    NTYPE NumVertices() const
    { return(num_vertices); }

    template <typename VTYPE2>
    NTYPE NumAdjacent(const VTYPE2 iv) const
    { return(this->ListLength(iv)); }

    /// Return list element containing j'th vertex adjacent to vertex iv.
    template <typename VTYPE, typename NTYPE2>
    const ETYPE & Element(const VTYPE iv, const NTYPE2 j) const
    { return(this->element[this->ElementIndex(iv,j)]); }

    /// Return j'th vertex adjacent to vertex iv.
    template <typename VTYPE2, typename NTYPE2>
    VTYPE AdjacentVertex(const VTYPE2 iv, const NTYPE2 j) const
    { return(Element(iv,j).Vertex()); }

    /// Return true if iv1 is adjacent to iv0.
    template <typename VTYPE0, typename VTYPE1>
    bool IsAdjacent(const VTYPE0 iv0, const VTYPE1 iv1) const;

    template <typename NTYPE2>
    void SetNumVertices(const NTYPE2 num_vertices);

    /// Set adjaceny list from mesh of 2D polygons.
    /// @param polymesh Mesh of 2D polygons.
    ///   - Polygon vertices are listed in clockwise or counter-clockwise
    ///     order around the polygon.
    template <typename VTYPE2, typename NTYPE2>
    void SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    /// Set adjacency list from mesh of (hyper) cubes.
    /// @param polymesh Mesh of (hyper) cubes.
    ///   - Cube vertices are listed in order:
    ///     (0,0,0), (1,0,0), (0,1,0), (1,1,0), 
    ///     (0,0,1), (1,0,1), (0,1,1), (1,1,1).
    template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
    void SetFromMeshOfCubes
    (const POLYMESH<VTYPE2,NTYPE2> & polymesh,
     const CUBE_TYPE & cube);

    /// Clear all lists.
    void Clear();
  };

  /// List of vertices adjacent to each vertex.
  /// @tparam VTYPE Vertex type.
  /// @tparam NTYPE Number type.
  template <typename VTYPE, typename NTYPE>
  class VERTEX_ADJACENCY_LIST:
    public VERTEX_ADJACENCY_LIST_BASE
    <VERTEX_ADJACENCY_LIST_ELEMENT<VTYPE>,NTYPE> {

  public:
    VERTEX_ADJACENCY_LIST() {};
  };


  // **************************************************
  // Class VERTEX_POLY_EDGE_INCIDENCE
  // **************************************************

  /// List of polytopes incident on each vertex.
  /// List of vertices adjacent to each vertex.
  /// For each poly ipoly incident on vertex iv, list of edges
  ///   incident on ipoly and iv.
  /// @tparam VP_ELEMENT_TYPE Element type for vertex poly incidence list.
  ///   Class VP_ELEMENT_TYPE must include member functions SetPolyIndex(),
  ///     first_edge and num_edges and type POLY_INDEX_TYPE.
  /// @tparam VA_ELEMENT_TYPE Element type for vertex adjacent list.
  ///   Class VA_ELEMENT_TYPE must include member functions SetVertex()
  ///     and type VERTEX_INDEX_TYPE.
  template <typename VP_ELEMENT_TYPE, typename VA_ELEMENT_TYPE, typename NTYPE>
  class VERTEX_POLY_EDGE_INCIDENCE {

  protected:
    LIST_OF_LISTS<VP_ELEMENT_TYPE,NTYPE> vertex_poly_incidence;
    LIST_OF_LISTS<VA_ELEMENT_TYPE,NTYPE> vertex_adjacency_list;

    /// For each vertex iv and incident poly ipoly, list of edges
    ///   incident on iv and ipoly.
    /// incident_edge_list.Element(ia,ib)=t is a reference to the t'th element 
    ///   of some list in vertex_adjacency_list.
    /// If (iv0,iv1) is incident on ipoly and 
    ///   vertex_poly_incidence.Element(iv0,j) = ipoly and
    ///   vertex_adjacency_list.Element(iv0,k)=iv1,
    ///   then edge_sublist.Element(IncidentEdgeListIndex(iv0,j),i)=k
    ///   for some i.
    LIST_OF_LISTS<NTYPE,NTYPE> incident_edge_list;

  protected:
    typedef typename VP_ELEMENT_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename VA_ELEMENT_TYPE::VERTEX_INDEX_TYPE VTYPE;

    /// Increment length of list of edges incident on vertex iv and
    ///   poly vertex_poly_incidence(iv,j).
    template <typename VTYPE2, typename NTYPE2>
    void IncrementIncidentEdgeListLength(const VTYPE2 iv, const NTYPE2 j)
    { incident_edge_list.list_length[IncidentEdgeListIndex(iv,j)]++; }

    /// Set k'th element of list of edges incident on vertex iv and
    ///   poly vertex_poly_incidence(iv,j) to x.
    template <typename VTYPE2, typename NTYPE2, typename NTYPE3,
              typename XTYPE>
    void SetIncidentEdgeListElement
    (const VTYPE2 iv, const NTYPE2 j, const NTYPE3 k, const XTYPE x)
    { 
      NTYPE list_index = IncidentEdgeListIndex(iv,j);
      incident_edge_list.ElementRef(list_index,k) = x;
    }


  public:
    VERTEX_POLY_EDGE_INCIDENCE() {};

    NTYPE NumVertices() const
    { return(vertex_poly_incidence.NumLists()); }

    template <typename VTYPE2, typename NTYPE2>
    void SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    template <typename VTYPE>
    NTYPE NumIncidentPoly(const VTYPE iv) const
    { return(vertex_poly_incidence.ListLength(iv)); }

    /// Return j'th poly containing vertex iv.
    template <typename VTYPE, typename NTYPE2>
    PTYPE IncidentPoly(const VTYPE iv, const NTYPE2 j) const
    { return(vertex_poly_incidence.ElementRefConst(iv,j).PolyIndex()); }

    template <typename VTYPE2>
    NTYPE NumAdjacent(const VTYPE2 iv) const
    { return(vertex_adjacency_list.ListLength(iv)); }

    /// Return j'th vertex adjacent to vertex iv.
    template <typename VTYPE2, typename NTYPE2>
    VTYPE AdjacentVertex(const VTYPE2 iv, const NTYPE2 j) const
    { return(vertex_adjacency_list.ElementRefConst(iv,j).Vertex()); }

    /// Return number of edges incident on vertex iv and 
    ///   poly IncidentPoly(iv,j).
    template <typename VTYPE2, typename NTYPE2>
    NTYPE NumIncidentEdges(const VTYPE2 iv, const NTYPE2 j) const
    { return(incident_edge_list.list_length[IncidentEdgeListIndex(iv,j)]); }; 

    /// Return k'th element in list of edges incident on vertex iv and
    ///   poly IncidentPoly(iv,j).
    template <typename VTYPE2, typename NTYPE2, typename NTYPE3>
    NTYPE IncidentEdge(const VTYPE2 iv, const NTYPE2 j, const NTYPE3 k) const
    { return(incident_edge_list.Element(IncidentEdgeListIndex(iv,j),k)); }; 

    /// Return index of list of edges incident on vertex iv and
    ///   poly vertex_poly_incidence(iv,j).
    template <typename VTYPE2, typename NTYPE2>
    NTYPE IncidentEdgeListIndex(const VTYPE2 iv, const NTYPE2 j) const
    { return(vertex_poly_incidence.ElementIndex(iv, j)); }

    /// Return index of k'th edge incident on vertex iv and iv
    ///   poly vertex_poly_incidence(iv,j).
    template <typename VTYPE2, typename NTYPE2, typename NTYPE3>
    NTYPE IncidentEdgeIndex
    (const VTYPE2 iv, const NTYPE2 j, const NTYPE3 k) const
    { return(incident_edge_list.ElementIndex
             (IncidentEdgeListIndex(iv, j), k)); }

    /// Return true if list of edges incident on vertex iv0 and 
    ///  poly vertex_poly_incidence(iv,j) contains edge (iv0, iv1).
    template <typename VTYPE0, typename VTYPE1, typename NTYPE2>
    bool DoesIncidenceListContainEdge
    (const VTYPE0 iv0, const VTYPE1 iv1, const NTYPE2 j) const;

    void Clear();
  };


  // **************************************************
  // POLYMESH compare
  // **************************************************

  /// function class for comparing polytopes in POLYMESH
  /// @pre: Vertices for each poly are sorted.
  template <typename MTYPE> 
  class POLYMESH_LESS_THAN:public LIST_LESS_THAN<MTYPE> {

  public:
    POLYMESH_LESS_THAN(const MTYPE * pmesh):
      LIST_LESS_THAN<MTYPE>(pmesh)
    {};
  };


  // **************************************************
  // Class POLYMESH member functions
  // **************************************************

  // Return true if 2D polygon ipoly contains edge (iv0,iv1).
  template <typename VTYPE, typename NTYPE>
  template <typename VTYPE0, typename VTYPE1>
  bool POLYMESH<VTYPE,NTYPE>::DoesPoly2DContainEdge
  (const NTYPE ipoly, const VTYPE0 iv0, const VTYPE1 iv1) const
  {
    if (NumPolyVert(ipoly) <= 1) { return(false); }

    for (NTYPE j0 = 0; j0 < NumPolyVert(ipoly); j0++) {
      VTYPE jv0 = Vertex(ipoly, j0);
      NTYPE j1 = (j0+1)%NumPolyVert(ipoly);
      VTYPE jv1 = Vertex(ipoly, j1);

      if (jv0 == iv0 && jv1 == iv1) { return(true); }
      if (jv1 == iv0 && jv0 == iv1) { return(true); }
    }

    return(false);
  }


  // **************************************************
  // Class VERTEX_POLY_INCIDENCE member functions
  // **************************************************

  // Return true if vertex iv is incident on polytope ipoly.
  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename PTYPE2, typename NTYPE2>
  bool VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::
  IsVertexIncidentOnPoly
  (const VTYPE2 iv, const PTYPE2 ipoly, NTYPE2 & iloc) const
  {
    const ETYPE * list = this->List(iv);
    for (iloc = 0; iloc < this->ListLength(iv); iloc++) {
      if (list[iloc].PolyIndex() == ipoly) { return(true); }
    }

    return(false);
  }

  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE0, typename VTYPE1, typename PTYPE2>
  void VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::
  GetPolyContaining
  (const VTYPE0 iv0, const VTYPE1 iv1, std::vector<PTYPE2> & poly_list)
  {
    poly_list.clear();
    for (NTYPE j = 0; j < NumIncidentPoly(iv0);j++) {
      PTYPE ipoly = IncidentPoly(iv0, j);
      if (IsVertexIncidentOnPoly(iv1, ipoly)) 
        { poly_list.push_back(ipoly); }
    }
  }

  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::
  AllocateLists(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    this->num_poly = polymesh.NumPoly();
    this->num_vertices = 0;
    if (this->num_poly > 0) { this->num_vertices = get_max_vert(polymesh)+1; }
    this->SetNumLists(NumVertices());
    count_num_incident_poly(polymesh, this->list_length);
    this->SetFirstElement();
    this->AllocArrayElement();
  }


  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::
  Set(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    std::vector<PTYPE> current_poly_index;
    std::vector<PTYPE> last_poly;
    IJK::PROCEDURE_ERROR error("VERTEX_POLY_INCIDENCE_BASE::Set");

    Clear();
    AllocateLists(polymesh);

    last_poly.resize(NumVertices());
    current_poly_index.resize(NumVertices(), 0);
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      for (NTYPE i = 0; i < polymesh.NumPolyVert(ipoly); i++) {
        const VTYPE2 iv = polymesh.Vertex(ipoly, i);

        if (current_poly_index[iv] != 0) {
          if (last_poly[iv] == ipoly) {
            // Poly ipoly already inserted in list of poly incident on iv.
            continue;
          }
        }

        this->ElementRef(iv, current_poly_index[iv]).SetPolyIndex(ipoly);
        current_poly_index[iv]++;
        last_poly[iv] = ipoly;
      }
    }

    // Check stored correct number of poly for each vertex.
    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (current_poly_index[iv] != this->ListLength(iv)) {
        error.AddMessage
          ("Programming error.  Incorrect storage of poly vertices.");
        error.AddMessage
          ("  num_incident_poly[", iv, "] = ", this->ListLength(iv), ".");
        error.AddMessage
          ("  Stored ", current_poly_index[iv], " incident poly.");
        throw error;
      }
    }


  }

  template <typename ETYPE, typename NTYPE>
  void VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::Clear()
  {
    num_vertices = 0;
    num_poly = 0;
    LIST_OF_LISTS<ETYPE,NTYPE>::Clear();
  }


  // **************************************************
  // Class VERTEX_ADJACENCY_LIST_BASE member functions
  // **************************************************

  // Return true if iv1 is adjacent to iv0.
  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE0, typename VTYPE1>
  bool VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  IsAdjacent(const VTYPE0 iv0, const VTYPE1 iv1) const
  {
    const ETYPE * list = this->List(iv0);
    for (NTYPE i = 0; i < this->ListLength(iv0); i++) {
      if (list[i].Vertex() == iv1) { return(true); }
    }

    return(false);
  }


  template <typename ETYPE, typename NTYPE>
  template <typename NTYPE2>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  SetNumVertices(const NTYPE2 num_vertices)
  {
    this->num_vertices = num_vertices;
    this->SetNumLists(num_vertices);
  }


  // Allocate arrays list_length[], first_element[] and element[].
  // Set number of vertices.
  // Set arrays list_length[] and first_element[].
  template <typename ETYPE, typename NTYPE>
  template <typename NTYPE2>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  AllocateLists(const std::vector<NTYPE2> & num_adjacent)
  {
    num_vertices = num_adjacent.size();

    for (NTYPE i = 0; i < NumVertices(); i++) 
      { this->list_length[i] = num_adjacent[i]; }

    this->SetFirstElement();
    this->AllocArrayElement();
  }


  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    LIST_OF_LISTS<NTYPE,NTYPE> adjacent;
    std::vector<NTYPE> current_element;
    std::vector<NTYPE> num_adjacent;
    IJK::PROCEDURE_ERROR error("VERTEX_ADJACENCY_LIST::SetFrom2DMesh");
 
    Clear();

    if (polymesh.NumPoly() > 0) { 
      const int numv = get_max_vert(polymesh)+1; 
      SetNumVertices(numv);
    }

    if (NumVertices() == 0) { return; };

    adjacent.SetNumLists(NumVertices());
    count_num_poly2D_edges_incident_on_vertices
      (polymesh, adjacent.list_length);

    adjacent.SetFirstElement();
    adjacent.AllocArrayElement();

    current_element.resize(NumVertices());
    for (NTYPE iv = 0; iv < NumVertices(); iv++) 
      { current_element[iv] = adjacent.first_element[iv]; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);
      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          VTYPE2 iv0 = polymesh.Vertex(ipoly, i0);
          NTYPE i1 = (i0+1)%num_poly_vert;
          VTYPE2 iv1 = polymesh.Vertex(ipoly, i1);

          if (iv0 != iv1) {
            adjacent.element[current_element[iv0]] = iv1;
            current_element[iv0]++;
            adjacent.element[current_element[iv1]] = iv0;
            current_element[iv1]++;
          }
        }
      }
    }

    // Check stored correct number of adjacent vertices for each vertex.
    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (current_element[iv] != 
          adjacent.first_element[iv]+adjacent.list_length[iv]) {
        error.AddMessage
          ("Programming error.  Problem computing vertices adjacent to vertex ",
           iv, ".");
        error.AddMessage
          ("  Expected ", adjacent.list_length[iv], 
           " vertices (including duplicates)");
        error.AddMessage
          ("  but computed ", current_element[iv]-adjacent.first_element[iv],
           " adjacent vertices.");
        throw error;
      }
    }

    // Sort adjacent vertices.
    adjacent.SortEachList();

    // Set num_adjacent[i] to number of distinct vertices in adjacent list i.
    num_adjacent.resize(NumVertices());
    adjacent.CountNumDistinct(num_adjacent);

    AllocateLists(num_adjacent);

    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (this->list_length[iv] > 0) {
        NTYPE n = 0;
        NTYPE k = adjacent.first_element[iv];
        this->ElementRef(iv, n).SetVertex(adjacent.element[k]);
        n++;
        for (NTYPE j = 1; j < adjacent.list_length[iv]; j++) {
          if (adjacent.element[k+j] != adjacent.element[k+j-1]) {
            this->ElementRef(iv, n).SetVertex(adjacent.element[k+j]);
            n++;
          }
        }
        if (n != NumAdjacent(iv)) {
          error.AddMessage
            ("Programming error.  Problem computing vertices adjacent to vertex ",
             iv, ".");
          error.AddMessage
            ("  Expected ", NumAdjacent(iv), " vertices (no duplicates)");
          error.AddMessage
            ("  but computed ", n, " adjacent vertices.");
          throw error;
        }
      }
    }

  }


  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  SetFromMeshOfCubes(const POLYMESH<VTYPE2,NTYPE2> & polymesh,
                     const CUBE_TYPE & cube)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = cube.Dimension();
    LIST_OF_LISTS<NTYPE,NTYPE> adjacent;
    std::vector<NTYPE> current_element;
    std::vector<NTYPE> num_adjacent;
    IJK::PROCEDURE_ERROR error("VERTEX_ADJACENCY_LIST::SetFromMeshOfCubes");
 
    Clear();

    if (polymesh.NumPoly() > 0) { 
      const int numv = get_max_vert(polymesh)+1; 
      SetNumVertices(numv);
    }

    if (NumVertices() == 0) { return; };

    adjacent.SetNumLists(NumVertices());
    count_num_cube_edges_incident_on_vertices
      (polymesh, cube, adjacent.list_length);

    adjacent.SetFirstElement();
    adjacent.AllocArrayElement();

    current_element.resize(NumVertices());
    for (NTYPE iv = 0; iv < NumVertices(); iv++) 
      { current_element[iv] = adjacent.first_element[iv]; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);
      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          const VTYPE iv0 = polymesh.Vertex(ipoly, i0);
          for (DTYPE d = 0; d < dimension; d++) {
            const NTYPE i1 = cube.VertexNeighbor(i0,d);
            const VTYPE iv1 = polymesh.Vertex(ipoly, i1);

            if (iv0 < iv1) {
              adjacent.element[current_element[iv0]] = iv1;
              current_element[iv0]++;
              adjacent.element[current_element[iv1]] = iv0;
              current_element[iv1]++;
            }
          }
        }
      }
    }

    // Check stored correct number of adjacent vertices for each vertex.
    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (current_element[iv] != 
          adjacent.first_element[iv]+adjacent.list_length[iv]) {
        error.AddMessage
          ("Programming error.  Problem computing vertices adjacent to vertex ",
           iv, ".");
        error.AddMessage
          ("  Expected ", adjacent.list_length[iv], 
           " vertices (including duplicates)");
        error.AddMessage
          ("  but computed ", current_element[iv]-adjacent.first_element[iv],
           " adjacent vertices.");
        throw error;
      }
    }

    // Sort adjacent vertices.
    adjacent.SortEachList();

    // Set num_adjacent[i] to number of distinct vertices in adjacent list i.
    num_adjacent.resize(NumVertices());
    adjacent.CountNumDistinct(num_adjacent);

    AllocateLists(num_adjacent);

    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (this->list_length[iv] > 0) {
        NTYPE n = 0;
        NTYPE k = adjacent.first_element[iv];
        this->ElementRef(iv, n).SetVertex(adjacent.element[k]);
        n++;
        for (NTYPE j = 1; j < adjacent.list_length[iv]; j++) {
          if (adjacent.element[k+j] != adjacent.element[k+j-1]) {
            this->ElementRef(iv, n).SetVertex(adjacent.element[k+j]);
            n++;
          }
        }
        if (n != NumAdjacent(iv)) {
          error.AddMessage
            ("Programming error.  Problem computing vertices adjacent to vertex ",
             iv, ".");
          error.AddMessage
            ("  Expected ", NumAdjacent(iv), " vertices (no duplicates)");
          error.AddMessage
            ("  but computed ", n, " adjacent vertices.");
          throw error;
        }
      }
    }

  }


  template <typename ETYPE, typename NTYPE>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::Clear()
  {
    num_vertices = 0;
    LIST_OF_LISTS<ETYPE,NTYPE>::Clear();
  }


  // **************************************************
  // Class VERTEX_POLY_EDGE_INCIDENCE member functions
  // **************************************************

  template <typename ETYPE0, typename ETYPE1, typename NTYPE>
  void VERTEX_POLY_EDGE_INCIDENCE<ETYPE0,ETYPE1,NTYPE>::Clear()
  {
    vertex_poly_incidence.Clear();
    vertex_adjacency_list.Clear();
    incident_edge_list.Clear();
  }


  template <typename ETYPE0, typename ETYPE1, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_POLY_EDGE_INCIDENCE<ETYPE0,ETYPE1,NTYPE>::
  SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    // Vertex, incident poly, index of poly in vertex incident list.
    typedef std::tuple<VTYPE,PTYPE,NTYPE> VERTEX_POLY_TUPLE;

    LIST_OF_LISTS<VERTEX_POLY_TUPLE,NTYPE> adjacent;
    std::vector<NTYPE> current_element;
    std::vector<PTYPE> last_poly;
    std::vector<PTYPE> current_poly_index;
    std::vector<NTYPE> current_incident_edge_index;
    NTYPE numv, nump;
    IJK::PROCEDURE_ERROR error("VERTEX_POLY_EDGE_INCIDENCE::SetFrom2DMesh");

    Clear();

    if (polymesh.NumPoly() == 0) { return; }

    nump = polymesh.NumPoly();
    numv = 0;
    if (nump > 0) { numv = get_max_vert(polymesh)+1; }
    vertex_poly_incidence.SetNumLists(numv);
    count_num_incident_poly(polymesh, vertex_poly_incidence.list_length);
    vertex_poly_incidence.SetFirstElement();
    vertex_poly_incidence.AllocArrayElement();

    adjacent.SetNumLists(numv);
    count_num_poly2D_edges_incident_on_vertices(polymesh, adjacent.list_length);
    adjacent.SetFirstElement();
    adjacent.AllocArrayElement();

    current_poly_index.resize(numv, 0);
    current_element.resize(numv);

    for (NTYPE iv = 0; iv < numv; iv++) 
      { current_element[iv] = adjacent.FirstElement(iv); }


    last_poly.resize(numv);
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);

      // Insert ipoly in incident lists for its vertices.
      for (NTYPE i  = 0; i  < num_poly_vert; i ++) {
        const VTYPE2 iv = polymesh.Vertex(ipoly, i);

        if (current_poly_index[iv] != 0) {
          if (last_poly[iv] == ipoly) {
            // Poly ipoly already inserted in list of poly incident on iv.
            continue;
          }
        }

        const NTYPE k = 
          vertex_poly_incidence.ElementIndex(iv, current_poly_index[iv]);
        vertex_poly_incidence.element[k].SetPolyIndex(ipoly);
        current_poly_index[iv]++;
        last_poly[iv] = ipoly;
      }

      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          VTYPE2 iv0 = polymesh.Vertex(ipoly, i0);
          NTYPE i1 = (i0+1)%num_poly_vert;
          VTYPE2 iv1 = polymesh.Vertex(ipoly, i1);

          if (iv0 != iv1) {
            NTYPE k0 = current_poly_index[iv0]-1;
            NTYPE k1 = current_poly_index[iv1]-1;
            adjacent.element[current_element[iv0]] = 
              VERTEX_POLY_TUPLE(iv1,ipoly,k0);
            current_element[iv0]++;
            adjacent.element[current_element[iv1]] =
              VERTEX_POLY_TUPLE(iv0,ipoly,k1);
            current_element[iv1]++;
          }
        }
      }
    }

    // Check stored correct number of adjacent vertices for each vertex.
    for (NTYPE iv = 0; iv < numv; iv++) {
      if (current_element[iv] != adjacent.FirstElement(iv) + adjacent.ListLength(iv)) {
        error.AddMessage
          ("Programming error.  Problem computing vertices adjacent to vertex ",
           iv, ".");
        error.AddMessage
          ("  Expected ", adjacent.ListLength(iv), 
           " vertices (including duplicates)");
        error.AddMessage
          ("  but computed ", current_element[iv]-adjacent.FirstElement(iv), " adjacent vertices.");
        throw error;
      }
    }

    // Sort adjacent vertices.
    adjacent.SortEachList();

    vertex_adjacency_list.SetNumLists(numv);
    count_num_distinct_tuple0(adjacent, vertex_adjacency_list.list_length);
    vertex_adjacency_list.SetFirstElement();
    vertex_adjacency_list.AllocArrayElement();

    incident_edge_list.SetNumLists(vertex_poly_incidence.element.size());

    // Set arrays vertex_adjacency_list.element[] and 
    //   incident_edge_list.list_length[].
    for (NTYPE iv0 = 0; iv0 < NumVertices(); iv0++) {
      if (vertex_adjacency_list.ListLength(iv0) > 0) {
        NTYPE n = 0;
        NTYPE ifirst = adjacent.first_element[iv0];
        VTYPE2 iv1 = std::get<0>(adjacent.element[ifirst]);
        vertex_adjacency_list.ElementRef(iv0,n).SetVertex(iv1);
        NTYPE j = std::get<2>(adjacent.element[ifirst]);
        IncrementIncidentEdgeListLength(iv0, j);
        n++;
        for (NTYPE i = 1; i < adjacent.list_length[iv0]; i++) {
          VTYPE2 iv1 = std::get<0>(adjacent.element[ifirst+i]);
          NTYPE j1 = std::get<2>(adjacent.element[ifirst+i]);
          VTYPE2 iv2 = std::get<0>(adjacent.element[ifirst+i-1]);
          NTYPE j2 = std::get<2>(adjacent.element[ifirst+i-1]);
          if (iv1 != iv2) {
            vertex_adjacency_list.ElementRef(iv0,n).SetVertex(iv1);
            n++;
          }

          if (iv1 != iv2 || j1 != j2) 
            { IncrementIncidentEdgeListLength(iv0, j1); }
        }
        if (n != vertex_adjacency_list.ListLength(iv0)) {
          error.AddMessage
            ("Programming error.  Problem computing vertices adjacent to vertex ",
             iv0, ".");
          error.AddMessage
            ("  Expected ", vertex_adjacency_list.ListLength(iv0), 
             " vertices (no duplicates)");
          error.AddMessage
            ("    but computed ", n, " adjacent vertices.");
          throw error;
        }
      }
    }

    incident_edge_list.SetFirstElement();
    incident_edge_list.AllocArrayElement();

    current_incident_edge_index.resize(incident_edge_list.NumLists(), 0);

    // Set array incident_edge_list.element[].
    for (NTYPE iv0 = 0; iv0 < NumVertices(); iv0++) {
      if (vertex_adjacency_list.ListLength(iv0) > 0) {
        NTYPE n = 0;
        NTYPE ifirst = adjacent.first_element[iv0];
        VTYPE2 iv1 = std::get<0>(adjacent.element[ifirst]);
        NTYPE j = std::get<2>(adjacent.element[ifirst]);
        NTYPE k = IncidentEdgeListIndex(iv0,j);
        SetIncidentEdgeListElement(iv0, j, current_incident_edge_index[k], n);
        current_incident_edge_index[k]++;
        for (NTYPE i = 1; i < adjacent.list_length[iv0]; i++) {
          VTYPE2 iv1 = std::get<0>(adjacent.element[ifirst+i]);
          NTYPE j1 = std::get<2>(adjacent.element[ifirst+i]);
          VTYPE2 iv2 = std::get<0>(adjacent.element[ifirst+i-1]);
          NTYPE j2 = std::get<2>(adjacent.element[ifirst+i-1]);
          if (iv1 != iv2) { n++; }
          if (iv1 != iv2 || j1 != j2) {
            NTYPE k = IncidentEdgeListIndex(iv0,j1);
            SetIncidentEdgeListElement
              (iv0, j1, current_incident_edge_index[k], n);
            current_incident_edge_index[k]++;
          }
        }
        if (n+1 != vertex_adjacency_list.ListLength(iv0)) {
          error.AddMessage
            ("Programming error.  Problem setting incident edges for vertex ",
             iv0, ".");
          error.AddMessage
            ("  Expected ", vertex_adjacency_list.ListLength(iv0), 
             " vertices (no duplicates)");
          error.AddMessage
            ("    but computed ", n+1, " adjacent vertices.");
          throw error;
        }
      }
    }

    // Check stored correct number of edges for each incident edge list.
    for (NTYPE k = 0; k < incident_edge_list.NumLists(); k++) {
      if (current_incident_edge_index[k] != incident_edge_list.ListLength(k)) {
        error.AddMessage
          ("Programming error.  Incorrect storage of incident edges.");
        error.AddMessage
          ("  num_incident_edges[", k, "] = ", 
           incident_edge_list.ListLength(k), ".");
        error.AddMessage
          ("  Stored ", current_poly_index[k], " incident edges.");
        throw error;
      }
    }

  }

  // Return true if list of edges incident on vertex iv0 and 
  //  poly vertex_poly_incidence(iv,j) contains edge (iv0, iv1).
  template <typename ETYPE0, typename ETYPE1, typename NTYPE>
  template <typename VTYPE0, typename VTYPE1, typename NTYPE2>
  bool VERTEX_POLY_EDGE_INCIDENCE<ETYPE0,ETYPE1,NTYPE>::
  DoesIncidenceListContainEdge
  (const VTYPE0 iv0, const VTYPE1 iv1, const NTYPE2 j) const
  {
    const NTYPE2 k = IncidentEdgeListIndex(iv0, j);
    for (NTYPE2 i = 0; i < incident_edge_list.ListLength(k); i++) {
      NTYPE2 n = incident_edge_list.Element(k, i);
      if (vertex_adjacency_list.Element(iv0, n).Vertex() == iv1) 
        { return(true); }
    }

    return(false);
  }


  // **************************************************
  // Functions on POLYMESH
  // **************************************************

  /// Return maximum poly vertex index.
  /// Return 0 if there are no polytope vertices.
  template <typename VTYPE, typename NTYPE>
  const VTYPE get_max_vert(const POLYMESH<VTYPE,NTYPE> & polymesh)
  {
    if (polymesh.element.size() > 0) {

      const VTYPE vmax = 
        *(std::max_element
          (polymesh.element.begin(), polymesh.element.end()));

      return(vmax);
    }
    else {
      return(0);
    }
  }

  /// Return sum of number of vertices in all polytopes.
  template <typename VTYPE, typename NTYPE>
  const NTYPE sum_num_poly_vert
  (const POLYMESH<VTYPE,NTYPE> & polymesh)
  {
    NTYPE sum = 0;
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) 
      { sum += polymesh.NumPolyVert(ipoly); }

    return(sum);
  }

  /// Count number of polytopes incident on each vertex.
  /// @param[out] num_incident_poly[] Array of number of incident polytopes.
  /// @pre num_incident_poly is preallocated to size = number of vertices.
  template <typename MESH_TYPE, typename NTYPE>
  void count_num_incident_poly
  (const MESH_TYPE & polymesh, std::vector<NTYPE> & num_incident_poly)
  {
    const NTYPE numv = num_incident_poly.size();
    std::vector<NTYPE> last_poly(numv);

    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < numv; iv++) 
      { num_incident_poly[iv] = 0; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      for (NTYPE i = 0; i < polymesh.NumPolyVert(ipoly); i++) {
        VTYPE iv = polymesh.Vertex(ipoly, i);

        if (num_incident_poly[iv] > 0) {
          if (last_poly[iv] == ipoly) {
            // Poly ipoly already counted as incident on iv.
            continue;
          }
        }

        num_incident_poly[iv]++;
        last_poly[iv] = ipoly;
      }
    }
  }


  /// Count number of edges incident on each vertex.
  /// If same edge is in k polygons, counts edge k times.
  /// @param polymesh Mesh of 2D polygons.  Vertices are listed in cyclic
  ///   order around each polygon.
  /// @param[out] num_incident_edges[] Array of number of incident edges.
  /// @pre num_incident_edges is preallocated to size = number of vertices.
  template <typename MESH_TYPE, typename NTYPE>
  void count_num_poly2D_edges_incident_on_vertices
  (const MESH_TYPE & polymesh, std::vector<NTYPE> & num_incident_edges)
  {
    const NTYPE numv = num_incident_edges.size();

    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < numv; iv++) 
      { num_incident_edges[iv] = 0; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);
      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          VTYPE iv0 = polymesh.Vertex(ipoly, i0);
          NTYPE i1 = (i0+1)%num_poly_vert;
          VTYPE iv1 = polymesh.Vertex(ipoly, i1);
          if (iv0 != iv1) {
            num_incident_edges[iv0]++;
            num_incident_edges[iv1]++;
          }
        }
      }
    }
  }

  /// Count number of (hyper)cube edges incident on each vertex.
  /// If same edge is in k polygons, counts edge k times.
  /// @param polymesh Mesh of (hyper) cubes.  
  ///   - Vertices are listed in order: 
  ///  (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1).
  /// @param[out] num_incident_edges[] Array of number of incident edges.
  /// @pre num_incident_edges is preallocated to size = number of vertices.
  template <typename MESH_TYPE, typename CUBE_TYPE, typename NTYPE>
  void count_num_cube_edges_incident_on_vertices
  (const MESH_TYPE & polymesh, const CUBE_TYPE & cube,
   std::vector<NTYPE> & num_incident_edges)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = cube.Dimension();
    const NTYPE numv = num_incident_edges.size();

    for (VTYPE iv = 0; iv < numv; iv++) 
      { num_incident_edges[iv] = 0; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);
      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          const VTYPE iv0 = polymesh.Vertex(ipoly, i0);
          for (NTYPE d = 0; d < dimension; d++) {
            const NTYPE i1 = cube.VertexNeighbor(i0,d);
            const VTYPE iv1 = polymesh.Vertex(ipoly, i1);

            if (iv0 < iv1) {
              num_incident_edges[iv0]++;
              num_incident_edges[iv1]++;
            }
          }
        }
      }
    }
  }


  /// Compute vertex link in mesh of cubes.
  /// @param cube Cube information.
  /// @pre cube has same dimension as cubes in mesh.
  /// @param[out] link_mesh Polygonal mesh of links.
  ///   Note that quadrilaterals are stored lower-left, lower-right,
  ///   upper-left, upper_right, NOT in clockwise or counter-clockwise order.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename VTYPE, typename CUBE_TYPE, typename POLYMESH2_TYPE>
  void compute_vertex_link_in_cube_mesh
  (const POLYMESH_TYPE & polymesh, 
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
   const VTYPE iv,
   const CUBE_TYPE & cube,
   POLYMESH2_TYPE & link_mesh)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE num_vertices_per_facet = cube.NumFacetVertices();
    std::vector<VTYPE> facet_vertex(num_vertices_per_facet);

    link_mesh.Clear();

    for (NTYPE j = 0; j < vertex_info.NumIncidentPoly(iv); j++) {
      const NTYPE jpoly = vertex_info.IncidentPoly(iv,j);

      for (NTYPE jfacet = 0; jfacet < cube.NumFacets(); jfacet++) {

        for (NTYPE k = 0; k < num_vertices_per_facet; k++) {
          const NTYPE k2 = cube.FacetVertex(jfacet, k);
          const VTYPE kv = polymesh.Vertex(jpoly, k2);
          facet_vertex[k] = kv;
        }

        if (!IJK::does_list_contain(facet_vertex, iv)) {
          link_mesh.AddPolytope(facet_vertex);
        }
      }
    }
  }

}

#endif
