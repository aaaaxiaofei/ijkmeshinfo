/// \file ijkmesh_faces.txx
/// ijk templates for processing polyhedral meshes faces.
/// - Version 0.2.1

/*
  IJK: Isosurface Jeneration Kode
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

#ifndef _IJKMESH_FACES_
#define _IJKMESH_FACES_

#include "ijk.txx"
#include "ijklist.txx"
#include "ijkmesh_datastruct.txx"

#include <vector>

namespace IJK {


  // **************************************************
  // ORIENTATION FUNCTIONS
  // **************************************************

  template <typename POLYMESH_TYPE, 
            typename FACET_INFO_TYPEA, typename FACET_INFO_TYPEB,
            typename CUBE_TYPE>
  bool do_hexahedra_orientations_match
  (const POLYMESH_TYPE & polymesh, 
   const FACET_INFO_TYPEA & facetA, const FACET_INFO_TYPEB & facetB,
   const CUBE_TYPE & cube)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename POLYMESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const NTYPE NUM_VERT_PER_CUBE_FACET(cube.NumFacetVertices());
    VTYPE facetA_list[NUM_VERT_PER_CUBE_FACET];
    VTYPE facetB_list[NUM_VERT_PER_CUBE_FACET];

    const NTYPE ipolyA = facetA.poly_containing_face;
    const NTYPE ipolyB = facetB.poly_containing_face;
    const NTYPE ifacetA = facetA.face_index;
    const NTYPE ifacetB = facetB.face_index;
    for (NTYPE i = 0; i < NUM_VERT_PER_CUBE_FACET; i++) {
      const NTYPE jA = cube.FacetVertex(ifacetA, i);
      const NTYPE jB = cube.FacetVertex(ifacetB, i);
      facetA_list[i] = polymesh.Vertex(ipolyA, jA);
      facetB_list[i] = polymesh.Vertex(ipolyB, jB);
    }

    std::swap(facetA_list[2], facetA_list[3]);
    std::swap(facetB_list[2], facetB_list[3]);


    // Location of facetB_list[0] in facetA_list.
    NTYPE jloc0;
    if (!does_list_contain
        (facetA_list, NUM_VERT_PER_CUBE_FACET, facetB_list[0], jloc0)) {
      // facetA and facetB are not the same.
      // Cannot compare orientation.
      return(true);
    }

    const NTYPE jloc1 = (jloc0+1)%NUM_VERT_PER_CUBE_FACET;

    if (facetA_list[jloc1] == facetB_list[1]) {
      // Polytopes have opposite orientation.
      return(false); 
    }
    else
      { return(true); }
  }


  /// Return true if orientation of tetrahedron facetA matches orientation
  ///   of tetrahedron facetB.
  /// @pre facetA and facetB have same set of vertices.
  /// @pre polymesh has at least one simplex.
  template <typename POLYMESH_TYPE, 
            typename FACET_INFO_TYPEA, typename FACET_INFO_TYPEB>
  bool do_tetrahedra_orientations_match
  (const POLYMESH_TYPE & polymesh, 
   const FACET_INFO_TYPEA & facetA, const FACET_INFO_TYPEB & facetB)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename POLYMESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const NTYPE NUM_VERT_PER_TET(4);
    const NTYPE NUM_VERT_PER_FACET(NUM_VERT_PER_TET-1);
    VTYPE facetA_list[NUM_VERT_PER_FACET];
    VTYPE facetB_list[NUM_VERT_PER_FACET];

    const NTYPE ipolyA = facetA.poly_containing_face;
    const NTYPE ipolyB = facetB.poly_containing_face;
    const NTYPE ifacetA = facetA.face_index;
    const NTYPE ifacetB = facetB.face_index;
    for (NTYPE i = 0; i < NUM_VERT_PER_FACET; i++) {
      const NTYPE jA = (ifacetA+i+1)%NUM_VERT_PER_TET;
      const NTYPE jB = (ifacetB+i+1)%NUM_VERT_PER_TET;
      facetA_list[i] = polymesh.Vertex(ipolyA, jA);
      facetB_list[i] = polymesh.Vertex(ipolyB, jB);
    }

    // @pre facetA_list[] and facetB_list[] contain the same set of vertices.
    if (facetA_list[0] == facetB_list[0]) {
      if (facetA_list[1] == facetB_list[1]) { return(true); }
      else { return(false); }
    }
    else if (facetA_list[0] == facetB_list[1]) {
      if (facetA_list[1] == facetB_list[2]) { return(true); }
      else { return(false); }
    }
    else {
      // facetA_list[0] == facetB_list[2]

      if (facetA_list[1] == facetB_list[0]) { return(true); }
      else { return(false); }
    }

  }


  // **************************************************
  // SUBROUTINES FOR PROCESSING FACES
  // **************************************************

  namespace MESH_FACES {

    /// Store vertices of facet kf in C++ vector facet_vert.
    /// Store facet kf in C++ vector facet_list.
    template <typename FACET_LIST_TYPE, typename FACET_INDEX_TYPE, 
              typename VTYPE, typename FACE_INFO_TYPE>
    void store_facet
    (const FACET_LIST_TYPE & facets, const FACET_INDEX_TYPE kf,
     std::vector<VTYPE> & facet_vert, 
     std::vector<FACE_INFO_TYPE> & facet_list)
    {
      typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
      typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;

      // Store facet vertices.
      add_list(facets.VertexList(kf), facets.NumFacetVert(kf), facet_vert);

      // Store facet.
      const PTYPE jpoly = facets.PolyContainingFacet(kf);
      const FTYPE jf = facets.FacetIndex(kf);
      facet_list.push_back(FACE_INFO_TYPE(jpoly, jf));
    }


    /// Store vertices of facet kf0 in C++ vector facet_vert.
    /// Store facets kf0 and kf1 in C++ vector facet_list.
    template <typename FACET_LIST_TYPE, typename FACET_INDEX_TYPE, 
              typename VTYPE, typename FACE_INFO_TYPE>
    void store_two_facets
    (const FACET_LIST_TYPE & facets, 
     const FACET_INDEX_TYPE kf0, const FACET_INDEX_TYPE kf1,
     std::vector<VTYPE> & facet_vert, 
     std::vector<FACE_INFO_TYPE> & facet_list)
    {
      typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
      typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;

      // Store vertices of facet kf0.
      add_list(facets.VertexList(kf0), facets.NumFacetVert(kf0), facet_vert);

      // Store facets.
      const PTYPE jpoly0 = facets.PolyContainingFacet(kf0);
      const FTYPE jf0 = facets.FacetIndex(kf0);
      facet_list.push_back(FACE_INFO_TYPE(jpoly0, jf0));
      const PTYPE jpoly1 = facets.PolyContainingFacet(kf1);
      const FTYPE jf1 = facets.FacetIndex(kf1);
      facet_list.push_back(FACE_INFO_TYPE(jpoly1, jf1));
    }


    /// Store vertices of facet sorted_facet[k0] in C++ vector facet_vert.
    /// Store facets sorted_facet[k0], sorted_facet[k0+1], ...
    ///   sorted_facet[k1-1] in C++ vector facet_list.
    /// @pre k0 < k1.
    template <typename FACET_LIST_TYPE, typename FACET_INDEX_TYPE,
              typename NTYPE0, typename NTYPE1,
              typename VTYPE, typename FACE_INFO_TYPE>
    void store_facets_in_sorted_range
    (const FACET_LIST_TYPE & facets,
     const std::vector<FACET_INDEX_TYPE> & sorted_facet,
     const NTYPE0 k0, const NTYPE1 k1,
     std::vector<VTYPE> & facet_vert, 
     std::vector<FACE_INFO_TYPE> & facet_list)
    {
      typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
      typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;

      const NTYPE0 kf0 = sorted_facet[k0];

      // Store vertices of facet kf0.
      add_list(facets.VertexList(kf0), facets.NumFacetVert(kf0), facet_vert);

      // Store facets.
      for (int m = k0; m < k1; m++) {
        const FACET_INDEX_TYPE mf = sorted_facet[m];
        const PTYPE jpoly = facets.PolyContainingFacet(mf);
        const FTYPE jf = facets.FacetIndex(mf);
        facet_list.push_back(FACE_INFO_TYPE(jpoly, jf));
      }
    }
  }


  // **************************************************
  // GET NON-MANIFOLD AND BOUNDARY FACETS
  // **************************************************

  /// Get non-manifold and boundary facets of hex mesh.
  template <typename POLYMESH_TYPE, typename CUBE_TYPE,
            typename VTYPE, typename FACE_INFO_TYPE>
  void get_non_manifold_and_boundary_facets_of_hex_mesh
    (const POLYMESH_TYPE & polymesh, const CUBE_TYPE & cube,
     std::vector<VTYPE> & non_manifold_facet_vert,
     std::vector<FACE_INFO_TYPE> & non_manifold_facets,
     std::vector<VTYPE> & boundary_facet_vert,
     std::vector<FACE_INFO_TYPE> & boundary_facets,
     std::vector<VTYPE> & orientation_mismatch_facet_vert,
     std::vector<FACE_INFO_TYPE> & orientation_mismatch_facets)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;
    typedef FACE_INFO_BASE<PTYPE,NTYPE> FINFO_TYPE;

    FACET_LIST_BASE<VTYPE,NTYPE,FINFO_TYPE> facets;
    std::vector<NTYPE> sorted_facet;

    using namespace IJK::MESH_FACES;

    facets.SetFromMeshOfCubes(polymesh, cube);
    facets.SortVert();
    facets.GetSortedPolytopeIndices(sorted_facet);

    NTYPE k0 = 0;
    while (k0 < facets.NumFacets()) {
      NTYPE kf0 = sorted_facet[k0];
      NTYPE k1 = k0+1;
      while (k1 < facets.NumFacets()) {
        const NTYPE kf1 = sorted_facet[k1];
        if (!facets.ArePolytopesEqual(kf0, kf1))
          { break; }
        k1++;
      }

      NTYPE num_duplicate = k1 - k0;
      if (num_duplicate > 2) {
        store_facets_in_sorted_range
          (facets, sorted_facet, k0, k1, non_manifold_facet_vert,
           non_manifold_facets);
      }
      else if (num_duplicate == 1) {
        store_facet(facets, kf0, boundary_facet_vert, boundary_facets);
      }
      else if (num_duplicate == 2) {

        const NTYPE kf0 = sorted_facet[k0];
        const NTYPE kf1 = sorted_facet[k0+1];
        if (!do_hexahedra_orientations_match
            (polymesh, facets.poly_data[kf0], facets.poly_data[kf1], cube)) {

          store_two_facets(facets, kf0, kf1, orientation_mismatch_facet_vert,
                           orientation_mismatch_facets);
        }
      }

      k0 = k1;
    }

  }


  /// Get non-manifold and boundary facets of tetrahedral mesh.
  template <typename POLYMESH_TYPE,
            typename VTYPE, typename FACE_INFO_TYPE>
  void get_non_manifold_and_boundary_facets_of_tet_mesh
    (const POLYMESH_TYPE & polymesh,
     std::vector<VTYPE> & non_manifold_facet_vert,
     std::vector<FACE_INFO_TYPE> & non_manifold_facets,
     std::vector<VTYPE> & boundary_facet_vert,
     std::vector<FACE_INFO_TYPE> & boundary_facets,
     std::vector<VTYPE> & orientation_mismatch_facet_vert,
     std::vector<FACE_INFO_TYPE> & orientation_mismatch_facets)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;
    typedef FACE_INFO_BASE<PTYPE,NTYPE> FINFO_TYPE;

    FACET_LIST_BASE<VTYPE,NTYPE,FINFO_TYPE> facets;
    std::vector<NTYPE> sorted_facet;

    using namespace IJK::MESH_FACES;

    facets.SetFromMeshOfSimplices(polymesh);
    facets.SortVert();
    facets.GetSortedPolytopeIndices(sorted_facet);

    NTYPE k0 = 0;
    while (k0 < facets.NumFacets()) {
      NTYPE kf0 = sorted_facet[k0];
      NTYPE k1 = k0+1;
      while (k1 < facets.NumFacets()) {
        const NTYPE kf1 = sorted_facet[k1];
        if (!facets.ArePolytopesEqual(kf0, kf1))
          { break; }
        k1++;
      }

      NTYPE num_duplicate = k1 - k0;
      if (num_duplicate > 2) {
        store_facets_in_sorted_range
          (facets, sorted_facet, k0, k1, non_manifold_facet_vert,
           non_manifold_facets);
      }
      else if (num_duplicate == 1) {
        store_facet(facets, kf0, boundary_facet_vert, boundary_facets);
      }
      else if (num_duplicate == 2) {

        const NTYPE kf0 = sorted_facet[k0];
        const NTYPE kf1 = sorted_facet[k0+1];
        if (!do_tetrahedra_orientations_match
            (polymesh, facets.poly_data[kf0], facets.poly_data[kf1])) {

          store_two_facets(facets, kf0, kf1, orientation_mismatch_facet_vert,
                           orientation_mismatch_facets);
        }
      }

      k0 = k1;
    }

  }


  // **************************************************
  // GET NON-MANIFOLD VERTICES
  // **************************************************

}

#endif

